/*!
====================================================================
|                     MicromorphicMaterial.cpp                     |
====================================================================
| The source file for a class which computes the cauchy stress and |
| the associated jacobians for a micromorphic material.            |
--------------------------------------------------------------------
| Notes: Relies on libraries from the micromorphic_element         |
|        repository available at bitbucket.org/NateAM2             |
====================================================================
*/


#include "MicromorphicMaterial.h"

registerMooseObject("tartigrade", MicromorphicMaterial);

template<>
InputParameters
validParams<MicromorphicMaterial>(){
    InputParameters params = validParams<Material>();

    // Vectors of material properties
    params.addRequiredParam<std::vector<Real>>(
        "material_fparameters", "The vector of floating point material parameters required for the stiffness matrices");

    params.addRequiredParam<std::string>(
        "model_name", "The material model name");

    params.addRequiredParam<int>(
        "number_ADD_DOF", 0, "The number of additional degrees of freedom beyond u and phi in the problem");

    params.addRequiredParam<int>(
        "number_ADD_TERMS", 0, "The number of additional balance equations being solved beyond the balance of linear momentum and first moment of momentum");

    params.addRequiredParam<int>(
        "number_ADD_JACOBIANS", 0, "The number of additional jacobians being provided beyond those of the stress measures");

    return params;
}

MicromorphicMaterial::MicromorphicMaterial(const InputParameters & parameters)
    : Material(parameters),
    // Declare that this material is going to provide Eigen matrices containing the cauchy stress and 
    // jacobians that Kernels can use.
    _cauchy(declareProperty<Vector_9>("cauchy")),
    _s(declareProperty<Vector_9>("s")),
    _m(declareProperty<Vector_27>("m")),
    _DcauchyDgrad_u(declareProperty<Matrix_9x9>("DcauchyDgrad_u")),
    _DcauchyDphi(declareProperty<Matrix_9x9>("DcauchyDphi")),
    _DcauchyDgrad_phi(declareProperty<Matrix_9x27>("DcauchyDgrad_phi"))
    _DsDgrad_u(declareProperty<Matrix_9x9>("DsDgrad_u")),
    _DsDphi(declareProperty<Matrix_9x9>("DsDphi")),
    _DsDgrad_phi(declareProperty<Matrix_9x27>("DsDgrad_phi")),
    _DmDgrad_u(declareProperty<Matrix_27x9>("DmDgrad_u")),
    _DmDphi(declareProperty<Matrix_27x9>("DmDphi")),
    _DmDgrad_phi(declareProperty<Matrix_27x27>("DmDgrad_phi")),
    _fparameters(getParam<std::vector<Real>>("material_fparameters")),
    _n_ADD_DOF(getParam<int>("number_ADD_DOF")),
    _n_ADD_TERMS(getParam<int>("number_ADD_TERMS")),
    _n_ADD_JACOBIANS(getParam<int>("number_ADD_JACOBIANS")),
    _model_name(getParam<std::string>("model_name")){}

void MicromorphicMaterial::computeQpProperties(){
    /*!
    =============================
    |    computeQpProperties    |
    =============================

    Evaluate the constitutive model at the quadrature point.

    */

    //Define required DOF values for the
    //balance of linear momentum and first moment of momentum
    double __grad_u[3][3];
    double __phi[9];
    double __grad_phi[9][3];
    
    //Copy over the gradient of u
    for (int i=0; i<3; i++){__grad_u[0][i] = _grad_u1[_qp][i];}
    for (int i=0; i<3; i++){__grad_u[1][i] = _grad_u2[_qp][i];}
    for (int i=0; i<3; i++){__grad_u[2][i] = _grad_u3[_qp][i];}
    
    //Copy over phi
    __phi[0] = _phi_11[_qp];
    __phi[1] = _phi_22[_qp];
    __phi[2] = _phi_33[_qp];
    __phi[3] = _phi_23[_qp];
    __phi[4] = _phi_13[_qp];
    __phi[5] = _phi_12[_qp];
    __phi[6] = _phi_32[_qp];
    __phi[7] = _phi_31[_qp];
    __phi[8] = _phi_21[_qp];

    //Copy over grad_phi
    for (int i=0; i<3; i++){__grad_phi[0][i] = _grad_phi_11[_qp][i];}
    for (int i=0; i<3; i++){__grad_phi[1][i] = _grad_phi_22[_qp][i];}
    for (int i=0; i<3; i++){__grad_phi[2][i] = _grad_phi_33[_qp][i];}
    for (int i=0; i<3; i++){__grad_phi[3][i] = _grad_phi_23[_qp][i];}
    for (int i=0; i<3; i++){__grad_phi[4][i] = _grad_phi_13[_qp][i];}
    for (int i=0; i<3; i++){__grad_phi[5][i] = _grad_phi_12[_qp][i];}
    for (int i=0; i<3; i++){__grad_phi[6][i] = _grad_phi_32[_qp][i];}
    for (int i=0; i<3; i++){__grad_phi[7][i] = _grad_phi_31[_qp][i];}
    for (int i=0; i<3; i++){__grad_phi[8][i] = _grad_phi_21[_qp][i];}

    //Evaluate the stresses and their jacobians
    auto &factory = micromorphic_material_library::MaterialFactory::Instance();
    auto material = factory.GetMaterial(_model_name);

    //TODO: Remove these hardcoded values
    std::vector<double> time;
    time.resize(2);
    std::vector<double> SDVS;
    SDVS.resize(0);

    std::vector<double> ADD_DOF;
    ADD_DOF.resize(_n_ADD_DOF);

    std::vector<std::vector<double>> ADD_grad_DOF;
    ADD_grad_DOF.resize(_n_ADD_DOF);
    //END of hardcoded values

    //Set the sizes of the 
    std::vector<Eigen::VectorXd> ADD_TERMS;
    ADD_terms.resize(_n_ADD_TERMS);

    std::vector<Eigen::VectorXd> ADD_JACOBIANS;
    ADD_JACOBIANS.resize(_n_ADD_JACOBIANS);

    material->evaluate_model(time, fparams, __grad_u, __phi, __grad_phi,  SDVS, ADD_DOF, ADD_grad_DOF, 
                             _cauchy[_qp], _s[_qp], _m[_qp],
                             _DcauchyDgrad_u[_qp], _DcauchyDphi[_qp], _DcauchyDgrad_phi[_qp],
                             _DsDgrad_u[_qp],      _DsDphi[_qp],      _DsDgrad_phi[_qp],
                             _DmDgrad_u[_qp],      _DmDphi[_qp],      _DmDgrad_phi[_qp]);
}

