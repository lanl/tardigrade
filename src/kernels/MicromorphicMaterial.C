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

registerMooseObject("tardigradeApp", MicromorphicMaterial);

template<>
InputParameters
validParams<MicromorphicMaterial>(){
    InputParameters params = validParams<Material>();

    // Vectors of material properties
    params.addRequiredParam<std::vector<Real>>(
        "material_fparameters", "The vector of floating point material parameters required for the stiffness matrices");

    params.addRequiredParam<std::string>(
        "model_name", "The material model name");

    params.addParam<int>(
        "number_ADD_DOF", 0, "The number of additional degrees of freedom beyond u and phi in the problem");

    params.addParam<int>(
        "number_ADD_TERMS", 0, "The number of additional balance equations being solved beyond the balance of linear momentum and first moment of momentum");

    params.addParam<int>(
        "number_ADD_JACOBIANS", 0, "The number of additional jacobians being provided beyond those of the stress measures");

    // Coupled variables
    params.addCoupledVar(
        "u1", "The displacement in the x direction.");

    params.addCoupledVar(
        "u2", "The displacement in the y direction.");

    params.addCoupledVar(
        "u3", "The displacement in the z direction.");

    params.addCoupledVar(
        "phi_11", "The xx component of the phi tensor.");

    params.addCoupledVar(
        "phi_22", "The yy component of the phi tensor.");

    params.addCoupledVar(
        "phi_33", "The zz component of the phi tensor.");

    params.addCoupledVar(
        "phi_23", "The yz component of the phi tensor.");

    params.addCoupledVar(
        "phi_13", "The xz component of the phi tensor.");

    params.addCoupledVar(
        "phi_12", "The xy component of the phi tensor.");

    params.addCoupledVar(
        "phi_32", "The zy component of the phi tensor.");

    params.addCoupledVar(
        "phi_31", "The zx component of the phi tensor.");

    params.addCoupledVar(
        "phi_21", "The yx component of the phi tensor.");

    return params;
}

MicromorphicMaterial::MicromorphicMaterial(const InputParameters & parameters)
    : Material(parameters),
    // Declare that this material is going to provide Eigen matrices containing the cauchy stress and 
    // jacobians that Kernels can use.
    _fparams(getParam<std::vector<Real>>("material_fparameters")),
    _n_ADD_DOF(getParam<int>("number_ADD_DOF")),
    _n_ADD_TERMS(getParam<int>("number_ADD_TERMS")),
    _n_ADD_JACOBIANS(getParam<int>("number_ADD_JACOBIANS")),
    _model_name(getParam<std::string>("model_name")),
    _grad_u1(isCoupled("u1") ? coupledGradient("u1")
                             : _grad_zero),
    _grad_u2(isCoupled("u2") ? coupledGradient("u2")
                             : _grad_zero),
    _grad_u3(isCoupled("u3") ? coupledGradient("u2")
                             : _grad_zero),
    _phi_11(coupledValue("phi_11")),
    _phi_22(coupledValue("phi_22")),
    _phi_33(coupledValue("phi_33")),
    _phi_23(coupledValue("phi_23")),
    _phi_13(coupledValue("phi_13")),
    _phi_12(coupledValue("phi_12")),
    _phi_32(coupledValue("phi_32")),
    _phi_31(coupledValue("phi_31")),
    _phi_21(coupledValue("phi_21")),
    _grad_phi_11(isCoupled("phi_11") ? coupledGradient("phi_11")
                                     : _grad_zero),
    _grad_phi_22(isCoupled("phi_22") ? coupledGradient("phi_22")
                                     : _grad_zero),
    _grad_phi_33(isCoupled("phi_33") ? coupledGradient("phi_33")
                                     : _grad_zero),
    _grad_phi_23(isCoupled("phi_23") ? coupledGradient("phi_23")
                                     : _grad_zero),
    _grad_phi_13(isCoupled("phi_13") ? coupledGradient("phi_13")
                                     : _grad_zero),
    _grad_phi_12(isCoupled("phi_12") ? coupledGradient("phi_12")
                                     : _grad_zero),
    _grad_phi_32(isCoupled("phi_32") ? coupledGradient("phi_32")
                                     : _grad_zero),
    _grad_phi_31(isCoupled("phi_31") ? coupledGradient("phi_31")
                                     : _grad_zero),
    _grad_phi_21(isCoupled("phi_21") ? coupledGradient("phi_21")
                                     : _grad_zero),
    _cauchy(declareProperty<std::vector<double>>("cauchy")),
    _s(declareProperty<std::vector<double>>("s")),
    _m(declareProperty<std::vector<double>>("m")),
    _DcauchyDgrad_u(declareProperty<std::vector<std::vector<double>>>("DcauchyDgrad_u")),
    _DcauchyDphi(declareProperty<std::vector<std::vector<double>>>("DcauchyDphi")),
    _DcauchyDgrad_phi(declareProperty<std::vector<std::vector<double>>>("DcauchyDgrad_phi")),
    _DsDgrad_u(declareProperty<std::vector<std::vector<double>>>("DsDgrad_u")),
    _DsDphi(declareProperty<std::vector<std::vector<double>>>("DsDphi")),
    _DsDgrad_phi(declareProperty<std::vector<std::vector<double>>>("DsDgrad_phi")),
    _DmDgrad_u(declareProperty<std::vector<std::vector<double>>>("DmDgrad_u")),
    _DmDphi(declareProperty<std::vector<std::vector<double>>>("DmDphi")),
    _DmDgrad_phi(declareProperty<std::vector<std::vector<double>>>("DmDgrad_phi")),
    _ADD_TERMS(declareProperty<std::vector<std::vector<double>>>("ADD_TERMS")){
    /*!
    ==============================
    |    MicromorphicMaterial    |
    ==============================

    The constructor for MicromorphicMaterial.
    */


}

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
    bool print_grad_u = false;
    for (int i=0; i<3; i++){__grad_u[0][i] = _grad_u1[_qp](i); if(1e-6<fabs(__grad_u[0][i])){print_grad_u=true;}}
    for (int i=0; i<3; i++){__grad_u[1][i] = _grad_u2[_qp](i); if(1e-6<fabs(__grad_u[1][i])){print_grad_u=true;}}
    for (int i=0; i<3; i++){__grad_u[2][i] = _grad_u3[_qp](i); if(1e-6<fabs(__grad_u[2][i])){print_grad_u=true;}}

    if(print_grad_u){
    std::cout << "__grad_u:\n";
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            std::cout << __grad_u[i][j] << " ";
        }
        std::cout << "\n";
    }
    }
    
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

    //std::cout << "__phi: ";
    //for (int i=0; i<9; i++){std::cout << __phi[i] << " ";}
    //std::cout << "\n";

    //Copy over grad_phi
    for (int i=0; i<3; i++){__grad_phi[0][i] = _grad_phi_11[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[1][i] = _grad_phi_22[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[2][i] = _grad_phi_33[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[3][i] = _grad_phi_23[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[4][i] = _grad_phi_13[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[5][i] = _grad_phi_12[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[6][i] = _grad_phi_32[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[7][i] = _grad_phi_31[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[8][i] = _grad_phi_21[_qp](i);}
 
    //std::cout << "__grad_phi:\n";
    //for (int i=0; i<9; i++){
    //    for (int j=0; j<3; j++){
    //        std::cout << __grad_phi[i][j] << " ";
    //    }
    //    std::cout << "\n";
    //}

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

    std::vector<std::vector<std::vector<double>>> ADD_JACOBIANS;
    ADD_JACOBIANS.resize(_n_ADD_JACOBIANS);
    //END of hardcoded values

    //Set the sizes of the additional term vectors
    _ADD_TERMS[_qp].resize(_n_ADD_TERMS);

    //Evaluate the model
    material->evaluate_model(time, _fparams, __grad_u, __phi, __grad_phi,
                              SDVS,                ADD_DOF,           ADD_grad_DOF,
                             _cauchy[_qp],         _s[_qp],           _m[_qp],
                             _DcauchyDgrad_u[_qp], _DcauchyDphi[_qp], _DcauchyDgrad_phi[_qp],
                             _DsDgrad_u[_qp],      _DsDphi[_qp],      _DsDgrad_phi[_qp],
                             _DmDgrad_u[_qp],      _DmDphi[_qp],      _DmDgrad_phi[_qp],
                             _ADD_TERMS[_qp],      ADD_JACOBIANS);

    
}

