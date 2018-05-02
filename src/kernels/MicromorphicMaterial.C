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

//MOOSE includes
#include "Function.h"

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

    params.addParam<bool>(
        "MMS", false, "Flag for whether to compute the method of manufactured solutions");

    // Coupled variables
    params.addRequiredCoupledVar(
        "u1", "The displacement in the 1 direction.");

    params.addRequiredCoupledVar(
        "u2", "The displacement in the 2 direction.");

    params.addRequiredCoupledVar(
        "u3", "The displacement in the 3 direction.");

    params.addRequiredCoupledVar(
        "phi_11", "The 11 component of the phi tensor.");

    params.addRequiredCoupledVar(
        "phi_22", "The 22 component of the phi tensor.");

    params.addRequiredCoupledVar(
        "phi_33", "The 33 component of the phi tensor.");

    params.addRequiredCoupledVar(
        "phi_23", "The 23 component of the phi tensor.");

    params.addRequiredCoupledVar(
        "phi_13", "The 13 component of the phi tensor.");

    params.addRequiredCoupledVar(
        "phi_12", "The 12 component of the phi tensor.");

    params.addRequiredCoupledVar(
        "phi_32", "The 32 component of the phi tensor.");

    params.addRequiredCoupledVar(
        "phi_31", "The 31 component of the phi tensor.");

    params.addRequiredCoupledVar(
        "phi_21", "The 21 component of the phi tensor.");

    // Functions for method of manufactured solutions
    params.addParam<FunctionName>(
        "u1_fxn", "0", "The function for displacement in the 1 direction.");

    params.addParam<FunctionName>(
        "u2_fxn", "0", "The function for displacement in the 2 direction.");

    params.addParam<FunctionName>(
        "u3_fxn", "0", "The function for displacement in the 3 direction.");

    params.addParam<FunctionName>(
        "phi_11_fxn", "0", "The function for the 11 component of the phi tensor.");

    params.addParam<FunctionName>(
        "phi_22_fxn", "0", "The function for the 22 component of the phi tensor.");

    params.addParam<FunctionName>(
        "phi_33_fxn", "0", "The function for the 33 component of the phi tensor.");

    params.addParam<FunctionName>(
        "phi_23_fxn", "0", "The function for the 23 component of the phi tensor.");

    params.addParam<FunctionName>(
        "phi_13_fxn", "0", "The function for the 13 component of the phi tensor.");

    params.addParam<FunctionName>(
        "phi_12_fxn", "0", "The function for the 12 component of the phi tensor.");

    params.addParam<FunctionName>(
        "phi_32_fxn", "0", "The function for the 32 component of the phi tensor.");

    params.addParam<FunctionName>(
        "phi_31_fxn", "0", "The function for the 31 component of the phi tensor.");

    params.addParam<FunctionName>(
        "phi_21_fxn", "0", "The function for the 21 component of the phi tensor.");

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
    _MMS(getParam<bool>("MMS")),
    _u1(coupledValue("u1")),
    _u2(coupledValue("u2")),
    _u3(coupledValue("u3")),
    _grad_u1(coupledGradient("u1")),
    _grad_u2(coupledGradient("u2")),
    _grad_u3(coupledGradient("u3")),
    _phi_11(coupledValue("phi_11")),
    _phi_22(coupledValue("phi_22")),
    _phi_33(coupledValue("phi_33")),
    _phi_23(coupledValue("phi_23")),
    _phi_13(coupledValue("phi_13")),
    _phi_12(coupledValue("phi_12")),
    _phi_32(coupledValue("phi_32")),
    _phi_31(coupledValue("phi_31")),
    _phi_21(coupledValue("phi_21")),
    _grad_phi_11(coupledGradient("phi_11")),
    _grad_phi_22(coupledGradient("phi_22")),
    _grad_phi_33(coupledGradient("phi_33")),
    _grad_phi_23(coupledGradient("phi_23")),
    _grad_phi_13(coupledGradient("phi_13")),
    _grad_phi_12(coupledGradient("phi_12")),
    _grad_phi_32(coupledGradient("phi_32")),
    _grad_phi_31(coupledGradient("phi_31")),
    _grad_phi_21(coupledGradient("phi_21")),
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
    _ADD_TERMS(declareProperty<std::vector<std::vector<double>>>("ADD_TERMS")),
    _ADD_JACOBIANS(declareProperty<std::vector<std::vector<std::vector<double>>>>("ADD_JACOBIANS")),
    _u1_fxn(getFunction("u1_fxn")),
    _u2_fxn(getFunction("u2_fxn")),
    _u3_fxn(getFunction("u3_fxn")),
    _phi_11_fxn(getFunction("phi_11_fxn")),
    _phi_22_fxn(getFunction("phi_22_fxn")),
    _phi_33_fxn(getFunction("phi_33_fxn")),
    _phi_23_fxn(getFunction("phi_23_fxn")),
    _phi_13_fxn(getFunction("phi_13_fxn")),
    _phi_12_fxn(getFunction("phi_12_fxn")),
    _phi_32_fxn(getFunction("phi_32_fxn")),
    _phi_31_fxn(getFunction("phi_31_fxn")),
    _phi_21_fxn(getFunction("phi_21_fxn")),
    _cauchy_MMS(declareProperty<std::vector<double>>("cauchy_MMS")),
    _s_MMS(declareProperty<std::vector<double>>("s_MMS")),
    _m_MMS(declareProperty<std::vector<double>>("m_MMS")),
    _ADD_TERMS_MMS(declareProperty<std::vector<std::vector<double>>>("ADD_TERMS_MMS"))

{
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

    RealVectorValue tmp_grad;

    //Copy over the gradient of u
    for (int i=0; i<3; i++){__grad_u[0][i] = _grad_u1[_qp](i);}
    for (int i=0; i<3; i++){__grad_u[1][i] = _grad_u2[_qp](i);}
    for (int i=0; i<3; i++){__grad_u[2][i] = _grad_u3[_qp](i);}
    
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
    for (int i=0; i<3; i++){__grad_phi[0][i] = _grad_phi_11[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[1][i] = _grad_phi_22[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[2][i] = _grad_phi_33[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[3][i] = _grad_phi_23[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[4][i] = _grad_phi_13[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[5][i] = _grad_phi_12[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[6][i] = _grad_phi_32[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[7][i] = _grad_phi_31[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[8][i] = _grad_phi_21[_qp](i);}

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

    //Set the sizes of the additional term vectors
    _ADD_TERMS[_qp].resize(_n_ADD_TERMS);
    _ADD_JACOBIANS[_qp].resize(_n_ADD_JACOBIANS);

    //Evaluate the model
    material->evaluate_model(time, _fparams, __grad_u, __phi, __grad_phi,
                              SDVS,                ADD_DOF,           ADD_grad_DOF,
                             _cauchy[_qp],         _s[_qp],           _m[_qp],
                             _DcauchyDgrad_u[_qp], _DcauchyDphi[_qp], _DcauchyDgrad_phi[_qp],
                             _DsDgrad_u[_qp],      _DsDphi[_qp],      _DsDgrad_phi[_qp],
                             _DmDgrad_u[_qp],      _DmDphi[_qp],      _DmDgrad_phi[_qp],
                             _ADD_TERMS[_qp],      _ADD_JACOBIANS[_qp]);

    //Evaluate method of manufactured solutions stresses
    if(_MMS){

        //Compute and copy the gradient of u
        tmp_grad = _u1_fxn.gradient(_t,_q_point[_qp]);
        for (int i=0; i<3; i++){__grad_u[0][i] = tmp_grad(i);}
        tmp_grad = _u2_fxn.gradient(_t,_q_point[_qp]);
        for (int i=0; i<3; i++){__grad_u[1][i] = tmp_grad(i);}
        tmp_grad = _u3_fxn.gradient(_t,_q_point[_qp]);
        for (int i=0; i<3; i++){__grad_u[2][i] = tmp_grad(i);}

        //Compute and copy phi
        __phi[0] = _phi_11_fxn.value(_t,_q_point[_qp]);
        __phi[1] = _phi_22_fxn.value(_t,_q_point[_qp]);
        __phi[2] = _phi_33_fxn.value(_t,_q_point[_qp]);
        __phi[3] = _phi_23_fxn.value(_t,_q_point[_qp]);
        __phi[4] = _phi_13_fxn.value(_t,_q_point[_qp]);
        __phi[5] = _phi_12_fxn.value(_t,_q_point[_qp]);
        __phi[6] = _phi_32_fxn.value(_t,_q_point[_qp]);
        __phi[7] = _phi_31_fxn.value(_t,_q_point[_qp]);
        __phi[8] = _phi_21_fxn.value(_t,_q_point[_qp]);

        //Compute and copy grad_phi
        tmp_grad = _phi_11_fxn.gradient(_t,_q_point[_qp]);
        for (int i=0; i<3; i++){__grad_phi[0][i] = tmp_grad(i);}
        tmp_grad = _phi_22_fxn.gradient(_t,_q_point[_qp]);
        for (int i=0; i<3; i++){__grad_phi[1][i] = tmp_grad(i);}
        tmp_grad = _phi_33_fxn.gradient(_t,_q_point[_qp]);
        for (int i=0; i<3; i++){__grad_phi[2][i] = tmp_grad(i);}
        tmp_grad = _phi_23_fxn.gradient(_t,_q_point[_qp]);
        for (int i=0; i<3; i++){__grad_phi[3][i] = tmp_grad(i);}
        tmp_grad = _phi_13_fxn.gradient(_t,_q_point[_qp]);
        for (int i=0; i<3; i++){__grad_phi[4][i] = tmp_grad(i);}
        tmp_grad = _phi_12_fxn.gradient(_t,_q_point[_qp]);
        for (int i=0; i<3; i++){__grad_phi[5][i] = tmp_grad(i);}
        tmp_grad = _phi_32_fxn.gradient(_t,_q_point[_qp]);
        for (int i=0; i<3; i++){__grad_phi[6][i] = tmp_grad(i);}
        tmp_grad = _phi_31_fxn.gradient(_t,_q_point[_qp]);
        for (int i=0; i<3; i++){__grad_phi[7][i] = tmp_grad(i);}
        tmp_grad = _phi_21_fxn.gradient(_t,_q_point[_qp]);
        for (int i=0; i<3; i++){__grad_phi[8][i] = tmp_grad(i);}

        //TODO: Add in function support for the additional DOF and their gradients.        

        //Evaluate the method of manufactured solutions stresses
        material->evaluate_model(time, _fparams, __grad_u, __phi, __grad_phi,
                                 SDVS,             ADD_DOF,     ADD_grad_DOF,
                                 _cauchy_MMS[_qp], _s_MMS[_qp], _m_MMS[_qp], _ADD_TERMS_MMS[_qp]);

    }

}

