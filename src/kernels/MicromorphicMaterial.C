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

    params.set<bool>("use_displaced_mesh") = true; //TODO: Note that this is hard-coded such that we will always be using the
                                                   //      gradient w.r.t. the current coordinates. We might want to add a 
                                                   //      more general approach which allows the user to implement a 
                                                   //      total-Lagrangian formulation if desired. This would make the 
                                                   //      implementation faster but potentially could add confusion.

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

double relative_norm(const double &v1, const double &v2){
    /*!
    =======================
    |    relative_norm    |
    =======================

    Compute the relative norm between two values.
    */

    double tmp;
    tmp = std::max(fabs(v1),fabs(v2));
    if(tmp<1e-7){tmp=1;}

    return (v1-v2)/tmp;
    
}

void print_vector(const std::vector<double> &v){
    /*!
    ======================
    |    print_vector    |
    ======================

    Print a std::vector to the terminal.

    */

    for (int i=0; i<v.size(); i++){
        std::cout << v[i] << " ";
    }
    std::cout << std::endl;
    return;
}

void print_matrix(const std::vector<std::vector<double>> &M){
    /*!
    ======================
    |    print_matrix    |
    ======================

    Print a std::vector matrix to the terminal.

    */

    for (int i=0; i<M.size(); i++){
        print_vector(M[i]);
    }
    return;
}

void compare_vectors(const std::vector<double> &V1, const std::vector<double> &V2, double &relative_error){
    /*!
    =========================
    |    compare_vectors    |
    =========================

    Compare two std::vectors to determine the relative error between them.
    */

    if (V1.size() != V2.size()){mooseError("Vectors of different sizes cannot be compared!");}
    for (int i=0; i<V1.size(); i++){
        std::cout << relative_norm(V1[i],V2[i]) << " ";
        relative_error += fabs(relative_norm(V1[i],V2[i]));
    }
    std::cout << "\n";
    return;
}

void compare_matrices(const std::vector<std::vector<double>> &M1, const std::vector<std::vector<double>> &M2){
    /*!
    ==========================
    |    compare_matrices    |
    ==========================

    Compare two matrices to determine the relative error between them.

    */

    double tmp = 0;

    if (M1.size() != M2.size()){mooseError("Matricies of different sizes cannot be compared!");}
    for (int i=0; i<M1.size(); i++){compare_vectors(M1[i],M2[i],tmp);}
    std::cout << "relative error: " << tmp << "\n";
    if (tmp>1e-3){
                  std::cout << "###################################################"
                            << "#    WARNING: POOR QUALITY JACOBIAN DISCOVERED    #"
                            << "###################################################\n";
                  mooseError("Poor quality stress jacobian discovered");
                 }
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

    //std::cout << "coords: " << _q_point[_qp](0) << " " << _q_point[_qp](1) << " " << _q_point[_qp](2) << "\n";
    //std::cout << "u:      " << _u1[_qp] << " " << _u2[_qp] << " " << _u3[_qp] << "\n";

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

    //Extract the time
    std::vector<double> time;
    time.resize(2);

    time[0] = _t;
    time[1] = _dt;

    //TODO: Remove these hardcoded values
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

    //Jacobian debugging routines to follow uncomment only if absolutely necessary!
/*    std::vector<double> __cauchy;
    std::vector<double> __s;
    std::vector<double> __m;
    std::vector<std::vector<double>> __DcauchyDgrad_u;
    std::vector<std::vector<double>> __DcauchyDphi;
    std::vector<std::vector<double>> __DcauchyDgrad_phi;
    std::vector<std::vector<double>> __DsDgrad_u;
    std::vector<std::vector<double>> __DsDphi;
    std::vector<std::vector<double>> __DsDgrad_phi;
    std::vector<std::vector<double>> __DmDgrad_u;
    std::vector<std::vector<double>> __DmDphi;
    std::vector<std::vector<double>> __DmDgrad_phi;
    std::vector<std::vector<double>> __ADD_TERMS;
    std::vector<std::vector<std::vector<double>>> __ADD_JACOBIANS;

    material->evaluate_model_numeric_gradients(time, _fparams, __grad_u, __phi, __grad_phi,
                              SDVS,                ADD_DOF,           ADD_grad_DOF,
                             __cauchy,         __s,           __m,
                             __DcauchyDgrad_u, __DcauchyDphi, __DcauchyDgrad_phi,
                             __DsDgrad_u,      __DsDphi,      __DsDgrad_phi,
                             __DmDgrad_u,      __DmDphi,      __DmDgrad_phi,
                             __ADD_TERMS,      __ADD_JACOBIANS);

    //Check the returned values
    std::cout << "_qp:          " << _qp << "\n";
    std::cout << "__grad_u:\n";
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            std::cout << __grad_u[i][j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "Differences in stresses (must be zero or something is really wrong!)\n";
    std::cout << "cauchy error: ";
    for (int i=0; i<9; i++){std::cout << _cauchy[_qp][i] - __cauchy[i] << " ";}
    std::cout << "\n";
    std::cout << "s error: ";
    for (int i=0; i<9; i++){std::cout << _s[_qp][i] - __s[i] << " ";}
    std::cout << "\n";
    std::cout << "m error: ";
    for (int i=0; i<27; i++){std::cout << _m[_qp][i] - __m[i] << " ";}
    std::cout << "\n";

    std::cout << "Differences in stress gradients (should be close to zero!)\n";

    std::cout << "DcauchyDgrad_u analytic:\n";
    print_matrix(_DcauchyDgrad_u[_qp]);
    std::cout << "DcauchyDgrad_u numeric:\n";
    print_matrix(__DcauchyDgrad_u);
    std::cout << "_DcauchyDgrad_u error:\n";
    compare_matrices(_DcauchyDgrad_u[_qp],__DcauchyDgrad_u);

    std::cout << "DcauchyDphi analytic:\n";
    print_matrix(_DcauchyDphi[_qp]);
    std::cout << "DcauchyDphi numeric:\n";
    print_matrix(__DcauchyDphi);
    std::cout << "_DcauchyDphi error:\n";
    compare_matrices(_DcauchyDphi[_qp],__DcauchyDphi);

    std::cout << "DcauchyDgrad_phi analytic:\n";
    print_matrix(_DcauchyDgrad_phi[_qp]);
    std::cout << "DcauchyDgrad_phi numeric:\n";
    print_matrix(__DcauchyDgrad_phi);
    std::cout << "_DcauchyDgrad_phi error:\n";
    compare_matrices(_DcauchyDgrad_phi[_qp],__DcauchyDgrad_phi);

    std::cout << "DsDgrad_u analytic:\n";
    print_matrix(_DsDgrad_u[_qp]);
    std::cout << "DsDgrad_u numeric:\n";
    print_matrix(__DsDgrad_u);
    std::cout << "_DsDgrad_u error:\n";
    compare_matrices(_DsDgrad_u[_qp],__DsDgrad_u);

    std::cout << "DsDphi analytic:\n";
    print_matrix(_DsDphi[_qp]);
    std::cout << "DsDphi numeric:\n";
    print_matrix(__DsDphi);
    std::cout << "_DsDphi error:\n";
    compare_matrices(_DsDphi[_qp],__DsDphi);

    std::cout << "DsDgrad_phi analytic:\n";
    print_matrix(_DsDgrad_phi[_qp]);
    std::cout << "DsDgrad_phi numeric:\n";
    print_matrix(__DsDgrad_phi);
    std::cout << "_DsDgrad_phi error:\n";
    compare_matrices(_DsDgrad_phi[_qp],__DsDgrad_phi);

    std::cout << "DmDgrad_u analytic:\n";
    print_matrix(_DmDgrad_u[_qp]);
    std::cout << "DmDgrad_u numeric:\n";
    print_matrix(__DmDgrad_u);
    std::cout << "_DmDgrad_u error:\n";
    compare_matrices(_DmDgrad_u[_qp],__DmDgrad_u);

    std::cout << "DmDphi analytic:\n";
    print_matrix(_DmDphi[_qp]);
    std::cout << "DmDphi numeric:\n";
    print_matrix(__DmDphi);
    std::cout << "_DmDphi error:\n";
    compare_matrices(_DmDphi[_qp],__DmDphi);

    std::cout << "DmDgrad_phi analytic:\n";
    print_matrix(_DmDgrad_phi[_qp]);
    std::cout << "DmDgrad_phi numeric:\n";
    print_matrix(__DmDgrad_phi);
    std::cout << "_DmDgrad_phi error:\n";
    compare_matrices(_DmDgrad_phi[_qp],__DmDgrad_phi);

    if(_qp>=7){    
        mooseError("Expected error!");
    }
*/
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
/*        std::cout << "_cauchy_MMS[_qp]: ";
        for (int i=0; i<9; i++){std::cout << _cauchy_MMS[_qp][i] << " ";}
        std::cout << "\n";
        std::cout << "_s_MMS[_qp]: ";
        for (int i=0; i<9; i++){std::cout << _s_MMS[_qp][i] << " ";}
        std::cout << "\n";
        std::cout << "_m_MMS[_qp]: ";
        for (int i=0; i<27; i++){std::cout << _m_MMS[_qp][i] << " ";}
        std::cout << "\n";*/
    }

}

