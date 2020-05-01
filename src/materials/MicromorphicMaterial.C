/*!
====================================================================
|                     MicromorphicMaterial.cpp                     |
====================================================================
| The source file for a class which computes the PK2 stress and    |
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

    // The state variable array
    params.addParam< int >(
        "number_SDVS", 0, "The number of solution-dependent state variables" );

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

void MicromorphicMaterial::initQpStatefulProperties(){
    /*!
     * Initialize the internal state variable array
     */

    _SDVS[ _qp ] = std::vector< Real >(_n_SDVS, 0 );

    return;
}

MicromorphicMaterial::MicromorphicMaterial(const InputParameters & parameters)
    : Material(parameters),
    // Declare that this material is going to provide Eigen matrices containing the PK2 stress and 
    // jacobians that Kernels can use.
    _fparams(getParam<std::vector<Real>>("material_fparameters")),
    _n_ADD_DOF(getParam<int>("number_ADD_DOF")),
    _n_ADD_TERMS(getParam<int>("number_ADD_TERMS")),
    _n_ADD_JACOBIANS(getParam<int>("number_ADD_JACOBIANS")),
    _model_name(getParam<std::string>("model_name")),
    _MMS(getParam<bool>("MMS")),
    _n_SDVS( getParam< int >( "number_SDVS" ) ),
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
    _old_u1(coupledValueOld("u1")),
    _old_u2(coupledValueOld("u2")),
    _old_u3(coupledValueOld("u3")),
    _old_grad_u1(coupledGradientOld("u1")),
    _old_grad_u2(coupledGradientOld("u2")),
    _old_grad_u3(coupledGradientOld("u3")),
    _old_phi_11(coupledValueOld("phi_11")),
    _old_phi_22(coupledValueOld("phi_22")),
    _old_phi_33(coupledValueOld("phi_33")),
    _old_phi_23(coupledValueOld("phi_23")),
    _old_phi_13(coupledValueOld("phi_13")),
    _old_phi_12(coupledValueOld("phi_12")),
    _old_phi_32(coupledValueOld("phi_32")),
    _old_phi_31(coupledValueOld("phi_31")),
    _old_phi_21(coupledValueOld("phi_21")),
    _old_grad_phi_11(coupledGradientOld("phi_11")),
    _old_grad_phi_22(coupledGradientOld("phi_22")),
    _old_grad_phi_33(coupledGradientOld("phi_33")),
    _old_grad_phi_23(coupledGradientOld("phi_23")),
    _old_grad_phi_13(coupledGradientOld("phi_13")),
    _old_grad_phi_12(coupledGradientOld("phi_12")),
    _old_grad_phi_32(coupledGradientOld("phi_32")),
    _old_grad_phi_31(coupledGradientOld("phi_31")),
    _old_grad_phi_21(coupledGradientOld("phi_21")),
    _deformation_gradient(declareProperty<std::vector<double>>("MM_deformation_gradient")),
    _micro_displacement(declareProperty<std::vector<double>>("micro_displacement")),
    _gradient_micro_displacement(declareProperty<std::vector<std::vector<double>>>("gradient_micro_displacement")),
    _cauchy(declareProperty<std::vector<double>>("cauchy")),
    _s(declareProperty<std::vector<double>>("s")),
    _m(declareProperty<std::vector<double>>("m")),
    _PK2(declareProperty<std::vector<double>>("PK2")),
    _SIGMA(declareProperty<std::vector<double>>("SIGMA")),
    _M(declareProperty<std::vector<double>>("M")),
    _DPK2Dgrad_u(declareProperty<std::vector<std::vector<double>>>("DPK2Dgrad_u")),
    _DPK2Dphi(declareProperty<std::vector<std::vector<double>>>("DPK2Dphi")),
    _DPK2Dgrad_phi(declareProperty<std::vector<std::vector<double>>>("DPK2Dgrad_phi")),
    _DSIGMADgrad_u(declareProperty<std::vector<std::vector<double>>>("DSIGMADgrad_u")),
    _DSIGMADphi(declareProperty<std::vector<std::vector<double>>>("DSIGMADphi")),
    _DSIGMADgrad_phi(declareProperty<std::vector<std::vector<double>>>("DSIGMADgrad_phi")),
    _DMDgrad_u(declareProperty<std::vector<std::vector<double>>>("DMDgrad_u")),
    _DMDphi(declareProperty<std::vector<std::vector<double>>>("DMDphi")),
    _DMDgrad_phi(declareProperty<std::vector<std::vector<double>>>("DMDgrad_phi")),
    _ADD_TERMS(declareProperty<std::vector<std::vector<double>>>("ADD_TERMS")),
    _ADD_JACOBIANS(declareProperty<std::vector<std::vector<std::vector<double>>>>("ADD_JACOBIANS")),
    _SDVS( declareProperty< std::vector< double > > ( "SDVS" ) ),
    _old_SDVS( getMaterialPropertyOld< std::vector< double > >( "SDVS" ) ),
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
    _PK2_MMS(declareProperty<std::vector<double>>("PK2_MMS")),
    _SIGMA_MMS(declareProperty<std::vector<double>>("SIGMA_MMS")),
    _M_MMS(declareProperty<std::vector<double>>("M_MMS")),
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

    for (unsigned int i=0; i<v.size(); i++){
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

    for (unsigned int i=0; i<M.size(); i++){
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
    for (unsigned int i=0; i<V1.size(); i++){
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
    for (unsigned int i=0; i<M1.size(); i++){compare_vectors(M1[i],M2[i],tmp);}
    std::cout << "relative error: " << tmp << "\n";
    if (tmp>1e-3){
                  std::cout << "###################################################\n"
                            << "#    WARNING: POOR QUALITY JACOBIAN DISCOVERED    #\n"
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

    double __grad_u[ 3 ][ 3 ], __old_grad_u[ 3 ][ 3 ];
    double __phi[ 9 ], __old_phi[ 9 ];
    double __grad_phi[ 9 ][ 3 ], __old_grad_phi[ 9 ][ 3 ];

    RealVectorValue tmp_grad;

    //std::cout << "coords: " << _q_point[_qp](0) << " " << _q_point[_qp](1) << " " << _q_point[_qp](2) << "\n";
    //std::cout << "u:      " << _u1[_qp] << " " << _u2[_qp] << " " << _u3[_qp] << "\n";

    //Copy over the gradient of u
    for (int i=0; i<3; i++){__grad_u[0][i] = _grad_u1[_qp](i);}
    for (int i=0; i<3; i++){__grad_u[1][i] = _grad_u2[_qp](i);}
    for (int i=0; i<3; i++){__grad_u[2][i] = _grad_u3[_qp](i);}

    //Copy over the old gradient of u
    for (int i=0; i<3; i++){__old_grad_u[0][i] = _old_grad_u1[_qp](i);}
    for (int i=0; i<3; i++){__old_grad_u[1][i] = _old_grad_u2[_qp](i);}
    for (int i=0; i<3; i++){__old_grad_u[2][i] = _old_grad_u3[_qp](i);}
    
    //Copy over phi
    __phi[0] = _phi_11[_qp];
    __phi[1] = _phi_12[_qp];
    __phi[2] = _phi_13[_qp];
    __phi[3] = _phi_21[_qp];
    __phi[4] = _phi_22[_qp];
    __phi[5] = _phi_23[_qp];
    __phi[6] = _phi_31[_qp];
    __phi[7] = _phi_32[_qp];
    __phi[8] = _phi_33[_qp];

    //Copy over the old phi
    __old_phi[0] = _old_phi_11[_qp];
    __old_phi[1] = _old_phi_12[_qp];
    __old_phi[2] = _old_phi_13[_qp];
    __old_phi[3] = _old_phi_21[_qp];
    __old_phi[4] = _old_phi_22[_qp];
    __old_phi[5] = _old_phi_23[_qp];
    __old_phi[6] = _old_phi_31[_qp];
    __old_phi[7] = _old_phi_32[_qp];
    __old_phi[8] = _old_phi_33[_qp];

    //Copy over grad_phi
    for (int i=0; i<3; i++){__grad_phi[0][i] = _grad_phi_11[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[1][i] = _grad_phi_12[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[2][i] = _grad_phi_13[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[3][i] = _grad_phi_21[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[4][i] = _grad_phi_22[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[5][i] = _grad_phi_23[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[6][i] = _grad_phi_31[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[7][i] = _grad_phi_32[_qp](i);}
    for (int i=0; i<3; i++){__grad_phi[8][i] = _grad_phi_33[_qp](i);}

    //Copy over the old grad_phi
    for (int i=0; i<3; i++){__old_grad_phi[0][i] = _old_grad_phi_11[_qp](i);}
    for (int i=0; i<3; i++){__old_grad_phi[1][i] = _old_grad_phi_12[_qp](i);}
    for (int i=0; i<3; i++){__old_grad_phi[2][i] = _old_grad_phi_13[_qp](i);}
    for (int i=0; i<3; i++){__old_grad_phi[3][i] = _old_grad_phi_21[_qp](i);}
    for (int i=0; i<3; i++){__old_grad_phi[4][i] = _old_grad_phi_22[_qp](i);}
    for (int i=0; i<3; i++){__old_grad_phi[5][i] = _old_grad_phi_23[_qp](i);}
    for (int i=0; i<3; i++){__old_grad_phi[6][i] = _old_grad_phi_31[_qp](i);}
    for (int i=0; i<3; i++){__old_grad_phi[7][i] = _old_grad_phi_32[_qp](i);}
    for (int i=0; i<3; i++){__old_grad_phi[8][i] = _old_grad_phi_33[_qp](i);}

    //Store deformation values for use by other kernels

    //Copy the deformation gradient
    _deformation_gradient[_qp].resize(9);
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            _deformation_gradient[_qp][ 3 * i + j ] = __grad_u[i][j];
        }
    }
    _deformation_gradient[ _qp ][ 0 ] += 1;
    _deformation_gradient[ _qp ][ 4 ] += 1;
    _deformation_gradient[ _qp ][ 8 ] += 1;

    //Copy the micro-displacement
    _micro_displacement[_qp].resize(9);
    for ( int i=0; i<9; i++ ){
        _micro_displacement[ _qp ][ i ] = __phi[ i ];
    }
    _micro_displacement[ _qp ][ 0 ] += 1;
    _micro_displacement[ _qp ][ 4 ] += 1;
    _micro_displacement[ _qp ][ 8 ] += 1;

    //Copy the gradient of the micro-displacement
    _gradient_micro_displacement[_qp].resize(9);
    for (int i=0; i<9; i++){
        _gradient_micro_displacement[_qp][i].resize(3);
        for (int j=0; j<3; j++){
            _gradient_micro_displacement[_qp][i][j] = __grad_phi[i][j];
        }
    }

    //Evaluate the stresses and their jacobians
    auto &factory = micromorphic_material_library::MaterialFactory::Instance();
    auto material = factory.GetMaterial(_model_name);

    //Extract the time
    std::vector<double> time;
    time.resize(2);

    time[0] = _t;
    time[1] = _dt;

    //TODO: Remove these hardcoded values
    if ( _n_ADD_DOF > 0 ){
        mooseError( "MicromorphicMaterial does not support additional degrees of freedom" );
    }

    std::vector<double> ADD_DOF, old_ADD_DOF;
    ADD_DOF.resize( _n_ADD_DOF );
    old_ADD_DOF.resize( _n_ADD_DOF );

    std::vector<std::vector<double>> ADD_grad_DOF, old_ADD_grad_DOF;
    ADD_grad_DOF.resize( _n_ADD_DOF );
    old_ADD_grad_DOF.resize( _n_ADD_DOF );

    //END of hardcoded values

    //Set the sizes of the additional term vectors
    _ADD_TERMS[_qp].resize(_n_ADD_TERMS);
    _ADD_JACOBIANS[_qp].resize(_n_ADD_JACOBIANS);

    //Evaluate the model
    std::string output_message;

    int errorCode = material->evaluate_model( time, _fparams,
                                              __grad_u, __phi, __grad_phi,
                                              __old_grad_u, __old_phi, __old_grad_phi,
                                              _SDVS[ _qp ],
                                              ADD_DOF,            ADD_grad_DOF,
                                              old_ADD_DOF,        old_ADD_grad_DOF,
                                              _PK2[_qp],           _SIGMA[_qp],        _M[_qp],
                                              _DPK2Dgrad_u[_qp],   _DPK2Dphi[_qp],     _DPK2Dgrad_phi[_qp],
                                              _DSIGMADgrad_u[_qp], _DSIGMADphi[_qp],   _DSIGMADgrad_phi[_qp],
                                              _DMDgrad_u[_qp],     _DMDphi[_qp],       _DMDgrad_phi[_qp],
                                              _ADD_TERMS[_qp],     _ADD_JACOBIANS[_qp], output_message );

    if ( errorCode == 1 ){
        std::string error_message = "Convergence not achieved in material model. Requesting timestep cutback.\n";
        error_message += output_message;
        mooseException( error_message.c_str() );
    }

    if ( errorCode == 2 ){
        std::string error_message = "FATAL ERROR IN MICROMORPHIC MATERIAL MODEL\n";
        error_message += output_message;
        mooseError( error_message.c_str() );
    }

/*    //Jacobian debugging routines to follow uncomment only if absolutely necessary!
    std::vector<double> __PK2;
    std::vector<double> __SIGMA;
    std::vector<double> __M;
    std::vector<std::vector<double>> __DPK2Dgrad_u;
    std::vector<std::vector<double>> __DPK2Dphi;
    std::vector<std::vector<double>> __DPK2Dgrad_phi;
    std::vector<std::vector<double>> __DSIGMADgrad_u;
    std::vector<std::vector<double>> __DSIGMADphi;
    std::vector<std::vector<double>> __DSIGMADgrad_phi;
    std::vector<std::vector<double>> __DMDgrad_u;
    std::vector<std::vector<double>> __DMDphi;
    std::vector<std::vector<double>> __DMDgrad_phi;
    std::vector<std::vector<double>> __ADD_TERMS;
    std::vector<std::vector<std::vector<double>>> __ADD_JACOBIANS;

    material->evaluate_model_numeric_gradients(time, _fparams, __grad_u, __phi, __grad_phi,
                              SDVS,           ADD_DOF,       ADD_grad_DOF,
                             __PK2,           __SIGMA,       __M,
                             __DPK2Dgrad_u,   __DPK2Dphi,    __DPK2Dgrad_phi,
                             __DSIGMADgrad_u, __DSIGMADphi,  __DSIGMADgrad_phi,
                             __DMDgrad_u,     __DMDphi,      __DMDgrad_phi,
                             __ADD_TERMS,     __ADD_JACOBIANS, 1e-6);

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
    std::cout << "PK2 error: ";
    for (int i=0; i<9; i++){std::cout << _PK2[_qp][i] - __PK2[i] << " ";}
    std::cout << "\n";
    std::cout << "SIGMA error: ";
    for (int i=0; i<9; i++){std::cout << _SIGMA[_qp][i] - __SIGMA[i] << " ";}
    std::cout << "\n";
    std::cout << "M error: ";
    for (int i=0; i<27; i++){std::cout << _M[_qp][i] - __M[i] << " ";}
    std::cout << "\n";

    std::cout << "Differences in stress gradients (should be close to zero!)\n";

    std::cout << "DPK2Dgrad_u analytic:\n";
    print_matrix(_DPK2Dgrad_u[_qp]);
    std::cout << "DPK2Dgrad_u numeric:\n";
    print_matrix(__DPK2Dgrad_u);
    std::cout << "_DPK2Dgrad_u error:\n";
    compare_matrices(_DPK2Dgrad_u[_qp],__DPK2Dgrad_u);

    std::cout << "DPK2Dphi analytic:\n";
    print_matrix(_DPK2Dphi[_qp]);
    std::cout << "DPK2Dphi numeric:\n";
    print_matrix(__DPK2Dphi);
    std::cout << "_DPK2Dphi error:\n";
    compare_matrices(_DPK2Dphi[_qp],__DPK2Dphi);

    std::cout << "DPK2Dgrad_phi analytic:\n";
    print_matrix(_DPK2Dgrad_phi[_qp]);
    std::cout << "DPK2Dgrad_phi numeric:\n";
    print_matrix(__DPK2Dgrad_phi);
    std::cout << "_DPK2Dgrad_phi error:\n";
    compare_matrices(_DPK2Dgrad_phi[_qp],__DPK2Dgrad_phi);

    std::cout << "DSIGMADgrad_u analytic:\n";
    print_matrix(_DSIGMADgrad_u[_qp]);
    std::cout << "DSIGMADgrad_u numeric:\n";
    print_matrix(__DSIGMADgrad_u);
    std::cout << "_DSIGMADgrad_u error:\n";
    compare_matrices(_DSIGMADgrad_u[_qp],__DSIGMADgrad_u);

    std::cout << "DSIGMADphi analytic:\n";
    print_matrix(_DSIGMADphi[_qp]);
    std::cout << "DSIGMADphi numeric:\n";
    print_matrix(__DSIGMADphi);
    std::cout << "_DSIGMADphi error:\n";
    compare_matrices(_DSIGMADphi[_qp],__DSIGMADphi);

    std::cout << "DSIGMADgrad_phi analytic:\n";
    print_matrix(_DSIGMADgrad_phi[_qp]);
    std::cout << "DSIGMADgrad_phi numeric:\n";
    print_matrix(__DSIGMADgrad_phi);
    std::cout << "_DSIGMADgrad_phi error:\n";
    compare_matrices(_DSIGMADgrad_phi[_qp],__DSIGMADgrad_phi);

    std::cout << "DMDgrad_u analytic:\n";
    print_matrix(_DMDgrad_u[_qp]);
    std::cout << "DMDgrad_u numeric:\n";
    print_matrix(__DMDgrad_u);
    std::cout << "_DMDgrad_u error:\n";
    compare_matrices(_DMDgrad_u[_qp],__DMDgrad_u);

    std::cout << "DMDphi analytic:\n";
    print_matrix(_DMDphi[_qp]);
    std::cout << "DMDphi numeric:\n";
    print_matrix(__DMDphi);
    std::cout << "_DMDphi error:\n";
    compare_matrices(_DMDphi[_qp],__DMDphi);

    std::cout << "DMDgrad_phi analytic:\n";
    print_matrix(_DMDgrad_phi[_qp]);
    std::cout << "DMDgrad_phi numeric:\n";
    print_matrix(__DMDgrad_phi);
    std::cout << "_DMDgrad_phi error:\n";
    compare_matrices(_DMDgrad_phi[_qp],__DMDgrad_phi);

    if(_qp>=7){    
        mooseError("Expected error!");
    }
*/
    //Evaluate method of manufactured solutions stresses
    double mms_grad_u[ 3 ][ 3 ], mms_old_grad_u[ 3 ][ 3 ];
    double mms_phi[ 9 ], mms_old_phi[ 9 ];
    double mms_grad_phi[ 9 ][ 3 ], mms_old_grad_phi[ 9 ][ 3 ];
    if(_MMS){

        //compute and copy the gradient of u
        tmp_grad = _u1_fxn.gradient( _t, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++ ){
            mms_grad_u[ 0 ][ i ] = tmp_grad( i );
        }

        tmp_grad = _u2_fxn.gradient( _t, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++ ){
            mms_grad_u[ 1 ][ i ] = tmp_grad( i );
        }

        tmp_grad = _u3_fxn.gradient( _t, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++ ){
            mms_grad_u[ 2 ][ i ] = tmp_grad( i );
        }

        //compute and copy phi
        mms_phi[0] = _phi_11_fxn.value( _t, _q_point[ _qp ] );
        mms_phi[1] = _phi_12_fxn.value( _t, _q_point[ _qp ] );
        mms_phi[2] = _phi_13_fxn.value( _t, _q_point[ _qp ] );
        mms_phi[3] = _phi_21_fxn.value( _t, _q_point[ _qp ] );
        mms_phi[4] = _phi_22_fxn.value( _t, _q_point[ _qp ] );
        mms_phi[5] = _phi_23_fxn.value( _t, _q_point[ _qp ] );
        mms_phi[6] = _phi_31_fxn.value( _t, _q_point[ _qp ] );
        mms_phi[7] = _phi_32_fxn.value( _t, _q_point[ _qp ] );
        mms_phi[8] = _phi_33_fxn.value( _t, _q_point[ _qp ] );

        //compute and copy grad_phi
        tmp_grad = _phi_11_fxn.gradient( _t, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++){
            mms_grad_phi[ 0 ][ i ] = tmp_grad( i );
        }
        tmp_grad = _phi_12_fxn.gradient( _t, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++){
            mms_grad_phi[ 1 ][ i ] = tmp_grad( i );
        }
        tmp_grad = _phi_13_fxn.gradient( _t, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++){
            mms_grad_phi[ 2 ][ i ] = tmp_grad( i );
        }
        tmp_grad = _phi_21_fxn.gradient( _t, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++){
            mms_grad_phi[ 3 ][ i ] = tmp_grad( i );
        }
        tmp_grad = _phi_22_fxn.gradient( _t, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++){
            mms_grad_phi[ 4 ][ i ] = tmp_grad( i );
        }
        tmp_grad = _phi_23_fxn.gradient( _t, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++){
            mms_grad_phi[ 5 ][ i ] = tmp_grad( i );
        }
        tmp_grad = _phi_31_fxn.gradient( _t, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++){
            mms_grad_phi[ 6 ][ i ] = tmp_grad( i );
        }
        tmp_grad = _phi_32_fxn.gradient( _t, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++){
            mms_grad_phi[ 7 ][ i ] = tmp_grad( i );
        }
        tmp_grad = _phi_33_fxn.gradient( _t, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++){
            mms_grad_phi[ 8 ][ i ] = tmp_grad( i );
        }

        //Compute the old values
        //compute and copy the gradient of u
        tmp_grad = _u1_fxn.gradient( _t - _dt, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++ ){
            mms_old_grad_u[ 0 ][ i ] = tmp_grad( i );
        }

        tmp_grad = _u2_fxn.gradient( _t - _dt, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++ ){
            mms_old_grad_u[ 1 ][ i ] = tmp_grad( i );
        }

        tmp_grad = _u3_fxn.gradient( _t - _dt, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++ ){
            mms_old_grad_u[ 2 ][ i ] = tmp_grad( i );
        }

        //compute and copy phi
        mms_old_phi[0] = _phi_11_fxn.value( _t - _dt, _q_point[ _qp ] );
        mms_old_phi[1] = _phi_12_fxn.value( _t - _dt, _q_point[ _qp ] );
        mms_old_phi[2] = _phi_13_fxn.value( _t - _dt, _q_point[ _qp ] );
        mms_old_phi[3] = _phi_21_fxn.value( _t - _dt, _q_point[ _qp ] );
        mms_old_phi[4] = _phi_22_fxn.value( _t - _dt, _q_point[ _qp ] );
        mms_old_phi[5] = _phi_23_fxn.value( _t - _dt, _q_point[ _qp ] );
        mms_old_phi[6] = _phi_31_fxn.value( _t - _dt, _q_point[ _qp ] );
        mms_old_phi[7] = _phi_32_fxn.value( _t - _dt, _q_point[ _qp ] );
        mms_old_phi[8] = _phi_33_fxn.value( _t - _dt, _q_point[ _qp ] );

        //compute and copy grad_phi
        tmp_grad = _phi_11_fxn.gradient( _t - _dt, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++){
            mms_old_grad_phi[ 0 ][ i ] = tmp_grad( i );
        }
        tmp_grad = _phi_12_fxn.gradient( _t - _dt, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++){
            mms_old_grad_phi[ 1 ][ i ] = tmp_grad( i );
        }
        tmp_grad = _phi_13_fxn.gradient( _t - _dt, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++){
            mms_old_grad_phi[ 2 ][ i ] = tmp_grad( i );
        }
        tmp_grad = _phi_21_fxn.gradient( _t - _dt, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++){
            mms_old_grad_phi[ 3 ][ i ] = tmp_grad( i );
        }
        tmp_grad = _phi_22_fxn.gradient( _t - _dt, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++){
            mms_old_grad_phi[ 4 ][ i ] = tmp_grad( i );
        }
        tmp_grad = _phi_23_fxn.gradient( _t - _dt, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++){
            mms_old_grad_phi[ 5 ][ i ] = tmp_grad( i );
        }
        tmp_grad = _phi_31_fxn.gradient( _t - _dt, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++){
            mms_old_grad_phi[ 6 ][ i ] = tmp_grad( i );
        }
        tmp_grad = _phi_32_fxn.gradient( _t - _dt, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++){
            mms_old_grad_phi[ 7 ][ i ] = tmp_grad( i );
        }
        tmp_grad = _phi_33_fxn.gradient( _t - _dt, _q_point[ _qp ] );
        for ( int i = 0; i < 3; i++){
            mms_old_grad_phi[ 8 ][ i ] = tmp_grad( i );
        }

        //TODO: Add in function support for the additional DOF and their gradients.        

        //Evaluate the method of manufactured solutions stresses
        std::vector< double > mms_SDVS = _old_SDVS[ _qp ];
        errorCode = material->evaluate_model( time, _fparams,
                                              mms_grad_u, mms_phi, mms_grad_phi,
                                              mms_old_grad_u, mms_old_phi, mms_old_grad_phi,
                                              mms_SDVS,
                                              ADD_DOF,     ADD_grad_DOF,
                                              old_ADD_DOF, old_ADD_grad_DOF,
                                              _PK2_MMS[_qp], _SIGMA_MMS[_qp], _M_MMS[_qp], _ADD_TERMS_MMS[_qp], output_message );

        if ( errorCode == 1 ){
            std::string error_message = "Convergence not achieved in method of manufactured solutions evaluation";
            error_message += output_message;
            mooseException( error_message.c_str() );
        }
/*        std::cout << "_PK2_MMS[_qp]: ";
        for (int i=0; i<9; i++){std::cout << _PK2_MMS[_qp][i] << " ";}
        std::cout << "\n";
        std::cout << "_s_MMS[_qp]: ";
        for (int i=0; i<9; i++){std::cout << _s_MMS[_qp][i] << " ";}
        std::cout << "\n";
        std::cout << "_m_MMS[_qp]: ";
        for (int i=0; i<27; i++){std::cout << _m_MMS[_qp][i] << " ";}
        std::cout << "\n";*/
    }

}

