/*!
====================================================================
|                     GradientEnhancedDamagedMicromorphicMaterial.cpp                     |
====================================================================
| The source file for a class which computes the PK2 stress and    |
| the associated jacobians for a micromorphic material.            |
--------------------------------------------------------------------
| Notes: Relies on libraries from the micromorphic_element         |
|        repository available at bitbucket.org/NateAM2             |
====================================================================
*/


#include "GradientEnhancedDamagedMicromorphicMaterial.h"

//MOOSE includes
#include "Function.h"

registerMooseObject("tardigradeApp", GradientEnhancedDamagedMicromorphicMaterial);

template<>
InputParameters
validParams<GradientEnhancedDamagedMicromorphicMaterial>(){
    InputParameters params = validParams<Material>();

    // Vectors of material properties
    params.addRequiredParam<std::vector<Real>>(
        "material_fparameters", "The vector of floating point material parameters required for the stiffness matrices");

    params.addRequiredParam<std::string>(
        "model_name", "The material model name");

    params.addRequiredCoupledVar( "displacements", "The 3 displacement components" );
    params.addRequiredCoupledVar( "micro_displacement_gradient", "The 9 components of the micro displacement gradient" );
    /* params.addRequiredCoupledVar( "nonlocal_damage", "The nonlocal damage variable" ); */

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

void GradientEnhancedDamagedMicromorphicMaterial::initQpStatefulProperties(){
    /*!
     * Initialize the internal state variable array
     */

    _SDVS[ _qp ] = std::vector< Real >(_n_SDVS, 0 );

    return;
}

GradientEnhancedDamagedMicromorphicMaterial::GradientEnhancedDamagedMicromorphicMaterial(const InputParameters & parameters)
    : Material(parameters),
    // Declare that this material is going to provide Eigen matrices containing the PK2 stress and 
    // jacobians that Kernels can use.
    _fparams(getParam<std::vector<Real>>("material_fparameters")),
    /* _n_ADD_DOF(getParam<int>("number_ADD_DOF")), */
    /* _n_ADD_TERMS(getParam<int>("number_ADD_TERMS")), */
    /* _n_ADD_JACOBIANS(getParam<int>("number_ADD_JACOBIANS")), */
    _model_name(getParam<std::string>("model_name")),
    /* _MMS(getParam<bool>("MMS")), */
    _n_SDVS( getParam< int >( "number_SDVS" ) ),

    _grad_disp( coupledGradients( "displacements" ) ),
    _grad_disp_old( coupledGradientsOld( "displacements" ) ),

    _micro_disp_gradient( coupledValues( "micro_displacement_gradient" ) ),
    _micro_disp_gradient_old( coupledValuesOld( "micro_displacement_gradient" ) ),

    _grad_micro_disp_gradient( coupledGradients( "micro_displacement_gradient" ) ),
    _grad_micro_disp_gradient_old( coupledGradientsOld( "micro_displacement_gradient" ) ),

    /* _k( coupledValue( "nonlocal_damage" ) ), */

    _deformation_gradient(declareProperty<std::vector<double>>("MM_deformation_gradient")),
    _micro_deformation(declareProperty<std::vector<double>>("micro_deformation")),
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
    |    GradientEnhancedDamagedMicromorphicMaterial    |
    ==============================

    The constructor for GradientEnhancedDamagedMicromorphicMaterial.
    */

}

void GradientEnhancedDamagedMicromorphicMaterial::computeQpProperties(){
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

    //Copy over the gradient of u
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        __grad_u[i][j] = (*_grad_disp[i])[_qp](j);

    //Copy over the old gradient of u
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        __old_grad_u[i][j] = (*_grad_disp_old[i])[_qp](j);

    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
          __phi[i+j*3] = (*_micro_disp_gradient[i+j*3])[_qp];

    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
          __old_phi[i+j*3] = (*_micro_disp_gradient_old[i+j*3])[_qp];

    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        for (int k=0; k<3; k++)
          __grad_phi[i+j*3][k] = (*_grad_micro_disp_gradient[i+j*3])[_qp](k);

    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        for (int k=0; k<3; k++)
          __old_grad_phi[i+j*3][k] = (*_grad_micro_disp_gradient_old[i+j*3])[_qp](k);


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

    //Copy the micro-deformation
    _micro_deformation[_qp].resize(9);
    for ( int i=0; i<9; i++ ){
        _micro_deformation[ _qp ][ i ] = __phi[ i ];
    }
    _micro_deformation[ _qp ][ 0 ] += 1;
    _micro_deformation[ _qp ][ 4 ] += 1;
    _micro_deformation[ _qp ][ 8 ] += 1;

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

    /* //TODO: Remove these hardcoded values */
    /* if ( _n_ADD_DOF > 0 ){ */
    /*     mooseError( "GradientEnhancedDamagedMicromorphicMaterial does not support additional degrees of freedom" ); */
    /* } */

    const int _n_ADD_DOF = 0;
    const int _n_ADD_TERMS = 0;
    const int _n_ADD_JACOBIANS = 0;

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

#ifdef DEBUG_MODE
    std::map<std::string, std::map<std::string, std::map<std::string, std::vector<double>>>> debug;
#endif

    //Set the state variables to the previously converged values
    _SDVS[ _qp ] = _old_SDVS[ _qp ];

//    if ( _qp == 0 ){
//        std::cout << "SDVS pre model evaluation:\n";
//        vectorTools::print( _SDVS[ _qp ] );
//    }

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
                                              _ADD_TERMS[_qp],     _ADD_JACOBIANS[_qp], output_message
#ifdef DEBUG_MODE
                                              , debug
#endif  
                    );

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


    /*
     *
     * :param std::vector< double > &SDVS: The previously converged values of the state variables
     *     [ previousMacroStrainISV, previousMicroStrainISV, previousMicroGradientStrainISV,
     *       previousMacroGamma, previousMicroGamma, previousMicroGradientGamma,
     *       previousPlasticDeformationGradient - eye, previousPlasticMicroDeformation - eye,
     *       previousPlasticMicroGradient ]
     *
     * */
    const double* ptr_Fp = &_SDVS[_qp][7];
    const double* ptr_dFp_dF = &_ADD_JACOBIANS[_qp][0][0];
   


    //TODO: Add in function support for the additional DOF and their gradients.        

}

void computeGradientDamage(double* omega, double* dOmega_dPlasticDeformationGradient, double* plasticDeformationGradient, double* gradientDamageMaterialParameters, double* time, double* gradientDamageStateVars)
{


}

