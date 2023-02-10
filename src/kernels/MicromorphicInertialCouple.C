/*!
=====================================================================
|                  MicromorphicInertialCouple.cpp                   |
---------------------------------------------------------------------
| The source file for the micromorphic inertial couple kernel. This |
| kernel draws upon the balance equations from the git repository   |
| micromorphic_element.                                             |
| and requires nine degrees of freedom to be defined for solid      |
| mechanics. These degrees of freedom are:                          |
|     phi_11, phi_12, phi_13, phi_21, phi_22, phi_23, phi_31,       |
|     phi_32, phi_33                                                |
| which are the components of the micro-displacement tensor.        |
=====================================================================
*/

#include<MicromorphicInertialCouple.h>

//We define the valid parameters for this kernel and their default values
registerMooseObject("tardigradeApp", MicromorphicInertialCouple);

InputParameters
MicromorphicInertialCouple::validParams()
{
    InputParameters params = Kernel::validParams();
    params.set< bool >( "use_displaced_mesh" ) = false;
    params.addRequiredParam< int >( "component_i", "The i component of the inertial couple tensor" );
    params.addRequiredParam< int >( "component_j", "The j component of the inertial couple tensor" );
    params.addRequiredParam< int >( "dof_num",   "The degree of freedom to use for the diagonal jacobian calculation" );
    params.addParam< bool >( "MMS", false,       "The flag for whether the run will be using the method of manufactured solutions" );
    params.addRequiredParam< Real >( "reference_density", "The density in the reference configuration" );
    params.addRequiredParam< Real >( "I11", "The 11 component of the reference mass moment of inertia" );
    params.addRequiredParam< Real >( "I12", "The 12 component of the reference mass moment of inertia" );
    params.addRequiredParam< Real >( "I13", "The 13 component of the reference mass moment of inertia" );
    params.addRequiredParam< Real >( "I22", "The 22 component of the reference mass moment of inertia" );
    params.addRequiredParam< Real >( "I23", "The 23 component of the reference mass moment of inertia" );
    params.addRequiredParam< Real >( "I33", "The 33 component of the reference mass moment of inertia" );
    params.addRequiredCoupledVar( "phi11", "Micro-displacement in the 11 direction" );
    params.addRequiredCoupledVar( "phi12", "Micro-displacement in the 12 direction" );
    params.addRequiredCoupledVar( "phi13", "Micro-displacement in the 13 direction" );
    params.addRequiredCoupledVar( "phi21", "Micro-displacement in the 21 direction" );
    params.addRequiredCoupledVar( "phi22", "Micro-displacement in the 22 direction" );
    params.addRequiredCoupledVar( "phi23", "Micro-displacement in the 23 direction" );
    params.addRequiredCoupledVar( "phi31", "Micro-displacement in the 31 direction" );
    params.addRequiredCoupledVar( "phi32", "Micro-displacement in the 32 direction" );
    params.addRequiredCoupledVar( "phi33", "Micro-displacement in the 33 direction" );
    return params;
}

MicromorphicInertialCouple::MicromorphicInertialCouple( const InputParameters & parameters )
    : // We have to call the constructor for the base class first
        Kernel( parameters ),
        _component_i( getParam< int >( "component_i" ) ),
        _component_j( getParam< int >( "component_j" ) ),
        _dof_num( getParam< int >( "dof_num" ) ),
        _MMS( getParam< bool >( "MMS" ) ),
        _phi11_int(isCoupled("phi11") ? coupled("phi11")
                                      : 100),
        _phi12_int(isCoupled("phi12") ? coupled("phi12")
                                      : 100),
        _phi13_int(isCoupled("phi13") ? coupled("phi13")
                                      : 100),
        _phi21_int(isCoupled("phi21") ? coupled("phi21")
                                      : 100),
        _phi22_int(isCoupled("phi22") ? coupled("phi22")
                                      : 100),
        _phi23_int(isCoupled("phi23") ? coupled("phi23")
                                      : 100),
        _phi31_int(isCoupled("phi31") ? coupled("phi31")
                                      : 100),
        _phi32_int(isCoupled("phi32") ? coupled("phi32")
                                      : 100),
        _phi33_int(isCoupled("phi33") ? coupled("phi33")
                                      : 100),
        _density( getParam< Real >( "reference_density" ) ),
        _I11( getParam< Real >( "I11" ) ),
        _I12( getParam< Real >( "I12" ) ),
        _I13( getParam< Real >( "I13" ) ),
        _I22( getParam< Real >( "I22" ) ),
        _I23( getParam< Real >( "I23" ) ),
        _I33( getParam< Real >( "I33" ) ),
        _phi11( coupledValue( "phi11" ) ),
        _phi12( coupledValue( "phi12" ) ),
        _phi13( coupledValue( "phi13" ) ),
        _phi21( coupledValue( "phi21" ) ),
        _phi22( coupledValue( "phi22" ) ),
        _phi23( coupledValue( "phi23" ) ),
        _phi31( coupledValue( "phi31" ) ),
        _phi32( coupledValue( "phi32" ) ),
        _phi33( coupledValue( "phi33" ) ),
        _dotDotChi11( coupledDotDot( "phi11" ) ),
        _dotDotChi12( coupledDotDot( "phi12" ) ),
        _dotDotChi13( coupledDotDot( "phi13" ) ),
        _dotDotChi21( coupledDotDot( "phi21" ) ),
        _dotDotChi22( coupledDotDot( "phi22" ) ),
        _dotDotChi23( coupledDotDot( "phi23" ) ),
        _dotDotChi31( coupledDotDot( "phi31" ) ),
        _dotDotChi32( coupledDotDot( "phi32" ) ),
        _dotDotChi33( coupledDotDot( "phi33" ) ),
        _dDotDotChi11Du( coupledDotDotDu( "phi11" ) ),
        _dDotDotChi12Du( coupledDotDotDu( "phi12" ) ),
        _dDotDotChi13Du( coupledDotDotDu( "phi13" ) ),
        _dDotDotChi21Du( coupledDotDotDu( "phi21" ) ),
        _dDotDotChi22Du( coupledDotDotDu( "phi22" ) ),
        _dDotDotChi23Du( coupledDotDotDu( "phi23" ) ),
        _dDotDotChi31Du( coupledDotDotDu( "phi31" ) ),
        _dDotDotChi32Du( coupledDotDotDu( "phi32" ) ),
        _dDotDotChi33Du( coupledDotDotDu( "phi33" ) )
    {
    /*!
    =====================
    |    Constructor    |
    =====================

    The constructor for the MicromorphicInertialCouple class.
    Note that this constructor is just variable 
    assignments.

    */

    if ( _MMS ){

        mooseError( "Method of manufactured solutions not implemented for MicromorphicInertialCouple" );

    }

}

Real MicromorphicInertialCouple::computeQpResidual(){
    /*!
    ===========================
    |    computeQpResidual    |
    ===========================

    Compute the residual at the quadrature point for 
    the indicated component.

    cint_i = -N \rho_0 I_{KL } \ddot{ \chi }_{kK} \chi_{lL}

    where i = _component

    */

    Real cinertia_ij;

    double _chi[ 9 ] = { 1 + _phi11[ _qp ],     _phi12[ _qp ],     _phi13[ _qp ],
                             _phi21[ _qp ], 1 + _phi22[ _qp ],     _phi23[ _qp ],
                             _phi31[ _qp ],     _phi32[ _qp ],     _phi33[ _qp ]  };

    double _dotDotChi[ 9 ] = { _dotDotChi11[ _qp ], _dotDotChi12[ _qp ], _dotDotChi13[ _qp ],
                               _dotDotChi21[ _qp ], _dotDotChi22[ _qp ], _dotDotChi23[ _qp ],
                               _dotDotChi31[ _qp ], _dotDotChi32[ _qp ], _dotDotChi33[ _qp ] };

    double _inertia[ 9 ] = { _I11, _I12, _I13,
                             _I12, _I22, _I23,
                             _I13, _I23, _I33 };
    
    int errorCode = balance_equations::compute_inertia_couple( _component_i, _component_j,
                                                               _test[ _i ][ _qp ],
                                                               _density, _chi, _dotDotChi, _inertia,
                                                               cinertia_ij );

    if ( errorCode != 0 ){
        mooseError( "Error code: " + std::to_string( errorCode ) + 
                    " returned from computation of the inertia couple for component " + 
                    std::to_string( _component_i ) + std::to_string( _component_j ) );
    }

    return cinertia_ij;
}

Real MicromorphicInertialCouple::computeQpJacobian(){
    /*!
    ===========================
    |    computeQpJacobian    |
    ===========================

    Compute the diagonal jacobian term.

    */
    Real dcinertiadU_ij;

    double _chi[ 9 ] = { 1 + _phi11[ _qp ],     _phi12[ _qp ],     _phi13[ _qp ],
                             _phi21[ _qp ], 1 + _phi22[ _qp ],     _phi23[ _qp ],
                             _phi31[ _qp ],     _phi32[ _qp ], 1 + _phi33[ _qp ]  };

    double _dotDotChi[ 9 ] = { _dotDotChi11[ _qp ], _dotDotChi12[ _qp ], _dotDotChi13[ _qp ],
                               _dotDotChi21[ _qp ], _dotDotChi22[ _qp ], _dotDotChi23[ _qp ],
                               _dotDotChi31[ _qp ], _dotDotChi32[ _qp ], _dotDotChi33[ _qp ] };

    std::vector< double > _dDotDotChiDu = { _dDotDotChi11Du[ _qp ], _dDotDotChi12Du[ _qp ], _dDotDotChi13Du[ _qp ],
                                            _dDotDotChi21Du[ _qp ], _dDotDotChi22Du[ _qp ], _dDotDotChi23Du[ _qp ],
                                            _dDotDotChi31Du[ _qp ], _dDotDotChi32Du[ _qp ], _dDotDotChi33Du[ _qp ] };

    double _inertia[ 9 ] = { _I11, _I12, _I13,
                             _I12, _I22, _I23,
                             _I13, _I23, _I33 };

    //Encode the indices to the ordering of the output
    int _I = 3 * _component_i + _component_j;
    int _J = _I;

    //Copy the test and interpolation functions so that the balance equation function can read it
    int errorCode = balance_equations::compute_inertia_couple_jacobian( _I, _J,
                                                                        _test[ _i ][ _qp ], _phi[ _j ][ _qp ],
                                                                        _density, _chi, _dotDotChi, _dDotDotChiDu,
                                                                        _inertia, dcinertiadU_ij );

    if ( errorCode != 0 ){
        mooseError( "Error code: " + std::to_string( errorCode ) + 
                    " returned from computation of the inertia couple jacobian for component " + 
                    std::to_string( _component_i ) + std::to_string( _component_j ) );
    }

    return dcinertiadU_ij;
}

Real MicromorphicInertialCouple::computeQpOffDiagJacobian(unsigned int jvar){
    /*!
    ==================================
    |    computeQpOffDiagJacobian    |
    ==================================

    Compute the off-diagonal terms of the jacobian
    */

    Real dcinertiadU_ij;

    double _chi[ 9 ] = { 1 + _phi11[ _qp ],     _phi12[ _qp ],     _phi13[ _qp ],
                             _phi21[ _qp ], 1 + _phi22[ _qp ],     _phi23[ _qp ],
                             _phi31[ _qp ],     _phi32[ _qp ], 1 + _phi33[ _qp ]  };

    double _dotDotChi[ 9 ] = { _dotDotChi11[ _qp ], _dotDotChi12[ _qp ], _dotDotChi13[ _qp ],
                               _dotDotChi21[ _qp ], _dotDotChi22[ _qp ], _dotDotChi23[ _qp ],
                               _dotDotChi31[ _qp ], _dotDotChi32[ _qp ], _dotDotChi33[ _qp ] };

    std::vector< double > _dDotDotChiDu = { _dDotDotChi11Du[ _qp ], _dDotDotChi12Du[ _qp ], _dDotDotChi13Du[ _qp ],
                                            _dDotDotChi21Du[ _qp ], _dDotDotChi22Du[ _qp ], _dDotDotChi23Du[ _qp ],
                                            _dDotDotChi31Du[ _qp ], _dDotDotChi32Du[ _qp ], _dDotDotChi33Du[ _qp ] };

    double _inertia[ 9 ] = { _I11, _I12, _I13,
                             _I12, _I22, _I23,
                             _I13, _I23, _I33 };

    //Encode the indices
    int  _I = 3 * _component_i + _component_j;
    int  _J = -1;

    if(jvar == _phi11_int){
        _J = 0;
    }
    else if(jvar == _phi12_int){
        _J = 1;
    }
    else if(jvar == _phi13_int){
        _J = 2;
    }
    else if(jvar == _phi21_int){
        _J = 3;
    }
    else if(jvar == _phi22_int){
        _J = 4;
    }
    else if(jvar == _phi23_int){
        _J = 5;
    }
    else if(jvar == _phi31_int){
        _J = 6;
    }
    else if(jvar == _phi32_int){
        _J = 7;
    }
    else if(jvar == _phi33_int){
        _J = 8;
    }

    //Copy the test and interpolation functions so that the balance equation function can read it
    if ( _J >= 0 ){

        int errorCode = balance_equations::compute_inertia_couple_jacobian( _I, _J,
                                                                            _test[ _i ][ _qp ], _phi[ _j ][ _qp ],
                                                                            _density, _chi, _dotDotChi, _dDotDotChiDu,
                                                                            _inertia, dcinertiadU_ij );

        if ( errorCode != 0 ){
            mooseError( "Error code: " + std::to_string( errorCode ) + 
                        " returned from computation of the inertia couple jacobian for component " + 
                        std::to_string( _component_i ) + std::to_string( _component_j ) );
        }

        return dcinertiadU_ij;

    }
    else{

        return 0;

    }

}
