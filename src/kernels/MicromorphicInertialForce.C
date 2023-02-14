/*!
=====================================================================
|                   MicromorphicInertialForce.cpp                   |
---------------------------------------------------------------------
| The source file for the micromorphic inertial force kernel. This  |
| kernel draws upon the balance equations from the git repository   |
| micromorphic_element                                              |
| and requires three degrees of freedom to be defined for solid     |
| mechanics. These degrees of freedom are:                          |
|     u_1, u_2, u_3                                                 |
| where u_1 -> u_3 are the macro deformations                       |
=====================================================================
*/

#include<MicromorphicInertialForce.h>

//We define the valid parameters for this kernel and their default values
registerMooseObject("tardigradeApp", MicromorphicInertialForce);

InputParameters
MicromorphicInertialForce::validParams()
{
    InputParameters params = Kernel::validParams();
    params.set< bool >( "use_displaced_mesh" ) = false;
    params.addRequiredParam< int >( "component", "The component of the inertial force vector" );
    params.addRequiredParam< int >( "dof_num",   "The degree of freedom to use for the diagonal jacobian calculation" );
    params.addParam< bool >( "MMS", false,       "The flag for whether the run will be using the method of manufactured solutions" );
    params.addRequiredParam< Real >( "reference_density", "The density in the reference configuration" );
    params.addRequiredCoupledVar( "u1", "Macro-displacement in the 1 direction" );
    params.addRequiredCoupledVar( "u2", "Macro-displacement in the 2 direction" );
    params.addRequiredCoupledVar( "u3", "Macro-displacement in the 3 direction" );
    return params;
}

MicromorphicInertialForce::MicromorphicInertialForce( const InputParameters & parameters )
    : // We have to call the constructor for the base class first
        Kernel( parameters ),
        _component( getParam< int >( "component" ) ),
        _dof_num( getParam< int >( "dof_num" ) ),
        _MMS( getParam< bool >( "MMS" ) ),
        _u1_int(isCoupled("u1") ? coupled("u1")
                                : 100),
        _u2_int(isCoupled("u2") ? coupled("u2")
                                : 100),
        _u3_int(isCoupled("u3") ? coupled("u3")
                                : 100),
        _density( getParam< Real >( "reference_density" ) ),
        _a1( coupledDotDot( "u1" ) ),
        _a2( coupledDotDot( "u2" ) ),
        _a3( coupledDotDot( "u3" ) ),
        _da1du( coupledDotDotDu( "u1" ) ),
        _da2du( coupledDotDotDu( "u2" ) ),
        _da3du( coupledDotDotDu( "u3" ) )
    {
    /*!
    =====================
    |    Constructor    |
    =====================

    The constructor for the MicromorphicInertialForce class.
    Note that this constructor is just variable 
    assignments.

    */

    if ( _MMS ){

        mooseError( "Method of manufactured solutions not implemented for MicromorphicInertialForce" );

    }

}

Real MicromorphicInertialForce::computeQpResidual(){
    /*!
    ===========================
    |    computeQpResidual    |
    ===========================

    Compute the residual at the quadrature point for 
    the indicated component.

    fint_i = -psi rho_0 a_i

    where i = _component

    To compute the tangent for our current residual, i.e., finertia,
    read the paragraph near "(see their definition)" at the link:
    https://www.dealii.org/current/doxygen/deal.II/step_8.html

    In other words, finertia implicitly contains a Kronecker delta
    function which evaluates to 1 (as _i shares the same physical
    component as _component, by definition).
    */

    Real finertia;

    double _acceleration[ 3 ] = { _a1[ _qp ], _a2[ _qp ], _a3[ _qp ] };
    
    //Copy the test function so that the balance equation function can read it
    int errorCode = balance_equations::compute_inertia_force( _component, _test[ _i ][ _qp ], _density, _acceleration, finertia );

    if ( errorCode != 0 ){
        mooseError( "Error code: " + std::to_string( errorCode ) + " returned from computation of the inertia force for component " + std::to_string( _component ) );
    }

    return finertia;
}

Real MicromorphicInertialForce::computeQpJacobian(){
    /*!
    ===========================
    |    computeQpJacobian    |
    ===========================

    Compute the block diagonal jacobian term.

    By "diagonal," dfinertiadU_ij implicitly contains two Kronecker 
    delta functions which both evaluate to 1, as _i and _j share the same 
    physical component as _component.
    */
    Real dfinertiadU_ij;

    double _acceleration[ 3 ] = { _a1[ _qp ], _a2[ _qp ], _a3[ _qp ] };
    std::vector< double > DaDu_i = { _da1du[ _qp ], _da2du[ _qp ], _da3du[ _qp ] };

    //Copy the test and interpolation functions so that the balance equation function can read it
    int errorCode = balance_equations::compute_inertia_force_jacobian( _component, _component,
                                                                       _test[ _i ][ _qp ], _phi[ _j ][ _qp ], _density,
                                                                       _acceleration, DaDu_i, dfinertiadU_ij ); 

    if ( errorCode != 0 ){
        mooseError( "Error code: " + std::to_string( errorCode ) +
                    " returned from computation of the inertia force jacobian for component " +
                    std::to_string( _component ) );
    }

    return dfinertiadU_ij;
}

Real MicromorphicInertialForce::computeQpOffDiagJacobian(unsigned int jvar){
    /*!
    ==================================
    |    computeQpOffDiagJacobian    |
    ==================================

    Compute the block off-diagonal terms of the jacobian

    By "off-diagonal," _j does not share the physical component of _i.
    Since current residual finertia only conjugates _test[ _i ] with the
    acceleration component of _i, all off-diagonal tangents are zero.
    */

    int  _off_diag_dof_num = -1;

    if(jvar == _u1_int){
        _off_diag_dof_num = 0;
    }
    else if(jvar == _u2_int){
        _off_diag_dof_num = 1;
    }
    else if(jvar == _u3_int){
        _off_diag_dof_num = 2;
    }

    Real dfinertiadU_ij;

    double _acceleration[ 3 ] = { _a1[ _qp ], _a2[ _qp ], _a3[ _qp ] };
    std::vector< double > DaDu_i = { _da1du[ _qp ], _da2du[ _qp ], _da3du[ _qp ] };

    //Copy the test and interpolation functions so that the balance equation function can read it
    if ( _off_diag_dof_num >= 0 ){

        int errorCode = balance_equations::compute_inertia_force_jacobian( _component, _off_diag_dof_num,
                                                                           _test[ _i ][ _qp ], _phi[ _j ][ _qp ], _density,
                                                                           _acceleration, DaDu_i, dfinertiadU_ij ); 
    
        if ( errorCode != 0 ){
            mooseError( "Error code: " + std::to_string( errorCode ) +
                        " returned from computation of the inertia force jacobian for component " +
                        std::to_string( _component ) );
        }

        return dfinertiadU_ij;

    }
    else{

        return 0;

    }

}
