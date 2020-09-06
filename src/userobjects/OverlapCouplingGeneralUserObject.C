/*=============================================================================
 *                     OverlapCouplingGeneralUserObject                       *
 *============================================================================*
 * Communicate with the overlap coupling libraries outside of MOOSE to        *
 * perform the overlap coupling. The general schme is the following:          *
 *   - Prior to first increment                                               *
 *     - OverlapCouplingNodeUserObject: Initialize overlap variables to empty *
 *                                      ( both macro and micro )              *
 *   - Timestep Begins                                                        *
 *     - OverlapCouplingNodeUserObject: Initialize overlap variables to       *
 *                                      empty ( both macro and micro )        *
 *     - Nonlinear iteration                                                  *
 *       - Linear iteration( s )                                              *
 *         - All balance equation kernels update                              *
 *         - CoupledPenaltyDirichletBC: Apply nodal BC values from the        *
 *                                      OverlapCouplingNodeUserObject         *
 *                                      ( both macro and micro )              *
 *                                                                            *
 *       - OverlapCouplingNodeUserObject: Run the external overlap coupling   *
 *                                        method ( macro only ) and update    *
 *                                        the degree of freedom maps ( both   *
 *                                        macro and micro )                   *
 *============================================================================*/

#include "OverlapCouplingGeneralUserObject.h"

#include "overlapCoupling.h"

registerMooseObject( "tardigradeApp", OverlapCouplingGeneralUserObject );

template<>
InputParameters
validParams< OverlapCouplingGeneralUserObject >( ){
    InputParameters params = validParams< GeneralUserObject >( );
    params.addRequiredParam< std::string >( "overlap_configuration_filename", "The overlap configuration YAML filename" );
    params.addRequiredParam< bool > ( "is_macro", "Flag which indicates if this is the macro-scale ( true, micromorphic ) or the micro scale ( false, classical continuum )" );
    return params;

}

OverlapCouplingGeneralUserObject::OverlapCouplingGeneralUserObject( const InputParameters &parameters )
    : GeneralUserObject( parameters ),
        _overlap_coupling_filename( getParam< std::string >( "overlap_configuration_filename" ) ),
        _isMacro( getParam< bool >( "is_macro" ) ){
}

void OverlapCouplingGeneralUserObject::initialize( ){
    /*!
     * Initialize the overlap coupling general user object
     */

    std::cout << "initialize\n";
    _dim = _fe_problem.mesh( ).dimension( );

    return;

}

void OverlapCouplingGeneralUserObject::execute( ){
    /*!
     * Execute the general user object
     */

    std::cout << "execute\n";

    return;

}

void OverlapCouplingGeneralUserObject::threadJoin( const UserObject &y ){
    /*!
     * Join the threads
     */

    std::cout << "thread join\n";

    return;

}

void OverlapCouplingGeneralUserObject::finalize( ){
    /*!
     * Finalize the user object
     */

    std::cout << "finalizing\n";

    return;

}

bool OverlapCouplingGeneralUserObject::getNodalValue( const dof_id_type nodeId, const dof_id_type dofIndex, Real &nodalDOFValue ){
    /*!
     * Get the nodal value for the given DOF
     *
     * :param const dof_id_type node_id: The node ID
     * :param const dof_id_type dofIndex: The index of the DOF
     *     ( 0 - 2: Displacements, 3-11: Micro-displacements )
     * :param Real &nodalDOFValue: The value of the degree of freedom at the node
     */

    auto it = std::find( _nodeIds.begin( ), _nodeIds.end( ),  nodeId );

    if ( it == _nodeIds.end( ) ){

        return false;

    }
    else{

        dof_id_type index;
        if ( _isMacro ){

            index = ( _dim + _dim * _dim ) * ( it - _nodeIds.begin( ) ) + dofIndex;

        }
        else{

            index = _dim * ( it - _nodeIds.begin( ) );

        }

        if ( index >= _dofValues.size( ) ){

            mooseError( "The index for the node with id " + std::to_string( nodeId ) + " is too large for the DOF value vector" ); 

        }

        nodalDOFValue = _dofValues[ index ];

        return true;

    }

    return false;
}
