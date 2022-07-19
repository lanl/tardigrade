#include<CouplingForce.h>

//Register the object
registerMooseObject( "tardigradeApp", CouplingForce );

InputParameters
CouplingForce::validParams( )
{
    InputParameters params = NodalKernel::validParams( );
    params.addRequiredParam< UserObjectName >( "overlap_coupling_object", "The overlap coupling user object" );
    params.addRequiredParam< unsigned int >( "component", "The component of the DOF vector being constrained (indexes from 0)" );
    params.addRequiredParam< bool > ( "is_macroscale", "Flag for whether the domain is the macroscale or the micro-scale" );
    return params;
}

CouplingForce::CouplingForce( const InputParameters &parameters )
    : // Call the constructor for the base class
        NodalKernel( parameters ),
        _overlap_coupling( getUserObject< OverlapCoupling >( "overlap_coupling_object" ) ),
        _component( getParam< unsigned int >( "component" ) ),
        _isMacro( getParam< bool >( "is_macroscale" ) ){

}

Real CouplingForce::computeQpResidual( ){
    /*!
     *
     * NOTE: IT IS ASSUMED THAT THE COUPLING FORCE IS APPLIED TO ALL ELEMENTS
     *       EVEN THOSE NOT BEING OVERLAPPED!
     */

    unsigned int _dim = _mesh.dimension( );

    auto elementCount = _mesh.nodeToElemMap( ).find( _current_node->id( ) );
    if ( elementCount == _mesh.nodeToElemMap( ).end( ) ){
        mooseError("Node " + std::to_string( _current_node->id( ) ) + " was not found in the node to element map");
    }
    if ( _isMacro ){

        auto localIndex = _overlap_coupling.getMacroGlobalLocalNodeMap( )->find( _current_node->id( ) );

        if ( localIndex == _overlap_coupling.getMacroGlobalLocalNodeMap( )->end( ) ){

            return 0;
            mooseError( "The current node is not found in the global-local node map" );

        }

        if ( _overlap_coupling.getMacroAugmentedLagrangianForce( )->size( ) <= ( ( _dim * _dim + _dim ) * localIndex->second + _component ) ){
            mooseError("macro oops");
        }

//        _console << "name: " << name( ) << " dim: " << _dim << " global node: " << localIndex->first << " local index: " << localIndex->second << " component: " << _component << " value: " << -_overlap_coupling.getMacroAugmentedLagrangianForce( )->at( ( _dim * _dim + _dim ) * localIndex->second + _component ) << "\n";

        return -_overlap_coupling.getMacroAugmentedLagrangianForce( )->at( ( _dim * _dim + _dim ) * localIndex->second + _component );

    }
    else{

        auto localIndex = _overlap_coupling.getMicroGlobalLocalNodeMap( )->find( _current_node->id( ) );

        if ( localIndex == _overlap_coupling.getMicroGlobalLocalNodeMap( )->end( ) ){

            return 0;
            mooseError( "The current node is not found in the global-local node map" );

        }

        if ( _overlap_coupling.getMicroAugmentedLagrangianForce( )->size( ) <= ( _dim * localIndex->second + _component ) ){

            mooseError("micro oops");

        }
//        _console << _current_node->id() << ", " << localIndex->second << ", " << _component << ": " << _overlap_coupling.getMicroAugmentedLagrangianForce( )->at( _dim * localIndex->second + _component ) << "\n";
        return _overlap_coupling.getMicroAugmentedLagrangianForce( )->at( _dim * localIndex->second + _component );

    }

    return 0;

}

Real CouplingForce::computeQpJacobian( ){
    return 0;
}
