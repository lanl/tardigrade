#include<CoupledDirichletBC.h>

//Register the object
registerMooseObject( "tardigradeApp", CoupledDirichletBC );

template<>
InputParameters
validParams< CoupledDirichletBC >( ){
    InputParameters params = validParams<DirichletBCBase>( );
    params.addRequiredParam< UserObjectName >( "overlap_coupling_object", "The overlap coupling user object" );
    params.addRequiredParam< unsigned int >( "component", "The component of the DOF vector being constrained" );
    params.addRequiredParam< bool > ( "is_macroscale", "Flag for whether the domain is the macroscale or the micro-scale" );
    return params;
}

CoupledDirichletBC::CoupledDirichletBC( const InputParameters &parameters )
    : // Call the constructor for the base class
        DirichletBCBase( parameters ),
        _overlap_coupling( getUserObject< OverlapCoupling >( "overlap_coupling_object" ) ),
        _component( getParam< unsigned int >( "component" ) ),
        _isMacro( getParam< bool >( "is_macroscale" ) ){

}

Real CoupledDirichletBC::computeQpValue( ){

    unsigned int _dim = _mesh.dimension( );

    if ( _isMacro ){

        auto localIndex = _overlap_coupling.getMacroGlobalLocalNodeMap( )->find( _current_node->id( ) );

        if ( localIndex == _overlap_coupling.getMacroGlobalLocalNodeMap( )->end( ) ){

            mooseError( "The current node is not found in the global-local node map" );

        }

        return _overlap_coupling.getMacroDisplacementDOF( )->at( ( _dim * _dim + _dim ) * localIndex->second + _component );


    }
    else{

        auto localIndex = _overlap_coupling.getMicroGlobalLocalNodeMap( )->find( _current_node->id( ) );

        if ( localIndex == _overlap_coupling.getMicroGlobalLocalNodeMap( )->end( ) ){

            mooseError( "The current node is not found in the global-local node map" );

        }

        return _overlap_coupling.getMicroDisplacementDOF( )->at( _dim * localIndex->second + _component );

    }

    return 0;

}

bool CoupledDirichletBC::shouldApply( ){

    //If the current node is in the DOF ID map then it should apply
    if ( _isMacro ){
        if ( _overlap_coupling.getMacroGlobalLocalNodeMap( )->find( _current_node->id( ) ) 
                 != _overlap_coupling.getMacroGlobalLocalNodeMap( )->end( ) ){

            return true;

        }
    }
    else {
        if ( _overlap_coupling.getMicroGlobalLocalNodeMap( )->find( _current_node->id( ) ) 
                 != _overlap_coupling.getMicroGlobalLocalNodeMap( )->end( ) ){

            return true;

        }
    }

    return false;

}
