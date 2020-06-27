#include "NodalVolumeAverages.h"

registerMooseObject( "tardigradeApp", NodalVolumeAverages );

template <>
InputParameters
validParams<NodalVolumeAverages>()
{
    InputParameters params = validParams<AuxKernel>();
    params.addRequiredParam<UserObjectName>( "ElementIntegrateUserObject", "The coupled element integrate user object" );
    params.addParam<bool>( "compute_nodal_volume", false, "Flag for whether the nodal volume should be returned or not" );
    return params;
}

NodalVolumeAverages::NodalVolumeAverages( const InputParameters & parameters )
    : AuxKernel( parameters ),
      _elementIntegrate( getUserObject< ElementIntegrateUserObject >( "ElementIntegrateUserObject" ) ),
      _compute_nodal_volume( getParam< bool >( "compute_nodal_volume" ) )
{
}

Real
NodalVolumeAverages::computeValue()
{
    if ( _var.isNodalDefined( ) ){

        if ( _compute_nodal_volume ){
            return _elementIntegrate.get_nodal_volume( _current_node->id( ) );
        }
        else{        
            return _elementIntegrate.get_nodal_density( _current_node->id( ) );
        }

    }

    mooseError( "NodalVolumeAverages must be a nodal AuxKernel" );
}

