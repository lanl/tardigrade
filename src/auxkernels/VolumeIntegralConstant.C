#include "VolumeIntegralConstant.h"

registerMooseObject("tardigradeApp", VolumeIntegralConstant );

InputParameters
VolumeIntegralConstant::validParams()
{
    InputParameters params = AuxKernel::validParams();
    params.addClassDescription( "Computes the volume weighted integral of a constant value" );
    params.addParam<Real>("value", 0.0, "Some constant value that can be read from the input file" );
    params.declareControllable("value");
    return params;
}

VolumeIntegralConstant::VolumeIntegralConstant( const InputParameters & parameters )
    : AuxKernel( parameters ), _value( getParam< Real >( "value" ) )
{
}

Real
VolumeIntegralConstant::computeValue()
{
    return _value * ( _bnd ? _current_side_volume : _current_elem_volume );
}
