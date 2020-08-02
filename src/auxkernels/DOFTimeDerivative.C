#include "DOFTimeDerivative.h"

registerMooseObject( "tardigradeApp", DOFTimeDerivative );

template <>
InputParameters
validParams<DOFTimeDerivative>()
{
    InputParameters params = validParams<AuxKernel>();
    params.addParam< int >( "derivative_order", 1, "The order of the temporal derivative 1 or 2 is supported ( defaults to 1 )" );
    params.addRequiredCoupledVar( "coupled", "The variable to compute the temporal derivative of" );
    return params;
}

DOFTimeDerivative::DOFTimeDerivative( const InputParameters & parameters )
    : AuxKernel( parameters ),
      _derivativeOrder( getParam< int >( "derivative_order" ) ),
      _dotc( coupledDot( "coupled" ) ),
      _dotDotc( coupledDotDot( "coupled" ) )
{

    if ( ( _derivativeOrder != 1 ) && ( _derivativeOrder != 2 ) ){

        mooseError( "The derivative order of " + std::to_string( _derivativeOrder ) + " is not supported" );

    }

}

Real
DOFTimeDerivative::computeValue()
{
    if ( _var.isNodalDefined( ) ){

        if ( _derivativeOrder == 1 ){

            return _dotc[ _qp ];

        }
        else if ( _derivativeOrder == 2 ){

            return _dotDotc[ _qp ];

        }

    }

    mooseError( "DOFTimeDerivative must be a nodal AuxKernel" );
}

