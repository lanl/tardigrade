#pragma once

#include "AuxKernel.h"

class DOFTimeDerivative;

template <>
InputParameters validParams<DOFTimeDerivative>();

class DOFTimeDerivative : public AuxKernel
{

    public:
        DOFTimeDerivative( const InputParameters & parameters );

    protected:
        virtual Real computeValue( ) override;
        const int _derivativeOrder;

};
