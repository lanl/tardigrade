#pragma once

#include "AuxKernel.h"

class VolumeIntegralConstant;

class VolumeIntegralConstant : public AuxKernel
{

    public:
        static InputParameters validParams();

        VolumeIntegralConstant( const InputParameters & parameters );

    protected:
        virtual Real computeValue( ) override;

        const Real & _value;

};
