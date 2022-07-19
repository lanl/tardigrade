#pragma once

#include "AuxKernel.h"

class DOFTimeDerivative;


class DOFTimeDerivative : public AuxKernel
{

    public:
        DOFTimeDerivative( const InputParameters & parameters );
        
        static InputParameters validParams();

    protected:
        virtual Real computeValue( ) override;
        const int _derivativeOrder;
        const VariableValue &_dotc;
        const VariableValue &_dotDotc;

};
