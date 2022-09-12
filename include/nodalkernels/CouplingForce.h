#pragma once

#include<NodalKernel.h>
#include<OverlapCoupling.h>

class CouplingForce;

class CouplingForce : public NodalKernel
{

    public:

        static InputParameters validParams();

        CouplingForce( const InputParameters &parameters );

    protected:

        virtual Real computeQpResidual( ) override;

        virtual Real computeQpJacobian( ) override;

//        virtual bool shouldApply( ) override;

    private:

        //The user object which defines the values
        const OverlapCoupling &_overlap_coupling;

        //The component of the variable ( 0 - 2 for displacement, 3 - 11 for micro displacement )
        const unsigned int &_component;

        //Whether the force is applying to the macro-scale or micro-scale
        const bool &_isMacro;

};
