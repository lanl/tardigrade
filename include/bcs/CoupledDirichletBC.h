
#ifndef COUPLEDDIRICHLETBC_H
#define COUPLEDDIRICHELTBC_H

#include<DirichletBCBase.h>
#include<OverlapCoupling.h>

class CoupledDirichletBC;

template <>
InputParameters validParams<CoupledDirichletBC>();

class CoupledDirichletBC : public DirichletBCBase{

    public:

        CoupledDirichletBC( const InputParameters &parameters );

    protected:

        virtual Real computeQpValue( ) override;

        virtual bool shouldApply( ) override;

        //The user object which defines the values
        const OverlapCoupling &_overlap_coupling;

        //The component of the variable ( 0 - 2 for displacement, 3 - 11 for micro displacement )
        const unsigned int &_component;

        //Whether the BC is applying to the macro-scale or micro-scale
        const bool &_isMacro;

};

#endif
