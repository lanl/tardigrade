/*!
=====================================================================
|                    MicromorphicInertialForce.h                    |
---------------------------------------------------------------------
| The header file for the micromorphic inertial force kernel. This  |
| kernel draws upon the balance equations from the git repository   |
| micromorphic_element                                              |
| and requires three degrees of freedom to be defined for solid     |
| mechanics. These degrees of freedom are:                          |
|     u_1, u_2, u_3                                                 |
| where u_1 -> u_3 are the macro deformations                       |
=====================================================================
*/

#ifndef MICROMORPHICINERTIALFORCE_H
#define MICROMORPHICINERTIALFORCE_H

#include "Kernel.h"
#include <balance_equations.h> //The Micromorphic balance equations

//Forward declarations
class MicromorphicInertialForce;

template <>
InputParameters validParams<MicromorphicInertialForce>();

class MicromorphicInertialForce : public Kernel{
    /*!
    -----------------------------------------------------
    |    MicromorphicInertialForce kernel definition    |
    -----------------------------------------------------

    We define the kernel for the inernal force to be

    -psi rho_0 a_i

    where psi is the test function and a is the acceleration. Note that 
    this is a vector definition so we will need to define multiple kernels 
    each obtaining a specific component of the residual.
    */
    public:
        //The constructor definition
        MicromorphicInertialForce(const InputParameters & parameters);

    protected:
        //The residual at a quadrature point.
        //Note: As stated above, a kernel is a scalar quantity so we have to define
        //      as many of these as dimensions in our problem.
        virtual Real computeQpResidual() override;
        
        //The diagonal members of the residual at a quadrature point.
        //Note: The parameter chosen as the "diagonal" member is totally arbitrary
        //      We just need to have each DOF associated with some residual term.
        virtual Real computeQpJacobian() override;

        //The off diagonal members of the residual at a quadrature point.
        virtual Real computeQpOffDiagJacobian(unsigned int jvar) override; 

        //Parameters
        const int _component;
        const int _dof_num;
        const bool _MMS;

        //Coupled degrees of freedom
        unsigned int _u1_int;
        unsigned int _u2_int;
        unsigned int _u3_int;

        const Real _density;
        const VariableValue &_a1;
        const VariableValue &_a2;
        const VariableValue &_a3;
        const VariableValue &_da1du;
        const VariableValue &_da2du;
        const VariableValue &_da3du;


};

#endif
