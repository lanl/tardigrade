/*!
=====================================================================
|                   MicromorphicInertialCouple.h                    |
---------------------------------------------------------------------
| The header file for the micromorphic inertial couple kernel. This |
| kernel draws upon the balance equations from the git repository   |
| micromorphic_element.                                             |
| and requires nine degrees of freedom to be defined for solid      |
| mechanics. These degrees of freedom are:                          |
|     phi_11, phi_12, phi_13, phi_21, phi_22, phi_23, phi_31,       |
|     phi_32, phi_33                                                |
| which are the components of the micro-displacement tensor.        |
=====================================================================
*/


#ifndef MICROMORPHICINERTIALCOUPLE_H
#define MICROMORPHICINERTIALCOUPLE_H

#include "Kernel.h"
#include <balance_equations.h> //The Micromorphic balance equations

//Forward declarations
class MicromorphicInertialCouple;

class MicromorphicInertialCouple : public Kernel
{
    /*!
    ------------------------------------------
    |    MicromorphicInertialCouple kernel definition    |
    ------------------------------------------

    We define the kernel for the inertial couple to be

    -psi \rho_0 \ddot{ \chi }_{iI } \chi_{jJ} I_{IJ}

    where psi is the test function, \rho_0 is the density in the reference
    configuration, \ddot{ \chi }_{iI} is the second temporal derivative of 
    the micro deformation, \chi_{jJ} is the micro-deformation, and I_{IJ}
    is the reference inertia.
    Note that this is a vector definition so we will need to define multiple 
    kernels each obtaining a specific component of the residual.
    */
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        //The constructor definition
        MicromorphicInertialCouple(const InputParameters & parameters);
        
        static InputParameters validParams();

    protected:
        //The residual at a quadrature point.
        //Note: As stated above, a kernel is a scalar quantity so we have to define
        //      as many of these as dimensions in our problem.
        virtual Real computeQpResidual() override;
        
        //The diagonal members of the jacobian at a quadrature point.
        //Note: The parameter chosen as the "diagonal" member is totally arbitrary
        //      We just need to have each DOF associated with some residual term.
        virtual Real computeQpJacobian() override;

        //The off-diagonal members of the jacobian at the quadrature point.
        virtual Real computeQpOffDiagJacobian(unsigned int jvar) override; 

        //Assume 3D
        const unsigned int _dim = 3;

        //Parameters
        const int _component_i;
        const int _component_j;
        const int _dof_num;
        const bool _MMS;

        //Coupled degrees of freedom
        unsigned int _phi11_int;
        unsigned int _phi22_int;
        unsigned int _phi33_int;
        unsigned int _phi23_int;
        unsigned int _phi13_int;
        unsigned int _phi12_int;
        unsigned int _phi32_int;
        unsigned int _phi31_int;
        unsigned int _phi21_int;

        const Real _density;
        const Real _I11;
        const Real _I12;
        const Real _I13;
        const Real _I22;
        const Real _I23;
        const Real _I33;

        const VariableValue &_phi11;
        const VariableValue &_phi12;
        const VariableValue &_phi13;
        const VariableValue &_phi21;
        const VariableValue &_phi22;
        const VariableValue &_phi23;
        const VariableValue &_phi31;
        const VariableValue &_phi32;
        const VariableValue &_phi33;

        const VariableValue &_dotDotChi11;
        const VariableValue &_dotDotChi12;
        const VariableValue &_dotDotChi13;
        const VariableValue &_dotDotChi21;
        const VariableValue &_dotDotChi22;
        const VariableValue &_dotDotChi23;
        const VariableValue &_dotDotChi31;
        const VariableValue &_dotDotChi32;
        const VariableValue &_dotDotChi33;

        const VariableValue &_dDotDotChi11Du;
        const VariableValue &_dDotDotChi12Du;
        const VariableValue &_dDotDotChi13Du;
        const VariableValue &_dDotDotChi21Du;
        const VariableValue &_dDotDotChi22Du;
        const VariableValue &_dDotDotChi23Du;
        const VariableValue &_dDotDotChi31Du;
        const VariableValue &_dDotDotChi32Du;
        const VariableValue &_dDotDotChi33Du;
};

#endif
