/*!
=====================================================================
|                          InternalForce.h                          |
---------------------------------------------------------------------
| The header file for the micromorphic internal force kernel. This  |
| kernel draws upon the balance equations from the git repository   |
| micromorphic_element available at:                                |
|     https://bitbucket.org/NateAM2/                                |
| and requires twelve degrees of freedom to be defined for solid    |
| mechanics. These degrees of freedom are:                          |
|     u_1, u_2, u_3, phi_11, phi_22, phi_33, phi_23, phi_13, phi,12 |
|     phi_32, phi_31, phi_21                                        |
| where u_1 -> u_3 are macro deformations, and phi_11 -> phi_21 are |
| the components of the micro-displacement tensor.                  |
=====================================================================
*/


#ifndef INTERNALFORCE_H
#define INTERNALFORCE_H

#include "Kernel.h"
#include <balance_equations.h> //The Micromorphic balance equations

//Forward declarations
class InternalForce;

template <>
InputParameters validParams<InternalForce>();

class InternalForce : public Kernel{
    /*!
    -----------------------------------------
    |    InternalForce kernel definition    |
    -----------------------------------------

    We define the kernel for the internal force to be

    -psi_{,j} sigma_{ji}

    where psi is the test function and sigma is the cauchy stress. Note that 
    this is a vector definition so we will need to define multiple kernels 
    each obtaining a specific component of the residual.
    */
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        //The constructor definition
        InternalForce(const InputParameters & parameters);

    protected:
        //The residual at a quadrature point.
        //Note: As stated above, a kernel is a scalar quantity so we have to define
        //      as many of these as dimensions in our problem.
        virtual Real computeQpResidual() override;
        
        //The diagonal members of the residual at a quadrature point.
        //Note: The parameter chosen as the "diagonal" member is totally arbitrary
        //      We just need to have each DOF associated with some residual term.
        virtual Real computeQpJacobian() override;

        //Parameters
        const int _component;
        const int _dof_num;

        const MaterialProperty<std::vector<double>>              &_cauchy;           //The cauchy stress
        const MaterialProperty<std::vector<std::vector<double>>> &_DcauchyDgrad_u;   //The gradient of the cauchy stress w.r.t. u
        const MaterialProperty<std::vector<std::vector<double>>> &_DcauchyDphi;      //The gradient of the cauchy stress w.r.t. the micro-deformation tensor
        const MaterialProperty<std::vector<std::vector<double>>> &_DcauchyDgrad_phi; //The gradient of the cauchy stress w.r.t. the spatial gradient of the micro-deformation tensor
};

#endif
