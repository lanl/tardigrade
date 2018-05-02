/*!
=====================================================================
|                         InternalCouple.h                          |
---------------------------------------------------------------------
| The header file for the micromorphic internal couple kernel. This |
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


#ifndef INTERNALCOUPLE_H
#define INTERNALCOUPLE_H

#include "Kernel.h"
#include <balance_equations.h> //The Micromorphic balance equations

//Forward declarations
class InternalCouple;

template <>
InputParameters validParams<InternalCouple>();

class InternalCouple : public Kernel{
    /*!
    ------------------------------------------
    |    InternalCouple kernel definition    |
    ------------------------------------------

    We define the kernel for the internal force to be

    psi (sigma_{ij} - s_{ij}) - -psi_{,k} m_{kji}

    where psi is the test function, sigma is the cauchy stress, s is the 
    symmetric stress, and m is the higher order couple stress. Note that 
    this is a vector definition so we will need to define multiple kernels 
    each obtaining a specific component of the residual.
    */
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        //The constructor definition
        InternalCouple(const InputParameters & parameters);

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

        //Parameters
        const int _component_i;
        const int _component_j;
        const int _dof_num;
        const bool _MMS;

        unsigned int _u1_int;
        unsigned int _u2_int;
        unsigned int _u3_int;

        unsigned int _phi_11_int;
        unsigned int _phi_22_int;
        unsigned int _phi_33_int;
        unsigned int _phi_23_int;
        unsigned int _phi_13_int;
        unsigned int _phi_12_int;
        unsigned int _phi_32_int;
        unsigned int _phi_31_int;
        unsigned int _phi_21_int;

        const MaterialProperty<std::vector<double>>              &_cauchy;           //The cauchy stress
        const MaterialProperty<std::vector<double>>              &_s;                //The symmetric stress
        const MaterialProperty<std::vector<double>>              &_m;                //The higher order couple stress
        const MaterialProperty<std::vector<std::vector<double>>> &_DcauchyDgrad_u;   //The gradient of the cauchy stress w.r.t. u
        const MaterialProperty<std::vector<std::vector<double>>> &_DcauchyDphi;      //The gradient of the cauchy stress w.r.t. the micro-deformation tensor
        const MaterialProperty<std::vector<std::vector<double>>> &_DcauchyDgrad_phi; //The gradient of the cauchy stress w.r.t. the spatial gradient of the micro-deformation tensor
        const MaterialProperty<std::vector<std::vector<double>>> &_DsDgrad_u;        //The gradient of the symmetric stress w.r.t. u
        const MaterialProperty<std::vector<std::vector<double>>> &_DsDphi;           //The gradient of the symmetric stress w.r.t. the micro-deformation tensor
        const MaterialProperty<std::vector<std::vector<double>>> &_DsDgrad_phi;      //The gradient of the symmetric stress w.r.t. the spatial gradient of the micro-deformation tensor
        const MaterialProperty<std::vector<std::vector<double>>> &_DmDgrad_u;        //The gradient of the higher order couple stress w.r.t. u
        const MaterialProperty<std::vector<std::vector<double>>> &_DmDphi;           //The gradient of the higher order couple stress w.r.t. the micro-deformation tensor
        const MaterialProperty<std::vector<std::vector<double>>> &_DmDgrad_phi;      //The gradient of the higher order couple stress w.r.t. the spatial gradient of the micro-deformation tensor

        const MaterialProperty<std::vector<double>>              &_cauchy_MMS;       //The cauchy stress
        const MaterialProperty<std::vector<double>>              &_s_MMS;            //The symmetric stress
        const MaterialProperty<std::vector<double>>              &_m_MMS;            //The higher order couple stress
};

#endif
