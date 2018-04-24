/*!
====================================================================
|                      MicromorphicMaterial.h                      |
====================================================================
| The header file for a class which computes the cauchy stress and |
| the associated jacobians for a micromorphic  material.           |
--------------------------------------------------------------------
| Notes: Relies on libraries from the micromorphic_element         |
|        repository available at bitbucket.org/NateAM2             |
====================================================================
*/

#ifndef MICROMORPHICMATERIAL_H
#define MICROMORPHICMATERIAL_H

#include "Material.h"
#include<balance_equations.h>

//Forward declarations
class MicromorphicMaterial;

template <>
InputParameters validParams<MicromorphicMaterial>();

class MicromorphicMaterial : public Material{
    /*!
    ======================================
    |    MicromorphicMaterial    |
    ======================================

    Translation of the micromorphic_linear_elasticity library 
    available in the micromorphic_element repository for use 
    in MOOSE. More complete details of this model can be found 
    there.

    */

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        MicromorphicMaterial(const InputParameters &parameters);

    protected:
        virtual void computeQpProperties() override;

    private:
        //Coupled variables (i.e. u_i,j, phi_ij, and phi_ij,k)
        //grad u
        const VariableGradient & _grad_u1;
        const VariableGradient & _grad_u2;
        const VariableGradient & _grad_u3;

        //phi
        const Variable &_phi_11;
        const Variable &_phi_22;
        const Variable &_phi_33;
        const Variable &_phi_23;
        const Variable &_phi_13;
        const Variable &_phi_32;
        const Variable &_phi_31;
        const Variable &_phi_21;

        //grad phi
        const VariableGradient &_grad_phi_11;
        const VariableGradient &_grad_phi_22;
        const VariableGradient &_grad_phi_33;
        const VariableGradient &_grad_phi_23;
        const VariableGradient &_grad_phi_13;
        const VariableGradient &_grad_phi_12;
        const VariableGradient &_grad_phi_32;
        const VariableGradient &_grad_phi_31;
        const VariableGradient &_grad_phi_21;

        //Stresses (Material Properties in MOOSE parlance)
        MaterialProperty<Vector_9>     &_cauchy;
        MaterialProperty<Vector_9>     &_s;
        MaterialProperty<Vector_27>    &_m;

        MaterialProperty<Matrix_9x9>   &_DcauchyDgrad_u;
        MaterialProperty<Matrix_9x9>   &_DcauchyDphi;
        MaterialProperty<Matrix_9x27>  &_DcauchyDgrad_phi;

        MaterialProperty<Matrix_9x9>   &_DsDgrad_u;
        MaterialProperty<Matrix_9x9>   &_DsDphi;
        MaterialProperty<Matrix_9x27>  &_DsDgrad_phi;

        MaterialProperty<Matrix_27x9>  &_DmDgrad_u;
        MaterialProperty<Matrix_27x9>  &_DmDphi;
        MaterialProperty<Matrix_27x27> &_DmDgrad_phi;

        //TODO: Add additional values

};

#endif
