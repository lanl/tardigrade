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
#include<micromorphic_material_library.h>

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

        //Parameters
        std::vector<Real> _fparams;
        int               _n_ADD_DOF;
        int               _n_ADD_TERMS;
        int               _n_ADD_JACOBIANS;
        std::string       _model_name;

        //Coupled variables (i.e. u_i,j, phi_ij, and phi_ij,k)
        //grad u
        const VariableGradient & _grad_u1;
        const VariableGradient & _grad_u2;
        const VariableGradient & _grad_u3;

        //phi
        const VariableValue &_phi_11;
        const VariableValue &_phi_22;
        const VariableValue &_phi_33;
        const VariableValue &_phi_23;
        const VariableValue &_phi_13;
        const VariableValue &_phi_12;
        const VariableValue &_phi_32;
        const VariableValue &_phi_31;
        const VariableValue &_phi_21;

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
        MaterialProperty<std::vector<double>>               &_cauchy;
        MaterialProperty<std::vector<double>>               &_s;
        MaterialProperty<std::vector<double>>               &_m;

        MaterialProperty<std::vector<std::vector<double>>>  &_DcauchyDgrad_u;
        MaterialProperty<std::vector<std::vector<double>>>  &_DcauchyDphi;
        MaterialProperty<std::vector<std::vector<double>>>  &_DcauchyDgrad_phi;

        MaterialProperty<std::vector<std::vector<double>>>  &_DsDgrad_u;
        MaterialProperty<std::vector<std::vector<double>>>  &_DsDphi;
        MaterialProperty<std::vector<std::vector<double>>>  &_DsDgrad_phi;

        MaterialProperty<std::vector<std::vector<double>>>  &_DmDgrad_u;
        MaterialProperty<std::vector<std::vector<double>>>  &_DmDphi;
        MaterialProperty<std::vector<std::vector<double>>>  &_DmDgrad_phi;

        //TODO: Add additional values
        MaterialProperty<std::vector<std::vector<double>>>  &_ADD_TERMS;

};

#endif
