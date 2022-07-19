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
#include "MooseException.h"
#include<balance_equations.h>
#include<micromorphic_material_library.h>

//Forward declarations
class MicromorphicMaterial;
class Function;

class MicromorphicMaterial : public Material
{
    /*!
    ==============================
    |    MicromorphicMaterial    |
    ==============================

    Translation of the micromorphic_linear_elasticity library 
    available in the micromorphic_element repository for use 
    in MOOSE. More complete details of this model can be found 
    there.

    */

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        MicromorphicMaterial(const InputParameters &parameters);
        
        static InputParameters validParams();

    protected:
        virtual void initQpStatefulProperties() override;
        virtual void computeQpProperties() override;

    private:

        //Parameters
        std::vector<Real> _fparams;
        int               _n_ADD_DOF;
        int               _n_ADD_TERMS;
        int               _n_ADD_JACOBIANS;
        int               _n_SDVS;
        std::string       _model_name;
        bool              _MMS;

        //Old variable values
        //grad u
        const VariableValue    & _old_u1;
        const VariableValue    & _old_u2;
        const VariableValue    & _old_u3;
        const VariableGradient & _old_grad_u1;
        const VariableGradient & _old_grad_u2;
        const VariableGradient & _old_grad_u3;

        //phi
        const VariableValue    & _old_phi_11;
        const VariableValue    & _old_phi_22;
        const VariableValue    & _old_phi_33;
        const VariableValue    & _old_phi_23;
        const VariableValue    & _old_phi_13;
        const VariableValue    & _old_phi_12;
        const VariableValue    & _old_phi_32;
        const VariableValue    & _old_phi_31;
        const VariableValue    & _old_phi_21;

        //grad phi
        const VariableGradient & _old_grad_phi_11;
        const VariableGradient & _old_grad_phi_22;
        const VariableGradient & _old_grad_phi_33;
        const VariableGradient & _old_grad_phi_23;
        const VariableGradient & _old_grad_phi_13;
        const VariableGradient & _old_grad_phi_12;
        const VariableGradient & _old_grad_phi_32;
        const VariableGradient & _old_grad_phi_31;
        const VariableGradient & _old_grad_phi_21;


        //Coupled variables (i.e. u_i,j, phi_ij, and phi_ij,k)
        //grad u
        const VariableValue    & _u1;
        const VariableValue    & _u2;
        const VariableValue    & _u3;
        const VariableGradient & _grad_u1;
        const VariableGradient & _grad_u2;
        const VariableGradient & _grad_u3;

        //phi
        const VariableValue    & _phi_11;
        const VariableValue    & _phi_22;
        const VariableValue    & _phi_33;
        const VariableValue    & _phi_23;
        const VariableValue    & _phi_13;
        const VariableValue    & _phi_12;
        const VariableValue    & _phi_32;
        const VariableValue    & _phi_31;
        const VariableValue    & _phi_21;

        //grad phi
        const VariableGradient & _grad_phi_11;
        const VariableGradient & _grad_phi_22;
        const VariableGradient & _grad_phi_33;
        const VariableGradient & _grad_phi_23;
        const VariableGradient & _grad_phi_13;
        const VariableGradient & _grad_phi_12;
        const VariableGradient & _grad_phi_32;
        const VariableGradient & _grad_phi_31;
        const VariableGradient & _grad_phi_21;

        //Fundamental deformation measures
        MaterialProperty<std::vector<double>> & _deformation_gradient;
        MaterialProperty<std::vector<double>> & _micro_deformation;
        MaterialProperty<std::vector<std::vector<double>>> & _gradient_micro_displacement;

        //Stresses (Material Properties in MOOSE parlance)
        MaterialProperty<std::vector<double>>               & _cauchy;
        MaterialProperty<std::vector<double>>               & _s;
        MaterialProperty<std::vector<double>>               & _m;

        MaterialProperty<std::vector<double>>               & _PK2;
        MaterialProperty<std::vector<double>>               & _SIGMA;
        MaterialProperty<std::vector<double>>               & _M;

        MaterialProperty<std::vector<std::vector<double>>>  & _DPK2Dgrad_u;
        MaterialProperty<std::vector<std::vector<double>>>  & _DPK2Dphi;
        MaterialProperty<std::vector<std::vector<double>>>  & _DPK2Dgrad_phi;

        MaterialProperty<std::vector<std::vector<double>>>  & _DSIGMADgrad_u;
        MaterialProperty<std::vector<std::vector<double>>>  & _DSIGMADphi;
        MaterialProperty<std::vector<std::vector<double>>>  & _DSIGMADgrad_phi;

        MaterialProperty<std::vector<std::vector<double>>>  & _DMDgrad_u;
        MaterialProperty<std::vector<std::vector<double>>>  & _DMDphi;
        MaterialProperty<std::vector<std::vector<double>>>  & _DMDgrad_phi;

        //State variables
        MaterialProperty< std::vector< double > > & _SDVS;
        const MaterialProperty< std::vector< double > > & _old_SDVS;

        //TODO: Add additional values
        MaterialProperty<std::vector<std::vector<double>>>  & _ADD_TERMS;
        MaterialProperty<std::vector<std::vector<std::vector<double>>>>  & _ADD_JACOBIANS;

        //MMS function values
        const Function & _u1_fxn;
        const Function & _u2_fxn;
        const Function & _u3_fxn;
        const Function & _phi_11_fxn;
        const Function & _phi_22_fxn;
        const Function & _phi_33_fxn;
        const Function & _phi_23_fxn;
        const Function & _phi_13_fxn;
        const Function & _phi_12_fxn;
        const Function & _phi_32_fxn;
        const Function & _phi_31_fxn;
        const Function & _phi_21_fxn;

        //MMS stress measures
        MaterialProperty<std::vector<double>> & _PK2_MMS;
        MaterialProperty<std::vector<double>> & _SIGMA_MMS;
        MaterialProperty<std::vector<double>> & _M_MMS;
        MaterialProperty<std::vector<std::vector<double>>> & _ADD_TERMS_MMS;
};

#endif
