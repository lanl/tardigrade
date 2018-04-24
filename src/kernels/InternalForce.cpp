/*!
=====================================================================
|                         InternalForce.cpp                         |
---------------------------------------------------------------------
| The source file for the micromorphic internal force kernel. This  |
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

#include<InternalForce.h>

//We define the valid parameters for this kernel and their default values
registerMooseObject("tardigrade", InternalForce);

template<>
InputParameters
validParams<InputParameters>(){
    InputParameters params = validParams<Kernel>();
    params.addRequiredParam<int>("component", "The component of the internal force vector");
    params.addRequiredParam<int>("dof_num",   "The degree of freedom to use for the diagonal jacobian calculation");
    return params;
}

InternalForce::InternalForce(const InputParameters & parameters)
    : // We have to call the constructor for the base class first
        Kernel(parameters),
        _component(getParam<int>("component")),
        _compoennt(getParam<int>("dof_num")),
        _cauchy(getMaterialProperty<Vector_9>("cauchy")),
        _DcauchyDgrad_u(getMaterialProperty<Matrix_9x9>("DcauchyDgrad_u")),
        _Dcauchydphi(getMaterialProperty<Matrix_9x9>("DcauchyDphi")),
        _DcauchyDgrad_phi(getMaterialProperty<Matrix_9x27>("DcauchyDgrad_phi"))
    {
    /*!
    =====================
    |    Constructor    |
    =====================

    The constructor for the InternalForce class.
    Note that this constructor is just variable 
    assignments.

    */
}

Real InternalForce::computeQpResidual(){
    /*!
    ===========================
    |    computeQpResidual    |
    ===========================

    Compute the residual at the quadrature point for 
    the indicated component.

    fint_i = -psi_j sigma_ji

    where i = _component

    */
    Real fint;
    balance_equations::compute_internal_force(_component, _grad_test[_i][_qp], _cauchy, fint);
    return fint;
}

Real InternalForce::computeQpJacobian(){
    /*!
    ===========================
    |    computeQpJacobian    |
    ===========================

    Compute the diagonal jacobian term.

    */

    Real dfdUint;
    balance_equations::compute_internal_force_jacobian(_component,           _dof_num, 
                                                       _test[_i][_qp],       _grad_test[_i][_qp], _phi[_j][_qp],          _grad_phi[_j][_qp],
                                                       _DcauchyDgrad_u[_qp], _Dcauchydphi[_qp],   _DcauchyDgrad_phi[_qp], dfdUint);
    return dfdUint;
}
