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
        _dof_num(getParam<int>("dof_num")),
        _cauchy(getMaterialProperty<Vector_9>("cauchy")),
        _DcauchyDgrad_u(getMaterialProperty<Matrix_9x9>("DcauchyDgrad_u")),
        _DcauchyDphi(getMaterialProperty<Matrix_9x9>("DcauchyDphi")),
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
    
    //Copy the test function so that the balance equation function can read it
    double dNdx[3];
    //const double *p;
    //p = &_grad_test[_i][_qp](0);
    for (int i=0; i<3; i++){dNdx[i] = *(&_grad_test[_i][_qp](i));}//p+i);}
    
    balance_equations::compute_internal_force(_component, dNdx, _cauchy[_qp], fint);
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

    //Copy the test and interpolation functions so that the balance equation function can read it
    double dNdx[3];
    double detadx[3];
    //const double *p;
    //const double *__p;
    //p   = &_grad_test[_i][_qp](0);
    //__p = &_grad_phi[_i][_qp](0);
    for (int i=0; i<3; i++){
        dNdx[i]   = *(&_grad_test[_i][_qp](i));//p+i);
        detadx[i] = *(&_grad_phi[_i][_qp](i));//__p+i);;
    }


    balance_equations::compute_internal_force_jacobian(_component,           _dof_num, 
                                                       _test[_i][_qp],       dNdx, _phi[_j][_qp],          detadx,
                                                       _DcauchyDgrad_u[_qp], _DcauchyDphi[_qp],   _DcauchyDgrad_phi[_qp], dfdUint);
    return dfdUint;
}
