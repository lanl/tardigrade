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

#include<InternalCouple.h>

//We define the valid parameters for this kernel and their default values
registerMooseObject("tardigradeApp", InternalCouple);

template<>
InputParameters
validParams<InternalCouple>(){
    InputParameters params = validParams<Kernel>();
    params.addRequiredParam<int>("component_i", "The i component of the internal couple tensor");
    params.addRequiredParam<int>("component_j", "The j component of the internal couple tensor");
    params.addRequiredParam<int>("dof_num",   "The degree of freedom to use for the diagonal jacobian calculation");
    return params;
}

InternalCouple::InternalCouple(const InputParameters & parameters)
    : // We have to call the constructor for the base class first
        Kernel(parameters),
        _component_i(getParam<int>("component_i")),
        _component_j(getParam<int>("component_j")),
        _dof_num(getParam<int>("dof_num")),
        _cauchy(getMaterialProperty<std::vector<double>>("cauchy")),
        _s(getMaterialProperty<std::vector<double>>("s")),
        _m(getMaterialProperty<std::vector<double>>("m")),
        _DcauchyDgrad_u(getMaterialProperty<std::vector<std::vector<double>>>("DcauchyDgrad_u")),
        _DcauchyDphi(getMaterialProperty<std::vector<std::vector<double>>>("DcauchyDphi")),
        _DcauchyDgrad_phi(getMaterialProperty<std::vector<std::vector<double>>>("DcauchyDgrad_phi")),
        _DsDgrad_u(getMaterialProperty<std::vector<std::vector<double>>>("DsDgrad_u")),
        _DsDphi(getMaterialProperty<std::vector<std::vector<double>>>("DsDphi")),
        _DsDgrad_phi(getMaterialProperty<std::vector<std::vector<double>>>("DsDgrad_phi")),
        _DmDgrad_u(getMaterialProperty<std::vector<std::vector<double>>>("DmDgrad_u")),
        _DmDphi(getMaterialProperty<std::vector<std::vector<double>>>("DmDphi")),
        _DmDgrad_phi(getMaterialProperty<std::vector<std::vector<double>>>("DmDgrad_phi"))
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

Real InternalCouple::computeQpResidual(){
    /*!
    ===========================
    |    computeQpResidual    |
    ===========================

    Compute the residual at the quadrature point for 
    the indicated component.

    cint_ij = pse(sigma_ij - s_ij) -psi_k m_kji

    where i = _component_i
          j = _component_j

    */
    Real cint;
    
    //Copy the test function so that the balance equation function can read it
    double dNdx[3];
    //const double *p;
    //p = &_grad_test[_i][_qp](0);
    for (int i=0; i<3; i++){dNdx[i] = *(&_grad_test[_i][_qp](i));}//p+i);}
    
    balance_equations::compute_internal_couple(_component_i, _component_j, _test[_i][_qp], dNdx, 
                                               _cauchy[_qp], _s[_qp], _m[_qp],
                                               cint);
    return cint;
}

Real InternalCouple::computeQpJacobian(){
    /*!
    ===========================
    |    computeQpJacobian    |
    ===========================

    Compute the diagonal jacobian term.

    */

    Real dcdUint;

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


    balance_equations::compute_internal_couple_jacobian(_component_i,  _component_j, _dof_num, 
                                                        _test[_i][_qp], dNdx, _phi[_j][_qp], detadx,
                                                        _DcauchyDgrad_u[_qp], _DcauchyDphi[_qp], _DcauchyDgrad_phi[_qp],
                                                        _DsDgrad_u[_qp], _DsDphi[_qp], _DsDgrad_phi[_qp],
                                                        _DmDgrad_u[_qp], _DmDphi[_qp], _DmDgrad_phi[_qp],
                                                        dcdUint);
    return dcdUint;
}
