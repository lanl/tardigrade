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
//    std::cout << "In InputParameters InternalCouple\n";
    InputParameters params = validParams<Kernel>();
    params.addRequiredParam<int>("component_i", "The i component of the internal couple tensor");
    params.addRequiredParam<int>("component_j", "The j component of the internal couple tensor");
    params.addRequiredParam<int>("dof_num",   "The degree of freedom to use for the diagonal jacobian calculation");
    params.addParam<bool>("MMS", false, "The flag indicating whether the method of manufactured solutions is being used");
    params.addCoupledVar("u1", "The degree of freedom in the x direction");
    params.addCoupledVar("u2", "The degree of freedom in the y direction");
    params.addCoupledVar("u3", "The degree of freedom in the z direction");
    params.addCoupledVar("phi_11", "The xx component of the phi tensor");
    params.addCoupledVar("phi_22", "The yy component of the phi tensor");
    params.addCoupledVar("phi_33", "The zz component of the phi tensor");
    params.addCoupledVar("phi_23", "The yz component of the phi tensor");
    params.addCoupledVar("phi_13", "The yz component of the phi tensor");
    params.addCoupledVar("phi_12", "The yz component of the phi tensor");
    params.addCoupledVar("phi_32", "The yz component of the phi tensor");
    params.addCoupledVar("phi_31", "The yz component of the phi tensor");
    params.addCoupledVar("phi_21", "The yz component of the phi tensor");
    return params;
}

InternalCouple::InternalCouple(const InputParameters & parameters)
    : // We have to call the constructor for the base class first
        Kernel(parameters),
        _component_i(getParam<int>("component_i")),
        _component_j(getParam<int>("component_j")),
        _dof_num(getParam<int>("dof_num")),
        _MMS(getParam<bool>("MMS")),
        _u1_int(isCoupled("u1") ? coupled("u1")
                                : 100),
        _u2_int(isCoupled("u2") ? coupled("u2")
                                : 100),
        _u3_int(isCoupled("u3") ? coupled("u3")
                                : 100),
        _phi_11_int(isCoupled("phi_11") ? coupled("phi_11")
                                        : 100),
        _phi_22_int(isCoupled("phi_22") ? coupled("phi_22")
                                        : 100),
        _phi_33_int(isCoupled("phi_33") ? coupled("phi_33")
                                        : 100),
        _phi_23_int(isCoupled("phi_23") ? coupled("phi_23")
                                        : 100),
        _phi_13_int(isCoupled("phi_13") ? coupled("phi_13")
                                        : 100),
        _phi_12_int(isCoupled("phi_12") ? coupled("phi_12")
                                        : 100),
        _phi_32_int(isCoupled("phi_32") ? coupled("phi_32")
                                        : 100),
        _phi_31_int(isCoupled("phi_31") ? coupled("phi_31")
                                        : 100),
        _phi_21_int(isCoupled("phi_21") ? coupled("phi_21")
                                        : 100),
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
        _DmDgrad_phi(getMaterialProperty<std::vector<std::vector<double>>>("DmDgrad_phi")),
        _cauchy_MMS(getMaterialProperty<std::vector<double>>("cauchy_MMS")),
        _s_MMS(getMaterialProperty<std::vector<double>>("s_MMS")),
        _m_MMS(getMaterialProperty<std::vector<double>>("m_MMS"))
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
    Real cint_MMS;
    
    //Copy the test function so that the balance equation function can read it
    double dNdx[3];
    for (int indx=0; indx<3; indx++){dNdx[indx] = _grad_test[_i][_qp](indx);}
    
    balance_equations::compute_internal_couple(_component_i, _component_j, _test[_i][_qp], dNdx, 
                                               _cauchy[_qp], _s[_qp],      _m[_qp],
                                               cint);

    if(_MMS){
        
        balance_equations::compute_internal_couple(_component_i,     _component_j,     _test[_i][_qp], dNdx, 
                                                   _cauchy_MMS[_qp], _s_MMS[_qp],      _m_MMS[_qp],
                                                   cint_MMS);
        cint -= cint_MMS;
//        std::cout << "cint - cint_MMS: " << cint << "\n";
    }
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
    //mooseError("fail");
    //Copy the test and interpolation functions so that the balance equation function can read it
    double dNdx[3];
    double detadx[3];
    for (int indx=0; indx<3; indx++){
        dNdx[indx]   = _grad_test[_i][_qp](indx);
        detadx[indx] = _grad_phi[_j][_qp](indx);
    }


    balance_equations::compute_internal_couple_jacobian(_component_i,         _component_j,      _dof_num, 
                                                        _test[_i][_qp],        dNdx,             _phi[_j][_qp],          detadx,
                                                        _DcauchyDgrad_u[_qp], _DcauchyDphi[_qp], _DcauchyDgrad_phi[_qp],
                                                        _DsDgrad_u[_qp],      _DsDphi[_qp],      _DsDgrad_phi[_qp],
                                                        _DmDgrad_u[_qp],      _DmDphi[_qp],      _DmDgrad_phi[_qp],
                                                        dcdUint);
    return dcdUint;
}

Real InternalCouple::computeQpOffDiagJacobian(unsigned int jvar){
    /*!
    ==================================
    |    computeQpOffDiagJacobian    |
    ==================================

    Compute the off-diagonal terms of the jacobian
    */


    Real dcdUint;
    int  _off_diag_dof_num = -1;

    if(jvar == _u1_int){
        _off_diag_dof_num = 0;
    }
    else if(jvar == _u2_int){
        _off_diag_dof_num = 1;
    }
    else if(jvar == _u3_int){
        _off_diag_dof_num = 2;
    }
    else if(jvar == _phi_11_int){
        _off_diag_dof_num = 3;
    }
    else if(jvar == _phi_22_int){
        _off_diag_dof_num = 4;
    }
    else if(jvar == _phi_33_int){
        _off_diag_dof_num = 5;
    }
    else if(jvar == _phi_23_int){
        _off_diag_dof_num = 6;
    }
    else if(jvar == _phi_13_int){
        _off_diag_dof_num = 7;
    }
    else if(jvar == _phi_12_int){
        _off_diag_dof_num = 8;
    }
    else if(jvar == _phi_32_int){
        _off_diag_dof_num = 9;
    }
    else if(jvar == _phi_31_int){
        _off_diag_dof_num = 10;
    }
    else if(jvar == _phi_21_int){
        _off_diag_dof_num = 11;
    }

    double dNdx[3];
    double detadx[3];
    for (int indx=0; indx<3; indx++){
        dNdx[indx]   = _grad_test[_i][_qp](indx);
        detadx[indx] = _grad_phi[_j][_qp](indx);
    }

    if(_off_diag_dof_num >= 0){
        balance_equations::compute_internal_couple_jacobian(_component_i,         _component_j,      _off_diag_dof_num,
                                                            _test[_i][_qp],        dNdx,             _phi[_j][_qp],          detadx,
                                                            _DcauchyDgrad_u[_qp], _DcauchyDphi[_qp], _DcauchyDgrad_phi[_qp],
                                                            _DsDgrad_u[_qp],      _DsDphi[_qp],      _DsDgrad_phi[_qp],
                                                            _DmDgrad_u[_qp],      _DmDphi[_qp],      _DmDgrad_phi[_qp],
                                                            dcdUint);
        return dcdUint;
    }
    else{
        return 0;
    }
}
