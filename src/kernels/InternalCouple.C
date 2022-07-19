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

InputParameters
InternalCouple::validParams()
{
//    std::cout << "In InputParameters InternalCouple\n";
    InputParameters params = Kernel::validParams();
    params.set< bool >( "use_displaced_mesh" ) = false;
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
    params.addCoupledVar("phi_13", "The xz component of the phi tensor");
    params.addCoupledVar("phi_12", "The xy component of the phi tensor");
    params.addCoupledVar("phi_32", "The zy component of the phi tensor");
    params.addCoupledVar("phi_31", "The zx component of the phi tensor");
    params.addCoupledVar("phi_21", "The yx component of the phi tensor");
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
        _deformation_gradient(getMaterialProperty<std::vector<double>>("MM_deformation_gradient")),
        _micro_deformation(getMaterialProperty<std::vector<double>>("micro_deformation")),
        _gradient_micro_displacement(getMaterialProperty<std::vector<std::vector<double>>>("gradient_micro_displacement")),
        _PK2(getMaterialProperty<std::vector<double>>("PK2")),
        _SIGMA(getMaterialProperty<std::vector<double>>("SIGMA")),
        _M(getMaterialProperty<std::vector<double>>("M")),
        _DPK2Dgrad_u(getMaterialProperty<std::vector<std::vector<double>>>("DPK2Dgrad_u")),
        _DPK2Dphi(getMaterialProperty<std::vector<std::vector<double>>>("DPK2Dphi")),
        _DPK2Dgrad_phi(getMaterialProperty<std::vector<std::vector<double>>>("DPK2Dgrad_phi")),
        _DSIGMADgrad_u(getMaterialProperty<std::vector<std::vector<double>>>("DSIGMADgrad_u")),
        _DSIGMADphi(getMaterialProperty<std::vector<std::vector<double>>>("DSIGMADphi")),
        _DSIGMADgrad_phi(getMaterialProperty<std::vector<std::vector<double>>>("DSIGMADgrad_phi")),
        _DMDgrad_u(getMaterialProperty<std::vector<std::vector<double>>>("DMDgrad_u")),
        _DMDphi(getMaterialProperty<std::vector<std::vector<double>>>("DMDphi")),
        _DMDgrad_phi(getMaterialProperty<std::vector<std::vector<double>>>("DMDgrad_phi")),
        _PK2_MMS(getMaterialProperty<std::vector<double>>("PK2_MMS")),
        _SIGMA_MMS(getMaterialProperty<std::vector<double>>("SIGMA_MMS")),
        _M_MMS(getMaterialProperty<std::vector<double>>("M_MMS"))
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

    cint_ij = psi (PK2_ij - SIGMA_ij) -psi_k M_kji

    where i = _component_i
          j = _component_j

    */
    Real cint;
    Real cint_MMS;
    
    //Copy the test function so that the balance equation function can read it
    double dNdX[3];
    for (int indx=0; indx<3; indx++){dNdX[indx] = _grad_test[_i][_qp](indx);}
    
    balance_equations::compute_internal_couple(_component_i,               _component_j,             _test[_i][_qp], dNdX,
                                               _deformation_gradient[_qp], _micro_deformation[_qp],   
                                               _PK2[_qp],                  _SIGMA[_qp],              _M[_qp],
                                               cint);

    if(_MMS){
        
        balance_equations::compute_internal_couple(_component_i,               _component_j,             _test[_i][_qp], dNdX,
                                                   _deformation_gradient[_qp], _micro_deformation[_qp],
                                                   _PK2_MMS[_qp],              _SIGMA_MMS[_qp],          _M_MMS[_qp],
                                                   cint_MMS);
        cint -= cint_MMS;
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
    //Copy the test and interpolation functions so that the balance equation function can read it
    double dNdX[3];
    double detadX[3];
    for (int indx=0; indx<3; indx++){
        dNdX[indx]   = _grad_test[_i][_qp](indx);
        detadX[indx] = _grad_phi[_j][_qp](indx);
    }

    balance_equations::compute_internal_couple_jacobian( _dim * _component_i + _component_j, _dof_num, 
                                                         _test[_i][_qp],             dNdX,
                                                         _phi[_j][_qp],              detadX,
                                                         _deformation_gradient[_qp], _micro_deformation[_qp],
                                                         _PK2[_qp],                  _SIGMA[_qp],              _M[_qp],
                                                         _DPK2Dgrad_u[_qp],          _DPK2Dphi[_qp],           _DPK2Dgrad_phi[_qp],
                                                         _DSIGMADgrad_u[_qp],        _DSIGMADphi[_qp],         _DSIGMADgrad_phi[_qp],
                                                         _DMDgrad_u[_qp],            _DMDphi[_qp],             _DMDgrad_phi[_qp],
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
    else if(jvar == _phi_12_int){
        _off_diag_dof_num = 4;
    }
    else if(jvar == _phi_13_int){
        _off_diag_dof_num = 5;
    }
    else if(jvar == _phi_21_int){
        _off_diag_dof_num = 6;
    }
    else if(jvar == _phi_22_int){
        _off_diag_dof_num = 7;
    }
    else if(jvar == _phi_23_int){
        _off_diag_dof_num = 8;
    }
    else if(jvar == _phi_31_int){
        _off_diag_dof_num = 9;
    }
    else if(jvar == _phi_32_int){
        _off_diag_dof_num = 10;
    }
    else if(jvar == _phi_33_int){
        _off_diag_dof_num = 11;
    }

    double dNdX[3];
    double detadX[3];
    for (int indx=0; indx<3; indx++){
        dNdX[indx]   = _grad_test[_i][_qp](indx);
        detadX[indx] = _grad_phi[_j][_qp](indx);
    }

    if(_off_diag_dof_num >= 0){
        balance_equations::compute_internal_couple_jacobian( _dim * _component_i +  _component_j,         _off_diag_dof_num, 
                                                             _test[_i][_qp],             dNdX,
                                                             _phi[_j][_qp],              detadX,
                                                             _deformation_gradient[_qp], _micro_deformation[_qp],
                                                             _PK2[_qp],                  _SIGMA[_qp],              _M[_qp],
                                                             _DPK2Dgrad_u[_qp],          _DPK2Dphi[_qp],           _DPK2Dgrad_phi[_qp],
                                                             _DSIGMADgrad_u[_qp],        _DSIGMADphi[_qp],         _DSIGMADgrad_phi[_qp],
                                                             _DMDgrad_u[_qp],            _DMDphi[_qp],             _DMDgrad_phi[_qp],
                                                             dcdUint);

        return dcdUint;
    }
    else{
        return 0;
    }
}
