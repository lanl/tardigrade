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
registerMooseObject("tardigradeApp", InternalForce);

template<>
InputParameters
validParams<InternalForce>(){
    InputParameters params = validParams<Kernel>();
    params.addRequiredParam<int>("component", "The component of the internal force vector");
    params.addRequiredParam<int>("dof_num",   "The degree of freedom to use for the diagonal jacobian calculation");
    params.addParam<bool>("MMS", false,       "The flag for whether the run will be using the method of manufactured solutions");
    params.addCoupledVar("u1", "The degree of freedom in the 1 direction");
    params.addCoupledVar("u2", "The degree of freedom in the 2 direction");
    params.addCoupledVar("u3", "The degree of freedom in the 3 direction");
    params.addCoupledVar("phi_11", "The 11 component of the phi tensor");
    params.addCoupledVar("phi_22", "The 22 component of the phi tensor");
    params.addCoupledVar("phi_33", "The 33 component of the phi tensor");
    params.addCoupledVar("phi_23", "The 23 component of the phi tensor");
    params.addCoupledVar("phi_13", "The 13 component of the phi tensor");
    params.addCoupledVar("phi_12", "The 12 component of the phi tensor");
    params.addCoupledVar("phi_32", "The 32 component of the phi tensor");
    params.addCoupledVar("phi_31", "The 31 component of the phi tensor");
    params.addCoupledVar("phi_21", "The 21 component of the phi tensor");
    return params;
}

InternalForce::InternalForce(const InputParameters & parameters)
    : // We have to call the constructor for the base class first
        Kernel(parameters),
        _component(getParam<int>("component")),
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
        _deformation_gradient(getMaterialProperty<std::vector<std::vector<double>>>("deformation_gradient")),
        _micro_displacement(getMaterialProperty<std::vector<double>>("micro_displacement")),
        _grad_micro_displacement(getMaterialProperty<std::vector<std::vector<double>>>("gradient_micro_displacement")),
        _PK2(getMaterialProperty<std::vector<double>>("PK2")),
        _DPK2Dgrad_u(getMaterialProperty<std::vector<std::vector<double>>>("DPK2Dgrad_u")),
        _DPK2Dphi(getMaterialProperty<std::vector<std::vector<double>>>("DPK2Dphi")),
        _DPK2Dgrad_phi(getMaterialProperty<std::vector<std::vector<double>>>("DPK2Dgrad_phi")),
        _PK2_MMS(getMaterialProperty<std::vector<double>>("PK2_MMS"))
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
    Real fint_MMS;
    
    //Copy the test function so that the balance equation function can read it
    double dNdX[3];
    for (int indx=0; indx<3; indx++){dNdX[indx] = _grad_test[_i][_qp](indx);}//p+i);}

    balance_equations::compute_internal_force(_component, dNdX, _deformation_gradient[_qp], _PK2[_qp], fint);
    //std::cout << "fint: " << fint << "\n";

    if(_MMS){
        balance_equations::compute_internal_force(_component, dNdX, _PK2_MMS[_qp], fint_MMS);
        fint -= fint_MMS;
        //std::cout << "fint - fint_MMS: " << fint << "\n";
    }

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
    double dNdX[3];
    double detadX[3];
    for (int indx=0; indx<3; indx++){
        dNdX[indx]   = _grad_test[_i][_qp](indx);
        detadX[indx] = _grad_phi[_j][_qp](indx);
    }

    balance_equations::compute_internal_force_jacobian(                _component,                 _dof_num, 
                                                                   _test[_i][_qp],                     dNdX,     _phi[_j][_qp],   detadX,
                                                       _deformation_gradient[_qp], _micro_displacement[_qp],
                                                                        _PK2[_qp],        _DPK2Dgrad_u[_qp],    _DPK2Dphi[_qp],   _DPK2Dgrad_phi[_qp],
                                                                         dfdUint);
    return dfdUint;
}

Real InternalForce::computeQpOffDiagJacobian(unsigned int jvar){
    /*!
    ==================================
    |    computeQpOffDiagJacobian    |
    ==================================

    Compute the off-diagonal terms of the jacobian
    */


    Real dfdUint;
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


    //Copy the test and interpolation functions so that the balance equation function can read it
    double dNdX[3];
    double detadX[3];
    for (int indx=0; indx<3; indx++){
        dNdX[indx]   = _grad_test[_i][_qp](indx);
        detadX[indx] = _grad_phi[_j][_qp](indx);
    }

    if(_off_diag_dof_num>=0){
        balance_equations::compute_internal_force_jacobian(                _component,        _off_diag_dof_num, 
                                                                       _test[_i][_qp],                     dNdX,     _phi[_j][_qp],   detadX,
                                                           _deformation_gradient[_qp], _micro_displacement[_qp],
                                                                            _PK2[_qp],        _DPK2Dgrad_u[_qp],    _DPK2Dphi[_qp],   _DPK2Dgrad_phi[_qp],
                                                                             dfdUint);
        return dfdUint;
    }
    else{
        return 0;
    }
}
