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

InternalForce::InternalForce(const InputParameters & parameters)
    : // We have to call the constructor for the base class first
        Kernel(parameters),
        _component(getParam<int>("component")),
        _dof_num(getParam<int>("dof_num")),
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
        _DcauchyDgrad_u(getMaterialProperty<std::vector<std::vector<double>>>("DcauchyDgrad_u")),
        _DcauchyDphi(getMaterialProperty<std::vector<std::vector<double>>>("DcauchyDphi")),
        _DcauchyDgrad_phi(getMaterialProperty<std::vector<std::vector<double>>>("DcauchyDgrad_phi"))

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
    for (int indx=0; indx<3; indx++){dNdx[indx] = _grad_test[_i][_qp](indx);}//p+i);}

    balance_equations::compute_internal_force(_component, dNdx, _cauchy[_qp], fint);
//    std::cout << "fint: " << fint << "\n";
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
    for (int indx=0; indx<3; indx++){
        dNdx[indx]   = _grad_test[_i][_qp](indx);
        detadx[indx] = _grad_phi[_i][_qp](indx);
    }

    balance_equations::compute_internal_force_jacobian(_component,           _dof_num, 
                                                       _test[_i][_qp],       dNdx, _phi[_j][_qp],          detadx,
                                                       _DcauchyDgrad_u[_qp], _DcauchyDphi[_qp],   _DcauchyDgrad_phi[_qp], dfdUint);

//    std::cout << "N:    " << _test[_i][_qp] << "\n";
//    std::cout << "dNdx: ";
//    for (int indx=0; indx<3; indx++){ std::cout << dNdx[indx] << " ";}
//    std::cout << "\n";
//
//    std::cout << "eta:    " << _phi[_i][_qp] << "\n";
//    std::cout << "detadx: ";
//    for (int indx=0; indx<3; indx++){ std::cout << detadx[indx] << " ";}
//    std::cout << "\n";
//
//    std::cout << "_DcauchyDgrad_u:\n";
//    for (int indx=0; indx<_DcauchyDgrad_u[_qp].size(); indx++){
//        for (int jndx=0; jndx<_DcauchyDgrad_u[_qp][indx].size(); jndx++){
//            std::cout << _DcauchyDgrad_u[_qp][indx][jndx] << " ";
//        }
//        std::cout << "\n";
//    }
//
//    std::cout << "_DcauchyDphi:\n";
//    for (int indx=0; indx<_DcauchyDphi[_qp].size(); indx++){
//        for (int jndx=0; jndx<_DcauchyDphi[_qp][indx].size(); jndx++){
//            std::cout << _DcauchyDphi[_qp][indx][jndx] << " ";
//        }
//        std::cout << "\n";
//    }
//
//    std::cout << "_DcauchyDgrad_phi:\n";
//    for (int indx=0; indx<_DcauchyDgrad_phi[_qp].size(); indx++){
//        for (int jndx=0; jndx<_DcauchyDgrad_phi[_qp][indx].size(); jndx++){
//            std::cout << _DcauchyDgrad_phi[_qp][indx][jndx] << " ";
//        }
//        std::cout << "\n";
//    }
//
//
//
//    std::cout << "(i,dof): (" << _component << ", " <<  _dof_num << ")\n";
//    std::cout << "N: " << _test[_i][_qp] << "\n";
//    std::cout << "dNdx: " << dNdx[0] << ", " << dNdx[1] << ", " << dNdx[2] << "\n";
//    std::cout << "eta: " << _phi[_j][_qp] << "\n";
//    std::cout << "detadx: " << detadx[0] << ", " << detadx[1] << ", " << detadx[2] << "\n";
//    std::cout << "cauchy: ";
//    for (int i=0; i<9; i++){std::cout << _cauchy[_qp][i] << " ";}
//    std::cout << "\n";
//    std::cout << "\nDcauchyDgrad_u[_qp]:\n";
//    for (int i=0; i<9; i++){
//        for (int j=0; j<9; j++){
//            std::cout << _DcauchyDgrad_u[_qp][i][j] << " ";
//        }
//        std::cout << "\n";
//    }
//    std::cout << "\nDcauchyDphi[_qp}:\n";
//    for (int i=0; i<9; i++){
//        for (int j=0; j<9; j++){
//            std::cout << _DcauchyDphi[_qp][i][j] << " ";
//        }
//        std::cout << "\n";
//    }
//
//    std::cout << "\nDcauchyDgrad_phi[_qp]:\n";
//    for (int i=0; i<9; i++){
//        for (int j=0; j<27; j++){
//            std::cout << _DcauchyDgrad_phi[_qp][i][j] << " ";
//        }
//        std::cout << "\n";
//    }
//    std::cout << "dfdUint: " << dfdUint;
//    mooseError("fail");
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
    double dNdx[3];
    double detadx[3];
    for (int indx=0; indx<3; indx++){
        dNdx[indx]   = _grad_test[_i][_qp](indx);
        detadx[indx] = _grad_phi[_i][_qp](indx);
    }
    if(_off_diag_dof_num>=0){
        balance_equations::compute_internal_force_jacobian(_component,           _off_diag_dof_num, 
                                                           _test[_i][_qp],       dNdx, _phi[_j][_qp],          detadx,
                                                           _DcauchyDgrad_u[_qp], _DcauchyDphi[_qp],   _DcauchyDgrad_phi[_qp], dfdUint);
        return dfdUint;
    }
    else{
        return 0;
    }
}
