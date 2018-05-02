###############################################################################
#                                const_mms.i                                  #
###############################################################################
# A constant value method of manufactured solutions solution. This should be  #
# a relatively simple solution.                                               #
###############################################################################

[Mesh]
  type = GeneratedMesh
  displacements = 'disp_x disp_y disp_z'
  dim = 3
  nx = 4
  ny = 4
  nz = 4
#  file = unit_cube.e
[]

[Variables]
  [./disp_x]
#    [./InitialCondition]
#      type = RandomIC
#      min  = -1
#      max  =  1
#    [../]
    scaling = 1e-3
  [../]
  [./disp_y]
#    [./InitialCondition]
#      type = RandomIC
#      min  = -1
#      max  =  1
#    [../]
    scaling = 1e-3
  [../]
  [./disp_z]
#    [./InitialCondition]
#      type = RandomIC
#      min  = -1
#      max  =  1
#    [../]
    scaling = 1e-3
  [../]
  [./phi_xx]
#    [./InitialCondition]
#      type = RandomIC
#      min  = -1
#      max  =  1
#    [../]
    scaling = 1e-3
  [../]
  [./phi_yy]
#    [./InitialCondition]
#      type = RandomIC
#      min  = -1
#      max  =  1
#    [../]
    scaling = 1e-3
  [../]
  [./phi_zz]
#    [./InitialCondition]
#      type = RandomIC
#      min  = -1
#      max  =  1
#    [../]
    scaling = 1e-3
  [../]
  [./phi_yz]
#    [./InitialCondition]
#      type = RandomIC
#      min  = -1
#      max  =  1
#    [../]
    scaling = 1e-3
  [../]
  [./phi_xz]
#    [./InitialCondition]
#      type = RandomIC
#      min  = -1
#      max  =  1
#    [../]
    scaling = 1e-3
  [../]
  [./phi_xy]
#    [./InitialCondition]
#      type = RandomIC
#      min  = -1
#      max  =  1
#    [../]
    scaling = 1e-3
  [../]
  [./phi_zy]
#    [./InitialCondition]
#      type = RandomIC
#      min  = -1
#      max  =  1
#    [../]
    scaling = 1e-3
  [../]
  [./phi_zx]
#    [./InitialCondition]
#      type = RandomIC
#      min  = -1
#      max  =  1
#    [../]
    scaling = 1e-3
  [../]
  [./phi_yx]
#    [./InitialCondition]
#      type = RandomIC
#      min  = -1
#      max  =  1
#    [../]
    scaling = 1e-3
  [../]
[]

[Kernels]
  #Define the internal force balance equations
  [./force_1]
    type = InternalForce
    component = 0
    dof_num   = 0
    variable  = disp_x
    MMS       = true

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./force_2]
    type = InternalForce
    component = 1
    dof_num   = 1
    variable  = disp_y
    MMS       = true

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./force_3]
    type = InternalForce
    component = 2
    dof_num   = 2
    variable  = disp_z
    MMS       = true

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  #Define the internal couple balance equations
  [./couple_11]
    type = InternalCouple
    component_i = 0
    component_j = 0
    dof_num     = 3
    variable    = phi_xx
    MMS       = true

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_22]
    type = InternalCouple
    component_i = 1
    component_j = 1
    dof_num     = 4
    variable    = phi_yy
    MMS       = true

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_33]
    type = InternalCouple
    component_i = 2
    component_j = 2
    dof_num     = 5
    variable    = phi_zz
    MMS       = true

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_23]
    type = InternalCouple
    component_i = 1
    component_j = 2
    dof_num     = 6
    variable    = phi_yz
    MMS       = true

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_13]
    type = InternalCouple
    component_i = 0
    component_j = 2
    dof_num     = 7
    variable    = phi_xz
    MMS       = true

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_12]
    type = InternalCouple
    component_i = 0
    component_j = 1
    dof_num     = 8
    variable    = phi_xy
    MMS       = true

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_32]
    type = InternalCouple
    component_i = 2
    component_j = 1
    dof_num     = 9
    variable    = phi_zy
    MMS       = true

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_31]
    type = InternalCouple
    component_i = 2
    component_j = 0
    dof_num     = 10
    variable    = phi_zx
    MMS       = true

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_21]
    type = InternalCouple
    component_i = 1
    component_j = 0
    dof_num     = 11
    variable    = phi_yx
    MMS       = true

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
[]

[BCs]
  [./all_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'left right bottom top back front'
    function = u1_fxn
  [../]
  [./all_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'left right bottom top back front'
    function = u2_fxn
  [../]
  [./all_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = 'left right bottom top back front'
    function = u3_fxn
  [../]
[]

[Materials]
  [./linear_elastic]
    type = MicromorphicMaterial
    material_fparameters = '8e3 11e3 2e3 1.538e3 -1e3 -1.39e3 -2.11e3 0. 0. 0. 0. 0. 0. 0.769 0. 0. 0. 0.'
    model_name = "LinearElasticity"
    MMS       = true

    #Coupled variables
    u1     = 'disp_x'
    u2     = 'disp_y'
    u3     = 'disp_z'
    phi_11 = 'phi_xx'
    phi_22 = 'phi_yy'
    phi_33 = 'phi_zz'
    phi_23 = 'phi_yz'
    phi_13 = 'phi_xz'
    phi_12 = 'phi_xy'
    phi_32 = 'phi_zy'
    phi_31 = 'phi_zx'
    phi_21 = 'phi_yx'

    #The manufactured solutions function
    u1_fxn     = 'u1_fxn'
    u2_fxn     = 'u2_fxn'
    u3_fxn     = 'u3_fxn'
    phi_11_fxn = 'phi_11_fxn'
    phi_22_fxn = 'phi_22_fxn'
#    phi_33_fxn = 'phi_33_fxn'
#    phi_23_fxn = 'phi_23_fxn'
#    phi_13_fxn = 'phi_13_fxn'
#    phi_12_fxn = 'phi_12_fxn'
#    phi_32_fxn = 'phi_32_fxn'
#    phi_31_fxn = 'phi_31_fxn'
#    phi_21_fxn = 'phi_21_fxn'
  [../]
[]

[Functions]
  [./u1_fxn]
    type  = ParsedFunction
    value = 1.721
  [../]
  [./u2_fxn]
    type  = ParsedFunction
    value = 0.321
  [../]
  [./u3_fxn]
    type  = ParsedFunction
    value = -0.981
  [../]
  [./phi_11_fxn]
    type = ParsedFunction
    value = -.2152
  [../]
  [./phi_22_fxn]
    type = ParsedFunction
    value = 1.31
  [../]
[]

#[Preconditioning]
#  [./SMP]
#    type = SMP
#    full = true
#  [../]
#[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK'
#  solve_type = 'NEWTON'
#  nl_rel_tol = 1e-8
#  nl_abs_tol = 1e-8
#  nl_max_its = 100
  #Terms for debugging
  petsc_options = 'snes_type_test -snes_test_display' 
#  petsc_options = '-snes_converged_reason -ksp_converged_reason'
#  l_max_its  = 10
#  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
#  petsc_options_value = 'hypre    boomeramg      100'
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = '100'
[]

[Outputs]
  exodus = true
[]
