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
  nx = 2
  ny = 2
  nz = 2
#  file = unit_cube.e
[]

[Variables]
  [./disp_x]
#    scaling = 1e-1
  [../]
  [./disp_y]
#    scaling = 1e-1
  [../]
  [./disp_z]
#    scaling = 1e-1
  [../]
  [./phi_xx]
#    scaling = 1e-1
  [../]
  [./phi_yy]
#    scaling = 1e-1
  [../]
  [./phi_zz]
#    scaling = 1e-1
  [../]
  [./phi_yz]
#    scaling = 1e-1
  [../]
  [./phi_xz]
#    scaling = 1e-1
  [../]
  [./phi_xy]
#    scaling = 1e-1
  [../]
  [./phi_zy]
#    scaling = 1e-1
  [../]
  [./phi_zx]
#    scaling = 1e-1
  [../]
  [./phi_yx]
#    scaling = 1e-1
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
    dof_num     = 7
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
    dof_num     = 11
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
    dof_num     = 8
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
    dof_num     = 5
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
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_12]
    type = InternalCouple
    component_i = 0
    component_j = 1
    dof_num     = 4
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
    dof_num     = 10
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
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_31]
    type = InternalCouple
    component_i = 2
    component_j = 0
    dof_num     = 9
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
    dof_num     = 6
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
    preset = true
    variable = disp_x
    boundary = 'left right bottom top back front'
    function = u1_fxn
  [../]
  [./all_y]
    type = FunctionDirichletBC
    preset = true
    variable = disp_y
    boundary = 'left right bottom top back front'
    function = u2_fxn
  [../]
  [./all_z]
    type = FunctionDirichletBC
    preset = true
    variable = disp_z
    boundary = 'left right bottom top back front'
    function = u3_fxn
  [../]
  [./all_phi_xx]
    type = FunctionDirichletBC
    preset = true
    variable = phi_xx
    boundary = 'left right bottom top back front'
    function = phi_11_fxn
  [../]
  [./all_phi_yy]
    type = FunctionDirichletBC
    preset = true
    variable = phi_yy
    boundary = 'left right bottom top back front'
    function = phi_22_fxn
  [../]
  [./all_phi_zz]
    type = FunctionDirichletBC
    preset = true
    variable = phi_zz
    boundary = 'left right bottom top back front'
    function = phi_33_fxn
  [../]
  [./all_phi_yz]
    type = FunctionDirichletBC
    preset = true
    variable = phi_yz
    boundary = 'left right bottom top back front'
    function = phi_23_fxn
  [../]
  [./all_phi_xz]
    type = FunctionDirichletBC
    preset = true
    variable = phi_xz
    boundary = 'left right bottom top back front'
    function = phi_13_fxn
  [../]
  [./all_phi_xy]
     type = FunctionDirichletBC
     preset = true
     variable = phi_xy
     boundary = 'left right bottom top back front'
     function = phi_12_fxn
  [../]
  [./all_phi_zy]
     type = FunctionDirichletBC
     preset = true
     variable = phi_zy
     boundary = 'left right bottom top back front'
     function = phi_32_fxn
  [../]
  [./all_phi_zx]
     type = FunctionDirichletBC
     preset = true
     variable = phi_zx
     boundary = 'left right bottom top back front'
     function = phi_31_fxn
  [../]
  [./all_phi_yx]
     type = FunctionDirichletBC
     preset = true
     variable = phi_yx
     boundary = 'left right bottom top back front'
     function = phi_21_fxn
  [../]
[]

[Materials]
  [./linear_elastic]
    type = MicromorphicMaterial
#    material_fparameters = '8e3 11e3 2e3 1.538e3 -1e3 -1.39e3 -2.11e3 0. 0. 0. 0. 0. 0. 0.769 0. 0. 0. 0.'
#    material_fparameters = '8e3 11e3 2e3 1.538e3 -1e3 -1.39e3 -2.11e3 0.12 0.51 0.72 0.84 0.443 0.62 0.769 0.945 0.47 0.63 0.58'
#    material_fparameters = '29. 7. 60. 10. 10. 8. 5. 0. 0. 0. 0. 0. 0. 8. 0. 0. 0. 0.'
    material_fparameters = '2. 696.47 65.84 5. -7.69 -51.92 38.61 -27.31 5.13 11. 1.85 -0.19 -1.08 -1.57 2.29 -0.61 5.97 -2.02 2.38 -0.32 -3.25 2. -51.92 5.13'
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
    phi_33_fxn = 'phi_33_fxn'
    phi_23_fxn = 'phi_23_fxn'
    phi_13_fxn = 'phi_13_fxn'
    phi_12_fxn = 'phi_12_fxn'
    phi_32_fxn = 'phi_32_fxn'
    phi_31_fxn = 'phi_31_fxn'
    phi_21_fxn = 'phi_21_fxn'
  [../]
[]

# Functions for the method of manufactured solutions
[Functions]
  [./u1_fxn]
    type  = ParsedFunction
    value = 1.721*t
  [../]
  [./u2_fxn]
    type  = ParsedFunction
    value = 0.321*t
  [../]
  [./u3_fxn]
    type  = ParsedFunction
    value = -0.981*t
  [../]
  [./phi_11_fxn]
    type = ParsedFunction
#    value = -.2152
    value = 0.054142*t
  [../]
  [./phi_22_fxn]
    type = ParsedFunction
#    value = 1.31
    value = 0.07059678*t
  [../]
  [./phi_33_fxn]
    type = ParsedFunction
    #value = 2.142
#    value = -.521
    value = 0.04161017*t
  [../]
  [./phi_23_fxn]
    type = ParsedFunction
#    value = -0.177
    value = -0.00516283*t
  [../]
  [./phi_13_fxn]
    type = ParsedFunction
#    value = 0.606
    value = -0.0056683*t
  [../]
  [./phi_12_fxn]
    type = ParsedFunction
#    value = 3.72
    value = 0.00955174*t
  [../]
  [./phi_32_fxn]
    type = ParsedFunction
#    value = .827
    value = 0.01006055*t
  [../]
  [./phi_31_fxn]
    type = ParsedFunction
#    value = .718
    value = 0.00635417*t
  [../]
  [./phi_21_fxn]
    type = ParsedFunction
#    value = 2.271
    value = -0.00823092*t
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    #type = FDP
    full = true
  [../]
[]

[Executioner]
#  type = Steady
  type = Transient
  num_steps = 2
  dt = 0.5
  solve_type = 'PJFNK'
#  solve_type = 'NEWTON'
#  nl_rel_tol = 1e-8
#  nl_abs_tol = 1e-8
#  nl_max_its = 100
  #Terms for debugging
#  petsc_options = 'snes_type_test -snes_test_display' 
#  petsc_options = '-snes_converged_reason -ksp_converged_reason'
#  l_max_its  = 10
#  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
#  petsc_options_value = 'hypre    boomeramg      100'
#  petsc_options_iname = '-ksp_gmres_restart'
#  petsc_options_value = '100'
  petsc_options = '-snes_ksp_ew -snes_converged_reason -ksp_converged_reason -ksp_monitor_true_residual -ksp_compute_singular_values'
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart'
  petsc_options_value = 'asm      lu           1               101               '
[]

[Outputs]
  exodus = true
  perf_graph = true
[]
