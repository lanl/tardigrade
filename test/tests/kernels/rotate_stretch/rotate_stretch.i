###############################################################################
#                             rotate_stretch.i                                #
###############################################################################
# A test which rotates and stretches a material at the same time to make sure #
# that the stress is invariant w.r.t. rotations.                              #
###############################################################################

[Mesh]
  type = GeneratedMesh
  displacements = 'disp_x disp_y disp_z'
  dim = 3
  nx = 1
  ny = 1
  nz = 1
[]

[Variables]
  [./disp_x]
#    [./InitialCondition]
#      type = RandomIC
#      min  = -0.1
#      max  =  0.1
#    [../]
#    scaling = 1e-4
  [../]
  [./disp_y]
#    [./InitialCondition]
#      type = RandomIC
#      min  = -0.1
#      max  =  0.1
#    [../]
#    scaling = 1e-4
  [../]
  [./disp_z]
#    [./InitialCondition]
#      type = RandomIC
#      min  = -0.1
#      max  =  0.1
#    [../]
#    scaling = 1e-4
  [../]
  [./phi_xx]
#    order = CONSTANT
#    family = MONOMIAL
#    [./InitialCondition]
#      type = RandomIC
#      min  = -0.1
#      max  =  0.1
#    [../]
#    scaling = 1e-5
  [../]
  [./phi_yy]
#    order = CONSTANT
#    family = MONOMIAL
#    [./InitialCondition]
#      type = RandomIC
#      min  = -0.1
#      max  =  0.1
#    [../]
#    scaling = 1e-5
  [../]
  [./phi_zz]
#    order = CONSTANT
#    family = MONOMIAL
#    [./InitialCondition]
#      type = RandomIC
#      min  = -0.1
#      max  =  0.1
#    [../]
#    scaling = 1e-5
  [../]
  [./phi_yz]
#    order = CONSTANT
#    family = MONOMIAL
#    [./InitialCondition]
#      type = RandomIC
#      min  = -0.1
#      max  =  0.1
#    [../]
#    scaling = 1e-1
  [../]
  [./phi_xz]
#    order = CONSTANT
#    family = MONOMIAL
#    [./InitialCondition]
#      type = RandomIC
#      min  = -0.1
#      max  =  0.1
#    [../]
#    scaling = 1e-1
  [../]
  [./phi_xy]
#    order = CONSTANT
#    family = MONOMIAL
#    [./InitialCondition]
#      type = RandomIC
#      min  = -0.1
#      max  =  0.1
#    [../]
#    scaling = 1e-1
  [../]
  [./phi_zy]
#    order = CONSTANT
#    family = MONOMIAL
#    [./InitialCondition]
#      type = RandomIC
#      min  = -0.1
#      max  =  0.1
#    [../]
#    scaling = 1e-1
  [../]
  [./phi_zx]
#    order = CONSTANT
#    family = MONOMIAL
#    [./InitialCondition]
#      type = RandomIC
#      min  = -0.1
#      max  =  0.1
#    [../]
#    scaling = 1e-1
  [../]
  [./phi_yx]
#    order = CONSTANT
#    family = MONOMIAL
#    [./InitialCondition]
#      type = RandomIC
#      min  = -0.1
#      max  =  0.1
#    [../]
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
  [./couple_32]
    type = InternalCouple
    component_i = 2
    component_j = 1
    dof_num     = 10
    variable    = phi_zy

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
  active = 'left_x left_y back_z bottom_x bottom_y top_x top_y'
  [./left_x]
    type     = FunctionDirichletBC
    variable = disp_x
    boundary = 'left'
    function = moving_x
  [../]
  [./left_y]
    type     = FunctionDirichletBC
    variable = disp_y
    boundary = 'left'
    function = moving_y
  [../]
  [./right_x]
    type     = FunctionDirichletBC
    variable = disp_x
    boundary = 'right'
    function = moving_x
  [../]
  [./right_y]
    type     = FunctionDirichletBC
    variable = disp_y
    boundary = 'right'
    function = moving_y
  [../]
  [./back_z]
    type = PresetBC
    variable = disp_z
    boundary = 'back'
    value = 0
  [../]
  [./front_z]
    type = PresetBC
    variable = disp_z
    boundary = 'front'
    value = 0
  [../]
  [./bottom_x]
    type     = FunctionDirichletBC
    variable = disp_x
    boundary = 'bottom'
    function = fixed_x
  [../]
  [./bottom_y]
    type     = FunctionDirichletBC
    variable = disp_y
    boundary = 'bottom'
    function = fixed_y
  [../]
  [./top_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'top'
    function = moving_x
  [../]
  [./top_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'top'
    function = moving_y
  [../]
[]

[Functions]
  [./fixed_x]
    type  = ParsedFunction
    value = 'x*(cos(pi*t)-1)-y*sin(pi*t)'
  [../]
  [./fixed_y]
    type  = ParsedFunction
    value = 'x*sin(pi*t)+y*(cos(pi*t)-1)'
  [../]
  [./moving_x]
    type  = ParsedFunction
    value = 'x*(cos(pi*t)-1)-(1+0.05*t)*y*sin(pi*t)'
  [../]
  [./moving_y]
    type  = ParsedFunction
    value = 'x*sin(pi*t)+(1+0.05*t)*y*(cos(pi*t)-1)'
  [../]
[]

[Materials]
  [./linear_elastic]
    type = MicromorphicMaterial
#    material_fparameters = '0. 0. 0. 0. 0. 0. 0. 15.4 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
#    material_fparameters = '8e3 11e3 2e3 1.538e3 -1e3 -1.39e3 -2.11e3 0. 0. 0. 0. 0. 0. 0.769 0. 0. 0. 0.'
#    material_fparameters = '8e3 11e3 2e3 1.538e3 -1e3 -1.39e3 -2.11e3 0.12 0.51 0.72 0.84 0.443 0.62 0.769 0.945 0.47 0.63 0.58'
#    material_fparameters = '29.48e3 25.48e3 1e3 0.4e3 -1.5e3 -1.4e3 -3e3 0 0 0 0 0 0 10e5 0. 0. 0. 0.'
#    material_fparameters = '29. 7. 60. 10. 10. 8. 5. 0. 0. 0. 0. 0. 0. 8. 0. 0. 0. 0.'
#    material_fparameters = '2. 696.47 65.84 5. -7.69 -51.92 38.61 -27.31 5.13 11. 1.85 -0.19 -1.08 -1.57 2.29 -0.61 5.97 -2.02 2.38 -0.32 -3.25 2. -51.92 5.13'
    material_fparameters = '2 29.48e3 25.48e3 5 1e3 0.4e3 -1.5e3 -1.4e3 -3e3 11 0 0 0 0 0 0 10e5 0 0 0 0 2 .4e3 -3e3'
    model_name = "LinearElasticity"

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
  end_time = 1.0
  dtmin    = 1e-4
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-9
#  solve_type = 'NEWTON'
  solve_type = 'PJFNK'
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt   = 0.1
#    dt   = 0.01
    cutback_factor     = 0.4
    growth_factor      = 1.2
    optimal_iterations = 100
  [../]
  petsc_options = '-ksp_monitor_true_residual'
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart -print_linear_residuals'
  petsc_options_value = 'asm      lu           1               101                false                  '
 
[]

[Outputs]
  exodus = true
  perf_graph = true
#  [./out]
#    type = Exodus
#    output_material_properties = true
#  [../]
[]
