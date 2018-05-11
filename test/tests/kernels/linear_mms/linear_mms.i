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
    phi_32 = phi_zy
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
  [./all_phi_xx]
    type = FunctionDirichletBC
    variable = phi_xx
    boundary = 'left right bottom top back front'
    function = phi_11_fxn
  [../]
  [./all_phi_yy]
    type = FunctionDirichletBC
    variable = phi_yy
    boundary = 'left right bottom top back front'
    function = phi_22_fxn
  [../]
  [./all_phi_zz]
    type = FunctionDirichletBC
    variable = phi_zz
    boundary = 'left right bottom top back front'
    function = phi_33_fxn
  [../]
  [./all_phi_yz]
    type = FunctionDirichletBC
    variable = phi_yz
    boundary = 'left right bottom top back front'
    function = phi_23_fxn
  [../]
  [./all_phi_xz]
    type = FunctionDirichletBC
    variable = phi_xz
    boundary = 'left right bottom top back front'
    function = phi_13_fxn
  [../]
  [./all_phi_xy]
     type = FunctionDirichletBC
     variable = phi_xy
     boundary = 'left right bottom top back front'
     function = phi_12_fxn
  [../]
  [./all_phi_zy]
     type = FunctionDirichletBC
     variable = phi_zy
     boundary = 'left right bottom top back front'
     function = phi_32_fxn
  [../]
  [./all_phi_zx]
     type = FunctionDirichletBC
     variable = phi_zx
     boundary = 'left right bottom top back front'
     function = phi_31_fxn
  [../]
  [./all_phi_yx]
     type = FunctionDirichletBC
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
    material_fparameters = '29. 7. 60. 10. 10. 8. 5. 0. 0. 0. 0. 0. 0. 8. 0. 0. 0. 0.'
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
    value = (-0.392938371196*x+0.427721330099*y+0.546297092872*z+-0.102629538166)*t
  [../]
  [./u2_fxn]
    type  = ParsedFunction
    value = (-0.438937939571*x+0.153787079751*y+-0.961528396769*z+-0.36965947717)*t
  [../]
  [./u3_fxn]
    type  = ParsedFunction
    value = (0.0381361970313*x+0.215764963612*y+0.313643967698*z+-0.458099414768)*t
  [../]
  [./phi_11_fxn]
    type = ParsedFunction
    value = (0.273269772553*x+0.234751402616*y+0.348710677051*z+0.192506607384)*t
  [../]
  [./phi_22_fxn]
    type = ParsedFunction
    value = (0.810451063497*x+1.3773727073*y+1.3320544823*z+0.885074307114)*t
  [../]
  [./phi_33_fxn]
    type = ParsedFunction
    value = (0.153709634128*x+0.207560231216*y+0.223427008827*z+0.291508288536)*t
  [../]
  [./phi_23_fxn]
    type = ParsedFunction
    value = (0.408060510166*x+0.551020333375*y+0.593142646679*z+0.773881361519)*t
  [../]
  [./phi_13_fxn]
    type = ParsedFunction
    value = (-0.718279686634*x+-0.969921623204*y+-1.04406651409*z+-1.36220792749)*t
  [../]
  [./phi_12_fxn]
    type = ParsedFunction
    value = (0.453315731782*x+0.770416308732*y+0.745068122696*z+0.495055315838)*t
  [../]
  [./phi_32_fxn]
    type = ParsedFunction
    value = (-0.0332158273198*x+-0.0564507544767*y+-0.0545934155157*z+-0.036274214045)*t
  [../]
  [./phi_31_fxn]
    type = ParsedFunction
    value = (1.51763494893*x+1.30371877428*y+1.93660464388*z+1.06910747038)*t
  [../]
  [./phi_21_fxn]
    type = ParsedFunction
    value = (-0.0906507373851*x+-0.0778731857187*y+-0.115676460347*z+-0.0638594812293)*t
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
  type = Transient
  end_time = 1
  dtmin   = 1e-3
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
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt   = 0.5
    cutback_factor     = 0.4
    growth_factor      = 1.2
    optimal_iterations = 5
  [../]
  petsc_options = '-snes_ksp_ew -snes_converged_reason -ksp_converged_reason -ksp_monitor_true_residual -ksp_compute_singular_values'
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart'
  petsc_options_value = 'asm      lu           1               101               '
[]

[Outputs]
  exodus = true
[]
