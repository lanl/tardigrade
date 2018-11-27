###############################################################################
#                           overlapping_columns.i                             #
###############################################################################
# A uniaxial stretch type problem where a micromorphic and DNS domain overlap #
# eachother. The test detects whether the nodes are being assigned to the     #
# correct element domain and whether the overlap coupling is being performed  #
# correctly for a simple elastic stretch problem.                             #
###############################################################################

[Mesh]
  file = overlapping_columns.e
  displacements = 'disp_x disp_y disp_z'

  #Define names for blocks
  block_id = '1 2'
  block_name = 'micro DNS'

  #Define names for boundaries
  boundary_id   = '1 2 3 4 5 6 7 8 9 10 11 12'
  boundary_name = 'M_z- M_z+ M_y- M_x+ M_x- M_y+ D_z- D_z+ D_y- D_x+ D_x- D_y+'
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
  [./phi_xx]
    block = 'micro'
  [../]
  [./phi_yy]
    block = 'micro'
  [../]
  [./phi_zz]
    block = 'micro'
  [../]
  [./phi_yz]
    block = 'micro'
  [../]
  [./phi_xz]
    block = 'micro'
  [../]
  [./phi_xy]
    block = 'micro'
  [../]
  [./phi_zy]
    block = 'micro'
  [../]
  [./phi_zx]
    block = 'micro'
  [../]
  [./phi_yx]
    block = 'micro'
  [../]
[]

[Kernels]
  #Define the internal force balance equations for micromorphic
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

    #Define the block
    block = 'micro'
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

    #Define the block
    block = 'micro'
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

    #Define the block
    block = 'micro'
  [../]
  #Define the internal couple balance equations for micromorphic
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

    #Define the block
    block = 'micro'
  [../]
  [./couple_22]
    type = InternalCouple
    component_i = 1
    component_j = 1
    dof_num     = 4
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

    #Define the block
    block = 'micro'
  [../]
  [./couple_33]
    type = InternalCouple
    component_i = 2
    component_j = 2
    dof_num     = 5
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

    #Define the block
    block = 'micro'
  [../]
  [./couple_23]
    type = InternalCouple
    component_i = 1
    component_j = 2
    dof_num     = 6
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

    #Define the block
    block = 'micro'
  [../]
  [./couple_13]
    type = InternalCouple
    component_i = 0
    component_j = 2
    dof_num     = 7
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

    #Define the block
    block = 'micro'
  [../]
  [./couple_12]
    type = InternalCouple
    component_i = 0
    component_j = 1
    dof_num     = 8
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

    #Define the block
    block = 'micro'
  [../]
  [./couple_32]
    type = InternalCouple
    component_i = 2
    component_j = 1
    dof_num     = 9
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

    #Define the block
    block = 'micro'
  [../]
  [./couple_31]
    type = InternalCouple
    component_i = 2
    component_j = 0
    dof_num     = 10
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

    #Define the block
    block = 'micro'
  [../]
  [./couple_21]
    type = InternalCouple
    component_i = 1
    component_j = 0
    dof_num     = 11
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

    #Define the block
    block = 'micro'
  [../]
  [./TensorMechanics]
    #Stress divergence kernels
    type = StressDivergenceTensors
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
#    displacements = 'disp_x disp_y disp_z'
    use_displaced_mesh = true
    block = 'DNS'
  [../]
[]

#[Modules/TensorMechanics/Master]
#  [./block2]
#    strain = FINITE
#    block = 'DNS'
#  [../]
#[]

[BCs]
  active = 'left_x back_z bottom_y top_y'
  [./left_x]
    #type = DirichletBC
    type = PresetBC
    variable = disp_x
    boundary = 'M_x- D_x-'
    value = 0
  [../]
  [./back_z]
    #type = DirichletBC
    type = PresetBC
    variable = disp_z
    boundary = 'M_z- D_z-'
    value = 0
  [../]
  [./bottom_y]
    #type = DirichletBC
    type = PresetBC
    variable = disp_y
    boundary = 'M_y-'
    value = 0
  [../]
  [./top_y]
    #type = DirichletBC
    #type = PresetBC
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'D_y+'
    function = top_bc
  [../]
[]

[Functions]
  [./top_bc]
    type  = ParsedFunction
    value = 0.1*t
  [../]
[]

[Materials]
  [./linear_elastic]
    type = MicromorphicMaterial
#    material_fparameters = '0. 0. 0. 0. 0. 0. 0. 15.4 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
#    material_fparameters = '8e3 11e3 2e3 1.538e3 -1e3 -1.39e3 -2.11e3 0. 0. 0. 0. 0. 0. 0.769 0. 0. 0. 0.'
#    material_fparameters = '8e3 11e3 2e3 1.538e3 -1e3 -1.39e3 -2.11e3 0.12 0.51 0.72 0.84 0.443 0.62 0.769 0.945 0.47 0.63 0.58'
#    material_fparameters = '29.48e3 25.48e3 1e3 0.4e3 -1.5e3 -1.4e3 -3e3 0 0 0 0 0 0 10e5 0. 0. 0. 0.'
    material_fparameters = '29. 7. 60. 10. 10. 8. 5. 0. 0. 0. 0. 0. 0. 8. 0. 0. 0. 0.'
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

    block = 'micro'
  [../]
  [./TensorMechanics/elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 68e3
    poissons_ratio = 0.32
    block = 'DNS'
  [../]
#  [./TensorMechanics/strain]
#    type = ComputeFiniteStrain
#    block = 'DNS'
#  [../]
  [./TensorMechanics/stress]
    type = ComputeFiniteStrainElasticStress
    block = 'DNS'
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
  num_steps = 10
  dt        = 0.1
  solve_type = 'PJFNK'
#  solve_type = 'NEWTON'
#  nl_rel_tol = 1e-8
#  nl_abs_tol = 1e-8
#  nl_max_its = 100
  #Terms for debugging
#  petsc_options = '-ksp_monitor_true_residual -ksp_compute_singularvalues' 
#  petsc_options = '-snes_converged_reason -ksp_converged_reason'
#  l_max_its  = 10
#  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
#  petsc_options_value = 'hypre    boomeramg      100'
#  petsc_options_iname = '-ksp_gmres_restart'
#  petsc_options_value = '100'
#  petsc_options = '-snes_ksp_ew -ksp_monitor_true_residual -ksp_compute_singularvalues'# -pc_svd_monitor'
  petsc_options = '-ksp_monitor_true_residual -ksp_compute_singularvalues'# -pc_svd_monitor'
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart -print_linear_residuals'# -ksp_view_mat'
  petsc_options_value = 'asm      lu           1               101                false                  '# binary'
[]

[Outputs]
  exodus = true
  print_perf_log = true
[]
