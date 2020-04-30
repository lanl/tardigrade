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
    type = InternalForceOverlap
    nodal_overlap_userobject = nodal_overlap
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
    type = InternalForceOverlap
    nodal_overlap_userobject = nodal_overlap
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
    type = InternalForceOverlap
    nodal_overlap_userobject = nodal_overlap
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
    type = InternalCoupleOverlap
    nodal_overlap_userobject = nodal_overlap
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
    type = InternalCoupleOverlap
    nodal_overlap_userobject = nodal_overlap
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
    type = InternalCoupleOverlap
    nodal_overlap_userobject = nodal_overlap
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
    type = InternalCoupleOverlap
    nodal_overlap_userobject = nodal_overlap
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
    type = InternalCoupleOverlap
    nodal_overlap_userobject = nodal_overlap
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
    type = InternalCoupleOverlap
    nodal_overlap_userobject = nodal_overlap
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
    type = InternalCoupleOverlap
    nodal_overlap_userobject = nodal_overlap
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
    type = InternalCoupleOverlap
    nodal_overlap_userobject = nodal_overlap
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
    type = InternalCoupleOverlap
    nodal_overlap_userobject = nodal_overlap
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
  [./DNS_force_1]
    type = StressDivergenceTensorsOverlap
    nodal_overlap_userobject = nodal_overlap
    variable = 'disp_x'
    component = 0
    displacements = 'disp_x disp_y disp_z'
    use_finite_deform_jacobian = True
    use_displaced_mesh = True
    block = 'DNS'
  [../]
  [./DNS_force_2]
    type = StressDivergenceTensorsOverlap
    nodal_overlap_userobject = nodal_overlap
    variable = 'disp_y'
    component = 1
    displacements = 'disp_x disp_y disp_z'
    use_finite_deform_jacobian = True
    use_displaced_mesh = True
    block = 'DNS'
  [../]
  [./DNS_force_3]
    type = StressDivergenceTensorsOverlap
    nodal_overlap_userobject = nodal_overlap
    variable = 'disp_z'
    component = 2
    displacements = 'disp_x disp_y disp_z'
    use_finite_deform_jacobian = True
    use_displaced_mesh = True
    block = 'DNS'
  [../]
#  [./TensorMechanics]
#    strain = FINITE
#    block = 'DNS'
#  [../]
[]

[NodalKernels]
  [./DNS_coupling_u1]
    type = ProjectedDOF
    variable = disp_x
    block = 'DNS'
    nodal_overlap_userobject = nodal_overlap
    DNS_dof_userobject = DNS_dof
#    macro_dof_userobject = micromorphic_dof
    dof_num = 0
    is_DNS = true
#    scale_factor= 500
  [../]
  [./DNS_coupling_u2]
    type = ProjectedDOF
    variable = disp_y
    block = 'DNS'
    nodal_overlap_userobject = nodal_overlap
    DNS_dof_userobject = DNS_dof
#    macro_dof_userobject = micromorphic_dof
    dof_num = 1
    is_DNS = true
#    scale_factor= 500
  [../]
  [./DNS_coupling_u3]
    type = ProjectedDOF
    variable = disp_z
    block = 'DNS'
    nodal_overlap_userobject = nodal_overlap
    DNS_dof_userobject = DNS_dof
#    macro_dof_userobject = micromorphic_dof
    dof_num = 2
    is_DNS = true
#    scale_factor= 500
  [../]
  [./micromorphic_coupling_u1]
    type = ProjectedDOF
    variable = disp_x
    block = 'micro'
    nodal_overlap_userobject = nodal_overlap
    DNS_dof_userobject = DNS_dof
#    macro_dof_userobject = micromorphic_dof
    dof_num = 0
    is_DNS = false
#    scale_factor= 500
  [../]
  [./micromorphic_coupling_u2]
    type = ProjectedDOF
    variable = disp_y
    block = 'micro'
    nodal_overlap_userobject = nodal_overlap
    DNS_dof_userobject = DNS_dof
#    macro_dof_userobject = micromorphic_dof
    dof_num = 1
    is_DNS = false
#    scale_factor= 500
  [../]
  [./micromorphic_coupling_u3]
    type = ProjectedDOF
    variable = disp_z
    block = 'micro'
    nodal_overlap_userobject = nodal_overlap
    DNS_dof_userobject = DNS_dof
#    macro_dof_userobject = micromorphic_dof
    dof_num = 2
    is_DNS = false
#    scale_factor= 500
  [../]
  [./micromorphic_coupling_phi11]
    type = ProjectedDOF
    variable = phi_xx
    block = 'micro'
    nodal_overlap_userobject = nodal_overlap
    DNS_dof_userobject = DNS_dof
#    macro_dof_userobject = micromorphic_dof
    dof_num = 3
    is_DNS = false
#    scale_factor= 500
  [../]
  [./micromorphic_coupling_phi22]
    type = ProjectedDOF
    variable = phi_yy
    block = 'micro'
    nodal_overlap_userobject = nodal_overlap
    DNS_dof_userobject = DNS_dof
#    macro_dof_userobject = micromorphic_dof
    dof_num = 4
    is_DNS = false
#    scale_factor= 500
  [../]
  [./micromorphic_coupling_phi33]
    type = ProjectedDOF
    variable = phi_zz
    block = 'micro'
    nodal_overlap_userobject = nodal_overlap
    DNS_dof_userobject = DNS_dof
#    macro_dof_userobject = micromorphic_dof
    dof_num = 5
    is_DNS = false
#    scale_factor= 500
  [../]
  [./micromorphic_coupling_phi23]
    type = ProjectedDOF
    variable = phi_yz
    block = 'micro'
    nodal_overlap_userobject = nodal_overlap
    DNS_dof_userobject = DNS_dof
#    macro_dof_userobject = micromorphic_dof
    dof_num = 6
    is_DNS = false
#    scale_factor= 500
  [../]
  [./micromorphic_coupling_phi13]
    type = ProjectedDOF
    variable = phi_xz
    block = 'micro'
    nodal_overlap_userobject = nodal_overlap
    DNS_dof_userobject = DNS_dof
#    macro_dof_userobject = micromorphic_dof
    dof_num = 7
    is_DNS = false
#    scale_factor= 500
  [../]
  [./micromorphic_coupling_phi12]
    type = ProjectedDOF
    variable = phi_xy
    block = 'micro'
    nodal_overlap_userobject = nodal_overlap
    DNS_dof_userobject = DNS_dof
#    macro_dof_userobject = micromorphic_dof
    dof_num = 8
    is_DNS = false
#    scale_factor= 500
  [../]
  [./micromorphic_coupling_phi32]
    type = ProjectedDOF
    variable = phi_zy
    block = 'micro'
    nodal_overlap_userobject = nodal_overlap
    DNS_dof_userobject = DNS_dof
#    macro_dof_userobject = micromorphic_dof
    dof_num = 9
    is_DNS = false
#    scale_factor= 500
  [../]
  [./micromorphic_coupling_phi31]
    type = ProjectedDOF
    variable = phi_zx
    block = 'micro'
    nodal_overlap_userobject = nodal_overlap
    DNS_dof_userobject = DNS_dof
#    macro_dof_userobject = micromorphic_dof
    dof_num = 10
    is_DNS = false
#    scale_factor= 500
  [../]
  [./micromorphic_coupling_phi21]
    type = ProjectedDOF
    variable = phi_yx
    block = 'micro'
    nodal_overlap_userobject = nodal_overlap
    DNS_dof_userobject = DNS_dof
#    macro_dof_userobject = micromorphic_dof
    dof_num = 11
    is_DNS = false
#    scale_factor= 500
  [../]
[]

#[Modules]
#  [./TensorMechanics]
#    [./Master]
#      [./block2]
#        strain = FINITE
#        block = 'DNS'
#      [../]
#    [../]
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
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 68
    poissons_ratio = 0.32
    block = 'DNS'
  [../]
  [./strain]
    type = ComputeFiniteStrain
    block = 'DNS'
  [../]
  [./stress]
    type = ComputeFiniteStrainElasticStress
    block = 'DNS'
  [../]
[]

[AuxVariables]
  [./av1] #Dummy AuxVariable for ordering
  [../]
  [./av2] #Dummy AuxVariable for ordering
  [../]
#  [./element_aux]
#   order = CONSTANT
#   family = MONOMIAL
#   block = 'micro'
#  [../]
[]

[AuxKernels]
  [./eval_nodaloverlap]
    type = NodalUOAux
    variable = av1
#    priorvar = av2
    block = 'DNS'
    execute_on = initial
#    nodal_overlap_userobject = nodal_overlap
    nodal_userobject = nodal_overlap
  [../]
  [./eval_elementintegrate]
    type = ElementUOAux
    variable = av2
#    priorvar = av1
    block = 'micro'
    execute_on = initial
    element_integrate_userobject = element_integrate
  [../]
  [./eval_DNS_dof]
    type = NodalUOAux
    variable = av1
    block = 'DNS'
    execute_on = linear
    nodal_userobject = DNS_dof
  [../]
  [./eval_micromorphic_dof]
    type = NodalUOAux
    variable = av1
    block = 'micro'
    execute_on = linear
    nodal_userobject = micromorphic_dof
  [../]
#  [./compute_overlap]
#    type = ComputeProjectorAux
#    variable = element_aux
#    block = 'micro'
#    nodal_overlap_userobject = nodal_overlap
#    element_integrate_userobject = element_integrate
#  [../]
[]

[UserObjects]
  [./nodal_overlap]
    type = NodalOverlapUserObject
    block = 'DNS'
    execute_on = initial
    macroscale_domain = 'micro'
    unique_node_execute = true
  [../]
  [./element_integrate]
    type = ElementIntegrateUserObject
    block = 'DNS'
    variable = 'disp_x'
    execute_on = initial
  [../]
  [./projector]
    type = ProjectorUserObject
    block = 'micro'
    execute_on = initial
    nodal_overlap_userobject = nodal_overlap
    element_integrate_userobject = element_integrate
  [../]
  [./micromorphic_dof]
    type = MicromorphicDOFUserObject
    block = 'micro'
    execute_on = linear
    nodal_overlap_userobject = nodal_overlap
    projector_userobject = projector
    u1 = disp_x
    u2 = disp_y
    u3 = disp_z
    phi11 = phi_xx
    phi22 = phi_yy
    phi33 = phi_zz
    phi23 = phi_yz
    phi13 = phi_xz
    phi12 = phi_xy
    phi32 = phi_zy
    phi31 = phi_zx
    phi21 = phi_yx
  [../]
  [./DNS_dof]
    type = DNSDOFUserObject
    block = 'DNS'
    execute_on = linear
    nodal_overlap_userobject = nodal_overlap
    projector_userobject = projector
    micromorphic_DOF_userobject = micromorphic_dof
    u1 = disp_x
    u2 = disp_y
    u3 = disp_z
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
  nl_max_its = 20
  l_max_its  = 5
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
  perf_graph = true
[]
