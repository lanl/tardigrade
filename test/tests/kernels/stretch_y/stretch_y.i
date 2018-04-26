[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 10
  ny = 10
  nz = 10
[]

[Variables]
  [./u_x]
  [../]
  [./u_y]
  [../]
  [./u_z]
  [../]
  [./phi_xx]
  [../]
  [./phi_yy]
  [../]
  [./phi_zz]
  [../]
  [./phi_yz]
  [../]
  [./phi_xz]
  [../]
  [./phi_xy]
  [../]
  [./phi_zy]
  [../]
  [./phi_zx]
  [../]
  [./phi_yx]
  [../]
[]

[Kernels]
  #Define the internal force balance equations
  [./force_1]
    type = InternalForce
    component = 0
    dof_num   = 0
    variable = u_x
  [../]
  [./force_2]
    type = InternalForce
    component = 1
    dof_num   = 1
    variable = u_y
  [../]
  [./force_3]
    type = InternalForce
    component = 2
    dof_num   = 2
    variable = u_z
  [../]
  #Define the internal couple balance equations
  [./couple_11]
    type = InternalForce
    component_i = 0
    component_j = 0
    dof_num   = 3
    variable = phi_xx
  [../]
  [./couple_22]
    type = InternalForce
    component_i = 1
    component_j = 1
    dof_num   = 4
    variable = phi_yy
  [../]
  [./couple_33]
    type = InternalForce
    component_i = 2
    component_j = 2
    dof_num   = 5
    variable = phi_zz
  [../]
  [./couple_23]
    type = InternalForce
    component_i = 1
    component_j = 2
    dof_num   = 6
    variable = phi_yz
  [../]
  [./couple_13]
    type = InternalForce
    component_i = 0
    component_j = 2
    dof_num   = 7
    variable = phi_xz
  [../]
  [./couple_12]
    type = InternalForce
    component_i = 0
    component_j = 1
    dof_num   = 7
    variable = phi_xy
  [../]
  [./couple_32]
    type = InternalForce
    component_i = 2
    component_j = 1
    dof_num   = 8
    variable = phi_zy
  [../]
  [./couple_31]
    type = InternalForce
    component_i = 2
    component_j = 0
    dof_num   = 9
    variable = phi_zx
  [../]
  [./couple_21]
    type = InternalForce
    component_i = 1
    component_j = 0
    dof_num   = 7
    variable = phi_yx
  [../]
[]

[BCs]
  [./left]
    type = DirichletBC
    variable = u_x
    boundary = left
    value = 0
  [../]
  [./back]
    type = DirichletBC
    variable = u_y
    boundary = back
    value = 0
  [../]
  [./bottom]
    type = DirichletBC
    variable = u_z
    boundary = bottom
    value = 0
  [../]
  [./top]
    type = DirichletBC
    variable = u_y
    boundary = front
    value = 1
  [../]
[]

[Materials]
  [./linear_elastic]
    type = MicromorphicMaterial
    material_fparameters = '1000., 8e9, 11e9, 2e9, 1.538e9, -1e9, -1.39e9, -2.11e9, 0., 0., 0., 0., 0., 0., 0.769e6, 0., 0., 0., 0.'
    model_name = "LinearElasticity"

    #Coupled variables
    u1     = 'u_x'
    u2     = 'u_y'
    u3     = 'u_z'
    phi_11 = 'phi_11'
    phi_22 = 'phi_22'
    phi_33 = 'phi_33'
    phi_23 = 'phi_23'
    phi_13 = 'phi_13'
    phi_12 = 'phi_12'
    phi_32 = 'phi_32'
    phi_31 = 'phi_31'
    phi_21 = 'phi_21'
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
[]
