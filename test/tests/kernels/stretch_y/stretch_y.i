[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
[]

#[Mesh]
#  file = unit_cube.e
#[]

[Variables]
  [./u_x]
  [../]
  [./u_y]
    initial_condition = 0.0
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
    variable  = u_x
  [../]
  [./force_2]
    type = InternalForce
    component = 1
    dof_num   = 1
    variable  = u_y
  [../]
  [./force_3]
    type = InternalForce
    component = 2
    dof_num   = 2
    variable  = u_z
  [../]
  #Define the internal couple balance equations
  [./couple_11]
    type = InternalCouple
    component_i = 0
    component_j = 0
    dof_num     = 3
    variable    = phi_xx
  [../]
  [./couple_22]
    type = InternalCouple
    component_i = 1
    component_j = 1
    dof_num     = 4
    variable    = phi_yy
  [../]
  [./couple_33]
    type = InternalCouple
    component_i = 2
    component_j = 2
    dof_num     = 5
    variable    = phi_zz
  [../]
  [./couple_23]
    type = InternalCouple
    component_i = 1
    component_j = 2
    dof_num     = 6
    variable    = phi_yz
  [../]
  [./couple_13]
    type = InternalCouple
    component_i = 0
    component_j = 2
    dof_num     = 7
    variable    = phi_xz
  [../]
  [./couple_12]
    type = InternalCouple
    component_i = 0
    component_j = 1
    dof_num     = 8
    variable    = phi_xy
  [../]
  [./couple_32]
    type = InternalCouple
    component_i = 2
    component_j = 1
    dof_num     = 9
    variable    = phi_zy
  [../]
  [./couple_31]
    type = InternalCouple
    component_i = 2
    component_j = 0
    dof_num     = 10
    variable    = phi_zx
  [../]
  [./couple_21]
    type = InternalCouple
    component_i = 1
    component_j = 0
    dof_num     = 11
    variable    = phi_yx
  [../]
[]

[BCs]
  [./left_x]
    type = DirichletBC
    variable = u_x
    boundary = left
    value = 0
  [../]
  [./back_z]
    type = DirichletBC
    variable = u_z
    boundary = back
    value = 0
  [../]
  [./bottom_y]
    type = DirichletBC
    variable = u_y
    boundary = bottom
    value = 0
  [../]
  [./top_y]
    type = DirichletBC
    variable = u_y
    boundary = top
    value = 0.1
  [../]
[]

[Materials]
  [./linear_elastic]
    type = MicromorphicMaterial
    material_fparameters = '8e9 11e9 2e9 1.538e9 -1e9 -1.39e9 -2.11e9 0. 0. 0. 0. 0. 0. 0.769e6 0. 0. 0. 0.'
    model_name = "LinearElasticity"

    #Coupled variables
    u1     = u_x
    u2     = u_y
    u3     = u_z
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

[Executioner]
  type = Steady
  solve_type = 'PJFNK'
  l_max_its  = 10
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]

[Outputs]
  exodus = true
[]
