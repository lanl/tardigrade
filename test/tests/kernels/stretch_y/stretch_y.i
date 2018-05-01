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
    initial_condition = 1.22
  [../]
  [./u_y]
    initial_condition = 2.1
  [../]
  [./u_z]
    initial_condition = 4.1
  [../]
  [./phi_xx]
    initial_condition = -2.3
  [../]
  [./phi_yy]
    initial_condition = 0.124
  [../]
  [./phi_zz]
    initial_condition = 7.2
  [../]
  [./phi_yz]
    initial_condition = -8.2
  [../]
  [./phi_xz]
    initial_condition = 0.28
  [../]
  [./phi_xy]
    initial_condition = 7.21
  [../]
  [./phi_zy]
    initial_condition = 2.1
  [../]
  [./phi_zx]
    initial_condition = -9.2
  [../]
  [./phi_yx]
    initial_condition = 3.1
  [../]
[]

[Kernels]
  #Define the internal force balance equations
  [./force_1]
    type = InternalForce
    component = 0
    dof_num   = 0
    variable  = u_x

    #Coupled variables
    u1     = u_x
    u2     = u_y
    u3     = u_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_zy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./force_2]
    type = InternalForce
    component = 1
    dof_num   = 1
    variable  = u_y

    #Coupled variables
    u1     = u_x
    u2     = u_y
    u3     = u_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_zy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./force_3]
    type = InternalForce
    component = 2
    dof_num   = 2
    variable  = u_z

    #Coupled variables
    u1     = u_x
    u2     = u_y
    u3     = u_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_zy
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
    u1     = u_x
    u2     = u_y
    u3     = u_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_zy
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

    #Coupled variables
    u1     = u_x
    u2     = u_y
    u3     = u_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_zy
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

    #Coupled variables
    u1     = u_x
    u2     = u_y
    u3     = u_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_zy
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

    #Coupled variables
    u1     = u_x
    u2     = u_y
    u3     = u_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_zy
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

    #Coupled variables
    u1     = u_x
    u2     = u_y
    u3     = u_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_zy
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

    #Coupled variables
    u1     = u_x
    u2     = u_y
    u3     = u_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_zy
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

    #Coupled variables
    u1     = u_x
    u2     = u_y
    u3     = u_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_zy
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

    #Coupled variables
    u1     = u_x
    u2     = u_y
    u3     = u_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_zy
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

    #Coupled variables
    u1     = u_x
    u2     = u_y
    u3     = u_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_zy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
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

#[Preconditioning]
#  [./SMP]
#    type = SMP
#    full = true
#  [../]
#[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK'
  l_max_its  = 5
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]

[Outputs]
  exodus = true
[]
