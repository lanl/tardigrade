[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  micro_displacement_gradient = 'phi_xx phi_xy phi_xz phi_yx phi_yy phi_yz phi_zx phi_zy phi_zz'
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  xmax=10
  ymax=100
  zmax=10
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
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

[MicromorphicContinuum]
  [all]
    model_name = "LinearElasticityDruckerPragerPlasticity"
    material_fparameters = '2 1e3 1e4 
                            2 2e3 1e4 
                            2 1e3 1e4 
                            2 0.56 0.2
                            2 0. 0. 
                            2 0. 0. 
                            2 0. 0. 
                            2 0. 0. 
                            2 0. 0. 
                            2 29480 25480 
                            5 1000 400 -1500 -1400 -3000 
                            11 0 0 0 0 0 0 1e+06 0 0 0 0 
                            2 400 -3000 
                            0.5 0.5 0.5 
                            1e-09 1e-09'

    number_SDVS = 55

    save_in_disp_y = 'force_y'
  []
[]

[AuxVariables]
  [force_y]
  []
[]

[Postprocessors]
  [sum_bottom_force_y]
    type = NodalSum
    variable = force_y
    boundary = bottom
  []
[]

[BCs]
  [./left_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    preset = true
    value = 0
  [../]
  [./back_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    preset = true
    value = 0
  [../]
  [./bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    preset = true
    value = 0
  [../]
  [./top_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'top'
    preset = true
    function = top_bc
  [../]
[]

[Functions]
  [./top_bc]
    type  = ParsedFunction
    value = 5*t
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    petsc_options_iname = '     -pc_type
                                -pc_hypre_type
                                -ksp_type
                                -ksp_gmres_restart
                                -pc_hypre_boomeramg_relax_type_all
                                -pc_hypre_boomeramg_strong_threshold
                                -pc_hypre_boomeramg_agg_nl
                                -pc_hypre_boomeramg_agg_num_paths
                                -pc_hypre_boomeramg_max_levels
                                -pc_hypre_boomeramg_coarsen_type
                                -pc_hypre_boomeramg_interp_type
                                -pc_hypre_boomeramg_P_max
                                -pc_hypre_boomeramg_truncfactor' 

    petsc_options_value = '     hypre
                                boomeramg
                                gmres
                                301
                                symmetric-SOR/Jacobi
                                0.75
                                4 
                                2
                                25
                                Falgout
                                ext+i
                                0
                                0.1 '
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'

  nl_rel_tol = 1e-12
  nl_abs_tol = 5e-8
  l_tol = 1e-5
  l_max_its = 150
  nl_max_its = 12
  nl_div_tol = 1e3

  automatic_scaling=true
  compute_scaling_once =true
  verbose=false

  line_search = none

  dtmin = 1e-7
  dtmax= 5e-2
  
  start_time = 0.0
  end_time = 1.0 

  num_steps = 2000

  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 8
    iteration_window = 3
    linear_iteration_ratio = 1000
    growth_factor=1.2
    cutback_factor=0.5
    dt = 5e-2
  []
  [Quadrature]
    order=SECOND
  []
  [Predictor]
    type = SimplePredictor
    scale = 1.0
    skip_after_failed_timestep = true
  []
[] 

[Outputs]
  interval = 1
  print_linear_residuals = false
  csv = true
  exodus = true 
  [pgraph]
    type = PerfGraphOutput
    execute_on = 'timestep_end final'  # Default is "final"
    level = 2             # Default is 1
  []
[]
