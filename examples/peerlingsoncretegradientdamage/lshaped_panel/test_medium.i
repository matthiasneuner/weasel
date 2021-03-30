[Mesh]
  type = FileMesh
  file = ../../meshes/lshape_medium.e
[]

[GlobalParams]
  volumetric_locking_correction = false
  displacements = 'disp_x disp_y disp_z'
  nonlocal_damage = 'nonlocal_equivalent_strain'
  order = FIRST
  family = LAGRANGE
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
  [nonlocal_equivalent_strain]
  []
[]

[Kernels]
  [div_sig_x]
    type = WeaselStressDivergence
    variable = disp_x
    component = 0
    save_in = force_x
  []
  [div_sig_y]
    type = WeaselStressDivergence
    variable = disp_y
    component = 1
    save_in = force_y
  []
  [div_sig_z]
    type = WeaselStressDivergence
    variable = disp_z
    component = 2
    save_in = force_z
  []
  [nonlocal_damage_kernel]
    type = WeaselGradientDamage
    variable = nonlocal_equivalent_strain
  []
[]

[AuxVariables]
  [von_mises]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_first_invariant]
    order = CONSTANT
    family = MONOMIAL
  []
  [equivalent_strain]
    order = CONSTANT
    family = MONOMIAL
  []
  [damage]
    order = CONSTANT
    family = MONOMIAL
  []
  [force_x]
  []
  [force_y]
  []
  [force_z]
  []
[]


[Postprocessors]
  [bot_react_y]
    type = NodalSum
    variable = force_y
    boundary = boundary_bottom
  []
[]

[AuxKernels]
  [von_mises_kernel]
    type = RankTwoScalarAux
    variable = von_mises
    rank_two_tensor = stress
    execute_on = timestep_end
    scalar_type = VonMisesStress
  []
  [first_invariant_kernel]
    type = RankTwoScalarAux
    variable = stress_first_invariant
    rank_two_tensor = stress
    execute_on = timestep_end
    scalar_type = FirstInvariant
  []
  [equivalent_strain_kernel]
    type = MaterialRealAux
    property = kappa
    variable = equivalent_strain
    execute_on = timestep_end
  []
  [damage_kernel]
    type = MaterialRealAux
    property = damage
    variable = damage
    execute_on = timestep_end
  []
[]


[Materials]
[strain]
    type = ComputeSmallStrain
[]
[von_mises]
    type = WeaselMaterialPeerlingsConcreteGradientDamage
    E = 30000
    nu = 0.2
    kappa_0 = 4e-4
    alpha = 0.995
    beta = 500
    k = 10
    nonlocal_radius = 5
[]
[]

[BCs]
  [bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = boundary_bottom
    value = 0
  []
  [bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = boundary_bottom
    value = 0
  []
  [bottom_z]
    type = DirichletBC
    variable = disp_z
    boundary = boundary_bottom 
    value = 0
  []
  [load]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = boundary_load
    function = '1 * t'
    preset = false
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'

  petsc_options_iname = '-pc_type   -pc_hypre_type    -ksp_type     -ksp_gmres_restart  -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_max_levels -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_interp_type -pc_hypre_boomeramg_P_max -pc_hypre_boomeramg_truncfactor'
  petsc_options_value = 'hypre      boomeramg         gmres         301                  0.6                                  4                          5                                 25                             Falgout                          ext+i                           1                         0.3'

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  l_tol = 1e-5
  l_max_its = 300
  nl_max_its = 20

  line_search = 'none'

  automatic_scaling=true
  compute_scaling_once =true
  verbose=false

  dtmin = 1e-5
  dtmax= 1e-2
  
  start_time = 0.0
  end_time = 1.0 

  num_steps = 1000
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 15
    iteration_window = 3
    linear_iteration_ratio = 100
    growth_factor=1.5
    cutback_factor=0.5
    dt = 1e-2
  []
  [Predictor]
    type = SimplePredictor
    scale = 1.0
    skip_after_failed_timestep = true
  []
[] 

[Outputs]
  interval = 1
  #execute_on = 'initial timestep_end'
  print_linear_residuals = false
  csv = true
  exodus = true
  [pgraph]
    type = PerfGraphOutput
    execute_on = 'final'  # Default is "final"
    level = 1             # Default is 1
  []
[]
