Config.obj: Tools_log.obj M_kracken.obj
SFC_data_types.obj: Config.obj Tools_log.obj Heat_Equation/Heat_Eq_data_types.obj Darcy/Darcy_data_types.obj SWE/SWE_data_types.obj Tests/Tests_data_types.obj Generic/Generic_data_types.obj #NUMA/NUMA_data_types.obj
SFC_edge_traversal.obj: SFC_data_types.obj
SFC_node_traversal.obj: SFC_edge_traversal.obj
SFC_traversal.obj: SFC_node_traversal.obj Heat_Equation/Heat_Eq.obj Darcy/Darcy.obj Tests/Tests.obj SWE/SWE.obj Generic/Generic.obj #NUMA/NUMA.obj
SFC_main.obj: Config.obj M_kracken.obj SFC_traversal.obj
Tools_local_function_space_base.obj: SFC_data_types.obj
Tools_noise.obj: SFC_data_types.obj
Tools_stack_base.obj: SFC_data_types.obj


Conformity/Conformity.obj: SFC_edge_traversal.obj

geoclaw/c_bind_riemannsolvers.obj: geoclaw/riemannsolvers.obj geoclaw/riemannsolvers_sp.obj geoclaw/riemannsolvers_qp.obj

Tests/Tests.obj: Tests/Tests_data_types.obj Tests/Tests_basis_functions.obj Tests/Tests_initialize.obj Tests/Tests_node_dummy_traversal.obj Tests/Tests_flops_traversal.obj Tests/Tests_memory_traversal.obj Tests/Tests_consistency_traversal.obj Samoa/Samoa.obj
Tests/Tests_basis_functions.obj: SFC_data_types.obj
Tests/Tests_initialize.obj: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.obj
Tests/Tests_node_dummy_traversal.obj: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.obj
Tests/Tests_consistency_traversal.obj: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.obj
Tests/Tests_flops_traversal.obj: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.obj
Tests/Tests_memory_traversal.obj: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.obj

Darcy/Darcy.obj: Darcy/Darcy_data_types.obj Solver/LinearSolver.obj Darcy/Darcy_initialize.obj Darcy/Darcy_node_dummy.obj Darcy/Darcy_output.obj Darcy/Darcy_xml_output.obj Darcy/Darcy_grad_p.obj Darcy/Darcy_laplace_jacobi.obj Darcy/Darcy_laplace_cg.obj Darcy/Darcy_transport_eq.obj Darcy/Darcy_permeability.obj Darcy/Darcy_adapt.obj Samoa/Samoa.obj
Darcy/Darcy_local_function_spaces.obj: SFC_data_types.obj Tools_grid_variable.f90 Tools_local_function_space.f90
Darcy/Darcy_basis.obj: SFC_data_types.obj Darcy/Darcy_local_function_spaces.obj Samoa/Samoa.obj
Darcy/Darcy_adapt.obj: SFC_generic_adaptive_traversal.f90 Conformity/Conformity.obj Darcy/Darcy_local_function_spaces.obj Darcy/Darcy_basis.obj Samoa/Samoa.obj Darcy/Darcy_initialize.obj
Darcy/Darcy_node_dummy.obj: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.obj
Darcy/Darcy_grad_p.obj: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.obj Darcy/Darcy_local_function_spaces.obj Darcy/Darcy_basis.obj Samoa/Samoa.obj
Darcy/Darcy_initialize.obj: SFC_generic_traversal_ringbuffer.f90 Tools_noise.obj SFC_node_traversal.obj Darcy/Darcy_basis.obj Samoa/Samoa.obj
Darcy/Darcy_laplace_cg.obj: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.obj Darcy/Darcy_basis.obj Samoa/Samoa.obj
Darcy/Darcy_laplace_jacobi.obj: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.obj Darcy/Darcy_basis.obj Samoa/Samoa.obj
Darcy/Darcy_output.obj: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.obj Darcy/Darcy_basis.obj Samoa/Samoa.obj LIB_VTK_IO.obj
Darcy/Darcy_xml_output.obj: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.obj Darcy/Darcy_basis.obj Samoa/Samoa.obj LIB_VTK_IO.obj
Darcy/Darcy_transport_eq.obj: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.obj Darcy/Darcy_basis.obj Samoa/Samoa.obj
Darcy/Darcy_permeability.obj: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.obj Darcy/Darcy_basis.obj Samoa/Samoa.obj

Heat_Equation/Heat_Eq.obj: Heat_Equation/Heat_Eq_data_types.obj Heat_Equation/Heat_Eq_basis.obj Heat_Equation/Heat_Eq_initialize.obj Heat_Equation/Heat_Eq_output.obj Heat_Equation/Heat_Eq_xml_output.obj Heat_Equation/Heat_Eq_euler_timestep.obj Heat_Equation/Heat_Eq_midpoint_timestep.obj Heat_Equation/Heat_Eq_heun_timestep.obj Heat_Equation/Heat_Eq_adapt.obj Samoa/Samoa.obj
Heat_Equation/Heat_Eq_local_function_spaces.obj: SFC_data_types.obj
Heat_Equation/Heat_Eq_basis.obj: SFC_data_types.obj Heat_Equation/Heat_Eq_local_function_spaces.obj Samoa/Samoa.obj
Heat_Equation/Heat_Eq_adapt.obj: SFC_generic_adaptive_traversal.f90 Conformity/Conformity.obj Heat_Equation/Heat_Eq_basis.obj Samoa/Samoa.obj
Heat_Equation/Heat_Eq_euler_timestep.obj: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.obj Heat_Equation/Heat_Eq_basis.obj Samoa/Samoa.obj
Heat_Equation/Heat_Eq_heun_timestep.obj: SFC_generic_traversal_ringbuffer.f90 Heat_Equation/Heat_Eq_euler_timestep.obj
Heat_Equation/Heat_Eq_initialize.obj: SFC_generic_traversal_ringbuffer.f90 Tools_noise.obj SFC_node_traversal.obj Heat_Equation/Heat_Eq_basis.obj Samoa/Samoa.obj Heat_Equation/Heat_Eq_local_function_spaces.obj
Heat_Equation/Heat_Eq_midpoint_timestep.obj: SFC_generic_traversal_ringbuffer.f90 Heat_Equation/Heat_Eq_euler_timestep.obj
Heat_Equation/Heat_Eq_output.obj: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.obj Heat_Equation/Heat_Eq_basis.obj Heat_Equation/Heat_Eq_local_function_spaces.obj Samoa/Samoa.obj LIB_VTK_IO.obj
Heat_Equation/Heat_Eq_xml_output.obj: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.obj Heat_Equation/Heat_Eq_basis.obj Heat_Equation/Heat_Eq_local_function_spaces.obj Samoa/Samoa.obj LIB_VTK_IO.obj

SWE/SWE.obj: SWE/SWE_data_types.obj SWE/SWE_displace.obj SWE/SWE_basis.obj SWE/SWE_initialize.obj SWE/SWE_output.obj SWE/SWE_xml_output.obj SWE/SWE_euler_timestep.obj SWE/SWE_adapt.obj Samoa/Samoa.obj SWE/SWE_ascii_output.obj
SWE/SWE_local_function_spaces.obj: SFC_data_types.obj
SWE/SWE_basis.obj: SFC_data_types.obj SWE/SWE_local_function_spaces.obj Samoa/Samoa.obj
SWE/SWE_adapt.obj: SFC_generic_adaptive_traversal.f90 Conformity/Conformity.obj SWE/SWE_basis.obj Samoa/Samoa.obj SWE/SWE_euler_timestep.obj SWE/SWE_initialize.obj
SWE/SWE_euler_timestep.obj: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.obj SWE/SWE_basis.obj Samoa/Samoa.obj geoclaw/c_bind_riemannsolvers.obj
SWE/SWE_initialize.obj: SFC_generic_traversal_ringbuffer.f90 SWE/SWE_euler_timestep.obj Tools_noise.obj SFC_node_traversal.obj SWE/SWE_basis.obj Samoa/Samoa.obj SWE/SWE_local_function_spaces.obj
SWE/SWE_displace.obj: SWE/SWE_initialize.obj SFC_generic_traversal_ringbuffer.f90 SWE/SWE_euler_timestep.obj Tools_noise.obj SFC_node_traversal.obj SWE/SWE_basis.obj Samoa/Samoa.obj SWE/SWE_local_function_spaces.obj
SWE/SWE_output.obj: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.obj SWE/SWE_basis.obj SWE/SWE_local_function_spaces.obj SWE/SWE_euler_timestep.obj Samoa/Samoa.obj LIB_VTK_IO.obj
SWE/SWE_xml_output.obj: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.obj SWE/SWE_basis.obj SWE/SWE_local_function_spaces.obj SWE/SWE_euler_timestep.obj Samoa/Samoa.obj LIB_VTK_IO.obj
SWE/SWE_ascii_output.obj: SWE/ascii_output.obj SFC_edge_traversal.obj LIB_VTK_IO.obj SFC_edge_traversal.f90 Samoa/Samoa.obj SWE/SWE_euler_timestep.obj

#NUMA/NUMA.obj: NUMA/NUMA_data_types.obj NUMA/NUMA_basis.obj NUMA/NUMA_initialize.obj NUMA/NUMA_output.obj NUMA/NUMA_xml_output.obj NUMA/NUMA_euler_timestep.obj NUMA/NUMA_adapt.obj Samoa/Samoa.obj
#NUMA/NUMA_local_function_spaces.obj: SFC_data_types.obj
#NUMA/NUMA_basis.obj: SFC_data_types.obj NUMA/NUMA_local_function_spaces.obj Samoa/Samoa.obj
#NUMA/NUMA_adapt.obj: SFC_generic_adaptive_traversal.f90 Conformity/Conformity.obj NUMA/NUMA_basis.obj Samoa/Samoa.obj NUMA/NUMA_euler_timestep.obj NUMA/NUMA_initialize.obj
#NUMA/NUMA_euler_timestep.obj: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.obj NUMA/NUMA_basis.obj NUMA/NUMA_constants.obj Samoa/Samoa.obj geoclaw/c_bind_riemannsolvers.obj
#NUMA/NUMA_initialize.obj: SFC_generic_traversal_ringbuffer.f90 NUMA/NUMA_euler_timestep.obj Tools_noise.obj SFC_node_traversal.obj NUMA/NUMA_basis.obj Samoa/Samoa.obj NUMA/NUMA_local_function_spaces.obj
#NUMA/NUMA_output.obj: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.obj NUMA/NUMA_basis.obj NUMA/NUMA_local_function_spaces.obj NUMA/NUMA_euler_timestep.obj Samoa/Samoa.obj LIB_VTK_IO.obj
#NUMA/NUMA_xml_output.obj: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.obj NUMA/NUMA_basis.obj NUMA/NUMA_local_function_spaces.obj NUMA/NUMA_euler_timestep.obj Samoa/Samoa.obj LIB_VTK_IO.obj
#NUMA/NUMA_constants.obj: NUMA/NUMA_data_types.obj

Generic/Generic.obj: Generic/Generic_data_types.obj  SFC_node_traversal.obj Generic/Generic_initialize.obj Generic/Generic_template.obj Generic/Generic_adapt_template.obj
Generic/Generic_initialize.obj: SFC_generic_traversal_ringbuffer.f90 SFC_edge_traversal.obj Samoa/Samoa.obj
Generic/Generic_template.obj: SFC_generic_traversal_ringbuffer.f90 SFC_edge_traversal.obj Samoa/Samoa.obj
Generic/Generic_adapt_template.obj: SFC_generic_adaptive_traversal.f90 Conformity/Conformity.obj Samoa/Samoa.obj

Samoa/Samoa.obj: SFC_data_types.obj Samoa/Tools_quadrature_rule_base.obj
Samoa/Tools_quadrature_rule_base.obj: SFC_data_types.obj

Solver/LinearSolver.obj: SFC_edge_traversal.obj
