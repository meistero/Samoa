# DO NOT EDIT --- auto-generated file
SFC_checkpointing.obj: SFC_data_types.obj SFC_node_traversal.obj
SFC_containers.obj: SFC_data_types.obj
SFC_data_types.obj: Tools_log.obj Heat_Equation/Heat_Eq_data_types.obj Darcy/Darcy_data_types.obj Tests/Tests_data_types.obj
SFC_fine_grid.obj: SFC_data_types.obj Tools_stack_base.obj
SFC_data_types.obj: SFC_data_types.obj
SFC_edge_traversal.obj: SFC_data_types.obj SFC_fine_grid.obj SFC_containers.obj
SFC_node_traversal.obj: SFC_edge_traversal.obj
SFC_traversal.obj: SFC_node_traversal.obj Heat_Equation/Heat_Eq.obj Darcy/Darcy.obj Tests/Tests.obj
SFC_main.obj: M_kracken.obj SFC_traversal.obj
Tools_local_function_space_base.obj: SFC_data_types.obj
Tools_noise.obj: SFC_data_types.obj
Tools_stack_base.obj: SFC_data_types.obj

Tests/Tests.obj: Tests/Tests_data_types.obj Tests/Tests_basis_functions.obj Tests/Tests_initialize.obj Tests/Tests_node_dummy_traversal.obj Tests/Tests_flops_traversal.obj Tests/Tests_memory_traversal.obj Tests/Tests_consistency_traversal.obj Samoa/Samoa.obj
Tests/Tests_initialize.obj: SFC_node_traversal.obj
Tests/Tests_node_dummy_traversal.obj: SFC_node_traversal.obj
Tests/Tests_consistency_traversal.obj: SFC_node_traversal.obj
Tests/Tests_flops_traversal.obj: SFC_node_traversal.obj
Tests/Tests_memory_traversal.obj: SFC_node_traversal.obj
Tests/Tests_basis_functions.obj: SFC_data_types.obj

Darcy/Darcy.obj: Darcy/Darcy_data_types.obj Darcy/Darcy_initialize.obj Darcy/Darcy_node_dummy.obj Darcy/Darcy_output.obj Darcy/Darcy_xml_output.obj Darcy/Darcy_grad_p.obj Darcy/Darcy_laplace_jacobi.obj Darcy/Darcy_laplace_cg.obj Darcy/Darcy_transport_eq.obj Darcy/Darcy_permeability.obj Darcy/Darcy_adapt.obj Samoa/Samoa.obj
Darcy/Darcy_adapt.obj: SFC_node_traversal.obj Darcy/Darcy_local_function_spaces.obj Samoa/Samoa_quadrature_rules.obj Darcy/Darcy_basis.obj Samoa/Samoa.obj
Darcy/Darcy_node_dummy.obj: SFC_node_traversal.obj
Darcy/Darcy_grad_p.obj: SFC_node_traversal.obj Darcy/Darcy_local_function_spaces.obj Samoa/Samoa_quadrature_rules.obj Darcy/Darcy_basis.obj Samoa/Samoa.obj
Darcy/Darcy_initialize.obj: Tools_noise.obj SFC_node_traversal.obj Samoa/Samoa_quadrature_rules.obj Darcy/Darcy_basis.obj Samoa/Samoa.obj
Darcy/Darcy_laplace_cg.obj: SFC_node_traversal.obj Samoa/Samoa_quadrature_rules.obj Darcy/Darcy_basis.obj Samoa/Samoa.obj
Darcy/Darcy_laplace_jacobi.obj: SFC_node_traversal.obj Samoa/Samoa_quadrature_rules.obj Darcy/Darcy_basis.obj Samoa/Samoa.obj
Darcy/Darcy_local_function_spaces.obj: Tools_local_function_space_base.obj SFC_data_types.obj
Darcy/Darcy_output.obj: SFC_node_traversal.obj Darcy/Darcy_basis.obj Samoa/Samoa.obj
Darcy/Darcy_xml_output.obj: SFC_node_traversal.obj Darcy/Darcy_basis.obj Samoa/Samoa.obj LIB_VTK_IO.obj
Darcy/Darcy_transport_eq.obj: SFC_node_traversal.obj Darcy/Darcy_basis.obj Samoa/Samoa.obj
Darcy/Darcy_permeability.obj: SFC_node_traversal.obj Darcy/Darcy_basis.obj Samoa/Samoa.obj
Darcy/Darcy_basis.obj: SFC_data_types.obj Darcy/Darcy_local_function_spaces.obj Samoa/Samoa.obj

Heat_Equation/Heat_Eq.obj: Heat_Equation/Heat_Eq_data_types.obj Heat_Equation/Heat_Eq_basis.obj Heat_Equation/Heat_Eq_initialize.obj Heat_Equation/Heat_Eq_output.obj Heat_Equation/Heat_Eq_xml_output.obj Heat_Equation/Heat_Eq_euler_timestep.obj Heat_Equation/Heat_Eq_midpoint_timestep.obj Heat_Equation/Heat_Eq_heun_timestep.obj Heat_Equation/Heat_Eq_adapt.obj Samoa/Samoa.obj
Heat_Equation/Heat_Eq_adapt.obj: SFC_node_traversal.obj Heat_Equation/Heat_Eq_basis.obj Samoa/Samoa.obj
Heat_Equation/Heat_Eq_euler_timestep.obj: SFC_node_traversal.obj Heat_Equation/Heat_Eq_basis.obj Samoa/Samoa.obj
Heat_Equation/Heat_Eq_heun_timestep.obj: Heat_Equation/Heat_Eq_euler_timestep.obj
Heat_Equation/Heat_Eq_initialize.obj: Tools_noise.obj SFC_node_traversal.obj Heat_Equation/Heat_Eq_basis.obj Samoa/Samoa.obj Heat_Equation/Heat_Eq_local_function_spaces.obj
Heat_Equation/Heat_Eq_local_function_spaces.obj: Tools_local_function_space_base.obj SFC_data_types.obj
Heat_Equation/Heat_Eq_midpoint_timestep.obj: Heat_Equation/Heat_Eq_euler_timestep.obj
Heat_Equation/Heat_Eq_output.obj: SFC_node_traversal.obj Heat_Equation/Heat_Eq_basis.obj Heat_Equation/Heat_Eq_local_function_spaces.obj Samoa/Samoa.obj
Heat_Equation/Heat_Eq_xml_output.obj: SFC_node_traversal.obj Heat_Equation/Heat_Eq_basis.obj Heat_Equation/Heat_Eq_local_function_spaces.obj Samoa/Samoa.obj LIB_VTK_IO.obj
Heat_Equation/Heat_Eq_basis.obj: SFC_data_types.obj Heat_Equation/Heat_Eq_local_function_spaces.obj Samoa/Samoa.obj

Shallow_Water/Shallow_Water.obj: Shallow_Water/Shallow_Water_data_types.obj Shallow_Water/Shallow_Water_basis.obj Shallow_Water/Shallow_Water_initialize.obj Shallow_Water/Shallow_Water_xml_output.obj Shallow_Water/Shallow_Water_euler_timestep.obj Shallow_Water/Shallow_Water_adapt.obj Samoa/Samoa.obj
Shallow_Water/Shallow_Water_adapt.obj: SFC_node_traversal.obj Shallow_Water/Shallow_Water_basis.obj Samoa/Samoa.obj
Shallow_Water/Shallow_Water_euler_timestep.obj: SFC_node_traversal.obj Shallow_Water/Shallow_Water_basis.obj Samoa/Samoa.obj
Shallow_Water/Shallow_Water_initialize.obj: Tools_noise.obj SFC_node_traversal.obj Shallow_Water/Shallow_Water_basis.obj Samoa/Samoa.obj Shallow_Water/Shallow_Water_local_function_spaces.obj
Shallow_Water/Shallow_Water_local_function_spaces.obj: Tools_local_function_space_base.obj SFC_data_types.obj
Shallow_Water/Shallow_Water_xml_output.obj: SFC_node_traversal.obj Shallow_Water/Shallow_Water_basis.obj Shallow_Water/Shallow_Water_local_function_spaces.obj Samoa/Samoa.obj LIB_VTK_IO.obj
Shallow_Water/Shallow_Water_basis.obj: SFC_data_types.obj Shallow_Water/Shallow_Water_local_function_spaces.obj Samoa/Samoa.obj

Samoa/Samoa.obj: SFC_data_types.obj Tools_local_function_space_base.obj Samoa/Tools_quadrature_rule_base.obj
Samoa/Samoa_quadrature_rules.obj: SFC_data_types.obj Samoa/Samoa.obj Darcy/Darcy_local_function_spaces.obj Heat_Equation/Heat_Eq_local_function_spaces.obj Samoa/Samoa.obj
Samoa/Tools_quadrature_rule_base.obj: SFC_data_types.obj
Samoa/Tools_quadrature_rule.obj: SFC_data_types.obj Samoa/Tools_quadrature_rule_base.obj
