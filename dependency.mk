SFC_data_types.o: Tools_log.o Heat_Equation/Heat_Eq_data_types.o Darcy/Darcy_data_types.o SWE/SWE_data_types.o Tests/Tests_data_types.o NUMA/NUMA_data_types.o PyOP2/PyOP2_data_types.o
SFC_edge_traversal.o: SFC_data_types.o
SFC_node_traversal.o: SFC_edge_traversal.o
SFC_traversal.o: SFC_node_traversal.o Heat_Equation/Heat_Eq.o Darcy/Darcy.o Tests/Tests.o SWE/SWE.o NUMA/NUMA.o PyOP2/PyOP2.o
SFC_main.o: M_kracken.o SFC_traversal.o
Tools_local_function_space_base.o: SFC_data_types.o
Tools_noise.o: SFC_data_types.o
Tools_stack_base.o: SFC_data_types.o

Conformity/Conformity.o: SFC_edge_traversal.o

geoclaw/c_bind_riemannsolvers.o: geoclaw/riemannsolvers.o

Tests/Tests.o: Tests/Tests_data_types.o Tests/Tests_basis_functions.o Tests/Tests_initialize.o Tests/Tests_node_dummy_traversal.o Tests/Tests_flops_traversal.o Tests/Tests_memory_traversal.o Tests/Tests_consistency_traversal.o Samoa/Samoa.o
Tests/Tests_basis_functions.o: SFC_data_types.o
Tests/Tests_initialize.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o
Tests/Tests_node_dummy_traversal.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o
Tests/Tests_consistency_traversal.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o
Tests/Tests_flops_traversal.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o
Tests/Tests_memory_traversal.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o

Darcy/Darcy.o: Darcy/Darcy_data_types.o Darcy/Darcy_initialize.o Darcy/Darcy_node_dummy.o Darcy/Darcy_output.o Darcy/Darcy_xml_output.o Darcy/Darcy_grad_p.o Darcy/Darcy_laplace_jacobi.o Darcy/Darcy_laplace_cg.o Darcy/Darcy_transport_eq.o Darcy/Darcy_permeability.o Darcy/Darcy_adapt.o Samoa/Samoa.o
Darcy/Darcy_local_function_spaces.o: SFC_data_types.o
Darcy/Darcy_basis.o: SFC_data_types.o Darcy/Darcy_local_function_spaces.o Samoa/Samoa.o
Darcy/Darcy_adapt.o: SFC_generic_adaptive_traversal.f90 Conformity/Conformity.o Darcy/Darcy_local_function_spaces.o Darcy/Darcy_basis.o Samoa/Samoa.o Darcy/Darcy_initialize.o
Darcy/Darcy_node_dummy.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o
Darcy/Darcy_grad_p.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Darcy/Darcy_local_function_spaces.o Darcy/Darcy_basis.o Samoa/Samoa.o
Darcy/Darcy_initialize.o: SFC_generic_traversal_ringbuffer.f90 Tools_noise.o SFC_node_traversal.o Darcy/Darcy_basis.o Samoa/Samoa.o
Darcy/Darcy_laplace_cg.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Darcy/Darcy_basis.o Samoa/Samoa.o
Darcy/Darcy_laplace_jacobi.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Darcy/Darcy_basis.o Samoa/Samoa.o
Darcy/Darcy_output.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Darcy/Darcy_basis.o Samoa/Samoa.o LIB_VTK_IO.o
Darcy/Darcy_xml_output.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Darcy/Darcy_basis.o Samoa/Samoa.o LIB_VTK_IO.o
Darcy/Darcy_transport_eq.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Darcy/Darcy_basis.o Samoa/Samoa.o
Darcy/Darcy_permeability.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Darcy/Darcy_basis.o Samoa/Samoa.o

Heat_Equation/Heat_Eq.o: Heat_Equation/Heat_Eq_data_types.o Heat_Equation/Heat_Eq_basis.o Heat_Equation/Heat_Eq_initialize.o Heat_Equation/Heat_Eq_output.o Heat_Equation/Heat_Eq_xml_output.o Heat_Equation/Heat_Eq_euler_timestep.o Heat_Equation/Heat_Eq_midpoint_timestep.o Heat_Equation/Heat_Eq_heun_timestep.o Heat_Equation/Heat_Eq_adapt.o Samoa/Samoa.o
Heat_Equation/Heat_Eq_local_function_spaces.o: SFC_data_types.o
Heat_Equation/Heat_Eq_basis.o: SFC_data_types.o Heat_Equation/Heat_Eq_local_function_spaces.o Samoa/Samoa.o
Heat_Equation/Heat_Eq_adapt.o: SFC_generic_adaptive_traversal.f90 Conformity/Conformity.o Heat_Equation/Heat_Eq_basis.o Samoa/Samoa.o
Heat_Equation/Heat_Eq_euler_timestep.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Heat_Equation/Heat_Eq_basis.o Samoa/Samoa.o
Heat_Equation/Heat_Eq_heun_timestep.o: SFC_generic_traversal_ringbuffer.f90 Heat_Equation/Heat_Eq_euler_timestep.o
Heat_Equation/Heat_Eq_initialize.o: SFC_generic_traversal_ringbuffer.f90 Tools_noise.o SFC_node_traversal.o Heat_Equation/Heat_Eq_basis.o Samoa/Samoa.o Heat_Equation/Heat_Eq_local_function_spaces.o
Heat_Equation/Heat_Eq_midpoint_timestep.o: SFC_generic_traversal_ringbuffer.f90 Heat_Equation/Heat_Eq_euler_timestep.o
Heat_Equation/Heat_Eq_output.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Heat_Equation/Heat_Eq_basis.o Heat_Equation/Heat_Eq_local_function_spaces.o Samoa/Samoa.o LIB_VTK_IO.o
Heat_Equation/Heat_Eq_xml_output.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Heat_Equation/Heat_Eq_basis.o Heat_Equation/Heat_Eq_local_function_spaces.o Samoa/Samoa.o LIB_VTK_IO.o

SWE/SWE.o: SWE/SWE_data_types.o SWE/SWE_basis.o SWE/SWE_initialize.o SWE/SWE_output.o SWE/SWE_xml_output.o SWE/SWE_euler_timestep.o SWE/SWE_adapt.o Samoa/Samoa.o
SWE/SWE_local_function_spaces.o: SFC_data_types.o
SWE/SWE_basis.o: SFC_data_types.o SWE/SWE_local_function_spaces.o Samoa/Samoa.o
SWE/SWE_adapt.o: SFC_generic_adaptive_traversal.f90 Conformity/Conformity.o SWE/SWE_basis.o Samoa/Samoa.o SWE/SWE_euler_timestep.o SWE/SWE_initialize.o
SWE/SWE_euler_timestep.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o SWE/SWE_basis.o Samoa/Samoa.o geoclaw/c_bind_riemannsolvers.o
SWE/SWE_initialize.o: SFC_generic_traversal_ringbuffer.f90 SWE/SWE_euler_timestep.o Tools_noise.o SFC_node_traversal.o SWE/SWE_basis.o Samoa/Samoa.o SWE/SWE_local_function_spaces.o
SWE/SWE_output.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o SWE/SWE_basis.o SWE/SWE_local_function_spaces.o SWE/SWE_euler_timestep.o Samoa/Samoa.o LIB_VTK_IO.o
SWE/SWE_xml_output.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o SWE/SWE_basis.o SWE/SWE_local_function_spaces.o SWE/SWE_euler_timestep.o Samoa/Samoa.o LIB_VTK_IO.o

NUMA/NUMA.o: NUMA/NUMA_data_types.o NUMA/NUMA_basis.o NUMA/NUMA_initialize.o NUMA/NUMA_output.o NUMA/NUMA_xml_output.o NUMA/NUMA_euler_timestep.o NUMA/NUMA_adapt.o Samoa/Samoa.o
NUMA/NUMA_local_function_spaces.o: SFC_data_types.o
NUMA/NUMA_basis.o: SFC_data_types.o NUMA/NUMA_local_function_spaces.o Samoa/Samoa.o
NUMA/NUMA_adapt.o: SFC_generic_adaptive_traversal.f90 Conformity/Conformity.o NUMA/NUMA_basis.o Samoa/Samoa.o NUMA/NUMA_euler_timestep.o NUMA/NUMA_initialize.o
NUMA/NUMA_euler_timestep.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o NUMA/NUMA_basis.o NUMA/NUMA_constants.o Samoa/Samoa.o geoclaw/c_bind_riemannsolvers.o
NUMA/NUMA_initialize.o: SFC_generic_traversal_ringbuffer.f90 NUMA/NUMA_euler_timestep.o Tools_noise.o SFC_node_traversal.o NUMA/NUMA_basis.o Samoa/Samoa.o NUMA/NUMA_local_function_spaces.o
NUMA/NUMA_output.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o NUMA/NUMA_basis.o NUMA/NUMA_local_function_spaces.o NUMA/NUMA_euler_timestep.o Samoa/Samoa.o LIB_VTK_IO.o
NUMA/NUMA_xml_output.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o NUMA/NUMA_basis.o NUMA/NUMA_local_function_spaces.o NUMA/NUMA_euler_timestep.o Samoa/Samoa.o LIB_VTK_IO.o
NUMA/NUMA_constants.o: NUMA/NUMA_data_types.o

PyOP2/PyOP2.o: PyOP2/PyOP2_data_types.o PyOP2/PyOP2_initialize.o PyOP2/PyOP2_template.o PyOP2/PyOP2_adapt_template.o
PyOP2/PyOP2_initialize.o: SFC_generic_traversal_ringbuffer.f90 SFC_edge_traversal.o SFC_node_traversal.o Samoa/Samoa.o
PyOP2/PyOP2_template.o: SFC_generic_traversal_ringbuffer.f90 SFC_edge_traversal.o SFC_node_traversal.o Samoa/Samoa.o PyOP2/PyOP2_initialize.o
PyOP2/PyOP2_adapt_template.o: SFC_generic_adaptive_traversal.f90 Conformity/Conformity.o Samoa/Samoa.o PyOP2/PyOP2_initialize.o

Samoa/Samoa.o: SFC_data_types.o Samoa/Tools_quadrature_rule_base.o
Samoa/Tools_quadrature_rule_base.o: SFC_data_types.o
