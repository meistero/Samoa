! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_HEAT_EQ)
	MODULE Heat_Eq_xml_output
		use LIB_VTK_IO

		use SFC_node_traversal

		use Samoa_heat_eq

		!> Output point data
		type t_output_point_data
			real (kind = GRID_SR), dimension(2)						:: coords		!< position
			real (kind = GRID_SR)									:: T			!< temperature
		END type t_output_point_data

		!> Output cell data
		type t_output_cell_data
			real (kind = GRID_SR)									:: heat_conductivity
			integer (kind = BYTE)										:: depth
			integer (kind = BYTE)										:: refinement
		END type t_output_cell_data


		logical, PARAMETER											:: l_second_order = (_HEAT_EQ_ORDER > 1)

		type(t_output_point_data), DIMENSION(:), ALLOCATABLE		:: point_data
		type(t_output_cell_data), DIMENSION(:), ALLOCATABLE			:: cell_data
		integer (kind = GRID_SI), DIMENSION(:), ALLOCATABLE			:: i_connectivity

		integer (kind = GRID_SI)									:: i_point_data_index
		integer (kind = GRID_SI)									:: i_cell_data_index

		type(heat_eq_gv_T)											:: gv_T
		type(heat_eq_gv_T_temp)										:: gv_T_temp

#		define	_GT_NAME							heat_eq_xml_output_traversal

#		if (_HEAT_EQ_EDGE_SIZE > 0)
#			define	_GT_EDGES
#		endif

#		define	_GT_NODES
#		define _GT_REFINEMENTS

#		define	_GT_PRE_TRAVERSAL_OP				pre_traversal_op
#		define	_GT_POST_TRAVERSAL_OP				post_traversal_op

#		define	_GT_ELEMENT_OP						element_op

#		define _GT_INNER_EDGE_FIRST_TOUCH_OP		edge_first_touch_op
#		define _GT_INNER_NODE_FIRST_TOUCH_OP		node_first_touch_op
#		define _GT_EDGE_FIRST_TOUCH_OP				edge_first_touch_op
#		define _GT_NODE_FIRST_TOUCH_OP				node_first_touch_op

#		include "SFC_generic_traversal_ringbuffer.f90"


		subroutine pre_traversal_op(traversal, section)
 			type(t_grid_section), intent(inout)									:: grid

			type(t_section_info)                                        :: grid_info
			integer (kind = GRID_SI)									:: i_error, i_cells, i_points

			_log_write(1, '(A, I0)') " Heat Eq: output step ", i_output_iteration

            grid_info = grid%get_info()
			i_cells = grid_info%i_cells

			if (l_second_order) then
				i_points = grid_info%i_crossed_edges + grid_info%i_color_edges + grid_info%i_nodes + sum(grid_info%i_boundary_edges + grid_info%i_boundary_nodes)
				allocate(i_connectivity(6 * i_cells), stat = i_error); assert_eq(i_error, 0)
			else
				i_points = grid_info%i_nodes + sum(grid_info%i_boundary_nodes)
				allocate(i_connectivity(3 * i_cells), stat = i_error); assert_eq(i_error, 0)
			end if

			allocate(cell_data(i_cells), stat = i_error); assert_eq(i_error, 0)
			allocate(point_data(i_points), stat = i_error); assert_eq(i_error, 0)

			i_cell_data_index = 1
			i_point_data_index = 0
		end subroutine

		subroutine post_traversal_op(traversal, section)
			type(t_grid_section), intent(inout)							:: grid

			integer (kind = GRID_SI), dimension(:), allocatable			:: i_offsets
			integer (1), dimension(:), allocatable						:: i_types
			real (kind = GRID_SR), dimension(:), allocatable			:: r_empty

			type(t_section_info)                                        :: grid_info
			integer (kind = GRID_SI)									:: i_error, i_cells, i_points
			character (len = 64)										:: s_file_name
			integer(4)													:: e_io, i

            grid_info = grid%get_info()
			i_cells = grid_info%i_cells

			if (l_second_order) then
				i_points = grid_info%i_crossed_edges + grid_info%i_color_edges + grid_info%i_nodes + sum(grid_info%i_boundary_edges + grid_info%i_boundary_nodes)
			else
				i_points = grid_info%i_nodes + sum(grid_info%i_boundary_nodes)
			end if

			allocate(i_offsets(i_cells), stat = i_error); assert_eq(i_error, 0)
			allocate(i_types(i_cells), stat = i_error); assert_eq(i_error, 0)
			allocate(r_empty(max(i_cells, i_points)), stat = i_error); assert_eq(i_error, 0)

			r_empty = 0.0_GRID_SR

			if (l_second_order) then
				i_types = 22_1

				forall (i = 1 : i_cells)
					i_offsets(i) = 6_GRID_SI * i
				end forall
			else
				i_types = 5_1

				forall (i = 1 : i_cells)
					i_offsets(i) = 3_GRID_SI * i
				end forall
			end if

			write (s_file_name, "(A, A, I0, A)") TRIM(grid%s_file_stamp), "_", i_output_iteration, ".vtu"

			e_io = VTK_INI_XML('binary', s_file_name, 'UnstructuredGrid')
				e_io = VTK_GEO_XML(i_points, i_cells, point_data%coords(1), point_data%coords(2), r_empty(1:i_points))

				e_io = VTK_CON_XML(i_cells, i_connectivity, i_offsets, i_types)

				e_io = VTK_DAT_XML('node', 'OPEN')
					e_io = VTK_VAR_XML(i_points, 'temperature', point_data%T)
				e_io = VTK_DAT_XML('node', 'CLOSE')

				e_io = VTK_DAT_XML('cell', 'OPEN')
					e_io = VTK_VAR_XML(i_cells, 'heat conductivity', cell_data%heat_conductivity)
					e_io = VTK_VAR_XML(i_cells, 'grid depth', cell_data%depth)
					e_io = VTK_VAR_XML(i_cells, 'refinement flag', cell_data%refinement)
				e_io = VTK_DAT_XML('cell', 'CLOSE')

				e_io = VTK_GEO_XML()
			e_io = VTK_END_XML()

			deallocate(i_offsets, stat = i_error); assert_eq(i_error, 0)
			deallocate(i_types, stat = i_error); assert_eq(i_error, 0)
			deallocate(i_connectivity, stat = i_error); assert_eq(i_error, 0)
			deallocate(r_empty, stat = i_error); assert_eq(i_error, 0)

			deallocate(cell_data, stat = i_error); assert_eq(i_error, 0)
			deallocate(point_data, stat = i_error); assert_eq(i_error, 0)

			i_output_iteration = i_output_iteration + 1
		end subroutine


		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, section, element)
			type(t_grid_section), intent(inout)							:: grid
			type(t_element_base), intent(inout)					:: element

			!local variables

			integer (kind = GRID_SI)							:: i

			real (kind = GRID_SR), dimension(_HEAT_EQ_SIZE)		:: T						!< Temperature DoFs
			real (kind = GRID_SR), dimension(_HEAT_EQ_SIZE)		:: r_point_data_indices		!< point data indices
			integer (kind = GRID_SR), dimension(_HEAT_EQ_SIZE)	:: point_data_indices		!< point data indices

			call gv_T%read(element, T)
			call gv_T_temp%read(element, r_point_data_indices)

			point_data_indices = r_point_data_indices
			T = samoa_basis_T_dofs_to_values(T)

			cell_data(i_cell_data_index)%heat_conductivity = element%cell%data_pers%heat_conductivity
			cell_data(i_cell_data_index)%depth = element%cell%geometry%i_depth
			cell_data(i_cell_data_index)%refinement = element%cell%geometry%refinement

			if (l_second_order) then
				forall (i = 1 : 6)
					point_data(point_data_indices(i))%coords = samoa_barycentric_to_world_point(element%transform_data, samoa_basis_T_get_dof_coords(i))
					point_data(point_data_indices(i))%T = T(i)
				end forall

				i_connectivity(6 * i_cell_data_index - 5 : 6 * i_cell_data_index) = point_data_indices((/ 1, 2, 3, 6, 4, 5 /)) - 1
			else
				forall (i = 1 : 3)
					point_data(point_data_indices(i))%coords = samoa_barycentric_to_world_point(element%transform_data, samoa_basis_T_get_dof_coords(i))
					point_data(point_data_indices(i))%T = T(i)
				end forall

				i_connectivity(3 * i_cell_data_index - 2 : 3 * i_cell_data_index) = point_data_indices(1:3) - 1
			end if

			i_cell_data_index = i_cell_data_index + 1
		end subroutine

		subroutine edge_first_touch_op(traversal, section, edge)
			type(t_grid_section), intent(inout)												:: grid
			type(t_edge_data), intent(inout)		:: edge

			if (l_second_order) then
				call pre_dof_op(edge%data_pers%T_temp)
			end if
		end subroutine

		subroutine node_first_touch_op(traversal, section, node)
			type(t_grid_section), intent(inout)												:: grid
			type(t_node_data), intent(inout)				:: node

			call pre_dof_op(node%data_pers%T_temp)
		end subroutine

		subroutine pre_dof_op(r)
			real(kind = GRID_SR), dimension(:), intent(out)				:: r

			integer														:: i

			forall (i = 1 : size(r))
				r(i) = i_point_data_index + i
			end forall

			i_point_data_index = i_point_data_index + size(r)
		end subroutine
	END MODULE
#endif
