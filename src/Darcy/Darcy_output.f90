! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_vtk_output
		use LIB_VTK_IO

		use SFC_edge_traversal

		use Samoa_darcy

		!> Output point data
		type t_output_point_data
			real (kind = GRID_SR), dimension(2)						:: coords
			real (kind = GRID_SR)									:: p
			real (kind = GRID_SR)									:: S
		END type t_output_point_data

		!> Output cell data
		type t_output_cell_data
			real (kind = GRID_SR)									:: permeability
			real (kind = GRID_SR), dimension(2)						:: u
			integer (kind = GRID_SI)								:: depth
			integer (kind = GRID_SI)								:: refinement
		END type t_output_cell_data

		logical, PARAMETER											:: l_second_order = .false.

        type num_traversal_data
            type(t_output_point_data), allocatable	                :: point_data(:)
            type(t_output_cell_data), allocatable		            :: cell_data(:)
            integer (kind = GRID_SI), allocatable		            :: i_connectivity(:)
            character(len=256)							            :: s_file_stamp

            integer (kind = GRID_SI)								:: i_output_iteration
            integer (kind = GRID_SI)								:: i_point_data_index
            integer (kind = GRID_SI)								:: i_cell_data_index
        end type

        interface node_first_touch_op
            module procedure node_first_touch_op_scalar
            module procedure node_first_touch_op_array
        end interface

        type(darcy_gv_p)											:: gv_p
        !type(darcy_gv_u)											:: gv_u
        type(darcy_gv_saturation)									:: gv_saturation
        type(darcy_gv_r)											:: gv_r

#		define	_GT_NAME							t_darcy_vtk_output_traversal

#		if (_DARCY_P_EDGE_SIZE > 0)
#			define _GT_EDGES
#		endif

#		define	_GT_NODES
#		define	_GT_REFINEMENTS

#		define	_GT_PRE_TRAVERSAL_OP				pre_traversal_op
#		define	_GT_POST_TRAVERSAL_OP				post_traversal_op

#		define	_GT_ELEMENT_OP						element_op

#		define _GT_NODE_FIRST_TOUCH_OP				node_first_touch_op

#		include "SFC_generic_traversal_ringbuffer.f90"


		subroutine pre_traversal_op(traversal, section)
 			type(t_darcy_vtk_output_traversal), intent(inout)				:: traversal
  			type(t_grid_section), intent(inout)							:: section

			integer (kind = GRID_SI)									:: i_error, i_cells, i_points
			type(t_section_info)                                           :: grid_info

            if (rank_MPI == 0) then
                _log_write(1, '(A, I0, A, I0)') " Darcy: output step: ", traversal%i_output_iteration
            end if

            grid_info = section%get_info()
			i_cells = grid_info%i_cells

			if (l_second_order) then
				i_points = grid_info%i_crossed_edges + grid_info%i_color_edges + grid_info%i_nodes + sum(grid_info%i_boundary_edges + grid_info%i_boundary_nodes)
				allocate(traversal%i_connectivity(7 * i_cells), stat = i_error); assert_eq(i_error, 0)
			else
				i_points = grid_info%i_nodes + sum(grid_info%i_boundary_nodes)
				allocate(traversal%i_connectivity(4 * i_cells), stat = i_error); assert_eq(i_error, 0)
			end if

			allocate(traversal%cell_data(i_cells), stat = i_error); assert_eq(i_error, 0)
			allocate(traversal%point_data(i_points), stat = i_error); assert_eq(i_error, 0)

			traversal%i_cell_data_index = 1
			traversal%i_point_data_index = 1
		end subroutine

		subroutine post_traversal_op(traversal, section)
 			type(t_darcy_vtk_output_traversal), intent(inout)		    :: traversal
 			type(t_grid_section), intent(inout)							:: section

			integer (kind = GRID_SI), dimension(:), allocatable			:: i_offsets
			integer (kind = GRID_SI), dimension(:), allocatable			:: i_types
			real (kind = GRID_SR), dimension(:), allocatable			:: r_empty
            type(t_vtk_writer)                                          :: vtk

			type(t_section_info)                                        :: grid_info
			integer (kind = GRID_SI)									:: i_error, i_cells, i_points
			character (len = 256)										:: s_file_name
			integer(4)													:: e_io, i

            grid_info = section%get_info()
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

			write (s_file_name, "(A, A, I0, A)") TRIM(traversal%s_file_stamp), "_", traversal%i_output_iteration, ".vtk"

#           if defined(_QUAD_PRECISION)
#               warning VTK output does not work for quad precision
#           else
                e_io = vtk%VTK_INI('ASCII', s_file_name, 'Darcy', 'UNSTRUCTURED_GRID')
                    e_io = vtk%VTK_GEO(i_points, traversal%point_data%coords(1), traversal%point_data%coords(2), r_empty(1:i_points))

                    e_io = vtk%VTK_CON(i_cells, traversal%i_connectivity, i_types)

                    e_io = vtk%VTK_DAT(i_points, 'node')
                    e_io = vtk%VTK_VAR(i_points, 'pressure', traversal%point_data%p)
                    e_io = vtk%VTK_VAR(i_points, 'saturation', traversal%point_data%S)

                    e_io = vtk%VTK_DAT(i_cells, 'cell')
                    e_io = vtk%VTK_VAR(i_cells, 'permeability', traversal%cell_data%permeability)
                    e_io = vtk%VTK_VAR('vect', i_cells, 'velocity', traversal%cell_data%u(1), traversal%cell_data%u(2), r_empty(1:i_cells))
                    e_io = vtk%VTK_VAR(i_cells, 'grid_depth', traversal%cell_data%depth)
                    e_io = vtk%VTK_VAR(i_cells, 'refinement_flag', traversal%cell_data%refinement)
                e_io = vtk%VTK_END()
#           endif

			deallocate(i_offsets, stat = i_error); assert_eq(i_error, 0)
			deallocate(i_types, stat = i_error); assert_eq(i_error, 0)
			deallocate(r_empty, stat = i_error); assert_eq(i_error, 0)

			deallocate(traversal%i_connectivity, stat = i_error); assert_eq(i_error, 0)
			deallocate(traversal%cell_data, stat = i_error); assert_eq(i_error, 0)
			deallocate(traversal%point_data, stat = i_error); assert_eq(i_error, 0)

			traversal%i_output_iteration = traversal%i_output_iteration + 1
		end subroutine


		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, section, element)
 			type(t_darcy_vtk_output_traversal), intent(inout)				:: traversal
 			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)					:: element

			!local variables

			integer (kind = GRID_SI)    :: i

			!local variables

			real (kind = GRID_SR)       :: p(3)
			real (kind = GRID_SR)       :: u(2)
			real (kind = GRID_SR)       :: saturation(3)

			real (kind = GRID_SR)       :: r_point_data_indices(3)		!< point data indices
			integer (kind = GRID_SI)    :: point_data_indices(3)		!< point data indices

			call gv_p%read(element, p)
			call gv_r%read(element, r_point_data_indices)
			!call gv_u%read_from_element(element, u)
			call gv_saturation%read(element, saturation)

			point_data_indices = r_point_data_indices
			p = p / 6.89e3_GRID_SR
			saturation = samoa_basis_flow_dofs_to_values(saturation)

#           if (_DARCY_LAYERS > 0)
                !TODO
#           else
                traversal%cell_data(traversal%i_cell_data_index)%permeability = element%cell%data_pers%base_permeability
                traversal%cell_data(traversal%i_cell_data_index)%depth = element%cell%geometry%i_depth
                traversal%cell_data(traversal%i_cell_data_index)%u = u
                traversal%cell_data(traversal%i_cell_data_index)%refinement = element%cell%geometry%refinement
#           endif

			if (l_second_order) then
				forall (i = 1 : 6)
					traversal%point_data(point_data_indices(i))%coords = samoa_barycentric_to_world_point(element%transform_data, samoa_basis_p_get_dof_coords(i))
					traversal%point_data(point_data_indices(i))%p = p(i)
					traversal%point_data(point_data_indices(i))%S = saturation(i)
				end forall

				traversal%i_connectivity(7 * traversal%i_cell_data_index - 6) = 6
				traversal%i_connectivity(7 * traversal%i_cell_data_index - 5 : 7 * traversal%i_cell_data_index) = point_data_indices([ 1, 2, 3, 6, 4, 5 ]) - 1
			else
				forall (i = 1 : 3)
					traversal%point_data(point_data_indices(i))%coords = samoa_barycentric_to_world_point(element%transform_data, samoa_basis_p_get_dof_coords(i))
					traversal%point_data(point_data_indices(i))%p = p(i)
					traversal%point_data(point_data_indices(i))%S = saturation(i)
				end forall

				traversal%i_connectivity(4 * traversal%i_cell_data_index - 3) = 3
				traversal%i_connectivity(4 * traversal%i_cell_data_index - 2 : 4 * traversal%i_cell_data_index) = point_data_indices(1:3) - 1
			end if

			traversal%i_cell_data_index = traversal%i_cell_data_index + 1
		end subroutine

		subroutine node_first_touch_op_array(traversal, section, nodes)
 			type(t_darcy_vtk_output_traversal), intent(inout)	:: traversal
 			type(t_grid_section), intent(in)				    :: section
			type(t_node_data), intent(inout)				    :: nodes(:)

			integer (kind = GRID_SI)						    :: i, j

            do j = 1, size(nodes)
                do i = 1, _DARCY_LAYERS + 1
                    call pre_dof_op(traversal%i_point_data_index, nodes(j)%data_pers%r(i))
                end do
            end do
		end subroutine

		subroutine node_first_touch_op_scalar(traversal, section, node)
 			type(t_darcy_vtk_output_traversal), intent(inout)	:: traversal
 			type(t_grid_section), intent(in)				    :: section
			type(t_node_data), intent(inout)				    :: node

			integer (kind = GRID_SI)						    :: i

			do i = 1, _DARCY_LAYERS + 1
				call pre_dof_op(traversal%i_point_data_index, node%data_pers%r(i))
			end do
		end subroutine

		elemental subroutine pre_dof_op(i_point_data_index, r)
 			integer(kind = GRID_SI), intent(inout)			:: i_point_data_index
			real(kind = GRID_SR), intent(out)				:: r

			r = real(i_point_data_index, GRID_SR)
			i_point_data_index = i_point_data_index + 1
		end subroutine

		pure subroutine scatter_traversals_op(grid_data, section_data)
            type(num_traversal_data), intent(in)	        :: grid_data
            type(num_traversal_data), intent(inout)		    :: section_data(:)

            integer (kind = GRID_SI)                        :: i_section

            do i_section = 1, size(section_data)
                write(section_data(i_section)%s_file_stamp, '(A, A, I0)') grid_data%s_file_stamp, "_", i_section
            end do

            section_data%i_output_iteration = grid_data%i_output_iteration
        end subroutine
	END MODULE
#endif
