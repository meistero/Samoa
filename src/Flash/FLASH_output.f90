#include "Compilation_control.f90"

#if defined(_FLASH)
	MODULE FLASH_output
		use LIB_VTK_IO

		use SFC_edge_traversal

		use Samoa_FLASH
		use FLASH_euler_timestep

		implicit none

		!> Output point data
		type t_output_point_data
			type(t_state)											:: Q
			real (kind = GRID_SR), dimension(2)						:: coords		!< position
		end type

		!> Output cell data
		type t_output_cell_data
			type(t_state)											:: Q
			integer (kind = GRID_SI)								:: depth
			integer (kind = GRID_SI)								:: refinement
		end type

        type num_traversal_data
            type(t_output_point_data), allocatable		            :: point_data(:)
            type(t_output_cell_data), allocatable			        :: cell_data(:)
            character(len=64)							            :: s_file_stamp

            integer (kind = GRID_SI)								:: i_output_iteration
            integer (kind = GRID_SI)								:: i_point_data_index
            integer (kind = GRID_SI)								:: i_cell_data_index
        end type

		integer, PARAMETER											:: i_element_order = 0

		type(t_gv_Q)												:: gv_Q

#		define _GT_NAME								t_FLASH_output_traversal

#		define _GT_EDGES
#		define _GT_REFINEMENTS

#		define _GT_PRE_TRAVERSAL_OP					pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP				post_traversal_op

#		define _GT_ELEMENT_OP						element_op

#		define _GT_CELL_TO_EDGE_OP				    cell_to_edge_op

#		include "SFC_generic_traversal_ringbuffer.f90"


		subroutine pre_traversal_op(traversal, section)
			type(t_FLASH_output_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: section

			type(t_section_info)                                           :: grid_info
			integer (kind = GRID_SI)									:: i_error, i_cells, i_points

			_log_write(1, '(A, I0)') " FLASH: output step ", traversal%i_output_iteration

            grid_info = section%get_info()
			i_cells = grid_info%i_cells

			if (i_element_order > 1) then
				i_points = 6 * i_cells
			else
				i_points = 3 * i_cells
			end if

			allocate(traversal%cell_data(i_cells), stat = i_error); assert_eq(i_error, 0)
			allocate(traversal%point_data(i_points), stat = i_error); assert_eq(i_error, 0)

			traversal%i_cell_data_index = 1
			traversal%i_point_data_index = 1
		end subroutine

		subroutine post_traversal_op(traversal, grid)
			type(t_FLASH_output_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: grid

			integer (kind = GRID_SI), dimension(:), allocatable			:: i_offsets
			integer (kind = GRID_SI), dimension(:), allocatable			:: i_types
			integer (kind = GRID_SI), dimension(:), allocatable			:: i_connectivity
			real (kind = GRID_SR), dimension(:, :), allocatable			:: r_velocity
			real (kind = GRID_SR), dimension(:), allocatable			:: r_empty
            type(t_vtk_writer)                                          :: vtk

			type(t_section_info)                                           :: grid_info
			integer (kind = GRID_SI)									:: i_error, i_cells, i_points
			character (len = 64)										:: s_file_name
			integer(4)													:: e_io, i

            grid_info = grid%get_info()
			i_cells = grid_info%i_cells

			if (i_element_order > 1) then
				i_points = 6 * i_cells
				allocate(i_connectivity(7 * i_cells), stat = i_error); assert_eq(i_error, 0)
			else
				i_points = 3 * i_cells
				allocate(i_connectivity(4 * i_cells), stat = i_error); assert_eq(i_error, 0)
			end if

			allocate(i_offsets(i_cells), stat = i_error); assert_eq(i_error, 0)
			allocate(i_types(i_cells), stat = i_error); assert_eq(i_error, 0)
			allocate(r_velocity(2, max(i_cells, i_points)), stat = i_error); assert_eq(i_error, 0)
			allocate(r_empty(max(i_cells, i_points)), stat = i_error); assert_eq(i_error, 0)

			r_empty = 0.0_GRID_SR

			if (i_element_order > 1) then
				i_types = 22_1

				forall (i = 1 : i_cells)
					i_offsets(i) = 6_GRID_SI * i
					i_connectivity(7 * i - 6 : 7 * i) = [6, 6 * i - 6, 6 * i - 5, 6 * i - 4, 6 * i - 3, 6 * i - 2, 6 * i - 1]
				end forall
			else
				i_types = 5_1

				forall (i = 1 : i_cells)
					i_offsets(i) = 3_GRID_SI * i
					i_connectivity(4 * i - 3 : 4 * i) = [3, 3 * i - 3, 3 * i - 2, 3 * i - 1]
				end forall
			end if

			write (s_file_name, "(A, A, I0, A)") TRIM(traversal%s_file_stamp), "_", traversal%i_output_iteration, ".vtk"

			e_io = vtk%VTK_INI('ASCII', s_file_name, 'FLASH', 'UNSTRUCTURED_GRID')
				e_io = vtk%VTK_GEO(i_points, traversal%point_data%coords(1), traversal%point_data%coords(2), r_empty(1:i_points))

				e_io = vtk%VTK_CON(i_cells, i_connectivity, i_types)

				e_io = vtk%VTK_DAT(i_points, 'node')

				if (i_element_order > 0) then
					e_io = vtk%VTK_VAR(i_points, 'water_height', traversal%point_data%Q%h)
					e_io = vtk%VTK_VAR(i_points, 'bathymetry', traversal%point_data%Q%b)

					r_velocity(1, 1:i_points) = traversal%point_data%Q%p(1) / (traversal%point_data%Q%h - traversal%point_data%Q%b)
					r_velocity(2, 1:i_points) = traversal%point_data%Q%p(2) / (traversal%point_data%Q%h - traversal%point_data%Q%b)
					e_io = vtk%VTK_VAR('vect', i_points, 'velocity', r_velocity(1, 1:i_points), r_velocity(2, 1:i_points), r_empty(1:i_points))
				endif

				e_io = vtk%VTK_DAT(i_cells, 'cell')

				if (i_element_order == 0) then
					e_io = vtk%VTK_VAR(i_cells, 'water_height', traversal%cell_data%Q%h)
					e_io = vtk%VTK_VAR(i_cells, 'bathymetry', traversal%cell_data%Q%b)

					r_velocity(1, 1:i_cells) = traversal%cell_data%Q%p(1) / (traversal%cell_data%Q%h - traversal%cell_data%Q%b)
					r_velocity(2, 1:i_cells) = traversal%cell_data%Q%p(2) / (traversal%cell_data%Q%h - traversal%cell_data%Q%b)
					e_io = vtk%VTK_VAR('vect', i_cells, 'velocity', r_velocity(1, 1:i_cells), r_velocity(2, 1:i_cells), r_empty(1:i_cells))
				endif

				e_io = vtk%VTK_VAR(i_cells, 'grid_depth', traversal%cell_data%depth)
				e_io = vtk%VTK_VAR(i_cells, 'refinement_flag', traversal%cell_data%refinement)
			e_io = vtk%VTK_END()

			deallocate(i_offsets, stat = i_error); assert_eq(i_error, 0)
			deallocate(i_types, stat = i_error); assert_eq(i_error, 0)
			deallocate(i_connectivity, stat = i_error); assert_eq(i_error, 0)
			deallocate(r_velocity, stat = i_error); assert_eq(i_error, 0)
			deallocate(r_empty, stat = i_error); assert_eq(i_error, 0)

			deallocate(traversal%cell_data, stat = i_error); assert_eq(i_error, 0)
			deallocate(traversal%point_data, stat = i_error); assert_eq(i_error, 0)

			traversal%i_output_iteration = traversal%i_output_iteration + 1
		end subroutine


		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, grid, element)
			type(t_FLASH_output_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: grid
			type(t_element_base), intent(inout)					:: element

			!local variables

			integer (kind = GRID_SI)							:: i
			real (kind = GRID_SR), parameter, dimension(2, 6)	:: r_test_points = reshape((/ 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.5, 0.5, 0.0, 0.5, 0.5 /), (/ 2, 6 /))
			real (kind = GRID_SR), parameter, dimension(2)		:: r_test_point0 = [1.0_GRID_SR/3.0_GRID_SR, 1.0_GRID_SR/3.0_GRID_SR]

			type(t_state), dimension(_FLASH_CELL_SIZE)			:: Q
			type(t_state), dimension(6)							:: Q_test

			call gv_Q%read(element, Q)

			traversal%cell_data(traversal%i_cell_data_index)%depth = element%cell%geometry%i_depth
			traversal%cell_data(traversal%i_cell_data_index)%refinement = element%cell%geometry%refinement

			select case (i_element_order)
				case (2)
					forall (i = 1 : 6)
						traversal%point_data(traversal%i_point_data_index + i - 1)%coords = samoa_barycentric_to_world_point(element%transform_data, r_test_points(:, i))
						traversal%point_data(traversal%i_point_data_index + i - 1)%Q%h = t_basis_Q_eval(r_test_points(:, i), Q%h)
						traversal%point_data(traversal%i_point_data_index + i - 1)%Q%b = t_basis_Q_eval(r_test_points(:, i), Q%b)
						traversal%point_data(traversal%i_point_data_index + i - 1)%Q%p(1) = t_basis_Q_eval(r_test_points(:, i), Q%p(1))
						traversal%point_data(traversal%i_point_data_index + i - 1)%Q%p(2) = t_basis_Q_eval(r_test_points(:, i), Q%p(2))
					end forall

					traversal%i_point_data_index = traversal%i_point_data_index + 6
				case (1)
					forall (i = 1 : 3)
						traversal%point_data(traversal%i_point_data_index + i - 1)%coords = samoa_barycentric_to_world_point(element%transform_data, r_test_points(:, i))
						traversal%point_data(traversal%i_point_data_index + i - 1)%Q%h = t_basis_Q_eval(r_test_points(:, i), Q%h)
						traversal%point_data(traversal%i_point_data_index + i - 1)%Q%b = t_basis_Q_eval(r_test_points(:, i), Q%b)
						traversal%point_data(traversal%i_point_data_index + i - 1)%Q%p(1) = t_basis_Q_eval(r_test_points(:, i), Q%p(1))
						traversal%point_data(traversal%i_point_data_index + i - 1)%Q%p(2) = t_basis_Q_eval(r_test_points(:, i), Q%p(2))
					end forall

					traversal%i_point_data_index = traversal%i_point_data_index + 3
				case (0)
					forall (i = 1 : 3)
						traversal%point_data(traversal%i_point_data_index + i - 1)%coords = samoa_barycentric_to_world_point(element%transform_data, r_test_points(:, i))
					end forall

					traversal%i_point_data_index = traversal%i_point_data_index + 3

					traversal%cell_data(traversal%i_cell_data_index)%Q%h = t_basis_Q_eval(r_test_point0, Q%h)
					traversal%cell_data(traversal%i_cell_data_index)%Q%b = t_basis_Q_eval(r_test_point0, Q%b)
					traversal%cell_data(traversal%i_cell_data_index)%Q%p(1) = t_basis_Q_eval(r_test_point0, Q%p(1))
					traversal%cell_data(traversal%i_cell_data_index)%Q%p(2) = t_basis_Q_eval(r_test_point0, Q%p(2))
			end select

			traversal%i_cell_data_index = traversal%i_cell_data_index + 1
		end subroutine
	END MODULE
#endif
