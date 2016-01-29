! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


!*****************************************************************
! MODULE SFC_node_traversal: initialize node datastructures from an SFC_data_types
!*****************************************************************

#include "Compilation_control.f90"

MODULE SFC_node_traversal
	use SFC_edge_traversal

	implicit none

	private
	public init_grid

	CONTAINS

	subroutine init_grid(grid, depth)
		type(t_grid), intent(inout)			    :: grid
		integer (kind = BYTE), intent(in), optional    :: depth
		integer (kind = GRID_SI)            		:: i_section, i_thread

		! local variables
		type(t_section_info)           	            :: section_descriptor
		type(t_section_info_list)           	    :: section_descriptors
        type(t_edge_data)			                :: last_crossed_edge_data
		integer                                     :: start_depth

        !set start depth to the passed argument if available, otherwise to 0
        start_depth = 0

		if (present(depth)) then
            assert_ge(depth, 0)
            start_depth = depth
        end if

		!the start grid belongs to rank 0 and will be distributed during runtime
		if (rank_MPI == 0) then
			section_descriptor = t_section_info(&
                index = 1, &
                i_load = 0, &
				i_cells = 2 ** (start_depth + 1), &
				i_stack_nodes = [2 ** ((start_depth + 4) / 2), 2 ** ((start_depth + 3) / 2)], &
				i_stack_edges = [2 ** ((start_depth + 4) / 2) - 1, 2 ** ((start_depth + 3) / 2) - 1], &
				i_boundary_edges = 0, &
				i_boundary_nodes = 0, &
				i_comms = 0)

			call section_descriptor%estimate_bounds()

            !add only section to the section desctiptor list
            call section_descriptors%add(section_descriptor)
		endif

        grid%start_distance = 0
        grid%min_distance = 0
		grid%end_distance = 0

		!$omp parallel
        call grid%create(section_descriptors, section_descriptor%i_stack_nodes)
		!$omp end parallel

		if (rank_MPI == 0) then
            grid%sections%elements%t_global_data = grid%t_global_data
            i_thread = 1 + omp_get_thread_num()

            do i_section = 1, size(grid%sections%elements)
                call init_section(grid%threads%elements(i_thread), grid%sections%elements(i_section), start_depth)

                call grid%threads%elements(i_thread)%edges_stack(RED)%pop_data(last_crossed_edge_data)
                call find_section_boundary_elements(grid%threads%elements(i_thread), grid%sections%elements(i_section), grid%sections%elements(i_section)%cells%i_current_element, last_crossed_edge_data)
            end do

            call update_distances(grid)
            call update_neighbors(grid, grid)
        end if

        call grid%reverse()
 	end subroutine

	subroutine init_section(thread, section, start_depth)
		type(t_grid_thread), intent(inout)			:: thread
		type(t_grid_section), intent(inout)			:: section
		integer                                     :: start_depth

		type(fine_triangle), parameter				:: first_cell = fine_triangle(i_entity_types = FIRST_OLD_BND, i_depth = 0, &
            refinement = 0, i_plotter_type = -2, i_turtle_type = K, i_color_edge_color = RED)

		type(fine_triangle)							:: second_cell = fine_triangle(i_entity_types = LAST_NEW_BND, i_depth = 0, &
            refinement = 0, i_plotter_type = -6, i_turtle_type = H, i_color_edge_color = RED)

        real (kind = GRID_SR), parameter            :: first_coords(2, 3) = reshape([[0.0_GRID_SR, 0.0_GRID_SR], [1.0_GRID_SR, 0.0_GRID_SR], [1.0_GRID_SR, 1.0_GRID_SR]], shape(first_coords))
        real (kind = GRID_SR), parameter            :: second_coords(2, 3) = reshape([[1.0_GRID_SR, 1.0_GRID_SR], [0.0_GRID_SR, 1.0_GRID_SR], [0.0_GRID_SR, 0.0_GRID_SR]], shape(second_coords))

        integer :: i

#	    if (_DEBUG_LEVEL > 3)
            _log_write(4, '(X, A)') "Initial section state :"
            call section%print()
#	    endif

        !use hardcoded cells to form a unit square

        call sfc_init_recursion(thread, section, first_cell, first_coords, start_depth)
        call sfc_init_recursion(thread, section, second_cell, second_coords, start_depth)

#	    if (_DEBUG_LEVEL > 5)
            _log_write(6, '(2X, A)') "initial output cells:"
            do i = lbound(section%cells%elements, 1), ubound(section%cells%elements, 1)
                _log_write(6, '(3X, I0, X, A)') i, section%cells%elements(i)%to_string()
            end do
            _log_write(6, '(A)') ""
#	    endif
	end subroutine

	recursive subroutine sfc_init_recursion(thread, section, cell, coords, levels)
		type(t_grid_thread), intent(inout)			                    :: thread
		type(t_grid_section), intent(inout)								:: section
		type(fine_triangle), intent(in)                                 :: cell
        real (kind = GRID_SR), intent(in)							    :: coords(2, 3)
		integer                                                         :: levels

		!local variables
		type(t_cell_stream_data), pointer								:: p_cell_data
		type(fine_triangle)								                :: first_cell, second_cell
		real (kind = GRID_SR)								            :: first_coords(2, 3), second_coords(2, 3)

		if (levels == 0) then
			p_cell_data => section%cells%next()
			p_cell_data%fine_triangle = cell

			call sfc_init(thread, section, p_cell_data, coords)
		else
			call create_child_cells(cell, first_cell, second_cell)
			first_coords = reshape([coords(:, 1), 0.5_GRID_SR * (coords(:, 1) + coords(:, 3)), coords(:, 2)], shape(first_coords))
			second_coords = reshape([coords(:, 2), 0.5_GRID_SR * (coords(:, 1) + coords(:, 3)), coords(:, 3)], shape(second_coords))

			call sfc_init_recursion(thread, section, first_cell, first_coords, levels - 1)
			call sfc_init_recursion(thread, section, second_cell, second_coords, levels - 1)
		endif
	end subroutine

	subroutine sfc_init(thread, section, cell_data, coords)
		type(t_grid_thread), intent(inout)			                    :: thread
		type(t_grid_section), intent(inout)								:: section
		type(t_cell_stream_data), intent(inout)							:: cell_data
		real (kind = GRID_SR), intent(in)								:: coords(2, 3)

		integer (kind = BYTE)                                              :: edge_depths(3)

		type(t_crossed_edge_stream_data), pointer						:: previous_edge, next_edge
		type(t_edge_data)												:: color_edge
		type(t_edge_data), pointer										:: boundary_edge
		type(t_node_data), pointer										:: color_node_in, color_node_out, transfer_node
		integer (kind = BYTE)												:: i_previous_edge_type, i_color_edge_type, i_next_edge_type
		integer (kind = BYTE)												:: i_previous_edge_index, i_color_edge_index, i_next_edge_index
		integer (kind = BYTE)												:: i_color_node_in_index, i_transfer_node_index, i_color_node_out_index
		integer (kind = BYTE)										        :: i_color_edge_color

		i_color_edge_color = cell_data%i_color_edge_color

		call cell_data%get_edge_types(i_previous_edge_type, i_color_edge_type, i_next_edge_type)
		call cell_data%get_edge_indices(i_previous_edge_index, i_color_edge_index, i_next_edge_index)
		call cell_data%get_node_indices(i_color_node_out_index, i_transfer_node_index, i_color_node_in_index)
		edge_depths = [cell_data%i_depth, cell_data%i_depth - 1_1, cell_data%i_depth]

		select case (i_previous_edge_type)
			case (OLD)
				previous_edge => section%crossed_edges_out%current()

				if (thread%nodes_stack(RED + GREEN - i_color_edge_color)%is_empty()) then
				    transfer_node => section%boundary_nodes(RED + GREEN - i_color_edge_color)%current()
				else
				    transfer_node => thread%nodes_stack(RED + GREEN - i_color_edge_color)%current()
				end if
			case (OLD_BND)
				boundary_edge => section%boundary_edges(RED)%next()
                boundary_edge%depth = edge_depths(i_previous_edge_index)
				previous_edge => boundary_edge%t_crossed_edge_stream_data

                color_node_out => section%boundary_nodes(i_color_edge_color)%next()
                transfer_node => section%boundary_nodes(RED + GREEN - i_color_edge_color)%next()
		end select

		select case (i_color_edge_type)
			case (OLD)
				call thread%edges_stack(i_color_edge_color)%pop_data(color_edge)

				color_node_out => thread%nodes_stack(i_color_edge_color)%pop()

				if (thread%nodes_stack(i_color_edge_color)%is_empty()) then
				    color_node_in => section%boundary_nodes(i_color_edge_color)%current()
				else
				    color_node_in => thread%nodes_stack(i_color_edge_color)%current()
				end if
			case (NEW)
				call thread%edges_stack(i_color_edge_color)%push_data(color_edge)

				if (thread%nodes_stack(i_color_edge_color)%is_empty()) then
				    color_node_out => section%boundary_nodes(i_color_edge_color)%current()
				else
				    color_node_out => thread%nodes_stack(i_color_edge_color)%current()
				end if

				color_node_in => thread%nodes_stack(i_color_edge_color)%push()
			case (OLD_BND)
                color_edge%depth = edge_depths(i_color_edge_index)

				color_node_out => section%boundary_nodes(i_color_edge_color)%current()
				color_node_in => section%boundary_nodes(i_color_edge_color)%next()
			case (NEW_BND)
                color_edge%depth = edge_depths(i_color_edge_index)
				call thread%edges_stack(i_color_edge_color)%push_data(color_edge)

				if (thread%nodes_stack(i_color_edge_color)%is_empty()) then
				    color_node_out => section%boundary_nodes(i_color_edge_color)%current()
				else
				    color_node_out => thread%nodes_stack(i_color_edge_color)%current()
				end if

				color_node_in => thread%nodes_stack(i_color_edge_color)%push()
		end select

		select case (i_next_edge_type)
			case (NEW)
				next_edge => section%crossed_edges_out%next()
			case (NEW_BND)
				boundary_edge => thread%edges_stack(RED)%push()
				boundary_edge%depth = edge_depths(i_next_edge_index)
				next_edge => boundary_edge%t_crossed_edge_stream_data
		end select

		transfer_node%position = coords(:, i_transfer_node_index)
		color_node_out%position = coords(:, i_color_node_out_index)
		color_node_in%position = coords(:, i_color_node_in_index)

		select case (i_color_edge_type)
			case (OLD)
				call section%color_edges_out%write(color_edge%t_color_edge_stream_data)
				call section%nodes_out%write(color_node_out%t_node_stream_data)
			case (OLD_BND)
				call section%boundary_edges(i_color_edge_color)%write(color_edge)
		end select

		call cell_data%reverse()
	end subroutine

	subroutine create_child_cells(parent_cell, first_child_cell, second_child_cell)
        type(fine_triangle), intent(in)						:: parent_cell
        type(fine_triangle), intent(inout)					:: first_child_cell, second_child_cell

        integer(kind = BYTE), dimension(3, 2), parameter		:: i_turtle_child_type = reshape([ H, H, V, V, K, K ], [ 3, 2 ])
        integer(kind = BYTE), dimension(-8:8), parameter 	    :: i_plotter_child_type = [ 3, 2, 1, 8, 7, 6, 5, 4,     0,  -6, -7, -8, -1, -2, -3, -4, -5 ]
        integer(kind = BYTE)			 						:: i_previous_edge_type, i_color_edge_type, i_next_edge_type
        integer(kind = BYTE)								 	:: i

        first_child_cell%i_turtle_type = i_turtle_child_type(parent_cell%i_turtle_type, 1)
        second_child_cell%i_turtle_type = i_turtle_child_type(parent_cell%i_turtle_type, 2)

        !check for correctness
        assert_eq(parent_cell%i_turtle_type, ieor(first_child_cell%i_turtle_type, second_child_cell%i_turtle_type))
        assert_gt(first_child_cell%i_turtle_type, second_child_cell%i_turtle_type)

        first_child_cell%i_plotter_type = i_plotter_child_type(parent_cell%i_plotter_type)
        first_child_cell%i_depth = parent_cell%i_depth + 1
        second_child_cell%i_plotter_type = -i_plotter_child_type(-parent_cell%i_plotter_type)
        second_child_cell%i_depth = parent_cell%i_depth + 1

        call parent_cell%get_edge_types(i_previous_edge_type, i_color_edge_type, i_next_edge_type)

        select case (parent_cell%i_turtle_type)
            case (K)
                call first_child_cell%set_edge_types(i_previous_edge_type, i_next_edge_type, int(NEW, 1))
                call second_child_cell%set_edge_types(int(OLD, 1), i_color_edge_type, i_next_edge_type)
                first_child_cell%i_color_edge_color = RED + GREEN - parent_cell%i_color_edge_color
                second_child_cell%i_color_edge_color = parent_cell%i_color_edge_color
            case (V)
                call first_child_cell%set_edge_types(i_previous_edge_type, i_color_edge_type, int(NEW, 1))
                call second_child_cell%set_edge_types(int(OLD, 1), i_color_edge_type, i_next_edge_type)
                first_child_cell%i_color_edge_color = parent_cell%i_color_edge_color
                second_child_cell%i_color_edge_color = parent_cell%i_color_edge_color
            case (H)
                call first_child_cell%set_edge_types(i_previous_edge_type, i_color_edge_type, int(NEW, 1))
                call second_child_cell%set_edge_types(int(OLD, 1), i_previous_edge_type, i_next_edge_type)
                first_child_cell%i_color_edge_color = parent_cell%i_color_edge_color
                second_child_cell%i_color_edge_color = RED + GREEN - parent_cell%i_color_edge_color
        end select
    end subroutine
end MODULE
