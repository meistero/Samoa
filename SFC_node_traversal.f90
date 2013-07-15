! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


!*****************************************************************
! MODULE SFC_node_traversal: initialize node datastructures from an SFC_grid
!*****************************************************************

#include "Compilation_control.f90"

MODULE SFC_node_traversal
	use SFC_edge_traversal

	implicit none

	private
	public recursive_traversal_init

	CONTAINS

	subroutine recursive_traversal_init(triangles, num_coarse_triangles, grid)
		type(triangle_tree), dimension(:), intent(inout), target	:: triangles		! initial triangle array
		integer (kind = GRID_SI), intent(in)		:: num_coarse_triangles
		type(t_grid), intent(inout)			        :: grid
		integer (kind = GRID_SI)            		:: i_section, i_thread

		! local variables
        type(t_edge_data)			                :: last_crossed_edge_data

		if (rank_MPI == 0) then
            grid%sections%elements%t_global_data = grid%t_global_data
            i_thread = 1 + omp_get_thread_num()

            do i_section = 1, size(grid%sections%elements)
                call recursive_traversal_init_section(triangles, num_coarse_triangles, grid%threads%elements(i_thread), grid%sections%elements(i_section))

                call grid%threads%elements(i_thread)%edges_stack(RED)%pop(last_crossed_edge_data)
                call find_section_boundary_elements(grid%threads%elements(i_thread), grid%sections%elements(i_section), grid%sections%elements(i_section)%cells%i_current_element, last_crossed_edge_data)
            end do

            call update_distances(grid)
            call update_neighbors(grid, grid)
        end if

        call grid%reverse()
 	end subroutine

	subroutine recursive_traversal_init_section(triangles, num_coarse_triangles, thread, section)
		type(triangle_tree), intent(inout), target	:: triangles(:)		! initial triangle array
		integer (kind = GRID_SI), intent(in)		:: num_coarse_triangles
		type(t_grid_thread), intent(inout)			:: thread
		type(t_grid_section), intent(inout)			:: section

		type(fem_triangle)							:: first_fem_triangle
		integer (kind = GRID_SI)					:: i, i_color
		integer (kind = GRID_SI)	                :: i_triang


#	    if (_DEBUG_LEVEL > 3)
            _log_write(4, '(X, A)') "Initial section state :"
            call section%print()
#	    endif

        i_triang = 1

		DO i = 1, num_coarse_triangles
			triangles(i_triang)%i_output_index = triangles(i_triang)%i_input_min - 1

			call read_first_fem_triangle(triangles(i_triang), first_fem_triangle)
 			call sfc_init_recursion(triangles(i_triang), first_fem_triangle, thread, section)

			i_triang = i_triang + 1
		end DO

#	    if (_DEBUG_LEVEL > 5)
            _log_write(6, '(2X, A)') "initial output cells:"
            do i = lbound(section%cells%elements, 1), ubound(section%cells%elements, 1)
                _log_write(6, '(3X, I0, X, A)') i, section%cells%elements(i)%to_string()
            end do
            _log_write(6, '(A)') ""
#	    endif
	end subroutine

	recursive subroutine sfc_init_recursion(p_tree, fem_tri, thread, section)
		type(triangle_tree), intent(inout), target						:: p_tree 		! intent(inout) (for the contents of the target)
		type(fem_triangle), intent(inout)								:: fem_tri
		type(t_grid_thread), intent(inout)			:: thread
		type(t_grid_section), intent(inout)								:: section

		!local variables
		type(fem_triangle)												:: fem_child
		type(t_cell_stream_data), pointer								:: p_cell_data

		if (is_leaf_triangle(fem_tri)) then		! initialize EDGE structures
			p_cell_data => section%cells%next()
			p_cell_data%fine_triangle = fem_tri

			call sfc_init(thread, section, p_cell_data, fem_tri)
		else
			call read_first_child(p_tree, fem_tri, fem_child)
			call sfc_init_recursion(p_tree, fem_child, thread, section)

			call read_second_child(p_tree, fem_tri, fem_child)
			call sfc_init_recursion(p_tree, fem_child, thread, section)
		endif

		call write_triangle_out(p_tree, fem_tri%triangle)
	end subroutine

	subroutine sfc_init(thread, section, cell_data, fem_tri)
		type(t_grid_thread), intent(inout)			                    :: thread
		type(t_grid_section), intent(inout)								:: section
		type(t_cell_stream_data), intent(inout)							:: cell_data
		type(fem_triangle), intent(in)									:: fem_tri

		integer (kind = 1)                                              :: edge_depths(3)

		type(t_crossed_edge_stream_data), pointer						:: previous_edge, next_edge
		type(t_edge_data)												:: color_edge
		type(t_edge_data), pointer										:: boundary_edge
		type(t_node_data), pointer										:: color_node_in, color_node_out, transfer_node
		integer (KIND = 1)												:: i_previous_edge_type, i_color_edge_type, i_next_edge_type
		integer (KIND = 1)												:: i_previous_edge_index, i_color_edge_index, i_next_edge_index
		integer (KIND = 1)												:: i_color_node_in_index, i_transfer_node_index, i_color_node_out_index
		logical (KIND = GRID_SL)										:: l_color

		l_color = cell_data%l_color_edge_color

		call cell_data%get_edge_types(i_previous_edge_type, i_color_edge_type, i_next_edge_type)
		call cell_data%get_edge_indices(i_previous_edge_index, i_color_edge_index, i_next_edge_index)
		call cell_data%get_node_indices(i_color_node_out_index, i_transfer_node_index, i_color_node_in_index)
		edge_depths = [cell_data%i_depth, cell_data%i_depth - 1_1, cell_data%i_depth]

		select case (i_previous_edge_type)
			case (OLD)
				previous_edge => section%crossed_edges_out%current()
				transfer_node => thread%nodes_stack(.not. l_color)%current()
			case (OLD_BND)
				boundary_edge => section%boundary_edges(RED)%next()
                boundary_edge%depth = edge_depths(i_previous_edge_index)
				previous_edge => boundary_edge%t_crossed_edge_stream_data

                color_node_out => thread%nodes_stack(l_color)%push()
                call section%boundary_nodes(l_color)%read(color_node_out)

                transfer_node => thread%nodes_stack(.not. l_color)%push()
                call section%boundary_nodes(.not. l_color)%read(transfer_node)
		end select

		select case (i_color_edge_type)
			case (OLD)
				call thread%edges_stack(l_color)%pop(color_edge)

				color_node_out => thread%nodes_stack(l_color)%pop()
				color_node_in => thread%nodes_stack(l_color)%current()
			case (NEW)
				call thread%edges_stack(l_color)%push(color_edge)

				color_node_out => thread%nodes_stack(l_color)%current()
				color_node_in => thread%nodes_stack(l_color)%push()
			case (OLD_BND)
                color_edge%depth = edge_depths(i_color_edge_index)

				color_node_out => section%boundary_nodes(l_color)%current()
				call thread%nodes_stack(l_color)%pop(color_node_out)

				color_node_in => thread%nodes_stack(l_color)%push()
				call section%boundary_nodes(l_color)%read(color_node_in)
			case (NEW_BND)
                color_edge%depth = edge_depths(i_color_edge_index)
				call thread%edges_stack(l_color)%push(color_edge)

				color_node_out => thread%nodes_stack(l_color)%current()
				color_node_in => thread%nodes_stack(l_color)%push()
		end select

		select case (i_next_edge_type)
			case (NEW)
				next_edge => section%crossed_edges_out%next()
			case (NEW_BND)
				boundary_edge => thread%edges_stack(RED)%push()
				boundary_edge%depth = edge_depths(i_next_edge_index)
				next_edge => boundary_edge%t_crossed_edge_stream_data
		end select

		transfer_node%position = fem_tri%r_coords(:, i_transfer_node_index)
		color_node_out%position = fem_tri%r_coords(:, i_color_node_out_index)
		color_node_in%position = fem_tri%r_coords(:, i_color_node_in_index)

		select case (i_color_edge_type)
			case (OLD)
				call section%color_edges_out%write(color_edge%t_color_edge_stream_data)
				call section%nodes_out%write(color_node_out%t_node_stream_data)
			case (OLD_BND)
				call section%boundary_edges(l_color)%write(color_edge)
		end select

		call cell_data%reverse()
	end subroutine
end MODULE
