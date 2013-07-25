! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


!> Generic traversal template
!> Warning: this a template module body and requires preprocessor commands before inclusion.
!> Usage: Inside your module definition, insert
!>
!> The resulting method is defined as _GT_NAME
!> @author Oliver Meister


!Macro hell incoming:

!Recursively include this file until all (reasonable) variants of edge types are covered

#if .not. defined(_GT_PREVIOUS_EDGE_TYPE)
#	define _GT_PREVIOUS_EDGE_TYPE		_UNDEFINED
#	include __FILE__
#	undef _GT_PREVIOUS_EDGE_TYPE

#	define _GT_PREVIOUS_EDGE_TYPE		_OLD
#	define _GT_PREVIOUS_EDGE_SUFFIX		o
#	include __FILE__
#	undef _GT_PREVIOUS_EDGE_SUFFIX
#	undef _GT_PREVIOUS_EDGE_TYPE

#	define _GT_PREVIOUS_EDGE_TYPE		_OLD_BND
#	define _GT_PREVIOUS_EDGE_SUFFIX		d
#	include __FILE__
#	undef _GT_PREVIOUS_EDGE_SUFFIX
#	undef _GT_PREVIOUS_EDGE_TYPE
#else

#if .not. defined(_GT_COLOR_EDGE_TYPE)
#	define _GT_COLOR_EDGE_TYPE			_UNDEFINED
#	include __FILE__
#	undef _GT_COLOR_EDGE_TYPE

#	define _GT_COLOR_EDGE_TYPE			_OLD
#	define _GT_COLOR_EDGE_SUFFIX		o
#	include __FILE__
#	undef _GT_COLOR_EDGE_SUFFIX
#	undef _GT_COLOR_EDGE_TYPE

#	define _GT_COLOR_EDGE_TYPE			_NEW
#	define _GT_COLOR_EDGE_SUFFIX		n
#	include __FILE__
#	undef _GT_COLOR_EDGE_SUFFIX
#	undef _GT_COLOR_EDGE_TYPE

#	define _GT_COLOR_EDGE_TYPE			_OLD_BND
#	define _GT_COLOR_EDGE_SUFFIX		d
#	include __FILE__
#	undef _GT_COLOR_EDGE_SUFFIX
#	undef _GT_COLOR_EDGE_TYPE

#	define _GT_COLOR_EDGE_TYPE			_NEW_BND
#	define _GT_COLOR_EDGE_SUFFIX		b
#	include __FILE__
#	undef _GT_COLOR_EDGE_SUFFIX
#	undef _GT_COLOR_EDGE_TYPE

#else

#if .not. defined(_GT_NEXT_EDGE_TYPE)
#	define _GT_NEXT_EDGE_TYPE			_UNDEFINED
#	include __FILE__
#	undef _GT_NEXT_EDGE_TYPE

#	define _GT_NEXT_EDGE_TYPE			_NEW
#	define _GT_NEXT_EDGE_SUFFIX			n
#	include __FILE__
#	undef _GT_NEXT_EDGE_SUFFIX
#	undef _GT_NEXT_EDGE_TYPE

#	define _GT_NEXT_EDGE_TYPE			_NEW_BND
#	define _GT_NEXT_EDGE_SUFFIX			b
#	include __FILE__
#	undef _GT_NEXT_EDGE_SUFFIX
#	undef _GT_NEXT_EDGE_TYPE
#else

#if .not. defined(_GT_PREFIX)
#	define _OP0(x)						x
#	define _OP1(x, a)					_SUFFIX1(x, a)
#	define _OP2(x, a, b)				_SUFFIX2(x, a, b)
#	define _OP3(x, a, b, c)				_SUFFIX3(x, a, b, c)
#else
#	define _OP0(x)						_SUFFIX0E(x, _GT_PREFIX)
#	define _OP1(x, a)					_SUFFIX1E(x, _GT_PREFIX, a)
#	define _OP2(x, a, b)				_SUFFIX2E(x, _GT_PREFIX, a, b)
#	define _OP3(x, a, b, c)				_SUFFIX3E(x, _GT_PREFIX, a, b, c)
#endif

#define _SUFFIX1(x, a)					x ## _ ## a
#define _SUFFIX2(x, a, b)				x ## _ ## a ## b
#define _SUFFIX3(x, a, b, c)			x ## _ ## a ## b ## c

#define _SUFFIX0E(x, y)					_CONCAT0E(x, y)
#define _SUFFIX1E(x, y, a)				_CONCAT1E(x, y, a)
#define _SUFFIX2E(x, y, a, b)			_CONCAT2E(x, y, a, b)
#define _SUFFIX3E(x, y, a, b, c)		_CONCAT3E(x, y, a, b, c)

#define _CONCAT0E(x, y)					x ## _ ## y
#define _CONCAT1E(x, y, a)				x ## _ ## y ## _ ## a
#define _CONCAT2E(x, y, a, b)			x ## _ ## y ## _ ## a ## b
#define _CONCAT3E(x, y, a, b, c)		x ## _ ## y ## _ ## a ## b ## c

#define _UNDEFINED						-1
#define _OLD							0
#define _NEW							1
#define _OLD_BND					    2
#define _NEW_BND						3

! Short forms for node macros (enable/disable code execution based on preprocessor directives)

#if defined(_GT_NODES) .or. .not. defined(_GT_NO_COORDS)
#	define _GT_N(x)			x
#else
#	define _GT_N(x)
#endif

#if defined(_GT_NODES) .and. defined(_GT_INPUT)
#	define _GT_NR(x)		x
#else
#	define _GT_NR(x)
#endif

#if defined(_GT_NODES) .and. defined(_GT_OUTPUT)
#	define _GT_NW(x)		x
#else
#	define _GT_NW(x)
#endif

! Short forms for edge macros (enable/disable code execution based on preprocessor directives)

#if defined(_GT_EDGES)
#	define _GT_E(x)		x
#else
#	define _GT_E(x)
#endif

#if defined(_GT_EDGES) .and. defined(_GT_INPUT)
#	define _GT_ER(x)		x
#else
#	define _GT_ER(x)
#endif

#if defined(_GT_EDGES) .and. defined(_GT_OUTPUT)
#	define _GT_EW(x)		x
#else
#	define _GT_EW(x)
#endif

!For each skeleton operator there must be a boundary skeleton operator defined
#if defined(_GT_SKELETON_OP) .and. .not. (defined(_GT_BND_SKELETON_OP) .and. defined(_GT_CELL_TO_EDGE_OP) .and. defined(_GT_CELL_UPDATE_OP))
#	error "If a skeleton operator exists, a boundary skeleton operator, a cell-to-edge operator and a cell update operator must also be defined"
#endif

!implement this function if no edge types are set (which means only once)
#if (_GT_PREVIOUS_EDGE_TYPE .eq. _UNDEFINED) .and. (_GT_COLOR_EDGE_TYPE .eq. _UNDEFINED) .and. (_GT_NEXT_EDGE_TYPE .eq. _UNDEFINED)
	subroutine _OP0(pre_traversal_grid)(traversal, grid)
		type(_GT), intent(inout)					                :: traversal
		type(t_grid), intent(inout)					        		:: grid

#		if defined(_GT_PRE_TRAVERSAL_GRID_OP)
        	call _GT_PRE_TRAVERSAL_GRID_OP(traversal, grid)
#		endif
	end subroutine

	subroutine _OP0(post_traversal_grid)(traversal, grid)
		type(_GT), intent(inout)					                :: traversal
		type(t_grid), intent(inout)					        		:: grid

#		if defined(_GT_POST_TRAVERSAL_GRID_OP)
		    call _GT_POST_TRAVERSAL_GRID_OP(traversal, grid)
#		endif
	end subroutine

	subroutine _OP0(boundary_skeleton)(traversal, section)
		type(_GT), intent(inout)					                :: traversal
		type(t_grid_section), intent(inout)					        :: section

		integer (kind = GRID_SI)							        :: i_color
        integer (kind = GRID_SI)							        :: i_comm
        type(t_comm_interface), pointer                             :: comm

#       if defined(_GT_EDGES) .and. defined(_GT_SKELETON_OP)
            do i_color = RED, GREEN
                do i_comm = 1, size(section%comms(i_color)%elements)
                    comm => section%comms(i_color)%elements(i_comm)

                    if (comm%neighbor_rank .ge. 0 .and. comm%i_edges > 0) then
                        !the skeleton operator must be called for locally owned edges only.
                        !it is sufficient to check the first edge, since all edges of the comm are flagged equally
                        if (comm%p_local_edges(1)%owned_locally) then
                            call _GT_SKELETON_OP(traversal, section, comm%p_local_edges, comm%p_local_edges%rep, comm%p_neighbor_edges(comm%i_edges : 1 : -1)%rep, comm%p_local_edges%update, comm%p_neighbor_edges(comm%i_edges : 1 : -1)%update)
                        end if
                    else
                        call _GT_BND_SKELETON_OP(traversal, section, comm%p_local_edges, comm%p_local_edges%rep, comm%p_local_edges%update)
                    end if
                end do
            end do
#       endif
	end subroutine

	subroutine _OP0(pre_traversal)(traversal, section)
		type(_GT), intent(inout)					                :: traversal
		type(t_grid_section), intent(inout)					        :: section

		integer (kind = GRID_SI)							        :: i_color

#		if defined(_GT_PRE_TRAVERSAL_OP)
			call _GT_PRE_TRAVERSAL_OP(traversal, section)
#		endif

		do i_color = RED, GREEN
#			if defined(_GT_NODES) .and. defined(_GT_NODE_FIRST_TOUCH_OP)
                call _GT_NODE_FIRST_TOUCH_OP(traversal, section, section%boundary_nodes(i_color)%elements)
#			endif

#			if defined(_GT_EDGES) .and. defined(_GT_EDGE_FIRST_TOUCH_OP)
                call _GT_EDGE_FIRST_TOUCH_OP(traversal, section, section%boundary_edges(i_color)%elements)
#			endif
		end do
	end subroutine

	subroutine _OP0(post_traversal)(traversal, section)
		type(_GT), intent(inout)					        :: traversal
		type(t_grid_section), intent(inout)					:: section

		integer (kind = GRID_SI)							:: i_color, i_node, i_edge

        _log_write(5, '(3X, A, I0)'), "section: ", section%index

		do i_color = RED, GREEN
            _log_write(5, '(4X, A)'), color_to_char(i_color)

#			if defined(_GT_EDGES)
                _log_write(5, '(5X, A)'), "boundary edges:"

#			    if defined(_GT_EDGE_LAST_TOUCH_OP)
                    call _GT_EDGE_LAST_TOUCH_OP(traversal, section, section%boundary_edges(i_color)%elements)
#			    endif

				do i_edge = 1, size(section%boundary_edges(i_color)%elements)
                    if (section%boundary_edges(i_color)%elements(i_edge)%owned_globally) then
#			            if defined(_GT_EDGE_REDUCE_OP)
                            call _GT_EDGE_REDUCE_OP(traversal, section, section%boundary_edges(i_color)%elements(i_edge))
#			            endif
                    endif
               end do
#			endif

#			if defined(_GT_NODES)
                _log_write(5, '(5X, A)'), "boundary nodes:"

#			    if defined(_GT_NODE_LAST_TOUCH_OP)
                    call _GT_NODE_LAST_TOUCH_OP(traversal, section, section%boundary_nodes(i_color)%elements)
#               endif

                do i_node = 1, size(section%boundary_nodes(i_color)%elements)
                    if (section%boundary_nodes(i_color)%elements(i_node)%owned_globally) then
#			            if defined(_GT_NODE_REDUCE_OP)
                            call _GT_NODE_REDUCE_OP(traversal, section, section%boundary_nodes(i_color)%elements(i_node))
#           			endif
                    endif
				end do
#			endif
		end do

#		if defined(_GT_POST_TRAVERSAL_OP)
			call _GT_POST_TRAVERSAL_OP(traversal, section)
#		endif
	end subroutine

    elemental subroutine _OP0(set_stats)(traversal, section)
        type(_GT), intent(inout)			        :: traversal
        type(t_grid_section), intent(inout)			:: section

        type(t_section_info)               	        :: info

        info = section%get_capacity()

        traversal%current_stats%i_traversals = 1
        traversal%current_stats%i_traversed_cells = info%i_cells
       	traversal%current_stats%i_traversed_memory = sizeof(section%cells%elements(1)) * info%i_cells

#   	if defined(_GT_NODES)
            traversal%current_stats%i_traversed_nodes = info%i_nodes + sum(info%i_boundary_nodes)
            traversal%current_stats%i_traversed_memory = traversal%current_stats%i_traversed_memory + sizeof(section%boundary_nodes(RED)%elements(1)%t_node_stream_data) * (info%i_nodes + sum(info%i_boundary_nodes))
#   	endif

#   	if defined(_GT_EDGES)
            traversal%current_stats%i_traversed_edges = info%i_crossed_edges + info%i_color_edges + sum(info%i_boundary_edges)
            traversal%current_stats%i_traversed_memory = traversal%current_stats%i_traversed_memory + sizeof(section%boundary_edges(RED)%elements(1)%t_crossed_edge_stream_data) * info%i_crossed_edges
            traversal%current_stats%i_traversed_memory = traversal%current_stats%i_traversed_memory + sizeof(section%boundary_edges(RED)%elements(1)%t_color_edge_stream_data) * (info%i_color_edges + sum(info%i_boundary_edges))
#   	endif

        traversal%stats = traversal%stats + traversal%current_stats
        section%stats = section%stats + traversal%current_stats
        call section%estimate_load()
    end subroutine

	!> initializes a ring buffer element by setting the cell pointers
	subroutine _OP0(init)(section, element)
		type(t_grid_section), intent(inout)					:: section
		type(t_traversal_element), intent(inout)			:: element
		type(t_cell_stream_data), pointer					:: p_cell_data

		p_cell_data => section%cells%next()
		element%cell%geometry => p_cell_data%fine_triangle
		element%cell%data_pers => p_cell_data%data_pers
	end subroutine

	!> sets the element transformation data
	subroutine _OP0(set_transform_data)(section, element)
		type(t_grid_section), intent(inout)					:: section
		type(t_traversal_element), intent(inout)			:: element						!< input element

#		if defined(_GT_NODES) .or. .not. defined(_GT_NO_COORDS)
			integer (kind = 1)								:: i_color_node_out, i_transfer_node, i_color_node_in
#		endif

#		if defined(_GT_EDGES)
			integer (kind = 1)								:: i, i_previous_edge, i_color_edge, i_next_edge
			type(t_cell_transform_data), pointer			:: p_plotter_data
#		endif

		assert_ne(element%cell%geometry%i_plotter_type, 0)
		assert_ge(element%cell%geometry%i_plotter_type, -8)
		assert_le(element%cell%geometry%i_plotter_type, 8)

		element%transform_data%plotter_data => ref_plotter_data(element%cell%geometry%i_plotter_type)
		element%transform_data%custom_data%scaling = element%cell%geometry%get_scaling()

#		if defined(_GT_NODES) .or. .not. defined(_GT_NO_COORDS)
			call element%cell%geometry%get_node_indices(i_color_node_out, i_transfer_node, i_color_node_in)

			element%nodes(i_color_node_out) = element%color_node_out
			element%nodes(i_transfer_node) = element%transfer_node
			element%nodes(i_color_node_in) = element%color_node_in

#           if .not. defined(_GT_NO_COORDS)
                element%transform_data%custom_data%offset => element%nodes(2)%ptr%position
#           endif
#		endif

#		if defined(_GT_EDGES)
			call element%cell%geometry%get_edge_indices(i_previous_edge, i_color_edge, i_next_edge)

			element%edges(i_previous_edge) = element%previous_edge
			element%edges(i_color_edge) = element%color_edge
			element%edges(i_next_edge) = element%next_edge

			element%edges(1)%ptr%depth = element%cell%geometry%i_depth
			element%edges(2)%ptr%depth = element%cell%geometry%i_depth - 1
			element%edges(3)%ptr%depth = element%cell%geometry%i_depth

			p_plotter_data => element%transform_data%plotter_data

			forall (i = 1 : 3)
                element%edges(i)%ptr%transform_data => p_plotter_data%edges(i)
            end forall
#		endif
	end subroutine

	subroutine _OP0(read) (traversal, thread, section, element)
		type(_GT), intent(inout)					            :: traversal
		type(t_grid_thread), intent(inout)					    :: thread
		type(t_grid_section), intent(inout)					    :: section

#       if defined(_ASSERT)
            type(t_traversal_element), target, intent(inout)    :: element
#       else
            type(t_traversal_element), intent(inout)	        :: element
#       endif

        !unroll all cases
		select case (element%cell%geometry%i_edge_types)
            case(INNER_OLD)
                call _OP3(read, o, o, n)(traversal, thread, section, element)
            case(INNER_NEW)
                call _OP3(read, o, n, n)(traversal, thread, section, element)
            case(INNER_OLD_BND)
                call _OP3(read, o, d, n)(traversal, thread, section, element)
            case(INNER_NEW_BND)
                call _OP3(read, o, b, n)(traversal, thread, section, element)
            case(FIRST_NEW)
                call _OP3(read, d, n, n)(traversal, thread, section, element)
            case(FIRST_OLD_BND)
                call _OP3(read, d, d, n)(traversal, thread, section, element)
            case(FIRST_NEW_BND)
                call _OP3(read, d, b, n)(traversal, thread, section, element)
            case(LAST_OLD)
                call _OP3(read, o, o, b)(traversal, thread, section, element)
            case(LAST_OLD_BND)
                call _OP3(read, o, d, b)(traversal, thread, section, element)
            case(LAST_NEW_BND)
                call _OP3(read, o, b, b)(traversal, thread, section, element)
            case(SINGLE_OLD_BND)
                call _OP3(read, d, d, b)(traversal, thread, section, element)
            case(SINGLE_NEW_BND)
                call _OP3(read, d, b, b)(traversal, thread, section, element)
			case default
				assert_eq(element%cell%geometry%i_edge_types, INNER_OLD)
        end select
	end subroutine

	subroutine _OP0(write) (traversal, thread, section, element)
		type(_GT), intent(inout)					            :: traversal
		type(t_grid_thread), intent(inout)					    :: thread
		type(t_grid_section), intent(inout)					    :: section

#       if defined(_ASSERT)
            type(t_traversal_element), target, intent(inout)    :: element
#       else
            type(t_traversal_element), intent(inout)	        :: element
#       endif

        !unroll all cases
		select case (element%cell%geometry%i_edge_types)
            case(INNER_OLD)
                call _OP3(write, o, o, n)(traversal, thread, section, element)
            case(INNER_NEW)
                call _OP3(write, o, n, n)(traversal, thread, section, element)
            case(INNER_OLD_BND)
                call _OP3(write, o, d, n)(traversal, thread, section, element)
            case(INNER_NEW_BND)
                call _OP3(write, o, b, n)(traversal, thread, section, element)
            case(FIRST_NEW)
                call _OP3(write, d, n, n)(traversal, thread, section, element)
            case(FIRST_OLD_BND)
                call _OP3(write, d, d, n)(traversal, thread, section, element)
            case(FIRST_NEW_BND)
                call _OP3(write, d, b, n)(traversal, thread, section, element)
            case(LAST_OLD)
                call _OP3(write, o, o, b)(traversal, thread, section, element)
            case(LAST_OLD_BND)
                call _OP3(write, o, d, b)(traversal, thread, section, element)
            case(LAST_NEW_BND)
                call _OP3(write, o, b, b)(traversal, thread, section, element)
            case(SINGLE_OLD_BND)
                call _OP3(write, d, d, b)(traversal, thread, section, element)
            case(SINGLE_NEW_BND)
                call _OP3(write, d, b, b)(traversal, thread, section, element)
			case default
				assert_eq(element%cell%geometry%i_edge_types, INNER_OLD)
        end select
	end subroutine
#endif

!implement these functions if all edge types are defined (which means there's a variant for each type combination)
#if (_GT_PREVIOUS_EDGE_TYPE .ne. _UNDEFINED) .and. (_GT_COLOR_EDGE_TYPE .ne. _UNDEFINED) .and. (_GT_NEXT_EDGE_TYPE .ne. _UNDEFINED)
	subroutine _OP3(read, _GT_PREVIOUS_EDGE_SUFFIX, _GT_COLOR_EDGE_SUFFIX, _GT_NEXT_EDGE_SUFFIX) (traversal, thread, section, element)
		type(_GT), intent(inout)					            :: traversal
		type(t_grid_thread), intent(inout)					    :: thread
		type(t_grid_section), intent(inout)					    :: section

#       if defined(_ASSERT)
            type(t_traversal_element), target, intent(inout)	:: element
#       else
            type(t_traversal_element), intent(inout)	        :: element
#       endif

		logical(kind = GRID_SL)								    :: l_color_edge_color
		real (kind = GRID_SI), dimension(2, 3), parameter	    :: node_offset = reshape([0.0, 1.0, -1.0, 1.0, -1.0, 0.0], [2, 3])

#		if defined(_GT_SKELETON_OP) .or. defined(_GT_CELL_UPDATE_OP)
            type(num_cell_update)                               :: color_edge_local_update, next_edge_local_update
#       endif

	 	l_color_edge_color = element%cell%geometry%l_color_edge_color

		!read previous crossed edge
#		if (_GT_PREVIOUS_EDGE_TYPE == _OLD)
			!read transfer node from stack
			_GT_N(element%transfer_node%ptr => thread%nodes_stack(.not. l_color_edge_color)%current())
#		elif (_GT_PREVIOUS_EDGE_TYPE == _OLD_BND)
			!access previous edge on the process edge stream
			_GT_E(element%previous_edge%ptr => section%boundary_edges(RED)%next())

			!read transfer nodes from the process node stream and push them on the stacks
			_GT_N(element%color_node_out%ptr => thread%nodes_stack(l_color_edge_color)%push())
			_GT_N(call section%boundary_nodes(l_color_edge_color)%read(element%color_node_out%ptr))

			_GT_N(element%transfer_node%ptr => thread%nodes_stack(.not. l_color_edge_color)%push())
			_GT_N(call section%boundary_nodes(.not. l_color_edge_color)%read(element%transfer_node%ptr))
#		endif

		!read color edge and color node
#		if (_GT_COLOR_EDGE_TYPE == _OLD)
			_GT_E(element%color_edge%ptr => thread%edges_stack(l_color_edge_color)%pop())

			_GT_N(element%color_node_out%ptr => thread%nodes_stack(l_color_edge_color)%pop())
			_GT_N(element%color_node_in%ptr => thread%nodes_stack(l_color_edge_color)%current())
#		elif (_GT_COLOR_EDGE_TYPE == _NEW)
			_GT_E(element%color_edge%ptr => thread%edges_stack(l_color_edge_color)%push())
			_GT_ER(call section%color_edges_in%read(element%color_edge%ptr%t_color_edge_stream_data))

			_GT_N(element%color_node_out%ptr => thread%nodes_stack(l_color_edge_color)%current())
			_GT_N(element%color_node_in%ptr => thread%nodes_stack(l_color_edge_color)%push())
			_GT_NR(call section%nodes_in%read(element%color_node_in%ptr%t_node_stream_data))
#		elif (_GT_COLOR_EDGE_TYPE == _OLD_BND) .or. (_GT_COLOR_EDGE_TYPE == _NEW_BND)
			_GT_E(element%color_edge%ptr => section%boundary_edges(l_color_edge_color)%next())

			_GT_N(element%color_node_out%ptr => section%boundary_nodes(l_color_edge_color)%current())
			_GT_N(call thread%nodes_stack(l_color_edge_color)%pop(element%color_node_out%ptr))

			_GT_N(element%color_node_in%ptr => thread%nodes_stack(l_color_edge_color)%push())
			_GT_N(call section%boundary_nodes(l_color_edge_color)%read(element%color_node_in%ptr))
#		endif

		!read next crossed edge
#		if (_GT_NEXT_EDGE_TYPE == _NEW)
			_GT_ER(call section%crossed_edges_in%read(element%next_edge%ptr%t_crossed_edge_stream_data))
#		elif (_GT_NEXT_EDGE_TYPE == _NEW_BND)
			_GT_E(element%next_edge%ptr => section%boundary_edges(RED)%next())

			!pop transfer nodes from the stacks and write them to the process node streams
			_GT_N(element%color_node_in%ptr => section%boundary_nodes(l_color_edge_color)%current())
			_GT_N(call thread%nodes_stack(l_color_edge_color)%pop(element%color_node_in%ptr))

			_GT_N(element%transfer_node%ptr => section%boundary_nodes(.not. l_color_edge_color)%current())
			_GT_N(call thread%nodes_stack(.not. l_color_edge_color)%pop(element%transfer_node%ptr))
#		endif

		!set element tranformation data (must be defined)
		call _OP0(set_transform_data) (section, element)

#		if .not. defined(_GT_NO_COORDS)
#			if .not. defined(_STORE_NODE_COORDS) .and. (_GT_COLOR_EDGE_TYPE == _NEW)
				element%color_node_in%ptr%position = element%color_node_out%ptr%position + element%transform_data%custom_data%scaling * matmul(element%transform_data%plotter_data%jacobian, node_offset(:, element%cell%geometry%i_turtle_type))
#			endif
#		endif

		!First Touch Hooks

#		if defined(_GT_CELL_FIRST_TOUCH_OP)
			call _GT_CELL_FIRST_TOUCH_OP(traversal, section, element%cell)
#		endif

#		if (_GT_COLOR_EDGE_TYPE == _NEW)
#			if defined(_GT_INNER_EDGE_FIRST_TOUCH_OP)
				_GT_E(call _GT_INNER_EDGE_FIRST_TOUCH_OP(traversal, section, element%color_edge%ptr))
#			endif

#			if defined(_GT_INNER_NODE_FIRST_TOUCH_OP)
				_GT_N(call _GT_INNER_NODE_FIRST_TOUCH_OP(traversal, section, element%color_node_in%ptr))
#			endif
#		endif

#		if (_GT_NEXT_EDGE_TYPE == _NEW)
#           if defined(_GT_EDGES)
                assert(associated(element%next_edge%ptr, element%next_edge_data))
#           endif

#			if defined(_GT_INNER_EDGE_FIRST_TOUCH_OP)
				_GT_E(call _GT_INNER_EDGE_FIRST_TOUCH_OP(traversal, section, element%next_edge%ptr))
#			endif
#		endif

		!Call skeleton operator on inner new edges
#		if defined(_GT_SKELETON_OP)
#			if (_GT_COLOR_EDGE_TYPE == _NEW)
				call _GT_SKELETON_OP(traversal, section, element%color_edge%ptr, _GT_CELL_TO_EDGE_OP(element%t_element_base, element%color_edge%ptr), element%color_edge%ptr%rep, color_edge_local_update, element%color_edge%ptr%update)
#			endif

#			if (_GT_NEXT_EDGE_TYPE == _NEW)
				call _GT_SKELETON_OP(traversal, section, element%next_edge%ptr, _GT_CELL_TO_EDGE_OP(element%t_element_base, element%next_edge%ptr), _GT_CELL_TO_EDGE_OP(element%next%t_element_base, element%next_edge%ptr), next_edge_local_update, element%next_edge%ptr%update)
#			elif (_GT_NEXT_EDGE_TYPE == _NEW_BND)
				next_edge_local_update = element%next_edge%ptr%update
#			endif
#		endif

		!Call cell update operator
#		if defined(_GT_CELL_UPDATE_OP)
#			if (_GT_COLOR_EDGE_TYPE == _OLD) .or. (_GT_COLOR_EDGE_TYPE == _OLD_BND) .or. (_GT_COLOR_EDGE_TYPE == _NEW_BND)
                select case (element%cell%geometry%i_turtle_type)
                    case (K)
                        call _GT_CELL_UPDATE_OP(traversal, section, element%t_element_base, element%color_edge%ptr%update, next_edge_local_update, element%previous_edge%ptr%update)
                    case (V)
                        call _GT_CELL_UPDATE_OP(traversal, section, element%t_element_base, next_edge_local_update, element%color_edge%ptr%update, element%previous_edge%ptr%update)
                    case (H)
                        call _GT_CELL_UPDATE_OP(traversal, section, element%t_element_base, next_edge_local_update, element%previous_edge%ptr%update, element%color_edge%ptr%update)
                end select
#			elif (_GT_COLOR_EDGE_TYPE == _NEW)
                select case (element%cell%geometry%i_turtle_type)
                    case (K)
                        call _GT_CELL_UPDATE_OP(traversal, section, element%t_element_base, color_edge_local_update, next_edge_local_update, element%previous_edge%ptr%update)
                    case (V)
                        call _GT_CELL_UPDATE_OP(traversal, section, element%t_element_base, next_edge_local_update, color_edge_local_update, element%previous_edge%ptr%update)
                    case (H)
                        call _GT_CELL_UPDATE_OP(traversal, section, element%t_element_base, next_edge_local_update, element%previous_edge%ptr%update, color_edge_local_update)
                end select
#			endif
#		endif
	end subroutine
#endif

!implement these functions if all edge types are defined (which means there's a variant for each type combination)
#if (_GT_PREVIOUS_EDGE_TYPE .ne. _UNDEFINED) .and. (_GT_COLOR_EDGE_TYPE .ne. _UNDEFINED) .and. (_GT_NEXT_EDGE_TYPE .ne. _UNDEFINED)
	subroutine _OP3(write, _GT_PREVIOUS_EDGE_SUFFIX, _GT_COLOR_EDGE_SUFFIX, _GT_NEXT_EDGE_SUFFIX) (traversal, thread, section, element)
		type(_GT), intent(inout)					            :: traversal
		type(t_grid_thread), intent(inout)					    :: thread
		type(t_grid_section), intent(inout)				        :: section

#       if defined(_ASSERT)
            type(t_traversal_element), target, intent(inout)	:: element
#       else
            type(t_traversal_element), intent(inout)	        :: element
#       endif

		logical(kind = GRID_SL)							        :: l_color_edge_color

		l_color_edge_color = element%cell%geometry%l_color_edge_color

#		if defined(_GT_CELL_TO_EDGE_OP)
			!Copy cell representations to edges

#			if (_GT_PREVIOUS_EDGE_TYPE == _OLD_BND)
                element%previous_edge%ptr%rep = _GT_CELL_TO_EDGE_OP(element%t_element_base, element%previous_edge%ptr)
#			endif

#			if (_GT_COLOR_EDGE_TYPE == _NEW) .or. (_GT_COLOR_EDGE_TYPE == _OLD_BND) .or. (_GT_COLOR_EDGE_TYPE == _NEW_BND)
				element%color_edge%ptr%rep = _GT_CELL_TO_EDGE_OP(element%t_element_base, element%color_edge%ptr)
#			endif

#			if (_GT_NEXT_EDGE_TYPE == _NEW_BND)
				element%next_edge%ptr%rep = _GT_CELL_TO_EDGE_OP(element%t_element_base, element%next_edge%ptr)
#			endif
#		endif

		!Last Touch Hooks

#		if defined(_GT_CELL_LAST_TOUCH_OP)
			call _GT_CELL_LAST_TOUCH_OP(traversal, section, element%cell)
#		endif

#		if defined(_GT_CELL_REDUCE_OP)
			call _GT_CELL_REDUCE_OP(traversal, section, element%cell)
#		endif

#		if (_GT_PREVIOUS_EDGE_TYPE == _OLD)
#           if defined(_GT_EDGES)
                assert(associated(element%previous_edge%ptr, element%previous%next_edge_data))
#           endif

#			if defined(_GT_INNER_EDGE_LAST_TOUCH_OP)
				_GT_E(call _GT_INNER_EDGE_LAST_TOUCH_OP(traversal, section, element%previous_edge%ptr))
#			endif

#			if defined(_GT_EDGE_REDUCE_OP)
                _GT_E(call _GT_INNER_EDGE_REDUCE_OP(traversal, section, element%previous_edge%ptr))
#			endif
#		endif

#		if (_GT_COLOR_EDGE_TYPE == _OLD)
#			if defined(_GT_INNER_EDGE_LAST_TOUCH_OP)
				_GT_E(call _GT_INNER_EDGE_LAST_TOUCH_OP(traversal, section, element%color_edge%ptr))
#			endif

#			if defined(_GT_EDGE_REDUCE_OP)
                _GT_E(call _GT_EDGE_REDUCE_OP(traversal, section, element%color_edge%ptr))
#			endif

#			if defined(_GT_INNER_NODE_LAST_TOUCH_OP)
				_GT_N(call _GT_INNER_NODE_LAST_TOUCH_OP(traversal, section, element%color_node_out%ptr))
#			endif

#			if defined(_GT_NODE_REDUCE_OP)
                _GT_N(call _GT_INNER_NODE_REDUCE_OP(traversal, section, element%color_node_out%ptr))
#           endif
#		endif

#		if (_GT_PREVIOUS_EDGE_TYPE == _OLD)
            _GT_EW(call section%crossed_edges_out%write(element%previous_edge%ptr%t_crossed_edge_stream_data))
#		elif (_GT_PREVIOUS_EDGE_TYPE == _OLD_BND)
            _GT_E(element%previous_edge%ptr => element%previous%next_edge_data)
#       endif

		!write color edge & transfer node to stream if the edge is old or on the boundary, otherwise on the stack
#		if (_GT_COLOR_EDGE_TYPE == _OLD)
			_GT_EW(call section%color_edges_out%write(element%color_edge%ptr%t_color_edge_stream_data))
			_GT_NW(call section%nodes_out%write(element%color_node_out%ptr%t_node_stream_data))
#		endif

#		if (_GT_NEXT_EDGE_TYPE == _NEW_BND)
            _GT_E(element%next_edge%ptr => element%next%previous_edge%ptr)
#		endif

#       if defined(_GT_OUTPUT)
#		    if (_GT_PREVIOUS_EDGE_TYPE == _OLD .and. _GT_NEXT_EDGE_TYPE == _NEW)
                call element%cell%geometry%reverse_inner()
#		    else
                call element%cell%geometry%reverse()
#		    endif
#       endif
	end subroutine
#endif

#endif
#endif
#endif

#undef _GT_N
#undef _GT_NT
#undef _GT_NNT
#undef _GT_NR
#undef _GT_NW
#undef _GT_E
#undef _GT_ET
#undef _GT_ENT
#undef _GT_ER
#undef _GT_EW
