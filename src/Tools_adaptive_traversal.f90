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

#if !defined(_GT_PREVIOUS_EDGE_TYPE)
#	define _GT_PREVIOUS_EDGE_TYPE		_UNDEFINED
#	include "Tools_adaptive_traversal.f90"
#	undef _GT_PREVIOUS_EDGE_TYPE

#	define _GT_PREVIOUS_EDGE_TYPE		_OLD
#	define _GT_PREVIOUS_EDGE_SUFFIX		o
#	include "Tools_adaptive_traversal.f90"
#	undef _GT_PREVIOUS_EDGE_SUFFIX
#	undef _GT_PREVIOUS_EDGE_TYPE

#	define _GT_PREVIOUS_EDGE_TYPE		_OLD_BND
#	define _GT_PREVIOUS_EDGE_SUFFIX		d
#	include "Tools_adaptive_traversal.f90"
#	undef _GT_PREVIOUS_EDGE_SUFFIX
#	undef _GT_PREVIOUS_EDGE_TYPE
#else

#if !defined(_GT_COLOR_EDGE_TYPE)
#	define _GT_COLOR_EDGE_TYPE			_UNDEFINED
#	include "Tools_adaptive_traversal.f90"
#	undef _GT_COLOR_EDGE_TYPE

#	define _GT_COLOR_EDGE_TYPE			_OLD
#	define _GT_COLOR_EDGE_SUFFIX		o
#	include "Tools_adaptive_traversal.f90"
#	undef _GT_COLOR_EDGE_SUFFIX
#	undef _GT_COLOR_EDGE_TYPE

#	define _GT_COLOR_EDGE_TYPE			_NEW
#	define _GT_COLOR_EDGE_SUFFIX		n
#	include "Tools_adaptive_traversal.f90"
#	undef _GT_COLOR_EDGE_SUFFIX
#	undef _GT_COLOR_EDGE_TYPE

#	define _GT_COLOR_EDGE_TYPE			_OLD_BND
#	define _GT_COLOR_EDGE_SUFFIX		d
#	include "Tools_adaptive_traversal.f90"
#	undef _GT_COLOR_EDGE_SUFFIX
#	undef _GT_COLOR_EDGE_TYPE

#	define _GT_COLOR_EDGE_TYPE			_NEW_BND
#	define _GT_COLOR_EDGE_SUFFIX		b
#	include "Tools_adaptive_traversal.f90"
#	undef _GT_COLOR_EDGE_SUFFIX
#	undef _GT_COLOR_EDGE_TYPE

#else

#if !defined(_GT_NEXT_EDGE_TYPE)
#	define _GT_NEXT_EDGE_TYPE			_UNDEFINED
#	include "Tools_adaptive_traversal.f90"
#	undef _GT_NEXT_EDGE_TYPE

#	define _GT_NEXT_EDGE_TYPE			_NEW
#	define _GT_NEXT_EDGE_SUFFIX			n
#	include "Tools_adaptive_traversal.f90"
#	undef _GT_NEXT_EDGE_SUFFIX
#	undef _GT_NEXT_EDGE_TYPE

#	define _GT_NEXT_EDGE_TYPE			_NEW_BND
#	define _GT_NEXT_EDGE_SUFFIX			b
#	include "Tools_adaptive_traversal.f90"
#	undef _GT_NEXT_EDGE_SUFFIX
#	undef _GT_NEXT_EDGE_TYPE
#else

#if !defined(_GT_PREFIX)
#	define _OP0(x)						x
#	define _OP1(x, a)					_SUFFIX1(x,a)
#	define _OP2(x, a, b)				_SUFFIX2(x,a,b)
#	define _OP3(x, a, b, c)				_SUFFIX3(x,a,b,c)
#else
#	define _OP0(x)						_SUFFIX0E(x,_GT_PREFIX)
#	define _OP1(x, a)					_SUFFIX1E(x,_GT_PREFIX,a)
#	define _OP2(x, a, b)				_SUFFIX2E(x,_GT_PREFIX,a,b)
#	define _OP3(x, a, b, c)				_SUFFIX3E(x,_GT_PREFIX,a,b,c)
#endif

#define _SUFFIX1(x, a)					_conc3(x,_,a)
#define _SUFFIX2(x, a, b)				_conc4(x,_,a,b)
#define _SUFFIX3(x, a, b, c)			_conc5(x,_,a,b,c)

#define _SUFFIX0E(x, y)					_conc3(x,_,y)
#define _SUFFIX1E(x, y, a)				_conc5(x,_,y,_,a)
#define _SUFFIX2E(x, y, a, b)			_conc6(x,_,y,_,a,b)
#define _SUFFIX3E(x, y, a, b, c)		_conc7(x,_,y,_,a,b,c)

#define _UNDEFINED						-1
#define _OLD							0
#define _NEW							1
#define _OLD_BND					    2
#define _NEW_BND						3

! Short forms for node macros (enable/disable code execution based on preprocessor directives)

#if defined(_GT_NODES) || !defined(_GT_NO_COORDS)
#	define _GT_N(x)			x
#else
#	define _GT_N(x)
#endif

#if defined(_GT_NODES) && defined(_GT_INPUT)
#	define _GT_NR(x)		x
#else
#	define _GT_NR(x)
#endif

#if defined(_GT_NODES) && defined(_GT_OUTPUT)
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

#if defined(_GT_EDGES) && defined(_GT_INPUT)
#	define _GT_ER(x)		x
#else
#	define _GT_ER(x)
#endif

#if defined(_GT_EDGES) && defined(_GT_OUTPUT)
#	define _GT_EW(x)		x
#else
#	define _GT_EW(x)
#endif

!For each skeleton operator there must be a boundary skeleton operator defined
#if defined(_GT_SKELETON_OP) && !(defined(_GT_BND_SKELETON_OP) && defined(_GT_CELL_TO_EDGE_OP) && defined(_GT_CELL_UPDATE_OP))
#	error "If a skeleton operator exists, a boundary skeleton operator, a cell-to-edge operator and a cell update operator must also be defined"
#endif

!implement this function if no edge types are set (which means only once)
#if (_GT_PREVIOUS_EDGE_TYPE == _UNDEFINED) && (_GT_COLOR_EDGE_TYPE == _UNDEFINED) && (_GT_NEXT_EDGE_TYPE == _UNDEFINED)
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

#       if defined(_GT_EDGES) && defined(_GT_SKELETON_OP)
	        integer (kind = GRID_SI)							    :: i_color
            integer (kind = GRID_SI)							    :: i_comm
            type(t_comm_interface), pointer                         :: comm

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

		integer (kind = GRID_SI)							        :: i_color, i_node, i_edge

#		if defined(_GT_PRE_TRAVERSAL_OP)
			call _GT_PRE_TRAVERSAL_OP(traversal, section)
#		endif

		do i_color = RED, GREEN
#			if defined(_GT_EDGES)
                _log_write(5, '(5X, A)'), "boundary edges:"

                do i_edge = 1, size(section%boundary_edges(i_color)%elements)
#			        if defined(_GT_EDGE_FIRST_TOUCH_OP)
                        if (section%boundary_nodes(i_color)%elements(i_node)%owned_locally) then
                            call _GT_EDGE_FIRST_TOUCH_OP(traversal, section, section%boundary_edges(i_color)%elements)
                        endif
#			        endif
                end do
#			endif

#			if defined(_GT_NODES)
                _log_write(5, '(5X, A)'), "boundary nodes:"

                do i_node = 1, size(section%boundary_nodes(i_color)%elements)
#    			    if defined(_GT_NODE_FIRST_TOUCH_OP)
                        if (section%boundary_nodes(i_color)%elements(i_node)%owned_locally) then
                            call _GT_NODE_FIRST_TOUCH_OP(traversal, section, section%boundary_nodes(i_color)%elements(i_node))
                        endif
#                   endif
                end do
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

                do i_edge = 1, size(section%boundary_edges(i_color)%elements)
#			        if defined(_GT_EDGE_LAST_TOUCH_OP)
                        if (section%boundary_nodes(i_color)%elements(i_node)%owned_locally) then
                            call _GT_EDGE_LAST_TOUCH_OP(traversal, section, section%boundary_edges(i_color)%elements)
                        endif
#			        endif

#			        if defined(_GT_EDGE_REDUCE_OP)
                        if (section%boundary_edges(i_color)%elements(i_edge)%owned_globally) then
                            call _GT_EDGE_REDUCE_OP(traversal, section, section%boundary_edges(i_color)%elements(i_edge))
                        endif
#			        endif
                end do
#			endif

#			if defined(_GT_NODES)
                _log_write(5, '(5X, A)'), "boundary nodes:"

                do i_node = 1, size(section%boundary_nodes(i_color)%elements)
#    			    if defined(_GT_NODE_LAST_TOUCH_OP)
                        if (section%boundary_nodes(i_color)%elements(i_node)%owned_locally) then
                            call _GT_NODE_LAST_TOUCH_OP(traversal, section, section%boundary_nodes(i_color)%elements(i_node))
                        endif
#                   endif

#    			    if defined(_GT_NODE_REDUCE_OP)
                        if (section%boundary_nodes(i_color)%elements(i_node)%owned_globally) then
                            call _GT_NODE_REDUCE_OP(traversal, section, section%boundary_nodes(i_color)%elements(i_node))
                        endif
#                   endif
                end do
#			endif
		end do

#		if defined(_GT_POST_TRAVERSAL_OP)
			call _GT_POST_TRAVERSAL_OP(traversal, section)
#		endif
	end subroutine

    subroutine _OP0(set_stats_counters)(stats, section)
        type(t_statistics), intent(inout)		    :: stats
        type(t_grid_section), intent(in)			:: section

        type(t_section_info)               	        :: info

        info = section%get_info()

        stats%i_traversals = 1
        stats%i_traversed_cells = info%i_cells
       	stats%i_traversed_memory = sizeof(section%cells%elements(1)) * info%i_cells

#   	if defined(_GT_NODES)
            stats%i_traversed_nodes = info%i_nodes + sum(info%i_boundary_nodes)
            stats%i_traversed_memory = stats%i_traversed_memory + sizeof(section%boundary_nodes(RED)%elements(1)%t_node_stream_data) * (info%i_nodes + sum(info%i_boundary_nodes))
#   	endif

#   	if defined(_GT_EDGES)
            stats%i_traversed_edges = info%i_crossed_edges + info%i_color_edges + sum(info%i_boundary_edges)
            stats%i_traversed_memory = stats%i_traversed_memory + sizeof(section%boundary_edges(RED)%elements(1)%t_crossed_edge_stream_data) * info%i_crossed_edges
            stats%i_traversed_memory = stats%i_traversed_memory + sizeof(section%boundary_edges(RED)%elements(1)%t_color_edge_stream_data) * (info%i_color_edges + sum(info%i_boundary_edges))
#   	endif
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

#		if defined(_GT_NODES) || !defined(_GT_NO_COORDS)
			integer (kind = BYTE)								:: i_color_node_out, i_transfer_node, i_color_node_in
#		endif

#		if defined(_GT_EDGES)
			integer (kind = BYTE)								:: i, i_previous_edge, i_color_edge, i_next_edge
			type(t_cell_transform_data), pointer			:: p_plotter_data
#		endif

		assert_ne(element%cell%geometry%i_plotter_type, 0)
		assert_ge(element%cell%geometry%i_plotter_type, -8)
		assert_le(element%cell%geometry%i_plotter_type, 8)

		element%transform_data%plotter_data => ref_plotter_data(element%cell%geometry%i_plotter_type)
		element%transform_data%custom_data%scaling = element%cell%geometry%get_scaling()

#		if defined(_GT_NODES) || !defined(_GT_NO_COORDS)
			call element%cell%geometry%get_node_indices(i_color_node_out, i_transfer_node, i_color_node_in)

			element%nodes(i_color_node_out) = element%color_node_out
			element%nodes(i_transfer_node) = element%transfer_node
			element%nodes(i_color_node_in) = element%color_node_in

#           if !defined(_GT_NO_COORDS)
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

			do i = 1, 3
                element%edges(i)%ptr%transform_data => p_plotter_data%edges(i)
            end do
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
		select case (element%cell%geometry%i_entity_types)
            case(INNER_OLD)
                call _OP3(read,o,o,n)(traversal, thread, section, element)
            case(INNER_NEW)
                call _OP3(read,o,n,n)(traversal, thread, section, element)
            case(INNER_OLD_BND)
                call _OP3(read,o,d,n)(traversal, thread, section, element)
            case(INNER_NEW_BND)
                call _OP3(read,o,b,n)(traversal, thread, section, element)
            case(FIRST_NEW)
                call _OP3(read,d,n,n)(traversal, thread, section, element)
            case(FIRST_OLD_BND)
                call _OP3(read,d,d,n)(traversal, thread, section, element)
            case(FIRST_NEW_BND)
                call _OP3(read,d,b,n)(traversal, thread, section, element)
            case(LAST_OLD)
                call _OP3(read,o,o,b)(traversal, thread, section, element)
            case(LAST_OLD_BND)
                call _OP3(read,o,d,b)(traversal, thread, section, element)
            case(LAST_NEW_BND)
                call _OP3(read,o,b,b)(traversal, thread, section, element)
            case(SINGLE_OLD_BND)
                call _OP3(read,d,d,b)(traversal, thread, section, element)
            case(SINGLE_NEW_BND)
                call _OP3(read,d,b,b)(traversal, thread, section, element)
			case default
				assert_eq(element%cell%geometry%i_entity_types, INNER_OLD)
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
		select case (element%cell%geometry%i_entity_types)
            case(INNER_OLD)
                call _OP3(write,o,o,n)(traversal, thread, section, element)
            case(INNER_NEW)
                call _OP3(write,o,n,n)(traversal, thread, section, element)
            case(INNER_OLD_BND)
                call _OP3(write,o,d,n)(traversal, thread, section, element)
            case(INNER_NEW_BND)
                call _OP3(write,o,b,n)(traversal, thread, section, element)
            case(FIRST_NEW)
                call _OP3(write,d,n,n)(traversal, thread, section, element)
            case(FIRST_OLD_BND)
                call _OP3(write,d,d,n)(traversal, thread, section, element)
            case(FIRST_NEW_BND)
                call _OP3(write,d,b,n)(traversal, thread, section, element)
            case(LAST_OLD)
                call _OP3(write,o,o,b)(traversal, thread, section, element)
            case(LAST_OLD_BND)
                call _OP3(write,o,d,b)(traversal, thread, section, element)
            case(LAST_NEW_BND)
                call _OP3(write,o,b,b)(traversal, thread, section, element)
            case(SINGLE_OLD_BND)
                call _OP3(write,d,d,b)(traversal, thread, section, element)
            case(SINGLE_NEW_BND)
                call _OP3(write,d,b,b)(traversal, thread, section, element)
			case default
				assert_eq(element%cell%geometry%i_entity_types, INNER_OLD)
        end select
	end subroutine
#endif

!implement these functions if all edge types are defined (which means there's a variant for each type combination)
#if (_GT_PREVIOUS_EDGE_TYPE != _UNDEFINED) && (_GT_COLOR_EDGE_TYPE != _UNDEFINED) && (_GT_NEXT_EDGE_TYPE != _UNDEFINED)
	subroutine _OP3(read,_GT_PREVIOUS_EDGE_SUFFIX,_GT_COLOR_EDGE_SUFFIX,_GT_NEXT_EDGE_SUFFIX) (traversal, thread, section, element)
		type(_GT), intent(inout)					            :: traversal
		type(t_grid_thread), intent(inout)					    :: thread
		type(t_grid_section), intent(inout)					    :: section

#       if defined(_ASSERT)
            type(t_traversal_element), target, intent(inout)	:: element
#       else
            type(t_traversal_element), intent(inout)	        :: element
#       endif

		integer(kind = BYTE)								    :: i_color_edge_color
		real (kind = GRID_SI), dimension(2, 3), parameter	    :: node_offset = reshape([0.0, 1.0, -1.0, 1.0, -1.0, 0.0], [2, 3])

#       if !defined(_GT_NO_COORDS) && !defined(_GT_PASS_COORDS) && !defined(_STORE_NODE_COORDS) && (_GT_COLOR_EDGE_TYPE == _NEW)
            real (kind = GRID_SI)                               :: new_position(2)
#       endif

#		if defined(_GT_SKELETON_OP) || defined(_GT_CELL_UPDATE_OP)
            type(num_cell_update)                               :: color_edge_local_update, next_edge_local_update
#       endif

	 	i_color_edge_color = element%cell%geometry%i_color_edge_color

		!read previous crossed edge
#		if (_GT_PREVIOUS_EDGE_TYPE == _OLD)
            if (thread%nodes_stack(RED + GREEN - i_color_edge_color)%is_empty()) then
			    _GT_N(element%transfer_node%ptr => section%boundary_nodes(RED + GREEN - i_color_edge_color)%current())
            else
			    _GT_N(element%transfer_node%ptr => thread%nodes_stack(RED + GREEN - i_color_edge_color)%current())
			end if
#		elif (_GT_PREVIOUS_EDGE_TYPE == _OLD_BND)
			!access previous edge on the process edge stream
			_GT_E(element%previous_edge%ptr => section%boundary_edges(RED)%next())

            _GT_N(element%color_node_out%ptr => section%boundary_nodes(i_color_edge_color)%next())
            _GT_N(element%transfer_node%ptr => section%boundary_nodes(RED + GREEN - i_color_edge_color)%next())
#		endif

		!read color edge and color node
#		if (_GT_COLOR_EDGE_TYPE == _OLD)
			_GT_E(element%color_edge%ptr => thread%edges_stack(i_color_edge_color)%pop())
			_GT_N(element%color_node_out%ptr => thread%nodes_stack(i_color_edge_color)%pop())

            if (thread%nodes_stack(i_color_edge_color)%is_empty()) then
			    _GT_N(element%color_node_in%ptr => section%boundary_nodes(i_color_edge_color)%current())
            else
			    _GT_N(element%color_node_in%ptr => thread%nodes_stack(i_color_edge_color)%current())
            end if
#		elif (_GT_COLOR_EDGE_TYPE == _NEW)
            if (thread%nodes_stack(i_color_edge_color)%is_empty()) then
			    _GT_N(element%color_node_out%ptr => section%boundary_nodes(i_color_edge_color)%current())
            else
			    _GT_N(element%color_node_out%ptr => thread%nodes_stack(i_color_edge_color)%current())
			end if

			_GT_N(element%color_node_in%ptr => thread%nodes_stack(i_color_edge_color)%push())
            _GT_NR(call section%nodes_in%read(element%color_node_in%ptr%t_node_stream_data))

			_GT_E(element%color_edge%ptr => thread%edges_stack(i_color_edge_color)%push())
            _GT_ER(call section%color_edges_in%read(element%color_edge%ptr%t_color_edge_stream_data))
#		elif (_GT_COLOR_EDGE_TYPE == _OLD_BND) || (_GT_COLOR_EDGE_TYPE == _NEW_BND)
			_GT_E(element%color_edge%ptr => section%boundary_edges(i_color_edge_color)%next())

            _GT_N(element%color_node_out%ptr => section%boundary_nodes(i_color_edge_color)%current())
            _GT_N(element%color_node_in%ptr => section%boundary_nodes(i_color_edge_color)%next())
#		endif

		!read next crossed edge
#		if (_GT_NEXT_EDGE_TYPE == _NEW)
            _GT_ER(call section%crossed_edges_in%read(element%next_edge%ptr%t_crossed_edge_stream_data))
#		elif (_GT_NEXT_EDGE_TYPE == _NEW_BND)
			_GT_E(element%next_edge%ptr => section%boundary_edges(RED)%next())

            _GT_N(element%color_node_in%ptr => section%boundary_nodes(i_color_edge_color)%current())
            _GT_N(element%transfer_node%ptr => section%boundary_nodes(RED + GREEN - i_color_edge_color)%current())
#		endif

		!set element tranformation data (must be defined)
		call _OP0(set_transform_data) (section, element)

#		if !defined(_GT_NO_COORDS)
#           if defined(_GT_PASS_COORDS)
#               if (_GT_COLOR_EDGE_TYPE == _OLD) && (_GT_PREVIOUS_EDGE_TYPE == _OLD)

#               elif (_GT_PREVIOUS_EDGE_TYPE == _OLD)
                    element%color_node_in%ptr%position = element%color_node_out%ptr%position + element%transform_data%custom_data%scaling * matmul(element%transform_data%plotter_data%jacobian, node_offset(:, element%cell%geometry%i_turtle_type))
#               else
                    element%nodes(1)%ptr%position(:) = element%coords(:, 1)
                    element%nodes(2)%ptr%position(:) = element%coords(:, 2)
                    element%nodes(3)%ptr%position(:) = element%coords(:, 3)
#               endif
#		    elif !defined(_STORE_NODE_COORDS)
#               if (_GT_COLOR_EDGE_TYPE == _NEW)
                    new_position = element%color_node_out%ptr%position + element%transform_data%custom_data%scaling * matmul(element%transform_data%plotter_data%jacobian, node_offset(:, element%cell%geometry%i_turtle_type))
                    element%color_node_in%ptr%position = new_position
#               endif
#           endif
#		endif

        assert_veq(element%nodes(1)%ptr%position, element%nodes(1)%ptr%position)
        assert_veq(element%nodes(2)%ptr%position, element%nodes(2)%ptr%position)
        assert_veq(element%nodes(3)%ptr%position, element%nodes(3)%ptr%position)

        !assert(.not. any(element%nodes(1)%ptr%position > 0.0_SR .and. element%nodes(1)%ptr%position < 1.0e-100_SR))
        !assert(.not. any(element%nodes(2)%ptr%position > 0.0_SR .and. element%nodes(2)%ptr%position < 1.0e-100_SR))
        !assert(.not. any(element%nodes(3)%ptr%position > 0.0_SR .and. element%nodes(3)%ptr%position < 1.0e-100_SR))

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
#			if (_GT_COLOR_EDGE_TYPE == _OLD) || (_GT_COLOR_EDGE_TYPE == _OLD_BND) || (_GT_COLOR_EDGE_TYPE == _NEW_BND)
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
#if (_GT_PREVIOUS_EDGE_TYPE != _UNDEFINED) && (_GT_COLOR_EDGE_TYPE != _UNDEFINED) && (_GT_NEXT_EDGE_TYPE != _UNDEFINED)
	subroutine _OP3(write,_GT_PREVIOUS_EDGE_SUFFIX,_GT_COLOR_EDGE_SUFFIX,_GT_NEXT_EDGE_SUFFIX) (traversal, thread, section, element)
		type(_GT), intent(inout)					            :: traversal
		type(t_grid_thread), intent(inout)					    :: thread
		type(t_grid_section), intent(inout)				        :: section

#       if defined(_ASSERT)
            type(t_traversal_element), target, intent(inout)	:: element
#       else
            type(t_traversal_element), intent(inout)	        :: element
#       endif

		integer (kind = BYTE)							        :: i_color_edge_color

		i_color_edge_color = element%cell%geometry%i_color_edge_color

#		if defined(_GT_CELL_TO_EDGE_OP)
			!Copy cell representations to edges

#			if (_GT_PREVIOUS_EDGE_TYPE == _OLD_BND)
                element%previous_edge%ptr%rep = _GT_CELL_TO_EDGE_OP(element%t_element_base, element%previous_edge%ptr)
#			endif

#			if (_GT_COLOR_EDGE_TYPE == _NEW) || (_GT_COLOR_EDGE_TYPE == _OLD_BND) || (_GT_COLOR_EDGE_TYPE == _NEW_BND)
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
#		    if (_GT_PREVIOUS_EDGE_TYPE == _OLD && _GT_NEXT_EDGE_TYPE == _NEW)
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
