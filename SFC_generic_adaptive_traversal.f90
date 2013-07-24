! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


!> Generic traversal template
!> Warning: this a template module body and requires preprocessor commands before inclusion.
!>
!> The resulting method is defined as _GT_NAME
!> @author Oliver Meister

!multiple levels of indirection are necessary to properly resolve the names
#define _CONC2(X, Y)			X ## _ ## Y
#define _PREFIX(P, X)			_CONC2(P, X)
#define _GT_(X)					_PREFIX(_GT_NAME, X)

#define _GT						_GT_NAME

!if no dedicated inner operators exists, use the default operators
#if defined(_GT_EDGE_FIRST_TOUCH_OP) .and. .not. defined(_GT_INNER_EDGE_FIRST_TOUCH_OP)
#	define _GT_INNER_EDGE_FIRST_TOUCH_OP	_GT_EDGE_FIRST_TOUCH_OP
#endif

#if defined(_GT_EDGE_LAST_TOUCH_OP) .and. .not. defined(_GT_INNER_EDGE_LAST_TOUCH_OP)
#	define _GT_INNER_EDGE_LAST_TOUCH_OP		_GT_EDGE_LAST_TOUCH_OP
#endif

#if defined(_GT_EDGE_REDUCE_OP) .and. .not. defined(_GT_INNER_EDGE_REDUCE_OP)
#	define _GT_INNER_EDGE_REDUCE_OP		    _GT_EDGE_REDUCE_OP
#endif

#if defined(_GT_NODE_FIRST_TOUCH_OP) .and. .not. defined(_GT_INNER_NODE_FIRST_TOUCH_OP)
#	define _GT_INNER_NODE_FIRST_TOUCH_OP	_GT_NODE_FIRST_TOUCH_OP
#endif

#if defined(_GT_NODE_LAST_TOUCH_OP) .and. .not. defined(_GT_INNER_NODE_LAST_TOUCH_OP)
#	define _GT_INNER_NODE_LAST_TOUCH_OP		_GT_NODE_LAST_TOUCH_OP
#endif

#if defined(_GT_NODE_REDUCE_OP) .and. .not. defined(_GT_INNER_NODE_REDUCE_OP)
#	define _GT_INNER_NODE_REDUCE_OP		    _GT_NODE_REDUCE_OP
#endif

#define _GT_EDGES

PRIVATE
PUBLIC _GT

!Module types:

!> Traversal element ring buffer structure that provides local storage for some of the grid data
type, extends(t_element_base) :: t_traversal_element
	integer (KIND = GRID_SI)                            :: i_cell
	type(num_cell_data_temp)							:: cell_data_temp							!< cell temporary data

#	if defined(_GT_EDGES)
		type(t_edge_data)				                :: next_edge_data            				!< next edge data
#	endif

	type(t_traversal_element), pointer					:: previous, next							!< pointer to previous and next traversal element in the ringbuffer
end type

!> Traversal element structure that provides local storage for some of the grid data
type, extends(t_element_base) :: t_refinement_element
	type(fine_triangle)									    :: cell_geometry_storage
end type

type, extends(num_traversal_data) :: t_section_data
	type(t_adaptive_statistics)								:: stats, current_stats
end type

type, extends(num_traversal_data) :: t_thread_data
	type(t_refinement_element), dimension(2)				:: refinement_elements					!< Temporary refinement elements

    type(t_traversal_element), dimension(8)			        :: src_elements							!< Input element ring buffer (must be 8, because transfer nodes can be referenced back up to 8 elements)
    type(t_traversal_element), dimension(11)        		:: dest_elements						!< Output element ring buffer (must be 11, because transfer nodes can be referenced back up to 8 + 3 new refinement elements)

    type(t_traversal_element), pointer						:: p_src_element => null(), p_dest_element => null()        !< Current source and destination element
end type

type, extends(t_section_data) :: _GT
    type(_GT), pointer                                      :: children_alloc(:) => null(), children(:) => null()			!< section data

    contains

    procedure, pass :: traverse => traverse_in_place
    final :: finalize
end type

contains

subroutine finalize(traversal)
    type(_GT)       :: traversal
	integer         :: i_error

    if (associated(traversal%children_alloc)) then
        deallocate(traversal%children_alloc, stat = i_error); assert_eq(i_error, 0)
    end if
end subroutine


function edge_merge_wrapper_op(local_edges, neighbor_edges) result(l_conform)
    type(t_edge_data), intent(inout)    :: local_edges
    type(t_edge_data), intent(in)       :: neighbor_edges
    logical                             :: l_conform

    assert_eq(local_edges%min_distance, neighbor_edges%min_distance)
    assert(local_edges%owned_locally)

#   if defined(_GT_EDGE_MERGE_OP)
        call _GT_EDGE_MERGE_OP(local_edges, neighbor_edges)
#   endif

    l_conform = .true.
end function

function node_merge_wrapper_op(local_nodes, neighbor_nodes) result(l_conform)
    type(t_node_data), intent(inout)    :: local_nodes
    type(t_node_data), intent(in)       :: neighbor_nodes
    logical                             :: l_conform

    assert_eq(local_nodes%distance, neighbor_nodes%distance)
    assert_veq(local_nodes%position, neighbor_nodes%position)
    assert(local_nodes%owned_locally)

#   if defined(_GT_NODE_MERGE_OP)
        call _GT_NODE_MERGE_OP(local_nodes, neighbor_nodes)
#   endif

    l_conform = .true.
end function

function edge_write_op(local_edges, neighbor_edges) result(l_conform)
    type(t_edge_data), intent(inout)    :: local_edges
    type(t_edge_data), intent(in)       :: neighbor_edges
    logical                             :: l_conform

    assert_eq(local_edges%min_distance, neighbor_edges%min_distance)
    assert(.not. local_edges%owned_locally)
    assert(neighbor_edges%owned_locally)

    local_edges%data_pers = neighbor_edges%data_pers
    local_edges%data_temp = neighbor_edges%data_temp

    l_conform = .true.
end function

function node_write_op(local_nodes, neighbor_nodes) result(l_conform)
    type(t_node_data), intent(inout)    :: local_nodes
    type(t_node_data), intent(in)       :: neighbor_nodes
    logical                             :: l_conform

    assert_eq(local_nodes%distance, neighbor_nodes%distance)
    assert_veq(local_nodes%position, neighbor_nodes%position)
    assert(.not. local_nodes%owned_locally)
    assert(neighbor_nodes%owned_locally)

    local_nodes%data_pers = neighbor_nodes%data_pers
    local_nodes%data_temp = neighbor_nodes%data_temp

    l_conform = .true.
end function

!*****************************************************************
! Generic adaptive traversal
!*****************************************************************

!> Generic adaptive traversal subroutine
!> @author Oliver Meister
subroutine traverse_in_place(traversal, grid)
	class(_GT), target, intent(inout)	                :: traversal
	type(t_grid), intent(inout)							:: grid

	type(t_grid), save							        :: grid_temp
	integer                                             :: imbalance, balance_steps

    !no traversal on empty grids
    !$omp barrier

    !$omp single
    traversal%current_stats%r_integrity_time = -omp_get_wtime()
    !$omp end single

    call conformity_check(grid)

    !$omp single
    traversal%current_stats%r_integrity_time = traversal%current_stats%r_integrity_time + omp_get_wtime()

    traversal%current_stats%r_load_balancing_time = -omp_get_wtime()
    !$omp end single

#	if .not. defined(_GT_INPUT_DEST)
	    !exchange grid sections with neighbors as long as the grid is not balanced
		do balance_steps = 0, huge(1)
	        imbalance = distribute_load(grid)

	        if (imbalance .le. 0) then
	            exit
	        end if
	    end do

   		_log_write(2, "(2X, I0, A)") balance_steps, " load balancing iteration(s) performed."
        !$omp barrier
#	endif

    !$omp single
    traversal%current_stats%r_load_balancing_time = traversal%current_stats%r_load_balancing_time + omp_get_wtime()

    traversal%current_stats%r_allocation_time = -omp_get_wtime()
    !$omp end single

    !$omp barrier

    call rebalance_grid(grid, grid_temp)

    !if necessary, reverse order to match source and destination grid
    if (.not. grid%sections%is_forward()) then
        !$omp barrier
        _log_write(4, '(X, A)') "Reverse destination grid.."
        call grid_temp%reverse()
    end if

    assert_eq(grid%sections%is_forward(), grid_temp%sections%is_forward())

    !$omp single
    traversal%current_stats%r_allocation_time = traversal%current_stats%r_allocation_time + omp_get_wtime()
    !$omp end single

    !refine grid
    call traverse_out_of_place(traversal, grid, grid_temp)

    !$omp barrier

    call grid%destroy()

    !$omp single
    call grid_temp%move(grid)
    !$omp end single
end subroutine

!> Generic adaptive grid traversal
!> @author Oliver Meister
subroutine traverse_out_of_place(traversal, src_grid, dest_grid)
	class(_GT), target, intent(inout)	                        :: traversal
	type(t_grid), intent(inout)									:: src_grid, dest_grid

#	if (_DEBUG_LEVEL > 3)
		_log_write(4, '(2X, A)') "source grid initial state :"
		call src_grid%print()
#	endif

	call traverse_grids(traversal, src_grid, dest_grid)

#	if (_DEBUG_LEVEL > 3)
		_log_write(4, '(2X, A)') "destination grid final state:"
		call dest_grid%print()
#	endif
end subroutine

subroutine traverse_grids(traversal, src_grid, dest_grid)
	class(_GT), intent(inout)	                            :: traversal
	type(t_grid), intent(inout)								:: src_grid, dest_grid

	type(t_grid_section), target							:: src_section
	type(t_grid_section), pointer							:: dest_section
    integer (kind = GRID_SI)			                    :: i_thread, i_src_section, i_dest_section, i_first_local_section, i_last_local_section
    integer (kind = GRID_SI)			                    :: i_src_cell
	integer                                                 :: i_error
	integer (kind = 1)										:: i_color
	type(t_thread_data), target, save						:: thread_traversal
	!$omp threadprivate(thread_traversal)

    if (.not. associated(traversal%children) .or. size(traversal%children) .ne. size(dest_grid%sections%elements)) then
		!$omp barrier    	

		!$omp single
        if (.not. associated(traversal%children_alloc) .or. size(traversal%children_alloc) < size(dest_grid%sections%elements)) then
            if (associated(traversal%children_alloc)) then
                deallocate(traversal%children_alloc, stat = i_error); assert_eq(i_error, 0)
            end if

            allocate(traversal%children_alloc(max(omp_get_num_threads() * dest_grid%i_sections_per_thread, size(dest_grid%sections%elements))), stat = i_error); assert_eq(i_error, 0)
        end if

        traversal%children => traversal%children_alloc(1 : size(dest_grid%sections%elements))
    	!$omp end single
    end if

	if (.not. associated(thread_traversal%p_src_element)) then
		assert(.not. associated(thread_traversal%p_dest_element))
	
    	!create ringbuffers and temporary elements
		call create_ringbuffer(thread_traversal%src_elements)
		call create_ringbuffer(thread_traversal%dest_elements)
		call create_refinement_element(thread_traversal%refinement_elements)

		thread_traversal%p_src_element => thread_traversal%src_elements(1)
		thread_traversal%p_dest_element => thread_traversal%dest_elements(1)
	end if

    i_thread = 1 + omp_get_thread_num()
    call dest_grid%get_local_sections_in_traversal_order(i_first_local_section, i_last_local_section)

    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_traversal_time = -omp_get_wtime()
    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_barrier_time = -omp_get_wtime()

    !$omp barrier

    !$omp single
    call pre_traversal_grid_dest(traversal, dest_grid)

    call prefix_sum(src_grid%sections%elements%last_dest_cell, src_grid%sections%elements%dest_cells)
    call prefix_sum(dest_grid%sections%elements%last_dest_cell, dest_grid%sections%elements%dest_cells)
    !$omp end single

    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_barrier_time = traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_barrier_time + omp_get_wtime()

    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_computation_time = -omp_get_wtime()

    !call pre traversal operator
 	do i_dest_section = i_first_local_section, i_last_local_section
        dest_section => dest_grid%sections%elements(i_dest_section)

        call pre_traversal_dest(traversal%children(i_dest_section), dest_section)
    end do

    !$omp barrier

    !traversal

    if (i_last_local_section .ge. i_first_local_section) then
        dest_section => dest_grid%sections%elements(i_first_local_section)

        !find the source section that contains the first cell of the first destination section
        i_src_section = 1 + ((i_first_local_section - 1) * size(src_grid%sections%elements)) / size(dest_grid%sections%elements)
        i_src_cell = src_grid%sections%elements(i_src_section)%last_dest_cell - src_grid%sections%elements(i_src_section)%dest_cells + 1

        do while (src_grid%sections%elements(i_src_section)%last_dest_cell - src_grid%sections%elements(i_src_section)%dest_cells .gt. dest_section%last_dest_cell - dest_section%dest_cells)
            i_src_section = i_src_section - 1
            assert_le(i_src_section, size(src_grid%sections%elements))
            i_src_cell = src_grid%sections%elements(i_src_section)%last_dest_cell - src_grid%sections%elements(i_src_section)%dest_cells + 1
        end do

        do while (src_grid%sections%elements(i_src_section)%last_dest_cell .le. dest_section%last_dest_cell - dest_section%dest_cells)
            i_src_section = i_src_section + 1
            assert_le(i_src_section, size(src_grid%sections%elements))
            i_src_cell = src_grid%sections%elements(i_src_section)%last_dest_cell - src_grid%sections%elements(i_src_section)%dest_cells + 1
        end do

        src_section = src_grid%sections%elements(i_src_section)

        !traverse all unnecessary elements of the first source section with an empty traversal
        do while (i_src_cell .le. dest_section%last_dest_cell - dest_section%dest_cells)
            i_src_cell = i_src_cell + empty_leaf(thread_traversal, traversal%children(i_first_local_section), src_grid%threads%elements(i_thread), src_section)
        end do
    end if

    !traverse all destination sections
    do i_dest_section = i_first_local_section, i_last_local_section
        dest_section => dest_grid%sections%elements(i_dest_section)

#       if _DEBUG_LEVEL > 4
            _log_write(5, '(2X, A)') "destination section initial state :"
            call dest_section%print()
#       endif

        !traverse all source sections that overlap with the destination section
        do while (i_src_cell .le. dest_section%last_dest_cell)
            assert_le(i_src_section, size(src_grid%sections%elements))
            assert_ge(i_src_cell, src_section%last_dest_cell - src_section%dest_cells + 1)

            !traverse all elements
            do while (i_src_cell .le. min(src_section%last_dest_cell, dest_section%last_dest_cell))
                i_src_cell = i_src_cell + leaf(thread_traversal, traversal%children(i_dest_section), src_grid%threads%elements(i_thread), dest_grid%threads%elements(i_thread), src_section, dest_section)
            end do

            if (i_src_cell .gt. src_section%last_dest_cell .and. i_src_section < size(src_grid%sections%elements)) then
                i_src_section = i_src_section + 1
                src_section = src_grid%sections%elements(i_src_section)
            end if
        end do

#		if defined(_GT_CELL_TO_EDGE_OP)
            thread_traversal%p_dest_element%previous%next_edge_data%rep = _GT_CELL_TO_EDGE_OP(thread_traversal%p_dest_element%previous%t_element_base, thread_traversal%p_dest_element%previous%next_edge_data)
#	    endif

        thread_traversal%p_dest_element%previous%next_edge_data%remove = .false.

        call find_section_boundary_elements(dest_grid%threads%elements(i_thread), dest_section, dest_section%cells%i_current_element, thread_traversal%p_dest_element%previous%next_edge_data)
    end do

    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_computation_time = traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_computation_time + omp_get_wtime()

    !$omp barrier

    call update_distances(dest_grid)

    !$omp barrier

    !$omp single
    assert_veq(decode_distance(dest_grid%t_global_data%min_distance), decode_distance(src_grid%t_global_data%min_distance))
    !$omp end single

	!update communication info
	call update_neighbors(src_grid, dest_grid)

    !$omp barrier

    do i_dest_section = i_first_local_section, i_last_local_section
        call recv_mpi_boundary(dest_grid%sections%elements(i_dest_section))
        call send_mpi_boundary(dest_grid%sections%elements(i_dest_section))
    end do

    !$omp barrier

    !sync and call post traversal operator
    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_sync_time = -omp_get_wtime()
	call sync_boundary(dest_grid, edge_merge_wrapper_op, node_merge_wrapper_op, edge_write_op, node_write_op)
    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_sync_time = traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_sync_time + omp_get_wtime()

    !$omp barrier

    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_computation_time = traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_computation_time - omp_get_wtime()

    do i_dest_section = i_first_local_section, i_last_local_section
        call post_traversal_dest(traversal%children(i_dest_section), dest_grid%sections%elements(i_dest_section))
	end do

    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_computation_time = traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_computation_time + omp_get_wtime()

    !$omp barrier

    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_barrier_time = traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_barrier_time - omp_get_wtime()

    !$omp single
    call post_traversal_grid_dest(traversal, dest_grid)
    !$omp end single

    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_barrier_time = traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_barrier_time + omp_get_wtime()
    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_traversal_time = traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_traversal_time + omp_get_wtime()

	call set_stats_dest(traversal%children(i_first_local_section : i_last_local_section), dest_grid%sections%elements(i_first_local_section : i_last_local_section))

    !$omp barrier

#   if defined(_GT_OUTPUT_SRC)
        call src_grid%reverse()
#   endif

	call dest_grid%reverse()

    !$omp single
        call traversal%current_stats%reduce(traversal%children%current_stats)
        traversal%stats = traversal%stats + traversal%current_stats
        dest_grid%stats = dest_grid%stats + traversal%current_stats
    !$omp end single
end subroutine

function leaf(thread_traversal, traversal, src_thread, dest_thread, src_section, dest_section) result(i_dest_cells)
    type(t_thread_data), intent(inout)          :: thread_traversal
    type(_GT), intent(inout)                    :: traversal
    type(t_grid_thread), intent(inout)          :: src_thread, dest_thread
    type(t_grid_section), intent(inout)         :: src_section, dest_section
    integer  (kind = GRID_SI)                   :: i_dest_cells

    type(t_traversal_element), pointer          :: p_dest_element_2, p_dest_element_3, p_dest_element_4
    integer                                     :: i

	call init_src_element(src_section, thread_traversal%p_src_element)
	call init_dest_element(dest_section, thread_traversal%p_dest_element)

	select case (thread_traversal%p_src_element%cell%geometry%refinement)
		case (-1)
			!coarsening
			!p_src_element_2 => thread_traversal%p_src_element%next

			call init_src_element(src_section, thread_traversal%p_src_element%next)

#			if .not. defined(_GT_INPUT_DEST)
				call create_parent_cell(thread_traversal%p_src_element%cell%geometry, thread_traversal%p_src_element%next%cell%geometry, thread_traversal%p_dest_element%cell%geometry)
#			endif

			call read_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

			call read_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

#			if .not. defined(_GT_INPUT_DEST)
				thread_traversal%p_dest_element%nodes(1)%ptr%position = thread_traversal%p_src_element%nodes(1)%ptr%position
				thread_traversal%p_dest_element%nodes(2)%ptr%position = thread_traversal%p_src_element%nodes(3)%ptr%position
				thread_traversal%p_dest_element%nodes(3)%ptr%position = 2.0_GRID_SR * thread_traversal%p_src_element%nodes(2)%ptr%position - thread_traversal%p_src_element%nodes(1)%ptr%position
#			endif

			call _GT_COARSEN_OP(traversal, dest_section, thread_traversal%p_src_element, thread_traversal%p_dest_element, [1])
			call write_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

			thread_traversal%p_src_element => thread_traversal%p_src_element%next

			call read_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

#			if .not. defined(_GT_INPUT_DEST)
				thread_traversal%p_dest_element%nodes(3)%ptr%position = thread_traversal%p_src_element%nodes(3)%ptr%position
#			endif

			call _GT_COARSEN_OP(traversal, dest_section, thread_traversal%p_src_element, thread_traversal%p_dest_element, [2])
			call write_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

			call write_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

			thread_traversal%p_src_element => thread_traversal%p_src_element%next
			thread_traversal%p_dest_element => thread_traversal%p_dest_element%next

			i_dest_cells = 1
		case (0)
			!transfer
#			if .not. defined(_GT_INPUT_DEST)
				call copy_cell(thread_traversal%p_src_element%cell%geometry, thread_traversal%p_dest_element%cell%geometry)
#			endif

			call read_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

			call read_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

#			if .not. defined(_GT_INPUT_DEST)
				do i = 1, 3
					thread_traversal%p_dest_element%nodes(i)%ptr%position = thread_traversal%p_src_element%nodes(i)%ptr%position
				end do
#			endif

			call _GT_TRANSFER_OP(traversal, dest_section, thread_traversal%p_src_element, thread_traversal%p_dest_element)
			call write_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

			call write_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

			thread_traversal%p_src_element => thread_traversal%p_src_element%next
			thread_traversal%p_dest_element => thread_traversal%p_dest_element%next

			i_dest_cells = 1
		case (1)
			!single refinement
			p_dest_element_2 => thread_traversal%p_dest_element%next

			call init_dest_element(dest_section, p_dest_element_2)

#			if .not. defined(_GT_INPUT_DEST)
				call create_child_cells(thread_traversal%p_src_element%cell%geometry, thread_traversal%p_dest_element%cell%geometry, p_dest_element_2%cell%geometry)
#			endif

			call read_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

			call read_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

#			if .not. defined(_GT_INPUT_DEST)
				thread_traversal%p_dest_element%nodes(1)%ptr%position = thread_traversal%p_src_element%nodes(1)%ptr%position
				thread_traversal%p_dest_element%nodes(2)%ptr%position = 0.5_GRID_SR * (thread_traversal%p_src_element%nodes(1)%ptr%position + thread_traversal%p_src_element%nodes(3)%ptr%position)
				thread_traversal%p_dest_element%nodes(3)%ptr%position = thread_traversal%p_src_element%nodes(2)%ptr%position
#			endif

			call _GT_REFINE_OP(traversal, dest_section, thread_traversal%p_src_element, thread_traversal%p_dest_element, [1])
			call write_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

			call read_dest_element(traversal, dest_thread, dest_section, p_dest_element_2)

#			if .not. defined(_GT_INPUT_DEST)
				p_dest_element_2%nodes(3)%ptr%position = thread_traversal%p_src_element%nodes(3)%ptr%position
#			endif

			call _GT_REFINE_OP(traversal, dest_section, thread_traversal%p_src_element, p_dest_element_2, [2])
			call write_dest_element(traversal, dest_thread, dest_section, p_dest_element_2)

			call write_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

			thread_traversal%p_src_element => thread_traversal%p_src_element%next
			thread_traversal%p_dest_element => p_dest_element_2%next

			i_dest_cells = 2
		case (2)
			!single refinement + first child refinement
			p_dest_element_2 => thread_traversal%p_dest_element%next
			p_dest_element_3 => p_dest_element_2%next

			call init_dest_element(dest_section, p_dest_element_2)
			call init_dest_element(dest_section, p_dest_element_3)

#			if .not. defined(_GT_INPUT_DEST)
				call create_child_cells(thread_traversal%p_src_element%cell%geometry, thread_traversal%refinement_elements(1)%cell%geometry, p_dest_element_3%cell%geometry)
				call create_child_cells(thread_traversal%refinement_elements(1)%cell%geometry, thread_traversal%p_dest_element%cell%geometry, p_dest_element_2%cell%geometry)
#			endif

			call read_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

			call read_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

#			if .not. defined(_GT_INPUT_DEST)
				thread_traversal%p_dest_element%nodes(1)%ptr%position = thread_traversal%p_src_element%nodes(1)%ptr%position
				thread_traversal%p_dest_element%nodes(2)%ptr%position = 0.5_GRID_SR * (thread_traversal%p_src_element%nodes(1)%ptr%position + thread_traversal%p_src_element%nodes(2)%ptr%position)
				thread_traversal%p_dest_element%nodes(3)%ptr%position = 0.5_GRID_SR * (thread_traversal%p_src_element%nodes(1)%ptr%position + thread_traversal%p_src_element%nodes(3)%ptr%position)
#			endif

			call _GT_REFINE_OP(traversal, dest_section, thread_traversal%p_src_element, thread_traversal%p_dest_element, [1, 1])
			call write_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

			call read_dest_element(traversal, dest_thread, dest_section, p_dest_element_2)

#			if .not. defined(_GT_INPUT_DEST)
				p_dest_element_2%nodes(3)%ptr%position = thread_traversal%p_src_element%nodes(2)%ptr%position
#			endif

			call _GT_REFINE_OP(traversal, dest_section, thread_traversal%p_src_element, p_dest_element_2, [1, 2])
			call write_dest_element(traversal, dest_thread, dest_section, p_dest_element_2)

			call read_dest_element(traversal, dest_thread, dest_section, p_dest_element_3)

#			if .not. defined(_GT_INPUT_DEST)
				p_dest_element_3%nodes(3)%ptr%position = thread_traversal%p_src_element%nodes(3)%ptr%position
#			endif

			call _GT_REFINE_OP(traversal, dest_section, thread_traversal%p_src_element, p_dest_element_3, [2])
			call write_dest_element(traversal, dest_thread, dest_section, p_dest_element_3)

			call write_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

			thread_traversal%p_src_element => thread_traversal%p_src_element%next
			thread_traversal%p_dest_element => p_dest_element_3%next

			i_dest_cells = 3
		case (3)
			!single refinement + second child refinement
			p_dest_element_2 => thread_traversal%p_dest_element%next
			p_dest_element_3 => p_dest_element_2%next

			call init_dest_element(dest_section, p_dest_element_2)
			call init_dest_element(dest_section, p_dest_element_3)

#			if .not. defined(_GT_INPUT_DEST)
				call create_child_cells(thread_traversal%p_src_element%cell%geometry, thread_traversal%p_dest_element%cell%geometry, thread_traversal%refinement_elements(1)%cell%geometry)
				call create_child_cells(thread_traversal%refinement_elements(1)%cell%geometry, p_dest_element_2%cell%geometry, p_dest_element_3%cell%geometry)
#			endif

			call read_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

			call read_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

#			if .not. defined(_GT_INPUT_DEST)
				thread_traversal%p_dest_element%nodes(1)%ptr%position = thread_traversal%p_src_element%nodes(1)%ptr%position
				thread_traversal%p_dest_element%nodes(2)%ptr%position = 0.5_GRID_SR * (thread_traversal%p_src_element%nodes(1)%ptr%position + thread_traversal%p_src_element%nodes(3)%ptr%position)
				thread_traversal%p_dest_element%nodes(3)%ptr%position = thread_traversal%p_src_element%nodes(2)%ptr%position
#			endif

			call _GT_REFINE_OP(traversal, dest_section, thread_traversal%p_src_element, thread_traversal%p_dest_element, [1])
			call write_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

			call read_dest_element(traversal, dest_thread, dest_section, p_dest_element_2)

#			if .not. defined(_GT_INPUT_DEST)
				p_dest_element_2%nodes(2)%ptr%position = 0.5_GRID_SR * (thread_traversal%p_src_element%nodes(2)%ptr%position + thread_traversal%p_src_element%nodes(3)%ptr%position)
#			endif

			call _GT_REFINE_OP(traversal, dest_section, thread_traversal%p_src_element, p_dest_element_2, [2, 1])
			call write_dest_element(traversal, dest_thread, dest_section, p_dest_element_2)

			call read_dest_element(traversal, dest_thread, dest_section, p_dest_element_3)

#			if .not. defined(_GT_INPUT_DEST)
				p_dest_element_3%nodes(3)%ptr%position = thread_traversal%p_src_element%nodes(3)%ptr%position
#			endif

			call _GT_REFINE_OP(traversal, dest_section, thread_traversal%p_src_element, p_dest_element_3, [2, 2])
			call write_dest_element(traversal, dest_thread, dest_section, p_dest_element_3)

			call write_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

			thread_traversal%p_src_element => thread_traversal%p_src_element%next
			thread_traversal%p_dest_element => p_dest_element_3%next

			i_dest_cells = 3
		case (4)
			!double refinement
			p_dest_element_2 => thread_traversal%p_dest_element%next
			p_dest_element_3 => p_dest_element_2%next
			p_dest_element_4 => p_dest_element_3%next

			call init_dest_element(dest_section, p_dest_element_2)
			call init_dest_element(dest_section, p_dest_element_3)
			call init_dest_element(dest_section, p_dest_element_4)

#			if .not. defined(_GT_INPUT_DEST)
				call create_child_cells(thread_traversal%p_src_element%cell%geometry, thread_traversal%refinement_elements(1)%cell%geometry, thread_traversal%refinement_elements(2)%cell%geometry)
				call create_child_cells(thread_traversal%refinement_elements(1)%cell%geometry, thread_traversal%p_dest_element%cell%geometry, p_dest_element_2%cell%geometry)
				call create_child_cells(thread_traversal%refinement_elements(2)%cell%geometry, p_dest_element_3%cell%geometry, p_dest_element_4%cell%geometry)
#			endif

			call read_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

			call read_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

#			if .not. defined(_GT_INPUT_DEST)
				thread_traversal%p_dest_element%nodes(1)%ptr%position = thread_traversal%p_src_element%nodes(1)%ptr%position
				thread_traversal%p_dest_element%nodes(2)%ptr%position = 0.5_GRID_SR * (thread_traversal%p_src_element%nodes(1)%ptr%position + thread_traversal%p_src_element%nodes(2)%ptr%position)
				thread_traversal%p_dest_element%nodes(3)%ptr%position = 0.5_GRID_SR * (thread_traversal%p_src_element%nodes(1)%ptr%position + thread_traversal%p_src_element%nodes(3)%ptr%position)
#			endif

			call _GT_REFINE_OP(traversal, dest_section, thread_traversal%p_src_element, thread_traversal%p_dest_element, [1, 1])
			call write_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

			call read_dest_element(traversal, dest_thread, dest_section, p_dest_element_2)

#			if .not. defined(_GT_INPUT_DEST)
				p_dest_element_2%nodes(3)%ptr%position = thread_traversal%p_src_element%nodes(2)%ptr%position
#			endif

			call _GT_REFINE_OP(traversal, dest_section, thread_traversal%p_src_element, p_dest_element_2, [1, 2])
			call write_dest_element(traversal, dest_thread, dest_section, p_dest_element_2)

			call read_dest_element(traversal, dest_thread, dest_section, p_dest_element_3)

#			if .not. defined(_GT_INPUT_DEST)
				p_dest_element_3%nodes(2)%ptr%position = 0.5_GRID_SR * (thread_traversal%p_src_element%nodes(2)%ptr%position + thread_traversal%p_src_element%nodes(3)%ptr%position)
#			endif

			call _GT_REFINE_OP(traversal, dest_section, thread_traversal%p_src_element, p_dest_element_3, [2, 1])
			call write_dest_element(traversal, dest_thread, dest_section, p_dest_element_3)

			call read_dest_element(traversal, dest_thread, dest_section, p_dest_element_4)

#			if .not. defined(_GT_INPUT_DEST)
				p_dest_element_4%nodes(3)%ptr%position = thread_traversal%p_src_element%nodes(3)%ptr%position
#			endif

			call _GT_REFINE_OP(traversal, dest_section, thread_traversal%p_src_element, p_dest_element_4, [2, 2])
			call write_dest_element(traversal, dest_thread, dest_section, p_dest_element_4)

			call write_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

			thread_traversal%p_src_element => thread_traversal%p_src_element%next
			thread_traversal%p_dest_element => p_dest_element_4%next

			i_dest_cells = 4
	end select
end function

function empty_leaf(thread_traversal, traversal, src_thread, src_section) result(i_dest_cells)
    type(t_thread_data), intent(inout)          :: thread_traversal
    type(_GT), intent(inout)                    :: traversal
    type(t_grid_thread), intent(inout)          :: src_thread
    type(t_grid_section), intent(inout)         :: src_section
    integer(kind = GRID_SI)                     :: i_dest_cells

	call init_src_element(src_section, thread_traversal%p_src_element)
	select case (thread_traversal%p_src_element%cell%geometry%i_edge_types)
		case (INNER_OLD)
			call read_src_oon(traversal, src_thread, src_section, thread_traversal%p_src_element)
			call write_src_oon(traversal, src_thread, src_section, thread_traversal%p_src_element)
		case (INNER_NEW)
			call read_src_onn(traversal, src_thread, src_section, thread_traversal%p_src_element)
			call write_src_onn(traversal, src_thread, src_section, thread_traversal%p_src_element)
		case default
			call read_src(traversal, src_thread, src_section, thread_traversal%p_src_element)
			call write_src(traversal, src_thread, src_section, thread_traversal%p_src_element)
	end select

    select case (thread_traversal%p_src_element%cell%geometry%refinement)
		case (-1)
            thread_traversal%p_src_element => thread_traversal%p_src_element%next

            call init_src_element(src_section, thread_traversal%p_src_element)
            select case (thread_traversal%p_src_element%cell%geometry%i_edge_types)
                case (INNER_OLD)
                    call read_src_oon(traversal, src_thread, src_section, thread_traversal%p_src_element)
                    call write_src_oon(traversal, src_thread, src_section, thread_traversal%p_src_element)
                case (INNER_NEW)
                    call read_src_onn(traversal, src_thread, src_section, thread_traversal%p_src_element)
                    call write_src_onn(traversal, src_thread, src_section, thread_traversal%p_src_element)
                case default
                    call read_src(traversal, src_thread, src_section, thread_traversal%p_src_element)
                    call write_src(traversal, src_thread, src_section, thread_traversal%p_src_element)
            end select

            i_dest_cells = 1
		case (0)
            i_dest_cells = 1
		case (1)
            i_dest_cells = 2
        case (2, 3)
            i_dest_cells = 3
		case (4)
            i_dest_cells = 4
	end select

    thread_traversal%p_src_element => thread_traversal%p_src_element%next
end function

subroutine copy_cell(src_cell, dest_cell)
	type(fine_triangle), intent(in)						:: src_cell
	type(fine_triangle), intent(inout)					:: dest_cell

	dest_cell = src_cell

	call dest_cell%set_edge_types(OLD, iand(1, src_cell%get_color_edge_type()), NEW)

#	if (_DEBUG_LEVEL > 5)
		_log_write(4, '(A, A)') "  src        (in): ", src_cell%to_string()
		_log_write(4, '(A, A)') "  dest      (out): ", dest_cell%to_string()
#	endif
end subroutine

subroutine create_parent_cell(first_child_cell, second_child_cell, parent_cell)
	type(fine_triangle), intent(in)						:: first_child_cell, second_child_cell
	type(fine_triangle), intent(inout)					:: parent_cell

	integer(kind = 1), dimension(-8:8), parameter 	    :: i_plotter_parent_type = [ 3, 2, 1, 8, 7, 6, 5, 4,     0,  -6, -7, -8, -1, -2, -3, -4, -5 ]
	integer(kind = 1)			 						:: i_color_edge_type
	integer(kind = 1)								 	:: i

	!the parent turtle grammar type can be computed by a simple xor (it's not a symmetrical operator though, so check if the order is correct)
	parent_cell%i_turtle_type = ieor(first_child_cell%i_turtle_type, second_child_cell%i_turtle_type)
	assert_gt(first_child_cell%i_turtle_type, second_child_cell%i_turtle_type)

	!the parent plotter grammar type is always the same as the first grandchild's type
	parent_cell%i_plotter_type = i_plotter_parent_type(first_child_cell%i_plotter_type)
	parent_cell%i_depth = first_child_cell%i_depth - 1

	select case (parent_cell%i_turtle_type)
		case (K)
			parent_cell%l_color_edge_color = second_child_cell%l_color_edge_color
			i_color_edge_type = iand(1, second_child_cell%get_color_edge_type())
		case (V, H)
			parent_cell%l_color_edge_color = first_child_cell%l_color_edge_color
			i_color_edge_type = iand(1, first_child_cell%get_color_edge_type())
	end select

	call parent_cell%set_edge_types(OLD, i_color_edge_type, NEW)

#	if (_DEBUG_LEVEL > 5)
		_log_write(5, '(A, A)') "  1st child  (in): ", first_child_cell%to_string()
		_log_write(5, '(A, A)') "  2nd child  (in): ", second_child_cell%to_string()
		_log_write(5, '(A, A)') "  parent    (out): ", parent_cell%to_string()
#	endif
end subroutine

subroutine create_child_cells(parent_cell, first_child_cell, second_child_cell)
	type(fine_triangle), intent(in)						:: parent_cell
	type(fine_triangle), intent(inout)					:: first_child_cell, second_child_cell

	integer(kind = 1), dimension(3, 2), parameter		:: i_turtle_child_type = reshape([ H, H, V, V, K, K ], [ 3, 2 ])
	integer(kind = 1), dimension(-8:8), parameter 	    :: i_plotter_child_type = [ 3, 2, 1, 8, 7, 6, 5, 4,     0,  -6, -7, -8, -1, -2, -3, -4, -5 ]
	integer(kind = 1)			 						:: i_color_edge_type
	integer(kind = 1)								 	:: i

	first_child_cell%i_turtle_type = i_turtle_child_type(parent_cell%i_turtle_type, 1)
	second_child_cell%i_turtle_type = i_turtle_child_type(parent_cell%i_turtle_type, 2)

	!check for correctness
	assert_eq(parent_cell%i_turtle_type, ieor(first_child_cell%i_turtle_type, second_child_cell%i_turtle_type))
	assert_gt(first_child_cell%i_turtle_type, second_child_cell%i_turtle_type)

	first_child_cell%i_plotter_type = i_plotter_child_type(parent_cell%i_plotter_type)
	first_child_cell%i_depth = parent_cell%i_depth + 1
	second_child_cell%i_plotter_type = -i_plotter_child_type(-parent_cell%i_plotter_type)
	second_child_cell%i_depth = parent_cell%i_depth + 1

	i_color_edge_type = iand(1, parent_cell%get_color_edge_type())

	select case (parent_cell%i_turtle_type)
		case (K)
			call first_child_cell%set_edge_types(OLD, NEW, NEW)
			call second_child_cell%set_edge_types(OLD, i_color_edge_type, NEW)
			first_child_cell%l_color_edge_color = .not. parent_cell%l_color_edge_color
			second_child_cell%l_color_edge_color = parent_cell%l_color_edge_color
	 	case (V)
			call first_child_cell%set_edge_types(OLD, i_color_edge_type, NEW)
			call second_child_cell%set_edge_types(OLD, i_color_edge_type, NEW)
			first_child_cell%l_color_edge_color = parent_cell%l_color_edge_color
			second_child_cell%l_color_edge_color = parent_cell%l_color_edge_color
	 	case (H)
			call first_child_cell%set_edge_types(OLD, i_color_edge_type, NEW)
			call second_child_cell%set_edge_types(OLD, OLD, NEW)
			first_child_cell%l_color_edge_color = parent_cell%l_color_edge_color
			second_child_cell%l_color_edge_color = .not. parent_cell%l_color_edge_color
	end select

#	if (_DEBUG_LEVEL > 5)
		_log_write(5, '(A, A)') "  parent     (in): ", parent_cell%to_string()
		_log_write(5, '(A, A)') "  1st child (out): ", first_child_cell%to_string()
		_log_write(5, '(A, A)') "  2nd child (out): ", second_child_cell%to_string()
#	endif
end subroutine

elemental subroutine create_refinement_element(element)
	type(t_refinement_element), target, intent(inout)		:: element

	element%cell%geometry => element%cell_geometry_storage
end subroutine

subroutine create_ringbuffer(elements)
	type(t_traversal_element), dimension(:), target, intent(inout)		:: elements

	integer (kind = GRID_SI)											:: i

	do i = 1, size(elements)
		elements(i)%previous => elements(mod(i + size(elements) - 2, size(elements)) + 1)
		elements(i)%next => elements(mod(i, size(elements)) + 1)

		nullify(elements(i)%cell%geometry)
		nullify(elements(i)%cell%data_pers)
		elements(i)%cell%data_temp => elements(i)%cell_data_temp

		nullify(elements(i)%color_node_out%ptr)
		nullify(elements(i)%transfer_node%ptr)
		nullify(elements(i)%color_node_in%ptr)

#		if defined(_GT_EDGES)
            nullify(elements(i)%color_edge%ptr)
            elements(i)%previous_edge%ptr => elements(i)%previous%next_edge_data
            elements(i)%next_edge%ptr => elements(i)%next_edge_data
#		endif
	end do
end subroutine

!************
!I/O routines
!************

!****************
!Destination grid operations
!****************

#undef _GT_NO_COORDS
#undef _GT_ELEMENT_TO_EDGE_OP
#undef _GT_SKELETON_OP

#if defined(_GT_INPUT_DEST)
#	define _GT_INPUT
#endif

#	define _GT_OUTPUT
#	define _GT_PREFIX						dest

#	include "Tools_adaptive_traversal.f90"

#	undef _GT_INPUT
#	undef _GT_OUTPUT
#	undef _GT_PREFIX

subroutine init_dest_element(section, element)
	type(t_grid_section), intent(inout)             :: section
	type(t_traversal_element), intent(inout)		:: element

	call init_dest(section, element)
	element%i_cell = section%cells%i_current_element
end subroutine

subroutine read_dest_element(traversal, thread, section, element)
	type(_GT), intent(inout)                        :: traversal
	type(t_grid_thread), intent(inout)              :: thread
	type(t_grid_section), intent(inout)             :: section
	type(t_traversal_element), intent(inout)		:: element

    integer (kind = GRID_SI)                        :: i_empty

    if (element%i_cell == 1) then
        call element%cell%geometry%set_previous_edge_type(OLD_BND)
    end if

    select case (element%cell%geometry%get_color_edge_type())
        case (OLD)
            if(thread%indices_stack(element%cell%geometry%l_color_edge_color)%is_empty()) then
                call element%cell%geometry%set_color_edge_type(OLD_BND)
            else
                call thread%indices_stack(element%cell%geometry%l_color_edge_color)%pop(i_empty)
            end if
        case (NEW)
            call thread%indices_stack(element%cell%geometry%l_color_edge_color)%push(element%i_cell)
    end select

	call read_dest(traversal, thread, section, element)
end subroutine

subroutine write_dest_element(traversal, thread, section, element)
	type(_GT), intent(inout)                        :: traversal
	type(t_grid_thread), intent(inout)              :: thread
	type(t_grid_section), intent(inout)             :: section
	type(t_traversal_element), intent(inout)		:: element

    element%next_edge%ptr%remove = .false.
    element%color_edge%ptr%remove = .false.

    select case (element%cell%geometry%i_turtle_type)
        case (K)
            element%previous_edge%ptr%remove = .true.
        case (H)
            element%next_edge%ptr%remove = .true.
    end select

    select case (element%cell%geometry%i_edge_types)
        case (INNER_OLD)
            call write_dest_oon(traversal, thread, section, element)
        case (INNER_NEW)
            call write_dest_onn(traversal, thread, section, element)
        case (INNER_OLD_BND)
            call write_dest_odn(traversal, thread, section, element)
        case (INNER_NEW_BND)
            call write_dest_obn(traversal, thread, section, element)
        case (FIRST_NEW, FIRST_OLD_BND, FIRST_NEW_BND)
            element%previous_edge%ptr%remove = .false.
            call write_dest(traversal, thread, section, element)
        case (LAST_OLD, LAST_OLD_BND, LAST_NEW_BND)
            element%next_edge%ptr%remove = .false.
            call write_dest(traversal, thread, section, element)
        case (SINGLE_OLD_BND, SINGLE_NEW_BND)
            element%previous_edge%ptr%remove = .false.
            element%next_edge%ptr%remove = .false.
            call write_dest(traversal, thread, section, element)
    end select
end subroutine

#undef _GT_CELL_FIRST_TOUCH_OP
#undef _GT_CELL_LAST_TOUCH_OP
#undef _GT_CELL_REDUCE_OP
#undef _GT_EDGE_FIRST_TOUCH_OP
#undef _GT_EDGE_LAST_TOUCH_OP
#undef _GT_EDGE_REDUCE_OP
#undef _GT_NODE_FIRST_TOUCH_OP
#undef _GT_NODE_LAST_TOUCH_OP
#undef _GT_NODE_REDUCE_OP
#undef _GT_INNER_EDGE_FIRST_TOUCH_OP
#undef _GT_INNER_EDGE_LAST_TOUCH_OP
#undef _GT_INNER_EDGE_REDUCE_OP
#undef _GT_INNER_NODE_FIRST_TOUCH_OP
#undef _GT_INNER_NODE_LAST_TOUCH_OP
#undef _GT_INNER_NODE_REDUCE_OP
#undef _GT_EDGE_MERGE_OP
#undef _GT_NODE_MERGE_OP
#undef _GT_PRE_TRAVERSAL_GRID_OP
#undef _GT_POST_TRAVERSAL_GRID_OP
#undef _GT_PRE_TRAVERSAL_OP
#undef _GT_POST_TRAVERSAL_OP

!****************
!Source grid operations
!****************

#if defined(_GT_OUTPUT_SRC)
#	define _GT_OUTPUT
#endif

#	define _GT_INPUT
#	define _GT_PREFIX						src

#	include "Tools_adaptive_traversal.f90"

#	undef _GT_PREFIX
#   undef  _GT_INPUT
#   undef  _GT_OUTPUT

subroutine init_src_element(section, element)
	type(t_grid_section), intent(inout)             :: section
	type(t_traversal_element), intent(inout)		:: element

	call init_src(section, element)
	element%i_cell = section%cells%i_current_element
end subroutine

subroutine read_src_element(traversal, thread, section, element)
	type(_GT), intent(inout)                        :: traversal
	type(t_grid_thread), intent(inout)              :: thread
	type(t_grid_section), intent(inout)             :: section
	type(t_traversal_element), intent(inout)		:: element

	call read_src(traversal, thread, section, element)
end subroutine

subroutine write_src_element(traversal, thread, section, element)
	type(_GT), intent(inout)                        :: traversal
	type(t_grid_thread), intent(inout)              :: thread
	type(t_grid_section), intent(inout)             :: section
	type(t_traversal_element), intent(inout)		:: element

	call write_src(traversal, thread, section, element)
end subroutine

!---

#undef _GT_NAME
#undef _GT_NODES
#undef _GT_EDGES

#undef _GT_TRANSFER_OP
#undef _GT_REFINE_OP
#undef _GT_COARSEN_OP
