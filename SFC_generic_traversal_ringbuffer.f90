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
#if defined(_GT_ELEMENT_OP) .and. .not. defined(_GT_INNER_ELEMENT_OP)
#	define _GT_INNER_ELEMENT_OP		        _GT_ELEMENT_OP
#endif

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

PRIVATE
PUBLIC _GT

!Module types:

!> Traversal element ring buffer structure that provides local storage for some of the grid data
type, extends(t_element_base) :: t_traversal_element
	type(num_cell_data_temp)							    :: cell_data_temp							!< cell temporary data

#	if defined(_GT_EDGES)
		type(t_edge_data)				                    :: next_edge_data            				!< next crossed edge temp + local data
#	endif

	type(t_traversal_element), pointer					    :: previous, next							!< pointer to previous and next traversal element in the ringbuffer
end type

type, extends(num_traversal_data) :: t_section_data
	type(t_statistics)										:: stats, current_stats
end type

type, extends(num_traversal_data) :: t_thread_data
    type(t_traversal_element), dimension(8)     			:: elements									!< Element ring buffer (must contain at least 8 elements, because transfer nodes can be referenced back up to 8 elements)
    type(t_traversal_element), pointer						:: p_current_element => null()				!< Current element
end type

type, extends(t_section_data) :: _GT
    type(_GT), pointer                                      :: children_alloc(:) => null(), children(:) => null()			!< section data

    contains

    procedure, pass :: traverse => traverse_grid
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
! Generic traversal
!*****************************************************************

!> Generic iterative traversal subroutine
!> @author Oliver Meister
subroutine traverse_grid(traversal, grid)
	class(_GT), intent(inout)	                        :: traversal
	class(t_grid), intent(inout)					    :: grid

	integer (kind = GRID_SI)                            :: i_section, i_thread, i_first_local_section, i_last_local_section
	integer (kind = GRID_SI)                            :: i_error
	type(t_thread_data), target, save					:: thread_traversal
	!$omp threadprivate(thread_traversal)

    if (.not. associated(traversal%children) .or. size(traversal%children) .ne. size(grid%sections%elements_alloc)) then
        !$omp barrier

        !$omp single
        if (.not. associated(traversal%children_alloc) .or. size(traversal%children_alloc) < size(grid%sections%elements_alloc)) then
            if (associated(traversal%children_alloc)) then
                deallocate(traversal%children_alloc, stat = i_error); assert_eq(i_error, 0)
            end if

            allocate(traversal%children_alloc(max(omp_get_num_threads() * grid%i_sections_per_thread, size(grid%sections%elements_alloc))), stat = i_error); assert_eq(i_error, 0)
        end if

        traversal%children => traversal%children_alloc(1 : size(grid%sections%elements_alloc))
        !$omp end single
    end if

    if (.not. associated(thread_traversal%p_current_element)) then
        call create_ringbuffer(thread_traversal%elements)
        thread_traversal%p_current_element => thread_traversal%elements(1)
    end if

    call grid%get_local_sections(i_first_local_section, i_last_local_section)

    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_traversal_time = -omp_get_wtime()

    !$omp single
    call pre_traversal_grid(traversal, grid)
    !$omp end single nowait

    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_computation_time = 0.0

#   if defined(_GT_SKELETON_OP)
        traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_computation_time = traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_computation_time - omp_get_wtime()

        do i_section = i_first_local_section, i_last_local_section
            assert_eq(i_section, grid%sections%elements_alloc(i_section)%index)

            call boundary_skeleton(traversal%children(i_section), grid%sections%elements_alloc(i_section))
        end do

        traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_computation_time = traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_computation_time + omp_get_wtime()
#   endif

    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_barrier_time = -omp_get_wtime()

    !$omp barrier

    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_barrier_time = traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_barrier_time + omp_get_wtime()

    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_computation_time = traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_computation_time - omp_get_wtime()

	do i_section = i_first_local_section, i_last_local_section
        call recv_mpi_boundary(grid%sections%elements_alloc(i_section))
    end do

	do i_section = i_first_local_section, i_last_local_section
        !$omp task default(shared) firstprivate(i_section) private(i_thread) mergeable
            i_thread = 1 + omp_get_thread_num()
            call pre_traversal(traversal%children(i_section), grid%sections%elements_alloc(i_section))
            call traverse_section(thread_traversal, traversal%children(i_section), grid%threads%elements(i_thread), grid%sections%elements_alloc(i_section))
            call send_mpi_boundary(grid%sections%elements_alloc(i_section))
        !$omp end task
    end do

    !$omp taskwait

    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_computation_time = traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_computation_time + omp_get_wtime()

    call grid%reverse()

    !sync and call post traversal operator
    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_sync_time = -omp_get_wtime()
	call sync_boundary(grid, edge_merge_wrapper_op, node_merge_wrapper_op, edge_write_op, node_write_op)
    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_sync_time = traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_sync_time + omp_get_wtime()

    !call post traversal operator
    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_computation_time = traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_computation_time - omp_get_wtime()

	do i_section = i_first_local_section, i_last_local_section
        assert_eq(i_section, grid%sections%elements_alloc(i_section)%index)
        call post_traversal(traversal%children(i_section), grid%sections%elements_alloc(i_section))
    end do

    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_computation_time = traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_computation_time + omp_get_wtime()

    !$omp barrier

    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_barrier_time = traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_barrier_time - omp_get_wtime()

    !$omp single
    call post_traversal_grid(traversal, grid)
    !$omp end single

    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_barrier_time = traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_barrier_time + omp_get_wtime()
    traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_traversal_time = traversal%children(i_first_local_section : i_last_local_section)%current_stats%r_traversal_time + omp_get_wtime()
	call set_stats(traversal%children(i_first_local_section : i_last_local_section), grid%sections%elements_alloc(i_first_local_section : i_last_local_section))

    !$omp barrier

    !$omp single
        call traversal%current_stats%reduce(traversal%children%current_stats)
        traversal%stats = traversal%stats + traversal%current_stats
        grid%stats = grid%stats + traversal%current_stats
    !$omp end single
end subroutine

!> Generic iterative traversal subroutine
!> @author Oliver Meister
subroutine traverse_section(thread_traversal, traversal, thread, section)
	type(t_thread_data), target, intent(inout)	        :: thread_traversal
	class(_GT), target, intent(inout)	                :: traversal
	type(t_grid_thread), intent(inout)					:: thread
	type(t_grid_section), intent(inout)					:: section

	! local variables
	integer (kind = GRID_SI)							:: i
	type(t_traversal_element), pointer					:: p_next_element

#	if (_DEBUG_LEVEL > 4)
		_log_write(5, '(2X, A)') "section input state :"
		call section%print()
#	endif

#	if (_DEBUG_LEVEL > 5)
		_log_write(6, '(2X, A)') "input cells:"
		do i = lbound(section%cells%elements, 1), ubound(section%cells%elements, 1)
			_log_write(6, '(3X, I0, X, A)') i, section%cells%elements(i)%to_string()
		end do
		_log_write(6, '(A)') ""
#	endif

	if (size(section%cells%elements) > 0) then
        assert_ge(size(section%cells%elements), 2)
		call init(section, thread_traversal%p_current_element)

		!process first element
        p_next_element => thread_traversal%p_current_element%next
        call init(section, p_next_element)
        call leaf(traversal, thread, section, thread_traversal%p_current_element)
        thread_traversal%p_current_element => p_next_element

		do
			select case (thread_traversal%p_current_element%cell%geometry%i_edge_types)
				case (INNER_OLD)
                    !init next element for the skeleton operator
                    p_next_element => thread_traversal%p_current_element%next
                    call init(section, p_next_element)
					call old_leaf(traversal, thread, section, thread_traversal%p_current_element)
				case (INNER_NEW)
                    !init next element for the skeleton operator
                    p_next_element => thread_traversal%p_current_element%next
                    call init(section, p_next_element)
					call new_leaf(traversal, thread, section, thread_traversal%p_current_element)
				case (INNER_OLD_BND)
                    !init next element for the skeleton operator
                    p_next_element => thread_traversal%p_current_element%next
                    call init(section, p_next_element)
                    call old_bnd_leaf(traversal, thread, section, thread_traversal%p_current_element)
				case (INNER_NEW_BND)
                    !init next element for the skeleton operator
                    p_next_element => thread_traversal%p_current_element%next
                    call init(section, p_next_element)
                    call new_bnd_leaf(traversal, thread, section, thread_traversal%p_current_element)
                case default
                    !this should happen only for the last element
                    exit
			end select

			thread_traversal%p_current_element => p_next_element
		end do

		!process last element
        call leaf(traversal, thread, section, thread_traversal%p_current_element)
	end if

#	if (_DEBUG_LEVEL > 4)
		_log_write(5, '(2X, A)') "section output state :"
		call section%print()
#	endif

#	if (_DEBUG_LEVEL > 5)
		_log_write(6, '(2X, A)') "output cells:"
		do i = lbound(section%cells%elements, 1), ubound(section%cells%elements, 1)
			_log_write(6, '(3X, I0, X, A)') i, section%cells%elements(i)%to_string()
		end do
		_log_write(6, '(A)') ""
#	endif
end subroutine

subroutine leaf(traversal, thread, section, current_element)
	type(_GT), intent(inout)	                        :: traversal
	type(t_grid_thread), intent(inout)					:: thread
	type(t_grid_section), intent(inout)					:: section
	type(t_traversal_element), intent(inout)	        :: current_element

	call read(traversal, thread, section, current_element)

#	if defined(_GT_ELEMENT_OP)
		call _GT_ELEMENT_OP(traversal, section, current_element%t_element_base)
#	endif

	call write(traversal, thread, section, current_element)
end subroutine

subroutine old_leaf(traversal, thread, section, current_element)
	type(_GT), intent(inout)	                        :: traversal
	type(t_grid_thread), intent(inout)					:: thread
	type(t_grid_section), intent(inout)					:: section
	type(t_traversal_element), intent(inout)	        :: current_element

	call read_oon(traversal, thread, section, current_element)

#	if defined(_GT_INNER_ELEMENT_OP)
		call _GT_INNER_ELEMENT_OP(traversal, section, current_element%t_element_base)
#	endif

	call write_oon(traversal, thread, section, current_element)
end subroutine

subroutine new_leaf(traversal, thread, section, current_element)
	type(_GT), intent(inout)	                        :: traversal
	type(t_grid_thread), intent(inout)					:: thread
	type(t_grid_section), intent(inout)					:: section
	type(t_traversal_element), intent(inout)	        :: current_element

	call read_onn(traversal, thread, section, current_element)

#	if defined(_GT_INNER_ELEMENT_OP)
		call _GT_INNER_ELEMENT_OP(traversal, section, current_element%t_element_base)
#	endif

	call write_onn(traversal, thread, section, current_element)
end subroutine

subroutine old_bnd_leaf(traversal, thread, section, current_element)
	type(_GT), intent(inout)	                        :: traversal
	type(t_grid_thread), intent(inout)					:: thread
	type(t_grid_section), intent(inout)					:: section
	type(t_traversal_element), intent(inout)	        :: current_element

	call read_odn(traversal, thread, section, current_element)

#	if defined(_GT_ELEMENT_OP)
		call _GT_ELEMENT_OP(traversal, section, current_element%t_element_base)
#	endif

	call write_odn(traversal, thread, section, current_element)
end subroutine

subroutine new_bnd_leaf(traversal, thread, section, current_element)
	type(_GT), intent(inout)	                        :: traversal
	type(t_grid_thread), intent(inout)					:: thread
	type(t_grid_section), intent(inout)					:: section
	type(t_traversal_element), intent(inout)	        :: current_element

	call read_obn(traversal, thread, section, current_element)

#	if defined(_GT_ELEMENT_OP)
		call _GT_ELEMENT_OP(traversal, section, current_element%t_element_base)
#	endif

	call write_obn(traversal, thread, section, current_element)
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

#define _GT_INPUT
#define _GT_OUTPUT

#include "Tools_adaptive_traversal.f90"

#undef _GT_INPUT
#undef _GT_OUTPUT

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
#undef _GT_ELEMENT_OP
#undef _GT_INNER_ELEMENT_OP
#undef _GT_EDGE_MERGE_OP
#undef _GT_NODE_MERGE_OP
#undef _GT_PRE_TRAVERSAL_OP
#undef _GT_POST_TRAVERSAL_OP
#undef _GT_PRE_TRAVERSAL_GRID_OP
#undef _GT_POST_TRAVERSAL_GRID_OP

#undef _GT_NODES
#undef _GT_EDGES
#undef _GT_NO_COORDS
#undef _GT_NAME

