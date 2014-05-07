! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


!> Generic traversal template
!> Warning: this a template module body and requires preprocessor commands before inclusion.
!>
!> The resulting method is defined as _GT_NAME
!> @author Oliver Meister

!multiple levels of indirection are necessary to properly resolve the names
#define _GT						_GT_NAME

!if no dedicated inner operators exists, use the default operators
#if defined(_GT_ELEMENT_OP) && !defined(_GT_INNER_ELEMENT_OP)
#	define _GT_INNER_ELEMENT_OP		        _GT_ELEMENT_OP
#endif

#if defined(_GT_EDGE_FIRST_TOUCH_OP) && !defined(_GT_INNER_EDGE_FIRST_TOUCH_OP)
#	define _GT_INNER_EDGE_FIRST_TOUCH_OP	_GT_EDGE_FIRST_TOUCH_OP
#endif

#if defined(_GT_EDGE_LAST_TOUCH_OP) && !defined(_GT_INNER_EDGE_LAST_TOUCH_OP)
#	define _GT_INNER_EDGE_LAST_TOUCH_OP		_GT_EDGE_LAST_TOUCH_OP
#endif

#if defined(_GT_EDGE_REDUCE_OP) && !defined(_GT_INNER_EDGE_REDUCE_OP)
#	define _GT_INNER_EDGE_REDUCE_OP		    _GT_EDGE_REDUCE_OP
#endif

#if defined(_GT_NODE_FIRST_TOUCH_OP) && !defined(_GT_INNER_NODE_FIRST_TOUCH_OP)
#	define _GT_INNER_NODE_FIRST_TOUCH_OP	_GT_NODE_FIRST_TOUCH_OP
#endif

#if defined(_GT_NODE_LAST_TOUCH_OP) && !defined(_GT_INNER_NODE_LAST_TOUCH_OP)
#	define _GT_INNER_NODE_LAST_TOUCH_OP		_GT_NODE_LAST_TOUCH_OP
#endif

#if defined(_GT_NODE_REDUCE_OP) && !defined(_GT_INNER_NODE_REDUCE_OP)
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

type, extends(num_traversal_data) :: t_thread_data
    type(t_traversal_element), dimension(8)     			:: elements									!< Element ring buffer (must contain at least 8 elements, because transfer nodes can be referenced back up to 8 elements)
    type(t_traversal_element), pointer						:: p_current_element => null()				!< Current element
	type(t_statistics)										:: stats
end type

type, extends(num_traversal_data) :: _GT
    type(_GT), pointer                                      :: children(:) => null()			!< section data
    type(t_thread_data), pointer                            :: threads(:) => null()             !< thread data
    type(t_statistics)                                      :: stats
    integer                                                 :: mpi_node_type, mpi_edge_type

    contains

    procedure, pass :: traverse
    procedure, pass :: create
    procedure, pass :: destroy
    procedure, pass :: reduce_stats
end type

contains

subroutine reduce_stats(traversal, mpi_op, global)
    class(_GT)              :: traversal
    integer, intent(in)     :: mpi_op
    logical, intent(in)     :: global

    if (associated(traversal%threads)) then
        call traversal%stats%reduce(traversal%threads(:)%stats, mpi_op)
    end if

    if (global) then
        call traversal%stats%reduce(mpi_op)
    end if
end subroutine

subroutine create(traversal)
    class(_GT)      :: traversal
	integer         :: i_error

#    if defined(_GT_NODE_MPI_TYPE)
        call create_node_mpi_type(traversal%mpi_node_type)
#    endif

#    if defined(_GT_EDGE_MPI_TYPE)
        call create_edge_mpi_type(traversal%mpi_edge_type)
#    endif
end subroutine

subroutine destroy(traversal)
    class(_GT)      :: traversal
	integer         :: i_error

#    if defined(_GT_NODE_MPI_TYPE) && defined(_MPI)
        call MPI_Type_free(traversal%mpi_node_type, i_error); assert_eq(i_error, 0)
#    endif

#    if defined(_GT_EDGE_MPI_TYPE) && defined(_MPI)
        call MPI_Type_free(traversal%mpi_edge_type, i_error); assert_eq(i_error, 0)
#    endif

    if (associated(traversal%children)) then
        deallocate(traversal%children, stat = i_error); assert_eq(i_error, 0)
    end if

    if (associated(traversal%threads)) then
        deallocate(traversal%threads, stat = i_error); assert_eq(i_error, 0)
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

function edge_write_wrapper_op(local_edges, neighbor_edges) result(l_conform)
    type(t_edge_data), intent(inout)    :: local_edges
    type(t_edge_data), intent(in)       :: neighbor_edges
    logical                             :: l_conform

    assert_eq(local_edges%min_distance, neighbor_edges%min_distance)
    assert(.not. local_edges%owned_locally)
    assert(neighbor_edges%owned_locally)

#   if defined(_GT_EDGE_WRITE_OP)
        call _GT_EDGE_WRITE_OP(local_edges, neighbor_edges)
#   else
        local_edges%data_pers = neighbor_edges%data_pers
        local_edges%data_temp = neighbor_edges%data_temp
#   endif

    l_conform = .true.
end function

function node_write_wrapper_op(local_nodes, neighbor_nodes) result(l_conform)
    type(t_node_data), intent(inout)    :: local_nodes
    type(t_node_data), intent(in)       :: neighbor_nodes
    logical                             :: l_conform

    assert_eq(local_nodes%distance, neighbor_nodes%distance)
    assert_veq(local_nodes%position, neighbor_nodes%position)
    assert(.not. local_nodes%owned_locally)
    assert(neighbor_nodes%owned_locally)

#   if defined(_GT_NODE_WRITE_OP)
        call _GT_NODE_WRITE_OP(local_nodes, neighbor_nodes)
#   else
        local_nodes%data_pers = neighbor_nodes%data_pers
        local_nodes%data_temp = neighbor_nodes%data_temp
#   endif

    l_conform = .true.
end function

!*****************************************************************
! Generic traversal
!*****************************************************************

!> Generic iterative traversal subroutine
!> @author Oliver Meister
subroutine traverse(traversal, grid)
	class(_GT), intent(inout)	                        :: traversal
	type(t_grid), intent(inout)					        :: grid

	integer (kind = GRID_SI)                            :: i_section, i_first_local_section, i_last_local_section
	integer (kind = GRID_SI)                            :: i_error
	type(t_statistics)                                  :: thread_stats

	integer (kind = GRID_SI), save                      :: i_thread
	!$omp threadprivate(i_thread)

    if (.not. associated(traversal%children) .or. size(traversal%children) .ne. grid%sections%get_size()) then
        !$omp barrier

        !$omp single
        if (associated(traversal%children)) then
            deallocate(traversal%children, stat = i_error); assert_eq(i_error, 0)
        end if

        allocate(traversal%children(grid%sections%get_size()), stat = i_error); assert_eq(i_error, 0)
        !$omp end single
    end if

    if (.not. associated(traversal%threads)) then
        !$omp barrier

        !$omp single
        allocate(traversal%threads(omp_get_max_threads()), stat = i_error); assert_eq(i_error, 0)
        !$omp end single

        i_thread = 1 + omp_get_thread_num()
        call create_ringbuffer(traversal%threads(i_thread)%elements)
        traversal%threads(i_thread)%p_current_element => traversal%threads(i_thread)%elements(1)
    end if

    assert_eq(size(traversal%threads), omp_get_max_threads())

    call grid%get_local_sections(i_first_local_section, i_last_local_section)

    thread_stats%r_traversal_time = -get_wtime()

#   if defined(_ASAGI_TIMING)
        grid%sections%elements_alloc(i_first_local_section : i_last_local_section)%stats%r_asagi_time = 0.0
#   endif

    thread_stats%r_barrier_time = -get_wtime()

    select type (traversal)
        type is (_GT)
            !$omp single
            call pre_traversal_grid(traversal, grid)
            !$omp end single nowait
        class default
            assert(.false.)
    end select

    thread_stats%r_barrier_time = thread_stats%r_barrier_time + get_wtime()

    thread_stats%r_computation_time = 0.0

#   if defined(_GT_SKELETON_OP)
        thread_stats%r_computation_time = -get_wtime()

        do i_section = i_first_local_section, i_last_local_section
            assert_eq(i_section, grid%sections%elements_alloc(i_section)%index)

            call boundary_skeleton(traversal%children(i_section), grid%sections%elements_alloc(i_section))
        end do

        thread_stats%r_computation_time = thread_stats%r_computation_time + get_wtime()
#   endif

    !$omp barrier

    !thread_stats%r_computation_time = thread_stats%r_computation_time - get_wtime()

    do i_section = i_first_local_section, i_last_local_section
#       if !defined(_GT_NODE_MPI_TYPE) && !defined(_GT_EDGE_MPI_TYPE)
            call recv_mpi_boundary(grid%sections%elements_alloc(i_section))
#       elif defined(_GT_NODE_MPI_TYPE)
            call recv_mpi_boundary(grid%sections%elements_alloc(i_section), mpi_node_type_optional=traversal%mpi_node_type)
#       elif defined(_GT_EDGE_MPI_TYPE)
            call recv_mpi_boundary(grid%sections%elements_alloc(i_section), mpi_edge_type_optional=traversal%mpi_edge_type)
#       else
            call recv_mpi_boundary(grid%sections%elements_alloc(i_section), mpi_edge_type_optional=traversal%mpi_edge_type, mpi_node_type_optional=traversal%mpi_node_type)
#       endif
    end do

    do i_section = i_first_local_section, i_last_local_section
#       if defined(_OPENMP_TASKS)
            !$omp task default(shared) firstprivate(i_section) mergeable
#       endif

        traversal%children(i_section)%stats%r_computation_time = -get_wtime()

        call pre_traversal(traversal%children(i_section), grid%sections%elements_alloc(i_section))
        call traverse_section(traversal%threads(i_thread), traversal%children(i_section), grid%threads%elements(i_thread), grid%sections%elements_alloc(i_section))

#       if !defined(_GT_NODE_MPI_TYPE) && !defined(_GT_EDGE_MPI_TYPE)
            call send_mpi_boundary(grid%sections%elements_alloc(i_section))
#       elif defined(_GT_NODE_MPI_TYPE)
            call send_mpi_boundary(grid%sections%elements_alloc(i_section), mpi_node_type_optional=traversal%mpi_node_type)
#       elif defined(_GT_EDGE_MPI_TYPE)
            call send_mpi_boundary(grid%sections%elements_alloc(i_section), mpi_edge_type_optional=traversal%mpi_edge_type)
#       else
            call send_mpi_boundary(grid%sections%elements_alloc(i_section), mpi_edge_type_optional=traversal%mpi_edge_type, mpi_node_type_optional=traversal%mpi_node_type)
#       endif

        traversal%children(i_section)%stats%r_computation_time = traversal%children(i_section)%stats%r_computation_time + get_wtime()

#       if defined(_OPENMP_TASKS)
            !$omp end task
#       endif
    end do

#   if defined(_OPENMP_TASKS)
        !$omp taskwait
#   endif

    !thread_stats%r_computation_time = thread_stats%r_computation_time + get_wtime()

    !sync and call post traversal operator
    thread_stats%r_sync_time = -get_wtime()
    call sync_boundary(grid, edge_merge_wrapper_op, node_merge_wrapper_op, edge_write_wrapper_op, node_write_wrapper_op)
    thread_stats%r_sync_time = thread_stats%r_sync_time + get_wtime()

    !call post traversal operator
    thread_stats%r_computation_time = thread_stats%r_computation_time - get_wtime()

    do i_section = i_first_local_section, i_last_local_section
        assert_eq(i_section, grid%sections%elements_alloc(i_section)%index)
        call post_traversal(traversal%children(i_section), grid%sections%elements_alloc(i_section))
    end do

    thread_stats%r_computation_time = thread_stats%r_computation_time + get_wtime()

    call grid%reverse()

    !$omp barrier

    thread_stats%r_barrier_time = thread_stats%r_barrier_time - get_wtime()

    select type (traversal)
        type is (_GT)
            !$omp single
            call post_traversal_grid(traversal, grid)
            !$omp end single
        class default
            assert(.false.)
    end select

    thread_stats%r_barrier_time = thread_stats%r_barrier_time + get_wtime()
    thread_stats%r_traversal_time = thread_stats%r_traversal_time + get_wtime()

#   if defined(_ASAGI_TIMING)
        !HACK: in lack of a better method, we reduce ASAGI timing data like this for now - should be changed in the long run, so that stats belongs to the section and not the traversal
        thread_stats%r_asagi_time = sum(grid%sections%elements_alloc(i_first_local_section : i_last_local_section)%stats%r_asagi_time)
#   endif

    do i_section = i_first_local_section, i_last_local_section
        call set_stats_counters(traversal%children(i_section)%stats, grid%sections%elements_alloc(i_section))
        grid%sections%elements_alloc(i_section)%stats%t_statistics = grid%sections%elements_alloc(i_section)%stats%t_statistics + traversal%children(i_section)%stats
        thread_stats = thread_stats + traversal%children(i_section)%stats
    end do

    traversal%threads(i_thread)%stats = traversal%threads(i_thread)%stats + thread_stats
    grid%threads%elements(i_thread)%stats%t_statistics = grid%threads%elements(i_thread)%stats%t_statistics + thread_stats
end subroutine

!> Generic iterative traversal subroutine
!> @author Oliver Meister
subroutine traverse_section(thread_traversal, traversal, thread, section)
	type(t_thread_data), target, intent(inout)	        :: thread_traversal
	type(_GT), target, intent(inout)	                :: traversal
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

	if (section%cells%get_size() > 0) then
        assert_ge(section%cells%get_size(), 2)
		call init(section, thread_traversal%p_current_element)

		!process first element
        p_next_element => thread_traversal%p_current_element%next
        call init(section, p_next_element)
        call leaf(traversal, thread, section, thread_traversal%p_current_element)
        thread_traversal%p_current_element => p_next_element

		do
			select case (thread_traversal%p_current_element%cell%geometry%i_entity_types)
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
#undef _GT_EDGE_MPI_TYPE
#undef _GT_EDGE_MERGE_OP
#undef _GT_EDGE_WRITE_OP
#undef _GT_NODE_FIRST_TOUCH_OP
#undef _GT_NODE_LAST_TOUCH_OP
#undef _GT_NODE_REDUCE_OP
#undef _GT_NODE_MPI_TYPE
#undef _GT_NODE_MERGE_OP
#undef _GT_NODE_WRITE_OP
#undef _GT_INNER_EDGE_FIRST_TOUCH_OP
#undef _GT_INNER_EDGE_LAST_TOUCH_OP
#undef _GT_INNER_EDGE_REDUCE_OP
#undef _GT_INNER_NODE_FIRST_TOUCH_OP
#undef _GT_INNER_NODE_LAST_TOUCH_OP
#undef _GT_INNER_NODE_REDUCE_OP
#undef _GT_ELEMENT_OP
#undef _GT_INNER_ELEMENT_OP
#undef _GT_PRE_TRAVERSAL_OP
#undef _GT_POST_TRAVERSAL_OP
#undef _GT_PRE_TRAVERSAL_GRID_OP
#undef _GT_POST_TRAVERSAL_GRID_OP

#undef _GT_NODES
#undef _GT_EDGES
#undef _GT_NO_COORDS
#undef _GT_NAME
