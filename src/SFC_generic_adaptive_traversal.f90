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

#if defined(_GT_INNER_EDGE_FIRST_TOUCH_OP) || defined(_GT_INNER_NODE_FIRST_TOUCH_OP)
#	error "No inner first touch operators are allowed for adaptive traversal!"
#endif

#if defined(_GT_EDGE_FIRST_TOUCH_OP)
#	define _GT_INNER_EDGE_FIRST_TOUCH_OP	_GT_EDGE_FIRST_TOUCH_OP
#endif

#if defined(_GT_NODE_FIRST_TOUCH_OP)
#	define _GT_INNER_NODE_FIRST_TOUCH_OP	_GT_NODE_FIRST_TOUCH_OP
#endif

!if no dedicated inner operators exists, use the default operators
#if defined(_GT_EDGE_LAST_TOUCH_OP) && !defined(_GT_INNER_EDGE_LAST_TOUCH_OP)
#	define _GT_INNER_EDGE_LAST_TOUCH_OP		_GT_EDGE_LAST_TOUCH_OP
#endif

#if defined(_GT_NODE_LAST_TOUCH_OP) && !defined(_GT_INNER_NODE_LAST_TOUCH_OP)
#	define _GT_INNER_NODE_LAST_TOUCH_OP		_GT_NODE_LAST_TOUCH_OP
#endif

#if defined(_GT_EDGE_REDUCE_OP) && !defined(_GT_INNER_EDGE_REDUCE_OP)
#	define _GT_INNER_EDGE_REDUCE_OP		    _GT_EDGE_REDUCE_OP
#endif

#if defined(_GT_NODE_REDUCE_OP) && !defined(_GT_INNER_NODE_REDUCE_OP)
#	define _GT_INNER_NODE_REDUCE_OP		    _GT_NODE_REDUCE_OP
#endif

#define _GT_EDGES

PRIVATE
PUBLIC _GT

!Module types:

!> Traversal element ring buffer structure that provides local storage for some of the grid data
type, extends(t_element_base) :: t_traversal_element
	integer (KIND = GRID_SI)                            :: i_cell
	real(kind = GRID_SR)				                :: coords(2,3)            				    !< coordinates
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

type, extends(num_traversal_data) :: t_thread_data
	type(t_refinement_element), dimension(2)				:: refinement_elements					!< Temporary refinement elements

    type(t_traversal_element), dimension(8)			        :: src_elements							!< Input element ring buffer (must be 8, because transfer nodes can be referenced back up to 8 elements)
    type(t_traversal_element), dimension(11)        		:: dest_elements						!< Output element ring buffer (must be 11, because transfer nodes can be referenced back up to 8 + 3 new refinement elements)

    type(t_traversal_element), pointer						:: p_src_element => null(), p_dest_element => null()        !< Current source and destination element

	type(t_adaptive_statistics)								:: stats
end type

type, extends(num_traversal_data) :: _GT
    type(_GT), pointer                                      :: children(:) => null()			!< section data
    type(t_thread_data), pointer                            :: threads(:) => null()             !< thread data
	type(t_adaptive_statistics)								:: stats
    type(t_conformity)                                      :: conformity
    integer                                                 :: mpi_node_type, mpi_edge_type

    contains

    procedure, pass :: create
    procedure, pass :: destroy
    procedure, pass :: traverse => traverse_in_place
    procedure, pass :: reduce_stats
    procedure, pass :: clear_stats
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

    call traversal%conformity%reduce_stats(mpi_op, global)
end subroutine

subroutine clear_stats(traversal)
    class(_GT)              :: traversal
    integer                 :: i

    if (associated(traversal%threads)) then
        do i = 1, size(traversal%threads)
            call traversal%threads(i)%stats%clear()
        end do
    end if

    call traversal%stats%clear()

    call traversal%conformity%clear_stats()
end subroutine

subroutine create(traversal)
    class(_GT)      :: traversal
	integer         :: i_error

	call traversal%conformity%create()

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

	call traversal%conformity%destroy()

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
! Generic adaptive traversal
!*****************************************************************

!> Generic adaptive traversal subroutine
!> @author Oliver Meister
subroutine traverse_in_place(traversal, grid)
	class(_GT), target, intent(inout)	                :: traversal
	type(t_grid), intent(inout)							:: grid

    integer (kind = GRID_SI)                            :: i_src_section, i_first_src_section, i_last_src_section
	integer                                             :: i_error
	type(t_grid), save							        :: grid_temp
	type(t_adaptive_statistics)                         :: thread_stats

    if (.not. associated(traversal%threads) .or. size(traversal%threads) .ne. cfg%i_threads) then
        !$omp barrier

        !$omp single
        allocate(traversal%threads(cfg%i_threads), stat = i_error); assert_eq(i_error, 0)
        !$omp end single

    	assert(.not. associated(traversal%threads(i_thread)%p_dest_element))

    	!create ringbuffers and temporary elements
		call create_ringbuffer(traversal%threads(i_thread)%src_elements)
		call create_ringbuffer(traversal%threads(i_thread)%dest_elements)
		call create_refinement_element(traversal%threads(i_thread)%refinement_elements)

		traversal%threads(i_thread)%p_src_element => traversal%threads(i_thread)%src_elements(1)
		traversal%threads(i_thread)%p_dest_element => traversal%threads(i_thread)%dest_elements(1)
    end if

    assert_eq(size(traversal%threads), omp_get_max_threads())

    thread_stats = t_adaptive_statistics()

    !$omp barrier

    thread_stats%r_integrity_time = -get_wtime()
    call traversal%conformity%check(grid)
    thread_stats%r_integrity_time = thread_stats%r_integrity_time + get_wtime()

    !$omp barrier

    call grid%get_local_sections(i_first_src_section, i_last_src_section)

    do i_src_section = i_first_src_section, i_last_src_section
        call grid%sections%elements_alloc(i_src_section)%estimate_load()
    end do

    !duplicate source grid boundary before load balancing
    thread_stats%r_sync_time = thread_stats%r_sync_time - get_wtime()
    call duplicate_boundary_data(grid, edge_write_wrapper_op, node_write_wrapper_op)
    thread_stats%r_sync_time = thread_stats%r_sync_time + get_wtime()

    !wait until all boundary data has been copied
    !$omp barrier

    !balance load BEFORE refinement if splitting is DISABLED
    if (.not. cfg%l_split_sections) then
#	    if !defined(_GT_INPUT_DEST)
	        !exchange source grid sections with neighbors
            thread_stats%r_load_balancing_time = -get_wtime()
		    call distribute_load(grid, 0.0)
            thread_stats%r_load_balancing_time = thread_stats%r_load_balancing_time + get_wtime()
#	    endif
    end if

    thread_stats%r_allocation_time = -get_wtime()
    call create_destination_grid(grid, grid_temp)

    !if necessary, reverse order to match source and destination grid
    if (.not. grid%sections%is_forward()) then
        !$omp barrier
        _log_write(4, '(X, A)') "Reverse destination grid.."
        call grid_temp%reverse()
        !$omp barrier
    end if

    assert_eqv(grid%sections%is_forward(), grid_temp%sections%is_forward())
    thread_stats%r_allocation_time = thread_stats%r_allocation_time + get_wtime()

    !refine grid
    call traverse_out_of_place(traversal, grid, grid_temp)

    !$omp barrier

    call grid%destroy()

    !$omp single
    call grid_temp%move(grid)
    !$omp end single

    !balance load AFTER refinement if splitting is ENABLED
    if (cfg%l_split_sections) then
#       if !defined(_GT_INPUT_DEST)
            !duplicate destination grid boundary before load balancing
            thread_stats%r_sync_time = thread_stats%r_sync_time - get_wtime()
            call duplicate_boundary_data(grid, edge_write_wrapper_op, node_write_wrapper_op)
            thread_stats%r_sync_time = thread_stats%r_sync_time + get_wtime()

            !wait until all boundary data has been copied
            !$omp barrier

	        !exchange destination grid sections with neighbors
            thread_stats%r_load_balancing_time = -get_wtime()
		    call distribute_load(grid, 0.0)
            thread_stats%r_load_balancing_time = thread_stats%r_load_balancing_time + get_wtime()
#	    endif
    end if

    traversal%threads(i_thread)%stats = traversal%threads(i_thread)%stats + thread_stats
    grid%threads%elements(i_thread)%stats = grid%threads%elements(i_thread)%stats + thread_stats
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

	type(t_grid_section), pointer							:: dest_section
    integer (kind = GRID_SI)			                    :: i_dest_section, i_first_local_section, i_last_local_section
	integer                                                 :: i_error
	integer (kind = BYTE)								    :: i_color
	type(t_adaptive_statistics)                             :: thread_stats

    integer (kind = GRID_SI), save			                :: i_src_section, i_src_cell
	type(t_grid_section), target, save					    :: src_section
	!$omp threadprivate(i_src_section, i_src_cell, src_section)

    if (.not. associated(traversal%children) .or. size(traversal%children) .ne. dest_grid%sections%get_size()) then
		!$omp barrier

		!$omp single
        if (associated(traversal%children)) then
            deallocate(traversal%children, stat = i_error); assert_eq(i_error, 0)
        end if

        allocate(traversal%children(dest_grid%sections%get_size()), stat = i_error); assert_eq(i_error, 0)
    	!$omp end single
    end if

    call dest_grid%get_local_sections_in_traversal_order(i_first_local_section, i_last_local_section)

    thread_stats%r_traversal_time = -get_wtime()

#   if defined(_ASAGI_TIMING)
        dest_grid%sections%elements_alloc(i_first_local_section : i_last_local_section)%stats%r_asagi_time = 0.0
        thread_stats%r_asagi_time = 0.0
#   endif

    thread_stats%r_barrier_time = -get_wtime()
    select type (traversal)
        type is (_GT)
            !$omp single
            call pre_traversal_grid_dest(traversal, dest_grid)

            call prefix_sum(src_grid%sections%elements%last_dest_cell, src_grid%sections%elements%dest_cells)
            call prefix_sum(dest_grid%sections%elements%last_dest_cell, dest_grid%sections%elements%dest_cells)
            !$omp end single
        class default
            assert(.false.)
    end select
    thread_stats%r_barrier_time = thread_stats%r_barrier_time + get_wtime()

    thread_stats%r_computation_time = thread_stats%r_computation_time - get_wtime()

    !call pre traversal operator
    thread_stats%r_pre_compute_time = -get_wtime()
    do i_dest_section = i_first_local_section, i_last_local_section
        call pre_traversal_dest(traversal%children(i_dest_section), dest_grid%sections%elements(i_dest_section))
    end do
    thread_stats%r_pre_compute_time = thread_stats%r_pre_compute_time + get_wtime()

    thread_stats%r_computation_time = thread_stats%r_computation_time + get_wtime()

    !traversal

    if (src_grid%sections%get_size() > 0) then
        !make an initial guess for the first source section and the first source cell
        i_src_section = max(1, min(src_grid%sections%get_size(), &
            1 + ((i_first_local_section - 1) * src_grid%sections%get_size()) / dest_grid%sections%get_size()))
        i_src_cell = 0
    end if

    !traverse all destination sections
    do i_dest_section = i_first_local_section, i_last_local_section
#       if defined(_OPENMP_TASKS_ADAPTIVITY)
            !$omp task default(shared) firstprivate(i_dest_section) private(dest_section) mergeable
#       endif

        thread_stats%r_computation_time = thread_stats%r_computation_time - get_wtime()

        dest_section => dest_grid%sections%elements(i_dest_section)

        if (i_src_cell .ne. dest_section%last_dest_cell - dest_section%dest_cells + 1) then
            !find the source section that contains the first cell of the first destination section
            if (i_src_cell .eq. 0 .or. &
                i_src_cell .ge. dest_section%last_dest_cell - dest_section%dest_cells + 5 .or. &
                src_grid%sections%elements(i_src_section)%last_dest_cell .le. dest_section%last_dest_cell - dest_section%dest_cells) then

                do while (src_grid%sections%elements(i_src_section)%last_dest_cell - src_grid%sections%elements(i_src_section)%dest_cells .gt. dest_section%last_dest_cell - dest_section%dest_cells)
                    i_src_section = i_src_section - 1
                    assert_ge(i_src_section, 1)
                end do

                do while (src_grid%sections%elements(i_src_section)%last_dest_cell .le. dest_section%last_dest_cell - dest_section%dest_cells)
                    i_src_section = i_src_section + 1
                    assert_le(i_src_section, src_grid%sections%get_size())
                end do

                src_section = src_grid%sections%elements(i_src_section)
                i_src_cell = src_section%last_dest_cell - src_section%dest_cells + 1
            end if

            !traverse all unnecessary elements of the source section with an empty traversal
            do while (i_src_cell .le. dest_section%last_dest_cell - dest_section%dest_cells)
                i_src_cell = i_src_cell + empty_leaf(traversal%threads(i_thread), traversal%children(i_dest_section), src_grid%threads%elements(i_thread), src_section)
            end do
        end if

#       if _DEBUG_LEVEL > 4
            _log_write(5, '(2X, A)') "destination section initial state :"
            call dest_section%print()
#       endif

        !traverse all source sections that overlap with the destination section
        do while (i_src_cell .le. dest_section%last_dest_cell)
            assert_le(i_src_section, src_grid%sections%get_size())
            assert_ge(i_src_cell, src_section%last_dest_cell - src_section%dest_cells + 1)

            !traverse all elements
            do while (i_src_cell .le. min(src_section%last_dest_cell, dest_section%last_dest_cell))
                i_src_cell = i_src_cell + leaf(traversal%threads(i_thread), traversal%children(i_dest_section), src_grid%threads%elements(i_thread), dest_grid%threads%elements(i_thread), src_section, dest_section)
            end do

            if (i_src_cell .gt. src_section%last_dest_cell .and. i_src_section < src_grid%sections%get_size()) then
                i_src_section = i_src_section + 1
                src_section = src_grid%sections%elements(i_src_section)
            end if
        end do

#	    if defined(_GT_CELL_TO_EDGE_OP)
            traversal%threads(i_thread)%p_dest_element%previous%next_edge_data%rep = _GT_CELL_TO_EDGE_OP(traversal%threads(i_thread)%p_dest_element%previous%t_element_base, traversal%threads(i_thread)%p_dest_element%previous%next_edge_data)
#       endif

        traversal%threads(i_thread)%p_dest_element%previous%next_edge_data%remove = .false.

        call find_section_boundary_elements(dest_grid%threads%elements(i_thread), dest_section, dest_section%cells%i_current_element, traversal%threads(i_thread)%p_dest_element%previous%next_edge_data)

        thread_stats%r_computation_time = thread_stats%r_computation_time + get_wtime()

#       if defined(_OPENMP_TASKS_ADAPTIVITY)
            !$omp end task
#       endif
    end do

#   if defined(_OPENMP_TASKS_ADAPTIVITY)
        !$omp taskwait
#   endif

    thread_stats%r_update_distances_time = -get_wtime()
    call update_distances(dest_grid)
    thread_stats%r_update_distances_time = thread_stats%r_update_distances_time + get_wtime()

    !$omp single
    assert_veq(decode_distance(dest_grid%t_global_data%min_distance), decode_distance(src_grid%t_global_data%min_distance))
    !$omp end single nowait

    !update communication info
    thread_stats%r_update_neighbors_time = -get_wtime()
    call update_neighbors(src_grid, dest_grid)
    thread_stats%r_update_neighbors_time = thread_stats%r_update_neighbors_time + get_wtime()

    !$omp barrier

    thread_stats%r_sync_time = thread_stats%r_sync_time - get_wtime()
    do i_dest_section = i_first_local_section, i_last_local_section
#       if !defined(_GT_NODE_MPI_TYPE) && !defined(_GT_EDGE_MPI_TYPE)
            call recv_mpi_boundary(dest_grid%sections%elements(i_dest_section))
            call send_mpi_boundary(dest_grid%sections%elements(i_dest_section))
#       elif defined(_GT_NODE_MPI_TYPE) && !defined(_GT_EDGE_MPI_TYPE)
            call recv_mpi_boundary(dest_grid%sections%elements(i_dest_section), mpi_node_type_optional=traversal%mpi_node_type)
            call send_mpi_boundary(dest_grid%sections%elements(i_dest_section), mpi_node_type_optional=traversal%mpi_node_type)
#       elif defined(_GT_EDGE_MPI_TYPE) && !defined(_GT_NODE_MPI_TYPE)
            call recv_mpi_boundary(dest_grid%sections%elements(i_dest_section), mpi_edge_type_optional=traversal%mpi_edge_type)
            call send_mpi_boundary(dest_grid%sections%elements(i_dest_section), mpi_edge_type_optional=traversal%mpi_edge_type)
#       else
            call recv_mpi_boundary(dest_grid%sections%elements(i_dest_section), mpi_edge_type_optional=traversal%mpi_edge_type, mpi_node_type_optional=traversal%mpi_node_type)
            call send_mpi_boundary(dest_grid%sections%elements(i_dest_section), mpi_edge_type_optional=traversal%mpi_edge_type, mpi_node_type_optional=traversal%mpi_node_type)
#       endif
    end do

    !sync and call post traversal operator
#   if !defined(_GT_NODE_MPI_TYPE) && !defined(_GT_EDGE_MPI_TYPE)
        call collect_boundary_data(dest_grid, edge_merge_wrapper_op, node_merge_wrapper_op)
#   elif defined(_GT_NODE_MPI_TYPE) && !defined(_GT_EDGE_MPI_TYPE)
        call collect_boundary_data(dest_grid, edge_merge_wrapper_op, node_merge_wrapper_op, mpi_node_type_optional=traversal%mpi_node_type)
#   elif defined(_GT_EDGE_MPI_TYPE) && !defined(_GT_NODE_MPI_TYPE)
        call collect_boundary_data(dest_grid, edge_merge_wrapper_op, node_merge_wrapper_op, mpi_edge_type_optional=traversal%mpi_edge_type)
#   else
        call collect_boundary_data(dest_grid, edge_merge_wrapper_op, node_merge_wrapper_op, mpi_node_type_optional=traversal%mpi_node_type, mpi_edge_type_optional=traversal%mpi_edge_type)
#   endif
    thread_stats%r_sync_time = thread_stats%r_sync_time + get_wtime()

    thread_stats%r_post_compute_time = -get_wtime()
    do i_dest_section = i_first_local_section, i_last_local_section
        call post_traversal_dest(traversal%children(i_dest_section), dest_grid%sections%elements(i_dest_section))
    end do
    thread_stats%r_post_compute_time = thread_stats%r_post_compute_time + get_wtime()

    thread_stats%r_computation_time = thread_stats%r_computation_time + thread_stats%r_post_compute_time

    !$omp barrier

    thread_stats%r_barrier_time = thread_stats%r_barrier_time - get_wtime()

    select type (traversal)
        type is (_GT)
            !$omp single
            call post_traversal_grid_dest(traversal, dest_grid)
            !$omp end single
        class default
            assert(.false.)
    end select

    thread_stats%r_barrier_time = thread_stats%r_barrier_time + get_wtime()
    thread_stats%r_traversal_time = thread_stats%r_traversal_time + get_wtime()

    !HACK: in lack of a better method, we reduce ASAGI timing data like this for now - should be changed in the long run, so that stats belongs to the section and not the traversal
#   if defined(_ASAGI_TIMING)
        thread_stats%r_asagi_time = sum(dest_grid%sections%elements_alloc(i_first_local_section : i_last_local_section)%stats%r_asagi_time)
#   endif

    do i_dest_section = i_first_local_section, i_last_local_section
        call set_stats_counters_dest(traversal%children(i_dest_section)%stats%t_statistics, dest_grid%sections%elements_alloc(i_dest_section))
        dest_grid%sections%elements_alloc(i_dest_section)%stats = dest_grid%sections%elements_alloc(i_dest_section)%stats + traversal%children(i_dest_section)%stats
        thread_stats = thread_stats + traversal%children(i_dest_section)%stats
    end do

    traversal%threads(i_thread)%stats = traversal%threads(i_thread)%stats + thread_stats
    dest_grid%threads%elements(i_thread)%stats = dest_grid%threads%elements(i_thread)%stats + thread_stats

#   if defined(_GT_OUTPUT_SRC)
        call src_grid%reverse()
#   endif

	call dest_grid%reverse()

    !$omp barrier
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

#			if !defined(_GT_INPUT_DEST)
				call create_parent_cell(thread_traversal%p_src_element%cell%geometry, thread_traversal%p_src_element%next%cell%geometry, thread_traversal%p_dest_element%cell%geometry)
#			endif

			call read_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

#			if !defined(_GT_INPUT_DEST)
				thread_traversal%p_dest_element%coords(:,1) = thread_traversal%p_src_element%nodes(1)%ptr%position
				thread_traversal%p_dest_element%coords(:,2) = thread_traversal%p_src_element%nodes(3)%ptr%position
				thread_traversal%p_dest_element%coords(:,3) = thread_traversal%p_src_element%nodes(2)%ptr%position + (thread_traversal%p_src_element%nodes(2)%ptr%position - thread_traversal%p_src_element%nodes(1)%ptr%position)

                assert_veq(thread_traversal%p_dest_element%coords, thread_traversal%p_dest_element%coords)
                assert_veq(thread_traversal%p_dest_element%coords, thread_traversal%p_dest_element%coords)
                assert_veq(thread_traversal%p_dest_element%coords, thread_traversal%p_dest_element%coords)
#			endif

			call read_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

			call _GT_COARSEN_OP(traversal, dest_section, thread_traversal%p_src_element, thread_traversal%p_dest_element, [1])
			call write_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

			thread_traversal%p_src_element => thread_traversal%p_src_element%next

			call read_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

			call _GT_COARSEN_OP(traversal, dest_section, thread_traversal%p_src_element, thread_traversal%p_dest_element, [2])

			call write_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

			call write_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

			thread_traversal%p_src_element => thread_traversal%p_src_element%next
			thread_traversal%p_dest_element => thread_traversal%p_dest_element%next

			i_dest_cells = 1
		case (0)
			!transfer
#			if !defined(_GT_INPUT_DEST)
				call copy_cell(thread_traversal%p_src_element%cell%geometry, thread_traversal%p_dest_element%cell%geometry)
#			endif

			call read_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

#			if !defined(_GT_INPUT_DEST)
                thread_traversal%p_dest_element%coords(:,1) = thread_traversal%p_src_element%nodes(1)%ptr%position
                thread_traversal%p_dest_element%coords(:,2) = thread_traversal%p_src_element%nodes(2)%ptr%position
                thread_traversal%p_dest_element%coords(:,3) = thread_traversal%p_src_element%nodes(3)%ptr%position
                assert_veq(thread_traversal%p_dest_element%coords, thread_traversal%p_dest_element%coords)
#			endif

			call read_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

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

#			if !defined(_GT_INPUT_DEST)
				call create_child_cells(thread_traversal%p_src_element%cell%geometry, thread_traversal%p_dest_element%cell%geometry, p_dest_element_2%cell%geometry)
#			endif

			call read_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

#			if !defined(_GT_INPUT_DEST)
				thread_traversal%p_dest_element%coords(:,1) = thread_traversal%p_src_element%nodes(1)%ptr%position
				thread_traversal%p_dest_element%coords(:,2) = 0.5_GRID_SR * (thread_traversal%p_src_element%nodes(1)%ptr%position + thread_traversal%p_src_element%nodes(3)%ptr%position)
				thread_traversal%p_dest_element%coords(:,3) = thread_traversal%p_src_element%nodes(2)%ptr%position
                assert_veq(thread_traversal%p_dest_element%coords, thread_traversal%p_dest_element%coords)
#			endif

			call read_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

			call _GT_REFINE_OP(traversal, dest_section, thread_traversal%p_src_element, thread_traversal%p_dest_element, [1])
			call write_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

#			if !defined(_GT_INPUT_DEST)
                p_dest_element_2%coords(:,1) = thread_traversal%p_dest_element%coords(:,3)
                p_dest_element_2%coords(:,2) = thread_traversal%p_dest_element%coords(:,2)
				p_dest_element_2%coords(:,3) = thread_traversal%p_src_element%nodes(3)%ptr%position
                assert_veq(p_dest_element_2%coords, p_dest_element_2%coords)
#			endif

			call read_dest_element(traversal, dest_thread, dest_section, p_dest_element_2)

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

#			if !defined(_GT_INPUT_DEST)
				call create_child_cells(thread_traversal%p_src_element%cell%geometry, thread_traversal%refinement_elements(1)%cell%geometry, p_dest_element_3%cell%geometry)
				call create_child_cells(thread_traversal%refinement_elements(1)%cell%geometry, thread_traversal%p_dest_element%cell%geometry, p_dest_element_2%cell%geometry)
#			endif

			call read_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

#			if !defined(_GT_INPUT_DEST)
				thread_traversal%p_dest_element%coords(:,1) = thread_traversal%p_src_element%nodes(1)%ptr%position
				thread_traversal%p_dest_element%coords(:,2) = 0.5_GRID_SR * (thread_traversal%p_src_element%nodes(1)%ptr%position + thread_traversal%p_src_element%nodes(2)%ptr%position)
				thread_traversal%p_dest_element%coords(:,3) = 0.5_GRID_SR * (thread_traversal%p_src_element%nodes(1)%ptr%position + thread_traversal%p_src_element%nodes(3)%ptr%position)
                assert_veq(thread_traversal%p_dest_element%coords, thread_traversal%p_dest_element%coords)
#			endif

			call read_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

			call _GT_REFINE_OP(traversal, dest_section, thread_traversal%p_src_element, thread_traversal%p_dest_element, [1, 1])
			call write_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

#			if !defined(_GT_INPUT_DEST)
                p_dest_element_2%coords(:,1) = thread_traversal%p_dest_element%coords(:,3)
                p_dest_element_2%coords(:,2) = thread_traversal%p_dest_element%coords(:,2)
				p_dest_element_2%coords(:,3) = thread_traversal%p_src_element%nodes(2)%ptr%position
                assert_veq(p_dest_element_2%coords, p_dest_element_2%coords)
#			endif

			call read_dest_element(traversal, dest_thread, dest_section, p_dest_element_2)

			call _GT_REFINE_OP(traversal, dest_section, thread_traversal%p_src_element, p_dest_element_2, [1, 2])
			call write_dest_element(traversal, dest_thread, dest_section, p_dest_element_2)

#			if !defined(_GT_INPUT_DEST)
                p_dest_element_3%coords(:,1) = p_dest_element_2%coords(:,3)
                p_dest_element_3%coords(:,2) = p_dest_element_2%coords(:,1)
				p_dest_element_3%coords(:,3) = thread_traversal%p_src_element%nodes(3)%ptr%position
                assert_veq(p_dest_element_3%coords, p_dest_element_3%coords)
#			endif

			call read_dest_element(traversal, dest_thread, dest_section, p_dest_element_3)

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

#			if !defined(_GT_INPUT_DEST)
				call create_child_cells(thread_traversal%p_src_element%cell%geometry, thread_traversal%p_dest_element%cell%geometry, thread_traversal%refinement_elements(1)%cell%geometry)
				call create_child_cells(thread_traversal%refinement_elements(1)%cell%geometry, p_dest_element_2%cell%geometry, p_dest_element_3%cell%geometry)
#			endif

			call read_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

#			if !defined(_GT_INPUT_DEST)
				thread_traversal%p_dest_element%coords(:,1) = thread_traversal%p_src_element%nodes(1)%ptr%position
				thread_traversal%p_dest_element%coords(:,2) = 0.5_GRID_SR * (thread_traversal%p_src_element%nodes(1)%ptr%position + thread_traversal%p_src_element%nodes(3)%ptr%position)
				thread_traversal%p_dest_element%coords(:,3) = thread_traversal%p_src_element%nodes(2)%ptr%position
                assert_veq(thread_traversal%p_dest_element%coords, thread_traversal%p_dest_element%coords)
#			endif

			call read_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

			call _GT_REFINE_OP(traversal, dest_section, thread_traversal%p_src_element, thread_traversal%p_dest_element, [1])
			call write_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

#			if !defined(_GT_INPUT_DEST)
                p_dest_element_2%coords(:,1) = thread_traversal%p_dest_element%coords(:,3)
				p_dest_element_2%coords(:,2) = 0.5_GRID_SR * (thread_traversal%p_src_element%nodes(2)%ptr%position + thread_traversal%p_src_element%nodes(3)%ptr%position)
                p_dest_element_2%coords(:,3) = thread_traversal%p_dest_element%coords(:,2)
                assert_veq(p_dest_element_2%coords, p_dest_element_2%coords)
#			endif

			call read_dest_element(traversal, dest_thread, dest_section, p_dest_element_2)

			call _GT_REFINE_OP(traversal, dest_section, thread_traversal%p_src_element, p_dest_element_2, [2, 1])
			call write_dest_element(traversal, dest_thread, dest_section, p_dest_element_2)

#			if !defined(_GT_INPUT_DEST)
                p_dest_element_3%coords(:,1) = p_dest_element_2%coords(:,3)
				p_dest_element_3%coords(:,2) = p_dest_element_2%coords(:,2)
				p_dest_element_3%coords(:,3) = thread_traversal%p_src_element%nodes(3)%ptr%position
                assert_veq(p_dest_element_3%coords, p_dest_element_3%coords)
#			endif

			call read_dest_element(traversal, dest_thread, dest_section, p_dest_element_3)

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

#			if !defined(_GT_INPUT_DEST)
				call create_child_cells(thread_traversal%p_src_element%cell%geometry, thread_traversal%refinement_elements(1)%cell%geometry, thread_traversal%refinement_elements(2)%cell%geometry)
				call create_child_cells(thread_traversal%refinement_elements(1)%cell%geometry, thread_traversal%p_dest_element%cell%geometry, p_dest_element_2%cell%geometry)
				call create_child_cells(thread_traversal%refinement_elements(2)%cell%geometry, p_dest_element_3%cell%geometry, p_dest_element_4%cell%geometry)
#			endif

			call read_src_element(traversal, src_thread, src_section, thread_traversal%p_src_element)

#			if !defined(_GT_INPUT_DEST)
				thread_traversal%p_dest_element%coords(:,1) = thread_traversal%p_src_element%nodes(1)%ptr%position
				thread_traversal%p_dest_element%coords(:,2) = 0.5_GRID_SR * (thread_traversal%p_src_element%nodes(1)%ptr%position + thread_traversal%p_src_element%nodes(2)%ptr%position)
				thread_traversal%p_dest_element%coords(:,3) = 0.5_GRID_SR * (thread_traversal%p_src_element%nodes(1)%ptr%position + thread_traversal%p_src_element%nodes(3)%ptr%position)
                assert_veq(thread_traversal%p_dest_element%coords, thread_traversal%p_dest_element%coords)
#			endif

			call read_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

			call _GT_REFINE_OP(traversal, dest_section, thread_traversal%p_src_element, thread_traversal%p_dest_element, [1, 1])
			call write_dest_element(traversal, dest_thread, dest_section, thread_traversal%p_dest_element)

#			if !defined(_GT_INPUT_DEST)
				p_dest_element_2%coords(:,1) = thread_traversal%p_dest_element%coords(:,3)
				p_dest_element_2%coords(:,2) = thread_traversal%p_dest_element%coords(:,2)
				p_dest_element_2%coords(:,3) = thread_traversal%p_src_element%nodes(2)%ptr%position
                assert_veq(p_dest_element_2%coords, p_dest_element_2%coords)
#			endif

			call read_dest_element(traversal, dest_thread, dest_section, p_dest_element_2)

			call _GT_REFINE_OP(traversal, dest_section, thread_traversal%p_src_element, p_dest_element_2, [1, 2])
			call write_dest_element(traversal, dest_thread, dest_section, p_dest_element_2)

#			if !defined(_GT_INPUT_DEST)
				p_dest_element_3%coords(:,1) = p_dest_element_2%coords(:,3)
				p_dest_element_3%coords(:,2) = 0.5_GRID_SR * (thread_traversal%p_src_element%nodes(2)%ptr%position + thread_traversal%p_src_element%nodes(3)%ptr%position)
				p_dest_element_3%coords(:,3) = p_dest_element_2%coords(:,1)
                assert_veq(p_dest_element_3%coords, p_dest_element_3%coords)
#			endif

			call read_dest_element(traversal, dest_thread, dest_section, p_dest_element_3)

			call _GT_REFINE_OP(traversal, dest_section, thread_traversal%p_src_element, p_dest_element_3, [2, 1])
			call write_dest_element(traversal, dest_thread, dest_section, p_dest_element_3)

#			if !defined(_GT_INPUT_DEST)
				p_dest_element_4%coords(:,1) = p_dest_element_3%coords(:,3)
				p_dest_element_4%coords(:,2) = p_dest_element_3%coords(:,2)
				p_dest_element_4%coords(:,3) = thread_traversal%p_src_element%nodes(3)%ptr%position
                assert_veq(p_dest_element_4%coords, p_dest_element_4%coords)
#			endif

			call read_dest_element(traversal, dest_thread, dest_section, p_dest_element_4)

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
	select case (thread_traversal%p_src_element%cell%geometry%i_entity_types)
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
            select case (thread_traversal%p_src_element%cell%geometry%i_entity_types)
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

	call dest_cell%set_edge_types(int(OLD, 1), iand(1_1, src_cell%get_color_edge_type()), int(NEW, 1))

#	if (_DEBUG_LEVEL > 5)
		_log_write(4, '(A, A)') "  src        (in): ", src_cell%to_string()
		_log_write(4, '(A, A)') "  dest      (out): ", dest_cell%to_string()
#	endif
end subroutine

subroutine create_parent_cell(first_child_cell, second_child_cell, parent_cell)
	type(fine_triangle), intent(in)						:: first_child_cell, second_child_cell
	type(fine_triangle), intent(inout)					:: parent_cell

	integer(kind = BYTE), dimension(-8:8), parameter 	    :: i_plotter_parent_type = [ 3, 2, 1, 8, 7, 6, 5, 4,     0,  -6, -7, -8, -1, -2, -3, -4, -5 ]
	integer(kind = BYTE)			 						:: i_color_edge_type
	integer(kind = BYTE)								 	:: i

	!the parent turtle grammar type can be computed by a simple xor (it's not a symmetrical operator though, so check if the order is correct)
	parent_cell%i_turtle_type = ieor(first_child_cell%i_turtle_type, second_child_cell%i_turtle_type)
	assert_gt(first_child_cell%i_turtle_type, second_child_cell%i_turtle_type)

	!the parent plotter grammar type is always the same as the first grandchild's type
	parent_cell%i_plotter_type = i_plotter_parent_type(first_child_cell%i_plotter_type)
	parent_cell%i_depth = first_child_cell%i_depth - 1

	select case (parent_cell%i_turtle_type)
		case (K)
			parent_cell%i_color_edge_color = second_child_cell%i_color_edge_color
			i_color_edge_type = iand(1_BYTE, second_child_cell%get_color_edge_type())
		case (V, H)
			parent_cell%i_color_edge_color = first_child_cell%i_color_edge_color
			i_color_edge_type = iand(1_BYTE, first_child_cell%get_color_edge_type())
	end select

	call parent_cell%set_edge_types(int(OLD, 1), i_color_edge_type, int(NEW, 1))

#	if (_DEBUG_LEVEL > 5)
		_log_write(5, '(A, A)') "  1st child  (in): ", first_child_cell%to_string()
		_log_write(5, '(A, A)') "  2nd child  (in): ", second_child_cell%to_string()
		_log_write(5, '(A, A)') "  parent    (out): ", parent_cell%to_string()
#	endif
end subroutine

subroutine create_child_cells(parent_cell, first_child_cell, second_child_cell)
	type(fine_triangle), intent(in)						:: parent_cell
	type(fine_triangle), intent(inout)					:: first_child_cell, second_child_cell

	integer(kind = BYTE), dimension(3, 2), parameter		:: i_turtle_child_type = reshape([ H, H, V, V, K, K ], [ 3, 2 ])
	integer(kind = BYTE), dimension(-8:8), parameter 	    :: i_plotter_child_type = [ 3, 2, 1, 8, 7, 6, 5, 4,     0,  -6, -7, -8, -1, -2, -3, -4, -5 ]
	integer(kind = BYTE)			 						:: i_color_edge_type
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

	i_color_edge_type = iand(1, parent_cell%get_color_edge_type())

	select case (parent_cell%i_turtle_type)
		case (K)
			call first_child_cell%set_edge_types(int(OLD, 1), int(NEW, 1), int(NEW, 1))
			call second_child_cell%set_edge_types(int(OLD, 1), i_color_edge_type, int(NEW, 1))
			first_child_cell%i_color_edge_color = RED + GREEN - parent_cell%i_color_edge_color
			second_child_cell%i_color_edge_color = parent_cell%i_color_edge_color
	 	case (V)
			call first_child_cell%set_edge_types(int(OLD, 1), i_color_edge_type, int(NEW, 1))
			call second_child_cell%set_edge_types(int(OLD, 1), i_color_edge_type, int(NEW, 1))
			first_child_cell%i_color_edge_color = parent_cell%i_color_edge_color
			second_child_cell%i_color_edge_color = parent_cell%i_color_edge_color
	 	case (H)
			call first_child_cell%set_edge_types(int(OLD, 1), i_color_edge_type, int(NEW, 1))
			call second_child_cell%set_edge_types(int(OLD, 1), int(OLD, 1), int(NEW, 1))
			first_child_cell%i_color_edge_color = parent_cell%i_color_edge_color
			second_child_cell%i_color_edge_color = RED + GREEN - parent_cell%i_color_edge_color
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

!coordinates are passed through from old to new grid and cannot be disabled
#undef _GT_NO_COORDS
#define _GT_PASS_COORDS

!Skeleton operators are not allowed in adaptive traversal
#if defined(_GT_SKELETON_OP)
#   error A skeleton operator cannot be defined in an adaptive traversal!
#endif

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
        call element%cell%geometry%set_previous_edge_type(int(OLD_BND, 1))
    end if

    select case (element%cell%geometry%get_color_edge_type())
        case (OLD)
            if(thread%indices_stack(element%cell%geometry%i_color_edge_color)%is_empty()) then
                !correct this case to INNER_OLD_BND
                call element%cell%geometry%set_color_edge_type(int(OLD_BND, 1))
            else
                !this is the default case INNER_OLD
                call thread%indices_stack(element%cell%geometry%i_color_edge_color)%pop_data(i_empty)
            end if
        case (NEW)
            call thread%indices_stack(element%cell%geometry%i_color_edge_color)%push_data(element%i_cell)
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

    select case (element%cell%geometry%i_entity_types)
        case (INNER_OLD)
            call write_dest_oon(traversal, thread, section, element)
        case (INNER_NEW)
            call write_dest_onn(traversal, thread, section, element)
        case (INNER_OLD_BND)
            call write_dest_odn(traversal, thread, section, element)
        case (FIRST_NEW, FIRST_OLD_BND)
            element%previous_edge%ptr%remove = .false.
            call write_dest(traversal, thread, section, element)
        case default
            assert(.false.)
    end select
end subroutine

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
#undef _GT_PRE_TRAVERSAL_GRID_OP
#undef _GT_POST_TRAVERSAL_GRID_OP
#undef _GT_PRE_TRAVERSAL_OP
#undef _GT_POST_TRAVERSAL_OP
#undef _GT_PASS_COORDS

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
