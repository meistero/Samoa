! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"


module Stats_list
    use SFC_edge_traversal

    implicit none

#	define _CNT_DATA_TYPE    		type(t_statistics)
#	define _CNT_TYPE_NAME    		t_statistics_list

#	include "Tools_list.f90"
end module

module Conformity
    use SFC_edge_traversal
    use Stats_list

	implicit none

    private
    public t_conformity

    type :: t_thread_data
        type(t_statistics)  :: stats
    end type

    type t_conformity
        type(t_thread_data), pointer    :: threads(:) => null()
        type(t_statistics)              :: stats
        integer                         :: mpi_node_type, mpi_edge_type

        contains

        procedure, pass :: create
        procedure, pass :: check => conformity_check
        procedure, pass :: reduce_stats
        procedure, pass :: clear_stats
        procedure, pass :: destroy
    end type

#   if defined(__GFORTRAN__)
#       define add_edges(x)                         merge(2,1,x)
#   else
        integer(kind = GRID_SI), parameter	        :: add_edges(.true. : .false.) = [2, 1]
#   endif

    contains

    subroutine reduce_stats(conformity, mpi_op, global)
        class(t_conformity)     :: conformity
        integer, intent(in)     :: mpi_op
        logical, intent(in)     :: global

        if (associated(conformity%threads)) then
            call conformity%stats%reduce(conformity%threads(:)%stats, mpi_op)
        end if

        if (global) then
            call conformity%stats%reduce(mpi_op)
        end if
    end subroutine

    subroutine clear_stats(conformity)
        class(t_conformity)     :: conformity
        integer :: i

        if (associated(conformity%threads)) then
            do i = 1, size(conformity%threads)
                call conformity%threads(i)%stats%clear()
            end do
        end if

        call conformity%stats%clear()
    end subroutine

    subroutine create_node_mpi_type(mpi_node_type)
        integer, intent(out)            :: mpi_node_type

#       if defined(_MPI)
            type(t_node_data)                       :: node
            integer                                 :: blocklengths(2), types(2), disps(2), type_size, i_error
            integer (kind = MPI_ADDRESS_KIND)       :: lb, ub

            blocklengths(1) = 1
            blocklengths(2) = 1

            disps(1) = 0
            disps(2) = sizeof(node)

            types(1) = MPI_LB
            types(2) = MPI_UB

            call MPI_Type_struct(2, blocklengths, disps, types, mpi_node_type, i_error); assert_eq(i_error, 0)
            call MPI_Type_commit(mpi_node_type, i_error); assert_eq(i_error, 0)

            call MPI_Type_size(mpi_node_type, type_size, i_error); assert_eq(i_error, 0)
            call MPI_Type_get_extent(mpi_node_type, lb, ub, i_error); assert_eq(i_error, 0)

            assert_eq(0, lb)
            assert_eq(0, type_size)
            assert_eq(sizeof(node), ub)
#       endif
    end subroutine

    subroutine create_edge_mpi_type(mpi_edge_type)
        integer, intent(out)            :: mpi_edge_type

#       if defined(_MPI)
            type(t_edge_data)                       :: edge
            integer                                 :: blocklengths(3), types(3), disps(3), type_size, i_error
            integer (kind = MPI_ADDRESS_KIND)       :: lb, ub

            blocklengths(1) = 1
            blocklengths(2) = 3
            blocklengths(3) = 1

            disps(1) = 0
            disps(2) = loc(edge%refine) - loc(edge)
            disps(3) = sizeof(edge)

            types(1) = MPI_LB
            types(2) = MPI_LOGICAL
            types(3) = MPI_UB

            call MPI_Type_struct(3, blocklengths, disps, types, mpi_edge_type, i_error); assert_eq(i_error, 0)
            call MPI_Type_commit(mpi_edge_type, i_error); assert_eq(i_error, 0)

            call MPI_Type_size(mpi_edge_type, type_size, i_error); assert_eq(i_error, 0)
            call MPI_Type_get_extent(mpi_edge_type, lb, ub, i_error); assert_eq(i_error, 0)

            assert_eq(0, lb)
            assert_eq(3*4, type_size)
            assert_eq(sizeof(edge), ub)
#       endif
    end subroutine

    subroutine create(conformity)
        class(t_conformity) :: conformity

        call create_node_mpi_type(conformity%mpi_node_type)
        call create_edge_mpi_type(conformity%mpi_edge_type)
    end subroutine

    subroutine destroy(conformity)
        class(t_conformity) :: conformity
        integer             :: i_error

#       if defined(_MPI)
            call MPI_Type_free(conformity%mpi_node_type, i_error); assert_eq(i_error, 0)
#       endif

#       if defined(_MPI)
            call MPI_Type_free(conformity%mpi_edge_type, i_error); assert_eq(i_error, 0)
#       endif

        if (associated(conformity%threads)) then
            deallocate(conformity%threads, stat = i_error); assert_eq(i_error, 0)
        end if
    end subroutine

	!*******************
	!conformity check
	!*******************

	!> Corrects refinement array for conformity and returns the number of cells in the destination grid
	subroutine conformity_check(conformity, grid)
		class(t_conformity), intent(inout)   :: conformity
		type(t_grid), intent(inout)         :: grid

		integer (kind = GRID_SI)		    :: i_traversals
        integer                             :: i_error

		_log_write(3, "(2X, A)") "Check conformity..."

        if (.not. associated(conformity%threads)) then
            !$omp barrier

            !$omp single
            allocate(conformity%threads(omp_get_max_threads()), stat = i_error); assert_eq(i_error, 0)
            !$omp end single
        end if

		!mark all edge refinement flags according to the cell refinement flag
		call initial_conformity_traversal(conformity, grid)
        i_traversals = 1

        !update edge refinement flags until the grid is conform
		do
            call update_conformity_traversal(conformity, grid)

            i_traversals = i_traversals + 1

            if (grid%l_conform) then
                exit
            end if
		end do

        !do empty traversals where necessary in order to align the traversal directions
        call empty_traversal(conformity, grid)

 		_log_write(2, "(2X, I0, A)") i_traversals, " conformity traversal(s) performed."
	end subroutine

    !************************
    !grid traversals
    !************************

    subroutine initial_conformity_traversal(conformity, grid)
		class(t_conformity), intent(inout)  :: conformity
        type(t_grid), intent(inout)	        :: grid

        integer (kind = GRID_SI)            :: i_section, i_first_local_section, i_last_local_section
        type(t_statistics)                  :: thread_stats

		_log_write(3, "(3X, A)") "Initial conformity traversal:"

        call grid%get_local_sections(i_first_local_section, i_last_local_section)

        thread_stats%r_traversal_time = -get_wtime()
        thread_stats%r_computation_time = -get_wtime()

        do i_section = i_first_local_section, i_last_local_section
#           if defined(_OPENMP_TASKS)
                !$omp task default(shared) firstprivate(i_section) mergeable
#           endif

            call initial_conformity_traversal_wrapper(grid%threads%elements(i_thread), grid%sections%elements_alloc(i_section))
            call grid%sections%elements_alloc(i_section)%reverse()

#           if defined(_OPENMP_TASKS)
                !$omp end task
#           endif
        end do

#       if defined(_OPENMP_TASKS)
            !$omp taskwait
#       endif

        thread_stats%r_computation_time = thread_stats%r_computation_time + get_wtime()
        thread_stats%r_traversal_time = thread_stats%r_traversal_time + get_wtime()

        conformity%threads(i_thread)%stats = conformity%threads(i_thread)%stats + thread_stats
    end subroutine

    subroutine update_conformity_traversal(conformity, grid)
		class(t_conformity), intent(inout)  :: conformity
        type(t_grid), intent(inout)         :: grid

        integer (kind = GRID_SI)            :: i_section, i_first_local_section, i_last_local_section
        type(t_statistics)                  :: thread_stats

 		_log_write(3, '(3X, A)') "Update conformity traversal:"

        call grid%get_local_sections(i_first_local_section, i_last_local_section)

        thread_stats%r_traversal_time = -get_wtime()
        thread_stats%r_computation_time = -get_wtime()

        do i_section = i_first_local_section, i_last_local_section
            assert_eq(i_section, grid%sections%elements_alloc(i_section)%index)
            call recv_mpi_boundary(grid%sections%elements_alloc(i_section), mpi_node_type_optional=conformity%mpi_node_type)
        end do

        do i_section = i_first_local_section, i_last_local_section
#           if defined(_OPENMP_TASKS)
                !$omp task default(shared) firstprivate(i_section) mergeable
#           endif

            !do a conformity traversal only if it is required
            do while (.not. grid%sections%elements_alloc(i_section)%l_conform)
                call update_conformity_traversal_wrapper(grid%threads%elements(i_thread), grid%sections%elements_alloc(i_section))
                call grid%sections%elements_alloc(i_section)%reverse()
            end do

            call send_mpi_boundary(grid%sections%elements_alloc(i_section), mpi_node_type_optional=conformity%mpi_node_type)

#           if defined(_OPENMP_TASKS)
                !$omp end task
#           endif
        end do

#       if defined(_OPENMP_TASKS)
            !$omp taskwait
#       endif

        thread_stats%r_computation_time = thread_stats%r_computation_time + get_wtime()

        !wait until all computation is done
        !$omp barrier

        thread_stats%r_sync_time = -get_wtime()
        call collect_boundary_data(grid, edge_merge_op_integrity, node_merge_op_integrity, mpi_node_type_optional=conformity%mpi_node_type)

        !wait until all boundary data has been collected
        !$omp barrier

        call duplicate_boundary_data(grid, edge_write_op_integrity, node_write_op_integrity)
        thread_stats%r_sync_time = thread_stats%r_sync_time + get_wtime()

        !wait until all boundary data has been copied
        !$omp barrier

        thread_stats%r_barrier_time = -get_wtime()

        !$omp single
        call reduce(grid%l_conform, grid%sections%elements_alloc%l_conform, MPI_LAND, .true.)
        !$omp end single

        thread_stats%r_barrier_time = thread_stats%r_barrier_time + get_wtime()

        thread_stats%r_traversal_time = thread_stats%r_traversal_time + get_wtime()

        conformity%threads(i_thread)%stats = conformity%threads(i_thread)%stats + thread_stats
    end subroutine

    subroutine empty_traversal(conformity, grid)
		class(t_conformity), intent(inout)  :: conformity
        type(t_grid), intent(inout)         :: grid

        integer (kind = GRID_SI)            :: i_section, i_first_local_section, i_last_local_section, i
        type(t_statistics)                  :: thread_stats

 		_log_write(3, '(3X, A)') "Empty traversal:"

        call grid%get_local_sections(i_first_local_section, i_last_local_section)

        thread_stats%r_traversal_time = -get_wtime()
        thread_stats%r_computation_time = -get_wtime()

        do i_section = i_first_local_section, i_last_local_section
#           if defined(_OPENMP_TASKS)
                !$omp task default(shared) firstprivate(i_section) mergeable
#           endif

            assert_eq(i_section, grid%sections%elements_alloc(i_section)%index)

            !do an empty traversal if the section is backwards
            if (grid%sections%elements_alloc(i_section)%cells%is_forward() .neqv. grid%sections%is_forward()) then
                call empty_traversal_wrapper(grid%threads%elements(i_thread), grid%sections%elements_alloc(i_section))
                call grid%sections%elements_alloc(i_section)%reverse()
            end if

#           if defined(_OPENMP_TASKS)
                !$omp end task
#           endif
        end do

#       if defined(_OPENMP_TASKS)
            !$omp taskwait
#       endif

        thread_stats%r_computation_time = thread_stats%r_computation_time + get_wtime()

        thread_stats%r_barrier_time = -get_wtime()

        !$omp barrier

        !$omp single
        call gather_integrity(grid%t_global_data, grid%sections%elements%t_global_data)
        !$omp end single

        thread_stats%r_barrier_time = thread_stats%r_barrier_time + get_wtime()

        thread_stats%r_traversal_time = thread_stats%r_traversal_time + get_wtime()

        conformity%threads(i_thread)%stats = conformity%threads(i_thread)%stats + thread_stats

        _log_write(3, '(4X, A, I0, A, 2(X, I0), A, 2(X, I0))') "Estimate for #cells: ", grid%dest_cells, ", stack red:", grid%min_dest_stack(RED), grid%max_dest_stack(RED), ", stack green:", grid%min_dest_stack(GREEN), grid%max_dest_stack(GREEN)
    end subroutine

    !************************
    !NUMA wrapper functions
    !************************

    subroutine initial_conformity_traversal_wrapper(thread, section)
		type(t_grid_thread), intent(inout)      :: thread
		type(t_grid_section), intent(inout)     :: section

		type(t_grid_thread)					    :: thread_local
		type(t_grid_section)				    :: section_local

        thread_local = thread
        section_local = section

        call initial_conformity_traversal_section(thread_local, section_local)

        thread = thread_local
        section = section_local
    end subroutine

    subroutine update_conformity_traversal_wrapper(thread, section)
		type(t_grid_thread), intent(inout)      :: thread
		type(t_grid_section), intent(inout)     :: section

		type(t_grid_thread)					    :: thread_local
		type(t_grid_section)				    :: section_local

        thread_local = thread
        section_local = section

        call update_conformity_traversal_section(thread_local, section_local)

        thread = thread_local
        section = section_local
    end subroutine

    subroutine empty_traversal_wrapper(thread, section)
		type(t_grid_thread), intent(inout)      :: thread
		type(t_grid_section), intent(inout)     :: section

		type(t_grid_thread)					    :: thread_local
		type(t_grid_section)				    :: section_local

        thread_local = thread
        section_local = section

        call empty_traversal_section(thread_local, section_local)

        thread = thread_local
        section = section_local
    end subroutine
    !************************
    !initial conformity traversal
    !************************

    subroutine initial_conformity_traversal_section(thread, section)
		type(t_grid_thread), intent(inout)					                :: thread
		type(t_grid_section), intent(inout)					                :: section

		! local variables
		integer (kind = GRID_SI)                                            :: i_dest_stack(RED : GREEN)
		type(t_cell_stream_data), pointer				                    :: p_cell_data

		section%l_conform = .true.
		section%min_dest_stack = 0
		section%max_dest_stack = 0
		section%start_dest_stack = 0
		section%dest_cells = 0

		i_dest_stack = section%start_dest_stack

        select case (section%cells%get_size())
            case (1)
                !process the only element
                p_cell_data => section%cells%next()
                call initial_conformity_traversal_leaf(thread, section, p_cell_data%fine_triangle, i_dest_stack)
            case (2:)
                !process first element
                p_cell_data => section%cells%next()
                call initial_conformity_traversal_leaf(thread, section, p_cell_data%fine_triangle, i_dest_stack)

                do
                    p_cell_data => section%cells%next()

                    select case (p_cell_data%fine_triangle%i_entity_types)
                        case (INNER_OLD)
                            call initial_conformity_traversal_old_leaf(thread, section, p_cell_data%fine_triangle, i_dest_stack)
                        case (INNER_NEW)
                            call initial_conformity_traversal_new_leaf(thread, section, p_cell_data%fine_triangle, i_dest_stack)
                        case (INNER_OLD_BND)
                            call initial_conformity_traversal_old_bnd_leaf(thread, section, p_cell_data%fine_triangle, i_dest_stack)
                        case (INNER_NEW_BND)
                            call initial_conformity_traversal_new_bnd_leaf(thread, section, p_cell_data%fine_triangle, i_dest_stack)
                        case default
                            exit
                    end select
                end do

                !process last element
                call initial_conformity_traversal_leaf(thread, section, p_cell_data%fine_triangle, i_dest_stack)
        end select

		!add +1 to the maximum for any crossed edges that might have been refined
		!we did not take account of those, since they produce at most one more element on the stack

		section%max_dest_stack = section%max_dest_stack + 1
		section%last_dest_cell = section%dest_cells

		!swap start and end dest stack
		section%end_dest_stack = section%start_dest_stack
		section%start_dest_stack = i_dest_stack
    end subroutine

	subroutine initial_conformity_traversal_old_leaf(thread, section, cell, i_dest_stack)
		type(t_grid_thread), intent(inout)			:: thread
		type(t_grid_section), intent(inout)         :: section
		type(fine_triangle), intent(inout)          :: cell
		integer (kind = GRID_SI), intent(inout)     :: i_dest_stack(RED : GREEN)

		!local variables
		type(t_crossed_edge_stream_data), pointer	:: previous_edge, next_edge
		type(t_edge_data), pointer					:: color_edge
		type(t_edge_geometry)	                    :: previous_edge_geometry, color_edge_geometry

		previous_edge => section%crossed_edges_in%current()
		next_edge => section%crossed_edges_in%next()

		color_edge => thread%edges_stack(cell%i_color_edge_color)%pop()

		next_edge%refine = .false.
		next_edge%coarsen = .true.

		previous_edge_geometry = previous_edge%t_edge_geometry
		color_edge_geometry = color_edge%t_edge_geometry

		select case(cell%i_turtle_type)
			case (K)
				call mark_edges(section, cell, previous_edge%t_edge_geometry, next_edge%t_edge_geometry, color_edge%t_edge_geometry)
			case (V)
				call mark_edges(section, cell, previous_edge%t_edge_geometry, color_edge%t_edge_geometry, next_edge%t_edge_geometry)
			case (H)
				call mark_edges(section, cell, color_edge%t_edge_geometry, previous_edge%t_edge_geometry, next_edge%t_edge_geometry)
		end select

		section%l_conform = section%l_conform &
            .and. (previous_edge_geometry%refine .eqv. previous_edge%t_edge_geometry%refine) &
            .and. (previous_edge_geometry%coarsen .eqv. previous_edge%t_edge_geometry%coarsen) &
            .and. (color_edge_geometry%refine .eqv. color_edge%t_edge_geometry%refine) &
            .and. (color_edge_geometry%coarsen .eqv. color_edge%t_edge_geometry%coarsen)

        i_dest_stack(cell%i_color_edge_color) = i_dest_stack(cell%i_color_edge_color) - add_edges(color_edge%refine)
        section%min_dest_stack = min(section%min_dest_stack, i_dest_stack)

		call section%color_edges_out%write(color_edge%t_color_edge_stream_data)

		call cell%reverse_inner()
	end subroutine

	subroutine initial_conformity_traversal_new_leaf(thread, section, cell, i_dest_stack)
		type(t_grid_thread), intent(inout)			:: thread
		type(t_grid_section), intent(inout)			:: section
		type(fine_triangle), intent(inout)			:: cell
		integer (kind = GRID_SI), intent(inout)     :: i_dest_stack(RED : GREEN)

		!local variables
		type(t_crossed_edge_stream_data), pointer	:: previous_edge, next_edge
		type(t_edge_data), pointer					:: color_edge
		type(t_edge_geometry)	                    :: previous_edge_geometry

		previous_edge => section%crossed_edges_in%current()
		next_edge => section%crossed_edges_in%next()

		color_edge => thread%edges_stack(cell%i_color_edge_color)%push()
		call section%color_edges_in%read(color_edge%t_color_edge_stream_data)

		next_edge%refine = .false.
		next_edge%coarsen = .true.

		color_edge%refine = .false.
		color_edge%coarsen = .true.

		previous_edge_geometry = previous_edge%t_edge_geometry

		select case(cell%i_turtle_type)
			case (K)
				call mark_edges(section, cell, previous_edge%t_edge_geometry, next_edge%t_edge_geometry, color_edge%t_edge_geometry)
			case (V)
				call mark_edges(section, cell, previous_edge%t_edge_geometry, color_edge%t_edge_geometry, next_edge%t_edge_geometry)
			case (H)
				call mark_edges(section, cell, color_edge%t_edge_geometry, previous_edge%t_edge_geometry, next_edge%t_edge_geometry)
		end select

		section%l_conform = section%l_conform &
            .and. (previous_edge_geometry%refine .eqv. previous_edge%t_edge_geometry%refine) &
            .and. (previous_edge_geometry%coarsen .eqv. previous_edge%t_edge_geometry%coarsen)

        i_dest_stack(cell%i_color_edge_color) = i_dest_stack(cell%i_color_edge_color) + add_edges(color_edge%refine)
        section%max_dest_stack = max(section%max_dest_stack, i_dest_stack)

		call cell%reverse_inner()
	end subroutine

	subroutine initial_conformity_traversal_old_bnd_leaf(thread, section, cell, i_dest_stack)
		type(t_grid_thread), intent(inout)			:: thread
		type(t_grid_section), intent(inout)			:: section
		type(fine_triangle), intent(inout)			:: cell
		integer (kind = GRID_SI), intent(inout)     :: i_dest_stack(RED : GREEN)

		!local variables
		type(t_edge_data), pointer					:: color_edge
		type(t_crossed_edge_stream_data), pointer	:: previous_edge, next_edge
		type(t_edge_geometry)	                    :: previous_edge_geometry

		previous_edge => section%crossed_edges_in%current()
		next_edge => section%crossed_edges_in%next()

        color_edge => section%boundary_edges(cell%i_color_edge_color)%next()

		next_edge%refine = .false.
		next_edge%coarsen = .true.

		color_edge%refine = .false.
		color_edge%coarsen = .true.

		previous_edge_geometry = previous_edge%t_edge_geometry

		select case(cell%i_turtle_type)
			case (K)
				call mark_edges(section, cell, previous_edge%t_edge_geometry, next_edge%t_edge_geometry, color_edge%t_edge_geometry)
			case (V)
				call mark_edges(section, cell, previous_edge%t_edge_geometry, color_edge%t_edge_geometry, next_edge%t_edge_geometry)
			case (H)
				call mark_edges(section, cell, color_edge%t_edge_geometry, previous_edge%t_edge_geometry, next_edge%t_edge_geometry)
		end select

		section%l_conform = section%l_conform &
            .and. (previous_edge_geometry%refine .eqv. previous_edge%t_edge_geometry%refine) &
            .and. (previous_edge_geometry%coarsen .eqv. previous_edge%t_edge_geometry%coarsen)

        i_dest_stack(cell%i_color_edge_color) = i_dest_stack(cell%i_color_edge_color) - add_edges(color_edge%refine)
        section%min_dest_stack = min(section%min_dest_stack, i_dest_stack)

		call cell%reverse_inner()
	end subroutine

	subroutine initial_conformity_traversal_new_bnd_leaf(thread, section, cell, i_dest_stack)
		type(t_grid_thread), intent(inout)			:: thread
		type(t_grid_section), intent(inout)			:: section
		type(fine_triangle), intent(inout)			:: cell
		integer (kind = GRID_SI), intent(inout)     :: i_dest_stack(RED : GREEN)

		!local variables
		type(t_edge_data), pointer					:: color_edge
		type(t_crossed_edge_stream_data), pointer	:: previous_edge, next_edge
		type(t_edge_geometry)	                    :: previous_edge_geometry

		previous_edge => section%crossed_edges_in%current()
		next_edge => section%crossed_edges_in%next()

		color_edge => section%boundary_edges(cell%i_color_edge_color)%next()

		next_edge%refine = .false.
		next_edge%coarsen = .true.

		color_edge%refine = .false.
		color_edge%coarsen = .true.

		previous_edge_geometry = previous_edge%t_edge_geometry

		select case(cell%i_turtle_type)
			case (K)
				call mark_edges(section, cell, previous_edge%t_edge_geometry, next_edge%t_edge_geometry, color_edge%t_edge_geometry)
			case (V)
				call mark_edges(section, cell, previous_edge%t_edge_geometry, color_edge%t_edge_geometry, next_edge%t_edge_geometry)
			case (H)
				call mark_edges(section, cell, color_edge%t_edge_geometry, previous_edge%t_edge_geometry, next_edge%t_edge_geometry)
		end select

		section%l_conform = section%l_conform &
            .and. (previous_edge_geometry%refine .eqv. previous_edge%t_edge_geometry%refine) &
            .and. (previous_edge_geometry%coarsen .eqv. previous_edge%t_edge_geometry%coarsen)

		i_dest_stack(cell%i_color_edge_color) = i_dest_stack(cell%i_color_edge_color) + add_edges(color_edge%refine)
        section%max_dest_stack = max(section%max_dest_stack, i_dest_stack)

		call cell%reverse_inner()
	end subroutine

	subroutine initial_conformity_traversal_leaf(thread, section, cell, i_dest_stack)
		type(t_grid_thread), intent(inout)			:: thread
		type(t_grid_section), intent(inout)			:: section
		type(fine_triangle), intent(inout)			:: cell
		integer (kind = GRID_SI), intent(inout)     :: i_dest_stack(RED : GREEN)

		!local variables
		integer(BYTE)							:: i_previous_edge_type, i_color_edge_type, i_next_edge_type
		type(t_edge_data), pointer					:: color_edge, p_boundary_edge
		type(t_crossed_edge_stream_data), pointer	:: previous_edge, next_edge
		type(t_edge_geometry)	                    :: previous_edge_geometry, color_edge_geometry

		call cell%get_edge_types(i_previous_edge_type, i_color_edge_type, i_next_edge_type)

		select case(i_previous_edge_type)
			case (OLD)
				previous_edge => section%crossed_edges_in%current()

                previous_edge_geometry = previous_edge%t_edge_geometry
			case (OLD_BND)
				p_boundary_edge => section%boundary_edges(RED)%next()
				previous_edge => p_boundary_edge%t_crossed_edge_stream_data

				previous_edge%refine = .false.
				previous_edge%coarsen = .true.
		end select

		select case(i_color_edge_type)
			case (OLD)
				color_edge => thread%edges_stack(cell%i_color_edge_color)%pop()

				color_edge_geometry = color_edge%t_edge_geometry
			case (NEW)
				color_edge => thread%edges_stack(cell%i_color_edge_color)%push()
				call section%color_edges_in%read(color_edge%t_color_edge_stream_data)

				color_edge%refine = .false.
				color_edge%coarsen = .true.
			case (OLD_BND)
				color_edge => section%boundary_edges(cell%i_color_edge_color)%next()

				color_edge%refine = .false.
				color_edge%coarsen = .true.
			case (NEW_BND)
				color_edge => section%boundary_edges(cell%i_color_edge_color)%next()

				color_edge%refine = .false.
				color_edge%coarsen = .true.
		end select

		select case(i_next_edge_type)
			case (NEW)
				next_edge => section%crossed_edges_in%next()
			case (NEW_BND)
				p_boundary_edge => section%boundary_edges(RED)%next()
				next_edge => p_boundary_edge%t_crossed_edge_stream_data
		end select

		next_edge%refine = .false.
		next_edge%coarsen = .true.

		select case(cell%i_turtle_type)
			case (K)
				call mark_edges(section, cell, previous_edge%t_edge_geometry, next_edge%t_edge_geometry, color_edge%t_edge_geometry)
			case (V)
				call mark_edges(section, cell, previous_edge%t_edge_geometry, color_edge%t_edge_geometry, next_edge%t_edge_geometry)
			case (H)
				call mark_edges(section, cell, color_edge%t_edge_geometry, previous_edge%t_edge_geometry, next_edge%t_edge_geometry)
		end select

		select case(i_previous_edge_type)
            case (OLD_BND)
                i_dest_stack = i_dest_stack - [1, 2]
                section%min_dest_stack = min(section%min_dest_stack, i_dest_stack)
			case (OLD)
                section%l_conform = section%l_conform &
                    .and. (previous_edge_geometry%refine .eqv. previous_edge%refine) &
                    .and. (previous_edge_geometry%coarsen .eqv. previous_edge%coarsen)
		end select

		select case(i_color_edge_type)
			case (OLD)
                section%l_conform = section%l_conform &
                    .and. (color_edge_geometry%refine .eqv. color_edge%refine) &
                    .and. (color_edge_geometry%coarsen .eqv. color_edge%coarsen)

				i_dest_stack(cell%i_color_edge_color) = i_dest_stack(cell%i_color_edge_color) - add_edges(color_edge%refine)
                section%min_dest_stack = min(section%min_dest_stack, i_dest_stack)

				call section%color_edges_out%write(color_edge%t_color_edge_stream_data)
            case (OLD_BND)
				i_dest_stack(cell%i_color_edge_color) = i_dest_stack(cell%i_color_edge_color) - add_edges(color_edge%refine)
                section%min_dest_stack = min(section%min_dest_stack, i_dest_stack)
			case (NEW, NEW_BND)
                i_dest_stack(cell%i_color_edge_color) = i_dest_stack(cell%i_color_edge_color) + add_edges(color_edge%refine)
                section%max_dest_stack = max(section%max_dest_stack, i_dest_stack)
		end select

		select case(i_next_edge_type)
			case (NEW_BND)
                i_dest_stack = i_dest_stack + [1, 2]
                section%max_dest_stack = max(i_dest_stack, section%max_dest_stack)
		end select

		call cell%reverse()
	end subroutine

    !************************
    !update conformity traversal
    !************************

	subroutine update_conformity_traversal_section(thread, section)
		type(t_grid_thread), intent(inout)			:: thread
		type(t_grid_section), intent(inout)			:: section

		! local variables
		integer (kind = GRID_SI)                    :: i_dest_stack(RED : GREEN)
		type(t_cell_stream_data), pointer			:: p_cell_data

		section%l_conform = .true.
		section%min_dest_stack = 0
		section%max_dest_stack = 0
		section%start_dest_stack = 0
		section%dest_cells = 0

		i_dest_stack = section%start_dest_stack

        select case (section%cells%get_size())
            case (1)
                !process the only element
                p_cell_data => section%cells%next()
                call update_conformity_traversal_leaf(thread, section, p_cell_data%fine_triangle, i_dest_stack)
            case (2:)
                !process first element
                p_cell_data => section%cells%next()
                call update_conformity_traversal_leaf(thread, section, p_cell_data%fine_triangle, i_dest_stack)

                do
                    p_cell_data => section%cells%next()

                    select case (p_cell_data%fine_triangle%i_entity_types)
                        case (INNER_OLD)
                            call update_conformity_traversal_old_leaf(thread, section, p_cell_data%fine_triangle, i_dest_stack)
                        case (INNER_NEW)
                            call update_conformity_traversal_new_leaf(thread, section, p_cell_data%fine_triangle, i_dest_stack)
                        case (INNER_OLD_BND)
                            call update_conformity_traversal_old_bnd_leaf(thread, section, p_cell_data%fine_triangle, i_dest_stack)
                        case (INNER_NEW_BND)
                            call update_conformity_traversal_new_bnd_leaf(thread, section, p_cell_data%fine_triangle, i_dest_stack)
                        case default
                            exit
                    end select
                end do

                !process last element
                call update_conformity_traversal_leaf(thread, section, p_cell_data%fine_triangle, i_dest_stack)
        end select

		!add +1 to the maximum for any crossed edges that might have been refined
		!we did not take account of those, since they produce at most one more element on the stack

		section%max_dest_stack = section%max_dest_stack + 1
		section%last_dest_cell = section%dest_cells

		!swap start and end dest stack
		section%end_dest_stack = section%start_dest_stack
		section%start_dest_stack = i_dest_stack
	end subroutine

	subroutine update_conformity_traversal_old_leaf(thread, section, cell, i_dest_stack)
		type(t_grid_thread), intent(inout)			:: thread
		type(t_grid_section), intent(inout)			:: section
		type(fine_triangle), intent(inout)			:: cell
		integer (kind = GRID_SI), intent(inout)     :: i_dest_stack(RED : GREEN)

		!local variables
		type(t_crossed_edge_stream_data), pointer	:: previous_edge, next_edge
		type(t_edge_data), pointer					:: color_edge
		type(t_edge_geometry)	                    :: previous_edge_geometry, color_edge_geometry

		previous_edge => section%crossed_edges_in%current()
		next_edge => section%crossed_edges_in%next()

		color_edge => thread%edges_stack(cell%i_color_edge_color)%pop()

		previous_edge_geometry = previous_edge%t_edge_geometry
		color_edge_geometry = color_edge%t_edge_geometry

		! check marked legs and mark hypotenuse if needed
		select case (cell%i_turtle_type)
			case (K)
				call mark_edges_integrity(section, cell, previous_edge%t_edge_geometry, next_edge%t_edge_geometry, color_edge%t_edge_geometry)
			case (V)
				call mark_edges_integrity(section, cell, previous_edge%t_edge_geometry, color_edge%t_edge_geometry, next_edge%t_edge_geometry)
			case (H)
				call mark_edges_integrity(section, cell, color_edge%t_edge_geometry, previous_edge%t_edge_geometry, next_edge%t_edge_geometry)
		end select

		section%l_conform = section%l_conform &
            .and. (previous_edge_geometry%refine .eqv. previous_edge%t_edge_geometry%refine) &
            .and. (previous_edge_geometry%coarsen .eqv. previous_edge%t_edge_geometry%coarsen) &
            .and. (color_edge_geometry%refine .eqv. color_edge%t_edge_geometry%refine) &
            .and. (color_edge_geometry%coarsen .eqv. color_edge%t_edge_geometry%coarsen)

        i_dest_stack(cell%i_color_edge_color) = i_dest_stack(cell%i_color_edge_color) - add_edges(color_edge%refine)
        section%min_dest_stack = min(i_dest_stack, section%min_dest_stack)

		call section%color_edges_out%write(color_edge%t_color_edge_stream_data)

		call cell%reverse_inner()
	end subroutine

	subroutine update_conformity_traversal_new_leaf(thread, section, cell, i_dest_stack)
		type(t_grid_thread), intent(inout)			:: thread
		type(t_grid_section), intent(inout)			:: section
		type(fine_triangle), intent(inout)			:: cell
		integer (kind = GRID_SI), intent(inout)     :: i_dest_stack(RED : GREEN)

		!local variables
		type(t_crossed_edge_stream_data), pointer	        :: previous_edge, next_edge
		type(t_edge_data), pointer					        :: color_edge
		type(t_edge_geometry)	                            :: previous_edge_geometry

		previous_edge => section%crossed_edges_in%current()
		next_edge => section%crossed_edges_in%next()

		color_edge => thread%edges_stack(cell%i_color_edge_color)%push()
		call section%color_edges_in%read(color_edge%t_color_edge_stream_data)

		previous_edge_geometry = previous_edge%t_edge_geometry

		! check marked legs and mark hypotenuse if needed
		select case (cell%i_turtle_type)
			case (K)
				call mark_edges_integrity(section, cell, previous_edge%t_edge_geometry, next_edge%t_edge_geometry, color_edge%t_edge_geometry)
			case (V)
				call mark_edges_integrity(section, cell, previous_edge%t_edge_geometry, color_edge%t_edge_geometry, next_edge%t_edge_geometry)
			case (H)
				call mark_edges_integrity(section, cell, color_edge%t_edge_geometry, previous_edge%t_edge_geometry, next_edge%t_edge_geometry)
		end select

		section%l_conform = section%l_conform &
            .and. (previous_edge_geometry%refine .eqv. previous_edge%t_edge_geometry%refine) &
            .and. (previous_edge_geometry%coarsen .eqv. previous_edge%t_edge_geometry%coarsen)

        i_dest_stack(cell%i_color_edge_color) = i_dest_stack(cell%i_color_edge_color) + add_edges(color_edge%refine)
        section%max_dest_stack = max(i_dest_stack, section%max_dest_stack)

		call cell%reverse_inner()
	end subroutine

	subroutine update_conformity_traversal_old_bnd_leaf(thread, section, cell, i_dest_stack)
		type(t_grid_thread), intent(inout)			:: thread
		type(t_grid_section), intent(inout)			:: section
		type(fine_triangle), intent(inout)			:: cell
		integer (kind = GRID_SI), intent(inout)     :: i_dest_stack(RED : GREEN)

		!local variables
		type(t_crossed_edge_stream_data), pointer	:: previous_edge, next_edge
		type(t_edge_data), pointer					:: color_edge
		type(t_edge_geometry)	                    :: previous_edge_geometry

		previous_edge => section%crossed_edges_in%current()
		next_edge => section%crossed_edges_in%next()

        color_edge => section%boundary_edges(cell%i_color_edge_color)%next()

		previous_edge_geometry = previous_edge%t_edge_geometry

		! check marked legs and mark hypotenuse if needed
		select case (cell%i_turtle_type)
			case (K)
				call mark_edges_integrity(section, cell, previous_edge%t_edge_geometry, next_edge%t_edge_geometry, color_edge%t_edge_geometry)
			case (V)
				call mark_edges_integrity(section, cell, previous_edge%t_edge_geometry, color_edge%t_edge_geometry, next_edge%t_edge_geometry)
			case (H)
				call mark_edges_integrity(section, cell, color_edge%t_edge_geometry, previous_edge%t_edge_geometry, next_edge%t_edge_geometry)
		end select

		section%l_conform = section%l_conform &
            .and. (previous_edge_geometry%refine .eqv. previous_edge%t_edge_geometry%refine) &
            .and. (previous_edge_geometry%coarsen .eqv. previous_edge%t_edge_geometry%coarsen)

        i_dest_stack(cell%i_color_edge_color) = i_dest_stack(cell%i_color_edge_color) - add_edges(color_edge%refine)
        section%min_dest_stack = min(i_dest_stack, section%min_dest_stack)

		call cell%reverse_inner()
	end subroutine

	subroutine update_conformity_traversal_new_bnd_leaf(thread, section, cell, i_dest_stack)
		type(t_grid_thread), intent(inout)			:: thread
		type(t_grid_section), intent(inout)			:: section
		type(fine_triangle), intent(inout)			:: cell
		integer (kind = GRID_SI), intent(inout)     :: i_dest_stack(RED : GREEN)

		!local variables
		type(t_crossed_edge_stream_data), pointer	    :: previous_edge, next_edge
		type(t_edge_data), pointer					    :: color_edge
		type(t_edge_geometry)	                        :: previous_edge_geometry

		previous_edge => section%crossed_edges_in%current()
		next_edge => section%crossed_edges_in%next()

		color_edge => section%boundary_edges(cell%i_color_edge_color)%next()

		previous_edge_geometry = previous_edge%t_edge_geometry

		! check marked legs and mark hypotenuse if needed
		select case (cell%i_turtle_type)
			case (K)
				call mark_edges_integrity(section, cell, previous_edge%t_edge_geometry, next_edge%t_edge_geometry, color_edge%t_edge_geometry)
			case (V)
				call mark_edges_integrity(section, cell, previous_edge%t_edge_geometry, color_edge%t_edge_geometry, next_edge%t_edge_geometry)
			case (H)
				call mark_edges_integrity(section, cell, color_edge%t_edge_geometry, previous_edge%t_edge_geometry, next_edge%t_edge_geometry)
		end select

		section%l_conform = section%l_conform &
            .and. (previous_edge_geometry%refine .eqv. previous_edge%t_edge_geometry%refine) &
            .and. (previous_edge_geometry%coarsen .eqv. previous_edge%t_edge_geometry%coarsen)

        i_dest_stack(cell%i_color_edge_color) = i_dest_stack(cell%i_color_edge_color) + add_edges(color_edge%refine)
        section%max_dest_stack = max(i_dest_stack, section%max_dest_stack)

		call cell%reverse_inner()
	end subroutine

	subroutine update_conformity_traversal_leaf(thread, section, cell, i_dest_stack)
		type(t_grid_thread), intent(inout)			:: thread
		type(t_grid_section), intent(inout)			:: section
		type(fine_triangle), intent(inout)			:: cell
		integer (kind = GRID_SI), intent(inout)     :: i_dest_stack(RED : GREEN)

		!local variables
		integer(BYTE)							:: i_previous_edge_type, i_color_edge_type, i_next_edge_type
		type(t_edge_data), pointer				    :: color_edge, p_boundary_edge
		type(t_crossed_edge_stream_data), pointer   :: previous_edge, next_edge
		type(t_edge_geometry)	                    :: previous_edge_geometry, color_edge_geometry

		call cell%get_edge_types(i_previous_edge_type, i_color_edge_type, i_next_edge_type)

		select case(i_previous_edge_type)
			case (OLD)
				previous_edge => section%crossed_edges_in%current()

                previous_edge_geometry = previous_edge%t_edge_geometry
			case (OLD_BND)
				p_boundary_edge => section%boundary_edges(RED)%next()
				previous_edge => p_boundary_edge%t_crossed_edge_stream_data
		end select

		select case(i_color_edge_type)
			case (OLD)
				color_edge => thread%edges_stack(cell%i_color_edge_color)%pop()

                color_edge_geometry = color_edge%t_edge_geometry
 			case (NEW)
				color_edge => thread%edges_stack(cell%i_color_edge_color)%push()
				call section%color_edges_in%read(color_edge%t_color_edge_stream_data)
			case (OLD_BND, NEW_BND)
				color_edge => section%boundary_edges(cell%i_color_edge_color)%next()
		end select

		select case(i_next_edge_type)
			case (NEW)
				next_edge => section%crossed_edges_in%next()
			case (NEW_BND)
				p_boundary_edge => section%boundary_edges(RED)%next()
				next_edge => p_boundary_edge%t_crossed_edge_stream_data
		end select

		! check marked legs and mark hypotenuse if needed
		select case (cell%i_turtle_type)
			case (K)
				call mark_edges_integrity(section, cell, previous_edge%t_edge_geometry, next_edge%t_edge_geometry, color_edge%t_edge_geometry)
			case (V)
				call mark_edges_integrity(section, cell, previous_edge%t_edge_geometry, color_edge%t_edge_geometry, next_edge%t_edge_geometry)
			case (H)
				call mark_edges_integrity(section, cell, color_edge%t_edge_geometry, previous_edge%t_edge_geometry, next_edge%t_edge_geometry)
		end select

		select case(i_previous_edge_type)
			case (OLD)
                section%l_conform = section%l_conform &
                    .and. (previous_edge_geometry%refine .eqv. previous_edge%t_edge_geometry%refine) &
                    .and. (previous_edge_geometry%coarsen .eqv. previous_edge%t_edge_geometry%coarsen)
            case (OLD_BND)
                i_dest_stack = i_dest_stack - [1, 2]
                section%min_dest_stack = min(i_dest_stack, section%min_dest_stack)
		end select

		select case(i_color_edge_type)
			case (OLD)
                section%l_conform = section%l_conform &
                    .and. (color_edge_geometry%refine .eqv. color_edge%t_edge_geometry%refine) &
                    .and. (color_edge_geometry%coarsen .eqv. color_edge%t_edge_geometry%coarsen)

                i_dest_stack(cell%i_color_edge_color) = i_dest_stack(cell%i_color_edge_color) - add_edges(color_edge%refine)
                section%min_dest_stack = min(i_dest_stack, section%min_dest_stack)

				call section%color_edges_out%write(color_edge%t_color_edge_stream_data)
            case (OLD_BND)
                i_dest_stack(cell%i_color_edge_color) = i_dest_stack(cell%i_color_edge_color) - add_edges(color_edge%refine)
                section%min_dest_stack = min(i_dest_stack, section%min_dest_stack)
			case (NEW, NEW_BND)
				i_dest_stack(cell%i_color_edge_color) = i_dest_stack(cell%i_color_edge_color) + add_edges(color_edge%refine)
                section%max_dest_stack = max(i_dest_stack, section%max_dest_stack)
		end select

		select case(i_next_edge_type)
			case (NEW_BND)
                i_dest_stack = i_dest_stack + [1, 2]
                section%max_dest_stack = max(i_dest_stack, section%max_dest_stack)
		end select

		call cell%reverse()
	end subroutine

    !************************
    !empty traversal
    !************************

	subroutine empty_traversal_section(thread, section)
		type(t_grid_thread), intent(inout)			:: thread
		type(t_grid_section), intent(inout)			:: section
		integer (kind = GRID_SI)                    :: i_dest_stack(RED : GREEN)

		! local variables
		type(t_cell_stream_data), pointer			:: p_cell_data

        select case (section%cells%get_size())
            case (1)
                !process the only element
                p_cell_data => section%cells%next()
                call empty_traversal_leaf(thread, section, p_cell_data%fine_triangle)
            case (2:)
                !process first element
                p_cell_data => section%cells%next()
                call empty_traversal_leaf(thread, section, p_cell_data%fine_triangle)

                do
                    p_cell_data => section%cells%next()

                    select case (p_cell_data%fine_triangle%i_entity_types)
                        case (INNER_OLD)
                            call empty_traversal_old_leaf(thread, section, p_cell_data%fine_triangle)
                        case (INNER_NEW)
                            call empty_traversal_new_leaf(thread, section, p_cell_data%fine_triangle)
                        case (INNER_OLD_BND)
                            call empty_traversal_old_bnd_leaf(thread, section, p_cell_data%fine_triangle)
                        case (INNER_NEW_BND)
                            call empty_traversal_new_bnd_leaf(thread, section, p_cell_data%fine_triangle)
                        case default
                            exit
                    end select
                end do

                !process last element
                call empty_traversal_leaf(thread, section, p_cell_data%fine_triangle)
        end select

		!swap start and end dest stack
		i_dest_stack = section%end_dest_stack
		section%end_dest_stack = section%start_dest_stack
		section%start_dest_stack = i_dest_stack
	end subroutine

	subroutine empty_traversal_old_leaf(thread, section, cell)
		type(t_grid_thread), intent(inout)			:: thread
		type(t_grid_section), intent(inout)			:: section
		type(fine_triangle), intent(inout)			:: cell

		!local variables
		type(t_crossed_edge_stream_data), pointer	:: next_edge
		type(t_edge_data), pointer					:: color_edge

		next_edge => section%crossed_edges_in%next()
		color_edge => thread%edges_stack(cell%i_color_edge_color)%pop()
		call section%color_edges_out%write(color_edge%t_color_edge_stream_data)

        call cell%reverse_refinement()
		call cell%reverse_inner()
	end subroutine

	subroutine empty_traversal_new_leaf(thread, section, cell)
		type(t_grid_thread), intent(inout)			:: thread
		type(t_grid_section), intent(inout)			:: section
		type(fine_triangle), intent(inout)			:: cell

		!local variables
		type(t_crossed_edge_stream_data), pointer	                        :: next_edge
		type(t_edge_data), pointer					                        :: color_edge

		next_edge => section%crossed_edges_in%next()
		color_edge => thread%edges_stack(cell%i_color_edge_color)%push()
		call section%color_edges_in%read(color_edge%t_color_edge_stream_data)

        call cell%reverse_refinement()
		call cell%reverse_inner()
	end subroutine

	subroutine empty_traversal_old_bnd_leaf(thread, section, cell)
		type(t_grid_thread), intent(inout)			:: thread
		type(t_grid_section), intent(inout)			:: section
		type(fine_triangle), intent(inout)			:: cell

		!local variables
		type(t_crossed_edge_stream_data), pointer	:: next_edge
		type(t_edge_data), pointer					:: color_edge

		next_edge => section%crossed_edges_in%next()
        color_edge => section%boundary_edges(cell%i_color_edge_color)%next()

        call cell%reverse_refinement()
		call cell%reverse_inner()
	end subroutine

	subroutine empty_traversal_new_bnd_leaf(thread, section, cell)
		type(t_grid_thread), intent(inout)			:: thread
		type(t_grid_section), intent(inout)			:: section
		type(fine_triangle), intent(inout)			:: cell

		!local variables
		type(t_crossed_edge_stream_data), pointer	                        :: next_edge
		type(t_edge_data), pointer					                        :: color_edge

		next_edge => section%crossed_edges_in%next()
		color_edge => section%boundary_edges(cell%i_color_edge_color)%next()

        call cell%reverse_refinement()
		call cell%reverse_inner()
	end subroutine

	subroutine empty_traversal_leaf(thread, section, cell)
		type(t_grid_thread), intent(inout)			:: thread
		type(t_grid_section), intent(inout)			:: section
		type(fine_triangle), intent(inout)			:: cell

		!local variables
		integer(BYTE)							:: i_previous_edge_type, i_color_edge_type, i_next_edge_type
		type(t_edge_data), pointer				    :: color_edge, p_boundary_edge
		type(t_crossed_edge_stream_data), pointer   :: next_edge

		call cell%get_edge_types(i_previous_edge_type, i_color_edge_type, i_next_edge_type)

		select case(i_previous_edge_type)
			case (OLD_BND)
				p_boundary_edge => section%boundary_edges(RED)%next()
		end select

		select case(i_color_edge_type)
			case (OLD)
				color_edge => thread%edges_stack(cell%i_color_edge_color)%pop()
				call section%color_edges_out%write(color_edge%t_color_edge_stream_data)
 			case (NEW)
				color_edge => thread%edges_stack(cell%i_color_edge_color)%push()
				call section%color_edges_in%read(color_edge%t_color_edge_stream_data)
			case (OLD_BND, NEW_BND)
				p_boundary_edge => section%boundary_edges(cell%i_color_edge_color)%next()
		end select

		select case(i_next_edge_type)
			case (NEW)
				next_edge => section%crossed_edges_in%next()
			case (NEW_BND)
				p_boundary_edge => section%boundary_edges(RED)%next()
		end select

        call cell%reverse_refinement()
		call cell%reverse()
	end subroutine

    !************************
    !integrity operators
    !************************

	subroutine mark_edges(section, cell, left_leg, hypotenuse, right_leg)
		type(t_grid_section), intent(inout)				:: section
		type(fine_triangle), intent(inout)				:: cell
		type(t_edge_geometry), intent(inout)			:: hypotenuse, left_leg, right_leg

		!the hypotenuse cannot be coarsened
		hypotenuse%coarsen = .false.

        !set edge flags
		select case (cell%refinement)
			case (0)
				left_leg%coarsen = .false.
			case (1)
				hypotenuse%refine = .true.
			case (2)
				left_leg%refine = .true.
			case (3)
				right_leg%refine = .true.
			case (4)
				left_leg%refine = .true.
				right_leg%refine = .true.
		end select

		call mark_edges_integrity(section, cell, left_leg, hypotenuse, right_leg)
	end subroutine

	! integrity check for refinement
	subroutine mark_edges_integrity(section, cell, left_leg, hypotenuse, right_leg)
		type(t_grid_section), intent(inout)				:: section
		type(fine_triangle), intent(inout)				:: cell
		type(t_edge_geometry), intent(inout)			:: hypotenuse, left_leg, right_leg

		hypotenuse%refine = hypotenuse%refine .or. left_leg%refine .or. right_leg%refine
		left_leg%coarsen = .not. hypotenuse%refine .and. left_leg%coarsen .and. right_leg%coarsen .and. (left_leg%remove .or. right_leg%remove)
		right_leg%coarsen = left_leg%coarsen

		if (hypotenuse%refine) then
			if (left_leg%refine) then
				if (right_leg%refine) then
					cell%refinement = 4

					!increase the number of cells
					section%dest_cells = section%dest_cells + 4
				else
					cell%refinement = 3

					!increase the number of cells
					section%dest_cells = section%dest_cells + 3
				endif
			else if (right_leg%refine) then
                cell%refinement = 2

                !increase the number of cells
                section%dest_cells = section%dest_cells + 3
			else
                cell%refinement = 1

                !increase the number of cells
                section%dest_cells = section%dest_cells + 2
			endif
		else if (left_leg%coarsen) then
            cell%refinement = -1

            if (left_leg%remove) then
                !increase the number of cells
                section%dest_cells = section%dest_cells + 1
            end if
		else
            cell%refinement = 0

            !increase the number of cells
            section%dest_cells = section%dest_cells + 1
		endif
	end subroutine

    function edge_merge_op_integrity(local_edges, neighbor_edges) result(l_conform)
        type(t_edge_data), intent(inout)    :: local_edges
        type(t_edge_data), intent(in)       :: neighbor_edges
        logical                             :: l_conform

        assert(.not. local_edges%remove)
        assert(.not. neighbor_edges%remove)

        l_conform = (local_edges%refine .or. .not. neighbor_edges%refine) .and. (.not. local_edges%coarsen .or. neighbor_edges%coarsen)

        local_edges%refine = local_edges%refine .or. neighbor_edges%refine
        local_edges%coarsen = local_edges%coarsen .and. neighbor_edges%coarsen
    end function

    function edge_write_op_integrity(local_edges, neighbor_edges) result(l_conform)
        type(t_edge_data), intent(inout)    :: local_edges
        type(t_edge_data), intent(in)       :: neighbor_edges
        logical                             :: l_conform

        assert(.not. local_edges%remove)
        assert(.not. neighbor_edges%remove)
        assert(.not. local_edges%refine .or. neighbor_edges%refine)
        assert(local_edges%coarsen .or. .not. neighbor_edges%coarsen)

        l_conform = (local_edges%refine .eqv. neighbor_edges%refine) .and. (local_edges%coarsen .eqv. neighbor_edges%coarsen)

        local_edges%refine = neighbor_edges%refine
        local_edges%coarsen = neighbor_edges%coarsen
    end function

    function node_merge_op_integrity(local_node, neighbor_node) result(l_conform)
        type(t_node_data), intent(inout)    :: local_node
        type(t_node_data), intent(in)       :: neighbor_node
        logical                             :: l_conform

        l_conform = .true.
    end function

    function node_write_op_integrity(local_node, neighbor_node) result(l_conform)
        type(t_node_data), intent(inout)    :: local_node
        type(t_node_data), intent(in)       :: neighbor_node
        logical                             :: l_conform

        l_conform = .true.
    end function

    subroutine gather_integrity(grid_data, section_data)
        type(t_global_data), intent(inout)	    :: grid_data
        type(t_global_data), intent(inout)		:: section_data(:)

		integer (BYTE)					        :: i_color

        call reduce(grid_data%dest_cells, section_data%dest_cells, MPI_SUM, .false.)

        do i_color = RED, GREEN
            section_data%start_dest_stack(i_color) = section_data%start_dest_stack(i_color) - section_data%end_dest_stack(i_color)
            section_data%min_dest_stack(i_color) = section_data%min_dest_stack(i_color) - section_data%end_dest_stack(i_color)
            section_data%max_dest_stack(i_color) = section_data%max_dest_stack(i_color) - section_data%end_dest_stack(i_color)
            section_data%end_dest_stack(i_color) = -section_data%start_dest_stack(i_color)

            call prefix_sum(section_data%end_dest_stack(i_color), section_data%end_dest_stack(i_color))

            section_data%start_dest_stack(i_color) = section_data%start_dest_stack(i_color) + section_data%end_dest_stack(i_color)
            section_data%min_dest_stack(i_color) = section_data%min_dest_stack(i_color) + section_data%end_dest_stack(i_color)
            section_data%max_dest_stack(i_color) = section_data%max_dest_stack(i_color) + section_data%end_dest_stack(i_color)

            grid_data%min_dest_stack(i_color) = minval(section_data%min_dest_stack(i_color))
            grid_data%max_dest_stack(i_color) = maxval(section_data%max_dest_stack(i_color))
        end do
    end subroutine
end module
