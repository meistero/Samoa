! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


!*****************************************************************
! MODULE SFC_edge_traversal: initialize edge data structures from an SFC_data_types
!*****************************************************************
#include "Compilation_control.f90"

#include "Tools_grid.f90"

module SFC_edge_traversal
	use SFC_data_types
    use Grid
    use Int64_list

	implicit none

	public

    interface
        function op_edge_merge(local_edge, neighbor_edge) result(l_conform)
            import
            type(t_edge_data), intent(inout)    :: local_edge
            type(t_edge_data), intent(in)       :: neighbor_edge
            logical                             :: l_conform
        end function

        function op_node_merge(local_node, neighbor_node) result(l_conform)
            import
            type(t_node_data), intent(inout)    :: local_node
            type(t_node_data), intent(in)       :: neighbor_node
            logical                             :: l_conform
        end function
    end interface

	contains

	!*******************
	!create a new grid and cut it into sections of uniform load
	!*******************
	subroutine create_destination_grid(src_grid, dest_grid)
		type(t_grid), intent(inout)                             :: src_grid
		type(t_grid), intent(inout)           	                :: dest_grid

		integer (kind = GRID_DI), parameter                     :: min_section_size = 4_GRID_DI
		type(t_grid_section), pointer                           :: section
		type(t_section_info_list), save                         :: section_descs
		integer (kind = GRID_DI)                                :: i_grid_load, i_total_load, i_grid_partial_load, i_src_section_partial_load, i_section_partial_load, i_section_partial_load_prev
		integer (kind = GRID_DI)                                :: i_src_section, i_dest_section, i_src_sections, i_dest_sections, i_sum_cells, i_sum_cells_prev, i_total_section
		integer (kind = GRID_DI)                                :: i_total_sections, i_eff_dest_cells

        _log_write(4, '(3X, A)') "create splitting"

        i_src_sections = src_grid%sections%get_size()

        !$omp single
        call reduce(src_grid%load, src_grid%sections%elements_alloc%load, MPI_SUM, .false.)
        call prefix_sum(src_grid%sections%elements_alloc%last_dest_cell, src_grid%sections%elements_alloc%dest_cells)

        i_grid_load = src_grid%load

        if (cfg%l_split_sections) then
            i_grid_partial_load = i_grid_load
            call prefix_sum(i_grid_partial_load, MPI_SUM)
            i_total_load = i_grid_load
            call reduce(i_total_load, MPI_SUM)

            i_total_sections = int(cfg%i_sections_per_thread, GRID_DI) * int(cfg%i_threads, GRID_DI) * int(size_MPI, GRID_DI)
        else
            i_grid_partial_load = i_grid_load
            i_total_load = max(1_DI, i_grid_load)

            i_total_sections = int(cfg%i_sections_per_thread, GRID_DI) * int(cfg%i_threads, GRID_DI)
        end if

        i_dest_sections = (i_total_sections * i_grid_partial_load + i_total_load - 1) / i_total_load &
            - (i_total_sections * (i_grid_partial_load - i_grid_load)) / i_total_load

        if (src_grid%dest_cells < i_dest_sections * min_section_size) then
            i_dest_sections = min(src_grid%dest_cells / min_section_size, i_dest_sections)

            if (src_grid%dest_cells > 0) then
                i_dest_sections = max(i_dest_sections, 1)
            end if

            _log_write(4, '(3X, A, I0)') "sections: ", i_dest_sections
            call section_descs%resize(int(i_dest_sections, GRID_SI))

            i_sum_cells_prev = 0
            i_eff_dest_cells = src_grid%dest_cells - i_dest_sections * min_section_size
            assert(i_eff_dest_cells .ge. 0 .or. i_dest_sections == 1)

            do i_dest_section = 1, i_dest_sections
                !Set the number of cells to the difference of partial sums. This guarantees, that there are exactly src_grid%dest_cells cells in total.

                !Add +4 to create enough space for additional refinement cells
                !64bit arithmetics are needed here, (i_dest_section * src_grid%dest_cells) can become very big!
                i_sum_cells = (i_eff_dest_cells * i_dest_section) / i_dest_sections
                section_descs%elements(i_dest_section)%i_cells = min_section_size + i_sum_cells - i_sum_cells_prev + 4_GRID_DI
                section_descs%elements(i_dest_section)%i_load = (i_grid_load * i_dest_section) / i_dest_sections - (i_grid_load * (i_dest_section - 1)) / i_dest_sections
                assert(section_descs%elements(i_dest_section)%i_cells - 4_GRID_DI .ge. min_section_size .or. i_dest_sections == 1)
                i_sum_cells_prev = i_sum_cells

                section_descs%elements(i_dest_section)%index = i_dest_section
                section_descs%elements(i_dest_section)%i_stack_nodes = src_grid%max_dest_stack - src_grid%min_dest_stack + 1
                section_descs%elements(i_dest_section)%i_stack_edges = src_grid%max_dest_stack - src_grid%min_dest_stack + 1

                call section_descs%elements(i_dest_section)%estimate_bounds()
                _log_write(4, '(4X, "Cells: ", I0)') section_descs%elements(i_dest_section)%i_cells - 4_GRID_DI
                _log_write(4, '(4X, "Load: ", I0)') section_descs%elements(i_dest_section)%i_load
            end do

            if (i_dest_sections > 0) then
                assert_eq(i_sum_cells, i_eff_dest_cells)
            end if
        else
            _log_write(4, '(3X, A, I0)') "sections: ", i_dest_sections
            call section_descs%resize(int(i_dest_sections, GRID_SI))

            i_sum_cells = 0
            i_sum_cells_prev = 0
            i_section_partial_load_prev = i_grid_partial_load - i_grid_load
            i_eff_dest_cells = src_grid%dest_cells - i_dest_sections * min_section_size
            assert_ge(i_eff_dest_cells, 0)
            i_total_section = (i_total_sections * (i_grid_partial_load - i_grid_load)) / i_total_load + 1_GRID_DI

            i_src_section = 1
            i_src_section_partial_load = i_grid_partial_load - i_grid_load

            do i_dest_section = 1, i_dest_sections
                i_section_partial_load = min(1_GRID_DI + (i_total_load * i_total_section - 1_GRID_DI) / i_total_sections, i_grid_partial_load)
                assert_gt(i_section_partial_load, i_grid_partial_load - i_grid_load)
                assert_lt((i_total_load * (i_total_section - 1_GRID_DI)) / i_total_sections, i_grid_partial_load)

                do
                    if (i_src_section_partial_load + src_grid%sections%elements_alloc(i_src_section)%load >= i_section_partial_load) then
                        exit
                    endif

                    i_src_section_partial_load = i_src_section_partial_load + src_grid%sections%elements_alloc(i_src_section)%load
                    i_src_section = i_src_section + 1
                end do

                i_sum_cells = src_grid%sections%elements_alloc(i_src_section)%last_dest_cell - src_grid%sections%elements_alloc(i_src_section)%dest_cells &
                 + ((i_section_partial_load - i_src_section_partial_load) * src_grid%sections%elements_alloc(i_src_section)%dest_cells) / &
                    src_grid%sections%elements_alloc(i_src_section)%load

                section_descs%elements(i_dest_section)%i_load = i_section_partial_load - i_section_partial_load_prev

                section_descs%elements(i_dest_section)%i_cells = min_section_size &
                     + (i_sum_cells * i_eff_dest_cells) / src_grid%dest_cells &
                    - (i_sum_cells_prev * i_eff_dest_cells) / src_grid%dest_cells + 4_GRID_DI

                assert_ge(section_descs%elements(i_dest_section)%i_cells - 4_GRID_DI, min_section_size)

                section_descs%elements(i_dest_section)%index = i_dest_section
                section_descs%elements(i_dest_section)%i_stack_nodes = src_grid%max_dest_stack - src_grid%min_dest_stack + 1
                section_descs%elements(i_dest_section)%i_stack_edges = src_grid%max_dest_stack - src_grid%min_dest_stack + 1

                call section_descs%elements(i_dest_section)%estimate_bounds()
                _log_write(4, '(4X, "Cells: ", I0, " Load: ", I0)') section_descs%elements(i_dest_section)%i_cells - 4_GRID_DI, section_descs%elements(i_dest_section)%i_load

                i_sum_cells_prev = i_sum_cells
                i_section_partial_load_prev = i_section_partial_load
                i_total_section = i_total_section + 1
            end do

            if (i_dest_sections > 0) then
                assert_eq(i_total_section, 1 + (i_total_sections * i_grid_partial_load + i_total_load - 1_GRID_DI) / i_total_load)
                assert_eq(i_sum_cells, src_grid%dest_cells)
            end if
        end if
        !$omp end single copyprivate(i_dest_sections)

        assert_le(i_dest_sections, ishft(1, 14) - 1)

        !$omp single
            dest_grid%stats = src_grid%stats
            dest_grid%t_global_data = src_grid%t_global_data
        !$omp end single

        !create new grid
        call dest_grid%create(section_descs, src_grid%max_dest_stack - src_grid%min_dest_stack + 1)

        dest_grid%threads%elements(i_thread)%stats = src_grid%threads%elements(i_thread)%stats
    end subroutine

    subroutine update_distances(grid)
        type(t_grid), intent(inout)						:: grid
        integer (kind = GRID_SI)                        :: i_section, i_first_local_section, i_last_local_section, i_color
        type(t_grid_section), pointer					:: section
        integer (kind = GRID_DI)                        :: min_distances(RED : GREEN)

        _log_write(3, '(3X, A)') "update distances:"

        _log_write(4, '(4X, A)') "grid distances (before):"
        _log_write(4, '(5X, A, 2(F0.4, X))') "start:", decode_distance(grid%start_distance)
        _log_write(4, '(5X, A, 2(F0.4, X))') "min  :", decode_distance(grid%min_distance)
        _log_write(4, '(5X, A, 2(F0.4, X))') "end  :", decode_distance(grid%end_distance)

        call grid%get_local_sections(i_first_local_section, i_last_local_section)

        _log_write(4, '(4X, A)') "section distances (not matched):"
        do i_section = i_first_local_section, i_last_local_section
            section => grid%sections%elements_alloc(i_section)
            _log_write(4, '(5X, A, I0)') "section:", i_section

            _log_write(4, '(6X, A, 2(F0.4, X))') "start:", decode_distance(section%start_distance)
            _log_write(4, '(6X, A, 2(F0.4, X))') "min  :", decode_distance(section%min_distance)
            _log_write(4, '(6X, A, 2(F0.4, X))') "end  :", decode_distance(section%end_distance)

            section%min_distance = section%min_distance - section%start_distance
            section%end_distance = section%end_distance - section%start_distance
            section%start_distance = section%end_distance
        end do

        !$omp barrier

        !$omp single
        if (grid%sections%get_size() > 0) then
            call prefix_sum(grid%sections%elements%end_distance(RED), grid%sections%elements%end_distance(RED))
            call prefix_sum(grid%sections%elements%end_distance(GREEN), grid%sections%elements%end_distance(GREEN))
        end if
        !$omp end single

        do i_section = i_first_local_section, i_last_local_section
            section => grid%sections%elements_alloc(i_section)

            section%min_distance = section%min_distance + section%end_distance - section%start_distance
            section%start_distance = section%end_distance - section%start_distance
        end do

        !$omp barrier

        !$omp single
        min_distances(RED) = grid%min_distance(RED) - minval(grid%sections%elements%min_distance(RED))
        min_distances(GREEN) = grid%min_distance(GREEN) - minval(grid%sections%elements%min_distance(GREEN))
        !$omp end single copyprivate(min_distances)

        _log_write(4, '(4X, A)') "section distances (matched):"
        do i_section = i_first_local_section, i_last_local_section
            section => grid%sections%elements_alloc(i_section)

            _log_write(4, '(5X, A, I0)') "section: ", i_section

            section%min_distance = section%min_distance + min_distances
            section%start_distance = section%start_distance + min_distances
            section%end_distance = section%end_distance + min_distances

            do i_color = RED, GREEN
                section%boundary_edges(i_color)%elements%min_distance = section%boundary_edges(i_color)%elements%min_distance + section%min_distance(i_color)
                section%boundary_nodes(i_color)%elements%distance = section%boundary_nodes(i_color)%elements%distance + section%min_distance(i_color)
            end do

            _log_write(4, '(6X, A, 2(F0.4, X))') "start:", decode_distance(section%start_distance)
            _log_write(4, '(6X, A, 2(F0.4, X))') "min  :", decode_distance(section%min_distance)
            _log_write(4, '(6X, A, 2(F0.4, X))') "end  :", decode_distance(section%end_distance)
        end do

        !$omp barrier

        !$omp single
        if (grid%sections%get_size() > 0) then
            grid%start_distance = grid%sections%elements(1)%start_distance
            grid%end_distance = grid%sections%elements(grid%sections%get_size())%end_distance
        end if

        _log_write(4, '(4X, A)') "grid distances (after):"
        _log_write(4, '(5X, A, 2(F0.4, X))') "start:", decode_distance(grid%start_distance)
        _log_write(4, '(5X, A, 2(F0.4, X))') "min  :", decode_distance(grid%min_distance)
        _log_write(4, '(5X, A, 2(F0.4, X))') "end  :", decode_distance(grid%end_distance)
        !$omp end single
    end subroutine

	!*******************
	!update neighbors
	!*******************

    subroutine update_neighbors(src_grid, dest_grid)
        type(t_grid), intent(in)						:: src_grid
        type(t_grid), intent(inout)						:: dest_grid

        type(t_int64_list), save, allocatable           :: neighbor_min_distances_red(:), neighbor_min_distances_green(:)
        type(t_integer_list), save                      :: src_neighbor_list_red, src_neighbor_list_green
        integer                                         :: i_comm, i_error

        _log_write(4, '(X, A)') "update neighbors:"

        !$omp single
        src_neighbor_list_red = t_integer_list()
        src_neighbor_list_green = t_integer_list()

        !find all source grid neighbors
        call get_grid_neighbors(src_grid, src_neighbor_list_red, int(RED, 1))
        call get_grid_neighbors(src_grid, src_neighbor_list_green, int(GREEN, 1))

        !collect minimum distances
        call collect_minimum_distances(dest_grid, src_neighbor_list_red, neighbor_min_distances_red, int(RED, 1))
        call collect_minimum_distances(dest_grid, src_neighbor_list_green, neighbor_min_distances_green, int(GREEN, 1))
        !$omp end single

        !barrier here, all destination sections must have found their boundary elements in order to continue

        !find all destination grid neighbors by checking if any overlap between local sections and neighbors exists
        call create_dest_neighbor_lists(dest_grid, src_neighbor_list_red, src_neighbor_list_green, neighbor_min_distances_red, neighbor_min_distances_green)

        !$omp single
#		if defined(_MPI)
            do i_comm = 1, size(neighbor_min_distances_red)
                call neighbor_min_distances_red(i_comm)%clear()
            end do

            do i_comm = 1, size(neighbor_min_distances_green)
                call neighbor_min_distances_green(i_comm)%clear()
            end do

        	deallocate(neighbor_min_distances_red, stat = i_error); assert_eq(i_error, 0)
        	deallocate(neighbor_min_distances_green, stat = i_error); assert_eq(i_error, 0)
#		endif

        call src_neighbor_list_red%clear()
        call src_neighbor_list_green%clear()
        !$omp end single
    end subroutine

    subroutine get_grid_neighbors(grid, rank_list, i_color)
        type(t_grid), intent(in)						:: grid
        type(t_integer_list), intent(out)               :: rank_list
        integer (kind = BYTE), intent(in)			        :: i_color

        type(t_grid_section), pointer                   :: section
        type(t_comm_interface), pointer                 :: comm
        integer (kind = GRID_SI)			            :: i_section, i_comm

        _log_write(4, '(3X, A, A)') "get neighbors: ", trim(color_to_char(i_color))

        !orient rank list along the traversal direction
        rank_list%forward = grid%sections%forward

        do i_section = 1, grid%sections%get_size()
            section => grid%sections%elements(i_section)

            _log_write(4, '(5X, A, I0)') "section:", i_section

            do i_comm = 1, section%comms(i_color)%get_size()
                comm => section%comms(i_color)%elements(i_comm)

                _log_write(4, '(6X, A, I0)') "comm: ", comm%neighbor_rank

                if (comm%neighbor_rank .ge. 0 .and. comm%neighbor_rank .ne. rank_MPI) then
                    if (rank_list%get_size() .eq. 0) then
                        call rank_list%add(comm%neighbor_rank)
                        _log_write(4, '(7X, A)') "added (because list is empty)"
                    else if (.not. any(rank_list%elements .eq. comm%neighbor_rank)) then
                        call rank_list%add(comm%neighbor_rank)
                        _log_write(4, '(7X, A)') "added (because not in list)"
                    end if
                end if
            end do
        end do

        if (.not. rank_list%is_forward()) then
            call rank_list%reverse()
        end if

        _log_write(4, '(4X, A, I0, 1X, I0)') "#neighbors: ", rank_list%get_size()
    end subroutine

subroutine collect_minimum_distances(grid, rank_list, neighbor_min_distances, i_color)
        type(t_grid), intent(inout)						    :: grid
        type(t_integer_list), intent(inout)                 :: rank_list
        type(t_int64_list), allocatable, intent(inout)      :: neighbor_min_distances(:)
        integer (BYTE), intent(in)				            :: i_color

        integer(kind = GRID_DI)                             :: local_min_distances(grid%sections%get_size())
        integer							                    :: requests(rank_list%get_size(), 2), i_neighbor_sections(rank_list%get_size())
        integer											    :: i_comm, i_neighbors, i_section, i_error, i_sections

        !Collect minimum distances from all sections of all neighbor processes
        _log_write(3, '(3X, A, A)') "collect minimum distances from sections: ", trim(color_to_char(i_color))

#		if defined(_MPI)
            i_sections = grid%sections%get_size()
            i_neighbors = rank_list%get_size()

		    requests = MPI_REQUEST_NULL
            local_min_distances(:) = grid%sections%elements_alloc(:)%min_distance(i_color)

            _log_write(4, '(4X, A, I0, A)') "rank: ", rank_MPI, " (local)"
            do i_section = 1, i_sections
                _log_write(4, '(5X, A, I0, A, F0.4)') "local section ", i_section, " distance: ", decode_distance(local_min_distances(i_section))
            end do

            allocate(neighbor_min_distances(i_neighbors), stat = i_error); assert_eq(i_error, 0)

            if (i_neighbors == 0) then
                return
            end if

            !send/receive distances
            assert_veq(requests, MPI_REQUEST_NULL)

		    do i_comm = 1, i_neighbors
                call mpi_irecv(i_neighbor_sections(i_comm), 1, MPI_INTEGER, rank_list%elements(i_comm), 0, MPI_COMM_WORLD, requests(i_comm, 2), i_error); assert_eq(i_error, 0)
                call mpi_isend(i_sections, 1, MPI_INTEGER, rank_list%elements(i_comm), 0, MPI_COMM_WORLD, requests(i_comm, 1), i_error); assert_eq(i_error, 0)
 		    end do

            call mpi_waitall(i_neighbors, requests(:, 1), MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)
            call mpi_waitall(i_neighbors, requests(:, 2), MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)

		    do i_comm = 1, i_neighbors
                call neighbor_min_distances(i_comm)%resize(i_neighbor_sections(i_comm))
 		    end do

            assert_veq(requests, MPI_REQUEST_NULL)

		    do i_comm = 1, i_neighbors
                call mpi_irecv(neighbor_min_distances(i_comm)%elements(1), i_neighbor_sections(i_comm), MPI_INTEGER8, rank_list%elements(i_comm), 0, MPI_COMM_WORLD, requests(i_comm, 2), i_error); assert_eq(i_error, 0)
                call mpi_isend(local_min_distances(1), i_sections, MPI_INTEGER8, rank_list%elements(i_comm), 0, MPI_COMM_WORLD, requests(i_comm, 1), i_error); assert_eq(i_error, 0)
 		    end do

            call mpi_waitall(i_neighbors, requests(:, 1), MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)
            call mpi_waitall(i_neighbors, requests(:, 2), MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)

            do i_comm = 1, i_neighbors
                _log_write(4, '(4X, A, I0)') "rank: ", rank_list%elements(i_comm)
                do i_section = 1, neighbor_min_distances(i_comm)%get_size()
                    _log_write(4, '(5X, A, I0, A, F0.4)') "section: ", i_section, " , distance: ", decode_distance(neighbor_min_distances(i_comm)%elements(i_section))
                end do
            end do
#		endif
    end subroutine

    subroutine create_dest_neighbor_lists(grid, src_neighbor_list_red, src_neighbor_list_green, neighbor_min_distances_red, neighbor_min_distances_green)
        type(t_grid), intent(inout)						:: grid
        type(t_integer_list), intent(in)                :: src_neighbor_list_red, src_neighbor_list_green
        type(t_int64_list), intent(in)                  :: neighbor_min_distances_red(:), neighbor_min_distances_green(:)

        integer (KIND = GRID_SI)			            :: i_section, i_first_local_section, i_last_local_section
        integer (kind = BYTE)			                :: i_color
        type(t_grid_section), pointer					:: section

        _log_write(4, '(3X, A, I0, X, I0)') "create destination comm list: #neighbors: ", src_neighbor_list_red%get_size(), src_neighbor_list_red%get_size()

        call grid%get_local_sections(i_first_local_section, i_last_local_section)

        do i_section = i_first_local_section, i_last_local_section
            section => grid%sections%elements_alloc(i_section)

            call set_comms_local_data(grid, section, src_neighbor_list_red, neighbor_min_distances_red, int(RED, 1))
            call set_comms_local_data(grid, section, src_neighbor_list_green, neighbor_min_distances_green, int(GREEN, 1))
        end do

        !$omp barrier

        do i_section = i_first_local_section, i_last_local_section
            assert_eq(grid%sections%elements_alloc(i_section)%index, i_section)

            _log_write(4, '(4X, A, I0)') "section set neighbor pointers: ", i_section

            do i_color = RED, GREEN
                call set_comm_neighbor_data(grid, grid%sections%elements_alloc(i_section), i_color)
            end do
        end do
    end subroutine

    subroutine set_comms_local_data(grid, section, src_neighbor_list, neighbor_min_distances, i_color)
        type(t_grid), intent(inout)						:: grid
        type(t_grid_section), pointer, intent(inout)    :: section
        type(t_integer_list), intent(in)                :: src_neighbor_list
        type(t_int64_list), intent(in)                  :: neighbor_min_distances(:)
        integer (kind = BYTE), intent(in)			    :: i_color

        integer                                         :: i_comm

        _log_write(4, '(4X, A, I0)') "set local data: section ", section%index

        call find_section_neighbors(grid, section, src_neighbor_list, neighbor_min_distances, i_color)

        !count number of boundary elements shared with each comm and set the owner of each element
        call count_section_boundary_elements(section, i_color)
        call set_comm_local_pointers(section, i_color)

        !test for correctness
        assert_eq(sum(section%comms(i_color)%elements%i_edges), section%boundary_edges(i_color)%get_size())
        assert_eq(sum(max(section%comms(i_color)%elements%i_nodes - 1, 0)) + 1, section%boundary_nodes(i_color)%get_size())

        _log_write(4, '(5X, A, A, A)') "Neighbors (", trim(color_to_char(i_color)), "):"
        do i_comm = 1, section%comms(i_color)%get_size()
            _log_write(4, '(6X, (A))') trim(section%comms(i_color)%elements(i_comm)%to_string())
        end do
    end subroutine

    subroutine find_section_neighbors(grid, section, src_neighbor_list, neighbor_min_distances, i_color)
        type(t_grid), intent(inout)						:: grid
        type(t_grid_section), pointer , intent(inout)	:: section
        type(t_integer_list), intent(in)                :: src_neighbor_list
        type(t_int64_list), intent(in)                  :: neighbor_min_distances(:)
        integer (kind = BYTE), intent(in)			    :: i_color

        integer (KIND = GRID_SI)						:: i_comm, i_section_2, i_comms_old, i_comms_new
        integer (KIND = GRID_DI)                        :: min_distance, max_distance
        type(t_grid_section), pointer                   :: section_2

        min_distance = section%min_distance(i_color)

        !clear comm list if it is not empty
        assert_eq(section%comms(i_color)%get_size(), 0)
        call section%comms_type(OLD, i_color)%unattach()
        call section%comms_type(NEW, i_color)%unattach()

        !we build the old communicators in forward traversal order
        !and the new communicators in reverse order
        section%comms_type(OLD, i_color)%forward = .true.
        section%comms_type(NEW, i_color)%forward = .false.

        !first, iterate over old boundary
        max_distance = max(section%start_distance(i_color), section%end_distance(i_color))

        !first, check local sections
        do i_section_2 = section%index - 1, 1, -1
            section_2 => grid%sections%elements_alloc(i_section_2)

            if (max_distance < min_distance) then
                exit
            end if

            assert_eq(i_section_2, section_2%index)

            if (section_2%min_distance(i_color) .le. max_distance) then
                call section%comms_type(OLD, i_color)%add(t_comm_interface(local_rank = rank_MPI, neighbor_rank = rank_MPI, local_section = section%index, neighbor_section = i_section_2, min_distance = section_2%min_distance(i_color)))
                max_distance = section_2%min_distance(i_color)
            end if
        end do

        !next, check process neighbors
        do i_comm = 1, src_neighbor_list%get_size()
            if (max_distance < min_distance .or. src_neighbor_list%elements(i_comm) > rank_MPI) then
                exit
            end if

            do i_section_2 = neighbor_min_distances(i_comm)%get_size(), 1, -1
                if (max_distance < min_distance) then
                    exit
                end if

                if (neighbor_min_distances(i_comm)%elements(i_section_2) .le. max_distance) then
                    call section%comms_type(OLD, i_color)%add(t_comm_interface(local_rank = rank_MPI, neighbor_rank = src_neighbor_list%elements(i_comm), local_section = section%index, neighbor_section = i_section_2, min_distance = neighbor_min_distances(i_comm)%elements(i_section_2)))
                    max_distance = neighbor_min_distances(i_comm)%elements(i_section_2)
                end if
            end do
        end do

        !add a domain boundary comm and serialize comm list for faster access
        call section%comms_type(OLD, i_color)%add(t_comm_interface(local_rank = rank_MPI, neighbor_rank = -1, local_section = section%index, neighbor_section = 0, min_distance = -huge(1_GRID_DI)))

        !second, iterate over new boundary
        max_distance = max(section%start_distance(i_color), section%end_distance(i_color))

        !first, check local sections
        do i_section_2 = section%index + 1, grid%sections%get_size()
            section_2 => grid%sections%elements_alloc(i_section_2)

            if (max_distance < min_distance) then
                exit
            end if

            assert_eq(i_section_2, section_2%index)

            if (section_2%min_distance(i_color) .le. max_distance) then
                call section%comms_type(NEW, i_color)%add(t_comm_interface(local_rank = rank_MPI, neighbor_rank = rank_MPI, local_section = section%index, neighbor_section = i_section_2, min_distance = section_2%min_distance(i_color)))
                max_distance = section_2%min_distance(i_color)
            end if
        end do

        !next, check process neighbors
        do i_comm = src_neighbor_list%get_size(), 1, -1
            if (max_distance < min_distance .or. src_neighbor_list%elements(i_comm) < rank_MPI) then
                exit
            end if

            do i_section_2 = 1, neighbor_min_distances(i_comm)%get_size()
                if (max_distance < min_distance) then
                    exit
                end if

                if (neighbor_min_distances(i_comm)%elements(i_section_2) .le. max_distance) then
                    call section%comms_type(NEW, i_color)%add(t_comm_interface(local_rank = rank_MPI, neighbor_rank = src_neighbor_list%elements(i_comm), local_section = section%index, neighbor_section = i_section_2, min_distance = neighbor_min_distances(i_comm)%elements(i_section_2)))
                    max_distance = neighbor_min_distances(i_comm)%elements(i_section_2)
                end if
            end do
        end do

        !add a domain boundary comm and serialize comm list for faster access
        call section%comms_type(NEW, i_color)%add(t_comm_interface(local_rank = rank_MPI, neighbor_rank = -1, local_section = section%index, neighbor_section = 0, min_distance = -huge(1_GRID_DI)))

        !merge old and new comm lists (new list must be reversed first)
        call section%comms_type(NEW, i_color)%reverse()
        i_comms_old = section%comms_type(OLD, i_color)%get_size()
        i_comms_new = section%comms_type(NEW, i_color)%get_size()
        call section%comms(i_color)%clear()
        section%comms(i_color) = section%comms(i_color)%merge(section%comms_type(OLD, i_color), section%comms_type(NEW, i_color))

        call section%comms_type(OLD, i_color)%clear()
        call section%comms_type(NEW, i_color)%clear()

        if (grid%sections%is_forward()) then
            call section%comms_type(OLD, i_color)%attach(section%comms(i_color)%elements(1 : i_comms_old), .true.)
            call section%comms_type(NEW, i_color)%attach(section%comms(i_color)%elements(i_comms_old + 1 : ), .true.)
        else
            call section%comms(i_color)%reverse()
            call section%comms_type(OLD, i_color)%attach(section%comms(i_color)%elements(1 : i_comms_new), .false.)
            call section%comms_type(NEW, i_color)%attach(section%comms(i_color)%elements(i_comms_new + 1 : ), .false.)
        end if
    end subroutine

    subroutine count_section_boundary_elements(section, i_color)
        type(t_grid_section), intent(inout)				:: section
        integer (kind = BYTE), intent(in)			        :: i_color

        integer (KIND = GRID_SI)                        :: i_current_edge, i_current_node, i_pass, i_comm
        type(t_edge_data), pointer                      :: current_edge
        type(t_node_data), pointer                      :: current_node
        type(t_comm_interface), pointer                 :: comm

        !match section boundaries with comm distances to find the number of edges and nodes
        !that each comm shares with each section

        _log_write(4, '(5X, A, I0)') "count section boundary elements: section ", section%index

        assert_eqv(section%comms(i_color)%is_forward(), section%cells%is_forward())
        assert_eqv(section%comms_type(OLD, i_color)%is_forward(), section%cells%is_forward())
        assert_eqv(section%comms_type(NEW, i_color)%is_forward(), section%cells%is_forward())

        !reverse new neighbors in order to compare decreasing distances
        call section%boundary_type_edges(NEW, i_color)%reverse()
        call section%boundary_type_nodes(NEW, i_color)%reverse()
        call section%comms_type(NEW, i_color)%reverse()

        do i_pass = OLD, NEW
            _log_write(4, '(6X, A, A)') "pass: ", edge_type_to_char(i_pass)
            _log_write(4, '(7X, A)') "compare edges:"

            assert_eqv(section%comms_type(i_pass, i_color)%is_forward(), section%boundary_type_nodes(i_pass, i_color)%is_forward())
            assert_eqv(section%comms_type(i_pass, i_color)%is_forward(), section%boundary_type_edges(i_pass, i_color)%is_forward())

            i_comm = 1
            assert_le(i_comm, section%comms_type(i_pass, i_color)%get_size())
            comm => section%comms_type(i_pass, i_color)%elements(i_comm)

            _log_write(4, '(8X, A, A)') "nb: ", trim(comm%to_string())

            do i_current_edge = 1, section%boundary_type_edges(i_pass, i_color)%get_size()
                current_edge => section%boundary_type_edges(i_pass, i_color)%elements(i_current_edge)

                do while (current_edge%min_distance .lt. comm%min_distance)
                    i_comm = i_comm + 1
                    assert_le(i_comm, section%comms_type(i_pass, i_color)%get_size())
                    comm => section%comms_type(i_pass, i_color)%elements(i_comm)

                    _log_write(4, '(8X, A, A)') "nb: ", trim(comm%to_string())
                end do

                comm%i_edges = comm%i_edges + 1
                _log_write(4, '(8X, A, F0.4, X, F0.4)') "edge: ", decode_distance(current_edge%min_distance), decode_distance(current_edge%min_distance + encode_edge_size(current_edge%depth))
            end do

            _log_write(4, '(7X, A)') "compare nodes:"

            i_comm = 1
            assert_le(i_comm, section%comms_type(i_pass, i_color)%get_size())
            comm => section%comms_type(i_pass, i_color)%elements(i_comm)

            _log_write(4, '(8X, A, A)') "nb  : ", trim(comm%to_string())

            do i_current_node = 1, section%boundary_type_nodes(i_pass, i_color)%get_size()
                current_node => section%boundary_type_nodes(i_pass, i_color)%elements(i_current_node)

                do while (current_node%distance .lt. comm%min_distance)
                    i_comm = i_comm + 1
                    assert_le(i_comm, section%comms_type(i_pass, i_color)%get_size())
                    comm => section%comms_type(i_pass, i_color)%elements(i_comm)

                    _log_write(4, '(8X, A, A)') "nb: ", trim(comm%to_string())
                end do

                do
                    comm%i_nodes = comm%i_nodes + 1
                    _log_write(4, '(8X, A, F0.4)') "node: ", decode_distance(current_node%distance)

                    if (current_node%distance .eq. comm%min_distance) then
                        i_comm = i_comm + 1
                        assert_le(i_comm, section%comms_type(i_pass, i_color)%get_size())
                        comm => section%comms_type(i_pass, i_color)%elements(i_comm)

                        _log_write(4, '(8X, A, A)') "nb  : ", trim(comm%to_string())
                    else
                        exit
                    end if
                end do
            end do
        end do

        !restore the correct order by reversing again.
        call section%boundary_type_edges(NEW, i_color)%reverse()
        call section%boundary_type_nodes(NEW, i_color)%reverse()
        call section%comms_type(NEW, i_color)%reverse()
    end subroutine

    subroutine set_comm_local_pointers(section, i_color)
        type(t_grid_section), intent(inout)				:: section
        integer (kind = BYTE), intent(in)			        :: i_color

        integer, parameter                              :: divisor = 4

        integer (KIND = GRID_SI)                        :: i_comm
        integer (KIND = GRID_SI)                        :: i_first_edge, i_first_node, i_last_edge, i_last_node
        logical                                         :: l_forward
        type(t_comm_interface), pointer                 :: comm

        _log_write(4, '(4X, A, I0)') "set local pointers: section ", section%index
        _log_write(4, '(4X, A, A, A)') "Neighbors (", trim(color_to_char(i_color)), "):"

        !set ownership to true initally and falsify each wrong entry afterwards
        section%boundary_edges(i_color)%elements%owned_locally = .true.
        section%boundary_edges(i_color)%elements%owned_globally = .true.
        section%boundary_nodes(i_color)%elements%owned_locally = .true.
        section%boundary_nodes(i_color)%elements%owned_globally = .true.

        i_first_edge = 0
        i_last_edge = 0
        i_first_node = 1
        i_last_node = 0

        l_forward = section%boundary_nodes(i_color)%is_forward()

        assert_eqv(section%comms(i_color)%is_forward(), section%boundary_edges(i_color)%is_forward())
        assert_eqv(section%comms(i_color)%is_forward(), section%boundary_nodes(i_color)%is_forward())

        do i_comm = 1, section%comms(i_color)%get_size()
            comm => section%comms(i_color)%elements(i_comm)

            _log_write(4, '(6X, (A))') trim(comm%to_string())

            assert_eq(comm%local_rank, rank_MPI)
            assert_eq(comm%local_section, section%index)

            i_first_edge = i_last_edge + 1
            i_first_node = max(i_first_node, i_last_node)
            i_last_edge = i_first_edge + comm%i_edges - 1
            i_last_node = i_first_node + comm%i_nodes - 1

            if (l_forward) then
                comm%p_local_edges => section%boundary_edges(i_color)%elements(i_first_edge : i_last_edge)
                comm%p_local_nodes => section%boundary_nodes(i_color)%elements(i_first_node : i_last_node)
            else
                comm%p_local_edges => section%boundary_edges(i_color)%elements(i_last_edge : i_first_edge : -1)
                comm%p_local_nodes => section%boundary_nodes(i_color)%elements(i_last_node : i_first_node : -1)
            end if

            !each entity with a neighbor of lower section index is not owned by the current section
            if (comm%neighbor_rank .ge. 0 .and. comm%neighbor_rank < rank_MPI) then
                comm%p_local_edges%owned_globally = .false.
                comm%p_local_nodes%owned_globally = .false.
            else if (comm%neighbor_rank == rank_MPI) then
                if (comm%neighbor_section < section%index) then
                    comm%p_local_edges%owned_locally = .false.
                    comm%p_local_edges%owned_globally = .false.
                    comm%p_local_nodes((comm%i_nodes - 1)/divisor + 2 : comm%i_nodes - 1 - (comm%i_nodes - 1)/divisor)%owned_locally = .false.
                    comm%p_local_nodes((comm%i_nodes - 1)/divisor + 2 : comm%i_nodes - 1 - (comm%i_nodes - 1)/divisor)%owned_globally = .false.
                else
                    comm%p_local_nodes(1:(comm%i_nodes - 1)/divisor + 1)%owned_locally = .false.
                    comm%p_local_nodes(1:(comm%i_nodes - 1)/divisor + 1)%owned_globally = .false.
                    comm%p_local_nodes(comm%i_nodes - (comm%i_nodes - 1)/divisor : comm%i_nodes)%owned_locally = .false.
                    comm%p_local_nodes(comm%i_nodes - (comm%i_nodes - 1)/divisor : comm%i_nodes)%owned_globally = .false.
                end if
            end if

            assert_eq(size(comm%p_local_edges), comm%i_edges)
            assert_eq(size(comm%p_local_nodes), comm%i_nodes)
        end do
    end subroutine

	!>finds the correct comm index for a certain rank and section
    function find_section_comm(comms, i_rank, i_section) result(comm)
        type(t_comm_interface), pointer, intent(in)         :: comms(:)
        type(t_comm_interface), pointer                     :: comm
        integer (kind = GRID_SI), intent(in)                :: i_rank
        integer (kind = GRID_SI), intent(in)                :: i_section

        integer (kind = GRID_SI)                            :: i_comm

        !find the correct comm index

        i_comm = minloc(abs(comms%neighbor_rank - i_rank) + abs(comms%neighbor_section - i_section), 1)
        comm => comms(i_comm)

        assert_eq(comm%neighbor_rank, i_rank)
        assert_eq(comm%neighbor_section, i_section)
    end function

    subroutine set_comm_neighbor_data(grid, section, i_color)
        type(t_grid), intent(inout)				        :: grid
        type(t_grid_section),target, intent(inout)	    :: section
        integer (kind = BYTE), intent(in)			    :: i_color

        integer (KIND = GRID_SI)                        :: i_comm, i_pass
        type(t_comm_interface), pointer                 :: comm, comm_2
        type(t_grid_section), pointer				    :: section_2

        _log_write(4, '(4X, A, I0)') "set neighbor data: section ", section%index
        _log_write(4, '(4X, A, A, A)') "Neighbors (", trim(color_to_char(i_color)), "):"

        do i_pass = OLD, NEW
            _log_write(4, '(5X, A, A)') trim(edge_type_to_char(i_pass)), ":"

            do i_comm = 1, section%comms_type(i_pass, i_color)%get_size()
                comm => section%comms_type(i_pass, i_color)%elements(i_comm)
                assert_eq(section%index, comm%local_section)

                _log_write(4, '(6X, (A))') trim(comm%to_string())

                assert_eq(comm%local_rank, rank_MPI)
                assert_eq(comm%local_section, section%index)

                if (comm%neighbor_rank == rank_MPI) then
                    section_2 => grid%sections%elements_alloc(comm%neighbor_section)
                    assert_eq(section_2%index, comm%neighbor_section)

                    comm_2 => find_section_comm(section_2%comms_type(1 - i_pass, i_color)%elements, rank_MPI, section%index)

                    assert_eq(comm%local_rank, comm_2%neighbor_rank)
                    assert_eq(comm%local_section, comm_2%neighbor_section)
                    assert_eq(comm_2%local_rank, comm%neighbor_rank)
                    assert_eq(comm_2%local_section, comm%neighbor_section)

                    comm%p_neighbor_edges => comm_2%p_local_edges
                    comm%p_neighbor_nodes => comm_2%p_local_nodes

                    assert_eq(comm%i_edges, comm_2%i_edges)
                    assert_eq(comm%i_nodes, comm_2%i_nodes)
                    assert_veq(decode_distance(comm%p_local_edges%min_distance), decode_distance(comm_2%p_local_edges(comm%i_edges : 1 : -1)%min_distance))
                    assert_veq(decode_distance(comm%p_local_nodes%distance), decode_distance(comm_2%p_local_nodes(comm%i_nodes : 1 : -1)%distance))
                    assert_veq(comm%p_local_nodes%position(1), comm_2%p_local_nodes(comm%i_nodes : 1 : -1)%position(1))
                    assert_veq(comm%p_local_nodes%position(2), comm_2%p_local_nodes(comm%i_nodes : 1 : -1)%position(2))
                else if (comm%neighbor_rank .ge. 0) then
                    call comm%create_buffer()

                    assert_eq(size(comm%p_neighbor_edges), comm%i_edges)
                    assert_eq(size(comm%p_neighbor_nodes), comm%i_nodes)
                end if
            end do
        end do
    end subroutine

    subroutine delete_comm_neighbor_data(section, i_color)
        type(t_grid_section),target, intent(inout)	    :: section
        integer (kind = BYTE), intent(in)			    :: i_color

        integer (KIND = GRID_SI)                        :: i_comm
        type(t_comm_interface), pointer                 :: comm

        _log_write(4, '(4X, A, I0)') "set neighbor data: section ", section%index
        _log_write(4, '(4X, A, A, A)') "Neighbors (", trim(color_to_char(i_color)), "):"

        do i_comm = 1, section%comms(i_color)%get_size()
            comm => section%comms(i_color)%elements(i_comm)
            assert_eq(section%index, comm%local_section)

            _log_write(4, '(6X, (A))') trim(comm%to_string())

            call comm%destroy_buffer()
        end do
    end subroutine

    subroutine send_mpi_boundary(section, mpi_edge_type_optional, mpi_node_type_optional)
        type(t_grid_section), intent(inout)				:: section
        integer, intent(in), optional                   :: mpi_edge_type_optional, mpi_node_type_optional

        integer (kind = GRID_SI)						:: i_comm
        integer                                         :: i_error, send_tag, recv_tag
        integer                                         :: mpi_edge_type, mpi_node_type, mpi_edge_extent, mpi_node_extent, mpi_edge_size, mpi_node_size
        integer (BYTE)							        :: i_color
        type(t_comm_interface), pointer			        :: comm

#        if defined(_MPI)
            if (present(mpi_edge_type_optional)) then
                mpi_edge_type = mpi_edge_type_optional
                call MPI_Type_size(mpi_edge_type, mpi_edge_size, i_error); assert_eq(i_error, 0)
                mpi_edge_extent = 1
            else
                mpi_edge_type = MPI_BYTE
                mpi_edge_size = 1
                mpi_edge_extent = sizeof(section%boundary_edges(RED)%elements(1))
            end if

            if (present(mpi_node_type_optional)) then
                mpi_node_type = mpi_node_type_optional
                call MPI_Type_size(mpi_node_type, mpi_node_size, i_error); assert_eq(i_error, 0)
                mpi_node_extent = 1
            else
                mpi_node_type = MPI_BYTE
                mpi_node_size = 1
                mpi_node_extent = sizeof(section%boundary_nodes(RED)%elements(1))
            end if

            _log_write(4, '(4X, A, I0)') "send mpi boundary: section ", section%index

            do i_color = RED, GREEN
                _log_write(4, '(5X, A, A)') trim(color_to_char(i_color)), ":"

                do i_comm = 1, section%comms(i_color)%get_size()
                    comm => section%comms(i_color)%elements(i_comm)

                    assert(associated(comm%p_local_edges))
                    assert(associated(comm%p_local_nodes))

                    if (comm%neighbor_rank .ge. 0 .and. comm%neighbor_rank .ne. rank_MPI) then
                        _log_write(4, '(6X, A)') trim(comm%to_string())
                        assert(associated(comm%p_neighbor_edges))
                        assert(associated(comm%p_neighbor_nodes))
						assert_lt(comm%local_section, ishft(1, 14))
						assert_lt(comm%neighbor_section, ishft(1, 14))

                        send_tag = ishft(comm%local_section, 14) + comm%neighbor_section
                        recv_tag = ishft(comm%neighbor_section, 14) + comm%local_section

                        assert_lt(send_tag, ishft(1, 28))
                        assert_lt(recv_tag, ishft(1, 28))

                        _log_write(5, '(7X, A, I0, X, I0, A, I0, X, I0, A, I0)') "send from: ", comm%local_rank, comm%local_section,  " to  : ", comm%neighbor_rank, comm%neighbor_section, " send tag: ", send_tag

                        if (comm%i_edges * mpi_edge_size > 0) then
                            assert_eq(comm%send_requests(1), MPI_REQUEST_NULL)
                            call mpi_isend(get_c_pointer(comm%p_local_edges), comm%i_edges * mpi_edge_extent, mpi_edge_type, comm%neighbor_rank, send_tag, MPI_COMM_WORLD, comm%send_requests(1), i_error); assert_eq(i_error, 0)
                            assert_ne(comm%send_requests(1), MPI_REQUEST_NULL)
                        end if

                        if (comm%i_nodes * mpi_node_size > 0) then
                            assert_eq(comm%send_requests(2), MPI_REQUEST_NULL)
                            call mpi_isend(get_c_pointer(comm%p_local_nodes), comm%i_nodes * mpi_node_extent, mpi_node_type, comm%neighbor_rank, send_tag, MPI_COMM_WORLD, comm%send_requests(2), i_error); assert_eq(i_error, 0)
                            assert_ne(comm%send_requests(2), MPI_REQUEST_NULL)
                        endif
                    end if
                end do
            end do
#       endif
    end subroutine

    subroutine recv_mpi_boundary(section, mpi_edge_type_optional, mpi_node_type_optional)
        type(t_grid_section), intent(inout)				:: section
        integer, intent(in), optional                   :: mpi_edge_type_optional, mpi_node_type_optional

        integer (kind = GRID_SI)						:: i_comm
        integer                                         :: i_error, send_tag, recv_tag
        integer                                         :: mpi_edge_type, mpi_node_type, mpi_edge_extent, mpi_node_extent, mpi_edge_size, mpi_node_size
        integer (BYTE)							        :: i_color
        type(t_comm_interface), pointer			        :: comm

#        if defined(_MPI)
            if (present(mpi_edge_type_optional)) then
                mpi_edge_type = mpi_edge_type_optional
                call MPI_Type_size(mpi_edge_type, mpi_edge_size, i_error); assert_eq(i_error, 0)
                mpi_edge_extent = 1
            else
                mpi_edge_type = MPI_BYTE
                mpi_edge_size = 1
                mpi_edge_extent = sizeof(section%boundary_edges(RED)%elements(1))
            end if

            if (present(mpi_node_type_optional)) then
                mpi_node_type = mpi_node_type_optional
                call MPI_Type_size(mpi_node_type, mpi_node_size, i_error); assert_eq(i_error, 0)
                mpi_node_extent = 1
            else
                mpi_node_type = MPI_BYTE
                mpi_node_size = 1
                mpi_node_extent = sizeof(section%boundary_nodes(RED)%elements(1))
            end if

            _log_write(4, '(4X, A, I0)') "recv mpi boundary: section ", section%index

             do i_color = RED, GREEN
                _log_write(4, '(5X, A, A)') trim(color_to_char(i_color)), ":"

                do i_comm = 1, section%comms(i_color)%get_size()
                    comm => section%comms(i_color)%elements(i_comm)

                    assert(associated(comm%p_local_edges))
                    assert(associated(comm%p_local_nodes))

                    if (comm%neighbor_rank .ge. 0 .and. comm%neighbor_rank .ne. rank_MPI) then
                        _log_write(4, '(6X, A)') trim(comm%to_string())
                        assert(associated(comm%p_neighbor_edges))
                        assert(associated(comm%p_neighbor_nodes))
						assert_lt(comm%local_section, ishft(1, 14))
						assert_lt(comm%neighbor_section, ishft(1, 14))

                        send_tag = ishft(comm%local_section, 14) + comm%neighbor_section
                        recv_tag = ishft(comm%neighbor_section, 14) + comm%local_section

                        assert_lt(send_tag, ishft(1, 28))
                        assert_lt(recv_tag, ishft(1, 28))

                        _log_write(5, '(7X, A, I0, X, I0, A, I0, X, I0, A, I0)') "recv to  : ", comm%local_rank, comm%local_section, " from: ", comm%neighbor_rank, comm%neighbor_section, " recv tag: ", recv_tag

                        if (comm%i_edges * mpi_edge_size > 0) then
                            assert_eq(comm%recv_requests(1), MPI_REQUEST_NULL)
                            call mpi_irecv(get_c_pointer(comm%p_neighbor_edges), comm%i_edges * mpi_edge_extent, mpi_edge_type, comm%neighbor_rank, recv_tag, MPI_COMM_WORLD, comm%recv_requests(1), i_error); assert_eq(i_error, 0)
                            assert_ne(comm%recv_requests(1), MPI_REQUEST_NULL)
                        end if

                        if (comm%i_nodes * mpi_node_size > 0) then
                            assert_eq(comm%recv_requests(2), MPI_REQUEST_NULL)
                            call mpi_irecv(get_c_pointer(comm%p_neighbor_nodes), comm%i_nodes * mpi_node_extent, mpi_node_type, comm%neighbor_rank, recv_tag, MPI_COMM_WORLD, comm%recv_requests(2), i_error); assert_eq(i_error, 0)
                            assert_ne(comm%recv_requests(2), MPI_REQUEST_NULL)
                        endif
                    end if
                end do
            end do
#       endif
    end subroutine

    subroutine finish_send(section, i_color, mpi_edge_size, mpi_node_size)
        type(t_grid_section), intent(inout)			    :: section
        integer (BYTE), intent(in)					    :: i_color
        integer, intent(in)                             :: mpi_edge_size, mpi_node_size

#       if defined(_MPI)
            integer (kind = GRID_SI)				    :: i_comm
            integer                                     :: i_error
            type(t_comm_interface), pointer			    :: comm

            do i_comm = 1, section%comms(i_color)%get_size()
                comm => section%comms(i_color)%elements(i_comm)

                if (comm%neighbor_rank .ne. rank_MPI .and. comm%neighbor_rank .ge. 0) then
                    _log_write(4, '(7X, (A))') trim(comm%to_string())

                    assert(associated(comm%p_local_edges))
                    assert(associated(comm%p_local_nodes))
                    assert(associated(comm%p_neighbor_edges))
                    assert(associated(comm%p_neighbor_nodes))

                    _log_write(5, '(8X, A, I0, X, I0, A, I0, X, I0)') "send from : ", comm%local_rank, comm%local_section, " to  : ", comm%neighbor_rank, comm%neighbor_section

                    if (comm%i_edges * mpi_edge_size > 0) then
                        assert_ne(comm%send_requests(1), MPI_REQUEST_NULL)
                        call mpi_wait(comm%send_requests(1), MPI_STATUS_IGNORE, i_error); assert_eq(i_error, 0)
                        assert_eq(comm%send_requests(1), MPI_REQUEST_NULL)
                    end if

                    if (comm%i_nodes * mpi_node_size > 0) then
                        assert_ne(comm%send_requests(2), MPI_REQUEST_NULL)
                        call mpi_wait(comm%send_requests(2), MPI_STATUS_IGNORE, i_error); assert_eq(i_error, 0)
                        assert_eq(comm%send_requests(2), MPI_REQUEST_NULL)
                    end if
                end if
            end do
#       endif
    end subroutine

    function merge_comm(comm, edge_merge_op, node_merge_op) result(l_comm_conform)
        type(t_comm_interface), intent(inout)			:: comm
        procedure(op_edge_merge)                        :: edge_merge_op
        procedure(op_node_merge)                        :: node_merge_op
        logical                                         :: l_comm_conform, l_conform

        integer                                         :: i_error, i
        integer (kind = GRID_SI)       					:: i_first_node, i_last_node

        l_comm_conform = .true.

        !gather data in section with lowest index

        assert(associated(comm%p_local_edges))
        assert(associated(comm%p_local_nodes))
        assert(associated(comm%p_neighbor_edges))
        assert(associated(comm%p_neighbor_nodes))
        assert_eq(size(comm%p_local_edges), size(comm%p_neighbor_edges))
        assert_eq(size(comm%p_local_nodes), size(comm%p_neighbor_nodes))
        assert_eq(size(comm%p_local_edges), comm%i_edges)
        assert_eq(size(comm%p_local_nodes), comm%i_nodes)

        _log_write(4, '(7X, (A))') trim(comm%to_string())

        !merge on local edges
        !only the owner may execute merge operations (a race condition might occur otherwise)

        assert_veq(decode_distance(comm%p_local_edges%min_distance), decode_distance(comm%p_neighbor_edges(comm%i_edges : 1 : -1)%min_distance))

        do i = 1, comm%i_edges
            if (comm%p_local_edges(i)%owned_locally) then
                assert(comm%neighbor_rank .ne. rank_MPI .or. .not. comm%p_neighbor_edges(comm%i_edges + 1 - i)%owned_locally)

                !Careful: the Fortran standard does neither define if logical expressions are executed nor the order of operations.
                !Hence, under the assumption that edge_merge_op() must always be executed, the code below is incorrect as it may short circuit:
                ! l_section_conform = l_section_conform .and. edge_merge_op(...)
                !This is the correct way:

                l_conform = edge_merge_op(comm%p_local_edges(i), comm%p_neighbor_edges(comm%i_edges + 1 - i))
                l_comm_conform = l_comm_conform .and. l_conform
            end if
        end do

        assert_veq(decode_distance(comm%p_local_nodes%distance), decode_distance(comm%p_neighbor_nodes(comm%i_nodes : 1 : -1)%distance))
        assert_veq(comm%p_local_nodes%position(1), comm%p_neighbor_nodes(comm%i_nodes : 1 : -1)%position(1))
        assert_veq(comm%p_local_nodes%position(2), comm%p_neighbor_nodes(comm%i_nodes : 1 : -1)%position(2))

        !merge on local nodes
        !only the owner may execute merge operations (which is the local section)
        !a race condition can occur otherwise

        do i = 1, comm%i_nodes
            if (comm%p_local_nodes(i)%owned_locally) then
                assert(comm%neighbor_rank .ne. rank_MPI .or. .not. comm%p_neighbor_nodes(comm%i_nodes + 1 - i)%owned_locally)

                l_conform = node_merge_op(comm%p_local_nodes(i), comm%p_neighbor_nodes(comm%i_nodes + 1 - i))
                l_comm_conform = l_comm_conform .and. l_conform
            end if
        end do
    end function

    subroutine finish_receive_and_merge(section, i_color, edge_merge_op, node_merge_op, mpi_edge_size, mpi_node_size)
        type(t_grid_section), intent(inout)			    :: section
        integer (BYTE), intent(in)					    :: i_color
        procedure(op_edge_merge)                        :: edge_merge_op
        procedure(op_node_merge)                        :: node_merge_op
        integer, intent(in)                             :: mpi_edge_size, mpi_node_size

        integer                                         :: i_error
        integer (kind = GRID_SI)						:: i_comm
        logical                                         :: l_conform, l_section_conform
        type(t_comm_interface), pointer			        :: comm

        _log_write(4, '(3X, A)') "sync boundary sections:"

        !receive all boundary data, then merge

        _log_write(4, '(3X, A, I0)') "wait for MPI receives, read and write merged data"

        l_section_conform = .true.

        do i_comm = 1, section%comms(i_color)%get_size()
            comm => section%comms(i_color)%elements(i_comm)
            !proceed with edge and node merging

            !The way we defined ownership implies that the current section cannot own
            !any edges or nodes, if the neighbor section is of lower index.
            !In this case we can skip communication with this particular section.

            if (comm%neighbor_rank .eq. rank_MPI) then
                l_conform = merge_comm(comm, edge_merge_op, node_merge_op)
                l_section_conform = l_section_conform .and. l_conform
            end if
        end do

#       if defined(_MPI)
            do i_comm = 1, section%comms(i_color)%get_size()
                comm => section%comms(i_color)%elements(i_comm)

                if (comm%neighbor_rank .ne. rank_MPI .and. comm%neighbor_rank .ge. 0) then
                    _log_write(4, '(7X, (A))') trim(comm%to_string())

                    assert(associated(comm%p_local_edges))
                    assert(associated(comm%p_local_nodes))
                    assert(associated(comm%p_neighbor_edges))
                    assert(associated(comm%p_neighbor_nodes))
                    assert_eq(size(comm%p_local_edges), size(comm%p_neighbor_edges))
                    assert_eq(size(comm%p_local_nodes), size(comm%p_neighbor_nodes))

                    _log_write(5, '(8X, A, I0, X, I0, A, I0, X, I0)') "receive at: ", comm%local_rank, comm%local_section, " from: ", comm%neighbor_rank, comm%neighbor_section

                    if (comm%i_edges * mpi_edge_size > 0) then
                        assert_ne(comm%recv_requests(1), MPI_REQUEST_NULL)
                        call mpi_wait(comm%recv_requests(1), MPI_STATUS_IGNORE, i_error); assert_eq(i_error, 0)
                        assert_eq(comm%recv_requests(1), MPI_REQUEST_NULL)
                    end if

                    if (comm%i_nodes * mpi_node_size > 0) then
                        assert_ne(comm%recv_requests(2), MPI_REQUEST_NULL)
                        call mpi_wait(comm%recv_requests(2), MPI_STATUS_IGNORE, i_error); assert_eq(i_error, 0)
                        assert_eq(comm%recv_requests(2), MPI_REQUEST_NULL)
                    end if

                    !proceed with edge and node merging

                    l_conform = merge_comm(comm, edge_merge_op, node_merge_op)
                    l_section_conform = l_section_conform .and. l_conform
                end if
            end do
#       endif

        if (.not. l_section_conform) then
            !Careful, this is concurrent write access, but should not cause a race condition
            section%l_conform = .false.
        end if
    end subroutine

    subroutine duplicate_section_boundary_data(grid, section, i_color, edge_write_op, node_write_op)
        type(t_grid), intent(inout)			            :: grid
        type(t_grid_section), intent(inout)			    :: section
        integer (BYTE), intent(in)					    :: i_color
        procedure(op_edge_merge)                        :: edge_write_op
        procedure(op_node_merge)                        :: node_write_op

        integer (kind = GRID_SI)						:: i_comm
        integer                                         :: i
        type(t_comm_interface), pointer			        :: comm
        integer (kind = GRID_SI)       					:: i_first_node, i_last_node
        logical                                         :: l_conform, l_section_conform

        do i_comm = 1, section%comms(i_color)%get_size()
            comm => section%comms(i_color)%elements(i_comm)

            if (comm%neighbor_rank .eq. rank_MPI) then
                _log_write(4, '(7X, (A))') trim(comm%to_string())

                l_section_conform = .true.

                !write to neighbor edges
                !only the owner may execute write operations, as he merged all the data

                assert_veq(decode_distance(comm%p_local_edges%min_distance), decode_distance(comm%p_neighbor_edges(comm%i_edges : 1 : -1)%min_distance))

                do i = 1, comm%i_edges
                    if (comm%p_local_edges(i)%owned_locally) then
                        assert(.not. comm%p_neighbor_edges(comm%i_edges + 1 - i)%owned_locally)

                        l_conform = edge_write_op(comm%p_neighbor_edges(comm%i_edges + 1 - i), comm%p_local_edges(i))
                        l_section_conform = l_section_conform .and. l_conform
                    end if
                end do

                !write to neighbor nodes
                !only the owner may execute write operations, as he merged all the data

                assert_veq(decode_distance(comm%p_local_nodes%distance), decode_distance(comm%p_neighbor_nodes(comm%i_nodes : 1 : -1)%distance))
                assert_veq(comm%p_local_nodes%position(1), comm%p_neighbor_nodes(comm%i_nodes : 1 : -1)%position(1))
                assert_veq(comm%p_local_nodes%position(2), comm%p_neighbor_nodes(comm%i_nodes : 1 : -1)%position(2))

                do i = 1, comm%i_nodes
                    if (comm%p_local_nodes(i)%owned_locally) then
                        assert(.not. comm%p_neighbor_nodes(comm%i_nodes + 1 - i)%owned_locally)

                        l_conform = node_write_op(comm%p_neighbor_nodes(comm%i_nodes + 1 - i), comm%p_local_nodes(i))
                        l_section_conform = l_section_conform .and. l_conform
                    end if
                end do

                if (.not. l_section_conform) then
                    !Careful, this is concurrent write access, but should not cause a race condition
                    grid%sections%elements_alloc(comm%neighbor_section)%l_conform = .false.
                end if
            end if
        end do
    end subroutine

    subroutine collect_boundary_data(grid, edge_merge_op, node_merge_op, mpi_edge_type_optional, mpi_node_type_optional)
        type(t_grid), intent(inout)						:: grid
        procedure(op_edge_merge)                        :: edge_merge_op
        procedure(op_node_merge)                        :: node_merge_op
        integer, intent(in), optional                   :: mpi_edge_type_optional, mpi_node_type_optional

        integer (kind = GRID_SI)						:: i_section, i_first_local_section, i_last_local_section, i_comm
        integer                                         :: i_error
        integer (BYTE)							        :: i_color
        integer                                         :: mpi_edge_size, mpi_node_size

        _log_write(4, '(3X, A)') "sync boundary sections:"

        call grid%get_local_sections(i_first_local_section, i_last_local_section)

#       if defined(_MPI)
            if (present(mpi_edge_type_optional)) then
                call MPI_Type_size(mpi_edge_type_optional, mpi_edge_size, i_error); assert_eq(i_error, 0)
            else
                mpi_edge_size = 1
            end if

            if (present(mpi_node_type_optional)) then
                call MPI_Type_size(mpi_node_type_optional, mpi_node_size, i_error); assert_eq(i_error, 0)
            else
                mpi_node_size = 1
            end if
#       endif

        _log_write(4, '(6X, "wait for MPI sends/receives and merge:")')

        do i_section = i_first_local_section, i_last_local_section
            _log_write(4, '(4X, A, I0)') "section: ", i_section

            do i_color = RED, GREEN
                _log_write(4, '(5X, A, A)') trim(color_to_char(i_color)), ":"

                !send all boundary data
                !This must be done for all communicators of the section before merging the boundary, as boundary data may overlap.
                !If one communicator started overwriting a vertex when another communicator has not finished sending,
                !the vertex may be overwritten before it is sent.

                call finish_send(grid%sections%elements_alloc(i_section), i_color, mpi_edge_size, mpi_node_size)

                !receive all boundary data, then merge
                call finish_receive_and_merge(grid%sections%elements_alloc(i_section), i_color, edge_merge_op, node_merge_op, mpi_edge_size, mpi_node_size)
            end do
        end do
    end subroutine

    subroutine duplicate_boundary_data(grid, edge_write_op, node_write_op)
        type(t_grid), intent(inout)						:: grid
        procedure(op_edge_merge)                        :: edge_write_op
        procedure(op_node_merge)                        :: node_write_op

        integer (kind = GRID_SI)						:: i_section, i_first_local_section, i_last_local_section
        integer (BYTE)							        :: i_color

        _log_write(4, '(3X, A)') "duplicate boundary data from owner sections:"

        call grid%get_local_sections(i_first_local_section, i_last_local_section)

        do i_section = i_first_local_section, i_last_local_section
            assert_eq(i_section, grid%sections%elements_alloc(i_section)%index)

            _log_write(4, '(4X, A, I0)') "section: ", i_section

            do i_color = RED, GREEN
                _log_write(4, '(5X, A, A)') trim(color_to_char(i_color)), ":"

                call duplicate_section_boundary_data(grid, grid%sections%elements_alloc(i_section), i_color, edge_write_op, node_write_op)
            end do
        end do
    end subroutine

    subroutine distribute_load(grid, r_max_imbalance)
        type(t_grid), intent(inout)						:: grid
		real, intent(in)               					:: r_max_imbalance		!< maximum allowed global imbalance (i.e. 0.1 = 10%)

        integer (kind = GRID_DI)                        :: i_total_load, i_max_load
        integer (kind = GRID_SI)						:: i_first_local_section, i_last_local_section, i_section
        integer, pointer                                :: i_rank_out(:), i_section_index_out(:), i_rank_in(:)
        type(t_grid), save						        :: grid_temp
        integer	(BYTE)  		                        :: i_color
        logical                                         :: l_early_exit

#		if defined(_MPI)
            _log_write(3, '(3X, A)') "distribute load"

            !$omp barrier

		    !$omp single
			!compute global imbalance by a prefix sum over the grid load

		    call prefix_sum(grid%sections%elements_alloc%partial_load, grid%sections%elements_alloc%load)
		    call reduce(grid%load, grid%sections%elements_alloc%load, MPI_SUM, .false.)

			!check imbalance
			if (r_max_imbalance > 0.0d0) then
                i_total_load = grid%load
                i_max_load = grid%load
                call reduce(i_total_load, MPI_SUM)
                call reduce(i_max_load, MPI_MAX)
            else
                i_max_load = 2_GRID_DI
                i_total_load = 1_GRID_DI
            end if
	        !$omp end single copyprivate(i_max_load, i_total_load)

			!exit early if the imbalance is small enough
	        if (dble(i_max_load) * size_MPI .le. (1.0d0 + r_max_imbalance) * i_total_load) then
                _log_write(2, '(4X, "load balancing: imbalance above threshold? no: ", F0.3, " < ", F0.3)') dble(i_max_load) * size_MPI / dble(i_total_load) - 1.0d0, r_max_imbalance
                return
           	end if

	        _log_write(2, '(4X, "load balancing: imbalance above threshold? yes: ", F0.3, " > ", F0.3)') dble(i_max_load) * size_MPI / dble(i_total_load) - 1.0d0, r_max_imbalance

            if (cfg%l_serial_lb) then
                call compute_partition_serial(grid, i_rank_out, i_section_index_out, i_rank_in, l_early_exit)
            else
                call compute_partition_distributed(grid, i_rank_out, i_section_index_out, i_rank_in, l_early_exit)
            end if

            if (l_early_exit) then
                return
            end if

			_log_write(3, '(4X, "send comm changes:")')

	        _log_write(4, '(4X, A)') "grid distances (before):"
	        _log_write(4, '(5X, A, 2(F0.4, X))') "start:", decode_distance(grid%start_distance)
	        _log_write(4, '(5X, A, 2(F0.4, X))') "min  :", decode_distance(grid%min_distance)
	        _log_write(4, '(5X, A, 2(F0.4, X))') "end  :", decode_distance(grid%end_distance)

            !HACK: this allows neighbor representations to be transferred from old to new grid
	        call grid%get_local_sections(i_first_local_section, i_last_local_section)

	        do i_section = i_first_local_section, i_last_local_section
	            assert_eq(grid%sections%elements_alloc(i_section)%index, i_section)
                call save_neighbor_rep(grid%sections%elements_alloc(i_section))
	        end do

            !Rank and section indices may have changed, so send the new information to the neighbors
            call send_recv_comm_changes(grid, i_rank_out, i_section_index_out, i_rank_in)

	        !$omp barrier

			_log_write(3, '(4X, "migrate sections:")')

	        !exit early if nothing changes on this rank

	        if (size(i_rank_out) == 0 .and. size(i_rank_in) == 0) then
	            !if we have nothing and we get nothing from others we are done
                _log_write(2, '(4X, "load balancing: exit, we do not have nor get anything ")')
	            return
            else if (size(i_rank_in) == 0) then
                !if we get nothing we cannot give anything to ourselves
                assert_vne(i_rank_out, rank_MPI)
            else if (size(i_rank_out) == 0) then
                !if we have nothing we cannot get anything from ourselves
                assert_vne(i_rank_in, rank_MPI)
            else if (i_rank_out(1) .eq. rank_MPI .and. i_rank_out(size(i_rank_out)) .eq. rank_MPI .and. &
                    i_rank_in(1) .eq. rank_MPI .and. i_rank_in(size(i_rank_in)) .eq. rank_MPI) then

                !if we give nothing to others and get nothing from others we are done
                _log_write(2, '(4X, "load balancing: exit, we have something, but do not give nor get anything")')
                return
	        end if

            !Delete all buffers
	        do i_section = i_first_local_section, i_last_local_section
	            do i_color = RED, GREEN
                    call delete_comm_neighbor_data(grid%sections%elements_alloc(i_section), i_color)
	            end do
	        end do

	        !$omp single
	        if (grid%sections%get_size() > 0) then
	            assert_veq(decode_distance(grid%start_distance), decode_distance(grid%sections%elements(1)%start_distance))
	            assert_veq(decode_distance(grid%end_distance), decode_distance(grid%sections%elements(grid%sections%get_size())%end_distance))
			end if

            grid_temp%sections = t_grid_section_list()
            call grid_temp%sections%resize(size(i_rank_in))
            !$omp end single

            !if necessary, reverse order to match source and destination grid
            if (.not. grid%sections%is_forward()) then
                !$omp barrier
                _log_write(4, '(X, A)') "Reverse destination grid.."
                call grid_temp%reverse()
                !$omp barrier
            end if

			call send_recv_section_infos(grid, grid_temp, i_rank_out, i_rank_in)
			call send_recv_section_data(grid, grid_temp, i_rank_out, i_rank_in)

            !$omp barrier

            !$omp single
            call grid%sections%clear()
			grid%sections = grid_temp%sections

	        !update distances again and check for correctness
	        if (grid%sections%get_size() > 0) then
	            grid%start_distance = grid%sections%elements(1)%start_distance
	            grid%end_distance = grid%sections%elements(grid%sections%get_size())%end_distance
	        else
	            grid%start_distance = 0
	            grid%end_distance = 0
	        end if

	        do i_color = RED, GREEN
	            call reduce(grid%min_distance(i_color), grid%sections%elements%min_distance(i_color), MPI_MIN, .false.)
	        end do

	    	call prefix_sum(grid%sections%elements_alloc%partial_load, grid%sections%elements_alloc%load)
	    	call reduce(grid%load, grid%sections%elements_alloc%load, MPI_SUM, .false.)
	        call reduce(grid%dest_cells, grid%sections%elements%dest_cells, MPI_SUM, .false.)

	        _log_write(4, '(4X, A)') "grid distances (after):"
	        _log_write(4, '(5X, A, 2(F0.4, X))') "start:", decode_distance(grid%start_distance)
	        _log_write(4, '(5X, A, 2(F0.4, X))') "min  :", decode_distance(grid%min_distance)
	        _log_write(4, '(5X, A, 2(F0.4, X))') "end  :", decode_distance(grid%end_distance)
			!$omp end single

	        !resize stacks to make sure they are big enough for the new sections
	        call grid%threads%elements(i_thread)%destroy()
	        call grid%threads%elements(i_thread)%create(grid%max_dest_stack - grid%min_dest_stack + 1)

	        call grid%get_local_sections(i_first_local_section, i_last_local_section)

	        do i_section = i_first_local_section, i_last_local_section
	            assert_eq(grid%sections%elements_alloc(i_section)%index, i_section)

	            _log_write(4, '(4X, A, I0)') "section set local pointers: ", i_section

	            do i_color = RED, GREEN
	                call set_comm_local_pointers(grid%sections%elements_alloc(i_section), i_color)
	            end do
	        end do

	        !$omp barrier

	        do i_section = i_first_local_section, i_last_local_section
	            assert_eq(grid%sections%elements_alloc(i_section)%index, i_section)

	            _log_write(4, '(4X, A, I0)') "section set neighbor pointers: ", i_section

	            do i_color = RED, GREEN
	                call set_comm_neighbor_data(grid, grid%sections%elements_alloc(i_section), i_color)
	            end do

                !HACK: this allows neighbor representations to be transferred from old to new grid
                call load_neighbor_rep(grid%sections%elements_alloc(i_section))
	        end do

	        !$omp barrier
#		endif
    end subroutine

#	if defined(_MPI)
		subroutine load_neighbor_rep(section)
		    type(t_grid_section), intent(inout)             :: section
		    integer						                    :: i_section, i_color, i_comm, i_edge

#		    if defined(_USE_SKELETON_OP)
                type(num_cell_rep)                              :: rep
                type(num_cell_update)                           :: update
                type(t_comm_interface), pointer                 :: comm

                do i_color = RED, GREEN
                    do i_comm = 1, section%comms(i_color)%get_size()
                        comm => section%comms(i_color)%elements(i_comm)

                        if (comm%neighbor_rank .ge. 0 .and. comm%neighbor_rank .ne. rank_MPI) then
                            do i_edge = 1, comm%i_edges
                                update%flux(:) = comm%p_local_edges(i_edge)%update%flux(:)
                                rep = transfer(update, rep)
                                comm%p_neighbor_edges(i_edge)%rep%Q(:) = rep%Q(:)
                            end do
                        end if
                    end do
                end do
#		    endif
		end subroutine

        subroutine save_neighbor_rep(section)
		    type(t_grid_section), intent(inout)             :: section
		    integer						                    :: i_section, i_color, i_comm, i_edge

#		    if defined(_USE_SKELETON_OP)
                type(num_cell_rep)                              :: rep
                type(num_cell_update)                           :: update
                type(t_comm_interface), pointer                 :: comm

                do i_color = RED, GREEN
                    do i_comm = 1, section%comms(i_color)%get_size()
                        comm => section%comms(i_color)%elements(i_comm)

                        if (comm%neighbor_rank .ge. 0) then
                            do i_edge = 1, comm%i_edges
                                rep%Q(:) = comm%p_neighbor_edges(i_edge)%rep%Q(:)
                                update = transfer(rep, update)
                                comm%p_local_edges(i_edge)%update%flux(:) = update%flux(:)
                            end do
                        end if
                    end do
                end do
#		    endif
		end subroutine

        !> uses a distributed algorithm to compute the new partitiion
        !> this approach scales well, but requires local methods and cannot achieve optimal load balance
        subroutine compute_partition_distributed(grid, i_rank_out_local, i_section_index_out_local, i_rank_in_local, l_early_exit)
		    type(t_grid), intent(inout)		                :: grid
		    integer, pointer, intent(inout)                 :: i_rank_out_local(:), i_section_index_out_local(:), i_rank_in_local(:)
		    logical, intent(out)                            :: l_early_exit

            integer, pointer, save                          :: i_rank_out(:) => null(), i_section_index_out(:) => null(), i_rank_in(:) => null()
            integer, allocatable, save                      :: i_sections_out(:), i_sections_in(:), i_partial_sections_in(:), i_partial_sections_out(:), i_delta_out(:), requests_out(:, :), requests_in(:)
            integer (kind = GRID_SI)						:: i_first_local_section, i_last_local_section, i_rank, i_section, i_comm
            integer									        :: i_first_rank_out, i_last_rank_out, i_first_rank_in, i_last_rank_in
            type(t_grid_section), pointer					:: section
            integer						                    :: i_error, i_sections, requests(2)
            integer	(BYTE)  		                        :: i_color
            integer (kind = GRID_DI)					    :: load, partial_load, total_load, rank_imbalance
            character (len=100)                             :: msg

	        l_early_exit = .false.
        	call grid%get_local_sections(i_first_local_section, i_last_local_section)

            !$omp single
			!switch to integer arithmetics from now on, we need exact arithmetics

			load = grid%load
            partial_load = load
			call prefix_sum(partial_load, MPI_SUM)
			total_load = load
			call reduce(total_load, MPI_SUM)

			assert_gt(total_load, 0)
			i_sections = grid%sections%get_size()

            if (i_sections > 0) then
                rank_imbalance = \
                    ((2_GRID_DI * partial_load - grid%sections%elements_alloc(i_sections)%load) * size_MPI) / (2_GRID_DI * total_load) - \
                    ((2_GRID_DI * (partial_load - load) + grid%sections%elements_alloc(1)%load) * size_MPI) / (2_GRID_DI * total_load)
            else
                rank_imbalance = 0
            end if

			call reduce(rank_imbalance, MPI_MAX)

			!allocate arrays on first call
			if (.not. associated(i_rank_out)) then
                allocate(i_rank_out(0), stat=i_error); assert_eq(i_error, 0)
                allocate(i_section_index_out(0), stat=i_error); assert_eq(i_error, 0)
            end if

            if (.not. allocated(i_sections_out)) then
                allocate(i_sections_out(0), stat=i_error); assert_eq(i_error, 0)
                allocate(i_partial_sections_out(0), stat=i_error); assert_eq(i_error, 0)
                allocate(i_delta_out(0), stat=i_error); assert_eq(i_error, 0)
                allocate(requests_out(0, 0), stat=i_error); assert_eq(i_error, 0)
            end if

            if (.not. allocated(i_sections_in)) then
                allocate(i_sections_in(0), stat=i_error); assert_eq(i_error, 0)
                allocate(i_partial_sections_in(0), stat=i_error); assert_eq(i_error, 0)
                allocate(requests_in(0), stat=i_error, errmsg=msg); if (i_error .ne. 0) then; print *, msg; end if; assert_eq(i_error, 0)
            end if

	        !$omp end single copyprivate(rank_imbalance, i_sections, load, partial_load, total_load)

	        !exit early if the imbalance cannot be improved
	        if (rank_imbalance .le. 0) then
                _log_write(2, '(4X, "load balancing: imbalance is improvable? no.")')
                l_early_exit = .true.
                return
	        end if

	        _log_write(2, '(4X, "load balancing: imbalance is improvable? yes.")')

			!$omp single
            i_first_rank_out = ((partial_load - load) * size_MPI) / total_load	!round down
            i_last_rank_out = (partial_load * size_MPI - 1) / total_load		!round up and subtract 1

            !special case: integers are rounded up below 0, so if partial_load = 0, the result is 0, but should be -1
            if (partial_load == 0) then
                i_last_rank_out = -1
            end if

			!allocate space for variables that are stored per output section
			if (size(i_rank_out) .ne. i_sections) then
                deallocate(i_rank_out, stat = i_error); assert_eq(i_error, 0)
                deallocate(i_section_index_out, stat = i_error); assert_eq(i_error, 0)

                allocate(i_rank_out(i_sections), stat=i_error); assert_eq(i_error, 0)
                allocate(i_section_index_out(i_sections), stat=i_error); assert_eq(i_error, 0)
            end if

			!allocate space for variables that are stored per output rank
            if (lbound(i_sections_out, 1) .ne. i_first_rank_out .or. ubound(i_sections_out, 1) .ne. i_last_rank_out) then
                deallocate(i_sections_out, stat=i_error); assert_eq(i_error, 0)
                deallocate(i_partial_sections_out, stat=i_error); assert_eq(i_error, 0)
                deallocate(i_delta_out, stat=i_error); assert_eq(i_error, 0)
                deallocate(requests_out, stat=i_error); assert_eq(i_error, 0)

                allocate(i_sections_out(i_first_rank_out : i_last_rank_out), stat=i_error); assert_eq(i_error, 0)
                allocate(i_partial_sections_out(i_first_rank_out : i_last_rank_out), stat=i_error); assert_eq(i_error, 0)
                allocate(i_delta_out(i_first_rank_out : i_last_rank_out), stat=i_error); assert_eq(i_error, 0)
                allocate(requests_out(i_first_rank_out : i_last_rank_out, 2), stat=i_error); assert_eq(i_error, 0)
            end if

			i_sections_out = 0
			requests_out = MPI_REQUEST_NULL
			requests = MPI_REQUEST_NULL
			!$omp end single copyprivate(i_first_rank_out, i_last_rank_out)

			_log_write(3, '(4X, "match ranks: out ranks ", I0, " to ", I0)') i_first_rank_out, i_last_rank_out

	        !compute the number of sections, that are sent to each rank
	        do i_section = i_first_local_section, i_last_local_section
	            section => grid%sections%elements_alloc(i_section)
	            assert_eq(section%index, i_section)

	            !assign this section to its respective ideal rank
	            i_rank = (((2_GRID_DI * (partial_load - load + section%partial_load) - section%load)) * size_MPI) / (2_GRID_DI * total_load)
	        	assert_ge(i_rank, i_first_rank_out)
	        	assert_le(i_rank, i_last_rank_out)

	        	i_rank_out(i_section) = i_rank

			 	!$omp atomic
                i_sections_out(i_rank) = i_sections_out(i_rank) + 1
	        end do

			!$omp barrier

	        !$omp single
	        assert_eq(sum(i_sections_out), i_sections)
			!receive 2 messages, containing the rank number of the first and the last source rank (they must be unique)

			call mpi_irecv(i_first_rank_in, 1, MPI_INTEGER, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, requests(1), i_error); assert_eq(i_error, 0)
			call mpi_irecv(i_last_rank_in, 1, MPI_INTEGER, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, requests(2), i_error); assert_eq(i_error, 0)

			!My rank must be the first input rank for all my output ranks except the first one
			!and the last input rank for all my output ranks except the last one

            do i_rank = i_first_rank_out + 1, i_last_rank_out
                call mpi_isend(rank_MPI, 1, MPI_INTEGER, i_rank, 1, MPI_COMM_WORLD, requests_out(i_rank, 1), i_error); assert_eq(i_error, 0)
                call mpi_isend(rank_MPI, 1, MPI_INTEGER, i_rank - 1, 2, MPI_COMM_WORLD, requests_out(i_rank - 1, 2), i_error); assert_eq(i_error, 0)

                _log_write(3, '("send: from: ", I0, " to: ", I0, " I am your first input rank")') rank_MPI, i_rank
                _log_write(3, '("send: from: ", I0, " to: ", I0, " I am your last input rank")') rank_MPI, i_rank - 1
            end do

			!check if i am first or last input rank for first and last output rank (meaning an exact match)
			if (load > 0 .and. (partial_load - load) * size_MPI == total_load * i_first_rank_out) then
                call mpi_isend(rank_MPI, 1, MPI_INTEGER, i_first_rank_out, 1, MPI_COMM_WORLD, requests_out(i_first_rank_out, 1), i_error); assert_eq(i_error, 0)
               _log_write(3, '("send: from: ", I0, " to: ", I0, " I am your first input rank")') rank_MPI, i_first_rank_out
			end if

			if (load > 0 .and. partial_load * size_MPI == total_load * (i_last_rank_out + 1)) then
                call mpi_isend(rank_MPI, 1, MPI_INTEGER, i_last_rank_out, 2, MPI_COMM_WORLD, requests_out(i_last_rank_out, 2), i_error); assert_eq(i_error, 0)
                _log_write(3, '("send: from: ", I0, " to: ", I0, " I am your last input rank")') rank_MPI, i_last_rank_out
			end if

			call mpi_waitall(i_last_rank_out - i_first_rank_out + 1, requests_out(i_first_rank_out:i_last_rank_out, 1), MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)
			call mpi_waitall(i_last_rank_out - i_first_rank_out + 1, requests_out(i_first_rank_out:i_last_rank_out, 2), MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)
			call mpi_waitall(2, requests, MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)

			_log_write(3, '(4X, "count sections: out ranks ", I0, " to ", I0, " in ranks: ", I0, " to ", I0)') i_first_rank_out, i_last_rank_out, i_first_rank_in, i_last_rank_in

			!now that I know which ranks I get data from, allocate an array that stores the number of sections I get from each source rank
            if (lbound(i_sections_in, 1) .ne. i_first_rank_in .or. ubound(i_sections_in, 1) .ne. i_last_rank_in) then
                deallocate(i_sections_in, stat=i_error); assert_eq(i_error, 0)
                deallocate(i_partial_sections_in, stat=i_error); assert_eq(i_error, 0)
                deallocate(requests_in, stat=i_error, errmsg=msg); if (i_error .ne. 0) then; print *, msg; end if; assert_eq(i_error, 0)

                allocate(i_sections_in(i_first_rank_in : i_last_rank_in), stat=i_error); assert_eq(i_error, 0)
                allocate(i_partial_sections_in(i_first_rank_in : i_last_rank_in), stat=i_error); assert_eq(i_error, 0)
                allocate(requests_in(i_first_rank_in : i_last_rank_in), stat=i_error, errmsg=msg); if (i_error .ne. 0) then; print *, msg; end if; assert_eq(i_error, 0)
            end if

			i_sections_in = 0
			requests_in = MPI_REQUEST_NULL
			!$omp end single copyprivate(i_first_rank_in, i_last_rank_in)

	        !$omp single
 		    !Communicate the number of sections sent to the destination ranks / received from the source ranks

			_log_write(3, '(4X, "send number of outgoing sections and receive number of incoming sections per rank")')

	        assert_veq(requests_in, MPI_REQUEST_NULL)
	        assert_veq(requests_out, MPI_REQUEST_NULL)

			do i_rank = i_first_rank_in, i_last_rank_in
				call mpi_irecv(i_sections_in(i_rank), 1, MPI_INTEGER, i_rank, 0, MPI_COMM_WORLD, requests_in(i_rank), i_error); assert_eq(i_error, 0)
        	end do

			do i_rank = i_first_rank_out, i_last_rank_out
                call mpi_isend(i_sections_out(i_rank), 1, MPI_INTEGER, i_rank, 0, MPI_COMM_WORLD, requests_out(i_rank, 1), i_error); assert_eq(i_error, 0)
        	end do

			call mpi_waitall(i_last_rank_out - i_first_rank_out + 1, requests_out(i_first_rank_out:i_last_rank_out, 1), MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)
 	    	call mpi_waitall(i_last_rank_in - i_first_rank_in + 1, requests_in, MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)

	        assert_veq(requests_out, MPI_REQUEST_NULL)
			assert_veq(requests_in, MPI_REQUEST_NULL)

			_log_write(3, '(4X, "communicate new index of first incoming section per rank")')

			!compute prefix sum over number of input sections and output sections to find out which rank gets which section indices assigned
	        call prefix_sum(i_partial_sections_in, i_sections_in)
	        call prefix_sum(i_partial_sections_out, i_sections_out)

			!communicate this info again
			do i_rank = i_first_rank_out, i_last_rank_out
                call mpi_irecv(i_delta_out(i_rank), 1, MPI_INTEGER, i_rank, 0, MPI_COMM_WORLD, requests_out(i_rank, 1), i_error); assert_eq(i_error, 0)
        	end do

			do i_rank = i_first_rank_in, i_last_rank_in
                call mpi_isend(i_partial_sections_in(i_rank), 1, MPI_INTEGER, i_rank, 0, MPI_COMM_WORLD, requests_in(i_rank), i_error); assert_eq(i_error, 0)
        	end do

            if (.not. associated(i_rank_in) .or. size(i_rank_in) .ne. sum(i_sections_in)) then
                if (associated(i_rank_in)) then
                    deallocate(i_rank_in, stat=i_error); assert_eq(i_error, 0)
                end if

                allocate(i_rank_in(sum(i_sections_in)), stat=i_error); assert_eq(i_error, 0)
            end if

            i_section = 0
            do i_rank = i_first_rank_in, i_last_rank_in
                if (i_sections_in(i_rank) > 0) then
                    i_rank_in(i_section + 1 : i_section + i_sections_in(i_rank)) = i_rank
                    i_section = i_section + i_sections_in(i_rank)
                end if
            end do

			call mpi_waitall(i_last_rank_out - i_first_rank_out + 1, requests_out(i_first_rank_out:i_last_rank_out, 1), MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)
	    	call mpi_waitall(i_last_rank_in - i_first_rank_in + 1, requests_in, MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)

			_log_write(3, '(4X, "compute new indices of all outgoing sections")')

        	i_delta_out = i_delta_out - i_partial_sections_out
	        !$omp end single

        	!compute the new section indices
	        do i_section = i_first_local_section, i_last_local_section
	            i_section_index_out(i_section) = i_delta_out(i_rank_out(i_section)) + i_section
	        end do

            !pass private copies of the array pointers back to caller function
            i_rank_out_local => i_rank_out
            i_section_index_out_local => i_section_index_out
	        i_rank_in_local => i_rank_in
        end subroutine

        subroutine midpoint_cutoff(l, s)
            integer(kind = GRID_DI), intent(in)     :: l(:)
            integer(kind = GRID_SI), intent(inout)  :: s(:)

            integer                                 :: n, i, j
            integer (kind = GRID_DI)                :: L_j, L_n

            n = size(l)
            L_n = sum(l)
            L_j = 0_GRID_DI

            do j = 1, n
                i = int((size_MPI * (2_GRID_DI * L_j + l(j))) / (2_GRID_DI * L_n))
                s(j) = i

                L_j = L_j + l(j)
            end do
        end subroutine

        subroutine iterative_binary(l, s)
            integer(kind = GRID_DI), intent(in)     :: l(:)
            integer(kind = GRID_SI), intent(inout)  :: s(:)

            integer(kind = GRID_SI), allocatable    :: s_test(:)
            integer                                 :: n, i_error, i, j
            integer (kind = GRID_DI)                :: test, test_max, current_min, current_max, load_i

            n = size(l)

            allocate(s_test(n), stat=i_error); assert_eq(i_error, 0)

            current_max = 0_GRID_DI
            i = 0
            load_i = 0_GRID_DI

            do j = 1, n
                if (s(j) > i) then
                    current_max = max(current_max, load_i)
                    i = s(j)
                    load_i = 0_GRID_DI
                end if

                load_i = load_i + l(j)
            end do

            current_max = max(current_max, load_i) + 1_GRID_DI
            current_min = max(sum(l) / size_MPI, maxval(l))

            test = (current_min + current_max) / 2

            do
                if (current_max <= current_min) then
                    exit
                end if

                test_max = 0_GRID_DI
                i = 0
                load_i = 0_GRID_DI

                do j = 1, n
                    if (load_i + l(j) >= test .and. i < size_MPI - 1) then
                        test_max = max(test_max, load_i)
                        i = i + 1
                        load_i = 0_GRID_DI
                    end if

                    load_i = load_i + l(j)
                    s_test(j) = i
                end do

                test_max = max(test_max, load_i)

                if (test_max < current_max) then
                    current_max = test_max
                    test = (current_min + current_max) / 2

                    s(:) = s_test(:)
                else
                    current_min = test
                    test = current_max
                end if
            end do

            deallocate(s_test, stat=i_error); assert_eq(i_error, 0)
        end subroutine

        !> uses a serial algorithm to compute the new partition
        !> this approach scales badly, but can apply global methods and can achieve optimal load balance
        subroutine compute_partition_serial(grid, i_rank_out_local, i_section_index_out_local, i_rank_in_local, l_early_exit)
		    type(t_grid), intent(inout)		                :: grid
		    integer, pointer, intent(inout)                 :: i_rank_out_local(:), i_section_index_out_local(:), i_rank_in_local(:)
		    logical, intent(out)                            :: l_early_exit

            integer, pointer, save                          :: i_rank_out(:) => null(), i_section_index_out(:) => null(), i_rank_in(:) => null()
            integer, allocatable, save                      :: all_sections(:), displacements(:), all_ranks(:), all_section_indices_out(:)
            integer, allocatable, save                      :: requests_in(:), requests_out(:)
            integer (kind = GRID_DI), allocatable, save     :: all_load(:), local_load(:)
            integer                                         :: total_sections, i_sections_out, i_sections_in, i_error, i, j, k

			l_early_exit = .false.

            !$omp single

            !gather section indices

            if (rank_MPI == 0) then
                allocate(all_sections(0 : size_MPI - 1), stat=i_error); assert_eq(i_error, 0)
                allocate(displacements(0 : size_MPI - 1), stat=i_error); assert_eq(i_error, 0)
            else
                allocate(all_sections(0), stat=i_error); assert_eq(i_error, 0)
                allocate(displacements(0), stat=i_error); assert_eq(i_error, 0)
            end if

            i_sections_out = grid%sections%get_size()

            call mpi_gather(i_sections_out, 1, MPI_INTEGER, all_sections, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)

            !gather load

            if (rank_MPI == 0) then
                call prefix_sum(displacements, all_sections)
                total_sections = displacements(size_MPI - 1)
                displacements(:) = displacements(:) - all_sections(:)

                allocate(all_load(total_sections), stat=i_error); assert_eq(i_error, 0)
                allocate(all_ranks(total_sections), stat=i_error); assert_eq(i_error, 0)
                allocate(all_section_indices_out(total_sections), stat=i_error); assert_eq(i_error, 0)
            else
                allocate(all_load(0), stat=i_error); assert_eq(i_error, 0)
                allocate(all_ranks(0), stat=i_error); assert_eq(i_error, 0)
                allocate(all_section_indices_out(0), stat=i_error); assert_eq(i_error, 0)
            end if

            allocate(local_load(i_sections_out), stat=i_error); assert_eq(i_error, 0)

            local_load(:) = grid%sections%elements_alloc(:)%load
            call mpi_gatherv(local_load, i_sections_out, MPI_INTEGER8, all_load, all_sections, displacements, MPI_INTEGER8, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)

            if (rank_MPI == 0) then
                j = 0
                do i = 0, size_MPI - 1
                    do k = 1, all_sections(i)
                        all_ranks(j + k) = i
                    end do

                    j = j + all_sections(i)
                end do

                call iterative_binary(all_load, all_ranks)

                all_sections(:) = 0

                do j = 1, total_sections
                    i = all_ranks(j)
                    all_sections(i) = all_sections(i) + 1
                    all_section_indices_out(j) = all_sections(i)
                end do
            end if

            !scatter new section count

            call mpi_scatter(all_sections, 1, MPI_INTEGER, i_sections_in, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)

            !scatter new ranks

            !allocate space for variables that are stored per input section
			if (.not. associated(i_rank_in) .or. size(i_rank_in) .ne. i_sections_in) then
                if (associated(i_rank_in)) then
                    deallocate(i_rank_in, stat = i_error); assert_eq(i_error, 0)
                    deallocate(requests_in, stat=i_error); assert_eq(i_error, 0)
                end if

                allocate(i_rank_in(i_sections_in), stat=i_error); assert_eq(i_error, 0)
                allocate(requests_in(i_sections_in), stat=i_error); assert_eq(i_error, 0)
            end if

			!allocate space for variables that are stored per output section
			if (.not. associated(i_rank_out) .or. size(i_rank_out) .ne. i_sections_out) then
                if (associated(i_rank_out)) then
                    deallocate(i_rank_out, stat = i_error); assert_eq(i_error, 0)
                    deallocate(i_section_index_out, stat = i_error); assert_eq(i_error, 0)
                    deallocate(requests_out, stat = i_error); assert_eq(i_error, 0)
                end if

                allocate(i_rank_out(i_sections_out), stat=i_error); assert_eq(i_error, 0)
                allocate(i_section_index_out(i_sections_out), stat=i_error); assert_eq(i_error, 0)
                allocate(requests_out(i_sections_out), stat=i_error); assert_eq(i_error, 0)
            end if

            if (rank_MPI == 0) then
                all_sections(0 : size_MPI - 2) = displacements(1 : size_MPI - 1) - displacements(0 : size_MPI - 2)
                all_sections(size_MPI - 1) = total_sections - displacements(size_MPI - 1)
            end if

            call mpi_scatterv(all_ranks, all_sections, displacements, MPI_INTEGER, i_rank_out, i_sections_out, MPI_INTEGER, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
            call mpi_scatterv(all_section_indices_out, all_sections, displacements, MPI_INTEGER, i_section_index_out, i_sections_out, MPI_INTEGER, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)

            deallocate(local_load, stat=i_error); assert_eq(i_error, 0)
            deallocate(all_load, stat=i_error); assert_eq(i_error, 0)
            deallocate(all_sections, stat=i_error); assert_eq(i_error, 0)
            deallocate(displacements, stat=i_error); assert_eq(i_error, 0)
            deallocate(all_ranks, stat=i_error); assert_eq(i_error, 0)
            deallocate(all_section_indices_out, stat=i_error); assert_eq(i_error, 0)

            do j = 1, i_sections_in
                call mpi_irecv(i_rank_in(j), 1, MPI_INTEGER, MPI_ANY_SOURCE, j, MPI_COMM_WORLD, requests_in(j), i_error); assert_eq(i_error, 0)
            end do

            do j = 1, i_sections_out
                call mpi_isend(rank_MPI, 1, MPI_INTEGER, i_rank_out(j), i_section_index_out(j), MPI_COMM_WORLD, requests_out(j), i_error); assert_eq(i_error, 0)
            end do

            call mpi_waitall(i_sections_in, requests_in, MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)
            call mpi_waitall(i_sections_out, requests_out, MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)
            !$omp end single

            !pass private copies of the array pointers back to caller function
            i_rank_out_local => i_rank_out
            i_section_index_out_local => i_section_index_out
            i_rank_in_local => i_rank_in
        end subroutine

        subroutine send_recv_comm_changes(src_grid, i_rank_out, i_section_index_out, i_rank_in)
		    type(t_grid), intent(inout)		                :: src_grid
		    integer, intent(in)                             :: i_rank_out(:), i_section_index_out(:), i_rank_in(:)

            type(t_grid_section), pointer					:: section
            type(t_comm_interface), pointer                 :: comm
            type(t_comm_interface)                 			:: old_comm
            integer						                    :: i_error, send_tag, recv_tag, i_sections, new_rank, new_section, requests(2)
            integer (kind = GRID_SI)						:: i_first_local_section, i_last_local_section, i_section, i_comm
            integer (BYTE)						            :: i_color

	        call src_grid%get_local_sections(i_first_local_section, i_last_local_section)

	        !Since the other processes cannot possibly know if something changed between two neighbors,
	        !we have to send this information to all neighbors

	        do i_section = i_first_local_section, i_last_local_section
	            section => src_grid%sections%elements_alloc(i_section)
	            assert_eq(section%index, i_section)
	            _log_write(4, '(4X, A, I0)') "section send/recv comm changes: ", i_section

				new_rank = i_rank_out(i_section)
				new_section = i_section_index_out(i_section)
	            section%index = new_section

	            do i_color = RED, GREEN
	                _log_write(4, '(5X, A, A)') trim(color_to_char(i_color)), ":"

	                do i_comm = 1, section%comms(i_color)%get_size()
	                    comm => section%comms(i_color)%elements(i_comm)
	                    _log_write(4, '(6X, A)') trim(comm%to_string())

						old_comm = comm

	                    comm%local_rank = new_rank
	                    comm%local_section = new_section

	                    if (comm%neighbor_rank .ge. 0) then
							assert_lt(old_comm%local_section, ishft(1, 15))
							assert_lt(old_comm%neighbor_section, ishft(1, 15))
					        send_tag = ishft(old_comm%local_section, 15) + old_comm%neighbor_section
					        recv_tag = ishft(old_comm%neighbor_section, 15) + old_comm%local_section

	                        _log_write(4, '(7X, A, I0, X, I0, A, I0, X, I0, A, I0)') "send from: ", comm%local_rank, comm%local_section,  " to  : ", comm%neighbor_rank, comm%neighbor_section, " send tag: ", send_tag
	                        _log_write(4, '(7X, A, I0, X, I0, A, I0, X, I0, A, I0)') "recv to  : ", comm%local_rank, comm%local_section, " from: ", comm%neighbor_rank, comm%neighbor_section, " recv tag: ", recv_tag

	                        assert_veq(comm%send_requests, MPI_REQUEST_NULL)
	                        assert_veq(comm%recv_requests, MPI_REQUEST_NULL)

	                        call mpi_irecv(comm%neighbor_rank,           1, MPI_INTEGER, old_comm%neighbor_rank, recv_tag, MPI_COMM_WORLD, comm%recv_requests(1), i_error); assert_eq(i_error, 0)
	                        call mpi_irecv(comm%neighbor_section,        1, MPI_INTEGER, old_comm%neighbor_rank, recv_tag, MPI_COMM_WORLD, comm%recv_requests(2), i_error); assert_eq(i_error, 0)
	                        call mpi_isend(comm%local_rank,              1, MPI_INTEGER, old_comm%neighbor_rank, send_tag, MPI_COMM_WORLD, comm%send_requests(1), i_error); assert_eq(i_error, 0)
	                        call mpi_isend(comm%local_section,           1, MPI_INTEGER, old_comm%neighbor_rank, send_tag, MPI_COMM_WORLD, comm%send_requests(2), i_error); assert_eq(i_error, 0)

	                        assert_vne(comm%send_requests, MPI_REQUEST_NULL)
	                        assert_vne(comm%recv_requests, MPI_REQUEST_NULL)
	                    end if
	                end do
	            end do
	        end do

	        !$omp barrier

	        !wait until all sections sent and received their communication changes
	        do i_section = i_first_local_section, i_last_local_section
	            section => src_grid%sections%elements_alloc(i_section)

	            _log_write(4, '(4X, A, I0)') "section wait for send/recv comm changes: ", i_section

	            do i_color = RED, GREEN
	                _log_write(4, '(5X, A, A)') trim(color_to_char(i_color)), ":"
	                do i_comm = 1, section%comms(i_color)%get_size()
	                    comm => section%comms(i_color)%elements(i_comm)

	                    _log_write(4, '(6X, A)') trim(comm%to_string())

	                    if (comm%neighbor_rank .ge. 0) then
	                        _log_write(4, '(7X, A, I0, X, I0, A, I0, X, I0)') "wait from: ", comm%local_rank, comm%local_section, " to  : ", comm%neighbor_rank, comm%neighbor_section
	                        _log_write(4, '(7X, A, I0, X, I0, A, I0, X, I0)') "wait to  : ", comm%local_rank, comm%local_section, " from: ", comm%neighbor_rank, comm%neighbor_section

	                        assert_vne(comm%send_requests, MPI_REQUEST_NULL)
	                        assert_vne(comm%recv_requests, MPI_REQUEST_NULL)

	                        call mpi_waitall(2, comm%send_requests, MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)
	                        call mpi_waitall(2, comm%recv_requests, MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)

	                        comm%send_requests = MPI_REQUEST_NULL
	                        comm%recv_requests = MPI_REQUEST_NULL
	                    end if
	                end do
	            end do
	        end do
        end subroutine

		subroutine send_recv_section_infos(src_grid, dest_grid, i_rank_out, i_rank_in)
		    type(t_grid), intent(inout)	                    :: src_grid, dest_grid
		    integer, intent(in)                             :: i_rank_out(:), i_rank_in(:)

		    integer						                    :: i_section, i_error, i_comm, i_color
            integer (kind = GRID_SI)						:: i_first_src_section, i_last_src_section, i_first_dest_section, i_last_dest_section
			type(t_section_info), allocatable               :: src_infos(:), dest_infos(:)
		    integer, allocatable 						    :: src_requests(:), dest_requests(:)
		    type(t_grid_section), pointer                   :: section

	        call src_grid%get_local_sections(i_first_src_section, i_last_src_section)
	        call dest_grid%get_local_sections(i_first_dest_section, i_last_dest_section)

			allocate(src_infos(i_first_src_section  : i_last_src_section), stat=i_error); assert_eq(i_error, 0)
			allocate(dest_infos(i_first_dest_section : i_last_dest_section), stat=i_error); assert_eq(i_error, 0)

			allocate(src_requests(i_first_src_section : i_last_src_section), stat=i_error); assert_eq(i_error, 0)
			allocate(dest_requests(i_first_dest_section : i_last_dest_section), stat=i_error); assert_eq(i_error, 0)

		    !exchange section infos
		    src_requests = MPI_REQUEST_NULL
		    dest_requests = MPI_REQUEST_NULL

		    do i_section = i_first_dest_section, i_last_dest_section
                if (i_rank_in(i_section) .ne. rank_MPI) then
                    _log_write(4, '("receiving: from: ", I0, " to: ", I0, " tag: ", I0)') i_rank_in(i_section), rank_MPI, i_section
                    call mpi_irecv(dest_infos(i_section), sizeof(dest_infos(i_section)), MPI_BYTE, i_rank_in(i_section), i_section, MPI_COMM_WORLD, dest_requests(i_section), i_error); assert_eq(i_error, 0)
                end if
 		    end do

		    do i_section = i_first_src_section, i_last_src_section
                if (i_rank_out(i_section) .ne. rank_MPI) then
                    section => src_grid%sections%elements_alloc(i_section)
                    src_infos(i_section) = section%get_info()

                    _log_write(4, '("sending:   from: ", I0, " to: ", I0, " tag: ", I0)') rank_MPI, i_rank_out(i_section), section%index
                    call mpi_isend(src_infos(i_section), sizeof(src_infos(i_section)), MPI_BYTE, i_rank_out(i_section), section%index, MPI_COMM_WORLD, src_requests(i_section), i_error); assert_eq(i_error, 0)
                end if
		    end do

		    do i_section = i_first_dest_section, i_last_dest_section
                if (i_rank_in(i_section) .ne. rank_MPI) then
                    call mpi_wait(dest_requests(i_section), MPI_STATUS_IGNORE, i_error); assert_eq(i_error, 0)

                    call dest_grid%sections%elements_alloc(i_section)%create(dest_infos(i_section))
                end if
 		    end do

		    call mpi_waitall(i_last_src_section - i_first_src_section + 1, src_requests, MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)

			deallocate(src_infos, stat=i_error); assert_eq(i_error, 0)
			deallocate(dest_infos, stat=i_error); assert_eq(i_error, 0)

			deallocate(src_requests, stat=i_error); assert_eq(i_error, 0)
			deallocate(dest_requests, stat=i_error); assert_eq(i_error, 0)
		end subroutine

		subroutine send_recv_section_data(src_grid, dest_grid, i_rank_out, i_rank_in)
		    type(t_grid), intent(inout)		                :: src_grid, dest_grid
		    integer, intent(in)                             :: i_rank_out(:), i_rank_in(:)

		    integer						                    :: i_section, i_comm, i_error, i_color
            integer (kind = GRID_SI)						:: i_first_src_section, i_last_src_section, i_first_dest_section, i_last_dest_section
		    type(t_grid_section), pointer					:: section
		    integer, allocatable 						    :: src_requests(:,:), dest_requests(:,:)
            integer (kind = GRID_DI)						:: tmp_distances(RED: GREEN)

	        call src_grid%get_local_sections(i_first_src_section, i_last_src_section)
	        call dest_grid%get_local_sections(i_first_dest_section, i_last_dest_section)

			allocate(src_requests(11, i_first_src_section : i_last_src_section), stat=i_error); assert_eq(i_error, 0)
			allocate(dest_requests(11, i_first_dest_section : i_last_dest_section), stat=i_error); assert_eq(i_error, 0)

		    !exchange sections
		    src_requests = MPI_REQUEST_NULL
		    dest_requests = MPI_REQUEST_NULL

		    do i_section = i_first_dest_section, i_last_dest_section
                select case (i_rank_in(i_section) - rank_MPI)
				case (0)
					!do nothing
				case default
                    section => dest_grid%sections%elements_alloc(i_section)

                    call mpi_irecv(section%t_global_data, sizeof(section%t_global_data), MPI_BYTE, i_rank_in(i_section), i_section, MPI_COMM_WORLD, dest_requests(1, i_section), i_error); assert_eq(i_error, 0)

                    call mpi_irecv(section%cells%get_c_pointer(),                  sizeof(section%cells%elements),            MPI_BYTE, i_rank_in(i_section), i_section, MPI_COMM_WORLD, dest_requests(2, i_section), i_error); assert_eq(i_error, 0)
                    call mpi_irecv(section%crossed_edges_in%get_c_pointer(),       sizeof(section%crossed_edges_in%elements), MPI_BYTE, i_rank_in(i_section), i_section, MPI_COMM_WORLD, dest_requests(3, i_section), i_error); assert_eq(i_error, 0)
                    call mpi_irecv(section%color_edges_in%get_c_pointer(),         sizeof(section%color_edges_in%elements),   MPI_BYTE, i_rank_in(i_section), i_section, MPI_COMM_WORLD, dest_requests(4, i_section), i_error); assert_eq(i_error, 0)
                    call mpi_irecv(section%nodes_in%get_c_pointer(),               sizeof(section%nodes_in%elements),         MPI_BYTE, i_rank_in(i_section), i_section, MPI_COMM_WORLD, dest_requests(5, i_section), i_error); assert_eq(i_error, 0)

                    do i_color = RED, GREEN
                        call mpi_irecv(section%boundary_edges(i_color)%get_c_pointer(),  sizeof(section%boundary_edges(i_color)%elements), MPI_BYTE, i_rank_in(i_section), i_section, MPI_COMM_WORLD, dest_requests(7 + i_color, i_section), i_error); assert_eq(i_error, 0)
                        call mpi_irecv(section%boundary_nodes(i_color)%get_c_pointer(),  sizeof(section%boundary_nodes(i_color)%elements), MPI_BYTE, i_rank_in(i_section), i_section, MPI_COMM_WORLD, dest_requests(9 + i_color, i_section), i_error); assert_eq(i_error, 0)
                        call mpi_irecv(section%comms(i_color)%get_c_pointer(),           sizeof(section%comms(i_color)%elements),          MPI_BYTE, i_rank_in(i_section), i_section, MPI_COMM_WORLD, dest_requests(11 + i_color, i_section), i_error); assert_eq(i_error, 0)
                    end do
                end select
		    end do

		    do i_section = i_first_src_section, i_last_src_section
                section => src_grid%sections%elements_alloc(i_section)

	            assert(section%is_forward() .eqv. src_grid%sections%is_forward())

                select case (i_rank_out(i_section) - rank_MPI)
				case (0)
				    !local sections are simply copied to the new list
                    dest_grid%sections%elements_alloc(section%index) = section
				case default
                    call mpi_isend(section%t_global_data, sizeof(section%t_global_data), MPI_BYTE, i_rank_out(i_section), section%index, MPI_COMM_WORLD, src_requests(1, i_section), i_error); assert_eq(i_error, 0)

                    call mpi_isend(section%cells%get_c_pointer(),             sizeof(section%cells%elements),                MPI_BYTE, i_rank_out(i_section), section%index, MPI_COMM_WORLD, src_requests(2, i_section), i_error); assert_eq(i_error, 0)
                    call mpi_isend(section%crossed_edges_in%get_c_pointer(),  sizeof(section%crossed_edges_in%elements),     MPI_BYTE, i_rank_out(i_section), section%index, MPI_COMM_WORLD, src_requests(3, i_section), i_error); assert_eq(i_error, 0)
                    call mpi_isend(section%color_edges_in%get_c_pointer(),    sizeof(section%color_edges_in%elements),       MPI_BYTE, i_rank_out(i_section), section%index, MPI_COMM_WORLD, src_requests(4, i_section), i_error); assert_eq(i_error, 0)
                    call mpi_isend(section%nodes_in%get_c_pointer(),          sizeof(section%nodes_in%elements),             MPI_BYTE, i_rank_out(i_section), section%index, MPI_COMM_WORLD, src_requests(5, i_section), i_error); assert_eq(i_error, 0)

                    do i_color = RED, GREEN
                        call mpi_isend(section%boundary_edges(i_color)%get_c_pointer(), sizeof(section%boundary_edges(i_color)%elements),  MPI_BYTE, i_rank_out(i_section), section%index, MPI_COMM_WORLD, src_requests(7 + i_color, i_section), i_error); assert_eq(i_error, 0)
                        call mpi_isend(section%boundary_nodes(i_color)%get_c_pointer(), sizeof(section%boundary_nodes(i_color)%elements),  MPI_BYTE, i_rank_out(i_section), section%index, MPI_COMM_WORLD, src_requests(9 + i_color, i_section), i_error); assert_eq(i_error, 0)
                        call mpi_isend(section%comms(i_color)%get_c_pointer(), sizeof(section%comms(i_color)%elements),                    MPI_BYTE, i_rank_out(i_section), section%index, MPI_COMM_WORLD, src_requests(11 + i_color, i_section), i_error); assert_eq(i_error, 0)
                    end do
		        end select
		    end do

		    !wait until the load has been distributed by all neighbor threads and processes
			!$omp barrier

		    do i_section = i_first_src_section, i_last_src_section
				section => src_grid%sections%elements_alloc(i_section)

                select case (i_rank_out(i_section) - rank_MPI)
                    case (0)
						!do nothing
                    case default
                        call mpi_waitall(11, src_requests(:, i_section), MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)

                        call section%destroy()
                end select
 		    end do

			do i_section = i_first_dest_section, i_last_dest_section
				section => dest_grid%sections%elements_alloc(i_section)

                select case (i_rank_in(i_section) - rank_MPI)
                    case (0)
						!do nothing
					case default
                        call mpi_waitall(11, dest_requests(:, i_section), MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)

                        !$omp atomic
                        src_grid%max_dest_stack(RED) = src_grid%max_dest_stack(RED) + section%max_dest_stack(RED) - section%min_dest_stack(RED) + 1
                        !$omp atomic
                        src_grid%max_dest_stack(GREEN) = src_grid%max_dest_stack(GREEN) + section%max_dest_stack(GREEN) - section%min_dest_stack(GREEN) + 1
                end select

	            assert(section%is_forward() .eqv. src_grid%sections%is_forward())
	            assert_eq(section%cells%elements(1)%get_previous_edge_type(), OLD_BND)
	        end do

			deallocate(src_requests, stat=i_error); assert_eq(i_error, 0)
			deallocate(dest_requests, stat=i_error); assert_eq(i_error, 0)
		end subroutine
#	endif

    subroutine find_section_boundary_elements(thread, section, last_cell_index, last_crossed_edge_data)
        type(t_grid_thread), intent(inout)	                        :: thread
        type(t_grid_section), intent(inout)	                        :: section
        integer (kind = GRID_DI), intent(in)		                :: last_cell_index
        type(t_edge_data), intent(in)	                            :: last_crossed_edge_data

        type(t_edge_data), pointer                                  :: p_edge, p_edges(:)
        type(t_node_data), pointer                                  :: p_nodes(:)
        integer(kind = GRID_SI)							            :: i_color, i_pass, i_edges, i_nodes, i

        _log_write(4, '(2X, A, I0)') "find section boundary elements: ", section%index

        !The next edge of the last cell is a new boundary edge.
        !As the cell is already reversed, we therefore have to set the PREVIOUS EDGE to an OLD BOUNDARY edge instead.

        assert_ge(last_cell_index, 1)
        assert_le(last_cell_index, section%cells%get_size())
        call section%cells%elements(last_cell_index)%set_previous_edge_type(int(OLD_BND, 1))

        !and move the last crossed edge to the red stack
        p_edge => thread%edges_stack(RED)%push()
        p_edge = last_crossed_edge_data

        !find additional process edges and nodes
        do i_color = RED, GREEN
            !All cell indices that remain on the index stack are of cells with new boundary color edges.
            !As the cells are already reversed, we therefore have to set the COLOR EDGE to an OLD BOUNDARY EDGE for each cell instead.

            do i = 1, thread%indices_stack(i_color)%i_current_element
                call section%cells%elements(thread%indices_stack(i_color)%elements(i))%set_color_edge_type(int(OLD_BND, 1))
            end do

            section%min_distance(i_color) = 0

            !all edges on the boundary stream are old boundary edges

            i_edges = section%boundary_edges(i_color)%i_current_element
            p_edges => section%boundary_edges(i_color)%elements(i_edges : 1 : -1)
            call section%boundary_type_edges(OLD, i_color)%attach(p_edges, .not. section%boundary_edges(i_color)%is_forward())

            !reverse the stream so it is aligned with the traversal direction
            call section%boundary_type_edges(OLD, i_color)%reverse()

            if (i_edges > 0) then
                p_edges(1)%min_distance = 0
                p_edges(2 : i_edges)%min_distance = encode_edge_size(p_edges(1 : i_edges - 1)%depth)
                call prefix_sum(p_edges%min_distance, p_edges%min_distance)

                section%start_distance(i_color) = p_edges(i_edges)%min_distance + encode_edge_size(p_edges(i_edges)%depth)
            else
                section%start_distance(i_color) = 0
            end if

            !all nodes on the boundary stream are old boundary nodes

            i_nodes = section%boundary_nodes(i_color)%i_current_element
            p_nodes => section%boundary_nodes(i_color)%elements(i_nodes : 1 : -1)
            call section%boundary_type_nodes(OLD, i_color)%attach(p_nodes, .not. section%boundary_nodes(i_color)%is_forward())

            !reverse the stream so it is aligned with the traversal direction
            call section%boundary_type_nodes(OLD, i_color)%reverse()

            if (i_color == GREEN) then
                assert_eq(i_nodes, i_edges + 1)
                p_nodes(1 : i_edges)%distance = p_edges%min_distance
                p_nodes(i_nodes)%distance = section%start_distance(i_color)
            else
                assert_eq(i_nodes, i_edges)
                p_nodes%distance = p_edges%min_distance
            end if

            !all remaining edges on the stack are new boundary edges

            i_edges = thread%edges_stack(i_color)%i_current_element
            p_edges => section%boundary_edges(i_color)%elements(section%boundary_edges(i_color)%i_current_element + 1 : section%boundary_edges(i_color)%i_current_element + i_edges)
            section%boundary_edges(i_color)%i_current_element = section%boundary_edges(i_color)%i_current_element + i_edges
            p_edges = thread%edges_stack(i_color)%elements(1 : i_edges)
            call section%boundary_type_edges(NEW, i_color)%attach(p_edges, section%boundary_edges(i_color)%is_forward())

            if (i_edges > 0) then
                p_edges(1)%min_distance = 0
                p_edges(2 : i_edges)%min_distance = encode_edge_size(p_edges(1 : i_edges - 1)%depth)
                call prefix_sum(p_edges%min_distance, p_edges%min_distance)

                section%end_distance(i_color) = p_edges(i_edges)%min_distance + encode_edge_size(p_edges(i_edges)%depth)
            else
                section%end_distance(i_color) = 0
            end if

            !all remaining nodes on the stack are new boundary nodes

            i_nodes = thread%nodes_stack(i_color)%i_current_element
            p_nodes => section%boundary_nodes(i_color)%elements(section%boundary_nodes(i_color)%i_current_element : section%boundary_nodes(i_color)%i_current_element + i_nodes)
            section%boundary_nodes(i_color)%i_current_element = section%boundary_nodes(i_color)%i_current_element + i_nodes
            p_nodes(2:) = thread%nodes_stack(i_color)%elements(1 : i_nodes)
            call section%boundary_type_nodes(NEW, i_color)%attach(p_nodes, section%boundary_nodes(i_color)%is_forward())

            if (i_color == GREEN) then
                assert_eq(i_nodes , i_edges)
                p_nodes(1 : i_nodes)%distance = p_edges%min_distance
                p_nodes(i_nodes + 1)%distance = section%end_distance(i_color)
            else
                assert_eq(i_nodes + 1, i_edges)
                p_nodes%distance = p_edges%min_distance
            end if

            thread%indices_stack(i_color)%i_current_element = 0
            thread%edges_stack(i_color)%i_current_element = 0
            thread%nodes_stack(i_color)%i_current_element = 0
        end do

        call section%cells%trim()
        call section%crossed_edges_out%trim()
        section%crossed_edges_in = section%crossed_edges_out

        call section%nodes_out%trim()
        section%nodes_in = section%nodes_out

        call section%color_edges_out%trim()
        section%color_edges_in = section%color_edges_out

        call section%boundary_nodes(RED)%trim()
        call section%boundary_nodes(GREEN)%trim()
        call section%boundary_edges(RED)%trim()
        call section%boundary_edges(GREEN)%trim()

#	    if (_DEBUG_LEVEL > 3)
            do i_color = RED, GREEN
                _log_write(4, '(3X, A, A)') trim(color_to_char(i_color)), ":"
                do i_pass = OLD, NEW
                    _log_write(4, '(4X, A, A)') trim(edge_type_to_char(i_pass)), ":"

                    p_edges => section%boundary_type_edges(i_pass, i_color)%elements
                    _log_write(4, '(5X, A, I0)') "boundary edges: ", size(p_edges)
                    do i = 1, size(p_edges)
                        _log_write(4, '(6X, F0.4, X, F0.4)') decode_distance(p_edges(i)%min_distance), decode_distance(p_edges(i)%min_distance + encode_edge_size(p_edges(i)%depth))
                    end do

                    p_nodes => section%boundary_type_nodes(i_pass, i_color)%elements
                    _log_write(4, '(5X, A, I0)') "boundary nodes: ", size(p_nodes)
                    do i = 1, size(p_nodes)
                        _log_write(4, '(6X, F0.4)') decode_distance(p_nodes(i)%distance)
                    end do
                end do
            end do
#	    endif

#	    if (_DEBUG_LEVEL > 4)
            _log_write(5, '(2X, A)') "destination section final state :"
            call section%print()
#	    endif
    end subroutine
end module
