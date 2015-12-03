! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
	MODULE SWE_point_output
		use LIB_VTK_IO

		use SFC_edge_traversal

		use Samoa_swe
		use SWE_euler_timestep

		implicit none

		!> Output point data
		type t_output_point_data
			type(t_state)					:: Q
			real (kind = GRID_SR), dimension(2)		:: coords		!< position
		end type t_output_point_data

		!> Output cell dat
		type t_output_cell_data
			type(t_state)					:: Q
			integer (kind = GRID_SI)			:: rank
			integer (kind = GRID_SI)			:: section_index
			integer (kind = BYTE)				:: depth
			integer (kind = BYTE)				:: refinement
		end type t_output_cell_data

	        type num_traversal_data
	            type(t_output_point_data), allocatable	        :: point_data(:)
	            type(t_output_cell_data), allocatable		:: cell_data(:)
	            character(len=64)					:: s_file_stamp

	            integer (kind = GRID_SI)				:: i_output_iteration=0
	            integer (kind = GRID_SI)				:: i_point_data_index
	            integer (kind = GRID_SI)				:: i_cell_data_index
	        end type

		integer, parameter		:: i_element_order = 0
		real (kind = GRID_SR), allocatable		:: r_testpoints(:,:)

        integer, parameter      :: out_unit = 20
        character (len = 256)   :: pout_file_name

		type(t_gv_Q)			:: gv_Q

#		define _GT_NAME					t_swe_point_output_traversal

#		define _GT_EDGES
#		define _GT_EDGES_TEMP

#		define _GT_PRE_TRAVERSAL_OP			pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP			post_traversal_op
#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op

#		define _GT_ELEMENT_OP				element_op

#		define _GT_CELL_TO_EDGE_OP			cell_to_edge_op

#		include "SFC_generic_traversal_ringbuffer.f90"

  subroutine pre_traversal_grid_op(traversal, grid)
		type(t_swe_point_output_traversal), intent(inout)		:: traversal
		type(t_grid), intent(inout)				:: grid

		integer 						:: i, erro

        if (rank_MPI == 0) then
            _log_write(1, '(A, I0)') " SWE: output step: ", traversal%i_output_iteration
        end if

        ! initialize test point array

		if (allocated(r_testpoints)) then
        	deallocate(r_testpoints)
    	end if
		allocate (r_testpoints(size(cfg%r_testpoints, dim=1), 8), stat = erro)

        r_testpoints(:,:) = 0

        ! load with 0-1 coordinates and check for correctness

		do i=1, size(r_testpoints, dim=1)
			r_testpoints(i,1) = (cfg%r_testpoints(i,1) - cfg%offset(1)) / cfg%scaling
			try((0.0 <= r_testpoints(i,1) .and. r_testpoints(i,1) <= 1.0), 'Submitted testpoints contain invalid x-coordinate')
			r_testpoints(i,2) = (cfg%r_testpoints(i,2) - cfg%offset(2)) / cfg%scaling
			try((0.0 <= r_testpoints(i,2) .and. r_testpoints(i,2) <= 1.0), 'Submitted testpoints contain invalid y-coordinate')
            !r_testpoints(i,8) = 0
            r_testpoints(i,7) = huge(1.0)
		end do

        call scatter(traversal%s_file_stamp, traversal%children%s_file_stamp)
        call scatter(traversal%i_output_iteration, traversal%children%i_output_iteration)
        call scatter(grid%r_time, grid%sections%elements_alloc(:)%r_time)
	end subroutine

        subroutine post_traversal_grid_op(traversal, grid)
		type(t_swe_point_output_traversal), intent(inout)		:: traversal
		type(t_grid), intent(inout)				:: grid


        integer                         :: i_error, i, j, k, alloc_err
		integer(4)			:: i_rank, i_section, e_io
		logical                         :: l_exists

    	real (kind = GRID_SR), pointer				:: dummy_points(:,:,:) => null() !groesse/dimension auf null setzen
    	real (kind = GRID_SR), pointer 				:: big_points_array(:,:,:) => null()

        integer                                 :: counter = 1

#    if defined(_MPI)
            call mpi_barrier(MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)

            !array of point arrays for gather
            if (rank_MPI == 0) then
                if (associated(big_points_array)) then
        	        deallocate(big_points_array)
    	        end if

                allocate (big_points_array(size(r_testpoints,dim=1),size(r_testpoints,dim=2),size_MPI), stat = alloc_err)
                if (alloc_err > 0) then
                    write(*,'(A)') "Error when trying to allocate rank-0-pointmatrix"
                end if
            else
                if (associated(dummy_points)) then
        	        deallocate(dummy_points)
              	end if
	            allocate (dummy_points(1,1,1), stat = alloc_err)
                if (alloc_err > 0) then
                    write(*,'(A)') "Error when trying to allocate rank-0-pointmatrix"
                end if
	            big_points_array => dummy_points
            end if

            call MPI_gather(r_testpoints(1,1), sizeof(r_testpoints), MPI_BYTE, &
              big_points_array(1,1,1), sizeof(r_testpoints), MPI_BYTE, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)

            if (rank_MPI == 0) then

                do j=1, size(r_testpoints, dim=1) !for every point
                     do i=2, (size_MPI) !for every rank
                         if (big_points_array(j,7,i) < r_testpoints(j,7)) then
                            r_testpoints(j,:) = big_points_array(j,:,i)
                        end if
                    end do
                end do

            end if

            ! sum up p, h, b, how often processed, time
            !if (rank_MPI == 0) then
            !    do j=1, size(r_testpoints, dim=1) !for every point
            !        do k=3, 8 !for every quantitiy where it makes sense: 3,4,5,6,7, 8=time simpler if mean computed
            !            do i=2, (size_MPI) !for every rank
            !                r_testpoints(j,k) = r_testpoints(j,k) + big_points_array(j,k,i)
            !            end do
            !        end do
            !    end do

                !compute means
                !do j=1, size(r_testpoints, dim=1)
                !    do i=3,6
                !        try((r_testpoints(j,7) > dble(0.0)), 'testpoints not calculated')
                !        r_testpoints(j,i) = r_testpoints(j,i) / r_testpoints(j,7)
                !    end do
                !    r_testpoints(j,8) = r_testpoints(j,8) / r_testpoints(j,7)
                !end do
            !end if

#     endif

        if (rank_MPI == 0) then
            ! file schreiben
            write(pout_file_name, "(A, A, I0, A, F10.3, A, F10.3, A, I0, A)") TRIM(traversal%s_file_stamp), "_d", cfg%i_max_depth, "_cou", cfg%courant_number, "_dry", cfg%dry_tolerance, "_pointoutput_", traversal%i_output_iteration, ".txt"

            open(unit=out_unit, file=pout_file_name, action="write", status="replace")
                write(out_unit, "(A)") "x, y, z, p(1), p(2), h, b, dist_to_cell_center, time"
                do i=1, size(r_testpoints, dim=1)
                    write(out_unit, "(8(F0.15,A),F0.15)") cfg%r_testpoints(i,1), ",", cfg%r_testpoints(i,2), ",", dble(0.0), ",", r_testpoints(i,3), ",", r_testpoints(i,4), ",", r_testpoints(i,5), ",", r_testpoints(i,6), ",", r_testpoints(i,7), ",", r_testpoints(i,8)
                end do
            close(out_unit)
        end if

        traversal%i_output_iteration = traversal%i_output_iteration + 1

    end subroutine

	subroutine pre_traversal_op(traversal, section)
		type(t_swe_point_output_traversal), intent(inout)		:: traversal
		type(t_grid_section), intent(inout)			:: section

		type(t_section_info)                         :: grid_info
		integer (kind = GRID_SI)		     :: i_error, i_cells, i_points

	end subroutine

	subroutine post_traversal_op(traversal, section)
		type(t_swe_point_output_traversal), intent(inout)				:: traversal
		type(t_grid_section), intent(inout)					:: section

		integer (kind = GRID_SI), dimension(:), allocatable			:: i_offsets
		integer (1), dimension(:), allocatable					:: i_types
		integer (kind = GRID_SI), dimension(:), allocatable			:: i_connectivity
		real (kind = GRID_SR), dimension(:, :), allocatable			:: r_velocity
		real (kind = GRID_SR), dimension(:), allocatable			:: r_empty
        	type(t_vtk_writer)                                       		:: vtk

		type(t_section_info)                                        :: grid_info
		integer (kind = GRID_SI)				    :: i_error, i_cells, i_points, i
		integer(4)						    :: e_io

		traversal%i_output_iteration = traversal%i_output_iteration + 1
	end subroutine


		!******************
		!Geometry operators
		!******************

	subroutine element_op(traversal, section, element)
		type(t_swe_point_output_traversal), intent(inout)				:: traversal
		type(t_grid_section), intent(inout)							:: section
		type(t_element_base), intent(inout)					:: element

		!local variables
		integer (kind = GRID_SI)							:: i,j
		type(t_state), dimension(_SWE_CELL_SIZE)			:: Q
		type(t_state), dimension(6)							:: Q_test
		real (kind = GRID_SR) 						        :: h, b, p(2), local_coord(2), epsvec(2), eps, distvec(2), dist

		call gv_Q%read(element, Q)

        epsvec = samoa_world_to_barycentric_vector(element%transform_data, [0.001_GRID_SR, 0.0_GRID_SR])
        eps = sqrt(dot_product(epsvec, epsvec))

		do i=1, size(r_testpoints, dim=1)
			!check ob koordinaten in aktueller zelle
			local_coord = samoa_world_to_barycentric_point(element%transform_data, [r_testpoints(i,1), r_testpoints(i,2)])
			if (local_coord(1) + eps >= 0.0_GRID_SR .and. local_coord(2) + eps >= 0.0_GRID_SR .and. local_coord(1) + local_coord(2) - eps <= 1.0_GRID_SR) then

                distvec = samoa_barycentric_to_world_point(element%transform_data, [1.0_GRID_SR/3.0_GRID_SR, 1.0_GRID_SR/3.0_GRID_SR]) - [r_testpoints(i,1), r_testpoints(i,2)]
                dist = sqrt(dot_product(distvec, distvec))

                if (dist < r_testpoints(i,7)) then
                    h = t_basis_Q_eval(local_coord, Q%h)			!insert transformed coordinates
                    b = t_basis_Q_eval(local_coord, Q%b)
                    p(1) = t_basis_Q_eval(local_coord, Q%p(1))	    !insert transformed coordinates, gives world coordinates
                    p(2) = t_basis_Q_eval(local_coord, Q%p(2))

                    r_testpoints(i,3) = p(1)
                    r_testpoints(i,4) = p(2)
                    r_testpoints(i,5) = h
                    r_testpoints(i,6) = b
                    r_testpoints(i,7) = dist
                    r_testpoints(i,8) = section%r_time / (24.0_SR * 60.0_SR * 60.0_SR) + 70.2422_SR !convert to days since start of year for comparison with buoy data
                end if
			end if
		end do

	end subroutine
	END MODULE
#endif
