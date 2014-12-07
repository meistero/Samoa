! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE_DAMBREAK_CLASSIC)
	MODULE SWE_compute_error
		use LIB_VTK_IO

		use SFC_edge_traversal

		use Samoa_swe
		use SWE_euler_timestep

		implicit none

		!> Norm data
		type t_norm_data
            real (kind = GRID_SR)                   :: error_l1, error_l2, error_max
		end type

		type t_output_point_data
			type(t_state)											:: Q
			real (kind = GRID_SR), dimension(2)						:: coords		!< position
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
            type(t_norm_data)                  :: norm_data_h, norm_data_u

            type(t_output_point_data), allocatable		            :: point_data(:)
            type(t_output_cell_data), allocatable			        :: cell_data(:)
            character(len=64)							            :: s_file_stamp

            integer (kind = GRID_SI)								:: i_output_iteration = 0
            integer (kind = GRID_SI)								:: i_point_data_index
            integer (kind = GRID_SI)								:: i_cell_data_index
        end type

		integer, parameter		:: i_element_order = 0
		!real (kind = GRID_SR), allocatable		:: r_testpoints(:,:)


        integer, parameter      :: out_unit = 20
        character (len = 128)   :: pout_file_name

		type(t_gv_Q)			:: gv_Q

#		define _GT_NAME					t_swe_compute_error_traversal

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
		type(t_swe_compute_error_traversal), intent(inout)		:: traversal
		type(t_grid), intent(inout)				:: grid

        call scatter(traversal%s_file_stamp, traversal%children%s_file_stamp)
        call scatter(traversal%i_output_iteration, traversal%children%i_output_iteration)
        call scatter(grid%r_time, grid%sections%elements_alloc(:)%r_time)
	end subroutine

        subroutine post_traversal_grid_op(traversal, grid)
		type(t_swe_compute_error_traversal), intent(inout)		:: traversal
		type(t_grid), intent(inout)				:: grid


        integer                         :: i_error, i, j, k, alloc_err
		integer(4)			:: i_rank, i_section, e_io
		logical                         :: l_exists

    	real (kind = GRID_SR), pointer				:: dummy_points(:,:,:) => null() !groesse/dimension auf null setzen
    	real (kind = GRID_SR), pointer 				:: big_points_array(:,:,:) => null()

        integer                                 :: counter = 1

#    if defined(_MPI)
            call mpi_barrier(MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)

            call reduce(traversal%norm_data_h%error_l1, traversal%children%norm_data_h%error_l1, MPI_SUM, .true.)
            call reduce(traversal%norm_data_h%error_l2, traversal%children%norm_data_h%error_l2, MPI_SUM, .true.)
            call reduce(traversal%norm_data_h%error_max, traversal%children%norm_data_h%error_max, MPI_MAX, .true.)
            call reduce(traversal%norm_data_u%error_l1, traversal%children%norm_data_u%error_l1, MPI_SUM, .true.)
            call reduce(traversal%norm_data_u%error_l2, traversal%children%norm_data_u%error_l2, MPI_SUM, .true.)
            call reduce(traversal%norm_data_u%error_max, traversal%children%norm_data_u%error_max, MPI_MAX, .true.)

#     endif

        traversal%norm_data_h%error_l2 = sqrt(traversal%norm_data_h%error_l2)
        traversal%norm_data_u%error_l2 = sqrt(traversal%norm_data_u%error_l2)

        traversal%i_output_iteration = traversal%i_output_iteration + 1

        if (rank_MPI == 0) then
!           _log_write(1, '(A, ES14.7)') " SWE 1D BENCHMARK: water height L1-Error: ", traversal%norm_data_h%error_l1
!           _log_write(1, '(A, ES14.7)') " SWE 1D BENCHMARK: water height L2-Error: ", traversal%norm_data_h%error_l2
!           _log_write(1, '(A, ES14.7)') " SWE 1D BENCHMARK: water height Max-Error: ", traversal%norm_data_h%error_max
!           _log_write(1, '(A, ES14.7)') " SWE 1D BENCHMARK: velocity L1-Error: ", traversal%norm_data_u%error_l1
!           _log_write(1, '(A, ES14.7)') " SWE 1D BENCHMARK: velocity L2-Error: ", traversal%norm_data_u%error_l2
!           _log_write(1, '(A, ES14.7)') " SWE 1D BENCHMARK: velocity Max-Error: ", traversal%norm_data_u%error_max

            ! file schreiben
            write(pout_file_name, "(A, A, A, I2, A, I2, A, F6.5, A, F6.5, A, I3, A)") "erroroutput", TRIM(traversal%s_file_stamp), "_dmin", cfg%i_min_depth, "_dmax", cfg%i_max_depth, "_cou", cfg%courant_number, "_dry", cfg%dry_tolerance, "_p", size_MPI, ".txt"

            open(unit=out_unit, file=pout_file_name, action="write", status="replace")
                write(out_unit, "(A)") "dmin, dmax, cou, dry_tol, processes, sim_time, h_error_l1, h_error_l2, h_error_max, u_error_l1, u_error_l2, u_error_max"
                write(out_unit, "(2(I0, A), 2(F6.5, A), I0, A, F12.4, 6(A, ES14.7))") cfg%i_min_depth, ", ", cfg%i_max_depth, ", ", cfg%courant_number, ", ", cfg%dry_tolerance, ", ", size_MPI, ", ", cfg%t_phase, ", ", traversal%norm_data_h%error_l1, ", ", traversal%norm_data_h%error_l2, ", ", traversal%norm_data_h%error_max, ", ", traversal%norm_data_u%error_l1, ", ", traversal%norm_data_u%error_l2, ", ", traversal%norm_data_u%error_max
            close(out_unit)

        end if

    end subroutine

	subroutine pre_traversal_op(traversal, section)
		type(t_swe_compute_error_traversal), intent(inout)		:: traversal
		type(t_grid_section), intent(inout)			:: section

		type(t_section_info)                         :: info
!		integer (kind = GRID_SI)		     :: i_error, i_cells, i_points
    
        !info = section%get_info()

        !traversal%cells = info%i_cells

        traversal%norm_data_h%error_l1 = 0.0_GRID_SR
        traversal%norm_data_h%error_l2 = 0.0_GRID_SR
        traversal%norm_data_h%error_max = 0.0_GRID_SR
        traversal%norm_data_u%error_l1 = 0.0_GRID_SR
        traversal%norm_data_u%error_l2 = 0.0_GRID_SR
        traversal%norm_data_u%error_max = 0.0_GRID_SR

	end subroutine

	subroutine post_traversal_op(traversal, section)
		type(t_swe_compute_error_traversal), intent(inout)				:: traversal
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
		type(t_swe_compute_error_traversal), intent(inout)				:: traversal
		type(t_grid_section), intent(inout)							:: section
		type(t_element_base), intent(inout)					:: element

		!local variables
!		integer (kind = GRID_SI)							:: i,j
		type(t_state), dimension(_SWE_CELL_SIZE)			:: Q
		type(t_state), dimension(6)							:: Q_test
!		real (kind = GRID_SR) 						        :: h, b, p(2), local_coord(2), epsvec(2), eps, distvec(2), dist

        real (kind = GRID_SR)    ::  g, hl, hr, h, u, h_diff, u_diff, time, cell_volume
        real (kind = GRID_SR)    ::  state_left, state_rarefaction, state_intermediate   
        real (kind = GRID_SR)    ::  domain_length   !length of domain
        real (kind = GRID_SR)    ::  domain_first    !first coordinate of domain
        real (kind = GRID_SR)    ::  shock_position  !shock_position in domain coordinates
        real (kind = GRID_SR)    ::  coords(2)       !cell center coordinates


		call gv_Q%read(element, Q)

        time = section%r_time

        !compute analytic solution 

        g = 9.80665_GRID_SR
        hl = 3.0_GRID_SR
        hr = 1.0_GRID_SR

        state_left = -5.42401603979929_GRID_SR
        state_rarefaction = -1.92518574923489_GRID_SR
        state_intermediate = 5.08133721790905_GRID_SR

        domain_length = cfg%scaling
        domain_first = cfg%offset(1)
        shock_position = 0.0_GRID_SR

        coords = samoa_barycentric_to_world_point(element%transform_data, [1.0_GRID_SR/3.0_GRID_SR, 1.0_GRID_SR/3.0_GRID_SR]) * cfg%scaling + cfg%offset

        !compute h and u:coords(1)
        if (coords(1) < state_left * time) then
            h = 3.0
            u = 0.0
        else if (coords(1) < state_rarefaction * time) then
            h = 0.0113301801441992_GRID_SR * (10.8480320795986_GRID_SR - coords(1) / time)**2
            u = -(2.0/3.0) * sqrt((10.8480320795986_GRID_SR - coords(1) / time)**2) + 10.8480320795986_GRID_SR
        else if (coords(1) < state_intermediate * time) then
            h = 1.84857660309676
            u = 2.33255352704294
        else
            h = 1.0
            u = 0.0
        end if 

        !compute l1 and l2 local l1 and l2 norm (adding everything up is done by mpi, sqrt is done afterwards)
        h_diff = abs(h - t_basis_Q_eval([1.0_GRID_SR/3.0_GRID_SR, 1.0_GRID_SR/3.0_GRID_SR], Q%h))
        u_diff = abs(u - t_basis_Q_eval([1.0_GRID_SR/3.0_GRID_SR, 1.0_GRID_SR/3.0_GRID_SR], Q%p(1)) / t_basis_Q_eval([1.0_GRID_SR/3.0_GRID_SR, 1.0_GRID_SR/3.0_GRID_SR], Q%h))

        cell_volume = element%cell%geometry%get_volume()

        traversal%norm_data_h%error_l1 = traversal%norm_data_h%error_l1 + h_diff * cell_volume
        traversal%norm_data_h%error_l2 = traversal%norm_data_h%error_l2 + h_diff * h_diff * cell_volume
        traversal%norm_data_u%error_l1 = traversal%norm_data_u%error_l1 + u_diff * cell_volume
        traversal%norm_data_u%error_l2 = traversal%norm_data_u%error_l2 + u_diff * u_diff * cell_volume

        traversal%norm_data_h%error_max = max(h_diff, traversal%norm_data_h%error_max)
        traversal%norm_data_u%error_max = max(u_diff, traversal%norm_data_u%error_max)

	end subroutine
	END MODULE
#endif
