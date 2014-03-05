#include "Compilation_control.f90"

#if defined(_FLASH)
	MODULE FLASH
		use SFC_edge_traversal
		use FLASH_data_types
		use FLASH_dg_element

		use FLASH_adapt
		use FLASH_initialize
		use FLASH_output
		use FLASH_xml_output
		use FLASH_euler_timestep

		use Samoa_FLASH

		implicit none

		PRIVATE
		PUBLIC t_FLASH

		type t_FLASH
            type(t_FLASH_init_traversal)              :: init
            type(t_FLASH_output_traversal)            :: output
            type(t_FLASH_xml_output_traversal)        :: xml_output
            type(t_FLASH_euler_timestep_traversal)    :: euler
            type(t_FLASH_adaption_traversal)          :: adaption

            contains

            procedure , pass :: create => FLASH_create
            procedure , pass :: run => FLASH_run
            procedure , pass :: destroy => FLASH_destroy
        end type

		contains

		!> Creates all required runtime objects for the scenario
		subroutine FLASH_create(FLASH, grid, l_log, i_asagi_mode)
            class(t_FLASH), intent(inout)                                 :: FLASH
			type(t_grid), intent(inout)									:: grid
			logical, intent(in)						                    :: l_log
			integer, intent(in)											:: i_asagi_mode

			!local variables
			character(64)												:: s_log_name, s_date, s_time


                character(*), parameter :: s_solver = "flash solver"

            !$omp master
            if (rank_MPI == 0) then
                _log_write(0, '(" FLASH: flux solver: ", A)') s_solver
            end if
            !$omp end master

			!open log file
			call date_and_time(s_date, s_time)
			write (FLASH%output%s_file_stamp, "(A, A, A8, A, A6)") "output/FLASH", "_", s_date, "_", s_time
			write (FLASH%xml_output%s_file_stamp, "(A, A, A8, A, A6)") "output/FLASH", "_", s_date, "_", s_time
			write (s_log_name, '(A, A)') TRIM(FLASH%xml_output%s_file_stamp), ".log"

			if (l_log) then
				_log_open_file(s_log_name)
			endif

			!create a quadrature rule
			call t_qr_create_dunavant_rule(qr_Q, max(1, 2 * _FLASH_ORDER))

			call load_scenario(grid, i_asagi_mode, "data/alaska/displ.nc", "data/alaska/bath.nc", 2.0d6, [-0.5d6, -1.5d6])
		end subroutine

		subroutine load_scenario(grid, i_asagi_mode, ncd_displ, ncd_bath, scaling, offset)
			type(t_grid), target, intent(inout)     :: grid
			integer, intent(in)						:: i_asagi_mode

    		        character(*), intent(in)                :: ncd_displ, ncd_bath
            		double precision, intent(in)            :: scaling, offset(2)
			integer                                 :: i_error, k
			integer, pointer						:: afh

#			if defined(_ASAGI)
				grid%afh_displacement = asagi_create(grid_type = GRID_FLOAT, hint = i_asagi_mode, levels = 1)
				grid%afh_bathymetry = asagi_create(grid_type = GRID_FLOAT, hint = i_asagi_mode, levels = 1)

                i_error = asagi_open(grid%afh_displacement, "data/displ.nc", 0); assert_eq(i_error, GRID_SUCCESS)
                i_error = asagi_open(grid%afh_bathymetry, "data/bath.nc", 0); assert_eq(i_error, GRID_SUCCESS)

                if (rank_MPI == 0) then
                    afh => grid%afh_displacement
                    _log_write(1, '(A, A, A, F0.2, A, F0.2, A, F0.2, A, F0.2, A)') " FLASH: loaded '", "data/displ.nc", "', coordinate system: [", grid_min_x(afh), ", ", grid_min_y(afh), "] x [", grid_max_x(afh), ", ", grid_max_y(afh), "]"

                    afh => grid%afh_bathymetry
                    _log_write(1, '(A, A, A, F0.2, A, F0.2, A, F0.2, A, F0.2, A)') " FLASH: loaded '", "data/bath.nc", "', coordinate system: [", grid_min_x(afh), ", ", grid_min_y(afh), "] x [", grid_max_x(afh), ", ", grid_max_y(afh), "]"
                end if
         
 	        grid%scaling = scaling
                grid%offset = offset
#			endif
		end subroutine

		!> Destroys all required runtime objects for the scenario
		subroutine FLASH_destroy(FLASH, grid, l_log)
            class(t_FLASH), intent(inout)                                  :: FLASH
			type(t_grid), intent(inout)									:: grid
			logical, intent(in)						                    :: l_log

#			if defined(_ASAGI)
				call asagi_close(grid%afh_displacement)
				call asagi_close(grid%afh_bathymetry)
#			endif

			if (l_log) then
				_log_close_file()
			endif

			call t_qr_destroy(qr_Q)
		end subroutine

		!*********************************
		! run()-method
		!*********************************

		!> Sets the initial values of the FLASH and runs the time steps
		subroutine FLASH_run(FLASH, grid, i_max_time_steps, r_max_time, r_output_step)
            class(t_FLASH), intent(inout)                                 :: FLASH
			type(t_grid), intent(inout)									:: grid
			integer (kind = GRID_SI), intent(inout)						:: i_max_time_steps
			real (kind = GRID_SR), intent(in)							:: r_max_time
			real (kind = GRID_SR), intent(in)							:: r_output_step

			type (t_adaptive_statistics)								:: adaption_stats_initial, adaption_stats_time_steps
			type (t_statistics)									        :: grid_stats_initial, grid_stats_time_steps
			real (kind = GRID_SR)										:: r_t1, r_t2, r_t3, r_t4
			real (kind = GRID_SR)										:: r_time_next_output
			type(t_section_info)           	                            :: grid_info

			!init parameters

			call dgelmt_init

			r_time_next_output = 0.0_GRID_SR

            !$omp master
            if (rank_MPI == 0) then
                _log_write(0, *) "FLASH: setting initial values and a priori refinement.."
                _log_write(0, *) ""
            end if
            !$omp end master

			r_t1 = omp_get_wtime()

			do
				!set numerics and check for refinement
				call FLASH%init%traverse(grid)

                grid_info = grid%get_capacity(.false.)

                !$omp master
                if (rank_MPI == 0) then
                    _log_write(1, "(A, I0, A, I0, A)") " FLASH: ", FLASH%adaption%stats%i_traversals, " adaptions, ", grid_info%i_cells, " cells"
                end if
                !$omp end master

                grid_info = grid%get_capacity(.true.)
				if (FLASH%init%i_refinements_issued .le. grid_info%i_cells / 100) then
					exit
				endif

				call FLASH%adaption%traverse(grid)
			end do

			r_t2 = omp_get_wtime()

            !$omp master
            if (rank_MPI == 0) then
                _log_write(0, *) "FLASH: done."
                _log_write(0, *) ""

                call grid_info%print()
			end if
            !$omp end master

			!output initial grid
			if (r_output_step >= 0.0_GRID_SR) then
				call FLASH%xml_output%traverse(grid)
				r_time_next_output = r_time_next_output + r_output_step
			end if

            !$omp master
			!copy counters
			grid_stats_initial = grid%stats
            adaption_stats_initial = FLASH%adaption%stats

            if (rank_MPI == 0) then
                _log_write(0, *) "FLASH: running time steps.."
                _log_write(0, *) ""
			end if
            !$omp end master

            r_t3 = omp_get_wtime()

			do
				if ((r_max_time >= 0.0 .and. grid%r_time >= r_max_time) .or. (i_max_time_steps >= 0 .and. FLASH%euler%stats%i_traversals >= i_max_time_steps)) then
					exit
				end if

				!if (FLASH%euler%i_refinements_issued > 0) then
					call FLASH%adaption%traverse(grid)
				!end if

				!do a timestep
				call FLASH%euler%traverse(grid)

                grid_info = grid%get_capacity(.false.)

                !$omp master
                if (rank_MPI == 0) then
                    _log_write(1, '(A, I0, A, ES14.7, A, ES14.7, A, I0)') " FLASH: time step: ", FLASH%euler%stats%i_traversals, ", sim. time:", grid%r_time, " s, dt:", grid%r_dt, " s, cells: ", grid_info%i_cells
                end if
                !$omp end master

				!output grid
				if (r_output_step >= 0.0_GRID_SR .and. grid%r_time >= r_time_next_output) then
					call FLASH%xml_output%traverse(grid)
					r_time_next_output = r_time_next_output + r_output_step
				end if
			end do

			r_t4 = omp_get_wtime()

            grid_info = grid%get_capacity(.true.)

            !$omp master
			grid_stats_time_steps = grid%stats - grid_stats_initial
            adaption_stats_time_steps = FLASH%adaption%stats - adaption_stats_initial

            call FLASH%init%stats%reduce()
            call adaption_stats_initial%reduce()
            call grid_stats_initial%reduce()

            call FLASH%euler%stats%reduce()
            call adaption_stats_time_steps%reduce()
            call grid_stats_time_steps%reduce()

            if (rank_MPI == 0) then
                _log_write(0, *) "FLASH: done."
                _log_write(0, *) ""

                _log_write(0, *) "Initialization:"
                _log_write(0, '(A, T34, A)') " Init: ", trim(FLASH%init%stats%to_string())
                _log_write(0, '(A, T34, A)') " Adaptions: ", trim(adaption_stats_initial%to_string())
                _log_write(0, '(A, T34, A)') " Grid: ", trim(grid_stats_initial%to_string())
                _log_write(0, '(A, T34, F10.4, A)') " Element throughput: ", 1.0e-6 * real(grid_stats_initial%i_traversed_cells, GRID_SR) / (r_t2 - r_t1), " M/s"
                _log_write(0, '(A, T34, F10.4, A)') " Memory throughput: ", real(grid_stats_initial%i_traversed_memory, GRID_SR) / ((1024 * 1024 * 1024) * (r_t2 - r_t1)), " GB/s"
                _log_write(0, '(A, T34, F10.4, A)') " Asagi time:", grid_stats_initial%r_asagi_time, " s"
                _log_write(0, '(A, T34, F10.4, A)') " Initialization time:", r_t2 - r_t1, " s"
                _log_write(0, *) ""
                _log_write(0, *) "Execution:"
                _log_write(0, '(A, T34, A)') " Time steps: ", trim(FLASH%euler%stats%to_string())
                _log_write(0, '(A, T34, A)') " Adaptions: ", trim(adaption_stats_time_steps%to_string())
                _log_write(0, '(A, T34, A)') " Grid: ", trim(grid_stats_time_steps%to_string())
                _log_write(0, '(A, T34, F10.4, A)') " Element throughput: ", 1.0e-6 * real(grid_stats_time_steps%i_traversed_cells, GRID_SR) / (r_t4 - r_t3), " M/s"
                _log_write(0, '(A, T34, F10.4, A)') " Memory throughput: ", real(grid_stats_time_steps%i_traversed_memory, GRID_SR) / ((1024 * 1024 * 1024) * (r_t4 - r_t3)), " GB/s"
                _log_write(0, '(A, T34, F10.4, A)') " Cell update throughput: ", 1.0e-6 * real(FLASH%euler%stats%i_traversed_cells, GRID_SR) / (r_t4 - r_t3), " M/s"
                _log_write(0, '(A, T34, F10.4, A)') " Flux solver throughput: ", 1.0e-6 * real(FLASH%euler%stats%i_traversed_edges, GRID_SR) / (r_t4 - r_t3), " M/s"
                _log_write(0, '(A, T34, F10.4, A)') " Asagi time:", grid_stats_time_steps%r_asagi_time, " s"
                _log_write(0, '(A, T34, F10.4, A)') " Execution time:", r_t4 - r_t3, " s"
                _log_write(0, *) ""
                _log_write(0, '(A, T34, F10.4, A)') " Total time:", (r_t2 - r_t1) + (r_t4 - r_t3), " s"
                _log_write(0, *) "---"

                call grid_info%print()
			end if
            !$omp end master
		end subroutine
	END MODULE FLASH
#endif

