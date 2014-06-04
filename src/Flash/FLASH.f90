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

			call load_scenario(grid, cfg%s_bathymetry_file, cfg%s_displacement_file)
		end subroutine

	subroutine load_scenario(grid, ncd_bath, ncd_displ, scaling, offset)
			type(t_grid), target, intent(inout)     :: grid
            character(*), intent(in)                :: ncd_bath, ncd_displ
            double precision, optional,intent(in)   :: scaling, offset(2)

			integer						            :: i_asagi_hints
			integer                                 :: i_error

#			if defined(_ASAGI)
                select case(cfg%i_asagi_mode)
                    case (0)
                        i_asagi_hints = GRID_NO_HINT
                    case (1)
                        i_asagi_hints = ieor(GRID_NOMPI, GRID_PASSTHROUGH)
                    case (2)
                        i_asagi_hints = GRID_NOMPI
                    case (3)
                        i_asagi_hints = ieor(GRID_NOMPI, SMALL_CACHE)
                    case (4)
                        i_asagi_hints = GRID_LARGE_GRID
                    case default
                        try(.false., "Invalid asagi mode, must be in range 0 to 4")
                end select

#               if defined(_ASAGI_NUMA)
                    cfg%afh_bathymetry = grid_create_for_numa(grid_type = GRID_FLOAT, hint = i_asagi_hints, levels = 1, tcount=omp_get_max_threads())
                    cfg%afh_displacement = grid_create_for_numa(grid_type = GRID_FLOAT, hint = ior(i_asagi_hints, GRID_HAS_TIME), levels = 1, tcount=omp_get_max_threads())

                    !$omp parallel private(i_error)
						i_error = grid_register_thread(cfg%afh_bathymetry); assert_eq(i_error, GRID_SUCCESS)
						i_error = grid_register_thread(cfg%afh_displacement); assert_eq(i_error, GRID_SUCCESS)
                        i_error = asagi_open(cfg%afh_bathymetry, trim(ncd_bath), 0); assert_eq(i_error, GRID_SUCCESS)
                        i_error = asagi_open(cfg%afh_displacement, trim(ncd_displ), 0); assert_eq(i_error, GRID_SUCCESS)
                    !$omp end parallel
#               else
                    cfg%afh_bathymetry = asagi_create(grid_type = GRID_FLOAT, hint = i_asagi_hints, levels = 1)
                    cfg%afh_displacement = asagi_create(grid_type = GRID_FLOAT, hint = ior(i_asagi_hints, GRID_HAS_TIME), levels = 1)

                    i_error = asagi_open(cfg%afh_bathymetry, trim(ncd_bath), 0); assert_eq(i_error, GRID_SUCCESS)
                    i_error = asagi_open(cfg%afh_displacement, trim(ncd_displ), 0); assert_eq(i_error, GRID_SUCCESS)
#               endif

                associate(afh_d => cfg%afh_displacement, afh_b => cfg%afh_bathymetry)
                    if (present(scaling)) then
                    else
                        cfg%scaling = max(grid_max_x(afh_b) - grid_min_x(afh_b), grid_max_y(afh_b) - grid_min_y(afh_b))
                    end if

                    if (present(offset)) then
                    else
                        cfg%offset = [0.5_GRID_SR * (grid_min_x(afh_d) + grid_max_x(afh_d)), 0.5_GRID_SR * (grid_min_y(afh_d) + grid_max_y(afh_d))] - 0.5_GRID_SR * cfg%scaling
                        cfg%offset = min(max(cfg%offset, [grid_min_x(afh_b), grid_min_y(afh_b)]), [grid_max_x(afh_b), grid_max_y(afh_b)] - cfg%scaling)
                    end if

                    if (rank_MPI == 0) then
                        _log_write(1, '(" FLASH: loaded ", A, ", domain: [", F0.2, ", ", F0.2, "] x [", F0.2, ", ", F0.2, "], time: [", F0.2, ", ", F0.2, "]")') &
                            trim(ncd_bath), grid_min_x(afh_b), grid_max_x(afh_b),  grid_min_y(afh_b), grid_max_y(afh_b),  grid_min_z(afh_b), grid_max_z(afh_b)

                        _log_write(1, '(" FLASH: loaded ", A, ", domain: [", F0.2, ", ", F0.2, "] x [", F0.2, ", ", F0.2, "], time: [", F0.2, ", ", F0.2, "]")') &
                            trim(ncd_displ), grid_min_x(afh_d), grid_max_x(afh_d),  grid_min_y(afh_d), grid_max_y(afh_d),  grid_min_z(afh_d), grid_max_z(afh_d)

                        _log_write(1, '(" FLASH: computational domain: [", F0.2, ", ", F0.2, "] x [", F0.2, ", ", F0.2, "]")'), cfg%offset(1), cfg%offset(1) + cfg%scaling, cfg%offset(2), cfg%offset(2) + cfg%scaling
                    end if
               end associate
#           else
                cfg%scaling = 1.0_GRID_SR
                cfg%offset = [0.0_GRID_SR, 0.0_GRID_SR]
#			endif
		end subroutine

		!> Destroys all required runtime objects for the scenario
		subroutine FLASH_destroy(FLASH, grid, l_log)
            class(t_FLASH), intent(inout)                                  :: FLASH
			type(t_grid), intent(inout)									:: grid
			logical, intent(in)						                    :: l_log

#			if defined(_ASAGI)
				call asagi_close(cfg%afh_displacement)
				call asagi_close(cfg%afh_bathymetry)
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
		subroutine FLASH_run(FLASH, grid)
            class(t_FLASH), intent(inout)                                 :: FLASH
			type(t_grid), intent(inout)									:: grid

			type (t_adaptive_statistics)								:: adaption_stats_initial, adaption_stats_time_steps
			type (t_adaptive_statistics)							    :: grid_stats_initial, grid_stats_time_steps
			real (kind = GRID_SR)										:: r_t1, r_t2, r_t3, r_t4
			real (kind = GRID_SR)										:: r_time_next_output
			type(t_grid_info)           	                            :: grid_info

			!init parameters

			call dgelmt_init

			r_time_next_output = 0.0_GRID_SR

            !$omp master
            if (rank_MPI == 0) then
                _log_write(0, *) "FLASH: setting initial values and a priori refinement.."
                _log_write(0, *) ""
            end if
            !$omp end master

			r_t1 = get_wtime()

			do
				!set numerics and check for refinement
				call FLASH%init%traverse(grid)

                grid_info = grid%get_info(MPI_SUM, .false.)

                !$omp master
                if (rank_MPI == 0) then
                    _log_write(1, "(A, I0, A, I0, A)") " FLASH: ", FLASH%adaption%stats%i_traversals, " adaptions, ", grid_info%i_cells, " cells"
                end if
                !$omp end master

                grid_info = grid%get_info(MPI_SUM, .true.)
				if (FLASH%init%i_refinements_issued .le. grid_info%i_cells / 100) then
					exit
				endif

				call FLASH%adaption%traverse(grid)
			end do

			r_t2 = get_wtime()

            !$omp master
            if (rank_MPI == 0) then
                _log_write(0, *) "FLASH: done."
                _log_write(0, *) ""

                call grid_info%print()
			end if
            !$omp end master

			!output initial grid
			if (cfg%r_output_time_step >= 0.0_GRID_SR) then
				call FLASH%xml_output%traverse(grid)
				r_time_next_output = r_time_next_output + cfg%r_output_time_step
			end if

            !$omp master
			!copy counters
            call FLASH%init%reduce_stats(MPI_SUM, .true.)
            call FLASH%adaption%reduce_stats(MPI_SUM, .true.)
            call grid%reduce_stats(MPI_SUM, .true.)

			grid_stats_initial = grid%stats
            adaption_stats_initial = FLASH%adaption%stats

            if (rank_MPI == 0) then
                _log_write(0, *) "FLASH: running time steps.."
                _log_write(0, *) ""
			end if
            !$omp end master

            r_t3 = get_wtime()

			do
				if ((cfg%r_max_time >= 0.0 .and. grid%r_time >= cfg%r_max_time) .or. (cfg%i_max_time_steps >= 0 .and. FLASH%euler%stats%i_traversals >= cfg%i_max_time_steps)) then
					exit
				end if

				!if (FLASH%euler%i_refinements_issued > 0) then
					call FLASH%adaption%traverse(grid)
				!end if

				!do a timestep
				call FLASH%euler%traverse(grid)

                grid_info = grid%get_info(MPI_SUM, .false.)

                !$omp master
                if (rank_MPI == 0) then
                    _log_write(1, '(A, I0, A, ES14.7, A, ES14.7, A, I0)') " FLASH: time step: ", FLASH%euler%stats%i_traversals, ", sim. time:", grid%r_time, " s, dt:", grid%r_dt, " s, cells: ", grid_info%i_cells
                end if
                !$omp end master

				!output grid
				if (cfg%r_output_time_step >= 0.0_GRID_SR .and. grid%r_time >= r_time_next_output) then
					call FLASH%xml_output%traverse(grid)
					r_time_next_output = r_time_next_output + cfg%r_output_time_step
				end if
			end do

			r_t4 = get_wtime()

            grid_info = grid%get_info(MPI_SUM, .true.)

            !$omp master
            call FLASH%euler%reduce_stats(MPI_SUM, .true.)
            call FLASH%adaption%reduce_stats(MPI_SUM, .true.)
            call grid%reduce_stats(MPI_SUM, .true.)

			grid_stats_time_steps = grid%stats - grid_stats_initial
            adaption_stats_time_steps = FLASH%adaption%stats - adaption_stats_initial

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

