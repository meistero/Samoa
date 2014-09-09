! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
	MODULE SWE
		use SFC_edge_traversal
		use SWE_data_types

		use SWE_adapt
		use SWE_initialize
		use SWE_displace
		use SWE_output
		use SWE_xml_output
		use SWE_ascii_output
		use SWE_point_output
		use SWE_euler_timestep

		use Samoa_swe

		implicit none

		PRIVATE
		PUBLIC t_swe

		type t_swe
            type(t_swe_init_traversal)              :: init
            type(t_swe_displace_traversal)          :: displace
            type(t_swe_output_traversal)            :: output
            type(t_swe_xml_output_traversal)        :: xml_output
            type(t_swe_ascii_output_traversal)      :: ascii_output                     !-------------------------
	    type(t_swe_point_output_traversal)	    :: point_output

            type(t_swe_euler_timestep_traversal)    :: euler
            type(t_swe_adaption_traversal)          :: adaption

            contains

            procedure, pass :: create => swe_create
            procedure, pass :: run => swe_run
            procedure, pass :: destroy => swe_destroy
        end type

		contains

		!> Creates all required runtime objects for the scenario
		subroutine swe_create(swe, grid, l_log, i_asagi_mode)
            class(t_swe), intent(inout)                                 :: swe
			type(t_grid), intent(inout)									:: grid
			logical, intent(in)						                    :: l_log
			integer, intent(in)											:: i_asagi_mode

			!local variables
			character(64)												:: s_log_name, s_date, s_time

			!open log file
			call date_and_time(s_date, s_time)
			write (swe%output%s_file_stamp, "(A, A, A8, A, A6)") "output/swe", "_", s_date, "_", s_time
			write (swe%xml_output%s_file_stamp, "(A, A, A8, A, A6)") "output/swe", "_", s_date, "_", s_time
            write (swe%point_output%s_file_stamp, "(A, A, A8, A, A6)") "output/swe", "_", s_date, "_", s_time
			write (s_log_name, '(A, A)') TRIM(swe%xml_output%s_file_stamp), ".log"

			if (l_log) then
				_log_open_file(s_log_name)
			endif

			call load_scenario(grid, cfg%s_bathymetry_file, cfg%s_displacement_file)

			call swe%init%create()
            call swe%displace%create()
            call swe%output%create()
            call swe%xml_output%create()
            call swe%ascii_output%create()
            call swe%euler%create()
            call swe%adaption%create()
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
                        _log_write(1, '(" SWE: loaded ", A, ", domain: [", F0.2, ", ", F0.2, "] x [", F0.2, ", ", F0.2, "], time: [", F0.2, ", ", F0.2, "]")') &
                            trim(ncd_bath), grid_min_x(afh_b), grid_max_x(afh_b),  grid_min_y(afh_b), grid_max_y(afh_b),  grid_min_z(afh_b), grid_max_z(afh_b)

                        _log_write(1, '(" SWE: loaded ", A, ", domain: [", F0.2, ", ", F0.2, "] x [", F0.2, ", ", F0.2, "], time: [", F0.2, ", ", F0.2, "]")') &
                            trim(ncd_displ), grid_min_x(afh_d), grid_max_x(afh_d),  grid_min_y(afh_d), grid_max_y(afh_d),  grid_min_z(afh_d), grid_max_z(afh_d)

                        _log_write(1, '(" SWE: computational domain: [", F0.2, ", ", F0.2, "] x [", F0.2, ", ", F0.2, "]")'), cfg%offset(1), cfg%offset(1) + cfg%scaling, cfg%offset(2), cfg%offset(2) + cfg%scaling
                    end if
               end associate
#           else
                cfg%scaling = 1.0_GRID_SR
                cfg%offset = [0.0_GRID_SR, 0.0_GRID_SR]
#			endif
		end subroutine

		!> Destroys all required runtime objects for the scenario
		subroutine swe_destroy(swe, grid, l_log)
            class(t_swe), intent(inout)                                  :: swe
			type(t_grid), intent(inout)									:: grid
			logical, intent(in)		:: l_log

			call swe%init%destroy()
            call swe%displace%destroy()
            call swe%output%destroy()
            call swe%xml_output%destroy()
            call swe%ascii_output%destroy()
            call swe%point_output%destroy()
            call swe%euler%destroy()
            call swe%adaption%destroy()

#			if defined(_ASAGI)
				call asagi_close(cfg%afh_displacement)
				call asagi_close(cfg%afh_bathymetry)
#			endif

			if (l_log) then
				_log_close_file()
			endif
		end subroutine

		!*********************************
		! run()-method
		!*********************************

		!> Sets the initial values of the SWE and runs the time steps
		subroutine swe_run(swe, grid, i_max_time_steps, r_max_time, r_output_step)
            class(t_swe), intent(inout)                                 :: swe
			type(t_grid), intent(inout)									:: grid
			integer (kind = GRID_SI), intent(inout)						:: i_max_time_steps
			real (kind = GRID_SR), intent(in)							:: r_max_time
			real (kind = GRID_SR), intent(in)							:: r_output_step

			double precision										    :: t_initial, t_time_steps
			real (kind = GRID_SR)										:: r_time_next_output
			type(t_grid_info)           	                            :: grid_info, grid_info_max
			integer (kind = GRID_SI)                                    :: i_initial_step, i_time_step
			type (t_adaptive_statistics)								:: adaption_stats_initial, adaption_stats_time_steps
			type (t_adaptive_statistics)							    :: grid_stats_initial, grid_stats_time_steps

			!init parameters
			r_time_next_output = 0.0_GRID_SR

            if (rank_MPI == 0) then
                !$omp master
                _log_write(0, *) "SWE: setting initial values and a priori refinement.."
                _log_write(0, *) ""
                !$omp end master
            end if

			t_initial = -get_wtime()
            i_initial_step = 0

			do
				!set numerics and check for refinement
				call swe%init%traverse(grid)

                if (rank_MPI == 0) then
                    grid_info%i_cells = grid%get_cells(MPI_SUM, .false.)

                    !$omp master
                    _log_write(1, "(A, I0, A, I0, A)") " SWE: ", i_initial_step, " adaptions, ", grid_info%i_cells, " cells"
                    !$omp end master
                end if

                grid_info%i_cells = grid%get_cells(MPI_SUM, .true.)
				if (swe%init%i_refinements_issued .le. grid_info%i_cells / 100_GRID_DI) then
					exit
				endif

				call swe%adaption%traverse(grid)

				i_initial_step = i_initial_step + 1
			end do

			t_initial = t_initial + get_wtime()

            grid_info = grid%get_info(MPI_SUM, .true.)

            if (rank_MPI == 0) then
                !$omp master
                _log_write(0, *) "SWE: done."
                _log_write(0, *) ""

                call grid_info%print()
                !$omp end master
			end if

			!output initial grid
			if (r_output_step >= 0.0_GRID_SR) then
                if (cfg%l_ascii_output) then
                    call swe%ascii_output%traverse(grid)
                end if
		if(cfg%l_gridoutput) then
                    call swe%xml_output%traverse(grid)
                end if
		if (cfg%l_pointoutput) then
		    call swe%point_output%traverse(grid)
		end if

				r_time_next_output = r_time_next_output + r_output_step
			end if

            !$omp master
            call swe%init%reduce_stats(MPI_SUM, .true.)
            call swe%adaption%reduce_stats(MPI_SUM, .true.)
            call grid%reduce_stats(MPI_SUM, .true.)

			!copy counters
            adaption_stats_initial = swe%adaption%stats
            grid_stats_initial = grid%stats

            if (rank_MPI == 0) then
                _log_write(0, *) "SWE: running time steps.."
                _log_write(0, *) ""
			end if
            !$omp end master

            t_time_steps = -get_wtime()
            i_time_step = 0

#           if defined(_ASAGI)
                ! during the earthquake, do small time steps that include a displacement

                do
                    if ((r_max_time >= 0.0 .and. grid%r_time > r_max_time) .or. (i_max_time_steps >= 0 .and. i_time_step >= i_max_time_steps)) then
                        exit
                    end if

                    if (grid%r_time > grid_max_z(cfg%afh_displacement)) then
                        exit
                    end if

                    !do an euler time step
                    call swe%adaption%traverse(grid)

                    call swe%euler%traverse(grid)
                    i_time_step = i_time_step + 1

                    !displace time-dependent bathymetry
                    call swe%displace%traverse(grid)

                    if (rank_MPI == 0) then
                        grid_info%i_cells = grid%get_cells(MPI_SUM, .false.)

                        !$omp master
                        _log_write(1, '(A, I0, A, ES14.7, A, ES14.7, A, I0)') " SWE: EQ time step: ", i_time_step, ", sim. time:", grid%r_time, " s, dt:", grid%r_dt, " s, cells: ", grid_info%i_cells
                        !$omp end master
                    end if

                    !output grid
                    if (r_output_step >= 0.0_GRID_SR .and. grid%r_time >= r_time_next_output) then
                        if (cfg%l_ascii_output) then
        	            call swe%ascii_output%traverse(grid)
        	        end if
			if(cfg%l_gridoutput) then
        	            call swe%xml_output%traverse(grid)
        	        end if
			if (cfg%l_pointoutput) then
			    call swe%point_output%traverse(grid)
			end if

                        r_time_next_output = r_time_next_output + r_output_step
                    end if
                end do
#           endif

            !regular tsunami time steps begin after the earthquake is over

			do
				if ((r_max_time >= 0.0 .and. grid%r_time > r_max_time) .or. (i_max_time_steps >= 0 .and. i_time_step >= i_max_time_steps)) then
					exit
				end if

				call swe%adaption%traverse(grid)

				!do a time step
				call swe%euler%traverse(grid)
				i_time_step = i_time_step + 1

                if (rank_MPI == 0) then
                    grid_info%i_cells = grid%get_cells(MPI_SUM, .false.)

                    !$omp master
                    _log_write(1, '(A, I0, A, ES14.7, A, ES14.7, A, I0)') " SWE: time step: ", i_time_step, ", sim. time:", grid%r_time, " s, dt:", grid%r_dt, " s, cells: ", grid_info%i_cells
                    !$omp end master
                end if

				!output grid
				if (r_output_step >= 0.0_GRID_SR .and. grid%r_time >= r_time_next_output) then
                    if (cfg%l_ascii_output) then
             	       call swe%ascii_output%traverse(grid)
               	    end if
		    if(cfg%l_gridoutput) then
                    	call swe%xml_output%traverse(grid)
                    end if
		    if (cfg%l_pointoutput) then
		    	call swe%point_output%traverse(grid)
		    end if

					r_time_next_output = r_time_next_output + r_output_step
				end if
			end do

			t_time_steps = t_time_steps + get_wtime()

            grid_info = grid%get_info(MPI_SUM, .true.)
            grid_info_max = grid%get_info(MPI_MAX, .true.)

            !$omp master
            call swe%euler%reduce_stats(MPI_SUM, .true.)
            call swe%displace%reduce_stats(MPI_SUM, .true.)
            call swe%adaption%reduce_stats(MPI_SUM, .true.)
            if (cfg%l_pointoutput == .true.) then
                call swe%point_output%reduce_stats(MPI_SUM, .true.)
            end if
            call grid%reduce_stats(MPI_SUM, .true.)

			grid_stats_time_steps = grid%stats - grid_stats_initial
            adaption_stats_time_steps = swe%adaption%stats - adaption_stats_initial

            if (rank_MPI == 0) then
                _log_write(0, *) "SWE: done."
                _log_write(0, *) ""

                _log_write(0, *) "Initialization phase:"
                _log_write(0, *) ""
                _log_write(0, '(A, T34, A)') " Init: ", trim(swe%init%stats%to_string())
                _log_write(0, '(A, T34, A)') " Adaptions: ", trim(adaption_stats_initial%to_string())
                _log_write(0, '(A, T34, A)') " Grid: ", trim(grid_stats_initial%to_string())
                _log_write(0, '(A, T34, F10.4, A)') " Element throughput: ", 1.0d-6 * dble(grid_stats_initial%i_traversed_cells) / t_initial, " M/s"
                _log_write(0, '(A, T34, F10.4, A)') " Memory throughput: ", dble(grid_stats_initial%i_traversed_memory) / ((1024 * 1024 * 1024) * t_initial), " GB/s"
                _log_write(0, '(A, T34, F10.4, A)') " Asagi time:", grid_stats_initial%r_asagi_time, " s"
                _log_write(0, '(A, T34, F10.4, A)') " Initialization phase time:", t_initial, " s"
                _log_write(0, *) ""
                _log_write(0, *) "Time step phase:"
                _log_write(0, *) ""
                _log_write(0, '(A, T34, A)') " Time steps: ", trim(swe%euler%stats%to_string())
                _log_write(0, '(A, T34, A)') " Displace: ", trim(swe%displace%stats%to_string())
                _log_write(0, '(A, T34, A)') " Adaptions: ", trim(adaption_stats_time_steps%to_string())
                if (cfg%l_pointoutput == .true.) then
                    _log_write(0, '(A, T34, A)') " Point output: ", trim(swe%point_output%stats%to_string())
                end if
                _log_write(0, '(A, T34, A)') " Grid: ", trim(grid_stats_time_steps%to_string())
                _log_write(0, '(A, T34, F10.4, A)') " Element throughput: ", 1.0d-6 * dble(grid_stats_time_steps%i_traversed_cells) / t_time_steps, " M/s"
                _log_write(0, '(A, T34, F10.4, A)') " Memory throughput: ", dble(grid_stats_time_steps%i_traversed_memory) / ((1024 * 1024 * 1024) * t_time_steps), " GB/s"
                _log_write(0, '(A, T34, F10.4, A)') " Cell update throughput: ", 1.0d-6 * dble(swe%euler%stats%i_traversed_cells) / t_time_steps, " M/s"
                _log_write(0, '(A, T34, F10.4, A)') " Flux solver throughput: ", 1.0d-6 * dble(swe%euler%stats%i_traversed_edges) / t_time_steps, " M/s"
                _log_write(0, '(A, T34, F10.4, A)') " Asagi time:", grid_stats_time_steps%r_asagi_time, " s"
                _log_write(0, '(A, T34, F10.4, A)') " Time step phase time:", t_time_steps, " s"
                _log_write(0, *) ""
                _log_write(0, '(A, T34, F10.4, A)') " Total time:", t_initial + t_time_steps, " s"
                _log_write(0, *) "---"

                _log_write(0, '("  Cells: avg: ", I0, " max: ", I0)') grid_info%i_cells / (omp_get_max_threads() * size_MPI), grid_info_max%i_cells
                _log_write(0, '()')

                call grid_info%print()
			end if
            !$omp end master
		end subroutine
	END MODULE SWE
#endif

