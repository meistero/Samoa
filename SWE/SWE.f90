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
            type(t_swe_euler_timestep_traversal)    :: euler
            type(t_swe_adaption_traversal)          :: adaption

            contains

            procedure , pass :: create => swe_create
            procedure , pass :: run => swe_run
            procedure , pass :: destroy => swe_destroy
        end type

		contains

		!> Creates all required runtime objects for the scenario
		subroutine swe_create(swe, grid, l_log, i_asagi_mode)
            class(t_swe), intent(inout)                                 :: swe
			type(t_grid), intent(inout)									:: grid
			logical (kind = GRID_SL), intent(in)						:: l_log
			integer, intent(in)											:: i_asagi_mode

			!local variables
			character(64)												:: s_log_name, s_date, s_time

#           if defined (_SWE_LF)
                character(*), parameter :: s_solver = "Lax Friedrichs"
#           elif defined (_SWE_LLF)
                character(*), parameter :: s_solver = "Local Lax Friedrichs"
#           elif defined(_SWE_LF_BATH)
                character(*), parameter :: s_solver = "Lax Friedrichs + Bathymetry"
#           elif defined(_SWE_LLF_BATH)
                character(*), parameter :: s_solver = "Local Lax Friedrichs + Bathymetry"
#           elif defined(_SWE_FWAVE)
                character(*), parameter :: s_solver = "FWave"
#           elif defined(_SWE_SSQ_FWAVE)
                character(*), parameter :: s_solver = "SSQ-FWave"
#           elif defined(_SWE_AUG_RIEMANN)
                character(*), parameter :: s_solver = "Augmented Riemann"
#           endif

            !$omp master
            if (rank_MPI == 0) then
                _log_write(0, '(" SWE: flux solver: ", A)') s_solver
            end if
            !$omp end master

			!open log file
			call date_and_time(s_date, s_time)
			write (swe%output%s_file_stamp, "(A, A, A8, A, A6)") "output/swe", "_", s_date, "_", s_time
			write (swe%xml_output%s_file_stamp, "(A, A, A8, A, A6)") "output/swe", "_", s_date, "_", s_time
			write (s_log_name, '(A, A)') TRIM(swe%xml_output%s_file_stamp), ".log"

			if (l_log) then
				_log_open_file(s_log_name)
			endif

			!create a quadrature rule
			call t_qr_create_dunavant_rule(qr_Q, max(1, 2 * _SWE_ORDER))

			call load_scenario(grid, i_asagi_mode, "data/displ.nc", "data/bath.nc", 2.0d6, [-0.5d6, -1.0d6])
            !call load_scenario(grid, i_asagi_mode, "data/seissol_displ.nc", "data/seissol_bathymetry.nc", 1.0d5, [-0.5d5, -0.5d5])
		end subroutine

		subroutine load_scenario(grid, i_asagi_mode, ncd_displ, ncd_bath, scaling, offset)
			type(t_grid), target, intent(inout)     :: grid
			integer, intent(in)						:: i_asagi_mode
            character(*), intent(in)                :: ncd_displ, ncd_bath
            double precision, intent(in)            :: scaling, offset(2)

			integer                                 :: i_error, k

#			if defined(_ASAGI)
#               if defined(_ASAGI_NUMA)
                    grid%afh_displacement = grid_create_for_numa(grid_type = GRID_FLOAT, hint = ior(i_asagi_mode, GRID_HAS_TIME), levels = 1, tcount=omp_get_max_threads())
                    grid%afh_bathymetry = grid_create_for_numa(grid_type = GRID_FLOAT, hint = i_asagi_mode, levels = 1, tcount=omp_get_max_threads())

                    !$omp parallel private(i_error)
						i_error = grid_register_thread(grid%afh_displacement); assert_eq(i_error, GRID_SUCCESS)
						i_error = grid_register_thread(grid%afh_bathymetry); assert_eq(i_error, GRID_SUCCESS)
                        i_error = asagi_open(grid%afh_displacement, ncd_displ, 0); assert_eq(i_error, GRID_SUCCESS)
                        i_error = asagi_open(grid%afh_bathymetry, ncd_bath, 0); assert_eq(i_error, GRID_SUCCESS)
                    !$omp end parallel
#               else
                    grid%afh_displacement = asagi_create(grid_type = GRID_FLOAT, hint = ior(i_asagi_mode, GRID_HAS_TIME), levels = 1)
                    grid%afh_bathymetry = asagi_create(grid_type = GRID_FLOAT, hint = i_asagi_mode, levels = 1)

                    i_error = asagi_open(grid%afh_displacement, ncd_displ, 0); assert_eq(i_error, GRID_SUCCESS)
                    i_error = asagi_open(grid%afh_bathymetry, ncd_bath, 0); assert_eq(i_error, GRID_SUCCESS)
#               endif

                if (rank_MPI == 0) then
                    associate(afh_d => grid%afh_displacement, afh_b => grid%afh_bathymetry)
                        _log_write(1, '(" SWE: loaded ", A, ", domain: [", F0.2, ", ", F0.2, "] x [", F0.2, ", ", F0.2, "], time: [", F0.2, ", ", F0.2, "]")') &
                            ncd_displ, grid_min_x(afh_d), grid_max_x(afh_d),  grid_min_y(afh_d), grid_max_y(afh_d),  grid_min_z(afh_d), grid_max_z(afh_d)

                        _log_write(1, '(" SWE: loaded ", A, ", domain: [", F0.2, ", ", F0.2, "] x [", F0.2, ", ", F0.2, "], time: [", F0.2, ", ", F0.2, "]")') &
                            ncd_bath, grid_min_x(afh_b), grid_max_x(afh_b),  grid_min_y(afh_b), grid_max_y(afh_b),  grid_min_z(afh_b), grid_max_z(afh_b)
                    end associate
                end if

                grid%scaling = scaling
                grid%offset = offset
#			endif
		end subroutine

		!> Destroys all required runtime objects for the scenario
		subroutine swe_destroy(swe, grid, l_log)
            class(t_swe), intent(inout)                                  :: swe
			type(t_grid), intent(inout)									:: grid
			logical (kind = GRID_SL)		:: l_log

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

		!> Sets the initial values of the SWE and runs the time steps
		subroutine swe_run(swe, grid, i_max_time_steps, r_max_time, r_output_step)
            class(t_swe), intent(inout)                                 :: swe
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
			r_time_next_output = 0.0_GRID_SR

            !$omp master
            if (rank_MPI == 0) then
                _log_write(0, *) "SWE: setting initial values and a priori refinement.."
                _log_write(0, *) ""
            end if
            !$omp end master

			r_t1 = omp_get_wtime()

			do
				!set numerics and check for refinement
				call swe%init%traverse(grid)

                grid_info = grid%get_capacity(.false.)

                !$omp master
                if (rank_MPI == 0) then
                    _log_write(1, "(A, I0, A, I0, A)") " SWE: ", swe%adaption%stats%i_traversals, " adaptions, ", grid_info%i_cells, " cells"
                end if
                !$omp end master

                grid_info = grid%get_capacity(.true.)
				if (swe%init%i_refinements_issued .le. grid_info%i_cells / 100) then
					exit
				endif

				call swe%adaption%traverse(grid)
			end do

			r_t2 = omp_get_wtime()

            !$omp master
            if (rank_MPI == 0) then
                _log_write(0, *) "SWE: done."
                _log_write(0, *) ""

                call grid_info%print()
			end if
            !$omp end master

			!output initial grid
			if (r_output_step >= 0.0_GRID_SR) then
				call swe%xml_output%traverse(grid)
				r_time_next_output = r_time_next_output + r_output_step
			end if

            !$omp master
			!copy counters
			grid_stats_initial = grid%stats
            adaption_stats_initial = swe%adaption%stats

            if (rank_MPI == 0) then
                _log_write(0, *) "SWE: running time steps.."
                _log_write(0, *) ""
			end if
            !$omp end master

            r_t3 = omp_get_wtime()

#           if defined(_ASAGI)
                ! during the earthquake, do small time steps that include a displacement

                do
                    if ((r_max_time >= 0.0 .and. grid%r_time > r_max_time) .or. (i_max_time_steps >= 0 .and. swe%euler%stats%i_traversals >= i_max_time_steps)) then
                        exit
                    end if

                    if (grid%r_time > grid_max_z(grid%afh_displacement)) then
                        exit
                    end if

                    !do an euler time step
                    call swe%adaption%traverse(grid)
                    call swe%euler%traverse(grid)

                    !displace time-dependent bathymetry
                    call swe%displace%traverse(grid)

                    grid_info = grid%get_capacity(.false.)

                    !$omp master
                    if (rank_MPI == 0) then
                        _log_write(1, '(A, I0, A, ES14.7, A, ES14.7, A, I0)') " SWE: EQ time step: ", swe%euler%stats%i_traversals, ", sim. time:", grid%r_time, " s, dt:", grid%r_dt, " s, cells: ", grid_info%i_cells
                    end if
                    !$omp end master

                    !output grid
                    if (r_output_step >= 0.0_GRID_SR .and. grid%r_time >= r_time_next_output) then
                        call swe%xml_output%traverse(grid)
                        r_time_next_output = r_time_next_output + r_output_step
                    end if
                end do
#           endif

            !regular tsunami time steps begin after the earthquake is over

			do
				if ((r_max_time >= 0.0 .and. grid%r_time > r_max_time) .or. (i_max_time_steps >= 0 .and. swe%euler%stats%i_traversals >= i_max_time_steps)) then
					exit
				end if

				call swe%adaption%traverse(grid)

				!do a time step
				call swe%euler%traverse(grid)

                grid_info = grid%get_capacity(.false.)

                !$omp master
                if (rank_MPI == 0) then
                    _log_write(1, '(A, I0, A, ES14.7, A, ES14.7, A, I0)') " SWE: time step: ", swe%euler%stats%i_traversals, ", sim. time:", grid%r_time, " s, dt:", grid%r_dt, " s, cells: ", grid_info%i_cells
                end if
                !$omp end master

				!output grid
				if (r_output_step >= 0.0_GRID_SR .and. grid%r_time >= r_time_next_output) then
					call swe%xml_output%traverse(grid)
					r_time_next_output = r_time_next_output + r_output_step
				end if
			end do

			r_t4 = omp_get_wtime()

            grid_info = grid%get_capacity(.true.)

            !$omp master
			grid_stats_time_steps = grid%stats - grid_stats_initial
            adaption_stats_time_steps = swe%adaption%stats - adaption_stats_initial

            call swe%init%stats%reduce()
            call adaption_stats_initial%reduce()
            call grid_stats_initial%reduce()

            call swe%euler%stats%reduce()
            call adaption_stats_time_steps%reduce()
            call grid_stats_time_steps%reduce()

            if (rank_MPI == 0) then
                _log_write(0, *) "SWE: done."
                _log_write(0, *) ""

                _log_write(0, *) "Initialization phase:"
                _log_write(0, *) ""
                _log_write(0, '(A, T34, A)') " Init: ", trim(swe%init%stats%to_string())
                _log_write(0, '(A, T34, A)') " Adaptions: ", trim(adaption_stats_initial%to_string())
                _log_write(0, '(A, T34, A)') " Grid: ", trim(grid_stats_initial%to_string())
                _log_write(0, '(A, T34, F10.4, A)') " Element throughput: ", 1.0e-6 * real(grid_stats_initial%i_traversed_cells, GRID_SR) / (r_t2 - r_t1), " M/s"
                _log_write(0, '(A, T34, F10.4, A)') " Memory throughput: ", real(grid_stats_initial%i_traversed_memory, GRID_SR) / ((1024 * 1024 * 1024) * (r_t2 - r_t1)), " GB/s"
                _log_write(0, '(A, T34, F10.4, A)') " Asagi time:", grid_stats_initial%r_asagi_time, " s"
                _log_write(0, '(A, T34, F10.4, A)') " Initialization phase time:", r_t2 - r_t1, " s"
                _log_write(0, *) ""
                _log_write(0, *) "Time Step phase:"
                _log_write(0, *) ""
                _log_write(0, '(A, T34, A)') " Time steps: ", trim(swe%euler%stats%to_string())
                _log_write(0, '(A, T34, A)') " Displace: ", trim(swe%displace%stats%to_string())
                _log_write(0, '(A, T34, A)') " Adaptions: ", trim(adaption_stats_time_steps%to_string())
                _log_write(0, '(A, T34, A)') " Grid: ", trim(grid_stats_time_steps%to_string())
                _log_write(0, '(A, T34, F10.4, A)') " Element throughput: ", 1.0e-6 * real(grid_stats_time_steps%i_traversed_cells, GRID_SR) / (r_t4 - r_t3), " M/s"
                _log_write(0, '(A, T34, F10.4, A)') " Memory throughput: ", real(grid_stats_time_steps%i_traversed_memory, GRID_SR) / ((1024 * 1024 * 1024) * (r_t4 - r_t3)), " GB/s"
                _log_write(0, '(A, T34, F10.4, A)') " Cell update throughput: ", 1.0e-6 * real(swe%euler%stats%i_traversed_cells, GRID_SR) / (r_t4 - r_t3), " M/s"
                _log_write(0, '(A, T34, F10.4, A)') " Flux solver throughput: ", 1.0e-6 * real(swe%euler%stats%i_traversed_edges, GRID_SR) / (r_t4 - r_t3), " M/s"
                _log_write(0, '(A, T34, F10.4, A)') " Asagi time:", grid_stats_time_steps%r_asagi_time, " s"
                _log_write(0, '(A, T34, F10.4, A)') " Time Step phase time:", r_t4 - r_t3, " s"
                _log_write(0, *) ""
                _log_write(0, '(A, T34, F10.4, A)') " Total time:", (r_t2 - r_t1) + (r_t4 - r_t3), " s"
                _log_write(0, *) "---"

                call grid_info%print()
			end if
            !$omp end master
		end subroutine
	END MODULE SWE
#endif

