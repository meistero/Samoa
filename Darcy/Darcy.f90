! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy
		use Darcy_data_types
		use Darcy_initialize_pressure
		use Darcy_initialize_saturation
		use Darcy_edge_dummy
		use Darcy_vtk_output
		use Darcy_xml_output
		use Darcy_grad_p
		use Darcy_pressure_solver_cg
		use Darcy_pressure_solver_jacobi
		use Darcy_transport_eq
		use Darcy_permeability
		use Darcy_adapt

		use linear_solver

		use Samoa_darcy

		implicit none

		type t_darcy
            type(t_darcy_init_pressure_traversal)           :: init_pressure
            type(t_darcy_init_saturation_traversal)         :: init_saturation
            type(t_darcy_vtk_output_traversal)              :: vtk_output
            type(t_darcy_xml_output_traversal)              :: xml_output
            type(t_darcy_grad_p_traversal)                  :: grad_p
            type(t_darcy_transport_eq_traversal)            :: transport_eq
            type(t_darcy_permeability_traversal)            :: permeability
            type(t_darcy_adaption_traversal)                :: adaption
            class(t_linear_solver), allocatable             :: pressure_solver

            contains

            procedure , pass :: create => darcy_create
            procedure , pass :: run => darcy_run
            procedure , pass :: destroy => darcy_destroy
        end type

		private
		public t_darcy

		contains

		!> Creates all required runtime objects for the scenario
		subroutine darcy_create(darcy, grid, l_log, i_asagi_mode)
            class(t_darcy)                                              :: darcy
 			type(t_grid), intent(inout)									:: grid
			logical, intent(in)						                    :: l_log
			integer, intent(in)											:: i_asagi_mode

			!local variables
			character (len = 64)										:: s_log_name, s_date, s_time
			integer                                                     :: i_error
            type(t_darcy_pressure_solver_cg)                            :: pressure_solver_cg
            type(t_darcy_pressure_solver_jacobi)                        :: pressure_solver_jacobi

            !allocate solver

 			grid%r_time = 0.0_GRID_SR
            pressure_solver_cg = t_darcy_pressure_solver_cg(cfg%r_epsilon * cfg%r_p0)
            pressure_solver_jacobi = t_darcy_pressure_solver_jacobi(cfg%r_epsilon * cfg%r_p0)

            allocate(darcy%pressure_solver, source=pressure_solver_cg, stat=i_error); assert_eq(i_error, 0)

			!open log file
			call date_and_time(s_date, s_time)
			write (darcy%vtk_output%s_file_stamp, "(A, A, A8, A, A6)") "output/darcy", "_", s_date, "_", s_time
			write (darcy%xml_output%s_file_stamp, "(A, A, A8, A, A6)") "output/darcy", "_", s_date, "_", s_time
			write (s_log_name, '(A, A)') TRIM(darcy%vtk_output%s_file_stamp), ".log"

			if (l_log) then
				_log_open_file(s_log_name)
			endif

			call load_permeability(grid, cfg%s_permeability_file)
		end subroutine

		subroutine load_permeability(grid, s_template)
			type(t_grid), target, intent(inout)		:: grid
			character(256)					        :: s_template

            integer                                 :: i_asagi_hints
			integer									:: i_error, i, j, i_ext_pos
			character(256)					        :: s_file_name

#			if defined(_ASAGI)
                !convert ASAGI mode to ASAGI hints

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
                    cfg%afh_permeability = grid_create_for_numa(grid_type = GRID_FLOAT, hint = i_asagi_hints, levels = 1, tcount=omp_get_max_threads())

					!$omp parallel private(i_error, i, j, i_ext_pos, s_file_name)
						i_error = grid_register_thread(cfg%afh_permeability); assert_eq(i_error, GRID_SUCCESS)

						do j = min(10, cfg%i_max_depth / 2), 0, -1
                            i_ext_pos = index(s_template, ".", .true.)
                            write (s_file_name, fmt = '(A, "_", I0, A)') s_template(: i_ext_pos - 1), 2 ** j, trim(s_template(i_ext_pos :))

                            i_error = asagi_open(cfg%afh_permeability, trim(s_file_name), 0)

                            if (i_error == GRID_SUCCESS) then
                                exit
                            endif
                        end do

	                    assert_eq(i_error, GRID_SUCCESS)

	                    if (rank_MPI == 0) then
	                		associate(afh => cfg%afh_permeability)
	                        	_log_write(1, '(" Darcy: loaded ", A, ", domain: [", F0.2, ", ", F0.2, "] x [", F0.2, ", ", F0.2, "], time: [", F0.2, ", ", F0.2, "]")') &
	                        	    trim(s_file_name), grid_min_x(afh), grid_max_x(afh),  grid_min_y(afh), grid_max_y(afh),  grid_min_z(afh), grid_max_z(afh)
	                		end associate
	                    end if
					!$omp end parallel
#               else
                    cfg%afh_permeability = asagi_create(grid_type = GRID_FLOAT, hint = i_asagi_hints, levels = cfg%i_max_depth / 2 + 1)

                    do i = 0, cfg%i_max_depth / 2
                        do j = i, 0, -1
                            i_ext_pos = index(s_template, ".", .true.)
                            write (s_file_name, fmt = '(A, "_", I0, A)') s_template(: i_ext_pos - 1), 2 ** j, trim(s_template(i_ext_pos :))

                            i_error = asagi_open(cfg%afh_permeability, trim(s_file_name), i)

                            if (i_error == GRID_SUCCESS) then
                                exit
                            endif
                        end do

                        assert_eq(i_error, GRID_SUCCESS)

                        if (rank_MPI == 0) then
                    		associate(afh => cfg%afh_permeability)
                            	_log_write(1, '(" Darcy: loaded ", A, ", domain: [", F0.2, ", ", F0.2, "] x [", F0.2, ", ", F0.2, "], time: [", F0.2, ", ", F0.2, "]")') &
                            	    trim(s_file_name), grid_min_x(afh), grid_max_x(afh),  grid_min_y(afh), grid_max_y(afh),  grid_min_z(afh), grid_max_z(afh)
                    		end associate
                        end if
                    end do
#               endif
#			endif

            cfg%scaling = 1.0_GRID_SR
            cfg%offset = [0.0_GRID_SR, 0.0_GRID_SR]
		end subroutine

		!> Destroys all required runtime objects for the scenario
		subroutine darcy_destroy(darcy, grid, l_log)
            class(t_darcy)                                               :: darcy
 			type(t_grid), intent(inout)							:: grid
			logical		                    :: l_log
			integer (kind = 1)				:: i

			if (l_log) then
				_log_close_file()
			endif

#			if defined(_ASAGI)
				call asagi_close(cfg%afh_permeability)
#			endif
		end subroutine

		!> Sets the initial values of the scenario and runs the time steps
		subroutine darcy_run(darcy, grid, i_max_time_steps, r_max_time, r_output_step)
            class(t_darcy)                                              :: darcy
 			type(t_grid), intent(inout)									:: grid
			integer (kind = GRID_SI), intent(in)						:: i_max_time_steps
			real (kind = GRID_SR), intent(in)							:: r_max_time
			real (kind = GRID_SR), intent(in)							:: r_output_step

			type (t_adaptive_statistics)								:: adaption_stats_initial, adaption_stats_time_steps
			type (t_statistics)									        :: pressure_solver_stats_initial, pressure_solver_stats_time_steps
			type (t_statistics)									        :: grid_stats_initial, grid_stats_time_steps
            integer (kind = GRID_SI)									:: i_lse_iterations, i_lse_iterations_initial
			double precision										    :: r_t1, r_t2, r_t3, r_t4
			real (kind = GRID_SR)										:: r_time_next_output
			type(t_section_info)           	                            :: grid_info

			!init parameters
			r_time_next_output = 0.0_GRID_SR

            !$omp master
            if (rank_MPI == 0) then
                _log_write(0, *) "Darcy: setting initial values and solving initial system.."
                _log_write(0, *) ""
            end if
            !$omp end master

			r_t1 = omp_get_wtime()

			!set pressure initial condition
			call darcy%init_pressure%traverse(grid)

			do
				!reset saturation to initial condition, compute permeability and set refinement flag
				call darcy%init_saturation%traverse(grid)

				!solve pressure equation
				i_lse_iterations = darcy%pressure_solver%solve(grid)

                grid_info = grid%get_capacity(.false.)

                !$omp master
                if (rank_MPI == 0) then
                    _log_write(1, "(A, I0, A, I0, A, I0, A)") " Darcy: ", darcy%adaption%stats%i_traversals, " adaptions, ", i_lse_iterations, " iterations, ", grid_info%i_cells, " cells"
                end if
                !$omp end master

                grid_info = grid%get_capacity(.true.)
				if (darcy%init_saturation%i_refinements_issued .le. grid_info%i_cells / 100_GRID_DI) then
					exit
				endif

				!refine grid
				call darcy%adaption%traverse(grid)
			end do

			r_t2 = omp_get_wtime()

            !$omp master
            if (rank_MPI == 0) then
                _log_write(0, *) "Darcy: done."
                _log_write(0, *) ""

                call grid_info%print()
            end if
			!$omp end master

			!output initial grid
			if (r_output_step >= 0.0_GRID_SR) then
				call darcy%grad_p%traverse(grid)
				call darcy%xml_output%traverse(grid)
				r_time_next_output = r_time_next_output + r_output_step
			end if

            !$omp master
			!copy counters
			adaption_stats_initial = darcy%adaption%stats
            pressure_solver_stats_initial = darcy%pressure_solver%stats
			grid_stats_initial = grid%stats

            if (rank_MPI == 0) then
                _log_write(0, *) "Darcy: running time steps.."
                _log_write(0, *) ""
            end if
			!$omp end master

			r_t3 = omp_get_wtime()

			do
				if ((r_max_time >= 0.0 .and. grid%r_time > r_max_time) .or. (i_max_time_steps >= 0 .and. darcy%transport_eq%stats%i_traversals >= i_max_time_steps)) then
					exit
				end if

				!refine grid if necessary
				!if (darcy%permeability%i_refinements_issued > 0) then
					!refine grid
					call darcy%adaption%traverse(grid)

					!recompute permeability + refinement flag
					call darcy%permeability%traverse(grid)
				!end if

				!solve pressure equation
				i_lse_iterations = darcy%pressure_solver%solve(grid)

				!compute velocity field
				call darcy%grad_p%traverse(grid)

				!transport equation time step
				call darcy%transport_eq%traverse(grid)

				!compute permeability field + refinement flag
				call darcy%permeability%traverse(grid)

                grid_info = grid%get_capacity(.false.)

                !$omp master
                if (rank_MPI == 0) then
                    _log_write(1, '(A, I0, A, ES14.7, A, ES14.7, A, I0, A, I0)') " Darcy: time step: ", darcy%transport_eq%stats%i_traversals, ", sim. time:", grid%r_time, " s, dt:", grid%r_dt, " s, cells: ", grid_info%i_cells, ", LSE iterations: ", i_lse_iterations
                end if
                !$omp end master

				!output grid
				if (r_output_step >= 0.0_GRID_SR .and. grid%r_time >= r_time_next_output) then
					call darcy%xml_output%traverse(grid)
					r_time_next_output = r_time_next_output + r_output_step
				end if
			end do

			r_t4 = omp_get_wtime()

            grid_info = grid%get_capacity(.true.)

            !$omp master
            adaption_stats_time_steps = darcy%adaption%stats - adaption_stats_initial
            pressure_solver_stats_time_steps = darcy%pressure_solver%stats - pressure_solver_stats_initial
			grid_stats_time_steps = grid%stats - grid_stats_initial

            call darcy%init_saturation%stats%reduce()
            call adaption_stats_initial%reduce()
            call pressure_solver_stats_initial%reduce()
            call grid_stats_initial%reduce()

            call darcy%grad_p%stats%reduce()
            call darcy%transport_eq%stats%reduce()
            call darcy%permeability%stats%reduce()
            call adaption_stats_time_steps%reduce()
            call pressure_solver_stats_time_steps%reduce()
            call grid_stats_time_steps%reduce()

            if (rank_MPI == 0) then
                _log_write(0, *) "Darcy: done."
                _log_write(0, *) ""

                _log_write(0, *) "Initialization phase:"
                _log_write(0, *) ""
                _log_write(0, '(A, T34, A)') " Init: ", trim(darcy%init_saturation%stats%to_string())
                _log_write(0, '(A, T34, A)') " Adaptions: ", trim(adaption_stats_initial%to_string())
                _log_write(0, '(A, T34, A)') " Pressure Solver: ", trim(pressure_solver_stats_initial%to_string())
                _log_write(0, '(A, T34, A)') " Grid: ", trim(grid_stats_initial%to_string())
                _log_write(0, '(A, T34, F10.4, A)') " Element throughput: ", 1.0e-6 * real(grid_stats_initial%i_traversed_cells, GRID_SR) / (r_t2 - r_t1), " M/s"
                _log_write(0, '(A, T34, F10.4, A)') " Memory throughput: ", real(grid_stats_initial%i_traversed_memory, GRID_SR) / ((1024 * 1024 * 1024) * (r_t2 - r_t1)), " GB/s"
                _log_write(0, '(A, T34, F10.4, A)') " Asagi time:", grid_stats_initial%r_asagi_time, " s"
                _log_write(0, '(A, T34, F10.4, A)') " Initialization phase time:", r_t2 - r_t1, " s"
                _log_write(0, *) ""
                _log_write(0, *) "Time step phase:"
                _log_write(0, *) ""
                _log_write(0, '(A, T34, A)') " Transport: ", trim(darcy%transport_eq%stats%to_string())
                _log_write(0, '(A, T34, A)') " Gradient: ", trim(darcy%grad_p%stats%to_string())
                _log_write(0, '(A, T34, A)') " Permeability: ", trim(darcy%permeability%stats%to_string())
                _log_write(0, '(A, T34, A)') " Adaptions: ", trim(adaption_stats_time_steps%to_string())
                _log_write(0, '(A, T34, A)') " Pressure Solver: ", trim(pressure_solver_stats_time_steps%to_string())
                _log_write(0, '(A, T34, A)') " Grid: ", trim(grid_stats_time_steps%to_string())
                _log_write(0, '(A, T34, F10.4, A)') " Element throughput: ", 1.0e-6 * real(grid_stats_time_steps%i_traversed_cells, GRID_SR) / (r_t4 - r_t3), " M/s"
                _log_write(0, '(A, T34, F10.4, A)') " Memory throughput: ", real(grid_stats_time_steps%i_traversed_memory, GRID_SR) / ((1024 * 1024 * 1024) * (r_t4 - r_t3)), " GB/s"
                _log_write(0, '(A, T34, F10.4, A)') " Asagi time:", grid_stats_time_steps%r_asagi_time, " s"
                _log_write(0, '(A, T34, F10.4, A)') " Time step phase time:", r_t4 - r_t3, " s"
                _log_write(0, *) ""
                _log_write(0, '(A, T34, F10.4, A)') " Total time:", (r_t2 - r_t1) + (r_t4 - r_t3), " s"
                _log_write(0, *) "---"

                call grid_info%print()
            end if
            !$omp end master
		end subroutine
	END MODULE Darcy
#endif

