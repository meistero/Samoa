! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_NUMA)
	MODULE NUMA
		use SFC_edge_traversal
		use NUMA_data_types

		use NUMA_adapt
		use NUMA_initialize
		use NUMA_output
		use NUMA_xml_output
		use NUMA_euler_timestep

		use Samoa_numa

		implicit none

		PRIVATE
		PUBLIC t_numa

		type t_numa
            type(t_numa_init_traversal)              :: init
            type(t_numa_output_traversal)            :: output
            type(t_numa_xml_output_traversal)        :: xml_output
            type(t_numa_euler_timestep_traversal)    :: euler
            type(t_numa_adaption_traversal)          :: adaption

            contains

            procedure , pass :: create => numa_create
            procedure , pass :: run => numa_run
            procedure , pass :: destroy => numa_destroy
        end type

		contains

		!> Creates all required runtime objects for the scenario
		subroutine numa_create(numa, grid, l_log)
            class(t_numa), intent(inout)                                  :: numa
			type(t_grid), intent(inout)									:: grid
			logical (kind = GRID_SL)									:: l_log

			!local variables
			character(64)												:: s_log_name, s_date, s_time

			!open log file
			call date_and_time(s_date, s_time)
			write (numa%output%s_file_stamp, "(A, A, A8, A, A6)") "output/numa", "_", s_date, "_", s_time
			write (numa%xml_output%s_file_stamp, "(A, A, A8, A, A6)") "output/numa", "_", s_date, "_", s_time
			write (s_log_name, '(A, A)') TRIM(numa%xml_output%s_file_stamp), ".log"

			if (l_log) then
				_log_open_file(s_log_name)
			endif

			!create a quadrature rule
			call t_qr_create_dunavant_rule(qr_Q, max(1, 2 * _NUMA_ORDER))

			call load_external_data(grid)
		end subroutine

		subroutine load_external_data(grid)
			type(t_grid), target, intent(inout)     :: grid
			integer                                 :: i_error, k
			integer, pointer						:: afh

#			if defined(_ASAGI)
				grid%afh_displacement = asagi_create(grid_type = GRID_FLOAT, hint = ieor(GRID_NOMPI, SMALL_CACHE), levels = omp_get_max_threads())
				grid%afh_bathymetry = asagi_create(grid_type = GRID_FLOAT, hint = ieor(GRID_NOMPI, SMALL_CACHE), levels = omp_get_max_threads())

				do k = 0, omp_get_max_threads() - 1
                    i_error = asagi_open(grid%afh_displacement, "data/displ.nc", k); assert_eq(i_error, GRID_SUCCESS)
                    i_error = asagi_open(grid%afh_bathymetry, "data/bath.nc", k); assert_eq(i_error, GRID_SUCCESS)
                end do

                if (rank_MPI == 0) then
                    afh => grid%afh_displacement
                    _log_write(1, '(A, A, A, F0.2, A, F0.2, A, F0.2, A, F0.2, A)') " NUMA: loaded '", "data/displ.nc", "', coordinate system: [", grid_min_x(afh), ", ", grid_min_y(afh), "] x [", grid_max_x(afh), ", ", grid_max_y(afh), "]"

                    afh => grid%afh_bathymetry
                    _log_write(1, '(A, A, A, F0.2, A, F0.2, A, F0.2, A, F0.2, A)') " NUMA: loaded '", "data/bath.nc", "', coordinate system: [", grid_min_x(afh), ", ", grid_min_y(afh), "] x [", grid_max_x(afh), ", ", grid_max_y(afh), "]"
                end if
#			endif
		end subroutine

		!> Destroys all required runtime objects for the scenario
		subroutine numa_destroy(numa, grid, l_log)
            class(t_numa), intent(inout)                                  :: numa
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

		!> Sets the initial values of the NUMA and runs the time steps
		subroutine numa_run(numa, grid, i_max_time_steps, r_max_time, r_output_step)
            class(t_numa), intent(inout)                                 :: numa
			type(t_grid), intent(inout)									:: grid
			integer (kind = GRID_SI), intent(inout)						:: i_max_time_steps
			real (kind = GRID_SR), intent(in)							:: r_max_time
			real (kind = GRID_SR), intent(in)							:: r_output_step

			type (t_adaptive_statistics)								:: adaption_stats_initial, adaption_stats_time_steps
			type (t_statistics)									        :: grid_stats_initial, grid_stats_time_steps
			real (kind = GRID_SR)										:: r_t1, r_t2, r_t3, r_t4
			real (kind = GRID_SR)										:: r_time_next_output
			type(t_section_info)           	                            :: grid_info
		!	print *, "I WAS HERE IN NUMA.f90-------------"
			!init parameters
			r_time_next_output = 0.0_GRID_SR

            !$omp master
            if (rank_MPI == 0) then
                _log_write(0, '(A, I0, A, I0, A, I0)') " NUMA: mpi ranks: ", size_MPI, ", openmp threads: ", omp_get_max_threads(), ", openmp procs: ", omp_get_num_procs()
                _log_write(0, '(A, I0, A, I0, A, I0)') " NUMA: min depth: ", grid%i_min_depth, ", max depth: ", grid%i_max_depth, ", sections per thread: ", grid%i_sections_per_thread
                _log_write(0, '(A, I0, A, ES9.2, A, ES9.2)') " NUMA: max time steps: ", i_max_time_steps, ", max sim. time: ", r_max_time, ", output step: ", r_output_step
                _log_write(0, *) ""

                _log_write(0, *) "NUMA: setting initial values and a priori refinement.."
                _log_write(0, *) ""
            end if
            !$omp end master

			r_t1 = omp_get_wtime()

			do
				!set numerics and check for refinement
				call numa%init%traverse(grid)

                grid_info = grid%get_capacity(.false.)

                !$omp master
                if (rank_MPI == 0) then
                    _log_write(1, "(A, I0, A, I0, A)") " NUMA: ", numa%adaption%stats%i_traversals, " adaptions, ", grid_info%i_cells, " cells"
                end if
                !$omp end master

                grid_info = grid%get_capacity(.true.)
				if (numa%init%i_refinements_issued .le. grid_info%i_cells / 100) then
					exit
				endif

				call numa%adaption%traverse(grid)
			end do

			r_t2 = omp_get_wtime()

            !$omp master
            if (rank_MPI == 0) then
                _log_write(0, *) "NUMA: done."
                _log_write(0, *) ""

                call grid_info%print()
			end if
            !$omp end master

			!output initial grid
			if (r_output_step >= 0.0_GRID_SR) then
				 _log_write(0, *) "-----OUTPUT TRAVERSAL CALLED-----"
				call numa%xml_output%traverse(grid)
				r_time_next_output = r_time_next_output + r_output_step
			end if

            !$omp master
			!copy counters
			grid_stats_initial = grid%stats
            adaption_stats_initial = numa%adaption%stats

            if (rank_MPI == 0) then
                _log_write(0, *) "NUMA: running time steps.."
                _log_write(0, *) ""
			end if
            !$omp end master

            r_t3 = omp_get_wtime()
			do

				if ((r_max_time >= 0.0 .and. grid%r_time >= r_max_time) .or. (i_max_time_steps >= 0 .and. numa%euler%stats%i_traversals >= i_max_time_steps)) then
				_log_write(0, *) " EXIT CALLED!"	
				_log_write(0, *) (r_max_time >= 0.0) 
				_log_write(0, *) (r_max_time)
				_log_write(0, *) (grid%r_time >= r_max_time)
				_log_write(0, *) (grid%r_time)
				_log_write(0, *) (r_max_time)
				_log_write(0, *) (i_max_time_steps >= 0)
				_log_write(0, *) (i_max_time_steps)
				_log_write(0, *) (numa%euler%stats%i_traversals >= i_max_time_steps)
				_log_write(0, *) (numa%euler%stats%i_traversals)
				_log_write(0, *) (i_max_time_steps)
					exit
				end if
				if (numa%euler%i_refinements_issued > 0) then
					call numa%adaption%traverse(grid)
				end if

				!do a timestep

				call numa%euler%traverse(grid)
                grid_info = grid%get_capacity(.false.)

                !$omp master
                if (rank_MPI == 0) then
                    _log_write(1, '(A, I0, A, ES14.7, A, ES14.7, A, I0)') " NUMA: time step: ", numa%euler%stats%i_traversals, ", sim. time:", grid%r_time, " s, dt:", grid%r_dt, " s, cells: ", grid_info%i_cells
                end if
                !$omp end master

				!output grid
			_log_write(0, *) "OUTPUT-------------------------->",r_output_step,grid%r_time ,r_time_next_output , (r_output_step >= 0.0_GRID_SR .and. grid%r_time >= r_time_next_output)
				if (r_output_step >= 0.0_GRID_SR .and. grid%r_time >= r_time_next_output) then
					call numa%xml_output%traverse(grid)
					r_time_next_output = r_time_next_output + r_output_step
				end if
			end do

			r_t4 = omp_get_wtime()

            grid_info = grid%get_capacity(.true.)

            !$omp master
			grid_stats_time_steps = grid%stats - grid_stats_initial
            adaption_stats_time_steps = numa%adaption%stats - adaption_stats_initial

            call numa%init%stats%reduce()
            call adaption_stats_initial%reduce()
            call grid_stats_initial%reduce()

            call numa%euler%stats%reduce()
            call adaption_stats_time_steps%reduce()
            call grid_stats_time_steps%reduce()

            if (rank_MPI == 0) then
                _log_write(0, *) "NUMA: done."
                _log_write(0, *) ""

                _log_write(0, *) "Initialization:"
                _log_write(0, '(A, T34, A)') " Init: ", trim(numa%init%stats%to_string())
                _log_write(0, '(A, T34, A)') " Adaptions: ", trim(adaption_stats_initial%to_string())
                _log_write(0, '(A, T34, A)') " Grid: ", trim(grid_stats_initial%to_string())
                _log_write(0, '(A, T34, F10.4, A)') " Element throughput: ", 1.0e-6 * real(grid_stats_initial%i_traversed_cells, GRID_SR) / (r_t2 - r_t1), " M/s"
                _log_write(0, '(A, T34, F10.4, A)') " Memory throughput: ", real(grid_stats_initial%i_traversed_memory, GRID_SR) / ((1024 * 1024 * 1024) * (r_t2 - r_t1)), " GB/s"
                _log_write(0, '(A, T34, F10.4, A)') " Initialization time:", r_t2 - r_t1, " s"
                _log_write(0, *) ""
                _log_write(0, *) "Execution:"
                _log_write(0, '(A, T34, A)') " Time steps: ", trim(numa%euler%stats%to_string())
                _log_write(0, '(A, T34, A)') " Adaptions: ", trim(adaption_stats_time_steps%to_string())
                _log_write(0, '(A, T34, A)') " Grid: ", trim(grid_stats_time_steps%to_string())
                _log_write(0, '(A, T34, F10.4, A)') " Element throughput: ", 1.0e-6 * real(grid_stats_time_steps%i_traversed_cells, GRID_SR) / (r_t4 - r_t3), " M/s"
                _log_write(0, '(A, T34, F10.4, A)') " Memory throughput: ", real(grid_stats_time_steps%i_traversed_memory, GRID_SR) / ((1024 * 1024 * 1024) * (r_t4 - r_t3)), " GB/s"
                _log_write(0, '(A, T34, F10.4, A)') " Cell update throughput: ", 1.0e-6 * real(numa%euler%stats%i_traversed_cells, GRID_SR) / (r_t4 - r_t3), " M/s"
                _log_write(0, '(A, T34, F10.4, A)') " Flux solver throughput: ", 1.0e-6 * real(numa%euler%stats%i_traversed_edges, GRID_SR) / (r_t4 - r_t3), " M/s"
                _log_write(0, '(A, T34, F10.4, A)') " Execution time:", r_t4 - r_t3, " s"
                _log_write(0, *) ""
                _log_write(0, '(A, T34, F10.4, A)') " Total time:", (r_t2 - r_t1) + (r_t4 - r_t3), " s"
                _log_write(0, *) "---"

                call grid_info%print()
			end if
            !$omp end master
		end subroutine
	END MODULE NUMA
#endif

