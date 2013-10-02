! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


!*****************************************************************
! MODULE SFC_traversal: traversal for multigrid
!*****************************************************************
#include "Compilation_control.f90"

MODULE SFC_traversal
	use SFC_data_types
    use SFC_edge_traversal
	use SFC_node_traversal

#	if defined(_HEAT_EQ)
		use Heat_Eq
#	elif defined(_DARCY)
		use Darcy
#	elif defined(_TESTS)
		use Tests
#	elif defined(_SWE)
		use SWE
#	elif defined(_NUMA)
		use NUMA
#	elif defined(_GENERIC)
		use Generic
#	endif

	implicit none

	PUBLIC

	integer (kind = GRID_SI)			:: i_sections_per_thread							!< number of sections per thread
	integer (kind = GRID_SI)			:: i_max_time_steps									!< number of simulation time steps
	real (kind = GRID_SR)				:: r_max_time, r_output_time_step					!< maximum simulation time, outpout time step
	logical (kind = GRID_SL)			:: l_log                                            !< if true, a log file is used
	integer (kind = 1)			        :: i_min_depth, i_max_depth			                !< number of simulation time steps
	integer			        			:: i_asagi_mode			                		!< ASAGI hints

	contains

	subroutine sfc_generic()

		! local variables
		type(t_grid)														:: grid

#       if defined(_HEAT_EQ)

#	    elif defined(_TESTS)

#	    elif defined(_DARCY)
            type(t_darcy)                                                   :: darcy
#	    elif defined(_SWE)
           type(t_swe)          											:: swe
#	    elif defined(_NUMA)
           type(t_numa)                                               		:: numa
#	    elif defined(_GENERIC)
           type(t_generic)                                               		:: generic
#	    endif

        if (rank_MPI == 0) then
             _log_write(0, ' (" sam(oa)²: Space filling curves and Adaptive Meshes for Oceanic and Other Applications")')
             _log_write(0, "")

#    	    if defined(_TESTS)
                _log_write(0, '(" Scenario: Tests")')
#    		elif defined (_HEAT_EQ)
                _log_write(0, '(" Scenario: Heat Equation")')
#    		elif defined(_DARCY)
                _log_write(0, '(" Scenario: Darcy")')
#    		elif defined(_SWE)
                _log_write(0, '(" Scenario: SWE")')
#    		elif defined(_NUMA)
                _log_write(0, '(" Scenario: NUMA")')
#    		elif defined(_GENERIC)
                _log_write(0, '(" Scenario: GENERIC")')
#    		endif

#    		if defined(_OMP)
                _log_write(0, '(" OpenMP: Yes, threads: ", I0, ", procs: ", I0)') omp_get_max_threads(), omp_get_num_procs()
#    		else
                _log_write(0, '(" OpenMP: No")')
#    		endif

#    		if defined(_MPI)
                _log_write(0, '(" MPI: Yes, ranks: ", I0)') size_MPI
#    		else
                _log_write(0, '(" MPI: No")')
#    		endif

#    		if defined(_ASAGI)
                _log_write(0, '(" ASAGI: Yes, mode: ", I0)') i_asagi_mode
#    		else
                _log_write(0, '(" ASAGI: No")')
#    		endif

#    		if defined(_DEBUG_LEVEL)
                _log_write(0, '(" Debug Level: ", I0)') _DEBUG_LEVEL
#    		endif

#    		if defined(_ASSERT)
                _log_write(0, '(" Assertions: Yes")')
#    		else
                _log_write(0, '(" Assertions: No")')
#    		endif

            _log_write(0, '(" Sections per thread: ", I0)') i_sections_per_thread
            _log_write(0, '(" Adaptivity: min depth: ", I0, ", max depth: ", I0)') i_min_depth, i_max_depth
            _log_write(0, '(" Simulation: max time steps: ", I0, ", max time: ", ES9.2, ", output step: ", ES9.2)'), i_max_time_steps, r_max_time, r_output_time_step
            _log_write(0, "")
        end if

		!convert ASAGI hints to ASAGI format

#		if defined(_ASAGI)
			select case(i_asagi_mode)
				case (0)
					i_asagi_mode = GRID_NO_HINT
				case (1)
					i_asagi_mode = ieor(GRID_NOMPI, GRID_PASSTHROUGH)
				case (2)
					i_asagi_mode = GRID_NOMPI
				case (3)
					i_asagi_mode = ieor(GRID_NOMPI, SMALL_CACHE)
				case (4)
					i_asagi_mode = GRID_LARGE_GRID
			end select
#		endif

		!create, run and destroy scenario

#		if defined(_TESTS)
            grid%i_sections_per_thread = i_sections_per_thread

            !TODO: tests should be able to execute in addition to one of the scenarios!

            !create initial grid
            call init_grid(grid)
			call tests_create(grid, l_log, i_asagi_mode)

			!$omp parallel
			call tests_run(grid, i_max_time_steps)
			!$omp end parallel

			call tests_destroy(grid, l_log)
            call grid%destroy()
#		elif defined (_HEAT_EQ)
            grid%i_min_depth = i_min_depth
            grid%i_max_depth = i_max_depth
            grid%i_sections_per_thread = i_sections_per_thread

            !create initial grid
            call init_grid(grid)
			call heat_eq_create(grid, l_log, i_asagi_mode)

			!$omp parallel
			call heat_eq_run(grid, i_max_time_steps, r_max_time, r_output_time_step)
			!$omp end parallel

			call heat_eq_destroy(grid, l_log)
            call grid%destroy()
#		elif defined(_DARCY)
            grid%i_min_depth = i_min_depth
            grid%i_max_depth = i_max_depth
            grid%i_sections_per_thread = i_sections_per_thread

            !create initial grid
            call init_grid(grid)
			call darcy%create(grid, l_log, i_asagi_mode)

            !$omp parallel
			call darcy%run(grid, i_max_time_steps, r_max_time, r_output_time_step)
			!$omp end parallel

			call darcy%destroy(grid, l_log)
            call grid%destroy()
#		elif defined(_SWE)
            grid%i_min_depth = i_min_depth
            grid%i_max_depth = i_max_depth
            grid%i_sections_per_thread = i_sections_per_thread

            !create initial grid
            call init_grid(grid)
			call swe%create(grid, l_log, i_asagi_mode)

            !$omp parallel
			call swe%run(grid, i_max_time_steps, r_max_time, r_output_time_step)
			!$omp end parallel

			call swe%destroy(grid, l_log)
            call grid%destroy()
#		elif defined(_NUMA)
            grid%i_min_depth = i_min_depth
            grid%i_max_depth = i_max_depth
            grid%i_sections_per_thread = i_sections_per_thread

            !create initial grid
            call init_grid(grid)
			call numa%create(grid, l_log)

            !$omp parallel
			call numa%run(grid, i_max_time_steps, r_max_time, r_output_time_step)
			!$omp end parallel

			call numa%destroy(grid, l_log)
            call grid%destroy()
#		elif defined(_GENERIC)
            !this scenario is a special case - grid and parallel execution are managed on its own

			call generic%create(l_log)
			call generic%run()
			call generic%destroy(l_log)
#		endif
	end subroutine sfc_generic
end MODULE SFC_traversal
