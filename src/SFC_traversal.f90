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
#	elif defined(_FLASH)
		use FLASH
#	elif defined(_NUMA)
		use NUMA
#	elif defined(_GENERIC)
		use Generic
#	endif

	implicit none

	PUBLIC

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
#	    elif defined(_FLASH)
           type(t_flash)          											:: flash
#	    elif defined(_NUMA)
           type(t_numa)                                               		:: numa
#	    elif defined(_GENERIC)
           type(t_generic)                                                  :: generic
#	    endif

		!create, run and destroy scenario

#		if defined(_TESTS)
            !TODO: tests should be able to execute in addition to one of the scenarios!

            !create initial grid
            call init_grid(grid)
			call tests_create(grid, cfg%l_log, cfg%i_asagi_mode)

			!$omp parallel copyin(cfg)
			call tests_run(grid, cfg%i_max_time_steps)
            call grid%destroy()
			!$omp end parallel

			call tests_destroy(grid, cfg%l_log)
#		elif defined (_HEAT_EQ)
            !create initial grid
            call init_grid(grid)
			call heat_eq_create(grid, cfg%l_log, cfg%i_asagi_mode)

			!$omp parallel copyin(cfg)
			call heat_eq_run(grid, cfg%i_max_time_steps, real(cfg%r_max_time, GRID_SR), real(cfg%r_output_time_step, GRID_SR))
            call grid%destroy()
			!$omp end parallel

			call heat_eq_destroy(grid, cfg%l_log)
#		elif defined(_DARCY)
            !create initial grid
            call init_grid(grid)
			call darcy%create(grid, cfg%l_log, cfg%i_asagi_mode)

            !$omp parallel copyin(cfg)
			call darcy%run(grid)
            call grid%destroy()
			!$omp end parallel

			call darcy%destroy(grid, cfg%l_log)
#		elif defined(_SWE)
            !create initial grid
            call init_grid(grid)
			call swe%create(grid, cfg%l_log, cfg%i_asagi_mode)

            !$omp parallel copyin(cfg)
			call swe%run(grid)
            call grid%destroy()
			!$omp end parallel

			call swe%destroy(grid, cfg%l_log)
#		elif defined(_FLASH)
            !create initial grid
            call init_grid(grid)
			call flash%create(grid, cfg%l_log, cfg%i_asagi_mode)

            !$omp parallel copyin(cfg)
			call flash%run(grid, cfg%i_max_time_steps, real(cfg%r_max_time, GRID_SR), real(cfg%r_output_time_step, GRID_SR))
            call grid%destroy()
			!$omp end parallel

			call flash%destroy(grid, cfg%l_log)
#		elif defined(_NUMA)
            !create initial grid
            call init_grid(grid)
			call numa%create(grid, cfg%l_log)

            !$omp parallel copyin(cfg)
			call numa%run(grid, cfg%i_max_time_steps, real(cfg%r_max_time, GRID_SR), real(cfg%r_output_time_step, GRID_SR))
            call grid%destroy()
			!$omp end parallel

			call numa%destroy(grid, cfg%l_log)
#		elif defined(_GENERIC)
            !this scenario is a special case - grid and parallel execution are managed on its own

			call generic%create(cfg%l_log)
			call generic%run()
			call generic%destroy(cfg%l_log)
#		endif
	end subroutine sfc_generic
end MODULE SFC_traversal
