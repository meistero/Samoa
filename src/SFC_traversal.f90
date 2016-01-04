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

        call set_signals()

		!create, run and destroy scenario

#		if defined(_TESTS)
            !TODO: tests should be able to execute in addition to one of the scenarios!

            !create initial grid
            call init_grid(grid, cfg%i_start_depth)
			call tests_create(grid, cfg%l_log, cfg%i_asagi_mode)

			!$omp parallel copyin(cfg)
			call tests_run(grid)
            call grid%destroy()
			!$omp end parallel

			call tests_destroy(grid, cfg%l_log)
#		elif defined (_HEAT_EQ)
            !create initial grid
            call init_grid(grid, cfg%i_start_depth)
			call heat_eq_create(grid, cfg%l_log, cfg%i_asagi_mode)

			!$omp parallel copyin(cfg)
			call heat_eq_run(grid, cfg%i_start_depth)
            call grid%destroy()
			!$omp end parallel

			call heat_eq_destroy(grid, cfg%l_log)
#		elif defined(_DARCY)
            !create initial grid
            call init_grid(grid, cfg%i_start_depth)
			call darcy%create(grid, cfg%l_log, cfg%i_asagi_mode)

            !$omp parallel copyin(cfg)
			call darcy%run(grid)
            call grid%destroy()
			!$omp end parallel

			call darcy%destroy(grid, cfg%l_log)
#		elif defined(_SWE)
            !create initial grid
            call init_grid(grid, cfg%i_start_depth)
			call swe%create(grid, cfg%l_log, cfg%i_asagi_mode)

            !$omp parallel copyin(cfg)
			call swe%run(grid)
            call grid%destroy()
			!$omp end parallel

			call swe%destroy(grid, cfg%l_log)
#		elif defined(_FLASH)
            !create initial grid
            call init_grid(grid, cfg%i_start_depth)
			call flash%create(grid, cfg%l_log, cfg%i_asagi_mode)

            !$omp parallel copyin(cfg)
			call flash%run(grid)
            call grid%destroy()
			!$omp end parallel

			call flash%destroy(grid, cfg%l_log)
#		elif defined(_NUMA)
            !create initial grid
            call init_grid(grid, cfg%i_start_depth)
			call numa%create(grid, cfg%l_log)

            !$omp parallel copyin(cfg)
			call numa%run(grid)
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

    subroutine set_signals()
#       if !defined(__GFORTRAN__)
            use ifport
#       endif

        integer :: i_status

#       if defined(__GFORTRAN__)
            !i_status = signal(2, signal_interrupt)
#       else
            !i_status = signal(2, signal_interrupt, -1)
#       endif
    end subroutine

    subroutine unset_signal(i_signal)
#       if !defined(__GFORTRAN__)
            use ifport
#       endif

        integer(4), intent(in) :: i_signal

        integer :: i_status

#       if defined(__GFORTRAN__)
            !i_status = signal(i_signal, 0)
#       else
            !i_status = signal(i_signal, signal_interrupt, 0)
#       endif
    end subroutine

	function signal_interrupt() result(rcode)
        integer(4) :: rcode

        print '("*****************************************************")'
        print '("*** Interrupt signal caught, initiating soft exit ***")'
        print '("*****************************************************")'

        !initiate a soft exit by setting all possible exit conditions to return immediately
        cfg%i_min_depth = 1
        cfg%i_max_depth = 1
        cfg%i_max_time_steps = 0
        cfg%r_max_time = 0.0_GRID_SR

#       if defined (_DARCY)
            cfg%i_max_iterations = 0
#       endif

        !if ctrl-c is pressed again, the program should terminate immediately
        call unset_signal(2)

        rcode = 0
    end function
end MODULE SFC_traversal
