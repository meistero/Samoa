! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_PYOP2)
	MODULE PyOP2
		use PyOP2_data_types
		use PyOP2_init_indices
		use SFC_edge_traversal

		implicit none

		type t_pyop2
            type(t_pyop2_init_indices_traversal)    :: init_indices

            contains

            procedure , pass :: create => pyop2_create
            procedure , pass :: run => pyop2_run
            procedure , pass :: destroy => pyop2_destroy
        end type

		private
		public t_pyop2

		contains

		!> Creates all required runtime objects for the scenario
		subroutine pyop2_create(pyop2, grid, l_log)
            class(t_pyop2)                                              :: pyop2
 			type(t_grid), intent(inout)									:: grid
			logical (kind = GRID_SL), intent(in)						:: l_log

			!local variables
			character (len = 64)										:: s_log_name, s_date, s_time

			!open log file
			call date_and_time(s_date, s_time)
			write (s_log_name, '(A, "output/pyop2_", A, "_", A, ".log")') trim(s_date), trim(s_time)

			if (l_log) then
				_log_open_file(s_log_name)
			endif
		end subroutine


		!> Destroys all required runtime objects for the scenario
		subroutine pyop2_destroy(pyop2, grid, l_log)
            class(t_pyop2)                  :: pyop2
 			type(t_grid), intent(inout)     :: grid
			logical (kind = GRID_SL)		:: l_log

			if (l_log) then
				_log_close_file()
			endif
		end subroutine

		!> Sets the initial values of the scenario and runs the time steps
		subroutine pyop2_run(pyop2, grid)
            class(t_pyop2)                                              :: pyop2
 			type(t_grid), intent(inout)									:: grid

            call pyop2%init_indices%traverse(grid)
		end subroutine
	END MODULE PyOP2
#endif

