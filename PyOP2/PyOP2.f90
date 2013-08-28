! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_PYOP2)
	MODULE PyOP2
		use PyOP2_data_types
		use PyOP2_init_indices
		use PyOP2_traversal
		use PyOP2_adaptive_traversal
		use SFC_edge_traversal

		implicit none

		type t_pyop2
            type(t_pyop2_init_indices_traversal)    :: init_indices
            type(t_pyop2_traversal)                 :: traversal
            type(t_pyop2_adaptive_traversal)        :: adaptive_traversal

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

            !test the wrapper (call must come from a non-parallel-region)
            !give a kernel as argument and let the wrapper generate a grid
            call run_f90_kernel(print_indices_kernel)
		end subroutine

		!> Sets the initial values of the scenario and runs the time steps
		subroutine pyop2_run(pyop2, grid)
            class(t_pyop2)                                              :: pyop2
 			type(t_grid), intent(inout)									:: grid

            call pyop2%init_indices%traverse(grid)

            !test the traversal
            !bind a kernel and call traversal directly

            pyop2%traversal%kernel => print_indices_kernel
            call pyop2%traversal%traverse(grid)
            call pyop2%adaptive_traversal%traverse(grid)
		end subroutine

		subroutine print_indices_kernel(cell_index, edge_indices, vertex_indices, coords)
            use, intrinsic :: iso_c_binding
            integer(kind=c_int), intent(in) :: cell_index
            integer(kind=c_int), intent(in) :: edge_indices(3)
            integer(kind=c_int), intent(in) :: vertex_indices(3)
            real(kind=c_double), intent(in) :: coords(6)

            _log_write(1, '("cell index: ", I0 , " edge indices: ", 3(I0, X) , " vertex indices: ", 3(I0, X))') cell_index, edge_indices, vertex_indices
            _log_write(1, '("coords: ", 3("( ", 2(F0.3, X), ") "))') coords
        end subroutine
	END MODULE PyOP2
#endif

