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
		use SFC_node_traversal

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

        type(t_grid), allocatable, save               :: grid

		private
		public t_pyop2

		contains

		!> Creates all required runtime objects for the scenario
		subroutine pyop2_create(pyop2, l_log)
            class(t_pyop2)                                              :: pyop2
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
		subroutine pyop2_destroy(pyop2, l_log)
            class(t_pyop2)                  :: pyop2
			logical (kind = GRID_SL)		:: l_log

			if (l_log) then
				_log_close_file()
			endif
		end subroutine

		subroutine print_indices_kernel(cell_index, edge_indices, vertex_indices, coords, refinement)
            use, intrinsic :: iso_c_binding
            integer(kind=c_int), value, intent(in) :: cell_index
            integer(kind=c_int), intent(in) :: edge_indices(3)
            integer(kind=c_int), intent(in) :: vertex_indices(3)
            real(kind=c_double), intent(in) :: coords(6)
            integer(kind=c_char), intent(inout) :: refinement

            _log_write(1, '("cell index: ", I0 , " edge indices: ", 3(I0, X) , " vertex indices: ", 3(I0, X))') cell_index, edge_indices, vertex_indices
            _log_write(1, '("coords: ", 3("( ", 2(F0.3, X), ") "))') coords
        end subroutine

		!> Sets the initial values of the scenario and runs the time steps
		subroutine pyop2_run(pyop2)
            class(t_pyop2)                                              :: pyop2

            !test the wrapper (call must come from a non-parallel-region)
            !give a kernel as argument and let the wrapper generate a grid
            call samoa_run_f90_kernel(print_indices_kernel)
		end subroutine

		subroutine samoa_run_c_kernel(c_kernel) bind(c)
            use, intrinsic :: iso_c_binding
            type(c_funptr), value, intent(in)              :: c_kernel

            procedure(pyop2_kernel), pointer          :: f90_kernel

            ! Convert C to Fortran procedure pointer.
            call c_f_procpointer(c_kernel, f90_kernel)

            call samoa_run_f90_kernel(f90_kernel)
        end subroutine

		subroutine samoa_run_f90_kernel(f90_kernel)
            procedure(pyop2_kernel), pointer, intent(in)        :: f90_kernel

            type(t_pyop2_adaptive_traversal), save              :: adaptive_traversal
            type(t_pyop2_traversal), save                       :: kernel_traversal
            type(t_pyop2_init_indices_traversal), save          :: init_indices

            if (.not. allocated(grid)) then
                call init_MPI()

                !init element transformation data
                call init_transform_data()

                !init grid

                allocate(grid)
                grid%i_min_depth = 1
                grid%i_max_depth = 14
                grid%i_sections_per_thread = 1

                call omp_set_num_threads(1)

                call init_grid(grid)

                !set entity indices
                call init_indices%traverse(grid)
            end if

            !$omp parallel default(shared)
                if (kernel_traversal%adapt) then
                    !$omp barrier
                    call adaptive_traversal%traverse(grid)

                    !set entity indices
                    call init_indices%traverse(grid)
                end if

                !attach kernel and run traversal
                !$omp single
                    kernel_traversal%kernel => f90_kernel
                !$omp end single

                call kernel_traversal%traverse(grid)

                !$omp barrier

                !$omp single
                    kernel_traversal%kernel => null()
                !$omp end single
            !$omp end parallel
        end subroutine
	END MODULE PyOP2
#endif

