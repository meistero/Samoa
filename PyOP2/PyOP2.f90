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

        type t_c_section
            integer (kind = c_long_long)    :: i_cells, i_edges, i_nodes
            type(c_ptr)                     :: cells_to_edges, cells_to_nodes, edges_to_nodes
            type(c_ptr)                     :: coords
        end type

        type(t_grid), allocatable, save             :: grid
        type(t_c_section), allocatable, save        :: c_sections(:)

		private
		public t_pyop2, samoa_run_kernel, samoa_get_grid

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

		subroutine print_indices_kernel(section_index, cell_index, refinement)
            use, intrinsic :: iso_c_binding
            integer(kind=c_int), value, intent(in)          :: section_index
            integer(kind=c_long_long), value, intent(in)    :: cell_index
            integer(kind=c_char), intent(inout)             :: refinement

            _log_write(1, '("section index: ", I0 , "cell index: ", I0)') section_index, cell_index
        end subroutine

        !***********
        !C interface
        !***********

        subroutine samoa_create(i_depth, i_sections) bind(c)
            integer (kind = c_char), value, intent(in)          :: i_depth
            integer (kind = c_int), value, intent(in)           :: i_sections

            type(t_pyop2_init_indices_traversal), save          :: init_indices

            assert(.not. allocated(grid))

            allocate(grid)

            !init MPI and element transformation data
            call init_MPI()
            call init_transform_data()

            !init grid
            grid%i_sections_per_thread = i_sections
            call init_grid(grid, i_depth)

            !set entity indices
            call init_indices%traverse(grid)
        end subroutine

        subroutine samoa_destroy() bind(c)
            assert(allocated(grid))

            call grid%destroy()
            deallocate(grid)

            call finalize_mpi()
        end subroutine

		subroutine samoa_run_c_kernel(c_kernel) bind(c,name="samoa_run_kernel")
            use, intrinsic :: iso_c_binding
            type(c_funptr), intent(in), value   :: c_kernel

            procedure(pyop2_kernel), pointer    :: f90_kernel

            ! Convert C to Fortran procedure pointer.
            call c_f_procpointer(c_kernel, f90_kernel)

            call samoa_run_kernel(f90_kernel)
        end subroutine

        subroutine samoa_get_grid(i_sections, c_grid) bind(c)
            integer (kind = c_int), intent(out) :: i_sections
            type(c_ptr), intent(out)            :: c_grid

            type(t_grid_section), pointer       :: section
            integer (kind= GRID_SI)             :: i_first_section, i_last_section, i_section

            if (.not. allocated(grid)) then
                !init the grid with default settings
                call omp_set_num_threads(1)
                call samoa_create(0, 1)
            end if

            i_sections = size(grid%sections%elements_alloc)

            if (allocated(c_sections)) then
                if (size(c_sections) .ne. i_sections) then
                    deallocate(c_sections)
                    allocate(c_sections(i_sections))
                end if
            else
                allocate(c_sections(i_sections))
            end if

            !$omp parallel

            call grid%get_local_sections(i_first_section, i_last_section)

            do i_section = i_first_section, i_last_section
                section => grid%sections%elements_alloc(i_section)

                c_sections(i_section)%i_cells = section%i_cells
                c_sections(i_section)%i_edges = section%i_edges
                c_sections(i_section)%i_nodes = section%i_nodes
                c_sections(i_section)%cells_to_edges = c_loc(section%cells_to_edges_map)
                c_sections(i_section)%cells_to_nodes = c_loc(section%cells_to_nodes_map)
                c_sections(i_section)%edges_to_nodes = c_loc(section%edges_to_nodes_map)
                c_sections(i_section)%coords = c_loc(section%coords)
            end do

            c_grid = c_loc(c_sections)

            !$omp end parallel
        end subroutine


        !*****************
        !Fortran interface
        !*****************

		!> Runs a test kernel
		subroutine pyop2_run(pyop2)
            class(t_pyop2)                                              :: pyop2

            !test the wrapper (call must come from a non-parallel-region)
            !give a kernel as argument and let the wrapper generate a grid
            call samoa_run_kernel(print_indices_kernel)
		end subroutine

		subroutine samoa_run_kernel(kernel)
            procedure(pyop2_kernel), pointer, intent(in)        :: kernel

            type(t_pyop2_adaptive_traversal), save              :: adaptive_traversal
            type(t_pyop2_traversal), save                       :: kernel_traversal
            type(t_pyop2_init_indices_traversal), save          :: init_indices

            if (.not. allocated(grid)) then
                !init the grid with default settings
                call omp_set_num_threads(1)
                call samoa_create(0, 1)
            end if

            !$omp parallel default(shared)
                !attach kernel and run traversal
                !$omp single
                    kernel_traversal%kernel => kernel
                !$omp end single

                call kernel_traversal%traverse(grid)

                !$omp barrier

                !$omp single
                    kernel_traversal%kernel => null()
                !$omp end single

                if (kernel_traversal%adapt) then
                    !$omp barrier
                    call adaptive_traversal%traverse(grid)

                    !set entity indices
                    call init_indices%traverse(grid)
                end if
            !$omp end parallel
        end subroutine
	END MODULE PyOP2
#endif

