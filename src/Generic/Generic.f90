! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_GENERIC)
	MODULE Generic
		use Generic_data_types
		use Generic_init_indices
		use Generic_traversal
		use Generic_adaptive_traversal
		use SFC_edge_traversal
		use SFC_node_traversal

        use, intrinsic :: iso_c_binding

		implicit none

		type t_generic
            type(t_generic_init_indices_traversal)    :: init_indices
            type(t_generic_traversal)                 :: traversal
            type(t_generic_adaptive_traversal)        :: adaptive_traversal

            contains

            procedure, pass :: create => generic_create
            procedure, pass :: run => generic_run
            procedure, pass :: destroy => generic_destroy
        end type

        type t_c_section
            integer (kind = c_long_long)    :: i_cells, i_edges, i_nodes
            integer (kind = c_long_long)    :: i_bcells, i_bedges, i_bnodes
            type(c_ptr)                     :: cells_to_edges, cells_to_nodes, edges_to_nodes
            type(c_ptr)                     :: bcells_to_cells, bedges_to_edges, bnodes_to_nodes
            type(c_ptr)                     :: coords
        end type

        type, extends(t_grid) :: t_c_grid
            type(t_c_section), allocatable  :: c_sections(:)
            logical :: in_use = .false.
        end type

        type(t_c_grid), allocatable, save   :: grids(:)


		private
		public t_generic, samoa_run_kernel, samoa_get_grid

		contains

		!> Creates all required runtime objects for the scenario
		subroutine generic_create(generic, l_log)
            class(t_generic)                                              :: generic
			logical , intent(in)						:: l_log

			!local variables
			character (len = 64)										:: s_log_name, s_date, s_time

			!open log file
			call date_and_time(s_date, s_time)
			write (s_log_name, '(A, "output/generic_", A, "_", A, ".log")') trim(s_date), trim(s_time)

			if (l_log) then
				_log_open_file(s_log_name)
			endif
		end subroutine


		!> Destroys all required runtime objects for the scenario
		subroutine generic_destroy(generic, l_log)
            class(t_generic)                  :: generic
			logical 		:: l_log

			if (l_log) then
				_log_close_file()
			endif
		end subroutine

		subroutine print_indices_kernel(section_index, cell_index, refinement, data)
            use, intrinsic :: iso_c_binding
            integer(kind=c_int), value, intent(in)          :: section_index
            integer(kind=c_long_long), value, intent(in)    :: cell_index
            integer(kind=c_char), intent(inout)             :: refinement
            type(c_ptr), value, intent(in)                  :: data

            _log_write(1, '("section index: ", I0 , " cell index: ", I0)') section_index, cell_index
        end subroutine

        !***********
        !C interface
        !***********

        function samoa_create_grid(i_sections, i_depth) result(handle) bind(c)
            integer (kind = c_int), value, intent(in)           :: i_sections
            integer (kind = c_char), value, intent(in)          :: i_depth
            integer (kind = c_int)                              :: handle

            type(t_c_grid), allocatable                         :: grids_tmp(:)
            type(t_generic_init_indices_traversal)                :: init_indices

            do handle = 1, size(grids)
                if (.not. grids(handle)%in_use) then
                    exit
                end if
            end do

            if (.not. allocated(grids)) then
                handle = 1
                allocate(grids(handle))
            else
                if (handle > size(grids)) then
                    allocate(grids_tmp(handle))

                    grids_tmp(1 : size(grids)) = grids
                    call move_alloc(grids_tmp, grids)
                end if
            end if

            grids(handle)%in_use = .true.

            !init MPI and element transformation data

            call omp_set_num_threads(1);
            call init_mpi()
            call init_transform_data()

            !init grid
            grids(handle)%i_sections_per_thread = i_sections
            call init_grid(grids(handle), i_depth)

            !set entity indices
            call init_indices%traverse(grids(handle))
        end function

        subroutine samoa_destroy_grid(handle) bind(c)
            integer (kind = c_int), value, intent(in)       :: handle
            assert_ge(size(grids), handle)

            assert(grids(handle)%in_use)

            call grids(handle)%destroy()

            if (allocated(grids(handle)%c_sections)) then
                deallocate(grids(handle)%c_sections)
            end if

            grids(handle)%in_use = .false.

            call finalize_mpi()
        end subroutine

		subroutine samoa_run_c_kernel(handle, data, c_kernel) bind(c,name="samoa_run_kernel")
            integer (kind = c_int), value, intent(in)   :: handle
            type(c_ptr), intent(in), value              :: data
            type(c_funptr), intent(in), value           :: c_kernel

            procedure(generic_kernel), pointer            :: f90_kernel

            assert_ge(size(grids), handle)

            ! Convert C to Fortran procedure pointer.
            call c_f_procpointer(c_kernel, f90_kernel)

            call samoa_run_kernel(handle, data, f90_kernel)
        end subroutine

        subroutine samoa_get_grid(handle, i_sections, c_sections) bind(c)
            integer (kind = c_int), value, intent(in)   :: handle
            integer (kind = c_int), intent(out)         :: i_sections
            type(c_ptr), intent(out)                    :: c_sections

            type(t_grid_section), pointer               :: section
            integer (kind= GRID_SI)                     :: i_first_section, i_last_section, i_section

            assert_ge(size(grids), handle)

            i_sections = grids(handle)%sections%get_size()

            if (allocated(grids(handle)%c_sections)) then
                if (size(grids(handle)%c_sections) .ne. i_sections) then
                    deallocate(grids(handle)%c_sections)
                    allocate(grids(handle)%c_sections(i_sections))
                end if
            else
                allocate(grids(handle)%c_sections(i_sections))
            end if

            !$omp parallel
            call grids(handle)%get_local_sections(i_first_section, i_last_section)

            do i_section = i_first_section, i_last_section
                section => grids(handle)%sections%elements_alloc(i_section)

                grids(handle)%c_sections(i_section)%i_cells = section%i_cells
                grids(handle)%c_sections(i_section)%i_edges = section%i_edges
                grids(handle)%c_sections(i_section)%i_nodes = section%i_nodes
                grids(handle)%c_sections(i_section)%cells_to_edges = c_loc(section%cells_to_edges_map)
                grids(handle)%c_sections(i_section)%cells_to_nodes = c_loc(section%cells_to_nodes_map)
                grids(handle)%c_sections(i_section)%edges_to_nodes = c_loc(section%edges_to_nodes_map)

                grids(handle)%c_sections(i_section)%i_bcells = section%i_bcells
                grids(handle)%c_sections(i_section)%i_bedges = section%i_bedges
                grids(handle)%c_sections(i_section)%i_bnodes = section%i_bnodes
                grids(handle)%c_sections(i_section)%bcells_to_cells = c_loc(section%bcells_to_cells_map)
                grids(handle)%c_sections(i_section)%bedges_to_edges = c_loc(section%bedges_to_edges_map)
                grids(handle)%c_sections(i_section)%bnodes_to_nodes = c_loc(section%bnodes_to_nodes_map)

                grids(handle)%c_sections(i_section)%coords = c_loc(section%coords)
            end do
            !$omp end parallel

            c_sections = c_loc(grids(handle)%c_sections)
        end subroutine


        !*****************
        !Fortran interface
        !*****************

		!> Runs a test kernel
		subroutine generic_run(generic)
            class(t_generic)                                              :: generic

            integer :: grid_handle

            !test the wrapper (call must come from a non-parallel-region)
            !give a kernel as argument and let the wrapper generate a grid
            grid_handle = samoa_create_grid(1, 3)

            call samoa_run_kernel(grid_handle, c_ptr(null()), print_indices_kernel)

            call samoa_destroy_grid(grid_handle)
		end subroutine

		subroutine samoa_run_kernel(handle, data, kernel)
            integer (kind = c_int), intent(in)                  :: handle
            type(c_ptr), value, intent(in)                      :: data
            procedure(generic_kernel), pointer, intent(in)        :: kernel

            type(t_generic_adaptive_traversal)                    :: adaptive_traversal
            type(t_generic_traversal)                             :: kernel_traversal
            type(t_generic_init_indices_traversal)                :: init_indices

            assert_ge(size(grids), handle)
            assert(grids(handle)%in_use)

            kernel_traversal%kernel => kernel
            kernel_traversal%data = data

            !$omp parallel default(shared)
                call kernel_traversal%traverse(grids(handle))
                !$omp barrier

                if (kernel_traversal%adapt) then
                    !$omp barrier
                    call adaptive_traversal%traverse(grids(handle)%t_grid)

                    !set entity indices
                    call init_indices%traverse(grids(handle))
                end if

                !$omp barrier
            !$omp end parallel
        end subroutine
	END MODULE Generic
#endif

