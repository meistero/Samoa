! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_PYOP2)
	MODULE pyop2_traversal
		use SFC_edge_traversal
		use SFC_node_traversal
		use Samoa
		use pyOP2_init_indices
        use, intrinsic :: iso_c_binding

        implicit none

        public pyop2_kernel, run_f90_kernel, run_c_kernel

        ! Define kernel wrapper interface
        abstract interface
            subroutine pyop2_kernel(cell_index, edge_indices, vertex_indices, coords)
                use, intrinsic :: iso_c_binding
                integer(kind=c_int), intent(in) :: cell_index
                integer(kind=c_int), intent(in) :: edge_indices(3)
                integer(kind=c_int), intent(in) :: vertex_indices(3)
                real(kind=c_double), intent(in) :: coords(6)
            end subroutine
        end interface

        type num_traversal_data
            procedure(pyop2_kernel), nopass, pointer :: kernel => null()
        end type

#		define _GT_NAME							t_pyop2_traversal

#		define _GT_EDGES
#		define _GT_NODES

#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define _GT_PRE_TRAVERSAL_OP			    pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP			post_traversal_op

#		define _GT_ELEMENT_OP					element_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine pre_traversal_grid_op(traversal, grid)
  			type(t_pyop2_traversal), intent(inout)	        :: traversal
 			type(t_grid), intent(inout)					    :: grid

			traversal%children%kernel => traversal%kernel
		end subroutine

 		subroutine post_traversal_grid_op(traversal, grid)
  			type(t_pyop2_traversal), intent(inout)	        :: traversal
 			type(t_grid), intent(inout)					    :: grid

		end subroutine

 		subroutine pre_traversal_op(traversal, section)
 			type(t_pyop2_traversal), intent(inout)	        :: traversal
  			type(t_grid_section), intent(inout)				:: section

		end subroutine

 		subroutine post_traversal_op(traversal, section)
 			type(t_pyop2_traversal), intent(inout)	        :: traversal
  			type(t_grid_section), intent(inout)				:: section

		end subroutine

		!*******************************
		!Geometry operators
		!*******************************

		!element

		subroutine element_op(traversal, section, element)
 			type(t_pyop2_traversal), intent(inout)	        :: traversal
 			type(t_grid_section), intent(inout)			    :: section
			type(t_element_base), intent(inout), target		:: element

            real (kind = c_double),  parameter  :: base_coords(2, 3) = [[1, 0, 0], [0, 0, 1]]
			real (kind = c_double)              :: coords(6)
			integer (kind = c_int)              :: vertex_indices(3), edge_indices(3), cell_index, i

			cell_index = element%cell%data_pers%index

            do i = 1, 3
                edge_indices(i) = element%edges(i)%ptr%data_pers%index
                vertex_indices(i) = element%nodes(i)%ptr%data_pers%index

                coords(2 * i - 1 : 2 * i) = samoa_barycentric_to_world_point(element%transform_data, base_coords(:, i))
            end do

			call traversal%kernel(cell_index, edge_indices, vertex_indices, coords)
		end subroutine

		subroutine run_c_kernel(c_kernel) bind(c)
            use, intrinsic :: iso_c_binding
            type(c_funptr), value, intent(in)              :: c_kernel

            procedure(pyop2_kernel), pointer          :: f90_kernel

            ! Convert C to Fortran procedure pointer.
            call c_f_procpointer(c_kernel, f90_kernel)

            call run_f90_kernel(f90_kernel)
        end subroutine

		subroutine run_f90_kernel(f90_kernel)
            procedure(pyop2_kernel), pointer, intent(in)        :: f90_kernel

            type(t_pyop2_traversal)                             :: pyop2_wrapper
            type(t_pyop2_init_indices_traversal)                :: init_indices
            type(t_grid), allocatable, save                     :: grid

            if (.not. allocated(grid)) then
                call init_MPI()

                !init element transformation data
                call init_transform_data()

                !init grid

                allocate(grid)
                grid%i_min_depth = 1
                grid%i_max_depth = 14
                grid%i_sections_per_thread = 8

                call init_grid(grid)

                !set entity indices
                call init_indices%traverse(grid)
            end if

            !attach kernel and run traversal
            pyop2_wrapper%kernel => f90_kernel

            !$omp parallel
            call pyop2_wrapper%traverse(grid)
            !$omp end parallel

            pyop2_wrapper%kernel => null()
        end subroutine
	END MODULE
#endif
