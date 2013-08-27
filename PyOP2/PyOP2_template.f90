! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_PYOP2)
	MODULE pyop2_traversal
		use SFC_edge_traversal
		use Samoa
        use, intrinsic :: iso_c_binding

        implicit none

        ! Define kernel wrapper interface
        abstract interface
            subroutine kernel_wrapper_op(cell_index, edge_indices, vertex_indices, coords)
                use, intrinsic :: iso_c_binding
                integer(kind=c_int), intent(in) :: cell_index
                integer(kind=c_int), intent(in) :: edge_indices(3)
                integer(kind=c_int), intent(in) :: vertex_indices(3)
                real(kind=c_double), intent(in) :: coords(2, 3)
            end subroutine
        end interface

        type num_traversal_data
            procedure(kernel_wrapper_op), nopass, pointer :: kernel_wrapper
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

			traversal%children%kernel_wrapper => traversal%kernel_wrapper
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
 			type(t_pyop2_traversal), intent(inout)	    :: traversal
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
			real (kind = c_double)              :: coords(2, 3)
			integer (kind = c_int)              :: vertex_indices(3), edge_indices(3), cell_index, i

			cell_index = element%cell%data_pers%index

            do i = 1, 3
                edge_indices(i) = element%edges(i)%ptr%data_pers%index
                vertex_indices(i) = element%nodes(i)%ptr%data_pers%index

                coords(:, i) = samoa_barycentric_to_world_point(element%transform_data, base_coords(:, i))
            end do

			call traversal%kernel_wrapper(cell_index, edge_indices, vertex_indices, coords)
		end subroutine

		subroutine run_kernel(c_kernel_wrapper) bind(c)
            type(c_funptr), intent(in)      :: c_kernel_wrapper

            type(t_pyop2_traversal)         :: traversal
            type(t_grid)                    :: grid !TODO: this should be a global variable that is reused during traversals

            ! Convert C to Fortran procedure pointer.
            call c_f_procpointer(c_kernel_wrapper, traversal%kernel_wrapper)

            call traversal%traverse(grid)
        end subroutine
	END MODULE
#endif
