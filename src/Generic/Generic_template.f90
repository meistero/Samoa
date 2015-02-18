! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_GENERIC)
	MODULE generic_traversal
		use SFC_edge_traversal
		use Samoa
        use, intrinsic :: iso_c_binding

        implicit none

        public generic_kernel

        ! Define kernel wrapper interface
        abstract interface
            subroutine generic_kernel(section_index, cell_index, refinement, data)
                use, intrinsic :: iso_c_binding
                integer(kind=c_int), value, intent(in)          :: section_index
                integer(kind=c_long_long), value, intent(in)    :: cell_index
                integer(kind=c_char), intent(inout)             :: refinement
                type(c_ptr), intent(in), value                  :: data
            end subroutine
        end interface

        type num_traversal_data
            procedure(generic_kernel), nopass, pointer    :: kernel => null()
            type(c_ptr)                                 :: data
            logical                                     :: adapt = .false.
        end type

#		define _GT_NAME							t_generic_traversal

#		define _GT_NO_COORDS

#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define _GT_PRE_TRAVERSAL_OP			    pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP			post_traversal_op

#		define _GT_ELEMENT_OP					element_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine pre_traversal_grid_op(traversal, grid)
  			type(t_generic_traversal), intent(inout)	        :: traversal
 			type(t_grid), intent(inout)					    :: grid

 			integer :: i

			traversal%adapt = .false.

			do i = 1, size(traversal%children)
                traversal%children(i)%kernel => traversal%kernel
                traversal%children(i)%data = traversal%data
                traversal%children(i)%adapt = .false.
			end do
		end subroutine

 		subroutine post_traversal_grid_op(traversal, grid)
  			type(t_generic_traversal), intent(inout)	        :: traversal
 			type(t_grid), intent(inout)					    :: grid

            call reduce(traversal%adapt, traversal%children%adapt, MPI_LOR, .true.)
		end subroutine

 		subroutine pre_traversal_op(traversal, section)
 			type(t_generic_traversal), intent(inout)	        :: traversal
  			type(t_grid_section), intent(inout)				:: section

			traversal%adapt = .false.
		end subroutine

 		subroutine post_traversal_op(traversal, section)
 			type(t_generic_traversal), intent(inout)	        :: traversal
  			type(t_grid_section), intent(inout)				:: section
		end subroutine

		!*******************************
		!Geometry operators
		!*******************************

		!element

		subroutine element_op(traversal, section, element)
 			type(t_generic_traversal), intent(inout)	    :: traversal
 			type(t_grid_section), intent(inout)			    :: section
			type(t_element_base), intent(inout)		        :: element

			integer (kind = c_long_long)        :: cell_index
			integer (kind = c_char)             :: refinement

			cell_index = element%cell%data_pers%index
            refinement = 0

			call traversal%kernel(section%index - 1, cell_index, refinement, traversal%data)

			element%cell%geometry%refinement = refinement
            traversal%adapt = traversal%adapt .or. (refinement .ne. 0)
		end subroutine
	END MODULE
#endif
