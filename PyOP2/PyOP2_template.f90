! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_PYOP2)
	MODULE mod_%(name)
		use SFC_edge_traversal

        implicit none

        type num_traversal_data
            !additional data definitions
        end type

#		define _GT_NAME							%(name)

#		define _GT_EDGES
#		define _GT_NODES

#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define _GT_PRE_TRAVERSAL_OP			    pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP			post_traversal_op

#		define _GT_ELEMENT_OP					element_op

#		define _GT_INNER_NODE_FIRST_TOUCH_OP	inner_node_first_touch_op
#		define _GT_INNER_NODE_LAST_TOUCH_OP		inner_node_last_touch_op
#		define _GT_NODE_REDUCE_OP		        node_reduce_op
#		define _GT_NODE_MERGE_OP		        node_merge_op

#		define _GT_INNER_EDGE_FIRST_TOUCH_OP	inner_edge_first_touch_op
#		define _GT_INNER_EDGE_LAST_TOUCH_OP		inner_edge_last_touch_op
#		define _GT_EDGE_REDUCE_OP		        edge_reduce_op
#		define _GT_EDGE_MERGE_OP		        edge_merge_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine pre_traversal_grid_op(traversal, grid)
  			type(%(name)), intent(inout)	    :: traversal
 			type(t_grid), intent(inout)					    :: grid


		end subroutine

 		subroutine post_traversal_grid_op(traversal, grid)
  			type(%(name)), intent(inout)	    :: traversal
 			type(t_grid), intent(inout)					    :: grid


		end subroutine

 		subroutine pre_traversal_op(traversal, section)
 			type(t_darcy_jacobi_solver), intent(inout)	    :: traversal
  			type(t_grid_section), intent(inout)				:: section


		end subroutine

 		subroutine post_traversal_op(traversal, section)
 			type(t_darcy_jacobi_solver), intent(inout)	    :: traversal
  			type(t_grid_section), intent(inout)				:: section


		end subroutine

		!*******************************
		!Geometry operators
		!*******************************

		!element

		subroutine element_op(traversal, section, element)
 			type(%(name)), intent(inout)	                :: traversal
 			type(t_grid_section), intent(inout)			    :: section
			type(t_element_base), intent(inout), target		:: element


		end subroutine

		! first touches

		subroutine inner_node_first_touch_op(traversal, section, node)
 			type(t_darcy_jacobi_solver), intent(in)		    :: traversal
 			type(t_grid_section), intent(in)			    :: section
			type(t_node_data), intent(inout)			    :: node


		end subroutine

		subroutine inner_edge_first_touch_op(traversal, section, node)
 			type(%(name)), intent(in)		                :: traversal
 			type(t_grid_section), intent(in)			    :: section
			type(t_node_data), intent(inout)			    :: node


		end subroutine

		!last touches

		subroutine inner_node_last_touch_op(traversal, section, node)
 			type(t_darcy_jacobi_solver), intent(in)		    :: traversal
 			type(t_grid_section), intent(in)			    :: section
			type(t_node_data), intent(inout)			    :: node


		end subroutine

		subroutine inner_edge_last_touch_op(traversal, section, edge)
 			type(%(name)), intent(in)		                :: traversal
 			type(t_grid_section), intent(in)			    :: section
			type(t_edge_data), intent(inout)			    :: edge


		end subroutine

		subroutine node_reduce_op(traversal, section, node)
 			type(%(name)), intent(inout)	                :: traversal
 			type(t_grid_section), intent(in)			    :: section
			type(t_node_data), intent(in)			        :: node


		end subroutine

		subroutine node_merge_op(local_node, neighbor_node)
 			type(t_node_data), intent(inout)			    :: local_node
			type(t_node_data), intent(in)				    :: neighbor_node


        end subroutine
	END MODULE
#endif
