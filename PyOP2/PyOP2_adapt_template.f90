! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_PYOP2)
	MODULE mod_%(name)
		use SFC_edge_traversal
		use Conformity

        implicit none

        type num_traversal_data
            !additional data definitions
        end type

#		define	_GT_NAME					    %(name)

#		define _GT_EDGES
#		define _GT_NODES

#		define _GT_TRANSFER_OP				    transfer_op
#		define _GT_REFINE_OP				    refine_op
#		define _GT_COARSEN_OP				    coarsen_op

#		define _GT_INNER_NODE_FIRST_TOUCH_OP	inner_node_first_touch_op
#		define _GT_INNER_NODE_LAST_TOUCH_OP		inner_node_last_touch_op
#		define _GT_NODE_REDUCE_OP		        node_reduce_op
#		define _GT_NODE_MERGE_OP		        node_merge_op

#		define _GT_INNER_EDGE_FIRST_TOUCH_OP	inner_edge_first_touch_op
#		define _GT_INNER_EDGE_LAST_TOUCH_OP		inner_edge_last_touch_op
#		define _GT_EDGE_REDUCE_OP		        edge_reduce_op
#		define _GT_EDGE_MERGE_OP		        edge_merge_op

#		include "SFC_generic_adaptive_traversal.f90"

		!******************
		!Geometry operators
		!******************

		subroutine transfer_op(traversal, section, src_element, dest_element)
 			type(%(name)), intent(inout)							:: traversal
 			type(t_grid_section), intent(inout)							            :: section
			type(t_traversal_element), intent(inout)									:: src_element
			type(t_traversal_element), intent(inout)									:: dest_element


		end subroutine

		subroutine refine_op(traversal, section, src_element, dest_element, refinement_path)
 			type(%(name)), intent(inout)							:: traversal
 			type(t_grid_section), intent(inout)										:: section
			type(t_traversal_element), intent(inout)								:: src_element
			type(t_traversal_element), intent(inout)								:: dest_element
			integer, dimension(:), intent(in)										:: refinement_path


		end subroutine

		subroutine coarsen_op(traversal, section, src_element, dest_element, refinement_path)
  			type(%(name)), intent(inout)							:: traversal
			type(t_grid_section), intent(inout)													:: section
			type(t_traversal_element), intent(inout)									:: src_element
			type(t_traversal_element), intent(inout)									:: dest_element
			integer, dimension(:), intent(in)											:: refinement_path


		end subroutine

		! first touches

		subroutine inner_node_first_touch_op(traversal, section, node)
 			type(t_darcy_jacobi_solver), intent(in)		    :: traversal
 			type(t_grid_section), intent(in)			    :: section
			type(t_node_data), intent(inout)			    :: node


		end subroutine

		subroutine inner_edge_first_touch_op(traversal, section, node)
 			type(%(name)), intent(in)		    :: traversal
 			type(t_grid_section), intent(in)			    :: section
			type(t_node_data), intent(inout)			    :: node


		end subroutine

		!last touches

		subroutine inner_node_last_touch_op(traversal, section, node)
 			type(t_darcy_jacobi_solver), intent(in)		    :: traversal
 			type(t_grid_section), intent(in)			    :: section
			type(t_node_data), intent(inout)			    :: node


		end subroutine

		subroutine inner_edge_last_touch_op(traversal, section, node)
 			type(%(name)), intent(in)		    :: traversal
 			type(t_grid_section), intent(in)			    :: section
			type(t_node_data), intent(inout)			    :: node


		end subroutine

		subroutine node_reduce_op(traversal, section, node)
 			type(%(name)), intent(inout)	    :: traversal
 			type(t_grid_section), intent(in)			    :: section
			type(t_node_data), intent(in)			        :: node


		end subroutine

		subroutine node_merge_op(local_node, neighbor_node)
 			type(t_node_data), intent(inout)			    :: local_node
			type(t_node_data), intent(in)				    :: neighbor_node


        end subroutine
	END MODULE
#endif
