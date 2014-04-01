! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_TESTS)
	MODULE Tests_flops_data
		use SFC_node_traversal

		integer (kind = GRID_SI)				:: i_traversed_flops
		real (kind = GRID_SR)					:: r_dummy

		!*******************************
		!Geometry operators
		!*******************************

	END MODULE

	MODULE Tests_flops_30_edges_and_nodes
		use Tests_flops_data

#		define _GT_POST_TRAVERSAL_OP			post_traversal_op
#		define _GT_ELEMENT_OP					element_op

#		define _GT_NAME							Tests_flops_30_traversal

#		define _GT_EDGES						.true.
#		define _GT_EDGES_TEMP					.true.
#		define _GT_NODES						.true.
#		define _GT_NODES_TEMP					.true.

#		include "SFC_generic_traversal_ringbuffer.f90"

#		undef _GT_EDGES
#		undef _GT_EDGES_TEMP
#		undef _GT_NODES
#		undef _GT_NODES_TEMP

#		undef _GT_NAME

#		undef _GT_ELEMENT_OP
#		undef _GT_POST_TRAVERSAL_OP

		subroutine post_traversal_op(traversal, section)
			implicit none
			type(triangle_tree), dimension(:), intent(inout), target	:: triangles				! initial triangle array
			integer (kind = GRID_SI), intent(in)						:: num_coarse_triangles

			i_traversed_flops = i_traversed_flops + 30
		end subroutine

		subroutine element_op(element)
			implicit none
			type(t_element_base), intent(inout), target		:: element

		end subroutine
	END MODULE
#endif
