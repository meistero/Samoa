! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_TESTS)
	MODULE Tests_initialize
		use SFC_node_traversal

#		define	_GT_NAME						tests_init_traversal
#		define	_GT_EDGES						.true.
#		define	_GT_EDGES_TEMP					.false.
#		define	_GT_NODES						.true.
#		define	_GT_NODES_TEMP					.false.

#		define _GT_INNER_EDGE_FIRST_TOUCH_OP			edge_first_touch_op
#		define _GT_INNER_NODE_FIRST_TOUCH_OP			node_first_touch_op
#		define _GT_EDGE_FIRST_TOUCH_OP		edge_first_touch_op
#		define _GT_NODE_FIRST_TOUCH_OP		node_first_touch_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		!undefine macros to avoid compiler warnings
#		undef	_GT_NAME
#		undef	_GT_CELL
#		undef	_GT_CELL_TEMP
#		undef	_GT_EDGES
#		undef	_GT_EDGES_TEMP
#		undef	_GT_NODES
#		undef	_GT_NODES_TEMP

#		undef	_GT_INNER_EDGE_FIRST_TOUCH_OP
#		undef	_GT_INNER_NODE_FIRST_TOUCH_OP
#		undef	_GT_EDGE_FIRST_TOUCH_OP
#		undef	_GT_NODE_FIRST_TOUCH_OP

		!******************
		!Geometry operators
		!******************

		!first touches

		subroutine edge_first_touch_op(edge)
			implicit none
			type(t_edge_data), intent(inout)				:: edge

			edge%data_pers%i_ref_cnt = 0
		end subroutine

		subroutine node_first_touch_op(node)
			implicit none
			type(t_node_data), intent(inout)				:: node

			node%data_pers%i_ref_cnt = 0
		end subroutine
	END MODULE
#endif
