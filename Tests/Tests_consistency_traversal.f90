! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_TESTS)
	MODULE Tests_consistency_ringbuffer
		use SFC_node_traversal

		integer (kind = GRID_SI)				:: i_node_ref_cnt, i_edge_ref_cnt
		integer (kind = GRID_SI)				:: i_bnd_node_ref_cnt, i_bnd_edge_ref_cnt

#		define _GT_NAME							tests_consistency_ringbuffer_traversal
#		define _GT_EDGES
#		define _GT_EDGES_TEMP
#		define _GT_NODES
#		define _GT_NODES_TEMP

#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP			post_traversal_op

#		define _GT_CELL_FIRST_TOUCH_OP			cell_first_touch_op
#		define _GT_EDGE_FIRST_TOUCH_OP			edge_first_touch_op
#		define _GT_NODE_FIRST_TOUCH_OP			node_first_touch_op
#		define _GT_BND_EDGE_FIRST_TOUCH_OP		bnd_edge_first_touch_op
#		define _GT_BND_NODE_FIRST_TOUCH_OP		bnd_node_first_touch_op

#		define _GT_CELL_LAST_TOUCH_OP			cell_last_touch_op
#		define _GT_EDGE_LAST_TOUCH_OP			edge_last_touch_op
#		define _GT_NODE_LAST_TOUCH_OP			node_last_touch_op
#		define _GT_BND_EDGE_LAST_TOUCH_OP		bnd_edge_last_touch_op
#		define _GT_BND_NODE_LAST_TOUCH_OP		bnd_node_last_touch_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine pre_traversal_op(traversal, section)
			implicit none
			type(triangle_tree), dimension(:), intent(inout), target	:: triangles				! initial triangle array
			integer (kind = GRID_SI), intent(in)						:: num_coarse_triangles

			_log_write(0,  *) "Tests: running traversal consistency test"

			i_node_ref_cnt = 0
			i_edge_ref_cnt = 0
			i_bnd_node_ref_cnt = 0
			i_bnd_edge_ref_cnt = 0
		end subroutine

		subroutine post_traversal_op(traversal, section)
			implicit none
			type(triangle_tree), dimension(:), intent(inout), target	:: triangles				! initial triangle array
			integer (kind = GRID_SI), intent(in)						:: num_coarse_triangles

			assert_eq(i_node_ref_cnt, 0)
			assert_eq(i_edge_ref_cnt, 0)
			assert_eq(i_bnd_node_ref_cnt, 0)
			assert_eq(i_bnd_edge_ref_cnt, 0)
		end subroutine

		!*******************************
		!Geometry operators
		!*******************************

		!first touches

		subroutine cell_first_touch_op(cell)
			implicit none
			type(t_cell_data_ptr), intent(inout)				:: cell
		end subroutine

		subroutine edge_first_touch_op(edge)
			implicit none
			type(t_edge_data), intent(inout)				:: edge

			edge%data_pers%i_ref_cnt = edge%data_pers%i_ref_cnt + 1
			i_edge_ref_cnt = i_edge_ref_cnt + 1
			assert_eq(edge%data_pers%i_ref_cnt, 1)
		end subroutine

		subroutine node_first_touch_op(node)
			implicit none
			type(t_node_data), intent(inout)				:: node

			node%data_pers%i_ref_cnt = node%data_pers%i_ref_cnt + 1
			i_node_ref_cnt = i_node_ref_cnt + 1
			assert_eq(node%data_pers%i_ref_cnt, 1)
		end subroutine

		subroutine bnd_edge_first_touch_op(edge)
			implicit none
			type(t_edge_data), intent(inout)				:: edge

			edge%data_pers%i_ref_cnt = edge%data_pers%i_ref_cnt + 1
			i_bnd_edge_ref_cnt = i_bnd_edge_ref_cnt + 1
			assert_eq(edge%data_pers%i_ref_cnt, 1)
		end subroutine

		subroutine bnd_node_first_touch_op(node)
			implicit none
			type(t_node_data), intent(inout)				:: node

			node%data_pers%i_ref_cnt = node%data_pers%i_ref_cnt + 1
			i_bnd_node_ref_cnt = i_bnd_node_ref_cnt + 1
			assert_eq(node%data_pers%i_ref_cnt, 1)
		end subroutine

		!last touches

		subroutine cell_last_touch_op(cell)
			implicit none
			type(t_cell_data_ptr), intent(inout)				:: cell
		end subroutine

		subroutine edge_last_touch_op(edge)
			implicit none
			type(t_edge_data), intent(inout)				:: edge

			edge%data_pers%i_ref_cnt = edge%data_pers%i_ref_cnt - 1
   			i_edge_ref_cnt = i_edge_ref_cnt - 1
			assert_eq(edge%data_pers%i_ref_cnt, 0)
		end subroutine

		subroutine node_last_touch_op(node)
			implicit none
			type(t_node_data), intent(inout)				:: node

			node%data_pers%i_ref_cnt = node%data_pers%i_ref_cnt - 1
   			i_node_ref_cnt = i_node_ref_cnt - 1
			assert_eq(node%data_pers%i_ref_cnt, 0)
		end subroutine

		subroutine bnd_edge_last_touch_op(edge)
			implicit none
			type(t_edge_data), intent(inout)				:: edge

			edge%data_pers%i_ref_cnt = edge%data_pers%i_ref_cnt - 1
   			i_bnd_edge_ref_cnt = i_bnd_edge_ref_cnt - 1
			assert_eq(edge%data_pers%i_ref_cnt, 0)
		end subroutine

		subroutine bnd_node_last_touch_op(node)
			implicit none
			type(t_node_data), intent(inout)				:: node

			node%data_pers%i_ref_cnt = node%data_pers%i_ref_cnt - 1
   			i_bnd_node_ref_cnt = i_bnd_node_ref_cnt - 1
			assert_eq(node%data_pers%i_ref_cnt, 0)
		end subroutine
	END MODULE

	MODULE Tests_consistency
		use SFC_node_traversal

		integer (kind = GRID_SI)				:: i_node_ref_cnt, i_edge_ref_cnt
		integer (kind = GRID_SI)				:: i_bnd_node_ref_cnt, i_bnd_edge_ref_cnt

#		define _GT_NAME							tests_consistency_traversal
#		define _GT_CELL							.true.
#		define _GT_CELL_TEMP					.true.
#		define _GT_EDGES						.true.
#		define _GT_EDGES_TEMP					.true.
#		define _GT_NODES						.true.
#		define _GT_NODES_TEMP					.true.

#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP			post_traversal_op

#		define _GT_CELL_FIRST_TOUCH_OP			cell_first_touch_op
#		define _GT_EDGE_FIRST_TOUCH_OP			edge_first_touch_op
#		define _GT_NODE_FIRST_TOUCH_OP			node_first_touch_op
#		define _GT_BND_EDGE_FIRST_TOUCH_OP		bnd_edge_first_touch_op
#		define _GT_BND_NODE_FIRST_TOUCH_OP		bnd_node_first_touch_op

#		define _GT_CELL_LAST_TOUCH_OP			cell_last_touch_op
#		define _GT_EDGE_LAST_TOUCH_OP			edge_last_touch_op
#		define _GT_NODE_LAST_TOUCH_OP			node_last_touch_op
#		define _GT_BND_EDGE_LAST_TOUCH_OP		bnd_edge_last_touch_op
#		define _GT_BND_NODE_LAST_TOUCH_OP		bnd_node_last_touch_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine pre_traversal_op(traversal, section)
			implicit none
			type(triangle_tree), dimension(:), intent(inout), target	:: triangles				! initial triangle array
			integer (kind = GRID_SI), intent(in)						:: num_coarse_triangles

			_log_write(0,  *) "Tests: running ringbuffer traversal consistency test.."

			i_node_ref_cnt = 0
			i_edge_ref_cnt = 0
			i_bnd_node_ref_cnt = 0
			i_bnd_edge_ref_cnt = 0
		end subroutine

		subroutine post_traversal_op(traversal, section)
			implicit none
			type(triangle_tree), dimension(:), intent(inout), target	:: triangles				! initial triangle array
			integer (kind = GRID_SI), intent(in)						:: num_coarse_triangles

			assert_eq(i_node_ref_cnt, 0)
			assert_eq(i_edge_ref_cnt, 0)
			assert_eq(i_bnd_node_ref_cnt, 0)
			assert_eq(i_bnd_edge_ref_cnt, 0)
		end subroutine

		!*******************************
		!Geometry operators
		!*******************************

		!first touches

		subroutine cell_first_touch_op(cell)
			implicit none
			type(t_cell_data_ptr), intent(inout)				:: cell
		end subroutine

		subroutine edge_first_touch_op(edge)
			implicit none
			type(t_edge_data), intent(inout)				:: edge

			edge%data_pers%i_ref_cnt = edge%data_pers%i_ref_cnt + 1
			i_edge_ref_cnt = i_edge_ref_cnt + 1
			assert_eq(edge%data_pers%i_ref_cnt, 1)
		end subroutine

		subroutine node_first_touch_op(node)
			implicit none
			type(t_node_data), intent(inout)				:: node

			node%data_pers%i_ref_cnt = node%data_pers%i_ref_cnt + 1
			i_node_ref_cnt = i_node_ref_cnt + 1
			assert_eq(node%data_pers%i_ref_cnt, 1)
		end subroutine

		subroutine bnd_edge_first_touch_op(edge)
			implicit none
			type(t_edge_data), intent(inout)				:: edge

			edge%data_pers%i_ref_cnt = edge%data_pers%i_ref_cnt + 1
			i_bnd_edge_ref_cnt = i_bnd_edge_ref_cnt + 1
			assert_eq(edge%data_pers%i_ref_cnt, 1)
		end subroutine

		subroutine bnd_node_first_touch_op(node)
			implicit none
			type(t_node_data), intent(inout)				:: node

			node%data_pers%i_ref_cnt = node%data_pers%i_ref_cnt + 1
			i_bnd_node_ref_cnt = i_bnd_node_ref_cnt + 1
			assert_eq(node%data_pers%i_ref_cnt, 1)
		end subroutine

		!last touches

		subroutine cell_last_touch_op(cell)
			implicit none
			type(t_cell_data_ptr), intent(inout)				:: cell
		end subroutine

		subroutine edge_last_touch_op(edge)
			implicit none
			type(t_edge_data), intent(inout)				:: edge

			edge%data_pers%i_ref_cnt = edge%data_pers%i_ref_cnt - 1
   			i_edge_ref_cnt = i_edge_ref_cnt - 1
			assert_eq(edge%data_pers%i_ref_cnt, 0)
		end subroutine

		subroutine node_last_touch_op(node)
			implicit none
			type(t_node_data), intent(inout)				:: node

			node%data_pers%i_ref_cnt = node%data_pers%i_ref_cnt - 1
   			i_node_ref_cnt = i_node_ref_cnt - 1
			assert_eq(node%data_pers%i_ref_cnt, 0)
		end subroutine

		subroutine bnd_edge_last_touch_op(edge)
			implicit none
			type(t_edge_data), intent(inout)				:: edge

			edge%data_pers%i_ref_cnt = edge%data_pers%i_ref_cnt - 1
   			i_bnd_edge_ref_cnt = i_bnd_edge_ref_cnt - 1
			assert_eq(edge%data_pers%i_ref_cnt, 0)
		end subroutine

		subroutine bnd_node_last_touch_op(node)
			implicit none
			type(t_node_data), intent(inout)				:: node

			node%data_pers%i_ref_cnt = node%data_pers%i_ref_cnt - 1
   			i_bnd_node_ref_cnt = i_bnd_node_ref_cnt - 1
			assert_eq(node%data_pers%i_ref_cnt, 0)
		end subroutine
	END MODULE
#endif
