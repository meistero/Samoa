! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_PYOP2)
	MODULE pyOP2_init_indices
		use SFC_edge_traversal
		use Samoa

		implicit none

        type num_traversal_data
            integer :: cell_index
            integer :: edge_index
            integer :: vertex_index
        end type

#		define	_GT_NAME						t_pyop2_init_indices_traversal

#       define _GT_EDGES
#		define _GT_NODES

#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_PRE_TRAVERSAL_OP		        pre_traversal_op
#		define _GT_NODE_FIRST_TOUCH_OP		    node_first_touch_vector_op
#		define _GT_INNER_NODE_FIRST_TOUCH_OP    node_first_touch_scalar_op
#		define _GT_EDGE_FIRST_TOUCH_OP		    edge_first_touch_vector_op
#		define _GT_INNER_EDGE_FIRST_TOUCH_OP    edge_first_touch_scalar_op
#		define _GT_CELL_FIRST_TOUCH_OP		    cell_first_touch_op

#		define	_GT_ELEMENT_OP					element_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine pre_traversal_grid_op(traversal, grid)
 			type(t_pyop2_init_indices_traversal), intent(inout)      	:: traversal
 			type(t_grid), intent(inout)							        :: grid

 			traversal%cell_index = 0
 			traversal%edge_index = 0
 			traversal%vertex_index = 0
		end subroutine


		subroutine pre_traversal_op(traversal, section)
 			type(t_pyop2_init_indices_traversal), intent(inout)      	:: traversal
 			type(t_grid_section), intent(inout)							:: section

 			traversal%cell_index = 0
 			traversal%edge_index = 0
 			traversal%vertex_index = 0
		end subroutine

		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, section, element)
 			type(t_pyop2_init_indices_traversal)                :: traversal
 			type(t_grid_section), intent(inout)				    :: section
			type(t_element_base), intent(inout)				    :: element

		end subroutine

		subroutine node_first_touch_vector_op(traversal, section, nodes)
 			type(t_pyop2_init_indices_traversal), intent(inout) :: traversal
 			type(t_grid_section), intent(in)				    :: section
			type(t_node_data), intent(inout)			        :: nodes(:)

            nodes%data_pers%index = traversal%vertex_index

            traversal%vertex_index = traversal%vertex_index + size(nodes)
		end subroutine

        subroutine node_first_touch_scalar_op(traversal, section, node)
 			type(t_pyop2_init_indices_traversal), intent(inout) :: traversal
 			type(t_grid_section), intent(in)				    :: section
			type(t_node_data), intent(inout)			        :: node

            node%data_pers%index = traversal%vertex_index

            traversal%vertex_index = traversal%vertex_index + 1
		end subroutine

		subroutine edge_first_touch_vector_op(traversal, section, edges)
 			type(t_pyop2_init_indices_traversal), intent(inout) :: traversal
 			type(t_grid_section), intent(in)				    :: section
			type(t_edge_data), intent(inout)			        :: edges(:)

            edges%data_pers%index = traversal%edge_index

            traversal%edge_index = traversal%edge_index + size(edges)
		end subroutine

		subroutine edge_first_touch_scalar_op(traversal, section, edge)
 			type(t_pyop2_init_indices_traversal), intent(inout) :: traversal
 			type(t_grid_section), intent(in)				    :: section
			type(t_edge_data), intent(inout)			        :: edge

            edge%data_pers%index = traversal%edge_index

            traversal%edge_index = traversal%edge_index + 1
		end subroutine

		elemental subroutine cell_first_touch_op(traversal, section, cell)
 			type(t_pyop2_init_indices_traversal), intent(inout) :: traversal
 			type(t_grid_section), intent(in)					:: section
			type(t_cell_data_ptr), intent(inout)			    :: cell

            cell%data_pers%index = traversal%cell_index

            traversal%cell_index = traversal%cell_index + 1
		end subroutine
	END MODULE
#endif

