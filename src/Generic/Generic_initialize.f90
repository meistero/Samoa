! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_GENERIC)
	MODULE generic_init_indices
		use SFC_edge_traversal
		use Samoa

		implicit none

        type num_traversal_data
            integer (kind=c_long_long)          :: cell_index, edge_index, node_index
            integer (kind=c_long_long)          :: bcell_index, bedge_index, bnode_index
        end type

#		define	_GT_NAME						t_generic_init_indices_traversal

#       define _GT_EDGES
#		define _GT_NODES

#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_PRE_TRAVERSAL_OP		        pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP		    post_traversal_op
#		define _GT_NODE_FIRST_TOUCH_OP		    node_first_touch_op
#		define _GT_INNER_NODE_FIRST_TOUCH_OP    inner_node_first_touch_op
#		define _GT_EDGE_FIRST_TOUCH_OP		    edge_first_touch_op
#		define _GT_INNER_EDGE_FIRST_TOUCH_OP    inner_edge_first_touch_op
#		define _GT_CELL_FIRST_TOUCH_OP		    cell_first_touch_op

#		define	_GT_ELEMENT_OP					element_op
#		define	_GT_INNER_ELEMENT_OP			inner_element_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine pre_traversal_grid_op(traversal, grid)
 			type(t_generic_init_indices_traversal), intent(inout)      	:: traversal
 			type(t_grid), intent(inout)							        :: grid

		end subroutine


		subroutine pre_traversal_op(traversal, section)
 			type(t_generic_init_indices_traversal), intent(inout)      	:: traversal
 			type(t_grid_section), intent(inout)							:: section

 			traversal%cell_index = 0
 			traversal%edge_index = 0
 			traversal%node_index = 0

 			traversal%bcell_index = 0
 			traversal%bedge_index = 0
 			traversal%bnode_index = 0

 			section%i_bnodes = section%boundary_nodes(RED)%get_size() + section%boundary_nodes(GREEN)%get_size()
 			section%i_bedges = section%boundary_edges(RED)%get_size() + section%boundary_edges(GREEN)%get_size()
 			section%i_bcells = section%i_bedges !this is an upper bound, there must be less boundary elements than boundary edges

 			section%i_nodes = section%nodes_in%get_size() + section%i_bnodes
 			section%i_edges = section%crossed_edges_in%get_size() + section%color_edges_in%get_size() + section%i_bedges
 			section%i_cells = section%cells%get_size()

            if (.not. allocated(section%cells_to_edges_map)) then
                allocate(section%cells_to_edges_map(3, 0 : section%i_cells - 1))
                allocate(section%cells_to_nodes_map(3, 0 : section%i_cells - 1))
                allocate(section%edges_to_nodes_map(2, 0 : section%i_edges - 1))
                allocate(section%bcells_to_cells_map(0 : section%i_bcells - 1))
                allocate(section%bedges_to_edges_map(0 : section%i_bedges - 1))
                allocate(section%bnodes_to_nodes_map(0 : section%i_bnodes - 1))
                allocate(section%coords(2, 0 : section%i_nodes))
            end if
		end subroutine

		subroutine post_traversal_op(traversal, section)
 			type(t_generic_init_indices_traversal), intent(inout)      	:: traversal
 			type(t_grid_section), intent(inout)							:: section

 			section%i_bcells = traversal%bcell_index
 			section%i_bedges = traversal%bedge_index
 			section%i_bnodes = traversal%bnode_index
		end subroutine

		!******************
		!Geometry operators
		!******************

		subroutine inner_element_op(traversal, section, element)
 			type(t_generic_init_indices_traversal)                :: traversal
 			type(t_grid_section), intent(inout)				    :: section
			type(t_element_base), intent(inout)				    :: element

            real (kind = c_double),  parameter  :: base_coords(2, 3) = [[1, 0, 0], [0, 0, 1]]
			integer                             :: node_indices(3), edge_indices(3), cell_index
            integer (kind = BYTE)                  :: edge_types(3)

			integer (kind = BYTE)                  :: i_previous_edge, i_color_edge, i_next_edge, i

			cell_index = element%cell%data_pers%index

            do i = 1, 3
                edge_indices(i) = element%edges(i)%ptr%data_pers%index
                node_indices(i) = element%nodes(i)%ptr%data_pers%index

                section%coords(:, node_indices(i)) = samoa_barycentric_to_world_point(element%transform_data, base_coords(:, i))
            end do

			call element%cell%geometry%get_edge_indices(i_previous_edge, i_color_edge, i_next_edge)
			call element%cell%geometry%get_edge_types(edge_types(i_previous_edge), edge_types(i_color_edge), edge_types(i_next_edge))

            call init_maps(section %i_cells, section%i_edges, &
                section%cells_to_edges_map, section%cells_to_nodes_map, section%edges_to_nodes_map, &
                cell_index, edge_indices, node_indices, edge_types)
		end subroutine

		subroutine element_op(traversal, section, element)
 			type(t_generic_init_indices_traversal)                :: traversal
 			type(t_grid_section), intent(inout)				    :: section
			type(t_element_base), intent(inout)				    :: element

            call inner_element_op(traversal, section, element)

			section%bcells_to_cells_map(traversal%bcell_index) = element%cell%data_pers%index
			traversal%bcell_index = traversal%bcell_index + 1
		end subroutine


		subroutine init_maps(i_cells, i_edges, cells_to_edges_map, cells_to_nodes_map, edges_to_nodes_map, i_c, i_e, i_v, edge_types)
 			integer (kind=c_long_long), intent(in)              :: i_cells, i_edges
 			integer (kind=c_long_long), intent(inout)           :: cells_to_edges_map(3, 0 : i_cells - 1)
 			integer (kind=c_long_long), intent(inout)           :: cells_to_nodes_map(3, 0 : i_cells - 1)
 			integer (kind=c_long_long), intent(inout)           :: edges_to_nodes_map(2, 0 : i_edges - 1)
 			integer, intent(in)                                 :: i_c, i_e(3), i_v(3)
 			integer (kind = BYTE), intent(in)                      :: edge_types(3)

            integer :: i

            do i = 1, 3
                cells_to_nodes_map(i, i_c) = i_v(i)
            end do

            do i = 1, 3
                cells_to_edges_map(i, i_c) = i_e(i)
            end do

            do i = 1, 3
                select case (edge_types(i))
                    case (NEW, NEW_BND)
                        edges_to_nodes_map(1, i_e(i)) = i_v(1 + mod(i, 3))
                        edges_to_nodes_map(2, i_e(i)) = i_v(1 + mod(i + 1, 3))
                    case (OLD_BND)
                        edges_to_nodes_map(1, i_e(i)) = i_v(1 + mod(i + 1, 3))
                        edges_to_nodes_map(2, i_e(i)) = i_v(1 + mod(i, 3))
                    case (OLD)
                        assert_eq(edges_to_nodes_map(1, i_e(i)), i_v(1 + mod(i + 1, 3)))
                        assert_eq(edges_to_nodes_map(2, i_e(i)), i_v(1 + mod(i, 3)))
                end select
            end do
        end subroutine

		subroutine node_first_touch_op(traversal, section, nodes)
 			type(t_generic_init_indices_traversal), intent(inout) :: traversal
 			type(t_grid_section), intent(inout)				    :: section
			type(t_node_data), intent(inout)			        :: nodes(:)

            integer :: i

            forall (i = 1 : size(nodes))
                nodes(i)%data_pers%index = traversal%node_index + i - 1
                section%bnodes_to_nodes_map(traversal%bnode_index + i - 1) = traversal%node_index + i - 1
            end forall

            traversal%node_index = traversal%node_index + size(nodes)
            traversal%bnode_index = traversal%bnode_index + size(nodes)
		end subroutine

        subroutine inner_node_first_touch_op(traversal, section, node)
 			type(t_generic_init_indices_traversal), intent(inout) :: traversal
 			type(t_grid_section), intent(in)				    :: section
			type(t_node_data), intent(inout)			        :: node

            node%data_pers%index = traversal%node_index

            traversal%node_index = traversal%node_index + 1
		end subroutine

		subroutine edge_first_touch_op(traversal, section, edges)
 			type(t_generic_init_indices_traversal), intent(inout) :: traversal
 			type(t_grid_section), intent(inout)				    :: section
			type(t_edge_data), intent(inout)			        :: edges(:)

            integer :: i

            forall (i = 1 : size(edges))
                edges(i)%data_pers%index = traversal%edge_index + i - 1
                section%bedges_to_edges_map(traversal%bedge_index + i - 1) = traversal%edge_index + i - 1
            end forall

            traversal%edge_index = traversal%edge_index + size(edges)
            traversal%bedge_index = traversal%bedge_index + size(edges)
		end subroutine

		subroutine inner_edge_first_touch_op(traversal, section, edge)
 			type(t_generic_init_indices_traversal), intent(inout) :: traversal
 			type(t_grid_section), intent(in)				    :: section
			type(t_edge_data), intent(inout)			        :: edge

            edge%data_pers%index = traversal%edge_index

            traversal%edge_index = traversal%edge_index + 1
		end subroutine

		elemental subroutine cell_first_touch_op(traversal, section, cell)
 			type(t_generic_init_indices_traversal), intent(inout) :: traversal
 			type(t_grid_section), intent(in)					:: section
			type(t_cell_data_ptr), intent(inout)			    :: cell

            cell%data_pers%index = traversal%cell_index

            traversal%cell_index = traversal%cell_index + 1
		end subroutine
	END MODULE
#endif

