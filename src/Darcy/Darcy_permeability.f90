! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_permeability
		use SFC_edge_traversal
		use Darcy_initialize_saturation

		use Samoa_darcy

        type num_traversal_data
        end type

		type(darcy_gv_saturation)				:: gv_saturation
		type(darcy_gv_p)						:: gv_p
		type(darcy_gv_rhs)						:: gv_rhs
		type(darcy_gv_is_dirichlet_boundary)    :: gv_is_dirichlet
		type(darcy_gv_volume)                   :: gv_inflow

#		define _GT_NAME							t_darcy_permeability_traversal

#		if (_DARCY_FLOW_EDGE_SIZE > 0)
#			define _GT_EDGES
#		endif

#		define _GT_NODES

#		define _GT_ELEMENT_OP					element_op

#		define _GT_NODE_FIRST_TOUCH_OP		    node_first_touch_op

#		define _GT_NODE_MERGE_OP		        node_merge_op
#		define _GT_NODE_LAST_TOUCH_OP		    node_last_touch_op

#		define _GT_EDGE_MPI_TYPE

#		include "SFC_generic_traversal_ringbuffer.f90"


        subroutine create_edge_mpi_type(mpi_edge_type)
            integer, intent(out)            :: mpi_edge_type

#           if defined(_MPI)
                type(t_edge_data)                       :: edge
                integer                                 :: blocklengths(2), types(2), disps(2), type_size, i_error
                integer (kind = MPI_ADDRESS_KIND)       :: lb, ub

                blocklengths(1) = 1
                blocklengths(2) = 1

                disps(1) = 0
                disps(2) = sizeof(edge)

                types(1) = MPI_LB
                types(2) = MPI_UB

                call MPI_Type_struct(2, blocklengths, disps, types, mpi_edge_type, i_error); assert_eq(i_error, 0)
                call MPI_Type_commit(mpi_edge_type, i_error); assert_eq(i_error, 0)

                call MPI_Type_size(mpi_edge_type, type_size, i_error); assert_eq(i_error, 0)
                call MPI_Type_get_extent(mpi_edge_type, lb, ub, i_error); assert_eq(i_error, 0)

                assert_eq(0, lb)
                assert_eq(0, type_size)
                assert_eq(sizeof(edge), ub)
#           endif
        end subroutine

		!*******************************
		!Geometry operators
		!*******************************

		elemental subroutine node_first_touch_op(traversal, section, node)
 			type(t_darcy_permeability_traversal), intent(in)                     :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			:: node

			call flow_pre_dof_op(node%position(1), node%position(2), node%data_pers%p, node%data_pers%saturation, node%data_pers%rhs, node%data_temp%is_dirichlet_boundary, node%data_temp%volume)
		end subroutine

		subroutine element_op(traversal, section, element)
 			type(t_darcy_permeability_traversal), intent(inout)		:: traversal
 			type(t_grid_section), intent(inout)						:: section
			type(t_element_base), intent(inout)				        :: element

			real (kind = GRID_SR)   :: saturation(_DARCY_LAYERS + 1, 3)
			real (kind = GRID_SR)   :: p(_DARCY_LAYERS + 1, 3)
			real (kind = GRID_SR)   :: rhs(_DARCY_LAYERS + 1, 3)

			call gv_saturation%read_from_element(element, saturation)
			call gv_p%read_from_element(element, p)

			!call element operator
#           if (_DARCY_LAYERS > 0)
                call initialize_rhs(element, saturation, p, rhs, element%cell%data_pers%base_permeability)
#           else
                call initialize_rhs(element, saturation(1,:), p(1,:), rhs(1,:), element%cell%data_pers%base_permeability)
#           endif

			call gv_rhs%add_to_element(element, rhs)
		end subroutine

		elemental subroutine node_merge_op(local_node, neighbor_node)
			type(t_node_data), intent(inout)		:: local_node
 			type(t_node_data), intent(in)		    :: neighbor_node

            if (any(neighbor_node%data_temp%is_dirichlet_boundary)) then
                local_node%data_pers%p = neighbor_node%data_pers%p
            end if

			local_node%data_pers%saturation = max(local_node%data_pers%saturation, neighbor_node%data_pers%saturation)
			local_node%data_temp%is_dirichlet_boundary = local_node%data_temp%is_dirichlet_boundary .or. neighbor_node%data_temp%is_dirichlet_boundary
			local_node%data_pers%rhs = local_node%data_pers%rhs + neighbor_node%data_pers%rhs
			local_node%data_temp%volume = local_node%data_temp%volume + neighbor_node%data_temp%volume
		end subroutine

		elemental subroutine node_last_touch_op(traversal, section, node)
 			type(t_darcy_permeability_traversal), intent(in)            :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			                :: node

			real (kind = GRID_SR) :: total_inflow

            total_inflow = tiny(1.0_SR) + sum(node%data_temp%volume)

			call flow_post_dof_op(node%data_pers%saturation, node%data_pers%rhs, node%data_temp%is_dirichlet_boundary, node%data_temp%volume, total_inflow)
		end subroutine

		elemental subroutine flow_pre_dof_op(pos_x, pos_y, p, saturation, rhs, is_dirichlet, inflow)
			real (kind = GRID_SR), intent(in)					:: pos_x, pos_y
			real (kind = GRID_SR), intent(inout)				:: p
			real (kind = GRID_SR), intent(inout)				:: saturation
			real (kind = GRID_SR), intent(out)					:: rhs
			logical, intent(out)			                    :: is_dirichlet
			real (kind = GRID_SR), intent(out)					:: inflow

			rhs = 0.0_SR
			is_dirichlet = .false.
			inflow = 0.0_SR

#           if !defined(_ASAGI)
                if (pos_x == 0.0_SR) then
                    saturation = 1.0_SR
                end if

                if (pos_x == 1.0_SR) then
                    is_dirichlet = .true.
                    p = cfg%r_p_prod
                end if
#           endif
		end subroutine

		elemental subroutine flow_post_dof_op(saturation, rhs, is_dirichlet, inflow, total_inflow)
			real (kind = GRID_SR), intent(inout)				:: saturation
			real (kind = GRID_SR), intent(inout)				:: rhs
			real (kind = GRID_SR), intent(in)				    :: inflow
			real (kind = GRID_SR), intent(in)				    :: total_inflow
			logical, intent(in)			                        :: is_dirichlet

            rhs = rhs + cfg%r_inflow * inflow / total_inflow

            !limit the accumulated saturation at the injection well
            saturation = min(saturation, 1.0_SR)

            if (is_dirichlet) then
                rhs = 0.0_SR
            end if
		end subroutine
	END MODULE
#endif
