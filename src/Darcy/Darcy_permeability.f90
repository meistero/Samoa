! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_permeability
		use SFC_edge_traversal

		use Samoa_darcy

        type num_traversal_data
        end type

        type(darcy_gv_saturation)				:: gv_saturation
        type(darcy_gv_rhs)				        :: gv_rhs

#		define _GT_NAME							t_darcy_permeability_traversal

#		if (_DARCY_FLOW_EDGE_SIZE > 0)
#			define _GT_EDGES
#		endif

#		define _GT_NODES

#		define _GT_ELEMENT_OP					element_op

#		define _GT_NODE_FIRST_TOUCH_OP		    node_first_touch_op

#		define _GT_NODE_MERGE_OP		        node_merge_op

#		define _GT_EDGE_MPI_TYPE

#		include "SFC_generic_traversal_ringbuffer.f90"


        subroutine create_edge_mpi_type(mpi_edge_type)
            integer, intent(out)            :: mpi_edge_type

            type(t_edge_data)               :: edge
            integer                         :: blocklengths(2), types(2), disps(2), i_error, extent

#           if defined(_MPI)
                blocklengths(1) = 1
                blocklengths(2) = 1

                disps(1) = 0
                disps(2) = sizeof(edge)

                types(1) = MPI_LB
                types(2) = MPI_UB

                call MPI_Type_struct(2, blocklengths, disps, types, mpi_edge_type, i_error); assert_eq(i_error, 0)
                call MPI_Type_commit(mpi_edge_type, i_error); assert_eq(i_error, 0)

                call MPI_Type_extent(mpi_edge_type, extent, i_error); assert_eq(i_error, 0)
                assert_eq(sizeof(edge), extent)

                call MPI_Type_size(mpi_edge_type, extent, i_error); assert_eq(i_error, 0)
                assert_eq(0, extent)
#           endif
        end subroutine

		!*******************************
		!Geometry operators
		!*******************************

		elemental subroutine node_first_touch_op(traversal, section, node)
 			type(t_darcy_permeability_traversal), intent(in)                     :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			:: node

			node%data_pers%rhs = 0.0d0
		end subroutine

		subroutine element_op(traversal, section, element)
 			type(t_darcy_permeability_traversal), intent(inout)		:: traversal
 			type(t_grid_section), intent(inout)						:: section
			type(t_element_base), intent(inout), target				:: element

			real (kind = GRID_SR)		:: saturation(_DARCY_FLOW_SIZE)
			real (kind = GRID_SR)		:: rhs(_DARCY_P_SIZE)

			call gv_saturation%read(element, saturation)

			!call element operator
			call alpha_volume_op(element, saturation, element%cell%data_pers%base_permeability, element%cell%data_pers%permeability, rhs)

			call gv_rhs%add(element, rhs)
		end subroutine

		elemental subroutine node_merge_op(local_node, neighbor_node)
			type(t_node_data), intent(inout)		:: local_node
 			type(t_node_data), intent(in)		    :: neighbor_node

			local_node%data_pers%rhs = local_node%data_pers%rhs + neighbor_node%data_pers%rhs
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine alpha_volume_op(element, saturation, base_permeability, permeability, rhs)
			type(t_element_base), intent(inout)				                    :: element
			real (kind = GRID_SR), intent(in)		                            :: saturation(:)
			real (kind = GRID_SR), intent(in)									:: base_permeability
			real (kind = GRID_SR), intent(out)									:: permeability
			real (kind = GRID_SR), intent(out)									:: rhs(:)

			real (kind = GRID_SR), parameter            :: Dx(3, 3) = reshape([1.0d0/8.0d0, 1.0d0/4.0d0, 1.0d0/8.0d0, -1.0d0/8.0d0, -1.0d0/4.0d0, -1.0d0/8.0d0, 0.0d0, 0.0d0, 0.0d0], [3, 3])
			real (kind = GRID_SR), parameter            :: Dy(3, 3) = reshape([0.0d0, 0.0d0, 0.0d0, -1.0d0/8.0d0, -1.0d0/4.0d0, -1.0d0/8.0d0, 1.0d0/8.0d0, 1.0d0/4.0d0, 1.0d0/8.0d0], [3, 3])
			real (kind = GRID_SR)					    :: g_local(2), x(2), r_lambda_w(3), r_lambda_n(3)
			integer                                     :: i

			r_lambda_w = (saturation * saturation) / cfg%r_nu_w
			r_lambda_n = (1.0_GRID_SR - saturation) * (1.0_GRID_SR - saturation) / cfg%r_nu_n

			permeability = base_permeability * dot_product([0.25_GRID_SR, 0.5_GRID_SR, 0.25_GRID_SR], r_lambda_w + r_lambda_n)
            g_local = samoa_world_to_barycentric_normal(element%transform_data, g)

            rhs = g_local(1) * matmul(cfg%r_rho_w * r_lambda_w + cfg%r_rho_n * r_lambda_n, Dx) + g_local(2) * matmul(cfg%r_rho_w * r_lambda_w + cfg%r_rho_n * r_lambda_n, Dy)
            rhs = base_permeability * rhs
		end subroutine
	END MODULE
#endif
