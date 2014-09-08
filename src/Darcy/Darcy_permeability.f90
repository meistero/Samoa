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

#		define _GT_NAME							t_darcy_permeability_traversal

#		if (_DARCY_FLOW_EDGE_SIZE > 0)
#			define _GT_EDGES
#		endif

#		define _GT_NODES

#		define _GT_ELEMENT_OP					element_op

#		define _GT_NODE_MPI_TYPE
#		define _GT_EDGE_MPI_TYPE

#		include "SFC_generic_traversal_ringbuffer.f90"

        subroutine create_node_mpi_type(mpi_node_type)
            integer, intent(out)            :: mpi_node_type

            type(t_node_data)               :: node
            integer                         :: blocklengths(2), types(2), disps(2), i_error, extent

#           if defined(_MPI)
                blocklengths(1) = 1
                blocklengths(2) = 1

                disps(1) = 0
                disps(2) = sizeof(node)

                types(1) = MPI_LB
                types(2) = MPI_UB

                call MPI_Type_struct(2, blocklengths, disps, types, mpi_node_type, i_error); assert_eq(i_error, 0)
                call MPI_Type_commit(mpi_node_type, i_error); assert_eq(i_error, 0)

                call MPI_Type_extent(mpi_node_type, extent, i_error); assert_eq(i_error, 0)
                assert_eq(sizeof(node), extent)

                call MPI_Type_size(mpi_node_type, extent, i_error); assert_eq(i_error, 0)
                assert_eq(0, extent)
#           endif
        end subroutine

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

		subroutine element_op(traversal, section, element)
 			type(t_darcy_permeability_traversal), intent(inout)		:: traversal
 			type(t_grid_section), intent(inout)						:: section
			type(t_element_base), intent(inout), target				:: element

			real (kind = GRID_SR), dimension(_DARCY_FLOW_SIZE)		:: saturation

			call gv_saturation%read(element, saturation)

			!call element operator
			call alpha_volume_op(saturation, element%cell%data_pers%base_permeability, element%cell%data_pers%permeability)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine alpha_volume_op(saturation, base_permeability, permeability)
			real (kind = GRID_SR), dimension(_DARCY_FLOW_SIZE), intent(in)		:: saturation
			real (kind = GRID_SR), intent(in)									:: base_permeability
			real (kind = GRID_SR), intent(out)									:: permeability

			real (kind = GRID_SR)												:: r_lambda_n, r_lambda_w

			r_lambda_n = sum([0.25_GRID_SR, 0.5_GRID_SR, 0.25_GRID_SR] * (1.0_GRID_SR - saturation) * (1.0_GRID_SR - saturation))
			r_lambda_w = sum([0.25_GRID_SR, 0.5_GRID_SR, 0.25_GRID_SR] * saturation * saturation)

			permeability = base_permeability * (r_lambda_n + cfg%r_rel_permeability * r_lambda_w)
		end subroutine
	END MODULE
#endif
