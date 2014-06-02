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

#		define _GT_INNER_NODE_FIRST_TOUCH_OP	inner_node_first_touch_op
#		define _GT_NODE_FIRST_TOUCH_OP		    node_first_touch_op

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

		elemental subroutine inner_node_first_touch_op(traversal, section, node)
 			type(t_darcy_permeability_traversal), intent(in)                     :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			:: node

			call inner_flow_pre_dof_op(node%data_pers%rhs)
		end subroutine

		elemental subroutine node_first_touch_op(traversal, section, node)
 			type(t_darcy_permeability_traversal), intent(in)                     :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			:: node

			integer :: i

			do i = 1, _DARCY_FLOW_NODE_SIZE
				call flow_pre_dof_op(node%position, node%data_pers%rhs(i))
			end do
		end subroutine

		elemental subroutine inner_flow_pre_dof_op(rhs)
			real (kind = GRID_SR), intent(out)					:: rhs

			rhs = 0.0_GRID_SR
		end subroutine

		pure subroutine flow_pre_dof_op(pos, rhs)
			real (kind = GRID_SR), dimension(2), intent(in)		:: pos
			real (kind = GRID_SR), intent(out)					:: rhs

			rhs = 0.0_GRID_SR
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

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine alpha_volume_op(element, saturation, base_permeability, permeability, rhs)
			type(t_element_base), intent(inout)				                    :: element
			real (kind = GRID_SR), dimension(_DARCY_FLOW_SIZE), intent(in)		:: saturation
			real (kind = GRID_SR), intent(in)									:: base_permeability
			real (kind = GRID_SR), intent(out)									:: permeability
			real (kind = GRID_SR), intent(out)									:: rhs(:)

			real (kind = GRID_SR), parameter                                    :: g(2) = [0.0d0, -9.81d0]
			real (kind = GRID_SR)												:: x(2), grad_psi(2), r_lambda_w, r_lambda_n
			integer     :: i

			r_lambda_w = sum([0.25_GRID_SR, 0.5_GRID_SR, 0.25_GRID_SR] * saturation * saturation) / cfg%r_nu_w
			r_lambda_n = sum([0.25_GRID_SR, 0.5_GRID_SR, 0.25_GRID_SR] * (1.0_GRID_SR - saturation) * (1.0_GRID_SR - saturation)) / cfg%r_nu_n

			permeability = base_permeability * (r_lambda_w + r_lambda_n)

            do i = 1, 3
                x = samoa_basis_p_get_dof_coords(i)
                grad_psi = [samoa_basis_p_d_dx(x), samoa_basis_p_d_dy(x)]
                grad_psi = samoa_barycentric_to_world_normal(element%transform_data, grad_psi)
                rhs(i) = (cfg%r_rho_w * r_lambda_w + cfg%r_rho_n * r_lambda_n) * dot_product(g, grad_psi) !this should be the integral over lambda_t * g * grad psi
            end do
		end subroutine
	END MODULE
#endif
