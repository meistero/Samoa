! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_grad_p
		use SFC_edge_traversal

		use Samoa_darcy

        type num_traversal_data
			real (kind = GRID_SR)				:: u_max					!< global maximum of u
        end type

		type(darcy_gv_p)						:: gv_p
		type(darcy_gv_u)						:: gv_u

#		define _GT_NAME							t_darcy_grad_p_traversal

#		if (_DARCY_P_EDGE_SIZE > 0 || _DARCY_FLOW_EDGE_SIZE > 0)
#			define _GT_EDGES
#		endif

#		define _GT_NODES

#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op
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

		subroutine post_traversal_grid_op(traversal, grid)
			type(t_darcy_grad_p_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							:: grid

			call reduce(grid%u_max, traversal%children%u_max, MPI_MAX, .true.)
			grid%u_max = sqrt(grid%u_max)
		end subroutine

		subroutine pre_traversal_op(traversal, section)
  			type(t_darcy_grad_p_traversal), intent(inout)		:: traversal
 			type(t_grid_section), intent(inout)				    :: section

			traversal%u_max = 0.0_GRID_SR
		end subroutine

		!*******************************
		!Geometry operators
		!*******************************

		subroutine element_op(traversal, section, element)
 			type(t_darcy_grad_p_traversal), intent(inout)				:: traversal
 			type(t_grid_section), intent(inout)				:: section
			type(t_element_base), intent(inout), target		:: element

			real (kind = GRID_SR), dimension(_DARCY_P_SIZE)	:: p
			real (kind = GRID_SR), dimension(2, _DARCY_U_SIZE)	:: u

			call gv_p%read(element, p)

			!call element operator
			call alpha_volume_op(traversal, element, p, u, element%cell%data_pers%base_permeability)

			call gv_u%write_to_element(element, u)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine alpha_volume_op(traversal, element, p, u, base_permeability)
 			type(t_darcy_grad_p_traversal), intent(inout)						:: traversal
			type(t_element_base), intent(inout)									:: element
			real (kind = GRID_SR), dimension(:), intent(in)						:: p
			real (kind = GRID_SR), dimension(:, :), intent(out)					:: u
			real (kind = GRID_SR), intent(in)									:: base_permeability

			real (kind = GRID_SR), parameter, dimension(2)						:: g(2) = [0.0, -9.81]

			integer (kind = GRID_SI)											:: i

			!define velocity by u = -k grad p
			do i = 1, _DARCY_U_SIZE
				u(:, i) = samoa_basis_p_gradient(samoa_basis_u_get_dof_coords(i), p)
				u(:, i) = -base_permeability * samoa_barycentric_to_world_normal(element%transform_data, u(:, i)) + cfg%r_rho_w * g
				traversal%u_max = max(traversal%u_max, DOT_PRODUCT(u(:, i), u(:, i)))
			end do

			!compute DoFs from the values of u
			u(1, :) = samoa_basis_u_values_to_dofs(u(1, :))
			u(2, :) = samoa_basis_u_values_to_dofs(u(2, :))
		end subroutine
	END MODULE
#endif
