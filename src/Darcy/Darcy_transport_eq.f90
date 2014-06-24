! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_transport_eq
		use SFC_edge_traversal

		use Samoa_darcy

        type num_traversal_data
        end type

		type(darcy_gv_u)						:: gv_u
		type(darcy_gv_saturation)				:: gv_saturation
		type(darcy_gv_volume)					:: gv_volume
		type(darcy_gv_flux)						:: gv_flux
		type(darcy_gv_d)						:: gv_diverg

#		define _GT_NAME							t_darcy_transport_eq_traversal

#		if (_DARCY_U_EDGE_SIZE > 0 || _DARCY_FLOW_EDGE_SIZE > 0)
#			define _GT_EDGES
#		endif

#		define _GT_NODES

#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op

#		define _GT_INNER_ELEMENT_OP				inner_element_op
#		define _GT_ELEMENT_OP					element_op

#		define _GT_NODE_FIRST_TOUCH_OP		    node_first_touch_op
#		define _GT_NODE_LAST_TOUCH_OP		    node_last_touch_op
#		define _GT_INNER_NODE_LAST_TOUCH_OP		inner_node_last_touch_op

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

		subroutine pre_traversal_grid_op(traversal, grid)
			type(t_darcy_transport_eq_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							    :: grid

            !Dual cells have a size of of edge_length, the maximum Eigenvalue of the system is (2 S_max u_max) / (\Phi nu) = (2 u_max) / (\Phi nu)
            !This gives an upper bound of \Delta t \leq (\Phi nu edge_length) / (2 u_max)
			call scatter(grid%r_dt, grid%sections%elements_alloc%r_dt)
		end subroutine

		subroutine post_traversal_grid_op(traversal, grid)
			type(t_darcy_transport_eq_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							    :: grid

			grid%r_time = grid%r_time + grid%r_dt
		end subroutine

		!element

		subroutine inner_element_op(traversal, section, element)
 			type(t_darcy_transport_eq_traversal), intent(inout)	    :: traversal
 			type(t_grid_section), intent(inout)						:: section
			type(t_element_base), intent(inout), target				:: element

			real(kind = GRID_SR), dimension(2, _DARCY_U_SIZE)		:: u
			real(kind = GRID_SR), dimension(_DARCY_FLOW_SIZE)		:: saturation
			real(kind = GRID_SR), dimension(_DARCY_FLOW_SIZE)		:: flux, diverg
			real(kind = GRID_SR), dimension(_DARCY_FLOW_SIZE)		:: volume

			call gv_u%read_from_element(element, u)
			call gv_saturation%read(element, saturation)

			!call volume operator
			call compute_fluxes(element, u, saturation, volume, flux, diverg)

			call gv_volume%add(element, volume)
			call gv_flux%add(element, flux)
			call gv_diverg%add(element, diverg)
		end subroutine

		subroutine element_op(traversal, section, element)
 			type(t_darcy_transport_eq_traversal), intent(inout)	    :: traversal
 			type(t_grid_section), intent(inout)						:: section
			type(t_element_base), intent(inout), target				:: element

			real(kind = GRID_SR), dimension(2, _DARCY_U_SIZE)		:: u
			real(kind = GRID_SR), dimension(_DARCY_FLOW_SIZE)		:: saturation
			real(kind = GRID_SR), dimension(_DARCY_FLOW_SIZE)		:: flux, diverg
			real(kind = GRID_SR), dimension(_DARCY_FLOW_SIZE)		:: volume

			call gv_u%read_from_element(element, u)
			call gv_saturation%read(element, saturation)

			!call volume operator
			call compute_fluxes(element, u, saturation, volume, flux, diverg)

			call gv_volume%add(element, volume)
			call gv_flux%add(element, flux)
			call gv_diverg%add(element, diverg)
		end subroutine

		!first touches

		elemental subroutine node_first_touch_op(traversal, section, node)
 			type(t_darcy_transport_eq_traversal), intent(in)	    :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)				:: node

			call pre_dof_op(node%data_temp%flux, node%data_temp%volume, node%data_pers%d)
		end subroutine

		!last touches

		elemental subroutine inner_node_last_touch_op(traversal, section, node)
 			type(t_darcy_transport_eq_traversal), intent(in)	    :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)				:: node

            !do saturation update without divergence correction
            call inner_post_dof_op(section%r_dt, node%data_pers%saturation, node%data_temp%flux, node%data_temp%volume)
		end subroutine

		elemental subroutine node_last_touch_op(traversal, section, node)
 			type(t_darcy_transport_eq_traversal), intent(in)	    :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)				:: node

            if (any(abs(node%position - 0.5_SR) == 0.5_SR)) then
                !do a saturation update with divergence correction
                call post_dof_op(section%r_dt, node%data_pers%saturation, node%data_temp%flux, node%data_temp%volume, node%data_pers%d)
            else
                !do a saturation update without divergence correction
                call inner_post_dof_op(section%r_dt, node%data_pers%saturation, node%data_temp%flux, node%data_temp%volume)
            end if
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		elemental subroutine pre_dof_op(flux, volume, diverg)
 			real (kind = GRID_SR), intent(out)		:: flux
			real (kind = GRID_SR), intent(out)		:: volume
			real (kind = GRID_SR), intent(out)		:: diverg

			flux = 0.0_GRID_SR
			volume = 0.0_GRID_SR
			diverg = 0.0_SR
		end subroutine

		subroutine compute_inner_fluxes(element, u_in, saturation, volume, flux)
			type(t_element_base), intent(inout)									:: element
			real (kind = GRID_SR), dimension(:, :), intent(in)					:: u_in
			real (kind = GRID_SR), dimension(:), intent(in)						:: saturation
			real (kind = GRID_SR), dimension(:), intent(out)					:: volume
			real (kind = GRID_SR), dimension(:), intent(out)					:: flux

			real (kind = GRID_SR)									:: u(2), n(2), r_dual_edge_length
            real (kind = GRID_SR)								    :: u_dot_n, lambda_w(3)
            real (kind = GRID_SR)               				    :: pos_prod(2), radius
            integer (kind = GRID_SI)							    :: i

			volume(:) = [0.25_GRID_SR, 0.50_GRID_SR, 0.25_GRID_SR] * cfg%scaling * cfg%scaling * element%cell%data_pers%porosity * element%cell%geometry%get_volume()

            lambda_w(:) = (saturation(:) * saturation(:)) / cfg%r_nu_w

			u(:) = u_in(:, 1)
			r_dual_edge_length = 0.5_GRID_SR * cfg%scaling * element%cell%geometry%get_leg_size()

			!compute an upwind flux

			n = samoa_barycentric_to_world_normal(element%transform_data, [1.0_GRID_SR, 0.0_GRID_SR])
			n = n / norm2(n)
			u_dot_n = dot_product(n, u)
			flux(1) = (lambda_w(2) * max(u_dot_n, 0.0_GRID_SR) + lambda_w(1) * min(u_dot_n, 0.0_GRID_SR)) * r_dual_edge_length

			n = samoa_barycentric_to_world_normal(element%transform_data, [0.0_GRID_SR, 1.0_GRID_SR])
			n = n / norm2(n)
			u_dot_n = dot_product(n, u)
			flux(3) = (lambda_w(2) * max(u_dot_n, 0.0_GRID_SR) + lambda_w(3) * min(u_dot_n, 0.0_GRID_SR)) * r_dual_edge_length

			flux(2) = -flux(1) - flux(3)
		end subroutine

        subroutine compute_fluxes(element, u_in, saturation, volume, flux, diverg)
			type(t_element_base), intent(inout)									:: element
			real (kind = GRID_SR), dimension(:, :), intent(in)					:: u_in
			real (kind = GRID_SR), dimension(:), intent(in)						:: saturation
			real (kind = GRID_SR), dimension(:), intent(out)	    :: volume
			real (kind = GRID_SR), dimension(:), intent(out)	    :: flux
			real (kind = GRID_SR), dimension(:), intent(out)	    :: diverg

			real (kind = GRID_SR)									:: u(2), n(2), r_dual_edge_length
            real (kind = GRID_SR)								    :: u_dot_n, lambda_w(3), lambda_n(3)
            real (kind = GRID_SR)               				    :: pos_prod(2), radius
            integer (kind = GRID_SI)							    :: i

			volume(:) = [0.25_GRID_SR, 0.50_GRID_SR, 0.25_GRID_SR] * cfg%scaling * cfg%scaling * element%cell%data_pers%porosity * element%cell%geometry%get_volume()

            lambda_w(:) = (saturation(:) * saturation(:)) / cfg%r_nu_w
            lambda_n(:) = (1.0_SR - saturation(:)) * (1.0_SR - saturation(:)) / cfg%r_nu_n

			u(:) = u_in(:, 1)
			r_dual_edge_length = 0.5_GRID_SR * cfg%scaling * element%cell%geometry%get_leg_size()

			!compute an upwind flux
			!also, accumulate the divergence on boundary dual cells

			n = samoa_barycentric_to_world_normal(element%transform_data, [1.0_GRID_SR, 0.0_GRID_SR])
			n = n / norm2(n)
			u_dot_n = dot_product(n, u)
			flux(1) = (lambda_w(2) * max(u_dot_n, 0.0_GRID_SR) + lambda_w(1) * min(u_dot_n, 0.0_GRID_SR)) * r_dual_edge_length
			diverg(1) = ((lambda_w(2) + lambda_n(2)) * max(u_dot_n, 0.0_GRID_SR) + (lambda_w(1) + lambda_n(1)) * min(u_dot_n, 0.0_GRID_SR)) * r_dual_edge_length

			n = samoa_barycentric_to_world_normal(element%transform_data, [0.0_GRID_SR, 1.0_GRID_SR])
			n = n / norm2(n)
			u_dot_n = dot_product(n, u)
			flux(3) = (lambda_w(2) * max(u_dot_n, 0.0_GRID_SR) + lambda_w(3) * min(u_dot_n, 0.0_GRID_SR)) * r_dual_edge_length
			diverg(3) = ((lambda_w(2) + lambda_n(2)) * max(u_dot_n, 0.0_GRID_SR) + (lambda_w(1) + lambda_n(1)) * min(u_dot_n, 0.0_GRID_SR)) * r_dual_edge_length

			flux(2) = -flux(1) - flux(3)
			diverg(2) = -diverg(1) - diverg(3)

            !Careful: the FEM pressure solution implies that u = 0 on the boundary.
            !If we defined an outflow condition in the FV step, we would get a non-zero divergence.
		end subroutine

		elemental subroutine inner_post_dof_op(dt, saturation, flux, volume)
			real (kind = GRID_SR), intent(in)		:: dt
			real (kind = GRID_SR), intent(inout)	:: saturation
			real (kind = GRID_SR), intent(in)		:: flux
			real (kind = GRID_SR), intent(in)		:: volume

            if (volume > 0.0_SR) then
                saturation = saturation + dt * flux / volume
            end if

            !assert_pure(saturation <= 1.0_SR)
		end subroutine

		elemental subroutine post_dof_op(dt, saturation, flux, volume, diverg)
			real (kind = GRID_SR), intent(in)		:: dt
 			real (kind = GRID_SR), intent(inout)	:: saturation
			real (kind = GRID_SR), intent(in)		:: flux
			real (kind = GRID_SR), intent(in)		:: volume
			real (kind = GRID_SR), intent(in)		:: diverg

            real (kind = GRID_SR)				    :: lambda_w, lambda_n

            if (volume > 0.0_SR) then
                lambda_w = (saturation * saturation) / cfg%r_nu_w
                lambda_n = (1.0_SR - saturation) * (1.0_SR - saturation) / cfg%r_nu_n

                saturation = saturation + dt * (flux - lambda_w / (lambda_w + lambda_n) * diverg) / volume
            end if

            !assert_pure(saturation <= 1.0_SR)
		end subroutine

		elemental subroutine node_merge_op(local_node, neighbor_node)
			type(t_node_data), intent(inout)		:: local_node
 			type(t_node_data), intent(in)		    :: neighbor_node

			local_node%data_temp%flux = local_node%data_temp%flux + neighbor_node%data_temp%flux
			local_node%data_temp%volume = local_node%data_temp%volume + neighbor_node%data_temp%volume
			local_node%data_pers%d = local_node%data_pers%d + neighbor_node%data_pers%d
		end subroutine
	END MODULE
#endif
