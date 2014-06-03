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
            grid%r_dt = cfg%courant_number * cfg%scaling * cfg%r_phi * cfg%r_nu_w * get_edge_size(cfg%i_max_depth) / (2.0_GRID_SR * grid%u_max)
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
			real(kind = GRID_SR), dimension(_DARCY_FLOW_SIZE)		:: flux
			real(kind = GRID_SR), dimension(_DARCY_FLOW_SIZE)		:: volume

			call gv_u%read_from_element(element, u)
			call gv_saturation%read(element, saturation)

			!call volume operator
			call alpha_volume_op(element, u, saturation, volume, flux)

			call gv_volume%add(element, volume)
			call gv_flux%add(element, flux)
		end subroutine

		subroutine element_op(traversal, section, element)
 			type(t_darcy_transport_eq_traversal), intent(inout)	    :: traversal
 			type(t_grid_section), intent(inout)						:: section
			type(t_element_base), intent(inout), target				:: element

			real(kind = GRID_SR), dimension(2, _DARCY_U_SIZE)		:: u
			real(kind = GRID_SR), dimension(_DARCY_FLOW_SIZE)		:: saturation
			real(kind = GRID_SR), dimension(_DARCY_FLOW_SIZE)		:: flux
			real(kind = GRID_SR), dimension(_DARCY_FLOW_SIZE)		:: volume

			integer :: i
			real(kind = GRID_SR) :: u_dot_n
			real(kind = GRID_SR) :: edge_length

			call gv_u%read_from_element(element, u)
			call gv_saturation%read(element, saturation)

			!call volume operator
			call alpha_volume_op(element, u, saturation, volume, flux)

            !outflow on the right
            if (element%nodes(2)%ptr%position(1) == 1.0 .and. element%nodes(1)%ptr%position(1) == 1.0) then
                flux(2) = flux(2) - u(1, 1) * (saturation(2) * saturation(2)) * 0.5_GRID_SR * element%cell%geometry%get_leg_size() / cfg%r_nu_w
                flux(1) = flux(1) - u(1, 1) * (saturation(1) * saturation(1)) * 0.5_GRID_SR * element%cell%geometry%get_leg_size() / cfg%r_nu_w
            else if (element%nodes(2)%ptr%position(1) == 1.0 .and. element%nodes(3)%ptr%position(1) == 1.0) then
                flux(2) = flux(2) - u(1, 1) * (saturation(2) * saturation(2)) * 0.5_GRID_SR * element%cell%geometry%get_leg_size() / cfg%r_nu_w
                flux(3) = flux(3) - u(1, 1) * (saturation(3) * saturation(3)) * 0.5_GRID_SR * element%cell%geometry%get_leg_size() / cfg%r_nu_w
            else if (element%nodes(1)%ptr%position(1) == 1.0 .and. element%nodes(3)%ptr%position(1) == 1.0) then
                flux(1) = flux(1) - u(1, 1) * (saturation(1) * saturation(1)) * 0.5_GRID_SR * element%cell%geometry%get_hypo_size() / cfg%r_nu_w
                flux(3) = flux(3) - u(1, 1) * (saturation(3) * saturation(3)) * 0.5_GRID_SR * element%cell%geometry%get_hypo_size() / cfg%r_nu_w
            end if

            !outflow on the top
            if (element%nodes(2)%ptr%position(2) == 1.0 .and. element%nodes(1)%ptr%position(2) == 1.0) then
                flux(2) = flux(2) - u(2, 1) * (saturation(2) * saturation(2)) * 0.5_GRID_SR * element%cell%geometry%get_leg_size() / cfg%r_nu_w
                flux(1) = flux(1) - u(2, 1) * (saturation(1) * saturation(1)) * 0.5_GRID_SR * element%cell%geometry%get_leg_size() / cfg%r_nu_w
            else if (element%nodes(2)%ptr%position(2) == 1.0 .and. element%nodes(3)%ptr%position(2) == 1.0) then
                flux(2) = flux(2) - u(2, 1) * (saturation(2) * saturation(2)) * 0.5_GRID_SR * element%cell%geometry%get_leg_size() / cfg%r_nu_w
                flux(3) = flux(3) - u(2, 1) * (saturation(3) * saturation(3)) * 0.5_GRID_SR * element%cell%geometry%get_leg_size() / cfg%r_nu_w
            else if (element%nodes(1)%ptr%position(2) == 1.0 .and. element%nodes(3)%ptr%position(2) == 1.0) then
                flux(1) = flux(1) - u(2, 1) * (saturation(1) * saturation(1)) * 0.5_GRID_SR * element%cell%geometry%get_hypo_size() / cfg%r_nu_w
                flux(3) = flux(3) - u(2, 1) * (saturation(3) * saturation(3)) * 0.5_GRID_SR * element%cell%geometry%get_hypo_size() / cfg%r_nu_w
            end if

			!left and bottom boundary have symmetry conditions, so the fluxes cancel out

			call gv_volume%add(element, volume)
			call gv_flux%add(element, flux)
		end subroutine
		!first touches

		elemental subroutine node_first_touch_op(traversal, section, node)
 			type(t_darcy_transport_eq_traversal), intent(in)	    :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)				:: node

			call pre_dof_op(node%data_temp%flux, node%data_temp%volume)
		end subroutine

		!last touches

		elemental subroutine inner_node_last_touch_op(traversal, section, node)
 			type(t_darcy_transport_eq_traversal), intent(in)	    :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)				:: node
			integer (kind = GRID_SI)					:: i

			do i = 1, _DARCY_FLOW_NODE_SIZE
				call post_dof_op(section, node%data_pers%saturation(i), node%data_temp%flux(i), node%data_temp%volume(i))
			end do
		end subroutine

		elemental subroutine node_last_touch_op(traversal, section, node)
 			type(t_darcy_transport_eq_traversal), intent(in)	    :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)				:: node

			integer (kind = GRID_SI)					:: i

			if (node%position(1) > 0.0_GRID_SR .or. node%position(2) > 0.0_GRID_SR) then
                call post_dof_op(section, node%data_pers%saturation, node%data_temp%flux, node%data_temp%volume)
			end if
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		elemental subroutine pre_dof_op(flux, volume)
 			real (kind = GRID_SR), intent(out)		:: flux
			real (kind = GRID_SR), intent(out)		:: volume

			flux = 0.0_GRID_SR
			volume = 0.0_GRID_SR
		end subroutine

		subroutine alpha_volume_op(element, u, saturation, volume, flux)
			type(t_element_base), intent(inout)									:: element
			real (kind = GRID_SR), dimension(:, :), intent(in)					:: u
			real (kind = GRID_SR), dimension(:), intent(in)						:: saturation
			real (kind = GRID_SR), dimension(:), intent(out)					:: volume
			real (kind = GRID_SR), dimension(:), intent(out)					:: flux

			real (kind = GRID_SR), dimension(2)									:: r_u, n
			real (kind = GRID_SR)												:: r_u_dot_n, r_flux
			real (kind = GRID_SR), dimension(3)									:: r_lambda_w
			real (kind = GRID_SR)               								:: r_dual_edge_length
			integer (kind = GRID_SI)											:: i

			volume(:) = [0.25_GRID_SR, 0.50_GRID_SR, 0.25_GRID_SR] * element%cell%geometry%get_volume()

            r_lambda_w(:) = (saturation(:) * saturation(:)) / cfg%r_nu_w

			r_u(:) = u(:, 1)
			r_dual_edge_length = 0.5_GRID_SR * element%cell%geometry%get_leg_size()

			!compute an upwind flux

			n = samoa_barycentric_to_world_normal(element%transform_data, [1.0_GRID_SR, 0.0_GRID_SR])
			n = n / sqrt(DOT_PRODUCT(n, n))
			r_u_dot_n = DOT_PRODUCT(n, r_u)
			r_flux = (r_lambda_w(2) * max(r_u_dot_n, 0.0_GRID_SR) + r_lambda_w(1) * min(r_u_dot_n, 0.0_GRID_SR)) * r_dual_edge_length
			flux(1) = r_flux

			n = samoa_barycentric_to_world_normal(element%transform_data, [0.0_GRID_SR, 1.0_GRID_SR])
			n = n / sqrt(DOT_PRODUCT(n, n))
			r_u_dot_n = DOT_PRODUCT(n, r_u)
			r_flux = (r_lambda_w(2) * max(r_u_dot_n, 0.0_GRID_SR) + r_lambda_w(3) * min(r_u_dot_n, 0.0_GRID_SR)) * r_dual_edge_length
			flux(3) = r_flux

			flux(2) = -flux(1) - flux(3)
		end subroutine

		elemental subroutine post_dof_op(section, saturation, flux, volume)
			type(t_grid_section), intent(in)		:: section
 			real (kind = GRID_SR), intent(inout)	:: saturation
			real (kind = GRID_SR), intent(in)		:: flux
			real (kind = GRID_SR), intent(in)		:: volume

			real (kind = GRID_SR)					:: r_ds_dt

			r_ds_dt = flux / (cfg%r_phi * volume)
			saturation = saturation + section%r_dt * r_ds_dt
		end subroutine

		elemental subroutine node_merge_op(local_node, neighbor_node)
			type(t_node_data), intent(inout)		:: local_node
 			type(t_node_data), intent(in)		    :: neighbor_node

			local_node%data_temp%flux = local_node%data_temp%flux + neighbor_node%data_temp%flux
			local_node%data_temp%volume = local_node%data_temp%volume + neighbor_node%data_temp%volume
		end subroutine
	END MODULE
#endif
