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
		type(darcy_gv_flux)						:: gv_flux_w
		type(darcy_gv_d)						:: gv_flux_n

#		define _GT_NAME							t_darcy_transport_eq_traversal

#		define _GT_NODES

#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op

#		define _GT_ELEMENT_OP					element_op

#		define _GT_NODE_FIRST_TOUCH_OP		    node_first_touch_op
#		define _GT_NODE_LAST_TOUCH_OP		    node_last_touch_op
!#		define _GT_INNER_NODE_LAST_TOUCH_OP		inner_node_last_touch_op

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

		subroutine element_op(traversal, section, element)
 			type(t_darcy_transport_eq_traversal), intent(inout)	    :: traversal
 			type(t_grid_section), intent(inout)						:: section
			type(t_element_base), intent(inout), target				:: element

			real(kind = GRID_SR)    :: u(2)
			real(kind = GRID_SR)    :: saturation(3)
			real(kind = GRID_SR)    :: flux_w(3)
			real(kind = GRID_SR)    :: flux_n(3)
			real(kind = GRID_SR)    :: volume(3)

			call gv_u%read_from_element(element, u)
			call gv_saturation%read_from_element(element, saturation)

			!call volume operator
			call compute_fluxes(element, u, saturation, volume, flux_w, flux_n)

			call gv_volume%add_to_element(element, volume)
			call gv_flux_w%add_to_element(element, flux_w)
			call gv_flux_n%add_to_element(element, flux_n)
		end subroutine

		!first touches

		elemental subroutine node_first_touch_op(traversal, section, node)
 			type(t_darcy_transport_eq_traversal), intent(in)	    :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)				:: node

			call pre_dof_op(node%data_temp%flux, node%data_pers%d, node%data_temp%volume)
		end subroutine

		!last touches

		elemental subroutine inner_node_last_touch_op(traversal, section, node)
 			type(t_darcy_transport_eq_traversal), intent(in)	    :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)				:: node

            call post_dof_op(section%r_dt, node%data_pers%saturation, node%data_temp%flux, node%data_pers%d, node%data_temp%volume)
		end subroutine

		elemental subroutine node_last_touch_op(traversal, section, node)
 			type(t_darcy_transport_eq_traversal), intent(in)    :: traversal
 			type(t_grid_section), intent(in)				    :: section
			type(t_node_data), intent(inout)				    :: node

            call post_dof_op_correction(section%r_dt, node%data_pers%saturation, node%data_temp%flux, node%data_pers%d, node%data_temp%volume)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		elemental subroutine pre_dof_op(flux_w, flux_n, volume)
 			real (kind = GRID_SR), intent(out)		:: flux_w
 			real (kind = GRID_SR), intent(out)		:: flux_n
			real (kind = GRID_SR), intent(out)		:: volume

			flux_w = 0.0_SR
			flux_n = 0.0_SR
			volume = 0.0_SR
		end subroutine

        subroutine compute_fluxes(element, u, saturation, volume, flux_w, flux_n)
			type(t_element_base), intent(inout)	    :: element
			real (kind = GRID_SR), intent(in)	    :: u(:), saturation(:)
			real (kind = GRID_SR), intent(out)	    :: volume(:), flux_w(:), flux_n(:)

            real (kind = GRID_SR)                   :: u_w(2), u_n(2), r_dual_edge_length
            real (kind = GRID_SR)                   :: lambda_w(3), lambda_n(3)
            integer                                 :: i

            volume(:) = [0.25_SR, 0.50_SR, 0.25_SR] * cfg%scaling * cfg%scaling * element%cell%data_pers%porosity * element%cell%geometry%get_volume()

            lambda_w = (saturation * saturation) / cfg%r_nu_w
            lambda_n = (1.0_SR - saturation) * (1.0_SR - saturation) / cfg%r_nu_n

			!compute fluxes

			u_w = samoa_world_to_barycentric_normal(element%transform_data, u)
			u_w = u_w / (element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))
			u_n = u_w
			assert_le(norm2(u), (1.0_SR + 1.0e2_SR * epsilon(1.0_SR)) * norm2(u_w))
			assert_le(norm2(u_w), (1.0_SR + 1.0e2_SR * epsilon(1.0_SR)) * norm2(u))

            flux_w(1) = lambda_w(2) * max(u_w(1), 0.0_SR) + lambda_w(1) * min(u_w(1), 0.0_SR)
            flux_w(3) = lambda_w(2) * max(u_w(2), 0.0_SR) + lambda_w(3) * min(u_w(2), 0.0_SR)
            flux_w(2) = -flux_w(1) - flux_w(3)

            flux_n(1) = lambda_n(2) * max(u_n(1), 0.0_SR) + lambda_n(1) * min(u_n(1), 0.0_SR)
            flux_n(3) = lambda_n(2) * max(u_n(2), 0.0_SR) + lambda_n(3) * min(u_n(2), 0.0_SR)
            flux_n(2) = -flux_n(1) - flux_n(3)

            !Careful: the FEM pressure solution implies that u = 0 on the boundary.
            !If we define an outflow condition in the FV step, we will get a non-zero divergence.
		end subroutine

        !> Update saturation
		elemental subroutine post_dof_op(dt, saturation, flux_w, flux_n, volume)
			real (kind = GRID_SR), intent(in)		:: dt
			real (kind = GRID_SR), intent(inout)	:: saturation
			real (kind = GRID_SR), intent(in)		:: flux_w
			real (kind = GRID_SR), intent(in)		:: flux_n
			real (kind = GRID_SR), intent(in)		:: volume

            if (volume > 0.0_SR) then
                saturation = saturation + dt / volume * flux_w
            end if

            !assert_pure(flux_w + flux_n == 1.0_SR)
		end subroutine

        !> Update saturation and apply a divergence correction
        !> EITHER: Eliminate divergence by removing excess mass from the saturation:
        !> S_w \leftarrow S_w' / (S_w' + S_n') = (S_w + dt/V * f_w) / ((S_w + S_n) + dt/V * (f_w + f_n))
        !> If f_w + f_n = 0, then S_w \leftarrow S_w'
        !> saturation = (saturation + dt / volume * flux_w) / (1.0_SR + dt / volume * (flux_w + flux_n))
        !>
        !> OR: Eliminate divergence by removing excess mass from the flux:
        !> S_w \leftarrow S_w + dt/V * (f_w - lambda_w / (lambda_w + lambda_n) * (f_w + f_n))
        !> S_n \leftarrow S_n + dt/V * (f_n - lambda_n / (lambda_w + lambda_n) * (f_w + f_n))
        !> \Rightarrow S_w + S_n \leftarrow S_w + S_n + dt/V * ((f_w + f_n) - (f_w + f_n))
        !>
        !> The last works better, maybe because it prevents an unphysical cell state instead of fixing it
        !> Saturation conservation is violated by both approaches however.
		elemental subroutine post_dof_op_correction(dt, saturation, flux_w, flux_n, volume)
			real (kind = GRID_SR), intent(in)		:: dt
			real (kind = GRID_SR), intent(inout)	:: saturation
			real (kind = GRID_SR), intent(in)		:: flux_w
			real (kind = GRID_SR), intent(in)		:: flux_n
			real (kind = GRID_SR), intent(in)		:: volume

            real (kind = GRID_SR)                   :: lambda_w, lambda_n

            if (volume > 0.0_SR) then
                lambda_w = (saturation * saturation) / cfg%r_nu_w
                lambda_n = (1.0_SR - saturation) * (1.0_SR - saturation) / cfg%r_nu_n

                saturation = saturation + dt / volume * (flux_w - lambda_w / (lambda_w + lambda_n) * (flux_w + flux_n))
            end if

            !assert_pure(saturation .le. 1.0_SR)
		end subroutine

		elemental subroutine node_merge_op(local_node, neighbor_node)
			type(t_node_data), intent(inout)		:: local_node
 			type(t_node_data), intent(in)		    :: neighbor_node

			local_node%data_temp%flux = local_node%data_temp%flux + neighbor_node%data_temp%flux
			local_node%data_pers%d = local_node%data_pers%d + neighbor_node%data_pers%d
			local_node%data_temp%volume = local_node%data_temp%volume + neighbor_node%data_temp%volume
		end subroutine
	END MODULE
#endif
