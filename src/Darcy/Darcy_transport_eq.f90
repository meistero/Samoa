! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_transport_eq
		use SFC_edge_traversal
		use Darcy_grad_p

		use Samoa_darcy

		public compute_upwind_flux

        type num_traversal_data
            real (kind = GRID_SR)               :: prod_w(4), prod_n(4)
        end type

		type(darcy_gv_p)				        :: gv_p
		type(darcy_gv_saturation)				:: gv_saturation
		type(darcy_gv_volume)					:: gv_volume
		type(darcy_gv_flux)						:: gv_flux_w
		type(darcy_gv_d)						:: gv_flux_n

#		define _GT_NAME							t_darcy_transport_eq_traversal

#		define _GT_NODES

#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_PRE_TRAVERSAL_OP		        pre_traversal_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op

#		define _GT_ELEMENT_OP					element_op

#		define _GT_NODE_FIRST_TOUCH_OP		    node_first_touch_op
#		define _GT_NODE_LAST_TOUCH_OP		    node_last_touch_op
#		define _GT_NODE_REDUCE_OP		        node_reduce_op
#		define _GT_INNER_NODE_LAST_TOUCH_OP		inner_node_last_touch_op
#		define _GT_INNER_NODE_REDUCE_OP		    inner_node_reduce_op

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

            traversal%prod_w = 0.0_SR
            traversal%prod_n = 0.0_SR
		end subroutine

		subroutine pre_traversal_op(traversal, section)
			type(t_darcy_transport_eq_traversal), intent(inout)		:: traversal
			type(t_grid_section), intent(inout)							    :: section

            traversal%prod_w = 0.0_SR
            traversal%prod_n = 0.0_SR
		end subroutine

		subroutine post_traversal_grid_op(traversal, grid)
			type(t_darcy_transport_eq_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							    :: grid

            integer                :: i

            do i = 1, 4
                call reduce(traversal%prod_w(i), traversal%children%prod_w(i), MPI_SUM, .false.)
                call reduce(traversal%prod_n(i), traversal%children%prod_n(i), MPI_SUM, .false.)
            end do

            !accumulated production in bbl += dt in s * um^3/s * (6.28981077 bbl/m^3) * (cfg%scaling m/um)^3
            grid%prod_w_acc = grid%prod_w_acc + traversal%prod_w * grid%r_dt * (6.28981077_SR) * (cfg%scaling ** 3)
            grid%prod_n_acc = grid%prod_n_acc + traversal%prod_n * grid%r_dt * (6.28981077_SR) * (cfg%scaling ** 3)

            !production rate in bbl/d = um^3/s * (6.28981077 bbl/m^3) * (cfg%scaling m/um)^3 * (86400 s/d)
            grid%prod_w = traversal%prod_w * (6.28981077_SR * 86400.0_SR) * (cfg%scaling ** 3)
            grid%prod_n = traversal%prod_n * (6.28981077_SR * 86400.0_SR) * (cfg%scaling ** 3)

			grid%r_time = grid%r_time + grid%r_dt
		end subroutine

		!element

		subroutine element_op(traversal, section, element)
 			type(t_darcy_transport_eq_traversal), intent(inout)	    :: traversal
 			type(t_grid_section), intent(inout)						:: section
			type(t_element_base), intent(inout)			            :: element

#           if (_DARCY_LAYERS > 0)
                real(kind = GRID_SR)    :: p(_DARCY_LAYERS + 1, 3)
                real(kind = GRID_SR)    :: saturation(_DARCY_LAYERS + 1, 3)
                real(kind = GRID_SR)    :: flux_w(_DARCY_LAYERS + 1, 3)
                real(kind = GRID_SR)    :: flux_n(_DARCY_LAYERS + 1, 3)
                real(kind = GRID_SR)    :: volume(_DARCY_LAYERS + 1, 3)
                real(kind = GRID_SR)    :: porosity(_DARCY_LAYERS + 1)

                !Set porosity for each intermediate layer to the weighted average of lower and upper layer.
                !The top and bottom layers have half the volume of the inner layers.
                !We account for that by multiplying the porosity with 0.5 in the outer layers.

                porosity = 0.0_SR
                porosity(1 : _DARCY_LAYERS) = porosity(1 : _DARCY_LAYERS) + 0.5_SR * element%cell%data_pers%porosity
                porosity(2 : _DARCY_LAYERS + 1) = porosity(2 : _DARCY_LAYERS + 1) + 0.5_SR * element%cell%data_pers%porosity
#           else
                real(kind = GRID_SR)    :: p(3)
                real(kind = GRID_SR)    :: saturation(3)
                real(kind = GRID_SR)    :: flux_w(3)
                real(kind = GRID_SR)    :: flux_n(3)
                real(kind = GRID_SR)    :: volume(3)
                real(kind = GRID_SR)    :: porosity

                porosity = element%cell%data_pers%porosity
#           endif

			call gv_p%read_from_element(element, p)
			call gv_saturation%read_from_element(element, saturation)

			!call volume operator
			call compute_fluxes(element, p, saturation, volume, element%cell%data_pers%base_permeability, porosity, flux_w, flux_n)

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

		subroutine inner_node_last_touch_op(traversal, section, node)
 			type(t_darcy_transport_eq_traversal), intent(in)    :: traversal
 			type(t_grid_section), intent(in)				    :: section
			type(t_node_data), intent(inout)				    :: node

            call post_dof_op(section%r_dt, node%data_pers%saturation, node%data_temp%flux, node%data_pers%d, node%data_temp%volume)
		end subroutine

		elemental subroutine node_last_touch_op(traversal, section, node)
 			type(t_darcy_transport_eq_traversal), intent(in)    :: traversal
 			type(t_grid_section), intent(in)				    :: section
			type(t_node_data), intent(inout)				    :: node

            if (any(node%data_temp%is_dirichlet_boundary)) then
                call post_dof_op_correction(section%r_dt, node%data_pers%saturation, node%data_temp%flux, node%data_pers%d, node%data_temp%volume)
            else
                call post_dof_op(section%r_dt, node%data_pers%saturation, node%data_temp%flux, node%data_pers%d, node%data_temp%volume)
            end if
		end subroutine


		subroutine inner_node_reduce_op(traversal, section, node)
 			type(t_darcy_transport_eq_traversal), intent(inout)	    :: traversal
 			type(t_grid_section), intent(in)						:: section
			type(t_node_data), intent(inout)				        :: node

			! do nothing
		end subroutine


		subroutine node_reduce_op(traversal, section, node)
 			type(t_darcy_transport_eq_traversal), intent(inout)     :: traversal
 			type(t_grid_section), intent(in)				        :: section
			type(t_node_data), intent(inout)				        :: node

			real (kind = GRID_SR)  :: prod_w, prod_n
			integer :: i

            if (any(node%data_temp%is_dirichlet_boundary)) then
                prod_w = 0.0_SR
                prod_n = 0.0_SR

                do i = 1, _DARCY_LAYERS + 1
                    call reduce_op(node%data_pers%saturation(i), node%data_temp%flux(i), node%data_pers%d(i), prod_w, prod_n)
                end do

                if (node%position(2) > 0.5_SR) then
                    if (node%position(1) < 0.5_SR) then
                        traversal%prod_w(1) = traversal%prod_w(1) + prod_w
                        traversal%prod_n(1) = traversal%prod_n(1) + prod_n
                    else if (node%position(1) > 0.5_SR) then
                        traversal%prod_w(2) = traversal%prod_w(2) + prod_w
                        traversal%prod_n(2) = traversal%prod_n(2) + prod_n
                    end if
                else if (node%position(2) < 0.5_SR) then
                    if (node%position(1) < 0.5_SR) then
                        traversal%prod_w(4) = traversal%prod_w(4) + prod_w
                        traversal%prod_n(4) = traversal%prod_n(4) + prod_n
                    else if (node%position(1) > 0.5_SR) then
                        traversal%prod_w(3) = traversal%prod_w(3) + prod_w
                        traversal%prod_n(3) = traversal%prod_n(3) + prod_n
                    end if
                end if
            end if
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

        subroutine compute_fluxes(element, p, saturation, volume, base_permeability, porosity, flux_w, flux_n)
			type(t_element_base), intent(inout)	    :: element

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR), intent(in)	    :: p(:, :)
                real (kind = GRID_SR), intent(in)	    :: base_permeability(:,:)
                real (kind = GRID_SR), intent(in)	    :: porosity(:)
                real (kind = GRID_SR), intent(in)	    :: saturation(:, :)
                real (kind = GRID_SR), intent(out)	    :: volume(:, :), flux_w(:, :), flux_n(:, :)

                real (kind = GRID_SR)                   :: u_w(7), u_n(7)
                real (kind = GRID_SR)                   :: lambda_w(_DARCY_LAYERS + 1, 3), lambda_n(_DARCY_LAYERS + 1, 3)

                real (kind = GRID_SR)				    :: g_local(3), edge_length, surface, dz
                integer                                 :: i

                volume(:, 1) = cfg%dz * 0.25_SR * porosity(:) * element%cell%geometry%get_volume()
                volume(:, 2) = cfg%dz * 0.50_SR * porosity(:) * element%cell%geometry%get_volume()
                volume(:, 3) = cfg%dz * 0.25_SR * porosity(:) * element%cell%geometry%get_volume()

                lambda_w = (saturation * saturation) / cfg%r_nu_w
                lambda_n = (1.0_SR - saturation) * (1.0_SR - saturation) / cfg%r_nu_n

                flux_w = 0.0_SR
                flux_n = 0.0_SR

                !rotate g so it points in the right direction (no scaling!)
                g_local = g
                g_local(1:2) = samoa_world_to_barycentric_normal(element%transform_data, g_local(1:2))
                g_local(1:2) = g_local(1:2) / (element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))

                edge_length = element%cell%geometry%get_leg_size()
                surface = element%cell%geometry%get_volume()
                dz = cfg%dz

                do i = 1, _DARCY_LAYERS
                    call compute_velocity_1D(edge_length, 0.25_SR * edge_length * dz, base_permeability(i, 1), p(i, 2), p(i, 1), u_w(1), u_n(1), g_local(1))
                    call compute_velocity_1D(edge_length, 0.25_SR * edge_length * dz, base_permeability(i, 1), p(i, 2), p(i, 3), u_w(2), u_n(2), g_local(2))

                    call compute_velocity_1D(dz, 0.25_SR * surface, base_permeability(i, 2), p(i, 1), p(i + 1, 1), u_w(3), u_n(3), g_local(3))
                    call compute_velocity_1D(dz, 0.50_SR * surface, base_permeability(i, 2), p(i, 2), p(i + 1, 2), u_w(4), u_n(4), g_local(3))
                    call compute_velocity_1D(dz, 0.25_SR * surface, base_permeability(i, 2), p(i, 3), p(i + 1, 3), u_w(5), u_n(5), g_local(3))

                    call compute_velocity_1D(edge_length, 0.25_SR * edge_length * dz, base_permeability(i, 1), p(i + 1, 2), p(i + 1, 1), u_w(6), u_n(6), g_local(1))
                    call compute_velocity_1D(edge_length, 0.25_SR * edge_length * dz, base_permeability(i, 1), p(i + 1, 2), p(i + 1, 3), u_w(7), u_n(7), g_local(2))

                    !compute fluxes

                    call compute_flux(u_w, lambda_w(i:i+1, :), flux_w(i:i+1, :))
                    call compute_flux(u_n, lambda_n(i:i+1, :), flux_n(i:i+1, :))
                end do
#           else
                real (kind = GRID_SR), intent(in)	    :: p(:)
                real (kind = GRID_SR), intent(in)	    :: base_permeability
                real (kind = GRID_SR), intent(in)	    :: porosity
                real (kind = GRID_SR), intent(in)	    :: saturation(:)
                real (kind = GRID_SR), intent(out)	    :: volume(:), flux_w(:), flux_n(:)

                real (kind = GRID_SR)                   :: lambda_w(3), lambda_n(3)
                real (kind = GRID_SR)					:: g_local(2), edge_length
                real (kind = SR)                        :: u_w(2), u_n(2), u_norm

                volume = [0.25_SR, 0.50_SR, 0.25_SR] * porosity * element%cell%geometry%get_volume()

                lambda_w = (saturation * saturation) / cfg%r_nu_w
                lambda_n = (1.0_SR - saturation) * (1.0_SR - saturation) / cfg%r_nu_n

                flux_w = 0.0_SR
                flux_n = 0.0_SR

                !rotate g so it points in the right direction (no scaling!)
                g_local = g
                g_local(1:2) = samoa_world_to_barycentric_normal(element%transform_data, g_local(1:2))
                g_local(1:2) = g_local(1:2) / (element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))

                edge_length = element%cell%geometry%get_leg_size()

                !compute velocities

                call compute_velocity_1D(edge_length, 0.5_SR * edge_length, base_permeability, p(2), p(1), u_w(1), u_n(1), g_local(1))
                call compute_velocity_1D(edge_length, 0.5_SR * edge_length, base_permeability, p(2), p(3), u_w(2), u_n(2), g_local(2))

                !compute fluxes

                call compute_flux(u_w, lambda_w, flux_w)
                call compute_flux(u_n, lambda_n, flux_n)
#           endif

            !Careful: the FEM pressure solution implies that u = 0 on the boundary.
            !If we define an outflow condition in the FV step, we will get a non-zero divergence.
		end subroutine

		subroutine compute_flux(u, lambda, flux)
            real (kind = GRID_SR), intent(inout)       :: u(:)

            !Upwind and F-Wave solver are identical except for the treatment of boundaries.

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR), intent(in)       :: lambda(:, :)
                real (kind = GRID_SR), intent(out)	    :: flux(:, :)

#               if defined(_UPWIND_FLUX)
                    call compute_upwind_flux(u(1), lambda(1, 2), lambda(1, 1), flux(1, 2), flux(1, 1))
                    call compute_upwind_flux(u(2), lambda(1, 2), lambda(1, 3), flux(1, 2), flux(1, 3))

                    call compute_upwind_flux(u(3), lambda(1, 1), lambda(2, 1), flux(1, 1), flux(2, 1))
                    call compute_upwind_flux(u(4), lambda(1, 2), lambda(2, 2), flux(1, 2), flux(2, 2))
                    call compute_upwind_flux(u(5), lambda(1, 3), lambda(2, 3), flux(1, 3), flux(2, 3))

                    call compute_upwind_flux(u(6), lambda(2, 2), lambda(2, 1), flux(2, 2), flux(2, 1))
                    call compute_upwind_flux(u(7), lambda(2, 2), lambda(2, 3), flux(2, 2), flux(2, 3))
#               elif defined(_FWAVE_FLUX)
                    !Yeah, no..
#                   error Not yet implemented!
#               endif
#           else
                real (kind = GRID_SR), intent(in)       :: lambda(:)
                real (kind = GRID_SR), intent(out)	    :: flux(:)

#               if defined(_UPWIND_FLUX)
                    call compute_upwind_flux(u(1), lambda(2), lambda(1), flux(2), flux(1))
                    call compute_upwind_flux(u(2), lambda(2), lambda(3), flux(2), flux(3))
#               elif defined(_FWAVE_FLUX)
                    call compute_fwave_flux_dlambda(u(1), lambda(2), lambda(1), flux(2), flux(1))
                    call compute_fwave_flux_dlambda(u(2), lambda(2), lambda(3), flux(2), flux(3))

                    call compute_fwave_flux_du(-u(2), 0.0_SR, lambda(1), flux(1))
                    call compute_fwave_flux_du(u(1) + u(2), 0.0_SR, lambda(1), flux(1))

                    call compute_fwave_flux_du(-u(1), 0.0_SR, lambda(3), flux(3))
                    call compute_fwave_flux_du(u(1) + u(2), 0.0_SR, lambda(3), flux(3))

                    call compute_fwave_flux_du(-u(1), 0.0_SR, lambda(2), flux(2))
                    call compute_fwave_flux_du(-u(2), 0.0_SR, lambda(2), flux(2))
#               endif
#           endif
        end subroutine

        subroutine compute_upwind_flux(u, lambdaL, lambdaR, fluxL, fluxR)
            real (kind = GRID_SR), intent(in)       :: u, lambdaL, lambdaR
            real (kind = GRID_SR), intent(out)	    :: fluxL, fluxR

            fluxL = fluxL + (lambdaL * max(u, 0.0_SR) + lambdaR * min(u, 0.0_SR))
            fluxR = fluxR - (lambdaL * max(u, 0.0_SR) + lambdaR * min(u, 0.0_SR))
        end subroutine

        subroutine compute_fwave_flux_dlambda(u, lambdaL, lambdaR, fluxL, fluxR)
            real (kind = GRID_SR), intent(in)       :: u, lambdaL, lambdaR
            real (kind = GRID_SR), intent(out)	    :: fluxL, fluxR

            fluxL = fluxL + (lambdaR - lambdaL) * min(u, 0.0_SR)
            fluxR = fluxR + (lambdaR - lambdaL) * max(u, 0.0_SR)
        end subroutine

        subroutine compute_fwave_flux_du(uL, uR, lambda, fluxR)
            real (kind = GRID_SR), intent(in)       :: uL, uR, lambda
            real (kind = GRID_SR), intent(out)	    :: fluxR

            !if the intermediate state velocity is > 0 then add the flux difference to the left cell
            !otherwise add it to the right cell

            !But: since left and right cell are actually identical, we don't need the branch and the intermediate state.
            !Exception: boundary cells.

            fluxR = fluxR + lambda * (uR - uL)
        end subroutine

        !> Update saturation
        !>
        !> Here we assume the total flux is divergence-free, so f_w + f_n = 0:
        !> S'_w := S_w + dt/V * f_w
        !> S'_n := S_n + dt/V * f_n
        !> \Rightarrow S'_w + S'_n = S_w + S_n + dt/V * (f_w + f_n) = S_w + S_n
		elemental subroutine post_dof_op(dt, saturation, flux_w, flux_n, volume)
			real (kind = GRID_SR), intent(in)		:: dt
			real (kind = GRID_SR), intent(inout)	:: saturation
			real (kind = GRID_SR), intent(in)		:: flux_w
			real (kind = GRID_SR), intent(in)		:: flux_n
			real (kind = GRID_SR), intent(in)		:: volume

            if (volume > 0.0_SR) then
                saturation = saturation - dt / volume * flux_w
            end if

            !assert_pure(flux_w + flux_n == 0.0_SR)
		end subroutine

        !> Update saturation and apply a divergence correction
        !>
        !> Eliminate divergence by removing excess mass from the flux:
        !> S'_w := S_w + dt/V * (f_w - lambda_w / (lambda_w + lambda_n) * (f_w + f_n))
        !> S'_n := S_n + dt/V * (f_n - lambda_n / (lambda_w + lambda_n) * (f_w + f_n))
        !> \Rightarrow S'_w + S'_n = S_w + S_n + dt/V * ((f_w + f_n) - (f_w + f_n)) = S_w + S_n
		elemental subroutine post_dof_op_correction(dt, saturation, flux_w, flux_n, volume)
			real (kind = GRID_SR), intent(in)		:: dt
			real (kind = GRID_SR), intent(inout)	:: saturation
			real (kind = GRID_SR), intent(in)		:: flux_w, flux_n
			real (kind = GRID_SR), intent(in)		:: volume

            real (kind = GRID_SR)                   :: lambda_w, lambda_n

            if (volume > 0.0_SR) then
                lambda_w = (saturation * saturation) / cfg%r_nu_w
                lambda_n = (1.0_SR - saturation) * (1.0_SR - saturation) / cfg%r_nu_n

                saturation = saturation - dt / volume * (flux_w - lambda_w / (lambda_w + lambda_n) * (flux_w + flux_n))
            end if

            !assert_pure(saturation .le. 1.0_SR)
		end subroutine

		!> Compute production rates at the wells
		pure subroutine reduce_op(saturation, flux_w, flux_n, prod_w, prod_n)
			real (kind = GRID_SR), intent(in)	    :: saturation
			real (kind = GRID_SR), intent(in)		:: flux_w, flux_n
			real (kind = GRID_SR), intent(inout)	:: prod_w, prod_n

            real (kind = GRID_SR)                   :: lambda_w, lambda_n

            lambda_w = (saturation * saturation) / cfg%r_nu_w
            lambda_n = (1.0_SR - saturation) * (1.0_SR - saturation) / cfg%r_nu_n

            !Water and oil production rate:
            prod_w = prod_w - lambda_w / (lambda_w + lambda_n) * (flux_w + flux_n)
            prod_n = prod_n - lambda_n / (lambda_w + lambda_n) * (flux_w + flux_n)
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
