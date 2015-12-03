! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_transport_eq
		use SFC_edge_traversal
		use Darcy_grad_p
		use Darcy_initialize_saturation

		use Samoa_darcy

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

            !In the 2D case we always assumed that the height of the domain is 1, when in fact it should be delta_z.
            !so multiply the rates by delta_z
#           if (_DARCY_LAYERS == 0)
                traversal%prod_w = traversal%prod_w * cfg%dz
                traversal%prod_n = traversal%prod_n * cfg%dz
#           endif

            !accumulated production in bbl
            grid%prod_w_acc = grid%prod_w_acc + (traversal%prod_w * grid%r_dt) / _BBL
            grid%prod_n_acc = grid%prod_n_acc + (traversal%prod_n * grid%r_dt) / _BBL

            !production rate in bbl/d
            grid%prod_w = traversal%prod_w / (_BBL / _D)
            grid%prod_n = traversal%prod_n / (_BBL / _D)

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

            if (any(node%data_temp%is_pressure_dirichlet_boundary)) then
                call post_dof_op_correction(section%r_dt, node%data_pers%saturation, node%data_temp%flux, node%data_pers%d, node%data_temp%volume)
            else
                where(node%data_temp%is_saturation_dirichlet_boundary)
                    node%data_temp%flux = 0.0_SR
                end where

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

            if (any(node%data_temp%is_pressure_dirichlet_boundary)) then
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

                real (kind = GRID_SR)                   :: u_w(_DARCY_LAYERS, 7), u_n(_DARCY_LAYERS, 7)

                real (kind = GRID_SR)				    :: g_local(3), edge_length, surface, dz
                integer                                 :: i

                volume(:, 1) = cfg%dz * 0.25_SR * porosity(:) * element%cell%geometry%get_volume()
                volume(:, 2) = cfg%dz * 0.50_SR * porosity(:) * element%cell%geometry%get_volume()
                volume(:, 3) = cfg%dz * 0.25_SR * porosity(:) * element%cell%geometry%get_volume()

                flux_w = 0.0_SR
                flux_n = 0.0_SR

                !rotate g so it points in the right direction (no scaling!)
                g_local = cfg%g
                g_local(1:2) = samoa_world_to_barycentric_normal(element%transform_data, g_local(1:2))
                g_local(1:2) = g_local(1:2) / (element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))

                edge_length = element%cell%geometry%get_leg_size()
                surface = element%cell%geometry%get_volume()
                dz = cfg%dz

                call compute_base_fluxes_3D(p, base_permeability, edge_length, dz, 0.5_SR * edge_length * dz, surface, g_local, u_w, u_n)
                call compute_fluxes_3D(saturation, u_w, u_n, flux_w, flux_n)
#           else
                real (kind = GRID_SR), intent(in)	    :: p(:)
                real (kind = GRID_SR), intent(in)	    :: base_permeability
                real (kind = GRID_SR), intent(in)	    :: porosity
                real (kind = GRID_SR), intent(in)	    :: saturation(:)
                real (kind = GRID_SR), intent(out)	    :: volume(:), flux_w(:), flux_n(:)

                real (kind = GRID_SR)					:: g_local(2), edge_length
                real (kind = SR)                        :: u_w(2), u_n(2), u_norm

                volume = [0.25_SR, 0.50_SR, 0.25_SR] * porosity * element%cell%geometry%get_volume()

                flux_w = 0.0_SR
                flux_n = 0.0_SR

                !rotate g so it points in the right direction (no scaling!)
                g_local = cfg%g(1:2)
                g_local(1:2) = samoa_world_to_barycentric_normal(element%transform_data, g_local(1:2))
                g_local(1:2) = g_local(1:2) / (element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))

                edge_length = element%cell%geometry%get_leg_size()

                !compute fluxes

                call compute_base_fluxes_2D(p, base_permeability, edge_length, 0.5_SR * edge_length, g_local(1:2), u_w, u_n)
                call compute_fluxes_2D(saturation, u_w, u_n, flux_w, flux_n)
#           endif

            !Careful: the FEM pressure solution implies that u = 0 on the boundary.
            !If we define an outflow condition in the FV step, we will get a non-zero divergence.
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
            else
                saturation = max(0.0_SR, min(1.0_SR, saturation - flux_w))
            end if

            assert_pure(saturation .ge. 0.0_SR)
            assert_pure(saturation .le. 1.0_SR)
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
                lambda_w = l_w(saturation)
                lambda_n = l_n(saturation)

                saturation = saturation - dt / volume * (flux_w - lambda_w / (lambda_w + lambda_n) * (flux_w + flux_n))
            else
                saturation = max(0.0_SR, min(1.0_SR, saturation - flux_w))
            end if

            assert_pure(saturation .ge. 0.0_SR)
            assert_pure(saturation .le. 1.0_SR)
		end subroutine

		!> Compute production rates at the wells
		pure subroutine reduce_op(saturation, flux_w, flux_n, prod_w, prod_n)
			real (kind = GRID_SR), intent(in)	    :: saturation
			real (kind = GRID_SR), intent(in)		:: flux_w, flux_n
			real (kind = GRID_SR), intent(inout)	:: prod_w, prod_n

            real (kind = GRID_SR)                   :: lambda_w, lambda_n

            lambda_w = l_w(saturation)
            lambda_n = l_n(saturation)

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
