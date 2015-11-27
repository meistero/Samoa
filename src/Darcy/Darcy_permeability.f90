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

		type(darcy_gv_saturation)				        :: gv_saturation
		type(darcy_gv_p)						        :: gv_p
		type(darcy_gv_rhs)						        :: gv_rhs
		type(darcy_gv_is_pressure_dirichlet_boundary)   :: gv_is_pressure_dirichlet
		type(darcy_gv_is_saturation_dirichlet_boundary) :: gv_is_saturation_dirichlet
		type(darcy_gv_volume)                           :: gv_inflow

		public initialize_rhs

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

			call flow_pre_dof_op(node%position(1), node%position(2), node%data_pers%p, node%data_pers%saturation, node%data_pers%rhs, node%data_temp%is_pressure_dirichlet_boundary, node%data_temp%is_saturation_dirichlet_boundary, node%data_temp%volume)
		end subroutine

		subroutine element_op(traversal, section, element)
 			type(t_darcy_permeability_traversal), intent(inout)		:: traversal
 			type(t_grid_section), intent(inout)						:: section
			type(t_element_base), intent(inout)				        :: element

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR)   :: saturation(_DARCY_LAYERS + 1, 3)
                real (kind = GRID_SR)   :: p(_DARCY_LAYERS + 1, 3)
                real (kind = GRID_SR)   :: rhs(_DARCY_LAYERS + 1, 3)
#           else
                real (kind = GRID_SR)   :: saturation(3)
                real (kind = GRID_SR)   :: p(3)
                real (kind = GRID_SR)   :: rhs(3)
#           endif

			call gv_saturation%read_from_element(element, saturation)
			call gv_p%read_from_element(element, p)

			!call element operator
            call initialize_rhs(element, saturation, p, rhs, element%cell%data_pers%base_permeability)

			call gv_rhs%add_to_element(element, rhs)
		end subroutine

		elemental subroutine node_merge_op(local_node, neighbor_node)
			type(t_node_data), intent(inout)		:: local_node
 			type(t_node_data), intent(in)		    :: neighbor_node

            if (any(neighbor_node%data_temp%is_pressure_dirichlet_boundary)) then
                local_node%data_pers%p = neighbor_node%data_pers%p
            end if

			local_node%data_pers%saturation = max(local_node%data_pers%saturation, neighbor_node%data_pers%saturation)
			local_node%data_temp%is_pressure_dirichlet_boundary = local_node%data_temp%is_pressure_dirichlet_boundary .or. neighbor_node%data_temp%is_pressure_dirichlet_boundary
			local_node%data_temp%is_saturation_dirichlet_boundary = local_node%data_temp%is_saturation_dirichlet_boundary .or. neighbor_node%data_temp%is_saturation_dirichlet_boundary
			local_node%data_pers%rhs = local_node%data_pers%rhs + neighbor_node%data_pers%rhs
			local_node%data_temp%volume = local_node%data_temp%volume + neighbor_node%data_temp%volume
		end subroutine

		elemental subroutine node_last_touch_op(traversal, section, node)
 			type(t_darcy_permeability_traversal), intent(in)            :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			                :: node

			real (kind = GRID_SR) :: total_inflow

            total_inflow = tiny(1.0_SR) + sum(node%data_temp%volume)

			call flow_post_dof_op(node%data_pers%saturation, node%data_pers%rhs, node%data_temp%is_pressure_dirichlet_boundary, node%data_temp%volume, total_inflow)
		end subroutine

		elemental subroutine flow_pre_dof_op(pos_x, pos_y, p, saturation, rhs, is_pressure_dirichlet, is_saturation_dirichlet, inflow)
			real (kind = GRID_SR), intent(in)					:: pos_x, pos_y
			real (kind = GRID_SR), intent(inout)				:: p
			real (kind = GRID_SR), intent(inout)				:: saturation
			real (kind = GRID_SR), intent(out)					:: rhs
			logical, intent(out)			                    :: is_pressure_dirichlet, is_saturation_dirichlet
			real (kind = GRID_SR), intent(out)					:: inflow

			rhs = 0.0_SR
			is_saturation_dirichlet = .false.
			is_pressure_dirichlet = .false.
			inflow = 0.0_SR

#           if !defined(_ASAGI)
                if (pos_x == 0.0_SR) then
                    is_saturation_dirichlet = .true.
                    saturation = 1.0_SR
                else if (pos_x == 1.0_SR) then
                    is_pressure_dirichlet = .true.
                    p = cfg%r_p_prod
                end if
#           endif
		end subroutine

		elemental subroutine flow_post_dof_op(saturation, rhs, is_pressure_dirichlet, inflow, total_inflow)
			real (kind = GRID_SR), intent(inout)				:: saturation
			real (kind = GRID_SR), intent(inout)				:: rhs
			real (kind = GRID_SR), intent(in)				    :: inflow
			real (kind = GRID_SR), intent(in)				    :: total_inflow
			logical, intent(in)			                        :: is_pressure_dirichlet

            rhs = rhs + cfg%r_inflow * inflow / total_inflow

            !limit the accumulated saturation at the injection well
            saturation = min(saturation, 1.0_SR)

            if (is_pressure_dirichlet) then
                rhs = 0.0_SR
            end if
		end subroutine

#       if (_DARCY_LAYERS > 0)
            subroutine initialize_rhs(element, saturation, p, rhs, base_permeability)
                type(t_element_base), intent(inout)				                    :: element
                real (kind = GRID_SR), intent(inout)		                        :: saturation(:, :)
                real (kind = GRID_SR), intent(inout)			                    :: p(:, :)
                real (kind = GRID_SR), intent(out)									:: rhs(:, :)
                real (kind = GRID_SR), intent(inout)                                :: base_permeability(:, :)

                real (kind = GRID_SR)					            :: coords(2, 3)
                real (kind = GRID_SR)					            :: g_local(3), pos_prod(2), pos_in(2), radius, weights(3), edge_length, surface, dz, permeability_sum
                real (kind = GRID_SR)                               :: lambda_w(_DARCY_LAYERS + 1, 3), lambda_n(_DARCY_LAYERS + 1, 3), lambda_t(_DARCY_LAYERS, 7)
                real (kind = GRID_SR)                               :: inflow(_DARCY_LAYERS + 1, 3)
                integer                                             :: i, layer
                logical		                                        :: is_dirichlet(_DARCY_LAYERS + 1, 3)

                rhs = 0.0_SR

                radius = cfg%r_well_radius / element%transform_data%custom_data%scaling

                pos_in = samoa_world_to_barycentric_point(element%transform_data, cfg%r_pos_in)

                lambda_w = l_w(saturation)
                lambda_n = l_n(saturation)

                if (norm2(pos_in) < 1.0_SR + epsilon(1.0_SR)) then
                    if (pos_in(1) > -epsilon(1.0_SR) .and. 1.0_SR - (pos_in(1) + pos_in(2)) > -epsilon(1.0_SR) .and. pos_in(2) > -epsilon(1.0_SR)) then
                        !injection well:
                        !set an inflow pressure condition and a constant saturation condition

                        !assume the whole well is filled with water
                        weights = [pos_in(1), 1.0_SR - (pos_in(1) + pos_in(2)), pos_in(2)]
                        saturation(:, 1) = max(0.0_SR, min(1.0_SR, saturation(:, 1) + weights(1)))
                        saturation(:, 2) = max(0.0_SR, min(1.0_SR, saturation(:, 2) + weights(2)))
                        saturation(:, 3) = max(0.0_SR, min(1.0_SR, saturation(:, 3) + weights(3)))
                        call gv_saturation%add_to_element(element, spread(weights, 1, _DARCY_LAYERS + 1))

                        lambda_w = l_w(saturation)
                        lambda_n = l_n(saturation)

                        !The inflow condition is given in um^3 / s
                        !If we devide this by the number of vertical layers, we obtain the 3D inflow for a vertical dual cell column
                        !Split the inflow over all primary cells in each layer that share the dual cell column

#                       define _DARCY_INJ_INFLOW
#                       if defined(_DARCY_INJ_INFLOW)
                            !Using Peaceman's well model we consider the well as an internal boundary and assume that
                            !near the well the following condition holds:
                            !The radial derivative p_r is constant over r and z.
                            !
                            !Thus the inflow q is q(r,phi,z) = lambda_w(S) K_r(r,phi,z) (-p_r)
                            !With \integral_{well boundary} q(r,phi,z) * dS = Q we obtain
                            ! \integral_{well boundary} lambda_w(S) K_r(r,phi,z) (-p_r) dOmega = Q
                            ! p_r = - Q / integral_(well boundary) lambda_w(S) K_r(r,phi,z) dOmega
                            ! q(r,phi,z) = Q K_r(r,phi,z) / integral_(well boundary) K_r(r,phi,z) dOmega

                            weights = 0.25_SR * [pos_in(1), 2.0_SR - 2.0_SR * (pos_in(1) + pos_in(2)), pos_in(2)]

                            is_dirichlet(:, 1) = weights(1) > epsilon(1.0_SR)
                            is_dirichlet(:, 2) = weights(2) > epsilon(1.0_SR)
                            is_dirichlet(:, 3) = weights(3) > epsilon(1.0_SR)

                            inflow = 0.0_SR

                            !split local inflows in primary layer and assign half contributions to dual layers

                            do i = 1, 3
                                inflow(1:_DARCY_LAYERS, i) = inflow(1:_DARCY_LAYERS, i) + 0.5_SR * base_permeability(:, 1) * weights(i)
                                inflow(2:_DARCY_LAYERS + 1, i) = inflow(2:_DARCY_LAYERS + 1, i) + 0.5_SR * base_permeability(:, 1) * weights(i)
                            end do

                            call gv_inflow%add_to_element(element, inflow)
                            call gv_is_saturation_dirichlet%add_to_element(element, is_dirichlet)
#                       elif defined(_DARCY_INJ_TOP_INFLOW)
                            weights = [pos_in(1), 2.0_SR - 2.0_SR * (pos_in(1) + pos_in(2)), pos_in(2)]

                            is_dirichlet = .false.
                            is_dirichlet(_DARCY_LAYERS + 1, 1) = weights(1) > epsilon(1.0_SR)
                            is_dirichlet(_DARCY_LAYERS + 1, 2) = weights(2) > epsilon(1.0_SR)
                            is_dirichlet(_DARCY_LAYERS + 1, 3) = weights(3) > epsilon(1.0_SR)

                            inflow = cfg%r_inflow / cfg%dz

                            !For this inflow condition, we need an asymmetric element matrix where contributions from neighbor dual cells to the well cells are 0,
                            !but contributions from the well cells to neighbor dual cells are nonzero.

                            !base_permeability(:, 1) = 0.0_SR
                            !base_permeability(:, 2) = 1.0_SR

                            !split local inflows and assign contributions to top layer
                            rhs(_DARCY_LAYERS + 1, :) = rhs(_DARCY_LAYERS + 1, :) + inflow / 8.0_SR * weights

                            call gv_is_saturation_dirichlet%add_to_element(element, is_dirichlet)
#                       elif defined(_DARCY_INJ_PRESSURE)
                            weights = [pos_in(1), 1.0_SR - (pos_in(1) + pos_in(2)), pos_in(2)]

                            is_dirichlet(:, 1) = weights(1) > epsilon(1.0_SR)
                            is_dirichlet(:, 2) = weights(2) > epsilon(1.0_SR)
                            is_dirichlet(:, 3) = weights(3) > epsilon(1.0_SR)

                            do i = 1, 3
                                if (is_dirichlet(_DARCY_LAYERS + 1, i)) then
                                    p(_DARCY_LAYERS + 1, i) = cfg%r_p_in

                                    do layer = _DARCY_LAYERS, 1, -1
                                        p(layer, i) = p(layer + 1, i) - &
                                            (lambda_w(layer, i) * cfg%r_rho_w + lambda_n(layer, i) * cfg%r_rho_n) &
                                            / (lambda_w(layer, i) + lambda_n(layer, i)) * cfg%dz * g(3)
                                    end do
                                end if
                            end do

                            call gv_p%write_to_element(element, p)
                            call gv_is_pressure_dirichlet%add_to_element(element, is_dirichlet)
#                       else
#                           error Injection condition must be defined!
#                       endif
                    end if
                end if

                pos_prod = sign(cfg%r_pos_prod - 0.5_SR, element%transform_data%custom_data%offset - 0.5_SR) + 0.5_SR
                pos_prod = samoa_world_to_barycentric_point(element%transform_data, pos_prod)

                if (norm2(pos_prod) < 1.0_SR + epsilon(1.0_SR)) then
                    if (pos_prod(1) > -epsilon(1.0_SR) .and. 1.0_SR - (pos_prod(1) + pos_prod(2)) > -epsilon(1.0_SR) .and. pos_prod(2) > -epsilon(1.0_SR)) then
                        !production well:
                        !set a constant pressure condition and an outflow saturation condition

                        weights = [pos_prod(1), 1.0_SR - (pos_prod(1) + pos_prod(2)), pos_prod(2)]

#                       define _DARCY_PROD_ALL_PRESSURE
#                       if defined(_DARCY_PROD_ALL_PRESSURE)
                            is_dirichlet(:, 1) = weights(1) > epsilon(1.0_SR)
                            is_dirichlet(:, 2) = weights(2) > epsilon(1.0_SR)
                            is_dirichlet(:, 3) = weights(3) > epsilon(1.0_SR)

                            do i = 1, 3
                                if (is_dirichlet(_DARCY_LAYERS + 1, i)) then
                                    p(_DARCY_LAYERS + 1, i) = cfg%r_p_prod

                                    do layer = _DARCY_LAYERS, 1, -1
                                        p(layer, i) = p(layer + 1, i) - &
                                            (lambda_w(layer, i) * cfg%r_rho_w + lambda_n(layer, i) * cfg%r_rho_n) &
                                            / (lambda_w(layer, i) + lambda_n(layer, i)) * cfg%dz * cfg%g(3)
                                    end do
                                end if
                            end do
#                       else
#                           error Production condition must be defined!
#                       endif

                        call gv_p%write_to_element(element, p)
                        call gv_is_pressure_dirichlet%add_to_element(element, is_dirichlet)
                    end if
                end if

                !rotate g so it points in the right direction (no scaling!)
                g_local = cfg%g
                g_local(1:2) = samoa_world_to_barycentric_normal(element%transform_data, g_local(1:2))
                g_local(1:2) = g_local(1:2) / (element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))

                edge_length = element%cell%geometry%get_leg_size()
                surface = element%cell%geometry%get_volume()
                dz = cfg%dz

                do i = 1, _DARCY_LAYERS
                    call compute_rhs_1D(edge_length, 0.25_SR * edge_length * dz, base_permeability(i, 1), p(i, 2), p(i, 1), lambda_w(i, 2), lambda_w(i, 1), lambda_n(i, 2), lambda_n(i, 1), g_local(1), rhs(i, 2), rhs(i, 1), lambda_t(i, 1))
                    call compute_rhs_1D(edge_length, 0.25_SR * edge_length * dz, base_permeability(i, 1), p(i, 2), p(i, 3), lambda_w(i, 2), lambda_w(i, 3), lambda_n(i, 2), lambda_n(i, 3), g_local(2), rhs(i, 2), rhs(i, 3), lambda_t(i, 2))

                    call compute_rhs_1D(dz, 0.25_SR * surface, base_permeability(i, 2), p(i, 1), p(i + 1, 1), lambda_w(i, 1), lambda_w(i + 1, 1), lambda_n(i, 1), lambda_n(i + 1, 1), g_local(3), rhs(i, 1), rhs(i + 1, 1), lambda_t(i, 3))
                    call compute_rhs_1D(dz, 0.50_SR * surface, base_permeability(i, 2), p(i, 2), p(i + 1, 2), lambda_w(i, 2), lambda_w(i + 1, 2), lambda_n(i, 2), lambda_n(i + 1, 2), g_local(3), rhs(i, 2), rhs(i + 1, 2), lambda_t(i, 4))
                    call compute_rhs_1D(dz, 0.25_SR * surface, base_permeability(i, 2), p(i, 3), p(i + 1, 3), lambda_w(i, 3), lambda_w(i + 1, 3), lambda_n(i, 3), lambda_n(i + 1, 3), g_local(3), rhs(i, 3), rhs(i + 1, 3), lambda_t(i, 5))

                    call compute_rhs_1D(edge_length, 0.25_SR * edge_length * dz, base_permeability(i, 1), p(i + 1, 2), p(i + 1, 1), lambda_w(i + 1, 2), lambda_w(i + 1, 1), lambda_n(i + 1, 2), lambda_n(i + 1, 1), g_local(1), rhs(i + 1, 2), rhs(i + 1, 1), lambda_t(i, 6))
                    call compute_rhs_1D(edge_length, 0.25_SR * edge_length * dz, base_permeability(i, 1), p(i + 1, 2), p(i + 1, 3), lambda_w(i + 1, 2), lambda_w(i + 1, 3), lambda_n(i + 1, 2), lambda_n(i + 1, 3), g_local(2), rhs(i + 1, 2), rhs(i + 1, 3), lambda_t(i, 7))
                end do

#               if !defined(_ASAGI)
                    coords(:, 1) = samoa_barycentric_to_world_point(element%transform_data, [1.0_SR, 0.0_SR])
                    coords(:, 2) = samoa_barycentric_to_world_point(element%transform_data, [0.0_SR, 0.0_SR])
                    coords(:, 3) = samoa_barycentric_to_world_point(element%transform_data, [0.0_SR, 1.0_SR])

                    if (coords(1, 1) + coords(1, 2) < epsilon(1.0_SR)) then
                        rhs(:, 1) = rhs(:, 1) + cfg%r_inflow / real(_DARCY_LAYERS + 1, SR) * 0.5_SR * element%cell%geometry%get_leg_size()
                        rhs(:, 2) = rhs(:, 2) + cfg%r_inflow / real(_DARCY_LAYERS + 1, SR) * 0.5_SR * element%cell%geometry%get_leg_size()
                    else if (coords(1, 1) + coords(1, 3) < epsilon(1.0_SR)) then
                        rhs(:, 1) = rhs(:, 1) + cfg%r_inflow / real(_DARCY_LAYERS + 1, SR) * 0.5_SR * element%cell%geometry%get_hypo_size()
                        rhs(:, 3) = rhs(:, 3) + cfg%r_inflow / real(_DARCY_LAYERS + 1, SR) * 0.5_SR * element%cell%geometry%get_hypo_size()
                    else if (coords(1, 2) + coords(1, 3) < epsilon(1.0_SR)) then
                        rhs(:, 2) = rhs(:, 2) + cfg%r_inflow / real(_DARCY_LAYERS + 1, SR) * 0.5_SR * element%cell%geometry%get_leg_size()
                        rhs(:, 3) = rhs(:, 3) + cfg%r_inflow / real(_DARCY_LAYERS + 1, SR) * 0.5_SR * element%cell%geometry%get_leg_size()
                    end if
#               endif

                if (element%transform_data%plotter_data%orientation > 0) then
                    element%cell%data_pers%lambda_t = lambda_t
                else
                    element%cell%data_pers%lambda_t(:, 1) = lambda_t(:, 2)
                    element%cell%data_pers%lambda_t(:, 2) = lambda_t(:, 1)
                    element%cell%data_pers%lambda_t(:, 3) = lambda_t(:, 5)
                    element%cell%data_pers%lambda_t(:, 4) = lambda_t(:, 4)
                    element%cell%data_pers%lambda_t(:, 5) = lambda_t(:, 3)
                    element%cell%data_pers%lambda_t(:, 6) = lambda_t(:, 7)
                    element%cell%data_pers%lambda_t(:, 7) = lambda_t(:, 6)
                end if
            end subroutine
#       else
            subroutine initialize_rhs(element, saturation, p, rhs, base_permeability)
                type(t_element_base), intent(inout)				                    :: element
                real (kind = GRID_SR), intent(inout)		                        :: saturation(:)
                real (kind = GRID_SR), intent(inout)			                    :: p(:)
                real (kind = GRID_SR), intent(out)									:: rhs(:)
                real (kind = GRID_SR), intent(in)                                   :: base_permeability

                real (kind = GRID_SR)					            :: coords(2, 3)
                real (kind = GRID_SR)					            :: pos_prod(2), pos_in(2), radius, inflow, edge_length
                real (kind = GRID_SR)                               :: g_local(2), lambda_w(3), lambda_n(3), lambda_t(2), weights(3)
                logical		                                        :: is_dirichlet(3)

                rhs = 0.0_SR

                radius = cfg%r_well_radius / element%transform_data%custom_data%scaling

                pos_in = samoa_world_to_barycentric_point(element%transform_data, cfg%r_pos_in)

                if (norm2(pos_in) < 1.0_SR + epsilon(1.0_SR)) then
                    if (norm2(pos_in - 0.5_SR) < sqrt(0.5_SR) + epsilon(1.0_SR)) then
                        !injection well:
                        !set an inflow pressure condition and a constant saturation condition

                        saturation = max(0.0_SR, min(1.0_SR, saturation + [pos_in(1), 1.0_SR - (pos_in(1) + pos_in(2)), pos_in(2)]))
                        call gv_saturation%add_to_element(element, [pos_in(1), 1.0_SR - (pos_in(1) + pos_in(2)), pos_in(2)])

                        weights = [pos_in(1), 2.0_SR - 2.0_SR * (pos_in(1) + pos_in(2)), pos_in(2)]
                        is_dirichlet = weights > epsilon(1.0_SR)

                        call gv_is_saturation_dirichlet%add_to_element(element, is_dirichlet)

                        if (base_permeability > 0) then
                            !The inflow condition is given in um^3 / s
                            !If we devide this by the height of the domain cfg%dz, we obtain the 2D inflow in um^2/s
                            !Split the inflow over all primary cells that share the dual cell

                            inflow = cfg%r_inflow / cfg%dz

                            rhs = rhs + inflow / 8.0_SR * weights
                        end if
                    end if
                end if

                pos_prod = sign(cfg%r_pos_prod - 0.5_SR, element%transform_data%custom_data%offset - 0.5_SR) + 0.5_SR
                pos_prod = samoa_world_to_barycentric_point(element%transform_data, pos_prod)

                if (norm2(pos_prod) < 1.0_SR + epsilon(1.0_SR)) then
                    if (pos_prod(1) > -epsilon(1.0_SR) .and. 1.0_SR - (pos_prod(1) + pos_prod(2)) > -epsilon(1.0_SR) .and. pos_prod(2) > -epsilon(1.0_SR)) then
                        !production well:
                        !set a constant pressure condition and an outflow saturation condition

                        is_dirichlet = [pos_prod(1) > epsilon(1.0_SR), 1.0_SR - (pos_prod(1) + pos_prod(2)) > epsilon(1.0_SR), pos_prod(2) > epsilon(1.0_SR)]

                        where (is_dirichlet)
                            p = cfg%r_p_prod
                        end where

                        call gv_p%write_to_element(element, p)
                        call gv_is_pressure_dirichlet%add_to_element(element, is_dirichlet)

                        rhs = 0.0_SR
                    end if
                end if

                lambda_w = l_w(saturation)
                lambda_n = l_n(saturation)

                !rotate g so it points in the right direction (no scaling!)
                g_local = cfg%g(1:2)
                g_local(1:2) = samoa_world_to_barycentric_normal(element%transform_data, g_local(1:2))
                g_local(1:2) = g_local(1:2) / (element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))

                edge_length = element%cell%geometry%get_leg_size()

                call compute_rhs_1D(edge_length, 0.5_SR * edge_length, base_permeability, p(2), p(1), lambda_w(2), lambda_w(1), lambda_n(2), lambda_n(1), g_local(1), rhs(2), rhs(1), lambda_t(1))
                call compute_rhs_1D(edge_length, 0.5_SR * edge_length, base_permeability, p(2), p(3), lambda_w(2), lambda_w(3), lambda_n(2), lambda_n(3), g_local(2), rhs(2), rhs(3), lambda_t(2))

#               if !defined(_ASAGI)
                    coords(:, 1) = samoa_barycentric_to_world_point(element%transform_data, [1.0_SR, 0.0_SR])
                    coords(:, 2) = samoa_barycentric_to_world_point(element%transform_data, [0.0_SR, 0.0_SR])
                    coords(:, 3) = samoa_barycentric_to_world_point(element%transform_data, [0.0_SR, 1.0_SR])

                    if (coords(1, 1) + coords(1, 2) < epsilon(1.0_SR)) then
                        rhs(1) = rhs(1) + cfg%r_inflow * 0.5_SR * element%cell%geometry%get_leg_size()
                        rhs(2) = rhs(2) + cfg%r_inflow * 0.5_SR * element%cell%geometry%get_leg_size()
                    else if (coords(1, 1) + coords(1, 3) < epsilon(1.0_SR)) then
                        rhs(1) = rhs(1) + cfg%r_inflow * 0.5_SR * element%cell%geometry%get_hypo_size()
                        rhs(3) = rhs(3) + cfg%r_inflow * 0.5_SR * element%cell%geometry%get_hypo_size()
                    else if (coords(1, 2) + coords(1, 3) < epsilon(1.0_SR)) then
                        rhs(2) = rhs(2) + cfg%r_inflow * 0.5_SR * element%cell%geometry%get_leg_size()
                        rhs(3) = rhs(3) + cfg%r_inflow * 0.5_SR * element%cell%geometry%get_leg_size()
                    end if
#               endif

                if (element%transform_data%plotter_data%orientation > 0) then
                    element%cell%data_pers%lambda_t = lambda_t
                else
                    element%cell%data_pers%lambda_t(1) = lambda_t(2)
                    element%cell%data_pers%lambda_t(2) = lambda_t(1)
                end if
            end subroutine
#       endif
	END MODULE
#endif
