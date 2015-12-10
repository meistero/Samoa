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

		public initialize_rhs, get_areas_and_lengths

#		define _GT_NAME							t_darcy_permeability_traversal

#		if (_DARCY_FLOW_EDGE_SIZE > 0)
#			define _GT_EDGES
#		endif

#		define _GT_NODES

#		define _GT_ELEMENT_OP					element_op

#		define _GT_NODE_FIRST_TOUCH_OP		    node_first_touch_op
#		define _GT_INNER_NODE_FIRST_TOUCH_OP	inner_node_first_touch_op

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
 			type(t_darcy_permeability_traversal), intent(in)    :: traversal
 			type(t_grid_section), intent(in)					:: section
			type(t_node_data), intent(inout)			        :: node

			call flow_pre_dof_op(node%position(1), node%position(2), node%data_pers%p, node%data_pers%saturation, node%data_pers%rhs, node%data_temp%volume, node%data_temp%is_pressure_dirichlet_boundary, node%data_temp%is_saturation_dirichlet_boundary)
		end subroutine

		elemental subroutine inner_node_first_touch_op(traversal, section, node)
 			type(t_darcy_permeability_traversal), intent(in)    :: traversal
 			type(t_grid_section), intent(in)					:: section
			type(t_node_data), intent(inout)			        :: node

            !ensure the temporary variables are set to .false., the may be invalid otherwise
            node%data_temp%is_saturation_dirichlet_boundary = .false.
            node%data_temp%is_pressure_dirichlet_boundary = .false.

			call flow_pre_dof_op(node%position(1), node%position(2), node%data_pers%p, node%data_pers%saturation, node%data_pers%rhs, node%data_temp%volume, node%data_temp%is_pressure_dirichlet_boundary, node%data_temp%is_saturation_dirichlet_boundary)
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

		elemental subroutine flow_pre_dof_op(pos_x, pos_y, p, saturation, rhs, inflow, is_pressure_dirichlet, is_saturation_dirichlet)
			real (kind = GRID_SR), intent(in)					:: pos_x, pos_y
			real (kind = GRID_SR), intent(inout)				:: p
			real (kind = GRID_SR), intent(inout)				:: saturation
			real (kind = GRID_SR), intent(out)					:: rhs
			real (kind = GRID_SR), intent(out)					:: inflow
            logical, intent(inout)                              :: is_pressure_dirichlet, is_saturation_dirichlet

            integer :: i

			rhs = 0.0_SR
			inflow = 0.0_SR

#           if !defined(_ASAGI)
                if (pos_x == 0.0_SR) then
                    is_saturation_dirichlet = .true.
                else if (pos_x == 1.0_SR) then
                    is_pressure_dirichlet = .true.
                end if
#           endif

            if (is_pressure_dirichlet) then
                p = cfg%r_p_prod
            endif
		end subroutine

		elemental subroutine flow_post_dof_op(saturation, rhs, is_pressure_dirichlet, inflow, total_inflow)
			real (kind = GRID_SR), intent(inout)				:: saturation
			real (kind = GRID_SR), intent(inout)				:: rhs
			real (kind = GRID_SR), intent(in)				    :: inflow
			real (kind = GRID_SR), intent(in)				    :: total_inflow
			logical, intent(in)			                        :: is_pressure_dirichlet

            rhs = rhs + cfg%r_inflow * inflow / total_inflow

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
                real (kind = GRID_SR)					            :: g_local(3), pos_prod(2), pos_in(2), weights(3), edge_length, dx, dy, dz, Ax, Ay, Az, permeability_sum
                real (kind = GRID_SR)                               :: lambda_t(_DARCY_LAYERS, 7)
                real (kind = GRID_SR)                               :: inflow(_DARCY_LAYERS + 1, 3)
                integer                                             :: i, layer
                logical		                                        :: is_dirichlet(_DARCY_LAYERS + 1, 3)

                rhs = 0.0_SR

#               if defined(_ASAGI)
                    call gv_is_saturation_dirichlet%read_from_element(element, is_dirichlet)

                    if (any(is_dirichlet)) then
                        pos_in = samoa_world_to_barycentric_point(element%transform_data, cfg%r_pos_in)

                        !injection well:

                        !The inflow condition is given in um^3 / s
                        !If we devide this by the number of vertical layers, we obtain the 3D inflow for a vertical dual cell column
                        !Split the inflow over all primary cells in each layer that share the dual cell column

#                       define _DARCY_INJ_INFLOW
#                       if defined(_DARCY_INJ_INFLOW)
                            !set a water inflow condition

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
                            inflow = 0.0_SR

                            !split local inflows in primary layer and assign half contributions to dual layers

                            do i = 1, 3
                                inflow(1:_DARCY_LAYERS, i) = inflow(1:_DARCY_LAYERS, i) + 0.5_SR * base_permeability(:, 1) * weights(i)
                                inflow(2:_DARCY_LAYERS + 1, i) = inflow(2:_DARCY_LAYERS + 1, i) + 0.5_SR * base_permeability(:, 1) * weights(i)
                            end do

                            call gv_inflow%add_to_element(element, inflow)
#                       elif defined(_DARCY_INJ_PRESSURE)
                            !set a pressure Dirichlet condition

                            weights = [pos_in(1), 1.0_SR - (pos_in(1) + pos_in(2)), pos_in(2)]

                            do i = 1, 3
                                if (is_dirichlet(_DARCY_LAYERS + 1, i)) then
                                    p(_DARCY_LAYERS + 1, i) = cfg%r_p_in

                                    do layer = _DARCY_LAYERS, 1, -1
                                        p(layer, i) = p(layer + 1, i) - &
                                            (l_w(saturation(layer, i)) * cfg%r_rho_w + l_n(saturation(layer, i)) * cfg%r_rho_n) &
                                            / (l_w(saturation(layer, i)) + l_n(saturation(layer, i))) * cfg%dz * g(3)
                                    end do
                                end if
                            end do

                            call gv_p%write_to_element(element, p)
#                       else
#                           error Injection condition must be defined!
#                       endif
                    end if

                    call gv_is_pressure_dirichlet%read_from_element(element, is_dirichlet)

                    if (any(is_dirichlet)) then
                        !switch to a different discretizaion near the well:
                        !the pressure becomes a logarithmic function p(r) - p(r0) = c log(r/r0),
                        !hence dp/dr(x) = c * r0 / x = (p(r) - p(r0)) / (log(r/r0) * x/r0).
                        !So dp/dr(r/2) = (p(r) - p(r0)) / (log(r/r0) * r/(2 r0))
                        !instead of the finite difference term (p(r) - p(r0)) / (r - r0).

#                       define _DARCY_PROD_ALL_PRESSURE
#                       if defined(_DARCY_PROD_ALL_PRESSURE)
                            !set a constant pressure condition

                            do i = 1, 3
                                if (is_dirichlet(_DARCY_LAYERS + 1, i)) then
                                    p(_DARCY_LAYERS + 1, i) = cfg%r_p_prod

                                    do layer = _DARCY_LAYERS, 1, -1
                                        p(layer, i) = p(layer + 1, i) - &
                                            (l_w(saturation(layer, i)) * cfg%r_rho_w + l_n(saturation(layer, i)) * cfg%r_rho_n) &
                                            / (l_w(saturation(layer, i)) + l_n(saturation(layer, i))) * cfg%dz * cfg%g(3)
                                    end do
                                end if
                            end do

                            call gv_p%write_to_element(element, p)
#                       else
#                           error Production condition must be defined!
#                       endif
                    end if
#               else
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

                !rotate g so it points in the right direction (no scaling!)
                g_local = cfg%g
                g_local(1:2) = samoa_world_to_barycentric_normal(element%transform_data, g_local(1:2))

                call get_areas_and_lengths(element, dx, dy, dz, Ax, Ay, Az)
                call setup_lse_3D(saturation, p, base_permeability, dx, dy, dz, Ax, Ay, Az, g_local, rhs, lambda_t)

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
                real (kind = GRID_SR)					            :: pos_prod(2), pos_in(2), inflow, dx, dy, Ax, Ay
                real (kind = GRID_SR)                               :: g_local(2), lambda_t(2), weights(3)
                logical		                                        :: is_dirichlet(3)

                rhs = 0.0_SR

#               if defined(_ASAGI)
                    call gv_is_saturation_dirichlet%read_from_element(element, is_dirichlet)

                    if (any(is_dirichlet)) then
                        !injection well:

#                       define _DARCY_INJ_INFLOW
#                       if defined(_DARCY_INJ_INFLOW)
                            !set a water inflow condition

                            pos_in = samoa_world_to_barycentric_point(element%transform_data, cfg%r_pos_in)

                            weights = [pos_in(1), 2.0_SR - 2.0_SR * (pos_in(1) + pos_in(2)), pos_in(2)]

                            if (base_permeability > 0) then
                                !The inflow condition is given in um^3 / s
                                !If we devide this by the height of the domain cfg%dz, we obtain the 2D inflow in um^2/s
                                !Split the inflow over all primary cells that share the dual cell

                                inflow = cfg%r_inflow / cfg%dz

                                rhs = rhs + inflow / 8.0_SR * weights
                           end if
#                       elif defined(_DARCY_INJ_PRESSURE)
                            !set a pressure Dirichlet condition

                            where (is_dirichlet)
                                p = cfg%r_p_in
                            end where

                            call gv_p%write_to_element(element, p)
#                       else
#                           error Injection condition must be defined!
#                       endif
                    end if

                    call gv_is_pressure_dirichlet%read_from_element(element, is_dirichlet)

                    if (any(is_dirichlet)) then
                        !production well:

#                       define _DARCY_PROD_ALL_PRESSURE
#                       if defined(_DARCY_PROD_ALL_PRESSURE)
                            !set a pressure Dirichlet condition

                            where (is_dirichlet)
                                p = cfg%r_p_prod
                            end where

                            call gv_p%write_to_element(element, p)
#                       else
#                           error Production condition must be defined!
#                       endif
                    end if
#               else
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

                !rotate g so it points in the right direction (no scaling!)
                g_local = cfg%g(1:2)
                g_local(1:2) = samoa_world_to_barycentric_normal(element%transform_data, g_local(1:2))

                call get_areas_and_lengths(element, dx, dy, Ax, Ay)
                call setup_lse_2D(saturation, p, base_permeability, dx, dy, Ax, Ay, g_local, rhs, lambda_t)

                if (element%transform_data%plotter_data%orientation > 0) then
                    element%cell%data_pers%lambda_t = lambda_t
                else
                    element%cell%data_pers%lambda_t(1) = lambda_t(2)
                    element%cell%data_pers%lambda_t(2) = lambda_t(1)
                end if
            end subroutine
#       endif

#       if (_DARCY_LAYERS > 0)
            subroutine get_areas_and_lengths(element, dx, dy, dz, Ax, Ay, Az)
                type(t_element_base), intent(inout)				:: element
                real (kind = GRID_SR), intent(inout)		    :: dx, dy, dz, Ax, Ay, Az

                logical		                                    :: is_dirichlet(_DARCY_LAYERS + 1, 3)

                dx = element%cell%geometry%get_leg_size()
                dy = dx
                dz = cfg%dz
                Ax = 0.5_SR * dx * dz
                Ay = Ax
                Az = element%cell%geometry%get_volume()

#               if defined(_ASAGI)
                    call gv_is_saturation_dirichlet%read_from_element(element, is_dirichlet)

                    if (any(is_dirichlet)) then
                        !switch to a different discretization near the well:
                        !the pressure becomes a logarithmic function p(r) - p(r0) = c log(r/r0),
                        !hence dp/dr(x) = c * r0 / x = (p(r) - p(r0)) / (log(r/r0) * x/r0).
                        !So dp/dr(r/2) = (p(r) - p(r0)) / (log(r/r0) * r/(2 r0))
                        !instead of the finite difference term (p(r) - p(r0)) / (r - r0).

                        if (is_dirichlet(1, 1) .or. is_dirichlet(1, 2)) then
                            dx = element%cell%geometry%get_leg_size()
                            dx = dx * 2.0_SR / PI * log(dx / cfg%r_well_radius)
                        end if

                        if (is_dirichlet(1, 3) .or. is_dirichlet(1, 2)) then
                            dy = element%cell%geometry%get_leg_size()
                            dy = dy * 2.0_SR / PI * log(dy / cfg%r_well_radius)
                        end if
                    end if

                    call gv_is_pressure_dirichlet%read_from_element(element, is_dirichlet)

                    if (any(is_dirichlet)) then
                        !switch to a different discretization near the well:
                        !the pressure becomes a logarithmic function p(r) - p(r0) = c log(r/r0),
                        !hence dp/dr(x) = c * r0 / x = (p(r) - p(r0)) / (log(r/r0) * x/r0).
                        !So dp/dr(r/2) = (p(r) - p(r0)) / (log(r/r0) * r/(2 r0))
                        !instead of the finite difference term (p(r) - p(r0)) / (r - r0).

                        if (is_dirichlet(1, 1) .or. is_dirichlet(1, 2)) then
                            dx = element%cell%geometry%get_leg_size()
                            dx = dx * 2.0_SR / PI * log(dx / cfg%r_well_radius)
                        end if

                        if (is_dirichlet(1, 3) .or. is_dirichlet(1, 2)) then
                            dy = element%cell%geometry%get_leg_size()
                            dy = dy * 2.0_SR / PI * log(dy / cfg%r_well_radius)
                        end if
                    end if
#               endif
            end subroutine
#       else
            subroutine get_areas_and_lengths(element, dx, dy, Ax, Ay)
                type(t_element_base), intent(inout)				:: element
                real (kind = GRID_SR), intent(inout)		    :: dx, dy, Ax, Ay

                logical		                                    :: is_dirichlet(3)

                dx = element%cell%geometry%get_leg_size()
                dy = dx
                Ax = 0.5_SR * dx
                Ay = Ax

#               if defined(_ASAGI)
                    call gv_is_saturation_dirichlet%read_from_element(element, is_dirichlet)

                    if (any(is_dirichlet)) then
                        !switch to a different discretization near the well:
                        !the pressure becomes a logarithmic function p(r) - p(r0) = c log(r/r0),
                        !hence dp/dr(x) = c * r0 / x = (p(r) - p(r0)) / (log(r/r0) * x/r0).
                        !So dp/dr(r/2) = (p(r) - p(r0)) / (log(r/r0) * r/(2 r0))
                        !instead of the finite difference term (p(r) - p(r0)) / (r - r0).

                        if (is_dirichlet(1) .or. is_dirichlet(2)) then
                            dx = element%cell%geometry%get_leg_size()
                            dx = dx * 2.0_SR / PI * log(dx / cfg%r_well_radius)
                        end if

                        if (is_dirichlet(3) .or. is_dirichlet(2)) then
                            dy = element%cell%geometry%get_leg_size()
                            dy = dy * 2.0_SR / PI * log(dy / cfg%r_well_radius)
                        end if
                    end if

                    call gv_is_pressure_dirichlet%read_from_element(element, is_dirichlet)

                    if (any(is_dirichlet)) then
                        !switch to a different discretization near the well:
                        !the pressure becomes a logarithmic function p(r) - p(r0) = c log(r/r0),
                        !hence dp/dr(x) = c * r0 / x = (p(r) - p(r0)) / (log(r/r0) * x/r0).
                        !So dp/dr(r/2) = (p(r) - p(r0)) / (log(r/r0) * r/(2 r0))
                        !instead of the finite difference term (p(r) - p(r0)) / (r - r0).

                        if (is_dirichlet(1) .or. is_dirichlet(2)) then
                            dx = element%cell%geometry%get_leg_size()
                            dx = dx * 2.0_SR / PI * log(dx / cfg%r_well_radius)
                        end if

                        if (is_dirichlet(3) .or. is_dirichlet(2)) then
                            dy = element%cell%geometry%get_leg_size()
                            dy = dy * 2.0_SR / PI * log(dy / cfg%r_well_radius)
                        end if
                    end if
#               endif
            end subroutine
#       endif

	END MODULE
#endif
