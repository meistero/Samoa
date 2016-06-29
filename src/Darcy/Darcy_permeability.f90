! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_permeability
		use SFC_edge_traversal
		use Samoa_darcy

        type num_traversal_data
            logical                                     :: is_matrix_modified
        end type

		type(darcy_gv_saturation)				        :: gv_saturation
		type(darcy_gv_p)						        :: gv_p
		type(darcy_gv_r)						        :: gv_r
		type(darcy_gv_rhs)						        :: gv_rhs
		type(darcy_gv_volume)                           :: gv_trace
		type(darcy_gv_volume)                           :: gv_inflow
		type(darcy_gv_boundary_condition)               :: gv_boundary_condition
		type(darcy_gm_A)                                :: gm_A

		public initialize_rhs, get_areas_and_lengths

#		define _GT_NAME							t_darcy_permeability_traversal

#		if (_DARCY_FLOW_EDGE_SIZE > 0)
#			define _GT_EDGES
#		endif

#		define _GT_NODES

#		define _GT_PRE_TRAVERSAL_OP		        pre_traversal_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op

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

		subroutine pre_traversal_op(traversal, section)
			type(t_darcy_permeability_traversal), intent(inout)		:: traversal
			type(t_grid_section), intent(inout)							    :: section

            traversal%is_matrix_modified = .false.
		end subroutine

		subroutine post_traversal_grid_op(traversal, grid)
			type(t_darcy_permeability_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							    :: grid

            call reduce(traversal%is_matrix_modified, traversal%sections%is_matrix_modified, MPI_LOR, .true.)

            !also reduce the injector pressure as it is used for the linear solver exit criterion
            call reduce(grid%p_bh, MPI_MAX)
		end subroutine

		!*******************************
		!Geometry operators
		!*******************************

		elemental subroutine node_first_touch_op(traversal, section, node)
 			type(t_darcy_permeability_traversal), intent(in)    :: traversal
 			type(t_grid_section), intent(in)					:: section
			type(t_node_data), intent(inout)			        :: node

			call flow_pre_dof_op(node%position(1), node%position(2), node%data_pers%p, node%data_pers%saturation, node%data_pers%r, node%data_pers%rhs, node%data_temp%volume)
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
            call initialize_rhs(element, saturation, p, rhs, element%cell%data_pers%base_permeability, traversal%is_matrix_modified)

			call gv_rhs%add_to_element(element, rhs)
		end subroutine

		elemental subroutine node_merge_op(local_node, neighbor_node)
			type(t_node_data), intent(inout)		:: local_node
 			type(t_node_data), intent(in)		    :: neighbor_node

			local_node%data_pers%r = local_node%data_pers%r + neighbor_node%data_pers%r
			local_node%data_pers%rhs = local_node%data_pers%rhs + neighbor_node%data_pers%rhs
			local_node%data_temp%volume = local_node%data_temp%volume + neighbor_node%data_temp%volume
		end subroutine

		elemental subroutine node_last_touch_op(traversal, section, node)
 			type(t_darcy_permeability_traversal), intent(in)            :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			                :: node

#           if defined(_ASAGI)
#               if defined(_DARCY_INJ_INFLOW)
                    if (any(node%data_pers%boundary_condition > 0)) then
                        call inflow_post_dof_op(node%data_pers%rhs, node%data_pers%r, node%data_temp%volume)
                    end if
#               elif defined(_DARCY_INJ_PRESSURE)
                    if (any(node%data_pers%boundary_condition > 0)) then
                        call pressure_post_dof_op(node%data_pers%rhs)
                    end if
#               endif
#           endif
		end subroutine

		elemental subroutine flow_pre_dof_op(pos_x, pos_y, p, saturation, r, rhs, inflow)
			real (kind = GRID_SR), intent(in)					:: pos_x, pos_y
			real (kind = GRID_SR), intent(inout)				:: p, saturation
			real (kind = GRID_SR), intent(out)					:: r, rhs, inflow

			r = 0.0_SR
			rhs = 0.0_SR
			inflow = 0.0_SR
		end subroutine

		pure subroutine inflow_post_dof_op(rhs, r, d)
			real (kind = GRID_SR), intent(inout)				:: rhs(:)
			real (kind = GRID_SR), intent(in)				    :: r(:), d(:)

            rhs = rhs + cfg%r_inflow / sum(d) * d
		end subroutine

		pure subroutine pressure_post_dof_op(rhs)
			real (kind = GRID_SR), intent(inout)				:: rhs(:)

            rhs = (sum(rhs) + cfg%r_inflow) / (_DARCY_LAYERS + 1)
		end subroutine

#       if (_DARCY_LAYERS > 0)
            subroutine initialize_rhs(element, saturation, p, rhs, base_permeability, is_matrix_modified)
                type(t_element_base), intent(inout)				                    :: element
                real (kind = GRID_SR), intent(inout)		                        :: saturation(:, :)
                real (kind = GRID_SR), intent(inout)			                    :: p(:, :)
                real (kind = GRID_SR), intent(out)									:: rhs(:, :)
                real (kind = GRID_SR), intent(inout)                                :: base_permeability(:, :)
                logical, intent(inout)                                              :: is_matrix_modified

                real (kind = GRID_SR)					            :: coords(2, 3)
                real (kind = GRID_SR)					            :: g_local(3), weights(3), edge_length, dx, dy, dz, Ax, Ay, Az
                real (kind = GRID_SR)                               :: lambda_t(_DARCY_LAYERS, 7)
                real (kind = GRID_SR)                               :: inflow((_DARCY_LAYERS + 1) * 3)
                integer                                             :: i, layer
                integer (kind = SI) 	                            :: boundary_condition(3)

                rhs = 0.0_SR

                !rotate g so it points in the right direction (no scaling!)
                g_local = cfg%g
                g_local(1:2) = samoa_world_to_barycentric_normal(element%transform_data, g_local(1:2))

                call get_areas_and_lengths(element, dx, dy, dz, Ax, Ay, Az)
                call setup_lse_3D(saturation, p, base_permeability, dx, dy, dz, Ax, Ay, Az, g_local, rhs, lambda_t)

#               if defined(_ASAGI)
                    call gv_boundary_condition%read_from_element(element, boundary_condition)

                    if (any(boundary_condition > 0)) then
                        !injection well:

#                       if defined(_DARCY_INJ_INFLOW)
                            !set a water inflow condition

                            inflow = 0.0_SR

                            call gm_A%get_trace(element, inflow)
                            call gv_trace%add_to_element(element, inflow)
#                       elif defined(_DARCY_INJ_PRESSURE)
                            do i = 1, 3
                                if (boundary_condition(i) > 0) then
                                    !the bottom hole pressure is determined by the linear solver,
                                    !so just make sure that the well pressure is analytical.

                                    do layer = 1, _DARCY_LAYERS
                                        p(layer, i) = p(_DARCY_LAYERS + 1, i) - cfg%r_rho_w * (_DARCY_LAYERS + 1 - layer) * cfg%dz * cfg%g(3)
                                    end do

                                    !The well is a Dirichlet condition (of sorts), so set
                                    !local vertical contributions to 0.

                                    lambda_t(:, 2 + i) = 0.0_SR
                                end if
                            end do

                            call gv_p%write_to_element(element, p)

#                       else
#                           error Injection condition must be defined!
#                       endif
                    end if

                    if (any(boundary_condition < 0)) then
#                       if defined(_DARCY_PROD_PRESSURE)
                            !set a constant pressure condition

                            do i = 1, 3
                                if (boundary_condition(i) < 0) then
                                    p(_DARCY_LAYERS + 1, i) = cfg%r_p_prod

                                    do layer = _DARCY_LAYERS, 1, -1
                                        p(layer, i) = p(layer + 1, i) - &
                                            (l_w(saturation(layer, i)) * cfg%r_rho_w + l_n(saturation(layer, i)) * cfg%r_rho_n) &
                                            / (l_w(saturation(layer, i)) + l_n(saturation(layer, i))) * cfg%dz * cfg%g(3)
                                    end do

                                    !Due to the Dirichlet condition the local vertical contributions are ignored,
                                    !we do not have to set them to 0 explicitly.
                                end if
                            end do

                            call gv_p%write_to_element(element, p)
#                       else
#                           error Production condition must be defined!
#                       endif
                    end if
#               else
                    call gv_boundary_condition%read_from_element(element, boundary_condition)

                    if (any(boundary_condition > 0)) then
                        if (boundary_condition(1) > 0 .and. boundary_condition(2) > 0) then
                            rhs(:, 1) = rhs(:, 1) + cfg%r_inflow / real(_DARCY_LAYERS + 1, SR) * 0.5_SR * element%cell%geometry%get_leg_size()
                            rhs(:, 2) = rhs(:, 2) + cfg%r_inflow / real(_DARCY_LAYERS + 1, SR) * 0.5_SR * element%cell%geometry%get_leg_size()
                        else if (boundary_condition(1) > 0 .and. boundary_condition(3) > 0) then
                            rhs(:, 1) = rhs(:, 1) + cfg%r_inflow / real(_DARCY_LAYERS + 1, SR) * 0.5_SR * element%cell%geometry%get_hypo_size()
                            rhs(:, 3) = rhs(:, 3) + cfg%r_inflow / real(_DARCY_LAYERS + 1, SR) * 0.5_SR * element%cell%geometry%get_hypo_size()
                        else if (boundary_condition(2) > 0 .and. boundary_condition(3) > 0) then
                            rhs(:, 2) = rhs(:, 2) + cfg%r_inflow / real(_DARCY_LAYERS + 1, SR) * 0.5_SR * element%cell%geometry%get_leg_size()
                            rhs(:, 3) = rhs(:, 3) + cfg%r_inflow / real(_DARCY_LAYERS + 1, SR) * 0.5_SR * element%cell%geometry%get_leg_size()
                        end if
                    end if

                    if (any(boundary_condition < 0)) then
                        !set a constant pressure condition

                        do i = 1, 3
                            if (boundary_condition(i) < 0) then
                                p(:, i) = cfg%r_p_prod
                            end if
                        end do

                        call gv_p%write_to_element(element, p)
                    end if
#               endif

                if (element%transform_data%plotter_data%orientation > 0) then
                    is_matrix_modified = is_matrix_modified .or. any(element%cell%data_pers%lambda_t .ne. lambda_t)

                    element%cell%data_pers%lambda_t = lambda_t
                else
                    is_matrix_modified = is_matrix_modified .or. any(element%cell%data_pers%lambda_t .ne. lambda_t(:, [2, 1, 5, 4, 3, 7, 6]))

                    element%cell%data_pers%lambda_t = lambda_t(:, [2, 1, 5, 4, 3, 7, 6])
                end if
            end subroutine
#       else
            subroutine initialize_rhs(element, saturation, p, rhs, base_permeability, is_matrix_modified)
                type(t_element_base), intent(inout)				                    :: element
                real (kind = GRID_SR), intent(inout)		                        :: saturation(:)
                real (kind = GRID_SR), intent(inout)			                    :: p(:)
                real (kind = GRID_SR), intent(out)									:: rhs(:)
                real (kind = GRID_SR), intent(in)                                   :: base_permeability
                logical, intent(inout)                                              :: is_matrix_modified

                real (kind = GRID_SR)					            :: coords(2, 3)
                real (kind = GRID_SR)					            :: inflow(3), r(3), dx, dy, Ax, Ay
                real (kind = GRID_SR)                               :: g_local(2), lambda_t(2), weights(3)
                integer                                             :: i
                integer (kind = SI)	                                :: boundary_condition(3)

                rhs = 0.0_SR

                !rotate g so it points in the right direction (no scaling!)
                g_local = cfg%g(1:2)
                g_local(1:2) = samoa_world_to_barycentric_normal(element%transform_data, g_local(1:2))

                call get_areas_and_lengths(element, dx, dy, Ax, Ay)
                call setup_lse_2D(saturation, p, base_permeability, dx, dy, Ax, Ay, g_local, rhs, lambda_t)

#               if defined(_ASAGI)
                    call gv_boundary_condition%read_from_element(element, boundary_condition)

                    if (any(boundary_condition > 0)) then
                        !injection well:

                        !in 2D, the two injection conditions should give the same output

#                       if defined(_DARCY_INJ_INFLOW)
                            !set an inflow injection condition

                            inflow = 0.0_SR

                            call gm_A%get_trace(element, inflow)
                            call gv_trace%add_to_element(element, inflow)
#                       elif defined(_DARCY_INJ_PRESSURE)
                            !do nothing, the pressure should be evaluated on its own
#                       else
#                           error Injection condition must be defined!
#                       endif
                    end if

                    if (any(boundary_condition < 0)) then
                        !production well:

#                       if defined(_DARCY_PROD_PRESSURE)
                            !set an internal pressure Dirichlet condition

                            where (boundary_condition < 0)
                                p = cfg%r_p_prod
                            end where

                            call gv_p%write_to_element(element, p)
#                       else
#                           error Production condition must be defined!
#                       endif
                    end if
#               else
                    call gv_boundary_condition%read_from_element(element, boundary_condition)

                    if (any(boundary_condition > 0)) then
                        if (boundary_condition(1) > 0 .and. boundary_condition(2) > 0) then
                            rhs(1) = rhs(1) + cfg%r_inflow * 0.5_SR * element%cell%geometry%get_leg_size()
                            rhs(2) = rhs(2) + cfg%r_inflow * 0.5_SR * element%cell%geometry%get_leg_size()
                        else if (boundary_condition(1) > 0 .and. boundary_condition(3) > 0) then
                            rhs(1) = rhs(1) + cfg%r_inflow * 0.5_SR * element%cell%geometry%get_hypo_size()
                            rhs(3) = rhs(3) + cfg%r_inflow * 0.5_SR * element%cell%geometry%get_hypo_size()
                        else if (boundary_condition(2) > 0 .and. boundary_condition(3) > 0) then
                            rhs(2) = rhs(2) + cfg%r_inflow * 0.5_SR * element%cell%geometry%get_leg_size()
                            rhs(3) = rhs(3) + cfg%r_inflow * 0.5_SR * element%cell%geometry%get_leg_size()
                        end if
                    end if

                    if (any(boundary_condition < 0)) then
                        !set a constant pressure condition

                        where (boundary_condition < 0)
                            p = cfg%r_p_prod
                        end where

                        call gv_p%write_to_element(element, p)
                    end if
#               endif

                !check if any values of lambda_t changed, indicating that the linear system must be solved again

                if (element%transform_data%plotter_data%orientation > 0) then
                    is_matrix_modified = is_matrix_modified .or. any(element%cell%data_pers%lambda_t .ne. lambda_t)

                    element%cell%data_pers%lambda_t = lambda_t
                else
                    is_matrix_modified = is_matrix_modified .or. any(element%cell%data_pers%lambda_t .ne. lambda_t([2, 1]))

                    element%cell%data_pers%lambda_t = lambda_t([2, 1])
                end if
            end subroutine
#       endif

        !According to the Peaceman well model, the pressure error is approximately
        ! err = Q / (2 pi lambda k) * (log(dx/r_w) + log(0.2)).

        !Using the pressure correction p'_w := p_w - err
        !we get (p(1) - p_w) / dx = (p(1) - p'_w) / dx +
        ! Q / (2 pi lambda k) * (log(dx/r_w) - gamma - log(4))) / dx
        !At the well the discretization gives -(p(1) - p'_w)/dx * dx/2 = Q / (8 lambda k), so
        ! (p(1) - p_w) / dx =
        ! = (p(1) - p'_w) / dx + (p(1) - p'_w) / 2 * (4 / pi) * (log(dx/r_w) + log(0.2))) / dx
        ! = (p(1) - p'_w) / dx * (1 + (2 / pi) * (log(dx/r_w) + log(0.2)))
        !Thus (p(1) - p'_w) / dx = (p(1) - p_w) / (dx * (1 + (2 / pi) * (log(dx/r_w) + log(0.2))))

#       if (_DARCY_LAYERS > 0)
            subroutine get_areas_and_lengths(element, dx, dy, dz, Ax, Ay, Az)
                type(t_element_base), intent(inout)				:: element
                real (kind = GRID_SR), intent(inout)		    :: dx, dy, dz, Ax, Ay, Az

                integer (kind = SI)                             :: boundary_condition(3)

                dx = element%cell%geometry%get_leg_size()
                dy = dx
                dz = cfg%dz
                Ax = 0.5_SR * dx * dz
                Ay = Ax
                Az = element%cell%geometry%get_volume()

#               if defined(_ASAGI)
                    call gv_boundary_condition%read_from_element(element, boundary_condition)

                    if (any(boundary_condition .ne. 0)) then
                        if (boundary_condition(1) .ne. 0 .or. boundary_condition(2) .ne. 0) then
                            dx = dx + dx * 2.0_SR / PI * (log(dx / cfg%r_well_radius) + log(0.2_SR))
                        end if

                        if (boundary_condition(3) .ne. 0 .or. boundary_condition(2) .ne. 0) then
                            dy = dy + dy * 2.0_SR / PI * (log(dy / cfg%r_well_radius) + log(0.2_SR))
                        end if
                    end if
#               endif
            end subroutine
#       else
            subroutine get_areas_and_lengths(element, dx, dy, Ax, Ay)
                type(t_element_base), intent(inout)				:: element
                real (kind = GRID_SR), intent(inout)		    :: dx, dy, Ax, Ay

                integer (kind = SI)                             :: boundary_condition(3)

                dx = element%cell%geometry%get_leg_size()
                dy = dx
                Ax = 0.5_SR * dx
                Ay = Ax

#               if defined(_ASAGI)
                    call gv_boundary_condition%read_from_element(element, boundary_condition)

                    if (any(boundary_condition .ne. 0)) then
                        if (boundary_condition(1) .ne. 0 .or. boundary_condition(2) .ne. 0) then
                            dx = dx + dx * 2.0_SR / PI * (log(dx / cfg%r_well_radius) + log(0.2_SR))
                        end if

                        if (boundary_condition(3) .ne. 0 .or. boundary_condition(2) .ne. 0) then
                            dy = dy + dy * 2.0_SR / PI * (log(dy / cfg%r_well_radius) + log(0.2_SR))
                        end if
                    end if
#               endif
            end subroutine
#       endif

	END MODULE
#endif
