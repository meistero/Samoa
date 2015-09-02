! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_grad_p
		use SFC_edge_traversal
		use Darcy_initialize_saturation

		use Samoa_darcy

        type num_traversal_data
			real (kind = GRID_SR)				:: r_dt					!< global minimum time step
        end type

		type(darcy_gv_p)						:: gv_p
		type(darcy_gv_saturation)				:: gv_saturation

        public compute_velocity_1D

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

#           if defined(_MPI)
                type(t_node_data)                       :: node
                integer                                 :: blocklengths(2), types(2), disps(2), type_size, i_error
                integer (kind = MPI_ADDRESS_KIND)       :: lb, ub

                blocklengths(1) = 1
                blocklengths(2) = 1

                disps(1) = 0
                disps(2) = sizeof(node)

                types(1) = MPI_LB
                types(2) = MPI_UB

                call MPI_Type_struct(2, blocklengths, disps, types, mpi_node_type, i_error); assert_eq(i_error, 0)
                call MPI_Type_commit(mpi_node_type, i_error); assert_eq(i_error, 0)

                call MPI_Type_size(mpi_node_type, type_size, i_error); assert_eq(i_error, 0)
                call MPI_Type_get_extent(mpi_node_type, lb, ub, i_error); assert_eq(i_error, 0)

                assert_eq(0, lb)
                assert_eq(0, type_size)
                assert_eq(sizeof(node), ub)
#           endif
        end subroutine

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

		subroutine post_traversal_grid_op(traversal, grid)
			type(t_darcy_grad_p_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							:: grid

			call reduce(traversal%r_dt, traversal%children%r_dt, MPI_MIN, .true.)
			grid%r_dt = cfg%courant_number * traversal%r_dt

			if (cfg%r_max_time > 0.0_SR) then
                grid%r_dt = min(cfg%r_max_time, grid%r_dt)
            end if

			if (cfg%r_output_time_step > 0.0_SR) then
                grid%r_dt = min(cfg%r_output_time_step, grid%r_dt)
            end if
		end subroutine

		subroutine pre_traversal_op(traversal, section)
  			type(t_darcy_grad_p_traversal), intent(inout)		:: traversal
 			type(t_grid_section), intent(inout)				    :: section

			traversal%r_dt = huge(1.0_SR)
		end subroutine

		!*******************************
		!Geometry operators
		!*******************************

		subroutine element_op(traversal, section, element)
 			type(t_darcy_grad_p_traversal), intent(inout)	:: traversal
 			type(t_grid_section), intent(inout)				:: section
			type(t_element_base), intent(inout)		        :: element

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR) :: S(_DARCY_LAYERS + 1, 3)
                real (kind = GRID_SR) :: p(_DARCY_LAYERS + 1, 3)

                integer :: i

                call gv_saturation%read_from_element(element, S)
                call gv_p%read_from_element(element, p)

                call alpha_volume_op(traversal, element, S, p, element%cell%data_pers%base_permeability, element%cell%data_pers%porosity)
#           else
                real (kind = GRID_SR)       :: S(3)
                real (kind = GRID_SR)       :: p(3)

                call gv_saturation%read_from_element(element, S)
                call gv_p%read_from_element(element, p)

                call alpha_volume_op(traversal, element, S, p, element%cell%data_pers%base_permeability, element%cell%data_pers%porosity)
#           endif
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

#       if (_DARCY_LAYERS > 0)
            subroutine alpha_volume_op(traversal, element, S, p, base_permeability, porosity)
                type(t_darcy_grad_p_traversal), intent(inout)						:: traversal
                type(t_element_base), intent(inout)									:: element
                real (kind = GRID_SR), intent(in)					                :: S(:,:)
                real (kind = GRID_SR), intent(in)					                :: p(:,:)
                real (kind = GRID_SR), intent(in)									:: base_permeability(:,:)
                real (kind = GRID_SR), intent(in)									:: porosity(:)

                real (kind = GRID_SR)					  :: g_local(3), edge_length, surface, dz
                real (kind = SR)                          :: u_w(7), u_n(7), dt_inv_max
                integer                                   :: i

                !rotate g so it points in the right direction (no scaling!)
                g_local = cfg%g
                g_local(1:2) = samoa_world_to_barycentric_normal(element%transform_data, g_local(1:2))
                g_local(1:2) = g_local(1:2) / (element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))

                edge_length = element%cell%geometry%get_leg_size()
                surface = element%cell%geometry%get_volume()
                dz = cfg%dz

                do i = 1, _DARCY_LAYERS
                    call compute_velocity_1D(edge_length, 1.0_SR, base_permeability(i, 1), p(i, 2), p(i, 1), u_w(1), u_n(1), g_local(1))
                    call compute_velocity_1D(edge_length, 1.0_SR, base_permeability(i, 1), p(i, 2), p(i, 3), u_w(2), u_n(2), g_local(2))

                    call compute_velocity_1D(dz, 1.0_SR, base_permeability(i, 2), p(i, 1), p(i + 1, 1), u_w(3), u_n(3), g_local(3))
                    call compute_velocity_1D(dz, 1.0_SR, base_permeability(i, 2), p(i, 2), p(i + 1, 2), u_w(4), u_n(4), g_local(3))
                    call compute_velocity_1D(dz, 1.0_SR, base_permeability(i, 2), p(i, 3), p(i + 1, 3), u_w(5), u_n(5), g_local(3))

                    call compute_velocity_1D(edge_length, 1.0_SR, base_permeability(i, 1), p(i + 1, 2), p(i + 1, 1), u_w(6), u_n(6), g_local(1))
                    call compute_velocity_1D(edge_length, 1.0_SR, base_permeability(i, 1), p(i + 1, 2), p(i + 1, 3), u_w(7), u_n(7), g_local(2))

                    dt_inv_max = sum([  0.25_SR * edge_length * dz * compute_max_wave_speed_1D(S(i, 2), S(i, 1), u_w(1), u_n(1), base_permeability(i, 1), g_local(1)), &
                                        0.25_SR * edge_length * dz * compute_max_wave_speed_1D(S(i, 2), S(i, 3), u_w(2), u_n(2), base_permeability(i, 1), g_local(2)), &
                                        0.25_SR * surface * compute_max_wave_speed_1D(S(i, 1), S(i + 1, 1), u_w(3), u_n(3), base_permeability(i, 2), g_local(3)), &
                                        0.50_SR * surface * compute_max_wave_speed_1D(S(i, 2), S(i + 1, 2), u_w(4), u_n(4), base_permeability(i, 2), g_local(3)), &
                                        0.25_SR * surface * compute_max_wave_speed_1D(S(i, 3), S(i + 1, 3), u_w(5), u_n(5), base_permeability(i, 2), g_local(3)), &
                                        0.25_SR * edge_length * dz * compute_max_wave_speed_1D(S(i + 1, 2), S(i + 1, 1), u_w(6), u_n(6), base_permeability(i, 1), g_local(1)), &
                                        0.25_SR * edge_length * dz * compute_max_wave_speed_1D(S(i + 1, 2), S(i + 1, 3), u_w(7), u_n(7), base_permeability(i, 1), g_local(2))])

                    if (dt_inv_max * porosity(i) > 1.0e5_SR * tiny(1.0_SR)) then
                        traversal%r_dt = min(traversal%r_dt, porosity(i) * surface * dz / dt_inv_max)
                    end if
                end do
            end subroutine
#       else
            subroutine alpha_volume_op(traversal, element, S, p, base_permeability, porosity)
                type(t_darcy_grad_p_traversal), intent(inout)						:: traversal
                type(t_element_base), intent(inout)									:: element
                real (kind = GRID_SR), intent(in)					                :: S(:)
                real (kind = GRID_SR), intent(in)					                :: p(:)
                real (kind = GRID_SR), intent(in)									:: base_permeability
                real (kind = GRID_SR), intent(in)									:: porosity

                real (kind = GRID_SR)					  :: g_local(2), edge_length, surface
                real (kind = SR)                          :: u_w(2), u_n(2), dt_inv_max

                !rotate g so it points in the right direction (no scaling!)
                g_local = cfg%g(1:2)
                g_local(1:2) = samoa_world_to_barycentric_normal(element%transform_data, g_local(1:2))
                g_local(1:2) = g_local(1:2) / (element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))

                edge_length = element%cell%geometry%get_leg_size()
                surface = element%cell%geometry%get_volume()

                call compute_velocity_1D(edge_length, 1.0_SR, base_permeability, p(2), p(1), u_w(1), u_n(1), g_local(1))
                call compute_velocity_1D(edge_length, 1.0_SR, base_permeability, p(2), p(3), u_w(2), u_n(2), g_local(2))

                dt_inv_max =    0.5_SR * edge_length * compute_max_wave_speed_1D(S(2), S(1), u_w(1), u_n(1), base_permeability, g_local(1)) + &
                                0.5_SR * edge_length * compute_max_wave_speed_1D(S(2), S(3), u_w(2), u_n(2), base_permeability, g_local(2))

                if (dt_inv_max * porosity > 1.0e5_SR * tiny(1.0_SR)) then
                    traversal%r_dt = min(traversal%r_dt, porosity * surface / dt_inv_max)
                end if
            end subroutine
#       endif

        function compute_max_wave_speed_1D(S_l, S_r, u_w_in, u_n_in, permeability, g_local) result(max_wave_speed)
            real (kind = GRID_SR), intent(in)       :: S_l, S_r, u_w_in, u_n_in, permeability, g_local
            real (kind = GRID_SR)                   :: max_wave_speed

            real (kind = GRID_SR)                   :: lambda_wl, lambda_wr, lambda_nl, lambda_nr, u_w, u_n, du_1, du_2

            !Upwind and F-Wave solver approximate the Riemann solution by two shock waves with velocities xi_1 and xi_2:
            !xi_1 = (u_w(S_r) - F_w) / (S_r - S_l) = (u_n(S_r) - F_n) / ((1 - S_r) - (1 - S_l))
            !xi_2 = (F_w - u_w(S_l)) / (S_r - S_l) = (F_n - u_n(S_l)) / ((1 - S_r) - (1 - S_l))

            lambda_wl = l_w(S_l)
            lambda_wr = l_w(S_r)
            lambda_nl = l_n(S_l)
            lambda_nr = l_n(S_r)

            if (u_w_in > 0.0) then
                u_w = lambda_wl * u_w_in
            else
                u_w = lambda_wr * u_w_in
            endif

            if (u_n_in > 0.0) then
                u_n = lambda_nl * u_n_in
            else
                u_n = lambda_nr * u_n_in
            endif

            du_1 = u_w - lambda_wl / (lambda_wl + lambda_nl) * ((u_w + u_n) + lambda_nl * (cfg%r_rho_w - cfg%r_rho_n) * permeability * g_local)
            du_2 = lambda_wr / (lambda_wr + lambda_nr) * ((u_w + u_n) + lambda_nr * (cfg%r_rho_w - cfg%r_rho_n) * permeability * g_local) - u_w

            if (S_r .ne. S_l) then
                max_wave_speed = max(abs(du_1 / (S_r - S_l)), abs(du_2 / (S_r - S_l)))
            else
                !when the saturations are equal, all signals have amplitude 0
                max_wave_speed = 0.0_SR
            end if
        end function

        elemental subroutine compute_velocity_1D(dx, area, base_permeability, pL, pR, u_w, u_n, g_local)
            real (kind = GRID_SR), intent(in)       :: dx, area, base_permeability, pL, pR, g_local
            real (kind = GRID_SR), intent(inout)    :: u_w, u_n

            u_w = area * base_permeability * (-(pR - pL) / dx + cfg%r_rho_w * g_local)
            u_n = area * base_permeability * (-(pR - pL) / dx + cfg%r_rho_n * g_local)
        end subroutine
	END MODULE
#endif
