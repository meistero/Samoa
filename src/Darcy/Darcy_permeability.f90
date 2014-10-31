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

		type(darcy_gv_saturation)				:: gv_saturation
		type(darcy_gv_p)						:: gv_p
		type(darcy_gv_rhs)						:: gv_rhs
		type(darcy_gv_is_dirichlet_boundary)    :: gv_is_dirichlet

#		define _GT_NAME							t_darcy_permeability_traversal

#		if (_DARCY_FLOW_EDGE_SIZE > 0)
#			define _GT_EDGES
#		endif

#		define _GT_NODES

#		define _GT_ELEMENT_OP					element_op

#		define _GT_NODE_FIRST_TOUCH_OP		    node_first_touch_op

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

		elemental subroutine node_first_touch_op(traversal, section, node)
 			type(t_darcy_permeability_traversal), intent(in)                     :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			:: node

			node%data_pers%rhs = 0.0_SR
			node%data_temp%is_dirichlet_boundary = .false.
		end subroutine

		subroutine element_op(traversal, section, element)
 			type(t_darcy_permeability_traversal), intent(inout)		:: traversal
 			type(t_grid_section), intent(inout)						:: section
			type(t_element_base), intent(inout), target				:: element

			real (kind = GRID_SR)   :: saturation(3)
			real (kind = GRID_SR)   :: p(3)
			real (kind = GRID_SR)   :: rhs(3)
			logical		            :: is_dirichlet(3)

			call gv_saturation%read_from_element(element, saturation)
			call gv_p%read_from_element(element, p)

			!call element operator
			call alpha_volume_op(element, saturation, p, rhs, is_dirichlet, element%cell%data_pers%base_permeability)

			call gv_rhs%add_to_element(element, rhs)
		end subroutine

		elemental subroutine node_merge_op(local_node, neighbor_node)
			type(t_node_data), intent(inout)		:: local_node
 			type(t_node_data), intent(in)		    :: neighbor_node

            if (any(neighbor_node%data_temp%is_dirichlet_boundary)) then
                local_node%data_pers%p = neighbor_node%data_pers%p
            end if

			local_node%data_pers%saturation = max(local_node%data_pers%saturation, neighbor_node%data_pers%saturation)
			local_node%data_temp%is_dirichlet_boundary = local_node%data_temp%is_dirichlet_boundary .or. neighbor_node%data_temp%is_dirichlet_boundary
			local_node%data_pers%rhs = local_node%data_pers%rhs + neighbor_node%data_pers%rhs
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine alpha_volume_op(element, saturation, p, rhs, is_dirichlet, base_permeability)
			type(t_element_base), intent(inout)				                    :: element
			real (kind = GRID_SR), intent(inout)		                        :: saturation(:)
			real (kind = GRID_SR), intent(inout)			                    :: p(:)
			real (kind = GRID_SR), intent(out)									:: rhs(:)
			logical, intent(out)			                                    :: is_dirichlet(:)
			real (kind = GRID_SR), intent(in)									:: base_permeability

			real (kind = GRID_SR), parameter            :: volumes(3) = [1.0_SR/8.0_SR, 1.0_SR/4.0_SR, 1.0_SR/8.0_SR]

			real (kind = GRID_SR)					    :: u_w(2), g_local(2), pos_prod(2), pos_in(2), lambda_w(3), lambda_n(3), radius, dual_edge_length
            real (kind = GRID_SR)                       :: lambda_t(2)

#           if (_DARCY_LAYERS > 1)
                real (kind = GRID_SR)                   :: K_base(2)

                K_base = element%cell%data_pers%base_permeability
#           else
                real (kind = GRID_SR)                   :: K_base

                K_base = element%cell%data_pers%base_permeability
#           endif

            rhs(:) = 0.0_SR
            radius = 0.0508_SR / (cfg%scaling * element%transform_data%custom_data%scaling)

            pos_prod = sign(cfg%r_pos_prod - 0.5_SR, element%transform_data%custom_data%offset - 0.5_SR) + 0.5_SR
            pos_prod = samoa_world_to_barycentric_point(element%transform_data, pos_prod)
            if (pos_prod(1) .ge. -radius .and. pos_prod(2) .ge. -radius .and. pos_prod(1) + pos_prod(2) .le. 1.0_SR + radius) then
                !production well:
                !set a constant pressure condition and an outflow saturation condition
                is_dirichlet = .true.
                p = cfg%r_p_prod

                call gv_p%write(element, p)
                call gv_is_dirichlet%add(element, is_dirichlet)
            end if

            pos_in = samoa_world_to_barycentric_point(element%transform_data, cfg%r_pos_in)
            if (base_permeability > 0.0_SR .and. pos_in(1) .ge. -radius .and. pos_in(2) .ge. -radius .and. pos_in(1) + pos_in(2) .le. 1.0_SR + radius) then
                !injection well:
                !set an inflow pressure condition and a constant saturation condition
                saturation = 1.0_SR
                call gv_saturation%write(element, saturation)

                rhs(:) = 0.009201_SR / 51.816_SR * min(1.0_SR, 1.0_SR / (0.0508_SR * 0.0508_SR * pi) * ((cfg%scaling * element%transform_data%custom_data%scaling) ** 2)) * [0.125_SR, 0.25_SR, 0.125_SR]
            end if

			lambda_w = (saturation * saturation) / cfg%r_nu_w
			lambda_n = (1.0_SR - saturation) * (1.0_SR - saturation) / cfg%r_nu_n

            g_local = samoa_world_to_barycentric_normal(element%transform_data, g)
            g_local = g_local / (cfg%scaling * element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))
            dual_edge_length = 0.5_SR * element%cell%geometry%get_leg_size()

#           if (_DARCY_LAYERS > 1)
                u_w = K_base(1) * (-[p(1) - p(2), p(3) - p(2)] + dual_edge_length * cfg%r_rho_w * g_local)
#           else
                u_w = K_base * (-[p(1) - p(2), p(3) - p(2)] + dual_edge_length * cfg%r_rho_w * g_local)
#           endif

            if (u_w(1) > 0) then
#               if (_DARCY_LAYERS > 1)
                    rhs(1:2) = rhs(1:2) + K_base(1) * (cfg%r_rho_w * lambda_w(2) + cfg%r_rho_n * lambda_n(2)) * dual_edge_length * g_local)
#               else
                    rhs(1:2) = rhs(1:2) + K_base * (cfg%r_rho_w * lambda_w(2) + cfg%r_rho_n * lambda_n(2)) * dual_edge_length * g_local(1)
#               endif

                lambda_t(1) = (lambda_w(2) + lambda_n(2))
            else
#               if (_DARCY_LAYERS > 1)
                    rhs(1:2) = rhs(1:2) + K_base(1) * (cfg%r_rho_w * lambda_w(1) + cfg%r_rho_n * lambda_n(1)) * dual_edge_length * g_local)
#               else
                    rhs(1:2) = rhs(1:2) + K_base * (cfg%r_rho_w * lambda_w(1) + cfg%r_rho_n * lambda_n(1)) * dual_edge_length * g_local(1)
#               endif

                lambda_t(1) = (lambda_w(1) + lambda_n(1))
            end if

            if (u_w(2) > 0) then
#               if (_DARCY_LAYERS > 1)
                    rhs(2:3) = rhs(2:3) + K_base(1) * (cfg%r_rho_w * lambda_w(2) + cfg%r_rho_n * lambda_n(2)) * dual_edge_length * g_local)
#               else
                    rhs(2:3) = rhs(2:3) + K_base * (cfg%r_rho_w * lambda_w(2) + cfg%r_rho_n * lambda_n(2)) * dual_edge_length * g_local(2)
#               endif

                lambda_t(2) = (lambda_w(2) + lambda_n(2))
            else
#               if (_DARCY_LAYERS > 1)
                    rhs(2:3) = rhs(2:3) + K_base(1) * (cfg%r_rho_w * lambda_w(3) + cfg%r_rho_n * lambda_n(3)) * dual_edge_length * g_local)
#               else
                    rhs(2:3) = rhs(2:3) + K_base * (cfg%r_rho_w * lambda_w(3) + cfg%r_rho_n * lambda_n(3)) * dual_edge_length * g_local(2)
#               endif

                lambda_t(2) = (lambda_w(3) + lambda_n(3))
            end if

            if (element%transform_data%plotter_data%orientation > 0) then
                element%cell%data_pers%lambda_t = lambda_t
            else
                element%cell%data_pers%lambda_t(2) = lambda_t(1)
                element%cell%data_pers%lambda_t(1) = lambda_t(2)
            end if
		end subroutine
	END MODULE
#endif
