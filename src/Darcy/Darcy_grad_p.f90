! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_grad_p
		use SFC_edge_traversal

		use Samoa_darcy

        type num_traversal_data
			real (kind = GRID_SR)				:: r_dt					!< global minimum time step
        end type

		type(darcy_gv_p)						:: gv_p
		type(darcy_gv_u)						:: gv_u

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

            type(t_node_data)               :: node
            integer                         :: blocklengths(2), types(2), disps(2), i_error, extent

#           if defined(_MPI)
                blocklengths(1) = 1
                blocklengths(2) = 1

                disps(1) = 0
                disps(2) = sizeof(node)

                types(1) = MPI_LB
                types(2) = MPI_UB

                call MPI_Type_struct(2, blocklengths, disps, types, mpi_node_type, i_error); assert_eq(i_error, 0)
                call MPI_Type_commit(mpi_node_type, i_error); assert_eq(i_error, 0)

                call MPI_Type_extent(mpi_node_type, extent, i_error); assert_eq(i_error, 0)
                assert_eq(sizeof(node), extent)

                call MPI_Type_size(mpi_node_type, extent, i_error); assert_eq(i_error, 0)
                assert_eq(0, extent)
#           endif
        end subroutine

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

		subroutine post_traversal_grid_op(traversal, grid)
			type(t_darcy_grad_p_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							:: grid

			call reduce(grid%r_dt, traversal%children%r_dt, MPI_MIN, .true.)
			grid%r_dt = cfg%courant_number * cfg%r_nu_w * cfg%scaling * grid%r_dt
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
 			type(t_darcy_grad_p_traversal), intent(inout)				:: traversal
 			type(t_grid_section), intent(inout)				:: section
			type(t_element_base), intent(inout), target		:: element

			real (kind = GRID_SR) :: p(_DARCY_LAYERS, 3)
			real (kind = GRID_SR) :: u(_DARCY_LAYERS, 2)
			integer (kind = GRID_SI)											:: i

			call gv_p%read_from_element(element, p)

			!call element operator
			do i = 1, _DARCY_LAYERS
                call alpha_volume_op(traversal, element, p(i, :), u(i, :), element%cell%data_pers%base_permeability, element%cell%data_pers%porosity)
                !call alpha_volume_op(traversal, element, p(i, :), u(i, :), element%cell%data_pers%base_permeability(i), element%cell%data_pers%porosity(i))
			end do

			call gv_u%write_to_element(element, u)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine alpha_volume_op(traversal, element, p, u, base_permeability, porosity)
 			type(t_darcy_grad_p_traversal), intent(inout)						:: traversal
			type(t_element_base), intent(inout)									:: element
			real (kind = GRID_SR), intent(in)					                :: p(:)
			real (kind = GRID_SR), intent(out)					                :: u(:)
			real (kind = GRID_SR), intent(in)									:: base_permeability
			real (kind = GRID_SR), intent(in)									:: porosity

			real (kind = SR)                                                    :: u_norm

			!define velocity by u = k (-grad p + rho g)
            u = samoa_barycentric_to_world_normal(element%transform_data, 0.5_SR * [p(1) - p(2), p(3) - p(2)])
            u = u * (element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))
            assert_le(norm2(u), (1.0_SR + 1.0e2_SR * epsilon(1.0_SR)) * norm2(0.5_SR * [p(1) - p(2), p(3) - p(2)]))
            assert_le(norm2(0.5_SR * [p(1) - p(2), p(3) - p(2)]), (1.0_SR + 1.0e2_SR * epsilon(1.0_SR)) * norm2(u))
            u = base_permeability * (-u + cfg%r_rho_w * g)

            u_norm = norm2(u)
            if (u_norm > 0.0_SR .and. porosity > 0.0_SR) then
                traversal%r_dt = min(traversal%r_dt, porosity * element%cell%geometry%get_leg_size() / (2.0_SR * u_norm))
            end if
		end subroutine
	END MODULE
#endif
