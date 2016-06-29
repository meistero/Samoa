! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_grad_p
		use SFC_edge_traversal
		use Darcy_permeability

		use Samoa_darcy

        type num_traversal_data
			real (kind = GRID_SR)				:: r_dt					!< global minimum time step
        end type

		type(darcy_gv_p)				        :: gv_p
		type(darcy_gv_saturation)				:: gv_saturation
		type(darcy_gv_volume)					:: gv_volume
		type(darcy_gv_flux)						:: gv_xi_w
		type(darcy_gv_d)						:: gv_xi_n

#		define _GT_NAME							t_darcy_grad_p_traversal

#		define _GT_NODES

#		define _GT_PRE_TRAVERSAL_OP		        pre_traversal_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op

#		define _GT_ELEMENT_OP					element_op

#		define _GT_NODE_FIRST_TOUCH_OP		    node_first_touch_op
#		define _GT_NODE_REDUCE_OP		        node_reduce_op

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

		subroutine pre_traversal_op(traversal, section)
  			type(t_darcy_grad_p_traversal), intent(inout)		:: traversal
 			type(t_grid_section), intent(inout)				    :: section

			traversal%r_dt = huge(1.0_SR)
		end subroutine

		subroutine post_traversal_grid_op(traversal, grid)
			type(t_darcy_grad_p_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							:: grid

			call reduce(traversal%r_dt, traversal%sections%r_dt, MPI_MIN, .true.)
			grid%r_dt = cfg%courant_number * traversal%r_dt

			if (cfg%r_max_time > 0.0_SR) then
                grid%r_dt = min(cfg%r_max_time, grid%r_dt)
            end if

			if (cfg%r_output_time_step > 0.0_SR) then
                grid%r_dt = min(cfg%r_output_time_step, grid%r_dt)
            end if
		end subroutine

		!*******************************
		!Geometry operators
		!*******************************

		subroutine element_op(traversal, section, element)
 			type(t_darcy_grad_p_traversal), intent(inout)	    :: traversal
 			type(t_grid_section), intent(inout)						:: section
			type(t_element_base), intent(inout)			            :: element

#           if (_DARCY_LAYERS > 0)
                real(kind = GRID_SR)    :: p(_DARCY_LAYERS + 1, 3)
                real(kind = GRID_SR)    :: saturation(_DARCY_LAYERS + 1, 3)
                real(kind = GRID_SR)    :: xi_w(_DARCY_LAYERS + 1, 3)
                real(kind = GRID_SR)    :: xi_n(_DARCY_LAYERS + 1, 3)
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
                real(kind = GRID_SR)    :: xi_w(3)
                real(kind = GRID_SR)    :: xi_n(3)
                real(kind = GRID_SR)    :: volume(3)
                real(kind = GRID_SR)    :: porosity

                porosity = element%cell%data_pers%porosity
#           endif

			call gv_p%read_from_element(element, p)
			call gv_saturation%read_from_element(element, saturation)

			!call volume operator
			call compute_wavespeeds(element, p, saturation, volume, element%cell%data_pers%base_permeability, porosity, xi_w, xi_n)

			call gv_volume%add_to_element(element, volume)
			call gv_xi_w%add_to_element(element, xi_w)
			call gv_xi_n%add_to_element(element, xi_n)
		end subroutine


		!first touches

		elemental subroutine node_first_touch_op(traversal, section, node)
 			type(t_darcy_grad_p_traversal), intent(in)	    :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)				:: node

			call pre_dof_op(node%data_temp%flux, node%data_pers%d, node%data_temp%volume)
		end subroutine

		!last touches

		subroutine node_reduce_op(traversal, section, node)
 			type(t_darcy_grad_p_traversal), intent(inout)       :: traversal
 			type(t_grid_section), intent(in)				    :: section
			type(t_node_data), intent(in)				        :: node

            integer :: i

            do i = 1, _DARCY_LAYERS + 1
                call reduce_dof_op(node%data_temp%flux(i), node%data_pers%d(i), node%data_temp%volume(i), traversal%r_dt)
            end do
		end subroutine

		elemental subroutine node_merge_op(local_node, neighbor_node)
			type(t_node_data), intent(inout)		:: local_node
 			type(t_node_data), intent(in)		    :: neighbor_node

			local_node%data_temp%flux = local_node%data_temp%flux + neighbor_node%data_temp%flux
			local_node%data_pers%d = local_node%data_pers%d + neighbor_node%data_pers%d
			local_node%data_temp%volume = local_node%data_temp%volume + neighbor_node%data_temp%volume
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		elemental subroutine pre_dof_op(xi_w, xi_n, volume)
 			real (kind = GRID_SR), intent(out)		:: xi_w, xi_n
			real (kind = GRID_SR), intent(out)		:: volume

			xi_w = 0.0_SR
			xi_n = 0.0_SR
			volume = 0.0_SR
		end subroutine

		pure subroutine reduce_dof_op(xi_w, xi_n, volume, dt)
			real (kind = GRID_SR), intent(in)		:: xi_w, xi_n
			real (kind = GRID_SR), intent(in)		:: volume
			real (kind = GRID_SR), intent(inout)	:: dt

            if ((xi_w + xi_n) * volume > 1.0e5_SR * tiny(1.0_SR)) then
                dt = min(dt, volume / (xi_w + xi_n))
            end if
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

        subroutine compute_wavespeeds(element, p, saturation, volume, base_permeability, porosity, xi_w, xi_n)
			type(t_element_base), intent(inout)	    :: element

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR), intent(in)	    :: p(:, :)
                real (kind = GRID_SR), intent(in)	    :: base_permeability(:,:)
                real (kind = GRID_SR), intent(in)	    :: porosity(:)
                real (kind = GRID_SR), intent(in)	    :: saturation(:, :)
                real (kind = GRID_SR), intent(out)	    :: volume(:, :), xi_w(:, :), xi_n(:, :)

                real (kind = GRID_SR)                   :: u_w(_DARCY_LAYERS, 7), u_n(_DARCY_LAYERS, 7)

                real (kind = GRID_SR)				    :: g_local(3), dx, dy, dz, Ax, Ay, Az

                volume(:, 1) = cfg%dz * 0.25_SR * porosity(:) * element%cell%geometry%get_volume()
                volume(:, 2) = cfg%dz * 0.50_SR * porosity(:) * element%cell%geometry%get_volume()
                volume(:, 3) = cfg%dz * 0.25_SR * porosity(:) * element%cell%geometry%get_volume()

                xi_w = 0.0_SR
                xi_n = 0.0_SR

                !rotate g so it points in the right direction (no scaling!)
                g_local = cfg%g
                g_local(1:2) = samoa_world_to_barycentric_normal(element%transform_data, g_local(1:2))

                call get_areas_and_lengths(element, dx, dy, dz, Ax, Ay, Az)
                call compute_base_fluxes_3D(p, base_permeability, dx, dy, dz, Ax, Ay, Az, g_local, u_w, u_n)
                call compute_wave_speeds_3D(saturation, u_w, u_n, xi_w, xi_n)
#           else
                real (kind = GRID_SR), intent(in)	    :: p(:)
                real (kind = GRID_SR), intent(in)	    :: base_permeability
                real (kind = GRID_SR), intent(in)	    :: porosity
                real (kind = GRID_SR), intent(in)	    :: saturation(:)
                real (kind = GRID_SR), intent(out)	    :: volume(:), xi_w(:), xi_n(:)

                real (kind = GRID_SR)					:: g_local(2), dx, dy, Ax, Ay
                real (kind = SR)                        :: u_w(2), u_n(2), u_norm

                volume = [0.25_SR, 0.50_SR, 0.25_SR] * porosity * element%cell%geometry%get_volume()

                xi_w = 0.0_SR
                xi_n = 0.0_SR

                !rotate g so it points in the right direction (no scaling!)
                g_local = cfg%g(1:2)
                g_local(1:2) = samoa_world_to_barycentric_normal(element%transform_data, g_local(1:2))

                !compute base fluxes

                call get_areas_and_lengths(element, dx, dy, Ax, Ay)
                call compute_base_fluxes_2D(p, base_permeability, dx, dy, Ax, Ay, g_local(1:2), u_w, u_n)
                call compute_wave_speeds_2D(saturation, u_w, u_n, xi_w, xi_n)
#           endif
		end subroutine
	END MODULE
#endif
