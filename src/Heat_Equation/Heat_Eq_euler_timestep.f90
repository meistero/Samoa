! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_HEAT_EQ)
	MODULE Heat_Eq_Euler_Timestep
		use SFC_node_traversal
		use Samoa_heat_eq

		implicit none

		PUBLIC timestep_pre_traversal_op, timestep_post_traversal_op, timestep_element_op, timestep_cell_first_touch_op, timestep_edge_first_touch_op, timestep_node_first_touch_op, timestep_bnd_edge_first_touch_op, timestep_bnd_node_first_touch_op, timestep_bnd_edge_last_touch_op, timestep_bnd_node_last_touch_op, timestep_count_issued_refinements

		type(heat_eq_gv_T)						:: gv_T
		type(heat_eq_gv_r)						:: gv_r
		type(heat_eq_gv_mat_mass_diagonal)		:: gv_mat_mass_diagonal

		real (kind = GRID_SR), dimension(2)		:: laser_pos
		integer (kind = GRID_DI)				:: i_refinements_issued

#		define _GT_NAME							euler_timestep_traversal

#		if (_HEAT_EQ_EDGE_SIZE > 0)
#			define _GT_EDGES
#			define _GT_EDGES_TEMP
#		endif

#		define _GT_NODES
#		define _GT_NODES_TEMP
#		define _GT_REFINEMENTS

#		define _GT_PRE_TRAVERSAL_OP				timestep_pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP			timestep_post_traversal_op

#		define _GT_ELEMENT_OP					timestep_element_op

#		define _GT_CELL_FIRST_TOUCH_OP			timestep_cell_first_touch_op
#		define _GT_INNER_EDGE_FIRST_TOUCH_OP	timestep_edge_first_touch_op
#		define _GT_INNER_NODE_FIRST_TOUCH_OP	timestep_node_first_touch_op
#		define _GT_EDGE_FIRST_TOUCH_OP			timestep_bnd_edge_first_touch_op
#		define _GT_NODE_FIRST_TOUCH_OP			timestep_bnd_node_first_touch_op

#		define _GT_CELL_LAST_TOUCH_OP			timestep_cell_last_touch_op
#		define _GT_INNER_EDGE_LAST_TOUCH_OP		timestep_edge_last_touch_op
#		define _GT_INNER_NODE_LAST_TOUCH_OP		timestep_node_last_touch_op
#		define _GT_EDGE_LAST_TOUCH_OP			timestep_bnd_edge_last_touch_op
#		define _GT_NODE_LAST_TOUCH_OP			timestep_bnd_node_last_touch_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		function timestep_count_issued_refinements()
			integer (kind = GRID_SI)			:: timestep_count_issued_refinements

			timestep_count_issued_refinements = i_refinements_issued
		end function

		!*******************************
		!Geometry operators
		!*******************************

		subroutine timestep_pre_traversal_op(traversal, section)
 			type(t_grid_section), intent(inout)							:: grid
            integer (kind = GRID_SI)                                    :: i, i_color

			laser_pos = [0.50125_GRID_SR, 0.50715_GRID_SR] - 0.31747_GRID_SR * [cos(2.0_GRID_SR * PI * grid%r_laser_rps * grid%r_time), sin(2.0_GRID_SR * PI * grid%r_laser_rps * grid%r_time)]

			!this variable will be incremented for each cell with a refinement request
			i_refinements_issued = 0
		end subroutine

		subroutine timestep_post_traversal_op(traversal, section)
 			type(t_grid_section), intent(inout)							:: grid
		end subroutine

		subroutine timestep_element_op(traversal, section, element)
			type(t_grid_section), intent(inout)						:: grid
			type(t_element_base), intent(inout)	        :: element

			!local variables

			real (kind = GRID_SR), dimension(_HEAT_EQ_SIZE)	:: T
			real (kind = GRID_SR), dimension(_HEAT_EQ_SIZE)	:: r
			real (kind = GRID_SR), dimension(_HEAT_EQ_SIZE)	:: mat_mass_diagonal

			call gv_T%read(element, T)

			!call element operator
			call timestep_alpha_volume(grid, element, T, r, mat_mass_diagonal, element%cell%data_pers%heat_conductivity)

			call gv_r%add(element, r)
			call gv_mat_mass_diagonal%add(element, mat_mass_diagonal)
		end subroutine

		subroutine timestep_cell_first_touch_op(traversal, section, cell)
			type(t_grid_section), intent(inout)												:: grid
			type(t_cell_data_ptr), intent(inout)				:: cell

			call timestep_pre_dof(cell%data_temp%r, cell%data_temp%mat_mass_diagonal)
		end subroutine

		subroutine timestep_edge_first_touch_op(traversal, section, edge)
			type(t_grid_section), intent(inout)												:: grid
			type(t_edge_data), intent(inout)		:: edge

			call timestep_pre_dof(edge%data_temp%r, edge%data_temp%mat_mass_diagonal)
		end subroutine

		subroutine timestep_node_first_touch_op(traversal, section, node)
			type(t_grid_section), intent(inout)												:: grid
			type(t_node_data), intent(inout)				:: node

			call timestep_pre_dof(node%data_temp%r, node%data_temp%mat_mass_diagonal)
		end subroutine

		subroutine timestep_bnd_edge_first_touch_op(traversal, section, edge)
			type(t_grid_section), intent(inout)												:: grid
			type(t_edge_data), intent(inout)		:: edge
		end subroutine

		subroutine timestep_bnd_node_first_touch_op(traversal, section, node)
			type(t_grid_section), intent(inout)												:: grid
			type(t_node_data), intent(inout)				:: node
		end subroutine

		subroutine timestep_cell_last_touch_op(traversal, section, cell)
			type(t_grid_section), intent(inout)												:: grid
			type(t_cell_data_ptr), intent(inout)				:: cell

			call euler_timestep_post_dof(grid, cell%data_pers%T, cell%data_temp%r, cell%data_temp%mat_mass_diagonal)
		end subroutine

		subroutine timestep_edge_last_touch_op(traversal, section, edge)
			type(t_grid_section), intent(inout)												:: grid
			type(t_edge_data), intent(inout)		:: edge

			call euler_timestep_post_dof(grid, edge%data_pers%T, edge%data_temp%r, edge%data_temp%mat_mass_diagonal)
		end subroutine

		subroutine timestep_bnd_edge_last_touch_op(traversal, section, edge)
			type(t_grid_section), intent(inout)												:: grid
			type(t_edge_data), intent(inout)		:: edge

			edge%data_pers%T_temp = 0
		end subroutine

		subroutine timestep_bnd_node_last_touch_op(traversal, section, node)
			type(t_grid_section), intent(inout)												:: grid
			type(t_node_data), intent(inout)				:: node

			node%data_pers%T_temp = 0
		end subroutine

		subroutine timestep_node_last_touch_op(traversal, section, node)
			type(t_grid_section), intent(inout)												:: grid
			type(t_node_data), intent(inout)				:: node

			call euler_timestep_post_dof(grid, node%data_pers%T, node%data_temp%r, node%data_temp%mat_mass_diagonal)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		elemental subroutine timestep_pre_dof(r, mat_mass_diagonal)
			real(kind = GRID_SR), intent(out)				:: r
			real(kind = GRID_SR), intent(out)				:: mat_mass_diagonal

			!init temp variables to 0

			r = 0.0_GRID_SR
			mat_mass_diagonal = tiny(1.0_GRID_SR)
		end subroutine

		subroutine timestep_alpha_volume(grid, element, T, r, mat_mass_diagonal, heat_conductivity)
			type(t_grid_section), intent(inout)						:: grid
			type(t_element_base), intent(inout)						:: element
			real (kind = GRID_SR), dimension(:), intent(inout)		:: T
			real (kind = GRID_SR), dimension(:), intent(inout)		:: r
			real (kind = GRID_SR), dimension(:), intent(inout)		:: mat_mass_diagonal
			real (kind = GRID_SR), intent(in)						:: heat_conductivity

			real (kind = GRID_SR), dimension(2)						:: laser_dist
			real (kind = GRID_SR)									:: T_norm
			integer (kind = GRID_SI)								:: i, depth

			!add up global mass matrix diagonal
			forall (i = 1 : _HEAT_EQ_SIZE)
				mat_mass_diagonal(i) = gm_mass_matrix(i, i)
			end forall

			mat_mass_diagonal = abs(element%transform_data%plotter_data%det_jacobian) * (element%transform_data%custom_data%scaling ** 2) * mat_mass_diagonal

			!set residual
			r = -heat_conductivity * MATMUL(gm_stiffness_matrix, T)

			!**set refinement info**

			!set default refinement to none
			depth = element%cell%geometry%i_depth
			T_norm = abs(T(3) - T(2)) + abs(T(1) - T(2))

			if (depth < cfg%i_max_depth .and. T_norm > 3.0e-2_GRID_SR) then
				_log_write(5, "(A, T30, A, I0)") "  refinement issued:", "depth ", depth

				element%cell%geometry%refinement = 1
				i_refinements_issued = i_refinements_issued + 1
			else if (depth > cfg%i_min_depth .and. T_norm < 1.0e-2_GRID_SR) then
				_log_write(5, "(A, T30, A, I0)") "  coarsening issued:", "depth ", depth
				element%cell%geometry%refinement = -1
			else
				element%cell%geometry%refinement = 0
			end if

			!**add external heat source**

			!get barycentric coordinates of the laser position rotating with a given angular frequency
			laser_dist = samoa_world_to_barycentric_point(element%transform_data, laser_pos)

			!check first if the laser is close enough to the cell to intersect it
			if (dot_product(laser_dist, laser_dist) <= 1.0_GRID_SR) then
				!now check if it is actually inside the cell
				if (laser_dist(1) >= 0.0_GRID_SR .and. laser_dist(2) >= 0.0_GRID_SR .and. laser_dist(1) + laser_dist(2) <= 1.0_GRID_SR) then
					!assuming a dirac-delta-function => element integral reduced to point-evaluation of the basis functions
					r = r + heat_conductivity * samoa_basis_T_at(laser_dist)

					!also, always refine in the heat source
					if (depth < cfg%i_max_depth) then
						_log_write(5, "(A, T30, A, I0)") "  refinement issued:", "depth ", depth

						element%cell%geometry%refinement = 1
						i_refinements_issued = i_refinements_issued + 1
					end if
				end if
			end if
		end subroutine

		subroutine euler_timestep_post_dof(grid, T, r, mat_mass_diagonal)
			type(t_grid_section), intent(inout)							:: grid

			real(kind = GRID_SR), dimension(:), intent(inout)	:: T
			real(kind = GRID_SR), dimension(:), intent(in)		:: r
			real(kind = GRID_SR), dimension(:), intent(in)		:: mat_mass_diagonal

			T = T + grid%r_dt * r / mat_mass_diagonal
		end subroutine
	END MODULE
#endif
