! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_HEAT_EQ)
	MODULE Heat_Eq_Heun_Timestep_1
		use SFC_node_traversal

		use Samoa_heat_eq
		use Heat_Eq_Euler_Timestep

		implicit none

		type(heat_eq_gv_T)						:: gv_T
		type(heat_eq_gv_T_temp)					:: gv_T_temp

#		define _GT_NAME							heun_timestep_1_traversal

#		if (_HEAT_EQ_EDGE_SIZE > 0)
#			define _GT_EDGES
#			define _GT_EDGES_TEMP
#		endif

#		define _GT_NODES
#		define _GT_NODES_TEMP
#		define _GT_REFINEMENTS

#		define _GT_PRE_TRAVERSAL_OP				timestep_pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP				timestep_post_traversal_op

#		define _GT_ELEMENT_OP					timestep_element_op

#		define _GT_CELL_FIRST_TOUCH_OP			timestep_cell_first_touch_op
#		define _GT_INNER_EDGE_FIRST_TOUCH_OP			timestep_edge_first_touch_op
#		define _GT_INNER_NODE_FIRST_TOUCH_OP			timestep_node_first_touch_op
#		define _GT_EDGE_FIRST_TOUCH_OP		timestep_bnd_edge_first_touch_op
#		define _GT_NODE_FIRST_TOUCH_OP		timestep_bnd_node_first_touch_op

#		define _GT_CELL_LAST_TOUCH_OP			heun_timestep_1_cell_last_touch_op
#		define _GT_INNER_EDGE_LAST_TOUCH_OP			heun_timestep_1_edge_last_touch_op
#		define _GT_INNER_NODE_LAST_TOUCH_OP			heun_timestep_1_node_last_touch_op
#		define _GT_EDGE_LAST_TOUCH_OP		timestep_bnd_edge_last_touch_op
#		define _GT_NODE_LAST_TOUCH_OP		timestep_bnd_node_last_touch_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		!*******************************
		!Geometry operators
		!*******************************

		subroutine heun_timestep_1_cell_last_touch_op(traversal, section, cell)
			type(t_grid_section), intent(inout)												:: grid

			type(t_cell_data_ptr), intent(inout)				:: cell

			call heun_timestep_1_post_dof(grid, cell%data_pers%T, cell%data_pers%T_temp, cell%data_temp%r, cell%data_temp%mat_mass_diagonal)
		end subroutine

		subroutine heun_timestep_1_edge_last_touch_op(traversal, section, edge)
			type(t_grid_section), intent(inout)												:: grid

			type(t_edge_data), intent(inout)		:: edge

			call heun_timestep_1_post_dof(grid, edge%data_pers%T, edge%data_pers%T_temp, edge%data_temp%r, edge%data_temp%mat_mass_diagonal)
		end subroutine

		subroutine heun_timestep_1_node_last_touch_op(traversal, section, node)
			type(t_grid_section), intent(inout)												:: grid

			type(t_node_data), intent(inout)				:: node

			call heun_timestep_1_post_dof(grid, node%data_pers%T, node%data_pers%T_temp, node%data_temp%r, node%data_temp%mat_mass_diagonal)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine heun_timestep_1_post_dof(grid, T, T_temp, r, mat_mass_diagonal)
			type(t_grid_section), intent(inout)									:: grid
			real(kind = GRID_SR), dimension(:), intent(inout)			:: T
			real(kind = GRID_SR), dimension(:), intent(out)				:: T_temp
			real(kind = GRID_SR), dimension(:), intent(in)				:: r
			real(kind = GRID_SR), dimension(:), intent(in)				:: mat_mass_diagonal

			T_temp = T + 0.5_GRID_SR * grid%r_dt * r / mat_mass_diagonal
			T = T + grid%r_dt * r / mat_mass_diagonal
		end subroutine
	END MODULE

	MODULE Heat_Eq_Heun_Timestep_2
		use SFC_node_traversal

		use Samoa_heat_eq
		use Heat_Eq_Euler_Timestep

		implicit none

		type(heat_eq_gv_T)						:: gv_T
		type(heat_eq_gv_T_temp)					:: gv_T_temp

#		define _GT_NAME							heun_timestep_2_traversal

#		if (_HEAT_EQ_EDGE_SIZE > 0)
#			define _GT_EDGES
#			define _GT_EDGES_TEMP
#		endif

#		define _GT_NODES
#		define _GT_NODES_TEMP
#		define _GT_REFINEMENTS

#		define _GT_PRE_TRAVERSAL_OP				timestep_pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP				timestep_post_traversal_op

#		define _GT_ELEMENT_OP					timestep_element_op

#		define _GT_CELL_FIRST_TOUCH_OP			timestep_cell_first_touch_op
#		define _GT_INNER_EDGE_FIRST_TOUCH_OP			timestep_edge_first_touch_op
#		define _GT_INNER_NODE_FIRST_TOUCH_OP			timestep_node_first_touch_op
#		define _GT_EDGE_FIRST_TOUCH_OP		timestep_bnd_edge_first_touch_op
#		define _GT_NODE_FIRST_TOUCH_OP		timestep_bnd_node_first_touch_op

#		define _GT_CELL_LAST_TOUCH_OP			heun_timestep_2_cell_last_touch_op
#		define _GT_INNER_EDGE_LAST_TOUCH_OP			heun_timestep_2_edge_last_touch_op
#		define _GT_INNER_NODE_LAST_TOUCH_OP			heun_timestep_2_node_last_touch_op
#		define _GT_EDGE_LAST_TOUCH_OP		timestep_bnd_edge_last_touch_op
#		define _GT_NODE_LAST_TOUCH_OP		timestep_bnd_node_last_touch_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		!*******************************
		!Geometry operators
		!*******************************

		subroutine heun_timestep_2_cell_last_touch_op(traversal, section, cell)
			type(t_grid_section), intent(inout)												:: grid

			type(t_cell_data_ptr), intent(inout)				:: cell

			call heun_timestep_2_post_dof(grid, cell%data_pers%T, cell%data_pers%T_temp, cell%data_temp%r, cell%data_temp%mat_mass_diagonal)
		end subroutine

		subroutine heun_timestep_2_edge_last_touch_op(traversal, section, edge)
			type(t_grid_section), intent(inout)												:: grid

			type(t_edge_data), intent(inout)		:: edge

			call heun_timestep_2_post_dof(grid, edge%data_pers%T, edge%data_pers%T_temp, edge%data_temp%r, edge%data_temp%mat_mass_diagonal)
		end subroutine

		subroutine heun_timestep_2_node_last_touch_op(traversal, section, node)
			type(t_grid_section), intent(inout)												:: grid

			type(t_node_data), intent(inout)				:: node

			call heun_timestep_2_post_dof(grid, node%data_pers%T, node%data_pers%T_temp, node%data_temp%r, node%data_temp%mat_mass_diagonal)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine heun_timestep_2_post_dof(grid, T, T_temp, r, mat_mass_diagonal)
			type(t_grid_section), intent(inout)						:: grid
			real(kind = GRID_SR), dimension(:), intent(out)	:: T
			real(kind = GRID_SR), dimension(:), intent(in)	:: T_temp
			real(kind = GRID_SR), dimension(:), intent(in)	:: r
			real(kind = GRID_SR), dimension(:), intent(in)	:: mat_mass_diagonal

			T = T_temp + 0.5_GRID_SR * grid%r_dt * r / mat_mass_diagonal
		end subroutine
	END MODULE

	MODULE Heat_Eq_Heun_Timestep
		use SFC_node_traversal

		use Samoa_heat_eq
		use Heat_Eq_Heun_Timestep_1
		use Heat_Eq_Heun_Timestep_2

		implicit none

		PUBLIC

		CONTAINS

		!> Wrapper method that calls the two traversal methods for the Heun timestep
		subroutine heun_timestep_traversal(grid)
			type(t_grid), intent(inout)							:: grid

			call heun_timestep_1_traversal(grid)
			grid%r_time = grid%r_time + grid%r_dt
			call heun_timestep_2_traversal(grid)
			grid%r_time = grid%r_time - grid%r_dt
		end subroutine
	END MODULE
#endif
