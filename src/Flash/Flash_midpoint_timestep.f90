! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_FLASH)
	MODULE FLASH_Midpoint_Timestep_1
		use SFC_edge_traversal

		use Samoa_FLASH
		use FLASH_euler_timestep

		implicit none

        	type num_traversal_data
            		integer (kind = GRID_SI)			:: i_refinements_issued
        	end type

		type(t_gv_Q)						:: gv_T
		type(t_lfs_flux)					:: lfs_flux

#		define _GT_NAME						midpoint_timestep_1_traversal
#		define _GT_EDGES

#		define _GT_PRE_TRAVERSAL_OP				timestep_pre_traversal_op

#		define _GT_POST_TRAVERSAL_GRID_OP			post_traversal_grid_op
#		define _GT_ELEMENT_OP					timestep_element_op

#		define _GT_CELL_FIRST_TOUCH_OP				timestep_cell_first_touch_op
#		define _GT_INNER_EDGE_FIRST_TOUCH_OP			timestep_edge_first_touch_op
#		define _GT_EDGE_FIRST_TOUCH_OP				timestep_bnd_edge_first_touch_op

#		define _GT_CELL_LAST_TOUCH_OP				midpoint_timestep_1_cell_last_touch_op
#		define _GT_INNER_EDGE_LAST_TOUCH_OP			midpoint_timestep_1_edge_last_touch_op
#		define _GT_EDGE_LAST_TOUCH_OP				timestep_bnd_edge_last_touch_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		!*******************************
		!Geometry operators
		!*******************************

		subroutine midpoint_timestep_1_cell_last_touch_op(traversal, section, cell)
 			type(midpoint_timestep_1_traversal), intent(inout)	                :: traversal
			type(t_grid_section), intent(inout)				:: section
			type(t_cell_data_ptr), intent(inout)				:: cell

			call midpoint_timestep_1_post_dof(section)
		end subroutine

		subroutine midpoint_timestep_1_edge_last_touch_op(traversal, section, edge)
 			type(midpoint_timestep_1_traversal), intent(inout)	                :: traversal			
			type(t_grid_section), intent(inout)		:: section
			type(t_edge_data), intent(inout)		:: edge

			call midpoint_timestep_1_post_dof(section)
		end subroutine



		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine midpoint_timestep_1_post_dof(grid)
			type(t_grid_section), intent(inout)					:: grid

		end subroutine
	END MODULE

	MODULE FLASH_Midpoint_Timestep_2
		use SFC_edge_traversal

		use Samoa_flash
		use FLASH_Euler_Timestep

		implicit none

        	type num_traversal_data
            		integer (kind = GRID_SI)			:: i_refinements_issued
        	end type

		type(t_gv_Q)						:: gv_T

#		define _GT_NAME							midpoint_timestep_2_traversal

#		if (_FLASH_EDGE_SIZE > 0)
#			define _GT_EDGES
#			define _GT_EDGES_TEMP
#		endif

#		define _GT_NODES
#		define _GT_NODES_TEMP
#		define _GT_REFINEMENTS

#		define _GT_PRE_TRAVERSAL_OP				timestep_pre_traversal_op

#		define _GT_ELEMENT_OP					timestep_element_op

#		define _GT_CELL_FIRST_TOUCH_OP			timestep_cell_first_touch_op
#		define _GT_INNER_EDGE_FIRST_TOUCH_OP		timestep_edge_first_touch_op
#		define _GT_EDGE_FIRST_TOUCH_OP			timestep_bnd_edge_first_touch_op

#		define _GT_CELL_LAST_TOUCH_OP			midpoint_timestep_2_cell_last_touch_op
#		define _GT_INNER_EDGE_LAST_TOUCH_OP		midpoint_timestep_2_edge_last_touch_op
#		define _GT_EDGE_LAST_TOUCH_OP			timestep_bnd_edge_last_touch_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		!*******************************
		!Geometry operators
		!*******************************

		subroutine midpoint_timestep_2_cell_last_touch_op(traversal, section, cell)
 			type(midpoint_timestep_2_traversal), intent(inout)	                :: traversal
			type(t_grid_section), intent(inout)				:: section
			type(t_cell_data_ptr), intent(inout)				:: cell

			call midpoint_timestep_2_post_dof(section)
		end subroutine

		subroutine midpoint_timestep_2_edge_last_touch_op(traversal, section, edge)
 			type(midpoint_timestep_2_traversal), intent(inout)	:: traversal
			type(t_grid_section), intent(inout)		:: section

			type(t_edge_data), intent(inout)		:: edge

			call midpoint_timestep_2_post_dof(section)
		end subroutine

		subroutine midpoint_timestep_2_node_last_touch_op(traversal, section, node)
 			type(midpoint_timestep_2_traversal), intent(inout)	                :: traversal
			type(t_grid_section), intent(inout)				:: section
			type(t_node_data), intent(inout)				:: node

			call midpoint_timestep_2_post_dof(section)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine midpoint_timestep_2_post_dof(section)
			type(t_grid_section), intent(inout)		:: section

			!T = T_temp + grid%r_dt * r / mat_mass_diagonal
		end subroutine
	END MODULE

	MODULE FLASH_Midpoint_Timestep
		use SFC_edge_traversal

		use Samoa_flash
		use FLASH_Midpoint_Timestep_1
		use FLASH_Midpoint_Timestep_2

		implicit none

		PUBLIC

		CONTAINS

		!> Wrapper method that calls the two traversal methods for the Midpoint timestep
		subroutine midpoint_timestep_traversal(grid)
			type(t_grid), intent(inout)							:: grid

			!call midpoint_timestep_1_traversal(grid)
			grid%r_time = grid%r_time + 0.5_GRID_SR * grid%r_dt
			!call midpoint_timestep_2_traversal(grid)
			grid%r_time = grid%r_time - 0.5_GRID_SR * grid%r_dt
		end subroutine
	END MODULE
#endif
