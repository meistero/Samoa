! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_HEAT_EQ)
	MODULE Heat_Eq_Initialize
		use Tools_noise

		use SFC_node_traversal

		use Samoa_heat_eq

		integer (kind = GRID_DI)												:: i_refinements_issued

		type(heat_eq_gv_T)														:: gv_T

		PUBLIC heat_eq_init_count_issued_refinements

#		define	_GT_NAME							heat_eq_init_traversal

#		if (_HEAT_EQ_EDGE_SIZE > 0)
#			define _GT_EDGES
#			define _GT_EDGES_TEMP
#		endif

#		define _GT_NODES
#		define _GT_NODES_TEMP
#		define _GT_REFINEMENTS

#		define	_GT_ELEMENT_OP						element_op

#		define	_GT_PRE_TRAVERSAL_OP				pre_traversal_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		function heat_eq_init_count_issued_refinements()
			integer (kind = GRID_SI)			:: heat_eq_init_count_issued_refinements

			heat_eq_init_count_issued_refinements = i_refinements_issued
		end function

		subroutine pre_traversal_op(traversal, section)
 			type(t_grid_section), intent(inout)							:: grid

			!this variable will be incremented for each cell with a refinement request
			i_refinements_issued = 0
		end subroutine

		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, section, element)
			type(t_grid_section), intent(inout)												:: grid
			type(t_element_base), intent(inout)					:: element

			real (kind = GRID_SR), dimension(_HEAT_EQ_SIZE)		:: T

			call alpha_volume_op(traversal, section, element, T, element%cell%data_pers%heat_conductivity)

			call gv_T%write(element, T)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine alpha_volume_op(traversal, section, element, T, heat_conductivity)
			type(t_grid_section), intent(inout)												:: grid
			type(t_element_base), intent(inout)						:: element
			real (kind = GRID_SR), dimension(:), intent(out)		:: T
			real (kind = GRID_SR), intent(out)						:: heat_conductivity

			!evaluate initial function values at dof positions and compute DoFs
			T = 0.0_GRID_SR

			heat_conductivity = get_heat_conductivity(grid, samoa_barycentric_to_world_point(element%transform_data, [1.0_GRID_SR / 3.0_GRID_SR, 1.0_GRID_SR / 3.0_GRID_SR]), element%cell%geometry%i_depth / 2_GRID_SI)

			!refine the cell if any of these conditions hold:
			! * not yet at min depth and in fluid phase

			if (element%cell%geometry%i_depth < cfg%i_min_depth .and. heat_conductivity > 0.0_GRID_SR) then
				element%cell%geometry%refinement = 1
				i_refinements_issued = i_refinements_issued + 1
			else
				element%cell%geometry%refinement = 0
			end if
		end subroutine

		function get_heat_conductivity(grid, x, lod) result(r_heat_conductivity)
 			type(t_grid_section), intent(inout)					:: grid
			real (kind = GRID_SR), dimension(:), intent(in)		:: x						!< position in world coordinates
			integer (kind = GRID_SI),intent(in)					:: lod						!< level of detail
			real (kind = GRID_SR)								:: r_heat_conductivity		!< heat conductivity

            assert_ge(x(1), 0.0); assert_ge(x(2), 0.0)
            assert_le(x(1), 1.0); assert_le(x(2), 1.0)

			r_heat_conductivity = -5.0_GRID_SR * (t_noise_2D(8.0_GRID_SR * x, lod, 0.2_GRID_SR) - 0.2_GRID_SR) + 0.5_GRID_SR
			r_heat_conductivity = max(0.0_GRID_SR, min(0.01_GRID_SR, r_heat_conductivity))
		end function
	END MODULE
#endif
