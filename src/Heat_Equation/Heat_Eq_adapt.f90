! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_HEAT_EQ)
	MODULE Heat_Eq_Adapt
		use SFC_edge_traversal
        use Conformity

		use Samoa_heat_eq
		use Tools_noise

		implicit none

        type num_traversal_data
			real (kind = GRID_SR), dimension(_HEAT_EQ_SIZE, 2)				    :: T_in
        end type

		type(heat_eq_gv_T)							:: gv_T

#		define	_GT_NAME							heat_eq_adaption_traversal

#		if (_HEAT_EQ_EDGE_SIZE > 0)
#			define _GT_EDGES
#		endif

#		define	_GT_NODES

#		define	_GT_TRANSFER_OP						transfer_op
#		define	_GT_REFINE_OP						refine_op
#		define	_GT_COARSEN_OP						coarsen_op

#		include "SFC_generic_adaptive_traversal.f90"

		!******************
		!Geometry operators
		!******************
		subroutine transfer_op(traversal, section, src_element, dest_element)
			type(t_grid_section), intent(inout)													:: grid
			type(t_traversal_element), intent(inout)									:: src_element
			type(t_traversal_element), intent(inout)									:: dest_element

			real (kind = GRID_SR), dimension(_HEAT_EQ_SIZE)								:: T

			call gv_T%read( src_element%t_element_base, T)
			call gv_T%write( dest_element%t_element_base, T)

			dest_element%cell%data_pers%heat_conductivity = src_element%cell%data_pers%heat_conductivity
		end subroutine

		subroutine refine_op(traversal, section, src_element, dest_element, refinement_path)
			type(t_grid_section), intent(inout)												:: grid
			type(t_traversal_element), intent(inout)								:: src_element
			type(t_traversal_element), intent(inout)								:: dest_element
			integer, dimension(:), intent(in)										:: refinement_path

			real (kind = GRID_SR), dimension(_HEAT_EQ_SIZE)							:: T_in
			real (kind = GRID_SR), dimension(_HEAT_EQ_SIZE, 2)						:: T_out
			integer																	:: i

			!temperature

			call gv_T%read( src_element%t_element_base, T_in)

			do i = 1, size(refinement_path)
				call samoa_basis_T_split(T_in, T_out(:, 1), T_out(:, 2))

				T_in = T_out(:, refinement_path(i))
			end do

			call gv_T%write( dest_element%t_element_base, T_in)

			!heat conductivity
			dest_element%cell%data_pers%heat_conductivity = get_heat_conductivity(grid, samoa_barycentric_to_world_point(dest_element%transform_data, [1.0_GRID_SR / 3.0_GRID_SR, 1.0_GRID_SR / 3.0_GRID_SR]), dest_element%cell%geometry%i_depth / 2_GRID_SI)
		end subroutine

		subroutine coarsen_op(traversal, section, src_element, dest_element, refinement_path)
			type(t_grid_section), intent(inout)												:: grid
			type(t_traversal_element), intent(inout)								:: src_element
			type(t_traversal_element), intent(inout)								:: dest_element
			integer, dimension(:), intent(in)										:: refinement_path

			real (kind = GRID_SR), dimension(_HEAT_EQ_SIZE)							:: T_out
			integer																	:: i

			!temperature
			i = refinement_path(1)
			call gv_T%read( src_element%t_element_base, T_in(:, i))

			if (i > 1) then
				call samoa_basis_T_merge(T_in(:, 1), T_in(:, 2), T_out)
				call gv_T%write( dest_element%t_element_base, T_out)

				!heat conductivity (compute the average)
				dest_element%cell%data_pers%heat_conductivity = dest_element%cell%data_pers%heat_conductivity + 0.5_GRID_SR * src_element%cell%data_pers%heat_conductivity
			else
				dest_element%cell%data_pers%heat_conductivity = 0.5_GRID_SR * src_element%cell%data_pers%heat_conductivity
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
