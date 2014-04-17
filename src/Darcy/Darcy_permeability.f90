! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_permeability
		use SFC_edge_traversal

		use Samoa_darcy

        type num_traversal_data
            integer (kind = GRID_DI)			:: i_refinements_issued
        end type

        type(darcy_gv_saturation)				:: gv_saturation
        type(darcy_gv_p)						:: gv_p

#		define _GT_NAME							t_darcy_permeability_traversal

#		if (_DARCY_FLOW_EDGE_SIZE > 0)
#			define _GT_EDGES
#		endif

#		define _GT_NODES
#		define _GT_REFINEMENTS

#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op
#		define _GT_ELEMENT_OP					element_op

#		include "SFC_generic_traversal_ringbuffer.f90"

 		subroutine post_traversal_grid_op(traversal, grid)
  			type(t_darcy_permeability_traversal), intent(inout)			:: traversal
 			type(t_grid), intent(inout)							        :: grid

			call reduce(traversal%i_refinements_issued, traversal%children%i_refinements_issued, MPI_SUM, .true.)
		end subroutine

		subroutine pre_traversal_op(traversal, section)
            type(t_darcy_permeability_traversal), intent(inout)         :: traversal
  			type(t_grid_section), intent(inout)							:: section

			!this variable will be incremented for each cell with a refinement request
			traversal%i_refinements_issued = 0
		end subroutine

		!*******************************
		!Geometry operators
		!*******************************

		subroutine element_op(traversal, section, element)
 			type(t_darcy_permeability_traversal), intent(inout)		:: traversal
 			type(t_grid_section), intent(inout)						:: section
			type(t_element_base), intent(inout), target				:: element

			real (kind = GRID_SR), dimension(_DARCY_FLOW_SIZE)		:: saturation
			real (kind = GRID_SR), dimension(_DARCY_P_SIZE)			:: p

			call gv_saturation%read(element, saturation)
			call gv_p%read(element, p)

			!call element operator
			call alpha_volume_op(traversal, section, element, saturation, p, element%cell%data_pers%base_permeability, element%cell%data_pers%permeability)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine alpha_volume_op(traversal, section, element, saturation, p, base_permeability, permeability)
            type(t_darcy_permeability_traversal)                                :: traversal
 			type(t_grid_section), intent(inout)							        :: section
			type(t_element_base), intent(inout)									:: element
			real (kind = GRID_SR), dimension(_DARCY_FLOW_SIZE), intent(in)		:: saturation
			real (kind = GRID_SR), dimension(_DARCY_P_SIZE), intent(in)			:: p
			real (kind = GRID_SR), intent(in)									:: base_permeability
			real (kind = GRID_SR), intent(out)									:: permeability

			real (kind = GRID_SR)												:: r_lambda_n, r_lambda_w
			real (kind = GRID_SR)												:: r_sat_norm, r_p_norm
			integer (kind = GRID_SI)											:: i_depth, i
			logical 											:: l_coarsen_p, l_coarsen_sat, l_refine_sat

			r_lambda_n = sum([0.25_GRID_SR, 0.5_GRID_SR, 0.25_GRID_SR] * (1.0_GRID_SR - saturation) * (1.0_GRID_SR - saturation))
			r_lambda_w = sum([0.25_GRID_SR, 0.5_GRID_SR, 0.25_GRID_SR] * saturation * saturation)

			permeability = base_permeability * (r_lambda_n + cfg%r_rel_permeability * r_lambda_w)

			r_sat_norm = max(abs(saturation(3) - saturation(2)), abs(saturation(1) - saturation(2)))
			r_p_norm = max(abs(p(3) - p(2)), abs(p(1) - p(2)))

			l_refine_sat = r_sat_norm > 0.1_GRID_SR
			l_coarsen_sat = r_sat_norm < 0.02_GRID_SR
			l_coarsen_p = r_p_norm < 0.003_GRID_SR * cfg%r_p0

			!* refine the cell if the saturation becomes too steep
			!* coarsen the cell if pressure and saturation are constant within a cell
			!* leave it as is otherwise

			i_depth = element%cell%geometry%i_depth

			if (i_depth < cfg%i_max_depth .and. base_permeability > 0.0_GRID_SR .and. (l_refine_sat .or. i_depth < cfg%i_min_depth)) then
				_log_write(5, "(A, T30, A, I0, A, L, A, ES9.2)") "  refinement issued:", "depth ", i_depth, ", sat ", l_refine_sat, ", perm ", base_permeability

				element%cell%geometry%refinement = 1
				traversal%i_refinements_issued = traversal%i_refinements_issued + 1_GRID_DI
			else if ((i_depth > cfg%i_min_depth .and. l_coarsen_p .and. l_coarsen_sat) .or. (base_permeability == 0.0_GRID_SR)) then
				_log_write(5, "(A, T30, A, I0, A, L, A, L, A, ES9.2)") "  coarsening issued:", "depth ", i_depth, ", p ", l_coarsen_p, ", sat ", l_coarsen_sat, ", perm ", base_permeability

				element%cell%geometry%refinement = -1
			else
				_log_write(5, "(A, T30, A, I0)") "  keep:", "depth ", i_depth

				element%cell%geometry%refinement = 0
			end if
		end subroutine
	END MODULE
#endif
