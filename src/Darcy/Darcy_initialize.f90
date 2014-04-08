! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_initialize_pressure
		use SFC_edge_traversal

		use Samoa_darcy
		use Tools_noise

		implicit none

        type num_traversal_data
        end type

		type(darcy_gv_p)						:: gv_p
		type(darcy_gv_saturation)				:: gv_saturation

		public get_base_permeability

#		define	_GT_NAME						t_darcy_init_pressure_traversal

#		if (_DARCY_P_EDGE_SIZE > 0)
#			define _GT_EDGES
#		endif

#		define _GT_NODES

#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_NODE_FIRST_TOUCH_OP		    node_first_touch_op
#		define _GT_INNER_NODE_FIRST_TOUCH_OP    inner_node_first_touch_op

#		define	_GT_ELEMENT_OP					element_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine pre_traversal_grid_op(traversal, grid)
 			type(t_darcy_init_pressure_traversal), intent(inout)      	:: traversal
 			type(t_grid), intent(inout)							        :: grid

			call scatter(grid%r_time, grid%sections%elements_alloc%r_time)
		end subroutine

		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, section, element)
 			type(t_darcy_init_pressure_traversal)                       :: traversal
 			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)											:: element

			real (kind = GRID_SR), dimension(_DARCY_P_SIZE)								:: p

			call alpha_volume_op(section, element, element%cell%data_pers%base_permeability)
		end subroutine

		elemental subroutine node_first_touch_op(traversal, section, node)
 			type(t_darcy_init_pressure_traversal), intent(in)                       :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			:: node

			if (node%position(1) > 0.0_GRID_SR .and. node%position(1) < 1.0_GRID_SR) then
				node%data_temp%is_dirichlet_boundary = .false.
			else
				node%data_temp%is_dirichlet_boundary = .true.
			end if

            call inner_node_first_touch_op(traversal, section, node)
		end subroutine

		elemental subroutine inner_node_first_touch_op(traversal, section, node)
 			type(t_darcy_init_pressure_traversal), intent(in)                       :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			:: node
			integer										:: i

			do i = 1, _DARCY_P_NODE_SIZE
				call pressure_pre_dof_op(real(cfg%r_p0, GRID_SR), node%position, node%data_pers%p(i), node%data_pers%r(i), node%data_pers%d(i), node%data_pers%A_d(i))
			end do
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine alpha_volume_op(section, element, base_permeability)
 			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(in)									:: element
			real (kind = GRID_SR), intent(out)									:: base_permeability

			real (kind = GRID_SR), dimension(2)									:: pos

			!set base permeability
			base_permeability = get_base_permeability(section, samoa_barycentric_to_world_point(element%transform_data, (/ 1.0_GRID_SR / 3.0_GRID_SR, 1.0_GRID_SR / 3.0_GRID_SR /)), element%cell%geometry%i_depth / 2_GRID_SI)
		end subroutine

		function get_base_permeability(section, x, lod) result(r_base_permeability)
 			type(t_grid_section), intent(inout)					:: section
			real (kind = GRID_SR), dimension(:), intent(in)		:: x						!< position in world coordinates
			integer (kind = GRID_SI), intent(in)				:: lod						!< level of detail
			real (kind = GRID_SR)								:: r_base_permeability		!< heat conductivity

            real (kind = GRID_SR)                               :: xs(2)

            xs = cfg%scaling * x + cfg%offset

            assert_ge(x(1), 0.0); assert_ge(x(2), 0.0)
            assert_le(x(1), 1.0); assert_le(x(2), 1.0)

#			if defined(_ASAGI)

#               if defined(_ASAGI_TIMING)
                    section%stats%r_asagi_time = section%stats%r_asagi_time - omp_get_wtime()
#               endif

#               if defined(_ASAGI_NUMA)
                	r_base_permeability = asagi_get_float(cfg%afh_permeability, dble(xs(1)), dble(xs(2)), 0)
#				else
                	r_base_permeability = asagi_get_float(cfg%afh_permeability, dble(xs(1)), dble(xs(2)), lod)
#				endif

#               if defined(_ASAGI_TIMING)
                    section%stats%r_asagi_time = section%stats%r_asagi_time + omp_get_wtime()
#               endif
#			else
				r_base_permeability = -7.0e-8_GRID_SR * (t_noise_2D((/ 10.0_GRID_SR * xs(1) - 2.0_GRID_SR, 10.0_GRID_SR * xs(2) /), lod, 0.2_GRID_SR) + 0.7_GRID_SR - 4.0_GRID_SR * xs(2) * (1.0_GRID_SR - xs(2))) + 0.5e-8_GRID_SR
				r_base_permeability = max(0.0_GRID_SR, min(1.0e-8_GRID_SR, r_base_permeability))
#			endif
		end function

		pure subroutine pressure_pre_dof_op(p0, pos, p, r, d, A_d)
 			real (kind = GRID_SR), intent(in)					:: p0
			real (kind = GRID_SR), intent(in)					:: pos(2)
			real (kind = GRID_SR), intent(out)					:: p
			real (kind = GRID_SR), intent(out)					:: r
			real (kind = GRID_SR), intent(out)					:: d
			real (kind = GRID_SR), intent(out)					:: A_d

			p = (1.0_GRID_SR - pos(1)) * p0
			r = 0.0_GRID_SR
			d = 0.0_GRID_SR
			A_d = 0.0_GRID_SR
		end subroutine
	END MODULE

	MODULE Darcy_initialize_saturation
		use SFC_edge_traversal
		use Samoa_darcy

		implicit none

        type num_traversal_data
            integer (kind = GRID_DI)			:: i_refinements_issued
        end type

		type(darcy_gv_saturation)				:: gv_saturation
		type(darcy_gv_p)						:: gv_p

#		define	_GT_NAME						t_darcy_init_saturation_traversal

#		if (_DARCY_P_EDGE_SIZE > 0 || _DARCY_FLOW_EDGE_SIZE > 0)
#			define _GT_EDGES
#		endif

#		define _GT_NODES
#		define _GT_REFINEMENTS

#		define _GT_INNER_NODE_FIRST_TOUCH_OP	inner_node_first_touch_op
#		define _GT_NODE_FIRST_TOUCH_OP		    node_first_touch_op

#		define	_GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define	_GT_PRE_TRAVERSAL_OP			pre_traversal_op
#		define	_GT_ELEMENT_OP					element_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine post_traversal_grid_op(traversal, grid)
 			type(t_darcy_init_saturation_traversal), intent(inout)      :: traversal
 			type(t_grid), intent(inout)							        :: grid

			call reduce(traversal%i_refinements_issued, traversal%children%i_refinements_issued, MPI_SUM, .true.)
		end subroutine

		subroutine pre_traversal_op(traversal, section)
 			type(t_darcy_init_saturation_traversal), intent(inout)      :: traversal
 			type(t_grid_section), intent(inout)							:: section

			!this variable will be incremented for each cell with a refinement request
			traversal%i_refinements_issued = 0
		end subroutine

		!******************
		!Geometry operators
		!******************

		elemental subroutine inner_node_first_touch_op(traversal, section, node)
 			type(t_darcy_init_saturation_traversal), intent(in)                     :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			:: node

			call inner_flow_pre_dof_op(node%data_pers%saturation)
		end subroutine

		elemental subroutine node_first_touch_op(traversal, section, node)
 			type(t_darcy_init_saturation_traversal), intent(in)                     :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			:: node

			integer :: i

			do i = 1, _DARCY_FLOW_NODE_SIZE
				call flow_pre_dof_op(node%position, node%data_pers%saturation(i))
			end do
		end subroutine

		subroutine element_op(traversal, section, element)
 			type(t_darcy_init_saturation_traversal)                 :: traversal
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

		elemental subroutine inner_flow_pre_dof_op(saturation)
			real (kind = GRID_SR), intent(out)					:: saturation

			saturation = 0.0_GRID_SR
		end subroutine

		pure subroutine flow_pre_dof_op(pos, saturation)
			real (kind = GRID_SR), dimension(2), intent(in)		:: pos
			real (kind = GRID_SR), intent(out)					:: saturation

			if (pos(1) == 0.0_GRID_SR .and. pos(2) > 3.1_GRID_SR / 8.0_GRID_SR .and. pos(2) < 4.9_GRID_SR/8.0_GRID_SR) then
				saturation = 1.0_GRID_SR
			else
				saturation = 0.0_GRID_SR
			end if
		end subroutine

		subroutine alpha_volume_op(traversal, section, element, saturation, p, base_permeability, permeability)
 			type(t_darcy_init_saturation_traversal)                 :: traversal
 			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)									:: element
			real (kind = GRID_SR), dimension(_DARCY_FLOW_SIZE), intent(in)		:: saturation
			real (kind = GRID_SR), dimension(_DARCY_P_SIZE), intent(in)			:: p
			real (kind = GRID_SR), intent(in)									:: base_permeability
			real (kind = GRID_SR), intent(out)									:: permeability

			real (kind = GRID_SR)												:: r_lambda_n, r_lambda_w
			integer (kind = GRID_SI)											:: i_depth
			logical (kind = GRID_SL)											:: l_refine_p, l_coarsen_p, l_coarsen_sat, l_refine_sat

			!compute permeability

			r_lambda_n = sum([0.25_GRID_SR, 0.5_GRID_SR, 0.25_GRID_SR] * (1.0_GRID_SR - saturation) * (1.0_GRID_SR - saturation))
			r_lambda_w = sum([0.25_GRID_SR, 0.5_GRID_SR, 0.25_GRID_SR] * saturation * saturation)

			permeability = base_permeability * (r_lambda_n + cfg%r_rel_permeability * r_lambda_w)

			l_refine_sat = max(abs(saturation(3) - saturation(2)), abs(saturation(1) - saturation(2))) > 0.1_GRID_SR
			l_refine_p = max(abs(p(3) - p(2)), abs(p(1) - p(2))) > 0.01_GRID_SR * cfg%r_p0

			!refine the cell if necessary (no coarsening in the initialization!)

			i_depth = element%cell%geometry%i_depth

			if (i_depth < cfg%i_max_depth .and. base_permeability > 0.0_GRID_SR .and. (i_depth < cfg%i_min_depth .or. l_refine_p .or. l_refine_sat)) then
				element%cell%geometry%refinement = 1
				traversal%i_refinements_issued = traversal%i_refinements_issued + 1_GRID_DI
			else
				element%cell%geometry%refinement = 0
			end if
		end subroutine
	END MODULE
#endif
