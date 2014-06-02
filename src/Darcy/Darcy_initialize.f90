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

			if (all(node%position == [0.0_GRID_SR, 0.0_GRID_SR]) .or. all(node%position == [1.0_GRID_SR, 1.0_GRID_SR])) then
				node%data_temp%is_dirichlet_boundary = .true.
			else
				node%data_temp%is_dirichlet_boundary = .false.
			end if

            call inner_node_first_touch_op(traversal, section, node)
		end subroutine

		elemental subroutine inner_node_first_touch_op(traversal, section, node)
 			type(t_darcy_init_pressure_traversal), intent(in)                       :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			:: node
			integer										:: i

			do i = 1, _DARCY_P_NODE_SIZE
				call pressure_pre_dof_op(real(cfg%r_p_in, GRID_SR), real(cfg%r_p_prod, GRID_SR), node%position, node%data_pers%p(i), node%data_pers%r(i), node%data_pers%d(i), node%data_pers%A_d(i))
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
			real (kind = GRID_SR)								:: r_base_permeability		!< permeability tensor

            real (kind = GRID_SR)                               :: xs(2)

            xs = cfg%scaling * x + cfg%offset

            assert_ge(x(1), 0.0); assert_ge(x(2), 0.0)
            assert_le(x(1), 1.0); assert_le(x(2), 1.0)

#			if defined(_ASAGI)
#               if defined(_ASAGI_TIMING)
                    section%stats%r_asagi_time = section%stats%r_asagi_time - get_wtime()
#               endif
                r_base_permeability = 1.0e-8_GRID_SR + asagi_get_float(cfg%afh_permeability, dble(xs(1)), dble(xs(2)), 0)
#               if defined(_ASAGI_TIMING)
                    section%stats%r_asagi_time = section%stats%r_asagi_time + get_wtime()
#               endif
#			else
				r_base_permeability = -7.0e-8_GRID_SR * (t_noise_2D((/ 10.0_GRID_SR * xs(1) - 2.0_GRID_SR, 10.0_GRID_SR * xs(2) /), lod, 0.2_GRID_SR) + 0.7_GRID_SR - 4.0_GRID_SR * xs(2) * (1.0_GRID_SR - xs(2))) + 0.5e-8_GRID_SR
				r_base_permeability = max(0.0_GRID_SR, min(1.0e-8_GRID_SR, r_base_permeability))
#			endif

			!r_base_permeability = 1.0e-8
		end function

		pure subroutine pressure_pre_dof_op(p_in, p_prod, pos, p, r, d, A_d)
 			real (kind = GRID_SR), intent(in)					:: p_in
 			real (kind = GRID_SR), intent(in)					:: p_prod
			real (kind = GRID_SR), intent(in)					:: pos(2)
			real (kind = GRID_SR), intent(out)					:: p
			real (kind = GRID_SR), intent(out)					:: r
			real (kind = GRID_SR), intent(out)					:: d
			real (kind = GRID_SR), intent(out)					:: A_d

			p = 0.5_GRID_SR * (2.0_GRID_SR - pos(1) - pos(2)) * p_in + 0.5_GRID_SR * (pos(1) + pos(2)) * p_prod
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
		type(darcy_gv_rhs)						:: gv_rhs

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

			call inner_flow_pre_dof_op(node%data_pers%saturation, node%data_pers%rhs)
		end subroutine

		elemental subroutine node_first_touch_op(traversal, section, node)
 			type(t_darcy_init_saturation_traversal), intent(in)                     :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			:: node

			integer :: i

			do i = 1, _DARCY_FLOW_NODE_SIZE
				call flow_pre_dof_op(node%position, node%data_pers%saturation(i), node%data_pers%rhs(i))
			end do
		end subroutine

		subroutine element_op(traversal, section, element)
 			type(t_darcy_init_saturation_traversal)                 :: traversal
 			type(t_grid_section), intent(inout)						:: section
			type(t_element_base), intent(inout), target				:: element

			real (kind = GRID_SR), dimension(_DARCY_FLOW_SIZE)		:: saturation
			real (kind = GRID_SR), dimension(_DARCY_P_SIZE)			:: p
			real (kind = GRID_SR)		                            :: rhs(_DARCY_P_SIZE)

			call gv_saturation%read(element, saturation)
			call gv_p%read(element, p)

			!call element operator
			call alpha_volume_op(traversal, section, element, saturation, p, element%cell%data_pers%base_permeability, element%cell%data_pers%permeability, rhs)

			call gv_rhs%add(element, rhs)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		elemental subroutine inner_flow_pre_dof_op(saturation, rhs)
			real (kind = GRID_SR), intent(out)					:: saturation
			real (kind = GRID_SR), intent(out)					:: rhs

			saturation = 0.0_GRID_SR
			rhs = 0.0_GRID_SR
		end subroutine

		pure subroutine flow_pre_dof_op(pos, saturation, rhs)
			real (kind = GRID_SR), dimension(2), intent(in)		:: pos
			real (kind = GRID_SR), intent(out)					:: saturation
			real (kind = GRID_SR), intent(out)					:: rhs

			if (pos(1) == 0.0_GRID_SR .and. pos(2) == 0.0_GRID_SR) then
				saturation = 1.0_GRID_SR
			else
				saturation = 0.0_GRID_SR
			end if

			rhs = 0.0_GRID_SR
		end subroutine

		subroutine alpha_volume_op(traversal, section, element, saturation, p, base_permeability, permeability, rhs)
 			type(t_darcy_init_saturation_traversal)                 :: traversal
 			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)									:: element
			real (kind = GRID_SR), dimension(_DARCY_FLOW_SIZE), intent(in)		:: saturation
			real (kind = GRID_SR), dimension(_DARCY_P_SIZE), intent(in)			:: p
			real (kind = GRID_SR), intent(in)									:: base_permeability
			real (kind = GRID_SR), intent(out)									:: permeability
			real (kind = GRID_SR), intent(out)									:: rhs(:)

			integer (kind = GRID_SI)											:: i_depth
			logical 											                :: l_refine_p, l_coarsen_p, l_coarsen_sat, l_refine_sat
			real (kind = GRID_SR), parameter                                    :: g(2) = [0.0d0, -9.81d0]
			real (kind = GRID_SR)												:: x(2), grad_psi(2), r_lambda_w, r_lambda_n
			integer     :: i

			!compute permeability

			r_lambda_w = sum([0.25_GRID_SR, 0.5_GRID_SR, 0.25_GRID_SR] * saturation * saturation) / cfg%r_nu_w
			r_lambda_n = sum([0.25_GRID_SR, 0.5_GRID_SR, 0.25_GRID_SR] * (1.0_GRID_SR - saturation) * (1.0_GRID_SR - saturation)) / cfg%r_nu_n

			permeability = base_permeability * (r_lambda_w + r_lambda_n)

            do i = 1, 3
                x = samoa_basis_p_get_dof_coords(i)
                grad_psi = [samoa_basis_p_d_dx(x), samoa_basis_p_d_dy(x)]
                grad_psi = samoa_barycentric_to_world_normal(element%transform_data, grad_psi)
                rhs(i) = (cfg%r_rho_w * r_lambda_w + cfg%r_rho_n * r_lambda_n) * dot_product(g, grad_psi) !this should be the integral over lambda_t * g * grad psi
            end do

			l_refine_sat = max(abs(saturation(3) - saturation(2)), abs(saturation(1) - saturation(2))) > 0.1_GRID_SR
			l_refine_p = max(abs(p(3) - p(2)), abs(p(1) - p(2))) > 0.01_GRID_SR * (cfg%r_p_in - cfg%r_p_prod)

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
