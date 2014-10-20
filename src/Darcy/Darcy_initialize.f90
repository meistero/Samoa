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

		public get_base_permeability, get_porosity

#		define	_GT_NAME						t_darcy_init_pressure_traversal

#		if (_DARCY_P_EDGE_SIZE > 0)
#			define _GT_EDGES
#		endif

#		define _GT_NODES

#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_NODE_FIRST_TOUCH_OP		    node_first_touch_op
#		define _GT_INNER_NODE_FIRST_TOUCH_OP	inner_node_first_touch_op

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

			call alpha_volume_op(section, element, element%cell%data_pers%base_permeability, element%cell%data_pers%porosity)
		end subroutine

		subroutine node_first_touch_op(traversal, section, nodes)
 			type(t_darcy_init_pressure_traversal), intent(in)   :: traversal
 			type(t_grid_section), intent(inout)				    :: section
			type(t_node_data), intent(inout)			        :: nodes(:)

            integer :: i

            do i = 1, size(nodes)
                call inner_node_first_touch_op(traversal, section, nodes(i))
            end do
		end subroutine

		subroutine inner_node_first_touch_op(traversal, section, node)
 			type(t_darcy_init_pressure_traversal), intent(in)                       :: traversal
 			type(t_grid_section), intent(inout)		    :: section
			type(t_node_data), intent(inout)			:: node
			integer										:: i

			call pressure_pre_dof_op(real(cfg%r_p_prod, GRID_SR), node%data_pers%p, node%data_pers%r, node%data_pers%d, node%data_pers%A_d)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine alpha_volume_op(section, element, base_permeability, porosity)
 			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(in)									:: element
			real (kind = GRID_SR), intent(out)									:: base_permeability
			real (kind = GRID_SR), intent(out)									:: porosity

			base_permeability = get_base_permeability(section, samoa_barycentric_to_world_point(element%transform_data, [1.0_SR / 3.0_SR, 1.0_SR / 3.0_SR]), element%cell%geometry%i_depth / 2_GRID_SI)
			porosity = get_porosity(section, samoa_barycentric_to_world_point(element%transform_data, [1.0_SR / 3.0_SR, 1.0_SR / 3.0_SR]))
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

#           if defined(_ASAGI_TIMING)
                section%stats%r_asagi_time = section%stats%r_asagi_time - get_wtime()
#           endif

            if (grid_min_x(cfg%afh_permeability) <= xs(1) .and. grid_min_y(cfg%afh_permeability) <= xs(2) &
                    .and. xs(1) <= grid_max_x(cfg%afh_permeability) .and. xs(2) <= grid_max_y(cfg%afh_permeability)) then

                r_base_permeability = grid_get_float_3d(cfg%afh_permeability, dble(xs(1)), dble(xs(2)), 0.1_SR, 0)

                !convert from mD to m^2
                r_base_permeability =  9.869233e-16_SR * r_base_permeability
            else
                r_base_permeability = 0.0_SR
            end if

#           if defined(_ASAGI_TIMING)
                section%stats%r_asagi_time = section%stats%r_asagi_time + get_wtime()
#           endif
		end function

		function get_porosity(section, x) result(porosity)
 			type(t_grid_section), intent(inout)					:: section
			real (kind = GRID_SR), dimension(:), intent(in)		:: x						!< position in world coordinates
			real (kind = GRID_SR)								:: porosity		        !< porosity

            real (kind = GRID_SR)                               :: xs(2)

            xs = cfg%scaling * x + cfg%offset

            assert_ge(x(1), 0.0); assert_ge(x(2), 0.0)
            assert_le(x(1), 1.0); assert_le(x(2), 1.0)

#           if defined(_ASAGI_TIMING)
                section%stats%r_asagi_time = section%stats%r_asagi_time - get_wtime()
#           endif

            if (grid_min_x(cfg%afh_porosity) <= xs(1) .and. grid_min_y(cfg%afh_porosity) <= xs(2) &
                    .and. xs(1) <= grid_max_x(cfg%afh_porosity) .and. xs(2) <= grid_max_y(cfg%afh_porosity)) then

                porosity = max(0.0d0, grid_get_float_3d(cfg%afh_porosity, dble(xs(1)), dble(xs(2)), 0.1_SR, 0))
            else
                porosity = 0.0d0
            end if

#           if defined(_ASAGI_TIMING)
                section%stats%r_asagi_time = section%stats%r_asagi_time + get_wtime()
#           endif
		end function

		elemental subroutine pressure_pre_dof_op(p_prod, p, r, d, A_d)
 			real (kind = GRID_SR), intent(in)					:: p_prod
			real (kind = GRID_SR), intent(out)					:: p
			real (kind = GRID_SR), intent(out)					:: r
			real (kind = GRID_SR), intent(out)					:: d
			real (kind = GRID_SR), intent(out)					:: A_d

			p = p_prod
			r = 0.0_SR
			d = 0.0_SR
			A_d = 0.0_SR
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
		type(darcy_gv_is_dirichlet_boundary)    :: gv_is_dirichlet

#		define	_GT_NAME						t_darcy_init_saturation_traversal

#		if (_DARCY_P_EDGE_SIZE > 0 || _DARCY_FLOW_EDGE_SIZE > 0)
#			define _GT_EDGES
#		endif

#		define _GT_NODES

#		define _GT_INNER_NODE_FIRST_TOUCH_OP	inner_node_first_touch_op
#		define _GT_NODE_FIRST_TOUCH_OP		    node_first_touch_op

#		define	_GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define	_GT_PRE_TRAVERSAL_OP			pre_traversal_op
#		define	_GT_ELEMENT_OP					element_op

#		define _GT_NODE_MERGE_OP		        node_merge_op

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

			call flow_pre_dof_op(node%data_pers%saturation, node%data_pers%rhs, node%data_temp%is_dirichlet_boundary)
		end subroutine

		elemental subroutine node_first_touch_op(traversal, section, node)
 			type(t_darcy_init_saturation_traversal), intent(in)                     :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			:: node

			call flow_pre_dof_op(node%data_pers%saturation, node%data_pers%rhs, node%data_temp%is_dirichlet_boundary)
		end subroutine

		subroutine element_op(traversal, section, element)
 			type(t_darcy_init_saturation_traversal)                 :: traversal
 			type(t_grid_section), intent(inout)						:: section
			type(t_element_base), intent(inout), target				:: element

			real (kind = GRID_SR)   :: saturation(_DARCY_FLOW_SIZE)
			real (kind = GRID_SR)   :: p(_DARCY_P_SIZE)
			real (kind = GRID_SR)   :: rhs(_DARCY_P_SIZE)
			logical		            :: is_dirichlet(_DARCY_P_SIZE)

			call gv_saturation%read(element, saturation)
			call gv_p%read(element, p)

			!call element operator
			call alpha_volume_op(traversal, element, saturation, p, rhs, is_dirichlet, element%cell%data_pers%base_permeability, element%cell%data_pers%permeability)

			call gv_rhs%add(element, rhs)
		end subroutine

		elemental subroutine node_merge_op(local_node, neighbor_node)
			type(t_node_data), intent(inout)		:: local_node
 			type(t_node_data), intent(in)		    :: neighbor_node

            if (any(neighbor_node%data_temp%is_dirichlet_boundary)) then
                local_node%data_pers%p = neighbor_node%data_pers%p
            end if

			local_node%data_pers%saturation = max(local_node%data_pers%saturation, neighbor_node%data_pers%saturation)
			local_node%data_temp%is_dirichlet_boundary = local_node%data_temp%is_dirichlet_boundary .or. neighbor_node%data_temp%is_dirichlet_boundary
			local_node%data_pers%rhs = local_node%data_pers%rhs + neighbor_node%data_pers%rhs
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		elemental subroutine flow_pre_dof_op(saturation, rhs, is_dirichlet)
			real (kind = GRID_SR), intent(out)					:: saturation
			real (kind = GRID_SR), intent(out)					:: rhs
			logical, intent(out)			                    :: is_dirichlet

			saturation = 0.0_SR
			rhs = 0.0_SR
			is_dirichlet = .false.
		end subroutine

		subroutine alpha_volume_op(traversal, element, saturation, p, rhs, is_dirichlet, base_permeability, permeability)
 			type(t_darcy_init_saturation_traversal)                             :: traversal
			type(t_element_base), intent(inout)									:: element
			real (kind = GRID_SR), intent(inout)		                        :: saturation(:)
			real (kind = GRID_SR), intent(inout)			                    :: p(:)
			real (kind = GRID_SR), intent(out)									:: rhs(:)
			logical, intent(out)			                                    :: is_dirichlet(:)
			real (kind = GRID_SR), intent(in)									:: base_permeability
			real (kind = GRID_SR), intent(out)									:: permeability

			real (kind = GRID_SR), parameter            :: dx(3) = [1.0_SR, -1.0_SR, 0.0_SR]
			real (kind = GRID_SR), parameter            :: dy(3) = [0.0_SR, -1.0_SR, 1.0_SR]
			real (kind = GRID_SR), parameter            :: volumes(3) = [1.0_SR/8.0_SR, 1.0_SR/4.0_SR, 1.0_SR/8.0_SR]

			real (kind = GRID_SR)					    :: g_local(2), pos_prod(2), pos_in(2), r_lambda_w(3), r_lambda_n(3), radius
			integer (kind = GRID_SI)	                :: i_depth
			logical 								    :: l_refine_initial, l_refine_solution

            !set boundary conditions and source terms

            rhs(:) = 0.0_SR
            radius = 0.0508_SR / (cfg%scaling * element%transform_data%custom_data%scaling)

            l_refine_initial = .false.

            pos_prod = sign(cfg%r_pos_prod - 0.5_SR, element%transform_data%custom_data%offset - 0.5_SR) + 0.5_SR
            pos_prod = samoa_world_to_barycentric_point(element%transform_data, pos_prod)
            if (pos_prod(1) .ge. -radius .and. pos_prod(2) .ge. -radius .and. pos_prod(1) + pos_prod(2) .le. 1.0_SR + radius) then
                l_refine_initial = .true.

                !production well:
                !set a constant pressure condition and an outflow saturation condition
                is_dirichlet = .true.
                p = cfg%r_p_prod

                call gv_p%write(element, p)
                call gv_is_dirichlet%add(element, is_dirichlet)
            end if

            !saturation = 1.0_SR
            !call gv_saturation%write(element, saturation)

            pos_in = samoa_world_to_barycentric_point(element%transform_data, cfg%r_pos_in)
            if (base_permeability > 0.0_SR .and. pos_in(1) .ge. -radius .and. pos_in(2) .ge. -radius .and. pos_in(1) + pos_in(2) .le. 1.0_SR + radius) then
                l_refine_initial = .true.

                !injection well:
                !set an inflow pressure condition and a constant saturation condition
                saturation = 1.0_SR
                call gv_saturation%write(element, saturation)

                rhs(:) = 0.009201_SR / 51.816_SR * min(1.0_SR, 1.0_SR / (0.0508_SR * 0.0508_SR * pi) * ((cfg%scaling * element%transform_data%custom_data%scaling) ** 2)) * [1.0_SR/6.0_SR, 1.0_SR/6.0_SR, 1.0_SR/6.0_SR]
            end if

			!compute permeability

			r_lambda_w = (saturation * saturation) / cfg%r_nu_w
			r_lambda_n = (1.0_SR - saturation) * (1.0_SR - saturation) / cfg%r_nu_n

			permeability = base_permeability * (dot_product([0.25_SR, 0.5_SR, 0.25_SR], r_lambda_w + r_lambda_n))
            g_local = samoa_world_to_barycentric_normal(element%transform_data, g / cfg%scaling)

            rhs = rhs + base_permeability * dot_product(cfg%r_rho_w * r_lambda_w + cfg%r_rho_n * r_lambda_n, volumes) * (g_local(1) * dx + g_local(2) * dy)

			l_refine_solution = max(abs(saturation(3) - saturation(2)), abs(saturation(1) - saturation(2))) > 0.1_SR
			l_refine_solution = l_refine_solution .or. max(abs(p(3) - p(2)), abs(p(1) - p(2))) > 0.01_SR * cfg%r_p_prod

			!refine the cell if necessary (no coarsening in the initialization!)

			i_depth = element%cell%geometry%i_depth

			if (i_depth < cfg%i_max_depth .and. (((base_permeability > 0.0_SR) .and. (i_depth < cfg%i_min_depth .or. l_refine_solution)) .or. l_refine_initial)) then
				element%cell%geometry%refinement = 1
				traversal%i_refinements_issued = traversal%i_refinements_issued + 1_DI
			else
				element%cell%geometry%refinement = 0
			end if
		end subroutine
	END MODULE
#endif
