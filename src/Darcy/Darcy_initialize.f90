! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_initialize_pressure
		use SFC_edge_traversal

        use iso_c_binding
		use Samoa_darcy
		use Tools_noise

		implicit none

        type num_traversal_data
        end type

		type(darcy_gv_p)						:: gv_p
		type(darcy_gv_saturation)				:: gv_saturation

		public get_base_permeability, get_porosity

#		define	_GT_NAME						t_darcy_init_pressure_traversal

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

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR), intent(out)									:: base_permeability(:, :)
                real (kind = GRID_SR), intent(out)									:: porosity(:)

                real (kind = GRID_SR)       :: x(3)
                integer (kind = GRID_SI)    :: i

                x(1:2) = samoa_barycentric_to_world_point(element%transform_data, [1.0_SR / 3.0_SR, 1.0_SR / 3.0_SR])

                do i = 1, _DARCY_LAYERS
                    x(3) = (i - 0.5_SR) / real(_DARCY_LAYERS, SR)

                    !horizontal permeability is assumed to be isotropic
                    base_permeability(i, :) = get_base_permeability(section, x, element%cell%geometry%i_depth / 2_GRID_SI)
                    porosity(i) = get_porosity(section, x)
                end do
#           else
                real (kind = GRID_SR), intent(out)									:: base_permeability
                real (kind = GRID_SR), intent(out)									:: porosity

                real (kind = GRID_SR)   :: x(3)

                x(1:2) = samoa_barycentric_to_world_point(element%transform_data, [1.0_SR / 3.0_SR, 1.0_SR / 3.0_SR])
                x(3) = 1.0e-5_SR

                base_permeability = get_base_permeability(section, x, element%cell%geometry%i_depth / 2_GRID_SI)
                porosity = get_porosity(section, x)
#           endif

            !set initial refinement to 0
            element%cell%geometry%refinement = 0
 		end subroutine

		function get_base_permeability(section, x, lod) result(r_base_permeability)
 			type(t_grid_section), intent(inout)					:: section
			real (kind = GRID_SR), intent(in)		            :: x(:)						!< position in world coordinates
			integer (kind = GRID_SI), intent(in)				:: lod						!< level of detail

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR)							:: r_base_permeability(2)	!< permeability tensor
#           else
                real (kind = GRID_SR)							:: r_base_permeability		!< permeability tensor
#           endif

            real (kind = GRID_SR)                               :: xs(3)

            xs(1:2) = cfg%scaling * x(1:2) + cfg%offset
            xs(3) = cfg%scaling * max(1, _DARCY_LAYERS) * cfg%dz * x(3)

            assert_ge(x(1), 0.0); assert_ge(x(2), 0.0)
            assert_le(x(1), 1.0); assert_le(x(2), 1.0)

#           if defined(_ASAGI_TIMING)
                section%stats%r_asagi_time = section%stats%r_asagi_time - get_wtime()
#           endif

#           if defined(_ASAGI)
                if (asagi_grid_min(cfg%afh_permeability_X, 0) <= xs(1) .and. asagi_grid_min(cfg%afh_permeability_X, 1) <= xs(2) &
                        .and. xs(1) <= asagi_grid_max(cfg%afh_permeability_X, 0) .and. xs(2) <= asagi_grid_max(cfg%afh_permeability_X, 1)) then

#                   if (_DARCY_LAYERS > 0)
                        r_base_permeability(1) = asagi_grid_get_float(cfg%afh_permeability_X, real(xs, c_double), 0)
                        !buffer(2) = asagi_grid_get_float(cfg%afh_permeability_Y, real(xs, c_double), 0)
                        r_base_permeability(2) = asagi_grid_get_float(cfg%afh_permeability_Z, real(xs, c_double), 0)

                        !assume horizontally isotropic permeability
                        !assert(abs(buffer(1) - buffer(2)) < epsilon(1.0_SR))
#                   else
                        r_base_permeability = asagi_grid_get_float(cfg%afh_permeability_X, real(xs, c_double), 0)
                        !buffer(2) = asagi_grid_get_float(cfg%afh_permeability_Y, real(xs, c_double), 0)
                        !buffer(3) = asagi_grid_get_float(cfg%afh_permeability_Z, real(xs, c_double), 0)
#                   endif

                    !convert from mD to m^2 to um^2
                    r_base_permeability = r_base_permeability * 9.869233e-16_SR / (cfg%scaling ** 2)
                else
                    r_base_permeability = 0.0_SR
                end if
#           else
                r_base_permeability = 5.0e-12_SR / (cfg%scaling ** 2)
#           endif

#           if defined(_ASAGI_TIMING)
                section%stats%r_asagi_time = section%stats%r_asagi_time + get_wtime()
#           endif
		end function

		function get_porosity(section, x) result(porosity)
 			type(t_grid_section), intent(inout)					:: section
			real (kind = GRID_SR), dimension(:), intent(in)		:: x						!< position in world coordinates
			real (kind = GRID_SR)								:: porosity		        !< porosity

            real (kind = GRID_SR)                               :: xs(3)

            xs(1:2) = cfg%scaling * x(1:2) + cfg%offset
            xs(3) = cfg%scaling * max(1, _DARCY_LAYERS) * cfg%dz * x(3)

            assert_ge(x(1), 0.0); assert_ge(x(2), 0.0)
            assert_le(x(1), 1.0); assert_le(x(2), 1.0)

#           if defined(_ASAGI_TIMING)
                section%stats%r_asagi_time = section%stats%r_asagi_time - get_wtime()
#           endif

#           if defined(_ASAGI)
                if (asagi_grid_min(cfg%afh_porosity, 0) <= xs(1) .and. asagi_grid_min(cfg%afh_porosity, 1) <= xs(2) &
                        .and. xs(1) <= asagi_grid_max(cfg%afh_porosity, 0) .and. xs(2) <= asagi_grid_max(cfg%afh_porosity, 1)) then

                    porosity = asagi_grid_get_float(cfg%afh_porosity, real(xs, c_double), 0)
                else
                    porosity = 0.0_SR
                end if
#           else
                porosity = 0.2_SR
#           endif

            !reduce porosity to account for residual saturations of wetting and non-wetting phase
            porosity = porosity * (1.0_SR - cfg%S_wr - cfg%S_nr)

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
		use Darcy_permeability
		use Darcy_error_estimate

		implicit none

        type num_traversal_data
            integer (kind = GRID_DI)			:: i_refinements_issued
        end type

		type(darcy_gv_saturation)				:: gv_saturation
		type(darcy_gv_p)						:: gv_p
		type(darcy_gv_rhs)						:: gv_rhs
		type(darcy_gv_is_pressure_dirichlet_boundary)    :: gv_is_pressure_dirichlet
		type(darcy_gv_is_saturation_dirichlet_boundary)    :: gv_is_saturation_dirichlet
		type(darcy_gv_volume)					:: gv_inflow

#		define	_GT_NAME						t_darcy_init_saturation_traversal

#		define _GT_NODES

#		define _GT_NODE_FIRST_TOUCH_OP		    node_first_touch_op
#		define _GT_NODE_LAST_TOUCH_OP		    node_last_touch_op

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

		elemental subroutine node_first_touch_op(traversal, section, node)
 			type(t_darcy_init_saturation_traversal), intent(in)                     :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			:: node

			call flow_pre_dof_op(node%position(1), node%position(2), node%data_pers%p, node%data_pers%saturation, node%data_pers%rhs, node%data_temp%is_pressure_dirichlet_boundary, node%data_temp%is_saturation_dirichlet_boundary, node%data_temp%volume)
		end subroutine

		subroutine element_op(traversal, section, element)
 			type(t_darcy_init_saturation_traversal)                 :: traversal
 			type(t_grid_section), intent(inout)						:: section
			type(t_element_base), intent(inout)				        :: element

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR)   :: saturation(_DARCY_LAYERS + 1, 3)
                real (kind = GRID_SR)   :: p(_DARCY_LAYERS + 1, 3)
                real (kind = GRID_SR)   :: rhs(_DARCY_LAYERS + 1, 3)
#           else
                real (kind = GRID_SR)   :: saturation(3)
                real (kind = GRID_SR)   :: p(3)
                real (kind = GRID_SR)   :: rhs(3)
#           endif

			call gv_saturation%read_from_element(element, saturation)
			call gv_p%read_from_element(element, p)

			!call element operator
			call alpha_volume_op(traversal, element, saturation, p, rhs, element%cell%data_pers%base_permeability)

			call gv_rhs%add_to_element(element, rhs)
		end subroutine

		elemental subroutine node_merge_op(local_node, neighbor_node)
			type(t_node_data), intent(inout)		:: local_node
 			type(t_node_data), intent(in)		    :: neighbor_node

            where (neighbor_node%data_temp%is_pressure_dirichlet_boundary)
                local_node%data_pers%p = neighbor_node%data_pers%p
            end where

			local_node%data_pers%saturation = max(local_node%data_pers%saturation, neighbor_node%data_pers%saturation)
			local_node%data_temp%is_pressure_dirichlet_boundary = local_node%data_temp%is_pressure_dirichlet_boundary .or. neighbor_node%data_temp%is_pressure_dirichlet_boundary
			local_node%data_temp%is_saturation_dirichlet_boundary = local_node%data_temp%is_saturation_dirichlet_boundary .or. neighbor_node%data_temp%is_saturation_dirichlet_boundary
			local_node%data_pers%rhs = local_node%data_pers%rhs + neighbor_node%data_pers%rhs
			local_node%data_temp%volume = local_node%data_temp%volume + neighbor_node%data_temp%volume
		end subroutine

		elemental subroutine node_last_touch_op(traversal, section, node)
 			type(t_darcy_init_saturation_traversal), intent(in)                     :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			:: node

			real (kind = GRID_SR) :: total_inflow

            total_inflow = tiny(1.0_SR) + sum(node%data_temp%volume)

			call flow_post_dof_op(node%data_pers%saturation, node%data_pers%rhs, node%data_temp%is_pressure_dirichlet_boundary, node%data_temp%volume, total_inflow)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		elemental subroutine flow_pre_dof_op(pos_x, pos_y, p, saturation, rhs, is_pressure_dirichlet, is_saturation_dirichlet, inflow)
			real (kind = GRID_SR), intent(in)					:: pos_x, pos_y
			real (kind = GRID_SR), intent(inout)				:: p
			real (kind = GRID_SR), intent(inout)				:: saturation
			real (kind = GRID_SR), intent(out)					:: rhs
			logical, intent(out)			                    :: is_pressure_dirichlet, is_saturation_dirichlet
			real (kind = GRID_SR), intent(out)					:: inflow

            saturation = 0.0_SR
			rhs = 0.0_SR
			is_pressure_dirichlet = .false.
			is_saturation_dirichlet = .false.
			inflow = 0.0_SR

#           if !defined(_ASAGI)
                if (pos_x < 0.5_SR) then
                    saturation = 1.0_SR
                else if (pos_x > 0.5_SR) then
                    saturation = 0.0_SR
                else
                    saturation = 0.5_SR
                end if

                if (pos_x == 0.0_SR) then
                    is_saturation_dirichlet = .true.
                    saturation = 1.0_SR
                else if (pos_x == 1.0_SR) then
                    is_pressure_dirichlet = .true.
                    p = cfg%r_p_prod
                end if
#           endif
		end subroutine

		elemental subroutine flow_post_dof_op(saturation, rhs, is_pressure_dirichlet, inflow, total_inflow)
			real (kind = GRID_SR), intent(inout)				:: saturation
			real (kind = GRID_SR), intent(inout)				:: rhs
			real (kind = GRID_SR), intent(in)				    :: inflow
			real (kind = GRID_SR), intent(in)				    :: total_inflow
			logical, intent(in)			                        :: is_pressure_dirichlet

            rhs = rhs + cfg%r_inflow * inflow / total_inflow

            !limit the accumulated saturation at the injection well
            saturation = min(saturation, 1.0_SR)

            if (is_pressure_dirichlet) then
                rhs = 0.0_SR
            end if
		end subroutine

		subroutine alpha_volume_op(traversal, element, saturation, p, rhs, base_permeability)
 			type(t_darcy_init_saturation_traversal)                             :: traversal
			type(t_element_base), intent(inout)									:: element

			real (kind = GRID_SR)					                :: pos_prod(2), pos_in(2)

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR), intent(inout)		        :: saturation(:, :)
                real (kind = GRID_SR), intent(inout)			    :: p(:, :)
                real (kind = GRID_SR), intent(inout)				:: rhs(:, :)
                real (kind = GRID_SR), intent(inout)                :: base_permeability(:, :)
#           else
                real (kind = GRID_SR), intent(inout)		        :: saturation(:)
                real (kind = GRID_SR), intent(inout)			    :: p(:)
                real (kind = GRID_SR), intent(inout)				:: rhs(:)
                real (kind = GRID_SR), intent(in)                   :: base_permeability
#           endif

            !set boundary conditions and source terms

            call initialize_rhs(element, saturation, p, rhs, base_permeability)

            !saturation = 1.0_SR
            !call gv_saturation%write_to_element(element, saturation)

			!check refinement indicator
            !make sure no coarsening happens during the initial phase
			call compute_refinement_indicator(element, traversal%i_refinements_issued, element%cell%geometry%i_depth, element%cell%geometry%refinement, saturation, p, element%cell%data_pers%base_permeability)
            element%cell%geometry%refinement = max(element%cell%geometry%refinement, 0)

            if (element%cell%geometry%i_depth < cfg%i_max_depth) then
                pos_in = samoa_world_to_barycentric_point(element%transform_data, cfg%r_pos_in)
                if (norm2(pos_in) < 1.0_SR + epsilon(1.0_SR)) then
                    if (pos_in(1) > -epsilon(1.0_SR) .and. 1.0_SR - (pos_in(1) + pos_in(2)) > -epsilon(1.0_SR) .and. pos_in(2) > -epsilon(1.0_SR)) then
                        element%cell%geometry%refinement = 1
                        traversal%i_refinements_issued = traversal%i_refinements_issued + 1
                    end if
                end if

                pos_prod = sign(cfg%r_pos_prod - 0.5_SR, element%transform_data%custom_data%offset - 0.5_SR) + 0.5_SR
                pos_prod = samoa_world_to_barycentric_point(element%transform_data, pos_prod)
                if (norm2(pos_prod) < 1.0_SR + epsilon(1.0_SR)) then
                    if (pos_prod(1) > -epsilon(1.0_SR) .and. 1.0_SR - (pos_prod(1) + pos_prod(2)) > -epsilon(1.0_SR) .and. pos_prod(2) > -epsilon(1.0_SR)) then
                        element%cell%geometry%refinement = 1
                        traversal%i_refinements_issued = traversal%i_refinements_issued + 1
                    end if
                end if
            end if
		end subroutine
	END MODULE
#endif
