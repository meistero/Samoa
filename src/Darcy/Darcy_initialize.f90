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

		public get_permeability_and_porosity_at_element, get_permeability_and_porosity_at_position, mean_transform, mean_invert

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

			call get_permeability_and_porosity_at_element(section, element, element%cell%data_pers%base_permeability, element%cell%data_pers%porosity)

            !set initial refinement to 0
            element%cell%geometry%refinement = 0
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

        !> converts the value p to a transformed space where the respective mean is additive
		elemental function mean_transform(p) result(p_trans)
            real (kind = SR), intent(in)    :: p
            real (kind = SR)                :: p_trans

#           if defined(_PERM_MEAN_ARITHMETIC)
                p_trans = p
#           elif defined(_PERM_MEAN_HARMONIC)
                !this will cause an intended floating point overflow if p = 0. Backtranformation will return the value 0.
                p_trans = 1.0_SR / p
#           elif defined(_PERM_MEAN_GEOMETRIC)
                !this will cause an intended floating point overflow if p = 0. Backtranformation will return the value 0.
                p_trans = log(p)
#           endif
		end function

        !> converts the value p_trans back to normal space
        elemental function mean_invert(p_trans) result(p)
            real (kind = SR), intent(in)    :: p_trans
            real (kind = SR)                :: p

#           if defined(_PERM_MEAN_ARITHMETIC)
                p = p_trans
#           elif defined(_PERM_MEAN_HARMONIC)
                p = 1.0_SR / p_trans
#           elif defined(_PERM_MEAN_GEOMETRIC)
                p = exp(p_trans)
#           endif
		end function

        recursive subroutine refine_2D_recursive(section, x1, x2, x3, perm, por, no_samples, depth, nz, x_min, x_max)
 			type(t_grid_section), intent(inout)     :: section
            real (kind = GRID_SR), intent(in)       :: x1(:), x2(:), x3(:)
            integer, intent(in)                     :: depth, nz
            real (kind = GRID_SR), intent(in)       :: x_min(:), x_max(:)

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR), intent(inout)    :: perm(:, :), por(:)
                integer, intent(inout)                  :: no_samples(:)
                real (kind = GRID_SR)                   :: perm_buffer(2), por_buffer
#           else
                real (kind = GRID_SR), intent(inout)    :: perm, por
                integer, intent(inout)                  :: no_samples
                real (kind = GRID_SR)                   :: perm_buffer, por_buffer
#           endif

            real (kind = GRID_SR) :: x(3)
            integer :: i, level, k, k_start, k_end

            if (depth > 0) then
                call refine_2D_recursive(section, x1, 0.5_SR * (x1 + x3), x2, perm, por, no_samples, depth - 1, nz, x_min, x_max)
                call refine_2D_recursive(section, x2, 0.5_SR * (x1 + x3), x3, perm, por, no_samples, depth - 1, nz, x_min, x_max)
            else
                x(1:2) = (x1 + x2 + x3) / 3.0_SR

                if (any(x(1:2) < x_min .or. x(1:2) > x_max)) then
                    return
                endif

#               if (_DARCY_LAYERS > 0)
                    do level = 1, _DARCY_LAYERS
                        !find the k-range in the source data that we have to read from
                        !and ensure that we read at least one value from the source data
                        k_start = ((level - 1) * nz) / _DARCY_LAYERS + 1
                        k_end = max(k_start, (level * nz) / _DARCY_LAYERS)

                        do k = k_start, k_end
                            x(3) = (real(k, SR) - 0.5_SR) / real(nz, SR)
                            call get_permeability_and_porosity_at_position(section, x, perm_buffer, por_buffer)
                            perm(level, :) = perm(level, :) + mean_transform(perm_buffer)
                            por(level) = por(level) + por_buffer
                        end do

                        no_samples(level) = no_samples(level) + (k_end - k_start + 1)
                    end do
#               else
                    k_start = 1
                    k_end = 1

                    do k = k_start, k_end
                        x(3) = (real(k, SR) - 0.5_SR) / real(nz, SR)
                        call get_permeability_and_porosity_at_position(section, x, perm_buffer, por_buffer)
                        perm = perm + mean_transform(perm_buffer)
                        por = por + por_buffer
                    end do

                    no_samples = no_samples + (k_end - k_start + 1)
#               endif
            end if
        end subroutine

        subroutine get_permeability_and_porosity_at_element(section, element, permeability, porosity)
			type(t_grid_section), intent(inout)     :: section
			type(t_element_base), intent(inout)     :: element

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR), intent(out)		:: permeability(:, :)	!< permeability tensor
                real (kind = GRID_SR), intent(out)		:: porosity(:)          !< porosity
#           else
                real (kind = GRID_SR), intent(out)		:: permeability		    !< permeability tensor
                real (kind = GRID_SR), intent(out)		:: porosity	            !< porosity
#           endif

            real (kind = GRID_SR)   :: x(3), x1(2), x2(2), x3(2), x_min(2), x_max(2)
            integer					:: level, i, ddepth, nz

#           if (_DARCY_LAYERS > 0)
                integer :: no_samples(_DARCY_LAYERS)
#           else
                integer :: no_samples
#           endif

#           if defined(_ADAPT_INTEGRATE)
#               if defined(_ASAGI)
                    ddepth = nint(log(1.0_SR / (asagi_grid_delta(cfg%afh_permeability_X, 0) * asagi_grid_delta(cfg%afh_permeability_X, 1) * _M * _M)) / log(2.0_SR)) - element%cell%geometry%i_depth
                    nz = max(1, nint(cfg%dz * real(max(1, _DARCY_LAYERS, SR)) / (asagi_grid_delta(cfg%afh_permeability_X, 2)  * _M)))
#               else
                    ddepth = cfg%i_max_depth - element%cell%geometry%i_depth
                    nz = 1
#               endif

                x1 = samoa_barycentric_to_world_point(element%transform_data, [1.0_SR, 0.0_SR])
                x2 = samoa_barycentric_to_world_point(element%transform_data, [0.0_SR, 0.0_SR])
                x3 = samoa_barycentric_to_world_point(element%transform_data, [0.0_SR, 1.0_SR])

                permeability = 0.0_SR
                porosity = 0.0_SR
                no_samples = 0

#               if defined(_ASAGI)
                    x_min = ([asagi_grid_min(cfg%afh_permeability_X, 0), asagi_grid_min(cfg%afh_permeability_X, 1)] - cfg%offset) / cfg%scaling
                    x_max = ([asagi_grid_max(cfg%afh_permeability_X, 0), asagi_grid_max(cfg%afh_permeability_X, 1)] - cfg%offset) / cfg%scaling

                    x_min = max(x_min, anint(x_min * 2.0_SR ** real(element%cell%geometry%i_depth / 2, SR)) * 0.5_SR ** real(element%cell%geometry%i_depth / 2, SR))
                    x_max = min(x_max, anint(x_max * 2.0_SR ** real(element%cell%geometry%i_depth / 2, SR)) * 0.5_SR ** real(element%cell%geometry%i_depth / 2, SR))
#               else
                    x_min = 0.0_SR
                    x_max = 1.0_SR
#               endif

                call refine_2D_recursive(section, x1, x2, x3, permeability, porosity, no_samples, 8, nz, x_min, x_max)

#               if (_DARCY_LAYERS > 0)
                    where(no_samples > 0)
                        permeability(:, 1) = mean_invert(permeability(:, 1) / no_samples)
                        permeability(:, 2) = mean_invert(permeability(:, 2) / no_samples)
                        porosity = porosity / no_samples
                    elsewhere
                        permeability(:, 1) = 0.0_SR
                        permeability(:, 2) = 0.0_SR
                        porosity = 0.0_SR
                    endwhere
#               else
                    if (no_samples > 0) then
                        permeability = mean_invert(permeability / no_samples)
                        porosity = porosity / no_samples
                    else
                        permeability = 0.0_SR
                        porosity = 0.0_SR
                    end if
#               endif
#           elif defined(_ADAPT_SAMPLE)
                x(1:2) = samoa_barycentric_to_world_point(dest_element%transform_data, [1.0_SR / 3.0_SR, 1.0_SR / 3.0_SR])

#               if (_DARCY_LAYERS > 0)
                    do level = 1, _DARCY_LAYERS
                        x(3) = (level - 0.5_SR) / real(_DARCY_LAYERS, SR)

                        call get_permeability_and_porosity_at_position(section, x, permeability(level, :), porosity(level))
                    end do
#               else
                    x(3) = 1.0e-5_SR

                    call get_permeability_and_porosity_at_position(section, x, permeability, porosity)
#               endif
#           endif
        end subroutine

		subroutine get_permeability_and_porosity_at_position(section, x, permeability, porosity)
 			type(t_grid_section), intent(inout)					:: section
			real (kind = GRID_SR), intent(in)		            :: x(:)				!< position in world coordinates
			real (kind = GRID_SR), intent(out)				    :: porosity	        !< porosity

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR)							:: permeability(2)	!< permeability tensor
#           else
                real (kind = GRID_SR)							:: permeability		!< permeability tensor
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
                        permeability(1) = asagi_grid_get_float(cfg%afh_permeability_X, real(xs, c_double), 0)
                        !buffer(2) = asagi_grid_get_float(cfg%afh_permeability_Y, real(xs, c_double), 0)
                        permeability(2) = asagi_grid_get_float(cfg%afh_permeability_Z, real(xs, c_double), 0)

                        !assume horizontally isotropic permeability
                        !assert(abs(buffer(1) - buffer(2)) < epsilon(1.0_SR))
#                   else
                        permeability = asagi_grid_get_float(cfg%afh_permeability_X, real(xs, c_double), 0)
                        !buffer(2) = asagi_grid_get_float(cfg%afh_permeability_Y, real(xs, c_double), 0)
                        !buffer(3) = asagi_grid_get_float(cfg%afh_permeability_Z, real(xs, c_double), 0)
#                   endif

                    !permeability is given in millidarcy
                    permeability = permeability * _MDY
                else
                    permeability = 0.0_SR
                end if

                if (asagi_grid_min(cfg%afh_porosity, 0) <= xs(1) .and. asagi_grid_min(cfg%afh_porosity, 1) <= xs(2) &
                        .and. xs(1) <= asagi_grid_max(cfg%afh_porosity, 0) .and. xs(2) <= asagi_grid_max(cfg%afh_porosity, 1)) then

                    porosity = asagi_grid_get_float(cfg%afh_porosity, real(xs, c_double), 0)
                else
                    porosity = 0.0_SR
                end if
#           else
                !set the permeability in m^2
                permeability = 5.0e-12_SR * (_M * _M)
                porosity = 0.2_SR
#           endif

            !reduce porosity to account for residual saturations of wetting and non-wetting phase
            porosity = porosity * (1.0_SR - cfg%S_wr - cfg%S_nr)

#           if defined(_ASAGI_TIMING)
                section%stats%r_asagi_time = section%stats%r_asagi_time + get_wtime()
#           endif
		end subroutine

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

#		define _GT_INNER_NODE_FIRST_TOUCH_OP	inner_node_first_touch_op
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
 			type(t_darcy_init_saturation_traversal), intent(in) :: traversal
 			type(t_grid_section), intent(in)					:: section
			type(t_node_data), intent(inout)			        :: node

            call flow_pre_dof_op(node%position(1), node%position(2), node%data_pers%p, node%data_pers%saturation, node%data_pers%rhs, node%data_temp%volume, node%data_temp%is_pressure_dirichlet_boundary, node%data_temp%is_saturation_dirichlet_boundary)
		end subroutine

		elemental subroutine inner_node_first_touch_op(traversal, section, node)
 			type(t_darcy_init_saturation_traversal), intent(in) :: traversal
 			type(t_grid_section), intent(in)					:: section
			type(t_node_data), intent(inout)			        :: node

            !ensure the temporary variables are set to .false., the may be invalid otherwise
            node%data_temp%is_saturation_dirichlet_boundary = .false.
            node%data_temp%is_pressure_dirichlet_boundary = .false.

			call flow_pre_dof_op(node%position(1), node%position(2), node%data_pers%p, node%data_pers%saturation, node%data_pers%rhs, node%data_temp%volume, node%data_temp%is_pressure_dirichlet_boundary, node%data_temp%is_saturation_dirichlet_boundary)
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

		elemental subroutine flow_pre_dof_op(pos_x, pos_y, p, saturation, rhs, inflow, is_pressure_dirichlet, is_saturation_dirichlet)
			real (kind = GRID_SR), intent(in)					:: pos_x, pos_y
			real (kind = GRID_SR), intent(inout)				:: p
			real (kind = GRID_SR), intent(inout)				:: saturation
			real (kind = GRID_SR), intent(out)					:: rhs
			real (kind = GRID_SR), intent(out)					:: inflow
            logical, intent(inout)                              :: is_pressure_dirichlet, is_saturation_dirichlet

            saturation = 0.0_SR
			rhs = 0.0_SR
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
                else if (pos_x == 1.0_SR) then
                    is_pressure_dirichlet = .true.
                end if
#           endif

            if (is_pressure_dirichlet) then
                p = cfg%r_p_prod
            endif
		end subroutine

		elemental subroutine flow_post_dof_op(saturation, rhs, is_pressure_dirichlet, inflow, total_inflow)
			real (kind = GRID_SR), intent(inout)				:: saturation
			real (kind = GRID_SR), intent(inout)				:: rhs
			real (kind = GRID_SR), intent(in)				    :: inflow
			real (kind = GRID_SR), intent(in)				    :: total_inflow
			logical, intent(in)			                        :: is_pressure_dirichlet

            rhs = rhs + cfg%r_inflow * inflow / total_inflow

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

			!check refinement indicator
            !make sure no coarsening happens during the initial phase
			call compute_refinement_indicator(element, traversal%i_refinements_issued, element%cell%geometry%i_depth, element%cell%geometry%refinement, saturation, p, element%cell%data_pers%base_permeability)
            element%cell%geometry%refinement = max(element%cell%geometry%refinement, 0)

            if (element%cell%geometry%i_depth < cfg%i_max_depth) then
                pos_in = samoa_world_to_barycentric_point(element%transform_data, cfg%r_pos_in)
                if (norm2(pos_in) < 1.0_SR + epsilon(1.0_SR)) then
                    if (pos_in(1) > -epsilon(1.0_SR) .and. 1.0_SR - (pos_in(1) + pos_in(2)) > -epsilon(1.0_SR) .and. pos_in(2) > - epsilon(1.0_SR)) then
                        element%cell%geometry%refinement = 1
                        traversal%i_refinements_issued = traversal%i_refinements_issued + 1
                    end if
                end if

                pos_prod = sign(cfg%r_pos_prod - 0.5_SR, element%transform_data%custom_data%offset - 0.5_SR) + 0.5_SR
                pos_prod = samoa_world_to_barycentric_point(element%transform_data, pos_prod)
                if (norm2(pos_prod) < 4.0_SR + epsilon(1.0_SR)) then
                    if (pos_prod(1) > -3.0_SR .and. 1.0_SR - (pos_prod(1) + pos_prod(2)) > -3.0_SR .and. pos_prod(2) > -3.0_SR) then
                        element%cell%geometry%refinement = 1
                        traversal%i_refinements_issued = traversal%i_refinements_issued + 1
                    end if
                end if
            end if
		end subroutine
	END MODULE
#endif
