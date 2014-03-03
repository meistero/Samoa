#include "Compilation_control.f90"

#if defined(_FLASH)
	MODULE FLASH_Initialize
 		use Tools_noise

		use SFC_edge_traversal
		use FLASH_euler_timestep

		use Samoa_FLASH

		implicit none

        type num_traversal_data
            integer (kind = GRID_SI)			:: i_refinements_issued
        end type

		type(t_gv_Q)							:: gv_Q
		type(t_lfs_flux)						:: lfs_flux

		PUBLIC get_bathymetry

#		define _GT_NAME							t_FLASH_init_traversal

#		define _GT_EDGES
#		define _GT_EDGES_TEMP

#		define _GT_PRE_TRAVERSAL_OP			pre_traversal_op
#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define _GT_ELEMENT_OP					element_op

#		define _GT_CELL_TO_EDGE_OP				cell_to_edge_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine pre_traversal_grid_op(traversal, grid)
			type(t_FLASH_init_traversal), intent(inout)		        :: traversal
			type(t_grid), intent(inout)							    :: grid

			grid%r_time = 0.0_GRID_SR
			grid%d_max = grid%i_max_depth
			grid%u_max = sqrt(g)

            call scatter(grid%r_time, grid%sections%elements_alloc%r_time)
            call scatter(grid%afh_displacement, grid%sections%elements_alloc%afh_displacement)
            call scatter(grid%afh_bathymetry, grid%sections%elements_alloc%afh_bathymetry)
            call scatter(grid%scaling, grid%sections%elements_alloc%scaling)
            call scatter(grid%offset(1), grid%sections%elements_alloc%offset(1))
            call scatter(grid%offset(2), grid%sections%elements_alloc%offset(2))
		end subroutine

		subroutine post_traversal_grid_op(traversal, grid)
			type(t_FLASH_init_traversal), intent(inout)		        :: traversal
			type(t_grid), intent(inout)							    :: grid

            call reduce(traversal%i_refinements_issued, traversal%children%i_refinements_issued, MPI_SUM, .true.)
            call reduce(grid%u_max, grid%sections%elements_alloc%u_max, MPI_MAX, .true.)
		end subroutine

		subroutine pre_traversal_op(traversal, section)
			type(t_FLASH_init_traversal), intent(inout)				    :: traversal
			type(t_grid_section), intent(inout)							:: section

			!this variable will be incremented for each cell with a refinement request
			traversal%i_refinements_issued = 0
			section%u_max = 0.0_GRID_SR
		end subroutine

		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, section, element)
			type(t_FLASH_init_traversal), intent(inout)				    :: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)					:: element

			type(t_state), dimension(_FLASH_CELL_SIZE)			:: Q

			call alpha_volume_op(traversal, section, element, Q)

			call gv_Q%write(element, Q)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine alpha_volume_op(traversal, section, element, Q)
			type(t_FLASH_init_traversal), intent(inout)				    :: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)						:: element
			type(t_state), dimension(_FLASH_CELL_SIZE), intent(out)	:: Q

			real (kind = GRID_SR), dimension(2)						:: pos
			integer (kind = GRID_SI)								:: i
			real (kind = GRID_SR), parameter, dimension(2, 3)		:: r_test_points = reshape([1.0, 0.0, 0.0, 0.0, 0.0, 1.0], [2, 3])
			type(t_state), dimension(3)								:: Q_test
			real (kind = GRID_SR), dimension(_FLASH_CELL_SIZE)		:: lambda

			!evaluate initial function values at dof positions and compute DoFs

			do i = 1, _FLASH_CELL_SIZE
				Q(i) = get_initial_state(section, samoa_barycentric_to_world_point(element%transform_data, t_basis_Q_get_dof_coords(i)), element%cell%geometry%i_depth / 2_GRID_SI)
			end do

			if (_FLASH_ORDER == 0) then
				do i = 1, 3
					Q_test(i) = get_initial_state(section, samoa_barycentric_to_world_point(element%transform_data, r_test_points(:, i)), element%cell%geometry%i_depth / 2_GRID_SI)
				end do
			else
				do i = 1, 3
					Q_test(i)%h = t_basis_Q_eval(r_test_points(:, i), Q%h)
				end do
			endif

			if (element%cell%geometry%i_depth < section%i_min_depth .or. (element%cell%geometry%i_depth < section%i_max_depth .and. ( &
					(maxval(Q_test%h - Q_test%b) > 0.0_GRID_SR .and. minval(Q_test%h - Q_test%b) < 0.0_GRID_SR) .or. &			!refine coast lines
					abs(Q_test(3)%h - Q_test(2)%h) + abs(Q_test(1)%h - Q_test(2)%h) > 0.0_GRID_SR))) then						!refine waves

				element%cell%geometry%refinement = 1
				traversal%i_refinements_issued = traversal%i_refinements_issued + 1
			else
				element%cell%geometry%refinement = 0
			end if

			!estimate initial u_max

			where (Q%h - Q%b > 0)
				lambda = sqrt(g * (Q%h - Q%b)) + sqrt((Q%p(1) * Q%p(1) + Q%p(2) * Q%p(2)) / ((Q%h - Q%b) * (Q%h - Q%b)))
			elsewhere
				lambda = 0.0_GRID_SR
			end where

			section%u_max = max(section%u_max, maxval(lambda))
		end subroutine

		function get_initial_state(section, x, lod) result(Q)
			type(t_grid_section), intent(inout)					:: section
			real (kind = GRID_SR), dimension(:), intent(in)		:: x
			integer (kind = GRID_SI), intent(in)				:: lod						!< level of detail
			type(t_state)										:: Q

			!init height, momentum and bathymetry
			Q%t_dof_state = get_initial_dof_state(section, x, lod)
			Q%b = get_bathymetry(section, x, 0.0_GRID_SR, lod)
			!Q%b = 0.0 !<-------HACK
		end function

		function get_initial_dof_state(section, x, lod) result(Q)
			type(t_grid_section), intent(inout)							:: section
			real (kind = GRID_SR), dimension(:), intent(in)		:: x						!< position in world coordinates
			integer (kind = GRID_SI), intent(in)				:: lod						!< level of detail
			type(t_dof_state)									:: Q						!< initial state

			real (kind = GRID_SR), dimension(2), parameter		:: dam_center = [0.5, 0.5]
			real (kind = GRID_SR), parameter					:: dam_radius = 0.02
			real (kind = GRID_SR), parameter					:: outer_height = 0.0
			real (kind = GRID_SR), parameter					:: inner_height = 10.0

#			if defined(_ASAGI)
				Q%h = 0
				Q%h_old = Q%h
#			else
				Q%h = 0.5_GRID_SR * (inner_height + outer_height) + (inner_height - outer_height) * sign(0.5_GRID_SR, (dam_radius ** 2) - dot_product(x - dam_center, x - dam_center))
				Q%h_old = Q%h
#			endif

			Q%p = 0.0_GRID_SR
			Q%p_old = Q%p
		end function

		function get_bathymetry(section, x, t, lod) result(bathymetry)
			type(t_grid_section), intent(inout)					:: section
			real (kind = GRID_SR), dimension(:), intent(in)		:: x						!< position in world coordinates
			real (kind = GRID_SR), intent(in)		            :: t						!< simulation time
			integer (kind = GRID_SI), intent(in)				:: lod						!< level of detail
			real (kind = GRID_SR)								:: bathymetry				!< bathymetry

#			if defined(_ASAGI)
				real (kind = GRID_SR) :: xs(2), ts

				xs = section%scaling * x + section%offset

#               if defined(_ASAGI_TIMING)
                    section%stats%r_asagi_time = section%stats%r_asagi_time - omp_get_wtime()
#               endif

				assert(grid_min_x(section%afh_bathymetry) <= xs(1) .and. grid_min_y(section%afh_bathymetry) <= xs(2) &
						.and. xs(1) <= grid_max_x(section%afh_bathymetry) .and. xs(2) <= grid_max_y(section%afh_bathymetry))

                bathymetry = asagi_get_float(section%afh_bathymetry, xs(1), xs(2), 0)

                if(grid_min_x(section%afh_displacement) <= xs(1) .and. grid_min_y(section%afh_displacement) <= xs(2) &
                        .and. xs(1) <= grid_max_x(section%afh_displacement) .and. xs(2) <= grid_max_y(section%afh_displacement) &
                        .and. grid_min_z(section%afh_displacement) < t) then

                    ts = min(t, grid_max_z(section%afh_displacement))
                    bathymetry = bathymetry + asagi_get_float(section%afh_displacement, xs(1) * cos(.52) + xs(2)*sin(.52), xs(2)*cos(.52)-xs(1)*sin(.52), ts, 0)
                end if

#               if defined(_ASAGI_TIMING)
                    section%stats%r_asagi_time = section%stats%r_asagi_time + omp_get_wtime()
#               endif
#			else
                real (kind = GRID_SR), dimension(2), parameter		:: dam_center = [0.5, 0.5]
                real (kind = GRID_SR), parameter					:: dam_radius = 0.1
                real (kind = GRID_SR), parameter					:: outer_height = -100.0
                real (kind = GRID_SR), parameter					:: inner_height = -5.0

				bathymetry = 0.5_GRID_SR * (inner_height + outer_height) + (inner_height - outer_height) * sign(0.5_GRID_SR, (dam_radius ** 2) - dot_product(x - dam_center, x - dam_center))
#			endif
		end function
	END MODULE
#endif
