! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
	MODULE SWE_Initialize_Bathymetry
        use iso_c_binding
		use SFC_edge_traversal
		use SWE_euler_timestep

		use Samoa_swe

		implicit none

        type num_traversal_data
        end type

		type(t_gv_Q)							:: gv_Q
		type(t_lfs_flux)						:: lfs_flux

		PUBLIC get_bathymetry_at_element, get_bathymetry_at_position

#		define _GT_NAME							t_swe_init_b_traversal

#		define _GT_EDGES
#		define _GT_EDGES_TEMP

#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op
#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define _GT_ELEMENT_OP					element_op

#		define _GT_CELL_TO_EDGE_OP				cell_to_edge_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine pre_traversal_grid_op(traversal, grid)
			type(t_swe_init_b_traversal), intent(inout)		        :: traversal
			type(t_grid), intent(inout)							    :: grid

			grid%r_time = 0.0_GRID_SR
			grid%r_dt = 0.0_GRID_SR
			grid%r_dt_new = 0.0_GRID_SR

            call scatter(grid%r_time, grid%sections%elements_alloc%r_time)
		end subroutine

		subroutine post_traversal_grid_op(traversal, grid)
			type(t_swe_init_b_traversal), intent(inout)		        :: traversal
			type(t_grid), intent(inout)							    :: grid

            call reduce(grid%r_dt_new, grid%sections%elements_alloc%r_dt_new, MPI_MIN, .true.)

            grid%r_dt = grid%r_dt_new
		end subroutine

		subroutine pre_traversal_op(traversal, section)
			type(t_swe_init_b_traversal), intent(inout)				    :: traversal
			type(t_grid_section), intent(inout)							:: section

			!this variable will be incremented for each cell with a refinement request
			section%r_dt_new = huge(1.0_GRID_SR)
		end subroutine

		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, section, element)
			type(t_swe_init_b_traversal), intent(inout)				    :: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)					:: element

			type(t_state), dimension(_SWE_CELL_SIZE)			:: Q

			call alpha_volume_op(traversal, section, element, Q)

			call gv_Q%write(element, Q)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine alpha_volume_op(traversal, section, element, Q)
			type(t_swe_init_b_traversal), intent(inout)				    :: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)						:: element
			type(t_state), dimension(_SWE_CELL_SIZE), intent(out)	:: Q

			real (kind = GRID_SR), dimension(2)						:: pos
			integer (kind = GRID_SI)								:: i
			real (kind = GRID_SR), parameter		                :: r_test_points(2, 3) = reshape([1.0, 0.0, 0.0, 0.0, 0.0, 1.0], [2, 3])
			real (kind = GRID_SR)                                   :: centroid_square(2), centroid_triangle(2)
			type(t_state), dimension(3)								:: Q_test
			real (kind = GRID_SR), dimension(_SWE_CELL_SIZE)		:: lambda

			!evaluate initial function values at dof positions and compute DoFs

			Q%b = get_bathymetry_at_element(section, element)
		end subroutine

        function get_bathymetry_at_element(section, element) result(bathymetry)
			type(t_grid_section), intent(inout)     :: section
			type(t_element_base), intent(inout)     :: element
            real (kind = GRID_SR)					:: bathymetry

            real (kind = GRID_SR)   :: x(2), p(2), v1(2), v2(2), alpha, beta
			integer					:: i, j, n

#define     _ADAPT_INTEGRATE
!#define     _ADAPT_SAMPLE
#           if defined(_ADAPT_INTEGRATE)
#               if defined(_ASAGI)
                    n = max(1, int(cfg%scaling * element%transform_data%custom_data%scaling / asagi_grid_delta(cfg%afh_bathymetry, 0)))
#               else
                    n = max(1, 512 * int(element%transform_data%custom_data%scaling))
#               endif

                p = samoa_barycentric_to_world_point(element%transform_data, [0.0_SR, 0.0_SR])
                v1 = samoa_barycentric_to_world_vector(element%transform_data, [0.5_SR, 0.5_SR])
                v2 = samoa_barycentric_to_world_vector(element%transform_data, [0.5_SR, -0.5_SR])

                bathymetry = 0.0_SR

                do i = 0, n - 1
                    alpha = (i + 0.5_SR) / real(n, SR)

                    do j = -i, i
                        beta = real(j, SR) / real(n, SR)
                        x = p + alpha * v1 + beta * v2

                        bathymetry = bathymetry + get_bathymetry_at_position(section, x, section%r_time)
                    end do
                end do

                bathymetry = bathymetry / (n * n)
#           elif defined(_ADAPT_SAMPLE)
                x = samoa_barycentric_to_world_point(element%transform_data, [1.0_SR / 3.0_SR, 1.0_SR / 3.0_SR])
                bathymetry = get_bathymetry_at_position(section, x, section%r_time)
#           endif
        end function

		function get_bathymetry_at_position(section, x, t) result(bathymetry)
			type(t_grid_section), intent(inout)					:: section
			real (kind = GRID_SR), dimension(:), intent(in)		:: x						!< position in world coordinates
			real (kind = GRID_SR), intent(in)		            :: t						!< simulation time
            real (kind = GRID_SR)								:: bathymetry				!< bathymetry

            real (kind = GRID_SR)                               :: xs(3)

            xs(1:2) = cfg%scaling * x + cfg%offset

#			if defined(_ASAGI)
#               if defined(_ASAGI_TIMING)
                    section%stats%r_asagi_time = section%stats%r_asagi_time - get_wtime()
#               endif

				if (asagi_grid_min(cfg%afh_bathymetry, 0) <= xs(1) .and. asagi_grid_min(cfg%afh_bathymetry, 1) <= xs(2) &
                        .and. xs(1) <= asagi_grid_max(cfg%afh_bathymetry, 0) .and. xs(2) <= asagi_grid_max(cfg%afh_bathymetry, 1)) then

                    xs(3) = 0.0
                    bathymetry = asagi_grid_get_float(cfg%afh_bathymetry, real(xs, c_double), 0)
                else
                    bathymetry = -5000.0 !we assume that the sea floor is constant here
                end if

                if (asagi_grid_min(cfg%afh_displacement, 0) <= xs(1) .and. asagi_grid_min(cfg%afh_displacement, 1) <= xs(2) &
                        .and. xs(1) <= asagi_grid_max(cfg%afh_displacement, 0) .and. xs(2) <= asagi_grid_max(cfg%afh_displacement, 1) &
                        .and. t > cfg%t_min_eq) then

                    xs(3) = min(t, cfg%t_max_eq)
                    bathymetry = bathymetry + asagi_grid_get_float(cfg%afh_displacement, real(xs, c_double), 0)
                end if

#               if defined(_ASAGI_TIMING)
                    section%stats%r_asagi_time = section%stats%r_asagi_time + get_wtime()
#               endif
#			else
				bathymetry = 0.0_SR
#			endif
		end function
	END MODULE

	MODULE SWE_Initialize_Dofs
 		use Tools_noise

        use iso_c_binding
		use SFC_edge_traversal
		use SWE_euler_timestep
		use SWE_Initialize_Bathymetry
		use Samoa_swe

		implicit none

        type num_traversal_data
            integer (kind = GRID_DI)			:: i_refinements_issued
        end type

		type(t_gv_Q)							:: gv_Q
		type(t_lfs_flux)						:: lfs_flux

#		define _GT_NAME							t_swe_init_dofs_traversal

#		define _GT_EDGES
#		define _GT_EDGES_TEMP

#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op
#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define _GT_ELEMENT_OP					element_op

#		define _GT_CELL_TO_EDGE_OP				cell_to_edge_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine pre_traversal_grid_op(traversal, grid)
			type(t_swe_init_dofs_traversal), intent(inout)		        :: traversal
			type(t_grid), intent(inout)							    :: grid

			grid%r_time = 0.0_GRID_SR
			grid%r_dt = 0.0_GRID_SR
			grid%r_dt_new = 0.0_GRID_SR

            call scatter(grid%r_time, grid%sections%elements_alloc%r_time)
		end subroutine

		subroutine post_traversal_grid_op(traversal, grid)
			type(t_swe_init_dofs_traversal), intent(inout)		        :: traversal
			type(t_grid), intent(inout)							    :: grid

            call reduce(traversal%i_refinements_issued, traversal%children%i_refinements_issued, MPI_SUM, .true.)
            call reduce(grid%r_dt_new, grid%sections%elements_alloc%r_dt_new, MPI_MIN, .true.)

            grid%r_dt = grid%r_dt_new
		end subroutine

		subroutine pre_traversal_op(traversal, section)
			type(t_swe_init_dofs_traversal), intent(inout)				    :: traversal
			type(t_grid_section), intent(inout)							:: section

			!this variable will be incremented for each cell with a refinement request
			traversal%i_refinements_issued = 0
			section%r_dt_new = huge(1.0_GRID_SR)
		end subroutine

		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, section, element)
			type(t_swe_init_dofs_traversal), intent(inout)				    :: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)					        :: element

			type(t_state)			                                    :: Q(_SWE_CELL_SIZE)

			call gv_Q%read(element, Q)

			call alpha_volume_op(traversal, section, element, Q)

			call gv_Q%write(element, Q)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine alpha_volume_op(traversal, section, element, Q)
			type(t_swe_init_dofs_traversal), intent(inout)				    :: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)						    :: element
			type(t_state), dimension(_SWE_CELL_SIZE), intent(inout)	    :: Q

			real (kind = GRID_SR), dimension(2)						:: pos
			integer (kind = GRID_SI)								:: i
			real (kind = GRID_SR), parameter		                :: r_test_points(2, 3) = reshape([1.0, 0.0, 0.0, 0.0, 0.0, 1.0], [2, 3])
			real (kind = GRID_SR)                                   :: centroid_square(2), centroid_triangle(2)
			type(t_state), dimension(3)								:: Q_test
			real (kind = GRID_SR), dimension(_SWE_CELL_SIZE)		:: lambda

			!evaluate initial DoFs (not the bathymetry!)

			Q%t_dof_state = get_initial_dof_state_at_element(section, element)

			element%cell%geometry%refinement = 0

			if (element%cell%geometry%i_depth < cfg%i_min_depth) then
                !refine if the minimum depth is not met

 				element%cell%geometry%refinement = 1
				traversal%i_refinements_issued = traversal%i_refinements_issued + 1
            else if (element%cell%geometry%i_depth < cfg%i_max_depth) then
                do i = 1, 3
                    Q_test(i)%t_dof_state = get_initial_dof_state_at_position(section, samoa_barycentric_to_world_point(element%transform_data, r_test_points(:, i)))
                    Q_test(i)%b = get_bathymetry_at_position(section, samoa_barycentric_to_world_point(element%transform_data, r_test_points(:, i)), 0.0_SR)
                end do

#               if defined (_ASAGI)
                    centroid_square = 0.5_GRID_SR * [asagi_grid_min(cfg%afh_displacement, 0) + asagi_grid_max(cfg%afh_displacement, 0), asagi_grid_min(cfg%afh_displacement, 1) + asagi_grid_max(cfg%afh_displacement, 1)]
                    centroid_square = 1.0_GRID_SR / cfg%scaling * (centroid_square - cfg%offset)
                    centroid_square = samoa_world_to_barycentric_point(element%transform_data, centroid_square)

                    centroid_triangle = [1.0_GRID_SR/3.0_GRID_SR, 1.0_GRID_SR/3.0_GRID_SR]
                    centroid_triangle = cfg%scaling * samoa_barycentric_to_world_point(element%transform_data, centroid_triangle) + cfg%offset
                    centroid_triangle = [ &
                        (centroid_triangle(1) - asagi_grid_min(cfg%afh_displacement, 0)) / (asagi_grid_max(cfg%afh_displacement, 0) - asagi_grid_min(cfg%afh_displacement, 0)), &
                        (centroid_triangle(2) - asagi_grid_min(cfg%afh_displacement, 1)) / (asagi_grid_max(cfg%afh_displacement, 1) - asagi_grid_min(cfg%afh_displacement, 1)) &
                    ]

                    if (maxval(Q_test%h - Q_test%b) > 0.0 .and. minval(Q_test%h - Q_test%b) <= 0.0) then
                        !refine coast lines

                        element%cell%geometry%refinement = 1
                        traversal%i_refinements_issued = traversal%i_refinements_issued + 1
                    elseif (centroid_square(1) >= 0.0 .and. centroid_square(2) >= 0.0 .and. centroid_square(1) + centroid_square(2) <= 1.0) then
                        !refine the triangle if it contains the centroid of the initial condition

                        element%cell%geometry%refinement = 1
                        traversal%i_refinements_issued = traversal%i_refinements_issued + 1
                    elseif (centroid_triangle(1) >= 0.0 .and. centroid_triangle(2) >= 0.0 .and. centroid_triangle(1) <= 1.0 .and. centroid_triangle(2) <= 1.0) then
                        !refine the triangle if its centroid is contained in the initial condition

                        element%cell%geometry%refinement = 1
                        traversal%i_refinements_issued = traversal%i_refinements_issued + 1
                    end if
#               else
                    if (maxval(Q_test%h - Q_test%b) > 0.0 .and. minval(Q_test%h - Q_test%b) <= 0.0 .or. any(Q_test%h .ne. 0.0)) then
                        !refine coast lines and initial displacement

                        element%cell%geometry%refinement = 1
                        traversal%i_refinements_issued = traversal%i_refinements_issued + 1
                    end if
#               endif
			end if

			!estimate initial time step

			where (Q%h - Q%b > 0.0_GRID_SR)
				lambda = sqrt(g * (Q%h - Q%b)) + sqrt((Q%p(1) * Q%p(1) + Q%p(2) * Q%p(2)) / ((Q%h - Q%b) * (Q%h - Q%b)))
			elsewhere
				lambda = 0.0_GRID_SR
			end where

			section%r_dt_new = min(section%r_dt_new, cfg%scaling * element%cell%geometry%get_volume() / (sum(element%cell%geometry%get_edge_sizes()) * maxval(lambda)))
		end subroutine

		function get_initial_dof_state_at_element(section, element) result(Q)
			type(t_grid_section), intent(inout)		:: section
			type(t_element_base), intent(inout)     :: element
			type(t_dof_state)						:: Q

			real (kind = GRID_SR)		            :: x(2)

            x = samoa_barycentric_to_world_point(element%transform_data, [1.0_SR / 3.0_SR, 1.0_SR / 3.0_SR])
            Q = get_initial_dof_state_at_position(section, x)
		end function

		function get_initial_dof_state_at_position(section, x) result(Q)
			type(t_grid_section), intent(inout)					:: section
			real (kind = GRID_SR), intent(in)		            :: x(:)            !< position in world coordinates
			type(t_dof_state)							        :: Q

            real (kind = GRID_SR), parameter		            :: hL = 1.0_SR, hR = 0.0_SR
            real (kind = GRID_SR)                               :: xs(2)

            xs = cfg%scaling * x + cfg%offset

#			if defined(_ASAGI)
				Q%h = 0.0_GRID_SR
#			else
				if (xs(1) < 0.0_SR) then
                    Q%h = hL
                else
                    Q%h = hR
                end if
#			endif

			Q%p = 0.0_GRID_SR
		end function
	END MODULE
#endif

