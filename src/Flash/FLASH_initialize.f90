! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_FLASH)
	MODULE FLASH_Initialize_Bathymetry
        use iso_c_binding
		use SFC_edge_traversal
		use FLASH_euler_timestep

		use Samoa_FLASH
#       if defined(_FLASH_PATCH)
            use FLASH_PATCH
#       endif

        ! No ASAGI -> Artificial scenario selector 
#       if !defined(_ASAGI)
            use FLASH_Scenario
#       endif 

		implicit none

        type num_traversal_data
        end type

		type(t_gv_Q)							:: gv_Q
		type(t_lfs_flux)						:: lfs_flux

		PUBLIC get_bathymetry_at_element, get_bathymetry_at_position
#       if defined (_FLASH_PATCH)
            PUBLIC get_bathymetry_at_patch
#       endif

#		define _GT_NAME							t_FLASH_init_b_traversal

#		define _GT_EDGES
#		define _GT_EDGES_TEMP

#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op
#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define _GT_ELEMENT_OP					element_op

#		define _GT_CELL_TO_EDGE_OP				cell_to_edge_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine pre_traversal_grid_op(traversal, grid)
			type(t_FLASH_init_b_traversal), intent(inout)		        :: traversal
			type(t_grid), intent(inout)							    :: grid

			grid%r_time = 0.0_GRID_SR
			grid%r_dt = 0.0_GRID_SR

            call scatter(grid%r_time, grid%sections%elements_alloc%r_time)
		end subroutine

		subroutine post_traversal_grid_op(traversal, grid)
			type(t_FLASH_init_b_traversal), intent(inout)		        :: traversal
			type(t_grid), intent(inout)							    :: grid

		end subroutine

		subroutine pre_traversal_op(traversal, section)
			type(t_FLASH_init_b_traversal), intent(inout)				    :: traversal
			type(t_grid_section), intent(inout)							:: section

			!this variable will be incremented for each cell with a refinement request
			!section%r_dt_new = huge(1.0_GRID_SR)
		end subroutine

		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, section, element)
			type(t_FLASH_init_b_traversal), intent(inout)				    :: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)					:: element
#           if defined(_FLASH_PATCH)
                type(t_state), dimension(_FLASH_PATCH_ORDER_SQUARE)   :: Q
#           else
			type(t_state), dimension(_FLASH_CELL_SIZE)			:: Q
#           endif

			call alpha_volume_op(traversal, section, element, Q)

#           if defined(_FLASH_PATCH)
                element%cell%data_pers%B = Q(:)%b
#           else
			call gv_Q%write(element, Q)
#           endif
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine alpha_volume_op(traversal, section, element, Q)
			type(t_FLASH_init_b_traversal), intent(inout)				    :: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)						:: element
#           if defined(_FLASH_PATCH)
                type(t_state), dimension(_FLASH_PATCH_ORDER_SQUARE), intent(out)  :: Q
#           else
                type(t_state), dimension(_FLASH_CELL_SIZE), intent(out)   :: Q
#           endif

			real (kind = GRID_SR), dimension(2)						:: pos
			integer (kind = GRID_SI)								:: i
			real (kind = GRID_SR), parameter		                :: r_test_points(2, 3) = reshape([1.0, 0.0, 0.0, 0.0, 0.0, 1.0], [2, 3])
			real (kind = GRID_SR)                                   :: centroid_square(2), centroid_triangle(2)
			type(t_state), dimension(3)								:: Q_test

			!evaluate initial function values at dof positions and compute DoFs

#           if defined (_FLASH_PATCH)
                Q%b = get_bathymetry_at_patch(section, element, section%r_time)
#           else
                Q%b = get_bathymetry_at_element(section, element, section%r_time)
#           endif
		end subroutine

        recursive function refine_2D_recursive(section, x1, x2, x3, t, depth) result(bath)
 			type(t_grid_section), intent(inout)     :: section
            real (kind = GRID_SR), intent(in)       :: x1(:), x2(:), x3(:), t
            real (kind = GRID_SR)                   :: bath
            integer, intent(in)                     :: depth

            if (depth > 0) then
                bath = 0.5_GRID_SR * (  refine_2D_recursive(section, x1, 0.5_GRID_SR * (x1 + x3), x2, t, depth - 1) &
                                 + refine_2D_recursive(section, x2, 0.5_GRID_SR * (x1 + x3), x3, t, depth - 1))
            else
                bath = get_bathymetry_at_position(section, (x1 + x2 + x3) / 3.0_GRID_SR, t)
            end if
        end function

        function get_bathymetry_at_element(section, element, t) result(bathymetry)
			type(t_grid_section), intent(inout)     :: section
			type(t_element_base), intent(inout)     :: element
            real (kind = GRID_SR), intent(in)		:: t
            real (kind = GRID_SR)					:: bathymetry

            real (kind = GRID_SR)   :: x(2), x1(2), x2(2), x3(2)
            integer                 :: ddepth, data_depth

#           if defined(_ADAPT_INTEGRATE)
                !limit to 16 refinement levels and maximum depth + 3
#               if defined(_ASAGI)
                    data_depth = nint(log(cfg%scaling ** 2 / (asagi_grid_delta(cfg%afh_bathymetry, 0) * asagi_grid_delta(cfg%afh_bathymetry, 1))) / log(2.0_GRID_SR))
                    ddepth = min(16, min(cfg%i_max_depth + 3, data_depth + 3) - element%cell%geometry%i_depth)
#               else
                    ddepth = min(16, cfg%i_max_depth - element%cell%geometry%i_depth)
#               endif

                x1 = samoa_barycentric_to_world_point(element%transform_data, [1.0_GRID_SR, 0.0_GRID_SR])
                x2 = samoa_barycentric_to_world_point(element%transform_data, [0.0_GRID_SR, 0.0_GRID_SR])
                x3 = samoa_barycentric_to_world_point(element%transform_data, [0.0_GRID_SR, 1.0_GRID_SR])
                
                bathymetry = refine_2D_recursive(section, x1, x2, x3, t, ddepth)
#           elif defined(_ADAPT_SAMPLE)
                x = samoa_barycentric_to_world_point(element%transform_data, [1.0_GRID_SR / 3.0_GRID_SR, 1.0_GRID_SR / 3.0_GRID_SR])
                bathymetry = get_bathymetry_at_position(section, x, t)
#           endif
        end function
        
#if defined (_FLASH_PATCH)
        function get_bathymetry_at_patch(section, element, t) result(bathymetry)
            type(t_grid_section), intent(inout)     :: section
            type(t_element_base), intent(inout)     :: element
            real (kind = GRID_SR), intent(in)       :: t
            real (kind = GRID_SR), dimension(_FLASH_PATCH_ORDER_SQUARE) :: bathymetry
            
            integer (kind = GRID_SI)                                        :: i, j, row, col, cell_id

            real (kind = GRID_SR)   :: x(2), x1(2), x2(2), x3(2)
            integer                 :: ddepth, data_depth
            
            !iterate through cells in patch
            row = 1
            col = 1
            do i=1, _FLASH_PATCH_ORDER_SQUARE
                ! if orientation is backwards, the plotter uses a transformation that mirrors the cell...
                ! this simple change solves the problem :)
                if (element%cell%geometry%i_plotter_type > 0) then 
                    cell_id = i
                else
                    cell_id = (row-1)*(row-1) + 2 * row - col
                end if
                
#              if defined(_ADAPT_INTEGRATE)
                    !limit to 16 refinement levels and maximum depth + 3
#                  if defined(_ASAGI)
                        data_depth = nint(log(cfg%scaling ** 2 / (asagi_grid_delta(cfg%afh_bathymetry, 0) * asagi_grid_delta(cfg%afh_bathymetry, 1))) / log(2.0_GRID_SR))
                        ddepth = min(16, min(cfg%i_max_depth + 3, data_depth + 3) - element%cell%geometry%i_depth)
#                  else
                        ddepth = min(16, cfg%i_max_depth - element%cell%geometry%i_depth)
#                  endif

                    if (mod(col,2) == 1) then 
                        x1 = samoa_barycentric_to_world_point(element%transform_data, FLASH_PATCH_geometry%coords(:,1,cell_id) )
                        x2 = samoa_barycentric_to_world_point(element%transform_data, FLASH_PATCH_geometry%coords(:,2,cell_id) )
                        x3 = samoa_barycentric_to_world_point(element%transform_data, FLASH_PATCH_geometry%coords(:,3,cell_id) )
                    else
                        x3 = samoa_barycentric_to_world_point(element%transform_data, FLASH_PATCH_geometry%coords(:,1,cell_id) )
                        x2 = samoa_barycentric_to_world_point(element%transform_data, FLASH_PATCH_geometry%coords(:,2,cell_id) )
                        x1 = samoa_barycentric_to_world_point(element%transform_data, FLASH_PATCH_geometry%coords(:,3,cell_id) )
                    end if
                    
                    bathymetry(i) = refine_2D_recursive(section, x1, x2, x3, t, ddepth)
#              elif defined(_ADAPT_SAMPLE)
                    ! x = average of the 3 vertices
                    x = 0
                    do j=1,3
                        x = x  + FLASH_PATCH_geometry%coords(:,j,cell_id) 
                    end do
                    x = x / 3
                    ! transform x to world coordinates
                    x = samoa_barycentric_to_world_point(element%transform_data, x)
                    bathymetry(i) = get_bathymetry_at_position(section, x, t)
#              endif

                col = col + 1
                if (col == 2*row) then
                    col = 1
                    row = row + 1
                end if
            end do
        end function
#endif

		function get_bathymetry_at_position(section, x, t) result(bathymetry)
			type(t_grid_section), intent(inout)					:: section
			real (kind = GRID_SR), dimension(:), intent(in)		:: x						!< position in world coordinates
			real (kind = GRID_SR), intent(in)		            :: t						!< simulation time
            real (kind = GRID_SR)								:: bathymetry				!< bathymetry

            real (kind = c_double)                              :: xs(3)

            xs(1:2) = real(cfg%scaling * x + cfg%offset, c_double)

#			if defined(_ASAGI)
#               if defined(_ASAGI_TIMING)
                    call section%stats%start_time(asagi_time)
#               endif

			if (grid_min_x(cfg%afh_bathymetry) <= xs(1) .and. grid_min_y(cfg%afh_bathymetry) <= xs(2) &
                        .and. xs(1) <= grid_max_x(cfg%afh_bathymetry) .and. xs(2) <= grid_max_y(cfg%afh_bathymetry)) then

                    xs(3) = 0.0
                    bathymetry = asagi_get_float(cfg%afh_bathymetry, xs, 0)
                else
                    bathymetry = -5000.0_GRID_SR !we assume that the sea floor is constant here
                end if

                if (asagi_grid_min(cfg%afh_displacement, 0) <= xs(1) .and. asagi_grid_min(cfg%afh_displacement, 1) <= xs(2) &
                        .and. xs(1) <= asagi_grid_max(cfg%afh_displacement, 0) .and. xs(2) <= asagi_grid_max(cfg%afh_displacement, 1) &
                        .and. t > cfg%t_min_eq) then

                    xs(3) = real(min(t, cfg%t_max_eq), c_double)
                    bathymetry = bathymetry + asagi_grid_get_float(cfg%afh_displacement, xs, 0)
                end if

#               if defined(_ASAGI_TIMING)
                    call section%stats%stop_time(asagi_time)
#               endif
#			else
				bathymetry = FLASH_Scenario_get_bathymetry(real(xs, GRID_SR))
#			endif
		end function
	END MODULE

	MODULE FLASH_Initialize_Dofs
 		use Tools_noise

        use iso_c_binding
		use SFC_edge_traversal
		use FLASH_euler_timestep
		use FLASH_Initialize_Bathymetry
		use Samoa_FLASH
#       if defined(_FLASH_PATCH)
            use FLASH_PATCH
#       endif

        ! No ASAGI -> Artificial scenario selector 
#       if !defined(_ASAGI)
            use FLASH_Scenario
#       endif 

		implicit none

        type num_traversal_data
            integer (kind = GRID_DI)			:: i_refinements_issued
        end type

		type(t_gv_Q)							:: gv_Q
		type(t_lfs_flux)						:: lfs_flux

#		define _GT_NAME							t_FLASH_init_dofs_traversal

#		define _GT_EDGES
#		define _GT_EDGES_TEMP

#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op
#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define _GT_ELEMENT_OP					element_op

#		define _GT_CELL_TO_EDGE_OP				cell_to_edge_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine pre_traversal_grid_op(traversal, grid)
			type(t_FLASH_init_dofs_traversal), intent(inout)		        :: traversal
			type(t_grid), intent(inout)							    :: grid

			grid%r_time = 0.0_GRID_SR
			grid%r_dt = 0.0_GRID_SR
			grid%r_dt_new = 0.0_GRID_SR

            call scatter(grid%r_time, grid%sections%elements_alloc%r_time)
		end subroutine

		subroutine post_traversal_grid_op(traversal, grid)
			type(t_FLASH_init_dofs_traversal), intent(inout)		        :: traversal
			type(t_grid), intent(inout)							    :: grid

            call reduce(traversal%i_refinements_issued, traversal%sections%i_refinements_issued, MPI_SUM, .true.)
            call reduce(grid%r_dt_new, grid%sections%elements_alloc%r_dt_new, MPI_MIN, .true.)

            grid%r_dt = cfg%courant_number * grid%r_dt_new
		end subroutine

		subroutine pre_traversal_op(traversal, section)
			type(t_FLASH_init_dofs_traversal), intent(inout)				    :: traversal
			type(t_grid_section), intent(inout)							:: section

			!this variable will be incremented for each cell with a refinement request
			traversal%i_refinements_issued = 0
			section%r_dt_new = huge(1.0_GRID_SR)
		end subroutine

		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, section, element)
			type(t_FLASH_init_dofs_traversal), intent(inout)				    :: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)					        :: element

#           if defined(_FLASH_PATCH)
                type(t_state), dimension(_FLASH_PATCH_ORDER_SQUARE)   :: Q
                
                Q(:)%b = element%cell%data_pers%B

                call alpha_volume_op(traversal, section, element, Q)
                
                element%cell%data_pers%H = Q(:)%h
                element%cell%data_pers%HU = Q(:)%p(1)
                element%cell%data_pers%HV = Q(:)%p(2)
                
#           else
            type(t_state), dimension(_FLASH_CELL_SIZE)            :: Q

			call gv_Q%read(element, Q)

			call alpha_volume_op(traversal, section, element, Q)

			call gv_Q%write(element, Q)
#           endif
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine alpha_volume_op(traversal, section, element, Q)
			type(t_FLASH_init_dofs_traversal), intent(inout)				    :: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)						    :: element
#           if defined(_FLASH_PATCH)
                type(t_state), dimension(_FLASH_PATCH_ORDER_SQUARE), intent(out)  :: Q
                real (kind = GRID_SR), dimension(_FLASH_PATCH_ORDER_SQUARE)       :: max_wave_speed
                real (kind = GRID_SR), DIMENSION(2)                             :: r_coords     !< cell coords within patch
                integer (kind = GRID_SI)                                        :: j, row, col, cell_id
#           else
                type(t_state), dimension(_FLASH_CELL_SIZE), intent(out)   :: Q
                real (kind = GRID_SR)               :: max_wave_speed(_FLASH_CELL_SIZE)
#           endif

			real (kind = GRID_SR)               :: x(2)
			real (kind = GRID_SR)               :: dQ_norm
			real (kind = GRID_SR), parameter    :: probes(2, 3) = reshape([1.0, 0.0, 0.0, 0.0, 0.0, 1.0], [2, 3])
			type(t_dof_state)	                :: Q_test(3)
			integer                             :: i



			!evaluate initial DoFs (not the bathymetry!)
			
#           if defined(_FLASH_PATCH)
                row = 1
                col = 1
                do i=1, _FLASH_PATCH_ORDER_SQUARE
                
                    ! if orientation is backwards, the plotter uses a transformation that mirrors the cell...
                    ! this simple change solves the problem :)
                    if (element%cell%geometry%i_plotter_type > 0) then 
                        cell_id = i
                    else
                        cell_id = (row-1)*(row-1) + 2 * row - col
                    end if
                
                    r_coords = [0_GRID_SR, 0_GRID_SR]
                    do j=1,3
                        r_coords(:) = r_coords(:) + FLASH_PATCH_geometry%coords(:,j,cell_id) 
                    end do
                    r_coords = r_coords / 3
                    Q(i)%t_dof_state = get_initial_dof_state_at_position(section, samoa_barycentric_to_world_point(element%transform_data, r_coords))

                    col = col + 1
                    if (col == 2*row) then
                        col = 1
                        row = row + 1
                    end if
                end do
                
                ! dry cells
                where (element%cell%data_pers%H < element%cell%data_pers%B + cfg%dry_tolerance) 
                    element%cell%data_pers%H  = element%cell%data_pers%B
                    element%cell%data_pers%HU = 0.0_GRID_SR
                    element%cell%data_pers%HV = 0.0_GRID_SR
                end where
#           else			
			Q%t_dof_state = get_initial_dof_state_at_element(section, element)
			
            ! dry cells
            if (Q(1)%h < Q(1)%b + cfg%dry_tolerance) then
                Q%h  = Q%b
                Q(1)%p(:) = [0.0_GRID_SR, 0.0_GRID_SR]
            end if
#           endif

			element%cell%geometry%refinement = 0

			if (element%cell%geometry%i_depth < cfg%i_min_depth) then
                !refine if the minimum depth is not met

 				element%cell%geometry%refinement = 1
				traversal%i_refinements_issued = traversal%i_refinements_issued + 1
            else if (element%cell%geometry%i_depth < cfg%i_max_depth) then
#               if defined (_ASAGI)
                    !refine the displacements

                    x = samoa_barycentric_to_world_point(element%transform_data, [1.0_GRID_SR/3.0_GRID_SR, 1.0_GRID_SR/3.0_GRID_SR])
#                   if defined (_FLASH_PATCH)
                        dQ_norm = maxval(abs(get_bathymetry_at_patch(section, element, real(cfg%t_max_eq + 1.0, GRID_SR) ) - element%cell%data_pers%B))
#                   else
                        dQ_norm = abs(get_bathymetry_at_element(section, element, real(cfg%t_max_eq + 1.0, GRID_SR) ) - Q(1)%b)
#                   endif

                    if (dQ_norm > 100.0_GRID_SR *cfg%dry_tolerance) then
                        element%cell%geometry%refinement = 1
                        traversal%i_refinements_issued = traversal%i_refinements_issued + 1
                    end if
#               else
                    !refine any slopes in the initial state

                    do i = 1, 3
                        x = samoa_barycentric_to_world_point(element%transform_data, probes(:, i))
                        Q_test(i) = get_initial_dof_state_at_position(section, x)
                    end do

                    if (maxval(Q_test%h) > minval(Q_test%h)) then
                        element%cell%geometry%refinement = 1
                        traversal%i_refinements_issued = traversal%i_refinements_issued + 1
                    end if
#               endif
			end if

			!estimate initial time step

			where (Q%h - Q%b > 0.0_GRID_SR)
				max_wave_speed = sqrt(g * (Q%h - Q%b)) + sqrt((Q%p(1) * Q%p(1) + Q%p(2) * Q%p(2)) / ((Q%h - Q%b) * (Q%h - Q%b)))
			elsewhere
				max_wave_speed = 0.0_GRID_SR
			end where

            !This will cause a division by zero if the wave speeds are 0.
            !Bue to the min operator, the error will not affect the time step.
#           if defined _FLASH_PATCH
                section%r_dt_new = min(section%r_dt_new, element%cell%geometry%get_volume() / (sum(element%cell%geometry%get_edge_sizes()) * maxval(max_wave_speed) * _FLASH_PATCH_ORDER))
#           else
                section%r_dt_new = min(section%r_dt_new, element%cell%geometry%get_volume() / (sum(element%cell%geometry%get_edge_sizes()) * maxval(max_wave_speed)))
#           endif
		end subroutine

		function get_initial_dof_state_at_element(section, element) result(Q)
			type(t_grid_section), intent(inout)		:: section
			type(t_element_base), intent(inout)     :: element
			type(t_dof_state)						:: Q

			real (kind = GRID_SR)		            :: x(2)

            x = samoa_barycentric_to_world_point(element%transform_data, [1.0_GRID_SR / 3.0_GRID_SR, 1.0_GRID_SR / 3.0_GRID_SR])
            Q = get_initial_dof_state_at_position(section, x)
		end function

		function get_initial_dof_state_at_position(section, x) result(Q)
			type(t_grid_section), intent(inout)					:: section
			real (kind = GRID_SR), intent(in)		            :: x(:)            !< position in world coordinates
			type(t_dof_state)							        :: Q

            real (kind = GRID_SR)                               :: xs(2)

            xs = cfg%scaling * x + cfg%offset

#			if defined(_ASAGI)
				Q%h = 0.0_GRID_SR
#			else
                Q = FLASH_Scenario_get_initial_Q(xs)
#			endif

			Q%p = 0.0_GRID_SR
		end function
	END MODULE
#endif

