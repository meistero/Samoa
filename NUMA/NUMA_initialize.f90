! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_NUMA)
	MODULE NUMA_Initialize
 		!use Tools_noise
		use SFC_edge_traversal
		use NUMA_euler_timestep

		use Samoa_NUMA

		implicit none

		type num_traversal_data
		    integer (kind = GRID_SI)			:: i_refinements_issued
		end type

		type(t_gv_Q)							:: gv_Q
		type(t_lfs_flux)						:: lfs_flux

#		define _GT_NAME							t_numa_init_traversal

#		define _GT_EDGES
#		define _GT_EDGES_TEMP

#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op
#		define _GT_ELEMENT_OP					element_op

#		define _GT_CELL_TO_EDGE_OP				cell_to_edge_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine pre_traversal_op(traversal, section)
			type(t_numa_init_traversal), intent(inout)				    :: traversal
			type(t_grid_section), intent(inout)							:: section

			!this variable will be incremented for each cell with a refinement request
			traversal%i_refinements_issued = 0
			section%u_max = 1.0_GRID_SR
		end subroutine

		subroutine post_traversal_grid_op(traversal, grid)
			type(t_numa_init_traversal), intent(inout)		        :: traversal
			type(t_grid), intent(inout)							    :: grid

		    call reduce(traversal%i_refinements_issued, traversal%children%i_refinements_issued, MPI_SUM, .true.)
		    call reduce(grid%u_max, grid%sections%elements_alloc%u_max, MPI_MAX, .true.)
		end subroutine

		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, section, element)
			type(t_numa_init_traversal), intent(inout)				    :: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)					:: element

			type(t_state), dimension(_NUMA_CELL_SIZE)			:: Q

			call alpha_volume_op(traversal, section, element, Q)

			call gv_Q%write(element, Q)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine alpha_volume_op(traversal, section, element, Q)
			type(t_numa_init_traversal), intent(inout)				    :: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)						:: element
			type(t_state), dimension(_NUMA_CELL_SIZE), intent(out)	:: Q

			real (kind = GRID_SR), dimension(2)						:: pos
			integer (kind = GRID_SI)								:: i
			real (kind = GRID_SR), parameter, dimension(2, 3)		:: r_test_points = reshape([1.0, 0.0, 0.0, 0.0, 0.0, 1.0], [2, 3])
			type(t_state), dimension(3)								:: Q_test
			real (kind = GRID_SR), dimension(_NUMA_CELL_SIZE)		:: lambda

			!evaluate initial function values at dof positions and compute DoFs

			do i = 1, _NUMA_CELL_SIZE
				Q(i) = get_initial_state(section, samoa_barycentric_to_world_point(element%transform_data, t_basis_Q_get_dof_coords(i)), element%cell%geometry%i_depth / 2_GRID_SI)
			!	print *,"IN THE BEGINNING: ", t_basis_Q_get_dof_coords(i)
			!	print *,"IN THE SECOND BEGINNING: ", samoa_barycentric_to_world_point(element%transform_data, t_basis_Q_get_dof_coords(i))
			end do

			if (_NUMA_ORDER == 0) then
				do i = 1, 3
			!		Q_test(i) = get_initial_state(grid, samoa_barycentric_to_world_point(element%transform_data, r_test_points(:, i)), element%cell%geometry%i_depth / 2_GRID_SI)
			!		print *,"IN THE test: ", r_test_points(:, i), samoa_barycentric_to_world_point(element%transform_data, r_test_points(:, i))
				end do
			else
				do i = 1, 3
					Q_test(i)%h = t_basis_Q_eval(r_test_points(:, i), Q%h)
				end do
			endif
!			print *, "ALPHA VOLUME OP BEFORE REFINEMENT"
			if (element%cell%geometry%i_depth < section%i_min_depth ) then!.or. (element%cell%geometry%i_depth < section%i_max_depth .and. ( &
!					(maxval(Q_test%h) > 0.0_GRID_SR .and. minval(Q_test%h) < 0.0_GRID_SR) .or. &			!refine coast lines
!					abs(Q_test(3)%h - Q_test(2)%h) + abs(Q_test(1)%h - Q_test(2)%h) > 0.0_GRID_SR))) then						!refine waves
!			print *, "ALPHA VOLUME OP AFTER REFINEMENT"
				element%cell%geometry%refinement = 1
				traversal%i_refinements_issued = traversal%i_refinements_issued + 1
			else
				element%cell%geometry%refinement = 0
			end if

			!estimate initial u_max

			where (Q%h> 0)
				lambda = sqrt(g * (Q%h )) + sqrt((Q%p(1) * Q%p(1) + Q%p(2) * Q%p(2)) / ((Q%h ) * (Q%h)))
			elsewhere
				lambda = 0.0_GRID_SR
			end where

			section%u_max = max(section%u_max, maxval(lambda))
			
		end subroutine

		function get_initial_state(grid, x, lod) result(Q)
			type(t_grid_section), intent(inout)							:: grid
			real (kind = GRID_SR), dimension(:), intent(in)		:: x
			integer (kind = GRID_SI), intent(in)				:: lod						!< level of detail
			type(t_state)										:: Q

			!init height, momentum and bathymetry
			Q%t_dof_state = get_initial_dof_state(grid, x, lod)
		end function

		function get_initial_dof_state(grid, x, lod) result(Q)

			use NUMA_constants, only: rgas, p00, cp, cv, gravity, gamma, xappa, omega, earth_radius

			type(t_grid_section), intent(inout)							:: grid
			real (kind = GRID_SR), dimension(:), intent(in)		:: x						!< position in world coordinates
			integer (kind = GRID_SI), intent(in)				:: lod						!< level of detail
			type(t_dof_state)									:: Q						!< initial state

			real (kind = GRID_SR), dimension(2), parameter		:: dam_center = [0.5, 0.5]
			real (kind = GRID_SR), parameter					:: dam_radius = 0.02
			real (kind = GRID_SR), parameter					:: outer_height = 0.0
			real (kind = GRID_SR), parameter					:: inner_height = 0.2
			  !Local Arrays
			  real pi, f0
			  real sigma, sigma2, u, v, w, r2
			  real xx, yy, zz
			  real temp0, c, c2, u0, v0, w0, bv, bv2, g2, ac, xl, yl, zl, xc, yc, zc, thetac
			  real u_k, v_k, w_k, pi_ref, rho_ref, theta_ref, dtheta, theta_k, rho_k
			  real theta0, rc, r, pi_k, xr, yr, zr
			  real rho00, hc, n2, l2, l, lz, p_k, rhofac, afac, hafac, delta_k
			  real ddelta_dx_k, ddelta_dy_k, ddelta_dz_k
			  real rho_exact, pi_exact, u_exact, v_exact, w_exact, theta_exact
			  real rho_pert, p_pert, t_pert, theta_pert, rho
			  real xf, yf, zf, radius, zradius
			  integer e, i, j, k
			  real gasr
			gasr= 287.17
			  !Constants
			  pi=4.0*atan(1.0)

			  !initialize

		! print *,"the solution for position ",x,"are: "

		     theta0=300
		     u0=0
		     w0=0
		     xc=0.50
		     zc=0.35
		     thetac=2
		     rc=.20
		    pi=4.0*atan(1.0)
		    earth_radius=160
		    omega=7.29e-5
		    gravity=9.80616
		    gamma=1.4
		    cp=1004.67
		    cv=717.5
		    rgas=287.17
		    xappa=0.286
		    p00=1e5
		    c=cv/gasr
		     c2=gasr/cp

			r = sqrt( (x(1)-xc)**2 + (x(2)-zc)**2 )
		 dtheta=0
                 if (r < rc) then
                     dtheta=thetac *(1 + cos(pi*r/rc) )
	!	     print *, "IF CONDITION HAPPENED!", dtheta
                 end if
		theta_k=theta0 + dtheta
      	        pi_k=1 - gravity/(cp*theta0)*x(2)
      	        rho_k=p00/(gasr*theta_k)*(pi_k)**c


		Q%h = rho_k
		Q%p(1) = rho_k * u0
		Q%p(2) = rho_k * w0
		Q%e = (theta_k)
	!	 if (Q%e < 348) then
	!		print *,"state: ", Q%e,Q%h, dtheta,gasr
	!	endif
  	        theta_ref=300
  	        pi_ref=1 - gravity/(cp*theta_ref)*x(2)
  	        rho_ref=p00/(gasr*theta_ref)*(pi_ref)**c
		Q%h_ref = rho_ref
		Q%p_ref (1) = rho_ref * u0 +10
		Q%p_ref (2) = rho_ref * w0 +2.5
		Q%e_ref  = (theta_ref)
		!print *, Q%e_ref,rho_ref,theta_ref
		!print *,"TERMS: ", Q%e_ref
		Q%e = (Q%e - Q%e_ref)*rho_k
		Q%h = -(Q%h - Q%h_ref)
		end function
	END MODULE
#endif
