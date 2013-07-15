! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_NUMA)
	MODULE NUMA_Euler_Timestep
		use SFC_edge_traversal
		use Samoa_NUMA


		implicit none

		type num_traversal_data
		    integer (kind = GRID_SI)			:: i_refinements_issued
		end type

        interface skeleton_op
            module procedure skeleton_array_op
            module procedure skeleton_scalar_op
        end interface

        interface bnd_skeleton_op
            module procedure bnd_skeleton_array_op
            module procedure bnd_skeleton_scalar_op
        end interface

		PUBLIC cell_to_edge_op

		type(t_gv_Q)				        :: gv_Q
		type(t_lfs_flux)				:: lfs_flux

#		define _GT_NAME					t_numa_euler_timestep_traversal

#		define _GT_EDGES
#		define _GT_EDGES_TEMP

#		define _GT_PRE_TRAVERSAL_OP			pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP			post_traversal_op
#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op

#		define _GT_CELL_TO_EDGE_OP			cell_to_edge_op
#		define _GT_SKELETON_OP				skeleton_op
#		define _GT_BND_SKELETON_OP			bnd_skeleton_op
#		define _GT_CELL_UPDATE_OP			cell_update_op
#		define _GT_CELL_LAST_TOUCH_OP			cell_last_touch_op

#		include "SFC_generic_traversal_ringbuffer.f90"


		!*******************************
		!Geometry operators
		!*******************************
		subroutine pre_traversal_grid_op(traversal, grid)
			type(t_numa_euler_timestep_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							    :: grid
			_log_write(0, *) "**************************************"
			_log_write(0, *) (grid%u_max)
			_log_write(0, *) "**************************************"
			!grid%r_dt = 0.45_GRID_SR * get_edge_size(grid%d_max) / ((2.0_GRID_SR + sqrt(2.0_GRID_SR)) * grid%u_max)
			!print *, " pre_traversal_grid_op: dt= ", grid%r_dt
			grid%u_max = 1.0_GRID_SR
		end subroutine
   		
		subroutine post_traversal_grid_op(traversal, grid)
			type(t_numa_euler_timestep_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							    :: grid

		            call reduce(traversal%i_refinements_issued, traversal%children%i_refinements_issued, MPI_SUM, .true.)
		            call reduce(grid%u_max, grid%sections%elements_alloc%u_max, MPI_MAX, .true.)
			grid%r_time = grid%r_time + grid%r_dt
		end subroutine

		subroutine pre_traversal_op(traversal, section)
			type(t_numa_euler_timestep_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: section

			!this variable will be incremented for each cell with a refinement request
			traversal%i_refinements_issued = 0
			section%u_max = 0.0_GRID_SR
		end subroutine

		subroutine post_traversal_op(traversal, section)
			type(t_numa_euler_timestep_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							    :: section
		end subroutine

		function cell_to_edge_op(element, edge) result(rep)
			implicit none
			type(t_element_base), intent(in)						:: element
			type(t_edge_data), intent(in)						:: edge
			type(num_cell_rep)										:: rep

			type(t_state), dimension(_NUMA_CELL_SIZE)				:: Q
			integer(kind = GRID_SI)									:: i, j, i_edge
			real(kind = GRID_SR), dimension(2, _NUMA_EDGE_SIZE)		:: dof_pos
			real(kind = GRID_SR), dimension(2, 3), parameter		:: edge_offsets = reshape([0.0, 0.0, 0.0, 1.0, 1.0, 0.0], [2, 3])
			real(kind = GRID_SR), dimension(2, 3), parameter		:: edge_vectors = reshape([0.0, 1.0, 1.0, -1.0, -1.0, 0.0], [2, 3])

			call gv_Q%read(element, Q)
			
			i_edge = edge%transform_data%index

			_log_write(6, '(A, I0, A, F0.3, X, F0.3, X, F0.3, X, F0.3)') "NUMA: el->edge ", i_edge, " Q in: ", Q

			if (_NUMA_CELL_SIZE > 1) then
				forall (i = 1 : _NUMA_EDGE_SIZE)
					dof_pos(:, i) = edge_offsets(:, i_edge) + t_basis_flux_get_dof_coords(i) * edge_vectors(:, i_edge)
				end forall


				forall (i = 1 : _NUMA_EDGE_SIZE)
					rep%Q(i)%h = t_basis_Q_eval(dof_pos(:, i), Q%h)
					rep%Q(i)%p(1) = t_basis_Q_eval(dof_pos(:, i), Q%p(1))
					rep%Q(i)%p(2) = t_basis_Q_eval(dof_pos(:, i), Q%p(2))
				end forall
			else
				rep%Q(1) = Q(1)
			end if

			_log_write(6, '(A, I0, A, F0.3, X, F0.3, X, F0.3, X, F0.3)') "NUMA: el->edge ", i_edge, " rep out: ", rep%Q
		end function
		subroutine skeleton_array_op(traversal, grid, edges, rep1, rep2, update1, update2)
			type(t_numa_euler_timestep_traversal), intent(in)				:: traversal
			type(t_grid_section), intent(in)							    :: grid
			type(t_edge_data), intent(in)								    :: edges(:)
			type(num_cell_rep), intent(in)									:: rep1(:), rep2(:)
			type(num_cell_update), intent(out)								:: update1(:), update2(:)

            integer (kind = GRID_SI)                                        :: i

            do i = 1, size(edges)
                call skeleton_scalar_op(traversal, grid, edges(i), rep1(i), rep2(i), update1(i), update2(i))
            end do
		end subroutine

		subroutine skeleton_scalar_op(traversal, grid, edge, rep1, rep2, update1, update2)
			type(t_numa_euler_timestep_traversal), intent(in)				:: traversal
			type(t_grid_section), intent(in)							:: grid
			type(t_edge_data), intent(in)								:: edge
			type(num_cell_rep), intent(in)									:: rep1, rep2
			type(num_cell_update), intent(out)								:: update1, update2

			real (kind = GRID_SR), dimension(2)				:: press
			_log_write(6, '(A, F0.3, X, F0.3, X, F0.3, X, F0.3)') "NUMA: skel rep1 in: ", rep1
			_log_write(6, '(A, F0.3, X, F0.3, X, F0.3, X, F0.3)') "NUMA: skel rep2 in: ", rep2

				!print *, "NORMAL OP CALLED!"
				call compute_pressure(press(1),rep1%Q(1)%e,rep1%Q(1)%e_ref)
				call compute_pressure(press(2),rep2%Q(1)%e,rep2%Q(1)%e_ref)
				!_log_write(0, *) "Pressure: ", press
				call rusanov_flux_weak_nc(edge%transform_data%normal, rep1%Q(1), rep2%Q(1), update1%flux(1), update2%flux(1), press)


			_log_write(6, '(A, F0.3, X, F0.3, X, F0.3, X, F0.3)') "NUMA: skel fluxL out: ", update1
			_log_write(6, '(A, F0.3, X, F0.3, X, F0.3, X, F0.3)') "NUMA: skel fluxR out: ", update2
		end subroutine
		subroutine bnd_skeleton_array_op(traversal, grid, edges, rep, update)
				type(t_numa_euler_timestep_traversal), intent(in)				:: traversal
				type(t_grid_section), intent(in)							    :: grid
				type(t_edge_data), intent(in)								    :: edges(:)
				type(num_cell_rep), intent(in)									:: rep(:)
				type(num_cell_update), intent(out)								:: update(:)

		    integer (kind = GRID_SI)                                        :: i

		    do i = 1, size(edges)
		        call bnd_skeleton_scalar_op(traversal, grid, edges(i), rep(i), update(i))
		    end do
			end subroutine
	
		subroutine bnd_skeleton_scalar_op(traversal, grid, edge, rep, update)
			implicit none
			type(t_numa_euler_timestep_traversal), intent(in)				:: traversal
			type(t_grid_section), intent(in)							:: grid
			type(t_edge_data), intent(in)								:: edge
			type(num_cell_rep), intent(in)									:: rep
			type(num_cell_update), intent(out)								:: update
			type(t_state)													:: bnd_rep
			type(t_update)													:: bnd_flux
			
			real (kind = GRID_SR), dimension(2)				:: press

			bnd_rep = rep%Q(1)!t_state(0.0, [0.0, 0.0], 0.0, 0.0, [0.0, 0.0], 0.0)

			_log_write(6, '(A, F0.3, X, F0.3, X, F0.3, X, F0.3)') "NUMA: bnd skel rep in : ", rep

				!print *, "BND OP CALLED!"
				call compute_pressure(press(1),rep%Q(1)%e,rep%Q(1)%e_ref)
				call compute_pressure(press(2),rep%Q(1)%e,rep%Q(1)%e_ref)
				call rusanov_flux_weak_nc(edge%transform_data%normal, rep%Q(1), bnd_rep, update%flux(1), bnd_flux,press)


			_log_write(6, '(A, F0.3, X, F0.3, X, F0.3, X, F0.3)') "NUMA: bnd skel flux out: ", update
		end subroutine

		subroutine cell_update_op(traversal, grid, element, update1, update2, update3)
			implicit none
			type(t_numa_euler_timestep_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: grid
			type(t_element_base), intent(inout)						:: element
			type(num_cell_update), intent(in)						:: update1, update2, update3

			!local variables

			type(t_state), dimension(_NUMA_CELL_SIZE)				:: dQ

			call volume_op(traversal, grid, element, dQ, [update1%flux, update2%flux, update3%flux])

			call gv_Q%add(element, dQ)
		end subroutine

		subroutine cell_last_touch_op(traversal, grid, cell)
			implicit none
			type(t_numa_euler_timestep_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: grid
			type(t_cell_data_ptr), intent(inout)				:: cell
			integer (kind = 1)								:: depth
			real(kind = GRID_SR)							:: b_norm

			depth = cell%geometry%i_depth
			b_norm = minval(abs(cell%data_pers%Q%h))

			!refine also on the coasts
			if (depth < grid%i_max_depth .and. b_norm < 0.01_GRID_SR) then
				cell%geometry%refinement = 1
				traversal%i_refinements_issued = traversal%i_refinements_issued + 1
			else if (b_norm < 0.05_GRID_SR) then
				cell%geometry%refinement = 0
			endif			
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine volume_op(traversal, grid, element, dQ, fluxes)
			implicit none
			type(t_numa_euler_timestep_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: grid
			type(t_element_base), intent(inout)								:: element
			type(t_state), dimension(:), intent(out)						:: dQ
			type(t_update), dimension(:), intent(in)						:: fluxes

			real(kind = GRID_SR)											:: volume, dQ_norm, edge_lengths(3)
			integer (kind = 1)												:: i, depth

			_log_write(5, '(A, 4(X, F0.3))') "NUMA: volume op edge 1 flux in:", fluxes(1)
			_log_write(5, '(A, 4(X, F0.3))') "NUMA: volume op edge 2 flux in:", fluxes(2)
			_log_write(5, '(A, 4(X, F0.3))') "NUMA: volume op edge 3 flux in:", fluxes(3)

			volume = element%cell%geometry%get_volume()
			edge_lengths = [element%cell%geometry%get_leg_size(), element%cell%geometry%get_hypo_size(), element%cell%geometry%get_leg_size()]
			!_log_write(0, *) "before: flux_h: ",fluxes%h
			!_log_write(0, *) "before: h: ",dQ%h,"e: ", dQ%e
			!_log_write(0, *) "--------------------------------------------"
			dQ%h = sum(edge_lengths * fluxes%h)
			dQ%p(1) = sum(edge_lengths * fluxes%p(1))
			dQ%p(2) = sum(edge_lengths * fluxes%p(2))
			dQ%e = sum(edge_lengths * fluxes%e)

			!set refinement condition

			element%cell%geometry%refinement = 0
			dQ_norm = dot_product(dQ%e , dQ%e)

			depth = element%cell%geometry%i_depth
			if (depth < grid%i_max_depth .and. dQ_norm > (1.0e-5_GRID_SR ** 2)) then
				element%cell%geometry%refinement = 1
				traversal%i_refinements_issued = traversal%i_refinements_issued + 1
			else if (depth > grid%i_min_depth .and. dQ_norm < (5.0e-7_GRID_SR ** 2)) then
				element%cell%geometry%refinement = -1
			endif			
			grid%u_max = max(grid%u_max, maxval(fluxes%max_wave_speed))
			dQ%t_dof_state = (-grid%r_dt / volume) * dQ%t_dof_state
			!print *, fluxes%p(1),fluxes%h
			!print *, "VOLUME OPERATOR!!!!"
			!_log_write(0, *) "after: h: ",dQ%h,"e: ", dQ%e
			!_log_write(0, *) "*********************************************"
			_log_write(6, '(A, 4(X, F0.3))') "NUMA: volume op out: ", dQ
		end subroutine

		!----------------- NUMA

		subroutine rusanov_flux_weak_nc(normal, QL, QR, fluxL, fluxR,r)
		  
		  
		  use NUMA_constants, only: gravity, p00, cv, rgas, gamma

		  implicit none
			type(t_state), intent(in)						:: QL, QR
			type(t_update), intent(out) 						:: fluxL, fluxR
			real(kind = GRID_SR), dimension(2), intent(in)				:: normal
			real (kind = GRID_SR), dimension(2)				::r		
		  !local
		  real (kind = GRID_SR) nxl,nzl,nxr,nzr
		  real (kind = GRID_SR) normx_r1,normz_r1,normx_r2,normz_r2
		  real (kind = GRID_SR) pi, temp
		  real (kind = GRID_SR) rl_k, ul_k, wl_k, tl_k, pl_k, al_k
		  real (kind = GRID_SR) rul_k, rwl_k, rtl_k
		  real (kind = GRID_SR) rr_k, ur_k, wr_k, tr_k, pr_k, ar_k
		  real (kind = GRID_SR) rur_k, rwr_k, rtr_k
		  real (kind = GRID_SR) unl, unr, utl, utr, claml, clamr, clam_glob, clam
		  real (kind = GRID_SR) fxl, fzl
		  real (kind = GRID_SR) fxr, fzr
		  real (kind = GRID_SR) flux_rl, flux_rr, diss_r
		  real (kind = GRID_SR) flux_ul, flux_ur, diss_u
		  real (kind = GRID_SR) flux_wl, flux_wr, diss_w
		  real (kind = GRID_SR) flux_tl, flux_tr, diss_t
		  integer i
		!  print *, "normal: ",normal
		!  print *, "QL: ",QL

		  !Project normal and tangent vectors

	!	print *, "FLUX FUNCTION RGAS: ",rgas

		  !Do Gauss-Lobatto Integration
	
		     !Store Normal Vectors
		     nxl=normal(1)
		     nzl=normal(2)
		     
		     nxr=-nxl
		     nzr=-nzl
	!	     print *, "EULER: ",rgas
		     !Interpolate onto Quadrature Points
		     rl_k=QL%h + QL%h_ref
		     ul_k=QL%p(1)/rl_k
		     wl_k=QL%p(2)/rl_k
		     tl_k=(QL%e + QL%e_ref)/rl_k
		     pl_k=r(1)
		     rul_k=rl_k*ul_k
		     rwl_k=rl_k*wl_k
		     rtl_k=rl_k*tl_k
		     
		     rr_k=QR%h + QR%h_ref
		 !    print *, rr_k ,QR%h , QR%h_ref
		     ur_k=QR%p(1)/rr_k
		     wr_k=QR%p(2)/rr_k
		     tr_k=(QR%e + QR%e_ref)/rr_k
		     pr_k=r(2)
		     rur_k=rr_k*ur_k
		     rwr_k=rr_k*wr_k
		     rtr_k=rr_k*tr_k
		!     print *, "P00: ",p00,cv
		     
		     !Compute Rusanov flux Constant
		     unl=nxl*ul_k + nzl*wl_k
		     unr=nxl*ur_k + nzl*wr_k
		     pi=(rl_k*rgas*tl_k/p00)**(rgas/cv)
		     temp=tl_k*pi
		     al_k=sqrt(gamma*rgas*temp)
		     pi=(rr_k*rgas*tr_k/p00)**(rgas/cv)
		     temp=tr_k*pi
		     !if (temp < 0) then
			!print *,  tr_k, QR%e ,QR%e_ref,rr_k,QR%h , QR%h_ref
			!stop
		     !endif
		     ar_k=sqrt(gamma*rgas*temp)
		     !print *,  ar_k,rgas,temp,gamma
		     claml=abs(unl) + al_k
		     clamr=abs(unr) + ar_k
		     clam=max(claml,clamr)
		     
		     !----Mass Equation-----!
		     !Flux Variables
		     fxl=rul_k
		     fzl=rwl_k
		     
		     fxr=rur_k
		     fzr=rwr_k
		     
		     !Normal Flux Component
		     flux_rl=( nxl*fxl + nzl*fzl )
		     flux_rr=( nxr*fxr + nzr*fzr )
		     fluxL%h=flux_rl - flux_rr
		     
		     !Dissipation Term
		     diss_r=clam*(rr_k - rl_k)
		     
		     !----U-Momentum Equation-----!
		     !Flux Variables
		     fxl=rul_k*ul_k + pl_k
		     fzl=rul_k*wl_k
		     
		     fxr=rur_k*ur_k + pr_k
		     fzr=rur_k*wr_k 
		     
		     !Normal Flux Component
		     flux_ul=( nxl*fxl + nzl*fzl )
		     flux_ur=( nxr*fxr + nzr*fzr )
		     fluxL%p(1)=flux_ul - flux_ur
		     
		     !Dissipation Term
		     diss_u=clam*(rur_k - rul_k)
		     !print *,"DISS: ",diss_u
		     !----W-Momentum Equation-----!
		     !Flux Variables
		     fxl=rwl_k*ul_k 
		     fzl=rwl_k*wl_k + pl_k
		     
		     fxr=rwr_k*ur_k 
		     fzr=rwr_k*wr_k + pr_k
		     
		     !Left Normal Flux Component
		     flux_wl=( nxl*fxl + nzl*fzl )
		     flux_wr=( nxr*fxr + nzr*fzr )
		     fluxL%p(2)=flux_wl - flux_wr
		     
		     !Dissipation Term
		     diss_w=clam*(rwr_k - rwl_k)
		     
		     !----Temperature Equation-----!
		     !Flux Variables
		     fxl=rtl_k*ul_k 
		     fzl=rtl_k*wl_k
		     
		     fxr=rtr_k*ur_k 
		     fzr=rtr_k*wr_k
		     
		     !Left Normal Flux Component
		     flux_tl=( nxl*fxl + nzl*fzl )
		     flux_tr=( nxr*fxr + nzr*fzr )
		     fluxL%e=flux_tl - flux_tr
		     
		     !Dissipation Term
		     diss_t=clam*(rtr_k - rtl_k)
		     
		     !Construct Rusanov Flux
		     fluxL%h=0.5*( fluxL%h- diss_r)
		     fluxL%p(1)=0.5*( fluxL%p(1)- diss_u)
		     fluxL%p(2)=0.5*(fluxL%p(2)- diss_w)
		     fluxL%e=0.5*(fluxL%e - diss_t)


		     fluxR = t_update_state_inv(fluxL)
		     !fluxR=fluxL
		!   print * , "FLUXR: ",fluxR
		  end subroutine rusanov_flux_weak_nc

		subroutine compute_pressure(press,e,e_ref)

			  use NUMA_constants, only: rgas, p00, cp, cv

			  implicit none
			  real (kind = GRID_SR)  press,e,e_ref
			  real (kind = GRID_SR)  c, rho, theta
			  real (kind = GRID_SR)  p_ref

			  !constant
			  c=cp/cv
			! _log_write(0, *) "*****************E_REF: ", e, e_ref
       		          theta=e + e_ref
			  call compute_pressure_ref(p_ref,e_ref)
			 ! print *, press, p_ref, c,p00 ,theta ,rgas, e , e_ref
			 !_log_write(0, *) "----------------P_REF: ", p_ref
			 !_log_write(0, *) "****************THETA: ", theta
			  press=p00*(rgas*theta/p00)**c - p_ref
		end subroutine compute_pressure

		subroutine compute_pressure_ref(press_ref,e_ref)

		  use NUMA_constants, only: rgas, p00, cp, cv

		  implicit none

		 real (kind = GRID_SR)  press_ref,e_ref
		  !local variables
		  real (kind = GRID_SR) c, rho, theta
		  integer i, k, e, ie

		  c=cp/cv
	!	  print *,rgas, p00, cp, cv
			   theta=e_ref
			   press_ref=p00*(rgas*theta/p00)**c

		end subroutine compute_pressure_ref

		!----------------- NUMA

		!> Augmented Riemann solver from geoclaw
	END MODULE
#endif
