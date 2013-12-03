! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
	MODULE SWE_Euler_Timestep
		use SFC_edge_traversal

		use Samoa_swe
		use c_bind_riemannsolvers

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

		type(t_gv_Q)							:: gv_Q
		type(t_lfs_flux)						:: lfs_flux

#		define _GT_NAME							t_swe_euler_timestep_traversal

#		define _GT_EDGES
#		define _GT_EDGES_TEMP

#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP			post_traversal_op
#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op

#		define _GT_CELL_TO_EDGE_OP				cell_to_edge_op
#		define _GT_SKELETON_OP					skeleton_op
#		define _GT_BND_SKELETON_OP				bnd_skeleton_op
#		define _GT_CELL_UPDATE_OP				cell_update_op
#		define _GT_CELL_LAST_TOUCH_OP			cell_last_touch_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		!*******************************
		!Geometry operators
		!*******************************

		subroutine pre_traversal_grid_op(traversal, grid)
			type(t_swe_euler_timestep_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							    :: grid

            grid%r_dt = 0.45_GRID_SR * grid%scaling * get_edge_size(grid%d_max) / ((2.0_GRID_SR + sqrt(2.0_GRID_SR)) * grid%u_max)

#           if defined(_ASAGI)
                if (grid%r_time < grid_max_z(grid%afh_bathymetry)) then
                    grid%r_dt = min(grid%r_dt, 0.1/15.0 * grid_max_z(grid%afh_bathymetry))
                end if
#           endif

			grid%u_max = 0.0_GRID_SR
		end subroutine

		subroutine post_traversal_grid_op(traversal, grid)
			type(t_swe_euler_timestep_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							    :: grid

            call reduce(traversal%i_refinements_issued, traversal%children%i_refinements_issued, MPI_SUM, .true.)
            call reduce(grid%u_max, grid%sections%elements_alloc%u_max, MPI_MAX, .true.)
			grid%r_time = grid%r_time + grid%r_dt

			grid%sections%elements%r_time = grid%r_time
		end subroutine

		subroutine pre_traversal_op(traversal, section)
			type(t_swe_euler_timestep_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: section

			!this variable will be incremented for each cell with a refinement request
			traversal%i_refinements_issued = 0
			section%u_max = 0.0_GRID_SR
		end subroutine

		subroutine post_traversal_op(traversal, section)
			type(t_swe_euler_timestep_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							    :: section
		end subroutine

		function cell_to_edge_op(element, edge) result(rep)
			type(t_element_base), intent(in)						:: element
			type(t_edge_data), intent(in)						    :: edge
			type(num_cell_rep)										:: rep

			type(t_state), dimension(_SWE_CELL_SIZE)				:: Q
			integer(kind = GRID_SI)									:: i, j, i_edge
			real(kind = GRID_SR), dimension(2, _SWE_EDGE_SIZE)		:: dof_pos
			real(kind = GRID_SR), dimension(2, 3), parameter		:: edge_offsets = reshape([0.0, 0.0, 0.0, 1.0, 1.0, 0.0], [2, 3])
			real(kind = GRID_SR), dimension(2, 3), parameter		:: edge_vectors = reshape([0.0, 1.0, 1.0, -1.0, -1.0, 0.0], [2, 3])

			call gv_Q%read(element, Q)

			_log_write(6, '(3X, A)') "swe cell to edge op:"
			_log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q in: ", Q
            _log_write(6, '(4X, A, F0.3, 1X, F0.3)') "normal in : ", edge%transform_data%normal

#           if (_SWE_CELL_SIZE > 1)
                i_edge = edge%transform_data%index
                _log_write(6, '(4X, A, I0)') "edge ", i_edge

				forall (i = 1 : _SWE_EDGE_SIZE)
					dof_pos(:, i) = edge_offsets(:, i_edge) + t_basis_flux_get_dof_coords(i) * edge_vectors(:, i_edge)
				end forall

				call lfs_flux%transform(edge%transform_data, dof_pos(1, :))
				call lfs_flux%transform(edge%transform_data, dof_pos(2, :))

				forall (i = 1 : _SWE_EDGE_SIZE)
					rep%Q(i)%h = t_basis_Q_eval(dof_pos(:, i), Q%h)
					rep%Q(i)%p(1) = t_basis_Q_eval(dof_pos(:, i), Q%p(1))
					rep%Q(i)%p(2) = t_basis_Q_eval(dof_pos(:, i), Q%p(2))
				end forall
#           else
				rep%Q(1) = Q(1)
#           endif

			_log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q out: ", rep%Q
		end function

		subroutine skeleton_array_op(traversal, grid, edges, rep1, rep2, update1, update2)
			type(t_swe_euler_timestep_traversal), intent(in)				:: traversal
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
			type(t_swe_euler_timestep_traversal), intent(in)				:: traversal
			type(t_grid_section), intent(in)							    :: grid
			type(t_edge_data), intent(in)								    :: edge
			type(num_cell_rep), intent(in)									:: rep1, rep2
			type(num_cell_update), intent(out)								:: update1, update2

			_log_write(6, '(3X, A)') "swe skeleton op:"
			_log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q 1 in: ", rep1%Q
			_log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q 2 in: ", rep2%Q

#			if defined (_SWE_LF) .or.  defined (_SWE_LF_BATH) .or.  defined (_SWE_LLF) .or.  defined (_SWE_LLF_BATH)
				call compute_lf_flux(edge%transform_data%normal, rep1%Q(1), rep2%Q(1), update1%flux(1), update2%flux(1))
#			else
				call compute_geoclaw_flux(edge%transform_data%normal, rep1%Q(1), rep2%Q(1), update1%flux(1), update2%flux(1))
#			endif

			_log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "flux 1 out: ", update1%flux
			_log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "flux 2 out: ", update2%flux
		end subroutine

		subroutine bnd_skeleton_array_op(traversal, grid, edges, rep, update)
			type(t_swe_euler_timestep_traversal), intent(in)				:: traversal
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
			type(t_swe_euler_timestep_traversal), intent(in)				:: traversal
			type(t_grid_section), intent(in)							    :: grid
			type(t_edge_data), intent(in)								    :: edge
			type(num_cell_rep), intent(in)									:: rep
			type(num_cell_update), intent(out)								:: update

			type(t_state)													:: bnd_rep
			type(t_update)													:: bnd_flux

			bnd_rep = t_state(0.0, [0.0, 0.0], rep%Q(1)%b)

#			if defined (_SWE_LF) .or.  defined (_SWE_LF_BATH) .or.  defined (_SWE_LLF) .or.  defined (_SWE_LLF_BATH)
				call compute_lf_flux(edge%transform_data%normal, rep%Q(1), bnd_rep, update%flux(1), bnd_flux)
#			else
				call compute_geoclaw_flux(edge%transform_data%normal, rep%Q(1), bnd_rep, update%flux(1), bnd_flux)
#			endif
		end subroutine

		subroutine cell_update_op(traversal, grid, element, update1, update2, update3)
			type(t_swe_euler_timestep_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: grid
			type(t_element_base), intent(inout)						:: element
			type(num_cell_update), intent(in)						:: update1, update2, update3

			!local variables

			type(t_state), dimension(_SWE_CELL_SIZE)				:: dQ

			call volume_op(traversal, grid, element, dQ, [update1%flux, update2%flux, update3%flux])

			call gv_Q%add(element, dQ)
		end subroutine

		subroutine cell_last_touch_op(traversal, grid, cell)
			type(t_swe_euler_timestep_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: grid
			type(t_cell_data_ptr), intent(inout)				:: cell
			integer (kind = 1)								:: depth
			real(kind = GRID_SR)							:: b_norm

			depth = cell%geometry%i_depth
			b_norm = minval(abs(cell%data_pers%Q%h - cell%data_pers%Q%b))

			!refine also on the coasts
			if (depth < grid%i_max_depth .and. b_norm < 100.0_GRID_SR) then
				cell%geometry%refinement = 1
				traversal%i_refinements_issued = traversal%i_refinements_issued + 1
			else if (b_norm < 300.0_GRID_SR) then
				cell%geometry%refinement = max(cell%geometry%refinement, 0)
			endif
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine volume_op(traversal, section, element, dQ, fluxes)
			type(t_swe_euler_timestep_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							    :: section
			type(t_element_base), intent(inout)								:: element
			type(t_state), dimension(:), intent(out)						:: dQ
			type(t_update), dimension(:), intent(in)						:: fluxes

			real(kind = GRID_SR)											:: volume, dQ_norm, edge_lengths(3)
			integer (kind = 1)												:: i, depth

			_log_write(6, '(3X, A)') "swe cell update op:"
			_log_write(6, '(4X, A, 4(X, F0.3))') "edge 1 flux in:", fluxes(1)
			_log_write(6, '(4X, A, 4(X, F0.3))') "edge 2 flux in:", fluxes(2)
			_log_write(6, '(4X, A, 4(X, F0.3))') "edge 3 flux in:", fluxes(3)

			volume = section%scaling * section%scaling * element%cell%geometry%get_volume()
			edge_lengths = section%scaling * element%cell%geometry%get_edge_sizes()

			dQ%h = sum(edge_lengths * fluxes%h)
			dQ%p(1) = sum(edge_lengths * fluxes%p(1))
			dQ%p(2) = sum(edge_lengths * fluxes%p(2))
			dQ%b = 0.0_GRID_SR

			!set refinement condition

			element%cell%geometry%refinement = 0
			dQ_norm = dot_product(dQ(1)%p, dQ(1)%p)

			depth = element%cell%geometry%i_depth
			if (depth < section%i_max_depth .and. dQ_norm > (section%scaling * 2.0_GRID_SR) ** 2) then
				element%cell%geometry%refinement = 1
				traversal%i_refinements_issued = traversal%i_refinements_issued + 1
			else if (depth > section%i_min_depth .and. dQ_norm < (section%scaling * 1.0_GRID_SR) ** 2) then
				element%cell%geometry%refinement = -1
			endif

			section%u_max = max(section%u_max, maxval(fluxes%max_wave_speed))

			dQ%t_dof_state = dQ%t_dof_state * (-section%r_dt / volume)

			_log_write(6, '(4X, A, 4(X, F0.3))') "dQ out: ", dQ
		end subroutine

		!> Lax Friedrichs flux. Depending on compiler flags, the function implements
		!> the global or local variant with or without bathymetry
		pure subroutine compute_lf_flux(normal, QL, QR, fluxL, fluxR)
			type(t_state), intent(in)							:: QL, QR
			type(t_update), intent(out) 						:: fluxL, fluxR
			real(kind = GRID_SR), intent(in)		            :: normal(2)

			real(kind = GRID_SR), parameter						:: dry_tol = 0.01_GRID_SR
			real(kind = GRID_SR)								:: vL, vR, alpha

#           if defined(_SWE_LF_BATH) .or. defined(_SWE_LLF_BATH)
                if (QL%h - QL%b < dry_tol) then
                    vL = 0.0_GRID_SR
                    fluxL%max_wave_speed = 0.0_GRID_SR
                else
                    vL = DOT_PRODUCT(normal, QL%p / (QL%h - QL%b))
                    fluxL%max_wave_speed = sqrt(g * (QL%h - QL%b)) + sqrt(vL * vL)
                end if

                if (QR%h - QR%b < dry_tol) then
                    vR = 0.0_GRID_SR
                    fluxR%max_wave_speed = 0.0_GRID_SR
                else
                    vR = DOT_PRODUCT(normal, QR%p / (QR%h - QR%b))
                    fluxR%max_wave_speed = sqrt(g * (QR%h - QR%b)) + sqrt(vR * vR)
                end if

#               if defined(_SWE_LLF_BATH)
                    alpha = max(fluxL%max_wave_speed, fluxR%max_wave_speed)
#               else
                    alpha = 100.0
#               endif

                fluxL%h = 0.5_GRID_SR * (vL * (QL%h - QL%b) + vR * (QR%h - QR%b) + alpha * (QL%h - QR%h))
                fluxR%h = -fluxL%h

                fluxL%p = 0.5_GRID_SR * ((vL + alpha) * QL%p + (vR - alpha) * QR%p) + 0.5_GRID_SR * g * (max(QR%h - QL%b, 0.0_GRID_SR) ** 2) * normal
                fluxR%p = -0.5_GRID_SR * ((vL + alpha) * QL%p + (vR - alpha) * QR%p) - 0.5_GRID_SR * g * (max(QL%h - QR%b, 0.0_GRID_SR) ** 2) * normal
#           else
                real(kind = GRID_SR), parameter					:: b = -1000.0_GRID_SR     !default constant bathymetry

                !use the height of the water pillars for computation
                vL = DOT_PRODUCT(normal, QL%p / (QL%h - b))
                vR = DOT_PRODUCT(normal, QR%p / (QR%h - b))

                fluxL%max_wave_speed = sqrt(g * (QL%h - b)) + sqrt(vL * vL)
                fluxR%max_wave_speed = sqrt(g * (QR%h - b)) + sqrt(vR * vR)

#               if defined(_SWE_LLF)
                    alpha = max(fluxL%max_wave_speed, fluxR%max_wave_speed)
#               else
                    alpha = 100.0
#               endif

                fluxL%h = 0.5_GRID_SR * (vL * (QL%h - b) + vR * (QR%h - b) + alpha * (QL%h - QR%h))
                fluxR%h = -fluxL%h

                fluxL%p = 0.5_GRID_SR * ((vL + alpha) * QL%p + (vR - alpha) * QR%p) + 0.5_GRID_SR * g * (QR%h - b) * (QR%h - b) * normal
                fluxR%p = -0.5_GRID_SR * ((vL + alpha) * QL%p + (vR - alpha) * QR%p) - 0.5_GRID_SR * g * (QL%h - b) * (QL%h - b) * normal
#           endif
		end subroutine

		!> Augmented Riemann solver from geoclaw
		subroutine compute_geoclaw_flux(normal, QL, QR, fluxL, fluxR)
			type(t_state), intent(in)           :: QL, QR
			type(t_update), intent(out)         :: fluxL, fluxR
			real(kind = GRID_SR), intent(in)    :: normal(2)

			real(kind = GRID_SR)				:: transform_matrix(2, 2)
			real(kind = GRID_SR)			    :: net_updatesL(3), net_updatesR(3), max_wave_speed
			real(kind = GRID_SR)                :: pL(2), pR(2), hL, hR, bL, bR

			transform_matrix(1, :) = normal
			transform_matrix(2, :) = [-normal(2), normal(1)]

			pL = matmul(transform_matrix, QL%p)
			pR = matmul(transform_matrix, QR%p)
			hL = QL%h - QL%b
			hR = QR%h - QR%b
			bL = QL%b
			bR = QR%b

#           if defined(_SWE_FWAVE)
                call c_bind_geoclaw_solver(GEOCLAW_FWAVE, 1, 3, hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, 0.01_8, g, net_updatesL, net_updatesR, max_wave_speed)
#           elif defined(_SWE_SSQ_FWAVE)
                call c_bind_geoclaw_solver(GEOCLAW_SSQ_FWAVE, 1, 3, hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, 0.01_8, g, net_updatesL, net_updatesR, max_wave_speed)
#           elif defined(_SWE_AUG_RIEMANN)
                call c_bind_geoclaw_solver(GEOCLAW_AUG_RIEMANN, 1, 3, hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, 0.01_8, g, net_updatesL, net_updatesR, max_wave_speed)
#           endif

			fluxL%h = net_updatesL(1)
			fluxL%p = matmul(net_updatesL(2:3), transform_matrix)
			fluxL%max_wave_speed = max_wave_speed

			fluxR%h = net_updatesR(1)
			fluxR%p = matmul(net_updatesR(2:3), transform_matrix)
			fluxR%max_wave_speed = max_wave_speed
		end subroutine
	END MODULE
#endif
