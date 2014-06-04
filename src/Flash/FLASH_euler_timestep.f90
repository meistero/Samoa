#include "Compilation_control.f90"

#if defined(_FLASH)
  MODULE FLASH_Euler_Timestep
    use SFC_edge_traversal

    use Samoa_FLASH

    use FLASH_dg_element

    implicit none

    type num_traversal_data
      integer (kind = GRID_SI)      :: i_refinements_issued
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

    type(t_gv_Q)                :: gv_Q
    type(t_lfs_flux)            :: lfs_flux

#   define _GT_NAME                     t_FLASH_euler_timestep_traversal

#   define _GT_EDGES
#   define _GT_EDGES_TEMP

#   define _GT_PRE_TRAVERSAL_OP         pre_traversal_op
#   define _GT_POST_TRAVERSAL_OP        post_traversal_op
#   define _GT_PRE_TRAVERSAL_GRID_OP    pre_traversal_grid_op
#   define _GT_POST_TRAVERSAL_GRID_OP   post_traversal_grid_op

#   define _GT_CELL_TO_EDGE_OP          cell_to_edge_op
#   define _GT_SKELETON_OP              skeleton_op
#   define _GT_BND_SKELETON_OP          bnd_skeleton_op
#   define _GT_CELL_UPDATE_OP           cell_update_op
#   define _GT_CELL_LAST_TOUCH_OP       cell_last_touch_op

#   include "SFC_generic_traversal_ringbuffer.f90"

    !*******************************
    !Geometry operators
    !*******************************

	subroutine pre_traversal_grid_op(traversal, grid)
			type(t_flash_euler_timestep_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							    :: grid

        grid%r_dt = cfg%courant_number * cfg%scaling * get_edge_size(grid%d_max) / ((2.0_GRID_SR + sqrt(2.0_GRID_SR)) * grid%u_max)
	!print *, "Calculated dt: ", grid%r_dt, cfg%scaling, grid%u_max,get_edge_size(grid%d_max)
	!print *, "calculatzed Dt: ",0.45_GRID_SR * grid%scaling * get_edge_size(grid%d_max) / ((2.0_GRID_SR + sqrt(2.0_GRID_SR)) * grid%u_max), grid%r_dt

#           if defined(_ASAGI)
                if (grid%r_time < grid_max_z(cfg%afh_bathymetry)) then
                    grid%r_dt = min(grid%r_dt, 0.1/15.0 * grid_max_z(cfg%afh_bathymetry))
	!	print *, "calculatzed Dt in if: ",min(grid%r_dt, 0.1/15.0 * grid_max_z(grid%afh_bathymetry)), grid%r_dt
                end if
#           endif

	grid%u_max = 0.0_GRID_SR
	end subroutine

	subroutine post_traversal_grid_op(traversal, grid)
			type(t_flash_euler_timestep_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							    :: grid

            call reduce(traversal%i_refinements_issued, traversal%children%i_refinements_issued, MPI_SUM, .true.)
            call reduce(grid%u_max, grid%sections%elements_alloc%u_max, MPI_MAX, .true.)
			grid%r_time = grid%r_time + grid%r_dt

			grid%sections%elements%r_time = grid%r_time
	end subroutine

    subroutine pre_traversal_op(traversal, section)
      type(t_FLASH_euler_timestep_traversal), intent(inout)     :: traversal
      type(t_grid_section), intent(inout)                       :: section

      !this variable will be incremented for each cell with a refinement request
      traversal%i_refinements_issued = 0
      section%u_max = 0.0_GRID_SR
    end subroutine

    subroutine post_traversal_op(traversal, section)
      type(t_FLASH_euler_timestep_traversal), intent(inout)     :: traversal
      type(t_grid_section), intent(inout)                       :: section
    end subroutine

    function cell_to_edge_op(element, edge) result(rep)
      type(t_element_base), intent(in)                          :: element
      type(t_edge_data), intent(in)                             :: edge
      type(num_cell_rep)                                        :: rep

      type(t_state), dimension(_FLASH_CELL_SIZE)                :: Q
      integer(kind = GRID_SI)                                   :: i, j, i_edge
      real(kind = GRID_SR), dimension(2, _FLASH_EDGE_SIZE)      :: dof_pos
      real(kind = GRID_SR), dimension(2, 3), parameter    :: edge_offsets = reshape([0.0, 0.0, 0.0, 1.0, 1.0, 0.0], [2, 3])
      real(kind = GRID_SR), dimension(2, 3), parameter    :: edge_vectors = reshape([0.0, 1.0, 1.0, -1.0, -1.0, 0.0], [2, 3])

      call gv_Q%read(element, Q)

      _log_write(6, '(3X, A)') "FLASH cell to edge op:"
      _log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q in: ", Q
            _log_write(6, '(4X, A, F0.3, 1X, F0.3)') "normal in : ", edge%transform_data%normal

!#     if (_FLASH_CELL_SIZE > 1)
      !  i_edge = edge%transform_data%index
  !      _log_write(6, '(4X, A, I0)') "edge ", i_edge
	!print *, "-------------------------------------------"
	!print *, "GLOBAL INDEXXXXXX: " , i_edge
!	print *, "LOCAL INDEXXXXXX: " ,  edge%transform_data%orientation
	!print *, "normal: ", edge%transform_data%normal
	!print *, "-------------------------------------------"
   !     forall (i = 1 : _FLASH_EDGE_SIZE)
    !      dof_pos(:, i) = edge_offsets(:, i_edge) + t_basis_flux_get_dof_coords(i) * edge_vectors(:, i_edge)
  !      end forall

       ! call lfs_flux%transform(edge%transform_data, dof_pos(1, :))
        !call lfs_flux%transform(edge%transform_data, dof_pos(2, :))

        !forall (i = 1 : _FLASH_EDGE_SIZE)
         ! rep%Q(i)%h = t_basis_Q_eval(dof_pos(:, i), Q%h)
          !rep%Q(i)%p(1) = t_basis_Q_eval(dof_pos(:, i), Q%p(1))
          !rep%Q(i)%p(2) = t_basis_Q_eval(dof_pos(:, i), Q%p(2))
        !end forall
!#     else
      !  rep%Q(1) = Q(1)
!#     endif
	rep%Q = Q





      _log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q out: ", rep%Q
    end function

    subroutine skeleton_array_op(traversal, grid, edges, rep1, rep2, update1, update2)
      type(t_FLASH_euler_timestep_traversal), intent(in)        :: traversal
      type(t_grid_section), intent(in)                          :: grid
      type(t_edge_data), intent(in)                             :: edges(:)
      type(num_cell_rep), intent(in)                            :: rep1(:), rep2(:)
      type(num_cell_update), intent(out)                        :: update1(:), update2(:)

      integer (kind = GRID_SI)                                  :: i

      do i = 1, size(edges)
        call skeleton_scalar_op(traversal, grid, edges(i), rep1(i), rep2(i), update1(i), update2(i))
      end do
    end subroutine

    subroutine skeleton_scalar_op(traversal, grid, edge, rep1, rep2, update1, update2)
      type(t_FLASH_euler_timestep_traversal), intent(in)        :: traversal
      type(t_grid_section), intent(in)                          :: grid
      type(t_edge_data), intent(in)                             :: edge
      type(num_cell_rep), intent(in)                            :: rep1, rep2
      type(num_cell_update), intent(out)                        :: update1, update2

      REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE,3)      :: r_rhs_l, r_rhs_r
      REAL (KIND = GRID_SR)					:: max_wave_speed
      REAL (KIND = GRID_SR)                                     :: r_minh_l, r_minh_r
      REAL (KIND = GRID_SR), DIMENSION(_FLASH_EDGE_SIZE)        :: r_h_l, r_hv_l, r_hu_l, &
                                                                   r_h_r, r_hv_r, r_hu_r

      r_minh_r = 1
      r_minh_l = 1
      r_h_l  = rep1%Q(:)%h
      r_hu_l = rep1%Q(:)%p(1)
      r_hv_l = rep1%Q(:)%p(2)
      r_h_r  = rep2%Q(:)%h
      r_hu_r = rep2%Q(:)%p(1)
      r_hv_r = rep2%Q(:)%p(2)

      !_log_write(0, *) "normal: ", edge%transform_data%normal
      _log_write(6, '(3X, A)') "FLASH skeleton op:"
      _log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q 1 in: ", rep1%Q
      _log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q 2 in: ", rep2%Q

      call compute_flash_flux(r_rhs_l, r_rhs_r,max_wave_speed, edge%transform_data%normal, r_minh_l, r_minh_r, &
                              _FLASH_CELL_SIZE, _FLASH_EDGE_SIZE, gquadwei, gMinvpsi, &
                              r_h_l, r_hu_l, r_hv_l, r_h_r, r_hu_r, r_hv_r,rep2%Q(1)%b)

      update1%flux(:)%h    = -r_rhs_l(:,1)
      update1%flux(:)%p(1) = -r_rhs_l(:,2)
      update1%flux(:)%p(2) = -r_rhs_l(:,3)

      update2%flux(:)%h    =  r_rhs_r(:,1)
      update2%flux(:)%p(1) =  r_rhs_r(:,2)
      update2%flux(:)%p(2) =  r_rhs_r(:,3)

	update1%flux(:)%max_wave_speed = max_wave_speed
	update2%flux(:)%max_wave_speed = max_wave_speed
      _log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "flux 1 out: ", update1%flux
      _log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "flux 2 out: ", update2%flux
    end subroutine

    subroutine bnd_skeleton_array_op(traversal, grid, edges, rep, update)
      type(t_FLASH_euler_timestep_traversal), intent(in)        :: traversal
      type(t_grid_section), intent(in)                          :: grid
      type(t_edge_data), intent(in)                             :: edges(:)
      type(num_cell_rep), intent(in)                            :: rep(:)
      type(num_cell_update), intent(out)                        :: update(:)

      integer (kind = GRID_SI)                                  :: i

      do i = 1, size(edges)
        call bnd_skeleton_scalar_op(traversal, grid, edges(i), rep(i), update(i))
      end do
    end subroutine

    subroutine bnd_skeleton_scalar_op(traversal, grid, edge, rep, update)
      type(t_FLASH_euler_timestep_traversal), intent(in)        :: traversal
      type(t_grid_section), intent(in)                          :: grid
      type(t_edge_data), intent(in)                             :: edge
      type(num_cell_rep), intent(in)                            :: rep
      type(num_cell_update), intent(out)                        :: update

      type(t_state)                                             :: bnd_rep
      type(t_update)                                            :: bnd_flux

      REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE,3)      :: r_rhs_l, r_rhs_r
      REAL (KIND = GRID_SR)					:: max_wave_speed
      REAL (KIND = GRID_SR)                                     :: r_minh_l, r_minh_r
      REAL (KIND = GRID_SR), DIMENSION(_FLASH_EDGE_SIZE)        :: r_h_l, r_hv_l, r_hu_l, &
                                                                   r_h_r, r_hv_r, r_hu_r

      bnd_rep = t_state(0.0, [0.0, 0.0],0.0, [0.0, 0.0], rep%Q(1)%b)
      r_minh_r = 1
      r_minh_l = 1
      r_h_l  = rep%Q(:)%h
      r_hu_l = rep%Q(:)%p(1)
      r_hv_l = rep%Q(:)%p(2)
      r_h_r  = rep%Q(:)%h
      r_hu_r = rep%Q(:)%p(1)
      r_hv_r = rep%Q(:)%p(2)

      !_log_write(0, *) "BOUNDARY normal: ", edge%transform_data%normal

      call compute_flash_flux(r_rhs_l, r_rhs_r, max_wave_speed,edge%transform_data%normal, r_minh_l, r_minh_r, &
                              _FLASH_CELL_SIZE, _FLASH_EDGE_SIZE, gquadwei, gMinvpsi, &
                              r_h_l, r_hu_l, r_hv_l, r_h_r, r_hu_r, r_hv_r,rep%Q(1)%b)

      update%flux(:)%h    = -r_rhs_l(:,1)
      update%flux(:)%p(1) = -r_rhs_l(:,2)
      update%flux(:)%p(2) = -r_rhs_l(:,3)

      update%flux(:)%max_wave_speed = max_wave_speed

    end subroutine

    subroutine cell_update_op(traversal, grid, element, update1, update2, update3)
      type(t_FLASH_euler_timestep_traversal), intent(inout)     :: traversal
      type(t_grid_section), intent(inout)                       :: grid
      type(t_element_base), intent(inout)                       :: element
      type(num_cell_update), intent(in)                         :: update1, update2, update3

      !local variables

      type(t_state), dimension(_FLASH_CELL_SIZE)                :: dQ

      call volume_op(traversal, grid, element, dQ, [update1%flux, update2%flux, update3%flux])

      call gv_Q%add(element, dQ)
    end subroutine

		subroutine cell_last_touch_op(traversal, grid, cell)
			type(t_FLASH_euler_timestep_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: grid
			type(t_cell_data_ptr), intent(inout)				:: cell
			integer (kind = 1)								:: depth
			real(kind = GRID_SR)							:: b_norm

			depth = cell%geometry%i_depth
			b_norm = minval(abs(cell%data_pers%Q%h - cell%data_pers%Q%b))

			!refine also on the coasts
			if (depth < cfg%i_max_depth .and. b_norm < 100.0_GRID_SR) then
				cell%geometry%refinement = 1
				traversal%i_refinements_issued = traversal%i_refinements_issued + 1_GRID_DI
			else if (b_norm < 300.0_GRID_SR) then
				cell%geometry%refinement = max(cell%geometry%refinement, 0)
			endif
		end subroutine

    !*******************************
    !Volume and DoF operators
    !*******************************

    subroutine volume_op(traversal, section, element, dQ, fluxes)
      type(t_FLASH_euler_timestep_traversal), intent(inout)     :: traversal
      type(t_grid_section), intent(inout)                       :: section
      type(t_element_base), intent(inout)                       :: element
      type(t_state), dimension(:), intent(out)                  :: dQ
      type(t_update), dimension(:), intent(in)                  :: fluxes

      real(kind = GRID_SR)                                      :: volume, dQ_norm, edge_lengths(3)
      integer (kind = 1)                                        :: i, depth

      _log_write(6, '(3X, A)') "FLASH cell update op:"
      _log_write(6, '(4X, A, 4(X, F0.3))') "edge 1 flux in:", fluxes(1)
      _log_write(6, '(4X, A, 4(X, F0.3))') "edge 2 flux in:", fluxes(2)
      _log_write(6, '(4X, A, 4(X, F0.3))') "edge 3 flux in:", fluxes(3)

      volume = element%cell%geometry%get_volume()
      edge_lengths = element%cell%geometry%get_edge_sizes()

      dQ%h    = sum(edge_lengths * fluxes%h)
      dQ%p(1) = sum(edge_lengths * fluxes%p(1))
      dQ%p(2) = sum(edge_lengths * fluxes%p(2))
      dQ%b    = 0.0_GRID_SR

      dQ%h_old    = sum(edge_lengths * fluxes%h)
      dQ%p_old(1) = sum(edge_lengths * fluxes%p(1))
      dQ%p_old(2) = sum(edge_lengths * fluxes%p(2))

      !set refinement condition

      element%cell%geometry%refinement = 0
      dQ_norm = dot_product(dQ(1)%p, dQ(1)%p)

      depth = element%cell%geometry%i_depth
      if (depth < cfg%i_max_depth .and. dQ_norm > (1.5_GRID_SR ** 2)) then
        element%cell%geometry%refinement = 1
        traversal%i_refinements_issued = traversal%i_refinements_issued + 1
      else if (depth > cfg%i_min_depth .and. dQ_norm < (1.45_GRID_SR ** 2)) then
        element%cell%geometry%refinement = -1
      endif

      section%u_max = max(section%u_max, maxval(fluxes%max_wave_speed))

      !_log_write(0, *) "U_max: ", section%u_max

        do i = 1, _FLASH_CELL_SIZE
            dQ(i)%t_dof_state = dQ(i)%t_dof_state * (-section%r_dt / volume)
        end do

      _log_write(6, '(4X, A, 4(X, F0.3))') "dQ out: ", dQ
    end subroutine

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine compute_flash_flux(r_rhs_l, r_rhs_r, max_wave_speed, r_normal, r_minh_l, r_minh_r, &
                      i_faceunknowns, i_gquadpts, r_gqwei, r_gMinvpsi, &
                      r_h_l, r_hu_l, r_hv_l, r_h_r, r_hu_r, r_hv_r,b)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(out):: r_rhs_l, r_rhs_r
    REAL (KIND = GRID_SR), INTENT(out):: max_wave_speed
    REAL (KIND = GRID_SR), INTENT(in)                 :: r_minh_l, r_minh_r
    REAL (KIND = GRID_SR), DIMENSION(2), INTENT(in)   :: r_normal
    INTEGER (KIND = GRID_SI), INTENT(in)              :: i_faceunknowns, i_gquadpts
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(in)   :: r_gqwei, &
                                                         r_h_l, r_hv_l, r_hu_l, &
                                                         r_h_r, r_hv_r, r_hu_r
    REAL (KIND = GRID_SR), DIMENSION(:,:), INTENT(in) :: r_gMinvpsi

!--- local declarations
    INTEGER (KIND = GRID_SI)                          :: i_quad, i_dof
    REAL (KIND = GRID_SR), DIMENSION(3)               :: r_intFstar
    REAL (KIND = GRID_SR), DIMENSION(3,2)             :: r_F_l, r_F_r
    REAL (KIND = GRID_SR), DIMENSION(3)               :: r_flux_l, r_flux_r, r_Fstar
    REAL (KIND = GRID_SR)			      :: b

    real(kind = GRID_SR)								:: vL

    r_rhs_l = 0._GRID_SR
    r_rhs_r = 0._GRID_SR
    r_intFstar(:) = 0._GRID_SR

!--- perform edge quadrature
    edge_quad_loop: DO i_quad=1,i_gquadpts

      r_F_l = flux(r_h_l(i_quad), r_hu_l(i_quad), r_hv_l(i_quad))
      r_F_r = flux(r_h_r(i_quad), r_hu_r(i_quad), r_hv_r(i_quad))

      r_flux_l = MATMUL(r_F_l, r_normal)
      r_flux_r = MATMUL(r_F_r, r_normal)

      r_Fstar  = riemannsolver(r_h_l(i_quad), r_h_r(i_quad), r_hu_l(i_quad), r_hu_r(i_quad), &
                               r_hv_l(i_quad), r_hv_r(i_quad), r_normal)

      ! note: this is not used right now but might be useful for wetting/drying...
      r_intFstar(:) = r_intFstar(:) + r_gqwei(i_quad)*r_Fstar

      edge_dof_loop: DO i_dof=1,i_faceunknowns
        r_rhs_l(i_dof,:) = r_rhs_l(i_dof,:) + &
                           r_gqwei(i_quad)*r_gMinvpsi(i_dof,i_quad)*(r_flux_l - r_Fstar)
        r_rhs_r(i_dof,:) = r_rhs_r(i_dof,:) + &
                           r_gqwei(i_quad)*r_gMinvpsi(i_dof,i_quad)*(r_flux_r - r_Fstar)


  !if (DOT_PRODUCT(r_rhs_l(1,:),r_rhs_l(1,:)) /= 0.0) then
  ! _log_write(0, *) "Fsatr: ", r_flux_l, r_flux_r, r_Fstar
  ! _log_write(0, *) "Fluxes: " , r_rhs_l, r_rhs_r
  ! _log_write(0, *) "h: ",r_h_l,r_h_r

  !endif

!-------------**************_____________****************_______________****************----------------------
!					YOU WERE HERE!
	max_wave_speed = 0
	vL = DOT_PRODUCT(r_normal, r_hu_l / (r_h_l - b))
	max_wave_speed = sqrt(g * (r_h_l(1) - b)) + sqrt(vL * vL)

      END DO edge_dof_loop
    END DO edge_quad_loop

    end subroutine
!*******************************************************************************
  FUNCTION flux(r_h, r_hu, r_hv) RESULT(r_flux)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), INTENT(in)     :: r_h, r_hu, r_hv
    REAL (KIND = GRID_SR), DIMENSION(3,2) :: r_flux

    REAL (KIND = GRID_SR)                 :: r_abs

    r_abs = heaviside(r_h)

    r_flux(1,1) = r_hu
    r_flux(1,2) = r_hv

    r_flux(2,1) = r_abs     * ((r_hu**2)+ 0.5*(r_h*r_h)) + &
                  (1-r_abs) * ((r_hu**2)/(r_h+10E-20) + 0.5*(r_h*r_h))
    r_flux(2,2) = r_abs     * (r_hu*r_hv) + &
                  (1-r_abs) * (r_hu*r_hv/(r_h+10E-20))

    r_flux(3,1) = r_abs     * (r_hu*r_hv) + &
                  (1-r_abs) * (r_hu*r_hv/(r_h+10E-20))
    r_flux(3,2) = r_abs     * ((r_hv**2)+ 0.5*(r_h*r_h)) + &
                  (1-r_abs) * ((r_hv**2)/(r_h+10E-20) + 0.5*(r_h*r_h))

  END FUNCTION flux

!*******************************************************************************

  FUNCTION riemannsolver(r_h_l, r_h_r, r_hu_l, r_hu_r, r_hv_l, r_hv_r, &
                         r_normal) RESULT(r_flux)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(3)               :: r_flux
    REAL (KIND = GRID_SR), INTENT(in)                 :: r_h_l, r_h_r, &
                                                         r_hu_l, r_hu_r, r_hv_l, r_hv_r
    REAL (KIND = GRID_SR), DIMENSION(:), INTENT(in) :: r_normal

    REAL (KIND = GRID_SR)                             :: r_u_l, r_u_r, r_v_l, r_v_r, &
                                                         r_a_l, r_a_r, r_u_norm_l, r_u_norm_r, &
                                                         r_abs_l, r_abs_r, r_lambda, r_tol
    REAL (KIND = GRID_SR), DIMENSION(3)               :: r_dq
    REAL (KIND = GRID_SR), DIMENSION(3,2)             :: r_flux_left, r_flux_right, &
                                                         r_flux_help, r_fluxrus

    r_tol = 10E-10

!--- compute velocity; if height reaches minimum, set velocity to zero
    IF(r_h_l < r_tol) THEN
      r_u_l = 0._GRID_SR
      r_v_l = 0._GRID_SR
    ELSE
      r_u_l = r_hu_l / r_h_l
      r_v_l = r_hv_l / r_h_l
    END IF
    IF(r_h_r < r_tol) THEN
      r_u_r = 0._GRID_SR
      r_v_r = 0._GRID_SR
    ELSE
      r_u_r = r_hu_r / r_h_r
      r_v_r = r_hv_r / r_h_r
    END IF

!--- compute normal component of velocity
    r_u_norm_l =  r_u_l*r_normal(1) + r_v_l*r_normal(2)
    r_u_norm_r =  r_u_r*r_normal(1) + r_v_r*r_normal(2)

!--- compute sound speeds
    r_a_l = sqrt(r_h_l)
    r_a_r = sqrt(r_h_r)

!--- compute maximal flux strength
    r_abs_l  = abs(r_u_norm_l) + r_a_l
    r_abs_r  = abs(r_u_norm_r) + r_a_r
    r_lambda = abs(max(r_abs_l, r_abs_r))

!--- compute left and right fluxes
    r_flux_left  = flux(r_h_l, r_hu_l, r_hv_l)
    r_flux_right = flux(r_h_r, r_hu_r, r_hv_r)

!--- compute q_r - q_l as three-dimensional vector
    r_dq(1) = r_h_r  - r_h_l
    r_dq(2) = r_hu_r - r_hu_l
    r_dq(3) = r_hv_r - r_hv_l

    r_flux_help = r_lambda*matmul(reshape(r_dq, (/ 3, 1 /)), reshape(r_normal, (/1,2/)))

    r_fluxrus = (r_flux_left + r_flux_right - r_flux_help)
    r_fluxrus = r_fluxrus/2

    r_flux = matmul(r_fluxrus, r_normal)

  END FUNCTION riemannsolver

!*******************************************************************************
  FUNCTION heaviside(r_real)

    IMPLICIT NONE

    REAL (KIND = GRID_SR)           :: r_real
    REAL (KIND = GRID_SR)           :: heaviside

!     heaviside = (1-SIGNUM(r_real))
    heaviside = MERGE(1._GRID_SR, 0._GRID_SR, r_real <= 0._GRID_SR)

  END FUNCTION heaviside

!*******************************************************************************


  END MODULE
#endif
