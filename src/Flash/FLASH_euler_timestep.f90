! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema, Stefan Vater
! This program is licensed under the GPL, for details see the file LICENSE

#include "Compilation_control.f90"


#if defined(_FLASH)
  MODULE FLASH_Euler_Timestep
    use SFC_edge_traversal

    use Samoa_FLASH

    use FLASH_dg_element
    USE DG_riemann_solver
    USE DG_equation

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

!     grid%r_dt = 0.45_GRID_SR * cfg%scaling * get_edge_size(grid%d_max) / ((2.0_GRID_SR + sqrt(2.0_GRID_SR)) * grid%u_max)
    grid%r_dt = 1e-3!0.45_GRID_SR * get_edge_size(grid%d_max) / ((2.0_GRID_SR + sqrt(2.0_GRID_SR)) * grid%u_max)
!     print *, "Calculated dt: ", grid%r_dt, cfg%scaling, grid%u_max,get_edge_size(grid%d_max)
!     print *, "calculatzed Dt: ",0.45_GRID_SR * grid%scaling * get_edge_size(grid%d_max) / ((2.0_GRID_SR + sqrt(2.0_GRID_SR)) * grid%u_max), grid%r_dt

# if defined(_ASAGI)
    if (grid%r_time < grid_max_z(cfg%afh_bathymetry)) then
      grid%r_dt = min(grid%r_dt, 0.1/15.0 * grid_max_z(cfg%afh_bathymetry))
!       print *, "calculatzed Dt in if: ",min(grid%r_dt, 0.1/15.0 * grid_max_z(grid%afh_bathymetry)), grid%r_dt
    end if
# endif

    grid%u_max = 0.0_GRID_SR
  end subroutine

!*******************************************************************************
  subroutine post_traversal_grid_op(traversal, grid)
    type(t_flash_euler_timestep_traversal), intent(inout)		:: traversal
    type(t_grid), intent(inout)							    :: grid

    call reduce(traversal%i_refinements_issued, traversal%children%i_refinements_issued, MPI_SUM, .true.)
    call reduce(grid%u_max, grid%sections%elements_alloc%u_max, MPI_MAX, .true.)
    grid%r_time = grid%r_time + grid%r_dt

    grid%sections%elements%r_time = grid%r_time

  end subroutine

!*******************************************************************************
  subroutine pre_traversal_op(traversal, section)
    type(t_FLASH_euler_timestep_traversal), intent(inout)     :: traversal
    type(t_grid_section), intent(inout)                       :: section

    !this variable will be incremented for each cell with a refinement request
    traversal%i_refinements_issued = 0
    section%u_max = 0.0_GRID_SR
  end subroutine

!*******************************************************************************
  subroutine post_traversal_op(traversal, section)
    type(t_FLASH_euler_timestep_traversal), intent(inout)     :: traversal
    type(t_grid_section), intent(inout)                       :: section
  end subroutine

!*******************************************************************************
  function cell_to_edge_op(element, edge) result(rep)
    type(t_element_base), intent(in)                          :: element
    type(t_edge_data), intent(in)                             :: edge
    type(num_cell_rep)                                        :: rep

    type(t_state), dimension(_FLASH_CELL_SIZE)                :: Q
    integer(kind = GRID_SI)                                   :: i, i_dof, i_edge
    real(kind = GRID_SR), dimension(2, _FLASH_EDGE_SIZE)      :: dof_pos
    real(kind = GRID_SR), dimension(2, 3), parameter    :: edge_offsets = reshape([0.0, 0.0, 0.0, 1.0, 1.0, 0.0], [2, 3])
    real(kind = GRID_SR), dimension(2, 3), parameter    :: edge_vectors = reshape([0.0, 1.0, 1.0, -1.0, -1.0, 0.0], [2, 3])

    call gv_Q%read(element, Q)

    _log_write(6, '(3X, A)') "FLASH cell to edge op:"
    _log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q in: ", Q
    _log_write(6, '(4X, A, F0.3, 1X, F0.3)') "normal in : ", edge%transform_data%normal

!# if (_FLASH_CELL_SIZE > 1)
!     i_edge = edge%transform_data%index
!     _log_write(6, '(4X, A, I0)') "edge ", i_edge
!     print *, "-------------------------------------------"
!     print *, "GLOBAL INDEXXXXXX: " , i_edge
!     print *, "LOCAL INDEXXXXXX: " ,  edge%transform_data%orientation
!     print *, "normal: ", edge%transform_data%normal
!     print *, "-------------------------------------------"
!     forall (i = 1 : _FLASH_EDGE_SIZE)
!       dof_pos(:, i) = edge_offsets(:, i_edge) + t_basis_flux_get_dof_coords(i) * edge_vectors(:, i_edge)
!     end forall
!
!     call lfs_flux%transform(edge%transform_data, dof_pos(1, :))
!     call lfs_flux%transform(edge%transform_data, dof_pos(2, :))
!
!     forall (i = 1 : _FLASH_EDGE_SIZE)
!       rep%Q(i)%h = t_basis_Q_eval(dof_pos(:, i), Q%h)
!       rep%Q(i)%p(1) = t_basis_Q_eval(dof_pos(:, i), Q%p(1))
!       rep%Q(i)%p(2) = t_basis_Q_eval(dof_pos(:, i), Q%p(2))
!     end forall
!# else
!     rep%Q(1) = Q(1)
!# endif
    forall (i = 1 : _FLASH_EDGE_SIZE)
      rep%Q(i) = Q(i)
    end forall
    ! set some markers which are later used for wetting and drying
    rep%iswete = 1.0_GRID_SR
    rep%iswetg = 1.0_GRID_SR

    _log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q out: ", rep%Q
  end function

!*******************************************************************************
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

!*******************************************************************************
  SUBROUTINE skeleton_scalar_op(traversal, grid, edge, rep1, rep2, update1, update2)

    IMPLICIT NONE

    type(t_FLASH_euler_timestep_traversal), INTENT(in)      :: traversal
    type(t_grid_section), INTENT(in)                        :: grid
    type(t_edge_data), INTENT(in)                           :: edge
    type(num_cell_rep), INTENT(in)                          :: rep1, rep2
    type(num_cell_update), intent(out)                      :: update1, update2

!--- local declarations
    INTEGER (KIND = GRID_SI)                                :: i_quad, i_dof
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_EDGE_SIZE)      :: r_h, r_hu, r_hv
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_EDGE_QUAD_SIZE) :: r_h_l, r_hu_l, r_hv_l, &
                                                               r_h_r, r_hu_r, r_hv_r
    REAL (KIND = GRID_SR), DIMENSION(3)                     :: r_Fstar, r_Fleft, r_Frght

    forall (i_dof = 1 : _FLASH_CELL_SIZE)
      update1%flux(i_dof)%h    = 0.0_GRID_SR
      update1%flux(i_dof)%p(1) = 0.0_GRID_SR
      update1%flux(i_dof)%p(2) = 0.0_GRID_SR
      update2%flux(i_dof)%h    = 0.0_GRID_SR
      update2%flux(i_dof)%p(1) = 0.0_GRID_SR
      update2%flux(i_dof)%p(2) = 0.0_GRID_SR
    end forall

    IF (MAXVAL(rep1%Q(:)%h) > r_wettol .OR. MAXVAL(rep2%Q(:)%h) > r_wettol) THEN

!-- compute values at quadrature points
      forall (i_dof = 1 : _FLASH_EDGE_SIZE)
        r_h(i_dof)  = rep1%Q(i_dof)%h
        r_hu(i_dof) = rep1%Q(i_dof)%p(1)
        r_hv(i_dof) = rep1%Q(i_dof)%p(2)
      end forall

      r_h_l  = MATMUL(r_h , r_gpsi)
      r_hu_l = MATMUL(r_hu, r_gpsi)
      r_hv_l = MATMUL(r_hv, r_gpsi)

      forall (i_dof = 1 : _FLASH_EDGE_SIZE)
        r_h(i_dof)  = rep2%Q(i_dof)%h
        r_hu(i_dof) = rep2%Q(i_dof)%p(1)
        r_hv(i_dof) = rep2%Q(i_dof)%p(2)
      end forall

      r_h_r  = MATMUL(r_h , r_gpsi)
      r_hu_r = MATMUL(r_hu, r_gpsi)
      r_hv_r = MATMUL(r_hv, r_gpsi)

!--- perform edge quadrature
      edge_quad_loop: DO i_quad=1,_FLASH_EDGE_QUAD_SIZE

        r_Fstar = riemannsolver(r_h_l(i_quad), r_h_r(i_quad), r_hu_l(i_quad), r_hu_r(i_quad), &
                                r_hv_l(i_quad), r_hv_r(i_quad), edge%transform_data%normal)

        r_Fleft = (1.0-rep1%iswete)*rep1%iswetg*MATMUL(fluxGrav(r_h_l(i_quad)), edge%transform_data%normal)
        r_Frght = (1.0-rep2%iswete)*rep2%iswetg*MATMUL(fluxGrav(r_h_r(i_quad)), edge%transform_data%normal)

        edge_dof_loop: DO i_dof=1,_FLASH_CELL_SIZE
          update1%flux(i_dof)%h = update1%flux(i_dof)%h + r_gqwei(i_quad)*r_gMinvpsi(i_dof,i_quad)*(r_Fstar(1)  -r_Fleft(1)  )
          update1%flux(i_dof)%p = update1%flux(i_dof)%p + r_gqwei(i_quad)*r_gMinvpsi(i_dof,i_quad)*(r_Fstar(2:3)-r_Fleft(2:3))
          update2%flux(i_dof)%h = update2%flux(i_dof)%h - r_gqwei(i_quad)*r_gMinvpsi(i_dof,i_quad)*(r_Fstar(1)  -r_Frght(1)  )
          update2%flux(i_dof)%p = update2%flux(i_dof)%p - r_gqwei(i_quad)*r_gMinvpsi(i_dof,i_quad)*(r_Fstar(2:3)-r_Frght(2:3))
        END DO edge_dof_loop

      END DO edge_quad_loop

    END IF

  END SUBROUTINE skeleton_scalar_op

!*******************************************************************************
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

!*******************************************************************************
  SUBROUTINE bnd_skeleton_scalar_op(traversal, grid, edge, rep, update)

    IMPLICIT NONE

    type(t_FLASH_euler_timestep_traversal), INTENT(in)      :: traversal
    type(t_grid_section), INTENT(in)                        :: grid
    type(t_edge_data), INTENT(in)                           :: edge
    type(num_cell_rep), INTENT(in)                          :: rep
    type(num_cell_update), intent(out)                      :: update

!--- local declarations
    INTEGER (KIND = GRID_SI)                                :: i_quad, i_dof
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_EDGE_SIZE)      :: r_h, r_hu, r_hv
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_EDGE_QUAD_SIZE) :: r_h_l, r_hu_l, r_hv_l
    REAL (KIND = GRID_SR), DIMENSION(3)                     :: r_Fstar, r_Fleft

    forall (i_dof = 1 : _FLASH_CELL_SIZE)
      update%flux(i_dof)%h    = 0.0_GRID_SR
      update%flux(i_dof)%p(1) = 0.0_GRID_SR
      update%flux(i_dof)%p(2) = 0.0_GRID_SR
    end forall

    IF (MAXVAL(rep%Q(:)%h) > r_wettol) THEN

!-- compute values at quadrature points
      forall (i_dof = 1 : _FLASH_EDGE_SIZE)
        r_h(i_dof)  = rep%Q(i_dof)%h
        r_hu(i_dof) = rep%Q(i_dof)%p(1)
        r_hv(i_dof) = rep%Q(i_dof)%p(2)
      end forall

      r_h_l  = MATMUL(r_h , r_gpsi)
      r_hu_l = MATMUL(r_hu, r_gpsi)
      r_hv_l = MATMUL(r_hv, r_gpsi)

!--- perform edge quadrature
      edge_quad_loop: DO i_quad=1,_FLASH_EDGE_QUAD_SIZE

        r_Fstar = riemannsolver(r_h_l(i_quad), r_h_l(i_quad), r_hu_l(i_quad), r_hu_l(i_quad), &
                                r_hv_l(i_quad), r_hv_l(i_quad), edge%transform_data%normal)

        r_Fleft = (1.0-rep%iswete)*rep%iswetg*MATMUL(fluxGrav(r_h_l(i_quad)), edge%transform_data%normal)

        edge_dof_loop: DO i_dof=1,_FLASH_CELL_SIZE
          update%flux(i_dof)%h = update%flux(i_dof)%h + r_gqwei(i_quad)*r_gMinvpsi(i_dof,i_quad)*(r_Fstar(1)  -r_Fleft(1)  )
          update%flux(i_dof)%p = update%flux(i_dof)%p + r_gqwei(i_quad)*r_gMinvpsi(i_dof,i_quad)*(r_Fstar(2:3)-r_Fleft(2:3))
        END DO edge_dof_loop

      END DO edge_quad_loop

    END IF

  END SUBROUTINE bnd_skeleton_scalar_op

!*******************************************************************************
  subroutine cell_update_op(traversal, section, element, update1, update2, update3)
    type(t_FLASH_euler_timestep_traversal), intent(inout)     :: traversal
    type(t_grid_section), intent(inout)                       :: section
    type(t_element_base), intent(inout)                       :: element
    type(num_cell_update), intent(in)                         :: update1, update2, update3

    !local variables

    type(t_update), dimension(3)                              :: fluxes
    type(t_update)                                            :: eflux
    type(t_state), dimension(_FLASH_CELL_SIZE)                :: dQ
    real(kind = GRID_SR)                                      :: volume, dQ_norm, edge_lengths(3)
    integer (kind = 1)                                        :: depth, i_dof

    fluxes = [update1%flux, update2%flux, update3%flux]
!
!       _log_write(6, '(3X, A)') "FLASH cell update op:"
!       _log_write(6, '(4X, A, 4(X, F0.3))') "edge 1 flux in:", fluxes(1)
!       _log_write(6, '(4X, A, 4(X, F0.3))') "edge 2 flux in:", fluxes(2)
!       _log_write(6, '(4X, A, 4(X, F0.3))') "edge 3 flux in:", fluxes(3)

    volume       = element%cell%geometry%get_volume()
    edge_lengths = element%cell%geometry%get_edge_sizes()

    dQ%h    = sum(edge_lengths * fluxes%h)
    dQ%p(1) = sum(edge_lengths * fluxes%p(1))
    dQ%p(2) = sum(edge_lengths * fluxes%p(2))
    dQ%b    = 0.0_GRID_SR

    dQ%h_old    = sum(edge_lengths * fluxes%h)
    dQ%p_old(1) = sum(edge_lengths * fluxes%p(1))
    dQ%p_old(2) = sum(edge_lengths * fluxes%p(2))

    ! set refinement condition

!     element%cell%geometry%refinement = 0
!     dQ_norm = dot_product(dQ(1)%p, dQ(1)%p)
!
!     depth = element%cell%geometry%i_depth
!     if (depth < cfg%i_max_depth .and. dQ_norm > (1.5_GRID_SR ** 2)) then
!       element%cell%geometry%refinement = 1
!       traversal%i_refinements_issued = traversal%i_refinements_issued + 1
!     else if (depth > cfg%i_min_depth .and. dQ_norm < (1.45_GRID_SR ** 2)) then
!       element%cell%geometry%refinement = -1
!     endif

!     section%u_max = max(section%u_max, maxval(fluxes%max_wave_speed))

    !_log_write(0, *) "U_max: ", section%u_max

    forall (i_dof = 1 : _FLASH_CELL_SIZE)
      dQ(i_dof)%t_dof_state = dQ(i_dof)%t_dof_state * (-section%r_dt / volume)
    end forall

!     _log_write(6, '(4X, A, 4(X, F0.3))') "dQ out: ", dQ

!       call volume_op(section, element, eflux)

    call gv_Q%add(element, dQ)

  end subroutine

!*******************************************************************************
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

!*******************************************************************************
  subroutine volume_op(section, element, update)
    type(t_grid_section), intent(inout)                       :: section
    type(t_element_base), intent(inout)                       :: element
    type(num_cell_update), intent(inout)                        :: update

    real(kind = GRID_SR)                                      :: volume

    _log_write(6, '(3X, A)') "FLASH volume op:"

    volume = element%cell%geometry%get_volume()

  end subroutine

  END MODULE

#endif
