#include "Compilation_control.f90"

#if defined(_FLASH)
  MODULE FLASH_DG_Element

    USE SFC_data_types

    IMPLICIT NONE

    PUBLIC

    !> barycentric coordinates of element/face DOFs
    REAL (KIND = GRID_SR), DIMENSION(3, _FLASH_CELL_SIZE)                     :: r_edofcoo
    !> barycentric coordinates wrt edge of edge quadrature points
    REAL (KIND = GRID_SR), DIMENSION(2, _FLASH_EDGE_QUAD_SIZE)                :: r_gquadcoo
    !> weights for quadrature rule corresponding to edge quadrature points
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_EDGE_QUAD_SIZE)                   :: r_gqwei
    !> basis function evaluations at edge quadrature points
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_EDGE_SIZE, _FLASH_EDGE_QUAD_SIZE) :: r_gpsi
    !> nonzero basis function indices for reference edge
    INTEGER (KIND = GRID_SI), DIMENSION(_FLASH_EDGE_SIZE)                     :: i_gpsiidx
    !> barycentric coordinates of element/face quadrature points
    REAL (KIND = GRID_SR), DIMENSION(3, _FLASH_CELL_QUAD_SIZE)                :: r_equadcoo
    !> weights for the quadrature rule corresponding to the element quadrature points
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_QUAD_SIZE)                   :: r_eqwei
    !> basis function evaluations at element quadrature points
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE, _FLASH_CELL_QUAD_SIZE) :: r_epsi
    !> derivative matrix in \f$\xi\f$ direction
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE, _FLASH_CELL_SIZE)      :: r_Dxi
    !> derivative matrix in \f$\eta\f$ direction
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE, _FLASH_CELL_SIZE)      :: r_Deta
    !> Vandermonde matrix
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE, _FLASH_CELL_SIZE)      :: r_evander
    !> inverse Vandermonde matrix
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE, _FLASH_CELL_SIZE)      :: r_einvvander
    !> index to get dof order according to cell orientation
    INTEGER (KIND = GRID_SI), DIMENSION(_FLASH_CELL_SIZE, 2)                  :: i_reflect


    !< COMPUTED MATICES:
    !> inverse mass matrix
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE, _FLASH_CELL_SIZE)      :: r_emassinv
    !> transformed basis functions einvmass*psi at element quadrature points
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE, _FLASH_EDGE_QUAD_SIZE) :: r_gMinvpsi
    !> transformed basis functions einvmass*psi at edge quadrature points
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE, _FLASH_CELL_QUAD_SIZE) :: r_eMinvpsi
    !> transformed basis functions einvmass*psi at edge quadrature points
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE, _FLASH_CELL_QUAD_SIZE) :: r_eMinvdpsidxi
    !> transformed basis functions einvmass*psi at edge quadrature points
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE, _FLASH_CELL_QUAD_SIZE) :: r_eMinvdpsideta
    !> auxiliary data
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE, _FLASH_EDGE_QUAD_SIZE) :: r_gpsitmp
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE, _FLASH_CELL_QUAD_SIZE) :: r_edpsidxi, r_edpsideta

    CONTAINS

    SUBROUTINE dgelmt_init
      IMPLICIT NONE
# if (_FLASH_ORDER == 0)
      r_edofcoo(:,1)    = (/ 0.33333333333333333, 0.33333333333333333, 0.33333333333333333 /)

      r_gquadcoo(:,1)   = (/ 0.5, 0.5 /)

      r_gqwei(:)        = (/ 2.0 /)

      r_gpsi(:,1)       = (/ 1.0 /)

      i_gpsiidx(:)      = (/ 1 /)

      r_equadcoo(:,1)   = (/ 0.33333333333333333, 0.33333333333333333, 0.33333333333333333 /)

      r_eqwei(:)        = (/ 2.0 /)

      r_epsi(:,1)       = (/ 1.0 /)

      r_Dxi(:,1)        = (/ 0.0 /)

      r_Deta(:,1)       = (/ 0.0 /)

      r_evander(:,1)    = (/ 1.0 /)

      r_einvvander(:,1) = (/ 1.0 /)

      i_reflect(:,1)    = (/ 1 /)
      i_reflect(:,2)    = (/ 1 /)

# endif
# if (_FLASH_ORDER == 1)

      r_edofcoo(:,1)    = (/ 1.00000000000000e+00 , 0.00000000000000e+00 , 0.00000000000000e+00 /)
      r_edofcoo(:,2)    = (/ 0.00000000000000e+00 , 1.00000000000000e+00 , 0.00000000000000e+00 /)
      r_edofcoo(:,3)    = (/ 0.00000000000000e+00 , 0.00000000000000e+00 , 1.00000000000000e+00 /)

      r_gquadcoo(:,1)   = (/ 8.87298334620742e-01 , 1.12701665379258e-01 /)
      r_gquadcoo(:,2)   = (/ 5.00000000000000e-01 , 5.00000000000000e-01 /)
      r_gquadcoo(:,3)   = (/ 1.12701665379258e-01 , 8.87298334620742e-01 /)

      r_gqwei(:)        = (/ 5.55555555555556e-01 , 8.88888888888889e-01 , 5.55555555555555e-01 /)

      r_gpsi(:,1)       = (/ 8.87298334620742e-01 , 1.12701665379258e-01 /)
      r_gpsi(:,2)       = (/ 5.00000000000000e-01 , 5.00000000000000e-01 /)
      r_gpsi(:,3)       = (/ 1.12701665379258e-01 , 8.87298334620742e-01 /)

      i_gpsiidx(:)      = (/ 2 , 3 /)

      r_equadcoo(:,1)   = (/ 6.66666666666667e-01 , 1.66666666666667e-01 , 1.66666666666667e-01 /)
      r_equadcoo(:,2)   = (/ 1.66666666666667e-01 , 6.66666666666667e-01 , 1.66666666666667e-01 /)
      r_equadcoo(:,3)   = (/ 1.66666666666667e-01 , 1.66666666666667e-01 , 6.66666666666667e-01 /)

      r_eqwei(:)        = (/ 6.66666666666667e-01 , 6.66666666666667e-01 , 6.66666666666667e-01 /)

      r_epsi(:,1)       = (/ 6.66666666666667e-01 , 1.66666666666666e-01 , 1.66666666666667e-01 /)
      r_epsi(:,2)       = (/ 1.66666666666667e-01 , 6.66666666666667e-01 , 1.66666666666667e-01 /)
      r_epsi(:,3)       = (/ 1.66666666666667e-01 , 1.66666666666667e-01 , 6.66666666666667e-01 /)

      r_Dxi(:,1)        = (/-5.00000000000000e-01 , 5.00000000000000e-01 , 0.00000000000000e+00 /)
      r_Dxi(:,2)        = (/-5.00000000000000e-01 , 5.00000000000000e-01 , 0.00000000000000e+00 /)
      r_Dxi(:,3)        = (/-5.00000000000000e-01 , 5.00000000000000e-01 , 0.00000000000000e+00 /)

      r_Deta(:,1)       = (/-5.00000000000000e-01 , 0.00000000000000e+00 , 5.00000000000000e-01 /)
      r_Deta(:,2)       = (/-5.00000000000000e-01 , 0.00000000000000e+00 , 5.00000000000000e-01 /)
      r_Deta(:,3)       = (/-5.00000000000000e-01 , 0.00000000000000e+00 , 5.00000000000000e-01 /)

      r_evander(:,1)    = (/ 7.07106781186547e-01 ,-1.00000000000000e+00 ,-1.73205080756888e+00 /)
      r_evander(:,2)    = (/ 7.07106781186547e-01 ,-1.00000000000000e+00 , 1.73205080756888e+00 /)
      r_evander(:,3)    = (/ 7.07106781186547e-01 , 2.00000000000000e+00 , 0.00000000000000e+00 /)

      r_einvvander(:,1) = (/ 4.71404520791032e-01 , 4.71404520791032e-01 , 4.71404520791032e-01 /)
      r_einvvander(:,2) = (/-1.66666666666667e-01 ,-1.66666666666667e-01 , 3.33333333333333e-01 /)
      r_einvvander(:,3) = (/-2.88675134594813e-01 , 2.88675134594813e-01 , 0.00000000000000e+00 /)

      i_reflect(:,1)    = (/ 1 , 2 , 3 /)
      i_reflect(:,2)    = (/ 1 , 3 , 2 /)

# endif

      ! some corrections since we use different reference elements in Stormflash2d and Samoa
      r_gqwei     = r_gqwei / 2.0_GRID_SR
!       r_eqwei     = r_eqwei / 2.0_GRID_SR
!       r_Dxi       = r_Dxi   * 2.0_GRID_SR
!       r_Deta      = r_Deta  * 2.0_GRID_SR
      ! Note Stefan: maybe r_evander and r_einvvander have also to adjusted

      ! compute transformed basis functions
      r_edpsidxi  = MATMUL(r_Dxi , r_epsi)
      r_edpsideta = MATMUL(r_Deta, r_epsi)
      r_emassinv  = MATMUL(r_evander, TRANSPOSE(r_evander))

      r_gpsitmp(:,:)         = 0.0_GRID_SR
      r_gpsitmp(i_gpsiidx,:) = r_gpsi

      r_gMinvpsi      = MATMUL(r_emassinv, r_gpsitmp)
      r_eMinvpsi      = MATMUL(r_emassinv, r_epsi)
      r_eMinvdpsidxi  = MATMUL(r_emassinv, r_edpsidxi)
      r_eMinvdpsideta = MATMUL(r_emassinv, r_edpsideta)

    END SUBROUTINE dgelmt_init

  END MODULE FLASH_DG_Element
#endif
