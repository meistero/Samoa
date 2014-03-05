#include "Compilation_control.f90"

#if defined(_FLASH)
  MODULE FLASH_DG_Element

    USE SFC_data_types

    IMPLICIT NONE

    PUBLIC

    !> number of basis functions nonzero on a specific edge
# if (_FLASH_ORDER == 0)
#define _GPSINONZERO  1
# else if (_FLASH_ORDER == 1)
#define _GPSINONZERO  2
# endif

    !> barycentric coordinates of element/face DOFs
    REAL (KIND = GRID_SR), DIMENSION(3, _FLASH_CELL_SIZE)                   :: edofcoo
    !> barycentric coordinates wrt edge of edge quadrature points
    REAL (KIND = GRID_SR), DIMENSION(2, _FLASH_EDGE_SIZE)                   :: gquadcoo
    !> weights for quadrature rule corresponding to edge quadrature points
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_EDGE_SIZE)                      :: gquadwei
    !> basis function evaluations at edge quadrature points
    REAL (KIND = GRID_SR), DIMENSION(_GPSINONZERO, _FLASH_EDGE_SIZE)        :: gpsiquad
    !> nonzero basis function indices for reference edge
    INTEGER (KIND = GRID_SI), DIMENSION(_GPSINONZERO)                       :: gpsiidx
    !> barycentric coordinates of element/face quadrature points
    REAL (KIND = GRID_SR), DIMENSION(3, _FLASH_CELL_SIZE)                   :: equadcoo
    !> weights for the quadrature rule corresponding to the element quadrature points
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE)                      :: equadwei
    !> basis function evaluations at element quadrature points
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE, _FLASH_CELL_SIZE)    :: epsiquad
    !> derivative matrix in \f$\xi\f$ direction
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE, _FLASH_CELL_SIZE)    :: dpsidxi
    !> derivative matrix in \f$\eta\f$ direction
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE, _FLASH_CELL_SIZE)    :: dpsideta
    !> Vandermonde matrix
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE, _FLASH_CELL_SIZE)    :: evander
    !> inverse Vandermonde matrix
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE, _FLASH_CELL_SIZE)    :: einvvander

    !< COMPUTED MATICES:
    !> inverse mass matrix
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE, _FLASH_CELL_SIZE)    :: emassinv
    !> transformed basis functions einvmass*psi at element quadrature points
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE, _FLASH_EDGE_SIZE)    :: gMinvpsi
    !> transformed basis functions einvmass*psi at edge quadrature points
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE, _FLASH_CELL_SIZE)    :: eMinvpsi
    !> auxiliary data
    REAL (KIND = GRID_SR), DIMENSION(_FLASH_CELL_SIZE, _FLASH_EDGE_SIZE)    :: gpsitmp

    CONTAINS

    SUBROUTINE dgelmt_init
      IMPLICIT NONE
# if (_FLASH_ORDER == 0)
      edofcoo(:,1)    = (/ 0.33333333333333333, 0.33333333333333333, 0.33333333333333333 /)

      gquadcoo(:,1)   = (/ 0.5, 0.5 /)

      gquadwei(:)     = (/ 1.0 /) !(/ 2.0 /)

      gpsiquad(:,1)   = (/ 1.0 /)

      gpsiidx(:)      = (/ 1 /)

      equadcoo(:,1)   = (/ 0.33333333333333333, 0.33333333333333333, 0.33333333333333333 /)

      equadwei(:)     = (/ 1.0 /) !(/ 2.0 /)

      epsiquad(:,1)   = (/ 1.0 /)

      dpsidxi(:,1)    = (/ 0.0 /)

      dpsideta(:,1)   = (/ 0.0 /)

      evander(:,1)    = (/ 1.0 /)

      einvvander(:,1) = (/ 1.0 /)

# else if (_FLASH_ORDER == 1)

    edofcoo(:,1)    = (/ 1.00000000000000e+00 , 0.00000000000000e+00 , 0.00000000000000e+00 /)
    edofcoo(:,2)    = (/ 0.00000000000000e+00 , 1.00000000000000e+00 , 0.00000000000000e+00 /)
    edofcoo(:,3)    = (/ 0.00000000000000e+00 , 0.00000000000000e+00 , 1.00000000000000e+00 /)

    gquadcoo(:,1)   = (/ 8.87298334620742e-01 , 1.12701665379258e-01 /)
    gquadcoo(:,2)   = (/ 5.00000000000000e-01 , 5.00000000000000e-01 /)
    gquadcoo(:,3)   = (/ 1.12701665379258e-01 , 8.87298334620742e-01 /)

    gquadwei(:)     = (/ 5.55555555555556e-01 , 8.88888888888889e-01 , 5.55555555555555e-01 /)

    gpsiquad(:,1)   = (/ 8.87298334620742e-01 , 1.12701665379258e-01 /)
    gpsiquad(:,2)   = (/ 5.00000000000000e-01 , 5.00000000000000e-01 /)
    gpsiquad(:,3)   = (/ 1.12701665379258e-01 , 8.87298334620742e-01 /)

    gpsiidx(:)      = (/ 2 , 3 /)

    equadcoo(:,1)   = (/ 6.66666666666667e-01 , 1.66666666666667e-01 , 1.66666666666667e-01 /)
    equadcoo(:,2)   = (/ 1.66666666666667e-01 , 6.66666666666667e-01 , 1.66666666666667e-01 /)
    equadcoo(:,3)   = (/ 1.66666666666667e-01 , 1.66666666666667e-01 , 6.66666666666667e-01 /)

    equadwei(:)     = (/ 6.66666666666667e-01 , 6.66666666666667e-01 , 6.66666666666667e-01 /)

    epsiquad(:,1)   = (/ 6.66666666666667e-01 , 1.66666666666666e-01 , 1.66666666666667e-01 /)
    epsiquad(:,2)   = (/ 1.66666666666667e-01 , 6.66666666666667e-01 , 1.66666666666667e-01 /)
    epsiquad(:,3)   = (/ 1.66666666666667e-01 , 1.66666666666667e-01 , 6.66666666666667e-01 /)

    dpsidxi(:,1)    = (/-5.00000000000000e-01 , 5.00000000000000e-01 , 0.00000000000000e+00 /)
    dpsidxi(:,2)    = (/-5.00000000000000e-01 , 5.00000000000000e-01 , 0.00000000000000e+00 /)
    dpsidxi(:,3)    = (/-5.00000000000000e-01 , 5.00000000000000e-01 , 0.00000000000000e+00 /)

    dpsideta(:,1)   = (/-5.00000000000000e-01 , 0.00000000000000e+00 , 5.00000000000000e-01 /)
    dpsideta(:,2)   = (/-5.00000000000000e-01 , 0.00000000000000e+00 , 5.00000000000000e-01 /)
    dpsideta(:,3)   = (/-5.00000000000000e-01 , 0.00000000000000e+00 , 5.00000000000000e-01 /)

    evander(:,1)    = (/ 7.07106781186547e-01 ,-1.00000000000000e+00 ,-1.73205080756888e+00 /)
    evander(:,2)    = (/ 7.07106781186547e-01 ,-1.00000000000000e+00 , 1.73205080756888e+00 /)
    evander(:,3)    = (/ 7.07106781186547e-01 , 2.00000000000000e+00 , 0.00000000000000e+00 /)

    einvvander(:,1) = (/ 4.71404520791032e-01 , 4.71404520791032e-01 , 4.71404520791032e-01 /)
    einvvander(:,2) = (/-1.66666666666667e-01 ,-1.66666666666667e-01 , 3.33333333333333e-01 /)
    einvvander(:,3) = (/-2.88675134594813e-01 , 2.88675134594813e-01 , 0.00000000000000e+00 /)

# endif

      emassinv = MATMUL(evander, TRANSPOSE(evander))

      gpsitmp(:,:)       = 0.0_GRID_SR
      gpsitmp(gpsiidx,:) = gpsiquad

      gMinvpsi = MATMUL(emassinv, gpsitmp)
      eMinvpsi = MATMUL(emassinv, epsiquad)

    END SUBROUTINE dgelmt_init

  END MODULE FLASH_DG_Element
#endif
