!*******************************************************************************
! MODULE NAME:
!       DG_equation
! FUNCTION:
!       Here the user has to specify the flux and the source term.
!       The divergence of the flux also has to be specified manualy.
!       Also for dual, if desired.
!
! VERSION(s):
! 1. first version            s.beckers        07/2013
! 2. include dummies for DWR  s.beckers        06/2014
!*******************************************************************************
#include "Compilation_control.f90"


#if defined(_FLASH)

MODULE DG_equation

  USE SFC_data_types
!   USE GRID_api

  PUBLIC :: r_wettol, velocity, divflux, flux, fluxGrav, source

  REAL (KIND = GRID_SR)                             :: r_wettol = 1.0E-8

  CONTAINS

!*******************************************************************************
  PURE FUNCTION velocity(r_h, r_hu) RESULT(r_u)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), INTENT(in)               :: r_h, r_hu
    REAL (KIND = GRID_SR)                           :: r_u

    r_u = MERGE(r_hu/r_h, 0.0_GRID_SR, r_h > r_wettol)

  END FUNCTION velocity

!*******************************************************************************
!> @brief Routine for flux computation. Has to be defined by the user.
!>
!*******************************************************************************
  FUNCTION flux(r_h, r_hu, r_hv, r_bathy, r_iswet) RESULT(r_flux)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), INTENT(in)               :: r_h, r_hu, r_hv
    REAL (KIND = GRID_SR), INTENT(in), OPTIONAL     :: r_bathy
    REAL (KIND = GRID_SR), INTENT(in), OPTIONAL     :: r_iswet
    REAL (KIND = GRID_SR), DIMENSION(3,2)           :: r_flux

    REAL (KIND = GRID_SR)                           :: r_u, r_v

    r_u = velocity(r_h, r_hu)
    r_v = velocity(r_h, r_hv)

    r_flux(1,1) = r_hu
    r_flux(1,2) = r_hv

    r_flux(2,1) = r_hu*r_u
    r_flux(2,2) = r_hu*r_v

    r_flux(3,1) = r_hv*r_u
    r_flux(3,2) = r_hv*r_v

    IF (.NOT. PRESENT(r_iswet) .OR. r_iswet==1.0) THEN
      r_flux(2,1) = r_flux(2,1) + 0.5*r_h**2
      r_flux(3,2) = r_flux(3,2) + 0.5*r_h**2
    END IF

  END FUNCTION flux

!*******************************************************************************
  FUNCTION fluxGrav(r_h) RESULT(r_flux)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), INTENT(in)               :: r_h
    REAL (KIND = GRID_SR), DIMENSION(3,2)           :: r_flux

    r_flux = 0.0_GRID_SR

    r_flux(2,1) = 0.5*r_h**2
    r_flux(3,2) = 0.5*r_h**2

  END FUNCTION fluxGrav

!*******************************************************************************
  FUNCTION source(r_h, r_hu, r_hv, r_wet, r_gradb) RESULT(r_source)

    IMPLICIT NONE

    REAL (KIND = GRID_SR)                  :: r_h, r_hu, r_hv
    REAL (KIND = GRID_SR), DIMENSION(:)    :: r_gradb
    REAL (KIND = GRID_SR)                  :: r_tol, r_wet

    REAL (KIND = GRID_SR), DIMENSION(3)    :: r_source

!------- compute the source terms
    r_source(1) = 0.0_GRID_SR
    r_source(2) = r_wet*r_h*r_gradb(1)
    r_source(3) = r_wet*r_h*r_gradb(2)

  END FUNCTION

!*******************************************************************************
! DESCRIPTION of [SUBROUTINE divflux]:
!> @brief computes the divergence and source terms of the equation in consistent form
!>
!> @param[inout]        r_rhs           REAL array of DIMENSION(3);
!>                                      sum of divergence and source term per quadrature point
!> @param[in]           i_face          INTEGER; local number of quadrature point
!> @param[in]           r_h_e           REAL array of DIMENSION(i_faceunknowns);
!>                                      height at element DOFs
!> @param[in]           r_bathy_e       REAL array of DIMENSION(i_faceunknowns);
!>                                      bathymetry at element DOFs
!> @param[in]           r_hu_e          REAL array of DIMENSION(i_faceunknowns);
!>                                      momentum in x direction at element DOFs
!> @param[in]           r_hv_e          REAL array of DIMENSION(i_faceunknowns);
!>                                      momentum in y direction at element DOFs
!> @param[in]           r_Jacobi        REAL array of DIMENSION(2,2);
!>                                      metric terms for the element
!> @param[in]           i_faceunknowns  INTEGER; number of DOFs inside the element
!> @param[in]           r_Dxi           REAL array of DIMENSION(i_faceunknowns, i_faceunknowns);
!>                                      differentiated basisfunctions in x direction at Lagrange points
!> @param[in]           r_Deta          REAL array of DIMENSION(i_faceunknowns, i_faceunknowns);
!>                                      differentiated basisfunctions in y direction at Lagrange points
!> @param[in]           r_Psi           REAL array of DIMENSION(i_faceunknowns, i_faceunknowns);
!>                                      matrix of basis function evaluation at element quadrature points
!> @param[in]           r_Minv          REAL array of DIMENSION(i_faceunknowns, i_faceunknowns);
!>                                      inverse mass matrix interpolated on element quadrature points
!> @param[in]           r_visc          REAL; viscosity parameter
!> @param[in]           r_linear        INTEGER; switch for linearization (=1). default is 0
!> @param[in]           r_gamma         REAL; bottom friction parameter
!> @param[in]           r_gammat        REAL; wind friction parameter
!> @param[in]           r_rho           REAL; uniform density of fluid
!> @param[in]           r_coriol_e      REAL array of DIMENSION(i_faceunknowns);
!>                                      Coriolis forcing at element DOFs
!> @param[in]           r_taux_e        REAL array of DIMENSION(i_faceunknowns);
!>                                      x-component of wind forcing at element DOFs
!> @param[in]           r_tauy_e        REAL array of DIMENSION(i_faceunknowns);
!>                                      y-component of wind forcing at element DOFs

  FUNCTION divflux(r_h, r_hu, r_hv, r_wet, r_gradh, r_gradhu, r_gradhv) RESULT(r_Div)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), INTENT(in)                 :: r_hu, r_hv, r_h, r_wet
    REAL (KIND = GRID_SR), DIMENSION(2), INTENT(in)   :: r_gradhu, r_gradhv, r_gradh
    REAL (KIND = GRID_SR), DIMENSION(3)               :: r_Div

!--- local declarations
    REAL (KIND = GRID_SR)                             :: r_u, r_v
    REAL (KIND = GRID_SR), DIMENSION(3,2)             :: r_div_tmp

    r_u = velocity(r_h, r_hu)
    r_v = velocity(r_h, r_hv)

    r_div_tmp(1,1) = r_gradhu(1)
    r_div_tmp(1,2) = r_gradhv(2)

    r_div_tmp(2,1) = - r_u*r_u * r_gradh(1) + 2.0 * r_gradhu(1) * r_u &
                     + r_wet*r_h*r_gradh(1)

    r_div_tmp(2,2) = - r_u*r_v * r_gradh(2) + r_gradhu(2)*r_v + r_gradhv(2)*r_u

    r_div_tmp(3,1) = - r_u*r_v * r_gradh(1) + r_gradhu(1)*r_v + r_gradhv(1)*r_u

    r_div_tmp(3,2) = - r_v*r_v * r_gradh(2) + 2.0 * r_gradhv(2) * r_v &
                     + r_wet*r_h*r_gradh(2)

!--- compute full RHS
      r_Div(:) = r_div_tmp(:,1) + r_div_tmp(:,2)

  END FUNCTION divflux

END MODULE DG_equation

#endif
