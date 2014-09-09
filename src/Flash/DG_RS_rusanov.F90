!*******************************************************************************
!
!> @file  DG_RS_rusanov.F90
!> @brief contains module DG_riemann_solver
!
!*******************************************************************************
#include "Compilation_control.f90"


#if defined(_FLASH)

MODULE DG_riemann_solver

#define GRID_DIMENSION 2
  USE SFC_data_types
!   USE GRID_api
  USE DG_equation

  PRIVATE
  PUBLIC  :: riemannsolver, riemannsolverAdv

  CONTAINS

!*******************************************************************************
! DESCRIPTION of [FUNCTION riemannsolver]:
!> @brief Rusanov Riemann Solver for Shallow Water Equations
!>
!> @param[in]     r_h_l           REAL array of DIMENSION(i_faceunknowns);
!>                                height at DOFs of left element
!> @param[in]     r_h_r           REAL array of DIMENSION(i_faceunknowns);
!>                                height at DOFs of right element
!> @param[in]     r_hu_l          REAL array of DIMENSION(i_faceunknowns);
!>                                momentum in x-direction at DOFs of left element
!> @param[in]     r_hu_r          REAL array of DIMENSION(i_faceunknowns);
!>                                momentum in x-direction at DOFs of right element
!> @param[in]     r_hv_l          REAL array of DIMENSION(i_faceunknowns);
!>                                momentum in x-direction at DOFs of left element
!> @param[in]     r_hv_r          REAL array of DIMENSION(i_faceunknowns);
!>                                momentum in y-direction at DOFs of right element
!> @param[in]     r_normal        REAL array of DIMENSION(2);
!>                                normal vector pointing from left to right element
!>
!> @return                        REAL array of DIMENSION(3); Solution of Riemann problem
!>
  FUNCTION riemannsolver(r_h_l, r_h_r, r_hu_l, r_hu_r, r_hv_l, r_hv_r, &
                         r_normal) RESULT(r_flux)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(3)               :: r_flux
    REAL (KIND = GRID_SR), INTENT(in)                 :: r_h_l, r_h_r, &
                                                         r_hu_l, r_hu_r, r_hv_l, r_hv_r
    REAL (KIND = GRID_SR), DIMENSION(GRID_DIMENSION), INTENT(in) :: r_normal

    REAL (KIND = GRID_SR)                             :: r_u_l, r_u_r, r_v_l, r_v_r, &
                                                         r_a_l, r_a_r, r_u_norm_l, r_u_norm_r, &
                                                         r_abs_l, r_abs_r, r_lambda
    REAL (KIND = GRID_SR), DIMENSION(3)               :: r_dq
    REAL (KIND = GRID_SR), DIMENSION(3,2)             :: r_flux_left, r_flux_right, &
                                                         r_flux_help, r_fluxrus

!--- compute velocities and their normal components
    r_u_l = velocity(r_h_l, r_hu_l)
    r_v_l = velocity(r_h_l, r_hv_l)
    r_u_r = velocity(r_h_r, r_hu_r)
    r_v_r = velocity(r_h_r, r_hv_r)

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

    r_flux_help = r_lambda*matmul(reshape(r_dq, (/ 3, 1 /)), reshape(r_normal, (/ 1, 2 /)))

    r_fluxrus = (r_flux_left + r_flux_right - r_flux_help)
    r_fluxrus = r_fluxrus/2.0_GRID_SR

    r_flux = matmul(r_fluxrus, r_normal)

  END FUNCTION riemannsolver

!*******************************************************************************
! DESCRIPTION of [FUNCTION riemannsolverAdv]:
!> @brief Rusanov Riemann Solver for Shallow Water Equations
!>
!> @param[in]     r_h_l           REAL array of DIMENSION(i_faceunknowns);
!>                                height at DOFs of left element
!> @param[in]     r_h_r           REAL array of DIMENSION(i_faceunknowns);
!>                                height at DOFs of right element
!> @param[in]     r_hu_l          REAL array of DIMENSION(i_faceunknowns);
!>                                momentum in x-direction at DOFs of left element
!> @param[in]     r_hu_r          REAL array of DIMENSION(i_faceunknowns);
!>                                momentum in x-direction at DOFs of right element
!> @param[in]     r_hv_l          REAL array of DIMENSION(i_faceunknowns);
!>                                momentum in x-direction at DOFs of left element
!> @param[in]     r_hv_r          REAL array of DIMENSION(i_faceunknowns);
!>                                momentum in y-direction at DOFs of right element
!> @param[in]     r_normal        REAL array of DIMENSION(2);
!>                                normal vector pointing from left to right element
!>
!> @return                        REAL array of DIMENSION(3); Solution of Riemann problem
!>
  FUNCTION riemannsolverAdv(r_h_l, r_h_r, r_hu_l, r_hu_r, r_hv_l, r_hv_r, &
                            r_normal) RESULT(r_flux)

    IMPLICIT NONE

    REAL (KIND = GRID_SR), DIMENSION(3)               :: r_flux
    REAL (KIND = GRID_SR), INTENT(in)                 :: r_h_l, r_h_r, &
                                                         r_hu_l, r_hu_r, r_hv_l, r_hv_r
    REAL (KIND = GRID_SR), DIMENSION(GRID_DIMENSION), INTENT(in) :: r_normal

    REAL (KIND = GRID_SR)                             :: r_u_l, r_u_r, r_v_l, r_v_r, &
                                                         r_u_norm_l, r_u_norm_r, &
                                                         r_lambda
    REAL (KIND = GRID_SR), DIMENSION(3)               :: r_dq
    REAL (KIND = GRID_SR), DIMENSION(3,2)             :: r_flux_left, r_flux_right, &
                                                         r_flux_help, r_fluxrus

!--- compute velocities and their normal components
    r_u_l = velocity(r_h_l, r_hu_l)
    r_v_l = velocity(r_h_l, r_hv_l)
    r_u_r = velocity(r_h_r, r_hu_r)
    r_v_r = velocity(r_h_r, r_hv_r)

    r_u_norm_l =  r_u_l*r_normal(1) + r_v_l*r_normal(2)
    r_u_norm_r =  r_u_r*r_normal(1) + r_v_r*r_normal(2)

!--- compute maximal flux strength
    r_lambda = MAX(ABS(r_u_norm_l), ABS(r_u_norm_r))

!--- compute left and right fluxes
    r_flux_left  = flux(r_h_l, r_hu_l, r_hv_l, r_iswet=0.0_GRID_SR)
    r_flux_right = flux(r_h_r, r_hu_r, r_hv_r, r_iswet=0.0_GRID_SR)

!--- compute q_r - q_l as three-dimensional vector
    r_dq(1) = r_h_r  - r_h_l
    r_dq(2) = r_hu_r - r_hu_l
    r_dq(3) = r_hv_r - r_hv_l

    r_flux_help = r_lambda*MATMUL(RESHAPE(r_dq, (/ 3, 1 /)), RESHAPE(r_normal, (/ 1, 2 /)))

    r_fluxrus = (r_flux_left + r_flux_right - r_flux_help)
    r_fluxrus = r_fluxrus/2.0_GRID_SR

    r_flux = MATMUL(r_fluxrus, r_normal)

  END FUNCTION riemannsolverAdv

!*******************************************************************************
END MODULE DG_riemann_solver

#endif
