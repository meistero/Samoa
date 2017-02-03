!*******************************************************************************
!
!> @file  DG_limiter_BJSV.F90
!> @brief contains module DG_limiter
!
!*******************************************************************************
MODULE DG_limiter

  USE GRID_api
  USE FLASH_parameters
  USE DG_equation

  PRIVATE
  PUBLIC  :: limiter

  CONTAINS

!*******************************************************************************
!> @brief flux limiter for DG method
!>
!>  Slope limiter of BArth and Jesperson combined with Wetting and Drying fix
!>
!> @param[in]     p_ghand             TYPE(grid_handle); grid handle
!> @param[inout]  r_Q                 REAL array of DIMENSION(i_numelmt, i_faceunknowns,3);
!>                                    read-in flux and limited version of the flux
!> @param[in]     i_elementedges      INTEGER array of DIMENSION(3,i_numelmt);
!>                                    all three edges per element
!> @param[in]     i_edgeinfo          INTEGER array of DIMENSION(2,i_numedge);
!>                                    number of left and right element of the edge
!> @param[in]     r_edgelength        REAL array of DIMENSION(i_numedge);
!>                                    edge length
!> @param[in]     i_numelmt           INTEGER; number of elements
!> @param[in]     r_coo               REAL array of DIMENSION(2,i_numdofs);
!>                                    DOF coordinates
!> @param[in]     i_edgenodes         INTEGER array of DIMENSION(i_numedge,2);
!>                                    nodes per edge
!> @param[in]     i_elementnodes      INTEGER array of DIMENSION(i_numelmt,3);
!>                                    nodes per element
!> @param[in]     i_numedge           INTEGER; number of edges
!> @param[in]     r_bathy             REAL array of DIMENSION(i_numelmt, i_faceunknowns)
!>                                    bathymetry of DOFs
!> @param[in]     i_faceunknowns      INTEGER; number of DOFs
!>
  SUBROUTINE limiter(r_Q, i_elementedges, i_iedgeinfo, i_bedgeinfo, &
                     r_edgelength, r_metrics_inv, i_numelmt, r_coo, r_coodofs, &
                     i_edgenodes, i_elementnodes, i_elementdofs, i_numedge, &
                     r_bathy, i_faceunknowns, r_TV)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)                            :: p_ghand
    REAL(KIND = GRID_SR), DIMENSION(:,:), INTENT(inout)       :: r_Q
    INTEGER(KIND = GRID_SI), DIMENSION(:,:), INTENT(in)       :: i_elementedges, i_iedgeinfo, &
                                                                 i_bedgeinfo, i_edgenodes, &
                                                                 i_elementnodes, i_elementdofs
    REAL(KIND = GRID_SR), DIMENSION(:), INTENT(in)            :: r_edgelength, r_bathy
    REAL(KIND = GRID_SR), DIMENSION(:,:,:), INTENT(in)        :: r_metrics_inv
    REAL(KIND = GRID_SR), DIMENSION(:,:), INTENT(in)          :: r_coo, r_coodofs
    INTEGER(KIND = GRID_SI), INTENT(in)                       :: i_numelmt, i_faceunknowns, &
                                                                 i_numedge
    REAL(KIND = GRID_SR), DIMENSION(:,:), INTENT(out)         :: r_TV

    INTEGER(KIND = GRID_SI)                                   :: i_alct, i_edge, i_elmt
    REAL(KIND = GRID_SR)                                      :: r_alpha, r_hext, r_utmp
    INTEGER(KIND = GRID_SI), DIMENSION(:), ALLOCATABLE        :: i_dofssorted
    REAL(KIND = GRID_SR), DIMENSION(:), ALLOCATABLE           :: r_corr, r_htmp
    REAL(KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE         :: r_evander, r_emean, r_emin, r_emax, &
                                                                 r_ulim, r_Qlim
    LOGICAL, DIMENSION(:), ALLOCATABLE                        :: l_loc

    r_TV = 0._GRID_SR

!--- allocate workspace
    ALLOCATE(i_dofssorted(i_faceunknowns), &
             r_corr(i_faceunknowns), &
             r_htmp(i_faceunknowns), &
             r_evander(i_faceunknowns, i_faceunknowns), &
             r_emean(6, i_numelmt), &
             r_emin(6, i_numelmt), &
             r_emax(6, i_numelmt), &
             r_ulim(i_faceunknowns, 3), &
             r_Qlim(3, i_numelmt*i_faceunknowns), &
             l_loc(i_faceunknowns), stat=i_alct)
    IF (i_alct /= 0) CALL grid_error(c_error='[DG_limiter_BJSV]: Could not allocate enough memory')

!--- get info from signature and initialize
    r_evander = GRID_femtypes%p_type(FEM_ZETA)%sig%r_evander

!--- compute mean, min and max values for each element
    DO i_elmt=1, i_numelmt
      r_emean(1,i_elmt) = DOT_PRODUCT(r_evander(:,1), r_Q(1,i_elementdofs(:,i_elmt))) / SUM(r_evander(:,1))
      r_emean(2,i_elmt) = DOT_PRODUCT(r_evander(:,1), r_Q(2,i_elementdofs(:,i_elmt))) / SUM(r_evander(:,1))
      r_emean(3,i_elmt) = DOT_PRODUCT(r_evander(:,1), r_Q(3,i_elementdofs(:,i_elmt))) / SUM(r_evander(:,1))
      r_emean(4,i_elmt) = DOT_PRODUCT(r_evander(:,1), r_Q(1,i_elementdofs(:,i_elmt)) + &
                                      r_bathy(i_elementdofs(:,i_elmt))) / SUM(r_evander(:,1))
      r_emean(5,i_elmt) = velocity(r_emean(1,i_elmt), r_emean(2,i_elmt))
      r_emean(6,i_elmt) = velocity(r_emean(1,i_elmt), r_emean(3,i_elmt))

      r_emin(1,i_elmt)  = r_emean(1,i_elmt)
      r_emin(2,i_elmt)  = r_emean(2,i_elmt)
      r_emin(3,i_elmt)  = r_emean(3,i_elmt)
      r_emin(4,i_elmt)  = r_emean(4,i_elmt)
      r_emin(5,i_elmt)  = r_emean(5,i_elmt)
      r_emin(6,i_elmt)  = r_emean(6,i_elmt)

      r_emax(1,i_elmt)  = r_emean(1,i_elmt)
      r_emax(2,i_elmt)  = r_emean(2,i_elmt)
      r_emax(3,i_elmt)  = r_emean(3,i_elmt)
      r_emax(4,i_elmt)  = r_emean(4,i_elmt)
      r_emax(5,i_elmt)  = r_emean(5,i_elmt)
      r_emax(6,i_elmt)  = r_emean(6,i_elmt)
    END DO


!--- loop over all inner edges to update min and max
    DO i_edge=1, SIZE(i_iedgeinfo, 2)
      r_emin(1,i_iedgeinfo(3,i_edge)) = MIN(r_emin(1,i_iedgeinfo(3,i_edge)), r_emean(1,i_iedgeinfo(4,i_edge)))
      r_emin(1,i_iedgeinfo(4,i_edge)) = MIN(r_emin(1,i_iedgeinfo(4,i_edge)), r_emean(1,i_iedgeinfo(3,i_edge)))
      r_emin(2,i_iedgeinfo(3,i_edge)) = MIN(r_emin(2,i_iedgeinfo(3,i_edge)), r_emean(2,i_iedgeinfo(4,i_edge)))
      r_emin(2,i_iedgeinfo(4,i_edge)) = MIN(r_emin(2,i_iedgeinfo(4,i_edge)), r_emean(2,i_iedgeinfo(3,i_edge)))
      r_emin(3,i_iedgeinfo(3,i_edge)) = MIN(r_emin(3,i_iedgeinfo(3,i_edge)), r_emean(3,i_iedgeinfo(4,i_edge)))
      r_emin(3,i_iedgeinfo(4,i_edge)) = MIN(r_emin(3,i_iedgeinfo(4,i_edge)), r_emean(3,i_iedgeinfo(3,i_edge)))
      r_emin(4,i_iedgeinfo(3,i_edge)) = MIN(r_emin(4,i_iedgeinfo(3,i_edge)), r_emean(4,i_iedgeinfo(4,i_edge)))
      r_emin(4,i_iedgeinfo(4,i_edge)) = MIN(r_emin(4,i_iedgeinfo(4,i_edge)), r_emean(4,i_iedgeinfo(3,i_edge)))
      r_emin(5,i_iedgeinfo(3,i_edge)) = MIN(r_emin(5,i_iedgeinfo(3,i_edge)), r_emean(5,i_iedgeinfo(4,i_edge)))
      r_emin(5,i_iedgeinfo(4,i_edge)) = MIN(r_emin(5,i_iedgeinfo(4,i_edge)), r_emean(5,i_iedgeinfo(3,i_edge)))
      r_emin(6,i_iedgeinfo(3,i_edge)) = MIN(r_emin(6,i_iedgeinfo(3,i_edge)), r_emean(6,i_iedgeinfo(4,i_edge)))
      r_emin(6,i_iedgeinfo(4,i_edge)) = MIN(r_emin(6,i_iedgeinfo(4,i_edge)), r_emean(6,i_iedgeinfo(3,i_edge)))

      r_emax(1,i_iedgeinfo(3,i_edge)) = MAX(r_emax(1,i_iedgeinfo(3,i_edge)), r_emean(1,i_iedgeinfo(4,i_edge)))
      r_emax(1,i_iedgeinfo(4,i_edge)) = MAX(r_emax(1,i_iedgeinfo(4,i_edge)), r_emean(1,i_iedgeinfo(3,i_edge)))
      r_emax(2,i_iedgeinfo(3,i_edge)) = MAX(r_emax(2,i_iedgeinfo(3,i_edge)), r_emean(2,i_iedgeinfo(4,i_edge)))
      r_emax(2,i_iedgeinfo(4,i_edge)) = MAX(r_emax(2,i_iedgeinfo(4,i_edge)), r_emean(2,i_iedgeinfo(3,i_edge)))
      r_emax(3,i_iedgeinfo(3,i_edge)) = MAX(r_emax(3,i_iedgeinfo(3,i_edge)), r_emean(3,i_iedgeinfo(4,i_edge)))
      r_emax(3,i_iedgeinfo(4,i_edge)) = MAX(r_emax(3,i_iedgeinfo(4,i_edge)), r_emean(3,i_iedgeinfo(3,i_edge)))
      r_emax(4,i_iedgeinfo(3,i_edge)) = MAX(r_emax(4,i_iedgeinfo(3,i_edge)), r_emean(4,i_iedgeinfo(4,i_edge)))
      r_emax(4,i_iedgeinfo(4,i_edge)) = MAX(r_emax(4,i_iedgeinfo(4,i_edge)), r_emean(4,i_iedgeinfo(3,i_edge)))
      r_emax(5,i_iedgeinfo(3,i_edge)) = MIN(r_emax(5,i_iedgeinfo(3,i_edge)), r_emean(5,i_iedgeinfo(4,i_edge)))
      r_emax(5,i_iedgeinfo(4,i_edge)) = MIN(r_emax(5,i_iedgeinfo(4,i_edge)), r_emean(5,i_iedgeinfo(3,i_edge)))
      r_emax(6,i_iedgeinfo(3,i_edge)) = MIN(r_emax(6,i_iedgeinfo(3,i_edge)), r_emean(6,i_iedgeinfo(4,i_edge)))
      r_emax(6,i_iedgeinfo(4,i_edge)) = MIN(r_emax(6,i_iedgeinfo(4,i_edge)), r_emean(6,i_iedgeinfo(3,i_edge)))
    END DO

!--- loop over all elements to do the limiting
    DO i_elmt=1, i_numelmt

!--- limit in surface elevation
      r_hext = MAXVAL(r_Q(1,i_elementdofs(:,i_elmt))+r_bathy(i_elementdofs(:,i_elmt)))
      r_alpha = 1.0_GRID_SR
      IF (r_hext > r_emean(4,i_elmt)+1e-16) &
        r_alpha = MIN(r_alpha, (r_emax(4,i_elmt) - r_emean(4,i_elmt)) / (r_hext - r_emean(4,i_elmt)))
      r_hext = MINVAL(r_Q(1,i_elementdofs(:,i_elmt))+r_bathy(i_elementdofs(:,i_elmt)))
      IF (r_hext < r_emean(4,i_elmt)-1e-16) &
        r_alpha = MIN(r_alpha, (r_emean(4,i_elmt) - r_emin(4,i_elmt)) / (r_emean(4,i_elmt) - r_hext))

      r_corr    = (r_alpha-1.0_GRID_SR) * (r_Q(1,i_elementdofs(:,i_elmt))+r_bathy(i_elementdofs(:,i_elmt)) - r_emean(4,i_elmt))
      r_corr(3) = -(r_corr(1) + r_corr(2))

      r_Qlim(1,i_elementdofs(:,i_elmt)) = r_Q(1,i_elementdofs(:,i_elmt)) + r_corr

!--- positivity preserving limiter ala Bunya et al. (2009)
      IF(MINVAL(r_Qlim(1,i_elementdofs(:,i_elmt))) < 0.0_GRID_SR) THEN
        l_loc = .TRUE.
        i_dofssorted(1) = MINLOC(r_Qlim(1,i_elementdofs(:,i_elmt)), 1, l_loc)
        l_loc(i_dofssorted(1)) = .FALSE.
        i_dofssorted(2) = MINLOC(r_Qlim(1,i_elementdofs(:,i_elmt)), 1, l_loc)
        l_loc(i_dofssorted(2)) = .FALSE.
        i_dofssorted(3) = MINLOC(r_Qlim(1,i_elementdofs(:,i_elmt)), 1, l_loc)

        r_htmp = r_Qlim(1,i_elementdofs(:,i_elmt))
        r_Qlim(1,i_elementdofs(i_dofssorted(1),i_elmt)) = 0.0_GRID_SR
        r_Qlim(1,i_elementdofs(i_dofssorted(2),i_elmt)) = &
          MAX(0.0_GRID_SR, r_htmp(i_dofssorted(2)) + r_htmp(i_dofssorted(1))/2.0_GRID_SR)
        r_Qlim(1,i_elementdofs(i_dofssorted(3),i_elmt)) = &
          MAX(0.0_GRID_SR, r_htmp(i_dofssorted(3)) + r_htmp(i_dofssorted(1)) - &     ! Note: the max-function is a fix to enforce positivity,
          (r_Qlim(1,i_elementdofs(i_dofssorted(2),i_elmt))-r_htmp(i_dofssorted(2)))) ! it might affect mass conservation
      END IF

!--- limit momentum based on velocity
!--- compute limited u-velocity
      r_utmp      = velocity(r_Q(1,i_elementdofs(1,i_elmt)), r_Q(2,i_elementdofs(1,i_elmt)))
      r_ulim(1,:) = MAX(MIN(r_utmp, r_emax(5, i_elmt)), r_emin(5, i_elmt))
      r_utmp      = velocity(r_Q(1,i_elementdofs(2,i_elmt)), r_Q(2,i_elementdofs(2,i_elmt)))
      r_ulim(2,:) = MAX(MIN(r_utmp, r_emax(5, i_elmt)), r_emin(5, i_elmt))
      r_utmp      = velocity(r_Q(1,i_elementdofs(3,i_elmt)), r_Q(2,i_elementdofs(3,i_elmt)))
      r_ulim(3,:) = MAX(MIN(r_utmp, r_emax(5, i_elmt)), r_emin(5, i_elmt))

!--- compute u-velocities at other nodes based on linear momentum assumption
      r_ulim(1,1) = velocity(r_Qlim(1,i_elementdofs(1,i_elmt)), &
                             3.0_GRID_SR*r_emean(2,i_elmt)-r_Qlim(1,i_elementdofs(2,i_elmt))*r_ulim(2,1) &
                                                          -r_Qlim(1,i_elementdofs(3,i_elmt))*r_ulim(3,1))
      r_ulim(2,2) = velocity(r_Qlim(1,i_elementdofs(2,i_elmt)), &
                             3.0_GRID_SR*r_emean(2,i_elmt)-r_Qlim(1,i_elementdofs(1,i_elmt))*r_ulim(1,2) &
                                                          -r_Qlim(1,i_elementdofs(3,i_elmt))*r_ulim(3,2))
      r_ulim(3,3) = velocity(r_Qlim(1,i_elementdofs(3,i_elmt)), &
                             3.0_GRID_SR*r_emean(2,i_elmt)-r_Qlim(1,i_elementdofs(1,i_elmt))*r_ulim(1,3) &
                                                          -r_Qlim(1,i_elementdofs(2,i_elmt))*r_ulim(2,3))

      IF((MAXVAL(r_ulim(:,1))-MINVAL(r_ulim(:,1)) <= MAXVAL(r_ulim(:,2))-MINVAL(r_ulim(:,2))) .AND. &
         (MAXVAL(r_ulim(:,1))-MINVAL(r_ulim(:,1)) <= MAXVAL(r_ulim(:,3))-MINVAL(r_ulim(:,3)))) THEN
        r_Qlim(2,i_elementdofs(:,i_elmt)) = r_Qlim(1,i_elementdofs(:,i_elmt))*r_ulim(:,1)
      ELSE
        IF(MAXVAL(r_ulim(:,2))-MINVAL(r_ulim(:,2)) <= MAXVAL(r_ulim(:,3))-MINVAL(r_ulim(:,3))) THEN
          r_Qlim(2,i_elementdofs(:,i_elmt)) = r_Qlim(1,i_elementdofs(:,i_elmt))*r_ulim(:,2)
        ELSE
          r_Qlim(2,i_elementdofs(:,i_elmt)) = r_Qlim(1,i_elementdofs(:,i_elmt))*r_ulim(:,3)
        END IF
      END IF

!--- compute limited v-velocity
      r_utmp      = velocity(r_Q(1,i_elementdofs(1,i_elmt)), r_Q(3,i_elementdofs(1,i_elmt)))
      r_ulim(1,:) = MAX(MIN(r_utmp, r_emax(6, i_elmt)), r_emin(6, i_elmt))
      r_utmp      = velocity(r_Q(1,i_elementdofs(2,i_elmt)), r_Q(3,i_elementdofs(2,i_elmt)))
      r_ulim(2,:) = MAX(MIN(r_utmp, r_emax(6, i_elmt)), r_emin(6, i_elmt))
      r_utmp      = velocity(r_Q(1,i_elementdofs(3,i_elmt)), r_Q(3,i_elementdofs(3,i_elmt)))
      r_ulim(3,:) = MAX(MIN(r_utmp, r_emax(6, i_elmt)), r_emin(6, i_elmt))

!--- compute v-velocities at other nodes based on linear momentum assumption
      r_ulim(1,1) = velocity(r_Qlim(1,i_elementdofs(1,i_elmt)), &
                             3.0_GRID_SR*r_emean(3,i_elmt)-r_Qlim(1,i_elementdofs(2,i_elmt))*r_ulim(2,1) &
                                                          -r_Qlim(1,i_elementdofs(3,i_elmt))*r_ulim(3,1))
      r_ulim(2,2) = velocity(r_Qlim(1,i_elementdofs(2,i_elmt)), &
                             3.0_GRID_SR*r_emean(3,i_elmt)-r_Qlim(1,i_elementdofs(1,i_elmt))*r_ulim(1,2) &
                                                          -r_Qlim(1,i_elementdofs(3,i_elmt))*r_ulim(3,2))
      r_ulim(3,3) = velocity(r_Qlim(1,i_elementdofs(3,i_elmt)), &
                             3.0_GRID_SR*r_emean(3,i_elmt)-r_Qlim(1,i_elementdofs(1,i_elmt))*r_ulim(1,3) &
                                                          -r_Qlim(1,i_elementdofs(2,i_elmt))*r_ulim(2,3))

      IF((MAXVAL(r_ulim(:,1))-MINVAL(r_ulim(:,1)) <= MAXVAL(r_ulim(:,2))-MINVAL(r_ulim(:,2))) .AND. &
         (MAXVAL(r_ulim(:,1))-MINVAL(r_ulim(:,1)) <= MAXVAL(r_ulim(:,3))-MINVAL(r_ulim(:,3)))) THEN
        r_Qlim(3,i_elementdofs(:,i_elmt)) = r_Qlim(1,i_elementdofs(:,i_elmt))*r_ulim(:,1)
      ELSE
        IF(MAXVAL(r_ulim(:,2))-MINVAL(r_ulim(:,2)) <= MAXVAL(r_ulim(:,3))-MINVAL(r_ulim(:,3))) THEN
          r_Qlim(3,i_elementdofs(:,i_elmt)) = r_Qlim(1,i_elementdofs(:,i_elmt))*r_ulim(:,2)
        ELSE
          r_Qlim(3,i_elementdofs(:,i_elmt)) = r_Qlim(1,i_elementdofs(:,i_elmt))*r_ulim(:,3)
        END IF
      END IF

    END DO

    r_Q(:,:) = r_Qlim(:,:)

!--- deallocate workspace
    DEALLOCATE(i_dofssorted, r_corr, r_htmp, r_evander, r_emean, r_emin, r_emax, r_ulim, r_Qlim, l_loc)

  END SUBROUTINE limiter

!*******************************************************************************
END MODULE DG_limiter
