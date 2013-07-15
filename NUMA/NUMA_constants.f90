! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


!----------------------------------------------------------------------!
! This module defines the Constants that were in PARAM.H
! F.X. Giraldo
! Department of Applied Mathematics
! Naval Postgraduate School
! Monterey, CA 93943-5216
!----------------------------------------------------------------------!
#if defined(_NUMA)
module NUMA_constants
  use NUMA_data_types
  public :: &
       mod_constants_create, &
       pi, tol, earth_radius, omega, gravity, gamma, cp, cv, rgas, xappa, p00, nnorm
  
  private
  !-----------------------------------------------------------------------
  !module variables and parameters
  real (kind = GRID_SR) pi, tol, earth_radius, omega, gravity, gamma, cp, cv, rgas, xappa, p00
  integer nnorm
  !-----------------------------------------------------------------------
  
contains
  
  !-----------------------------------------------------------------------
  subroutine mod_constants_create()
    
    implicit none

    !Initialize
    pi=4.0*atan(1.0)
    tol=1e-10
!    earth_radius=6.37122e6
    earth_radius=160
    omega=7.29e-5
    gravity=9.80616
    gamma=1.4
    cp=1004.67
    cv=717.5
    rgas=287.17
    xappa=0.286
    p00=1e5
    nnorm=1000
    
  end subroutine mod_constants_create
  
end module NUMA_constants
#endif
