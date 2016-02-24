module c_bind_riemannsolvers
  implicit none

  enum, bind(c)
    enumerator :: GEOCLAW_FWAVE = 1, GEOCLAW_SSQ_FWAVE, GEOCLAW_AUG_RIEMANN
  end enum

  interface c_bind_geoclaw_solver
    module procedure c_bind_geoclaw_solver_qp
    module procedure c_bind_geoclaw_solver_dp
    module procedure c_bind_geoclaw_solver_sp
  end interface

  contains

  subroutine c_bind_geoclaw_solver_dp(i_solver, i_maxIter, i_numberOfFWaves, i_hL, i_hR, i_huL, i_huR, i_hvL, i_hvR, i_bL, i_bR, i_dryTol, i_g, o_netUpdatesLeft, o_netUpdatesRight, o_maxWaveSpeed )
  !variable declaration
    !input
      integer, intent(in)                                 :: i_solver, i_maxIter, i_numberOfFWaves

      double precision, intent(inout)                     :: i_hL,i_hR, i_huL,i_huR, i_hvL,i_hvR,i_bL,i_bR;
      double precision, intent(in)                        :: i_dryTol,i_g

    !output
      double precision, dimension(3), intent(out)         :: o_netUpdatesLeft, o_netUpdatesRight;
      double precision, intent(out)                       :: o_maxWaveSpeed;

    !local
      integer                                             :: waveNumber, equationNumber;
      double precision, dimension(i_numberOfFWaves)       :: waveSpeeds;
      double precision, dimension(3, i_numberOfFWaves)    :: fWaves;
      double precision, dimension(3)                      :: wall;
      double precision                                    :: hstar, hstartest, uL, uR, vL, vR, uhat, chat, delphi;
      double precision                                    :: sL, sR, sRoe1, sRoe2, sE1, sE2, s1m, s2m;
      logical                                             :: rare1,rare2;


  !******************************************************************
  !* necessary (changed) part of the GeoClaw subroutine rpn - start *
  !******************************************************************

  !reset max wave speed
  o_maxWaveSpeed = 0.d0

  !reset net updates
  o_netUpdatesLeft = 0.d0
  o_netUpdatesRight = 0.d0


  !Initialize Riemann problem for grid interface
  waveSpeeds=0.d0
  fWaves=0.d0

  !skip problem if in a completely dry area
  if (i_hL .le. i_dryTol .and. i_hR .le. i_dryTol) then
    return
  endif

  !check for wet/dry boundary
  if (i_hR.gt.i_dryTol) then
    uR=i_huR/i_hR
    vR=i_hvR/i_hR
  else
    i_hR = 0.d0
    i_huR = 0.d0
    i_hvR = 0.d0
    uR = 0.d0
    vR = 0.d0
  endif

  if (i_hL.gt.i_dryTol) then
    uL=i_huL/i_hL
    vL=i_hvL/i_hL
  else
    i_hL=0.d0
    i_huL=0.d0
    i_hvL=0.d0
    uL=0.d0
    vL=0.d0
  endif

  !per default there is no wall
  wall(1) = 1.d0
  wall(2) = 1.d0
  wall(3) = 1.d0

  if (i_hR.le.i_dryTol) then
    call riemanntype(i_hL,i_hL,uL,-uL,hstar,s1m,s2m,rare1,rare2,1,i_dryTol,i_g)
    hstartest=max(i_hL,hstar)
    if (hstartest+i_bL.lt.i_bR) then !right state should become ghost values that mirror left for wall problem
      wall(2)=0.d0
      wall(3)=0.d0
      i_hR=i_hL
      i_huR=-i_huL
      i_bR=i_bL
      uR=-uL
      vR=vL
    elseif (i_hL+i_bL.lt.i_bR) then
      i_bR=i_hL+i_bL
    endif
  elseif (i_hL.le.i_dryTol) then ! right surface is lower than left topo
    call riemanntype(i_hR,i_hR,-uR,uR,hstar,s1m,s2m,rare1,rare2,1,i_dryTol,i_g)
    hstartest=max(i_hR,hstar)
    if (hstartest+i_bR.lt.i_bL) then  !left state should become ghost values that mirror right
      wall(1)=0.d0
      wall(2)=0.d0
      i_hL=i_hR
      i_huL=-i_huR
      i_bL=i_bR
      uL=-uR
      vL=vR
    elseif (i_hR+i_bR.lt.i_bL) then
      i_bL=i_hR+i_bR
    endif
  endif

  !BUGFIX:
  !Problem: loss of significance may occur in phiR-phiL, causing divergence of the steady state.
  !Action:  Compute delphi=phiR-phiL explicitly. delphi is arithmetically equivalent to phiR-phiL, but with far smaller numerical loss.
  delphi = (i_huR - i_huL)*(uL + uR) - uL*uR*(i_hR - i_hL) + (0.5d0*i_g*(i_bR + i_hR - i_bL - i_hL)*(i_hR + i_hL)) - 0.5d0*i_g*(i_hR + i_hL)*(i_bR - i_bL)

  !determine wave speeds
  sL=uL-sqrt(i_g*i_hL) ! 1 wave speed of left state
  sR=uR+sqrt(i_g*i_hR) ! 2 wave speed of right state

  uhat=(sqrt(i_g*i_hL)*uL + sqrt(i_g*i_hR)*uR)/(sqrt(i_g*i_hR)+sqrt(i_g*i_hL)) ! Roe average
  chat=sqrt(i_g*0.5d0*(i_hR+i_hL)) ! Roe average
  sRoe1=uhat-chat ! Roe wave speed 1 wave
  sRoe2=uhat+chat ! Roe wave speed 2 wave

  sE1 = min(sL,sRoe1) ! Eindfeldt speed 1 wave
  sE2 = max(sR,sRoe2) ! Eindfeldt speed 2 wave

  !*******************
  !* call the solver *
  !*******************

  select case (i_solver)
    case (GEOCLAW_FWAVE)
        call riemann_fwave(3,i_numberOfFWaves,i_hL,i_hR,i_huL,i_huR,i_hvL,i_hvR,i_bL,i_bR,uL,uR,vL,vR,delphi,sE1,sE2,i_dryTol,i_g,waveSpeeds,fWaves)
    case (GEOCLAW_SSQ_FWAVE)
        call riemann_ssqfwave(i_maxIter,3,i_numberOfFWaves,i_hL,i_hR,i_huL,i_huR,i_hvL,i_hvR,i_bL,i_bR,uL,uR,vL,vR,delphi,sE1,sE2,i_dryTol,i_g,waveSpeeds,fWaves)
    case (GEOCLAW_AUG_RIEMANN)
        call riemann_aug_JCP(i_maxIter,3,i_numberOfFWaves,i_hL,i_hR,i_huL,i_huR,i_hvL,i_hvR,i_bL,i_bR,uL,uR,vL,vR,delphi,sE1,sE2,i_dryTol,i_g,waveSpeeds,fWaves)
  end select

  !eliminate ghost fluxes for wall
  do waveNumber=1,i_numberOfFWaves
    waveSpeeds(waveNumber)=waveSpeeds(waveNumber)*wall(waveNumber)
    do equationNumber=1,3
      fWaves(equationNumber,waveNumber)=fWaves(equationNumber,waveNumber)*wall(waveNumber)
    enddo
  enddo

  !compute net updates
  do equationNumber=1,3
    do  waveNumber=1,i_numberOfFWaves
      if (waveSpeeds(waveNumber).lt.0.d0) then
       o_netUpdatesLeft(equationNumber)=o_netUpdatesLeft(equationNumber) + fWaves(equationNumber,waveNumber);
      elseif (waveSpeeds(waveNumber).gt.0.d0) then
       o_netUpdatesRight(equationNumber)=o_netUpdatesRight(equationNumber) + fWaves(equationNumber,waveNumber);
      else
        o_netUpdatesLeft(equationNumber)=o_netUpdatesLeft(equationNumber) + .5d0*fWaves(equationNumber,waveNumber);
        o_netUpdatesRight(equationNumber)=o_netUpdatesRight(equationNumber) + .5d0*fWaves(equationNumber,waveNumber);
      endif
    enddo
  enddo

  !******************************************************
  !* necessary part of the GeoClaw subroutine rpn - end *
  !******************************************************

  !compute maximum wave speed
  waveSpeeds = abs(waveSpeeds)
  o_maxWaveSpeed = maxVal(waveSpeeds)
  end subroutine c_bind_geoclaw_solver_dp

  subroutine c_bind_geoclaw_solver_sp(i_solver, i_maxIter, i_numberOfFWaves, i_hL, i_hR, i_huL, i_huR, i_hvL, i_hvR, i_bL, i_bR, i_dryTol, i_g, o_netUpdatesLeft, o_netUpdatesRight, o_maxWaveSpeed )
  !variable declaration
    !input
      integer, intent(in)                                 :: i_solver, i_maxIter, i_numberOfFWaves

      real, intent(inout)                     :: i_hL,i_hR, i_huL,i_huR, i_hvL,i_hvR,i_bL,i_bR;
      real, intent(in)                        :: i_dryTol,i_g

    !output
      real, dimension(3), intent(out)         :: o_netUpdatesLeft, o_netUpdatesRight;
      real, intent(out)                       :: o_maxWaveSpeed;

    !local
      integer                                             :: waveNumber, equationNumber;
      real, dimension(i_numberOfFWaves)       :: waveSpeeds;
      real, dimension(3, i_numberOfFWaves)    :: fWaves;
      real, dimension(3)                      :: wall;
      real                                    :: hstar, hstartest, uL, uR, vL, vR, uhat, chat, delphi;
      real                                    :: sL, sR, sRoe1, sRoe2, sE1, sE2, s1m, s2m;
      logical                                             :: rare1,rare2;


  !******************************************************************
  !* necessary (changed) part of the GeoClaw subroutine rpn - start *
  !******************************************************************

  !reset max wave speed
  o_maxWaveSpeed = 0.e0;

  !reset net updates
  o_netUpdatesLeft = 0.e0;
  o_netUpdatesRight = 0.e0;


  !Initialize Riemann problem for grid interface
  waveSpeeds=0.e0
  fWaves=0.e0

  !skip problem if in a completely dry area
  if (i_hL .le. i_dryTol .and. i_hR .le. i_dryTol) then
    return;
  endif

  !BUGFIX:
  !Problem: loss of significance may occur in phiR-phiL, causing divergence of the steady state.
  !Action:  Compute delphi=phiR-phiL explicitly. delphi is arithmetically equivalent to phiR-phiL, but with far smaller numerical loss.
  delphi=0.0e0

  !check for wet/dry boundary
  if (i_hR.gt.i_dryTol) then
    uR=i_huR/i_hR
    vR=i_hvR/i_hR
  else
    i_hR = 0.e0
    i_huR = 0.e0
    i_hvR = 0.e0
    uR = 0.e0
    vR = 0.e0
  endif

  if (i_hL.gt.i_dryTol) then
    uL=i_huL/i_hL
    vL=i_hvL/i_hL
  else
    i_hL=0.e0
    i_huL=0.e0
    i_hvL=0.e0
    uL=0.e0
    vL=0.e0
  endif

  !per default there is no wall
  wall(1) = 1.e0
  wall(2) = 1.e0
  wall(3) = 1.e0

  if (i_hR.le.i_dryTol) then
    call riemanntype_sp(i_hL,i_hL,uL,-uL,hstar,s1m,s2m,rare1,rare2,1,i_dryTol,i_g)
    hstartest=max(i_hL,hstar)
    if (hstartest+i_bL.lt.i_bR) then !right state should become ghost values that mirror left for wall problem
      wall(2)=0.e0
      wall(3)=0.e0
      i_hR=i_hL
      i_huR=-i_huL
      i_bR=i_bL
      uR=-uL
      vR=vL
    elseif (i_hL+i_bL.lt.i_bR) then
      i_bR=i_hL+i_bL
    endif
  elseif (i_hL.le.i_dryTol) then ! right surface is lower than left topo
    call riemanntype_sp(i_hR,i_hR,-uR,uR,hstar,s1m,s2m,rare1,rare2,1,i_dryTol,i_g)
    hstartest=max(i_hR,hstar)
    if (hstartest+i_bR.lt.i_bL) then  !left state should become ghost values that mirror right
      wall(1)=0.e0
      wall(2)=0.e0
      i_hL=i_hR
      i_huL=-i_huR
      i_bL=i_bR
      uL=-uR
      vL=vR
    elseif (i_hR+i_bR.lt.i_bL) then
      i_bL=i_hR+i_bR
    endif
  endif

  !BUGFIX:
  !Problem: loss of significance may occur in phiR-phiL, causing divergence of the steady state.
  !Action:  Compute delphi=phiR-phiL explicitly. delphi is arithmetically equivalent to phiR-phiL, but with far smaller numerical loss.
  delphi = (i_huR - i_huL)*(uL + uR) - uL*uR*(i_hR - i_hL) + (0.5e0*i_g*(i_bR + i_hR - i_bL - i_hL)*(i_hR + i_hL)) - 0.5e0*i_g*(i_hR + i_hL)*(i_bR - i_bL)

  !determine wave speeds
  sL=uL-sqrt(i_g*i_hL) ! 1 wave speed of left state
  sR=uR+sqrt(i_g*i_hR) ! 2 wave speed of right state

  uhat=(sqrt(i_g*i_hL)*uL + sqrt(i_g*i_hR)*uR)/(sqrt(i_g*i_hR)+sqrt(i_g*i_hL)) ! Roe average
  chat=sqrt(i_g*0.5e0*(i_hR+i_hL)) ! Roe average
  sRoe1=uhat-chat ! Roe wave speed 1 wave
  sRoe2=uhat+chat ! Roe wave speed 2 wave

  sE1 = min(sL,sRoe1) ! Eindfeldt speed 1 wave
  sE2 = max(sR,sRoe2) ! Eindfeldt speed 2 wave

  !*******************
  !* call the solver *
  !*******************

  select case (i_solver)
    case (GEOCLAW_FWAVE)
        call riemann_fwave_sp(3,i_numberOfFWaves,i_hL,i_hR,i_huL,i_huR,i_hvL,i_hvR,i_bL,i_bR,uL,uR,vL,vR,delphi,sE1,sE2,i_dryTol,i_g,waveSpeeds,fWaves)
    case (GEOCLAW_SSQ_FWAVE)
        call riemann_ssqfwave_sp(i_maxIter,3,i_numberOfFWaves,i_hL,i_hR,i_huL,i_huR,i_hvL,i_hvR,i_bL,i_bR,uL,uR,vL,vR,delphi,sE1,sE2,i_dryTol,i_g,waveSpeeds,fWaves)
    case (GEOCLAW_AUG_RIEMANN)
        call riemann_aug_JCP_sp(i_maxIter,3,i_numberOfFWaves,i_hL,i_hR,i_huL,i_huR,i_hvL,i_hvR,i_bL,i_bR,uL,uR,vL,vR,delphi,sE1,sE2,i_dryTol,i_g,waveSpeeds,fWaves)
  end select

  !eliminate ghost fluxes for wall
  do waveNumber=1,i_numberOfFWaves
    waveSpeeds(waveNumber)=waveSpeeds(waveNumber)*wall(waveNumber)
    do equationNumber=1,3
      fWaves(equationNumber,waveNumber)=fWaves(equationNumber,waveNumber)*wall(waveNumber)
    enddo
  enddo

  !compute net updates
  do equationNumber=1,3
    do  waveNumber=1,i_numberOfFWaves
      if (waveSpeeds(waveNumber).lt.0.e0) then
       o_netUpdatesLeft(equationNumber)=o_netUpdatesLeft(equationNumber) + fWaves(equationNumber,waveNumber);
      elseif (waveSpeeds(waveNumber).gt.0.e0) then
       o_netUpdatesRight(equationNumber)=o_netUpdatesRight(equationNumber) + fWaves(equationNumber,waveNumber);
      else
        o_netUpdatesLeft(equationNumber)=o_netUpdatesLeft(equationNumber) + .5e0*fWaves(equationNumber,waveNumber);
        o_netUpdatesRight(equationNumber)=o_netUpdatesRight(equationNumber) + .5e0*fWaves(equationNumber,waveNumber);
      endif
    enddo
  enddo

  !******************************************************
  !* necessary part of the GeoClaw subroutine rpn - end *
  !******************************************************

  !compute maximum wave speed
  waveSpeeds = abs(waveSpeeds)
  o_maxWaveSpeed = maxVal(waveSpeeds)
  end subroutine c_bind_geoclaw_solver_sp

  subroutine c_bind_geoclaw_solver_qp(i_solver, i_maxIter, i_numberOfFWaves, i_hL, i_hR, i_huL, i_huR, i_hvL, i_hvR, i_bL, i_bR, i_dryTol, i_g, o_netUpdatesLeft, o_netUpdatesRight, o_maxWaveSpeed )
  !variable declaration
    !input
      integer, parameter  :: QP = kind(1.0q0)
      integer, intent(in)                                 :: i_solver, i_maxIter, i_numberOfFWaves

      real (kind = QP), intent(inout)                     :: i_hL,i_hR, i_huL,i_huR, i_hvL,i_hvR,i_bL,i_bR;
      real (kind = QP), intent(in)                        :: i_dryTol,i_g

    !output
      real (kind = QP), dimension(3), intent(out)         :: o_netUpdatesLeft, o_netUpdatesRight;
      real (kind = QP), intent(out)                       :: o_maxWaveSpeed;

    !local
      integer                                             :: waveNumber, equationNumber;
      real (kind = QP), dimension(i_numberOfFWaves)       :: waveSpeeds;
      real (kind = QP), dimension(3, i_numberOfFWaves)    :: fWaves;
      real (kind = QP), dimension(3)                      :: wall;
      real (kind = QP)                                    :: hstar, hstartest, uL, uR, vL, vR, uhat, chat, delphi;
      real (kind = QP)                                    :: sL, sR, sRoe1, sRoe2, sE1, sE2, s1m, s2m;
      logical                                             :: rare1,rare2;


  !******************************************************************
  !* necessary (changed) part of the GeoClaw subroutine rpn - start *
  !******************************************************************

  !reset max wave speed
  o_maxWaveSpeed = 0.q0;

  !reset net updates
  o_netUpdatesLeft = 0.q0;
  o_netUpdatesRight = 0.q0;


  !Initialize Riemann problem for grid interface
  waveSpeeds=0.q0
  fWaves=0.q0

  !skip problem if in a completely dry area
  if (i_hL .le. i_dryTol .and. i_hR .le. i_dryTol) then
    return;
  endif


  !BUGFIX:
  !Problem: loss of significance may occur in phiR-phiL, causing divergence of the steady state.
  !Action:  Compute delphi=phiR-phiL explicitly. delphi is arithmetically equivalent to phiR-phiL, but with far smaller numerical loss.
  delphi=0.0q0

  !check for wet/dry boundary
  if (i_hR.gt.i_dryTol) then
    uR=i_huR/i_hR
    vR=i_hvR/i_hR
  else
    i_hR = 0.q0
    i_huR = 0.q0
    i_hvR = 0.q0
    uR = 0.q0
    vR = 0.q0
  endif

  if (i_hL.gt.i_dryTol) then
    uL=i_huL/i_hL
    vL=i_hvL/i_hL
  else
    i_hL=0.q0
    i_huL=0.q0
    i_hvL=0.q0
    uL=0.q0
    vL=0.q0
  endif

  !per default there is no wall
  wall(1) = 1.q0
  wall(2) = 1.q0
  wall(3) = 1.q0

  if (i_hR.le.i_dryTol) then
    call riemanntype_qp(i_hL,i_hL,uL,-uL,hstar,s1m,s2m,rare1,rare2,1,i_dryTol,i_g)
    hstartest=max(i_hL,hstar)
    if (hstartest+i_bL.lt.i_bR) then !right state should become ghost values that mirror left for wall problem
      wall(2)=0.q0
      wall(3)=0.q0
      i_hR=i_hL
      i_huR=-i_huL
      i_bR=i_bL
      uR=-uL
      vR=vL
    elseif (i_hL+i_bL.lt.i_bR) then
      i_bR=i_hL+i_bL
    endif
  elseif (i_hL.le.i_dryTol) then ! right surface is lower than left topo
    call riemanntype_qp(i_hR,i_hR,-uR,uR,hstar,s1m,s2m,rare1,rare2,1,i_dryTol,i_g)
    hstartest=max(i_hR,hstar)
    if (hstartest+i_bR.lt.i_bL) then  !left state should become ghost values that mirror right
      wall(1)=0.q0
      wall(2)=0.q0
      i_hL=i_hR
      i_huL=-i_huR
      i_bL=i_bR
      uL=-uR
      vL=vR
    elseif (i_hR+i_bR.lt.i_bL) then
      i_bL=i_hR+i_bR
    endif
  endif

  !BUGFIX:
  !Problem: loss of significance may occur in phiR-phiL, causing divergence of the steady state.
  !Action:  Compute delphi=phiR-phiL explicitly. delphi is arithmetically equivalent to phiR-phiL, but with far smaller numerical loss.
  delphi = (i_huR - i_huL)*(uL + uR) - uL*uR*(i_hR - i_hL) + (0.5q0*i_g*(i_bR + i_hR - i_bL - i_hL)*(i_hR + i_hL)) - 0.5q0*i_g*(i_hR + i_hL)*(i_bR - i_bL)

  !determine wave speeds
  sL=uL-sqrt(i_g*i_hL) ! 1 wave speed of left state
  sR=uR+sqrt(i_g*i_hR) ! 2 wave speed of right state

  uhat=(sqrt(i_g*i_hL)*uL + sqrt(i_g*i_hR)*uR)/(sqrt(i_g*i_hR)+sqrt(i_g*i_hL)) ! Roe average
  chat=sqrt(i_g*0.5q0*(i_hR+i_hL)) ! Roe average
  sRoe1=uhat-chat ! Roe wave speed 1 wave
  sRoe2=uhat+chat ! Roe wave speed 2 wave

  sE1 = min(sL,sRoe1) ! Eindfeldt speed 1 wave
  sE2 = max(sR,sRoe2) ! Eindfeldt speed 2 wave

  !*******************
  !* call the solver *
  !*******************

  select case (i_solver)
    case (GEOCLAW_FWAVE)
        call riemann_fwave_qp(3,i_numberOfFWaves,i_hL,i_hR,i_huL,i_huR,i_hvL,i_hvR,i_bL,i_bR,uL,uR,vL,vR,delphi,sE1,sE2,i_dryTol,i_g,waveSpeeds,fWaves)
    case (GEOCLAW_SSQ_FWAVE)
        call riemann_ssqfwave_qp(i_maxIter,3,i_numberOfFWaves,i_hL,i_hR,i_huL,i_huR,i_hvL,i_hvR,i_bL,i_bR,uL,uR,vL,vR,delphi,sE1,sE2,i_dryTol,i_g,waveSpeeds,fWaves)
    case (GEOCLAW_AUG_RIEMANN)
        call riemann_aug_JCP_qp(i_maxIter,3,i_numberOfFWaves,i_hL,i_hR,i_huL,i_huR,i_hvL,i_hvR,i_bL,i_bR,uL,uR,vL,vR,delphi,sE1,sE2,i_dryTol,i_g,waveSpeeds,fWaves)
  end select

  !eliminate ghost fluxes for wall
  do waveNumber=1,i_numberOfFWaves
    waveSpeeds(waveNumber)=waveSpeeds(waveNumber)*wall(waveNumber)
    do equationNumber=1,3
      fWaves(equationNumber,waveNumber)=fWaves(equationNumber,waveNumber)*wall(waveNumber)
    enddo
  enddo

  !compute net updates
  do equationNumber=1,3
    do  waveNumber=1,i_numberOfFWaves
      if (waveSpeeds(waveNumber).lt.0.q0) then
       o_netUpdatesLeft(equationNumber)=o_netUpdatesLeft(equationNumber) + fWaves(equationNumber,waveNumber);
      elseif (waveSpeeds(waveNumber).gt.0.q0) then
       o_netUpdatesRight(equationNumber)=o_netUpdatesRight(equationNumber) + fWaves(equationNumber,waveNumber);
      else
        o_netUpdatesLeft(equationNumber)=o_netUpdatesLeft(equationNumber) + .5q0*fWaves(equationNumber,waveNumber);
        o_netUpdatesRight(equationNumber)=o_netUpdatesRight(equationNumber) + .5q0*fWaves(equationNumber,waveNumber);
      endif
    enddo
  enddo

  !******************************************************
  !* necessary part of the GeoClaw subroutine rpn - end *
  !******************************************************

  !compute maximum wave speed
  waveSpeeds = abs(waveSpeeds)
  o_maxWaveSpeed = maxVal(waveSpeeds)
  end subroutine c_bind_geoclaw_solver_qp
end module
