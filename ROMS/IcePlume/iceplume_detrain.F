#include "cppdefs.h"
#ifdef ICEPLUME
!
! ==================================================================!
!                                                                   !
! Use this function to calculate the velocity of detrainment.       !
!                                                                   !
! ==================================================================!
!
SUBROUTINE ICEPLUME_DETRAIN(ng, I,                            &
                          & iceDepthK, plumeDepthK, osDepthK, &
                          & lc, detr)
  USE mod_iceplume
  implicit none
!
! ==================================================================!
!                                                                   !
! Input/Output variables:                                           !
!                                                                   !
! ==================================================================!
!                                                                   !
! ng          - grid identifier                                     !
! I           - point source index identifier                       !
! iceDepthK   - glacier grounding line depth vertical grid index    !
! plumeDepthK - plume detrainment depth vertical grid index         !
! osDepthK    - plume overshoot depth vertical grid index           !
! lc          - plume/ocean contact length [m]                      !
! detr        - detrainment volume flux [m^3 s^-1]                  !
!                                                                   !
! ==================================================================!
!                                                                   !
! Local variables:                                                  !
!                                                                   !
! ==================================================================!
!                                                                   !
! detrDz      - detrainment layer total thickness [m]               !
! minDz       - minimum detrainment layer thickness [m]             !
! maxVel      - maximum velocity based on Richardson number         !
!               criterion [m s^-1]                                  !
! Fr          - Froude number                                       !
! rhoP        - plume density [kg m^-3]                             !
! rho1        - upper layer density [kg m^-3]                       !
! rho2        - lower layer density [kg m^-3]                       !
! h1          - upper layer thickness [m]                           !
! h2          - lower layer thickness [m]                           !
! gRed        - reduced gravity [m s^-2]                            !
! isSurface   - logical, if the plume discharge at surface          !
!                                                                   !
! udSwitch    - switch to determine search one layer up(1)/down(-1) !
! potE        - potential energy of plume/ambient water             !
!                                                                   !
! detrVel     - detrainment velocity [m s^-1]                       !
! detrWeight  - detrainment weight function                         !
! detrSum     - detrainment volume flux before normalize [m^3 s^-1] !
!                                                                   !
! ==================================================================!
!
! In/out variables
!
  integer, intent(in) :: ng, I
  integer, intent(in) :: iceDepthK, plumeDepthK, osDepthK
  real(r8), intent(in) :: lc, detr
!
! Local variables declaration
!
  real(r8) :: detrDz, minDz, maxVel, Fr
  real(r8) :: rhoP, rho1, rho2, h1, h2, gRed, potE
  integer  :: udSwitch
  logical  :: isSurface = .FALSE.
!
! For detrainment weight function
!
  real(r8) :: detrVel, detrWeight, detrSum
!
! Other local variables
!
  integer  :: KI, K, counter
  real(r8) :: rhoH, h01, h02
!
! ==================================================================!
!                                                                   !
! PART I - Find minimum detrainment depth                           !
!                                                                   !
! ==================================================================!
!
  detrDz = 0.0
  rhoP = PLUME(ng) % rho(I, osDepthK)
!
! First, check if plume detrains in surface layer.
!
  IF (plumeDepthK .EQ. N(ng)) isSurface = .TRUE.
!
  IF (isSurface) THEN
!
! Use Froude number to limit detrainment flow speed.
!
    rhoH = 0.0
    h1 = 0.0
    DO K = N(ng), 1, -1
      PLUME(ng) % detI(I, K) = 1
      detrDz = detrDz + PLUME(ng) % dz(I, K)
      h1 = h1+PLUME(ng) % dz(I, K)
      rhoH = rhoH+PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
      rho1 = rhoH/h1
      gRed = MAX(g*(rho1-rhoP)/rhoRef, 0.0)
      Fr = detr/(lc*h1)/SQRT(gRed*h1)
      IF (Fr .LE. 1.0/RiB) EXIT
    ENDDO
  ELSE  ! (.NOT. isSurface)
!
! Use the Richardson number criteria to determine if to distribute
! the detrainment in several layers. First, calculate BV frequency.
!
    PLUME(ng) % detI(I, plumeDepthK) = 1
    detrDz = detrDz + PLUME(ng) % dz(I, plumeDepthK)
!
! Get density of the upper and lower layer.
!
    rhoH = 0.0
    h1 = 0.0
    DO K = plumeDepthK, N(ng)
      IF (K .EQ. plumeDepthK) THEN
        h1 = h1+0.5*PLUME(ng) % dz(I, K)
        rhoH = rhoH+0.5*PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
      ELSE
        h1 = h1+PLUME(ng) % dz(I, K)
        rhoH = rhoH+PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
      ENDIF
    ENDDO
    rho1 = rhoH/h1
    rhoH = 0.0
    h2 = 0.0
    DO K = 1, plumeDepthK
      IF (K .EQ. plumeDepthK) THEN
        h2 = h2+0.5*PLUME(ng) % dz(I, K)
        rhoH = rhoH+0.5*PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
      ELSE
        h2 = h2+PLUME(ng) % dz(I, K)
        rhoH = rhoH+PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
      ENDIF
    ENDDO
    rho2 = rhoH/h2
!
! Compute minimum detrainment thickness.
!
    gRed = MAX(-g*(rho1-rho2)/rhoRef, gRedBkg)
    maxVel = SQRT(gRed*h1*h2/(h1+h2))/RiB
    minDz = detr/maxVel/lc
!
! Add layers around plumeDepthK until it reaches critical thickness.
!
    counter = 0
# ifdef ICEPLUME_DETRAIN_HALF
    udSwitch = 1
# else
    KI = plumeDepthK
    h01 = 0.5*PLUME(ng) % dz(I, KI)
    h02 = 0.5*PLUME(ng) % dz(I, KI)
    potE = g*PLUME(ng) % dz(I, KI)*(PLUME(ng) % rhoAm(I, KI)-rhoP)
    IF (potE .GT. 0.0) THEN
      udSwitch = 1
    ELSE
      udSwitch = -1
    ENDIF
# endif
    DO WHILE ( (detrDz .LT. minDz) .AND. (counter .LT. 100) )
      counter = counter + 1
      IF (udSwitch .EQ. 1) THEN
!
! Search one layer up
!
        DO K = 1, N(ng)
          IF ((PLUME(ng) % detI(I, K) .EQ. 0) .AND.   &
            & (PLUME(ng) % detI(I, K-1) .EQ. 1)) THEN
            KI = K
            h01 = h01 + PLUME(ng) % dz(I, KI)
          ENDIF
        ENDDO
      ELSEIF (udSwitch .EQ. -1) THEN
!
! Search one layer down
!
        DO K = 1, N(ng)-1
          IF ((PLUME(ng) % detI(I, K+1) .EQ. 1) .AND. &
            & (PLUME(ng) % detI(I, K) .EQ. 0)) THEN
            KI = K
            h02 = h02 + PLUME(ng) % dz(I, KI)
          ENDIF
        ENDDO
      ENDIF
# ifdef ICEPLUME_DETRAIN_HALF
!
! Search one layer up then one layer down.
!
      udSwitch = udSwitch*(-1)
# else
!
! Search layer to minimize potential energy anomaly.
!
      potE = potE + &
        & g*PLUME(ng) % dz(I, KI)*(PLUME(ng) % rhoAm(I, KI)-rhoP)
      IF (potE .GT. 0.0) THEN
        udSwitch = 1
      ELSE
        udSwitch = -1
      ENDIF
# endif
!
! Update detrainment layer flag and detrainment depth.
!
      PLUME(ng) % detI(I, KI) = 1
      detrDz = detrDz + PLUME(ng) % dz(I, KI)
!
! IF layer reached surface/bottom, force udSwitch to switch sign
!
      IF     ( (PLUME(ng) % detI(I, N(ng)) .EQ. 1) .AND.    &
        &      (PLUME(ng) % detI(I, 1)     .EQ. 0) ) THEN
        udSwitch = -1
      ELSEIF ( (PLUME(ng) % detI(I, 1)     .EQ. 1) .AND.    &
        &      (PLUME(ng) % detI(I, N(ng)) .EQ. 0)) THEN
        udSwitch = 1
      ELSEIF ( (PLUME(ng) % detI(I, 1)     .EQ. 1) .AND.    &
        &      (PLUME(ng) % detI(I, N(ng)) .EQ. 1)) THEN
        udSwitch = 0
      ENDIF
!
! This is another exit loop criteria.
!
      IF ( (detrDz .EQ. (detrDz+1.0)) .OR. (udSwitch .EQ. 0) ) EXIT
    ENDDO
  ENDIF
!
! ==================================================================!
!                                                                   !
! PART II - Distribute detrainmnet in several layers                !
!                                                                   !
! Update log                                                        !
! Use a Gause function to smooth the distribution.                  !
! 2018/06/08 Chuning Wang                                           !
!                                                                   !
! ==================================================================!
!
! Update detrainment volume flux
!
  detrVel = detr / (lc * detrDz)
  DO K = 1,N(ng)
    IF (PLUME(ng) % detI(I, K) .EQ. 1) THEN
!
! First, calculate weight using Gaussian function.
!
# ifdef ICEPLUME_DETRAIN_HALF
      detrWeight = EXP(-1.0 *                   &
        & ((PLUME(ng) % zR(I, K)-               &
        &   PLUME(ng) % zR(I, plumeDepthK)) /   &
        &  detrDz)**2)
# else
!       IF (K .GT. plumeDepthK) THEN
!         detrWeight = EXP(-0.5 *                 &
!           & ((PLUME(ng) % zR(I, K)-             &
!           &   PLUME(ng) % zR(I, plumeDepthK)) / &
!           &  h01)**2)
!       ELSE
!         detrWeight = EXP(-0.5 *                 &
!           & ((PLUME(ng) % zR(I, K)-             &
!           &   PLUME(ng) % zR(I, plumeDepthK)) / &
!           &  h02)**2)
!       ENDIF
      detrWeight = 1.0
# endif
!
! Distribute the detrainment.
!
      PLUME(ng) % det(I, K) = detrWeight * detrVel * lc * &
        & PLUME(ng) % dz(I, K)
    ENDIF
  ENDDO
!
! Normalize
!
  detrSum = SUM(PLUME(ng) % det(I, :))
  DO K = 1,N(ng)
    PLUME(ng) % det(I, K) = PLUME(ng) % det(I, K)*detr/detrSum
  ENDDO
END SUBROUTINE ICEPLUME_DETRAIN
# ifdef ICEPLUME_DETRAIN_NEUTRAL
!
! ==================================================================!
!                                                                   !
! Another method to calculate detrainment in each layer. It         !
! computes tracer concentration in each layer to make detrainment   !
! have similar density to ambient water.                            !
!                                                                   !
! ==================================================================!
!
SUBROUTINE ICEPLUME_DETRAIN2(ng, I,                           &
                          & iceDepthK, plumeDepthK, osDepthK, &
                          & lc, det)
  USE mod_iceplume
  implicit none
!
! ==================================================================!
!                                                                   !
! Input/Output variables:                                           !
!                                                                   !
! ==================================================================!
!                                                                   !
! ng          - grid identifier                                     !
! I           - point source index identifier                       !
! iceDepthK   - glacier grounding line depth vertical grid index    !
! plumeDepthK - plume detrainment depth vertical grid index         !
! osDepthK    - plume overshoot depth vertical grid index           !
! lc          - plume/ocean contact length [m]                      !
! det         - detrainment volume flux [m^3 s^-1]                  !
!                                                                   !
! ==================================================================!
!                                                                   !
! Local variables:                                                  !
!                                                                   !
! ==================================================================!
!                                                                   !
! detDz       - detrainment layer total thickness [m]               !
! minDz       - minimum detrainment layer thickness [m]             !
!               criterion [m s^-1]                                  !
! potE        - potential energy of plume/ambient water             !
! udSwitch    - switch to determine search one layer up(1)/down(-1) !
! isSurface   - logical, if the plume discharge at surface          !
! isBottom    - logical, if the plume discharge at bottom           !
! detVel      - detrainment velocity [m s^-1]                       !
!                                                                   !
! volR        - detrainment volume of runoff [m^3]                  !
! volM        - detrainment volume of meltwater                     !
! volF        - detrainment volume of total freshwater              !
! volE        - detrainment volume of entrainment                   !
! volP        - detrainment volume of plume                         !
!                                                                   !
! tF          - temperature of total freshwater [degC]              !
! tE          - temperature of entrainment                          !
! tP          - temperature of plume                                !
! tGade       - temperature of meltwater (Gade line)                !
!                                                                   !
! sF          - salinity of total freshwater [PSU]                  !
! sE          - salinity of entrainment                             !
! sP          - salinity of plume                                   !
!                                                                   !
! rhoF        - density of total freshwater [kg m^-3]               !
! rhoE        - density of entrainment                              !
! rhoP        - density of plume                                    !
!                                                                   !
! ==================================================================!
!
! In/out variables
!
  integer, intent(in) :: ng, I
  integer, intent(in) :: iceDepthK, plumeDepthK, osDepthK
  real(r8), intent(in) :: lc, det
!
! Local variables declaration
!
  real(r8) :: detDz, minDz, potE
  integer  :: udSwitch
  logical  :: isSurface = .FALSE.
  logical  :: isBottom = .FALSE.
  real(r8) :: detVel
!
  real(r8) :: volR, volM, volF, volE, volP
  real(r8) :: tF, tE, tP, tGade, sF, sE, sP, rhoF, rhoE, rhoP
  real(r8) :: RHO
  integer  :: KI, K1, K2, K, itrc, counter
  real(r8) :: cff, cff1, cff2
!
! Initialize.
!
  DO K = 1, N(ng)
    PLUME(ng) % detFrac(I, K) = 0.0
    PLUME(ng) % detF(I, K) = 0.0
    PLUME(ng) % detE(I, K) = 0.0
  ENDDO
  DO itrc = 1, NT(ng)
    DO K = 1, N(ng)
      PLUME(ng) % detTrc(I, K, itrc) = 0.0
    ENDDO
  ENDDO
!
! ==================================================================!
!                                                                   !
! PART I - Find minimum detrainment depth.                          !
!                                                                   !
! ==================================================================!
!
  IF (plumeDepthK .EQ. N(ng)) isSurface = .TRUE.
  IF (plumeDepthK .EQ. 1) isBottom = .TRUE.
  detDz = 0.0
  rhoP = PLUME(ng) % rho(I, osDepthK)
!
! Compute minimum detrainment thickness.
!
  minDz = det/lc
!
! Add layers around plumeDepthK until it reaches critical thickness.
!
  counter = 0
  KI = plumeDepthK
  PLUME(ng) % detI(I, KI) = 1
  detDz = PLUME(ng) % dz(I, KI)
  PLUME(ng) % detFrac(I, KI) = 1.0
  K1 = KI
  K2 = KI
  potE = g*PLUME(ng) % dz(I, KI)*(PLUME(ng) % rhoAm(I, KI)-rhoP)
  IF (potE .GT. 0.0) THEN
    udSwitch = 1
  ELSE
    udSwitch = -1
  ENDIF
  DO WHILE ( (detDz .LT. minDz) .AND. (counter .LT. 100) )
    counter = counter + 1
    IF (udSwitch .EQ. 1) THEN
!
! Search one layer up
!
      DO K = 1, N(ng)
        IF ((PLUME(ng) % detI(I, K) .EQ. 0) .AND.   &
          & (PLUME(ng) % detI(I, K-1) .EQ. 1)) THEN
          KI = K
        ENDIF
      ENDDO
      K1 = KI
    ELSEIF (udSwitch .EQ. -1) THEN
!
! Search one layer down
!
      DO K = 1, N(ng)-1
        IF ((PLUME(ng) % detI(I, K+1) .EQ. 1) .AND. &
          & (PLUME(ng) % detI(I, K) .EQ. 0)) THEN
          KI = K
        ENDIF
      ENDDO
      K2 = KI
    ENDIF
    IF (KI .EQ. N(ng)) isSurface = .TRUE.
    IF (KI .EQ. 1) isBottom = .TRUE.
!
! Search layer to minimize potential energy anomaly.
!
    potE = potE + &
      & g*PLUME(ng) % dz(I, KI)*(PLUME(ng) % rhoAm(I, KI)-rhoP)
    IF (potE .GT. 0.0) THEN
      udSwitch = 1
    ELSE
      udSwitch = -1
    ENDIF
!
! Update detrainment layer flag and detrainment depth.
!
    PLUME(ng) % detI(I, KI) = 1
    detDz = detDz + PLUME(ng) % dz(I, KI)
    PLUME(ng) % detFrac(I, KI) = 1.0
!
! IF layer reached surface/bottom, force udSwitch to switch sign
!
    IF     ( (       isSurface ) .AND. ( .NOT. isBottom ) ) THEN
      udSwitch = -1
    ELSEIF ( ( .NOT. isSurface ) .AND. (       isBottom ) ) THEN
      udSwitch = 1
    ELSEIF ( (       isSurface ) .AND. (       isBottom ) ) THEN
      udSwitch = 0
    ENDIF
!
! Other exit loop criteria.
!
    IF ( (detDz .EQ. (detDz+1.0)) .OR. (udSwitch .EQ. 0) ) EXIT
  ENDDO
!
! Only take partial grid from top/bottom layer to satisfy potE=0
!
  IF ( (.NOT. isSurface) .AND. (counter .GT. 0) ) THEN
    IF ( potE .GT. 0.0 ) THEN
      KI = K2
    ELSEIF ( potE .LT. 0.0 ) THEN
      KI = K1
    ENDIF
      cff = PLUME(ng) % dz(I, KI) - &
        & potE/(g*(PLUME(ng) % rhoAm(I, KI)-rhoP))
      PLUME(ng) % detFrac(I, KI) = cff / PLUME(ng) % dz(I, KI)
      detDz = detDz - (PLUME(ng) % dz(I, KI) - cff)
  ENDIF
!
! ==================================================================!
!                                                                   !
! PART II - Distribute detrainmnet in several layers.               !
!                                                                   !
! ==================================================================!
!
! Update detrainment volume flux
!
  detVel = det / (lc * detDz)
!
! Distribute the detrainment.
!
  DO K = 1, N(ng)
    PLUME(ng) % det(I, K) = detVel * lc *                 &
      & PLUME(ng) % detFrac(I, K) * PLUME(ng) % dz(I, K)
  ENDDO
!
! Calculate proportion of each water mass.
! There are three main water masses -
! 1. Runoff (R),
! 2. Melt (M),
! 3. Entrainment (E).
! R and M combined forms a new water mass,
! 4. Freshwater (F).
! M and F combined froms the final water mass,
! 5. Plume (P).
!
  volP = PLUME(ng) % f(I, osDepthK)
  volR = PLUME(ng) % f(I, iceDepthK)
  volM = PLUME(ng) % m(I, osDepthK)
  volF = PLUME(ng) % m(I, osDepthK) + PLUME(ng) % f(I, iceDepthK)
  volE = volP - volF
!
! Temperature. Note that tGade is effective temperature of Gade line.
!
  tP = PLUME(ng) % trc(I, itemp)
  tGade = (L - cI*tIce)/cW
  tF = (tGade*volM + PLUME(ng) % t(I, iceDepthK)*volR)/volF
!
! Salinity.
!
  sP = PLUME(ng) % trc(I, isalt)
  sF = (sIce*volM  + PLUME(ng) % s(I, iceDepthK)*volR)/volF
!
! Density.
!
  rhoF = RHO(tF, sF, PLUME(ng) % zR(I, plumeDepthK))
  rhoE = (rhoP*volP-rhoF*volF)/volE
!
! Calculate proportion from mixing line.
!
  DO K = 1, N(ng)
    IF (PLUME(ng) % detI(I, K) .EQ. 1) THEN
      cff1 = MIN(MAX((PLUME(ng) % rhoAm(I, K)-rhoE)/(rhoF-rhoE), &
        &            0.0), 1.0)
      cff2 = 1-cff1
      PLUME(ng) % detF(I, K) = PLUME(ng) % det(I, K)*cff1
      PLUME(ng) % detE(I, K) = PLUME(ng) % det(I, K)*cff2
    ENDIF
  ENDDO
!
! Make corrections to make sure each water mass sums to its previous
! volume.
!
  cff1 = volF - SUM(PLUME(ng) % detF(I, :))
  cff2 = volE - SUM(PLUME(ng) % detE(I, :))
  DO K = 1, N(ng)
    IF (PLUME(ng) % detI(I, K) .EQ. 1) THEN
      cff = PLUME(ng) % dz(I, K)*PLUME(ng) % detFrac(I, K)/detDz
      PLUME(ng) % detF(I, K) = PLUME(ng) % detF(I, K) + cff*cff1
      PLUME(ng) % detE(I, K) = PLUME(ng) % detE(I, K) + cff*cff2
    ENDIF
  ENDDO
!
! ==================================================================!
!                                                                   !
! Update tracer concentration in detrainment model.                 !
!                                                                   !
! ==================================================================!
!
! Active tracers.
!
  tE = (tP*volP-tF*volF)/volE
  sE = (sP*volP-sF*volF)/volE
  DO K = 1, N(ng)
    IF (PLUME(ng) % detI(I, K) .EQ. 1) THEN
      PLUME(ng) % detTrc(I, K, itemp) =                           &
        & (tF*PLUME(ng) % detF(I, K)+tE*PLUME(ng) % detE(I, K)) / &
        & (PLUME(ng) % detF(I, K)+PLUME(ng) % detE(I, K))
      PLUME(ng) % detTrc(I, K, isalt) =                           &
        & (sF*PLUME(ng) % detF(I, K)+sE*PLUME(ng) % detE(I, K)) / &
        & (PLUME(ng) % detF(I, K)+PLUME(ng) % detE(I, K))
    ENDIF
  ENDDO
!
! Passive tracers.
!
  DO itrc = 3, NT(ng)
    cff1 = PLUME(ng) % trcIni(I, itrc)*volR/volF
#  ifdef ICEPLUME_MELT_TRACER
    IF ( itrc .EQ. NT(ng)-1 ) THEN
      cff1 = PLUME(ng) % trcIni(I, itrc)*volM/volF
    ELSEIF ( itrc .EQ. NT(ng) ) THEN
      cff1 = 0.0
    ENDIF
#  endif
    cff2 = (PLUME(ng) % trc(I, itrc)*volP - cff1*volF)/volE
    DO K = 1, N(ng)
      IF (PLUME(ng) % detI(I, K) .EQ. 1) THEN
        PLUME(ng) % detTrc(I, K, itrc) =    &
          & (cff1*PLUME(ng) % detF(I, K) +  &
          &  cff2*PLUME(ng) % detE(I, K)) / &
          & (PLUME(ng) % detF(I, K) +       &
          &  PLUME(ng) % detE(I, K))
      ENDIF
    ENDDO
  ENDDO
END SUBROUTINE ICEPLUME_DETRAIN2
!
SUBROUTINE ICEPLUME_DETRAIN3(ng, I,                           &
                          & iceDepthK, plumeDepthK, osDepthK, &
                          & lc, det)
  USE mod_iceplume
  implicit none
!
! ==================================================================!
!                                                                   !
! Input/Output variables:                                           !
!                                                                   !
! ==================================================================!
!                                                                   !
! ng          - grid identifier                                     !
! I           - point source index identifier                       !
! iceDepthK   - glacier grounding line depth vertical grid index    !
! plumeDepthK - plume detrainment depth vertical grid index         !
! osDepthK    - plume overshoot depth vertical grid index           !
! lc          - plume/ocean contact length [m]                      !
! det         - detrainment volume flux [m^3 s^-1]                  !
!                                                                   !
! ==================================================================!
!                                                                   !
! Local variables:                                                  !
!                                                                   !
! ==================================================================!
!                                                                   !
! detDz       - detrainment layer total thickness [m]               !
! minDz       - minimum detrainment layer thickness [m]             !
!               criterion [m s^-1]                                  !
! potE        - potential energy of plume/ambient water             !
! udSwitch    - switch to determine search one layer up(1)/down(-1) !
! isSurface   - logical, if the plume discharge at surface          !
! isBottom    - logical, if the plume discharge at bottom           !
! detVel      - detrainment velocity [m s^-1]                       !
!                                                                   !
! volR        - detrainment volume of runoff [m^3]                  !
! volM        - detrainment volume of meltwater                     !
! volF        - detrainment volume of total freshwater              !
! volE        - detrainment volume of entrainment                   !
! volP        - detrainment volume of plume                         !
!                                                                   !
! tF          - temperature of total freshwater [degC]              !
! tE          - temperature of entrainment                          !
! tP          - temperature of plume                                !
! tGade       - temperature of meltwater (Gade line)                !
!                                                                   !
! sF          - salinity of total freshwater [PSU]                  !
! sE          - salinity of entrainment                             !
! sP          - salinity of plume                                   !
!                                                                   !
! rhoF        - density of total freshwater [kg m^-3]               !
! rhoE        - density of entrainment                              !
! rhoP        - density of plume                                    !
!                                                                   !
! ==================================================================!
!
! In/out variables
!
  integer, intent(in) :: ng, I
  integer, intent(in) :: iceDepthK, plumeDepthK, osDepthK
  real(r8), intent(in) :: lc, det
!
! Local variables declaration
!
  real(r8) :: detDz, minDz, potE
  integer  :: udSwitch
  logical  :: isSurface = .FALSE.
  logical  :: isBottom = .FALSE.
  real(r8) :: detVel
!
  real(r8) :: volR, volM, volF, volE, volP
  real(r8) :: tF, tE, tP, tGade, sF, sE, sP, rhoF, rhoE, rhoP
  real(r8) :: RHO
  integer  :: KI, K1, K2, K, itrc, counter
  real(r8) :: cff, cff1, cff2
  real(r8) :: rhoH, rho1, rho2, h1, h2, gRed, maxVel, Fr
!
! Initialize.
!
  DO K = 1, N(ng)
    PLUME(ng) % detFrac(I, K) = 0.0
    PLUME(ng) % detF(I, K) = 0.0
    PLUME(ng) % detE(I, K) = 0.0
  ENDDO
  DO itrc = 1, NT(ng)
    DO K = 1, N(ng)
      PLUME(ng) % detTrc(I, K, itrc) = 0.0
    ENDDO
  ENDDO
!
! ==================================================================!
!                                                                   !
! PART I - Find minimum detrainment depth.                          !
!                                                                   !
! ==================================================================!
!
  IF (plumeDepthK .EQ. N(ng)) isSurface = .TRUE.
  IF (plumeDepthK .EQ. 1) isBottom = .TRUE.
  detDz = 0.0
  rhoP = PLUME(ng) % rho(I, osDepthK)
!
! Compute minimum detrainment thickness.
!
!   minDz = det/lc
!
! Get density of the upper and lower layer.
!
  rhoH = 0.0
  h1 = 0.0
  DO K = plumeDepthK, N(ng)
    IF (K .EQ. plumeDepthK) THEN
      h1 = h1+0.5*PLUME(ng) % dz(I, K)
      rhoH = rhoH+0.5*PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
    ELSE
      h1 = h1+PLUME(ng) % dz(I, K)
      rhoH = rhoH+PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
    ENDIF
  ENDDO
  rho1 = rhoH/h1
  rhoH = 0.0
  h2 = 0.0
  DO K = 1, plumeDepthK
    IF (K .EQ. plumeDepthK) THEN
      h2 = h2+0.5*PLUME(ng) % dz(I, K)
      rhoH = rhoH+0.5*PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
    ELSE
      h2 = h2+PLUME(ng) % dz(I, K)
      rhoH = rhoH+PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
    ENDIF
  ENDDO
  rho2 = rhoH/h2
!
! Compute minimum detrainment thickness.
!
  gRed = MAX(-g*(rho1-rho2)/rhoRef, gRedBkg)
  maxVel = SQRT(gRed*h1*h2/(h1+h2))/RiB
  minDz = det/maxVel/lc
!
! Add layers around plumeDepthK until it reaches critical thickness.
!
  counter = 0
  KI = plumeDepthK
  PLUME(ng) % detI(I, KI) = 1
  detDz = PLUME(ng) % dz(I, KI)
  PLUME(ng) % detFrac(I, KI) = 1.0
  K1 = KI
  K2 = KI
  potE = g*PLUME(ng) % dz(I, KI)*(PLUME(ng) % rhoAm(I, KI)-rhoP)
  IF (potE .GT. 0.0) THEN
    udSwitch = 1
  ELSE
    udSwitch = -1
  ENDIF
  DO WHILE ( (detDz .LT. minDz) .AND. (counter .LT. 100) )
    counter = counter + 1
    IF (udSwitch .EQ. 1) THEN
!
! Search one layer up
!
      DO K = 1, N(ng)
        IF ((PLUME(ng) % detI(I, K) .EQ. 0) .AND.   &
          & (PLUME(ng) % detI(I, K-1) .EQ. 1)) THEN
          KI = K
        ENDIF
      ENDDO
      K1 = KI
    ELSEIF (udSwitch .EQ. -1) THEN
!
! Search one layer down
!
      DO K = 1, N(ng)-1
        IF ((PLUME(ng) % detI(I, K+1) .EQ. 1) .AND. &
          & (PLUME(ng) % detI(I, K) .EQ. 0)) THEN
          KI = K
        ENDIF
      ENDDO
      K2 = KI
    ENDIF
    IF (KI .EQ. N(ng)) isSurface = .TRUE.
    IF (KI .EQ. 1) isBottom = .TRUE.
!
! Search layer to minimize potential energy anomaly.
!
    potE = potE + &
      & g*PLUME(ng) % dz(I, KI)*(PLUME(ng) % rhoAm(I, KI)-rhoP)
    IF (potE .GT. 0.0) THEN
      udSwitch = 1
    ELSE
      udSwitch = -1
    ENDIF
!
! Update detrainment layer flag and detrainment depth.
!
    PLUME(ng) % detI(I, KI) = 1
    detDz = detDz + PLUME(ng) % dz(I, KI)
    PLUME(ng) % detFrac(I, KI) = 1.0
!
! IF layer reached surface/bottom, force udSwitch to switch sign
!
    IF     ( (       isSurface ) .AND. ( .NOT. isBottom ) ) THEN
      udSwitch = -1
    ELSEIF ( ( .NOT. isSurface ) .AND. (       isBottom ) ) THEN
      udSwitch = 1
    ELSEIF ( (       isSurface ) .AND. (       isBottom ) ) THEN
      udSwitch = 0
    ENDIF
!
! Other exit loop criteria.
!
    IF ( (detDz .EQ. (detDz+1.0)) .OR. (udSwitch .EQ. 0) ) EXIT
  ENDDO
!
! Only take partial grid from top/bottom layer to satisfy potE=0
!
  IF ( (.NOT. isSurface) .AND. (counter .GT. 0) ) THEN
    IF ( potE .GT. 0.0 ) THEN
      KI = K2
    ELSEIF ( potE .LT. 0.0 ) THEN
      KI = K1
    ENDIF
      cff = PLUME(ng) % dz(I, KI) - &
        & potE/(g*(PLUME(ng) % rhoAm(I, KI)-rhoP))
      PLUME(ng) % detFrac(I, KI) = cff / PLUME(ng) % dz(I, KI)
      detDz = detDz - (PLUME(ng) % dz(I, KI) - cff)
  ENDIF
  IF (isSurface) THEN
!
! If the plume reach surface, use Froude number to limit detrainment
! flow speed, rewrite detrainment.
!
    DO K = 1, N(ng)
      PLUME(ng) % detI(I, K) = 0
      PLUME(ng) % detFrac(I, K) = 0.0
    ENDDO
    detDz = 0.0
    rhoH = 0.0
    h1 = 0.0
    DO K = N(ng), 1, -1
      PLUME(ng) % detI(I, K) = 1
      PLUME(ng) % detFrac(I, K) = 1.0
      detDz = detDz + PLUME(ng) % dz(I, K)
      h1 = h1+PLUME(ng) % dz(I, K)
      rhoH = rhoH+PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
      rho1 = rhoH/h1
      gRed = MAX(g*(rho1-rhoP)/rhoRef, gRedBkg)
      Fr = det/(lc*h1)/SQRT(gRed*h1)
      IF (Fr .LE. 1.0/RiB) EXIT
    ENDDO
  ENDIF
!
! ==================================================================!
!                                                                   !
! PART II - Distribute detrainmnet in several layers.               !
!                                                                   !
! ==================================================================!
!
! Update detrainment volume flux
!
  detVel = det / (lc * detDz)
!
! Distribute the detrainment.
!
  DO K = 1, N(ng)
    PLUME(ng) % det(I, K) = detVel * lc *                 &
      & PLUME(ng) % detFrac(I, K) * PLUME(ng) % dz(I, K)
  ENDDO
!
! Calculate proportion of each water mass.
! There are three main water masses -
! 1. Runoff (R),
! 2. Melt (M),
! 3. Entrainment (E).
! R and M combined forms a new water mass,
! 4. Freshwater (F).
! M and F combined froms the final water mass,
! 5. Plume (P).
!
  volP = PLUME(ng) % f(I, osDepthK)
  volR = PLUME(ng) % f(I, iceDepthK)
  volM = PLUME(ng) % m(I, osDepthK)
  volF = PLUME(ng) % m(I, osDepthK) + PLUME(ng) % f(I, iceDepthK)
  volE = volP - volF
!
! Temperature. Note that tGade is effective temperature of Gade line.
!
  tP = PLUME(ng) % trc(I, itemp)
  tGade = (L - cI*tIce)/cW
  tF = (tGade*volM + PLUME(ng) % t(I, iceDepthK)*volR)/volF
!
! Salinity.
!
  sP = PLUME(ng) % trc(I, isalt)
  sF = (sIce*volM  + PLUME(ng) % s(I, iceDepthK)*volR)/volF
!
! Density.
!
  rhoF = RHO(tF, sF, PLUME(ng) % zR(I, plumeDepthK))
  rhoE = (rhoP*volP-rhoF*volF)/volE
!
! Calculate proportion from mixing line.
!
  DO K = 1, N(ng)
    IF (PLUME(ng) % detI(I, K) .EQ. 1) THEN
      cff1 = MIN(MAX((PLUME(ng) % rhoAm(I, K)-rhoE)/(rhoF-rhoE), &
        &            0.0), 1.0)
      cff2 = 1-cff1
      PLUME(ng) % detF(I, K) = PLUME(ng) % det(I, K)*cff1
      PLUME(ng) % detE(I, K) = PLUME(ng) % det(I, K)*cff2
    ENDIF
  ENDDO
!
! Make corrections to make sure each water mass sums to its previous
! volume.
!
  cff1 = volF - SUM(PLUME(ng) % detF(I, :))
  cff2 = volE - SUM(PLUME(ng) % detE(I, :))
  DO K = 1, N(ng)
    IF (PLUME(ng) % detI(I, K) .EQ. 1) THEN
      cff = PLUME(ng) % dz(I, K)*PLUME(ng) % detFrac(I, K)/detDz
      PLUME(ng) % detF(I, K) = PLUME(ng) % detF(I, K) + cff*cff1
      PLUME(ng) % detE(I, K) = PLUME(ng) % detE(I, K) + cff*cff2
    ENDIF
  ENDDO
!
! ==================================================================!
!                                                                   !
! Update tracer concentration in detrainment model.                 !
!                                                                   !
! ==================================================================!
!
! Active tracers.
!
  tE = (tP*volP-tF*volF)/volE
  sE = (sP*volP-sF*volF)/volE
  DO K = 1, N(ng)
    IF (PLUME(ng) % detI(I, K) .EQ. 1) THEN
      PLUME(ng) % detTrc(I, K, itemp) =                           &
        & (tF*PLUME(ng) % detF(I, K)+tE*PLUME(ng) % detE(I, K)) / &
        & (PLUME(ng) % detF(I, K)+PLUME(ng) % detE(I, K))
      PLUME(ng) % detTrc(I, K, isalt) =                           &
        & (sF*PLUME(ng) % detF(I, K)+sE*PLUME(ng) % detE(I, K)) / &
        & (PLUME(ng) % detF(I, K)+PLUME(ng) % detE(I, K))
    ENDIF
  ENDDO
!
! Passive tracers.
!
  DO itrc = 3, NT(ng)
    cff1 = PLUME(ng) % trcIni(I, itrc)*volR/volF
#  ifdef ICEPLUME_MELT_TRACER
    IF ( itrc .EQ. NT(ng)-1 ) THEN
      cff1 = PLUME(ng) % trcIni(I, itrc)*volM/volF
    ELSEIF ( itrc .EQ. NT(ng) ) THEN
      cff1 = 0.0
    ENDIF
#  endif
    cff2 = (PLUME(ng) % trc(I, itrc)*volP - cff1*volF)/volE
    DO K = 1, N(ng)
      IF (PLUME(ng) % detI(I, K) .EQ. 1) THEN
        PLUME(ng) % detTrc(I, K, itrc) =    &
          & (cff1*PLUME(ng) % detF(I, K) +  &
          &  cff2*PLUME(ng) % detE(I, K)) / &
          & (PLUME(ng) % detF(I, K) +       &
          &  PLUME(ng) % detE(I, K))
      ENDIF
    ENDDO
  ENDDO
END SUBROUTINE ICEPLUME_DETRAIN3
# endif
#endif
