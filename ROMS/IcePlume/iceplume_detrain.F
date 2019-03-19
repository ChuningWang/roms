#include "cppdefs.h"
#ifdef ICEPLUME
!
! ==================================================================!
!                                                                   !
! Use this function to calculate the velocity of detrainment.       !
!                                                                   !
! ==================================================================!
!
SUBROUTINE ICEPLUME_DETRAIN_HALF(ng, I,                   &
                              &  iceDepthK, plumeDepthK,  &
                              &  osDepthK,                &
                              &  lc, detr, detrDz)
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
! osDepthK    - maximum overshooting depth vertical grid index      !
! lc          - plume/ocean contact length [m]                      !
! detr        - detrainment volume flux [m^3 s^-1]                  !
! detrDz      - initial detrainment layer thickness [m]             !
!                                                                   !
! ==================================================================!
!                                                                   !
! Local variables:                                                  !
!                                                                   !
! ==================================================================!
!                                                                   !
! maxVel      - maximum velocity based on Richardson number         !
!               criterion [m s^-1]                                  !
! meanVel     - ice tangental velocity [m s^-1]                     !
! Fr          - Froude number                                       !
! minDz       - minimum detrainment thickness [m]                   !
! rhoP        - plume density [kg m^-3]                             !
! rho1        - upper layer density [kg m^-3]                       !
! rho2        - lower layer density [kg m^-3]                       !
! h1          - upper layer thickness [m]                           !
! h2          - lower layer thickness [m]                           !
! gRed        - reduced gravity [m s^-2]                            !
!                                                                   !
! ==================================================================!
!
! In/out variables
!
  integer, intent(in) :: ng, I
  integer, intent(in) :: iceDepthK, plumeDepthK, osDepthK
  real(r8), intent(in) :: lc, detr
  real(r8), intent(out) :: detrDz
!
! Local variables declaration
!
  real(r8) :: maxVel, meanVel, minDz, Fr
  real(r8) :: rhoP, rho1, rho2, h1, h2, gRed
  logical :: isSurface = .FALSE.
!
! For detrainment weight function
!
  real(r8) :: detrVel, detrWeight, detrWeightSum
  integer  :: K1, K2, KI, searchSwitch
  integer :: K, counter
  real(r8) :: cpp
!
! ==================================================================!
!                                                                   !
! PART I - Find minimum detrainment depth                           !
!                                                                   !
! ==================================================================!
!
  detrDz = 0.0
!
! First, check if plume detrains in surface layer.
!
  IF (plumeDepthK .EQ. N(ng)) isSurface = .TRUE.
!
  IF (isSurface) THEN
!
! Use Froude number to limit detrainment flow speed.
!
    cpp = 0.0
    h1 = 0.0
    rhoP = PLUME(ng) % rho(I, plumeDepthK)
    DO K = N(ng):1:-1
      PLUME(ng) % detI(I, K) = 1
      detrDz = detrDz + PLUME(ng) % dz(I, K)
      h1 = h1+PLUME(ng) % dz(I, K)
      cpp = cpp+PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
      rho1 = cpp/h1
      gRed = MAX(g*(rho1-rhoP)/rhoRef, 0.0)
      Fr = detr/(lc*h1)/SQRT(gRed*h1)
      IF (Fr .LE. 1.0) EXIT
    ENDDO
  ELSE  ! (.NOT. isSurface)
!
! Check the Richardson number criteria to determine if to distribute
! the detrainment in several layers. First, calculate BV frequency
!
    PLUME(ng) % detI(I, plumeDepthK) = 1
    detrDz = detrDz + PLUME(ng) % dz(I, plumeDepthK)
!
    cpp = 0.0
    h1 = 0.0
    DO K = plumeDepthK:N(ng)
      IF (K .EQ. plumeDpethK) THEN
        h1 = h1+0.5*PLUME(ng) % dz(I, K)
        cpp = cpp+0.5*PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
      ELSE
        h1 = h1+PLUME(ng) % dz(I, K)
        cpp = cpp+PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
      ENDIF
    ENDDO
    rho1 = cpp/h1
!
    cpp = 0.0
    h2 = 0.0
    DO K = 1:plumeDepthK
      IF (K .EQ. plumeDpethK) THEN
        h2 = h2+0.5*PLUME(ng) % dz(I, K)
        cpp = cpp+0.5*PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
      ELSE
        h2 = h2+PLUME(ng) % dz(I, K)
        cpp = cpp+PLUME(ng) % dz(I, K)*PLUME(ng) % rhoAm(I, K)
      ENDIF
    ENDDO
    rho2 = cpp/h2
!
    gRed = MAX(-g*(rho1-rho2)/rhoRef, 0.0)
    maxVel = RiBMin*MAX(SQRT(gRed*h1*h2/(h1+h2)), 1.0E-5)
    minDz = detr/maxVel/lc
    searchSwitch = 1
    counter = 0
!
    DO WHILE ( (detrDz .LT. minDz) .AND. (counter .LT. 100) )
      counter = counter + 1
!
! Search one layer up then one layer down. If the plume has reached
! surface/bottom, set searchswitch to -1/1
!
      IF     ( (PLUME(ng) % detI(I, N(ng)) .EQ. 1) .AND.  &
        &      (PLUME(ng) % detI(I, 0) .EQ. 1) ) THEN
        searchSwitch = 0
        EXIT
      ELSEIF ( (PLUME(ng) % detI(I, N(ng)) .EQ. 1) .AND.  &
        &      (PLUME(ng) % detI(I, 0) .EQ. 0) ) THEN
        searchSwitch = -1
      ELSEIF ( (PLUME(ng) % detI(I, N(ng)) .EQ. 0) .AND.  &
        &      (PLUME(ng) % detI(I, 0) .EQ. 1) ) THEN
        searchSwitch = 1
      ENDIF
!
      IF (searchSwitch .EQ. 1) THEN
!
! Search one layer up
!
        DO K = 1,N(ng)
          IF ((PLUME(ng) % detI(I, K) .EQ. 0) .AND. &
            & (PLUME(ng) % detI(I, K-1) .EQ. 1)) THEN
            KI = K
          ENDIF
        ENDDO
      ELSEIF (searchSwitch .EQ. -1) THEN
!
! Search one layer down
!
        DO K = 1,N(ng)-1
          IF ((PLUME(ng) % detI(I, K+1) .EQ. 1) .AND. &
            & (PLUME(ng) % detI(I, K) .EQ. 0)) THEN
            KI = K
          ENDIF
        ENDDO
      ENDIF
      searchSwitch = searchSwitch*(-1)
!
! Update detrainment layer flag and detrainment depth
!
      PLUME(ng) % detI(I, KI) = 1
      detrDz = detrDz + PLUME(ng) % dz(I, KI)
!
! Another exit loop criteria
!
      IF (detrDz .EQ. (detrDz+1.0)) EXIT
    ENDDO
  ENDIF
!
!
! ==================================================================!
!                                                                   !
! PART II - Distribute detrainmnet in several layers                !
!                                                                   !
! Update log                                                        !
! Use a Gause function to smooth the distribution                   !
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
      detrWeight = EXP(-0.5 *                 &
        & ((PLUME(ng) % zR(I, K)-             &
        &   PLUME(ng) % zR(I, plumeDepthK)) / &
        &  detrDz)**2)
      PLUME(ng) % det(I, K) = detrWeight * detrVel * lc * &
        & PLUME(ng) % dz(I, K)
    ENDIF
  ENDDO
!
! Normalize
!
  detrWeightSum = SUM(PLUME(ng) % det(I, :))
  DO K = 1,N(ng)
    PLUME(ng) % det(I, K) = PLUME(ng) % det(I, K) * &
      & detr / detrWeightSum
  ENDDO
END SUBROUTINE ICEPLUME_DETRAIN_HALF
#endif
