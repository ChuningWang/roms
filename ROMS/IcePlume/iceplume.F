#include "cppdefs.h"
MODULE iceplume_mod
#ifdef ICEPLUME
!
!=======================================================================
!                                                                      !
!  This subroutine is a wraper for the ICEPLUME module.                !
!                                                                      !
!=======================================================================
!
implicit none
!
PRIVATE
PUBLIC :: iceplume
!
CONTAINS
!
!***********************************************************************
SUBROUTINE iceplume (ng, tile)
!***********************************************************************
!
  USE mod_grid
  USE mod_ocean
  USE mod_stepping
!
! Imported variable declarations.
!
  integer, intent(in) :: ng, tile
!
! Local variable declarations.
!
# include "tile.h"
!
  CALL iceplume_tile (ng, tile,                                    &
    &                 LBi, UBi, LBj, UBj,                          &
    &                 IminS, ImaxS, JminS, JmaxS,                  &
    &                 nrhs(ng),                                    &
    &                 GRID(ng) % rmask,                            &
    &                 GRID(ng) % om_r,                             &
    &                 GRID(ng) % on_r,                             &
    &                 GRID(ng) % h,                                &
    &                 GRID(ng) % z_w,                              &
    &                 OCEAN(ng) % u,                               &
    &                 OCEAN(ng) % v,                               &
    &                 OCEAN(ng) % wvel,                            &
    &                 OCEAN(ng) % t)
!
  RETURN
END SUBROUTINE iceplume
!
!***********************************************************************
SUBROUTINE iceplume_tile(ng, tile,                                 &
    &                    LBi, UBi, LBj, UBj,                       &
    &                    IminS, ImaxS, JminS, JmaxS,               &
    &                    nrhs,                                     &
    &                    rmask,                                    &
    &                    om_r,                                     &
    &                    on_r,                                     &
    &                    h,                                        &
    &                    z_w,                                      &
    &                    u,                                        &
    &                    v,                                        &
    &                    wvel,                                     &
    &                    t)
!***********************************************************************
!
! Load public modules
!
  USE mod_scalars, ONLY : dt, LuvSrc
  USE mod_sources
!
! Load private module
!
  USE mod_iceplume
!
!  Imported variable declarations.
!
  integer, intent(in) :: ng, tile
  integer, intent(in) :: LBi, UBi, LBj, UBj
  integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
  integer, intent(in) :: nrhs
!
# ifdef ASSUMED_SHAPE
  real(r8), intent(in) :: rmask(LBi:,LBj:)
  real(r8), intent(in) :: om_r(LBi:,LBj:)
  real(r8), intent(in) :: on_r(LBi:,LBj:)
  real(r8), intent(in) :: h(LBi:,LBj:)
  real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
  real(r8), intent(in) :: u(LBi:,LBj:,:,:)
  real(r8), intent(in) :: v(LBi:,LBj:,:,:)
  real(r8), intent(in) :: wvel(LBi:,LBj:,0:)
  real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
# else
  real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
  real(r8), intent(in) :: om_r(LBi:UBi,LBj:UBj)
  real(r8), intent(in) :: on_r(LBi:UBi,LBj:UBj)
  real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
  real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
  real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
  real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
  real(r8), intent(in) :: wvel(LBi:UBi,LBj:UBj,0:N(ng))
  real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
# endif
!
! Local variable declarations.
!
  integer :: is, Isrc, Jsrc, Dsrc, D2src
  integer :: K, iTracer
  real(r8) :: pr, prRef
  real(r8) :: dx, dy
  real(r8) :: fIni, tIni, sIni
  real(r8) :: Qsrc
  real(r8) :: cpp1, cpp2
!
# include "set_bounds.h"
!
! ==================================================================
!
  IF (LuvSrc(ng)) THEN
    DO is=1,Nsrc(ng)
      Isrc = SOURCES(ng) % Isrc(is)
      Jsrc = SOURCES(ng) % Jsrc(is)
      IF ( ((IstrR .LE. Isrc) .AND. (Isrc .LE. IendR)) .AND. &
        &  ((JstrR .LE. Jsrc) .AND. (Jsrc .LE. JendR)) ) THEN
!
! ==================================================================!
! Find the right Rho coordinates                                    !
! ==================================================================!
!
        Dsrc = INT(SOURCES(ng) % Dsrc(is))
        IF (Dsrc .EQ. 0) THEN
          IF ( (rmask(Isrc  , Jsrc) .EQ. 0) .AND. &
            &  (rmask(Isrc-1, Jsrc) .EQ. 1) ) THEN
            Isrc = Isrc-1
            D2src = -1
          ELSEIF ( (rmask(Isrc  , Jsrc) .EQ. 1) .AND. &
                &  (rmask(Isrc-1, Jsrc) .EQ. 0) ) THEN
            D2src = 1
          ELSE
            D2src = 0
          ENDIF
        ELSEIF (Dsrc .EQ. 1) THEN
          IF ((rmask(Isrc, Jsrc  ) .EQ. 0) .AND. &
            & (rmask(Isrc, Jsrc-1) .EQ. 1)) THEN
            Jsrc = Jsrc-1
            D2src = -1
          ELSEIF ((rmask(Isrc, Jsrc  ) .EQ. 0) .AND. &
                & (rmask(Isrc, Jsrc-1) .EQ. 1)) THEN
            D2src = 1
          ELSE
            D2src = 0
          ENDIF
        ENDIF
!
! Get grid dimensions
!
        IF (Dsrc .EQ. 0) THEN
          dx = ABS(om_r(Isrc, Jsrc))
          dy = ABS(on_r(Isrc, Jsrc))
        ELSE
          dx = ABS(on_r(Isrc, Jsrc))
          dy = ABS(om_r(Isrc, Jsrc))
        ENDIF
!
! Get initial condition for iceplume model
!
        fIni = ABS(SOURCES(ng) % SGtra(is))
        tIni = &
# ifdef ONE_TRACER_SOURCE
          & SOURCES(ng) % SGtrc(itemp)
# else
          & SOURCES(ng) % SGtrc(is, itemp)
# endif
        sIni = &
# ifdef ONE_TRACER_SOURCE
          & SOURCES(ng) % SGtrc(isalt)
# else
          & SOURCES(ng) % SGtsrc(is, isalt)
# endif
!
! ==================================================================!
! Read in profiles from OCEAN.                                      !
! ==================================================================!
!
        DO K = 0, N(ng)
          PLUME(ng) % zW(is, K) = z_w(Isrc, Jsrc, K)
        ENDDO
!
        DO K = 1, N(ng)
!
! temp and salt
!
          PLUME(ng) % sAm(is, K) = t(Isrc, Jsrc, K, nrhs, isalt)
          PLUME(ng) % tpAm(is, K) = t(Isrc, Jsrc, K, nrhs, itemp)
!
! convert potential temp to in-situ temp
!
          prRef = 101.d3*1.d-4
          pr = prRef + &
            & (ABS(PLUME(ng) % zW(is, K))*rho_ref*g)*1.d-4  ! [dbar]
          CALL SW_TEMP(PLUME(ng) % sAm(is, K),  &
                     & PLUME(ng) % tpAm(is, K), &
                     & pr,                     &
                     & prRef,                  &
                     & PLUME(ng) % tAm(is, K))
!
! tracers
!
          IF (useTracers) THEN
            DO iTracer = 3, NT(ng)
              PLUME(ng) % trcAm(is, K, iTracer) = &
                & t(Isrc, Jsrc, K, nrhs, iTracer)
            ENDDO
          ENDIF
!
! u/v, w
!
          IF (Dsrc .EQ. 0) THEN
            PLUME(ng) % vAm(is, K) =     &
              & 0.5d0 * (               &
              & v(Isrc  , Jsrc, K, nrhs) +  &
              & v(Isrc+1, Jsrc, K, nrhs)    &
              & )
          ELSEIF (Dsrc .EQ. 1) THEN
            PLUME(ng) % vAm(is, K) =     &
              & 0.5d0 * (               &
              & u(Isrc, Jsrc  , K, nrhs) +  &
              & u(Isrc, Jsrc+1, K, nrhs)    &
              & )
          ENDIF
!
          PLUME(ng) % wAm(is, K) =       &
            & 0.5d0 * (                 &
            & wvel(Isrc, Jsrc, K-1) +       &
            & wvel(Isrc, Jsrc, K)           &
            & )
        ENDDO
!
! Discharge tracer concentration
!
        IF (useTracers .and. useInputTracers) THEN
          DO iTracer = 3, NT(ng)
            PLUME(ng) % trcIni(is, iTracer) = &
# ifdef ONE_TRACER_SOURCE
              & SOURCES(ng) % SGtrc(iTracer)
# else
              & SOURCES(ng) % SGtrc(is, iTracer)
# endif
          ENDDO
        ELSE
          DO iTracer = 3, NT(ng)
            PLUME(ng) % trcIni(is, iTracer) = 0.d0
          ENDDO
        ENDIF
!
! ==================================================================
!
        CALL iceplume_calc(ng, is, dx, dy, &
                         & fIni, tIni, sIni)
!
! ==================================================================!
!                                                                   !
! Rewrite Qsrc, T & S in SOURCES                                    !
!                                                                   !
! Update log                                                        !
! Instead read T/S directly from PLUME(ng)%t/s, read them from      !
! tracer values. This is done to accomodate the new checkCLF        !
! switch.                                                           !
! 2018/04/20 Chuning Wang                                           !
!                                                                   !
! Since entrainment / detrainment are separated in PLUME, here I    !
! updte the volume flux Qsrc and tracer flux Tsrc accordingly.      !
! 2018/05/21 Chuning Wang                                           !
!                                                                   !
! ==================================================================!
!
!         DO K = 1, N(ng)
!           Qsrc = (PLUME(ng) % det(is, K) +        &
!             &     PLUME(ng) % ent(is, K) +        &
!             &     PLUME(ng) % mAm(is, K) +        &
!             &     ABS(SOURCES(ng) % Qsrc(is, K)))
! !
! ! T, S, Tracer concentration
! !
!           IF (Qsrc .NE. 0.0d0) THEN
!             SOURCES(ng) % Tsrc(is, K, isalt) =            &
!               & (PLUME(ng) % det(is, K) *                 &
!               &  PLUME(ng) % trc(is, isalt) +             &
!               &  PLUME(ng) % ent(is, K) *                 &
!               &  PLUME(ng) % sAm(is, K) +                 &
!               &  PLUME(ng) % mAm(is, K) *                 &
!               &  PLUME(ng) % sB(is, K) +                  &
!               &  ABS(SOURCES(ng) % Qsrc(is, K)) *         &
!               &  SOURCES(ng) % Tsrc(is, K, isalt)) / Qsrc
!             SOURCES(ng) % Tsrc(is, K, itemp) =            &
!               & (PLUME(ng) % det(is, K) *                 &
!               &  PLUME(ng) % trc(is, itemp) +             &
!               &  PLUME(ng) % ent(is, K) *                 &
!               &  PLUME(ng) % tAm(is, K) +                 &
!               &  PLUME(ng) % mAm(is, K) *                 &
!               &  PLUME(ng) % tB(is, K) +                  &
!               &  ABS(SOURCES(ng) % Qsrc(is, K)) *         &
!               &  SOURCES(ng) % Tsrc(is, K, itemp)) / Qsrc
!             DO iTracer = 3, NT(ng)
!               SOURCES(ng) % Tsrc(is, K, iTracer) =            &
!                 & (PLUME(ng) % det(is, K) *                   &
!                 &  PLUME(ng) % trc(is, iTracer) +             &
!                 &  PLUME(ng) % ent(is, K) *                   &
!                 &  PLUME(ng) % trcAm(is, K, iTracer) +        &
!                 &  ABS(SOURCES(ng) % Qsrc(is, K)) *           &
!                 &  SOURCES(ng) % Tsrc(is, K, iTracer)) / Qsrc
!             ENDDO
! !
! ! If use background melt tracer, rewrite tracer concentration
! ! of the last tracer type.
! !
! !            IF (useBkgMelt .AND. useTracers .AND. &
! !              & useMeltTracers) THEN
! !              IF (D2src .GT. 0.d0) THEN
! !                SOURCES(ng) % Tsrc(is, K, NT(ng)) = &
! !                  & (PLUME(ng) % mAm(is, K) *       &
! !                  &  PLUME(ng) % trc(is, NT(ng))) / &
! !                  & Qsrc
! !              ELSEIF (D2src .LT. 0.d0) THEN
! !                SOURCES(ng) % Tsrc(is, K, NT(ng)) = &
! !                  & (PLUME(ng) % mAm(is, K) *       &
! !                  &  PLUME(ng) % trc(is, NT(ng))) / &
! !                  & (-Qsrc)
! !              ELSE
! !                SOURCES(ng) % Tsrc(is, K, NT(ng)) = 0.d0
! !              ENDIF
! !            ENDIF
!           ELSE  ! Qsrc .EQ. 0.0d0
!             SOURCES(ng) % Tsrc(is, K, isalt) = 0.d0
!             SOURCES(ng) % Tsrc(is, K, itemp) = 0.d0
!             DO iTracer = 3, NT(ng)
!               SOURCES(ng) % Tsrc(is, K, iTracer) = 0.d0
!             ENDDO
!           ENDIF
! !
! ! Qsrc
! !
!           IF (D2src .GT. 0.d0) THEN
!             SOURCES(ng) % Qsrc(is, K) = Qsrc
!           ELSEIF (D2src .LT. 0.d0) THEN
!             SOURCES(ng) % Qsrc(is, K) = -Qsrc
!           ELSE
!             SOURCES(ng) % Qsrc(is, K) = 0.d0
!           ENDIF
!         ENDDO
      ENDIF
    ENDDO
  ENDIF
END SUBROUTINE iceplume_tile
#endif
END MODULE iceplume_mod
