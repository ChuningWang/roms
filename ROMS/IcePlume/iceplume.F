#include "cppdefs.h"
MODULE iceplume_mod
#ifdef ICEPLUME
!
!=======================================================================
!                                                                      !
!  This subroutine is a wraaper for the ICEPLUME module.               !
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
  integer :: I, II, JI, DIR, DIR2
  integer :: K, iTracer
  real(r8) :: pr, prRef
  real(r8) :: dx, dy
  real(r8) :: QIni, TIni, SIni
  real(r8) :: Qsrc
!
# include "set_bounds.h"
!
! ==================================================================
!
  IF (LuvSrc(ng)) THEN
    DO I=1,Nsrc(ng)
      II = SOURCES(ng) % Isrc(I)
      JI = SOURCES(ng) % Jsrc(I)
      IF ( ((IstrR .LE. II) .AND. (II .LE. IendR)) .AND. &
        &  ((JstrR .LE. JI) .AND. (JI .LE. JendR)) ) THEN
!
! ==================================================================
! Find the right Rho coordinates
! ==================================================================
!
        DIR = SOURCES(ng) % Dsrc(I)
        IF (DIR .EQ. 0) THEN
          IF ( (rmask(II  , JI) .EQ. 0) .AND. &
            &  (rmask(II-1, JI) .EQ. 1) ) THEN
            II = II-1
            DIR2 = -1
          ELSE
            DIR2 = 1
          ENDIF
        ELSEIF (DIR .EQ. 1) THEN
          IF ((rmask(II, JI  ) .EQ. 0) .AND. &
            & (rmask(II, JI-1) .EQ. 1)) THEN
            JI = JI-1
            DIR2 = -1
          ELSE
            DIR2 = 1
          ENDIF
        ENDIF
!
! Get grid dimensions
!
        IF ((SOURCES(ng) % Dsrc(I)) .EQ. 0) THEN  ! read in v direction
          dx = ABS(om_r(II, JI))
          dy = ABS(on_r(II, JI))
        ELSE
          dx = ABS(on_r(II, JI))
          dy = ABS(om_r(II, JI))
        ENDIF
!
! Get initial condition for iceplume model
!
        QIni = ABS(SOURCES(ng) % Qsrc(I, 1))
        TIni = &
# ifdef ONE_TRACER_SOURCE
          & SOURCES(ng) % Tsrc(itemp)
# elif defined TWO_D_TRACER_SOURCE
          & SOURCES(ng) % Tsrc(I, itemp)
# else
          & SOURCES(ng) % Tsrc(I, 1, itemp)
# endif
        SIni = &
# ifdef ONE_TRACER_SOURCE
          & SOURCES(ng) % Tsrc(isalt)
# elif defined TWO_D_TRACER_SOURCE
          & SOURCES(ng) % Tsrc(I, isalt)
# else
          & SOURCES(ng) % Tsrc(I, 1, isalt)
# endif
        SOURCES(ng) % Qsrc(I, 1) = 0.d0
!
! ==================================================================
! Read in profiles from OCEAN.
! ==================================================================
!
        DO K = 0, N(ng)
          PLUME(ng) % zW(I, K) = z_w(II, JI, K)
        ENDDO
!
        DO K = 1, N(ng)
!
! temp and salt
!
          PLUME(ng) % sAm(I, K) = t(II, JI, K, nrhs, isalt)
          PLUME(ng) % tpAm(I, K) = t(II, JI, K, nrhs, itemp)
!
! convert potential temp to in-situ temp
!
          prRef = 101.d3*1.d-4
          pr = prRef + &
            & (ABS(PLUME(ng) % zW(I, K))*rho_ref*g)*1.d-4  ! [dbar]
          CALL SW_TEMP(PLUME(ng) % sAm(I, K),  &
                     & PLUME(ng) % tpAm(I, K), &
                     & pr,                     &
                     & prRef,                  &
                     & PLUME(ng) % tAm(I, K))
!
! tracers
!
          IF (useTracers) THEN
            DO iTracer = 3, NT(ng)
              PLUME(ng) % trcAm(I, K, iTracer) = &
                & t(II, JI, K, nrhs, iTracer)
            ENDDO
          ENDIF
!
! u/v, w
!
          IF (DIR .EQ. 0) THEN
            PLUME(ng) % vAm(I, K) =     &
              & 0.5d0 * (               &
              & v(II  , JI, K, nrhs) +  &
              & v(II+1, JI, K, nrhs)    &
              & )
          ELSEIF (DIR .EQ. 1) THEN
            PLUME(ng) % vAm(I, K) =     &
              & 0.5d0 * (               &
              & u(II, JI  , K, nrhs) +  &
              & u(II, JI+1, K, nrhs)    &
              & )
          ENDIF
!
          PLUME(ng) % wAm(I, K) =       &
            & 0.5d0 * (                 &
            & wvel(II, JI, K-1) +       &
            & wvel(II, JI, K)           &
            & )
        ENDDO
!
! Discharge tracer concentration
!
        IF (useTracers .and. useInputTracers) THEN
          DO iTracer = 3, NT(ng)
            PLUME(ng) % trcIni(I, iTracer) = &
# ifdef ONE_TRACER_SOURCE
              & SOURCES(ng) % Tsrc(iTracer)
# elif defined TWO_D_TRACER_SOURCE
              & SOURCES(ng) % Tsrc(I, iTracer)
# else
              & SOURCES(ng) % Tsrc(I, 1, iTracer)
# endif
          ENDDO
        ELSE
          DO iTracer = 3, NT(ng)
            PLUME(ng) % trcIni(I, iTracer) = 0.d0
          ENDDO
        ENDIF
!
! ==================================================================
!
        CALL iceplume_calc(ng, I, dx, dy, &
                         & QIni, TIni, SIni)
!
! ==================================================================
! Rewrite Qsrc, T & S in SOURCES
!
! Update log
! Instead read T/S directly from PLUME(ng)%t/s, read them from tracer
! values. This is done to accomodate the new checkCLF switch.
! 2018/04/20 Chuning Wang
!
! Since entrainment / detrainment are separated in PLUME, here I
! updte the volume flux Qsrc and tracer flux Tsrc accordingly.
! 2018/05/21 Chuning Wang
! ==================================================================
!
        DO K = 1, N(ng)
          Qsrc = (PLUME(ng) % det(I, K) + &
               &  PLUME(ng) % ent(I, K) + &
               &  ABS(SOURCES(ng) % Qsrc(I, K)))
!
! T, S, Tracer concentration
!
          IF (Qsrc .NE. 0.0d0) THEN
            SOURCES(ng) % Tsrc(I, K, isalt) = &
              & (PLUME(ng) % det(I, K) * &
              &  PLUME(ng) % trc(I, isalt) + &
              &  PLUME(ng) % ent(I, K) * &
              &  PLUME(ng) % sAm(I, K) + &
              &  ABS(SOURCES(ng) % Qsrc(I, K)) * &
              &  SOURCES(ng) % Tsrc(I, K, isalt)) / Qsrc
            SOURCES(ng) % Tsrc(I, K, itemp) = &
              & (PLUME(ng) % det(I, K) * &
              &  PLUME(ng) % trc(I, itemp) + &
              &  PLUME(ng) % ent(I, K) * &
              &  PLUME(ng) % tAm(I, K) + &
              &  ABS(SOURCES(ng) % Qsrc(I, K)) * &
              &  SOURCES(ng) % Tsrc(I, K, itemp)) / Qsrc
            DO iTracer = 3, NT(ng)
              SOURCES(ng) % Tsrc(I, K, iTracer) = &
                & (PLUME(ng) % det(I, K) * &
                &  PLUME(ng) % trc(I, iTracer) + &
                &  PLUME(ng) % ent(I, K) * &
                &  PLUME(ng) % trcAm(I, K, iTracer) + &
                &  ABS(SOURCES(ng) % Qsrc(I, K)) * &
                &  SOURCES(ng) % Tsrc(I, K, iTracer)) / Qsrc
            ENDDO
          ELSE
            SOURCES(ng) % Tsrc(I, K, isalt) = 0.d0
            SOURCES(ng) % Tsrc(I, K, itemp) = 0.d0
            DO iTracer = 3, NT(ng)
              SOURCES(ng) % Tsrc(I, K, iTracer) = 0.d0
            ENDDO
          ENDIF
!
! Qsrc
!
          IF (SOURCES(ng) % Qbar(I) .GT. 0.d0) THEN
            SOURCES(ng) % Qsrc(I, K) = Qsrc
          ELSEIF (SOURCES(ng) % Qbar(I) .LT. 0.d0) THEN
            SOURCES(ng) % Qsrc(I, K) = -Qsrc
          ELSE
            SOURCES(ng) % Qsrc(I, K) = 0.d0
          ENDIF
        ENDDO
!
! Use the tendency terms to update T, S & tracers in OCEAN
!
        IF (useBkgMelt) THEN
          DO K = 1, N(ng)
            DO iTracer = 1, NT(ng)
              t(II, JI, K, nrhs, iTracer) = &
                & t(II, JI, K, nrhs, iTracer) + &
                & PLUME(ng) % tendT(I, K, iTracer)
            ENDDO
            IF (useBkgTracers) THEN
!
! Hard code meltwater concentration to at least 0.
!
              t(II, JI, K, nrhs, NT(ng)) = &
                MAX(t(II, JI, K, nrhs, NT(ng)), 0.0d0)
            ENDIF
          ENDDO
        ENDIF
!
! Write to file for diagonse
!
!        IF (QIni .GT. 0.0d0) THEN
!          write(*, *)  II, JI
!          write(*, '(A3, 99 A12)')  'lev', 'zW', &
!            & 't', 's', 'r', 'w', 'a', 'mInt', 'volFlux'
!          DO K = 0, N(ng)
!            write(*, '(I4, 99 E12.4)')  K, &
!              & PLUME(ng) % zW(I, K), &
!              & PLUME(ng) % t(I, K), PLUME(ng) % s(I, K), &
!              & PLUME(ng) % r(I, K), PLUME(ng) % w(I, K), &
!              & PLUME(ng) % a(I, K), PLUME(ng) % mInt(I, K), &
!              & PLUME(ng) % volFlux(I, K)
!          ENDDO
!          write(*, '(A3, 99 A12)')  'lev', 'zR', &
!            & 'tAm', 'sAm', 'vAm', 'wAm', 'ent', 'det', &
!            & 'fwFlux', 'heatFlux'
!          DO K = 1, N(ng)
!            write(*, '(I4, 99 E12.4)')  K, PLUME(ng) % zR(I, K), &
!              & PLUME(ng) % tAm(I, K), PLUME(ng) % sAm(I, K), &
!              & PLUME(ng) % vAm(I, K), PLUME(ng) % wAm(I, K), &
!              & PLUME(ng) % ent(I, K), PLUME(ng) % det(I, K), &
!              & PLUME(ng) % fwFlux(I, K), PLUME(ng) % heatFlux(I, K)
!          ENDDO
!        ENDIF
!
      ENDIF
    ENDDO
  ENDIF
END SUBROUTINE iceplume_tile
#endif
END MODULE iceplume_mod
