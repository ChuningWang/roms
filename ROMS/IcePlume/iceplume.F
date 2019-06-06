#include "cppdefs.h"
#ifdef ICEPLUME
!
!=======================================================================
!                                                                      !
!  This subroutine is a wraper for the iceplume module.                !
!                                                                      !
!=======================================================================
!
!***********************************************************************
SUBROUTINE iceplume (ng, tile)
!***********************************************************************
!
  USE mod_param
!
!  Imported variable declarations.
!
  integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# ifdef PROFILE
  CALL wclock_on (ng, iNLM, 4, __LINE__, __FILE__)
# endif
  CALL iceplume_tile (ng, tile)
# ifdef PROFILE
  CALL wclock_off (ng, iNLM, 4, __LINE__, __FILE__)
# endif

  RETURN
END SUBROUTINE iceplume
!
!***********************************************************************
SUBROUTINE iceplume_tile (ng, tile)
!***********************************************************************
!
  USE mod_grid, ONLY : GRID
  USE mod_ocean, ONLY : OCEAN
  USE mod_stepping, ONLY : nstp
  USE mod_scalars, ONLY : dt, LuvSrc
  USE mod_sources
  USE mod_iceplume
!
!  Imported variable declarations.
!
  integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
  integer :: is, Isrc, Jsrc, Dsrc, D2src
  integer :: K, K2, itrc
  real(r8) :: pr, prRef = 101.d3*1.d-4
  real(r8) :: dx, dy
  real(r8) :: cff, cff1, cff2, cff3, cff4
  real(r8) :: RHO
!
# include "set_bounds.h"
!
  IF (LuvSrc(ng)) THEN
    DO is=1,Nsrc(ng)
      Isrc = SOURCES(ng) % Isrc(is)
      Jsrc = SOURCES(ng) % Jsrc(is)
      IF ( ((IstrR .LE. Isrc) .AND. (Isrc .LE. IendR)) .AND.            &
     &     ((JstrR .LE. Jsrc) .AND. (Jsrc .LE. JendR)) ) THEN
!
! ==================================================================!
!                                                                   !
! Find the right Rho coordinates.                                   !
!                                                                   !
! ==================================================================!
!
        Dsrc = NINT(SOURCES(ng) % Dsrc(is))
        IF (Dsrc .EQ. 0) THEN
          IF ( (GRID(ng) % rmask(Isrc  , Jsrc) .EQ. 0) .AND.            &
     &         (GRID(ng) % rmask(Isrc-1, Jsrc) .EQ. 1) ) THEN
            Isrc = Isrc-1
            D2src = -1
          ELSEIF ( (GRID(ng) % rmask(Isrc  , Jsrc) .EQ. 1) .AND.        &
     &             (GRID(ng) % rmask(Isrc-1, Jsrc) .EQ. 0) ) THEN
            D2src = 1
          ELSE
            D2src = 0
          ENDIF
        ELSEIF (Dsrc .EQ. 1) THEN
          IF ( (GRID(ng) % rmask(Isrc, Jsrc  ) .EQ. 0) .AND.            &
     &         (GRID(ng) % rmask(Isrc, Jsrc-1) .EQ. 1)) THEN
            Jsrc = Jsrc-1
            D2src = -1
          ELSEIF ((GRID(ng) % rmask(Isrc, Jsrc  ) .EQ. 1) .AND.         &
     &            (GRID(ng) % rmask(Isrc, Jsrc-1) .EQ. 0)) THEN
            D2src = 1
          ELSE
            D2src = 0
          ENDIF
        ENDIF
        PLUME(ng) % dir(is) = REAL(D2src)
!
! Get grid dimensions
!
        IF (Dsrc .EQ. 0) THEN
          dx = ABS(GRID(ng) % om_r(Isrc, Jsrc))
          dy = ABS(GRID(ng) % on_r(Isrc, Jsrc))
        ELSE
          dx = ABS(GRID(ng) % on_r(Isrc, Jsrc))
          dy = ABS(GRID(ng) % om_r(Isrc, Jsrc))
        ENDIF
!
! ==================================================================!
!                                                                   !
! Read in profiles from OCEAN.                                      !
!                                                                   !
! ==================================================================!
!
        DO K = 0, N(ng)
          PLUME(ng) % zW(is, K) = GRID(ng) % z_w(Isrc, Jsrc, K)
        ENDDO
!
        DO K = 1, N(ng)
          PLUME(ng) % zR(is, K) =                                       &
     &        0.5d0 * (PLUME(ng) % zW(is, K-1) + PLUME(ng) % zW(is, K))
          PLUME(ng) % dz(is, K) =                                       &
     &        PLUME(ng) % zW(is, K) - PLUME(ng) % zW(is, K-1)
!
! temp, salt and density
!
          PLUME(ng) % sAm(is, K) =                                      &
     &        OCEAN(ng) % t(Isrc, Jsrc, K, nstp(ng), isalt)
          PLUME(ng) % tpAm(is, K) =                                     &
     &        OCEAN(ng) % t(Isrc, Jsrc, K, nstp(ng), itemp)
!
! convert potential temp to in-situ temp
!
          pr = prRef+(ABS(PLUME(ng) % zR(is, K))*rhoRef*g)*1.d-4
          CALL SW_TEMP(PLUME(ng) % sAm(is, K), PLUME(ng) % tpAm(is, K), &
     &                 pr, prRef, PLUME(ng) % tAm(is, K))
!
! tracers
!
          DO itrc = 1, NT(ng)
            PLUME(ng) % trcAm(is, K, itrc) =                            &
     &          OCEAN(ng) % t(Isrc, Jsrc, K, nstp(ng), itrc)
          ENDDO
!
! u/v, w
!
          IF (Dsrc .EQ. 0) THEN
            PLUME(ng) % vAm(is, K) = 0.5d0 * (                          &
     &          OCEAN(ng) % v(Isrc  , Jsrc, K, nstp(ng)) +              &
     &          OCEAN(ng) % v(Isrc+1, Jsrc, K, nstp(ng)) )
          ELSEIF (Dsrc .EQ. 1) THEN
            PLUME(ng) % vAm(is, K) = 0.5d0 * (                          &
     &          OCEAN(ng) % u(Isrc, Jsrc  , K, nstp(ng)) +              &
     &          OCEAN(ng) % u(Isrc, Jsrc+1, K, nstp(ng)) )
          ENDIF
          PLUME(ng) % wAm(is, K) = 0.5d0 * (                            &
     &        OCEAN(ng) % wvel(Isrc, Jsrc, K-1) +                       &
     &        OCEAN(ng) % wvel(Isrc, Jsrc, K  ) )
        ENDDO
!
! Calculate rho-layer depth, thickness, and ambient density
!
        IF (iic(ng) .EQ. ntstart(ng)) THEN
          PLUME(ng) % rhoAm(is, N(ng)+1) = rhoAir
        ENDIF
# ifdef ICEPLUME_DETRAIN_AVERAGE
        cff4 = REAL(                                                    &
     &      (SOURCES(ng) % SGIrange(is, 2)-                             &
     &       SOURCES(ng) % SGIrange(is, 1)+1) *                         &
     &      (SOURCES(ng) % SGJrange(is, 2)-                             &
     &       SOURCES(ng) % SGJrange(is, 1)+1))
        IF (iic(ng) .EQ. ntstart(ng)) THEN
          PLUME(ng) % wAvg = MAX(MIN(PLUME(ng) % wAvg, 1.0), 0.0)
        ENDIF
# endif
        DO K = 1, N(ng)
# ifdef ICEPLUME_DETRAIN_AVERAGE
!
! Average density profiles in a box to surpress some oscillation
!
          cff1 = SUM( OCEAN(ng) % t(                                    &
     & SOURCES(ng) % SGIrange(is, 1):SOURCES(ng) % SGIrange(is, 2),     &
     & SOURCES(ng) % SGJrange(is, 1):SOURCES(ng) % SGJrange(is, 2),     &
     & K, nstp(ng), isalt) ) / cff4
          cff3 = SUM( OCEAN(ng) % t(                                    &
     & SOURCES(ng) % SGIrange(is, 1):SOURCES(ng) % SGIrange(is, 2),     &
     & SOURCES(ng) % SGJrange(is, 1):SOURCES(ng) % SGJrange(is, 2),     &
     & K, nstp(ng), itemp) ) / cff4
          pr = prRef+(ABS(PLUME(ng) % zR(is, K))*rhoRef*g)*1.d-4
          CALL SW_TEMP(cff1, cff3, pr, prRef, cff2)
          cff = RHO(cff2, cff1, PLUME(ng) % zR(is, K))
          IF (iic(ng) .EQ. ntstart(ng)) THEN
            PLUME(ng) % rhoAm(is, K) = cff
          ELSE
            PLUME(ng) % rhoAm(is, K) =                                  &
     &          PLUME(ng) % rhoAm(is, K)*(1.0-PLUME(ng) % wAvg) +       &
     &          cff*PLUME(ng) % wAvg
          ENDIF
# else
          PLUME(ng) % rhoAm(is, K) =                                    &
     &        RHO(PLUME(ng) % tAm(is, K),                               &
     &            PLUME(ng) % sAm(is, K),                               &
     &            PLUME(ng) % zR(is, K))
# endif
        ENDDO
!
! Discharge tracer concentration
!
# ifdef ICEPLUME_TRACER
        PLUME(ng) % trcIni(is, itemp) = 0.0
        PLUME(ng) % trcIni(is, isalt) = 0.0
        DO itrc = 3, NT(ng)
          PLUME(ng) % trcIni(is, itrc) =                                &
#  ifdef ONE_TRACER_SOURCE
     &        SOURCES(ng) % SGtrc(itrc)
#  else
     &        SOURCES(ng) % SGtrc(is, itrc)
#  endif
        ENDDO
# else
        DO itrc = 1, NT(ng)
          PLUME(ng) % trcIni(is, itrc) = 0.d0
        ENDDO
# endif
!
! ==================================================================!
!                                                                   !
! Call the main function here.                                      !
!                                                                   !
! ==================================================================!
!
        IF (NINT(SOURCES(ng) % SGtyp(is)) .NE. 1) THEN
          CALL iceplume_calc(ng, is, dx, dy,                            &
     &                       ABS(SOURCES(ng) % SGbar(is)),              &
# ifdef ONE_TRACER_SOURCE
     &                       SOURCES(ng) % SGtrc(itemp),                &
     &                       SOURCES(ng) % SGtrc(isalt),                &
# else
     &                       SOURCES(ng) % SGtrc(is, itemp),            &
     &                       SOURCES(ng) % SGtrc(is, isalt),            &
# endif
     &                       NINT(SOURCES(ng) % SGtyp(is)),             &
     &                       SOURCES(ng) % SGdep(is),                   &
     &                       SOURCES(ng) % SGlen(is))
        ENDIF
# ifdef ICEPLUME_VIRTUAL_MIX
!
! ==================================================================!
!                                                                   !
! Calculate virtual vertical mixing if grid is not adjoint to land. !
!                                                                   !
! ==================================================================!
!
! Calculate new box vertical boundary.
!
        PLUME(ng) % zLNew(is, 0) = PLUME(ng) % zW(is, 0)
        DO K = 1, N(ng)
          PLUME(ng) % dzNew(is, K) = PLUME(ng) % dz(is, K) +            &
     &               (PLUME(ng) % det(is, K) +                          &
     &                PLUME(ng) % ent(is, K) +                          &
     &                PLUME(ng) % mB(is, K))*dt(ng)/(dx*dy)
          PLUME(ng) % zLNew(is, K) =                                    &
     &        PLUME(ng) % zLNew(is, K-1) +                              &
     &        PLUME(ng) % dzNew(is, K)
        ENDDO
!
! Calculate tracer concentration.
!
        DO K = 1, N(ng)
          cff  = dx*dy*PLUME(ng) % dz(is, K)
          cff1 = (PLUME(ng) % det(is, K) +                              &
     &            PLUME(ng) % ent(is, K) +                              &
     &            PLUME(ng) % mB(is, K))*dt(ng)
          DO itrc = 1, NT(ng)
            cff2 = (PLUME(ng) % det(is, K)*                             &
     &              PLUME(ng) % trc(is, itrc)+                          &
     &              PLUME(ng) % ent(is, K)*                             &
     &              PLUME(ng) % trcAm(is, K, itrc)+                     &
     &              PLUME(ng) % mB(is,k)*                               &
     &              PLUME(ng)%trcB(is,k,itrc)+                          &
     &              PLUME(ng) % trcAmToB(is,k,itrc))*dt(ng)
            PLUME(ng) % trcL(is, K, itrc) =                             &
     &          (PLUME(ng) % trcAm(is, K, itrc)*cff + cff2) /           &
     &          (cff + cff1)
          ENDDO
        ENDDO
!
! Correction for grid stretching.
!
        cff1 = PLUME(ng) % zW(is, N(ng)) -                              &
     &         PLUME(ng) % zW(is, 0)
        cff2 = PLUME(ng) % zLNew(is, N(ng)) -                           &
     &         PLUME(ng) % zW(is, 0)
        DO K = 0, N(ng)
          PLUME(ng) % zWNew(is, K) =                                    &
     &        (PLUME(ng) % zW(is, K) - PLUME(ng) %zWNew(is, 0))*        &
     &        cff2/cff1 + PLUME(ng) % zWNew(is, 0)
        ENDDO
!
! Calculate weight function.
!
        DO K = 1, N(ng)
          cff = PLUME(ng) % zWNew(is, K) -                              &
     &          PLUME(ng) % zWNew(is, K-1)
          DO K2 = 1, N(ng)
            cff1 = PLUME(ng) % zWNew(is, K) -                           &
     &             PLUME(ng) % zLNew(is, K2)
            cff1 = MIN(MAX(cff1, 0.0_r8), cff)
            cff2 = PLUME(ng) % zWNew(is, K-1) -                         &
     &             PLUME(ng) % zLNew(is, K2-1)
            cff2 = MIN(-MIN(cff2, 0.0_r8), cff)
            PLUME(ng) % dzWei(is, K2) = (cff-cff1-cff2)/cff
          ENDDO
          DO itrc = 1, NT(ng)
            PLUME(ng) % trcNew(is, K, itrc) = 0.0_r8
            DO K2 = 1, N(ng)
              PLUME(ng) % trcNew(is, K, itrc) =                         &
     &            PLUME(ng) % trcNew(is, K, itrc) +                     &
     &            PLUME(ng) % dzWei(is, K2)*                            &
     &            PLUME(ng) % trcL(is, K2, itrc)
            ENDDO
          ENDDO
        ENDDO
!
! Scale tracers with original omega levels.
!
        cff = (PLUME(ng) % zWNew(is, N(ng)) -                           &
     &         PLUME(ng) % zWNew(is, 0)) /                              &
     &        (PLUME(ng) % zW(is, N(ng)) - PLUME(ng) % zW(is, 0))
        DO K = 1, N(ng)
          DO itrc = 1, NT(ng)
            PLUME(ng) % trcNew(is, K, itrc) =                           &
     &          PLUME(ng) % trcNew(is, K, itrc)*cff
          ENDDO
        ENDDO
!
! Substitute new tracer values back to OCEAN.
!
        DO K = 1, N(ng)
          DO itrc = 1, NT(ng)
            OCEAN(ng) % t(Isrc, Jsrc, K, nstp(ng), itrc) =              &
     &          PLUME(ng) % trcNew(is, K, itrc)
          ENDDO
        ENDDO
# endif
      ENDIF
    ENDDO
  ENDIF
  RETURN
END SUBROUTINE iceplume_tile
#endif
