#include "cppdefs.h"
#ifdef ICEPLUME
!
!=======================================================================
!                                                                      !
!  This subroutine is a wraper for the ICEPLUME module.                !
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
  integer :: K, itrc
  real(r8) :: pr, prRef
  real(r8) :: dx, dy
  real(r8) :: cff, cff1, cff2, cff3, cff4
# ifdef ICEPLUME_VIRTUAL_MIX
  real(r8), dimension(0:N(ng)) :: dzNew, zWNew, zLNew
# endif
!
# include "set_bounds.h"
!
  IF (LuvSrc(ng)) THEN
    DO is=1,Nsrc(ng)
      Isrc = SOURCES(ng) % Isrc(is)
      Jsrc = SOURCES(ng) % Jsrc(is)
      IF ( ((IstrR .LE. Isrc) .AND. (Isrc .LE. IendR)) .AND. &
        &  ((JstrR .LE. Jsrc) .AND. (Jsrc .LE. JendR)) ) THEN
!
! ==================================================================!
!                                                                   !
! Find the right Rho coordinates                                    !
!                                                                   !
! ==================================================================!
!
        Dsrc = INT(SOURCES(ng) % Dsrc(is))
        IF (Dsrc .EQ. 0) THEN
          IF ( (GRID(ng) % rmask(Isrc  , Jsrc) .EQ. 0) .AND.      &
            &  (GRID(ng) % rmask(Isrc-1, Jsrc) .EQ. 1) ) THEN
            Isrc = Isrc-1
            D2src = -1
          ELSEIF ( (GRID(ng) % rmask(Isrc  , Jsrc) .EQ. 1) .AND.  &
                &  (GRID(ng) % rmask(Isrc-1, Jsrc) .EQ. 0) ) THEN
            D2src = 1
          ELSE
            D2src = 0
          ENDIF
        ELSEIF (Dsrc .EQ. 1) THEN
          IF ( (GRID(ng) % rmask(Isrc, Jsrc  ) .EQ. 0) .AND.      &
            &  (GRID(ng) % rmask(Isrc, Jsrc-1) .EQ. 1)) THEN
            Jsrc = Jsrc-1
            D2src = -1
          ELSEIF ((GRID(ng) % rmask(Isrc, Jsrc  ) .EQ. 1) .AND.   &
                & (GRID(ng) % rmask(Isrc, Jsrc-1) .EQ. 0)) THEN
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
!
! temp and salt
!
          PLUME(ng) % sAm(is, K) = &
            & OCEAN(ng) % t(Isrc, Jsrc, K, nstp(ng), isalt)
          PLUME(ng) % tpAm(is, K) = &
            & OCEAN(ng) % t(Isrc, Jsrc, K, nstp(ng), itemp)
!
! convert potential temp to in-situ temp
!
          prRef = 101.d3*1.d-4
          pr = prRef + &
            & (ABS(PLUME(ng) % zW(is, K))*rho_ref*g)*1.d-4  ! [dbar]
          CALL SW_TEMP(PLUME(ng) % sAm(is, K),  &
                     & PLUME(ng) % tpAm(is, K), &
                     & pr,                      &
                     & prRef,                   &
                     & PLUME(ng) % tAm(is, K))
!
! tracers
!
          DO itrc = 1, NT(ng)
            PLUME(ng) % trcAm(is, K, itrc) = &
              & OCEAN(ng) % t(Isrc, Jsrc, K, nstp(ng), itrc)
          ENDDO
!
! u/v, w
!
          IF (Dsrc .EQ. 0) THEN
            PLUME(ng) % vAm(is, K) = 0.5d0 * (              &
              & OCEAN(ng) % v(Isrc  , Jsrc, K, nstp(ng)) +  &
              & OCEAN(ng) % v(Isrc+1, Jsrc, K, nstp(ng)) )
          ELSEIF (Dsrc .EQ. 1) THEN
            PLUME(ng) % vAm(is, K) = 0.5d0 * (              &
              & OCEAN(ng) % u(Isrc, Jsrc  , K, nstp(ng)) +  &
              & OCEAN(ng) % u(Isrc, Jsrc+1, K, nstp(ng)) )
          ENDIF
          PLUME(ng) % wAm(is, K) = 0.5d0 * (                &
            & OCEAN(ng) % wvel(Isrc, Jsrc, K-1) +           &
            & OCEAN(ng) % wvel(Isrc, Jsrc, K  ) )
        ENDDO
!
! Discharge tracer concentration
!
        IF (useTracers .and. useInputTracers) THEN
          DO itrc = 3, NT(ng)
            PLUME(ng) % trcIni(is, itrc) = &
# ifdef ONE_TRACER_SOURCE
              & SOURCES(ng) % SGtrc(itrc)
# else
              & SOURCES(ng) % SGtrc(is, itrc)
# endif
          ENDDO
        ELSE
          DO itrc = 3, NT(ng)
            PLUME(ng) % trcIni(is, itrc) = 0.d0
          ENDDO
        ENDIF
!
! ==================================================================!
!                                                                   !
! Call the main function here.                                      !
!                                                                   !
! ==================================================================!
!
        IF (INT(SOURCES(ng) % SGtyp(is)) .NE. 1) THEN
          CALL iceplume_calc(ng, is, dx, dy,                  &
                           & ABS(SOURCES(ng) % SGbar(is)),    &
# ifdef ONE_TRACER_SOURCE
                           & SOURCES(ng) % SGtrc(itemp),      &
                           & SOURCES(ng) % SGtrc(isalt),      &
# else
                           & SOURCES(ng) % SGtrc(is, itemp),  &
                           & SOURCES(ng) % SGtrc(is, isalt),  &
# endif
                           & INT(SOURCES(ng) % SGtyp(is)),    &
                           & SOURCES(ng) % SGdep(is),         &
                           & SOURCES(ng) % SGlen(is))
        ENDIF
# ifdef ICEPLUME_VIRTUAL_MIX
!
! ==================================================================!
!                                                                   !
! Calculate virtual vertical mixing if grid is not adjoint to land. !
!                                                                   !
! ==================================================================!
!
        IF (INT(PLUME(ng) % dir(is)) .EQ. 0) THEN
!
! Calculate new box vertical boundary.
!
          zLNew(0) = PLUME(ng) % zW(is, 0)
          DO K = 1, N(ng)
            dzNew(K) = PLUME(ng) % dz(is, K) +                  &
              &        (PLUME(ng) % det(is, K) +                &
              &         PLUME(ng) % ent(is, K) +                &
              &         PLUME(ng) % mAm(is, K))*dt(ng)/(dx*dy)
            zLNew(K) = zLNew(K-1) + dzNew(K)
          ENDDO
!
! Calculate tracer concentration.
!
          DO K = 1, N(ng)
            cff  = dx*dy*PLUME(ng) % dz(is, K)
            cff1 = (PLUME(ng) % det(is, K) +        &
              &     PLUME(ng) % ent(is, K) +        &
              &     PLUME(ng) % mAm(is, K))*dt(ng)
            DO itrc = 1, NT(ng)
              cff2 = (PLUME(ng) % det(is, K)*                 &
                &     PLUME(ng) % trc(is, itrc)+              &
                &     PLUME(ng) % ent(is, K)*                 &
                &     PLUME(ng) % trcAm(is, K, itrc)+         &
                &     PLUME(ng) % mAm(is, K)*                 &
                &     PLUME(ng) % trcB(is, K, itrc))*dt(ng)
              PLUME(ng) % trcNew(is, K, itrc) =                   &
                & (PLUME(ng) % trcAm(is, K, itrc)*cff + cff2) /   &
                & (cff + cff1)
            ENDDO
          ENDDO
!
! Correction for grid stretching.
!
          cff1 = PLUME(ng) % zW(is, N(ng)) - PLUME(ng) % zW(is, 0)
          cff2 = zLNew(N(ng))              - PLUME(ng) % zW(is, 0)
          DO K = 0, N(ng)
            zWNew(K) = &
              & (PLUME(ng) % zW(is, K)-zWNew(0))*cff2/cff1+zWNew(0)
          ENDDO
!
          cff  = zWNew(1) - zWNew(0)
          cff1 = MAX(zWNew(1) - zLNew(1), 0.0d0)
          cff3 = cff - cff1
          DO itrc = 1, NT(ng)
            PLUME(ng) % trcNew(is, 1, itrc) =                  &
              & (PLUME(ng) % trcNew(is, 1, itrc)*cff3 +        &
              &  PLUME(ng) % trcNew(is, 2, itrc)*cff1) / cff
          ENDDO
!
          cff  = zWNew(N(ng)) - zWNew(N(ng)-1)
          cff2 = -MIN(zWNew(N(ng)-1)-zLNew(N(ng)-1), 0.0d0)
          cff3 = cff - cff2
          DO itrc = 1, NT(ng)
            PLUME(ng) % trcNew(is, N(ng), itrc) =                   &
              & (PLUME(ng) % trcNew(is, N(ng)  , itrc)*cff3 +       &
              &  PLUME(ng) % trcNew(is, N(ng)-1, itrc)*cff2) / cff
          ENDDO
!
          DO K = 2, N(ng)-1
            cff  = zWNew(K) - zWNew(K-1)
            cff1 = MAX(zWNew(K) - zLNew(K), 0.0d0)
            cff2 = -MIN(zWNew(K-1)-zLNew(K-1), 0.0d0)
            cff3 = cff - cff1 - cff2
            DO itrc = 1, NT(ng)
              PLUME(ng) % trcNew(is, K, itrc) =               &
                & (PLUME(ng) % trcNew(is, K+1, itrc)*cff1 +   &
                &  PLUME(ng) % trcNew(is, K-1, itrc)*cff2 +   &
                &  PLUME(ng) % trcNew(is, K  , itrc)*cff3)/   &
                & cff
            ENDDO
          ENDDO
!
! Scale tracers with original omega levels.
!
          cff = (zWNew(N(ng))-zWnew(0)) / &
            &   (PLUME(ng) % zW(is, N(ng))-PLUME(ng) % zW(is, 0))
          DO K = 1, N(ng)
            DO itrc = 1, NT(ng)
              PLUME(ng) % trcNew(is, K, itrc) =       &
                & PLUME(ng) % trcNew(is, K, itrc)*cff
            ENDDO
          ENDDO
!
! Substitute new tracer values back to OCEAN.
!
          DO K = 1, N(ng)
            DO itrc = 1, NT(ng)
              OCEAN(ng) % t(Isrc, Jsrc, K, nstp(ng), itrc) =  &
                & PLUME(ng) % trcNew(is, K, itrc)
            ENDDO
          ENDDO
        ENDIF
# endif
      ENDIF
    ENDDO
  ENDIF
  RETURN
END SUBROUTINE iceplume_tile
#endif
