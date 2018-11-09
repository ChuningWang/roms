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
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 4, __LINE__, __FILE__)
# endif
  CALL iceplume_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                IminS, ImaxS, JminS, JmaxS)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 4, __LINE__, __FILE__)
# endif

  RETURN
END SUBROUTINE iceplume
!
!***********************************************************************
SUBROUTINE iceplume_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS)
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
  integer, intent(in) :: LBi, UBi, LBj, UBj
  integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
  integer :: is, Isrc, Jsrc, Dsrc, D2src
  integer :: K, itrc
  real(r8) :: pr, prRef
  real(r8) :: dx, dy
  real(r8) :: cff, cff1, cff2, cff3, cff4
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
! Calculate virtual vertical mixing.                                !
!                                                                   !
! ==================================================================!
!
      PLUME(ng) % fW(is, 0) = 0.0d0
      DO K = 1, N(ng)
        PLUME(ng) % fW(is, K) = SUM( &
          & PLUME(ng) % det(is, 1:K) + PLUME(ng) % ent(is, 1:K)+ &
          & PLUME(ng) % mAm(is, 1:K))
      ENDDO
!
! Subtract the vertical velocities of the moving S-coordinates.
!
      cff1 = PLUME(ng) % fW(is, N(ng))/ &
        & (PLUME(ng) % zW(is, N(ng))-PLUME(ng) % zW(is, 0))
      DO K = 1, N(ng)
        PLUME(ng) % fW(is, K) = PLUME(ng) % fW(is, K) - &
          & cff1*(PLUME(ng) % zW(is, K)-PLUME(ng) % zW(is, 0))
      ENDDO
      if (is.eq.2) write(cff1)
# endif
      ENDIF
    ENDDO
  ENDIF
  RETURN
END SUBROUTINE iceplume_tile
#endif
