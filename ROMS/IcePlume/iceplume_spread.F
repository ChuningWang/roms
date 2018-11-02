#include "cppdefs.h"
MODULE iceplume_spread_mod
#ifdef ICEPLUME_spread
!
!=======================================================================
!                                                                      !
!  This subroutine spreads the iceplume detrainment in several grids.  !
!                                                                      !
!=======================================================================
!
implicit none
!
PRIVATE
PUBLIC :: iceplume_spread
!
CONTAINS
!
!***********************************************************************
SUBROUTINE iceplume_spread (ng, tile)
!***********************************************************************
!
! Imported variable declarations.
!
  integer, intent(in) :: ng, tile
!
  CALL iceplume_spread_tile (ng, tile)
  RETURN
END SUBROUTINE iceplume_spread
!
!***********************************************************************
SUBROUTINE iceplume_spread_tile(ng, tile)
!***********************************************************************
!
! Load public modules
!
  USE mod_scalars, ONLY : LuvSrc
  USE mod_sources
  USE mod_param, ONLY : BOUNDS
!
! Load private module
!
  USE mod_iceplume
!
!  Imported variable declarations.
!
  integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
  integer :: is, Isrc, Jsrc, Dsrc, is2, Isrc2, Jsrc2
  integer :: Gsrc, Gsrc2
  integer :: K, iTracer
  real(r8) :: cff, cff1, cff2, cff3, cff4
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
        Dsrc = INT(SOURCES(ng) % Dsrc(is))
        Gsrc = INT(SOURCES(ng) % SGgid(is))
        IF ( SOURCES(ng) % SGbar(is) .GT. 0.0) THEN
          cff = 0.0d0
          DO is2=1,Nsrc(ng)
            Isrc2 = SOURCES(ng) % Isrc(is2)
            Jsrc2 = SOURCES(ng) % Jsrc(is2)
            Gsrc2 = INT(SOURCES(ng) % SGgid(is2))
!
!=======================================================================!
!                                                                       !
!  Find the grids next to the runoff grid.                              !
!                                                                       !
!=======================================================================!
!
  IF ( ((Dsrc .EQ. 0) .AND. (Isrc .EQ. Isrc2) .AND.               &
    &   (ABS(Jsrc-Jsrc2) .EQ. 1) .AND. (Gsrc .EQ. Gsrc2)) .OR.    &
    &  ((Dsrc .EQ. 1) .AND. (Jsrc .EQ. Jsrc2) .AND.               &
    &   (ABS(Isrc-Isrc2) .EQ. 1) .AND. (Gsrc .EQ. Gsrc2)) ) THEN
    cff = cff+1.0
    DO K = 1, N(ng)
      PLUME(ng) % det(is2, K) = 0.25d0*PLUME(ng) % det(is, K)
    ENDDO
    DO iTracer = 1, NT(ng)
      PLUME(ng) % trc(is2, iTracer) = PLUME(ng) % trc(is, iTracer)
    ENDDO
  ENDIF
!
!=======================================================================!
!
          ENDDO
          IF (cff .GT. 0) THEN
            DO K = 1, N(ng)
              PLUME(ng) % det(is, K) = &
                & (1.0d0-cff*0.25)*PLUME(ng) % det(is, K)
            ENDDO
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  ENDIF
END SUBROUTINE iceplume_spread_tile
#endif
END MODULE iceplume_spread_mod
