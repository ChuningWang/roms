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
    &                 nstp(ng),                                    &
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
    &                    nstp,                                     &
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
  integer, intent(in) :: nstp
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
  integer :: K, itrc
  real(r8) :: pr, prRef
  real(r8) :: dx, dy
  real(r8) :: fIni, tIni, sIni
  real(r8) :: Qsrc
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
!
! ==================================================================!
!                                                                   !
! Find the right Rho coordinates                                    !
!                                                                   !
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
        PLUME(ng) % dir = real(D2src)
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
        fIni = ABS(SOURCES(ng) % SGbar(is))
        tIni = SOURCES(ng) % SGtrc(is, itemp)
        sIni = SOURCES(ng) % SGtrc(is, isalt)
!
! ==================================================================!
!                                                                   !
! Read in profiles from OCEAN.                                      !
!                                                                   !
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
          PLUME(ng) % sAm(is, K) = t(Isrc, Jsrc, K, nstp, isalt)
          PLUME(ng) % tpAm(is, K) = t(Isrc, Jsrc, K, nstp, itemp)
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
              & t(Isrc, Jsrc, K, nstp, itrc)
          ENDDO
!
! u/v, w
!
          IF (Dsrc .EQ. 0) THEN
            PLUME(ng) % vAm(is, K) =        &
              & 0.5d0 * (                   &
              & v(Isrc  , Jsrc, K, nstp) +  &
              & v(Isrc+1, Jsrc, K, nstp)    &
              & )
          ELSEIF (Dsrc .EQ. 1) THEN
            PLUME(ng) % vAm(is, K) =        &
              & 0.5d0 * (                   &
              & u(Isrc, Jsrc  , K, nstp) +  &
              & u(Isrc, Jsrc+1, K, nstp)    &
              & )
          ENDIF
!
          PLUME(ng) % wAm(is, K) =          &
            & 0.5d0 * (                     &
            & wvel(Isrc, Jsrc, K-1) +       &
            & wvel(Isrc, Jsrc, K)           &
            & )
        ENDDO
!
! Discharge tracer concentration
!
        IF (useTracers .and. useInputTracers) THEN
          DO itrc = 3, NT(ng)
            PLUME(ng) % trcIni(is, itrc) = &
              & SOURCES(ng) % SGtrc(is, itrc)
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
        CALL iceplume_calc(ng, is, dx, dy, &
                         & fIni, tIni, sIni)
!
! ==================================================================!
!                                                                   !
! Update SOURCES and OCEAN.                                         !
!                                                                   !
! ==================================================================!
!
!        DO K = 1, N(ng)
!          Qsrc = (PLUME(ng) % det(is, K) +         &
!            &     PLUME(ng) % ent(is, K) +         &
!            &     ABS(SOURCES(ng) % Qsrc(is, K)))
!!
!! T, S, Tracer concentration
!!
!          IF (Qsrc .NE. 0.0d0) THEN
!            DO itrc = 1, NT(ng)
!              SOURCES(ng) % Tsrc(is,K,itrc) =                   &
!                &  (PLUME(ng) % det(is, K) *                    &
!                &   PLUME(ng) % trc(is, itrc) +                 &
!                &   PLUME(ng) % ent(is, K) *                    &
!                &   PLUME(ng) % trcAm(is, K, itrc) +            &
!                &   ABS(SOURCES(ng) % Qsrc(is, K)) *            &
!                &   SOURCES(ng) % Tsrc(is,K,itrc)               &
!                &  ) / Qsrc
!            ENDDO
!          ELSE
!            DO itrc = 1, NT(ng)
!              SOURCES(ng) % Tsrc(is, K, itrc) = 0.d0
!            ENDDO
!          ENDIF
!!
!! Qsrc
!!
!          IF (D2src .GT. 0) THEN
!            SOURCES(ng) % Qsrc(is, K) = Qsrc
!          ELSEIF (D2src .LT. 0) THEN
!            SOURCES(ng) % Qsrc(is, K) = -Qsrc
!          ELSE
!            SOURCES(ng) % Qsrc(is, K) = 0.d0
!          ENDIF
!        ENDDO
!!
!! Use the tendency terms to update T, S & tracers in OCEAN
!!
!        IF (useBkgMelt) THEN
!          DO K = 1, N(ng)
!          cff1 = dx*dy*PLUME(ng) % dz(is, K)
!          cff2 = dt(ng)*PLUME(ng) % mAm(is, K)
!            t(Isrc, Jsrc, K, nstp, isalt) =               &
!              &  (cff1*t(Isrc, Jsrc, K, nstp, isalt) +    &
!              &   cff2*PLUME(ng) % sB(is, K)) /           &
!              &  (cff1 + cff2)
!            t(Isrc, Jsrc, K, nstp, itemp) =               &
!              &  (cff1*t(Isrc, Jsrc, K, nstp, itemp) +    &
!              &   cff2*PLUME(ng) % tB(is, K)) /           &
!              &  (cff1 + cff2)
!            DO itrc = 3, NT(ng)
!              t(Isrc, Jsrc, K, nstp, itrc) =                &
!                &  (cff1*t(Isrc, Jsrc, K, nstp, itrc)) /    &
!                &  (cff1 + cff2)
!            ENDDO
!            IF (useMeltTracers) THEN
!              SOURCES(ng) % Tsrc(is, K, NT(ng)) = 0.d0
!              t(Isrc, Jsrc, K, nstp, NT(ng)) =              &
!                &  (cff1*t(Isrc, Jsrc, K, nstp, NT(ng)) +   &
!                &   cff2*PLUME(ng) % trc(is, NT(ng))) /     &
!                &  (cff1 + cff2)
!!
!! Hard code meltwater concentration to at least 0.
!!
!              t(Isrc, Jsrc, K, nstp, NT(ng)) = &
!                MAX(t(Isrc, Jsrc, K, nstp, NT(ng)), 0.0d0)
!            ENDIF
!          ENDDO
!        ENDIF
!
      ENDIF
    ENDDO
  ENDIF
END SUBROUTINE iceplume_tile
#endif
END MODULE iceplume_mod
