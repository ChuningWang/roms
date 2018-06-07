#include "cppdefs.h"
MODULE iceplume_mod
#ifdef ICEPLUME
!
implicit none
!
PRIVATE
PUBLIC :: iceplume
!
CONTAINS
!
SUBROUTINE iceplume (ng, tile)
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
# ifdef MASKING
    &                 GRID(ng) % rmask,                            &
# endif
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

SUBROUTINE iceplume_tile(ng, tile,                                 &
    &                    LBi, UBi, LBj, UBj,                       &
    &                    IminS, ImaxS, JminS, JmaxS,               &
    &                    nrhs,                                     &
# ifdef MASKING
    &                    rmask,                                    &
# endif
    &                    om_r, on_r, h, z_w,                       &
    &                    u, v, wvel, t)
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
#  ifdef MASKING
  real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
  real(r8), intent(in) :: om_r(LBi:,LBj:)
  real(r8), intent(in) :: on_r(LBi:,LBj:)
  real(r8), intent(in) :: h(LBi:,LBj:)
  real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
  real(r8), intent(in) :: u(LBi:,LBj:,:,:)
  real(r8), intent(in) :: v(LBi:,LBj:,:,:)
  real(r8), intent(in) :: wvel(LBi:,LBj:,:)
  real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
# else
#  ifdef MASKING
  real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
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
  real(r8) :: pr, prRef
!
# include "set_bounds.h"
!
! ==================================================================
! Read in some scalar parameters
  ngr = ng
  Nr = N(ng)
!   iceDepth = -abs(GRID(ng) % h(II, JI))
! If checkCFL is activited, need to read in slow-mode time step
  IF (checkCFL) THEN
    dtn = dt(ng)
  ENDIF
!
  IF (LuvSrc(ng)) THEN
    DO I=1,Nsrc(ng)
      II=SOURCES(ng)%Isrc(I)
      JI=SOURCES(ng)%Jsrc(I)
      IF (((IstrR.le.II).and.(II.le.IendR)).and. &
        & ((JstrR.le.JI).and.(JI.le.JendR))) THEN
!
! ==================================================================
! Find the right Rho coordinates
! ==================================================================
!
        DIR=SOURCES(ng)%Dsrc(I)
        IF (DIR .EQ. 0) THEN
          IF ((rmask(II  , JI) .EQ. 0) .AND. &
            & (rmask(II-1, JI) .EQ. 1)) THEN
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

        IF ((SOURCES(ng) % Dsrc(I)) .EQ. 0) THEN  ! read in v direction
          dx = abs(om_r(II, JI))
          dy = abs(on_r(II, JI))
        ELSE
          dx = abs(on_r(II, JI))
          dy = abs(om_r(II, JI))
        ENDIF
!
! ==================================================================
! Read in profiles from OCEAN.
! ==================================================================
!
        DO K = 1, N(ng)+1
          PLUME(ng) % zW(K) = z_w(II, JI, K-1)
        ENDDO
!
        DO K = 1, N(ng)
! temp and salt
!     PLUME(ng) % tAm(K) = OCEAN(ng) % t(II, JI, K, nrhs, itemp)
          PLUME(ng) % sAm(K) = t(II, JI, K, nrhs, isalt)
          PLUME(ng) % tpAm(K) = t(II, JI, K, nrhs, itemp)
! convert potential temp to in-situ temp
          prRef = 101.d3*1.d-4
          pr = prRef + &
            & (abs(PLUME(ng) % zW(K))*rho_ref*g)*1.d-4  ! [dbar]
          CALL SW_TEMP(PLUME(ng) % sAm(K),  &
                     & PLUME(ng) % tpAm(K), &
                     & pr, prRef, PLUME(ng) % tAm(K))
! tracers
          IF (useTracers) THEN
            DO iTracer = 3, NT(ng)
              PLUME(ng) % trcAm(K, iTracer) = &
                & t(II, JI, K, nrhs, iTracer)
            ENDDO
          ENDIF
! u/v, w
          IF (DIR .EQ. 0) THEN
            PLUME(ng) % vAm(K) = &
              & 0.25d0 * ( &
              & v(II, JI,   K, nrhs) + &
              & v(II, JI+1, K, nrhs) &
              & )
          ELSEIF (DIR .EQ. 1) THEN
            PLUME(ng) % vAm(K) = &
              & 0.25d0 * ( &
              & u(II,   JI, K, nrhs) + &
              & u(II+1, JI, K, nrhs) &
              & )
          ENDIF

          PLUME(ng) % wAm(K) = &
            & 0.5 * ( &
            & wvel(II, JI, K-1) + &
            & wvel(II, JI, K) &
            & )
        ENDDO
! Discharge tracer concentration
        IF (useTracers .and. useInputTracers) THEN
          DO iTracer = 3, NT(ng)
            PLUME(ng) % trcIni(iTracer) = &
              & SOURCES(ng) % Tsrc(I, 1, iTracer)
          ENDDO
        ELSE
          DO iTracer = 3, NT(ng)
            PLUME(ng) % trcIni(iTracer) = 0.d0
          ENDDO
        ENDIF

! ==================================================================
        CALL iceplume_calc(ng, &
                         & abs(SOURCES(ng) % Qbar(I)),      &
                         & SOURCES(ng) % Tsrc(I, 1, itemp), &
                         & SOURCES(ng) % Tsrc(I, 1, isalt))
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
!
! ==================================================================
        DO K = 1, Nr
! Qsrc
          IF (SOURCES(ng) % Qbar(I) .GT. 0.d0) THEN
            SOURCES(ng) % Qsrc(I, K) = &
              & (PLUME(ng) % ent(K) + PLUME(ng) % det(K))
          ELSEIF (SOURCES(ng) % Qbar(I) .LT. 0.d0) THEN
            SOURCES(ng) % Qsrc(I, K) = &
              & -(PLUME(ng) % ent(K) + PLUME(ng) % det(K))
          ELSE
            SOURCES(ng) % Qsrc(I, K) = 0.d0
          ENDIF
! T, S, Tracer concentration
          IF (SOURCES(ng) % Qsrc(I, K) .NE. 0.d0) THEN
            SOURCES(ng) % Tsrc(I, K, isalt) = &
              & (PLUME(ng) % det(K)*PLUME(ng) % trc(isalt) + &
              &  PLUME(ng) % ent(K)*PLUME(ng) % sAm(K)) / &
              & (PLUME(ng) % det(K)+PLUME(ng) % ent(K))
            SOURCES(ng) % Tsrc(I, K, itemp) = &
              & (PLUME(ng) % det(K)*PLUME(ng) % trc(itemp) + &
              &  PLUME(ng) % ent(K)*PLUME(ng) % tAm(K)) / &
              & (PLUME(ng) % det(K)+PLUME(ng) % ent(K))
            DO iTracer = 3, NT(ng)
              SOURCES(ng) % Tsrc(I, K, iTracer) = &
                & (PLUME(ng) % det(K)*PLUME(ng) % trc(iTracer) + &
                &  PLUME(ng) % ent(K)*PLUME(ng) % trcAm(K,iTracer)) / &
                & (PLUME(ng) % det(K)+PLUME(ng) % ent(K))
            ENDDO
          ELSE
            SOURCES(ng) % Tsrc(I, K, isalt) = 0.d0
            SOURCES(ng) % Tsrc(I, K, itemp) = 0.d0
            DO iTracer = 3, NT(ng)
              SOURCES(ng) % Tsrc(I, K, iTracer) = 0.d0
            ENDDO
          ENDIF
        ENDDO
! Use the tendency terms to update T & S in OCEAN
        DO K = 1, Nr
          t(II, JI, K, nrhs, isalt) = &
            & t(II, JI, K, nrhs, isalt) + &
            & dt(ng) * PLUME(ng) % tendS(K)
          t(II, JI, K, nrhs, itemp) = &
            & t(II, JI, K, nrhs, itemp) + &
            & dt(ng) * PLUME(ng) % tendT(K)
        ENDDO

        !DO K = 1, Nr
        !  write(*, *)  K, PLUME(ng) % sAm(K), &
        !    & PLUME(ng) % ent(K), &
        !    & PLUME(ng) % det(K), &
        !    & SOURCES(ng) % Tsrc(I, K, isalt)
        !ENDDO

! Write to file for diagonse
        IF (useDebug) THEN
          ! IF (iic(ng) .EQ. ntimes(ng)) THEN
          IF (iic(ng) .EQ. 1) THEN
          ! IF (mod(icc, 120) .EQ. 0)
            open(unit=15, file='./plume_debug.txt')
            write(15, '(I3)')  I
            write(15, '(99 E12.4)')  PLUME(ng) % zW
            write(15, '(99 E12.4)')  PLUME(ng) % tAm
            write(15, '(99 E12.4)')  PLUME(ng) % sAm
            write(15, '(99 E12.4)')  PLUME(ng) % vAm
            write(15, '(99 E12.4)')  PLUME(ng) % wAm
            write(15, *)  ' '
            write(15, '(99 E12.4)')  PLUME(ng) % ent
            write(15, '(99 E12.4)')  PLUME(ng) % s
            write(15, '(99 E12.4)')  PLUME(ng) % t
            write(15, '(99 E12.4)')  PLUME(ng) % r
            write(15, '(99 E12.4)')  PLUME(ng) % a
            write(15, '(99 E12.4)')  PLUME(ng) % w
            write(15, '(99 E12.4)')  PLUME(ng) % mInt
            write(15, *)  ' '
            write(15, '(99 E12.4)')  SOURCES(ng) % Tsrc(I, :, 1)
            write(15, '(99 E12.4)')  SOURCES(ng) % Tsrc(I, :, 2)
            write(15, '(99 E12.4)')  SOURCES(ng) % Tsrc(I, :, 3)
            write(15, *)  ' '
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  ENDIF

END SUBROUTINE iceplume_tile

#endif
END MODULE iceplume_mod
