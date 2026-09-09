#include "cppdefs.h"
      MODULE isf_masks_mod
#ifdef ICESHELF
!
!git $Id$
!====================================================== Chuning Wang ===
!  Copyright (c) 2002-2025 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  These routines set internal ice shelf/sheet masking arrays.         !
# ifdef ICESHEET
!                                                                      !
!  If using ice flow model, activate wetting and drying. The grounding !
!  line is computed dynamically from the flotation criterion (using    !
!  mod_scalars rho0 as seawater density):                              !
!    grounded:  rho_i * H >  rho0 * (zeta+h-Dcrit-eps)                 !
!    floating:  rho_i * H <= rho0 * (zeta+h-Dcrit-eps)  (and H > 0)    !
# endif
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: isf_masks
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE isf_masks (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
# ifdef ICESHEET
      USE mod_stepping, ONLY : knew, lsnew
      USE mod_ocean
      USE mod_isf,      ONLY : isIhis, ISF
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 42, __LINE__, MyFile)
# endif
      CALL isf_masks_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     GRID(ng) % h,                                &
     &                     GRID(ng) % zice,                             &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
# ifdef ICESHEET
     &                     OCEAN(ng) % zeta(:,:,knew(ng)),              &
     &                     ISF(ng) % Si(:,:,lsnew(ng),isIhis),          &
     &                     GRID(ng) % rmask_ist,                        &
     &                     GRID(ng) % umask_ist,                        &
     &                     GRID(ng) % vmask_ist,                        &
     &                     GRID(ng) % pmask_ist,                        &
# endif
     &                     GRID(ng) % rmask_isf)
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE isf_masks
!
!***********************************************************************
      SUBROUTINE isf_masks_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           h, zice,                               &
# ifdef MASKING
     &                           rmask,                                 &
# endif
# ifdef ICESHEET
     &                           zeta, hice,                            &
     &                           rmask_ist, umask_ist,                  &
     &                           vmask_ist, pmask_ist,                  &
# endif
     &                           rmask_isf)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
# ifdef ICESHEET
      USE mod_isf,    ONLY : rho_i
# endif
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: h(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
#  ifdef ICESHEET
      real(r8), intent(inout) :: zice(LBi:,LBj:)
      real(r8), intent(in) :: zeta(LBi:,LBj:)
      real(r8), intent(in) :: hice(LBi:,LBj:)
      real(r8), intent(out) :: rmask_ist(LBi:,LBj:)
      real(r8), intent(out) :: umask_ist(LBi:,LBj:)
      real(r8), intent(out) :: vmask_ist(LBi:,LBj:)
      real(r8), intent(out) :: pmask_ist(LBi:,LBj:)
#  else
      real(r8), intent(in) :: zice(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: rmask_isf(LBi:,LBj:)
# else
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ICESHEET
      real(r8), intent(inout) :: zice(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: hice(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rmask_ist(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: umask_ist(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vmask_ist(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: pmask_ist(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: rmask_isf(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8), parameter :: eps = 1.0E-10_r8
!
# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set internal full ice shelf/sheet masks.
!-----------------------------------------------------------------------
!
# ifdef ICESHEET
!  Dynamic grounding line: compute masks from flotation criterion.
!
!    grounded:  rho_i * H > rho_w * h  (ice fills ocean column)
!    floating:  rho_i * H <= rho_w * h  (ice shelf cavity exists)
!    ice-free:  H <= 0
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          IF (hice(i,j).gt.0.0_r8) THEN
            IF (rho_i*hice(i,j).gt.                                     &
     &          rho0*(zeta(i,j)+h(i,j)-Dcrit(ng)-eps)) THEN
!
!  Grounded ice: no cavity. set water column depth to Dcrit(ng)+eps
!
              rmask_isf(i,j)=0.0_r8
              zice(i,j)=-zeta(i,j)-h(i,j)+Dcrit(ng)+eps
            ELSE
!
!  Floating ice: cavity exists, zice from hydrostatic equilibrium.
!
              rmask_isf(i,j)=1.0_r8
              zice(i,j)=-(rho_i/rho0)*hice(i,j)
            END IF
          ELSE
!
!  Ice-free.
!
            rmask_isf(i,j)=0.0_r8
            zice(i,j)=0.0_r8
          END IF
#  ifdef MASKING
          rmask_isf(i,j)=rmask_isf(i,j)*rmask(i,j)
#  endif
        END DO
      END DO
!
! Compute local ice sheet mask at RHO-points in terms of ice sheet
! thickness.
!
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
          IF (hice(i,j).gt.0.0_r8) THEN
            rmask_ist(i,j)=1.0_r8
          ELSE
            rmask_ist(i,j)=0.0_r8
          END IF
        END DO
      END DO
!
!  ice sheet mask at PSI-points.
!
      DO j=Jstr,JendR
        DO i=Istr,IendR
          umask_ist(i,j)=rmask_ist(i-1,j)+rmask_ist(i,j)
          IF (umask_ist(i,j).eq.1.0_r8) THEN
            umask_ist(i,j)=rmask_ist(i-1,j)-rmask_ist(i,j)
          END IF
!
          vmask_ist(i,j)=rmask_ist(i,j-1)+rmask_ist(i,j)
          IF (vmask_ist(i,j).eq.1.0_r8) THEN
            vmask_ist(i,j)=rmask_ist(i,j-1)-rmask_ist(i,j)
          END IF
!
          IF ((rmask_ist(i-1,j-1).eq.1.0_r8).and.                       &
     &        (rmask_ist(i  ,j-1).eq.1.0_r8).and.                       &
     &        (rmask_ist(i-1,j  ).eq.1.0_r8).and.                       &
     &        (rmask_ist(i  ,j  ).eq.1.0_r8)) THEN
            pmask_ist(i,j)=1.0_r8
          ELSE
            pmask_ist(i,j)=0.0_r8
          END IF
        END DO
      END DO
# else
!  Static ice shelf: use zice from grid file.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          IF (ABS(zice(i,j)).gt.eps) THEN
            rmask_isf(i,j)=1.0_r8
          ELSE
            rmask_isf(i,j)=0.0_r8
          END IF
#  ifdef MASKING
          rmask_isf(i,j)=rmask_isf(i,j)*rmask(i,j)
#  endif
        END DO
      END DO
# endif
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_isf)
# ifdef ICESHEET
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          zice)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_ist)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask_ist)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask_ist)
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask_ist)
# endif
      END IF

# ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rmask_isf)
#  ifdef ICESHEET
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    zice)
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rmask_ist, umask_ist, vmask_ist, pmask_ist)
#  endif
# endif

      RETURN
      END SUBROUTINE isf_masks_tile
#endif
      END MODULE isf_masks_mod
