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
# ifdef WET_DRY
!                                                                      !
!  If wetting and drying, masks array are time dependent and changed   !
!  at every time-step.                                                 !
# endif
!                                                                      !
!  Under ICESHEET_SSA, the grounding line is computed dynamically      !
!  from the flotation criterion:                                       !
!    grounded:  rho_i * H > rho_w * h                                  !
!    floating:  rho_i * H <= rho_w * h  (and H > 0)                    !
!                                                                      !
!=======================================================================
!
# ifdef WET_DRY
      USE wetdry_mod,         ONLY : wetdry_mask_tile
# endif
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
# ifdef ICESHEET_SSA
      USE mod_isf, ONLY : ISF, lsnew
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
# ifdef ICESHEET_SSA
     &                     ISF(ng) % Si(:,:,lsnew(ng),1),               &
# endif
# ifdef WET_DRY
     &                     GRID(ng) % pmask_wet,                        &
     &                     GRID(ng) % rmask_wet,                        &
     &                     GRID(ng) % umask_wet,                        &
     &                     GRID(ng) % vmask_wet,                        &
# elif defined MASKING
     &                     GRID(ng) % pmask,                            &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
# endif
# ifdef WET_DRY
     &                     GRID(ng) % pmask_ist,                        &
     &                     GRID(ng) % rmask_ist,                        &
     &                     GRID(ng) % umask_ist,                        &
     &                     GRID(ng) % vmask_ist,                        &
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
# ifdef ICESHEET_SSA
     &                           H_ice,                                 &
# endif
# if defined WET_DRY || defined MASKING
     &                           pmask, rmask,                          &
     &                           umask, vmask,                          &
# endif
# ifdef WET_DRY
     &                           pmask_ist, rmask_ist,                  &
     &                           umask_ist, vmask_ist,                  &
# endif
     &                           rmask_isf)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sources
# ifdef ICESHEET_SSA
      USE mod_isf, ONLY : rho_i, min_cavity_thickness
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
      real(r8), intent(inout) :: zice(LBi:,LBj:)
#  ifdef ICESHEET_SSA
      real(r8), intent(in) :: H_ice(LBi:,LBj:)
#  endif
#  if defined WET_DRY || defined MASKING
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(out) :: pmask_ist(LBi:,LBj:)
      real(r8), intent(out) :: rmask_ist(LBi:,LBj:)
      real(r8), intent(out) :: umask_ist(LBi:,LBj:)
      real(r8), intent(out) :: vmask_ist(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: rmask_isf(LBi:,LBj:)
# else
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: zice(LBi:UBi,LBj:UBj)
#  ifdef ICESHEET_SSA
      real(r8), intent(in) :: H_ice(LBi:UBi,LBj:UBj)
#  endif
#  if defined WET_DRY || defined MASKING
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(out) :: pmask_ist(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rmask_ist(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: umask_ist(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vmask_ist(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: rmask_isf(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
!
# ifdef ICESHEET_SSA
      real(r8) :: H_rho, rho_w_local
# endif
# ifdef WET_DRY
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wetdry
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set internal full ice shelf/sheet masks.
!-----------------------------------------------------------------------
!
# ifdef ICESHEET_SSA
!  Dynamic grounding line: compute masks from flotation criterion.
!    grounded:  rho_i * H > rho_w * h  (ice fills ocean column)
!    floating:  rho_i * H <= rho_w * h  (ice shelf cavity exists)
!    ice-free:  H <= 0
!
!  rmask_isf = 1 where floating ice (cavity exists)
!  rmask_ist = 1 where grounded ice (under WET_DRY)
!
      rho_w_local=REAL(rho0,r8)
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          H_rho=H_ice(i,j)
          IF (H_rho.gt.0.0_r8) THEN
            IF (rho_i*H_rho.gt.rho_w_local*h(i,j)) THEN
!  Grounded ice: no cavity, zice=0.
              rmask_isf(i,j)=0.0_r8
              zice(i,j)=0.0_r8
            ELSE
!  Floating ice: cavity exists, zice from hydrostatic equilibrium.
!  Enforce minimum cavity thickness to avoid numerical instability
!  in very thin water columns near the grounding line.
              rmask_isf(i,j)=1.0_r8
              zice(i,j)=-(rho_i/rho_w_local)*H_rho
              IF (h(i,j)+zice(i,j).lt.min_cavity_thickness(ng)) THEN
                rmask_isf(i,j)=0.0_r8
                zice(i,j)=0.0_r8
              END IF
            END IF
          ELSE
!  Ice-free: use original zice criterion (static ice shelf).
            IF (zice(i,j).ne.0.0_r8) THEN
              rmask_isf(i,j)=1.0_r8
            ELSE
              rmask_isf(i,j)=0.0_r8
            END IF
          END IF
#  if defined WET_DRY || defined MASKING
          rmask_isf(i,j)=rmask_isf(i,j)*rmask(i,j)
#  endif
        END DO
      END DO
# else
!  Static ice shelf: use zice from grid file.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          IF (zice(i,j).ne.0.0_r8) THEN
            rmask_isf(i,j)=1.0_r8
          ELSE
            rmask_isf(i,j)=0.0_r8
          END IF
#  if defined WET_DRY || defined MASKING
          rmask_isf(i,j)=rmask_isf(i,j)*rmask(i,j)
#  endif
        END DO
      END DO
# endif
!
!  Exchange halo for rmask_isf and zice BEFORE wet/dry computation.
!  This ensures rmask_isf and zice are defined on the extended range
!  (Jstr-1:JendR, Istr-1:IendR) needed by wetdry_mask_tile.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_isf)
# ifdef ICESHEET_SSA
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          zice)
# endif
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rmask_isf)
#  ifdef ICESHEET_SSA
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    zice)
#  endif
# endif
# ifdef WET_DRY
!
!  Use rmask_isf (computed above) for wet/dry mask.
!  wetdry=1 where ice shelf cavity exists (floating ice).
!
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
          wetdry(i,j)=rmask_isf(i,j)
        END DO
      END DO
!
! Compute wet/dry mask arrays.
!
      CALL wetdry_mask_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       wetdry,                                    &
     &                       pmask_ist, rmask_ist,                      &
     &                       umask_ist, vmask_ist)
# endif
!
!  Exchange halo for wet/dry mask arrays.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
# ifdef WET_DRY
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask_ist)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_ist)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask_ist)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask_ist)
# endif
      END IF
!
# ifdef DISTRIBUTE
#  ifdef WET_DRY
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pmask_ist, rmask_ist, umask_ist, vmask_ist)
#  endif
# endif

      RETURN
      END SUBROUTINE isf_masks_tile
#endif
      END MODULE isf_masks_mod
