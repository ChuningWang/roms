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
     &                     GRID(ng) % zice,                             &
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
     &                           zice,                                  &
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
      real(r8), intent(in) :: zice(LBi:,LBj:)
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
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
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

# ifdef WET_DRY
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wetdry
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set internal full ice shelf/sheet masks.
!-----------------------------------------------------------------------
!
!  Mask at RHO-points.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          IF (zice(i,j).ne.0.0_r8) THEN
            rmask_isf(i,j)=1.0_r8
          ELSE
            rmask_isf(i,j)=0.0_r8
          END IF
# if defined WET_DRY || defined MASKING
          rmask_isf(i,j)=rmask_isf(i,j)*rmask(i,j)
# endif
        END DO
      END DO
# ifdef WET_DRY
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
          IF (zice(i,j).ne.0.0_r8) THEN
            wetdry(i,j)=1.0_r8
          ELSE
            wetdry(i,j)=0.0_r8
          END IF
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
!      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!# ifdef WET_DRY
!        CALL exchange_p2d_tile (ng, tile,                               &
!     &                          LBi, UBi, LBj, UBj,                     &
!     &                          pmask_ist)
!        CALL exchange_r2d_tile (ng, tile,                               &
!     &                          LBi, UBi, LBj, UBj,                     &
!     &                          rmask_ist)
!        CALL exchange_u2d_tile (ng, tile,                               &
!     &                          LBi, UBi, LBj, UBj,                     &
!     &                          umask_ist)
!        CALL exchange_v2d_tile (ng, tile,                               &
!     &                          LBi, UBi, LBj, UBj,                     &
!     &                          vmask_ist)
!# endif
!        CALL exchange_r2d_tile (ng, tile,                               &
!     &                          LBi, UBi, LBj, UBj,                     &
!     &                          rmask_isf)
!      END IF
!
!# ifdef DISTRIBUTE
!#  ifdef WET_DRY
!      CALL mp_exchange2d (ng, tile, model, 4,                           &
!     &                    LBi, UBi, LBj, UBj,                           &
!     &                    NghostPoints,                                 &
!     &                    EWperiodic(ng), NSperiodic(ng),               &
!     &                    pmask_ist, rmask_ist, umask_ist, vmask_ist)
!#  endif
!      CALL mp_exchange2d (ng, tile, model, 1,                           &
!     &                    LBi, UBi, LBj, UBj,                           &
!     &                    NghostPoints,                                 &
!     &                    EWperiodic(ng), NSperiodic(ng),               &
!     &                    rmask_isf)
!# endif

      RETURN
      END SUBROUTINE isf_masks_tile
#endif
      END MODULE isf_masks_mod
