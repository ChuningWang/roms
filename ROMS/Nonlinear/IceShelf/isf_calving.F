#include "cppdefs.h"
      MODULE isf_calving_mod
#ifdef ICESHEET_CALVING
!
!git $Id$
!====================================================== Chuning Wang ===
!  Copyright (c) 2002-2025 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module implements thickness threshold calving for the ice      !
!  sheet/shelf model. Ice at the calving front that is thinner than    !
!  a critical threshold H_calve is removed.                            !
!                                                                      !
!  Reference:                                                          !
!    PISM CalvingAtThickness.cc                                        !
!    MISOMIP1 protocol: H_calve = 100 m (Asay-Davis et al., 2016)      !
!                                                                      !
!  Algorithm:                                                          !
!    1. Detect ice front cells: H > 0 but at least one neighbor H = 0  !
!    2. If H < H_calve at ice front, remove ice (set H = 0)            !
!    3. Update masks and zice accordingly                              !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
      PRIVATE
      PUBLIC :: isf_calving
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE isf_calving (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_isf, ONLY : ISF, lsnew
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
      CALL wclock_on (ng, model, 44, __LINE__, MyFile)
# endif
      CALL isf_calving_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
     &                       GRID(ng) % h,                              &
     &                       GRID(ng) % zice,                           &
     &                       GRID(ng) % rmask_isf,                      &
     &                       ISF(ng) % Si(:,:,lsnew(ng),1))
# ifdef PROFILE
      CALL wclock_off (ng, model, 44, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE isf_calving
!
!***********************************************************************
      SUBROUTINE isf_calving_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             pm, pn, h, zice, rmask_isf, H_ice)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_isf, ONLY : H_calve, rho_i
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
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(inout) :: zice(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_isf(LBi:,LBj:)
      real(r8), intent(inout) :: H_ice(LBi:,LBj:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: zice(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask_isf(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: H_ice(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
      logical :: is_front
      real(r8) :: rho_w_local
!
# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Thickness threshold calving.
!-----------------------------------------------------------------------
!
!  Remove ice at the calving front that is thinner than H_calve.
!  Ice front detection: H > 0 but at least one neighbor has H = 0.
!
!  First, exchange halo for H_ice so neighbor checks work correctly.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          H_ice)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    H_ice)
# endif
!
      rho_w_local=REAL(rho0,r8)
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          IF (H_ice(i,j).gt.0.0_r8) THEN
!  Check if this is an ice front cell (has at least one ice-free neighbor).
            is_front=.FALSE.
            IF (H_ice(i-1,j).le.0.0_r8) is_front=.TRUE.
            IF (H_ice(i+1,j).le.0.0_r8) is_front=.TRUE.
            IF (H_ice(i,j-1).le.0.0_r8) is_front=.TRUE.
            IF (H_ice(i,j+1).le.0.0_r8) is_front=.TRUE.
!
!  Apply calving if at front and thickness below threshold.
!  Only remove ice thickness here.  The mask and zice update is
!  handled by isf_masks (requires WET_DRY) to ensure a consistent
!  transition of the ocean vertical coordinate system.
!
            IF (is_front.and.H_ice(i,j).lt.H_calve(ng)) THEN
              H_ice(i,j)=0.0_r8
# ifdef WET_DRY
              rmask_isf(i,j)=0.0_r8
              zice(i,j)=0.0_r8
# endif
            END IF
          END IF
        END DO
      END DO
!
!  Exchange halo for modified fields.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          H_ice)
# ifdef WET_DRY
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_isf)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          zice)
# endif
      END IF
# ifdef DISTRIBUTE
#  ifdef WET_DRY
      CALL mp_exchange2d (ng, tile, model, 3,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    H_ice, rmask_isf, zice)
#  else
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    H_ice)
#  endif
# endif
!
      RETURN
      END SUBROUTINE isf_calving_tile
#endif
      END MODULE isf_calving_mod
