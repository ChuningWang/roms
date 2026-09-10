#include "cppdefs.h"
      MODULE isf_calving_mod
#ifdef ICESHEET_CALVING
!
!git $Id$
!====================================================== Chuning Wang ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module implements pre-defined mask calving for the ice         !
!  sheet/shelf model.  A static mask (GRID%calving_mask, read from     !
!  the grid NetCDF file) defines the region where ice is allowed       !
!  to exist.  The ice shelf flows outward normally, but any ice        !
!  volume advected into a cell with mask = 0 calves immediately:       !
!  thickness, enthalpy and temperature are set to zero there.          !
!                                                                      !
!  Semantics:                                                          !
!    calving_mask = 0  ->  H = 0, E = 0, T = 0 (ice removed)           !
!    calving_mask != 0 ->  no action                                   !
!                                                                      !
!  The calving front is therefore pinned at the mask boundary.         !
!  Velocity is NOT touched: the next SSA solve relaxes it to the       !
!  (vanishing) driving stress.  Cavity geometry (zice, rmask_isf)      !
!  is NOT touched here: isf_masks recomputes it from the flotation     !
!  criterion and handles the H = 0 (ice-free) branch afterwards.       !
!                                                                      !
!  This is the first-stage deterministic scheme.  Future, more         !
!  physical parameterizations (stress/strain-rate or ocean-forced      !
!  calving) can be added as sibling routines selected by CPP flags.    !
!  The previous thickness-threshold scheme (H_calve) is archived in    !
!  .agents/attic/isf_calving_threshold.F.                              !
!                                                                      !
!  Reference:                                                          !
!    MISOMIP2 protocol (De Rydt et al., 2024): fixed calving front     !
!    via pre-defined mask is a standard coupled-experiment setup.      !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
      PRIVATE
      PUBLIC :: isf_calving_mask
      PUBLIC :: isf_calving_mask_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE isf_calving_mask (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_stepping, ONLY : lsnew
      USE mod_isf,      ONLY : ISF, isIhis, isIeis, ifItis
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
      CALL isf_calving_mask_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            GRID(ng) % calving_mask,              &
     &                            ISF(ng) % Si(:,:,lsnew(ng),isIhis),   &
     &                            ISF(ng) % Si(:,:,lsnew(ng),isIeis),   &
     &                            ISF(ng) % Fi(:,:,ifItis))
# ifdef PROFILE
      CALL wclock_off (ng, model, 44, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE isf_calving_mask
!
!***********************************************************************
      SUBROUTINE isf_calving_mask_tile (ng, tile, model,                &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
     &                                  calving_mask,                   &
     &                                  h_ice, e_ice, t_ice)
!***********************************************************************
!
!  Applies the pre-defined calving mask to the current time level
!  (lsnew) of the ice state.  Called from iceshelf.F after thickness
!  advection, every SSA coupling cycle:
!
!    calving_mask = 0  ->  h_ice = 0, e_ice = 0, t_ice = 0
!    calving_mask != 0 ->  no action
!
!  Zeroing enthalpy and temperature together with thickness keeps the
!  E = H*T state consistent (a stale E would pollute T = E/H when ice
!  re-advects into the cell).
!
!  Imported variable declarations.
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: calving_mask(LBi:,LBj:)
      real(r8), intent(inout) :: h_ice(LBi:,LBj:)
      real(r8), intent(inout) :: e_ice(LBi:,LBj:)
      real(r8), intent(inout) :: t_ice(LBi:,LBj:)
# else
      real(r8), intent(in) :: calving_mask(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: h_ice(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: e_ice(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: t_ice(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
!
# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Remove ice where the pre-defined calving mask forbids it.
!-----------------------------------------------------------------------
!
!  The loop covers IstrT:IendT, JstrT:JendT (interior + ghost points)
!  so boundary ghost cells are masked as well; halo exchange below
!  restores inter-tile consistency.  The 0.5 tolerance accepts integer
!  or scaled 0/1 input read from the grid file.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          IF (calving_mask(i,j).lt.0.5_r8) THEN
            h_ice(i,j)=0.0_r8
            e_ice(i,j)=0.0_r8
            t_ice(i,j)=0.0_r8
          END IF
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          h_ice)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          e_ice)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          t_ice)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 3,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    h_ice, e_ice, t_ice)
# endif
!
      RETURN
      END SUBROUTINE isf_calving_mask_tile
#endif
      END MODULE isf_calving_mod
