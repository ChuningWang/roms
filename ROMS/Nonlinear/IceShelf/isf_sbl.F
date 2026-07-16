#include "cppdefs.h"
      MODULE isf_sbl_mod
#if defined ICESHELF_THERMO && !defined ICESHELF_NO_SBL
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2024 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This routine computes the sub-boundary layer (SBL) vertical weights !
!  used to distribute ice shelf basal fluxes over a prescribed         !
!  boundary layer thickness (isf_Hsbl) instead of injecting all flux   !
!  into the topmost cell.                                              !
!                                                                      !
!  The weight Oi(i,j,k,ioHsbl) represents the fraction of the total    !
!  ice shelf flux assigned to vertical layer k.  Weights are           !
!  proportional to cell thickness within the boundary layer and        !
!  normalized so that the vertical sum equals 1 for each water         !
!  column.                                                             !
!                                                                      !
!  Three cases are handled:                                            !
!    1. isf_Hsbl <= 0:   fallback to top-cell only (Oi=1 at k=N, 0     !
!                        elsewhere), equivalent to original behavior.  !
!    2. Hwater <= Hsbl:  entire water column is within the boundary    !
!                        layer; weights proportional to Hz.            !
!    3. Hwater >  Hsbl:  boundary layer extends from the ice base      !
!                        (k=N) downward by isf_Hsbl meters; cells      !
!                        fully within the layer get full Hz, the       !
!                        straddling cell gets the remaining fraction.  !
!                                                                      !
!  ROMS vertical indexing: k=N(ng) is the topmost cell (adjacent to    !
!  the ice base), k=1 is the bottommost cell.                          !
!                                                                      !
!  For static ice shelf configurations, Oi depends only on geometry    !
!  (h, zice, Hz) and is computed once on the first time step.  For     !
!  dynamic ice sheet configurations, this routine must be called       !
!  every baroclinic step to update the weights as geometry evolves.    !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: isf_sbl
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE isf_sbl (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_isf
!
      integer, intent(in) :: ng, tile
!
# include "tile.h"
!
      CALL isf_sbl_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % zice,                               &
     &                   GRID(ng) % h,                                  &
     &                   GRID(ng) % rmask_isf,                          &
     &                   ISF(ng) % Oi)
!
      RETURN
      END SUBROUTINE isf_sbl
!
!***********************************************************************
      SUBROUTINE isf_sbl_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Hz, zice, h, rmask_isf, Oi)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_isf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: zice(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: rmask_isf(LBi:,LBj:)
      real(r8), intent(out) :: Oi(LBi:,LBj:,:,:)
# else
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask_isf(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Oi(LBi:UBi,LBj:UBj,N(ng),nIsfO)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8) :: cff, depth_remaining, Hsbl, Hwater
!
# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute SBL vertical weights.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (rmask_isf(i,j).eq.1.0_r8) THEN
            Hsbl = isf_Hsbl(ng)
            Hwater = h(i,j) + zice(i,j)
!
!  Hsbl <= 0: fallback to top-cell only (weight=1 at k=N, 0 elsewhere).
!  This reproduces the original single-cell flux injection behavior.
!
            IF (Hsbl.le.0.0_r8) THEN
              DO k=1,N(ng)
                Oi(i,j,k,ioHsbl) = 0.0_r8
              END DO
              Oi(i,j,N(ng),ioHsbl) = 1.0_r8
!
!  Water column thinner than Hsbl: distribute over entire column
!  proportional to cell thickness.  This ensures uniform flux
!  distribution when the cavity is shallower than the boundary layer.
!
            ELSE IF (Hwater.le.Hsbl) THEN
              cff = 0.0_r8
              DO k=1,N(ng)
                Oi(i,j,k,ioHsbl) = Hz(i,j,k)
                cff = cff + Hz(i,j,k)
              END DO
              IF (cff.gt.0.0_r8) THEN
                DO k=1,N(ng)
                  Oi(i,j,k,ioHsbl) = Oi(i,j,k,ioHsbl) / cff
                END DO
              ELSE
                DO k=1,N(ng)
                  Oi(i,j,k,ioHsbl) = 0.0_r8
                END DO
                Oi(i,j,N(ng),ioHsbl) = 1.0_r8
              END IF
!
!  Normal case: boundary layer is thinner than the water column.
!  Accumulate cell thicknesses from the top (k=N) downward until
!  the accumulated thickness reaches Hsbl.  Cells fully within the
!  boundary layer receive their full Hz as weight; the straddling
!  cell receives only the remaining thickness (depth_remaining).
!  Cells below the boundary layer receive zero weight.
!
            ELSE
              depth_remaining = Hsbl
              cff = 0.0_r8
              DO k=N(ng),1,-1
                IF (depth_remaining.ge.Hz(i,j,k)) THEN
                  Oi(i,j,k,ioHsbl) = Hz(i,j,k)
                  depth_remaining = depth_remaining - Hz(i,j,k)
                ELSE IF (depth_remaining.gt.0.0_r8) THEN
                  Oi(i,j,k,ioHsbl) = depth_remaining
                  depth_remaining = 0.0_r8
                ELSE
                  Oi(i,j,k,ioHsbl) = 0.0_r8
                END IF
                cff = cff + Oi(i,j,k,ioHsbl)
              END DO
              IF (cff.gt.0.0_r8) THEN
                DO k=1,N(ng)
                  Oi(i,j,k,ioHsbl) = Oi(i,j,k,ioHsbl) / cff
                END DO
              ELSE
!  Degenerate column: assign all weight to top cell
                DO k=1,N(ng)
                  Oi(i,j,k,ioHsbl) = 0.0_r8
                END DO
                Oi(i,j,N(ng),ioHsbl) = 1.0_r8
              END IF
            END IF
          ELSE
!  Non-ice-shelf cells: zero weight (no flux redistribution)
            DO k=1,N(ng)
              Oi(i,j,k,ioHsbl) = 0.0_r8
            END DO
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE isf_sbl_tile
#endif
      END MODULE isf_sbl_mod
