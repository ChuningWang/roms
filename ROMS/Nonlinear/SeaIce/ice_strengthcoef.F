#include "cppdefs.h"

      MODULE ice_strengthcoef_mod

#ifdef ICE_STRENGTH_QUAD
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2025 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initializes pstar_grid based on the grid spacing.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ice
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
      implicit none
!
      PRIVATE
      PUBLIC  :: ice_strengthcoef
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_strengthcoef (ng, tile, model)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
#include "tile.h"
!
      CALL ice_strengthcoef_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            GRID(ng) % grdscl,                    &
     &                            ICE(ng) % Fi)

      RETURN
      END SUBROUTINE ice_strengthcoef
!
!***********************************************************************
      SUBROUTINE ice_strengthcoef_tile (ng, tile, model,                &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
     &                                  grdscl,                         &
     &                                  Fi)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: grdscl(LBi:,LBj:)
      real(r8), intent(out) :: Fi(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: grdscl(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Fi(LBi:UBi,LBj:UBj,nIceF)
# endif
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax

      integer :: i, j

# include "set_bounds.h"

!
!-----------------------------------------------------------------------
!  Scale ice strength according to the grid size. Scale by 1 km.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Fi(i,j,icPgrd)=Pstar_e(ng)*0.001_r8*grdscl(i,j)
        END DO
      END DO
!
      RETURN
      END SUBROUTINE ice_strengthcoef_tile
#endif
      END MODULE ice_strengthcoef_mod
