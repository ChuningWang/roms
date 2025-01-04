#include "cppdefs.h"
      MODULE isf_mom_mod
#ifdef ICESHELF_MOM
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2024 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!========================================== Benjamin K. Galton-Fenzi ===
!                                                                      !
! This module computes iceshelf momemtum equation                      !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: isf_mom
!
      CONTAINS

!
!***********************************************************************
      SUBROUTINE isf_mom (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_ocean
      USE mod_stepping
      USE mod_isf, ONLY : ISF
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
      CALL wclock_on (ng, model, 6, __LINE__, MyFile)
# endif
      CALL isf_mom_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nrhs(ng),                                      &
     &                   ISF(ng) % Si)
# ifdef PROFILE
      CALL wclock_off (ng, model, 6, __LINE__, MyFile)
# endif
      RETURN
      END SUBROUTINE isf_mom
!
!***********************************************************************
      SUBROUTINE isf_mom_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs,                                    &
     &                         Si)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_isf
!
      USE bc_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: Si(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: Si(LBi:UBi,LBj:UBj,2,nIsfS)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
# include "set_bounds.h"
!
!  TBC
!
      RETURN
      END SUBROUTINE isf_mom_tile
#endif ICESHELF_MOM
      END MODULE isf_mom_mod
