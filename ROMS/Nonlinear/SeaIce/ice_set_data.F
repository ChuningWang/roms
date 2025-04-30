#include "cppdefs.h"

      MODULE ice_set_data_mod

#if defined ICE_MODEL && defined ICE_FIXED
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2024 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                            Chuning Wang       !
!================================================== Hernan G. Arango ===
!                                                                      !
!  Set sea ice state variables. This module replaces the dynamic sea   !
!  ice model; sea ice thermodynamics and bulk stress terms are still   !
!  calculated.                                                         !
!                                                                      !
!======================================================================!
!
      USE mod_param
      USE mod_forces
      USE mod_ice
      USE mod_scalars
!
      USE ice_bc2d_mod,    ONLY : ice_bc2d_tile
      USE ice_tibc_mod,    ONLY : ice_tibc_tile
      USE ice_uibc_mod,    ONLY : ice_uibc_tile
      USE ice_vibc_mod,    ONLY : ice_vibc_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
!
      implicit none
!
      PRIVATE
      PUBLIC  ice_set_data
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_set_data (ng, tile, model)
!***********************************************************************
!
      USE mod_stepping
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
      CALL ice_set_data_tile (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        liold(ng), linew(ng),                     &
     &                        FORCES(ng) % aice,                        &
     &                        FORCES(ng) % hice,                        &
     &                        FORCES(ng) % hsno,                        &
     &                        FORCES(ng) % hmel,                        &
     &                        FORCES(ng) % tice,                        &
     &                        FORCES(ng) % uice,                        &
     &                        FORCES(ng) % vice,                        &
     &                        ICE(ng) % Si)
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE ice_set_data
!
!***********************************************************************
      SUBROUTINE ice_set_data_tile (ng, tile, model,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              liold, linew,                       &
     &                              aice, hice, hsno, hmel, tice,       &
     &                              uice, vice,                         &
     &                              Si)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(inout) :: liold, linew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: aice(LBi:,LBj:)
      real(r8), intent(in) :: hice(LBi:,LBj:)
      real(r8), intent(in) :: hsno(LBi:,LBj:)
      real(r8), intent(in) :: hmel(LBi:,LBj:)
      real(r8), intent(in) :: tice(LBi:,LBj:)
      real(r8), intent(in) :: uice(LBi:,LBj:)
      real(r8), intent(in) :: vice(LBi:,LBj:)
      real(r8), intent(out) :: Si(LBi:,LBj:,:,:)
# else
      real(r8), intent(in) :: aice(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: hice(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: hsno(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: hmel(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tice(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: uice(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vice(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Si(LBi:UBi,LBj:UBj,2,nIceS)
# endif
!
! Local variable definitions
!
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Impose positive concentration and thicknesses.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Si(i,j,liold,isAice)=aice(i,j)
          Si(i,j,liold,isHice)=hice(i,j)
          Si(i,j,liold,isHsno)=hsno(i,j)
          Si(i,j,liold,isHmel)=hmel(i,j)
          Si(i,j,liold,isTice)=tice(i,j)
          Si(i,j,liold,isEnth)=tice(i,j)*hice(i,j)
          Si(i,j,linew,isAice)=aice(i,j)
          Si(i,j,linew,isHice)=hice(i,j)
          Si(i,j,linew,isHsno)=hsno(i,j)
          Si(i,j,linew,isHmel)=hmel(i,j)
          Si(i,j,linew,isTice)=tice(i,j)
          Si(i,j,linew,isEnth)=tice(i,j)*hice(i,j)
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          Si(i,j,liold,isUice)=uice(i,j)
          Si(i,j,linew,isUice)=uice(i,j)
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          Si(i,j,liold,isVice)=vice(i,j)
          Si(i,j,linew,isVice)=vice(i,j)
        END DO
      END DO
!
!  Set lateral boundary conditions.
!
      CALL ice_bc2d_tile (ng, tile, model, isAice,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, liold,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isAice),                             &
     &                    LBC(:,ibICE(isAice),ng))
!
      CALL ice_bc2d_tile (ng, tile, model, isHice,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, liold,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isHice),                             &
     &                    LBC(:,ibICE(isHice),ng))
!
      CALL ice_bc2d_tile (ng, tile, model, isHsno,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, liold,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isHsno),                             &
     &                    LBC(:,ibICE(isHsno),ng))
!
      CALL ice_bc2d_tile (ng, tile, model, isHmel,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, liold,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isHmel),                             &
     &                    LBC(:,ibICE(isHmel),ng))
!
      CALL ice_tibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    liold, liold,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isHice),                             &
     &                    Si(:,:,:,isTice),                             &
     &                    Si(:,:,:,isEnth))
!
      CALL ice_uibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, liold,                                 &
     &                    Si(:,:,:,isUice))
!
      CALL ice_vibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, liold,                                 &
     &                    Si(:,:,:,isVice))
!
      CALL ice_bc2d_tile (ng, tile, model, isAice,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isAice),                             &
     &                    LBC(:,ibICE(isAice),ng))
!
      CALL ice_bc2d_tile (ng, tile, model, isHice,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isHice),                             &
     &                    LBC(:,ibICE(isHice),ng))
!
      CALL ice_bc2d_tile (ng, tile, model, isHsno,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isHsno),                             &
     &                    LBC(:,ibICE(isHsno),ng))
!
      CALL ice_bc2d_tile (ng, tile, model, isHmel,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isHmel),                             &
     &                    LBC(:,ibICE(isHmel),ng))
!
      CALL ice_tibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isHice),                             &
     &                    Si(:,:,:,isTice),                             &
     &                    Si(:,:,:,isEnth))
!
      CALL ice_uibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,isUice))
!
      CALL ice_vibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,isVice))
# ifdef DISTRIBUTE
!
      CALL mp_exchange3d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj, 1, 2,                     &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Si(:,:,:,isAice),                             &
     &                    Si(:,:,:,isHice),                             &
     &                    Si(:,:,:,isHsno),                             &
     &                    Si(:,:,:,isHmel))
!
      CALL mp_exchange3d (ng, tile, model, 3,                           &
     &                    LBi, UBi, LBj, UBj, 1, 2,                     &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Si(:,:,:,isTice),                             &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice))
# endif
!
      RETURN
      END SUBROUTINE ice_set_data_tile
#endif
      END MODULE ice_set_data_mod
