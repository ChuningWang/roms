#include "cppdefs.h"
      MODULE ice_ini_fields_mod
#ifdef NONLINEAR
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2024 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license               Chuning Wang   !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This routine initializes other time levels for ice fields.          !
!                                                                      !
!=======================================================================
!
      USE mod_ice
!
      implicit none
!
      PRIVATE
      PUBLIC :: ice_ini_fields
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_ini_fields (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
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
      CALL wclock_on (ng, iNLM, 2, __LINE__, MyFile)
# endif
      CALL ice_ini_fields_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          liold(ng), liuol(ng),                   &
# ifdef ICE_BIO
     &                          isIphy, isINO3, isINH4,                 &
# endif
     &                          ibICE,                                  &
# ifdef MASKING
     &                          GRID(ng) % rmask,                       &
     &                          GRID(ng) % umask,                       &
     &                          GRID(ng) % vmask,                       &
# endif
# ifdef WET_DRY
     &                          GRID(ng) % rmask_wet,                   &
     &                          GRID(ng) % umask_wet,                   &
     &                          GRID(ng) % vmask_wet,                   &
# endif
     &                          ICE(ng) % Si,                           &
     &                          ICE(ng) % Fi,                           &
     &                          ICE(ng) % Ti)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 2, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE ice_ini_fields
!
!***********************************************************************
      SUBROUTINE ice_ini_fields_tile (ng, tile, model,                  &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                liold, liuol,                     &
# ifdef ICE_BIO
     &                                isIphy, isINO3, isINH4,           &
# endif
     &                                ibICE,                            &
# ifdef MASKING
     &                                rmask, umask, vmask,              &
# endif
# ifdef WET_DRY
     &                                rmask_wet, umask_wet, vmask_wet,  &
# endif
     &                                Si, Fi, Ti)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      USE bc_3d_mod,       ONLY : bc_r3d_tile
      USE ice_uibc_mod,    ONLY : ice_uibc_tile
      USE ice_vibc_mod,    ONLY : ice_vibc_tile
      USE ice_bc2d_mod,    ONLY : ice_bc2d_tile
      USE exchange_3d_mod, ONLY : exchange_r3d_tile, exchange_u3d_tile, &
     &                            exchange_v3d_tile
      USE exchange_4d_mod, ONLY : exchange_r4d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: liold, liuol
# ifdef ICE_BIO
      integer, intent(in) :: isIphy, isINO3, isINH4
# endif
      integer, intent(in) :: ibICE(nIceS)
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: Fi(LBi:,LBj:,:)
      real(r8), intent(inout) :: Si(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: Ti(LBi:,LBj:,:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: Fi(LBi:UBi,LBj:UBj,nIceF)
      real(r8), intent(inout) :: Si(LBi:UBi,LBj:UBj,2,nIceS)
      real(r8), intent(inout) :: Ti(LBi:UBi,LBj:UBj,2,NTICE)
# endif
!
!  Local variable declarations.
!
      integer :: linew, liunw
      integer :: i, j, nf, ic, itrc

      real(r8) :: cff

# include "set_bounds.h"

!
!-----------------------------------------------------------------------
!  If not perfect restart, initialize other time levels for ice state
!  variables.
!-----------------------------------------------------------------------
!
      linew=3-liold
      liunw=3-liuol
      IF (.not.PerfectRST(ng)) THEN
        DO j=JstrB,JendB
          DO i=IstrB,IendB
            cff=Si(i,j,liold,isAice)
# ifdef MASKING
            cff=cff*rmask(i,j)
# endif
# ifdef WET_DRY
            cff=cff*rmask_wet(i,j)
# endif
            Si(i,j,liold,isAice)=cff
            Si(i,j,linew,isAice)=cff
!
            cff=Si(i,j,liold,isHice)
# ifdef MASKING
            cff=cff*rmask(i,j)
# endif
# ifdef WET_DRY
            cff=cff*rmask_wet(i,j)
# endif
            Si(i,j,liold,isHice)=cff
            Si(i,j,linew,isHice)=cff
!
            cff=Si(i,j,liold,isHmel)
# ifdef MASKING
            cff=cff*rmask(i,j)
# endif
# ifdef WET_DRY
            cff=cff*rmask_wet(i,j)
# endif
            Si(i,j,liold,isHmel)=cff
            Si(i,j,linew,isHmel)=cff
!
            cff=Si(i,j,liold,isHsno)
# ifdef MASKING
            cff=cff*rmask(i,j)
# endif
# ifdef WET_DRY
            cff=cff*rmask_wet(i,j)
# endif
            Si(i,j,liold,isHsno)=cff
            Si(i,j,linew,isHsno)=cff
!
            cff=Si(i,j,liold,isIage)
# ifdef MASKING
            cff=cff*rmask(i,j)
# endif
# ifdef WET_DRY
            cff=cff*rmask_wet(i,j)
# endif
            Si(i,j,liold,isIage)=cff
            Si(i,j,linew,isIage)=cff
!
            cff=Si(i,j,liold,isISxx)
# ifdef MASKING
            cff=cff*rmask(i,j)
# endif
# ifdef WET_DRY
            cff=cff*rmask_wet(i,j)
# endif
            Si(i,j,liold,isISxx)=cff
            Si(i,j,linew,isISxx)=cff
!
            cff=Si(i,j,liold,isISxy)
# ifdef MASKING
            cff=cff*rmask(i,j)
# endif
# ifdef WET_DRY
            cff=cff*rmask_wet(i,j)
# endif
            Si(i,j,liold,isISxy)=cff
            Si(i,j,linew,isISxy)=cff
!
            cff=Si(i,j,liold,isISyy)
# ifdef MASKING
            cff=cff*rmask(i,j)
# endif
# ifdef WET_DRY
            cff=cff*rmask_wet(i,j)
# endif
            Si(i,j,liold,isISyy)=cff
            Si(i,j,linew,isISyy)=cff
!
            cff=Si(i,j,liold,isTice)
# ifdef MASKING
            cff=cff*rmask(i,j)
# endif
# ifdef WET_DRY
            cff=cff*rmask_wet(i,j)
# endif
            Si(i,j,liold,isTice)=cff
            Si(i,j,linew,isTice)=cff
# ifdef ICE_BIO
!
            cff=Si(i,j,liold,isIphy)
#  ifdef MASKING
            cff=cff*rmask(i,j)
#  endif
#  ifdef WET_DRY
            cff=cff*rmask_wet(i,j)
#  endif
            Si(i,j,liold,isIphy)=cff
            Si(i,j,linew,isIphy)=cff
!
            cff=Si(i,j,liold,isINO3)
#  ifdef MASKING
            cff=cff*rmask(i,j)
#  endif
#  ifdef WET_DRY
            cff=cff*rmask_wet(i,j)
#  endif
            Si(i,j,liold,isINO3)=cff
            Si(i,j,linew,isINO3)=cff
!
            cff=Si(i,j,liold,isINH4)
#  ifdef MASKING
            cff=cff*rmask(i,j)
#  endif
#  ifdef WET_DRY
            cff=cff*rmask_wet(i,j)
#  endif
            Si(i,j,liold,isINH4)=cff
            Si(i,j,linew,isINH4)=cff
# endif
!
            IF (Licetrc(ng)) THEN
              DO ic=1,NTICE(ng)
                cff=Ti(i,j,liold,ic)
# ifdef MASKING
                cff=cff*rmask(i,j)
# endif
# ifdef WET_DRY
                cff=cff*rmask_wet(i,j)
# endif
                Ti(i,j,liold,ic)=cff
                Ti(i,j,linew,ic)=cff
              END DO
            END IF
          END DO
        END DO
!
        DO j=JstrB,JendB
          DO i=IstrM,IendB
            cff=Si(i,j,liuol,isUice)
# ifdef MASKING
            cff=cff*umask(i,j)
# endif
# ifdef WET_DRY
            cff=cff*umask_wet(i,j)
# endif
            Si(i,j,liuol,isUice)=cff
            Si(i,j,liunw,isUice)=cff
          END DO
        END DO
!
        DO j=JstrM,JendB
          DO i=IstrB,IendB
            cff=Si(i,j,liuol,isVice)
# ifdef MASKING
            cff=cff*vmask(i,j)
# endif
# ifdef WET_DRY
            cff=cff*vmask_wet(i,j)
# endif
            Si(i,j,liuol,isVice)=cff
            Si(i,j,liunw,isVice)=cff
          END DO
        END DO
!
!  Apply boundary conditions.
!
        CALL ice_uibc_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liuol, liuol,                               &
     &                      Si(:,:,:,isUice))
!
        CALL ice_vibc_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liuol, liuol,                               &
     &                      Si(:,:,:,isVice))
!
        CALL ice_uibc_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liuol, liunw,                               &
     &                      Si(:,:,:,isUice))
!
        CALL ice_vibc_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liuol, liunw,                               &
     &                      Si(:,:,:,isVice))
!
        CALL ice_bc2d_tile (ng, tile, model, isAice,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, liold,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isAice),                           &
     &                      LBC(:,ibICE(isAice),ng))
!
        CALL ice_bc2d_tile (ng, tile, model, isHice,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, liold,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isHice),                           &
     &                      LBC(:,ibICE(isHice),ng))
!
        CALL ice_bc2d_tile (ng, tile, model, isHmel,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, liold,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isHmel),                           &
     &                      LBC(:,ibICE(isHmel),ng))
!
        CALL ice_bc2d_tile (ng, tile, model, isHsno,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, liold,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isHsno),                           &
     &                      LBC(:,ibICE(isHsno),ng))
!
        CALL ice_bc2d_tile (ng, tile, model, isIage,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, liold,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isIage),                           &
     &                      LBC(:,ibICE(isIage),ng))
!
        CALL ice_bc2d_tile (ng, tile, model, isISxx,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, liold,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isISxx),                           &
     &                      LBC(:,ibICE(isISxx),ng))
!
        CALL ice_bc2d_tile (ng, tile, model, isISxy,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, liold,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isISxy),                           &
     &                      LBC(:,ibICE(isISxy),ng))
!
        CALL ice_bc2d_tile (ng, tile, model, isISyy,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, liold,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isISyy),                           &
     &                      LBC(:,ibICE(isISyy),ng))
!
        CALL ice_bc2d_tile (ng, tile, model, isTice,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, liold,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isTice),                           &
     &                      LBC(:,ibICE(isTice),ng))
!
# ifdef ICE_BIO
        CALL ice_bc2d_tile (ng, tile, model, isIphy,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, liold,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isIphy),                           &
     &                      LBC(:,ibICE(isIphy),ng))
!
        CALL ice_bc2d_tile (ng, tile, model, isINO3,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, liold,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isINO3),                           &
     &                      LBC(:,ibICE(isINO3),ng))
!
        CALL ice_bc2d_tile (ng, tile, model, isINH4,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, liold,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isINH4),                           &
     &                      LBC(:,ibICE(isINH4),ng))
# endif
!
        IF (Licetrc(ng)) THEN
          CALL bc_r3d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, 1, NTICE(ng),           &
     &                      Ti(:,:,liold,:))
        END IF
!
        CALL ice_bc2d_tile (ng, tile, model, isAice,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, linew,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isAice),                           &
     &                      LBC(:,ibICE(isAice),ng))
!
        CALL ice_bc2d_tile (ng, tile, model, isHice,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, linew,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isHice),                           &
     &                      LBC(:,ibICE(isHice),ng))
!
        CALL ice_bc2d_tile (ng, tile, model, isHmel,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, linew,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isHmel),                           &
     &                      LBC(:,ibICE(isHmel),ng))
!
        CALL ice_bc2d_tile (ng, tile, model, isHsno,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, linew,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isHsno),                           &
     &                      LBC(:,ibICE(isHsno),ng))
!
        CALL ice_bc2d_tile (ng, tile, model, isIage,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, linew,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isIage),                           &
     &                      LBC(:,ibICE(isIage),ng))
!
        CALL ice_bc2d_tile (ng, tile, model, isISxx,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, linew,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isISxx),                           &
     &                      LBC(:,ibICE(isISxx),ng))
!
        CALL ice_bc2d_tile (ng, tile, model, isISxy,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, linew,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isISxy),                           &
     &                      LBC(:,ibICE(isISxy),ng))
!
        CALL ice_bc2d_tile (ng, tile, model, isISyy,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, linew,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isISyy),                           &
     &                      LBC(:,ibICE(isISyy),ng))
!
        CALL ice_bc2d_tile (ng, tile, model, isTice,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, linew,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isTice),                           &
     &                      LBC(:,ibICE(isTice),ng))
!
# ifdef ICE_BIO
        CALL ice_bc2d_tile (ng, tile, model, isIphy,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, linew,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isIphy),                           &
     &                      LBC(:,ibICE(isIphy),ng))
!
        CALL ice_bc2d_tile (ng, tile, model, isINO3,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, linew,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isINO3),                           &
     &                      LBC(:,ibICE(isINO3),ng))
!
        CALL ice_bc2d_tile (ng, tile, model, isINH4,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, linew,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isINH4),                           &
     &                      LBC(:,ibICE(isINH4),ng))
# endif
!
        IF (Licetrc(ng)) THEN
          CALL bc_r3d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, 1, NTICE(ng),           &
     &                      Ti(:,:,linew,:))
        END IF
      END IF
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, 2,               &
     &                          Si(:, :,:,isAice))
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, 2,               &
     &                          Si(:, :,:,isHice))
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, 2,               &
     &                          Si(:, :,:,isHmel))
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, 2,               &
     &                          Si(:, :,:,isHsno))
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, 2,               &
     &                          Si(:, :,:,isIage))
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, 2,               &
     &                          Si(:, :,:,isISxx))
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, 2,               &
     &                          Si(:, :,:,isISxy))
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, 2,               &
     &                          Si(:, :,:,isISyy))
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, 2,               &
     &                          Si(:, :,:,isTice))
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, 2,               &
     &                          Si(:, :,:,isUice))
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, 2,               &
     &                          Si(:, :,:,isVice))
        IF (Licetrc(ng)) THEN
          CALL exchange_r4d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            1, 2, 1, NTICE(ng),                   &
     &                            Ti)
        END IF
      END IF
# ifdef DISTRIBUTE
!
      CALL mp_exchange3d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj, 1, 2,                     &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Si(:,:,:,isAice),                             &
     &                    Si(:,:,:,isHice),                             &
     &                    Si(:,:,:,isHmel),                             &
     &                    Si(:,:,:,isHsno))
      CALL mp_exchange3d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj, 1, 2,                     &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Si(:,:,:,isIage),                             &
     &                    Si(:,:,:,isISxx),                             &
     &                    Si(:,:,:,isISxy),                             &
     &                    Si(:,:,:,isISyy))
      CALL mp_exchange3d (ng, tile, model, 3,                           &
     &                    LBi, UBi, LBj, UBj, 1, 2,                     &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Si(:,:,:,isTice),                             &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice))
      IF (Licetrc(ng)) THEN
        CALL mp_exchange4d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj, 1, 2, 1, NTICE(ng),     &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      Ti)
      END IF
# endif
!
      RETURN
      END SUBROUTINE ice_ini_fields_tile
#endif
      END MODULE ice_ini_fields_mod
