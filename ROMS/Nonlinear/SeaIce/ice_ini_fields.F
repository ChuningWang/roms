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
      USE mod_ice
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
     &                          linew(ng), liold(ng),                   &
     &                          nIceS, nIceF, NTICE(ng), Licetrc(ng),   &
     &                          isUice, isVice,                         &
     &                          isAice, isHice, isHsno, isHmel, isIage, &
# ifdef ICE_BIO
     &                          isIphy, isINO3, isINH4,                 &
# endif
     &                          ibICE,                                  &
# ifdef MASKING
     &                          GRID(ng) % rmask,                       &
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
     &                                linew, liold,                     &
     &                                nIceS, nIceF, NTICE, Licetrc,     &
     &                                isUice, isVice,                   &
     &                                isAice, isHice, isHsno, isHmel,   &
     &                                isIage,                           &
# ifdef ICE_BIO
     &                                isIphy, isINO3, isINH4,           &
# endif
     &                                ibICE,                            &
# ifdef MASKING
     &                                rmask,                            &
# endif
     &                                Si, Fi, Ti)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      USE bc_4d_mod,       ONLY : bc_r4d_tile
      USE ice_bc2d_mod,    ONLY : ice_bc2d_tile
      USE ice_tibc_mod,    ONLY : ice_tibc_tile
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
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
      integer, intent(in) :: linew, liold
      integer, intent(in) :: nIceS, nIceF, NTICE
      logical, intent(in) :: Licetrc
      integer, intent(in) :: isUice, isVice
      integer, intent(in) :: isAice, isHice, isHsno, isHmel, isIage
# ifdef ICE_BIO
      integer, intent(in) :: isIphy, isINO3, isINH4
# endif
      integer, intent(in) :: ibICE(nIceS)
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: Fi(LBi:,LBj:,:)
      real(r8), intent(inout) :: Si(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: Ti(LBi:,LBj:,:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: Fi(LBi:UBi,LBj:UBj,nIceF)
      real(r8), intent(inout) :: Si(LBi:UBi,LBj:UBj,2,nIceS)
      real(r8), intent(inout) :: Ti(LBi:UBi,LBj:UBj,2,NTICE)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, nf, ic, itrc

      real(r8) :: cff

# include "set_bounds.h"

!
!-----------------------------------------------------------------------
!  If not perfect restart, initialize other time levels for ice state
!  variables.
!-----------------------------------------------------------------------
!
      IF (.not.PerfectRST(ng)) THEN
        DO nf=1,nIceS
          DO j=JstrB,JendB
            DO i=IstrB,IendB
              cff=Si(i,j,liold,nf)
# ifdef MASKING
              cff=cff1*rmask(i,j)
# endif
              Si(i,j,liold,nf)=cff
              Si(i,j,linew,nf)=cff
            END DO
          END DO
        END DO
!
        DO nf=1,NTICE
          DO j=JstrB,JendB
            DO i=IstrB,IendB
              cff=Ti(i,j,liold,nf)
# ifdef MASKING
              cff=cff1*rmask(i,j)
# endif
              Ti(i,j,liold,nf)=cff
              Ti(i,j,linew,nf)=cff
            END DO
          END DO
        END DO
!
!  Apply boundary conditions.
!
        CALL ice_bc2d_tile (ng, tile, model, isAice,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, linew,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isAice),                           &
     &                      LBC(:,ibICE(isAice),ng))

        CALL ice_bc2d_tile (ng, tile, model, isHice,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, linew,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isHice),                           &
     &                      LBC(:,ibICE(isHice),ng))

        CALL ice_bc2d_tile (ng, tile, model, isHsno,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, linew,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isHsno),                           &
     &                      LBC(:,ibICE(isHsno),ng))

        CALL ice_bc2d_tile (ng, tile, model, isHmel,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, linew,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isHmel),                           &
     &                      LBC(:,ibICE(isHmel),ng))

        CALL ice_bc2d_tile (ng, tile, model, isIage,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, linew,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isIage),                           &
     &                      LBC(:,ibICE(isIage),ng))

# ifdef ICE_BIO
        CALL ice_bc2d_tile (ng, tile, model, isIphy,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, linew,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isIphy),                           &
     &                      LBC(:,ibICE(isIphy),ng))

        CALL ice_bc2d_tile (ng, tile, model, isINO3,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, linew,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isINO3),                           &
     &                      LBC(:,ibICE(isINO3),ng))

        CALL ice_bc2d_tile (ng, tile, model, isINH4,                    &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, linew,                               &
     &                      Si(:,:,:,isUice),                           &
     &                      Si(:,:,:,isVice),                           &
     &                      Si(:,:,:,isINH4),                           &
     &                      LBC(:,ibICE(isINH4),ng))
# endif

        IF (Licetrc) THEN
          CALL bc_r4d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, 1, 2, 1, NTICE,         &
     &                      Ti)
        END IF
!
      END IF
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r4d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, 2, 1, nIceS,     &
     &                          Si)
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, nIceF,           &
     &                          Fi)
        IF (Licetrc) THEN
          CALL exchange_r4d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, 2, 1, NTICE,   &
     &                            Ti)
        END IF
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, 2, 1, nIceS,           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Si)
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, nIceF,                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Fi)
      IF (Licetrc) THEN
        CALL mp_exchange4d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj, 1, 2, 1, NTICE,         &
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
