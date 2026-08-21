#include "cppdefs.h"

      MODULE isf_ini_fields_mod

#ifdef ICESHEET
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2024 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           W. Paul Budgell    !
!    See License_ROMS.md                            Katherine Hedstrom !
!================================================== Hernan G. Arango ===
!                                                       Chuning Wang   !
!  Sets the initial fields for the iceshelf state variables and        !
!  applies their lateral boundary conditions.                          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_isf
      USE mod_ncparam
      USE mod_scalars
      USE mod_stepping
!
      USE isf_bc2d_mod,     ONLY : isf_bc2d_tile
      USE isf_pibc_mod,     ONLY : isf_pibc_tile
      USE isf_tibc_mod,     ONLY : isf_tibc_tile
!
      USE exchange_2d_mod,  ONLY : exchange_r2d_tile, exchange_p2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod,  ONLY : mp_exchange2d
# endif
!
      implicit none
!
      PRIVATE
      PUBLIC  :: isf_ini_fields
      PUBLIC  :: isf_ini_fields_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE isf_ini_fields (ng, tile, model)
!***********************************************************************
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
      CALL isf_ini_fields_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          lsold(ng),                              &
# ifdef MASKING
     &                          GRID(ng) % rmask,                       &
     &                          GRID(ng) % pumask_ist,                  &
     &                          GRID(ng) % pvmask_ist,                  &
# endif
     &                          ISF(ng) % Si,                           &
     &                          ISF(ng) % Fi)
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE isf_ini_fields
!
!***********************************************************************
      SUBROUTINE isf_ini_fields_tile (ng, tile, model,                  &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                lsold,                            &
# ifdef MASKING
     &                                rmask, pumask_ist, pvmask_ist,    &
# endif
     &                                Si, Fi)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: lsold
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: pumask_ist(LBi:,LBj:)
      real(r8), intent(in) :: pvmask_ist(LBi:,LBj:)
# endif
      real(r8), intent(inout) :: Si(LBi:UBi,LBj:UBj,2,nIsfS)
      real(r8), intent(inout) :: Fi(LBi:UBi,LBj:UBj,nIsfF)
!
!  Local variable declarations.
!
      integer :: i, j, lsnew
!
# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  If not a perfect restart, initialize the "lsnew" time level of the
!  iceshelf state variables to the "lsold" values.
!-----------------------------------------------------------------------
!
      lsnew=3-lsold
!
      IF (.not.PerfectRST(ng)) THEN
        DO j=JstrB,JendB
          DO i=IstrB,IendB
            Si(i,j,lsnew,isIhis)=Si(i,j,lsold,isIhis)
            Si(i,j,lsnew,isIeis)=Si(i,j,lsold,isIeis)
# ifdef MASKING
            Si(i,j,lsold,isIhis)=Si(i,j,lsold,isIhis)*rmask(i,j)
            Si(i,j,lsnew,isIhis)=Si(i,j,lsnew,isIhis)*rmask(i,j)
            Si(i,j,lsold,isIeis)=Si(i,j,lsold,isIeis)*rmask(i,j)
            Si(i,j,lsnew,isIeis)=Si(i,j,lsnew,isIeis)*rmask(i,j)
# endif
          END DO
        END DO
!
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            Si(i,j,lsnew,isIuis)=Si(i,j,lsold,isIuis)
            Si(i,j,lsnew,isIvis)=Si(i,j,lsold,isIvis)
# ifdef MASKING
            Si(i,j,lsold,isIuis)=Si(i,j,lsold,isIuis)*pumask_ist(i,j)
            Si(i,j,lsnew,isIuis)=Si(i,j,lsnew,isIuis)*pumask_ist(i,j)
            Si(i,j,lsold,isIvis)=Si(i,j,lsold,isIvis)*pvmask_ist(i,j)
            Si(i,j,lsnew,isIvis)=Si(i,j,lsnew,isIvis)*pvmask_ist(i,j)
# endif
          END DO
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions to the iceshelf state variables
!  for both time levels.
!-----------------------------------------------------------------------
!
      CALL isf_pibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    Si(:,:,lsold,isIuis),                         &
     &                    Si(:,:,lsold,isIvis))
!
      CALL isf_pibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    Si(:,:,lsnew,isIuis),                         &
     &                    Si(:,:,lsnew,isIvis))
!
!  Thickness must be set before enthalpy: isf_tibc_tile computes the
!  enthalpy boundary value as h*T and divides by the boundary thickness.
!
      CALL isf_bc2d_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    Si(:,:,lsold,isIhis))
!
      CALL isf_bc2d_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    Si(:,:,lsnew,isIhis))
!
      CALL isf_tibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Si(:,:,lsold,isIhis),                         &
     &                    Si(:,:,lsold,isIeis),                         &
     &                    Fi(:,:,icItis))
!
      CALL isf_tibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Si(:,:,lsnew,isIhis),                         &
     &                    Si(:,:,lsnew,isIeis),                         &
     &                    Fi(:,:,icItis))
!
!-----------------------------------------------------------------------
!  Exchange halo data in the periodic or distributed-memory case.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile, LBi, UBi, LBj, UBj,           &
     &                          Si(:,:,lsold,isIhis))
        CALL exchange_r2d_tile (ng, tile, LBi, UBi, LBj, UBj,           &
     &                          Si(:,:,lsnew,isIhis))
        CALL exchange_r2d_tile (ng, tile, LBi, UBi, LBj, UBj,           &
     &                          Si(:,:,lsold,isIeis))
        CALL exchange_r2d_tile (ng, tile, LBi, UBi, LBj, UBj,           &
     &                          Si(:,:,lsnew,isIeis))
        CALL exchange_p2d_tile (ng, tile, LBi, UBi, LBj, UBj,           &
     &                          Si(:,:,lsold,isIuis))
        CALL exchange_p2d_tile (ng, tile, LBi, UBi, LBj, UBj,           &
     &                          Si(:,:,lsnew,isIuis))
        CALL exchange_p2d_tile (ng, tile, LBi, UBi, LBj, UBj,           &
     &                          Si(:,:,lsold,isIvis))
        CALL exchange_p2d_tile (ng, tile, LBi, UBi, LBj, UBj,           &
     &                          Si(:,:,lsnew,isIvis))
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 4, LBi, UBi, LBj, UBj,        &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Si(:,:,lsold,isIhis),                         &
     &                    Si(:,:,lsnew,isIhis),                         &
     &                    Si(:,:,lsold,isIeis),                         &
     &                    Si(:,:,lsnew,isIeis))
      CALL mp_exchange2d (ng, tile, iNLM, 4, LBi, UBi, LBj, UBj,        &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Si(:,:,lsold,isIuis),                         &
     &                    Si(:,:,lsnew,isIuis),                         &
     &                    Si(:,:,lsold,isIvis),                         &
     &                    Si(:,:,lsnew,isIvis))
# endif
!
      RETURN
      END SUBROUTINE isf_ini_fields_tile
#endif
      END MODULE isf_ini_fields_mod
