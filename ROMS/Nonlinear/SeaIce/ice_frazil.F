#include "cppdefs.h"
      MODULE ice_frazil_mod

#if defined ICE_MODEL && defined ICE_THERMO
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2024 The ROMS/TOMS Group      Paul Budgell       !
!    Licensed under a MIT/X style license           Katherine Hedstrom !
!    See License_ROMS.md                            Scott M. Durski    !
!=======================================================================
!                                                                      !
!  This routine computes the frazil ice growth in the water when the   !
!  water temperature gets below freezing. It adjusts both the water    !
!  temperature and salinity accordingly.                               !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Steele, M., G.L. Mellor, and M.G. Mcphee, 1989: Role of the       !
!      molecular sublayer in the melting or freezing of sea ice,       !
!      J. Phys. Oceanogr., 19, 139-147.                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ice
      USE mod_ocean
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
# endif
!
      implicit none
!
      PRIVATE
      PUBLIC  :: ice_frazil
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_frazil (ng, tile, model)
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
      CALL ice_frazil_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), linew(ng),                        &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
# endif
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % z_r,                             &
     &                      OCEAN(ng) % rho,                            &
     &                      OCEAN(ng) % t,                              &
     &                      ICE(ng) % Fi, ICE(ng) % Ti)
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE ice_frazil
!
!***********************************************************************
      SUBROUTINE ice_frazil_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nrhs, linew,                          &
# ifdef MASKING
     &                            rmask,                                &
# endif
# ifdef WET_DRY
     &                            rmask_wet,                            &
# endif
     &                            Hz, z_r, rho, t,                      &
     &                            Fi, Ti)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, linew

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: Fi(LBi:,LBj:,:)
      real(r8), intent(inout) :: Ti(LBi:,LBj:,:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: Fi(LBi:UBi,LBj:UBj,nIceF)
      real(r8), intent(inout) :: Ti(LBi:UBi,LBj:UBj,2,NTICE(ng))
# endif
!
! Local variable definitions
!
      logical :: overLand
!
      integer :: i, j, k, ic, itrc
!
      real(r8), parameter :: eps      = 1.0E-4_r8       ! zero division
      real(r8), parameter :: r        = 0.52_r8         ! Cp_i/Cp_w
      real(r8), parameter :: Lhat     = 83.88_r8        ! hfus/Cp_w [K]
      real(r8), parameter :: frln     = -0.0543_r8      ! [psu C-1]
      real(r8), parameter :: sice_ref = 3.2_r8          ! [psu]
!
      real(r8), dimension(NAT+1:NT(ng)) :: trc
!
      real(r8) :: delta_wfr, gamma_k, orhoi, t_fr, sice
      real(r8) :: cff, cff1, cff2

# ifdef DISTRIBUTE
      real(r8), allocatable           :: buffer(:)
      character (len=3), allocatable  :: op_handle(:)
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute frazil ice growth.
!-----------------------------------------------------------------------
!
!  Initialize rate of frazil ice growth (m3/s).
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Fi(i,j,icW_fr)=0.0_r8
        END DO
      END DO
!
!  Original formulation.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
# ifdef MASKING
          overLand=rmask(i,j).lt.1.0_r8
#  ifdef WET_DRY
          overLand=overland.or.(rmask_wet(i,j).lt.1.0_r8)
#  endif
# else
          overLand=.FALSE.
# endif
          IF (.not.overLand) THEN
            IF (Licetrc(ng)) THEN
              DO itrc=NAT+1,NT(ng)
                IF (LtracerIce(itrc,ng)) THEN
                  trc(itrc)=0.0_r8
                END IF
              END DO
            END IF
            DO k=1,N(ng)
              orhoi = 1.0_r8/IceRho(ng)
              t_fr = frln*t(i,j,k,nrhs,isalt) + 0.000759_r8*z_r(i,j,k)
              delta_wfr = 0._r8
              IF (t(i,j,k,nrhs,itemp).lt.t_fr) THEN
                gamma_k = (t_fr - t(i,j,k,nrhs,itemp)) /                &
     &                    (Lhat + t(i,j,k,nrhs,itemp)*(1.0_r8-r) - t_fr)
                delta_wfr = gamma_k*Hz(i,j,k)*(rhow+rho(i,j,k))*orhoi
                Fi(i,j,icW_fr) = Fi(i,j,icW_fr) + delta_wfr
                t(i,j,k,nrhs,itemp) = t(i,j,k,nrhs,itemp) +             &
     &            gamma_k*(Lhat + t(i,j,k,nrhs,itemp)*(1.0_r8-r))
                sice = MIN(sice_ref, t(i,j,k,nrhs,isalt))
                t(i,j,k,nrhs,isalt) = t(i,j,k,nrhs,isalt) *             &
     &            (1.0_r8+gamma_k)
              ELSE IF ((Fi(i,j,icW_fr).gt.0.0_r8).and.                  &
     &                 (t(i,j,k,nrhs,itemp).gt.t_fr)) THEN
                gamma_k = (t_fr - t(i,j,k,nrhs,itemp)) /                &
     &                    (Lhat + t(i,j,k,nrhs,itemp)*(1.0_r8-r) - t_fr)
                delta_wfr = gamma_k*Hz(i,j,k)*(rhow+rho(i,j,k))*orhoi
                IF ((Fi(i,j,icW_fr)+delta_wfr).gt.0.0_r8) THEN
                  Fi(i,j,icW_fr) = Fi(i,j,icW_fr) + delta_wfr
                ELSE
                  gamma_k = -Fi(i,j,icW_fr)*                            &
     &                      IceRho(ng)/((rhow+rho(i,j,k))*Hz(i,j,k))
                  delta_wfr = -Fi(i,j,icW_fr)
                  Fi(i,j,icW_fr) = 0.0_r8
                ENDIF
                t(i,j,k,nrhs,itemp) = t(i,j,k,nrhs,itemp) +             &
     &            gamma_k*(Lhat + t(i,j,k,nrhs,itemp)*(1.0_r8-r))
                sice = MIN(sice_ref, t(i,j,k,nrhs,isalt))
                t(i,j,k,nrhs,isalt) = t(i,j,k,nrhs,isalt) *             &
     &            (1.0_r8+gamma_k)
              END IF
              IF (Licetrc(ng)) THEN
                cff = delta_wfr*IceRho(ng)/(rhow+rho(i,j,k))
                DO itrc=NAT+1,NT(ng)
                  IF (LtracerIce(itrc,ng)) THEN
                    IF (delta_wfr.lt.0.0_r8) THEN
                      cff1=trc(itrc)/                                   &
     &  (Fi(i,j,icW_fr)*(IceRho(ng)/(rhow+rho(i,j,k)))+eps)
                    ELSE
                      cff1=t(i,j,k,nrhs,itrc)
                    ENDIF
                    t(i,j,k,nrhs,itrc) = t(i,j,k,nrhs,itrc) -           &
     &                cff/Hz(i,j,k)*(cff1-t(i,j,k,nrhs,itrc))
                    trc(itrc) = trc(itrc) + cff*cff1
                    trc(itrc) = MAX(trc(itrc), 0.0_r8)
                  END IF
                END DO
              END IF
              IF (Lfrzsrc(ng)) THEN
                DO itrc=NAT+1,NT(ng)
                  IF (LsourceFrz(itrc,ng)) THEN
                    t(i,j,k,nrhs,itrc) = t(i,j,k,nrhs,itrc) +           &
     &                MAX(delta_wfr, 0.0_r8)/Hz(i,j,k)*trcFrz(itrc,ng)
                  END IF
                END DO
              END IF
            END DO
!
!  Compute rate of frazil ice growth.
!
            Fi(i,j,icW_fr) = Fi(i,j,icW_fr) / dt(ng)
!
!  Update tracer concentration.
!
            IF (Licetrc(ng)) THEN
              ic=0
              DO itrc=NAT+1,NT(ng)
                IF (LtracerIce(itrc,ng)) THEN
                  ic=ic+1
                  Ti(i,j,linew,ic) = Ti(i,j,linew,ic)+trc(itrc)
                END IF
              END DO
            END IF
          END IF
        END DO
      END DO
!
!  Apply periodic boundary conditions.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          t(:,:,:,nrhs,itemp))
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          t(:,:,:,nrhs,isalt))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fi(:,:,icW_fr))
        IF (Licetrc(ng)) THEN
          DO itrc=NAT+1,NT(ng)
            IF (LtracerIce(itrc,ng)) THEN
              CALL exchange_r3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                t(:,:,:,nrhs,itrc))
            END IF
          END DO
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, NTICE(ng),     &
     &                            Ti(:,:,linew,:))
        END IF
      END IF
!
# ifdef DISTRIBUTE
!
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    t(:,:,:,nrhs,itemp))
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    t(:,:,:,nrhs,isalt))
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Fi(:,:,icW_fr))

      IF (Licetrc(ng)) THEN
        DO itrc=NAT+1,NT(ng)
          IF (LtracerIce(itrc,ng)) THEN
            CALL mp_exchange3d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          t(:,:,:,nrhs,itrc))
          END IF
        END DO
        CALL mp_exchange3d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj, 1, NTICE(ng),           &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      Ti(:,:,linew,:))
      END IF
# endif
!
      RETURN
      END SUBROUTINE ice_frazil_tile
#endif
      END MODULE ice_frazil_mod

