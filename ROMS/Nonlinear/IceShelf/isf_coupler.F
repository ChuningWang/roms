#include "cppdefs.h"

      MODULE isf_coupler_mod

#if defined ICESHEET_SSA && defined ICESHELF_THERMO
!
!git $Id$
!====================================================== Chuning Wang ===
!  Copyright (c) 2002-2025 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  Ice-ocean coupling routines for the dynamic ice sheet model.        !
!                                                                      !
!  This module provides:                                               !
!    1. Ice-ocean momentum coupling: computes the drag stress at the   !
!       ice-ocean interface from the relative velocity between ice     !
!       (SSA) and ocean (top layer), and writes it to sustr/svstr      !
!       so the ocean momentum equation sees the correct top BC.        !
!    2. Friction velocity update: computes u* = sqrt(Cd_isf)*|du| at       !
!       RHO points for use by the 3EQN thermodynamics (isf_vbc).       !
!    3. Basal melt rate temporal filtering: applies an exponential     !
!       moving average (EMA) to smooth high-frequency melt variability !
!       before it feeds back to ice thickness evolution.               !
!                                                                      !
!  References:                                                         !
!    Holland & Jenkins (1999), J. Phys. Oceanogr. 29, 1787-1800.      !
!    McPhee (1987), J. Geophys. Res. 92, 7017-7032.                    !
!    MITgcm SHELFICE: shelfice_v_drag_coeff.F                          !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC :: isf_coupler
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE isf_coupler (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
      USE mod_isf
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
      CALL isf_coupler_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nrhs(ng),                                  &
     &                       GRID(ng) % rmask_isf,                      &
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
     &                       ISF(ng) % Si,                              &
     &                       ISF(ng) % Fi,                              &
     &                       FORCES(ng) % sustr,                        &
     &                       FORCES(ng) % svstr)
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE isf_coupler
!
!***********************************************************************
      SUBROUTINE isf_coupler_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nrhs,                                &
     &                             rmask_isf,                           &
     &                             u, v,                                &
     &                             Si, Fi,                              &
     &                             sustr, svstr)
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
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: rmask_isf(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: Si(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: Fi(LBi:,LBj:,:)
      real(r8), intent(inout) :: sustr(LBi:,LBj:)
      real(r8), intent(inout) :: svstr(LBi:,LBj:)
# else
      real(r8), intent(in) :: rmask_isf(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: Si(LBi:UBi,LBj:UBj,2,nIsfS)
      real(r8), intent(inout) :: Fi(LBi:UBi,LBj:UBj,nIsfF)
      real(r8), intent(inout) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: svstr(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, know
!
      real(r8) :: alpha, dt_step, du, dv, spd
      real(r8) :: ui_u, vi_v, ustar_rho
      real(r8) :: cff, cff1, cff2, melt_inst, melt_prev
!
# include "set_bounds.h"
!
!=======================================================================
!  Step 1: Ice-ocean momentum coupling.
!
!  Compute the drag stress at the ice-ocean interface:
!    tau_io = Cd_isf * |u_ice - u_ocean| * (u_ice - u_ocean)
!
!  Ice velocity is at PSI points; interpolate to U/V points.
!  Ocean velocity is at U/V points (top layer k=N).
!  Write stress to sustr/svstr (ocean surface stress).
!  Compute friction velocity u* at RHO points for 3EQN.
!=======================================================================
!
      know=lsnew(ng)
      dt_step=dt(ng)
!
!  U-stress at U-points.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
!
!  Check if ice shelf cavity exists at both adjacent RHO points.
!  Same criterion as set_vbc.F (rmask_isf*rmask_isf.eq.1) so the
!  drag stress is only applied where both cells are cavity water.
!
          IF (rmask_isf(i-1,j)*rmask_isf(i,j).eq.1.0_r8) THEN
!
!  Interpolate ice velocity from PSI to U-point.
!  PSI(i,j) is at the corner between RHO(i-1,j-1), RHO(i,j-1),
!  RHO(i-1,j), RHO(i,j).  U(i,j) is between RHO(i-1,j) and RHO(i,j).
!  Average the two adjacent PSI values in y:
!    u_ice(U) = 0.5*(Si(i,j) + Si(i,j+1))
!
            ui_u=0.5_r8*(Si(i,j  ,know,isIuis)+                         &
     &                   Si(i,j+1,know,isIuis))
!
!  Relative velocity: ice - ocean.
!  Ocean v is at V-points; interpolate to U-point (i,j) with the
!  same 4-point average as set_vbc.F:
!    v(U) = 0.25*(v(i-1,j) + v(i,j) + v(i-1,j+1) + v(i,j+1))
!
            du=ui_u-u(i,j,N(ng),nrhs)
            dv=0.5_r8*(Si(i,j  ,know,isIvis)+                           &
     &                 Si(i,j+1,know,isIvis))-                          &
     &         0.25_r8*(v(i-1,j  ,N(ng),nrhs)+                          &
     &                  v(i  ,j  ,N(ng),nrhs)+                          &
     &                  v(i-1,j+1,N(ng),nrhs)+                          &
     &                  v(i  ,j+1,N(ng),nrhs))
            spd=SQRT(du*du+dv*dv)
!
!  Kinematic stress (divided by rho0): tau/rho0 = Cd_isf*|du|*du.
!
            sustr(i,j)=Cd_isf(ng)*spd*du
          ELSE
!  Outside ice shelf cavity: zero stress so migrating grounding lines
!  do not leave stale drag values on newly-exposed ocean cells.
            sustr(i,j)=0.0_r8
          END IF
        END DO
      END DO
!
!  V-stress at V-points.
!
      DO j=JstrV,Jend
        DO i=Istr,Iend
          IF (rmask_isf(i,j-1)*rmask_isf(i,j).eq.1.0_r8) THEN
!
!  Interpolate ice velocity from PSI to V-point.
!  V(i,j) is between RHO(i,j-1) and RHO(i,j).
!  Average the two adjacent PSI values in x:
!    v_ice(V) = 0.5*(Si(i,j) + Si(i+1,j))
!
            vi_v=0.5_r8*(Si(i  ,j,know,isIvis)+                         &
     &                   Si(i+1,j,know,isIvis))
!
!  Relative velocity: ice - ocean.
!  Ocean u is at U-points; interpolate to V-point (i,j) with the
!  same 4-point average as set_vbc.F:
!    u(V) = 0.25*(u(i,j-1) + u(i+1,j-1) + u(i,j) + u(i+1,j))
!
            dv=vi_v-v(i,j,N(ng),nrhs)
            du=0.5_r8*(Si(i  ,j,know,isIuis)+                           &
     &                 Si(i+1,j,know,isIuis))-                          &
     &         0.25_r8*(u(i  ,j-1,N(ng),nrhs)+                          &
     &                  u(i+1,j-1,N(ng),nrhs)+                          &
     &                  u(i  ,j  ,N(ng),nrhs)+                          &
     &                  u(i+1,j  ,N(ng),nrhs))
            spd=SQRT(du*du+dv*dv)
!
            svstr(i,j)=Cd_isf(ng)*spd*dv
          ELSE
            svstr(i,j)=0.0_r8
          END IF
        END DO
      END DO
!
!  Exchange halo of the ice-shelf drag stress so that Step 2 (friction
!  velocity at RHO points) reads consistent U/V-point values at tile
!  boundaries.  Mirrors the boundary-condition treatment at the end of
!  the ice shelf cavity stress block in set_vbc.F.
!
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  sustr)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  svstr)
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    sustr, svstr)
# endif
!
!=======================================================================
!  Step 2: Compute friction velocity u* at RHO points.
!
!  u* = sqrt(Cd_isf) * |u_ice - u_ocean| at RHO points.
!  Interpolate stress from U/V to RHO and compute:
!    u* = (tau_x^2 + tau_y^2)^{1/4}
!
!  This replaces the wind-stress-based u* in isf_vbc under ICESHEET_SSA.
!=======================================================================
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (rmask_isf(i,j).eq.1.0_r8) THEN
!
!  Average sustr from surrounding U-points and svstr from V-points.
!  u* = (tau_x^2 + tau_y^2)^{1/4} = sqrt(|tau|), where tau = Cd_isf*spd*du.
!
            cff1=0.5_r8*(sustr(i,j)+sustr(i+1,j))
            cff2=0.5_r8*(svstr(i,j)+svstr(i,j+1))
            ustar_rho=SQRT(SQRT(cff1*cff1+cff2*cff2))
            Fi(i,j,icIust)=MAX(1.0E-6_r8, ustar_rho)
          ELSE
            Fi(i,j,icIust)=0.0_r8
          END IF
        END DO
      END DO
!
!=======================================================================
!  Step 3: Basal melt rate EMA temporal filtering.
!
!  Apply exponential moving average to smooth high-frequency melt
!  variability before it feeds back to ice thickness evolution:
!    melt_filtered = alpha*melt_instant + (1-alpha)*melt_prev
!  where alpha = dt / (tau_melt + dt).
!
!  Fi(i,j,icImel) holds the instantaneous melt rate from isf_vbc.
!  Fi(i,j,icImef) stores the EMA state (persistent between steps).
!  isf_thickness_adv reads Fi(i,j,icImef) for the thickness source.
!=======================================================================
!
      IF (tau_melt(ng).gt.0.0_r8) THEN
        alpha=dt_step/(tau_melt(ng)+dt_step)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            IF (rmask_isf(i,j).eq.1.0_r8) THEN
              melt_inst=Fi(i,j,icImel)
              melt_prev=Fi(i,j,icImef)
              Fi(i,j,icImef)=alpha*melt_inst+                           &
     &                       (1.0_r8-alpha)*melt_prev
            ELSE
              Fi(i,j,icImef)=0.0_r8
            END IF
          END DO
        END DO
      ELSE
!  No filtering: pass through instantaneous melt rate.
        DO j=Jstr,Jend
          DO i=Istr,Iend
            Fi(i,j,icImef)=Fi(i,j,icImel)
          END DO
        END DO
      END IF
!
      RETURN
      END SUBROUTINE isf_coupler_tile
!
#endif
      END MODULE isf_coupler_mod
