#include "cppdefs.h"

      MODULE isf_advection_mod

#ifdef ICESHEET
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2024 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module advects the icesheet scalar state: ice thickness        !
!  (mass conservation) and ice enthalpy (heat content),                !
!                                                                      !
!    dH/dt = -div(u*H) + SMB - basal_melt                              !
!    dE/dt = -div(u*E) + SMB*T - basal_melt*T                          !
!                                                                      !
!  where u is the depth-averaged ice velocity from the SSA solver      !
!  (PSI points) interpolated to RHO cell faces, H is ice thickness,    !
!  E = H*T is the ice enthalpy, and T is the ice temperature           !
!  (diagnostic, Fi(:,:,icItis)).  The melt and SMB enthalpy source     !
!  terms act at the column temperature T = E/H so that mass and heat   !
!  sources do not alter the ice temperature.                           !
!                                                                      !
!  Numerical scheme (first-order upwind, flux form):                   !
!    - PSI velocity interpolated to RHO cell faces (2-point average)   !
!    - First-order upwind flux for monotonicity                        !
!    - CFL-based sub-stepping for stability                            !
!    - Positivity preservation (H >= 0, E = 0 where H = 0)             !
!                                                                      !
!  Reference: STREAMICE streamice_advect_thickness.F                   !
!             PISM mass/transport.cc                                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
      PRIVATE
      PUBLIC  :: isf_advection
!
!  Global dot product shared variable (for MPI reduction).
!
      real(r8) :: isf_dot_prod
!$OMP THREADPRIVATE (isf_dot_prod)
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE isf_advection (ng, tile, model)
!***********************************************************************
!
      USE mod_grid
      USE mod_stepping, ONLY : lsold, lsnew
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
      CALL wclock_on (ng, model, 43, __LINE__, MyFile)
# endif
      CALL isf_advection_tile (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         GRID(ng) % pm,                           &
     &                         GRID(ng) % pn,                           &
     &                         ISF(ng) % Si,                            &
     &                         ISF(ng) % Fi)
# ifdef PROFILE
      CALL wclock_off (ng, model, 43, __LINE__, MyFile)
# endif
      RETURN
      END SUBROUTINE isf_advection
!
!***********************************************************************
      SUBROUTINE isf_advection_tile (ng, tile, model,                   &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               pm, pn,                            &
     &                               Si, Fi)
!***********************************************************************
!
!  Solves the ice thickness (mass) and enthalpy (heat content)
!  evolution equations in flux form:
!
!    dH/dt = -div(u*H) + SMB - melt_rate
!    dE/dt = -div(u*E) + SMB*T - melt_rate*T
!
!  where u is the depth-averaged ice velocity (from SSA solver),
!  H is ice thickness, E = H*T is the ice enthalpy, T is the ice
!  temperature, and melt_rate is the basal melt rate from the
!  3-equation thermodynamics (positive = melting).
!
!  Numerical scheme:
!    - PSI velocity interpolated to RHO cell faces (2-point average)
!    - First-order upwind flux for monotonicity
!    - CFL-based sub-stepping for stability
!    - Positivity preservation (H >= 0, E = 0 where H = 0)
!
!  The enthalpy fluxes take the upwind state from the same cell as
!  the thickness fluxes so that the advected temperature T = E/H is
!  preserved under pure advection.
!
!  Grid layout (B-grid / PSI velocity):
!    PSI(i,j) is the SW corner of RHO cell (i,j).
!    x-face between RHO(i-1,j) and RHO(i,j) has velocity:
!      u_face(i,j) = 0.5*(u_psi(i,j) + u_psi(i,j+1))
!    y-face between RHO(i,j-1) and RHO(i,j) has velocity:
!      v_face(i,j) = 0.5*(v_psi(i,j) + v_psi(i+1,j))
!
!  Reference: STREAMICE streamice_advect_thickness.F,
!             PISM mass transport (upwind FD)
!
      USE mod_scalars
      USE mod_parallel
      USE mod_stepping
      USE mod_iounits
      USE mod_isf
!
# ifdef DISTRIBUTE
      USE distribute_mod,  ONLY : mp_reduce
# endif
      USE mp_exchange_mod, ONLY : mp_exchange2d
      USE exchange_2d_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)    :: pm(LBi:,LBj:)
      real(r8), intent(in)    :: pn(LBi:,LBj:)
      real(r8), intent(inout) :: Si(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: Fi(LBi:,LBj:,:)
# else
      real(r8), intent(in)    :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Si(LBi:UBi,LBj:UBj,2,nIsfS)
      real(r8), intent(inout) :: Fi(LBi:UBi,LBj:UBj,nIsfF)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
      integer :: ksa_old, ksa_new
!
      real(r8) :: dt_full, dt_rem, dt_sub
      real(r8) :: max_icfl, loc_icfl, cfl_u, cfl_v
      real(r8) :: H_old, E_old, T_old
      real(r8) :: H_left, H_right, H_below, H_above
      real(r8) :: E_left, E_right, E_below, E_above
      real(r8) :: Fx_left, Fx_right, Fy_bottom, Fy_top
      real(r8) :: Gx_left, Gx_right, Gy_bottom, Gy_top
!
!  Work arrays.
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: H_new
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: E_new
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: u_xface
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: v_yface
!
# ifdef DISTRIBUTE
      character(len=3) :: op_handle
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", isf_advection_tile"
!
# include "set_bounds.h"
!
      SourceFile=MyFile
!
!  Set time level indices.  The lsold/lsnew swap was performed in
!  iceshelf.F at the SSA solve cadence; the velocity used here is
!  Si(:,:,ksa_new,isIuis/isIvis) from the most recent SSA solve.
!
      ksa_old=lsold(ng)
      ksa_new=lsnew(ng)
!
!-----------------------------------------------------------------------
!  Advect ice thickness and enthalpy with CFL sub-stepping.
!-----------------------------------------------------------------------
!
!  Initialize work arrays.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          H_new(i,j)=0.0_r8
          E_new(i,j)=0.0_r8
        END DO
      END DO
!
!  Compute face velocities from PSI-point velocity.
!  u_xface(i,j): x-direction velocity on the x-face between
!                RHO(i-1,j) and RHO(i,j).
!  v_yface(i,j): y-direction velocity on the y-face between
!                RHO(i,j-1) and RHO(i,j).
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          u_xface(i,j)=0.5_r8*                                          &
     &        (Si(i,j,ksa_new,isIuis)+Si(i,j+1,ksa_new,isIuis))
          v_yface(i,j)=0.5_r8*                                          &
     &        (Si(i,j,ksa_new,isIvis)+Si(i+1,j,ksa_new,isIvis))
        END DO
      END DO
!
!  Compute maximum CFL number for sub-stepping.
!  icfl = max(|u_face|/dx + |v_face|/dy) over all active cells.
!
      max_icfl=1.0E-20_r8
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (Si(i,j,ksa_old,isIhis).gt.0.0_r8) THEN
            cfl_u=MAX(ABS(u_xface(i,j)),                              &
     &                ABS(u_xface(i+1,j)))*pm(i,j)
            cfl_v=MAX(ABS(v_yface(i,j)),                              &
     &                ABS(v_yface(i,j+1)))*pn(i,j)
            loc_icfl=cfl_u+cfl_v
            IF (loc_icfl.gt.max_icfl) max_icfl=loc_icfl
          END IF
        END DO
      END DO
!
!  Global maximum for distributed memory.
!
# ifdef DISTRIBUTE
      isf_dot_prod=max_icfl
      op_handle='MAX'
      CALL mp_reduce(ng,iNLM,1,isf_dot_prod,op_handle)
      max_icfl=isf_dot_prod
# endif
!
!  CFL sub-stepping loop.
!
      dt_full=dt(ng)*REAL(N_ssa(ng)*N_thick(ng),r8)
      dt_rem=dt_full
!
      DO WHILE (dt_rem.gt.1.0E-10_r8)
!
        dt_sub=MIN(cfl_factor(ng)/max_icfl, dt_rem)
!
!  Compute upwind fluxes and update thickness and enthalpy.
!  Thickness flux:  Fx = u_face * H_upwind
!  Enthalpy flux:   Gx = u_face * E_upwind
!
!  H_new = H_old - dt_sub/area * (Fx_out - Fx_in + Fy_out - Fy_in)
!                 - dt_sub * melt_rate + dt_sub * SMB
!  E_new = E_old - dt_sub/area * (Gx_out - Gx_in + Gy_out - Gy_in)
!                 - dt_sub * melt_rate*T + dt_sub * SMB*T
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            H_old=Si(i,j,ksa_old,isIhis)
            E_old=Si(i,j,ksa_old,isIeis)
!
            IF (H_old.gt.0.0_r8) THEN
!
!  Column ice temperature driving the melt/SMB enthalpy sources.
!  E = 0 wherever H = 0 (maintained by initialization, boundary
!  conditions and the clipping below), so the enthalpy upwind
!  states need no separate guard.
!
              T_old=E_old/MAX(H_old,eps)
!
!  X-direction fluxes.
!  Left x-face (i): between RHO(i-1,j) and RHO(i,j).
!  Right x-face (i+1): between RHO(i,j) and RHO(i+1,j).
!
              H_left=Si(i-1,j,ksa_old,isIhis)
              H_right=Si(i+1,j,ksa_old,isIhis)
              E_left=Si(i-1,j,ksa_old,isIeis)
              E_right=Si(i+1,j,ksa_old,isIeis)
!
!  Upwind: if u > 0, the state comes from left (i-1); if u < 0, from
!  right (i).
!
              IF (u_xface(i,j).ge.0.0_r8) THEN
                Fx_left=u_xface(i,j)*MAX(H_left,0.0_r8)
                Gx_left=u_xface(i,j)*E_left
              ELSE
                Fx_left=u_xface(i,j)*H_old
                Gx_left=u_xface(i,j)*E_old
              END IF
!
              IF (u_xface(i+1,j).ge.0.0_r8) THEN
                Fx_right=u_xface(i+1,j)*H_old
                Gx_right=u_xface(i+1,j)*E_old
              ELSE
                Fx_right=u_xface(i+1,j)*MAX(H_right,0.0_r8)
                Gx_right=u_xface(i+1,j)*E_right
              END IF
!
!  Y-direction fluxes.
!  Bottom y-face (j): between RHO(i,j-1) and RHO(i,j).
!  Top y-face (j+1): between RHO(i,j) and RHO(i,j+1).
!
              H_below=Si(i,j-1,ksa_old,isIhis)
              H_above=Si(i,j+1,ksa_old,isIhis)
              E_below=Si(i,j-1,ksa_old,isIeis)
              E_above=Si(i,j+1,ksa_old,isIeis)
!
              IF (v_yface(i,j).ge.0.0_r8) THEN
                Fy_bottom=v_yface(i,j)*MAX(H_below,0.0_r8)
                Gy_bottom=v_yface(i,j)*E_below
              ELSE
                Fy_bottom=v_yface(i,j)*H_old
                Gy_bottom=v_yface(i,j)*E_old
              END IF
!
              IF (v_yface(i,j+1).ge.0.0_r8) THEN
                Fy_top=v_yface(i,j+1)*H_old
                Gy_top=v_yface(i,j+1)*E_old
              ELSE
                Fy_top=v_yface(i,j+1)*MAX(H_above,0.0_r8)
                Gy_top=v_yface(i,j+1)*E_above
              END IF
!
!  Thickness and enthalpy update.  On a uniform grid:
!  div(uH) ~ (Fx_right - Fx_left)*pm + (Fy_top - Fy_bottom)*pn.
!
!  Melt rate source: instantaneous basal melt (icImel) from isf_vbc
!  [m/s, positive = melting].  No EMA-filtered field is maintained.
!
!  Surface mass balance (SMB): uniform rate smb_rate(ng) [m/s ice
!  equivalent, positive = accumulation].  Applied everywhere ice
!  exists (H > 0), both grounded and floating.
!
              H_new(i,j)=H_old                                         &
     &          -dt_sub*((Fx_right-Fx_left)*pm(i,j)                   &
     &                  +(Fy_top-Fy_bottom)*pn(i,j))                  &
# ifdef ICESHELF_THERMO
     &          -dt_sub*Fi(i,j,icImel)                                &
# endif
     &          +dt_sub*smb_rate(ng)
!
              E_new(i,j)=E_old                                         &
     &          -dt_sub*((Gx_right-Gx_left)*pm(i,j)                   &
     &                  +(Gy_top-Gy_bottom)*pn(i,j))                  &
# ifdef ICESHELF_THERMO
     &          -dt_sub*Fi(i,j,icImel)*T_old                          &
# endif
     &          +dt_sub*smb_rate(ng)*T_old
!
!  Positivity preservation and safety cap.  Where H is clipped, E is
!  zeroed or rescaled so the invariant E = H*T is retained.
!
              IF (H_new(i,j).lt.0.0_r8) THEN
                H_new(i,j)=0.0_r8
                E_new(i,j)=0.0_r8
              ELSE IF (H_new(i,j).gt.max_thickness(ng)) THEN
                E_new(i,j)=E_new(i,j)*max_thickness(ng)/H_new(i,j)
                H_new(i,j)=max_thickness(ng)
              END IF
!
            ELSE
              H_new(i,j)=0.0_r8
              E_new(i,j)=0.0_r8
            END IF
          END DO
        END DO
!
!  Update H and E for next sub-step (copy back to ksa_old for
!  multi-substep advection).
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            Si(i,j,ksa_old,isIhis)=H_new(i,j)
            Si(i,j,ksa_old,isIeis)=E_new(i,j)
          END DO
        END DO
!
!  Exchange halo for multi-substep advection (only if more sub-steps
!  remain).
!
        IF (dt_rem-dt_sub.gt.1.0E-10_r8) THEN
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile(ng,tile,LBi,UBi,LBj,UBj,          &
     &                             Si(:,:,ksa_old,isIhis))
            CALL exchange_r2d_tile(ng,tile,LBi,UBi,LBj,UBj,          &
     &                             Si(:,:,ksa_old,isIeis))
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange2d(ng,tile,iNLM,2,LBi,UBi,LBj,UBj,         &
     &                       NghostPoints,EWperiodic(ng),              &
     &                       NSperiodic(ng),                           &
     &                       Si(:,:,ksa_old,isIhis),                   &
     &                       Si(:,:,ksa_old,isIeis))
# endif
        END IF
!
        dt_rem=dt_rem-dt_sub
!
      END DO  ! sub-stepping
!
!  Store final thickness and enthalpy at ksa_new and update the ice
!  temperature diagnostic (T = E/H).
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          Si(i,j,ksa_new,isIhis)=H_new(i,j)
          Si(i,j,ksa_new,isIeis)=E_new(i,j)
          IF (H_new(i,j).gt.0.0_r8) THEN
            Fi(i,j,icItis)=E_new(i,j)/H_new(i,j)
          ELSE
            Fi(i,j,icItis)=0.0_r8
          END IF
        END DO
      END DO
!
!  Exchange halo for updated thickness and enthalpy.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile(ng,tile,LBi,UBi,LBj,UBj,              &
     &                         Si(:,:,ksa_new,isIhis))
        CALL exchange_r2d_tile(ng,tile,LBi,UBi,LBj,UBj,              &
     &                         Si(:,:,ksa_new,isIeis))
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d(ng,tile,iNLM,2,LBi,UBi,LBj,UBj,             &
     &                   NghostPoints,EWperiodic(ng),NSperiodic(ng),    &
     &                   Si(:,:,ksa_new,isIhis),                       &
     &                   Si(:,:,ksa_new,isIeis))
# endif
!
      RETURN
      END SUBROUTINE isf_advection_tile
#endif
      END MODULE isf_advection_mod
