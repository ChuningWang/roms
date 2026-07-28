#include "cppdefs.h"

      MODULE isf_mom_mod

#ifdef ICESHELF_MOM
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2024 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!========================================== Benjamin K. Galton-Fenzi ===
!                                                                      !
!  This module computes the ice sheet/shelf momentum equation using    !
!  the Shallow Shelf Approximation (SSA).                              !
!                                                                      !
!  The SSA solver uses a Q1 bilinear finite element discretization     !
!  on a B-grid (velocity at PSI points, scalars at RHO points) with    !
!  matrix-free preconditioned conjugate gradient (PCG) linear solve    !
!  and Picard nonlinear iteration.                                     !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
      PRIVATE
      PUBLIC  :: isf_mom
!
!  CG global dot product shared variables (for OMP critical reduction).
!
      integer :: isf_tile_count
      real(r8) :: isf_dot_prod
!$OMP THREADPRIVATE (isf_tile_count, isf_dot_prod)
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
      CALL wclock_on (ng, model, 42, __LINE__, MyFile)
# endif
      CALL isf_mom_tile (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nrhs(ng),                                      &
     &                   GRID(ng) % pm,                                 &
     &                   GRID(ng) % pn,                                 &
     &                   GRID(ng) % zice,                               &
     &                   GRID(ng) % h,                                  &
     &                   GRID(ng) % rmask_isf,                          &
     &                   GRID(ng) % pmask,                              &
     &                   GRID(ng) % rmask,                              &
     &                   ISF(ng) % Si,                                  &
     &                   ISF(ng) % Fi)
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
      RETURN
      END SUBROUTINE isf_mom
!
!***********************************************************************
      SUBROUTINE isf_mom_tile (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs_val,                                &
     &                         pm, pn,                                  &
     &                         zice, h,                                 &
     &                         rmask_isf, pmask, rmask,                 &
     &                         Si, Fi)
!***********************************************************************
!
      USE mod_kinds
      USE mod_param
      USE mod_scalars
      USE mod_parallel
      USE mod_stepping
      USE mod_iounits
      USE mod_isf
!
      USE bc_2d_mod
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_reduce
# endif
      USE mp_exchange_mod, ONLY : mp_exchange2d
      USE exchange_2d_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs_val
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)    :: pm(LBi:,LBj:)
      real(r8), intent(in)    :: pn(LBi:,LBj:)
      real(r8), intent(in)    :: zice(LBi:,LBj:)
      real(r8), intent(in)    :: h(LBi:,LBj:)
      real(r8), intent(in)    :: rmask_isf(LBi:,LBj:)
      real(r8), intent(in)    :: pmask(LBi:,LBj:)
      real(r8), intent(in)    :: rmask(LBi:,LBj:)
      real(r8), intent(inout) :: Si(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: Fi(LBi:,LBj:,:)
# else
      real(r8), intent(in)    :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: zice(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: rmask_isf(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Si(LBi:UBi,LBj:UBj,2,nIsfS)
      real(r8), intent(inout) :: Fi(LBi:UBi,LBj:UBj,nIsfF)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, iter
      integer :: ksa_old, ksa_new
!
      real(r8) :: cff, eps_e, H_rho, u_mid, v_mid
      real(r8) :: rho_w_local
      real(r8) :: rhoi_g
      real(r8) :: B_glen_sq, C_basal_sq, u_slide
      real(r8) :: T_star, RT_star, A_glen
      real(r8) :: eps_xx, eps_yy, eps_xy, eps_e2
      real(r8) :: u1, u2, u3, u4, v1, v2, v3, v4
      real(r8) :: jac
      real(r8), parameter :: gravity = 9.81_r8
!
!  SSA solver scratch arrays (PSI-point and RHO-point).
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rmask_ist_scr
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: pmask_ist_scr
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: pmask_ssa
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: float_frac
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: taudx_psi
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: taudy_psi
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: nuH_rho
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: beta_rho
!
!  CG solver work arrays (PSI-point).
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Ru, Rv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Zu, Zv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Du, Dv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Au, Av
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: diag_u, diag_v
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Ru_old, Rv_old
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Zu_old, Zv_old
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: u_scr, v_scr
!
!  CG scalar variables.
!
      real(r8) :: resid, resid_0, alpha, beta_cg
      real(r8) :: dot_zr, dot_zr_old, dot_dad
      real(r8) :: cgtol, nl_resid, nl_resid_0
      real(r8) :: vel_change, vel_norm
      integer :: cg_iter, nl_iter
      logical :: converged
!
!  FE shape function derivatives and quadrature data.
!
      real(r8), dimension(4,4,2) :: DPhi
      real(r8), dimension(4,4) :: phi_val
      real(r8), dimension(4) :: Xquad, jacq
      real(r8) :: dx_rho, dy_rho, uq, vq, exx, eyy, exy
      real(r8) :: dphi_x, dphi_y, phiv
      real(r8) :: dSdx_q, dSdy_q, H_q
      real(r8) :: zs_corner(4), H_corner(4), zs_sum, zs_mean
      integer :: nq, k_node, iq, jq, n_ice
      integer :: i1, i2, j1, j2
!  Smoothing scratch for driving stress regularization.
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zs_smooth
      integer :: n_smooth, is_iter
      real(r8), parameter :: sqrt3_inv = 1.0_r8/SQRT(3.0_r8)
# ifdef DISTRIBUTE
      character(len=3) :: op_handle
# endif
!
!  Pointer aliases for velocity arrays (for halo exchange calls).
!
      real(r8) :: dummy_r8
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", isf_mom_tile"
!
# include "set_bounds.h"
!
      SourceFile=MyFile
!
!  Set time level indices.  The lsold/lsnew swap was performed in
!  iceshelf.F before calling isf_mom.  ksa_old = previous (initial guess),
!  ksa_new = current (solution to be computed).
!
      ksa_old=lsold(ng)
      ksa_new=lsnew(ng)
!
!  Use rho0 as seawater density for flotation criterion.
!
      rho_w_local=REAL(rho0,r8)
!
!-----------------------------------------------------------------------
!  Compute ice geometry, grounding mask, and SSA active mask.
!-----------------------------------------------------------------------
!
!  RHO-point grounded ice mask (rmask_ist) from flotation criterion:
!    grounded if rho_i * H > rho_w * D  (D = water depth = -h)
!    floating if rho_i * H <= rho_w * D and H > 0
!    ice-free if H <= 0
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          H_rho=Si(i,j,ksa_old,isIhis)
          IF ((H_rho.gt.0.0_r8).and.                                  &
     &        (rho_i*H_rho.gt.rho_w_local*h(i,j))) THEN
            rmask_ist_scr(i,j)=1.0_r8         ! grounded
          ELSE
            rmask_ist_scr(i,j)=0.0_r8         ! floating or ice-free
          END IF
        END DO
      END DO
!
!  PSI-point grounded ice mask (pmask_ist): a PSI node is grounded if
!  any of its 4 surrounding RHO cells is grounded.
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          pmask_ist_scr(i,j)=MAX(rmask_ist_scr(i,j),                   &
     &                           rmask_ist_scr(i-1,j),                 &
     &                           rmask_ist_scr(i,j-1),                 &
     &                           rmask_ist_scr(i-1,j-1))
        END DO
      END DO
!
!  SSA active mask (pmask_ssa): a PSI node is active if any of its 4
!  surrounding RHO cells has ice (H > 0).  Inactive nodes are Dirichlet
!  (u = v = 0).
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          IF (Si(i,j,ksa_old,isIhis).gt.0.0_r8 .or.                       &
     &        Si(i-1,j,ksa_old,isIhis).gt.0.0_r8 .or.                     &
     &        Si(i,j-1,ksa_old,isIhis).gt.0.0_r8 .or.                     &
     &        Si(i-1,j-1,ksa_old,isIhis).gt.0.0_r8) THEN
            pmask_ssa(i,j)=1.0_r8
          ELSE
            pmask_ssa(i,j)=0.0_r8
          END IF
        END DO
      END DO
!
!  Float fraction at RHO points: 0 = grounded, 1 = floating.
!  Used to modulate basal drag (floating ice has no basal drag).
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          H_rho=Si(i,j,ksa_old,isIhis)
          IF (H_rho.gt.0.0_r8) THEN
            IF (rmask_ist_scr(i,j).eq.1.0_r8) THEN
              float_frac(i,j)=0.0_r8           ! grounded
            ELSE
              float_frac(i,j)=1.0_r8           ! floating
            END IF
          ELSE
            float_frac(i,j)=0.0_r8             ! ice-free
          END IF
        END DO
      END DO
!
!  Compute ice surface elevation (icIzsf) and bed elevation (icIzbd)
!  at RHO points.
!
!    grounded: S = -h + H,  B = -h
!    floating: S = (1 - rho_i/rho_w) * H,  B = -(rho_i/rho_w) * H
!    ice-free: S = 0, B = -h
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          H_rho=Si(i,j,ksa_old,isIhis)
          IF (H_rho.gt.0.0_r8) THEN
            IF (rmask_ist_scr(i,j).eq.1.0_r8) THEN
              Fi(i,j,icIzsf)=-h(i,j)+H_rho
              Fi(i,j,icIzbd)=-h(i,j)
            ELSE
              Fi(i,j,icIzsf)=(1.0_r8-rho_i/rho_w_local)*H_rho
              Fi(i,j,icIzbd)=-(rho_i/rho_w_local)*H_rho
            END IF
          ELSE
            Fi(i,j,icIzsf)=0.0_r8
            Fi(i,j,icIzbd)=-h(i,j)
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Precompute Q1 bilinear shape function derivatives at 2x2 Gauss pts.
!  Node numbering: 1=(i,j), 2=(i+1,j), 3=(i,j+1), 4=(i+1,j+1)
!  Reference element: xi in [-1,1], eta in [-1,1]
!  Used by both driving stress FE assembly and CG stiffness assembly.
!-----------------------------------------------------------------------
!
      Xquad(1)=0.5_r8*(1.0_r8-sqrt3_inv)
      Xquad(2)=0.5_r8*(1.0_r8+sqrt3_inv)
!
      DO jq=1,2
        DO iq=1,2
          nq=2*(jq-1)+iq
          jacq(nq)=0.25_r8
!
!  Node 1: (xi=-1, eta=-1)
!
          DPhi(1,nq,1)=-0.25_r8*(1.0_r8-Xquad(jq))*2.0_r8
          DPhi(1,nq,2)=-0.25_r8*(1.0_r8-Xquad(iq))*2.0_r8
          phi_val(1,nq)=(1.0_r8-Xquad(iq))*(1.0_r8-Xquad(jq))
!
!  Node 2: (xi=+1, eta=-1)
!
          DPhi(2,nq,1)= 0.25_r8*(1.0_r8-Xquad(jq))*2.0_r8
          DPhi(2,nq,2)=-0.25_r8*Xquad(iq)*2.0_r8
          phi_val(2,nq)=Xquad(iq)*(1.0_r8-Xquad(jq))
!
!  Node 3: (xi=-1, eta=+1)
!
          DPhi(3,nq,1)=-0.25_r8*Xquad(jq)*2.0_r8
          DPhi(3,nq,2)= 0.25_r8*(1.0_r8-Xquad(iq))*2.0_r8
          phi_val(3,nq)=(1.0_r8-Xquad(iq))*Xquad(jq)
!
!  Node 4: (xi=+1, eta=+1)
!
          DPhi(4,nq,1)= 0.25_r8*Xquad(jq)*2.0_r8
          DPhi(4,nq,2)= 0.25_r8*Xquad(iq)*2.0_r8
          phi_val(4,nq)=Xquad(iq)*Xquad(jq)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute driving stress RHS via FE assembly (taudx_psi, taudy_psi).
!-----------------------------------------------------------------------
!
# ifdef ICESHEET_SSA
#   include <isf_driving_stress.h>
#  endif
!
!-----------------------------------------------------------------------
!  Compute effective viscosity (nuH) and basal drag (beta) at RHO pts.
!-----------------------------------------------------------------------
!
#  ifdef ICESHEET_SSA
#   include <isf_viscosity.h>
#  endif
!
!-----------------------------------------------------------------------
!  SSA Picard iteration: solve A(u)*u = taud using matrix-free PCG.
!-----------------------------------------------------------------------
!
#  ifdef ICESHEET_SSA
#   include <isf_cg_solve.h>
#  endif
!
!  Copy solution from ksa_new to output (lsnew was set by driver).
!  The CG solver writes directly into Si(:,:,ksa_new,isIuis/isIvis).
!
      RETURN
      END SUBROUTINE isf_mom_tile
#endif
      END MODULE isf_mom_mod
