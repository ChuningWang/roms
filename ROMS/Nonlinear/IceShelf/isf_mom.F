#include "cppdefs.h"

      MODULE isf_mom_mod

#ifdef ICESHEET
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
# ifdef SSA_MMS_TESTING
      USE isf_mms_source_mod, ONLY : isf_read_mms_source
# endif
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
      integer :: i, j, iter, ii, jj
      integer :: ksa_old, ksa_new
!
      real(r8) :: cff, eps_e, H_rho, u_mid, v_mid, avg_ff
# ifdef ISF_DEBUG
      real(r8) :: cff1, cff2
# endif
      real(r8) :: u_cg_max
      real(r8) :: rho_w_local
      real(r8) :: rhoi_g
      real(r8) :: B_hardness, C_basal_sq, u_slide
      real(r8) :: gl_half, H_float
      real(r8), parameter :: pi_const = 3.141592653589793_r8
      real(r8) :: T_star, RT_star, A_pat
      real(r8) :: eps_xx, eps_yy, eps_xy, eps_e2
      real(r8) :: u1, u2, u3, u4, v1, v2, v3, v4
      real(r8) :: jac
      real(r8), parameter :: gravity = 9.81_r8
!
!  SSA solver scratch arrays (PSI-point and RHO-point).
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rmask_ist_scr
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: pmask_ist_scr
      real(r8), dimension(LBi:UBi,LBj:UBj) :: pmask_ssa
      real(r8), dimension(LBi:UBi,LBj:UBj) :: float_frac
!
      real(r8), dimension(LBi:UBi,LBj:UBj) :: taudx_psi
      real(r8), dimension(LBi:UBi,LBj:UBj) :: taudy_psi
!
!  NOTE: nuH_rho/beta_rho are passed to mp_exchange2d /
!  exchange_p2d_tile (halo fill for the extended cell assembly in
!  isf_cg_solve.h), which index them on the GLOBAL bounds
!  (LBi:UBi, LBj:UBj); keep them declared on those bounds.
      real(r8), dimension(LBi:UBi,LBj:UBj) :: nuH_rho
      real(r8), dimension(LBi:UBi,LBj:UBj) :: beta_rho
!
!  CG solver work arrays (PSI-point).
!  NOTE: Au/Av/diag_u/diag_v/u_scr/v_scr are passed to mp_exchange2d
!  and exchange_p2d_tile, which index them on the GLOBAL bounds
!  (LBi:UBi, LBj:UBj).  Declaring them on the tile-local IminS:ImaxS
!  range overruns the array on multi-tile MPI decompositions
!  (IminS=Istr-4, ImaxS=Iend+3), corrupting the halo exchange and
!  making the SSA solve diverge under MPI while passing on 1 tile.
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Ru, Rv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Zu, Zv
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Du, Dv
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Au, Av
      real(r8), dimension(LBi:UBi,LBj:UBj) :: diag_u, diag_v
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Ru_old, Rv_old
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Zu_old, Zv_old
      real(r8), dimension(LBi:UBi,LBj:UBj) :: u_scr, v_scr
!
!  CG scalar variables.
!
      real(r8) :: resid, resid_0, alpha, beta_cg
      real(r8) :: dot_zr, dot_zr_old, dot_dad
      real(r8) :: cgtol, nl_resid, nl_resid_0
      real(r8) :: vel_change, vel_norm
      real(r8) :: resid_norm, rhs_norm
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
!  Extended cell-assembly bounds (Istr..Iend + one halo-cell ring).
      integer :: iCstr, iCend, jCstr, jCend
!  MMS debug dump file name (per tile).
      character (len=24) :: dump_fname
!  Smoothing scratch for driving stress regularization.
!  Global bounds: passed to mp_exchange2d / exchange_p2d_tile.
      real(r8), dimension(LBi:UBi,LBj:UBj) :: zs_smooth
      integer :: is_iter
      integer :: dist_margin, di, dj
      real(r8) :: taper_w
      real(r8) :: H11, H21, H12, H22, H_psi
      real(r8) :: dSdx, dSdy, jac_avg
      real(r8), parameter :: sqrt3_inv = 1.0_r8/SQRT(3.0_r8)
# ifdef SSA_MMS_TESTING
!  MMS source terms are written by isf_read_mms_source with GLOBAL
!  indices (LBi:UBi, LBj:UBj).  The arrays must be declared on the
!  global bounds, NOT the tile-local IminS:ImaxS range, otherwise the
!  subroutine overruns the array on multi-tile MPI decompositions
!  (memory corruption -> Picard divergence).
      real(r8), dimension(LBi:UBi,LBj:UBj) :: mms_fu_psi
      real(r8), dimension(LBi:UBi,LBj:UBj) :: mms_fv_psi
# endif
# ifdef DISTRIBUTE
      character(len=3) :: op_handle
# endif
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
!  Use reference seawater density for the flotation criterion.
!  Reuse mod_scalars rho0 (mean Boussinesq density, default 1025.0)
!  for consistency with the isf_masks.F grounding-line criterion.
!
      rho_w_local=rho0
!
!-----------------------------------------------------------------------
!  Prepare the scalar ice state (thickness, temperature, enthalpy)
!  BEFORE any geometry/mask computation.
!
!  The lsold/lsnew swap alternates the two time layers, and the ini
!  file populates layer 1 only (get_state Tindex=1).  Two issues:
!
!  1. The layer now appearing as ksa_new was last filled by the
!     previous call's carry-forward copy, which only touches the tile
!     interior; its ghost ring still holds the zero initial value.
!     pmask_ssa, float_frac and the viscosity all read H from ksa_old
!     including halo cells (the 4 RHO corners of boundary PSI nodes),
!     so the ksa_old halo must be exchanged first; otherwise seam PSI
!     nodes are deactivated (MIN corner H = 0 < min_thickness) and the
!     assembled operator is corrupted under MPI.
!
!  2. ksa_new must be carried forward from ksa_old so that any
!     component NOT updated every SSA call (thickness advection
!     skipped at cold start, or N_thick > 1) does not leave the new
!     layer at zero.  When thickness advection runs it overwrites
!     ksa_new afterwards, so this copy does not disturb the coupling.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile(ng,tile,LBi,UBi,LBj,UBj,                &
     &                         Si(:,:,ksa_old,isIhis))
        CALL exchange_r2d_tile(ng,tile,LBi,UBi,LBj,UBj,                &
     &                         Si(:,:,ksa_old,isIeis))
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d(ng,tile,iNLM,2,LBi,UBi,LBj,UBj,               &
     &                   NghostPoints,EWperiodic(ng),NSperiodic(ng),     &
     &                   Si(:,:,ksa_old,isIhis),                         &
     &                   Si(:,:,ksa_old,isIeis))
# endif
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          Si(i,j,ksa_new,isIhis)=Si(i,j,ksa_old,isIhis)
          Si(i,j,ksa_new,isIeis)=Si(i,j,ksa_old,isIeis)
        END DO
      END DO
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile(ng,tile,LBi,UBi,LBj,UBj,                &
     &                         Si(:,:,ksa_new,isIhis))
        CALL exchange_r2d_tile(ng,tile,LBi,UBi,LBj,UBj,                &
     &                         Si(:,:,ksa_new,isIeis))
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d(ng,tile,iNLM,2,LBi,UBi,LBj,UBj,               &
     &                   NghostPoints,EWperiodic(ng),NSperiodic(ng),     &
     &                   Si(:,:,ksa_new,isIhis),                         &
     &                   Si(:,:,ksa_new,isIeis))
# endif
!
!-----------------------------------------------------------------------
!  Compute ice geometry, grounding mask, and SSA active mask.
!-----------------------------------------------------------------------
!
!  RHO-point grounded ice mask (rmask_ist) from flotation criterion:
!    grounded if rho_i * H > rho_w * h  (h = water depth, positive downward)
!    floating if rho_i * H <= rho_w * h and H > 0
!    ice-free if H <= 0
!
!  NOTE: The SSA solver uses H from ksa_old (the previous SSA state),
!  while isf_masks.F uses H from lsnew (after thickness advection).
!  This intentional asynchrony is part of the sequential coupling
!  design: solve velocity from old geometry, then advect thickness,
!  then update ocean-side masks from new geometry.
!
!  The flotation density rho_w_local is mod_scalars rho0 (mean
!  Boussinesq density, overridable via RHO0 in the ocean input
!  script) -- the same rho0 used by isf_masks.F.  Unlike isf_masks.F,
!  the SSA-side criterion uses the static water depth h (no zeta,
!  Dcrit, eps terms): the SSA solve runs before the barotropic step
!  and must not chase high-frequency sea level oscillations in the
!  grounding-line position.  The ocean-side wet/dry masks are the
!  authority for the cavity (isf_masks.F); SSA float_frac only
!  modulates basal drag in the momentum equation.
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
!  Safety cap on ice thickness (prevent garbage values from
!  initialization or numerical artifacts from affecting the solve).
!  Configurable via max_thickness in isf.in (default 5000 m).
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          IF (Si(i,j,ksa_old,isIhis).gt.max_thickness(ng)) THEN
            Si(i,j,ksa_old,isIhis)=max_thickness(ng)
          END IF
        END DO
      END DO
!
!  SSA active mask (pmask_ssa): a PSI node is active if any of its 4
!  surrounding RHO cells has ice (H > 0).  Inactive nodes are Dirichlet
!  (u = v = 0).
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          IF (Si(i,j,ksa_old,isIhis).gt.0.0_r8 .or.                    &
     &        Si(i-1,j,ksa_old,isIhis).gt.0.0_r8 .or.                  &
     &        Si(i,j-1,ksa_old,isIhis).gt.0.0_r8 .or.                  &
     &        Si(i-1,j-1,ksa_old,isIhis).gt.0.0_r8) THEN
            pmask_ssa(i,j)=1.0_r8
          ELSE
            pmask_ssa(i,j)=0.0_r8
          END IF
        END DO
      END DO
!
!  Deactivate ice-margin PSI nodes where the minimum surrounding H
!  is below min_thickness.  These nodes have poorly-defined stress
!  balance and would produce garbage velocity.  By setting pmask_ssa=0
!  BEFORE the Picard solve they become Dirichlet (u=v=0) throughout
!  the iteration, making the solved problem consistent with the
!  returned solution (previously this was a post-hoc cleanup at
!  isf_cg_solve.h end, which fixed the output but not the solve).
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          IF (pmask_ssa(i,j).eq.1.0_r8) THEN
            IF (MIN(Si(i,j,ksa_old,isIhis),                           &
     &              Si(i-1,j,ksa_old,isIhis),                           &
     &              Si(i,j-1,ksa_old,isIhis),                           &
     &              Si(i-1,j-1,ksa_old,isIhis)).lt.                    &
     &            min_thickness(ng)) THEN
              pmask_ssa(i,j)=0.0_r8
            END IF
          END IF
        END DO
      END DO
!
!  Exchange halo of the SSA active mask.  The element assembly in
!  isf_cg_solve.h scatters each cell's contribution to its 4 corner PSI
!  nodes only when pmask_ssa(corner) == 1; at tile boundaries these
!  corners are the neighbouring tile's halo.  Without an exchange the
!  halo holds uninitialized values, so cross-tile boundary nodes lose
!  their contributions (Au = 0) and the solve diverges under MPI.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,pmask_ssa)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d(ng,tile,iNLM,1,LBi,UBi,LBj,UBj,               &
     &                   NghostPoints,EWperiodic(ng),NSperiodic(ng),     &
     &                   pmask_ssa)
# endif
!
!  Float fraction at RHO points: 0 = grounded, 1 = floating.
!  Used to modulate basal drag (floating ice has no basal drag).
!
!  Sub-grid cosine smoothing near the grounding line (Schoof 2007,
!  Gladstone et al. 2010): float_frac transitions smoothly from 0
!  (grounded) to 1 (floating) over a thickness range defined by
!  gl_float_width.  This improves Picard convergence near the grounding
!  line and gives more physically correct migration rates.
!
!  The transition is centered on the flotation thickness H_float =
!  rho_w * h / rho_i, with half-width dH = gl_float_width / 2.
!  float_frac = 0.5 * (1 + cos(pi * (H - H_float + dH) / (2*dH)))
!  clamped to [0, 1].
!
!  NOTE: no thin-cavity clamp here.  isf_masks.F no longer treats thin
!  cavities as grounded (min_cavity_thickness removed), so the SSA
!  float_frac and the ocean-side masks share the same flotation-only
!  criterion.
!
      gl_half=0.5_r8*gl_float_width(ng)
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          H_rho=Si(i,j,ksa_old,isIhis)
          IF (H_rho.gt.0.0_r8) THEN
            H_float=rho_w_local*h(i,j)/rho_i
            IF (H_rho.ge.H_float+gl_half) THEN
              float_frac(i,j)=0.0_r8         ! fully grounded
            ELSE IF (H_rho.le.H_float-gl_half) THEN
              float_frac(i,j)=1.0_r8         ! floating
            ELSE
!  Transition zone: cosine interpolation from floating (1, thin ice)
!  to grounded (0, thick ice).
              float_frac(i,j)=0.5_r8*(1.0_r8+                          &
     &          COS(pi_const*(H_rho-H_float+gl_half)/                 &
     &              (2.0_r8*gl_half)))
            END IF
          ELSE
            float_frac(i,j)=0.0_r8           ! ice-free
          END IF
        END DO
      END DO
!
!  Exchange halo for float_frac (needed by avg_ff in isf_cg_solve.h
!  which reads neighbor float_frac values at tile boundaries).
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile(ng,tile,LBi,UBi,LBj,UBj,                 &
     &                         float_frac)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d(ng,tile,iNLM,1,LBi,UBi,LBj,UBj,               &
     &                   NghostPoints,EWperiodic(ng),NSperiodic(ng),     &
     &                   float_frac)
# endif
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
# ifdef SSA_MMS_TESTING
!  Read MMS source terms from NetCDF (one-time per SSA call).
      CALL isf_read_mms_source (ng, tile, LBi, UBi, LBj, UBj,          &
     &                          mms_fu_psi, mms_fv_psi)
#  ifdef ISF_DEBUG
!  Report initial velocity (ksa_old - the layer read from the ini file),
!  mms and thickness magnitudes (tile 0).
      IF (tile.eq.0) THEN
        u_slide=0.0_r8; cff=0.0_r8; cff1=0.0_r8
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            cff=MAX(cff,ABS(Si(i,j,ksa_old,isIhis)))
          END DO
        END DO
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            IF (pmask_ssa(i,j).eq.1.0_r8) THEN
              u_slide=MAX(u_slide,SQRT(Si(i,j,ksa_old,isIuis)**2+      &
     &                                 Si(i,j,ksa_old,isIvis)**2))
              cff1=MAX(cff1,ABS(mms_fu_psi(i,j)))
            END IF
          END DO
        END DO
        WRITE(*,'(A,ES12.4,A,ES12.4,A,ES12.4)')                        &
     &      '  MMS init: max|u_old|=',u_slide,' max|mms_fu|=',cff1,   &
     &      ' max|H|=',cff
      END IF
#  endif
# endif
!
!-----------------------------------------------------------------------
!  SSA Picard iteration: solve A(u)*u = taud using matrix-free PCG.
!  Active under ICESHEET_SSA (full ice-sheet) or SSA_MMS_TESTING
!  (standalone solver verification).
!-----------------------------------------------------------------------
!
#  if defined ICESHEET_SSA || defined SSA_MMS_TESTING
!  NOTE: the halo exchange of the ksa_new velocity estimate happens
!  INSIDE isf_cg_solve.h, AFTER the initial-guess copy (ksa_old ->
!  ksa_new) and the Dirichlet BC enforcement.  Exchanging before the
!  copy leaves the halo stale (zero) because the copy only touches the
!  tile interior, corrupting Au/Av at tile boundaries under MPI.
#   include <isf_cg_solve.h>
#  endif
!
!  Exchange halo of updated velocity so that downstream routines
!  (thickness advection, output) see consistent tile-boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,              &
     &                         Si(:,:,ksa_new,isIuis))
        CALL exchange_p2d_tile(ng,tile,LBi,UBi,LBj,UBj,              &
     &                         Si(:,:,ksa_new,isIvis))
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d(ng,tile,iNLM,2,LBi,UBi,LBj,UBj,             &
     &                   NghostPoints,EWperiodic(ng),NSperiodic(ng),   &
     &                   Si(:,:,ksa_new,isIuis),                       &
     &                   Si(:,:,ksa_new,isIvis))
# endif
!
      RETURN
      END SUBROUTINE isf_mom_tile
#endif
      END MODULE isf_mom_mod
