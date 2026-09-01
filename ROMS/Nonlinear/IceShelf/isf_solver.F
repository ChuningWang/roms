#include "cppdefs.h"

      MODULE isf_solver_mod
#ifdef ICESHEET_SSA
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module solves the ice sheet/shelf momentum equation using      !
!  the Shallow Shelf Approximation (SSA).  It is a self-contained      !
!  solver designed to be called from the ice momentum driver           !
!  (isf_mom):  the driver computes the ice geometry, masks and         !
!  boundary fields, this module computes the horizontal velocity.      !
!                                                                      !
!  Discretization:  Q1 bilinear finite elements on the ROMS B-grid     !
!  (velocity at PSI points, scalars at RHO points) with 2x2 Gauss      !
!  quadrature, following MITgcm STREAMICE.  The nonlinear Glen-law     !
!  viscosity is resolved by Picard iteration; each linear solve is     !
!  a matrix-free preconditioned conjugate gradient (PCG) with a        !
!  Jacobi (diagonal) preconditioner.                                   !
!                                                                      !
!  The basal drag coefficient beta is an INPUT:  it is taken from      !
!  Fi(:,:,icIbta), computed by isf_stress from the ocean state         !
!  averaged over the SSA coupling window (Weertman or regularized      !
!  Coulomb sliding law).  Freezing beta at the coupling-window value   !
!  is the intended frozen-coefficient (Picard) linearization of the    !
!  drag term; only the viscosity is updated inside the Picard loop.    !
!                                                                      !
!  Input contract (caller must fill before calling isf_ssa):           !
!                                                                      !
!    Si(:,:,ksa_old,isIhis)        ice thickness [m], halo exchanged   !
!    Si(:,:,ksa_old,isIuis/isIvis) velocity initial guess [m s-1]      !
!    Fi(:,:,icItis)                ice temperature [degC]              !
!    Fi(:,:,icIzsf)                ice surface elevation [m]           !
!    Fi(:,:,icIbta)                basal drag coefficient (isf_stress) !
!                                                                      !
!  Output:                                                             !
!                                                                      !
!    Si(:,:,ksa_new,isIuis/isIvis) solved SSA velocity [m s-1],        !
!                                  halo exchanged                      !
!    Fi(:,:,icInuh)                depth-integrated viscosity          !
!    Fi(:,:,icItdx/icItdy)         driving stress diagnostics          !
!                                                                      !
!  Lateral boundaries are LBC-aware:  closed edges impose zero         !
!  velocity, clamped edges fix the velocity to ISF_LOBC data,          !
!  gradient edges are free.  Dirichlet enforcement acts on the         !
!  normal component of each edge (u on west, v on south/north);        !
!  the eastern (calving front) edge is always left free.  The          !
!  ice-front pressure boundary term is future work.                    !
!                                                                      !
!  MPI notes:  PSI tile ranges do not overlap (get_bounds.F:           !
!  IstrP=my_Istr everywhere, IendP=my_Iend except on the eastern       !
!  domain edge), so a seam PSI node is owned by exactly one tile       !
!  while the RHO cells feeding its FE stencil lie on both sides of     !
!  the seam.  The element assembly therefore integrates one extra      !
!  ring of halo cells along the western/southern tile seams (bounds    !
!  iCstr/jCstr below; the legacy isf_cg_solve.h documents this         !
!  strategy but assembles Istr..Iend only, leaving seam node rows      !
!  incomplete under multi-tile MPI) and scatters only to owned PSI     !
!  nodes (ownership filter at every scatter).  mp_exchange2d fills     !
!  corner ghost points (E/W segments span the full J range), so the    !
!  halo ring reads valid exchanged data.                               !
!                                                                      !
!  Parallel decomposition:  correct in serial and MPI (DISTRIBUTE).    !
!  Pure OpenMP is NOT supported:  the CG iteration state lives in      !
!  tile-local automatic arrays which cannot be merged across           !
!  shared-memory tiles (same limitation as the legacy solver).         !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    MacAyeal (1989), J. Geophys. Res. 94, 4071-4087.         (SSA)    !
!    Goldberg (2011), J. Comput. Phys. 230, 6513-6525.                 !
!    Bueler & Brown (2009), J. Glaciol. 55, 2511-2529.                 !
!      (Paterson temperature-dependent Glen rate factor)               !
!    MITgcm STREAMICE streamice_cg_solve.F / streamice_visc_beta.F     !
!      (Q1 FE assembly, matrix-free PCG, Picard iteration).            !
!    PISM compute_nuH_everywhere (nuH regularization, thin-ice         !
!      strength extension).                                            !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
      PRIVATE
      PUBLIC  :: isf_ssa
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE isf_ssa (ng, tile, model)
!***********************************************************************
!
!  Outer interface of the SSA solver.  Sets tile bounds and calls
!  the tile computational routine.
!
      USE mod_param
      USE mod_grid
      USE mod_stepping
      USE mod_isf, ONLY : ISF
!
      implicit none
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
      CALL isf_ssa_tile (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   lsold(ng), lsnew(ng),                          &
     &                   GRID(ng) % pm,                                 &
     &                   GRID(ng) % pn,                                 &
     &                   ISF(ng) % Si,                                  &
     &                   ISF(ng) % Fi)
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE isf_ssa
!
!***********************************************************************
      SUBROUTINE isf_ssa_tile (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         ksa_old, ksa_new,                        &
     &                         pm, pn,                                  &
     &                         Si, Fi)
!***********************************************************************
!
!  Tile-level SSA solve.  Orchestrates:  active-node mask, initial
!  guess (copy + basal-drag warm start), driving stress RHS, and the
!  Picard/PCG iteration.  All kernels (viscosity, matrix action,
!  Jacobi diagonal, driving stress, LBC enforcement) live in the
!  private subroutines below;  this routine owns halo exchanges,
!  global reductions and convergence control.
!
      USE mod_kinds
      USE mod_param
      USE mod_scalars
      USE mod_parallel
      USE mod_iounits
      USE mod_isf
!
      USE exchange_2d_mod, ONLY : exchange_p2d_tile, exchange_r2d_tile
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_reduce
# endif
      USE mp_exchange_mod, ONLY : mp_exchange2d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: ksa_old, ksa_new
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
      integer :: i, j, cg_iter, nl_iter
!
!  LBC-aware lateral boundary flags:  *_fix marks a Dirichlet
!  constrained edge (closed or clamped);  *_clamp marks an edge with
!  a prescribed boundary value from ISF_LOBC (clamped).  Gradient
!  edges are free.  See ssa_lbc_velocity / ssa_lbc_zero.
!
      logical :: Lu_west_fix, Lu_west_clamp
      logical :: Lv_south_fix, Lv_south_clamp
      logical :: Lv_north_fix, Lv_north_clamp
      logical :: converged
!
      real(r8) :: cff, u_slide, beta_avg, jac_avg
# ifdef ISF_DEBUG
      real(r8) :: cff1
# endif
!
!  CG scalar variables.
!
      real(r8) :: resid, resid_0, alpha, beta_cg
      real(r8) :: dot_zr, dot_zr_old, dot_dad
      real(r8) :: cgtol, nl_resid
      real(r8) :: vel_change, vel_norm
      real(r8) :: resid_norm, rhs_norm
!
!  FE shape function tables (Q1 bilinear, 2x2 Gauss quadrature).
!  Node numbering:  1=(i,j), 2=(i+1,j), 3=(i,j+1), 4=(i+1,j+1).
!  Reference element [0,1]x[0,1];  physical derivatives carry the
!  metric factor 2/dx (2/dy) and jac = dx*dy = 1/(pm*pn).
!
      real(r8), dimension(4,4,2) :: DPhi
      real(r8), dimension(4,4)   :: phi_tab
      real(r8), dimension(4)     :: Xquad
      real(r8), parameter :: sqrt3_inv = 1.0_r8/SQRT(3.0_r8)
      integer :: iq, jq
!
!  SSA solver fields.  Arrays passed to exchange / mp_exchange are
!  declared on the GLOBAL bounds (LBi:UBi, LBj:UBj).  NOTE: Ru/Rv
!  (and Zu/Zv) are also declared on the global bounds, matching the
!  explicit-shape dummy declarations of ssa_lbc_zero:  the tile
!  scratch bounds (IminS:ImaxS, JminS:JmaxS) differ from LBi:UBi by
!  the ghost margin, and sequence association of an explicit-shape
!  dummy with a narrower-width actual would map the fixed-node
!  zeroing onto the wrong elements.
!
      real(r8), dimension(LBi:UBi,LBj:UBj) :: pmask_ssa
      real(r8), dimension(LBi:UBi,LBj:UBj) :: taudx_psi
      real(r8), dimension(LBi:UBi,LBj:UBj) :: taudy_psi
      real(r8), dimension(LBi:UBi,LBj:UBj) :: nuH_rho
      real(r8), dimension(LBi:UBi,LBj:UBj) :: beta_rho
!
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Ru, Rv
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Zu, Zv
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Du, Dv
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Au, Av
      real(r8), dimension(LBi:UBi,LBj:UBj) :: diag_u, diag_v
      real(r8), dimension(LBi:UBi,LBj:UBj) :: u_scr, v_scr
!
!  Extended cell-assembly bounds (see module header):  own RHO cells
!  Istr..Iend plus the western/southern halo ring that feeds the
!  seam PSI nodes owned by this tile.  The ring is skipped along
!  physical (non-periodic) domain edges where the halo cells do not
!  exist.
!
      integer :: iCstr, iCend, jCstr, jCend
!
# ifdef DISTRIBUTE
      character (len=3) :: op_handle
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", isf_ssa_tile"
!
# include "set_bounds.h"
!
      SourceFile=MyFile
!
!=======================================================================
!  Stage 1:  SSA active-node mask (PSI points).
!=======================================================================
!
!  Exchange the thickness halo first:  the mask stencil reads the 4
!  RHO cells surrounding each PSI node, including halo cells at MPI
!  tile seams (idempotent if the caller already exchanged it).
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,ksa_old,isIhis))
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Si(:,:,ksa_old,isIhis))
# endif
!
!  A PSI node is active if any of its 4 surrounding RHO cells has
!  ice (H > 0).  Margin nodes whose minimum surrounding H is below
!  min_thickness are deactivated:  their stress balance is poorly
!  defined, and setting pmask_ssa=0 BEFORE the Picard solve makes
!  them Dirichlet (u=v=0) throughout the iteration so that the
!  solved problem is consistent with the returned solution.
!
      DO j=LBj,UBj
        DO i=LBi,UBi
          pmask_ssa(i,j)=0.0_r8
        END DO
      END DO
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          IF (Si(i,j,ksa_old,isIhis).gt.0.0_r8 .or.                     &
     &        Si(i-1,j,ksa_old,isIhis).gt.0.0_r8 .or.                   &
     &        Si(i,j-1,ksa_old,isIhis).gt.0.0_r8 .or.                   &
     &        Si(i-1,j-1,ksa_old,isIhis).gt.0.0_r8) THEN
            pmask_ssa(i,j)=1.0_r8
          END IF
        END DO
      END DO
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          IF (pmask_ssa(i,j).eq.1.0_r8) THEN
            IF (MIN(Si(i,j,ksa_old,isIhis),                             &
     &              Si(i-1,j,ksa_old,isIhis),                           &
     &              Si(i,j-1,ksa_old,isIhis),                           &
     &              Si(i-1,j-1,ksa_old,isIhis)).lt.                     &
     &            min_thickness(ng)) THEN
              pmask_ssa(i,j)=0.0_r8
            END IF
          END IF
        END DO
      END DO
!
!  Exchange the mask halo:  the element assembly below tests
!  pmask_ssa at every scatter target, and the CG boundary treatment
!  requires consistent values at tile seams.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, pmask_ssa)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pmask_ssa)
# endif
!
!=======================================================================
!  Stage 2:  LBC-aware Dirichlet flags.
!=======================================================================
!
!  Mixed is not a supported ice-sheet boundary type (Cla/Clo/Gra
!  only).  The eastern (calving front) edge is always free.
!
      Lu_west_fix   = LBC(iwest ,ibISF(isIuis),ng)%closed .or.          &
     &                LBC(iwest ,ibISF(isIuis),ng)%clamped
      Lu_west_clamp = LBC(iwest ,ibISF(isIuis),ng)%clamped
      Lv_south_fix  = LBC(isouth,ibISF(isIvis),ng)%closed .or.          &
     &                LBC(isouth,ibISF(isIvis),ng)%clamped
      Lv_south_clamp= LBC(isouth,ibISF(isIvis),ng)%clamped
      Lv_north_fix  = LBC(inorth,ibISF(isIvis),ng)%closed .or.          &
     &                LBC(inorth,ibISF(isIvis),ng)%clamped
      Lv_north_clamp= LBC(inorth,ibISF(isIvis),ng)%clamped
!
!=======================================================================
!  Stage 3:  basal drag coefficient (input from isf_stress).
!=======================================================================
!
!  Copy beta from Fi(:,:,icIbta) into the local work array and
!  exchange its halo:  the element assembly and the warm start read
!  beta at the RHO cells surrounding PSI nodes, including halo ring
!  cells.  The zero-initialization also protects halo corner cells
!  that the exchanges may not fill.
!
      DO j=LBj,UBj
        DO i=LBi,UBi
          beta_rho(i,j)=0.0_r8
          nuH_rho(i,j)=0.0_r8
        END DO
      END DO
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          beta_rho(i,j)=Fi(i,j,icIbta)
        END DO
      END DO
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, beta_rho)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    beta_rho)
# endif
!
!=======================================================================
!  Stage 4:  Q1 shape function tables at the 2x2 Gauss points.
!=======================================================================
!
      Xquad(1)=0.5_r8*(1.0_r8-sqrt3_inv)
      Xquad(2)=0.5_r8*(1.0_r8+sqrt3_inv)
!
      DO jq=1,2
        DO iq=1,2
!
!  Node 1:  (xi=0, eta=0)
!
          DPhi(1,2*(jq-1)+iq,1)=-0.5_r8*(1.0_r8-Xquad(jq))
          DPhi(1,2*(jq-1)+iq,2)=-0.5_r8*(1.0_r8-Xquad(iq))
          phi_tab(1,2*(jq-1)+iq)=(1.0_r8-Xquad(iq))*(1.0_r8-Xquad(jq))
!
!  Node 2:  (xi=1, eta=0)
!
          DPhi(2,2*(jq-1)+iq,1)= 0.5_r8*(1.0_r8-Xquad(jq))
          DPhi(2,2*(jq-1)+iq,2)=-0.5_r8*Xquad(iq)
          phi_tab(2,2*(jq-1)+iq)=Xquad(iq)*(1.0_r8-Xquad(jq))
!
!  Node 3:  (xi=0, eta=1)
!
          DPhi(3,2*(jq-1)+iq,1)=-0.5_r8*Xquad(jq)
          DPhi(3,2*(jq-1)+iq,2)= 0.5_r8*(1.0_r8-Xquad(iq))
          phi_tab(3,2*(jq-1)+iq)=(1.0_r8-Xquad(iq))*Xquad(jq)
!
!  Node 4:  (xi=1, eta=1)
!
          DPhi(4,2*(jq-1)+iq,1)= 0.5_r8*Xquad(jq)
          DPhi(4,2*(jq-1)+iq,2)= 0.5_r8*Xquad(iq)
          phi_tab(4,2*(jq-1)+iq)=Xquad(iq)*Xquad(jq)
        END DO
      END DO
!
!  Extended cell-assembly bounds (see module header).
!
      iCstr=Istr
      IF (.not.DOMAIN(ng)%Western_Edge(tile).or.EWperiodic(ng)) THEN
        iCstr=Istr-1
      END IF
      jCstr=Jstr
      IF (.not.DOMAIN(ng)%Southern_Edge(tile).or.NSperiodic(ng)) THEN
        jCstr=Jstr-1
      END IF
      iCend=Iend
      jCend=Jend
!
!=======================================================================
!  Stage 5:  initial velocity guess.
!=======================================================================
!
!  Copy the previous SSA velocity (ksa_old) as the initial guess,
!  zero the inactive nodes and enforce the Dirichlet lateral
!  boundaries.
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          Si(i,j,ksa_new,isIuis)=Si(i,j,ksa_old,isIuis)
          Si(i,j,ksa_new,isIvis)=Si(i,j,ksa_old,isIvis)
          IF (pmask_ssa(i,j).eq.0.0_r8) THEN
            Si(i,j,ksa_new,isIuis)=0.0_r8
            Si(i,j,ksa_new,isIvis)=0.0_r8
          END IF
        END DO
      END DO
!
      CALL ssa_lbc_velocity (ng, ksa_new,                               &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Istr, Jstr, Jend,                          &
     &                       IstrP, IendP, JstrP, JendP,                &
     &                       Lu_west_fix, Lu_west_clamp,                &
     &                       Lv_south_fix, Lv_south_clamp,              &
     &                       Lv_north_fix, Lv_north_clamp,              &
     &                       Si)
!
!  Exchange the halo of the initial guess AFTER the copy and the
!  Dirichlet enforcement:  the element assembly reads all 4 PSI
!  corners of every integrated cell, including halo nodes.
!
      CALL ssa_exchange_uv (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IstrP, IendP, JstrP, JendP,                 &
     &                      ksa_new, Si)
!
!=======================================================================
!  Stage 6:  driving stress RHS (FE-consistent nodal forces).
!=======================================================================
!
      CALL ssa_driving_stress (ng, tile, model,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Istr, Iend, Jstr, Jend,                   &
     &                        IstrP, IendP, JstrP, JendP,               &
     &                        ksa_old,                                  &
     &                        pm, pn, pmask_ssa,                        &
     &                        Si, Fi,                                   &
     &                        taudx_psi, taudy_psi)
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          Fi(i,j,icItdx)=taudx_psi(i,j)
          Fi(i,j,icItdy)=taudy_psi(i,j)
        END DO
      END DO
!
!=======================================================================
!  Stage 7:  basal-drag warm start.
!=======================================================================
!
!  If the initial velocity is near zero, invert the local basal drag
!  balance (tau_d = beta * u) to provide a physically reasonable
!  first guess.  This avoids the Picard "high-viscosity trap"
!  (zero velocity -> eps_e = eps_glen_min -> huge nuH -> tiny
!  velocity -> convergence to the wrong fixed point).  Because beta
!  is the frozen coupling-window coefficient, the inversion is
!  valid for any sliding exponent and for both the Weertman and
!  Coulomb laws.  On floating ice (beta = 0) the warm start is
!  skipped and the previous velocity field is kept.
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          IF (pmask_ssa(i,j).eq.1.0_r8) THEN
            u_slide=SQRT(Si(i,j,ksa_new,isIuis)**2+                     &
     &                   Si(i,j,ksa_new,isIvis)**2)
            IF (u_slide.lt.1.0E-5_r8) THEN
              beta_avg=0.25_r8*(beta_rho(i,j)+                          &
     &                          beta_rho(i-1,j)+                        &
     &                          beta_rho(i,j-1)+                        &
     &                          beta_rho(i-1,j-1))
              IF (beta_avg.gt.1.0E-10_r8) THEN
                jac_avg=0.25_r8*(                                       &
     &            1.0_r8/(pm(i-1,j-1)*pn(i-1,j-1))+                     &
     &            1.0_r8/(pm(i,  j-1)*pn(i,  j-1))+                     &
     &            1.0_r8/(pm(i-1,j  )*pn(i-1,j  ))+                     &
     &            1.0_r8/(pm(i,  j  )*pn(i,  j  )))
                cff=1.0_r8/(beta_avg*jac_avg*0.25_r8)
                Si(i,j,ksa_new,isIuis)=taudx_psi(i,j)*cff
                Si(i,j,ksa_new,isIvis)=taudy_psi(i,j)*cff
              END IF
            END IF
          END IF
        END DO
      END DO
!
!  Re-enforce the Dirichlet boundaries and refresh the velocity
!  halo after the warm start (it modifies interior nodes only).
!
      CALL ssa_lbc_velocity (ng, ksa_new,                               &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Istr, Jstr, Jend,                          &
     &                       IstrP, IendP, JstrP, JendP,                &
     &                       Lu_west_fix, Lu_west_clamp,                &
     &                       Lv_south_fix, Lv_south_clamp,              &
     &                       Lv_north_fix, Lv_north_clamp,              &
     &                       Si)
!
      CALL ssa_exchange_uv (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IstrP, IendP, JstrP, JendP,                 &
     &                      ksa_new, Si)
!
!=======================================================================
!  Stage 8:  initialize CG work arrays.
!=======================================================================
!
      DO j=LBj,UBj
        DO i=LBi,UBi
          Au(i,j)=0.0_r8
          Av(i,j)=0.0_r8
          diag_u(i,j)=0.0_r8
          diag_v(i,j)=0.0_r8
          u_scr(i,j)=0.0_r8
          v_scr(i,j)=0.0_r8
          Du(i,j)=0.0_r8
          Dv(i,j)=0.0_r8
          Ru(i,j)=0.0_r8
          Rv(i,j)=0.0_r8
          Zu(i,j)=0.0_r8
          Zv(i,j)=0.0_r8
        END DO
      END DO
!
!=======================================================================
!  Stage 9:  Picard iteration.
!=======================================================================
!
!  Each iteration freezes the viscosity at the current velocity,
!  solves the linearized SSA with matrix-free PCG, applies the
!  robustness measures (speed cap, relaxation) and checks the dual
!  convergence criterion.
!
      cgtol=cg_tol(ng)
      converged=.FALSE.
!
      DO nl_iter=1,nl_max_iter(ng)
!
!  Effective viscosity at RHO points from the current velocity
!  estimate (Glen law + Paterson rate factor + regularization).
!
        CALL ssa_viscosity (ng,                                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Istr, Iend, Jstr, Jend,                     &
     &                      ksa_old, ksa_new,                           &
     &                      pm, pn,                                     &
     &                      Si, Fi,                                     &
     &                      nuH_rho)
!
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Fi(i,j,icInuh)=nuH_rho(i,j)
          END DO
        END DO
!
!  Exchange the nuH halo:  the assembly below integrates the halo
!  ring cells (iCstr/jCstr), whose viscosity was computed on the
!  neighbouring tile.
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_p2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, nuH_rho)
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      nuH_rho)
# endif
!
!  Jacobi preconditioner:  diagonal of the assembled operator at
!  the active PSI nodes.
!
        CALL ssa_jacobi_diag (ng,                                       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        iCstr, iCend, jCstr, jCend,               &
     &                        IstrP, IendP, JstrP, JendP,               &
     &                        pm, pn, pmask_ssa,                        &
     &                        nuH_rho, beta_rho,                        &
     &                        DPhi, phi_tab,                            &
     &                        diag_u, diag_v)
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_p2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, diag_u)
          CALL exchange_p2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, diag_v)
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, model, 2,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      diag_u, diag_v)
# endif
!
!  Initial residual r = b - A*u for the current velocity estimate.
!
        DO j=LBj,UBj
          DO i=LBi,UBi
            u_scr(i,j)=Si(i,j,ksa_new,isIuis)
            v_scr(i,j)=Si(i,j,ksa_new,isIvis)
          END DO
        END DO
!
        CALL ssa_matvec (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   iCstr, iCend, jCstr, jCend,                    &
     &                   IstrP, IendP, JstrP, JendP,                    &
     &                   pm, pn, pmask_ssa,                             &
     &                   nuH_rho, beta_rho,                             &
     &                   DPhi, phi_tab,                                 &
     &                   u_scr, v_scr, Au, Av)
!
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            IF (pmask_ssa(i,j).eq.1.0_r8) THEN
              Ru(i,j)=taudx_psi(i,j)-Au(i,j)
              Rv(i,j)=taudy_psi(i,j)-Av(i,j)
            ELSE
              Ru(i,j)=0.0_r8
              Rv(i,j)=0.0_r8
            END IF
          END DO
        END DO
!
!  Zero the residual on Dirichlet-constrained components so CG
!  never updates those nodes.
!
        CALL ssa_lbc_zero (ng,                                          &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Istr, Jstr, Jend,                            &
     &                     IstrP, IendP, JstrP, JendP,                  &
     &                     Lu_west_fix, Lv_south_fix,                   &
     &                     Lv_north_fix,                                &
     &                     Ru, Rv)
!
!  Preconditioned residual and initial search direction.
!
        dot_zr=0.0_r8
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            IF (pmask_ssa(i,j).eq.1.0_r8) THEN
              IF (diag_u(i,j).gt.eps) THEN
                Zu(i,j)=Ru(i,j)/diag_u(i,j)
              ELSE
                Zu(i,j)=0.0_r8
              END IF
              IF (diag_v(i,j).gt.eps) THEN
                Zv(i,j)=Rv(i,j)/diag_v(i,j)
              ELSE
                Zv(i,j)=0.0_r8
              END IF
              Du(i,j)=Zu(i,j)
              Dv(i,j)=Zv(i,j)
              dot_zr=dot_zr+Zu(i,j)*Ru(i,j)+Zv(i,j)*Rv(i,j)
            ELSE
              Zu(i,j)=0.0_r8
              Zv(i,j)=0.0_r8
              Du(i,j)=0.0_r8
              Dv(i,j)=0.0_r8
            END IF
          END DO
        END DO
!
        CALL ssa_lbc_zero (ng,                                          &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Istr, Jstr, Jend,                            &
     &                     IstrP, IendP, JstrP, JendP,                  &
     &                     Lu_west_fix, Lv_south_fix,                   &
     &                     Lv_north_fix,                                &
     &                     Du, Dv)
!
# ifdef DISTRIBUTE
        op_handle='SUM'
        CALL mp_reduce (ng, model, 1, dot_zr, op_handle)
# endif
        resid_0=SQRT(ABS(dot_zr))
        resid=resid_0
!
!-----------------------------------------------------------------------
!  Inner PCG loop.
!-----------------------------------------------------------------------
!
        DO cg_iter=1,cg_max_iter(ng)
          IF (resid.le.cgtol*resid_0 .or. resid.lt.eps) EXIT
!
!  Matrix action on the search direction, A*d.
!
          CALL ssa_matvec (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     iCstr, iCend, jCstr, jCend,                  &
     &                     IstrP, IendP, JstrP, JendP,                  &
     &                     pm, pn, pmask_ssa,                           &
     &                     nuH_rho, beta_rho,                           &
     &                     DPhi, phi_tab,                               &
     &                     Du, Dv, Au, Av)
!
!  alpha = (z,r) / (d,A*d);  both inner products are reduced in a
!  single pass over the active nodes.
!
          dot_dad=0.0_r8
          dot_zr_old=0.0_r8
          DO j=JstrP,JendP
            DO i=IstrP,IendP
              IF (pmask_ssa(i,j).eq.1.0_r8) THEN
                dot_dad=dot_dad+Du(i,j)*Au(i,j)+Dv(i,j)*Av(i,j)
                dot_zr_old=dot_zr_old+Zu(i,j)*Ru(i,j)+Zv(i,j)*Rv(i,j)
              END IF
            END DO
          END DO
!
# ifdef DISTRIBUTE
          op_handle='SUM'
          CALL mp_reduce (ng, model, 1, dot_dad, op_handle)
          CALL mp_reduce (ng, model, 1, dot_zr_old, op_handle)
# endif
!
          IF (ABS(dot_dad).lt.eps) EXIT
          alpha=dot_zr_old/dot_dad
!
          IF (alpha.ne.alpha.or.alpha.gt.1.0e30_r8) THEN
            IF (Master) WRITE(stdout,'(a,i3,3e12.4)')                   &
     &          ' SSA CG NaN: iter=',cg_iter,' dot_dad=',dot_dad,       &
     &          ' dot_zr_old=',dot_zr_old
            EXIT
          END IF
!
!  Update the solution and the residual.
!
          DO j=JstrP,JendP
            DO i=IstrP,IendP
              IF (pmask_ssa(i,j).eq.1.0_r8) THEN
                Si(i,j,ksa_new,isIuis)=Si(i,j,ksa_new,isIuis)+          &
     &                                alpha*Du(i,j)
                Si(i,j,ksa_new,isIvis)=Si(i,j,ksa_new,isIvis)+          &
     &                                alpha*Dv(i,j)
                Ru(i,j)=Ru(i,j)-alpha*Au(i,j)
                Rv(i,j)=Rv(i,j)-alpha*Av(i,j)
              END IF
            END DO
          END DO
!
          CALL ssa_lbc_zero (ng,                                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Istr, Jstr, Jend,                          &
     &                       IstrP, IendP, JstrP, JendP,                &
     &                       Lu_west_fix, Lv_south_fix,                 &
     &                       Lv_north_fix,                              &
     &                       Ru, Rv)
!
!  Preconditioned residual update and (z,r) accumulation
!  (Fletcher-Reeves form:  beta = (z_k+1,r_k+1)/(z_k,r_k)).
!
          dot_zr=0.0_r8
          DO j=JstrP,JendP
            DO i=IstrP,IendP
              IF (pmask_ssa(i,j).eq.1.0_r8) THEN
                IF (diag_u(i,j).gt.eps) THEN
                  Zu(i,j)=Ru(i,j)/diag_u(i,j)
                ELSE
                  Zu(i,j)=0.0_r8
                END IF
                IF (diag_v(i,j).gt.eps) THEN
                  Zv(i,j)=Rv(i,j)/diag_v(i,j)
                ELSE
                  Zv(i,j)=0.0_r8
                END IF
                dot_zr=dot_zr+Zu(i,j)*Ru(i,j)+Zv(i,j)*Rv(i,j)
              END IF
            END DO
          END DO
!
# ifdef DISTRIBUTE
          op_handle='SUM'
          CALL mp_reduce (ng, model, 1, dot_zr, op_handle)
# endif
!
          IF (ABS(dot_zr_old).lt.eps) EXIT
          beta_cg=dot_zr/dot_zr_old
!
!  New search direction.
!
          DO j=JstrP,JendP
            DO i=IstrP,IendP
              IF (pmask_ssa(i,j).eq.1.0_r8) THEN
                Du(i,j)=Zu(i,j)+beta_cg*Du(i,j)
                Dv(i,j)=Zv(i,j)+beta_cg*Dv(i,j)
              END IF
            END DO
          END DO
!
          CALL ssa_lbc_zero (ng,                                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Istr, Jstr, Jend,                          &
     &                       IstrP, IendP, JstrP, JendP,                &
     &                       Lu_west_fix, Lv_south_fix,                 &
     &                       Lv_north_fix,                              &
     &                       Du, Dv)
!
          resid=SQRT(ABS(dot_zr))
!
        END DO  ! cg_iter
!
!-----------------------------------------------------------------------
!  Post-CG robustness measures:  halo refresh + Dirichlet BCs,
!  speed cap, Picard relaxation.
!-----------------------------------------------------------------------
!
        CALL ssa_exchange_uv (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IstrP, IendP, JstrP, JendP,               &
     &                        ksa_new, Si)
!
        CALL ssa_lbc_velocity (ng, ksa_new,                             &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Istr, Jstr, Jend,                        &
     &                         IstrP, IendP, JstrP, JendP,              &
     &                         Lu_west_fix, Lu_west_clamp,              &
     &                         Lv_south_fix, Lv_south_clamp,            &
     &                         Lv_north_fix, Lv_north_clamp,            &
     &                         Si)
!
!  PISM-style velocity cap.
!
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            IF (pmask_ssa(i,j).eq.1.0_r8) THEN
              u_slide=SQRT(Si(i,j,ksa_new,isIuis)**2+                   &
     &                     Si(i,j,ksa_new,isIvis)**2)
              IF (u_slide.gt.ssa_max_speed(ng)) THEN
                cff=ssa_max_speed(ng)/u_slide
                Si(i,j,ksa_new,isIuis)=Si(i,j,ksa_new,isIuis)*cff
                Si(i,j,ksa_new,isIvis)=Si(i,j,ksa_new,isIvis)*cff
              END IF
            END IF
          END DO
        END DO
!
!  Picard relaxation u = alpha*u_new + (1-alpha)*u_old.  Prevents
!  oscillation between high/low viscosity regimes.  Skipped on the
!  first iteration, where ksa_old holds the (physical) initial
!  guess rather than a Picard iterate.
!
        IF (nl_iter.gt.1) THEN
          DO j=JstrP,JendP
            DO i=IstrP,IendP
              IF (pmask_ssa(i,j).eq.1.0_r8) THEN
                Si(i,j,ksa_new,isIuis)=picard_alpha(ng)*                &
     &            Si(i,j,ksa_new,isIuis)+                               &
     &            (1.0_r8-picard_alpha(ng))*                            &
     &            Si(i,j,ksa_old,isIuis)
                Si(i,j,ksa_new,isIvis)=picard_alpha(ng)*                &
     &            Si(i,j,ksa_new,isIvis)+                               &
     &            (1.0_r8-picard_alpha(ng))*                            &
     &            Si(i,j,ksa_old,isIvis)
              END IF
            END DO
          END DO
        END IF
!
!  Refresh the velocity halo and the Dirichlet BCs after the cap
!  and relaxation, then recompute the viscosity so nuH is
!  consistent with the relaxed velocity used in the next Picard
!  iteration and in the convergence check below.
!
        CALL ssa_exchange_uv (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IstrP, IendP, JstrP, JendP,               &
     &                        ksa_new, Si)
!
        CALL ssa_lbc_velocity (ng, ksa_new,                             &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Istr, Jstr, Jend,                        &
     &                         IstrP, IendP, JstrP, JendP,              &
     &                         Lu_west_fix, Lu_west_clamp,              &
     &                         Lv_south_fix, Lv_south_clamp,            &
     &                         Lv_north_fix, Lv_north_clamp,            &
     &                         Si)
!
        CALL ssa_viscosity (ng,                                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Istr, Iend, Jstr, Jend,                     &
     &                      ksa_old, ksa_new,                           &
     &                      pm, pn,                                     &
     &                      Si, Fi,                                     &
     &                      nuH_rho)
!
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Fi(i,j,icInuh)=nuH_rho(i,j)
          END DO
        END DO
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_p2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, nuH_rho)
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      nuH_rho)
# endif
!
!  True nonlinear residual ||A(u)u - b|| with the relaxed velocity
!  and the updated viscosity (one extra matrix action per Picard
!  iteration;  the legacy solver reused the last CG search
!  direction action, which does not measure the SSA residual).
!
        DO j=LBj,UBj
          DO i=LBi,UBi
            u_scr(i,j)=Si(i,j,ksa_new,isIuis)
            v_scr(i,j)=Si(i,j,ksa_new,isIvis)
          END DO
        END DO
!
        CALL ssa_matvec (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   iCstr, iCend, jCstr, jCend,                    &
     &                   IstrP, IendP, JstrP, JendP,                    &
     &                   pm, pn, pmask_ssa,                             &
     &                   nuH_rho, beta_rho,                             &
     &                   DPhi, phi_tab,                                 &
     &                   u_scr, v_scr, Au, Av)
!
!-----------------------------------------------------------------------
!  Picard convergence check:  dual criterion.
!    1. Relative velocity change
!           ||u_relaxed - u_old|| / ||u_relaxed||  <  nl_tol
!    2. Relative residual norm
!           ||A(u)u - b|| / ||b||                     <  nl_tol
!  Either criterion passing triggers convergence.  The residual
!  criterion catches stagnation at a non-solution fixed point.
!-----------------------------------------------------------------------
!
        vel_change=0.0_r8
        vel_norm=0.0_r8
        resid_norm=0.0_r8
        rhs_norm=0.0_r8
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            IF (pmask_ssa(i,j).eq.1.0_r8) THEN
              vel_change=vel_change+                                    &
     &          (Si(i,j,ksa_new,isIuis)-Si(i,j,ksa_old,isIuis))**2+     &
     &          (Si(i,j,ksa_new,isIvis)-Si(i,j,ksa_old,isIvis))**2
              vel_norm=vel_norm+                                        &
     &          Si(i,j,ksa_new,isIuis)**2+                              &
     &          Si(i,j,ksa_new,isIvis)**2
              resid_norm=resid_norm+                                    &
     &          (Au(i,j)-taudx_psi(i,j))**2+                            &
     &          (Av(i,j)-taudy_psi(i,j))**2
              rhs_norm=rhs_norm+                                        &
     &          taudx_psi(i,j)**2+taudy_psi(i,j)**2
            END IF
          END DO
        END DO
# ifdef DISTRIBUTE
        op_handle='SUM'
        CALL mp_reduce (ng, model, 1, vel_change, op_handle)
        CALL mp_reduce (ng, model, 1, vel_norm, op_handle)
        CALL mp_reduce (ng, model, 1, resid_norm, op_handle)
        CALL mp_reduce (ng, model, 1, rhs_norm, op_handle)
# endif
!
        IF (vel_norm.gt.eps) THEN
          nl_resid=SQRT(vel_change/vel_norm)
        ELSE
          nl_resid=0.0_r8
        END IF
!
        IF (nl_resid.le.nl_tol(ng)) THEN
          converged=.TRUE.
        END IF
!
        IF (.not.converged .and. rhs_norm.gt.eps) THEN
          IF (SQRT(resid_norm/rhs_norm).le.nl_tol(ng)) THEN
            converged=.TRUE.
          END IF
        END IF
!
!  DEBUG:  Picard iteration diagnostics.  Entirely gated by
!  ISF_DEBUG so production runs pay no output cost.
!
# ifdef ISF_DEBUG
        IF (Master .and. tile.eq.0) THEN
          IF (converged.or.nl_iter.eq.1.or.                             &
     &        MOD(nl_iter,10).eq.0) THEN
            u_slide=0.0_r8
            DO j=JstrP,JendP
              DO i=IstrP,IendP
                IF (pmask_ssa(i,j).eq.1.0_r8) THEN
                  cff1=SQRT(Si(i,j,ksa_new,isIuis)**2+                  &
     &                     Si(i,j,ksa_new,isIvis)**2)
                  IF (cff1.gt.u_slide) u_slide=cff1
                END IF
              END DO
            END DO
            WRITE(stdout,'(A,I3,A,ES12.4,A,ES12.4,A,L1)')               &
     &          '  Picard iter ',nl_iter,' nl_resid=',nl_resid,         &
     &          ' max|u|=',u_slide,' conv=',converged
          END IF
        END IF
# endif
!
        IF (converged) EXIT
!
!  Use the current velocity as the initial guess of the next Picard
!  iteration.
!
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            Si(i,j,ksa_old,isIuis)=Si(i,j,ksa_new,isIuis)
            Si(i,j,ksa_old,isIvis)=Si(i,j,ksa_new,isIvis)
          END DO
        END DO
!
      END DO  ! nl_iter (Picard)
!
!=======================================================================
!  Stage 10:  NaN guard and final halo refresh.
!=======================================================================
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          IF (pmask_ssa(i,j).eq.1.0_r8) THEN
            IF (Si(i,j,ksa_new,isIuis).ne.Si(i,j,ksa_new,isIuis).or.    &
     &          Si(i,j,ksa_new,isIvis).ne.Si(i,j,ksa_new,isIvis)) THEN
              IF (Master) WRITE(stdout,'(a,2i4)')                       &
     &            ' SSA NaN at PSI node:',i,j
              Si(i,j,ksa_new,isIuis)=0.0_r8
              Si(i,j,ksa_new,isIvis)=0.0_r8
            END IF
          END IF
        END DO
      END DO
!
!  Final velocity halo exchange so downstream routines (thickness
!  advection, output) see consistent tile-boundary data.
!
      CALL ssa_exchange_uv (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IstrP, IendP, JstrP, JendP,                 &
     &                      ksa_new, Si)
!
      CALL ssa_lbc_velocity (ng, ksa_new,                               &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Istr, Jstr, Jend,                          &
     &                       IstrP, IendP, JstrP, JendP,                &
     &                       Lu_west_fix, Lu_west_clamp,                &
     &                       Lv_south_fix, Lv_south_clamp,              &
     &                       Lv_north_fix, Lv_north_clamp,              &
     &                       Si)
!
      RETURN
      END SUBROUTINE isf_ssa_tile
!
!***********************************************************************
      SUBROUTINE ssa_exchange_uv (ng, tile, model,                      &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IstrP, IendP, JstrP, JendP,            &
     &                           ksa, Si)
!***********************************************************************
!
!  Refresh the halo of the two PSI-point velocity components of
!  the state array Si(:,:,ksa,isIuis/isIvis).  The components are
!  staged through internal scratch arrays so the exchanges operate
!  on contiguous arrays, and the FULL range (including halo) is
!  copied back afterwards:  downstream kernels (ssa_viscosity,
!  thickness advection) read the Si velocity halo cells at tile
!  seams, so the state array halo must carry the exchanged values.
!
      USE mod_kinds
      USE mod_param, ONLY : NghostPoints
      USE mod_scalars
      USE mod_parallel
      USE exchange_2d_mod, ONLY : exchange_p2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
      USE mod_isf, ONLY : nIsfS, isIuis, isIvis
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IstrP, IendP, JstrP, JendP
      integer, intent(in) :: ksa
!
      real(r8), intent(inout) :: Si(LBi:UBi,LBj:UBj,2,nIsfS)
!
!  Local variable declarations.
!
      integer :: i, j
!
      real(r8), dimension(LBi:UBi,LBj:UBj) :: u_scr, v_scr
!
!-----------------------------------------------------------------------
!  Stage, exchange and copy back the velocity.  The copy-back covers
!  the full declared range so the halo of Si is refreshed with the
!  exchanged values (mp_exchange2d fills the scratch halo from the
!  neighbouring tiles' interiors).
!-----------------------------------------------------------------------
!
      DO j=LBj,UBj
        DO i=LBi,UBi
          u_scr(i,j)=Si(i,j,ksa,isIuis)
          v_scr(i,j)=Si(i,j,ksa,isIvis)
        END DO
      END DO
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, u_scr)
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, v_scr)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    u_scr, v_scr)
# endif
!
      DO j=LBj,UBj
        DO i=LBi,UBi
          Si(i,j,ksa,isIuis)=u_scr(i,j)
          Si(i,j,ksa,isIvis)=v_scr(i,j)
        END DO
      END DO
!
      RETURN
      END SUBROUTINE ssa_exchange_uv
!
!***********************************************************************
      SUBROUTINE ssa_lbc_velocity (ng, ksa,                             &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             Istr, Jstr, Jend,                    &
     &                             IstrP, IendP, JstrP, JendP,          &
     &                             Lu_west_fix, Lu_west_clamp,          &
     &                             Lv_south_fix, Lv_south_clamp,        &
     &                             Lv_north_fix, Lv_north_clamp,        &
     &                             Si)
!***********************************************************************
!
!  Enforce the Dirichlet lateral boundary conditions on the SSA
!  velocity state:
!
!    closed   -> u = 0      (no-slip / ice divide)
!    clamped  -> u = boundary data from ISF_LOBC
!    gradient -> free (no enforcement)
!
!  Dirichlet enforcement acts on the normal component of each edge:
!  u on the western edge, v on the southern/northern edges.  The
!  eastern (calving front) edge is always left free.
!
      USE mod_kinds
      USE mod_param, ONLY : Mm
      USE mod_isf
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, ksa
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Istr, Jstr, Jend
      integer, intent(in) :: IstrP, IendP, JstrP, JendP
      logical, intent(in) :: Lu_west_fix, Lu_west_clamp
      logical, intent(in) :: Lv_south_fix, Lv_south_clamp
      logical, intent(in) :: Lv_north_fix, Lv_north_clamp
!
      real(r8), intent(inout) :: Si(LBi:UBi,LBj:UBj,2,nIsfS)
!
!  Local variable declarations.
!
      integer :: i, j
!
!-----------------------------------------------------------------------
!  Apply Dirichlet lateral boundary conditions.
!-----------------------------------------------------------------------
!
      IF (Jstr.eq.1.and.Lv_south_fix) THEN
        DO i=IstrP,IendP
          IF (Lv_south_clamp) THEN
            Si(i,JstrP,ksa,isIvis)=ISF_LOBC(isIvis,ng)%isf_south(i)
          ELSE
            Si(i,JstrP,ksa,isIvis)=0.0_r8
          END IF
        END DO
      END IF
!
      IF (Jend.eq.Mm(ng).and.Lv_north_fix) THEN
        DO i=IstrP,IendP
          IF (Lv_north_clamp) THEN
            Si(i,JendP,ksa,isIvis)=ISF_LOBC(isIvis,ng)%isf_north(i)
          ELSE
            Si(i,JendP,ksa,isIvis)=0.0_r8
          END IF
        END DO
      END IF
!
      IF (Istr.eq.1.and.Lu_west_fix) THEN
        DO j=JstrP,JendP
          IF (Lu_west_clamp) THEN
            Si(IstrP,j,ksa,isIuis)=ISF_LOBC(isIuis,ng)%isf_west(j)
          ELSE
            Si(IstrP,j,ksa,isIuis)=0.0_r8
          END IF
        END DO
      END IF
!
      RETURN
      END SUBROUTINE ssa_lbc_velocity
!
!***********************************************************************
      SUBROUTINE ssa_lbc_zero (ng,                                      &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Istr, Jstr, Jend,                        &
     &                         IstrP, IendP, JstrP, JendP,              &
     &                         Lu_west_fix, Lv_south_fix,               &
     &                         Lv_north_fix,                            &
     &                         uf, vf)
!***********************************************************************
!
!  Zero the Dirichlet-constrained components of a general PSI-point
!  field pair.  Used to keep the CG residual and search direction
!  zero on fixed nodes so the iteration never updates them.  The
!  initial guess already satisfies the boundary values (enforced by
!  ssa_lbc_velocity), hence zero residual means the solution keeps
!  them unchanged.
!
      USE mod_kinds
      USE mod_param, ONLY : Mm
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Istr, Jstr, Jend
      integer, intent(in) :: IstrP, IendP, JstrP, JendP
      logical, intent(in) :: Lu_west_fix
      logical, intent(in) :: Lv_south_fix
      logical, intent(in) :: Lv_north_fix
!
      real(r8), intent(inout) :: uf(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vf(LBi:UBi,LBj:UBj)
!
!  Local variable declarations.
!
      integer :: i, j
!
!-----------------------------------------------------------------------
!  Zero fixed components.
!-----------------------------------------------------------------------
!
      IF (Jstr.eq.1.and.Lv_south_fix) THEN
        DO i=IstrP,IendP
          vf(i,JstrP)=0.0_r8
        END DO
      END IF
!
      IF (Jend.eq.Mm(ng).and.Lv_north_fix) THEN
        DO i=IstrP,IendP
          vf(i,JendP)=0.0_r8
        END DO
      END IF
!
      IF (Istr.eq.1.and.Lu_west_fix) THEN
        DO j=JstrP,JendP
          uf(IstrP,j)=0.0_r8
        END DO
      END IF
!
      RETURN
      END SUBROUTINE ssa_lbc_zero
!
!***********************************************************************
      SUBROUTINE ssa_driving_stress (ng, tile, model,                   &
     &                               LBi, UBi, LBj, UBj,                &
     &                               Istr, Iend, Jstr, Jend,            &
     &                               IstrP, IendP, JstrP, JendP,        &
     &                               ksa_old,                           &
     &                               pm, pn, pmask_ssa,                 &
     &                               Si, Fi,                            &
     &                               taudx, taudy)
!***********************************************************************
!
!  Compute the driving stress RHS at PSI points:
!
!    b_k = -rho_i * g * H_psi * grad(S) * (nodal area)
!
!  This is the mass-lumped Q1 nodal force vector:  a finite
!  difference gradient of the smoothed surface elevation on the 4
!  RHO cells surrounding each PSI node, multiplied by the lumped
!  nodal area (jac/4;  4 elements share each Q1 node).
!
!  The surface elevation S (Fi(:,:,icIzsf)) is smoothed with
!  n_smooth iterations of a 5-point Laplacian filter before the
!  gradient is taken.  This regularizes the grounding-line kink
!  where S transitions between the grounded (S = -h + H) and
!  floating (S = (1-rho_i/rho_w) H) branches.
!
!  Ice margin taper:  the driving stress is linearly tapered to
!  zero within n_taper cells of the ice margin (any of the 8
!  neighbouring RHO cells ice-free).  The missing resistance is to
!  be provided by an ice-front pressure boundary term (future
!  work).  Margin PSI nodes with any ice-free surrounding cell are
!  skipped entirely.
!
!  Reference:  PISM compute_nuH_everywhere (slope regularization),
!  MITgcm STREAMICE streamice_driving_stress.F.
!
      USE mod_kinds
      USE mod_param, ONLY : NghostPoints
      USE mod_scalars
      USE mod_parallel
      USE exchange_2d_mod, ONLY : exchange_p2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
      USE mod_isf
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Istr, Iend, Jstr, Jend
      integer, intent(in) :: IstrP, IendP, JstrP, JendP
      integer, intent(in) :: ksa_old
!
      real(r8), intent(in)    :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: pmask_ssa(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: Si(LBi:UBi,LBj:UBj,2,nIsfS)
      real(r8), intent(in)    :: Fi(LBi:UBi,LBj:UBj,nIsfF)
      real(r8), intent(inout) :: taudx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: taudy(LBi:UBi,LBj:UBj)
!
!  Local variable declarations.
!
      integer :: i, j, di, dj, is_iter, dist_margin
!
      real(r8) :: H11, H21, H12, H22, H_psi
      real(r8) :: dSdx, dSdy, dx_rho, dy_rho, jac_avg
      real(r8) :: zs_sum, taper_w, rhoi_g
      real(r8), parameter :: gravity = 9.81_r8
!
      real(r8), dimension(LBi:UBi,LBj:UBj) :: zs_smooth
!
!-----------------------------------------------------------------------
!  Initialize the nodal forces.
!-----------------------------------------------------------------------
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          taudx(i,j)=0.0_r8
          taudy(i,j)=0.0_r8
        END DO
      END DO
!
      rhoi_g=rho_i*gravity
!
!-----------------------------------------------------------------------
!  Smooth the surface elevation for the driving stress gradient.
!-----------------------------------------------------------------------
!
!  Initialize the smoothed field (including the halo for
!  tile-boundary elements).  Halo values are refreshed from the
!  neighbouring tile interiors by the exchange at the top of each
!  smoothing iteration.
!
      DO j=LBj,UBj
        DO i=LBi,UBi
          IF (Si(i,j,ksa_old,isIhis).gt.0.0_r8) THEN
            zs_smooth(i,j)=Fi(i,j,icIzsf)
          ELSE
            zs_smooth(i,j)=0.0_r8
          END IF
        END DO
      END DO
!
!  Apply the 5-point Laplacian smoothing (n_smooth iterations).
!  Only ice-covered cells are smoothed;  ice-free neighbours are
!  replaced by the centre value (no-flux treatment at the margin).
!  The interior-only loop keeps the writes inside the array bounds.
!
      DO is_iter=1,n_smooth(ng)
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      zs_smooth)
# endif
        DO j=Jstr,Jend
          DO i=Istr,Iend
            IF (Si(i,j,ksa_old,isIhis).gt.0.0_r8) THEN
              zs_sum=zs_smooth(i,j)
              IF (Si(i+1,j,ksa_old,isIhis).gt.0.0_r8) THEN
                zs_sum=zs_sum+zs_smooth(i+1,j)
              ELSE
                zs_sum=zs_sum+zs_smooth(i,j)
              END IF
              IF (Si(i-1,j,ksa_old,isIhis).gt.0.0_r8) THEN
                zs_sum=zs_sum+zs_smooth(i-1,j)
              ELSE
                zs_sum=zs_sum+zs_smooth(i,j)
              END IF
              IF (Si(i,j+1,ksa_old,isIhis).gt.0.0_r8) THEN
                zs_sum=zs_sum+zs_smooth(i,j+1)
              ELSE
                zs_sum=zs_sum+zs_smooth(i,j)
              END IF
              IF (Si(i,j-1,ksa_old,isIhis).gt.0.0_r8) THEN
                zs_sum=zs_sum+zs_smooth(i,j-1)
              ELSE
                zs_sum=zs_sum+zs_smooth(i,j)
              END IF
              zs_smooth(i,j)=zs_sum/5.0_r8
            END IF
          END DO
        END DO
      END DO
!
!  Exchange the halo of the smoothed field (needed by the
!  tile-boundary PSI stencils below).
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, zs_smooth)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    zs_smooth)
# endif
!
!-----------------------------------------------------------------------
!  Assemble the driving stress at the active PSI nodes.
!-----------------------------------------------------------------------
!
!  At PSI point (i,j) the 4 surrounding RHO cells are
!  (i-1,j-1), (i,j-1), (i-1,j), (i,j):
!
!    dS/dx = [(zs(i,j-1)+zs(i,j)) - (zs(i-1,j-1)+zs(i-1,j))] / (2 dx)
!    dS/dy = [(zs(i-1,j)+zs(i,j)) - (zs(i-1,j-1)+zs(i,j-1))] / (2 dy)
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          IF (pmask_ssa(i,j).eq.1.0_r8) THEN
!
            H11=Si(i-1,j-1,ksa_old,isIhis)
            H21=Si(i,  j-1,ksa_old,isIhis)
            H12=Si(i-1,j,  ksa_old,isIhis)
            H22=Si(i,  j,  ksa_old,isIhis)
!
!  Skip margin nodes with any ice-free surrounding cell.
!
            IF (H11.gt.0.0_r8.and.H21.gt.0.0_r8.and.                    &
     &          H12.gt.0.0_r8.and.H22.gt.0.0_r8) THEN
!
!  Ice margin taper:  distance to the nearest ice-free cell among
!  the 8 neighbours, normalized by n_taper and clamped to 1.
!
              dist_margin=n_taper(ng)
              DO dj=-1,1
                DO di=-1,1
                  IF (di.ne.0.or.dj.ne.0) THEN
                    IF (Si(i+di,j+dj,ksa_old,isIhis).le.0.0_r8) THEN
                      dist_margin=MIN(dist_margin,                      &
     &                    MAX(ABS(di),ABS(dj)))
                    END IF
                  END IF
                END DO
              END DO
              taper_w=MIN(REAL(dist_margin,r8)/                         &
     &                    REAL(n_taper(ng),r8),1.0_r8)
!
!  FD gradient of the smoothed surface at the PSI node.
!
              dx_rho=0.5_r8*(1.0_r8/pm(i-1,j)+1.0_r8/pm(i,j))
              dy_rho=0.5_r8*(1.0_r8/pn(i,j-1)+1.0_r8/pn(i,j))
              jac_avg=0.25_r8*(1.0_r8/(pm(i-1,j-1)*pn(i-1,j-1))+        &
     &                         1.0_r8/(pm(i,  j-1)*pn(i,  j-1))+        &
     &                         1.0_r8/(pm(i-1,j  )*pn(i-1,j  ))+        &
     &                         1.0_r8/(pm(i,  j  )*pn(i,  j  )))
!
              dSdx=((zs_smooth(i,j-1)+zs_smooth(i,j))-                  &
     &              (zs_smooth(i-1,j-1)+zs_smooth(i-1,j)))/             &
     &             (2.0_r8*dx_rho)
              dSdy=((zs_smooth(i-1,j)+zs_smooth(i,j))-                  &
     &              (zs_smooth(i-1,j-1)+zs_smooth(i,j-1)))/             &
     &             (2.0_r8*dy_rho)
!
              H_psi=0.25_r8*(H11+H21+H12+H22)
!
!  FE-consistent nodal force (lumped nodal area = jac_avg/4).
!
              taudx(i,j)=-taper_w*rhoi_g*H_psi*dSdx*                    &
     &                       jac_avg*0.25_r8
              taudy(i,j)=-taper_w*rhoi_g*H_psi*dSdy*                    &
     &                       jac_avg*0.25_r8
            END IF
          END IF
        END DO
      END DO
!
!  Exchange the halo of the assembled RHS.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, taudx)
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, taudy)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    taudx, taudy)
# endif
!
      RETURN
      END SUBROUTINE ssa_driving_stress
!
!***********************************************************************
      SUBROUTINE ssa_viscosity (ng,                                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Istr, Iend, Jstr, Jend,                 &
     &                          ksa_old, ksa_new,                       &
     &                          pm, pn,                                 &
     &                          Si, Fi,                                 &
     &                          nuH)
!***********************************************************************
!
!  Compute the depth-integrated effective viscosity at RHO points:
!
!    gamma = eps_xx^2 + eps_yy^2 + eps_xx*eps_yy + eps_xy^2
!           + eps_glen_min^2            (squared second invariant)
!    nuH   = 0.5 * B * gamma^((1-n)/(2n)) * H
!
!  where B = A^(-1/n) is the Glen hardness.  The exponent (1-n)/(2n)
!  acts on the SQUARED invariant gamma:  this is the PISM convention
!  (FlowLaw::effective_viscosity, FlowLaw.cc) and is equivalent to
!  the STREAMICE form  nuH = 0.5*B*eps_e^((1-n)/n)*H  with
!  eps_e = sqrt(gamma)  (streamice_visc_beta.F);  both reduce to the
!  MacAyeal (1989) Glen-law effective viscosity nu = 0.5*B*
!  gamma^((1-n)/(2n)).  The Glen rate factor
!  A is either the constant A_glen (isothermal override, used for
!  ISOMIP-style tests) or the Paterson (1994) temperature-dependent
!  formula with pressure-adjusted temperature
!  T* = T + 273.15 + 8.66e-4 * 0.5 * H  (Bueler & Brown 2009):
!
!    A = 3.61e-13 * exp(-60000/(R T*)),   T* <= 263.15 K
!    A = 1.73e3   * exp(-139000/(R T*)),  T* >  263.15 K
!
!  Strain rates are centered differences of the 4 PSI-corner
!  velocities of each RHO cell.  PISM robustness measures:
!    - strength extension for thin ice (H < min_thickness)
!    - additive nuH regularization (nuH += nuH_reg)
!
!  Reference:  Paterson (1994) Table 5.2, Bueler & Brown (2009),
!  MITgcm STREAMICE streamice_visc_beta.F, PISM compute_nuH_everywhere.
!
      USE mod_kinds
      USE mod_isf
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Istr, Iend, Jstr, Jend
      integer, intent(in) :: ksa_old, ksa_new
!
      real(r8), intent(in)    :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: Si(LBi:UBi,LBj:UBj,2,nIsfS)
      real(r8), intent(in)    :: Fi(LBi:UBi,LBj:UBj,nIsfF)
      real(r8), intent(inout) :: nuH(LBi:UBi,LBj:UBj)
!
!  Local variable declarations.
!
      integer :: i, j
!
      real(r8) :: H_rho, B_hardness, T_star, RT_star, A_pat
      real(r8) :: cff, eps_xx, eps_yy, eps_xy
!
!-----------------------------------------------------------------------
!  Compute viscosity at the RHO cells with ice (H > 0).  Cell (i,j)
!  has 4 PSI corner nodes:  (i,j), (i+1,j), (i,j+1), (i+1,j+1).
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          nuH(i,j)=0.0_r8
          H_rho=Si(i,j,ksa_old,isIhis)
          IF (H_rho.gt.0.0_r8) THEN
!
!  Glen rate factor:  constant (isothermal) or Paterson.
!
            IF (A_glen(ng).gt.0.0_r8) THEN
              B_hardness=A_glen(ng)**(-1.0_r8/n_glen(ng))
            ELSE
              T_star=Fi(i,j,icItis)+273.15_r8+                          &
     &               8.66E-4_r8*0.5_r8*H_rho
              RT_star=8.314_r8*T_star
              IF (T_star.le.263.15_r8) THEN
                A_pat=3.61E-13_r8*EXP(-6.0E4_r8/RT_star)
              ELSE
                A_pat=1.73E3_r8*EXP(-13.9E4_r8/RT_star)
              END IF
              B_hardness=A_pat**(-1.0_r8/n_glen(ng))
            END IF
!
!  Centered strain rates from the 4 PSI-corner velocities:
!  nodes 1=(i,j), 2=(i+1,j), 3=(i,j+1), 4=(i+1,j+1).
!
            eps_xx=0.5_r8*((Si(i+1,j,  ksa_new,isIuis)-                 &
     &                      Si(i,  j,  ksa_new,isIuis))+                &
     &                     (Si(i+1,j+1,ksa_new,isIuis)-                 &
     &                      Si(i,  j+1,ksa_new,isIuis)))*pm(i,j)
            eps_yy=0.5_r8*((Si(i,  j+1,ksa_new,isIvis)-                 &
     &                      Si(i,  j,  ksa_new,isIvis))+                &
     &                     (Si(i+1,j+1,ksa_new,isIvis)-                 &
     &                      Si(i+1,j,  ksa_new,isIvis)))*pn(i,j)
            eps_xy=0.5_r8*(                                             &
     &          0.5_r8*((Si(i+1,j,  ksa_new,isIvis)-                    &
     &                   Si(i,  j,  ksa_new,isIvis))+                   &
     &                  (Si(i+1,j+1,ksa_new,isIvis)-                    &
     &                   Si(i,  j+1,ksa_new,isIvis)))*pm(i,j)+          &
     &          0.5_r8*((Si(i,  j+1,ksa_new,isIuis)-                    &
     &                   Si(i,  j,  ksa_new,isIuis))+                   &
     &                  (Si(i+1,j+1,ksa_new,isIuis)-                    &
     &                   Si(i+1,j,  ksa_new,isIuis)))*pn(i,j))
!
!  Glen-law effective viscosity from the squared second invariant,
!  with strength extension for thin ice and additive regularization.
!  The bracket is strictly positive through the eps_glen_min^2 term.
!
            cff=0.5_r8*B_hardness*                                      &
     &          (eps_xx*eps_xx+eps_yy*eps_yy+eps_xx*eps_yy+             &
     &           eps_xy*eps_xy+                                         &
     &           eps_glen_min(ng)*eps_glen_min(ng))**                   &
     &          ((1.0_r8-n_glen(ng))/(2.0_r8*n_glen(ng)))
!
            IF (H_rho.lt.min_thickness(ng)) THEN
              nuH(i,j)=cff*min_thickness(ng)
            ELSE
              nuH(i,j)=cff*H_rho
            END IF
!
            nuH(i,j)=nuH(i,j)+nuH_reg(ng)
!
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE ssa_viscosity
!
!***********************************************************************
      SUBROUTINE ssa_jacobi_diag (ng,                                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            iCstr, iCend, jCstr, jCend,           &
     &                            IstrP, IendP, JstrP, JendP,           &
     &                            pm, pn, pmask_ssa,                    &
     &                            nuH, beta,                            &
     &                            DPhi, phi_tab,                        &
     &                            diag_u, diag_v)
!***********************************************************************
!
!  Assemble the Jacobi preconditioner:  the diagonal of the SSA
!  operator at the active PSI nodes.  For each Q1 node k the
!  diagonal is the element self-coupling with a unit velocity at
!  node k and zero elsewhere:
!
!    diag_u(k) = int [ nuH*(4*dphi_x^2 + dphi_y^2) + beta*phi^2 ] dA
!    diag_v(k) = int [ nuH*(4*dphi_y^2 + dphi_x^2) + beta*phi^2 ] dA
!
!  (set u_k = 1:  exx = dphi_x, eyy = 0, exy = 0.5*dphi_y, hence
!  u-row:  dphi_x*(4*dphi_x) + dphi_y*(2*0.5*dphi_y);  the v-row
!  follows by symmetry).
!
!  The cell loop covers the own cells plus the western/southern
!  halo ring (iCstr..iCend);  scatters land only on owned PSI
!  nodes (ownership filter).
!
!  Reference:  MITgcm STREAMICE CG_ADIAG.
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: iCstr, iCend, jCstr, jCend
      integer, intent(in) :: IstrP, IendP, JstrP, JendP
!
      real(r8), intent(in)    :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: pmask_ssa(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: nuH(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: beta(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: DPhi(4,4,2)
      real(r8), intent(in)    :: phi_tab(4,4)
      real(r8), intent(inout) :: diag_u(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: diag_v(LBi:UBi,LBj:UBj)
!
!  Local variable declarations.
!
      integer :: i, j, i1, j1, nq, k_node
!
      real(r8) :: dx_rho, dy_rho, jac
      real(r8) :: dphi_x, dphi_y, phv
!
!  Quadrature weights (2x2 Gauss on the reference [0,1]x[0,1]).
!
      real(r8), parameter :: wquad(4) = 0.25_r8
!
!-----------------------------------------------------------------------
!  Assemble the diagonal contributions.
!-----------------------------------------------------------------------
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          diag_u(i,j)=0.0_r8
          diag_v(i,j)=0.0_r8
        END DO
      END DO
!
      DO j=jCstr,jCend
        DO i=iCstr,iCend
          IF (nuH(i,j).gt.0.0_r8) THEN
            dx_rho=1.0_r8/pm(i,j)
            dy_rho=1.0_r8/pn(i,j)
            jac=1.0_r8/(pm(i,j)*pn(i,j))
!
            DO nq=1,4
              DO k_node=1,4
                IF (k_node.eq.1) THEN
                  i1=i;   j1=j
                ELSE IF (k_node.eq.2) THEN
                  i1=i+1; j1=j
                ELSE IF (k_node.eq.3) THEN
                  i1=i;   j1=j+1
                ELSE
                  i1=i+1; j1=j+1
                END IF
!
!  Ownership filter:  keep only contributions landing on PSI nodes
!  owned by this tile (seam nodes are fed by the halo ring cells
!  integrated here, and by the own cells integrated on the
!  neighbouring tile).
!
                IF (pmask_ssa(i1,j1).eq.1.0_r8.and.                     &
     &              i1.ge.IstrP.and.i1.le.IendP.and.                    &
     &              j1.ge.JstrP.and.j1.le.JendP) THEN
                  dphi_x=DPhi(k_node,nq,1)*2.0_r8/dx_rho
                  dphi_y=DPhi(k_node,nq,2)*2.0_r8/dy_rho
                  phv=phi_tab(k_node,nq)
!
                  diag_u(i1,j1)=diag_u(i1,j1)+                          &
     &              wquad(nq)*jac*nuH(i,j)*                             &
     &              (4.0_r8*dphi_x*dphi_x+dphi_y*dphi_y)+               &
     &              wquad(nq)*jac*beta(i,j)*phv*phv
!
                  diag_v(i1,j1)=diag_v(i1,j1)+                          &
     &              wquad(nq)*jac*nuH(i,j)*                             &
     &              (4.0_r8*dphi_y*dphi_y+dphi_x*dphi_x)+               &
     &              wquad(nq)*jac*beta(i,j)*phv*phv
                END IF
              END DO
            END DO
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE ssa_jacobi_diag
!
!***********************************************************************
      SUBROUTINE ssa_matvec (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       iCstr, iCend, jCstr, jCend,                &
     &                       IstrP, IendP, JstrP, JendP,                &
     &                       pm, pn, pmask_ssa,                         &
     &                       nuH, beta,                                 &
     &                       DPhi, phi_tab,                             &
     &                       xu, xv, Au, Av)
!***********************************************************************
!
!  Matrix-free action of the SSA operator on a PSI-point field
!  pair:  (Au,Av) = A * (xu,xv).  Q1 element assembly over the own
!  cells plus the western/southern halo ring, with 2x2 Gauss
!  quadrature.  At each quadrature point the velocity is
!  interpolated, strain rates are formed from the shape function
!  derivatives, and membrane stress divergence plus basal drag are
!  scattered to the 4 corner nodes:
!
!    (Au)_m += int [ nuH*(dphi_m/dx*(4 exx + 2 eyy)
!                         + dphi_m/dy*(2 exy))
!                    + beta*phi_m*u_q ] dA
!
!  with the analogous v-row.  The (4,2) coefficients follow from
!  the depth-integrated membrane stresses of the incompressible
!  SSA:  N_xx = 2 nuH (2 exx + eyy),  N_xy = 2 nuH exy.
!
!  Input/output halos are exchanged inside this routine:  the
!  assembly reads the PSI corners of every integrated cell (halo
!  included).
!
!  Reference:  MITgcm STREAMICE streamice_cg_functions.F (CG_ACTION).
!
      USE mod_kinds
      USE mod_param, ONLY : NghostPoints
      USE mod_scalars
      USE mod_parallel
      USE exchange_2d_mod, ONLY : exchange_p2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: iCstr, iCend, jCstr, jCend
      integer, intent(in) :: IstrP, IendP, JstrP, JendP
!
      real(r8), intent(in)    :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: pmask_ssa(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: nuH(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: beta(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: DPhi(4,4,2)
      real(r8), intent(in)    :: phi_tab(4,4)
      real(r8), intent(inout) :: xu(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: xv(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Au(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Av(LBi:UBi,LBj:UBj)
!
!  Local variable declarations.
!
      integer :: i, j, i1, j1, nq, k_node
!
      real(r8) :: dx_rho, dy_rho, jac
      real(r8) :: dphi_x, dphi_y, phv
      real(r8) :: uq, vq, exx, eyy, exy
!
!  Quadrature weights (2x2 Gauss on the reference [0,1]x[0,1]).
!
      real(r8), parameter :: wquad(4) = 0.25_r8
!
!-----------------------------------------------------------------------
!  Exchange the input halos:  the assembly reads all 4 PSI corners
!  of every integrated cell, including halo nodes owned by the
!  neighbouring tiles.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, xu)
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, xv)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    xu, xv)
# endif
!
!-----------------------------------------------------------------------
!  Element assembly.
!-----------------------------------------------------------------------
!
      DO j=JstrP-1,JendP+1
        DO i=IstrP-1,IendP+1
          Au(i,j)=0.0_r8
          Av(i,j)=0.0_r8
        END DO
      END DO
!
      DO j=jCstr,jCend
        DO i=iCstr,iCend
          IF (nuH(i,j).gt.0.0_r8) THEN
            dx_rho=1.0_r8/pm(i,j)
            dy_rho=1.0_r8/pn(i,j)
            jac=1.0_r8/(pm(i,j)*pn(i,j))
!
            DO nq=1,4
!
!  Interpolate the velocity and form the strain rates at the
!  quadrature point.
!
              uq=0.0_r8
              vq=0.0_r8
              exx=0.0_r8
              eyy=0.0_r8
              exy=0.0_r8
              DO k_node=1,4
                IF (k_node.eq.1) THEN
                  i1=i;   j1=j
                ELSE IF (k_node.eq.2) THEN
                  i1=i+1; j1=j
                ELSE IF (k_node.eq.3) THEN
                  i1=i;   j1=j+1
                ELSE
                  i1=i+1; j1=j+1
                END IF
                phv=phi_tab(k_node,nq)
                dphi_x=DPhi(k_node,nq,1)*2.0_r8/dx_rho
                dphi_y=DPhi(k_node,nq,2)*2.0_r8/dy_rho
                uq=uq+xu(i1,j1)*phv
                vq=vq+xv(i1,j1)*phv
                exx=exx+xu(i1,j1)*dphi_x
                eyy=eyy+xv(i1,j1)*dphi_y
                exy=exy+0.5_r8*(xu(i1,j1)*dphi_y+xv(i1,j1)*dphi_x)
              END DO
!
!  Scatter to the 4 corner nodes (ownership filter).
!
              DO k_node=1,4
                IF (k_node.eq.1) THEN
                  i1=i;   j1=j
                ELSE IF (k_node.eq.2) THEN
                  i1=i+1; j1=j
                ELSE IF (k_node.eq.3) THEN
                  i1=i;   j1=j+1
                ELSE
                  i1=i+1; j1=j+1
                END IF
                IF (pmask_ssa(i1,j1).eq.1.0_r8.and.                     &
     &              i1.ge.IstrP.and.i1.le.IendP.and.                    &
     &              j1.ge.JstrP.and.j1.le.JendP) THEN
                  dphi_x=DPhi(k_node,nq,1)*2.0_r8/dx_rho
                  dphi_y=DPhi(k_node,nq,2)*2.0_r8/dy_rho
                  phv=phi_tab(k_node,nq)
!
                  Au(i1,j1)=Au(i1,j1)+                                  &
     &              wquad(nq)*jac*nuH(i,j)*                             &
     &              (dphi_x*(4.0_r8*exx+2.0_r8*eyy)+                    &
     &               dphi_y*(2.0_r8*exy))+                              &
     &              wquad(nq)*jac*beta(i,j)*phv*uq
!
                  Av(i1,j1)=Av(i1,j1)+                                  &
     &              wquad(nq)*jac*nuH(i,j)*                             &
     &              (dphi_y*(4.0_r8*eyy+2.0_r8*exx)+                    &
     &               dphi_x*(2.0_r8*exy))+                              &
     &              wquad(nq)*jac*beta(i,j)*phv*vq
                END IF
              END DO
            END DO
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Exchange the output halos:  shared seam PSI node values computed
!  on the neighbouring tiles arrive through the halo refresh.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, Au)
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, Av)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Au, Av)
# endif
!
      RETURN
      END SUBROUTINE ssa_matvec
#endif
      END MODULE isf_solver_mod
