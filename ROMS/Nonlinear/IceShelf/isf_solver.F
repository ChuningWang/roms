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
!  Fi(:,:,ifIbta), computed by isf_stress from the ocean state         !
!  averaged over the SSA coupling window (Weertman or regularized      !
!  Coulomb sliding law).  Freezing beta at the coupling-window value   !
!  is the intended frozen-coefficient (Picard) linearization of the    !
!  drag term; only the viscosity is updated inside the Picard loop.    !
!                                                                      !
!  Input contract (caller must fill before calling isf_ssa):           !
!                                                                      !
!    Si(:,:,lsold,isIhis)          ice thickness [m], halo exchanged   !
!    Si(:,:,lsold,isIuis/isIvis)   velocity initial guess [m s-1]      !
!    Fi(:,:,ifItis)                ice temperature [degC]              !
!    Fi(:,:,ifIzsf)                ice surface elevation [m]           !
!    Fi(:,:,ifIbta)                basal drag coefficient (isf_stress) !
!    GRID%pmask_act                active solve mask at PSI points,    !
!                                  halo exchanged (isf_masks): OR of   !
!                                  the 4 surrounding RHO cells having  !
!                                  ice, so the ice-margin ring solves  !
!    GRID%rmask_ist                grounded/ice mask at RHO points,    !
!                                  halo exchanged (isf_masks)          !
!                                                                      !
!  Output:                                                             !
!                                                                      !
!    Si(:,:,lsnew,isIuis/isIvis)   solved SSA velocity [m s-1],        !
!                                  halo exchanged                      !
!    Fi(:,:,ifInuh)                depth-integrated viscosity          !
!    Fi(:,:,ifItdx/ifItdy)         driving stress diagnostics          !
!                                                                      !
!  Lateral boundaries are fully LBC-configurable per edge and per      !
!  velocity component (Cla / Clo / Gra, same flags as the ice state    !
!  LBC):  closed edges impose zero velocity, clamped edges fix the     !
!  velocity to ISF_LOBC data, gradient edges are free.  Dirichlet      !
!  enforcement (isf_pibc_tile for values, ssa_lbc_zero for the CG      !
!  residual/search directions) acts on every edge/component flagged    !
!  closed or clamped -- which edges carry which condition is an        !
!  application choice (e.g. MISOMIP uses a clamped western inflow      !
!  and a free eastern calving front, but other setups may differ).     !
!  Ice-front pressure boundary load:  front faces (ice cell against an !
!  ice-free neighbour, free gradient edges included) carry the net     !
!  hydrostatic push of the seawater as a Q1 line load (ssa_front_load),!
!  with the floating / grounded-in-water / land-margin cases of        !
!  Morland (1987) / MacAyeal (1989);  the ice-margin ring solves       !
!  under the OR active mask GRID%pmask_act.                            !
!                                                                      !
!  MPI notes:  PSI tile ranges do not overlap (get_bounds.F:           !
!  IstrP=my_Istr everywhere, IendP=my_Iend except on the eastern       !
!  domain edge), so a seam PSI node is owned by exactly one tile       !
!  while the RHO cells feeding its FE stencil lie on both sides of     !
!  the seam.  The element assembly therefore integrates one extra      !
!  ring of halo cells along the western/southern tile seams (bounds    !
!  IstrC/JstrC below; the legacy isf_cg_solve.h documents this         !
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
!    MacAyeal (1989), J. Geophys. Res. 94, 4071-4087.                  !
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
      USE mod_grid
      USE mod_iounits
      USE mod_isf
      USE mod_parallel
      USE mod_param
      USE mod_scalars
      USE mod_stepping
      USE isf_pibc_mod,    ONLY : isf_pibc_tile
!
      USE exchange_2d_mod, ONLY : exchange_p2d_tile, exchange_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
      USE distribute_mod,  ONLY : mp_reduce
# endif
!
      implicit none
!
      PUBLIC  :: isf_ssa
      PRIVATE
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE isf_ssa (ng, tile, model)
!***********************************************************************
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
     &                   GRID(ng) % h,                                  &
     &                   ISF(ng) % Ci(:,:,icIzta),                      &
     &                   GRID(ng) % rmask_ist,                          &
     &                   GRID(ng) % pmask_act,                          &
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
     &                         lsold, lsnew,                            &
     &                         pm, pn, h, zeta_avg,                     &
     &                         rmask_ist, pmask_act,                    &
     &                         Si, Fi)
!***********************************************************************
!
!  Tile-level SSA solver.  Orchestrates:  active-node mask, initial
!  guess (copy + basal-drag warm start), driving stress RHS, and the
!  Picard/PCG iteration.  All kernels (viscosity, matrix action, Jacobi
!  diagonal, driving stress, LBC enforcement) live in the private
!  subroutines below;  this routine owns halo exchanges, global
!  reductions and convergence control.
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: lsold, lsnew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)    :: pm(LBi:,LBj:)
      real(r8), intent(in)    :: pn(LBi:,LBj:)
      real(r8), intent(in)    :: h(LBi:,LBj:)
      real(r8), intent(in)    :: zeta_avg(LBi:,LBj:)
      real(r8), intent(in)    :: pmask_act(LBi:,LBj:)
      real(r8), intent(in)    :: rmask_ist(LBi:,LBj:)
      real(r8), intent(inout) :: Si(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: Fi(LBi:,LBj:,:)
# else
      real(r8), intent(in)    :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: zeta_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: pmask_act(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: rmask_ist(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Si(LBi:UBi,LBj:UBj,2,nIsfS)
      real(r8), intent(inout) :: Fi(LBi:UBi,LBj:UBj,nIsfF)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, cg_iter, nl_iter
!
      logical :: converged
!
      real(r8) :: cff, cff1, u_slide, beta_avg, jac_avg
!
!  Glen hardness field B = A^(-1/n) at RHO points, computed once
!  per solve by ssa_hardness (Stage 4):  it depends only on the ice
!  temperature and the lsold thickness, both frozen during the
!  Picard iteration, so the Paterson exponentials are not
!  re-evaluated by every ssa_viscosity call.
!
      real(r8), dimension(LBi:UBi,LBj:UBj) :: B_hard
!
!  Land-margin PSI indicator for the driving-stress taper
!  (ssa_land_margin):  1 at nodes adjacent to a land-margin front
!  cell;  the taper distance is measured to these nodes only
!  (plan option A -- marine fronts carry the full front load).
!
      real(r8), dimension(LBi:UBi,LBj:UBj) :: pmarg
!
!  Global-reduction buffer (batched mp_reduce).
!
      real(r8), dimension(4) :: norms
      character (len=3), dimension(4) :: op_norms
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
!  phi_tab holds the shape function VALUES at the [0,1] Gauss
!  abscissae gxa/gxb;  DPhi holds the shape function derivatives on
!  the reference element in a [-1,1]-convention scaling, so the
!  physical derivative is dphi/dx = DPhi*2/dx (and likewise in y),
!  as used by the assembly routines below.
!
      real(r8), parameter :: sqrt3_inv = 1.0_r8/SQRT(3.0_r8)
      real(r8), parameter :: gxa = 0.5_r8*(1.0_r8-sqrt3_inv)
      real(r8), parameter :: gxb = 0.5_r8*(1.0_r8+sqrt3_inv)
      real(r8), parameter :: phi_tab(4,4) = RESHAPE((/                  &
     &    (1.0_r8-gxa)*(1.0_r8-gxa), gxa*(1.0_r8-gxa),                  &
     &    (1.0_r8-gxa)*gxa,          gxa*gxa,                           &
     &    (1.0_r8-gxb)*(1.0_r8-gxa), gxb*(1.0_r8-gxa),                  &
     &    (1.0_r8-gxb)*gxa,          gxb*gxa,                           &
     &    (1.0_r8-gxa)*(1.0_r8-gxb), gxa*(1.0_r8-gxb),                  &
     &    (1.0_r8-gxa)*gxb,          gxa*gxb,                           &
     &    (1.0_r8-gxb)*(1.0_r8-gxb), gxb*(1.0_r8-gxb),                  &
     &    (1.0_r8-gxb)*gxb,          gxb*gxb                            &
     &  /),(/4,4/))
      real(r8), parameter :: DPhi(4,4,2) = RESHAPE((/                   &
     &   -0.25_r8*(1.0_r8-gxa)*2.0_r8,  0.25_r8*(1.0_r8-gxa)*2.0_r8,    &
     &   -0.25_r8*gxa*2.0_r8,           0.25_r8*gxa*2.0_r8,             &
     &   -0.25_r8*(1.0_r8-gxa)*2.0_r8,  0.25_r8*(1.0_r8-gxa)*2.0_r8,    &
     &   -0.25_r8*gxa*2.0_r8,           0.25_r8*gxa*2.0_r8,             &
     &   -0.25_r8*(1.0_r8-gxb)*2.0_r8,  0.25_r8*(1.0_r8-gxb)*2.0_r8,    &
     &   -0.25_r8*gxb*2.0_r8,           0.25_r8*gxb*2.0_r8,             &
     &   -0.25_r8*(1.0_r8-gxb)*2.0_r8,  0.25_r8*(1.0_r8-gxb)*2.0_r8,    &
     &   -0.25_r8*gxb*2.0_r8,           0.25_r8*gxb*2.0_r8,             &
     &   -0.25_r8*(1.0_r8-gxa)*2.0_r8, -0.25_r8*gxa*2.0_r8,             &
     &    0.25_r8*(1.0_r8-gxa)*2.0_r8,  0.25_r8*gxa*2.0_r8,             &
     &   -0.25_r8*(1.0_r8-gxb)*2.0_r8, -0.25_r8*gxb*2.0_r8,             &
     &    0.25_r8*(1.0_r8-gxb)*2.0_r8,  0.25_r8*gxb*2.0_r8,             &
     &   -0.25_r8*(1.0_r8-gxa)*2.0_r8, -0.25_r8*gxa*2.0_r8,             &
     &    0.25_r8*(1.0_r8-gxa)*2.0_r8,  0.25_r8*gxa*2.0_r8,             &
     &   -0.25_r8*(1.0_r8-gxb)*2.0_r8, -0.25_r8*gxb*2.0_r8,             &
     &    0.25_r8*(1.0_r8-gxb)*2.0_r8,  0.25_r8*gxb*2.0_r8              &
     &  /),(/4,4,2/))
!
!  SSA solver fields.  All work arrays below are declared on the GLOBAL
!  bounds (LBi:UBi, LBj:UBj) because they are passed to mp_exchange2d /
!  exchange_p2d_tile, which index them on those bounds (the tile scratch
!  bounds IminS:ImaxS, JminS:JmaxS are wider by the ghost margin).
!
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Ru, Rv
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Zu, Zv
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Du, Dv
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Au, Av
      real(r8), dimension(LBi:UBi,LBj:UBj) :: diag_u, diag_v
!
!  Extended cell-assembly bounds (see module header):  own RHO
!  cells Istr..Iend plus the western/southern halo ring that feeds
!  the seam PSI nodes owned by this tile.  The ring is skipped
!  along physical (non-periodic) domain edges, where the halo
!  cells beyond the boundary do not exist.
!
      integer :: IstrC, IendC, JstrC, JendC
!
# ifdef DISTRIBUTE
      character (len=3) :: op_handle
# endif
!
# include "set_bounds.h"
!
!  Extended cell-assembly bounds (see module header).
!
      IstrC=Istr
      IF (.not.DOMAIN(ng)%Western_Edge(tile).or.EWperiodic(ng)) THEN
        IstrC=Istr-1
      END IF
      JstrC=Jstr
      IF (.not.DOMAIN(ng)%Southern_Edge(tile).or.NSperiodic(ng)) THEN
        JstrC=Jstr-1
      END IF
      IendC=Iend
      JendC=Jend
!
!=======================================================================
!  Stage 1:  initial velocity guess.
!=======================================================================
!
!  Safety cap on the lsold ice thickness (robustness against garbage
!  values from initialization or numerical artifacts feeding the
!  geometry, the viscosity and the driving stress).  Configurable
!  via max_thickness in isf.in (default 5000 m).  The capped value
!  propagates to the new layer through the next advection step.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          IF (Si(i,j,lsold,isIhis).gt.max_thickness(ng)) THEN
            Si(i,j,lsold,isIhis)=max_thickness(ng)
          END IF
        END DO
      END DO
!
!  Copy the previous SSA velocity (lsold) as the initial guess and
!  zero the inactive nodes;  the Dirichlet lateral boundary values
!  are enforced later by isf_pibc_tile (after the warm start).
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          IF (pmask_act(i,j).eq.0.0_r8) THEN
            Si(i,j,lsnew,isIuis)=0.0_r8
            Si(i,j,lsnew,isIvis)=0.0_r8
          ELSE
            Si(i,j,lsnew,isIuis)=Si(i,j,lsold,isIuis)
            Si(i,j,lsnew,isIvis)=Si(i,j,lsold,isIvis)
          END IF
        END DO
      END DO
!
!=======================================================================
!  Stage 2:  driving stress RHS (FE-consistent nodal forces).
!=======================================================================
!
      CALL ssa_land_margin (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Istr, Iend, Jstr, Jend,                    &
     &                       IstrC, IendC, JstrC, JendC,                &
     &                       IstrP, IendP, JstrP, JendP,                &
     &                       h, zeta_avg, rmask_ist,                    &
     &                       Si(:,:,lsold,isIhis), pmarg)
!
      CALL ssa_driving_stress (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         pm, pn,                                  &
     &                         rmask_ist, pmask_act, pmarg,             &
     &                         Si(:,:,lsold,isIhis),                    &
     &                         Fi(:,:,ifIzsf),                          &
     &                         Fi(:,:,ifItdx),                          &
     &                         Fi(:,:,ifItdy))
!
      CALL ssa_front_load (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Istr, Iend, Jstr, Jend,                      &
     &                     IstrC, IendC, JstrC, JendC,                  &
     &                     IstrP, IendP, JstrP, JendP,                  &
     &                     pm, pn, h, zeta_avg,                         &
     &                     rmask_ist, pmask_act,                        &
     &                     Si(:,:,lsold,isIhis),                        &
     &                     Fi(:,:,ifItdx),                              &
     &                     Fi(:,:,ifItdy))
!
!=======================================================================
!  Stage 3:  basal-drag warm start.
!=======================================================================
!
!  If the initial velocity is near zero, invert the local basal drag
!  balance (tau_d = beta * u * A_node) to provide a first guess of
!  physically reasonable magnitude.  This avoids the Picard
!  "high-viscosity trap" (zero velocity -> eps_e = eps_glen_min ->
!  huge nuH -> tiny velocity -> convergence to the wrong fixed
!  point).  The inversion treats the frozen beta linearly, so for
!  sliding exponents other than the linear law it is a heuristic
!  initial estimate, not an exact balance.  On floating ice
!  (beta = 0) the warm start is skipped and the previous velocity
!  field is kept.  The node balance beta*A_node*u = b_k uses the
!  same lumped nodal area A_node = jac_avg as the driving stress.
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          IF (pmask_act(i,j).eq.1.0_r8) THEN
            u_slide=SQRT(Si(i,j,lsnew,isIuis)**2+                       &
     &                   Si(i,j,lsnew,isIvis)**2)
            IF (u_slide.lt.u_warm_min) THEN
              beta_avg=0.25_r8*(Fi(i  ,j  ,ifIbta)+                     &
     &                          Fi(i-1,j  ,ifIbta)+                     &
     &                          Fi(  i,j-1,ifIbta)+                     &
     &                          Fi(i-1,j-1,ifIbta))
              IF (beta_avg.gt.1.0E-10_r8) THEN
                jac_avg=0.25_r8*(                                       &
     &            1.0_r8/(pm(i-1,j-1)*pn(i-1,j-1))+                     &
     &            1.0_r8/(pm(i,  j-1)*pn(i,  j-1))+                     &
     &            1.0_r8/(pm(i-1,j  )*pn(i-1,j  ))+                     &
     &            1.0_r8/(pm(i,  j  )*pn(i,  j  )))
!  Node balance beta_eff * A_node * u = b_k with the same lumped
!  nodal area A_node = jac_avg as the driving stress assembly.
!
                cff=1.0_r8/(beta_avg*jac_avg)
                Si(i,j,lsnew,isIuis)=Fi(i,j,ifItdx)*cff
                Si(i,j,lsnew,isIvis)=Fi(i,j,ifItdy)*cff
              END IF
            END IF
          END IF
        END DO
      END DO
!
      CALL isf_pibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    Si(:,:,lsnew,isIuis),                         &
     &                    Si(:,:,lsnew,isIvis))
!
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Si(:,:,lsnew,isIuis),                         &
     &                    Si(:,:,lsnew,isIvis))
# endif
!
!=======================================================================
!  Stage 4:  Glen hardness field (constant during the solve).
!=======================================================================
!
!  B = A^(-1/n) depends only on the ice temperature (frozen at the
!  coupling-window diagnostic) and the lsold thickness, so it is
!  computed once per solve;  ssa_viscosity then only combines it
!  with the current strain-rate invariant (ssa_hardness below).
!
      CALL ssa_hardness (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Istr, Iend, Jstr, Jend,                        &
     &                   Si(:,:,lsold,isIhis),                          &
     &                   Fi(:,:,ifItis),                                &
     &                   B_hard)
!
!=======================================================================
!  Stage 5:  initialize CG work arrays.
!=======================================================================
!
      DO j=LBj,UBj
        DO i=LBi,UBi
          Au(i,j)=0.0_r8
          Av(i,j)=0.0_r8
          diag_u(i,j)=0.0_r8
          diag_v(i,j)=0.0_r8
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
!  Stage 6:  Picard iteration.
!=======================================================================
!
!  Each iteration freezes the viscosity at the current velocity,
!  solves the linearized SSA with matrix-free PCG, applies the
!  robustness measures (speed cap, relaxation) and checks the dual
!  convergence criterion;  a warning is printed if nl_max_iter is
!  exhausted without either criterion being met.
!
      cgtol=cg_tol(ng)
      converged=.FALSE.
!
      DO nl_iter=1,nl_max_iter(ng)
!
!  Effective viscosity at RHO points from the current velocity
!  estimate (Glen law with the Stage-4 cached hardness field,
!  thin-ice strength extension and additive regularization).
!
        CALL ssa_viscosity (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Istr, Iend, Jstr, Jend,                     &
     &                      pm, pn,                                     &
     &                      Si(:,:,lsold,isIhis),                       &
     &                      Si(:,:,lsnew,isIuis),                       &
     &                      Si(:,:,lsnew,isIvis),                       &
     &                      B_hard,                                     &
     &                      Fi(:,:,ifInuh))
!
!  Jacobi preconditioner:  the diagonal of the assembled SSA
!  operator at the active PSI nodes (ssa_jacobi_diag).
!
        CALL ssa_jacobi_diag (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IstrC, IendC, JstrC, JendC,               &
     &                        IstrP, IendP, JstrP, JendP,               &
     &                        DPhi, phi_tab,                            &
     &                        pm, pn,                                   &
     &                        pmask_act,                                &
     &                        Fi(:,:,ifInuh),                           &
     &                        Fi(:,:,ifIbta),                           &
     &                        diag_u, diag_v)
!
        CALL ssa_matvec (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IstrC, IendC, JstrC, JendC,                    &
     &                   IstrP, IendP, JstrP, JendP,                    &
     &                   DPhi, phi_tab,                                 &
     &                   pm, pn,                                        &
     &                   pmask_act,                                     &
     &                   Fi(:,:,ifInuh),                                &
     &                   Fi(:,:,ifIbta),                                &
     &                   Si(:,:,lsnew,isIuis),                          &
     &                   Si(:,:,lsnew,isIvis),                          &
     &                   Au, Av)
!
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            IF (pmask_act(i,j).eq.1.0_r8) THEN
              Ru(i,j)=Fi(i,j,ifItdx)-Au(i,j)
              Rv(i,j)=Fi(i,j,ifItdy)-Av(i,j)
            ELSE
              Ru(i,j)=0.0_r8
              Rv(i,j)=0.0_r8
            END IF
          END DO
        END DO
!
!  Zero the residual on Dirichlet-constrained components (edges
!  flagged closed or clamped in the ice-state LBC) so the CG
!  update never touches those nodes.
!
        CALL ssa_lbc_zero (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     Ru, Rv)
!
!  Preconditioned residual and initial search direction.
!
        dot_zr=0.0_r8
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            IF (pmask_act(i,j).eq.1.0_r8) THEN
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
        CALL ssa_lbc_zero (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
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
!  Matrix action on the SEARCH DIRECTION, A*d.  (The action on the
!  velocity itself is only needed for the initial and the
!  true-residual evaluations, which use separate matvec calls.)
!
          CALL ssa_matvec (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IstrC, IendC, JstrC, JendC,                  &
     &                     IstrP, IendP, JstrP, JendP,                  &
     &                     DPhi, phi_tab,                               &
     &                     pm, pn,                                      &
     &                     pmask_act,                                   &
     &                     Fi(:,:,ifInuh),                              &
     &                     Fi(:,:,ifIbta),                              &
     &                     Du, Dv,                                      &
     &                     Au, Av)
!
!  alpha = (z,r) / (d,A*d):  dot_zr is carried from the previous
!  iteration (or the initial evaluation) instead of being recomputed,
!  so only dot_dad needs the fused pass over the active nodes.
!
          dot_dad=0.0_r8
          dot_zr_old=dot_zr
          DO j=JstrP,JendP
            DO i=IstrP,IendP
              IF (pmask_act(i,j).eq.1.0_r8) THEN
                dot_dad=dot_dad+Du(i,j)*Au(i,j)+Dv(i,j)*Av(i,j)
              END IF
            END DO
          END DO
!
# ifdef DISTRIBUTE
          op_handle='SUM'
          CALL mp_reduce (ng, model, 1, dot_dad, op_handle)
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
              IF (pmask_act(i,j).eq.1.0_r8) THEN
                Si(i,j,lsnew,isIuis)=Si(i,j,lsnew,isIuis)+alpha*Du(i,j)
                Si(i,j,lsnew,isIvis)=Si(i,j,lsnew,isIvis)+alpha*Dv(i,j)
                Ru(i,j)=Ru(i,j)-alpha*Au(i,j)
                Rv(i,j)=Rv(i,j)-alpha*Av(i,j)
              END IF
            END DO
          END DO
!
          CALL ssa_lbc_zero (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       Ru, Rv)
!
!  Preconditioned residual update and (z,r) accumulation
!  (Fletcher-Reeves form:  beta = (z_k+1,r_k+1)/(z_k,r_k)).
!
          dot_zr=0.0_r8
          DO j=JstrP,JendP
            DO i=IstrP,IendP
              IF (pmask_act(i,j).eq.1.0_r8) THEN
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
              IF (pmask_act(i,j).eq.1.0_r8) THEN
                Du(i,j)=Zu(i,j)+beta_cg*Du(i,j)
                Dv(i,j)=Zv(i,j)+beta_cg*Dv(i,j)
              END IF
            END DO
          END DO
!
          CALL ssa_lbc_zero (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       Du, Dv)
!
          resid=SQRT(ABS(dot_zr))
!
        END DO  ! cg_iter
!
!-----------------------------------------------------------------------
!  Post-CG robustness measures:  speed cap, Picard relaxation, then
!  a single pibc + halo refresh (below) before the consistency
!  viscosity and the true-residual matvec.  No exchange is needed
!  between the CG loop and the cap:  the cap and the relaxation
!  read and write interior nodes only, and CG never moves the
!  Dirichlet nodes (ssa_lbc_zero keeps their search direction
!  zero), so the boundary values are still in place.
!-----------------------------------------------------------------------
!
!  PISM-style velocity cap.
!
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            IF (pmask_act(i,j).eq.1.0_r8) THEN
              u_slide=SQRT(Si(i,j,lsnew,isIuis)**2+                     &
     &                     Si(i,j,lsnew,isIvis)**2)
              IF (u_slide.gt.ssa_max_speed(ng)) THEN
                cff=ssa_max_speed(ng)/u_slide
                Si(i,j,lsnew,isIuis)=Si(i,j,lsnew,isIuis)*cff
                Si(i,j,lsnew,isIvis)=Si(i,j,lsnew,isIvis)*cff
              END IF
            END IF
          END DO
        END DO
!
!  Picard relaxation u = alpha*u_new + (1-alpha)*u_old.  Prevents
!  oscillation between high/low viscosity regimes.  Skipped on
!  the first iteration, where lsold holds the (physical) initial
!  guess rather than a Picard iterate; later iterations anchor
!  to the previous iterate written back into lsold.
!
        IF (nl_iter.gt.1) THEN
          DO j=JstrP,JendP
            DO i=IstrP,IendP
              IF (pmask_act(i,j).eq.1.0_r8) THEN
                Si(i,j,lsnew,isIuis)=                                   &
     &                    picard_alpha(ng) *Si(i,j,lsnew,isIuis)+       &
     &            (1.0_r8-picard_alpha(ng))*Si(i,j,lsold,isIuis)
                Si(i,j,lsnew,isIvis)=                                   &
     &                    picard_alpha(ng) *Si(i,j,lsnew,isIvis)+       &
     &            (1.0_r8-picard_alpha(ng))*Si(i,j,lsold,isIvis)
              END IF
            END DO
          END DO
        END IF
!
        CALL isf_pibc_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      Si(:,:,lsnew,isIuis),                       &
     &                      Si(:,:,lsnew,isIvis))
!
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, model, 2,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      Si(:,:,lsnew,isIuis),                       &
     &                      Si(:,:,lsnew,isIvis))
# endif
!
!  Refresh the velocity halo and the Dirichlet BCs after the cap
!  and relaxation, then recompute the viscosity so nuH is
!  consistent with the relaxed velocity used in the next Picard
!  iteration and in the convergence check below.
!
        CALL ssa_viscosity (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Istr, Iend, Jstr, Jend,                     &
     &                      pm, pn,                                     &
     &                      Si(:,:,lsold,isIhis),                       &
     &                      Si(:,:,lsnew,isIuis),                       &
     &                      Si(:,:,lsnew,isIvis),                       &
     &                      B_hard,                                     &
     &                      Fi(:,:,ifInuh))
!
!  True nonlinear residual ||A(u)u - b|| with the relaxed velocity
!  and the updated viscosity (one extra matrix action per Picard
!  iteration;  the legacy solver reused the last CG search
!  direction action, which does not measure the SSA residual).
!
        CALL ssa_matvec (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IstrC, IendC, JstrC, JendC,                    &
     &                   IstrP, IendP, JstrP, JendP,                    &
     &                   DPhi, phi_tab,                                 &
     &                   pm, pn,                                        &
     &                   pmask_act,                                     &
     &                   Fi(:,:,ifInuh),                                &
     &                   Fi(:,:,ifIbta),                                &
     &                   Si(:,:,lsnew,isIuis),                          &
     &                   Si(:,:,lsnew,isIvis),                          &
     &                   Au, Av)
!
!-----------------------------------------------------------------------
!  Picard convergence check:  dual criterion.
!    1. Relative velocity change
!           ||u_relaxed - u_old|| / ||u_relaxed||  <  nl_tol
!    2. Relative residual norm
!                   ||A(u)u - b|| / ||b||          <  nl_tol
!  Either criterion passing triggers convergence.  The residual
!  criterion catches stagnation at a non-solution fixed point.
!-----------------------------------------------------------------------
!
        norms(1)=0.0_r8
        norms(2)=0.0_r8
        norms(3)=0.0_r8
        norms(4)=0.0_r8
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            IF (pmask_act(i,j).eq.1.0_r8) THEN
              norms(1)=norms(1)+                                        &
     &          (Si(i,j,lsnew,isIuis)-Si(i,j,lsold,isIuis))**2+         &
     &          (Si(i,j,lsnew,isIvis)-Si(i,j,lsold,isIvis))**2
              norms(2)=norms(2)+                                        &
     &          Si(i,j,lsnew,isIuis)**2+Si(i,j,lsnew,isIvis)**2
              norms(3)=norms(3)+                                        &
     &          (Au(i,j)-Fi(i,j,ifItdx))**2+(Av(i,j)-Fi(i,j,ifItdy))**2
              norms(4)=norms(4)+                                        &
     &          Fi(i,j,ifItdx)**2+Fi(i,j,ifItdy)**2
            END IF
          END DO
        END DO
# ifdef DISTRIBUTE
        op_norms='SUM'
        CALL mp_reduce (ng, model, 4, norms, op_norms)
# endif
        vel_change=norms(1)
        vel_norm=norms(2)
        resid_norm=norms(3)
        rhs_norm=norms(4)
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
# ifdef ISF_DEBUG
!
!  DEBUG:  Picard iteration diagnostics.  Entirely gated by
!  ISF_DEBUG so production runs pay no output cost.
!
        IF (Master) THEN
          IF (converged.or.(nl_iter.eq.1).or.MOD(nl_iter,10).eq.0) THEN
            u_slide=0.0_r8
            DO j=JstrP,JendP
              DO i=IstrP,IendP
                IF (pmask_act(i,j).eq.1.0_r8) THEN
                  cff1=SQRT(Si(i,j,lsnew,isIuis)**2+                    &
     &                      Si(i,j,lsnew,isIvis)**2)
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
!  iteration.  NOTE: this overwrites the lsold VELOCITY layer, so
!  after the solve it holds the penultimate Picard iterate rather
!  than the previous SSA solution.  That is safe here:  every
!  downstream consumer reads the solved lsnew layer (advection
!  carry-forward, isf_stress beta, output), and after the next
!  time-level swap the overwritten layer is refreshed before any
!  read;  under PERFECT_RESTART both layers are written, so even
!  the Picard anchor chain is reproduced exactly after a restart.
!
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            Si(i,j,lsold,isIuis)=Si(i,j,lsnew,isIuis)
            Si(i,j,lsold,isIvis)=Si(i,j,lsnew,isIvis)
          END DO
        END DO
!
        CALL isf_pibc_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      Si(:,:,lsold,isIuis),                       &
     &                      Si(:,:,lsold,isIvis))
!
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, model, 2,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      Si(:,:,lsold,isIuis),                       &
     &                      Si(:,:,lsold,isIvis))
# endif
!
      END DO  ! nl_iter (Picard)
!
!  Warn when the Picard iteration is exhausted without meeting
!  either convergence criterion:  the (unconverged) last iterate
!  is returned, and the ocean coupling must not treat it as a
!  silently converged solution.
!
      IF (.not.converged .and. Master) THEN
        WRITE(stdout,'(a,i4,a,es9.2)')                                  &
     &    ' SSA WARNING: Picard not converged in ',nl_max_iter(ng),     &
     &    ' iterations; last nl_resid=',nl_resid
      END IF
!
!=======================================================================
!  Stage 7:  NaN guard and final halo refresh.
!=======================================================================
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          IF (pmask_act(i,j).eq.1.0_r8) THEN
            IF (Si(i,j,lsnew,isIuis).ne.Si(i,j,lsnew,isIuis).or.        &
     &          Si(i,j,lsnew,isIvis).ne.Si(i,j,lsnew,isIvis)) THEN
              WRITE(stdout,'(a,i3,a,2i4)')                              &
     &            ' SSA NaN zeroed: rank=',MyRank,                      &
     &            '  tile-local PSI node (i,j)=',i,j
              Si(i,j,lsnew,isIuis)=0.0_r8
              Si(i,j,lsnew,isIvis)=0.0_r8
            END IF
          END IF
        END DO
      END DO
!
      CALL isf_pibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    Si(:,:,lsnew,isIuis),                         &
     &                    Si(:,:,lsnew,isIvis))
!
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Si(:,:,lsnew,isIuis),                         &
     &                    Si(:,:,lsnew,isIvis))
# endif
!
      RETURN
      END SUBROUTINE isf_ssa_tile
!
!***********************************************************************
      SUBROUTINE ssa_lbc_zero (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         uf, vf)
!***********************************************************************
!
!  Zero the Dirichlet-constrained components of a general PSI-point
!  field pair.  Used to keep the CG residual and search direction
!  zero on fixed nodes so the iteration never updates them.  The
!  initial guess already satisfies the boundary values (written
!  by isf_pibc_tile); zero residual then keeps the solution at
!  those values throughout the iteration.
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: uf(LBi:,LBj:)
      real(r8), intent(inout) :: vf(LBi:,LBj:)
# else
      real(r8), intent(inout) :: uf(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vf(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Zero fixed components.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
!
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
          IF (LBC(iwest,ibISF(isIuis),ng)%clamped .or.                  &
     &        LBC(iwest,ibISF(isIuis),ng)%closed) THEN
            DO j=JstrP,JendP
              uf(IstrP,j)=0.0_r8
            END DO
          END IF
!
          IF (LBC(iwest,ibISF(isIvis),ng)%clamped .or.                  &
     &        LBC(iwest,ibISF(isIvis),ng)%closed) THEN
            DO j=JstrP,JendP
              vf(IstrP,j)=0.0_r8
            END DO
          END IF
        END IF
!
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
!
          IF (LBC(ieast,ibISF(isIuis),ng)%clamped .or.                  &
     &        LBC(ieast,ibISF(isIuis),ng)%closed) THEN
            DO j=JstrP,JendP
              uf(IendP,j)=0.0_r8
            END DO
          END IF
!
          IF (LBC(ieast,ibISF(isIvis),ng)%clamped .or.                  &
      &       LBC(ieast,ibISF(isIvis),ng)%closed) THEN
            DO j=JstrP,JendP
              vf(IendP,j)=0.0_r8
            END DO
          END IF
        END IF
      END IF

      IF (.not.NSperiodic(ng)) THEN
!
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!
          IF (LBC(isouth,ibISF(isIuis),ng)%clamped .or.                 &
     &        LBC(isouth,ibISF(isIuis),ng)%closed) THEN
            DO i=IstrP,IendP
              uf(i,JstrP)=0.0_r8
            END DO
          END IF
!
          IF (LBC(isouth,ibISF(isIvis),ng)%clamped .or.                 &
     &        LBC(isouth,ibISF(isIvis),ng)%closed) THEN
            DO i=IstrP,IendP
              vf(i,JstrP)=0.0_r8
            END DO
          END IF
        END IF
!
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
!
          IF (LBC(inorth,ibISF(isIuis),ng)%clamped .or.                 &
     &        LBC(inorth,ibISF(isIuis),ng)%closed) THEN
            DO i=IstrP,IendP
              uf(i,JendP)=0.0_r8
            END DO
          END IF
!
          IF (LBC(inorth,ibISF(isIvis),ng)%clamped .or.                 &
     &        LBC(inorth,ibISF(isIvis),ng)%closed) THEN
            DO i=IstrP,IendP
              vf(i,JendP)=0.0_r8
            END DO
          END IF
        END IF
      END IF
!
      RETURN
      END SUBROUTINE ssa_lbc_zero
!
!***********************************************************************
      SUBROUTINE ssa_driving_stress (ng, tile, model,                   &
     &                               LBi, UBi, LBj, UBj,                &
     &                               pm, pn,                            &
     &                               rmask_ist, pmask_act, pmarg,       &
     &                               hi, zs,                            &
     &                               tdx, tdy)
!***********************************************************************
!
!  Compute the driving stress RHS at PSI points:
!
!    b_k = -rho_i * g * H_psi * grad(S) * (nodal area)
!
!  This is the mass-lumped Q1 nodal force vector:  a finite
!  difference gradient of the smoothed surface elevation on the 4
!  RHO cells surrounding each PSI node, multiplied by the lumped
!  nodal area jac_avg = (A_1+A_2+A_3+A_4)/4 = sum_e (A_e/4), the
!  exact integral of the Q1 hat function over its 4-element
!  support.  This matches the consistent-load convention of PISM
!  SSAFEM (int phi_k tau_d dOmega by quadrature) and of the
!  MITgcm STREAMICE driving stress.  NOTE: the legacy
!  isf_driving_stress.h carries an extra 0.25 factor (a 4x too
!  small nodal load); correcting it here makes the new solver's
!  velocities larger than the legacy solver's by up to 4x -- an
!  expected A/B difference until the legacy factor is corrected
!  or the legacy path is retired.
!
!  The surface elevation S (Fi(:,:,ifIzsf)) is smoothed with
!  n_smooth iterations of a 5-point Laplacian filter before the
!  gradient is taken.  This regularizes the grounding-line kink
!  where S transitions between the grounded (S = -h + H) and
!  floating (S = (1-rho_i/rho_w) H) branches.
!
!  Land-margin taper (plan option A):  the driving stress is
!  linearly tapered to zero within n_taper cells of a LAND margin
!  (grounded ice with bed at or above sea level terminating on
!  ice-free cells;  indicator pmarg from ssa_land_margin).  Marine
!  fronts (floating or grounded in water) carry the full driving
!  stress plus the front-pressure load and are NOT tapered -- the
!  cliff force is represented by that load alone (see the constant
!  extension below).  n_taper <= 0 disables the taper (guard added:
!  the previous 0/0 produced NaN for the default n_taper = 0).
!
!  Reference:  PISM compute_nuH_everywhere (slope regularization),
!  MITgcm STREAMICE streamice_driving_stress.F.
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: pm(LBi:,LBj:)
      real(r8), intent(in)  :: pn(LBi:,LBj:)
      real(r8), intent(in)  :: rmask_ist(LBi:,LBj:)
      real(r8), intent(in)  :: pmask_act(LBi:,LBj:)
      real(r8), intent(in)  :: pmarg(LBi:,LBj:)
      real(r8), intent(in)  :: hi(LBi:,LBj:)
      real(r8), intent(in)  :: zs(LBi:,LBj:)
      real(r8), intent(out) :: tdx(LBi:,LBj:)
      real(r8), intent(out) :: tdy(LBi:,LBj:)
# else
      real(r8), intent(in)  :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: rmask_ist(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pmask_act(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pmarg(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: hi(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: zs(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tdx(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tdy(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, di, dj, is_iter, dist_margin
!
      real(r8) :: hi_psi
      real(r8) :: dSdx, dSdy, dx_rho, dy_rho, jac_avg
      real(r8) :: zs_sum, taper_w
      real(r8) :: cff, cff1
!
      real(r8), dimension(LBi:UBi,LBj:UBj) :: zs_smooth
      real(r8), dimension(LBi:UBi,LBj:UBj) :: zs_smooth_temp
!
#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize the nodal forces.
!-----------------------------------------------------------------------
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          tdx(i,j)=0.0_r8
          tdy(i,j)=0.0_r8
        END DO
      END DO
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
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          IF (rmask_ist(i,j).eq.1.0_r8) THEN
            zs_smooth(i,j)=zs(i,j)
          ELSE
            zs_smooth(i,j)=0.0_r8
          END IF
        END DO
      END DO

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    zs_smooth)
# endif
!
!  Constant extension of the surface elevation into ice-free cells
!  adjacent to ice.  The ice-front cliff is represented by the
!  front-pressure boundary load ONLY (ssa_front_load);  without this
!  extension the FD gradient at margin nodes would see the zs = 0
!  cliff drop of the ice-free side and add a spurious nodal force
!  of order 2*rho_i*g*(1-rho_i/rho0)*H**2*dy -- twice the front load
!  it replaces.  See .agents/docs/plans/ssa_front_load_gradient_
!  rationale.md for the full derivation.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          IF (rmask_ist(i,j).ne.1.0_r8) THEN
            cff=0.0_r8
            cff1=0.0_r8
!
!  Same bounds guards as the smoothing gates above:  the +/-1
!  neighbour reads must stay inside the array on domain-edge
!  tiles (out-of-bounds = no such neighbour = skip).
!
            IF (i.gt.LBi.and.rmask_ist(i-1,j).eq.1.0_r8) THEN
              cff=cff+zs_smooth(i-1,j)
              cff1=cff1+1.0_r8
            END IF
            IF (i.lt.UBi.and.rmask_ist(i+1,j).eq.1.0_r8) THEN
              cff=cff+zs_smooth(i+1,j)
              cff1=cff1+1.0_r8
            END IF
            IF (j.gt.LBj.and.rmask_ist(i,j-1).eq.1.0_r8) THEN
              cff=cff+zs_smooth(i,j-1)
              cff1=cff1+1.0_r8
            END IF
            IF (j.lt.UBj.and.rmask_ist(i,j+1).eq.1.0_r8) THEN
              cff=cff+zs_smooth(i,j+1)
              cff1=cff1+1.0_r8
            END IF
            IF (cff1.gt.0.0_r8) THEN
              zs_smooth(i,j)=cff/cff1
            END IF
          END IF
        END DO
      END DO
!
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    zs_smooth)
# endif
!
!  Apply the 5-point Laplacian smoothing (n_smooth iterations).
!  Only ice-covered cells are smoothed;  ice-free neighbours are
!  replaced by the centre value (no-flux treatment at the margin).
!  The interior-only loop keeps the writes inside the array bounds.
!
      DO is_iter=1,n_smooth(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            IF (rmask_ist(i,j).eq.1.0_r8) THEN
              zs_sum=zs_smooth(i,j)
!
!  Gate guards:  the +/-1 neighbours reach one cell beyond the
!  JstrR..JendR / IstrR..IendR loop, which on domain-edge tiles
!  sits OUTSIDE the array bounds (e.g. LBj = 0 with JstrR = 0:
!  the south gate would read rmask_ist(i,-1)).  Those reads are
!  undefined behaviour fetching adjacent heap/stack memory;  when
!  the fetched bits happened to equal 1.0_r8 the gate opened and
!  an out-of-bounds zs_smooth value entered the average (the L11
!  nondeterminism).  An out-of-bounds neighbour is treated as
!  ice-free (no-flux), which is exactly the boundary semantics.
!
              IF (i.lt.UBi.and.                                         &
     &            rmask_ist(i+1,j  ).eq.1.0_r8) THEN
                zs_sum=zs_sum+zs_smooth(i+1,j  )
              ELSE
                zs_sum=zs_sum+zs_smooth(i  ,j  )
              END IF
              IF (i.gt.LBi.and.                                         &
     &            rmask_ist(i-1,j  ).eq.1.0_r8) THEN
                zs_sum=zs_sum+zs_smooth(i-1,j  )
              ELSE
                zs_sum=zs_sum+zs_smooth(i  ,j  )
              END IF
              IF (j.lt.UBj.and.                                         &
     &            rmask_ist(i  ,j+1).eq.1.0_r8) THEN
                zs_sum=zs_sum+zs_smooth(i  ,j+1)
              ELSE
                zs_sum=zs_sum+zs_smooth(i  ,j  )
              END IF
              IF (j.gt.LBj.and.                                         &
     &            rmask_ist(i  ,j-1).eq.1.0_r8) THEN
                zs_sum=zs_sum+zs_smooth(i  ,j-1)
              ELSE
                zs_sum=zs_sum+zs_smooth(i  ,j  )
              END IF
              zs_smooth_temp(i,j)=zs_sum/5.0_r8
            ELSE
!
!  Keep the constant-extension value (above) instead of resetting
!  to zero:  the gradient stencil at margin nodes reads these cells.
!
              zs_smooth_temp(i,j)=zs_smooth(i,j)
            END IF
          END DO
        END DO
!
!  Reset zs_smooth
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            zs_smooth(i,j)=zs_smooth_temp(i,j)
          END DO
        END DO
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      zs_smooth)
# endif
      END DO
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
          IF (pmask_act(i,j).eq.1.0_r8) THEN
!
!  Ice margin taper:  distance to the nearest ice-free cell among
!  the 8 neighbours, normalized by n_taper and clamped to 1.
!
            IF (n_taper(ng).le.0) THEN
              taper_w=1.0_r8
            ELSE
              dist_margin=n_taper(ng)
              DO dj=-1,1
                DO di=-1,1
                  IF (di.ne.0.or.dj.ne.0) THEN
!
!  Bounds guard:  JendP can equal UBj on the northern domain
!  edge, where the +1 scan would read pmarg(i,UBj+1) out of
!  bounds (same undefined-memory hazard as the gates above).
!
                    IF (i+di.ge.LBi.and.i+di.le.UBi.and.               &
     &                  j+dj.ge.LBj.and.j+dj.le.UBj) THEN
                      IF (pmarg(i+di,j+dj).eq.1.0_r8) THEN
                        dist_margin=MIN(dist_margin,                    &
     &                                  MAX(ABS(di),ABS(dj)))
                      END IF
                    END IF
                  END IF
                END DO
              END DO
              taper_w=MIN(REAL(dist_margin,r8)/REAL(n_taper(ng),r8),    &
     &                    1.0_r8)
            END IF
!
!  FD gradient of the smoothed surface at the PSI node.
!
            dx_rho=0.5_r8*(1.0_r8/pm(i-1,j  )+1.0_r8/pm(i  ,j  ))
            dy_rho=0.5_r8*(1.0_r8/pn(i  ,j-1)+1.0_r8/pn(i  ,j  ))
            jac_avg=0.25_r8*(1.0_r8/(pm(i-1,j-1)*pn(i-1,j-1))+          &
     &                       1.0_r8/(pm(i,  j-1)*pn(i,  j-1))+          &
     &                       1.0_r8/(pm(i-1,j  )*pn(i-1,j  ))+          &
     &                       1.0_r8/(pm(i,  j  )*pn(i,  j  )))
!
            dSdx=((zs_smooth(i  ,j-1)+zs_smooth(i  ,j  ))-              &
     &            (zs_smooth(i-1,j-1)+zs_smooth(i-1,j  )))/             &
     &           (2.0_r8*dx_rho)
            dSdy=((zs_smooth(i-1,j  )+zs_smooth(i  ,j  ))-              &
     &            (zs_smooth(i-1,j-1)+zs_smooth(i  ,j-1)))/             &
     &           (2.0_r8*dy_rho)
!
            hi_psi=0.25_r8*(hi(i-1,j-1)+hi(i,j-1)+hi(i-1,j)+hi(i,j))
!
!  FE-consistent nodal force:  tau_d times the lumped nodal area
!  jac_avg (= sum of the surrounding A_e/4;  see header).
!
            cff=-taper_w*rho_i*g*hi_psi*jac_avg
            tdx(i,j)=cff*dSdx
            tdy(i,j)=cff*dSdy
          END IF
        END DO
      END DO
!
!  No halo exchange of the assembled RHS:  every consumer (warm
!  start, CG residual, convergence norms) reads tdx/tdy at
!  owned PSI nodes only.
!
      RETURN
      END SUBROUTINE ssa_driving_stress
!
!***********************************************************************
      SUBROUTINE ssa_front_load (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           IstrC, IendC, JstrC, JendC,            &
     &                           IstrP, IendP, JstrP, JendP,            &
     &                           pm, pn, h, zeta_avg,                   &
     &                           rmask_ist, pmask_act, hi,              &
     &                           tdx, tdy)
!***********************************************************************
!
!  Ice-front pressure boundary load:  the depth-integrated net
!  hydrostatic push acting on each ice-front face, scattered to the
!  face's two PSI nodes as a Q1 line load.  Each endpoint receives
!  p_hat*L_face/2 along the OUTWARD normal (seaward); the tangential
!  component is zero (N_nt = 0 for uniform H along the front, the
!  same approximation as MITgcm STREAMICE).
!
!  A front face separates an ice cell (rmask_ist = 1) from an
!  ice-free neighbour cell;  a free (gradient) lateral boundary
!  counts as ice-free (the ice continues into the open domain),
!  while closed/clamped edges are walls and carry no load.  The
!  load case depends only on the ICE cell's own state, evaluated
!  with the same window-averaged criterion as the basal drag
!  (isf_stress),  D_avg = h + zeta_avg - Dcrit:
!
!    floating, rho_i*H <= rho0*D_avg :
!        p_hat = 0.5*rho_i*g*(1-rho_i/rho0)*H**2
!    grounded in water, rho_i*H > rho0*D_avg > 0 :
!        p_hat = 0.5*g*(rho_i*H**2 - rho0*D_avg**2)
!    grounded land margin, D_avg <= 0 :
!        p_hat = 0.5*rho_i*g*H**2
!
!  A consistently classified grounded cell always satisfies
!  rho_i*H**2 > rho0*D_avg**2 (grounding implies rho_i*H > rho0*D_avg
!  and H/rho_i > 0), so p_hat > 0 without any clamp.  The load
!  accumulates on the driving-stress RHS (tdx/tdy) assembled by
!  ssa_driving_stress, so the warm start and the Picard residual
!  norms consume it automatically.
!
!  The cell loop covers the own cells plus the western/southern
!  halo ring (IstrC/JstrC), and scatters only to owned PSI nodes
!  (ownership filter), so the assembly is complete under multi-tile
!  MPI by the same argument as ssa_matvec.  No halo exchange of the
!  output:  every consumer reads tdx/tdy at owned PSI nodes only.
!
!  References:  Morland (1987), MacAyeal (1989), Schoof (2007);
!  MITgcm STREAMICE streamice_driving_stress.F:124-211 (float_frac=1
!  is the grounded branch there; see ssa_front_load_plan erratum).
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Istr, Iend, Jstr, Jend
      integer, intent(in) :: IstrC, IendC, JstrC, JendC
      integer, intent(in) :: IstrP, IendP, JstrP, JendP
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: pm(LBi:,LBj:)
      real(r8), intent(in)  :: pn(LBi:,LBj:)
      real(r8), intent(in)  :: h(LBi:,LBj:)
      real(r8), intent(in)  :: zeta_avg(LBi:,LBj:)
      real(r8), intent(in)  :: rmask_ist(LBi:,LBj:)
      real(r8), intent(in)  :: pmask_act(LBi:,LBj:)
      real(r8), intent(in)  :: hi(LBi:,LBj:)
      real(r8), intent(inout) :: tdx(LBi:,LBj:)
      real(r8), intent(inout) :: tdy(LBi:,LBj:)
# else
      real(r8), intent(in)  :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: zeta_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: rmask_ist(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pmask_act(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: hi(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tdx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tdy(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
!
      logical :: front_e, front_w, front_n, front_s
!
      real(r8) :: H_ice, D_avg, phat
      real(r8) :: dyf, dxf, cff
!
!-----------------------------------------------------------------------
!  Assemble the front-pressure line loads.
!-----------------------------------------------------------------------
!
      DO j=JstrC,JendC
        DO i=IstrC,IendC
          IF (rmask_ist(i,j).eq.1.0_r8) THEN
            H_ice=hi(i,j)
            D_avg=h(i,j)+zeta_avg(i,j)-Dcrit(ng)
!
!  Net depth-integrated front pressure of this ice cell.
!
            IF (rho_i*H_ice.le.rho0*D_avg) THEN
              phat=0.5_r8*rho_i*g*(1.0_r8-rho_i/rho0)*H_ice**2
            ELSE IF (D_avg.gt.0.0_r8) THEN
              phat=0.5_r8*g*(rho_i*H_ice**2-rho0*D_avg**2)
            ELSE
              phat=0.5_r8*rho_i*g*H_ice**2
            END IF
!
!  Front-face detection:  ice-free neighbour, with a free (gradient)
!  lateral boundary counting as ice-free (shared predicate).
!
            CALL ssa_front_faces (ng, tile, i, j,                       &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Istr, Iend, Jstr, Jend,               &
     &                            rmask_ist,                            &
     &                            front_e, front_w, front_n, front_s)
!
!  Q1 line-load scatter:  each endpoint node of the face receives
!  p_hat*L_face/2 along the outward normal.
!
            IF (front_e) THEN
              dyf=0.5_r8*(1.0_r8/pn(i,j)+1.0_r8/pn(i+1,j))
              cff=0.5_r8*phat*dyf
              IF (pmask_act(i+1,j).eq.1.0_r8.and.                       &
     &            i+1.ge.IstrP.and.i+1.le.IendP.and.                    &
     &            j.ge.JstrP.and.j.le.JendP) THEN
                tdx(i+1,j)=tdx(i+1,j)+cff
              END IF
              IF (pmask_act(i+1,j+1).eq.1.0_r8.and.                     &
     &            i+1.ge.IstrP.and.i+1.le.IendP.and.                    &
     &            j+1.ge.JstrP.and.j+1.le.JendP) THEN
                tdx(i+1,j+1)=tdx(i+1,j+1)+cff
              END IF
            END IF
!
            IF (front_w) THEN
              dyf=0.5_r8*(1.0_r8/pn(i-1,j)+1.0_r8/pn(i,j))
              cff=0.5_r8*phat*dyf
              IF (pmask_act(i,j).eq.1.0_r8.and.                         &
     &            i.ge.IstrP.and.i.le.IendP.and.                        &
     &            j.ge.JstrP.and.j.le.JendP) THEN
                tdx(i,j)=tdx(i,j)-cff
              END IF
              IF (pmask_act(i,j+1).eq.1.0_r8.and.                       &
     &            i.ge.IstrP.and.i.le.IendP.and.                        &
     &            j+1.ge.JstrP.and.j+1.le.JendP) THEN
                tdx(i,j+1)=tdx(i,j+1)-cff
              END IF
            END IF
!
            IF (front_n) THEN
              dxf=0.5_r8*(1.0_r8/pm(i,j)+1.0_r8/pm(i,j+1))
              cff=0.5_r8*phat*dxf
              IF (pmask_act(i,j+1).eq.1.0_r8.and.                       &
     &            i.ge.IstrP.and.i.le.IendP.and.                        &
     &            j+1.ge.JstrP.and.j+1.le.JendP) THEN
                tdy(i,j+1)=tdy(i,j+1)+cff
              END IF
              IF (pmask_act(i+1,j+1).eq.1.0_r8.and.                     &
     &            i+1.ge.IstrP.and.i+1.le.IendP.and.                    &
     &            j+1.ge.JstrP.and.j+1.le.JendP) THEN
                tdy(i+1,j+1)=tdy(i+1,j+1)+cff
              END IF
            END IF
!
            IF (front_s) THEN
              dxf=0.5_r8*(1.0_r8/pm(i,j-1)+1.0_r8/pm(i,j))
              cff=0.5_r8*phat*dxf
              IF (pmask_act(i,j).eq.1.0_r8.and.                         &
     &            i.ge.IstrP.and.i.le.IendP.and.                        &
     &            j.ge.JstrP.and.j.le.JendP) THEN
                tdy(i,j)=tdy(i,j)-cff
              END IF
              IF (pmask_act(i+1,j).eq.1.0_r8.and.                       &
     &            i+1.ge.IstrP.and.i+1.le.IendP.and.                    &
     &            j.ge.JstrP.and.j.le.JendP) THEN
                tdy(i+1,j)=tdy(i+1,j)-cff
              END IF
            END IF
!
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE ssa_front_load
!
!***********************************************************************
      SUBROUTINE ssa_front_faces (ng, tile, i, j,                       &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Istr, Iend, Jstr, Jend,               &
     &                            rmask_ist,                            &
     &                            front_e, front_w, front_n, front_s)
!***********************************************************************
!
!  Front-face predicate shared by ssa_front_load and ssa_land_margin:
!  the east/west/north/south face of cell (i,j) is a front face if
!  the corresponding neighbour cell is ice-free.  A free (gradient)
!  lateral boundary counts as ice-free (the ice continues into the
!  open domain);  closed and clamped edges are walls and never
!  fronts.  Domain-edge ghost cells are outside the caller cell
!  loops (IstrC/JstrC skip them), so only the override reaches
!  beyond the physical domain.
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, i, j
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Istr, Iend, Jstr, Jend
!
      logical, intent(out) :: front_e, front_w, front_n, front_s
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: rmask_ist(LBi:,LBj:)
# else
      real(r8), intent(in) :: rmask_ist(LBi:UBi,LBj:UBj)
# endif
!
      front_e=(rmask_ist(i+1,j  ).eq.0.0_r8)
      front_w=(rmask_ist(i-1,j  ).eq.0.0_r8)
      front_n=(rmask_ist(i  ,j+1).eq.0.0_r8)
      front_s=(rmask_ist(i  ,j-1).eq.0.0_r8)
!
!  At a physical domain edge the neighbour cell does not exist and
!  the ghost mask is not consulted:  the edge status comes from the
!  LBC alone.  A free (gradient) edge is open ocean (a front);  a
!  closed or clamped edge is a wall (never a front).  Overriding in
!  BOTH directions matters:  leaving the zero-ghost mask value in
!  place at closed/clamped edges declared every ice cell at the
!  wall a front -- the spurious wall load landed on the Dirichlet
!  boundary nodes (inert for the velocity solution, but polluting
!  the residual norms and the is_taudx/is_taudy diagnostics), and
!  ssa_land_margin flagged a spurious land margin (driving-stress
!  taper) along the wall.
!
      IF (.not.EWperiodic(ng)) THEN
        IF (i.eq.Iend) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            front_e=LBC(ieast,ibISF(isIuis),ng)%gradient
          END IF
        END IF
        IF (i.eq.Istr) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            front_w=LBC(iwest,ibISF(isIuis),ng)%gradient
          END IF
        END IF
      END IF
      IF (.not.NSperiodic(ng)) THEN
        IF (j.eq.Jend) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            front_n=LBC(inorth,ibISF(isIvis),ng)%gradient
          END IF
        END IF
        IF (j.eq.Jstr) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            front_s=LBC(isouth,ibISF(isIvis),ng)%gradient
          END IF
        END IF
      END IF
!
      RETURN
      END SUBROUTINE ssa_front_faces
!
!***********************************************************************
      SUBROUTINE ssa_land_margin (ng, tile, model,                      &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           IstrC, IendC, JstrC, JendC,            &
     &                           IstrP, IendP, JstrP, JendP,            &
     &                           h, zeta_avg, rmask_ist, hi,            &
     &                           pmarg)
!***********************************************************************
!
!  Land-margin indicator for the driving-stress taper (plan option
!  A).  A RHO cell is a land-margin FRONT cell when it holds
!  grounded ice with its bed at or above sea level (D_avg <= 0,
!  same window-averaged criterion as the basal drag and
!  ssa_front_load) AND at least one of its faces is a front face
!  (ssa_front_faces).  pmarg = 1 at PSI nodes one of whose 4 RHO
!  cells is such a cell;  the taper distance in ssa_driving_stress
!  is measured to these nodes ONLY -- marine fronts (floating or
!  grounded in water) carry the full front-pressure load and keep
!  the full driving stress.
!
!  The cell loop runs over the assembly ring (IstrC..IendC) so the
!  seam nodes of the neighbouring tile are classified identically
!  from both sides;  pmarg is then halo-exchanged for the 8-neighbour
!  taper scan at owned nodes.
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Istr, Iend, Jstr, Jend
      integer, intent(in) :: IstrC, IendC, JstrC, JendC
      integer, intent(in) :: IstrP, IendP, JstrP, JendP
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: h(LBi:,LBj:)
      real(r8), intent(in)  :: zeta_avg(LBi:,LBj:)
      real(r8), intent(in)  :: rmask_ist(LBi:,LBj:)
      real(r8), intent(in)  :: hi(LBi:,LBj:)
      real(r8), intent(out) :: pmarg(LBi:,LBj:)
# else
      real(r8), intent(in)  :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: zeta_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: rmask_ist(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: hi(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: pmarg(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
!
      logical :: front_e, front_w, front_n, front_s
!
      real(r8) :: D_avg
!
      real(r8), dimension(LBi:UBi,LBj:UBj) :: lcell
!
!-----------------------------------------------------------------------
!  Classify the land-margin front cells and project to PSI nodes.
!-----------------------------------------------------------------------
!
      DO j=LBj,UBj
        DO i=LBi,UBi
          lcell(i,j)=0.0_r8
        END DO
      END DO
!
      DO j=JstrC,JendC
        DO i=IstrC,IendC
          IF (rmask_ist(i,j).eq.1.0_r8) THEN
            D_avg=h(i,j)+zeta_avg(i,j)-Dcrit(ng)
!
!  D_avg <= 0 (bed above the window-averaged sea level) implies
!  rho_i*H > rho0*D_avg for every H > 0 (rmask gate above), so the
!  single comparison suffices.
!
            IF (D_avg.le.0.0_r8) THEN
              CALL ssa_front_faces (ng, tile, i, j,                     &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Istr, Iend, Jstr, Jend,             &
     &                              rmask_ist,                          &
     &                              front_e, front_w,                   &
     &                              front_n, front_s)
              IF (front_e.or.front_w.or.front_n.or.front_s) THEN
                lcell(i,j)=1.0_r8
              END IF
            END IF
          END IF
        END DO
      END DO
!
!  Zero the FULL array first:  the taper scan in ssa_driving_stress
!  reads pmarg one node beyond the owned range, and on single-tile
!  runs that halo is never exchanged (an uninitialized automatic
!  array made the taper nondeterministic there).
!
      DO j=LBj,UBj
        DO i=LBi,UBi
          pmarg(i,j)=0.0_r8
        END DO
      END DO
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          pmarg(i,j)=MAX(lcell(i-1,j-1),lcell(i  ,j-1),                 &
     &                   lcell(i-1,j  ),lcell(i  ,j  ))
        END DO
      END DO
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmarg)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pmarg)
# endif
!
      RETURN
      END SUBROUTINE ssa_land_margin
!
!***********************************************************************
      SUBROUTINE ssa_hardness (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         hi, ti,                                  &
     &                         B_hard)
!***********************************************************************
!
!  Compute the Glen hardness field B = A^(-1/n) at RHO points.
!  B depends only on the ice temperature and the ice thickness,
!  both frozen during a solve (lsold values), so this routine is
!  called once per solve (Stage 4) and ssa_viscosity only combines
!  the cached field with the current strain-rate invariant.
!
!  The Glen rate factor A is either the constant A_glen
!  (isothermal override, used for ISOMIP-style tests) or the
!  Paterson (1994) temperature-dependent formula with
!  pressure-adjusted temperature
!  T* = T + 273.15 + 8.66e-4 * 0.5 * H  (Bueler & Brown 2009):
!
!    A = 3.61e-13 * exp(-60000/(R T*)),   T* <= 263.15 K
!    A = 1.73e3   * exp(-139000/(R T*)),  T* >  263.15 K
!
!  Cells without ice get B = 0 (their viscosity is zeroed by
!  ssa_viscosity anyway).  No halo exchange:  B_hard is consumed
!  by ssa_viscosity over the same interior range only.
!
!  Reference:  Paterson (1994) Table 5.2, Bueler & Brown (2009).
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Istr, Iend, Jstr, Jend
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: hi(LBi:,LBj:)
      real(r8), intent(in)  :: ti(LBi:,LBj:)
      real(r8), intent(out) :: B_hard(LBi:,LBj:)
# else
      real(r8), intent(in)  :: hi(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: ti(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: B_hard(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
!
      real(r8) :: T_star, RT_star, A_pat
!
!-----------------------------------------------------------------------
!  Assemble the hardness field.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (hi(i,j).gt.0.0_r8) THEN
            IF (A_glen(ng).gt.0.0_r8) THEN
              B_hard(i,j)=A_glen(ng)**(-1.0_r8/n_glen(ng))
            ELSE
              T_star=ti(i,j)+273.15_r8+8.66E-4_r8*0.5_r8*hi(i,j)
              RT_star=8.314_r8*T_star
              IF (T_star.le.263.15_r8) THEN
                A_pat=3.61E-13_r8*EXP(-6.0E4_r8/RT_star)
              ELSE
                A_pat=1.73E3_r8*EXP(-13.9E4_r8/RT_star)
              END IF
              B_hard(i,j)=A_pat**(-1.0_r8/n_glen(ng))
            END IF
          ELSE
            B_hard(i,j)=0.0_r8
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE ssa_hardness
!
!***********************************************************************
      SUBROUTINE ssa_viscosity (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Istr, Iend, Jstr, Jend,                 &
     &                          pm, pn,                                 &
     &                          hi, ui, vi, B_hard,                     &
     &                          nuh)
!***********************************************************************
!
!  Compute the depth-integrated effective viscosity at RHO points:
!
!    gamma = eps_xx^2 + eps_yy^2 + eps_xx*eps_yy + eps_xy^2 +
!            eps_glen_min^2  (squared second invariant)
!    nuH   = 0.5 * B * gamma^((1-n)/(2n)) * H
!
!  where B = A^(-1/n) is the Glen hardness, PRECOMPUTED by the
!  ssa_hardness (isf_ssa_tile Stage 4):  B depends only on the ice
!  temperature and the lsold thickness, both frozen during the
!  Picard iteration, so the Paterson (1994) rate-factor evaluation
!  runs once per solve instead of once per viscosity call.
!
!  The exponent (1-n)/(2n) acts on the SQUARED invariant gamma:
!  this is the PISM convention (FlowLaw::effective_viscosity,
!  FlowLaw.cc) and is equivalent to the STREAMICE form
!  nuH = 0.5*B*eps_e^((1-n)/n)*H with eps_e = sqrt(gamma)
!  (streamice_visc_beta.F);  both reduce to the MacAyeal (1989)
!  Glen-law effective viscosity nu = 0.5*B*gamma^((1-n)/(2n)).
!
!  Strain rates are centered differences of the 4 PSI-corner
!  velocities of each RHO cell.  PISM robustness measures:
!    - strength extension for thin ice (H_eff = max(H, min_thick-
!      ness)), so margin cells keep membrane stiffness
!    - additive nuH regularization (nuH += nuH_reg)
!
!  Reference:  Paterson (1994) Table 5.2, Bueler & Brown (2009),
!  MITgcm STREAMICE streamice_visc_beta.F, PISM compute_nuH_everywhere.
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Istr, Iend, Jstr, Jend
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: pm(LBi:,LBj:)
      real(r8), intent(in)  :: pn(LBi:,LBj:)
      real(r8), intent(in)  :: hi(LBi:,LBj:)
      real(r8), intent(in)  :: ui(LBi:,LBj:)
      real(r8), intent(in)  :: vi(LBi:,LBj:)
      real(r8), intent(in)  :: B_hard(LBi:,LBj:)
      real(r8), intent(out) :: nuh(LBi:,LBj:)
# else
      real(r8), intent(in)  :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: hi(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: ui(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: vi(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: B_hard(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: nuh(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
!
      real(r8) :: cff, eps_xx, eps_yy, eps_xy
!
!-----------------------------------------------------------------------
!  Compute viscosity at the RHO cells with ice (H > 0).  Cell (i,j)
!  has 4 PSI corner nodes:  (i,j), (i+1,j), (i,j+1), (i+1,j+1).
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          nuh(i,j)=0.0_r8
          IF (hi(i,j).gt.0.0_r8) THEN
!
!  Centered strain rates from the 4 PSI-corner velocities:
!  nodes 1=(i,j), 2=(i+1,j), 3=(i,j+1), 4=(i+1,j+1).
!
            eps_xx=0.5_r8*((ui(i+1,j  )-ui(i,  j  ))+                   &
     &                     (ui(i+1,j+1)-ui(i,  j+1)))*pm(i,j)
            eps_yy=0.5_r8*((vi(i,  j+1)-vi(i,  j  ))+                   &
     &                     (vi(i+1,j+1)-vi(i+1,j  )))*pn(i,j)
            eps_xy=0.5_r8*(0.5_r8*((vi(i+1,j  )-vi(i,  j  ))+           &
     &                             (vi(i+1,j+1)-vi(i,  j+1)))*pm(i,j)+  &
     &                     0.5_r8*((ui(i,  j+1)-ui(i,  j  ))+           &
     &                             (ui(i+1,j+1)-ui(i+1,j  )))*pn(i,j))
!
!  Glen-law effective viscosity from the squared second invariant,
!  with strength extension for thin ice and additive regularization.
!  The bracket is strictly positive through the eps_glen_min^2 term.
!
            cff=0.5_r8*B_hard(i,j)*                                     &
     &          (eps_xx**2+eps_yy**2+eps_xx*eps_yy+eps_xy**2+           &
     &           eps_glen_min(ng)**2)**                                 &
     &          ((1.0_r8-n_glen(ng))/(2.0_r8*n_glen(ng)))
!
!  Strength extension (PISM):  thin ice is treated with the minimum
!  thickness in the H factor so margin cells retain membrane
!  stiffness (H < min_thickness would give nuH -> 0).
!
            nuh(i,j)=cff*MAX(hi(i,j),min_thickness(ng))
            nuh(i,j)=nuh(i,j)+nuH_reg(ng)
          END IF
        END DO
      END DO
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          nuh)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    nuh)
# endif
!
      RETURN
      END SUBROUTINE ssa_viscosity
!
!***********************************************************************
      SUBROUTINE ssa_jacobi_diag (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IstrC, IendC, JstrC, JendC,           &
     &                            IstrP, IendP, JstrP, JendP,           &
     &                            DPhi, phi_tab,                        &
     &                            pm, pn,                               &
     &                            pmask_act,                            &
     &                            nuH, beta,                            &
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
!  halo ring (IstrC..IendC);  scatters land only on owned PSI
!  nodes (ownership filter).
!
!  Reference:  MITgcm STREAMICE CG_ADIAG.
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IstrC, IendC, JstrC, JendC
      integer, intent(in) :: IstrP, IendP, JstrP, JendP
!
      real(r8), intent(in)  :: DPhi(4,4,2)
      real(r8), intent(in)  :: phi_tab(4,4)
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: pm(LBi:,LBj:)
      real(r8), intent(in)  :: pn(LBi:,LBj:)
      real(r8), intent(in)  :: pmask_act(LBi:,LBj:)
      real(r8), intent(in)  :: nuH(LBi:,LBj:)
      real(r8), intent(in)  :: beta(LBi:,LBj:)
      real(r8), intent(out) :: diag_u(LBi:,LBj:)
      real(r8), intent(out) :: diag_v(LBi:,LBj:)
# else
      real(r8), intent(in)  :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pmask_act(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: nuH(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: beta(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: diag_u(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: diag_v(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, i1, j1, nq, k_node
!
      real(r8) :: dx_rho, dy_rho, jac
      real(r8) :: dphi_x, dphi_y, phv
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
      DO j=JstrC,JendC
        DO i=IstrC,IendC
          IF (nuH(i,j).gt.0.0_r8) THEN
            dx_rho=1.0_r8/pm(i,j)
            dy_rho=1.0_r8/pn(i,j)
            jac=1.0_r8/(pm(i,j)*pn(i,j))
!
            DO nq=1,4
              DO k_node=1,4
                IF (k_node.eq.1) THEN
                  i1=i
                  j1=j
                ELSE IF (k_node.eq.2) THEN
                  i1=i+1
                  j1=j
                ELSE IF (k_node.eq.3) THEN
                  i1=i
                  j1=j+1
                ELSE
                  i1=i+1
                  j1=j+1
                END IF
!
!  Ownership filter:  keep only contributions landing on PSI nodes
!  owned by this tile (seam nodes are fed by the halo ring cells
!  integrated here, and by the own cells integrated on the
!  neighbouring tile).
!
                IF (pmask_act(i1,j1).eq.1.0_r8.and.                     &
     &              i1.ge.IstrP.and.i1.le.IendP.and.                    &
     &              j1.ge.JstrP.and.j1.le.JendP) THEN
                  dphi_x=DPhi(k_node,nq,1)*2.0_r8/dx_rho
                  dphi_y=DPhi(k_node,nq,2)*2.0_r8/dy_rho
                  phv=phi_tab(k_node,nq)
!
                  diag_u(i1,j1)=diag_u(i1,j1)+                          &
     &                          0.25_r8*jac*nuH(i,j)*                   &
     &                          (4.0_r8*dphi_x*dphi_x+dphi_y*dphi_y)+   &
     &                          0.25_r8*jac*beta(i,j)*phv**2
                  diag_v(i1,j1)=diag_v(i1,j1)+                          &
     &                          0.25_r8*jac*nuH(i,j)*                   &
     &                          (4.0_r8*dphi_y*dphi_y+dphi_x*dphi_x)+   &
     &                          0.25_r8*jac*beta(i,j)*phv**2
                END IF
              END DO
            END DO
          END IF
        END DO
      END DO
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, diag_u)
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, diag_v)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    diag_u, diag_v)
# endif
!
      RETURN
      END SUBROUTINE ssa_jacobi_diag
!
!***********************************************************************
      SUBROUTINE ssa_matvec (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IstrC, IendC, JstrC, JendC,                &
     &                       IstrP, IendP, JstrP, JendP,                &
     &                       DPhi, phi_tab,                             &
     &                       pm, pn,                                    &
     &                       pmask_act,                                 &
     &                       nuH, beta,                                 &
     &                       xu, xv,                                    &
     &                       Au, Av)
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
!  The INPUT halos are exchanged inside this routine:  the ring
!  cells' assembly reads PSI corners owned by the neighbouring
!  tiles, so (xu,xv) carry their halo values on entry (hence the
!  intent(inout) declaration).  No output exchange:  every
!  consumer of (Au,Av) reads owned PSI nodes only.
!
!  Reference:  MITgcm STREAMICE streamice_cg_functions.F (CG_ACTION).
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IstrC, IendC, JstrC, JendC
      integer, intent(in) :: IstrP, IendP, JstrP, JendP
!
      real(r8), intent(in)    :: DPhi(4,4,2)
      real(r8), intent(in)    :: phi_tab(4,4)
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: pm(LBi:,LBj:)
      real(r8), intent(in)  :: pn(LBi:,LBj:)
      real(r8), intent(in)  :: pmask_act(LBi:,LBj:)
      real(r8), intent(in)  :: nuH(LBi:,LBj:)
      real(r8), intent(in)  :: beta(LBi:,LBj:)
      real(r8), intent(inout) :: xu(LBi:,LBj:)
      real(r8), intent(inout) :: xv(LBi:,LBj:)
      real(r8), intent(out) :: Au(LBi:,LBj:)
      real(r8), intent(out) :: Av(LBi:,LBj:)
# else
      real(r8), intent(in)  :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pmask_act(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: nuH(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: beta(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: xu(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: xv(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Au(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Av(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, i1, j1, nq, k_node
!
      real(r8) :: dx_rho, dy_rho, jac
      real(r8) :: dphi_x, dphi_y, phv
      real(r8) :: uq, vq, exx, eyy, exy
!
!-----------------------------------------------------------------------
!  Exchange the INPUT halos before assembly:  the ring cells
!  (IstrC/JstrC) integrate elements whose PSI corners include nodes
!  owned by the neighbouring tiles, so the quadrature-point
!  interpolation must read the neighbours' current values of the
!  input pair.  Only the halo points are written; the interior is
!  returned unchanged (hence intent(inout) for xu/xv).
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
      DO j=JstrC,JendC
        DO i=IstrC,IendC
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
                  i1=i
                  j1=j
                ELSE IF (k_node.eq.2) THEN
                  i1=i+1
                  j1=j
                ELSE IF (k_node.eq.3) THEN
                  i1=i
                  j1=j+1
                ELSE
                  i1=i+1
                  j1=j+1
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
                IF (pmask_act(i1,j1).eq.1.0_r8.and.                     &
     &              i1.ge.IstrP.and.i1.le.IendP.and.                    &
     &              j1.ge.JstrP.and.j1.le.JendP) THEN
                  dphi_x=DPhi(k_node,nq,1)*2.0_r8/dx_rho
                  dphi_y=DPhi(k_node,nq,2)*2.0_r8/dy_rho
                  phv=phi_tab(k_node,nq)
!
                  Au(i1,j1)=Au(i1,j1)+                                  &
     &              0.25_r8*jac*nuH(i,j)*                               &
     &              (dphi_x*(4.0_r8*exx+2.0_r8*eyy)+                    &
     &               dphi_y*(2.0_r8*exy))+                              &
     &              0.25_r8*jac*beta(i,j)*phv*uq
!
                  Av(i1,j1)=Av(i1,j1)+                                  &
     &              0.25_r8*jac*nuH(i,j)*                               &
     &              (dphi_y*(4.0_r8*eyy+2.0_r8*exx)+                    &
     &               dphi_x*(2.0_r8*exy))+                              &
     &              0.25_r8*jac*beta(i,j)*phv*vq
                END IF
              END DO
            END DO
          END IF
        END DO
      END DO
!
!  No output exchange:  every consumer of (Au,Av) -- the CG inner
!  products, the residual update and the convergence norms -- reads
!  them at owned PSI nodes only, and each owned node receives the
!  complete contribution of all its elements through the ring
!  assembly above (seam elements are integrated by the tile that
!  owns the node).
!
      RETURN
      END SUBROUTINE ssa_matvec
#endif
      END MODULE isf_solver_mod
