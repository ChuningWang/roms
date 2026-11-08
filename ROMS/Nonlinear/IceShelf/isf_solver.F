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
# ifdef ICESHEET_PETSC
#  include <petsc/finclude/petscksp.h>
      USE petscksp
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
      integer :: i, j, cg_iter, nl_iter, ic
      logical :: u_dW, u_dE, u_dS, u_dN, v_dW, v_dE, v_dS, v_dN
!
      logical :: converged
!
      real(r8) :: cff, cff1, u_slide, beta_avg, jac_avg
      real(r8) :: u_prev, u_allow
      integer :: nb
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
!  FD solver (n_ssa_fd > 0) work arrays:  the FD discretization
!  (PISM SSAFD recipe) defines the equations at the RHO cell
!  centers, so the solver state lives on RHO points;  the staggered
!  viscosities live on the U/V faces (nuhx = x-face between
!  (i,j)-(i+1,j), nuhy = y-face between (i,j)-(i,j+1)) and the
!  driving forces on the cell centers (fdx/fdy).
!
      real(r8), dimension(LBi:UBi,LBj:UBj) :: nuhx
      real(r8), dimension(LBi:UBi,LBj:UBj) :: nuhy
      real(r8), dimension(LBi:UBi,LBj:UBj) :: fdx
      real(r8), dimension(LBi:UBi,LBj:UBj) :: fdy
      real(r8), dimension(LBi:UBi,LBj:UBj) :: urho
      real(r8), dimension(LBi:UBi,LBj:UBj) :: vrho
      real(r8), dimension(LBi:UBi,LBj:UBj) :: aurho
      real(r8), dimension(LBi:UBi,LBj:UBj) :: avrho
      real(r8), dimension(LBi:UBi,LBj:UBj) :: beta_rho
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
!  Static per-cell Q1 metric factors for the assembly kernels (filled
!  once per solve from the fixed grid metrics pm/pn; see the fill loop
!  below):  fac_dx = 2*pm = 2/dx (the Q1 reference cell spans [-1,1]),
!  fac_dy = 2*pn = 2/dy, fac_jac = 1/(pm*pn) (cell area).  Passing
!  these instead of pm/pn removes ~35 divisions per cell from every
!  ssa_matvec / ssa_jacobi_diag call -- matvec runs once per CG
!  iteration, the hottest kernel of the solve.
!
      real(r8), dimension(LBi:UBi,LBj:UBj) :: fac_dx, fac_dy, fac_jac
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
!  Carry the refreeze reservoirs Vi/Ti across the time-level swap
!  (review finding M-6):  isf_vbc accumulates them in the new level
!  every step, but the driver swaps the layer roles just before this
!  solve without copying them, so under N_ssa > 1 the reservoir chain
!  would resume from one-to-two windows ago.  Interior ranges only:
!  isf_vbc consumes interior values and the reservoirs are not halo-
!  exchanged (output masks the halo).
!
      IF (Lisftrc(ng)) THEN
        DO j=Jstr,Jend
          DO i=Istr,Iend
            ISF(ng) % Vi(i,j,lsnew)=ISF(ng) % Vi(i,j,lsold)
            DO ic=1,NTISF(ng)
              ISF(ng) % Ti(i,j,lsnew,ic)=ISF(ng) % Ti(i,j,lsold,ic)
            END DO
          END DO
        END DO
      END IF
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
!  Precompute the static per-cell metric factors over the assembly
!  ring (the exact range the assembly kernels evaluate them on; the
!  pm/pn halo values feeding the seam ring are exchanged grid fields).
!  Replacing DPhi*2/dx_rho (two roundings and an FDIV per quadrature
!  point) by DPhi*fac_dx changes results only at roundoff level.
!
      DO j=JstrC,JendC
        DO i=IstrC,IendC
          fac_dx(i,j)=2.0_r8*pm(i,j)
          fac_dy(i,j)=2.0_r8*pn(i,j)
          fac_jac(i,j)=1.0_r8/(pm(i,j)*pn(i,j))
        END DO
      END DO
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
      CALL ssa_land_margin (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Istr, Iend, Jstr, Jend,                     &
     &                      IstrC, IendC, JstrC, JendC,                 &
     &                      IstrP, IendP, JstrP, JendP,                 &
     &                      h, zeta_avg, rmask_ist,                     &
     &                      Si(:,:,lsold,isIhis), pmarg)
!
      IF (n_ds_fem(ng).gt.0) THEN
!
!  B5.8:  PISM-FEM style driving stress (smooth bed/H gradients +
!  pointwise flotation branch at quadrature points;  replaces the
!  node-based FD gradient of the branch-mixed surface field).
!
        CALL ssa_ds_fem (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Istr, Iend, Jstr, Jend,                       &
     &                    IstrP, IendP, JstrP, JendP,                   &
     &                    DPhi, phi_tab,                                &
     &                    pm, pn, h, zeta_avg,                          &
     &                    rmask_ist, pmask_act,                         &
     &                    Si(:,:,lsold,isIhis),                         &
     &                    Fi(:,:,ifItdx),                                &
     &                    Fi(:,:,ifItdy))
      ELSE
        CALL ssa_driving_stress (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           pm, pn, h, zeta_avg,                   &
     &                           rmask_ist, pmask_act, pmarg,           &
     &                           Si(:,:,lsold,isIhis),                  &
     &                           Fi(:,:,ifIzsf),                        &
     &                           Fi(:,:,ifItdx),                        &
     &                           Fi(:,:,ifItdy))
      END IF
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
!  FD path (n_ssa_fd > 0):  replace the nodal RHS with the mapped
!  cell-centered force of ssa_fd_driving (branch-aware gradient +
!  cfbc delta_p).  Linear in the state, so it is assembled once
!  here;  the is_taudx/is_taudy diagnostics then carry the FD node
!  forces in this mode.
!
        IF (n_ssa_fd(ng).gt.0) THEN
!
!  Map the current PSI state to cell centers (FD state mask).  The
!  4-corner average beta_rho is the cell-centered drag consumed by
!  the FD drivers (PISM samples tauc at cell centers);  before the
!  2026-10-31 rewiring it was computed but never used -- the FD
!  operators sampled the PSI field at cell centers instead, a
!  half-cell offset of the drag.
!
          DO j=LBj,UBj
            DO i=LBi,UBi
              IF (Si(i,j,lsold,isIhis).gt.0.0_r8) THEN
                urho(i,j)=0.25_r8*(Si(i,j,lsnew,isIuis)+                 &
     &                             Si(i+1,j,lsnew,isIuis)+               &
     &                             Si(i,j+1,lsnew,isIuis)+               &
     &                             Si(i+1,j+1,lsnew,isIuis))
                vrho(i,j)=0.25_r8*(Si(i,j,lsnew,isIvis)+                 &
     &                             Si(i+1,j,lsnew,isIvis)+               &
     &                             Si(i,j+1,lsnew,isIvis)+               &
     &                             Si(i+1,j+1,lsnew,isIvis))
                beta_rho(i,j)=0.25_r8*(Fi(i,j,ifIbta)+                   &
     &                                 Fi(i+1,j,ifIbta)+                 &
     &                                 Fi(i,j+1,ifIbta)+                 &
     &                                 Fi(i+1,j+1,ifIbta))
              ELSE
                urho(i,j)=0.0_r8
                vrho(i,j)=0.0_r8
                beta_rho(i,j)=0.0_r8
              END IF
            END DO
          END DO
!
          CALL ssa_fd_driving (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         pm, pn, h, zeta_avg,                     &
     &                         Si(:,:,lsold,isIhis),                    &
     &                         fdx, fdy)
!
!  Transpose map to the PSI nodes (consumed by the shared PCG).
!
          DO j=JstrP,JendP
            DO i=IstrP,IendP
              nb=0
              cff=0.0_r8
              cff1=0.0_r8
              IF (i-1.ge.Istr.and.j-1.ge.Jstr) THEN
                cff=cff+fdx(i-1,j-1)
                cff1=cff1+fdy(i-1,j-1)
                nb=nb+1
              END IF
              IF (i.le.Iend.and.j-1.ge.Jstr) THEN
                cff=cff+fdx(i,j-1)
                cff1=cff1+fdy(i,j-1)
                nb=nb+1
              END IF
              IF (i-1.ge.Istr.and.j.le.Jend) THEN
                cff=cff+fdx(i-1,j)
                cff1=cff1+fdy(i-1,j)
                nb=nb+1
              END IF
              IF (i.le.Iend.and.j.le.Jend) THEN
                cff=cff+fdx(i,j)
                cff1=cff1+fdy(i,j)
                nb=nb+1
              END IF
              IF (nb.gt.0) THEN
                Fi(i,j,ifItdx)=cff/REAL(nb,r8)
                Fi(i,j,ifItdy)=cff1/REAL(nb,r8)
              ELSE
                Fi(i,j,ifItdx)=0.0_r8
                Fi(i,j,ifItdy)=0.0_r8
              END IF
            END DO
          END DO
#  ifdef ISF_DEBUG
          cff=0.0_r8
          DO j=Jstr,Jend
            DO i=Istr,Iend
              cff=MAX(cff,ABS(Fi(i,j,ifItdx)))
            END DO
          END DO
          WRITE (*,*) 'FD-RHS probe: max|Fi(ifItdx)| =', cff,           &
     &                ' max|fdx(cell)| probe:'
          cff=0.0_r8
          DO j=MAX(LBj,Jstr-1),MIN(UBj,Jend)
            DO i=MAX(LBi,Istr-1),MIN(UBi,Iend)
              cff=MAX(cff,ABS(fdx(i,j)))
            END DO
          END DO
          WRITE (*,*) 'FD-RHS probe: max|fdx| =', cff
#  endif
        END IF
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
      IF (n_ssa_fd(ng).gt.0) THEN
!
!-----------------------------------------------------------------------
!  FD path (PISM SSAFD recipe):  the Picard/PCG runs entirely in the
!  RHO space (ssa_fd_picard;  the averaged-projection variant
!  diverged through the checkerboard null space of the PSI-space
!  operator).  The cell-centered RHS lives in fdx/fdy from Stage 2,
!  and the converged cell state maps back with the transpose
!  average.
!-----------------------------------------------------------------------
!
        DO j=LBj,UBj
          DO i=LBi,UBi
            IF (Si(i,j,lsold,isIhis).gt.0.0_r8) THEN
              urho(i,j)=0.25_r8*(Si(i,j,lsnew,isIuis)+                   &
     &                           Si(i+1,j,lsnew,isIuis)+                 &
     &                           Si(i,j+1,lsnew,isIuis)+                 &
     &                           Si(i+1,j+1,lsnew,isIuis))
              vrho(i,j)=0.25_r8*(Si(i,j,lsnew,isIvis)+                   &
     &                           Si(i+1,j,lsnew,isIvis)+                 &
     &                           Si(i,j+1,lsnew,isIvis)+                 &
     &                           Si(i+1,j+1,lsnew,isIvis))
            ELSE
              urho(i,j)=0.0_r8
              vrho(i,j)=0.0_r8
            END IF
          END DO
        END DO
!
        IF (fd_ksp(ng).gt.0) THEN
# ifdef ICESHEET_PETSC
          CALL ssa_fd_picard_ksp (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Istr, Iend, Jstr, Jend,               &
     &                            pm, pn,                               &
     &                            Si(:,:,lsold,isIhis),                 &
     &                            nuhx, nuhy, beta_rho,           &
     &                            fdx, fdy, B_hard,                     &
     &                            urho, vrho,                           &
     &                            Ru, Rv, Zu, Zv, Du, Dv,               &
     &                            Au, Av, diag_u, diag_v)
# else
          IF (Master) WRITE(stdout,'(a)')                               &
     &    ' SSA-FD WARNING: fd_ksp>0 but ICESHEET_PETSC not compiled;'// &
     &    ' falling back to the internal PCG.'
          CALL ssa_fd_picard (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Istr, Iend, Jstr, Jend,                   &
     &                        pm, pn,                                   &
     &                        Si(:,:,lsold,isIhis),                     &
     &                        nuhx, nuhy, beta_rho,               &
     &                        fdx, fdy, B_hard,                         &
     &                        urho, vrho,                               &
     &                        Ru, Rv, Zu, Zv, Du, Dv,                   &
     &                        Au, Av, diag_u, diag_v)
# endif
        ELSE
          CALL ssa_fd_picard (ng, tile, model,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Istr, Iend, Jstr, Jend,                     &
     &                      pm, pn,                                     &
     &                      Si(:,:,lsold,isIhis),                       &
     &                      nuhx, nuhy, beta_rho,                 &
     &                      fdx, fdy, B_hard,                           &
     &                      urho, vrho,                                 &
     &                      Ru, Rv, Zu, Zv, Du, Dv,                     &
     &                      Au, Av, diag_u, diag_v)
        END IF
!
!  Map the converged cell state back to the PSI nodes (transpose
!  average;  boundary nodes use the available cells).
!
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            nb=0
            cff=0.0_r8
            cff1=0.0_r8
            IF (i-1.ge.Istr.and.j-1.ge.Jstr) THEN
              cff=cff+urho(i-1,j-1)
              cff1=cff1+vrho(i-1,j-1)
              nb=nb+1
            END IF
            IF (i.le.Iend.and.j-1.ge.Jstr) THEN
              cff=cff+urho(i,j-1)
              cff1=cff1+vrho(i,j-1)
              nb=nb+1
            END IF
            IF (i-1.ge.Istr.and.j.le.Jend) THEN
              cff=cff+urho(i-1,j)
              cff1=cff1+vrho(i-1,j)
              nb=nb+1
            END IF
            IF (i.le.Iend.and.j.le.Jend) THEN
              cff=cff+urho(i,j)
              cff1=cff1+vrho(i,j)
              nb=nb+1
            END IF
            IF (nb.gt.0) THEN
              Si(i,j,lsnew,isIuis)=cff/REAL(nb,r8)
              Si(i,j,lsnew,isIvis)=cff1/REAL(nb,r8)
            END IF
          END DO
        END DO
!
!  FD diagnostics:  the cell-centred effective viscosity for the
!  is_nuH output (the mean of the cell's four face values from the
!  last Picard iteration;  the FE path fills this in ssa_viscosity).
!
        DO j=Jstr,Iend
          DO i=Istr,Iend
            IF (Si(i,j,lsold,isIhis).gt.0.0_r8) THEN
              Fi(i,j,ifInuh)=0.25_r8*(nuhx(i-1,j)+nuhx(i,j)+            &
     &                                nuhy(i,j-1)+nuhy(i,j))
            ELSE
              Fi(i,j,ifInuh)=0.0_r8
            END IF
          END DO
        END DO
!
      ELSE
!
      DO nl_iter=1,nl_max_iter(ng)
!
!  Effective viscosity at RHO points from the current velocity
!  estimate (Glen law with the Stage-4 cached hardness field,
!  thin-ice strength extension and additive regularization).
!
        CALL ssa_viscosity (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Istr, Iend, Jstr, Jend,                   &
     &                        pm, pn,                                   &
     &                        Si(:,:,lsold,isIhis),                     &
     &                        Si(:,:,lsnew,isIuis),                     &
     &                        Si(:,:,lsnew,isIvis),                     &
     &                        B_hard,                                   &
     &                        Fi(:,:,ifInuh))
!
          CALL ssa_jacobi_diag (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IstrC, IendC, JstrC, JendC,             &
     &                          IstrP, IendP, JstrP, JendP,             &
     &                          DPhi, phi_tab,                          &
     &                          fac_dx, fac_dy, fac_jac,                &
     &                          pmask_act,                              &
     &                          Fi(:,:,ifInuh),                         &
     &                          Fi(:,:,ifIbta),                         &
     &                          diag_u, diag_v)
!
          CALL ssa_matvec (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Istr, Iend, Jstr, Jend,                    &
     &                     IstrC, IendC, JstrC, JendC,                  &
     &                     IstrP, IendP, JstrP, JendP,                  &
     &                     DPhi, phi_tab,                               &
     &                     fac_dx, fac_dy, fac_jac,                     &
     &                     pmask_act, rmask_ist,                        &
     &                     Fi(:,:,ifInuh),                              &
     &                     Fi(:,:,ifIbta),                              &
     &                     Si(:,:,lsnew,isIuis),                        &
     &                     Si(:,:,lsnew,isIvis),                        &
     &                     Au, Av)
!
!  Initial residual R = b - A(u) of the Picard linearization (the
!  c942bde1 Stage-6 restructure accidentally dropped this block:
!  with Stage-5-zeroed Ru the first dot_zr is zero, the CG loop
!  exits immediately and the FE solve degenerates to echoing the
!  input velocity -- regression found and fixed 2026-10-28).
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
        dot_zr_old=dot_zr
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
     &                     Istr, Iend, Jstr, Jend,                    &
     &                     IstrC, IendC, JstrC, JendC,                  &
     &                     IstrP, IendP, JstrP, JendP,                  &
     &                     DPhi, phi_tab,                               &
     &                     fac_dx, fac_dy, fac_jac,                     &
     &                     pmask_act, rmask_ist,                        &
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
            IF (Master) WRITE(stdout,'(a,i3,a,1p,e12.4,a,e12.4)')      &
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
!  Velocity cap as a Picard-update damper (2026-10-24 redesign).
!
!  The former post-solve truncation (u -> u*cap/|u| whenever |u|
!  exceeded the cap) is incompatible with the convergence criteria
!  and the transport:  a binding cap pins the iterate, the two time
!  levels become identical, nl_resid evaluates to 0/0 = 0 and Picard
!  "converges" on the first iteration with the full cap velocity
!  released into the advection (CFL blow-up, the 2026-10-24 4 km
!  scan failure).  Instead the cap now bounds the UPDATE relative to
!  the previous iterate:  where the relaxed update would push the
!  node past the cap, the step is shortened so the node approaches
!  the cap magnitude along the update direction but never gets
!  clamped onto it.  Nodes whose solution is below the cap -- the
!  converged regime everywhere in normal physics -- are untouched,
!  so the fixed point is exactly the uncapped solution.
!
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            IF (pmask_act(i,j).eq.1.0_r8) THEN
!
!  Picard relaxation u = alpha*u_new + (1-alpha)*u_old.  Prevents
!  oscillation between high/low viscosity regimes.  Skipped on
!  the first iteration, where lsold holds the (physical) initial
!  guess rather than a Picard iterate.
!
              IF (nl_iter.gt.1) THEN
                Si(i,j,lsnew,isIuis)=                                   &
     &                    picard_alpha(ng) *Si(i,j,lsnew,isIuis)+       &
     &            (1.0_r8-picard_alpha(ng))*Si(i,j,lsold,isIuis)
                Si(i,j,lsnew,isIvis)=                                   &
     &                    picard_alpha(ng) *Si(i,j,lsnew,isIvis)+       &
     &            (1.0_r8-picard_alpha(ng))*Si(i,j,lsold,isIvis)
              END IF
!
!  Cap the relaxed UPDATE:  the allowed magnitude grows by at most
!  50 percent of the cap per Picard iteration beyond the cap --
!  overshooting nodes approach the cap region geometrically instead
!  of being clamped onto it, and a diverging iterate (u_prev far
!  above the cap) is pulled BACK toward the cap each iteration.
!  Nodes below the allowed magnitude are untouched, so the fixed
!  point is exactly the uncapped solution.
!
              u_slide=SQRT(Si(i,j,lsnew,isIuis)**2+                     &
     &                     Si(i,j,lsnew,isIvis)**2)
              u_allow=1.5_r8*ssa_max_speed(ng)
              IF (u_slide.gt.u_allow) THEN
                u_prev=SQRT(Si(i,j,lsold,isIuis)**2+                    &
     &                      Si(i,j,lsold,isIvis)**2)
                IF (u_prev.gt.u_allow) THEN
                  u_allow=MAX(u_allow, 0.9_r8*u_prev)
                END IF
                IF (u_slide.gt.u_allow) THEN
                  cff=u_allow/u_slide
                  Si(i,j,lsnew,isIuis)=cff*Si(i,j,lsnew,isIuis)
                  Si(i,j,lsnew,isIvis)=cff*Si(i,j,lsnew,isIvis)
                END IF
              END IF
            END IF
          END DO
        END DO
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
     &                     Istr, Iend, Jstr, Jend,                    &
     &                   IstrC, IendC, JstrC, JendC,                    &
     &                   IstrP, IendP, JstrP, JendP,                    &
     &                   DPhi, phi_tab,                                 &
     &                   fac_dx, fac_dy, fac_jac,                       &
     &                   pmask_act, rmask_ist,                          &
     &                   Fi(:,:,ifInuh),                                &
     &                   Fi(:,:,ifIbta),                                &
     &                   Si(:,:,lsnew,isIuis),                          &
     &                   Si(:,:,lsnew,isIvis),                          &
     &                   Au, Av)
!
!-----------------------------------------------------------------------
!  Picard convergence check:  dual criterion.        END IF
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
!  Dirichlet-constrained components are excluded from the residual
!  norms (review finding M-8):  at closed/clamped boundary nodes
!  ssa_lbc_zero keeps the CG residual at zero, so the local true
!  residual b - Au is the constraint reaction, O(driving stress),
!  and never decays -- including it gives criterion 2 an irreducible
!  floor and blinds the stagnation detector for clamped-boundary
!  configurations.  The node sets mirror ssa_lbc_zero.
!
        u_dW=(.not.EWperiodic(ng)).and.DOMAIN(ng)%Western_Edge(tile)   &
     &      .and.(LBC(iwest,ibISF(isIuis),ng)%closed .or.              &
     &            LBC(iwest,ibISF(isIuis),ng)%clamped)
        u_dE=(.not.EWperiodic(ng)).and.DOMAIN(ng)%Eastern_Edge(tile)   &
     &      .and.(LBC(ieast,ibISF(isIuis),ng)%closed .or.              &
     &            LBC(ieast,ibISF(isIuis),ng)%clamped)
        u_dS=(.not.NSperiodic(ng)).and.DOMAIN(ng)%Southern_Edge(tile)  &
     &      .and.(LBC(isouth,ibISF(isIuis),ng)%closed .or.             &
     &            LBC(isouth,ibISF(isIuis),ng)%clamped)
        u_dN=(.not.NSperiodic(ng)).and.DOMAIN(ng)%Northern_Edge(tile)  &
     &      .and.(LBC(inorth,ibISF(isIuis),ng)%closed .or.             &
     &            LBC(inorth,ibISF(isIuis),ng)%clamped)
        v_dW=(.not.EWperiodic(ng)).and.DOMAIN(ng)%Western_Edge(tile)   &
     &      .and.(LBC(iwest,ibISF(isIvis),ng)%closed .or.              &
     &            LBC(iwest,ibISF(isIvis),ng)%clamped)
        v_dE=(.not.EWperiodic(ng)).and.DOMAIN(ng)%Eastern_Edge(tile)   &
     &      .and.(LBC(ieast,ibISF(isIvis),ng)%closed .or.              &
     &            LBC(ieast,ibISF(isIvis),ng)%clamped)
        v_dS=(.not.NSperiodic(ng)).and.DOMAIN(ng)%Southern_Edge(tile)  &
     &      .and.(LBC(isouth,ibISF(isIvis),ng)%closed .or.             &
     &            LBC(isouth,ibISF(isIvis),ng)%clamped)
        v_dN=(.not.NSperiodic(ng)).and.DOMAIN(ng)%Northern_Edge(tile)  &
     &      .and.(LBC(inorth,ibISF(isIvis),ng)%closed .or.             &
     &            LBC(inorth,ibISF(isIvis),ng)%clamped)
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
              IF (.not.((u_dW.and.i.eq.IstrP).or.                       &
     &                  (u_dE.and.i.eq.IendP).or.                       &
     &                  (u_dS.and.j.eq.JstrP).or.                       &
     &                  (u_dN.and.j.eq.JendP))) THEN
                norms(3)=norms(3)+(Au(i,j)-Fi(i,j,ifItdx))**2
                norms(4)=norms(4)+Fi(i,j,ifItdx)**2
              END IF
              IF (.not.((v_dW.and.i.eq.IstrP).or.                       &
     &                  (v_dE.and.i.eq.IendP).or.                       &
     &                  (v_dS.and.j.eq.JstrP).or.                       &
     &                  (v_dN.and.j.eq.JendP))) THEN
                norms(3)=norms(3)+(Av(i,j)-Fi(i,j,ifItdy))**2
                norms(4)=norms(4)+Fi(i,j,ifItdy)**2
              END IF
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
            IF (rhs_norm.gt.eps) THEN
              WRITE(stdout,'(A,ES12.4)')                                &
     &          '  Picard true-resid ||Au-b||/||b|| =',                 &
     &          SQRT(resid_norm/rhs_norm)
            END IF
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
      END IF  ! n_ssa_fd FE/FD branch
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
     &                               pm, pn, h, zeta_avg,               &
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
      real(r8), intent(in)  :: h(LBi:,LBj:)
      real(r8), intent(in)  :: zeta_avg(LBi:,LBj:)
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
      real(r8), intent(in)  :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: zeta_avg(LBi:UBi,LBj:UBj)
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
      integer :: iw, ie, js, jn
      integer :: cw, ce, cs, cn, cc
      integer :: IstrC, IendC, JstrC, JendC
!
      real(r8) :: hi_psi
      real(r8) :: dSdx, dSdy, dx_rho, dy_rho, jac_avg
      real(r8) :: zs_sum, taper_w
      real(r8) :: cff, cff1
      real(r8) :: D_c, hx_c, hy_c
!
      real(r8), dimension(LBi:UBi,LBj:UBj) :: zs_smooth
      real(r8), dimension(LBi:UBi,LBj:UBj) :: zs_smooth_temp
      real(r8), dimension(LBi:UBi,LBj:UBj) :: taperfld
!
#include "set_bounds.h"
!
!  Extended cell-assembly bounds (the same convention as the other
!  element-assembly routines:  the western/southern halo ring is
!  integrated so seam nodes receive complete contributions).
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
!-----------------------------------------------------------------------
!  Branch-aware driving stress (n_dsbranch > 0, B5.5 ablation A'):
!  PISM SSAFD recipe -- the surface gradient at each CELL uses
!  one-sided differences within its own flotation branch, never
!  differencing across a grounded/floating interface (where the
!  surface elevation switches between the S = -h + H and
!  S = (1-rho_i/rho_w) H branches).  The audit (2026-10-27) showed
!  the centered cross-branch node gradient puts ~5.7x PISM's
!  driving stress on the first floating cell behind the grounding
!  line (~18 kPa excess), the surviving candidate for the ~85 km
!  westward equilibrium bias;  the ice front is unaffected (the
!  constant extension above already flattens the cross-front
!  gradient and the front-load line load is carried unchanged).
!  The cell stresses are scattered to the 4 corner PSI nodes with
!  the Q1 hat weight A_e/4 (the same lumping as jac_avg), which
!  reproduces the centered node gradient to O(dx^2) in the
!  interior, so only the GL band differs from the default path.
!
      IF (n_dsbranch(ng).gt.0) THEN
!
!  Nodal taper weights (same rule as the default path above).
!
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            IF (n_taper(ng).le.0) THEN
              taperfld(i,j)=1.0_r8
            ELSE
              dist_margin=n_taper(ng)
              DO dj=-1,1
                DO di=-1,1
                  IF (di.ne.0.or.dj.ne.0) THEN
                    IF (i+di.ge.LBi.and.i+di.le.UBi.and.             &
     &                  j+dj.ge.LBj.and.j+dj.le.UBj) THEN
                      IF (pmarg(i+di,j+dj).eq.1.0_r8) THEN
                        dist_margin=MIN(dist_margin,                  &
     &                                  MAX(ABS(di),ABS(dj)))
                      END IF
                    END IF
                  END IF
                END DO
              END DO
              taperfld(i,j)=MIN(REAL(dist_margin,r8)/                 &
     &                          REAL(n_taper(ng),r8),1.0_r8)
            END IF
          END DO
        END DO
!
!  Overwrite the nodal forces with the branch-aware assembly.
!
        DO j=JstrP,JendP
          DO i=IstrP,IendP
            tdx(i,j)=0.0_r8
            tdy(i,j)=0.0_r8
          END DO
        END DO
!
!  Cell loop over the assembly ring (the same reach as ssa_matvec)
!  so seam nodes receive complete contributions from both tiles.
!  Branch code:  0 = ice-free, 1 = grounded, 2 = floating, using
!  the window-mean criterion shared with isf_stress and
!  ssa_front_load (D_avg = h + zeta_avg - Dcrit).
!
        DO j=JstrC,JendC
          DO i=IstrC,IendC
            IF (hi(i,j).gt.0.0_r8) THEN
!
!  Branch codes of the cell and its 4 neighbors.
!
              D_c=h(i,j)+zeta_avg(i,j)-Dcrit(ng)
              IF (rho_i*hi(i,j).gt.rho0*MAX(D_c,0.0_r8)) THEN
                cc=1
              ELSE
                cc=2
              END IF
              IF (i-1.ge.LBi) THEN
                D_c=h(i-1,j)+zeta_avg(i-1,j)-Dcrit(ng)
                IF (hi(i-1,j).le.0.0_r8) THEN
                  cw=0
                ELSE IF (rho_i*hi(i-1,j).gt.rho0*MAX(D_c,0.0_r8)) THEN
                  cw=1
                ELSE
                  cw=2
                END IF
              ELSE
                cw=0
              END IF
              IF (i+1.le.UBi) THEN
                D_c=h(i+1,j)+zeta_avg(i+1,j)-Dcrit(ng)
                IF (hi(i+1,j).le.0.0_r8) THEN
                  ce=0
                ELSE IF (rho_i*hi(i+1,j).gt.rho0*MAX(D_c,0.0_r8)) THEN
                  ce=1
                ELSE
                  ce=2
                END IF
              ELSE
                ce=0
              END IF
              IF (j-1.ge.LBj) THEN
                D_c=h(i,j-1)+zeta_avg(i,j-1)-Dcrit(ng)
                IF (hi(i,j-1).le.0.0_r8) THEN
                  cs=0
                ELSE IF (rho_i*hi(i,j-1).gt.rho0*MAX(D_c,0.0_r8)) THEN
                  cs=1
                ELSE
                  cs=2
                END IF
              ELSE
                cs=0
              END IF
              IF (j+1.le.UBj) THEN
                D_c=h(i,j+1)+zeta_avg(i,j+1)-Dcrit(ng)
                IF (hi(i,j+1).le.0.0_r8) THEN
                  cn=0
                ELSE IF (rho_i*hi(i,j+1).gt.rho0*MAX(D_c,0.0_r8)) THEN
                  cn=1
                ELSE
                  cn=2
                END IF
              ELSE
                cn=0
              END IF
!
!  Branch-aware one-sided weights:  zero across a grounded/floating
!  interface only (ice-free neighbours keep the constant-extension
!  gradient;  the front load is the sole cliff force there).
!
              IF ((cc.eq.1.and.cw.eq.2).or.(cc.eq.2.and.cw.eq.1)) THEN
                iw=0
              ELSE
                iw=1
              END IF
              IF ((cc.eq.1.and.ce.eq.2).or.(cc.eq.2.and.ce.eq.1)) THEN
                ie=0
              ELSE
                ie=1
              END IF
              IF ((cc.eq.1.and.cs.eq.2).or.(cc.eq.2.and.cs.eq.1)) THEN
                js=0
              ELSE
                js=1
              END IF
              IF ((cc.eq.1.and.cn.eq.2).or.(cc.eq.2.and.cn.eq.1)) THEN
                jn=0
              ELSE
                jn=1
              END IF
!
!  Weighted one-sided gradient of the (extended) surface.
!
              IF (iw+ie.gt.0) THEN
                hx_c=(REAL(iw,r8)*(zs_smooth(i,j)-zs_smooth(i-1,j))+   &
     &                REAL(ie,r8)*(zs_smooth(i+1,j)-zs_smooth(i,j)))/  &
     &               (REAL(iw+ie,r8)/pm(i,j))
              ELSE
                hx_c=0.0_r8
              END IF
              IF (js+jn.gt.0) THEN
                hy_c=(REAL(js,r8)*(zs_smooth(i,j)-zs_smooth(i,j-1))+   &
     &                REAL(jn,r8)*(zs_smooth(i,j+1)-zs_smooth(i,j)))/  &
     &               (REAL(js+jn,r8)/pn(i,j))
              ELSE
                hy_c=0.0_r8
              END IF
!
!  Cell driving stress scattered to the 4 corner PSI nodes with
!  the Q1 hat weight A_e/4 (lumped-area consistent).
!
              cff1=-rho_i*g*hi(i,j)*0.25_r8/(pm(i,j)*pn(i,j))
!
              IF (pmask_act(i,j).eq.1.0_r8.and.                         &
     &            i.ge.IstrP.and.i.le.IendP.and.                        &
     &            j.ge.JstrP.and.j.le.JendP) THEN
                tdx(i,j)=tdx(i,j)+cff1*hx_c*taperfld(i,j)
                tdy(i,j)=tdy(i,j)+cff1*hy_c*taperfld(i,j)
              END IF
              IF (pmask_act(i+1,j).eq.1.0_r8.and.                       &
     &            i+1.ge.IstrP.and.i+1.le.IendP.and.                    &
     &            j.ge.JstrP.and.j.le.JendP) THEN
                tdx(i+1,j)=tdx(i+1,j)+cff1*hx_c*taperfld(i+1,j)
                tdy(i+1,j)=tdy(i+1,j)+cff1*hy_c*taperfld(i+1,j)
              END IF
              IF (pmask_act(i,j+1).eq.1.0_r8.and.                       &
     &            i.ge.IstrP.and.i.le.IendP.and.                        &
     &            j+1.ge.JstrP.and.j+1.le.JendP) THEN
                tdx(i,j+1)=tdx(i,j+1)+cff1*hx_c*taperfld(i,j+1)
                tdy(i,j+1)=tdy(i,j+1)+cff1*hy_c*taperfld(i,j+1)
              END IF
              IF (pmask_act(i+1,j+1).eq.1.0_r8.and.                     &
     &            i+1.ge.IstrP.and.i+1.le.IendP.and.                    &
     &            j+1.ge.JstrP.and.j+1.le.JendP) THEN
                tdx(i+1,j+1)=tdx(i+1,j+1)+cff1*hx_c*taperfld(i+1,j+1)
                tdy(i+1,j+1)=tdy(i+1,j+1)+cff1*hy_c*taperfld(i+1,j+1)
              END IF
            END IF
          END DO
        END DO
      END IF
!
      RETURN
      END SUBROUTINE ssa_driving_stress
!
!***********************************************************************
      SUBROUTINE ssa_ds_fem (ng, tile, model,                           &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Istr, Iend, Jstr, Jend,                   &
     &                        IstrP, IendP, JstrP, JendP,               &
     &                        DPhi, phi_tab,                            &
     &                        pm, pn, h, zeta_avg,                      &
     &                        rmask_ist, pmask_act,                     &
     &                        hi,                                       &
     &                        tdx, tdy)
!***********************************************************************
!
!  B5.8:  PISM-FEM style driving stress (SSAFEM::driving_stress):
!  element-wise quadrature-point evaluation from the SMOOTH bed and
!  thickness fields, with the flotation branch determined POINTWISE
!  at each quadrature point.  The surface gradient is never taken
!  across the grounded/floating branch boundary (the B5.5 audit
!  measured a 5.7x local RHS difference from the centered
!  differencing of the branch-mixed surface field, and the B5.7
!  three-operator comparison attributed the ROMS-FE 12-44% velocity
!  excess primarily to this driving-stress construction difference).
!
!  At each quadrature point q of element (i,j):
!    b_q, H_q, sl_q   = interpolated bed, thickness, sea level
!    b_x, b_y, H_x, H_y = shape-function gradients (smooth fields)
!    grounded_q        = rho_i*H_q > rho0*max(sl_q - b_q, 0)
!    h_x               = grounded_q ? b_x + H_x : alpha*H_x
!    tau_d,q           = -rho_i*g*H_q*(h_x, h_y)
!  scattered to the 4 corner PSI nodes with test-function weights,
!  same assembly convention as ssa_matvec.  ssa_front_load remains
!  the cliff/boundary term (B5.3 verified equivalent to PISM cfbc).
!
      implicit none
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Istr, Iend, Jstr, Jend
      integer, intent(in) :: IstrP, IendP, JstrP, JendP
!
      real(r8), intent(in)  :: DPhi(4,4,2)
      real(r8), intent(in)  :: phi_tab(4,4)
!
      real(r8), intent(in)  :: pm(LBi:,LBj:)
      real(r8), intent(in)  :: pn(LBi:,LBj:)
      real(r8), intent(in)  :: h(LBi:,LBj:)
      real(r8), intent(in)  :: zeta_avg(LBi:,LBj:)
      real(r8), intent(in)  :: rmask_ist(LBi:,LBj:)
      real(r8), intent(in)  :: pmask_act(LBi:,LBj:)
      real(r8), intent(in)  :: hi(LBi:,LBj:)
      real(r8), intent(out) :: tdx(LBi:,LBj:)
      real(r8), intent(out) :: tdy(LBi:,LBj:)
!
      integer :: i, j, nq, k, i1, j1
      integer :: IstrC, IendC, JstrC, JendC
!
      real(r8) :: alpha_c
      real(r8) :: bed_k(4), thk_k(4), zeta_k(4)
      real(r8) :: b_q, H_q, sl_q, dep_q
      real(r8) :: b_x, b_y, H_x, H_y
      real(r8) :: hx_q, hy_q, taud_x, taud_y
      real(r8) :: phv, dpx, dpy
      real(r8) :: f_dx, f_dy, f_jac
!
!
!  Assembly-ring bounds (the same convention as ssa_matvec).
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
      alpha_c=1.0_r8-rho_i/rho0
!
!  Clear output.
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          tdx(i,j)=0.0_r8
          tdy(i,j)=0.0_r8
        END DO
      END DO
!
      DO j=JstrC,JendC
        DO i=IstrC,IendC
          IF (hi(i,j).gt.0.0_r8) THEN
            f_dx=2.0_r8*pm(i,j)
            f_dy=2.0_r8*pn(i,j)
            f_jac=1.0_r8/(pm(i,j)*pn(i,j))
!
            bed_k(1)=-h(i,j)
            bed_k(2)=-h(i+1,j)
            bed_k(3)=-h(i,j+1)
            bed_k(4)=-h(i+1,j+1)
            thk_k(1)=hi(i,j)
            thk_k(2)=hi(i+1,j)
            thk_k(3)=hi(i,j+1)
            thk_k(4)=hi(i+1,j+1)
            zeta_k(1)=zeta_avg(i,j)
            zeta_k(2)=zeta_avg(i+1,j)
            zeta_k(3)=zeta_avg(i,j+1)
            zeta_k(4)=zeta_avg(i+1,j+1)
!
            DO nq=1,4
              b_q=0.0_r8
              H_q=0.0_r8
              sl_q=0.0_r8
              b_x=0.0_r8
              b_y=0.0_r8
              H_x=0.0_r8
              H_y=0.0_r8
              DO k=1,4
                phv=phi_tab(k,nq)
                dpx=DPhi(k,nq,1)*f_dx
                dpy=DPhi(k,nq,2)*f_dy
                b_q=b_q+phv*bed_k(k)
                H_q=H_q+phv*thk_k(k)
                sl_q=sl_q+phv*zeta_k(k)
                b_x=b_x+dpx*bed_k(k)
                b_y=b_y+dpy*bed_k(k)
                H_x=H_x+dpx*thk_k(k)
                H_y=H_y+dpy*thk_k(k)
              END DO
!
!  Pointwise flotation branch at the quadrature point.
!
              dep_q=MAX(sl_q-b_q,0.0_r8)
              IF (rho_i*H_q.gt.rho0*dep_q) THEN
                hx_q=b_x+H_x
                hy_q=b_y+H_y
              ELSE
                hx_q=alpha_c*H_x
                hy_q=alpha_c*H_y
              END IF
!
              taud_x=-rho_i*g*H_q*hx_q
              taud_y=-rho_i*g*H_q*hy_q
!
!  Scatter to the 4 corner PSI nodes (test-function weights,
!  ownership filter;  same convention as ssa_matvec).
!
              DO k=1,4
                IF (k.eq.1) THEN
                  i1=i
                  j1=j
                ELSE IF (k.eq.2) THEN
                  i1=i+1
                  j1=j
                ELSE IF (k.eq.3) THEN
                  i1=i
                  j1=j+1
                ELSE
                  i1=i+1
                  j1=j+1
                END IF
                IF (pmask_act(i1,j1).eq.1.0_r8.and.                    &
     &              i1.ge.IstrP.and.i1.le.IendP.and.                    &
     &              j1.ge.JstrP.and.j1.le.JendP) THEN
                  tdx(i1,j1)=tdx(i1,j1)+                               &
     &                0.25_r8*f_jac*phi_tab(k,nq)*taud_x
                  tdy(i1,j1)=tdy(i1,j1)+                               &
     &                0.25_r8*f_jac*phi_tab(k,nq)*taud_y
                END IF
              END DO
            END DO
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE ssa_ds_fem
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
!  closed or clamped edge is a wall (never a front).  (2026-11-08
!  correction:  an experiment replacing the gradient override with
!  the ghost-ring mask value -- intended to let ice-continuation
!  edges stay load-free -- killed the MISMIP+ east-edge front load,
!  because the restart's ghost ring holds generator-filled ice at
!  the calving front.  The LBC rule stands.  Ice-continuation edges
!  would need a distinct LBC type plus the membrane-traction
!  operator term;  see work_plan 2026-11-08.)
!
!  edge is an open front (the LBC-alone rule;  MISMIP+ east edge).
!  0 -- the ghost-ring mask decides (ghost ice = continuation, no
!  front;  ghost ocean = open front), and ssa_matvec supplies the
!  membrane traction of the continuing ice (the SSA_TEST front
!  channel symmetry edges).
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
!  Floor at 100 K:  a garbage ice temperature at or below absolute
!  zero would overflow EXP(-Q/RT) or divide by zero below.
              T_star=MAX(ti(i,j)+273.15_r8+8.66E-4_r8*0.5_r8*hi(i,j), &
     &                   100.0_r8)
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
     &                            fac_dx, fac_dy, fac_jac,              &
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
      real(r8), intent(in)  :: fac_dx(LBi:,LBj:)
      real(r8), intent(in)  :: fac_dy(LBi:,LBj:)
      real(r8), intent(in)  :: fac_jac(LBi:,LBj:)
      real(r8), intent(in)  :: pmask_act(LBi:,LBj:)
      real(r8), intent(in)  :: nuH(LBi:,LBj:)
      real(r8), intent(in)  :: beta(LBi:,LBj:)
      real(r8), intent(out) :: diag_u(LBi:,LBj:)
      real(r8), intent(out) :: diag_v(LBi:,LBj:)
# else
      real(r8), intent(in)  :: fac_dx(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: fac_dy(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: fac_jac(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pmask_act(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: rmask_ist(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: nuH(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: beta(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: diag_u(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: diag_v(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, i1, j1, i2, j2, kk, nq, k_node
!
      real(r8) :: dphi_x, dphi_y, phv, nuhq
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
                  dphi_x=DPhi(k_node,nq,1)*fac_dx(i,j)
                  dphi_y=DPhi(k_node,nq,2)*fac_dy(i,j)
                  phv=phi_tab(k_node,nq)
!
!  B5.6 variant:  nuH at the quadrature point by bilinear
!  interpolation (matching ssa_matvec;  the diagonal must stay the
!  operator diagonal of the same discretization).
!
                  IF (n_nuh_interp(ng).gt.0) THEN
                    nuhq=0.0_r8
                    DO kk=1,4
                      IF (kk.eq.1) THEN
                        i2=i
                        j2=j
                      ELSE IF (kk.eq.2) THEN
                        i2=i+1
                        j2=j
                      ELSE IF (kk.eq.3) THEN
                        i2=i
                        j2=j+1
                      ELSE
                        i2=i+1
                        j2=j+1
                      END IF
                      nuhq=nuhq+nuH(i2,j2)*phi_tab(kk,nq)
                    END DO
                  ELSE
                    nuhq=nuH(i,j)
                  END IF
!
                  diag_u(i1,j1)=diag_u(i1,j1)+                          &
     &                          0.25_r8*fac_jac(i,j)*nuhq*              &
     &                          (4.0_r8*dphi_x*dphi_x+dphi_y*dphi_y)+   &
     &                          0.25_r8*fac_jac(i,j)*beta(i,j)*phv**2
                  diag_v(i1,j1)=diag_v(i1,j1)+                          &
     &                          0.25_r8*fac_jac(i,j)*nuhq*              &
     &                          (4.0_r8*dphi_y*dphi_y+dphi_x*dphi_x)+   &
     &                          0.25_r8*fac_jac(i,j)*beta(i,j)*phv**2
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
     &                       Istr, Iend, Jstr, Jend,                    &
     &                       IstrC, IendC, JstrC, JendC,                &
     &                       IstrP, IendP, JstrP, JendP,                &
     &                       DPhi, phi_tab,                             &
     &                       fac_dx, fac_dy, fac_jac,                   &
     &                       pmask_act, rmask_ist,                      &
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
      integer, intent(in) :: Istr, Iend, Jstr, Jend
      integer, intent(in) :: IstrC, IendC, JstrC, JendC
      integer, intent(in) :: IstrP, IendP, JstrP, JendP
!
      real(r8), intent(in)    :: DPhi(4,4,2)
      real(r8), intent(in)    :: phi_tab(4,4)
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: fac_dx(LBi:,LBj:)
      real(r8), intent(in)  :: fac_dy(LBi:,LBj:)
      real(r8), intent(in)  :: fac_jac(LBi:,LBj:)
      real(r8), intent(in)  :: pmask_act(LBi:,LBj:)
      real(r8), intent(in)  :: rmask_ist(LBi:,LBj:)
      real(r8), intent(in)  :: nuH(LBi:,LBj:)
      real(r8), intent(in)  :: beta(LBi:,LBj:)
      real(r8), intent(inout) :: xu(LBi:,LBj:)
      real(r8), intent(inout) :: xv(LBi:,LBj:)
      real(r8), intent(out) :: Au(LBi:,LBj:)
      real(r8), intent(out) :: Av(LBi:,LBj:)
# else
      real(r8), intent(in)  :: fac_dx(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: fac_dy(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: fac_jac(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pmask_act(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: rmask_ist(LBi:UBi,LBj:UBj)
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
      real(r8), dimension(4) :: dpx, dpy, phv4
      real(r8) :: uq, vq, exx, eyy, exy, nuhq
      real(r8) :: exx_e, eyy_e, Nxx_e, Nyy_e
      logical :: cont_s, cont_n, cont_w, cont_e
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
!  Clear the output pair over exactly the reach of the scatters below:
!  cell corners of the assembly ring, [IstrC..IendC+1]x[JstrC..JendC+1].
!  This always lies inside the array bounds (IendC=Iend, so the upper
!  edge is Iend+1 <= UBi with equality on eastern/northern domain-edge
!  tiles), whereas the original PSI-ring loop DO j=JstrP-1,JendP+1 ran
!  to IendP+1/JendP+1, which exceed UBi/UBj on domain-edge tiles
!  (IendP=IendR=Iend+1 there) and wrote one full row/column past the
!  Au/Av allocations on EVERY call -- the heap corruption behind the
!  L11 crashes (smashed malloc metadata surfacing as the set_zeta
!  null-descriptor SIGSEGV, netCDF/HDF5 malloc aborts, and freelist
!  surgery inside live allocations).  The range below also covers all
!  owned PSI nodes (IstrC <= IstrP, IendC+1 >= IendP), so every node
!  the consumers read is either scattered or cleared.
!
      DO j=JstrC,JendC+1
        DO i=IstrC,IendC+1
          Au(i,j)=0.0_r8
          Av(i,j)=0.0_r8
        END DO
      END DO
!
      DO j=JstrC,JendC
        DO i=IstrC,IendC
          IF (nuH(i,j).gt.0.0_r8) THEN
!
!  Continuation-edge membrane traction (LBC type Cnt):  at a Cnt
!  domain edge the sheet continues beyond the domain (declared by
!  the user, no ghost-state test), and the row convention
!  +int N:grad(phi) = b carries the boundary integral +oint (N.n) phi
!  in the OPERATOR.  Inert unless the edge LBC is Cnt (a gradient
!  edge is an open front and carries the ssa_front_load pressure
!  instead;  a closed edge is a wall).
!
            cont_s=.FALSE.; cont_n=.FALSE.
            cont_w=.FALSE.; cont_e=.FALSE.
            IF ((.not.EWperiodic(ng)).and.(.not.NSperiodic(ng))) THEN
              IF ((j.eq.Jstr).and.DOMAIN(ng)%Southern_Edge(tile)) THEN
                cont_s=LBC(isouth,ibISF(isIvis),ng)%continuation
              END IF
              IF ((j.eq.Jend).and.DOMAIN(ng)%Northern_Edge(tile)) THEN
                cont_n=LBC(inorth,ibISF(isIvis),ng)%continuation
              END IF
              IF ((i.eq.Istr).and.DOMAIN(ng)%Western_Edge(tile)) THEN
                cont_w=LBC(iwest,ibISF(isIuis),ng)%continuation
              END IF
              IF ((i.eq.Iend).and.DOMAIN(ng)%Eastern_Edge(tile)) THEN
                cont_e=LBC(ieast,ibISF(isIuis),ng)%continuation
              END IF
            END IF
            IF (cont_s.or.cont_n.or.cont_w.or.cont_e) THEN
              exx_e=0.5_r8*((xu(i+1,j  )-xu(i,  j  ))+                  &
     &                      (xu(i+1,j+1)-xu(i,  j+1)))*fac_dx(i,j)*0.5_r8
              eyy_e=0.5_r8*((xv(i,  j+1)-xv(i,  j  ))+                  &
     &                      (xv(i+1,j+1)-xv(i+1,j  )))*fac_dy(i,j)*0.5_r8
              Nxx_e=2.0_r8*nuH(i,j)*(2.0_r8*exx_e+eyy_e)
              Nyy_e=2.0_r8*nuH(i,j)*(2.0_r8*eyy_e+exx_e)
!
!  South:  (N.n)_y = -N_yy -> the operator gains +N_yy*dx/2 per node.
!  North:  -N_yy*dx/2.  West:  +N_xx*dy/2.  East:  -N_xx*dy/2.
!
              IF (cont_s) THEN
                IF (pmask_act(i,j).eq.1.0_r8)                           &
     &            Av(i,j)=Av(i,j)+Nyy_e/fac_dy(i,j)
                IF (pmask_act(i+1,j).eq.1.0_r8)                         &
     &            Av(i+1,j)=Av(i+1,j)+Nyy_e/fac_dy(i,j)
              END IF
              IF (cont_n) THEN
                IF (pmask_act(i,j+1).eq.1.0_r8)                         &
     &            Av(i,j+1)=Av(i,j+1)-Nyy_e/fac_dy(i,j)
                IF (pmask_act(i+1,j+1).eq.1.0_r8)                       &
     &            Av(i+1,j+1)=Av(i+1,j+1)-Nyy_e/fac_dy(i,j)
              END IF
              IF (cont_w) THEN
                IF (pmask_act(i,j).eq.1.0_r8)                           &
     &            Au(i,j)=Au(i,j)+Nxx_e/fac_dx(i,j)
                IF (pmask_act(i,j+1).eq.1.0_r8)                         &
     &            Au(i,j+1)=Au(i,j+1)+Nxx_e/fac_dx(i,j)
              END IF
              IF (cont_e) THEN
                IF (pmask_act(i+1,j).eq.1.0_r8)                         &
     &            Au(i+1,j)=Au(i+1,j)-Nxx_e/fac_dx(i,j)
                IF (pmask_act(i+1,j+1).eq.1.0_r8)                       &
     &            Au(i+1,j+1)=Au(i+1,j+1)-Nxx_e/fac_dx(i,j)
              END IF
            END IF
!
            DO nq=1,4
!
!  Interpolate the velocity and form the strain rates at the
!  quadrature point.  The shape-function values and derivatives at
!  this quadrature point (one per corner node) are formed once and
!  reused by the scatter below.
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
                phv4(k_node)=phi_tab(k_node,nq)
                dpx(k_node)=DPhi(k_node,nq,1)*fac_dx(i,j)
                dpy(k_node)=DPhi(k_node,nq,2)*fac_dy(i,j)
                uq=uq+xu(i1,j1)*phv4(k_node)
                vq=vq+xv(i1,j1)*phv4(k_node)
                exx=exx+xu(i1,j1)*dpx(k_node)
                eyy=eyy+xv(i1,j1)*dpy(k_node)
                exy=exy+0.5_r8*(xu(i1,j1)*dpy(k_node)+                  &
     &                           xv(i1,j1)*dpx(k_node))
              END DO
!
!  B5.6 membrane variant (n_nuh_interp > 0):  the effective
!  viscosity is bilinearly interpolated to the quadrature point from
!  the cell-centred nuH field (a continuous Q1 nuH over the mesh,
!  PISM SSAFEM-style) instead of the element-constant value --
!  capturing the sub-element nuH variation that the cell-centre
!  constant drops and that the B5.6.1 operator audit measured as a
!  60-140% membrane-divergence bias.
!
              IF (n_nuh_interp(ng).gt.0) THEN
                nuhq=0.0_r8
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
                  nuhq=nuhq+nuH(i1,j1)*phv4(k_node)
                END DO
              ELSE
                nuhq=nuH(i,j)
              END IF
!
!  Scatter to the 4 corner nodes (ownership filter).
!
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
                IF (pmask_act(i1,j1).eq.1.0_r8.and.                     &
     &              i1.ge.IstrP.and.i1.le.IendP.and.                    &
     &              j1.ge.JstrP.and.j1.le.JendP) THEN
!
                  Au(i1,j1)=Au(i1,j1)+                                  &
     &              0.25_r8*fac_jac(i,j)*nuhq*                           &
     &              (dpx(k_node)*(4.0_r8*exx+2.0_r8*eyy)+               &
     &               dpy(k_node)*(2.0_r8*exy))+                         &
     &              0.25_r8*fac_jac(i,j)*beta(i,j)*phv4(k_node)*uq
!
                  Av(i1,j1)=Av(i1,j1)+                                  &
     &              0.25_r8*fac_jac(i,j)*nuhq*                           &
     &              (dpy(k_node)*(4.0_r8*eyy+2.0_r8*exx)+               &
     &               dpx(k_node)*(2.0_r8*exy))+                         &
     &              0.25_r8*fac_jac(i,j)*beta(i,j)*phv4(k_node)*vq
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
!
      RETURN
      END SUBROUTINE ssa_matvec
!
!***********************************************************************
!  FD solver (n_ssa_fd > 0):  finite-difference discretization of
!  the SSA following the PISM SSAFD recipe.  Every piece below is
!  the one-for-one FD counterpart of an FE piece above:
!
!    ssa_fd_viscosity  <-> ssa_viscosity
!        staggered face viscosities nuhx (x-face between RHO(i,j)
!        and RHO(i+1,j)) and nuhy (y-face between RHO(i,j) and
!        RHO(i,j+1));  strain rates from face differences of the
!        RHO-point velocity with the 4-point cross term;  same Glen
!        invariant, hardness field, thin-ice strength extension and
!        (unused here) additive regularization as the FE path.
!
!    ssa_fd_driving    <-> ssa_driving_stress (+ ssa_front_load)
!        cell-centered driving stress with the branch-aware weighted
!        one-sided surface gradient (PISM compute_driving_stress)
!        and the marginal-cell front pressure difference
!        delta_p = H*(P_ice - P_water) folded into the cell RHS
!        (PISM cfbc), replacing the FE pair of the nodal gradient
!        plus the Q1 line load.
!
!    ssa_fd_matvec     <-> ssa_matvec
!        the 13-point FD stencil action (PISM fd_operator Maxima
!        coefficients) on the RHO-point velocity, with the marginal
!        switch weights (centered -> one-sided differences at
!        ice-adjacent-to-ice-free interfaces) and the cell-centered
!        beta on the diagonal.
!
!    ssa_fd_diag       <-> ssa_jacobi_diag
!        stencil diagonal + beta.
!
!  The solver state (search directions, residuals, preconditioned
!  vectors) lives on RHO points;  velocity enters the solve as the
!  4-corner PSI average and leaves as the 4-cell RHO average
!  (the two maps are transposes of each other, so interior velocity
!  magnitude is preserved to O(dx^2);  the PSI boundary values are
!  handled by isf_pibc as before).  LBC semantics:  the FE path
!  zeroes residual/search rows at closed/clamped PSI edge nodes
!  (ssa_lbc_zero);  the FD path carries the same intent by fixing
!  the RHO cells adjacent to closed/clamped edges (identity rows),
!  the closest cell-centered analogue.
!***********************************************************************
!
      SUBROUTINE ssa_fd_viscosity (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             Istr, Iend, Jstr, Jend,              &
     &                             pm, pn,                              &
     &                             hi, B_hard,                          &
     &                             urho, vrho,                          &
     &                             nuhx, nuhy)
!***********************************************************************
!
      implicit none
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Istr, Iend, Jstr, Jend
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: pm(LBi:,LBj:)
      real(r8), intent(in)  :: pn(LBi:,LBj:)
      real(r8), intent(in)  :: hi(LBi:,LBj:)
      real(r8), intent(in)  :: B_hard(LBi:,LBj:)
      real(r8), intent(in)  :: urho(LBi:,LBj:)
      real(r8), intent(in)  :: vrho(LBi:,LBj:)
      real(r8), intent(out) :: nuhx(LBi:,LBj:)
      real(r8), intent(out) :: nuhy(LBi:,LBj:)
# else
      real(r8), intent(in)  :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: hi(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: B_hard(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: urho(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: vrho(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: nuhx(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: nuhy(LBi:UBi,LBj:UBj)
# endif
!
      integer :: i, j
!
      real(r8) :: Hf, gam, ux, uy, vx, vy, wz
!
!
!  x-face viscosities (PISM o=0):  u_x from the 2-point face
!  difference, u_y from the 4-point average of the two straddled
!  faces;  face thickness is the mean of the two cells.
!
      DO j=MAX(LBj,Jstr-1),MIN(UBj,Jend)
        DO i=MAX(LBi,Istr-1),MIN(UBi-1,Iend)
          Hf=0.5_r8*(hi(i,j)+hi(i+1,j))
          IF (Hf.gt.0.0_r8) THEN
            ux=(urho(i+1,j)-urho(i,j))*pm(i,j)
            uy=0.25_r8*((urho(i,j+1)+urho(i+1,j+1)-                     &
     &                   urho(i,j-1)-urho(i+1,j-1))*pn(i,j))
            vx=(vrho(i+1,j)-vrho(i,j))*pm(i,j)
            vy=0.25_r8*((vrho(i,j+1)+vrho(i+1,j+1)-                     &
     &                   vrho(i,j-1)-vrho(i+1,j-1))*pn(i,j))
            wz=-(ux+vy)
            gam=0.5_r8*(ux*ux+vy*vy+wz*wz+0.5_r8*(uy+vx)*(uy+vx))+      &
     &          eps_glen_min(ng)**2
            IF (i+1.le.UBi) THEN
              nuhx(i,j)=0.25_r8*(B_hard(i,j)+B_hard(i+1,j))*            &
     &                  gam**((1.0_r8-n_glen(ng))/(2.0_r8*n_glen(ng)))* &
     &                  MAX(Hf,min_thickness(ng))
            ELSE
              nuhx(i,j)=0.5_r8*B_hard(i,j)*                             &
     &                  gam**((1.0_r8-n_glen(ng))/(2.0_r8*n_glen(ng)))* &
     &                  MAX(Hf,min_thickness(ng))
            END IF
          ELSE
            nuhx(i,j)=0.0_r8
          END IF
        END DO
        nuhx(MIN(UBi-1,Iend),j)=nuhx(MIN(UBi-2,Iend-1),j)
      END DO
!
!  y-face viscosities (PISM o=1).
!
      DO j=MAX(LBj,Jstr-1),MIN(UBj-1,Jend)
        DO i=MAX(LBi,Istr-1),MIN(UBi,Iend)
          Hf=0.5_r8*(hi(i,j)+hi(i,j+1))
          IF (Hf.gt.0.0_r8) THEN
            ux=0.25_r8*((urho(i+1,j)+urho(i+1,j+1)-                     &
     &                   urho(i-1,j)-urho(i-1,j+1))*pm(i,j))
            uy=(vrho(i,j+1)-vrho(i,j))*pn(i,j)
            vx=0.25_r8*((vrho(i+1,j)+vrho(i+1,j+1)-                     &
     &                   vrho(i-1,j)-vrho(i-1,j+1))*pm(i,j))
            vy=(vrho(i,j+1)-vrho(i,j))*pn(i,j)
            wz=-(ux+vy)
            gam=0.5_r8*(ux*ux+vy*vy+wz*wz+0.5_r8*(uy+vx)*(uy+vx))+      &
     &          eps_glen_min(ng)**2
            IF (j+1.le.UBj) THEN
              nuhy(i,j)=0.25_r8*(B_hard(i,j)+B_hard(i,j+1))*            &
     &                  gam**((1.0_r8-n_glen(ng))/(2.0_r8*n_glen(ng)))* &
     &                  MAX(Hf,min_thickness(ng))
            ELSE
              nuhy(i,j)=0.5_r8*B_hard(i,j)*                             &
     &                  gam**((1.0_r8-n_glen(ng))/(2.0_r8*n_glen(ng)))* &
     &                  MAX(Hf,min_thickness(ng))
            END IF
          ELSE
            nuhy(i,j)=0.0_r8
          END IF
        END DO
      END DO
      DO i=MAX(LBi,Istr-1),MIN(UBi,Iend)
        nuhy(i,MIN(UBj-1,Jend))=nuhy(i,MIN(UBj-2,Jend-1))
      END DO
!
      RETURN
      END SUBROUTINE ssa_fd_viscosity
!
!***********************************************************************
      SUBROUTINE ssa_fd_driving (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           Istr, Iend, Jstr, Jend,                &
     &                           pm, pn, h, zeta_avg,                   &
     &                           hi,                                    &
     &                           fdx, fdy)
!***********************************************************************
!
!  Cell-centered driving stress (PISM compute_driving_stress) plus
!  the marginal front-pressure difference (PISM assemble_rhs cfbc
!  block), folded into one cell RHS:  fdx = taud_x + (E-W)*delta_p/dx.
!  Surface gradients are branch-aware one-sided differences (the
!  audit-established PISM semantics);  ice-free cells carry a zero
!  RHS (their velocity rows are identity in ssa_fd_matvec).
!
      implicit none
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Istr, Iend, Jstr, Jend
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: pm(LBi:,LBj:)
      real(r8), intent(in)  :: pn(LBi:,LBj:)
      real(r8), intent(in)  :: h(LBi:,LBj:)
      real(r8), intent(in)  :: zeta_avg(LBi:,LBj:)
      real(r8), intent(in)  :: hi(LBi:,LBj:)
      real(r8), intent(out) :: fdx(LBi:,LBj:)
      real(r8), intent(out) :: fdy(LBi:,LBj:)
# else
      real(r8), intent(in)  :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: zeta_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: hi(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: fdx(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: fdy(LBi:UBi,LBj:UBj)
# endif
!
      integer :: i, j, ii, jj
      integer :: cw, ce, cs, cn, cc
      real(r8) :: cff
!
      real(r8) :: D_c, D_n
      real(r8) :: wsumx, wsumy
      real(r8) :: hx, hy, Pice, Pwat, dp, wc
      real(r8) :: wsum
      real(r8), dimension(LBi:UBi,LBj:UBj) :: sfld
!
!
      DO j=MAX(LBj,Jstr-1),MIN(UBj,Jend)
        DO i=MAX(LBi,Istr-1),MIN(UBi,Iend)
          fdx(i,j)=0.0_r8
          fdy(i,j)=0.0_r8
!
          IF (hi(i,j).gt.0.0_r8) THEN
!
!  Branch codes:  0 = ice-free, 1 = grounded, 2 = floating (window
!  criterion shared with isf_stress / ssa_front_load).
!
            D_c=h(i,j)+zeta_avg(i,j)-Dcrit(ng)
            IF (rho_i*hi(i,j).gt.rho0*MAX(D_c,0.0_r8)) THEN
              cc=1
            ELSE
              cc=2
            END IF
            IF (i.gt.LBi) THEN
              D_n=h(i-1,j)+zeta_avg(i-1,j)-Dcrit(ng)
              IF (hi(i-1,j).le.0.0_r8) THEN
                cw=0
              ELSE IF (rho_i*hi(i-1,j).gt.rho0*MAX(D_n,0.0_r8)) THEN
                cw=1
              ELSE
                cw=2
              END IF
            ELSE
              cw=0
            END IF
            IF (i.lt.UBi) THEN
              D_n=h(i+1,j)+zeta_avg(i+1,j)-Dcrit(ng)
              IF (hi(i+1,j).le.0.0_r8) THEN
                ce=0
              ELSE IF (rho_i*hi(i+1,j).gt.rho0*MAX(D_n,0.0_r8)) THEN
                ce=1
              ELSE
                ce=2
              END IF
            ELSE
              ce=0
            END IF
            IF (j.gt.LBj) THEN
              D_n=h(i,j-1)+zeta_avg(i,j-1)-Dcrit(ng)
              IF (hi(i,j-1).le.0.0_r8) THEN
                cs=0
              ELSE IF (rho_i*hi(i,j-1).gt.rho0*MAX(D_n,0.0_r8)) THEN
                cs=1
              ELSE
                cs=2
              END IF
            ELSE
              cs=0
            END IF
            IF (j.lt.UBj) THEN
              D_n=h(i,j+1)+zeta_avg(i,j+1)-Dcrit(ng)
              IF (hi(i,j+1).le.0.0_r8) THEN
                cn=0
              ELSE IF (rho_i*hi(i,j+1).gt.rho0*MAX(D_n,0.0_r8)) THEN
                cn=1
              ELSE
                cn=2
              END IF
            ELSE
              cn=0
            END IF
!
!  Surface elevation from the branch codes (s = H - h grounded,
!  s = (1-rho_i/rho_0) H floating):  self-contained so the RHS does
!  not depend on the isf_masks refresh order (Fi(ifIzsf) is not part
!  of the restart state).
!
          DO jj=LBj,UBj
            DO ii=LBi,UBi
              IF (hi(ii,jj).gt.0.0_r8) THEN
                D_n=h(ii,jj)+zeta_avg(ii,jj)-Dcrit(ng)
                IF (rho_i*hi(ii,jj).gt.rho0*MAX(D_n,0.0_r8)) THEN
                  sfld(ii,jj)=hi(ii,jj)-h(ii,jj)
                ELSE
                  sfld(ii,jj)=(1.0_r8-rho_i/rho0)*hi(ii,jj)
                END IF
              ELSE
                sfld(ii,jj)=0.0_r8
              END IF
            END DO
          END DO
!
!  Branch-aware weighted one-sided surface gradient:  grounded --
!  centered over the two grounded neighbors;  floating -- centered
!  over the two floating neighbors;  GL/marginal transitions use
!  the single in-branch one-sided difference (PISM weight()).
!  Ice-free neighbors are excluded (zs=0 extension would fake a
!  cliff;  the cliff force is the delta_p term below).
!
            wsumx=0.0_r8
            wsumy=0.0_r8
            hx=0.0_r8
            hy=0.0_r8
!
!  PISM SSAFD driving-stress gradient:  the two one-sided differences
!  enter with the SAME sign, (h.c-h.w)/dx + (h.e-h.c)/dx over
!  (west+east), summing to the centered difference.  The negative sign
!  on the west term below (found 2026-11-05 by the Python per-cell
!  recompute) turned the sum into a second difference divided by a
!  first-difference denominator:  fdx vanished on uniform-slope
!  interiors (the 2 km equilibrium FD collapse) while staying correct
!  on zero-slope geometry (the front test).
!
!
!  The x and y gradient normalizations are SEPARATE:  PISM divides
!  h_x by (west+east) and h_y by (south+north).  Sharing one wsum
!  (the 2026-11-05 verification find) divided each interior gradient
!  by 4 instead of 2, halving both driving-stress components.
!
            IF (cw.eq.cc.and.cc.ne.0) THEN
              hx=hx+(sfld(i,j)-sfld(i-1,j))*pm(i,j)
              wsumx=wsumx+1.0_r8
            END IF
            IF (ce.eq.cc.and.cc.ne.0) THEN
              hx=hx+(sfld(i+1,j)-sfld(i,j))*pm(i,j)
              wsumx=wsumx+1.0_r8
            END IF
            IF (cs.eq.cc.and.cc.ne.0) THEN
              hy=hy+(sfld(i,j)-sfld(i,j-1))*pn(i,j)
              wsumy=wsumy+1.0_r8
            END IF
            IF (cn.eq.cc.and.cc.ne.0) THEN
              hy=hy+(sfld(i,j+1)-sfld(i,j))*pn(i,j)
              wsumy=wsumy+1.0_r8
            END IF
            IF (wsumx.gt.0.0_r8) THEN
              hx=hx/wsumx
            END IF
            IF (wsumy.gt.0.0_r8) THEN
              hy=hy/wsumy
            END IF
!
!  Floating cells at the ice front:  halve the gradient (PISM's
!  constant-extrapolation half-weight at ocean fronts).
!
            IF (cc.eq.2.and.(ce.eq.0.or.cw.eq.0)) hx=0.5_r8*hx
            IF (cc.eq.2.and.(cn.eq.0.or.cs.eq.0)) hy=0.5_r8*hy
!
            fdx(i,j)=-rho_i*g*hi(i,j)*hx
            fdy(i,j)=-rho_i*g*hi(i,j)*hy
!
!  Marginal-cell front pressure difference (cfbc):  any ice-free
!  ocean neighbor in the 8-neighborhood (PISM is_marginal) adds
!  delta_p = H*(P_ice - P_water) per open face, along the outward
!  normal.  Fjord-wall override:  ice-free neighbor with bed above
!  the ice surface gets no pressure difference.
!  PISM average_water_column_pressure:  the water column at the
!  interface runs from the ice bottom (max(bed, flotation depth))
!  up to sea level -- NOT from the sea bed (a floating front would
!  otherwise carry the deep-basin pressure below its keel).
!
            IF (cw.eq.0.and.h(i-1,j).gt.0.0_r8) THEN
              IF (-h(i-1,j).lt.sfld(i,j)) THEN
                Pice=0.5_r8*rho_i*g*hi(i,j)
                wc=MAX(zeta_avg(i-1,j)                                   &
     &                -MAX(-h(i-1,j),                                   &
     &                     -rho_i/rho0*hi(i,j)),0.0_r8)
                Pwat=0.5_r8*rho0*g*wc*wc/MAX(hi(i,j),1.0_r8)
                dp=hi(i,j)*(Pice-Pwat)
                fdx(i,j)=fdx(i,j)-dp*pm(i,j)
              END IF
            END IF
            IF (ce.eq.0.and.h(i+1,j).gt.0.0_r8) THEN
              IF (-h(i+1,j).lt.sfld(i,j)) THEN
                Pice=0.5_r8*rho_i*g*hi(i,j)
                wc=MAX(zeta_avg(i+1,j)                                   &
     &                -MAX(-h(i+1,j),                                   &
     &                     -rho_i/rho0*hi(i,j)),0.0_r8)
                Pwat=0.5_r8*rho0*g*wc*wc/MAX(hi(i,j),1.0_r8)
                dp=hi(i,j)*(Pice-Pwat)
                fdx(i,j)=fdx(i,j)+dp*pm(i,j)
              END IF
            END IF
            IF (cs.eq.0.and.h(i,j-1).gt.0.0_r8) THEN
              IF (-h(i,j-1).lt.sfld(i,j)) THEN
                Pice=0.5_r8*rho_i*g*hi(i,j)
                wc=MAX(zeta_avg(i,j-1)                                   &
     &                -MAX(-h(i,j-1),                                   &
     &                     -rho_i/rho0*hi(i,j)),0.0_r8)
                Pwat=0.5_r8*rho0*g*wc*wc/MAX(hi(i,j),1.0_r8)
                dp=hi(i,j)*(Pice-Pwat)
                fdy(i,j)=fdy(i,j)-dp*pn(i,j)
              END IF
            END IF
            IF (cn.eq.0.and.h(i,j+1).gt.0.0_r8) THEN
              IF (-h(i,j+1).lt.sfld(i,j)) THEN
                Pice=0.5_r8*rho_i*g*hi(i,j)
                wc=MAX(zeta_avg(i,j+1)                                   &
     &                -MAX(-h(i,j+1),                                   &
     &                     -rho_i/rho0*hi(i,j)),0.0_r8)
                Pwat=0.5_r8*rho0*g*wc*wc/MAX(hi(i,j),1.0_r8)
                dp=hi(i,j)*(Pice-Pwat)
                fdy(i,j)=fdy(i,j)+dp*pn(i,j)
              END IF
            END IF
          END IF
        END DO
      END DO
!
#  ifdef ISF_DEBUG
      cff=0.0_r8
      DO j=MAX(LBj,Jstr-1),MIN(UBj,Jend)
        DO i=MAX(LBi,Istr-1),MIN(UBi,Iend)
          cff=MAX(cff,ABS(fdx(i,j)))
        END DO
      END DO
      IF (Master) WRITE (*,*) 'FD-DRIVE probe: max|fdx| =', cff,        &
     &    ' at iic=', iic(ng)
#  endif
!
      RETURN
      END SUBROUTINE ssa_fd_driving
!
!***********************************************************************
      SUBROUTINE ssa_fd_matvec (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Istr, Iend, Jstr, Jend,                 &
     &                          pm, pn, hi,                             &
     &                          nuhx, nuhy, beta,                       &
     &                          xu, xv,                                 &
     &                          Au, Av)
!***********************************************************************
!
!  FD stencil action on the RHO-point velocity pair (PISM SSAFD
!  fd_operator, Maxima coefficients verbatim).  The FD system lives
!  entirely in the RHO space:  no PSI mapping here (the checkerboard null
!  space of the averaged-projection form made the shared PCG
!  diverge on real geometries;  the FD Picard driver ssa_fd_picard
!  runs its own PCG on the cell-centered state instead).
!
!  Marginal switch weights:  an ice-free direct neighbor zeroes that
!  face;  an ice-free diagonal neighbor zeroes the two adjacent
!  interface switches (centered -> one-sided, PISM cfbc rule).
!  Cell-centered beta enters the diagonals (diag_u=eq1(5),
!  diag_v=eq2(14)).  Zero-diagonal cells fall back to identity
!  (PISM replace_zero_diagonal_entries).
!
      implicit none
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Istr, Iend, Jstr, Jend
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)    :: pm(LBi:,LBj:)
      real(r8), intent(in)    :: pn(LBi:,LBj:)
      real(r8), intent(in)    :: hi(LBi:,LBj:)
      real(r8), intent(in)    :: nuhx(LBi:,LBj:)
      real(r8), intent(in)    :: nuhy(LBi:,LBj:)
      real(r8), intent(in)    :: beta(LBi:,LBj:)
      real(r8), intent(in)    :: xu(LBi:,LBj:)
      real(r8), intent(in)    :: xv(LBi:,LBj:)
      real(r8), intent(out)   :: Au(LBi:,LBj:)
      real(r8), intent(out)   :: Av(LBi:,LBj:)
# else
      real(r8), intent(in)    :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: hi(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: nuhx(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: nuhy(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: beta(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: xu(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: xv(LBi:UBi,LBj:UBj)
      real(r8), intent(out)   :: Au(LBi:UBi,LBj:UBj)
      real(r8), intent(out)   :: Av(LBi:UBi,LBj:UBj)
# endif
!
      integer :: i, j, kk, ic, jc
      integer :: Nf, Ef, Sf, Wf
      integer :: NNE, NNW, ENE, ESE, SSW, SSE, WSW, WNW
!
      real(r8) :: cw, ce, cn, cs
      real(r8) :: dx, dy, dx2, dy2, d4, d2
      real(r8) :: bu, bv
!
      real(r8), dimension(18) :: eq1
      real(r8), dimension(18) :: eq2
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (hi(i,j).gt.0.0_r8) THEN
!
            IF (i.gt.LBi) THEN
              cw=nuhx(i-1,j)
            ELSE
              cw=0.0_r8
            END IF
            IF (i.lt.UBi) THEN
              ce=nuhx(i,j)
            ELSE
              ce=0.0_r8
            END IF
            IF (j.lt.UBj) THEN
              cn=nuhy(i,j)
            ELSE
              cn=0.0_r8
            END IF
            IF (j.gt.LBj) THEN
              cs=nuhy(i,j-1)
            ELSE
              cs=0.0_r8
            END IF
!
            Nf=1; Ef=1; Sf=1; Wf=1
            NNE=1; NNW=1; ENE=1; ESE=1; SSW=1; SSE=1; WSW=1; WNW=1
            IF (j.ge.UBj.or.hi(i,j+1).le.0.0_r8) THEN
              Nf=0; NNE=0; NNW=0
            END IF
            IF (j.le.LBj.or.hi(i,j-1).le.0.0_r8) THEN
              Sf=0; SSW=0; SSE=0
            END IF
            IF (i.ge.UBi.or.hi(i+1,j).le.0.0_r8) THEN
              Ef=0; ENE=0; ESE=0
            END IF
            IF (i.le.LBi.or.hi(i-1,j).le.0.0_r8) THEN
              Wf=0; WSW=0; WNW=0
            END IF
            IF (i.lt.UBi.and.j.lt.UBj) THEN
              IF (hi(i+1,j+1).le.0.0_r8) THEN
                NNE=0; ENE=0
              END IF
            END IF
            IF (i.lt.UBi.and.j.gt.LBj) THEN
              IF (hi(i+1,j-1).le.0.0_r8) THEN
                ESE=0; SSE=0
              END IF
            END IF
            IF (i.gt.LBi.and.j.gt.LBj) THEN
              IF (hi(i-1,j-1).le.0.0_r8) THEN
                SSW=0; WSW=0
              END IF
            END IF
            IF (i.gt.LBi.and.j.lt.UBj) THEN
              IF (hi(i-1,j+1).le.0.0_r8) THEN
                WNW=0; NNW=0
              END IF
            END IF
!
            dx=1.0_r8/pm(i,j)
            dy=1.0_r8/pn(i,j)
            dx2=dx*dx
            dy2=dy*dy
            d4=4.0_r8*dx*dy
            d2=2.0_r8*dx*dy
!
!  ---- u-equation (eq1), PISM Maxima coefficients ----
!
            eq1(1)=0.0_r8
            eq1(2)=-cn*Nf/dy2
            eq1(3)=0.0_r8
            eq1(4)=-4.0_r8*cw*Wf/dx2
            eq1(5)=(cn*Nf+cs*Sf)/dy2 + 4.0_r8*(ce*Ef+cw*Wf)/dx2
            eq1(6)=-4.0_r8*ce*Ef/dx2
            eq1(7)=0.0_r8
            eq1(8)=-cs*Sf/dy2
            eq1(9)=0.0_r8
            eq1(10)= cw*Wf*WNW/d2 + cn*Nf*NNW/d4
            eq1(11)=(cn*Nf*NNE-cn*Nf*NNW)*Nf/d4                           &
     &             +(cw*Wf-ce*Ef)*Nf/d2
            eq1(12)=-ce*Ef*ENE/d2 - cn*Nf*NNE/d4
            eq1(13)=(cw*Wf*WSW-cw*Wf*WNW)*Wf/d2                           &
     &             +(cn*Nf-cs*Sf)*Wf/d4
            eq1(14)=(cn*Nf*Ef*Nf-cn*Nf*Wf*Nf                             &
     &              -cs*Sf*Ef*Sf+cs*Sf*Wf*Sf)/d4                         &
     &             +(ce*Ef*Nf-cw*Wf*Nf-ce*Ef*Sf+cw*Wf*Sf)/d2
            eq1(15)=(ce*Ef*ENE-ce*Ef*ESE)*Ef/d2                           &
     &             +(cs*Sf*Ef*Sf-cn*Nf*Ef*Nf)/d4
            eq1(16)=-cw*Wf*WSW/d2 - cs*Sf*SSW/d4
            eq1(17)=(cs*Sf*SSW-cs*Sf*SSE)*Sf/d4                           &
     &             +(ce*Ef-cw*Wf)*Sf/d2
            eq1(18)= ce*Ef*ESE/d2 + cs*Sf*SSE/d4
!
!  ---- v-equation (eq2) ----
!
            eq2(1)= cw*Wf*WNW/d4 + cn*Nf*NNW/d2
            eq2(2)=(cn*Nf*NNE-cn*Nf*NNW)*Nf/d2                            &
     &            +(cw*Wf-ce*Ef)*Nf/d4
            eq2(3)=-ce*Ef*ENE/d4 - cn*Nf*NNE/d2
            eq2(4)=(cw*Wf*WSW-cw*Wf*WNW)*Wf/d4                            &
     &            +(cn*Nf-cs*Sf)*Wf/d2
!
!  PISM SSAFD eq2 center cross term:  the first factor is the
!  (E-W) SWITCH difference, not the face-coefficient difference.
!  The transcription below had (c.e*E-c.w*W) there, producing a
!  spurious term quadratic in nuH (O(1e16) entries across the GL
!  band) -- the source of the FD-PCG divergence and of the LU
!  conditioning catastrophe on realistic nuH fields.  Invisible
!  with uniform nuH (c.e=c.w), which is why the analytic tests
!  never caught it.
!
            eq2(5)=(Ef-Wf)*(cn*Nf-cs*Sf)/d2                               &
     &            +(ce*Ef-cw*Wf)*(Nf-Sf)/d4
            eq2(6)=(ce*Ef*ENE-ce*Ef*ESE)*Ef/d4                            &
     &            +(cs*Sf*Ef*Sf-cn*Nf*Ef*Nf)/d2
            eq2(7)=-cw*Wf*WSW/d4 - cs*Sf*SSW/d2
            eq2(8)=(cs*Sf*SSW-cs*Sf*SSE)*Sf/d2                            &
     &            +(ce*Ef-cw*Wf)*Sf/d4
            eq2(9)= ce*Ef*ESE/d4 + cs*Sf*SSE/d2
            eq2(10)=0.0_r8
            eq2(11)=-4.0_r8*cn*Nf/dy2
            eq2(12)=0.0_r8
            eq2(13)=-cw*Wf/dx2
            eq2(14)=4.0_r8*(cn*Nf+cs*Sf)/dy2 + (ce*Ef+cw*Wf)/dx2
            eq2(15)=-ce*Ef/dx2
            eq2(16)=0.0_r8
            eq2(17)=-4.0_r8*cs*Sf/dy2
            eq2(18)=0.0_r8
!
            bu=beta(i,j)
            bv=beta(i,j)
            eq1(5)=eq1(5)+bu
            eq2(14)=eq2(14)+bv
!
            Au(i,j)=0.0_r8
            Av(i,j)=0.0_r8
!
!  The 9-point coefficients are ordered NW,N,NE,W,C,E,SW,S,SE
!  (PISM SSAFD J-array convention, J(k)=j+1 for k=1..3 and j-1 for
!  k=7..9), so the j offset must DESCEND with kk:  jc=1-(kk-1)/3.
!  The previous jc=(kk-1)/3-1 flipped the S and N rows of the
!  stencil, which makes the assembled matrix NON-symmetric
!  (A[(i,j),(i,j+1)] = -nuh(i,j-1) vs A[(i,j+1),(i,j)] =
!  -nuh(i,j)) -- the root cause of the FD-PCG divergence in the
!  grounded region.  Invisible in y-symmetric test geometries.
!
            DO kk=1,9
              ic=MOD(kk-1,3)-1
              jc=1-(kk-1)/3
              Au(i,j)=Au(i,j)+eq1(kk)*xu(i+ic,j+jc)
              Av(i,j)=Av(i,j)+eq2(kk+9)*xv(i+ic,j+jc)
            END DO
            DO kk=10,18
              ic=MOD(kk-10,3)-1
              jc=1-(kk-10)/3
              Au(i,j)=Au(i,j)+eq1(kk)*xv(i+ic,j+jc)
              Av(i,j)=Av(i,j)+eq2(kk-9)*xu(i+ic,j+jc)
            END DO
!
            IF (ABS(eq1(5)).lt.1.0E-16_r8) Au(i,j)=xu(i,j)
            IF (ABS(eq2(14)).lt.1.0E-16_r8) Av(i,j)=xv(i,j)
          ELSE
!
!  Ice-free cell:  identity row (the FD state is zero there).
!
            Au(i,j)=xu(i,j)
            Av(i,j)=xv(i,j)
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE ssa_fd_matvec
!
      SUBROUTINE ssa_fd_diag (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Istr, Iend, Jstr, Jend,                   &
     &                        IstrP, IendP, JstrP, JendP,               &
     &                        pm, pn,                                   &
     &                        hi, nuhx, nuhy, beta,                     &
     &                        diag_u, diag_v)
!***********************************************************************
!
!  Stencil diagonal (eq1(5) / eq2(14) with the same marginal switch
!  weights) plus the cell-centered beta:  the Jacobi preconditioner
!  of the FD system.  Ice-free cells:  1 (identity rows).
!
      implicit none
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Istr, Iend, Jstr, Jend
      integer, intent(in) :: IstrP, IendP, JstrP, JendP
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: pm(LBi:,LBj:)
      real(r8), intent(in)  :: pn(LBi:,LBj:)
      real(r8), intent(in)  :: hi(LBi:,LBj:)
      real(r8), intent(in)  :: nuhx(LBi:,LBj:)
      real(r8), intent(in)  :: nuhy(LBi:,LBj:)
      real(r8), intent(in)  :: beta(LBi:,LBj:)
      real(r8), intent(out) :: diag_u(LBi:,LBj:)
      real(r8), intent(out) :: diag_v(LBi:,LBj:)
# else
      real(r8), intent(in)  :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: hi(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: nuhx(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: nuhy(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: beta(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: diag_u(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: diag_v(LBi:UBi,LBj:UBj)
# endif
!
      integer :: i, j, nb
!
      real(r8) :: cw, ce, cn, cs
      real(r8) :: dx, dy
      real(r8), dimension(LBi:UBi,LBj:UBj) :: wrkx
      real(r8), dimension(LBi:UBi,LBj:UBj) :: wrky
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (hi(i,j).gt.0.0_r8) THEN
            IF (i.gt.LBi.and.hi(i-1,j).gt.0.0_r8) THEN
              cw=nuhx(i-1,j)
            ELSE
              cw=0.0_r8
            END IF
            IF (i.lt.UBi.and.hi(i+1,j).gt.0.0_r8) THEN
              ce=nuhx(i,j)
            ELSE
              ce=0.0_r8
            END IF
            IF (j.lt.UBj.and.hi(i,j+1).gt.0.0_r8) THEN
              cn=nuhy(i,j)
            ELSE
              cn=0.0_r8
            END IF
            IF (j.gt.LBj.and.hi(i,j-1).gt.0.0_r8) THEN
              cs=nuhy(i,j-1)
            ELSE
              cs=0.0_r8
            END IF
            dx=1.0_r8/pm(i,j)
            dy=1.0_r8/pn(i,j)
            wrkx(i,j)=(cn+cs)/dy**2 + 4.0_r8*(cw+ce)/dx**2               &
     &              +beta(i,j)
            wrky(i,j)=4.0_r8*(cn+cs)/dy**2 + (cw+ce)/dx**2               &
     &              +beta(i,j)
            IF (ABS(wrkx(i,j)).lt.1.0E-16_r8) wrkx(i,j)=1.0_r8
            IF (ABS(wrky(i,j)).lt.1.0E-16_r8) wrky(i,j)=1.0_r8
          ELSE
            wrkx(i,j)=1.0_r8
            wrky(i,j)=1.0_r8
          END IF
        END DO
      END DO
!
!  Transpose map to the PSI nodes (Jacobi diagonals are consumed at
!  PSI nodes by the shared PCG driver).
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          nb=0
          diag_u(i,j)=0.0_r8
          diag_v(i,j)=0.0_r8
          IF (i-1.ge.Istr.and.j-1.ge.Jstr) THEN
            diag_u(i,j)=diag_u(i,j)+wrkx(i-1,j-1)
            diag_v(i,j)=diag_v(i,j)+wrky(i-1,j-1)
            nb=nb+1
          END IF
          IF (i.le.Iend.and.j-1.ge.Jstr) THEN
            diag_u(i,j)=diag_u(i,j)+wrkx(i,j-1)
            diag_v(i,j)=diag_v(i,j)+wrky(i,j-1)
            nb=nb+1
          END IF
          IF (i-1.ge.Istr.and.j.le.Jend) THEN
            diag_u(i,j)=diag_u(i,j)+wrkx(i-1,j)
            diag_v(i,j)=diag_v(i,j)+wrky(i-1,j)
            nb=nb+1
          END IF
          IF (i.le.Iend.and.j.le.Jend) THEN
            diag_u(i,j)=diag_u(i,j)+wrkx(i,j)
            diag_v(i,j)=diag_v(i,j)+wrky(i,j)
            nb=nb+1
          END IF
          IF (nb.gt.0) THEN
            diag_u(i,j)=diag_u(i,j)/REAL(nb,r8)
            diag_v(i,j)=diag_v(i,j)/REAL(nb,r8)
          ELSE
            diag_u(i,j)=1.0_r8
            diag_v(i,j)=1.0_r8
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE ssa_fd_diag
!
!***********************************************************************
      SUBROUTINE ssa_fd_picard (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Istr, Iend, Jstr, Jend,                 &
     &                          pm, pn, hi,                             &
     &                          nuhx, nuhy, beta,                       &
     &                          fdx, fdy, B_hard,                       &
     &                          urho, vrho,                             &
     &                          Ru, Rv, Zu, Zv, Du, Dv,                 &
     &                          Au, Av, diag_u, diag_v)
!***********************************************************************
!
!  FD Picard/PCG driver (PISM SSAFD semantics, RHO space).  One-for-
!  one counterpart of the FE Picard block in isf_ssa_tile Stage 6:
!  same dual convergence criterion (relative velocity change +
!  relative residual), same Jacobi-preconditioned Fletcher-Reeves
!  PCG, same speed-cap damper and Picard relaxation -- every vector
!  and operator in the cell-centered space.
!
      implicit none
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Istr, Iend, Jstr, Jend
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)    :: pm(LBi:,LBj:)
      real(r8), intent(in)    :: pn(LBi:,LBj:)
      real(r8), intent(in)    :: hi(LBi:,LBj:)
      real(r8), intent(inout) :: nuhx(LBi:,LBj:)
      real(r8), intent(inout) :: nuhy(LBi:,LBj:)
      real(r8), intent(in)    :: beta(LBi:,LBj:)
      real(r8), intent(in)    :: fdx(LBi:,LBj:)
      real(r8), intent(in)    :: fdy(LBi:,LBj:)
      real(r8), intent(in)    :: B_hard(LBi:,LBj:)
      real(r8), intent(inout) :: urho(LBi:,LBj:)
      real(r8), intent(inout) :: vrho(LBi:,LBj:)
      real(r8), intent(inout) :: Ru(LBi:,LBj:)
      real(r8), intent(inout) :: Rv(LBi:,LBj:)
      real(r8), intent(inout) :: Zu(LBi:,LBj:)
      real(r8), intent(inout) :: Zv(LBi:,LBj:)
      real(r8), intent(inout) :: Du(LBi:,LBj:)
      real(r8), intent(inout) :: Dv(LBi:,LBj:)
      real(r8), intent(inout) :: Au(LBi:,LBj:)
      real(r8), intent(inout) :: Av(LBi:,LBj:)
      real(r8), intent(inout) :: diag_u(LBi:,LBj:)
      real(r8), intent(inout) :: diag_v(LBi:,LBj:)
# else
      real(r8), intent(in)    :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: hi(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: nuhx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: nuhy(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: beta(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: fdx(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: fdy(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: B_hard(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: urho(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vrho(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Ru(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Rv(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Zu(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Zv(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Du(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Dv(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Au(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Av(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: diag_u(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: diag_v(LBi:UBi,LBj:UBj)
# endif
!
      integer :: i, j, ii, jj, nl_iter, cg_iter
      integer :: cw, ce, cs, cn
      logical :: converged
!
      real(r8) :: cgtol, eps, alpha, beta_cg, cff
      real(r8) :: dot_zr, dot_zr_old, dot_dad
      real(r8) :: resid, resid_0, nl_resid
      real(r8) :: u_slide, u_prev, u_allow
      real(r8) :: vel_norm, vel_change, rhs_norm, resid_norm
      real(r8) :: uo_n, un_n
      real(r8) :: dx2, dy2
      real(r8) :: norms(4)
      character (len=3) :: op_norms(4), op_handle
!
      eps=1.0E-16_r8
      cgtol=cg_tol(ng)
      converged=.FALSE.
!
!  Anchored ring:  the outermost computational cells carry identity
!  rows (velocity held at the mapped-in value).  This anchors the
!  rigid-translation null space of the beta=0 floating ice the same
!  way the FE path is anchored by ssa_lbc_zero at the domain edges
!  (and matches PISM, where the ice-sheet interior supplies the
!  reference frame).
!
      DO j=LBj,UBj
        DO i=LBi,UBi
          IF ((i.eq.Istr.or.i.eq.Iend.or.j.eq.Jstr.or.j.eq.Jend).and.   &
     &        hi(i,j).gt.0.0_r8) THEN
            urho(i,j)=0.0_r8
            vrho(i,j)=0.0_r8
          END IF
        END DO
      END DO
!
      DO nl_iter=1,nl_max_iter(ng)
!
        CALL ssa_fd_viscosity (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         pm, pn, hi, B_hard,                      &
     &                         urho, vrho, nuhx, nuhy)
!
!  RHO-space Jacobi diagonal (eq1(5)/eq2(14) + beta):  the FD PCG
!  runs in the cell space, so no PSI mapping here (the mapped
!  diagonal belongs to the shared-driver variant).
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            IF (hi(i,j).gt.0.0_r8) THEN
              IF (i.gt.LBi.and.hi(i-1,j).gt.0.0_r8) THEN
                cw=nuhx(i-1,j)
              ELSE
                cw=0.0_r8
              END IF
              IF (i.lt.UBi.and.hi(i+1,j).gt.0.0_r8) THEN
                ce=nuhx(i,j)
              ELSE
                ce=0.0_r8
              END IF
              IF (j.lt.UBj.and.hi(i,j+1).gt.0.0_r8) THEN
                cn=nuhy(i,j)
              ELSE
                cn=0.0_r8
              END IF
              IF (j.gt.LBj.and.hi(i,j-1).gt.0.0_r8) THEN
                cs=nuhy(i,j-1)
              ELSE
                cs=0.0_r8
              END IF
              dx2=(1.0_r8/pm(i,j))**2
              dy2=(1.0_r8/pn(i,j))**2
              diag_u(i,j)=(cn+cs)/dy2+4.0_r8*(cw+ce)/dx2+beta(i,j)
              diag_v(i,j)=4.0_r8*(cn+cs)/dy2+(cw+ce)/dx2+beta(i,j)
              IF (ABS(diag_u(i,j)).lt.eps) diag_u(i,j)=1.0_r8
              IF (ABS(diag_v(i,j)).lt.eps) diag_v(i,j)=1.0_r8
            ELSE
              diag_u(i,j)=1.0_r8
              diag_v(i,j)=1.0_r8
            END IF
          END DO
        END DO
!
        CALL ssa_fd_matvec (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Istr, Iend, Jstr, Jend,                     &
     &                      pm, pn, hi, nuhx, nuhy, beta,               &
     &                      urho, vrho, Au, Av)
!
        rhs_norm=0.0_r8
        resid_norm=0.0_r8
        DO j=Jstr,Jend
          DO i=Istr,Iend
            IF (hi(i,j).gt.0.0_r8.and.                                &
     &                 .not.(i.eq.Istr.or.i.eq.Iend.or.                     &
     &                      j.eq.Jstr.or.j.eq.Jend)) THEN
              Ru(i,j)=fdx(i,j)-Au(i,j)
              Rv(i,j)=fdy(i,j)-Av(i,j)
              rhs_norm=rhs_norm+fdx(i,j)*fdx(i,j)+fdy(i,j)*fdy(i,j)
              resid_norm=resid_norm+Ru(i,j)*Ru(i,j)+Rv(i,j)*Rv(i,j)
            ELSE
              Ru(i,j)=0.0_r8
              Rv(i,j)=0.0_r8
            END IF
          END DO
        END DO
!
        dot_zr=0.0_r8
        DO j=Jstr,Jend
          DO i=Istr,Iend
            IF (hi(i,j).gt.0.0_r8.and.                                &
     &                 .not.(i.eq.Istr.or.i.eq.Iend.or.                     &
     &                      j.eq.Jstr.or.j.eq.Jend)) THEN
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
# ifdef DISTRIBUTE
        op_handle='SUM'
        CALL mp_reduce (ng, model, 1, dot_zr, op_handle)
# endif
        resid_0=SQRT(ABS(dot_zr))
        resid=resid_0
        dot_zr_old=dot_zr
!
#  ifdef ISF_DEBUG
        cff=1.0E30_r8
        DO jj=Jstr,Jend
          DO ii=Istr,Iend
            IF (hi(ii,jj).gt.0.0_r8) cff=MIN(cff,diag_u(ii,jj))
          END DO
        END DO
        IF (Master) WRITE (*,*) 'FD-PCG probe: rhs_norm=',             &
     &      SQRT(rhs_norm),' resid_0=',resid_0,' diag_min=',cff,       &
     &      ' nuhx_max=',MAXVAL(nuhx)
#  endif
!
        DO cg_iter=1,cg_max_iter(ng)
          IF (resid.le.cgtol*resid_0 .or. resid.lt.eps) EXIT
!
          CALL ssa_fd_matvec (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Istr, Iend, Jstr, Jend,                   &
     &                        pm, pn, hi, nuhx, nuhy, beta,             &
     &                        Du, Dv, Au, Av)
!
          dot_dad=0.0_r8
          DO j=Jstr,Jend
            DO i=Istr,Iend
              IF (hi(i,j).gt.0.0_r8.and.                                &
     &                 .not.(i.eq.Istr.or.i.eq.Iend.or.                     &
     &                      j.eq.Jstr.or.j.eq.Jend)) THEN
                dot_dad=dot_dad+Du(i,j)*Au(i,j)+Dv(i,j)*Av(i,j)
              END IF
            END DO
          END DO
!
# ifdef DISTRIBUTE
          CALL mp_reduce (ng, model, 1, dot_dad, op_handle)
# endif
!
          IF (ABS(dot_dad).lt.eps) THEN
            IF (Master) WRITE(stdout,'(a,i3,a,1p,e12.4)')               &
     &          ' SSA-FD CG stagnation: iter=',cg_iter,' dot_dad=',dot_dad
            EXIT
          END IF
          alpha=dot_zr/dot_dad
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
              IF (hi(i,j).gt.0.0_r8.and.                                &
     &                 .not.(i.eq.Istr.or.i.eq.Iend.or.                     &
     &                      j.eq.Jstr.or.j.eq.Jend)) THEN
                urho(i,j)=urho(i,j)+alpha*Du(i,j)
                vrho(i,j)=vrho(i,j)+alpha*Dv(i,j)
                Ru(i,j)=Ru(i,j)-alpha*Au(i,j)
                Rv(i,j)=Rv(i,j)-alpha*Av(i,j)
              END IF
            END DO
          END DO
!
          dot_zr=0.0_r8
          DO j=Jstr,Jend
            DO i=Istr,Iend
              IF (hi(i,j).gt.0.0_r8.and.                                &
     &                 .not.(i.eq.Istr.or.i.eq.Iend.or.                     &
     &                      j.eq.Jstr.or.j.eq.Jend)) THEN
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
          CALL mp_reduce (ng, model, 1, dot_zr, op_handle)
# endif
!
          IF (ABS(dot_zr_old).lt.eps) EXIT
          beta_cg=dot_zr/dot_zr_old
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
              IF (hi(i,j).gt.0.0_r8.and.                                &
     &                 .not.(i.eq.Istr.or.i.eq.Iend.or.                     &
     &                      j.eq.Jstr.or.j.eq.Jend)) THEN
                Du(i,j)=Zu(i,j)+beta_cg*Du(i,j)
                Dv(i,j)=Zv(i,j)+beta_cg*Dv(i,j)
              END IF
            END DO
          END DO
!
          resid=SQRT(ABS(dot_zr))
          dot_zr_old=dot_zr
!
        END DO  ! cg_iter
!
!  Speed-cap damper on the cell-centered state (FE semantics).
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            IF (hi(i,j).gt.0.0_r8.and.                                &
     &                 .not.(i.eq.Istr.or.i.eq.Iend.or.                     &
     &                      j.eq.Jstr.or.j.eq.Jend)) THEN
              u_slide=SQRT(urho(i,j)**2+vrho(i,j)**2)
              u_allow=1.5_r8*ssa_max_speed(ng)
              IF (u_slide.gt.u_allow) THEN
                u_prev=SQRT(MAX(urho(i,j)-alpha*Du(i,j),0.0_r8)**2+      &
     &                      MAX(vrho(i,j)-alpha*Dv(i,j),0.0_r8)**2)
                IF (u_prev.gt.u_allow) THEN
                  u_allow=MAX(u_allow,0.9_r8*u_prev)
                END IF
                IF (u_slide.gt.u_allow) THEN
                  cff=u_allow/u_slide
                  urho(i,j)=cff*urho(i,j)
                  vrho(i,j)=cff*vrho(i,j)
                END IF
              END IF
            END IF
          END DO
        END DO
!
!  Picard relaxation:  u = alpha*u_new + (1-alpha)*u_old.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            IF (hi(i,j).gt.0.0_r8.and.                                &
     &                 .not.(i.eq.Istr.or.i.eq.Iend.or.                     &
     &                      j.eq.Jstr.or.j.eq.Jend)) THEN
              urho(i,j)=picard_alpha(ng)*urho(i,j)                       &
     &                 +(1.0_r8-picard_alpha(ng))                        &
     &                  *(urho(i,j)-alpha*Du(i,j))
              vrho(i,j)=picard_alpha(ng)*vrho(i,j)                       &
     &                 +(1.0_r8-picard_alpha(ng))                        &
     &                  *(vrho(i,j)-alpha*Dv(i,j))
            END IF
          END DO
        END DO
!
!  Dual convergence:  relative velocity change and relative residual.
!
        uo_n=0.0_r8
        un_n=0.0_r8
        DO j=Jstr,Jend
          DO i=Istr,Iend
            IF (hi(i,j).gt.0.0_r8.and.                                &
     &                 .not.(i.eq.Istr.or.i.eq.Iend.or.                     &
     &                      j.eq.Jstr.or.j.eq.Jend)) THEN
              uo_n=uo_n+MAX(urho(i,j)-alpha*Du(i,j),0.0_r8)**2           &
     &                  +MAX(vrho(i,j)-alpha*Dv(i,j),0.0_r8)**2
              un_n=un_n+urho(i,j)**2+vrho(i,j)**2
            END IF
          END DO
        END DO
# ifdef DISTRIBUTE
        norms(1)=uo_n
        norms(2)=un_n
        norms(3)=resid_norm
        norms(4)=rhs_norm
        op_norms(1)='SUM'
        op_norms(2)='SUM'
        op_norms(3)='SUM'
        op_norms(4)='SUM'
        CALL mp_reduce (ng, model, 4, norms, op_norms)
        uo_n=norms(1)
        un_n=norms(2)
        resid_norm=norms(3)
        rhs_norm=norms(4)
# endif
        vel_norm=SQRT(un_n)
        vel_change=SQRT(ABS(un_n-uo_n))
        IF (vel_norm.gt.eps) THEN
          nl_resid=vel_change/vel_norm
        ELSE
          nl_resid=0.0_r8
        END IF
        IF (rhs_norm.gt.eps) THEN
          resid_norm=SQRT(resid_norm)/SQRT(rhs_norm)
          nl_resid=MAX(nl_resid,resid_norm)
        END IF
!
        IF (nl_resid.lt.nl_tol(ng)) converged=.TRUE.
!
        IF (Master) WRITE(stdout,'(a,i3,a,1p,e9.2,a,e9.2,a,l1)')        &
     &      ' SSA-FD Picard iter',nl_iter,' nl_resid=',nl_resid,        &
     &      ' max|u|=',SQRT(un_n),' conv=',converged
!
        IF (converged) EXIT
!
      END DO  ! nl_iter
!
      IF (.not.converged .and. Master) THEN
        WRITE(stdout,'(a,i4,a,es9.2)')                                  &
     &      ' SSA-FD WARNING: Picard not converged in ',                &
     &      nl_max_iter(ng),' iterations; last nl_resid=',nl_resid
      END IF
!
      RETURN
      END SUBROUTINE ssa_fd_picard
!
# ifdef ICESHEET_PETSC
      SUBROUTINE ssa_fd_picard_ksp (ng, tile, model,                     &
     &                          LBi, UBi, LBj, UBj,                      &
     &                          Istr, Iend, Jstr, Jend,                  &
     &                          pm, pn, hi,                              &
     &                          nuhx, nuhy, beta,                        &
     &                          fdx, fdy, B_hard,                        &
     &                          urho, vrho,                              &
     &                          Ru, Rv, Zu, Zv, Du, Dv,                  &
     &                          Au, Av, diag_u, diag_v)
!***********************************************************************
!
!  FD Picard driver with the linear solve routed through PETSc KSP
!  (ICESHEET_PETSC, fd_ksp=1).  Counterpart of ssa_fd_picard:  the
!  Picard loop, viscosity, speed-cap damper, relaxation and dual
!  convergence criterion are the same;  the built-in Jacobi PCG is
!  replaced by an explicitly assembled global matrix solved on rank
!  0 (SeqAIJ, KSP preonly + LU direct by default;  -ksp_type and
!  -pc_type command line flags override).  This removes the PCG
!  divergence seen on the wall-corner MISMIP+ geometry and, as a
!  side effect, is also cleaner in MPI:  the matrix is assembled
!  for the full domain with anchoring only on the LBC-flagged
!  edges, so tile seams carry no artificial zero-velocity pinning,
!  unlike the per-tile anchoring of the PCG path.
!
!  Data flow:  tile interiors are gathered to rank 0 (fd_gather),
!  the matrix/RHS are assembled from the same PISM Maxima stencil
!  coefficients as ssa_fd_matvec (eq1/eq2 transcribed verbatim),
!  the solved, damped and relaxed cell state is scattered back to
!  every tile interior (fd_scatter) each Picard iteration.  The
!  anchoring is LBC-aware and per component (like ssa_lbc_zero on
!  the FE path):  domain edges flagged closed/clamped in the
!  ice-state LBC carry identity rows with zero right-hand sides,
!  and one reference ice cell is pinned when a component has no
!  Dirichlet edge at all (rigid-translation null space of a free
!  component under a direct solve).  Unconditional edge anchoring
!  would pin free (gradient) edges and crush y-uniform flow through
!  the membrane stiffness -- the 250x front-test discrepancy.
!  Cells with a vanishing diagonal fall back to identity rows
!  driven by the local RHS, matching the replace_zero_diagonal
!  semantics of the matvec.
!
      implicit none
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Istr, Iend, Jstr, Jend
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in)    :: pm(LBi:,LBj:)
      real(r8), intent(in)    :: pn(LBi:,LBj:)
      real(r8), intent(in)    :: hi(LBi:,LBj:)
      real(r8), intent(inout) :: nuhx(LBi:,LBj:)
      real(r8), intent(inout) :: nuhy(LBi:,LBj:)
      real(r8), intent(in)    :: beta(LBi:,LBj:)
      real(r8), intent(in)    :: fdx(LBi:,LBj:)
      real(r8), intent(in)    :: fdy(LBi:,LBj:)
      real(r8), intent(in)    :: B_hard(LBi:,LBj:)
      real(r8), intent(inout) :: urho(LBi:,LBj:)
      real(r8), intent(inout) :: vrho(LBi:,LBj:)
      real(r8), intent(inout) :: Ru(LBi:,LBj:)
      real(r8), intent(inout) :: Rv(LBi:,LBj:)
      real(r8), intent(inout) :: Zu(LBi:,LBj:)
      real(r8), intent(inout) :: Zv(LBi:,LBj:)
      real(r8), intent(inout) :: Du(LBi:,LBj:)
      real(r8), intent(inout) :: Dv(LBi:,LBj:)
      real(r8), intent(inout) :: Au(LBi:,LBj:)
      real(r8), intent(inout) :: Av(LBi:,LBj:)
      real(r8), intent(inout) :: diag_u(LBi:,LBj:)
      real(r8), intent(inout) :: diag_v(LBi:,LBj:)
#  else
      real(r8), intent(in)    :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: hi(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: nuhx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: nuhy(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: beta(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: fdx(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: fdy(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: B_hard(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: urho(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vrho(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Ru(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Rv(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Zu(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Zv(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Du(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Dv(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Au(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Av(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: diag_u(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: diag_v(LBi:UBi,LBj:UBj)
#  endif
!
!  Local variable declarations.
!
#  ifdef DISTRIBUTE
      include 'mpif.h'
#  endif
!
      logical :: converged
!
      logical :: ancu_w, ancu_e, ancu_s, ancu_n
      logical :: ancv_w, ancv_e, ancv_s, ancv_n
      logical :: hasu, hasv
!
      integer :: i, j, k, kk, ic, jc, nl_iter
      integer :: Lm_g, Mm_g, nslab, ncell2
      integer :: gi, gn
      integer :: k2
      integer :: fi_i, fi_j
      integer :: iu
      integer :: ii, jj
      integer :: Nf, Ef, Sf, Wf
      integer :: NNE, NNW, ENE, ESE, SSW, SSE, WSW, WNW
      integer :: iconv
!
!  mrank/nrpc/mpi_err are declared unconditionally:  the serial
!  branch (mrank=0) and the rank-0 PETSc object blocks reference
!  them outside the DISTRIBUTE guards.  nrpc/mpi_err stay unused
!  in serial builds.
!
      integer :: nrpc, mrank, mpi_err
#  ifdef DISTRIBUTE
      integer :: myrect(4)
      integer, allocatable :: rects(:)
      integer, allocatable :: counts(:)
      integer, allocatable :: displs(:)
#  endif
!
      real(r8) :: eps, cff
      real(r8) :: rm_row
      real(r8) :: u_slide, u_prev, u_allow
      real(r8) :: vel_norm, vel_change
      real(r8) :: nl_resid
      real(r8) :: uo_n, un_n
      real(r8) :: cwv, cev, cnv, csv
      real(r8) :: dx, dy, dx2, dy2, d4, d2
      real(r8) :: bu, bv
      real(r8), dimension(18) :: eq1
      real(r8), dimension(18) :: eq2
!
      real(r8), allocatable, dimension(:,:) :: g_hi
      real(r8), allocatable, dimension(:,:) :: g_beta
      real(r8), allocatable, dimension(:,:) :: g_fdx
      real(r8), allocatable, dimension(:,:) :: g_fdy
      real(r8), allocatable, dimension(:,:) :: g_pm
      real(r8), allocatable, dimension(:,:) :: g_pn
      real(r8), allocatable, dimension(:,:) :: g_nuhx
      real(r8), allocatable, dimension(:,:) :: g_nuhy
      real(r8), allocatable, dimension(:,:) :: g_uold
      real(r8), allocatable, dimension(:,:) :: g_vold
      real(r8), allocatable, dimension(:,:) :: g_unew
      real(r8), allocatable, dimension(:,:) :: g_vnew
      real(r8), allocatable, dimension(:) :: sbuf
      real(r8), allocatable, dimension(:) :: rscale
      real(r8), allocatable, dimension(:) :: rbuf
      logical, allocatable, dimension(:,:) :: pin_u
      logical, allocatable, dimension(:,:) :: pin_v
!
!  PETSc objects and scalars.
!
      Mat  :: pA
      Vec  :: pb, px, pw
      KSP  :: pksp
      PC   :: ppc
      PetscInt :: ndof, i1(1), icol(18)
      PetscScalar :: vcol(18), one
      PetscScalar, pointer :: xv(:)
      PetscErrorCode :: perr
      PetscReal :: resn, rhsn
      KSPConvergedReason :: preason
!
!-----------------------------------------------------------------------
!  Setup:  global field buffers, MPI tile layout.
!-----------------------------------------------------------------------
!
      eps=1.0E-16_r8
      converged=.FALSE.
      Lm_g=Lm(ng)
      Mm_g=Mm(ng)
      nslab=(Iend-Istr+1)*(Jend-Jstr+1)
      ncell2=2*Lm_g*Mm_g
      one=1.0d0
!
!  Serial (non-MPI) builds:  the KSP path requires a single tile
!  covering the full domain;  otherwise fall back to the internal
!  PCG (which anchors per tile).  MPI builds:  any decomposition.
!
#  ifndef DISTRIBUTE
      IF (Istr.ne.1.or.Iend.ne.Lm_g.or.                           &
     &    Jstr.ne.1.or.Jend.ne.Mm_g) THEN
        IF (Master) WRITE(stdout,'(a)')                            &
     &    ' SSA-FD-KSP: serial multi-tile run; using internal PCG.'
        CALL ssa_fd_picard (ng, tile, model,                       &
     &                      LBi, UBi, LBj, UBj,                    &
     &                      Istr, Iend, Jstr, Jend,                &
     &                      pm, pn, hi, nuhx, nuhy, beta,          &
     &                      fdx, fdy, B_hard, urho, vrho,          &
     &                      Ru, Rv, Zu, Zv, Du, Dv,                &
     &                      Au, Av, diag_u, diag_v)
        RETURN
      END IF
#  endif
!
      allocate ( g_hi(0:Lm_g+1,0:Mm_g+1) )
      allocate ( g_beta(0:Lm_g+1,0:Mm_g+1) )
      allocate ( g_fdx(0:Lm_g+1,0:Mm_g+1) )
      allocate ( g_fdy(0:Lm_g+1,0:Mm_g+1) )
      allocate ( g_pm(0:Lm_g+1,0:Mm_g+1) )
      allocate ( g_pn(0:Lm_g+1,0:Mm_g+1) )
      allocate ( g_nuhx(0:Lm_g+1,0:Mm_g+1) )
      allocate ( g_nuhy(0:Lm_g+1,0:Mm_g+1) )
      allocate ( g_uold(0:Lm_g+1,0:Mm_g+1) )
      allocate ( g_vold(0:Lm_g+1,0:Mm_g+1) )
      allocate ( g_unew(0:Lm_g+1,0:Mm_g+1) )
      allocate ( g_vnew(0:Lm_g+1,0:Mm_g+1) )
      allocate ( sbuf(nslab) )
      allocate ( rbuf(Lm_g*Mm_g) )
      allocate ( pin_u(Lm_g,Mm_g) )
      allocate ( pin_v(Lm_g,Mm_g) )
      allocate ( rscale(ncell2) )
      g_hi=0.0_r8
      g_beta=0.0_r8
      g_fdx=0.0_r8
      g_fdy=0.0_r8
      g_pm=0.0_r8
      g_pn=0.0_r8
      g_nuhx=0.0_r8
      g_nuhy=0.0_r8
      g_uold=0.0_r8
      g_vold=0.0_r8
      g_unew=0.0_r8
      g_vnew=0.0_r8
!
#  ifdef DISTRIBUTE
      mrank=MyRank
      myrect(1)=Istr
      myrect(2)=Iend
      myrect(3)=Jstr
      myrect(4)=Jend
      CALL mpi_comm_size (OCN_COMM_WORLD, nrpc, mpi_err)
      allocate ( rects(4*nrpc) )
      allocate ( counts(nrpc) )
      allocate ( displs(nrpc) )
      CALL mpi_allgather (myrect, 4, MPI_INTEGER,                   &
     &                    rects, 4, MPI_INTEGER,                    &
     &                    OCN_COMM_WORLD, mpi_err)
      displs(1)=0
      DO k=1,nrpc
        counts(k)=(rects(4*(k-1)+2)-rects(4*(k-1)+1)+1)*           &
     &            (rects(4*(k-1)+4)-rects(4*(k-1)+3)+1)
        IF (k.gt.1) displs(k)=displs(k-1)+counts(k-1)
      END DO
#  else
      mrank=0
#  endif
!
!  Gather the solve-static fields once (rank 0 consumes them).
!
      CALL fd_gather (hi, g_hi)
      CALL fd_gather (beta, g_beta)
      CALL fd_gather (fdx, g_fdx)
      CALL fd_gather (fdy, g_fdy)
      CALL fd_gather (pm, g_pm)
      CALL fd_gather (pn, g_pn)
!
!  PETSc lifecycle:  initialize on every rank (collective), keep
!  all objects sequential on rank 0 (PETSC_COMM_SELF), finalize on
!  every rank at the end of the solve.  Init/finalize cycling is
!  supported and keeps this routine self-contained in the ROMS
!  driver (no hooks needed in the ROMS startup/shutdown path).
!
      CALL PetscInitialize (perr)
!
      IF (mrank.eq.0) THEN
        ndof=ncell2
        CALL MatCreateSeqAIJ (PETSC_COMM_SELF, ndof, ndof, 18,     &
     &                        PETSC_NULL_INTEGER_ARRAY, pA, perr)
        CALL VecCreateSeq (PETSC_COMM_SELF, ndof, pb, perr)
        CALL VecDuplicate (pb, px, perr)
        CALL VecDuplicate (pb, pw, perr)
        CALL KSPCreate (PETSC_COMM_SELF, pksp, perr)
        CALL KSPSetType (pksp, KSPPREONLY, perr)
        CALL KSPGetPC (pksp, ppc, perr)
        CALL PCSetType (ppc, PCLU, perr)
        CALL KSPSetFromOptions (pksp, perr)
      END IF
!
!-----------------------------------------------------------------------
!  LBC-aware anchoring (rank 0).  Identity rows are placed on the
!  edges flagged closed/clamped in the ice-state LBC, per velocity
!  component -- the same per-edge/per-component logic as ssa_lbc_zero
!  on the FE path.  Unconditional domain-edge anchoring would pin
!  free (gradient) edges and crush any y-uniform flow through the
!  membrane stiffness (found as the 250x front-test discrepancy).
!  If a component has no closed/clamped edge at all, one reference
!  ice cell is pinned to zero to remove the rigid-translation null
!  space (the direct solve cannot tolerate a singular matrix;  the
!  physical solution has that component zero at the pinned cell).
!-----------------------------------------------------------------------
!
      IF (mrank.eq.0) THEN
        IF (fd_ksp(ng).eq.2) THEN
!
!         Legacy semantics (fd_ksp=2):  unconditional domain-edge
!         anchoring, both components -- the pre-LBC-aware behavior.
!         Use for MISMIP-type diagnostic runs where the LU solve is
!         known to behave;  the free-edge crushing of y-uniform flow
!         (the 250x front-test artifact) is part of this mode.
!
          ancu_w=.TRUE.
          ancu_e=.TRUE.
          ancu_s=.TRUE.
          ancu_n=.TRUE.
          ancv_w=.TRUE.
          ancv_e=.TRUE.
          ancv_s=.TRUE.
          ancv_n=.TRUE.
          hasu=.TRUE.
          hasv=.TRUE.
        ELSE
!
!         LBC-aware anchoring (fd_ksp=1):  read the ice-state LBC.
!
        ancu_w=LBC(iwest ,ibISF(isIuis),ng)%clamped.or.               &
     &         LBC(iwest ,ibISF(isIuis),ng)%closed
        ancu_e=LBC(ieast ,ibISF(isIuis),ng)%clamped.or.               &
     &         LBC(ieast ,ibISF(isIuis),ng)%closed
        ancu_s=LBC(isouth,ibISF(isIuis),ng)%clamped.or.               &
     &         LBC(isouth,ibISF(isIuis),ng)%closed
        ancu_n=LBC(inorth,ibISF(isIuis),ng)%clamped.or.               &
     &         LBC(inorth,ibISF(isIuis),ng)%closed
        ancv_w=LBC(iwest ,ibISF(isIvis),ng)%clamped.or.               &
     &         LBC(iwest ,ibISF(isIvis),ng)%closed
        ancv_e=LBC(ieast ,ibISF(isIvis),ng)%clamped.or.               &
     &         LBC(ieast ,ibISF(isIvis),ng)%closed
        ancv_s=LBC(isouth,ibISF(isIvis),ng)%clamped.or.               &
     &         LBC(isouth,ibISF(isIvis),ng)%closed
        ancv_n=LBC(inorth,ibISF(isIvis),ng)%clamped.or.               &
     &         LBC(inorth,ibISF(isIvis),ng)%closed
        hasu=ancu_w.or.ancu_e.or.ancu_s.or.ancu_n
        hasv=ancv_w.or.ancv_e.or.ancv_s.or.ancv_n
        END IF
!
#  ifdef ISF_DEBUG
        WRITE(stdout,'(a,4i4)') ' FD-KSP LBC slots (u,v):',         &
     &      ibISF(isIuis), ibISF(isIvis)
        WRITE(stdout,'(a,8l2)') ' FD-KSP anchor flags u(w,e,s,n)'//  &
     &      ' v(w,e,s,n):',                                          &
     &      ancu_w,ancu_e,ancu_s,ancu_n,ancv_w,ancv_e,ancv_s,ancv_n
        DO kk=1,2
          IF (kk.eq.1) THEN
            iu=ibISF(isIuis)
          ELSE
            iu=ibISF(isIvis)
          END IF
          WRITE(stdout,'(a,i3,4(a,2l2))') ' FD-KSP slot',iu,        &
     &       ' W(c,g)',LBC(iwest,iu,ng)%closed,                     &
     &                   LBC(iwest,iu,ng)%gradient,                 &
     &       ' E(c,g)',LBC(ieast,iu,ng)%closed,                     &
     &                   LBC(ieast,iu,ng)%gradient,                 &
     &       ' S(c,g)',LBC(isouth,iu,ng)%closed,                    &
     &                   LBC(isouth,iu,ng)%gradient,                &
     &       ' N(c,g)',LBC(inorth,iu,ng)%closed,                    &
     &                   LBC(inorth,iu,ng)%gradient
        END DO
#  endif
!
!  Reference cell for a component without any Dirichlet edge:  the
!  first ice cell in scan order (deterministic).
!
        fi_i=0
        fi_j=0
        DO j=1,Mm_g
          DO i=1,Lm_g
            IF (fi_i.eq.0.and.g_hi(i,j).gt.0.0_r8) THEN
              fi_i=i
              fi_j=j
            END IF
          END DO
        END DO
        DO j=1,Mm_g
          DO i=1,Lm_g
            pin_u(i,j)=(g_hi(i,j).le.0.0_r8).or.                      &
     &                 (i.eq.1    .and.ancu_w).or.                    &
     &                 (i.eq.Lm_g.and.ancu_e).or.                    &
     &                 (j.eq.1    .and.ancu_s).or.                    &
     &                 (j.eq.Mm_g.and.ancu_n).or.                    &
     &                 (.not.hasu.and.i.eq.fi_i.and.j.eq.fi_j)
            pin_v(i,j)=(g_hi(i,j).le.0.0_r8).or.                      &
     &                 (i.eq.1    .and.ancv_w).or.                    &
     &                 (i.eq.Lm_g.and.ancv_e).or.                    &
     &                 (j.eq.1    .and.ancv_s).or.                    &
     &                 (j.eq.Mm_g.and.ancv_n).or.                    &
     &                 (.not.hasv.and.i.eq.fi_i.and.j.eq.fi_j)
          END DO
        END DO
      END IF
!
!  One-shot system dump for offline matrix analysis (first SSA solve
!  only):  row/col/val triplets exactly as inserted, plus the RHS.
!
      IF (mrank.eq.0) THEN
        OPEN (371, file='ksp_dump_A.dat', status='replace',           &
     &        action='write')
        OPEN (372, file='ksp_dump_b.dat', status='replace',           &
     &        action='write')
      END IF
!
!-----------------------------------------------------------------------
!  Picard iteration.
!-----------------------------------------------------------------------
!
      DO nl_iter=1,nl_max_iter(ng)
!
!  Face viscosity at the current tile velocity (tile-local, then
!  gathered;  the assembly on rank 0 reads global faces so tile
!  halos never enter the operator).
!
        CALL ssa_fd_viscosity (ng, tile, model,                     &
     &                         LBi, UBi, LBj, UBj,                  &
     &                         Istr, Iend, Jstr, Jend,              &
     &                         pm, pn, hi, B_hard,                  &
     &                         urho, vrho, nuhx, nuhy)
!
        CALL fd_gather (nuhx, g_nuhx)
        CALL fd_gather (nuhy, g_nuhy)
        CALL fd_gather (urho, g_uold)
        CALL fd_gather (vrho, g_vold)
!
#  ifdef ISF_DEBUG
        IF (nl_iter.eq.1) THEN
          cff=0.0_r8
          DO jj=1,Mm_g
            DO ii=1,Lm_g
              cff=MAX(cff,ABS(g_beta(ii,jj)))
            END DO
          END DO
          WRITE(stdout,'(a,1p,e10.3,a,1p,e10.3,a,1p,e10.3)')        &
     &      ' FD-KSP inputs: max|g_beta|=',cff,                     &
     &      ' max|g_nuhx|=',MAXVAL(g_nuhx),                         &
     &      ' max|g_fdx|=',MAXVAL(ABS(g_fdx))
        END IF
#  endif
!
        IF (mrank.eq.0) THEN
!
!  Assemble the global operator:  eq1/eq2 coefficients verbatim
!  from ssa_fd_matvec, global bounds replacing the tile-array
!  bounds (cw=ce=cn=cs=0 and all switches off at the domain edge,
!  i.e. the zero-filled ghost ring plays the no-ice neighbor).
!
          CALL MatZeroEntries (pA, perr)
          CALL VecSet (pb, 0.0d0, perr)
          CALL VecSet (px, 0.0d0, perr)
!
          DO j=1,Mm_g
            DO i=1,Lm_g
              gi=2*((j-1)*Lm_g+i-1)
!
              IF (g_hi(i,j).le.0.0_r8) THEN
!
!  Identity rows (ice-free cell):  velocity zero.
!
                i1(1)=gi
                icol(1)=gi
                vcol(1)=one
                CALL MatSetValues (pA, 1, i1, 1, icol, vcol,        &
     &                             INSERT_VALUES, perr)
                i1(1)=gi+1
                icol(1)=gi+1
                CALL MatSetValues (pA, 1, i1, 1, icol, vcol,        &
     &                             INSERT_VALUES, perr)
                rscale(gi+1)=1.0_r8
                rscale(gi+2)=1.0_r8
                IF (nl_iter.eq.1) THEN
                  WRITE(371,'(i7,1x,i2)') gi, 1
                  WRITE(371,'(2i7,1x,1p,e23.16)') gi, gi, one
                  WRITE(371,'(i7,1x,i2)') gi+1, 1
                  WRITE(371,'(2i7,1x,1p,e23.16)') gi+1, gi+1, one
                END IF
              ELSE
!
                IF (i.gt.1) THEN
                  cwv=g_nuhx(i-1,j)
                ELSE
                  cwv=0.0_r8
                END IF
                IF (i.lt.Lm_g) THEN
                  cev=g_nuhx(i,j)
                ELSE
                  cev=0.0_r8
                END IF
                IF (j.lt.Mm_g) THEN
                  cnv=g_nuhy(i,j)
                ELSE
                  cnv=0.0_r8
                END IF
                IF (j.gt.1) THEN
                  csv=g_nuhy(i,j-1)
                ELSE
                  csv=0.0_r8
                END IF
!
                Nf=1; Ef=1; Sf=1; Wf=1
                NNE=1; NNW=1
                ENE=1; ESE=1
                SSW=1; SSE=1
                WSW=1; WNW=1
                IF (j.ge.Mm_g.or.g_hi(i,j+1).le.0.0_r8) THEN
                  Nf=0; NNE=0; NNW=0
                END IF
                IF (j.le.1.or.g_hi(i,j-1).le.0.0_r8) THEN
                  Sf=0; SSW=0; SSE=0
                END IF
                IF (i.ge.Lm_g.or.g_hi(i+1,j).le.0.0_r8) THEN
                  Ef=0; ENE=0; ESE=0
                END IF
                IF (i.le.1.or.g_hi(i-1,j).le.0.0_r8) THEN
                  Wf=0; WSW=0; WNW=0
                END IF
                IF (i.lt.Lm_g.and.j.lt.Mm_g) THEN
                  IF (g_hi(i+1,j+1).le.0.0_r8) THEN
                    NNE=0; ENE=0
                  END IF
                END IF
                IF (i.lt.Lm_g.and.j.gt.1) THEN
                  IF (g_hi(i+1,j-1).le.0.0_r8) THEN
                    ESE=0; SSE=0
                  END IF
                END IF
                IF (i.gt.1.and.j.gt.1) THEN
                  IF (g_hi(i-1,j-1).le.0.0_r8) THEN
                    SSW=0; WSW=0
                  END IF
                END IF
                IF (i.gt.1.and.j.lt.Mm_g) THEN
                  IF (g_hi(i-1,j+1).le.0.0_r8) THEN
                    WNW=0; NNW=0
                  END IF
                END IF
!
                dx=1.0_r8/g_pm(i,j)
                dy=1.0_r8/g_pn(i,j)
                dx2=dx*dx
                dy2=dy*dy
                d4=4.0_r8*dx*dy
                d2=2.0_r8*dx*dy
!
                eq1(1)=0.0_r8
                eq1(2)=-cnv*Nf/dy2
                eq1(3)=0.0_r8
                eq1(4)=-4.0_r8*cwv*Wf/dx2
                eq1(5)=(cnv*Nf+csv*Sf)/dy2+                         &
     &                 4.0_r8*(cev*Ef+cwv*Wf)/dx2
                eq1(6)=-4.0_r8*cev*Ef/dx2
                eq1(7)=0.0_r8
                eq1(8)=-csv*Sf/dy2
                eq1(9)=0.0_r8
                eq1(10)= cwv*Wf*WNW/d2+cnv*Nf*NNW/d4
                eq1(11)=(cnv*Nf*NNE-cnv*Nf*NNW)*Nf/d4+              &
     &                  (cwv*Wf-cev*Ef)*Nf/d2
                eq1(12)=-cev*Ef*ENE/d2-cnv*Nf*NNE/d4
                eq1(13)=(cwv*Wf*WSW-cwv*Wf*WNW)*Wf/d2+              &
     &                  (cnv*Nf-csv*Sf)*Wf/d4
                eq1(14)=(cnv*Nf*Ef*Nf-cnv*Nf*Wf*Nf-                 &
     &                   csv*Sf*Ef*Sf+csv*Sf*Wf*Sf)/d4+             &
     &                  (cev*Ef*Nf-cwv*Wf*Nf-cev*Ef*Sf+             &
     &                   cwv*Wf*Sf)/d2
                eq1(15)=(cev*Ef*ENE-cev*Ef*ESE)*Ef/d2+              &
     &                  (csv*Sf*Ef*Sf-cnv*Nf*Ef*Nf)/d4
                eq1(16)=-cwv*Wf*WSW/d2-csv*Sf*SSW/d4
                eq1(17)=(csv*Sf*SSW-csv*Sf*SSE)*Sf/d4+              &
     &                  (cev*Ef-cwv*Wf)*Sf/d2
                eq1(18)= cev*Ef*ESE/d2+csv*Sf*SSE/d4
!
                eq2(1)= cwv*Wf*WNW/d4+cnv*Nf*NNW/d2
                eq2(2)=(cnv*Nf*NNE-cnv*Nf*NNW)*Nf/d2+               &
     &                 (cwv*Wf-cev*Ef)*Nf/d4
                eq2(3)=-cev*Ef*ENE/d4-cnv*Nf*NNE/d2
                eq2(4)=(cwv*Wf*WSW-cwv*Wf*WNW)*Wf/d4+               &
     &                 (cnv*Nf-csv*Sf)*Wf/d2
                eq2(5)=(Ef-Wf)*(cnv*Nf-csv*Sf)/d2+                  &
     &                 (cev*Ef-cwv*Wf)*(Nf-Sf)/d4
                eq2(6)=(cev*Ef*ENE-cev*Ef*ESE)*Ef/d4+               &
     &                 (csv*Sf*Ef*Sf-cnv*Nf*Ef*Nf)/d2
                eq2(7)=-cwv*Wf*WSW/d4-csv*Sf*SSW/d2
                eq2(8)=(csv*Sf*SSW-csv*Sf*SSE)*Sf/d2+               &
     &                 (cev*Ef-cwv*Wf)*Sf/d4
                eq2(9)= cev*Ef*ESE/d4+csv*Sf*SSE/d2
                eq2(10)=0.0_r8
                eq2(11)=-4.0_r8*cnv*Nf/dy2
                eq2(12)=0.0_r8
                eq2(13)=-cwv*Wf/dx2
                eq2(14)=4.0_r8*(cnv*Nf+csv*Sf)/dy2+                 &
     &                 (cev*Ef+cwv*Wf)/dx2
                eq2(15)=-cev*Ef/dx2
                eq2(16)=0.0_r8
                eq2(17)=-4.0_r8*csv*Sf/dy2
                eq2(18)=0.0_r8
!
                bu=g_beta(i,j)
                bv=g_beta(i,j)
                eq1(5)=eq1(5)+bu
                eq2(14)=eq2(14)+bv
!
!  u-row:  9 u-neighbour columns then 9 v-neighbour columns.
!  LBC-pinned rows (and the translation reference cell) become
!  identity rows with zero right-hand side.
!
                IF (pin_u(i,j).or.ABS(eq1(5)).lt.eps) THEN
                  i1(1)=gi
                  icol(1)=gi
                  vcol(1)=one
                  CALL MatSetValues (pA, 1, i1, 1, icol, vcol,      &
     &                               INSERT_VALUES, perr)
                  rscale(gi+1)=1.0_r8
                  IF (nl_iter.eq.1) THEN
                    WRITE(371,'(i7,1x,i2)') gi, 1
                    WRITE(371,'(2i7,1x,1p,e23.16)') gi, gi, one
                  END IF
                ELSE
                  rm_row=0.0_r8
                  DO kk=1,18
                    rm_row=MAX(rm_row,ABS(eq1(kk)))
                  END DO
                  rscale(gi+1)=1.0_r8/MAX(rm_row,1.0d-300)
!
!  Free-edge rows reference neighbours outside the global index
!  range on the far side of the edge;  those columns are skipped
!  (the ghost velocity is a Dirichlet zero in PISM SSAFD, so the
!  column contributes nothing).  Without the skip, MatSetValues
!  rejects the whole row and the assembled matrix loses the edge
!  row entirely (the silent LU failure of fd_ksp=1 on MISMIP).
!
                  kk=0
                  k2=0
                  DO jc=1,-1,-1
                    DO ic=-1,1
                      gn=2*((j+jc-1)*Lm_g+i+ic-1)
                      kk=kk+1
                      IF (gn.ge.0.and.gn.lt.ndof) THEN
                        k2=k2+1
                        icol(k2)=gn
                        vcol(k2)=eq1(kk)*rscale(gi+1)
                      END IF
                    END DO
                  END DO
                  DO jc=1,-1,-1
                    DO ic=-1,1
                      gn=2*((j+jc-1)*Lm_g+i+ic-1)+1
                      kk=kk+1
                      IF (gn.ge.0.and.gn.lt.ndof) THEN
                        k2=k2+1
                        icol(k2)=gn
                        vcol(k2)=eq1(kk)*rscale(gi+1)
                      END IF
                    END DO
                  END DO
                  i1(1)=gi
                  CALL MatSetValues (pA, 1, i1, k2, icol, vcol,     &
     &                               INSERT_VALUES, perr)
                  IF (nl_iter.eq.1) THEN
                    WRITE(371,'(i7,1x,i2)') gi, k2
                    DO k=1,k2
                      WRITE(371,'(2i7,1x,1p,e23.16)') gi, icol(k),  &
     &                      vcol(k)
                    END DO
                  END IF
                END IF
!
!  v-row:  9 v-neighbour columns then 9 u-neighbour columns.
!
                IF (pin_v(i,j).or.ABS(eq2(14)).lt.eps) THEN
                  i1(1)=gi+1
                  icol(1)=gi+1
                  vcol(1)=one
                  CALL MatSetValues (pA, 1, i1, 1, icol, vcol,      &
     &                               INSERT_VALUES, perr)
                  rscale(gi+2)=1.0_r8
                  IF (nl_iter.eq.1) THEN
                    WRITE(371,'(i7,1x,i2)') gi+1, 1
                    WRITE(371,'(2i7,1x,1p,e23.16)') gi+1, gi+1, one
                  END IF
                ELSE
                  rm_row=0.0_r8
                  DO kk=10,18
                    rm_row=MAX(rm_row,ABS(eq2(kk)))
                  END DO
                  rscale(gi+2)=1.0_r8/MAX(rm_row,1.0d-300)
                  kk=0
                  k2=0
                  DO jc=1,-1,-1
                    DO ic=-1,1
                      gn=2*((j+jc-1)*Lm_g+i+ic-1)+1
                      kk=kk+1
                      IF (gn.ge.0.and.gn.lt.ndof) THEN
                        k2=k2+1
                        icol(k2)=gn
                        vcol(k2)=eq2(kk+9)*rscale(gi+2)
                      END IF
                    END DO
                  END DO
                  DO jc=1,-1,-1
                    DO ic=-1,1
                      gn=2*((j+jc-1)*Lm_g+i+ic-1)
                      kk=kk+1
                      IF (gn.ge.0.and.gn.lt.ndof) THEN
                        k2=k2+1
                        icol(k2)=gn
                        vcol(k2)=eq2(kk-9)*rscale(gi+2)
                      END IF
                    END DO
                  END DO
                  i1(1)=gi+1
                  CALL MatSetValues (pA, 1, i1, k2, icol, vcol,     &
     &                               INSERT_VALUES, perr)
                  IF (nl_iter.eq.1) THEN
                    WRITE(371,'(i7,1x,i2)') gi+1, k2
                    DO k=1,k2
                      WRITE(371,'(2i7,1x,1p,e23.16)') gi+1, icol(k),&
     &                      vcol(k)
                    END DO
                  END IF
                END IF
!
!  Right-hand side and initial guess (LBC-pinned rows carry zero
!  right-hand side and zero initial guess;  zero-diagonal fallback
!  rows keep the local-RHS identity semantics).
!
                IF (.not.pin_u(i,j)) THEN
                  CALL VecSetValue (pb, gi, g_fdx(i,j)*rscale(gi+1),&
     &                              INSERT_VALUES, perr)
                  CALL VecSetValue (px, gi, g_uold(i,j),            &
     &                              INSERT_VALUES, perr)
                  IF (nl_iter.eq.1) WRITE(372,'(i7,1x,1p,e23.16)')  &
     &              gi, g_fdx(i,j)*rscale(gi+1)
                END IF
                IF (.not.pin_v(i,j)) THEN
                  CALL VecSetValue (pb, gi+1, g_fdy(i,j)            &
     &                              *rscale(gi+2),                  &
     &                              INSERT_VALUES, perr)
                  CALL VecSetValue (px, gi+1, g_vold(i,j),          &
     &                              INSERT_VALUES, perr)
                  IF (nl_iter.eq.1) WRITE(372,'(i7,1x,1p,e23.16)')  &
     &              gi+1, g_fdy(i,j)*rscale(gi+2)
                END IF
              END IF
            END DO
          END DO
!
          CALL MatAssemblyBegin (pA, MAT_FINAL_ASSEMBLY, perr)
          CALL MatAssemblyEnd (pA, MAT_FINAL_ASSEMBLY, perr)
          CALL VecAssemblyBegin (pb, perr)
          CALL VecAssemblyEnd (pb, perr)
          CALL VecAssemblyBegin (px, perr)
          CALL VecAssemblyEnd (px, perr)
          IF (nl_iter.eq.1) THEN
            CLOSE (371)
            CLOSE (372)
          END IF
!
!  Zero-diagonal rows were already replaced by identity rows inside
!  the assembly loop (same ABS(diag) < eps test as ssa_fd_matvec,
!  beta included), so the RHS-driven identity semantics carry over.
!
          CALL KSPSetOperators (pksp, pA, pA, perr)
          CALL KSPSolve (pksp, pb, px, perr)
          CALL KSPGetConvergedReason (pksp, preason, perr)
          IF (preason.lt.0) THEN
            WRITE(stdout,'(a,i5)')                                  &
     &       ' SSA-FD-KSP WARNING: KSP diverged, reason=',preason
          END IF
!
#  ifdef ISF_DEBUG
          CALL VecNorm (px, NORM_INFINITY, resn, perr)
          WRITE(stdout,'(a,i4,a,i8,a,1p,e10.3)')                    &
     &      ' FD-KSP solve reason=',preason,' perr=',perr,          &
     &      ' max|x|=',resn
#  endif
!
!
#  ifdef ISF_DEBUG
!  Linear-solve exactness probe:  ||b - A x_LU|| / ||b|| must be at
!  roundoff level if the assembly and the direct solve are correct.
!
          CALL VecNorm (pb, NORM_2, rhsn, perr)
          CALL MatMult (pA, px, pw, perr)
          CALL VecAXPY (pw, -1.0d0, pb, perr)
          CALL VecNorm (pw, NORM_2, resn, perr)
          IF (rhsn.gt.eps) THEN
            WRITE(stdout,'(a,i3,a,1p,e10.3,a,e10.3)')               &
     &        ' FD-KSP lin probe iter',nl_iter,' ||b||=',rhsn,      &
     &        ' rel_lin_resid=',resn/rhsn
          ELSE
            WRITE(stdout,'(a,i3,a,1p,e10.3)')                       &
     &        ' FD-KSP lin probe iter',nl_iter,' ||b||=',rhsn
          END IF
#  endif
!
!  Extract, apply the speed-cap damper and the Picard relaxation.
!
          CALL VecGetArrayF90 (px, xv, perr)
          DO j=1,Mm_g
            DO i=1,Lm_g
              gi=2*((j-1)*Lm_g+i-1)
              g_unew(i,j)=xv(gi+1)
              g_vnew(i,j)=xv(gi+2)
            END DO
          END DO
          CALL VecRestoreArrayF90 (px, xv, perr)
!
          DO j=1,Mm_g
            DO i=1,Lm_g
              IF (g_hi(i,j).gt.0.0_r8) THEN
                u_slide=SQRT(g_unew(i,j)**2+g_vnew(i,j)**2)
                u_allow=1.5_r8*ssa_max_speed(ng)
                IF (u_slide.gt.u_allow) THEN
                  u_prev=SQRT(g_uold(i,j)**2+g_vold(i,j)**2)
                  IF (u_prev.gt.u_allow) THEN
                    u_allow=MAX(u_allow,0.9_r8*u_prev)
                  END IF
                  IF (u_slide.gt.u_allow) THEN
                    cff=u_allow/u_slide
                    g_unew(i,j)=cff*g_unew(i,j)
                    g_vnew(i,j)=cff*g_vnew(i,j)
                  END IF
                END IF
                IF (.not.pin_u(i,j)) THEN
                  g_unew(i,j)=picard_alpha(ng)*g_unew(i,j)+         &
     &                       (1.0_r8-picard_alpha(ng))*g_uold(i,j)
                ELSE
                  g_unew(i,j)=0.0_r8
                END IF
                IF (.not.pin_v(i,j)) THEN
                  g_vnew(i,j)=picard_alpha(ng)*g_vnew(i,j)+         &
     &                       (1.0_r8-picard_alpha(ng))*g_vold(i,j)
                ELSE
                  g_vnew(i,j)=0.0_r8
                END IF
              ELSE
                g_unew(i,j)=0.0_r8
                g_vnew(i,j)=0.0_r8
              END IF
            END DO
          END DO
!
!  Dual convergence criterion:  relative velocity change and the
!  true relative residual of the frozen-viscosity system at the
!  relaxed state.
!
          uo_n=0.0_r8
          un_n=0.0_r8
          vel_change=0.0_r8
          DO j=1,Mm_g
            DO i=1,Lm_g
              IF (g_hi(i,j).gt.0.0_r8) THEN
                IF (.not.pin_u(i,j)) THEN
                  uo_n=uo_n+g_uold(i,j)**2
                  un_n=un_n+g_unew(i,j)**2
                  vel_change=vel_change+                            &
     &              (g_unew(i,j)-g_uold(i,j))**2
                END IF
                IF (.not.pin_v(i,j)) THEN
                  uo_n=uo_n+g_vold(i,j)**2
                  un_n=un_n+g_vnew(i,j)**2
                  vel_change=vel_change+                            &
     &              (g_vnew(i,j)-g_vold(i,j))**2
                END IF
              END IF
            END DO
          END DO
          vel_norm=SQRT(un_n)
          vel_change=SQRT(vel_change)
          IF (vel_norm.gt.eps) THEN
            nl_resid=vel_change/vel_norm
          ELSE
            nl_resid=0.0_r8
          END IF
!
          DO j=1,Mm_g
            DO i=1,Lm_g
              gi=2*((j-1)*Lm_g+i-1)
              CALL VecSetValue (px, gi, g_unew(i,j),                &
     &                          INSERT_VALUES, perr)
              CALL VecSetValue (px, gi+1, g_vnew(i,j),              &
     &                          INSERT_VALUES, perr)
            END DO
          END DO
          CALL VecAssemblyBegin (px, perr)
          CALL VecAssemblyEnd (px, perr)
          CALL VecNorm (pb, NORM_2, rhsn, perr)
          CALL MatMult (pA, px, pw, perr)
          CALL VecAXPY (pw, -1.0d0, pb, perr)
          CALL VecNorm (pw, NORM_2, resn, perr)
          IF (rhsn.gt.eps) THEN
            nl_resid=MAX(nl_resid,resn/rhsn)
          END IF
!
          IF (nl_resid.lt.nl_tol(ng)) converged=.TRUE.
!
          WRITE(stdout,'(a,i3,a,1p,e9.2,a,e9.2,a,l1)')              &
     &      ' SSA-FD-KSP Picard iter',nl_iter,' nl_resid=',         &
     &      nl_resid,' max|u|=',vel_norm,' conv=',converged
!
!  The relaxed state becomes the tile state for the next
!  iteration (scatter copies it into every tile interior).
!
          DO j=1,Mm_g
            DO i=1,Lm_g
              g_uold(i,j)=g_unew(i,j)
              g_vold(i,j)=g_vnew(i,j)
            END DO
          END DO
!
        END IF  ! mrank.eq.0
!
        CALL fd_scatter (g_uold, urho)
        CALL fd_scatter (g_vold, vrho)
!
#  ifdef DISTRIBUTE
!  Refresh the velocity halos:  the westernmost/southernmost face
!  viscosities of each tile read the neighbour velocity from the
!  halo in the next Picard iteration (the eastern/northern seam
!  faces are copied one inward by ssa_fd_viscosity and do not).
!  Without this exchange the seam-face nuH would lag one Picard
!  iteration behind.  Interior values are untouched, so single-
!  tile runs are bit-identical with or without it.
!
        CALL mp_exchange2d (ng, tile, model, 2,                     &
     &                      LBi, UBi, LBj, UBj,                     &
     &                      NghostPoints,                           &
     &                      EWperiodic(ng), NSperiodic(ng),         &
     &                      urho, vrho)
#  endif
!
        iconv=0
        IF (converged) iconv=1
#  ifdef DISTRIBUTE
        CALL mpi_bcast (iconv, 1, MPI_INTEGER, 0,                   &
     &                  OCN_COMM_WORLD, mpi_err)
#  endif
        IF (iconv.eq.1) converged=.TRUE.
!
        IF (converged) EXIT
!
      END DO  ! nl_iter
!
      IF (.not.converged.and.Master) THEN
        WRITE(stdout,'(a,i4,a,es9.2)')                               &
     &      ' SSA-FD-KSP WARNING: Picard not converged in ',        &
     &      nl_max_iter(ng),' iterations; last nl_resid=',nl_resid
      END IF
!
      IF (mrank.eq.0) THEN
        CALL KSPDestroy (pksp, perr)
        CALL VecDestroy (px, perr)
        CALL VecDestroy (pb, perr)
        CALL VecDestroy (pw, perr)
        CALL MatDestroy (pA, perr)
      END IF
      CALL PetscFinalize (perr)
!
      deallocate ( g_hi, g_beta, g_fdx, g_fdy, g_pm, g_pn )
      deallocate ( g_nuhx, g_nuhy )
      deallocate ( g_uold, g_vold, g_unew, g_vnew )
      deallocate ( sbuf, rbuf )
      deallocate ( pin_u, pin_v )
      deallocate ( rscale )
#  ifdef DISTRIBUTE
      deallocate ( rects, counts, displs )
#  endif
!
      RETURN
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE fd_gather (fld, gfld)
!***********************************************************************
!
!  Gather one RHO-space tile-interior field into the rank-0 global
!  buffer (ghost ring stays zero).  MPI_GATHERV in distributed
!  runs;  direct copy in serial single-tile runs.
!
      real(r8), intent(in) :: fld(LBi:,LBj:)
      real(r8), intent(out) :: gfld(0:,0:)
!
      integer :: ii, jj, kk, p
!
#  ifdef DISTRIBUTE
      p=0
      DO jj=Jstr,Jend
        DO ii=Istr,Iend
          p=p+1
          sbuf(p)=fld(ii,jj)
        END DO
      END DO
      CALL mpi_gatherv (sbuf, nslab, MPI_DOUBLE_PRECISION,          &
     &                  rbuf, counts, displs,                        &
     &                  MPI_DOUBLE_PRECISION, 0,                     &
     &                  OCN_COMM_WORLD, mpi_err)
      IF (mrank.eq.0) THEN
        DO kk=1,nrpc
          DO jj=rects(4*(kk-1)+3),rects(4*(kk-1)+4)
            DO ii=rects(4*(kk-1)+1),rects(4*(kk-1)+2)
              gfld(ii,jj)=rbuf(displs(kk)+(jj-rects(4*(kk-1)+3))*   &
     &                   (rects(4*(kk-1)+2)-rects(4*(kk-1)+1)+1)+   &
     &                   ii-rects(4*(kk-1)+1)+1)
            END DO
          END DO
        END DO
      END IF
#  else
      DO jj=Jstr,Jend
        DO ii=Istr,Iend
          gfld(ii,jj)=fld(ii,jj)
        END DO
      END DO
#  endif
!
      RETURN
      END SUBROUTINE fd_gather
!
!***********************************************************************
      SUBROUTINE fd_scatter (gfld, fld)
!***********************************************************************
!
!  Scatter one rank-0 global field back into the tile interiors.
!
      real(r8), intent(in) :: gfld(0:,0:)
      real(r8), intent(out) :: fld(LBi:,LBj:)
!
      integer :: ii, jj, kk, p
!
#  ifdef DISTRIBUTE
      IF (mrank.eq.0) THEN
        DO kk=1,nrpc
          p=displs(kk)
          DO jj=rects(4*(kk-1)+3),rects(4*(kk-1)+4)
            DO ii=rects(4*(kk-1)+1),rects(4*(kk-1)+2)
              p=p+1
              rbuf(p)=gfld(ii,jj)
            END DO
          END DO
        END DO
      END IF
      CALL mpi_scatterv (rbuf, counts, displs,                      &
     &                   MPI_DOUBLE_PRECISION, sbuf, nslab,          &
     &                   MPI_DOUBLE_PRECISION, 0,                    &
     &                   OCN_COMM_WORLD, mpi_err)
      p=0
      DO jj=Jstr,Jend
        DO ii=Istr,Iend
          p=p+1
          fld(ii,jj)=sbuf(p)
        END DO
      END DO
#  else
      DO jj=Jstr,Jend
        DO ii=Istr,Iend
          fld(ii,jj)=gfld(ii,jj)
        END DO
      END DO
#  endif
!
      RETURN
      END SUBROUTINE fd_scatter
!
      END SUBROUTINE ssa_fd_picard_ksp
# endif
#endif
      END MODULE isf_solver_mod
