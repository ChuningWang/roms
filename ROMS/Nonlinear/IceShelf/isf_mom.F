#include "cppdefs.h"

      MODULE isf_mom_mod

#ifdef ICESHEET
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!========================================== Benjamin K. Galton-Fenzi ===
!                                                                      !
!  This module is the wrapper for the ICESHEET momentum solver.        !
!                                                                      !
!  It prepares the scalar ice state required by the SSA solver        !
!  contract (isf_solver_mod) and delegates the velocity solve:        !
!                                                                      !
!    1. exchange the lsold thickness/enthalpy halos and carry        !
!       them forward into the lsnew layer (components not            !
!       updated by the solve itself must stay valid);                !
!    2. cap the thickness at max_thickness (robustness);             !
!    3. compute the input diagnostics of the solver contract:        !
!         Fi(:,:,ifItis)  ice temperature T = E/H                    !
!         Fi(:,:,ifIzsf)  ice surface elevation (flotation branch)   !
!         Fi(:,:,ifIzbd)  ice bed elevation (flotation branch)       !
!    4. CALL isf_ssa (isf_solver_mod): Q1 FE Picard/PCG velocity    !
!       solve (LBC-aware, matrix-free, MPI seam-safe).               !
!                                                                      !
!  The masks, driving stress, viscosity, Glen hardness and the       !
!  entire Picard/PCG iteration live in isf_solver_mod; the           !
!  dynamic masks (rmask_ist, pmask_ist) come from isf_masks.F        !
!  through the GRID structure.  The legacy inline solver            !
!  (isf_driving_stress.h / isf_viscosity.h / isf_cg_solve.h) was    !
!  removed;  archived copies live in .agents/attic/.                !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
      USE isf_solver_mod, ONLY : isf_ssa
!
      implicit none
!
      PRIVATE
      PUBLIC  :: isf_mom
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE isf_mom (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
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
     &                   lsold(ng), lsnew(ng),                          &
     &                   GRID(ng) % h,                                  &
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
     &                         lsold, lsnew,                            &
     &                         h,                                       &
     &                         Si, Fi)
!***********************************************************************
!
      USE mod_kinds
      USE mod_param
      USE mod_scalars
      USE mod_parallel
      USE mod_isf
!
      USE exchange_2d_mod
      USE mp_exchange_mod, ONLY : mp_exchange2d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: lsold, lsnew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)    :: h(LBi:,LBj:)
      real(r8), intent(inout) :: Si(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: Fi(LBi:,LBj:,:)
# else
      real(r8), intent(in)    :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Si(LBi:UBi,LBj:UBj,2,nIsfS)
      real(r8), intent(inout) :: Fi(LBi:UBi,LBj:UBj,nIsfF)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
      integer :: ksa_old, ksa_new
!
      real(r8) :: H_rho
      real(r8) :: rho_w_local
!
!  RHO-point grounded-ice scratch mask (flotation criterion with the
!  static water depth), used only for the surface/bed elevation
!  branches below.  Tile-local, no halo exchange needed.
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rmask_ist_scr
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", isf_mom_tile"
!
# include "set_bounds.h"
!
!  Set time level indices.  The lsold/lsnew swap was performed in
!  iceshelf.F before calling isf_mom.  ksa_old = previous (input),
!  ksa_new = current (solution to be computed).
!
      ksa_old=lsold
      ksa_new=lsnew
!
!  Use reference seawater density for the flotation criterion.
!  Reuse mod_scalars rho0 (mean Boussinesq density, default 1025.0)
!  for consistency with the isf_masks.F grounding-line criterion.
!
      rho_w_local=rho0
!
!-----------------------------------------------------------------------
!  Exchange the lsold scalar-state halos and carry the state forward
!  into the lsnew layer.
!
!  The lsold/lsnew swap alternates the two time layers, and the ini
!  file populates layer 1 only (get_state Tindex=1).  Two issues:
!
!  1. The layer now appearing as ksa_new was last filled by the
!     previous call's carry-forward copy, which only touches the tile
!     interior; its ghost ring still holds the zero initial value.
!     The SSA solver reads H from ksa_old including halo cells (the
!     4 RHO corners of boundary PSI nodes), so the ksa_old halo must
!     be exchanged first; otherwise seam PSI nodes are deactivated
!     and the assembled operator is corrupted under MPI.
!
!  2. ksa_new must be carried forward from ksa_old so that any
!     component NOT updated by the solve (thickness advection is
!     skipped at cold start) does not leave the new layer at zero.
!     When thickness advection runs it overwrites ksa_new afterwards,
!     so this copy does not disturb the coupling.
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
!  Safety cap on ice thickness (prevent garbage values from
!  initialization or numerical artifacts from affecting the solve).
!  Configurable via max_thickness in isf.in (default 5000 m).
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          IF (Si(i,j,ksa_old,isIhis).gt.max_thickness(ng)) THEN
            Si(i,j,ksa_old,isIhis)=max_thickness(ng)
            Si(i,j,ksa_new,isIhis)=max_thickness(ng)
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  RHO-point grounded-ice mask from the flotation criterion (static
!  water depth, no zeta/Dcrit terms: the SSA-side criterion must not
!  chase high-frequency sea level oscillations).  Used below for the
!  surface/bed elevation branches;  the ocean-side masks (isf_masks.F)
!  are the authority for the cavity.
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
!  Solver-contract input diagnostics at RHO points:
!
!    Fi(:,:,ifItis)  ice temperature T = E/H
!    Fi(:,:,ifIzsf)  surface elevation:  grounded S = -h + H,
!                    floating S = (1 - rho_i/rho_w) H
!    Fi(:,:,ifIzbd)  bed elevation:      grounded B = -h,
!                    floating B = -(rho_i/rho_w) H
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          H_rho=Si(i,j,ksa_old,isIhis)
          IF (H_rho.gt.0.0_r8) THEN
            Fi(i,j,ifItis)=Si(i,j,ksa_old,isIeis)/H_rho
            IF (rmask_ist_scr(i,j).eq.1.0_r8) THEN
              Fi(i,j,ifIzsf)=-h(i,j)+H_rho
              Fi(i,j,ifIzbd)=-h(i,j)
            ELSE
              Fi(i,j,ifIzsf)=(1.0_r8-rho_i/rho_w_local)*H_rho
              Fi(i,j,ifIzbd)=-(rho_i/rho_w_local)*H_rho
            END IF
          ELSE
            Fi(i,j,ifItis)=0.0_r8
            Fi(i,j,ifIzsf)=0.0_r8
            Fi(i,j,ifIzbd)=-h(i,j)
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Solve the SSA velocity field (isf_solver_mod):  Q1 FE Picard/PCG,
!  LBC-aware, matrix-free.  All masks, driving stress, viscosity and
!  iteration control live there;  the solved velocity returns in
!  Si(:,:,lsnew,isIuis/isIvis) with halos exchanged.
!-----------------------------------------------------------------------
!
# ifdef ICESHEET_SSA
      CALL isf_ssa (ng, tile, model)
# endif
!
      RETURN
      END SUBROUTINE isf_mom_tile
#endif
      END MODULE isf_mom_mod
