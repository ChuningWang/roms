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
!  (diagnostic, Fi(:,:,ifItis)).  The melt and SMB enthalpy source     !
!  terms act at the column temperature T = E/H so that mass and heat   !
!  sources do not alter the ice temperature.                           !
!                                                                      !
!  Cadence: called from iceshelf.F every ocean timestep, decoupled     !
!  from the SSA solve (every N_ssa steps).  The velocity field is      !
!  the latest SSA solution: it stays in the new time level between     !
!  solves (the layer roles swap only on SSA steps, and the solver      !
!  overwrites the new level at each solve).  Integrating transport     !
!  at the ocean dt (instead of lumping it at the SSA cadence)          !
!  removes the O(N_ssa) staircase in the ice thickness seen by the     !
!  ocean-side coupling, so thermodynamics and mass-exchange            !
!  feedbacks vary smoothly in time.                                    !
!                                                                      !
!  Old-level write-back (L9 fix, 2026-09-27):  after the transported   !
!  state is written to the new level it is also copied back into the   !
!  old level (thickness and enthalpy only, not velocities).  Without   !
!  the write-back the old level accumulates only the source terms      !
!  between solves and the transport divergence enters the state only   !
!  at the window-end layer exchange, under-applying the transport by   !
!  a factor of N_ssa for N_ssa > 1 (diagnostic: transport increment    !
!  ratio N_ssa=10 vs N_ssa=1 = 0.1032).  With the write-back the old   !
!  level carries the running state between solves; the old-level halo  !
!  is refreshed by the source-term exchange at the start of the next   !
!  step, before the flux stencils read it.                             !
!                                                                      !
!  Numerical scheme (first-order upwind, flux form):                   !
!    - PSI velocity interpolated to RHO cell faces (2-point average)   !
!    - First-order upwind flux for monotonicity                        !
!    - Single explicit step at the ocean dt (ice-flow CFL << 1)        !
!    - Source terms added at the old time level before advection       !
!                                                                      !
!  Reference: STREAMICE streamice_advect_thickness.F                   !
!             PISM mass/transport.cc                                   !
!                                                                      !
!=======================================================================
!
      USE mod_grid
      USE mod_isf
      USE mod_param
      USE mod_scalars
      USE mod_stepping
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile, exchange_p2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
      implicit none
!
      PUBLIC  :: isf_advection
      PRIVATE
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE isf_advection (ng, tile, model)
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
      CALL wclock_on (ng, model, 43, __LINE__, MyFile)
# endif
      CALL isf_advection_tile (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         lsold(ng), lsnew(ng),                    &
     &                         GRID(ng) % pm,                           &
     &                         GRID(ng) % pn,                           &
     &                         GRID(ng) % rmask_ist,                    &
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
     &                               lsold, lsnew,                      &
     &                               pm, pn,                            &
     &                               rmask_ist,                         &
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
!    - Single explicit step at the ocean dt (ice-flow CFL << 1)
!    - Source terms (basal melt, SMB) added at the old time level
!
!  Cadence: called every ocean timestep (see module header).  The
!  velocity stays in the new level between solves (the latest SSA
!  solution; the layer roles swap only on SSA steps).  The transported
!  state written to the new level is copied back into the old level
!  below (L9 write-back) so the transport feeds back every step.
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
      real(r8), intent(in)    :: rmask_ist(LBi:,LBj:)
      real(r8), intent(inout) :: Si(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: Fi(LBi:,LBj:,:)
# else
      real(r8), intent(in)    :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: rmask_ist(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Si(LBi:UBi,LBj:UBj,2,nIsfS)
      real(r8), intent(inout) :: Fi(LBi:UBi,LBj:UBj,nIsfF)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
!
      real(r8) :: Fx, Fe, cff
!
!  Work arrays.
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ux
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ve
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", isf_advection_tile"
!
# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Advect ice thickness and enthalpy (single explicit upwind step).
!-----------------------------------------------------------------------
!
!  Add surface and bottom fluxes from surface mass balance and
!  ocean-driven melting.
!
!  H_new = H_old - dt*melt   + dt*SMB
!  E_new = E_old - dt*melt*T + dt*SMB*T
!
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (rmask_ist(i,j).gt.0.5_r8) THEN
!  Cap the melt by the available thickness (review finding P3):  an
!  uncapped melt on a thin advected front drives H negative, and the
!  negative "anti-ice" would keep transporting until divergence
!  removes it (isf_vbc already caps the tracer side the same way).
!
# if defined ICESHELF_THERMO || defined ICESHEET_PRESCRIBED_MELT
            cff=MIN(dt(ng)*Fi(i,j,ifImel),                                 &
     &              Si(i,j,lsold,isIhis))
            Si(i,j,lsold,isIhis)=Si(i,j,lsold,isIhis)                      &
     &                           -cff                                     &
     &                           +dt(ng)*smb_rate(ng)
            Si(i,j,lsold,isIeis)=Si(i,j,lsold,isIeis)                      &
     &                           -cff*Fi(i,j,ifItis)                       &
     &                           +dt(ng)*smb_rate(ng)*Fi(i,j,ifItis)
# else
            Si(i,j,lsold,isIhis)=Si(i,j,lsold,isIhis)                      &
     &                           +dt(ng)*smb_rate(ng)
            Si(i,j,lsold,isIeis)=Si(i,j,lsold,isIeis)                      &
     &                           +dt(ng)*smb_rate(ng)*Fi(i,j,ifItis)
# endif
          END IF
        END DO
      END DO
!
!  Exchange the modified old-level halo: the upwind fluxes below read
!  neighbor values at lsold, which the source terms above changed in
!  the tile interior.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile(ng,tile,                                 &
     &                         LBi,UBi,LBj,UBj,                         &
     &                         Si(:,:,lsold,isIhis))
        CALL exchange_r2d_tile(ng,tile,                                 &
     &                         LBi,UBi,LBj,UBj,                         &
     &                         Si(:,:,lsold,isIeis))
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d(ng,tile,model,2,                               &
     &                   LBi,UBi,LBj,UBj,                               &
     &                   NghostPoints,EWperiodic(ng),                   &
     &                   NSperiodic(ng),                                &
     &                   Si(:,:,lsold,isIhis),                          &
     &                   Si(:,:,lsold,isIeis))
# endif
!
!  Compute face velocities from PSI-point velocity.
!  ux(i,j): x-direction velocity on the x-face between
!           RHO(i-1,j) and RHO(i,j).
!  ve(i,j): y-direction velocity on the y-face between
!           RHO(i,j-1) and RHO(i,j).
!
!  Only the consumed ranges are computed (review finding P3):  the
!  flux loop reads ux(i..i+1, j) over i=Istr..Iend+1, j=Jstr..Jend and
!  ve(i, j..j+1) over i=Istr..Iend, j=Jstr..Jend+1;  the wider loop
!  above read PSI points one ring beyond the domain that no boundary
!  routine ever writes.
!
      DO j=Jstr,Jend
        DO i=Istr-1,Iend+1
          ux(i,j)=0.5_r8*(Si(i,j,lsnew,isIuis)+Si(i,j+1,lsnew,isIuis))
        END DO
      END DO
      DO j=Jstr-1,Jend+1
        DO i=Istr,Iend
          ve(i,j)=0.5_r8*(Si(i,j,lsnew,isIvis)+Si(i+1,j,lsnew,isIvis))
        END DO
      END DO
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
!
!  Xi and Eta direction upwind volume fluxes.
!  Left   Xi-face  (i  ): between RHO(i-1,j  ) and RHO(i  ,j  ).
!  Right  Xi-face  (i+1): between RHO(i  ,j  ) and RHO(i+1,j  ).
!  Bottom Eta-face (j  ): between RHO(i  ,j-1) and RHO(i  ,j  ).
!  Top    Eta-face (j+1): between RHO(i  ,j  ) and RHO(i  ,j+1).
!
          Fx=(MAX(ux(i  ,j),0.0_r8)*Si(i-1,j,lsold,isIhis)+             &
     &        MIN(ux(i  ,j),0.0_r8)*Si(i  ,j,lsold,isIhis)) -           &
     &       (MAX(ux(i+1,j),0.0_r8)*Si(i  ,j,lsold,isIhis)+             &
     &        MIN(ux(i+1,j),0.0_r8)*Si(i+1,j,lsold,isIhis))
          Fe=(MAX(ve(i,j  ),0.0_r8)*Si(i,j-1,lsold,isIhis)+             &
     &        MIN(ve(i,j  ),0.0_r8)*Si(i,j  ,lsold,isIhis)) -           &
     &       (MAX(ve(i,j+1),0.0_r8)*Si(i,j  ,lsold,isIhis)+             &
     &        MIN(ve(i,j+1),0.0_r8)*Si(i,j+1,lsold,isIhis))
          Si(i,j,lsnew,isIhis)=Si(i,j,lsold,isIhis) +                   &
     &                         dt(ng)*(Fx*pm(i,j)+Fe*pn(i,j))
!
!  Non-negativity floor on the transported thickness:  the explicit
!  upwind step is monotone but a diverging flux (e.g. under a capped
!  speed next to an over-thick pile) can take H below zero in one
!  step, and a negative state then corrupts every consumer (masks,
!  flotation, SSA, output -- the 2026-10-26 8 km crash).  Clamping at
!  zero keeps the damage local:  rmask_ist reclassifies the cell as
!  ice-free on the next isf_masks call.  The melt source above is
!  already clamped by MIN(dt*melt,H);  the floor is a no-op for
!  healthy states (all runs to date transport monotone positive).
!
          Si(i,j,lsnew,isIhis)=MAX(Si(i,j,lsnew,isIhis),0.0_r8)
!
!  Xi and Eta direction upwind enthalpy fluxes.
!  Left   Xi-face  (i  ): between RHO(i-1,j  ) and RHO(i  ,j  ).
!  Right  Xi-face  (i+1): between RHO(i  ,j  ) and RHO(i+1,j  ).
!  Bottom Eta-face (j  ): between RHO(i  ,j-1) and RHO(i  ,j  ).
!  Top    Eta-face (j+1): between RHO(i  ,j  ) and RHO(i  ,j+1).
!
          Fx=(MAX(ux(i  ,j),0.0_r8)*Si(i-1,j,lsold,isIeis)+             &
     &        MIN(ux(i  ,j),0.0_r8)*Si(i  ,j,lsold,isIeis)) -           &
     &       (MAX(ux(i+1,j),0.0_r8)*Si(i  ,j,lsold,isIeis)+             &
     &        MIN(ux(i+1,j),0.0_r8)*Si(i+1,j,lsold,isIeis))
          Fe=(MAX(ve(i,j  ),0.0_r8)*Si(i,j-1,lsold,isIeis)+             &
     &        MIN(ve(i,j  ),0.0_r8)*Si(i,j  ,lsold,isIeis)) -           &
     &       (MAX(ve(i,j+1),0.0_r8)*Si(i,j  ,lsold,isIeis)+             &
     &        MIN(ve(i,j+1),0.0_r8)*Si(i,j+1,lsold,isIeis))
          Si(i,j,lsnew,isIeis)=Si(i,j,lsold,isIeis) +                   &
     &                         dt(ng)*(Fx*pm(i,j)+Fe*pn(i,j))
!
!  Consistent enthalpy for floored (ice-free) cells:  E = h*T carries
!  no meaning at h = 0.
!
          IF (Si(i,j,lsnew,isIhis).le.0.0_r8) THEN
            Si(i,j,lsnew,isIeis)=0.0_r8
          END IF
!
!  Update ice sheet temperature T=E/H.
!
          IF (Si(i,j,lsnew,isIhis).gt.0.0_r8) THEN
            Fi(i,j,ifItis)=Si(i,j,lsnew,isIeis)/Si(i,j,lsnew,isIhis)
          ELSE
            Fi(i,j,ifItis)=0.0_r8
          END IF
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile(ng,tile,                                 &
     &                         LBi,UBi,LBj,UBj,                         &
     &                         Si(:,:,lsnew,isIhis))
        CALL exchange_r2d_tile(ng,tile,                                 &
     &                         LBi,UBi,LBj,UBj,                         &
     &                         Si(:,:,lsnew,isIeis))
        CALL exchange_r2d_tile(ng,tile,                                 &
     &                         LBi,UBi,LBj,UBj,                         &
     &                         Fi(:,:,ifItis))
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d(ng,tile,model,3,                               &
     &                   LBi,UBi,LBj,UBj,                               &
     &                   NghostPoints,EWperiodic(ng),                   &
     &                   NSperiodic(ng),                                &
     &                   Si(:,:,lsnew,isIhis),                          &
     &                   Si(:,:,lsnew,isIeis),                          &
     &                   Fi(:,:,ifItis))
# endif
!
!-----------------------------------------------------------------------
!  Write the transported state back into the old time level (L9 fix,
!  2026-09-27):  thickness and enthalpy only, NOT the velocities (the
!  solver overwrites the velocities in the new level at each solve and
!  they stay frozen there between solves).  Without this feedback the
!  old level accumulates only the source terms between solves and the
!  transport enters the state just once per window, at the layer
!  exchange, under-applying it by a factor of N_ssa (N_ssa > 1).  The
!  interior-only copy leaves the old-level halo stale for one step;
!  the exchange after the source terms at the top of the next call
!  refreshes it before the flux stencils read it.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Si(i,j,lsold,isIhis)=Si(i,j,lsnew,isIhis)
          Si(i,j,lsold,isIeis)=Si(i,j,lsnew,isIeis)
        END DO
      END DO
!
      RETURN
      END SUBROUTINE isf_advection_tile
#endif
      END MODULE isf_advection_mod
