#include "cppdefs.h"

      MODULE iceshelf_mod

#ifdef ICESHELF_DYNAM
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2024 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           W. Paul Budgell    !
!    See License_ROMS.md                            Katherine Hedstrom !
!================================================== Hernan G. Arango ===
!                                                       Chuning Wang   !
!  This module is the primary driver routine for the ROMS nonlinear    !
!  iceshelf model.                                                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
      USE mod_stepping
      USE mod_isf
      USE isf_stress_mod,  ONLY : isf_stress
!
# ifdef ICESHELF_THERMO
      USE isf_vbc_mod,      ONLY : isf_vbc
# endif
# if defined ICESHELF_THERMO && !defined ICESHELF_NO_SBL
      USE isf_sbl_mod,      ONLY : isf_sbl
# endif
# ifdef ICESHEET
      USE isf_mom_mod,      ONLY : isf_mom
# endif
# ifdef ICESHEET
      USE isf_advection_mod, ONLY : isf_advection
# endif
# ifdef ICESHEET
      USE isf_bc2d_mod,     ONLY : isf_bc2d
      USE isf_pibc_mod,     ONLY : isf_pibc
      USE isf_tibc_mod,     ONLY : isf_tibc
# endif
# ifdef ICESHEET_CALVING
      USE isf_calving_mod,  ONLY : isf_calving_mask
# endif
      USE strings_mod,      ONLY : FoundError, uppercase
# if defined WET_DRY
      USE isf_masks_mod,    ONLY : isf_masks
# endif
!
      implicit none
!
      PUBLIC  :: iceshelf
      PRIVATE
!
      CONTAINS
!
      SUBROUTINE iceshelf (model)
!
!=======================================================================
!                                                                      !
!  Iceshelf model timestepping driver.                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: model
!
!  Local variable declarations.
!
      integer :: tile
      integer :: ig, ng, nl
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!=======================================================================
!  Timestep iceshelf model equations.
!=======================================================================
!
!  In nesting applications, the number of nesting layers (NestLayers) is
!  used to facilitate refinement grids and composite/refinament grids
!  combinations. Otherwise, the solution it is looped once for a single
!  grid application (NestLayers = 1).
!
      NEST_LAYER : DO nl=1,NestLayers
!
!-----------------------------------------------------------------------
!  Ice-ocean momentum coupling and surface flux masking.
!
!  Computes the drag stress at the ice-ocean interface (sustr/svstr),
!  zeroes atmospheric surface tracer fluxes over cavity cells, and
!  (under ICESHELF_THERMO) computes the friction velocity u* for the
!  3EQN thermodynamics.  Must run before isf_vbc, which reads sustr/
!  svstr and Fi(ifIust).
!-----------------------------------------------------------------------
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL isf_stress (ng, tile, iNLM)
          END DO
!$OMP BARRIER
        END DO

# ifdef ICESHELF_THERMO
!
!-----------------------------------------------------------------------
!  Compute the ice-ocean thermodynamics.
!-----------------------------------------------------------------------
!
#  if !defined ICESHELF_NO_SBL
!
!  Compute sub-boundary layer weights on the first time step only.
!  Oi depends only on static geometry (h, zice, Hz), so the result
!  is the same whether computed from cold start or restart.
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          IF (iic(ng).eq.ntstart(ng)) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL isf_sbl (ng, tile)
            END DO
!$OMP BARRIER
          END IF
        END DO
#  endif
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL isf_vbc (ng, tile, model)
          END DO
!$OMP BARRIER
        END DO
# endif

# ifdef ICESHEET
!
!-----------------------------------------------------------------------
!  Compute the iceshelf momentum (SSA).  With ICESHEET, solved
!  every N_ssa ocean timesteps for multi-timescale coupling.  The
!  scalar state (thickness, enthalpy) is advected below every ocean
!  timestep, decoupled from the SSA cadence.
!-----------------------------------------------------------------------
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
!
!  Advance the iceshelf time levels every ocean timestep.  Since the
!  advection below integrates the scalar state over each ocean dt,
!  the lsold/lsnew swap is no longer tied to N_ssa: it runs on SSA
!  steps (the solver expects a fresh target level) and on every
!  step that advects.  At the initialization timestep neither runs
!  and the swap is skipped so the restart/current layer is kept.
!
          IF (MOD(iic(ng)-1,N_ssa(ng)).eq.0) THEN
            lsold(ng)=lsnew(ng)
            lsnew(ng)=3-lsold(ng)
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL isf_mom (ng, tile, model)
            END DO
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL isf_pibc (ng, tile, model)
            END DO
!$OMP BARRIER
          END IF
!
          IF (MOD(iic(ng)-1,N_ssa(ng)).eq.0) THEN
!
!  Apply the velocity lateral boundary conditions at the SSA cadence:
!  the velocity only changes when the solver runs.  Order matters:
!  the velocity boundary (isf_pibc) and the thickness boundary
!  (isf_bc2d) must be applied before the enthalpy boundary
!  (isf_tibc), which computes the boundary enthalpy as h*T and
!  divides by the boundary thickness.
!
!$OMP BARRIER
          END IF
!
!  Advect ice thickness and enthalpy (mass and heat conservation)
!  every ocean timestep, with the SSA velocity held constant between
!  solves.  Integrating transport at the ocean dt removes the
!  O(N_ssa) staircase in the thickness seen by the ocean coupling:
!  the basal melt source uses the instantaneous rate Fi(:,:,ifImel)
!  from isf_vbc (no window accumulation), and the advected geometry
!  propagates to the ocean through isf_masks below, which also runs
!  every timestep and keeps zice consistent with the new thickness
!  for set_depth (the reason advection was previously disabled).
!
          IF (iic(ng).gt.ntstart(ng)) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL isf_advection (ng, tile, model)
            END DO
!$OMP BARRIER
          END IF
!
!  Apply the scalar lateral boundary conditions every timestep: the
!  advection above overwrites the new time level each step, so its
!  ghost ring and boundary values must be refreshed at the same
!  cadence (the velocity boundary above runs only on SSA steps).
!
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL isf_bc2d (ng, tile, model)
          END DO
!$OMP BARRIER
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL isf_tibc (ng, tile, model)
          END DO
!$OMP BARRIER
# ifdef ICESHEET_CALVING
!
!  Apply calving (pre-defined mask), after thickness advection: ice
!  advected into cells with calving_mask = 0 is removed (H, E, T
!  zeroed).  Runs every timestep: the mask is a static constraint
!  and enforcing it is idempotent, so boundary input into the
!  forbidden region is removed immediately.
!
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL isf_calving_mask (ng, tile, model)
          END DO
!$OMP BARRIER
# endif
        END DO
# endif

#if defined WET_DRY || defined ICESHEET
!
!-----------------------------------------------------------------------
!  Update grounding line position and ice shelf masks.
!  With ICESHEET: executed every ocean timestep, after thickness
!  advection, so zice and the cavity masks track the advected ice
!  geometry each step and the ocean (set_depth, isf_vbc, isf_stress)
!  sees a continuously updated ice state instead of a staircase at
!  the SSA cadence.  The previous grounding-line cadence parameter
!  N_gl is no longer referenced here (input retained, unused).
!  With WET_DRY only (no ICESHEET): executed every ocean timestep.
!
!  Under ICESHEET the masks must be recomputed even without WET_DRY
!  so that the dynamic zice/rmask_isf (flotation criterion) propagate to
!  set_depth.F and isf_vbc.F.  Without WET_DRY the wet/dry *_ist masks
!  are not updated (diagnostic only); the ocean grid mask (rmask) stays
!  static, so a migrating grounding line does not dry out ocean cells
!  (known limitation, see work_plan stage 4a.5).
!-----------------------------------------------------------------------
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
# ifdef ICESHEET
!  Cold-start skip as above; initial.F already ran isf_masks during
!  initialization, so the step-1 recomputation is redundant.
          IF (iic(ng).gt.ntstart(ng)) THEN
# endif
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL isf_masks (ng, tile, iNLM)
            END DO
!$OMP BARRIER
# ifdef ICESHEET
          END IF
# endif
        END DO
#endif

      END DO NEST_LAYER
!
      RETURN
      END SUBROUTINE iceshelf
!
#endif
      END MODULE iceshelf_mod
