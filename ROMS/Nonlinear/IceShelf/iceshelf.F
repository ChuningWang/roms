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
      USE mod_forces
!
# ifdef ICESHELF_THERMO
      USE isf_vbc_mod,      ONLY : isf_vbc
# endif
# if defined ICESHELF_THERMO && !defined ICESHELF_NO_SBL
      USE isf_sbl_mod,      ONLY : isf_sbl
# endif
# if defined ICESHEET && defined ICESHELF_THERMO
      USE isf_coupler_mod,  ONLY : isf_coupler
# endif
# ifdef ICESHELF_MOM
      USE isf_mom_mod,      ONLY : isf_mom
# endif
# ifdef ICESHEET
      USE isf_thickness_adv_mod, ONLY : isf_thickness_adv
# endif
# ifdef ICESHEET_CALVING
      USE isf_calving_mod,  ONLY : isf_calving
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
      real(r8), parameter :: dt_large = 1.0E+23_r8
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
!  MMS testing: skip ice-ocean thermodynamics.  isf_vbc (3EQN) melts
!  or freezes ice thickness between SSA calls, corrupting the exact
!  manufactured-solution state that the MMS operator check assumes.
!
#  ifndef SSA_MMS_TESTING
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL isf_vbc (ng, tile, model)
          END DO
!$OMP BARRIER
        END DO
#  endif
# endif

# if defined ICESHEET_SSA && defined ICESHELF_THERMO
!
!-----------------------------------------------------------------------
!  Ice-ocean coupling: momentum exchange (ice-ocean drag stress)
!  and basal melt rate temporal filtering (EMA).
!  Called after isf_vbc so that Fi(icImel) has the instantaneous melt.
!  Overwrites sustr/svstr under ice shelf, updates u*, filters melt.
!-----------------------------------------------------------------------
!
#  ifndef SSA_MMS_TESTING
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL isf_coupler (ng, tile, model)
          END DO
!$OMP BARRIER
        END DO
#  endif
# endif

# ifdef ICESHELF_MOM
!
!-----------------------------------------------------------------------
!  Compute the iceshelf momentum (SSA).  With ICESHEET_SSA, solved
!  every N_ssa ocean timesteps for multi-timescale coupling.
!-----------------------------------------------------------------------
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
#  ifdef ICESHEET_SSA
          IF (MOD(iic(ng)-1,N_ssa(ng)).eq.0) THEN
#  endif
            lsold(ng)=lsnew(ng)
            lsnew(ng)=3-lsold(ng)
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL isf_mom (ng, tile, model)
            END DO
!$OMP BARRIER
#  ifdef ICESHEET_SSA
          END IF
#  endif
        END DO
# endif

# ifdef ICESHEET_SSA
!
!-----------------------------------------------------------------------
!  Compute ice thickness advection (mass conservation).
!  Executed every N_ssa*N_thick ocean timesteps.
!-----------------------------------------------------------------------
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
!  Skip the cold-start first step: MOD(iic-1,N)=0 also triggers at
!  iic=ntstart=1 where no time has elapsed, but the update would apply
!  the full accumulated interval dt*N_ssa*N_thick, producing a spurious
!  thickness change (e.g. large 3EQN freeze-on at a cold start).
!  On a restart (ntstart>1) the first call IS needed: it covers the
!  interval since the last recorded update.
!  TODO: Re-enable thickness advection after zice update is implemented.
!  Currently disabled because isf_thickness_adv updates Si(isIhis) but
!  does not update zice, causing hwater=h+zice inconsistency in set_depth.
          IF (.FALSE.) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL isf_thickness_adv (ng, tile, model)
            END DO
!$OMP BARRIER
          END IF
        END DO
# endif

# ifdef ICESHEET_CALVING
!
!-----------------------------------------------------------------------
!  Apply calving (thickness threshold).
!  Executed together with thickness advection.
!-----------------------------------------------------------------------
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
!  Same cold-start skip as thickness advection above.
          IF (.FALSE.) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL isf_calving (ng, tile, model)
            END DO
!$OMP BARRIER
          END IF
        END DO
# endif

#if defined WET_DRY || defined ICESHEET_SSA
!
!-----------------------------------------------------------------------
!  Update grounding line position and ice shelf masks.
!  With ICESHEET_SSA: executed every N_ssa*N_gl ocean timesteps.
!  With WET_DRY only (no ICESHEET_SSA): executed every ocean timestep.
!
!  Under ICESHEET_SSA the masks must be recomputed even without WET_DRY
!  so that the dynamic zice/rmask_isf (flotation criterion) propagate to
!  set_depth.F and isf_vbc.F.  Without WET_DRY the wet/dry *_ist masks
!  are not updated (diagnostic only); the ocean grid mask (rmask) stays
!  static, so a migrating grounding line does not dry out ocean cells
!  (known limitation, see work_plan stage 4a.5).
!-----------------------------------------------------------------------
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
# ifdef ICESHEET_SSA
!  Cold-start skip as above; initial.F already ran isf_masks during
!  initialization, so the step-1 recomputation is redundant.
          IF (MOD(iic(ng)-1,N_ssa(ng)*N_gl(ng)).eq.0 .and.             &
     &        iic(ng).gt.ntstart(ng)) THEN
# endif
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL isf_masks (ng, tile, iNLM)
            END DO
!$OMP BARRIER
# ifdef ICESHEET_SSA
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
