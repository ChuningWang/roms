#include "cppdefs.h"
      MODULE isf_coupling_mod
#ifdef ICESHEET
!
!git $Id$
!====================================================== Chuning Wang ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  SSA coupling window accumulator for the ice sheet model.            !
!                                                                      !
!  This module is the ONLY place where the Ci coupling variables       !
!  are accumulated and processed.  It is called at the top of the      !
!  iceshelf.F driver, before every other ice computation, so the       !
!  window state is up to date before its consumers read it on the      !
!  same timestep.  One field is maintained:                            !
!                                                                      !
!    Ci(:,:,icIzta)  free surface zeta accumulated over the SSA        !
!                    coupling window [m]; at the last step of each     !
!                    window it holds the window MEAN, which the        !
!                    basal drag (isf_stress) reads at that step to     !
!                    build Fi(:,:,ifIbta) for the SSA solve.           !
!                                                                      !
!  Window phase (window_phase = MOD(iic,N_ssa), keyed to the           !
!  absolute timestep iic so partial windows continue seamlessly        !
!  across restart; the accumulator is part of the restart state,       !
!  netCDF variable is_zta):                                            !
!                                                                      !
!    N_ssa <= 1          degenerate single-step window: Ci = zeta      !
!                        (the instantaneous sample, explicit branch)   !
!    window_phase = 1    window first step, including the cold-        !
!                        start first step: Ci = zeta (clear + seed     !
!                        with the variable's initial value)            !
!    window_phase >= 2   Ci = Ci + zeta (accumulate)                   !
!    window_phase = 0    window last step (MOD(iic,N_ssa)=0, the       !
!                        SSA solve step): the accumulate above adds    !
!                        the final sample, then Ci = Ci/N_ssa          !
!                        finalizes the mean of the cycle just          !
!                        ended; the next window's first step           !
!                        clears it and starts a new cycle.             !
!                                                                      !
!  The per-cell update needs no halo exchange:                         !
!  each tile accumulates its ghost points from the exchanged zeta      !
!  identically to their owner tiles.                                   !
!                                                                      !
!  Future coupling accumulators (e.g. a cavity-mask floating           !
!  fraction for window-mean geometry) extend here: add the icI         !
!  index in isf_mod.h, allocate the Ci slot in mod_isf.F, and          !
!  add the field to the loops below.                                   !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
      PRIVATE
      PUBLIC :: isf_coupling
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE isf_coupling (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping, ONLY : knew
      USE mod_isf,      ONLY : ISF
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
      CALL isf_coupling_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        OCEAN(ng) % zeta(:,:,knew(ng)),           &
     &                        ISF(ng) % Ci)
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE isf_coupling
!
!***********************************************************************
      SUBROUTINE isf_coupling_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              zeta, Ci)
!***********************************************************************
!
!  Advance the SSA coupling window accumulators by one ocean step.
!  See the module header for the phase scheme; isf_stress reads the
!  finalized window mean on the SSA solve steps (MOD(iic,N_ssa)=0).
!
      USE mod_param
      USE mod_scalars
      USE mod_isf,      ONLY : N_ssa, nIsfC, icIzta
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: zeta(LBi:,LBj:)
      real(r8), intent(inout) :: Ci(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Ci(LBi:UBi,LBj:UBj,nIsfC)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, window_phase
!
      real(r8) :: cff1
!
# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Accumulate the free-surface sample of the current ocean step and
!  finalize or reseed the SSA coupling window.
!-----------------------------------------------------------------------
!
      IF (N_ssa(ng).le.1) THEN
!
!  Degenerate single-step window (N_ssa = 1): the accumulator holds
!  the instantaneous sample.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Ci(i,j,icIzta)=zeta(i,j)
          END DO
        END DO
      ELSE
        window_phase=MOD(iic(ng),N_ssa(ng))
        cff1=1.0_r8/REAL(N_ssa(ng),r8)
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            IF (window_phase.eq.1) THEN
              Ci(i,j,icIzta)=zeta(i,j)
            ELSE
              Ci(i,j,icIzta)=Ci(i,j,icIzta)+zeta(i,j)
            END IF
            IF (window_phase.eq.0) THEN
              Ci(i,j,icIzta)=Ci(i,j,icIzta)*cff1
            END IF
          END DO
        END DO
      END IF
!
      RETURN
      END SUBROUTINE isf_coupling_tile
!
#endif
      END MODULE isf_coupling_mod
