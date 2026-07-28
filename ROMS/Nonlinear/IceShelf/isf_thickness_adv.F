#include "cppdefs.h"

      MODULE isf_thickness_adv_mod

#ifdef ICESHEET_SSA
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2024 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!                                                                      !
!  This module computes ice thickness evolution via mass conservation: !
!                                                                      !
!    dH/dt = -div(u*H) + SMB - basal_melt                              !
!                                                                      !
!  using first-order upwind fluxes with CFL-based sub-stepping.        !
!  The ice velocity is taken from the SSA solver (PSI points) and      !
!  interpolated to RHO cell faces for flux computation.                !
!                                                                      !
!  Reference: STREAMICE streamice_advect_thickness.F                   !
!             PISM mass/transport.cc                                   !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
      PRIVATE
      PUBLIC  :: isf_thickness_adv
!
!  Global dot product shared variable (for MPI reduction).
!
      real(r8) :: isf_dot_prod
!$OMP THREADPRIVATE (isf_dot_prod)
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE isf_thickness_adv (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_ocean
      USE mod_scalars
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
      CALL wclock_on (ng, model, 43, __LINE__, MyFile)
# endif
      CALL isf_thickness_adv_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             GRID(ng) % pm,                       &
     &                             GRID(ng) % pn,                       &
     &                             ISF(ng) % Si,                        &
     &                             ISF(ng) % Fi)
# ifdef PROFILE
      CALL wclock_off (ng, model, 43, __LINE__, MyFile)
# endif
      RETURN
      END SUBROUTINE isf_thickness_adv
!
!***********************************************************************
      SUBROUTINE isf_thickness_adv_tile (ng, tile, model,               &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   pm, pn,                        &
     &                                   Si, Fi)
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
      integer :: i, j
      integer :: ksa_old, ksa_new
!
      real(r8) :: dt_full, dt_rem, dt_sub
      real(r8) :: max_icfl, loc_icfl, cfl_u, cfl_v
      real(r8) :: H_old, H_left, H_right, H_below, H_above
      real(r8) :: Fx_left, Fx_right, Fy_bottom, Fy_top
!
!  Work arrays.
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: H_new
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: u_xface
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: v_yface
!
# ifdef DISTRIBUTE
      character(len=3) :: op_handle
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", isf_thickness_adv_tile"
!
# include "set_bounds.h"
!
      SourceFile=MyFile
!
!  Set time level indices.
!
      ksa_old=lsold(ng)
      ksa_new=lsnew(ng)
!
!-----------------------------------------------------------------------
!  Compute ice thickness advection with CFL sub-stepping.
!-----------------------------------------------------------------------
!
# include <isf_thickness_adv.h>
!
      RETURN
      END SUBROUTINE isf_thickness_adv_tile
#endif
      END MODULE isf_thickness_adv_mod
