#include "cppdefs.h"

      MODULE isf_coupler_mod

#if defined ICESHELF && defined ICESHELF_THERMO
!
!git $Id$
!====================================================== Chuning Wang ===
!  Copyright (c) 2002-2025 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  Ice-ocean coupling routines for the ice shelf/sheet model.          !
!                                                                      !
!  This module provides:                                               !
!    1. Ice-ocean momentum coupling: computes the drag stress at the   !
!       ice-ocean interface and writes it to sustr/svstr so the ocean  !
!       momentum equation sees the correct top BC.  Three drag         !
!       parameterizations mirror the bottom drag options:              !
!         ICESHELF_LOGDRAG  logarithmic law (ZoSur_isf + von Karman)    !
!         ICESHELF_QDRAG    quadratic drag (tau = rdrag2_isf |du| du)   !
!         ICESHELF_LDRAG    linear drag (tau = rdrag_isf du)            !
!       With ICESHEET the drag acts on the ice-ocean relative velocity !
!       (ice at PSI points interpolated to U/V points); for a static   !
!       ice shelf (no ICESHEET) the ice velocity is zero and the drag  !
!       acts on the ocean top velocity alone.                          !
!    2. Friction velocity update: computes u* = |tau|^(1/2) at RHO     !
!       points for use by the 3EQN thermodynamics (isf_vbc).           !
!                                                                      !
!  References:                                                         !
!    Holland & Jenkins (1999), J. Phys. Oceanogr. 29, 1787-1800.       !
!    McPhee (1987), J. Geophys. Res. 92, 7017-7032.                    !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC :: isf_coupler
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE isf_coupler (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
      USE mod_isf
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
      CALL isf_coupler_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nrhs(ng),                                  &
     &                       GRID(ng) % rmask_isf,                      &
# if defined ICESHELF_LOGDRAG
     &                       GRID(ng) % ZoSur_isf,                      &
     &                       GRID(ng) % z_r,                            &
     &                       GRID(ng) % z_w,                            &
# endif
# if defined ICESHELF_LDRAG
     &                       GRID(ng) % rdrag_isf,                      &
# endif
# if defined ICESHELF_QDRAG
     &                       GRID(ng) % rdrag2_isf,                     &
# endif
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
     &                       ISF(ng) % Si,                              &
     &                       ISF(ng) % Fi,                              &
     &                       FORCES(ng) % sustr,                        &
     &                       FORCES(ng) % svstr)
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE isf_coupler
!
!***********************************************************************
      SUBROUTINE isf_coupler_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nrhs,                                &
     &                             rmask_isf,                           &
# if defined ICESHELF_LOGDRAG
     &                             ZoSur_isf, z_r, z_w,                 &
# endif
# if defined ICESHELF_LDRAG
     &                             rdrag_isf,                           &
# endif
# if defined ICESHELF_QDRAG
     &                             rdrag2_isf,                          &
# endif
     &                             u, v,                                &
     &                             Si, Fi,                              &
     &                             sustr, svstr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_isf
!
      USE bc_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: rmask_isf(LBi:,LBj:)
#  if defined ICESHELF_LOGDRAG
      real(r8), intent(in) :: ZoSur_isf(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,:)
#  endif
#  if defined ICESHELF_LDRAG
      real(r8), intent(in) :: rdrag_isf(LBi:,LBj:)
#  endif
#  if defined ICESHELF_QDRAG
      real(r8), intent(in) :: rdrag2_isf(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: Si(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: Fi(LBi:,LBj:,:)
      real(r8), intent(inout) :: sustr(LBi:,LBj:)
      real(r8), intent(inout) :: svstr(LBi:,LBj:)
# else
      real(r8), intent(in) :: rmask_isf(LBi:UBi,LBj:UBj)
#  if defined ICESHELF_LOGDRAG
      real(r8), intent(in) :: ZoSur_isf(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
#  if defined ICESHELF_LDRAG
      real(r8), intent(in) :: rdrag_isf(LBi:UBi,LBj:UBj)
#  endif
#  if defined ICESHELF_QDRAG
      real(r8), intent(in) :: rdrag2_isf(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: Si(LBi:UBi,LBj:UBj,2,nIsfS)
      real(r8), intent(inout) :: Fi(LBi:UBi,LBj:UBj,nIsfF)
      real(r8), intent(inout) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: svstr(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, know
!
      real(r8) :: du, dv, spd
      real(r8) :: ui_u, vi_v, ustar_rho
      real(r8) :: cff, cff1, cff2
# if defined ICESHELF_LOGDRAG
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wrk
# endif
!
# include "set_bounds.h"
!
!=======================================================================
!  Step 1: Ice-ocean momentum coupling.
!
!  Compute the drag stress at the ice-ocean interface:
!    LOGDRAG: tau = Cd_log |du| du,  Cd_log = (kappa/ln(dz/z0_sur))^2
!    QDRAG:   tau = rdrag2_isf |du| du
!    LDRAG:   tau = rdrag_isf du
!  where du = u_ice - u_ocean is the relative velocity (ice velocity is
!  zero for a static ice shelf without ICESHEET).  The drag parameters
!  rdrag_isf/rdrag2_isf/ZoSur_isf are grid arrays (defaulted to the
!  scalar rdrg_isf/rdrg2_isf/ZoSur_isf, overridden by the grid file
!  under ICESHELF_DRAG_GRID).
!
!  Ice velocity is at PSI points; interpolate to U/V points.
!  Ocean velocity is at U/V points (top layer k=N).
!  Write stress to sustr/svstr (ocean surface stress).
!  Compute friction velocity u* at RHO points for 3EQN.
!=======================================================================
!
      know=lsnew(ng)
!
# if defined ICESHELF_LOGDRAG
!  Precompute logarithmic drag coefficient at RHO points using the top
!  layer thickness dz = z_w(N) - z_r(N) and the ice shelf surface
!  roughness length, clamped to [Cdb_min, Cdb_max] as in set_vbc.F.
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          cff1=1.0_r8/LOG((z_w(i,j,N(ng))-z_r(i,j,N(ng)))/              &
     &                     ZoSur_isf(i,j))
          cff2=vonKar*vonKar*cff1*cff1
          wrk(i,j)=MIN(Cdb_max,MAX(Cdb_min,cff2))
        END DO
      END DO
# endif
!
!  U-stress at U-points.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
!
!  Check if ice shelf cavity exists at both adjacent RHO points.
!  Same criterion as set_vbc.F (rmask_isf*rmask_isf.eq.1) so the
!  drag stress is only applied where both cells are cavity water.
!
          IF (rmask_isf(i-1,j)*rmask_isf(i,j).eq.1.0_r8) THEN
!
!  Interpolate ice velocity from PSI to U-point.
!  PSI(i,j) is at the corner between RHO(i-1,j-1), RHO(i,j-1),
!  RHO(i-1,j), RHO(i,j).  U(i,j) is between RHO(i-1,j) and RHO(i,j).
!  Average the two adjacent PSI values in y:
!    u_ice(U) = 0.5*(Si(i,j) + Si(i,j+1))
!  For a static ice shelf (no ICESHEET) the ice velocity is zero.
!
# ifdef ICESHEET
            ui_u=0.5_r8*(Si(i,j  ,know,isIuis)+                         &
     &                   Si(i,j+1,know,isIuis))
# else
            ui_u=0.0_r8
# endif
!
!  Relative velocity: ice - ocean.
!  Ocean v is at V-points; interpolate to U-point (i,j) with the
!  same 4-point average as set_vbc.F:
!    v(U) = 0.25*(v(i-1,j) + v(i,j) + v(i-1,j+1) + v(i,j+1))
!
            du=ui_u-u(i,j,N(ng),nrhs)
# ifdef ICESHEET
            dv=0.5_r8*(Si(i,j  ,know,isIvis)+                           &
     &                 Si(i,j+1,know,isIvis))-                          &
# else
            dv=-                                                        &
# endif
     &         0.25_r8*(v(i-1,j  ,N(ng),nrhs)+                          &
     &                  v(i  ,j  ,N(ng),nrhs)+                          &
     &                  v(i-1,j+1,N(ng),nrhs)+                          &
     &                  v(i  ,j+1,N(ng),nrhs))
            spd=SQRT(du*du+dv*dv)
!
# if defined ICESHELF_LOGDRAG
            sustr(i,j)=0.5_r8*(wrk(i-1,j)+wrk(i,j))*du*spd
# elif defined ICESHELF_LDRAG
            sustr(i,j)=0.5_r8*(rdrag_isf(i-1,j)+rdrag_isf(i,j))*du
# else
            sustr(i,j)=0.5_r8*(rdrag2_isf(i-1,j)+rdrag2_isf(i,j))*      &
     &                 spd*du
# endif
          ELSE
!  Outside ice shelf cavity: zero stress so migrating grounding lines
!  do not leave stale drag values on newly-exposed ocean cells.
            sustr(i,j)=0.0_r8
          END IF
        END DO
      END DO
!
!  V-stress at V-points.
!
      DO j=JstrV,Jend
        DO i=Istr,Iend
          IF (rmask_isf(i,j-1)*rmask_isf(i,j).eq.1.0_r8) THEN
!
!  Interpolate ice velocity from PSI to V-point.
!  V(i,j) is between RHO(i,j-1) and RHO(i,j).
!  Average the two adjacent PSI values in x:
!    v_ice(V) = 0.5*(Si(i,j) + Si(i+1,j))
!
# ifdef ICESHEET
            vi_v=0.5_r8*(Si(i  ,j,know,isIvis)+                         &
     &                   Si(i+1,j,know,isIvis))
# else
            vi_v=0.0_r8
# endif
!
!  Relative velocity: ice - ocean.
!  Ocean u is at U-points; interpolate to V-point (i,j) with the
!  same 4-point average as set_vbc.F:
!    u(V) = 0.25*(u(i,j-1) + u(i+1,j-1) + u(i,j) + u(i+1,j))
!
            dv=vi_v-v(i,j,N(ng),nrhs)
# ifdef ICESHEET
            du=0.5_r8*(Si(i  ,j,know,isIuis)+                           &
     &                 Si(i+1,j,know,isIuis))-                          &
# else
            du=-                                                        &
# endif
     &         0.25_r8*(u(i  ,j-1,N(ng),nrhs)+                          &
     &                  u(i+1,j-1,N(ng),nrhs)+                          &
     &                  u(i  ,j  ,N(ng),nrhs)+                          &
     &                  u(i+1,j  ,N(ng),nrhs))
            spd=SQRT(du*du+dv*dv)
!
# if defined ICESHELF_LOGDRAG
            svstr(i,j)=0.5_r8*(wrk(i,j-1)+wrk(i,j))*dv*spd
# elif defined ICESHELF_LDRAG
            svstr(i,j)=0.5_r8*(rdrag_isf(i,j-1)+rdrag_isf(i,j))*dv
# else
            svstr(i,j)=0.5_r8*(rdrag2_isf(i,j-1)+rdrag2_isf(i,j))*      &
     &                 spd*dv
# endif
          ELSE
            svstr(i,j)=0.0_r8
          END IF
        END DO
      END DO
!
!  Exchange halo of the ice-shelf drag stress so that downstream
!  consumers (isf_vbc static u* path, output) read consistent U/V-point
!  values at tile boundaries.
!
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  sustr)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  svstr)
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    sustr, svstr)
# endif
!
!=======================================================================
!  Step 2: Compute friction velocity u* at RHO points (ICESHEET only).
!
!  u* = |tau|^(1/2) = (tau_x^2 + tau_y^2)^(1/4), valid for all three
!  drag schemes: QDRAG gives u* = sqrt(rdrag2_isf) |du|, LOGDRAG gives
!  u* = sqrt(Cd_log) |du|, and LDRAG gives u* = sqrt(rdrag_isf |du|).
!
!  Under ICESHEET this replaces the gamU-based u* in isf_vbc.  For a
!  static ice shelf (no ICESHEET) Step 2 is skipped: isf_vbc computes
!  u* from sustr/svstr with the gamU scaling factor instead.
!=======================================================================
!
# ifdef ICESHEET
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (rmask_isf(i,j).eq.1.0_r8) THEN
!
!  Average sustr from surrounding U-points and svstr from V-points.
!
            cff1=0.5_r8*(sustr(i,j)+sustr(i+1,j))
            cff2=0.5_r8*(svstr(i,j)+svstr(i,j+1))
            ustar_rho=SQRT(SQRT(cff1*cff1+cff2*cff2))
            Fi(i,j,icIust)=MAX(1.0E-6_r8, ustar_rho)
          ELSE
            Fi(i,j,icIust)=0.0_r8
          END IF
        END DO
      END DO
# endif
!
!  Basal melt rate passes through to thickness evolution unfiltered:
!  isf_thickness_adv reads Fi(i,j,icImel) (instantaneous melt rate from
!  isf_vbc) as the thickness source.  No EMA filtering is applied here.
!
      RETURN
      END SUBROUTINE isf_coupler_tile
!
#endif
      END MODULE isf_coupler_mod
