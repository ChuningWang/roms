#include "cppdefs.h"

      MODULE isf_stress_mod

#ifdef ICESHELF
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
!         ICESHELF_LOGDRAG  logarithmic law (ZoSur_isf + von Karman)   !
!         ICESHELF_QDRAG    quadratic drag (tau = rdrag2_isf |du| du)  !
!         ICESHELF_LDRAG    linear drag (tau = rdrag_isf du)           !
!       With ICESHEET the drag acts on the ice-ocean relative velocity !
!       (ice at PSI points interpolated to U/V points); for a static   !
!       ice shelf (no ICESHEET) the ice velocity is zero and the drag  !
!       acts on the ocean top velocity alone.                          !
!    2. Surface tracer flux masking: zeroes stflx (all tracers) and    !
!       srflx (SHORTWAVE) over ice shelf cavity cells (rmask_isf=1)    !
!       so atmospheric fluxes are replaced by ice-ocean fluxes from    !
!       the 3EQN thermodynamics in isf_vbc.                            !
!    3. Friction velocity update (ICESHELF_THERMO only): computes u* = !
!       gamU*|tau|^(1/2) at RHO points for use by the 3EQN             !
!       thermodynamics (isf_vbc).                                      !
!    4. Ice sheet basal drag coefficient (ICESHEET only): computes the !
!       SSA bottom drag beta at RHO points from the window-mean free   !
!       surface (Ci(:,:,icIzta), owned and finalized by isf_coupling   !
!       at the solve step) and stores it in Fi(:,:,ifIbta) for the     !
!       SSA momentum solve (isf_ssa).  Two sliding laws, CPP-selected  !
!       with Weertman as the default:                                  !
!         default                   Weertman power law                 !
!                                       beta = C^2 |u|^(m-1)           !
!         ICESHEET_COULOMB  regularized Coulomb (Tsai et al.           !
!                                       2015) with effective pressure  !
!                                       from hydraulic connectivity.   !
!       The float state follows the isf_masks.F flotation criterion    !
!       evaluated with the window-averaged free surface.               !
!       With ICESHEET_SMOOTH_GL the drag of both sliding laws is       !
!       scaled by a cosine-ramped grounded fraction (STREAMICE         !
!       float_frac) about the flotation thickness; the ramp width is   !
!       gl_float_width [m ice thickness].                              !
!                                                                      !
!  References:                                                         !
!    Holland & Jenkins (1999), J. Phys. Oceanogr. 29, 1787-1800.       !
!    McPhee (1987), J. Geophys. Res. 92, 7017-7032.                    !
!    Weertman (1957), J. Glaciol. 3, 33-38.                            !
!    Tsai et al. (2015), J. Glaciol. 61, 649-655.                      !
!    MITgcm STREAMICE streamice_upd_ffrac_uncoupled.F (float_frac).    !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC :: isf_stress
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE isf_stress (ng, tile, model)
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
      CALL isf_stress_tile (ng, tile, model,                            &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nrhs(ng), lsnew(ng),                       &
     &                       GRID(ng) % rmask_isf,                      &
# if defined ICESHELF_LOGDRAG
     &                       GRID(ng) % ZoSur_isf,                      &
     &                       GRID(ng) % z_r,                            &
     &                       GRID(ng) % z_w,                            &
# elif defined ICESHELF_QDRAG
     &                       GRID(ng) % rdrag2_isf,                     &
# elif defined ICESHELF_LDRAG
     &                       GRID(ng) % rdrag_isf,                      &
# endif
# ifdef ICESHEET
     &                       GRID(ng) % h,                              &
# endif
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
     &                       ISF(ng) % Si,                              &
     &                       ISF(ng) % Fi,                              &
     &                       ISF(ng) % Ci,                              &
     &                       FORCES(ng) % stflx,                        &
# ifdef SHORTWAVE
     &                       FORCES(ng) % srflx,                        &
# endif
     &                       FORCES(ng) % sustr,                        &
     &                       FORCES(ng) % svstr)
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE isf_stress
!
!***********************************************************************
      SUBROUTINE isf_stress_tile (ng, tile, model,                      &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nrhs, lsnew,                         &
     &                             rmask_isf,                           &
# if defined ICESHELF_LOGDRAG
     &                             ZoSur_isf, z_r, z_w,                 &
# elif defined ICESHELF_QDRAG
     &                             rdrag2_isf,                          &
# elif defined ICESHELF_LDRAG
     &                             rdrag_isf,                           &
# endif
# ifdef ICESHEET
     &                             h,                                   &
# endif
     &                             u, v,                                &
     &                             Si, Fi, Ci,                          &
     &                             stflx,                               &
# ifdef SHORTWAVE
     &                             srflx,                               &
# endif
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
      integer, intent(in) :: nrhs, lsnew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: rmask_isf(LBi:,LBj:)
#  if defined ICESHELF_LOGDRAG
      real(r8), intent(in) :: ZoSur_isf(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,:)
#  elif defined ICESHELF_QDRAG
      real(r8), intent(in) :: rdrag2_isf(LBi:,LBj:)
#  elif defined ICESHELF_LDRAG
      real(r8), intent(in) :: rdrag_isf(LBi:,LBj:)
#  endif
#  ifdef ICESHEET
      real(r8), intent(in) :: h(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: Si(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: Fi(LBi:,LBj:,:)
      real(r8), intent(in) :: Ci(LBi:,LBj:,:)
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
#  ifdef SHORTWAVE
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: sustr(LBi:,LBj:)
      real(r8), intent(inout) :: svstr(LBi:,LBj:)
# else
      real(r8), intent(in) :: rmask_isf(LBi:UBi,LBj:UBj)
#  if defined ICESHELF_LOGDRAG
      real(r8), intent(in) :: ZoSur_isf(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  elif defined ICESHELF_QDRAG
      real(r8), intent(in) :: rdrag2_isf(LBi:UBi,LBj:UBj)
#  elif defined ICESHELF_LDRAG
      real(r8), intent(in) :: rdrag_isf(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ICESHEET
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: Si(LBi:UBi,LBj:UBj,2,nIsfS)
      real(r8), intent(inout) :: Fi(LBi:UBi,LBj:UBj,nIsfF)
      real(r8), intent(in) :: Ci(LBi:UBi,LBj:UBj,nIsfC)
      real(r8), intent(inout) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
#  ifdef SHORTWAVE
      real(r8), intent(inout) :: srflx(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: svstr(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, itrc
!
      real(r8) :: du, dv, spd
      real(r8) :: ui_u, vi_v, ustar_rho
      real(r8) :: cff, cff1, cff2
# if defined ICESHELF_LOGDRAG
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wrk
# endif
# ifdef ICESHEET
      real(r8) :: u_mid, v_mid, u_slide
      real(r8) :: C_sq, zeta_avg, D_avg
      real(r8) :: H_eff, H_float
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
!  scalar rdrg_isf/rdrg2_isf/ZoSr_isf, overridden by the grid file
!  under ICESHELF_DRAG_GRID).
!
!  Ice velocity is at PSI points; interpolate to U/V points.
!  Ocean velocity is at U/V points (top layer k=N).
!  Write stress to sustr/svstr (ocean surface stress).
!=======================================================================
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
!  Drag stress is only applied where both cells are cavity water.
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
            ui_u=0.5_r8*(Si(i,j  ,lsnew,isIuis)+                        &
     &                   Si(i,j+1,lsnew,isIuis))
# else
            ui_u=0.0_r8
# endif
!
!  Relative velocity: ice - ocean.
!  Ocean v is at V-points; interpolate to U-point (i,j) with the
!  4-point average:
!    v(U) = 0.25*(v(i-1,j) + v(i,j) + v(i-1,j+1) + v(i,j+1))
!
            du=ui_u-u(i,j,N(ng),nrhs)
# ifdef ICESHEET
            dv=0.5_r8*(Si(i,j  ,lsnew,isIvis)+                          &
     &                 Si(i,j+1,lsnew,isIvis))-                         &
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
# elif defined ICESHELF_QDRAG
            sustr(i,j)=0.5_r8*(rdrag2_isf(i-1,j)+rdrag2_isf(i,j))*      &
     &                 spd*du
# elif defined ICESHELF_LDRAG
            sustr(i,j)=0.5_r8*(rdrag_isf(i-1,j)+rdrag_isf(i,j))*du
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
            vi_v=0.5_r8*(Si(i  ,j,lsnew,isIvis)+                        &
     &                   Si(i+1,j,lsnew,isIvis))
# else
            vi_v=0.0_r8
# endif
!
!  Relative velocity: ice - ocean.
!  Ocean u is at U-points; interpolate to V-point (i,j) with the
!  4-point average:
!    u(V) = 0.25*(u(i,j-1) + u(i+1,j-1) + u(i,j) + u(i+1,j))
!
            dv=vi_v-v(i,j,N(ng),nrhs)
# ifdef ICESHEET
            du=0.5_r8*(Si(i  ,j,lsnew,isIuis)+                          &
     &                 Si(i+1,j,lsnew,isIuis))-                         &
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
# elif defined ICESHELF_QDRAG
            svstr(i,j)=0.5_r8*(rdrag2_isf(i,j-1)+rdrag2_isf(i,j))*      &
     &                 spd*dv
# elif defined ICESHELF_LDRAG
            svstr(i,j)=0.5_r8*(rdrag_isf(i,j-1)+rdrag_isf(i,j))*dv
# endif
          ELSE
            svstr(i,j)=0.0_r8
          END IF
        END DO
      END DO
!
!  Exchange halo of the ice-shelf drag stress so that downstream
!  consumers read consistent U/V-point values at tile boundaries.
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
!  Step 2: Zero out surface tracer fluxes over the ice shelf cavity.
!
!  Over cavity water (rmask_isf=1) the ocean surface is the ice base,
!  not the atmosphere, so atmospheric heat, salt, and passive-tracer
!  fluxes must be suppressed.  The 3EQN thermodynamics in isf_vbc
!  overwrites stflx(:,:,itemp/isalt) for cavity cells; other tracers
!  remain zeroed here.  Shortwave radiation is similarly masked.
!=======================================================================
!
      DO itrc=1,NT(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            stflx(i,j,itrc)=stflx(i,j,itrc)*(1.0_r8-rmask_isf(i,j))
          END DO
        END DO
      END DO
# ifdef SHORTWAVE
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          srflx(i,j)=srflx(i,j)*(1.0_r8-rmask_isf(i,j))
        END DO
      END DO
# endif
!
# ifdef ICESHELF_THERMO
!=======================================================================
!  Step 3: Compute friction velocity u* at RHO points.
!
!  u* = gamU*|tau|^(1/2) = gamU*(tau_x^2 + tau_y^2)^(1/4), valid for all
!  drag schemes: QDRAG gives u* = sqrt(rdrag2_isf) |du|, LOGDRAG gives
!  u* = sqrt(Cd_log) |du|, and LDRAG gives u* = sqrt(rdrag_isf |du|).
!=======================================================================
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (rmask_isf(i,j).eq.1.0_r8) THEN
!
!  Average sustr from surrounding U-points and svstr from V-points.
!
            cff1=0.5_r8*(sustr(i,j)+sustr(i+1,j))
            cff2=0.5_r8*(svstr(i,j)+svstr(i,j+1))
            ustar_rho=gamU(ng)*SQRT(SQRT(cff1*cff1+cff2*cff2))
            Fi(i,j,ifIust)=MAX(min_ustar(ng), ustar_rho)
          ELSE
            Fi(i,j,ifIust)=0.0_r8
          END IF
        END DO
      END DO
# endif
# ifdef ICESHEET
!
!=======================================================================
!  Step 4: Ice sheet basal drag coefficient (SSA coupling window).
!
!  The SSA momentum equation is solved every N_ssa ocean timesteps
!  (iceshelf.F), so the basal drag it needs is evaluated from the
!  ocean state averaged over the coupling window.  The accumulator
!  Ci(:,:,icIzta) is owned by isf_coupling, called at the top of
!  iceshelf.F: at the last step of each window (MOD(iic,N_ssa)=0,
!  the same cadence as the solver call) it holds the finalized
!  window mean; on the cold-start first step (iic=1, where the
!  driver also solves for the initial velocity) it holds the seed
!  sample, i.e. the instantaneous free surface.  Both are read
!  directly here.  This routine runs before the lsold/lsnew swap
!  and the isf_ssa call in iceshelf.F, so lsnew is the previous
!  SSA solution: sliding velocity and thickness come from that
!  time level, matching the solver's frozen-coefficient (Picard)
!  linearization of the drag term.
!
!  Float state follows the isf_masks.F flotation criterion evaluated
!  with the window-averaged free surface:
!    grounded:  rho_i * H  >  rho0 * (zeta_avg + h - Dcrit)
!    floating:  rho_i * H <= rho0 * (zeta_avg + h - Dcrit)
!  so the drag, the effective pressure (Coulomb law) and the float
!  state all see the same averaged ocean state.  The window phase is
!  keyed to the absolute timestep, so restart files carry partial
!  windows across; the first window of a cold start is seeded from
!  its own first sample in isf_coupling.
!
!  Sliding laws (CPP-selected, Weertman is the default):
!    default            Weertman power law
!                       beta = C^2 |u|^(m-1)
!    ICESHEET_COULOMB   regularized Coulomb (Tsai et al. 2015):
!                       H_eff = rho_i*H - rho0*D_avg,
!                       beta = C^2 |u|^(m-1) * f_c*H_eff /
!                              (C^(2/m)|u| + (f_c*H_eff)^(1/m))^m
!  with C=C_basal, m=n_basal, |u| the 4-corner PSI-average sliding
!  speed regularized by eps_u_min, and f_c=f_coulomb.  H_eff is the
!  effective pressure scaled by 1/g (the g factors cancel in beta);
!  the Coulomb form reduces to Weertman as |u|->0 and saturates at
!  the Coulomb limit tau_b = f_c*g*H_eff as |u|->infinity.  Over
!  floating and ice-free cells the drag is zeroed: the Coulomb
!  formula clamps H_eff at zero (the drag vanishes continuously as
!  H_eff -> 0+), and without ICESHEET_SMOOTH_GL the Weertman value
!  is zeroed a posteriori wherever H_eff <= 0.
!
!  With ICESHEET_SMOOTH_GL the drag is scaled by a cosine-ramped
!  grounded fraction gfrac (STREAMICE float_frac,
!  phi_gl_streamice.F):
!
!    H_float = (rho0/rho_i) * D_avg  flotation thickness
!    x       = H - H_float + 0.5*W,  W = gl_float_width
!    gfrac   = 0                        x <= 0
!            = 1                        x >= W
!            = 0.5*(1-cos(pi*x/W))      otherwise
!
!  gfrac is 1 well grounded, 0 floating, 0.5 exactly at the binary
!  grounding line, and W=0 degenerates to the binary criterion
!  (bit-identical behavior).  The Weertman law is evaluated over all
!  ice cells and scaled by gfrac, so the drag ramps continuously
!  from zero at the floating edge of the band to the full law at
!  H_float+0.5*W (STREAMICE float_frac semantics); the Coulomb drag
!  additionally vanishes with H_eff on the grounded side.  The gfrac
!  diagnostic is stored for all ice cells in Fi(:,:,ifIgrf) (netCDF
!  is_gfrac) at the SSA cadence.
!=======================================================================
!
      IF ((iic(ng).eq.1).or.(MOD(iic(ng),N_ssa(ng)).eq.0)) THEN
!
        C_sq=C_basal(ng)*C_basal(ng)
!
!  Compute beta at RHO points from the window-averaged ocean state.
!  Interior loop only: the 4-corner velocity stencil needs the
!  i+1/j+1 neighbors.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            zeta_avg=Ci(i,j,icIzta)
            D_avg=h(i,j)+zeta_avg-Dcrit(ng)
!
            Fi(i,j,ifIbta)=0.0_r8
#  ifdef ICESHEET_SMOOTH_GL
            Fi(i,j,ifIgrf)=0.0_r8
#  endif
            IF (Si(i,j,lsnew,isIhis).gt.0.0_r8) THEN
!
!  Sliding speed at the RHO point from the previous SSA solution
!  (4-corner PSI average, the same stencil the solver uses for
!  beta).
!
              u_mid=0.25_r8*(Si(i  ,j  ,lsnew,isIuis)+                  &
     &                       Si(i+1,j  ,lsnew,isIuis)+                  &
     &                       Si(i  ,j+1,lsnew,isIuis)+                  &
     &                       Si(i+1,j+1,lsnew,isIuis))
              v_mid=0.25_r8*(Si(i  ,j  ,lsnew,isIvis)+                  &
     &                       Si(i+1,j  ,lsnew,isIvis)+                  &
     &                       Si(i  ,j+1,lsnew,isIvis)+                  &
     &                       Si(i+1,j+1,lsnew,isIvis))
              u_slide=SQRT(u_mid**2+v_mid**2+eps_u_min(ng)**2)
!
              H_eff=rho_i*Si(i,j,lsnew,isIhis)-rho0*D_avg
!
#  ifdef ICESHEET_COULOMB
!
!  Effective pressure from full hydraulic connectivity to the ocean.
!  Negative H_eff (floating and ice-free cells) is clamped at zero,
!  which zeroes the drag and matches the continuous limit of the
!  formula as H_eff -> 0+.  With ICESHEET_SMOOTH_GL the drag is
!  additionally scaled by the grounded fraction gfrac.
!
              cff2=f_coulomb(ng)*MAX(H_eff, 0.0_r8)
              Fi(i,j,ifIbta)=C_sq*cff2*u_slide**(n_basal(ng)-1.0_r8)/   &
     &                       ((C_sq**(1.0_r8/n_basal(ng)))*u_slide+     &
     &                        cff2**(1.0_r8/n_basal(ng)))**n_basal(ng)
#  else
!
!  Weertman power law, evaluated over all ice cells.  With
!  ICESHEET_SMOOTH_GL the value is scaled below by the grounded
!  fraction gfrac (continuous ramp through the grounding line);
!  without it the value is zeroed below wherever H_eff <= 0.
!
              Fi(i,j,ifIbta)=C_sq*u_slide**(n_basal(ng)-1.0_r8)
#  endif
#  ifdef ICESHEET_SMOOTH_GL
!
!  Grounded fraction:  cosine ramp of H about the flotation thickness
!  H_float=(rho0/rho_i)*D_avg; gfrac=0.5 coincides with the binary
!  grounding line and W=0 degenerates to the binary criterion.
!
              H_float=(rho0/rho_i)*D_avg
              cff=Si(i,j,lsnew,isIhis)-                                 &
     &            H_float+0.5_r8*gl_float_width(ng)
              IF (cff.le.0.0_r8) THEN
                Fi(i,j,ifIgrf)=0.0_r8
              ELSE IF (cff.ge.gl_float_width(ng)) THEN
                Fi(i,j,ifIgrf)=1.0_r8
              ELSE
                Fi(i,j,ifIgrf)=0.5_r8*                                  &
     &                         (1.0_r8-COS(pi*cff/gl_float_width(ng)))
              END IF
!
              Fi(i,j,ifIbta)=Fi(i,j,ifIbta)*Fi(i,j,ifIgrf)
#  else
!
!  Binary flotation switch (no grounding line smoothing):  zero the
!  drag over floating and ice-free cells (for the Coulomb law this
!  is already guaranteed by the clamp above).
!
              IF (H_eff.le.0.0_r8) THEN
                Fi(i,j,ifIbta)=0.0_r8
              END IF
#  endif
            END IF
          END DO
        END DO
      END IF
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fi(:,:,ifIbta))
!
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Fi(:,:,ifIbta))
#  endif
# endif
!
      RETURN
      END SUBROUTINE isf_stress_tile
!
#endif
      END MODULE isf_stress_mod
