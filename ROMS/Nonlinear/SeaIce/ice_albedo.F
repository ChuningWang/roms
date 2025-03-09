#include "cppdefs.h"
      MODULE ice_albedo_mod

#if defined ICE_MODEL && defined ICE_ALBEDO
#define ICE_ALB_EC92
#define ICE_ALB_SANGLE
#define ICE_ALB_GREGORIAN
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2024 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module computes the surface albedo over seawater, snow, or ice !
!  using selected formulation.                                         !
!                                                                      !
!  Additionally, it computes the fraction of shortwave radtion         !
!    penetrating ice base (frswp), so that shortwave radiation         !
!    reaching ice base is:                                             !
!                                                                      !
!      SW_ice_base = SW*(1-albedo_ice)*frswp                           !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Briegleb, B.P., P. Minnis, V. Ramanathan, and E. Harrison,        !
!      1986: Comparison of Regional Clear-Sky Albedos Inferred         !
!      from Satellite Observations and Model Computations, J.          !
!      Climate and Appled Meteor., 25, 214-226.                        !
!                                                                      !
!    Briegleb, B.P., et al., 2004: Scientific Description of the Sea   !
!      Ice Component in the Community Climate System Model, Version    !
!      Three                                                           !
!                                                                      !
!    Ebert, E.E. and J.A. Curry, 1992: A parameterization of ice       !
!      cloud optical properties for climate models, J. Geophys.        !
!      Res., 97, 3831-3836, doi: doi.org/10.1029/91jD02472.            !
!                                                                      !
!    Iqbal, M.,1983 An Introduction to Solar Radiation [Book]          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ice
      USE mod_scalars
!
      USE dateclock_mod,   ONLY : caldate
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
      implicit none
!
      PUBLIC  :: ice_albedo
      PRIVATE
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_albedo (ng, tile, model)
!***********************************************************************
!
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 42, __LINE__, MyFile)
# endif
      CALL ice_albedo_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold(ng), linew(ng),                       &
# ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
# endif
# ifdef ICE_ALB_SANGLE
     &                      GRID(ng) % latr,                            &
     &                      GRID(ng) % lonr,                            &
#  ifdef CLOUDS
     &                      FORCES(ng) % cloud,                         &
#  endif
# endif
     &                      ICE(ng) % Fi,                               &
     &                      ICE(ng) % Si,                               &
     &                      FORCES(ng) % srflx,                         &
     &                      FORCES(ng) % albedo_ice,                    &
     &                      FORCES(ng) % albedo_ocn,                    &
     &                      FORCES(ng) % frswp,                         &
# ifdef ICE_LMD_SHORTWAVE
     &                      FORCES(ng) % rfaci,                         &
# endif
     &                      FORCES(ng) % srflx_ics,                     &
     &                      FORCES(ng) % srflx_icb,                     &
     &                      FORCES(ng) % srflx_ocn)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 17, __LINE__, __FILE__)
# endif
      RETURN
      END SUBROUTINE ice_albedo
!
!***********************************************************************
      SUBROUTINE ice_albedo_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            liold, linew,                         &
# ifdef ICESHELF
     &                            zice,                                 &
# endif
# ifdef ICE_ALB_SANGLE
     &                            latr, lonr,                           &
#  ifdef CLOUDS
     &                            cloud,                                &
#  endif
     &                            Fi, Si, srflx,                        &
     &                            albedo_ice, albedo_ocn, frswp,        &
# ifdef ICE_LMD_SHORTWAVE
     &                            rfaci,                                &
# endif
     &                            srflx_ics, srflx_icb, srflx_ocn)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: liold, linew
!
# ifdef ASSUMED_SHAPE
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#  endif
#  ifdef ICE_ALB_SANGLE
      real(r8), intent(in) :: latr(LBi:,LBj:)
      real(r8), intent(in) :: lonr(LBi:,LBj:)
#   ifdef CLOUDS
      real(r8), intent(in) :: cloud(LBi:,LBj:)
#   endif
#  endif
      real(r8), intent(in) :: Fi(LBi:,LBj:,:)
      real(r8), intent(in) :: Si(LBi:,LBj:,:,:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
      real(r8), intent(out) :: albedo_ice(LBi:,LBj:)
      real(r8), intent(out) :: albedo_ocn(LBi:,LBj:)
      real(r8), intent(out) :: frswp(LBi:,LBj:)
#  ifdef ICE_LMD_SHORTWAVE
      real(r8), intent(out) :: rfaci(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: srflx_ics(LBi:,LBj:)
      real(r8), intent(out) :: srflx_icb(LBi:,LBj:)
      real(r8), intent(out) :: srflx_ocn(LBi:,LBj:)

# else
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ICE_ALB_SANGLE
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: lonr(LBi:UBi,LBj:UBj)
#   ifdef CLOUDS
      real(r8), intent(in) :: cloud(LBi:UBi,LBj:UBj)
#   endif
#  endif
      real(r8), intent(in) :: Fi(LBi:UBi,LBj:UBj,nIceF)
      real(r8), intent(in) :: Si(LBi:UBi,LBj:UBj,2,nIceS)
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: albedo_ice(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: albedo_ocn(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: frswp(LBi:UBi,LBj:UBj)
#  ifdef ICE_LMD_SHORTWAVE
      real(r8), intent(out) :: rfaci(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: srflx_ics(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: srflx_icb(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: srflx_ocn(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, li_stp
!
!  Default ocean surface albedo
!
      real(r8), parameter :: alb_o  = 0.060_r8
      real(r8), parameter :: alb_oz = 0.065_r8
!
# ifdef ICE_ALB_CSIM
!
! DuVivier et al., 2018
!
!   Visible ice albedo for h > ahmax; CESM2/CICE5.1 default
!     (cf. 0.73 in Briegleb et al., 2004)
!
      real(r8), parameter :: albicev   = 0.78_r8
!
!   Near-ir ice albedo for h > ahmax; CESM2/CICE5.1 default
!     (cf. 0.33 in Briegleb et al., 2004)
!
      real(r8), parameter :: albicei   = 0.36_r8
!
!   Cold snow albedo, visible; CESM2/CICE5.1 default
!     (cf. 0.96 in Briegleb et al., 2004)
!
      real(r8), parameter :: albsnowv  = 0.98_r8
!
!   Cold snow albedo, near IR; CESM2/CICE5.1 default
!     (cf. 0.68 in Briegleb et al., 2004)
!
      real(r8), parameter :: albsnowi  = 0.70_r8
!
!   Thickness above which ice albedo is const (m); CESM2/CICE5.1 default
!     (cf. 0.5 m in Briegleb et al., 2004)
!
      real(r8), parameter :: ahmax     = 0.3_r8
!
!   Parameter for fractional snow area (m); CESM2/CICE5.1 default
!     (cf. 0.02 m in Briegleb et al., 2004)
!
      real(r8), parameter :: snowpatch = 0.02_r8
!
!   Change in temp to give dalb_mlt albedo (C); CESM2/CICE5.1 default
!     (cf. 1.0 C in Briegleb et al., 2004)
!
      real(r8), parameter :: dT_melt   = 1.0_r8
!
!   Albedo change per dT_melt in temp for ice; CESM2/CICE5.1 default
!     (cf. -0.075 in Briegleb et al., 2004)
!
      real(r8), parameter :: dalb_mlt  = -0.075_r8
!
!   Albedo vis change per dT_melt in temp for snow; CESM2/CICE5.1 default
!     (cf. -0.1 in Briegleb et al., 2004)
!
      real(r8), parameter :: dalb_mltv = -0.1_r8
!
!   Albedo ir change per dT_melt in temp for snow; CESM2/CICE5.1 default
!     (cf. -0.15 in Briegleb et al., 2004)
!
      real(r8), parameter :: dalb_mlti = -0.15_r8
!
!   Melting temperature, ice top surface (C)
!
      real(r8), parameter :: Timelt    = 0.0_r8
!
!   Melting temperature, snow top surface (C)
!
      real(r8), parameter :: Tsmelt    = 0.0_r8
      real(r8), parameter :: puny      = 1.0e-11_r8
!
!   Visual fraction of total shortwave (Aumont et al., 2015)
!
      real(r8), parameter :: fv        = 0.43_r8
!
!   Fraction of penetrating solar rad (visible); CESM2/CICE5.1 default
!     (cf. 0.7 in Briegleb et al., 2004)
!
      real(r8), parameter :: i0vis     = 0.70_r8

      real(r8) :: fh, alb, alvni, alini, dTs, fT
      real(r8) :: alvns, alins, asnow, alvn, alin
# elif defined ICE_ALB_EC92
      real(r8) :: alb_s = 0.83_r8
# else
      real(r8), parameter :: alb_i_dry = 0.65_r8
      real(r8), parameter :: alb_i_wet = 0.60_r8
      real(r8), parameter :: alb_s_dry = 0.85_r8
      real(r8), parameter :: alb_s_wet = 0.72_r8
# endif

      real(r8) :: i_thick, s_thick
      real(r8) :: hours, days, dangle, hangle, langle
      real(r8) :: sin2, cos2, cosz, alb_sr
# ifdef ICE_ALB_DMEAN
      integer,  parameter :: ifrac     = 24
      real(r8), parameter :: absh2o    = 0.09_r8
      real(r8), parameter :: fraci     = 1._r8/REAL(ifrac,r8)
      integer  :: npart
      real(r8) :: stot, sads, ssurf, sdir, sdif, altdeg, cfac
# endif
      real(r8) :: cff, cff1, cff2, sfc_temp

# include "set_bounds.h"
!
!=======================================================================
!  Compute the surface albedo over seawater, snow or ice
!=======================================================================
!
!  Set ice model time level.
!
      IF (PerfectRST(ng) .and. iic(ng).eq.ntstart(ng)) THEN
        li_stp = liold
      ELSE
        li_stp = linew
      END IF

# if defined ICE_ALB_SANGLE || defined ICE_ALB_DMEAN
!
!  Calculate solar angle.
!    Iqbal (1992) Pg 7
!
      CALL caldate (tdays(ng), yd_dp=days, h_dp=hours)
      days   = MOD(days-1._r8,365._r8)
      dangle = 2.0_r8*pi*days/365._r8
      dangle = .006918_r8 +                                             &
     &  .070257_r8*SIN(      dangle) - .399912_r8*COS(      dangle) +   &
     &  .000907_r8*SIN(2._r8*dangle) - .006758_r8*COS(2._r8*dangle) +   &
     &  .001480_r8*SIN(3._r8*dangle) - .002697_r8*COS(3._r8*dangle)
      hangle = pi*(hours-12.0_r8)/12.0_r8
# endif

      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
# ifdef ICESHELF
          IF (zice(i,j).ne.0.0_r8) THEN
!
!  Zero out fluxes under iceshelf
!
            albedo_ice(i,j) = 0.0_r8
            albedo_ocn(i,j) = 0.0_r8
            frswp(i,j)      = 0.0_r8
#  ifdef ICE_LMD_SHORTWAVE
            rfaci(i,j)      = 0.0_r8
#  endif
            srflx_ics(i,j)  = 0.0_r8
            srflx_icb(i,j)  = 0.0_r8
            srflx_ocn(i,j)  = 0.0_r8
          ELSE
# endif
!
! Calculate ice/snow surface albedo
!
          cff      = 1.0_r8 / (Si(i,j,li_stp,isAice) + 0.000001_r8)
          i_thick  = cff    *  Si(i,j,li_stp,isHice)
          s_thick  = cff    *  Si(i,j,li_stp,isHsno)
          sfc_temp = Fi(i,j,icIsst)
# ifdef ICE_ALB_CSIM
!
!  Follow ice_shortwave.F90/compute_albedos subroutine (CESM2/CICE5.1)
!
          IF (Si(i,j,li_stp,isAice) .gt. min_ai(ng)) THEN
!
!  Bare ice, thickness dependence
!    Briegleb (2004) Eqs 30, 31
!
            fh    = MIN(1.0_r8, ATAN(i_thick*4.0_r8)/ATAN(ahmax*4.0_r8))
            alb   = alb_o*(1.0_r8-fh)
            alvni = alb + albicev*fh
            alini = alb + albicei*fh
!
!  Bare ice, temperature dependence
!    Briegleb (2004) Eq 29
!
            dTs   = Timelt - sfc_temp
!
!  If dTs >= dT_melt, fT = 0; if dTs == 0, fT=-1;
!    else follow Briegleb (2004) Eqs 29
!
            fT    = MIN(0.0_r8, dTs/dT_melt-1.0_r8)
!
!  Ice albedo decreases with warming, capped by ocean albedo
!    Briegleb (2004) Eq 33
!
!  Note: the CESM2/CICE5.1 code (ice_shortwave.F90/compute_albedos) does
!        not impose lower cap here, but it seems safer to do so for the
!        case of tis>>Timelt (=> fT<<0)
!
            alvni = MAX(alb_o, alvni - dalb_mlt*fT)
            alini = MAX(alb_o, alini - dalb_mlt*fT)
!
            IF (s_thick > puny) THEN
!
!  Snow on ice, temperature dependence, capped by ocean albedo
!    Briegleb (2004) Eq 29
!
              alvns = MAX(alb_o, albsnowv - dalb_mltv*fT)
              alins = MAX(alb_o, albsnowi - dalb_mlti*fT)
!
!  Combined ice/snow visual and near IR albedo
!    Briegleb (2004) Eqs 34, 35
!
              asnow = s_thick / (s_thick + snowpatch)
              alvn  = alvni*(1.0_r8-asnow) + alvns*asnow
              alin  = alini*(1.0_r8-asnow) + alins*asnow
            ELSE
              asnow = 0.0_r8
              alvn  = alvni
              alin  = alini
            END IF
!
!  Combined ice/snow albedo with respect to total incident SW
!    (absorbed SW = (1-albedo_i)*SW)
!
            albedo_ice(i,j) = fv*alvn + (1.0_r8-fv)*alin
          ELSE
            albedo_ice(i,j) = alb_o
          END IF
# elif defined ICE_ALB_EC92
!
!  Ice and snow albedo is calculated from Ebert and Curry (1992)
!
          IF (Si(i,j,li_stp,isAice) .gt. min_ai(ng)) THEN
            IF (s_thick .gt. 0._r8) THEN
              albedo_ice(i,j) = alb_s
            ELSE
!
!  Approximated values for ice alb (depend on COSZ, but small variation)
!
              i_thick = MAX(i_thick, 0.000001_r8)
              IF (i_thick .ge. 2._r8) THEN
                albedo_ice(i,j) = 0.561632_r8
              ELSE IF (i_thick .ge. 1._r8) THEN
                albedo_ice(i,j) = 0.07616_r8 *    i_thick  + 0.414492_r8
              ELSE
                albedo_ice(i,j) = 0.082409_r8*LOG(i_thick) + 0.485472_r8
                albedo_ice(i,j) = MAX(alb_o, albi(i,j))
              END IF
            END IF
          ELSE
            albedo_ice(i,j) = alb_o
          END IF
# else
          cff1 = alb_s_wet - alb_s_dry
          cff2 = alb_i_wet - alb_i_dry
          IF (Si(i,j,li_stp,isAice) .gt. min_ai(ng)) THEN
            IF (s_thick .gt. 0._r8) THEN
              IF (sfc_temp .gt. -1.0_r8) THEN
                albedo_ice(i,j) = alb_s_dry + cff1*(sfc_temp+1.0_r8)
              ELSE
                albedo_ice(i,j) = alb_s_dry
              END IF
            ELSE
              IF (sfc_temp .gt. -1.0_r8) THEN
                albedo_ice(i,j) = alb_i_dry + cff2*(sfc_temp+1.0_r8)
              ELSE
                albedo_ice(i,j) = alb_i_dry
              END IF
            END IF
          ELSE
            albedo_ice(i,j) = alb_o
          END IF
# endif

# ifdef ICE_ALB_DMEAN
!
! ----------------------------------------------------------------------
!  compute 24 hrs time-averaged albedo at ocean surface
! ----------------------------------------------------------------------
!
          sin2 = SIN(latr(i,j)*deg2rad)*SIN(dangle)
          cos2 = COS(latr(i,j)*deg2rad)*COS(dangle)
          stot = 0._r8
          sads = 0._r8
          DO npart=1,ifrac
!
!  Hour angle, in radians
!
            hangle = pi*(REAL(npart,r8)/12._r8-1._r8)
!
!  Cosine of the zenith angle
!
            cosz   = MIN(1._r8, MAX(0._r8, sin2+cos2*COS(hangle)))
!
!  Direct radiation component
!
            sdir   = 0.7_r8**MIN(100._r8, 1._r8/(cosz+tol))
!
!  Diffusive radiation component
!
            sdif   = ((1._r8-absh2o)-sdir)*.5_r8
!
!  Solar noon altitude in degrees
!
            altdeg = MAX(0._r8, ASIN(sin2+cos2))*rad2deg
!
!  Cloudiness correction
!
            cfac   = 1._r8 - 0.62_r8*cloud(i,j) + 0.0019_r8*altdeg
            ssurf  = (sdir + sdif)*cfac+tol
            stot   = stot + ssurf
!
!  Albedo for direct light
!
            alb_sr = 0.03_r8*EXP(0.742_r8*ACOS(cosz)**2.866_r8)
!
!  Radiation weighted sum of direct albedo
!
            sads   = sads + alb_sr*ssurf
          END DO
!
!  Daily mean co-albedo over ocean surface
!
          alb_sr = sads/stot
!
!  Daily mean albedo weighted by cloud fraction
!
          albedo_ocn(i,j) = 1._r8 -                                     &
     &                      ((1.0_r8-cloud(i,j))*(1.0_r8-alb_sr) +      &
     &                               cloud(i,j) *(1.0_r8-alb_oz))
# elif defined ICE_ALB_SANGLE
!
! ----------------------------------------------------------------------
!  Compute ocean surface albedo: Iqbal (1983) Pg 15
! ----------------------------------------------------------------------
!
          sin2   = SIN(latr(i,j)*deg2rad)*SIN(dangle)
          cos2   = COS(latr(i,j)*deg2rad)*COS(dangle)
#  ifdef ICE_ALB_GREGORIAN
!
!  Correction for local time if Gregorian reference time is used
!
          langle = MOD(lonr(i,j), 360._r8)*rad2deg
#  else
          langle = 0.0_r8
#  endif
!
!  Cosine of zenith angle
!
          cosz   = MIN(1._r8, MAX(0._r8, sin2+cos2*COS(hangle+langle)))
!
!  Albedo of direct sun light
!
          alb_sr = 0.03_r8*EXP(0.742_r8*ACOS(cosz)**2.866_r8)
!
!  Albedo weighted by cloud fraction
!
          albedo_ocn(i,j) = 1._r8 -                                     &
     &                      ((1.0_r8-cloud(i,j))*(1.0_r8-alb_sr) +      &
     &                               cloud(i,j) *(1.0_r8-alb_oz))
# else
          albedo_ocn(i,j) = alb_o
# endif
!
! ----------------------------------------------------------------------
!  Calculate shortwave ice penetration factor
! ----------------------------------------------------------------------
!
          IF (Si(i,j,li_stp,isAice) .gt. min_ai(ng)) THEN
# ifdef ICE_ALB_CSIM
!
!  Note: We assume that only the visual part incident on snow-uncovered
!        ice can penetrate through the surface layer (if not reflected)
!
!  CICE 5.1 cice/source/ice_shortwave.F90 Ln 996
!
            frswp(i,j) = (1.0_r8-alvni)*(1.0_r8-asnow)*i0vis*fv /       &
     &                   (1.0_r8-albedo_ice(i,j))
# else
            IF (i_thick .le. 0.1_r8) THEN
              cff = i_thick*5._r8
            ELSE
              cff = 0.4_r8 + i_thick
            END IF
            cff = cff + s_thick*20._r8
            frswp(i,j) = fv*EXP(-cff)
# endif
          ELSE
            frswp(i,j) = 1.0_r8
          END IF
!
! ----------------------------------------------------------------------
!  Calculate net downward short-wave radiation at ocean surface, ice
!    surface, and ice base
! ----------------------------------------------------------------------
!
          srflx_ocn(i,j) = srflx(i,j)*(1.0_r8-albedo_ocn(i,j))
          srflx_ics(i,j) = srflx(i,j)*(1.0_r8-albedo_ice(i,j))
          srflx_icb(i,j) = srflx_ics(i,j)*frswp(i,j)
# ifdef ICE_LMD_SHORTWAVE
!
! ----------------------------------------------------------------------
!  Calculate correction factor for lmd_r1
! ----------------------------------------------------------------------
!
          cff1       =  1.0_r8 - albedo_ocn(i,j)
          cff2       = (1.0_r8 - albedo_ice(i,j))*frswp(i,j)
          rfaci(i,j) =           (1.0_r8 - ai(i,j,li_stp))*cff1 /       &
     &                 MAX(tol, ((1.0_r8 - ai(i,j,li_stp))*cff1 +       &
     &                                     ai(i,j,li_stp) *cff2))
# endif
# ifdef ICESHELF
          END IF
# endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          albedo_ocn)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          albedo_ice)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          frswp)
# ifdef ICE_LMD_SHORTWAVE
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rfaci)
# endif
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          srflx_ics)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          srflx_icb)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          srflx_ocn)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    albedo_ocn, albedo_ice, frswp)
#  ifdef ICE_LMD_SHORTWAVE
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rfaci)
#  endif
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    srflx_ics, srflx_icb, srflx_ocn)
# endif
      RETURN
!
      END SUBROUTINE ice_albedo_tile
#endif
      END module ice_albedo_mod
