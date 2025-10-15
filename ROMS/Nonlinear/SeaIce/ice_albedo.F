#include "cppdefs.h"

      MODULE ice_albedo_mod

#if defined ICE_MODEL && defined ICE_ALBEDO
# define ICE_ALB_EC92
# define ICE_ALB_SANGLE
# define ICE_ALB_BRIEGLEB
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2025 The ROMS Group               Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module computes the surface albedo over seawater, snow, or ice !
!  using selected formulation.                                         !
!                                                                      !
!  Additionally, it computes the fraction of shortwave radtion         !
!    penetrating ice base (Frsw), so that shortwave radiation          !
!    reaching ice base is:                                             !
!                                                                      !
!      SW_ice_base = SW*(1-albedo_ice)*Frsw                            !
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
!    Iqbal, M. ,1983 An Introduction to Solar Radiation [Book]         !
!                                                                      !
!    Taylor, J.P., et al., 1996: Studies with a flexible new radiation !
!      code: Comarisions with aircraft short-wave observations         !
!                                                                      !
!    Enomoto T., 2007: Ocean surface albedo in AFES                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
# ifdef ICE_LMD_SHORTWAVE
      USE mod_mixing
# endif
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
      CALL wclock_on (ng, model, 42, __LINE__, MyFile)
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
     &                      FORCES(ng) % srflx,                         &
     &                      ICE(ng) % Si,                               &
     &                      ICE(ng) % Fi)
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
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
# endif
     &                            srflx, Si, Fi)
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
      real(r8), intent(in) :: srflx(LBi:,LBj:)
      real(r8), intent(in) :: Si(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: Fi(LBi:,LBj:,:)
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
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Si(LBi:UBi,LBj:UBj,2,nIceS)
      real(r8), intent(inout) :: Fi(LBi:UBi,LBj:UBj,nIceF)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, li_stp
!
!  Default ocean surface albedo
!
      real(r8), parameter :: alb_o  = 0.060_r8
      real(r8), parameter :: alb_od = 0.065_r8
!
!   Visual fraction of total shortwave (Aumont et al., 2015)
!
      real(r8), parameter :: fv     = 0.43_r8
!
      real(r8), parameter :: tol    = 0.000001_r8
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
!   Fraction of penetrating solar rad (visible); CESM2/CICE5.1 default
!     (cf. 0.7 in Briegleb et al., 2004)
!
      real(r8), parameter :: i0vis     = 0.70_r8

      real(r8) :: fh, alb, alvni, alini, dTs, fT
      real(r8) :: alvns, alins, asnow, alvn, alin
# elif defined ICE_ALB_EC92
      real(r8) :: alb_ice, alb_snow
# else
      real(r8), parameter :: alb_i_dry = 0.65_r8
      real(r8), parameter :: alb_i_wet = 0.60_r8
      real(r8), parameter :: alb_s_dry = 0.85_r8
      real(r8), parameter :: alb_s_wet = 0.72_r8
# endif

      real(r8) :: i_thick, s_thick
      real(r8) :: hour, yday, dangle, hangle
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
      CALL caldate (tdays(ng), yd_dp=yday, h_dp=hour)
      dangle = 2.0_r8*pi*yday/365._r8
      dangle = .006918_r8 +                                             &
     &  .070257_r8*SIN(      dangle) - .399912_r8*COS(      dangle) +   &
     &  .000907_r8*SIN(2._r8*dangle) - .006758_r8*COS(2._r8*dangle) +   &
     &  .001480_r8*SIN(3._r8*dangle) - .002697_r8*COS(3._r8*dangle)
      hangle = pi*(12.0_r8-hour)/12.0_r8
# endif

      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
# ifdef ICESHELF
          IF (zice(i,j).ne.0.0_r8) THEN
!
!  Zero out fluxes under iceshelf
!
            Fi(i,j,icAlbi) = 0.0_r8
            Fi(i,j,icAlbo) = 0.0_r8
            Fi(i,j,icFrsw) = 0.0_r8
            Fi(i,j,icSics) = 0.0_r8
            Fi(i,j,icSicb) = 0.0_r8
            Fi(i,j,icSocn) = 0.0_r8
#  ifdef ICE_LMD_SHORTWAVE
            Fi(i,j,icRfac) = 0.0_r8
#  endif
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
!    (absorbed SW = (1-albedo_ice)*SW)
!
            Fi(i,j,icAlbi) = fv*alvn + (1.0_r8-fv)*alin
          ELSE
            Fi(i,j,icAlbi) = alb_o
          END IF
# elif defined ICE_ALB_EC92
!
!  Ice and snow albedo is calculated from Ebert and Curry (1992)
!
          IF (Si(i,j,li_stp,isAice) .gt. min_ai(ng)) THEN
!
!  The threshold test can still lead to a negative albedo, so we modify
!  the minimum ice thickness for this estimate to something more
!  reasonable (SMD).
!
            i_thick=MAX(i_thick, 0.01_r8)
            IF (i_thick.ge.2.0_r8) THEN
              alb_ice =                              0.561632_r8
            ELSE IF (i_thick.ge.1.0_r8) THEN
              alb_ice = 0.07616_r8  *     i_thick  + 0.414492_r8
            ELSE
              alb_ice = 0.082409_r8 * LOG(i_thick) + 0.485472_r8
            END IF
!
!  Approximated values for snow albedo depending on COSZ, but small
!  variation.
!
            IF (Si(i,j,li_stp,isHmel).gt.0.0_r8) THEN
              IF (s_thick.ge.0.1_r8) THEN
                alb_snow = 0.701009_r8
              ELSE
                alb_snow = alb_ice+(0.701009_r8-alb_ice)*s_thick/0.1_r8
              END IF
            ELSE
              alb_snow = 0.83_r8
            END IF
!
!  Ebert and Curry estimate the melt pond albedo to be between 0.14 and
!  0.26, depending on the thickness of the meltwater ponds.  But of
!  course, the meltwater ponds do not cover all the ice surface, except
!  maybe effectively when the ice consists of small pieces with little
!  freeboard. So let's try the idea that the effective albedo of the
!  ice/melt pond/snow is that of ice for bare ice, that of snow
!  (potentially modified by meltwater if snow is present), and the
!  average of the high-end melt pond (0.26) and the low-end ice bare
!  ice (0.56) when melt ponds of more than a few centimeters deep are
!  present (SMD).
!
            IF (Si(i,j,li_stp,isHsno).gt.0.0_r8) THEN
              Fi(i,j,icAlbi) = alb_snow
            ELSE IF (Si(i,j,li_stp,isHmel).gt.0.02_r8) THEN
              Fi(i,j,icAlbi) = 0.42_r8
            ELSE
              Fi(i,j,icAlbi) = alb_ice
            END IF
          ELSE
            Fi(i,j,icAlbi) = alb_o
          END IF
# else
          cff1 = alb_s_wet - alb_s_dry
          cff2 = alb_i_wet - alb_i_dry
          IF (Si(i,j,li_stp,isAice) .gt. min_ai(ng)) THEN
            IF (s_thick .gt. 0._r8) THEN
              IF (sfc_temp .gt. -1.0_r8) THEN
                Fi(i,j,icAlbi) = alb_s_dry + cff1*(sfc_temp+1.0_r8)
              ELSE
                Fi(i,j,icAlbi) = alb_s_dry
              END IF
            ELSE
              IF (sfc_temp .gt. -1.0_r8) THEN
                Fi(i,j,icAlbi) = alb_i_dry + cff2*(sfc_temp+1.0_r8)
              ELSE
                Fi(i,j,icAlbi) = alb_i_dry
              END IF
            END IF
          ELSE
            Fi(i,j,icAlbi) = alb_o
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
          Fi(i,j,icAlbo) = 1._r8 -                                      &
     &                     ((1.0_r8-cloud(i,j))*(1.0_r8-alb_sr) +       &
     &                              cloud(i,j) *(1.0_r8-alb_od))
# elif defined ICE_ALB_SANGLE
!
! ----------------------------------------------------------------------
!  Compute ocean surface albedo
! ----------------------------------------------------------------------
!
!  Calculate Cosine of zenith angle: Iqbal (1983) Pg 15
!
          sin2   = SIN(latr(i,j)*deg2rad)*SIN(dangle)
          cos2   = COS(latr(i,j)*deg2rad)*COS(dangle)
!
!  Correction for local time if time is given in UTC
!
          hangle = hangle - lonr(i,j)*deg2rad
!
!  Cosine of zenith angle
!
          cosz   = MIN(1._r8, MAX(0._r8, sin2+cos2*COS(hangle)))
#  ifdef ICE_ALB_BRIEGLEB
!
!  Albedo of direct sun light: Briegleb et al. (1986) (default)
!
          alb_sr = 0.026_r8 / (       cosz**1.7_r8 + 0.065_r8) +        &
     &             0.150_r8 * (cosz-0.1_r8)*(cosz-0.5_r8)*(cosz-1.0_r8)
#  elif defined ICE_ALB_ENOMOTO
!
!  Albedo of direct sun light: Enomoto et al. (2007)
!
          alb_sr = 0.026_r8 / (1.1_r8*cosz**1.7_r8 + 0.065_r8) +        &
     &             0.150_r8 * (cosz-0.1_r8)*(cosz-0.5_r8)*(cosz-1.0_r8)
#  elif defined ICE_ALB_TAYLOR
!
!  Albedo of direct sun light: Taylor et al. (1996)
!
          alb_sr = 0.037_r8 / (1.1_r8*cosz**1.4_r8 + 0.150_r8)
#  endif
!
!  Albedo weighted by cloud fraction
!
          Fi(i,j,icAlbo) = 1._r8 -                                      &
     &                     ((1.0_r8-cloud(i,j))*(1.0_r8-alb_sr) +       &
     &                              cloud(i,j) *(1.0_r8-alb_od))
# else
          Fi(i,j,icAlbo) = alb_o
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
            Fi(i,j,icFrsw) = (1.0_r8-alvni)*(1.0_r8-asnow)*i0vis*fv /   &
     &                       (1.0_r8-Fi(i,j,icAlbi))
# else
            IF (i_thick .le. 0.1_r8) THEN
              cff = i_thick*5._r8
            ELSE
              cff = 0.4_r8 + i_thick
            END IF
            cff = cff + s_thick*20._r8
!
!  Assume only visual part penetrates (Briegleb et al., 2004)
!
            Fi(i,j,icFrsw) = fv*EXP(-cff)
# endif
          ELSE
            Fi(i,j,icFrsw) = 0.0_r8
          END IF
!
! ----------------------------------------------------------------------
!  Calculate net downward short-wave radiation at ocean surface, ice
!    surface, and ice base
! ----------------------------------------------------------------------
!
          Fi(i,j,icSocn) = srflx(i,j)*(1.0_r8-Fi(i,j,icAlbo))
          Fi(i,j,icSics) = srflx(i,j)*(1.0_r8-Fi(i,j,icAlbi))
          Fi(i,j,icSicb) = Fi(i,j,icSics)*Fi(i,j,icFrsw)
# ifdef ICE_LMD_SHORTWAVE
!
! ----------------------------------------------------------------------
!  Calculate correction factor for lmd_r1
! ----------------------------------------------------------------------
!
          cff1           =  1.0_r8 - Fi(i,j,icAlbo)
          cff2           = (1.0_r8 - Fi(i,j,icAlbi))*Fi(i,j,icFrsw)
          Fi(i,j,icRfac) =                                              &
     &      (1.0_r8-Si(i,j,li_stp,isAice))*cff1 /                       &
     &      MAX(tol, ((1.0_r8-Si(i,j,li_stp,isAice))*cff1 +             &
     &                        Si(i,j,li_stp,isAice) *cff2 *             &
     &                        INT(MIXING(ng)%Jwtype(i,j))))
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
     &                          Fi(:,:,icAlbo))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fi(:,:,icAlbi))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fi(:,:,icFrsw))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fi(:,:,icSics))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fi(:,:,icSicb))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fi(:,:,icSocn))
# ifdef ICE_LMD_SHORTWAVE
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fi(:,:,icRfac))
# endif
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 3,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Fi(:,:,icAlbo), Fi(:,:,icAlbi),               &
     &                    Fi(:,:,icFrsw))
      CALL mp_exchange2d (ng, tile, model, 3,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Fi(:,:,icSics), Fi(:,:,icSicb),               &
     &                    Fi(:,:,icSocn))
#  ifdef ICE_LMD_SHORTWAVE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Fi(:,:,icRfac))
#  endif
# endif
!
      RETURN
      END SUBROUTINE ice_albedo_tile
#endif
      END module ice_albedo_mod
