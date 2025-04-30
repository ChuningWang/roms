      MODULE ice_thermo_mod
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2024 The ROMS/TOMS Group      Paul Budgell       !
!    Licensed under a MIT/X style license           Katherine Hedstrom !
!    See License_ROMS.md                            Scott M. Durski    !
!================================================== Hernan G. Arango ===
!                                                                      !
!  It computes the ice thermodynamic growth and decay term based on    !
!  Mellor and Kantha (1989) and Parkinson and Washington (1979).       !
!                                                                      !
!  It computes heat fluxes and ice production rates:                   !
!                                                                      !
!    Fi(i,j,icW_ai) = -(qai(i,j) - qi2(i,j)) /(hfus1(i,j)*RhoSW)       !
!                                                                      !
!  and updates the internal ice temperature:                           !
!                                                                      !
!    Si(i,j,linew,isTice) = Si(i,j,linew,isTice) + dtice*RHS           !
!                                                                      !
!  Required atmospheric fields:                                        !
!                                                                      !
!    Qnet_ai(:,:)          net heat flux at the air-ice interface      !
!    rain(:,:)             rain fall rate                              !
!    sustr(:,:)            surface u-wind stress                       !
!    svstr(:,:)            surface v-wind stress                       !
!                                                                      !
!  Required ocean fields:                                              !
!                                                                      !
!    dtice                 ice kernel timestep                         !
!    rmask(:,:)            land/sea mask                               !
!    t(:,:,N(ng),:,:)      surface ocean temperature and salinity      !
!                                                                      !
!  Required ice variables:                                             !
!                                                                      !
!    Fi(:,:,icSnof)        snow fall rate                              !
!    Fi(:,:,icQcon)        gradient coefficient for heat conductivity  !
!    Fi(:,:,icQrhs)        downward heat conductivity term             !
!    Fi(:,:,icS0mk)        salinity of molecular sublayer under ice    !
!    Fi(:,:,icT0mk)        temperature of molecular sublayer under ice !
!    Fi(:,:,icIsst)        temperature at snow-air interface           !
!    Fi(:,:,icW_ai)        rate of melt/freeze at air-ice interface    !
!    Fi(:,:,icW_ao)        rate of melt/freeze at air-ocean interface  !
!    Fi(:,:,icW_fr)        rate of ice accretion by frazil growth      !
!    Fi(:,:,icW_io)        rate of melt/freeze at ice-ocean interface  !
!    Fi(:,:,icW_ro)        rate of melt/freeze runoff into ocean       !
!    Fi(:,:,icIOmf)        ice-ocean mass flux                         !
!                                                                      !
!    Si(:,:,linew,isAice)  ice concentration                           !
!    Si(:,:,linew,isIage)  ice age                                     !
!    Si(:,:,linew,isEnth)  scaled perturbation ice heat content        !
!    Si(:,:,linew,isHice)  ice thickness, ice mass (divided by area)   !
!    Si(:,:,linew,isHsno)  snow thickness, mass snow per area          !
!    Si(:,:,linew,isHmel)  melt water thickness on ice                 !
!    Si(:,:,linew,isTice)  ice interior temperature (ice layer middle) !
!                                                                      !
!  Relevant Internal variables:                                        !
!                                                                      !
!    brnfr(:,:)            brine fraction                              !
!    hfus1(:,:)            latent heat of fusion (L_o or L_3)          !
!    ice_thick(:,:)        ice thickness                               !
!    qai(:,:)              upward air-ice heat flux                    !
!    qio(:,:)              upward ice-ocean heat flux                  !
!    qi2(:,:)              heat flux in ice                            !
!    sice(:,:)             ice salinity                                !
!    snow_thick(:,:)       snow thickness                              !
!    t2(:,:)               temperature at ice/snow interface           !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Mellor, G.L. and L. Kantha, 1989: An Ice-Ocean Coupled Model,     !
!      J. Geophys. Res., 94, 10937-10954.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_boundary
#ifdef ICE_SHOREFAST
      USE mod_coupling
#endif
      USE mod_forces
      USE mod_grid
      USE mod_ice
      USE mod_ocean
      USE mod_scalars
!
      USE bc_2d_mod,       ONLY : bc_r2d_tile
      USE bc_3d_mod,       ONLY : bc_r3d_tile
      USE ice_bc2d_mod,    ONLY : ice_bc2d_tile
      USE ice_tibc_mod,    ONLY : ice_tibc_tile
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
#endif
!
      implicit none
!
      PUBLIC  :: ice_thermo
      PRIVATE :: ice_thermo_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_thermo (ng, tile, model)
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
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, model, 42, __LINE__, MyFile)
#endif
      CALL ice_thermo_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), liold(ng), linew(ng),             &
#ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#endif
#ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
#endif
#ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
#endif
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % z_w,                             &
#ifdef ICE_SHOREFAST
     &                      GRID(ng) % h,                               &
     &                      COUPLING(ng) % Zt_avg1,                     &
#endif
     &                      OCEAN(ng) % t,                              &
     &                      FORCES(ng) % sustr,                         &
     &                      FORCES(ng) % svstr,                         &
     &                      FORCES(ng) % rain,                          &
     &                      FORCES(ng) % Qnet_ai,                       &
     &                      FORCES(ng) % stflx,                         &
     &                      ICE(ng) % Ti,                               &
     &                      ICE(ng) % Fi,                               &
     &                      ICE(ng) % Si)
#ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
#endif
!
      RETURN
      END SUBROUTINE ice_thermo
!
!***********************************************************************
      SUBROUTINE ice_thermo_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nrhs, liold, linew,                   &
#ifdef MASKING
     &                            rmask,                                &
#endif
#ifdef WET_DRY
     &                            rmask_wet,                            &
#endif
#ifdef ICESHELF
     &                            zice,                                 &
#endif
     &                            z_r, z_w,                             &
#ifdef ICE_SHOREFAST
     &                            h, Zt_avg1,                           &
#endif
     &                            t,                                    &
     &                            sustr, svstr,                         &
     &                            rain,                                 &
     &                            Qnet_ai,                              &
     &                            stflx,                                &
     &                            Ti, Fi, Si)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, liold, linew
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
# endif
# ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
# endif
# ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
# endif
# ifdef ICE_SHOREFAST
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)
# endif
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: rain(LBi:,LBj:)
      real(r8), intent(in) :: Qnet_ai(LBi:,LBj:)
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
      real(r8), intent(inout) :: Ti(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: Fi(LBi:,LBj:,:)
      real(r8), intent(inout) :: Si(LBi:,LBj:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
# ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
# endif
# ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
# endif
# ifdef ICE_SHOREFAST
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rain(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Qnet_ai(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(inout) :: Ti(LBi:UBi,LBj:UBj,2,NTICE(ng))
      real(r8), intent(inout) :: Fi(LBi:UBi,LBj:UBj,nIceF)
      real(r8), intent(inout) :: Si(LBi:UBi,LBj:UBj,2,nIceS)
#endif

! Local variable definitions
!
      integer :: i, j, ic, itrc
!
      real(r8), parameter :: eps      = 1.0E-4_r8     ! zero division
      real(r8), parameter :: sice_ref = 3.2_r8        ! [psu]
      real(r8), parameter :: MelRho   = 1003.1_r8     ! [kg m-3]
      real(r8), parameter :: AlphIc   = 2.034_r8      ! [W m-1 K-1]
      real(r8), parameter :: AlphSn   = 0.31_r8       ! [W m-1 K-1]
      real(r8), parameter :: Cp_i     = 2093.0_r8     ! [J kg-1 K-1]
      real(r8), parameter :: Cp_w     = 3990.0_r8     ! [J kg-1 K-1]
      real(r8), parameter :: frln     = -0.0543_r8    ! [psu C-1]
      real(r8), parameter :: hfus     = 3.347E+5_r8   ! [J kg-1]
      real(r8), parameter :: kappa    = 0.4_r8        ! von Karman
      real(r8), parameter :: nu       = 1.8E-6_r8     ! m2/s
      real(r8), parameter :: prs      = 2432.0_r8     ! S Schmidt Num.
      real(r8), parameter :: prt      = 13.0_r8       ! T Prandtl Num.
      real(r8), parameter :: RhoSW    = 1026.0_r8     ! [kg m-3]
      real(r8), parameter :: RhoCp    = RhoSW*Cp_w    ! [kg m-1 s-2 K-1]
      real(r8), parameter :: RhoCpr   = 1.0_r8/rhocp  ! 1/(Rho*Cp)
      real(r8), parameter :: tpr      = 0.85_r8       ! Turb. Prandtl
      real(r8), parameter :: ykf      = 3.14          ! Yaglom/Kader
!
!  (CW) In MK89, z0ii is 0.05/3 (Table 1)
!       We also add ocean roughness following MK89 PG 10,941
!
      real(r8), parameter :: z0ii = 0.01666666666667_r8 ! ice roughness
      real(r8), parameter :: z0oi = 1.317633460955_r8   ! sea roughness
!
      real(r8) :: cff, cff1, cff2, cff3
      real(r8) :: dztop, corfac, cot, delta_mi
      real(r8) :: hicehinv, hstar, ai_old, hi_old, phi
      real(r8) :: Qsur, rno, termt, terms, tfrz
      real(r8) :: xtot, z0i, z0o, z0, zdz0
      real(r8) :: hh, clear, fac_sf
!
      logical, dimension(IminS:ImaxS,JminS:JmaxS) :: overOcean
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: alph
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: brnfr
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: b2d
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: chs
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: cht
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Coa
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: hfus1
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ice_thick
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Qai
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Qio
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Qao
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Qi2
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: salt_top
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: sice
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: snow_thick
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: temp_top
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: t2
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: utau

#include "set_bounds.h"
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          overOcean(i,j)=.TRUE.
#ifdef MASKING
          overOcean(i,j)=overOcean(i,j).and.(rmask(i,j).gt.0.0_r8)
# ifdef WET_DRY
          overOcean(i,j)=overOcean(i,j).and.(rmask_wet(i,j).gt.0.0_r8)
# endif
#endif
#ifdef ICESHELF
          overOcean(i,j)=overOcean(i,j).and.(zice(i,j).ge.0.0_r8)
#endif
        END DO
      END DO
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (overOcean(i,j)) THEN
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
!  Extract ocean surface temperature, salinity, heat flux.
!
            temp_top(i,j)=t(i,j,N(ng),nrhs,itemp)
            salt_top(i,j)=t(i,j,N(ng),nrhs,isalt)
            salt_top(i,j)=MIN(MAX(0.0_r8, salt_top(i,j)), 40.0_r8)
            Qao(i,j)=-stflx(i,j,itemp)*RhoCp
!
!  Compute squared-root of surface wind stress magnitude.
!
            utau(i,j)=SQRT(SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+  &
     &                          (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2))
            utau(i,j)=MAX(utau(i,j), eps)
!
!  Compute snow thickness, ice thickness, freezing temperature, brine
!  fraction, and thermal conductivity.
!
            ice_thick(i,j)=0.05_r8+                                     &
     &        Si(i,j,linew,isHice)/MAX(Si(i,j,linew,isAice), eps)
            snow_thick(i,j)=                                            &
     &        Si(i,j,linew,isHsno)/MAX(Si(i,j,linew,isAice), eps)
            sice(i,j)=MIN(sice_ref, salt_top(i,j))
            tfrz=frln*sice(i,j)
            brnfr(i,j)=tfrz/MIN(Si(i,j,linew,isTice), -eps)
            brnfr(i,j)=MIN(MAX(brnfr(i,j), 0.0_r8), 0.2_r8)
            alph(i,j)=AlphIc*MAX(1.0_r8-1.2_r8*brnfr(i,j), 0.25_r8)
            cff=Si(i,j,linew,isHice)**2
            corfac=2.0_r8/(1.0_r8+EXP(-cff))
            alph(i,j)=alph(i,j)*corfac
            Coa(i,j)=2.0_r8*alph(i,j)*snow_thick(i,j)/                  &
     &               (AlphSn*ice_thick(i,j))
!
!  Gradient coefficient for heat conductivity term.
!
            b2d(i,j)=2.0_r8*alph(i,j)/(ice_thick(i,j)*(1.0_r8+Coa(i,j)))
            Fi(i,j,icQcon)=Fi(i,j,icQcon)+b2d(i,j)
!
!-----------------------------------------------------------------------
!  Solve for temperature at the top of the ice layer.
!-----------------------------------------------------------------------
!
            IF (Si(i,j,linew,isAice).gt.min_ai(ng)) THEN
!
!  Downward conductivity term, assuming the ocean at the freezing point
!  (need to convert ice temperature to Kelvin).
!
              Fi(i,j,icQrhs)=Fi(i,j,icQrhs)+                            &
     &                       b2d(i,j)*(Si(i,j,linew,isTice)+273.15_r8)
!
!  Compute temperature at the snow/ice interface (convert to Celsius).
!
              Fi(i,j,icIsst)=(Fi(i,j,icQrhs)/Fi(i,j,icQcon))-273.15_r8
!
!  Bound value at zero Celsius for stability. It can occassionably be
!  unstable and take "icIsst" in the wrong direction.
!
              Fi(i,j,icIsst)=MIN(MAX(Fi(i,j,icIsst), -45.0_r8), 0.0_r8)
!
!  Calculate new interior ice temperature.
!
!  The "cot" calculation derives from an assumption of a linear
!  relationship between ice temperature and salinity.  As "isTice"
!  approaches zero Celsius, "cot" goes to infinity and any change in
!  "isTice" in the timestepping below is stymied, so keep the
!  denominator here below zero (SMD).
!
              cot=-tfrz*hfus/                                           &
     &            (MIN(Si(i,j,linew,isTice), frln*sice_ref))**2+Cp_i
              cff1=IceRho(ng)*cot*ice_thick(i,j)**2
              cff2=Fi(i,j,icIsst)-(2.0_r8+Coa(i,j))*Si(i,j,linew,isTice)
              cff3=1.0_r8+Coa(i,j)
              Si(i,j,linew,isTice)=Si(i,j,linew,isTice)+                &
     &                             dtice(ng)*                           &
     &                             (2.0_r8*alph(i,j)/cff1*              &
     &                              (Fi(i,j,icT0mk)+cff2/cff3))
!
!  Ensure that "isTice" remains below "frln*sice_ref" whenever ice is
!  present. Also ensure that it remains below the maximum of either
!  the ice surface temperature or the water temperature below. This
!  is imperfect as warmer ice can advect from elsewhere and alter
!  its heat content, or thick ice might actually have a warmer
!  interior temperature than either the surface or bottom. But in
!  general it does not make sense that ice forms at a warmer
!  temperature than exists in any of its surroundings (SMD).
!
              Si(i,j,linew,isTice)=MAX(Si(i,j,linew,isTice), -35.0_r8)
              Si(i,j,linew,isTice)=MIN(Si(i,j,linew,isTice), tfrz)
              Si(i,j,linew,isTice)=MIN(Si(i,j,linew,isTice),            &
     &                                 MAX(Fi(i,j,icIsst),              &
     &                                     Fi(i,j,icT0mk)))
!
!  Compute temperature at ice/snow boundary.
!
              t2(i,j)=(Fi(i,j,icIsst)+Coa(i,j)*Si(i,j,linew,isTice))/   &
     &                (1.0_r8+Coa(i,j))
!
!  Compute heat fluxes.
!
              hicehinv=2.0_r8/ice_thick(i,j)
!
!  Net heat flux from ice interior to snow (MK89 Eq 9b).
!
              Qi2(i,j)=alph(i,j)*hicehinv*                              &
     &                 (Si(i,j,linew,isTice)-t2(i,j))
!
!  Net heat flux from ocean to ice interior (MK89 Eq 11).
!
              Qio(i,j)=alph(i,j)*hicehinv*                              &
     &                 (Fi(i,j,icT0mk)-Si(i,j,linew,isTice))
!
!  Set net heat flux from ice to atmosphere.
!
              Qai(i,j)=Qnet_ai(i,j)
            ELSE
!
!  If no ice, set state variable to default values.
!
              Fi(i,j,icIsst)=Fi(i,j,icT0mk)
              t2(i,j)=Fi(i,j,icT0mk)
              Si(i,j,linew,isTice)=-2.0_r8
              Qi2(i,j)=0.0_r8
              Qio(i,j)=0.0_r8
              Qai(i,j)=0.0_r8
            END IF
!
!-----------------------------------------------------------------------
!  Solve for ice suface phase exchange (ice melt/refreeze).
!-----------------------------------------------------------------------
!
!  Compute meltpond thickness.
!  (CW) Add notes to clarify the decision tree of Durski et al. (2019).
!
            IF (Si(i,j,linew,isAice).gt.min_ai(ng)) THEN
              cff=1.0_r8-brnfr(i,j)
              hfus1(i,j)=hfus*cff+                                      &
     &                   Fi(i,j,icIsst)*Cp_w-                           &
     &                   (cff*Cp_i+brnfr(i,j)*Cp_w)*Si(i,j,linew,isTice)
!
!  Recompute net heat flux from ice interior to snow (MK89 Eq 9a, b).
!
              Qi2(i,j)=b2d(i,j)*(Si(i,j,linew,isTice)-Fi(i,j,icIsst))
!
!  Compute net heat flux through snow layer.
!
              IF ((Si(i,j,linew,isHsno).le.eps).and.                    &
     &            (Si(i,j,linew,isHmel).le.eps)) THEN
                Qsur=-(Qai(i,j)-Qi2(i,j))/(hfus1(i,j)*RhoSW)
              ELSE IF ((Si(i,j,linew,isHsno).le.eps).and.               &
     &                 (Si(i,j,linew,isHmel).gt.eps)) THEN
                Qsur=-(Qai(i,j)-Qi2(i,j))/(hfus1(i,j)*MelRho)
              ELSE
                Qsur=-(Qai(i,j)-Qi2(i,j))/(hfus*SnowWetRho(ng))
              END IF
!
!  Compute snow melt/meltpond refreeze rate using the decision tree of
!  Durski et al. (2019).
!
              Fi(i,j,icW_ai)=0.0_r8
              Fi(i,j,icW_ro)=0.0_r8
!
              IF ((Si(i,j,linew,isHsno).gt.eps).and.                    &
     &            (Fi(i,j,icIsst).ge.0.0_r8)) THEN
!
!  (CW) snow>0 -> sfw -> Tis>0 -> Qsurf>0, snow melt (1, 7)
!       Alternatively, Qsurf<0, no snow melt (2, 8)
!       Note that because it is snow melting into meltpond, the density
!       scaling factor should be SnowWetRho/MelRho instead of
!       SnowWetRho/RhoSW.
!
                Si(i,j,linew,isHsno)=Si(i,j,linew,isHsno)-              &
     &                               Si(i,j,linew,isAice)*              &
     &                               MAX(Qsur, 0.0_r8)*dtice(ng)
                Si(i,j,linew,isHmel)=Si(i,j,linew,isHmel)+              &
     &                               Si(i,j,linew,isAice)*              &
     &                               MAX(Qsur, 0.0_r8)*dtice(ng)*       &
     &                               SnowWetRho(ng)/MelRho
!
              ELSE IF ((Si(i,j,linew,isHmel).gt.eps).and.               &
     &                 (Fi(i,j,icIsst).le.tfrz)) THEN
!
!  (CW) snow -> sfw>0 -> Tis<Tfrz -> Qsurf<0, meltpond refreeze (6, 15)
!       Alternatively, Qsurf>0, no meltpond refreeze (5, 14).
!
                Fi(i,j,icW_ai)=MIN(Qsur, 0.0_r8)
                Si(i,j,linew,isHmel)=Si(i,j,linew,isHmel)+              &
     &                               Si(i,j,linew,isAice)*              &
     &                               MIN(Qsur, 0.0_r8)*dtice(ng)
!
              ELSE IF ((Si(i,j,linew,isHsno).le.eps).and.               &
     &                 (Si(i,j,linew,isHmel).ge.eps).and.               &
     &                 (Fi(i,j,icIsst).gt.tfrz)) THEN
!
!  (CW) snow=0 -> sfw>0 -> Tis>Tfrz -> Qsurf>0, ice melt (10, 12)
!       Alternatively, Qsurf<0, no ice melt (11, 13).
!
                Fi(i,j,icW_ai)=MAX(Qsur, 0.0_r8)
                Si(i,j,linew,isHmel)=Si(i,j,linew,isHmel)+              &
     &                               Si(i,j,linew,isAice)*              &
     &                               MAX(Qsur, 0.0_r8)*dtice(ng)
!
              ELSE IF ((Si(i,j,linew,isHsno).lt.eps).and.               &
     &                 (Si(i,j,linew,isHmel).lt.eps).and.               &
     &                 (Fi(i,j,icIsst).gt.tfrz)) THEN
!
!  (CW) snow=0 -> sfw=0 -> Tis>Tfrz -> Qsurf>0, ice melt (16, 18)
!       Alternatively, Qsurf<0, no ice melt (17, 19).
!
                Fi(i,j,icW_ai)=MAX(Qsur, 0.0_r8)
                Si(i,j,linew,isHmel)=Si(i,j,linew,isHmel)+              &
     &                               Si(i,j,linew,isAice)*              &
     &                               MAX(Qsur, 0.0_r8)*dtice(ng)
!
!  (CW) Otherwise no surface phase change, only ice cooling.
!       (3, 4, 9, 20, 21)
!
              END IF
!
!  (CW) Update snow thickness with snowfall rate.
!
              IF (Fi(i,j,icSnof).gt.0.0_r8) THEN
                Si(i,j,linew,isHsno)=Si(i,j,linew,isHsno)+              &
     &                               Si(i,j,linew,isAice)*              &
     &                               Fi(i,j,icSnof)/SnowDryRho(ng)*     &
     &                               dtice(ng)
              END IF
!
!  (CW) Update meltpond thickness with rainfall rate.
!
              IF (rain(i,j).gt.0.0_r8) THEN
                Si(i,j,linew,isHmel)=Si(i,j,linew,isHmel)+              &
     &                               Si(i,j,linew,isAice)*              &
     &                               rain(i,j)/RhoSW*dtice(ng)
              END IF
!
!  Limit the amount of surface water by the smaller of a max limit and
!  ice thickness (SMD).
!
              cff=MIN(max_hmelt(ng), Si(i,j,linew,isHice))
              IF (Si(i,j,linew,isHmel).gt.cff) THEN
                Fi(i,j,icW_ro)=(Si(i,j,linew,isHmel)-cff)/dtice(ng)
                Si(i,j,linew,isHmel)=cff
              END IF
            ELSE
!
!  If no ice, set state variable to default values.
!
              Fi(i,j,icW_ai)=0.0_r8
              Fi(i,j,icW_ro)=0.0_r8
              Si(i,j,linew,isHsno)=0.0_r8
              Si(i,j,linew,isHmel)=0.0_r8
            END IF
!
!-----------------------------------------------------------------------
!  Solve for ice-ocean boundary layer exchange.
!-----------------------------------------------------------------------
!
!  Yaglom and Kader (1974) formulation for turbulent roughness length
!  scales "z0t" and "z0s".
!
!  (CW) Revert algorithm back to MK89 Eq 19.
!
            z0o=MAX(z0ii*ice_thick(i,j), 0.00001_r8)
            z0i=MAX(z0oi*utau(i,j)*utau(i,j), 0.00001_r8)
            z0=EXP(       Si(i,j,linew,isAice) *LOG(z0i) +              &
     &             (1._r8-Si(i,j,linew,isAice))*LOG(z0o))
            z0=MIN(MAX(z0, 0.01_r8), 0.1_r8)
            dztop=z_w(i,j,N(ng))-z_r(i,j,N(ng))
            zdz0=dztop/z0
            zdz0=MAX(zdz0, 3.0_r8)
!
!  (CW) rno = utau*z0/nu (MK89 Eq 18b).
!
            rno=utau(i,j)*z0/nu
            termt=ykf*SQRT(rno)*prt**0.666667_r8
            terms=ykf*SQRT(rno)*prs**0.666667_r8
            cht(i,j)=utau(i,j)/(tpr*(LOG(zdz0)/kappa+termt))
            chs(i,j)=utau(i,j)/(tpr*(LOG(zdz0)/kappa+terms))
!
!  Temperature and salinity of molecular sublayer under ice.
!
            cff=1.0_r8-brnfr(i,j)
            hfus1(i,j)=hfus*cff+                                        &
     &                 Fi(i,j,icT0mk)*Cp_w-                             &
     &                 (cff*Cp_i+brnfr(i,j)*Cp_w)*Si(i,j,linew,isTice)
!
!  If heat flux is large enough to bring surface temperature below
!  freezing point, freeze surface water.
!
            tfrz=frln*salt_top(i,j)
!
!  (CW) See explaination of W_ao term in MK89 Pg 10,939
!
            Fi(i,j,icW_ao)=0.0_r8
            IF (((temp_top(i,j).le.tfrz).and.(Qao(i,j).gt.0.0_r8)).or.  &
     &          ((temp_top(i,j).ge.tfrz).and.(Qao(i,j).lt.0.0_r8).and.  &
     &           (Si(i,j,linew,isAice).gt.0.0_r8))) THEN
              Fi(i,j,icW_ao)=Qao(i,j)/(hfus1(i,j)*RhoSW)
            END IF
!
            IF ((Si(i,j,linew,isAice).gt.min_ai(ng)).and.               &
     &          (Si(i,j,linew,isHice).gt.min_hi(ng))) THEN
!
!  Calculate s0mk, t0mk and total freeze/melt of ice layer
!    Note s0mk and t0mk are used in the next time step, thus the first
!    time step in a simulation uses predefined values only.
!
              Fi(i,j,icW_io)=(Qio(i,j)/RhoSW+                           &
     &                        Cp_w*cht(i,j)*(Fi(i,j,icT0mk)-            &
     &                                       temp_top(i,j)))/hfus1(i,j)
!
!  Compute total mass exchange of ice-covered and ice-free zones (MK89
!  Eq 16b).
!
              xtot=Si(i,j,linew,isAice)*Fi(i,j,icW_io)+                 &
     &             (1.0_r8-Si(i,j,linew,isAice))*Fi(i,j,icW_ao)
!
!  (CW) Revert algorithm back to MK89 Eqs 20, 21.
!
              Fi(i,j,icS0mk) =                                          &
     &  (chs(i,j)*salt_top(i,j) -                                       &
     &   (xtot-Si(i,j,linew,isAice)*Fi(i,j,icW_ro))*sice(i,j)) /        &
     &  (chs(i,j) -                                                     &
     &   (xtot-Si(i,j,linew,isAice)*Fi(i,j,icW_ro)) -                   &
     &   (1._r8-Si(i,j,linew,isAice))*stflx(i,j,isalt)/salt_top(i,j))
!
!  Bound values.
!
              Fi(i,j,icS0mk)=MAX(Fi(i,j,icS0mk), 0.0_r8)
              Fi(i,j,icS0mk)=MIN(Fi(i,j,icS0mk), 40.0_r8)
!
!  MK89 Eq 23b.
!
              Fi(i,j,icT0mk)=frln*Fi(i,j,icS0mk)
            ELSE
!
!  If no ice, set state variable to default values.
!
              Fi(i,j,icW_io)=0.0_r8
              xtot=(1.0_r8-Si(i,j,linew,isAice))*Fi(i,j,icW_ao)
              Fi(i,j,icS0mk)=salt_top(i,j)
              Fi(i,j,icT0mk)=temp_top(i,j)
            END IF
!
!  Compute total ice-ocean mass flux.
!
            Fi(i,j,icIOmf)=xtot-                                        &
     &                     Si(i,j,linew,isAice)*Fi(i,j,icW_ro)+         &
     &                     Fi(i,j,icW_fr)
!
!-----------------------------------------------------------------------
!  Recompute surface heat and salt fluxes.
!-----------------------------------------------------------------------
!
            IF (Si(i,j,linew,isAice).gt.min_ai(ng)) THEN
#ifdef ICE_SHOREFAST
              hh=h(i,j)+Zt_avg1(i,j)
              clear=MAX(hh-0.9_r8*Si(i,j,linew,isHice), 0.0_r8)
              fac_sf=MIN(MAX(clear-0.5_r8, 0.0_r8)/1.0_r8, 1.0_r8)
#else
              fac_sf=1.0_r8
#endif
!
!  Adjust surface heat flux (MK89 Eq 16).
!
              stflx(i,j,itemp)=                                         &
     &          (1.0_r8-Si(i,j,linew,isAice))*Qao(i,j)+                 &
     &          (Si(i,j,linew,isAice)*Qio(i,j)-xtot*hfus1(i,j)*RhoSW)*  &
     &          fac_sf
#ifdef ICE_ALBEDO
!
!  Add short-wave radiation flux that penetrats sea ice.
!
              stflx(i,j,itemp)=stflx(i,j,itemp)-                        &
     &          Si(i,j,linew,isAice)*Fi(i,j,icSicb)*RhoCp
#endif
!
!  Change surface heat flux back to ROMS convention (Celsius m/s).
!
              stflx(i,j,itemp)=-stflx(i,j,itemp)*RhoCpr
!
!  Adjust surface freshwater flux (MK89 Eq 20).
!
              stflx(i,j,isalt)=                                         &
     &          (1.0_r8-Si(i,j,linew,isAice))*stflx(i,j,isalt)-         &
     &          (xtot-Si(i,j,linew,isAice)*Fi(i,j,icW_ro))*             &
     &          (sice(i,j)-salt_top(i,j))*fac_sf
            END IF
          ELSE
!
!  Under ice shelves, surface fluxes are not computed, but frazil ice
!  formation is still allowed.
!
            Fi(i,j,icIOmf)=Fi(i,j,icW_fr)
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute tracer fluxes in sea ice.
!-----------------------------------------------------------------------
!
      IF (Licetrc(ng)) THEN
        ic=0
        DO itrc=NAT+1,NT(ng)
          IF (LtracerIce(itrc,ng)) THEN
            ic=ic+1
            DO j=Jstr,Jend
              DO i=Istr,Iend
                IF (overOcean(i,j)) THEN
                  IF (Si(i,j,linew,isAice).gt.min_ai(ng) .and.          &
     &                Si(i,j,linew,isHice).gt.min_hi(ng)) THEN
                    cff1=Fi(i,j,icIOmf)-Fi(i,j,icW_fr)
                    IF (cff1.lt.0.0_r8) THEN
                      cff2=Ti(i,j,linew,ic)/(Si(i,j,linew,isHice)+eps)
                    ELSE
                      cff2=t(i,j,N(ng),nrhs,itrc)
                    ENDIF
                    stflx(i,j,itrc)=stflx(i,j,itrc)-                    &
     &                              cff1*(cff2-t(i,j,N(ng),nrhs,itrc))
                    Ti(i,j,linew,ic)=Ti(i,j,linew,ic)+                  &
     &                               cff1*cff2*dtice(ng)
                    Ti(i,j,linew,ic)=MAX(Ti(i,j,linew,ic), 0.0_r8)
                  ELSE
                    Ti(i,j,linew,ic)=0.0_r8
                  END IF
                END IF
              END DO
            END DO
          END IF
        END DO
      END IF
!
!  Add external tracer from ice melt.
!
      IF (Licesrc(ng)) THEN
        DO itrc=NAT+1,NT(ng)
          IF (LsourceIce(itrc,ng)) THEN
            DO j=Jstr,Jend
              DO i=Istr,Iend
                IF (overOcean(i,j)) THEN
                  cff=Si(i,j,linew,isAice)*Fi(i,j,icW_ro)+              &
     &                Si(i,j,linew,isAice)*                             &
     &                MAX(-Fi(i,j,icW_io), 0.0_r8)+                     &
     &                (1.0_r8-Si(i,j,linew,isAice))*                    &
     &                MAX(-Fi(i,j,icW_ao), 0.0_r8)
                  stflx(i,j,itrc)=stflx(i,j,itrc)+cff*trcIce(itrc,ng)
                END IF
              END DO
            END DO
          END IF
        END DO
      END IF
!
!  Add external tracer from ice freeze.
!
      IF (Lfrzsrc(ng)) THEN
        DO itrc=NAT+1,NT(ng)
          IF (LsourceFrz(itrc,ng)) THEN
            DO j=Jstr,Jend
              DO i=Istr,Iend
                IF (overOcean(i,j)) THEN
                  cff=Si(i,j,linew,isAice)*                             &
     &                MAX(Fi(i,j,icW_io), 0.0_r8)+                      &
     &                (1.0_r8-Si(i,j,linew,isAice))*                    &
     &                MAX(Fi(i,j,icW_ao), 0.0_r8)
                  stflx(i,j,itrc)=stflx(i,j,itrc)+cff*trcFrz(itrc,ng)
                END IF
              END DO
            END DO
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Update ice properties.
!-----------------------------------------------------------------------
!
!  Track the amount of new ice produced thermodynamically to calculate
!  average ice age.
!
!  (CW) Revert algorithm back to MK89.
!
      cff=dtice(ng)*RhoSw/IceRho(ng)
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (overOcean(i,j)) THEN
!
!  Store old ice mass and thickness.
!
            hi_old=Si(i,j,linew,isHice)
            ai_old=Si(i,j,linew,isAice)
!
!  MK89 Eq 3.
!
            Si(i,j,linew,isHice)=Si(i,j,linew,isHice)+                  &
     &        (Si(i,j,linew,isAice)*(Fi(i,j,icW_io)-Fi(i,j,icW_ai))+    &
     &         (1.0_r8-Si(i,j,linew,isAice))*Fi(i,j,icW_ao)+            &
     &         Fi(i,j,icW_fr))*cff
!
!  MK89 Eq 4, implicit formulation.
!
            IF (Fi(i,j,icW_ao).lt.0.0_r8) THEN
              phi=0.5_r8
            ELSE
              phi=4.0_r8
            END IF
            cff1=cff*(phi*Fi(i,j,icW_ao)+Fi(i,j,icW_fr))/ice_thick(i,j)
            Si(i,j,linew,isAice)=                                      &
     &        (cff1+Si(i,j,linew,isAice))/(1.0_r8+cff1)
            Si(i,j,linew,isAice)=MIN(Si(i,j,linew,isAice), max_ai(ng))
!
!  If sea ice concentration decreases, reduce snow and meltpond masses
!  accordingly; modify heat and salt flux to account the change in snow
!  and meltpond mass.
!
            IF (Si(i,j,linew,isAice).lt.ai_old) THEN
              cff1=Si(i,j,linew,isHsno)*SnowWetRho(ng)*                 &
     &             (1.0_r8-Si(i,j,linew,isAice)/MAX(ai_old, eps))
              cff2=Si(i,j,linew,isHmel)*MelRho*                         &
     &             (1.0_r8-Si(i,j,linew,isAice)/MAX(ai_old, eps))
!
              Fi(i,j,icW_ro)=Fi(i,j,icW_ro)+                            &
     &                       (cff1+cff2)/(RhoSW*dtice(ng)*              &
     &                                    (Si(i,j,linew,isAice)+eps))
              Fi(i,j,icIOmf)=Fi(i,j,icIOmf)-                            &
     &                       (cff1+cff2)/(RhoSW*dtice(ng)*              &
     &                                    (Si(i,j,linew,isAice)+eps))
!
              stflx(i,j,itemp)=stflx(i,j,itemp)-                        &
     &                         hfus*cff1/dtice(ng)*RhoCpr
              stflx(i,j,isalt)=stflx(i,j,isalt)-                        &
     &                         (cff1+cff2)/(RhoSW*dtice(ng))*           &
     &                         (salt_top(i,j)-sice(i,j))
!
              Si(i,j,linew,isHsno)=Si(i,j,linew,isHsno)*                &
     &                             Si(i,j,linew,isAice)/MAX(ai_old, eps)
              Si(i,j,linew,isHmel)=Si(i,j,linew,isHmel)*                &
     &                             Si(i,j,linew,isAice)/MAX(ai_old, eps)
            END IF

#ifdef ICE_CONVSNOW
!
!  If snow base is below sea level, then raise the snow base to sea
!  level by converting some snow to ice (N.B. "hstar" is also weighted
!  by "isAice" like "isHsno" and "isHice").
!
            hstar=Si(i,j,linew,isHsno)-                                 &
     &            Si(i,j,linew,isHice)*(RhoSW-IceRho(ng))/SnowDryRho(ng)
            IF (hstar.gt.0.0_r8) THEN
              cff1=hstar/RhoSW
              Si(i,j,linew,isHsno)=Si(i,j,linew,isHsno)-                &
     &                             cff1*IceRho(ng)
              Si(i,j,linew,isHice)=Si(i,j,linew,isHice)+                &
     &                             cff1*SnowDryRho(ng)
            END IF
#endif
!
!  Determine age of the sea ice. Any new ice production reduces the
!  overall age of the ice parcel.
!
            IF ((Si(i,j,linew,isIage).le.0.0_r8).and.                   &
     &          (Si(i,j,linew,isHice).gt.min_hi(ng))) THEN
!
!  New ice.
!
              Si(i,j,linew,isIage)=dtice(ng)*sec2day
            ELSE IF((Si(i,j,linew,isIage).gt.0.0_r8).and.               &
     &              (Si(i,j,linew,isHice).gt.min_hi(ng))) THEN
!
!  Ice age increases.
!
              delta_mi=MIN(MAX(Si(i,j,linew,isHice)-hi_old, 0.0_r8)/    &
     &                     Si(i,j,linew,isHice), 1.0_r8)
              Si(i,j,linew,isIage)=Si(i,j,linew,isIage)+                &
     &                             dtice(ng)*sec2day-                   &
     &                             Si(i,j,linew,isIage)*delta_mi
            ELSE
!
!  Ice melted.
!
              Si(i,j,linew,isIage)=0.0_r8
            ENDIF
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Bound values.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Si(i,j,linew,isAice)=MIN(Si(i,j,linew,isAice), max_ai(ng))
          Si(i,j,linew,isAice)=MAX(Si(i,j,linew,isAice), 0.0_r8)
          Si(i,j,linew,isHice)=MAX(Si(i,j,linew,isHice), 0.0_r8)
          Si(i,j,linew,isHsno)=MAX(Si(i,j,linew,isHsno), 0.0_r8)
          Si(i,j,linew,isHmel)=MAX(Si(i,j,linew,isHmel), 0.0_r8)
          Si(i,j,linew,isTice)=MAX(Si(i,j,linew,isTice), -70.0_r8)
!
          IF (Si(i,j,linew,isHice).le.0.0_r8)                           &
     &      Si(i,j,linew,isAice)=0.0_r8
          IF (Si(i,j,linew,isAice).le.0.0_r8)                           &
     &      Si(i,j,linew,isHice)=0.0_r8
        END DO
      END DO
      ic=0
!
      IF (Licetrc(ng)) THEN
        DO itrc=NAT+1,NT(ng)
          IF (LtracerIce(itrc,ng)) THEN
            ic=ic+1
            DO j=Jstr,Jend
              DO i=Istr,Iend
                IF (Si(i,j,linew,isAice).le.min_ai(ng) .or.             &
     &              Si(i,j,linew,isHice).le.min_hi(ng)) THEN
                  Ti(i,j,linew,ic)=0.0_r8
                END IF
              END DO
            END DO
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions.
!-----------------------------------------------------------------------
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fi(:,:,icIsst))
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fi(:,:,icQcon))
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fi(:,:,icQrhs))
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  stflx(:,:,isalt))
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  stflx(:,:,itemp))
!
      CALL ice_bc2d_tile (ng, tile, model, isAice,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isAice),                             &
     &                    LBC(:,ibICE(isAice),ng))
!
      CALL ice_bc2d_tile (ng, tile, model, isHice,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isHice),                             &
     &                    LBC(:,ibICE(isHice),ng))
!
      CALL ice_bc2d_tile (ng, tile, model, isHsno,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isHsno),                             &
     &                    LBC(:,ibICE(isHsno),ng))
!
      CALL ice_bc2d_tile (ng, tile, model, isHmel,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isHmel),                             &
     &                    LBC(:,ibICE(isHmel),ng))
!
      CALL ice_bc2d_tile (ng, tile, model, isIage,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isIage),                             &
     &                    LBC(:,ibICE(isIage),ng))
!
      CALL ice_tibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isHice),                             &
     &                    Si(:,:,:,isTice),                             &
     &                    Si(:,:,:,isEnth))
!
      IF (Licetrc(ng)) THEN
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, NTICE(ng),             &
     &                    Ti(:,:,linew,:))
      END IF
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isAice))
!
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isHice))
!
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isHsno))
!
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isHmel))
!
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isIage))
!
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isEnth))
!
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isTice))
!
        IF (Licetrc(ng)) THEN
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, NTICE(ng),     &
     &                            Ti(:,:,linew,:))
        END IF
      END IF
#ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Si(:,:,linew,isAice),                         &
     &                    Si(:,:,linew,isHice),                         &
     &                    Si(:,:,linew,isHsno),                         &
     &                    Si(:,:,linew,isHmel))
!
      CALL mp_exchange2d (ng, tile, model, 3,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Si(:,:,linew,isIage),                         &
     &                    Si(:,:,linew,isEnth),                         &
     &                    Si(:,:,linew,isTice))
!
      IF (Licetrc(ng)) THEN
        CALL mp_exchange3d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj, 1, NTICE(ng),           &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      Ti(:,:,linew,:))
      END IF
#endif
!
      RETURN
      END SUBROUTINE ice_thermo_tile
      END MODULE ice_thermo_mod
