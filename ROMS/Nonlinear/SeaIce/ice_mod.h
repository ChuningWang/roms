      MODULE mod_ice
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2024 The ROMS/TOMS Group      Paul Budgell       !
!    Licensed under a MIT/X style license           Katherine Hedstrom !
!    See License_ROMS.md                            Scott M. Durski    !
!=======================================================================
!                                                                      !
!  The sea-ice model was originally written and tested by Paul Budgell !
!  (2005). It has been updated and improved by Kate Hedstrom and Scott !
!  Durski.                                                             !
!                                                                      !
!  It uses the elastic-viscous-plastic (EVP) rheology of Hunke and     !
!  Dukowicz (1997) and Hunke (2001). It has a simple one-layer ice and !
!  snow thermodynamics with a molecular sublayer under the ice (Mellor !
!  and Kantha, 1989).                                                  !
!                                                                      !
!  Variables:                                                          !
!                                                                      !
!  The ice model prognostic 'state' variable are declares in compact   !
!  form in the derived-type structure as:                              !
!                                                                      !
!    ICE(ng) % state(:,:,TimeLevel,StateIndex)                         !
!                                                                      !
!  where the 'StateIndex' are as follows:                              !
!                                                                      !
!    isAice    ice concentration (grid cell fraction)                  !
!    isEnth    enthalpy of the ice/brine system, ice heat content      !
!    isHage    thickness associated with age of ice (m)                !
!    isHice    average ice thickness (m), ice mass divided by area     !
!    isHmel    surface meltwater thickness on ice (m)                  !
!    isHsno    average thickness of snow coverage (m), mass snow       !
!    isIage    age of ice(s)                                           !
!    isISxx    internal ice stress, xx-component (N/m)                 !
!    isISxy    internal ice stress, xy-component (N/m)                 !
!    isISyy    internal ice stress, yy-component (N/m)                 !
!    isTice    ice interior temperature (Celsius)                      !
!    isUice    ice U-velocity (m/s)                                    !
!    isUevp    elastic-viscous-plastic ice U-velocity (m/s)            !
!    isVice    ice V-velocity (m/s)                                    !
!    isVevp    elastic-viscous-plastic ice V-velocity (m/s)            !
#if defined ICE_BIO
!    isIphy    ice biology, algae                                      !
!    isINO3    ice biology, nitrate                                    !
!    isINH4    ice biology, ammonia                                    !
!    isIlog    ice biology log                                         !
#endif
!                                                                      !
!  The ice model internal 'field' arrays are declares in compact form  !
!  in the derived-type structure as:                                   !
!                                                                      !
!    ICE(ng) % field(:,:,FieldIndex)                                   !
!                                                                      !
!  where the 'FieldIndex' are as follows:                              !
!                                                                      !
!    icAIus    surface Air-Ice U-stress (N/m2)                         !
!    icAIvs    surface Air-Ice V-stress (N/m2)                         !
!    icBvis    ice bulk viscosity                                      !
!    icHsse    sea surface elevation (m)                               !
!    icPgrd    gridded ice strength parameter (unitless)               !
!    icPice    ice pressure or strength (N/m2)                         !
!    icQcon    gradient heat conductivity over ice/snow (W/m2/K)       !
!    icQrhs    RHS surface net heat flux over ice/snow (W/m2)          !
!    icSvis    ice shear viscosity                                     !
!    icS0mk    salinity of molecular sublayer under ice (unitless)     !
!    icT0mk    temperature of molecular sublayer under ice (Celsius)   !
!    icIsst    temperature at the snow/atmosphere interface (Celsius)  !
!    icUavg    vertically averaged mixed-layer U-velocity (m/s)        !
!    icVavg    vertically averaged mixed-layer V-velocity (m/s)        !
!    icWdiv    rate of ice divergence (m3/s)                           !
!    icW_ai    rate of melt/freeze at Air/Ice edge (m3/s)              !
!    icW_ao    rate of melt/freeze at Air/Ocean edge (m3/s)            !
!    icW_fr    rate of ice accretion due to frazil growth (m3/s)       !
!    icW_io    rate of melt/freeze at Ice/Ocean edge (m3/s)            !
!    icW_ro    rate of melt/freeze runoff into ocean (m3/s)            !
#if defined ICE_BULK_FLUXES
!    icAlbi    sea ice albedo                                          !
!    icAlbo    sea surface albedo                                      !
!    icFrsw    short-wave penetration fration in ice                   !
!    icSics    net short-wave flux at sea ice surface (W/m2)           !
!    icSicb    net short-wave flux at sea ice base (W/m2)              !
!    icSocn    net short-wave flux at sea surface (W/m2)               !
!    icRfac    lmd_r1 scaling factor                                   !
!    icSnof    snowfall rate (kg/m2/s)                                 !
#endif
!    icIOmf    Ice-Ocean mass flux (m/s)                               !
!    icIOvs    Ice-Ocean velocity shear magnitude (m/s)                !
!    icIOfv    Ice-Ocean friction velocity (m/s)                       !
!    icIOmt    Ice-Ocean momentum transfer coefficient (m/s)           !
!                                                                      !
!  The ice model internal 'tracer' arrays:                             !
!                                                                      !
!    Ti        tracer concentration in sea ice (m TU)                  !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
      PUBLIC :: allocate_ice
      PUBLIC :: deallocate_ice
      PUBLIC :: initialize_ice
!
!-----------------------------------------------------------------------
!  Ice model identification indices.
!-----------------------------------------------------------------------
!
      integer :: idAice       ! ice concentration
      integer :: idAiCL       ! ice concentration climatology
      integer :: idHage       ! thickness associated with age of ice
      integer :: idHice       ! ice thickness
      integer :: idHiCL       ! ice thickness climatology
      integer :: idHmel       ! surface meltwater thickness
      integer :: idHsno       ! snow cover thickness
      integer :: idIage       ! ice age
      integer :: idIOfv       ! ice-ocean friction velocity
      integer :: idIOmf       ! ice-ocean mass flux
      integer :: idIOmt       ! ice-ocean momentum transfer coefficient
      integer :: idIsst       ! ice/snow surface temperature
      integer :: idISxx       ! internal ice stress xx-component
      integer :: idISxy       ! internal ice stress xy-component
      integer :: idISyy       ! internal ice stress yy-component
      integer :: idIphy       ! ice biology, algae
      integer :: idINO3       ! ice biology, nitrate
      integer :: idINH4       ! ice biology, ammonia
      integer :: idIlog       ! ice biology log
      integer :: idS0mk       ! under ice molecular sublayer salinity
      integer :: idT0mk       ! under ice molecular sublayer temperature
      integer :: idTice       ! ice interior temperature
      integer :: idUice       ! ice U-velocity
      integer :: idUiCL       ! ice U-velocity climatology
      integer :: idUiER       ! ice U-eastward  at RHO-points
      integer :: idVice       ! ice V-velocity
      integer :: idViCL       ! ice V-velocity climatology
      integer :: idViNR       ! ice V-northward at RHO-points
      integer :: idWdiv       ! rate of ice divergence
      integer :: idW_ai       ! rate of melt/freeze at air/ice edge
      integer :: idW_ao       ! rate of melt/freeze at air/ocean edge
      integer :: idW_fr       ! rate of ice accretion by frazil growth
      integer :: idW_io       ! rate of melt/freeze at ice/ocean edge
      integer :: idW_ro       ! rate of melt/freeze runoff into ocean
      integer :: idAlbi       ! sea ice albedo
      integer :: idAlbo       ! ocean sufface albedo
      integer :: idFrsw       ! short-wave penetration fraction
      integer :: idSics       ! net short-wave flux at sea ice surface
      integer :: idSicb       ! net short-wave flux at sea ice base
      integer :: idSocn       ! net short-wave flux at sea surface
      integer :: idRfac       ! lmd_r1 scaling factor
      integer :: idSnof       ! snowfall rate
!
      integer, allocatable :: idStrc(:)  ! Tracer variables
!
!  Ice model state prognostic variables indices.
!
#if defined ICE_BIO
      integer, parameter :: nIceS = 19  ! number of ice state variables
#else
      integer, parameter :: nIceS = 15  ! number of ice state variables
#endif
!
      integer :: iSice(nIceS)           ! state I/O indices
!
      integer, parameter :: isAice =  1 ! ice concentration
      integer, parameter :: isHice =  2 ! ice thickness
      integer, parameter :: isHmel =  3 ! meltwater thickness on ice
      integer, parameter :: isHsno =  4 ! snow thickness
      integer, parameter :: isIage =  5 ! ice age
      integer, parameter :: isISxx =  6 ! internal ice xx-stress
      integer, parameter :: isISxy =  7 ! internal ice xy-stress
      integer, parameter :: isISyy =  8 ! internal ice yy-stress
      integer, parameter :: isTice =  9 ! ice interior temperature
      integer, parameter :: isUice = 10 ! ice U-velocity
      integer, parameter :: isVice = 11 ! ice V-velocity
      integer, parameter :: isEnth = 12 ! ice/brine enthalpy
      integer, parameter :: isHage = 13 ! thickness linked with ice age
      integer, parameter :: isUevp = 14 ! EVP ice U-velocity
      integer, parameter :: isVevp = 15 ! EVP ice V-velocity
#if defined ICE_BIO
      integer, parameter :: isIphy = 16 ! ice biology, algae
      integer, parameter :: isINO3 = 17 ! ice biology, nitrate
      integer, parameter :: isINH4 = 18 ! ice biology, ammonia
      integer, parameter :: isIlog = 19 ! ice biology log
#endif
!
      integer :: ice_ic_brylast  ! Boundary layer salinity
!
!  Ice model state lateral boundary conditions indices.
!
      integer :: ibICE(nIceS)           ! indices to LBC switch
      integer :: iceOBC(4,nIceS)        ! I/O metadata indices
!
!  Ice model internal variables indices.
!
#if defined ICE_BULK_FLUXES
      integer, parameter :: nIceF = 32  ! number of ice field variables
#else
      integer, parameter :: nIceF = 24  ! number of ice field variables
#endif
!
      integer :: iFice(nIceF)           ! internal fields I/O indices
!
      integer, parameter :: icAIus =  1 ! surface Air-Ice U-stress
      integer, parameter :: icAIvs =  2 ! surface Air-Ice V-stress
      integer, parameter :: icBvis =  3 ! ice bulk viscosity
      integer, parameter :: icHsse =  4 ! sea surface elevation
      integer, parameter :: icIOfv =  5 ! Ice-Ocean friction velocity
      integer, parameter :: icIOmf =  6 ! Ice-Ocean mass flux
      integer, parameter :: icIOmt =  7 ! Ice-Ocean momentum transfer
      integer, parameter :: icIOvs =  8 ! Ice-Ocean velocity shear
      integer, parameter :: icIsst =  9 ! ice/snow surface temperature
      integer, parameter :: icPgrd = 10 ! gridded ice strength
      integer, parameter :: icPice = 11 ! ice pressure or strength
      integer, parameter :: icQcon = 12 ! ice/snow heat conductivity
      integer, parameter :: icQrhs = 13 ! RHS heat flux over ice/snow
      integer, parameter :: icSvis = 14 ! ice shear viscosity
      integer, parameter :: icS0mk = 15 ! molecular sublayer salinity
      integer, parameter :: icT0mk = 16 ! molecular sublayer temperature
      integer, parameter :: icUavg = 17 ! average mixed-layer U-velocity
      integer, parameter :: icVavg = 18 ! average mixed-layer V-velocity
      integer, parameter :: icWdiv = 19 ! ice divergence rate
      integer, parameter :: icW_ai = 20 ! melt/freeze rate at Air/Ice
      integer, parameter :: icW_ao = 21 ! melt/freeze rate at Air/Ocean
      integer, parameter :: icW_fr = 22 ! ice accretion rate by frazil
      integer, parameter :: icW_io = 23 ! melt/freeze rate at Ice/Ocean
      integer, parameter :: icW_ro = 24 ! melt/freeze rate runoff
#if defined ICE_BULK_FLUXES
      integer, parameter :: icAlbi = 25 ! sea ice albedo
      integer, parameter :: icAlbo = 26 ! ocean surface albedo
      integer, parameter :: icFrsw = 27 ! SW penetration fraction in ice
      integer, parameter :: icSics = 28 ! net SW flux at sea ice surface
      integer, parameter :: icSicb = 29 ! net SW flux at sea ice base
      integer, parameter :: icSocn = 30 ! net SW flux at sea surface
      integer, parameter :: icRfac = 31 ! lmd_r1 scaling factor
      integer, parameter :: icSnof = 32 ! snowfall rate
#endif
!
!-----------------------------------------------------------------------
!  Ice model parameters.
!-----------------------------------------------------------------------

#ifdef AVERAGES
!
!  Switches to process time-averaged ice model state and internal
!  variables.
!
      logical, allocatable :: LiceFavg(:,:)     ! internal variables
      logical, allocatable :: LiceSavg(:,:)     ! state variables
      logical, allocatable :: LiceTavg(:,:)     ! tracer
#endif
!
!  Counter and number of Elastic-Viscous-Plastic (EVP) rheology
!  equations timesteps to resolve elastic dynamics.
!
      integer, allocatable :: iEVP(:)           ! timestep counter
      integer, allocatable :: nEVP(:)           ! number of timesteps
!
!  Ice equations time step (s).
!
      integer, allocatable :: dtice(:)          ! viscous
      integer, allocatable :: dtevp(:)          ! elastic
!
!  Density parameters (kg/m3).
!
      real(r8), allocatable :: AirRho(:)        ! air density
      real(r8), allocatable :: IceRho(:)        ! ice density
      real(r8), allocatable :: SnowDryRho(:)    ! dry snow density
      real(r8), allocatable :: SnowWetRho(:)    ! wet snow density
!
!  Boundary layer bulk drag coefficients (nondimensional).
!
      real(r8), allocatable :: Cd_ai(:)         ! air-ice boundary
      real(r8), allocatable :: Cd_io(:)         ! ice-ocean boundary
!
!  Ice strength exponential weighting coefficient on concentration
!  (nondimensional).
!
      real(r8), allocatable :: Astrength(:)
      real(r8), allocatable :: Pstar(:)
!
!  Minimum and maximum average fractional ice concetration coverage
!  (nomdimensional) limiters.
!
      real(r8), allocatable :: min_ai(:)        ! minimum limiter
      real(r8), allocatable :: max_ai(:)        ! maximun limiter
!
!  Minimum average ice thickness (m) limiter.
!
      real(r8), allocatable :: min_hi(:)        ! minimum limiter
!
!  Maximum surface melt water thickness (m) limiter.
!
      real(r8), allocatable :: max_hmelt(:)
!
!  Minimum and maximum ice shear strength (N/m2) limiters.
!
      real(r8), allocatable :: zetaMin(:)       ! minimum limiter
      real(r8), allocatable :: zetaMax(:)       ! maximum limiter
!
!  Turning angle for ice-water drag (radians).
!
      real(r8), allocatable :: stressAng(:)
!
!  Squared ellipticity/eccentricity of the yield curve (nondimensional).
!
      real(r8), allocatable :: ellip_sq(:)
!
!  Model formulation constants.
!
      real(r8) :: ice_emiss       ! ice emissivity (unitless)
      real(r8) :: spec_heat_air   ! specific heat of air (J/kg/K)
      real(r8) :: trans_coeff     ! heat transfer coefficient (unitless)
      real(r8) :: sublimation     ! latent heat of sublimation (J/kg)
!
!  Logical switches
!
      integer, allocatable :: NTICE(:)
      logical, allocatable :: Licetrc(:)
      logical, allocatable :: LtracerIce(:,:)   ! tracer exchange switch
      logical, allocatable :: Licesrc(:)
      logical, allocatable :: LsourceIce(:,:)   ! tracer source switch
      logical, allocatable :: Lfrzsrc(:)
      logical, allocatable :: LsourceFrz(:,:)   ! tracer source switch
!
!  Ice melt water tracer concentration
!
      real(r8), allocatable :: trcIce(:,:)      ! tracer concentration
      real(r8), allocatable :: trcFrz(:,:)      ! tracer concentration
!
!-----------------------------------------------------------------------
!  Define derived-type structure ice model state and internal arrays.
!-----------------------------------------------------------------------
!
      TYPE T_ICE

        real(r8), pointer :: Si(:,:,:,:)             ! [i,j,1:2,1:nIceS]
        real(r8), pointer :: Fi(:,:,:)               ! [i,j,1:nIceF]
        real(r8), pointer :: Ti(:,:,:,:)             ! [i,j,1:2,1:NTICE]

      END TYPE T_ICE
!
      TYPE (T_ICE), allocatable :: ICE(:)            ! [Ngrids]
!
!-----------------------------------------------------------------------
!  Define derived-type structure ice model lateral boundary variables.
!-----------------------------------------------------------------------
!
      TYPE T_ICE_LOBC

        real(r8), pointer :: ice_west (:)
        real(r8), pointer :: ice_east (:)
        real(r8), pointer :: ice_south(:)
        real(r8), pointer :: ice_north(:)
!
        real(r8), pointer :: iceG_west (:,:)
        real(r8), pointer :: iceG_east (:,:)
        real(r8), pointer :: iceG_south(:,:)
        real(r8), pointer :: iceG_north(:,:)

      END TYPE T_ICE_LOBC
!
      TYPE (T_ICE_LOBC), allocatable :: ICE_LOBC(:,:)   ! [nIceS,Ngrids]

#ifdef AVERAGES
!
!-----------------------------------------------------------------------
!  Define derived-type structure ice model state and internal arrays
!  time-averaged variables. Notice that only the requested arrays are
!  allocated and processed.
!-----------------------------------------------------------------------
!
      TYPE T_ICE_AVG

        real(r8), pointer :: var(:,:)                ! [i,j]

      END TYPE T_ICE_AVG
!
      TYPE (T_ICE_AVG), allocatable :: ICE_FAVG(:,:) ! [nIceF,Ngrids]
      TYPE (T_ICE_AVG), allocatable :: ICE_SAVG(:,:) ! [nIceS,Ngrids]
      TYPE (T_ICE_AVG), allocatable :: ICE_TAVG(:,:) ! [MT,Ngrids]
#endif
!
      CONTAINS
!
      SUBROUTINE allocate_ice (ng, LBi, UBi, LBj, UBj, ice_kernel)
!
!=======================================================================
!                                                                      !
!  This routine allocates either the derived-type ice model variables  !
!  (ice_kernel=TRUE) or module parameters (ice_kernel=FALSE).          !
!                                                                      !
!=======================================================================
!
      USE mod_param,   ONLY : Dmem, LBC, Ngrids, NT, MT, NAT
      USE mod_scalars, ONLY : iwest, ieast, isouth, inorth
!
!  Imported variable declarations.
!
      logical, intent(in) :: ice_kernel
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: i, ic, itrc
!
      real(r8), parameter :: IniVal = 0.0_r8
!
      real(r8) :: size2d, Xsize, Ysize
!
      size2d=REAL((UBi-LBi+1)*(UBj-LBj+1),r8)
      Xsize =REAL(UBi-LBi,r8)
      Ysize =REAL(UBj-LBj,r8)

!
!-----------------------------------------------------------------------
!  Allocate ice model parameters.
!-----------------------------------------------------------------------
!
      IF (.not.ice_kernel) THEN

        IF (.not.allocated(idStrc))                                     &
     &    allocate ( idStrc(MT) )
        Dmem(1)=Dmem(1)+REAL(MT,r8)
        idStrc(1:MT)=0

#ifdef AVERAGES

        IF (.not.allocated(LiceFavg))                                   &
          allocate ( LiceFavg(nIceF,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(nIceF*Ngrids,r8)
        LiceFavg(1:nIceF,1:Ngrids)=.FALSE.

        IF (.not.allocated(LiceSavg))                                   &
          allocate ( LiceSavg(nIceS,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(nIceS*Ngrids,r8)
        LiceSavg(1:nIceS,1:Ngrids)=.FALSE.

        IF (.not.allocated(LiceTavg))                                   &
          allocate ( LiceTavg(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
        LiceTavg(1:MT,1:Ngrids)=.FALSE.
#endif
        IF (.not.allocated(iEVP))                                       &
     &    allocate ( iEVP(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        iEVP(1:Ngrids)=0

        IF (.not.allocated(nEVP))                                       &
     &    allocate ( nEVP(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        nEVP(1:Ngrids)=0

        IF (.not.allocated(dtice))                                      &
     &    allocate ( dtice(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        dtice(1:Ngrids)=IniVal

        IF (.not.allocated(dtevp))                                      &
     &    allocate ( dtevp(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        dtevp(1:Ngrids)=IniVal

        IF (.not.allocated(AirRho))                                     &
     &    allocate ( AirRho(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        AirRho(1:Ngrids)=IniVal

        IF (.not.allocated(IceRho))                                     &
     &    allocate ( IceRho(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        IceRho(1:Ngrids)=IniVal

        IF (.not.allocated(SnowDryRho))                                &
     &    allocate ( SnowDryRho(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        SnowDryRho(1:Ngrids)=IniVal

        IF (.not.allocated(SnowWetRho))                                &
     &    allocate ( SnowWetRho(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        SnowWetRho(1:Ngrids)=IniVal

        IF (.not.allocated(Cd_ai))                                      &
     &    allocate ( Cd_ai(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        Cd_ai(1:Ngrids)=IniVal

        IF (.not.allocated(Cd_io))                                      &
     &    allocate ( Cd_io(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        Cd_io(1:Ngrids)=IniVal

        IF (.not.allocated(Astrength))                                  &
     &    allocate ( Astrength(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        Astrength(1:Ngrids)=IniVal

        IF (.not.allocated(Pstar))                                      &
     &    allocate ( Pstar(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        Pstar(1:Ngrids)=IniVal

        IF (.not.allocated(min_ai))                                     &
     &    allocate ( min_ai(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        min_ai(1:Ngrids)=IniVal

        IF (.not.allocated(max_ai))                                     &
     &    allocate ( max_ai(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        max_ai(1:Ngrids)=IniVal

        IF (.not.allocated(min_hi))                                     &
     &    allocate ( min_hi(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        min_hi(1:Ngrids)=IniVal

        IF (.not.allocated(max_hmelt))                                  &
     &    allocate ( max_hmelt(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        max_hmelt(1:Ngrids)=IniVal

        IF (.not.allocated(zetaMin))                                    &
     &    allocate ( zetaMin(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        zetaMin(1:Ngrids)=IniVal

        IF (.not.allocated(zetaMax))                                    &
     &    allocate ( zetaMax(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        zetaMax(1:Ngrids)=IniVal

        IF (.not.allocated(stressAng))                                  &
     &    allocate ( stressAng(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        stressAng(1:Ngrids)=IniVal

        IF (.not.allocated(ellip_sq))                                   &
     &    allocate ( ellip_sq(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        ellip_sq(1:Ngrids)=IniVal

        IF (.not.allocated(NTICE))                                      &
     &    allocate ( NTICE(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        NTICE(1:Ngrids)=0

        IF (.not.allocated(Licetrc))                                    &
     &    allocate ( Licetrc(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        Licetrc(1:Ngrids)=.FALSE.

        IF (.not.allocated(LtracerIce))                                 &
     &    allocate ( LtracerIce(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
        LtracerIce(1:MT,1:Ngrids)=.FALSE.

        IF (.not.allocated(Licesrc))                                    &
     &    allocate ( Licesrc(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        Licesrc(1:Ngrids)=.FALSE.

        IF (.not.allocated(LsourceIce))                                 &
     &    allocate ( LsourceIce(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
        LsourceIce(1:MT,1:Ngrids)=.FALSE.

        IF (.not.allocated(Lfrzsrc))                                    &
     &    allocate ( Lfrzsrc(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
        Lfrzsrc(1:Ngrids)=.FALSE.

        IF (.not.allocated(LsourceFrz))                                 &
     &    allocate ( LsourceFrz(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
        LsourceFrz(1:MT,1:Ngrids)=.FALSE.

        IF (.not.allocated(trcIce))                                     &
     &    allocate ( trcIce(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
        trcIce(1:MT,1:Ngrids)=IniVal

        IF (.not.allocated(trcFrz))                                     &
     &    allocate ( trcFrz(MT,Ngrids) )
        Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
        trcFrz(1:MT,1:Ngrids)=IniVal

      ELSE
!
!-----------------------------------------------------------------------
!  Allocate derived-type structure ice model kernel variables.
!-----------------------------------------------------------------------
!
        IF (ng.eq.1) allocate ( ICE(Ngrids) )
!
!  Nonlinear ice model state variables (Si) and internal kernel field
!  (Fi) arrays.
!
        allocate ( ICE(ng) % Si(LBi:UBi,LBj:UBj,2,nIceS) )
        Dmem(ng)=Dmem(ng)+2.0_r8*REAL(nIceS,r8)*size2d
!
        allocate ( ICE(ng) % Fi(LBi:UBi,LBj:UBj,nIceF) )
        Dmem(ng)=Dmem(ng)+REAL(nIceF,r8)*size2d
!
        IF (Licetrc(ng)) THEN
          allocate ( ICE(ng) % Ti(LBi:UBi,LBj:UBj,2,NTICE(ng)) )
          Dmem(ng)=Dmem(ng)+2.0_r8*REAL(NTICE(ng),r8)*size2d
        END IF
!
!-----------------------------------------------------------------------
!  Allocate derived-type structure ice model lateral boundary variables.
!-----------------------------------------------------------------------
!
        IF (ng.eq.1) allocate ( ICE_LOBC(nIceS,Ngrids) )
!
!  Nonlinear ice model 'state' lateral boundary arrays.
!
        DO i=1,nIceS
          IF (LBC(iwest,ibICE(i),ng)%acquire) THEN
            allocate ( ICE_LOBC(i,ng) % ice_west(LBj:UBj) )
            Dmem(ng)=Dmem(ng)+Ysize

            allocate ( ICE_LOBC(i,ng) % iceG_west(LBj:UBj,2) )
            Dmem(ng)=Dmem(ng)+2.0_r8*Ysize
          END IF
!
          IF (LBC(ieast,ibICE(i),ng)%acquire) THEN
            allocate ( ICE_LOBC(i,ng) % ice_east(LBj:UBj) )
            Dmem(ng)=Dmem(ng)+Ysize

            allocate ( ICE_LOBC(i,ng) % iceG_east(LBj:UBj,2) )
            Dmem(ng)=Dmem(ng)+2.0_r8*Ysize
          END IF
!
          IF (LBC(isouth,ibICE(i),ng)%acquire) THEN
            allocate ( ICE_LOBC(i,ng) % ice_south(LBi:UBi) )
            Dmem(ng)=Dmem(ng)+Xsize

            allocate ( ICE_LOBC(i,ng) % iceG_south(LBi:UBi,2) )
            Dmem(ng)=Dmem(ng)+2.0_r8*Xsize
          END IF
!
          IF (LBC(inorth,ibICE(i),ng)%acquire) THEN
            allocate ( ICE_LOBC(i,ng) % ice_north(LBi:UBi) )
            Dmem(ng)=Dmem(ng)+Xsize

            allocate ( ICE_LOBC(i,ng) % iceG_north(LBi:UBi,2) )
            Dmem(ng)=Dmem(ng)+2.0_r8*Xsize
          END IF
        END DO
#ifdef AVERAGES
!
!-----------------------------------------------------------------------
!  Allocate derived-type structure ice model state and internal arrays
!  time-averaged variables.
!-----------------------------------------------------------------------
!
        IF (ng.eq.1) THEN
          allocate ( ICE_FAVG(nIceF,Ngrids) )
          allocate ( ICE_SAVG(nIceS,Ngrids) )
        END IF
!
!  Time-averaged internal fields. Only fields with metadata (iFice > 0)
!  are allocated and processed.
!
        DO i=1,nIceF
          IF (iFice(i).gt.0) THEN
            IF (LiceFavg(i,ng)) THEN
              allocate ( ICE_FAVG(i,ng) % var(LBi:UBi,LBj:UBj) )
              Dmem(ng)=Dmem(ng)+size2d
            END IF
          END IF
        END DO
!
!  Time-averaged state fields. Only fields with metadata (iSice > 0)
!  are allocated and processed.
!
        DO i=1,nIceS
          IF (iSice(i).gt.0) THEN
            IF (LiceSavg(i,ng)) THEN
              allocate ( ICE_SAVG(i,ng) % var(LBi:UBi,LBj:UBj) )
              Dmem(ng)=Dmem(ng)+size2d
            END IF
          END IF
        END DO
!
!  Time-averaged tracer fields.
!
        IF (ng.eq.1) THEN
          allocate ( ICE_TAVG(MT,Ngrids) )
        END IF
!
        IF (Licetrc(ng)) THEN
          DO itrc=NAT+1,NT(ng)
            IF (LiceTavg(itrc,ng)) THEN
              allocate ( ICE_TAVG(itrc,ng) % var(LBi:UBi,LBj:UBj) )
              Dmem(ng)=Dmem(ng)+size2d
            END IF
          END DO
        END IF
      END IF
#endif
!
      RETURN
      END SUBROUTINE allocate_ice
!
      SUBROUTINE deallocate_ice (ng)
!
!=======================================================================
!                                                                      !
!  This routine deallocates all variables in the module for all nested !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param,   ONLY : Ngrids
#ifdef SUBOBJECT_DEALLOCATION
      USE mod_param,   ONLY : LBC
      USE mod_scalars, ONLY : iwest, ieast, isouth, inorth
!
      USE destroy_mod, ONLY : destroy
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
#ifdef SUBOBJECT_DEALLOCATION
      integer :: i
!
#endif
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", deallocate_ice"

#ifdef SUBOBJECT_DEALLOCATION
!
!-----------------------------------------------------------------------
!  Deallocate each variable in the derived-type T_ICE structure
!  separately.
!-----------------------------------------------------------------------
!
!  Nonlinear ice model state variables 'Si' and internal kernel arrays
!  'Fi'.
!
      IF (.not.destroy(ng, ICE(ng)%Si, MyFile,                          &
     &                 __LINE__, 'ICE(ng)%Si')) RETURN
!
      IF (.not.destroy(ng, ICE(ng)%Fi, MyFile,                          &
     &                 __LINE__, 'ICE(ng)%Fi')) RETURN
!
      IF (.not.destroy(ng, ICE(ng)%Ti, MyFile,                          &
     &                 __LINE__, 'ICE(ng)%Ti')) RETURN
!
!  Nonlinear ice model state lateral boundary arrays.
!
      DO i=1,nIceS
        IF (LBC(iwest,ibICE(i),ng)%acquire) THEN
          IF (.not.destroy(ng, ICE_LOBC(i,ng)%ice_west, MyFile,         &
     &                     __LINE__, 'ICE_LOBC(i,ng)%ice_west')) RETURN
          IF (.not.destroy(ng, ICE_LOBC(i,ng)%iceG_west, MyFile,        &
     &                     __LINE__, 'ICE_LOBC(i,ng)%iceG_west')) RETURN
        END IF
!
        IF (LBC(ieast,ibICE(i),ng)%acquire) THEN
          IF (.not.destroy(ng, ICE_LOBC(i,ng)%ice_east, MyFile,         &
     &                     __LINE__, 'ICE_LOBC(i,ng)%ice_west')) RETURN
          IF (.not.destroy(ng, ICE_LOBC(i,ng)%iceG_east, MyFile,        &
     &                     __LINE__, 'ICE_LOBC(i,ng)%iceG_west')) RETURN
        END IF
!
        IF (LBC(isouth,ibICE(i),ng)%acquire) THEN
          IF (.not.destroy(ng, ICE_LOBC(i,ng)%ice_south, MyFile,        &
     &                     __LINE__, 'ICE_LOBC(i,ng)%ice_south')) RETURN
          IF (.not.destroy(ng, ICE_LOBC(i,ng)%iceG_south, MyFile,       &
     &                     __LINE__,'ICE_LOBC(i,ng)%iceG_south')) RETURN
        END IF
!
        IF (LBC(inorth,ibICE(i),ng)%acquire) THEN
          IF (.not.destroy(ng, ICE_LOBC(i,ng)%ice_north, MyFile,        &
     &                     __LINE__, 'ICE_LOBC(i,ng)%ice_north')) RETURN
          IF (.not.destroy(ng, ICE_LOBC(i,ng)%iceG_south, MyFile,       &
     &                     __LINE__,'ICE_LOBC(i,ng)%iceG_north')) RETURN
        END IF
      END DO

# ifdef AVERAGES
!
!-----------------------------------------------------------------------
!  Time-averaged ice model state and internal variables.
!-----------------------------------------------------------------------
!
      DO i=1,nIceF
        IF (LiceFavg(i,ng)) THEN
          IF (.not.destroy(ng, ICE_FAVG(i,ng)%var, MyFile,              &
     &                     __LINE__, 'ICE_FAVG(i,ng)%var')) RETURN
        END IF
      END DO
!
      DO i=1,nIceS
        IF (LiceSavg(i,ng)) THEN
          IF (.not.destroy(ng, ICE_SAVG(i,ng)%var, MyFile,              &
     &                     __LINE__, 'ICE_SAVG(i,ng)%var')) RETURN
        END IF
      END DO
!
      IF (Licetrc(ng)) THEN
        DO itrc=1,NT(ng)
          IF (LiceTavg(itrc,ng)) THEN
            IF (.not.destroy(ng, ICE_TAVG(itrc,ng)%var, MyFile,         &
     &                       __LINE__,'ICE_TAVG(itrc,ng)%var')) RETURN
          END IF
        END DO
      END IF
# endif
#endif
      IF (allocated(idStrc))   deallocate ( idStrc )
!
!-----------------------------------------------------------------------
!  Deallocate derived-type ICE and ICE_LOBC structures.
!-----------------------------------------------------------------------
!
      IF (allocated(ICE))      deallocate ( ICE )
      IF (allocated(ICE_LOBC)) deallocate ( ICE_LOBC )
#ifdef AVERAGES
      IF (allocated(ICE_FAVG)) deallocate ( ICE_FAVG )
      IF (allocated(ICE_SAVG)) deallocate ( ICE_SAVG )
      IF (allocated(ICE_TAVG)) deallocate ( ICE_TAVG )
#endif
!
!-----------------------------------------------------------------------
!  Deallocate other variables in modules.
!-----------------------------------------------------------------------
!
#ifdef AVERAGES
      IF (allocated(LiceFavg)) deallocate ( LiceFavg )
      IF (allocated(LiceSavg)) deallocate ( LiceSavg )
      IF (allocated(LiceTavg)) deallocate ( LiceTavg )
#endif
!
      RETURN
      END SUBROUTINE deallocate_ice
!
      SUBROUTINE initialize_ice (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param,   ONLY : BOUNDS, DOMAIN, LBC, iNLM, NT, MT
      USE mod_scalars, ONLY : iwest, ieast, isouth, inorth
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: i, j, nf, ns, ic, itrc
      integer :: Imin, Imax, Jmin, Jmax

      real(r8), parameter :: IniVal = 0.0_r8
!
#include "tile.h"
!
!  Set array initialization range.
!
#ifdef _OPENMP
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
#else
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
#endif
!
!-----------------------------------------------------------------------
!  Initialize ice model state (Si) and internal field (Fi) arrays.
!-----------------------------------------------------------------------
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          DO ns=1,nIceS
            ICE(ng) % Si(i,j,1,ns) = IniVal
            ICE(ng) % Si(i,j,2,ns) = IniVal
          END DO
!
          DO nf=1,nIceF
            ICE(ng) % Fi(i,j,nf) = IniVal
          END DO
        END DO
      END DO
!
      IF (Licetrc(ng)) THEN
        ic=0
        DO itrc=1,NT(ng)
          IF (LtracerIce(itrc,ng)) THEN
            ic=ic+1
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                ICE(ng) % Ti(i,j,1,ic) = IniVal
                ICE(ng) % Ti(i,j,2,ic) = IniVal
              END DO
            END DO
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Initialize ice model lateral boundary arrays.
!-----------------------------------------------------------------------
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN

        DO i=1,nIceS
          IF (DOMAIN(ng)%NorthWest_Test(tile)) THEN
            IF (LBC(iwest,ibICE(i),ng)%acquire) THEN
              ICE_LOBC(i,ng) % ice_west  = IniVal
              ICE_LOBC(i,ng) % iceG_west = IniVal
            END IF
          END IF
!
          IF (DOMAIN(ng)%SouthEast_Test(tile)) THEN
            IF (LBC(ieast,ibICE(i),ng)%acquire) THEN
              ICE_LOBC(i,ng) % ice_east  = IniVal
              ICE_LOBC(i,ng) % iceG_east = IniVal
            END IF
          END IF
!
          IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
            IF (LBC(isouth,ibICE(i),ng)%acquire) THEN
              ICE_LOBC(i,ng) % ice_south  = IniVal
              ICE_LOBC(i,ng) % iceG_south = IniVal
            END IF
          END IF
!
          IF (DOMAIN(ng)%NorthEast_Test(tile)) THEN
            IF (LBC(inorth,ibICE(i),ng)%acquire) THEN
              ICE_LOBC(i,ng) % ice_north  = IniVal
              ICE_LOBC(i,ng) % iceG_north = IniVal
            END IF
          END IF
        END DO

      END IF

#ifdef AVERAGES
!
!-----------------------------------------------------------------------
!  Initialize time-averaged variables.
!-----------------------------------------------------------------------
!
!  Ice model internal fields.
!
      DO nf=1,nIceF
        IF (iFice(nf).gt.0) THEN
          IF (LiceFavg(nf,ng)) THEN
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                ICE_FAVG(nf,ng) % var(i,j) = IniVal
              END DO
            END DO
          END IF
        END IF
      END DO
!
!  Time-averaged state fields.
!
      DO ns=1,nIceS
        IF (iSice(ns).gt.0) THEN
          IF (LiceSavg(ns,ng)) THEN
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                ICE_SAVG(ns,ng) % var(i,j) = IniVal
              END DO
            END DO
          END IF
        END IF
      END DO
!
!  Time-averaged tracer fields.
!
      IF (Licetrc(ng)) THEN
        DO itrc=1,NT(ng)
          IF (LiceTavg(itrc,ng)) THEN
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                ICE_TAVG(itrc,ng) % var(i,j) = IniVal
              END DO
            END DO
          END IF
        END DO
      END IF
#endif
!
      RETURN
      END SUBROUTINE initialize_ice
!
      END MODULE mod_ice
