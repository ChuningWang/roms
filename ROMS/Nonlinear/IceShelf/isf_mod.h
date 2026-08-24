      MODULE mod_isf
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2024 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!========================================== Benjamin K. Galton-Fenzi ===
!                                                                      !
!  Ice Shelf Model Kernel Variables:                                   !
!                                                                      !
!  This module declares Ice Shelf Model Kernel constants and structure !
!  containing several variables describing ice draft and thickness     !
!  properties.                                                         !
!                                                                      !
!  Variables:                                                          !
#ifdef ICESHEET
!                                                                      !
!  The iceshelf model prognostic 'state' variable are declares in      !
!  compact form in the derived-type structure as:                      !
!                                                                      !
!    ISF(ng) % state(:,:,TimeLevel,StateIndex)                         !
!                                                                      !
!  where the 'StateIndex' are as follows:                              !
!                                                                      !
!    isIhis    iceshelf thickness [m]                                  !
!    isIuis    iceshelf U-velocity [m s-1]                             !
!    isIvis    iceshelf V-velocity [m s-1]                             !
!    isIeis    iceshelf enthalpy [degC m]                              !
#endif
#ifdef ICESHELF_THERMO
!                                                                      !
!  The iceshelf model internal 'field' arrays are declares in compact  !
!  form in the derived-type structure as:                              !
!                                                                      !
!    ISF(ng) % field(:,:,FieldIndex)                                   !
!                                                                      !
!  where the 'FieldIndex' are as follows:                              !
!                                                                      !
!    icImel    iceshelf basal melt rate [m s-1]                        !
!    icIgat    heat turbulent exchange coefficient [m s-1]             !
!    icIgas    salt turbulent exchange coefficient [m s-1]             !
!    icItst    temperature diff between ambient/boundary layer [degC]  !
!    icIust    friction velocity in ice/water boundary layer [m s-1]   !
!    icItbo    temperature in laminar boundary layer [degC]            !
!    icIsbo    salinity in laminar boundary layer [PSU]                !
!                                                                      !
!  The iceshelf model refreezed volume and tracer concentration:       !
!                                                                      !
!    Vi        iceshelf freeezed volume (m)                            !
!    Ti        tracer concentration in refreezed volume (m TU)         !
#endif
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none

      PUBLIC :: allocate_isf
      PUBLIC :: deallocate_isf
      PUBLIC :: initialize_isf
!
!-----------------------------------------------------------------------
!  Iceshelf model identification indices.
!-----------------------------------------------------------------------
!
#ifdef ICESHEET
      integer :: idIhis       ! Iceshelf thickness
      integer :: idIuis       ! Iceshelf u-velocity
      integer :: idIvis       ! Iceshelf v-velocity
      integer :: idIeis       ! Iceshelf enthalpy
      integer :: idIuer       ! Iceshelf eastward velocity (output)
      integer :: idIvnr       ! Iceshelf northward velocity (output)
!
      integer, parameter :: nIsfS = 4    ! number of state variables
!
      integer, parameter :: isIhis =  1  ! Iceshelf thickness
      integer, parameter :: isIuis =  2  ! Iceshelf u-velocity
      integer, parameter :: isIvis =  3  ! Iceshelf v-velocity
      integer, parameter :: isIeis =  4  ! Iceshelf enthalpy
#else
      integer, parameter :: nIsfS = 0    ! number of ice state variables
#endif
!
!  Ice model state prognostic variables indices.
!
      integer :: iSisf(nIsfS)            ! state I/O indices
!
!  Iceshelf model state lateral boundary conditions indices.
!
      integer :: ibISF(nIsfS)            ! indices to LBC switch
      integer :: isfOBC(4,nIsfS)         ! I/O metadata indices
!
#ifdef ICESHELF_THERMO
!
      integer :: idImel       ! Iceshelf melt rate
      integer :: idIgat       ! gammaT
      integer :: idIgas       ! gammaS
      integer :: idItst       ! T*
      integer :: idIust       ! U*
      integer :: idItbo       ! Boundary layer temperature
      integer :: idIsbo       ! Boundary layer salinity
      integer :: idItba       ! Ambient boundary layer temperature
      integer :: idIsba       ! Ambient boundary layer salinity
      integer :: idIsbl       ! SBL vertical weight (for restart)
!
      integer, parameter :: nIsfF_thermo = 9  ! thermo field variables
!
      integer, parameter :: icImel =  1  ! Iceshelf melt rate
      integer, parameter :: icIgat =  2  ! gammaT
      integer, parameter :: icIgas =  3  ! gammaS
      integer, parameter :: icItst =  4  ! T*
      integer, parameter :: icIust =  5  ! U*
      integer, parameter :: icItbo =  6  ! Boundary layer temperature
      integer, parameter :: icIsbo =  7  ! Boundary layer salinity
      integer, parameter :: icItba =  8  ! Ambient T in boundary layer
      integer, parameter :: icIsba =  9  ! Ambient S in boundary layer
!
!  Iceshelf model ocean variable indices (SBL vertical weights).
!  Oi(i,j,k,ioHsbl) stores the fraction of ice shelf basal flux
!  assigned to vertical layer k within the sub-boundary layer.
!
      integer, parameter :: nIsfO = 1    ! number of ocean variables
!
      integer, parameter :: ioHsbl = 1   ! SBL vertical weight
#else
!
      integer, parameter :: nIsfF_thermo = 0  ! thermo field variables
!  Placeholders so that ICESHEET applications without ICESHELF_THERMO
!  (no SBL weights) still compile: iOisf/Oi become zero-sized.
      integer, parameter :: nIsfO      = 0    ! number of ocean variables
      integer, parameter :: ioHsbl     = 1    ! SBL vertical weight
#endif
#ifdef ICESHEET
!
      integer :: idItis       ! Iceshelf temperature
      integer :: idIzsf       ! Ice sheet surface elevation
      integer :: idIzbd       ! Ice sheet bed elevation
      integer :: idInuh       ! Depth-integrated effective viscosity
      integer :: idItdx       ! Driving stress x-component
      integer :: idItdy       ! Driving stress y-component
      integer :: idIbta       ! Basal drag coefficient
!
      integer, parameter :: nIsfF_ssa = 7    ! SSA diag field variables
!
# ifdef ICESHELF_THERMO
      integer, parameter :: icItis = 10  ! Iceshelf temperature
      integer, parameter :: icIzsf = 11  ! Ice sheet surface elevation
      integer, parameter :: icIzbd = 12  ! Ice sheet bed elevation
      integer, parameter :: icInuh = 13  ! Depth-inte. eff. viscosity
      integer, parameter :: icItdx = 14  ! Driving stress x-component
      integer, parameter :: icItdy = 15  ! Driving stress y-component
      integer, parameter :: icIbta = 16  ! Basal drag coefficient
# else
      integer, parameter :: icItis =  1  ! Iceshelf temperature
      integer, parameter :: icIzsf =  2  ! Ice sheet surface elevation
      integer, parameter :: icIzbd =  3  ! Ice sheet bed elevation
      integer, parameter :: icInuh =  4  ! Depth-inte. eff. viscosity
      integer, parameter :: icItdx =  5  ! Driving stress x-component
      integer, parameter :: icItdy =  6  ! Driving stress y-component
      integer, parameter :: icIbta =  7  ! Basal drag coefficient
# endif
#else
!
      integer, parameter :: nIsfF_ssa = 0    ! SSA diag field variables
#endif
!
      integer, parameter :: nIsfF = nIsfF_thermo + nIsfF_ssa
!
      integer :: iFisf(nIsfF)            ! internal fields I/O indices
!
      integer :: iOisf(nIsfO)            ! ocean variables I/O indices
!
      integer :: isf_ic_brylast          ! Boundary layer salinity
!
      integer              :: idIvol     ! Volume of refreezed water
      integer, allocatable :: idItrc(:)  ! Tracer variables
      integer, allocatable :: idIamb(:)  ! Ambient tracer variables
!
!-----------------------------------------------------------------------
!  Iceshelf model parameters.
!                                                                      !
!   a, b, c  - Freezing point EOS constants [degC PSU-1, N/A, degC m-1]!
!   L        - Latent heat of melting [J kg-1]                         !
#ifdef ICESHELF_THERMO
!   Pr, Sc   - Prandtl number and Schmidt number                       !
!   Cd       - Drag coefficient                                        !
!   visc     - Background viscosity                                    !
!   L        - Latent heat of melting [J kg-1]                         !
!   rho_i    - Ice density [kg m-3]                                    !
!   Tis, Sis - Ice temperature/salinity                                !
#endif
!-----------------------------------------------------------------------
!
#ifdef AVERAGES
!
!  Switches to process time-averaged iceshelf model state and internal
!  variables.
!
      logical, allocatable :: LisfFavg(:,:)     ! internal variables
      logical, allocatable :: LisfSavg(:,:)     ! state variables
      logical, allocatable :: LisfVavg(:)       ! refreezed volume
      logical, allocatable :: LisfTavg(:,:)     ! refreezed tracer
      logical, allocatable :: LisfAavg(:,:)     ! ambient tracer
#endif
!
!  Gamma constraining factors
!
      real(r8), allocatable :: isf_drdz(:)      ! Drho/Dz value
      real(r8), allocatable :: gamU(:)          ! gamma scaling factor
      real(r8), allocatable :: min_ustar(:)     ! minimum ustar value
      real(r8), allocatable :: isf_Hsbl(:)      ! SBL thickness [m]
!
!  Ice-ocean coupling parameters
!
      real(r8), allocatable :: tau_melt(:)      ! MR EMA filter TS [s]
#ifdef ICESHELF
      real(r8), allocatable :: rdrg_isf(:)      ! linear drag coeff. [m/s]
      real(r8), allocatable :: rdrg2_isf(:)     ! quadratic drag coeff. [-]
      real(r8), allocatable :: ZoSr_isf(:)     ! surface roughness [m]
#endif
!
!  Logical switches
!
      integer, allocatable :: NTISF(:)
      logical, allocatable :: Lisftrc(:)
      logical, allocatable :: LtracerIsf(:,:)   ! tracer exchange switch
      logical, allocatable :: Lisfsrc(:)
      logical, allocatable :: LsourceIsf(:,:)   ! tracer source switch
!
!  Iceshelf melt water tracer concentration
!
      real(r8), allocatable :: trcIsf(:,:)      ! tracer concentration
!
#ifdef ICESHEET
!  SSA physical parameters
!
      real(r8), allocatable :: n_glen(:)        ! Glen flow law exponent
      real(r8), allocatable :: A_glen(:)        ! Glen rate factor
      real(r8), allocatable :: eps_glen_min(:)  ! minimum strain rate
      real(r8), allocatable :: C_basal(:)       ! basal drag coefficient
      real(r8), allocatable :: n_basal(:)       ! basal sliding exponent
      real(r8), allocatable :: eps_u_min(:)     ! minimum sliding vel.
!
!  SSA solver parameters
!
      integer, allocatable :: nl_max_iter(:)    ! max Picard iterations
      real(r8), allocatable :: nl_tol(:)        ! Picard relative tol.
      real(r8), allocatable :: nl_tol_fp(:)     ! Picard vel. tol.
      integer, allocatable :: cg_max_iter(:)    ! max CG iterations
      real(r8), allocatable :: cg_tol(:)        ! CG relative tolerance
      real(r8), allocatable :: nuH_reg(:)       ! nuH regularization
      real(r8), allocatable :: min_thickness(:) ! minimum ice thickness
      real(r8), allocatable :: ssa_max_speed(:) ! velocity cap [m/s]
!
!  Driving stress regularization parameters
!
      integer, allocatable :: n_smooth(:)       ! elev. smo. iter.
      integer, allocatable :: n_taper(:)        ! ice margin taper width
!
!  Picard iteration parameters
!
      real(r8), allocatable :: picard_alpha(:)  ! Picard relax. factor
!
!  Thickness advection parameters
!
      real(r8), allocatable :: cfl_factor(:)    ! CFL safety factor
!
!  Calving parameters
!
      real(r8), allocatable :: H_calve(:)       ! Thick. thre. calv. [m]
!
!  Multi-timescale parameters
!
      integer, allocatable :: N_ssa(:)          ! SSA solve freq.
      integer, allocatable :: N_gl(:)           ! GL update freq.
!
!  Surface mass balance parameters
!
      real(r8), allocatable :: smb_rate(:)      ! Uniform SMB rate [m/s]
!
!  SSA geometry parameters
!
      real(r8), allocatable :: max_thickness(:) ! thick. safety cap [m]
      real(r8), allocatable :: gl_float_width(:)! GL floatfrac width [m]
#endif
!
!  Other constant parameters
!
      real(r8), parameter :: min_his = 10.0_r8
      real(r8), parameter :: a       = -0.057_r8
      real(r8), parameter :: b       = 0.0939_r8
      real(r8), parameter :: c       = -7.61e-4_r8
      real(r8), parameter :: L       = 3.33e5_r8
      real(r8), parameter :: rho_i   = 916.7_r8
#ifdef ICESHELF_THERMO
      real(r8), parameter :: Pr      = 13.8_r8
      real(r8), parameter :: Sc      = 2432.2_r8
      real(r8), parameter :: Cd      = 5.0e-3_r8
      real(r8), parameter :: visc    = 1.95e-6_r8
!     real(r8), parameter :: small   = 1.0e-3_r8 ! Diffusion + tidal
      real(r8), parameter :: small   = 8.0e-6_r8 ! Diffusion only
      real(r8), parameter :: Tis     = -20.0_r8
      real(r8), parameter :: Sis     = 0.0_r8
#endif
      real(r8), parameter :: eps     = 1.0E-14_r8
!
!-----------------------------------------------------------------------
!  Define derived-type structure iceshelf model state and internal
!  arrays.
!-----------------------------------------------------------------------
!
      TYPE T_ISF

        real(r8), pointer :: Fi(:,:,:)    ! [i, j,         1:nIsfF]
        real(r8), pointer :: Si(:,:,:,:)  ! [i, j,    1:2, 1:nIsfS]
        real(r8), pointer :: Vi(:,:,:)    ! [i, j,    1:2         ]
        real(r8), pointer :: Ti(:,:,:,:)  ! [i, j,    1:2, 1:NTISF]
        real(r8), pointer :: Oi(:,:,:,:)  ! [i, j, k,      1:nIsfO]
        real(r8), pointer :: Ai(:,:,:)    ! [i, j,         1:NTISF]

      END TYPE T_ISF
!
      TYPE (T_ISF), allocatable :: ISF(:)            ! [Ngrids]
#ifdef ICESHEET
!
!-----------------------------------------------------------------------
!  Define derived-type structure iceshelf model lateral boundary
!  variables.
!-----------------------------------------------------------------------
!
      TYPE T_ISF_LOBC

        real(r8), pointer :: isf_west (:)
        real(r8), pointer :: isf_east (:)
        real(r8), pointer :: isf_south(:)
        real(r8), pointer :: isf_north(:)
!
        real(r8), pointer :: isfG_west (:,:)
        real(r8), pointer :: isfG_east (:,:)
        real(r8), pointer :: isfG_south(:,:)
        real(r8), pointer :: isfG_north(:,:)

      END TYPE T_ISF_LOBC
!
      TYPE (T_ISF_LOBC), allocatable :: ISF_LOBC(:,:)   ! [nIsfS,Ngrids]
#endif

#ifdef AVERAGES
!
!-----------------------------------------------------------------------
!  Define derived-type structure iceself model state and internal arrays
!  time-averaged variables. Notice that only the requested arrays are
!  allocated and processed.
!-----------------------------------------------------------------------
!
      TYPE T_ISF_AVG

        real(r8), pointer :: var(:,:)                   ! [i,j]

      END TYPE T_ISF_AVG
!
      TYPE (T_ISF_AVG), allocatable :: ISF_FAVG(:,:)    ! [nIsfF,Ngrids]
      TYPE (T_ISF_AVG), allocatable :: ISF_SAVG(:,:)    ! [nIsfS,Ngrids]
      TYPE (T_ISF_AVG), allocatable :: ISF_VAVG(:)      ! [Ngrids]
      TYPE (T_ISF_AVG), allocatable :: ISF_TAVG(:,:)    ! [MT,Ngrids]
      TYPE (T_ISF_AVG), allocatable :: ISF_AAVG(:,:)    ! [MT,Ngrids]
#endif
!
      CONTAINS
!
      SUBROUTINE allocate_isf (ng, LBi, UBi, LBj, UBj, isf_kernel)
!
!=======================================================================
!                                                                      !
!  This routine allocates either the derived-type iceshelf model       !
!  variables (isf_kernel=TRUE) or module parameters (isf_kernel=FALSE).!
!                                                                      !
!=======================================================================
!
      USE mod_param,   ONLY : BOUNDS, Dmem, LBC, Ngrids, NT, MT, NAT, N
      USE mod_scalars, ONLY : iwest, ieast, isouth, inorth
!
!  Imported variable declarations.
!
      logical, intent(in) :: isf_kernel
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: LBiG, UBiG, LBjG, UBjG, my_tile
      integer :: i, ic, itrc

      real(r8), parameter :: IniVal = 0.0_r8
!
      real(r8) :: size2d, Xsize, Ysize
!
      size2d=REAL((UBi-LBi+1)*(UBj-LBj+1),r8)
      Xsize =REAL(UBi-LBi,r8)
      Ysize =REAL(UBj-LBj,r8)

!
!-----------------------------------------------------------------------
!  Allocate iceshelf model parameters.
!-----------------------------------------------------------------------
!
      IF (.not.isf_kernel) THEN

        IF (.not.allocated(idItrc)) THEN
          allocate ( idItrc(MT) )
          Dmem(1)=Dmem(1)+REAL(MT,r8)
          idItrc(1:MT)=0
        END IF

        IF (.not.allocated(idIamb)) THEN
          allocate ( idIamb(MT) )
          Dmem(1)=Dmem(1)+REAL(MT,r8)
          idIamb(1:MT)=0
        END IF

#ifdef AVERAGES

        IF (.not.allocated(LisfFavg)) THEN
          allocate ( LisfFavg(nIsfF,Ngrids) )
          Dmem(1)=Dmem(1)+REAL(nIsfF*Ngrids,r8)
          LisfFavg(1:nIsfF,1:Ngrids)=.FALSE.
        END IF

        IF (.not.allocated(LisfSavg)) THEN
          allocate ( LisfSavg(nIsfS,Ngrids) )
          Dmem(1)=Dmem(1)+REAL(nIsfS*Ngrids,r8)
          LisfSavg(1:nIsfS,1:Ngrids)=.FALSE.
        END IF

        IF (.not.allocated(LisfVavg)) THEN
          allocate ( LisfVavg(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          LisfVavg(1:Ngrids)=.FALSE.
        END IF

        IF (.not.allocated(LisfTavg)) THEN
          allocate ( LisfTavg(MT,Ngrids) )
          Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
          LisfTavg(1:MT,1:Ngrids)=.FALSE.
        END IF

        IF (.not.allocated(LisfAavg)) THEN
          allocate ( LisfAavg(MT,Ngrids) )
          Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
          LisfAavg(1:MT,1:Ngrids)=.FALSE.
        END IF
#endif

        IF (.not.allocated(isf_drdz)) THEN
          allocate ( isf_drdz(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          isf_drdz(1:Ngrids)=IniVal
        END IF

        IF (.not.allocated(gamU)) THEN
          allocate ( gamU(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          gamU(1:Ngrids)=IniVal
        END IF

        IF (.not.allocated(min_ustar)) THEN
          allocate ( min_ustar(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          min_ustar(1:Ngrids)=IniVal
        END IF

        IF (.not.allocated(isf_Hsbl)) THEN
          allocate ( isf_Hsbl(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          isf_Hsbl(1:Ngrids)=IniVal
        END IF
        IF (.not.allocated(tau_melt)) THEN
          allocate ( tau_melt(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          tau_melt(1:Ngrids)=IniVal
        END IF
#ifdef ICESHELF
        IF (.not.allocated(rdrg_isf)) THEN
          allocate ( rdrg_isf(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          rdrg_isf(1:Ngrids)=IniVal
        END IF
        IF (.not.allocated(rdrg2_isf)) THEN
          allocate ( rdrg2_isf(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          rdrg2_isf(1:Ngrids)=IniVal
        END IF
        IF (.not.allocated(ZoSr_isf)) THEN
          allocate ( ZoSr_isf(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          ZoSr_isf(1:Ngrids)=IniVal
        END IF
#endif

        IF (.not.allocated(NTISF)) THEN
          allocate ( NTISF(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          NTISF(1:Ngrids)=0
        END IF

        IF (.not.allocated(Lisftrc)) THEN
          allocate ( Lisftrc(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          Lisftrc(1:Ngrids)=.FALSE.
        END IF

        IF (.not.allocated(LtracerIsf)) THEN
          allocate ( LtracerIsf(MT,Ngrids) )
          Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
          LtracerIsf(1:MT,1:Ngrids)=.FALSE.
        END IF

        IF (.not.allocated(Lisfsrc)) THEN
          allocate ( Lisfsrc(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          Lisfsrc(1:Ngrids)=.FALSE.
        END IF

        IF (.not.allocated(LsourceIsf)) THEN
          allocate ( LsourceIsf(MT,Ngrids) )
          Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
          LsourceIsf(1:MT,1:Ngrids)=.FALSE.
        END IF

        IF (.not.allocated(trcIsf)) THEN
          allocate ( trcIsf(MT,Ngrids) )
          Dmem(1)=Dmem(1)+REAL(MT*Ngrids,r8)
          trcIsf(1:MT,1:Ngrids)=IniVal
        END IF

#ifdef ICESHEET
!
!  Allocate SSA physical and solver parameters.
!
        IF (.not.allocated(n_glen)) THEN
          allocate ( n_glen(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          n_glen(1:Ngrids)=IniVal
        END IF
        IF (.not.allocated(A_glen)) THEN
          allocate ( A_glen(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          A_glen(1:Ngrids)=IniVal
        END IF
        IF (.not.allocated(eps_glen_min)) THEN
          allocate ( eps_glen_min(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          eps_glen_min(1:Ngrids)=IniVal
        END IF
        IF (.not.allocated(C_basal)) THEN
          allocate ( C_basal(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          C_basal(1:Ngrids)=IniVal
        END IF
        IF (.not.allocated(n_basal)) THEN
          allocate ( n_basal(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          n_basal(1:Ngrids)=IniVal
        END IF
        IF (.not.allocated(eps_u_min)) THEN
          allocate ( eps_u_min(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          eps_u_min(1:Ngrids)=IniVal
        END IF
        IF (.not.allocated(nl_max_iter)) THEN
          allocate ( nl_max_iter(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          nl_max_iter(1:Ngrids)=0
        END IF
        IF (.not.allocated(nl_tol)) THEN
          allocate ( nl_tol(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          nl_tol(1:Ngrids)=IniVal
        END IF
        IF (.not.allocated(nl_tol_fp)) THEN
          allocate ( nl_tol_fp(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          nl_tol_fp(1:Ngrids)=IniVal
        END IF
        IF (.not.allocated(cg_max_iter)) THEN
          allocate ( cg_max_iter(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          cg_max_iter(1:Ngrids)=0
        END IF
        IF (.not.allocated(cg_tol)) THEN
          allocate ( cg_tol(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          cg_tol(1:Ngrids)=IniVal
        END IF
        IF (.not.allocated(nuH_reg)) THEN
          allocate ( nuH_reg(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          nuH_reg(1:Ngrids)=IniVal
        END IF
        IF (.not.allocated(min_thickness)) THEN
          allocate ( min_thickness(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          min_thickness(1:Ngrids)=IniVal
        END IF
        IF (.not.allocated(ssa_max_speed)) THEN
          allocate ( ssa_max_speed(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          ssa_max_speed(1:Ngrids)=IniVal
        END IF
        IF (.not.allocated(n_smooth)) THEN
          allocate ( n_smooth(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          n_smooth(1:Ngrids)=0
        END IF
        IF (.not.allocated(n_taper)) THEN
          allocate ( n_taper(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          n_taper(1:Ngrids)=0
        END IF
        IF (.not.allocated(picard_alpha)) THEN
          allocate ( picard_alpha(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          picard_alpha(1:Ngrids)=IniVal
        END IF
        IF (.not.allocated(cfl_factor)) THEN
          allocate ( cfl_factor(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          cfl_factor(1:Ngrids)=IniVal
        END IF
        IF (.not.allocated(H_calve)) THEN
          allocate ( H_calve(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          H_calve(1:Ngrids)=IniVal
        END IF
        IF (.not.allocated(N_ssa)) THEN
          allocate ( N_ssa(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          N_ssa(1:Ngrids)=0
        END IF
        IF (.not.allocated(N_gl)) THEN
          allocate ( N_gl(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          N_gl(1:Ngrids)=0
        END IF
        IF (.not.allocated(smb_rate)) THEN
          allocate ( smb_rate(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          smb_rate(1:Ngrids)=IniVal
        END IF
!
!  SSA geometry parameters
!
        IF (.not.allocated(max_thickness)) THEN
          allocate ( max_thickness(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          max_thickness(1:Ngrids)=IniVal
        END IF
        IF (.not.allocated(gl_float_width)) THEN
          allocate ( gl_float_width(Ngrids) )
          Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
          gl_float_width(1:Ngrids)=IniVal
        END IF
#endif

      ELSE
!
!-----------------------------------------------------------------------
!  Allocate derived-type structure iceshelf model kernel variables.
!-----------------------------------------------------------------------
!
        IF (ng.eq.1) allocate ( ISF(Ngrids) )
!
!  Nonlinear iceshelf model state variables (Si) and internal kernel
!  field (Fi) arrays.
!
        allocate ( ISF(ng) % Si(LBi:UBi,LBj:UBj,2,nIsfS) )
        Dmem(ng)=Dmem(ng)+2.0_r8*REAL(nIsfS,r8)*size2d
!
        allocate ( ISF(ng) % Fi(LBi:UBi,LBj:UBj,nIsfF) )
        Dmem(ng)=Dmem(ng)+REAL(nIsfF,r8)*size2d
!
        allocate ( ISF(ng) % Oi(LBi:UBi,LBj:UBj,N(ng),nIsfO) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng)*nIsfO,r8)*size2d
!
        IF (Lisftrc(ng)) THEN
          allocate ( ISF(ng) % Vi(LBi:UBi,LBj:UBj,2) )
          Dmem(ng)=Dmem(ng)+2.0_r8*size2d
!
          allocate ( ISF(ng) % Ti(LBi:UBi,LBj:UBj,2,NTISF(ng)) )
          Dmem(ng)=Dmem(ng)+2.0_r8*REAL(NTISF(ng))*size2d
!
          allocate ( ISF(ng) % Ai(LBi:UBi,LBj:UBj,NTISF(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(NTISF(ng),r8)*size2d
        END IF
#ifdef ICESHEET
!
!-----------------------------------------------------------------------
!  Allocate derived-type structure iceshelf model lateral boundary
!  variables.
!-----------------------------------------------------------------------
!
        my_tile=-1                           ! for global values
        LBiG=BOUNDS(ng)%LBi(my_tile)
        UBiG=BOUNDS(ng)%UBi(my_tile)
        LBjG=BOUNDS(ng)%LBj(my_tile)
        UBjG=BOUNDS(ng)%UBj(my_tile)
!
        IF (ng.eq.1) allocate ( ISF_LOBC(nIsfS,Ngrids) )
!
!  Nonlinear iceself model 'state' lateral boundary arrays.
!
        DO i=1,nIsfS
          IF (LBC(iwest,ibISF(i),ng)%acquire) THEN
            allocate ( ISF_LOBC(i,ng) % isf_west(LBjG:UBjG) )
            Dmem(ng)=Dmem(ng)+Ysize

            allocate ( ISF_LOBC(i,ng) % isfG_west(LBjG:UBjG,2) )
            Dmem(ng)=Dmem(ng)+2.0_r8*Ysize
          END IF
!
          IF (LBC(ieast,ibISF(i),ng)%acquire) THEN
            allocate ( ISF_LOBC(i,ng) % isf_east(LBjG:UBjG) )
            Dmem(ng)=Dmem(ng)+Ysize

            allocate ( ISF_LOBC(i,ng) % isfG_east(LBjG:UBjG,2) )
            Dmem(ng)=Dmem(ng)+2.0_r8*Ysize
          END IF
!
          IF (LBC(isouth,ibISF(i),ng)%acquire) THEN
            allocate ( ISF_LOBC(i,ng) % isf_south(LBiG:UBiG) )
            Dmem(ng)=Dmem(ng)+Xsize

            allocate ( ISF_LOBC(i,ng) % isfG_south(LBiG:UBiG,2) )
            Dmem(ng)=Dmem(ng)+2.0_r8*Xsize
          END IF
!
          IF (LBC(inorth,ibISF(i),ng)%acquire) THEN
            allocate ( ISF_LOBC(i,ng) % isf_north(LBiG:UBiG) )
            Dmem(ng)=Dmem(ng)+Xsize

            allocate ( ISF_LOBC(i,ng) % isfG_north(LBiG:UBiG,2) )
            Dmem(ng)=Dmem(ng)+2.0_r8*Xsize
          END IF
        END DO
#endif
#ifdef AVERAGES
!
!-----------------------------------------------------------------------
!  Allocate derived-type structure iceshelf model state and internal
!  arrays time-averaged variables.
!-----------------------------------------------------------------------
!
        IF (ng.eq.1) THEN
          allocate ( ISF_FAVG(nIsfF,Ngrids) )
          allocate ( ISF_SAVG(nIsfS,Ngrids) )
        END IF
!
!  Time-averaged internal fields. Only fields with metadata (iFisf > 0)
!  are allocated and processed.
!
        DO i=1,nIsfF
          IF (iFisf(i).gt.0) THEN
            IF (LisfFavg(i,ng)) THEN
              allocate ( ISF_FAVG(i,ng) % var(LBi:UBi,LBj:UBj) )
              Dmem(ng)=Dmem(ng)+size2d
            END IF
          END IF
        END DO
!
!  Time-averaged state fields. Only fields with metadata (iSisf > 0)
!  are allocated and processed.
!
        DO i=1,nIsfS
          IF (iSisf(i).gt.0) THEN
            IF (LisfSavg(i,ng)) THEN
              allocate ( ISF_SAVG(i,ng) % var(LBi:UBi,LBj:UBj) )
              Dmem(ng)=Dmem(ng)+size2d
            END IF
          END IF
        END DO
!
!  Time-averaged tracer fields.
!
        IF (ng.eq.1) THEN
          allocate ( ISF_VAVG(Ngrids) )
          allocate ( ISF_TAVG(MT,Ngrids) )
          allocate ( ISF_AAVG(MT,Ngrids) )
        END IF
!
        IF (LisfVavg(ng)) THEN
          allocate ( ISF_VAVG(ng) % var(LBi:UBi,LBj:UBj) )
          Dmem(ng)=Dmem(ng)+size2d
        END IF
!
        IF (Lisftrc(ng)) THEN
          DO itrc=NAT+1,NT(ng)
            IF (LisfTavg(itrc,ng)) THEN
              allocate ( ISF_TAVG(itrc,ng) % var(LBi:UBi,LBj:UBj) )
              Dmem(ng)=Dmem(ng)+size2d
            END IF
            IF (LisfAavg(itrc,ng)) THEN
              allocate ( ISF_AAVG(itrc,ng) % var(LBi:UBi,LBj:UBj) )
              Dmem(ng)=Dmem(ng)+size2d
            END IF
          END DO
        END IF
#endif
      END IF
!
      RETURN
      END SUBROUTINE allocate_isf
!
      SUBROUTINE deallocate_isf (ng)
!
!=======================================================================
!                                                                      !
!  This routine deallocates all variables in the module for all nested !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param,   ONLY : Ngrids, MT
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
      integer :: i, itrc
!
#endif
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", deallocate_isf"
#ifdef SUBOBJECT_DEALLOCATION
!
!-----------------------------------------------------------------------
!  Deallocate each variable in the derived-type ISF structure
!  separately.
!-----------------------------------------------------------------------
!
!  Nonlinear iceshelf model state variables 'Si' and internal kernel
!  arrays 'Fi'.
!
      IF (.not.destroy(ng, ISF(ng)%Si, MyFile,                          &
     &                 __LINE__, 'ISF(ng)%Si')) RETURN
!
      IF (.not.destroy(ng, ISF(ng)%Fi, MyFile,                          &
     &                 __LINE__, 'ISF(ng)%Fi')) RETURN
!
      IF (.not.destroy(ng, ISF(ng)%Oi, MyFile,                          &
     &                 __LINE__, 'ISF(ng)%Oi')) RETURN
!
      IF (Lisftrc(ng)) THEN
        IF (.not.destroy(ng, ISF(ng)%Vi, MyFile,                        &
     &                   __LINE__, 'ISF(ng)%Vi')) RETURN
!
        IF (.not.destroy(ng, ISF(ng)%Ti, MyFile,                        &
     &                   __LINE__, 'ISF(ng)%Ti')) RETURN
!
        IF (.not.destroy(ng, ISF(ng)%Ai, MyFile,                        &
     &                   __LINE__, 'ISF(ng)%Ai')) RETURN
      END IF
# ifdef ICESHEET
!
!  Nonlinear iceshelf model state lateral boundary arrays.
!
      DO i=1,nIsfS
        IF (LBC(iwest,ibISF(i),ng)%acquire) THEN
          IF (.not.destroy(ng, ISF_LOBC(i,ng)%isf_west, MyFile,         &
     &                     __LINE__, 'ISF_LOBC(i,ng)%isf_west')) RETURN
          IF (.not.destroy(ng, ISF_LOBC(i,ng)%isfG_west, MyFile,        &
     &                     __LINE__, 'ISF_LOBC(i,ng)%isfG_west')) RETURN
        END IF
!
        IF (LBC(ieast,ibISF(i),ng)%acquire) THEN
          IF (.not.destroy(ng, ISF_LOBC(i,ng)%isf_east, MyFile,         &
     &                     __LINE__, 'ISF_LOBC(i,ng)%isf_west')) RETURN
          IF (.not.destroy(ng, ISF_LOBC(i,ng)%isfG_east, MyFile,        &
     &                     __LINE__, 'ISF_LOBC(i,ng)%isfG_west')) RETURN
        END IF
!
        IF (LBC(isouth,ibISF(i),ng)%acquire) THEN
          IF (.not.destroy(ng, ISF_LOBC(i,ng)%isf_south, MyFile,        &
     &                     __LINE__, 'ISF_LOBC(i,ng)%isf_south')) RETURN
          IF (.not.destroy(ng, ISF_LOBC(i,ng)%isfG_south, MyFile,       &
     &                     __LINE__,'ISF_LOBC(i,ng)%isfG_south')) RETURN
        END IF
!
        IF (LBC(inorth,ibISF(i),ng)%acquire) THEN
          IF (.not.destroy(ng, ISF_LOBC(i,ng)%isf_north, MyFile,        &
     &                     __LINE__, 'ISF_LOBC(i,ng)%isf_north')) RETURN
          IF (.not.destroy(ng, ISF_LOBC(i,ng)%isfG_south, MyFile,       &
     &                     __LINE__,'ISF_LOBC(i,ng)%isfG_north')) RETURN
        END IF
      END DO
# endif

# ifdef AVERAGES
!
!-----------------------------------------------------------------------
!  Time-averaged iceshelf model state and internal variables.
!-----------------------------------------------------------------------
!
      DO i=1,nIsfF
        IF (LisfFavg(i,ng)) THEN
          IF (.not.destroy(ng, ISF_FAVG(i,ng)%var, MyFile,              &
     &                     __LINE__, 'ISF_FAVG(i,ng)%var')) RETURN
        END IF
      END DO
!
      DO i=1,nIsfS
        IF (LisfSavg(i,ng)) THEN
          IF (.not.destroy(ng, ISF_SAVG(i,ng)%var, MyFile,              &
     &                     __LINE__, 'ISF_SAVG(i,ng)%var')) RETURN
        END IF
      END DO
!
      IF (LisfVavg(ng)) THEN
        IF (.not.destroy(ng, ISF_VAVG(ng)%var, MyFile,                  &
     &                   __LINE__, 'ISF_VAVG(ng)%var')) RETURN
      END IF
!
      IF (Lisftrc(ng)) THEN
        DO itrc=1,NT(ng)
          IF (LisfTavg(itrc,ng)) THEN
            IF (.not.destroy(ng, ISF_TAVG(itrc,ng)%var, MyFile,         &
     &                       __LINE__,'ISF_TAVG(itrc,ng)%var')) RETURN
          END IF
          IF (LisfAavg(itrc,ng)) THEN
            IF (.not.destroy(ng, ISF_AAVG(itrc,ng)%var, MyFile,         &
     &                       __LINE__,'ISF_AAVG(itrc,ng)%var')) RETURN
          END IF
        END DO
      END IF
# endif
#endif
      IF (allocated(idItrc))   deallocate ( idItrc )
      IF (allocated(idIamb))   deallocate ( idIamb )
!
!-----------------------------------------------------------------------
!  Deallocate derived-type ISF and ISF_LOBC structures.
!-----------------------------------------------------------------------
!
      IF (allocated(ISF))      deallocate ( ISF )
#ifdef ICESHEET
      IF (allocated(ISF_LOBC)) deallocate ( ISF_LOBC )
#endif
#ifdef AVERAGES
      IF (allocated(ISF_FAVG)) deallocate ( ISF_FAVG )
      IF (allocated(ISF_SAVG)) deallocate ( ISF_SAVG )
      IF (allocated(ISF_VAVG)) deallocate ( ISF_VAVG )
      IF (allocated(ISF_TAVG)) deallocate ( ISF_TAVG )
      IF (allocated(ISF_AAVG)) deallocate ( ISF_AAVG )
#endif
!
!-----------------------------------------------------------------------
!  Deallocate other variables in modules.
!-----------------------------------------------------------------------
!
#ifdef AVERAGES
      IF (allocated(LisfFavg)) deallocate ( LisfFavg )
      IF (allocated(LisfSavg)) deallocate ( LisfSavg )
      IF (allocated(LisfVavg)) deallocate ( LisfVavg )
      IF (allocated(LisfTavg)) deallocate ( LisfTavg )
      IF (allocated(LisfAavg)) deallocate ( LisfAavg )
#endif
!
      IF (allocated(isf_Hsbl))  deallocate ( isf_Hsbl )
      IF (allocated(tau_melt))  deallocate ( tau_melt )
#ifdef ICESHELF
      IF (allocated(rdrg_isf))  deallocate ( rdrg_isf )
      IF (allocated(rdrg2_isf)) deallocate ( rdrg2_isf )
      IF (allocated(ZoSr_isf)) deallocate ( ZoSr_isf )
#endif
!
#ifdef ICESHEET
      IF (allocated(n_glen))         deallocate ( n_glen )
      IF (allocated(A_glen))         deallocate ( A_glen )
      IF (allocated(eps_glen_min))   deallocate ( eps_glen_min )
      IF (allocated(C_basal))        deallocate ( C_basal )
      IF (allocated(n_basal))        deallocate ( n_basal )
      IF (allocated(eps_u_min))      deallocate ( eps_u_min )
      IF (allocated(nl_max_iter))    deallocate ( nl_max_iter )
      IF (allocated(nl_tol))         deallocate ( nl_tol )
      IF (allocated(nl_tol_fp))      deallocate ( nl_tol_fp )
      IF (allocated(cg_max_iter))    deallocate ( cg_max_iter )
      IF (allocated(cg_tol))         deallocate ( cg_tol )
      IF (allocated(nuH_reg))        deallocate ( nuH_reg )
      IF (allocated(min_thickness))  deallocate ( min_thickness )
      IF (allocated(ssa_max_speed))  deallocate ( ssa_max_speed )
      IF (allocated(n_smooth))       deallocate ( n_smooth )
      IF (allocated(n_taper))        deallocate ( n_taper )
      IF (allocated(picard_alpha))   deallocate ( picard_alpha )
      IF (allocated(cfl_factor))     deallocate ( cfl_factor )
      IF (allocated(H_calve))        deallocate ( H_calve )
      IF (allocated(N_ssa))          deallocate ( N_ssa )
      IF (allocated(N_gl))           deallocate ( N_gl )
      IF (allocated(smb_rate))       deallocate ( smb_rate )
      IF (allocated(max_thickness))   deallocate ( max_thickness )
      IF (allocated(gl_float_width))  deallocate ( gl_float_width )
#endif
!
      RETURN
      END SUBROUTINE deallocate_isf
!
      SUBROUTINE initialize_isf (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine initialize structure variables in the module using     !
!  first touch distribution policy. In shared-memory configuration,    !
!  this operation actually performs the propagation of the "shared     !
!  arrays" across the cluster,  unless another policy is specified     !
!  to  override the default.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param,   ONLY : BOUNDS, DOMAIN, LBC, iNLM, NT, MT, N
      USE mod_scalars, ONLY : iwest, ieast, isouth, inorth
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: i, j, nf, ns, ic, itrc, k, no
      integer :: Imin, Imax, Jmin, Jmax

      real(r8), parameter :: IniVal = 0.0_r8

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
!  Initialize iceshelf model state (Si) and internal field (Fi) arrays.
!-----------------------------------------------------------------------
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          DO ns=1,nIsfS
            ISF(ng) % Si(i,j,1,ns) = IniVal
            ISF(ng) % Si(i,j,2,ns) = IniVal
          END DO
!
          DO nf=1,nIsfF
            ISF(ng) % Fi(i,j,nf) = IniVal
          END DO
!
          DO no=1,nIsfO
            DO k=1,N(ng)
              ISF(ng) % Oi(i,j,k,no) = IniVal
            END DO
          END DO
        END DO
      END DO
	!
      IF (Lisftrc(ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            ISF(ng) % Vi(i,j,1) = IniVal
            ISF(ng) % Vi(i,j,2) = IniVal
          END DO
        END DO
        ic=0
        DO itrc=1,NT(ng)
          IF (LtracerIsf(itrc,ng)) THEN
            ic=ic+1
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                ISF(ng) % Ti(i,j,1,ic) = IniVal
                ISF(ng) % Ti(i,j,2,ic) = IniVal
                ISF(ng) % Ai(i,j,  ic) = IniVal
              END DO
            END DO
          END IF
        END DO
      END IF
#ifdef ICESHEET
!
!-----------------------------------------------------------------------
!  Initialize iceshelf model lateral boundary arrays.
!-----------------------------------------------------------------------
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN

        DO i=1,nIsfS
          IF (DOMAIN(ng)%NorthWest_Test(tile)) THEN
            IF (LBC(iwest,ibISF(i),ng)%acquire) THEN
              ISF_LOBC(i,ng) % isf_west  = IniVal
              ISF_LOBC(i,ng) % isfG_west = IniVal
            END IF
          END IF
!
          IF (DOMAIN(ng)%SouthEast_Test(tile)) THEN
            IF (LBC(ieast,ibISF(i),ng)%acquire) THEN
              ISF_LOBC(i,ng) % isf_east  = IniVal
              ISF_LOBC(i,ng) % isfG_east = IniVal
            END IF
          END IF
!
          IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
            IF (LBC(isouth,ibISF(i),ng)%acquire) THEN
              ISF_LOBC(i,ng) % isf_south  = IniVal
              ISF_LOBC(i,ng) % isfG_south = IniVal
            END IF
          END IF
!
          IF (DOMAIN(ng)%NorthEast_Test(tile)) THEN
            IF (LBC(inorth,ibISF(i),ng)%acquire) THEN
              ISF_LOBC(i,ng) % isf_north  = IniVal
              ISF_LOBC(i,ng) % isfG_north = IniVal
            END IF
          END IF
        END DO

      END IF
#endif

#ifdef AVERAGES
!
!-----------------------------------------------------------------------
!  Initialize time-averaged variables.
!-----------------------------------------------------------------------
!
!  Iceshelf model internal fields.
!
      DO nf=1,nIsfF
        IF (iFisf(nf).gt.0) THEN
          IF (LisfFavg(nf,ng)) THEN
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                ISF_FAVG(nf,ng) % var(i,j) = IniVal
              END DO
            END DO
          END IF
        END IF
      END DO
!
!  Time-averaged state fields.
!
      DO ns=1,nIsfS
        IF (iSisf(ns).gt.0) THEN
          IF (LisfSavg(ns,ng)) THEN
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                ISF_SAVG(ns,ng) % var(i,j) = IniVal
              END DO
            END DO
          END IF
        END IF
      END DO
!
!  Time-averaged tracer fields.
!
      IF (LisfVavg(ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            ISF_VAVG(ng) % var(i,j) = IniVal
          END DO
        END DO
      END IF
!
      IF (Lisftrc(ng)) THEN
        DO itrc=1,NT(ng)
          IF (LisfTavg(itrc,ng)) THEN
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                ISF_TAVG(itrc,ng) % var(i,j) = IniVal
              END DO
            END DO
          END IF
          IF (LisfAavg(itrc,ng)) THEN
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                ISF_AAVG(itrc,ng) % var(i,j) = IniVal
              END DO
            END DO
          END IF
        END DO
      END IF
#endif
!
      RETURN
      END SUBROUTINE initialize_isf
!
      END MODULE mod_isf
