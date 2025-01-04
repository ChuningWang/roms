      MODULE mod_isf
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2024 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!====================================== === Benjamin K. Galton-Fenzi ===
!                                                                      !
!  Ice Shelf Model Kernel Variables:                                   !
!  This module declares Ice Shelf Model Kernel constants and structure !
!  containing several variables describing ice draft and thickness     !
!  properties.                                                         !
!                                                                      !
!  Variables:                                                          !
#ifdef ICESHELF_MOM
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
!    isItis    iceshelf temperature [degC]                             !
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
# ifdef ICESHELF_TRACER
!  Iceshelf tracer model variables:                                    !
!                                                                      !
!    icIvol    iceshelf freeezed volume [m]                            !
!    icIvar    tracer concentration in refreezed volume [TU]           !
!                                                                      !
# endif
#endif
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
      integer :: idIhis       ! Iceshelf thickness
      integer :: idIuis       ! Iceshelf u-velocity
      integer :: idIvis       ! Iceshelf v-velocity
      integer :: idIuer       ! iceshelf eastward velocity
      integer :: idIvnr       ! iceshelf northward velocity
      integer :: idItis       ! Iceshelf temperature
      integer :: idIeis       ! Iceshelf enthalpy
      integer :: idImel       ! Iceshelf melt rate
      integer :: idIgat       ! gammaT
      integer :: idIgas       ! gammaS
      integer :: idItst       ! T*
      integer :: idIust       ! U*
      integer :: idItbo       ! Boundary layer temperature
      integer :: idIsbo       ! Boundary layer salinity
#ifdef ICESHELF_TRACER
      integer  :: idIvol                 ! Volume of refreezed water
      integer, allocatable :: idItrc(:)  ! Tracer variables
#endif
!
!  Ice model state prognostic variables indices.
!
#ifdef ICESHELF_MOM
      integer, parameter :: nIsfS = 5    ! number of ice state variables
#else
      integer, parameter :: nIsfS = 0    ! number of ice state variables
#endif
!
      integer :: iSisf(nIsfS)            ! state I/O indices
!
      integer, parameter :: isIhis =  1  ! Iceshelf thickness
      integer, parameter :: isIuis =  2  ! Iceshelf u-velocity
      integer, parameter :: isIvis =  3  ! Iceshelf v-velocity
      integer, parameter :: isItis =  4  ! Iceshelf temperature
      integer, parameter :: isIeis =  5  ! Iceshelf enthalpy

#ifdef ICESHELF_MOM
!
!  Iceshelf model state lateral boundary conditions indices.
!
      integer :: ibISF(nIsfS)            ! indices to LBC switch
      integer :: isfOBC(4,nIsfS)         ! I/O metadata indices
#endif
!
#if defined ICESHELF_3EQN_VBC
      integer, parameter :: nIsfF = 7    ! number of ice field variables
#elif defined ICESHELF_2EQN_VBC
      integer, parameter :: nIsfF = 1    ! number of ice field variables
#else
      integer, parameter :: nIsfF = 0    ! number of ice field variables
#endif
!
      integer :: iFisf(nIsfF)            ! internal fields I/O indices
!
      integer, parameter :: icImel =  1  ! Iceshelf melt rate
      integer, parameter :: icIgat =  2  ! gammaT
      integer, parameter :: icIgas =  3  ! gammaS
      integer, parameter :: icItst =  4  ! T*
      integer, parameter :: icIust =  5  ! U*
      integer, parameter :: icItbo =  6  ! Boundary layer temperature
      integer, parameter :: icIsbo =  7  ! Boundary layer salinity
!
      integer :: isf_ic_brylast  ! Boundary layer salinity
!
!-----------------------------------------------------------------------
!  Iceshelf model parameters.
!                                                                      !
!   a, b, c  - Freezing point EOS constants [degC PSU-1, N/A, degC m-1]!
!   L        - Latent heat of melting [J kg-1]                         !
#if defined ICESHELF_2EQN_VBC
!   gammaTS  - 2-EQ formula turbulent exchange coefficient             !
!   refSalt  - Reference salinity [PSU]                                !
#elif defined ICESHELF_3EQN_VBC
!   Pr, Sc   - Prandtl number and Schmidt number                       !
!   Cd       - Drag coefficient                                        !
!   visc     - Background viscosity                                    !
!   L        - Latent heat of melting [J kg-1]                         !
!   rho_i    - Ice density [kg m-3]                                    !
!   Tis, Sis - Ice temperature/salinity                                !
#endif
!-----------------------------------------------------------------------
!
      integer, allocatable :: lsnew(:)
      integer, allocatable :: lsold(:)
!$OMP THREADPRIVATE (lsnew, lsold)
#ifdef AVERAGES
!
!  Switches to process time-averaged iceshelf model state and internal
!  variables.
!
      logical, allocatable :: LisfFavg(:,:)     ! internal variables
      logical, allocatable :: LisfSavg(:,:)     ! state variables
# ifdef ICESHELF_TRACER
      logical, allocatable :: LisfVavg(:)       ! refreezed volume
      logical, allocatable :: LisfTavg(:,:)     ! refreezed tracer
# endif
#endif
!
!  Gamma constraining factors
!
      real(r8), allocatable :: gamU(:)          ! gamma scaling factor
      real(r8), allocatable :: min_ustar(:)     ! minimum ustar value
#ifdef ICESHELF_TRACER
!
!  Logical switches
!
      logical, allocatable :: LtracerIsf(:,:)   ! tracer source switch
#endif
!
!  Other constant parameters
!
      real(r8), parameter :: min_his = 1.0_r8
      real(r8), parameter :: a       = -0.057_r8
      real(r8), parameter :: b       = 0.0939_r8
      real(r8), parameter :: c       = -7.61e-4_r8
      real(r8), parameter :: L       = 3.33e5_r8
#if defined ICESHELF_2EQN_VBC
      real(r8), parameter :: gammaTS = 0.0001_r8
      real(r8), parameter :: refSalt = 34.4_r8
#elif defined ICESHELF_3EQN_VBC
      real(r8), parameter :: Pr      = 13.8_r8
      real(r8), parameter :: Sc      = 2432.2_r8
      real(r8), parameter :: Cd      = 5.0e-3_r8
      real(r8), parameter :: visc    = 1.95e-6_r8
!     real(r8), parameter :: small   = 1.0e-3_r8 ! Diffusion + tidal
      real(r8), parameter :: small   = 8.0e-6_r8 ! Diffusion only
      real(r8), parameter :: rho_i   = 916.7_r8
      real(r8), parameter :: Tis     = -20.0_r8
      real(r8), parameter :: Sis     = 0.0_r8
#endif
      real(r8), parameter :: eps     = 1.0E-14_r8
      real(r8), parameter :: trcIsf  = 1.0_r8
!
!-----------------------------------------------------------------------
!  Define derived-type structure iceshelf model state and internal
!  arrays.
!-----------------------------------------------------------------------
!
      TYPE T_ISF

        real(r8), pointer :: Fi(:,:,:)               ! [i,j,1:nIsfF]
        real(r8), pointer :: Si(:,:,:,:)             ! [i,j,1:2,1:nIsfS]
#ifdef ICESHELF_TRACER
        real(r8), pointer :: vol(:,:)                ! [i,j]
        real(r8), pointer :: trc(:,:,:)              ! [i,j,1:NPT]
#endif

      END TYPE T_ISF
!
      TYPE (T_ISF), allocatable :: ISF(:)            ! [Ngrids]
#ifdef ICESHELF_MOM
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

        real(r8), pointer :: var(:,:)                ! [i,j]

      END TYPE T_ISF_AVG
!
      TYPE (T_ISF_AVG), allocatable :: ISF_FAVG(:,:) ! [nIsfF,Ngrids]
      TYPE (T_ISF_AVG), allocatable :: ISF_SAVG(:,:) ! [nIsfS,Ngrids]
# ifdef ICESHELF_TRACER
      TYPE (T_ISF_AVG), allocatable :: ISF_VAVG(:)   ! [Ngrids]
      TYPE (T_ISF_AVG), allocatable :: ISF_TAVG(:,:) ! [NPT,Ngrids]
# endif
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
      USE mod_param,   ONLY : Dmem, LBC, Ngrids, NPT
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
      integer :: i, itrc
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
#ifdef ICESHELF_TRACER

        IF (.not.allocated(idItrc))                                     &
     &    allocate ( idItrc(NPT) )
#endif

        IF (.not.allocated(lsnew))                                      &
          allocate ( lsnew(Ngrids) )
        lsnew(1:Ngrids) = 1

        IF (.not.allocated(lsold))                                      &
          allocate ( lsold(Ngrids) )
        lsold(1:Ngrids) = 1

#ifdef AVERAGES

        IF (.not.allocated(LisfFavg))                                   &
          allocate ( LisfFavg(nIsfF,Ngrids) )

        IF (.not.allocated(LisfSavg))                                   &
          allocate ( LisfSavg(nIsfS,Ngrids) )
# ifdef ICESHELF_TRACER

        IF (.not.allocated(LisfVavg))                                   &
          allocate ( LisfVavg(Ngrids) )

        IF (.not.allocated(LisfTavg))                                   &
          allocate ( LisfTavg(NPT,Ngrids) )
# endif
#endif

        IF (.not.allocated(gamU))                                       &
     &    allocate ( gamU(Ngrids) )

        IF (.not.allocated(min_ustar))                                  &
     &    allocate ( min_ustar(Ngrids) )
#ifdef ICESHELF_TRACER

        IF (.not.allocated(LtracerIsf))                                 &
     &    allocate ( LtracerIsf(NPT,Ngrids) )
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
#ifdef ICESHELF_TRACER
!
        allocate ( ISF(ng) % vol(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
!
        allocate ( ISF(ng) % trc(LBi:UBi,LBj:UBj,NPT) )
        Dmem(ng)=Dmem(ng)+size2d
#endif
#ifdef ICESHELF_MOM
!
!-----------------------------------------------------------------------
!  Allocate derived-type structure iceshelf model lateral boundary
!  variables.
!-----------------------------------------------------------------------
!
        IF (ng.eq.1) allocate ( ISF_LOBC(nIsfS,Ngrids) )
!
!  Nonlinear iceself model 'state' lateral boundary arrays.
!
        DO i=1,nIsfS
          IF (LBC(iwest,ibISF(i),ng)%acquire) THEN
            allocate ( ISF_LOBC(i,ng) % isf_west(LBj:UBj) )
            Dmem(ng)=Dmem(ng)+Ysize

            allocate ( ISF_LOBC(i,ng) % isfG_west(LBj:UBj,2) )
            Dmem(ng)=Dmem(ng)+2.0_r8*Ysize
          END IF
!
          IF (LBC(ieast,ibISF(i),ng)%acquire) THEN
            allocate ( ISF_LOBC(i,ng) % isf_east(LBj:UBj) )
            Dmem(ng)=Dmem(ng)+Ysize

            allocate ( ISF_LOBC(i,ng) % isfG_east(LBj:UBj,2) )
            Dmem(ng)=Dmem(ng)+2.0_r8*Ysize
          END IF
!
          IF (LBC(isouth,ibISF(i),ng)%acquire) THEN
            allocate ( ISF_LOBC(i,ng) % isf_south(LBi:UBi) )
            Dmem(ng)=Dmem(ng)+Xsize

            allocate ( ISF_LOBC(i,ng) % isfG_south(LBi:UBi,2) )
            Dmem(ng)=Dmem(ng)+2.0_r8*Xsize
          END IF
!
          IF (LBC(inorth,ibISF(i),ng)%acquire) THEN
            allocate ( ISF_LOBC(i,ng) % isf_north(LBi:UBi) )
            Dmem(ng)=Dmem(ng)+Xsize

            allocate ( ISF_LOBC(i,ng) % isfG_north(LBi:UBi,2) )
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
# ifdef ICESHELF_TRACER
!
        IF (ng.eq.1) THEN
          allocate ( ISF_VAVG(Ngrids) )
          allocate ( ISF_TAVG(NPT,Ngrids) )
        END IF
!
        IF (LisfVavg(ng)) THEN
          allocate ( ISF_VAVG(ng) % var(LBi:UBi,LBj:UBj) )
          Dmem(ng)=Dmem(ng)+size2d
        END IF
!
        DO itrc=1,NPT
          IF (LisfTavg(itrc,ng)) THEN
            allocate ( ISF_TAVG(itrc,ng) % var(LBi:UBi,LBj:UBj) )
            Dmem(ng)=Dmem(ng)+size2d
          END IF
        END DO
# endif
      END IF
#endif
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
      USE mod_param,   ONLY : Ngrids, NPT
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
# ifdef ICESHELF_TRACER
!
      IF (.not.destroy(ng, ISF(ng)%vol, MyFile,                         &
     &                 __LINE__, 'ISF(ng)%vol')) RETURN
!
      IF (.not.destroy(ng, ISF(ng)%trc, MyFile,                         &
     &                 __LINE__, 'ISF(ng)%trc')) RETURN
# endif
# ifdef ICESHELF_MOM
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
#  ifdef ICESHELF_TRACER
!
      IF (LisfVavg(ng)) THEN
        IF (.not.destroy(ng, ISF_VAVG(ng)%var, MyFile,                  &
     &                   __LINE__, 'ISF_VAVG(ng)%var')) RETURN
      END IF
!
      DO itrc=1,NPT
        IF (LisfTavg(itrc,ng)) THEN
          IF (.not.destroy(ng, ISF_TAVG(itrc,ng)%var, MyFile,           &
     &                     __LINE__,'ISF_TAVG(itrc,ng)%var')) RETURN
        END IF
      END DO
#  endif
# endif
#endif
!
!-----------------------------------------------------------------------
!  Deallocate derived-type ISF and ISF_LOBC structures.
!-----------------------------------------------------------------------
!
      IF (allocated(ISF))      deallocate ( ISF )
#ifdef ICESHELF_MOM
      IF (allocated(ISF_LOBC)) deallocate ( ISF_LOBC )
#endif
#ifdef AVERAGES
      IF (allocated(ISF_FAVG)) deallocate ( ISF_FAVG )
      IF (allocated(ISF_SAVG)) deallocate ( ISF_SAVG )
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
      USE mod_param,   ONLY : BOUNDS, DOMAIN, LBC, iNLM, NPT
      USE mod_scalars, ONLY : iwest, ieast, isouth, inorth
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: i, j, nf, ns, itrc
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
        END DO
      END DO
#ifdef ICESHELF_TRACER
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          ISF(ng) % vol(i,j) = IniVal
          DO itrc=1,NPT
            ISF(ng) % trc(i,j,itrc) = IniVal
          END DO
        END DO
      END DO
#endif
#ifdef ICESHELF_MOM
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
# ifdef ICESHELF_TRACER
      IF (LisfVavg(ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            ISF_VAVG(ng) % var(i,j) = IniVal
          END DO
        END DO
      END IF
!
      DO itrc=1,NPT
        IF (LisfTavg(itrc,ng)) THEN
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ISF_TAVG(itrc,ng) % var(i,j) = IniVal
            END DO
          END DO
        END IF
      END DO
# endif
#endif
!
      RETURN
      END SUBROUTINE initialize_isf
!
      END MODULE mod_isf
