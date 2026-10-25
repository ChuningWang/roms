#include "cppdefs.h"
      MODULE isf_melt_prescribed_mod
#ifdef ICESHEET_PRESCRIBED_MELT
!
!git $Id$
!====================================================== Chuning Wang ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  Prescribed (analytical) basal melt profiles for the MISMIP+         !
!  ice-only experiments (Asay-Davis et al. 2016, Sect. 2.1/2.2;       !
!  protocol freeze table .agents/docs/plans/                           !
!  mismip_plus_protocol_freeze.md, section 5).                         !
!                                                                      !
!  In no-thermodynamics builds there is no 3EQN source writing        !
!  Fi(:,:,ifImel); this module fills it from the protocol's ad hoc    !
!  melt parameterizations so isf_advection consumes a prescribed       !
!  field, exactly as it consumes the instantaneous 3EQN rate under    !
!  ICESHELF_THERMO.  The experiment is selected by the input flag     !
!  mismip_melt (isf.in):                                               !
!                                                                      !
!    0        none (spin-up, Ice0, Ice1ra/Ice2ra re-advance)           !
!    1  Ice1  m = 0.2/a * tanh(Hc/75 m) * max(-100 m - zd, 0)          !
!              (cavity-thickness taper + linear thermal driving;       !
!              zero at the grounding line, ~75-100 m/a at deep draft)  !
!    2  Ice2  m = 100 m/a where x > 480 km (calving-event analogue)   !
!                                                                      !
!  Floating ice only (rmask_isf = 1); grounded and ice-free cells      !
!  get zero.  Rates are ice-equivalent [m/s] (M-9 convention),         !
!  converted with the MISMIP+ year of 31 556 926 s.  The 100-a         !
!  application window of Ice1r/Ice2r is realized by run boundaries     !
!  (the re-advance experiments restart with mismip_melt = 0), not by   !
!  internal time logic.                                                !
!                                                                      !
!  Geometry: zd = zice (ice draft, negative; for floating ice          !
!  zice = -(rho_i/rho0)*H from isf_masks), Hc = zd - zb is the        !
!  ocean-cavity thickness Hc = zd - zb = zice - (-h) = zice + h.       !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
      PRIVATE
      PUBLIC :: isf_melt_prescribed
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE isf_melt_prescribed (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_grid
      USE mod_isf,      ONLY : ISF, ifImel
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 42, __LINE__, MyFile)
# endif
      CALL isf_melt_prescribed_tile (ng, tile, model,                   &
     &                               LBi, UBi, LBj, UBj,               &
     &                               IminS, ImaxS, JminS, JmaxS,       &
     &                               GRID(ng) % h,                      &
     &                               GRID(ng) % zice,                   &
     &                               GRID(ng) % xr,                     &
     &                               GRID(ng) % rmask_isf,              &
     &                               ISF(ng) % Fi(:,:,ifImel))
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE isf_melt_prescribed
!
!***********************************************************************
      SUBROUTINE isf_melt_prescribed_tile (ng, tile, model,             &
     &                                    LBi, UBi, LBj, UBj,          &
     &                                    IminS, ImaxS, JminS, JmaxS,  &
     &                                    h, zice, xr, rmask_isf,      &
     &                                    melt)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_isf,      ONLY : mismip_melt
      USE mod_iounits,  ONLY : stdout
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)    :: h(LBi:,LBj:)
      real(r8), intent(in)    :: zice(LBi:,LBj:)
      real(r8), intent(in)    :: xr(LBi:,LBj:)
      real(r8), intent(in)    :: rmask_isf(LBi:,LBj:)
      real(r8), intent(inout) :: melt(LBi:,LBj:)
# else
      real(r8), intent(in)    :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: zice(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: rmask_isf(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: melt(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8) :: cff, h_cav, sec_year
!
!  MISMIP+ protocol constants (freeze table section 5).
!
      real(r8), parameter :: gamma_melt = 0.2_r8        ! [1/a]
      real(r8), parameter :: Hc0        = 75.0_r8       ! [m]
      real(r8), parameter :: z0         = -100.0_r8     ! [m]
      real(r8), parameter :: m_ice2     = 100.0_r8      ! [m/a]
      real(r8), parameter :: x_ice2     = 480.0E3_r8    ! [m]
!
# include "set_bounds.h"
!
      sec_year=31556926.0_r8
!
!-----------------------------------------------------------------------
!  Fill the prescribed melt rate on floating ice; zero elsewhere.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          IF (rmask_isf(i,j).eq.1.0_r8) THEN
            IF (mismip_melt(ng).eq.1) THEN
!
!  Ice1 profile: cavity-thickness taper times linear thermal driving.
!
              h_cav=zice(i,j)+h(i,j)
              cff=gamma_melt*TANH(h_cav/Hc0)*MAX(z0-zice(i,j), 0.0_r8)
              melt(i,j)=cff/sec_year
            ELSE IF (mismip_melt(ng).eq.2) THEN
!
!  Ice2 profile: uniform ablation downstream of x = 480 km.
!
              IF (xr(i,j).gt.x_ice2) THEN
                melt(i,j)=m_ice2/sec_year
              ELSE
                melt(i,j)=0.0_r8
              END IF
            ELSE
              melt(i,j)=0.0_r8
            END IF
          ELSE
            melt(i,j)=0.0_r8
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE isf_melt_prescribed_tile
#endif
      END MODULE isf_melt_prescribed_mod
