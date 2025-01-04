#include "cppdefs.h"
      MODULE iceshelf_mod
#ifdef ICESHELF_DYNAM
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2024 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           W. Paul Budgell    !
!    See License_ROMS.md                            Katherine Hedstrom !
!================================================== Hernan G. Arango ===
!                                                       Chuning Wang   !
!  This module is the primary driver routine for the ROMS nonlinear    !
!  iceshelf model.                                                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
      USE mod_stepping
      USE mod_isf
      USE mod_forces
!
# ifdef ICESHELF_THERMO
      USE isf_vbc_mod,      ONLY : isf_vbc
# endif
# ifdef ICESHELF_MOM
      USE isf_mom_mod,      ONLY : isf_mom
# endif
      USE strings_mod,      ONLY : FoundError, uppercase
!
      implicit none
!
      PUBLIC  :: iceshelf
      PRIVATE
!
      CONTAINS
!
      SUBROUTINE iceshelf (model)
!
!=======================================================================
!                                                                      !
!  Iceshelf model timestepping driver.                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: model
!
!  Local variable declarations.
!
      integer :: tile
      integer :: ig, ng, nl
!
      real(r8), parameter :: dt_large = 1.0E+23_r8
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!=======================================================================
!  Timestep iceshelf model equations.
!=======================================================================
!
!  In nesting applications, the number of nesting layers (NestLayers) is
!  used to facilitate refinement grids and composite/refinament grids
!  combinations. Otherwise, the solution it is looped once for a single
!  grid application (NestLayers = 1).
!
      NEST_LAYER : DO nl=1,NestLayers
# ifdef ICESHELF_THERMO
!
!-----------------------------------------------------------------------
!  Compute the ice-ocean thermodynamics.
!-----------------------------------------------------------------------
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL isf_vbc (ng, tile, model)
          END DO
!$OMP BARRIER
        END DO
# endif

# ifdef ICESHELF_MOM
!
!-----------------------------------------------------------------------
!  Compute the iceshelf momemtum.
!-----------------------------------------------------------------------
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          lsold(ng)=lsnew(ng)
          lsnew(ng)=3-lsold(ng)
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL isf_mom (ng, tile, model)
          END DO
!$OMP BARRIER
        END DO
# endif

      END DO NEST_LAYER
!
      RETURN
      END SUBROUTINE iceshelf
!
#endif
      END MODULE iceshelf_mod
