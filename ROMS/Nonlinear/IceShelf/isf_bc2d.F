#include "cppdefs.h"

      MODULE isf_bc2d_mod

#ifdef ICESHEET
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2024 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           W. Paul Budgell    !
!    See License_ROMS.md                            Katherine Hedstrom !
!================================================== Hernan G. Arango ===
!                                                       Chuning Wang   !
!  Sets the lateral boundary conditions on the iceshelf thickness      !
!  (RHO points).                                                       !
!                                                                      !
!  Supported boundary condition types (see isf.in):                    !
!     clamped    - set the ghost ring to the boundary data             !
!     gradient   - copy the interior value to the boundary             !
!     closed     - copy the interior value to the boundary             !
!                                                                      !
!  The ice-sheet model does NOT support radiation or nudging options.  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_isf
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
!
      implicit none
!
      PRIVATE
      PUBLIC  :: isf_bc2d
      PUBLIC  :: isf_bc2d_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE isf_bc2d (ng, tile, model)
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
# ifdef PROFILE
      CALL wclock_on (ng, model, 42, __LINE__, MyFile)
# endif
      CALL isf_bc2d_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    ISF(ng) % Si(:,:,lsnew(ng),isIhis))
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE isf_bc2d
!
!***********************************************************************
      SUBROUTINE isf_bc2d_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          hi)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: hi(LBi:,LBj:)
# else
      real(r8), intent(inout) :: hi(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
!
#include "set_bounds.h"
      IF (.not.EWperiodic(ng)) THEN
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          IF (LBC(iwest,ibISF(isIhis),ng)%clamped) THEN
            DO j=Jstr,Jend
              hi(Istr-1,j)=ISF_LOBC(isIhis,ng)%isf_west(j)
# ifdef MASKING
              hi(Istr-1,j)=hi(Istr-1,j)*GRID(ng)%rmask(Istr-1,j)
# endif
            END DO
          ELSE IF (LBC(iwest,ibISF(isIhis),ng)%gradient .or.            &
     &             LBC(iwest,ibISF(isIhis),ng)%closed) THEN
            DO j=Jstr,Jend
              hi(Istr-1,j)=hi(Istr,j)
# ifdef MASKING
              hi(Istr-1,j)=hi(Istr-1,j)*GRID(ng)%rmask(Istr-1,j)
# endif
            END DO
          END IF
        END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          IF (LBC(ieast,ibISF(isIhis),ng)%clamped) THEN
            DO j=Jstr,Jend
              hi(Iend+1,j)=ISF_LOBC(isIhis,ng)%isf_east(j)
# ifdef MASKING
              hi(Iend+1,j)=hi(Iend+1,j)*GRID(ng)%rmask(Iend+1,j)
# endif
            END DO
          ELSE IF (LBC(ieast,ibISF(isIhis),ng)%gradient .or.            &
     &             LBC(ieast,ibISF(isIhis),ng)%closed) THEN
            DO j=Jstr,Jend
              hi(Iend+1,j)=hi(Iend,j)
# ifdef MASKING
              hi(Iend+1,j)=hi(Iend+1,j)*GRID(ng)%rmask(Iend+1,j)
# endif
            END DO
          END IF
        END IF
      END IF

      IF (.not.NSperiodic(ng)) THEN
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          IF (LBC(isouth,ibISF(isIhis),ng)%clamped) THEN
            DO i=Istr,Iend
              hi(i,Jstr-1)=ISF_LOBC(isIhis,ng)%isf_south(i)
# ifdef MASKING
              hi(i,Jstr-1)=hi(i,Jstr-1)*GRID(ng)%rmask(i,Jstr-1)
# endif
            END DO
          ELSE IF (LBC(isouth,ibISF(isIhis),ng)%gradient .or.           &
     &             LBC(isouth,ibISF(isIhis),ng)%closed) THEN
            DO i=Istr,Iend
              hi(i,Jstr-1)=hi(i,Jstr)
# ifdef MASKING
              hi(i,Jstr-1)=hi(i,Jstr-1)*GRID(ng)%rmask(i,Jstr-1)
# endif
            END DO
          END IF
        END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          IF (LBC(inorth,ibISF(isIhis),ng)%clamped) THEN
            DO i=Istr,Iend
              hi(i,Jend+1)=ISF_LOBC(isIhis,ng)%isf_north(i)
# ifdef MASKING
              hi(i,Jend+1)=hi(i,Jend+1)*GRID(ng)%rmask(i,Jend+1)
# endif
            END DO
          ELSE IF (LBC(inorth,ibISF(isIhis),ng)%gradient .or.           &
     &             LBC(inorth,ibISF(isIhis),ng)%closed) THEN
            DO i=Istr,Iend
              hi(i,Jend+1)=hi(i,Jend)
# ifdef MASKING
              hi(i,Jend+1)=hi(i,Jend+1)*GRID(ng)%rmask(i,Jend+1)
# endif
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          hi(Istr-1,Jstr-1)=0.5_r8*(hi(Istr,Jstr-1)+hi(Istr-1,Jstr))
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          hi(Iend+1,Jstr-1)=0.5_r8*(hi(Iend+1,Jstr)+hi(Iend,Jstr-1))
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          hi(Istr-1,Jend+1)=0.5_r8*(hi(Istr-1,Jend)+hi(Istr,Jend+1))
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          hi(Iend+1,Jend+1)=0.5_r8*(hi(Iend+1,Jend)+hi(Iend,Jend+1))
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Periodic boundary conditions: exchange the RHO-point halos.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          hi)
      END IF
!
      RETURN
      END SUBROUTINE isf_bc2d_tile
#endif
      END MODULE isf_bc2d_mod
