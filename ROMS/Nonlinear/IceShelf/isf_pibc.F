#include "cppdefs.h"

      MODULE isf_pibc_mod

#ifdef ICESHEET
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2024 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           W. Paul Budgell    !
!    See License_ROMS.md                            Katherine Hedstrom !
!================================================== Hernan G. Arango ===
!                                                       Chuning Wang   !
!  Sets the lateral boundary conditions on the icesheet U- and         !
!  V-velocities.  Both components are stored at PSI points (p2dvar),   !
!  so they share the same boundary ring and are updated together.      !
!                                                                      !
!  Supported boundary condition types (see isf.in):                    !
!     clamped    - set the boundary ring to the boundary data          !
!     gradient   - copy the interior value to the boundary             !
!     closed     - zero                                                !
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
      USE exchange_2d_mod, ONLY : exchange_p2d_tile
!
      implicit none
!
      PRIVATE
      PUBLIC  :: isf_pibc
      PUBLIC  :: isf_pibc_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE isf_pibc (ng, tile, model)
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
      CALL isf_pibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    ISF(ng) % Si(:,:,lsnew(ng),isIuis),           &
     &                    ISF(ng) % Si(:,:,lsnew(ng),isIvis))
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE isf_pibc
!
!***********************************************************************
      SUBROUTINE isf_pibc_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          ui, vi)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ui(LBi:,LBj:)
      real(r8), intent(inout) :: vi(LBi:,LBj:)
# else
      real(r8), intent(inout) :: ui(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vi(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
!
!  PSI-point boundary ring.  The PSI grid spans [1:Lm+1]x[1:Mm+1]:
!  IstrP..IendP and JstrP..JendP are the PSI tile bounds including the
!  boundary nodes themselves (at the western/southern physical edge
!  IstrP=Istr, JstrP=Jstr; at the eastern/northern edge IendP=Iend+1,
!  JendP=Jend+1).  Each physical edge writes its boundary node only
!  (IstrP/IendP column, JstrP/JendP row) over the full PSI extent in
!  the tangential direction; the four corner nodes are written by both
!  adjacent edges with the south/north edges applied last.  The ghost
!  ring outside the boundary node (IstrP-1/IendP+1, ...) is NOT filled:
!  the SSA solver exchanges its own halos and owns every node it reads.
!
#include "set_bounds.h"
      IF (.not.EWperiodic(ng)) THEN
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
          IF (LBC(iwest,ibISF(isIuis),ng)%clamped) THEN
            DO j=JstrP,JendP
              ui(IstrP,j)=ISF_LOBC(isIuis,ng)%isf_west(j)
# ifdef MASKING
              ui(IstrP,j)=ui(IstrP,j)*GRID(ng)%pumask_ist(IstrP,j)
# endif
            END DO
          ELSE IF (LBC(iwest,ibISF(isIuis),ng)%gradient) THEN
            DO j=JstrP,JendP
              ui(IstrP,j)=ui(IstrP+1,j)
# ifdef MASKING
              ui(IstrP,j)=ui(IstrP,j)*GRID(ng)%pumask_ist(IstrP,j)
# endif
            END DO
          ELSE IF (LBC(iwest,ibISF(isIuis),ng)%closed) THEN
            DO j=JstrP,JendP
              ui(IstrP,j)=0.0_r8
            END DO
          END IF
!
          IF (LBC(iwest,ibISF(isIvis),ng)%clamped) THEN
            DO j=JstrP,JendP
              vi(IstrP,j)=ISF_LOBC(isIvis,ng)%isf_west(j)
# ifdef MASKING
              vi(IstrP,j)=vi(IstrP,j)*GRID(ng)%pvmask_ist(IstrP,j)
# endif
            END DO
          ELSE IF (LBC(iwest,ibISF(isIvis),ng)%gradient) THEN
            DO j=JstrP,JendP
              vi(IstrP,j)=vi(IstrP+1,j)
# ifdef MASKING
              vi(IstrP,j)=vi(IstrP,j)*GRID(ng)%pvmask_ist(IstrP,j)
# endif
            END DO
          ELSE IF (LBC(iwest,ibISF(isIvis),ng)%closed) THEN
            DO j=JstrP,JendP
              vi(IstrP,j)=0.0_r8
            END DO
          END IF
        END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
!
          IF (LBC(ieast,ibISF(isIuis),ng)%clamped) THEN
            DO j=JstrP,JendP
              ui(IendP,j)=ISF_LOBC(isIuis,ng)%isf_east(j)
# ifdef MASKING
              ui(IendP,j)=ui(IendP,j)*GRID(ng)%pumask_ist(IendP,j)
# endif
            END DO
          ELSE IF (LBC(ieast,ibISF(isIuis),ng)%gradient) THEN
            DO j=JstrP,JendP
              ui(IendP,j)=ui(IendP-1,j)
# ifdef MASKING
              ui(IendP,j)=ui(IendP,j)*GRID(ng)%pumask_ist(IendP,j)
# endif
            END DO
          ELSE IF (LBC(ieast,ibISF(isIuis),ng)%closed) THEN
            DO j=JstrP,JendP
              ui(IendP,j)=0.0_r8
            END DO
          END IF
!
          IF (LBC(ieast,ibISF(isIvis),ng)%clamped) THEN
            DO j=JstrP,JendP
              vi(IendP,j)=ISF_LOBC(isIvis,ng)%isf_east(j)
# ifdef MASKING
              vi(IendP,j)=vi(IendP,j)*GRID(ng)%pvmask_ist(IendP,j)
# endif
            END DO
          ELSE IF (LBC(ieast,ibISF(isIvis),ng)%gradient) THEN
            DO j=JstrP,JendP
              vi(IendP,j)=vi(IendP-1,j)
# ifdef MASKING
              vi(IendP,j)=vi(IendP,j)*GRID(ng)%pvmask_ist(IendP,j)
# endif
            END DO
          ELSE IF (LBC(ieast,ibISF(isIvis),ng)%closed) THEN
            DO j=JstrP,JendP
              vi(IendP,j)=0.0_r8
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
!
          IF (LBC(isouth,ibISF(isIuis),ng)%clamped) THEN
            DO i=IstrP,IendP
              ui(i,JstrP)=ISF_LOBC(isIuis,ng)%isf_south(i)
# ifdef MASKING
              ui(i,JstrP)=ui(i,JstrP)*GRID(ng)%pumask_ist(i,JstrP)
# endif
            END DO
          ELSE IF (LBC(isouth,ibISF(isIuis),ng)%gradient) THEN
            DO i=IstrP,IendP
              ui(i,JstrP)=ui(i,JstrP+1)
# ifdef MASKING
              ui(i,JstrP)=ui(i,JstrP)*GRID(ng)%pumask_ist(i,JstrP)
# endif
            END DO
          ELSE IF (LBC(isouth,ibISF(isIuis),ng)%closed) THEN
            DO i=IstrP,IendP
              ui(i,JstrP)=0.0_r8
            END DO
          END IF
!
          IF (LBC(isouth,ibISF(isIvis),ng)%clamped) THEN
            DO i=IstrP,IendP
              vi(i,JstrP)=ISF_LOBC(isIvis,ng)%isf_south(i)
# ifdef MASKING
              vi(i,JstrP)=vi(i,JstrP)*GRID(ng)%pvmask_ist(i,JstrP)
# endif
            END DO
          ELSE IF (LBC(isouth,ibISF(isIvis),ng)%gradient) THEN
            DO i=IstrP,IendP
              vi(i,JstrP)=vi(i,JstrP+1)
# ifdef MASKING
              vi(i,JstrP)=vi(i,JstrP)*GRID(ng)%pvmask_ist(i,JstrP)
# endif
            END DO
          ELSE IF (LBC(isouth,ibISF(isIvis),ng)%closed) THEN
            DO i=IstrP,IendP
              vi(i,JstrP)=0.0_r8
            END DO
          END IF
        END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
!
          IF (LBC(inorth,ibISF(isIuis),ng)%clamped) THEN
            DO i=IstrP,IendP
              ui(i,JendP)=ISF_LOBC(isIuis,ng)%isf_north(i)
# ifdef MASKING
              ui(i,JendP)=ui(i,JendP)*GRID(ng)%pumask_ist(i,JendP)
# endif
            END DO
          ELSE IF (LBC(inorth,ibISF(isIuis),ng)%gradient) THEN
            DO i=IstrP,IendP
              ui(i,JendP)=ui(i,JendP-1)
# ifdef MASKING
              ui(i,JendP)=ui(i,JendP)*GRID(ng)%pumask_ist(i,JendP)
# endif
            END DO
          ELSE IF (LBC(inorth,ibISF(isIuis),ng)%closed) THEN
            DO i=IstrP,IendP
              ui(i,JendP)=0.0_r8
            END DO
          END IF
!
          IF (LBC(inorth,ibISF(isIvis),ng)%clamped) THEN
            DO i=IstrP,IendP
              vi(i,JendP)=ISF_LOBC(isIvis,ng)%isf_north(i)
# ifdef MASKING
              vi(i,JendP)=vi(i,JendP)*GRID(ng)%pvmask_ist(i,JendP)
# endif
            END DO
          ELSE IF (LBC(inorth,ibISF(isIvis),ng)%gradient) THEN
            DO i=IstrP,IendP
              vi(i,JendP)=vi(i,JendP-1)
# ifdef MASKING
              vi(i,JendP)=vi(i,JendP)*GRID(ng)%pvmask_ist(i,JendP)
# endif
            END DO
          ELSE IF (LBC(inorth,ibISF(isIvis),ng)%closed) THEN
            DO i=IstrP,IendP
              vi(i,JendP)=0.0_r8
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Periodic boundary conditions: exchange the PSI-point halos.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ui)
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vi)
      END IF
!
      RETURN
      END SUBROUTINE isf_pibc_tile
#endif
      END MODULE isf_pibc_mod
