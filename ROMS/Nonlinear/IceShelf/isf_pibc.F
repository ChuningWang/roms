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
!     mixed      - clamp on inflow, gradient on outflow                !
!     gradient   - copy the interior value to the boundary             !
!     closed     - zero (or gamma2 free-slip on the S/N edges)         !
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
     &                    lsold(ng), lsnew(ng),                         &
     &                    ISF(ng) % Si(:,:,:,isIuis),                   &
     &                    ISF(ng) % Si(:,:,:,isIvis))
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
     &                          lsold, lsnew,                           &
     &                          ui, vi)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: lsold, lsnew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ui(LBi:,LBj:,:)
      real(r8), intent(inout) :: vi(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: vi(LBi:UBi,LBj:UBj,2)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
!
!  PSI-point boundary ring.  The PSI grid spans [1:Lm+1]x[1:Mm+1]
!  (one point shorter than RHO in each direction); IstrP..IendP and
!  JstrP..JendP are the PSI interior bounds.  The ring immediately
!  outside the interior (IstrP-1/IendP+1, JstrP-1/JendP+1) is the
!  ghost layer read by the SSA element assembly and the thickness
!  advection, so it is filled with the same boundary value.
!
#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        IF (LBC(iwest,ibISF(isIuis),ng)%clamped) THEN
          DO j=JstrP,JendP
            ui(IstrP  ,j,lsnew)=ISF_LOBC(isIuis,ng)%isf_west(j)
            ui(IstrP-1,j,lsnew)=ui(IstrP,j,lsnew)
            vi(IstrP  ,j,lsnew)=ISF_LOBC(isIvis,ng)%isf_west(j)
            vi(IstrP-1,j,lsnew)=vi(IstrP,j,lsnew)
# ifdef MASKING
            ui(IstrP,j,lsnew)=ui(IstrP,j,lsnew)*                        &
     &                        GRID(ng)%pmask(IstrP,j)
            ui(IstrP-1,j,lsnew)=ui(IstrP-1,j,lsnew)*                    &
     &                          GRID(ng)%pmask(IstrP-1,j)
            vi(IstrP,j,lsnew)=vi(IstrP,j,lsnew)*                        &
     &                        GRID(ng)%pmask(IstrP,j)
            vi(IstrP-1,j,lsnew)=vi(IstrP-1,j,lsnew)*                    &
     &                          GRID(ng)%pmask(IstrP-1,j)
# endif
          END DO
!
        ELSE IF (LBC(iwest,ibISF(isIuis),ng)%mixed) THEN
          DO j=JstrP,JendP
            IF (ui(IstrP,j,lsnew).ge.0.0_r8) THEN
              ui(IstrP,j,lsnew)=ISF_LOBC(isIuis,ng)%isf_west(j)
            ELSE
              ui(IstrP,j,lsnew)=ui(IstrP+1,j,lsnew)
            END IF
            IF (vi(IstrP,j,lsnew).ge.0.0_r8) THEN
              vi(IstrP,j,lsnew)=ISF_LOBC(isIvis,ng)%isf_west(j)
            ELSE
              vi(IstrP,j,lsnew)=vi(IstrP+1,j,lsnew)
            END IF
            ui(IstrP-1,j,lsnew)=ui(IstrP,j,lsnew)
            vi(IstrP-1,j,lsnew)=vi(IstrP,j,lsnew)
# ifdef MASKING
            ui(IstrP,j,lsnew)=ui(IstrP,j,lsnew)*                        &
     &                        GRID(ng)%pmask(IstrP,j)
            ui(IstrP-1,j,lsnew)=ui(IstrP-1,j,lsnew)*                    &
     &                          GRID(ng)%pmask(IstrP-1,j)
            vi(IstrP,j,lsnew)=vi(IstrP,j,lsnew)*                        &
     &                        GRID(ng)%pmask(IstrP,j)
            vi(IstrP-1,j,lsnew)=vi(IstrP-1,j,lsnew)*                    &
     &                          GRID(ng)%pmask(IstrP-1,j)
# endif
          END DO
!
        ELSE IF (LBC(iwest,ibISF(isIuis),ng)%gradient) THEN
          DO j=JstrP,JendP
            ui(IstrP,j,lsnew)=ui(IstrP+1,j,lsnew)
            ui(IstrP-1,j,lsnew)=ui(IstrP,j,lsnew)
            vi(IstrP,j,lsnew)=vi(IstrP+1,j,lsnew)
            vi(IstrP-1,j,lsnew)=vi(IstrP,j,lsnew)
# ifdef MASKING
            ui(IstrP,j,lsnew)=ui(IstrP,j,lsnew)*                        &
     &                        GRID(ng)%pmask(IstrP,j)
            ui(IstrP-1,j,lsnew)=ui(IstrP-1,j,lsnew)*                    &
     &                          GRID(ng)%pmask(IstrP-1,j)
            vi(IstrP,j,lsnew)=vi(IstrP,j,lsnew)*                        &
     &                        GRID(ng)%pmask(IstrP,j)
            vi(IstrP-1,j,lsnew)=vi(IstrP-1,j,lsnew)*                    &
     &                          GRID(ng)%pmask(IstrP-1,j)
# endif
          END DO
!
        ELSE IF (LBC(iwest,ibISF(isIuis),ng)%closed) THEN
          DO j=JstrP,JendP
            ui(IstrP  ,j,lsnew)=0.0_r8
            ui(IstrP-1,j,lsnew)=0.0_r8
            vi(IstrP  ,j,lsnew)=0.0_r8
            vi(IstrP-1,j,lsnew)=0.0_r8
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        IF (LBC(ieast,ibISF(isIuis),ng)%clamped) THEN
          DO j=JstrP,JendP
            ui(IendP  ,j,lsnew)=ISF_LOBC(isIuis,ng)%isf_east(j)
            ui(IendP+1,j,lsnew)=ui(IendP,j,lsnew)
            vi(IendP  ,j,lsnew)=ISF_LOBC(isIvis,ng)%isf_east(j)
            vi(IendP+1,j,lsnew)=vi(IendP,j,lsnew)
# ifdef MASKING
            ui(IendP,j,lsnew)=ui(IendP,j,lsnew)*                        &
     &                        GRID(ng)%pmask(IendP,j)
            ui(IendP+1,j,lsnew)=ui(IendP+1,j,lsnew)*                    &
     &                          GRID(ng)%pmask(IendP+1,j)
            vi(IendP,j,lsnew)=vi(IendP,j,lsnew)*                        &
     &                        GRID(ng)%pmask(IendP,j)
            vi(IendP+1,j,lsnew)=vi(IendP+1,j,lsnew)*                    &
     &                          GRID(ng)%pmask(IendP+1,j)
# endif
          END DO
!
        ELSE IF (LBC(ieast,ibISF(isIuis),ng)%mixed) THEN
          DO j=JstrP,JendP
            IF (ui(IendP+1,j,lsnew).le.0.0_r8) THEN
              ui(IendP+1,j,lsnew)=ISF_LOBC(isIuis,ng)%isf_east(j)
            ELSE
              ui(IendP+1,j,lsnew)=ui(IendP,j,lsnew)
            END IF
            IF (vi(IendP+1,j,lsnew).le.0.0_r8) THEN
              vi(IendP+1,j,lsnew)=ISF_LOBC(isIvis,ng)%isf_east(j)
            ELSE
              vi(IendP+1,j,lsnew)=vi(IendP,j,lsnew)
            END IF
            ui(IendP  ,j,lsnew)=ui(IendP+1,j,lsnew)
            vi(IendP  ,j,lsnew)=vi(IendP+1,j,lsnew)
# ifdef MASKING
            ui(IendP,j,lsnew)=ui(IendP,j,lsnew)*                        &
     &                        GRID(ng)%pmask(IendP,j)
            ui(IendP+1,j,lsnew)=ui(IendP+1,j,lsnew)*                    &
     &                          GRID(ng)%pmask(IendP+1,j)
            vi(IendP,j,lsnew)=vi(IendP,j,lsnew)*                        &
     &                        GRID(ng)%pmask(IendP,j)
            vi(IendP+1,j,lsnew)=vi(IendP+1,j,lsnew)*                    &
     &                          GRID(ng)%pmask(IendP+1,j)
# endif
          END DO
!
        ELSE IF (LBC(ieast,ibISF(isIuis),ng)%gradient) THEN
          DO j=JstrP,JendP
            ui(IendP+1,j,lsnew)=ui(IendP,j,lsnew)
            ui(IendP  ,j,lsnew)=ui(IendP+1,j,lsnew)
            vi(IendP+1,j,lsnew)=vi(IendP,j,lsnew)
            vi(IendP  ,j,lsnew)=vi(IendP+1,j,lsnew)
# ifdef MASKING
            ui(IendP,j,lsnew)=ui(IendP,j,lsnew)*                        &
     &                        GRID(ng)%pmask(IendP,j)
            ui(IendP+1,j,lsnew)=ui(IendP+1,j,lsnew)*                    &
     &                          GRID(ng)%pmask(IendP+1,j)
            vi(IendP,j,lsnew)=vi(IendP,j,lsnew)*                        &
     &                        GRID(ng)%pmask(IendP,j)
            vi(IendP+1,j,lsnew)=vi(IendP+1,j,lsnew)*                    &
     &                          GRID(ng)%pmask(IendP+1,j)
# endif
          END DO
!
        ELSE IF (LBC(ieast,ibISF(isIuis),ng)%closed) THEN
          DO j=JstrP,JendP
            ui(IendP  ,j,lsnew)=0.0_r8
            ui(IendP+1,j,lsnew)=0.0_r8
            vi(IendP  ,j,lsnew)=0.0_r8
            vi(IendP+1,j,lsnew)=0.0_r8
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        IF (LBC(isouth,ibISF(isIuis),ng)%clamped) THEN
          DO i=IstrP,IendP
            ui(i,JstrP  ,lsnew)=ISF_LOBC(isIuis,ng)%isf_south(i)
            ui(i,JstrP-1,lsnew)=ui(i,JstrP,lsnew)
            vi(i,JstrP  ,lsnew)=ISF_LOBC(isIvis,ng)%isf_south(i)
            vi(i,JstrP-1,lsnew)=vi(i,JstrP,lsnew)
# ifdef MASKING
            ui(i,JstrP,lsnew)=ui(i,JstrP,lsnew)*                        &
     &                        GRID(ng)%pmask(i,JstrP)
            ui(i,JstrP-1,lsnew)=ui(i,JstrP-1,lsnew)*                    &
     &                          GRID(ng)%pmask(i,JstrP-1)
            vi(i,JstrP,lsnew)=vi(i,JstrP,lsnew)*                        &
     &                        GRID(ng)%pmask(i,JstrP)
            vi(i,JstrP-1,lsnew)=vi(i,JstrP-1,lsnew)*                    &
     &                          GRID(ng)%pmask(i,JstrP-1)
# endif
          END DO
!
        ELSE IF (LBC(isouth,ibISF(isIuis),ng)%mixed) THEN
          DO i=IstrP,IendP
            IF (vi(i,JstrP,lsnew).ge.0.0_r8) THEN
              vi(i,JstrP,lsnew)=ISF_LOBC(isIvis,ng)%isf_south(i)
            ELSE
              vi(i,JstrP,lsnew)=vi(i,JstrP+1,lsnew)
            END IF
            IF (ui(i,JstrP,lsnew).ge.0.0_r8) THEN
              ui(i,JstrP,lsnew)=ISF_LOBC(isIuis,ng)%isf_south(i)
            ELSE
              ui(i,JstrP,lsnew)=ui(i,JstrP+1,lsnew)
            END IF
            ui(i,JstrP-1,lsnew)=ui(i,JstrP,lsnew)
            vi(i,JstrP-1,lsnew)=vi(i,JstrP,lsnew)
# ifdef MASKING
            ui(i,JstrP,lsnew)=ui(i,JstrP,lsnew)*                        &
     &                        GRID(ng)%pmask(i,JstrP)
            ui(i,JstrP-1,lsnew)=ui(i,JstrP-1,lsnew)*                    &
     &                          GRID(ng)%pmask(i,JstrP-1)
            vi(i,JstrP,lsnew)=vi(i,JstrP,lsnew)*                        &
     &                        GRID(ng)%pmask(i,JstrP)
            vi(i,JstrP-1,lsnew)=vi(i,JstrP-1,lsnew)*                    &
     &                          GRID(ng)%pmask(i,JstrP-1)
# endif
          END DO
!
        ELSE IF (LBC(isouth,ibISF(isIuis),ng)%gradient) THEN
          DO i=IstrP,IendP
            ui(i,JstrP,lsnew)=ui(i,JstrP+1,lsnew)
            ui(i,JstrP-1,lsnew)=ui(i,JstrP,lsnew)
            vi(i,JstrP,lsnew)=vi(i,JstrP+1,lsnew)
            vi(i,JstrP-1,lsnew)=vi(i,JstrP,lsnew)
# ifdef MASKING
            ui(i,JstrP,lsnew)=ui(i,JstrP,lsnew)*                        &
     &                        GRID(ng)%pmask(i,JstrP)
            ui(i,JstrP-1,lsnew)=ui(i,JstrP-1,lsnew)*                    &
     &                          GRID(ng)%pmask(i,JstrP-1)
            vi(i,JstrP,lsnew)=vi(i,JstrP,lsnew)*                        &
     &                        GRID(ng)%pmask(i,JstrP)
            vi(i,JstrP-1,lsnew)=vi(i,JstrP-1,lsnew)*                    &
     &                          GRID(ng)%pmask(i,JstrP-1)
# endif
          END DO
!
        ELSE IF (LBC(isouth,ibISF(isIuis),ng)%closed) THEN
          DO i=IstrP,IendP
            ui(i,JstrP  ,lsnew)=0.0_r8
            ui(i,JstrP-1,lsnew)=0.0_r8
            vi(i,JstrP  ,lsnew)=0.0_r8
            vi(i,JstrP-1,lsnew)=0.0_r8
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        IF (LBC(inorth,ibISF(isIuis),ng)%clamped) THEN
          DO i=IstrP,IendP
            ui(i,JendP  ,lsnew)=ISF_LOBC(isIuis,ng)%isf_north(i)
            ui(i,JendP+1,lsnew)=ui(i,JendP,lsnew)
            vi(i,JendP  ,lsnew)=ISF_LOBC(isIvis,ng)%isf_north(i)
            vi(i,JendP+1,lsnew)=vi(i,JendP,lsnew)
# ifdef MASKING
            ui(i,JendP,lsnew)=ui(i,JendP,lsnew)*                        &
     &                        GRID(ng)%pmask(i,JendP)
            ui(i,JendP+1,lsnew)=ui(i,JendP+1,lsnew)*                    &
     &                          GRID(ng)%pmask(i,JendP+1)
            vi(i,JendP,lsnew)=vi(i,JendP,lsnew)*                        &
     &                        GRID(ng)%pmask(i,JendP)
            vi(i,JendP+1,lsnew)=vi(i,JendP+1,lsnew)*                    &
     &                          GRID(ng)%pmask(i,JendP+1)
# endif
          END DO
!
        ELSE IF (LBC(inorth,ibISF(isIuis),ng)%mixed) THEN
          DO i=IstrP,IendP
            IF (vi(i,JendP+1,lsnew).le.0.0_r8) THEN
              vi(i,JendP+1,lsnew)=ISF_LOBC(isIvis,ng)%isf_north(i)
            ELSE
              vi(i,JendP+1,lsnew)=vi(i,JendP,lsnew)
            END IF
            IF (ui(i,JendP+1,lsnew).le.0.0_r8) THEN
              ui(i,JendP+1,lsnew)=ISF_LOBC(isIuis,ng)%isf_north(i)
            ELSE
              ui(i,JendP+1,lsnew)=ui(i,JendP,lsnew)
            END IF
            ui(i,JendP  ,lsnew)=ui(i,JendP+1,lsnew)
            vi(i,JendP  ,lsnew)=vi(i,JendP+1,lsnew)
# ifdef MASKING
            ui(i,JendP,lsnew)=ui(i,JendP,lsnew)*                        &
     &                        GRID(ng)%pmask(i,JendP)
            ui(i,JendP+1,lsnew)=ui(i,JendP+1,lsnew)*                    &
     &                          GRID(ng)%pmask(i,JendP+1)
            vi(i,JendP,lsnew)=vi(i,JendP,lsnew)*                        &
     &                        GRID(ng)%pmask(i,JendP)
            vi(i,JendP+1,lsnew)=vi(i,JendP+1,lsnew)*                    &
     &                          GRID(ng)%pmask(i,JendP+1)
# endif
          END DO
!
        ELSE IF (LBC(inorth,ibISF(isIuis),ng)%gradient) THEN
          DO i=IstrP,IendP
            ui(i,JendP+1,lsnew)=ui(i,JendP,lsnew)
            ui(i,JendP  ,lsnew)=ui(i,JendP+1,lsnew)
            vi(i,JendP+1,lsnew)=vi(i,JendP,lsnew)
            vi(i,JendP  ,lsnew)=vi(i,JendP+1,lsnew)
# ifdef MASKING
            ui(i,JendP,lsnew)=ui(i,JendP,lsnew)*                        &
     &                        GRID(ng)%pmask(i,JendP)
            ui(i,JendP+1,lsnew)=ui(i,JendP+1,lsnew)*                    &
     &                          GRID(ng)%pmask(i,JendP+1)
            vi(i,JendP,lsnew)=vi(i,JendP,lsnew)*                        &
     &                        GRID(ng)%pmask(i,JendP)
            vi(i,JendP+1,lsnew)=vi(i,JendP+1,lsnew)*                    &
     &                          GRID(ng)%pmask(i,JendP+1)
# endif
          END DO
!
        ELSE IF (LBC(inorth,ibISF(isIuis),ng)%closed) THEN
          DO i=IstrP,IendP
            ui(i,JendP  ,lsnew)=0.0_r8
            ui(i,JendP+1,lsnew)=0.0_r8
            vi(i,JendP  ,lsnew)=0.0_r8
            vi(i,JendP+1,lsnew)=0.0_r8
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners (only where the SSA assembly reads them).
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          ui(IstrP,JstrP-1,lsnew)=0.5_r8*                                &
     &                            (ui(IstrP  ,JstrP,lsnew)+              &
     &                             ui(IstrP+1,JstrP-1,lsnew))
          vi(IstrP,JstrP-1,lsnew)=0.5_r8*                                &
     &                            (vi(IstrP  ,JstrP,lsnew)+              &
     &                             vi(IstrP+1,JstrP-1,lsnew))
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          ui(IendP+1,JstrP-1,lsnew)=0.5_r8*                              &
     &                              (ui(IendP+1,JstrP,lsnew)+            &
     &                               ui(IendP  ,JstrP-1,lsnew))
          vi(IendP+1,JstrP-1,lsnew)=0.5_r8*                              &
     &                              (vi(IendP+1,JstrP,lsnew)+            &
     &                               vi(IendP  ,JstrP-1,lsnew))
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          ui(IstrP,JendP+1,lsnew)=0.5_r8*                                &
     &                            (ui(IstrP  ,JendP,lsnew)+              &
     &                             ui(IstrP+1,JendP+1,lsnew))
          vi(IstrP,JendP+1,lsnew)=0.5_r8*                                &
     &                            (vi(IstrP  ,JendP,lsnew)+              &
     &                             vi(IstrP+1,JendP+1,lsnew))
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          ui(IendP+1,JendP+1,lsnew)=0.5_r8*                              &
     &                              (ui(IendP+1,JendP,lsnew)+            &
     &                               ui(IendP  ,JendP+1,lsnew))
          vi(IendP+1,JendP+1,lsnew)=0.5_r8*                              &
     &                              (vi(IendP+1,JendP,lsnew)+            &
     &                               vi(IendP  ,JendP+1,lsnew))
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Periodic boundary conditions: exchange the PSI-point halos.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile, LBi, UBi, LBj, UBj,            &
     &                          ui(:,:,lsnew))
        CALL exchange_p2d_tile (ng, tile, LBi, UBi, LBj, UBj,            &
     &                          vi(:,:,lsnew))
      END IF
!
      RETURN
      END SUBROUTINE isf_pibc_tile
#endif
      END MODULE isf_pibc_mod
