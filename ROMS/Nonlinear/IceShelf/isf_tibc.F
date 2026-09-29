#include "cppdefs.h"

      MODULE isf_tibc_mod

#ifdef ICESHEET
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2024 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           W. Paul Budgell    !
!    See License_ROMS.md                            Katherine Hedstrom !
!================================================== Hernan G. Arango ===
!                                                       Chuning Wang   !
!  Sets the lateral boundary conditions on the internal iceshelf       !
!  enthalpy (Si(:,:,:,isIeis)).                                        !
!                                                                      !
!  The isIeis boundary file provides the ice TEMPERATURE, not the      !
!  enthalpy; the enthalpy boundary value is computed as                !
!                                                                      !
!     E = h * T                                                        !
!                                                                      !
!  using the thickness boundary value (isIhis) and the temperature     !
!  boundary value (isIeis).  The boundary temperature is also stored   !
!  in the diagnostic field Fi(:,:,ifItis) so that downstream physics   !
!  (e.g. the Glen rate factor) sees a consistent ice temperature.      !
!                                                                      !
!  Supported boundary condition types (see isf.in):                    !
!     clamped    - E = h_bnd * T_bnd                                   !
!     gradient   - copy the interior value to the boundary             !
!     closed     - copy the interior value to the boundary             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_isf
      USE mod_scalars
!
      implicit none
!
      PRIVATE
      PUBLIC  :: isf_tibc
      PUBLIC  :: isf_tibc_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE isf_tibc (ng, tile, model)
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
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 42, __LINE__, MyFile)
# endif
!  Apply to BOTH time levels (same rationale as isf_bc2d):  each
!  layer's enthalpy ghost is constructed from that layer's thickness
!  ghost (E = h*T), keeping the carrier level's boundary state
!  coherent with its interior between solves.
!
      CALL isf_tibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    ISF(ng) % Si(:,:,lsnew(ng),isIhis),           &
     &                    ISF(ng) % Si(:,:,lsnew(ng),isIeis),           &
     &                    ISF(ng) % Fi(:,:,ifItis))
      CALL isf_tibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    ISF(ng) % Si(:,:,lsold(ng),isIhis),           &
     &                    ISF(ng) % Si(:,:,lsold(ng),isIeis),           &
     &                    ISF(ng) % Fi(:,:,ifItis))
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE isf_tibc
!
!***********************************************************************
      SUBROUTINE isf_tibc_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          hi, ei, ti)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)    :: hi(LBi:,LBj:)
      real(r8), intent(inout) :: ei(LBi:,LBj:)
      real(r8), intent(inout) :: ti(LBi:,LBj:)
# else
      real(r8), intent(in)    :: hi(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ei(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ti(LBi:UBi,LBj:UBj)
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
          IF (LBC(iwest,ibISF(isIeis),ng)%clamped) THEN
            DO j=Jstr,Jend
              ei(Istr-1,j)=ISF_LOBC(isIhis,ng)%isf_west(j)*             &
     &                     ISF_LOBC(isIeis,ng)%isf_west(j)
              IF (hi(Istr-1,j).le.min_his(ng)) THEN
                ei(Istr-1,j)=0.0_r8
                ti(Istr-1,j)=0.0_r8
              ELSE
                ti(Istr-1,j)=ISF_LOBC(isIeis,ng)%isf_west(j)
              END IF
# ifdef MASKING
              ei(Istr-1,j)=ei(Istr-1,j)*GRID(ng)%rmask(Istr-1,j)
              ti(Istr-1,j)=ti(Istr-1,j)*GRID(ng)%rmask(Istr-1,j)
# endif
            END DO
          ELSE IF (LBC(iwest,ibISF(isIeis),ng)%gradient .or.            &
     &             LBC(iwest,ibISF(isIeis),ng)%closed) THEN
            DO j=Jstr,Jend
              ei(Istr-1,j)=hi(Istr,j)*ti(Istr,j)
              IF (hi(Istr-1,j).le.min_his(ng)) THEN
                ei(Istr-1,j)=0.0_r8
                ti(Istr-1,j)=0.0_r8
              ELSE
                ti(Istr-1,j)=ti(Istr,j)
              END IF
# ifdef MASKING
              ei(Istr-1,j)=ei(Istr-1,j)*GRID(ng)%rmask(Istr-1,j)
              ti(Istr-1,j)=ti(Istr-1,j)*GRID(ng)%rmask(Istr-1,j)
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
          IF (LBC(ieast,ibISF(isIeis),ng)%clamped) THEN
            DO j=Jstr,Jend
              ei(Iend+1,j)=ISF_LOBC(isIhis,ng)%isf_east(j)*             &
     &                     ISF_LOBC(isIeis,ng)%isf_east(j)
              IF (hi(Iend+1,j).le.min_his(ng)) THEN
                ei(Iend+1,j)=0.0_r8
                ti(Iend+1,j)=0.0_r8
              ELSE
                ti(Iend+1,j)=ISF_LOBC(isIeis,ng)%isf_east(j)
              END IF
# ifdef MASKING
              ei(Iend+1,j)=ei(Iend+1,j)*GRID(ng)%rmask(Iend+1,j)
              ti(Iend+1,j)=ti(Iend+1,j)*GRID(ng)%rmask(Iend+1,j)
# endif
            END DO
          ELSE IF (LBC(ieast,ibISF(isIeis),ng)%gradient .or.            &
     &             LBC(ieast,ibISF(isIeis),ng)%closed) THEN
            DO j=Jstr,Jend
              ei(Iend+1,j)=hi(Iend,j)*ti(Iend,j)
              IF (hi(Iend+1,j).le.min_his(ng)) THEN
                ei(Iend+1,j)=0.0_r8
                ti(Iend+1,j)=0.0_r8
              ELSE
                ti(Iend+1,j)=ti(Iend,j)
              END IF
# ifdef MASKING
              ei(Iend+1,j)=ei(Iend+1,j)*GRID(ng)%rmask(Iend+1,j)
              ti(Iend+1,j)=ti(Iend+1,j)*GRID(ng)%rmask(Iend+1,j)
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
          IF (LBC(isouth,ibISF(isIeis),ng)%clamped) THEN
            DO i=Istr,Iend
              ei(i,Jstr-1)=ISF_LOBC(isIhis,ng)%isf_south(i)*            &
     &                     ISF_LOBC(isIeis,ng)%isf_south(i)
              IF (hi(i,Jstr-1).le.min_his(ng)) THEN
                ei(i,Jstr-1)=0.0_r8
                ti(i,Jstr-1)=0.0_r8
              ELSE
                ti(i,Jstr-1)=ISF_LOBC(isIeis,ng)%isf_south(i)
              END IF
# ifdef MASKING
              ei(i,Jstr-1)=ei(i,Jstr-1)*GRID(ng)%rmask(i,Jstr-1)
              ti(i,Jstr-1)=ti(i,Jstr-1)*GRID(ng)%rmask(i,Jstr-1)
# endif
            END DO
          ELSE IF (LBC(isouth,ibISF(isIeis),ng)%gradient .or.           &
     &             LBC(isouth,ibISF(isIeis),ng)%closed) THEN
            DO i=Istr,Iend
              ei(i,Jstr-1)=hi(i,Jstr)*ti(i,Jstr)
              IF (hi(i,Jstr-1).le.min_his(ng)) THEN
                ei(i,Jstr-1)=0.0_r8
                ti(i,Jstr-1)=0.0_r8
              ELSE
                ti(i,Jstr-1)=ti(i,Jstr)
              END IF
# ifdef MASKING
              ei(i,Jstr-1)=ei(i,Jstr-1)*GRID(ng)%rmask(i,Jstr-1)
              ti(i,Jstr-1)=ti(i,Jstr-1)*GRID(ng)%rmask(i,Jstr-1)
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
          IF (LBC(inorth,ibISF(isIeis),ng)%clamped) THEN
            DO i=Istr,Iend
              ei(i,Jend+1)=ISF_LOBC(isIhis,ng)%isf_north(i)*            &
     &                     ISF_LOBC(isIeis,ng)%isf_north(i)
              IF (hi(i,Jend+1).le.min_his(ng)) THEN
                ei(i,Jend+1)=0.0_r8
                ti(i,Jend+1)=0.0_r8
              ELSE
                ti(i,Jend+1)=ISF_LOBC(isIeis,ng)%isf_north(i)
              END IF
# ifdef MASKING
              ei(i,Jend+1)=ei(i,Jend+1)*GRID(ng)%rmask(i,Jend+1)
              ti(i,Jend+1)=ti(i,Jend+1)*GRID(ng)%rmask(i,Jend+1)
# endif
            END DO
          ELSE IF (LBC(inorth,ibISF(isIeis),ng)%gradient .or.           &
     &             LBC(inorth,ibISF(isIeis),ng)%closed) THEN
            DO i=Istr,Iend
              ei(i,Jend+1)=hi(i,Jend)*ti(i,Jend)
              IF (hi(i,Jend+1).le.min_his(ng)) THEN
                ei(i,Jend+1)=0.0_r8
                ti(i,Jend+1)=0.0_r8
              ELSE
                ti(i,Jend+1)=ti(i,Jend)
              END IF
# ifdef MASKING
              ei(i,Jend+1)=ei(i,Jend+1)*GRID(ng)%rmask(i,Jend+1)
              ti(i,Jend+1)=ti(i,Jend+1)*GRID(ng)%rmask(i,Jend+1)
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
          ei(Istr-1,Jstr-1)=0.5_r8*(hi(Istr,Jstr-1)*ti(Istr,Jstr-1) +   &
     &                              hi(Istr-1,Jstr)*ti(Istr-1,Jstr))
          IF (hi(Istr-1,Jstr-1).le.min_his(ng)) THEN
            ei(Istr-1,Jstr-1)=0.0_r8
            ti(Istr-1,Jstr-1)=0.0_r8
          ELSE
            ti(Istr-1,Jstr-1)=ei(Istr-1,Jstr-1) /                       &
     &                        MAX(hi(Istr-1,Jstr-1),min_his(ng))
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          ei(Iend+1,Jstr-1)=0.5_r8*(hi(Iend+1,Jstr)*ti(Iend+1,Jstr) +   &
     &                              hi(Iend,Jstr-1)*ti(Iend,Jstr-1))
          IF (hi(Iend+1,Jstr-1).le.min_his(ng)) THEN
            ei(Iend+1,Jstr-1)=0.0_r8
            ti(Iend+1,Jstr-1)=0.0_r8
          ELSE
            ti(Iend+1,Jstr-1)=ei(Iend+1,Jstr-1) /                       &
     &                        MAX(hi(Iend+1,Jstr-1),min_his(ng))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          ei(Istr-1,Jend+1)=0.5_r8*(hi(Istr-1,Jend)*ti(Istr-1,Jend) +   &
     &                              hi(Istr,Jend+1)*ti(Istr,Jend+1))
          IF (hi(Istr-1,Jend+1).le.min_his(ng)) THEN
            ei(Istr-1,Jend+1)=0.0_r8
            ti(Istr-1,Jend+1)=0.0_r8
          ELSE
            ti(Istr-1,Jend+1)=ei(Istr-1,Jend+1) /                       &
     &                        MAX(hi(Istr-1,Jend+1),min_his(ng))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          ei(Iend+1,Jend+1)=0.5_r8*(hi(Iend+1,Jend)*ti(Iend+1,Jend) +   &
     &                              hi(Iend,Jend+1)*ti(Iend,Jend+1))
          IF (hi(Iend+1,Jend+1).le.min_his(ng)) THEN
            ei(Iend+1,Jend+1)=0.0_r8
            ti(Iend+1,Jend+1)=0.0_r8
          ELSE
            ti(Iend+1,Jend+1)=ei(Iend+1,Jend+1) /                       &
     &                        MAX(hi(Iend+1,Jend+1),min_his(ng))
          END IF
        END IF
      END IF
!
      RETURN
      END SUBROUTINE isf_tibc_tile
#endif
      END MODULE isf_tibc_mod
