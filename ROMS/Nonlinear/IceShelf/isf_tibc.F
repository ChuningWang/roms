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
!  in the diagnostic field Fi(:,:,icItis) so that downstream physics   !
!  (e.g. the Glen rate factor) sees a consistent ice temperature.      !
!                                                                      !
!  Supported boundary condition types (see isf.in):                    !
!     clamped    - E = h_bnd * T_bnd                                   !
!     mixed      - clamp on inflow, gradient on outflow                !
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
      CALL isf_tibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    lsold(ng), lsnew(ng),                         &
     &                    ISF(ng) % Si(:,:,:,isIuis),                   &
     &                    ISF(ng) % Si(:,:,:,isIvis),                   &
     &                    ISF(ng) % Si(:,:,:,isIhis),                   &
     &                    ISF(ng) % Si(:,:,:,isIeis),                   &
     &                    ISF(ng) % Fi(:,:,:))
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
     &                          lsold, lsnew,                           &
     &                          ui, vi, hi, enthalpy, Fi)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: lsold, lsnew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)    :: ui(LBi:,LBj:,:)
      real(r8), intent(in)    :: vi(LBi:,LBj:,:)
      real(r8), intent(in)    :: hi(LBi:,LBj:,:)
      real(r8), intent(inout) :: enthalpy(LBi:,LBj:,:)
      real(r8), intent(inout) :: Fi(LBi:,LBj:,:)
# else
      real(r8), intent(in)    :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: hi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: enthalpy(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: Fi(LBi:UBi,LBj:UBj,nIsfF)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
!
      real(r8) :: T_bnd
!
      real(r8), parameter :: eps=1.0E-6_r8

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        IF (LBC(iwest,ibISF(isIeis),ng)%clamped) THEN
          DO j=Jstr,Jend
            T_bnd=ISF_LOBC(isIeis,ng)%isf_west(j)
            enthalpy(Istr-1,j,lsnew)=ISF_LOBC(isIhis,ng)%isf_west(j)*   &
     &                               T_bnd
# ifdef MASKING
            enthalpy(Istr-1,j,lsnew)=enthalpy(Istr-1,j,lsnew)*          &
     &                               GRID(ng)%rmask(Istr-1,j)
# endif
# ifdef WET_DRY
            enthalpy(Istr-1,j,lsnew)=enthalpy(Istr-1,j,lsnew)*          &
     &                               GRID(ng)%rmask_wet(Istr-1,j)
# endif
            Fi(Istr-1,j,icItis)=T_bnd
            IF (hi(Istr-1,j,lsnew).le.min_his) THEN
              enthalpy(Istr-1,j,lsnew)=0.0_r8
              Fi(Istr-1,j,icItis)=0.0_r8
            END IF
          END DO
!
        ELSE IF (LBC(iwest,ibISF(isIeis),ng)%mixed) THEN
          DO j=Jstr,Jend
            IF (ui(Istr,j,lsnew).ge.0.0_r8) THEN
              T_bnd=ISF_LOBC(isIeis,ng)%isf_west(j)
              enthalpy(Istr-1,j,lsnew)=ISF_LOBC(isIhis,ng)%isf_west(j)* &
     &                                 T_bnd
              Fi(Istr-1,j,icItis)=T_bnd
            ELSE
              enthalpy(Istr-1,j,lsnew)=enthalpy(Istr,j,lsold)
              Fi(Istr-1,j,icItis)=enthalpy(Istr-1,j,lsnew)/             &
     &                            MAX(hi(Istr-1,j,lsnew),eps)
            END IF
# ifdef MASKING
            enthalpy(Istr-1,j,lsnew)=enthalpy(Istr-1,j,lsnew)*          &
     &                               GRID(ng)%rmask(Istr-1,j)
# endif
# ifdef WET_DRY
            enthalpy(Istr-1,j,lsnew)=enthalpy(Istr-1,j,lsnew)*          &
     &                               GRID(ng)%rmask_wet(Istr-1,j)
# endif
            IF (hi(Istr-1,j,lsnew).le.min_his) THEN
              enthalpy(Istr-1,j,lsnew)=0.0_r8
              Fi(Istr-1,j,icItis)=0.0_r8
            END IF
          END DO
!
        ELSE IF (LBC(iwest,ibISF(isIeis),ng)%gradient) THEN
          DO j=Jstr,Jend
            enthalpy(Istr-1,j,lsnew)=enthalpy(Istr,j,lsnew)
# ifdef MASKING
            enthalpy(Istr-1,j,lsnew)=enthalpy(Istr-1,j,lsnew)*          &
     &                               GRID(ng)%rmask(Istr-1,j)
# endif
# ifdef WET_DRY
            enthalpy(Istr-1,j,lsnew)=enthalpy(Istr-1,j,lsnew)*          &
     &                               GRID(ng)%rmask_wet(Istr-1,j)
# endif
            Fi(Istr-1,j,icItis)=enthalpy(Istr-1,j,lsnew)/               &
     &                          MAX(hi(Istr-1,j,lsnew),eps)
            IF (hi(Istr-1,j,lsnew).le.min_his) THEN
              enthalpy(Istr-1,j,lsnew)=0.0_r8
              Fi(Istr-1,j,icItis)=0.0_r8
            END IF
          END DO
!
        ELSE IF (LBC(iwest,ibISF(isIeis),ng)%closed) THEN
          DO j=Jstr,Jend
            enthalpy(Istr-1,j,lsnew)=enthalpy(Istr,j,lsnew)
# ifdef MASKING
            enthalpy(Istr-1,j,lsnew)=enthalpy(Istr-1,j,lsnew)*          &
     &                               GRID(ng)%rmask(Istr-1,j)
# endif
# ifdef WET_DRY
            enthalpy(Istr-1,j,lsnew)=enthalpy(Istr-1,j,lsnew)*          &
     &                               GRID(ng)%rmask_wet(Istr-1,j)
# endif
            Fi(Istr-1,j,icItis)=enthalpy(Istr-1,j,lsnew)/               &
     &                          MAX(hi(Istr-1,j,lsnew),eps)
            IF (hi(Istr-1,j,lsnew).le.min_his) THEN
              enthalpy(Istr-1,j,lsnew)=0.0_r8
              Fi(Istr-1,j,icItis)=0.0_r8
            END IF
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
            T_bnd=ISF_LOBC(isIeis,ng)%isf_east(j)
            enthalpy(Iend+1,j,lsnew)=ISF_LOBC(isIhis,ng)%isf_east(j)*   &
     &                               T_bnd
# ifdef MASKING
            enthalpy(Iend+1,j,lsnew)=enthalpy(Iend+1,j,lsnew)*          &
     &                               GRID(ng)%rmask(Iend+1,j)
# endif
# ifdef WET_DRY
            enthalpy(Iend+1,j,lsnew)=enthalpy(Iend+1,j,lsnew)*          &
     &                               GRID(ng)%rmask_wet(Iend+1,j)
# endif
            Fi(Iend+1,j,icItis)=T_bnd
            IF (hi(Iend+1,j,lsnew).le.min_his) THEN
              enthalpy(Iend+1,j,lsnew)=0.0_r8
              Fi(Iend+1,j,icItis)=0.0_r8
            END IF
          END DO
!
        ELSE IF (LBC(ieast,ibISF(isIeis),ng)%mixed) THEN
          DO j=Jstr,Jend
            IF (ui(Iend+1,j,lsnew).le.0.0_r8) THEN
              T_bnd=ISF_LOBC(isIeis,ng)%isf_east(j)
              enthalpy(Iend+1,j,lsnew)=ISF_LOBC(isIhis,ng)%isf_east(j)* &
     &                                 T_bnd
              Fi(Iend+1,j,icItis)=T_bnd
            ELSE
              enthalpy(Iend+1,j,lsnew)=enthalpy(Iend,j,lsold)
              Fi(Iend+1,j,icItis)=enthalpy(Iend+1,j,lsnew)/             &
     &                            MAX(hi(Iend+1,j,lsnew),eps)
            END IF
# ifdef MASKING
            enthalpy(Iend+1,j,lsnew)=enthalpy(Iend+1,j,lsnew)*          &
     &                               GRID(ng)%rmask(Iend+1,j)
# endif
# ifdef WET_DRY
            enthalpy(Iend+1,j,lsnew)=enthalpy(Iend+1,j,lsnew)*          &
     &                               GRID(ng)%rmask_wet(Iend+1,j)
# endif
            IF (hi(Iend+1,j,lsnew).le.min_his) THEN
              enthalpy(Iend+1,j,lsnew)=0.0_r8
              Fi(Iend+1,j,icItis)=0.0_r8
            END IF
          END DO
!
        ELSE IF (LBC(ieast,ibISF(isIeis),ng)%gradient) THEN
          DO j=Jstr,Jend
            enthalpy(Iend+1,j,lsnew)=enthalpy(Iend,j,lsnew)
# ifdef MASKING
            enthalpy(Iend+1,j,lsnew)=enthalpy(Iend+1,j,lsnew)*          &
     &                               GRID(ng)%rmask(Iend+1,j)
# endif
# ifdef WET_DRY
            enthalpy(Iend+1,j,lsnew)=enthalpy(Iend+1,j,lsnew)*          &
     &                               GRID(ng)%rmask_wet(Iend+1,j)
# endif
            Fi(Iend+1,j,icItis)=enthalpy(Iend+1,j,lsnew)/               &
     &                          MAX(hi(Iend+1,j,lsnew),eps)
            IF (hi(Iend+1,j,lsnew).le.min_his) THEN
              enthalpy(Iend+1,j,lsnew)=0.0_r8
              Fi(Iend+1,j,icItis)=0.0_r8
            END IF
          END DO
!
        ELSE IF (LBC(ieast,ibISF(isIeis),ng)%closed) THEN
          DO j=Jstr,Jend
            enthalpy(Iend+1,j,lsnew)=enthalpy(Iend,j,lsnew)
# ifdef MASKING
            enthalpy(Iend+1,j,lsnew)=enthalpy(Iend+1,j,lsnew)*          &
     &                               GRID(ng)%rmask(Iend+1,j)
# endif
# ifdef WET_DRY
            enthalpy(Iend+1,j,lsnew)=enthalpy(Iend+1,j,lsnew)*          &
     &                               GRID(ng)%rmask_wet(Iend+1,j)
# endif
            Fi(Iend+1,j,icItis)=enthalpy(Iend+1,j,lsnew)/               &
     &                          MAX(hi(Iend+1,j,lsnew),eps)
            IF (hi(Iend+1,j,lsnew).le.min_his) THEN
              enthalpy(Iend+1,j,lsnew)=0.0_r8
              Fi(Iend+1,j,icItis)=0.0_r8
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        IF (LBC(isouth,ibISF(isIeis),ng)%clamped) THEN
          DO i=Istr,Iend
            T_bnd=ISF_LOBC(isIeis,ng)%isf_south(i)
            enthalpy(i,Jstr-1,lsnew)=ISF_LOBC(isIhis,ng)%isf_south(i)*  &
     &                               T_bnd
# ifdef MASKING
            enthalpy(i,Jstr-1,lsnew)=enthalpy(i,Jstr-1,lsnew)*          &
     &                               GRID(ng)%rmask(i,Jstr-1)
# endif
# ifdef WET_DRY
            enthalpy(i,Jstr-1,lsnew)=enthalpy(i,Jstr-1,lsnew)*          &
     &                               GRID(ng)%rmask_wet(i,Jstr-1)
# endif
            Fi(i,Jstr-1,icItis)=T_bnd
            IF (hi(i,Jstr-1,lsnew).le.min_his) THEN
              enthalpy(i,Jstr-1,lsnew)=0.0_r8
              Fi(i,Jstr-1,icItis)=0.0_r8
            END IF
          END DO
!
        ELSE IF (LBC(isouth,ibISF(isIeis),ng)%mixed) THEN
          DO i=Istr,Iend
            IF (vi(i,Jstr,lsnew).ge.0.0_r8) THEN
              T_bnd=ISF_LOBC(isIeis,ng)%isf_south(i)
              enthalpy(i,Jstr-1,lsnew)=ISF_LOBC(isIhis,ng)%isf_south(i)*&
     &                                 T_bnd
              Fi(i,Jstr-1,icItis)=T_bnd
            ELSE
              enthalpy(i,Jstr-1,lsnew)=enthalpy(i,Jstr,lsold)
              Fi(i,Jstr-1,icItis)=enthalpy(i,Jstr-1,lsnew)/             &
     &                            MAX(hi(i,Jstr-1,lsnew),eps)
            END IF
# ifdef MASKING
            enthalpy(i,Jstr-1,lsnew)=enthalpy(i,Jstr-1,lsnew)*          &
     &                               GRID(ng)%rmask(i,Jstr-1)
# endif
# ifdef WET_DRY
            enthalpy(i,Jstr-1,lsnew)=enthalpy(i,Jstr-1,lsnew)*          &
     &                               GRID(ng)%rmask_wet(i,Jstr-1)
# endif
            IF (hi(i,Jstr-1,lsnew).le.min_his) THEN
              enthalpy(i,Jstr-1,lsnew)=0.0_r8
              Fi(i,Jstr-1,icItis)=0.0_r8
            END IF
          END DO
!
        ELSE IF (LBC(isouth,ibISF(isIeis),ng)%gradient) THEN
          DO i=Istr,Iend
            enthalpy(i,Jstr-1,lsnew)=enthalpy(i,Jstr,lsnew)
# ifdef MASKING
            enthalpy(i,Jstr-1,lsnew)=enthalpy(i,Jstr-1,lsnew)*          &
     &                               GRID(ng)%rmask(i,Jstr-1)
# endif
# ifdef WET_DRY
            enthalpy(i,Jstr-1,lsnew)=enthalpy(i,Jstr-1,lsnew)*          &
     &                               GRID(ng)%rmask_wet(i,Jstr-1)
# endif
            Fi(i,Jstr-1,icItis)=enthalpy(i,Jstr-1,lsnew)/               &
     &                          MAX(hi(i,Jstr-1,lsnew),eps)
            IF (hi(i,Jstr-1,lsnew).le.min_his) THEN
              enthalpy(i,Jstr-1,lsnew)=0.0_r8
              Fi(i,Jstr-1,icItis)=0.0_r8
            END IF
          END DO
!
        ELSE IF (LBC(isouth,ibISF(isIeis),ng)%closed) THEN
          DO i=Istr,Iend
            enthalpy(i,Jstr-1,lsnew)=enthalpy(i,Jstr,lsnew)
# ifdef MASKING
            enthalpy(i,Jstr-1,lsnew)=enthalpy(i,Jstr-1,lsnew)*          &
     &                               GRID(ng)%rmask(i,Jstr-1)
# endif
# ifdef WET_DRY
            enthalpy(i,Jstr-1,lsnew)=enthalpy(i,Jstr-1,lsnew)*          &
     &                               GRID(ng)%rmask_wet(i,Jstr-1)
# endif
            Fi(i,Jstr-1,icItis)=enthalpy(i,Jstr-1,lsnew)/               &
     &                          MAX(hi(i,Jstr-1,lsnew),eps)
            IF (hi(i,Jstr-1,lsnew).le.min_his) THEN
              enthalpy(i,Jstr-1,lsnew)=0.0_r8
              Fi(i,Jstr-1,icItis)=0.0_r8
            END IF
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
            T_bnd=ISF_LOBC(isIeis,ng)%isf_north(i)
            enthalpy(i,Jend+1,lsnew)=ISF_LOBC(isIhis,ng)%isf_north(i)*  &
     &                               T_bnd
# ifdef MASKING
            enthalpy(i,Jend+1,lsnew)=enthalpy(i,Jend+1,lsnew)*          &
     &                               GRID(ng)%rmask(i,Jend+1)
# endif
# ifdef WET_DRY
            enthalpy(i,Jend+1,lsnew)=enthalpy(i,Jend+1,lsnew)*          &
     &                               GRID(ng)%rmask_wet(i,Jend+1)
# endif
            Fi(i,Jend+1,icItis)=T_bnd
            IF (hi(i,Jend+1,lsnew).le.min_his) THEN
              enthalpy(i,Jend+1,lsnew)=0.0_r8
              Fi(i,Jend+1,icItis)=0.0_r8
            END IF
          END DO
!
        ELSE IF (LBC(inorth,ibISF(isIeis),ng)%mixed) THEN
          DO i=Istr,Iend
            IF (vi(i,Jend+1,lsnew).le.0.0_r8) THEN
              T_bnd=ISF_LOBC(isIeis,ng)%isf_north(i)
              enthalpy(i,Jend+1,lsnew)=ISF_LOBC(isIhis,ng)%isf_north(i)*&
     &                                 T_bnd
              Fi(i,Jend+1,icItis)=T_bnd
            ELSE
              enthalpy(i,Jend+1,lsnew)=enthalpy(i,Jend,lsold)
              Fi(i,Jend+1,icItis)=enthalpy(i,Jend+1,lsnew)/             &
     &                            MAX(hi(i,Jend+1,lsnew),eps)
            END IF
# ifdef MASKING
            enthalpy(i,Jend+1,lsnew)=enthalpy(i,Jend+1,lsnew)*          &
     &                               GRID(ng)%rmask(i,Jend+1)
# endif
# ifdef WET_DRY
            enthalpy(i,Jend+1,lsnew)=enthalpy(i,Jend+1,lsnew)*          &
     &                               GRID(ng)%rmask_wet(i,Jend+1)
# endif
            IF (hi(i,Jend+1,lsnew).le.min_his) THEN
              enthalpy(i,Jend+1,lsnew)=0.0_r8
              Fi(i,Jend+1,icItis)=0.0_r8
            END IF
          END DO
!
        ELSE IF (LBC(inorth,ibISF(isIeis),ng)%gradient) THEN
          DO i=Istr,Iend
            enthalpy(i,Jend+1,lsnew)=enthalpy(i,Jend,lsnew)
# ifdef MASKING
            enthalpy(i,Jend+1,lsnew)=enthalpy(i,Jend+1,lsnew)*          &
     &                               GRID(ng)%rmask(i,Jend+1)
# endif
# ifdef WET_DRY
            enthalpy(i,Jend+1,lsnew)=enthalpy(i,Jend+1,lsnew)*          &
     &                               GRID(ng)%rmask_wet(i,Jend+1)
# endif
            Fi(i,Jend+1,icItis)=enthalpy(i,Jend+1,lsnew)/               &
     &                          MAX(hi(i,Jend+1,lsnew),eps)
            IF (hi(i,Jend+1,lsnew).le.min_his) THEN
              enthalpy(i,Jend+1,lsnew)=0.0_r8
              Fi(i,Jend+1,icItis)=0.0_r8
            END IF
          END DO
!
        ELSE IF (LBC(inorth,ibISF(isIeis),ng)%closed) THEN
          DO i=Istr,Iend
            enthalpy(i,Jend+1,lsnew)=enthalpy(i,Jend,lsnew)
# ifdef MASKING
            enthalpy(i,Jend+1,lsnew)=enthalpy(i,Jend+1,lsnew)*          &
     &                               GRID(ng)%rmask(i,Jend+1)
# endif
# ifdef WET_DRY
            enthalpy(i,Jend+1,lsnew)=enthalpy(i,Jend+1,lsnew)*          &
     &                               GRID(ng)%rmask_wet(i,Jend+1)
# endif
            Fi(i,Jend+1,icItis)=enthalpy(i,Jend+1,lsnew)/               &
     &                          MAX(hi(i,Jend+1,lsnew),eps)
            IF (hi(i,Jend+1,lsnew).le.min_his) THEN
              enthalpy(i,Jend+1,lsnew)=0.0_r8
              Fi(i,Jend+1,icItis)=0.0_r8
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          enthalpy(Istr-1,Jstr-1,lsnew)=0.5_r8*                         &
     &                                  (enthalpy(Istr  ,Jstr-1,lsnew)+ &
     &                                   enthalpy(Istr-1,Jstr  ,lsnew))
# ifdef MASKING
          enthalpy(Istr-1,Jstr-1,lsnew)=enthalpy(Istr-1,Jstr-1,lsnew)*  &
     &                                  GRID(ng)%rmask(Istr-1,Jstr-1)
# endif
# ifdef WET_DRY
          enthalpy(Istr-1,Jstr-1,lsnew)=enthalpy(Istr-1,Jstr-1,lsnew)*  &
     &                                 GRID(ng)%rmask_wet(Istr-1,Jstr-1)
# endif
          Fi(Istr-1,Jstr-1,icItis)=enthalpy(Istr-1,Jstr-1,lsnew)/       &
     &                             MAX(hi(Istr-1,Jstr-1,lsnew),eps)
          IF (hi(Istr-1,Jstr-1,lsnew).le.min_his) THEN
            enthalpy(Istr-1,Jstr-1,lsnew)=0.0_r8
            Fi(Istr-1,Jstr-1,icItis)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          enthalpy(Iend+1,Jstr-1,lsnew)=0.5_r8*                         &
     &                                  (enthalpy(Iend+1,Jstr  ,lsnew)+ &
     &                                   enthalpy(Iend  ,Jstr-1,lsnew))
# ifdef MASKING
          enthalpy(Iend+1,Jstr-1,lsnew)=enthalpy(Iend+1,Jstr-1,lsnew)*  &
     &                                  GRID(ng)%rmask(Iend+1,Jstr-1)
# endif
# ifdef WET_DRY
          enthalpy(Iend+1,Jstr-1,lsnew)=enthalpy(Iend+1,Jstr-1,lsnew)*  &
     &                                 GRID(ng)%rmask_wet(Iend+1,Jstr-1)
# endif
          Fi(Iend+1,Jstr-1,icItis)=enthalpy(Iend+1,Jstr-1,lsnew)/       &
     &                             MAX(hi(Iend+1,Jstr-1,lsnew),eps)
          IF (hi(Iend+1,Jstr-1,lsnew).LE.min_his) THEN
            enthalpy(Iend+1,Jstr-1,lsnew)=0.0_r8
            Fi(Iend+1,Jstr-1,icItis)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          enthalpy(Istr-1,Jend+1,lsnew)=0.5_r8*                         &
     &                                  (enthalpy(Istr-1,Jend  ,lsnew)+ &
     &                                   enthalpy(Istr  ,Jend+1,lsnew))
# ifdef MASKING
          enthalpy(Istr-1,Jend+1,lsnew)=enthalpy(Istr-1,Jend+1,lsnew)*  &
     &                                  GRID(ng)%rmask(Istr-1,Jend+1)
# endif
# ifdef WET_DRY
          enthalpy(Istr-1,Jend+1,lsnew)=enthalpy(Istr-1,Jend+1,lsnew)*  &
     &                                 GRID(ng)%rmask_wet(Istr-1,Jend+1)
# endif
          Fi(Istr-1,Jend+1,icItis)=enthalpy(Istr-1,Jend+1,lsnew)/       &
     &                             MAX(hi(Istr-1,Jend+1,lsnew),eps)
          IF (hi(Istr-1,Jend+1,lsnew).le.min_his) THEN
            enthalpy(Istr-1,Jend+1,lsnew)=0.0_r8
            Fi(Istr-1,Jend+1,icItis)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          enthalpy(Iend+1,Jend+1,lsnew)=0.5_r8*                         &
     &                                  (enthalpy(Iend+1,Jend  ,lsnew)+ &
     &                                   enthalpy(Iend  ,Jend+1,lsnew))
# ifdef MASKING
          enthalpy(Iend+1,Jend+1,lsnew)=enthalpy(Iend+1,Jend+1,lsnew)*  &
     &                                  GRID(ng)%rmask(Iend+1,Jend+1)
# endif
# ifdef WET_DRY
          enthalpy(Iend+1,Jend+1,lsnew)=enthalpy(Iend+1,Jend+1,lsnew)*  &
     &                                 GRID(ng)%rmask_wet(Iend+1,Jend+1)
# endif
          Fi(Iend+1,Jend+1,icItis)=enthalpy(Iend+1,Jend+1,lsnew)/       &
     &                             MAX(hi(Iend+1,Jend+1,lsnew),eps)
          IF (hi(Iend+1,Jend+1,lsnew).le.min_his) THEN
            enthalpy(Iend+1,Jend+1,lsnew)=0.0_r8
            Fi(Iend+1,Jend+1,icItis)=0.0_r8
          END IF
        END IF
      END IF
!
      RETURN
      END SUBROUTINE isf_tibc_tile
#endif
      END MODULE isf_tibc_mod
