#include "cppdefs.h"
      MODULE isf_mms_source_mod
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2024 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module reads MMS (Method of Manufactured Solutions) source     !
!  terms from a NetCDF file for SSA solver verification testing.       !
!  The source terms are added to the CG solver RHS to force the        !
!  solution toward a known manufactured solution.                      !
!                                                                      !
!  Only active when SSA_MMS_TESTING is defined.                        !
!                                                                      !
!  The NetCDF file is read once per process (first SSA call for a      !
!  given grid) and cached in module-level arrays; subsequent tile      !
!  calls copy their local range from memory.  This avoids redundant    !
!  file I/O from every tile on every SSA solve.                        !
!                                                                      !
!=======================================================================
!
#ifdef SSA_MMS_TESTING
      USE mod_kinds
      USE mod_param,      ONLY : Ngrids, Lm, Mm
      USE mod_parallel,   ONLY : Master
      USE mod_scalars,    ONLY : exit_flag, NoError
!
      USE netcdf
!
      implicit none
!
      PRIVATE
      PUBLIC :: isf_read_mms_source
!
!  Per-process cache of the MMS source fields (PSI-point, global
!  dimensions (Lm+1) x (Mm+1)), one entry per nested grid so that
!  multi-grid configurations never re-read the file after the first
!  call for each grid.
!
!  THREAD SAFETY: the cache and the NetCDF read are protected by an
!  OMP critical section below.  Under OpenMP, multiple tiles in the
!  same process call this routine concurrently; without protection the
!  NF90_GET_VAR calls and the allocate/deallocate of the module-level
!  arrays would race, corrupting the cache.
!
      TYPE T_MMS_CACHE
        logical :: cached=.FALSE.
        real(r8), allocatable :: fu(:,:)
        real(r8), allocatable :: fv(:,:)
      END TYPE T_MMS_CACHE
!
      type (T_MMS_CACHE), allocatable, save :: mms_cache(:)
!
      CONTAINS
!
      SUBROUTINE isf_read_mms_source (ng, tile, LBi, UBi, LBj, UBj,    &
     &                                mms_fu_psi, mms_fv_psi)
!
!=======================================================================
!                                                                      !
!  This routine returns the MMS source terms (mms_fu_psi, mms_fv_psi) !
!  on PSI-points for the local tile.  On first call for a given grid   !
!  the fields are read from the NetCDF file specified by               !
!  mms_source_file in isf.in (dimensions (xi_psi, eta_psi) =           !
!  (Lm+1, Mm+1)); subsequent calls copy from the module-level cache.   !
!                                                                      !
!  On Input:                                                           !
!     ng            Nested grid number                                 !
!     tile          Domain partition tile                              !
!     LBi, UBi      I-dimension bounds                                 !
!     LBj, UBj      J-dimension bounds                                 !
!                                                                      !
!  On Output:                                                          !
!     mms_fu_psi    MMS u-source term on PSI-points [Pa m]            !
!     mms_fv_psi    MMS v-source term on PSI-points [Pa m]            !
!                                                                      !
!=======================================================================
!
      USE mod_isf,        ONLY : mms_source_file
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, LBi, UBi, LBj, UBj
!
      real(r8), intent(out) :: mms_fu_psi(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: mms_fv_psi(LBi:UBi,LBj:UBj)
!
!  Local variable declarations.
!
      integer :: ncid, varid_u, varid_v
      integer :: i, j, status
      integer :: start(2), count(2)
      integer :: nx_psi, ny_psi
!
      character (len=256) :: fname
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", isf_read_mms_source"
!
!-----------------------------------------------------------------------
!  Initialize output arrays to zero.
!-----------------------------------------------------------------------
!
      DO j=LBj,UBj
        DO i=LBi,UBi
          mms_fu_psi(i,j)=0.0_r8
          mms_fv_psi(i,j)=0.0_r8
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  If the source file for this grid is already cached, copy the local
!  tile range from memory and return.
!
!  The entire cache-check + NetCDF-read path is wrapped in an OMP
!  critical section to prevent concurrent access from multiple tiles
!  under OpenMP.  The critical section is only entered on the first
!  call per grid; subsequent calls (cache hit) execute quickly.
!-----------------------------------------------------------------------
!
!$OMP CRITICAL (mms_cache)
      IF (.not.allocated(mms_cache)) THEN
        allocate ( mms_cache(Ngrids) )
      END IF
      IF (mms_cache(ng) % cached) THEN
        nx_psi=SIZE(mms_cache(ng) % fu,1)
        ny_psi=SIZE(mms_cache(ng) % fu,2)
        DO j=LBj,UBj
          DO i=LBi,UBi
            IF (i.ge.1.and.i.le.nx_psi.and.                           &
     &          j.ge.1.and.j.le.ny_psi) THEN
              mms_fu_psi(i,j)=mms_cache(ng) % fu(i,j)
              mms_fv_psi(i,j)=mms_cache(ng) % fv(i,j)
            END IF
          END DO
        END DO
        GO TO 100
      END IF
!
!-----------------------------------------------------------------------
!  First call for this grid: read the full field from the NetCDF file.
!  PSI-point dimensions: xi_psi = Lm(ng)+1, eta_psi = Mm(ng)+1.
!-----------------------------------------------------------------------
!
      fname=TRIM(mms_source_file(ng))
!
      status=NF90_OPEN(TRIM(fname), NF90_NOWRITE, ncid)
      IF (status.ne.NF90_NOERR) THEN
        IF (Master) WRITE (*,10) TRIM(fname),                           &
     &    TRIM(NF90_STRERROR(status))
        exit_flag=2
        GO TO 100
      END IF
!
!  Inquire variable IDs.
!
      status=NF90_INQ_VARID(ncid, 'mms_fu_psi', varid_u)
      IF (status.ne.NF90_NOERR) THEN
        IF (Master) WRITE (*,20) 'mms_fu_psi', TRIM(fname),            &
     &    TRIM(NF90_STRERROR(status))
        exit_flag=2
        GO TO 100
      END IF
!
      status=NF90_INQ_VARID(ncid, 'mms_fv_psi', varid_v)
      IF (status.ne.NF90_NOERR) THEN
        IF (Master) WRITE (*,20) 'mms_fv_psi', TRIM(fname),            &
     &    TRIM(NF90_STRERROR(status))
        exit_flag=2
        GO TO 100
      END IF
!
      nx_psi=Lm(ng)+1
      ny_psi=Mm(ng)+1
!
      IF (.not.allocated(mms_cache(ng) % fu)) THEN
        allocate ( mms_cache(ng) % fu(nx_psi,ny_psi) )
        allocate ( mms_cache(ng) % fv(nx_psi,ny_psi) )
      END IF
!
      start=(/ 1, 1 /)
      count=(/ nx_psi, ny_psi /)
!
      status=NF90_GET_VAR(ncid, varid_u, mms_cache(ng) % fu,           &
     &                    start=start, count=count)
      IF (status.ne.NF90_NOERR) THEN
        IF (Master) WRITE (*,30) 'mms_fu_psi',                          &
     &    TRIM(NF90_STRERROR(status))
        exit_flag=2
        GO TO 100
      END IF
!
      status=NF90_GET_VAR(ncid, varid_v, mms_cache(ng) % fv,           &
     &                    start=start, count=count)
      IF (status.ne.NF90_NOERR) THEN
        IF (Master) WRITE (*,30) 'mms_fv_psi',                          &
     &    TRIM(NF90_STRERROR(status))
        exit_flag=2
        GO TO 100
      END IF
!
      mms_cache(ng) % cached=.TRUE.
!
!-----------------------------------------------------------------------
!  Clean up.
!-----------------------------------------------------------------------
!
      status=NF90_CLOSE(ncid)
!
      IF (Master) WRITE (*,40) TRIM(fname)
!
!  Copy the full field to the local tile range (1-based NetCDF global
!  indexing to ROMS array bounds).
!
      DO j=LBj,UBj
        DO i=LBi,UBi
          IF (i.ge.1.and.i.le.nx_psi.and.                               &
     &        j.ge.1.and.j.le.ny_psi) THEN
            mms_fu_psi(i,j)=mms_cache(ng) % fu(i,j)
            mms_fv_psi(i,j)=mms_cache(ng) % fv(i,j)
          END IF
        END DO
      END DO
!
  100 CONTINUE
!$OMP END CRITICAL (mms_cache)
!
  10  FORMAT (/,' ISF_READ_MMS_SOURCE - unable to open MMS source',    &
     &        ' file: ',a,/,19x,'NF90 error: ',a)
  20  FORMAT (/,' ISF_READ_MMS_SOURCE - unable to find variable: ',a,  &
     &        /,19x,'in file: ',a,/,19x,'NF90 error: ',a)
  30  FORMAT (/,' ISF_READ_MMS_SOURCE - unable to read variable: ',a,  &
     &        /,19x,'NF90 error: ',a)
  40  FORMAT (/,' ISF_READ_MMS_SOURCE - read MMS source terms from: ', &
     &        a)
!
      RETURN
      END SUBROUTINE isf_read_mms_source
!
#endif
      END MODULE isf_mms_source_mod
