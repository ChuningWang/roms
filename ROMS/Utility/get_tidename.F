#include "cppdefs.h"
      SUBROUTINE get_tidename (ng, model, ifield, ncid,                 &
#if defined PIO_LIB && defined DISTRIBUTE
     &                         pioFile,                                 &
#endif
     &                         nfiles, S,                               &
     &                         ntide,                                   &
     &                         Cout)
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2025 The ROMS Group               Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This routine reads in name of tidal harmonic constituents from      !
!  NetCDF file.                                                        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     ifield     Field ID.                                             !
!     ncid       NetCDF file ID.                                       !
# if defined PIO_LIB && defined DISTRIBUTE
!     pioFile    PIO file descriptor structure, TYPE(file_desc_t)      !
!                  pioFile%fh         file handler                     !
!                  pioFile%iosystem   IO system descriptor (struct)    !
# endif
!     nfiles     Number of input NetCDF files.                         !
!     S          I/O derived type structure, TYPE(T_IO).               !
!     ntide      Number of tides to read.                              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Cout       Tidal constituent names.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_netcdf
#if defined PIO_LIB && defined DISTRIBUTE
      USE mod_pio_netcdf
#endif
!
      USE strings_mod, ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ifield, ncid, nfiles, ntide
!
#if defined PIO_LIB && defined DISTRIBUTE
      TYPE (File_desc_t), intent(in) :: pioFile
#endif
      TYPE(T_IO), intent(in) :: S(nfiles)
!
      character (len=4), intent(out) :: Cout(ntide)
!
!  Local variable declarations.
!
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Read in requested 2D field according to IO type.
!-----------------------------------------------------------------------
!
      ncname=Cinfo(ifield,ng)
      SELECT CASE (S(1)%IOtype)
        CASE (io_nf90)
          CALL netcdf_get_svar (ng, model, ncname, 'tide_name',         &
     &                          Cout,                                   &
     &                          ncid = ncid,                            &
     &                          start = (/1, 1/),                       &
     &                          total = (/4, ntide/))
#if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL pio_netcdf_get_svar (ng, model, ncname, 'tide_name',     &
     &                              Cout,                               &
     &                              pioFile = pioFile,                  &
     &                              start = (/1, 1/),                   &
     &                              total = (/4, ntide/))
#endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) S(1)%IOtype
          exit_flag=2
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      IF (Master) THEN
        WRITE (stdout,60) TRIM(Vname(2,ifield)), ng
      END IF
!
  10  FORMAT (' GET_TIDENAME - Illegal input file type, io_type = ',i0, &
     &        /,13x,'Check KeyWord ''INP_LIB'' in ''roms.in''.')
  60  FORMAT (2x,'GET_TIDENAME     - ',a,/,22x,'(Grid = ',i2.2,44x,')')
!
      RETURN
      END SUBROUTINE get_tidename
