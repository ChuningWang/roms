#include "cppdefs.h"

      MODULE isf_output_mod

#ifdef ICESHELF_DYNAM
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2024 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module defines/writes the isf model prognostic variables into  !
!  output NetCDF files.                                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_isf
      USE mod_ncparam
      USE mod_stepping
!
      USE def_var_mod,     ONLY : def_var
# ifdef STATIONS
      USE extract_sta_mod, ONLY : extract_sta2d
# endif
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# ifdef PERFECT_RESTART
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
# endif
      USE uv_rotate_mod,   ONLY : uv_rotate2d
      USE strings_mod,     ONLY : FoundError
!
      implicit none
!
      PUBLIC :: isf_def_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: isf_def_pio
# endif
# ifdef STATIONS
      PUBLIC :: isf_def_station_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: isf_def_station_pio
#  endif
# endif
      PUBLIC :: isf_wrt_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: isf_wrt_pio
# endif
# ifdef STATIONS
      PUBLIC :: isf_wrt_station_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: isf_wrt_station_pio
#  endif
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE isf_def_nf90 (ng, model, ldef, VarOut, S,              &
     &                         t2dgrd, u2dgrd, v2dgrd, p2dgrd,          &
     &                         sr2dgrd, ru2dgrd, rv2dgrd, rp2dgrd)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
      TYPE(T_IO), intent(inout) :: S(Ngrids)
      integer, intent(in), optional, target ::                          &
     &  t2dgrd(:), u2dgrd(:), v2dgrd(:), p2dgrd(:),                    &
     &  sr2dgrd(:), ru2dgrd(:), rv2dgrd(:), rp2dgrd(:)
!
!  Local variable declarations.
!
      logical :: LdefVar, got_var(NV)
!
      integer, parameter :: Natt = 25

      integer :: i, j, ifield, nf, nvd3, nvd4, nvd, status, vtype
      integer :: ic, itrc
!
      integer, pointer :: isfgrd(:)
!
      real(r8) :: Aval(6)
!
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", isf_def_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!
      DEFINE : IF (ldef) THEN
!
!  Set number of dimensions for output variables.
!
# if defined WRITE_WATER && defined MASKING
#  ifdef PERFECT_RESTART
        IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
          nvd3=3
          nvd4=4
        ELSE
          nvd3=2
          nvd4=2
        END IF
#  else
        nvd3=2
        nvd4=2
#  endif
# else
        nvd3=3
        nvd4=4
# endif
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,LEN(Vinfo(1))
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
        IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
          vtype=NF_FRST
        ELSE
          vtype=NF_FOUT
        END IF
# ifdef ICESHEET
!
!-----------------------------------------------------------------------
!  Define isf model state prognostic variables.
!-----------------------------------------------------------------------
!
        DO nf=1,nIsfS
          IF (iSisf(nf).gt.0) THEN
            ifield=iSisf(nf)
            IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
              LdefVar=.TRUE.
#  ifdef PERFECT_RESTART
              nvd=nvd4
#  else
              nvd=nvd3
#  endif
            ELSE
              LdefVar=VarOut(ifield,ng)
              nvd=nvd3
            END IF
            IF (LdefVar) THEN
              Vinfo( 1)=Vname(1,ifield)
              Vinfo( 2)=Vname(2,ifield)
              Vinfo( 3)=Vname(3,ifield)
              Vinfo(14)=Vname(4,ifield)
              Vinfo(16)=Vname(1,idtime)
              Vinfo(21)=Vname(6,ifield)
              Vinfo(22)='coordinates'
              Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
!
              SELECT CASE (nf)
                CASE (isIuis)
#  ifdef PERFECT_RESTART
                  IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
                    isfgrd => rp2dgrd
                  ELSE
                    isfgrd => p2dgrd
                  END IF
#  else
                  isfgrd => p2dgrd
#  endif
#  if defined WRITE_WATER && defined MASKING
#   if !defined WET_DRY && defined PERFECT_RESTART
                  IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
                    Vinfo(24)='_FillValue'
                    Aval(6)=spval
                  END IF
#   else
                  Vinfo(20)='mask_psi'
#   endif
#  endif
                CASE (isIvis)
#  ifdef PERFECT_RESTART
                  IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
                    isfgrd => rp2dgrd
                  ELSE
                    isfgrd => p2dgrd
                  END IF
#  else
                  isfgrd => p2dgrd
#  endif
#  if defined WRITE_WATER && defined MASKING
#   if !defined WET_DRY && defined PERFECT_RESTART
                  IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
                    Vinfo(24)='_FillValue'
                    Aval(6)=spval
                  END IF
#   else
                  Vinfo(20)='mask_psi'
#   endif
#  endif
                CASE DEFAULT
                  isfgrd => t2dgrd
#  if defined WRITE_WATER && defined MASKING
#   if !defined WET_DRY && defined PERFECT_RESTART
                  IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
                    Vinfo(24)='_FillValue'
                    Aval(6)=spval
                  END IF
#   else
                  Vinfo(20)='mask_rho'
#   endif
#  endif
              END SELECT
!
              status=def_var(ng, model, S(ng)%ncid,                     &
     &                       S(ng)%Vid(ifield),                         &
     &                       vtype, nvd, isfgrd, Aval, Vinfo, ncname)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
!
!
!-----------------------------------------------------------------------
!  Define the ice temperature in restart files.  Restart files carry
!  BOTH the exact enthalpy state (above) and its diagnostic
!  temperature (T = E/h, single 2D record): a restart restores the
!  prognostic enthalpy exactly and the temperature consistently.
!-----------------------------------------------------------------------
!
        IF ((S(ng)%ncid.eq.RST(ng)%ncid).and.(idItis.gt.0)) THEN
          DO i=1,Natt
            DO j=1,LEN(Vinfo(1))
              Vinfo(i)(j:j)=' '
            END DO
          END DO
          Aval(6)=0.0_r8
          Vinfo( 1)=Vname(1,idItis)
          Vinfo( 2)=Vname(2,idItis)
          Vinfo( 3)=Vname(3,idItis)
          Vinfo(14)=Vname(4,idItis)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idItis)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idItis,ng),r8)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
#  ifdef PERFECT_RESTART
          isfgrd => sr2dgrd
#  else
          isfgrd => t2dgrd
#  endif
          status=def_var(ng, model, S(ng)%ncid,                       &
     &                   S(ng)%Vid(idItis),                           &
     &                   vtype, nvd3, isfgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError,                          &
     &                   __LINE__, MyFile)) RETURN
        END IF
        IF (S(ng)%ncid.ne.RST(ng)%ncid) THEN
!
!  Reset isfgrid dimension info.
!
          isfgrd => t2dgrd
          nvd=nvd3
!
!  Define isf EASTward velocity at RHO-points.
!
          IF (VarOut(idIuer,ng)) THEN
            Vinfo( 1)=Vname(1,idIuer)
            Vinfo( 2)=Vname(2,idIuer)
            Vinfo( 3)=Vname(3,idIuer)
            Vinfo(14)=Vname(4,idIuer)
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(21)=Vname(6,idIuer)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idIuer,ng),r8)
            status=def_var(ng, model, S(ng)%ncid,                       &
     &                     S(ng)%Vid(idIuer),                           &
     &                     vtype, nvd, isfgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!  Define isf NORTHward velocity at RHO-points.
!
          IF (VarOut(idIvnr,ng)) THEN
            Vinfo( 1)=Vname(1,idIvnr)
            Vinfo( 2)=Vname(2,idIvnr)
            Vinfo( 3)=Vname(3,idIvnr)
            Vinfo(14)=Vname(4,idIvnr)
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(21)=Vname(6,idIvnr)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idIvnr,ng),r8)
            status=def_var(ng, model, S(ng)%ncid,                       &
     &                     S(ng)%Vid(idIvnr),                           &
     &                     vtype, nvd, isfgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
# endif
# if defined ICESHEET || defined ICESHELF_THERMO
!
!-----------------------------------------------------------------------
!  Define isf model internal variables.
!-----------------------------------------------------------------------
!
        IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
          isfgrd => sr2dgrd
        ELSE
          isfgrd => t2dgrd
        END IF
        nvd=nvd3
!
        DO nf=1,nIsfF
          IF (iFisf(nf).gt.0) THEN
            ifield=iFisf(nf)
!
!  The ice temperature in restart files is defined by the dedicated
!  block above; skip it here to avoid a duplicate variable.
!
#  ifdef ICESHEET
            IF ((nf.eq.icItis).and.                                     &
     &          (S(ng)%ncid.eq.RST(ng)%ncid)) THEN
              CYCLE
            END IF
#  endif
            IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
#  ifdef PERFECT_RESTART
              LdefVar=.TRUE.
#  else
              LdefVar=.FALSE.
#  endif
            ELSE
              LdefVar=VarOut(ifield,ng)
            END IF
            IF (LdefVar) THEN
              Vinfo( 1)=Vname(1,ifield)
              Vinfo( 2)=Vname(2,ifield)
              Vinfo( 3)=Vname(3,ifield)
              Vinfo(14)=Vname(4,ifield)
              Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#   if !defined WET_DRY && defined PERFECT_RESTART
              IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
                Vinfo(24)='_FillValue'
                Aval(6)=spval
              END IF
#   else
              Vinfo(20)='mask_rho'
#   endif
#  endif
              Vinfo(21)=Vname(6,ifield)
              Vinfo(22)='coordinates'
              Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
!
              status=def_var(ng, model, S(ng)%ncid,                     &
     &                       S(ng)%Vid(ifield),                         &
     &                       vtype, nvd, isfgrd, Aval, Vinfo, ncname)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Define isf tracer amounts.
!-----------------------------------------------------------------------
!
!  Reset isfgrid dimension info.
!
        isfgrd => t2dgrd
        IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
#  ifdef PERFECT_RESTART
          nvd=nvd4
#  else
          nvd=nvd3
#  endif
        ELSE
          nvd=nvd3
        END IF
        IF (Lisftrc(ng)) THEN
          ifield=idIvol
          IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
            LdefVar=.TRUE.
          ELSE
            LdefVar=VarOut(ifield,ng)
          END IF
          IF (Ldefvar) THEN
            Vinfo( 1)=Vname(1,ifield)
            Vinfo( 2)=Vname(2,ifield)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#   if !defined WET_DRY && defined PERFECT_RESTART
            IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
              Vinfo(24)='_FillValue'
              Aval(6)=spval
            END IF
#   else
            Vinfo(20)='mask_rho'
#   endif
#  endif
            Vinfo(21)=Vname(6,ifield)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            status=def_var(ng, model, S(ng)%ncid,                       &
     &                     S(ng)%Vid(ifield),                           &
     &                     vtype, nvd, isfgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          DO itrc=NAT+1,NT(ng)
            IF (LtracerIsf(itrc,ng)) THEN
              ifield=idItrc(itrc)
              IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
                LdefVar=.TRUE.
              ELSE
                LdefVar=VarOut(ifield,ng)
              END IF
              IF (LdefVar) THEN
                Vinfo( 1)=Vname(1,ifield)
                Vinfo( 2)=Vname(2,ifield)
                Vinfo( 3)=Vname(3,ifield)
                Vinfo(14)=Vname(4,ifield)
                Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
#   if !defined WET_DRY && defined PERFECT_RESTART
                IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
                  Vinfo(24)='_FillValue'
                  Aval(6)=spval
                END IF
#   else
                Vinfo(20)='mask_rho'
#   endif
#  endif
                Vinfo(21)=Vname(6,ifield)
                Vinfo(22)='coordinates'
                Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
                status=def_var(ng, model, S(ng)%ncid,                   &
     &                         S(ng)%Vid(ifield),                       &
     &                         vtype, nvd, isfgrd, Aval, Vinfo, ncname)
                IF (FoundError(exit_flag, NoError, __LINE__, MyFile))   &
     &            RETURN
              END IF
            END IF
          END DO
!
!  Define ambient passive tracer concentration at RHO-points.
!
          DO itrc=NAT+1,NT(ng)
            IF (LtracerIsf(itrc,ng)) THEN
              ifield=idIamb(itrc)
              IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
                LdefVar=.FALSE.
              ELSE
                LdefVar=VarOut(ifield,ng)
              END IF
              IF (LdefVar) THEN
                Vinfo( 1)=Vname(1,ifield)
                Vinfo( 2)=Vname(2,ifield)
                Vinfo( 3)=Vname(3,ifield)
                Vinfo(14)=Vname(4,ifield)
                Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
                Vinfo(20)='mask_rho'
#  endif
                Vinfo(21)=Vname(6,ifield)
                Vinfo(22)='coordinates'
                Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
                status=def_var(ng, model, S(ng)%ncid,                     &
     &                         S(ng)%Vid(ifield),                         &
     &                         vtype, nvd, isfgrd, Aval, Vinfo, ncname)
                IF (FoundError(exit_flag, NoError, __LINE__, MyFile))     &
     &            RETURN
              END IF
            END IF
          END DO
        END IF
# endif
      END IF DEFINE
!
!-----------------------------------------------------------------------
!  Otherwise, check existing output file and prepare for appending
!  data.
!-----------------------------------------------------------------------
!
      QUERY : IF (.not.ldef) THEN
!
!  Initialize locallogical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  ice shelf variables. Get variable IDs.
!
        DO i=1,n_var
# ifdef ICESHEET
          DO nf=1,nIsfS
            IF (iSisf(nf).gt.0) THEN
              ifield=iSisf(nf)
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                got_var(ifield)=.TRUE.
                S(ng)%Vid(ifield)=var_id(i)
                CYCLE
              END IF
            END IF
          END DO
!
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idIuer))) THEN
            got_var(idIuer)=.TRUE.
            S(ng)%Vid(idIuer)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idIvnr))) THEN
            got_var(idIvnr)=.TRUE.
            S(ng)%Vid(idIvnr)=var_id(i)
          END IF
# endif
# if defined ICESHEET || defined ICESHELF_THERMO
!
          DO nf=1,nIsfF
            IF (iFisf(nf).gt.0) THEN
              ifield=iFisf(nf)
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                got_var(ifield)=.TRUE.
                S(ng)%Vid(ifield)=var_id(i)
                CYCLE
              END IF
            END IF
          END DO
!
          IF (Lisftrc(ng)) THEN
            ifield=idIvol
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
              got_var(ifield)=.TRUE.
              S(ng)%Vid(ifield)=var_id(i)
              CYCLE
            END IF
            DO itrc=NAT+1,NT(ng)
              IF (LtracerIsf(itrc, ng)) THEN
                ifield=idItrc(itrc)
                IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                  got_var(ifield)=.TRUE.
                  S(ng)%Vid(ifield)=var_id(i)
                  CYCLE
                END IF
              END IF
            END DO
          END IF
# endif
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
# ifdef ICESHEET
        DO nf=1,nIsfS
          IF (iSisf(nf).gt.0) THEN
            ifield=iSisf(nf)
            IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              RETURN
            END IF
          END IF
        END DO
!
        IF (.not.got_var(idIuer).and.VarOut(idIuer,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idIuer)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idIvnr).and.VarOut(idIvnr,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idIvnr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined ICESHEET || defined ICESHELF_THERMO
!
        DO nf=1,nIsfF
          IF (iFisf(nf).gt.0) THEN
            ifield=iFisf(nf)
!
!  The ice temperature in restart files is defined by the dedicated
!  block above; skip it here to avoid a duplicate variable.
!
#  ifdef ICESHEET
            IF ((nf.eq.icItis).and.                                     &
     &          (S(ng)%ncid.eq.RST(ng)%ncid)) THEN
              CYCLE
            END IF
#  endif
            IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
#  ifdef PERFECT_RESTART
              LdefVar=.TRUE.
#  else
              LdefVar=.FALSE.
#  endif
            ELSE
              LdefVar=VarOut(ifield,ng)
            END IF
            IF (.not.got_var(ifield).and.LdefVar) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              RETURN
            END IF
          END IF
        END DO
!
        IF (Lisftrc(ng)) THEN
          ifield=idIvol
          IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),    &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          DO itrc=NAT+1,NT(ng)
            IF (LtracerIsf(itrc, ng)) THEN
              ifield=idItrc(itrc)
              IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                RETURN
              END IF
            END IF
          END DO
        END IF
# endif
      END IF QUERY
!
  10  FORMAT (/,' ISF_DEF_NF90 - unable to find variable: ',a,2x,       &
     &        ' in output NetCDF file: ',a)
      RETURN
      END SUBROUTINE isf_def_nf90
!
!***********************************************************************
     SUBROUTINE isf_wrt_nf90 (ng, model, tile,                          &
    &                         LBi, UBi, LBj, UBj,                       &
    &                         VarOut, S)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      logical, intent(in) :: VarOut(NV,Ngrids)
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: LwrtVar
!
      integer :: gfactor, gtype, status
      integer :: ifield, ifld, ic, itrc, no
!
      real(dp) :: scale
!
      real(r8), pointer :: isfField2d(:,:)
      real(r8), pointer :: isfField3d(:,:,:)
      real(r8), pointer :: isfMask(:,:)
!
      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
      real(r8), allocatable, target :: Tmp2d(:,:)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", isf_wrt_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out isf prognostic variables into specified output NetCDF file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
# if defined WRITE_WATER && defined MASKING
#  ifdef PERFECT_RESTART
      IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
        gfactor=1
      ELSE
        gfactor=-1
      END IF
#  else
      gfactor=-1
#  endif
# else
      gfactor=1
# endif
# ifdef ICESHEET
!
!-----------------------------------------------------------------------
!  write out ice shelf model state prognostic variables.
!-----------------------------------------------------------------------
!
      DO ifld=1,nIsfS
        IF (iSisf(ifld).gt.0) THEN
          ifield=iSisf(ifld)
          IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
            LwrtVar=.TRUE.
          ELSE
            LwrtVar=VarOut(ifield,ng)
          END IF
!
          IF (LwrtVar) THEN
            IF ((model.eq.iNLM).and.                                    &
                ((S(ng)%ncid.eq.HIS(ng)%ncid).or.                       &
     &           (S(ng)%ncid.eq.QCK(ng)%ncid))) THEN
              isfField2d => ISF(ng) % Si(LBi:UBi,LBj:UBj,lsnew(ng),ifld)
            ELSE IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
#  ifdef PERFECT_RESTART
              isfField3d => ISF(ng) % Si(LBi:UBi,LBj:UBj,1:2,ifld)
#  else
              isfField2d => ISF(ng) % Si(LBi:UBi,LBj:UBj,lsnew(ng),ifld)
#  endif
#  ifdef AVERAGES
            ELSE IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              isfField2d => ISF_SAVG(ifld,ng) % var(LBi:UBi,LBj:UBj)
#  endif
            END IF
!
!  The ice state must be masked by the STATIC land/sea mask only:
!  the ocean wet/dry mask (*_full) is zero over grounded-ice cells
!  (Dcrit water film), which would fill the restart file with
!  _FillValue exactly where the ice is and destroy the state on
!  the restart read.
!
            SELECT CASE (ifld)
              CASE (isIuis)
                IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
#  ifdef PERFECT_RESTART
                  gtype=gfactor*p3dvar
#  else
                  gtype=gfactor*p2dvar
#  endif
                ELSE
                  gtype=gfactor*p2dvar
                END IF
#  ifdef MASKING
                isfMask => GRID(ng) % pmask
#  endif
              CASE (isIvis)
                IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
#  ifdef PERFECT_RESTART
                  gtype=gfactor*p3dvar
#  else
                  gtype=gfactor*p2dvar
#  endif
                ELSE
                  gtype=gfactor*p2dvar
                END IF
#  ifdef MASKING
                isfMask => GRID(ng) % pmask
#  endif
              CASE DEFAULT
                IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
#  ifdef PERFECT_RESTART
                  gtype=gfactor*r3dvar
#  else
                  gtype=gfactor*r2dvar
#  endif
                ELSE
                  gtype=gfactor*r2dvar
                END IF
#  ifdef MASKING
                isfMask => GRID(ng) % rmask
#  endif
            END SELECT
            scale=1.0_dp
!
#  ifdef PERFECT_RESTART
            IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
              status=nf_fwrite3d(ng, model, S(ng)%ncid, ifield,         &
     &                           S(ng)%Vid(ifield),                     &
     &                           S(ng)%Rindex, gtype,                   &
     &                           LBi, UBi, LBj, UBj, 1, 2, scale,       &
#   ifdef MASKING
     &                           isfMask,                               &
#   endif
     &                           isfField3d)
            ELSE
              status=nf_fwrite2d(ng, model, S(ng)%ncid, ifield,         &
     &                           S(ng)%Vid(ifield),                     &
     &                           S(ng)%Rindex, gtype,                   &
     &                           LBi, UBi, LBj, UBj, scale,             &
#   ifdef MASKING
     &                           isfMask,                               &
#   endif
     &                           isfField2d)
            END IF
#  else
            status=nf_fwrite2d(ng, model, S(ng)%ncid, ifield,           &
     &                         S(ng)%Vid(ifield),                       &
     &                         S(ng)%Rindex, gtype,                     &
     &                         LBi, UBi, LBj, UBj, scale,               &
#   ifdef MASKING
     &                         isfMask,                                 &
#   endif
     &                         isfField2d)
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      S(ng)%Rindex,               &
     &                                      TRIM(S(ng)%name)
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          END IF
        END IF
      END DO
!
!
!-----------------------------------------------------------------------
!  Write out the ice temperature (T = E/h) into restart files so a
!  restart restores both the exact enthalpy state (above) and its
!  diagnostic temperature.  Masked by the STATIC land/sea mask, like
!  the state fields above.
!-----------------------------------------------------------------------
!
        IF ((S(ng)%ncid.eq.RST(ng)%ncid).and.(idItis.gt.0)) THEN
          IF (.not.allocated(Tmp2d)) THEN
            allocate (Tmp2d(LBi:UBi,LBj:UBj))
          END IF
          WHERE (ISF(ng)%Si(LBi:UBi,LBj:UBj,lsnew(ng),isIhis)         &
     &           .gt.0.0_r8)
            Tmp2d=ISF(ng)%Si(LBi:UBi,LBj:UBj,lsnew(ng),isIeis)/       &
     &            ISF(ng)%Si(LBi:UBi,LBj:UBj,lsnew(ng),isIhis)
          ELSEWHERE
            Tmp2d=0.0_r8
          END WHERE
          scale=1.0_dp
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idItis,           &
     &                       S(ng)%Vid(idItis),                       &
     &                       S(ng)%Rindex, gfactor*r2dvar,            &
     &                       LBi, UBi, LBj, UBj, scale,               &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                        &
#  endif
     &                       Tmp2d)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,idItis)),      &
     &                                    S(ng)%Rindex,               &
     &                                    TRIM(S(ng)%name)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!  Write isf EASTward and NORTHward velocity component (m/s) at
!  RHO-points.
!
      IF (VarOut(idIuer,ng).and.VarOut(idIvnr,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
          Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
          Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
#  ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
#  endif
     &                    ISF(ng) % Si(:,:,lsnew(ng),isIuis),           &
     &                    ISF(ng) % Si(:,:,lsnew(ng),isIvis),           &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idIuer,               &
     &                     S(ng)%Vid(idIuer),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     Ur2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idIuer)),          &
     &                                  S(ng)%Rindex, TRIM(S(ng)%name)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idIvnr,               &
     &                     S(ng)%Vid(idIvnr),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     Vr2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idIvnr)),          &
     &                                  S(ng)%Rindex, TRIM(S(ng)%name)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
# endif
# if defined ICESHEET || defined ICESHELF_THERMO
!
!-----------------------------------------------------------------------
!  write out isf model internal variables.
!-----------------------------------------------------------------------
!
      DO ifld=1,nIsfF
        IF (iFisf(ifld).gt.0) THEN
          ifield=iFisf(ifld)
!
!  The ice temperature in restart files is written by the dedicated
!  block above; skip it here to avoid a duplicate write.
!
#  ifdef ICESHEET
          IF ((ifld.eq.icItis).and.                                     &
     &        (S(ng)%ncid.eq.RST(ng)%ncid)) THEN
            CYCLE
          END IF
#  endif
          IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
#  ifdef PERFECT_RESTART
            LwrtVar=.TRUE.
#  else
            LwrtVar=.FALSE.
#  endif
          ELSE
            LwrtVar=VarOut(ifield,ng)
          END IF
!
          IF (LwrtVar) THEN
            IF ((model.eq.iNLM).and.                                    &
                ((S(ng)%ncid.eq.HIS(ng)%ncid).or.                       &
     &           (S(ng)%ncid.eq.QCK(ng)%ncid).or.                       &
     &           (S(ng)%ncid.eq.RST(ng)%ncid))) THEN
              isfField2d => ISF(ng) % Fi(LBi:UBi,LBj:UBj,ifld)
#  ifdef AVERAGES
            ELSE IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              isfField2d => ISF_FAVG(ifld,ng) % var(LBi:UBi,LBj:UBj)
#  endif
            END IF
!
            scale=1.0_dp
            gtype=gfactor*r2dvar
            status=nf_fwrite2d(ng, model, S(ng)%ncid, ifield,           &
     &                         S(ng)%Vid(ifield),                       &
     &                         S(ng)%Rindex, gtype,                     &
     &                         LBi, UBi, LBj, UBj, scale,               &
#  ifdef MASKING
     &                         GRID(ng) % rmask_full,                   &
#  endif
     &                         isfField2d)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      S(ng)%Rindex,               &
     &                                      TRIM(S(ng)%name)
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          END IF
        END IF
      END DO
#  if defined ICESHELF_THERMO && !defined ICESHELF_NO_SBL
#   ifdef PERFECT_RESTART
!
!  Write out SBL weights (Oi) to restart file.
!
      IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
        DO no=1,nIsfO
          IF (iOisf(no).gt.0) THEN
            ifield=iOisf(no)
            isfField3d => ISF(ng) % Oi(LBi:UBi,LBj:UBj,:,no)
            scale=1.0_dp
            gtype=gfactor*r3dvar
            status=nf_fwrite3d(ng, model, S(ng)%ncid, ifield,           &
     &                         S(ng)%Vid(ifield),                       &
     &                         S(ng)%Rindex, gtype,                     &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#   ifdef MASKING
     &                         GRID(ng) % rmask_full,                   &
#   endif
     &                         isfField3d)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      S(ng)%Rindex,               &
     &                                      TRIM(S(ng)%name)
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          END IF
        END DO
      END IF
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  write out isf tracer amount at RHO-points.
!-----------------------------------------------------------------------
!
      IF (Lisftrc(ng)) THEN
        ifield=idIvol
        IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
          LwrtVar=.TRUE.
        ELSE
          LwrtVar=VarOut(ifield,ng)
        END IF
!
        IF (LwrtVar) THEN
          IF ((model.eq.iNLM).and.                                  &
              ((S(ng)%ncid.eq.HIS(ng)%ncid).or.                     &
     &         (S(ng)%ncid.eq.QCK(ng)%ncid))) THEN
            isfField2d => ISF(ng) % Vi(LBi:UBi,LBj:UBj,lsnew(ng))
          ELSE IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
#  ifdef PERFECT_RESTART
            isfField3d => ISF(ng) % Vi(LBi:UBi,LBj:UBj,1:2)
#  else
            isfField2d => ISF(ng) % Vi(LBi:UBi,LBj:UBj,lsnew(ng))
#  endif
#  ifdef AVERAGES
          ELSE IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            isfField2d => ISF_VAVG(ng) % var(LBi:UBi,LBj:UBj)
#  endif
          END IF
!
          scale=1.0_dp
#  ifdef PERFECT_RESTART
          IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
            gtype=gfactor*r3dvar
          ELSE
            gtype=gfactor*r2dvar
          END IF
#  else
          gtype=gfactor*r2dvar
#  endif
#  ifdef PERFECT_RESTART
          IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
            status=nf_fwrite3d(ng, model, S(ng)%ncid, ifield,           &
     &                         S(ng)%Vid(ifield),                       &
     &                         S(ng)%Rindex, gtype,                     &
     &                         LBi, UBi, LBj, UBj, 1, 2, scale,         &
#   ifdef MASKING
     &                         GRID(ng) % rmask_full,                   &
#   endif
     &                         isfField3d)
          ELSE
            status=nf_fwrite2d(ng, model, S(ng)%ncid, ifield,           &
     &                         S(ng)%Vid(ifield),                       &
     &                         S(ng)%Rindex, gtype,                     &
     &                         LBi, UBi, LBj, UBj, scale,               &
#   ifdef MASKING
     &                         GRID(ng) % rmask_full,                   &
#   endif
     &                         isfField2d)
          END IF
#  else
          status=nf_fwrite2d(ng, model, S(ng)%ncid, ifield,             &
     &                       S(ng)%Vid(ifield),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask_full,                     &
#   endif
     &                       isfField2d)
#  endif
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),        &
     &                                    S(ng)%Rindex,                 &
     &                                    TRIM(S(ng)%name)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
        ic=0
        DO itrc=NAT+1,NT(ng)
          IF (LtracerIsf(itrc, ng)) THEN
            ic=ic+1
            ifield=idItrc(itrc)
            IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
              LwrtVar=.TRUE.
            ELSE
              LwrtVar=VarOut(ifield,ng)
            END IF
!
            IF (LwrtVar) THEN
              IF ((model.eq.iNLM).and.                                  &
                  ((S(ng)%ncid.eq.HIS(ng)%ncid).or.                     &
     &             (S(ng)%ncid.eq.QCK(ng)%ncid))) THEN
                isfField2d => ISF(ng) % Ti(LBi:UBi,LBj:UBj,lsnew(ng),ic)
              ELSE IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
#  ifdef PERFECT_RESTART
                isfField3d => ISF(ng) % Ti(LBi:UBi,LBj:UBj,1:2,ic)
#  else
                isfField2d => ISF(ng) % Ti(LBi:UBi,LBj:UBj,lsnew(ng),ic)
#  endif
#  ifdef AVERAGES
              ELSE IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
                isfField2d => ISF_TAVG(itrc,ng) % var(LBi:UBi,LBj:UBj)
#  endif
              END IF
!
              scale=1.0_dp
#  ifdef PERFECT_RESTART
              IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
                gtype=gfactor*r3dvar
              ELSE
                gtype=gfactor*r2dvar
              END IF
#  else
              gtype=gfactor*r2dvar
#  endif
#  ifdef PERFECT_RESTART
              IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
                status=nf_fwrite3d(ng, model, S(ng)%ncid, ifield,       &
     &                             S(ng)%Vid(ifield),                   &
     &                             S(ng)%Rindex, gtype,                 &
     &                             LBi, UBi, LBj, UBj, 1, 2, scale,     &
#   ifdef MASKING
     &                             GRID(ng) % rmask_full,               &
#   endif
     &                             isfField3d)
              ELSE
                status=nf_fwrite2d(ng, model, S(ng)%ncid, ifield,       &
     &                             S(ng)%Vid(ifield),                   &
     &                             S(ng)%Rindex, gtype,                 &
     &                             LBi, UBi, LBj, UBj, scale,           &
#   ifdef MASKING
     &                             GRID(ng) % rmask_full,               &
#   endif
     &                             isfField2d)
              END IF
#  else
              status=nf_fwrite2d(ng, model, S(ng)%ncid, ifield,         &
     &                           S(ng)%Vid(ifield),                     &
     &                           S(ng)%Rindex, gtype,                   &
     &                           LBi, UBi, LBj, UBj, scale,             &
#   ifdef MASKING
     &                           GRID(ng) % rmask_full,                 &
#   endif
     &                           isfField2d)
#  endif
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),    &
     &                                        S(ng)%Rindex,             &
     &                                        TRIM(S(ng)%name)
                exit_flag=3
                ioerror=status
                RETURN
              END IF
            END IF
          END IF
        END DO
!
!  Write out ambient passive tracer concentration at RHO-points.
!
        ic=0
        DO itrc=NAT+1,NT(ng)
          IF (LtracerIsf(itrc, ng)) THEN
            ic=ic+1
            ifield=idIamb(itrc)
            IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
              LwrtVar=.FALSE.
            ELSE
              LwrtVar=VarOut(ifield,ng)
            END IF
!
            IF (LwrtVar) THEN
              IF ((model.eq.iNLM).and.                                  &
                  ((S(ng)%ncid.eq.HIS(ng)%ncid).or.                     &
     &             (S(ng)%ncid.eq.QCK(ng)%ncid))) THEN
                isfField2d => ISF(ng) % Ai(LBi:UBi,LBj:UBj,ic)
#  ifdef AVERAGES
              ELSE IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
                isfField2d => ISF_AAVG(itrc,ng) % var(LBi:UBi,LBj:UBj)
#  endif
              END IF
!
              scale=1.0_dp
              gtype=gfactor*r2dvar
              status=nf_fwrite2d(ng, model, S(ng)%ncid, ifield,         &
     &                           S(ng)%Vid(ifield),                     &
     &                           S(ng)%Rindex, gtype,                   &
     &                           LBi, UBi, LBj, UBj, scale,             &
#   ifdef MASKING
     &                           GRID(ng) % rmask_full,                 &
#   endif
     &                           isfField2d)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),    &
     &                                        S(ng)%Rindex,             &
     &                                        TRIM(S(ng)%name)
                exit_flag=3
                ioerror=status
                RETURN
              END IF
            END IF
          END IF
        END DO
      END IF
# endif
!
  10  FORMAT (/," ISF_WRT_NF90 - error while writing variable '",a,     &
     &        "', time record = ",i0,/,11x,'into file: ',a)
!
      RETURN
      END SUBROUTINE isf_wrt_nf90

# ifdef STATIONS
!
!***********************************************************************
      SUBROUTINE isf_def_station_nf90 (ng, model, ldef, VarOut, S,      &
     &                                 pgrd)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: pgrd(:)
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer, parameter :: Natt = 25

      integer :: i, ifield, j, nf, status, ic, itrc
!
      real(r8) :: Aval(6)
!
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", isf_def_station_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define prognostic isf variables into specified output NetCDF.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!
!  Initialize local information variable arrays.
!
      DO i=1,Natt
        DO j=1,LEN(Vinfo(1))
          Vinfo(i)(j:j)=' '
        END DO
      END DO
      DO i=1,6
        Aval(i)=0.0_r8
      END DO
!
!-----------------------------------------------------------------------
!  Define isf model state prognostic variables.
!-----------------------------------------------------------------------
!
      DEFINE : IF (ldef) THEN

#  ifdef ICESHEET
        DO nf=1,nIsfS
          IF (iSisf(nf).gt.0) THEN
            ifield=iSisf(nf)
            IF (VarOut(ifield,ng)) THEN
              Vinfo( 1)=Vname(1,ifield)
              Vinfo( 2)=Vname(2,ifield)
              Vinfo( 3)=Vname(3,ifield)
              Vinfo(14)=Vname(4,ifield)
              Vinfo(16)=Vname(1,idtime)
!
              status=def_var(ng, model, S(ng)%ncid,                     &
     &                       S(ng)%Vid(ifield),                         &
     &                       NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,     &
     &                       SetFillVal = .TRUE.,                       &
     &                       SetParAccess = .TRUE.)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
!
        IF (VarOut(idIuer,ng)) THEN
          Vinfo( 1)=Vname(1,idIuer)
          Vinfo( 2)=Vname(2,idIuer)
          Vinfo( 3)=Vname(3,idIuer)
          Vinfo(14)=Vname(4,idIuer)
          Vinfo(16)=Vname(1,idtime)
!
          status=def_var(ng, model, S(ng)%ncid,                         &
     &                   S(ng)%Vid(idIuer),                             &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (VarOut(idIvnr,ng)) THEN
          Vinfo( 1)=Vname(1,idIvnr)
          Vinfo( 2)=Vname(2,idIvnr)
          Vinfo( 3)=Vname(3,idIvnr)
          Vinfo(14)=Vname(4,idIvnr)
          Vinfo(16)=Vname(1,idtime)
!
          status=def_var(ng, model, S(ng)%ncid,                         &
     &                   S(ng)%Vid(idIvnr),                             &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  if defined ICESHEET || defined ICESHELF_THERMO
!
!-----------------------------------------------------------------------
!  Define isf model internal variables.
!-----------------------------------------------------------------------
!
        DO nf=1,nIsfF
          IF (iFisf(nf).gt.0) THEN
            ifield=iFisf(nf)
            IF (VarOut(ifield,ng)) THEN
              Vinfo( 1)=Vname(1,ifield)
              Vinfo( 2)=Vname(2,ifield)
              Vinfo( 3)=Vname(3,ifield)
              Vinfo(14)=Vname(4,ifield)
              Vinfo(16)=Vname(1,idtime)
!
              status=def_var(ng, model, S(ng)%ncid,                     &
     &                       S(ng)%Vid(ifield),                         &
     &                       NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,     &
     &                       SetFillVal = .TRUE.,                       &
     &                       SetParAccess = .TRUE.)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Define isf tracer amounts.
!-----------------------------------------------------------------------
!
        IF (Lisftrc(ng)) THEN
          ifield=idIvol
          IF (VarOut(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            Vinfo( 2)=Vname(2,ifield)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
!
            status=def_var(ng, model, S(ng)%ncid,                       &
     &                     S(ng)%Vid(ifield),                           &
     &                     NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,       &
     &                     SetFillVal = .TRUE.,                         &
     &                     SetParAccess = .TRUE.)
            IF (FoundError(exit_flag, NoError,                          &
     &                     __LINE__, MyFile)) RETURN
          END IF
          DO itrc=NAT+1,NT(ng)
            IF (LtracerIsf(itrc, ng)) THEN
              ifield=idItrc(itrc)
              IF (VarOut(ifield,ng)) THEN
                Vinfo( 1)=Vname(1,ifield)
                Vinfo( 2)=Vname(2,ifield)
                Vinfo( 3)=Vname(3,ifield)
                Vinfo(14)=Vname(4,ifield)
                Vinfo(16)=Vname(1,idtime)
!
                status=def_var(ng, model, S(ng)%ncid,                   &
     &                         S(ng)%Vid(ifield),                       &
     &                         NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,   &
     &                         SetFillVal = .TRUE.,                     &
     &                         SetParAccess = .TRUE.)
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
              END IF
            END IF
          END DO
!
!  Define ambient passive tracer concentration.
!
          DO itrc=NAT+1,NT(ng)
            IF (LtracerIsf(itrc, ng)) THEN
              ifield=idIamb(itrc)
              IF (VarOut(ifield,ng)) THEN
                Vinfo( 1)=Vname(1,ifield)
                Vinfo( 2)=Vname(2,ifield)
                Vinfo( 3)=Vname(3,ifield)
                Vinfo(14)=Vname(4,ifield)
                Vinfo(16)=Vname(1,idtime)
!
                status=def_var(ng, model, S(ng)%ncid,                   &
     &                         S(ng)%Vid(ifield),                       &
     &                         NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,   &
     &                         SetFillVal = .TRUE.,                     &
     &                         SetParAccess = .TRUE.)
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
              END IF
            END IF
          END DO
        END IF
#  endif
      END IF DEFINE
!
!-----------------------------------------------------------------------
!  Otherwise, check existing output file and prepare for appending
!  data.
!-----------------------------------------------------------------------
!
      QUERY : IF (.not.ldef) THEN
!
!  Initialize locallogical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  ice shelf variables. Get variable IDs.
!
        DO i=1,n_var
#  ifdef ICESHEET
          DO nf=1,nIsfS
            IF (iSisf(nf).gt.0) THEN
              ifield=iSisf(nf)
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                got_var(ifield)=.TRUE.
                S(ng)%Vid(ifield)=var_id(i)
                CYCLE
              END IF
            END IF
          END DO
!
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idIuer))) THEN
            got_var(idIuer)=.TRUE.
            S(ng)%Vid(idIuer)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idIvnr))) THEN
            got_var(idIvnr)=.TRUE.
            S(ng)%Vid(idIvnr)=var_id(i)
          END IF
#  endif
#  if defined ICESHEET || defined ICESHELF_THERMO
!
          DO nf=1,nIsfF
            IF (iFisf(nf).gt.0) THEN
              ifield=iFisf(nf)
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                got_var(ifield)=.TRUE.
                S(ng)%Vid(ifield)=var_id(i)
                CYCLE
              END IF
            END IF
          END DO
!
          IF (Lisftrc(ng)) THEN
            ifield=idIvol
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
              got_var(ifield)=.TRUE.
              S(ng)%Vid(ifield)=var_id(i)
              CYCLE
            END IF
            DO itrc=NAT+1,NT(ng)
              IF (LtracerIsf(itrc, ng)) THEN
                ifield=idItrc(itrc)
                IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                  got_var(ifield)=.TRUE.
                  S(ng)%Vid(ifield)=var_id(i)
                  CYCLE
                END IF
              END IF
            END DO
          END IF
#  endif
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
#  ifdef ICESHEET
        DO nf=1,nIsfS
          IF (iSisf(nf).gt.0) THEN
            ifield=iSisf(nf)
            IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              RETURN
            END IF
          END IF
        END DO
!
        IF (.not.got_var(idIuer).and.VarOut(idIuer,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idIuer)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idIvnr).and.VarOut(idIvnr,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idIvnr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined ICESHEET || defined ICESHELF_THERMO
!
        DO nf=1,nIsfF
          IF (iFisf(nf).gt.0) THEN
            ifield=iFisf(nf)
            IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              RETURN
            END IF
          END IF
        END DO
!
        IF (Lisftrc(ng)) THEN
          ifield=idIvol
          IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),    &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          DO itrc=NAT+1,NT(ng)
            IF (LtracerIsf(itrc, ng)) THEN
              ifield=idItrc(itrc)
              IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                RETURN
              END IF
            END IF
          END DO
        END IF
#  endif
      END IF QUERY
!
  10  FORMAT (/,' ISF_DEF_STATION_NF90 - unable to find variable: ',    &
     &        a,2x,' in output NetCDF file: ',a)
      RETURN
      END SUBROUTINE isf_def_station_nf90
!
!***********************************************************************
      SUBROUTINE isf_wrt_station_nf90 (ng, model, tile,                 &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 VarOut, S)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: Cgrid
!
      integer :: i, ifield, ifld, ic, itrc
!
      real(dp) :: scale
!
      real(r8), dimension(Nstation(ng)) :: Xpos, Ypos, psta
!
      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", isf_wrt_station_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out isf prognostic variables into specified output NetCDF file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set switch to extract station data at native C-grid position (TRUE)
!  or at RHO-points (FALSE).
!
#  ifdef STATIONS_CGRID
      Cgrid=.TRUE.
#  else
      Cgrid=.FALSE.
#  endif
!
!  Set positions for generic extraction routine.
!
      DO i=1,Nstation(ng)
        Xpos(i)=SCALARS(ng)%SposX(i)
        Ypos(i)=SCALARS(ng)%SposY(i)
      END DO
#  ifdef ICESHEET
!
!-----------------------------------------------------------------------
!  write out ice shelf model state prognostic variables.
!-----------------------------------------------------------------------
!
      DO ifld=1,nIsfS
        IF (iSisf(ifld).gt.0) THEN
          ifield=iSisf(ifld)
          IF (VarOut(ifield,ng)) THEN
            scale=1.0_dp
            CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          scale, ISF(ng)%Si(:,:,lsnew(ng),ifld),  &
     &                          Nstation(ng), Xpos, Ypos, psta)
!
            CALL netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,ifield)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            ncid = S(ng)%ncid,                    &
     &                            varid = S(ng)%Vid(ifield))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
      END DO
!
!  Write isf EASTward and NORTHward velocity component (m/s) at
!  RHO-points.
!
      IF (VarOut(idIuer,ng).and.VarOut(idIvnr,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
          Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
          Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
#   ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
#   endif
     &                    ISF(ng) % Si(:,:,lsnew(ng),isIuis),           &
     &                    ISF(ng) % Si(:,:,lsnew(ng),isIvis),           &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idIuer, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Ur2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idIuer)), psta,              &
     &                        (/1,S(ng)%Rindex/),                       &
     &                        (/Nstation(ng),1/),                       &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idIuer))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL extract_sta2d (ng, model, Cgrid, idIvnr, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Vr2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idIvnr)), psta,              &
     &                        (/1,S(ng)%Rindex/),                       &
     &                        (/Nstation(ng),1/),                       &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idIvnr))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
#  endif
#  if defined ICESHEET || defined ICESHELF_THERMO
!
!-----------------------------------------------------------------------
!  write out isf model internal variables.
!-----------------------------------------------------------------------
!
      DO ifld=1,nIsfF
        IF (iFisf(ifld).gt.0) THEN
          ifield=iFisf(ifld)
          IF (VarOut(ifield,ng)) THEN
            scale=1.0_dp
            CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          scale, ISF(ng)%Fi(:,:,ifld),            &
     &                          Nstation(ng), Xpos, Ypos, psta)
!
            CALL netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,ifield)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            ncid = S(ng)%ncid,                    &
     &                            varid = S(ng)%Vid(ifield))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  write out isf tracer amount.
!-----------------------------------------------------------------------
!
      IF (Lisftrc(ng)) THEN
        ifield=idIvol
        IF (VarOut(ifield,ng)) THEN
          scale=1.0_dp
          CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        scale, ISF(ng)%Vi(:,:,lsnew(ng)),         &
     &                        Nstation(ng), Xpos, Ypos, psta)
!
          CALL netcdf_put_fvar (ng, model, S(ng)%name,                  &
     &                          TRIM(Vname(1,ifield)), psta,            &
     &                          (/1,S(ng)%Rindex/),                     &
     &                          (/Nstation(ng),1/),                     &
     &                          ncid = S(ng)%ncid,                      &
     &                          varid = S(ng)%Vid(ifield))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile))         &
     &      RETURN
        END IF
        ic=0
        DO itrc=NAT+1,NT(ng)
          IF (LtracerIsf(itrc, ng)) THEN
            ic=ic+1
            ifield=idItrc(itrc)
            IF (VarOut(ifield,ng)) THEN
              scale=1.0_dp
              CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            scale, ISF(ng)%Ti(:,:,lsnew(ng),ic),  &
     &                            Nstation(ng), Xpos, Ypos, psta)
!
              CALL netcdf_put_fvar (ng, model, S(ng)%name,              &
     &                              TRIM(Vname(1,ifield)), psta,        &
     &                              (/1,S(ng)%Rindex/),                 &
     &                              (/Nstation(ng),1/),                 &
     &                              ncid = S(ng)%ncid,                  &
     &                              varid = S(ng)%Vid(ifield))
              IF (FoundError(exit_flag, NoError, __LINE__, MyFile))     &
     &          RETURN
            END IF
          END IF
        END DO
!
!  Write out ambient passive tracer concentration.
!
        ic=0
        DO itrc=NAT+1,NT(ng)
          IF (LtracerIsf(itrc, ng)) THEN
            ic=ic+1
            ifield=idIamb(itrc)
            IF (VarOut(ifield,ng)) THEN
              scale=1.0_dp
              CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            scale, ISF(ng)%Ai(:,:,ic),            &
     &                            Nstation(ng), Xpos, Ypos, psta)
!
              CALL netcdf_put_fvar (ng, model, S(ng)%name,              &
     &                              TRIM(Vname(1,ifield)), psta,        &
     &                              (/1,S(ng)%Rindex/),                 &
     &                              (/Nstation(ng),1/),                 &
     &                              ncid = S(ng)%ncid,                  &
     &                              varid = S(ng)%Vid(ifield))
              IF (FoundError(exit_flag, NoError, __LINE__, MyFile))     &
     &          RETURN
            END IF
          END IF
        END DO
      END IF
#  endif
!
  10  FORMAT (/,' ISF_WRT_STATION_NF90 - error while writing ',         &
     &        "variable '",a,"', time record = ",i0,/,11x,              &
     &        'into file: ',a)
!
      RETURN
      END SUBROUTINE isf_wrt_station_nf90
# endif

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE isf_def_pio (ng, model, ldef, VarOut, S,               &
     &                        t2dgrd, u2dgrd, v2dgrd, p2dgrd,           &
     &                        sr2dgrd, ru2dgrd, rv2dgrd, rp2dgrd)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
      TYPE(T_IO), intent(inout) :: S(Ngrids)
      integer, intent(in), optional, target ::                          &
     &  t2dgrd(:), u2dgrd(:), v2dgrd(:), p2dgrd(:),                    &
     &  sr2dgrd(:), ru2dgrd(:), rv2dgrd(:), rp2dgrd(:)
!
!  Local variable declarations.
!
      logical :: LdefVar, got_var(NV)
!
      integer, parameter :: Natt = 25
!
      integer :: FileH, i, j, ifield, nf, nvd3, nvd4, nvd, status, vtype
      integer :: ic, itrc
!
      integer, pointer :: isfgrd(:)
!
      real(r8) :: Aval(6)
!
      character (len=256) :: ncname
      character (len=MaxLen) :: Vinfo(Natt)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", isf_def_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
      FileH=ABS(S(ng)%pioFile%fh)
!
      DEFINE : IF (ldef) THEN
!
!  Set number of dimensions for output variables.
!
#  if defined WRITE_WATER && defined MASKING
#   ifdef PERFECT_RESTART
        IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
          nvd3=3
          nvd4=4
        ELSE
          nvd3=2
          nvd4=2
        END IF
#   else
        nvd3=2
        nvd4=2
#   endif
#  else
        nvd3=3
        nvd4=4
#  endif
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,LEN(Vinfo(1))
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
        IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
          vtype=PIO_FRST
        ELSE
          vtype=PIO_FOUT
        END IF
#  ifdef ICESHEET
!
!-----------------------------------------------------------------------
!  Define isf model state prognostic variables.
!-----------------------------------------------------------------------
!
        DO nf=1,nIsfS
          IF (iSisf(nf).gt.0) THEN
            ifield=iSisf(nf)
            IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
              LdefVar=.TRUE.
#   ifdef PERFECT_RESTART
              nvd=nvd4
#   else
              nvd=nvd3
#   endif
            ELSE
              LdefVar=VarOut(ifield,ng)
              nvd=nvd3
            END IF
            IF (LdefVar) THEN
              Vinfo( 1)=Vname(1,ifield)
              Vinfo( 2)=Vname(2,ifield)
              Vinfo( 3)=Vname(3,ifield)
              Vinfo(14)=Vname(4,ifield)
              Vinfo(16)=Vname(1,idtime)
              Vinfo(21)=Vname(6,ifield)
              Vinfo(22)='coordinates'
              Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
              S(ng)%pioVar(ifield)%dkind=PIO_FOUT
!
              SELECT CASE (nf)
                CASE (isIuis)
#   ifdef PERFECT_RESTART
                  IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
                    isfgrd => rp2dgrd
                  ELSE
                    isfgrd => p2dgrd
                  END IF
#   else
                  isfgrd => p2dgrd
#   endif
#   if defined WRITE_WATER && defined MASKING
#    if !defined WET_DRY && defined PERFECT_RESTART
                  IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
                    Vinfo(24)='_FillValue'
                    Aval(6)=spval
                  END IF
#    else
                  Vinfo(20)='mask_psi'
#    endif
#   endif
                  S(ng)%pioVar(ifield)%gtype=p2dvar
                CASE (isIvis)
#   ifdef PERFECT_RESTART
                  IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
                    isfgrd => rp2dgrd
                  ELSE
                    isfgrd => p2dgrd
                  END IF
#   else
                  isfgrd => p2dgrd
#   endif
#   if defined WRITE_WATER && defined MASKING
#    if !defined WET_DRY && defined PERFECT_RESTART
                  IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
                    Vinfo(24)='_FillValue'
                    Aval(6)=spval
                  END IF
#    else
                  Vinfo(20)='mask_psi'
#    endif
#   endif
                  S(ng)%pioVar(ifield)%gtype=p2dvar
                CASE DEFAULT
                  isfgrd => t2dgrd
#   if defined WRITE_WATER && defined MASKING
#    if !defined WET_DRY && defined PERFECT_RESTART
                  IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
                    Vinfo(24)='_FillValue'
                    Aval(6)=spval
                  END IF
#    else
                  Vinfo(20)='mask_rho'
#    endif
#   endif
                  S(ng)%pioVar(ifield)%gtype=r2dvar
              END SELECT
!
              status=def_var(ng, model, S(ng)%pioFile,                  &
     &                       S(ng)%pioVar(ifield)%vd,                   &
     &                       vtype, nvd, isfgrd, Aval, Vinfo, ncname)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
!
!
!-----------------------------------------------------------------------
!  Define the ice temperature in restart files (see the NetCDF branch
!  for the rationale): single 2D record, T = E/h.
!-----------------------------------------------------------------------
!
        IF ((FileH.eq.ABS(RST(ng)%pioFile%fh)).and.                    &
     &      (idItis.gt.0)) THEN
          DO i=1,Natt
            DO j=1,LEN(Vinfo(1))
              Vinfo(i)(j:j)=' '
            END DO
          END DO
          Aval(6)=0.0_r8
          Vinfo( 1)=Vname(1,idItis)
          Vinfo( 2)=Vname(2,idItis)
          Vinfo( 3)=Vname(3,idItis)
          Vinfo(14)=Vname(4,idItis)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idItis)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idItis,ng),r8)
          S(ng)%pioVar(idItis)%dkind=PIO_FOUT
          S(ng)%pioVar(idItis)%gtype=r2dvar
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
#   ifdef PERFECT_RESTART
          isfgrd => sr2dgrd
#   else
          isfgrd => t2dgrd
#   endif
!
          status=def_var(ng, model, S(ng)%pioFile,                     &
     &                   S(ng)%pioVar(idItis)%vd,                      &
     &                   vtype, nvd3, isfgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError,                          &
     &                   __LINE__, MyFile)) RETURN
        END IF
        IF (FileH.ne.ABS(RST(ng)%pioFile%fh)) THEN
!
!  Reset isfgrid dimension info.
!
          isfgrd => t2dgrd
          nvd=nvd3
!
!  Define isf EASTward velocity at RHO-points.
!
          IF (VarOut(idIuer,ng)) THEN
            Vinfo( 1)=Vname(1,idIuer)
            Vinfo( 2)=Vname(2,idIuer)
            Vinfo( 3)=Vname(3,idIuer)
            Vinfo(14)=Vname(4,idIuer)
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#   endif
            Vinfo(21)=Vname(6,idIuer)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idIuer,ng),r8)
            S(ng)%pioVar(idIuer)%dkind=PIO_FOUT
            S(ng)%pioVar(idIuer)%gtype=r2dvar
!
            status=def_var(ng, model, S(ng)%pioFile,                    &
     &                     S(ng)%pioVar(idIuer)%vd,                     &
     &                     vtype, nvd, isfgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!  Define isf NORTHward velocity at RHO-points.
!
          IF (VarOut(idIvnr,ng)) THEN
            Vinfo( 1)=Vname(1,idIvnr)
            Vinfo( 2)=Vname(2,idIvnr)
            Vinfo( 3)=Vname(3,idIvnr)
            Vinfo(14)=Vname(4,idIvnr)
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#   endif
            Vinfo(21)=Vname(6,idIvnr)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idIvnr,ng),r8)
            S(ng)%pioVar(idIvnr)%dkind=vtype
            S(ng)%pioVar(idIvnr)%gtype=r2dvar
!
            status=def_var(ng, model, S(ng)%pioFile,                    &
     &                     S(ng)%pioVar(idIvnr)%vd,                     &
     &                     vtype, nvd, isfgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
#  endif
#  if defined ICESHEET || defined ICESHELF_THERMO
!
!-----------------------------------------------------------------------
!  Define isf model internal variables.
!-----------------------------------------------------------------------
!
        IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
          isfgrd => sr2dgrd
        ELSE
          isfgrd => t2dgrd
        END IF
        nvd=nvd3
!
        DO nf=1,nIsfF
          IF (iFisf(nf).gt.0) THEN
            ifield=iFisf(nf)
!
!  The ice temperature in restart files is defined by the dedicated
!  block above; skip it here to avoid a duplicate variable.
!
#  ifdef ICESHEET
            IF ((nf.eq.icItis).and.                                     &
     &          (FileH.eq.ABS(RST(ng)%pioFile%fh))) THEN
              CYCLE
            END IF
#  endif
            IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
#   ifdef PERFECT_RESTART
              LdefVar=.TRUE.
#   else
              LdefVar=.FALSE.
#   endif
            ELSE
              LdefVar=VarOut(ifield,ng)
            END IF
            IF (LdefVar) THEN
              Vinfo( 1)=Vname(1,ifield)
              Vinfo( 2)=Vname(2,ifield)
              Vinfo( 3)=Vname(3,ifield)
              Vinfo(14)=Vname(4,ifield)
              Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
#    if !defined WET_DRY && defined PERFECT_RESTART
              IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
                Vinfo(24)='_FillValue'
                Aval(6)=spval
              END IF
#    else
              Vinfo(20)='mask_rho'
#    endif
#   endif
              Vinfo(21)=Vname(6,ifield)
              Vinfo(22)='coordinates'
              Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
              S(ng)%pioVar(ifield)%dkind=vtype
              S(ng)%pioVar(ifield)%gtype=r2dvar
!
              status=def_var(ng, model, S(ng)%pioFile,                  &
     &                       S(ng)%pioVar(ifield)%vd,                   &
     &                       vtype, nvd, isfgrd, Aval, Vinfo, ncname)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Define isf tracer amounts.
!-----------------------------------------------------------------------
!
!  Reset isfgrid dimension info.
!
        isfgrd => t2dgrd
        IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
#   ifdef PERFECT_RESTART
          nvd=nvd4
#   else
          nvd=nvd3
#   endif
        ELSE
          nvd=nvd3
        END IF
        IF (Lisftrc(ng)) THEN
          ifield=idIvol
          IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
            LdefVar=.TRUE.
          ELSE
            LdefVar=VarOut(ifield,ng)
          END IF
          IF (LdefVar) THEN
            Vinfo( 1)=Vname(1,ifield)
            Vinfo( 2)=Vname(2,ifield)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
#    if !defined WET_DRY && defined PERFECT_RESTART
            IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
              Vinfo(24)='_FillValue'
              Aval(6)=spval
            END IF
#    else
            Vinfo(20)='mask_rho'
#    endif
#   endif
            Vinfo(21)=Vname(6,ifield)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            S(ng)%pioVar(ifield)%dkind=vtype
            S(ng)%pioVar(ifield)%gtype=r2dvar
!
            status=def_var(ng, model, S(ng)%pioFile,                    &
     &                     S(ng)%pioVar(ifield)%vd,                     &
     &                     vtype, nvd, isfgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError,                          &
     &                     __LINE__, MyFile)) RETURN
          END IF
          DO itrc=NAT+1,NT(ng)
            IF (LtracerIsf(itrc, ng)) THEN
              ifield=idItrc(itrc)
              IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
                LdefVar=.TRUE.
              ELSE
                LdefVar=VarOut(ifield,ng)
              END IF
              IF (LdefVar) THEN
                Vinfo( 1)=Vname(1,ifield)
                Vinfo( 2)=Vname(2,ifield)
                Vinfo( 3)=Vname(3,ifield)
                Vinfo(14)=Vname(4,ifield)
                Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
#    if !defined WET_DRY && defined PERFECT_RESTART
                IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
                  Vinfo(24)='_FillValue'
                  Aval(6)=spval
                END IF
#    else
                Vinfo(20)='mask_rho'
#    endif
#   endif
                Vinfo(21)=Vname(6,ifield)
                Vinfo(22)='coordinates'
                Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
                S(ng)%pioVar(ifield)%dkind=vtype
                S(ng)%pioVar(ifield)%gtype=r2dvar
!
                status=def_var(ng, model, S(ng)%pioFile,                &
     &                         S(ng)%pioVar(ifield)%vd,                 &
     &                         vtype, nvd, isfgrd, Aval, Vinfo, ncname)
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
              END IF
            END IF
          END DO
!
!  Define ambient passive tracer concentration at RHO-points.
!
          DO itrc=NAT+1,NT(ng)
            IF (LtracerIsf(itrc, ng)) THEN
              ifield=idIamb(itrc)
              IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
                LdefVar=.FALSE.
              ELSE
                LdefVar=VarOut(ifield,ng)
              END IF
              IF (LdefVar) THEN
                Vinfo( 1)=Vname(1,ifield)
                Vinfo( 2)=Vname(2,ifield)
                Vinfo( 3)=Vname(3,ifield)
                Vinfo(14)=Vname(4,ifield)
                Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
                Vinfo(20)='mask_rho'
#   endif
                Vinfo(21)=Vname(6,ifield)
                Vinfo(22)='coordinates'
                Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
                S(ng)%pioVar(ifield)%dkind=vtype
                S(ng)%pioVar(ifield)%gtype=r2dvar
!
                status=def_var(ng, model, S(ng)%pioFile,                &
     &                         S(ng)%pioVar(ifield)%vd,                 &
     &                         vtype, nvd, isfgrd, Aval, Vinfo, ncname)
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
              END IF
            END IF
          END DO
        END IF
#  endif
      END IF DEFINE
!
!-----------------------------------------------------------------------
!  Otherwise, check existing output file and prepare for appending
!  data.
!-----------------------------------------------------------------------
!
      QUERY : IF (.not.ldef) THEN
!
!  Initialize locallogical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  ice shelf variables. Get variable IDs.
!
        IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
          vtype=PIO_FRST
        ELSE
          vtype=PIO_FOUT
        END IF
!
        DO i=1,n_var
#  ifdef ICESHEET
          DO nf=1,nIsfS
            IF (iSisf(nf).gt.0) THEN
              ifield=iSisf(nf)
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                got_var(ifield)=.TRUE.
                S(ng)%pioVar(ifield)%vd=var_desc(i)
                S(ng)%pioVar(ifield)%dkind=vtype
!
                SELECT CASE (nf)
                  CASE (isIuis)
                    S(ng)%pioVar(ifield)%gtype=p2dvar
                  CASE (isIvis)
                    S(ng)%pioVar(ifield)%gtype=p2dvar
                  CASE DEFAULT
                    S(ng)%pioVar(ifield)%gtype=r2dvar
                END SELECT
!
                CYCLE
              END IF
            END IF
          END DO
!
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idIuer))) THEN
            got_var(idIuer)=.TRUE.
            S(ng)%pioVar(idIuer)%vd=var_desc(i)
            S(ng)%pioVar(idIuer)%dkind=PIO_FOUT
            S(ng)%pioVar(idIuer)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idIvnr))) THEN
            got_var(idIvnr)=.TRUE.
            S(ng)%pioVar(idIvnr)%vd=var_desc(i)
            S(ng)%pioVar(idIvnr)%dkind=PIO_FOUT
            S(ng)%pioVar(idIvnr)%gtype=r2dvar
          END IF
#  endif
#  if defined ICESHEET || defined ICESHELF_THERMO
!
          DO nf=1,nIsfF
            IF (iFisf(nf).gt.0) THEN
              ifield=iFisf(nf)
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                got_var(ifield)=.TRUE.
                S(ng)%pioVar(ifield)%vd=var_desc(i)
                S(ng)%pioVar(ifield)%dkind=vtype
                S(ng)%pioVar(ifield)%gtype=r2dvar
                CYCLE
              END IF
            END IF
          END DO
!
          IF (Lisftrc(ng)) THEN
            ifield=idIvol
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
              got_var(ifield)=.TRUE.
              S(ng)%pioVar(ifield)%vd=var_desc(i)
              S(ng)%pioVar(ifield)%dkind=vtype
              S(ng)%pioVar(ifield)%gtype=r2dvar
              CYCLE
            END IF
            DO itrc=NAT+1,NT(ng)
              IF (LtracerIsf(itrc, ng)) THEN
                ifield=idItrc(itrc)
                IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                  got_var(ifield)=.TRUE.
                  S(ng)%pioVar(ifield)%vd=var_desc(i)
                  S(ng)%pioVar(ifield)%dkind=vtype
                  S(ng)%pioVar(ifield)%gtype=r2dvar
                  CYCLE
                END IF
              END IF
            END DO
          END IF
#  endif
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
#  ifdef ICESHEET
        DO nf=1,nIsfS
          IF (iSisf(nf).gt.0) THEN
            ifield=iSisf(nf)
            IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              RETURN
            END IF
          END IF
        END DO
!
        IF (.not.got_var(idIuer).and.VarOut(idIuer,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idIuer)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idIvnr).and.VarOut(idIvnr,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idIvnr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined ICESHEET || defined ICESHELF_THERMO
!
        DO nf=1,nIsfF
          IF (iFisf(nf).gt.0) THEN
            ifield=iFisf(nf)
!
!  The ice temperature in restart files is defined by the dedicated
!  block above; skip it here to avoid a duplicate variable.
!
#  ifdef ICESHEET
            IF ((nf.eq.icItis).and.                                     &
     &          (FileH.eq.ABS(RST(ng)%pioFile%fh))) THEN
              CYCLE
            END IF
#  endif
            IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
#   ifdef PERFECT_RESTART
              LdefVar=.TRUE.
#   else
              LdefVar=.FALSE.
#   endif
            ELSE
              LdefVar=VarOut(ifield,ng)
            END IF
            IF (.not.got_var(ifield).and.LdefVar) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              RETURN
            END IF
          END IF
        END DO
!
        IF (Lisftrc(ng)) THEN
          ifield=idIvol
          IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),    &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          DO itrc=NAT+1,NT(ng)
            IF (LtracerIsf(itrc, ng)) THEN
              ifield=idItrc(itrc)
              IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                RETURN
              END IF
            END IF
          END DO
        END IF
#  endif
      END IF QUERY
!
  10  FORMAT (/,' ISF_DEF_PIO - unable to find variable: ',a,2x,        &
     &        ' in output NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE isf_def_pio
!
!***********************************************************************
     SUBROUTINE isf_wrt_pio (ng, model, tile,                           &
    &                        LBi, UBi, LBj, UBj,                        &
    &                        VarOut, S)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      logical, intent(in) :: VarOut(NV,Ngrids)
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: LwrtVar
!
      integer :: gfactor, gtype, status
      integer :: FileH, ifield, ifld, ic, itrc, no
!
      real(dp) :: scale
!
      real(r8), pointer :: isfField2d(:,:)
      real(r8), pointer :: isfField3d(:,:,:)
      real(r8), pointer :: isfMask(:,:)
!
      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
      real(r8), allocatable, target :: Tmp2d(:,:)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", isf_wrt_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out isf prognostic variables into specified output NetCDF file.
!-----------------------------------------------------------------------
!
      FileH=ABS(S(ng)%pioFile%fh)

      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
#  if defined WRITE_WATER && defined MASKING
#   ifdef PERFECT_RESTART
      IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
        gfactor=1
      ELSE
        gfactor=-1
      END IF
#   else
      gfactor=-1
#   endif
#  else
      gfactor=1
#  endif
#  ifdef ICESHEET
!
!-----------------------------------------------------------------------
!  write out ice shelf model state prognostic variables.
!-----------------------------------------------------------------------
!
      DO ifld=1,nIsfS
        IF (iSisf(ifld).gt.0) THEN
          ifield=iSisf(ifld)
          IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
            LwrtVar=.TRUE.
          ELSE
            LwrtVar=VarOut(ifield,ng)
          END IF
!
          IF (LwrtVar) THEN
            IF ((model.eq.iNLM).and.                                    &
     &          ((FileH.eq.ABS(HIS(ng)%pioFile%fh)).or.                 &
     &           (FileH.eq.ABS(QCK(ng)%pioFile%fh)))) THEN
              isfField2d => ISF(ng) % Si(LBi:UBi,LBj:UBj,lsnew(ng),ifld)
            ELSE IF ((model.eq.iNLM).and.                               &
     &               (FileH.eq.ABS(RST(ng)%pioFile%fh))) THEN
#   ifdef PERFECT_RESTART
              isfField3d => ISF(ng) % Si(LBi:UBi,LBj:UBj,1:2,ifld)
#   else
              isfField2d => ISF(ng) % Si(LBi:UBi,LBj:UBj,lsnew(ng),ifld)
#   endif
#   ifdef AVERAGES
            ELSE IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
              isfField2d => ISF_SAVG(ifld,ng) % var(LBi:UBi,LBj:UBj)
#   endif
            END IF
!
!  The ice state must be masked by the STATIC land/sea mask only
!  (see the NetCDF branch for the reason).
!
            SELECT CASE (ifld)
              CASE (isIuis)
#   ifdef PERFECT_RESTART
                IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
                  IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                    ioDesc => ioDesc_dp_rpbar(ng)
                  ELSE
                    ioDesc => ioDesc_sp_rpbar(ng)
                  END IF
                ELSE
                  IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                    ioDesc => ioDesc_dp_p2dvar(ng)
                  ELSE
                    ioDesc => ioDesc_sp_p2dvar(ng)
                  END IF
                END IF
#   else
                IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                  ioDesc => ioDesc_dp_p2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_p2dvar(ng)
                END IF
#   endif
#   ifdef MASKING
                isfMask => GRID(ng) % pmask
#   endif
              CASE (isIvis)
#   ifdef PERFECT_RESTART
                IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
                  IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                    ioDesc => ioDesc_dp_rpbar(ng)
                  ELSE
                    ioDesc => ioDesc_sp_rpbar(ng)
                  END IF
                ELSE
                  IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                    ioDesc => ioDesc_dp_p2dvar(ng)
                  ELSE
                    ioDesc => ioDesc_sp_p2dvar(ng)
                  END IF
                END IF
#   else
                IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                  ioDesc => ioDesc_dp_p2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_p2dvar(ng)
                END IF
#   endif
#   ifdef MASKING
                isfMask => GRID(ng) % pmask
#   endif
              CASE DEFAULT
#   ifdef PERFECT_RESTART
                IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
                  IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                    ioDesc => ioDesc_dp_rzeta(ng)
                  ELSE
                    ioDesc => ioDesc_sp_rzeta(ng)
                  END IF
                ELSE
                  IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                    ioDesc => ioDesc_dp_r2dvar(ng)
                  ELSE
                    ioDesc => ioDesc_sp_r2dvar(ng)
                  END IF
                END IF
#   else
                IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                  ioDesc => ioDesc_dp_r2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_r2dvar(ng)
                END IF
#   endif
#   ifdef MASKING
                isfMask => GRID(ng) % rmask
#   endif
            END SELECT
            scale=1.0_dp
!
#   ifdef PERFECT_RESTART
            IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
              status=nf_fwrite3d(ng, model, S(ng)%pioFile, ifield,      &
     &                           S(ng)%pioVar(ifield), S(ng)%Rindex,    &
     &                           ioDesc,                                &
     &                           LBi, UBi, LBj, UBj, 1, 2, scale,       &
#    ifdef MASKING
     &                           isfMask,                               &
#    endif
     &                           isfField3d)
            ELSE
              status=nf_fwrite2d(ng, model, S(ng)%pioFile, ifield,      &
     &                           S(ng)%pioVar(ifield), S(ng)%Rindex,    &
     &                           ioDesc,                                &
     &                           LBi, UBi, LBj, UBj, scale,             &
#    ifdef MASKING
     &                           isfMask,                               &
#    endif
     &                           isfField2d)
            END IF
#   else
            status=nf_fwrite2d(ng, model, S(ng)%pioFile, ifield,        &
     &                         S(ng)%pioVar(ifield), S(ng)%Rindex,      &
     &                         ioDesc,                                  &
     &                         LBi, UBi, LBj, UBj, scale,               &
#    ifdef MASKING
     &                         isfMask,                                 &
#    endif
     &                         isfField2d)
#   endif
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      S(ng)%Rindex,               &
     &                                      TRIM(S(ng)%name)
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          END IF
        END IF
      END DO
!
!
!-----------------------------------------------------------------------
!  Write out the ice temperature (T = E/h) into restart files so a
!  restart restores both the exact enthalpy state (above) and its
!  diagnostic temperature (see the NetCDF branch).
!-----------------------------------------------------------------------
!
        IF ((FileH.eq.ABS(RST(ng)%pioFile%fh)).and.                    &
     &      (idItis.gt.0)) THEN
          IF (.not.allocated(Tmp2d)) THEN
            allocate (Tmp2d(LBi:UBi,LBj:UBj))
          END IF
          WHERE (ISF(ng)%Si(LBi:UBi,LBj:UBj,lsnew(ng),isIhis)         &
     &           .gt.0.0_r8)
            Tmp2d=ISF(ng)%Si(LBi:UBi,LBj:UBj,lsnew(ng),isIeis)/       &
     &            ISF(ng)%Si(LBi:UBi,LBj:UBj,lsnew(ng),isIhis)
          ELSEWHERE
            Tmp2d=0.0_r8
          END WHERE
          IF (S(ng)%pioVar(idItis)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dvar(ng)
          END IF
          scale=1.0_dp
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idItis,        &
     &                       S(ng)%pioVar(idItis), S(ng)%Rindex,      &
     &                       ioDesc,                                  &
     &                       LBi, UBi, LBj, UBj, scale,               &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                        &
#   endif
     &                       Tmp2d)
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,idItis)),      &
     &                                    S(ng)%Rindex,               &
     &                                    TRIM(S(ng)%name)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!  Write isf EASTward and NORTHward velocity component (m/s) at
!  RHO-points.
!
      IF (VarOut(idIuer,ng).and.VarOut(idIvnr,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
          Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
          Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
#   ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
#   endif
     &                    ISF(ng) % Si(:,:,lsnew(ng),isIuis),           &
     &                    ISF(ng) % Si(:,:,lsnew(ng),isIvis),           &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        IF (S(ng)%pioVar(idIuer)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idIuer,            &
     &                     S(ng)%pioVAR(idIuer), S(ng)%Rindex,          &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     Ur2d)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idIuer)),          &
     &                                  S(ng)%Rindex, TRIM(S(ng)%name)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        IF (S(ng)%pioVar(idIvnr)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idIvnr,            &
     &                     S(ng)%pioVar(idIvnr), S(ng)%Rindex,          &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     Vr2d)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idIvnr)),          &
     &                                  S(ng)%Rindex, TRIM(S(ng)%name)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
#  endif
#  if defined ICESHEET || defined ICESHELF_THERMO
!
!-----------------------------------------------------------------------
!  write out isf model internal variables.
!-----------------------------------------------------------------------
!
      DO ifld=1,nIsfF
        IF (iFisf(ifld).gt.0) THEN
          ifield=iFisf(ifld)
!
!  The ice temperature in restart files is written by the dedicated
!  block above; skip it here to avoid a duplicate write.
!
#  ifdef ICESHEET
          IF ((ifld.eq.icItis).and.                                     &
     &        (FileH.eq.ABS(RST(ng)%pioFile%fh))) THEN
            CYCLE
          END IF
#  endif
          IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
#   ifdef PERFECT_RESTART
            LwrtVar=.TRUE.
#   else
            LwrtVar=.FALSE.
#   endif
          ELSE
            LwrtVar=VarOut(ifield,ng)
          END IF
!
          IF (LwrtVar) THEN
            IF ((model.eq.iNLM).and.                                    &
     &          ((FileH.eq.ABS(HIS(ng)%pioFile%fh)).or.                 &
     &           (FileH.eq.ABS(QCK(ng)%pioFile%fh)).or.                 &
     &           (FileH.eq.ABS(RST(ng)%pioFile%fh)))) THEN
              isfField2d => ISF(ng) % Fi(LBi:UBi,LBj:UBj,ifld)
#   ifdef AVERAGES
            ELSE IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
              isfField2d => ISF_FAVG(ifld,ng) % var(LBi:UBi,LBj:UBj)
#   endif
            END IF
!
            scale=1.0_dp
            IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
            status=nf_fwrite2d(ng, model, S(ng)%pioFile, ifield,        &
     &                         S(ng)%pioVar(ifield), S(ng)%Rindex,      &
     &                         ioDesc,                                  &
     &                         LBi, UBi, LBj, UBj, scale,               &
#   ifdef MASKING
     &                         GRID(ng) % rmask_full,                   &
#   endif
     &                         isfField2d)
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      S(ng)%Rindex,               &
     &                                      TRIM(S(ng)%name)
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          END IF
        END IF
      END DO
#  if defined ICESHELF_THERMO && !defined ICESHELF_NO_SBL
#   ifdef PERFECT_RESTART
!
!  Write out SBL weights (Oi) to restart file.
!
      IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
        DO no=1,nIsfO
          IF (iOisf(no).gt.0) THEN
            ifield=iOisf(no)
            isfField3d => ISF(ng) % Oi(LBi:UBi,LBj:UBj,:,no)
            scale=1.0_dp
            IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
              ioDesc => ioDesc_dp_r3dvar(ng)
            ELSE
              ioDesc => ioDesc_sp_r3dvar(ng)
            END IF
            status=nf_fwrite3d(ng, model, S(ng)%pioFile, ifield,        &
     &                         S(ng)%pioVar(ifield), S(ng)%Rindex,      &
     &                         ioDesc,                                  &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#   ifdef MASKING
     &                         GRID(ng) % rmask_full,                   &
#   endif
     &                         isfField3d)
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      S(ng)%Rindex,               &
     &                                      TRIM(S(ng)%name)
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          END IF
        END DO
      END IF
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  write out isf tracer amount at RHO-points.
!-----------------------------------------------------------------------
!
      IF (Lisftrc(ng)) THEN
        ifield=idIvol
        IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
          LwrtVar=.TRUE.
        ELSE
          LwrtVar=VarOut(ifield,ng)
        END IF
!
        IF (LwrtVar) THEN
          IF ((model.eq.iNLM).and.                                      &
     &        ((FileH.eq.ABS(HIS(ng)%pioFile%fh)).or.                   &
     &         (FileH.eq.ABS(QCK(ng)%pioFile%fh)))) THEN
            isfField2d => ISF(ng) % Vi(LBi:UBi,LBj:UBj,lsnew(ng))
          ELSE IF ((model.eq.iNLM).and.                                 &
     &             (FileH.eq.ABS(RST(ng)%pioFile%fh))) THEN
#   ifdef PERFECT_RESTART
            isfField3d => ISF(ng) % Vi(LBi:UBi,LBj:UBj,1:2)
#   else
            isfField2d => ISF(ng) % Vi(LBi:UBi,LBj:UBj,lsnew(ng))
#   endif
#   ifdef AVERAGES
          ELSE IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
            isfField2d => ISF_VAVG(ng) % var(LBi:UBi,LBj:UBj)
#   endif
          END IF
!
          scale=1.0_dp
#   ifdef PERFECT_RESTART
          IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
            IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
              ioDesc => ioDesc_dp_rzeta(ng)
            ELSE
              ioDesc => ioDesc_sp_rzeta(ng)
            END IF
          ELSE
            IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
          END IF
#   else
          IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dvar(ng)
          END IF
#   endif
#   ifdef PERFECT_RESTART
          IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
            status=nf_fwrite3d(ng, model, S(ng)%pioFile, ifield,        &
     &                         S(ng)%pioVar(ifield), S(ng)%Rindex,      &
     &                         ioDesc,                                  &
     &                         LBi, UBi, LBj, UBj, 1, 2, scale,         &
#    ifdef MASKING
     &                         GRID(ng) % rmask_full,                   &
#    endif
     &                         isfField3d)
          ELSE
            status=nf_fwrite2d(ng, model, S(ng)%pioFile, ifield,        &
     &                         S(ng)%pioVar(ifield), S(ng)%Rindex,      &
     &                         ioDesc,                                  &
     &                         LBi, UBi, LBj, UBj, scale,               &
#    ifdef MASKING
     &                         GRID(ng) % rmask_full,                   &
#    endif
     &                         isfField2d)
          END IF
#   else
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, ifield,          &
     &                       S(ng)%pioVar(ifield), S(ng)%Rindex,        &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask_full,                     &
#    endif
     &                       isfField2d)
#   endif
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),        &
     &                                    S(ng)%Rindex,                 &
     &                                    TRIM(S(ng)%name)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
        ic=0
        DO itrc=NAT+1,NT(ng)
          IF (LtracerIsf(itrc, ng)) THEN
            ic=ic+1
            ifield=idItrc(itrc)
            IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
              LwrtVar=.TRUE.
            ELSE
              LwrtVar=VarOut(ifield,ng)
            END IF
!
            IF (LwrtVar) THEN
              IF ((model.eq.iNLM).and.                                  &
     &            ((FileH.eq.ABS(HIS(ng)%pioFile%fh)).or.               &
     &             (FileH.eq.ABS(QCK(ng)%pioFile%fh)))) THEN
                isfField2d => ISF(ng) % Ti(LBi:UBi,LBj:UBj,lsnew(ng),ic)
              ELSE IF ((model.eq.iNLM).and.                             &
     &                 (FileH.eq.ABS(RST(ng)%pioFile%fh))) THEN
#   ifdef PERFECT_RESTART
                isfField3d => ISF(ng) % Ti(LBi:UBi,LBj:UBj,1:2,ic)
#   else
                isfField2d => ISF(ng) % Ti(LBi:UBi,LBj:UBj,lsnew(ng),ic)
#   endif
#   ifdef AVERAGES
              ELSE IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
                isfField2d => ISF_TAVG(itrc,ng) % var(LBi:UBi,LBj:UBj)
#   endif
              END IF
!
              scale=1.0_dp
              IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
#   ifdef PERFECT_RESTART
                IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                  ioDesc => ioDesc_dp_r3dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_r3dvar(ng)
                END IF
#   else
                IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                  ioDesc => ioDesc_dp_r2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_r2dvar(ng)
                END IF
#   endif
              ELSE
                IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                  ioDesc => ioDesc_dp_r2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_r2dvar(ng)
                END IF
              END IF
#   ifdef PERFECT_RESTART
              IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
                status=nf_fwrite3d(ng, model, S(ng)%pioFile, ifield,    &
     &                             S(ng)%pioVar(ifield), S(ng)%Rindex,  &
     &                             ioDesc,                              &
     &                             LBi, UBi, LBj, UBj, 1, 2, scale,     &
#    ifdef MASKING
     &                             GRID(ng) % rmask_full,               &
#    endif
     &                             isfField3d)
              ELSE
                status=nf_fwrite2d(ng, model, S(ng)%pioFile, ifield,    &
     &                             S(ng)%pioVar(ifield), S(ng)%Rindex,  &
     &                             ioDesc,                              &
     &                             LBi, UBi, LBj, UBj, scale,           &
#    ifdef MASKING
     &                             GRID(ng) % rmask_full,               &
#    endif
     &                             isfField2d)
              END IF
#   else
              status=nf_fwrite2d(ng, model, S(ng)%pioFile, ifield,      &
     &                           S(ng)%pioVar(ifield), S(ng)%Rindex,    &
     &                           ioDesc,                                &
     &                           LBi, UBi, LBj, UBj, scale,             &
#    ifdef MASKING
     &                           GRID(ng) % rmask_full,                 &
#    endif
     &                           isfField2d)
#   endif
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),    &
     &                                        S(ng)%Rindex,             &
     &                                        TRIM(S(ng)%name)
                exit_flag=3
                ioerror=status
                RETURN
              END IF
            END IF
          END IF
        END DO
!
!  Write out ambient passive tracer concentration at RHO-points.
!
        ic=0
        DO itrc=NAT+1,NT(ng)
          IF (LtracerIsf(itrc, ng)) THEN
            ic=ic+1
            ifield=idIamb(itrc)
            IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
              LwrtVar=.FALSE.
            ELSE
              LwrtVar=VarOut(ifield,ng)
            END IF
!
            IF (LwrtVar) THEN
              IF ((model.eq.iNLM).and.                                  &
     &            ((FileH.eq.ABS(HIS(ng)%pioFile%fh)).or.               &
     &             (FileH.eq.ABS(QCK(ng)%pioFile%fh)))) THEN
                isfField2d => ISF(ng) % Ai(LBi:UBi,LBj:UBj,ic)
#   ifdef AVERAGES
              ELSE IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
                isfField2d => ISF_AAVG(itrc,ng) % var(LBi:UBi,LBj:UBj)
#   endif
              END IF
!
              scale=1.0_dp
              IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, S(ng)%pioFile, ifield,      &
     &                           S(ng)%pioVar(ifield), S(ng)%Rindex,    &
     &                           ioDesc,                                &
     &                           LBi, UBi, LBj, UBj, scale,             &
#    ifdef MASKING
     &                           GRID(ng) % rmask_full,                 &
#    endif
     &                           isfField2d)
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),    &
     &                                        S(ng)%Rindex,             &
     &                                        TRIM(S(ng)%name)
                exit_flag=3
                ioerror=status
                RETURN
              END IF
            END IF
          END IF
        END DO
      END IF
#  endif
!
  10  FORMAT (/," ISF_WRT - error while writing variable '",a,          &
     &        "', time record = ",i0,/,11x,'into file: ',a)
!
      RETURN
      END SUBROUTINE isf_wrt_pio
#  ifdef STATIONS
!
!***********************************************************************
      SUBROUTINE isf_def_station_pio (ng, model, ldef, VarOut, S,       &
     &                                pgrd)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: pgrd(:)
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer, parameter :: Natt = 25

      integer :: i, ifield, j, nf, status, ic, itrc
!
      real(r8) :: Aval(6)
!
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", isf_def_station_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define prognostic isf variables into specified output NetCDF.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!
!  Initialize local information variable arrays.
!
      DO i=1,Natt
        DO j=1,LEN(Vinfo(1))
          Vinfo(i)(j:j)=' '
        END DO
      END DO
      DO i=1,6
        Aval(i)=0.0_r8
      END DO
!
!-----------------------------------------------------------------------
!  Define isf model state prognostic variables.
!-----------------------------------------------------------------------
!
      DEFINE : IF (ldef) THEN

#   ifdef ICESHEET
        DO nf=1,nIsfS
          IF (iSisf(nf).gt.0) THEN
            ifield=iSisf(nf)
            IF (VarOut(ifield,ng)) THEN
              Vinfo( 1)=Vname(1,ifield)
              Vinfo( 2)=Vname(2,ifield)
              Vinfo( 3)=Vname(3,ifield)
              Vinfo(14)=Vname(4,ifield)
              Vinfo(16)=Vname(1,idtime)
              S(ng)%pioVar(ifield)%dkind=PIO_FOUT
              S(ng)%pioVar(ifield)%gtype=0
!
              status=def_var(ng, model, S(ng)%pioFile,                  &
     &                       S(ng)%pioVar(ifield)%vd,                   &
     &                       PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,    &
     &                       SetFillVal = .TRUE.,                       &
     &                       SetParAccess = .TRUE.)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
!
        IF (VarOut(idIuer,ng)) THEN
          Vinfo( 1)=Vname(1,idIuer)
          Vinfo( 2)=Vname(2,idIuer)
          Vinfo( 3)=Vname(3,idIuer)
          Vinfo(14)=Vname(4,idIuer)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idIuer)%dkind=PIO_FOUT
          S(ng)%pioVar(idIuer)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idIuer)%vd,                       &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (VarOut(idIvnr,ng)) THEN
          Vinfo( 1)=Vname(1,idIvnr)
          Vinfo( 2)=Vname(2,idIvnr)
          Vinfo( 3)=Vname(3,idIvnr)
          Vinfo(14)=Vname(4,idIvnr)
          Vinfo(16)=Vname(1,idtime)
          S(ng)%pioVar(idIvnr)%dkind=PIO_FOUT
          S(ng)%pioVar(idIvnr)%gtype=0
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idIvnr)%vd,                       &
     &                   PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,        &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#   if defined ICESHEET || defined ICESHELF_THERMO
!
!-----------------------------------------------------------------------
!  Define isf model internal variables.
!-----------------------------------------------------------------------
!
        DO nf=1,nIsfF
          IF (iFisf(nf).gt.0) THEN
            ifield=iFisf(nf)
            IF (VarOut(ifield,ng)) THEN
              Vinfo( 1)=Vname(1,ifield)
              Vinfo( 2)=Vname(2,ifield)
              Vinfo( 3)=Vname(3,ifield)
              Vinfo(14)=Vname(4,ifield)
              Vinfo(16)=Vname(1,idtime)
              S(ng)%pioVar(ifield)%dkind=PIO_FOUT
              S(ng)%pioVar(ifield)%gtype=0
!
              status=def_var(ng, model, S(ng)%pioFile,                  &
     &                       S(ng)%pioVar(ifield)%vd,                   &
     &                       PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,    &
     &                       SetFillVal = .TRUE.,                       &
     &                       SetParAccess = .TRUE.)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Define isf tracer amounts.
!-----------------------------------------------------------------------
!
        IF (Lisftrc(ng)) THEN
          ifield=idIvol
          IF (VarOut(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            Vinfo( 2)=Vname(2,ifield)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
            S(ng)%pioVar(ifield)%dkind=PIO_FOUT
            S(ng)%pioVar(ifield)%gtype=0
!
            status=def_var(ng, model, S(ng)%pioFile,                    &
     &                     S(ng)%pioVar(ifield)%vd,                     &
     &                     PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,      &
     &                     SetFillVal = .TRUE.,                         &
     &                     SetParAccess = .TRUE.)
            IF (FoundError(exit_flag, NoError,                          &
     &                     __LINE__, MyFile)) RETURN
          END IF
          DO itrc=NAT+1,NT(ng)
            IF (LtracerIsf(itrc, ng)) THEN
              ifield=idStrc(itrc)
              IF (VarOut(ifield,ng)) THEN
                Vinfo( 1)=Vname(1,ifield)
                Vinfo( 2)=Vname(2,ifield)
                Vinfo( 3)=Vname(3,ifield)
                Vinfo(14)=Vname(4,ifield)
                Vinfo(16)=Vname(1,idtime)
                S(ng)%pioVar(ifield)%dkind=PIO_FOUT
                S(ng)%pioVar(ifield)%gtype=0
!
                status=def_var(ng, model, S(ng)%pioFile,                &
     &                         S(ng)%pioVar(ifield)%vd,                 &
     &                         PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,  &
     &                         SetFillVal = .TRUE.,                     &
     &                         SetParAccess = .TRUE.)
                IF (FoundError(exit_flag, NoError,                      &
     &                     __LINE__, MyFile)) RETURN
              END IF
            END IF
          END DO
!
!  Define ambient passive tracer concentration.
!
          DO itrc=NAT+1,NT(ng)
            IF (LtracerIsf(itrc, ng)) THEN
              ifield=idIamb(itrc)
              IF (VarOut(ifield,ng)) THEN
                Vinfo( 1)=Vname(1,ifield)
                Vinfo( 2)=Vname(2,ifield)
                Vinfo( 3)=Vname(3,ifield)
                Vinfo(14)=Vname(4,ifield)
                Vinfo(16)=Vname(1,idtime)
                S(ng)%pioVar(ifield)%dkind=PIO_FOUT
                S(ng)%pioVar(ifield)%gtype=0
!
                status=def_var(ng, model, S(ng)%pioFile,                &
     &                         S(ng)%pioVar(ifield)%vd,                 &
     &                         PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,  &
     &                         SetFillVal = .TRUE.,                     &
     &                         SetParAccess = .TRUE.)
                IF (FoundError(exit_flag, NoError,                      &
     &                     __LINE__, MyFile)) RETURN
              END IF
            END IF
          END DO
        END IF
#   endif
      END IF DEFINE
!
!-----------------------------------------------------------------------
!  Otherwise, check existing output file and prepare for appending
!  data.
!-----------------------------------------------------------------------
!
      QUERY : IF (.not.ldef) THEN
!
!  Initialize locallogical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  ice shelf variables. Get variable IDs.
!
        DO i=1,n_var
#   ifdef ICESHEET
          DO nf=1,nIsfS
            IF (iSisf(nf).gt.0) THEN
              ifield=iSisf(nf)
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                got_var(ifield)=.TRUE.
                S(ng)%pioVar(ifield)%vd=var_desc(i)
                S(ng)%pioVar(ifield)%dkind=PIO_FOUT
                S(ng)%pioVar(ifield)%gtype=0
                CYCLE
              END IF
            END IF
          END DO
!
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idIuer))) THEN
            got_var(idIuer)=.TRUE.
            S(ng)%pioVar(idIuer)%vd=var_desc(i)
            S(ng)%pioVar(idIuer)%dkind=PIO_FOUT
            S(ng)%pioVar(idIuer)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idIvnr))) THEN
            got_var(idIvnr)=.TRUE.
            S(ng)%pioVar(idIvnr)%vd=var_desc(i)
            S(ng)%pioVar(idIvnr)%dkind=PIO_FOUT
            S(ng)%pioVar(idIvnr)%gtype=0
          END IF
#   endif
#   if defined ICESHEET || defined ICESHELF_THERMO
!
          DO nf=1,nIsfF
            IF (iFisf(nf).gt.0) THEN
              ifield=iFisf(nf)
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                got_var(ifield)=.TRUE.
                S(ng)%pioVar(ifield)%vd=var_desc(i)
                S(ng)%pioVar(ifield)%dkind=PIO_FOUT
                S(ng)%pioVar(ifield)%gtype=0
                CYCLE
              END IF
            END IF
          END DO
!
          IF (Lisftrc(ng)) THEN
            ifield=idIvol
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
              got_var(ifield)=.TRUE.
              S(ng)%pioVar(ifield)%vd=var_desc(i)
              S(ng)%pioVar(ifield)%dkind=PIO_FOUT
              S(ng)%pioVar(ifield)%gtype=0
              CYCLE
            END IF
            DO itrc=NAT+1,NT(ng)
              IF (LtracerIsf(itrc, ng)) THEN
                ifield=idStrc(itrc)
                IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                  got_var(ifield)=.TRUE.
                  S(ng)%pioVar(ifield)%vd=var_desc(i)
                  S(ng)%pioVar(ifield)%dkind=PIO_FOUT
                  S(ng)%pioVar(ifield)%gtype=0
                  CYCLE
                END IF
              END IF
            END DO
          END IF
#   endif
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
#   ifdef ICESHEET
        DO nf=1,nIsfS
          IF (iSisf(nf).gt.0) THEN
            ifield=iSisf(nf)
            IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              RETURN
            END IF
          END IF
        END DO
!
        IF (.not.got_var(idIuer).and.VarOut(idIuer,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idIuer)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idIvnr).and.VarOut(idIvnr,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idIvnr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   if defined ICESHEET || defined ICESHELF_THERMO
!
        DO nf=1,nIsfF
          IF (iFisf(nf).gt.0) THEN
            ifield=iFisf(nf)
            IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              RETURN
            END IF
          END IF
        END DO
!
        IF (Lisftrc(ng)) THEN
          ifield=idIvol
          IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),    &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          DO itrc=NAT+1,NT(ng)
            IF (LtracerIsf(itrc, ng)) THEN
              ifield=idStrc(itrc)
              IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                RETURN
              END IF
            END IF
          END DO
        END IF
#   endif
      END IF QUERY
!
  10  FORMAT (/,' ISF_DEF_STATION_PIO - unable to find variable: ',     &
     &        a,2x,' in output NetCDF file: ',a)
      RETURN
      END SUBROUTINE isf_def_station_pio
!
!***********************************************************************
      SUBROUTINE isf_wrt_station_pio (ng, model, tile,                  &
     &                                LBi, UBi, LBj, UBj,               &
     &                                VarOut, S)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: Cgrid
!
      integer :: i, ifield, ifld, ic, itrc
!
      real(dp) :: scale
!
      real(r8), dimension(Nstation(ng)) :: Xpos, Ypos, psta
!
      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", isf_wrt_station_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out isf prognostic variables into specified output NetCDF file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set switch to extract station data at native C-grid position (TRUE)
!  or at RHO-points (FALSE).
!
#   ifdef STATIONS_CGRID
      Cgrid=.TRUE.
#   else
      Cgrid=.FALSE.
#   endif
!
!  Set positions for generic extraction routine.
!
      DO i=1,Nstation(ng)
        Xpos(i)=SCALARS(ng)%SposX(i)
        Ypos(i)=SCALARS(ng)%SposY(i)
      END DO
#   ifdef ICESHEET
!
!-----------------------------------------------------------------------
!  write out ice shelf model state prognostic variables.
!-----------------------------------------------------------------------
!
      DO ifld=1,nIsfS
        IF (iSisf(ifld).gt.0) THEN
          ifield=iSisf(ifld)
          IF (VarOut(ifield,ng)) THEN
            scale=1.0_dp
            CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          scale, ISF(ng)%Si(:,:,lsnew(ng),ifld),  &
     &                          Nstation(ng), Xpos, Ypos, psta)
!
            CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,            &
     &                                TRIM(Vname(1,ifield)), psta,      &
     &                                (/1,S(ng)%Rindex/),               &
     &                                (/Nstation(ng),1/),               &
     &                                pioFile = S(ng)%pioFile,          &
     &                                pioVar = S(ng)%pioVar(ifield)%vd)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
      END DO
!
!  Write isf EASTward and NORTHward velocity component (m/s) at
!  RHO-points.
!
      IF (VarOut(idIuer,ng).and.VarOut(idIvnr,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
          Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
          Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
#    ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
#    endif
     &                    ISF(ng) % Si(:,:,lsnew(ng),isIuis),           &
     &                    ISF(ng) % Si(:,:,lsnew(ng),isIvis),           &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idIuer, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Ur2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idIuer)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idIuer)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL extract_sta2d (ng, model, Cgrid, idIvnr, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Vr2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                &
     &                            TRIM(Vname(1,idIvnr)), psta,          &
     &                            (/1,S(ng)%Rindex/),                   &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = S(ng)%pioFile,              &
     &                            pioVar = S(ng)%pioVar(idIvnr)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
#   endif
#   if defined ICESHEET || defined ICESHELF_THERMO
!
!-----------------------------------------------------------------------
!  write out isf model internal variables.
!-----------------------------------------------------------------------
!
      DO ifld=1,nIsfF
        IF (iFisf(ifld).gt.0) THEN
          ifield=iFisf(ifld)
          IF (VarOut(ifield,ng)) THEN
            scale=1.0_dp
            CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          scale, ISF(ng)%Fi(:,:,ifld),            &
     &                          Nstation(ng), Xpos, Ypos, psta)
!
            CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,            &
     &                                TRIM(Vname(1,ifield)), psta,      &
     &                                (/1,S(ng)%Rindex/),               &
     &                                (/Nstation(ng),1/),               &
     &                                pioFile = S(ng)%pioFile,          &
     &                                pioVar = S(ng)%pioVar(ifield)%vd)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  write out isf tracer amount.
!-----------------------------------------------------------------------
!
      IF (Lisftrc(ng)) THEN
        ifield=idIvol
        IF (VarOut(ifield,ng)) THEN
          scale=1.0_dp
          CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        scale, ISF(ng)%Vi(:,:,lsnew(ng)),         &
     &                        Nstation(ng), Xpos, Ypos, psta)
!
          CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,              &
     &                              TRIM(Vname(1,ifield)), psta,        &
     &                              (/1,S(ng)%Rindex/),                 &
     &                              (/Nstation(ng),1/),                 &
     &                              pioFile=S(ng)%pioFile,              &
     &                              pioVar=S(ng)%pioVar(ifield)%vd)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile))         &
     &      RETURN
        END IF
        ic=0
        DO itrc=NAT+1,NT(ng)
          IF (LtracerIsf(itrc, ng)) THEN
            ic=ic+1
            ifield=idStrc(itrc)
            IF (VarOut(ifield,ng)) THEN
              scale=1.0_dp
              CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            scale, ISF(ng)%Ti(:,:,lsnew(ng),ic),  &
     &                            Nstation(ng), Xpos, Ypos, psta)
!
              CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,          &
     &                                  TRIM(Vname(1,ifield)), psta,    &
     &                                  (/1,S(ng)%Rindex/),             &
     &                                  (/Nstation(ng),1/),             &
     &                                  pioFile=S(ng)%pioFile,          &
     &                                  pioVar=S(ng)%pioVar(ifield)%vd)
              IF (FoundError(exit_flag, NoError, __LINE__, MyFile))     &
     &          RETURN
            END IF
          END IF
        END DO
!
!  Write out ambient passive tracer concentration.
!
        ic=0
        DO itrc=NAT+1,NT(ng)
          IF (LtracerIsf(itrc, ng)) THEN
            ic=ic+1
            ifield=idIamb(itrc)
            IF (VarOut(ifield,ng)) THEN
              scale=1.0_dp
              CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            scale, ISF(ng)%Ai(:,:,ic),            &
     &                            Nstation(ng), Xpos, Ypos, psta)
!
              CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,          &
     &                                  TRIM(Vname(1,ifield)), psta,    &
     &                                  (/1,S(ng)%Rindex/),             &
     &                                  (/Nstation(ng),1/),             &
     &                                  pioFile=S(ng)%pioFile,          &
     &                                  pioVar=S(ng)%pioVar(ifield)%vd)
              IF (FoundError(exit_flag, NoError, __LINE__, MyFile))     &
     &          RETURN
            END IF
          END IF
        END DO
      END IF
#   endif
!
  10  FORMAT (/,' ISF_WRT_STATION_PIO - error while writing ',          &
     &        "variable '",a,"', time record = ",i0,/,11x,              &
     &        'into file: ',a)
!
      RETURN
      END SUBROUTINE isf_wrt_station_pio
#  endif
# endif
#endif
      END MODULE isf_output_mod
