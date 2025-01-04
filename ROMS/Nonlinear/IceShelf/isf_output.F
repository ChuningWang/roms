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
!  This module defines/writes the iceshelf model prognostic variables  !
!  into output NetCDF files.                                           !
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
     &                         t2dgrd, u2dgrd, v2dgrd)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: t2dgrd(:), u2dgrd(:), v2dgrd(:)
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: got_var(NV)
      logical :: LdefVar
!
      integer, parameter :: Natt = 25

      integer :: i, ifield, j, nf, nvd3, status, vtype, itrc

      integer :: isfgrd(3)
!

      real(r8) :: Aval(6)
!
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", isf_def_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define prognostic iceshelf variables into specified output NetCDF.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!
!-----------------------------------------------------------------------
!  Define iceshelf model state prognostic variables.
!-----------------------------------------------------------------------
!
      DEFINE : IF (ldef) THEN
!
!  Set number of dimensions for output variables.
!
# if defined WRITE_WATER && defined MASKING
        nvd3=2
# else
        nvd3=3
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

#ifdef ICESHELF_MOM
        DO nf=1,nIsfS
          IF (iSisf(nf).gt.0) THEN
            ifield=iSisf(nf)
            IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
              LdefVar=.TRUE.
              vtype=NF_FRST
            ELSE
              LdefVar=VarOut(ifield,ng)
              vtype=NF_FOUT
            END IF
            IF (Ldefvar) THEN
              Vinfo( 1)=Vname(1,ifield)
              Vinfo( 2)=Vname(2,ifield)
              vtype=NF_FOUT
              Vinfo( 3)=Vname(3,ifield)
              Vinfo(14)=Vname(4,ifield)
              Vinfo(16)=Vname(1,idtime)
              Vinfo(21)=Vname(6,ifield)
              Vinfo(22)='coordinates'
              Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
!
              SELECT CASE (nf)
                CASE (isIuis)
                  isfgrd(1:3)=u2dgrd(1:3)
# if defined WRITE_WATER && defined MASKING
                  Vinfo(20)='mask_u'
# endif
                CASE (isIvis)
                  isfgrd(1:3)=v2dgrd(1:3)
# if defined WRITE_WATER && defined MASKING
                  Vinfo(20)='mask_v'
# endif
                CASE DEFAULT
                  isfgrd(1:3)=t2dgrd(1:3)
# if defined WRITE_WATER && defined MASKING
                  Vinfo(20)='mask_rho'
# endif
              END SELECT
!
              status=def_var(ng, model, S(ng)%ncid,                     &
     &                       S(ng)%Vid(ifield),                         &
     &                       vtype, nvd3, isfgrd, Aval, Vinfo, ncname)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
!
!  Define iceshelf EASTward velocity at RHO-points.
!
        IF (VarOut(idIuer,ng)) THEN
          vtype=NF_FOUT
          Vinfo( 1)=Vname(1,idIuer)
          Vinfo( 2)=Vname(2,idIuer)
          Vinfo( 3)=Vname(3,idIuer)
          Vinfo(14)=Vname(4,idIuer)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idIuer)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idIuer,ng),r8)
          status=def_var(ng, model, S(ng)%ncid,                         &
     &                   S(ng)%Vid(idIuer),                             &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define iceshelf NORTHward velocity at RHO-points.
!
        IF (VarOut(idIvnr,ng)) THEN
          vtype=NF_FOUT
          Vinfo( 1)=Vname(1,idIvnr)
          Vinfo( 2)=Vname(2,idIvnr)
          Vinfo( 3)=Vname(3,idIvnr)
          Vinfo(14)=Vname(4,idIvnr)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idIvnr)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idIvnr,ng),r8)
          status=def_var(ng, model, S(ng)%ncid,                         &
     &                   S(ng)%Vid(idIvnr),                             &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#endif
#ifdef ICESHELF_THERMO
!
!-----------------------------------------------------------------------
!  Define iceshelf model internal variables.
!-----------------------------------------------------------------------
!
        DO nf=1,nIsfF
          IF (iFisf(nf).gt.0) THEN
            ifield=iFisf(nf)
            IF (VarOut(ifield,ng)) THEN
              vtype=NF_FOUT
              Vinfo( 1)=Vname(1,ifield)
              Vinfo( 2)=Vname(2,ifield)
              Vinfo( 3)=Vname(3,ifield)
              Vinfo(14)=Vname(4,ifield)
              Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
              Vinfo(20)='mask_rho'
# endif
              Vinfo(21)=Vname(6,ifield)
              Vinfo(22)='coordinates'
              Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
!
              status=def_var(ng, model, S(ng)%ncid,                     &
     &                       S(ng)%Vid(ifield),                         &
     &                       vtype, nvd3, t2dgrd, Aval, Vinfo, ncname)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
# ifdef ICESHELF_TRACER
!
!-----------------------------------------------------------------------
!  Define iceshelf refreezed volume.
!-----------------------------------------------------------------------
!
        ifield=idIvol
        IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
          LdefVar=.TRUE.
          vtype=NF_FRST
        ELSE
          LdefVar=VarOut(ifield,ng)
          vtype=NF_FOUT
        END IF
        IF (Ldefvar) THEN
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
          status=def_var(ng, model, S(ng)%ncid,                         &
     &                   S(ng)%Vid(ifield),                             &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define iceshelf refreezed tracer amounts.
!-----------------------------------------------------------------------
!
        DO itrc=1,NPT
          ifield=idItrc(itrc)
          IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
            LdefVar=.TRUE.
            vtype=NF_FRST
          ELSE
            LdefVar=VarOut(ifield,ng)
            vtype=NF_FOUT
          END IF
          IF (Ldefvar) THEN
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
            status=def_var(ng, model, S(ng)%ncid,                       &
     &                     S(ng)%Vid(ifield),                           &
     &                     vtype, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
# endif
#endif
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
!  iceshelf variables. Get variable IDs.
!
        DO i=1,n_var
#ifdef ICESHELF_MOM
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
#endif
#ifdef ICESHELF_THERMO
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
# ifdef ICESHELF_TRACER
!
          ifield=idIvol
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
            got_var(ifield)=.TRUE.
            S(ng)%Vid(ifield)=var_id(i)
            CYCLE
          END IF
!
          DO itrc=1,NPT
            ifield=idItrc(itrc)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
              got_var(ifield)=.TRUE.
              S(ng)%Vid(ifield)=var_id(i)
              CYCLE
            END IF
          END DO
# endif
#endif
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
#ifdef ICESHELF_MOM
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
#endif
#ifdef ICESHELF_THERMO
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
# ifdef ICESHELF_TRACER
!
        ifield=idIvol
        IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
!
        DO itrc=1,NPT
          ifield=idItrc(itrc)
          IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),        &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
#endif
      END IF QUERY
!
  10  FORMAT (/,' ICESHELF_DEF_NF90 - unable to find variable: ',a,2x,  &
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
      logical, intent(in) :: VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: LwrtVar
!
      integer :: ifield, ifld, itrc
      integer :: gfactor, gtype, status
!
      real(dp) :: scale
!
      real(r8), pointer :: isfField(:,:)
      real(r8), pointer :: isfMask(:,:)
!
      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", isf_wrt_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out iceshelf prognostic variables into specified output NetCDF
!  file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
# if defined WRITE_WATER && defined MASKING
      gfactor=-1
# else
      gfactor=1
# endif
#ifdef ICESHELF_MOM
!
!-----------------------------------------------------------------------
!  write out iceshelf model state prognostic variables.
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
     &           (S(ng)%ncid.eq.QCK(ng)%ncid).or.                       &
     &           (S(ng)%ncid.eq.RST(ng)%ncid))) THEN
              isfField => ISF(ng) % Si(LBi:UBi,LBj:UBj,lsnew(ng),ifld)
# ifdef AVERAGES
            ELSE IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              isfField => ISF_SAVG(ifld,ng) % var(LBi:UBi,LBj:UBj)
# endif
            END IF
!
            SELECT CASE (ifld)
              CASE (isIuis)
                gtype=gfactor*u2dvar
# ifdef MASKING
                isfMask => GRID(ng) % umask_full
# endif
              CASE (isIvis)
                gtype=gfactor*v2dvar
# ifdef MASKING
                isfMask => GRID(ng) % vmask_full
# endif
              CASE DEFAULT
                gtype=gfactor*r2dvar
# ifdef MASKING
                isfMask => GRID(ng) % rmask_full
# endif
            END SELECT
            scale=1.0_dp
!
            status=nf_fwrite2d(ng, model, S(ng)%ncid, ifield,           &
     &                         S(ng)%Vid(ifield),                       &
     &                         S(ng)%Rindex, gtype,                     &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         isfMask,                                 &
# endif
     &                         isfField)
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
!  Write iceshelf EASTward and NORTHward velocity component (m/s) at
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
# ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
# endif
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
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
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
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
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
#endif
#ifdef ICESHELF_THERMO
!
!-----------------------------------------------------------------------
!  write out iceshelf model internal variables.
!-----------------------------------------------------------------------
!
      DO ifld=1,nIsfF
        IF (iFisf(ifld).gt.0) THEN
          ifield=iFisf(ifld)
          IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
            LwrtVar=.FALSE.
          ELSE
            LwrtVar=VarOut(ifield,ng)
          END IF
!
          IF (LwrtVar) THEN
            IF ((model.eq.iNLM).and.                                    &
                ((S(ng)%ncid.eq.HIS(ng)%ncid).or.                       &
     &           (S(ng)%ncid.eq.QCK(ng)%ncid).or.                       &
     &           (S(ng)%ncid.eq.RST(ng)%ncid))) THEN
              isfField => ISF(ng) % Fi(LBi:UBi,LBj:UBj,ifld)
# ifdef AVERAGES
            ELSE IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              isfField => ISF_FAVG(ifld,ng) % var(LBi:UBi,LBj:UBj)
# endif
            END IF
!
            scale=1.0_dp
            gtype=gfactor*r2dvar
            status=nf_fwrite2d(ng, model, S(ng)%ncid, ifield,           &
     &                         S(ng)%Vid(ifield),                       &
     &                         S(ng)%Rindex, gtype,                     &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         GRID(ng) % rmask_full,                   &
# endif
     &                         isfField)
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
# ifdef ICESHELF_TRACER
!
!-----------------------------------------------------------------------
!  write out iceshelf refreezed volume at RHO-points.
!-----------------------------------------------------------------------
!
      ifield=idIvol
      IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
        LwrtVar=.TRUE.
      ELSE
        LwrtVar=VarOut(ifield,ng)
      END IF
!
      IF (LwrtVar) THEN
        IF ((model.eq.iNLM).and.                                        &
            ((S(ng)%ncid.eq.HIS(ng)%ncid).or.                           &
     &       (S(ng)%ncid.eq.QCK(ng)%ncid).or.                           &
     &       (S(ng)%ncid.eq.RST(ng)%ncid))) THEN
          isfField => ISF(ng) % vol(LBi:UBi,LBj:UBj)
#  ifdef AVERAGES
        ELSE IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
          isfField => ISF_VAVG(ng) % var(LBi:UBi,LBj:UBj)
#  endif
        END IF
!
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, S(ng)%ncid, ifield,               &
     &                     S(ng)%Vid(ifield),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     isfField)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),          &
     &                                  S(ng)%Rindex,                   &
     &                                  TRIM(S(ng)%name)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  write out iceshelf refreezed tracer amount at RHO-points.
!-----------------------------------------------------------------------
!
      DO itrc=1,NPT
        ifield=idItrc(itrc)
        IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
          LwrtVar=.TRUE.
        ELSE
          LwrtVar=VarOut(ifield,ng)
        END IF
!
        IF (LwrtVar) THEN
          IF ((model.eq.iNLM).and.                                      &
              ((S(ng)%ncid.eq.HIS(ng)%ncid).or.                         &
     &         (S(ng)%ncid.eq.QCK(ng)%ncid).or.                         &
     &         (S(ng)%ncid.eq.RST(ng)%ncid))) THEN
            isfField => ISF(ng) % trc(LBi:UBi,LBj:UBj,itrc)
#  ifdef AVERAGES
          ELSE IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            isfField => ISF_TAVG(itrc,ng) % var(LBi:UBi,LBj:UBj)
#  endif
          END IF
!
          scale=1.0_dp
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, model, S(ng)%ncid, ifield,             &
     &                       S(ng)%Vid(ifield),                         &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask_full,                     &
#  endif
     &                       isfField)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),        &
     &                                    S(ng)%Rindex,                 &
     &                                    TRIM(S(ng)%name)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
#endif
!
  10  FORMAT (/," ICESHELF_WRT_NF90 - error while writing variable '",a,&
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

      integer :: i, ifield, j, nf, status, itrc
!
      real(r8) :: Aval(6)
!
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", isf_def_stations_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define prognostic iceshelf variables into specified output NetCDF.
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
!  Define iceshelf model state prognostic variables.
!-----------------------------------------------------------------------
!
      DEFINE : IF (ldef) THEN

#ifdef ICESHELF_MOM
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
#endif
#ifdef ICESHELF_THERMO
!
!-----------------------------------------------------------------------
!  Define iceshelf model internal variables.
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
# ifdef ICESHELF_TRACER
!
!-----------------------------------------------------------------------
!  Define iceshelf refreezed volume.
!-----------------------------------------------------------------------
!
        ifield=idIvol
        IF (VarOut(ifield,ng)) THEN
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
!
          status=def_var(ng, model, S(ng)%ncid,                         &
     &                   S(ng)%Vid(ifield),                             &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError,                            &
     &                   __LINE__, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define iceshelf refreezed tracer amounts.
!-----------------------------------------------------------------------
!
        DO itrc=1,NPT
          ifield=idItrc(itrc)
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
        END DO
# endif
#endif
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
!  iceshelf variables. Get variable IDs.
!
        DO i=1,n_var
#ifdef ICESHELF_MOM
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
#endif
#ifdef ICESHELF_THERMO
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
# ifdef ICESHELF_TRACER
!
          ifield=idIvol
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
            got_var(ifield)=.TRUE.
            S(ng)%Vid(ifield)=var_id(i)
            CYCLE
          END IF
!
          DO itrc=1,NPT
            ifield=idItrc(itrc)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
              got_var(ifield)=.TRUE.
              S(ng)%Vid(ifield)=var_id(i)
              CYCLE
            END IF
          END DO
# endif
#endif
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
#ifdef ICESHELF_MOM
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
#endif
#ifdef ICESHELF_THERMO
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
# ifdef ICESHELF_TRACER
!
        ifield=idIvol
        IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
!
        DO itrc=1,NPT
          ifield=idItrc(itrc)
          IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),        &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
#endif
      END IF QUERY
!
  10  FORMAT (/,                                                        &
     &        ' ICESHELF_DEF_STATION_NF90 - unable to find variable: ', &
     &        a,2x,' in output NetCDF file: ',a)
      RETURN
      END SUBROUTINE isf_def_station_nf90
!
!***********************************************************************
     SUBROUTINE isf_wrt_station_nf90 (ng, model, tile,                  &
    &                                 LBi, UBi, LBj, UBj,               &
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
      integer :: i, ifield, ifld, itrc
!
      real(dp) :: scale
!
      real(r8), dimension(Nstation(ng)) :: Xpos, Ypos, psta
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", isf_wrt_station_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out iceshelf prognostic variables into specified output NetCDF
!  file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set switch to extract station data at native C-grid position (TRUE)
!  or at RHO-points (FALSE).
!
# ifdef STATIONS_CGRID
      Cgrid=.TRUE.
# else
      Cgrid=.FALSE.
# endif
!
!  Set positions for generic extraction routine.
!
      DO i=1,Nstation(ng)
        Xpos(i)=SCALARS(ng)%SposX(i)
        Ypos(i)=SCALARS(ng)%SposY(i)
      END DO
#ifdef ICESHELF_MOM
!
!-----------------------------------------------------------------------
!  write out iceshelf model state prognostic variables.
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
#endif
#ifdef ICESHELF_THERMO
!
!-----------------------------------------------------------------------
!  write out iceshelf model internal variables.
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
# ifdef ICESHELF_TRACER
!
!-----------------------------------------------------------------------
!  write out iceshelf refreezed volume.
!-----------------------------------------------------------------------
!
      ifield=idIvol
      IF (VarOut(ifield,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ISF(ng)%vol(:,:),                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
!
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,ifield)), psta,              &
     &                        (/1,S(ng)%Rindex/),                       &
     &                        (/Nstation(ng),1/),                       &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(ifield))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!-----------------------------------------------------------------------
!  write out iceshelf refreezed tracer amount.
!-----------------------------------------------------------------------
!
      DO itrc=1,NPT
        ifield=idItrc(itrc)
        IF (VarOut(ifield,ng)) THEN
          scale=1.0_dp
          CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        scale, ISF(ng)%trc(:,:,itrc),             &
     &                        Nstation(ng), Xpos, Ypos, psta)
!
          CALL netcdf_put_fvar (ng, model, S(ng)%name,                  &
     &                          TRIM(Vname(1,ifield)), psta,            &
     &                          (/1,S(ng)%Rindex/),                     &
     &                          (/Nstation(ng),1/),                     &
     &                          ncid = S(ng)%ncid,                      &
     &                          varid = S(ng)%Vid(ifield))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
# endif
#endif
!
  10  FORMAT (/,' ICESHELF_WRT_STATION_NF90 - error while writing ',    &
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
     &                        t2dgrd, u2dgrd, v2dgrd)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: t2dgrd(:), u2dgrd(:), v2dgrd(:)
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: LdefVar
!
      integer, parameter :: Natt = 25

      integer :: i, ifield, j, nf, nvd3, status, vtype, itrc

      integer :: isfgrd(3)
!
      real(r8) :: Aval(6)
!
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", isf_def_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define iceshelf model state prognostic variables.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!!
      DEFINE : IF (ldef) THEN
!
!  Set number of dimensions for output variables.
!
# if defined WRITE_WATER && defined MASKING
        nvd3=2
# else
        nvd3=3
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
#ifdef ICESHELF_MOM
!
!  Define iceshelf model state prognostic variables.
!
        DO nf=1,nIsfS
          IF (iSisf(nf).gt.0) THEN
            ifield=iSisf(nf)
            IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
              LdefVar=.TRUE.
              vtype=PIO_FRST
            ELSE
              LdefVar=VarOut(ifield,ng)
              vtype=PIO_FOUT
            END IF

            IF (Ldefvar) THEN
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
                  isfgrd(1:3)=u2dgrd(1:3)
# if defined WRITE_WATER && defined MASKING
                  Vinfo(20)='mask_u'
# endif
                  S(ng)%pioVar(ifield)%gtype=u2dvar
                CASE (isIvis)
                  isfgrd(1:3)=v2dgrd(1:3)
# if defined WRITE_WATER && defined MASKING
                  Vinfo(20)='mask_v'
# endif
                  S(ng)%pioVar(ifield)%gtype=v2dvar
                CASE DEFAULT
                  isfgrd(1:3)=t2dgrd(1:3)
# if defined WRITE_WATER && defined MASKING
                  Vinfo(20)='mask_rho'
# endif
                  S(ng)%pioVar(ifield)%gtype=r2dvar
              END SELECT
!
              status=def_var(ng, model, S(ng)%pioFile,                  &
     &                       S(ng)%pioVar(ifield)%vd,                   &
     &                       vtype, nvd3, isfgrd,                       &
     &                       Aval, Vinfo, ncname)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
!
!  Define iceshelf EASTward velocity at RHO-points.
!
        IF (VarOut(idIuer,ng)) THEN
          vtype=PIO_FOUT
          Vinfo( 1)=Vname(1,idIuer)
          Vinfo( 2)=Vname(2,idIuer)
          Vinfo( 3)=Vname(3,idIuer)
          Vinfo(14)=Vname(4,idIuer)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idIuer)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idIuer,ng),r8)
          S(ng)%pioVar(idIuer)%dkind=PIO_FOUT
          S(ng)%pioVar(idIuer)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idIuer)%vd,                       &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define iceshelf NORTHward velocity at RHO-points.
!
        IF (VarOut(idIvnr,ng)) THEN
          vtype=PIO_FOUT
          Vinfo( 1)=Vname(1,idIvnr)
          Vinfo( 2)=Vname(2,idIvnr)
          Vinfo( 3)=Vname(3,idIvnr)
          Vinfo(14)=Vname(4,idIvnr)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idIvnr)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idIvnr,ng),r8)
          S(ng)%pioVar(idIvnr)%dkind=PIO_FOUT
          S(ng)%pioVar(idIvnr)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idIvnr)%vd,                       &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#endif
#ifdef ICESHELF_THERMO
!
!-----------------------------------------------------------------------
!  Define iceshelf model internal variables.
!-----------------------------------------------------------------------
!
        DO nf=1,nIsfF
          IF (iFisf(nf).gt.0) THEN
            ifield=iFisf(nf)
            IF (VarOut(ifield,ng)) THEN
              vtype=PIO_FOUT
              Vinfo( 1)=Vname(1,ifield)
              Vinfo( 2)=Vname(2,ifield)
              Vinfo( 3)=Vname(3,ifield)
              Vinfo(14)=Vname(4,ifield)
              Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
              Vinfo(20)='mask_rho'
# endif
              Vinfo(21)=Vname(6,ifield)
              Vinfo(22)='coordinates'
              Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
              S(ng)%pioVar(ifield)%dkind=PIO_FOUT
              S(ng)%pioVar(ifield)%gtype=r2dvar
!
              status=def_var(ng, model, S(ng)%pioFile,                  &
     &                       S(ng)%pioVar(ifield)%vd,                   &
     &                       vtype, nvd3, t2dgrd,                       &
     &                       Aval, Vinfo, ncname)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
# ifdef ICESHELF_TRACER
!
!-----------------------------------------------------------------------
!  Define iceshelf refreezed volume.
!-----------------------------------------------------------------------
!
        ifield=idIvol
        IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
          LdefVar=.TRUE.
          vtype=PIO_FRST
        ELSE
          LdefVar=VarOut(ifield,ng)
          vtype=PIO_FOUT
        END IF
        IF (Ldefvar) THEN
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
          S(ng)%pioVar(ifield)%dkind=PIO_FOUT
          S(ng)%pioVar(ifield)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(ifield)%vd,                       &
     &                   vtype, nvd3, t2dgrd,                           &
     &                   Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError,                            &
     &                   __LINE__, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define iceshelf refreezed tracer amounts.
!-----------------------------------------------------------------------
!
        DO itrc=1,NPT
          ifield=idItrc(itrc)
          IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
            LdefVar=.TRUE.
            vtype=PIO_FRST
          ELSE
            LdefVar=VarOut(ifield,ng)
            vtype=PIO_FOUT
          END IF
          IF (Ldefvar) THEN
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
            S(ng)%pioVar(ifield)%dkind=PIO_FOUT
            S(ng)%pioVar(ifield)%gtype=r2dvar
!
            status=def_var(ng, model, S(ng)%pioFile,                    &
     &                     S(ng)%pioVar(ifield)%vd,                     &
     &                     vtype, nvd3, t2dgrd,                         &
     &                     Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError,                          &
     &                     __LINE__, MyFile)) RETURN
          END IF
        END DO
# endif
#endif
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
!  iceshelf variables. Get variable IDs.
!
        IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
          vtype=PIO_FRST
        ELSE
          vtype=PIO_FOUT
        END IF
!
        DO i=1,n_var
#ifdef ICESHELF_MOM
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
                    S(ng)%pioVar(ifield)%gtype=u2dvar
                  CASE (isIvis)
                    S(ng)%pioVar(ifield)%gtype=v2dvar
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
#endif
#ifdef ICESHELF_THERMO
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
# ifdef ICESHELF_TRACER
!
          ifield=idIvol
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
            got_var(ifield)=.TRUE.
            S(ng)%pioVar(ifield)%vd=var_desc(i)
            S(ng)%pioVar(ifield)%dkind=vtype
            S(ng)%pioVar(ifield)%gtype=r2dvar
            CYCLE
          END IF
!
          DO itrc=1,NPT
            ifield=idItrc(itrc)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
              got_var(ifield)=.TRUE.
              S(ng)%pioVar(ifield)%vd=var_desc(i)
              S(ng)%pioVar(ifield)%dkind=vtype
              S(ng)%pioVar(ifield)%gtype=r2dvar
              CYCLE
            END IF
          END DO
# endif
#endif
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
#ifdef ICESHELF_MOM
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
#endif
#ifdef ICESHELF_THERMO
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
# ifdef ICESHELF_TRACER
!
        ifield=idIvol
        IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
!
        DO itrc=1,NPT
          ifield=idItrc(itrc)
          IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),        &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
#endif
      END IF QUERY
!
  10  FORMAT (/,' ICESHELF_DEF_PIO - unable to find variable: ',a,2x,   &
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
      logical, intent(in) :: VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: LwrtVar
!
      integer :: ifield, ifld, itrc
      integer :: gfactor, gtype, status
!
      real(dp) :: scale
!
      real(r8), pointer :: isfField(:,:)
      real(r8), pointer :: isfMask(:,:)
!
      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", isf_wrt_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out iceshelf prognostic variables into specified output NetCDF
!  file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
# if defined WRITE_WATER && defined MASKING
      gfactor=-1
# else
      gfactor=1
# endif
!
!-----------------------------------------------------------------------
!  write out iceshelf model state prognostic variables.
!-----------------------------------------------------------------------
!
#ifdef ICESHELF_MOM
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
     &           (S(ng)%ncid.eq.QCK(ng)%ncid).or.                       &
     &           (S(ng)%ncid.eq.RST(ng)%ncid))) THEN
              isfField => ISF(ng) % Si(LBi:UBi,LBj:UBj,lsnew(ng),ifld)
# ifdef AVERAGES
            ELSE IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              isfField => ISF_SAVG(ifld,ng) % var(LBi:UBi,LBj:UBj)
# endif
            END IF
!
            SELECT CASE (ifld)
              CASE (isIuis)
                IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                  ioDesc => ioDesc_dp_u2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_u2dvar(ng)
                END IF
# ifdef MASKING
                isfMask => GRID(ng) % umask_full
# endif
              CASE (isIvis)
                IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                  ioDesc => ioDesc_dp_v2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_v2dvar(ng)
                END IF
# ifdef MASKING
                isfMask => GRID(ng) % vmask_full
# endif
              CASE DEFAULT
                IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                  ioDesc => ioDesc_dp_r2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_r2dvar(ng)
                END IF
# ifdef MASKING
                isfMask => GRID(ng) % rmask_full
# endif
            END SELECT
            scale=1.0_dp
!
            status=nf_fwrite2d(ng, model, S(ng)%pioFile, ifield,        &
     &                         S(ng)%pioVar(ifield), S(ng)%Rindex,      &
     &                         ioDesc,                                  &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         isfMask,                                 &
# endif
     &                         isfField)
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
!  Write iceshelf EASTward and NORTHward velocity component (m/s) at
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
# ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
# endif
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
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
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
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
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
#endif
#ifdef ICESHELF_THERMO
!
!-----------------------------------------------------------------------
!  write out iceshelf model internal variables.
!-----------------------------------------------------------------------
!
      DO ifld=1,nIsfF
        IF (iFisf(ifld).gt.0) THEN
          ifield=iFisf(ifld)
          IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
            LwrtVar=.FALSE.
          ELSE
            LwrtVar=VarOut(ifield,ng)
          END IF
!
          IF (LwrtVar) THEN
            IF ((model.eq.iNLM).and.                                    &
                ((S(ng)%ncid.eq.HIS(ng)%ncid).or.                       &
     &           (S(ng)%ncid.eq.QCK(ng)%ncid).or.                       &
     &           (S(ng)%ncid.eq.RST(ng)%ncid))) THEN
              isfField => ISF(ng) % Fi(LBi:UBi,LBj:UBj,ifld)
# ifdef AVERAGES
            ELSE IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              isfField => ISF_FAVG(ifld,ng) % var(LBi:UBi,LBj:UBj)
# endif
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
# ifdef MASKING
     &                         GRID(ng) % rmask_full,                   &
# endif
     &                         isfField)
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
# ifdef ICESHELF_TRACER
!
!-----------------------------------------------------------------------
!  write out iceshelf refreezed volume at RHO-points.
!-----------------------------------------------------------------------
!
      ifield=idIvol
      IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
        LwrtVar=.TRUE.
      ELSE
        LwrtVar=VarOut(ifield,ng)
      END IF
!
      IF (LwrtVar) THEN
        IF ((model.eq.iNLM).and.                                        &
            ((S(ng)%ncid.eq.HIS(ng)%ncid).or.                           &
     &       (S(ng)%ncid.eq.QCK(ng)%ncid).or.                           &
     &       (S(ng)%ncid.eq.RST(ng)%ncid))) THEN
          isfField => ISF(ng) % vol(LBi:UBi,LBj:UBj)
#  ifdef AVERAGES
        ELSE IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
          isfField => ISF_VAVG(ng) % var(LBi:UBi,LBj:UBj)
#  endif
        END IF
!
        scale=1.0_dp
        IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, ifield,            &
     &                     S(ng)%pioVar(ifield), S(ng)%Rindex,          &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     isfField)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),          &
     &                                  S(ng)%Rindex,                   &
     &                                  TRIM(S(ng)%name)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  write out iceshelf refreezed tracer amount at RHO-points.
!-----------------------------------------------------------------------
!
      DO itrc=1,NPT
        ifield=idItrc(itrc)
        IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
          LwrtVar=.TRUE.
        ELSE
          LwrtVar=VarOut(ifield,ng)
        END IF
!
        IF (LwrtVar) THEN
          IF ((model.eq.iNLM).and.                                      &
              ((S(ng)%ncid.eq.HIS(ng)%ncid).or.                         &
     &         (S(ng)%ncid.eq.QCK(ng)%ncid).or.                         &
     &         (S(ng)%ncid.eq.RST(ng)%ncid))) THEN
            isfField => ISF(ng) % trc(LBi:UBi,LBj:UBj,itrc)
#  ifdef AVERAGES
          ELSE IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
            isfField => ISF_TAVG(itrc,ng) % var(LBi:UBi,LBj:UBj)
#  endif
          END IF
!
          scale=1.0_dp
          IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dvar(ng)
          END IF
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, ifield,          &
     &                       S(ng)%pioVar(ifield), S(ng)%Rindex,        &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask_full,                     &
#  endif
     &                       isfField)
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),        &
     &                                    S(ng)%Rindex,                 &
     &                                    TRIM(S(ng)%name)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
#endif
!
  10  FORMAT (/," ICESHELF_WRT - error while writing variable '",a,     &
     &        "', time record = ",i0,/,11x,'into file: ',a)
!
      RETURN
      END SUBROUTINE isf_wrt_pio
# endif
!
#endif
      END MODULE isf_output_mod
