#include "cppdefs.h"

      MODULE ice_output_mod

#ifdef ICE_MODEL
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2025 The ROMS Group               Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module defines/writes the ice model prognostic variables into  !
!  output NetCDF files.                                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_ice
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
      PUBLIC :: ice_def_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: ice_def_pio
# endif
# ifdef STATIONS
      PUBLIC :: ice_def_station_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: ice_def_station_pio
#  endif
# endif
      PUBLIC :: ice_wrt_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: ice_wrt_pio
# endif
# ifdef STATIONS
      PUBLIC :: ice_wrt_station_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: ice_wrt_station_pio
#  endif
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_def_nf90 (ng, model, ldef, VarOut, S,              &
     &                         t2dgrd, u2dgrd, v2dgrd,                  &
     &                         sr2dgrd, ru2dgrd, rv2dgrd)
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
     &  t2dgrd(:), u2dgrd(:), v2dgrd(:),                                &
     &  sr2dgrd(:), ru2dgrd(:), rv2dgrd(:)
!
!  Local variable declarations.
!
      logical :: LdefVar, got_var(NV)
!
      integer, parameter :: Natt = 25

      integer :: i, j, ifield, nf, nvd3, nvd4, nvd, status, vtype
      integer :: ic, itrc
!
      integer, pointer :: icegrd(:)
!
      real(r8) :: Aval(6)
!
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ice_def_nf90"
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
!
!-----------------------------------------------------------------------
!  Define ice model state prognostic variables.
!-----------------------------------------------------------------------
!
        DO nf=1,nIceS
          IF (iSice(nf).gt.0) THEN
            ifield=iSice(nf)
            IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
              LdefVar=.TRUE.
# ifdef PERFECT_RESTART
              nvd=nvd4
# else
              nvd=nvd3
# endif
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
                CASE (isUice)
# ifdef PERFECT_RESTART
                  IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
                    icegrd => ru2dgrd
                  ELSE
                    icegrd => u2dgrd
                  END IF
# else
                  icegrd => u2dgrd
# endif
# if defined WRITE_WATER && defined MASKING
#  if !defined WET_DRY && defined PERFECT_RESTART
                  IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
                    Vinfo(24)='_FillValue'
                    Aval(6)=spval
                  END IF
#  else
                  Vinfo(20)='mask_u'
#  endif
# endif
                CASE (isVice)
# ifdef PERFECT_RESTART
                  IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
                    icegrd => rv2dgrd
                  ELSE
                    icegrd => v2dgrd
                  END IF
# else
                  icegrd => v2dgrd
# endif
# if defined WRITE_WATER && defined MASKING
#  if !defined WET_DRY && defined PERFECT_RESTART
                  IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
                    Vinfo(24)='_FillValue'
                    Aval(6)=spval
                  END IF
#  else
                  Vinfo(20)='mask_v'
#  endif
# endif
                CASE DEFAULT
                  icegrd => t2dgrd
# if defined WRITE_WATER && defined MASKING
#  if !defined WET_DRY && defined PERFECT_RESTART
                  IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
                    Vinfo(24)='_FillValue'
                    Aval(6)=spval
                  END IF
#  else
                  Vinfo(20)='mask_rho'
#  endif
# endif
              END SELECT
!
              status=def_var(ng, model, S(ng)%ncid,                     &
     &                       S(ng)%Vid(ifield),                         &
     &                       vtype, nvd, icegrd, Aval, Vinfo, ncname)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
!
        IF (S(ng)%ncid.ne.RST(ng)%ncid) THEN
!
!  Reset icegrid dimension info.
!
          icegrd => t2dgrd
          nvd=nvd3
!
!  Define ice EASTward velocity at RHO-points.
!
          IF (VarOut(idUiER,ng)) THEN
            Vinfo( 1)=Vname(1,idUiER)
            Vinfo( 2)=Vname(2,idUiER)
            Vinfo( 3)=Vname(3,idUiER)
            Vinfo(14)=Vname(4,idUiER)
            Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
# endif
            Vinfo(21)=Vname(6,idUiER)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idUiER,ng),r8)
            status=def_var(ng, model, S(ng)%ncid,                       &
     &                     S(ng)%Vid(idUiER),                           &
     &                     vtype, nvd, icegrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!  Define ice NORTHward velocity at RHO-points.
!
          IF (VarOut(idViNR,ng)) THEN
            Vinfo( 1)=Vname(1,idViNR)
            Vinfo( 2)=Vname(2,idViNR)
            Vinfo( 3)=Vname(3,idViNR)
            Vinfo(14)=Vname(4,idViNR)
            Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
# endif
            Vinfo(21)=Vname(6,idViNR)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idViNR,ng),r8)
            status=def_var(ng, model, S(ng)%ncid,                       &
     &                     S(ng)%Vid(idViNR),                           &
     &                     vtype, nvd, icegrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
!
!-----------------------------------------------------------------------
!  Define ice model internal variables.
!-----------------------------------------------------------------------
!
        IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
          icegrd => sr2dgrd
        ELSE
          icegrd => t2dgrd
        END IF
        nvd=nvd3
!
        DO nf=1,nIceF
          IF (iFice(nf).gt.0) THEN
            ifield=iFice(nf)
            IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
# ifdef PERFECT_RESTART
              LdefVar=.TRUE.
# else
              LdefVar=.FALSE.
# endif
            ELSE
              LdefVar=VarOut(ifield,ng)
            END IF
            IF (LdefVar) THEN
              Vinfo( 1)=Vname(1,ifield)
              Vinfo( 2)=Vname(2,ifield)
              Vinfo( 3)=Vname(3,ifield)
              Vinfo(14)=Vname(4,ifield)
              Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
#  if !defined WET_DRY && defined PERFECT_RESTART
              IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
                Vinfo(24)='_FillValue'
                Aval(6)=spval
              END IF
#  else
              Vinfo(20)='mask_rho'
#  endif
# endif
              Vinfo(21)=Vname(6,ifield)
              Vinfo(22)='coordinates'
              Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
!
              SELECT CASE (nf)
                CASE (icWdiv)
                  Vinfo(11)='increase ice thickness'
                  Vinfo(12)='decrease ice concentration'
                CASE (icW_ai)
                  Vinfo(11)='freezing'
                  Vinfo(12)='melting'
                CASE (icW_ao, icW_io)
                  Vinfo(11)='melting'
                  Vinfo(12)='freezing'
              END SELECT
!
              status=def_var(ng, model, S(ng)%ncid,                     &
     &                       S(ng)%Vid(ifield),                         &
     &                       vtype, nvd, icegrd, Aval, Vinfo, ncname)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Define ice tracer amounts.
!-----------------------------------------------------------------------
!
!  Reset icegrid dimension info.
!
        icegrd => t2dgrd
        IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
# ifdef PERFECT_RESTART
          nvd=nvd4
# else
          nvd=nvd3
# endif
        ELSE
          nvd=nvd3
        END IF
        IF (Licetrc(ng)) THEN
          DO itrc=NAT+1,NT(ng)
            IF (LtracerIce(itrc,ng)) THEN
              ifield=idStrc(itrc)
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
# if defined WRITE_WATER && defined MASKING
#  if !defined WET_DRY && defined PERFECT_RESTART
                IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
                  Vinfo(24)='_FillValue'
                  Aval(6)=spval
                END IF
#  else
                Vinfo(20)='mask_rho'
#  endif
# endif
                Vinfo(21)=Vname(6,ifield)
                Vinfo(22)='coordinates'
                Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
                status=def_var(ng, model, S(ng)%ncid,                   &
     &                         S(ng)%Vid(ifield),                       &
     &                         vtype, nvd, icegrd, Aval, Vinfo, ncname)
                IF (FoundError(exit_flag, NoError, __LINE__, MyFile))   &
     &            RETURN
              END IF
            END IF
          END DO
        END IF
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
!  sea-ice variables. Get variable IDs.
!
        DO i=1,n_var
          DO nf=1,nIceS
            IF (iSice(nf).gt.0) THEN
              ifield=iSice(nf)
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                got_var(ifield)=.TRUE.
                S(ng)%Vid(ifield)=var_id(i)
                CYCLE
              END IF
            END IF
          END DO
!
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUiER))) THEN
            got_var(idUiER)=.TRUE.
            S(ng)%Vid(idUiER)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idViNR))) THEN
            got_var(idViNR)=.TRUE.
            S(ng)%Vid(idViNR)=var_id(i)
          END IF
!
          DO nf=1,nIceF
            IF (iFice(nf).gt.0) THEN
              ifield=iFice(nf)
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                got_var(ifield)=.TRUE.
                S(ng)%Vid(ifield)=var_id(i)
                CYCLE
              END IF
            END IF
          END DO
!
          IF (Licetrc(ng)) THEN
            DO itrc=NAT+1,NT(ng)
              IF (LtracerIce(itrc, ng)) THEN
                ifield=idStrc(itrc)
                IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                  got_var(ifield)=.TRUE.
                  S(ng)%Vid(ifield)=var_id(i)
                  CYCLE
                END IF
              END IF
            END DO
          END IF
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
        DO nf=1,nIceS
          IF (iSice(nf).gt.0) THEN
            ifield=iSice(nf)
            IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              RETURN
            END IF
          END IF
        END DO
!
        IF (.not.got_var(idUiER).and.VarOut(idUiER,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUiER)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idViNR).and.VarOut(idViNR,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idViNR)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
!
        DO nf=1,nIceF
          IF (iFice(nf).gt.0) THEN
            ifield=iFice(nf)
            IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
# ifdef PERFECT_RESTART
              LdefVar=.TRUE.
# else
              LdefVar=.FALSE.
# endif
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
        IF (Licetrc(ng)) THEN
          DO itrc=NAT+1,NT(ng)
            IF (LtracerIce(itrc, ng)) THEN
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
      END IF QUERY
!
  10  FORMAT (/,' ICE_DEF_NF90 - unable to find variable: ',a,2x,       &
     &        ' in output NetCDF file: ',a)
      RETURN
      END SUBROUTINE ice_def_nf90
!
!***********************************************************************
     SUBROUTINE ice_wrt_nf90 (ng, model, tile,                          &
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
      integer :: ifield, ifld, ic, itrc
!
      real(dp) :: scale
!
      real(r8), pointer :: iceField2d(:,:)
      real(r8), pointer :: iceField3d(:,:,:)
      real(r8), pointer :: iceMask(:,:)
!
      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ice_wrt_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out ice prognostic variables into specified output NetCDF file.
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
!
!-----------------------------------------------------------------------
!  write out sea-ice model state prognostic variables.
!-----------------------------------------------------------------------
!
      DO ifld=1,nIceS
        IF (iSice(ifld).gt.0) THEN
          ifield=iSice(ifld)
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
              iceField2d => ICE(ng) % Si(LBi:UBi,LBj:UBj,IOUT,ifld)
            ELSE IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
# ifdef PERFECT_RESTART
              iceField3d => ICE(ng) % Si(LBi:UBi,LBj:UBj,1:2,ifld)
# else
              iceField2d => ICE(ng) % Si(LBi:UBi,LBj:UBj,IOUT,ifld)
# endif
# ifdef AVERAGES
            ELSE IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              iceField2d => ICE_SAVG(ifld,ng) % var(LBi:UBi,LBj:UBj)
# endif
            END IF
!
            SELECT CASE (ifld)
              CASE (isUice)
                IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
# ifdef PERFECT_RESTART
                  gtype=gfactor*u3dvar
# else
                  gtype=gfactor*u2dvar
# endif
                ELSE
                  gtype=gfactor*u2dvar
                END IF
# ifdef MASKING
                iceMask => GRID(ng) % umask_full
# endif
              CASE (isVice)
                IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
# ifdef PERFECT_RESTART
                  gtype=gfactor*v3dvar
# else
                  gtype=gfactor*v2dvar
# endif
                ELSE
                  gtype=gfactor*v2dvar
                END IF
# ifdef MASKING
                iceMask => GRID(ng) % vmask_full
# endif
              CASE DEFAULT
                IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
# ifdef PERFECT_RESTART
                  gtype=gfactor*r3dvar
# else
                  gtype=gfactor*r2dvar
# endif
                ELSE
                  gtype=gfactor*r2dvar
                END IF
# ifdef MASKING
                iceMask => GRID(ng) % rmask_full
# endif
            END SELECT
            scale=1.0_dp
!
# ifdef PERFECT_RESTART
            IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
              status=nf_fwrite3d(ng, model, S(ng)%ncid, ifield,         &
     &                           S(ng)%Vid(ifield),                     &
     &                           S(ng)%Rindex, gtype,                   &
     &                           LBi, UBi, LBj, UBj, 1, 2, scale,       &
#  ifdef MASKING
     &                           iceMask,                               &
#  endif
     &                           iceField3d)
            ELSE
              status=nf_fwrite2d(ng, model, S(ng)%ncid, ifield,         &
     &                           S(ng)%Vid(ifield),                     &
     &                           S(ng)%Rindex, gtype,                   &
     &                           LBi, UBi, LBj, UBj, scale,             &
#  ifdef MASKING
     &                           iceMask,                               &
#  endif
     &                           iceField2d)
            END IF
# else
            status=nf_fwrite2d(ng, model, S(ng)%ncid, ifield,           &
     &                         S(ng)%Vid(ifield),                       &
     &                         S(ng)%Rindex, gtype,                     &
     &                         LBi, UBi, LBj, UBj, scale,               &
#  ifdef MASKING
     &                         iceMask,                                 &
#  endif
     &                         iceField2d)
# endif
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
!  Write ice EASTward and NORTHward velocity component (m/s) at
!  RHO-points.
!
      IF (VarOut(idUiER,ng).and.VarOut(idViNR,ng)) THEN
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
     &                    ICE(ng) % Si(:,:,IOUT,isUice),                &
     &                    ICE(ng) % Si(:,:,IOUT,isVice),                &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idUiER,               &
     &                     S(ng)%Vid(idUiER),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Ur2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUiER)),          &
     &                                  S(ng)%Rindex, TRIM(S(ng)%name)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idViNR,               &
     &                     S(ng)%Vid(idViNR),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Vr2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idViNR)),          &
     &                                  S(ng)%Rindex, TRIM(S(ng)%name)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
!
!-----------------------------------------------------------------------
!  write out ice model internal variables.
!-----------------------------------------------------------------------
!
      DO ifld=1,nIceF
        IF (iFice(ifld).gt.0) THEN
          ifield=iFice(ifld)
          IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
# ifdef PERFECT_RESTART
            LwrtVar=.TRUE.
# else
            LwrtVar=.FALSE.
# endif
          ELSE
            LwrtVar=VarOut(ifield,ng)
          END IF
!
          IF (LwrtVar) THEN
            IF ((model.eq.iNLM).and.                                    &
                ((S(ng)%ncid.eq.HIS(ng)%ncid).or.                       &
     &           (S(ng)%ncid.eq.QCK(ng)%ncid).or.                       &
     &           (S(ng)%ncid.eq.RST(ng)%ncid))) THEN
              iceField2d => ICE(ng) % Fi(LBi:UBi,LBj:UBj,ifld)
# ifdef AVERAGES
            ELSE IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              iceField2d => ICE_FAVG(ifld,ng) % var(LBi:UBi,LBj:UBj)
# endif
            END IF
!
# ifdef ICE_BULK_FLUXES
            IF ((ifld.eq.icSics).or.(ifld.eq.icSicb).or.                &
     &          (ifld.eq.icSocn)) THEN
              scale=rho0*Cp
            ELSE
              scale=1.0_dp
            END IF
# else
            scale=1.0_dp
# endif
            gtype=gfactor*r2dvar
            status=nf_fwrite2d(ng, model, S(ng)%ncid, ifield,           &
     &                         S(ng)%Vid(ifield),                       &
     &                         S(ng)%Rindex, gtype,                     &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         GRID(ng) % rmask_full,                   &
# endif
     &                         iceField2d)
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
!-----------------------------------------------------------------------
!  write out ice tracer amount at RHO-points.
!-----------------------------------------------------------------------
!
      IF (Licetrc(ng)) THEN
        ic=0
        DO itrc=NAT+1,NT(ng)
          IF (LtracerIce(itrc, ng)) THEN
            ic=ic+1
            ifield=idStrc(itrc)
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
                iceField2d => ICE(ng) % Ti(LBi:UBi,LBj:UBj,IOUT,ic)
              ELSE IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
# ifdef PERFECT_RESTART
                iceField3d => ICE(ng) % Ti(LBi:UBi,LBj:UBj,1:2,ic)
# else
                iceField2d => ICE(ng) % Ti(LBi:UBi,LBj:UBj,IOUT,ic)
# endif
# ifdef AVERAGES
              ELSE IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
                iceField2d => ICE_TAVG(itrc,ng) % var(LBi:UBi,LBj:UBj)
# endif
              END IF
!
              scale=1.0_dp
# ifdef PERFECT_RESTART
              IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
                gtype=gfactor*r3dvar
              ELSE
                gtype=gfactor*r2dvar
              END IF
# else
              gtype=gfactor*r2dvar
# endif
# ifdef PERFECT_RESTART
              IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
                status=nf_fwrite3d(ng, model, S(ng)%ncid, ifield,       &
     &                             S(ng)%Vid(ifield),                   &
     &                             S(ng)%Rindex, gtype,                 &
     &                             LBi, UBi, LBj, UBj, 1, 2, scale,     &
#  ifdef MASKING
     &                             GRID(ng) % rmask_full,               &
#  endif
     &                             iceField3d)
              ELSE
                status=nf_fwrite2d(ng, model, S(ng)%ncid, ifield,       &
     &                             S(ng)%Vid(ifield),                   &
     &                             S(ng)%Rindex, gtype,                 &
     &                             LBi, UBi, LBj, UBj, scale,           &
#  ifdef MASKING
     &                             GRID(ng) % rmask_full,               &
#  endif
     &                             iceField2d)
              END IF
# else
              status=nf_fwrite2d(ng, model, S(ng)%ncid, ifield,         &
     &                           S(ng)%Vid(ifield),                     &
     &                           S(ng)%Rindex, gtype,                   &
     &                           LBi, UBi, LBj, UBj, scale,             &
#  ifdef MASKING
     &                           GRID(ng) % rmask_full,                 &
#  endif
     &                           iceField2d)
# endif
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
!
  10  FORMAT (/," ICE_WRT_NF90 - error while writing variable '",a,     &
     &        "', time record = ",i0,/,11x,'into file: ',a)
!
      RETURN
      END SUBROUTINE ice_wrt_nf90

# ifdef STATIONS
!
!***********************************************************************
      SUBROUTINE ice_def_station_nf90 (ng, model, ldef, VarOut, S,      &
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
     &  __FILE__//", ice_def_stations_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define prognostic ice variables into specified output NetCDF.
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
!  Define ice model state prognostic variables.
!-----------------------------------------------------------------------
!
      DEFINE : IF (ldef) THEN

        DO nf=1,nIceS
          IF (iSice(nf).gt.0) THEN
            ifield=iSice(nf)
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
        IF (VarOut(idUiER,ng)) THEN
          Vinfo( 1)=Vname(1,idUiER)
          Vinfo( 2)=Vname(2,idUiER)
          Vinfo( 3)=Vname(3,idUiER)
          Vinfo(14)=Vname(4,idUiER)
          Vinfo(16)=Vname(1,idtime)
!
          status=def_var(ng, model, S(ng)%ncid,                         &
     &                   S(ng)%Vid(idUiER),                             &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define isf NORTHward velocity at RHO-points.
!
        IF (VarOut(idViNR,ng)) THEN
          Vinfo( 1)=Vname(1,idViNR)
          Vinfo( 2)=Vname(2,idViNR)
          Vinfo( 3)=Vname(3,idViNR)
          Vinfo(14)=Vname(4,idViNR)
          Vinfo(16)=Vname(1,idtime)
!
          status=def_var(ng, model, S(ng)%ncid,                         &
     &                   S(ng)%Vid(idViNR),                             &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .TRUE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define ice model internal variables.
!-----------------------------------------------------------------------
!
        DO nf=1,nIceF
          IF (iFice(nf).gt.0) THEN
            ifield=iFice(nf)
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
!  Define ice tracer amounts.
!-----------------------------------------------------------------------
!
        IF (Licetrc(ng)) THEN
          DO itrc=NAT+1,NT(ng)
            IF (LtracerIce(itrc, ng)) THEN
              ifield=idStrc(itrc)
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
!  sea-ice variables. Get variable IDs.
!
        DO i=1,n_var
          DO nf=1,nIceS
            IF (iSice(nf).gt.0) THEN
              ifield=iSice(nf)
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                got_var(ifield)=.TRUE.
                S(ng)%Vid(ifield)=var_id(i)
                CYCLE
              END IF
            END IF
          END DO
!
          DO nf=1,nIceF
            IF (iFice(nf).gt.0) THEN
              ifield=iFice(nf)
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                got_var(ifield)=.TRUE.
                S(ng)%Vid(ifield)=var_id(i)
                CYCLE
              END IF
            END IF
          END DO
!
          IF (Licetrc(ng)) THEN
            DO itrc=NAT+1,NT(ng)
              IF (LtracerIce(itrc, ng)) THEN
                ifield=idStrc(itrc)
                IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                  got_var(ifield)=.TRUE.
                  S(ng)%Vid(ifield)=var_id(i)
                  CYCLE
                END IF
              END IF
            END DO
          END IF
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
        DO nf=1,nIceS
          IF (iSice(nf).gt.0) THEN
            ifield=iSice(nf)
            IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              RETURN
            END IF
          END IF
        END DO
!
        DO nf=1,nIceF
          IF (iFice(nf).gt.0) THEN
            ifield=iFice(nf)
            IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              RETURN
            END IF
          END IF
        END DO
!
        IF (Licetrc(ng)) THEN
          DO itrc=NAT+1,NT(ng)
            IF (LtracerIce(itrc, ng)) THEN
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
      END IF QUERY
!
  10  FORMAT (/,' ICE_DEF_STATION_NF90 - unable to find variable: ',    &
     &        a,2x,' in output NetCDF file: ',a)
      RETURN
      END SUBROUTINE ice_def_station_nf90
!
!***********************************************************************
      SUBROUTINE ice_wrt_station_nf90 (ng, model, tile,                 &
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
     &  __FILE__//", ice_wrt_station_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out ice prognostic variables into specified output NetCDF file.
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
!
!-----------------------------------------------------------------------
!  write out sea-ice model state prognostic variables.
!-----------------------------------------------------------------------
!
      DO ifld=1,nIceS
        IF (iSice(ifld).gt.0) THEN
          ifield=iSice(ifld)
          IF (VarOut(ifield,ng)) THEN
            scale=1.0_dp
            CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          scale, ICE(ng)%Si(:,:,IOUT,ifld),       &
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
!  Write ice EASTward and NORTHward velocity component (m/s) at
!  RHO-points.
!
      IF (VarOut(idUiER,ng).and.VarOut(idViNR,ng)) THEN
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
     &                    ICE(ng) % Si(:,:,IOUT,isUice),                &
     &                    ICE(ng) % Si(:,:,IOUT,isVice),                &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUiER, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Ur2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idUiER)), psta,              &
     &                        (/1,S(ng)%Rindex/),                       &
     &                        (/Nstation(ng),1/),                       &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idUiER))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL extract_sta2d (ng, model, Cgrid, idViNR, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Vr2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, S(ng)%name,                    &
     &                        TRIM(Vname(1,idViNR)), psta,              &
     &                        (/1,S(ng)%Rindex/),                       &
     &                        (/Nstation(ng),1/),                       &
     &                        ncid = S(ng)%ncid,                        &
     &                        varid = S(ng)%Vid(idViNR))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
!
!-----------------------------------------------------------------------
!  write out ice model internal variables.
!-----------------------------------------------------------------------
!
      DO ifld=1,nIceF
        IF (iFice(ifld).gt.0) THEN
          ifield=iFice(ifld)
          IF (VarOut(ifield,ng)) THEN
#  ifdef ICE_BULK_FLUXES
            IF ((ifld.eq.icSics).or.(ifld.eq.icSicb).or.                &
     &          (ifld.eq.icSocn)) THEN
              scale=rho0*Cp
            ELSE
              scale=1.0_dp
            END IF
#  else
            scale=1.0_dp
#  endif
            CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          scale, ICE(ng)%Fi(:,:,ifld),            &
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
!  write out ice tracer amount.
!-----------------------------------------------------------------------
!
      IF (Licetrc(ng)) THEN
        ic=0
        DO itrc=NAT+1,NT(ng)
          IF (LtracerIce(itrc, ng)) THEN
            ic=ic+1
            ifield=idStrc(itrc)
            IF (VarOut(ifield,ng)) THEN
              scale=1.0_dp
              CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            scale, ICE(ng)%Ti(:,:,IOUT,ic),       &
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
!
  10  FORMAT (/,' ICE_WRT_STATION_NF90 - error while writing ',         &
     &        "variable '",a,"', time record = ",i0,/,11x,              &
     &        'into file: ',a)
!
      RETURN
      END SUBROUTINE ice_wrt_station_nf90
# endif

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE ice_def_pio (ng, model, ldef, VarOut, S,               &
     &                        t2dgrd, u2dgrd, v2dgrd,                   &
     &                        sr2dgrd, ru2dgrd, rv2dgrd)
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
     &  t2dgrd(:), u2dgrd(:), v2dgrd(:),                                &
     &  sr2dgrd(:), ru2dgrd(:), rv2dgrd(:)
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
      integer, pointer :: icegrd(:)
!
      real(r8) :: Aval(6)
!
      character (len=256) :: ncname
      character (len=MaxLen) :: Vinfo(Natt)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ice_def_pio"
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
!
!-----------------------------------------------------------------------
!  Define ice model state prognostic variables.
!-----------------------------------------------------------------------
!
        DO nf=1,nIceS
          IF (iSice(nf).gt.0) THEN
            ifield=iSice(nf)
            IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
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
              S(ng)%pioVar(ifield)%dkind=vtype
!
              SELECT CASE (nf)
                CASE (isUice)
#  ifdef PERFECT_RESTART
                  IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
                    icegrd => ru2dgrd
                  ELSE
                    icegrd => u2dgrd
                  END IF
#  else
                  icegrd => u2dgrd
#  endif
#  if defined WRITE_WATER && defined MASKING
#   if !defined WET_DRY && defined PERFECT_RESTART
                  IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
                    Vinfo(24)='_FillValue'
                    Aval(6)=spval
                  END IF
#   else
                  Vinfo(20)='mask_u'
#   endif
#  endif
                  S(ng)%pioVar(ifield)%gtype=u2dvar
                CASE (isVice)
#  ifdef PERFECT_RESTART
                  IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
                    icegrd => rv2dgrd
                  ELSE
                    icegrd => v2dgrd
                  END IF
#  else
                  icegrd => v2dgrd
#  endif
#  if defined WRITE_WATER && defined MASKING
#   if !defined WET_DRY && defined PERFECT_RESTART
                  IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
                    Vinfo(24)='_FillValue'
                    Aval(6)=spval
                  END IF
#   else
                  Vinfo(20)='mask_v'
#   endif
#  endif
                  S(ng)%pioVar(ifield)%gtype=v2dvar
                CASE DEFAULT
                  icegrd => t2dgrd
#  if defined WRITE_WATER && defined MASKING
#   if !defined WET_DRY && defined PERFECT_RESTART
                  IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
                    Vinfo(24)='_FillValue'
                    Aval(6)=spval
                  END IF
#   else
                  Vinfo(20)='mask_rho'
#   endif
#  endif
                  S(ng)%pioVar(ifield)%gtype=r2dvar
              END SELECT
!
              status=def_var(ng, model, S(ng)%pioFile,                  &
     &                       S(ng)%pioVar(ifield)%vd,                   &
     &                       vtype, nvd, icegrd, Aval, Vinfo, ncname)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
!
        IF (FileH.ne.ABS(RST(ng)%pioFile%fh)) THEN
!
!  Reset icegrid dimension info.
!
          icegrd => t2dgrd
          nvd=nvd3
!
!  Define ice EASTward velocity at RHO-points.
!
          IF (VarOut(idUiER,ng)) THEN
            Vinfo( 1)=Vname(1,idUiER)
            Vinfo( 2)=Vname(2,idUiER)
            Vinfo( 3)=Vname(3,idUiER)
            Vinfo(14)=Vname(4,idUiER)
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(21)=Vname(6,idUiER)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idUiER,ng),r8)
            S(ng)%pioVar(idUiER)%dkind=PIO_FOUT
            S(ng)%pioVar(idUiER)%gtype=r2dvar
!
            status=def_var(ng, model, S(ng)%pioFile,                    &
     &                     S(ng)%pioVar(idUiER)%vd,                     &
     &                     vtype, nvd, icegrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!  Define ice NORTHward velocity at RHO-points.
!
          IF (VarOut(idViNR,ng)) THEN
            Vinfo( 1)=Vname(1,idViNR)
            Vinfo( 2)=Vname(2,idViNR)
            Vinfo( 3)=Vname(3,idViNR)
            Vinfo(14)=Vname(4,idViNR)
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(21)=Vname(6,idViNR)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idViNR,ng),r8)
            S(ng)%pioVar(idViNR)%dkind=PIO_FOUT
            S(ng)%pioVar(idViNR)%gtype=r2dvar
!
            status=def_var(ng, model, S(ng)%pioFile,                    &
     &                     S(ng)%pioVar(idViNR)%vd,                     &
     &                     vtype, nvd, icegrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
!
!-----------------------------------------------------------------------
!  Define ice model internal variables.
!-----------------------------------------------------------------------
!
        IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
          icegrd => sr2dgrd
        ELSE
          icegrd => t2dgrd
        END IF
        nvd=nvd3
!
        DO nf=1,nIceF
          IF (iFice(nf).gt.0) THEN
            ifield=iFice(nf)
            IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
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
              IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
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
              S(ng)%pioVar(ifield)%dkind=vtype
              S(ng)%pioVar(ifield)%gtype=r2dvar
!
              SELECT CASE (nf)
                CASE (icWdiv)
                  Vinfo(11)='increase ice thickness'
                  Vinfo(12)='decrease ice concentration'
                CASE (icW_ai)
                  Vinfo(11)='freezing'
                  Vinfo(12)='melting'
                CASE (icW_ao, icW_io)
                  Vinfo(11)='melting'
                  Vinfo(12)='freezing'
              END SELECT
!
              status=def_var(ng, model, S(ng)%pioFile,                  &
     &                       S(ng)%pioVar(ifield)%vd,                   &
     &                       vtype, nvd, icegrd, Aval, Vinfo, ncname)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Define ice tracer amounts.
!-----------------------------------------------------------------------
!
!  Reset icegrid dimension info.
!
        icegrd => t2dgrd
        IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
#  ifdef PERFECT_RESTART
          nvd=nvd4
#  else
          nvd=nvd3
#  endif
        ELSE
          nvd=nvd3
        END IF
        IF (Licetrc(ng)) THEN
          DO itrc=NAT+1,NT(ng)
            IF (LtracerIce(itrc, ng)) THEN
              ifield=idStrc(itrc)
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
#  if defined WRITE_WATER && defined MASKING
#   if !defined WET_DRY && defined PERFECT_RESTART
                IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
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
                S(ng)%pioVar(ifield)%dkind=vtype
                S(ng)%pioVar(ifield)%gtype=r2dvar
!
                status=def_var(ng, model, S(ng)%pioFile,                &
     &                         S(ng)%pioVar(ifield)%vd,                 &
     &                         vtype, nvd, icegrd, Aval, Vinfo, ncname)
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
              END IF
            END IF
          END DO
        END IF
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
!  sea-ice variables. Get variable IDs.
!
        IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
          vtype=PIO_FRST
        ELSE
          vtype=PIO_FOUT
        END IF
!
        DO i=1,n_var
          DO nf=1,nIceS
            IF (iSice(nf).gt.0) THEN
              ifield=iSice(nf)
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                got_var(ifield)=.TRUE.
                S(ng)%pioVar(ifield)%vd=var_desc(i)
                S(ng)%pioVar(ifield)%dkind=vtype
!
                SELECT CASE (nf)
                  CASE (isUice)
                    S(ng)%pioVar(ifield)%gtype=u2dvar
                  CASE (isVice)
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
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUiER))) THEN
            got_var(idUiER)=.TRUE.
            S(ng)%pioVar(idUiER)%vd=var_desc(i)
            S(ng)%pioVar(idUiER)%dkind=PIO_FOUT
            S(ng)%pioVar(idUiER)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idViNR))) THEN
            got_var(idViNR)=.TRUE.
            S(ng)%pioVar(idViNR)%vd=var_desc(i)
            S(ng)%pioVar(idViNR)%dkind=PIO_FOUT
            S(ng)%pioVar(idViNR)%gtype=r2dvar
          END IF
!
          DO nf=1,nIceF
            IF (iFice(nf).gt.0) THEN
              ifield=iFice(nf)
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
          IF (Licetrc(ng)) THEN
            DO itrc=NAT+1,NT(ng)
              IF (LtracerIce(itrc, ng)) THEN
                ifield=idStrc(itrc)
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
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
        DO nf=1,nIceS
          IF (iSice(nf).gt.0) THEN
            ifield=iSice(nf)
            IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              RETURN
            END IF
          END IF
        END DO
!
        IF (.not.got_var(idUiER).and.VarOut(idUiER,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUiER)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idViNR).and.VarOut(idViNR,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idViNR)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
!
        DO nf=1,nIceF
          IF (iFice(nf).gt.0) THEN
            ifield=iFice(nf)
            IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
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
        IF (Licetrc(ng)) THEN
          DO itrc=NAT+1,NT(ng)
            IF (LtracerIce(itrc, ng)) THEN
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
      END IF QUERY
!
  10  FORMAT (/,' ICE_DEF_PIO - unable to find variable: ',a,2x,        &
     &        ' in output NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE ice_def_pio
!
!***********************************************************************
     SUBROUTINE ice_wrt_pio (ng, model, tile,                           &
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
      integer :: FileH, ifield, ifld, ic, itrc
!
      real(dp) :: scale
!
      real(r8), pointer :: iceField2d(:,:)
      real(r8), pointer :: iceField3d(:,:,:)
      real(r8), pointer :: iceMask(:,:)
!
      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ice_wrt_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out ice prognostic variables into specified output NetCDF file.
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
!
!-----------------------------------------------------------------------
!  write out sea-ice model state prognostic variables.
!-----------------------------------------------------------------------
!
      DO ifld=1,nIceS
        IF (iSice(ifld).gt.0) THEN
          ifield=iSice(ifld)
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
              iceField2d => ICE(ng) % Si(LBi:UBi,LBj:UBj,IOUT,ifld)
            ELSE IF ((model.eq.iNLM).and.                               &
     &               (FileH.eq.ABS(RST(ng)%pioFile%fh))) THEN
#  ifdef PERFECT_RESTART
              iceField3d => ICE(ng) % Si(LBi:UBi,LBj:UBj,1:2,ifld)
#  else
              iceField2d => ICE(ng) % Si(LBi:UBi,LBj:UBj,IOUT,ifld)
#  endif
#  ifdef AVERAGES
            ELSE IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
              iceField2d => ICE_SAVG(ifld,ng) % var(LBi:UBi,LBj:UBj)
#  endif
            END IF
!
            SELECT CASE (ifld)
              CASE (isUice)
#  ifdef PERFECT_RESTART
                IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
                  IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                    ioDesc => ioDesc_dp_rubar(ng)
                  ELSE
                    ioDesc => ioDesc_sp_rubar(ng)
                  END IF
                ELSE
                  IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                    ioDesc => ioDesc_dp_u2dvar(ng)
                  ELSE
                    ioDesc => ioDesc_sp_u2dvar(ng)
                  END IF
                END IF
#  else
                IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                  ioDesc => ioDesc_dp_u2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_u2dvar(ng)
                END IF
#  endif
#  ifdef MASKING
                iceMask => GRID(ng) % umask_full
#  endif
              CASE (isVice)
#  ifdef PERFECT_RESTART
                IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
                  IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                    ioDesc => ioDesc_dp_rvbar(ng)
                  ELSE
                    ioDesc => ioDesc_sp_rvbar(ng)
                  END IF
                ELSE
                  IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                    ioDesc => ioDesc_dp_v2dvar(ng)
                  ELSE
                    ioDesc => ioDesc_sp_v2dvar(ng)
                  END IF
                END IF
#  else
                IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                  ioDesc => ioDesc_dp_v2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_v2dvar(ng)
                END IF
#  endif
#  ifdef MASKING
                iceMask => GRID(ng) % vmask_full
#  endif
              CASE DEFAULT
#  ifdef PERFECT_RESTART
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
#  else
                IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                  ioDesc => ioDesc_dp_r2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_r2dvar(ng)
                END IF
#  endif
#  ifdef MASKING
                iceMask => GRID(ng) % rmask_full
#  endif
            END SELECT
            scale=1.0_dp
!
#  ifdef PERFECT_RESTART
            IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
              status=nf_fwrite3d(ng, model, S(ng)%pioFile, ifield,      &
     &                           S(ng)%pioVar(ifield), S(ng)%Rindex,    &
     &                           ioDesc,                                &
     &                           LBi, UBi, LBj, UBj, 1, 2, scale,       &
#   ifdef MASKING
     &                           iceMask,                               &
#   endif
     &                           iceField3d)
            ELSE
              status=nf_fwrite2d(ng, model, S(ng)%pioFile, ifield,      &
     &                           S(ng)%pioVar(ifield), S(ng)%Rindex,    &
     &                           ioDesc,                                &
     &                           LBi, UBi, LBj, UBj, scale,             &
#   ifdef MASKING
     &                           iceMask,                               &
#   endif
     &                           iceField2d)
            END IF
#  else
            status=nf_fwrite2d(ng, model, S(ng)%pioFile, ifield,        &
     &                         S(ng)%pioVar(ifield), S(ng)%Rindex,      &
     &                         ioDesc,                                  &
     &                         LBi, UBi, LBj, UBj, scale,               &
#   ifdef MASKING
     &                         iceMask,                                 &
#   endif
     &                         iceField2d)
#  endif
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
!  Write ice EASTward and NORTHward velocity component (m/s) at
!  RHO-points.
!
      IF (VarOut(idUiER,ng).and.VarOut(idViNR,ng)) THEN
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
     &                    ICE(ng) % Si(:,:,IOUT,isUice),                &
     &                    ICE(ng) % Si(:,:,IOUT,isVice),                &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        IF (S(ng)%pioVar(idUiER)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idUiER,            &
     &                     S(ng)%pioVAR(idUiER), S(ng)%Rindex,          &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     Ur2d)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUiER)),          &
     &                                  S(ng)%Rindex, TRIM(S(ng)%name)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        IF (S(ng)%pioVar(idViNR)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idViNR,            &
     &                     S(ng)%pioVar(idViNR), S(ng)%Rindex,          &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     Vr2d)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idViNR)),          &
     &                                  S(ng)%Rindex, TRIM(S(ng)%name)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
!
!-----------------------------------------------------------------------
!  write out ice model internal variables.
!-----------------------------------------------------------------------
!
      DO ifld=1,nIceF
        IF (iFice(ifld).gt.0) THEN
          ifield=iFice(ifld)
          IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
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
     &          ((FileH.eq.ABS(HIS(ng)%pioFile%fh)).or.                 &
     &           (FileH.eq.ABS(QCK(ng)%pioFile%fh)).or.                 &
     &           (FileH.eq.ABS(RST(ng)%pioFile%fh)))) THEN
              iceField2d => ICE(ng) % Fi(LBi:UBi,LBj:UBj,ifld)
#  ifdef AVERAGES
            ELSE IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
              iceField2d => ICE_FAVG(ifld,ng) % var(LBi:UBi,LBj:UBj)
#  endif
            END IF
!
#  ifdef ICE_BULK_FLUXES
            IF ((ifld.eq.icSics).or.(ifld.eq.icSicb).or.                &
     &          (ifld.eq.icSocn)) THEN
              scale=rho0*Cp
            ELSE
              scale=1.0_dp
            END IF
#  else
            scale=1.0_dp
#  endif
            IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
            status=nf_fwrite2d(ng, model, S(ng)%pioFile, ifield,        &
     &                         S(ng)%pioVar(ifield), S(ng)%Rindex,      &
     &                         ioDesc,                                  &
     &                         LBi, UBi, LBj, UBj, scale,               &
#  ifdef MASKING
     &                         GRID(ng) % rmask_full,                   &
#  endif
     &                         iceField2d)
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
!-----------------------------------------------------------------------
!  write out ice tracer amount at RHO-points.
!-----------------------------------------------------------------------
!
      IF (Licetrc(ng)) THEN
        ic=0
        DO itrc=NAT+1,NT(ng)
          IF (LtracerIce(itrc, ng)) THEN
            ic=ic+1
            ifield=idStrc(itrc)
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
                iceField2d => ICE(ng) % Ti(LBi:UBi,LBj:UBj,IOUT,ic)
              ELSE IF ((model.eq.iNLM).and.                             &
     &                 (FileH.eq.ABS(RST(ng)%pioFile%fh))) THEN
#  ifdef PERFECT_RESTART
                iceField3d => ICE(ng) % Ti(LBi:UBi,LBj:UBj,1:2,ic)
#  else
                iceField2d => ICE(ng) % Ti(LBi:UBi,LBj:UBj,IOUT,ic)
#  endif
#  ifdef AVERAGES
              ELSE IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
                iceField2d => ICE_TAVG(itrc,ng) % var(LBi:UBi,LBj:UBj)
#  endif
              END IF
!
              scale=1.0_dp
#  ifdef PERFECT_RESTART
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
#  else
              IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
#  endif
#  ifdef PERFECT_RESTART
              IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
                status=nf_fwrite3d(ng, model, S(ng)%pioFile, ifield,    &
     &                             S(ng)%pioVar(ifield), S(ng)%Rindex,  &
     &                             ioDesc,                              &
     &                             LBi, UBi, LBj, UBj, 1, 2, scale,     &
#   ifdef MASKING
     &                             GRID(ng) % rmask_full,               &
#   endif
     &                             iceField3d)
              ELSE
                status=nf_fwrite2d(ng, model, S(ng)%pioFile, ifield,    &
     &                             S(ng)%pioVar(ifield), S(ng)%Rindex,  &
     &                             ioDesc,                              &
     &                             LBi, UBi, LBj, UBj, scale,           &
#   ifdef MASKING
     &                             GRID(ng) % rmask_full,               &
#   endif
     &                             iceField2d)
              END IF
#  else
              status=nf_fwrite2d(ng, model, S(ng)%pioFile, ifield,      &
     &                           S(ng)%pioVar(ifield), S(ng)%Rindex,    &
     &                           ioDesc,                                &
     &                           LBi, UBi, LBj, UBj, scale,             &
#   ifdef MASKING
     &                           GRID(ng) % rmask_full,                 &
#   endif
     &                           iceField2d)
#  endif
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
!
  10  FORMAT (/," ICE_WRT - error while writing variable '",a,          &
     &        "', time record = ",i0,/,11x,'into file: ',a)
!
      RETURN
      END SUBROUTINE ice_wrt_pio
# endif
#endif
      END MODULE ice_output_mod
