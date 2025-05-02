      SUBROUTINE read_IcePar (model, inp, out, Lwrite)
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2024 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This routine reads and reports ice model input parameters.          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ice
      USE mod_ncparam
      USE mod_scalars
      USE inp_decode_mod
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
!
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      integer :: Npts, Nval, ng, status
      integer :: ioff, nline, igrid
#ifdef AVERAGES
      integer :: nf, ns
#endif
      integer :: i, ic, itrc
      logical, allocatable :: Lswitch(:)
#ifdef T_PASSIVE
      real(r8), allocatable :: Rtracer(:,:)
      logical, allocatable :: Linert(:,:)
#endif
!
      real(r8) :: Rvalue(1)
      real(r8), dimension(200) :: Rval
!
      character (len=40) :: KeyWord
      character (len=256) :: line
      character (len=256), dimension(200) :: Cval

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
      igrid=1                            ! nested grid counter
      nline=0                            ! LBC multi-line counter
!
!-----------------------------------------------------------------------
!  Allocate ice model parameters that depend on Ngrids.
!-----------------------------------------------------------------------
!
      IF (.not.allocated(Lswitch)) THEN
        allocate (Lswitch(Ngrids))
      END IF
#ifdef T_PASSIVE
      IF (.not.allocated(Rtracer)) THEN
        allocate (Rtracer(NPT,Ngrids))
      END IF
      IF (.not.allocated(Linert)) THEN
        allocate (Linert(NPT,Ngrids))
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Read in ice model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.true.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('LBC(isAice)')
              ioff = ice_ic_brylast + isAice
              Npts=load_lbc(Nval, Cval, line, nline, ioff, igrid,       &
     &                      0, 0, Vname(1,idAice), LBC)
            CASE ('LBC(isHice)')
              ioff = ice_ic_brylast + isHice
              Npts=load_lbc(Nval, Cval, line, nline, ioff, igrid,       &
     &                      0, 0, Vname(1,idHice), LBC)
            CASE ('LBC(isHsno)')
              ioff = ice_ic_brylast + isHsno
              Npts=load_lbc(Nval, Cval, line, nline, ioff, igrid,       &
     &                      0, 0, Vname(1,idHsno), LBC)
            CASE ('LBC(isTice)')
              ioff = ice_ic_brylast + isTice
              Npts=load_lbc(Nval, Cval, line, nline, ioff, igrid,       &
     &                      0, 0, Vname(1,idTice), LBC)
            CASE ('LBC(isHmel)')
              ioff = ice_ic_brylast + isHmel
              Npts=load_lbc(Nval, Cval, line, nline, ioff, igrid,       &
     &                      0, 0, Vname(1,idHmel), LBC)
            CASE ('LBC(isIage)')
              ioff = ice_ic_brylast + isIage
              Npts=load_lbc(Nval, Cval, line, nline, ioff, igrid,       &
     &                      0, 0, Vname(1,idIage), LBC)
            CASE ('LBC(isISxx)')
              ioff = ice_ic_brylast + isISxx
              Npts=load_lbc(Nval, Cval, line, nline, ioff, igrid,       &
     &                      0, 0, Vname(1,idISxx), LBC)
            CASE ('LBC(isISxy)')
              ioff = ice_ic_brylast + isISxy
              Npts=load_lbc(Nval, Cval, line, nline, ioff, igrid,       &
     &                      0, 0, Vname(1,idISxy), LBC)
            CASE ('LBC(isISyy)')
              ioff = ice_ic_brylast + isISyy
              Npts=load_lbc(Nval, Cval, line, nline, ioff, igrid,       &
     &                      0, 0, Vname(1,idISyy), LBC)
            CASE ('LBC(isUice)')
              ioff = ice_ic_brylast + isUice
              Npts=load_lbc(Nval, Cval, line, nline, ioff, igrid,       &
     &                      0, 0, Vname(1,idUice), LBC)
            CASE ('LBC(isVice)')
              ioff = ice_ic_brylast + isVice
              Npts=load_lbc(Nval, Cval, line, nline, ioff, igrid,       &
     &                      0, 0, Vname(1,idVice), LBC)
            CASE ('nEVP')
              Npts=load_i(Nval, Rval, Ngrids, nEVP)
            CASE ('AirRho')
              Npts=load_r(Nval, Rval, Ngrids, AirRho)
            CASE ('IceRho')
              Npts=load_r(Nval, Rval, Ngrids, IceRho)
            CASE ('SnowDryRho')
              Npts=load_r(Nval, Rval, Ngrids, SnowDryRho)
            CASE ('SnowWetRho')
              Npts=load_r(Nval, Rval, Ngrids, SnowWetRho)
            CASE ('Cd_ai')
              Npts=load_r(Nval, Rval, Ngrids, Cd_ai)
            CASE ('Cd_io')
              Npts=load_r(Nval, Rval, Ngrids, Cd_io)
            CASE ('Astrength')
              Npts=load_r(Nval, Rval, Ngrids, Astrength)
            CASE ('Pstar')
              Npts=load_r(Nval, Rval, Ngrids, Pstar)
            CASE ('Pstar_e')
              Npts=load_r(Nval, Rval, Ngrids, Pstar_e)
            CASE ('zetaMin')
              Npts=load_r(Nval, Rval, Ngrids, zetaMin)
            CASE ('zetaMax')
              Npts=load_r(Nval, Rval, Ngrids, zetaMax)
            CASE ('ellip_sq')
              Npts=load_r(Nval, Rval, Ngrids, ellip_sq)
            CASE ('min_ai')
              Npts=load_r(Nval, Rval, Ngrids, min_ai)
            CASE ('max_ai')
              Npts=load_r(Nval, Rval, Ngrids, max_ai)
            CASE ('min_hi')
              Npts=load_r(Nval, Rval, Ngrids, min_hi)
            CASE ('max_hmelt')
              Npts=load_r(Nval, Rval, Ngrids, max_hmelt)
            CASE ('stressAng')
              Npts=load_r(Nval, Rval, Ngrids, stressAng)
              DO ng=1,Ngrids
                 stressAng(ng) = stressAng(ng)*deg2rad
              END DO
            CASE ('ice_emiss')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              ice_emiss=Rvalue(1)
            CASE ('spec_heat_air')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              spec_heat_air=Rvalue(1)
            CASE ('trans_coeff')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              trans_coeff=Rvalue(1)
            CASE ('sublimation')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              sublimation=Rvalue(1)
#ifdef T_PASSIVE
            CASE ('LtracerIce')
              Npts=load_l(Nval, Cval, NPT, Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=inert(i)
                  LtracerIce(itrc,ng)=Linert(i,ng)
                END DO
              END DO
              DO ng=1,Ngrids
                IF (ANY(LtracerIce(:,ng))) THEN
                  Licetrc(ng)=.TRUE.
                  ic=0
                  DO itrc=NAT+1,NT(ng)
                    IF (LtracerIce(itrc,ng)) THEN
                      ic=ic+1
                    END IF
                  END DO
                  NTICE(ng)=ic
                END IF
              END DO
            CASE ('LsourceIce')
              Npts=load_l(Nval, Cval, NPT, Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=inert(i)
                  LsourceIce(itrc,ng)=Linert(i,ng)
                END DO
              END DO
              DO ng=1,Ngrids
                IF (ANY(LsourceIce(:,ng))) THEN
                  Licesrc(ng)=.TRUE.
                END IF
              END DO
            CASE ('trcIce')
              Npts=load_r(Nval, Rval, NPT, Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=inert(i)
                  trcIce(itrc,ng)=Rtracer(i,ng)
                END DO
              END DO
            CASE ('LsourceFrz')
              Npts=load_l(Nval, Cval, NPT, Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=inert(i)
                  LsourceFrz(itrc,ng)=Linert(i,ng)
                END DO
              END DO
              DO ng=1,Ngrids
                IF (ANY(LsourceFrz(:,ng))) THEN
                  Lfrzsrc(ng)=.TRUE.
                END IF
              END DO
            CASE ('trcFrz')
              Npts=load_r(Nval, Rval, NPT, Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=inert(i)
                  trcFrz(itrc,ng)=Rtracer(i,ng)
                END DO
              END DO
#endif
            CASE ('Hout(idUice)')
              IF (idUice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idUice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idUice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idVice)')
              IF (idVice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idVice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idVice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idUiER)')
              IF (idUiER.eq.0) THEN
                IF (Master) WRITE (out,80) 'idUiER'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idUiER,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idViNR)')
              IF (idViNR.eq.0) THEN
                IF (Master) WRITE (out,80) 'idViNR'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idViNR,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idAice)')
              IF (idAice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idAice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idAice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idIage)')
              IF (idIage.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIage'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idIage,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idHice)')
              IF (idHice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idHice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idHice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idHmel)')
              IF (idHmel.eq.0) THEN
                IF (Master) WRITE (out,80) 'idHmel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idHmel,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idHsno)')
              IF (idHsno.eq.0) THEN
                IF (Master) WRITE (out,80) 'idHsno'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idHsno,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idTice)')
              IF (idTice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idTice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idTice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idISxx)')
              IF (idISxx.eq.0) THEN
                IF (Master) WRITE (out,80) 'idISxx'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idISxx,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idISxy)')
              IF (idISxy.eq.0) THEN
                IF (Master) WRITE (out,80) 'idISxy'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idISxy,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idISyy)')
              IF (idISyy.eq.0) THEN
                IF (Master) WRITE (out,80) 'idISyy'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idISyy,1:Ngrids)=Lswitch(1:Ngrids)
#if defined ICE_BIO
            CASE ('Hout(idIphy)')
              IF (idIphy.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIphy'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idIphy,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idINO3)')
              IF (idINO3.eq.0) THEN
                IF (Master) WRITE (out,80) 'idINO3'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idINO3,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idINH4)')
              IF (idINH4.eq.0) THEN
                IF (Master) WRITE (out,80) 'idINH4'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idINH4,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idIlog)')
              IF (idIlog.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIlog'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idIlog,1:Ngrids)=Lswitch(1:Ngrids)
#endif
            CASE ('Hout(idIsst)')
              IF (idIsst.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIsst'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idIsst,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idIOmf)')
              IF (idIOmf.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIOmf'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idIOmf,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idIOfv)')
              IF (idIOfv.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIOfv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idIOfv,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idIOmt)')
              IF (idIOmt.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIOmt'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idIOmt,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idS0mk)')
              IF (idS0mk.eq.0) THEN
                IF (Master) WRITE (out,80) 'idS0mk'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idS0mk,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idT0mk)')
              IF (idT0mk.eq.0) THEN
                IF (Master) WRITE (out,80) 'idT0mk'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idT0mk,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWdiv)')
              IF (idWdiv.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWdiv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWdiv,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idW_fr)')
              IF (idW_fr.eq.0) THEN
                IF (Master) WRITE (out,80) 'idW_fr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idW_fr,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idW_ai)')
              IF (idW_ai.eq.0) THEN
                IF (Master) WRITE (out,80) 'idW_ai'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idW_ai,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idW_ao)')
              IF (idW_ao.eq.0) THEN
                IF (Master) WRITE (out,80) 'idW_ao'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idW_ao,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idW_io)')
              IF (idW_io.eq.0) THEN
                IF (Master) WRITE (out,80) 'idW_io'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idW_io,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idW_ro)')
              IF (idW_ro.eq.0) THEN
                IF (Master) WRITE (out,80) 'idW_ro'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idW_ro,1:Ngrids)=Lswitch(1:Ngrids)
#if defined ICE_BULK_FLUXES
            CASE ('Hout(idAlbi)')
              IF (idAlbi.eq.0) THEN
                IF (Master) WRITE (out,80) 'idAlbi'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idAlbi,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idAlbo)')
              IF (idAlbo.eq.0) THEN
                IF (Master) WRITE (out,80) 'idAlbo'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idAlbo,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idFrsw)')
              IF (idFrsw.eq.0) THEN
                IF (Master) WRITE (out,80) 'idFrsw'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idFrsw,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idSics)')
              IF (idSics.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSics'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idSics,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idSicb)')
              IF (idSicb.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSicb'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idSicb,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idSocn)')
              IF (idSocn.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSocn'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idSocn,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idRfac)')
              IF (idRfac.eq.0) THEN
                IF (Master) WRITE (out,80) 'idRfac'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idRfac,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idSnof)')
              IF (idSnof.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSnof'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idSnof,1:Ngrids)=Lswitch(1:Ngrids)
#endif
#ifdef T_PASSIVE
            CASE ('Hout(idStrc)')
              IF (MAXVAL(idStrc).eq.0) THEN
                IF (Master) WRITE (out,80) 'idStrc'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NPT, Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=inert(i)
                  IF (LtracerIce(itrc,ng)) THEN
                    Hout(idStrc(itrc),ng)=Linert(i,ng)
                  END IF
                END DO
              END DO
#endif
            CASE ('Qout(idUice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idUice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idVice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idVice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idUiER)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idUiER,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idViNR)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idViNR,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idAice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idAice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idIage)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idIage,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idHice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idHice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idHmel)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idHmel,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idHsno)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idHsno,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idTice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idTice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idISxx)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idISxx,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idISxy)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idISxy,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idISyy)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idISyy,1:Ngrids)=Lswitch(1:Ngrids)
#if defined ICE_BIO
            CASE ('Qout(idIphy)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idIphy,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idINO3)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idINO3,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idINH4)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idINH4,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idIlog)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idIlog,1:Ngrids)=Lswitch(1:Ngrids)
#endif
            CASE ('Qout(idIsst)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idIsst,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idIOmf)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idIOmf,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idIOfv)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idIOfv,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idIOmt)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idIOmt,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idS0mk)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idS0mk,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idT0mk)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idT0mk,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idWdiv)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idWdiv,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idW_fr)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idW_fr,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idW_ai)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idW_ai,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idW_ao)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idW_ao,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idW_io)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idW_io,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idW_ro)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idW_ro,1:Ngrids)=Lswitch(1:Ngrids)
#if defined ICE_BULK_FLUXES
            CASE ('Qout(idAlbi)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idAlbi,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idAlbo)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idAlbo,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idFrsw)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idFrsw,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idSics)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idSics,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idSicb)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idSicb,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idSocn)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idSocn,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idRfac)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idRfac,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idSnof)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idSnof,1:Ngrids)=Lswitch(1:Ngrids)
#endif
#ifdef T_PASSIVE
            CASE ('Qout(idStrc)')
              Npts=load_l(Nval, Cval, NPT, Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=inert(i)
                  IF (LtracerIce(itrc,ng)) THEN
                    Qout(idStrc(itrc),ng)=Linert(i,ng)
                  END IF
                END DO
              END DO
#endif
#ifdef AVERAGES
            CASE ('Aout(idUice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idUice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idVice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idVice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idUiER)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idUiER,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idViNR)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idViNR,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idAice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idAice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idIage)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idIage,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idHice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idHice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idHmel)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idHmel,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idHsno)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idHsno,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idTice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idTice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idISxx)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idISxx,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idISxy)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idISxy,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idISyy)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idISyy,1:Ngrids)=Lswitch(1:Ngrids)
# if defined ICE_BIO
            CASE ('Aout(idIphy)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idIphy,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idINO3)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idINO3,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idINH4)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idINH4,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idIlog)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idIlog,1:Ngrids)=Lswitch(1:Ngrids)
# endif
            CASE ('Aout(idIsst)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idIsst,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idIOmf)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idIOmf,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idIOfv)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idIOfv,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idIOmt)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idIOmt,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idS0mk)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idS0mk,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idT0mk)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idT0mk,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idWdiv)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idWdiv,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idW_fr)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idW_fr,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idW_ai)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idW_ai,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idW_ao)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idW_ao,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idW_io)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idW_io,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idW_ro)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idW_ro,1:Ngrids)=Lswitch(1:Ngrids)
# if defined ICE_BULK_FLUXES
            CASE ('Aout(idAlbi)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idAlbi,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idAlbo)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idAlbo,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idFrsw)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idFrsw,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idSics)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idSics,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idSicb)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idSicb,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idSocn)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idSocn,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idRfac)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idRfac,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idSnof)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idSnof,1:Ngrids)=Lswitch(1:Ngrids)
# endif
# ifdef T_PASSIVE
            CASE ('Aout(idStrc)')
              Npts=load_l(Nval, Cval, NPT, Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=inert(i)
                  IF (LtracerIce(itrc,ng)) THEN
                    Aout(idStrc(itrc),ng)=Linert(i,ng)
                  END IF
                END DO
              END DO
# endif
#endif
#ifdef STATIONS
           CASE ('Sout(idUice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idUice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idVice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idVice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idUiER)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idUiER,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idViNR)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idViNR,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idAice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idAice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idIage)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idIage,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idHice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idHice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idHmel)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idHmel,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idHsno)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idHsno,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idTice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idTice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idISxx)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idISxx,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idISxy)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idISxy,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idISyy)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idISyy,1:Ngrids)=Lswitch(1:Ngrids)
# if defined ICE_BIO
            CASE ('Sout(idIphy)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idIphy,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idINO3)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idINO3,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idINH4)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idINH4,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idIlog)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idIlog,1:Ngrids)=Lswitch(1:Ngrids)
# endif
            CASE ('Sout(idIsst)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idIsst,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idIOmf)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idIOmf,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idIOfv)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idIOfv,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idIOmt)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idIOmt,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idS0mk)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idS0mk,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idT0mk)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idT0mk,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idWdiv)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idWdiv,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idW_fr)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idW_fr,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idW_ai)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idW_ai,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idW_ao)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idW_ao,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idW_io)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idW_io,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idW_ro)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idW_ro,1:Ngrids)=Lswitch(1:Ngrids)
# if defined ICE_BULK_FLUXES
            CASE ('Sout(idAlbi)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idAlbi,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idAlbo)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idAlbo,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idFrsw)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idFrsw,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idSics)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idSics,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idSicb)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idSicb,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idSocn)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idSocn,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idRfac)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idRfac,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idSnof)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idSnof,1:Ngrids)=Lswitch(1:Ngrids)
# endif
# ifdef T_PASSIVE
            CASE ('Sout(idStrc)')
              Npts=load_l(Nval, Cval, NPT, Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=inert(i)
                  IF (LtracerIce(itrc,ng)) THEN
                    Sout(idStrc(itrc),ng)=Linert(i,ng)
                  END IF
                END DO
              END DO
# endif
#endif
          END SELECT
        END IF
      END DO
  10  IF (Master) WRITE (out,30) line
      exit_flag=4
      RETURN
  20  CLOSE (inp)

! Set ice time step to ocean time step

      DO ng = 1,Ngrids
        dtice(ng) = dt(ng)
      END DO
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          WRITE (out,40) ng
          WRITE (out,50) nEVP(ng), 'nEVP',                              &
     &          'Number of elastic steps per plastic step in EVP.'
          WRITE (out,60) AirRho(ng), 'AirRho',                          &
     &          'Air density (kg/m3).'
          WRITE (out,60) IceRho(ng), 'IceRho',                          &
     &          'Density of sea ice (kg/m3).'
          WRITE (out,60) SnowDryRho(ng), 'SnowDryRho',                  &
     &          'Dry snow density (kg/m3).'
          WRITE (out,60) SnowWetRho(ng), 'SnowWetRho',                  &
     &          'Wet snow density (kg/m3).'
          WRITE (out,60) Cd_ai(ng), 'Cd_ai',                            &
     &          'Air-Ice drag coefficient (nondimensional).'
          WRITE (out,60) Cd_io(ng), 'Cd_io',                            &
     &          'Ice-Ocean drag coefficient (nondimensional).'
          WRITE (out,60) Astrength(ng), 'Astrength',                    &
     &          'Ice strength exponential weighting (nondimensional).'
          WRITE (out,60) Pstar(ng), 'Pstar',                            &
     &          'Ice strength coefficient (N/m2).'
          WRITE (out,60) Pstar_e(ng), 'Pstar_e',                        &
     &          'Ice strength coefficient in quadratic form (N/m).'
          WRITE (out,60) ZetaMin(ng), 'zetaMin',                        &
     &          'Minimum ice shear strength (N/m2) limiter.'
          WRITE (out,60) ZetaMax(ng), 'zetaMax',                        &
     &          'Maximum ice shear strength (N/m2) limiter.'
          WRITE (out,60) ellip_sq(ng), 'ellip_sq',                      &
     &          'Ellipticity squared of yield curve (nondimensional).'
          WRITE (out,60) min_ai(ng), 'min_ai',                          &
     &          'Minimum ice concentration (nondimensional) limiter.'
          WRITE (out,60) max_ai(ng), 'max_ai',                          &
     &          'Maximum ice concentration (nondimensional) limiter.'
          WRITE (out,60) min_hi(ng), 'min_hi',                          &
     &          'Minimum average ice thickness (m) limiter.'
          WRITE (out,60) max_hmelt(ng), 'max_hmelt',                    &
     &          'Maximum surface melt water thickness (m) limiter.'
          WRITE (out,60) stressAng(ng)*rad2deg, 'stressAng',            &
     &          'Turning angle for ice-water drag (degrees).'
          WRITE (out,60) ice_emiss, 'ice_emiss',                        &
     &          'Ice emissivity (nondimensional).'
          WRITE (out,60) spec_heat_air, 'spec_heat_air',                &
     &          'Specific heat of air (J/kg/K).'
          WRITE (out,60) trans_coeff, 'trans_coeff',                    &
     &          'Heat transfer coefficient (nondimensional).'
          WRITE (out,60) sublimation, 'sublimation',                    &
     &          'Latent heat of sublimation (J/kg).'
#ifdef T_PASSIVE
          DO i=1,NPT
            itrc=inert(i)
            IF (LtracerIce(itrc,ng)) THEN
              WRITE (out,185) LtracerIce(itrc,ng), 'LtracerIce', itrc,  &
     &            'Turning ON  ice-ocean tracer exchange on tracer ',   &
     &            itrc, TRIM(Vname(1,idStrc(itrc)))
            ELSE
              WRITE (out,185) LtracerIce(itrc,ng), 'LtracerIce', itrc,  &
     &            'Turning OFF ice-ocean tracer exchange on tracer ',   &
     &            itrc, TRIM(Vname(1,idStrc(itrc)))
            END IF
          END DO
#endif
!
          IF (nHIS(ng).gt.0) THEN
            WRITE (out,'(1x)')
            IF (Hout(idUice,ng)) WRITE (out,70) Hout(idUice,ng),        &
     &         'Hout(idUice)',                                          &
     &         'Write out ice U-velocity component.'
            IF (Hout(idVice,ng)) WRITE (out,70) Hout(idVice,ng),        &
     &         'Hout(idVice)',                                          &
     &         'Write out ice V-velocity component.'
            IF (Hout(idUiER,ng)) WRITE (out,70) Hout(idUiER,ng),        &
     &         'Hout(idUiER)',                                          &
     &         'Write out ice Eastward velocity at RHO-points.'
            IF (Hout(idViNR,ng)) WRITE (out,70) Hout(idViNR,ng),        &
     &         'Hout(idViNR)',                                          &
     &         'Write out ice Northward velocity at RHO-points.'
            IF (Hout(idAice,ng)) WRITE (out,70) Hout(idAice,ng),        &
     &         'Hout(idAice)',                                          &
     &         'Write out ice concentration (fractional area coverage).'
            IF (Hout(idIage,ng)) WRITE (out,70) Hout(idIage,ng),        &
     &         'Hout(idIage)',                                          &
     &         'Write out age of ice.'
            IF (Hout(idHice,ng)) WRITE (out,70) Hout(idHice,ng),        &
     &         'Hout(idHice)',                                          &
     &         'Write out average ice thickness (ice mass per area).'
            IF (Hout(idHmel,ng)) WRITE (out,70) Hout(idHmel,ng),        &
     &         'Hout(idHmel)',                                          &
     &         'Write out melt pond water thickness on ice.'
            IF (Hout(idHsno,ng)) WRITE (out,70) Hout(idHsno,ng),        &
     &         'Hout(idHsno)',                                          &
     &         'Write out average snow coverage thickness.'
            IF (Hout(idTice,ng)) WRITE (out,70) Hout(idTice,ng),        &
     &         'Hout(idTice)',                                          &
     &         'Write out interior ice temperature.'
            IF (Hout(idISxx,ng)) WRITE (out,70) Hout(idISxx,ng),        &
     &         'Hout(idISxx)',                                          &
     &         'Write out internal ice stress tensor, xx-component.'
            IF (Hout(idISxy,ng)) WRITE (out,70) Hout(idISxy,ng),        &
     &         'Hout(idISxy)',                                          &
     &         'Write out internal ice stress tensor, xy-component.'
            IF (Hout(idISyy,ng)) WRITE (out,70) Hout(idISyy,ng),        &
     &         'Hout(idISyy)',                                          &
     &         'Write out internal ice stress tensor, yy-component.'
#if defined ICE_BIO
            IF (Hout(idIphy,ng)) WRITE (out,70) Hout(idIphy,ng),        &
     &         'Hout(idIphy)',                                          &
     &         'Write out ice biology, algae.'
            IF (Hout(idINO3,ng)) WRITE (out,70) Hout(idINO3,ng),        &
     &         'Hout(idINO3)',                                          &
     &         'Write out ice biology, nitrate.'
            IF (Hout(idINH4,ng)) WRITE (out,70) Hout(idINH4,ng),        &
     &         'Hout(idINH4)',                                          &
     &         'Write out ice biology, ammonia.'
            IF (Hout(idIlog,ng)) WRITE (out,70) Hout(idIlog,ng),        &
     &         'Hout(idIlog)',                                          &
     &         'Write out ice biology log.'
#endif
            IF (Hout(idIsst,ng)) WRITE (out,70) Hout(idIsst,ng),        &
     &         'Hout(idIsst)',                                          &
     &         'Write out ice/snow surface temperature.'
            IF (Hout(idIOmf,ng)) WRITE (out,70) Hout(idIOmf,ng),        &
     &         'Hout(idIOmf)',                                          &
     &         'Write out ice-ocean mass flux.'
            IF (Hout(idIOfv,ng)) WRITE (out,70) Hout(idIOfv,ng),        &
     &         'Hout(idIOfv)',                                          &
     &         'Write out ice-ocean friction velocity.'
            IF (Hout(idIOmt,ng)) WRITE (out,70) Hout(idIOmt,ng),        &
     &         'Hout(idIOmt)',                                          &
     &         'Write out ice-ocean momentum transfer coefficient.'
            IF (Hout(idS0mk,ng)) WRITE (out,70) Hout(idS0mk,ng),        &
     &         'Hout(idS0mk)',                                          &
     &         'Write out salinity of molecular sublayer under ice.'
            IF (Hout(idT0mk,ng)) WRITE (out,70) Hout(idT0mk,ng),        &
     &         'Hout(idT0mk)',                                          &
     &         'Write out temperature of molecular sublayer under ice.'
            IF (Hout(idWdiv,ng)) WRITE (out,70) Hout(idWdiv,ng),        &
     &         'Hout(idWdiv)',                                          &
     &         'Write out rate of ice divergence.'
            IF (Hout(idW_fr,ng)) WRITE (out,70) Hout(idW_fr,ng),        &
     &         'Hout(idW_fr)',                                          &
     &         'Write out rate of ice accretion by frazil ice growth.'
            IF (Hout(idW_ai,ng)) WRITE (out,70) Hout(idW_ai,ng),        &
     &         'Hout(idW_ai)',                                          &
     &         'Write out rate of melt/freeze at air/ice interface.'
            IF (Hout(idW_ao,ng)) WRITE (out,70) Hout(idW_ao,ng),        &
     &         'Hout(idW_ao)',                                          &
     &         'Write out rate of melt/freeze at air/ocean interface.'
            IF (Hout(idW_io,ng)) WRITE (out,70) Hout(idW_io,ng),        &
     &         'Hout(idW_io)',                                          &
     &         'Write out rate of melt/freeze at ice/ocean interface.'
            IF (Hout(idW_ro,ng)) WRITE (out,70) Hout(idW_ro,ng),        &
     &         'Hout(idW_ro)',                                          &
     &         'Write out rate of melt/freeze runoff into ocean.'
#if defined ICE_BULK_FLUXES
            IF (Hout(idAlbi,ng)) WRITE (out,70) Hout(idAlbi,ng),        &
     &         'Hout(idAlbi)',                                          &
     &         'Write out sea ice albedo.'
            IF (Hout(idAlbo,ng)) WRITE (out,70) Hout(idAlbo,ng),        &
     &         'Hout(idAlbo)',                                          &
     &         'Write out sea surface albedo.'
            IF (Hout(idFrsw,ng)) WRITE (out,70) Hout(idFrsw,ng),        &
     &         'Hout(idFrsw)',                                          &
     &         'Write out short-wave penetration fraction in sea ice.'
            IF (Hout(idSics,ng)) WRITE (out,70) Hout(idSics,ng),        &
     &         'Hout(idSics)',                                          &
     &         'Write out net short-wave radiation at ice surafce.'
            IF (Hout(idSicb,ng)) WRITE (out,70) Hout(idSicb,ng),        &
     &         'Hout(idSicb)',                                          &
     &         'Write out net short-wave radiation at ice base.'
            IF (Hout(idSocn,ng)) WRITE (out,70) Hout(idSocn,ng),        &
     &         'Hout(idSocn)',                                          &
     &         'Write out net short-wave radiation at sea surface.'
            IF (Hout(idRfac,ng)) WRITE (out,70) Hout(idRfac,ng),        &
     &         'Hout(idRfac)',                                          &
     &         'Write out lmd_r1 scaling factor.'
            IF (Hout(idSnof,ng)) WRITE (out,70) Hout(idSnof,ng),        &
     &         'Hout(idSnof)',                                          &
     &         'Write out snowfall rate.'
#endif
#ifdef T_PASSIVE
            DO i=1,NPT
              itrc=inert(i)
              IF (Hout(idStrc(itrc),ng)) WRITE (out,180)                &
     &            Hout(idStrc(itrc),ng), 'Hout(idStrc)',                &
     &            'Write out sea ice passive tracer ', i,               &
     &            TRIM(Vname(1,idStrc(itrc)))
            END DO
#endif
          END IF
!
          IF (nQCK(ng).gt.0) THEN
            WRITE (out,'(1x)')
            IF (Qout(idUice,ng)) WRITE (out,70) Qout(idUice,ng),        &
     &         'Qout(idUice)',                                          &
     &         'Write out ice U-velocity component.'
            IF (Qout(idVice,ng)) WRITE (out,70) Qout(idVice,ng),        &
     &         'Qout(idVice)',                                          &
     &         'Write out ice V-velocity component.'
            IF (Qout(idUiER,ng)) WRITE (out,70) Qout(idUiER,ng),        &
     &         'Qout(idUiER)',                                          &
     &         'Write out ice Eastward velocity at RHO-points.'
            IF (Qout(idViNR,ng)) WRITE (out,70) Qout(idViNR,ng),        &
     &         'Qout(idViNR)',                                          &
     &         'Write out ice Northward velocity at RHO-points.'
            IF (Qout(idAice,ng)) WRITE (out,70) Qout(idAice,ng),        &
     &         'Qout(idAice)',                                          &
     &         'Write out ice concentration (fractional area coverage).'
            IF (Qout(idIage,ng)) WRITE (out,70) Qout(idIage,ng),        &
     &         'Qout(idIage)',                                          &
     &         'Write out age of ice.'
            IF (Qout(idHice,ng)) WRITE (out,70) Qout(idHice,ng),        &
     &         'Qout(idHice)',                                          &
     &         'Write out average ice thickness (ice mass per area).'
            IF (Qout(idHmel,ng)) WRITE (out,70) Qout(idHmel,ng),        &
     &         'Qout(idHmel)',                                          &
     &         'Write out surface melt water thickness on ice.'
            IF (Qout(idHsno,ng)) WRITE (out,70) Qout(idHsno,ng),        &
     &         'Qout(idHsno)',                                          &
     &         'Write out average snow coverage thickness.'
            IF (Qout(idTice,ng)) WRITE (out,70) Qout(idTice,ng),        &
     &         'Qout(idTice)',                                          &
     &         'Write out interior ice temperature.'
            IF (Qout(idISxx,ng)) WRITE (out,70) Qout(idISxx,ng),        &
     &         'Qout(idISxx)',                                          &
     &         'Write out internal ice stress tensor, xx-component.'
            IF (Qout(idISxy,ng)) WRITE (out,70) Qout(idISxy,ng),        &
     &         'Qout(idISxy)',                                          &
     &         'Write out internal ice stress tensor, xy-component.'
            IF (Qout(idISyy,ng)) WRITE (out,70) Qout(idISyy,ng),        &
     &         'Qout(idISyy)',                                          &
     &         'Write out internal ice stress tensor, yy-component.'
#if defined ICE_BIO
            IF (Qout(idIphy,ng)) WRITE (out,70) Qout(idIphy,ng),        &
     &         'Qout(idIphy)',                                          &
     &         'Write out ice biology, algae.'
            IF (Qout(idINO3,ng)) WRITE (out,70) Qout(idINO3,ng),        &
     &         'Qout(idINO3)',                                          &
     &         'Write out ice biology, nitrate.'
            IF (Qout(idINH4,ng)) WRITE (out,70) Qout(idINH4,ng),        &
     &         'Qout(idINH4)',                                          &
     &         'Write out ice biology, ammonia.'
            IF (Qout(idIlog,ng)) WRITE (out,70) Qout(idIlog,ng),        &
     &         'Qout(idIlog)',                                          &
     &         'Write out ice biology log.'
#endif
            IF (Qout(idIsst,ng)) WRITE (out,70) Qout(idIsst,ng),        &
     &         'Qout(idIsst)',                                          &
     &         'Write out ice/snow surface temperature.'
            IF (Qout(idIOmf,ng)) WRITE (out,70) Qout(idIOmf,ng),        &
     &         'Qout(idIOmf)',                                          &
     &         'Write out ice-ocean mass flux.'
            IF (Qout(idIOfv,ng)) WRITE (out,70) Qout(idIOfv,ng),        &
     &         'Qout(idIOfv)',                                          &
     &         'Write out ice-ocean friction velocity.'
            IF (Qout(idIOmt,ng)) WRITE (out,70) Qout(idIOmt,ng),        &
     &         'Qout(idIOmt)',                                          &
     &         'Write out ice-ocean momentum transfer coefficient.'
            IF (Qout(idS0mk,ng)) WRITE (out,70) Qout(idS0mk,ng),        &
     &         'Qout(idS0mk)',                                          &
     &         'Write out salinity of molecular sublayer under ice.'
            IF (Qout(idT0mk,ng)) WRITE (out,70) Qout(idT0mk,ng),        &
     &         'Qout(idT0mk)',                                          &
     &         'Write out temperature of molecular sublayer under ice.'
            IF (Qout(idWdiv,ng)) WRITE (out,70) Qout(idWdiv,ng),        &
     &         'Qout(idWdiv)',                                          &
     &         'Write out rate of ice divergence.'
            IF (Qout(idW_fr,ng)) WRITE (out,70) Qout(idW_fr,ng),        &
     &         'Qout(idW_fr)',                                          &
     &         'Write out rate of ice accretion by frazil ice growth.'
            IF (Qout(idW_ai,ng)) WRITE (out,70) Qout(idW_ai,ng),        &
     &         'Qout(idW_ai)',                                          &
     &         'Write out rate of melt/freeze at air/ice interface.'
            IF (Qout(idW_ao,ng)) WRITE (out,70) Qout(idW_ao,ng),        &
     &         'Qout(idW_ao)',                                          &
     &         'Write out rate of melt/freeze at air/ocean interface.'
            IF (Qout(idW_io,ng)) WRITE (out,70) Qout(idW_io,ng),        &
     &         'Qout(idW_io)',                                          &
     &         'Write out rate of melt/freeze at ice/ocean interface.'
            IF (Qout(idW_ro,ng)) WRITE (out,70) Qout(idW_ro,ng),        &
     &         'Qout(idW_ro)',                                          &
     &         'Write out rate of melt/freeze runoff into ocean.'
#if defined ICE_BULK_FLUXES
            IF (Qout(idAlbi,ng)) WRITE (out,70) Qout(idAlbi,ng),        &
     &         'Qout(idAlbi)',                                          &
     &         'Write out sea ice albedo.'
            IF (Qout(idAlbo,ng)) WRITE (out,70) Qout(idAlbo,ng),        &
     &         'Qout(idAlbo)',                                          &
     &         'Write out sea surface albedo.'
            IF (Qout(idFrsw,ng)) WRITE (out,70) Qout(idFrsw,ng),        &
     &         'Qout(idFrsw)',                                          &
     &         'Write out short-wave penetration fraction in sea ice.'
            IF (Qout(idSics,ng)) WRITE (out,70) Qout(idSics,ng),        &
     &         'Qout(idSics)',                                          &
     &         'Write out net short-wave radiation at ice surafce.'
            IF (Qout(idSicb,ng)) WRITE (out,70) Qout(idSicb,ng),        &
     &         'Qout(idSicb)',                                          &
     &         'Write out net short-wave radiation at ice base.'
            IF (Qout(idSocn,ng)) WRITE (out,70) Qout(idSocn,ng),        &
     &         'Qout(idSocn)',                                          &
     &         'Write out net short-wave radiation at sea surface.'
            IF (Qout(idRfac,ng)) WRITE (out,70) Qout(idRfac,ng),        &
     &         'Qout(idRfac)',                                          &
     &         'Write out lmd_r1 scaling factor.'
            IF (Qout(idSnof,ng)) WRITE (out,70) Qout(idSnof,ng),        &
     &         'Qout(idSnof)',                                          &
     &         'Write out snowfall rate.'
#endif
#ifdef T_PASSIVE
            DO i=1,NPT
              itrc=inert(i)
              IF (Qout(idStrc(itrc),ng)) WRITE (out,180)                &
     &            Qout(idStrc(itrc),ng), 'Qout(idStrc)',                &
     &            'Write out sea ice passive tracer ', i,               &
     &            TRIM(Vname(1,idStrc(itrc)))
            END DO
#endif
          END IF
#ifdef AVERAGES
!
          IF (nAVG(ng).gt.0) THEN
            WRITE (out,'(1x)')
            IF (Aout(idUice,ng)) WRITE (out,70) Aout(idUice,ng),        &
     &         'Aout(idUice)',                                          &
     &         'Write out ice U-velocity component.'
            IF (Aout(idVice,ng)) WRITE (out,70) Aout(idVice,ng),        &
     &         'Aout(idVice)',                                          &
     &         'Write out ice V-velocity component.'
            IF (Aout(idUiER,ng)) WRITE (out,70) Aout(idUiER,ng),        &
     &         'Aout(idUiER)',                                          &
     &         'Write out ice Eastward velocity at RHO-points.'
            IF (Aout(idViNR,ng)) WRITE (out,70) Aout(idViNR,ng),        &
     &         'Aout(idViNR)',                                          &
     &         'Write out ice Northward velocity at RHO-points.'
            IF (Aout(idAice,ng)) WRITE (out,70) Aout(idAice,ng),        &
     &         'Aout(idAice)',                                          &
     &         'Write out ice concentration (fractional area coverage).'
            IF (Aout(idIage,ng)) WRITE (out,70) Aout(idIage,ng),        &
     &         'Aout(idIage)',                                          &
     &         'Write out age of ice.'
            IF (Aout(idHice,ng)) WRITE (out,70) Aout(idHice,ng),        &
     &         'Aout(idHice)',                                          &
     &         'Write out average ice thickness (ice mass per area).'
            IF (Aout(idHmel,ng)) WRITE (out,70) Aout(idHmel,ng),        &
     &         'Aout(idHmel)',                                          &
     &         'Write out surface melt water thickness on ice.'
            IF (Aout(idHsno,ng)) WRITE (out,70) Aout(idHsno,ng),        &
     &         'Aout(idHsno)',                                          &
     &         'Write out average snow coverage thickness.'
            IF (Aout(idTice,ng)) WRITE (out,70) Aout(idTice,ng),        &
     &         'Aout(idTice)',                                          &
     &         'Write out interior ice temperature.'
            IF (Aout(idISxx,ng)) WRITE (out,70) Aout(idISxx,ng),        &
     &         'Aout(idISxx)',                                          &
     &         'Write out internal ice stress tensor, xx-component.'
            IF (Aout(idISxy,ng)) WRITE (out,70) Aout(idISxy,ng),        &
     &         'Aout(idISxy)',                                          &
     &         'Write out internal ice stress tensor, xy-component.'
            IF (Aout(idISyy,ng)) WRITE (out,70) Aout(idISyy,ng),        &
     &         'Aout(idISyy)',                                          &
     &         'Write out internal ice stress tensor, yy-component.'
# if defined ICE_BIO
            IF (Aout(idIphy,ng)) WRITE (out,70) Aout(idIphy,ng),        &
     &         'Aout(idIphy)',                                          &
     &         'Write out ice biology, algae.'
            IF (Aout(idINO3,ng)) WRITE (out,70) Aout(idINO3,ng),        &
     &         'Aout(idINO3)',                                          &
     &         'Write out ice biology, nitrate.'
            IF (Aout(idINH4,ng)) WRITE (out,70) Aout(idINH4,ng),        &
     &         'Aout(idINH4)',                                          &
     &         'Write out ice biology, ammonia.'
            IF (Aout(idIlog,ng)) WRITE (out,70) Aout(idIlog,ng),        &
     &         'Aout(idIlog)',                                          &
     &         'Write out ice biology log.'
# endif
            IF (Aout(idIsst,ng)) WRITE (out,70) Aout(idisst,ng),        &
     &         'Aout(idIsst)',                                          &
     &         'Write out ice/snow surface temperature.'
            IF (Aout(idIOmf,ng)) WRITE (out,70) Aout(idIOmf,ng),        &
     &         'Aout(idIOmf)',                                          &
     &         'Write out ice-ocean mass flux.'
            IF (Aout(idIOfv,ng)) WRITE (out,70) Aout(idIOfv,ng),        &
     &         'Aout(idIOfv)',                                          &
     &         'Write out ice-ocean friction velocity.'
            IF (Aout(idIOmt,ng)) WRITE (out,70) Aout(idIOmt,ng),        &
     &         'Aout(idIOmt)',                                          &
     &         'Write out ice-ocean momentum transfer coefficient.'
            IF (Aout(idS0mk,ng)) WRITE (out,70) Aout(idS0mk,ng),        &
     &         'Aout(idS0mk)',                                          &
     &         'Write out salinity of molecular sublayer under ice.'
            IF (Aout(idT0mk,ng)) WRITE (out,70) Aout(idT0mk,ng),        &
     &         'Aout(idT0mk)',                                          &
     &         'Write out temperature of molecular sublayer under ice.'
            IF (Aout(idWdiv,ng)) WRITE (out,70) Aout(idWdiv,ng),        &
     &         'Aout(idWdiv)',                                          &
     &         'Write out rate of ice divergence.'
            IF (Aout(idW_fr,ng)) WRITE (out,70) Aout(idW_fr,ng),        &
     &         'Aout(idW_fr)',                                          &
     &         'Write out rate of ice accretion by frazil ice growth.'
            IF (Aout(idW_ai,ng)) WRITE (out,70) Aout(idW_ai,ng),        &
     &         'Aout(idW_ai)',                                          &
     &         'Write out rate of melt/freeze at air/ice interface.'
            IF (Aout(idW_ao,ng)) WRITE (out,70) Aout(idW_ao,ng),        &
     &         'Aout(idW_ao)',                                          &
     &         'Write out rate of melt/freeze at air/ocean interface.'
            IF (Aout(idW_io,ng)) WRITE (out,70) Aout(idW_io,ng),        &
     &         'Aout(idW_io)',                                          &
     &         'Write out rate of melt/freeze at ice/ocean interface.'
            IF (Aout(idW_ro,ng)) WRITE (out,70) Aout(idW_ro,ng),        &
     &         'Aout(idW_ro)',                                          &
     &         'Write out rate of melt/freeze runoff into ocean.'
# if defined ICE_BULK_FLUXES
            IF (Aout(idAlbi,ng)) WRITE (out,70) Aout(idAlbi,ng),        &
     &         'Aout(idAlbi)',                                          &
     &         'Write out sea ice albedo.'
            IF (Aout(idAlbo,ng)) WRITE (out,70) Aout(idAlbo,ng),        &
     &         'Aout(idAlbo)',                                          &
     &         'Write out sea surface albedo.'
            IF (Aout(idFrsw,ng)) WRITE (out,70) Aout(idFrsw,ng),        &
     &         'Aout(idFrsw)',                                          &
     &         'Write out short-wave penetration fraction in sea ice.'
            IF (Aout(idSics,ng)) WRITE (out,70) Aout(idSics,ng),        &
     &         'Aout(idSics)',                                          &
     &         'Write out net short-wave radiation at ice surafce.'
            IF (Aout(idSicb,ng)) WRITE (out,70) Aout(idSicb,ng),        &
     &         'Aout(idSicb)',                                          &
     &         'Write out net short-wave radiation at ice base.'
            IF (Aout(idSocn,ng)) WRITE (out,70) Aout(idSocn,ng),        &
     &         'Aout(idSocn)',                                          &
     &         'Write out net short-wave radiation at sea surface.'
            IF (Aout(idRfac,ng)) WRITE (out,70) Aout(idRfac,ng),        &
     &         'Aout(idRfac)',                                          &
     &         'Write out lmd_r1 scaling factor.'
            IF (Aout(idSnof,ng)) WRITE (out,70) Aout(idSnof,ng),        &
     &         'Aout(idSnof)',                                          &
     &         'Write out snowfall rate.'
# endif
# ifdef T_PASSIVE
            DO i=1,NPT
              itrc=inert(i)
              IF (Aout(idStrc(itrc),ng)) WRITE (out,180)                &
     &            Aout(idStrc(itrc),ng), 'Aout(idStrc)',                &
     &            'Write out sea ice passive tracer ', i,               &
     &            TRIM(Vname(1,idStrc(itrc)))
            END DO
# endif
          END IF
#endif
#ifdef STATIONS
          IF (Lstations(ng).and.nSTA(ng).gt.0) THEN
            WRITE (out,'(1x)')
            IF (Sout(idUice,ng)) WRITE (out,70) Sout(idUice,ng),        &
     &         'Sout(idUice)',                                          &
     &         'Write out ice U-velocity component.'
            IF (Sout(idVice,ng)) WRITE (out,70) Sout(idVice,ng),        &
     &         'Sout(idVice)',                                          &
     &         'Write out ice V-velocity component.'
            IF (Sout(idUiER,ng)) WRITE (out,70) Sout(idUiER,ng),        &
     &         'Sout(idUiER)',                                          &
     &         'Write out ice Eastward velocity at RHO-points.'
            IF (Sout(idViNR,ng)) WRITE (out,70) Sout(idViNR,ng),        &
     &         'Sout(idViNR)',                                          &
     &         'Write out ice Northward velocity at RHO-points.'
            IF (Sout(idAice,ng)) WRITE (out,70) Sout(idAice,ng),        &
     &         'Sout(idAice)',                                          &
     &         'Write out ice concentration (fractional area coverage).'
            IF (Sout(idIage,ng)) WRITE (out,70) Sout(idIage,ng),        &
     &         'Sout(idIage)',                                          &
     &         'Write out age of ice.'
            IF (Sout(idHice,ng)) WRITE (out,70) Sout(idHice,ng),        &
     &         'Sout(idHice)',                                          &
     &         'Write out average ice thickness (ice mass per area).'
            IF (Sout(idHmel,ng)) WRITE (out,70) Sout(idHmel,ng),        &
     &         'Sout(idHmel)',                                          &
     &         'Write out surface melt water thickness on ice.'
            IF (Sout(idHsno,ng)) WRITE (out,70) Sout(idHsno,ng),        &
     &         'Sout(idHsno)',                                          &
     &         'Write out average snow coverage thickness.'
            IF (Sout(idTice,ng)) WRITE (out,70) Sout(idTice,ng),        &
     &         'Sout(idTice)',                                          &
     &         'Write out interior ice temperature.'
            IF (Sout(idISxx,ng)) WRITE (out,70) Sout(idISxx,ng),        &
     &         'Sout(idISxx)',                                          &
     &         'Write out internal ice stress tensor, xx-component.'
            IF (Sout(idISxy,ng)) WRITE (out,70) Sout(idISxy,ng),        &
     &         'Sout(idISxy)',                                          &
     &         'Write out internal ice stress tensor, xy-component.'
            IF (Sout(idISyy,ng)) WRITE (out,70) Sout(idISyy,ng),        &
     &         'Sout(idISyy)',                                          &
     &         'Write out internal ice stress tensor, yy-component.'
# if defined ICE_BIO
            IF (Sout(idIphy,ng)) WRITE (out,70) Sout(idIphy,ng),        &
     &         'Sout(idIphy)',                                          &
     &         'Write out ice biology, algae.'
            IF (Sout(idINO3,ng)) WRITE (out,70) Sout(idINO3,ng),        &
     &         'Sout(idINO3)',                                          &
     &         'Write out ice biology, nitrate.'
            IF (Sout(idINH4,ng)) WRITE (out,70) Sout(idINH4,ng),        &
     &         'Sout(idINH4)',                                          &
     &         'Write out ice biology, ammonia.'
            IF (Sout(idIlog,ng)) WRITE (out,70) Sout(idIlog,ng),        &
     &         'Sout(idIlog)',                                          &
     &         'Write out ice biology log.'
# endif
            IF (Sout(idIsst,ng)) WRITE (out,70) Sout(idIsst,ng),        &
     &         'Sout(idIsst)',                                          &
     &         'Write out ice/snow surface temperature.'
            IF (Sout(idIOmf,ng)) WRITE (out,70) Sout(idIOmf,ng),        &
     &         'Sout(idIOmf)',                                          &
     &         'Write out ice-ocean mass flux.'
            IF (Sout(idIOfv,ng)) WRITE (out,70) Sout(idIOfv,ng),        &
     &         'Sout(idIOfv)',                                          &
     &         'Write out ice-ocean friction velocity.'
            IF (Sout(idIOmt,ng)) WRITE (out,70) Sout(idIOmt,ng),        &
     &         'Sout(idIOmt)',                                          &
     &         'Write out ice-ocean momentum transfer coefficient.'
            IF (Sout(idS0mk,ng)) WRITE (out,70) Sout(idS0mk,ng),        &
     &         'Sout(idS0mk)',                                          &
     &         'Write out salinity of molecular sublayer under ice.'
            IF (Sout(idT0mk,ng)) WRITE (out,70) Sout(idT0mk,ng),        &
     &         'Sout(idT0mk)',                                          &
     &         'Write out temperature of molecular sublayer under ice.'
            IF (Sout(idWdiv,ng)) WRITE (out,70) Sout(idWdiv,ng),        &
     &         'Sout(idWdiv)',                                          &
     &         'Write out rate of ice divergence.'
            IF (Sout(idW_fr,ng)) WRITE (out,70) Sout(idW_fr,ng),        &
     &         'Sout(idW_fr)',                                          &
     &         'Write out rate of ice accretion by frazil ice growth.'
            IF (Sout(idW_ai,ng)) WRITE (out,70) Sout(idW_ai,ng),        &
     &         'Sout(idW_ai)',                                          &
     &         'Write out rate of melt/freeze at air/ice interface.'
            IF (Sout(idW_ao,ng)) WRITE (out,70) Sout(idW_ao,ng),        &
     &         'Sout(idW_ao)',                                          &
     &         'Write out rate of melt/freeze at air/ocean interface.'
            IF (Sout(idW_io,ng)) WRITE (out,70) Sout(idW_io,ng),        &
     &         'Sout(idW_io)',                                          &
     &         'Write out rate of melt/freeze at ice/ocean interface.'
            IF (Sout(idW_ro,ng)) WRITE (out,70) Sout(idW_ro,ng),        &
     &         'Sout(idW_ro)',                                          &
     &         'Write out rate of melt/freeze runoff into ocean.'
# if defined ICE_BULK_FLUXES
            IF (Sout(idAlbi,ng)) WRITE (out,70) Sout(idAlbi,ng),        &
     &         'Sout(idAlbi)',                                          &
     &         'Write out sea ice albedo.'
            IF (Sout(idAlbo,ng)) WRITE (out,70) Sout(idAlbo,ng),        &
     &         'Sout(idAlbo)',                                          &
     &         'Write out ocean surface albedo.'
            IF (Sout(idFrsw,ng)) WRITE (out,70) Sout(idFrsw,ng),        &
     &         'Sout(idFrsw)',                                          &
     &         'Write out short-wave penetration fraction in sea ice.'
            IF (Sout(idSics,ng)) WRITE (out,70) Sout(idSics,ng),        &
     &         'Sout(idSics)',                                          &
     &         'Write out net short-wave flux at sea ice surface.'
            IF (Sout(idSicb,ng)) WRITE (out,70) Sout(idSicb,ng),        &
     &         'Sout(idSicb)',                                          &
     &         'Write out net short-wave flux at sea ice base.'
            IF (Sout(idSocn,ng)) WRITE (out,70) Sout(idSocn,ng),        &
     &         'Sout(idSocn)',                                          &
     &         'Write out net short-wave flux at sea surface.'
            IF (Sout(idRfac,ng)) WRITE (out,70) Sout(idRfac,ng),        &
     &         'Sout(idRfac)',                                          &
     &         'Write out lmd_r1 scaling factor.'
            IF (Sout(idSnof,ng)) WRITE (out,70) Sout(idSnof,ng),        &
     &         'Sout(idSnof)',                                          &
     &         'Write out snowfall rate.'
# endif
# ifdef T_PASSIVE
            DO i=1,NPT
              itrc=inert(i)
              IF (Sout(idStrc(itrc),ng)) WRITE (out,180)                &
     &            Sout(idStrc(itrc),ng), 'Sout(idStrc)',                &
     &            'Write out sea ice passive tracer ', i,               &
     &            TRIM(Vname(1,idStrc(itrc)))
            END DO
# endif
          END IF
#endif
        END DO
      END IF

#ifdef AVERAGES
!
!  Turn out switches to process time-averaged ice model state and
!  internal variables. It is necessary to avoid recursive dependency
!  between "mod_ncparam" and "mod_ice".
!
      LiceFavg=.FALSE.
      LiceSavg=.FALSE.
!
      DO ng=1,Ngrids
        DO nf=1,nIceF
          IF (iFice(nf).gt.0) THEN
            IF (Aout(iFice(nf),ng)) THEN
              LiceFavg(nf,ng)=.TRUE.
            END IF
          END IF
        END DO
!
        DO ns=1,nIceS
          IF (iSice(ns).gt.0) THEN
            IF (Aout(iSice(ns),ng)) THEN
              LiceSavg(ns,ng)=.TRUE.
            END IF
          END IF
        END DO
      END DO
# ifdef T_PASSIVE
!
      DO ng=1,Ngrids
        DO i=1,NPT
          itrc=inert(i)
          IF (Aout(idStrc(itrc),ng)) THEN
            LiceTavg(itrc,ng)=.TRUE.
          END IF
        END DO
      END DO
# endif
#endif
!
  30  FORMAT (/,' READ_IcePar - Error while processing line: ',/,a)
  40  FORMAT (/,/,' Ice Parameters, Grid: ',i2.2,                       &
     &        /,  ' ========================',/)
  50  FORMAT (1x,i10,2x,a,t32,a)
  60  FORMAT (1p,e11.4,2x,a,t32,a)
  70  FORMAT (10x,l1,2x,a,t32,a)
  80  FORMAT (/,' READ_IcePar - variable index not yet loaded, ', a)
 180  FORMAT (10x,l1,2x,a,t32,a,i2.2,':',1x,a)
 185  FORMAT (10x,l1,2x,a,'(',i2.2,')',t32,a,i2.2,':',1x,a)
!
      RETURN
      END SUBROUTINE read_IcePar
