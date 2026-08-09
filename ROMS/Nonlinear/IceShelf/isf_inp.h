      SUBROUTINE read_IsfPar (model, inp, out, Lwrite)
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2024 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This routine reads and reports iceshelf model input parameters.     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_isf
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
      integer :: soff, igrid, nline
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
!  Allocate iceshelf model parameters that depend on Ngrids.
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
!  Read in iceshelf model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.true.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
#ifdef ICESHEET
            CASE ('LBC(isIhis)')
              soff = isf_ic_brylast + isIhis
              Npts=load_lbc(Nval, Cval, line, nline, soff, igrid,       &
     &                      0, 0, Vname(1,idIhis), LBC)
            CASE ('LBC(isIuis)')
              soff = isf_ic_brylast + isIuis
              Npts=load_lbc(Nval, Cval, line, nline, soff, igrid,       &
     &                      0, 0, Vname(1,idIuis), LBC)
            CASE ('LBC(isIvis)')
              soff = isf_ic_brylast + isIvis
              Npts=load_lbc(Nval, Cval, line, nline, soff, igrid,       &
     &                      0, 0, Vname(1,idIvis), LBC)
            CASE ('LBC(isIeis)')
              soff = isf_ic_brylast + isIeis
              Npts=load_lbc(Nval, Cval, line, nline, soff, igrid,       &
     &                      0, 0, Vname(1,idIeis), LBC)
#endif
            CASE ('isf_drdz')
              Npts=load_r(Nval, Rval, Ngrids, isf_drdz)
#ifdef ICESHELF
# ifdef ICESHELF_LDRAG
            CASE ('RDRG_ISF')
              Npts=load_r(Nval, Rval, Ngrids, rdrg_isf)
# endif
# ifdef ICESHELF_QDRAG
            CASE ('RDRG2_ISF')
              Npts=load_r(Nval, Rval, Ngrids, rdrg2_isf)
# endif
# ifdef ICESHELF_LOGDRAG
            CASE ('ZOSUR_ISF')
              Npts=load_r(Nval, Rval, Ngrids, ZoSr_isf)
# endif
#endif
#ifdef ICESHELF_THERMO
            CASE ('gamU')
              Npts=load_r(Nval, Rval, Ngrids, gamU)
            CASE ('min_ustar')
              Npts=load_r(Nval, Rval, Ngrids, min_ustar)
            CASE ('isf_Hsbl')
              Npts=load_r(Nval, Rval, Ngrids, isf_Hsbl)
            CASE ('tau_melt')
              Npts=load_r(Nval, Rval, Ngrids, tau_melt)
# ifdef T_PASSIVE
            CASE ('LtracerIsf')
              Npts=load_l(Nval, Cval, NPT, Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=inert(i)
                  LtracerIsf(itrc,ng)=Linert(i,ng)
                END DO
              END DO
              DO ng=1,Ngrids
                IF (ANY(LtracerIsf(:,ng))) THEN
                  Lisftrc(ng)=.TRUE.
                  ic=0
                  DO itrc=NAT+1,NT(ng)
                    IF (LtracerIsf(itrc,ng)) THEN
                      ic=ic+1
                    END IF
                  END DO
                  NTISF(ng)=ic
                END IF
              END DO
            CASE ('LsourceIsf')
              Npts=load_l(Nval, Cval, NPT, Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=inert(i)
                  LsourceIsf(itrc,ng)=Linert(i,ng)
                END DO
              END DO
              DO ng=1,Ngrids
                IF (ANY(LsourceIsf(:,ng))) THEN
                  Lisfsrc(ng)=.TRUE.
                END IF
              END DO
            CASE ('trcIsf')
              Npts=load_r(Nval, Rval, NPT, Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=inert(i)
                  trcIsf(itrc,ng)=Rtracer(i,ng)
                END DO
              END DO
# endif
#endif
#ifdef ICESHEET
!
!  SSA physical parameters
!
            CASE ('n_glen')
              Npts=load_r(Nval, Rval, Ngrids, n_glen)
            CASE ('A_glen')
              Npts=load_r(Nval, Rval, Ngrids, A_glen)
            CASE ('eps_glen_min')
              Npts=load_r(Nval, Rval, Ngrids, eps_glen_min)
            CASE ('C_basal')
              Npts=load_r(Nval, Rval, Ngrids, C_basal)
            CASE ('n_basal')
              Npts=load_r(Nval, Rval, Ngrids, n_basal)
            CASE ('eps_u_min')
              Npts=load_r(Nval, Rval, Ngrids, eps_u_min)
!
!  SSA solver parameters
!
            CASE ('nl_max_iter')
              Npts=load_i(Nval, Rval, Ngrids, nl_max_iter)
            CASE ('nl_tol')
              Npts=load_r(Nval, Rval, Ngrids, nl_tol)
            CASE ('nl_tol_fp')
              Npts=load_r(Nval, Rval, Ngrids, nl_tol_fp)
            CASE ('cg_max_iter')
              Npts=load_i(Nval, Rval, Ngrids, cg_max_iter)
            CASE ('cg_tol')
              Npts=load_r(Nval, Rval, Ngrids, cg_tol)
            CASE ('nuH_reg')
              Npts=load_r(Nval, Rval, Ngrids, nuH_reg)
            CASE ('min_thickness')
              Npts=load_r(Nval, Rval, Ngrids, min_thickness)
            CASE ('ssa_max_speed')
              Npts=load_r(Nval, Rval, Ngrids, ssa_max_speed)
!
!  Driving stress regularization parameters
!
            CASE ('n_smooth')
              Npts=load_i(Nval, Rval, Ngrids, n_smooth)
            CASE ('n_taper')
              Npts=load_i(Nval, Rval, Ngrids, n_taper)
!
!  Picard iteration parameters
!
            CASE ('picard_alpha')
              Npts=load_r(Nval, Rval, Ngrids, picard_alpha)
!
!  Thickness advection parameters
!
            CASE ('cfl_factor')
              Npts=load_r(Nval, Rval, Ngrids, cfl_factor)
!
!  Calving parameters
!
            CASE ('H_calve')
              Npts=load_r(Nval, Rval, Ngrids, H_calve)
!
!  Multi-timescale parameters
!
            CASE ('N_ssa')
              Npts=load_i(Nval, Rval, Ngrids, N_ssa)
            CASE ('N_thick')
              Npts=load_i(Nval, Rval, Ngrids, N_thick)
            CASE ('N_gl')
              Npts=load_i(Nval, Rval, Ngrids, N_gl)
!
!  Surface mass balance parameters
!
            CASE ('smb_rate')
              Npts=load_r(Nval, Rval, Ngrids, smb_rate)
!
!  SSA geometry parameters (available under ICESHEET)
!
            CASE ('max_thickness')
              Npts=load_r(Nval, Rval, Ngrids, max_thickness)
            CASE ('gl_float_width')
              Npts=load_r(Nval, Rval, Ngrids, gl_float_width)
#endif
#ifdef ICESHEET
            CASE ('Hout(idIhis)')
              IF (idIhis.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIhis'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idIhis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idIuis)')
              IF (idIuis.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIuis'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idIuis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idIvis)')
              IF (idIvis.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIvis'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idIvis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idIuer)')
              IF (idIuer.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIuer'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idIuer,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idIvnr)')
              IF (idIvnr.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIvnr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idIvnr,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idItis)')
              IF (idItis.eq.0) THEN
                IF (Master) WRITE (out,80) 'idItis'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idItis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idIeis)')
              IF (idIeis.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIeis'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idIeis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idIzsf)')
              IF (idIzsf.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIzsf'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idIzsf,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idIzbd)')
              IF (idIzbd.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIzbd'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idIzbd,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idInuh)')
              IF (idInuh.eq.0) THEN
                IF (Master) WRITE (out,80) 'idInuh'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idInuh,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idItdx)')
              IF (idItdx.eq.0) THEN
                IF (Master) WRITE (out,80) 'idItdx'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idItdx,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idItdy)')
              IF (idItdy.eq.0) THEN
                IF (Master) WRITE (out,80) 'idItdy'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idItdy,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idIbta)')
              IF (idIbta.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIbta'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idIbta,1:Ngrids)=Lswitch(1:Ngrids)
#endif
#ifdef ICESHELF_THERMO
            CASE ('Hout(idImel)')
              IF (idImel.eq.0) THEN
                IF (Master) WRITE (out,80) 'idImel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idImel,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idIgat)')
              IF (idIgat.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIgat'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idIgat,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idIgas)')
              IF (idIgas.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIgas'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idIgas,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idItst)')
              IF (idItst.eq.0) THEN
                IF (Master) WRITE (out,80) 'idItst'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idItst,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idIust)')
              IF (idIust.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIust'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idIust,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idItbo)')
              IF (idItbo.eq.0) THEN
                IF (Master) WRITE (out,80) 'idItbo'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idItbo,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idIsbo)')
              IF (idIsbo.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIsbo'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idIsbo,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idItba)')
              IF (idItba.eq.0) THEN
                IF (Master) WRITE (out,80) 'idItba'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idItba,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idIsba)')
              IF (idIsba.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIsba'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idIsba,1:Ngrids)=Lswitch(1:Ngrids)
# ifdef T_PASSIVE
            CASE ('Hout(idIvol)')
              IF (idIvol.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIvol'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              DO ng=1,Ngrids
                IF (Lisftrc(ng)) THEN
                  Hout(idIvol,ng)=Lswitch(ng)
                END IF
              END DO
            CASE ('Hout(idItrc)')
              IF (MAXVAL(idItrc).eq.0) THEN
                IF (Master) WRITE (out,80) 'idItrc'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NPT, Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=inert(i)
                  IF (LtracerIsf(itrc,ng)) THEN
                    Hout(idItrc(itrc),ng)=Linert(i,ng)
                  END IF
                END DO
              END DO
            CASE ('Hout(idIamb)')
              IF (MAXVAL(idIamb).eq.0) THEN
                IF (Master) WRITE (out,80) 'idIamb'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NPT, Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=inert(i)
                  IF (LtracerIsf(itrc,ng)) THEN
                    Hout(idIamb(itrc),ng)=Linert(i,ng)
                  END IF
                END DO
              END DO
# endif
#endif
#ifdef ICESHEET
            CASE ('Qout(idIhis)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idIhis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idIuis)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idIuis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idIvis)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idIvis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idIuer)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idIuer,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idIvnr)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idIvnr,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idItis)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idItis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idIeis)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idIeis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idIzsf)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idIzsf,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idIzbd)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idIzbd,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idInuh)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idInuh,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idItdx)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idItdx,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idItdy)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idItdy,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idIbta)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idIbta,1:Ngrids)=Lswitch(1:Ngrids)
#endif
#ifdef ICESHELF_THERMO
            CASE ('Qout(idImel)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idImel,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idIgat)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idIgat,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idIgas)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idIgas,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idItst)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idItst,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idIust)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idIust,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idItbo)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idItbo,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idIsbo)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idIsbo,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idItba)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idItba,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idIsba)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idIsba,1:Ngrids)=Lswitch(1:Ngrids)
# ifdef T_PASSIVE
            CASE ('Qout(idIvol)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              DO ng=1,Ngrids
                IF (Lisftrc(ng)) THEN
                  Qout(idIvol,ng)=Lswitch(ng)
                END IF
              END DO
            CASE ('Qout(idItrc)')
              Npts=load_l(Nval, Cval, NPT, Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=inert(i)
                  IF (LtracerIsf(itrc,ng)) THEN
                    Qout(idItrc(itrc),ng)=Linert(i,ng)
                  END IF
                END DO
              END DO
            CASE ('Qout(idIamb)')
              Npts=load_l(Nval, Cval, NPT, Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=inert(i)
                  IF (LtracerIsf(itrc,ng)) THEN
                    Qout(idIamb(itrc),ng)=Linert(i,ng)
                  END IF
                END DO
              END DO
# endif
#endif
#ifdef AVERAGES
# ifdef ICESHEET
            CASE ('Aout(idIhis)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idIhis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idIuis)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idIuis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idIvis)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idIvis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idIuer)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idIuer,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idIvnr)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idIvnr,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idItis)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idItis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idIeis)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idIeis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idIzsf)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idIzsf,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idIzbd)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idIzbd,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idInuh)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idInuh,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idItdx)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idItdx,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idItdy)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idItdy,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idIbta)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idIbta,1:Ngrids)=Lswitch(1:Ngrids)
# endif
# ifdef ICESHELF_THERMO
            CASE ('Aout(idImel)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idImel,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idIgat)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idIgat,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idIgas)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idIgas,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idItst)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idItst,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idIust)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idIust,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idItbo)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idItbo,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idIsbo)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idIsbo,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idItba)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idItba,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idIsba)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idIsba,1:Ngrids)=Lswitch(1:Ngrids)
#  ifdef T_PASSIVE
            CASE ('Aout(idIvol)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              DO ng=1,Ngrids
                IF (Lisftrc(ng)) THEN
                  Aout(idIvol,ng)=Lswitch(ng)
                END IF
              END DO
            CASE ('Aout(idItrc)')
              Npts=load_l(Nval, Cval, NPT, Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=inert(i)
                  IF (LtracerIsf(itrc,ng)) THEN
                    Aout(idItrc(itrc),ng)=Linert(i,ng)
                  END IF
                END DO
              END DO
            CASE ('Aout(idIamb)')
              Npts=load_l(Nval, Cval, NPT, Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=inert(i)
                  IF (LtracerIsf(itrc,ng)) THEN
                    Aout(idIamb(itrc),ng)=Linert(i,ng)
                  END IF
                END DO
              END DO
#  endif
# endif
#endif
#ifdef STATIONS
# ifdef ICESHEET
            CASE ('Sout(idIhis)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idIhis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idIuis)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idIuis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idIvis)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idIvis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idIuer)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idIuer,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idIvnr)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idIvnr,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idItis)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idItis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idIeis)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idIeis,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idIzsf)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idIzsf,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idIzbd)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idIzbd,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idInuh)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idInuh,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idItdx)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idItdx,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idItdy)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idItdy,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idIbta)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idIbta,1:Ngrids)=Lswitch(1:Ngrids)
# endif
# ifdef ICESHELF_THERMO
            CASE ('Sout(idImel)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idImel,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idIgat)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idIgat,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idIgas)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idIgas,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idItst)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idItst,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idIust)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idIust,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idItbo)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idItbo,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idIsbo)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idIsbo,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idItba)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idItba,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Sout(idIsba)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Sout(idIsba,1:Ngrids)=Lswitch(1:Ngrids)
#  ifdef T_PASSIVE
            CASE ('Sout(idIvol)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              DO ng=1,Ngrids
                IF (Lisftrc(ng)) THEN
                  Sout(idIvol,ng)=Lswitch(ng)
                END IF
              END DO
            CASE ('Sout(idItrc)')
              Npts=load_l(Nval, Cval, NPT, Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=inert(i)
                  IF (LtracerIsf(itrc,ng)) THEN
                    Sout(idItrc(itrc),ng)=Linert(i,ng)
                  END IF
                END DO
              END DO
            CASE ('Sout(idIamb)')
              Npts=load_l(Nval, Cval, NPT, Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=inert(i)
                  IF (LtracerIsf(itrc,ng)) THEN
                    Sout(idIamb(itrc),ng)=Linert(i,ng)
                  END IF
                END DO
              END DO
#  endif
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

!     DO ng = 1,Ngrids
!       dtice(ng) = dt(ng)
!     END DO
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          WRITE (out,40) ng
#ifdef ICESHELF
# ifdef ICESHELF_LDRAG
          WRITE (out,60) rdrg_isf(ng), 'rdrg_isf',                      &
     &          'Linear ice shelf drag coefficient (m/s).'
# endif
# ifdef ICESHELF_QDRAG
          WRITE (out,60) rdrg2_isf(ng), 'rdrg2_isf',                    &
     &          'Quadratic ice shelf drag coefficient.'
# endif
# ifdef ICESHELF_LOGDRAG
          WRITE (out,60) ZoSr_isf(ng), 'ZoSr_isf',                      &
     &          'Ice shelf surface roughness (m).'
# endif
#endif
#ifdef ICESHELF_THERMO
          WRITE (out,60) gamU(ng), 'gamU',                              &
     &          'Iceshelf Ustar scaling factor (nondimensional).'
          WRITE (out,60) min_ustar(ng), 'min_ustar',                    &
     &          'Iceshelf minimum Ustar value.'
          WRITE (out,60) isf_Hsbl(ng), 'isf_Hsbl',                      &
     &          'Iceshelf sub-boundary layer thickness [m].'
          WRITE (out,60) tau_melt(ng), 'tau_melt',                      &
     &          'Melt rate EMA filter timescale [s].'
# ifdef T_PASSIVE
          DO i=1,NPT
            itrc=inert(i)
            IF (LtracerIsf(itrc,ng)) THEN
              WRITE (out,185) LtracerIsf(itrc,ng), 'LtracerIsf', itrc,  &
     &            'Turning ON  point Sources/Sinks on iceshelf tracer ',&
     &            i, TRIM(Vname(1,idItrc(itrc)))
            ELSE
              WRITE (out,185) LtracerIsf(itrc,ng), 'LtracerIsf', itrc,  &
     &            'Turning OFF point Sources/Sinks on iceshelf tracer ',&
     &            i, TRIM(Vname(1,idItrc(itrc)))
            END IF
          END DO
# endif
#endif
#ifdef ICESHEET
          WRITE (out,60) n_glen(ng), 'n_glen',                          &
     &          'SSA Glen flow law exponent.'
          WRITE (out,60) A_glen(ng), 'A_glen',                          &
     &          'SSA Glen rate factor [Pa^-3 s^-1].'
          WRITE (out,60) eps_glen_min(ng), 'eps_glen_min',              &
     &          'SSA minimum strain rate regularization [1/s].'
          WRITE (out,60) C_basal(ng), 'C_basal',                        &
     &          'SSA basal drag coefficient [Pa (s/m)^(1/m)].'
          WRITE (out,60) n_basal(ng), 'n_basal',                        &
     &          'SSA basal sliding exponent.'
          WRITE (out,60) eps_u_min(ng), 'eps_u_min',                    &
     &          'SSA minimum sliding velocity [m/s].'
          WRITE (out,50) nl_max_iter(ng), 'nl_max_iter',                &
     &          'SSA max Picard iterations.'
          WRITE (out,60) nl_tol(ng), 'nl_tol',                          &
     &          'SSA Picard relative tolerance.'
          WRITE (out,60) nl_tol_fp(ng), 'nl_tol_fp',                    &
     &          'SSA Picard absolute velocity tolerance [m/s].'
          WRITE (out,50) cg_max_iter(ng), 'cg_max_iter',                &
     &          'SSA max CG iterations.'
          WRITE (out,60) cg_tol(ng), 'cg_tol',                          &
     &          'SSA CG relative tolerance.'
          WRITE (out,60) nuH_reg(ng), 'nuH_reg',                        &
     &          'SSA nuH regularization [Pa s m].'
          WRITE (out,60) min_thickness(ng), 'min_thickness',            &
     &          'SSA minimum ice thickness [m].'
          WRITE (out,60) ssa_max_speed(ng), 'ssa_max_speed',            &
     &          'SSA velocity cap [m/s].'
          WRITE (out,50) n_smooth(ng), 'n_smooth',                      &
     &          'Surface elevation smoothing iterations.'
          WRITE (out,50) n_taper(ng), 'n_taper',                        &
     &          'Ice margin taper width [cells].'
          WRITE (out,60) picard_alpha(ng), 'picard_alpha',              &
     &          'Picard relaxation factor (0 < alpha <= 1).'
          WRITE (out,60) cfl_factor(ng), 'cfl_factor',                  &
     &          'Thickness advection CFL safety factor [-].'
          WRITE (out,60) H_calve(ng), 'H_calve',                        &
     &          'Calving thickness threshold [m].'
          WRITE (out,50) N_ssa(ng), 'N_ssa',                            &
     &          'SSA solve frequency [ocean steps].'
          WRITE (out,50) N_thick(ng), 'N_thick',                        &
     &          'Thickness update frequency [SSA steps].'
          WRITE (out,50) N_gl(ng), 'N_gl',                              &
     &          'Grounding line update frequency [SSA steps].'
          WRITE (out,60) smb_rate(ng), 'smb_rate',                      &
     &          'Surface mass balance rate [m/s ice equiv, +=accum].'
          WRITE (out,60) gl_float_width(ng), 'gl_float_width',          &
     &          'Grounding line float_frac transition width [m].'
          WRITE (out,60) max_thickness(ng), 'max_thickness',            &
     &          'Maximum ice thickness safety cap [m].'
#endif
!
          IF (nHIS(ng).gt.0) THEN
            WRITE (out,'(1x)')
#ifdef ICESHEET
            IF (Hout(idIhis,ng)) WRITE (out,70) Hout(idIhis,ng),        &
     &         'Hout(idIhis)',                                          &
     &         'Write out iceshelf thickness.'
            IF (Hout(idIuis,ng)) WRITE (out,70) Hout(idIuis,ng),        &
     &         'Hout(idIuis)',                                          &
     &         'Write out iceshelf U-velocity.'
            IF (Hout(idIvis,ng)) WRITE (out,70) Hout(idIvis,ng),        &
     &         'Hout(idIvis)',                                          &
     &         'Write out iceshelf V-velocity.'
            IF (Hout(idIuer,ng)) WRITE (out,70) Hout(idIuer,ng),        &
     &         'Hout(idIuer)',                                          &
     &         'Write out iceshelf eastward velocity.'
            IF (Hout(idIvnr,ng)) WRITE (out,70) Hout(idIvnr,ng),        &
     &         'Hout(idIvnr)',                                          &
     &         'Write out iceshelf northward velocity.'
            IF (Hout(idItis,ng)) WRITE (out,70) Hout(idItis,ng),        &
     &         'Hout(idItis)',                                          &
     &         'Write out iceshelf temperature.'
            IF (Hout(idIeis,ng)) WRITE (out,70) Hout(idIeis,ng),        &
     &         'Hout(idIeis)',                                          &
     &         'Write out iceshelf enthalpy.'
            IF (Hout(idIzsf,ng)) WRITE (out,70) Hout(idIzsf,ng),        &
     &         'Hout(idIzsf)',                                          &
     &         'Write out ice sheet surface elevation.'
            IF (Hout(idIzbd,ng)) WRITE (out,70) Hout(idIzbd,ng),        &
     &         'Hout(idIzbd)',                                          &
     &         'Write out ice sheet bed elevation.'
            IF (Hout(idInuh,ng)) WRITE (out,70) Hout(idInuh,ng),        &
     &         'Hout(idInuh)',                                          &
     &         'Write out depth-integrated ice effective viscosity.'
            IF (Hout(idItdx,ng)) WRITE (out,70) Hout(idItdx,ng),        &
     &         'Hout(idItdx)',                                          &
     &         'Write out ice driving stress x-component.'
            IF (Hout(idItdy,ng)) WRITE (out,70) Hout(idItdy,ng),        &
     &         'Hout(idItdy)',                                          &
     &         'Write out ice driving stress y-component.'
            IF (Hout(idIbta,ng)) WRITE (out,70) Hout(idIbta,ng),        &
     &         'Hout(idIbta)',                                          &
     &         'Write out ice basal drag coefficient.'
#endif
#ifdef ICESHELF_THERMO
            IF (Hout(idImel,ng)) WRITE (out,70) Hout(idImel,ng),        &
     &         'Hout(idImel)',                                          &
     &         'Write out iceshelf basal melt rate.'
            IF (Hout(idIgat,ng)) WRITE (out,70) Hout(idIgat,ng),        &
     &         'Hout(idIgat)',                                          &
     &         'Write out iceshelf gammaT value.'
            IF (Hout(idIgas,ng)) WRITE (out,70) Hout(idIgas,ng),        &
     &         'Hout(idIgas)',                                          &
     &         'Write out iceshelf gammaS value.'
            IF (Hout(idItst,ng)) WRITE (out,70) Hout(idItst,ng),        &
     &         'Hout(idItst)',                                          &
     &         'Write out iceshelf Tstar value.'
            IF (Hout(idIust,ng)) WRITE (out,70) Hout(idIust,ng),        &
     &         'Hout(idIust)',                                          &
     &         'Write out iceshelf Ustar value.'
            IF (Hout(idItbo,ng)) WRITE (out,70) Hout(idItbo,ng),        &
     &         'Hout(idItbo)',                                          &
     &         'Write out iceshelf laminar boundary layer temperature.'
            IF (Hout(idIsbo,ng)) WRITE (out,70) Hout(idIsbo,ng),        &
     &         'Hout(idIsbo)',                                          &
     &         'Write out iceshelf laminar boundary layer salinity.'
# ifdef T_PASSIVE
            IF (Hout(idIvol,ng)) WRITE (out,70) Hout(idIvol,ng),        &
     &         'Hout(idIvol)',                                          &
     &         'Write out iceshelf refreezed volume.'
            DO i=1,NPT
              itrc=inert(i)
              IF (Hout(idItrc(itrc),ng)) WRITE (out,180)                &
     &            Hout(idItrc(itrc),ng), 'Hout(idItrc)',                &
     &            'Write out iceshelf refreezed passive tracer ', i,    &
     &            TRIM(Vname(1,idItrc(itrc)))
            END DO
            DO i=1,NPT
              itrc=inert(i)
              IF (Hout(idIamb(itrc),ng)) WRITE (out,180)                &
     &            Hout(idIamb(itrc),ng), 'Hout(idIamb)',                &
     &            'Write out iceshelf ambient passive tracer ', i,      &
     &            TRIM(Vname(1,idIamb(itrc)))
            END DO
# endif
#endif
          END IF
!
          IF (nQCK(ng).gt.0) THEN
            WRITE (out,'(1x)')
#ifdef ICESHEET
            IF (Qout(idIhis,ng)) WRITE (out,70) Qout(idIhis,ng),        &
     &         'Qout(idIhis)',                                          &
     &         'Write out iceshelf thickness.'
            IF (Qout(idIuis,ng)) WRITE (out,70) Qout(idIuis,ng),        &
     &         'Qout(idIuis)',                                          &
     &         'Write out iceshelf U-velocity.'
            IF (Qout(idIvis,ng)) WRITE (out,70) Qout(idIvis,ng),        &
     &         'Qout(idIvis)',                                          &
     &         'Write out iceshelf V-velocity.'
            IF (Qout(idIuer,ng)) WRITE (out,70) Qout(idIuer,ng),        &
     &         'Qout(idIuer)',                                          &
     &         'Write out iceshelf eastward velocity.'
            IF (Qout(idIvnr,ng)) WRITE (out,70) Qout(idIvnr,ng),        &
     &         'Qout(idIvnr)',                                          &
     &         'Write out iceshelf northward velocity.'
            IF (Qout(idIeis,ng)) WRITE (out,70) Qout(idIeis,ng),        &
     &         'Qout(idIeis)',                                          &
     &         'Write out iceshelf enthalpy.'
            IF (Qout(idItis,ng)) WRITE (out,70) Qout(idItis,ng),        &
     &         'Qout(idItis)',                                          &
     &         'Write out ice temperature.'
            IF (Qout(idIzsf,ng)) WRITE (out,70) Qout(idIzsf,ng),        &
     &         'Qout(idIzsf)',                                          &
     &         'Write out ice sheet surface elevation.'
            IF (Qout(idIzbd,ng)) WRITE (out,70) Qout(idIzbd,ng),        &
     &         'Qout(idIzbd)',                                          &
     &         'Write out ice sheet bed elevation.'
            IF (Qout(idInuh,ng)) WRITE (out,70) Qout(idInuh,ng),        &
     &         'Qout(idInuh)',                                          &
     &         'Write out depth-integrated ice effective viscosity.'
            IF (Qout(idItdx,ng)) WRITE (out,70) Qout(idItdx,ng),        &
     &         'Qout(idItdx)',                                          &
     &         'Write out ice driving stress x-component.'
            IF (Qout(idItdy,ng)) WRITE (out,70) Qout(idItdy,ng),        &
     &         'Qout(idItdy)',                                          &
     &         'Write out ice driving stress y-component.'
            IF (Qout(idIbta,ng)) WRITE (out,70) Qout(idIbta,ng),        &
     &         'Qout(idIbta)',                                          &
     &         'Write out ice basal drag coefficient.'
#endif
#ifdef ICESHELF_THERMO
            IF (Qout(idImel,ng)) WRITE (out,70) Qout(idImel,ng),        &
     &         'Qout(idImel)',                                          &
     &         'Write out iceshelf basal melt rate.'
            IF (Qout(idIgat,ng)) WRITE (out,70) Qout(idIgat,ng),        &
     &         'Qout(idIgat)',                                          &
     &         'Write out iceshelf gammaT value.'
            IF (Qout(idIgas,ng)) WRITE (out,70) Qout(idIgas,ng),        &
     &         'Qout(idIgas)',                                          &
     &         'Write out iceshelf gammaS value.'
            IF (Qout(idItst,ng)) WRITE (out,70) Qout(idItst,ng),        &
     &         'Qout(idItst)',                                          &
     &         'Write out iceshelf Tstar value.'
            IF (Qout(idIust,ng)) WRITE (out,70) Qout(idIust,ng),        &
     &         'Qout(idIust)',                                          &
     &         'Write out iceshelf Ustar value.'
            IF (Qout(idItbo,ng)) WRITE (out,70) Qout(idItbo,ng),        &
     &         'Qout(idItbo)',                                          &
     &         'Write out iceshelf laminar boundary layer temperature.'
            IF (Qout(idIsbo,ng)) WRITE (out,70) Qout(idIsbo,ng),        &
     &         'Qout(idIsbo)',                                          &
     &         'Write out iceshelf laminar boundary layer salinity.'
# ifdef T_PASSIVE
            IF (Qout(idIvol,ng)) WRITE (out,70) Qout(idIvol,ng),        &
     &         'Qout(idIvol)',                                          &
     &         'Write out iceshelf refreezed volume.'
            DO i=1,NPT
              itrc=inert(i)
              IF (Qout(idItrc(itrc),ng)) WRITE (out,180)                &
     &            Qout(idItrc(itrc),ng), 'Qout(idItrc)',                &
     &            'Write out iceshelf refreezed passive tracer ', i,    &
     &            TRIM(Vname(1,idItrc(itrc)))
            END DO
            DO i=1,NPT
              itrc=inert(i)
              IF (Qout(idIamb(itrc),ng)) WRITE (out,180)                &
     &            Qout(idIamb(itrc),ng), 'Qout(idIamb)',                &
     &            'Write out iceshelf ambient passive tracer ', i,      &
     &            TRIM(Vname(1,idIamb(itrc)))
            END DO
# endif
#endif
          END IF
#ifdef AVERAGES
!
          IF (nAVG(ng).gt.0) THEN
            WRITE (out,'(1x)')
# ifdef ICESHEET
            IF (Aout(idIhis,ng)) WRITE (out,70) Aout(idIhis,ng),        &
     &         'Aout(idIhis)',                                          &
     &         'Write out iceshelf thickness.'
            IF (Aout(idIuis,ng)) WRITE (out,70) Aout(idIuis,ng),        &
     &         'Aout(idIuis)',                                          &
     &         'Write out iceshelf U-velocity.'
            IF (Aout(idIvis,ng)) WRITE (out,70) Aout(idIvis,ng),        &
     &         'Aout(idIvis)',                                          &
     &         'Write out iceshelf V-velocity.'
            IF (Aout(idItis,ng)) WRITE (out,70) Aout(idItis,ng),        &
     &         'Aout(idItis)',                                          &
     &         'Write out iceshelf temperature.'
            IF (Aout(idIuer,ng)) WRITE (out,70) Aout(idIuer,ng),        &
     &         'Aout(idIuer)',                                          &
     &         'Write out iceshelf eastward velocity.'
            IF (Aout(idIvnr,ng)) WRITE (out,70) Aout(idIvnr,ng),        &
     &         'Aout(idIvnr)',                                          &
     &         'Write out iceshelf northward velocity.'
            IF (Aout(idIeis,ng)) WRITE (out,70) Aout(idIeis,ng),        &
     &         'Aout(idIeis)',                                          &
     &         'Write out iceshelf enthalpy.'
            IF (Aout(idIzsf,ng)) WRITE (out,70) Aout(idIzsf,ng),        &
     &         'Aout(idIzsf)',                                          &
     &         'Write out ice sheet surface elevation.'
            IF (Aout(idIzbd,ng)) WRITE (out,70) Aout(idIzbd,ng),        &
     &         'Aout(idIzbd)',                                          &
     &         'Write out ice sheet bed elevation.'
            IF (Aout(idInuh,ng)) WRITE (out,70) Aout(idInuh,ng),        &
     &         'Aout(idInuh)',                                          &
     &         'Write out depth-integrated ice effective viscosity.'
            IF (Aout(idItdx,ng)) WRITE (out,70) Aout(idItdx,ng),        &
     &         'Aout(idItdx)',                                          &
     &         'Write out ice driving stress x-component.'
            IF (Aout(idItdy,ng)) WRITE (out,70) Aout(idItdy,ng),        &
     &         'Aout(idItdy)',                                          &
     &         'Write out ice driving stress y-component.'
            IF (Aout(idIbta,ng)) WRITE (out,70) Aout(idIbta,ng),        &
     &         'Aout(idIbta)',                                          &
     &         'Write out ice basal drag coefficient.'
# endif
# ifdef ICESHELF_THERMO
            IF (Aout(idImel,ng)) WRITE (out,70) Aout(idImel,ng),        &
     &         'Aout(idImel)',                                          &
     &         'Write out iceshelf basal melt rate.'
            IF (Aout(idIgat,ng)) WRITE (out,70) Aout(idIgat,ng),        &
     &         'Aout(idIgat)',                                          &
     &         'Write out iceshelf gammaT value.'
            IF (Aout(idIgas,ng)) WRITE (out,70) Aout(idIgas,ng),        &
     &         'Aout(idIgas)',                                          &
     &         'Write out iceshelf gammaS value.'
            IF (Aout(idItst,ng)) WRITE (out,70) Aout(idItst,ng),        &
     &         'Aout(idItst)',                                          &
     &         'Write out iceshelf Tstar value.'
            IF (Aout(idIust,ng)) WRITE (out,70) Aout(idIust,ng),        &
     &         'Aout(idIust)',                                          &
     &         'Write out iceshelf Ustar value.'
            IF (Aout(idItbo,ng)) WRITE (out,70) Aout(idItbo,ng),        &
     &         'Aout(idItbo)',                                          &
     &         'Write out iceshelf laminar boundary layer temperature.'
            IF (Aout(idIsbo,ng)) WRITE (out,70) Aout(idIsbo,ng),        &
     &         'Aout(idIsbo)',                                          &
     &         'Write out iceshelf laminar boundary layer salinity.'
#  ifdef T_PASSIVE
            IF (Aout(idIvol,ng)) WRITE (out,70) Aout(idIvol,ng),        &
     &         'Aout(idIvol)',                                          &
     &         'Write out iceshelf refreezed volume.'
            DO i=1,NPT
              itrc=inert(i)
              IF (Aout(idItrc(itrc),ng)) WRITE (out,180)                &
     &            Aout(idItrc(itrc),ng), 'Aout(idItrc)',                &
     &            'Write out iceshelf refreezed passive tracer ', i,    &
     &            TRIM(Vname(1,idItrc(itrc)))
            END DO
            DO i=1,NPT
              itrc=inert(i)
              IF (Aout(idIamb(itrc),ng)) WRITE (out,180)                &
     &            Aout(idIamb(itrc),ng), 'Aout(idIamb)',                &
     &            'Write out iceshelf ambient passive tracer ', i,      &
     &            TRIM(Vname(1,idIamb(itrc)))
            END DO
#  endif
# endif
          END IF
#endif
#ifdef STATIONS
          IF (Lstations(ng).and.nSTA(ng).gt.0) THEN
            WRITE (out,'(1x)')
# ifdef ICESHEET
            IF (Sout(idIhis,ng)) WRITE (out,70) Sout(idIhis,ng),        &
     &         'Sout(idIhis)',                                          &
     &         'Write out iceshelf thickness.'
            IF (Sout(idIuis,ng)) WRITE (out,70) Sout(idIuis,ng),        &
     &         'Sout(idIuis)',                                          &
     &         'Write out iceshelf U-velocity.'
            IF (Sout(idIvis,ng)) WRITE (out,70) Sout(idIvis,ng),        &
     &         'Sout(idIvis)',                                          &
     &         'Write out iceshelf V-velocity.'
            IF (Sout(idIuer,ng)) WRITE (out,70) Sout(idIuer,ng),        &
     &         'Sout(idIuer)',                                          &
     &         'Write out iceshelf eastward velocity.'
            IF (Sout(idIvnr,ng)) WRITE (out,70) Sout(idIvnr,ng),        &
     &         'Sout(idIvnr)',                                          &
     &         'Write out iceshelf northward velocity.'
            IF (Sout(idItis,ng)) WRITE (out,70) Sout(idItis,ng),        &
     &         'Sout(idItis)',                                          &
     &         'Write out iceshelf temperature.'
            IF (Sout(idIeis,ng)) WRITE (out,70) Sout(idIeis,ng),        &
     &         'Sout(idIeis)',                                          &
     &         'Write out iceshelf enthalpy.'
            IF (Sout(idIzsf,ng)) WRITE (out,70) Sout(idIzsf,ng),        &
     &         'Sout(idIzsf)',                                          &
     &         'Write out ice sheet surface elevation.'
            IF (Sout(idIzbd,ng)) WRITE (out,70) Sout(idIzbd,ng),        &
     &         'Sout(idIzbd)',                                          &
     &         'Write out ice sheet bed elevation.'
            IF (Sout(idInuh,ng)) WRITE (out,70) Sout(idInuh,ng),        &
     &         'Sout(idInuh)',                                          &
     &         'Write out depth-integrated ice effective viscosity.'
            IF (Sout(idItdx,ng)) WRITE (out,70) Sout(idItdx,ng),        &
     &         'Sout(idItdx)',                                          &
     &         'Write out ice driving stress x-component.'
            IF (Sout(idItdy,ng)) WRITE (out,70) Sout(idItdy,ng),        &
     &         'Sout(idItdy)',                                          &
     &         'Write out ice driving stress y-component.'
            IF (Sout(idIbta,ng)) WRITE (out,70) Sout(idIbta,ng),        &
     &         'Sout(idIbta)',                                          &
     &         'Write out ice basal drag coefficient.'
# endif
# ifdef ICESHELF_THERMO
            IF (Sout(idImel,ng)) WRITE (out,70) Sout(idImel,ng),        &
     &         'Sout(idImel)',                                          &
     &         'Write out iceshelf basal melt rate.'
            IF (Sout(idIgat,ng)) WRITE (out,70) Sout(idIgat,ng),        &
     &         'Sout(idIgat)',                                          &
     &         'Write out iceshelf gammaT value.'
            IF (Sout(idIgas,ng)) WRITE (out,70) Sout(idIgas,ng),        &
     &         'Sout(idIgas)',                                          &
     &         'Write out iceshelf gammaS value.'
            IF (Sout(idItst,ng)) WRITE (out,70) Sout(idItst,ng),        &
     &         'Sout(idItst)',                                          &
     &         'Write out iceshelf Tstar value.'
            IF (Sout(idIust,ng)) WRITE (out,70) Sout(idIust,ng),        &
     &         'Sout(idIust)',                                          &
     &         'Write out iceshelf Ustar value.'
            IF (Sout(idItbo,ng)) WRITE (out,70) Sout(idItbo,ng),        &
     &         'Sout(idItbo)',                                          &
     &         'Write out iceshelf laminar boundary layer temperature.'
            IF (Sout(idIsbo,ng)) WRITE (out,70) Sout(idIsbo,ng),        &
     &         'Sout(idIsbo)',                                          &
     &         'Write out iceshelf laminar boundary layer salinity.'
#  ifdef T_PASSIVE
            IF (Sout(idIvol,ng)) WRITE (out,70) Sout(idIvol,ng),        &
     &         'Sout(idIvol)',                                          &
     &         'Write out iceshelf refreezed volume.'
            DO i=1,NPT
              itrc=inert(i)
              IF (Sout(idItrc(itrc),ng)) WRITE (out,180)                &
     &            Sout(idItrc(itrc),ng), 'Sout(idItrc)',                &
     &            'Write out iceshelf refreezed passive tracer ', i,    &
     &            TRIM(Vname(1,idItrc(itrc)))
            END DO
            DO i=1,NPT
              itrc=inert(i)
              IF (Sout(idIamb(itrc),ng)) WRITE (out,180)                &
     &            Sout(idIamb(itrc),ng), 'Sout(idIamb)',                &
     &            'Write out iceshelf ambient passive tracer ', i,      &
     &            TRIM(Vname(1,idIamb(itrc)))
            END DO
#  endif
# endif
          END IF
#endif
        END DO
      END IF

#ifdef AVERAGES
!
!  Turn out switches to process time-averaged iceshelf model state and
!  internal variables. It is necessary to avoid recursive dependency
!  between "mod_ncparam" and "mod_isf".
!
      LisfFavg=.FALSE.
      LisfSavg=.FALSE.
!
      DO ng=1,Ngrids
        DO nf=1,nIsfF
          IF (iFisf(nf).gt.0) THEN
            IF (Aout(iFisf(nf),ng)) THEN
              LisfFavg(nf,ng)=.TRUE.
            END IF
          END IF
        END DO
!
        DO ns=1,nIsfS
          IF (iSisf(ns).gt.0) THEN
            IF (Aout(iSisf(ns),ng)) THEN
              LisfSavg(ns,ng)=.TRUE.
            END IF
          END IF
        END DO
      END DO
# ifdef T_PASSIVE
!
      LisfVavg=.FALSE.
      LisfTavg=.FALSE.
      LisfAavg=.FALSE.
!
      DO ng=1,Ngrids
        IF (Aout(idIvol,ng)) THEN
          LisfVavg(ng)=.TRUE.
        END IF
        DO i=1,NPT
          itrc=inert(i)
          IF (Aout(idItrc(itrc),ng)) THEN
            LisfTavg(itrc,ng)=.TRUE.
          END IF
          IF (Aout(idIamb(itrc),ng)) THEN
            LisfAavg(itrc,ng)=.TRUE.
          END IF
        END DO
      END DO
# endif
#endif
!
  30  FORMAT (/,' READ_IsfPar - Error while processing line: ',/,a)
  40  FORMAT (/,/,' Iceshelf Parameters, Grid: ',i2.2,                  &
     &        /,  ' ========================',/)
  50  FORMAT (1x,i10,2x,a,t32,a)
  60  FORMAT (1p,e11.4,2x,a,t32,a)
  70  FORMAT (10x,l1,2x,a,t32,a)
  80  FORMAT (/,' READ_IsfPar - variable index not yet loaded, ', a)
  90  FORMAT (/,' READ_IsfPar - wrong number of values: ',i3,           &
      &        ' for string parameter: ',a)
 180  FORMAT (10x,l1,2x,a,t32,a,i2.2,':',1x,a)
 185  FORMAT (10x,l1,2x,a,'(',i2.2,')',t32,a,i2.2,':',1x,a)
!
      RETURN
      END SUBROUTINE read_IsfPar
