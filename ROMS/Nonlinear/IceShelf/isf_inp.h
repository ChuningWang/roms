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
      integer :: Npts, Nval, ng, status, soff
      integer :: igrid, nline
#ifdef AVERAGES
      integer :: nf, ns
#endif
      integer :: itrc
#ifdef ICESHELF_TRACER
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
!  Allocate ice model parameters that depend on Ngrids/
!-----------------------------------------------------------------------
!
#ifdef ICESHELF_TRACER
      IF (.not.allocated(Linert)) THEN
        allocate (Linert(NPT,Ngrids))
      END IF
#endif
!
      igrid=1                            ! nested grid counter
      nline=0                            ! LBC multi-line counter
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
#ifdef ICESHELF_MOM
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
            CASE ('LBC(isItis)')
              soff = isf_ic_brylast + isItis
              Npts=load_lbc(Nval, Cval, line, nline, soff, igrid,       &
     &                      0, 0, Vname(1,idItis), LBC)
#endif
#ifdef ICESHELF_THERMO
            CASE ('gamU')
              Npts=load_r(Nval, Rval, Ngrids, gamU)
            CASE ('min_ustar')
              Npts=load_r(Nval, Rval, Ngrids, min_ustar)
# ifdef ICESHELF_TRACER
            CASE ('LtracerIsf')
              Npts=load_l(Nval, Cval, NPT, Ngrids, Linert)
              DO ng=1,Ngrids
                DO itrc=1,NPT
                  LtracerIsf(itrc,ng)=Linert(itrc,ng)
                END DO
              END DO
# endif
#endif
#ifdef ICESHELF_MOM
            CASE ('Hout(idIhis)')
              IF (idIhis.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIhis'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idIhis,:))
            CASE ('Hout(idIuis)')
              IF (idIuis.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIuis'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idIuis,:))
            CASE ('Hout(idIvis)')
              IF (idIvis.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIvis'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idIvis,:))
            CASE ('Hout(idIuer)')
              IF (idIuer.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIuer'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idIuer,:))
            CASE ('Hout(idIvnr)')
              IF (idIvnr.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIvnr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idIvnr,:))
            CASE ('Hout(idItis)')
              IF (idItis.eq.0) THEN
                IF (Master) WRITE (out,80) 'idItis'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idItis,:))
            CASE ('Hout(idIeis)')
              IF (idIeis.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIeis'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idIeis,:))
#endif
#ifdef ICESHELF_THERMO
            CASE ('Hout(idImel)')
              IF (idImel.eq.0) THEN
                IF (Master) WRITE (out,80) 'idImel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idImel,:))
            CASE ('Hout(idIgat)')
              IF (idIgat.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIgat'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idIgat,:))
            CASE ('Hout(idIgas)')
              IF (idIgas.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIgas'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idIgas,:))
            CASE ('Hout(idItst)')
              IF (idItst.eq.0) THEN
                IF (Master) WRITE (out,80) 'idItst'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idItst,:))
            CASE ('Hout(idIust)')
              IF (idIust.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIust'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idIust,:))
            CASE ('Hout(idItbo)')
              IF (idItbo.eq.0) THEN
                IF (Master) WRITE (out,80) 'idItbo'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idItbo,:))
            CASE ('Hout(idIsbo)')
              IF (idIsbo.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIsbo'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idIsbo,:))
# ifdef ICESHELF_TRACER
            CASE ('Hout(idIvol)')
              IF (idIvol.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIvol'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idIvol,:))
            CASE ('Hout(idItrc)')
              IF (MAXVAL(idItrc).eq.0) THEN
                IF (Master) WRITE (out,280) 'idItrc'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NPT, Ngrids, Linert)
              DO ng=1,Ngrids
                DO itrc=1,NPT
                  Hout(itrc,ng)=Linert(itrc,ng)
                END DO
              END DO
# endif
#endif
#ifdef ICESHELF_MOM
            CASE ('Qout(idIhis)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idIhis,:))
            CASE ('Qout(idIuis)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idIuis,:))
            CASE ('Qout(idIvis)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idIvis,:))
            CASE ('Qout(idIuer)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idIuer,:))
            CASE ('Qout(idIvnr)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idIvnr,:))
            CASE ('Qout(idItis)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idItis,:))
            CASE ('Qout(idIeis)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idIeis,:))
#endif
#ifdef ICESHELF_THERMO
            CASE ('Qout(idImel)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idImel,:))
            CASE ('Qout(idIgat)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idIgat,:))
            CASE ('Qout(idIgas)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idIgas,:))
            CASE ('Qout(idItst)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idItst,:))
            CASE ('Qout(idIust)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idIust,:))
            CASE ('Qout(idItbo)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idItbo,:))
            CASE ('Qout(idIsbo)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idIsbo,:))
# ifdef ICESHELF_TRACER
            CASE ('Qout(idIvol)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idIvol,:))
            CASE ('Qout(idItrc)')
              Npts=load_l(Nval, Cval, NPT, Ngrids, Linert)
              DO ng=1,Ngrids
                DO itrc=1,NPT
                  Qout(itrc,ng)=Linert(itrc,ng)
                END DO
              END DO
# endif
#endif
#ifdef AVERAGES
# ifdef ICESHELF_MOM
            CASE ('Aout(idIhis)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idIhis,:))
            CASE ('Aout(idIuis)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idIuis,:))
            CASE ('Aout(idIvis)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idIvis,:))
            CASE ('Aout(idIuer)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idIuer,:))
            CASE ('Aout(idIvnr)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idIvnr,:))
            CASE ('Aout(idItis)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idItis,:))
            CASE ('Aout(idIeis)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idIeis,:))
# endif
# ifdef ICESHELF_THERMO
            CASE ('Aout(idImel)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idImel,:))
            CASE ('Aout(idIgat)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idIgat,:))
            CASE ('Aout(idIgas)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idIgas,:))
            CASE ('Aout(idItst)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idItst,:))
            CASE ('Aout(idIust)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idIust,:))
            CASE ('Aout(idItbo)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idItbo,:))
            CASE ('Aout(idIsbo)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idIsbo,:))
#  ifdef ICESHELF_TRACER
            CASE ('Aout(idIvol)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idIvol,:))
            CASE ('Aout(idItrc)')
              Npts=load_l(Nval, Cval, NPT, Ngrids, Linert)
              DO ng=1,Ngrids
                DO itrc=1,NPT
                  Aout(itrc,ng)=Linert(itrc,ng)
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
#ifdef ICESHELF_THERMO
          WRITE (out,60) gamU(ng), 'gamU',                              &
     &          'Iceshelf Ustar scaling factor (nondimensional).'
          WRITE (out,60) min_ustar(ng), 'min_ustar',                    &
     &          'Iceshelf minimum Ustar value.'
# ifdef ICESHELF_TRACER
          DO itrc=1,NPT
            IF (LtracerIsf(itrc,ng)) THEN
              WRITE (out,185) LtracerIsf(itrc,ng), 'LtracerIsf', itrc,  &
     &            'Turning ON  point Sources/Sinks on iceshelf tracer ',&
     &            itrc, TRIM(Vname(1,idItrc(itrc)))
            ELSE
              WRITE (out,185) LtracerIsf(itrc,ng), 'LtracerIsf', itrc,  &
     &            'Turning OFF point Sources/Sinks on iceshelf tracer ',&
     &            itrc, TRIM(Vname(1,idItrc(itrc)))
            END IF
          END DO
# endif
#endif
!
#ifdef ICESHELF_MOM
          IF (Hout(idIhis,ng)) WRITE (out,70) Hout(idIhis,ng),          &
     &       'Hout(idIhis)',                                            &
     &       'Write out iceshelf thickness.'
          IF (Hout(idIuis,ng)) WRITE (out,70) Hout(idIuis,ng),          &
     &       'Hout(idIuis)',                                            &
     &       'Write out iceshelf U-velocity.'
          IF (Hout(idIvis,ng)) WRITE (out,70) Hout(idIvis,ng),          &
     &       'Hout(idIvis)',                                            &
     &       'Write out iceshelf V-velocity.'
          IF (Hout(idIuer,ng)) WRITE (out,70) Hout(idIuer,ng),          &
     &       'Hout(idIuer)',                                            &
     &       'Write out iceshelf eastward velocity.'
          IF (Hout(idIvnr,ng)) WRITE (out,70) Hout(idIvnr,ng),          &
     &       'Hout(idIvnr)',                                            &
     &       'Write out iceshelf northward velocity.'
          IF (Hout(idItis,ng)) WRITE (out,70) Hout(idItis,ng),          &
     &       'Hout(idItis)',                                            &
     &       'Write out iceshelf temperature.'
          IF (Hout(idIeis,ng)) WRITE (out,70) Hout(idIeis,ng),          &
     &       'Hout(idIeis)',                                            &
     &       'Write out iceshelf enthalpy.'
#endif
#ifdef ICESHELF_THERMO
          IF (Hout(idImel,ng)) WRITE (out,70) Hout(idImel,ng),          &
     &       'Hout(idImel)',                                            &
     &       'Write out iceshelf basal melt rate.'
          IF (Hout(idIgat,ng)) WRITE (out,70) Hout(idIgat,ng),          &
     &       'Hout(idIgat)',                                            &
     &       'Write out iceshelf gammaT value.'
          IF (Hout(idIgas,ng)) WRITE (out,70) Hout(idIgas,ng),          &
     &       'Hout(idIgas)',                                            &
     &       'Write out iceshelf gammaS value.'
          IF (Hout(idItst,ng)) WRITE (out,70) Hout(idItst,ng),          &
     &       'Hout(idItst)',                                            &
     &       'Write out iceshelf Tstar value.'
          IF (Hout(idIust,ng)) WRITE (out,70) Hout(idIust,ng),          &
     &       'Hout(idIust)',                                            &
     &       'Write out iceshelf Ustar value.'
          IF (Hout(idItbo,ng)) WRITE (out,70) Hout(idItbo,ng),          &
     &       'Hout(idItbo)',                                            &
     &       'Write out iceshelf laminar boundary layer temperature.'
          IF (Hout(idIsbo,ng)) WRITE (out,70) Hout(idIsbo,ng),          &
     &       'Hout(idIsbo)',                                            &
     &       'Write out iceshelf laminar boundary layer salinity.'
# ifdef ICESHELF_TRACER
          IF (Hout(idIvol,ng)) WRITE (out,70) Hout(idIvol,ng),          &
     &       'Hout(idIvol)',                                            &
     &       'Write out iceshelf refreezed volume.'
          DO itrc=1,NPT
            IF (Hout(idItrc(itrc),ng)) WRITE (out,180)                  &
     &         Hout(idItrc(itrc),ng), 'Hout(idItrc)',                   &
     &         'Write out iceshelf refreezed passive tracer ', itrc,    &
     &         TRIM(Vname(1,idItrc(itrc)))                       
          END DO
# endif
#endif
!
#ifdef ICESHELF_MOM
          IF (Qout(idIhis,ng)) WRITE (out,70) Qout(idIhis,ng),          &
     &       'Qout(idIhis)',                                            &
     &       'Write out iceshelf thickness.'
          IF (Qout(idIuis,ng)) WRITE (out,70) Qout(idIuis,ng),          &
     &       'Qout(idIuis)',                                            &
     &       'Write out iceshelf U-velocity.'
          IF (Qout(idIvis,ng)) WRITE (out,70) Qout(idIvis,ng),          &
     &       'Qout(idIvis)',                                            &
     &       'Write out iceshelf V-velocity.'
          IF (Qout(idIuer,ng)) WRITE (out,70) Qout(idIuer,ng),          &
     &       'Qout(idIuer)',                                            &
     &       'Write out iceshelf eastward velocity.'
          IF (Qout(idIvnr,ng)) WRITE (out,70) Qout(idIvnr,ng),          &
     &       'Qout(idIvnr)',                                            &
     &       'Write out iceshelf northward velocity.'
          IF (Qout(idItis,ng)) WRITE (out,70) Qout(idItis,ng),          &
     &       'Qout(idItis)',                                            &
     &       'Write out iceshelf temperature.'
          IF (Qout(idIeis,ng)) WRITE (out,70) Qout(idIeis,ng),          &
     &       'Qout(idIeis)',                                            &
     &       'Write out iceshelf enthalpy.'
#endif
#ifdef ICESHELF_THERMO
          IF (Qout(idImel,ng)) WRITE (out,70) Qout(idImel,ng),          &
     &       'Qout(idImel)',                                            &
     &       'Write out iceshelf basal melt rate.'
          IF (Qout(idIgat,ng)) WRITE (out,70) Qout(idIgat,ng),          &
     &       'Qout(idIgat)',                                            &
     &       'Write out iceshelf gammaT value.'
          IF (Qout(idIgas,ng)) WRITE (out,70) Qout(idIgas,ng),          &
     &       'Qout(idIgas)',                                            &
     &       'Write out iceshelf gammaS value.'
          IF (Qout(idItst,ng)) WRITE (out,70) Qout(idItst,ng),          &
     &       'Qout(idItst)',                                            &
     &       'Write out iceshelf Tstar value.'
          IF (Qout(idIust,ng)) WRITE (out,70) Qout(idIust,ng),          &
     &       'Qout(idIust)',                                            &
     &       'Write out iceshelf Ustar value.'
          IF (Qout(idItbo,ng)) WRITE (out,70) Qout(idItbo,ng),          &
     &       'Qout(idItbo)',                                            &
     &       'Write out iceshelf laminar boundary layer temperature.'
          IF (Qout(idIsbo,ng)) WRITE (out,70) Qout(idIsbo,ng),          &
     &       'Qout(idIsbo)',                                            &
     &       'Write out iceshelf laminar boundary layer salinity.'
# ifdef ICESHELF_TRACER
          IF (Qout(idIvol,ng)) WRITE (out,70) Qout(idIvol,ng),          &
     &       'Qout(idIvol)',                                            &
     &       'Write out iceshelf refreezed volume.'
          DO itrc=1,NPT
            IF (Qout(idItrc(itrc),ng)) WRITE (out,180)                  &
     &          Qout(idItrc(itrc),ng), 'Qout(idItrc)',                  &
     &          'Write out iceshelf refreezed passive tracer ', itrc,   &
     &          TRIM(Vname(1,idItrc(itrc)))
          END DO
# endif
#endif
#ifdef AVERAGES
!
# ifdef ICESHELF_MOM
          IF (Aout(idIhis,ng)) WRITE (out,70) Aout(idIhis,ng),          &
     &       'Aout(idIhis)',                                            &
     &       'Write out iceshelf thickness.'
          IF (Aout(idIuis,ng)) WRITE (out,70) Aout(idIuis,ng),          &
     &       'Aout(idIuis)',                                            &
     &       'Write out iceshelf U-velocity.'
          IF (Aout(idIvis,ng)) WRITE (out,70) Aout(idIvis,ng),          &
     &       'Aout(idIvis)',                                            &
     &       'Write out iceshelf V-velocity.'
          IF (Aout(idItis,ng)) WRITE (out,70) Aout(idItis,ng),          &
     &       'Aout(idItis)',                                            &
     &       'Write out iceshelf temperature.'
          IF (Aout(idIuer,ng)) WRITE (out,70) Aout(idIuer,ng),          &
     &       'Aout(idIuer)',                                            &
     &       'Write out iceshelf eastward velocity.'
          IF (Aout(idIvnr,ng)) WRITE (out,70) Aout(idIvnr,ng),          &
     &       'Aout(idIvnr)',                                            &
     &       'Write out iceshelf northward velocity.'
          IF (Aout(idIeis,ng)) WRITE (out,70) Aout(idIeis,ng),          &
     &       'Aout(idIeis)',                                            &
     &       'Write out iceshelf enthalpy.'
# endif
# ifdef ICESHELF_THERMO
          IF (Aout(idImel,ng)) WRITE (out,70) Aout(idImel,ng),          &
     &       'Aout(idImel)',                                            &
     &       'Write out iceshelf basal melt rate.'
          IF (Aout(idIgat,ng)) WRITE (out,70) Aout(idIgat,ng),          &
     &       'Aout(idIgat)',                                            &
     &       'Write out iceshelf gammaT value.'
          IF (Aout(idIgas,ng)) WRITE (out,70) Aout(idIgas,ng),          &
     &       'Aout(idIgas)',                                            &
     &       'Write out iceshelf gammaS value.'
          IF (Aout(idItst,ng)) WRITE (out,70) Aout(idItst,ng),          &
     &       'Aout(idItst)',                                            &
     &       'Write out iceshelf Tstar value.'
          IF (Aout(idIust,ng)) WRITE (out,70) Aout(idIust,ng),          &
     &       'Aout(idIust)',                                            &
     &       'Write out iceshelf Ustar value.'
          IF (Aout(idItbo,ng)) WRITE (out,70) Aout(idItbo,ng),          &
     &       'Aout(idItbo)',                                            &
     &       'Write out iceshelf laminar boundary layer temperature.'
          IF (Aout(idIsbo,ng)) WRITE (out,70) Aout(idIsbo,ng),          &
     &       'Aout(idIsbo)',                                            &
     &       'Write out iceshelf laminar boundary layer salinity.'
#  ifdef ICESHELF_TRACER
          IF (Aout(idIvol,ng)) WRITE (out,70) Aout(idIvol,ng),          &
     &       'Aout(idIvol)',                                            &
     &       'Write out iceshelf refreezed volume.'
          DO itrc=1,NPT
            IF (Aout(idItrc(itrc),ng)) WRITE (out,180)                  &
     &          Aout(idItrc(itrc),ng), 'Aout(idItrc)',                  &
     &          'Write out iceshelf refreezed passive tracer ', itrc,   &
     &          TRIM(Vname(1,idItrc(itrc)))
          END DO
#  endif
# endif
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
# ifdef ICESHELF_TRACER
!
      LisfVavg=.FALSE.
      LisfTavg=.FALSE.
!
      DO ng=1,Ngrids
        IF (Aout(idIvol,ng)) THEN
          LisfVavg(ng)=.TRUE.
        END IF
!
        DO itrc=1,NPT
          IF (Aout(idItrc(itrc),ng)) THEN
            LisfTavg(itrc,ng)=.TRUE.
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
 180  FORMAT (10x,l1,2x,a,t32,a,i2.2,':',1x,a)
 185  FORMAT (10x,l1,2x,a,'(',i2.2,')',t32,a,i2.2,':',1x,a)
 280  FORMAT (/,' READ_PHYPAR - Variable index not yet loaded, ', a)
!
      RETURN
      END SUBROUTINE read_IsfPar
