#    ifdef ICESHEET
!
!  Read iceshelf model state variables.
!
        DO i=1,nIsfS
          IF (iSisf(i).gt.0) THEN
            ifield=iSisf(i)
            IF (get_var(ifield)) THEN
              foundit=find_string(var_name, n_var,                      &
     &                            TRIM(Vname(1,ifield)), vindex)
              IF (foundit) THEN
                my_pioVar%vd=var_desc(vindex)
                SELECT CASE (i)
                  CASE (isIuis, isIvis)
                    IF (Perfect2D) THEN
                      my_pioVar%gtype=p3dvar
                    ELSE
                      my_pioVar%gtype=p2dvar
                    END IF
                    IF (KIND(ISF(ng)%Si).eq.8) THEN
                      my_pioVar%dkind=PIO_double
                      IF (Perfect2D) THEN
                        ioDesc => ioDesc_dp_p3dvar(ng)
                      ELSE
                        ioDesc => ioDesc_dp_p2dvar(ng)
                      END IF
                    ELSE
                      my_pioVar%dkind=PIO_real
                      IF (Perfect2D) THEN
                        ioDesc => ioDesc_sp_p3dvar(ng)
                      ELSE
                        ioDesc => ioDesc_sp_p2dvar(ng)
                      END IF
                    END IF
!
                    IF (Perfect2D) THEN
                      status=nf_fread3d(ng, IDmod, ncname, pioFile,     &
     &                                  Vname(1,ifield), my_pioVar,     &
     &                                  InpRec, ioDesc, Vsize,          &
     &                                  LBi, UBi, LBj, UBj, 1, 2,       &
     &                                  Fscl, Fmin, Fmax,               &
#     ifdef MASKING
     &                                  GRID(ng) % pmask,               &
#     endif
#     ifdef CHECKSUM
     &                                  ISF(ng) % Si(:,:,:,i),          &
     &                                  checksum = Fhash)
#     else
     &                                  ISF(ng) % Si(:,:,:,i))
#     endif
                    ELSE
                      status=nf_fread2d(ng, IDmod, ncname, pioFile,     &
     &                                  Vname(1,ifield), my_pioVar,     &
     &                                  InpRec, ioDesc, Vsize,          &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  Fscl, Fmin, Fmax,               &
#     ifdef MASKING
     &                                  GRID(ng) % pmask,               &
#     endif
#     ifdef CHECKSUM
     &                                  ISF(ng) % Si(:,:,Tindex,i),     &
     &                                  checksum = Fhash)
#     else
     &                                  ISF(ng) % Si(:,:,Tindex,i))
#     endif
                    END IF
                  CASE DEFAULT
                    IF (Perfect2D) THEN
                      my_pioVar%gtype=r3dvar
                    ELSE
                      my_pioVar%gtype=r2dvar
                    END IF
                    IF (KIND(ISF(ng)%Si).eq.8) THEN
                      my_pioVar%dkind=PIO_double
                      IF (Perfect2D) THEN
                        ioDesc => ioDesc_dp_r3dvar(ng)
                      ELSE
                        ioDesc => ioDesc_dp_r2dvar(ng)
                      END IF
                    ELSE
                      my_pioVar%dkind=PIO_real
                      IF (Perfect2D) THEN
                        ioDesc => ioDesc_sp_r3dvar(ng)
                      ELSE
                        ioDesc => ioDesc_sp_r2dvar(ng)
                      END IF
                    END IF
!
                    IF (Perfect2D) THEN
                      status=nf_fread3d(ng, IDmod, ncname, pioFile,     &
     &                                  Vname(1,ifield), my_pioVar,     &
     &                                  InpRec, ioDesc, Vsize,          &
     &                                  LBi, UBi, LBj, UBj, 1, 2,       &
     &                                  Fscl, Fmin, Fmax,               &
#     ifdef MASKING
     &                                  GRID(ng) % rmask,               &
#     endif
#     ifdef CHECKSUM
     &                                  ISF(ng) % Si(:,:,:,i),          &
     &                                  checksum = Fhash)
#     else
     &                                  ISF(ng) % Si(:,:,:,i))
#     endif
                    ELSE
                      status=nf_fread2d(ng, IDmod, ncname, pioFile,     &
     &                                  Vname(1,ifield), my_pioVar,     &
     &                                  InpRec, ioDesc, Vsize,          &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  Fscl, Fmin, Fmax,               &
#     ifdef MASKING
     &                                  GRID(ng) % rmask,               &
#     endif
#     ifdef CHECKSUM
     &                                  ISF(ng) % Si(:,:,Tindex,i),     &
     &                                  checksum = Fhash)
#     else
     &                                  ISF(ng) % Si(:,:,Tindex,i))
#     endif
                    END IF
                END SELECT
!
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) THEN
                    WRITE (stdout,60) string, TRIM(Vname(1,ifield)),    &
     &                                InpRec, TRIM(ncname)
                  END IF
                  exit_flag=2
                  ioerror=status
                  RETURN
                ELSE
                  IF (Master) THEN
#     ifdef CHECKSUM
                    WRITE (stdout,70) TRIM(Vname(2,ifield)),            &
     &                                Fmin, Fmax, Fhash
#     else
                    WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
#     endif
                  END IF
                END IF
              ELSE
                IF (Master) THEN
                  WRITE (stdout,80) string, TRIM(Vname(1,ifield)),      &
     &                              TRIM(ncname)
                END IF
                exit_flag=4
                IF (FoundError(exit_flag, PIO_noerr,                    &
     &                         __LINE__, MyFile)) THEN
                  RETURN
                END IF
              END IF
            END IF
          END IF
        END DO
!
!  Read the ice temperature (is_temp).  It is only flagged by
!  checkvars when the file does NOT carry the ice enthalpy
!  (temperature-only IC files), so the enthalpy state is always
!  derived here as E = h*T; files that carry the enthalpy read it
!  in the state loop above and never reach this block.
!
        IF (get_var(idItis)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idItis)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=r2dvar
            IF (KIND(ISF(ng)%Si).eq.8) THEN
              my_pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              my_pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
            status=nf_fread2d(ng, IDmod, ncname, pioFile,               &
     &                        Vname(1,idItis), my_pioVar,               &
     &                        InpRec, ioDesc, Vsize,                    &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#     ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#     endif
#     ifdef CHECKSUM
     &                        ISF(ng) % Fi(:,:,ifItis),                 &
     &                        checksum = Fhash)
#     else
     &                        ISF(ng) % Fi(:,:,ifItis))
#     endif
!
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idItis)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
#     ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idItis)), Fmin, Fmax,    &
     &                            Fhash
#     else
                WRITE (stdout,70) TRIM(Vname(2,idItis)), Fmin, Fmax
#     endif
              END IF
!
!  Convert temperature to enthalpy (E = h*T): both time levels on
!  perfect restart, else the loaded time level only (isf_ini_fields
!  copies it to the other level).
!
              IF (Perfect2D) THEN
                DO i=1,2
                  ISF(ng) % Si(:,:,i,isIeis) =                         &
     &                ISF(ng) % Si(:,:,i,isIhis) *                     &
     &                ISF(ng) % Fi(:,:,ifItis)
                END DO
              ELSE
                ISF(ng) % Si(:,:,Tindex,isIeis) =                       &
     &              ISF(ng) % Si(:,:,Tindex,isIhis) *                   &
     &              ISF(ng) % Fi(:,:,ifItis)
              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idItis)),          &
     &                          TRIM(ncname)
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, PIO_noerr,                        &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
#    endif
!
!  Read sea ice model state variables.
!
        DO i=1,nIsfF
#   ifdef ICESHEET
          IF (iFisf(i).gt.0 .and. i.ne.ifItis) THEN
#   else
          IF (iFisf(i).gt.0) THEN
#   endif
            ifield=iFisf(i)
            IF (get_var(ifield).and.Perfect2D) THEN
              foundit=find_string(var_name, n_var,                      &
     &                            TRIM(Vname(1,ifield)), vindex)
              IF (foundit) THEN
                my_pioVar%vd=var_desc(vindex)
                my_pioVar%gtype=r2dvar
                IF (KIND(ISF(ng)%Si).eq.8) THEN
                  my_pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_r2dvar(ng)
                ELSE
                  my_pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_r2dvar(ng)
                END IF
!
                status=nf_fread2d(ng, IDmod, ncname, pioFile,           &
     &                            Vname(1,ifield), my_pioVar,           &
     &                            InpRec, ioDesc, Vsize,                &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
#    ifdef MASKING
     &                            GRID(ng) % rmask,                     &
#    endif
#    ifdef CHECKSUM
     &                            ISF(ng) % Fi(:,:,i),                  &
     &                            checksum = Fhash)
#    else
     &                            ISF(ng) % Fi(:,:,i))
#    endif
!
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) THEN
                    WRITE (stdout,60) string, TRIM(Vname(1,ifield)),    &
     &                                InpRec, TRIM(ncname)
                  END IF
                  exit_flag=2
                  ioerror=status
                  RETURN
                ELSE
                  IF (Master) THEN
#    ifdef CHECKSUM
                    WRITE (stdout,70) TRIM(Vname(2,ifield)),            &
     &                                Fmin, Fmax, Fhash
#    else
                    WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
#    endif
                  END IF
                END IF
              ELSE
                IF (Master) THEN
                  WRITE (stdout,80) string, TRIM(Vname(1,ifield)),      &
     &                              TRIM(ncname)
                END IF
                exit_flag=4
                IF (FoundError(exit_flag, PIO_noerr,                    &
     &                         __LINE__, MyFile)) THEN
                  RETURN
                END IF
              END IF
            END IF
          END IF
        END DO
!
!  Read ice-ocean coupling (Ci) accumulators from restart file.
!
#  ifdef ICESHEET
        DO nc=1,nIsfC
          IF (iCisf(nc).gt.0) THEN
            ifield=iCisf(nc)
            IF (get_var(ifield).and.Perfect2D) THEN
              foundit=find_string(var_name, n_var,                      &
     &                            TRIM(Vname(1,ifield)), vindex)
              IF (foundit) THEN
                my_pioVar%vd=var_desc(vindex)
                my_pioVar%gtype=r2dvar
                IF (KIND(ISF(ng)%Ci).eq.8) THEN
                  my_pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_r2dvar(ng)
                ELSE
                  my_pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_r2dvar(ng)
                END IF
!
                status=nf_fread2d(ng, IDmod, ncname, pioFile,           &
     &                            Vname(1,ifield), my_pioVar,           &
     &                            InpRec, ioDesc, Vsize,                &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
#   ifdef MASKING
     &                            GRID(ng) % rmask,                     &
#   endif
#   ifdef CHECKSUM
     &                            ISF(ng) % Ci(:,:,nc),                 &
     &                            checksum = Fhash)
#   else
     &                            ISF(ng) % Ci(:,:,nc))
#   endif
!
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) THEN
                    WRITE (stdout,60) string, TRIM(Vname(1,ifield)),    &
     &                                InpRec, TRIM(ncname)
                  END IF
                  exit_flag=2
                  ioerror=status
                  RETURN
                ELSE
                  IF (Master) THEN
#   ifdef CHECKSUM
                    WRITE (stdout,70) TRIM(Vname(2,ifield)),            &
     &                                Fmin, Fmax, Fhash
#   else
                    WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
#   endif
                  END IF
                END IF
              END IF
            END IF
          END IF
        END DO
#  endif
!
!  Read SBL weights (Oi) from restart file (PERFECT_RESTART only).
!
#  if defined ICESHELF_THERMO && !defined ICESHELF_NO_SBL
#   ifdef PERFECT_RESTART
        DO no=1,nIsfO
          IF (iOisf(no).gt.0) THEN
            ifield=iOisf(no)
            IF (get_var(ifield).and.Perfect2D) THEN
              foundit=find_string(var_name, n_var,                      &
     &                            TRIM(Vname(1,ifield)), vindex)
              IF (foundit) THEN
                my_pioVar%vd=var_desc(vindex)
                my_pioVar%gtype=r3dvar
                IF (KIND(ISF(ng)%Oi).eq.8) THEN
                  my_pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_r3dvar(ng)
                ELSE
                  my_pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_r3dvar(ng)
                END IF
!
                status=nf_fread3d(ng, IDmod, ncname, pioFile,           &
     &                            Vname(1,ifield), my_pioVar,           &
     &                            InpRec, ioDesc, Vsize,                &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            Fscl, Fmin, Fmax,                     &
#   ifdef MASKING
     &                            GRID(ng) % rmask,                     &
#   endif
#   ifdef CHECKSUM
     &                            ISF(ng) % Oi(:,:,:,no),                &
     &                            checksum = Fhash)
#   else
     &                            ISF(ng) % Oi(:,:,:,no))
#   endif
!
                IF (FoundError(status, PIO_noerr, __LINE__, MyFile))    &
     &          THEN
                  IF (Master) THEN
                    WRITE (stdout,60) string, TRIM(Vname(1,ifield)),    &
     &                                InpRec, TRIM(ncname)
                  END IF
                  exit_flag=2
                  ioerror=status
                  RETURN
                ELSE
                  IF (Master) THEN
#   ifdef CHECKSUM
                    WRITE (stdout,70) TRIM(Vname(2,ifield)),            &
     &                                Fmin, Fmax, Fhash
#   else
                    WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
#   endif
                  END IF
                END IF
              END IF
            END IF
          END IF
        END DO
#   endif
#  endif
!
!  Read in iceshelf refreezed volume (m).  Only present under
!  ICESHELF_THERMO (refreezing is a thermodynamic process).
!
#  if defined ICESHELF_THERMO
        ifield=idIvol
        IF (get_var(ifield)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            IF (Perfect2D) THEN
              my_pioVar%gtype=r3dvar
            ELSE
              my_pioVar%gtype=r2dvar
            END IF
            IF (KIND(ISF(ng)%Vi).eq.8) THEN
              my_pioVar%dkind=PIO_double
              IF (Perfect2D) THEN
                ioDesc => ioDesc_dp_r3dvar(ng)
              ELSE
                ioDesc => ioDesc_dp_r2dvar(ng)
              END IF
            ELSE
              my_pioVar%dkind=PIO_real
              IF (Perfect2D) THEN
                ioDesc => ioDesc_sp_r3dvar(ng)
              ELSE
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
            END IF
            IF (Perfect2D) THEN
              status=nf_fread3d(ng, IDmod, ncname, pioFile,             &
     &                          Vname(1,ifield), my_pioVar,             &
     &                          InpRec, ioDesc, Vsize,                  &
     &                          LBi, UBi, LBj, UBj, 1, 2,               &
     &                          Fscl, Fmin, Fmax,                       &
#    ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#    endif
#    ifdef CHECKSUM
     &                          ISF(ng) % Vi,                           &
     &                          checksum = Fhash)
#    else
     &                          ISF(ng) % Vi)
#    endif
            ELSE
              status=nf_fread2d(ng, IDmod, ncname, pioFile,             &
     &                          Vname(1,ifield), my_pioVar,             &
     &                          InpRec, ioDesc, Vsize,                  &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
#    ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#    endif
#    ifdef CHECKSUM
     &                          ISF(ng) % Vi(:,:,Tindex),               &
     &                          checksum = Fhash)
#    else
     &                          ISF(ng) % Vi(:,:,Tindex))
#    endif
              END IF
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,ifield)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
#    ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax,    &
     &                            Fhash
#    else
                WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
#    endif

              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,ifield)),          &
     &                          TRIM(ncname)
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, PIO_noerr,                        &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
#  endif
!
!  Read in integrated tracer (m Tunits).  Only present under
!  ICESHELF_THERMO (ice tracers are thermodynamic tracers).
!
#  if defined ICESHELF_THERMO
        itrc=0
        DO i=NAT+1,NT(ng)
          ifield=idItrc(i)
          IF (get_var(ifield)) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,ifield)), vindex)
            IF (foundit) THEN
              itrc=itrc+1
              scale=1.0_dp
              my_pioVar%vd=var_desc(vindex)
              IF (Perfect2D) THEN
                my_pioVar%gtype=r3dvar
              ELSE
                my_pioVar%gtype=r2dvar
              END IF
              IF (KIND(ISF(ng)%Ti).eq.8) THEN
                my_pioVar%dkind=PIO_double
                IF (Perfect2D) THEN
                  ioDesc => ioDesc_dp_r3dvar(ng)
                ELSE
                  ioDesc => ioDesc_dp_r2dvar(ng)
                END IF
              ELSE
                my_pioVar%dkind=PIO_real
                IF (Perfect2D) THEN
                  ioDesc => ioDesc_sp_r3dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_r2dvar(ng)
                END IF
              END IF
!
              IF (Perfect2D) THEN
                status=nf_fread3d(ng, IDmod, ncname, pioFile,           &
     &                            Vname(1,ifield), my_pioVar,           &
     &                            InpRec, ioDesc, Vsize,                &
     &                            LBi, UBi, LBj, UBj, 1, 2,             &
     &                            scale, Fmin, Fmax,                    &
#    ifdef MASKING
     &                            GRID(ng) % rmask,                     &
#    endif
#    ifdef CHECKSUM
     &                            ISF(ng) % Ti(:,:,:,itrc),             &
     &                            checksum = Fhash)
#    else
     &                            ISF(ng) % Ti(:,:,:,itrc))
#    endif
              ELSE
                status=nf_fread2d(ng, IDmod, ncname, pioFile,           &
     &                            Vname(1,ifield), my_pioVar,           &
     &                            InpRec, ioDesc, Vsize,                &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            scale, Fmin, Fmax,                    &
#    ifdef MASKING
     &                            GRID(ng) % rmask,                     &
#    endif
#    ifdef CHECKSUM
     &                            ISF(ng) % Ti(:,:,Tindex,itrc),        &
     &                            checksum = Fhash)
#    else
     &                            ISF(ng) % Ti(:,:,Tindex,itrc))
#    endif
              END IF
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile))      &
     &          THEN
                IF (Master) THEN
                  WRITE (stdout,60) string,                             &
     &                              TRIM(Vname(1,ifield)),              &
     &                              InpRec, TRIM(ncname)
                END IF
                exit_flag=2
                ioerror=status
                RETURN
              ELSE
                IF (Master) THEN
#    ifdef CHECKSUM
                  WRITE (stdout,70) TRIM(Vname(2,ifield)),              &
     &                              Fmin, Fmax, Fhash
#    else
                  WRITE (stdout,70) TRIM(Vname(2,ifield)),              &
     &                              Fmin, Fmax
#    endif
                END IF
              END IF
            ELSE
              IF (Master) THEN
                WRITE (stdout,80) string, TRIM(Vname(1,ifield)),        &
     &                            TRIM(ncname)
              END IF
              exit_flag=4
              IF (FoundError(exit_flag, PIO_noerr,                      &
     &                       __LINE__, MyFile)) THEN
                RETURN
              END IF
            END IF
          END IF
        END DO
#  endif
