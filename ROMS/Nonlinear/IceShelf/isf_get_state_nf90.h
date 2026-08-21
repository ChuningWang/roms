#   ifdef ICESHEET
!
!  Read iceshelf model state variables.
!
        DO i=1,nIsfS
          IF (iSisf(i).gt.0) THEN
            ifield=iSisf(i)
            IF (get_var(ifield)) THEN
              foundit=find_string(var_name, n_var,                      &
     &                            TRIM(Vname(1,ifield)), varid)
              IF (foundit) THEN
                SELECT CASE (i)
                  CASE (isIuis, isIvis)
                    IF (Perfect2D) THEN
                      gtype=var_flag(varid)*p3dvar
                    ELSE
                      gtype=var_flag(varid)*p2dvar
                    END IF
                    IF (Perfect2D) THEN
                      status=nf_fread3d(ng, IDmod, ncname, ncINPid,     &
     &                                  Vname(1,ifield), varid,         &
     &                                  InpRec, gtype, Vsize,           &
     &                                  LBi, UBi, LBj, UBj, 1, 2,       &
     &                                  Fscl, Fmin, Fmax,               &
#    ifdef MASKING
     &                                  GRID(ng) % pmask,               &
#    endif
#    ifdef CHECKSUM
     &                                  ISF(ng) % Si(:,:,:,i),          &
     &                                  checksum = Fhash)
#    else
     &                                  ISF(ng) % Si(:,:,:,i))
#    endif
                    ELSE
                      status=nf_fread2d(ng, IDmod, ncname, ncINPid,     &
     &                                  Vname(1,ifield), varid,         &
     &                                  InpRec, gtype, Vsize,           &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  Fscl, Fmin, Fmax,               &
#    ifdef MASKING
     &                                  GRID(ng) % pmask,               &
#    endif
#    ifdef CHECKSUM
     &                                  ISF(ng) % Si(:,:,Tindex,i),     &
     &                                  checksum = Fhash)
#    else
     &                                  ISF(ng) % Si(:,:,Tindex,i))
#    endif
                    END IF
                  CASE DEFAULT
                    IF (Perfect2D) THEN
                      gtype=var_flag(varid)*r3dvar
                    ELSE
                      gtype=var_flag(varid)*r2dvar
                    END IF
                    IF (Perfect2D) THEN
                      status=nf_fread3d(ng, IDmod, ncname, ncINPid,     &
     &                                  Vname(1,ifield), varid,         &
     &                                  InpRec, gtype, Vsize,           &
     &                                  LBi, UBi, LBj, UBj, 1, 2,       &
     &                                  Fscl, Fmin, Fmax,               &
#    ifdef MASKING
     &                                  GRID(ng) % rmask,               &
#    endif
#    ifdef CHECKSUM
     &                                  ISF(ng) % Si(:,:,:,i),          &
     &                                  checksum = Fhash)
#    else
     &                                  ISF(ng) % Si(:,:,:,i))
#    endif
                    ELSE
                      status=nf_fread2d(ng, IDmod, ncname, ncINPid,     &
     &                                  Vname(1,ifield), varid,         &
     &                                  InpRec, gtype, Vsize,           &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  Fscl, Fmin, Fmax,               &
#    ifdef MASKING
     &                                  GRID(ng) % rmask,               &
#    endif
#    ifdef CHECKSUM
     &                                  ISF(ng) % Si(:,:,Tindex,i),     &
     &                                  checksum = Fhash)
#    else
     &                                  ISF(ng) % Si(:,:,Tindex,i))
#    endif
                    END IF
                END SELECT
!
                IF (FoundError(status, nf90_noerr, __LINE__,            &
     &              MyFile)) THEN
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
                IF (FoundError(exit_flag, nf90_noerr,                   &
     &                         __LINE__, MyFile)) THEN
                  RETURN
                END IF
              END IF
            END IF
          END IF
        END DO
!
!  Read the ice temperature (is_temp) from the initial conditions
!  file and derive the ice enthalpy state variable, E = h*T.  IC
!  files carry the ice temperature instead of the enthalpy; restart
!  files carry the enthalpy state directly (nrrec.ne.0).
!
        IF (get_var(idItis).and.(nrrec(ng).eq.0)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idItis)),   &
     &                        varid)
          IF (foundit) THEN
            gtype=var_flag(varid)*r2dvar
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idItis), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#    ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#    endif
#    ifdef CHECKSUM
     &                        ISF(ng) % Fi(:,:,icItis),                 &
     &                        checksum = Fhash)
#    else
     &                        ISF(ng) % Fi(:,:,icItis))
#    endif
!
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idItis)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
#    ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,idItis)), Fmin, Fmax,    &
     &                            Fhash
#    else
                WRITE (stdout,70) TRIM(Vname(2,idItis)), Fmin, Fmax
#    endif
              END IF
!
!  Convert temperature to enthalpy (E = h*T) on the loaded time
!  level; isf_ini_fields copies it to the other time level.
!
              ISF(ng) % Si(:,:,Tindex,isIeis) =                         &
     &            ISF(ng) % Si(:,:,Tindex,isIhis) *                     &
     &            ISF(ng) % Fi(:,:,icItis)
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,idItis)),          &
     &                          TRIM(ncname)
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, nf90_noerr,                       &
     &                     __LINE__, MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
#   endif
!
!  Read iceshelf model state variables.
!
        DO i=1,nIsfF
          IF (iFisf(i).gt.0) THEN
            ifield=iFisf(i)
            IF (get_var(ifield).and.Perfect2D) THEN
              foundit=find_string(var_name, n_var,                      &
     &                            TRIM(Vname(1,ifield)), varid)
              IF (foundit) THEN
                gtype=var_flag(varid)*r2dvar
                status=nf_fread2d(ng, IDmod, ncname, ncINPid,           &
     &                            Vname(1,ifield), varid,               &
     &                            InpRec, gtype, Vsize,                 &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
#   ifdef MASKING
     &                            GRID(ng) % rmask,                     &
#   endif
#   ifdef CHECKSUM
     &                            ISF(ng) % Fi(:,:,i),                  &
     &                            checksum = Fhash)
#   else
     &                            ISF(ng) % Fi(:,:,i))
#   endif
!
                IF (FoundError(status, nf90_noerr, __LINE__,            &
     &              MyFile)) THEN
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
              ELSE
                IF (Master) THEN
                  WRITE (stdout,80) string, TRIM(Vname(1,ifield)),      &
     &                              TRIM(ncname)
                END IF
                exit_flag=4
                IF (FoundError(exit_flag, nf90_noerr,                   &
     &                         __LINE__, MyFile)) THEN
                  RETURN
                END IF
              END IF
            END IF
          END IF
        END DO
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
     &                            TRIM(Vname(1,ifield)), varid)
              IF (foundit) THEN
                gtype=var_flag(varid)*r3dvar
                status=nf_fread3d(ng, IDmod, ncname, ncINPid,           &
     &                            Vname(1,ifield), varid,               &
     &                            InpRec, gtype, Vsize,                 &
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
                IF (FoundError(status, nf90_noerr, __LINE__,            &
     &              MyFile)) THEN
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
     &                        varid)
          IF (foundit) THEN
            IF (Perfect2D) THEN
              gtype=var_flag(varid)*r3dvar
            ELSE
              gtype=var_flag(varid)*r2dvar
            END IF
            IF (Perfect2D) THEN
              status=nf_fread3d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,ifield), varid,                 &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj, 1, 2,               &
     &                          Fscl, Fmin, Fmax,                       &
#   ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#   endif
#   ifdef CHECKSUM
     &                          ISF(ng) % Vi,                           &
     &                          checksum = Fhash)
#   else
     &                          ISF(ng) % Vi)
#   endif
            ELSE
              status=nf_fread2d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,ifield), varid,                 &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
#   ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#   endif
#   ifdef CHECKSUM
     &                          ISF(ng) % Vi(:,:,Tindex),               &
     &                          checksum = Fhash)
#   else
     &                          ISF(ng) % Vi(:,:,Tindex))
#   endif
            END IF
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,ifield)),        &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
#   ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax,    &
     &                            Fhash
#   else
                WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
#   endif

              END IF
            END IF
          ELSE
            IF (Master) THEN
              WRITE (stdout,80) string, TRIM(Vname(1,ifield)),          &
     &                          TRIM(ncname)
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, nf90_noerr,                       &
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
     &                          TRIM(Vname(1,ifield)), varid)
            IF (foundit) THEN
              itrc=itrc+1
              IF (Perfect2D) THEN
                gtype=var_flag(varid)*r3dvar
              ELSE
                gtype=var_flag(varid)*r2dvar
              END IF
              IF (Perfect2D) THEN
                status=nf_fread3d(ng, IDmod, ncname, ncINPid,           &
     &                            Vname(1,ifield), varid,               &
     &                            InpRec, gtype, Vsize,                 &
     &                            LBi, UBi, LBj, UBj, 1, 2,             &
     &                            Fscl, Fmin, Fmax,                     &
#   ifdef MASKING
     &                            GRID(ng) % rmask,                     &
#   endif
#   ifdef CHECKSUM
     &                            ISF(ng) % Ti(:,:,:,itrc),             &
     &                            checksum = Fhash)
#   else
     &                            ISF(ng) % Ti(:,:,:,itrc))
#   endif
              ELSE
                status=nf_fread2d(ng, IDmod, ncname, ncINPid,           &
     &                            Vname(1,ifield), varid,               &
     &                            InpRec, gtype, Vsize,                 &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscl, Fmin, Fmax,                     &
#   ifdef MASKING
     &                            GRID(ng) % rmask,                     &
#   endif
#   ifdef CHECKSUM
     &                            ISF(ng) % Ti(:,:,Tindex,itrc),        &
     &                            checksum = Fhash)
#   else
     &                            ISF(ng) % Ti(:,:,Tindex,itrc))
#   endif
              END IF
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile))     &
     &          THEN
                IF (Master) THEN
                  WRITE (stdout,60) string, TRIM(Vname(1,ifield)),      &
     &                              InpRec, TRIM(ncname)
                END IF
                exit_flag=2
                ioerror=status
                RETURN
              ELSE
                IF (Master) THEN
#   ifdef CHECKSUM
                  WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax,  &
     &                              Fhash
#   else
                  WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
#   endif
                END IF
              END IF
            ELSE
              IF (Master) THEN
                WRITE (stdout,80) string, TRIM(Vname(1,ifield)),        &
     &                            TRIM(ncname)
              END IF
              exit_flag=4
              IF (FoundError(exit_flag, nf90_noerr,                     &
     &                       __LINE__, MyFile)) THEN
                RETURN
              END IF
            END IF
          END IF
        END DO
#  endif
