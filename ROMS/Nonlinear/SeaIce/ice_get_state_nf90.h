!
!  Read sea ice model state variables.
!
        DO i=1,nIceS
          IF (iSice(i).gt.0) THEN
            ifield=iSice(i)
            IF (get_var(ifield)) THEN
              foundit=find_string(var_name, n_var,                      &
     &                            TRIM(Vname(1,ifield)), varid)
              IF (foundit) THEN
                SELECT CASE (i)
                  CASE (isUice, isUevp)
                    IF (Perfect2D) THEN
                      gtype=var_flag(varid)*u3dvar
                    ELSE
                      gtype=var_flag(varid)*u2dvar
                    END IF
                    IF (Perfect2D) THEN
                      status=nf_fread3d(ng, IDmod, ncname, ncINPid,     &
     &                                  Vname(1,ifield), varid,         &
     &                                  InpRec, gtype, Vsize,           &
     &                                  LBi, UBi, LBj, UBj, 1, 2,       &
     &                                  Fscl, Fmin, Fmax,               &
#   ifdef MASKING
     &                                  GRID(ng) % umask,               &
#   endif
#   ifdef CHECKSUM
     &                                  ICE(ng) % Si(:,:,:,i),          &
     &                                  checksum = Fhash)
#   else
     &                                  ICE(ng) % Si(:,:,:,i))
#   endif
                    ELSE
                      status=nf_fread2d(ng, IDmod, ncname, ncINPid,     &
     &                                  Vname(1,ifield), varid,         &
     &                                  InpRec, gtype, Vsize,           &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  Fscl, Fmin, Fmax,               &
#   ifdef MASKING
     &                                  GRID(ng) % umask,               &
#   endif
#   ifdef CHECKSUM
     &                                  ICE(ng) % Si(:,:,Tindex,i),     &
     &                                  checksum = Fhash)
#   else
     &                                  ICE(ng) % Si(:,:,Tindex,i))
#   endif
                    END IF
                  CASE (isVice, isVevp)
                    IF (Perfect2D) THEN
                      gtype=var_flag(varid)*v3dvar
                    ELSE
                      gtype=var_flag(varid)*v2dvar
                    END IF
                    IF (Perfect2D) THEN
                      status=nf_fread3d(ng, IDmod, ncname, ncINPid,     &
     &                                  Vname(1,ifield), varid,         &
     &                                  InpRec, gtype, Vsize,           &
     &                                  LBi, UBi, LBj, UBj, 1, 2,       &
     &                                  Fscl, Fmin, Fmax,               &
#   ifdef MASKING
     &                                  GRID(ng) % vmask,               &
#   endif
#   ifdef CHECKSUM
     &                                  ICE(ng) % Si(:,:,:,i),          &
     &                                  checksum = Fhash)
#   else
     &                                  ICE(ng) % Si(:,:,:,i))
#   endif
                    ELSE
                      status=nf_fread2d(ng, IDmod, ncname, ncINPid,     &
     &                                  Vname(1,ifield), varid,         &
     &                                  InpRec, gtype, Vsize,           &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  Fscl, Fmin, Fmax,               &
#   ifdef MASKING
     &                                  GRID(ng) % vmask,               &
#   endif
#   ifdef CHECKSUM
     &                                  ICE(ng) % Si(:,:,Tindex,i),     &
     &                                  checksum = Fhash)
#   else
     &                                  ICE(ng) % Si(:,:,Tindex,i))
#   endif
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
#   ifdef MASKING
     &                                  GRID(ng) % rmask,               &
#   endif
#   ifdef CHECKSUM
     &                                  ICE(ng) % Si(:,:,:,i),          &
     &                                  checksum = Fhash)
#   else
     &                                  ICE(ng) % Si(:,:,:,i))
#   endif
                    ELSE
                      status=nf_fread2d(ng, IDmod, ncname, ncINPid,     &
     &                                  Vname(1,ifield), varid,         &
     &                                  InpRec, gtype, Vsize,           &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  Fscl, Fmin, Fmax,               &
#   ifdef MASKING
     &                                  GRID(ng) % rmask,               &
#   endif
#   ifdef CHECKSUM
     &                                  ICE(ng) % Si(:,:,Tindex,i),     &
     &                                  checksum = Fhash)
#   else
     &                                  ICE(ng) % Si(:,:,Tindex,i))
#   endif
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
#   ifdef CHECKSUM
                    WRITE (stdout,70) TRIM(Vname(2,ifield)),            &
     &                                Fmin, Fmax, Fhash
#   else
                    WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
#   endif
                  END IF
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
!
!  Read iceshelf model state variables.
!
        DO i=1,nIceF
          IF (iFice(i).gt.0) THEN
            ifield=iFice(i)
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
     &                            ICE(ng) % Fi(:,:,i),                  &
     &                            checksum = Fhash)
#   else
     &                            ICE(ng) % Fi(:,:,i))
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
!
!  Read in integrated tracer (m Tunits).
!
        itrc=0
        DO i=NAT+1,NT(ng)
          ifield=idStrc(i)
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
     &                            ICE(ng) % Ti(:,:,:,itrc),             &
     &                            checksum = Fhash)
#   else
     &                            ICE(ng) % Ti(:,:,:,itrc))
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
     &                            ICE(ng) % Ti(:,:,Tindex,itrc),        &
     &                            checksum = Fhash)
#   else
     &                            ICE(ng) % Ti(:,:,Tindex,itrc))
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
