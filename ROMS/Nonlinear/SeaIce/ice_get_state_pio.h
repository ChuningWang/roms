!
!  Read sea ice model state variables.
!
        DO i=1,nIceS
          IF (iSice(i).gt.0) THEN
            ifield=iSice(i)
            IF (get_var(ifield)) THEN
              foundit=find_string(var_name, n_var,                      &
     &                            TRIM(Vname(1,ifield)), vindex)
              IF (foundit) THEN
                my_pioVar%vd=var_desc(vindex)
                SELECT CASE (i)
                  CASE (isUice, isUevp)
                    IF (Perfect2D) THEN
                      my_pioVar%gtype=u3dvar
                    ELSE
                      my_pioVar%gtype=u2dvar
                    END IF
                    IF (KIND(ICE(ng)%Si).eq.8) THEN
                      my_pioVar%dkind=PIO_double
                      IF (Perfect2D) THEN
                        ioDesc => ioDesc_dp_u3dvar(ng)
                      ELSE
                        ioDesc => ioDesc_dp_u2dvar(ng)
                      END IF
                    ELSE
                      my_pioVar%dkind=PIO_real
                      IF (Perfect2D) THEN
                        ioDesc => ioDesc_sp_u3dvar(ng)
                      ELSE
                        ioDesc => ioDesc_sp_u2dvar(ng)
                      END IF
                    END IF
!
                    IF (Perfect2D) THEN
                      status=nf_fread3d(ng, IDmod, ncname, pioFile,     &
     &                                  Vname(1,ifield), my_pioVar,     &
     &                                  InpRec, ioDesc, Vsize,          &
     &                                  LBi, UBi, LBj, UBj, 1, 2,       &
     &                                  Fscl, Fmin, Fmax,               &
#    ifdef MASKING
     &                                  GRID(ng) % umask,               &
#    endif
#    ifdef CHECKSUM
     &                                  ICE(ng) % Si(:,:,:,i),          &
     &                                  checksum = Fhash)
#    else
     &                                  ICE(ng) % Si(:,:,:,i))
#    endif
                    ELSE
                      status=nf_fread2d(ng, IDmod, ncname, pioFile,     &
     &                                  Vname(1,ifield), my_pioVar,     &
     &                                  InpRec, ioDesc, Vsize,          &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  Fscl, Fmin, Fmax,               &
#    ifdef MASKING
     &                                  GRID(ng) % umask,               &
#    endif
#    ifdef CHECKSUM
     &                                  ICE(ng) % Si(:,:,Tindex,i),     &
     &                                  checksum = Fhash)
#    else
     &                                  ICE(ng) % Si(:,:,Tindex,i))
#    endif
                    END IF
                  CASE (isVice, isVevp)
                    IF (Perfect2D) THEN
                      my_pioVar%gtype=v3dvar
                    ELSE
                      my_pioVar%gtype=v2dvar
                    END IF
                    IF (KIND(ICE(ng)%Si).eq.8) THEN
                      my_pioVar%dkind=PIO_double
                      IF (Perfect2D) THEN
                        ioDesc => ioDesc_dp_v3dvar(ng)
                      ELSE
                        ioDesc => ioDesc_dp_v2dvar(ng)
                      END IF
                    ELSE
                      my_pioVar%dkind=PIO_real
                      IF (Perfect2D) THEN
                        ioDesc => ioDesc_sp_v3dvar(ng)
                      ELSE
                        ioDesc => ioDesc_sp_v2dvar(ng)
                      END IF
                    END IF
!
                    IF (Perfect2D) THEN
                      status=nf_fread3d(ng, IDmod, ncname, pioFile,     &
     &                                  Vname(1,ifield), my_pioVar,     &
     &                                  InpRec, ioDesc, Vsize,          &
     &                                  LBi, UBi, LBj, UBj, 1, 2,       &
     &                                  Fscl, Fmin, Fmax,               &
#    ifdef MASKING
     &                                  GRID(ng) % vmask,               &
#    endif
#    ifdef CHECKSUM
     &                                  ICE(ng) % Si(:,:,:,i),          &
     &                                  checksum = Fhash)
#    else
     &                                  ICE(ng) % Si(:,:,:,i))
#    endif
                    ELSE
                      status=nf_fread2d(ng, IDmod, ncname, pioFile,     &
     &                                  Vname(1,ifield), my_pioVar,     &
     &                                  InpRec, ioDesc, Vsize,          &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  Fscl, Fmin, Fmax,               &
#    ifdef MASKING
     &                                  GRID(ng) % vmask,               &
#    endif
#    ifdef CHECKSUM
     &                                  ICE(ng) % Si(:,:,Tindex,i),     &
     &                                  checksum = Fhash)
#    else
     &                                  ICE(ng) % Si(:,:,Tindex,i))
#    endif
                    END IF
                  CASE DEFAULT
                    IF (Perfect2D) THEN
                      my_pioVar%gtype=r3dvar
                    ELSE
                      my_pioVar%gtype=r2dvar
                    END IF
                    IF (KIND(ICE(ng)%Si).eq.8) THEN
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
#    ifdef MASKING
     &                                  GRID(ng) % rmask,               &
#    endif
#    ifdef CHECKSUM
     &                                  ICE(ng) % Si(:,:,:,i),          &
     &                                  checksum = Fhash)
#    else
     &                                  ICE(ng) % Si(:,:,:,i))
#    endif
                    ELSE
                      status=nf_fread2d(ng, IDmod, ncname, pioFile,     &
     &                                  Vname(1,ifield), my_pioVar,     &
     &                                  InpRec, ioDesc, Vsize,          &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  Fscl, Fmin, Fmax,               &
#    ifdef MASKING
     &                                  GRID(ng) % rmask,               &
#    endif
#    ifdef CHECKSUM
     &                                  ICE(ng) % Si(:,:,Tindex,i),     &
     &                                  checksum = Fhash)
#    else
     &                                  ICE(ng) % Si(:,:,Tindex,i))
#    endif
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
!  Read sea ice model state variables.
!
        DO i=1,nIceF
          IF (iFice(i).gt.0) THEN
            ifield=iFice(i)
            IF (get_var(ifield).and.Perfect2D) THEN
              foundit=find_string(var_name, n_var,                      &
     &                            TRIM(Vname(1,ifield)), vindex)
              IF (foundit) THEN
                my_pioVar%vd=var_desc(vindex)
                my_pioVar%gtype=r2dvar
                IF (KIND(ICE(ng)%Si).eq.8) THEN
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
     &                            ICE(ng) % Fi(:,:,i),                  &
     &                            checksum = Fhash)
#    else
     &                            ICE(ng) % Fi(:,:,i))
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
!  Read in integrated tracer (m Tunits).
!
        itrc=0
        DO i=NAT+1,NT(ng)
          ifield=idStrc(i)
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
              IF (KIND(ICE(ng)%Ti).eq.8) THEN
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
     &                            ICE(ng) % Ti(:,:,:,itrc),             &
     &                            checksum = Fhash)
#    else
     &                            ICE(ng) % Ti(:,:,:,itrc))
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
     &                            ICE(ng) % Ti(:,:,Tindex,itrc),        &
     &                            checksum = Fhash)
#    else
     &                            ICE(ng) % Ti(:,:,Tindex,itrc))
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
