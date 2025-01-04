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
                  CASE (isIuis)
                    my_pioVar%gtype=u2dvar
                    IF (KIND(ISF(ng)%Si).eq.8) THEN
                      my_pioVar%dkind=PIO_double
                      ioDesc => ioDesc_dp_u2dvar(ng)
                    ELSE
                      my_pioVar%dkind=PIO_real
                      ioDesc => ioDesc_sp_u2dvar(ng)
                    END IF
!
                    status=nf_fread2d(ng, IDmod, ncname, pioFile,       &
     &                                Vname(1,ifield), my_pioVar,       &
     &                                InpRec, ioDesc, Vsize,            &
     &                                LBi, UBi, LBj, UBj,               &
     &                                Fscl, Fmin, Fmax,                 &
#    ifdef MASKING
     &                                GRID(ng) % umask,                 &
#    endif
#    ifdef CHECKSUM
     &                                ISF(ng) % Si(:,:,Tindex,i),       &
     &                                checksum = Fhash)
#    else
     &                                ISF(ng) % Si(:,:,Tindex,i))
#    endif
                  CASE (isIvis)
                    my_pioVar%gtype=v2dvar
                    IF (KIND(ISF(ng)%Si).eq.8) THEN
                      my_pioVar%dkind=PIO_double
                      ioDesc => ioDesc_dp_v2dvar(ng)
                    ELSE
                      my_pioVar%dkind=PIO_real
                      ioDesc => ioDesc_sp_v2dvar(ng)
                    END IF
!
                    status=nf_fread2d(ng, IDmod, ncname, pioFile,       &
     &                                Vname(1,ifield), my_pioVar,       &
     &                                InpRec, ioDesc, Vsize,            &
     &                                LBi, UBi, LBj, UBj,               &
     &                                Fscl, Fmin, Fmax,                 &
#    ifdef MASKING
     &                                GRID(ng) % vmask,                 &
#    endif
#    ifdef CHECKSUM
     &                                ISF(ng) % Si(:,:,Tindex,i),       &
     &                                checksum = Fhash)
#    else
     &                                ISF(ng) % Si(:,:,Tindex,i))
#    endif
                  CASE DEFAULT
                    my_pioVar%gtype=r2dvar
                    IF (KIND(ISF(ng)%Si).eq.8) THEN
                      my_pioVar%dkind=PIO_double
                      ioDesc => ioDesc_dp_r2dvar(ng)
                    ELSE
                      my_pioVar%dkind=PIO_real
                      ioDesc => ioDesc_sp_r2dvar(ng)
                    END IF
!
                    status=nf_fread2d(ng, IDmod, ncname, pioFile,       &
     &                                Vname(1,ifield), my_pioVar,       &
     &                                InpRec, ioDesc, Vsize,            &
     &                                LBi, UBi, LBj, UBj,               &
     &                                Fscl, Fmin, Fmax,                 &
#    ifdef MASKING
     &                                GRID(ng) % rmask,                 &
#    endif
#    ifdef CHECKSUM
     &                                ISF(ng) % Si(:,:,Tindex,i),       &
     &                                checksum = Fhash)
#    else
     &                                ISF(ng) % Si(:,:,Tindex,i))
#    endif
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
#    ifdef ICESHELF_TRACER
!
!  Read in iceshelf volume (m).
!
        ifield=idIvol
        IF (get_var(ifield)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,ifield)),   &
     &                        vindex)
          IF (foundit) THEN
            my_pioVar%vd=var_desc(vindex)
            my_pioVar%gtype=r2dvar
            IF (KIND(ISF(ng)%vol).eq.8) THEN
              my_pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              my_pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
            status=nf_fread2d(ng, IDmod, ncname, pioFile,               &
     &                        Vname(1,ifield), my_pioVar,               &
     &                        InpRec, ioDesc, Vsize,                    &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#     ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#     endif
#     ifdef CHECKSUM
     &                        ISF(ng) % vol,                            &
     &                        checksum = Fhash)
#     else
     &                        ISF(ng) % vol)
#     endif
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
#     ifdef CHECKSUM
                WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax,    &
     &                            Fhash
#     else
                WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
#     endif

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
!
!  Read in integrated tracer (m Tunits).
!
        DO itrc=1,NPT
          ifield=idItrc(itrc)
          IF (get_var(ifield)) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,ifield)), vindex)
            IF (foundit) THEN
              scale=1.0_dp
              my_pioVar%vd=var_desc(vindex)
              my_pioVar%gtype=r2dvar
              IF (KIND(ISF(ng)%trc).eq.8) THEN
                my_pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                my_pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, IDmod, ncname, pioFile,             &
     &                          Vname(1,ifield), my_pioVar,             &
     &                          InpRec, ioDesc, Vsize,                  &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          scale, Fmin, Fmax,                      &
#     ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#     endif
#     ifdef CHECKSUM
     &                          ISF(ng)% trc(:,:,itrc),                 &
     &                          checksum = Fhash)
#     else
     &                          ISF(ng)% trc(:,:,itrc))
#     endif
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
#     ifdef CHECKSUM
                  WRITE (stdout,70) TRIM(Vname(2,ifield)),              &
     &                              Fmin, Fmax, Fhash
#     else
                  WRITE (stdout,70) TRIM(Vname(2,ifield)),              &
     &                              Fmin, Fmax
#     endif
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
#    endif
