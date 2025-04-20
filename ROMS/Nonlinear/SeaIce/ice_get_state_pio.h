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
                  CASE (isUice)
                    my_pioVar%gtype=u2dvar
                    IF (KIND(ICE(ng)%Si).eq.8) THEN
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
     &                                ICE(ng) % Si(:,:,Tindex,i),       &
     &                                checksum = Fhash)
#    else
     &                                ICE(ng) % Si(:,:,Tindex,i))
#    endif
                  CASE (isVice)
                    my_pioVar%gtype=v2dvar
                    IF (KIND(ICE(ng)%Si).eq.8) THEN
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
     &                                ICE(ng) % Si(:,:,Tindex,i),       &
     &                                checksum = Fhash)
#    else
     &                                ICE(ng) % Si(:,:,Tindex,i))
#    endif
                  CASE DEFAULT
                    my_pioVar%gtype=r2dvar
                    IF (KIND(ICE(ng)%Si).eq.8) THEN
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
     &                                ICE(ng) % Si(:,:,Tindex,i),       &
     &                                checksum = Fhash)
#    else
     &                                ICE(ng) % Si(:,:,Tindex,i))
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
              my_pioVar%gtype=r2dvar
              IF (KIND(ICE(ng)%Ti).eq.8) THEN
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
#    ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#    endif
#    ifdef CHECKSUM
     &                          ICE(ng) % Ti(:,:,Tindex,itrc),          &
     &                          checksum = Fhash)
#    else
     &                          ICE(ng) % Ti(:,:,Tindex,itrc))
#    endif
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
