/*
** git $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2024 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.md                                              **
************************************************************************
**                                                                    **
**  Assigns metadata indices for the iceshelf model variables that    **
**  are used in input and output NetCDF files.  The metadata          **
**  information is read from "varinfo.yaml".                          **
**                                                                    **
**  This file is included in file "mod_ncparam.F", routine            **
**  "initialize_ncparm".                                              **
**                                                                    **
************************************************************************
*/

/*
**  Iceshelf state variables.
*/

#ifdef ICESHELF_TRACER
            CASE ('idItrc(inert(i))')
              IF (NPT.gt.0) THEN
                varid=varid-1
# ifdef AGE_MEAN
                ic=0
                DO i=1,NPT,2
                  varid=varid+1
                  ic=ic+1
                  idItrc(i)=varid
                  DO ng=1,Ngrids
                    Fscale(varid,ng)=scale
                    Iinfo(1,varid,ng)=gtype
                  END DO
                  WRITE (Vname(1,varid),'(a,i2.2)')                     &
     &                  TRIM(ADJUSTL(Vinfo(1))), ic
                  WRITE (Vname(2,varid),'(a,a,i2.2)')                   &
     &                  TRIM(ADJUSTL(Vinfo(2))), ', type ', ic
                  WRITE (Vname(3,varid),'(a)')                          &
     &                  TRIM(ADJUSTL(Vinfo(3)))
                  WRITE (Vname(4,varid),'(a,i2.2)')                     &
     &                  TRIM(ADJUSTL(Vinfo(4))), ic
                  WRITE (Vname(5,varid),'(a)')                          &
     &                  TRIM(ADJUSTL(Vinfo(5)))
                  WRITE (Vname(6,varid),'(a,i2.2)')                     &
     &                  TRIM(ADJUSTL(Vinfo(6))), ic
                END DO
!
                ic=0
                DO i=2,NPT,2
                  varid=varid+1
                  ic=ic+1
                  idItrc(i)=varid
                  DO ng=1,Ngrids
                    Fscale(varid,ng)=scale
                    Iinfo(1,varid,ng)=gtype
                  END DO
                  WRITE (Vname(1,varid),'(a,i2.2,a)')                   &
     &                  TRIM(ADJUSTL(Vinfo(1))), ic, '_age'
                  WRITE (Vname(2,varid),'(a,i2.2)')                     &
     &                  'age concentration, type ', ic
                  WRITE (Vname(3,varid),'(a)')                          &
     &                  'second kilogram meter-3'
                  WRITE (Vname(4,varid),'(a,i2.2,a)')                   &
     &                  TRIM(ADJUSTL(Vinfo(4))), ic, '_age'
                  WRITE (Vname(5,varid),'(a)')                          &
     &                  TRIM(ADJUSTL(Vinfo(5)))
                  WRITE (Vname(6,varid),'(a,i2.2,a)')                   &
     &                  TRIM(ADJUSTL(Vinfo(6))), ic, '_age'
                END DO
# else
                DO i=1,NPT
                  varid=varid+1
                  idItrc(i)=varid
                  DO ng=1,Ngrids
                    Fscale(varid,ng)=scale
                    Iinfo(1,varid,ng)=gtype
                  END DO
                  WRITE (Vname(1,varid),'(a,i2.2)')                     &
     &                  TRIM(ADJUSTL(Vinfo(1))), i
                  WRITE (Vname(2,varid),'(a,a,i2.2)')                   &
     &                  TRIM(ADJUSTL(Vinfo(2))), ', type ', i
                  WRITE (Vname(3,varid),'(a)')                          &
     &                  TRIM(ADJUSTL(Vinfo(3)))
                  WRITE (Vname(4,varid),'(a,i2.2)')                     &
     &                  TRIM(ADJUSTL(Vinfo(4))), i
                  WRITE (Vname(5,varid),'(a)')                          &
     &                  TRIM(ADJUSTL(Vinfo(5)))
                  WRITE (Vname(6,varid),'(a,i2.2)')                     &
     &                  TRIM(ADJUSTL(Vinfo(6))), i
                END DO
# endif
              END IF
#endif
