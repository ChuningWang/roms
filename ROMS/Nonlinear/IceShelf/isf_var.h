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

#ifdef ICESHEET
          CASE ('idIhis(iwest)')
            isfOBC(iwest,isIhis)=varid
          CASE ('idIhis(ieast)')
            isfOBC(ieast,isIhis)=varid
          CASE ('idIhis(isouth)')
            isfOBC(isouth,isIhis)=varid
          CASE ('idIhis(inorth)')
            isfOBC(inorth,isIhis)=varid
          CASE ('idIuis(iwest)')
            isfOBC(iwest,isIuis)=varid
          CASE ('idIuis(ieast)')
            isfOBC(ieast,isIuis)=varid
          CASE ('idIuis(isouth)')
            isfOBC(isouth,isIuis)=varid
          CASE ('idIuis(inorth)')
            isfOBC(inorth,isIuis)=varid
          CASE ('idIvis(iwest)')
            isfOBC(iwest,isIvis)=varid
          CASE ('idIvis(ieast)')
            isfOBC(ieast,isIvis)=varid
          CASE ('idIvis(isouth)')
            isfOBC(isouth,isIvis)=varid
          CASE ('idIvis(inorth)')
            isfOBC(inorth,isIvis)=varid
!
!  The isIeis boundary file carries the ice TEMPERATURE (not the
!  enthalpy); the enthalpy boundary value is computed as h*T in
!  isf_tibc.
!
          CASE ('idIeis(iwest)')
            isfOBC(iwest,isIeis)=varid
          CASE ('idIeis(ieast)')
            isfOBC(ieast,isIeis)=varid
          CASE ('idIeis(isouth)')
            isfOBC(isouth,isIeis)=varid
          CASE ('idIeis(inorth)')
            isfOBC(inorth,isIeis)=varid
          CASE ('idIhis')
            idIhis=varid
            iSisf(isIhis)=idIhis
          CASE ('idIuis')
            idIuis=varid
            iSisf(isIuis)=idIuis
          CASE ('idIvis')
            idIvis=varid
            iSisf(isIvis)=idIvis
          CASE ('idIuer')
            idIuer=varid
          CASE ('idIvnr')
            idIvnr=varid
          CASE ('idIeis')
            idIeis=varid
            iSisf(isIeis)=idIeis
          CASE ('idItis')
            idItis=varid
            iFisf(ifItis)=idItis
          CASE ('idIzsf')
            idIzsf=varid
            iFisf(ifIzsf)=idIzsf
          CASE ('idIzbd')
            idIzbd=varid
            iFisf(ifIzbd)=idIzbd
          CASE ('idInuh')
            idInuh=varid
            iFisf(ifInuh)=idInuh
          CASE ('idItdx')
            idItdx=varid
            iFisf(ifItdx)=idItdx
          CASE ('idItdy')
            idItdy=varid
            iFisf(ifItdy)=idItdy
          CASE ('idIbta')
            idIbta=varid
            iFisf(ifIbta)=idIbta
          CASE ('idCzta')
            idCzta=varid
            iCisf(icIzta)=idCzta
#endif
#ifdef ICESHELF_THERMO
          CASE ('idImel')
            idImel=varid
            iFisf(ifImel)=idImel
          CASE ('idIgat')
            idIgat=varid
            iFisf(ifIgat)=idIgat
          CASE ('idIgas')
            idIgas=varid
            iFisf(ifIgas)=idIgas
          CASE ('idItst')
            idItst=varid
            iFisf(ifItst)=idItst
          CASE ('idIust')
            idIust=varid
            iFisf(ifIust)=idIust
          CASE ('idItbo')
            idItbo=varid
            iFisf(ifItbo)=idItbo
          CASE ('idIsbo')
            idIsbo=varid
            iFisf(ifIsbo)=idIsbo
          CASE ('idItba')
            idItba=varid
            iFisf(ifItba)=idItba
          CASE ('idIsba')
            idIsba=varid
            iFisf(ifIsba)=idIsba
          CASE ('idIsbl')
            idIsbl=varid
            iOisf(ioHsbl)=idIsbl
          CASE ('idIvol')
            idIvol=varid
          CASE ('idItrc(inert(i))')
            load=.TRUE.
          CASE ('idIamb(inert(i))')
            load=.TRUE.
#endif
