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
          CASE ('isIhis(iwest)')
            isfOBC(iwest,isIhis)=varid
          CASE ('isIhis(ieast)')
            isfOBC(ieast,isIhis)=varid
          CASE ('isIhis(isouth)')
            isfOBC(isouth,isIhis)=varid
          CASE ('isIhis(inorth)')
            isfOBC(inorth,isIhis)=varid
          CASE ('isIuis(iwest)')
            isfOBC(iwest,isIuis)=varid
          CASE ('isIuis(ieast)')
            isfOBC(ieast,isIuis)=varid
          CASE ('isIuis(isouth)')
            isfOBC(isouth,isIuis)=varid
          CASE ('isIuis(inorth)')
            isfOBC(inorth,isIuis)=varid
          CASE ('isIvis(iwest)')
            isfOBC(iwest,isIvis)=varid
          CASE ('isIvis(ieast)')
            isfOBC(ieast,isIvis)=varid
          CASE ('isIvis(isouth)')
            isfOBC(isouth,isIvis)=varid
          CASE ('isIvis(inorth)')
            isfOBC(inorth,isIvis)=varid
          CASE ('isIeis(iwest)')
            isfOBC(iwest,isIeis)=varid
          CASE ('isIeis(ieast)')
            isfOBC(ieast,isIeis)=varid
          CASE ('isIeis(isouth)')
            isfOBC(isouth,isIeis)=varid
          CASE ('isIeis(inorth)')
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
            iFisf(icItis)=idItis
          CASE ('idIzsf')
            idIzsf=varid
            iFisf(icIzsf)=idIzsf
          CASE ('idIzbd')
            idIzbd=varid
            iFisf(icIzbd)=idIzbd
          CASE ('idInuh')
            idInuh=varid
            iFisf(icInuh)=idInuh
          CASE ('idItdx')
            idItdx=varid
            iFisf(icItdx)=idItdx
          CASE ('idItdy')
            idItdy=varid
            iFisf(icItdy)=idItdy
          CASE ('idIbta')
            idIbta=varid
            iFisf(icIbta)=idIbta
#endif
#ifdef ICESHELF_THERMO
          CASE ('idImel')
            idImel=varid
            iFisf(icImel)=idImel
          CASE ('idIgat')
            idIgat=varid
            iFisf(icIgat)=idIgat
          CASE ('idIgas')
            idIgas=varid
            iFisf(icIgas)=idIgas
          CASE ('idItst')
            idItst=varid
            iFisf(icItst)=idItst
          CASE ('idIust')
            idIust=varid
            iFisf(icIust)=idIust
          CASE ('idItbo')
            idItbo=varid
            iFisf(icItbo)=idItbo
          CASE ('idIsbo')
            idIsbo=varid
            iFisf(icIsbo)=idIsbo
          CASE ('idItba')
            idItba=varid
            iFisf(icItba)=idItba
          CASE ('idIsba')
            idIsba=varid
            iFisf(icIsba)=idIsba
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
