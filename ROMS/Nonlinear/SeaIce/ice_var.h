/*
** git $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2025 The ROMS Group                             **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.md                                              **
************************************************************************
**                                                                    **
**  Assigns metadata indices for the sea ice model (Budgell, 2005;    **
**  Durski and Kurapov, 2019, 2020) variables that are used in input  **
**  and output NetCDF files.  The metadata information is read from   **
**  "varinfo.yaml".                                                   **
**                                                                    **
**  This file is included in file "mod_ncparam.F", routine            **
**  "initialize_ncparm".                                              **
**                                                                    **
************************************************************************
*/

/*
**  Ice state variables.
*/

          CASE ('idsnow')
            idsnow=varid
          CASE ('idAice')
            idAice=varid
            iSice(isAice)=idAice
          CASE ('idHage')
            idHage=varid
            iSice(isHage)=idHage
          CASE ('idHice')
            idHice=varid
            iSice(isHice)=idHice
          CASE ('idHmel')
            idHmel=varid
            iSice(isHmel)=idHmel
          CASE ('idHsno')
            idHsno=varid
            iSice(isHsno)=idHsno
          CASE ('idIage')
            idIage=varid
            iSice(isIage)=idIage
          CASE ('idISxx')
            idISxx=varid
            iSice(isISxx)=idISxx
          CASE ('idISxy')
            idISxy=varid
            iSice(isISxy)=idISxy
          CASE ('idISyy')
            idISyy=varid
            iSice(isISyy)=idISyy
          CASE ('idTice')
            idTice=varid
            iSice(isTice)=idTice
          CASE ('idUice')
            idUice=varid
            iSice(isUice)=idUice
          CASE ('idVice')
            idVice=varid
            iSice(isVice)=idVice
#if defined ICE_BIO
          CASE ('idIphy')
            isIphy=varid
            iSice(isIphy)=idIphy
          CASE ('idINO3')
            isINO3=varid
            iSice(isINO3)=idINO3
          CASE ('idINH4')
            isINH4=varid
            iSice(isINH4)=idINH4
          CASE ('idIlog')
            isIlog=varid
            iSice(isIlog)=idIlog
#endif
          CASE ('idAiCL')
            idAiCL=varid
          CASE ('idHiCL')
            idHiCL=varid
          CASE ('idIOfv')
            idIOfv=varid
            iFice(icIOfv)=idIOfv
          CASE ('idIOmf')
            idIOmf=varid
            iFice(icIOmf)=idIOmf
          CASE ('idIOmt')
            idIOmt=varid
            iFice(icIOmt)=idIOmt
          CASE ('idIsst')
            idIsst=varid
            iFice(icIsst)=idIsst
          CASE ('idS0mk')
            idS0mk=varid
            iFice(icS0mk)=idS0mk
          CASE ('idT0mk')
            idT0mk=varid
            iFice(icT0mk)=idT0mk
          CASE ('idUiCL')
            idUiCL=varid
          CASE ('idUiER')
            idUiER=varid
          CASE ('idViCL')
            idViCL=varid
          CASE ('idViNR')
            idViNR=varid
          CASE ('idWdiv')
            idWdiv=varid
            iFice(icWdiv)=idWdiv
          CASE ('idW_ai')
            idW_ai=varid
            iFice(icW_ai)=idW_ai
          CASE ('idW_ao')
            idW_ao=varid
            iFice(icW_ao)=idW_ao
          CASE ('idW_fr')
            idW_fr=varid
            iFice(icW_fr)=idW_fr
          CASE ('idW_io')
            idW_io=varid
            iFice(icW_io)=idW_io
          CASE ('idW_ro')
            idW_ro=varid
            iFice(icW_ro)=idW_ro
#if defined ICE_BULK_FLUXES
          CASE ('idAlbi')
            idAlbi=varid
            iFice(icAlbi)=idAlbi
          CASE ('idAlbo')
            idAlbo=varid
            iFice(icAlbo)=idAlbo
          CASE ('idFrsw')
            idFrsw=varid
            iFice(icFrsw)=idFrsw
          CASE ('idSics')
            idSics=varid
            iFice(icSics)=idSics
          CASE ('idSicb')
            idSicb=varid
            iFice(icSicb)=idSicb
          CASE ('idSocn')
            idSocn=varid
            iFice(icSocn)=idSocn
          CASE ('idRfac')
            idRfac=varid
            iFice(icRfac)=idRfac
          CASE ('idSnof')
            idSnof=varid
            iFice(icSnof)=idSnof
#endif
          CASE ('idStrc(inert(i))')
            load=.TRUE.
/*
**  Ice model open boundary conditions.
*/

          CASE ('idAice(iwest)')
            iceOBC(iwest,isAice)=varid
          CASE ('idAice(ieast)')
            iceOBC(ieast,isAice)=varid
          CASE ('idAice(isouth)')
            iceOBC(isouth,isAice)=varid
          CASE ('idAice(inorth)')
            iceOBC(inorth,isAice)=varid

          CASE ('idHice(iwest)')
            iceOBC(iwest,isHice)=varid
          CASE ('idHice(ieast)')
            iceOBC(ieast,isHice)=varid
          CASE ('idHice(isouth)')
            iceOBC(isouth,isHice)=varid
          CASE ('idHice(inorth)')
            iceOBC(inorth,isHice)=varid

          CASE ('idHmel(iwest)')
            iceOBC(iwest,isHmel)=varid
          CASE ('idHmel(ieast)')
            iceOBC(ieast,isHmel)=varid
          CASE ('idHmel(isouth)')
            iceOBC(isouth,isHmel)=varid
          CASE ('idHmel(inorth)')
            iceOBC(inorth,isHmel)=varid

          CASE ('idHsno(iwest)')
            iceOBC(iwest,isHsno)=varid
          CASE ('idHsno(ieast)')
            iceOBC(ieast,isHsno)=varid
          CASE ('idHsno(isouth)')
            iceOBC(isouth,isHsno)=varid
          CASE ('idHsno(inorth)')
            iceOBC(inorth,isHsno)=varid

          CASE ('idIage(iwest)')
            iceOBC(iwest,isIage)=varid
          CASE ('idIage(ieast)')
            iceOBC(ieast,isIage)=varid
          CASE ('idIage(isouth)')
            iceOBC(isouth,isIage)=varid
          CASE ('idIage(inorth)')
            iceOBC(inorth,isIage)=varid

          CASE ('idISxx(iwest)')
            iceOBC(iwest,isISxx)=varid
          CASE ('idISxx(ieast)')
            iceOBC(ieast,isISxx)=varid
          CASE ('idISxx(isouth)')
            iceOBC(isouth,isISxx)=varid
          CASE ('idISxx(inorth)')
            iceOBC(inorth,isISxx)=varid

          CASE ('idISxy(iwest)')
            iceOBC(iwest,isISxy)=varid
          CASE ('idISxy(ieast)')
            iceOBC(ieast,isISxy)=varid
          CASE ('idISxy(isouth)')
            iceOBC(isouth,isISxy)=varid
          CASE ('idISxy(inorth)')
            iceOBC(inorth,isISxy)=varid

          CASE ('idISyy(iwest)')
            iceOBC(iwest,isISyy)=varid
          CASE ('idISyy(ieast)')
            iceOBC(ieast,isISyy)=varid
          CASE ('idISyy(isouth)')
            iceOBC(isouth,isISyy)=varid
          CASE ('idISyy(inorth)')
            iceOBC(inorth,isISyy)=varid

          CASE ('idTice(iwest)')
            iceOBC(iwest,isTice)=varid
          CASE ('idTice(ieast)')
            iceOBC(ieast,isTice)=varid
          CASE ('idTice(isouth)')
            iceOBC(isouth,isTice)=varid
          CASE ('idTice(inorth)')
            iceOBC(inorth,isTice)=varid

          CASE ('idUice(iwest)')
            iceOBC(iwest,isUice)=varid
          CASE ('idUice(ieast)')
            iceOBC(ieast,isUice)=varid
          CASE ('idUice(isouth)')
            iceOBC(isouth,isUice)=varid
          CASE ('idUice(inorth)')
            iceOBC(inorth,isUice)=varid

          CASE ('idVice(iwest)')
            iceOBC(iwest,isVice)=varid
          CASE ('idVice(ieast)')
            iceOBC(ieast,isVice)=varid
          CASE ('idVice(isouth)')
            iceOBC(isouth,isVice)=varid
          CASE ('idVice(inorth)')
            iceOBC(inorth,isVice)=varid
