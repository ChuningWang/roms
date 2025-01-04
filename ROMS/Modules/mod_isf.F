#include "cppdefs.h"
#ifdef ICESHELF_DYNAM
/*
** git $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2024 The ROMS/TOMS Group          Chuning Wang  **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.md                                              **
************************************************************************
**                                                                    **
** This module declares iceshelf model variables and parameteres.     **
**                                                                    **
** The current design allows the user to have a lot of latitude for   **
** customizing or adding any particular sea-ice model.                **
**                                                                    **
** Each ice model is composed of these files:                         **
**                                                                    **
**   (1) Internal model parameters declaration:                       **
**                                                                    **
**         isf_mod.h                                                  **
**                                                                    **
**   (2) Model parameters standard input script:                      **
**                                                                    **
**         isf.in                                                     **
**                                                                    **
**   (3) Code to read input model parameters:                         **
**                                                                    **
**         isf_inp.h                                                  **
**                                                                    **
**   (4) Code to assign indices to model variables during the         **
**       reading of metadata information from "varinfo.yaml":         **
**                                                                    **
**         isf_var.h                                                  **
**         isf_var_tracer.h                                           **
**                                                                    **
** Note that all the *.h files are located in ROMS/Nonlinear/IceShelf **
** and included within  <...>  to allow the user to customize any of  **
** them in the project directory  while keeping the distributed code  **
** intact (check the build script for details).                       **
**                                                                    **
************************************************************************
*/
# include <isf_mod.h>
#else
      MODULE mod_isf
      END MODULE mod_isf
#endif
