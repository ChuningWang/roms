#include "cppdefs.h"
#ifdef ICESHELF_DYNAM
!!
!!git $Id$
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2024 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.md                                               !
!!======================================================================
!!                                                                     !
!!  This routine reads and reports iceshelf model input parameters.    !
!!                                                                     !
!!======================================================================
!!
# include <isf_inp.h>
#else
      SUBROUTINE read_IsfPar
      END SUBROUTINE read_IsfPar
#endif
