# git $Id$
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2025 The ROMS Group                                :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Nonlinear/IceShelf"

list( APPEND _files
      ROMS/Nonlinear/IceShelf/iceshelf.F
      ROMS/Nonlinear/IceShelf/isf_stress.F
      ROMS/Nonlinear/IceShelf/isf_masks.F
      ROMS/Nonlinear/IceShelf/isf_mms_source.F
      ROMS/Nonlinear/IceShelf/isf_mom.F
      ROMS/Nonlinear/IceShelf/isf_output.F
      ROMS/Nonlinear/IceShelf/isf_sbl.F
      ROMS/Nonlinear/IceShelf/isf_set_avg.F
      ROMS/Nonlinear/IceShelf/isf_advection.F
      ROMS/Nonlinear/IceShelf/isf_calving.F
      ROMS/Nonlinear/IceShelf/isf_bc2d.F
      ROMS/Nonlinear/IceShelf/isf_ini_fields.F
      ROMS/Nonlinear/IceShelf/isf_pibc.F
      ROMS/Nonlinear/IceShelf/isf_tibc.F
      ROMS/Nonlinear/IceShelf/isf_vbc.F
)

set ( ROMS_Nonlinear_IceShelf_files
      ${_files}
      PARENT_SCOPE
)
