#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Module containing data corresponding to the NLSCATM namelist

Module Rcf_CntlAtm_Mod

! Description:
!   Module to contain the CNTLATM deck
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.2   17/11/00   Tracer, Vegetation and MOSES II logicals
!                    P.Selwood and M. Best.
!   5.3   22/10/01   Coastal tiling and fixed leads temperature
!                    added.                        Nic Gedney.
!   5.3   19/06/01   Add Support for Tropopause-based Ozone.  Dave Tan
!   5.3   29/08/01   Sulphate and sea-salt logicals split into separate
!                    versions for different processes.  A. Jones
!   5.3   29/11/01   Add logical switches for spatial degradation of
!                    radiation, and radiative effects of murk. S.Cusack
!   5.3   14/11/01   put problem_number in namelist  A. Malcolm
!
!   5.3   29/11/01   Add changes for methane oxidation. David Jackson
!   5.4   15/08/02   Use only the UM version

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Implicit None

#include "cntlatm.h"

End Module Rcf_CntlAtm_Mod
#endif
