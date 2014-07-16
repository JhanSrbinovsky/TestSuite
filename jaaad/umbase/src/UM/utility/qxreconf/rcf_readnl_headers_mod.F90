#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  Read the HEADERS namelist

Module rcf_readnl_headers_mod

!  Subroutine Rcf_Readnl_Headers - Read the HEADERS namelist
!
! Description:
!   Read the HEADERS namelist for header overrides.
!
! Method:
!   Variables read into Rcf_Headers_Mod module.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine rcf_readnl_headers( nft )

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Oper

Use Rcf_Parvars_Mod, Only : &
    mype

Use rcf_headers_mod, Only : &
    Fixhd,                  &
    Inthd,                  &
    Relhd,                  &
    headers

Implicit None

! Arguments
Integer, Intent (In)    :: nft

! Comdecks
#include "c_mdi.h"

! Set defaults = MDI
FixHd(:) = IMDI
IntHd(:) = IMDI
RelHd(:) = RMDI

! Read Namelist
Read( Unit=nft, Nml = headers )

! Write out namelist for diagnostic - Just during development use
! a low output level
If (PrintStatus >= PrStatus_Oper) Then
  If (mype == 0) Then
    Write (Unit=6, Nml=headers)
  End If
End If

Return
End Subroutine rcf_readnl_headers
End Module rcf_readnl_headers_mod
#endif
