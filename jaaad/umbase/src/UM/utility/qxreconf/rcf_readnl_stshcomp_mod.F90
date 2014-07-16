#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Reads the STSHCOMP namelist

MODULE Rcf_readnl_stshcomp_Mod

!  Subroutine Rcf_Readnl_stshcomp - reads the STSHCOMP namelist.
!
! Description:
!   Read the STSHCOMP namelist for controlling section choices etc.
!
! Method:
!   Variables read into Rcf_Model_Mod module.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   28/11/01   Support tropopause-based ozone.  D.Tan
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CONTAINS

!+  Read in the STSHCOMP namelists.

Subroutine rcf_readnl_stshcomp (nft)

Use Rcf_Model_Mod, Only : &
    STSHCOMP,                                         &
    Indep_SR,     Atmos_SR,     Ocean_SR,             &
    Slab_SR,      Wave_SR,                            &
    BSPMSL,       CCEW,         IDO,                  &
    LOSSM,        MESO,         MLMO,                 &
    OAFLD,        OCARB,        OSFC,                 &
    SCAL,         UPD175,                             &
    ZONAVOZONE,   ZONAVTPPSOZONE

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Oper

Use Rcf_Parvars_Mod, Only : &
    mype

Implicit None

Integer nft

Indep_SR = ' '
Atmos_SR = ' '
Ocean_SR = ' '
Slab_SR  = ' '
Wave_SR  = ' '

BSPMSL = ' '
CCEW   = ' '
IDO    = ' '
LOSSM  = ' '
MESO   = ' '
MLMO   = ' '
OAFLD  = ' '
OCARB  = ' '
OSFC   = ' '
SCAL   = ' '
UPD175 = ' '

READ(Unit = nft, Nml = STSHCOMP)
If (PrintStatus >= PrStatus_Oper .AND. mype == 0) Then
  write(Unit = 6, Nml = STSHCOMP)
End If

! Set values for tropopause-based ozone.
! In future developments, obtain directly from namelist.
ZonAvTppsOzone=ZonAvOzone

Return

END SUBROUTINE  Rcf_readnl_stshcomp

END MODULE Rcf_readnl_stshcomp_Mod
#endif
