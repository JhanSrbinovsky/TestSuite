#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Read the NLSTCATM namelist

MODULE Rcf_readnl_nlstcatm_Mod

!  Subroutine Rcf_Readnl_Nlstcatm - Read the NLSTCATM namelist
!
! Description:
!   Read the NLSTCATM namelist of Atmosphere model control variables.
!
! Method:
!   Reads variables into the Rcf_CntlAtm_Mod module.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   28/11/01   Support tropopause-based ozone.  D.Tan
!   5.5   21/01/03   L_USE_TPPS_OZONE moved to NLSTCATM. D.Robinson
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CONTAINS

Subroutine rcf_readnl_nlstcatm (nft)

Use Rcf_CntlAtm_Mod

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Oper

Use Rcf_Parvars_Mod, Only : &
    mype

Implicit None

integer nft

H_Sect = ' '
!
! Set defaults for tropopause-based ozone
I_Tpps_Ozone_Opts =       0 ! No options for tropopause-based ozone
! Set defaults for calculating ozone tracer and including ozone tracer
! in radiation calculation
L_USE_CARIOLLE=.FALSE.
L_USE_OZONEINRAD=.FALSE.

READ(Unit = nft, Nml = nlstcatm)
If (PrintStatus >= PrStatus_Oper .AND. mype == 0) Then
  write(Unit = 6, Nml = nlstcatm)
End If

Return
END SUBROUTINE  Rcf_readnl_nlstcatm
END MODULE Rcf_readnl_nlstcatm_Mod
#endif
