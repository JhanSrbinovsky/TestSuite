
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  Read the UANCNUM namelist

Module Rcf_readnl_Uancnum_Mod

!  Subroutine Rcf_Readnl_uancnum - reads the uancnum namelist
!
! Description:
!    Reads the uancnum namelist for user ancilmaster processing.
!
! Method:
!    Reads the variables for the Ancil_Mod module.
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

Subroutine Rcf_readnl_Uancnum( nft )

Use Ancil_mod, Only :    &
    n_uancil,            &
    uancfils,            &
    uancnum

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Oper

Use Rcf_Parvars_Mod, Only : &
    mype

Implicit None

Integer nft               ! Unit number
Integer i                 ! Looper

! Initialise and read UANCNUM namelist : number of user
!   ancilmaster files and total no. of user ancil records

! Initialisation
N_UANCIL     = 0
Do I = 1,20
  UANCFILS(I)='        '
End Do

! Read namelist
Read( Unit = nft, Nml = UANCNUM )

If (PrintStatus >= PrStatus_Oper .AND. mype == 0) then
  Write ( unit = 6 , Nml = UANCNUM )
Endif

Return
End Subroutine Rcf_readnl_Uancnum
End Module Rcf_readnl_Uancnum_Mod
