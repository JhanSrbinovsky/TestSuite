
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Reads the USTSNUM namelist

Module Rcf_readnl_Ustsnum_Mod

!  Subroutine Rcf_Readnl_Ustsnum - reads the ustsnum namelist
!
! Description:
!   Reads information for user stashmaster processing.
!
! Method:
!   Read the namelist for variables in rcf_ppx_info_mod
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

Subroutine Rcf_readnl_Ustsnum( nft )

Use Rcf_Ppx_Info_mod, Only : &
    n_ustash,            &
    nrecs_ustash,        &
    ustsfils,            &
    ustsnum

Use Rcf_PrintStatus_Mod, ONly : &
    PrintStatus,            &
    PrStatus_Oper

Use Rcf_Parvars_Mod, Only : &
    mype

Implicit None

Integer nft               ! Unit number
Integer i                 ! Looper

! Open stash control file and read USTSNUM namelist: number of user
!   stash files and total no. of user stash records

! Initialisation
N_USTASH     = 0
NRECS_USTASH = 0
Do I = 1,20
  USTSFILS(I)='        '
End Do

! Read namelist
Read( Unit = nft, Nml = USTSNUM )
If (PrintStatus >= PrStatus_Oper .AND. mype == 0) Then
  Write( 6, Nml = USTSNUM )
End If

Return
End Subroutine Rcf_readnl_Ustsnum
End Module Rcf_readnl_Ustsnum_Mod
