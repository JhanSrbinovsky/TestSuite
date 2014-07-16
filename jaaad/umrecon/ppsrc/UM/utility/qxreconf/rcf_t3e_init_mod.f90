
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ T3E specific initialisation

Module Rcf_T3E_Init_Mod

!  Subroutine Rcf_T3E_Init - T3E initialisation tasks
!
! Description:
! This module does T3E specific initialisation - specifically
! specifying resolution for the timer and splitting the output
! streams - 1 per pe
!
! Method:
!   Similar to code in UM_SHELL from the UM.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   6.0   07/01/04   (re)Move section to send unit 6 output to a unique
!                    file for each PE. R.Sharp
!   6.2   11/08/05   Replace quotes in includes with angle brackets.
!                    P.Selwood.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_T3E_Init()

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Parvars_Mod, Only : &
    nproc,              &
    mype

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Normal

Implicit None

Return
End Subroutine Rcf_T3E_Init
End Module Rcf_T3E_Init_Mod
