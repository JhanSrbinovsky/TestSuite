
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Wrapper for architecture specific initialisation

Module Rcf_Arch_Initialise_Mod

!  Subroutine Rcf_Arch_Initialise - architecture specifics
!
! Description:
!   A wrapper for architecture specpific initialisation routines
!
! Method:
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

Subroutine Rcf_Arch_Initialise()

Use Rcf_T3E_Init_Mod, Only : &
    Rcf_T3E_Init

Implicit None






Return
End Subroutine Rcf_Arch_Initialise

End Module Rcf_Arch_Initialise_Mod
