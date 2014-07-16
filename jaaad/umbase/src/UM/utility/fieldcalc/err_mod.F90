#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module containing error codes for Fieldcalc

MODULE Err_Mod

! Description:
!
! Method:
!
! Owner: Dave Robinson
!
! History:
! Version Date     Comment
! ------- ----     -------
! 1.0     02/05/03 Original Code.  Sara James
! 6.0     12/09/03 Code implemented into UM. Dave Robinson
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

IMPLICIT None

! Error codes
INTEGER, PARAMETER :: StatusOK      =  0
INTEGER, PARAMETER :: StatusWarning = -9
INTEGER, PARAMETER :: StatusFatal   =  9
INTEGER, PARAMETER :: EndofFile     = -1

END MODULE Err_Mod
#endif
