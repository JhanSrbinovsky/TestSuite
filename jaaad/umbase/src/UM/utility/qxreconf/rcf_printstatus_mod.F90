#if defined(RECON) || defined(VAROPSVER)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Define variables for level of output

Module Rcf_PrintStatus_Mod

! Description:
!   VAriables to control timer and level of standard output
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   03/10/01   Addition of Implicit None. R.Sharp
!   5.3   25/09/01   Added L_IO_Timer. P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


! PrintStatus_Mod
! This module holds some parameters for the new Standard Output
! level information

Implicit None

Integer, Parameter   :: PrStatus_Min    = 1  ! Minimum output
Integer, Parameter   :: PrStatus_Normal = 2  ! Short informative output
Integer, Parameter   :: PrStatus_Oper   = 3  ! Full informative output
Integer, Parameter   :: PrStatus_Diag   = 4  ! Extra Diagnostic output

! The variable that controls message output - default setting
! is to 2 - Full informative output.

Integer, Save        :: PrintStatus     = PrStatus_Normal


! The variables that controls timing (and thus output levels to
! of timing information) are logicals

Logical, Save        :: LTimer       ! Subroutine timings
Logical, Save        :: L_IO_Timer   ! IO timings

End Module Rcf_PrintStatus_Mod
#endif
