#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Magic numbers  defining decompositions

Module Rcf_DecompTP_Mod

! Description:
!    This data module contains magic numbers defining decompositions
!    for the MPP reconfiguration
!
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   03/10/01   Addition of Implicit None and removal of unused
!                    header template.                            R.Sharp
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Implicit None

Integer, Parameter :: max_decomps  = 2
Integer, Parameter :: decomp_unset = -1
Integer, Parameter :: decomp_rcf_input  = 1
Integer, Parameter :: decomp_rcf_output = 2

End Module Rcf_DecompTP_Mod
#endif
