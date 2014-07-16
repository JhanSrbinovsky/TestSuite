
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  Define headers namelist

Module rcf_headers_mod

! Description:
!    Define headers namelist
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   03/10/01   Addition of Implicit None. R.Sharp
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Implicit None

Integer, Save :: FixHd(256)
Integer, Save :: IntHd(100)
Real, Save    :: RelHd(100)

Namelist /Headers/ FixHd,IntHd,RelHd
End Module rcf_headers_mod
