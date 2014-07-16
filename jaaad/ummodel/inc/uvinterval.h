#if defined(A01_3C) || defined(A01_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!    Define the upper and lower boundaries of the interval in which the
!    fluxes should be calculated.
!
! Method:
!
! Current Code Owner: Jean-Claude Thelen
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 6.2       13/02/06  Original code.  Jean-Claude Thelen
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v7.4 programming standards.
!
! Declarations:
!
! Global variables (#include statements etc):
!

!
!    Define the upper and lower boundaries of the interval in which the
!    fluxes should be calculated.
!
        Real, Parameter :: UV_INTERVAL_SHORT = 2.0E-07
        Real, Parameter :: UV_INTERVAL_LONG  = 3.2E-07
#endif
