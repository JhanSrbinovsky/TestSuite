#if defined(A70_1C) || defined(A70_1Z) \
 || defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to set the maximum number of calls to the radiation.
!
! Description:
!   This module declares the maximum number of calls to the
!   radiation code permitted on a single timestep.
!
! Current Code Owner: J.-C. Thelen
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   6.2   20/08/04   Original code.
!                        J.-C. Thelen
!
! Code Description:
!   Language: FORTRAN 90
!
!- End of header
!
!
MODULE max_calls
!
!
  INTEGER, Parameter :: npd_swcall=2
!   Size allocated for arrays concerned with the SW call
  INTEGER, Parameter :: npd_lwcall=2
!   Size allocated for arrays concerned with the LW call
!
END MODULE max_calls
#endif
