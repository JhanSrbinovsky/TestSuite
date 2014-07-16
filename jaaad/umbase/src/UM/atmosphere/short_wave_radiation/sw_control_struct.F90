#if defined(A70_1C) || defined(A70_1Z) \
 || defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   This module declares the controlling structures for SW radiation
!   and performs suitable default initialization.
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
MODULE sw_control_struct
!
!
  USE max_calls
  USE control_struc
!
!
  TYPE (control_option), Dimension(npd_swcall) :: sw_control
!
  INTEGER  n_swcall
!   The number of SW calls to the radiation
!
!
END MODULE sw_control_struct
!
!   Subroutine to set the default values of the control structure.
!
!END MODULE sw_control_struct
#endif
