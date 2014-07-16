#if defined(A70_1C) || defined(A70_1Z) \
 || defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Description:
!   This module defines the LW spectrum for each call to radiation.
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
MODULE spec_sw_lw
!
!
! Define the type of the spectrum

  USE max_calls
  USE dec_spec
!
!
  TYPE (spectrum), Dimension(npd_swcall) :: sw_spectrum
  TYPE (spectrum), Dimension(npd_lwcall) :: lw_spectrum
!
!
!
END MODULE spec_sw_lw
#endif
