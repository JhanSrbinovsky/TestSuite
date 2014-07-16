
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Module containing the O3INTP deck

Module Rcf_o3intp_Mod

! Description:
!   Module to contain the O3INTP deck
!
! Current Code Owner: J. M. Edwards
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.5   13/03/03   Original code.  Jean-Claude Thelen

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Implicit None

!+ ---------------------------------------------------------------------
!  Module to specify allowed methods of interpolating from the
!  ancillary file.
!
!  Current Code Owner: J. M. Edwards
!
!  History:
!
!  Version  Date      Comment.
!  5.2      14/11/00  Original code.
!                     (J. M. Edwards)
!
!- ---------------------------------------------------------------------
!
      Integer, parameter :: IO3_3DSPEC = 1
!       Ozone is provided as a full 3D field.
      Integer, parameter :: IO3_2DSPEC = 2
!       Ozone is expanded from a 2D field by direct copying.
      Integer, parameter :: IO3_2DMASSCON = 3
!       Ozone is expanded from a 2D field with conservation of mass.
      Integer, parameter :: IO3_TROP_MAP = 4
!       Ozone mixing ratios are set by mapping each height at each
!       grid-point in the real profile to a height in the ancillary
!       profile and using the mixing ratio there. The mapping is
!       set using the height of the tropopause
      Integer, parameter :: IO3_TROP_MAP_MASSCON = 5
!       Ozone mixing ratios are set as above, but scaled so as to
!       preserve the vertically integrated column ozone
!
! ----------------------------------------------------------------------

End Module Rcf_o3intp_Mod
