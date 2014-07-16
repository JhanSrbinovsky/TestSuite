#if defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Global data module for varying solar constant.

MODULE volcts_Mod

  IMPLICIT NONE
  SAVE

! Description:
!   Global data necessary for varying volcanic forcing
!
! Method:
!   reads data from text file for use in climate change runs.
!
! Current Code Owner: Chris Jones
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 6.6.2     05/06/09  Original code.  Chris Jones
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v7.4 programming standards.
!
! Declarations:
!
! Global variables (#include statements etc):

      REAL VOLCTS(4,12,1850:2300)

! SCVARY: Array of variable solar constant, from 1700 to 2300.
!
!- End of header

END MODULE volcts_Mod
#endif
