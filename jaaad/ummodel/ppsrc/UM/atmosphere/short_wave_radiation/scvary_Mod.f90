
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Global data module for varying solar constant.

MODULE scvary_Mod

  IMPLICIT NONE
  SAVE

! Description:
!   Global data necessary for varying solar constant
!
! Method:
!   reads data from text file for use in climate change runs.
!
! Current Code Owner: Chris Jones
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 6.6.2     04/06/09  Original code.  Chris Jones
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v7.4 programming standards.
!
! Declarations:
!
! Global variables (#include statements etc):

      REAL SCVARY (601)

! SCVARY: Array of variable solar constant, from 1700 to 2300.
!
!- End of header

END MODULE scvary_Mod
