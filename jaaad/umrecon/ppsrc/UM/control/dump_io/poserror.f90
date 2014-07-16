

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE POSERROR---------------------------------------
!LL
!LL  Purpose:
!LL           Prints out a message when position of a data block as
!LL           pointed to by fixed length header differs from actual
!LL           position in model dump.
!LL
!LL  Written by A. Dickinson 29/12/89
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL   4.4    15/10/97 Added code to print the error message to
!LL                   stderr, and call abort in case all the
!LL                   PE's have not detected the error condition.
!LL                     Author: Bob Carruthers, Cray Research
!LL
!LL  Programming standard:
!LL           Unified Model Documentation Paper No 3
!LL           Version No 1 15/1/90
!LL
!LL  System component: E4
!LL
!LL  System task: F3
!LL
!LL  Documentation:
!LL           None
!LL------------------------------------------------------------
!*L Arguments:-------------------------------------------------
      SUBROUTINE POSERROR(STRING,START_BLOCK,HEAD_POS,HEAD_ADDRESS)

      IMPLICIT NONE

      INTEGER                                                           &
     & START_BLOCK                                                      &
                    !IN Actual position of data block
     &,HEAD_POS                                                         &
                    !IN Position in FIXHD of pointer
     &,HEAD_ADDRESS !IN Position in file pointed to by FIXHD(HEAD_POS)




      CHARACTER*(80) STRING  !IN Description of block


! -------------------------------------------------------------
! Workspace usage:---------------------------------------------
! None
! -------------------------------------------------------------
!*L External subroutines called:-------------------------------
! None
!*-------------------------------------------------------------

!L Internal structure: none

      WRITE(6,'('' ******FATAL ERROR WHEN READING MODEL DUMP******'')')
      WRITE(6,'('' Conflict between start position of '',A)')STRING
      WRITE(6,'('' block and pointer in fixed length header: FIXHD('',  &
     &I3,'') ='',I9)')HEAD_POS,HEAD_ADDRESS
      WRITE(6,'('' Current position in file ='',I9,'' words in'')')     &
     &START_BLOCK
      WRITE(6,'('' ***********************************************'')')

      RETURN
      END SUBROUTINE POSERROR

