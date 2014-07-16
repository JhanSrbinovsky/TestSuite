#if defined(C70_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine: DEL_HIST -------------------------------------------
!LL
!LL  Purpose: delete a history file -called if problems writing
!LL           out partial sums
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 6.1.5A
!LL
!LL  Author:   R A Stratton
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL   3.1  05/02/93    Portable Fortran unit no assigns
!LL                    Author: A. Dickinson    Reviewer: R. Stratton
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: H
!LL
!LL  Project task: H
!LL
!LL  External documentation:
!LL
!LLEND------------------------------------------------------------------
!
!*L  Interface and arguments: ------------------------------------------
      SUBROUTINE DEL_HIST(HUNIT)
!
      IMPLICIT NONE
!
      INTEGER                                                           &
     &       HUNIT    ! IN  - history file unit number
!
! NOTE no error code and error message are passed back from this routine
! as this routine is only called in the event of a non-zero error
! code in U_MODEL.
!*----------------------------------------------------------------------
!  Common blocks
!
!  Subroutines called
!
      EXTERNAL GET_FILE
!
!  Local variables
!
      INTEGER   ICODE        ! error from open and close
      CHARACTER*80 CMESSAGE  ! error message this routine only
      CHARACTER*80 FILENAME
!
!L----------------------------------------------------------------------
!L 1. Open  history file
!L
      CALL GET_FILE(HUNIT,FILENAME,80,ICODE)
        OPEN(HUNIT,FILE=FILENAME,FORM='UNFORMATTED',IOSTAT=ICODE)
      IF (ICODE /= 0) THEN
        CMESSAGE='DELHIST: failure to open history file prior to its del&
     &etion'
        WRITE(6,*)CMESSAGE
      ENDIF
!L
!L 2. Close history file deleting it in the process.
!L
      CLOSE(HUNIT,IOSTAT=ICODE,STATUS='DELETE')
      IF (ICODE /= 0) THEN
        CMESSAGE='DELHIST: failure to delete history file'
        WRITE(6,*)CMESSAGE
      ENDIF

      RETURN
!L----------------------------------------------------------------------
      END SUBROUTINE DEL_HIST
!
#endif
