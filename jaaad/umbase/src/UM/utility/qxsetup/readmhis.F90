#if defined(SETUP)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Read master history namelist input file.
!
! Subroutine Interface:
      SUBROUTINE READMHIS( UNITHIST,ICODE,CMESSAGE )

      IMPLICIT NONE
!
! Description:
!  Read master history namelist input file to set up history
!   common block variables
!
! Method:
!
! Current Code Owner: R.T.H.Barnes.
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  3.5  23/05/95  Original code for submodels stage 1. Based on old
!                  routine READMCTL.  RTHBarnes.
!  4.0  18/10/95  Add ICODE error return to GET_FILE call. RTHBarnes
!  5.3  19/06/01  Added Filename to error message. Dave Tan
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered: <appropriate code>
! System Task:              <appropriate code>
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):
#include "csubmodl.h"
#include "chsunits.h"
#include "chistory.h"

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER      UNITHIST  ! IN  - Master history file unit no.
!   Array  arguments with intent(in):
!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
!   Scalar arguments with intent(out):
      INTEGER       ICODE    ! Out - Return code from routine
      CHARACTER*80 CMESSAGE ! Out - Return message if failure occured
!   Array  arguments with intent(out):

!   ErrorStatus <Delete if ErrorStatus not used>
      INTEGER      ErrorStatus          ! Error flag (0 = OK)

! Local parameters:

! Local scalars:
      CHARACTER *80 FILENAME

! Local dynamic arrays:

! Function & Subroutine calls:
      External GET_FILE

!- End of header
!L
!L 1. Read Master Control file namelist information
!L
!L  NLIHISTO : Integer overall model history data
!L  NLCHISTO : Character overall model history data
!L  NLIHISTG : Integer generic model history data
!L  NLCHISTG : Character generic model history data
!L  NLCFILES : Character variables used for logical filenames
!L
      CALL GET_FILE(UNITHIST,FILENAME,80,ICODE)
      OPEN(UNITHIST,FILE=FILENAME,IOSTAT=ICODE)
!
! Check for error
!
      IF(ICODE  >  0)THEN
        CMESSAGE='READMHIS  : Failed in OPEN of Master History File'
        GOTO 999
      ELSEIF(ICODE  <   0)THEN
        WRITE(6,*)'READMHIS  :                                          &
     &           Warning message on OPEN of Master History File'
        WRITE(6,*)'IOSTAT= ',ICODE
      ENDIF
!
      CMESSAGE='READMHIS: Problem reading namelist NLIHISTO '//FILENAME
      READ(UNITHIST,NLIHISTO,ERR=999)
!
      CMESSAGE='READMHIS: Problem reading namelist NLCHISTO '//FILENAME
      READ(UNITHIST,NLCHISTO,ERR=999)
!
      CMESSAGE='READMHIS: Problem reading namelist NLIHISTG '//FILENAME
      READ(UNITHIST,NLIHISTG,ERR=999)
!
      CMESSAGE='READMHIS: Problem reading namelist NLCHISTG '//FILENAME
      READ(UNITHIST,NLCHISTG,ERR=999)
!
      CMESSAGE='READMHIS: Problem reading namelist NLCFILES '//FILENAME
      READ(UNITHIST,NLCFILES,ERR=999)
!
!     Normal return
!
      ICODE=0
      CMESSAGE='READMHIS: Normal return'
      RETURN
!
!     Error return
!
 999  ICODE=1
      RETURN
      END SUBROUTINE READMHIS
#endif
