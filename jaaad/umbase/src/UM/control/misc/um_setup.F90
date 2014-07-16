#if defined(C70_1A) && !defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Reads History & control files; also interim for CRUN and HK for Op.
!
! Subroutine Interface:
      SUBROUTINE UM_SETUP(ICODE,CMESSAGE)

      IMPLICIT NONE
!
! Description:
!   Reads History and control namelist files for all runs.
!   For CRUNs also reads interim control file of values changed by user.
!   For operational runs reads housekeeping file.
!
! Method: as above
!
! Current Code Owner: RTHBarnes.
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered: C0
! System Task:              C0
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):
#include "csubmodl.h"
#include "chsunits.h"
#include "chistory.h"
#include "ccontrol.h"

! Subroutine arguments
!   Scalar arguments with intent(in):
!   Array  arguments with intent(in):
!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
!   Scalar arguments with intent(out):
      Integer ICODE ! Error code
      Character*80 CMESSAGE ! Error message
!   Array  arguments with intent(out):

! Local parameters:

! Local scalars:
      INTEGER I              ! Loop variable
      INTEGER UNITICTL       ! Unit no. for Interim Control file
      LOGICAL CRUN           ! T if continuation run
      INTEGER J              ! dummy for use with FORT_GET_ENV
      CHARACTER*4  TYPE      ! Type of run - NRUN or CRUN
      CHARACTER*(*) RoutineName
      PARAMETER (RoutineName='UM_SETUP')

! Local dynamic arrays:

!- End of header
!
! ----------------------------------------------------------------------
!  0. Initialise ICODE
!
      ICODE=0
! ----------------------------------------------------------------------
!  1. Read History file.
!
! DEPENDS ON: readhist
      CALL READHIST ( IHIST_UNIT,ICODE,CMESSAGE )
      IF (ICODE  >   0) GO TO 9999

! ----------------------------------------------------------------------
!  2. Read Control file on standard input.
!
! DEPENDS ON: readcntl
      CALL READCNTL ( 5,ICODE,CMESSAGE )
      IF (ICODE  >   0) GO TO 9999

! ----------------------------------------------------------------------
!  2.5 Call TIMER now that LTIMER has been read
!
      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('UM_SETUP',3)
      END IF

! ----------------------------------------------------------------------
!  3. Read Interim Control file for CRUNs only.
!
      CALL FORT_GET_ENV('TYPE',4,TYPE,4,J)
      CRUN = .false.
      DO  I = 1,N_INTERNAL_MODEL_MAX
        IF (H_STEPim(I)  >  0) THEN
          CRUN = .true.
          GO TO 300
        END IF
      END DO
  300 CONTINUE
      IF (CRUN .or. TYPE == 'CRUN') THEN
        WRITE(6,*)' UMSETUP; CRUN, read CONTCNTL from unit 9'
!  Open and read Interim Control file
        UNITICTL = 9
! DEPENDS ON: readcntl
        CALL READCNTL ( UNITICTL,ICODE,CMESSAGE )
        IF (ICODE  >   0) GO TO 9999
        CLOSE (UNITICTL)
      END IF

! ----------------------------------------------------------------------
!  4. Read Housekeeping file for operational runs only.
!
      IF (MODEL_STATUS  ==  'Operational') THEN
! DEPENDS ON: readhk
        CALL READHK(HKFILE_UNIT,ICODE,CMESSAGE)
        IF (ICODE  >   0) GO TO 9999
      END IF

! ----------------------------------------------------------------------
!  5. Output any error/warning message and stop if error.
!
 9999 CONTINUE
! DEPENDS ON: ereport
      IF (ICODE  /=  0) CALL EREPORT(RoutineName,ICODE,CMESSAGE)
      IF (ICODE  >   0) CALL ABORT

! ----------------------------------------------------------------------
!  6. Transfer end of TIMER section to UM_INDEX to get full SETUP time
!

      RETURN
      END SUBROUTINE UM_SETUP
#endif
