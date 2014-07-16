#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
      SUBROUTINE READCNTL(UNIT,ICODE,CMESSAGE)

      IMPLICIT NONE
!
! Description:
!  Reads overall, generic and model-specific control variables.
!
! Method:
!  Reads namelist containing all overall, generic and model-specific
!  control variables not in History file.  Control variables are then
!  available from COMMON block storage. Declarations, namelist and
!  common are all in included comdecks.
!
! Current Code Owner: R.T.H.Barnes
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
#include "ccontrol.h"
#include "mppconf.h"

! Subroutine arguments
!   Scalar arguments with intent(in):
!   Array  arguments with intent(in):
!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
!   Scalar arguments with intent(out):
      INTEGER UNIT  ! namelist input unit no.
      INTEGER ICODE ! return code
      CHARACTER*80 CMESSAGE ! error message
!   Array  arguments with intent(out):

! Local parameters:

! Local scalars:
      CHARACTER*80 FILENAME ! physical filename

! Local dynamic arrays:

! Function & Subroutine calls:
      External GET_FILE

!- End of header

!  Find filename and open for formatted i/o

      CALL GET_FILE(UNIT,FILENAME,80,ICODE)
      OPEN(UNIT,FILE=FILENAME,FORM='FORMATTED',IOSTAT=ICODE)

!  Read overall control data into COMMON/CNTLCALL/
      PP_LEN2_LOOK(:) = 0
      READ(UNIT,NLSTCALL)

!  Read control information for WRITD1 diagnostics
!  but only from unit 5, not for continuation runs from unit 9

      IF (UNIT  ==  5) THEN
! DEPENDS ON: readwritd
      CALL READWRITD(UNIT,ICODE,CMESSAGE)
      END IF

!  Read generic control data into COMMON/CNTLCGEN/

      READ(UNIT,NLSTCGEN)

! Read MPP configuration data into COMMON /COM_MPPCONF/
      GCOM_COLL_LIMIT = 1     ! Default value
      READ (UNIT,NLST_MPP)
!  Read atmospheric control data into COMMON/CNTLCATM/

#if defined(ATMOS)
!     Iteration Count output, L_MR_PHYSICS2 and STPH_RP defaults
      l_icount = .false.
      L_MR_PHYSICS2 = .FALSE.
      L_RPSEED_READ  = .FALSE.
      L_RPSEED_WRITE = .TRUE.
!     Height threshold default for pmsl geostrophic wind calculation
      npmsl_height = 500.0
! Natural forcing
      L_SCVARY=.FALSE.
      L_VOLCTS=.FALSE.
      READ(UNIT,NLSTCATM)
#endif

!  Read ocean control data into COMMON/CNTLCOCN/

#if defined(OCEAN)
      READ(UNIT,NLSTCOCN)
#endif

!  Read slab control data into COMMON/CNTLCSLB/

#if defined(SLAB)
      READ(UNIT,NLSTCSLB)
      READ(UNIT,NLSTCICES)
#endif

!  Read wave control data into COMMON/CNTLCWAV/

      REWIND(UNIT)

      RETURN
      END SUBROUTINE READCNTL

#endif
