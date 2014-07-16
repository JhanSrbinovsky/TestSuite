#if defined(CONTROL) || defined(FLDOP)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Convert multiple of specified time period to no. of timesteps
! Function Interface:

      INTEGER FUNCTION TOTIMP(PERIOD,UNIT,MDL)
      IMPLICIT NONE

! Description: Get no. of timesteps from time period information taken
!  from STASH profiles to convert to no. of timesteps required to
!  generate STASH list contents for controlling diagnostic output times.
!
! Method: Simple conversion of specified time periods. Illegal
!  combinations are returned as -999 to be trapped by calling routine.
!
! Current code owner:  UM System Team
!
! History:
! Version   Date       Comment
! =======   ====       =======
!   3.5     Apr. 95    Original code.  S.J.Swarbrick
!   4.4     Oct. 97    Changed error return value from -1 to -999
!                          Shaun de Witt
!   4.5     18/08/98   Added DEF,FLDOP   (A Van der Wal)
!LL  5.1  22/02/00  Add PARPARM for TYPSIZE                 P.Burton
!   5.3   23/10/01  Trap negative period specified and minor tidy
!                   R Rawlins
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
! Global variables:
#include "csubmmax.h"
#include "cntlgen.h"


! Function arguments:
!   Scalar arguments with intent(in):
      INTEGER, INTENT(IN)::     PERIOD ! multiples of time period
      CHARACTER*2, INTENT(IN):: UNIT   ! descriptor for time period
      INTEGER,  INTENT(IN)::    MDL    ! internal model

! Local scalars:
      INTEGER, PARAMETER::      TOTIMP_ERROR=-999 ! error marker
      REAL FAC

!- End of Header ----------------------------------------------------

      IF( PERIOD  >=  0 ) THEN            ! check for valid multiples

        FAC= REAL(STEPS_PER_PERIODim(MDL))/REAL(SECS_PER_PERIODim(MDL))

        IF      (UNIT == 'T ') THEN       ! timesteps
        TOTIMP=PERIOD
        ELSE IF (UNIT == 'H ') THEN       ! hours
        TOTIMP=PERIOD*FAC*3600.0 +0.5
        ELSE IF (UNIT == 'DA') THEN       ! days
        TOTIMP=PERIOD*FAC*3600.0*24.0 +0.5
        ELSE IF (UNIT == 'DU') THEN       ! dump periods
        IF (DUMPFREQim(MDL) == 0) THEN
          WRITE(6,*)'TOTIMP:IRREGULAR DUMPS FOR DUMP FREQUENCY'
            TOTIMP= TOTIMP_ERROR
        ELSE
          TOTIMP=PERIOD*DUMPFREQim(MDL)
        END IF
        ELSE                              ! illegal unit
          TOTIMP= TOTIMP_ERROR
          WRITE(6,*)'TOTIMP: UNEXPECTED TIME UNIT=',UNIT
        END IF
! Note: the special case of period  ==  -1 (indefinite) is trapped by
! the calling routine, otherwise it would be necessary to include lines
! ELSE IF(PERIOD  ==  -1) THEN TOTIMP= -1 here.

      ELSE ! trap illegal negative periods

        WRITE(6,*) 'TOTIMP: ILLEGAL TIME PERIOD=',PERIOD
        TOTIMP= TOTIMP_ERROR

      END IF

      END FUNCTION TOTIMP

!- End of Function code ---------------------------------------------
#endif
