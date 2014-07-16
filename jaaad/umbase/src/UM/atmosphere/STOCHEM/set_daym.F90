#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE SET_DAYM(DAYM,YEAR)
! ----------------------------------------------------------------------
!- Purpose:
!-   To set DAYM array (days in each month).
!-
!
! Current Owner of Code: C.E. Johnson
!
! History:
! Version   Date                    Comment
!  5.0    28/10/00  Created.  C.E. Johnson
!  5.2    25/07/01  Changed for new UM version. C.E. Johnson
!  5.5    13/02/04  Now calls SET_DAYM_ALL to calculate number of
!                   days in a year (360 or 365). K. Ketelsen
!  6.1    20/10/04  No change.
!
! ----------------------------------------------------------------------
      USE      IN_STOCHEM_GRD, only        :  set_daym_all
      INTEGER, DIMENSION(12), INTENT(OUT)  :: DAYM
      INTEGER,                INTENT(IN)   :: YEAR

      INTEGER, DIMENSION(12), PARAMETER ::                              &
     &  DAYM360=(/30,30,30,30,30,30,30,30,30,30,30,30/),                &
                                                          ! 360 days
     &  DAYM365=(/31,28,31,30,31,30,31,31,30,31,30,31/)   ! 365 days

#include "chsunits.h"
#include "cntlall.h"

! Set days of month array
      IF (LCAL360) THEN
        DAYM = DAYM360
      ELSE
        DAYM = DAYM365
      END IF
      CALL SET_DAYM_ALL              !kk

! Check for leap year and reset February length (365 day years)
      IF (.NOT. LCAL360) THEN
        IF (MOD(YEAR,4) == 0 .AND. MOD(YEAR,100) /= 0 .OR.              &
     &      MOD(YEAR,1000) == 0) THEN    ! leap year
          DAYM(2)=29
        ELSE
          DAYM(2)=28                     ! other years
        END IF
      END IF

      RETURN

      END SUBROUTINE SET_DAYM
#endif
