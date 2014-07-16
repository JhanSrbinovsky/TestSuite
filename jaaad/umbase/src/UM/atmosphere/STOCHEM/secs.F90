#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      FUNCTION SECS(DAY,IMONTH,IYEAR)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Returns number of seconds since
!-                        0000Z June 21st 0000, where IYEAR is the year,
!-                         no allowance for leap years.
!-   Inputs  : DAY,IMONTH,IYEAR
!-   Outputs :
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  3.4    31/08/94  Created. W.J. Collins
!  5.5    04/09/00  360/365 day automatic selection. C.E. Johnson
!  5.5    13/02/04  Uses daym_all to avoid summing daym repeatedly.
!                   K. Ketelsen
!  6.1    21/10/04  No change
!
!-
!VVV  V2.6  SECS 4/IX/99 360/365 selection.
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!----------------------------------------------------------------------
      REAL, INTENT(IN)    :: DAY
      INTEGER, INTENT(IN) :: IMONTH,IYEAR

      INTEGER :: I
      REAL :: SECS

      SECS=(REAL(daym_all)*(iyear-1)+                                   &
     &      1.0*SUM(DAYM((/(MOD(I-1,12)+1,I=6,IMONTH+12-1)/)))+         &
     &      DAY-21.0)*86400.0

      END FUNCTION SECS
#endif
