
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      FUNCTION COAG_QI(A, B, RMEDAIT, RMEDACC, SAITSQ, SACCSQ)
!
!---------------------------------------------------------------------
! Purpose: To calculate coagulation coefficients for sulphate aerosol.
!          Called by SULPHR
!
! Current owners of code:                         D Roberts, M Woodage
!
! History:
! Version    Date     Comment
! -------    ----     -------
!  5.4    05/09/02    New deck                                M Woodage
!
!  6.2    20/10/05    DEF for a17_2b added.                   A Jones
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!
! System components covered:
!
! System task:
!
! Documentation: Not yet available
!---------------------------------------------------------------------
!
      IMPLICIT NONE
!
! Comdecks:
!*L------------------COMDECK C_PI---------------------------------------
!LL
!LL 4.0 19/09/95  New value for PI. Old value incorrect
!LL               from 12th decimal place. D. Robinson
!LL 5.1 7/03/00   Fixed/Free format P.Selwood
!LL

      ! Pi
      Real, Parameter :: Pi                 = 3.14159265358979323846

      ! Conversion factor degrees to radians
      Real, Parameter :: Pi_Over_180        = Pi/180.0

      ! Conversion factor radians to degrees
      Real, Parameter :: Recip_Pi_Over_180  = 180.0/Pi

!*----------------------------------------------------------------------
!
! Arguments with intent IN:
      REAL                                                              &
     &    A                                                             &
     &,   B                                                             &
     &,   RMEDAIT                                                       &
     &,   RMEDACC                                                       &
     &,   SAITSQ                                                        &
     &,   SACCSQ
!
! Arguments with intent OUT:
!
      REAL COAG_QI      !Returned function required for coag caln
!
! Local variables:
!
      REAL                                                              &
     &  POWER1                                                          &
     &, POWER3                                                          &
     &, POWER4                                                          &
     &, Q1                                                              &
     &, Q3                                                              &
     &, Q4                                                              &
     &, Q5                                                              &
     &, SPREAD
!
      POWER1 = A + B - 2.0
      POWER3 = 3.0*(A - 1.0)
      POWER4 = 3.0*(B - 1.0)
!
      Q1 = (1.333333*PI)**POWER1
      Q3 = RMEDAIT**POWER3
      Q4 = RMEDACC**POWER4
!
      SPREAD = SAITSQ*(A*A - 1.0) + SACCSQ*(B*B - 1.0)
      SPREAD = SPREAD*4.5
      Q5 = EXP(SPREAD)
!
      COAG_QI = Q1*Q3*Q4*Q5
!
      RETURN
      END FUNCTION COAG_QI
