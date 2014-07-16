
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      SUBROUTINE HYGRO_FACT(HUMIDITY, ALPHA)
!
!---------------------------------------------------------------------
! Purpose: To calculate hygroscopic growth factor for sulphate aerosol.
!          Called by SULPHR
!
! Current owners of code:                         D Roberts, M Woodage
!
! History:
! Version    Date     Comment
! -------    ----     -------
!  5.4    05/09/02    New deck                                M Woodage
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
! Arguments with intent IN:
      REAL                                                              &
     &  HUMIDITY              !RELATIVE HUMIDITY (0-1)
!
! Arguments with intent OUT:
      REAL                                                              &
     &  ALPHA                 !FITZGERALD'S ALPHA PARAMETER
!
! Local variables:
      REAL                                                              &
     & ALPHA81                !VALUE OF ALPHA AT HUMIDITY 0.81
!
!     STATEMENT FUNCTION SPECIFICATIONS.
!
      REAL                                                              &
     & H                                                                &
                              !DUMMY VARIABLE FOR HUMIDITY
     &,AFUNC1
!
      AFUNC1(H) = 1.2E+00*EXP( (0.066E+00*H)/(1.058E+00 - H) )
!
! This routine calculates growth factors for AMMONIUM SULPHATE
! (based on routine GROW_PARTICLES by D. L. Roberts (1995)).
! Note no growth takes place below a humidity of 0.3.
! Above the deliquescence point, taken as 0.81, the scheme is the
! one due to Fitzgerald (1975).
! This is a simplified version for use with RH < 0.97, assuming that
! BETA=1, and only ALPHA needs to be calculated.
!
      ALPHA = 1.0E+00
!
      IF ( HUMIDITY  >=  0.3E+00 ) THEN
!
        IF ( HUMIDITY  <   0.81E+00 ) THEN
!
! Calculate ALPHA
! We have to be careful here to use the actual value of alpha
! at 0.81 otherwise the function would become discontinuous at 0.81.
!
          ALPHA81 = AFUNC1(0.81E+00)
          ALPHA = 1.0E+00 + (HUMIDITY-0.3E+00)                          &
     &                       *(ALPHA81-1.0E+00)/0.51E+00
!
        ELSE IF ( HUMIDITY  <=  0.97E+00 ) THEN
!
          ALPHA = AFUNC1(HUMIDITY)
!
        END IF
!
      END IF           ! End RH > 0.3 condition
!
      RETURN
      END SUBROUTINE HYGRO_FACT
