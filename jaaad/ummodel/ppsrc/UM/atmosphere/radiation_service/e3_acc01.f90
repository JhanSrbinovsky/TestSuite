


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to calculate third expoential integral.
!
! Method:
!       For small arguments a power series is used. For larger
!       arguments a Pade approximant derived from the asymptotic
!       expansion is used.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77 with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      FUNCTION E3_ACC01(X)
!
!
      IMPLICIT NONE
!
!
!     DUMMY ARGUMENTS
      REAL                                                              &
     &     X                                                            &
!             POINT OF EVALUATION
     &   , E3_ACC01
!             NAME OF FUNCTION
!
!     LOCAL VARIABLES
      REAL                                                              &
     &     EULER
!             EULER'S CONSTANT
!
      PARAMETER(EULER=0.5772156E+00)
!
!
      IF (X <  1.0E-06) THEN
         E3_ACC01=0.5E+00
      ELSE IF (X <  2.0E+00) THEN
         E3_ACC01=-0.5E+00*X*X*LOG(X)+0.5E+00                           &
     &      +X*(-1.0E+00+X*(0.75E+00-0.5E+00*EULER                      &
     &      +X*(1.0E+00/6.0E+00-X*(1.0E+00/48.0E+00                     &
     &      -X*(1.0E+00/3.60E+02-X*(1.0E+00/2.880E+03                   &
     &      -X*(1.0E+00/2.5200E+04)))))))
      ELSE
!        WE USE A DOUBLY CUBIC PADE APPROXIMANT DERIVED FROM THE
!        ASYMPTOTIC EXPRESSION.
         E3_ACC01=(EXP(-X)/X)                                           &
     &      *(6.0E+00+X*(48.0E+00+X*(15.0E+00+X)))                      &
     &      /(120.0E+00+X*(90.0E+00+X*(18.0E+00+X)))
      ENDIF
!
!
      RETURN
      END FUNCTION E3_ACC01
