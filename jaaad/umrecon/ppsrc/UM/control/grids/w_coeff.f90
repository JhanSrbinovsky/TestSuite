
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine W_COEFF------------------------------------------------
!LL
!LL  Purpose:
!LL          Calculates coefficients used to translate u and v compo-
!LL          nents of wind between equatorial (eq) latitude-longitude
!LL          grid and standard latitude-longitude grid (or visa versa).
!LL          Input latitudes and longitudes are in degrees.
!LL
!LL  Written by A. Dickinson
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL  5.1 17/04/00 Include c_pi for PI related variables. D. Robinson
!LL  5.3 24/09/01 Add reconfiguration to #defines. P.Selwood.
!LL  6.0 9/02/04  Add makebc to #defines. R. Sempers
!LL  6.3 18/12/06 Sign of COEFF2 reversed for phi_pole>90degs. See 
!LL               ticket #796. L. Jones
!LL
!LL  Documentation: The transformation formulae are described in
!LL                 unified model on-line documentation paper S1.
!LL
!LL  ------------------------------------------------------------------
!
!*L  Arguments:--------------------------------------------------------
      SUBROUTINE W_COEFF(                                               &
     & COEFF1,COEFF2,LAMBDA,LAMBDA_EQ,PHI_POLE,LAMBDA_POLE,POINTS)

      IMPLICIT NONE

      INTEGER                                                           &
     & POINTS            !IN  Number of points to be processed

      REAL                                                              &
     & COEFF1(POINTS)                                                   &
                         !OUT Coefficient of rotation no 1
     &,COEFF2(POINTS)                                                   &
                         !OUT Coefficient of rotation no 2
     &,LAMBDA(POINTS)                                                   &
                         !IN  Longitude
     &,LAMBDA_EQ(POINTS)                                                &
                         !IN  Longitude in equatorial lat-lon coords
     &,PHI_POLE                                                         &
                         !IN  Latitude of equatorial lat-lon pole
     &,LAMBDA_POLE       !IN  Longitude of equatorial lat-lon pole
! Workspace usage:-----------------------------------------------------
! None
!----------------------------------------------------------------------
! External subroutines called:-----------------------------------------
! None
!*---------------------------------------------------------------------
! Define local varables:-----------------------------------------------
      REAL A_LAMBDA,E_LAMBDA,SIN_E_LAMBDA,SIN_PHI_POLE                  &
     &    ,COS_PHI_POLE,C1,C2,LAMBDA_ZERO
      INTEGER I
!----------------------------------------------------------------------
! Constants from comdecks:---------------------------------------------
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
!----------------------------------------------------------------------

!L 1. Initialise local constants
!
! Longitude of zeroth meridian
      LAMBDA_ZERO=LAMBDA_POLE+180.

! Sine and cosine of latitude of eq pole

      SIN_PHI_POLE=SIN(PI_OVER_180*PHI_POLE)
      COS_PHI_POLE=COS(PI_OVER_180*PHI_POLE)

!L 2. Evaluate translation coefficients
!
      DO 200 I=1,POINTS

! Actual longitude converted to radians

      A_LAMBDA=PI_OVER_180*(LAMBDA(I)-LAMBDA_ZERO)

! Convert eq longitude to radians and take sine

      E_LAMBDA=LAMBDA_EQ(I)*PI_OVER_180
      SIN_E_LAMBDA=SIN(E_LAMBDA)

! Formulae used are from eqs (4.19) and (4.21)

      C1=SIN(A_LAMBDA)*SIN_E_LAMBDA*SIN_PHI_POLE                        &
     &           +COS(A_LAMBDA)*COS(E_LAMBDA)
      COEFF1(I)=C1
! avoid rounding error problems
      If (C1  >   0.9999) Then
        C1 = 1.0
        C2 = 0.0
      Else
        C2=SQRT(1.0-C1*C1)
      End IF

! Set the sign of C2. See ticket #796 for an explanation.
      COEFF2(I)=SIGN(C2,SIN_E_LAMBDA*COS_PHI_POLE)

200   CONTINUE

      RETURN
      END SUBROUTINE W_COEFF

