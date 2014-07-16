#if defined(C92_2A) || defined(RECON) || defined(VAROPSVER) \
 || defined(MAKEBC) || defined(FRAMES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine EQTOLL-------------------------------------------------
!LL
!LL  Purpose:  Calculates latitude and longitude on standard grid
!LL            from input arrays of latitude and longitude on
!LL            equatorial latitude-longitude (eq) grid used
!LL            in regional models. Both input and output latitudes
!LL            and longitudes are in degrees.
!LL
!LL  Written by A. Dickinson
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL    5.4  11/09/02  Cater for negative pole latitude. D.Robinson
!LL    6.0  05/09/03  Added new def for use with makebc. R Sempers
!LL
!LL  Documentation: The transformation formulae are described in
!LL                 unified model on-line documentation paper S1.
!LL
!LL Logical components covered : S131
!LL
!LL Project task :
!LL
!LL External documentation:
!LL
!LLEND-----------------------------------------------------------------
!
!*L  Arguments:--------------------------------------------------------
      SUBROUTINE EQTOLL                                                 &
     &(PHI_EQ,LAMBDA_EQ,PHI,LAMBDA,PHI_POLE,LAMBDA_POLE,POINTS)

      IMPLICIT NONE

      INTEGER                                                           &
     & POINTS            !IN  Number of points to be processed

      REAL                                                              &
     & PHI(POINTS)                                                      &
                         !OUT Latitude
     &,LAMBDA(POINTS)                                                   &
                         !OUT Longitude (0 =< LON < 360)
     &,LAMBDA_EQ(POINTS)                                                &
                         !IN  Longitude in equatorial lat-lon coords
     &,PHI_EQ(POINTS)                                                   &
                         !IN  Latitude in equatorial lat-lon coords
     &,PHI_POLE                                                         &
                         !IN  Latitude of equatorial lat-lon pole
     &,LAMBDA_POLE       !IN  Longitude of equatorial lat-lon pole

! Workspace usage:-----------------------------------------------------
! None
!----------------------------------------------------------------------
! External subroutines called:-----------------------------------------
! None
!*---------------------------------------------------------------------
! Local varables:------------------------------------------------------
      REAL E_LAMBDA,E_PHI,A_LAMBDA,ARG,A_PHI,SIN_PHI_POLE,COS_PHI_POLE
      REAL TERM1,TERM2,SMALL,LAMBDA_ZERO
      INTEGER I
      PARAMETER(SMALL=1.0E-6)
!----------------------------------------------------------------------
! Constants from comdecks:---------------------------------------------
#include "c_pi.h"
!----------------------------------------------------------------------

!L 1. Initialise local constants
!
! Latitude of zeroth meridian
      LAMBDA_ZERO=LAMBDA_POLE+180.
! Sine and cosine of latitude of eq pole
      IF (PHI_POLE >= 0.0) THEN
        SIN_PHI_POLE =  SIN(PI_OVER_180*PHI_POLE)
        COS_PHI_POLE =  COS(PI_OVER_180*PHI_POLE)
      ELSE
        SIN_PHI_POLE = -SIN(PI_OVER_180*PHI_POLE)
        COS_PHI_POLE = -COS(PI_OVER_180*PHI_POLE)
      ENDIF

!L 2. Transform from equatorial to standard latitude-longitude

      DO 200 I= 1,POINTS

! Scale eq longitude to range -180 to +180 degs

      E_LAMBDA=LAMBDA_EQ(I)
      IF(E_LAMBDA >   180.0) E_LAMBDA=E_LAMBDA-360.0
      IF(E_LAMBDA <  -180.0) E_LAMBDA=E_LAMBDA+360.0

! Convert eq latitude & longitude to radians

      E_LAMBDA=PI_OVER_180*E_LAMBDA
      E_PHI=PI_OVER_180*PHI_EQ(I)

! Compute latitude using equation (4.7)

      ARG=COS_PHI_POLE*COS(E_LAMBDA)*COS(E_PHI)                         &
     &                   +SIN(E_PHI)*SIN_PHI_POLE
      ARG=MIN(ARG, 1.0)
      ARG=MAX(ARG,-1.0)
      A_PHI=ASIN(ARG)
      PHI(I)=RECIP_PI_OVER_180*A_PHI

! Compute longitude using equation (4.8)

      TERM1 =(COS(E_PHI)*COS(E_LAMBDA)*SIN_PHI_POLE                     &
     &       -SIN(E_PHI)*COS_PHI_POLE)
      TERM2=COS(A_PHI)
      IF(TERM2 <  SMALL) THEN
        A_LAMBDA=0.0
      ELSE
        ARG=TERM1/TERM2
        ARG=MIN(ARG, 1.0)
        ARG=MAX(ARG,-1.0)
        A_LAMBDA=RECIP_PI_OVER_180*ACOS(ARG)
        A_LAMBDA=SIGN(A_LAMBDA,E_LAMBDA)
        A_LAMBDA=A_LAMBDA+LAMBDA_ZERO
      END IF

! Scale longitude to range 0 to 360 degs

      IF(A_LAMBDA >= 360.0) A_LAMBDA=A_LAMBDA-360.0
      IF(A_LAMBDA <  0.0) A_LAMBDA=A_LAMBDA+360.0
      LAMBDA(I)=A_LAMBDA

200   CONTINUE

      RETURN
      END SUBROUTINE EQTOLL
#endif
