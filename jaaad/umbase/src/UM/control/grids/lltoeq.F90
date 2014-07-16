#if defined(C92_2A) || defined(RECON) || defined(VAROPSVER) \
 || defined(MAKEBC) || defined(FRAMES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine LLTOEQ-------------------------------------------------
!LL
!LL  Purpose:  Calculates latitude and longitude on equatorial
!LL            latitude-longitude (eq) grid used in regional
!LL            models from input arrays of latitude and
!LL            longitude on standard grid. Both input and output
!LL            latitudes and longitudes are in degrees.
!LL
!LL  Written by A. Dickinson
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL    4.5   18/08/98  Added DEF,FLDOP   (A Van der Wal)
!LL    5.4   11/09/02  Cater for negative pole latitude. D.Robinson
!      6.0   05/09/03  Added def for enabling use in makebc. R.Sempers
!      6.1   19/07/04  Cater for Lamda_pole greater 180 G.Greed
!LL
!LL Programming standard :
!LL
!LL Logical components covered : S132
!LL
!LL Project task :
!LL
!LL  Documentation: The transformation formulae are described in
!LL                 unified model on-line documentation paper S1.
!LLEND -----------------------------------------------------------------
!
!*L  Arguments:--------------------------------------------------------
      SUBROUTINE LLTOEQ                                                 &
     &(PHI,LAMBDA,PHI_EQ,LAMBDA_EQ,PHI_POLE,LAMBDA_POLE,POINTS)

      IMPLICIT NONE

      INTEGER                                                           &
     & POINTS            !IN  Number of points to be processed

      REAL                                                              &
     & PHI(POINTS)                                                      &
                         !IN  Latitude
     &,LAMBDA(POINTS)                                                   &
                         !IN  Longitude
     &,LAMBDA_EQ(POINTS)                                                &
                         !OUT Longitude in equatorial lat-lon coords
     &,PHI_EQ(POINTS)                                                   &
                         !OUT Latitude in equatorial lat-lon coords
     &,PHI_POLE                                                         &
                         !IN  Latitude of equatorial lat-lon pole
     &,LAMBDA_POLE       !INOUT  Longitude of equatorial lat-lon pole
! Workspace usage:-----------------------------------------------------
! None
! ---------------------------------------------------------------------
! External subroutines called:-----------------------------------------
! None
!*---------------------------------------------------------------------
! Define local varables:-----------------------------------------------
      REAL A_LAMBDA,A_PHI,E_LAMBDA,ARG,E_PHI,SIN_PHI_POLE,COS_PHI_POLE
      REAL TERM1,TERM2,SMALL,LAMBDA_ZERO,LAMBDA_POLE_KEEP
      INTEGER I
      PARAMETER(SMALL=1.0E-6)
! ---------------------------------------------------------------------
! Constants from comdecks:---------------------------------------------
#include "c_pi.h"
! ---------------------------------------------------------------------

!L 1. Initialise local constants
! Scale lambda pole to range -180 to 180 degs
      LAMBDA_POLE_KEEP=LAMBDA_POLE
      IF (LAMBDA_POLE >   180.0) then
          LAMBDA_POLE=LAMBDA_POLE-360.0
      ENDIF
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

!L 2. Transform from standard to equatorial latitude-longitude

      DO 200 I= 1,POINTS

! Scale longitude to range -180 to +180 degs

      A_LAMBDA=LAMBDA(I)-LAMBDA_ZERO
      IF(A_LAMBDA >   180.0)A_LAMBDA=A_LAMBDA-360.
      IF(A_LAMBDA <= -180.0)A_LAMBDA=A_LAMBDA+360.

! Convert latitude & longitude to radians

      A_LAMBDA=PI_OVER_180*A_LAMBDA
      A_PHI=PI_OVER_180*PHI(I)

! Compute eq latitude using equation (4.4)

      ARG=-COS_PHI_POLE*COS(A_LAMBDA)*COS(A_PHI)                        &
     &                   +SIN(A_PHI)*SIN_PHI_POLE
      ARG=MIN(ARG, 1.0)
      ARG=MAX(ARG,-1.0)
      E_PHI=ASIN(ARG)
      PHI_EQ(I)=RECIP_PI_OVER_180*E_PHI

! Compute eq longitude using equation (4.6)

      TERM1 =(COS(A_PHI)*COS(A_LAMBDA)*SIN_PHI_POLE                     &
     &       +SIN(A_PHI)*COS_PHI_POLE)
      TERM2=COS(E_PHI)
      IF(TERM2 <  SMALL) THEN
        E_LAMBDA=0.0
      ELSE
        ARG=TERM1/TERM2
        ARG=MIN(ARG, 1.0)
        ARG=MAX(ARG,-1.0)
        E_LAMBDA=RECIP_PI_OVER_180*ACOS(ARG)
        E_LAMBDA=SIGN(E_LAMBDA,A_LAMBDA)
      ENDIF

! Scale longitude to range 0 to 360 degs

      IF(E_LAMBDA >= 360.0) E_LAMBDA=E_LAMBDA-360.0
      IF(E_LAMBDA <  0.0) E_LAMBDA=E_LAMBDA+360.0
      LAMBDA_EQ(I)=E_LAMBDA

200   CONTINUE

! Reset Lambda pole to the setting on entry to subroutine
      LAMBDA_POLE=LAMBDA_POLE_KEEP

      RETURN
      END SUBROUTINE LLTOEQ
#endif
