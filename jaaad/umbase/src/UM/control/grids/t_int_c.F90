#if defined(C92_2A) || defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE T_INT_C-------------------------------------------------
!LL
!LL  Purpose:
!LL    Carries out linear interpolation in time between two fields at
!LL    times T1 and T2. If the missing data indicator is present at one
!LL    of the times, the value at the other time is used. The interpolat
!LL    is controlled by a field ZI. A prescribed value is inserted where
!LL    If ZI changes between 0 and non-zero in the period T1 - T2, then
!LL    the field is linearly interpolated between its value at the time
!LL    ZI is non-zero and the prescibed value at the time when ZI become
!LL    The fractional time at which ZI changes between 0 and non-zero in
!LL    period T1 - T2 must be provided as input.
!LL
!LL  Written by A. Dickinson 30/03/90
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   3.2  17/03/93 : Correct for rounding problem, ie case where alpha
!LL                   should be exactly equal to frac_time.
!LL                   Author: R.A Stratton
!LL
!LL  Programming standard: Unified Model Documentation Paper No 3
!LL                        Version No 1 15/1/90
!LL
!LL  System component:S191
!LL
!LL  System task: S1
!LL
!LL  Documentation:  The interpolation formulae are described in
!LL            unified model on-line documentation paper S1.
!LL
!LL  -------------------------------------------------------------------
!*L  Arguments:---------------------------------------------------------

      SUBROUTINE T_INT_C(DATA_T1,T1,DATA_T2,T2,DATA_T3,T3,POINTS        &
     &,FRAC_TIME,ZI_T1,PRES_VALUE)

      IMPLICIT NONE

      INTEGER                                                           &
     & POINTS             !IN No of points to be processed

      REAL                                                              &
     & DATA_T1(POINTS)                                                  &
                          !IN Data at T1
     &,DATA_T2(POINTS)                                                  &
                          !IN Data at T2
     &,DATA_T3(POINTS)                                                  &
                          !OUT D_ta at T3
     &,ZI_T1(POINTS)                                                    &
                          !IN Value of controlling fieled at T1
     &,PRES_VALUE(POINTS)                                               &
                          !IN Prescribed value of Data when ZI=0
     &,FRAC_TIME(POINTS)                                                &
                          !IN Fractional time at which ZI changes betwee
                          !zero and non-zero in this time range
     &,T1                                                               &
          !IN Time of first data field
     &,T2                                                               &
          !In Time of second data field
     &,T3 !IN Time at which new field is required T1<=T3<=T2


! Local arrays:---------------------------------------------------------
       REAL INT_TIME(POINTS)
! ----------------------------------------------------------------------
!*L External subroutines called:----------------------------------------
      EXTERNAL T_INT
!*----------------------------------------------------------------------
! Local variables:------------------------------------------------------
      REAL                                                              &
     & ALPHA                                                            &
             !Fractional time
     &,EPSILON                                                          &
                   ! rounding error
     &,ALPHA_PLUS                                                       &
                   ! add rounding error to alpha
     &,ALPHA_MINUS ! alpha minus rounding error

      INTEGER                                                           &
     & I     !Loop index
! ----------------------------------------------------------------------
#include "c_mdi.h"

! set rounding error
      EPSILON=1.0E-6

! DEPENDS ON: t_int
      CALL T_INT(DATA_T1,T1,DATA_T2,T2,DATA_T3,T3,POINTS)

      ALPHA=(T3-T1)/(T2-T1)
      ALPHA_PLUS=ALPHA+EPSILON
      ALPHA_MINUS=ALPHA-EPSILON

      DO 100 I=1,POINTS

      IF(FRAC_TIME(I) /= RMDI)THEN
        IF(ZI_T1(I) == 0.0)THEN
           IF(ALPHA_MINUS <  FRAC_TIME(I))THEN
             DATA_T3(I)=PRES_VALUE(I)
           ELSE
             INT_TIME(I)=(ALPHA-FRAC_TIME(I))/(1.-FRAC_TIME(I))
             DATA_T3(I)=PRES_VALUE(I)*(1.-INT_TIME(I))                  &
     &                 +DATA_T2(I)*INT_TIME(I)
           ENDIF
        ELSE
           IF(ALPHA_PLUS >  FRAC_TIME(I))THEN
             DATA_T3(I)=PRES_VALUE(I)
           ELSE
             INT_TIME(I)=(FRAC_TIME(I)-ALPHA)/(FRAC_TIME(I))
             DATA_T3(I)=PRES_VALUE(I)*(1.-INT_TIME(I))                  &
     &                 +DATA_T1(I)*INT_TIME(I)
           ENDIF
        ENDIF
      ENDIF

100   CONTINUE




      RETURN
      END SUBROUTINE T_INT_C

#endif
