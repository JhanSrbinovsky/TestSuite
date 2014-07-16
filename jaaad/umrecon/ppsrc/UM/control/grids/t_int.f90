
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE T_INT:--------------------------------------------------
!LL
!LL  Purpose:
!LL       Carries out linear interpolation in time between two fields.
!LL       If the missing data indicator is present at one of the
!LL       times, the value at the other time is used.
!LL
!LL  Written by A. Dickinson 30/03/90
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL
!LL  Programming standard:
!LL           Unified Model Documentation Paper No 3
!LL           Version No 1 15/1/90
!LL
!LL  System component:S190
!LL
!LL  System task: S1
!LL
!LL  Documentation:
!LL       The interpolation formulae are described in
!LL       unified model on-line documentation paper S1.
!LL
!LL  -------------------------------------------------------------------
!*L  Arguments:---------------------------------------------------------

      SUBROUTINE T_INT(DATA_T1,T1,DATA_T2,T2,DATA_T3,T3,POINTS)

      IMPLICIT NONE

      INTEGER                                                           &
     & POINTS  !IN No of points to be processed

      REAL                                                              &
     & DATA_T1(POINTS)                                                  &
                       !IN Data at T1
     &,DATA_T2(POINTS)                                                  &
                       !IN Data at T2
     &,DATA_T3(POINTS)                                                  &
                       !OUT Data at T3
     &,T1                                                               &
          !IN Time of first data field
     &,T2                                                               &
          !IN Time of second data field
     &,T3 !IN Time at which new field is required T1<=T3<=T2


! Local arrays:---------------------------------------------------------
! None
! ----------------------------------------------------------------------
!*L External subroutines called:----------------------------------------
! None
!*----------------------------------------------------------------------
! Local variables:------------------------------------------------------
      REAL                                                              &
     & ALPHA !Fractional time

      INTEGER                                                           &
     & I     !Loop index
! ----------------------------------------------------------------------
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------


      ALPHA=(T3-T1)/(T2-T1)
      DO 100 I=1,POINTS
      DATA_T3(I)=DATA_T2(I)*ALPHA+DATA_T1(I)*(1-ALPHA)
      IF(DATA_T1(I) == RMDI)DATA_T3(I)=DATA_T2(I)
      IF(DATA_T2(I) == RMDI)DATA_T3(I)=DATA_T1(I)
100   CONTINUE

      RETURN
      END SUBROUTINE T_INT

