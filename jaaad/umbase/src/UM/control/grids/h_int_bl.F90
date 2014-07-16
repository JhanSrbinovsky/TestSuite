#if defined(C92_2A) || defined(RECON) || defined(VAROPSVER) \
 || defined(MAKEBC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Performs Bi-linear horizitontal interpolation
!
! Subroutine Interface:
      SUBROUTINE H_INT_BL(ROWS_IN,ROW_LENGTH_IN,LEN_FIELD_OUT           &
     &,                   INDEX_B_L,INDEX_B_R,DATA_IN                   &
     &,                   WEIGHT_B_L,WEIGHT_B_R,WEIGHT_T_L,WEIGHT_T_R   &
     &,                   DATA_OUT)

!LL  System component: S121
!LL
!LL  System task: S1
!LL
!LL  Purpose:
!LL
!LL  Documentation:
!LL            The interpolation formulae are described in
!LL            unified model on-line documentation paper S1.
!LL
      IMPLICIT NONE
!
! Description:
!   Carries out bi-linear horizontal interpolation using coefficients
!   and gather indices calculated in subroutine H_INT_CO
!
! Method:
!   See UMDP S1 for full desciption
!
! Current Code Owner: D.M. Goddard
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 3.0      ??????   Original code. A.Dickinson
! 4.0      160395   Renamed and brought up to new standard D.M. Goddard
!LL  5.1     10/04/00   South->North ordering. P.Selwood.
!LL  5.5     17/02/03 Include missing data for Wave model use
!LL                                             D.Holmes-Bell
!    6.0     05/09/03  Added new def for use with makebc. R.Sempers
!    6.2     23/11/05  Removed all references to the wavemodel.
!                      T.Edwards
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered: S121
! System Task:              S1
!
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER  ROWS_IN              !No of P rows on source grid
      INTEGER  ROW_LENGTH_IN        !No of pts per row on source grid
      INTEGER  LEN_FIELD_OUT        !No of points on target grid

!   Array  arguments with intent(in):
      INTEGER  INDEX_B_L(LEN_FIELD_OUT)
                                     !Index of bottom lefthand corner
                                     !  of source gridbox
      INTEGER  INDEX_B_R(LEN_FIELD_OUT)
                                     !Index of bottom righthand corner
                                     !  of source gridbox
      REAL     DATA_IN(ROWS_IN*ROW_LENGTH_IN)
                                      !Data before interpolation
      REAL     WEIGHT_B_L(LEN_FIELD_OUT)
                                     !Weight applied to value at bottom
                                     !lefthand corner of source gridbox
      REAL     WEIGHT_B_R(LEN_FIELD_OUT)
                                     !Weight applied to value at bottom
                                     !righthand corner of source gridbox
      REAL     WEIGHT_T_L(LEN_FIELD_OUT)
                                     !Weight applied to value at top
                                     !lefthand corner of source gridbox
      REAL     WEIGHT_T_R(LEN_FIELD_OUT)
                                     !Weight applied to value at top
                                     !righthand corner of source gridbox

!   Array  arguments with intent(out):
      REAL     DATA_OUT(LEN_FIELD_OUT) !Data after interpolation

! Local scalars:
      INTEGER      I
      Real, Parameter :: W_min_in = 1.0E-19
      Real, Parameter :: W_min_ou = 1.0E-5

#include "c_mdi.h"
! Function & Subroutine calls:
!     External None

!- End of header

!     1. Carry out horizontal interpolation using equation (2.1)

      DO I=1,LEN_FIELD_OUT

      DATA_OUT(I)=WEIGHT_B_L(I)*DATA_IN(INDEX_B_L(I))                   &
     &           +WEIGHT_B_R(I)*DATA_IN(INDEX_B_R(I))                   &
     &           +WEIGHT_T_L(I)*DATA_IN(INDEX_B_L(I)+ROW_LENGTH_IN)     &
     &           +WEIGHT_T_R(I)*DATA_IN(INDEX_B_R(I)+ROW_LENGTH_IN)

      END DO

      RETURN
      END SUBROUTINE H_INT_BL

#endif
