
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Performs Bi-linear horizitontal interpolation
!
! Subroutine Interface:
      SUBROUTINE H_INT_BL_O(ROWS_IN,ROW_LENGTH_IN,LEN_FIELD_OUT         &
     &,                   INDEX_B_L,INDEX_B_R,DATA_IN                   &
     &,                   WEIGHT_B_L,WEIGHT_B_R,WEIGHT_T_L,WEIGHT_T_R   &
     &,                   levs, i_field, levn                           &
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
!   and gather indices calculated in subroutine H_INT_CO. Interpolates
!   using values at ocean points (i.e. avoiding values at land points).
!
! Method:
!   See UMDP S1 for full desciption
!
! Current Code Owner: M. J. Bell
!
! History:
! Version   Date     Comment
! -------   ----     -------
!LL 5.2    30/07/00 New deck. Based on HINT_BL1 but taking into account
!LL                 land/sea mask to interpolate using only values at
!LL                 sea points. M. J. Bell
!LL  5.3  24/09/01  Portability changes.    Z. Gardner
!LL 5.3    13/11/01 Adjustments to take account of changes
!                   to interpolation indices following change
!                   from N -> S to S -> N atmosphere grid. M J Bell
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
      REAL     levs(ROW_LENGTH_IN*ROWS_IN)  ! IN levels data
      INTEGER  i_field   ! IN field type
      INTEGER  levn  ! IN level number

!   Array  arguments with intent(out):
      REAL     DATA_OUT(LEN_FIELD_OUT) !Data after interpolation

! Global parameter
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
! Local scalars:
      INTEGER      I
      integer  inw, ine, isw, ise  ! indices for each corner
      real     nw, ne, sw, se      ! values at each corner

! Function & Subroutine calls:
!     External None

!- End of header
!------------------------------------------------------------------
!     1. Carry out horizontal interpolation using equation (2.1)

      DO I=1,LEN_FIELD_OUT

        isw = INDEX_B_L(I)
        ise = INDEX_B_R(I)
        inw = INDEX_B_L(I)+ROW_LENGTH_IN
        ine = INDEX_B_R(I)+ROW_LENGTH_IN
        nw = DATA_IN(inw)
        ne = DATA_IN(ine)
        sw = DATA_IN(isw)
        se = DATA_IN(ise)

        IF ( ( levn  <=  levs(inw) .and. levn  <=  levs(ine) .and.      &
     &       levn  <=  levs(isw) .and. levn  <=  levs(ise) )            &
     &       .or. i_field  >=  2 ) then
! all values are sea points or field is
! either velocity or streamfn type; just interpolate
          DATA_OUT(I) = WEIGHT_B_L(I)*sw + WEIGHT_B_R(I)*se             &
     &                + WEIGHT_T_L(I)*nw + WEIGHT_T_R(I)*ne

        ELSE IF ( levn  >   levs(inw) .and. levn  >   levs(ine) .and.   &
     &            levn  >   levs(isw) .and. levn  >   levs(ise) ) then
! all tracer points are land points; set value to missing data

          DATA_OUT(I)= RMDI

        ELSE
! between 1 and 3 of surrounding tracer points are land

! work round all four points setting values from neighbouring ones

          IF ( levn  >   levs(inw) ) then
            IF ( levn  <=  levs(isw) .and. levn  <=  levs(ine) ) then
                nw = 0.5 * ( sw + ne )
            ELSE IF ( levn  <=  levs(isw) ) then
                nw = sw
            ELSE IF ( levn  <=  levs(ine) ) then
                nw = ne
            ELSE
                nw = RMDI
            END IF
          END IF

          IF ( levn  >   levs(isw) ) then
          IF ( levn  <=  levs(inw) .and. levn  <=  levs(ise) ) then
                sw = 0.5 * ( nw + se )
              ELSE IF ( levn  <=  levs(inw) ) then
                sw = nw
              ELSE IF ( levn  <=  levs(ise) ) then
                sw = se
              ELSE
                sw = RMDI
              END IF
          END IF

          IF ( levn  >   levs(ine) ) then
          IF ( levn  <=  levs(inw) .and. levn  <=  levs(ise) ) then
                ne = 0.5 * ( nw + se )
              ELSE IF ( levn  <=  levs(inw) ) then
                ne = nw
              ELSE IF ( levn  <=  levs(ise) ) then
                ne = se
              ELSE
                ne = RMDI
              END IF
          END IF

          IF ( levn  >   levs(ise) ) then
          IF ( levn  <=  levs(isw) .and. levn  <=  levs(ine) ) then
                se = 0.5 * ( sw + ne )
              ELSE IF ( levn  <=  levs(isw) ) then
                se = sw
              ELSE IF ( levn  <=  levs(ine) ) then
                se = ne
              ELSE
                se = RMDI
              END IF
          END IF

! finally fill in any point which is still not set (3 of 4 points
! will definitely have been set already)

          IF ( nw  ==  RMDI) nw = 0.5 * ( sw + ne )
          IF ( sw  ==  RMDI) sw = 0.5 * ( nw + se )
          IF ( ne  ==  RMDI) ne = 0.5 * ( nw + se )
          IF ( se  ==  RMDI) se = 0.5 * ( sw + ne )

          DATA_OUT(I) = WEIGHT_B_L(I)*sw + WEIGHT_B_R(I)*se             &
     &                + WEIGHT_T_L(I)*nw + WEIGHT_T_R(I)*ne

        END IF ! between 1 and 3 of surrounding points are land

      END DO ! I

      RETURN
      END SUBROUTINE H_INT_BL_O
