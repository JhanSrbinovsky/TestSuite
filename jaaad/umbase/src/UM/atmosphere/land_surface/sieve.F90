#if defined(A08_7A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE SIEVE--------------------------------------------------
!LL
!LL  PURPOSE : TO CALCULATE THE THROUGHFALL OF WATER FALLING
!LL            THROUGH THE SURFACE CANOPY
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 5.2:
!LL VERSION  DATE
!LL  5.2   15/11/00   New Deck         M. Best
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
!LL  VERSION NO. 1 18/1/90
!LL
!LL  LOGICAL COMPONENTS COVERED: P252
!LL
!LL  SYSTEM TASK : P252
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER NO 25
!CL                  SECTION (3B(II)), EQN(P252.9)
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE SIEVE (                                                &
     & NPNTS,TILE_PTS,TILE_INDEX,AREA,CAN_CPY,R,FRAC,TIMESTEP,          &
     & CAN_WCNT,TOT_TFALL                                               &
     & )

      IMPLICIT NONE

      INTEGER                                                           &
     & NPNTS                                                            &
                            ! IN Total number of land points.
     &,TILE_PTS                                                         &
                            ! IN Number of tile points.
     &,TILE_INDEX(NPNTS)    ! IN Index of tile points.

      REAL                                                              &
     & AREA                                                             &
                            ! IN Fractional area of gridbox over which
                            !    water falls (%).
     &,CAN_CPY(NPNTS)                                                   &
                            ! IN Canopy capacity (kg/m2).
     &,R(NPNTS)                                                         &
                            ! IN Water fall rate (kg/m2/s).
     &,FRAC(NPNTS)                                                      &
                            ! IN Tile fraction.
     &,TIMESTEP             ! IN Timestep (s).

      REAL                                                              &
     & CAN_WCNT(NPNTS)                                                  &
                            ! INOUT Canopy water content (kg/m2).
     &,TOT_TFALL(NPNTS)     ! INOUT Cummulative canopy throughfall
!                           !       (kg/m2/s).

!  Workspace --------------------------------------------------------
      REAL                                                              &
     & AEXP                                                             &
                            ! Used in calculation of exponential
                            ! in throughfall formula.
     &,CAN_RATIO                                                        &
                            ! CAN_WCNT / CAN_CPY
     &,TFALL(NPNTS)                                                     &
                            ! Local throughfall (kg/m2/s).
     &,SMALLESTP            ! Smallest +ve real which can be represented

      INTEGER                                                           &
     & I                                                                &
                            ! Land point index.
     &,J                    ! Counter for loop over tile points.

      SMALLESTP = TINY(1.0)

      DO J=1,TILE_PTS
        I = TILE_INDEX(J)
        IF (CAN_CPY(I)  >   0.0 .AND. R(I)  >  SMALLESTP) THEN
           AEXP = AREA*CAN_CPY(I)/(R(I)*TIMESTEP)
           ! Only calculate if AEXP is small enough to avoid underflow
           IF (AEXP < -LOG(SMALLESTP)) THEN
              AEXP = EXP(-AEXP)
           ELSE
              AEXP = 0.0
           ENDIF
           CAN_RATIO = CAN_WCNT(I) / CAN_CPY(I)
           CAN_RATIO = MIN(CAN_RATIO,1.0)
           TFALL(I) = R(I) * ((1.0-CAN_RATIO)*AEXP + CAN_RATIO)
        ELSE
           TFALL(I) = R(I)
        END IF
        CAN_WCNT(I) = CAN_WCNT(I) + (R(I) - TFALL(I))*TIMESTEP
        TOT_TFALL(I) = TOT_TFALL(I) + FRAC(I)*TFALL(I)
      ENDDO

      RETURN
      END SUBROUTINE SIEVE
#endif
