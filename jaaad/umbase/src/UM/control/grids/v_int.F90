#if defined(C92_2A) || defined(MAKEBC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE V_INT-------------------------------------------------
!LL
!LL  Purpose:  Performs vertical interpolation from one arbitrary set
!LL            of pressure levels to another. The technique used is
!LL            linear interpolation in log(p). When interpolating
!LL            wind components there is an option (controlled by
!LL            MAX_WIND) for including data from max wind modelling.
!LL
!LL  Written by A. Dickinson
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   3.1  23/02/93    DO ALL directive inserted before loop 240
!LL                    Author: A. Dickinson    Reviewer: F. Rawlins
!LL
!LL   4.2  01/07/96    Revised for CRAY T3E. Vector gather replaced
!LL                    by algorithm in which each term in the
!LL                    interpolation formula is first collected into
!LL                    a separate array and the interpolation
!LL                    calculation carried out after the loop over level
!LL                    New arguments START and END introduced to
!LL                    facilitate the removal of duplicate calculations
!LL                    when using domain decomposition in MPP mode.
!LL                    Author: A. Dickinson    Reviewer: F. Rawlins
!LL   4.5  14/04/98    Use assumption that neighbouring points are
!LL                    likely to be on or near same level. Jump out
!LL                    of loop-over-levels once level found. Results
!LL                    in a 40 percent speedup on 19 levels for
!LL                    non-vector machines. S.D.Mullerworth
!LL  5.3  24/09/01  Portability changes.    Z. Gardner
!    5.5  17/04/03     Replace reference to obsolete section
!                      C92_1A with C92_2A. T.White
!    6.0  08/09/03     Added new def to allow use with makebc. R.Sempers
!LL
!LL
!LL Programming standard :
!LL
!LL Logical components covered : S111
!LL
!LL Project task :
!LL
!LL  Documentation: The interpolation formulae are described in
!LL                 unified model on-line documentation paper S1.
!LL
!LLEND -----------------------------------------------------------------
!
!*L  ARGUMENTS:-------------------------------------------------------
      SUBROUTINE V_INT(P_IN,P_OUT,DATA_IN,DATA_OUT,POINTS,LEVELS        &
     &               ,DATA_MAXW,P_MAXW,MAX_WIND                         &
     &               ,START,END)

      IMPLICIT NONE

      INTEGER                                                           &
     & POINTS                                                           &
              ! Number of points to be processed.
     &,LEVELS                                                           &
              ! Number of levels in source data.
     &,START                                                            &
              ! Start position at each level
     &,END    ! Last point to be processed at each level

      REAL                                                              &
     & P_IN(POINTS,LEVELS)                                              &
                             !IN 3-D field of pressures at which
                             ! source data is stored.
     &,P_OUT(POINTS)                                                    &
                             !IN Array of pressure values to be
                             ! interpolated to.
     &,DATA_IN(POINTS,LEVELS)                                           &
                             !IN Source data as 3-D field.
     &,DATA_OUT(POINTS)                                                 &
                             !OUT Result of interpolation.
     &,DATA_MAXW(POINTS)                                                &
                             !IN Max wind data.
     &,P_MAXW(POINTS)        !IN Pressure of max wind data.

      LOGICAL                                                           &
     & MAX_WIND !IN Switch to include max winds if required.

! Workspace usage:-----------------------------------------------------
      REAL                                                              &
     & P1(POINTS)                                                       &
                             ! Upper input pressure \                  .
     &,P2(POINTS)                                                       &
                             ! Lower input pressure  \ Used in interp-
     &,D1(POINTS)                                                       &
                             ! Upper input data      / olation formula
     &,D2(POINTS)            ! Lower input data     /
! External subroutines called:-----------------------------------------
! None
!*---------------------------------------------------------------------
! Define local variables:----------------------------------------------
      INTEGER I,J                                                       &
     &  ,LAST                   ! Stores level of preceding point
      REAL ALPHA
!----------------------------------------------------------------------

! Initialise LAST to any value between 1 and LEVELS
      LAST=2
      DO I=START,END

! Start from same level as last point. First check whether this point
! is above or below, then continue search in appropriate direction
        IF(P_OUT(I) >= P_IN(I,LAST))THEN

! These next two loops exit immediately once level found.
! GOTO cuts out needless looping once level is found, reducing the
! cost of the routine by about 40 percent for 19 level runs.
          DO J=LAST,2,-1
            IF(P_OUT(I) <  P_IN(I,J-1))THEN
              GOTO 240
            ENDIF
          ENDDO
        ELSE
          DO J=LAST+1,LEVELS
            IF(P_OUT(I) >= P_IN(I,J))THEN
              GOTO 240
            ENDIF
          ENDDO
        ENDIF
 240    CONTINUE

! At this point, J is:
!    1         for below bottom level.
!    LEVELS+1  for above top level
!    Otherwise J is the level just above the point

        IF (J >  1.AND.J <= LEVELS)THEN
! Between top and bottom level
          P1(I)=P_IN(I,J)
          P2(I)=P_IN(I,J-1)
          D1(I)=DATA_IN(I,J)
          D2(I)=DATA_IN(I,J-1)
          LAST=J
        ELSE
! Special case; above top or below bottom.
! Set output field to top/bottom-most input field
          IF(J == LEVELS+1)J=LEVELS
          P1(I)=P_OUT(I)
          P2(I)=1.0
          D1(I)=DATA_IN(I,J)
          D2(I)=0.0
          LAST=J
        ENDIF
      ENDDO ! DO I=START,END

! If there is an extra level of winds from max wind modelling, include
! these in the interpolation. Repeat the level-finding logic because
! there are no calls with MAX_WIND=.TRUE. in UM so do not want to slow
! down the above loop by including the MAX_WIND test in the above.

      IF (MAX_WIND)THEN
        DO I=START,END

! If max wind level between current levels, redo interpolation
! incorporating max wind info.

! Start from same level as last point. First check whether this point
! is above or below, then check all levels above/below in turn
          IF(P_OUT(I) >= P_IN(I,LAST))THEN
! Below LAST level.
! These loops exit immediately once level found.
! GOTO cuts out needless looping once level is found, reducing the
! cost of the routine by about 40 percent for 19 level runs.
            DO J=LAST,2,-1
              IF(P_OUT(I) <  P_IN(I,J-1))THEN
                GOTO 340
              ENDIF
            ENDDO
          ELSE
            DO J=LAST+1,LEVELS
              IF(P_OUT(I) >= P_IN(I,J))THEN
                GOTO 340
              ENDIF
            ENDDO
          ENDIF
 340      CONTINUE

          IF(J >  1.AND.J <= LEVELS)THEN
            IF(P_MAXW(I) <  P_IN(I,J-1).AND.P_MAXW(I) >= P_IN(I,J))THEN

              IF(P_OUT(I) <  P_MAXW(I))THEN

! (i)  p(maxwind) > p(out) >= p(j)

                P2(I)=P_MAXW(I)
                D2(I)=DATA_MAXW(I)

              ELSE

! (ii) p(j-1) > p(out) >= p(maxwind)

                P1(I)=P_MAXW(I)
                D1(I)=DATA_MAXW(I)

              ENDIF
            ENDIF
          ENDIF

        ENDDO                   ! DO I=START,END

      ENDIF

!L 3. Compute equation (3.3)

#if defined(VECTLIB)
      CALL ONEOVER_V(END-START+1,P2(START),P2(START))
      DO I=START,END
        P1(I)=P1(I)*P2(I)
        P2(I)=P_OUT(I)*P2(I)
      ENDDO
      CALL ALOG_V(END-START+1,P1(START),P1(START))
      CALL ONEOVER_V(END-START+1,P1(START),P1(START))
      CALL ALOG_V(END-START+1,P2(START),P2(START))
      DO I=START,END
        ALPHA=P1(I)*P2(I)
        DATA_OUT(I)=ALPHA*D1(I)+(1.-ALPHA)*D2(I)
      ENDDO
#else
! Compute alpha, the interpolation weight given by equation (3.4)
      DO I=START,END
          ALPHA=ALOG(P_OUT(I)/P2(I))                                    &
     &         /ALOG(P1(I)/P2(I))
! Then apply equation (3.3)
          DATA_OUT(I)=ALPHA*D1(I)+(1.-ALPHA)*D2(I)
      ENDDO
#endif

      RETURN
      END SUBROUTINE V_INT
#endif
