#if defined(A08_7A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!LL  SUBROUTINE FRUNOFF------------------------------------------------
!LL
!LL  PURPOSE : TO CALCULATE SURFACE RUNOFF
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 5.2:
!LL VERSION  DATE
!!!  5.2   15/11/00   New Deck         M. Best
!    6.0   11/09/03   Weaken test on R(I) to avoid rounding problems
!                     on IBM.                                P.Dando
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
!LL
!LL  LOGICAL COMPONENTS : P252
!LL
!LL  SYSTEM TASK :
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER NO 25
!LL                  SECTION (3B(II)), EQN(P252.14)
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE FRUNOFF (                                              &
     & NPNTS,TILE_PTS,TILE_INDEX,AREA,                                  &
     & CAN_CPY,CAN_WCNT,INFIL,R,FRAC,TIMESTEP,                          &
     & SURF_ROFF                                                        &
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
     &,CAN_WCNT(NPNTS)                                                  &
                            ! IN Canopy water content (kg/m2).
     &,INFIL(NPNTS)                                                     &
                            ! IN Infiltration rate (kg/m2/s).
     &,R(NPNTS)                                                         &
                            ! IN Water fall rate (kg/m2/s).
     &,FRAC(NPNTS)                                                      &
                            ! IN Tile fraction.
     &,TIMESTEP             ! IN Timestep (s).

      REAL                                                              &
     & SURF_ROFF(NPNTS)     ! OUT Cummulative surface runoff (kg/m2/s).

!  Workspace --------------------------------------------------------
      REAL                                                              &
     & AEXP                                                             &
                            ! Used in the calculation of exponential
     &,AEXP1                                                            &
                            ! terms in the surface runoff formula.
     &,AEXP2                                                            &
                            !
     &,CM                                                               &
                            ! (CAN_CPY - CAN_WCNT)/TIMESTEP
     &,CAN_RATIO                                                        &
                            ! CAN_WCNT / CAN_CPY
     &,RUNOFF                                                           &
                            ! Local runoff.
     &,SMALLESTP            ! Smallest +ve number that can be represented

      INTEGER                                                           &
     & I                                                                &
                            ! Land point index.
     &,J                    ! Counter for loop over tile points.

      SMALLESTP = TINY(1.0)

!CDIR NODEP
      DO J=1,TILE_PTS
        I = TILE_INDEX(J)
        RUNOFF = 0.
#if defined(IBM)
        IF (R(I) >  0.00000000001) THEN
#else
        IF (R(I) >  0.0) THEN
#endif
          IF ( INFIL(I)*TIMESTEP <= CAN_WCNT(I)                         &
     &                                   .AND. CAN_CPY(I) >  0.0 ) THEN
! Infiltration in timestep < or = canopy water content
             AEXP = AREA*CAN_CPY(I)/R(I)
             IF (CAN_WCNT(I)  >   0.0) THEN
               AEXP1 = EXP( -AEXP*INFIL(I)/CAN_WCNT(I))
             ELSE
               AEXP1 = 0.0
             END IF
             AEXP2 = EXP( -AEXP/TIMESTEP)
             CAN_RATIO = CAN_WCNT(I)/CAN_CPY(I)
             CAN_RATIO = MIN(CAN_RATIO,1.0)
             RUNOFF = R(I) * ( CAN_RATIO*AEXP1 +                        &
     &                                       (1. - CAN_RATIO)*AEXP2 )
!                                                        ... P252.14A
          ELSE
! Infiltration in timestep > canopy water content
             CM = (CAN_CPY(I)-CAN_WCNT(I))/TIMESTEP
             CM = MAX(CM,0.0)
             ! Only compute AEXP if will not generate an underflow error
             IF ( AREA*(INFIL(I)+CM)/R(I) < - LOG(SMALLESTP) ) THEN
                AEXP = EXP( -AREA*(INFIL(I)+CM)/R(I))
             ELSE
                AEXP = 0.0
             ENDIF
             RUNOFF = R(I)*AEXP                    !     ... P252.14B
          ENDIF
        ENDIF
        SURF_ROFF(I) = SURF_ROFF(I) + FRAC(I)*RUNOFF
      ENDDO

      RETURN
      END SUBROUTINE FRUNOFF
#endif
