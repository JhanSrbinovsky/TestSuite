#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE LHN_SEARCH --------------------------------------------
!LL
!LL  Purpose : Search for suitable nearby latent heating profiles for
!LL            use in the LHN sceme.
!LL
!LL  For use on Cray C90
!LL  For Cray - Global  ; Enable defs GLOBAL
!LL
!LL             <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.4:
!LL version  Date
!LL  4.1  6/6/96 : Add "number of points searched" diagnostic (C Jones)
!LL  4.3  24/2/97 : T3E mods + remove timer call.  Stuart Bell
!    4.4 4/7/97 : Correct MPP version on LHN_SEARCH.  Deborah Salmond

!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL  Documentation :  FR  WP 171
!LL
!LLEND------------------------------------------------------------------
!L*  ARGUMENTS
      SUBROUTINE LHN_SEARCH(POINT,NEAR,RANGE                            &
     &                           ,SEARCH,IROWS,ICOLS,L_FIRST_SEARCH     &
     &                           ,TOT_PR,ANAL_PR,P_FIELD,RADIUS)

      IMPLICIT NONE

!-----AC common blocks
#include "acparm.h"
#include "comacp.h"

!-----DECLARE VARIABLES
      INTEGER       POINT,NEAR,P_FIELD,RANGE
      INTEGER       SEARCH(4*RANGE*(RANGE+1),2)
      INTEGER       IROWS, ICOLS
      INTEGER       RADIUS(5)
      REAL          TOT_PR(P_FIELD),ANAL_PR
      LOGICAL       L_FIRST_SEARCH

!-INTENT=IN---------------------------------------------------
!     P_FIELD                    - No. of points
!     POINT                      - Current point
!     RADIUS(5)                  - Diagnostic for breakdown of searches
!     RANGE                      - Search range in grid points
!     SEARCH                     - Template holding relative positions
!                                    of search points from POINT
!     TOT_PR(P_FIELD)            - Model total precip rates
!     ANAL_PR                    - Analysed ppn rate at POINT
!     IROWS                      - Number of rows in the model grid
!     ICOLS                      -    "    " columns  "    "     "
!-INTENT=INOUT-----------------------------------------------
!     L_FIRST_SEARCH             - True if this is the first search
!                                    this timestep.
!-INTENT=OUT-------------------------------------------------
!     NEAR                       - Nearest point with suitable profile
!*-----------------------------------------------------------

! Local arrays and variables
      INTEGER       XPOINT       ! X coord of POINT
      INTEGER       YPOINT       ! Y   "    "   "
      INTEGER       X            ! X coord of point to check
      INTEGER       Y            ! Y   "    "   "    "   "
      INTEGER       P            ! Number of point at (X,Y)
      INTEGER       JR , JN      ! Loop counters
      REAL          RATIO        ! ratio of model to obs ppn
      REAL          BEST         ! closest RATIO to 1
      LOGICAL       L_FOUND      ! false until a suitable point is found

!
!*----------------------------------------------------------------------
!
!
!C 1.0     Set up search template if this is the first time
!C         round this timestep.
!
      IF (L_FIRST_SEARCH) THEN
        DO JR = 1 , RANGE
          DO JN = 1 , 2*JR
            SEARCH(4*(JR-1)*JR+JN,1) = -JR-1+JN
            SEARCH(4*(JR-1)*JR+JN,2) = -JR
            SEARCH(4*(JR-1)*JR+JN+2*JR,1) = JR
            SEARCH(4*(JR-1)*JR+JN+2*JR,2) = -JR-1+JN
            SEARCH(4*(JR-1)*JR+JN+4*JR,1) = JR+1-JN
            SEARCH(4*(JR-1)*JR+JN+4*JR,2) = JR
            SEARCH(4*(JR-1)*JR+JN+6*JR,1) = -JR
            SEARCH(4*(JR-1)*JR+JN+6*JR,2) = JR+1-JN
          ENDDO     ! JN
        ENDDO       ! JR
        L_FIRST_SEARCH=.FALSE.
      ENDIF

!
!C 2.0     Loop round nearest points in ever increasing radius,
!C         looking for best fit to observed rain at POINT
!
!C 2.1     Initialise variables
!
      L_FOUND = .FALSE.
      BEST = 0.0
      RATIO = 0.0
      YPOINT = INT ((POINT-1)/ICOLS) + 1
      XPOINT = POINT - ( (YPOINT-1) * ICOLS )
!
!C 2.2     Loop over possible ranges
!
      DO JR = 1 , RANGE
        IF (.NOT. L_FOUND) THEN
!
!C  If not found one yet then loop over points at this range
!
          DO JN = 1 , 8*JR
!
!C  Make sure search point is within the grid.
!C    for global code this means P from 1 to number of pts,
!C    for limited area code, also need to make sure X within limits.
!
            X = XPOINT + SEARCH( 4*JR*(JR-1) +JN , 1)
            Y = YPOINT + SEARCH( 4*JR*(JR-1) +JN , 2)
#if defined(GLOBAL)
            P = (Y-1) * ICOLS + X
            IF ( P  >=  1 .AND. P  <=  (IROWS*ICOLS) ) THEN
#else
            IF ( X  >=  1 .AND. X  <=  ICOLS .AND.                      &
     &           Y  >=  1 .AND. Y  <=  IROWS )      THEN
              P = (Y-1) * ICOLS + X
#endif
! count points
              RADIUS(5)=RADIUS(5)+1
!
!C  Test model ppn at point, P
!
              IF ( TOT_PR(P)  >=  (EPSILON_LHN * ANAL_PR) ) THEN
                RATIO = TOT_PR(P) / ANAL_PR
                IF (RATIO  >   1.0) RATIO = 1.0/RATIO
!
!C  Keep record of best match at this range
!
                IF (RATIO  >   BEST) THEN
                  BEST    = RATIO
                  NEAR    = P
                  IF (.NOT.L_FOUND .AND. LHN_DIAG) THEN
                    IF (JR == 1) RADIUS(2)=RADIUS(2)+1
                    IF (JR == 2) RADIUS(3)=RADIUS(3)+1
                    IF (JR >= 3) RADIUS(4)=RADIUS(4)+1
                  ENDIF
                  L_FOUND = .TRUE.
                ENDIF       ! Ratio test
! Diagnostics of where pt is found
              ENDIF         ! Ppn rate test
            ENDIF           ! Within domain test
          ENDDO      ! JN
        ENDIF             ! Found one test
      ENDDO          ! JR
      RETURN
      END SUBROUTINE LHN_SEARCH
#endif
