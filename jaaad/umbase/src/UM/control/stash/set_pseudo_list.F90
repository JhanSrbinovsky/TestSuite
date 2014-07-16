#if defined(C84_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine SET_PSEUDO_LIST -----------------------------------------
!LL
!LL Purpose : To set up a list of pseudo levels at which a diagnostic
!LL           is required, using information in the STASH list.
!LL
!LL Service routine
!LL
!LL Written by D Robinson 7/10/92
!LL Copy of Subroutine SET_LEVELS_LIST (Deck SETLST1) taken and
!LL adapted for pseudo levels.
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL 3.1    12/1/93   : More error checking . Dave Robinson
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL
!LL Programming Standard : Unified Model Documentation paper number 3
!LL                      : Version no 4, dated 5/2/92
!LL
!LL Logical components covered : D3
!LL
!LL System task : P0
!LL
!LL Documentation: U.M. Documentation paper number P0,C4
!LL
!LLEND -----------------------------------------------------------------
!
!*L Arguments

      SUBROUTINE SET_PSEUDO_LIST                                        &
     &      (N_LEVELS,LEN_STLIST,STLIST,PSEUDO_LIST,                    &
     &      STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,ICODE,CMESSAGE)

      IMPLICIT NONE

      INTEGER                                                           &
     &   N_LEVELS                                                       &
                             ! IN Number of possible pseudo levels
     &,  LEN_STLIST                                                     &
                             ! IN Dimension of STLIST
     &,  STLIST(LEN_STLIST)                                             &
                             ! IN STASH list
     &,  NUM_STASH_PSEUDO                                               &
                             ! IN Dimension for STASH_PSEUDO_LEVELS
     &,  STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO+1,*)                      &
                                                   ! IN Pseudo levels
     &,  ICODE               ! OUT Return code

      LOGICAL                                                           &
     &   PSEUDO_LIST(N_LEVELS) ! OUT List of pseudo levels required.

      CHARACTER*(80) CMESSAGE ! Error message

!* ---------------------------------------------------------------------

!*L   Workspace Usage :- None

!*L   External Subroutines called :- None

#include "stparam.h"
!* ---------------------------------------------------------------------
!  Local variables

      INTEGER                                                           &
     &      JLEV                                                        &
                       ! Loop counter over levels
     &,     LEVEL_NO                                                    &
                       ! Level no in pseudo list
     &,     LIST_NO    ! Pseudo level list number

!* ---------------------------------------------------------------------

!L Initialise pseudo levels list to false

      DO JLEV=1,N_LEVELS
        PSEUDO_LIST(JLEV)= .FALSE.
      END DO

!L Get pseudo list number

      LIST_NO = STLIST(ST_PSEUDO_IN)

!L Check that Pseudo list number is valid (should be GE 0)

      IF (LIST_NO <  0) THEN

!       Illegal control data

        ICODE=1
        CMESSAGE = 'SET_PSEUDO_LIST: Illegal control data'
        WRITE(6,*) 'SET_PSEUDO_LIST: Illegal control data'
        WRITE(6,*) 'ST_PSEUDO_IN         = ',ST_PSEUDO_IN
        WRITE(6,*) 'STLIST(ST_PSEUDO_IN) = ',STLIST(ST_PSEUDO_IN)
        WRITE(6,*) 'Section and item numbers ',STLIST(2),STLIST(1)
        GO TO 999  !  Return

      ENDIF

!  Set logical array list to identify pseudo levels required.

      IF (LIST_NO >  0) THEN

        DO JLEV=2,STASH_PSEUDO_LEVELS(1,LIST_NO)+1
          LEVEL_NO = STASH_PSEUDO_LEVELS(JLEV,LIST_NO)
          IF (LEVEL_NO >= 1 .AND. LEVEL_NO <= N_LEVELS) THEN

!           Level is within range
            PSEUDO_LIST(LEVEL_NO) =.TRUE.

          ELSE

!           Level is out of range
            ICODE=2
            CMESSAGE=  ' SET_PSEUDO_LIST : level out of range'
            WRITE(6,*) ' SET_PSEUDO_LIST : level out of range'
            WRITE(6,*) ' pseudo list no = ',LIST_NO
            WRITE(6,*) ' level = ',LEVEL_NO
            WRITE(6,*) ' Section, Item = ',STLIST(2),STLIST(1)
            GO TO 999   !  Return

          END IF
        END DO

      END IF

 999  RETURN
      END SUBROUTINE SET_PSEUDO_LIST
#endif
