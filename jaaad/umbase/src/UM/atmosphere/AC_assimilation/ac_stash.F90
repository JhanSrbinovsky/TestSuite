#if defined(A18_2A) || defined(O35_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE AC_STASH -----------------------------------------------
!LL
!LL  Purpose : Stash processing for AC Scheme
!LL
!LL  For Cray - Global and Limited area
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  3.1  12/1/93 : Remove fix introduced in 2.8. Initialise
!LL               : N_LEVELS_LIST. Add extra checks. D Robinson.
!LL  3.2  8/7/93  :   Eliminate QA FORTRAN complaints    S Bell
!LL  5.2 20/8/00  :   Turn off STASH messages        A. Hines
!    6.0 11/09/03 :   Removed double ? for IBM cpp.  P.Dando
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LLEND------------------------------------------------------------------
      SUBROUTINE AC_STASH (STASH_ITEM_NO,LEVEL_NO,N_LEVELS,             &
     &                     GROUP_NO,N_GROUPS,TIMESTEP_NO,               &
     &                     STINDEX,STLIST,LEN_STLIST,SI,SF,             &
     &                     STASHWORK,STASH_LEVELS,NUM_STASH_LEVELS,     &
     &                     STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,        &
     &                     FIELD,LEN_FLD,LABEL,                         &
     &                     ICODE,CMESSAGE)

      IMPLICIT NONE

      INTEGER                                                           &
     &   STASH_ITEM_NO                                                  &
                              ! Stash Item Number
     &,  LEVEL_NO                                                       &
                              ! Level number
     &,  N_LEVELS                                                       &
                              ! No of levels
     &,  GROUP_NO                                                       &
                              ! Group number
     &,  N_GROUPS                                                       &
                              ! No of groups
     &,  TIMESTEP_NO                                                    &
                              ! Timestep number
     &,  LEN_FLD                                                        &
                              ! Length of field to be stashed
     &,  LEN_STLIST                                                     &
                              ! Dimension of STLIST
     &,  STINDEX(2,*)                                                   &
                              ! Start and no of items in STLIST
     &,  STLIST(LEN_STLIST,*)                                           &
                              ! Stash List of items to be output
     &,  SI(*)                                                          &
                              ! Address of Item in STASHWORK
     &,  NUM_STASH_LEVELS                                               &
                              ! Number of levels lists
     &,  NUM_STASH_PSEUDO                                               &
                              ! Number of pseudo lists
     &,  STASH_LEVELS(NUM_STASH_LEVELS+1,*)                             &
                                                   ! Levels lists
     &,  STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO+1,*) ! Pseudo lists

      REAL                                                              &
     &   STASHWORK(*)                                                   &
                              ! Work array for stashed data
     &,  FIELD(LEN_FLD)       ! Field to be stashed

      LOGICAL SF(*)           ! Stash Flags

      CHARACTER*(*) LABEL     !  Label to indicate field being stashed

      INTEGER ICODE           !  Return code
      CHARACTER*256 CMESSAGE  !  Error message

#include "stparam.h"

!*L   Dynamic allocated arrays

      LOGICAL LEVELS_LIST(N_LEVELS)  ! Expanded levels list
      LOGICAL PSEUDO_LIST(N_GROUPS)  ! Expanded pseudo list

!*L   External Subroutines called

      EXTERNAL SET_LEVELS_LIST, SET_PSEUDO_LIST

!     Local variables

      LOGICAL                                                           &
     &   L_SINGLE_LEV                                                   &
     &,  L_LEVELS_LIST                                                  &
     &,  L_PSEUDO_LIST                                                  &
     &,  LSTASH

      INTEGER                                                           &
     &   J,JLEV,JGRP                                                    &
                          !  Loop counters over levels/groups
     &,  J0                                                             &
                          !  Pointer in STASHWORK
     &,  IPOS                                                           &
                          !  Position in STLIST for this STASH_ITEM_NO
     &,  N_LEVELS_LIST                                                  &
                          !  No of levels in levels list
     &,  FLD_NO                                                         &
                          !  Field number in STASHWORK
     &,  GRP_NO           !  Group number in STASHWORK

!     Check that LEVEL_NO le N_LEVELS
      IF (LEVEL_NO >  N_LEVELS) THEN
        ICODE    = 1
        CMESSAGE = ' AC_STASH : LEVEL_NO gt N_LEVELS ?'
        WRITE (6,*) 'AC_STASH : LEVEL_NO must be LE to N_LEVELS'
        WRITE (6,*) 'Stash Item No = ',STASH_ITEM_NO
        WRITE (6,*) 'LEVEL_NO = ',LEVEL_NO,' N_LEVELS = ',N_LEVELS
        GO TO 999   !  Return
      ENDIF

!     Check that GROUP_NO le N_GROUPS
      IF (GROUP_NO >  N_GROUPS) THEN
        ICODE    = 1
        CMESSAGE = ' AC_STASH : GROUP_NO gt N_GROUPS ?'
        WRITE (6,*) 'AC_STASH : GROUP_NO must be LE to N_GROUPS'
        WRITE (6,*) 'Stash Item No = ',STASH_ITEM_NO
        WRITE (6,*) 'GROUP_NO = ',GROUP_NO,' N_GROUPS = ',N_GROUPS
        GO TO 999   !  Return
      ENDIF

      DO JLEV = 1,N_LEVELS
        LEVELS_LIST(JLEV) = .FALSE.
      ENDDO
      DO JGRP = 1,N_GROUPS
        PSEUDO_LIST(JGRP) = .FALSE.
      ENDDO

!     Get position in STLIST
      IPOS = STINDEX(1,STASH_ITEM_NO)

!     Determine if levels list used (Entry 10 in STLIST = negative)
      L_LEVELS_LIST = STLIST(ST_INPUT_BOTTOM,IPOS) <  0

!     Determine if single level (Entry 10 in STLIST = 100)
      L_SINGLE_LEV  = STLIST(ST_INPUT_BOTTOM,IPOS) == 100

!     Determine if pseudo list used (Entry 26 in STLIST = positive)
      L_PSEUDO_LIST = STLIST(ST_PSEUDO_IN,IPOS) >  0

      N_LEVELS_LIST = 1

      IF (L_LEVELS_LIST) THEN

        N_LEVELS_LIST = STASH_LEVELS(1,-STLIST(ST_INPUT_BOTTOM,IPOS))

!----   Get levels required for this field
!----   (Sets up LEVELS_LIST from STASH_LEVELS)
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST (N_LEVELS,LEN_STLIST,STLIST(1,IPOS),       &
     &       LEVELS_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,               &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) GO TO 999   !  Return

      ELSEIF (.NOT. L_SINGLE_LEV) THEN

        ICODE = 1
        CMESSAGE = 'AC_STASH ; No levels list ?'
        GO TO 999   !  Return

      ELSE
      ENDIF

      IF (L_PSEUDO_LIST) THEN

!----   Get pseudo levels list required for this field
!----   (Sets up PSEUDO_LIST from STASH_PSEUDO_LEVELS)
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST (N_GROUPS,LEN_STLIST,STLIST(1,IPOS),       &
     &       PSEUDO_LIST,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,          &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) GO TO 999   !  Return

      ENDIF

!     Determine if this field is to stashed
      LSTASH = .FALSE.
      IF ( L_SINGLE_LEV .OR.                                            &
     &    (L_LEVELS_LIST .AND. LEVELS_LIST(LEVEL_NO)) ) THEN
        IF (L_PSEUDO_LIST) THEN
          IF (PSEUDO_LIST(GROUP_NO)) LSTASH = .TRUE.
        ELSE
          LSTASH = .TRUE.
        ENDIF
      ENDIF

      IF (LSTASH) THEN

!----   Determine position in STASHWORK for this field
        FLD_NO = 0
        IF (L_SINGLE_LEV) THEN
          FLD_NO = 1
!         Could have pseudo levels - look at later
        ELSE
          DO JLEV=1,LEVEL_NO
            IF (LEVELS_LIST(JLEV)) THEN
              FLD_NO = FLD_NO+1
            ENDIF
          ENDDO
        ENDIF

        GRP_NO = 0
        IF (L_PSEUDO_LIST) THEN
          DO JLEV=1,GROUP_NO
            IF (PSEUDO_LIST(JLEV)) THEN
              GRP_NO = GRP_NO+1
            ENDIF
          ENDDO
        ELSE
          GRP_NO = 1
        ENDIF

!----   Check FLD_NO
        IF (FLD_NO == 0) THEN
          ICODE    = 1
          CMESSAGE = ' AC_STASH : FLD_NO = 0 ?'
          WRITE (6,*) 'AC_STASH : FLD_NO must be GT than 0'
          WRITE (6,*) 'Stash Item No = ',STASH_ITEM_NO
          WRITE (6,*) 'Level,Group,FLD_NO = ',LEVEL_NO,GROUP_NO,FLD_NO
          GO TO 999   !  Return
        ENDIF

!----   Check GRP_NO
        IF (GRP_NO == 0) THEN
          ICODE    = 1
          CMESSAGE = ' AC_STASH : GRP_NO = 0 ?'
          WRITE (6,*) 'AC_STASH : GRP_NO must be GT than 0'
          WRITE (6,*) 'Stash Item No = ',STASH_ITEM_NO
          WRITE (6,*) 'Level,Group,FLD_NO = ',LEVEL_NO,GROUP_NO,FLD_NO
          GO TO 999   !  Return
        ENDIF

!----   Set up pointer for this field in STASHWORK
        J0 = SI(STASH_ITEM_NO) - 1                                      &
     &     + (GRP_NO-1)*(N_LEVELS_LIST*LEN_FLD)                         &
     &     + (FLD_NO-1)*LEN_FLD

!----   Copy field into work space
        DO J=1,LEN_FLD
          STASHWORK(J0 + J) = FIELD (J)
        ENDDO


      ENDIF

 999  CONTINUE
      RETURN
      END SUBROUTINE AC_STASH
#endif
