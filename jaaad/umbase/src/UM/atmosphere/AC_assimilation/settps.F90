#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE SETTPS--------------------------------------------------
!LL
!LL  Purpose : Sets up the list of AC Observation Types in the
!LL            order they are to be processed in the assimilation.
!LL            This routine is called each time the AC Observation
!LL            files are read in.
!LL
!LL  For use on Cray
!LL
!LL D.Robinson  <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   4.1    9/5/96   Clearer error processing S Bell
!    4.2 25/11/96:T3E mods Stuart Bell
!    6.0 11/09/03: Removed double ? for IBM cpp.         P.Dando
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LLEND------------------------------------------------------------------
!*L  Arguments----------------------------------------------------------
      SUBROUTINE SETTPS (ICODE,CMESSAGE)

!     Disable multi-tasking
!FPP$ NOCONCUR R

      IMPLICIT NONE
      INTEGER       ICODE
      CHARACTER*256 CMESSAGE

!    INTENT=OUT--------------------------------------------------------
!     ICODE        : Return Code
!     CMESSAGE     : Reason for failure
!*   ------------------------------------------------------------------

!     The variables/arrays (all in comdeck COMACP) set up are :-

!     NACT  : Number of AC Observation Types to be processed.
!     LACT  : List of the Obs Types in the order to be processed.
!     N_GROUPS : Number of groups for processing in AC.
!     GROUP_NO : Group in which each type is to be processed.

!     The order of processing is controlled by the array AC_ORDER.

!     The array AC_OBS_TYPES in the ACP namelist is used to control
!     which observation types are to be processed.

!-----AC common blocks
#include "acparm.h"
#include "comacp.h"
#include "comag.h"
#if defined(MPP)
#include "parvars.h"
#else
      INTEGER  mype
      PARAMETER (mype = 0 ) ! always zero in non-MPP code
#endif

!*L   Workspace Usage:--------------------------------------------------
!     None
!*    ------------------------------------------------------------------

!*L   External Subroutine Calls:----------------------------------------
!     None
!*    ------------------------------------------------------------------

!     Local variables

      INTEGER THIS_TYPE,J,JOBT,JTYPE
      INTEGER LAST_GROUP,THIS_GROUP
      LOGICAL LUSE

!     THIS_TYPE  : Obs type in AC_ORDER
!     J          : Loop counter for obs types in LACT
!     JOBT       : Loop counter for obs types in AC_ORDER
!     JTYPE      : Loop counter for obs types in AC_ORDER
!     THIS_GROUP : Indicator from AC_ORDER of grouping (current type)
!     LAST_GROUP : Indicator from AC_ORDER of grouping (previous type)

!     ------------------------------------------------------------------
!L    1. Initialise arrays and variables set up by SETTPS
!     ------------------------------------------------------------------
      NACT  = 0
      N_GROUPS = 0
      LAST_GROUP =0
      DO JOBT=1,NOBTYPMX
        LACT (JOBT) = 0
        GROUP_NO(JOBT) = 0
        GROUP_FIRST(JOBT) = 0
        GROUP_LAST(JOBT) = 0
      ENDDO
!     ----------------------------------------------------------------
!L    2. Set up order of processing in LACT
!     -----------------------------------------------------------------
!     Loop over all AC Obs types known to AC Scheme

      DO JTYPE=1,NOBTYPMX
      THIS_TYPE  = MOD(DEF_AC_ORDER(JTYPE),1000)
      THIS_GROUP = (DEF_AC_ORDER(JTYPE)-THIS_TYPE)/1000

!     This loop determines whether the observation type - THIS_TYPE -
!     is to be used or not from the AC_OBS_TYPES array.

      IF (THIS_TYPE >  0) THEN

!       Use observation type if in namelist array AC_OBS_TYPES

        LUSE = .FALSE.
        DO JOBT=1,NOBTYPMX
          IF (THIS_TYPE == AC_OBS_TYPES(JOBT)) THEN
            LUSE = .TRUE.
          ENDIF
        ENDDO

        IF (LUSE) THEN

!         Set up to process this observation type
          NACT = NACT+1
          LACT(NACT) = THIS_TYPE

!         Group observation types ; Set up GROUP_NO and N_GROUPS

          IF (NACT == 1 .OR. THIS_GROUP /= LAST_GROUP) THEN

!           Start a new group.
            N_GROUPS = N_GROUPS+1
            GROUP_INDEX(N_GROUPS) = THIS_GROUP
            GROUP_FIRST(N_GROUPS) = NACT

          ENDIF
          LAST_GROUP = THIS_GROUP
          GROUP_NO(NACT) = N_GROUPS
          GROUP_LAST(N_GROUPS) = NACT

!         Find this type in MASTER_AC_TYPES ; Abort if not found.
          TYPE_INDEX(NACT)=0
          DO JOBT=1,NOBTYPMX
            IF (THIS_TYPE == MASTER_AC_TYPES(JOBT)) THEN
              TYPE_INDEX(NACT)=JOBT
            ENDIF
          ENDDO
          IF (TYPE_INDEX(NACT) == 0) THEN
            ICODE = 1
            CMESSAGE = 'SETTPS : Observation Type not in Master List ?'
            if(mype == 0)                                               &
     &      PRINT *, ' Observation Type ',THIS_TYPE,                    &
     &                  ' not in Master List ?'
            GO TO 999
          ENDIF

        ENDIF

      ENDIF
      ENDDO   !   End of JTYPE loop.

!     -----------------------------------------------------------------
!     3. Print out list of AC Obs types to be processed
!     -----------------------------------------------------------------

      IF (NACT >  0.AND.mype == 0) THEN

        PRINT '(/,A,/)', ' AC Obs Types to be processed this run'
        PRINT '(A,(T12,13I5))', ' Type  No ',                           &
     &        (LACT(J),J=1,NACT)
        PRINT '(A,(T12,13I5))', ' Group No ',                           &
     &        (GROUP_NO(J),J=1,NACT)
!       PRINT '(A,15I5)', ' Position in Obs Type List    ',
!    +  (TYPE_INDEX(J),J=1,NACT)

        PRINT '(/,A,15I5)', ' Group Number                ',            &
     &  (J,J=1,N_GROUPS)
        PRINT '(A,15I5)', ' Group Index                 ',              &
     &  (GROUP_INDEX(J),J=1,N_GROUPS)
        PRINT '(A,15I5)', ' First Type in Group         ',              &
     &  (GROUP_FIRST(J),J=1,N_GROUPS)
        PRINT '(A,15I5)', ' Last Type in Group          ',              &
     &  (GROUP_LAST (J),J=1,N_GROUPS)
        PRINT '(A,15I5)', ' No of iterations            ',              &
     &  (DEF_NO_ITERATIONS(GROUP_INDEX(J)),J=1,N_GROUPS)
        PRINT '(A,15I5)', ' Interval between Iterations ',              &
     &  (DEF_INTERVAL_ITER(GROUP_INDEX(J)),J=1,N_GROUPS)
        PRINT '(A,15I5)', ' Ratio of MG Rows to AG Rows ',              &
     &  (DEF_AGRES_ROWS(GROUP_INDEX(J)),J=1,N_GROUPS)
        PRINT '(A,15I5)', ' Ratio of MG Pts  to AG Pts  ',              &
     &  (DEF_AGRES_PTS(GROUP_INDEX(J)),J=1,N_GROUPS)
        PRINT '(A,15I5)', ' No of analysis levels       ',              &
     &  (DEF_NO_ANAL_LEVS(GROUP_INDEX(J)),J=1,N_GROUPS)
        PRINT '(A,15I5)', ' No of weight levels         ',              &
     &  (DEF_NO_WT_LEVS(GROUP_INDEX(J)),J=1,N_GROUPS)
        PRINT '(A,15I5)', ' Horizontal Analysis Mode    ',              &
     &  (DEF_MODE_HANAL(GROUP_INDEX(J)),J=1,N_GROUPS)

        PRINT '(/,A,15I5)', ' Group Dep scaling FACTORS in FI  '
        PRINT '(A,15I11)',' Group No ',(J,J=1,N_GROUPS)
        PRINT '(A,15E11.4)', '          ',                              &
     &  (DEF_FI_VAR_FACTOR(GROUP_INDEX(J)),J=1,N_GROUPS)

        PRINT '(/,A,15I5)', ' Nudging Coefficients '
        PRINT '(A,15I11)',' Group No ',(J,J=1,N_GROUPS)
#if defined(GLOBAL)
        PRINT '(A,15E11.4)', ' NH       ',                              &
     &  (DEF_NUDGE_NH(GROUP_INDEX(J)),J=1,N_GROUPS)
        PRINT '(A,15E11.4)', ' TR       ',                              &
     &  (DEF_NUDGE_TR(GROUP_INDEX(J)),J=1,N_GROUPS)
        PRINT '(A,15E11.4)', ' SH       ',                              &
     &  (DEF_NUDGE_SH(GROUP_INDEX(J)),J=1,N_GROUPS)
#else
        PRINT '(A,15E11.4)', '          ',                              &
     &  (DEF_NUDGE_LAM(GROUP_INDEX(J)),J=1,N_GROUPS)
#endif

      ELSEIF (NACT == 0.AND.mype == 0) THEN

        PRINT *,' SETTPS : No observation types to process ?'
        ICODE = 1
        CMESSAGE = 'SETTPS : No obs types to process ?'
        GO TO 999

      ENDIF

 999  CONTINUE
      RETURN
      END SUBROUTINE SETTPS
#endif
