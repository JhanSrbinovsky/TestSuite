#if defined(A18_2A)







! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE TYPE_DEP_VAR -------------------------------------------
!LL
!LL  Purpose : Process &ACP Namelist arrays which are
!LL            Observation Type Dependent.
!LL
!LL  For use on Cray
!LL
!LL S.Bell      <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.2  8/7/93      Add OBTHIN                         S Bell
!    4.2 25/11/96: T3E mods Stuart Bell
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LLEND------------------------------------------------------------------
      SUBROUTINE TYPE_DEP_VAR(TIMEB,TIMEA,TGETOBB,TGETOBA,RADINF,OBTHIN,&
     &                         CSCALE_START,CSCALE_OBTIME,CSCALE_END,   &
     &                         ICODE,CMESSAGE)

      IMPLICIT NONE

#include "acparm.h"

      INTEGER TIMEA   (NOBTYPMX)         !IN time window after ob time
      INTEGER TIMEB   (NOBTYPMX)         !IN time window before ob time
      INTEGER TGETOBA (NOBTYPMX)         !IN window for getting obs (a)
      INTEGER TGETOBB (NOBTYPMX)         !IN window for getting obs (b)
      INTEGER RADINF  (NOBTYPMX)         !IN horiz infl radius
      INTEGER OBTHIN  (NOBTYPMX)         !IN obs thinning factors
      INTEGER CSCALE_START  (NOBTYPMX)   !IN horiz cor scale at start
      INTEGER CSCALE_OBTIME (NOBTYPMX)   !IN   "    "    "   "  ob time
      INTEGER CSCALE_END    (NOBTYPMX)   !IN   "    "    "   "   end
      INTEGER ICODE                      !OUT error code and message
      CHARACTER*256 CMESSAGE

#include "comacp.h"
#if defined(MPP)
#include "parvars.h"
#else
      INTEGER  mype
      PARAMETER (mype = 0 ) ! always zero in non-MPP code
#endif

!     Local arrays/variables.

      INTEGER ISCFACT(NOBTYPMX)
      INTEGER JOBT,JOBT2,AC_TYPE,J
      INTEGER FIRST_TYPE,LAST_TYPE,N_TYPES
      LOGICAL LFOUND

  30  FORMAT (A,F5.2,A,I4)

!L    Process TIMEB
!     -------------
      DO JOBT=1,NOBTYPMX
        IF (TIMEB(JOBT) >  0) THEN
          AC_TYPE = TIMEB(JOBT)/1000
          TIMEB(JOBT) = MOD(TIMEB(JOBT),1000)
          IF (TIMEB(JOBT) <  1 .OR. TIMEB(JOBT) >  999) THEN
            ICODE=1
            CMESSAGE = 'TYPE_DEP_VAR : Invalid value given for TIMEB'
            GO TO 999
          ENDIF
          IF (AC_TYPE == 999) THEN
!           Use new value for ALL Obs types
            DO JOBT2=1,NOBTYPMX
              DEF_TIMEB(JOBT2)=TIMEB(JOBT)
            ENDDO
            if(mype == 0)                                               &
     &      PRINT *, 'TIMEB : Defaults changed to ',                    &
     &                  TIMEB(JOBT),' mins for all obs types.'
          ELSE
            LFOUND = .FALSE.
            DO JOBT2=1,NOBTYPMX
              IF (AC_TYPE == MASTER_AC_TYPES(JOBT2)) THEN
                LFOUND = .TRUE.
                DEF_TIMEB(JOBT2)=TIMEB(JOBT)
                if(mype == 0)                                           &
     &          PRINT *, 'TIMEB : Default changed to ',                 &
     &                TIMEB(JOBT),' mins for obs type ',AC_TYPE
              ENDIF
            ENDDO
            IF (.NOT.LFOUND) THEN
              ICODE=1
              CMESSAGE =                                                &
     &        'TYPE_DEP_VAR : Invalid AC Obs type used in TIMEB'
              GO TO 999
            ENDIF
          ENDIF
        ENDIF
      ENDDO

!L    Process TIMEA
!     -------------
      DO JOBT=1,NOBTYPMX
        IF (TIMEA(JOBT) >  0) THEN
          AC_TYPE = TIMEA(JOBT)/1000
          TIMEA(JOBT) = MOD(TIMEA(JOBT),1000)
          IF (TIMEA(JOBT) <  1 .OR. TIMEA(JOBT) >  999) THEN
            ICODE=1
            CMESSAGE = 'TYPE_DEP_VAR : Invalid value given for TIMEA'
            GO TO 999
          ENDIF
          IF (AC_TYPE == 999) THEN
!           Use new value for ALL Obs types
            DO JOBT2=1,NOBTYPMX
              DEF_TIMEA(JOBT2)=TIMEA(JOBT)
            ENDDO
            if(mype == 0)                                               &
     &      PRINT *, 'TIMEA : Defaults changed to ',                    &
     &                  TIMEA(JOBT),' mins for all obs types.'
          ELSE
            LFOUND = .FALSE.
            DO JOBT2=1,NOBTYPMX
              IF (AC_TYPE == MASTER_AC_TYPES(JOBT2)) THEN
                LFOUND = .TRUE.
                DEF_TIMEA(JOBT2)=TIMEA(JOBT)
                if(mype == 0)                                           &
     &          PRINT *, 'TIMEA : Default changed to ',                 &
     &                TIMEA(JOBT),' mins for obs type ',AC_TYPE
              ENDIF
            ENDDO
            IF (.NOT.LFOUND) THEN
              ICODE=1
              CMESSAGE =                                                &
     &        'TYPE_DEP_VAR : Invalid AC Obs type used in TIMEA'
              GO TO 999
            ENDIF
          ENDIF
        ENDIF
      ENDDO

!L    Process TGETOBB
!     ---------------
      DO JOBT=1,NOBTYPMX
        IF (TGETOBB(JOBT) >  0) THEN
          AC_TYPE = TGETOBB(JOBT)/1000
          TGETOBB(JOBT) = MOD(TGETOBB(JOBT),1000)
          IF (TGETOBB(JOBT) <  1 .OR. TGETOBB(JOBT) >  999) THEN
            ICODE=1
            CMESSAGE = 'TYPE_DEP_VAR : Invalid value given in TGETOBB'
            GO TO 999
          ENDIF
          IF (AC_TYPE == 999) THEN
!           Use new value for ALL Obs types
            DO JOBT2=1,NOBTYPMX
              DEF_TGETOBB(JOBT2)=TGETOBB(JOBT)
            ENDDO
            if(mype == 0)                                               &
     &      PRINT *, 'TGETOBB : Default changed to ',                   &
     &                   TGETOBB(JOBT),' mins for all obs types'
          ELSE
            LFOUND = .FALSE.
            DO JOBT2=1,NOBTYPMX
              IF (AC_TYPE == MASTER_AC_TYPES(JOBT2)) THEN
                LFOUND = .TRUE.
                DEF_TGETOBB(JOBT2)=TGETOBB(JOBT)
                if(mype == 0)                                           &
     &          PRINT *, 'TGETOBB : Default changed to ',               &
     &                TGETOBB(JOBT),' mins for obs type ',AC_TYPE
              ENDIF
            ENDDO
            IF (.NOT.LFOUND) THEN
              ICODE=1
              CMESSAGE =                                                &
     &        'TYPE_DEP_VAR : Invalid AC Obs type used in TGETOBB'
              GO TO 999
            ENDIF
          ENDIF
        ENDIF
      ENDDO

!L    Process TGETOBA
!     ---------------
      DO JOBT=1,NOBTYPMX
        IF (TGETOBA(JOBT) >  0) THEN
          AC_TYPE = TGETOBA(JOBT)/1000
          TGETOBA(JOBT) = MOD(TGETOBA(JOBT),1000)
          IF (TGETOBA(JOBT) <  1 .OR. TGETOBA(JOBT) >  999) THEN
            ICODE=1
            CMESSAGE = 'TYPE_DEP_VAR : Invalid value given in TGETOBA'
            GO TO 999
          ENDIF
          IF (AC_TYPE == 999) THEN
!           Use new value for ALL Obs types
            DO JOBT2=1,NOBTYPMX
              DEF_TGETOBA(JOBT2)=TGETOBA(JOBT)
            ENDDO
            if(mype == 0)                                               &
     &      PRINT *, 'TGETOBA : Default changed to ',                   &
     &                   TGETOBA(JOBT),' mins for all obs types'
          ELSE
            LFOUND = .FALSE.
            DO JOBT2=1,NOBTYPMX
              IF (AC_TYPE == MASTER_AC_TYPES(JOBT2)) THEN
                LFOUND = .TRUE.
                DEF_TGETOBA(JOBT2)=TGETOBA(JOBT)
                if(mype == 0)                                           &
     &          PRINT *, 'TGETOBA : Default changed to ',               &
     &                       TGETOBA(JOBT),' mins for obs type',AC_TYPE
              ENDIF
            ENDDO
            IF (.NOT.LFOUND) THEN
              ICODE=1
              CMESSAGE =                                                &
     &        'TYPE_DEP_VAR : Invalid AC Obs type used in TGETOBA'
              GO TO 999
            ENDIF
          ENDIF
        ENDIF
      ENDDO

!L    Process RADINF
!     --------------
      DO JOBT=1,NOBTYPMX
        IF (RADINF(JOBT) >  0) THEN
          AC_TYPE = RADINF(JOBT)/1000
          RADINF(JOBT) = MOD(RADINF(JOBT),1000)
          IF (RADINF(JOBT) <  1 .OR. RADINF(JOBT) >  999) THEN
            ICODE=1
            CMESSAGE = 'TYPE_DEP_VAR : Invalid value given in RADINF'
            GO TO 999
          ENDIF
          IF (AC_TYPE == 999) THEN
!           Use new value for ALL Obs types
            DO JOBT2=1,NOBTYPMX
              DEF_RADINF(JOBT2)=RADINF(JOBT)/100.0
            ENDDO
            if(mype == 0)                                               &
     &      PRINT 30, ' RADINF  : Defaults changed to ',                &
     &                  DEF_RADINF(1),' for all obs types'
          ELSE
            LFOUND = .FALSE.
            DO JOBT2=1,NOBTYPMX
              IF (AC_TYPE == MASTER_AC_TYPES(JOBT2)) THEN
                LFOUND = .TRUE.
                DEF_RADINF(JOBT2)=RADINF(JOBT)/100.0
                if(mype == 0)                                           &
     &          PRINT 30, ' RADINF  : Default changed to ',             &
     &                      DEF_RADINF(JOBT2),' for obs type ',AC_TYPE
              ENDIF
            ENDDO
            IF (.NOT.LFOUND) THEN
              ICODE=1
              CMESSAGE =                                                &
     &        'TYPE_DEP_VAR : Invalid AC Obs type used in RADINF'
              GO TO 999
            ENDIF
          ENDIF
        ENDIF
      ENDDO


!L    Process OBTHIN
!     --------------
      DO JOBT=1,NOBTYPMX
        IF (OBTHIN(JOBT) >  0) THEN
          AC_TYPE = OBTHIN(JOBT)/1000
          OBTHIN(JOBT) = MOD(OBTHIN(JOBT),1000)
          IF (OBTHIN(JOBT) <  1 .OR. OBTHIN(JOBT) >  999) THEN
            ICODE=1
            CMESSAGE = 'TYPE_DEP_VAR : Invalid value given in OBTHIN'
            GO TO 999
          ENDIF
          IF (AC_TYPE == 999) THEN
!           Use new value for ALL Obs types
            DO JOBT2=1,NOBTYPMX
              DEF_OBTHIN(JOBT2)=OBTHIN(JOBT)
            ENDDO
            if(mype == 0)                                               &
     &      PRINT *, ' OBTHIN  : Defaults changed to ',                 &
     &                  DEF_OBTHIN(1),' for all obs types'
          ELSE
            LFOUND = .FALSE.
            DO JOBT2=1,NOBTYPMX
              IF (AC_TYPE == MASTER_AC_TYPES(JOBT2)) THEN
                LFOUND = .TRUE.
                DEF_OBTHIN(JOBT2)=OBTHIN(JOBT)
                if(mype == 0)                                           &
     &          PRINT *, ' OBTHIN  : Default changed to ',              &
     &                      DEF_OBTHIN(JOBT2),' for obs type ',AC_TYPE
              ENDIF
            ENDDO
            IF (.NOT.LFOUND) THEN
              ICODE=1
              CMESSAGE =                                                &
     &        'TYPE_DEP_VAR : Invalid AC Obs type used in OBTHIN'
              GO TO 999
            ENDIF
          ENDIF
        ENDIF
      ENDDO

!L    Process CSCALE_START
!     --------------------
      DO JOBT=1,NOBTYPMX
        IF (CSCALE_START(JOBT) >  0) THEN
          AC_TYPE = CSCALE_START(JOBT)/1000
          CSCALE_START(JOBT) = MOD(CSCALE_START(JOBT),1000)
          IF (CSCALE_START(JOBT) <  1 .OR.                              &
     &        CSCALE_START(JOBT) >  999) THEN
            ICODE=1
            CMESSAGE =                                                  &
     &      'TYPE_DEP_VAR : Invalid value given in CSCALE_START'
            GO TO 999
          ENDIF
          IF (AC_TYPE == 999) THEN
!           Use new value for ALL Obs types
            DO JOBT2=1,NOBTYPMX
              DEF_CSCALE_START(JOBT2)=CSCALE_START(JOBT)
            ENDDO
            if(mype == 0)                                               &
     &               PRINT *, 'CSCALE_START  : Defaults changed to ',   &
     &                   CSCALE_START(JOBT),' km for all obs types.'
          ELSE
            LFOUND = .FALSE.
            DO JOBT2=1,NOBTYPMX
              IF (AC_TYPE == MASTER_AC_TYPES(JOBT2)) THEN
                LFOUND = .TRUE.
                DEF_CSCALE_START(JOBT2)=CSCALE_START(JOBT)
                if(mype == 0)                                           &
     &          PRINT *, 'CSCALE_START  : Default changed to ',         &
     &          CSCALE_START(JOBT),' km for obs type ',AC_TYPE
              ENDIF
            ENDDO
            IF (.NOT.LFOUND) THEN
              ICODE=1
              CMESSAGE =                                                &
     &        'TYPE_DEP_VAR : Invalid AC Obs type used in CSCALE_START'
              GO TO 999
            ENDIF
          ENDIF
        ENDIF
      ENDDO

!L    Process CSCALE_OBTIME
!     ---------------------
      DO JOBT=1,NOBTYPMX
        IF (CSCALE_OBTIME(JOBT) >  0) THEN
          AC_TYPE = CSCALE_OBTIME(JOBT)/1000
          CSCALE_OBTIME(JOBT) = MOD(CSCALE_OBTIME(JOBT),1000)
          IF (CSCALE_OBTIME(JOBT) <  1 .OR.                             &
     &        CSCALE_OBTIME(JOBT) >  999) THEN
            ICODE=1
            CMESSAGE =                                                  &
     &      'TYPE_DEP_VAR : Invalid value given in CSCALE_OBTIME'
            GO TO 999
          ENDIF
          IF (AC_TYPE == 999) THEN
!           Use new value for ALL Obs types
            DO JOBT2=1,NOBTYPMX
              DEF_CSCALE_OBTIME(JOBT2)=CSCALE_OBTIME(JOBT)
            ENDDO
            if(mype == 0)                                               &
     &      PRINT *, 'CSCALE_OBTIME : Defaults changed to ',            &
     &                  CSCALE_OBTIME(JOBT),' km for all obs types.'
          ELSE
            LFOUND = .FALSE.
            DO JOBT2=1,NOBTYPMX
              IF (AC_TYPE == MASTER_AC_TYPES(JOBT2)) THEN
                LFOUND = .TRUE.
                DEF_CSCALE_OBTIME(JOBT2)=CSCALE_OBTIME(JOBT)
                if(mype == 0)                                           &
     &          PRINT *, 'CSCALE_OBTIME : Default changed to ',         &
     &          CSCALE_OBTIME(JOBT),' km for obs type ',AC_TYPE
              ENDIF
            ENDDO
            IF (.NOT.LFOUND) THEN
              ICODE=1
              CMESSAGE =                                                &
     &        'TYPE_DEP_VAR : Invalid Obs type used in CSCALE_OBTIME'
              GO TO 999
            ENDIF
          ENDIF
        ENDIF
      ENDDO

!L    Process CSCALE_END
!     ------------------
      DO JOBT=1,NOBTYPMX
        IF (CSCALE_END(JOBT) >  0) THEN
          AC_TYPE = CSCALE_END(JOBT)/1000
          CSCALE_END(JOBT) = MOD(CSCALE_END(JOBT),1000)
          IF (CSCALE_END(JOBT) <  1 .OR.                                &
     &        CSCALE_END(JOBT) >  999) THEN
            ICODE=1
            CMESSAGE =                                                  &
     &      'TYPE_DEP_VAR : Invalid value given in CSCALE_END'
            GO TO 999
          ENDIF
          IF (AC_TYPE == 999) THEN
!           Use new value for ALL Obs types
            DO JOBT2=1,NOBTYPMX
              DEF_CSCALE_END(JOBT2)=CSCALE_END(JOBT)
            ENDDO
            if(mype == 0)                                               &
     &      PRINT *, 'CSCALE_END    : Defaults changed to ',            &
     &                   CSCALE_END(JOBT),' km for all obs types.'
          ELSE
            LFOUND = .FALSE.
            DO JOBT2=1,NOBTYPMX
              IF (AC_TYPE == MASTER_AC_TYPES(JOBT2)) THEN
                LFOUND = .TRUE.
                DEF_CSCALE_END(JOBT2)=CSCALE_END(JOBT)
                if(mype == 0)                                           &
     &          PRINT *, 'CSCALE_END    : Default changed to ',         &
     &          CSCALE_END(JOBT),' km for obs type ',AC_TYPE
              ENDIF
            ENDDO
            IF (.NOT.LFOUND) THEN
              ICODE=1
              CMESSAGE =                                                &
     &        'TYPE_DEP_VAR : Invalid AC Obs type used in CSCALE_END'
              GO TO 999
            ENDIF
          ENDIF
        ENDIF
      ENDDO

!L    Process NO_SCFACT
!     -----------------
      DO JOBT=1,NOBTYPMX
        ISCFACT(JOBT)   = NO_SCFACT(JOBT)
        NO_SCFACT(JOBT) = 1   !  Apply scale factor to all types.
      ENDDO
      DO JOBT=1,NOBTYPMX
        IF (ISCFACT(JOBT) >  0) THEN
          LFOUND = .FALSE.
          DO JOBT2=1,NOBTYPMX
            IF (ISCFACT(JOBT) == MASTER_AC_TYPES(JOBT2)) THEN
                LFOUND = .TRUE.
                NO_SCFACT(JOBT2) = 0   !  No scale factor for this type
            ENDIF
          ENDDO
          IF (.NOT.LFOUND) THEN
            ICODE=1
            CMESSAGE =                                                  &
     &      'TYPE_DEP_VAR : Invalid Obs type used in NO_SCFACT'
            GO TO 999
          ENDIF
        ENDIF
      ENDDO

      N_TYPES = 0
      DO JOBT=1,NOBTYPMX
        IF (MASTER_AC_TYPES(JOBT) >  0) THEN
          N_TYPES = N_TYPES+1
        ENDIF
      ENDDO

      DO J=1,(N_TYPES+6)/7
        FIRST_TYPE=(J-1)*7+1
        LAST_TYPE =MIN(J*7,N_TYPES)
          if(mype == 0)then
        PRINT *, ' '
        PRINT '(A,7I8)',   ' AC Obs Types      ',                       &
     &   (MASTER_AC_TYPES(JOBT),JOBT=FIRST_TYPE,LAST_TYPE)
        PRINT '(A,7F8.1)', ' TIMEB (mins)      ',                       &
     &   (DEF_TIMEB(JOBT),JOBT=FIRST_TYPE,LAST_TYPE)
        PRINT '(A,7F8.1)', ' TIMEA (mins)      ',                       &
     &   (DEF_TIMEA(JOBT),JOBT=FIRST_TYPE,LAST_TYPE)
        PRINT '(A,7F8.1)', ' TGETOBB (mins)    ',                       &
     &   (DEF_TGETOBB(JOBT),JOBT=FIRST_TYPE,LAST_TYPE)
        PRINT '(A,7F8.1)', ' TGETOBA (mins)    ',                       &
     &   (DEF_TGETOBA(JOBT),JOBT=FIRST_TYPE,LAST_TYPE)
        PRINT '(A,7F8.2)', ' RADINF            ',                       &
     &   (DEF_RADINF(JOBT),JOBT=FIRST_TYPE,LAST_TYPE)
        PRINT '(A,7I8)', ' OBTHIN            ',                         &
     &   (DEF_OBTHIN(JOBT),JOBT=FIRST_TYPE,LAST_TYPE)
        PRINT '(A,7F8.1)', ' CSCALE_START (km) ',                       &
     &   (DEF_CSCALE_START(JOBT),JOBT=FIRST_TYPE,LAST_TYPE)
        PRINT '(A,7F8.1)', ' CSCALE_OBTIME (km)',                       &
     &   (DEF_CSCALE_OBTIME(JOBT),JOBT=FIRST_TYPE,LAST_TYPE)
        PRINT '(A,7F8.1)', ' CSCALE_END (km)   ',                       &
     &   (DEF_CSCALE_END(JOBT),JOBT=FIRST_TYPE,LAST_TYPE)
        PRINT '(A,7I8)',   ' NO_SCFACT         ',                       &
     &   (NO_SCFACT(JOBT),JOBT=FIRST_TYPE,LAST_TYPE)
          endif
      ENDDO

!L    Convert Horizontal Correlation Scales from Km to metres.
      DO JOBT=1,NOBTYPMX
        DEF_CSCALE_START(JOBT)  = DEF_CSCALE_START(JOBT)  * 1000.0
        DEF_CSCALE_OBTIME(JOBT) = DEF_CSCALE_OBTIME(JOBT) * 1000.0
        DEF_CSCALE_END(JOBT)    = DEF_CSCALE_END(JOBT)    * 1000.0
      ENDDO

 999  CONTINUE
      RETURN
      END SUBROUTINE TYPE_DEP_VAR
#endif
