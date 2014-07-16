#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine Interface:
      Subroutine Oper_Emergency

      IMPLICIT NONE
!
!    Method:  Reset UMUI determined namelist items to alternate
!             predetermined values dependent on value of Environment
!             variables known to the operational suite
!
! Current Code Owner:  Tim Westmacott
!
! History:
! version  date      comment
! -------  ----      -------
!  4.5    3/09/98   New.  Stuart Bell
!  5.0   21/06/99   Amendments for 'C-P C dynamics' grid upgrade; wind
!                   limit test no longer supported. R Rawlins
!  5.1   13/04/00   New (temporary) RUN_ASSIM_MODE for FASTRUNs.
!                   Adam Clayton.
!  5.1   22/02/00   Move PARVARS for TYPSIZE                 P.Burton
!  5.3   29/11/00   Define run length changes wrt analysis time for
!  5.3              global runs. Adam Clayton
!  5.4   12/08/02   ONLYTO48 added. ONLYTO6 renamed ONLYTO9. A Clayton
!  5.4   22/03/02   Remove comment on same line as #include
!                                                S. Carroll
!  6.2   05/01/06   Move ONLYTO9 outside #if defined(GLOBAL)
!                   pre-processor directive M.Bush
!  6.2   05/01/06   Enable ONLYTO54 switch for Mesoscale ALTERNATE runs
!                   J.Bornemann

!
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered: Control
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):
#include "cmaxsize.h"
#include "csubmodl.h"
#include "chsunits.h"
#include "ccontrol.h"
#include "parvars.h"
! Print status information
#include "cprintst.h"

! Local scalars:
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='OperatorEmergency')
      INTEGER ICODE            ! return code from FORT_GET_ENV
      CHARACTER*80  ONLYTO3    ! value of EV MES SHORT RUN
      CHARACTER*80  ONLYTO9    ! value of EV SHORT RUN
      CHARACTER*80  ONLYTO12   ! value of EV MOGREPS SHORT RUN
      CHARACTER*80  ONLYTO15   ! value of EV MOGREPS SHORT RUN
      CHARACTER*80  ONLYTO18   ! value of EV MOGREPS SHORT RUN
      CHARACTER*80  ONLYTO21   ! value of EV MOGREPS SHORT RUN
      CHARACTER*80  ONLYTO24   ! value of EV FOAM SHORT RUN
      CHARACTER*80  ONLYTO42   ! value of EV EXTENDED UK4
      CHARACTER*80  ONLYTO48   ! value of EV GL 2DAY RUN
      CHARACTER*80  ONLYTO54   ! value of EV MES ALTERNATE RUN
      CHARACTER*80  ONLYTO72   ! value of EV GL 3DAY RUN
      CHARACTER*80  FASTRUN    ! value of EV FASTRUN
      CHARACTER*80  SHORTSTEP  ! value of EV SHORTSTEP

! Function & Subroutine calls:
      External FORT_GET_ENV

! FASTRUN
       CALL FORT_GET_ENV('FASTRUN',7,FASTRUN,80,ICODE)
       if(FASTRUN == 'true'.AND.ICODE == 0)then
        RUN_ASSIM_MODE='FastRun   '
       IF(PrintStatus >= PrStatus_Oper) THEN
        if(mype == 0)WRITE(6,*)                                         &
     &    "FASTRUN=true, setting RUN_ASSIM_MODE=",RUN_ASSIM_MODE
       ENDIF  ! PrintStatus test
       endif

#if defined(ATMOS)
! SHORTSTEP
       CALL FORT_GET_ENV('SHORTSTEP',9,SHORTSTEP,80,ICODE)
       if(SHORTSTEP == 'true'.AND.ICODE == 0)then
        IF(PrintStatus >= PrStatus_Oper) THEN
         WRITE(6,*) RoutineName,':Warning, Wind limit re-setting '      &
     &             ,'for SHORTSTEP no longer supported'
        ENDIF  ! PrintStatus test
       endif

#if defined(GLOBAL)
! UM6.5 : MODEL_ANALYSIS_HRS replaced by MODEL_ANALYSIS_MINS
! ONLYTO72 for GL ATMOS plus 6 hour assm
       CALL FORT_GET_ENV('ONLYTO72',8,ONLYTO72,80,ICODE)
       if(ONLYTO72 == 'true'.AND.ICODE == 0)then
       RUN_TARGET_END(3)=3                    ! days
       RUN_TARGET_END(4)=REAL(MODEL_ANALYSIS_MINS)/60.0 +3 ! hours
       RUN_TARGET_END(5)=REAL(MODEL_ANALYSIS_MINS)+(60.0*72)-           &
     &                   60.0*RUN_TARGET_END(4)          ! minutes 
       IF(PrintStatus >= PrStatus_Oper) THEN
        if(mype == 0)WRITE(6,*)                                         &
     &     "ONLYTO72=true, setting RUN_TARGET_END=",RUN_TARGET_END
       ENDIF  ! PrintStatus test
       endif

! ONLYTO48 for GL ATMOS
       CALL FORT_GET_ENV('ONLYTO48',8,ONLYTO48,80,ICODE)
       if(ONLYTO48 == 'true'.AND.ICODE == 0)then
       RUN_TARGET_END(3)=2                  ! days
       RUN_TARGET_END(4)=REAL(MODEL_ANALYSIS_MINS)/60.0 ! hours
       RUN_TARGET_END(5)=REAL(MODEL_ANALYSIS_MINS)+(60.0*48)-           &
     &                   60.0*RUN_TARGET_END(4)          ! minutes 
       IF(PrintStatus >= PrStatus_Oper) THEN
        if(mype == 0)WRITE(6,*)                                         &
     &     "ONLYTO48=true, setting RUN_TARGET_END=",RUN_TARGET_END
       ENDIF  ! PrintStatus test
       endif

! ONLYTO21 for MOGREPS GL ATMOS
       CALL FORT_GET_ENV('ONLYTO21',8,ONLYTO21,80,ICODE)
       if(ONLYTO21 == 'true'.AND.ICODE == 0)then
       RUN_TARGET_END(3)=0                  ! days
       RUN_TARGET_END(4)=REAL(MODEL_ANALYSIS_MINS)/60.0+21 ! hours
       RUN_TARGET_END(5)=REAL(MODEL_ANALYSIS_MINS)+(60.0*21)-           &
     &                   60.0*RUN_TARGET_END(4)          ! minutes 
       IF(PrintStatus >= PrStatus_Oper) THEN
        if(mype == 0)WRITE(6,*)                                         &
     &     "ONLYTO21=true, setting RUN_TARGET_END=",RUN_TARGET_END
       ENDIF  ! PrintStatus test
       endif

! ONLYTO18 for MOGREPS GL ATMOS
       CALL FORT_GET_ENV('ONLYTO18',8,ONLYTO18,80,ICODE)
       if(ONLYTO18 == 'true'.AND.ICODE == 0)then
       RUN_TARGET_END(3)=0                  ! days
       RUN_TARGET_END(4)=REAL(MODEL_ANALYSIS_MINS)/60.0+18 ! hours
       RUN_TARGET_END(5)=REAL(MODEL_ANALYSIS_MINS)+(60.0*18)-           &
     &                   60.0*RUN_TARGET_END(4)          ! minutes 
       IF(PrintStatus >= PrStatus_Oper) THEN
        if(mype == 0)WRITE(6,*)                                         &
     &     "ONLYTO18=true, setting RUN_TARGET_END=",RUN_TARGET_END
       ENDIF  ! PrintStatus test
       endif

#else
! ONLYTO3 for MES ATMOS plus 3 hours assm
       CALL FORT_GET_ENV('ONLYTO3',7,ONLYTO3,80,ICODE)
       if(ONLYTO3 == 'true'.AND.ICODE == 0)then
       RUN_TARGET_END(3)=0  !days
       RUN_TARGET_END(4)=REAL(MODEL_ANALYSIS_MINS)/60.0+4 ! hours
       RUN_TARGET_END(5)=REAL(MODEL_ANALYSIS_MINS)+(60.0*4)-            &
     &                   60.0*RUN_TARGET_END(4)          ! minutes 
       IF(PrintStatus >= PrStatus_Oper) THEN
        if(mype == 0)WRITE(6,*)                                         &
     &     "ONLYTO3=true, setting RUN_TARGET_END=",RUN_TARGET_END
       ENDIF  ! PrintStatus test
       endif

! ONLYTO42 for EXTENDED UK4 ATMOS
       CALL FORT_GET_ENV('ONLYTO42',8,ONLYTO42,80,ICODE)
       if (ONLYTO42 == 'true' .AND. ICODE == 0) then
       RUN_TARGET_END(3)=1  !days
       RUN_TARGET_END(4)=REAL(MODEL_ANALYSIS_MINS)/60.0+18 ! hours
       RUN_TARGET_END(5)=REAL(MODEL_ANALYSIS_MINS)+(60.0*42)-           &
     &                   60.0*RUN_TARGET_END(4)          ! minutes 
       IF(PrintStatus >= PrStatus_Oper) THEN
        if(mype == 0)WRITE(6,*)                                         &
           "ONLYTO42=true, setting RUN_TARGET_END=",RUN_TARGET_END
       ENDIF  ! PrintStatus test
       endif

! ONLYTO54 for ALTERNATE MES ATMOS
       CALL FORT_GET_ENV('ONLYTO54',8,ONLYTO54,80,ICODE)
       if(ONLYTO54 == 'true'.AND.ICODE == 0)then
       RUN_TARGET_END(3)=2  !days
       RUN_TARGET_END(4)=REAL(MODEL_ANALYSIS_MINS)/60.0+6 ! hours
       RUN_TARGET_END(5)=REAL(MODEL_ANALYSIS_MINS)+(60.0*54)-           &
     &                   60.0*RUN_TARGET_END(4)          ! minutes 
       IF(PrintStatus >= PrStatus_Oper) THEN
        if(mype == 0)WRITE(6,*)                                         &
     &     "ONLYTO54=true, setting RUN_TARGET_END=",RUN_TARGET_END
       ENDIF  ! PrintStatus test
       endif
#endif
!  ONLYTO9 for ATMOS
      CALL FORT_GET_ENV('ONLYTO9',7,ONLYTO9,80,ICODE)
      if(ONLYTO9 == 'true'.AND.ICODE == 0)then
      RUN_TARGET_END(3)=0                    ! days
      RUN_TARGET_END(4)=REAL(MODEL_ANALYSIS_MINS)/60.0+9 ! hours
      RUN_TARGET_END(5)=REAL(MODEL_ANALYSIS_MINS)+(60.0*9)-             &
     &                  60.0*RUN_TARGET_END(4)           ! minutes 
      IF(PrintStatus >= PrStatus_Oper) THEN
      if(mype == 0)WRITE(6,*)                                           &
     &   "ONLYTO9=true, setting RUN_TARGET_END=",RUN_TARGET_END
      ENDIF  ! PrintStatus test
      endif
!  ONLYTO12 for ATMOS
      CALL FORT_GET_ENV('ONLYTO12',8,ONLYTO12,80,ICODE)
      if(ONLYTO12 == 'true'.AND.ICODE == 0)then
      RUN_TARGET_END(3)=0                    ! days
      RUN_TARGET_END(4)=REAL(MODEL_ANALYSIS_MINS)/60.0+12 ! hours
      RUN_TARGET_END(5)=REAL(MODEL_ANALYSIS_MINS)+(60.0*12)-            &
     &                  60.0*RUN_TARGET_END(4)            ! minutes 
      IF(PrintStatus >= PrStatus_Oper) THEN
      if(mype == 0)WRITE(6,*)                                           &
     &   "ONLYTO12=true, setting RUN_TARGET_END=",RUN_TARGET_END
      ENDIF  ! PrintStatus test
      endif
!  ONLYTO15 for ATMOS
      CALL FORT_GET_ENV('ONLYTO15',8,ONLYTO15,80,ICODE)
      if(ONLYTO15 == 'true'.AND.ICODE == 0)then
      RUN_TARGET_END(3)=0                    ! days
      RUN_TARGET_END(4)=REAL(MODEL_ANALYSIS_MINS)/60.0+15 ! hours
      RUN_TARGET_END(5)=REAL(MODEL_ANALYSIS_MINS)+(60.0*15)-            &
     &                  60.0*RUN_TARGET_END(4)            ! minutes 
      IF(PrintStatus >= PrStatus_Oper) THEN
      if(mype == 0)WRITE(6,*)                                           &
     &   "ONLYTO15=true, setting RUN_TARGET_END=",RUN_TARGET_END
      ENDIF  ! PrintStatus test
      endif
#endif

      RETURN
      END SUBROUTINE Oper_Emergency
#endif
