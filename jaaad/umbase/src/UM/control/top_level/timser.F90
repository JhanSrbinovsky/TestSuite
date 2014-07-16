#if defined(CONTROL) || defined(FLDOP)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Process time-series domain data (if any)
! Subroutine interface:
      SUBROUTINE TIMSER(CMESSAGE,ErrorStatus)
      IMPLICIT NONE

! Description:
!
! Method:
!
! Current code owner:  S.J.Swarbrick
!
! History:
! Version   Date       Comment
! =======   ====       =======
!   4.0     Sept. 95   Original code.  S.J.Swarbrick
!   4.1     Apr.  96   Rationalise MDI   S.J.Swarbrick
!LL  5.1  22/02/00  Add PARPARM for TYPSIZE                 P.Burton
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
! Global variables:
#include "lenfil.h"
#include "csubmodl.h"
#include "version.h"
#include "parparm.h"
#include "typsize.h"
#include "model.h"
#include "cstash.h"
#include "stextend.h"
#include "c_mdi.h"

! Subroutine arguments
!   Scalar arguments with intent(out):
      CHARACTER*80 CMESSAGE      ! Error return message

! Error status:
      INTEGER        ErrorStatus ! Error return code

! Local variables:
      INTEGER BlkId           !Time series block identifier
      INTEGER BlkSt           !Start position of ts block data
      INTEGER Nrecs_prev      !No of recs in previous time ser block
      INTEGER IDP             !Domain profile loop counter
      INTEGER IPOS            !Position in ts limits arrays
      INTEGER ISBLIM,ISTLIM   !Used for converting vertical ts
      INTEGER IB,IT,IL,ILVL   !  domain limits to sequence nos.

!- End of Header ------------------------------------------------------

!Loop over domain profiles
      BlkSt =1
      DO IDP=1,NDPROF
        IF(NPOS_TS(IDP) >  0) THEN
!  This domain profile has a time series
!    Identify TS block using pointer array
          BlkId = NPOS_TS (IDP)
!    Find start position (in LIM_TS arrays) of data for this block
          IF (BlkId >  1) THEN
            BlkSt=BlkSt+Nrecs_prev
          END IF
!  Loop over records in ts block corresponding to domain profile IDP.
!  Adjust the TS records for domain profiles with vertical or horiz
!    averaging.
!  Convert the ts domain vertical limits to sequence nos.
!    in the domain profile levels range/levels list.
          DO IPOS=BlkSt,BlkSt+NRECS_TS(NPOS_TS(IDP))-1
!    Vertical levels
            IF(IOPL_D(IDP) == 1 .OR. IOPL_D(IDP) == 2 .OR.              &
     &                               IOPL_D(IDP) == 6) THEN
!           Model levels
              IF(IMN_D(IDP) == 1) THEN
!             Vertical mean
                BLIM_TS(IPOS)=1
                TLIM_TS(IPOS)=1
              ELSE
!             No vertical mean
                IF(LEVB_D(IDP) >= 0) THEN
!               Range of model levels
                  IF(BLIM_TS(IPOS) <  LEVB_D(IDP) .OR.                  &
     &              TLIM_TS(IPOS) >  LEVT_D(IDP)) THEN
                    WRITE(6,*) 'ERROR, TIMSER: ',                       &
     &            ' TS_DOMAIN LEVEL LIMIT OUT OF RANGE; ',              &
     &            ' DOM PROF: ',IDP,                                    &
     &            ' TS RECORD: ',IPOS
                    ErrorStatus=1
                    CMESSAGE='TS DOMAIN LEVEL LIMIT OUT OF RANGE'
                    GO TO 999
                  END IF
                  BLIM_TS(IPOS)=BLIM_TS(IPOS)-LEVB_D(IDP)+1
                  TLIM_TS(IPOS)=TLIM_TS(IPOS)-LEVB_D(IDP)+1
                ELSE
!               List of selected model levels;
!               LEVT_D(IDP)=no. of levels in list
                  ISBLIM=IMDI
                  ISTLIM=IMDI
                  DO IL=1,LEVT_D(IDP)
                    IF(BLIM_TS(IPOS) == LEVLST_D(IL,IDP)) ISBLIM=IL
                    IF(TLIM_TS(IPOS) == LEVLST_D(IL,IDP)) ISTLIM=IL
                  END DO
                  IF((ISTLIM == IMDI).OR.                               &
     &               (ISBLIM == IMDI)) THEN
                    WRITE(6,*)                                          &
     &             'ERROR TIMSER:T-SERIES INTEGER LEVEL NOT IN ',       &
     &             'LEVELS LIST; DOM PROF: ',IDP,' TS RECORD: ',IPOS
                    WRITE(6,*) 'SPECIFIED TS LEVELS LIMITS: ',          &
     &              BLIM_TS(IPOS),TLIM_TS(IPOS)
                    ErrorStatus = 1
                    CMESSAGE=                                           &
     &             'ERROR TIMSER:T-SERIES LEVEL NOT IN LEVELS LIST'
                    GO TO 999
                  END IF
!                 Store seq. nos. of ts domain level limits
                  BLIM_TS(IPOS)=ISBLIM
                  TLIM_TS(IPOS)=ISTLIM
                END IF
              END IF
!           List of specified real levels
            ELSE IF((IOPL_D(IDP) /= 5).AND.(IOPL_D(IDP) <= 9)) THEN
              IF(IMN_D(IDP) == 1) THEN
                BLIM_TS(IPOS)=1
                TLIM_TS(IPOS)=1
              ELSE
!             Determine sequence nos. of top & bottom ts domain
!             levels in real levels list (ISBLIM, ISTLIM), by
!             representing real level values as integers.
                ISBLIM=IMDI
                ISTLIM=IMDI
                IB=(BLIMR_TS(IPOS)*1000.+0.5)
                IT=(TLIMR_TS(IPOS)*1000.+0.5)
                DO IL=1,LEVT_D(IDP)
                  ILVL=(RLEVLST_D(IL,IDP)*1000.+0.5)
                  IF(IB == ILVL) ISBLIM=IL
                  IF(IT == ILVL) ISTLIM=IL
                END DO
                IF((ISTLIM == IMDI).OR.                                 &
     &             (ISBLIM == IMDI)) THEN
                  WRITE(6,*)                                            &
     &           'ERROR TIMSER:T-SERIES REAL LEVEL NOT IN ',            &
     &           'LEVELS LIST; DOM PROF: ',IDP,' TS RECORD: ',IPOS
                  WRITE(6,*) 'SPECIFIED TS LEVELS LIMITS: ',            &
     &            BLIMR_TS(IPOS),TLIMR_TS(IPOS)
                  ErrorStatus = 1
                  CMESSAGE=                                             &
     &           'ERROR TIMSER:T-SERIES LEVEL NOT IN LEVELS LIST'
                END IF
!               Store seq. nos. of ts domain level limits
                BLIM_TS(IPOS)=ISBLIM
                TLIM_TS(IPOS)=ISTLIM
              END IF
            ELSE IF(IOPL_D(IDP) == 5) THEN
!           Single level
              BLIM_TS(IPOS)=1
              TLIM_TS(IPOS)=1
            ELSE
              WRITE(6,*)                                                &
     &       'ERROR TIMSER: UNEXPECTED LEVEL TYPE CODE',IOPL_D(IDP)
              ErrorStatus=1
              GO TO 999
            END IF
!    Horizontal area
            IF(IMN_D(IDP) == 2) THEN
              ELIM_TS(IPOS)=1
              WLIM_TS(IPOS)=1
            ELSE IF(IMN_D(IDP) == 3) THEN
              NLIM_TS(IPOS)=1
              SLIM_TS(IPOS)=1
            ELSE IF(IMN_D(IDP) == 4) THEN
              ELIM_TS(IPOS)=1
              WLIM_TS(IPOS)=1
              NLIM_TS(IPOS)=1
              SLIM_TS(IPOS)=1
            END IF
            IG_TS =0  ! These constants are left-overs from the
            I1_TS =1  !  pre-vn3.5 TIMSER routine: they are used
            I51_TS=51 !  in the UM time-series routines.
          END DO      ! IPOS loop
        Nrecs_prev=NRECS_TS(NPOS_TS(IDP)) ! For next TS block
        END IF        ! TS(IDP) == 'Y'
      END DO          ! IDP loop

 999  RETURN
      END SUBROUTINE TIMSER
#endif
