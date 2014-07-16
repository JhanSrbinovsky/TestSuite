#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routines to initialise and write to PEER output files which can then
! be analysed with the PEER utility.
!
!


! Subroutine interface:
      SUBROUTINE PEER_FILE_DETAILS(                                     &
     &  PUNIT,PEER_FILENAME,SUNIT,SUMM_FILENAME                         &
     &  )

! Description:
!   Returns unit number and filename for sending output
!
! Current Code Owner: S.D.Mullerworth
!
! History:
!  Model    Date     Modification history from model version 3.5
!  version
!    4.5+   9/9/99   New DECK created. S.D.Mullerworth
!
! Subroutine Arguments:

      IMPLICIT NONE

      INTEGER                                                           &
     &  PUNIT                                                           &
                                ! OUT: Unit number for peer file
     &  ,SUNIT                  ! OUT: Unit number for summary file

      CHARACTER*30                                                      &
     &  PEER_FILENAME                                                   &
                                ! OUT: Peer filename
     &  ,SUMM_FILENAME          ! OUT: Summary filename

#include "chsunits.h"
#include "cntlall.h"
#include "parvars.h"

! Create filename for output from this pe
      WRITE(PEER_FILENAME,10)mype
      PEER_FILENAME(1:4)=EXPT_ID
      PEER_FILENAME(5:5)=JOB_ID
      WRITE(SUMM_FILENAME,20)mype
      SUMM_FILENAME(1:4)=EXPT_ID
      SUMM_FILENAME(5:5)=JOB_ID
 10   FORMAT('......peer.pe',I3.3)
 20   FORMAT('......summ.pe',I3.3)

      PUNIT=58
      SUNIT=59

      RETURN
      END SUBROUTINE PEER_FILE_DETAILS

! Subroutine interface:
#endif
