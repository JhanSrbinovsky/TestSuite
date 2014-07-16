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
      SUBROUTINE PEER_INITIALISE(                                       &
     &  ICODE,CMESSAGE)

! Description:
!   Creates output file to hold arrays and writes simple header.
!
! Current Code Owner: S.D.Mullerworth
!
! History:
!  Model    Date     Modification history from model version 3.5
!  version
!    4.5+   9/9/99   New DECK created. S.D.Mullerworth
!    5.3    05/03/01 Updates for new version. Now dependent on logical
!                    S.D.Mullerworth
!
! Subroutine Arguments:
      IMPLICIT NONE

      INTEGER                                                           &
     &  ICODE                  ! INOUT: Error return code

      CHARACTER*80 CMESSAGE     ! INOUT: Error message

! COMDECKS and common blocks

#include "parvars.h"
#include "c_writd.h"

      INTEGER                                                           &
     &  PEER_ADDRESS            ! Address pointer of peer file

      COMMON /PEER_ADDRESS / PEER_ADDRESS

! Local variables:

      INTEGER                                                           &
     &  HEADER(4)               ! Header for output file

      REAL                                                              &
     &  IOSTAT                  ! Return code from I/O routines

      INTEGER                                                           &
     &  PUNIT                                                           &
                                ! Unit number for output data file
     &  ,SUNIT                                                          &
                                ! Unit number for output summary file
     &  ,IOREQ                                                          &
                                ! Stores number of words to write out
     &  ,LENIO                  ! Length actually written out

      CHARACTER*30                                                      &
     &  PEER_FILENAME                                                   &
                                ! Filename for data
     &  ,SUMM_FILENAME          ! Filename for summary


      IF (L_PEER) THEN
! Initialise pointer
        PEER_ADDRESS=1

! Get file names and unit numbers, and open files for writing
! DEPENDS ON: peer_file_details
        CALL PEER_FILE_DETAILS(PUNIT,PEER_FILENAME,SUNIT,SUMM_FILENAME)
! File that holds the data
        CALL OPEN_SINGLE(PUNIT,PEER_FILENAME,20,1,1,ICODE)
        IF (ICODE /= 0)THEN
          WRITE(6,*)                                                    &
     &      'PEER: Error opening file ',PEER_FILENAME,' on unit ',PUNIT
          WRITE(6,*)'Error code ',ICODE
          CMESSAGE='PEER: Error opening file on one of the PEs'
          GOTO 999
        ENDIF
! File that holds details on each field
        OPEN(UNIT=SUNIT,FILE=SUMM_FILENAME,FORM='FORMATTED',            &
     &    IOSTAT=ICODE)
        IF (ICODE /= 0)THEN
          WRITE(6,*)                                                    &
     &      'PEER: Error opening file ',SUMM_FILENAME,' on unit ',SUNIT
          WRITE(6,*)'Error code ',ICODE
          CMESSAGE='PEER: Error opening file on one of the PEs'
          GOTO 999
        ENDIF
! Set the write position to the start of the file.

        WRITE(SUNIT,20)PEER_VN,nproc_x,nproc_y,mype
 20     FORMAT(4I5)
      ENDIF ! IF (L_PEER)

 999  CONTINUE
      RETURN
      END SUBROUTINE PEER_INITIALISE


! Subroutine interface:

! Subroutine interface:
#endif
