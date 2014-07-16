
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

! Subroutine interface:
      SUBROUTINE OUTPUT_ARRAY(                                          &
     &  FIELD                                                           &
     &  ,FIELD_NAME,FIELD_NUMBER,XPT,YPT,NPT                            &
     &  ,START_POINT,Xhalo,YHalo                                        &
     &  ,ICODE,CMESSAGE)

! Description:
!  Output array with header to file.
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
     &  FIELD_NUMBER                                                    &
                                ! IN: ID of field
     &  ,XPT                                                            &
                                ! IN: Nominal X dim of array inc halo
     &  ,YPT                                                            &
                                ! IN: Nominal Y dim of array inc halo
     &  ,NPT                                                            &
                                ! IN: Actual number of point to output
     &  ,START_POINT                                                    &
                                ! IN: Start position of data relative
                                !     to full XPT by YPT grid
     &  ,XHalo                                                          &
                                ! IN: Size of X halo
     &  ,YHalo                                                          &
                                ! IN: Size of Y halo
     &  ,ICODE                  ! INOUT: Error return code

      REAL                                                              &
     &  FIELD(NPT)              ! IN: Field to write out

      CHARACTER*(*)                                                     &
     &  FIELD_NAME              ! IN: Name of field - for info

      CHARACTER*80 CMESSAGE     ! INOUT: Error message

! Common blocks
!
! This Comdeck declares and stores the logical and integer
! variables used in the time-step control for writing general
! data.
!
! Code owner: S.J.Swarbrick
!
! Switch which activates output of arrays for Peer utility
      LOGICAL L_PEER

! Switches which activate dump writing
      LOGICAL                                                           &
     &  L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                  &
     &  ,L_WRIT_INIT

! Timesteps for dump writing
      INTEGER                                                           &
     &  T_WRITD1_START                                                  &
                              ! First timestep
     &  ,T_WRITD1_END                                                   &
                              ! Last timestep
     &  ,T_WRITD1_INT         ! Timestep interval between dumps

      INTEGER                                                           &
     &  PEER_VN                  ! Version of PEER utility

      NAMELIST/NLSTWRITEDATA/                                           &
     &  L_PEER,PEER_VN                                                  &
     &  ,T_WRITD1_START  ,T_WRITD1_END  ,T_WRITD1_INT                   &
     &  ,L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                 &
     &  ,L_WRIT_INIT

      COMMON/WRITEDATA/                                                 &
     &  L_PEER,PEER_VN                                                  &
     &  ,T_WRITD1_START  ,T_WRITD1_END  ,T_WRITD1_INT                   &
     &  ,L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                 &
     &  ,L_WRIT_INIT
!
      INTEGER                                                           &
     &  PEER_ADDRESS            ! Address pointer of peer file

      COMMON /PEER_ADDRESS / PEER_ADDRESS

! Local variables
      INTEGER                                                           &
     &  HEADER(7)               ! Header for field

      INTEGER                                                           &
     &  PUNIT                                                           &
                                ! Unit number for data file
     &  ,SUNIT                                                          &
                                ! Unit number for summary file
     &  ,IOREQ                                                          &
                                ! Stores number of words to write out
     &  ,LENIO                  ! Number of words actually written

      REAL                                                              &
     &  IOSTAT                  ! Error code from BUFFOUT. -1 means OK

      CHARACTER*80                                                      &
     &  FORMATTED_NAME          ! Ensure all names are 80 characters

      CHARACTER*30                                                      &
     &  PEER_FILENAME                                                   &
                                ! Filename for data
     &  ,SUMM_FILENAME          ! Filename for summary

      IF (L_PEER) THEN
! Get filename and unit number, and open file for writing
! DEPENDS ON: peer_file_details
        CALL PEER_FILE_DETAILS(PUNIT,PEER_FILENAME,SUNIT,SUMM_FILENAME)

! Set the write position to the appropriate part of the file.
!      CALL SETPOS_SINGLE(PUNIT,PEER_ADDRESS,ICODE)
!      IF (ICODE /= 0)THEN
!        WRITE(6,*)
!     &    'PEER: Error in SET_POS: file ',PEER_FILENAME,' unit ',PUNIT
!        WRITE(6,*)'Error code ',ICODE
!        CMESSAGE='PEER: Error accessing file on one of the pes'
!        GOTO 999
!      ENDIF

        HEADER(1)=PEER_ADDRESS
        HEADER(2)=NPT
        HEADER(3)=XPT
        HEADER(4)=YPT
        HEADER(5)=START_POINT
        HEADER(6)=XHalo
        HEADER(7)=YHalo

        FORMATTED_NAME=FIELD_NAME
        WRITE(SUNIT,*)FIELD_NUMBER,' ',FORMATTED_NAME
        WRITE(SUNIT,10)HEADER
 10     FORMAT(2I14,2I7,I14,2I7)

        CALL BUFFOUT_SINGLE(PUNIT,FIELD,NPT,LENIO,IOSTAT)
        IF (LENIO /= NPT.OR.IOSTAT /= -1.0)THEN
! DEPENDS ON: ioerror
          CALL IOERROR(                                                 &
     &      'OUTPUT_ARRAY: Failed to write field to PEER file'          &
     &      ,IOSTAT,NPT,LENIO)
          GOTO 999
        ENDIF
        PEER_ADDRESS=PEER_ADDRESS+NPT
      ENDIF ! IF (L_PEER)

 999  CONTINUE
      RETURN
      END SUBROUTINE OUTPUT_ARRAY
