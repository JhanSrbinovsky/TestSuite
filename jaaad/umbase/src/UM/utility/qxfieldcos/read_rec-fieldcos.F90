#if defined(FLDC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!LL  Subroutine: READ_REC
!LL
!LL  Purpose: To read a data record from a  pp file
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: ...
!LL
!LL  Project task: ...
!LL
!LL  External documentation:
!LL
!
      SUBROUTINE READ_REC(FIELD,NUM_CRAY_WORDS,IWA,PPUNIT,              &
     &                    ICODE,CMESSAGE)
      IMPLICIT NONE
!     arguments
      CHARACTER CMESSAGE*(*)      !OUT error message
      INTEGER                                                           &
     &     NUM_CRAY_WORDS                                               &
                                  !IN  No of CRAY words holding the data
     &    ,PPUNIT                                                       &
                                  !IN  unit no of the PP FILE
     &    ,IWA                                                          &
                                  !IN  WORD address of field to be read
     &    ,ICODE                  !OUT error code
      REAL                                                              &
     &     FIELD(NUM_CRAY_WORDS)  !OUT array holding data
!*---------------------------------------------------------------------
!     Called routines
      EXTERNAL SETPOS,BUFFIN
!*---------------------------------------------------------------------
!     arguments for called routines
      INTEGER                                                           &
     &     LEN_IO                 ! length of data read by BUFFIN
      REAL                                                              &
     &     A_IO                   ! return code from BUFFIN
!    LOCAL VARIABLES
      INTEGER                                                           &
     &     I                                                            &
                                  ! local counter
     &    ,J                                                            &
                                  ! local counter
     &    ,IX                     ! used in the UNIT command

! DEPENDS ON: setpos
      CALL SETPOS(PPUNIT,IWA,ICODE)
! DEPENDS ON: buffin
      CALL BUFFIN(PPUNIT,FIELD,NUM_CRAY_WORDS,LEN_IO,A_IO)

      RETURN
      END SUBROUTINE READ_REC

#endif


