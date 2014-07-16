

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE READ_FLH --------------------------------------
!LL
!LL  Written by D. Robinson 17/06/92
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL
!LL  Programming standard: Unified Model Documentation Paper No 3
!LL                        Version No 4  5/2/92
!LL
!LL  System component: R30
!LL
!LL  System task: F3
!LL
!LL  Purpose:
!LL           Reads in the fixed length header from file attached to
!LL           unit NFTIN.
!LL
!LL  Documentation:
!LL           Unified Model Documentation Paper No F3
!LL           Version No 5 9/2/90
!LL
!LL------------------------------------------------------------
!*L Arguments:-------------------------------------------------
      SUBROUTINE READ_FLH (NFTIN,FIXHD,LEN_FIXHD,                       &
     &           ICODE,CMESSAGE)

      IMPLICIT NONE

      INTEGER                                                           &
     & NFTIN                                                            &
                        ! IN  Unit no of dump
     &,LEN_FIXHD                                                        &
                        ! IN  Length of fixed length header
     &,FIXHD(LEN_FIXHD) ! OUT Fixed length header

      INTEGER  ICODE  ! OUT Return code; successful=0, error > 0
      CHARACTER*(80)                                                    &
     & CMESSAGE       !OUT Error message if ICODE > 0

! Local arrays:------------------------------------------------
! None
! -------------------------------------------------------------
! External subroutines called:---------------------------------
      EXTERNAL IOERROR,BUFFIN
! Local variables:---------------------------------------------
      INTEGER LEN_IO
      REAL A
! -------------------------------------------------------------

      ICODE=0
      CMESSAGE=' '

!L 1. Buffer in fixed length header record

! DEPENDS ON: buffin
      CALL BUFFIN (NFTIN,FIXHD(1),LEN_FIXHD,LEN_IO,A)

!L 2. Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= LEN_FIXHD)THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of fixed length header',A,LEN_IO        &
     &               ,LEN_FIXHD)
        CMESSAGE='READ_FLH: I/O error'
        ICODE=1
      ENDIF

      RETURN
      END SUBROUTINE READ_FLH
