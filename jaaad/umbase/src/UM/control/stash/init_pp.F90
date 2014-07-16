#if defined(C84_1A) || defined(FLDCALC) || defined(FLDMOD)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: INIT_PP  -------------------------------------------------
!LL
!LL  Purpose: Initialises direct access PP files at the start of
!LL           the run.  NB: Sequential PP files need no initialisation.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: D401
!LL
!LL  Project task:
!LL
!LL  External documentation: On-line UM document C61 - Zonal mean
!LL                          calculations.
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE INIT_PP ( FTN_UNIT,FILE_TYPE_LETTER,                   &
     &                     LEN1_LOOKUP,PP_LEN2_LOOKUP,FIXHD,            &
     &                     INTHD,REALHD,LEVDEPC,ROWDEPC,COLDEPC,        &
     &                     LEN_FIXHD,LEN_INTHD,                         &
     &                     LEN_REALHD,LEN1_LEVDEPC,LEN2_LEVDEPC,        &
     &                     LEN1_ROWDEPC,LEN2_ROWDEPC,                   &
     &                     LEN1_COLDEPC,LEN2_COLDEPC,                   &
     &                     PP_LEN_INTHD, PP_LEN_REALHD,                 &
     &                     ICODE,CMESSAGE)

#if defined(C84_1A)
      USE FIELD_BUFF_MOD, ONLY :                                        &
     &    INIT_FXH,                                                     &
     &    ATTACH_FXH,                                                   &
     &    INIT_IPPLOOK,                                                 &
     &    ATTACH_IPPLOOK
#endif

      IMPLICIT NONE
!
      CHARACTER*1                                                       &
     &    FILE_TYPE_LETTER    ! IN  - File type (p-PP, b-bndry)
      INTEGER                                                           &
     &    FTN_UNIT                                                      &
                              ! IN  - Fortran unit number
     &,   LEN1_LOOKUP                                                   &
                              ! IN  - Size of PP header
     &,   PP_LEN2_LOOKUP                                                &
                              ! IN  - Max allowable fields
     &,   LEN_FIXHD                                                     &
                              ! IN    LENGTH OF FIXED CONSTANTS
     &,   LEN_INTHD                                                     &
                              ! IN    LENGTH OF INTEGER CONSTANTS
     &,   LEN_REALHD                                                    &
                              ! IN    LENGTH OF REAL CONSTANTS
     &,   LEN1_LEVDEPC                                                  &
                              ! IN    LENGTH OF 1st Dim of lev depndt
     &,   LEN2_LEVDEPC                                                  &
                              ! IN    LENGTH OF 2nd Dim of lev depndt
     &,   LEN1_ROWDEPC                                                  &
                              ! IN    LENGTH OF 1st Dim of row depndt
     &,   LEN2_ROWDEPC                                                  &
                              ! IN    LENGTH OF 2nd Dim of row depndt 
     &,   LEN1_COLDEPC                                                  &
                              ! IN    LENGTH OF 1st Dim of col depndt
     &,   LEN2_COLDEPC                                                  &
                              ! IN    LENGTH OF 2nd Dim of col depndt
     &,   ICODE                                                         &
                              ! OUT - Error exit code
     &,   PP_LEN_INTHD                                                  &
                              ! IN - Length of PP FILE integer header
     &,   PP_LEN_REALHD       ! IN - Length of PP FILE real header
!
!
      INTEGER                                                           &
     &    FIXHD(LEN_FIXHD)                                              &
                                    ! IN    ARRAY OF FIXED CONSTANTS
     &,   INTHD(LEN_INTHD)                                              &
                                    ! IN    ARRAY OF integer CONSTANTS
     &,   LEVDEPC(LEN1_LEVDEPC*LEN2_LEVDEPC)                            &
                                              ! IN LEV DEP CONSTANTS
     &,   ROWDEPC(LEN1_ROWDEPC*LEN2_ROWDEPC)                            &
                                              ! IN ROW DEP CONSTANTS
     &,   COLDEPC(LEN1_COLDEPC*LEN2_COLDEPC)                            &
                                              ! IN COL DEP CONSTANT                                                                                  
     &,   PP_INTHD(PP_LEN_INTHD)                                        &
                                    ! OUT   ARRAY of integer constants
     &,   PP_LEVDEPC(LEN1_LEVDEPC*LEN2_LEVDEPC)                         &
                                                 ! OUT Level dep cts
     &,   PP_ROWDEPC(LEN1_ROWDEPC*LEN2_ROWDEPC)                         &
                                                 ! OUT Row dep cts
     &,   PP_COLDEPC(LEN1_COLDEPC*LEN2_COLDEPC)  ! OUT Col dep cts

#if defined(C84_1A)
      INTEGER, POINTER :: PP_FIXHD(:)
#else
      INTEGER                                                           &
     &    PP_FIXHD(LEN_FIXHD)       ! OUT   ARRAY of fixed constants
#endif

!
      REAL                                                              &
     &    REALHD(LEN_REALHD)                                            &
                                    ! IN    ARRAY OF REAL CONSTANTS
     &,   PP_REALHD(PP_LEN_REALHD)  ! OUT   ARRAY OF REAL CONSTANTS
!
      CHARACTER*80                                                      &
     &    CMESSAGE            ! OUT - Error message
!
!*----------------------------------------------------------------------
!
!  External subroutines
!
      EXTERNAL SETPOS,IOERROR,POSERROR,BUFFOUT,FLUSH_BUFFER
!
!  Local variables
!
#if defined(C84_1A)
      INTEGER, POINTER   :: IPPLOOK(:,:)
      INTEGER, PARAMETER :: current_io_pe=0

      INTEGER            :: DUMMY
      INTEGER            :: STEP
#else
      INTEGER IPPLOOK(LEN1_LOOKUP,PP_LEN2_LOOKUP)
#endif
!
!dir$ cache_align pp_fixhd, pp_inthd, pp_realhd, pp_levdepc, 
!     pp_rowdepc, pp_coldepc, ipplook
#include "cntl_io.h"
      INTEGER                                                           &
     &       II,JJ,IWA,IX,LEN_IO,START_BLOCK  !
      REAL A_IO
#include "c_mdi.h"
#if defined(C84_1A)
#include "parvars.h"
#endif

!L----------------------------------------------------------------------
!L 1. Reserve space
!L
#if defined(C84_1A)
       NULLIFY(IPPLOOK)
       STEP = 1
       CALL INIT_IPPLOOK(IPPLOOK, FTN_UNIT, LEN1_LOOKUP,                &
     &                   PP_LEN2_LOOKUP, DUMMY, STEP)
#else
      DO 1 II=1,PP_LEN2_LOOKUP
      DO 2 JJ=1,LEN1_LOOKUP
      IPPLOOK(JJ,II)=-99
    2 CONTINUE
    1 CONTINUE
#endif

!L----------------------------------------------------------------------
!L 1.1 Set up FIXED header record for the PP FILE
!L
#if defined(C84_1A)
! Attach fixed length header
      CALL INIT_FXH(FTN_UNIT,LEN_FIXHD)
      CALL ATTACH_FXH(PP_FIXHD,FTN_UNIT)
#endif
      DO 3 II=1,LEN_FIXHD
      PP_FIXHD(II)=FIXHD(II)
    3 CONTINUE
      IF (FILE_TYPE_LETTER == 'p' .OR.                                  &
     &    FILE_TYPE_LETTER == 'c') THEN
        PP_FIXHD(5)=3
      ELSEIF (FILE_TYPE_LETTER == 'b') THEN
        PP_FIXHD(5)=5
      ELSE
        ICODE=100
        CMESSAGE='INIT_PP  : Unknown output file type letter'
        RETURN
      ENDIF
      PP_FIXHD(101)=PP_LEN_INTHD
      PP_FIXHD(105)=PP_FIXHD(100)+PP_FIXHD(101)
      PP_FIXHD(106)=PP_LEN_REALHD
      PP_FIXHD(110)=PP_FIXHD(105)+PP_FIXHD(106)
      PP_FIXHD(111)=LEN1_LEVDEPC
      PP_FIXHD(112)=LEN2_LEVDEPC
      PP_FIXHD(115)=0
      PP_FIXHD(116)=rmdi
      PP_FIXHD(117)=rmdi
      PP_FIXHD(120)=0
      PP_FIXHD(121)=rmdi
      PP_FIXHD(122)=rmdi
      PP_FIXHD(125)=0
      PP_FIXHD(126)=rmdi
      PP_FIXHD(127)=rmdi
      PP_FIXHD(130)=0
      PP_FIXHD(131)=rmdi
      PP_FIXHD(135)=0
      PP_FIXHD(136)=rmdi
      PP_FIXHD(140)=0
      PP_FIXHD(141)=rmdi
      PP_FIXHD(142)=0
      PP_FIXHD(143)=rmdi
      PP_FIXHD(144)=0
      PP_FIXHD(145)=rmdi
      PP_FIXHD(150)=PP_FIXHD(110)+ PP_FIXHD(111)*PP_FIXHD(112)
      IF (LEN2_ROWDEPC > 0) THEN
        PP_FIXHD(115)=PP_FIXHD(110)+ PP_FIXHD(111)*PP_FIXHD(112)
        PP_FIXHD(116)=LEN1_ROWDEPC 
        PP_FIXHD(117)=LEN2_ROWDEPC
        PP_FIXHD(150)=PP_FIXHD(115)+ PP_FIXHD(116)*PP_FIXHD(117)  
      END IF 
      IF (LEN2_COLDEPC > 0) THEN 
        PP_FIXHD(120)=PP_FIXHD(115)+ PP_FIXHD(116)*PP_FIXHD(117)
        PP_FIXHD(121)=LEN1_COLDEPC
        PP_FIXHD(122)=LEN2_COLDEPC
        PP_FIXHD(150)=PP_FIXHD(120)+ PP_FIXHD(121)*PP_FIXHD(122) 
      END IF
      PP_FIXHD(151)=LEN1_LOOKUP
      PP_FIXHD(152)=PP_LEN2_LOOKUP
      pp_fixhd(160)=                                                    &
                         ! make sure the data starts on a sector bndry
     & ((pp_fixhd(150)+pp_len2_lookup*len1_lookup-1+um_sector_size-1)/  &
     & um_sector_size)*um_sector_size+1
!L----------------------------------------------------------------------
!L 1.2 Set up INTEGER constants record for the PP FILE
!L
      IF(PP_FIXHD(5) <= 2) THEN !  set all values initially to MDI
        DO II=1,PP_LEN_INTHD
          PP_INTHD(II)=INTHD(21)
        ENDDO
      ELSE
        DO II=1,PP_LEN_INTHD
          PP_INTHD(II)=IMDI
        ENDDO
      ENDIF

      PP_INTHD(6)=INTHD(6)
      PP_INTHD(7)=INTHD(7)
      PP_INTHD(8)=INTHD(8)
      PP_INTHD(9)=INTHD(9)
      PP_INTHD(10)=INTHD(10)
      PP_INTHD(13)=INTHD(13)
      PP_INTHD(17)=INTHD(17)
      PP_INTHD(24)=INTHD(24)
      PP_INTHD(25)=INTHD(25)
      PP_INTHD(28)=INTHD(28)
!L----------------------------------------------------------------------
!L 1.3 Set up REAL constants record for the PP FILE
!L
      DO II = 1, PP_LEN_REALHD
        PP_REALHD(II) = RMDI   ! Set all values to RMDI initially
      END DO
      PP_REALHD(1)=REALHD(1)
      PP_REALHD(2)=REALHD(2)
      PP_REALHD(3)=REALHD(3)
      PP_REALHD(4)=REALHD(4)
! Set to RMDI for VR      
      IF (LEN2_ROWDEPC > 0 .AND. LEN2_COLDEPC > 0) THEN
        PP_REALHD(1) = RMDI
        PP_REALHD(2) = RMDI     
        PP_REALHD(3) = RMDI    
        PP_REALHD(4) = RMDI       
      ENDIF      
      PP_REALHD(5)=REALHD(5)
      PP_REALHD(6)=REALHD(6)
#if defined(FLDMOD)
      IF(PP_LEN_REALHD >= 17)THEN
#endif
      PP_REALHD(16)=REALHD(16)
      PP_REALHD(17)=REALHD(17)
#if defined(FLDMOD)
      ENDIF
#endif

!L----------------------------------------------------------------------
!L 1.4 Set up LEVEL/ROW/COL DEPENDANT constants record for the PP FILE
!L
      DO 5 II=1,LEN1_LEVDEPC*LEN2_LEVDEPC
      PP_LEVDEPC(II)=LEVDEPC(II)
    5 CONTINUE
      
      DO II=1,LEN1_ROWDEPC*LEN2_ROWDEPC
         PP_ROWDEPC(II)=ROWDEPC(II)
      END DO
      
      DO II=1,LEN1_COLDEPC*LEN2_COLDEPC
         PP_COLDEPC(II)=COLDEPC(II)
      END DO
       
!L----------------------------------------------------------------------
!L 2.1 BUFFER OUT Header Records starting with the FIXED LENGTH
!L
! DEPENDS ON: buffout
      CALL BUFFOUT(FTN_UNIT,PP_FIXHD(1),LEN_FIXHD,LEN_IO,A_IO)
        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_FIXHD) THEN
! DEPENDS ON: ioerror
           CALL IOERROR('bufferout of fixed length header',A_IO,LEN_IO, &
     &                    LEN_FIXHD)
           CMESSAGE='INIT_PP:I/O error'
           ICODE=1
           RETURN
        ENDIF
      START_BLOCK=LEN_FIXHD+1
!L----------------------------------------------------------------------
!L 2.2 BUFFER OUT Integer Constants
!L

      IF(FIXHD(100) >  0) THEN  ! Any integer constants to output ?

! Check for error in file pointers

!        WRITE(6,*)  'START_BLOCK FIXHD(100)'
!        WRITE(6,*)   START_BLOCK
!        WRITE(6,*)   FIXHD(100)
!        WRITE(6,*)   FTN_UNIT
         IF(FIXHD(100) /= START_BLOCK) THEN  ! Check start address
! DEPENDS ON: poserror
            CALL POSERROR('integer constants',START_BLOCK,100,          &
     &      PP_FIXHD(100))
            CMESSAGE='INIT_PP:  Addressing conflict'
            ICODE=2
            RETURN
         END IF

! DEPENDS ON: buffout
         CALL BUFFOUT (FTN_UNIT,PP_INTHD(1),PP_FIXHD(101),LEN_IO,A_IO)

! Check for I/O errors

         IF(A_IO /= -1.0.OR.LEN_IO /= PP_FIXHD(101)) THEN
! DEPENDS ON: ioerror
            CALL IOERROR('buffer out of integer constants',A_IO,LEN_IO  &
     &     ,PP_FIXHD(101))
            CMESSAGE='INIT_PP: I/O Error'
            ICODE=3
            RETURN
         END IF

         START_BLOCK=START_BLOCK+PP_FIXHD(101)

      END IF

!L----------------------------------------------------------------------
!L 2.3 BUFFER OUT Real Constants
!L

      IF(PP_FIXHD(105) >  0) THEN   ! Any real constants to output ?

! Check for error in file pointers

        IF(PP_FIXHD(105) /= START_BLOCK) THEN
! DEPENDS ON: poserror
          CALL POSERROR('real constants',START_BLOCK,100,PP_FIXHD(105))
          CMESSAGE='INIT_PP: Addressing conflict'
          ICODE=4
          RETURN
        END IF

! DEPENDS ON: buffout
        CALL BUFFOUT(FTN_UNIT,PP_REALHD(1),PP_FIXHD(106),LEN_IO,A_IO)

! Check for I/O errors

        IF(A_IO /= -1.0.OR.LEN_IO /= PP_FIXHD(106)) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of real constants',A_IO,LEN_IO       &
     &                 ,PP_FIXHD(106))
          CMESSAGE='INIT_PP: I/O Error'
          ICODE=5
          RETURN
        END IF

        START_BLOCK=START_BLOCK+PP_FIXHD(106)

      END IF

!L----------------------------------------------------------------------
!L 2.4.1 BUFFER OUT Level Dependant Constants.
!L

      IF(PP_FIXHD(112) >  0) THEN ! Any level dependant constants ?

! Check for error in file pointers

         IF(PP_FIXHD(110) /= START_BLOCK) THEN
! DEPENDS ON: poserror
            CALL POSERROR('real constants',START_BLOCK,100,             &
     &                     PP_FIXHD(110))
            CMESSAGE='INIT_PP: Addressing conflict'
            ICODE=6
            RETURN
         END IF

! DEPENDS ON: buffout
         CALL BUFFOUT (FTN_UNIT,PP_LEVDEPC(1)                           &
     &              ,PP_FIXHD(111)*PP_FIXHD(112),LEN_IO,A_IO)

! Check for I/O errors

         IF(A_IO /= -1.0.OR.LEN_IO /= (PP_FIXHD(111)*PP_FIXHD(112)      &
     &        ))THEN
! DEPENDS ON: ioerror
           CALL IOERROR('buffer out of lev dep constants',A_IO,LEN_IO   &
     &            ,PP_FIXHD(111))
           CMESSAGE='INIT_PP: I/O Error'
           ICODE=7
           RETURN
         END IF

         START_BLOCK=START_BLOCK+ PP_FIXHD(111)*PP_FIXHD(112)

      END IF
!L----------------------------------------------------------------------
!L 2.4.2 BUFFER OUT Row Dependant Constants.
!L

      IF(PP_FIXHD(115) >  0) THEN ! Any row dependant constants ?

! Check for error in file pointers

        IF(PP_FIXHD(115) /= START_BLOCK) THEN
! DEPENDS ON: poserror
          CALL POSERROR('real constants',START_BLOCK,100,               &
     &                   PP_FIXHD(115))
          CMESSAGE='INIT_PP: Addressing conflict'
          ICODE=6
          RETURN
        END IF

! DEPENDS ON: buffout
        CALL BUFFOUT (FTN_UNIT,PP_ROWDEPC(1),                           &
     &                PP_FIXHD(116)*PP_FIXHD(117),LEN_IO,A_IO)

! Check for I/O errors

        IF(A_IO /= -1.0.OR.LEN_IO /= (PP_FIXHD(116)*PP_FIXHD(117)       &
     &       ))THEN

! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of row dep constants',A_IO,LEN_IO,   &
     &                  PP_FIXHD(116))
          CMESSAGE='INIT_PP: I/O Error'
          ICODE=7
          RETURN
        END IF

        START_BLOCK=START_BLOCK+ PP_FIXHD(116)*PP_FIXHD(117)

      END IF
!L----------------------------------------------------------------------
!L 2.4.3 BUFFER OUT Col Dependant Constants.
!L

      IF(PP_FIXHD(120) >  0) THEN ! Any col dependant constants ?

! Check for error in file pointers

        IF(PP_FIXHD(120) /= START_BLOCK) THEN
! DEPENDS ON: poserror
          CALL POSERROR('real constants',START_BLOCK,100,               &
     &                   PP_FIXHD(120))
          CMESSAGE='INIT_PP: Addressing conflict'
          ICODE=6
          RETURN
        END IF

! DEPENDS ON: buffout
        CALL BUFFOUT (FTN_UNIT,PP_COLDEPC(1),                           &
     &                PP_FIXHD(121)*PP_FIXHD(122),LEN_IO,A_IO)

! Check for I/O errors

        IF(A_IO /= -1.0.OR.LEN_IO /= (PP_FIXHD(121)*PP_FIXHD(122)       &
     &       ))THEN

! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of col dep constants',A_IO,LEN_IO ,  &
     &                  PP_FIXHD(121))
          CMESSAGE='INIT_PP: I/O Error'
          ICODE=7
          RETURN
        END IF

        START_BLOCK=START_BLOCK+ PP_FIXHD(121)*PP_FIXHD(122)

      END IF      
!L----------------------------------------------------------------------
!L 2.5 BUFFER OUT Lookup Table
!L
!     IWA= 0
!     CALL SETPOS(FTN_UNIT,3,IWA,ICODE)
           IF(PP_FIXHD(152) >  0) THEN

! Check for error in file pointers

             IF(PP_FIXHD(150) /= START_BLOCK) THEN
! DEPENDS ON: poserror
               CALL POSERROR('lookup table',START_BLOCK,100,            &
     &              PP_FIXHD(150))
               CMESSAGE='INIT_PP: Addressing conflict'
               ICODE=8
               RETURN
             END IF

! DEPENDS ON: buffout
      CALL BUFFOUT (FTN_UNIT,                                           &
     &              IPPLOOK,LEN1_LOOKUP*PP_LEN2_LOOKUP,LEN_IO,A_IO)

!
! Check for I/O errors

            IF(A_IO /= -1.0.OR.LEN_IO /= (PP_FIXHD(151)*PP_FIXHD(152))) &
     &          THEN
! DEPENDS ON: ioerror
              CALL IOERROR('buffer out of PP LOOKUP TABLE ',A_IO,LEN_IO &
     &            ,PP_FIXHD(152))
              CMESSAGE='INIT_PP: I/O Error'
              ICODE=9
              RETURN
            END IF
!
! Clear file buffer : force last buffer to be written to file
!  to avoid problems with continuation runs following hard failures.
!
      CALL FLUSH_BUFFER(FTN_UNIT,ICODE)
      IF(ICODE /= 0) THEN
         CMESSAGE='INIT_PP: Problem flushing buffer'
         ICODE=10
         RETURN
      ENDIF
!
            START_BLOCK=START_BLOCK+(PP_FIXHD(151)*PP_FIXHD(152))
#if defined(C84_1A)
!
! If we are the current I/O PE, we need to update our copy
! of the LOOKUP Table disk address
!
      STEP = 2
      CALL INIT_IPPLOOK(IPPLOOK, FTN_UNIT, DUMMY, DUMMY,                &
     &                  PP_FIXHD(150) - 1, STEP)
      NULLIFY(PP_FIXHD)
#endif

          END IF
      RETURN
      END SUBROUTINE INIT_PP
#endif
