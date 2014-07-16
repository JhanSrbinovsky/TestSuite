#if defined(CONTROL) || defined(UTILIO) || defined(FLDOP)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL   SUBROUTINE HDPPXRF ---------------------------------------------
!LL
!LL   PROGRAM TO READ THE HEADER RECORD OF THE PPXREF FILE
!LL   CHECK THE VALUES AND RETURN THE FILE DIMENSIONS
!LL
!LL   AUTHOR            M.J.CARTER
!LL
!LL   TESTED UNDER CFT77 ON OS 5.1
!LL
!LL
!LL   PROGRAMMING STANDARD  UMDP 4
!LL
!LL   LOGICAL COMPONENT  R911
!LL
!LL   PROJECT TASK: C4
!LL
!LL   EXTERNAL DOCUMENT C4
!LL
!LLEND---------------------------------------------------------------
      SUBROUTINE HDPPXRF(NFT,StmsrNam,ppxRecs,ICODE,CMESSAGE)
      IMPLICIT NONE
      INTEGER NFT,NFTU               !IN:  UNIT NUMBER FOR FILE
      CHARACTER*(80) CMESSAGE        !OUT: ERROR RETURN MESSAGE
      INTEGER ICODE                  !OUT: ERROR RETURN CODE

#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "cstash.h"
#include "lenfil.h"

! Local parameters:
      INTEGER, PARAMETER ::                                             &
     &  no_version = -9999 ! Value returned by functio get_umversion if
                           ! environment variable VN not set

! Local Scalars
      INTEGER :: get_um_version
      INTEGER LEN_IO
      INTEGER IU          ! Local - unit no. for stash control file
      INTEGER I
      INTEGER Int_Model_No
      INTEGER FirstBlank
      CHARACTER*13  StmsrNam
      CHARACTER*256 STASH_MSTR
      CHARACTER*1   CHAR1
      INTEGER       IOStatus

      character*8 c_um_version  !UM version as string
      character*8 c_stm_version !STASHmaster version string
      integer um_version,                                               &
                                !Version of UM
     &          um_revision     !Revision of UM
      integer stm_version,                                              &
                                !Version of STASHmaster file
     &  stm_revision            !Revision of STASHmaster file
      integer   ocode           !Copy of the input value of ICODE
      logical found_version     !Indicates presence of STM version
      REAL STATUS
      INTEGER RECORD(PPX_RECORDLEN)
#if defined(T3E)
      integer shmem_n_pes, msg, info, nproc, shmem_my_pe, mype
      common/shmem_hdppxrf/ IOStatus, found_version, stm_version
!dir$ cache_align /shmem_hdppxrf/
#endif
      IOStatus=0
!    Check if an error has already been encountered, and get out
!    if it has.
      ocode = 0
      IF (icode  >   0) then
         goto 9999
      ELSE IF (icode  <   0)then
         ocode = icode
         icode = 0
      END IF


#if defined(T3E)
      mype=shmem_my_pe()
      nproc=shmem_n_pes()
#endif
      IF (NFT == 22) THEN
#if defined(T3E)
        if(mype == 0) then
        stash_mstr='empty '
#endif
!Open STASHmaster file for current internal model
!  Get directory name for STASHmaster & append rest of filename
        CALL GET_FILE(NFT,STASH_MSTR,256,ICODE)
        FirstBlank = 0
        DO I = 1,256
          IF (STASH_MSTR(I:I) == ' '.AND.FirstBlank == 0)               &
     &                                   FirstBlank=I
        END DO
        STASH_MSTR(FirstBlank:FirstBlank)='/'
        STASH_MSTR(FirstBlank+1:FirstBlank+13)=StmsrNam
        OPEN(UNIT=NFT,FILE=STASH_MSTR,IOSTAT=IOStatus)
        write(6,*) '!!!! STASH_MSTR ',STASH_MSTR
#if defined(T3E)
        endif
!
        msg=7001
        info=0
        call gc_ibcast(msg, 1, 0, nproc, info, IOStatus)
#endif

        IF (IOStatus /= 0) THEN
          CMESSAGE=                                                     &
     &   'Error opening STASHmaster file, routine HDPPXRF'
          WRITE(6,*)                                                    &
     &   'HDPPXRF: Fortran Error Response = ',IOStatus,                 &
     &   ' Opening STASHmaster file ',StmsrNam
          ICODE=100
          GO TO 9999
        ENDIF

!    Get the UM version from the environment variable $VN.
! DEPENDS ON: get_um_version
      um_version = get_um_version()

      IF ( um_version  /=  no_version ) THEN


!     Now check through the header section of the STASHmaster
!     file looking for H3
      found_version = .false.
#if defined(T3E)
        if (mype == 0) then
#endif
      READ (nft, '(A1)') char1
      DO WHILE (char1  ==  'H' .or. char1  ==  '#')
         IF (char1  ==  'H') THEN
            BACKSPACE nft
            READ (nft, '(1X, A1)') char1
            IF (char1  ==  '3') THEN
!     This line starts with H3 and should
!     indicate the STASHmaster version. The line should look like
!     H3| UM_VERSION=4.3
               found_version = .true.
               BACKSPACE nft
               READ (nft, '(15x,a8)') c_stm_version
               READ (c_stm_version, '(i1,1x,i1)')                       &
     &              stm_version, stm_revision
               stm_version = stm_version*100 + stm_revision
!     Now perform the check against the UM version
               IF (stm_version  /=  um_version) then
#if defined(T3E)
!--in MPP mode, defer setting the variables until all PE's can
                  icode=1
                  go to 9997
#else
                  write (cmessage,*)                                    &
     & 'HDPPXRF : UM version and STASHmaster version differ'
                  write (6,*) 'Version of STASHmaster file ('           &
     &                 ,stm_version,                                    &
     &                 ') does not match UM version ('                  &
     &                 ,um_version,') in file ',StmsrNam
                  icode = 1
                  goto 9999
#endif
               END IF  ! version check
            END IF  ! char1 == '3'
         END IF                 ! char1 == 'H'
         READ (nft, '(a1)') char1
      END DO

#if defined(T3E)
      endif ! if(mype  ==  0)
!
!--in MPP Mode, get the Value of 'icode', 'stm_version',
!  and 'found_version'
 9997 continue
      msg=7007
      iostatus=icode
      call gc_ibcast (msg,3,0,nproc,info,iostatus)
      icode=iostatus
!--check if we generated a failure
      if(icode /= 0) then
        write (cmessage,*)                                              &
     &   'HDPPXRF: UM version and STASHmaster version differ'
        write (6,*) 'Version of STASHmaster file ('                     &
     &   ,stm_version,') does not match UM version ('                   &
     &   ,um_version,') in file ',StmsrNam
        go to 9999
      endif
#endif

      IF (.not. found_version) THEN
         write (6,*)                                                    &
     &        'HDPPXRF : No STASHmaster version available; Unable to'
         write (6,*)                                                    &
     &        'check against UM version for file ',StmsrNam
         cmessage = 'HDPPXRF : No STASHmaster version available'
         icode = -1
      END IF
#if defined(T3E)
      if(mype == 0) then
#endif
!     For safety, rewind to the start of the STASH file.
      rewind (nft)
#if defined(T3E)
      endif ! if(mype  ==  0)
#endif

      END IF
  100   continue
#if defined(T3E)
        if(mype == 0) then
#endif
!Count records - ppxRecs is counter
        READ(NFT,'(A1)') CHAR1
        IF (CHAR1 == '1') THEN
          BACKSPACE NFT
          READ(NFT,'(2X,I5)') Int_Model_No
          IF (Int_Model_No == -1) THEN
!End of file reached
!ppxRecs initialised to 1 before HDPPXRF - so subtract 1 now
            IF (StmsrNam(13:) == 'A') THEN
              IF (INTERNAL_MODEL_INDEX(A_IM) == 1) THEN
                ppxRecs=ppxRecs-1
              END IF
            END IF
            IF (StmsrNam(13:) == 'O') THEN
              IF (INTERNAL_MODEL_INDEX(O_IM) == 1) THEN
                ppxRecs=ppxRecs-1
              END IF
            END IF
            IF (StmsrNam(13:) == 'S') THEN
              IF (INTERNAL_MODEL_INDEX(S_IM) == 1) THEN
                ppxRecs=ppxRecs-1
              END IF
            END IF
            IF (StmsrNam(13:) == 'W') THEN
              IF (INTERNAL_MODEL_INDEX(W_IM) == 1) THEN
                ppxRecs=ppxRecs-1
              END IF
            END IF
            CLOSE(UNIT=NFT)
#if defined(T3E)
            GO TO 9998
#else
            GO TO 9999
#endif
          END IF
          ppxRecs = ppxRecs + 1
          GO TO 100
        ELSE
          GO TO 100
        END IF
#if defined(T3E)
        endif

9998    continue
        iostatus=ppxrecs
        msg=7002
        call gc_ibcast(msg, 1, 0, nproc, info, iostatus)
        ppxrecs=iostatus
        goto 9999

#endif
      ELSE

! Open stash control file and read USTSNUM namelist: number of user
!   stash files and total no. of user stash records
#if defined(CONVIEEE) || defined(CONVPP)           \
 || defined(CUMF) || defined(PUMF) || defined(MERGE) || defined(FLDOP)
! Read USTNUM namelist from unit 5
        IU = 5
#else
! Read USTNUM namelist from unit 4
        IU = 4                               ! Unit number
#if defined(T3E)
        if(mype == 0) then
        file='empty '
#endif
#if defined(IBM)
        file='empty '
#endif
        CALL GET_FILE(IU,FILE,80,icode)      ! Get name for stash file
        OPEN(IU,FILE=FILE,IOSTAT=icode)      ! Open stash file
#if defined(T3E)
        endif
!
        msg=7003
        info=0
        iostatus=icode
        call gc_ibcast(msg, 1, 0, nproc, info, iostatus)
        icode=iostatus
!
#endif
        IF(icode >  0)THEN                   ! Error check
          WRITE(6,*)'HDPPXRF : Failed in OPEN of Stash Control File'
          GOTO 9999
        ELSEIF(icode <  0)THEN
          WRITE(6,*)'HDPPXRF :                                          &
     &           Warning message on OPEN of Stash Control File'
          WRITE(6,*)'IOSTAT= ',icode
        ENDIF
#endif
!Initialisation
#if defined(PUMF) || defined(CUMF) || defined(CONVIEEE)                \
 || defined(MERGE) || defined(CONVPP)
        DO I = 1,20
          NRECS_USTASH(I)=0
        END DO
#elif defined(FLDOP)
        DO I = 1,20
          NRECS_USTASH(I)=0
        END DO
#else
        N_USTASH     = 0
#endif
        NRECS_USTASH = 0
        DO I = 1,20
          USTSFILS(I)='        '
        END DO
! Read namelist
#if defined(T3E)
        if(mype == 0) then
#endif
        READ(IU,USTSNUM)
#if defined(T3E)
        endif
!
        msg=7004
        call gc_ibcast(msg, 1, 0, nproc, info, N_USTASH)
        msg=7005
        call gc_ibcast(msg, 1, 0, nproc, info, NRECS_USTASH)
        msg=7006
        call gc_cbcast(msg, 160, 0, nproc, info, USTSFILS)
#endif
! Add no. of user stash records to ppxRecs
#if defined(PUMF) || defined(CUMF) || defined(CONVIEEE)                \
 || defined(MERGE) || defined(CONVPP)
        DO I=1,N_USTASH
          ppxRecs = ppxRecs + NRECS_USTASH(I)
        END DO
#elif defined(FLDOP)
        DO I=1,N_USTASH
          ppxRecs = ppxRecs + NRECS_USTASH(I)
        END DO
#else
        ppxRecs = ppxRecs + NRECS_USTASH
#endif
      END IF

 9999 CONTINUE
!     If we have found an error, leave it in icode.  If no error
!     occurred then check if the original input value of icode was]
!     non-zero (a previous untrapped error/warning), and copy this
!     back into ICODE before eaving the routine.
      IF (icode  ==  0 .and. ocode  /=  0) then
         icode = ocode
      END IF
      RETURN
      END SUBROUTINE HDPPXRF

#endif
