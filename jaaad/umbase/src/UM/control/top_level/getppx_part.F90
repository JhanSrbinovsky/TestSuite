#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Reads required portion of PPXREF file into "look-up" arrays
!
!  Subroutine Interface:

      SUBROUTINE GETPPX_PART(NFT,NFTU,StmsrNam,Im_ident,RowNumber,      &
#include "argppx.h"
     &                       ErrorStatus,CMESSAGE)
      IMPLICIT NONE
!
!  Description:
!    Reads records from PPXREF file into arrays PPXI (for integer data)
!    and PPXC (for character data, i.e. name of diagnostic/prognostic).
!    Only those ppxref records corresponding to entries in the STASH
!    addresses array IN_S are read in. Also set up pointer array PPXPTR.
!
!  Method:
!    Uses routines SETPOS and BUFFIN - these employ Cray-specific code
!
!  Current code owner: S.J.Swarbrick
!
!  History:
!  Version   Date       Comment
!  =======   ====       =======
!    3.5     Mar 95     Original code.  S.J.Swarbrick
!    4.0     Sept 95                    S.J.Swarbrick
!    4.0     Dec. 95   Check for ppxRecs LE (NDIAGP or NUM_DIAG_MAX)
!                                       (N Farnon)
!    4.1     Apr. 96   Changes associated with new STASHmaster format
!                      S.J.Swarbrick
!    4.5    30/10/97   Read stash data on PE 0 for the T3E
!                      and distribute it.
!                        Author: Bob Carruthers, Cray Research
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
!  Global Variables:
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "c_mdi.h"
#include "cstash.h"
! Declares IN_S
#include "stextend.h"

!  Subroutine arguments
!    Scalar arguments with intent(in):
      INTEGER      NFT,NFTU      ! Unit nos. for STASHmaster files
      CHARACTER*13  Stmsrnam      ! Names of stash master files

!    Array arguments with intent(out):
      CHARACTER*80   CMESSAGE    ! Error return message

!    Error status:
      INTEGER        ErrorStatus ! Error return code

!  Local scalars:
      INTEGER       I,J,K,Model   ! Loop counters
      CHARACTER*256 STASH_MSTR    ! File name for STASH master
      INTEGER       Im_index      ! Internal model index (run dependent)
      INTEGER       Im_ident      ! Internal model identifier (absolute)
      INTEGER       Section,Sec   ! section no.
      INTEGER       Item,Itm      ! item no.
      INTEGER       RowNumber     ! Row no. counter for PPXI, PPXC arrays
      INTEGER       RowNum_U      ! Do. for PPXI_U, PPXC_U (user diags.)
      CHARACTER*36  NAME
      CHARACTER*1   CHAR1
      INTEGER       FirstBlank
      INTEGER       IOStatus

!  Local arrays:
!  WARNING: must have PPXREF_CHARLEN=4*PPX_CHARWORD
!           to avoid overwriting
      CHARACTER DNAM (PPXREF_CHARLEN) ! For char part of ppx record
      INTEGER   CODES(PPXREF_CODELEN) ! For integer part of ppx record
      INTEGER   IMASK(20)
#if defined(MPP) && defined(T3E)

      common/shmem_getppx_c1/ dnam
!dir$ cache_align /shmem_getppx_c1/
      common/shmem_getppx_c2/ char1
!dir$ cache_align /shmem_getppx_c2/
      common/shmem_getppx/ codes, iostatus, Model, Sec, Itm
!dir$ cache_align /shmem_getppx/
!
      integer shmem_n_pes, msg, info, nproc, shmem_my_pe, mype
!
#endif

!  Function and subroutine calls:
      EXTERNAL READSTM
!
!- End of header -------------------------------------------------------
!
      ErrorStatus = 0
      IOStatus=0
!----------------------------------------------------------------------
! Check that the no. of requested diagnostics does not exceed max
! defined in comdecks VERSION and PPXLOOK.
!
      IF ( (ppxRecs  >   NDIAGP) .OR. (ppxRecs  >   NUM_DIAG_MAX) )     &
     &THEN
        WRITE(6,*) 'ERROR: no. of diags. reqested exceeds max'
        WRITE(6,*) 'ppxRecs=',ppxRecs,' NDIAGP=',NDIAGP,                &
     &             ' NUM_DIAG_MAX=',NUM_DIAG_MAX
        Errorstatus=104
        CMESSAGE= 'GTPPXPT1: ppxRecs GT (NDIAGP or NUM_DIAG_MAX)'
        GO TO 9999
      END IF
!----------------------------------------------------------------------
!Open STASHmaster file for current internal model
!  Get directory name for STASHmaster & append rest of filename
#if defined(MPP) && defined(T3E)
      mype=shmem_my_pe()
      nproc=shmem_n_pes()
#endif
#if defined(MPP) && defined(T3E)
      stash_mstr='empty '
      if(mype == 0) CALL GET_FILE(NFT,STASH_MSTR,256,ErrorStatus)
#else
      CALL GET_FILE(NFT,STASH_MSTR,256,ErrorStatus)
#endif
      FirstBlank = 0
      DO I = 1,256
        IF (STASH_MSTR(I:I) == ' '.AND.FirstBlank == 0)                 &
     &                                   FirstBlank=I
      END DO
      STASH_MSTR(FirstBlank:FirstBlank)='/'
      STASH_MSTR(FirstBlank+1:FirstBlank+13)=StmsrNam
#if defined(MPP) && defined(T3E)
      if(mype == 0) OPEN(UNIT=NFT,FILE=STASH_MSTR,IOSTAT=IOStatus)

      msg=7060
      info=0
      call gc_ibcast(msg, 1, 0, nproc, info, IOStatus)

#else
      OPEN(UNIT=NFT,FILE=STASH_MSTR,IOSTAT=IOStatus)
#endif
      IF(IOStatus /= 0) THEN
        WRITE(6,*) 'ERROR in routine GETPPX_PART'
        WRITE(6,*)                                                      &
     & 'CANNOT OPEN STASHmaster FILE, IOSTATUS=',IOStatus
        WRITE(6,*) 'UNIT=',NFT,' FILE=',STASH_MSTR
        ErrorStatus=100
        CMESSAGE=' GETPPX_PART: ERROR OPENING STASHmaster'
        GOTO 9999
      END IF

! Read the required ppxref records into PPXI, PPXC
      Im_index    = INTERNAL_MODEL_INDEX(Im_ident)
      DO Section  = 0,PPXREF_SECTIONS
        DO Item   = 1,PPXREF_ITEMS

! Check whether there is a stash entry
          IF (IN_S(1,Im_ident,Section,Item)  /=  0) THEN
! Assign pointer value
            PPXPTR(Im_index,Section,Item) = RowNumber

!  OriginFlag was compressed down at end of STASH_PROC,
!  to contain only those items requested.

            IF (OriginFlag(RowNumber) == 'U') THEN
!  Record is from user STASHmaster

!  GETPPX saved all userSTASHmaster records, not just
!  those requested, so search for correct record.
              DO  I = 1,NUM_USR_DIAG_MAX
                IF (PPXI_U(I,1) == Im_ident .and.                       &
     &              PPXI_U(I,2) == Section .and.                        &
     &              PPXI_U(I,3) == Item) THEN
!  Correct record found
                  RowNum_U = I
                END IF
              END DO
! Read user ppxref record from transfer arrays
              DO I=1,PPXREF_CHARLEN
                PPXC(RowNumber,I)=PPXC_U(RowNum_U,I)
              END DO
              DO I=1,PPXREF_CODELEN
                PPXI(RowNumber,I)=PPXI_U(RowNum_U,I)
              END DO
              IF ((PPXI(RowNumber,1) /= Im_ident).OR.                   &
     &            (PPXI(RowNumber,2) /= Section ).OR.                   &
     &            (PPXI(RowNumber,3) /= Item    )) THEN
                WRITE(6,*) 'ERROR, GETPPX_PART: '
                WRITE(6,*) 'Inconsistency in user ppxref transfer'
                WRITE(6,*) 'Model,Section,Item: ',                      &
     &                      Im_ident,Section,Item
                ErrorStatus=115
                GO TO 9999
              END IF

            ELSE IF (OriginFlag(RowNumber) == 'P') THEN

! Find appropriate record in STASHmaster file and read it in
#if defined(MPP) && defined(T3E)
 100          continue
              if(mype == 0) READ(NFT,'(A1)') CHAR1
!
              msg=7063
              info=0
              call gc_cbcast(msg, 1, 0, nproc, info, char1)
!
#else
 100          READ(NFT,'(A1)') CHAR1
#endif
              IF (CHAR1 == '1') THEN
#if defined(MPP) && defined(T3E)
                if(mype == 0) then
                  BACKSPACE NFT
                  READ(NFT,'(2X,3(I5,2X))') Model,Sec,Itm
                endif
!
                msg=7066
                info=0
                call gc_ibcast(msg, 3, 0, nproc,info, model)
!
#else
                BACKSPACE NFT
                READ(NFT,'(2X,3(I5,2X))') Model,Sec,Itm
#endif
                IF (Model == -1) THEN
                  WRITE(6,*)                                            &
     &           'GETPPX_PART: End of STASHmaster file ',               &
     &            StmsrNam,' reached'
                  GO TO 1100
                END IF
                IF (Sec == Section .AND. Itm == Item) THEN
!   Correct record found
#if defined(MPP) && defined(T3E)
                  if(mype == 0) then
                    BACKSPACE NFT
! DEPENDS ON: readstm
                    CALL READSTM                                        &
     &               (IMASK,DNAM,CODES,NFT,ErrorStatus,CMESSAGE)
                  endif

                  msg=7061
                  info=0
                  call gc_ibcast(msg, ppxref_codelen, 0, nproc,         &
     &             info, codes)
                  msg=7062
                  info=0
                  call gc_cbcast(msg, ppxref_charlen, 0, nproc,         &
     &             info, dnam)
#else
                  BACKSPACE NFT
! DEPENDS ON: readstm
                  CALL READSTM                                          &
     &           (IMASK,DNAM,CODES,NFT,ErrorStatus,CMESSAGE)
#endif
!   Transfer STASHmaster record to look-up arrays
                  DO I=1,PPXREF_CHARLEN
                    PPXC(RowNumber,I)=DNAM(I)
                  END DO
                  DO I=1,PPXREF_CODELEN
                    PPXI(RowNumber,I)=CODES(I)
                  END DO
                ELSE
                  GO TO 100
                END IF
              ELSE
                GO TO 100
              END IF
            ELSE IF (OriginFlag(RowNumber) /= ' ') THEN
              WRITE(6,*) 'ERROR, GETPPX_PART: INVALID OriginFlag'
              WRITE(6,*) 'Row number, Flag'
              WRITE(6,*) RowNumber, OriginFlag(RowNumber)
                ErrorStatus=135
                GO TO 9999
            END IF

            RowNumber = RowNumber + 1
 1100       CONTINUE

            IF ((RowNumber-1)  >   ppxRecs) THEN
              WRITE(6,*) 'Error in GETPPX_PART:'
              WRITE(6,*)                                                &
     &       ' PPXI row number exceeds total no. of ppx records'
              GO TO 9999
            END IF

          END IF   ! Stash entries
        END DO     ! Items
      END DO     ! Sections

 9999 CONTINUE

#if defined(MPP) && defined(T3E)
      if(mype == 0) CLOSE(UNIT=NFT)
#else
      CLOSE(UNIT=NFT)
#endif
      RETURN
      END SUBROUTINE GETPPX_PART
!- End of Subroutine code ---------------------------------------------
#endif
