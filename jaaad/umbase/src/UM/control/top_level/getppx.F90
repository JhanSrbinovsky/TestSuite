#if defined(CONTROL) || defined(UTILIO) || defined(FLDOP)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Reads PPXREF file into "look-up" arrays
!
!  Subroutine Interface:

      SUBROUTINE GETPPX(NFTPPXREF,NFTSTMSTU,StmsrNam,RowNumber,         &
#include "argppx.h"
     &                       ErrorStatus,CMESSAGE)
      IMPLICIT NONE
!
!  Description:
!    Reads records from PPXREF file into arrays PPXI (for integer data)
!    and PPXC (for character data, i.e. name of diagnostic/prognostic).
!    The entire PPXREF file is read in (non-null records only).
!
!  Method:
!    Uses routines SETPOS and BUFFIN - these employ Cray-specific code
!
!  Current code owner: S.J.Swarbrick
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

!  Subroutine arguments

!    Scalar arguments with intent(in):
      INTEGER        NFTPPXREF   ! Unit no. for PPXREF file
      INTEGER        NFTSTMSTU   ! Unit no. for user ppxref files
      CHARACTER*13    StmsrNam    ! Names of stash master files

!    Array arguments with intent(out):
      CHARACTER*80 CMESSAGE    ! Error return message

!    Error status:
      INTEGER      ErrorStatus ! Error return code

!  Local scalars:
      INTEGER      I,J,IE,ID,II  ! Loop counters
      INTEGER      hashcount
      INTEGER      IFIL,IREC     ! Do.
      INTEGER      LEN_IO        ! No. words read on each CALL BUFFIN
      REAL         STATUS        ! Error return code from BUFFIN
      INTEGER      IOSTATUS
      CHARACTER*80 UpsmFile      ! Full pathname for user psm files
      CHARACTER*256 STASH_MSTR    ! Do. STASH master files
      CHARACTER*(*)Routine       ! Subroutine name
      PARAMETER   (Routine = 'GETPPX')
      CHARACTER*1  CHAR1
      INTEGER      Im_index      !
      INTEGER      Im_ident      !
      INTEGER      Section       !
      INTEGER      Item          !
      INTEGER      LModel  ,DM
      INTEGER      LSection,DS
      INTEGER      LItem   ,DI
      INTEGER      USTrow
      INTEGER      RowNumber     ! Row no. counter for PPXI, PPXC arrays
      INTEGER      FirstBlank    ! Used to append Upsm file name to dir
      INTEGER      RI            ! Row index
      INTEGER      NU_recs       ! No. of records in a user psm file
      LOGICAL      OVERWRITE ! Set T if a system stash master record
                             !  is being overwritten by a user rec.
!  Local arrays:
!  WARNING: must have PPXREF_CHARLEN=4*PPX_CHARWORD
!           to avoid overwriting
      CHARACTER DNAM (PPXREF_CHARLEN) ! For character part of ppx rec
      INTEGER   CODES(PPXREF_CODELEN) ! For integer part of ppx record
      INTEGER   IMASK(20)             ! For ver mask in user psm
#if defined(T3E)

      common/shmem_getppx_c1/ dnam
!dir$ cache_align /shmem_getppx_c1/
      common/shmem_getppx_c2/ char1
!dir$ cache_align /shmem_getppx_c2/
      common/shmem_getppx/ codes, iostatus, nu_recs
!dir$ cache_align /shmem_getppx/
!
      integer shmem_n_pes, msg, info, nproc, shmem_my_pe, mype
!
#endif

!  Function and subroutine calls:
      EXTERNAL  READSTM
!
!- End of header -------------------------------------------------------
!
      ErrorStatus = 0
      NU_recs     = 0
      IOStatus   =0
!----------------------------------------------------------------------
! Check that the no. of requested diagnostics does not exceed max
! defined in comdecks VERSION and PPXLOOK.
!
      IF ( (ppxRecs  >   NDIAGP) .OR. (ppxRecs  >   NUM_DIAG_MAX) )     &
     &THEN
        WRITE(6,*) 'ERROR: no. of diags. requested exceeds max'
        WRITE(6,*) 'ppxRecs=',ppxRecs,' NDIAGP=',NDIAGP,                &
     &             ' NUM_DIAG_MAX=',NUM_DIAG_MAX
        Errorstatus=104
        CMESSAGE= 'GETPPX: ppxRecs GT (NDIAGP or NUM_DIAG_MAX)'
        GO TO 9999
      END IF
!----------------------------------------------------------------------

#if defined(T3E)
      mype=shmem_my_pe()
      nproc=shmem_n_pes()
#endif
      IF (NFTPPXREF == 22) THEN
!----------------------------------------------------------------------
!Read in records from STASHmaster for current internal model
!----------------------------------------------------------------------
!Open STASHmaster file for current internal model
!  Get directory name for STASHmaster & append rest of filename
#if defined(T3E)
        stash_mstr='empty '
        if(mype == 0) CALL GET_FILE(NFTPPXREF,STASH_MSTR,256,            &
     &   ErrorStatus)
#else
        CALL GET_FILE(NFTPPXREF,STASH_MSTR,256,ErrorStatus)
#endif
        FirstBlank = 0
        DO I = 1,256
          IF (STASH_MSTR(I:I) == ' '.AND.FirstBlank == 0)               &
     &                                   FirstBlank=I
        END DO
        STASH_MSTR(FirstBlank:FirstBlank)='/'
        STASH_MSTR(FirstBlank+1:FirstBlank+13)=StmsrNam
#if defined(T3E)
        if(mype == 0) OPEN(UNIT=NFTPPXREF,FILE=STASH_MSTR,              &
     &   IOSTAT=IOStatus)

        msg=7030
        info=0
        call gc_ibcast(msg, 1, 0, nproc, info, IOStatus)

#else
        OPEN(UNIT=NFTPPXREF,FILE=STASH_MSTR,IOSTAT=IOStatus)
#endif
        IF(IOStatus /= 0) THEN
          WRITE(6,*) 'ERROR in routine GETPPX'
          WRITE(6,*)                                                    &
     &   'CANNOT OPEN STASHmaster FILE, IOSTATUS=',IOStatus
          WRITE(6,*) 'UNIT=',NFTPPXREF,' FILE=',STASH_MSTR
          ErrorStatus=100
          CMESSAGE=' GETPPX: ERROR OPENING STASHmaster'
          GOTO 9999
        END IF


#if defined(T3E)
 100    continue
        if(mype == 0) READ(NFTPPXREF,'(A1)') CHAR1
!
        msg=7033
        info=0
        call gc_cbcast(msg, 1, 0, nproc, info, char1)
!
#else
 100    READ(NFTPPXREF,'(A1)') CHAR1
#endif
        IF (CHAR1 == '1') THEN
!Read block of records
#if defined(T3E)
          if(mype == 0) then
            BACKSPACE NFTPPXREF
! DEPENDS ON: readstm
            CALL READSTM(IMASK,DNAM,CODES,NFTPPXREF,                    &
     &       ErrorStatus,CMESSAGE)
          endif

          msg=7031
          info=0
          call gc_ibcast(msg, ppxref_codelen, 0, nproc, info, codes)
          msg=7032
          info=0
          call gc_cbcast(msg, ppxref_charlen, 0, nproc, info, dnam)
#else
          BACKSPACE NFTPPXREF
! DEPENDS ON: readstm
          CALL READSTM(IMASK,DNAM,CODES,NFTPPXREF,ErrorStatus,CMESSAGE)
#endif
          Im_ident = CODES(ppx_model_number)
          Section  = CODES(ppx_section_number)
          Item     = CODES(ppx_item_number)
          IF (Im_ident == -1) THEN
!End of file reached
#if defined(T3E)
            if(mype == 0) CLOSE(UNIT=NFTPPXREF)
#else
            CLOSE(UNIT=NFTPPXREF)
#endif
            GO TO 9999
          END IF
          Im_index= INTERNAL_MODEL_INDEX(Im_ident)
!   Increment row number
          RowNumber = RowNumber + 1
! Assign value to PPXPTR element corresponding to this record
#if defined(CONVIEEE) || defined(CONVPP)           \
 || defined(CUMF) || defined(PUMF) || defined(MERGE)
          PPXPTR(Im_ident,Section,Item) = RowNumber
#else
          PPXPTR(Im_index,Section,Item) = RowNumber
#endif
!   Transfer data from ppx record to look-up arrays
          DO I=1,PPXREF_CHARLEN
            PPXC(RowNumber,I)=DNAM(I)
          END DO
          DO I=1,PPXREF_CODELEN
            PPXI(RowNumber,I)=CODES(I)
          END DO
!   Set row index - indicates values of model,sec,item for this row
          RowIndex  (RowNumber)=  Im_ident*100000                       &
     &                          + Section *1000                         &
     &                          + Item
!   Set flag to indicate record originated from ppxref file
          OriginFlag(RowNumber)='P'
          IF (RowNumber  >   ppxRecs) THEN
            WRITE(6,*) 'Error in GETPPX:'
            WRITE(6,*)                                                  &
     &    ' PPXI row number exceeds total no. of ppx records ',         &
     &      RowNumber
            GO TO 9999
          END IF
          GO TO 100  ! Back to READ
        ELSE
! Skip to next line
          GO TO 100
        END IF
      ELSE         ! NFTPPXREF /= 1
! ----------------------------------------------------------
! Insert user-defined diagnostics into ppxref look-up arrays
! ----------------------------------------------------------

#if defined(PUMF) || defined(CUMF) || defined(CONVIEEE)                \
 || defined(MERGE) || defined(CONVPP)
      IF (NRECS_USTASH(1) >  0) THEN
#elif defined(FLDOP)
      IF (NRECS_USTASH(1) >  0) THEN
#else
      IF (NRECS_USTASH >  0) THEN
#endif
! There are user diagnostic records
      ErrorStatus=0
      IOStatus   =0
#if defined(PUMF) || defined(CUMF) || defined(CONVIEEE)                \
 || defined(MERGE) || defined(CONVPP)
#elif defined(FLDOP)
#else
! Get directory name for Upsm files
#if defined(T3E)
      upsmfile='empty '
      if(mype == 0) CALL GET_FILE(NFTSTMSTU,UpsmFile,80,                &
     & ErrorStatus)
#else
      CALL GET_FILE(NFTSTMSTU,UpsmFile,80,ErrorStatus)
#endif
      FirstBlank = 0
      DO I = 1,80
        IF (UpsmFile(I:I) == ' '.AND.FirstBlank == 0) FirstBlank=I
      END DO
#endif

! Loop over user pre-stash master files
      DO IFIL = 1,N_USTASH
#if defined(PUMF) || defined(CUMF) || defined(CONVIEEE)                \
 || defined(MERGE) || defined(CONVPP)
      UpsmFile=USTSFILS(IFIL)
#elif defined(FLDOP)
      UpsmFile=USTSFILS(IFIL)
#else
        UpsmFile(FirstBlank  :FirstBlank  )='.'
        UpsmFile(FirstBlank+1:FirstBlank+8)=USTSFILS(IFIL)

#endif
!   Open user stash master file
#if defined(T3E)
        if(mype == 0) OPEN(UNIT=NFTSTMSTU,FILE=UpsmFile,                &
     &   IOSTAT=IOStatus)

        msg=7040
        info=0
        call gc_ibcast(msg, 1, 0, nproc, info, IOStatus)

#else
        OPEN(NFTSTMSTU,FILE=UpsmFile,IOSTAT=IOStatus)
#endif
        IF(IOStatus /= 0) THEN
          WRITE(6,*) 'CANNOT OPEN USER PPXREF FILE.IOSTATUS=',          &
     &                                             IOStatus
          WRITE(6,*) 'UNIT=',NFTSTMSTU,' FILE=',UpsmFile
          ErrorStatus=100
          CMESSAGE=' GETPPX: ERROR OPENING USER PPXREF'
          GOTO 9999
        END IF

#if defined(PUMF) || defined(CUMF) || defined(CONVIEEE)                \
 || defined(MERGE) || defined(CONVPP)
        NU_recs = NRECS_USTASH(IFIL)
#elif defined(FLDOP)
        NU_recs = NRECS_USTASH(IFIL)
#else
!   Read number of records in this file
#if defined(T3E)
        if(mype == 0) READ(NFTSTMSTU,'(I3)') NU_recs

        msg=7050
        info=0
        call gc_ibcast(msg, 1, 0, nproc, info, nu_recs)

#else
        READ(NFTSTMSTU,'(I3)') NU_recs
#endif
#endif

!   Read in records from user pre-stash master file
        DO IREC = 1,NU_recs
!   Initialise OVERWRITE switch
        OVERWRITE=.FALSE.
        hashcount=0
#if defined(T3E)
 200    continue
        if(mype == 0) READ(NFTSTMSTU,'(A1)') CHAR1
!
        msg=7043
        info=0
        call gc_cbcast(msg, 1, 0, nproc, info, char1)
!
#else
 200    READ(NFTSTMSTU,'(A1)') CHAR1
#endif
        IF (CHAR1 /= '1') THEN
          hashcount=hashcount+1
          IF (hashcount >  20) THEN
            Errorstatus=100
            CMESSAGE='INCORRECT FORMAT IN USER STASHmaster FILE'
            WRITE(6,*) 'INCORRECT FORMAT IN USER STASHmaster FILE'
            WRITE(6,*) 'GAP BETWEEN RECORDS TOO LARGE?'
            GO TO 9999
          ELSE
            GO TO 200
          END IF
        ELSE
!Read block of records
#if defined(T3E)
          if(mype == 0) then
            BACKSPACE NFTSTMSTU
! DEPENDS ON: readstm
            CALL READSTM                                                &
     &       (IMASK,DNAM,CODES,NFTSTMSTU,ErrorStatus,CMESSAGE)
          endif

          msg=7041
          info=0
          call gc_ibcast(msg, ppxref_codelen, 0, nproc, info, codes)
          msg=7042
          info=0
          call gc_cbcast(msg, ppxref_charlen, 0, nproc, info, dnam)
#else
          BACKSPACE NFTSTMSTU
! DEPENDS ON: readstm
          CALL READSTM                                                  &
     &   (IMASK,DNAM,CODES,NFTSTMSTU,ErrorStatus,CMESSAGE)
#endif
          Im_ident = CODES(ppx_model_number)
          Section  = CODES(ppx_section_number)
          Item     = CODES(ppx_item_number)

!   Transfer data from ppx record to look-up arrays
!   No. of records extracted from STASHmaster file(s)= RowNumber.
          USTrow    =   0
          DO I=1,RowNumber
            RI      =   RowIndex(I)
            IF(RI <= 0) THEN   ! Check valid Rowindex
               Errorstatus=-1
               Cmessage= Routine//                                      &
     & ':Warning, invalid Rowindex for user STASHmaster record'
               write(6,*) Cmessage
               write(6,*) 'im_ident section item Rowindex=',            &
     &                     im_ident,section,item,RI
            ENDIF   ! Check valid Rowindex

!     Determine values of model,section,item for this row
            IF (USTrow >  0) THEN  ! Position of record found so
                exit               ! exit from I loop over RowNumber
            ENDIF

              LModel  =     RI/100000
              LSection=(RI-(RI/100000)*100000)/1000
              LItem   =(RI-(RI/1000  )*1000  )
!     Check whether previous item is being overwritten
              IF (Im_ident == LModel  .AND.                             &
     &            Section  == LSection.AND.                             &
     &            Item     == LItem        ) THEN
                IF      (OriginFlag(I) == 'P') THEN
                  OVERWRITE=.TRUE.
                  WRITE(6,*) 'MESSAGE FROM ROUTINE GETPPX:'
                  WRITE(6,*)                                            &
     &           'The following PPXREF record has been overwritten by'
                  WRITE(6,*)                                            &
     &           'a record read from a user-STASH master file: '
                  WRITE(6,*) 'Internal Model ',Im_ident,                &
     &           ' Section ',Section,' Item ',Item
                ELSE IF (OriginFlag(I) == 'U') THEN
                  WRITE(6,*) 'ERROR, GETPPX: '
                  WRITE(6,*) 'User diagnostic duplicated'
                  WRITE(6,*) 'Model,Section,Item ',                     &
     &                        Im_ident,Section,Item
                  ErrorStatus=100
                  CMESSAGE='ERROR,GETPPX:user diag duplicated'
                  GO TO 9999
                END IF
              END IF
!     Determine appropriate row number
              IF (LModel   == Im_ident.AND.                             &
     &            LSection == Section .AND.                             &
     &            LItem    == Item    .AND.USTrow == 0) THEN
                USTrow=I    ! Row number found
!     This record will overwrite a pre-existing record
!     Insert new record
                DO IE=1,PPXREF_CHARLEN
                  PPXC(USTrow,IE)=DNAM(IE)
                END DO
                DO IE=1,PPXREF_CODELEN
                  PPXI(USTRow,IE)=CODES(IE)
                END DO
!     Set flag to indicate record originated from user psm file
                OriginFlag(USTrow)='U'
              ELSE IF((LModel   >  Im_ident.AND.USTrow == 0) .OR.       &
     &                (LModel   == Im_ident.AND.                        &
     &                 LSection >  Section .AND.USTrow == 0) .OR.       &
     &                (LModel   == Im_ident.AND.                        &
     &                 LSection == Section .AND.                        &
     &                 LItem    >  Item    .AND.USTrow == 0)) THEN
                USTrow=I    ! Row number found
!     This record will be inserted between two pre-existing records
!     Create spare row - move all subsequent records up by one row
                DO ID = RowNumber+1,USTrow+1,-1
                  DO IE=1,PPXREF_CHARLEN
                    PPXC(ID,IE)=PPXC(ID-1,IE)
                  END DO
                  DO IE=1,PPXREF_CODELEN
                    PPXI(ID,IE)=PPXI(ID-1,IE)
                  END DO
                  RI            =RowIndex  (ID-1)
                  RowIndex  (ID)=RowIndex  (ID-1)
                  OriginFlag(ID)=OriginFlag(ID-1)
!     Determine values of model,section,item for this row
                  DM=     RI/100000
                  DS=(RI-(RI/100000)*100000)/1000
                  DI=(RI-(RI/1000  )*1000  )
!     Increment PPXPTR for record moved up
#if defined(CONVIEEE) || defined(CONVPP)           \
 || defined(CUMF) || defined(PUMF) || defined(MERGE)
                  PPXPTR(DM,DS,DI)=PPXPTR(DM,DS,DI)+1
#else
                  Im_index=INTERNAL_MODEL_INDEX(DM)
                  PPXPTR(Im_index,DS,DI)=PPXPTR(Im_index,DS,DI)+1
#endif
                END DO
!     Insert new record
                DO IE=1,PPXREF_CHARLEN
                  PPXC(USTrow,IE)=DNAM(IE)
                END DO
                DO IE=1,PPXREF_CODELEN
                  PPXI(USTRow,IE)=CODES(IE)
                END DO
!     Set row index - indicates model,sec,item for this row
                RowIndex  (USTrow)=  Im_ident*100000                    &
     &                             + Section *1000                      &
     &                             + Item
!     Set flag to indicate record originated from user psm file
                OriginFlag(USTrow)='U'
!     Set PPXPTR for the new record
#if defined(CONVIEEE) || defined(CONVPP)           \
 || defined(CUMF) || defined(PUMF) || defined(MERGE)
                PPXPTR(Im_ident,Section,Item)=USTrow
#else
                Im_index=INTERNAL_MODEL_INDEX(Im_ident)
                PPXPTR(Im_index,Section,Item)=USTrow
#endif


            ELSE IF (I == RowNumber) THEN ! last record and no match
                                          ! so append
!     This record will be added after all pre-existing records
              USTrow = I + 1      ! Set extra entry beyond current last
!     Add new record
              DO IE=1,PPXREF_CHARLEN
                PPXC(USTrow,IE)=DNAM(IE)
              END DO
              DO IE=1,PPXREF_CODELEN
                PPXI(USTrow,IE)=CODES(IE)
              END DO
!     Set row index - indicates model,sec,item for this row
              RowIndex  (USTrow)=  Im_ident*100000                      &
     &                           + Section *1000                        &
     &                           + Item
!     Set flag to indicate record originated from user psm file
              OriginFlag(USTrow)='U'
!     Set PPXPTR for the new record
#if defined(CONVIEEE) || defined(CONVPP)           \
 || defined(CUMF) || defined(PUMF) || defined(MERGE)
              PPXPTR(Im_ident,Section,Item)=USTrow
#else
              Im_index=INTERNAL_MODEL_INDEX(Im_ident)
              PPXPTR(Im_index,Section,Item)=USTrow
#endif
            END IF

          END DO  ! I=1,RowNumber : loop over current set of records

!     Increment RowNumber as UserSTASH record has been added.
!     don't increment it if a standard record has been overwritten.
        IF (.NOT.OVERWRITE) THEN
        RowNumber = RowNumber + 1
        END IF
        END IF        ! hashcount
        END DO        ! Loop over IREC recs in upsm file
      END DO          ! Loop over user psm files
      END IF          ! NRECS_USTASH >  0
#if !defined(UTILIO) && !defined(FLDOP)
! Copy user pre-stash master records to storage arrays -
!   for passing into to U_MODEL
! Note: OriginFlag will be compressed to requested items only
!  at the end of routine STASH_PROC (before used in GETPPX_PART)
      IF (NRECS_USTASH >  0) THEN
      RowNumber = 1
      DO I = 1,ppxRecs
        IF (OriginFlag(I) == 'U') THEN
          DO IE=1,PPXREF_CHARLEN
            PPXC_U(RowNumber,IE)=PPXC(I,IE)
          END DO
          DO IE=1,PPXREF_CODELEN
            PPXI_U(RowNumber,IE)=PPXI(I,IE)
          END DO
          RowNumber=RowNumber+1
        END IF
      END DO
      END IF
#endif

      END IF  !NFT == 22 (Standard STASHmstr or user STASHmstr)

 9999 CONTINUE
      RETURN
      END SUBROUTINE GETPPX
#endif
