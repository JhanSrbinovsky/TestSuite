#if defined(C80_1A) && !defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE WRITDUMP---------------------------------------
!LL
!LL AD, TJ      <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.1  19/02/93  Use FIXHD(12) not FIXHD(1) as Version no in P21BITS
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL   3.2  24/07/93  CHECK ON THE ERROR STATUS FROM BUFFOUT
!LL   3.2  25/05/93  Skip DIAG81 diagnostics for obs files. D. Robinson
!LL   3.3  08/04/94  Check that BUFLEN is long enough
!LL   3.3  22/11/93  Prevent dynamic allocation of zero dimension for
!LL                  BUF. (Possible for obs files) Call PR_LFLD to print
!LL                  logical fields. D Robinson.
!LL   3.5  28/03/95  MPP code : New code for parallel I/O
!LL                                              P.Burton
!LL   4.1  22/05/96  Fixes to MPP code    P.Burton
!     4.1  18/06/96  Changes to cope with changes in STASH addressing
!                    Author D.M. Goddard.
!LL   4.3  17/03/97  Changed name to UM_WRITDUMP and added
!LL                  D1_ADDRESSING arguments, which are passed
!LL                  to write_multi.                     P.Burton
!    4.4    23/07/97  Correct change_decomp error message   P.Burton
!     4.4  25/04/97  Changes to write well-formed records if the
!                    input dumpfile is in that format (almost PP file
!                    format)
!                      Author: Bob Carruthers, Cray Research
!     4.5  5/11/98   Removed check that field size is less than
!                    MaxFieldSize as the arrays are now dynamically
!                    allocated to the required size.         P.Burton
!    4.5    28/10/98  Introduce Single Column Model. J-C Thil.
!    5.0    14/09/99 Rewritten UM_WRITDUMP for vn5.x         P.Burton
!    5.3    08/06/01 Declare LEN_DATA before uisng.  A van der Wal
!    5.5    02/08/00 Modification for parallelisation of WAM.
!                    Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
!    6.0    13/11/03 Buffering of levels for improved SX-6 efficiency.
!                    Also tidy up of non-MPP code. P.Selwood.
!    6.1    17/08/04 Use buffout for big endian systems to utilise
!                    C buffering. JC Rioual, P Selwood.
!    6.2    06/12/05 Removed DIAG80 CPP define. T.Edwards
!LL
!LL  Programming standard: Unified Model Documentation Paper No 3
!LL                        Version No 1 15/1/90
!LL
!LL  Logical component: R30
!LL
!LL  Project task: F3
!LL
!LL  Purpose: Writes out model dump on unit NFTOUT and checks model
!LL           and dump dimensions for consistency.
!LL
!LL  Documentation: Unified Model Documentation Paper No F3
!LL                 Version No 5 9/2/90
!LLEND---------------------------------------------------------
!
!*L Arguments:-------------------------------------------------
      SUBROUTINE UM_WRITDUMP(NFTOUT,FIXHD,LEN_FIXHD,                    &
     &                    INTHD,LEN_INTHD,                              &
     &                    REALHD,LEN_REALHD,                            &
     &                    LEVDEPC,LEN1_LEVDEPC,LEN2_LEVDEPC,            &
     &                    ROWDEPC,LEN1_ROWDEPC,LEN2_ROWDEPC,            &
     &                    COLDEPC,LEN1_COLDEPC,LEN2_COLDEPC,            &
     &                    FLDDEPC,LEN1_FLDDEPC,LEN2_FLDDEPC,            &
     &                    EXTCNST,LEN_EXTCNST,                          &
     &                    DUMPHIST,LEN_DUMPHIST,                        &
     &                    CFI1,LEN_CFI1,                                &
     &                    CFI2,LEN_CFI2,                                &
     &                    CFI3,LEN_CFI3,                                &
     &                    LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP,               &
     &                    MPP_LOOKUP,MPP_LEN1_LOOKUP,                   &
     &                    BUFLEN,                                       &
#include "argppx.h"
     &                    SUBMODEL_ID,                                  &
     &                    N_OBJS_D1,D1_ADDR,                            &
     &                    LEN_DATA,D1)

      IMPLICIT NONE

      INTEGER                                                           &
     & NFTOUT                                                           &
                     !IN Unit no of dump
     &,LEN_FIXHD                                                        &
                     !IN Length of fixed length header
     &,LEN_INTHD                                                        &
                     !IN Length of integer header
     &,LEN_REALHD                                                       &
                     !IN Length of real header
     &,LEN1_LEVDEPC                                                     &
                     !IN 1st dim of level dep consts
     &,LEN2_LEVDEPC                                                     &
                     !IN 2ndt dim of level dep consts
     &,LEN1_ROWDEPC                                                     &
                     !IN 1st dim of row dep consts
     &,LEN2_ROWDEPC                                                     &
                     !IN 2nd dim of row dep consts
     &,LEN1_COLDEPC                                                     &
                     !IN 1st dim of column dep consts
     &,LEN2_COLDEPC                                                     &
                     !IN 2nd dim of column dep consts
     &,LEN1_FLDDEPC                                                     &
                     !IN 1st dim of field dep consts
     &,LEN2_FLDDEPC                                                     &
                     !IN 2nd dim of field dep consts
     &,LEN_EXTCNST                                                      &
                     !IN Length of extra constants
     &,LEN_DUMPHIST                                                     &
                     !IN Length of history block
     &,LEN_CFI1                                                         &
                     !IN Length of comp field index 1
     &,LEN_CFI2                                                         &
                     !IN Length of comp field index 2
     &,LEN_CFI3                                                         &
                     !IN Length of comp field index 3
     &,LEN1_LOOKUP                                                      &
                     !IN 1st dim of lookup
     &,LEN2_LOOKUP                                                      &
                     !IN 2nd dim of lookup
     &,MPP_LEN1_LOOKUP                                                  &
                       !IN 1st dim of MPP lookup

     &,SUBMODEL_ID                                                      &
                     !IN submodel of dump
     &,N_OBJS_D1     !IN number of objects (3D fields) in D1

      INTEGER                                                           &
     & FIXHD(LEN_FIXHD)                                                 &
                        !IN Fixed length header
     &,INTHD(LEN_INTHD)                                                 &
                        !IN Integer header
     &,LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)                                  &
                                       !IN PP lookup tables
     &,MPP_LOOKUP(MPP_LEN1_LOOKUP,LEN2_LOOKUP)                          &
                                               ! OUT
     &,CFI1(LEN_CFI1+1)                                                 &
                        !IN Compressed field index no 1
     &,CFI2(LEN_CFI2+1)                                                 &
                        !IN Compressed field index no 2
     &,CFI3(LEN_CFI3+1) !IN Compressed field index no 3

      INTEGER                                                           &
     & LEN_DATA       !IN Length of real data
      REAL                                                              &
     & REALHD(LEN_REALHD)                                               &
                          !IN Real header
     &,LEVDEPC(1+LEN1_LEVDEPC*LEN2_LEVDEPC)                             &
                                            !IN Lev dep consts
     &,ROWDEPC(1+LEN1_ROWDEPC*LEN2_ROWDEPC)                             &
                                            !IN Row dep consts
     &,COLDEPC(1+LEN1_COLDEPC*LEN2_COLDEPC)                             &
                                            !IN Col dep consts
     &,FLDDEPC(1+LEN1_FLDDEPC*LEN2_FLDDEPC)                             &
                                            !IN Field dep consts
     &,EXTCNST(LEN_EXTCNST+1)                                           &
                                !IN Extra constants
     &,DUMPHIST(LEN_DUMPHIST+1)                                         &
                                !IN History block
     &,D1(LEN_DATA)     !IN Real equivalence of data block

! Parameters required for dimensioning the D1_ADDR array
#include "d1_addr.h"

      INTEGER                                                           &
     &  D1_ADDR(D1_LIST_LEN,N_OBJS_D1)  ! IN D1 addressing info

      INTEGER                                                           &
     & BUFLEN         !IN Maximum length of single field in dump

! Comdecks/common blocks/parameters

#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "clookadd.h"
#include "c_mdi.h"
#include "decomptp.h"
#include "parvars.h"
#include "cntl_io.h"

! Diagnostic message switch
#include "cprintst.h"

! Local variables
      INTEGER                                                           &
     &  start_block                                                     &
                           ! first word of field data in dump (unused)
     &, object_index                                                    &
                           ! pointer to entry in D1_ADDR
     &, level                                                           &
                           ! level number of multi-level field
     &, number_of_fields                                                &
                           ! total number of fields to write out
     &, field_start                                                     &
                           ! start address of a field in the file
     &, data_size                                                       &
                           ! number of words of data on disk for a field
     &, data_write_size                                                 &
                           ! total number of words to read for a field
     &, len_io                                                          &
                           ! number of words actually written
     &, orig_decomp                                                     &
                           ! original decomposition type
     &, local_len                                                       &
                           ! length of local field from buffout
     &, address                                                         &
                           ! address of field in local D1 array
     &, k,i                                                             &
                           ! loop indicies
     &, number_of_data_words_in_memory                                  &
                                       ! unused return argument
     &, number_of_data_words_on_disk                                    &
                                       ! unused return argument
     &, disk_address       ! unused return argument

      INTEGER :: buffer_start     ! start field in output buffer
      INTEGER :: buffer_pos       ! current position in output buffer
      INTEGER :: buffered_fields  ! count of buffered fields

      LOGICAL                                                           &
     &  packed_field       ! TRUE if a field has been packed
                           !      to 32 bits

      REAL    :: IOSTAT    ! return from buffout_single

! Error reporting
      INTEGER       ICODE       ! =0 normal exit; >0 error exit
      CHARACTER*256 Cmessage    ! Error message
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='UM_WRITDUMP')

#if defined(NEC)
! NEC works best with lots of fields buffered up together for
! output - smaller default for other machines to save memory.
! May be worth fiddling with this on a platform-by-platform basis.
! Note that due to not forcing a flush of the buffer on change of
! data type (I/L/R), level_blocking needs to be 1 on little-endian
! systems when packing is done.
      INTEGER, PARAMETER :: level_blocking = 20
#else
      INTEGER, PARAMETER :: level_blocking = 1
#endif

! Local arrays
      REAL                                                              &
     &  buf((((buflen+um_sector_size)/um_sector_size)                   &
     &  * um_sector_size) * level_blocking)
!DIR$ CACHE_ALIGN buf

      REAL temp_buf(buflen)   ! For 32 bit packing

!--------------------------------------------------------------

      ICODE=0
      CMESSAGE=''

      IF (mype  ==  0) THEN
        WRITE(6,'(/,'' WRITING UNIFIED MODEL DUMP ON UNIT'',I3)')NFTOUT
        WRITE(6,'('' #####################################'',/)')
      ENDIF

! Select the relevant decomposition type for this dump

      orig_decomp=current_decomp_type

      IF (SUBMODEL_ID  ==  A_IM) THEN
        IF (current_decomp_type  /=  decomp_standard_atmos)             &
! DEPENDS ON: change_decomposition
     &  CALL CHANGE_DECOMPOSITION(decomp_standard_atmos,ICODE)

      ELSEIF (SUBMODEL_ID  ==  O_IM) THEN
        IF (current_decomp_type  /=  decomp_standard_ocean)             &
! DEPENDS ON: change_decomposition
     &  CALL CHANGE_DECOMPOSITION(decomp_standard_ocean,ICODE)

      ELSEIF (SUBMODEL_ID  ==  W_IM) THEN
        IF (current_decomp_type  /=  decomp_standard_wave)              &
! DEPENDS ON: change_decomposition
     &  CALL CHANGE_DECOMPOSITION(decomp_standard_wave,ICODE)

      ELSE  ! unsupported decomposition type
        WRITE(6,*)                                                      &
     &    'UM_WRITEDUMP : Could not change to decomposition required ', &
     &    'for submodel type ',SUBMODEL_ID
        ICODE=1
        CMESSAGE='Unsupported submodel for MPP code'
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ICODE,Cmessage)
      ENDIF

      IF (ICODE  /=  0) THEN
        ICODE=2
        WRITE(6,*)                                                      &
     &    'UM_WRITDUMP : Error - Could not set decomposition ',         &
     &    'for selected submodel.'
        CMESSAGE='Unsupported decomposition selected for MPP code'
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ICODE,Cmessage)
      ENDIF

! Reset the disk addresses and lengths for well-formed I/O
! DEPENDS ON: set_dumpfile_address
      CALL SET_DUMPFILE_ADDRESS(FIXHD, LEN_FIXHD,                       &
     &                          LOOKUP, LEN1_LOOKUP,                    &
     &                          LEN2_LOOKUP,                            &
     &                          number_of_data_words_in_memory,         &
     &                          number_of_data_words_on_disk,           &
     &                          disk_address)

! Write out the header records
! DEPENDS ON: writhead
      CALL WRITHEAD(NFTOUT,FIXHD,LEN_FIXHD,                             &
     &              INTHD,LEN_INTHD,                                    &
     &              REALHD,LEN_REALHD,                                  &
     &              LEVDEPC,LEN1_LEVDEPC,LEN2_LEVDEPC,                  &
     &              ROWDEPC,LEN1_ROWDEPC,LEN2_ROWDEPC,                  &
     &              COLDEPC,LEN1_COLDEPC,LEN2_COLDEPC,                  &
     &              FLDDEPC,LEN1_FLDDEPC,LEN2_FLDDEPC,                  &
     &              EXTCNST,LEN_EXTCNST,                                &
     &              DUMPHIST,LEN_DUMPHIST,                              &
     &              CFI1,LEN_CFI1,                                      &
     &              CFI2,LEN_CFI2,                                      &
     &              CFI3,LEN_CFI3,                                      &
     &              LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA,            &
#include "argppx.h"
     &              START_BLOCK,ICODE,CMESSAGE)

      IF (ICODE  /=  0) THEN
        WRITE(6,*) 'UM_WRITDUMP : Error writing dump header ',          &
     &             'on unit ',NFTOUT
        WRITE(6,*) 'Return code from WRITHEAD was ',ICODE,              &
     &             ' and error message was ',CMESSAGE
        ICODE=3
        CMESSAGE='Error writing dump header'
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ICODE,Cmessage)
      ENDIF


      IF (FIXHD(160)  >   0) THEN ! If there is data to write

! Loop over fields and write from D1

        number_of_fields= FIXHD(152)
        address         = 1
        object_index    = 1
        level           = 1
        buffered_fields = 0
        buffer_start    = 0
        buffer_pos      = 1

        DO K=1,number_of_fields  ! loop over fields to write out

          MPP_LOOKUP(P_LBLREC,K)=0
          MPP_LOOKUP(P_NADDR,K)=address

          IF (LOOKUP(LBLREC,K)  >   0) THEN  ! if there's data in
                                             ! the field

! Check that DATA_TYPE is valid no: +/-1 to +/-3
            IF (( ABS(LOOKUP(DATA_TYPE,K))  >=  1) .AND.                &
     &          ( ABS(LOOKUP(DATA_TYPE,K))  <=  3)) THEN

! OK - we're gonna write this one out. Setup buffering parameters
              buffered_fields = buffered_fields + 1

              If (buffer_start == 0) Then
                buffer_start = k
              End If

! Check that the "buf" array is big enough for this field

              IF (LOOKUP(LBLREC,K)  >   BUFLEN) THEN
                WRITE(6,*)                                              &
     &            'UM_WRITDUMP : Field length longer than buffer'
                WRITE(6,*) 'Field :',K
                WRITE(6,*) 'Field length on disk ',LOOKUP(LBLREC,K)
                WRITE(6,*) 'Buffer length ',BUFLEN
                ICODE=4
                CMESSAGE='Field length longer than buffer'
! DEPENDS ON: ereport
                CALL Ereport(RoutineName,ICODE,Cmessage)
              ENDIF

! Set "packed_field" to .TRUE. if 32bit packing has been used
              packed_field=(MOD((LOOKUP(LBPACK,K)),10) == 2)

! Set up the location of the field on disk and how much data
! needs to be written
              field_start=LOOKUP(LBEGIN,buffer_start)

! data_size contains the number of words to data used to store
! the field on disk
              IF (packed_field) THEN
                data_size=(LOOKUP(LBLREC,K)+1)/2
              ELSE
                data_size=LOOKUP(LBLREC,K)
              ENDIF

! data_write_size contains the number of words to data that need to
! be written out for a field. Each field has extra words of dummy data
! added at the end to ensure each field starts on a disk sector
! boundary.
              data_write_size=LOOKUP(LBNREC,K)

! All ocean diagnostics have no wrap-around points
              IF (SUBMODEL_ID  ==  O_IM) THEN
                IF (D1_ADDR(d1_object_type,object_index)  ==            &
     &              diagnostic) THEN
! DEPENDS ON: change_decomposition
                  CALL CHANGE_DECOMPOSITION(decomp_nowrap_ocean,ICODE)
                ELSE
! DEPENDS ON: change_decomposition
                  CALL CHANGE_DECOMPOSITION(decomp_standard_ocean,      &
     &                                      ICODE)
                ENDIF
              ENDIF


! Gather the field (into tempbuf) onto PE 0 for packing and output
! DEPENDS ON: general_gather_field
              CALL GENERAL_GATHER_FIELD(                                &
     &             D1(address),temp_buf,LOCAL_LEN,LOOKUP(LBLREC,K),1,   &
     &             D1_ADDR(1,object_index),0,                           &
     &             ICODE,CMESSAGE)

              IF (ICODE  /=  0) THEN
                WRITE(6,*) 'WRITEDUMP: Call to GENERAL_GATHER_FIELD ',  &
     &                     'failed'
                WRITE(6,*) 'Return code was ',ICODE
                WRITE(6,*) 'Error message was ',CMESSAGE
                WRITE(6,*) 'Field number ',LOOKUP(ITEM_CODE,K)
                WRITE(6,*) 'Dimensions ',LOOKUP(LBNPT,K),               &
     &                                 ' x ',LOOKUP(LBROW,K)
                WRITE(6,*) 'Grid type ',                                &
     &                      D1_ADDR(d1_grid_type,object_index)
                WRITE(6,*) 'Field was not written out'

                ICODE=300
                CMESSAGE='Failure to gather field'
! DEPENDS ON: ereport
                CALL Ereport(RoutineName,ICODE,Cmessage)
              ENDIF

!  Do the packing of the data

              IF (mype  ==  0) THEN
!             Does this field need to be compressed?
                IF( packed_field ) THEN
                  IF(LOOKUP(DATA_TYPE,K)  ==  1) THEN
! DEPENDS ON: pack21
                    CALL PACK21(LOOKUP(LBLREC,K),temp_buf,              &
                                buf(buffer_pos))
                  ENDIF
                ELSE ! no compression required - just do a copy
                  DO i=1,LOOKUP(LBLREC,k)
                    buf(buffer_pos + i - 1) = temp_buf(i)
                  ENDDO
                ENDIF
              ENDIF ! am I PE 0 ?


! If the field was compressed for writing on disk, we need to compress
! and expand the field in memory. This ensures the same field exists in
! memory that would exist if this dump was read back in.

              IF( packed_field ) THEN
                IF(LOOKUP(DATA_TYPE,K)  ==  1) THEN
! DEPENDS ON: pack21
                  CALL PACK21(LOCAL_LEN,D1(address),temp_buf)
! DEPENDS ON: expand21
                  CALL EXPAND21(LOCAL_LEN,temp_buf,D1(address))
                ENDIF
              ENDIF

! Increment position in buffer
              buffer_pos = buffer_pos + data_write_size

! Now write out the buffered fields if we've reached the limit
              If (buffered_fields == level_blocking .OR.                &
     &            k == number_of_fields) Then

! Set the position in the dump
! DEPENDS ON: setpos
                CALL SETPOS(NFTOUT,field_start,ICODE)
                IF (ICODE  /=  0) THEN
                  WRITE(6,*)                                            &
     &            'UM_WRITDUMP - SETPOS failed to move file pointer ',  &
     &            ' to ',field_start,' on unit ',NFTOUT
                  WRITE(6,*) 'SETPOS returned error code ',ICODE
                  ICODE=6
                  CMESSAGE='SETPOS failed while writing dump.'
! DEPENDS ON: ereport
                  CALL Ereport(RoutineName,ICODE,Cmessage)
                ENDIF

! Write out the data on PE 0
                IF (mype == 0) THEN
#if defined(LITTLE_END)
                  IF( packed_field ) THEN
                    IF(LOOKUP(DATA_TYPE,K)  ==  1) THEN
! Data is packed using CRAY 32 bit method - note that we need to write
! out 2*ISIZE 32 bit words using BUFFO32
                      CALL buffo32_single(NFTOUT,buf,2*(buffer_pos - 1),&
     &                                    LEN_IO,IOSTAT)
! And then halve LEN_IO to satisfy tests against ISIZE
                      LEN_IO = LEN_IO/2
                    ENDIF
                  ELSE
! For non-packed data
                    CALL buffout_single(NFTOUT,buf,buffer_pos - 1,      &
     &                                  LEN_IO,IOSTAT)
                  ENDIF
#else
                  CALL buffout_single(NFTOUT,buf,buffer_pos - 1,        &
     &                                  LEN_IO,IOSTAT)
#endif


                  IF ((IOSTAT  /=  -1.0) .OR.                           &
     &               (LEN_IO  /=  (buffer_pos - 1 ))) THEN
                    WRITE(6,*) 'WRITEDUMP: Error in call to ',          &
     &                         'BUFFOUT_SINGLE'
                    WRITE(6,*) 'Field : ', k
                    WRITE(6,*) 'LEN_IO : ',LEN_IO
                    WRITE(6,*) 'IOSTAT : ',IOSTAT

                    ICODE=400
                    CMESSAGE='Failure writing out field'
! DEPENDS ON: ereport
                    CALL Ereport(RoutineName,ICODE,Cmessage)
                  ENDIF
                ENDIF  ! mype == 0

! Reset data for indexing in buffer
                buffer_pos = 1
                buffer_start = 0
                buffered_fields = 0
              END IF ! Writing out data


              MPP_LOOKUP(P_LBLREC,K)=local_len
              address=address+local_len

              IF (ICODE  /=  0) THEN
                WRITE(6,*)                                              &
     &            'UM_WRITDUMP - Error while attempting to ',           &
     &            'write field ',K,' of ',number_of_fields,             &
     &            ' from unit ',NFTOUT
                WRITE(6,*) 'Return code from WRITE_MULTI was ',ICODE,   &
     &            'and error message was ',CMESSAGE
                WRITE(6,*) 'Field Information: '
                WRITE(6,*) 'Section ',D1_ADDR(d1_section,object_index), &
     &                     ' Item ',D1_ADDR(d1_item,object_index)
                WRITE(6,*) 'Disk address : ',field_start
                WRITE(6,*) 'D1 address : ',address
                WRITE(6,*) 'Number of words requested : ',              &
     &                     data_write_size
                WRITE(6,*) 'Number of words returned : ',len_io

                ICODE=7
                CMESSAGE='Error writing field to dump'
! DEPENDS ON: ereport
                CALL Ereport(RoutineName,ICODE,Cmessage)

              ENDIF ! If an error was detected reading the field

            ELSE ! invalid data type

              IF (( FIXHD(5)  <   6) .OR.                               &
                                           ! Not AC, Var or Cx
     &            ( FIXHD(5)  >   8)) THEN
! DEPENDS ON: pr_look
                CALL PR_LOOK(                                           &
#include "argppx.h"
     &                  LOOKUP,LOOKUP,LEN1_LOOKUP,K)
              ENDIF

              WRITE(6,*) 'UM_WRITDUMP : Failure for field ',K,' of ',   &
     &                   number_of_fields
              WRITE(6,*) 'LOOKUP(DATA_TYPE,K)= ',LOOKUP(DATA_TYPE,K)
              ICODE=8
              CMESSAGE='Invalid data type ( LOOKUP(DATA_TYPE,K) )'
! DEPENDS ON: ereport
              CALL Ereport(RoutineName,ICODE,Cmessage)

            ENDIF ! Check the LOOKUP(DATA_TYPE,K) is valid

#if defined (UTILIO)
          IF(PrintStatus >= PrStatus_Diag) THEN
            IF ((FIXHD(5)  <   6) .OR.                                  &
                                        ! Not AC/Var
     &          (FIXHD(5)  >   8)) THEN !     Obs/Cx

! Print out header and summary of data field
! DEPENDS ON: pr_look
              CALL PR_LOOK(                                             &
#include "argppx.h"
     &          LOOKUP,LOOKUP,LEN1_LOOKUP,K)
              IF (FIXHD(5)  /= 5 ) THEN   !  Skip if boundary dataset
                IF (LOOKUP(DATA_TYPE,K) == 1) THEN  !  Real
! DEPENDS ON: pr_rfld
                  CALL PR_RFLD(LOOKUP,LOOKUP,D1(LOOKUP(NADDR,K)),K)
                ELSE IF(LOOKUP(DATA_TYPE,K) == 2) THEN  !  Integer
! DEPENDS ON: pr_ifld
                  CALL PR_IFLD(LOOKUP,LOOKUP,D1(LOOKUP(NADDR,K)),K)
                ELSE IF(LOOKUP(DATA_TYPE,K) == 3) THEN  !  Logical
! DEPENDS ON: pr_lfld
                  CALL PR_LFLD(LOOKUP,LOOKUP,LEN1_LOOKUP,               &
     &                 D1(LOOKUP(NADDR,K)),K)
                ENDIF ! Test for different data types
              ENDIF ! test for boundary dataset
            ENDIF ! test for field type
          ENDIF ! PrintStatus >= PrStatus_Diag
#endif
          ENDIF ! If there was data in the field

          level=level+1
          IF (level  >   D1_ADDR(d1_no_levels,object_index)) THEN
            level=1
            object_index=object_index+1
          ENDIF

        ENDDO ! K : loop over fields to write out

        IF (PrintStatus  >=  PrStatus_Normal) THEN
          IF (mype  ==  0) THEN
            WRITE(6,*) 'Data successfully written'
            WRITE(6,*) FIXHD(161),' words written to unit ',NFTOUT
            IF ((FIXHD(5)  >=  6) .and.                                 &
                                         ! AC/Var
     &          (FIXHD(5)  <=  8)) THEN ! Obs/Cx
              WRITE(6,*) '(Observational data)'
            ELSE
              WRITE(6,*) '(Model data)'
            ENDIF
          ENDIF ! IF (mype  ==  0)
        ENDIF ! IF (PrintStatus  >=  PrStatus_Normal)

      ENDIF ! IF (FIXHD(160)  >   0)

! Reset to original decomposition type
! DEPENDS ON: change_decomposition
      CALL CHANGE_DECOMPOSITION(orig_decomp,ICODE)

      RETURN
      END SUBROUTINE UM_WRITDUMP
#endif
