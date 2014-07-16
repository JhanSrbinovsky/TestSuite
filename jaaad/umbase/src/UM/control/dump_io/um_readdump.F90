#if defined(C80_1A) || defined(MAKEBC)
#if !defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE UM_READDUMP---------------------------------------
!LL
!LL  Purpose: Reads in model dump on unit NFTIN and checks model
!LL           and dump dimensions for consistency.
!LL
!LL Rewritten and tidied up for vn5.x
!LL
!LL Model   Date      Modification history from model version 5.0
!LL version
!LL 5.0     3/6/99    Rewritten UM_READDUMP for vn5.x       P.Burton
!LL 5.2     18/12/00  Add levels argument to read_multi  P.Burton
!LL 5.3     30/01/02  Remove MPP defs. Pass size of unpacked field
!LL                   to read_multi to match change in READFLDS.
!LL                   S.D.Mullerworth
!    5.5    02/08/00  Modification for parallelisation of WAM.
!                     Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
!    6.0    02/07/03  Remove MPP defs completely.   E.Leung
!    6.0    14/11/03  Fix to avoid problem with reading WF block
!                     into obs array - only select valid data,
!                     discarding packing area. M. Huddleston
!   6.1   22/10/04  Correct logic when printing file type. R. Hill
!    6.2    23/11/05  Removed all references to the wavemodel.
!                     T.Edwards
!    6.2    26/05/06  Correct External statements. P.Selwood
!    6.1    22/10/04  Correct logic when printing file type. R. Hill
!    6.2    06/12/05  Removed DIAG80 CPP define. T.Edwards
!LL
!LL
!LL

! Subroutine Interface

      SUBROUTINE UM_READDUMP(NFTIN,FIXHD,LEN_FIXHD                      &
     & ,INTHD,LEN_INTHD                                                 &
     & ,REALHD,LEN_REALHD                                               &
     & ,LEVDEPC,LEN1_LEVDEPC,LEN2_LEVDEPC                               &
     & ,ROWDEPC,LEN1_ROWDEPC,LEN2_ROWDEPC                               &
     & ,COLDEPC,LEN1_COLDEPC,LEN2_COLDEPC                               &
     & ,FLDDEPC,LEN1_FLDDEPC,LEN2_FLDDEPC                               &
     & ,EXTCNST,LEN_EXTCNST                                             &
     & ,DUMPHIST,LEN_DUMPHIST                                           &
     & ,CFI1,LEN_CFI1                                                   &
     & ,CFI2,LEN_CFI2                                                   &
     & ,CFI3,LEN_CFI3                                                   &
     & ,LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP                                  &
     & ,MPP_LOOKUP,MPP_LEN1_LOOKUP                                      &
     & ,SUBMODEL_ID,N_OBJS_D1,D1_ADDR                                   &
     &      ,LEN_DATA,D1,                                               &
#include "argppx.h"
     &  READ_HEADER                                                     &
     &  )

      IMPLICIT NONE

       INTEGER                                                          &
     &  NFTIN                                                           &
                        !IN Unit number of dump
     &, LEN_FIXHD                                                       &
                        !IN Length of fixed length header
     &, LEN_INTHD                                                       &
                        !IN Length of integer header
     &, LEN_REALHD                                                      &
                        !IN Length of real header
     &, LEN1_LEVDEPC                                                    &
                        !IN 1st dim of level dep consts
     &, LEN2_LEVDEPC                                                    &
                        !IN 2nd dim of level dep consts
     &, LEN1_ROWDEPC                                                    &
                        !IN 1st dim of row dep consts
     &, LEN2_ROWDEPC                                                    &
                        !IN 2nd dim of row dep consts
     &, LEN1_COLDEPC                                                    &
                        !IN 1st dim of column dep consts
     &, LEN2_COLDEPC                                                    &
                        !IN 2nd dim of column dep consts
     &, LEN1_FLDDEPC                                                    &
                        !IN 1st dim of field dep consts
     &, LEN2_FLDDEPC                                                    &
                        !IN 2nd dim of field dep consts
     &, LEN_EXTCNST                                                     &
                        !IN Length of extra constants
     &, LEN_DUMPHIST                                                    &
                        !IN Length of history block
     &, LEN_CFI1                                                        &
                        !IN Length of comp field index 1
     &, LEN_CFI2                                                        &
                        !IN Length of comp field index 2
     &, LEN_CFI3                                                        &
                        !IN Length of comp field index 3
     &, LEN1_LOOKUP                                                     &
                        !IN 1st dim of lookup
     &, LEN2_LOOKUP                                                     &
                        !IN 2nd dim of lookup
     &, MPP_LEN1_LOOKUP                                                 &
                        !IN 1st dim of MPP lookup
     &, SUBMODEL_ID                                                     &
                        !IN submodel of dump
     &, N_OBJS_D1                                                       &
                        !IN number of objects (3D fields) in D1
     &, LEN_DATA        !IN length of model data

      INTEGER                                                           &
     &  FIXHD(LEN_FIXHD)                                                &
                           !IN Fixed length header
     &, INTHD(LEN_INTHD)                                                &
                           !IN Integer header
     &, LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)                                 &
                           !IN PP lookup tables
     &, CFI1(LEN_CFI1+1)                                                &
                           !IN Compressed field index no 1
     &, CFI2(LEN_CFI2+1)                                                &
                           !IN Compressed field index no 2
     &, CFI3(LEN_CFI3+1)                                                &
                           !IN Compressed field index no 3

     &, MPP_LOOKUP(MPP_LEN1_LOOKUP,LEN2_LOOKUP)
                           !OUT Local processor lookup

      REAL                                                              &
     &  REALHD(LEN_REALHD)                                              &
                           !IN Real header
     &, LEVDEPC(1+LEN1_LEVDEPC*LEN2_LEVDEPC)                            &
                                             !IN Lev dep consts
     &, ROWDEPC(1+LEN1_ROWDEPC*LEN2_ROWDEPC)                            &
                                             !IN Row dep consts
     &, COLDEPC(1+LEN1_COLDEPC*LEN2_COLDEPC)                            &
                                             !IN Col dep consts
     &, FLDDEPC(1+LEN1_FLDDEPC*LEN2_FLDDEPC)                            &
                                             !IN Field dep consts
     &, EXTCNST(LEN_EXTCNST+1)                                          &
                                 !IN Extra constants
     &, DUMPHIST(LEN_DUMPHIST+1)                                        &
                                 !IN History block

     &, D1(LEN_DATA)       !OUT Local subdomain of dump

      LOGICAL                                                           &
     & READ_HEADER         !IN  True if header is to be read in

! Parameters required for dimensioning the D1_ADDR array
#include "d1_addr.h"

      INTEGER                                                           &
     &  D1_ADDR(D1_LIST_LEN,N_OBJS_D1)
                           ! IN D1 addressing info.

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
     &, d1_item_code                                                    &
                      ! sec/item in d1_addr converted into single code
     &, number_of_fields                                                &
                      ! total number of fields to read in
     &, field_start                                                     &
                      ! start address of a field in the file
     &, data_size                                                       &
                      ! number of words of data on disk for a field
     &, data_read_size                                                  &
                      ! total number of words to read for a field
     &, data_full_size                                                  &
                      ! total number of words after any unpacking
     &, len_io                                                          &
                      ! number of words of data successfully read
     &, K                                                               &
                      ! loop counter over fields
     &, orig_decomp                                                     &
                      ! current decomposition on entry
     &, address                                                         &
                      ! start address of field in D1
     &, local_len     ! number of words of data put into D1 on this
                      ! processor


      LOGICAL                                                           &
     &  packed_field  ! TRUE if a field has been packed to 32 bits

! Error reporting
      INTEGER       ICODE       ! =0 normal exit; >0 error exit
      CHARACTER*256 Cmessage    ! Error message
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='UM_READDUMP')



!--------------------------------------------------------------

      ICODE=0
      CMESSAGE=''

      IF (mype  ==  0) THEN
        WRITE(6,'(/,'' READING UNIFIED MODEL DUMP ON UNIT'',I3)')NFTIN
        WRITE(6,'('' #####################################'',/)')
      ENDIF

! Change to the relevant decomposition type for this dump

      orig_decomp=current_decomp_type

      IF (SUBMODEL_ID  ==  A_IM) THEN
        IF (current_decomp_type  /=  decomp_standard_atmos)             &
! DEPENDS ON: change_decomposition
     &  CALL CHANGE_DECOMPOSITION(decomp_standard_atmos,ICODE)

      ELSEIF (SUBMODEL_ID  ==  O_IM) THEN
        IF (current_decomp_type  /=  decomp_standard_ocean)             &
! DEPENDS ON: change_decomposition
     &  CALL CHANGE_DECOMPOSITION(decomp_standard_ocean,ICODE)

      ELSEIF (SUBMODEL_ID  ==  X_IM) THEN
        IF (current_decomp_type  /=  decomp_smexe)                      &
! DEPENDS ON: change_decomposition
     &  CALL CHANGE_DECOMPOSITION(decomp_smexe,ICODE)
      ELSE  ! unsupported decomposition type
        WRITE(6,*)                                                      &
     &    'UM_READDUMP : Could not change to decomposition required ',  &
     &    'for submodel type ',SUBMODEL_ID
        ICODE=1
        CMESSAGE='Unsupported submodel for MPP code'
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ICODE,Cmessage)
      ENDIF

      IF (ICODE  /=  0) THEN
        ICODE=2
        WRITE(6,*) 'UM_READDUMP : Error - Could not set decomposition ',&
     &             'for selected submodel.'
        CMESSAGE='Unsupported decomposition selected for MPP code'
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ICODE,Cmessage)
      ENDIF

! Read in the header records, and do a consistency check

      IF (READ_HEADER) THEN

! DEPENDS ON: readhead
        CALL READHEAD(NFTIN,FIXHD,LEN_FIXHD,                            &
     &                INTHD,LEN_INTHD,                                  &
     &                REALHD,LEN_REALHD,                                &
     &                LEVDEPC,LEN1_LEVDEPC,LEN2_LEVDEPC,                &
     &                ROWDEPC,LEN1_ROWDEPC,LEN2_ROWDEPC,                &
     &                COLDEPC,LEN1_COLDEPC,LEN2_COLDEPC,                &
     &                FLDDEPC,LEN1_FLDDEPC,LEN2_FLDDEPC,                &
     &                EXTCNST,LEN_EXTCNST,                              &
     &                DUMPHIST,LEN_DUMPHIST,                            &
     &                CFI1,LEN_CFI1,                                    &
     &                CFI2,LEN_CFI2,                                    &
     &                CFI3,LEN_CFI3,                                    &
     &                LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP,                   &
     &                LEN_DATA,                                         &
#include "argppx.h"
     &                START_BLOCK,ICODE,CMESSAGE)

        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'UM_READDUMP : Error reading dump header ',        &
     &               'on unit ',NFTIN
          WRITE(6,*) 'Return code from READHEAD was ',ICODE,            &
     &               ' and error message was ',CMESSAGE
          ICODE=3
          CMESSAGE='Error reading dump header'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF

      ENDIF  ! IF (READ_HEADER)


      IF (FIXHD(160)  >   0) THEN ! If there is data to read

! Loop over fields and read into D1

        number_of_fields=FIXHD(152)
        address=1
        object_index=1
        level=1

        DO K=1,number_of_fields  ! loop over fields to read in

          MPP_LOOKUP(P_LBLREC,K)=0
          MPP_LOOKUP(P_NADDR,K)=address

          IF (LOOKUP(LBLREC,K)  >   0) THEN ! If there's data in
                                            ! the field

! Check that DATA_TYPE is valid no: +/-1 to +/-3
            IF (( ABS(LOOKUP(DATA_TYPE,K))  >=  1) .AND.                &
     &          ( ABS(LOOKUP(DATA_TYPE,K))  <=  3)) THEN

! Set "packed_field" to .TRUE. if 32bit packing has been used
              packed_field=(MOD((LOOKUP(LBPACK,K)),10) == 2)

! Check that the diagnostic in the dump matches that expected
! from D1_ADDR

              IF (D1_ADDR(d1_object_type,object_index)  ==  diagnostic) &
     &          THEN

                d1_item_code= (D1_ADDR(d1_section,object_index)*1000) + &
     &                         D1_ADDR(d1_item,object_index)
                IF (LOOKUP(ITEM_CODE,K)  /=  d1_item_code) THEN
                  WRITE(6,*)                                            &
     &              'UM_READDUMP : Dump field ',K,                      &
     &              ' does not match STASH request for item ',          &
     &              D1_ADDR(d1_item,object_index),                      &
     &              ' section ',D1_ADDR(d1_section,object_index)
                  WRITE(6,*) 'Expected code ',LOOKUP(ITEM_CODE,K)
                  CMESSAGE='UM_READDUMP Dump does not match STASH list'
                  ICODE=4
! DEPENDS ON: ereport
                  CALL Ereport(RoutineName,ICODE,Cmessage)
                ENDIF ! IF (LOOKUP(ITEM_CODE,K)  /=  d1_item_code)
              ENDIF ! IF (D1_ADDR(d1_object_type,object_index)  == 
                    !     diagnostic)

! Set up the location of the field on disk and how much data
! needs to be read in
              field_start=LOOKUP(LBEGIN,K) ! position of field in file

! data_size contains the number of words to data used to store
! the field on disk
              IF (packed_field) THEN
                data_size=(LOOKUP(LBLREC,K)+1)/2
              ELSE
                data_size=LOOKUP(LBLREC,K)
              ENDIF

! data_read_size contains the number of words to data that need to
! be read in for a field. Each field has extra words of dummy data
! added at the end to ensure each field starts on a disk sector
! boundary. The last field on a dump does not have these extra words
! added
              IF (K  /=  number_of_fields) THEN
                data_read_size=LOOKUP(LBNREC,K)
              ELSE
                data_read_size=data_size
              ENDIF

! This is the max of number of words required to store the field in
! memory after any unpacking is done, and number of words required
! to read in the data.

              data_full_size=max(LOOKUP(LBLREC,K),data_read_size)

! Move file pointer to the start of the field
! DEPENDS ON: setpos
              CALL SETPOS(NFTIN,field_start,ICODE)
              IF (ICODE  /=  0) THEN
                WRITE(6,*)                                              &
     &          'UM_READDUMP - SETPOS failed to move file pointer to ', &
     &          field_start,' on unit ',NFTIN
                WRITE(6,*) 'SETPOS returned error code ',ICODE
                ICODE=5
                CMESSAGE='SETPOS failed while reading dump. See output.'
! DEPENDS ON: ereport
                CALL Ereport(RoutineName,ICODE,Cmessage)
              ENDIF

! DEPENDS ON: read_multi
              CALL READ_MULTI(NFTIN,D1(address),data_read_size,         &
     &                        data_full_size,                           &
     &                        len_io,local_len,LOOKUP(1,K),FIXHD(12),   &
     &                        D1_ADDR(1,object_index),1,                &
     &                        ICODE,CMESSAGE)

              MPP_LOOKUP(P_LBLREC,K)=local_len
              address=address+local_len

              IF (ICODE  /=  0) THEN
                WRITE(6,*)                                              &
     &            'UM_READDUMP - Error while attempting to read field ',&
     &            K,' of ',number_of_fields,' from unit ',              &
     &            NFTIN
                WRITE(6,*) 'Return code from READ_MULTI was ',ICODE,    &
     &            'and error message was ',CMESSAGE
                WRITE(6,*) 'Field Information: '
                WRITE(6,*) 'Section ',D1_ADDR(d1_section,object_index), &
     &                     ' Item ',D1_ADDR(d1_item,object_index)
                WRITE(6,*) 'Disk address : ',field_start
                WRITE(6,*) 'D1 address : ',address
                WRITE(6,*) 'Number of words requested : ',data_read_size
                WRITE(6,*) 'Number of words returned : ',len_io

                ICODE=6
                CMESSAGE='Error reading field from dump'
! DEPENDS ON: ereport
                CALL Ereport(RoutineName,ICODE,Cmessage)

              ENDIF ! If an error was detected reading the field

            ELSE ! Error in LOOKUP(DATA_TYPE,K)

              IF (( FIXHD(5)  <   6) .OR.                               &
                                           ! Not AC, Var or Cx
     &            ( FIXHD(5)  >   8)) THEN
! DEPENDS ON: pr_look
                CALL PR_LOOK(                                           &
#include "argppx.h"
     &                  LOOKUP,LOOKUP,LEN1_LOOKUP,K)
              ENDIF

              WRITE(6,*) 'UM_READDUMP : Failure for field ',K,' of ',   &
     &                   number_of_fields
              WRITE(6,*) 'LOOKUP(DATA_TYPE,K)= ',LOOKUP(DATA_TYPE,K)
              ICODE=7
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
     &                         D1(LOOKUP(NADDR,K)),K)
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

        ENDDO ! K : loop over fields to read in

        IF (PrintStatus  >=  PrStatus_Normal) THEN
          IF (mype  ==  0) THEN
            WRITE(6,*) 'Data successfully read'
            WRITE(6,*) FIXHD(161),' words read from unit ',NFTIN
            IF ((FIXHD(5)  >=  6) .AND.                                 &
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
      END SUBROUTINE UM_READDUMP
!LL  SUBROUTINE READDUMP---------------------------------------
!LL
!LL  Purpose: Reads in model obs dump on unit NFTIN and checks model
!LL           and dump dimensions for consistency.
!LL
!LL  Code mostly copied from original READDUMP
!LL
!LL  Model            Modification history from model version 4.3:
!LL version  Date
!LL   4.3  19/3/97   New deck introduced                    P.Burton
!     6.1  27/07/04  Correction to avoid problems when obs file only
!                    has a small number of obs. A. Hines.
!     6.2  01/06/06  Replace BUFFIN_ACOBS and BUFFIN_SHMEM with
!                    BUFFIN. P.Selwood
!LL
!LL  Programming standard: Unified Model Documentation Paper No 3
!LL                        Version No 1 15/1/90
!LL
!LL  Logical component: R30
!LL
!LL  System task: F3
!LL
!LL  Documentation: Unified Model Documentation Paper No F3
!LL                 Version No 5 9/2/90
!LLEND---------------------------------------------------------
!
!*L Arguments:-------------------------------------------------
#endif
#endif
