#if defined(C80_1A) || defined(UTILIO) || defined(VAROPSVER)
#if !defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Reads in a number of fields from UM format file
!
! Subroutine Interface:
      SUBROUTINE READFLDS(NFTIN, NUMBER_OF_FIELDS,                      &
     &                    FIRST_FIELD,LOOKUP,LEN1_LOOKUP_ARG,           &
     &                    D1,DISUSED,FIXHD,                             &
#include "argppx.h"
#if defined(CONVIEEE) || defined(CUMF) || defined(PUMF)                \
 || defined(VAROPSVER) || defined(MAKEBC)
     &                    EXPAND,                                       &
#endif
     &                    ICODE,CMESSAGE)

      IMPLICIT NONE
!
! Description:
!  Reads in NUMBER_OF_FIELDS fields from file on unit NFTIN,
!  starting at field number FIRST_FIELD. The data is returned
!  in the D1 array.
!
! Current code owner: Paul Burton
!
! Subroutine Arguments:

      INTEGER                                                           &
     &  NFTIN                                                           &
                               ! IN: unit number to read data from
     &, NUMBER_OF_FIELDS                                                &
                               ! IN: number of fields to read in
     &, FIRST_FIELD                                                     &
                               ! IN: first field to read in
     &, LEN1_LOOKUP_ARG                                                 &
                               ! IN: first dimension of LOOKUP table
     &, LOOKUP(LEN1_LOOKUP_ARG,*)                                       &
                                   ! IN: lookup table starting
                                   !     at field 1
     &, DISUSED                                                         &
                               ! IN: Not used since move to MPP
     &, FIXHD(*)                                                        &
                               ! IN: fixed length header
#if defined(CONVIEEE) || defined(CUMF) || defined(PUMF)                \
 || defined(VAROPSVER) || defined(MAKEBC)
     &, EXPAND                                                          &
                         ! IN: (=1 if WGDOS or RLE packed data
                               !      is to be expanded)
#endif
     &, ICODE                  ! OUT: return code

      REAL                                                              &
     &  D1(*)                  ! OUT: array to return the data in

      CHARACTER*80                                                      &
     &  CMESSAGE               ! OUT: Error message if ICODE <> 0

! COMMON blocks and PARAMETERs

#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "c_mdi.h"
#include "parvars.h"
#include "d1_addr.h"
#include "typsize.h"
#include "clookadd.h"
#include "cprintst.h"

! Local variables

      INTEGER                                                           &
     &  K                                                               &
                               ! loop over fields to read in
     &, d1_off                                                          &
                               ! local offset into D1 for this field
     &, pack_code                                                       &
                               ! packing code for field
     &, field_start                                                     &
                               ! location of field on disk
     &, data_size                                                       &
                               ! number of words of data on disk
                               ! (including padding for WFIO)
     &, data_read_size                                                  &
                               ! number of words to read from disk
     &, data_full_size                                                  &
                               ! number of words after any unpacking
     &, len_io                                                          &
                               ! number of words read from disk
     &, num_cray_words                                                  &
     &, num_unpack_values                                               &
     &, len_full_word                                                   &
#if defined(CONVIEEE) || defined(CUMF) || defined(PUMF)                \
 || defined(VAROPSVER)
     &, dimx                                                            &
     &, dimy                                                            &
     &, idum                                                            &
#endif
     &, field_item                                                      &
                               ! Item number of field
     &, field_sect                                                      &
                               ! Section number of field
     &, field_model                                                     &
                               ! Model ID of field
     &, grid_type                                                       &
                               ! grid type code
     &, fld_type                                                        &
                               ! field type (P,U or V)
     &, halo_type                                                       &
                               ! halo type code
     &, i                                                               &
                               ! loop index
     &, local_len              ! size of field section on this PE

      REAL                                                              &
     &  A_IO                   ! Return code from BUFFIN

      INTEGER                                                           &
     &  fake_D1_ADDR(D1_LIST_LEN)    ! Fake D1_ADDR record to
                                     ! pass into read_multi
      INTEGER unset                  ! unset values
      PARAMETER (unset=-1)

! Functions

      INTEGER GET_FLD_TYPE,EXPPXI

!--------------------------------------------------------------

#if !defined(VAROPSVER)
      IF (FIXHD(12) <  403) THEN
        WRITE(6,*) 'READFLDS: file created by UM version ', FIXHD(12)
        ICODE=1
        CMESSAGE='READFLDS: Cannot read files from before vn4.3'
        GOTO 9999
      ENDIF
#endif

      d1_off=0

      DO K=FIRST_FIELD,FIRST_FIELD+NUMBER_OF_FIELDS-1

        pack_code=MOD((LOOKUP(LBPACK,K)),10)

! Set up the location of the field on disk
        field_start=LOOKUP(LBEGIN,K) ! position of field in file
        IF (field_start <= 0) THEN
          WRITE(6,*) 'READFLDS: start address =',field_start
          ICODE=1
          CMESSAGE='READFLDS: start address of field not given'
          GOTO 9999
        ENDIF

! data_size contains the number of words of data used to store the
! field on disk (needs to be halved if 32 bit packing has been used)
        IF (pack_code  ==  2) THEN
          data_size=(LOOKUP(LBLREC,K)+1)/2
        ELSE
          data_size=LOOKUP(LBLREC,K)
        ENDIF

! data_full_size is the number of words required to store the field
! in memory after any unpacking is done.
        ! This is to give buf the correct size in RDUNPCK, as
        ! buf will be the final expanded size of whole field
        ! including extra data
        ! WARNING LBEXT - may be -32768 MISSING VALUE !
        if ((pack_code == 4).and.(lookup(lbext, k) >  0)) then
          data_full_size=max(lookup(lbrow, k)*lookup(lbnpt, k)          &
     &      +lookup(lbext, k) ,lookup(lblrec,k))
        else
          data_full_size=max(lookup(lbrow, k)*lookup(lbnpt, k)          &
     &      ,lookup(lblrec,k))
        endif

        if ((lookup(lbrow,k) <  0).or.                                  &
     &      (lookup(lbnpt,k) <  0)) then
          data_full_size=lookup(lblrec,k)
        endif

! data_read_size contains the number of words to data that need to
! be read in for a field. Each field has extra words of dummy data
! added at the end to ensure each field starts on a disk sector
! boundary. The last field on a dump does not have these extra words
! added
        IF (K  /=  (FIRST_FIELD+NUMBER_OF_FIELDS-1)) THEN
          data_read_size=LOOKUP(LBNREC,K)
        ELSE
          data_read_size=data_size
        ENDIF
        IF (data_read_size <  0) THEN
          WRITE(6,*) 'READFLDS: number of words to read =',             &
     &     data_read_size
          ICODE=1
          CMESSAGE='READFLDS: number of words to read not given'
          GOTO 9999
        ENDIF

! data_full_size needs to be at least as big as data_read_size since
! it is used to dimension the BUF array in READ_MULTI.

        data_full_size = max(data_full_size, data_read_size)


! Move file pointer to the start of the field
! DEPENDS ON: setpos
        CALL SETPOS(NFTIN,field_start,ICODE)
        IF (ICODE  /=  0) THEN
          WRITE(6,*)                                                    &
     &     'READFLDS - SETPOS failed to move file pointer to ',         &
     &     field_start,' on unit ',NFTIN
          WRITE(6,*) 'SETPOS returned error code ',ICODE
          ICODE=100
          CMESSAGE='SETPOS failed in READFLDS. See output.'
          GOTO 9999
        ENDIF

! Get some information about this field

        field_item=MOD(LOOKUP(42,K),1000)
        field_sect=(LOOKUP(42,K)-field_item)/1000
        field_model=LOOKUP(45,K)
        if (FIXHD(5)  >=  6 .and. FIXHD(5)  <=  9) then
! Set grid_type and halo_type for ACobs and VARobs, Cx and CovStats
! (FIXHD(5)=6, 7, 8 and 9, respectively).
          grid_type = 1
          halo_type = halo_type_no_halo
        else
! DEPENDS ON: exppxi
          grid_type=EXPPXI(field_model,field_sect,field_item,           &
     &                          ppx_grid_type,                          &
#include "argppx.h"
     &                          ICODE,CMESSAGE)
! DEPENDS ON: exppxi
          halo_type=EXPPXI(field_model,field_sect,field_item,           &
     &                          ppx_halo_type,                          &
#include "argppx.h"
     &                          ICODE,CMESSAGE)
          IF (ICODE  /=  0) THEN
            WRITE(6,*)                                                  &
     &        'READFLDS - EXPPXI failed to get PPXREF information ',    &
     &        'for field: '
            WRITE(6,*) 'Model ID : ',field_model
            WRITE(6,*) 'Section  : ',field_sect
            WRITE(6,*) 'Item     : ',field_item
            WRITE(6,*) 'Error code was ',ICODE
            WRITE(6,*) 'Error message was ',CMESSAGE
            ICODE=1
            CMESSAGE='READFLDS failed to get PPXREF information'
            GOTO 9999
          ENDIF
        endif
#if defined(UTILIO)
        IF ((GRID_TYPE <  0).OR.(GRID_TYPE >  100)) THEN
!         WRITE(6,*)'READFL1A: CANNOT GET GRID_TYPE INFO'
!         WRITE(6,*)'          FROM STASHMASTER FOR '
!         WRITE(6,*)'          SECTION ',FIELD_SECT,' ITEM ',FIELD_ITEM
!         WRITE(6,*)'          GRID_TYPE SET TO 1 (NORMAL GRID)'
        ENDIF
        IF (FIXHD(12) >= 500) then
          IF ((HALO_TYPE <  1).OR.(HALO_TYPE >  3)) THEN
!           WRITE(6,*)'READFL1A: CANNOT GET HALO_TYPE INFO'
!           WRITE(6,*)'          FROM STASHMASTER FOR '
!           WRITE(6,*)'          SECTION ',FIELD_SECT,
!    &                           ' ITEM ',FIELD_ITEM
!           WRITE(6,*)'          HALO_TYPE SET TO 3 (NO HALO)'
          ENDIF
        ELSE
          HALO_TYPE=halo_type_no_halo
        ENDIF
#if !defined (PUMF) && !defined (CUMF)
        grid_type=1
#endif
        IF (LOOKUP(LBHEM,K) <  99) THEN
          halo_type=3
        ENDIF
#endif
! DEPENDS ON: get_fld_type
        fld_type=GET_FLD_TYPE(grid_type) ! field type P,U or V

! Set up fake D1_ADDR record to describe data to be read in
! Only set those items actually required by read_multi
! Assume that no diagnostic type fields will be read.

        DO i=1,D1_LIST_LEN
          fake_D1_ADDR(i)=unset
        ENDDO

        fake_D1_ADDR(d1_object_type)=prognostic
        fake_D1_ADDR(d1_imodl)=field_model
        fake_D1_ADDR(d1_section)=field_sect
        fake_D1_ADDR(d1_item)=field_item
        fake_D1_ADDR(d1_halo_type)=halo_type
        fake_D1_ADDR(d1_length)=lasize(1,fld_type,halo_type)*           &
     &                          lasize(2,fld_type,halo_type)
        fake_D1_ADDR(d1_no_levels)=1

! Grid type - for LBCs we need some special logic...
        IF ( (LOOKUP(LBHEM,K)  >=  99) .AND.                            &
                                               ! This is a LBC field
     &       (LOOKUP(LBHEM,K)  <   1000)) THEN

          IF (field_model  ==  ATMOS_IM) THEN
            IF (LOOKUP(LBHEM,K)  ==  99) THEN ! Old style LBCs
              fake_d1_ADDR(d1_grid_type)=ppx_atm_rim
            ELSE ! New style LBCs with different field types

              fake_d1_ADDR(d1_grid_type)=grid_type

              IF (grid_type  ==  ppx_atm_lbc_orog) THEN
                fake_D1_ADDR(d1_length)=                                &
     &            LENRIMA(fld_type,halo_type,rima_type_orog)*           &
     &            (LOOKUP(LBHEM,K)-100)
              ELSE
                fake_D1_ADDR(d1_length)=                                &
     &            LENRIMA(fld_type,halo_type,rima_type_norm)*           &
     &            (LOOKUP(LBHEM,K)-100)
              ENDIF ! IF (grid_type  ==  ppx_atm_lbc_orog)
              fake_D1_ADDR(d1_no_levels)=LOOKUP(LBHEM,K)-100
            ENDIF ! IF (LOOKUP(LBHEM,K)  ==  99)

          ELSE IF (field_model  ==  OCEAN_IM) THEN

            IF (LOOKUP(LBHEM,K)  >   99) THEN ! new style LBCs only
              fake_d1_ADDR(d1_grid_type)= lookup(lbcode,k)
              fake_d1_ADDR(d1_no_levels)=LOOKUP(LBHEM,K)-100
            END IF ! LOOKUP(LBHEM,K)

          ELSE IF (field_model  ==  WAVE_IM) THEN
            fake_D1_ADDR(d1_grid_type)=ppx_wam_rim

          ELSE
            ICODE=2
            write(6,*) 'READFLDS: Cannot process LBC for model type ',  &
     &                 field_model
            CMESSAGE='READFLDS : Cannot read LBCS for this model type'
            GOTO 9999
          ENDIF

        ELSE IF (FIXHD(5) == 4  .and.                                   &
                                                 ! Ancillary File
     &          (Mod( LOOKUP(LBPACK,K)/10, 10 ) == 0 )  .and.           &
     &           grid_type == ppx_atm_compressed ) THEN

          ! Compressed in stashmaster but uncompressed in header
          SELECT CASE( fld_type )
            CASE( fld_type_p )
              fake_D1_ADDR(d1_grid_type) = ppx_atm_tall

            CASE( fld_type_u )
              fake_D1_ADDR(d1_grid_type) = ppx_atm_cuall

            CASE( fld_type_v )
              fake_D1_ADDR(d1_grid_type) = ppx_atm_cvall

          END SELECT

        ELSE ! not an LBC
          fake_D1_ADDR(d1_grid_type)=grid_type
        ENDIF

! Read the field from disk and distribute it over the processors

! DEPENDS ON: read_multi
        CALL READ_MULTI(NFTIN,D1(d1_off+1),data_read_size,              &
     &                  data_full_size,len_io,local_len,                &
     &                  LOOKUP(1,K),FIXHD(12),fake_D1_ADDR,             &
     &                  fake_D1_ADDR(d1_no_levels),                     &
#if defined(CONVIEEE) || defined(CUMF) || defined(PUMF)      \
 || defined(VAROPSVER) || defined(MAKEBC)
     &                  expand,                                         &
#endif
     &                  ICODE,CMESSAGE)

! Check for I/O errors

        IF ((ICODE  /=  0) .OR. (len_io  /=  data_read_size)) THEN
          WRITE(6,*) 'READFLDS : Error reading in field ',K
          WRITE(6,*) 'UNIT ',NFTIN
          WRITE(6,*) 'MODEL ID ',field_model
          WRITE(6,*) 'SECTION ',field_sect
          WRITE(6,*) 'ITEM ',field_item
          IF (FIXHD(5) <  6 .OR. FIXHD(5) >  10) THEN

! DEPENDS ON: pr_look
            CALL PR_LOOK(                                               &
#include "argppx.h"
     &               LOOKUP,LOOKUP,LEN1_LOOKUP_ARG,K)

          ENDIF
          CMESSAGE='READFLDS: I/O error'
          GOTO 9999
        ENDIF
#if defined(UTILIO)
       IF(PrintStatus >= PrStatus_Diag) THEN
! Write out information about the field just read in
        IF (mype  ==  0) THEN
          IF (FIXHD(5) <  6 .OR. FIXHD(5) >  10) THEN

! DEPENDS ON: pr_look
            CALL PR_LOOK(                                               &
#include "argppx.h"
     &               LOOKUP,LOOKUP,LEN1_LOOKUP_ARG,K)

#if defined(CONVIEEE) || defined(CUMF) || defined(PUMF)                \
 || defined(VAROPSVER) || defined(MAKEBC)
            IF ( ((pack_code  ==  1).OR.(pack_code  ==  4)) .AND.       &
     &                                  (EXPAND  /=  1)) THEN
#else
            IF (pack_code  ==  1) THEN
#endif
              WRITE(6,*) 'WGDOS packing not supported .',               &
     &                   'Field summary omitted.'
            ELSEIF (pack_code  ==  3) THEN
              WRITE(6,*) 'GRIB compression not supported .',            &
     &                   'Field summary omitted.'
              WRITE(6,*) 'RLE packing not supported .',                 &
     &                   'Field summary omitted.'
            ELSEIF (FIXHD(5)  ==  5) THEN
              WRITE(6,*) 'Boundary dataset .',                          &
     &                   'Field summary omitted.'
            ELSEIF (LOOKUP(DATA_TYPE,K)  <   0) THEN
              WRITE(6,*) 'Time series field. ',                         &
     &                   'Field summary omitted.'
            ELSE
! Write out summary of field

              IF ( (FIXHD(2)  ==  1) .AND.                              &
     &             (LOOKUP(ITEM_CODE,K) == 30) ) THEN ! Land-sea mask

! DEPENDS ON: pr_lfld
                CALL PR_LFLD(LOOKUP,LOOKUP,64,D1(d1_off+1),K)

              ELSE IF(LOOKUP(DATA_TYPE,K) == 1) THEN  !  Real
! DEPENDS ON: pr_rfld
                CALL PR_RFLD(LOOKUP,LOOKUP,D1(d1_off+1),K)

              ELSE IF(LOOKUP(DATA_TYPE,K) == 2) THEN  !  Integer
! DEPENDS ON: pr_ifld
                CALL PR_IFLD(LOOKUP,LOOKUP,D1(d1_off+1),K)

              ELSE IF(LOOKUP(DATA_TYPE,K) == 3) THEN  !  Logical
! DEPENDS ON: pr_lfld
                CALL PR_LFLD(LOOKUP,LOOKUP,64,D1(d1_off+1),K)

              ENDIF ! type of field

            ENDIF ! if this field can have a summary

          ENDIF ! if this file is suitable for summaries

        ENDIF ! IF (mype  ==  0)
       ENDIF ! PrintStatus >= PrStatus_Diag
#endif

        d1_off=d1_off+local_len
        if(local_len == 0)then
          d1_off=d1_off+LOOKUP(LBLREC,K)
        endif

      ENDDO ! K : loop over fields

 9999 CONTINUE

      RETURN
      END SUBROUTINE READFLDS

#endif
#endif
