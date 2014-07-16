#if defined(MERGE)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL
!LL  MAIN PROGRAM FOR MERGE UTILITY---------------------------------
!LL
!LL  Written by D.M. Goddard 07/09/93
!LL
!LL  Reviewed by
!LL
!LL  Modification History:
!LL
!LL   3.4  13/06/94    Fixes bug that causes merge to fail for
!LL                    boundary files if  the second file starts at
!LL                    at an earlier hour but later day/month/year.
!LL                    Author D.M.Goddard
!LL   3.5  24/03/95    Changed OPEN to FILE_OPEN  P.Burton
!     4.0   02/10/95  Submodels project. Open PPXREF file and
!                     user preSTASHmaster file for use in PRLOOK
!                     Author D.M. Goddard
!     4.1  18/06/96    Changes to cope with changes in STASH addressing
!                      Author D.M. Goddard.
!LL  4.4   Oct. 1997 Changed error handling from routine HDPPXRF
!LL                  so only fatal (+ve) errors are handled.
!LL                                             Shaun de Witt
!     4.5  30/07/97    Code modified to output well-formed IO.
!                      Author D.M. Goddard.
!     5.1  11/05/00    Added call to InitPrintStatus. D.P.Matthews
!     5.3  22/11/01    Enable MPP as the only option for
!                      small executables         E.Leung
!     5.4  03/09/02    Arguments are added to DECOMPOSE_SMEXE.
!                      Check for landpt only and LBC fields
!                      - fields not supported         E.Leung
!     5.5  07/04/03    Use INTHD 6,7 instead of LOOKUP 18,19 as GLSIZE
!                                                              E.Leung
!     6.2  30/03/06    Add call to FILE_CLOSE at the end of
!                      subroutine merge, such that the IO buffer is
!                      flushed and the file closed before program
!                      exits.
!                      R. Sempers (frpz)
!     6.2  10/04/05    Removed calls to ABORT.  J. Gill
!
!LL  Purpose:
!LL          This program was primarily written to merge two boundary
!LL          datasets for use with the mesoscale model model in .
!LL          test mode. It has been extended to cope with the merging
!LL          of any two SEQUENTIAL datasets in unified model format,
!LL          provided they are of the same type and resolution.
!LL          A namelist allows the user to specify the point of merging
!LL          and in the case of time series to merge the datasets at the
!LL          point the times overlap.
!LL
!LL          MAIN_MERGE reads in fixed length ,integer and lookup
!LL          headers of UM files to be merged, extracts dimensions
!LL          of each file, sets the dimensions of the merged file
!LL          to that of the first input file and then passes these
!LL          values to subroutine MERGE.
!LL
!LL  Documentation:
!LL          UM Doc Paper F?
!LL
!LL  -----------------------------------------------------------------
!LL  SUBROUTINE MERGE-----------------------------------------------
!LL
!LL  Written by D. Goddard 14/07/93
!LL
!LL  Reviewed by
!LL
!LL  Modification History:
!LL
!LL  Purpose:
!LL   w      This program was primarily written to merge two boundary
!LL          datasets for use with the mesoscale model model in .
!LL          test mode. It has been extended to cope with the merging
!LL          of any two SEQUENTIAL datasets in unified model format,
!LL          provided they are of the same type and resolution.
!LL          A namelist allows the user to specify the point of merging
!LL          and in the case of time series to merge the datasets at the
!LL          point the times overlap.
!LL
!LL          MERGE reads in headers from files on NFTIN1 and NFTIN2,
!LL          comparing values. If differences occur where they are not
!LL          expected then the program aborts. Then the user decides if
!LL          the files are to be merged at a stated point or for time
!LL          series where they overlap temporally. This is done through
!LL          namelist CONTROL. If files are to be merged temporally, the
!LL          lookup table from file 1 is scanned for the first record in
!LL          the lookup table of file 2, when a common record is found
!LL          its number is used to set IDIFF.  Otherwise IDIFF is taken
!LL          from the namelist. The new merged file is then produced
!LL          by taking the first IDIFF records from file 1 then the
!LL          whole of file 2.
!LL
!LL
!LL  Documentation:
!LL          UM Doc Paper F?
!LL
!LL  -----------------------------------------------------------------
!*L  Arguments:-------------------------------------------------------
      SUBROUTINE MERGE(LEN_FIXHD1,LEN_INTHD1,LEN_REALHD1,               &
     &  P_ROWS1,P_ROWS2,P_ROWS3,                                        &
     &  ROW_LENGTH1,ROW_LENGTH2,ROW_LENGTH3,                            &
     &  LEN1_LEVDEPC1,LEN2_LEVDEPC1,LEN1_ROWDEPC1,                      &
     &  LEN2_ROWDEPC1,LEN1_COLDEPC1,LEN2_COLDEPC1,                      &
     &  LEN1_FLDDEPC1,LEN2_FLDDEPC1,LEN_EXTCNST1,                       &
     &  LEN_DUMPHIST1,LEN_CFI11,LEN_CFI21,LEN_CFI31,                    &
     &  LEN1_LOOKUP1,LEN2_LOOKUP1,LEN_DATA1,P_FIELD1,                   &
     &  LEN_FIXHD2,LEN_INTHD2,LEN_REALHD2,                              &
     &  LEN1_LEVDEPC2,LEN2_LEVDEPC2,LEN1_ROWDEPC2,                      &
     &  LEN2_ROWDEPC2,LEN1_COLDEPC2,LEN2_COLDEPC2,                      &
     &  LEN1_FLDDEPC2,LEN2_FLDDEPC2,LEN_EXTCNST2,                       &
     &  LEN_DUMPHIST2,LEN_CFI12,LEN_CFI22,LEN_CFI32,                    &
     &  LEN1_LOOKUP2,LEN2_LOOKUP2,LEN_DATA2,P_FIELD2,                   &
     &  LEN_FIXHD3,LEN_INTHD3,LEN_REALHD3,                              &
     &  LEN1_LEVDEPC3,LEN2_LEVDEPC3,LEN1_ROWDEPC3,                      &
     &  LEN2_ROWDEPC3,LEN1_COLDEPC3,LEN2_COLDEPC3,                      &
     &  LEN1_FLDDEPC3,LEN2_FLDDEPC3,LEN_EXTCNST3,                       &
     &  LEN_DUMPHIST3,LEN_CFI13,LEN_CFI23,LEN_CFI33,                    &
     &  LEN1_LOOKUP3,LEN2_LOOKUP3,LEN_DATA3,P_FIELD3                    &
     & ,NFTIN1,NFTIN2)

      IMPLICIT NONE

      INTEGER                                                           &
     & LEN_FIXHD1                                                       &
                    !IN Length of fixed length header on file 1
     &,LEN_INTHD1                                                       &
                    !IN Length of integer header on file 1
     &,LEN_REALHD1                                                      &
                    !IN Length of real header on file 1
     &,LEN1_LEVDEPC1                                                    &
                    !IN 1st dim of lev dependent consts on file 1
     &,LEN2_LEVDEPC1                                                    &
                    !IN 2nd dim of lev dependent consts on file 1
     &,LEN1_ROWDEPC1                                                    &
                    !IN 1st dim of row dependent consts on file 1
     &,LEN2_ROWDEPC1                                                    &
                    !IN 2nd dim of row dependent consts on file 1
     &,LEN1_COLDEPC1                                                    &
                    !IN 1st dim of col dependent consts on file 1
     &,LEN2_COLDEPC1                                                    &
                    !IN 2nd dim of col dependent consts on file 1
     &,LEN1_FLDDEPC1                                                    &
                    !IN 1st dim of field dependent consts on file 1
     &,LEN2_FLDDEPC1                                                    &
                    !IN 2nd dim of field dependent consts on file 1
     &,LEN_EXTCNST1                                                     &
                    !IN Length of extra consts on file 1
     &,LEN_DUMPHIST1                                                    &
                    !IN Length of history header on file 1
     &,LEN_CFI11                                                        &
                    !IN Length of index1 on file 1
     &,LEN_CFI21                                                        &
                    !IN Length of index2 on file 1
     &,LEN_CFI31                                                        &
                    !IN Length of index3 on file 1
     &,LEN1_LOOKUP1                                                     &
                    !IN 1st dim of LOOKUP on file 1
     &,LEN2_LOOKUP1                                                     &
                    !IN 2nd dim of LOOKUP on file 1
     &,LEN_DATA1                                                        &
                    !IN Length of data on file 1
     &,P_FIELD1                                                         &
                    !IN No of p-points per level on file 1
     &,P_ROWS1                                                          &
     &,ROW_LENGTH1

      INTEGER                                                           &
     & LEN_FIXHD2                                                       &
                    !IN Length of fixed length header on file 2
     &,LEN_INTHD2                                                       &
                    !IN Length of integer header on file 2
     &,LEN_REALHD2                                                      &
                    !IN Length of real header on file 2
     &,LEN1_LEVDEPC2                                                    &
                    !IN 1st dim of lev dependent consts on file 2
     &,LEN2_LEVDEPC2                                                    &
                    !IN 2nd dim of lev dependent consts on file 2
     &,LEN1_ROWDEPC2                                                    &
                    !IN 1st dim of row dependent consts on file 2
     &,LEN2_ROWDEPC2                                                    &
                    !IN 2nd dim of row dependent consts on file 2
     &,LEN1_COLDEPC2                                                    &
                    !IN 1st dim of col dependent consts on file 2
     &,LEN2_COLDEPC2                                                    &
                    !IN 2nd dim of col dependent consts on file 2
     &,LEN1_FLDDEPC2                                                    &
                    !IN 1st dim of field dependent consts on file 2
     &,LEN2_FLDDEPC2                                                    &
                    !IN 2nd dim of field dependent consts on file 2
     &,LEN_EXTCNST2                                                     &
                    !IN Length of extra consts on file 2
     &,LEN_DUMPHIST2                                                    &
                    !IN Length of history header on file 2
     &,LEN_CFI12                                                        &
                    !IN Length of index1 on file 2
     &,LEN_CFI22                                                        &
                    !IN Length of index2 on file 2
     &,LEN_CFI32                                                        &
                    !IN Length of index3 on file 2
     &,LEN1_LOOKUP2                                                     &
                    !IN 1st dim of LOOKUP on file 2
     &,LEN2_LOOKUP2                                                     &
                    !IN 2nd dim of LOOKUP on file 2
     &,LEN_DATA2                                                        &
                    !IN Length of data on file 2
     &,P_FIELD2                                                         &
                    !IN No of p-points per level on file 2
     &,P_ROWS2                                                          &
     &,ROW_LENGTH2

      INTEGER                                                           &
     & LEN_FIXHD3                                                       &
                    ! OUT Length of fixed length header on file 3
     &,LEN_INTHD3                                                       &
                    ! OUT Length of teger header on file 3
     &,LEN_REALHD3                                                      &
                    ! OUT Length of real header on file 3
     &,LEN1_LEVDEPC3                                                    &
                    ! OUT 1st dim of lev dependent consts on file 3
     &,LEN2_LEVDEPC3                                                    &
                    ! OUT 2nd dim of lev dependent consts on file 3
     &,LEN1_ROWDEPC3                                                    &
                    ! OUT 1st dim of row dependent consts on file 3
     &,LEN2_ROWDEPC3                                                    &
                    ! OUT 2nd dim of row dependent consts on file 3
     &,LEN1_COLDEPC3                                                    &
                    ! OUT 1st dim of col dependent consts on file 3
     &,LEN2_COLDEPC3                                                    &
                    ! OUT 2nd dim of col dependent consts on file 3
     &,LEN1_FLDDEPC3                                                    &
                    ! OUT 1st dim of field dependent consts on file 3
     &,LEN2_FLDDEPC3                                                    &
                    ! OUT 2nd dim of field dependent consts on file 3
     &,LEN_EXTCNST3                                                     &
                    ! OUT Length of extra consts on file 3
     &,LEN_DUMPHIST3                                                    &
                    ! OUT Length of history header on file 3
     &,LEN_CFI13                                                        &
                    ! OUT Length of index1 on file 3
     &,LEN_CFI23                                                        &
                    ! OUT Length of index2 on file 3
     &,LEN_CFI33                                                        &
                    ! OUT Length of index3 on file 3
     &,LEN1_LOOKUP3                                                     &
                    ! OUT 1st dim of LOOKUP on file 3
     &,LEN2_LOOKUP3                                                     &
                    ! OUT 2nd dim of LOOKUP on file 3
     &,LEN_DATA3                                                        &
                    ! OUT Length of data on file 3
     &,P_FIELD3                                                         &
                    ! OUT No of p-points per level on file 3
     &,P_ROWS3                                                          &
     &,ROW_LENGTH3

      INTEGER                                                           &
     & NFTIN1                                                           &
                    !IN Unit number for file 1
     &,NFTIN2       !IN Unit number for file 2


! Comdecks: ------------------------------------------------------------
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "clookadd.h"
#include "cntl_io.h"
#include "cstash.h"
#include "decomptp.h"

! Local arrays:---------------------------------------------------------
      INTEGER                                                           &
     & FIXHD1(LEN_FIXHD1),                                              &
                                                 !
     & INTHD1(LEN_INTHD1),                                              &
                                                 !\                    .
     & CFI11(LEN_CFI11+1),CFI21(LEN_CFI21+1),                           &
                                                 ! > file 1 headers
     & CFI31(LEN_CFI31+1),                                              &
                                                 !/
     & LOOKUP1(LEN1_LOOKUP1,LEN2_LOOKUP1)        !

      INTEGER                                                           &
     & FIXHD2(LEN_FIXHD2),                                              &
                                                 !
     & INTHD2(LEN_INTHD2),                                              &
                                                 !\                    .
     & CFI12(LEN_CFI12+1),CFI22(LEN_CFI22+1),                           &
                                                 ! > file 2 headers
     & CFI32(LEN_CFI32+1),                                              &
                                                 !/
     & LOOKUP2(LEN1_LOOKUP2,LEN2_LOOKUP2)        !

      INTEGER                                                           &
     & FIXHD3(256),                                                     &
                                                 !
     & INTHD3(100),                                                     &
                                                 !\                    .
     & CFI13(LEN_CFI13+1),CFI23(LEN_CFI23+1),                           &
                                                 ! > file 3 headers
     & CFI33(LEN_CFI33+1),                                              &
                                                 !/
     & LOOKUP3(LEN1_LOOKUP3,LEN2_LOOKUP3)        !

      REAL                                                              &
     & REALHD1(LEN_REALHD1),                                            &
                                                 !
     & LEVDEPC1(1+LEN1_LEVDEPC1*LEN2_LEVDEPC1),                         &
                                                 !
     & ROWDEPC1(1+LEN1_ROWDEPC1*LEN2_ROWDEPC1),                         &
                                                 !\                    .
     & COLDEPC1(1+LEN1_COLDEPC1*LEN2_COLDEPC1),                         &
                                                 ! > file 1 headers
     & FLDDEPC1(1+LEN1_FLDDEPC1*LEN2_FLDDEPC1),                         &
                                                 !/
     & EXTCNST1(LEN_EXTCNST1+1),                                        &
                                                 !
     & DUMPHIST1(LEN_DUMPHIST1+1),                                      &
                                                 !
     & D1(P_FIELD1)  ! Data array used to read in each field on file 1

      REAL                                                              &
     & REALHD2(LEN_REALHD2),                                            &
                                                 !
     & LEVDEPC2(1+LEN1_LEVDEPC2*LEN2_LEVDEPC2),                         &
                                                 !
     & ROWDEPC2(1+LEN1_ROWDEPC2*LEN2_ROWDEPC2),                         &
                                                 !\                    .
     & COLDEPC2(1+LEN1_COLDEPC2*LEN2_COLDEPC2),                         &
                                                 ! > file 2 headers
     & FLDDEPC2(1+LEN1_FLDDEPC2*LEN2_FLDDEPC2),                         &
                                                 !/
     & EXTCNST2(LEN_EXTCNST2+1),                                        &
                                                 !
     & DUMPHIST2(LEN_DUMPHIST2+1),                                      &
                                                 !
     & D2(P_FIELD2)  ! Data array used to read in each field on file 2

      REAL                                                              &
     & REALHD3(LEN_REALHD3),                                            &
                                                 !
     & LEVDEPC3(1+LEN1_LEVDEPC3*LEN2_LEVDEPC3),                         &
                                                 !
     & ROWDEPC3(1+LEN1_ROWDEPC3*LEN2_ROWDEPC3),                         &
                                                 !\                    .
     & COLDEPC3(1+LEN1_COLDEPC3*LEN2_COLDEPC3),                         &
                                                 ! > file 3 headers
     & FLDDEPC3(1+LEN1_FLDDEPC3*LEN2_FLDDEPC3),                         &
                                                 !/
     & EXTCNST3(LEN_EXTCNST3+1),                                        &
                                                 !
     & DUMPHIST3(LEN_DUMPHIST3+1),                                      &
                                                 !
     & D3(P_FIELD3)  ! Data array used to read in each field on file 3

      INTEGER                                                           &
     & PP_XREF(PPXREF_CODELEN)  !PPXREF codes for a given section/item

! External subroutines called:------------------------------------------
      EXTERNAL ABORT,ABORT_IO,READHEAD,READFLDS,WRITHEAD,WRITFLDS
      EXTERNAL HDPPXRF,GETPPX
!*----------------------------------------------------------------------
!*L  Local variables:---------------------------------------------------
      REAL                                                              &
     & MAX_DIFF  ! Maximum difference between two fields

      INTEGER                                                           &
     & ICODE                                                            &
                    ! Error return code from subroutines
     &,START_BLOCK                                                      &
                    ! READHEAD argument (not used)
     &,I,J,K,L                                                          &
                    ! Loop indices
     &,JMIN                                                             &
                    ! Minimum length of two headers
     &,SECTION                                                          &
                    ! STASH section number
     &,MAX_J                                                            &
                    ! Point number showing max difference in field
     &,IDIFF                                                            &
                    ! Number of records passed through before match
     &,NDIFF                                                            &
                    ! Number of differences between two header records
     &,NRECF1                                                           &
                    ! Number of records to be copied from file 1
     &,LEN_TIMESTEP                                                     &
                    !Combined length of data for a specified time
     &,LEN_BUF                                                          &
     &,MAX_LEN_BUF                                                      &
     &,POS

      INTEGER IROWDEPC1
      INTEGER IROWDEPC2

      INTEGER      disk_address ! Current rounded disk address
      INTEGER      number_of_data_words_on_disk
                                 ! Number of data words on disk
      INTEGER      number_of_data_words_in_memory                       &

     &,NFTOUT                                                           &
                    ! Unit number for file 3
     &,ERROR        ! Return code from subroutine OPEN

! Input arguments for decompose_smexe
      INTEGER                                                           &

     &  global_row_len1,                                                &
                           ! IN  :number of E-W points of entire model
     &  global_n_rows1,                                                 &
                           ! IN  :number of P rows of entire mode
     &  global_row_len2,                                                &
                           ! IN  :number of E-W points of entire model
     &  global_n_rows2,                                                 &
                           ! IN  :number of P rows of entire mode
     &  field_item,                                                     &
     &  field_sect,                                                     &
     &  field_model,                                                    &
     &  grid_type,                                                      &
     &  tot_levels         ! IN  :total number of levels


      CHARACTER                                                         &
     & CMESSAGE*100                                                     &
                    ! Character string returned if ICODE  /=  0
     &,PHRASE*(PPXREF_CHARLEN) ! Name of field
      INTEGER RowNumber

      INTEGER EXPPXI
      INTEGER NFT1,NFT2
      PARAMETER (NFT1=22, NFT2=2)
!*----------------------------------------------------------------------

      NAMELIST/CONTROL/NRECF1

!L 0. Read in PPXREF

      ppxRecs=1
      RowNumber=0
      cmessage = ' '
      ICODE = 0
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_A',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_A'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('MERGE', ICODE, 'Error reading STASHmaster_A')

      END IF
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_O',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_O'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('MERGE', ICODE, 'Error reading STASHmaster_O')

      END IF
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_S',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_S'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('MERGE', ICODE, 'Error reading STASHmaster_S')

      END IF
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_W',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_W'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('MERGE', ICODE, 'Error reading STASHmaster_W')

      ENDIF
      IF(ICODE >  0)THEN
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('MERGE', ICODE, 'Generic error - aborting')

      ENDIF

! DEPENDS ON: getppx
      CALL GETPPX(NFT1,NFT2,'STASHmaster_A',RowNumber,                  &
#include "argppx.h"
     &            ICODE,CMESSAGE)
! DEPENDS ON: getppx
      CALL GETPPX(NFT1,NFT2,'STASHmaster_O',RowNumber,                  &
#include "argppx.h"
     &            ICODE,CMESSAGE)
! DEPENDS ON: getppx
      CALL GETPPX(NFT1,NFT2,'STASHmaster_S',RowNumber,                  &
#include "argppx.h"
     &            ICODE,CMESSAGE)
! DEPENDS ON: getppx
      CALL GETPPX(NFT1,NFT2,'STASHmaster_W',RowNumber,                  &
#include "argppx.h"
     &            ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('MERGE', ICODE, CMESSAGE)

      ENDIF

!User STASHmaster
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(0,' ',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('MERGE', ICODE, CMESSAGE)

      ENDIF

! DEPENDS ON: getppx
      CALL GETPPX(0,NFT2,' ',RowNumber,                                 &
#include "argppx.h"
     &            ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('MERGE', ICODE, CMESSAGE)

      ENDIF
!L 1. Read in file 1 header

      WRITE(6,*)' '
      WRITE(6,*)'          FILE 1'
      WRITE(6,*)'          ------'
! DEPENDS ON: readhead
      CALL READHEAD(NFTIN1,FIXHD1,LEN_FIXHD1,                           &
     &                INTHD1,LEN_INTHD1,                                &
     &                REALHD1,LEN_REALHD1,                              &
     &                LEVDEPC1,LEN1_LEVDEPC1,LEN2_LEVDEPC1,             &
     &                ROWDEPC1,LEN1_ROWDEPC1,LEN2_ROWDEPC1,             &
     &                COLDEPC1,LEN1_COLDEPC1,LEN2_COLDEPC1,             &
     &                FLDDEPC1,LEN1_FLDDEPC1,LEN2_FLDDEPC1,             &
     &                EXTCNST1,LEN_EXTCNST1,                            &
     &                DUMPHIST1,LEN_DUMPHIST1,                          &
     &                CFI11,LEN_CFI11,                                  &
     &                CFI21,LEN_CFI21,                                  &
     &                CFI31,LEN_CFI31,                                  &
     &                LOOKUP1,LEN1_LOOKUP1,LEN2_LOOKUP1,                &
     &                LEN_DATA1,                                        &
#include "argppx.h"
     &                START_BLOCK,ICODE,CMESSAGE)

      IF(ICODE /= 0)THEN
        WRITE(6,*)CMESSAGE,ICODE
! DEPENDS ON: ereport
        CALL EREPORT('MERGE', ICODE, CMESSAGE)

      ENDIF

!L 2. Read in file 2 header

      WRITE(6,*)' '
      WRITE(6,*)'          FILE 2'
      WRITE(6,*)'          ------'
! DEPENDS ON: readhead
      CALL READHEAD(NFTIN2,FIXHD2,LEN_FIXHD2,                           &
     &                INTHD2,LEN_INTHD2,                                &
     &                REALHD2,LEN_REALHD2,                              &
     &                LEVDEPC2,LEN1_LEVDEPC2,LEN2_LEVDEPC2,             &
     &                ROWDEPC2,LEN1_ROWDEPC2,LEN2_ROWDEPC2,             &
     &                COLDEPC2,LEN1_COLDEPC2,LEN2_COLDEPC2,             &
     &                FLDDEPC2,LEN1_FLDDEPC2,LEN2_FLDDEPC2,             &
     &                EXTCNST2,LEN_EXTCNST2,                            &
     &                DUMPHIST2,LEN_DUMPHIST2,                          &
     &                CFI12,LEN_CFI12,                                  &
     &                CFI22,LEN_CFI22,                                  &
     &                CFI32,LEN_CFI32,                                  &
     &                LOOKUP2,LEN1_LOOKUP2,LEN2_LOOKUP2,                &
     &                LEN_DATA2,                                        &
#include "argppx.h"
     &                START_BLOCK,ICODE,CMESSAGE)


      IF(ICODE /= 0)THEN
        WRITE(6,*)CMESSAGE,ICODE
! DEPENDS ON: ereport
        CALL EREPORT('MERGE', ICODE, CMESSAGE)

      ENDIF

!L 3. Compare fixed length headers and substitute the value of
!L    file1 in file3.

      WRITE(6,*)' '
      WRITE(6,*)'FIXED LENGTH HEADER:'
      IF(LEN_FIXHD1 /= LEN_FIXHD2)THEN
        WRITE(6,*)'ERROR: LEN1=',LEN_FIXHD1,' LEN2=',LEN_FIXHD2
        WRITE(6,*)'Files are incompatable and cannot be merged.'
! DEPENDS ON: ereport
        CALL EREPORT('MERGE', 1010,                                     &
     &     'Files are incompatable and cannot be merged.')

      ELSE
        LEN_FIXHD3=LEN_FIXHD1
      ENDIF
      IF(FIXHD1(5) == 3.OR.FIXHD1(5) >  5)THEN
! Data type not supported. Abort.
        WRITE(6,*) 'ERROR Data type not supported'
! DEPENDS ON: ereport
        CALL EREPORT('MERGE', 1011,                                     &
     &     'Data type not supported.')

      ENDIF
      DO I=1,LEN_FIXHD1
        IF(FIXHD1(I) /= FIXHD2(I))THEN
          IF(I >= 2.AND.I <  6)THEN
            WRITE(6,*)'ERROR: FIXHD1(I)=', FIXHD1(I),                   &
     &                       'FIXHD2(I)=', FIXHD2(I)
            WRITE(6,*) 'Files are incompatable and cannot be merged.'
! DEPENDS ON: ereport
            CALL EREPORT('MERGE', 1012,                                 &
     &       'Files are incompatable and cannot be merged')

          ELSE IF(I == 101)THEN
            WRITE(6,*) 'ERROR: integer constant arrays have different ',&
     &                 'lengths'
            WRITE(6,*) 'File 1 = ',FIXHD1(I),' File 2 = ',FIXHD2(I)
! DEPENDS ON: ereport
            CALL EREPORT('MERGE', 1013,                                 &
     &       'integer constant arrays have different lengths')

          ELSE IF(I == 106)THEN
            WRITE(6,*) 'ERROR: real constant arrays have different ',   &
     &                 'lengths'
            WRITE(6,*) 'File 1 = ',FIXHD1(I),' File 2 = ',FIXHD2(I)
! DEPENDS ON: ereport
            CALL EREPORT('MERGE', 1014,                                 &
     &       'real constant arrays have different lengths')

          ELSE IF(I == 111)THEN
            WRITE(6,*) 'ERROR: level dependant constant arrays have ',  &
     &                 ' different lengths'
            WRITE(6,*) 'File 1 = ',FIXHD1(I),' File 2 = ',FIXHD2(I)
! DEPENDS ON: ereport
            CALL EREPORT('MERGE', 1015,                                 &
     &       'level dependant constant arrays have different lengths')

          ELSE IF(I == 116)THEN
            WRITE(6,*) 'ERROR: row dependant constant arrays have ',    &
     &                 ' different lengths'
            WRITE(6,*) 'File 1 = ',FIXHD1(I),' File 2 = ',FIXHD2(I)
! DEPENDS ON: ereport
            CALL EREPORT('MERGE', 1016,                                 &
     &       'row dependant constant arrays have different lengths')

          ELSE IF(I == 121)THEN
            WRITE(6,*) 'ERROR: column dependant constant arrays have ', &
     &                 ' different lengths'
            WRITE(6,*) 'File 1 = ',FIXHD1(I),' File 2 = ',FIXHD2(I)
! DEPENDS ON: ereport
            CALL EREPORT('MERGE', 1017,                                 &
     &       'column dependant constant arrays have different lengths')

          ELSE IF(I == 126)THEN
            WRITE(6,*) 'ERROR: field of constant arrays have ',         &
     &                 ' different lengths'
            WRITE(6,*) 'File 1 = ',FIXHD1(I),' File 2 = ',FIXHD2(I)
! DEPENDS ON: ereport
            CALL EREPORT('MERGE', 1018,                                 &
     &       'field of constant arras have different lengths')

          ELSE IF(I == 127)THEN
            WRITE(6,*) 'ERROR: field of constant arrays have ',         &
     &                 ' different lengths'
            WRITE(6,*) 'File 1 = ',FIXHD1(I),' File 2 = ',FIXHD2(I)
! DEPENDS ON: ereport
            CALL EREPORT('MERGE', 1019,                                 &
     &       'field of constant arrays have different lengths')

          ELSE IF(I == 131)THEN
            WRITE(6,*) 'ERROR: extra consatant arrays have ',           &
     &                 ' different lengths'
            WRITE(6,*) 'File 1 = ',FIXHD1(I),' File 2 = ',FIXHD2(I)
! DEPENDS ON: ereport
            CALL EREPORT('MERGE', 1020,                                 &
     &       'extra constant arrays have different lengths')

          ELSE IF(I == 136)THEN
            WRITE(6,*) 'ERROR: temp historyfile arrays have ',          &
     &                 ' different lengths'
            WRITE(6,*) 'File 1 = ',FIXHD1(I),' File 2 = ',FIXHD2(I)
! DEPENDS ON: ereport
            CALL EREPORT('MERGE', 1021,                                 &
     &       'temp historyfile arrays have different lengths')

          ELSE IF(I == 141)THEN
            WRITE(6,*) 'ERROR: compressed field index 1 arrays have ',  &
     &                 ' different lengths'
            WRITE(6,*) 'File 1 = ',FIXHD1(I),' File 2 = ',FIXHD2(I)
! DEPENDS ON: ereport
            CALL EREPORT('MERGE', 1022,                                 &
     &       'compressed field index 1 arrays have different lengths')

          ELSE IF(I == 143)THEN
            WRITE(6,*) 'ERROR: compressed field index 1 arrays have ',  &
     &                 ' different lengths'
            WRITE(6,*) 'File 1 = ',FIXHD1(I),' File 2 = ',FIXHD2(I)
! DEPENDS ON: ereport
            CALL EREPORT('MERGE', 1023,                                 &
     &       'compressed field index 1 arrays have different lengths')

          ELSE IF(I == 145)THEN
            WRITE(6,*) 'ERROR: compressed field index 1 arrays have ',  &
     &                 ' different lengths'
            WRITE(6,*) 'File 1 = ',FIXHD1(I),' File 2 = ',FIXHD2(I)
! DEPENDS ON: ereport
            CALL EREPORT('MERGE', 1024,                                 &
     &       'compressed field index 1 arrays have different lengths')

          END IF
        END IF
        FIXHD3(I)=FIXHD1(I)
      END DO

!L 4. Compare integer headers and substitute the value of
!L    file1 in file3.

      IF(LEN_INTHD1 >  0.OR.LEN_INTHD2 >  0)THEN
        WRITE(6,*)' '
        WRITE(6,*)'INTEGER HEADER:'
        DO I=1,LEN_INTHD1
          IF(FIXHD1(5) == 5.AND.FIXHD1(12) <= 303)THEN
            INTHD1(15)=INTHD1(13)
          END IF
          IF(FIXHD2(5) == 5.AND.FIXHD2(12) <= 303)THEN
            INTHD2(15)=INTHD2(13)
          END IF
          IF(INTHD1(I) /= INTHD2(I))THEN
            IF(I == 6)THEN
              WRITE(6,*) 'ERROR: Different number of points in row'
              WRITE(6,*) 'File 1 = ',INTHD1(I),' File 2 = ',INTHD2(I)
! DEPENDS ON: ereport
              CALL EREPORT('MERGE', 1025,                               &
     &         'Different number of points in row')

            ELSE IF(I == 7)THEN
              WRITE(6,*) 'ERROR: Different number of points in column'
              WRITE(6,*) 'File 1 = ',INTHD1(I),' File 2 = ',INTHD2(I)
! DEPENDS ON: ereport
              CALL EREPORT('MERGE', 1026,                               &
     &         'Different number of points in column')

            ELSE IF(I == 8)THEN
              WRITE(6,*) 'ERROR: Different number of levels'
              WRITE(6,*) 'File 1 = ',INTHD1(I),' File 2 = ',INTHD2(I)
! DEPENDS ON: ereport
              CALL EREPORT('MERGE', 1027,                               &
     &         'Different number of levels')

            ELSE IF(I == 9)THEN
              IF(FIXHD1(2) == 1)THEN
                IF(FIXHD1(5) == 1.OR.FIXHD1(5) == 2.OR.FIXHD1(5) == 5)  &
     &          THEN
                  WRITE(6,*) 'ERROR: Different number of wet levels'
                  WRITE(6,*) 'File 1 = ',INTHD1(I),                     &
     &                      ' File 2 = ',INTHD2(I)
! DEPENDS ON: ereport
                  CALL EREPORT('MERGE', 1028,                           &
     &             'Different number of wet levels')

                END IF
              END IF
            ELSE IF(I == 10)THEN
              IF(FIXHD1(2) == 1)THEN
                IF(FIXHD1(5) == 1.OR.FIXHD1(5) == 2)THEN
                  WRITE(6,*) 'ERROR: Different number of soil levels'
                  WRITE(6,*) 'File 1 = ',INTHD1(I),                     &
     &                      ' File 2 = ',INTHD2(I)
! DEPENDS ON: ereport
                  CALL EREPORT('MERGE', 1029,                           &
     &             'Different number of soil levels')

                END IF
              END IF
            ELSE IF(I == 12)THEN
              IF(FIXHD1(2) == 1)THEN
                IF(FIXHD1(5) == 1.OR.FIXHD1(5) == 2)THEN
                  WRITE(6,*) 'ERROR: Different number of tracers'
                  WRITE(6,*) 'File 1 = ',INTHD1(I),                     &
     &                      ' File 2 = ',INTHD2(I)
! DEPENDS ON: ereport
                  CALL EREPORT('MERGE', 1030,                           &
     &             'Different number of tracers')

                END IF
              END IF
            ELSE IF(I == 13)THEN
              IF(FIXHD1(2) == 1)THEN
                IF(FIXHD1(5) == 1.OR.FIXHD1(5) == 2)THEN
                  WRITE(6,*) 'ERROR: Different number of boundary ',    &
     &                       'layer  levels'
                  WRITE(6,*) 'File 1 = ',INTHD1(I),                     &
     &                      ' File 2 = ',INTHD2(I)
! DEPENDS ON: ereport
                  CALL EREPORT('MERGE', 1031,                           &
     &             'Different number of boundary layer levels')

                END IF
              END IF
            ELSE IF(I == 15)THEN
              IF(FIXHD1(2) == 1)THEN
                IF(FIXHD1(5) == 5)THEN
                  WRITE(6,*) 'ERROR: Different number of field types'
                  WRITE(6,*) 'File 1 = ',INTHD1(I),                     &
     &                      ' File 2 = ',INTHD2(I)
! DEPENDS ON: ereport
                  CALL EREPORT('MERGE', 1032,                           &
     &             'Different number of field types')

                END IF
              END IF
            ELSE IF(I == 25)THEN
              IF(FIXHD1(2) == 1)THEN
                IF(FIXHD1(5) == 1.OR.FIXHD1(5) == 2)THEN
                  WRITE(6,*) 'ERROR: Different number of land points'
                  WRITE(6,*) 'File 1 = ',INTHD1(I),                     &
     &                      ' File 2 = ',INTHD2(I)
! DEPENDS ON: ereport
                  CALL EREPORT('MERGE', 1033,                           &
     &             'Different number of land points')

                END IF
              END IF
            ELSE IF(I == 26)THEN
              IF(FIXHD1(2) == 1)THEN
                IF(FIXHD1(5) == 1.OR.FIXHD1(5) == 2)THEN
                  WRITE(6,*) 'ERROR: Different number of ozone levels'
                  WRITE(6,*) 'File 1 = ',INTHD1(I),                     &
     &                      ' File 2 = ',INTHD2(I)
! DEPENDS ON: ereport
                  CALL EREPORT('MERGE', 1034,                           &
     &             'Different number of ozone levels')

                END IF
              END IF
            END IF
          END IF
          INTHD3(I)=INTHD1(I)
        END DO
      END IF

!L 5. Compare real headers and substitute the value of
!L    file1 in file3 if elements the same:

      IF(LEN_REALHD1 >  0.OR.LEN_REALHD2 >  0)THEN
        WRITE(6,*)' '
        WRITE(6,*)'REAL HEADER:'
        IF(LEN_REALHD1 /= LEN_REALHD2)THEN
          WRITE(6,*)'WARNING LEN1=',LEN_REALHD1,' LEN2=',LEN_REALHD2
        ELSE
          LEN_REALHD3=LEN_REALHD1
        ENDIF
        DO I=1,LEN_REALHD1
#if defined(T3E)
          IF(XOR(REALHD1(I),REALHD2(I)) /= 0) THEN
#else
          IF(REALHD1(I) /= REALHD2(I))THEN
#endif
            IF(I == 1)THEN
              WRITE(6,*) 'ERROR: Different row spacing'
              WRITE(6,*) 'File 1 = ',REALHD1(I),                        &
     &                  ' File 2 = ',REALHD2(I)
! DEPENDS ON: ereport
              CALL EREPORT('MERGE', 1035,                               &
     &         'Different row spacing')

            ELSE IF(I == 2)THEN
              WRITE(6,*) 'ERROR: Different column spacing'
              WRITE(6,*) 'File 1 = ',REALHD1(I),                        &
     &                  ' File 2 = ',REALHD2(I)
! DEPENDS ON: ereport
              CALL EREPORT('MERGE', 1036,                               &
     &         'Different column spacing')

            ELSE IF(I == 3)THEN
              WRITE(6,*) 'ERROR: Different latitude of 1st row'
              WRITE(6,*) 'File 1 = ',REALHD1(I),                        &
     &                  ' File 2 = ',REALHD2(I)
! DEPENDS ON: ereport
              CALL EREPORT('MERGE', 1037,                               &
     &         'Different lattitude of 1st row')

            ELSE IF(I == 4)THEN
              WRITE(6,*) 'ERROR: Different longitude of 1st row'
              WRITE(6,*) 'File 1 = ',REALHD1(I),                        &
     &                  ' File 2 = ',REALHD2(I)
! DEPENDS ON: ereport
              CALL EREPORT('MERGE', 1038,                               &
     &         'Different longtitude of 1st row')

            ELSE IF(I == 5)THEN
              WRITE(6,*) 'ERROR: Different latitude of pseudo north ',  &
     &                   'pole'
              WRITE(6,*) 'File 1 = ',REALHD1(I),                        &
     &                  ' File 2 = ',REALHD2(I)
! DEPENDS ON: ereport
              CALL EREPORT('MERGE', 1039,                               &
     &         'Different latitude of pseudo north')

            ELSE IF(I == 6)THEN
              WRITE(6,*) 'ERROR: Different longitude of pseudo north ', &
     &                   'pole'
              WRITE(6,*) 'File 1 = ',REALHD1(I),                        &
     &                  ' File 2 = ',REALHD2(I)
! DEPENDS ON: ereport
              CALL EREPORT('MERGE', 1040,                               &
     &         'Different longtitude of pseduo north')

            END IF
          END IF
          REALHD3(I)=REALHD1(I)
        END DO
      END IF

!L 6. Compare level dependent constants

      IF(LEN1_LEVDEPC1 /= LEN1_LEVDEPC2)THEN
        WRITE(6,*)'ERROR different number of levels'
        WRITE(6,*)'LEV1=',LEN1_LEVDEPC1,' LEV2=',LEN1_LEVDEPC2
! DEPENDS ON: ereport
        CALL EREPORT('MERGE', 1041, 'ERROR different number of levels')

      ELSE
        LEN1_LEVDEPC3=LEN1_LEVDEPC1
      ENDIF
      IF(LEN2_LEVDEPC1 >  0.OR.LEN2_LEVDEPC2 >  0)THEN
        WRITE(6,*)' '
        WRITE(6,*)'LEVEL DEPENDENT CONSTS:'
        IF(LEN2_LEVDEPC1 /= LEN2_LEVDEPC2)THEN
          WRITE(6,*)'WARNING LEN1=',LEN2_LEVDEPC1,' LEN2=',LEN2_LEVDEPC2
        ELSE
          LEN2_LEVDEPC3=LEN2_LEVDEPC1
        ENDIF
        DO I=1,LEN2_LEVDEPC1
          DO J=1,LEN1_LEVDEPC1
#if defined(T3E)
            IF(XOR(LEVDEPC1((I-1)*LEN1_LEVDEPC1+J),                     &
     &        LEVDEPC2((I-1)*LEN1_LEVDEPC1+J)) /= 0)THEN
#else
            IF(LEVDEPC1((I-1)*LEN1_LEVDEPC1+J) /=                       &
     &        LEVDEPC2((I-1)*LEN1_LEVDEPC1+J))THEN
#endif
              WRITE(6,*) 'ERROR: Level dependant constants are ',       &
     &                   'different'
              WRITE(6,*) 'Level = ',J,' Item = ',I,                     &
     &                  ' File 1 = ',LEVDEPC1((I-1)*LEN1_LEVDEPC1+J),   &
     &                  ' File 2 = ',LEVDEPC2((I-1)*LEN1_LEVDEPC2+J)
! DEPENDS ON: ereport
              CALL EREPORT('MERGE', 1042,                               &
     &         'Level dependant constants are different')

            END IF
            LEVDEPC3((I-1)*LEN1_LEVDEPC1+J)=                            &
     &      LEVDEPC1((I-1)*LEN1_LEVDEPC1+J)
          END DO
        END DO
      ENDIF

!L 7. Compare row dependent constants

      IF(LEN1_ROWDEPC1 /= LEN1_ROWDEPC2)THEN
        WRITE(6,*)'ERROR different number of rows'
        WRITE(6,*)'ROW1=',LEN1_ROWDEPC1,' ROW2=',LEN1_ROWDEPC2
! DEPENDS ON: ereport
        CALL EREPORT('MERGE', 1043, 'Different number of rows')

      ELSE
        LEN1_ROWDEPC3=LEN1_ROWDEPC1
      ENDIF
      IF(LEN2_ROWDEPC1 >  0.OR.LEN2_ROWDEPC2 >  0)THEN
        WRITE(6,*)' '
        WRITE(6,*)'ROW DEPENDENT CONSTS:'
        IF(LEN2_ROWDEPC1 /= LEN2_ROWDEPC2)THEN
          WRITE(6,*)'WARNING LEN1=',LEN2_ROWDEPC1,' LEN2=',LEN2_ROWDEPC2
        ELSE
          LEN2_ROWDEPC3=LEN2_ROWDEPC1
        ENDIF
! Row dependent constants may be of different data types,
! so comparsion is skipped.
        DO I=1,LEN2_ROWDEPC1
          DO J=1,LEN1_ROWDEPC1
            ROWDEPC3((I-1)*LEN1_ROWDEPC1+J)=                            &
     &      ROWDEPC1((I-1)*LEN1_ROWDEPC1+J)
          END DO
        END DO
      ENDIF
!
!L 8. Compare column dependent constants

      IF(LEN1_COLDEPC1 /= LEN1_COLDEPC2)THEN
        WRITE(6,*)'ERROR different number of columns'
        WRITE(6,*)'COL1=',LEN1_COLDEPC1,' ROW2=',LEN1_ROWDEPC2
! DEPENDS ON: ereport
        CALL EREPORT('MERGE', 1044, 'Different number of columns')

      ELSE
        LEN1_COLDEPC3=LEN1_COLDEPC1
      ENDIF
      IF(LEN2_COLDEPC1 >  0.OR.LEN2_COLDEPC2 >  0)THEN
        WRITE(6,*)' '
        WRITE(6,*)'COLUMN DEPENDENT CONSTS:'
        IF(LEN2_COLDEPC1 /= LEN2_COLDEPC2)THEN
          WRITE(6,*)'WARNING LEN1=',LEN2_COLDEPC1,' LEN2=',LEN2_COLDEPC2
        ELSE
          LEN2_COLDEPC3=LEN2_COLDEPC1
        ENDIF
        DO I=1,LEN2_COLDEPC1
          DO J=1,LEN1_COLDEPC1
#if defined(T3E)
            IF(XOR(COLDEPC1((I-1)*LEN1_COLDEPC1+J),                     &
     &        COLDEPC2((I-1)*LEN1_COLDEPC1+J)) /= 0) THEN
#else
            IF(COLDEPC1((I-1)*LEN1_COLDEPC1+J) /=                       &
     &        COLDEPC2((I-1)*LEN1_COLDEPC1+J))THEN
#endif
              WRITE(6,*) 'ERROR: column dependant constants are ',      &
     &                   'different'
              WRITE(6,*) 'Column = ',J,' Item = ',I,                    &
     &                  ' File 1 = ',COLDEPC1((I-1)*LEN1_COLDEPC1+J),   &
     &                  ' File 2 = ',COLDEPC2((I-1)*LEN1_COLDEPC2+J)
! DEPENDS ON: ereport
              CALL EREPORT('MERGE', 1045,                               &
     &         'Column dependant constants are different')

            END IF
            COLDEPC3((I-1)*LEN1_COLDEPC1+J)=                            &
     &      COLDEPC1((I-1)*LEN1_COLDEPC1+J)
          END DO
        END DO
      ENDIF

!L 9. Compare field dependent constants

      IF(LEN1_FLDDEPC1 /= LEN1_FLDDEPC2)THEN
        WRITE(6,*)'ERROR different number of fields'
        WRITE(6,*)'FLD1=',LEN1_FLDDEPC1,' FLD2=',LEN1_FLDDEPC2
! DEPENDS ON: ereport
        CALL EREPORT('MERGE', 1046, 'different number of fields')

      ELSE
        LEN1_FLDDEPC3=LEN1_FLDDEPC1
      ENDIF
      IF(LEN2_FLDDEPC1 >  0.OR.LEN2_FLDDEPC2 >  0)THEN
        WRITE(6,*)' '
        WRITE(6,*)'FIELD DEPENDENT CONSTS:'
        IF(LEN2_FLDDEPC1 /= LEN2_FLDDEPC2)THEN
          WRITE(6,*)'WARNING LEN1=',LEN2_FLDDEPC1,' LEN2=',LEN2_FLDDEPC2
        ELSE
          LEN2_FLDDEPC3=LEN2_FLDDEPC1
        ENDIF
        DO I=1,LEN2_FLDDEPC1
          DO J=1,LEN1_FLDDEPC1
#if defined(T3E)
            IF(XOR(FLDDEPC1((I-1)*LEN1_FLDDEPC1+J),                     &
     &        FLDDEPC2((I-1)*LEN1_FLDDEPC1+J)) /= 0) THEN
#else
            IF(FLDDEPC1((I-1)*LEN1_FLDDEPC1+J) /=                       &
     &        FLDDEPC2((I-1)*LEN1_FLDDEPC1+J))THEN
#endif
              WRITE(6,*) 'ERROR: field dependant constants are ',       &
     &                   'different'
              WRITE(6,*) 'Field = ',J,' Item = ',I,                     &
     &                  ' File 1 = ',FLDDEPC1((I-1)*LEN1_FLDDEPC1+J),   &
     &                  ' File 2 = ',FLDDEPC2((I-1)*LEN1_FLDDEPC2+J)
! DEPENDS ON: ereport
              CALL EREPORT('MERGE', 1047,                               &
     &         'field dependant constants are different')

            END IF
            FLDDEPC3((I-1)*LEN1_FLDDEPC1+J)=                            &
     &      FLDDEPC1((I-1)*LEN1_FLDDEPC1+J)
          END DO
        END DO
      ENDIF

!L 10. Compare extra constants

      IF(LEN_EXTCNST1 >  0.OR.LEN_EXTCNST2 >  0)THEN
        WRITE(6,*)' '
        WRITE(6,*)'EXTRA CONSTANTS:'
        IF(LEN_EXTCNST1 /= LEN_EXTCNST2)THEN
          WRITE(6,*)'WARNING LEN1=',LEN_EXTCNST1,' LEN2=',LEN_EXTCNST2
        ELSE
          LEN_EXTCNST3=LEN_EXTCNST1
        ENDIF
        DO I=1,LEN_EXTCNST1
#if defined(T3E)
          IF(XOR(EXTCNST1(I),EXTCNST2(I)) /= 0) THEN
#else
          IF(EXTCNST1(I) /= EXTCNST2(I))THEN
#endif
            WRITE(6,*) 'ERROR: extra constants are different'
            WRITE(6,*) 'Item = ',I,                                     &
     &                  ' File 1 = ',EXTCNST1(I),                       &
     &                  ' File 2 = ',EXTCNST2(I)
! DEPENDS ON: ereport
            CALL EREPORT('MERGE', 1048,                                 &
     &       'extra constants are different')

          END IF
          EXTCNST3(I)=EXTCNST1(I)
        END DO
      ENDIF

!L 11. Compare dump history

      IF(LEN_DUMPHIST1 >  0.OR.LEN_DUMPHIST2 >  0)THEN
        WRITE(6,*)' '
        WRITE(6,*)'HISTORY BLOCK:'
        IF(LEN_DUMPHIST1 /= LEN_DUMPHIST2)THEN
          WRITE(6,*)'WARNING LEN1=',LEN_DUMPHIST1,' LEN2=',LEN_DUMPHIST2
        ELSE
          LEN_DUMPHIST3=LEN_DUMPHIST1
        ENDIF
        DO I=1,LEN_DUMPHIST1
#if defined(T3E)
          IF(XOR(DUMPHIST1(I),DUMPHIST2(I)) /= 0) THEN
#else
          IF(DUMPHIST1(I) /= DUMPHIST2(I))THEN
#endif
            WRITE(6,*) 'ERROR: dump histories are different'
            WRITE(6,*) 'Item = ',I,                                     &
     &                  ' File 1 = ',DUMPHIST1(I),                      &
     &                  ' File 2 = ',DUMPHIST2(I)
! DEPENDS ON: ereport
            CALL EREPORT('MERGE', 1049,                                 &
     &       'dump histories are different')

          END IF
          DUMPHIST3(I)=DUMPHIST1(I)
        END DO
      ENDIF

!L 12. Compare compressed index 1

      IF(LEN_CFI11 >  0.OR.LEN_CFI12 >  0)THEN
        WRITE(6,*)' '
        WRITE(6,*)'COMPRESSED INDEX 1:'
        IF(LEN_CFI11 /= LEN_CFI12)THEN
          WRITE(6,*)'WARNING LEN1=',LEN_CFI11,' LEN2=',LEN_CFI12
        ELSE
          LEN_CFI13=LEN_CFI11
        ENDIF
        DO I=1,LEN_CFI11
#if defined(T3E)
          IF(XOR(CFI11(I),CFI12(I)) /= 0) THEN
#else
          IF(CFI11(I) /= CFI12(I))THEN
#endif
            WRITE(6,*) 'ERROR: compressed index 1 is different'
            WRITE(6,*) 'Item = ',I,                                     &
     &                  ' File 1 = ',CFI11(I),                          &
     &                  ' File 2 = ',CFI12(I)
! DEPENDS ON: ereport
            CALL EREPORT('MERGE', 1050,                                 &
     &       'compressed index 1 is different')

          END IF
          CFI13(I)=CFI11(I)
        END DO
      ENDIF

!L 13. Compare compressed index 2

      IF(LEN_CFI21 >  0.OR.LEN_CFI22 >  0)THEN
        WRITE(6,*)' '
        WRITE(6,*)'COMPRESSED INDEX 2:'
        IF(LEN_CFI21 /= LEN_CFI22)THEN
          WRITE(6,*)'WARNING LEN1=',LEN_CFI21,' LEN2=',LEN_CFI22
        ELSE
          LEN_CFI23=LEN_CFI21
        ENDIF
        DO I=1,LEN_CFI21
#if defined(T3E)
          IF(XOR(CFI21(I),CFI22(I)) /= 0) THEN
#else
          IF(CFI21(I) /= CFI22(I))THEN
#endif
            WRITE(6,*) 'ERROR: compressed index 2 is different'
            WRITE(6,*) 'Item = ',I,                                     &
     &                  ' File 1 = ',CFI21(I),                          &
     &                  ' File 2 = ',CFI22(I)
! DEPENDS ON: ereport
            CALL EREPORT('MERGE', 1051,                                 &
     &       'compressed index 2 is different')

          END IF
          CFI23(I)=CFI21(I)
        END DO
      ENDIF

!L 14. Compare compressed index 3

      IF(LEN_CFI31 >  0.OR.LEN_CFI32 >  0)THEN
        WRITE(6,*)' '
        WRITE(6,*)'COMPRESSED INDEX 3:'
        IF(LEN_CFI31 /= LEN_CFI32)THEN
          WRITE(6,*)'WARNING LEN1=',LEN_CFI31,' LEN2=',LEN_CFI32
        ELSE
          LEN_CFI33=LEN_CFI31
        ENDIF
        DO I=1,LEN_CFI31
#if defined(T3E)
          IF(XOR(CFI31(I),CFI32(I)) /= 0) THEN
#else
          IF(CFI31(I) /= CFI32(I))THEN
#endif
            WRITE(6,*) 'ERROR: compressed index 3 is different'
            WRITE(6,*) 'Item = ',I,                                     &
     &                  ' File 1 = ',CFI31(I),                          &
     &                  ' File 2 = ',CFI32(I)
! DEPENDS ON: ereport
            CALL EREPORT('MERGE', 1052,                                 &
     &       'compressed index 3 is different')

          END IF
          CFI33(I)=CFI31(I)
        END DO
      ENDIF

!L 15. Compare lookup tables

      IF(LEN1_LOOKUP1 /= LEN1_LOOKUP2)THEN
        WRITE(6,*)'ERROR lookup tables of different length'
        WRITE(6,*)'LEN1=',LEN1_LOOKUP1,' LEN2=',LEN1_LOOKUP2
! DEPENDS ON: ereport
        CALL EREPORT('MERGE', 1053,                                     &
     &    'lookup tables of different length')

      ENDIF
      IF(LEN2_LOOKUP1 >  0.OR.LEN2_LOOKUP2 >  0)THEN
        WRITE(6,*)' '
        WRITE(6,*)'LOOKUP:'
        JMIN=MIN0(LEN2_LOOKUP1,LEN2_LOOKUP2)
        IDIFF=0
        NDIFF=0

! Read in namelist.                                   .
! NRECF1>=0 If file 2 is to be appended to file 1 after NRECF1
!           records.
! NRECF1<0  If the files are time series and the output file is a
!           time series. The point of overlap is calculated
!           automatically. This is the setting for merging       .
!           boundary datasets
        READ(5,CONTROL)
        IF(NRECF1 >  LEN2_LOOKUP1)THEN
          WRITE(6,*)'ERROR: NRECF1 is larger than LEN2_LOOKUP1'
          WRITE(6,*)' NRECF1 = ',NRECF1,' LEN2_LOOKUP1 = ',LEN2_LOOKUP1
! DEPENDS ON: ereport
          CALL EREPORT('MERGE', 1054,                                   &
     &     'NRECF1 is larger than LEN2_LOOKUP1')

        ELSE IF(NRECF1 >= 0)THEN
          IF(FIXHD1(5) /= 5)THEN
            IDIFF = NRECF1
          ELSE IF(MOD(NRECF1,INTHD1(15)) == 0)THEN
            IDIFF = NRECF1
          ELSE
            WRITE(6,*) 'ERROR: Files are time series.'
            WRITE(6,*) 'NRECF1 must be a multiple of ',INTHD1(15)
            WRITE(6,*) 'NRECF1 = ',NRECF1
! DEPENDS ON: ereport
            CALL EREPORT('MERGE', 1055,                                 &
     &       'Files are time serices, NRECF1 must be multiple')

          END IF
        ELSE
          IF(FIXHD1(10) /= 1)THEN
            WRITE(6,*)'ERROR: File 1 not a time series'
! DEPENDS ON: ereport
            CALL EREPORT('MERGE', 1056,                                 &
     &       'File 1 not a time series')

          END IF
          IF(FIXHD2(10) /= 1)THEN
            WRITE(6,*)'ERROR: File 2 not a time series'
! DEPENDS ON: ereport
            CALL EREPORT('MERGE', 1057,                                 &
     &       'File 2 not a time series')

          END IF

! Compare each lookup record in file 1 with the first looup record
! in file 2. When match is found set IDIFF.
          DO I=1,LEN2_LOOKUP1
            DO J=1,LEN1_LOOKUP1
              IF(LOOKUP1(J,I) /= LOOKUP2(J,1)                           &
     &.AND.(J <= 6.OR.J == 23.OR.J == 26))THEN
                NDIFF=NDIFF+1
              ENDIF
            ENDDO
            IF((NDIFF == 0).AND.(IDIFF == 0))THEN
              WRITE(6,*)                                                &
     &' File 1 lookup record ',I,' matched with File 2 record 1'
              IDIFF=I-1
            ELSE
              NDIFF=0
            ENDIF
          ENDDO

! If first lookup record in file 2 not found in file 1. Abort with
! error message
          IF(IDIFF == 0)THEN
            WRITE(6,*)                                                  &
     &'ERROR First lookup record in file 2 not found in file 1'
            WRITE(6,*) 'Cannot merge files'
! DEPENDS ON: ereport
            CALL EREPORT('MERGE', 1058,                                 &
     &       'First lookup record in file 2 not found in file 1')

          ENDIF
        ENDIF
      ENDIF

! Copy the first IDIFF records from file 1 and the remainder from
! file 2.
      LEN2_LOOKUP3=LEN2_LOOKUP2+IDIFF
      DO I=1,LEN2_LOOKUP3
        DO J=1,LEN1_LOOKUP1
          IF(I <= IDIFF)THEN
            LOOKUP3(J,I)=LOOKUP1(J,I)
          ELSE
            LOOKUP3(J,I)=LOOKUP2(J,I-IDIFF)
          ENDIF
        ENDDO
      ENDDO

!L 16 Ammend header information

! Check and correct fixed header.
      DO J=1,7
        IF(FIXHD3(5) == 5)THEN
          FIXHD3(20+J)=FIXHD1(20+J)   ! First validity time from file 1
          FIXHD3(27+J)=FIXHD2(27+J)   ! Last validity time from file 2
          IF(FIXHD1(20+J) >  FIXHD2(20+J))THEN
            IF(FIXHD1(20+J-1) >= (FIXHD2(20+J-1)))THEN
              WRITE(6,*) 'ERROR: File 2 is earlier than file 1  ',      &
     &                   FIXHD1(20+J),FIXHD2(20+J)
! DEPENDS ON: ereport
              CALL EREPORT('MERGE', 1059,                               &
     &         'File 2 is earliers than file1')

            ENDIF
          ENDIF
        ELSE
          FIXHD3(20+J)=FIXHD1(20+J)
          FIXHD3(27+J)=FIXHD1(27+J)
          IF(FIXHD1(20+J) /= FIXHD2(20+J))THEN
            WRITE(6,*) 'WARNING: Initial data time differs',            &
     &                 FIXHD1(20+J),FIXHD2(20+J)
          ENDIF
          IF(FIXHD1(27+J) /= FIXHD2(27+J))THEN
            WRITE(6,*) 'WARNING: Validity time differs',                &
     &                 FIXHD1(27+J),FIXHD2(27+J)
          ENDIF
        ENDIF
      ENDDO
      FIXHD3(152)=FIXHD2(152)+IDIFF
      FIXHD3(160)=FIXHD3(150)+FIXHD3(151)*FIXHD3(152)

!L 17 Calculate addressing and length of DATA in file 3

! Atmospheric dump dataset or Ancillary dataset
      IF((FIXHD3(2) == 1.OR.FIXHD3(2) == 2).AND.                        &
     &   (FIXHD3(5) <= 2.OR.FIXHD3(5) == 4))THEN
        LEN_DATA3=0
        DO I=1,LEN2_LOOKUP3
          LOOKUP3(NADDR,I)=LEN_DATA3+1
          LEN_DATA3=LEN_DATA3+LOOKUP3(LBLREC,I)
        ENDDO

  ! Call SET_DUMPFILE_ADDRESS to calculate start address
! DEPENDS ON: set_dumpfile_address
        CALL SET_DUMPFILE_ADDRESS(FIXHD3,LEN_FIXHD3,                    &
     &                            LOOKUP3,LEN1_LOOKUP3,LEN2_LOOKUP3,    &
     &                            NUMBER_OF_DATA_WORDS_IN_MEMORY,       &
     &                            NUMBER_OF_DATA_WORDS_ON_DISK,         &
     &                            DISK_ADDRESS)

! Boundary dataset
      ELSEIF(FIXHD3(2) == 1.AND.FIXHD3(5) == 5)THEN
  !   Calcuate start address from header and round it up
  !   to ensure we start on a sector boundary
        DISK_ADDRESS=FIXHD3(160)-1
        DISK_ADDRESS=((DISK_ADDRESS+UM_SECTOR_SIZE-1)/                  &
     &                UM_SECTOR_SIZE)*UM_SECTOR_SIZE
        FIXHD3(160)=DISK_ADDRESS+1

  ! Loop over number of times for which data is present in dataset
        INTHD3(3)=LEN2_LOOKUP3/INTHD3(15)
        LEN_DATA3=0
        DO J=1,INTHD3(3)
          LEN_BUF=0
          MAX_LEN_BUF=0
          DO I=1,INTHD3(15)
            POS=(J-1)*INTHD3(15)+I
            LOOKUP3(LBEGIN,POS)=DISK_ADDRESS+LEN_BUF
            LOOKUP3(LBNREC,POS)=LOOKUP3(LBLREC,POS)
            LOOKUP3(NADDR,POS)=LEN_DATA3+1
            LEN_BUF=LEN_BUF+LOOKUP3(LBLREC,POS)
          END DO
          MAX_LEN_BUF=MAX0(LEN_BUF,MAX_LEN_BUF)
  ! Update disk address and ensure that next time starts
  ! on a sector boundary
          DISK_ADDRESS=DISK_ADDRESS+LEN_BUF
          DISK_ADDRESS=((DISK_ADDRESS+UM_SECTOR_SIZE-1)/                &
     &                  UM_SECTOR_SIZE)*UM_SECTOR_SIZE
          IF(FIXHD3(12) <= 303)THEN
            LEN_DATA3=LEN_DATA3+LEN_BUF/2
          ELSE
            LEN_DATA3=LEN_DATA3+LEN_BUF
          END IF
        ENDDO
      END IF
      FIXHD3(161)=LEN_DATA3
!L 18. Print out header for file 3 and check for consistency

! DEPENDS ON: pr_fixhd
      CALL PR_FIXHD(FIXHD3,LEN_FIXHD3,LEN_INTHD3,LEN_REALHD3            &
     &,LEN1_LEVDEPC3,LEN2_LEVDEPC3,LEN1_ROWDEPC3,LEN2_ROWDEPC3          &
     &,LEN1_COLDEPC3,LEN2_COLDEPC3,LEN1_FLDDEPC3,LEN2_FLDDEPC3          &
     &,LEN_EXTCNST3,LEN_DUMPHIST3,LEN_CFI13,LEN_CFI23,LEN_CFI33         &
     &,LEN1_LOOKUP3,LEN2_LOOKUP3,LEN_DATA3                              &
     &,ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*)CMESSAGE,ICODE
! DEPENDS ON: ereport
        CALL EREPORT('MERGE', ICODE, CMESSAGE)

      ENDIF
! DEPENDS ON: chk_look
      CALL CHK_LOOK(FIXHD3,LOOKUP3,LEN1_LOOKUP3,LEN_DATA3,              &
#include "argppx.h"
     &              ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*)CMESSAGE,ICODE
! DEPENDS ON: ereport
        CALL EREPORT('MERGE', ICODE, CMESSAGE)

      ENDIF

!L 19. OPEN output file and write out header

      NFTOUT=22
! DEPENDS ON: file_open
      CALL FILE_OPEN(NFTOUT,'FILE3',5,1,0,ERROR)
      IF(ERROR /= 0)THEN
        WRITE(6,*) 'Error opening output file'
! DEPENDS ON: ereport
        CALL EREPORT('MERGE', ICODE, CMESSAGE)

      ENDIF
! DEPENDS ON: writhead
      CALL WRITHEAD(NFTOUT,FIXHD3,LEN_FIXHD3,                           &
     &                INTHD3,LEN_INTHD3,                                &
     &                REALHD3,LEN_REALHD3,                              &
     &                LEVDEPC3,LEN1_LEVDEPC3,LEN2_LEVDEPC3,             &
     &                ROWDEPC3,LEN1_ROWDEPC3,LEN2_ROWDEPC3,             &
     &                COLDEPC3,LEN1_COLDEPC3,LEN2_COLDEPC3,             &
     &                FLDDEPC3,LEN1_FLDDEPC3,LEN2_FLDDEPC3,             &
     &                EXTCNST3,LEN_EXTCNST3,                            &
     &                DUMPHIST3,LEN_DUMPHIST3,                          &
     &                CFI13,LEN_CFI13,                                  &
     &                CFI23,LEN_CFI23,                                  &
     &                CFI33,LEN_CFI33,                                  &
     &                LOOKUP3,LEN1_LOOKUP3,LEN2_LOOKUP3,                &
     &                LEN_DATA3,                                        &
#include "argppx.h"
     &                START_BLOCK,ICODE,CMESSAGE)

!L 19. Write data fields

      WRITE(6,*)' '
      WRITE(6,*)'DATA FIELDS:'
      JMIN=MIN0(LEN2_LOOKUP1,LEN2_LOOKUP2)

      DO I=1,IDIFF

! Read first fields from file 1 and write them to field 3

        IF (LOOKUP1(1,I) /= -99) THEN
        FIELD_ITEM=MOD(LOOKUP1(42,I),1000)
        FIELD_SECT=(LOOKUP1(42,I)-FIELD_ITEM)/1000
        FIELD_MODEL=LOOKUP1(45,I)
! DEPENDS ON: exppxi
        GRID_TYPE=EXPPXI(FIELD_MODEL,FIELD_SECT,FIELD_ITEM,             &
     &                          ppx_grid_type,                          &
#include "argppx.h"
     &                          ICODE,CMESSAGE)

! DEPENDS ON: decompose_smexe
        CALL DECOMPOSE_SMEXE(ROW_LENGTH1, P_ROWS1,                      &
     &                       0,0,TOT_LEVELS)
! DEPENDS ON: change_decomposition
        CALL CHANGE_DECOMPOSITION(decomp_smexe,ICODE)
! DEPENDS ON: readflds
        CALL READFLDS(NFTIN1,1,I,LOOKUP1,LEN1_LOOKUP1,                  &
     &                D1,P_FIELD1,FIXHD1,                               &
#include "argppx.h"
     &                ICODE,CMESSAGE)
! DEPENDS ON: abort_io
        IF(ICODE /= 0)CALL ABORT_IO('MERGE',CMESSAGE,ICODE,NFTIN1)
! DEPENDS ON: writflds
        CALL WRITFLDS(NFTOUT,1,I,LOOKUP3,LEN1_LOOKUP3,                  &
     &                D1,P_FIELD3,FIXHD3,                               &
#include "argppx.h"
     &                ICODE,CMESSAGE)
! DEPENDS ON: abort_io
        IF(ICODE /= 0)CALL ABORT_IO('MERGE',CMESSAGE,ICODE,NFTIN1)
        ENDIF
      ENDDO

! Read remaining fields from file 2 and write them to file 3
      DO I=1,LEN2_LOOKUP2

        IF ((LOOKUP2(1,I) /= -99).AND.(LOOKUP1(1,I) /= -99)) THEN
        FIELD_ITEM=MOD(LOOKUP2(42,I),1000)
        FIELD_SECT=(LOOKUP2(42,I)-FIELD_ITEM)/1000
        FIELD_MODEL=LOOKUP2(45,I)
! DEPENDS ON: exppxi
        GRID_TYPE=EXPPXI(FIELD_MODEL,FIELD_SECT,FIELD_ITEM,             &
     &                          ppx_grid_type,                          &
#include "argppx.h"
     &                          ICODE,CMESSAGE)

!       A check is performed as some grid types are not
!       supported in this version
        IF ((GRID_TYPE >  50.AND.GRID_TYPE <  54).OR.                   &
     &      (GRID_TYPE >  24.AND.GRID_TYPE <  30).OR.                   &
     &      (GRID_TYPE == 21)) THEN
          WRITE(CMESSAGE,*)'FIELD ',I,' WITH GRID_TYPE ',GRID_TYPE,     &
     &              ' NOT SUPPORTED.'
! DEPENDS ON: ereport
          CALL EREPORT('MERGE', 1063, CMESSAGE)

        ENDIF

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN2,1,I,LOOKUP2,LEN1_LOOKUP2,                  &
     &                D2,P_FIELD2,FIXHD2,                               &
#include "argppx.h"
     &                ICODE,CMESSAGE)
! DEPENDS ON: abort_io
        IF(ICODE /= 0)CALL ABORT_IO('MERGE',CMESSAGE,ICODE,NFTIN2)
! DEPENDS ON: writflds
        CALL WRITFLDS(NFTOUT,1,I+IDIFF,LOOKUP3,LEN1_LOOKUP3,            &
     &                D2,P_FIELD3,FIXHD3,                               &
#include "argppx.h"
     &                ICODE,CMESSAGE)
! DEPENDS ON: abort_io
        IF(ICODE /= 0)CALL ABORT_IO('MERGE',CMESSAGE,ICODE,NFTIN2)

        ENDIF
      ENDDO
! Call file close to flush io buffers and close output file
! DEPENDS ON: file_close
      CALL FILE_CLOSE(NFTOUT,'FILE3',5,1,0,ERROR)
      RETURN
      END SUBROUTINE MERGE
#endif
