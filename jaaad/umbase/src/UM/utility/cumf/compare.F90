#if defined(CUMF)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Program  MAIN_COMPARE and Subroutine COMPARE
!LL
!LL  Purpose: Compares two UM atmosphere, ocean, or ancillary files.
!LL           MAIN_COMPARE reads in fixed length and integer
!LL           headers of UM files to be compared, extracts dimensions
!LL           of each file and then passes these values to
!LL           subroutine COMPARE.
!LL
!LL            COMPARE subroutine:
!LL          Compares two UM atmosphere, ocean, or ancillary files.
!LL          COMPARE reads in headers and data fields from files on
!LL          NFTIN1 and NFTIN2, comparing values.
!LL          UNIT 6: If an exact compare is found the message 'OK'
!LL          is written out, otherwise
!LL          i)  if header, all differring values are printed
!LL          ii) if field, 1st 10 differring values are printed plus
!LL              the maximum difference between the fields.
!LL          iii) if field only present in one file, a warning message
!LL               is displayed
!LL          UNIT 7: Number of differences displayed for each header.
!LL                  Number of fields with differences is also
!LL                  displayed along with the number of differences
!LL                  for each field which has differences
!LL
!LL  Written by A. Dickinson 20/03/92
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL  3.3   31/10/93   Dimension of data array set to maximum value
!LL                   Author: A. Dickinson      Reviewer: P.Burton
!LL
!LL   3.3   22/11/93  Compare logical fields correctly. Print integer
!LL                   and logical differences. Do not compare data
!LL                   section for obs files. D. Robinson
!LL   3.3   15/12/93  Skip comparing fields if lookup record
!LL                   contains -99's. Allow compare to continue for
!LL                   files with different no of fields. Do not compare
!LL                   fields packed/compressed via WGDOS/GRIB method.
!LL                   Author: D.M.Goddard     Reviewer: D. Robinson
!LL
!LL   3.3   08/12/93  Extra argument for READFLDS. D. Robinson.
!LL
!LL   3.4   08/09/94  Print real values for LOOKUP 46-64 differences.
!LL                   Compare arrays only if both exist. D. Robinson.
!LL
!LL   3.4   12/12/94  Compare fields if LOOKUP(39) is -1 -2 -3
!LL                   ie Timeseries
!LL   3.5  24/03/95    Changed OPEN to FILE_OPEN  P.Burton
!       3.5  27/06/95  Submodels project. Replace call to RDPPXRF by
!                      function EXPPXC to extract name of diagnostic
!                      item.
!                      Author D.M.Goddard    Reviewer S Swarbrick
!     4.0   06/09/95  Allows comparison of pre-vn4.0 and vn4.0 dumps
!                     contain u and v currents as grid type for these
!                     fields as corrected at vn4.0 from 3 to 13.
!                     Author D.M. Goddard
!     4.0  18/09/95    Changes for submodel project
!   4.1  18/06/96   Changes to cope with changes in STASH addressing
!                   Author D.M. Goddard.
!     4.1  21/03/96    Fields read into correctly typed arrays
!                      Added more detailed output:
!                       - Deviation charts
!                       - Basic statistical analysis
!                            P.Burton
!     4.2  10/05/96    Added some checks to avoid FPE's by
!                      checking for NaN's and using xor for
!                      comparisons. UDG2F402
!                      Author: Bob Carruthers
!     4.2  10/05/96    Extension to process WGDOS packed fields UBC3F402
!                      Author: Bob Carruthers
!     4.3  12/03/97    Correct comparsion of integers
!          24/04/97    Corrections for comparing packed fieldsfiles
!                      Write out position of maximum difference
!                      Author: D.M. Goddard and Richard Barnes
!LL  4.4   Oct. 1997 Changed error handling from routine HDPPXRF
!LL                  so only fatal (+ve) errors are handled.
!LL                                             Shaun de Witt
!     4.4  11/06/97    Changes in print statements to reflect the
!                      well-formed Dumpfile I/O.
!                        Author: Bob Carruthers, Cray Research.
!   4.4  24/10/97   Initialise ICODE as it is no longer
!                   initialised in HDPPXRF
!                   Author D.M. Goddard
!                   + extra write statement for statistics. R.Rawlins
!     4.5  14/07/98    Replaced 'xor' and 'and' bitwise operators for
!                      workstations due to non-portability
!                      (A Van der Wal)
!   4.5  10/11/98   General upgrade to program.
!                   1) Files with different sets of fields can now
!                      be compared.
!                   2) Summary file now contains more information.
!                   Author D.M Goddard
!   5.1  31/03/00   Set value of PACK_CODE2 (bug introduced at vn4.5).
!                   Added call to InitPrintStatus.
!                   Added comma missing from write format. D.P.Matthews
!   5.1  11/05/00   Increase format size. D Robinson.
!LL  5.2  7/11/99  Enable run length encoding of ocean fieldsfiles to
!LL                compress the sequences of mdi values that represent
!LL                the land points. Ian Edmond
!   5.2  05/12/00   Produce difference-map displays for
!                   land-only fields.
!                   E.Leung
!LL  5.3  24/09/01  Portability changes.    Z. Gardner
!   5.3  22/11/01   Enable MPP as the only option for
!                   small executables         E.Leung
!   5.4  03/09/02   Allow landpt only field and LBC to be read in
!                   correctly after MPP removal           E.Leung
!   5.5  07/02/03   Use INTHD 6,7 instead of LOOKUP 18,19 as GLSIZE
!                                                           E.Leung
!   6.0  30/12/03   Set Model Code from 10 to 1 for FieldCalc
!                   diagnostics to use Atmos StashMaster file.
!                   D. Robinson
!   6.0  29/01/04   Print out filenames being compared.
!                   D. Robinson
!   6.0  11/09/03   Add call to gc_exit                     P.Dando
!   6.0  20/01/04   Read land mask only if one is in dump
!                   Output boundary diff maps on arbitrary grid
!                   S.D. Mullerworth
!   6.0  18/06/03   Fix PRINT_DIF_MAP input argument   E.Leung
!   6.1  18/08/04   Include file sx_dt.h moved to correct location
!                   for data initialisation (i.e. in BLKDATA) and remove
!                   repeated declaration of TOT_LEVELS.  P.Dando
!   6.2  22/08/05   Improve T3E construct for FCM. P.Selwood
!   6.2  10/4/2005  Removed calls to ABORT.  J. Gill
!   6.2  18/01/06   Print out field descriptor and level for fields.
!                   D.Robinson
!LL  Programming standard:
!LL
!LL  Logical components covered:
!LL
!LL  System Tasks: F3,F4,F6
!LL
!LL  Documentation: UM Doc Paper F5
!LL
!LL  -----------------------------------------------------------------
!*L  Arguments:-------------------------------------------------------
      SUBROUTINE COMPARE(LEN_FIXHD1,LEN_INTHD1,LEN_REALHD1,             &
     &  LEN1_LEVDEPC1,LEN2_LEVDEPC1,LEN1_ROWDEPC1,                      &
     &  LEN2_ROWDEPC1,LEN1_COLDEPC1,LEN2_COLDEPC1,                      &
     &  LEN1_FLDDEPC1,LEN2_FLDDEPC1,LEN_EXTCNST1,                       &
     &  LEN_DUMPHIST1,LEN_CFI11,LEN_CFI21,LEN_CFI31,                    &
     &  LEN1_LOOKUP1,LEN2_LOOKUP1,LEN_DATA1,P_FIELD1,                   &
     &  P_ROWS1,P_ROWS2,ROW_LENGTH1,ROW_LENGTH2,                        &
     &  LEN_FIXHD2,LEN_INTHD2,LEN_REALHD2,                              &
     &  LEN1_LEVDEPC2,LEN2_LEVDEPC2,LEN1_ROWDEPC2,                      &
     &  LEN2_ROWDEPC2,LEN1_COLDEPC2,LEN2_COLDEPC2,                      &
     &  LEN1_FLDDEPC2,LEN2_FLDDEPC2,LEN_EXTCNST2,                       &
     &  LEN_DUMPHIST2,LEN_CFI12,LEN_CFI22,LEN_CFI32,                    &
     &  LEN1_LOOKUP2,LEN2_LOOKUP2,LEN_DATA2,P_FIELD2                    &
     & ,NFTIN1,NFTIN2,MAX_FIELD_SIZE1,MAX_FIELD_SIZE2,                  &
#include "argppx.h"
     & expand)

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
     &,ROW_LENGTH1                                                      &
     &,MAX_FIELD_SIZE1                                                  &
                       !IN Maximum field size on file 1
     &,expand ! IN set to 1 to expand WGDOS and run length encoded
              ! Fields for comparison.
!
      integer lblrec_1, lblrec_2, length_changed

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
     &,ROW_LENGTH2                                                      &
     &,MAX_FIELD_SIZE2 !IN Maximum field size on file 2

      INTEGER                                                           &
     & NFTIN1                                                           &
                    !IN Unit number for file 1
     &,NFTIN2       !IN Unit number for file 2


! Comdecks: ------------------------------------------------------------
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "clookadd.h"
#include "c_mdi.h"
#include "cstash.h"
#include "parvars.h"
#include "decomptp.h"
#include "sx_size.h"
#include "atm_lsm.h"

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
     & R_D1(MAX_FIELD_SIZE1) ! REAL Array for field on file 1

      INTEGER                                                           &
     & I_D1(MAX_FIELD_SIZE1) ! INTEGER Array for field on file 1

      LOGICAL                                                           &
     & L_D1(MAX_FIELD_SIZE1) ! LOGICAL Array for field on file 1

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
     & R_D2(MAX_FIELD_SIZE2) ! REAL Array for field on file 2

      INTEGER                                                           &
     & I_D2(MAX_FIELD_SIZE1) ! INTEGER Array for field on file 2

      LOGICAL                                                           &
     & L_D2(MAX_FIELD_SIZE1) ! LOGICAL Array for field on file 2

      INTEGER                                                           &
     & PP_XREF(PPXREF_CODELEN)  !PPXREF codes for a given section/item

      LOGICAL                                                           &
     &  LAND_MASK_FOUND  ! Is there a land mask in the dump

! External subroutines called:------------------------------------------
      EXTERNAL ABORT,ABORT_IO,READHEAD,READFLDS,HDPPXRF,GETPPX,         &
     &         DERVSIZE
!*----------------------------------------------------------------------
!*L  Local variables:---------------------------------------------------
      REAL                                                              &
     & MAX_DIFF                                                         &
                 ! Maximum difference between two real fields
     &,RD1,RD2   ! Real variables to equivalent with LD1/ID1 & LD2/ID2
      REAL DIFF_PER,RMS_F1,RMS_F2,RMS_DIFF,A,RCODE

#if defined(T3E)

      integer jrc_nan
!
      integer jrc_mask
      integer deb_mask
!
      data jrc_mask/X'7FF0000000000000'/
      data deb_mask/X'FFF0000000000000'/
!
#endif
      INTEGER                                                           &
     & ICODE                                                            &
                    ! Error return code from subroutines
     &,START_BLOCK                                                      &
                    ! READHEAD argument (not used)
     &,I,J,K,L,M,N                                                      &
                      ! Loop indices
     &,JMIN                                                             &
                    ! Minimum length of two headers
     &,S_ITEM_CODE                                                      &
                    ! STASH item code
     &,SECTION                                                          &
                    ! STASH section number
     &,ID1,ID2                                                          &
                    ! Integer variables to equivalent with RD1 and RD2
     &,N_DIFF                                                           &
                    ! No of differences to be listed
     &,IMAX_DIFF                                                        &
                    ! Maximum difference between two integer fields
     &,PACK_CODE1                                                       &
                    ! Packing code for LOOKUP table 1
     &,PACK_CODE2                                                       &
                    ! Packing code for LOOKUP table 2
     &,MAX_J        ! Location of max.diff
       INTEGER RowNumber
      INTEGER MODEL              !Internal model number
      INTEGER LEN_FIELD          !Number of points in field to be
                                 !compared
      INTEGER N1,N2
      INTEGER OFFSET1
      INTEGER OFFSET2
      INTEGER NUMREC1
      INTEGER NUMREC2
      INTEGER NMISSING1
      INTEGER NMISSING2
      INTEGER IROWDEPC1
      INTEGER IROWDEPC2
      INTEGER IMASK,IMASK_DUMMY,LEN_IO
      INTEGER HALOX1,HALOY1,HALOX2,HALOY2
      INTEGER LAND_POINTS
      INTEGER OLEN_FIELD
      INTEGER FIELD_ITEM1,FIELD_SECT1,FIELD_MODEL1
      INTEGER FIELD_ITEM2,FIELD_SECT2,FIELD_MODEL2
      INTEGER GRID_TYPE1,GRID_TYPE2
      INTEGER INDEX(LEN2_LOOKUP1)
      INTEGER NDIFFER(LEN2_LOOKUP1)
      LOGICAL LMISSING1(LEN2_LOOKUP1)
      LOGICAL LMISSING2(LEN2_LOOKUP2)
       INTEGER      EXPPXI
       CHARACTER*36 EXPPXC

      LOGICAL                                                           &
     & LD1,LD2      ! Logical variables to equivalent with RD1 and RD2

      CHARACTER                                                         &
     & CMESSAGE*100                                                     &
                    ! Character string returned if ICODE  /=  0
     &,PHRASE*(PPXREF_CHARLEN) ! Name of field
      CHARACTER*1 DIFF(MAX_FIELD_SIZE1)
      CHARACTER*200 KEY
      CHARACTER*80 FILENAME ! Name of output file

      EQUIVALENCE (RD1,ID1,LD1) , (RD2,ID2,LD2)

      PARAMETER (N_DIFF=10)
      INTEGER NFT1,NFT2
      PARAMETER (NFT1=22, NFT2=2)
!*----------------------------------------------------------------------

!  0. Open PPXREF file

      ppxRecs=1
      RowNumber=0
      cmessage = ' '
      ICODE=0
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_A',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_A'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('COMPARE', ICODE,                                  &
     &   'Error reading STASHmaster_A')

      END IF
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_O',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_O'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('COMPARE', ICODE,                                  &
     &   'Error reading STASHmaster_O')

      END IF
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_S',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_S'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('COMPARE', ICODE,                                  &
     &   'Error reading STASHmaster_S')

      END IF
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_W',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_W'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('COMPARE', ICODE,                                  &
     &   'Error reading STASHmaster_W')

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
        CALL EREPORT('COMPARE', ICODE,                                  &
     &   CMESSAGE)

      ENDIF

!User STASHmaster
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(0,' ',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('COMPARE', ICODE,                                  &
     &   CMESSAGE)

      ENDIF
! DEPENDS ON: getppx
      CALL GETPPX(0,NFT2,' ',RowNumber,                                 &
#include "argppx.h"
     &            ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('COMPARE', ICODE,                                  &
     &   CMESSAGE)

      ENDIF
! 1: Open output files
!
! Open up unit 7: Summary part one
      CALL GET_FILE(7,FILENAME,80,ICODE)
      OPEN(7,FILE=FILENAME,STATUS='NEW',IOSTAT=ICODE)
      IF (ICODE /= 0) THEN
        WRITE(6,*) 'Can not write to ',FILENAME
      ELSE
        WRITE(6,*) 'OPEN: 7:',FILENAME,'has been created'
      ENDIF
      WRITE(7,*)' COMPARE - SUMMARY MODE'
      WRITE(7,*)'-----------------------'
      WRITE(7,*)' '

! Open up unit 8: Summary part two
      CALL GET_FILE(8,FILENAME,80,ICODE)
      OPEN(8,FILE=FILENAME,STATUS='NEW',IOSTAT=ICODE)
      IF (ICODE /= 0) THEN
        WRITE(6,*) 'Can not write to ',FILENAME
      ELSE
        WRITE(6,*) 'OPEN: 8:',FILENAME,'has been created'
      ENDIF

 ! Open up unit 10
      CALL GET_FILE(10,FILENAME,80,ICODE)
      OPEN(10,FILE=FILENAME,STATUS='NEW',IOSTAT=ICODE)
      IF (ICODE /= 0) THEN
        WRITE(6,*) 'Can not write to ',FILENAME
      ELSE
        WRITE(6,*) 'OPEN: 10:',FILENAME,'has been created'
      ENDIF
      WRITE(10,*)' COMPARE - DIFFERENCE CHARTS'
      WRITE(10,*)'----------------------------'
      WRITE(10,*)' '

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
        CALL EREPORT('COMPARE', ICODE,                                  &
     &   CMESSAGE)

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
        CALL EREPORT('COMPARE', ICODE,                                  &
     &   CMESSAGE)

      ENDIF

!L 3. Compare fixed length headers

      IF(FIXHD1(5) /= FIXHD2(5))THEN
        WRITE(6,'(''WARNING: FIXHD1(5)  = '',I3,'' FIXHD2(5)  = '',I3)')&
     &    FIXHD1(5),FIXHD2(5)
        WRITE(7,'(''WARNING: FIXHD1(5)  = '',I3,'' FIXHD2(5)  = '',I3)')&
     &    FIXHD1(5),FIXHD2(5)
        WRITE(6,'(''         File types are different'')')
        WRITE(7,'(''         File types are different'')')
      END IF
      WRITE(6,*)' '
      WRITE(6,*)'FIXED LENGTH HEADER:'

! Check length of fixed length headers
      JMIN=MIN0(LEN_FIXHD1,LEN_FIXHD2)
      IF(LEN_FIXHD1 /= LEN_FIXHD2)THEN
        WRITE(6,'(''WARNING: LEN_FIXHD1 = '',I3,'' LEN_FIXHD2 = '',I3)')&
     &    LEN_FIXHD1,LEN_FIXHD2
        WRITE(7,'(''WARNING: LEN_FIXHD1 = '',I3,'' LEN_FIXHD2 = '',I3)')&
     &    LEN_FIXHD1,LEN_FIXHD2
        WRITE(6,'(''         Fixed length headers have different '',    &
     &            ''lengths'')')
        WRITE(7,'(''         Fixed length headers have different '',    &
     &            ''lengths'')')
        WRITE(6,'(''         Comparing first '',I3,''elements only'')') &
     &    JMIN
        WRITE(7,'(''         Comparing first '',I3,''elements only'')') &
     &    JMIN
      END IF

! Check fixed length header
      IF(FIXHD1(152) == FIXHD2(152))THEN
        IF(FIXHD1(160) /= FIXHD2(160))THEN
          WRITE(6,'(''WARNING: LEN1 = '',i9,'' and LEN2 = '',i9)')      &
     &      FIXHD1(160),FIXHD2(160)
          WRITE(7,'(''WARNING: LEN1 = '',i9,'' and LEN2 = '',i9)')      &
     &      FIXHD1(160),FIXHD2(160)
          WRITE(6,'(''         Data start address differs'')')
          WRITE(7,'(''         Data start address differs'')')
          WRITE(6,'(''         Possibly due to comparing old and new '',&
     &              ''format UM dumps or fieldsfiles'')')
          WRITE(7,'(''         Possibly due to comparing old and new '',&
     &              ''format UM dumps or fieldsfiles'')')
        ELSE IF(FIXHD1(161) /= FIXHD2(161))THEN
          WRITE(6,'(''WARNING: LEN1 = '',i9,'' and LEN2 = '',i9)')      &
     &      FIXHD1(161),FIXHD2(161)
          WRITE(7,'(''WARNING: LEN1 = '',i9,'' and LEN2 = '',i9)')      &
     &      FIXHD1(161),FIXHD2(161)
          WRITE(6,'(''         Length of data differs'')')
          WRITE(7,'(''         Length of data differs'')')
          WRITE(6,'(''         Possibly due to comparing old and new '',&
     &              ''format UM dumps or fieldsfiles'')')
          WRITE(7,'(''         Possibly due to comparing old and new '',&
     &              ''format UM dumps or fieldsfiles'')')
        END IF
      END IF

      K = 0
      length_changed=0
      DO I=1,JMIN
        IF(FIXHD1(I) /= FIXHD2(I))THEN
        WRITE(6,'(''ITEM = '',i4,''  Values = '',i11,'' and '',i11)')   &
     &    i, fixhd1(i), fixhd2(i)
        K = K + 1
        ENDIF
      ENDDO

      IF(K == 0) WRITE(6,*) 'OK'
      WRITE(8,*) 'FIXED LENGTH HEADER:        ',                        &
     &           'Number of differences = ',K

!L 4. Compare integer headers

      IF(LEN_INTHD1 >  0.OR.LEN_INTHD2 >  0)THEN
        WRITE(6,*)' '
        WRITE(6,*)'INTEGER HEADER:'
        IF(LEN_INTHD1 /= LEN_INTHD2)THEN
          WRITE(6,*)'WARNING LEN1=',LEN_INTHD1,' LEN2=',LEN_INTHD2
        ENDIF
        JMIN=MIN0(LEN_INTHD1,LEN_INTHD2)
        K=0
        DO I=1,JMIN
          IF(INTHD1(I) /= INTHD2(I))THEN
            K=K+1
            WRITE(6,*)'ITEM=',I,INTHD1(I),INTHD2(I)
          ENDIF
        ENDDO
      ENDIF

      IF(K == 0) WRITE(6,*) 'OK'
      WRITE(8,*) 'INTEGER HEADER:             ',                        &
     &           'Number of differences = ',K
      L=K

!L 5. Compare real headers

      IF(LEN_REALHD1 >  0.OR.LEN_REALHD2 >  0)THEN
        WRITE(6,*)' '
        WRITE(6,*)'REAL HEADER:'
        IF(LEN_REALHD1 /= LEN_REALHD2)THEN
          WRITE(6,*)'WARNING LEN1=',LEN_REALHD1,' LEN2=',LEN_REALHD2
        ENDIF
        JMIN=MIN0(LEN_REALHD1,LEN_REALHD2)
        K=0
        DO I=1,JMIN
#if defined(T3E)
          IF(XOR(REALHD1(I),REALHD2(I)) /= 0) THEN
#else
          IF(REALHD1(I) /= REALHD2(I))THEN
#endif
            K=K+1
            WRITE(6,*)'ITEM=',I,REALHD1(I),REALHD2(I)
          ENDIF
        ENDDO
      ENDIF

      IF(K == 0) WRITE(6,*) 'OK'
      WRITE(8,*) 'REAL HEADER:                ',                        &
     &           'Number of differences = ',K
      L=L+K

!L 6. Compare level dependent constants

      WRITE(6,*)' '
      WRITE(6,*)'LEVEL DEPENDENT CONSTS:'
      IF(FIXHD1(110) >  0 .AND. FIXHD2(110) >  0) THEN
      IF(LEN1_LEVDEPC1 /= LEN1_LEVDEPC2)THEN
        WRITE(6,*)'ERROR : different number of levels'
        WRITE(6,*)'LEV1=',LEN1_LEVDEPC1,' LEV2=',LEN1_LEVDEPC2
! DEPENDS ON: ereport
        CALL EREPORT('COMPARE', 1009,                                   &
     &   'Different number of levels')

      ELSEIF(LEN2_LEVDEPC1 >  0.OR.LEN2_LEVDEPC2 >  0)THEN
        IF(LEN2_LEVDEPC1 /= LEN2_LEVDEPC2)THEN
          WRITE(6,*)'WARNING LEN1=',LEN2_LEVDEPC1,' LEN2=',LEN2_LEVDEPC2
        ENDIF
        JMIN=MIN0(LEN2_LEVDEPC1,LEN2_LEVDEPC2)
        K=0
        DO I=1,JMIN
          DO J=1,LEN1_LEVDEPC1
#if defined(T3E)
            IF(XOR(LEVDEPC1((I-1)*LEN1_LEVDEPC1+J),                     &
     &        LEVDEPC2((I-1)*LEN1_LEVDEPC1+J)) /= 0)THEN
#else
            IF(LEVDEPC1((I-1)*LEN1_LEVDEPC1+J) /=                       &
     &        LEVDEPC2((I-1)*LEN1_LEVDEPC1+J))THEN
#endif
              K=K+1
              WRITE(6,*)'LEVEL=',J,'ITEM=',I,                           &
     &        LEVDEPC1((I-1)*LEN1_LEVDEPC1+J),                          &
     &        LEVDEPC2((I-1)*LEN1_LEVDEPC1+J)
           ENDIF
          ENDDO
        ENDDO

        IF(K == 0) WRITE(6,*) 'OK'
        WRITE(8,*) 'LEVEL DEPENDENT CONSTANTS:  ',                      &
     &             'Number of differences = ',K
        L=L+K
      ENDIF
      ELSE
        WRITE(6,*)'No comparison done'
        IF (FIXHD1(110) <= 0) WRITE(6,*)'No array in FILE 1'
        IF (FIXHD2(110) <= 0) WRITE(6,*)'No array in FILE 2'
      ENDIF

!L 7. Compare row dependent constants

      WRITE(6,*)' '
      WRITE(6,*)'ROW DEPENDENT CONSTS:'
      IF(FIXHD1(115) >  0 .AND. FIXHD2(115) >  0) THEN
      IF(LEN1_ROWDEPC1 /= LEN1_ROWDEPC2)THEN
        WRITE(6,*)'ERROR : different number of rows'
        WRITE(6,*)'ROW1=',LEN1_ROWDEPC1,' ROW2=',LEN1_ROWDEPC2
! DEPENDS ON: ereport
        CALL EREPORT('COMPARE', 1010,                                   &
     &   'Different number of rows')

      ELSEIF(LEN2_ROWDEPC1 >  0.OR.LEN2_ROWDEPC2 >  0)THEN
        IF(LEN2_ROWDEPC1 /= LEN2_ROWDEPC2)THEN
          WRITE(6,*)'WARNING different second dimension'
          WRITE(6,*)'LEN1=',LEN2_ROWDEPC1,' LEN2=',LEN2_ROWDEPC2
        ENDIF
        JMIN=MIN0(LEN2_ROWDEPC1,LEN2_ROWDEPC2)
        K=0
        DO I=1,JMIN
          DO J=1,LEN1_ROWDEPC1
#if defined(T3E)
            IF(XOR(ROWDEPC1((I-1)*LEN1_ROWDEPC1+J),                     &
     &         ROWDEPC2((I-1)*LEN1_ROWDEPC1+J)) /= 0)THEN
#else
            IF(ROWDEPC1((I-1)*LEN1_ROWDEPC1+J) /=                       &
     &         ROWDEPC2((I-1)*LEN1_ROWDEPC1+J))THEN
#endif
              K=K+1
              WRITE(6,*)'ROW=',I,'ITEM=',J,                             &
     &        ROWDEPC1((I-1)*LEN1_ROWDEPC1+J),                          &
     &        ROWDEPC2((I-1)*LEN1_ROWDEPC1+J)
            ENDIF
          ENDDO
        ENDDO

        IF(K == 0) WRITE(6,*) 'OK'
        WRITE(8,*) 'ROW DEPENDENT CONSTANTS:    ',                      &
     &             'Number of differences = ',K
        L=L+K
      ENDIF
      ELSE
        WRITE(6,*)'No comparison done'
        IF (FIXHD1(115) <= 0) WRITE(6,*)'No array in FILE 1'
        IF (FIXHD2(115) <= 0) WRITE(6,*)'No array in FILE 2'
      ENDIF

!L 8. Compare column dependent constants

      WRITE(6,*)' '
      WRITE(6,*)'COLUMN DEPENDENT CONSTS:'
      IF(FIXHD1(120) >  0 .AND. FIXHD2(120) >  0) THEN
      IF(LEN1_COLDEPC1 /= LEN1_COLDEPC2)THEN
        WRITE(6,*)'ERROR : different number of columns.'
        WRITE(6,*)'COL1=',LEN1_COLDEPC1,' COL2=',LEN1_COLDEPC2
! DEPENDS ON: ereport
        CALL EREPORT('COMPARE', 1011,                                   &
     &   'Different number of columns')

      ELSEIF(LEN2_COLDEPC1 >  0.OR.LEN2_COLDEPC2 >  0)THEN
        IF(LEN2_COLDEPC1 /= LEN2_COLDEPC2)THEN
          WRITE(6,*)'WARNING LEN1=',LEN2_COLDEPC1,' LEN2=',LEN2_COLDEPC2
        ENDIF
        JMIN=MIN0(LEN2_COLDEPC1,LEN2_COLDEPC2)
        K=0
        DO I=1,JMIN
          DO J=1,LEN1_COLDEPC1
#if defined(T3E)
            IF(XOR(COLDEPC1((I-1)*LEN1_COLDEPC1+J),                     &
     &        COLDEPC2((I-1)*LEN1_COLDEPC1+J)) /= 0) THEN
#else
            IF(COLDEPC1((I-1)*LEN1_COLDEPC1+J) /=                       &
     &        COLDEPC2((I-1)*LEN1_COLDEPC1+J))THEN
#endif
              K=K+1
              WRITE(6,*)'COL=',I,'ITEM=',J,                             &
     &        COLDEPC1((I-1)*LEN1_COLDEPC1+J),                          &
     &        COLDEPC2((I-1)*LEN1_COLDEPC1+J)
            ENDIF
          ENDDO
        ENDDO

        IF(K == 0) WRITE(6,*) 'OK'
        WRITE(8,*) 'COLUMN DEPENDENT CONSTANTS: ',                      &
     &             'Number of differences = ',K
        L=L+K
      ENDIF
      ELSE
        WRITE(6,*)'No comparison done'
        IF (FIXHD1(120) <= 0) WRITE(6,*)'No array in FILE 1'
        IF (FIXHD2(120) <= 0) WRITE(6,*)'No array in FILE 2'
      ENDIF

!L 9. Compare field dependent constants

      WRITE(6,*)' '
      WRITE(6,*)'FIELD DEPENDENT CONSTS:'
      IF(FIXHD1(125) >  0 .AND. FIXHD2(125) >  0) THEN
      IF(LEN1_FLDDEPC1 /= LEN1_FLDDEPC2)THEN
        WRITE(6,*)'ERROR : different number of fields.'
        WRITE(6,*)'FLD1=',LEN1_FLDDEPC1,' FLD2=',LEN1_FLDDEPC2
! DEPENDS ON: ereport
        CALL EREPORT('COMPARE', 1012,                                   &
     &   'Different number of fields')

      ELSEIF(LEN2_FLDDEPC1 >  0.OR.LEN2_FLDDEPC2 >  0)THEN
        IF(LEN2_FLDDEPC1 /= LEN2_FLDDEPC2)THEN
          WRITE(6,*)'WARNING LEN1=',LEN2_FLDDEPC1,' LEN2=',LEN2_FLDDEPC2
        ENDIF
        JMIN=MIN0(LEN2_FLDDEPC1,LEN2_FLDDEPC2)
        K=0
        DO I=1,JMIN
          DO J=1,LEN1_FLDDEPC1
#if defined(T3E)
            IF(XOR(FLDDEPC1((I-1)*LEN1_FLDDEPC1+J),                     &
     &        FLDDEPC2((I-1)*LEN1_FLDDEPC1+J)) /= 0) THEN
#else
            IF(FLDDEPC1((I-1)*LEN1_FLDDEPC1+J) /=                       &
     &        FLDDEPC2((I-1)*LEN1_FLDDEPC1+J))THEN
#endif
             K=K+1
             WRITE(6,*)'FIELD=',J,'ITEM=',I,                            &
     &       FLDDEPC1((I-1)*LEN1_FLDDEPC1+J),                           &
     &       FLDDEPC2((I-1)*LEN1_FLDDEPC1+J)
           ENDIF
          ENDDO
        ENDDO

        IF(K == 0) WRITE(6,*) 'OK'
        WRITE(8,*) 'FIELD DEPENDENT CONSTANTS:  ',                      &
     &             'Number of differences = ',K

        L=L+K
      ENDIF
      ELSE
        WRITE(6,*)'No comparison done'
        IF (FIXHD1(125) <= 0) WRITE(6,*)'No array in FILE 1'
        IF (FIXHD2(125) <= 0) WRITE(6,*)'No array in FILE 2'
      ENDIF

!L 10. Compare extra constants

      WRITE(6,*)' '
      WRITE(6,*)'EXTRA CONSTANTS:'
      IF(FIXHD1(130) >  0 .AND. FIXHD2(130) >  0) THEN
      IF(LEN_EXTCNST1 >  0.OR.LEN_EXTCNST2 >  0)THEN
        IF(LEN_EXTCNST1 /= LEN_EXTCNST2)THEN
          WRITE(6,*)'WARNING LEN1=',LEN_EXTCNST1,' LEN2=',LEN_EXTCNST2
        ENDIF
        JMIN=MIN0(LEN_EXTCNST1,LEN_EXTCNST2)
        K=0
        DO I=1,JMIN
#if defined(T3E)
          IF(XOR(EXTCNST1(I),EXTCNST2(I)) /= 0) THEN
#else
          IF(EXTCNST1(I) /= EXTCNST2(I))THEN
#endif
            K=K+1
            WRITE(6,*)'ITEM=',I,EXTCNST1(I),EXTCNST2(I)
          ENDIF
        ENDDO

        IF(K == 0) WRITE(6,*) 'OK'
        WRITE(8,*) 'EXTRA CONSTANTS:            ',                      &
     &             'Number of differences = ',K

        L=L+K
      ENDIF
      ELSE
        WRITE(6,*)'No comparison done'
        IF (FIXHD1(130) <= 0) WRITE(6,*)'No array in FILE 1'
        IF (FIXHD2(130) <= 0) WRITE(6,*)'No array in FILE 2'
      ENDIF

!L 11. Compare dump history

      WRITE(6,*)' '
      WRITE(6,*)'HISTORY BLOCK:'
      IF(FIXHD1(135) >  0 .AND. FIXHD2(135) >  0) THEN
      IF(LEN_DUMPHIST1 >  0.OR.LEN_DUMPHIST2 >  0)THEN
        IF(LEN_DUMPHIST1 /= LEN_DUMPHIST2)THEN
          WRITE(6,*)'WARNING LEN1=',LEN_DUMPHIST1,' LEN2=',LEN_DUMPHIST2
        ENDIF
        JMIN=MIN0(LEN_DUMPHIST1,LEN_DUMPHIST2)
        K=0
        DO I=1,JMIN
#if defined(T3E)
          IF(XOR(DUMPHIST1(I),DUMPHIST2(I)) /= 0) THEN
#else
          IF(DUMPHIST1(I) /= DUMPHIST2(I))THEN
#endif
            K=K+1
            WRITE(6,*)'ITEM=',I,DUMPHIST1(I),DUMPHIST2(I)
          ENDIF
        ENDDO

        IF(K == 0) WRITE(6,*) 'OK'
        WRITE(8,*) 'HISTORY BLOCK:              ',                      &
     &             'Number of differences = ',K
        L=L+K
      ENDIF
      ELSE
        WRITE(6,*)'No comparison done'
        IF (FIXHD1(135) <= 0) WRITE(6,*)'No array in FILE 1'
        IF (FIXHD2(135) <= 0) WRITE(6,*)'No array in FILE 2'
      ENDIF

!L 12. Compare compressed index 1

      WRITE(6,*)' '
      WRITE(6,*)'COMPRESSED INDEX 1:'
      IF(FIXHD1(140) >  0 .AND. FIXHD2(140) >  0) THEN
      IF(LEN_CFI11 >  0.OR.LEN_CFI12 >  0)THEN
        IF(LEN_CFI11 /= LEN_CFI12)THEN
          WRITE(6,*)'WARNING LEN1=',LEN_CFI11,' LEN2=',LEN_CFI12
        ENDIF
        JMIN=MIN0(LEN_CFI11,LEN_CFI12)
        K=0
        DO I=1,JMIN
#if defined(T3E)
          IF(XOR(CFI11(I),CFI12(I)) /= 0) THEN
#else
          IF(CFI11(I) /= CFI12(I))THEN
#endif
            K=K+1
            WRITE(6,*)'ITEM=',I,CFI11(I),CFI12(I)
          ENDIF
        ENDDO

        IF(K == 0) WRITE(6,*) 'OK'
        WRITE(8,*) 'COMPRESSED INDEX 1:         ',                      &
     &             'Number of differences = ',K
        L=L+K
      ENDIF
      ELSE
        WRITE(6,*)'No comparison done'
        IF (FIXHD1(140) <= 0) WRITE(6,*)'No array in FILE 1'
        IF (FIXHD2(140) <= 0) WRITE(6,*)'No array in FILE 2'
      ENDIF

!L 13. Compare compressed index 2

      WRITE(6,*)' '
      WRITE(6,*)'COMPRESSED INDEX 2:'
      IF(FIXHD1(142) >  0 .AND. FIXHD2(142) >  0) THEN
      IF(LEN_CFI21 >  0.OR.LEN_CFI22 >  0)THEN
        IF(LEN_CFI21 /= LEN_CFI22)THEN
          WRITE(6,*)'WARNING LEN1=',LEN_CFI21,' LEN2=',LEN_CFI22
        ENDIF
        JMIN=MIN0(LEN_CFI21,LEN_CFI22)
        K=0
        DO I=1,JMIN
#if defined(T3E)
          IF(XOR(CFI21(I),CFI22(I)) /= 0) THEN
#else
          IF(CFI21(I) /= CFI22(I))THEN
#endif
            K=K+1
            WRITE(6,*)'ITEM=',I,CFI21(I),CFI22(I)
          ENDIF
        ENDDO

        IF(K == 0) WRITE(6,*) 'OK'
        WRITE(8,*) 'COMPRESSED INDEX 2:         ',                      &
     &             'Number of differences = ',K
        L=L+K
      ENDIF
      ELSE
        WRITE(6,*)'No comparison done'
        IF (FIXHD1(142) <= 0) WRITE(6,*)'No array in FILE 1'
        IF (FIXHD2(142) <= 0) WRITE(6,*)'No array in FILE 2'
      ENDIF

!L 14. Compare compressed index 3

      WRITE(6,*)' '
      WRITE(6,*)'COMPRESSED INDEX 3:'
      IF(FIXHD1(144) >  0 .AND. FIXHD2(144) >  0) THEN
      IF(LEN_CFI31 >  0.OR.LEN_CFI32 >  0)THEN
        IF(LEN_CFI31 /= LEN_CFI32)THEN
          WRITE(6,*)'WARNING LEN1=',LEN_CFI31,' LEN2=',LEN_CFI32
        ENDIF
        JMIN=MIN0(LEN_CFI31,LEN_CFI32)
        K=0
        DO I=1,JMIN
#if defined(T3E)
          IF(XOR(CFI31(I),CFI32(I)) /= 0) THEN
#else
          IF(CFI31(I) /= CFI32(I))THEN
#endif
            K=K+1
            WRITE(6,*)'ITEM=',I,CFI31(I),CFI32(I)
          ENDIF
        ENDDO

        IF(K == 0) WRITE(6,*) 'OK'
        WRITE(8,*) 'COMPRESSED INDEX 3:         ',                      &
     &             'Number of differences = ',K
        L=L+K
      ENDIF
      ELSE
        WRITE(6,*)'No comparison done'
        IF (FIXHD1(144) <= 0) WRITE(6,*)'No array in FILE 1'
        IF (FIXHD2(144) <= 0) WRITE(6,*)'No array in FILE 2'
      ENDIF

!L 15. Compare lookup tables

      IF(LEN1_LOOKUP1 /= LEN1_LOOKUP2)THEN
        WRITE(6,*)'ERROR first dimensions of lookup tables different'
        WRITE(6,*)'LEN1=',LEN1_LOOKUP1,' LEN2=',LEN1_LOOKUP2
! DEPENDS ON: ereport
        CALL EREPORT('COMPARE', 1013,                                   &
     &   'First dimensions of lookup tables different')

      ENDIF
      IF(LEN2_LOOKUP1 >  0.OR.LEN2_LOOKUP2 >  0)THEN
        WRITE(6,*)' '
        WRITE(6,*)'LOOKUP:'

! Check length of lookup tables
        IF(FIXHD1(5) == 3)THEN
          DO I=1,LEN2_LOOKUP1
            IF(LOOKUP1(1,I) /= -99)NUMREC1=I
          END DO
          DO I=1,LEN2_LOOKUP2
            IF(LOOKUP2(1,I) /= -99)NUMREC2=I
          END DO
        ELSE
          NUMREC1=LEN2_LOOKUP1
          NUMREC2=LEN2_LOOKUP2
        END IF

        IF(LEN2_LOOKUP1 /= LEN2_LOOKUP2)THEN
           WRITE(6,'(''WARNING LEN1 = '',i9,'' and LEN2 = '',i9)')      &
     &       len2_lookup1, len2_lookup2
          IF(FIXHD1(5) == 3)THEN
            WRITE(6,'(''Fieldsfile file1 contains '',i5,'' fields '',   &
     &                ''and '',I5,'' empty records'')')                 &
     &               NUMREC1,LEN2_LOOKUP1-NUMREC1
            WRITE(7,'(''Fieldsfile file1 contains '',i5,'' fields '',   &
     &                ''and '',I5,'' empty records'')')                 &
     &               NUMREC1,LEN2_LOOKUP1-NUMREC1
            WRITE(6,'(''Fieldsfile file2 contains '',i5,'' fields '',   &
     &                ''and '',I5,'' empty records'')')                 &
     &               NUMREC2,LEN2_LOOKUP2-NUMREC2
            WRITE(7,'(''Fieldsfile file2 contains '',i5,'' fields '',   &
     &                ''and '',I5,'' empty records'')')                 &
     &               NUMREC2,LEN2_LOOKUP2-NUMREC2
            IF( NUMREC1 == NUMREC2)THEN
              WRITE(6,*) 'Files contain same number of fields'
              WRITE(7,*) 'Files contain same number of fields'
              IF(LEN2_LOOKUP1 == NUMREC1)THEN
                WRITE(6,*) 'Empty records at the end of file1 ',        &
     &                     'have probably been removed by convieee'
                WRITE(7,*) 'Empty records at the end of file1 ',        &
     &                     'have probably been removed by convieee'
              ELSE IF(LEN2_LOOKUP2 == NUMREC2)THEN
                WRITE(6,*) 'Empty records at the end of file2 ',        &
     &                     'have probably been removed by convieee'
                WRITE(7,*) 'Empty records at the end of file2 ',        &
     &                     'have probably been removed by convieee'
              END IF
            END IF
          END IF
        END IF

! Build cross reference index
        OFFSET1=0
        OFFSET2=0
        DO I=1,LEN2_LOOKUP1
          INDEX(I)    = 0
          LMISSING1(I) = .TRUE.
        END DO
        DO I=1,LEN2_LOOKUP2
          LMISSING2(I) = .TRUE.
        END DO
        DO I=1,NUMREC1
          N1=I+OFFSET1
          N2=I+OFFSET2
          IF(LOOKUP1(ITEM_CODE,N1) == LOOKUP2(ITEM_CODE,N2))THEN
            INDEX(I)      =  N2
            LMISSING1(I)  = .FALSE.
            LMISSING2(N2) = .FALSE.
          ELSE
            DO J=N2+1,NUMREC2
              IF(INDEX(I) == 0)THEN
                IF(LOOKUP1(ITEM_CODE,N1) ==                             &
     &             LOOKUP2(ITEM_CODE,J))THEN
                  OFFSET2      =  OFFSET2+J-N2
                  INDEX(I)     =  J
                  LMISSING1(I) = .FALSE.
                  LMISSING2(J) = .FALSE.
                END IF
              END IF
            END DO
            IF(INDEX(I) == 0)THEN
              OFFSET2=OFFSET2-1
            END IF
          END IF
        END DO

        NMISSING1=0
        DO I=1,LEN2_LOOKUP1
          IF(LMISSING1(I).AND.LOOKUP1(1,I) /= -99)THEN
            NMISSING1=NMISSING1+1
            WRITE(6,'(''WARNING: Field '',I5,'' of file1 '',            &
     &                            ''has no match in file2'')') I
          END IF
        END DO
        NMISSING2=0
        DO I=1,LEN2_LOOKUP2
          IF(LMISSING2(I).AND.LOOKUP2(1,I) /= -99)THEN
            NMISSING2=NMISSING2+1
            WRITE(6,'(''WARNING: Field '',I5,'' of file2 '',            &
     &                            ''has no match in file1'')') I
          END IF
        END DO

        K=0
        DO I=1,NUMREC1
          IF(.NOT.LMISSING1(I).AND.LOOKUP1(1,I) /= -99)THEN
            DO J=1,LEN1_LOOKUP1
              IF(LOOKUP1(J,I) /= LOOKUP2(J,INDEX(I)))THEN
                K=K+1
                ID1=LOOKUP1(J,I)
                ID2=LOOKUP2(J,INDEX(I))
                IF (J >= 46 .AND. J <= 64) THEN
                  WRITE(6,'(''Header1: '',I5,'' Header2: '',I5,         &
     &                      '' Item: '',I3,'' Values: '',F12.5,F12.5)') &
     &                   I,INDEX(I),J,RD1,RD2
                ELSE
                  WRITE(6,'(''Header1: '',I5,'' Header2: '',I5,         &
     &                    '' Item: '',I3,'' Values: '',I8,I8)')         &
     &                   I,INDEX(I),J,ID1,ID2
                END IF
              END IF
            END DO
          END IF
        END DO

        IF(K == 0) WRITE(6,*) 'OK'
        WRITE(7,'(''Number of fields in file 1 = '',I5)') NUMREC1
        WRITE(7,'(''Number of fields in file 2 = '',I5)') NUMREC2
        WRITE(7,'(''Number of fields compared  = '',I5)')               &
     &        NUMREC1-NMISSING1
        IF(NMISSING1 /= 0)THEN
          WRITE(7,'(''Number of fields from file 1 omitted from '',     &
     &              ''comparison = '',I5)') NMISSING1
        END IF
        IF(NMISSING2 /= 0)THEN
          WRITE(7,'(''Number of fields from file 2 omitted from '',     &
     &              ''comparison = '',I5)') NMISSING2
        END IF
        WRITE(8,*) 'LOOKUP:                     ',                      &
     &             'Number of differences = ',K
        L=L+K
      END IF

!L 16. Compare data fields

!L Get decompostion information
!L ----------------------------
!L The two files should have same resolution
!L TOT_LEVELS not used in SX
! DEPENDS ON: decompose_smexe
        CALL DECOMPOSE_SMEXE(ROW_LENGTH1, P_ROWS1,                      &
     &                       0,0,TOT_LEVELS)
! DEPENDS ON: change_decomposition
        CALL CHANGE_DECOMPOSITION(decomp_smexe,ICODE)

! If it is an instantenous atmospheric dump
      IF(((FIXHD1(5) == 1).AND.(FIXHD2(5) == 1)).AND.                   &
     &   ((FIXHD1(2) == 1).AND.(FIXHD2(2) == 1)))THEN

!L Pulling out land-sea mask
!  -------------------------
! This assumes land-sea mask of 2 files are identical
        LAND_MASK_FOUND=.FALSE.
!       Check whether one exists before trying to read it
        DO I=1,LEN2_LOOKUP1
          IF (LOOKUP1(ITEM_CODE,i)  ==  30) THEN
            LAND_MASK_FOUND=.TRUE.
          END IF
        END DO
        IF (LAND_MASK_FOUND) THEN
! DEPENDS ON: read_land_sea
          CALL READ_LAND_SEA(NFTIN1,RCODE,LOOKUP1,                      &
     &      LEN1_LOOKUP1,LEN2_LOOKUP1,FIXHD1,LEN_FIXHD1)


!         Check error code from read_land_sea
          IF (RCODE /= -1.0) THEN
            WRITE (6,*) 'Error in READ_LAND_SEA.'
            WRITE (6,*) 'Return code from READ_LAND_SEA ',rcode
            ICODE = 200
            WRITE (CMESSAGE,*) 'DRLANDF1 : Error in READ_LAND_SEA.'
            GO TO 9999
          ENDIF
        END IF
      ENDIF
! Extract information of grid type for each file from STASHmaster

      WRITE(6,*)' '
      WRITE(6,*)'DATA FIELDS:'

      M=0
      N=0
      NDIFFER(:) = 0
      DO I=1,NUMREC1   ! Begin loop over number of fields in file1

        S_ITEM_CODE=MOD(LOOKUP1(42,I),1000)
        SECTION=(LOOKUP1(42,I)-S_ITEM_CODE)/1000
        IF(FIXHD1(12) >= 305)THEN
          MODEL=LOOKUP1(45,I)
        ELSEIF(S_ITEM_CODE <= 100.OR.                                   &
     &        (S_ITEM_CODE >= 200.AND.S_ITEM_CODE <= 205))THEN
          MODEL=1
        ELSEIF((S_ITEM_CODE >  100.AND.S_ITEM_CODE <= 176).OR.          &
     &         (S_ITEM_CODE >= 180.AND.S_ITEM_CODE <= 200))THEN
          MODEL=2
        ELSEIF((S_ITEM_CODE >= 177.AND.S_ITEM_CODE <= 179).OR.          &
     &         (S_ITEM_CODE >= 210.AND.S_ITEM_CODE <= 212))THEN
          MODEL=3
        END IF

!       All diagnostics under model code of 10 are in section 20
!       of Atmos StashMaster file.
        If (MODEL == 10) Then
          MODEL = 1
        End If

! DEPENDS ON: exppxc
        PHRASE=EXPPXC(MODEL,SECTION,S_ITEM_CODE,                        &
#include "argppx.h"
     &                ICODE,CMESSAGE)
        IF(ICODE /= 0)THEN
           PHRASE='NON-STANDARD FIELD'
        END IF

        IF(.NOT.LMISSING1(I))THEN
           M=INDEX(I)

        IF((LOOKUP1(42,I) == 28.OR.LOOKUP1(42,I) == 29).AND.            &
     &     (FIXHD1(12) /= FIXHD2(12).AND.                               &
     &     (FIXHD1(12) >= 400.OR.FIXHD2(12) >= 400)))THEN
         LEN_FIELD=MIN0(LOOKUP1(15,I),LOOKUP2(15,M))
        ELSE
          LEN_FIELD=LOOKUP1(15,I)
        END IF
        IF(FIXHD1(12) <  0.AND.FIXHD1(5) /= 3)LOOKUP1(30,I)=0
        IF(FIXHD2(12) <  0.AND.FIXHD2(5) /= 3)LOOKUP2(30,M)=0
        IF (LOOKUP1(1,I) /= -99 .AND. LOOKUP2(1,M) /= -99) THEN

        FIELD_ITEM1=MOD(LOOKUP1(42,I),1000)
        FIELD_SECT1=(LOOKUP1(42,I)-FIELD_ITEM1)/1000
        FIELD_MODEL1=LOOKUP1(45,I)

        If (FIELD_MODEL1 == 10) Then
          FIELD_MODEL1 = 1
        End If

! DEPENDS ON: exppxi
        GRID_TYPE1=EXPPXI(FIELD_MODEL1,FIELD_SECT1,FIELD_ITEM1,         &
     &                          ppx_grid_type,                          &
#include "argppx.h"
     &                          ICODE,CMESSAGE)
          IF ((GRID_TYPE1 <  0).OR.(GRID_TYPE1 >  100)) THEN
            GRID_TYPE1=1
!           WRITE(6,*)'COMPARE: CANNOT GET GRID_TYPE1 INFO'
!           WRITE(6,*)'         FROM STASHMASTER FOR '
!           WRITE(6,*)'         SECTION ',FIELD_SECT1,
!    &                ' ITEM ',FIELD_ITEM1
!           WRITE(6,*)'         GRID_TYPE1 SET TO 1 (NORMAL GRID)'
          ENDIF
        FIELD_ITEM2=MOD(LOOKUP2(42,M),1000)
        FIELD_SECT2=(LOOKUP2(42,M)-FIELD_ITEM2)/1000
        FIELD_MODEL2=LOOKUP2(45,M)

        If (FIELD_MODEL2 == 10) Then
          FIELD_MODEL2 = 1
        End If

! DEPENDS ON: exppxi
        GRID_TYPE2=EXPPXI(FIELD_MODEL2,FIELD_SECT2,FIELD_ITEM2,         &
     &                          ppx_grid_type,                          &
#include "argppx.h"
     &                          ICODE,CMESSAGE)
          IF ((GRID_TYPE2 <  0).OR.(GRID_TYPE2 >  100)) THEN
            GRID_TYPE2=1
!           WRITE(6,*)'COMPARE: CANNOT GET GRID_TYPE2 INFO'
!           WRITE(6,*)'         FROM STASHMASTER FOR '
!           WRITE(6,*)'         SECTION ',FIELD_SECT2,
!    &                ' ITEM ',FIELD_ITEM2
!           WRITE(6,*)'         GRID_TYPE2 SET TO 1 (NORMAL GRID)'
          ENDIF

        PACK_CODE1 = MOD(LOOKUP1(21,I),10)
        PACK_CODE2 = MOD(LOOKUP2(21,M),10)

       lblrec_1=lookup1(15, i)
       lblrec_2=lookup2(15, m)

        if (((pack_code1 == 1 .or. pack_code2 == 1) .or.                &
     &       (pack_code1 == 4 .or. pack_code2 == 4)) .and.              &
     &  expand /= 1) then

        ELSEIF (PACK_CODE1 == 3 .OR. PACK_CODE2 == 3) THEN

          WRITE(6,*)                                                    &
     &    'Field No ',I,' not compared. GRIB data not supported.'

        ELSE

!         Since, for LBC data, header definition for pre Vn5.0
!         and Vn5.2- is different, a check is to be done.
          IF ((GRID_TYPE1 == ppx_atm_lbc_theta).OR.                     &
     &        (GRID_TYPE1 == ppx_atm_lbc_u).or.                         &
     &        (GRID_TYPE1 == ppx_atm_lbc_v).or.                         &
     &        (GRID_TYPE1 == ppx_atm_lbc_orog).or.                      &
     &        (GRID_TYPE1 == ppx_ocn_lbc_theta).or.                     &
     &        (GRID_TYPE1 == ppx_ocn_lbc_u)) then

            IF ((FIXHD1(12) >= 502).AND.                                &
     &          (FIXHD2(12) >= 502)) THEN
              IF ((LOOKUP1(41,i)/10000) /=                              &
     &            (LOOKUP2(41,i)/10000)) THEN
                write(6,*)'COMPARE: Fields with different rims',        &
     &                    '         not yet supported'
! DEPENDS ON: ereport
                CALL EREPORT('COMPARE', 1014,                           &
     &           'Fields with different rims not yet supported')

              ENDIF

              RIMWIDTHA=LOOKUP1(41,i)/10000
              haloY1=MOD(LOOKUP1(41,i),10000)
              haloY1=haloY1/100
              haloX1=MOD(haloY1,100)

              haloY2=MOD(LOOKUP2(41,i),10000)
              haloY2=haloY2/100
              haloX2=MOD(haloY2,100)

            ENDIF
          ENDIF
      IF((LOOKUP1(39,I) ==  1 .AND. LOOKUP2(39,M) ==  1).OR.            &
     &   (LOOKUP1(39,I) == -1 .AND. LOOKUP2(39,M) == -1))THEN
! This is a REAL field


!     field is atmos/ocean LB field
!     Decomposition is compulsory as each field may have
!     different halo size
      IF ((GRID_TYPE1 == ppx_atm_lbc_theta).OR.                         &
     &    (GRID_TYPE1 == ppx_atm_lbc_u).OR.                             &
     &    (GRID_TYPE1 == ppx_atm_lbc_v).OR.                             &
     &    (GRID_TYPE1 == ppx_atm_lbc_orog).OR.                          &
     &    (GRID_TYPE1 == ppx_ocn_lbc_theta).OR.                         &
     &    (GRID_TYPE1 == ppx_ocn_lbc_u)) THEN

        CURRENT_DECOMP_TYPE=-1

! DEPENDS ON: decompose_smexe
        CALL DECOMPOSE_SMEXE(ROW_LENGTH1,P_ROWS1,                       &
     &                       haloX1, haloY1,                            &
     &                       TOT_LEVELS)
! DEPENDS ON: change_decomposition
        CALL CHANGE_DECOMPOSITION(decomp_smexe,ICODE)
! DEPENDS ON: dervsize
        CALL DERVSIZE(ICODE,CMESSAGE)

      ENDIF


! DEPENDS ON: readflds
        CALL READFLDS(NFTIN1,1,I,LOOKUP1,LEN1_LOOKUP1,                  &
     &               R_D1,MAX_FIELD_SIZE1,FIXHD1,                       &
#include "argppx.h"
     &               expand,icode,cmessage)
! DEPENDS ON: abort_io
        IF(ICODE /= 0)CALL ABORT_IO('COMPARE',CMESSAGE,ICODE,NFTIN1)


!       field is atmos/ocean LB field
        IF ((GRID_TYPE2 == ppx_atm_lbc_theta).OR.                       &
     &      (GRID_TYPE2 == ppx_atm_lbc_u).or.                           &
     &      (GRID_TYPE2 == ppx_atm_lbc_v).or.                           &
     &      (GRID_TYPE2 == ppx_atm_lbc_orog).or.                        &
     &      (GRID_TYPE2 == ppx_ocn_lbc_theta).or.                       &
     &      (GRID_TYPE2 == ppx_ocn_lbc_u)) then

          CURRENT_DECOMP_TYPE=-1

! DEPENDS ON: decompose_smexe
          CALL DECOMPOSE_SMEXE(ROW_LENGTH2,P_ROWS2,                     &
     &                         haloX2, haloY2,                          &
     &                         TOT_LEVELS)
! DEPENDS ON: change_decomposition
          CALL CHANGE_DECOMPOSITION(decomp_smexe,ICODE)
! DEPENDS ON: dervsize
          CALL DERVSIZE(ICODE,CMESSAGE)
        ENDIF

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN2,1,M,LOOKUP2,LEN1_LOOKUP2,                  &
     &                R_D2,MAX_FIELD_SIZE2,FIXHD2,                      &
#include "argppx.h"
     &               expand,icode,cmessage)
! DEPENDS ON: abort_io
        IF(ICODE /= 0)CALL ABORT_IO('COMPARE',CMESSAGE,ICODE,NFTIN1)

      ELSE IF((LOOKUP1(39,I) ==  2 .AND. LOOKUP2(39,M) ==  2).OR.       &
     &        (LOOKUP1(39,I) == -2 .AND. LOOKUP2(39,M) == -2))THEN
! This is an INTEGER field

!     field is atmos/ocean LB field
      IF ((GRID_TYPE1 == ppx_atm_lbc_theta).OR.                         &
     &    (GRID_TYPE1 == ppx_atm_lbc_u).OR.                             &
     &    (GRID_TYPE1 == ppx_atm_lbc_v).OR.                             &
     &    (GRID_TYPE1 == ppx_atm_lbc_orog).OR.                          &
     &    (GRID_TYPE1 == ppx_ocn_lbc_theta).OR.                         &
     &    (GRID_TYPE1 == ppx_ocn_lbc_u)) THEN

        CURRENT_DECOMP_TYPE=-1

! DEPENDS ON: decompose_smexe
        CALL DECOMPOSE_SMEXE(ROW_LENGTH1,P_ROWS1,                       &
     &                       haloX1, haloY1,                            &
     &                       TOT_LEVELS)
! DEPENDS ON: change_decomposition
        CALL CHANGE_DECOMPOSITION(4,ICODE)
! DEPENDS ON: change_decomposition
        CALL CHANGE_DECOMPOSITION(decomp_smexe,ICODE)
! DEPENDS ON: dervsize
        CALL DERVSIZE(ICODE,CMESSAGE)

      ENDIF


! DEPENDS ON: readflds
        CALL READFLDS(NFTIN1,1,I,LOOKUP1,LEN1_LOOKUP1,                  &
     &                I_D1,MAX_FIELD_SIZE1,FIXHD1,                      &
#include "argppx.h"
     &               expand,icode,cmessage)
! DEPENDS ON: abort_io
        IF(ICODE /= 0)CALL ABORT_IO('COMPARE',CMESSAGE,ICODE,NFTIN1)

!       field is atmos/ocean LB field
        IF ((GRID_TYPE2 == ppx_atm_lbc_theta).OR.                       &
     &      (GRID_TYPE2 == ppx_atm_lbc_u).or.                           &
     &      (GRID_TYPE2 == ppx_atm_lbc_v).or.                           &
     &      (GRID_TYPE2 == ppx_atm_lbc_orog).or.                        &
     &      (GRID_TYPE2 == ppx_ocn_lbc_theta).or.                       &
     &      (GRID_TYPE2 == ppx_ocn_lbc_u)) then

          CURRENT_DECOMP_TYPE=-1

! DEPENDS ON: decompose_smexe
          CALL DECOMPOSE_SMEXE(ROW_LENGTH2,P_ROWS2,                     &
     &                         haloX2, haloY2,                          &
     &                         TOT_LEVELS)
! DEPENDS ON: change_decomposition
          CALL CHANGE_DECOMPOSITION(4,ICODE)
! DEPENDS ON: dervsize
          CALL DERVSIZE(ICODE,CMESSAGE)
        ENDIF

! DEPENDS ON: readflds
       CALL READFLDS(NFTIN2,1,M,LOOKUP2,LEN1_LOOKUP2,                   &
     &                I_D2,MAX_FIELD_SIZE2,FIXHD2,                      &
#include "argppx.h"
     &               expand,icode,cmessage)
! DEPENDS ON: abort_io
        IF(ICODE /= 0)CALL ABORT_IO('COMPARE',CMESSAGE,ICODE,NFTIN1)

      ELSE IF((LOOKUP1(39,I) ==  3 .AND. LOOKUP2(39,M) ==  3).OR.       &
     &        (LOOKUP1(39,I) == -3 .AND. LOOKUP2(39,M) == -3))THEN
! This is an LOGICAL field

!     field is atmos/ocean LB field
      IF ((GRID_TYPE1 == ppx_atm_lbc_theta).OR.                         &
     &    (GRID_TYPE1 == ppx_atm_lbc_u).OR.                             &
     &    (GRID_TYPE1 == ppx_atm_lbc_v).OR.                             &
     &    (GRID_TYPE1 == ppx_atm_lbc_orog).OR.                          &
     &    (GRID_TYPE1 == ppx_ocn_lbc_theta).OR.                         &
     &    (GRID_TYPE1 == ppx_ocn_lbc_u)) THEN

        CURRENT_DECOMP_TYPE=-1

! DEPENDS ON: decompose_smexe
        CALL DECOMPOSE_SMEXE(ROW_LENGTH1,P_ROWS1,                       &
     &                       haloX1, haloY1,                            &
     &                       TOT_LEVELS)
! DEPENDS ON: change_decomposition
        CALL CHANGE_DECOMPOSITION(decomp_smexe,ICODE)
! DEPENDS ON: dervsize
        CALL DERVSIZE(ICODE,CMESSAGE)
      ENDIF


! DEPENDS ON: readflds
        CALL READFLDS(NFTIN1,1,I,LOOKUP1,LEN1_LOOKUP1,                  &
     &                L_D1,MAX_FIELD_SIZE1,FIXHD1,                      &
#include "argppx.h"
     &               expand,icode,cmessage)
! DEPENDS ON: abort_io
        IF(ICODE /= 0)CALL ABORT_IO('COMPARE',CMESSAGE,ICODE,NFTIN1)

!       field is atmos/ocean LB field
        IF ((GRID_TYPE2 == ppx_atm_lbc_theta).OR.                       &
     &      (GRID_TYPE2 == ppx_atm_lbc_u).or.                           &
     &      (GRID_TYPE2 == ppx_atm_lbc_v).or.                           &
     &      (GRID_TYPE2 == ppx_atm_lbc_orog).or.                        &
     &      (GRID_TYPE2 == ppx_ocn_lbc_theta).or.                       &
     &      (GRID_TYPE2 == ppx_ocn_lbc_u)) then

          CURRENT_DECOMP_TYPE=-1

! DEPENDS ON: decompose_smexe
          CALL DECOMPOSE_SMEXE(ROW_LENGTH2,P_ROWS2,                     &
     &                         haloX2, haloY2,                          &
     &                         TOT_LEVELS)
! DEPENDS ON: change_decomposition
          CALL CHANGE_DECOMPOSITION(decomp_smexe,ICODE)
! DEPENDS ON: dervsize
          CALL DERVSIZE(ICODE,CMESSAGE)
        ENDIF

! DEPENDS ON: readflds
       CALL READFLDS(NFTIN2,1,M,LOOKUP2,LEN1_LOOKUP2,                   &
     &                L_D2,MAX_FIELD_SIZE2,FIXHD2,                      &
#include "argppx.h"
     &               expand,icode,cmessage)
! DEPENDS ON: abort_io
        IF(ICODE /= 0)CALL ABORT_IO('COMPARE',CMESSAGE,ICODE,NFTIN1)

      ELSE
! This is an unrecognized field
          WRITE(6,*)                                                    &
     &     'Field No ',I,' not compared. Unrecognized type.'

      ENDIF
       if (((pack_code1 == 1 .or. pack_code2 == 1).or.                  &
     &      (pack_code1 == 4 .or. pack_code2 == 4)).and.                &
     &  expand == 1) then

         LEN_FIELD=LOOKUP1(15,I)
       endif
       lookup1(15, i)=lblrec_1
       lookup2(15, M)=lblrec_2


        WRITE(6,*)' '
        write(6,'(/''Field '',i5,'' : Stash Code '',i5,                 &
     &   '' : '',a,'' : Level '',i4/)')                                 &
     &   i, lookup1(42,i), phrase, lookup1(33,i)
        write(10,'(/''Field '',i5,'' : Stash Code '',i5,                &
     &   '' : '',a)') i, lookup1(42,i), phrase
        IF(GRID_TYPE1 >  100)THEN
          WRITE(10,*)'Grid type = information not available'
        ELSE
          WRITE(10,*)'Grid type = ',GRID_TYPE1
        ENDIF

        RMS_F1=0.0
        RMS_F2=0.0
        RMS_DIFF=0.0
        K=0
#if defined(T3E)
        jrc_nan=0
#endif
!       Real
        IF (LOOKUP1(39,I) == 1 .AND. LOOKUP2(39,M) == 1) THEN
          MAX_DIFF=0.
          DO J=1,LEN_FIELD
            DIFF(J)='.'
#if defined(T3E)
            IF(XOR(R_D1(J),R_D2(J)) /= 0) THEN
#else
            IF(R_D1(J) /= R_D2(J))THEN
#endif
              k=k+1
              if(k <= 10) then
                write(6,'(a,i6,2(e25.15,'' ('',z16,'')''))')            &
     &            'ITEM=',j,r_d1(j),r_d1(j),r_d2(j),r_d2(j)
              endif
#if defined(T3E)
              if((xor(and(r_d1(j),jrc_mask),jrc_mask) /= 0) .and.       &
     &         (xor(and(r_d2(j),jrc_mask),jrc_mask) /= 0)   .and.       &
     &         (and(r_d1(j),deb_mask) /= 0.or.                          &
     &          xor(r_d1(j),0) == 0).and.                               &
     &         (and(r_d2(j),deb_mask) /= 0.or.                          &
     &          xor(r_d2(j),0) == 0))then
#endif
                RD1=R_D1(J)
                RD2=R_D2(J)
                MAX_DIFF=MAX(MAX_DIFF,ABS(RD1-RD2))
                IF(MAX_DIFF == ABS(RD1-RD2)) MAX_J = J
                if(rd1 == 0.) then
                  if(rd2 == 0.) then
                    diff_per=0.
                  else
                    diff_per=(abs(rd1-rd2)/abs(rd2))*100
                  endif
                else
                  diff_per=(abs(rd1-rd2)/abs(rd1))*100
                endif
              IF (DIFF_PER  >   10.0) DIFF(J)="#"
              IF (DIFF_PER  <   10.0) DIFF(J)="X"
              IF (DIFF_PER  <   1.0) DIFF(J)="O"
              IF (DIFF_PER  <   0.1) DIFF(J)="o"
              IF (DIFF_PER  <   0.01) DIFF(J)=":"
                RMS_F1=RMS_F1+(R_D1(J)*R_D1(J))
                RMS_F2=RMS_F2+(R_D2(J)*R_D2(J))
                RMS_DIFF=RMS_DIFF+(R_D1(J)-R_D2(J))*(R_D1(J)-R_D2(J))
#if defined(T3E)
              else
                jrc_nan=jrc_nan+1
              endif
#endif
            else
#if defined(T3E)
              if((xor(and(r_d1(j),jrc_mask),jrc_mask) /= 0) .and.       &
     &         (xor(and(r_d2(j),jrc_mask),jrc_mask) /= 0)   .and.       &
     &         (and(r_d1(j),deb_mask) /= 0.or.                          &
     &       xor(r_d1(j),0) == 0).and.                                  &
     &         (and(r_d2(j),deb_mask) /= 0.or.                          &
     &       xor(r_d2(j),0) == 0))then
#endif
                RMS_F1=RMS_F1+(R_D1(J)*R_D1(J))
                RMS_F2=RMS_F2+(R_D2(J)*R_D2(J))
                RMS_DIFF=RMS_DIFF+(R_D1(J)-R_D2(J))*(R_D1(J)-R_D2(J))
#if defined(T3E)
              else
                jrc_nan=jrc_nan+1
              endif
#endif
            endif
          ENDDO
          IF (K /= 0) THEN
            WRITE(6,*)'NUMBER OF DIFFERENT VALUES = ',K
            WRITE(6,*)'MAXIMUM DIFFERENCE= ',MAX_DIFF,' AT PT. ',MAX_J
            RMS_F1=SQRT(RMS_F1/LEN_FIELD)
            RMS_F2=SQRT(RMS_F2/LEN_FIELD)
            RMS_DIFF=SQRT(RMS_DIFF/LEN_FIELD)
            WRITE(6,*) 'RMS FIELD1 : ',RMS_F1
            WRITE(6,*) 'RMS FIELD2 : ',RMS_F2
            DIFF_PER=ABS(RMS_F1-RMS_F2)
            WRITE(6,*) 'Difference: ',DIFF_PER                          &
     &                ,' RMS_difference: ',RMS_DIFF
            rd1=diff_per
            if(rms_f1 /= 0.) then
              diff_per=(diff_per/rms_f1)*100
              write(6,'(''Field '',i5,'' has a Difference between'',    &
     &         '' the RMS Values of '',e10.5,'' which is '',f10.3,      &
     &         '' Percent of Field 1, whose RMS Value is '',e10.5)')    &
     &         i, rd1, diff_per, rms_f1
              write(6,*) 'Difference as % of RMS FIELD1= ',DIFF_PER
            else if(rms_f2 /= 0.) then
              diff_per=(diff_per/rms_f2)*100
              write(6,'(''Field '',i5,'' has a Difference between'',    &
     &         '' the RMS Values of '',e10.5,'' which is '',f10.3,      &
     &         '' Percent of Field 2, whose RMS Value is '',e10.5)')    &
     &         i, rd1, diff_per, rms_f2
              write(6,*) 'Difference as % of RMS FIELD2= ',DIFF_PER
#if defined(T3E)
            else
              if(jrc_nan == 0) write(6,'(''Field '',i5,                 &
     &         '' - the Fields in Both Files have RMS Values of Zero''  &
     &         )') i
#endif
            endif
!
#if defined(T3E)
            if (diff_per  >   5 .or. jrc_nan /= 0) THEN
#else
            if (diff_per  >   5) THEN
#endif
              WRITE(6,*)
              WRITE(6,*)
#if defined(T3E)
              if(jrc_nan /= 0) then
                write(6,*) '********** NaN Values Detected **********', &
     &           '**'
              endif
#endif
              WRITE(6,*) '************** WARNING ********************'
              WRITE(6,*) '***** LARGE DIFFERENCE ENCOUNTERED ********'
              WRITE(6,*) '*******************************************'
              WRITE(6,*)
            ENDIF

            KEY='# d>10%  ;  X 10%>d>1%  ;  O 1%>d>0.1%  ; '//          &
     &          'o 0.1%>d>0.01%  ;  : d<0.01%  ;  . d=0%  ; '//         &
     &          '~ no data (sea-points)'

! In some cases, such as non-standard fields, grid type is not
! specified. Random number result.  In these cases, diff map will
! be printed if field size is (73*96).
! Grid-type 21 will have its new converted diff map of size
! (nrow*ncol) which is the same as land-mask it's been applied to
! Grid-type other then 21 have its map dimension lookup(18)
! *lookup(19)
            IF(((GRID_TYPE1 <  60).AND.                                 &
     &         (GRID_TYPE1 /= 31).AND.                                  &
     &         (GRID_TYPE1 /= 32).AND.                                  &
     &         (GRID_TYPE1 /= 47).AND.                                  &
     &         (GRID_TYPE1 /= 51))                                      &
     &      .OR.                                                        &
     &        ((GRID_TYPE1 >  100).AND.                                 &
     &         (LOOKUP1(18,I) >  0).AND.(LOOKUP1(18,I) <  300).AND.     &
     &         (LOOKUP1(19,I) >  0).AND.(LOOKUP1(19,I) <  300)))THEN
! Only certain grid types are suitable for difference maps
! For real data
              IF((GRID_TYPE1 == 21).AND.(GRID_TYPE2 == 21)) THEN
! DEPENDS ON: print_dif_map
                CALL PRINT_DIF_MAP(DIFF,P_ROWS1,ROW_LENGTH1,KEY,        &
     &            GRID_TYPE1,GRID_TYPE2)
              ELSE IF(GRID_TYPE1 == GRID_TYPE2.and.                     &
     &            ((GRID_TYPE2 == ppx_atm_lbc_theta).OR.                &
     &            (GRID_TYPE2 == ppx_atm_lbc_u).or.                     &
     &            (GRID_TYPE2 == ppx_atm_lbc_v).or.                     &
     &            (GRID_TYPE2 == ppx_atm_lbc_orog).or.                  &
     &            (GRID_TYPE2 == ppx_ocn_lbc_theta).or.                 &
     &            (GRID_TYPE2 == ppx_ocn_lbc_u))) then

                write(10,*)'LBC field: use arbitrary grid dimension'
                write(10,*)'80 columns wide to display differences'
! DEPENDS ON: print_dif_map
               CALL PRINT_DIF_MAP(DIFF,LEN_FIELD/80,80,KEY,             &
     &           GRID_TYPE1,GRID_TYPE2)
              ELSE
! DEPENDS ON: print_dif_map
               CALL PRINT_DIF_MAP(DIFF,LOOKUP1(18,I),LOOKUP1(19,I),KEY, &
     &           GRID_TYPE1,GRID_TYPE2)
              ENDIF
            ELSE
              write(6,*) 'Difference map not printed'
              write(6,*) 'Grid Type not suitable for difference maps'
              write(6,*) 'Grid Type = ',LOOKUP1(16,I)
            ENDIF

          ELSE
            WRITE(6,*)'OK'
            WRITE(10,*)'OK'
          ENDIF
!       Integer
        ELSE IF (LOOKUP1(39,I) == 2 .AND. LOOKUP2(39,M) == 2) THEN
          IMAX_DIFF=0
          DO J=1,LEN_FIELD
            DIFF(J)='.'
#if defined(T3E)
            IF (XOR(I_D1(J),I_D2(J)) /= 0) THEN
#else
            IF(I_D1(J) /= I_D2(J))THEN
#endif
              K=K+1
              if (k <= n_diff) write(6,*)'item=',j,i_d1(j),i_d2(j)
#if defined(T3E)
              if((xor(and(r_d1(j),jrc_mask),jrc_mask) /= 0) .and.       &
     &         (xor(and(r_d2(j),jrc_mask),jrc_mask) /= 0)) then
#endif
                ID1=I_D1(J)
                ID2=I_D2(J)
                IMAX_DIFF=MAX(IMAX_DIFF,ABS(ID1-ID2))

                IF (ID1  ==  0) THEN
                  IF (ID1  ==  ID2) THEN
                    DIFF_PER=0.0
                  ELSE
                    DIFF_PER=100.0
                  ENDIF
                ELSE
                  DIFF_PER=(REAL(ABS(ID1-ID2))/REAL(ABS(ID1)))*100.0
                ENDIF
              IF (DIFF_PER  >   10.0) DIFF(J)="#"
              IF (DIFF_PER  <   10.0) DIFF(J)="X"
              IF (DIFF_PER  <   1.0) DIFF(J)="O"
              IF (DIFF_PER  <   0.1) DIFF(J)="o"
              IF (DIFF_PER  <   0.01) DIFF(J)=":"
#if defined(T3E)
              else
                jrc_nan=jrc_nan+1
              endif
#endif
            ENDIF
          ENDDO
          IF (K /= 0) THEN
            write(6,'(''Field '',i5,'' has '',i5,                       &
     &       '' INTEGER Differences'',                                  &
     &       '' with a Maximum Difference of '',i20)') i, k, imax_diff
#if defined(T3E)
            WRITE(6,*)'NUMBER OF DIFFERENT VALUES = ',K                 &
     &      ,' (',jrc_nan,') NaN Values Detected)'
#else
            WRITE(6,*)'NUMBER OF DIFFERENT VALUES = ',K
#endif
            WRITE(6,*)'MAXIMUM DIFFERENCE= ',IMAX_DIFF
            KEY='# d>10%  ;  X 10%>d>1%  ;  O 1%>d>0.1%  ; '//          &
     &          'o 0.1%>d>0.01%  ;  : d<0.01%  ;  . d=0%  ; '//         &
     &          '~ no data (sea-points)'

! In some cases, such as non-standard fields, grid type is not
! specified. Random number result.  In these cases, diff map will
! be printed if field size is (73*96).
! Grid-type 21 will have its new converted diff map of size
! (nrow*ncol) which is the same as land-mask it's been applied to
! Grid-type other then 21 have its map dimension lookup(18)
! *lookup(19)
            IF(((GRID_TYPE1 <  60).AND.                                 &
     &         (GRID_TYPE1 /= 31).AND.                                  &
     &         (GRID_TYPE1 /= 32).AND.                                  &
     &         (GRID_TYPE1 /= 47).AND.                                  &
     &         (GRID_TYPE1 /= 51))                                      &
     &      .OR.                                                        &
     &        ((GRID_TYPE1 >  100).AND.                                 &
     &         (LOOKUP1(18,I) >  0).AND.(LOOKUP1(18,I) <  300).AND.     &
     &         (LOOKUP1(19,I) >  0).AND.(LOOKUP1(19,I) <  300)))THEN
! Only certain grid types are suitable for difference maps
! For integer data
              IF((GRID_TYPE1 == 21).AND.(GRID_TYPE2 == 21)) THEN
! DEPENDS ON: print_dif_map
                CALL PRINT_DIF_MAP(DIFF,P_ROWS1,ROW_LENGTH1,KEY,        &
     &            GRID_TYPE1,GRID_TYPE2)
              ELSE
! DEPENDS ON: print_dif_map
                CALL PRINT_DIF_MAP(DIFF,LOOKUP1(18,I),LOOKUP1(19,I),KEY,&
     &            GRID_TYPE1,GRID_TYPE2)
              ENDIF
            ENDIF

          ELSE
           write(6,'(''Field '',i5,                                     &
     &       '' has '',i5,'' INTEGER Differences'')') i, k
            write (6,*) 'OK'
            WRITE(6,*)'OK'
            WRITE(10,*)'OK'
          ENDIF
!       Logical
        ELSE IF (LOOKUP1(39,I) == 3 .AND. LOOKUP2(39,M) == 3) THEN
          DO J=1,LEN_FIELD
            DIFF(J)='.'
            IF (L_D1(J).NEQV.L_D2(J)) THEN
              K=K+1
              LD1=L_D1(J)
              LD2=L_D2(J)
              IF (K <= N_DIFF) WRITE(6,*)'ITEM=',J,LD1,LD2
              DIFF(J)="#"
            ENDIF
          ENDDO
          IF (K /= 0) THEN
            write(6,'(''Field '',i5,'' has '',i5,                       &
     &       '' LOGICAL Differences'')') i, k
            WRITE(6,*)'NUMBER OF DIFFERENT VALUES = ',K
            KEY='# Different values  ;  . identical  ; '//              &
     &          '~ no data (sea-points)'

! In some cases, such as non-standard fields, grid type is not
! specified. Random number result.  In these cases, diff map will
! be printed if field size is (73*96).
! Grid-type 21 will have its new converted diff map of size
! (nrow*ncol) which is the same as land-mask it's been applied to
! Grid-type other then 21 have its map dimension lookup(18)
! *lookup(19)
            IF(((GRID_TYPE1 <  60).AND.                                 &
     &         (GRID_TYPE1 /= 31).AND.                                  &
     &         (GRID_TYPE1 /= 32).AND.                                  &
     &         (GRID_TYPE1 /= 47).AND.                                  &
     &         (GRID_TYPE1 /= 51))                                      &
     &      .OR.                                                        &
     &        ((GRID_TYPE1 >  100).AND.                                 &
     &         (LOOKUP1(18,I) >  0).AND.(LOOKUP1(18,I) <  300).AND.     &
     &         (LOOKUP1(19,I) >  0).AND.(LOOKUP1(19,I) <  300)))THEN
! Only certain grid types are suitable for difference maps
! For logical data
              IF((GRID_TYPE1 == 21).AND.(GRID_TYPE2 == 21)) THEN
! DEPENDS ON: print_dif_map
                CALL PRINT_DIF_MAP(DIFF,P_ROWS1,ROW_LENGTH1,KEY,        &
     &            GRID_TYPE1,GRID_TYPE2)
              ELSE
! DEPENDS ON: print_dif_map
                CALL PRINT_DIF_MAP(DIFF,LOOKUP1(18,I),LOOKUP1(19,I),KEY,&
     &            GRID_TYPE1,GRID_TYPE2)
              ENDIF
            ENDIF
          ELSE
           write(6,'(''Field '',i5,'' has '',i5,                        &
     &       '' LOGICAL Differences'')') i, k
            write (6,*) 'OK'
            WRITE(6,*)'OK'
            WRITE(10,*)'OK'
          ENDIF
!       Real Timeseries
        ELSE IF (LOOKUP1(39,I) == -1 .AND. LOOKUP2(39,M) == -1) THEN
          MAX_DIFF=0.
          DO J=1,LEN_FIELD
#if defined(T3E)
            IF(XOR(R_D1(J),R_D2(J)) /= 0) THEN
#else
            IF(R_D1(J) /= R_D2(J))THEN
#endif
              MAX_DIFF=AMAX1(MAX_DIFF,ABS(R_D1(J)-R_D2(J)))
              K=K+1
              IF(K <= 10)WRITE(6,*)'ITEM=',J,R_D1(J),R_D2(J)
            ENDIF
          ENDDO
          IF (K /= 0) THEN
            WRITE(6,*)'NUMBER OF DIFFERENT VALUES = ',K
            WRITE(6,*)'MAXIMUM DIFFERENCE= ',MAX_DIFF
          ELSE
            WRITE(6,*)'OK'
          ENDIF
!       Integer Timeseries
        ELSE IF (LOOKUP1(39,I) == -2 .AND. LOOKUP2(39,M) == -2) THEN
          IMAX_DIFF=0
          DO J=1,LEN_FIELD
            IF (I_D1(J) /= I_D2(J)) THEN
              K=K+1
              ID1=I_D1(J)
              ID2=I_D2(J)
              IMAX_DIFF=MAX(IMAX_DIFF,IABS(ID1-ID2))
              IF (K <= N_DIFF) WRITE(6,*)'ITEM=',J,ID1,ID2
            ENDIF
          ENDDO
          IF (K /= 0) THEN
            WRITE(6,*)'NUMBER OF DIFFERENT VALUES = ',K
            WRITE(6,*)'MAXIMUM DIFFERENCE= ',IMAX_DIFF
          ELSE
            WRITE(6,*)'OK'
          ENDIF
!       Logical Timeseries
        ELSE IF (LOOKUP1(39,I) == -3 .AND. LOOKUP2(39,M) == -3) THEN
          DO J=1,LEN_FIELD
            IF (L_D1(J).NEQV.L_D2(J)) THEN
              K=K+1
              LD1=L_D1(J)
              LD2=L_D2(J)
              IF (K <= N_DIFF) WRITE(6,*)'ITEM=',J,LD1,LD2
            ENDIF
          ENDDO
          IF (K /= 0) THEN
            WRITE(6,*)'NUMBER OF DIFFERENT VALUES = ',K
          ELSE
            WRITE(6,*)'OK'
          ENDIF
        ELSE
          WRITE(6,*)                                                    &
     &    'Field No ',I,' not compared. Different Data Type Numbers ?'
        ENDIF
        WRITE(6,*)' '
        IF(K /= 0)THEN
          NDIFFER(I)=K
          N=N+1
        END IF
        L=L+K
        END IF

      ENDIF
      ENDIF
      ENDDO     !End loop over number of fields

! Output remainder of summary information
      WRITE(8,*) 'DATA FIELDS:                ',                        &
     &             'Number of fields with differences = ',N
      DO I = 1,NUMREC1   ! Begin loop over number of fields in file1
        IF(LOOKUP1(1,I) /= -99)THEN
          S_ITEM_CODE=MOD(LOOKUP1(42,I),1000)
          SECTION=(LOOKUP1(42,I)-S_ITEM_CODE)/1000
          IF(FIXHD2(12) >= 305)THEN
            MODEL=LOOKUP1(45,I)
          ELSEIF(S_ITEM_CODE <= 100.OR.                                 &
     &          (S_ITEM_CODE >= 200.AND.S_ITEM_CODE <= 205))THEN
            MODEL=1
          ELSEIF((S_ITEM_CODE >  100.AND.S_ITEM_CODE <= 176).OR.        &
     &           (S_ITEM_CODE >= 180.AND.S_ITEM_CODE <= 200))THEN
            MODEL=2
          ELSEIF((S_ITEM_CODE >= 177.AND.S_ITEM_CODE <= 179).OR.        &
     &           (S_ITEM_CODE >= 210.AND.S_ITEM_CODE <= 212))THEN
            MODEL=3
          END IF

!         All diagnostics under model code of 10 are in section 20
!         of Atmos StashMaster file.
          If (MODEL == 10) Then
            MODEL = 1
          End If

! DEPENDS ON: exppxc
          PHRASE=EXPPXC(MODEL,SECTION,S_ITEM_CODE,                      &
#include "argppx.h"
     &                  ICODE,CMESSAGE)
          IF(ICODE /= 0)THEN
             PHRASE='NON-STANDARD FIELD'
          END IF
          IF(LMISSING1(I))THEN
            WRITE(8,'(/''Field '',i5,'' : Stash Code '',i5,'' : '',a,   &
     &                 '' : No equivalent in file2'')')                 &
     &        I,LOOKUP1(42,I),PHRASE
          ELSE IF(NDIFFER(I) /= 0)THEN
            WRITE(8,'(/''Field '',I5,'' : Stash Code '',I5,             &
     &                 '' : '',A,'' : Number of differences = '',I8)')  &
     &        I, LOOKUP1(42,I), PHRASE, NDIFFER(I)
          END IF
        END IF
      END DO
      DO I = 1,NUMREC2   ! Begin loop over number of fields in file2
        IF(LMISSING2(I).AND.LOOKUP2(1,I) /= -99)THEN
          S_ITEM_CODE=MOD(LOOKUP2(42,I),1000)
          SECTION=(LOOKUP2(42,I)-S_ITEM_CODE)/1000
          IF(FIXHD2(12) >= 305)THEN
            MODEL=LOOKUP2(45,I)
          ELSEIF(S_ITEM_CODE <= 100.OR.                                 &
     &          (S_ITEM_CODE >= 200.AND.S_ITEM_CODE <= 205))THEN
            MODEL=1
          ELSEIF((S_ITEM_CODE >  100.AND.S_ITEM_CODE <= 176).OR.        &
     &           (S_ITEM_CODE >= 180.AND.S_ITEM_CODE <= 200))THEN
            MODEL=2
          ELSEIF((S_ITEM_CODE >= 177.AND.S_ITEM_CODE <= 179).OR.        &
     &           (S_ITEM_CODE >= 210.AND.S_ITEM_CODE <= 212))THEN
            MODEL=3
          END IF

          If (MODEL == 10) Then
            MODEL = 1
          End If

! DEPENDS ON: exppxc
          PHRASE=EXPPXC(MODEL,SECTION,S_ITEM_CODE,                      &
#include "argppx.h"
     &                  ICODE,CMESSAGE)
          IF(ICODE /= 0)THEN
             PHRASE='NON-STANDARD FIELD'
          END IF
          WRITE(8,'(/''Field '',i5,'' : Stash Code '',i5,'' : '',a,     &
     &               '' : No equivalent in file1'')')                   &
     &      I,LOOKUP2(42,I),PHRASE
        ENDIF
      END DO
      CLOSE(10)
      IF(L == 0)THEN
        WRITE(8,*)' files compare (ignoring Fixed Length Header)'
      ELSE
        WRITE(8,*)' files DO NOT compare'
      ENDIF
      WRITE(7,*)' '
      CLOSE(7)

 9999 CONTINUE
      RETURN
      END SUBROUTINE COMPARE

#endif
