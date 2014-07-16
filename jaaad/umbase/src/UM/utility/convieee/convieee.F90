#if defined(CONVIEEE)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL   PROGRAM MAIN_CONVIEEE and SUBROUTINE CONVIEEE ------------------
!LL
!LL  Written by A. Dickinson 05/05/92
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.2  06/04/93   Correct use of packing indicator as per vn2.8
!LL                   Author: A.Dickinson        Reviewer: P.Burton
!LL
!LL   3.2  06/05/93    Extend code to recognise PP type files
!LL                    Author: A. Dickinson    Reviewer: D. Richardson
!LL   3.3  08/12/93    Extra argument for READFLDS.
!LL                    Author: A. Dickinson    Reviewer: D. Richardson
!LL
!LL  3.3   31/10/93   Dimension of data array set to maximum value
!LL                   Author: A. Dickinson      Reviewer: P.Burton
!LL  3.4   11/10/94  Part of modset which makes sure that LOGICAL's are
!LL                  set correctly for IEEE machines covered by CONVIEEE
!LL                  Necessary to port model to a T3D.
!LL                  Author D.M. Goddard
!LL   3.5  24/03/95    Changed OPEN to FILE_OPEN  P.Burton
!    4.0  21/12/95  Timeseries now catered for
!                   Author D.M. Goddard.
!     4.0  18/06/96   Changes to cope with changes in STASH addressing
!                     Author D.M. Goddard.
!     4.3   02/04/97  Remove surplus definition of GETARG D.M.Goddard
!     4.3  15/04/97   Extra argument for READFLDS to select 32-64 bit
!                     expansion routine EXPAND21 or C90_EXPAND21
!                     D.M.Goddard
!     4.3  06/05/97   Prevents program crashing if a number, which is
!                     unrepresentable in 32 bits, is present.
!                     Unrepresentable number is replaced by RMDI
!                     Author: D.M. Goddard
!LL  4.4   Oct. 1997 Changed error handling from routine HDPPXRF
!LL                  so only fatal (+ve) errors are handled.
!LL                                             Shaun de Witt
!     4.4  23/09/97   Produce correct well-formed 32-bit
!                     dumpfiles.
!                       Author:  Bob Carruthhers, Cray Research
!     4.4  17/07/97  Introduce conversion from ieee to Cray PVP
!                    numbers and reintroduce functionality for
!                    PVP machines
!                    Author: D.M. Goddard
!   4.4  24/10/97   Initialise ICODE as it is no longer
!                   initialised in HDPPXRF
!                   Author D.M. Goddard
!   4.5  01/04/98   Removed SETPOS32 subroutine as it is now available
!                   in C PORTIO2A.                           P.Burton
!     4.5  14/10/97   Sets correct most significant number
!                     in packing indicator.
!                     Either 2 for CRAY format or 3 for IEEE format
!                     Author D.M. Goddard
!     4.5  14/01/98   Conversion of data moved into new subroutines
!                     ATMOS_CONVIEEE and OCEAN_CONVIEEE depending on
!                     whether an atmosphere or ocean dataset is being
!                     processed. The ocean subroutine contains
!                     additional code to expand ocean compressed fields
!                     before conversion.
!                     Author D.M. Goddard
!     4.5  13/07/98   In boundary datasets the entire block of data for
!                     a given time is rounded up to a sector boundary
!                     in well formed datasets rather than individual
!                     fields. Subroutine set_dumpfile_address is
!                     skipped for boundary files and the addressing
!                     caluclated in subroutine CONVIEEE.
!                     Author D.M. Goddard
!     5.1  09/05/00   Updated to work with new version of READFLDS.
!                     Fixed bugs in code to handle WGDOS packed fields.
!                     Reduced memory usage for ocean datasets which do
!                     not contain compressed fields.
!                     Added call to InitPrintStatus.
!                     Removed some unused variables. D.P.Matthews
!LL  5.2  7/11/99  Enable run length encoding of ocean fieldsfiles to
!LL                compress the sequences of mdi values that represent
!LL                the land points. Ian Edmond
!    5.3  22/11/01  Enable MPP as the only option for
!                   small executables         E.Leung
!    5.4  02/05/02  Included extra-data recognition algorithm to
!                   allow extra-data to be processed     E.Leung
!    5.4  03/09/02  Arguments are added to DECOMPOSE_SMEXE
!                                                  E.Leung
!    5.5  07/04/03  Use INTHD 6,7 instead of LOOKUP 18,19 as GLSIZE
!                                                           E.Leung
!    6.0  09/07/03  Porting to NEC   E.Leung
!    6.0  10/02/04  Correct LOOKUP values for LBLREC and LBNREC in
!                   WGDOS packed fields when field is not being
!                   expanded.        P.Dando
!    6.1  18/08/04  Correct problems converting EXTRA data on
!                   little-endian platforms and remove resultant
!                   redundant array.                     P.Dando
!    6.2  21/11/05  Add level and forecast time to print out in
!                   ATMOS_CONVIEEE. D.Robinson
!    6.2  22/08/05  Improve for free format conversion. P.Selwood
!    6.2  10/4/2005 Removed calls to ABORT.  J. Gill
!    6.2  28/10/05  Correct 64e behaviour and removed redundant
!                   code. P.Selwood
!
!    Purpose: Converts a dump, ancillary or fieldsfile
!             (atmosphere or ocean) between different
!             into 32-bit or 64-bit IEEE format or vice-versa.
!             The following conversions are supported:-
!               On a IEEE machine
!                 64-bit IEEE to 32-bit IEEE
!                 64-bit IEEE to 64-bit IEEE
!             In either case, WGDOS data will be unpacked if
!             requested.
!
!             MAIN_CONVIEEE reads in fixed length and integer
!             headers of UM file to be converted, extracts dimensions
!             of file and then passes these values to
!             subroutine CONVIEEE.
!
!            CONVIEEE reads in headers and data fields from unit NFTIN
!            converts them to IEEE format and writes them to NFTOUT.
!
!    Documentation: UM Doc Paper F5
!
!LLEND----------------------------------------------------------------
!*L  Arguments:-------------------------------------------------------
      SUBROUTINE CONVIEEE                                               &
     &  (LEN_FIXHD,LEN_INTHD,LEN_REALHD,                                &
     &  LEN1_LEVDEPC,LEN2_LEVDEPC,LEN1_ROWDEPC,                         &
     &  LEN2_ROWDEPC,LEN1_COLDEPC,LEN2_COLDEPC,                         &
     &  LEN1_FLDDEPC,LEN2_FLDDEPC,LEN_EXTCNST,                          &
     &  LEN_DUMPHIST,LEN_CFI1,LEN_CFI2,LEN_CFI3,                        &
     &  LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA,P_FIELD,                       &
     &  P_ROWS,ROW_LENGTH,                                              &
     &  NFTIN,NFTOUT,IEEE_TYPE,                                         &
     &  MAX_FIELD_SIZE, WGDOS_EXPAND)

      IMPLICIT NONE

      INTEGER                                                           &

     & LEN_FIXHD                                                        &
                    !IN Length of fixed length header on input file
     &,LEN_INTHD                                                        &
                    !IN Length of integer header on input file
     &,LEN_REALHD                                                       &
                    !IN Length of real header on input file
     &,LEN1_LEVDEPC                                                     &
                    !IN 1st dim of lev dependent consts on input file
     &,LEN2_LEVDEPC                                                     &
                    !IN 2nd dim of lev dependent consts on input file
     &,LEN1_ROWDEPC                                                     &
                    !IN 1st dim of row dependent consts on input file
     &,LEN2_ROWDEPC                                                     &
                    !IN 2nd dim of row dependent consts on input file
     &,LEN1_COLDEPC                                                     &
                    !IN 1st dim of col dependent consts on input file
     &,LEN2_COLDEPC                                                     &
                    !IN 2nd dim of col dependent consts on input file
     &,LEN1_FLDDEPC                                                     &
                    !IN 1st dim of field dependent consts on input fi
     &,LEN2_FLDDEPC                                                     &
                    !IN 2nd dim of field dependent consts on input fi
     &,LEN_EXTCNST                                                      &
                    !IN Length of extra consts on input file
     &,LEN_DUMPHIST                                                     &
                    !IN Length of history header on input file
     &,LEN_CFI1                                                         &
                    !IN Length of index1 on input file
     &,LEN_CFI2                                                         &
                    !IN Length of index2 on input file
     &,LEN_CFI3                                                         &
                    !IN Length of index3 on input file
     &,LEN1_LOOKUP                                                      &
                    !IN 1st dim of LOOKUP on input file
     &,LEN2_LOOKUP                                                      &
                    !IN 2nd dim of LOOKUP on input file
     &,LEN_DATA                                                         &
                    !IN Length of data on input file
     &,P_FIELD                                                          &
                    !IN No of p-points per level on input file
     &,P_ROWS                                                           &
     &,ROW_LENGTH                                                       &
     &,MAX_FIELD_SIZE !Maximum field size on file
      integer wgdos_expand  ! set to 1 for expansion of WGDOS Fields

      INTEGER                                                           &
     & NFTIN                                                            &
                !IN Unit number of input UM dump
     &,NFTOUT                                                           &
                !IN Unit number of output IEEE dump
     &,IEEE_TYPE ! Output file precision

! Local arrays:---------------------------------------------------------
      INTEGER                                                           &
     & FIXHD(LEN_FIXHD),                                                &
                                                 !
     & INTHD(LEN_INTHD),                                                &
                                                 !\  integer
     & CFI1(LEN_CFI1+1),CFI2(LEN_CFI2+1),                               &
                                                 ! > file headers
     & CFI3(LEN_CFI3+1),                                                &
                                                 !/
     & LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP),                                 &
                                                 !
     & LOOKUP_21(LEN2_LOOKUP)                                           &
                                ! Holds values of input LOOKUP(21,K)
     &,LOOKUP_LBNREC(LEN2_LOOKUP)                                       &
     &,lookup_lblrec(len2_lookup)                                       &
     &,lookup_lbegin(len2_lookup)                                       &
                                       ! Old value of lbegin in lookup
     &,lookup_lblrec_new(len2_lookup)                                   &
                                       ! New value of lblrec in lookup
     &,lookup_lbnrec_new(len2_lookup)                                   &
                                       ! New value of lbnrec in lookup
     &,lookup_lbegin_new(len2_lookup)                                   &
                                       ! New value of lbegin in lookup
     &,disk_address                                                     &
                                       ! Current rounded disk address
     &,number_of_data_words_on_disk                                     &
                                       ! Number of data words on disk
     &,number_of_data_words_in_memory                                   &
                                       ! Number of Data Words in memory
     &,old_fixhd_160                                                    &
                                       ! Input value of FIXHD(160)
     &,new_fixhd_160                   ! Output value of FIXHD(160)

      REAL                                                              &
     & REALHD(LEN_REALHD),                                              &
     & LEVDEPC(1+LEN1_LEVDEPC*LEN2_LEVDEPC),                            &
                                                 !
     & ROWDEPC(1+LEN1_ROWDEPC*LEN2_ROWDEPC),                            &
                                                 !
     & COLDEPC(1+LEN1_COLDEPC*LEN2_COLDEPC),                            &
                                                 !\  real
     & FLDDEPC(1+LEN1_FLDDEPC*LEN2_FLDDEPC),                            &
                                                 ! > file headers
     & EXTCNST(LEN_EXTCNST+1),                                          &
                                                 !/
     & DUMPHIST(LEN_DUMPHIST+1)
      INTEGER                                                           &
     & D1(MAX_FIELD_SIZE)  ! Data array used to read in each field

! External subroutines called:------------------------------------------
      EXTERNAL READHEAD,WRITHEAD,ABORT,READFLDS,ABORT_IO,HDPPXRF,GETPPX
!*----------------------------------------------------------------------
!*L  Local variables:---------------------------------------------------

      INTEGER                                                           &
     & ICODE                                                            &
                    ! Error return code from subroutines
     &,START_BLOCK                                                      &
                    ! READHEAD argument (not used)
     &,I,J,K,L                                                          &
                    ! Loop indices
     &,LEN_IO                                                           &
                    ! I/O length
     &,ITYPE                                                            &
                    ! Conversion type
     &,MODEL                                                            &
                        ! Internal model number
     &,SECTION                                                          &
                        ! Section number
     &,ITEM                                                             &
                        ! Item code
     &,JOC_NO_SEAPTS                                                    &
                        ! Number of points in compressed ocean field
     &,LEN_OCFLD                                                        &
                        ! Length of uncompressed ocean field
     &,INIT_FIXHD_161                                                   &
                        ! Initialised value of FIXHD(161)
     &,PPXREF_GRID_TYPE                                                 &
                        ! Grid type form ppxref
     &,LEN_BUF                                                          &
                        ! Record length of boundary dataset
     &,MAX_LEN_BUF                                                      &
                        ! Maximum record length of boundary dataset
     &,POS              ! Position of field in file

! Input arguments for decompose_smexe
      INTEGER                                                           &

     &  global_row_len,                                                 &
                          ! IN  :number of E-W points of entire model
     &  global_n_rows     ! IN  :number of P rows of entire mode

      INTEGER EXPPXI
      EXTERNAL EXPPXI
      INTEGER RowNumber


      REAL A        !Return code from BUFFIN; -1.0 = O.K.

      CHARACTER                                                         &
     & CMESSAGE*100 ! Character string returned if ICODE  /=  0
      INTEGER NFT1,NFT2
      PARAMETER (NFT1=22, NFT2=2)
!*----------------------------------------------------------------------
#include "clookadd.h"
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "cstash.h"
#include "c_mdi.h"
#include "cntl_io.h"
#include "decomptp.h"
!L 0. Read in PPXREF
      cmessage = ' '
      ppxRecs=1
      RowNumber=0
      ICODE=0
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_A',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_A'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('CONVIEEE', ICODE,                                 &
     &   'Error reading STASHmaster_A')

      END IF
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_O',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_O'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('CONVIEEE', ICODE,                                 &
     &   'Error reading STASHmaster_O')

      END IF
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_S',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_S'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('CONVIEEE', ICODE,                                 &
     &   'Error reading STASHmaster_S')

      END IF
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_W',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_W'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('CONVIEEE', ICODE,                                 &
     &   'Error reading STASHmaster_W')

      END IF

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
        CALL EREPORT('CONVIEEE', ICODE,                                 &
     &   CMESSAGE)

      END IF

!User STASHmaster
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(0,' ',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('CONVIEEE', ICODE,                                 &
     &   CMESSAGE)

      END IF

! DEPENDS ON: getppx
      CALL GETPPX(0,NFT2,' ',RowNumber,                                 &
#include "argppx.h"
     &            ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('CONVIEEE', ICODE,                                 &
     &   CMESSAGE)

      END IF

!L 1. Read in file header

! DEPENDS ON: readhead
      CALL READHEAD(NFTIN,FIXHD,LEN_FIXHD,                              &
     &              INTHD,LEN_INTHD,REALHD,LEN_REALHD,                  &
     &              LEVDEPC,LEN1_LEVDEPC,LEN2_LEVDEPC,                  &
     &              ROWDEPC,LEN1_ROWDEPC,LEN2_ROWDEPC,                  &
     &              COLDEPC,LEN1_COLDEPC,LEN2_COLDEPC,                  &
     &              FLDDEPC,LEN1_FLDDEPC,LEN2_FLDDEPC,                  &
     &              EXTCNST,LEN_EXTCNST,DUMPHIST,LEN_DUMPHIST,          &
     &              CFI1,LEN_CFI1,CFI2,LEN_CFI2,CFI3,LEN_CFI3,          &
     &              LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA,            &
#include "argppx.h"
     &              START_BLOCK,ICODE,CMESSAGE)

      IF(ICODE /= 0)THEN
        WRITE(6,*)CMESSAGE,ICODE
! DEPENDS ON: ereport
        CALL EREPORT('CONVIEEE', ICODE,                                 &
     &   CMESSAGE)

      ENDIF

! 2: Check for PP format dataset if field to be expanded

! Preserve the original length values for re-use
      DO K=1,LEN2_LOOKUP
        LOOKUP_LBLREC(K)=LOOKUP(LBLREC,K)
        LOOKUP_LBEGIN(K)=LOOKUP(LBEGIN,K)
        LOOKUP_LBNREC(K)=LOOKUP(LBNREC,K)
      END DO

! Get decompostion information
! ----------------------------

! DEPENDS ON: decompose_smexe
        CALL DECOMPOSE_SMEXE(ROW_LENGTH, P_ROWS,0,0,-99)
! DEPENDS ON: change_decomposition
        CALL CHANGE_DECOMPOSITION(decomp_smexe,ICODE)
      IF(LOOKUP(LBNREC,1) >  0.AND.FIXHD(12) >  0)THEN
! Check for WGDOS expansion
        IF (WGDOS_EXPAND == 1)THEN
! Issue a message on why we are doing this
          write(6,'(//''***** Initial Scan for PP Format Dataset'',     &
     &     '' *****''/)')
          DO I=1,LEN2_LOOKUP
            IF(LOOKUP(1,I) == -99) GOTO 195
! DEPENDS ON: readflds
            CALL READFLDS(NFTIN,1,I,LOOKUP,LEN1_LOOKUP,                 &
     &                    D1,MAX_FIELD_SIZE,FIXHD,                      &
#include "argppx.h"
     &                    WGDOS_EXPAND,ICODE,CMESSAGE)
! DEPENDS ON: abort_io
            IF(ICODE /= 0)CALL ABORT_IO('CONVIEEE',CMESSAGE,ICODE,NFTIN)
          END DO
195       CONTINUE
        END IF
      END IF

! 3: Reset LOOKUP and FIXHD

      INIT_FIXHD_161=0
      OLD_FIXHD_160=FIXHD(160)
      DO K=1,LEN2_LOOKUP

        IF(LOOKUP(1,K) == -99)GOTO 200
! Set LOOKUP(LBNREC) = 0  in old dumps where UM version number
!   not in fixed length header
        IF(FIXHD(12) <  0.AND.FIXHD(5) /= 3)THEN
          LOOKUP(LBNREC,K)=0
        END IF

! Packing code = -2 now obselete, reset packing code to 2
        IF(LOOKUP(LBPACK,K) == -2)LOOKUP(LBPACK,K)=2

! Store values of packing indicator and set least significant
!   number in LOOKUP(LBPACK,K) to 0 to indicate no packing
        LOOKUP_21(K)=LOOKUP(LBPACK,K)
        LOOKUP(LBPACK,K)=MOD(LOOKUP(LBPACK,K),1000)
          IF((MOD(LOOKUP(LBPACK,K),10) /= 1) .AND.                      &
     &                     (MOD(LOOKUP(LBPACK,K),10) /= 4))THEN
            LOOKUP(LBPACK,K)=(LOOKUP(LBPACK,K)/10)*10
          ELSE IF(WGDOS_EXPAND == 1) THEN
            LOOKUP(LBPACK,K)=(LOOKUP(LBPACK,K)/10)*10
          END IF
! For option IEEE_TYPE=32 but the field is WGDOS packed but NOT
! being expanded then need to double the size of LOOKUP(LBLREC,K)
! before the call to set_dumpfile_address so that the output field
! size is calculated correctly
          IF (IEEE_TYPE  ==  32 .AND. WGDOS_EXPAND  ==  0               &
     &                  .AND. (MOD(LOOKUP(LBPACK,K),10)  ==  1)) THEN
            LOOKUP(LBLREC,K) = 2*LOOKUP(LBLREC,K)
          ENDIF
! Process compressed fields
          IF(MOD(LOOKUP(LBPACK, K),1000) == 110)THEN
            IF(K <= (INTHD(14)+2)*INTHD(8))THEN
! Calculate expanded field lengths for ocean compressed fields
              MODEL=LOOKUP(MODEL_CODE, K)
              ITEM=MOD(LOOKUP(ITEM_CODE, K),1000)
              SECTION=(LOOKUP(ITEM_CODE, K)-ITEM)/1000
! DEPENDS ON: exppxi
              PPXREF_GRID_TYPE=EXPPXI(MODEL,SECTION,ITEM,PPX_GRID_TYPE, &
#include "argppx.h"
     &                                ICODE,CMESSAGE)
              IF(PPXREF_GRID_TYPE == 36)THEN
! Ocean mass points.
                LOOKUP(LBNPT,K)   = INTHD(6)
                LOOKUP(LBROW,K)   = INTHD(7)
                LOOKUP(LBLREC, K) = INTHD(6)*INTHD(7)
              ELSEIF(PPXREF_GRID_TYPE == 37)THEN
! Ocean velocity points. One less row.
                LOOKUP(LBNPT,K)   = INTHD(6)
                LOOKUP(LBROW,K)   = INTHD(7)-1
                LOOKUP(LBLREC, K) = INTHD(6)*(INTHD(7)-1)
              END IF
              LOOKUP(LBPACK, K) = 0

            ELSE
! Field not compressed onto sea points. Correct packing code
              LOOKUP(LBPACK, K)=MOD(LOOKUP(LBPACK, K),10)
            END IF

          END IF
! Add to length of data
          INIT_FIXHD_161=INIT_FIXHD_161+LOOKUP(LBLREC,K)

        END DO

200     CONTINUE
        FIXHD(160)=FIXHD(150)+FIXHD(151)*FIXHD(152)
        FIXHD(161)=INIT_FIXHD_161
        LEN_DATA=INIT_FIXHD_161

        DO K=1,LEN2_LOOKUP
!  indicate output format.
          LOOKUP(LBPACK,K)=LOOKUP(LBPACK,K)+3000
        END DO

      IF(FIXHD(12) <  208)FIXHD(12)=208

! Boundary datasets are structured differently.
! Skip call to set_dumpfile_address for boundary datasets and
! Calculate addressing for well formed  boundary datasets explicitly.
      IF (FIXHD(5) /= 5)THEN

! Not a boundary dataset. Call set_dumpfile_address
!
!--reset the 32/64 bit lookup headers after packing, etc
!  has been removed
! DEPENDS ON: set_dumpfile_address
      call set_dumpfile_address(fixhd, len_fixhd,                       &
     &                          lookup, len1_lookup,                    &
     &                          len2_lookup,                            &
     &                          number_of_data_words_in_memory,         &
     &                          number_of_data_words_on_disk,           &
     &                          disk_address)
      ELSE

! Boundary  dataset. Calcuate start address from header and round it up
! to ensure we start on a sector boundary
        DISK_ADDRESS=FIXHD(160)-1
        DISK_ADDRESS=((DISK_ADDRESS+UM_SECTOR_SIZE-1)/                  &
     &                UM_SECTOR_SIZE)*UM_SECTOR_SIZE

! Loop over number of times for which data is present in dataset
        DO K=1,INTHD(3)

! Loop over number of different field types present
          LEN_BUF=0
          MAX_LEN_BUF=0
          DO I=1,INTHD(15)
            POS=(K-1)*INTHD(15)+I
            LOOKUP(LBEGIN,POS)=DISK_ADDRESS+LEN_BUF
            LOOKUP(LBNREC,POS)=LOOKUP(LBLREC,POS)
            LEN_BUF=LEN_BUF+LOOKUP(LBLREC,POS)
          END DO
          MAX_LEN_BUF=MAX0(LEN_BUF,MAX_LEN_BUF)

! Update disk address and ensure that next time starts
! on a sector boundary
          DISK_ADDRESS=DISK_ADDRESS+LEN_BUF
          DISK_ADDRESS=((DISK_ADDRESS+UM_SECTOR_SIZE-1)/                &
     &                  UM_SECTOR_SIZE)*UM_SECTOR_SIZE

        END DO

      END IF
!--preserve the new length values for re-use
      new_fixhd_160=fixhd(160)
      do k=1,len2_lookup
        lookup_lblrec_new(k)=lookup(lblrec, k)
        lookup_lbnrec_new(k)=lookup(lbnrec, k)
        lookup_lbegin_new(k)=lookup(lbegin, k)
      end do
!L 1. Write out file header

! DEPENDS ON: writhead
      CALL WRITHEAD(NFTOUT,FIXHD,LEN_FIXHD,                             &
     &              INTHD,LEN_INTHD,REALHD,LEN_REALHD,                  &
     &              LEVDEPC,LEN1_LEVDEPC,LEN2_LEVDEPC,                  &
     &              ROWDEPC,LEN1_ROWDEPC,LEN2_ROWDEPC,                  &
     &              COLDEPC,LEN1_COLDEPC,LEN2_COLDEPC,                  &
     &              FLDDEPC,LEN1_FLDDEPC,LEN2_FLDDEPC,                  &
     &              EXTCNST,LEN_EXTCNST,DUMPHIST,LEN_DUMPHIST,          &
     &              CFI1,LEN_CFI1,CFI2,LEN_CFI2,CFI3,LEN_CFI3,          &
     &              LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA,            &
#if defined(IEEE)
     &              IEEE_TYPE,                                          &
#endif
#include "argppx.h"
     &              START_BLOCK,ICODE,CMESSAGE)

      IF(ICODE /= 0)THEN
        WRITE(6,*)CMESSAGE,ICODE
! DEPENDS ON: ereport
        CALL EREPORT('CONVIEEE', ICODE,CMESSAGE)

      ENDIF

! Reset PP file indicator
      IF(LOOKUP_LBNREC(1) >  0)THEN
        DO K=1,LEN2_LOOKUP
          LOOKUP(LBNREC,K)=LOOKUP_LBNREC(K)
          lookup(lblrec,k)=lookup_lblrec(k)
          lookup(lbegin,k)=lookup_lbegin(k)
        ENDDO
      ENDIF
! Restore value of packing indicator

      DO K=1,LEN2_LOOKUP
        LOOKUP(21,K)=LOOKUP_21(K)
      ENDDO

!L 3. Read in each field, convert to IEEE format and write out
!L    results to new file

      IF (FIXHD(2) == 1)THEN

! Atmosphere file
! DEPENDS ON: atmos_convieee
        CALL ATMOS_CONVIEEE(NFTIN,NFTOUT,IEEE_TYPE,MAX_FIELD_SIZE       &
     &,                         LEN_FIXHD,LEN1_LOOKUP,LEN2_LOOKUP       &
     &,                         OLD_FIXHD_160,NEW_FIXHD_160             &
     &,                         LOOKUP_LBLREC                           &
     &,                         LOOKUP_LBEGIN,LOOKUP_LBNREC             &
     &,                         LOOKUP_LBEGIN_NEW,LOOKUP_LBNREC_NEW,    &
#include "argppx.h"
     &                          FIXHD,LOOKUP,WGDOS_EXPAND)

      ELSEIF (FIXHD(2) == 2)THEN

! Ocean file

! Calculate sizes of compressed and uncompressed ocean fields:
! First decide whether there are any compressed fields (use LBPACK
! to determine whether the first field contains sea points only)
        IF( MOD(LOOKUP(LBPACK, 1)/10,10)  ==  0) THEN
          JOC_NO_SEAPTS=1
          LEN_OCFLD    =1
        ELSE
        JOC_NO_SEAPTS=INTHD(11)
        LEN_OCFLD    =INTHD(6)*INTHD(7)*INTHD(8)
        END IF

! DEPENDS ON: ocean_convieee
        CALL OCEAN_CONVIEEE(NFTIN,NFTOUT,IEEE_TYPE,MAX_FIELD_SIZE       &
     &,                         LEN_FIXHD,LEN_INTHD                     &
     &,                         LEN_CFI1,LEN_CFI2,LEN_CFI3              &
     &,                         LEN1_LOOKUP,LEN2_LOOKUP                 &
     &,                         JOC_NO_SEAPTS,LEN_OCFLD                 &
     &,                         OLD_FIXHD_160,NEW_FIXHD_160             &
     &,                         LOOKUP_LBLREC,LOOKUP_LBLREC_NEW         &
     &,                         LOOKUP_LBEGIN,LOOKUP_LBEGIN_NEW         &
     &,                         LOOKUP_LBNREC,LOOKUP_LBNREC_NEW         &
     &,                         FIXHD,INTHD,LOOKUP,CFI1,CFI2,CFI3,      &
#include "argppx.h"
     &                          WGDOS_EXPAND)

      END IF

      WRITE(6,'(I4,'' fields have been converted'')') LEN2_LOOKUP

      RETURN
      END SUBROUTINE CONVIEEE


!----------------------------------------------------------------------

#endif
