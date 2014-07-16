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


!----------------------------------------------------------------------

      SUBROUTINE ATMOS_CONVIEEE(NFTIN,NFTOUT,IEEE_TYPE,MAX_FIELD_SIZE   &
     &,                         LEN_FIXHD,LEN1_LOOKUP,LEN2_LOOKUP       &
     &,                         OLD_FIXHD_160,NEW_FIXHD_160             &
     &,                         LOOKUP_LBLREC                           &
     &,                         LOOKUP_LBEGIN,LOOKUP_LBNREC             &
     &,                         LOOKUP_LBEGIN_NEW,LOOKUP_LBNREC_NEW,    &
#include "argppx.h"
     &                          FIXHD,LOOKUP,WGDOS_EXPAND)
      IMPLICIT NONE
      INTEGER IEEE_TYPE
      INTEGER LEN_FIXHD
      INTEGER LEN1_LOOKUP
      INTEGER LEN2_LOOKUP
      INTEGER MAX_FIELD_SIZE
      INTEGER NEW_FIXHD_160
      INTEGER NFTIN
      INTEGER NFTOUT
      INTEGER OLD_FIXHD_160
      INTEGER WGDOS_EXPAND

      INTEGER FIXHD(LEN_FIXHD)
      INTEGER LOOKUP_LBLREC(LEN2_LOOKUP)
      INTEGER LOOKUP_LBEGIN(LEN2_LOOKUP)
      INTEGER LOOKUP_LBNREC(LEN2_LOOKUP)
      INTEGER LOOKUP_LBEGIN_NEW(LEN2_LOOKUP)
      INTEGER LOOKUP_LBNREC_NEW(LEN2_LOOKUP)
      INTEGER LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)
! Local arrays:--------------------------------------------------------
      INTEGER D1(MAX_FIELD_SIZE) ! Data array used to read in each field
      REAL IEEE_64(MAX_FIELD_SIZE) !Array containing 64 bit IEEE data
! External subroutines called:-----------------------------------------
      EXTERNAL ABORT,READFLDS,ABORT_IO,BUFFO32,PR_LOOK
#if defined(CRAY)
      EXTERNAL CRI2IEG,CRAY2CRI,CRI2CRAY
#endif
      EXTERNAL INT_FROM_REAL
!----------------------------------------------------------------------
! Local variables:-----------------------------------------------------
      INTEGER      I,J      ! Loop variables
      INTEGER      K        ! Return code from CRAY intrinsic functions
      INTEGER      ICODE    ! Error return code from READFLDS
      INTEGER      ITYPE    ! Conversion type
      INTEGER      LEN_IO   ! I/O length
      INTEGER      TOT_VALUES   ! Normal data + extra data (if applied)
      INTEGER      IEXTRAW      ! Extra data length
      INTEGER      NET_DATA_LEN ! Normal data length
      INTEGER      ADDR         ! Start address for extra data
      INTEGER      IEEE_ADDR    ! Start address for extra data in IEEE
      INTEGER      BIT_OFF      ! Bit offset in IEEE array
      INTEGER      CODE         ! Encoded info for extra data set
      INTEGER      DATA_VALUES  ! Decoded real extra data length
      INTEGER      PACK_CODE    ! Packing code
      INTEGER      IEEE2IEG
      INTEGER      INT_FROM_REAL ! Function to convert real to int
      REAL         A        ! Return code from BUFFIN; -1.0 = O.K.
      CHARACTER*80 CMESSAGE ! Character string returned if ICODE  /=  0
!----------------------------------------------------------------------
! Input arguments for decompose_smexe
      INTEGER                                                           &

     &  global_row_len,                                                 &
                          ! IN  :number of E-W points of entire model
     &  global_n_rows    ! IN  :number of P rows of entire mode


#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "cstash.h"
#include "clookadd.h"
#include "c_mdi.h"
#include "cprintst.h"

! Loop over all fields
      DO I=1,LEN2_LOOKUP

! Check for the end of a PP format lookup table
        IF(LOOKUP(1,I) == -99) GOTO 2000

! Reset the headers in case WDGOS packing has altered them
        IF(WGDOS_EXPAND == 1) THEN
        LOOKUP(LBLREC,I)=LOOKUP_LBLREC(I)
        LOOKUP(LBEGIN,I)=LOOKUP_LBEGIN(I)
        LOOKUP(LBNREC,I)=LOOKUP_LBNREC(I)
        FIXHD(160)=OLD_FIXHD_160
        ENDIF

! Check if this field has already been converted - WGDOS only
        IF(MOD(LOOKUP(21,I),10) == 1 .AND. WGDOS_EXPAND /= 1) THEN
! Read in field
          IF(IEEE_TYPE == 32)THEN
! DEPENDS ON: readflds
            CALL READFLDS(NFTIN,1,I,LOOKUP,LEN1_LOOKUP,                 &
     &                    IEEE_64,MAX_FIELD_SIZE,FIXHD,                 &
#include "argppx.h"
     &                    WGDOS_EXPAND,ICODE,CMESSAGE)
          ELSEIF(IEEE_TYPE == 64)THEN
! Read in field
! DEPENDS ON: readflds
            CALL READFLDS(NFTIN,1,I,LOOKUP,LEN1_LOOKUP,                 &
     &                    IEEE_64,MAX_FIELD_SIZE,FIXHD,                 &
#include "argppx.h"
     &                    WGDOS_EXPAND,ICODE,CMESSAGE)
          END IF
! DEPENDS ON: abort_io
          IF(ICODE /= 0)CALL ABORT_IO('CONVIEEE',CMESSAGE,ICODE,NFTIN)

        ELSE IF(MOD(LOOKUP(21,I),10) == 1 .AND. WGDOS_EXPAND == 1 .AND. &
     &          IEEE_TYPE == 64 ) THEN
! Read in field
! DEPENDS ON: readflds
          CALL READFLDS(NFTIN,1,I,LOOKUP,LEN1_LOOKUP,                   &
     &                  IEEE_64,MAX_FIELD_SIZE,FIXHD,                   &
#include "argppx.h"
     &                  WGDOS_EXPAND,ICODE,CMESSAGE)
! DEPENDS ON: abort_io
          IF(ICODE /= 0)CALL ABORT_IO('CONVIEEE',CMESSAGE,ICODE,NFTIN)

        ELSE
! Read in field
! DEPENDS ON: readflds
          CALL READFLDS(NFTIN,1,I,LOOKUP,LEN1_LOOKUP,                   &
     &                  D1,MAX_FIELD_SIZE,FIXHD,                        &
#include "argppx.h"
     &                  WGDOS_EXPAND,ICODE,CMESSAGE)
! DEPENDS ON: abort_io
          IF(ICODE /= 0)CALL ABORT_IO('CONVIEEE',CMESSAGE,ICODE,NFTIN)

! Set data type
          IF(ABS(LOOKUP(DATA_TYPE,I)) == 1) THEN
! Type real
            IF(IEEE_TYPE == 32)THEN
              ITYPE=3
            ELSEIF(IEEE_TYPE == 64)THEN
              ITYPE=2
            ENDIF
          ELSE IF(ABS(LOOKUP(DATA_TYPE,I)) == 2) THEN
! Type integer
            IF(IEEE_TYPE == 32)THEN
              ITYPE=2
            ELSEIF(IEEE_TYPE == 64)THEN
              ITYPE=1
            ENDIF
          ELSE IF(ABS(LOOKUP(DATA_TYPE,I)) == 3) THEN
! Type logical
            ITYPE=5
          ELSE
! DEPENDS ON: pr_look
            CALL PR_LOOK(                                               &
#include "argppx.h"
     &                   LOOKUP,LOOKUP,LEN1_LOOKUP,I)
            ICODE=3
            CMESSAGE='CONVIEEE: Invalid code in LOOKUP(39,K)'
            RETURN
          ENDIF

! Convert to IEEE format and write to disk
          IF(ITYPE >= 0)THEN
            IF(IEEE_TYPE == 32)THEN

              TOT_VALUES=LOOKUP(LBLREC,I)
              IEXTRAW=0
              IF(LOOKUP(LBEXT,I) >  0) THEN ! got some extra data
                IEXTRAW=LOOKUP(LBEXT,I)

                ! Check that no pack_code 1 file left packed if we have
                ! extra data. (files with pack_code 2 & 3 are forced
                ! unpacked).  This version supports pack_code 4 data

                PACK_CODE=MOD((LOOKUP(LBPACK,I)),10)
                IF((LOOKUP(LBROW,I)*LOOKUP(LBNPT,I)+LOOKUP(LBEXT,I)  /= &
     &              LOOKUP(LBLREC,I)).AND.                              &
     &              (PACK_CODE /= 4)) THEN
                  CMESSAGE=                                             &
     &            'CONVIEE1 : Packing of extra data not supported'
                  WRITE(6,*)'Please use expand option'
                  ICODE=1
! DEPENDS ON: abort_io
                  CALL ABORT_IO('CONVIEEE',CMESSAGE,ICODE,NFTIN)
                ENDIF
              ENDIF

              ! Converting data into IEEE format consists of 3 stages
              ! a) convert normal grid data into IEEE format
              ! b) convert integer header of extra data vector into IEEE
              ! c) convert rest of extra data vector into IEEE format

              ! Process normal data

              NET_DATA_LEN=LOOKUP(LBLREC,I)-IEXTRAW
              IF (NET_DATA_LEN  ==  0) THEN
! In some cases (e.g., CovStats files), fields contain no data.  Instead
! of giving an error from the data conversion  routines, output a
! message and skip conversion.
                WRITE(6,*) 'Data length = 0 for field ',i,              &
     &                     ' - nothing to convert'
              ELSE
              K=IEEE2IEG(ITYPE,NET_DATA_LEN,IEEE_64,0,                  &
     &                 D1,1,64,IEEE_TYPE)
              IF(K /= 0)THEN
                WRITE(6,*)'Conversion Error - Return Code is ',K
                DO J=1,NET_DATA_LEN
                    WRITE(6,'(''Error converting field '',i5,           &
     &                        '' : Stash Code '',i5,                    &
     &                        '' : Point No. '',i5)')                   &
     &                        I, LOOKUP(ITEM_CODE,I),J
                    WRITE(6,*) 'Number unconvertable reset to RMDI'
                    IEEE_64(J)=RMDI
                END DO
              END IF
              END IF
              ! Process extra data

              ! About BIT OFFSET
              !
              ! 1         2         3         4         5 (addr)
              ! |---------|---------|---------|---------|  FIELD
              ! .        .         .
              ! .      .       .
              ! .    .    .
              ! |----|----|----|----|----|----|----|----|  IEEE_FIELD
              ! 1         2         3         4         5 (ieee_addr)
              !                     |    |
              ! <--------->         |    |
              !  a "word"           | bit_off=32
              !                 bit_off=0
              ! Example:
              !  if ADDR=2, IEEE_ADDR=3/2=1
              !  IEEE_ADDR*2 eq 1;  so BIT_OFF=32
              !


              IF (IEXTRAW >  0) THEN ! process extra data as got some
                ! init values for while loop
                ! start address in field for extra data
                ADDR=TOT_VALUES-IEXTRAW+1
                IEEE_ADDR=(ADDR+1)/2
                IF (IEEE_ADDR*2 == ADDR) THEN
                  BIT_OFF=32
                ELSE
                  BIT_OFF=0
                ENDIF

                DO WHILE (ADDR <  TOT_VALUES)
                ! CODE is integer header for extra data vector
                ! which contains encoded info for extra data -
                ! vector length & data type
                ! Decode CODE: data_values will be vector length
                ! Details about extra data, see Paper F3
                ! NB. integer header for extra data vector is converted
                !     to its real EQUIVALENCE during model run.  Hence,
                !     INT_FROM_REAL serves to convert it back to INTEGER

! DEPENDS ON: int_from_real
                  CODE=INT_FROM_REAL(D1(ADDR))
! DEPENDS ON: check_extra
                  CALL CHECK_EXTRA(CODE,DATA_VALUES,ICODE,CMESSAGE)
                  IF (ICODE /= 0) THEN
                    write(6,*)'Fail in CHECK_EXTRA'
                    RETURN
                  ENDIF

                  ! Convert integer extra_data head into IEEE
                  K=IEEE2IEG(2,1,IEEE_64(IEEE_ADDR),BIT_OFF,            &
     &                      D1(ADDR),1,64,IEEE_TYPE)

                  IF (K /= 0) THEN
                    ICODE=1
                    CMESSAGE=                                           &
     &              'CONVIEE1: failed in integer conv of extra data'
                    RETURN
                  ENDIF

                  ! update bit_off, addr and ibm_addr
                  IF (BIT_OFF == 0) THEN
                    BIT_OFF=32
                  ELSE
                    BIT_OFF=0
                    IEEE_ADDR=IEEE_ADDR+1 ! GONE ON ANOTHER WORD..
                  ENDIF
                  ADDR=ADDR+1           ! INCREMENT ADDRESS

                  ! Convert REAL vector to IEEE format.
                  K=IEEE2IEG(3,DATA_VALUES,IEEE_64(IEEE_ADDR),          &
     &                      BIT_OFF,D1(ADDR),1,64,IEEE_TYPE)

                  IF (K /= 0) THEN
                    ICODE=1
                    CMESSAGE=                                           &
     &              'CONVIEE1: FAILED IN REAL CONV OF EXTRA DATA'
                    RETURN
                  ENDIF

                  ! Update loop variables.
                  ADDR=ADDR+DATA_VALUES
                  IEEE_ADDR=IEEE_ADDR+DATA_VALUES/2
                  ! Odd no. of values.
                  IF ((DATA_VALUES/2)*2 /= DATA_VALUES) THEN
                    IF (BIT_OFF == 0) THEN
                      BIT_OFF=32
                    ELSE
                      BIT_OFF=0
                      IEEE_ADDR=IEEE_ADDR+1 ! GONE ON ANOTHER WORD..
                    ENDIF
                  ENDIF
                ENDDO                 ! continue until run out of data

                ! Verify addr and ieee_addr have correct values at end
                ! of whileloop. First check that addr is ok
                IF (ADDR /= TOT_VALUES+1) THEN
                  WRITE(CMESSAGE,109)ADDR,TOT_VALUES+1
 109              FORMAT('CONVIEE1: addr',i5,1x,'<> tot_values+1',i5)
                  ICODE=1
                  RETURN
                ENDIF
                ! and so is ieee_addr
                IF (BIT_OFF == 0) IEEE_ADDR=IEEE_ADDR-1
                IF (IEEE_ADDR /= (TOT_VALUES+1)/2) THEN
                  WRITE(CMESSAGE,110)IEEE_ADDR,(TOT_VALUES+1)/2
 110              FORMAT('CONVIEE1: ieee_addr ',i5,1x,                  &
     &                   ' <> (tot_values+1)/2',i5)
                  ICODE=1
                  RETURN
                ENDIF
              ENDIF ! end processing of extra data
            END IF

! 64 bit case - just copy into ieee_64 with appropriate type change
            IF (IEEE_TYPE == 64) THEN
              DO K = 1, LOOKUP(LBLREC,I)
                IEEE_64(K) = TRANSFER(D1(K),IEEE_64(1))
              END DO
            END IF

          ELSE
            DO K=1,LOOKUP(LBLREC,I)
              IEEE_64(k)=IAND(D1(K),1)
            END DO
          ENDIF
        ENDIF

! Write out field
        FIXHD(160)=NEW_FIXHD_160
        IF(IEEE_TYPE == 32)THEN
          CALL SETPOS32(NFTOUT, LOOKUP_LBEGIN_NEW(I), K)
! DEPENDS ON: buffo32
          CALL BUFFO32(NFTOUT, IEEE_64, LOOKUP_LBNREC_NEW(I), LEN_IO, A)
        ELSEIF(IEEE_TYPE == 64)THEN
! DEPENDS ON: setpos
          CALL SETPOS(NFTOUT, LOOKUP_LBEGIN_NEW(I), K)
! DEPENDS ON: buffout
          CALL BUFFOUT(NFTOUT, IEEE_64, LOOKUP_LBNREC_NEW(I), LEN_IO, A)
        ENDIF

! Check for I/O errors
        if(A /= -1.0.OR.LEN_IO /= LOOKUP_LBNREC_NEW(I)) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of data field',                      &
     &                 A,LEN_IO,LOOKUP(15,I))
! DEPENDS ON: ereport
          CALL EREPORT('ATMOS_CONVIEEE', 1002,                          &
     &     'Buffer out of data field wrong size')

        ENDIF

        If ( PrintStatus >= PrStatus_Diag ) Then
          WRITE(6,'(''Field '',i5,'' : Stash Code '',i5,                &
     &     '' : Level '',i4,                                            &
     &     '' : T + '',i3,                                              &
     &     '' : has been converted'')') I, LOOKUP(42,I)                 &
     &          , LOOKUP(33,I), LOOKUP(14,I)
        End If

! Reset the headers in case WDGOS packing has altered them
        IF(WGDOS_EXPAND == 1) THEN
        LOOKUP(LBLREC,I)=LOOKUP_LBLREC(I)
        LOOKUP(LBEGIN,I)=LOOKUP_LBEGIN(I)
        LOOKUP(LBNREC,I)=LOOKUP_LBNREC(I)
        FIXHD(160)=OLD_FIXHD_160
        ENDIF

      END DO
2000  CONTINUE

      RETURN
      END SUBROUTINE ATMOS_CONVIEEE
#endif
