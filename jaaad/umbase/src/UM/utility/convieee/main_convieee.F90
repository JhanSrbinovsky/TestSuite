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
      PROGRAM MAIN_CONVIEEE

      IMPLICIT NONE

      INTEGER                                                           &
     & FIXHD(256)                                                       &
                      !Space for fixed length header
     &,INTHD(100)     !Space for integer header

      INTEGER                                                           &
     & LEN_FIXHD                                                        &
                      !Length of fixed length header on input file
     &,LEN_INTHD                                                        &
                      !Length of integer header on input file
     &,LEN_REALHD                                                       &
                      !Length of real header on input file
     &,LEN1_LEVDEPC                                                     &
                      !1st dim of lev dependent consts on input file
     &,LEN2_LEVDEPC                                                     &
                      !2nd dim of lev dependent consts on input file
     &,LEN1_ROWDEPC                                                     &
                      !1st dim of row dependent consts on input file
     &,LEN2_ROWDEPC                                                     &
                      !2nd dim of row dependent consts on input file
     &,LEN1_COLDEPC                                                     &
                      !1st dim of col dependent consts on input file
     &,LEN2_COLDEPC                                                     &
                      !2nd dim of col dependent consts on input file
     &,LEN1_FLDDEPC                                                     &
                      !1st dim of field dependent consts on input file
     &,LEN2_FLDDEPC                                                     &
                      !2nd dim of field dependent consts on input file
     &,LEN_EXTCNST                                                      &
                      !Length of extra consts on input file
     &,LEN_DUMPHIST                                                     &
                      !Length of history header on input file
     &,LEN_CFI1                                                         &
                      !Length of index1 on input file
     &,LEN_CFI2                                                         &
                      !Length of index2 on input file
     &,LEN_CFI3                                                         &
                      !Length of index3 on input file
     &,LEN1_LOOKUP                                                      &
                      !1st dim of LOOKUP on input file
     &,LEN2_LOOKUP                                                      &
                      !2nd dim of LOOKUP on input file
     &,LEN_DATA                                                         &
                      !Length of data on input file
     &,ROW_LENGTH                                                       &
                      !No of points E-W on input file
     &,P_ROWS                                                           &
                      !No of p-rows on input file
     &,P_FIELD                                                          &
                      !No of p-points per level on input file
     &,MAX_FIELD_SIZE !Maximum field size on file

      INTEGER                                                           &
     & LEN_IO                                                           &
                !Length of I/O returned by BUFFER IN
     &,I                                                                &
                !Loop index
     &,NFTIN                                                            &
                !Unit number of input UM dump
     &,NFTOUT                                                           &
                !Unit number of output IEEE dump
     &,ERR                                                              &
                !Return code from OPEN
     &,IEEE_TYPE                                                        &
                 ! Output file precision
     &,ICODE                                                            &
                !Return code from setpos
     &,LEN                                                              &
                !Length of string returned by PXFGETARG
     &,IERR     !Return code from PXFGETARG

      CHARACTER*5 INPREC
      CHARACTER *80                                                     &
     & STRING    ! Character string holding command line arg

      REAL A    !Return code from BUFFIN; -1.0 = O.K.

      integer wgdos_expand

      CHARACTER*100 PAREXE_ENV   ! hold name of the // exec script
      INTEGER ME_GC,NPROC_GC
#include "cntl_io.h"

! External subroutines called:------------------------------------------
      EXTERNAL IOERROR,ABORT_IO,BUFFIN,FILE_OPEN,                       &
     &         SETPOS,ABORT,CONVIEEE,InitPrintStatus                    &
     &        ,GC_INIT
#if defined(CRAY)
      EXTERNAL :: PXFGETARG
#endif
!*----------------------------------------------------------------------

      PAREXE_ENV=' '

      CALL GC_INIT(PAREXE_ENV,ME_GC,NPROC_GC)
! Initialise print status for standard output
! DEPENDS ON: initprintstatus
      CALL InitPrintStatus

!--select no WGDOS expansion
      wgdos_expand=0

!L 0. Read in precision of output file
      CALL FORT_GET_ENV("PRECISION",9,INPREC,5,IERR)
      READ(INPREC,'(a5)') string
      LEN = LEN_TRIM(string)
      IF(LEN /= 2.OR.IERR /= 0)THEN
        IEEE_TYPE=32
        if(len == 3 .and. ierr == 0) then
          if(string == '32e' .or. string == '32E') then
            wgdos_expand=1
          else if(string == '64e' .or. string == '64E') then
            ieee_type=64
            wgdos_expand=1
          else
            WRITE(6,*)'Unsupported word length ',STRING
! DEPENDS ON: ereport
            CALL EREPORT('MAIN_CONVIEEE', 1000,                         &
     &       'Unsupported word length')

          endif
        endif
      ELSE
        IF(STRING == '32')THEN
          IEEE_TYPE=32
        ELSEIF(STRING == '64')THEN
          IEEE_TYPE=64
        ELSE
          WRITE(6,*)'Unsupported word length ',STRING
! DEPENDS ON: ereport
          CALL EREPORT('MAIN_CONVIEEE', 1002,                           &
     &     'Unsupported word length')

        ENDIF
      ENDIF
!
      IF(WGDOS_EXPAND == 0) THEN
        WRITE(6,'(/''Conversion to IEEE  '',i2,''-bit Format'',         &
     &    '' with no expansion of WGDOS Fields''/)') ieee_type
      ELSE
        WRITE(6,'(/''Conversion to IEEE  '',i2,''-bit Format'',         &
     &    '' with expansion of WGDOS Fields''/)') ieee_type
      END IF


!L 1. Assign unit numbers

      NFTIN=20
      NFTOUT=21

      WRITE(6,'(20x,''FILE STATUS'')')
      WRITE(6,'(20x,''==========='')')
!     CALL OPEN(1,'PPXREF',6,0,0,ERR)
! DEPENDS ON: file_open
      CALL FILE_OPEN(20,'FILE1',5,0,0,ERR)
! DEPENDS ON: file_open
      CALL FILE_OPEN(21,'FILE2',5,1,0,ERR)


!L 2. Buffer in fixed length header record

! DEPENDS ON: buffin
      CALL BUFFIN(NFTIN,FIXHD,256,LEN_IO,A)

! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= 256)THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of fixed length header of input dump',  &
     &  A,LEN_IO,256)
! DEPENDS ON: ereport
      CALL EREPORT('MAIN_CONVIEEE', 1003,                               &
     & 'Buffer in of fixed length header of input dump wrong size')

      ENDIF

! Set missing data indicator to zero
      DO  I=1,256
        IF(FIXHD(I) <  0)FIXHD(I)=0
      ENDDO

! Input file dimensions
      LEN_FIXHD=256
      LEN_INTHD=FIXHD(101)
      LEN_REALHD=FIXHD(106)
      LEN1_LEVDEPC=FIXHD(111)
      LEN2_LEVDEPC=FIXHD(112)
      LEN1_ROWDEPC=FIXHD(116)
      LEN2_ROWDEPC=FIXHD(117)
      LEN1_COLDEPC=FIXHD(121)
      LEN2_COLDEPC=FIXHD(122)
      LEN1_FLDDEPC=FIXHD(126)
      LEN2_FLDDEPC=FIXHD(127)
      LEN_EXTCNST=FIXHD(131)
      LEN_DUMPHIST=FIXHD(136)
      LEN_CFI1=FIXHD(141)
      LEN_CFI2=FIXHD(143)
      LEN_CFI3=FIXHD(145)
      LEN1_LOOKUP=FIXHD(151)
      LEN2_LOOKUP=FIXHD(152)
      LEN_DATA=FIXHD(161)


!L 3. Buffer in integer constants from dump

! DEPENDS ON: buffin
       CALL BUFFIN(NFTIN,INTHD,FIXHD(101),LEN_IO,A)

! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= FIXHD(101))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of integer constants in input dump',    &
     &  A,LEN_IO,FIXHD(101))
! DEPENDS ON: ereport
      CALL EREPORT('MAIN_CONVIEEE', 1004,                               &
     & 'Buffer in of integer constants in input dump wrong size')

      ENDIF

! Set missing data indicator to zero
      DO  I=1,FIXHD(101)
        IF(INTHD(I) <  0)INTHD(I)=0
      ENDDO

       ROW_LENGTH=INTHD(6)
       P_ROWS=INTHD(7)
       P_FIELD=ROW_LENGTH*P_ROWS

!L Extract maximum field size from LOOKUP header
! DEPENDS ON: find_max_field_size
      CALL FIND_MAX_FIELD_SIZE(                                         &
     &      NFTIN,FIXHD(151),FIXHD(152),FIXHD,MAX_FIELD_SIZE,           &
     &      wgdos_expand)
! Rewind file
! DEPENDS ON: setpos
      CALL SETPOS(NFTIN,0,ICODE)

!L 4. Call CONVIEEE

! DEPENDS ON: convieee
      CALL CONVIEEE(LEN_FIXHD,LEN_INTHD,LEN_REALHD,                     &
     &  LEN1_LEVDEPC,LEN2_LEVDEPC,LEN1_ROWDEPC,                         &
     &  LEN2_ROWDEPC,LEN1_COLDEPC,LEN2_COLDEPC,                         &
     &  LEN1_FLDDEPC,LEN2_FLDDEPC,LEN_EXTCNST,                          &
     &  LEN_DUMPHIST,LEN_CFI1,LEN_CFI2,LEN_CFI3,                        &
     &  LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA,P_FIELD,                       &
     &  P_ROWS,ROW_LENGTH,                                              &
     &  NFTIN,NFTOUT,IEEE_TYPE,                                         &
     &  MAX_FIELD_SIZE, WGDOS_EXPAND)

      call gc_exit()

      STOP
      END PROGRAM MAIN_CONVIEEE
!*L  Arguments:-------------------------------------------------------


!----------------------------------------------------------------------

#endif
