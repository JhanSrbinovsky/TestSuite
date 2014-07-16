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
      PROGRAM MAIN_COMPARE

      IMPLICIT NONE

!  Global Variables:
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "c_mdi.h"
#include "cstash.h"
      INTEGER                                                           &
     & FIXHD1(256)                                                      &
                         !Space for fixed length header file 1
     &,INTHD1(100)       !Space for integer header file 1

      INTEGER                                                           &
     & FIXHD2(256)                                                      &
                         !Space for fixed length header file 2
     &,INTHD2(100)       !Space for integer header file 2

      INTEGER                                                           &
     & LEN_FIXHD1                                                       &
                      !Length of fixed length header on file 1
     &,LEN_INTHD1                                                       &
                      !Length of integer header on file 1
     &,LEN_REALHD1                                                      &
                      !Length of real header on file 1
     &,LEN1_LEVDEPC1                                                    &
                      !1st dim of lev dependent consts on file 1
     &,LEN2_LEVDEPC1                                                    &
                      !2nd dim of lev dependent consts on file 1
     &,LEN1_ROWDEPC1                                                    &
                      !1st dim of row dependent consts on file 1
     &,LEN2_ROWDEPC1                                                    &
                      !2nd dim of row dependent consts on file 1
     &,LEN1_COLDEPC1                                                    &
                      !1st dim of col dependent consts on file 1
     &,LEN2_COLDEPC1                                                    &
                      !2nd dim of col dependent consts on file 1
     &,LEN1_FLDDEPC1                                                    &
                      !1st dim of field dependent consts on file 1
     &,LEN2_FLDDEPC1                                                    &
                      !2nd dim of field dependent consts on file 1
     &,LEN_EXTCNST1                                                     &
                      !Length of extra consts on file 1
     &,LEN_DUMPHIST1                                                    &
                      !Length of history header on file 1
     &,LEN_CFI11                                                        &
                      !Length of index1 on file 1
     &,LEN_CFI21                                                        &
                      !Length of index2 on file 1
     &,LEN_CFI31                                                        &
                      !Length of index3 on file 1
     &,LEN1_LOOKUP1                                                     &
                      !1st dim of LOOKUP on file 1
     &,LEN2_LOOKUP1                                                     &
                      !2nd dim of LOOKUP on file 1
     &,LEN_DATA1                                                        &
                      !Length of data on file 1
     &,ROW_LENGTH1                                                      &
                      !No of points E-W on file 1
     &,P_ROWS1                                                          &
                      !No of p-rows on file 1
     &,P_FIELD1                                                         &
                      !No of p-points per level on file 1
     &,MAX_FIELD_SIZE1 !Maximum field size on file 1

      INTEGER                                                           &
     & LEN_FIXHD2                                                       &
                      !Length of fixed length header on file 2
     &,LEN_INTHD2                                                       &
                      !Length of integer header on file 2
     &,LEN_REALHD2                                                      &
                      !Length of real header on file 2
     &,LEN1_LEVDEPC2                                                    &
                      !1st dim of lev dependent consts on file 2
     &,LEN2_LEVDEPC2                                                    &
                      !2nd dim of lev dependent consts on file 2
     &,LEN1_ROWDEPC2                                                    &
                      !1st dim of row dependent consts on file 2
     &,LEN2_ROWDEPC2                                                    &
                      !2nd dim of row dependent consts on file 2
     &,LEN1_COLDEPC2                                                    &
                      !1st dim of col dependent consts on file 2
     &,LEN2_COLDEPC2                                                    &
                      !2nd dim of col dependent consts on file 2
     &,LEN1_FLDDEPC2                                                    &
                      !1st dim of field dependent consts on file 2
     &,LEN2_FLDDEPC2                                                    &
                      !2nd dim of field dependent consts on file 2
     &,LEN_EXTCNST2                                                     &
                      !Length of extra consts on file 2
     &,LEN_DUMPHIST2                                                    &
                      !Length of history header on file 2
     &,LEN_CFI12                                                        &
                      !Length of index1 on file 2
     &,LEN_CFI22                                                        &
                      !Length of index2 on file 2
     &,LEN_CFI32                                                        &
                      !Length of index3 on file 2
     &,LEN1_LOOKUP2                                                     &
                      !1st dim of LOOKUP on file 2
     &,LEN2_LOOKUP2                                                     &
                      !2nd dim of LOOKUP on file 2
     &,LEN_DATA2                                                        &
                      !Length of data on file 2
     &,ROW_LENGTH2                                                      &
                      !No of points E-W on file 2
     &,P_ROWS2                                                          &
                      !No of p-rows on file 2
     &,P_FIELD2                                                         &
                      !No of p-points per level on file 2
     &,MAX_FIELD_SIZE2 !Maximum field size on file 2


      INTEGER                                                           &
     & LEN_IO                                                           &
                !Length of I/O returned by BUFFER IN
     &,I                                                                &
                !Loop index
     &,NFTIN1                                                           &
                !Unit number of input UM file 1
     &,NFTIN2                                                           &
                !Unit number of input UM file 2

     &,ERR                                                              &
                !Return code from OPEN
     &,ICODE    !Return code from setpos
      REAL A    !BUFFER IN UNIT function
!
      integer expand
      CHARACTER*100 PAREXE_ENV  ! hold name of the // exec script
      INTEGER ME_GC,NPROC_GC

      Integer, Parameter :: Max_Filename_Len = 150
      Character (Len=Max_Filename_Len) :: FileName1
      Character (Len=Max_Filename_Len) :: FileName2


! External subroutines called:------------------------------------------
      EXTERNAL IOERROR,ABORT_IO,BUFFIN,FILE_OPEN,                       &
     &         SETPOS,ABORT,COMPARE,InitPrintStatus,GC_INIT
!*----------------------------------------------------------------------

      PAREXE_ENV=' '

      CALL GC_INIT(PAREXE_ENV,ME_GC,NPROC_GC)
! Initialise print status for standard output
! DEPENDS ON: initprintstatus
      CALL InitPrintStatus
!
      expand=1
!L 1. Assign unit numbers

      NFTIN1=20
      NFTIN2=21

      WRITE(6,*)' COMPARE - FULL MODE'
      WRITE(6,*)' -------------------'
      WRITE(6,*)' '

      Call Fort_Get_Env ('FILE1',5,FileName1,Max_Filename_Len,Err)
      Call Fort_Get_Env ('FILE2',5,FileName2,Max_Filename_Len,Err)

      write (6,*) ' Files being compared '
      write (6,*) ' -------------------- '
      write (6,*) ' File 1 : ',FileName1( 1:len_trim(FileName1) )
      write (6,*) ' File 2 : ',FileName2( 1:len_trim(FileName2) )
      write (6,*) ' '

      WRITE(6,'(20x,''FILE STATUS'')')
      WRITE(6,'(20x,''==========='')')
!     CALL OPEN(1,'PPXREF',6,0,0,ERR)
! DEPENDS ON: file_open
      CALL FILE_OPEN(NFTIN1,'FILE1',5,0,0,ERR)
! DEPENDS ON: file_open
      CALL FILE_OPEN(NFTIN2,'FILE2',5,0,0,ERR)

!L 2. Buffer in fixed length header record from file 1

! DEPENDS ON: buffin
      CALL BUFFIN(NFTIN1,FIXHD1,256,LEN_IO,A)

! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= 256)THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of fixed length header of input file',  &
     &  A,LEN_IO,256)
! DEPENDS ON: ereport
      CALL EREPORT('MAIN_COMPARE', 1000,                                &
     & 'Buffer in of fixed length header wrong size')

      ENDIF

! Set missing data indicator to zero
      DO  I=1,256
        IF(FIXHD1(I) <  0)FIXHD1(I)=0
      ENDDO

! Input file dimensions
      LEN_FIXHD1=256
      LEN_INTHD1=FIXHD1(101)
      LEN_REALHD1=FIXHD1(106)
      LEN1_LEVDEPC1=FIXHD1(111)
      LEN2_LEVDEPC1=FIXHD1(112)
      LEN1_ROWDEPC1=FIXHD1(116)
      LEN2_ROWDEPC1=FIXHD1(117)
      LEN1_COLDEPC1=FIXHD1(121)
      LEN2_COLDEPC1=FIXHD1(122)
      LEN1_FLDDEPC1=FIXHD1(126)
      LEN2_FLDDEPC1=FIXHD1(127)
      LEN_EXTCNST1=FIXHD1(131)
      LEN_DUMPHIST1=FIXHD1(136)
      LEN_CFI11=FIXHD1(141)
      LEN_CFI21=FIXHD1(143)
      LEN_CFI31=FIXHD1(145)
      LEN1_LOOKUP1=FIXHD1(151)
      LEN2_LOOKUP1=FIXHD1(152)
      LEN_DATA1=FIXHD1(161)

!L 3. Buffer in fixed length header record from file 2

! DEPENDS ON: buffin
      CALL BUFFIN(NFTIN2,FIXHD2,256,LEN_IO,A)

! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= 256)THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of fixed length header of input file',  &
     &  A,LEN_IO,256)
! DEPENDS ON: ereport
      CALL EREPORT('MAIN_COMPARE', 1001,                                &
     & 'Buffer in of fixed length header wrong size')

      ENDIF

! Set missing data indicator to zero
      DO  I=1,256
        IF(FIXHD2(I) <  0)FIXHD2(I)=0
      ENDDO

! Input file dimensions
      LEN_FIXHD2=256
      LEN_INTHD2=FIXHD2(101)
      LEN_REALHD2=FIXHD2(106)
      LEN1_LEVDEPC2=FIXHD2(111)
      LEN2_LEVDEPC2=FIXHD2(112)
      LEN1_ROWDEPC2=FIXHD2(116)
      LEN2_ROWDEPC2=FIXHD2(117)
      LEN1_COLDEPC2=FIXHD2(121)
      LEN2_COLDEPC2=FIXHD2(122)
      LEN1_FLDDEPC2=FIXHD2(126)
      LEN2_FLDDEPC2=FIXHD2(127)
      LEN_EXTCNST2=FIXHD2(131)
      LEN_DUMPHIST2=FIXHD2(136)
      LEN_CFI12=FIXHD2(141)
      LEN_CFI22=FIXHD2(143)
      LEN_CFI32=FIXHD2(145)
      LEN1_LOOKUP2=FIXHD2(151)
      LEN2_LOOKUP2=FIXHD2(152)
      LEN_DATA2=FIXHD2(161)


!L 4. Buffer in integer constants from file 1

! DEPENDS ON: buffin
       CALL BUFFIN(NFTIN1,INTHD1,FIXHD1(101),LEN_IO,A)

! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= FIXHD1(101))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of integer constants in input file 1',  &
     &  A,LEN_IO,FIXHD1(101))
! DEPENDS ON: ereport
      CALL EREPORT('MAIN_COMPARE', 1002,                                &
     & 'Buffer in of integer constants wrong length')

      ENDIF

! Set missing data indicator to zero
      DO  I=1,FIXHD1(101)
        IF(INTHD1(I) <  0)INTHD1(I)=0
      ENDDO

       ROW_LENGTH1=INTHD1(6)
       P_ROWS1=INTHD1(7)
       P_FIELD1=ROW_LENGTH1*P_ROWS1

!L Extract maximum field size from LOOKUP header
! DEPENDS ON: find_max_field_size
      CALL FIND_MAX_FIELD_SIZE(NFTIN1,FIXHD1(151),FIXHD1(152),FIXHD1    &
     &    ,max_field_size1, expand)

!L 5. Buffer in integer constants from file 2

! DEPENDS ON: buffin
       CALL BUFFIN(NFTIN2,INTHD2,FIXHD2(101),LEN_IO,A)

! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= FIXHD2(101))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of integer constants in input file 2',  &
     &  A,LEN_IO,FIXHD2(101))
! DEPENDS ON: ereport
      CALL EREPORT('MAIN_COMPARE', 1003,                                &
     & 'Buffer integer constants wrong length')

      ENDIF

! Set missing data indicator to zero
      DO  I=1,FIXHD2(101)
        IF(INTHD2(I) <  0)INTHD2(I)=0
      ENDDO

!L 6. Cause abort if files obviously different

      ROW_LENGTH2=INTHD2(6)
      P_ROWS2=INTHD2(7)
      P_FIELD2=ROW_LENGTH2*P_ROWS2

!L Extract maximum field size from LOOKUP header
! DEPENDS ON: find_max_field_size
      CALL FIND_MAX_FIELD_SIZE(NFTIN2,FIXHD2(151),FIXHD2(152),FIXHD2    &
     &    ,max_field_size2, expand)

      IF(P_FIELD1 /= P_FIELD2)THEN
       WRITE(6,*)'COMPARE: ERROR Dumps are at different resolutions'
! DEPENDS ON: ereport
       CALL EREPORT('MAIN_COMPARE', 1004,                               &
     &  'Dumps are at difference resolutions')

      ENDIF
      IF(LEN2_LOOKUP1 /= LEN2_LOOKUP2)THEN
       WRITE(6,*)                                                       &
     & 'COMPARE: WARNING Dumps have different number of fields'
      ENDIF

! Rewind files
! DEPENDS ON: setpos
      CALL SETPOS(NFTIN1,0,ICODE)
! DEPENDS ON: setpos
      CALL SETPOS(NFTIN2,0,ICODE)

!L 7. Call COMPARE

! DEPENDS ON: compare
      CALL COMPARE(LEN_FIXHD1,LEN_INTHD1,LEN_REALHD1,                   &
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

      call gc_exit()

      STOP
      END PROGRAM MAIN_COMPARE
!*L  Arguments:-------------------------------------------------------

#endif
