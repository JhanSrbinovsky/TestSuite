#if defined(CONVPP)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  PROGRAM MAIN_CONVPP --------------------------------------------
!LL
!LL  Purpose: Converts a UM file into PP format.
!LL
!LL  Written by A. Dickinson 05/07/93
!LL
!LL  Model            Modification history:
!LL version  Date
!LL
!LL  3.3   31/10/93   Dimension of data array set to maximum value
!LL                   Author: A. Dickinson      Reviewer: P.Burton
!LL
!LL   3.3   15/12/93  Rename subroutine PRINTDUMP to CONVPP. D Robinson
!LL
!LL   3.3   08/12/93  Extra argument for READFLDS. D. Robinson
!LL
!LL   3.4   23/09/94  Extended to process ocean dumps. Alternative
!LL                   subroutine introduced
!LL                   Author D.M.Goddard
!LL   3.5  24/03/95    Changed OPEN to FILE_OPEN  P.Burton
!LL   4.4  24/10/96   Ocean data is written out without wrap points.
!LL                   Catherine Jones
!LL   4.4  23/04/97   Compressed fields are uncompressed using the
!LL                   subroutine UNPACK
!LL                   Catherine Jones


!LL   4.2  Oct. 96    DEF CRAY replaced by DEF T3E
!LL                             S.J.Swarbrick
!LL  4.4   Oct. 1997 Changed error handling from routine HDPPXRF
!LL                  so only fatal (+ve) errors are handled.
!LL                                             Shaun de Witt
!   4.4  23/04/97   Corrections to processing of Land-sea mask and
!                   Land compressed fields
!                   D.M. Goddard
!   4.4  24/10/97   Initialise ICODE as it is no longer
!                   initialised in HDPPXRF
!                   Author D.M. Goddard
!     4.5  14/10/97   Sets most significant number in packing indicator
!                     to zero (Native format) to enable PP package to be
!                     used on fieldsfile output
!                     Author D.M. Goddard
!     5.1  11/05/00   Added call to InitPrintStatus. D.P.Matthews
!LL  5.2  7/11/99  Enable run length encoding of ocean fieldsfiles to
!LL                compress the sequences of mdi values that represent
!LL                the land points. Ian Edmond
!    5.3  22/11/01  Enable MPP as the only option for
!                   small executables         E.Leung
!    5.4  03/09/02  Arguments are added to DECOMPOSE_SMEXE
!                                                  E.Leung
!    5.5  31/01/03  - Use INTHD 6,7 instead of LOOKUP 18,19 as GLSIZE
!                   - Enable convpp to be built on HP
!                                                             E.Leung
!    6.0  11/09/03  Add call to gc_exit                       P.Dando
!    6.2  10/04/05  Removed calls to ABORT.  J. Gill
!LL
!LL  Programming standards:
!LL
!LL  Logical components covered:
!LL
!LL  System Tasks: F3,F4,F6
!LL
!LL  Documentation: UM Doc Paper F5
!LL
!LL  -----------------------------------------------------------------
      PROGRAM MAIN_CONVPP

      IMPLICIT NONE


      CHARACTER*80 ARG1,ARG2  ! Filenames


      INTEGER                                                           &
     & FIXHD(256)                                                       &
                         !Space for fixed length header
     &,INTHD(100)        !Space for integer header

      INTEGER                                                           &
     & LEN_FIXHD                                                        &
                      !Length of fixed length header on input file
     &,LEN_INTHD                                                        &
                      !Length of integer header on input file
     &,JOC_NO_SEAPTS                                                    &
                      !Number of points in compressed array
     &,LEN_OCFLD                                                        &
                      !Length of uncompressed ocean field
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
                !Unit number of input UM dump 1
     &,ERR                                                              &
                !Return code from OPEN

     &,ICODE    !Return code from setpos
      REAL A    !BUFFER IN UNIT function

      CHARACTER*100 PAREXE_ENV  ! hold name of the // exec script
      INTEGER ME_GC,NPROC_GC

! External subroutines called:------------------------------------------
      EXTERNAL IOERROR,ABORT_IO,BUFFIN,FILE_OPEN,SETPOS,ABORT,          &
     &         ATMOS_CONVPP,OCEAN_CONVPP,InitPrintStatus,GC_INIT
!*----------------------------------------------------------------------

      PAREXE_ENV=' '

      CALL GC_INIT(PAREXE_ENV,ME_GC,NPROC_GC)
! Initialise print status for standard output
! DEPENDS ON: initprintstatus
      CALL InitPrintStatus

!L 1. Assign unit numbers

      NFTIN=20

      WRITE(6,'(20x,''FILE STATUS'')')
      WRITE(6,'(20x,''==========='')')

! DEPENDS ON: file_open
      CALL FILE_OPEN(20,'FILE1',5,0,0,ERR)
      CALL GET_FILE(10,ARG2,80,ICODE)
      OPEN(10,FILE=ARG2,FORM='UNFORMATTED')

!L 2. Buffer in fixed length header record

! DEPENDS ON: buffin
      CALL BUFFIN(NFTIN,FIXHD,256,LEN_IO,A)

! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= 256)THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of fixed length header of input dump',  &
     &  A,LEN_IO,256)
! DEPENDS ON: ereport
      CALL EREPORT('MAIN_CONVPP', 1000,                                 &
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
      CALL EREPORT('MAIN_CONVPP', 1001,                                 &
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
       CALL FIND_MAX_FIELD_SIZE(                                        &
     &      NFTIN,FIXHD(151),FIXHD(152),FIXHD,MAX_FIELD_SIZE)

! Calculate sizes of compressed and uncompressed ocean fields
      JOC_NO_SEAPTS=INTHD(11)
      IF(FIXHD(2) == 2)THEN
        LEN_OCFLD=INTHD(6)*INTHD(7)*INTHD(8)
      ELSE
        LEN_OCFLD=0
      ENDIF
! Rewind file
! DEPENDS ON: setpos
      CALL SETPOS(NFTIN,0,ICODE)

      IF(FIXHD(2) == 1)THEN

! Atmospheric dump
! DEPENDS ON: atmos_convpp
      CALL ATMOS_CONVPP (LEN_FIXHD,LEN_INTHD,LEN_REALHD,                &

     &  LEN1_LEVDEPC,LEN2_LEVDEPC,LEN1_ROWDEPC,                         &
     &  LEN2_ROWDEPC,LEN1_COLDEPC,LEN2_COLDEPC,                         &
     &  LEN1_FLDDEPC,LEN2_FLDDEPC,LEN_EXTCNST,                          &
     &  LEN_DUMPHIST,LEN_CFI1,LEN_CFI2,LEN_CFI3,                        &
     &  LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA,P_FIELD,                       &
     &  NFTIN,MAX_FIELD_SIZE)

      ELSEIF (FIXHD(2) == 2)THEN

! Oceanic dump
! DEPENDS ON: ocean_convpp
      CALL OCEAN_CONVPP (LEN_FIXHD,LEN_INTHD,LEN_REALHD,                &
     &  LEN1_LEVDEPC,LEN2_LEVDEPC,LEN1_ROWDEPC,                         &
     &  LEN2_ROWDEPC,LEN1_COLDEPC,LEN2_COLDEPC,                         &
     &  LEN1_FLDDEPC,LEN2_FLDDEPC,LEN_EXTCNST,                          &
     &  LEN_DUMPHIST,LEN_CFI1,LEN_CFI2,LEN_CFI3,                        &
     &  LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA,P_FIELD,                       &
     &  NFTIN,MAX_FIELD_SIZE,JOC_NO_SEAPTS,LEN_OCFLD)
      ENDIF

      call gc_exit()

      STOP
      END PROGRAM MAIN_CONVPP
!LL  SUBROUTINE ATMOS_CONVPP -----------------------------------------
!LL
!LL Purpose: Converts UM file to PP format.
!LL
!LL  Written by A. Dickinson
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL   3.4   23/09/94  New output lookup table array introduced because
!LL                   element 21 set to  0 for output file would need
!LL                   to be reset before attempting to read next record
!LL                   Routine renamed because separate routine for
!LL                   ocean dump introduced.
!LL                   Author D.M.Goddard
!     4.1  18/06/96  Changes to cope with changes in STASH addressing
!                    Author D.M. Goddard.
!LL
!LL  Documentation: UM Doc Paper F5
!LL
!LL  System Tasks: F3,F4,F6
!LL
!LL  -----------------------------------------------------------------
!*L  Arguments:-------------------------------------------------------
!LL  SUBROUTINE OCEAN_CONVPP-----------------------------------------
!LL
!LL Purpose: Converts UM ocean file to PP format.
!LL
!LL  Written by D.M. Goddard
!LL
!LL  Model            Modification history from model version 3.4:
!LL version  Date
!LL
!LL   3.4   23/09/94  New routine at version 3.4
!LL
!LL  Documentation: UM Doc Paper F5
!LL
!LL  System Tasks: F3,F4,F6
!LL
!LL  -----------------------------------------------------------------
!*L  Arguments:-------------------------------------------------------

#endif
