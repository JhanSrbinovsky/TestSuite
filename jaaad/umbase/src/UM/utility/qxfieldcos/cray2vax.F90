#if defined(FLDC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: FIELDCOS -------------------------------------------------
!LL
!LL  Purpose:
!LL  To read a model dump or direct access fieldsfile and convert it to
!LL  a sequential PP file ready for transfer to a different platform.
!LL
!LL   A general note on fieldcos -- When doing a bit compare on
!LL   the output of fieldcos half words may disagree. This is caused
!LL   by the extra half word after an odd number of words in a field
!LL   and is nothing to worry about. (Simon Tett 13/5/92)
!LL
!LL       16/10/92 Added routines for conversion to VAX or IEEE
!LL  data formats, changes to LBPACK/LBUSER1 codes.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  3.1   19/02/93 Use FIXHD(12) not FIXHD(1) as Version no in P21BITS
!LL  3.1   29/01/93 Reset LBLREC when unpacking data
!LL  3.2   25/03/93 use COMDECK CHSUNITS for size of FLAG_IO
!LL  3.2   31/03/93 check dumps indicator in fixed header,
!LL                 correct data lengths for model dump conversions
!LL                 correct INTENT comments for subroutine arguments
!LL                 add fix to put correct m08 code on max/min temps
!LL                 Correct OPEN statement for UNICOS 7.0     PJS
!LL                 Code for real missing data indicator from PPHEADER
!LL  3.3   08/02/94 Modify calls to TIME2SEC/SEC2TIME to output/input
!LL                 elapsed times in days & secs, for portability, TCJ
!LL                 and correct day number when oper=.true.        RR
!LL  3.3   19/04/94 Check and correct invalid LBREL
!LL                 correct error in END_SECOND usage. P.Smith
!LL  3.4   18/05/94 Add processing of Logical data
!LL                 with fix for Land/Sea mask         P.Smith
!LL  3.4   09/09/94 Add GRIB decoder                   P.Smith
!LL  3.4   17/06/94 *CALL CCONTROL inserted (declares logical switches
!LL                  which replace *DEFs - LCAL360 replaces CAL360)
!LL                 Argument LCAL360 passed to S/R's READ_WRITE,
!LL                  CRAY_IBM, CRAY_VAX, CRAY_IBM and passed on to
!LL                  S/R's SEC2TIM, TIME2SEC
!LL                                               S.J.Swarbrick
!    4.0   30/03/95 Add new format option - GRIB to strip grib output
!                   from the model of its pp headers and output as
!                   pure binary grib. Also allow conversion of stash
!                   codes to standard grib code table 2 values or
!                   a user set of codes. R A Stratton
!LL  3.5  13/06/95  Remove comdeck CCONTROL and replace with locally
!LL                 declared LCAL360.                  RTHBarnes.
!    4.2   25/02/97  In order to remove the need for "assign" in the
!    calling script the C I/O routines GET_FILE and FILE_OPEN are used
!    in place of the FORTRAN OPEN statement. This results in a calling
!    script with unit declarations i.e.
!    export UNIT07="Diagnostic filename"
!    export UNIT10="Input filename"
!    export UNIT11="Output filename"
!    Also data conversion routines CRAY2IBM and CRAY2IEG changed to
!    Cray IEEE CRI2IBM and CRI2IEG conversion routines.    Ian Edmond
!    4.3 17/4/97 Cray 32 unpacking functionality added again  IEdmond
!    4.4 17/7/97 Fix to subroutine READFF to read wfio dumpfiles. IE
!LL   4.5    18/09/98  Corrected non-standard FORMAT statments
!LL                                                  P.Burton
!LL  5.0  06/05/99  Pass RoutineName to EREPORT. M Gallani
!LL  5.2  7/11/99  Enable run length encoding of ocean fieldsfiles to
!LL                compress the sequences of mdi values that represent
!LL                the land points. Ian Edmond
!LL  5.2  15/12/00  A check in initial data time needed to recalculate
!LL                 data time for accum/time mean fields is added.
!LL                 E.Leung
!LL  5.4  17/04/02  Allow all currently valid extra data
!LL                 vector types and simplify validation. R. Hill
!LL  5.4  02/05/02  Updated extra-data type and improved comments
!LL                 E.Leung
!LL  5.5  25/04/03  Grib data format not supported on non-CRAY
!LL                 platform                           E.Leung
!LL  5.5  28/02/03  Insert code for portable data conversion routines
!LL                 to replace Cray-specific CRI2IBM etc.     P.Dando
!    6.0  10/09/03  Conversion of portable data conversion routines
!                   (IEEE2IBM etc) into functions with error return
!                   codes matching those of CRAY routines.    P.Dando
!LL  6.0  02/07/03  Move subroutine CHECK_EXTRA and function
!LL                 INT_FROM_REAL in deck FIELDCOS to deck PREXTRA
!LL                 (minimize exec_xref decks on NEC)     E.Leung
!LL  6.0  18/06/03  Grib data format not supported     E.Leung
!LL  6.0  21/01/04  Call to DEGRIB put back in but NECSX6 DEF
!LL                 added around it as libgrib was then only
!LL                 available on NECSX6 W Roseblade
!    6.2  10/04/05  Removed calls to ABORT.  J. Gill
!LL
!    6.0  11/09/03  Replaced ABORT call with call to ABORT_IO. P.Dando
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: C41
!LL
!LL  Project task: C4
!LL
!LL  External documentation: UM documentation paper Y8
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
!LL  Routine: READ_WRITE -----------------------------------------------
!LL
!LL  Purpose: To read a   direct access PP file  and convert it to a
!LL  sequential file read to be passed across to the IBM
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: C41
!LL
!LL  Project task: C4
!LL
!LL  External documentation: UM Documentation paper C4
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
!LL  Routine: CRAY_IBM-------------------------------------------------
!LL
!LL  Purpose: To read a   direct access PP file  and convert it to a
!LL  sequential file read to be passed across to the IBM
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  4.5 24/7/98 Change to output land sea mask as a real field (as a
!LL              special case). Rick Rawlins
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: C41
!LL
!LL  Project task: C4
!LL
!LL  External documentation:
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!

!LL  Routine: CRAY_VAX-------------------------------------------------
!LL
!LL  Purpose: To read a   direct access PP file  and convert it to a
!LL  sequential file read to be passed in VAX format
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version:
!LL
!LL  Author:   P.J .Smith         Date: 26 June  1992
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: C41
!LL
!LL  Project task: C4
!LL
!LL  External documentation:
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!

!LL  Routine: CRAY_IEEE-------------------------------------------------
!LL
!LL  Purpose: To read a   direct access PP file  and convert it to a
!LL  sequential file read to be passed in IEEE format
!LL
!LL  Author:   P.J .Smith         Date: 26 June  1992
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version:
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: C41
!LL
!LL  Project task: C4
!LL
!LL  External documentation:
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!

#if defined(CRAY)
!=======================================================================
!    Routine: CRAY_GRIB
!
!    Purpose: To read a   direct access PP file  and convert it to a
!    pure grib file ready to be passed to HDS or workstation.
!
!    Tested under compiler:   cft77
!    Tested under OS version: UNICOS 7 & 8
!
!    Model            Modification history from model version 3.3:
!   version  Date
!    4.0    31/03/95  : Added to FIELDCOS
!
!    Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
!    Logical components covered: C41
!
!    Project task: C4
!
!    External documentation:
!
!-----------------------------------------------------------------------
!    Interface and arguments: ------------------------------------------
!
! =====================================================================
#endif
!LL  Routine: READFF---------------------------------------------------
!LL
!LL  Purpose: To read a   direct access PP file.
!LL
!LL  Author:   P.Trevelyan
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL  3.4   29/06/94 Correct unpacking of 32 bit data. Error affected
!LL                 odd-length fields. PP and Stash codes added to
!LL                 output. D. Robinson
!LL  6.2   22/08/05 Remove redundent space. P.Selwood
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered:
!LL
!LL  Project task:
!LL
!LL  External documentation:
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------

!LL  Routine: READ_REC--------------------------------------------------
!LL
!LL  Purpose: To read a data record from a  pp file
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: ...
!LL
!LL  Project task: ...
!LL
!LL  External documentation:
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
!LL  Routine: UN_PACK  -------------------------------------------------
!LL
!LL  Purpose: To unpack data from the input array FIELD and return
!LL  the data in FIELD.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered:
!LL
!LL  Project task:
!LL
!LL  External documentation:
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!LL  Routine: LOGICAL_TO_REAL ------------------------------------------
!LL
!LL  Purpose: To convert logical data within FIELD to real data.
!LL  the data in FIELD.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered:
!LL
!LL  Project task:
!LL
!LL  External documentation:
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!LL  Routine: INTEGER_TO_REAL ------------------------------------------
!LL
!LL  Purpose: To convert logical data within FIELD to real data.
!LL  the data in FIELD.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Model            Modification history:
!LL version  Date
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered:
!LL
!LL  Project task:
!LL
!LL  External documentation:
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
#endif


      INTEGER FUNCTION CRAY2VAX(I,LEN_RLABEL,                           &
     &                          VAX_LABEL,BIT_OFF,RLABEL)
      INTEGER I,LEN_RLABEL,BIT_OFF,VAX_LABEL
      REAL RLABEL
      CRAY2VAX=-1
      RETURN
      END FUNCTION CRAY2VAX
