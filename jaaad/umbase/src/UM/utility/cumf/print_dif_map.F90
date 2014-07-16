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

      SUBROUTINE PRINT_DIF_MAP(DIFF,ROWS,COLS,KEY,                      &
     &  GRID_TYPE1,GRID_TYPE2)
! Declarations
!LL Writes out a map of the differences between two fields - with one
!LL character per point. This allows points in two fields which are
!LL different to be quickly identified.
!LL Writes to UNIT10 - opened in COMPARE - filename must be supplied
!LL by UNIT10 environment variable via the cumf script

      IMPLICIT NONE

      INTEGER                                                           &
     &  ROWS                                                            &
              ! IN : number of rows in field
     &, COLS                                                            &
              ! IN : number of cols in field
     &, GRID_TYPE1                                                      &
                      !IN grid type for file 1
     &, GRID_TYPE2    !IN grid type for file 2

      CHARACTER*1                                                       &
     &  DIFF(ROWS*COLS)  ! IN : difference map field to be output

      CHARACTER*(*)                                                     &
     &  KEY  ! IN : key to difference map

! Local variables
      INTEGER X,Y,Z
      integer i ,j                                                      &
     &, LAND_POINTS                                                     &
                      !no of land points
     &, POINTS        !total no of land & sea points to be processed

! Constants from comdecks:----------------------
#include "c_mdi.h"
#include "parvars.h"
#include "atm_lsm.h"
! External subroutines called:--------------------------
      EXTERNAL FROM_LAND_POINTS

      REAL WORK_DIFF(ROWS*COLS),NWORK_DIFF(ROWS*COLS)
      character*1 numb(10), blank

      data numb/'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/
      data blank /' '/

      WRITE(10,'(/a/)') KEY


      POINTS=ROWS*COLS

!L 1. Converting symbols to real numbers (input to FROM_LAND_POINTS)
!  -----------------------------------------------------------------
          IF((GRID_TYPE1 == 21).AND.(GRID_TYPE2 == 21)) THEN
            do I=1,rows*cols
              if (diff(I) == "#") work_diff(I)=5.0
              if (diff(I) == "X") work_diff(I)=4.0
              if (diff(I) == "O") work_diff(I)=3.0
              if (diff(I) == "o") work_diff(I)=2.0
              if (diff(I) == ":") work_diff(I)=1.0
              if (diff(I) == ".") work_diff(I)=0.0
            enddo

! Pass work_diff to FROM_LAND_POINTS
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS(NWORK_DIFF,WORK_DIFF,ATMOS_LANDMASK,  &
     &                            POINTS,LAND_POINTS)

!L 2. Converting real numbers to symbols (to be output)
!  ----------------------------------------------------
            do I=1,points
              if (nwork_diff(I) == 5.0) diff(I)="#"
              if (nwork_diff(I) == 4.0) diff(I)="X"
              if (nwork_diff(I) == 3.0) diff(I)="O"
              if (nwork_diff(I) == 2.0) diff(I)="o"
              if (nwork_diff(I) == 1.0) diff(I)=":"
              if (nwork_diff(I) == 0.0) diff(I)="."
              if (nwork_diff(I) == rmdi) diff(I)="~"
            enddo
          ENDIF
#if defined(LFOK)

      write(10,'(6x,120a1)') ((blank, j=1,9),                           &
     & numb(mod((i+10)/10, 10)+1), i=1, cols, 10)
!
      write(10,'(6x,120a1)') (numb(mod(i, 10)+1), i=1,cols)

      do y=1,rows
        z=(y-1)*cols
        if(cols == 120) then
          write(10,123) y,(diff(x+z),x=1,cols)
 123      format(1x,i3,'->',120a1)
        else
          write(10,124) y,(diff(x+z),x=1,cols)
 124      format(1x,i3,'->',120a1/(6x,120a1))
        endif
      enddo
#else
      WRITE(6,*) 'Difference maps not supported on this platform'
#endif

      RETURN
      END SUBROUTINE PRINT_DIF_MAP
#endif
