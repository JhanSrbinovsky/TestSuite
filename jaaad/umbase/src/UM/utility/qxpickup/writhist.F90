#if defined(PICK)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: WRITHIST
!LL
!LL  Purpose: To prepend the contents of History common blocks to
!LL           the beginning of the history file input.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.0
!LL
!LL  Author:   A.Sangster
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL
!LL   3.1  05/02/93    Portable Fortran unit no assigns
!LL                    Author: A. Dickinson    Reviewer: R. Stratton
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL  3.5  06/04/95  Sub-Models stage 1: revise History and Control file
!LL                 contents.  RTHBarnes.
!LL  4.0  18/10/95  Add ICODE error return to GET_FILE call. RTHBarnes
!LL  4.1  01/05/96  Correct *IF DEF - remove CONTROL & COMB. RTHBarnes
!LL  4.5  05/06/98  Add DELIM='APOSTROPHE' to OPEN statement so that
!LL                 History file is written in correct form from
!LL                 Fujitsu.                       RBarnes@ecmwf.int
!    5.5  20/02/03  Add DELIM='APOSTROPHE' to OPEN statement for LINUX
!                   and NEC platforms.                         P.Dando
!    6.0  11/09/03  Add DELIM='APOSTROPHE' to OPEN statement for IBM
!                   and SGI platforms.                         P.Dando
!    6.1  18/08/04  Add DELIM='APOSTROPHE' in open statement for ALTIX
!                   (supplied by Z. Chaplin).                  P.Dando
!    6.2  31/01/06   Changed defs from specific compilers to a general
!                    GNU/Linux option. T. Edwards
!    6.2  16/05/05  Add DELIM=APOSTROPHE in open statement for HP.
!                   Simon Wilson
!LL
!LL  Programming standard: UM Doc Paper 3, draft version 3 (15/1/90)
!LL
!LL  Logical components covered: H4
!LL
!LL  Project task: H
!LL
!LL  Documentation:  Unified Model Documentation Paper
!LL                  H- History Bricks
!LL
!
!*L  Interface and arguments:
!
      SUBROUTINE WRITHIST                                               &
     &         ( UNITHIST,UNITCOPY,ICODE,CMESSAGE )
!
      IMPLICIT NONE
!
      INTEGER   UNITHIST  ! In  - History file unit
      INTEGER   UNITCOPY  ! In  - unit no. for copy of old history file
      INTEGER       ICODE    ! Out - Return code from routine
      CHARACTER*80  CMESSAGE ! Out - Return message if failure occured
!*
!
!L Common blocks
!
#include "csubmodl.h"
#include "chsunits.h"
#include "chistory.h"
!
      CHARACTER*80 FILENAME
      CHARACTER*8  NLNAME
!*L EXTERNAL subroutines called
      EXTERNAL GET_FILE
!*
!L
!L 1. Open history file and rewind
!L
      CALL GET_FILE(UNITHIST,FILENAME,80,ICODE)
#if defined(LINUX) \
 || defined(NEC) || defined(IBM) || defined(SGI) || defined(ALTIX)   \
 || defined(HP)
      OPEN(UNITHIST,FILE=FILENAME,FORM='FORMATTED',                     &
     &  DELIM='APOSTROPHE',IOSTAT=ICODE)
#else
      OPEN(UNITHIST,FILE=FILENAME,FORM='FORMATTED',IOSTAT=ICODE)
#endif
!
! Check for error
!
      IF(ICODE  >  0)THEN
        CMESSAGE='WRITHIST: Failed in OPEN of permanent history file'
        GOTO 999
      ELSEIF(ICODE  <   0)THEN
      WRITE(6,*)'WRITHIST: Warning message on OPEN of permanent history &
     & file'
        WRITE(6,*)'IOSTAT= ',ICODE
      ENDIF
!
      REWIND(UNITHIST)
!L
!L 2. Write new record as first record of history file
!L
      NLNAME = 'NLIHISTO'
      WRITE(UNITHIST,NLIHISTO,ERR=200)
      NLNAME = 'NLCHISTO'
      WRITE(UNITHIST,NLCHISTO,ERR=200)
      NLNAME = 'NLIHISTG'
      WRITE(UNITHIST,NLIHISTG,ERR=200)
      NLNAME = 'NLCHISTG'
      WRITE(UNITHIST,NLCHISTG,ERR=200)
      NLNAME = 'NLCFILES'
      WRITE(UNITHIST,NLCFILES,ERR=200)
      go to 199
!
! Check for error
!
! Write error
  200 continue
      ICODE = 2
      CMESSAGE='WRITHIST: Write ERROR on history file for namelist'//   &
     & NLNAME
      go to 999
  199 continue
!L
!L 2. Open copy of old history file and rewind
!L
      CALL GET_FILE(UNITCOPY,FILENAME,80,ICODE)
#if defined(LINUX) \
 || defined(NEC) || defined(IBM) || defined(SGI) || defined(ALTIX)
       OPEN(UNITCOPY,FILE=FILENAME,FORM='FORMATTED',                    &
     &   DELIM='APOSTROPHE',IOSTAT=ICODE)
#else
       OPEN(UNITCOPY,FILE=FILENAME,FORM='FORMATTED',IOSTAT=ICODE)
#endif
!
! Check for error
!
      IF(ICODE  >  0)THEN
        CMESSAGE='WRITHIST: Failed in OPEN of copy of old history file'
        GOTO 999
      ELSEIF(ICODE  <   0)THEN
      WRITE(6,*)'WRITHIST:Warning message on OPEN of copy of old history&
     &y file'
        WRITE(6,*)'IOSTAT= ',ICODE
      ENDIF
!
      REWIND(UNITCOPY)
!L
!L 3. Read each record of old and write to new permanent history file
!L
  250 continue
      NLNAME = 'NLIHISTO'
      READ (UNITCOPY,NLIHISTO,END=100,ERR=300)
      WRITE(UNITHIST,NLIHISTO,ERR=400)
      NLNAME = 'NLCHISTO'
      READ (UNITCOPY,NLCHISTO,END=100,ERR=300)
      WRITE(UNITHIST,NLCHISTO,ERR=400)
      NLNAME = 'NLIHISTG'
      READ (UNITCOPY,NLIHISTG,END=100,ERR=300)
      WRITE(UNITHIST,NLIHISTG,ERR=400)
      NLNAME = 'NLCHISTG'
      READ (UNITCOPY,NLCHISTG,END=100,ERR=300)
      WRITE(UNITHIST,NLCHISTG,ERR=400)
      NLNAME = 'NLCFILES'
      READ (UNITCOPY,NLCFILES,END=100,ERR=300)
      WRITE(UNITHIST,NLCFILES,ERR=400)
      go to 250
!
! Check for error
!
! End-of-file
  100 continue
      IF (NLNAME  ==  'NLIHISTO')  THEN   ! expected end-of-file
      WRITE(6,*)'Copied old history records to new phist file completed'
      go to 999
      ELSE   ! unexpected end-of-file
      ICODE = 1
      CMESSAGE='WRITHIST: End of file in READ from history file for name&
     &list '//NLNAME
      go to 999
      END IF
! Read error
  300 continue
      ICODE = 3
      CMESSAGE='WRITHIST: Read ERROR on history file for namelist'//    &
     & NLNAME
      go to 999
! Write error
  400 continue
      ICODE = 4
      CMESSAGE='WRITHIST: Write ERROR on history file for namelist'//   &
     & NLNAME
!
 999  CONTINUE
!L
!L 3. Close and return
!L
      CLOSE(UNITHIST)
      RETURN
      END SUBROUTINE WRITHIST
#endif
