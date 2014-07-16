#if defined(CONTROL) || defined(SETUP) || defined(COMB)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Routine : TEMPHIST -------------------------------------------------
!LL
!LL Purpose :Write current contents of history common block to temporary
!LL          or interim history file - overwriting previous record
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.0
!LL
!LL  Author:   A.B.SANGSTER       Date:           20 January 1990
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
!LL  4.5  05/05/98  Add DELIM='APOSTROPHE' to OPEN statement so that
!LL                 History file is written in correct form from
!LL                 Fujitsu.                       RBarnes@ecmwf.int
!    5.4  09/04/02   Need to specify DELIM='APOSTROPHE' in OPEN for
!                    namelists for Linux compilers
!                                                        S. Carroll
!    5.5  20/02/03  Modified to include NEC cpp        A.A.Dickinson
!    6.0  11/09/03  Modified to include IBM and SGI pre-processor
!                   directives (supplied by Zoe Chaplin)   P.Dando
!    6.1  18/08/04  Add DELIM=APOSTROPHE in open statement for ALTIX
!                   (supplied by Z. Chaplin).              P.Dando
!    6.2  16/05/05  Add DELIM=APOSTROPHE in open statement for HP.
!                   Simon Wilson
!    6.2  31/01/06   Changed defs from specific compilers to a general
!                    GNU/Linux option. T. Edwards
!LL
!LL  Programming standard: UM Doc Paper 3, draft version 3 (15/1/90)
!LL
!LL  Logical components covered: H3,H40
!LL
!LL  Project task: H
!LL
!LL  Documentation:  Unified Model Documentation Paper
!LL                  H- History Bricks
!LL
!*L  Interface and arguments:
!
      SUBROUTINE TEMPHIST                                               &
     &         ( UNITHIST,ICODE,CMESSAGE )
!
      IMPLICIT NONE
!
      INTEGER       UNITHIST ! In  - Temporary history file unit
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
       CHARACTER *80 FILENAME
       CHARACTER *8  NLNAME
!*L EXTERNAL subroutines called
      EXTERNAL GET_FILE
!*
!
!L
!L 1. Open, rewind and write a new record
!L
      CALL GET_FILE(UNITHIST,FILENAME,80,ICODE)
#if defined(FUJITSU) || \
    defined(LINUX) || defined(NEC) || defined(IBM) \
   || defined(SGI) || defined(ALTIX) || defined(HP)
       OPEN(UNITHIST,FILE=FILENAME,FORM='FORMATTED',                    &
     &  DELIM='APOSTROPHE',IOSTAT=ICODE)
#else
       OPEN(UNITHIST,FILE=FILENAME,FORM='FORMATTED',IOSTAT=ICODE)
#endif
!
! Check for error
!
      IF(ICODE  >  0)THEN
        CMESSAGE='TEMPHIST: Failed in OPEN of history file'
        GOTO 999
      ELSEIF(ICODE  <   0)THEN
        WRITE(6,*)'TEMPHIST: Warning message on OPEN of history file'
        WRITE(6,*)'IOSTAT= ',ICODE
      ENDIF
!
      REWIND(UNITHIST)

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

      go to 999
!
! Check for error
!
! Write error
  200 continue
      ICODE = 2
      CMESSAGE='TEMPHIST: Write ERROR on history file for namelist'//   &
     & NLNAME
!
 999  CONTINUE
!L
!L 2. Close and return
!L
      CLOSE(UNITHIST)
      RETURN
      END SUBROUTINE TEMPHIST
#endif
