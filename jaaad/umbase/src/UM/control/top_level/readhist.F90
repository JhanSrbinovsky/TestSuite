#if defined(CONTROL) || defined(COMB) || defined(PICK) || defined(HPRT)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: READHIST
!LL
!LL  Purpose: To initialise History common block from most recent
!LL           record in history file input
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
!    5.4  09/04/02   Work-around for some Linux compilers to make sure
!                    the read position is in the right spot for the
!                    required namelist.                     S. Carroll
!    5.5  20/02/03   Modified to include NEC cpp        A.A.Dickinson
!    6.1  18/08/04   Fix namelist read on Linux platforms - also
!                    works on NEC.  P.Dando
!    6.2  31/01/06   Changed defs from specific compilers to a general
!                    GNU/Linux option. T. Edwards
!LL
!LL  Programming standard: UM Doc Paper 3, draft version 3 (15/1/90)
!LL
!LL  Logical components covered: H2,H20
!LL
!LL  Project task: H
!LL
!LL  Documentation:  Unified Model Documentation Paper
!LL                  H- History Bricks
!LL
!
!*L  Interface and arguments:
!
      SUBROUTINE READHIST                                               &
     &         ( UNITHIST,ICODE,CMESSAGE )
!
      IMPLICIT NONE
!
      INTEGER       UNITHIST ! In  - History file unit
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
#if defined(HP)
      CHARACTER*80 TMP_LINE
#endif
      CHARACTER*80 FILENAME
      CHARACTER*8  NLNAME
!*L EXTERNAL subroutines called
      EXTERNAL INITCHST,GET_FILE
!*
!L
!L 1. Set common block area to zero or blank
!L
! DEPENDS ON: initchst
      CALL INITCHST
      ICODE = 0
!L
!L 2. Open history file and rewind
!L
      CALL GET_FILE(UNITHIST,FILENAME,80,ICODE)
        OPEN(UNITHIST,FILE=FILENAME,FORM='FORMATTED',IOSTAT=ICODE)
!
! Check for error
!
      IF(ICODE  >  0)THEN
        CMESSAGE='READHIST: Failed in OPEN of permanent history file'
        GOTO 999
      ELSEIF(ICODE  <   0)THEN
      WRITE(6,*)'READHIST: Warning message on OPEN of permanent history &
     & file'
        WRITE(6,*)'IOSTAT= ',ICODE
      ENDIF
!
      REWIND(UNITHIST)
!L
!L 3. Read most recent records
!L

      NLNAME = 'NLIHISTO'
#if defined(LINUX) \
 || defined(NEC)
      REWIND(UNITHIST)
#elif defined(HP)
!  Code to position read on the correct namelist
!  Rewind file to the start
      REWIND(UNITHIST)
!  Read line into local variable TMP_LINE
      READ(UNITHIST,*,END=100,ERR=200) TMP_LINE
!  Loop through file looking for "&[NLNAME]"
!  then stop and go back one line
      DO WHILE( TMP_LINE  /=  "&"//NLNAME)
        READ( UNITHIST,*,END=100,ERR=200) TMP_LINE
      END DO
      BACKSPACE UNITHIST
#endif
      READ(UNITHIST,NLIHISTO,END=100,ERR=200)

!  ----------

      NLNAME = 'NLCHISTO'
#if defined(LINUX) \
 || defined(NEC)
      REWIND(UNITHIST)
#elif defined(HP)
!  Code to position read on the correct namelist
!  Rewind file to the start
      REWIND(UNITHIST)
!  Read line into local variable TMP_LINE
      READ(UNITHIST,*,END=100,ERR=200) TMP_LINE
!  Loop through file looking for "&[NLNAME]"
!  then stop and go back one line
      DO WHILE( TMP_LINE  /=  "&"//NLNAME)
        READ( UNITHIST,*,END=100,ERR=200) TMP_LINE
      END DO
      BACKSPACE UNITHIST
#endif
      READ(UNITHIST,NLCHISTO,END=100,ERR=200)

!  ----------

      NLNAME = 'NLIHISTG'
#if defined(LINUX) \
 || defined(NEC)
      REWIND(UNITHIST)
#elif defined(HP)
!  Code to position read on the correct namelist
!  Rewind file to the start
      REWIND(UNITHIST)
!  Read line into local variable TMP_LINE
      READ(UNITHIST,*,END=100,ERR=200) TMP_LINE
!  Loop through file looking for "&[NLNAME]"
!  then stop and go back one line
      DO WHILE( TMP_LINE  /=  "&"//NLNAME)
        READ( UNITHIST,*,END=100,ERR=200) TMP_LINE
      END DO
      BACKSPACE UNITHIST
#endif
      READ(UNITHIST,NLIHISTG,END=100,ERR=200)

!  ----------

      NLNAME = 'NLCHISTG'
#if defined(LINUX) \
 || defined(NEC)
      REWIND(UNITHIST)
#elif defined(HP)
!  Code to position read on the correct namelist
!  Rewind file to the start
      REWIND(UNITHIST)
!  Read line into local variable TMP_LINE
      READ(UNITHIST,*,END=100,ERR=200) TMP_LINE
!  Loop through file looking for "&[NLNAME]"
!  then stop and go back one line
      DO WHILE( TMP_LINE  /=  "&"//NLNAME)
        READ( UNITHIST,*,END=100,ERR=200) TMP_LINE
      END DO
      BACKSPACE UNITHIST
#endif
      READ(UNITHIST,NLCHISTG,END=100,ERR=200)

!  ----------

      NLNAME = 'NLCFILES'
#if defined(LINUX) \
 || defined(NEC)
      REWIND(UNITHIST)
#elif defined(HP)
!  Code to position read on the correct namelist
!  Rewind file to the start
      REWIND(UNITHIST)
!  Read line into local variable TMP_LINE
      READ(UNITHIST,*,END=100,ERR=200) TMP_LINE
!  Loop through file looking for "&[NLNAME]"
!  then stop and go back one line
      DO WHILE( TMP_LINE  /=  "&NLCFILES")
        READ( UNITHIST,*,END=100,ERR=200) TMP_LINE
      END DO
      BACKSPACE UNITHIST
#endif
      READ(UNITHIST,NLCFILES,END=100,ERR=200)

      go to 999
!
! Check for error
!
! End-of-file
  100 continue
      ICODE = 1
      CMESSAGE='READHIST: End of file in READ from history file for name&
     &list '//NLNAME
      go to 999
! Read error
  200 continue
      ICODE = 2
      CMESSAGE='READHIST: Read ERROR on history file for namelist '//   &
     & NLNAME
!
 999  CONTINUE
!L
!L 4. Close and return
!L
      CLOSE(UNITHIST)
      RETURN
      END SUBROUTINE READHIST
#endif
