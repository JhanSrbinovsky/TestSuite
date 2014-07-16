#if defined(C80_1A) || defined(UTILIO) || defined(FLDOP)               \
 || defined(FLDC) || defined(VAROPSVER)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine ABORT_IO-----------------------------------------------
!LL
!LL  Purpose:  Prints out message and stops execution of program.
!LL            Called if ICODE  /=  0
!LL
!LL  Written by A. Dickinson
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author: Tracey Smith
!LL   4.1    23/05/96 Added MPP abort code - one PE calling ABORT
!LL                   will cause all others to abort     P.Burton
!LL   4.4    15/10/97 Added code to print the error messages
!LL                   to stderr as well as the unit 6
!LL                     Author: Bob Carruthers, Cray Research
!LL   4.5    08/07/98 Print only the leading non-blank
!LL                   characters in 'cmessage'
!LL                     Author: Bob Carruthers, Cray Research
!     5.3    22/11/01 Enable MPP as the only option for
!                     small executables         E.Leung
!     6.0    24/04/03 Use Fortran intrinsic len_trim instead of
!                     get_char_len. Tidy up writes. T.White
!     6.0    31/11/03 Enable ABORT_IO for FIELDCOS. P.Dando
!     6.2    10/04/05 Removed calls to ABORT.  J. Gill
!LL
!LL  Logical component number: E5
!LL
!LL  External Documentation: None
!LL
!LLEND
!*L  Arguments:--------------------------------------------------------

      SUBROUTINE ABORT_IO(STRING,CMESSAGE,ICODE,NFT)

      IMPLICIT NONE

      INTEGER                                                           &
     & ICODE                                                            &
               !IN Code returned by UM routines
     &,NFT     !IN Unit no being processed

      character*(*)                                                     &
     & STRING                                                           &
               !IN Subroutine name and position
     &,CMESSAGE!IN Message returned by UM routines

      INTEGER, PARAMETER :: STDOUT = 0  ! Standard error unit number.
      INTEGER, PARAMETER :: STDERR = 6  ! Standard error unit number.
#include "parvars.h"
!----------------------------------------------------------------------

!L Internal structure: None

      Write(STDOUT,*) 'Processor ',mype,' calling ABORT'
      Write(STDERR,*) 'Processor ',mype,' calling ABORT'

      Write(STDOUT,'('' Error detected in subroutine '',A)') string
      Write(STDERR,'('' Error detected in subroutine '',A)') string

      IF(NFT /= 0)THEN
        Write(STDOUT,'('' while doing I/O on unit'',I3)') nft
        Write(STDERR,'('' while doing I/O on unit'',I3)') nft
      ENDIF
      Write(STDOUT,'(A)')Trim(cmessage)
      Write(STDERR,'(A)')Trim(cmessage)

      Write(STDOUT,'('' ICODE='',I6)') icode
      CALL GC_ABORT(mype,nproc,CMESSAGE)
! DEPENDS ON: ereport
      CALL EREPORT('ABORT_IO', ICODE,                                   &
     & CMESSAGE)

      STOP
      END SUBROUTINE ABORT_IO

#endif
