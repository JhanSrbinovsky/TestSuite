#if defined(C80_1A) || defined(UTILIO)                                 \
 || defined(FLDIO) || defined(RECON) || defined(VAROPSVER)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE IOERROR----------------------------------------
!LL
!LL  Purpose: Prints out a message after using buffer in/out when
!LL           either a return code < 0.0 is encountered
!LL           by UNIT function or value returned by LENGTH
!LL           differs from length of I/O request.
!LL
!LL  Written by A. Dickinson
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL   4.1    12/06/96 Break up write statement. D. Robinson.
!LL   4.4    15/10/97 Added code to print the error message to
!LL                   stderr, and call abort in case all the
!LL                   PE's have not detected the error condition.
!LL                     Author: Bob Carruthers, Cray Research
!LL   4.5    08/07/98 Print only the leading non-blank
!LL                   characters in 'string'
!LL                     Author: Bob Carruthers, Cray Research
!LL   5.0    05/03/99 Remove DEF,C98_1A. D. Robinson.
!     6.0    24/04/03  Use Fortran intrinsic len_trim instead of
!                      get_char_len. T.White
!     6.2    12/08/05 Fix continuations for free format. P.Selwood
!     6.2    12/05/05 Remove call to abort.     P.Dando
!LL
!LL  Programming standard: Unified Model Documentation Paper No 3
!LL                        Version No 1 15/1/90
!LL
!LL  Logical component: E4
!LL
!LL  System task: F3
!LL
!LL  Documentation: CFT77 reference manual SR-0018 C  Page 9-3
!LL------------------------------------------------------------
!*L Arguments:-------------------------------------------------
      SUBROUTINE IOERROR(STRING,ERROR,LEN_IO,LEN_IO_REQ)

      IMPLICIT NONE

      INTEGER  :: LEN_IO      ! Number of 64-bit words transferred as
                              !  registered  by LENGTH function
      INTEGER  :: LEN_IO_REQ  ! Number of 64-bit words requested for
                              ! transfer via BUFFER IN/OUT

      character*(*) string

      REAL                                                              &
     & ERROR   ! Error code returned by UNIT function

! -------------------------------------------------------------
! Workspace usage:---------------------------------------------
! None
! -------------------------------------------------------------
!*L External subroutines called:-------------------------------
! None
!*-------------------------------------------------------------

!L Internal structure: none

      WRITE(6,'('' **FATAL ERROR WHEN READING/WRITING MODEL DUMP**'')')
      WRITE(6,'('' '',A)') TRIM(STRING)
      WRITE(6,'('' Error code = '',F6.2)') ERROR
      WRITE(6,'('' Length requested            = '',I9)') LEN_IO_REQ
      WRITE(6,'('' Length actually transferred = '',I9)') LEN_IO
      WRITE(6,'(''  Fatal error codes are as follows:'')')
      WRITE(6,'('' -1.0 Mismatch between actual and requested data'',   &
     &          '' length'')')
      WRITE(6,'(''  0.0 End-of-file was read'')')
      WRITE(6,'(''  1.0 Error occurred during read'')')
      WRITE(6,'(''  2.0 Other disk malfunction'')')
      WRITE(6,'('' 3.0 File does not exist'')')
      WRITE(6,'('' ***********************************************'')')
!
      write(0,'(//)')
      WRITE(0,'('' **FATAL ERROR WHEN READING/WRITING MODEL DUMP**'')')
      WRITE(0,'('' '',A)') TRIM(STRING)
      WRITE(0,'('' Error code = '',F6.2)') ERROR
      WRITE(0,'('' Length requested            = '',I9)') LEN_IO_REQ
      WRITE(0,'('' Length actually transferred = '',I9)') LEN_IO
      WRITE(0,'(''  Fatal error codes are as follows:'')')
      WRITE(0,'('' -1.0 Mismatch between actual and requested data'',   &
     &          '' length'')')
      WRITE(0,'(''  0.0 End-of-file was read'')')
      WRITE(0,'(''  1.0 Error occurred during read'')')
      WRITE(0,'(''  2.0 Other disk malfunction'')')
      WRITE(0,'('' 3.0 File does not exist'')')
      WRITE(0,'('' ***********************************************'')')

      RETURN
      END SUBROUTINE IOERROR
#endif
