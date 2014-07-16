#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
 || defined(UTILIO) || defined(FLDIO) || defined(MAKEBC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ A copy from an input 32-bit buffer to an output 32-bit buffer

      Subroutine copy_buffer_32( in, out, count, off_out, off_in )

      Implicit None
!
! Description:
!   Copies input buffer (32 bit data) to output buffer (32 bit data)
!   possibly with offsets in either input or output.
!
! Method:
!   Very simple copy of input to output with offsets.
!   Only required as a seperate deck to ensure 32 bit sizes are
!   respected (on NEC).
!
! Current Code Owner: Paul Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   6.1   13/07/04   Original code. Paul Selwood.
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables (#include statements etc):
#include "c_kinds.h"

! Arguments
      Integer (kind=integer64), Intent(In) ::                           &
     &         count,                                                   &
                                    ! size of data to copy
     &         off_in,                                                  &
                                    ! input offset
     &         off_out              ! output offset

      Integer (kind=integer32), Intent(InOut) ::                        &
     &         in(*),                                                   &
                                    ! input data
     &         out(*)               ! output data

      Integer :: i                  ! Looper

! End of header

      Do i = 1, count - off_in
        out(off_out + i) = in(i+off_in)
      End Do

      Return
      End Subroutine copy_buffer_32
#endif
