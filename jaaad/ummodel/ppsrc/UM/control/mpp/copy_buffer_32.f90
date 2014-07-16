

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
! Start C_KINDS

! Description:
!   Contains parameters defining kinds for 32 and 64 integers
!   and reals.
!
! Current Code Owner: Paul Selwood
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   6.1  10/03/04   Original code. Paul Selwood

! Parameters for 32 and 64 bit kinds

! Precision and range for 64 bit real
      Integer, Parameter :: prec64  = 15
      Integer, Parameter :: range64 = 307

! Precision and range for 32 bit real
      Integer, Parameter :: prec32  = 6
      Integer, Parameter :: range32 = 37

! Range for integers
      Integer, Parameter :: irange64=15
      Integer, Parameter :: irange32=9

! Kind for 64 bit real
      Integer, Parameter :: real64  = selected_real_kind(prec64,range64)
! Kind for 32 bit real
      Integer, Parameter :: real32  = selected_real_kind(prec32,range32)
! Kind for 64 bit integer
      Integer, Parameter :: integer64 = selected_int_kind(irange64)
! Kind for 32 bit integer
      Integer, Parameter :: integer32 = selected_int_kind(irange32)

! Kinds for 64 and 32 bit logicals. Note that there is no
! "selected_logical_kind", but using the equivalent integer kind is a
! workaround that works on every platform we've tested.
      Integer, Parameter :: logical64 = integer64
      Integer, Parameter :: logical32 = integer32

! End C_KINDS

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
