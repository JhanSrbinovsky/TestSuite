

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  SUBROUTINE EXPAND21
!
!  Purpose: Unpacks IEEE 32-bit data into IEEE 64-bit data.
!
!  Author:  Bob Carruthers, Cray Research
!
!  Tested under compiler:   f90
!  Tested under OS version: UNICOS/mk 1.4.1
!
!  Model            Modification history
! version  date
!
!   4.3  18/03/97   Original version
!   5.2  31/07/00   Change to only 4 arguments as last one is never
!                   used and confuses -Ra compiler option  P.Selwood
!   5.5  20/02/03   Mod to enable routine to work on non-CRAY
!                   platforms.                        P.Dando
!   6.1  22/03/04   Replace kinds definitions by c_kinds.
!                   P.Selwood
!
subroutine expand21(n, in, out)
!--expands input array 'in' from 32-bit into 'out' in
!  64-bit
!
!  n       the number of floating point words to convert
!  in      the input array of 32-bit numbers
!  out     the output array of 64-bit numbers
!
!--
implicit none
!
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

! Argument variables
integer (kind=integer64), intent(in) :: n
real (kind=real32), intent(in)       :: in(1:n)
real (kind=real64), intent(out)      :: out(1:n)

! Local integer variables
integer (kind=integer64) :: i

do i = 1,n
  out(i) = in(i)
end do

return
END SUBROUTINE expand21
