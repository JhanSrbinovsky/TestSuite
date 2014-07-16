#if defined(C80_1A) || defined(UTILIO)                                 \
 || defined(RECON) || defined(FLDIO) || defined(VAROPSVER)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  SUBROUTINE PACK21
!
!  Purpose: Packs IEEE 64-bit data into IEEE 32-bit data.
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
!   5.5  17/04/03   Remove reference to obsolete section
!                   C91_2A. T.White
!   5.5  20/02/03   Mod to enable routine to work on non-CRAY
!                   platforms.                        P.Dando
!   6.1  22/03/04   Replace kinds definitions by c_kinds.
!                   P.Selwood
!
SUBROUTINE pack21(n, in, out)
!--compresses input array 'in' from 64-bit into 'out' in
!  32-bit
!
!  n       the number of floating point words to convert
!  in      the input array of 64-bit numbers
!  out     the output array of 32-bit numbers
!
!--
implicit none

#include "c_kinds.h"

! Argument variables
integer (kind=integer64), intent(in) :: n
real (kind=real64), intent(in)       :: in(1:n)
real (kind=real32), intent(out)      :: out(1:n)

! Local real parameter
real, parameter :: tiny32=tiny(out(1))
! Local variables 
integer (kind=integer64) :: i

do i = 1,n
  if (ABS(in(i))  <   tiny32) then
! Prevent 'denormalized' numbers
    out(i) = 0.0
  else
    out(i) = in(i)
  endif
end do

return
END SUBROUTINE pack21
#endif
