

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   -------------------- SUBROUTINE UM_FORT_FLUSH ---------------------
!  
!   Purpose: A wrapper script for the flush intrinsic
!
!   -------------------------------------------------------------------

      SUBROUTINE UM_FORT_FLUSH(lunit,icode)

!     Required if using the NAG compiler




      Implicit None

!     The subroutine's arguments, whatever the compiler
      INTEGER, INTENT(IN)  :: lunit
      INTEGER, INTENT(OUT) :: icode

!     If on NAG or X1, require two 32 bit arguments to flush.
      CALL flush(lunit)
      icode = 0
      END SUBROUTINE UM_FORT_FLUSH
