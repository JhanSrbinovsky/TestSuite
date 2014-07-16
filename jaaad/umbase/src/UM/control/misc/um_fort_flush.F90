#if defined(C70_1A) || defined(FLDIO) || defined(UTILHIST) \
 || defined(UTILIO) || defined(RECON)
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
#if defined(LINUX_NAG_COMPILER)
      USE F90_UNIX_IO,ONLY:FLUSH
#endif

      Implicit None

!     The subroutine's arguments, whatever the compiler
      INTEGER, INTENT(IN)  :: lunit
      INTEGER, INTENT(OUT) :: icode

!     If on NAG or X1, require two 32 bit arguments to flush.
#if defined(LINUX_NAG_COMPILER) || defined(_X1) || defined(XD1) \
 || defined(XT3)
      INTEGER(KIND=4) :: icode1
      INTEGER(KIND=4) :: lunit1
      lunit1 = lunit
      CALL flush(lunit1,icode1)
      icode = icode1

!     If on the T3E we require two 64 bit arguments
#elif defined(T3E)
      CALL flush(lunit,icode)

!     All others use one 64 bit argument

#elif defined(IBM)
      flush(lunit)
#else
      CALL flush(lunit)
      icode = 0
#endif
      END SUBROUTINE UM_FORT_FLUSH
#endif
