#if defined(CONTROL) && defined(ATMOS) && defined(A35_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine RPOpenOutput
! -----------------------
! Opens the file stream constructing a filename from
! the the RPSEED variables that has been set in a
! namelist. If this is blank or just spaces, then
! construct a name from the $DATAW and $RUNID
! environment variables.
!
! Subroutine RPCloseOutput
! ------------------------
! Closes the output stream if it is open
!
! Subroutine RPWriteEntry(value)
! ------------------------------
! This routine is passed the random seed in its
! argument and it then writes this to the output file if
! it is open.

SUBROUTINE RPwriteEntry(seed,seedsize)
IMPLICIT NONE

#include "cntlatm.h"
! Arguments
 INTEGER, INTENT(IN) :: seedsize       ! The size of the random seed
 INTEGER, INTENT(IN) :: seed(seedsize) ! The random seed

! Local Variables
 INTEGER :: i                          ! loop variable

! If the unit is open the write the value to the file stream
 IF (L_RPSEED_WRITE) then
    DO i=1,seedsize
       WRITE(149,*) seed(i)
    ENDDO
 ENDIF
 END SUBROUTINE RPwriteEntry
#endif
