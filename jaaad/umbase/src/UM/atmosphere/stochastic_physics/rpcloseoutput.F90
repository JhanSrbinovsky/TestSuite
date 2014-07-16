#if defined(CONTROL) && defined(ATMOS) && defined(A35_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine RPCloseOutput
! ------------------------
! Closes the RPSEED output stream if it is open

 SUBROUTINE RPCloseOutput()
 IMPLICIT NONE
#include "cntlatm.h"

! Local Variables
 INTEGER :: icode ! Error return code

 IF (L_RPSEED_WRITE) THEN
   CLOSE(unit=149,iostat=icode)
   L_RPSEED_WRITE=.false.
   IF(icode > 0) THEN
     WRITE(6,*) "RPSEED OUT : Failed to close output file"
   ENDIF
 ENDIF
 END SUBROUTINE RPCloseOutput
#endif
