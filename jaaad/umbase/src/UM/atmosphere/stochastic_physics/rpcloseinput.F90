#if defined(CONTROL) && defined(ATMOS) && defined(A35_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine RPCloseInput
! ------------------------
! Closes the input RPSEED stream if it is open

 SUBROUTINE RPCloseInput()
 IMPLICIT NONE
#include "cntlatm.h"

! Local Variables
 INTEGER :: icode ! Error return code
 icode=0
 
 IF (L_RPSEED_READ) THEN
       CLOSE(unit=149,iostat=icode)
   IF(icode > 0) THEN
     WRITE(6,*) "RPSEED IN : Failed to close input file"
   ELSE
     WRITE(6,*) "RPSEED IN : Input file closed"
   ENDIF
 ENDIF
 END SUBROUTINE RPCloseInput

#endif
