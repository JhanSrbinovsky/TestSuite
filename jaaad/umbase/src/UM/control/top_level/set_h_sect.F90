#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Read run-time control information from namelists for atmos model
!
! Subroutine Interface:

      SUBROUTINE SET_H_SECT(H_SECT,MAXSECTS)
!  Temporary routine introduced at vn3.5 to set H_SECT in COMDECK
!  CCONTROL to values in ATMOS_SR in COMDECK MODEL.
!  Necessary bc. CCONTROL & MODEL contain some variable names the same
!  - this needs to be tidied at a future version.
      IMPLICIT NONE
      INTEGER MAXSECTS
      CHARACTER*3 H_SECT(0:MAXSECTS)
#include "csubmodl.h"
#include "version.h"
#include "model.h"
!  Local variables
      INTEGER I ! loop variable

      DO  I = 0,NSECTP
        H_SECT(I)(1:1) = '0'
        H_SECT(I)(2:3) = ATMOS_SR(I)
      END DO

      RETURN
      END SUBROUTINE SET_H_SECT
#endif
