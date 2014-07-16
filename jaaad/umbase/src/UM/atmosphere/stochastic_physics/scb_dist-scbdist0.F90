#if defined(A35_0A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Dummy routine to resolve externals
      REAL FUNCTION SCB_DIST()
! History:
!   Vn    Date     Comment
!  ----   ----     -------
!  6.2  23/05/06  Original code.  P. Selwood
!  6.4  15/12/06  Change Subroutine to Function.  P. Selwood

      IMPLICIT NONE

      INTEGER                       ::  ICODE
      CHARACTER (Len=80)            ::  CMESSAGE
      CHARACTER (Len=* ), Parameter ::  RoutineName='SCB_DIST'

#include "c_mdi.h"

      CMESSAGE = 'Routine should not be callable'
      ICODE = 1
! DEPENDS ON: ereport
      CALL EReport(RoutineName,ICODE,CMESSAGE)

      SCB_DIST = -RMDI

      RETURN
      END FUNCTION SCB_DIST
#endif
