
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Dummy routine to resolve externals
      SUBROUTINE UPDATE_LBC_ITERSL()
! History:
!   Vn    Date     Comment
!  ----   ----     -------
!  6.4  13/12/06  Original code.        Michail Diamantakis  

      IMPLICIT NONE

      INTEGER                       ::  ICODE
      CHARACTER (Len=80)            ::  CMESSAGE
      CHARACTER (Len=* ), Parameter ::  RoutineName='UPDATE_LBC_ITERSL'

      CMESSAGE = 'Routine should not be callable'
      ICODE = 1
! DEPENDS ON: ereport
      CALL EReport(RoutineName,ICODE,CMESSAGE)

      RETURN
      END SUBROUTINE UPDATE_LBC_ITERSL
