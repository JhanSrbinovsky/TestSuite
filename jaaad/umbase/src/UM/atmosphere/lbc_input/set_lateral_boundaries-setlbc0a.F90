#if defined(A31_0A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Dummy routine to resolve externals
      SUBROUTINE SET_LATERAL_BOUNDARIES()
! History:
!   Vn    Date     Comment
!  ----   ----     -------
!  6.2  23/05/06  Original code.  P. Selwood

      IMPLICIT NONE

      INTEGER                       ::  ICODE
      CHARACTER (Len=80)            ::  CMESSAGE
      CHARACTER (Len=* ), Parameter ::                                  &
     &                   RoutineName='SET_LATERAL_BOUNDARIES'

      CMESSAGE = 'Routine should not be callable'
      ICODE = 1
! DEPENDS ON: ereport
      CALL EReport(RoutineName,ICODE,CMESSAGE)

      RETURN
      END SUBROUTINE SET_LATERAL_BOUNDARIES


#endif
