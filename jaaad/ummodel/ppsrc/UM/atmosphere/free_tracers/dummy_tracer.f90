
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Dummy routine so there is something to compile for section a33_1a
      SUBROUTINE DUMMY_TRACER(ICODE,CMESSAGE)
! History:
!   Vn    Date     Comment
!  ----   ----     -------
!  6.1  10/09/04  Original code.  R Barnes

      INTEGER                       ::  ICODE
      CHARACTER (Len=*)             ::  CMESSAGE
      CHARACTER (Len=*), Parameter  ::  RoutineName='DUMMY_TRACER'

      Write(6,*)'**ERROR**: DUMMY_TRACER should not be callable'
      CMESSAGE = 'Routine DUMMY_TRACER should not be callable'
      ICODE = 1
! DEPENDS ON: ereport
      CALL EReport(RoutineName,ICODE,CMESSAGE)

      RETURN
      END SUBROUTINE DUMMY_TRACER
