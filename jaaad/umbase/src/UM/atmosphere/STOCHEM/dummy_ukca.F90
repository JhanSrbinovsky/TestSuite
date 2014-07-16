#if defined(A25_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Dummy routine so there is something to compile for section A25_2A
      SUBROUTINE DUMMY_UKCA(ICODE,CMESSAGE)
! History:
!   Vn    Date     Comment
!  ----   ----     -------
!  6.1  10/09/04  Original code.  R Barnes

      INTEGER                       ::  ICODE
      CHARACTER (Len=*)             ::  CMESSAGE
      CHARACTER (Len=*), Parameter  ::  RoutineName='DUMMY_UKCA'

      Write(6,*)'**ERROR**: DUMMY_UKCA should not be callable'
      CMESSAGE = 'Routine DUMMY_UKCA should not be callable'
      ICODE = 1
! DEPENDS ON: ereport
      CALL EReport(RoutineName,ICODE,CMESSAGE)

      RETURN
      END SUBROUTINE DUMMY_UKCA
#endif
