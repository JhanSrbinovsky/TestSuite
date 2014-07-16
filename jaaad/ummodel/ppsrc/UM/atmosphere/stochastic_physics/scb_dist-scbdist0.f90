
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

! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

      CMESSAGE = 'Routine should not be callable'
      ICODE = 1
! DEPENDS ON: ereport
      CALL EReport(RoutineName,ICODE,CMESSAGE)

      SCB_DIST = -RMDI

      RETURN
      END FUNCTION SCB_DIST
