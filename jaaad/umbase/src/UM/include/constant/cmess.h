!----------------------------------------------------------------------
! comdeck: CMESS
! Purpose: declares and stores values used for warning, error and
!          standard log messages
! History:
! version  date         change
! 4.5      21/09/98     New code
! Author:  M. J. Bell
!----------------------------------------------------------------------
! declarations:
! unit numbers for debugging, error, warning and standard log output
      integer OutUnitDbg   ! output unit for debugging information
      integer UnWarn
      integer UnErr
      integer UnStd
      PARAMETER (OutUnitDbg = 90, UnErr = 91, UnWarn = 92, UnStd = 93 )

! parameters for start of message
      Character*9 CWarn
      Character*7 CErr
      Character*1 CStd
      PARAMETER (CWarn  = 'Warning: ')
      PARAMETER (CErr   = 'ERROR: ')
      PARAMETER (CStd   = ' ')

! local character variable for subroutine name
      Character*20 CSub
!----------------------------------------------------------------------
