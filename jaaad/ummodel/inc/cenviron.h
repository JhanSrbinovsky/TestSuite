!----------------------------------------------------------------------
! comdeck: CENVIRON
! Purpose: defines environment variables for units
!          which have to be opened by open_file
! History:
! version  date         change
! 4.5      21/09/98     New code
! Author:  M. J. Bell
!----------------------------------------------------------------------
! declaration of parameters
#include "chsunits.h"

! declarations of common blocks
      common / LEnviron /    LEnv
      common / CEnviron /    CEnv

! declarations of variables
      integer LEnv(NUnits)      ! lengths of environment variable names
      character*15 CEnv(NUnits) ! names of environment variables
!----------------------------------------------------------------------
