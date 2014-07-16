#if !defined(OASIS3) && !defined(OASIS4)
      Subroutine oasis_inita2o(                                        &
#include "argd1.h"
#include "argsts.h"
#include "argduma.h"
#include "argptra.h"
     & icode,                                                           &
     & cmessage)


      Implicit none

!
! Description: Stub routine for oasis_inita2o. 
!              Run time controls should mean this routine is never 
!              actually called. 
!
! Author: R. Hill
! Current Code Owner : R. Hill 
!
!---------------------------------------------------------------------
#include "csubmodl.h"
#include "cntlatm.h"
#include "parparm.h"
#include "typsize.h"
#include "typd1.h"
#include "typsts.h"
#include "typduma.h"
#include "typptra.h"
!
      Integer       icode       ! OUT - Error return code
      Character*(*) cmessage    ! OUT - Error return message
!*----------------------------------------------------------------------
      CHARACTER(Len=52)   :: Message
      CHARACTER(Len=20)   :: RoutineName
      INTEGER             :: ErrorStat        ! Return code:
                                              !   0 = Normal exit
                                              ! +ve = Fatal Error
                                              ! -ve = Warning

      ErrorStat   = 1
      RoutineName = 'oasis_inita2o_stub'
      Message     = 'OASIS Routines unavailable - see output.'

      WRITE (6,*) '**ERROR**: oasis_inita2o unavailable.'
      WRITE (6,*) 'Check OASIS cpp key is set'

! DEPENDS ON: ereport
      CALL ereport(RoutineName, ErrorStat, Message)

      End Subroutine oasis_inita2o
#endif
