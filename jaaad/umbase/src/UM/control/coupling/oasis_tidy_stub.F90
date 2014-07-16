#if !defined(OASIS3) && !defined(OASIS4)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
      Subroutine oasis_tidy(                                        &
#include "argd1.h"
#include "argsts.h"
#include "argduma.h"
#include "argptra.h"
     & icode,cmessage)

      Implicit none

!
! Description: Stub routine for oasis_tidy.
!              Run-time logic should mean this is never called. 
!
! Author: R. Hill
! Current Code Owner : R. Hill 
!
!--------------------------------------------------------------------
#include "csubmodl.h"
#include "cntlatm.h"
#include "parvars.h"
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
      RoutineName = 'oasis_tidy_stub'
      Message     = 'OASIS3 Routines unavailable - see output.'

      WRITE (6,*) '**ERROR**: oasis_tidy unavailable.'
      WRITE (6,*) 'Check OASIS3 cpp key is set'

! DEPENDS ON: ereport
      CALL ereport(RoutineName, ErrorStat, Message)


      End Subroutine oasis_tidy
#endif
