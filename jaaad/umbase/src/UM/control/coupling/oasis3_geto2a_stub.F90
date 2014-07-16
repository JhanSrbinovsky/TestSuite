#if !defined(OASIS3)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE oasis3_geto2a(                                         &
#include "argd1.h"
#include "arg_atm_fields.h"
  &  GET_STEP,cmessage)

  IMPLICIT NONE

  !
  ! Description:
  ! Dummy version of oasis3_geto2a - should never be called.
  !
  ! Author: R. Hill
  !
  ! Current Code Owner : R. Hill
  !
  !
  !=====================================================================

#include "parvars.h"
#include "typsize.h"
#include "typd1.h"
#include "ctracera.h"
#include "typ_atm_fields.h"

  !     Subroutine arguments
  LOGICAL :: GET_STEP       ! Whether to transfer incoming
  ! data to UM atmos D1 or just ignore it.

  CHARACTER*(*) :: cmessage ! OUT - Error return message

  CHARACTER(Len=52)   :: Message
  CHARACTER(Len=20)   :: RoutineName
  INTEGER             :: ErrorStat        ! Return code:
  !   0 = Normal exit
  ! +ve = Fatal Error
  ! -ve = Warning

  ErrorStat   = 1
  RoutineName = 'oasis3_geto2a_stub'
  Message     = 'OASIS3 Routines unavailable - see output.'

  WRITE (6,*) '**ERROR**: oasis3_geto2a called but is unavailable.'
  WRITE (6,*) 'Check OASIS3 cpp key is set'


  ! DEPENDS ON: ereport
  CALL ereport(RoutineName, ErrorStat, Message)



END SUBROUTINE oasis3_geto2a
#endif
