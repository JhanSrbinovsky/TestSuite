#if !defined(OASIS3) && !defined(OASIS4)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
      Subroutine oasis_updatecpl(                                       &                       
#include "argd1.h" 
#include "arg_atm_fields.h"
     &  cmessage)

      Implicit none
!
!
! Description:
! Stub routine for oasis_updatecpl to allow compilation in 
! stand-alone atmosphere jobs. Run time logic should
! mean this version never actually gets called. 
!
! Author: R. Hill
! Current Code Owner : Richard Hill
!
!=====================================================================
#include "parvars.h"   
#include "decomptp.h" 
#include "decompdb.h" 
#include "atm_lsm.h"  

#include "cmaxsize.h" 
#include "typsize.h"  
#include "typd1.h"    
#include "typptra.h"  
#include "typcona.h"  
#include "typlndm.h"  
#include "ctracera.h" 
#include "typ_atm_fields.h"

#include "caoptr.h"   


      Character*(*) :: cmessage ! OUT - Error return message

      CHARACTER(Len=52)   :: Message
      CHARACTER(Len=20)   :: RoutineName
      INTEGER             :: ErrorStat        ! Return code:
                                              !   0 = Normal exit
                                              ! +ve = Fatal Error
                                              ! -ve = Warning

      ErrorStat   = 1
      RoutineName = 'oasis_updatecpl_stub'
      Message     = 'OASIS Routines unavailable - see output.'

      WRITE (6,*) '**ERROR**: oasis_updatecpl is unavailable.'
      WRITE (6,*) 'Check OASIS3 or OASIS4 cpp key is set'

! DEPENDS ON: ereport
      CALL ereport(RoutineName, ErrorStat, Message)

      End Subroutine oasis_updatecpl
#endif
