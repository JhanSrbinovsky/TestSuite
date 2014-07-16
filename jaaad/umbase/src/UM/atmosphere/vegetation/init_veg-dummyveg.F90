#if defined(A19_0A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Dummy version of routine INIT_VEG
      SUBROUTINE INIT_VEG(A_STEP,                                       &
#include "argd1.h"
#include "argduma.h"
#include "argptra.h"
#include "arglndm.h"
     &              ICODE,CMESSAGE)

! History:
! Version   Date     Comment
! -------   ----     -------
!   6.0    11/09/03  Replaced ABORT call with call to EREPORT. P.Dando
!
      INTEGER                                                           &
     & A_STEP             ! IN Current timestep in atmosphere model

      INTEGER                          ::  ICODE
      CHARACTER (Len=*)                ::  CMESSAGE
      CHARACTER (Len=*),  Parameter    ::  RoutineName='INIT_VEG'

      WRITE (6,*) '**ERROR**: INIT_VEG has been called but is'
      WRITE (6,*) 'unavailable.  Either set L_VEG_FRACS to .FALSE. or'
      WRITE (6,*) 'select section A19_1A or A19_2A.'

      CMESSAGE = 'Routine unavailable - see output for details'
      ICODE = 1
! DEPENDS ON: ereport
      CALL EReport(RoutineName,ICODE,CMESSAGE)

      RETURN
      END SUBROUTINE INIT_VEG
#endif
