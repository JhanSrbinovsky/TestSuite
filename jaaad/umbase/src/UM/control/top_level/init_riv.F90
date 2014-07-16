#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      SUBROUTINE INIT_RIV(                                              &
#include "argd1.h"
#include "argduma.h"
#include "argptra.h"
     &              ICODE,CMESSAGE)

! Purpose: To initialise variables for river routing. This will be
! extended later.
! Method: Sets up variables in the dump

      IMPLICIT NONE
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 5.5      28/03/03  Original code. C. Bunton
! 6.2      21/03/06  Included nstypes.h J Ridley.
!
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!

      INTEGER ICODE                 ! OUT Internal return code
      CHARACTER*80 CMESSAGE         ! OUT Internal error message

#include "nstypes.h"
#include "cmaxsize.h"
#include "csubmodl.h"
#include "parvars.h"
#include "typsize.h"
#include "typd1.h"
#include "typduma.h"
#include "typptra.h"
#include "chsunits.h"
#include "ccontrol.h"
#include "cruntimc.h"
#include "ctime.h"

! A_INTHD(16) holds the number of atmosphere timesteps since the last
! call to River routing


      ICODE=0
      CMESSAGE=""
      A_INTHD(16)=0



      RETURN
      END SUBROUTINE INIT_RIV


#endif
