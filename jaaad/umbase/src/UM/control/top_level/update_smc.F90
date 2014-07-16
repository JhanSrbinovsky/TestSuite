#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine UPDATE_SMC
!
!   Model
!  version  Date
!   6.2   23/02/06  New Deck.  Clive Jones
!
!  Programming standard : UM documentation paper no3,
!                         version no.1, dated 15/01/90
!
!  Purpose: To update the partitioning between unfrozen and
!           frozen soil moisture after the total soil moisture in a
!           layer has been updated.

      Subroutine UPDATE_SMC (                                           &
#include "argd1.h"
#include "argptra.h"
     &                   DZ_SOIL)

      IMPLICIT NONE

!*L Arguments

!L

!*L Comdecks
#include "parparm.h"
#include "typsize.h"
#include "typd1.h"
#include "typptra.h"
#include "cntlatm.h"
#include "nstypes.h"
#include "cancila.h"

      REAL               DZ_SOIL(SM_LEVELS) !IN soil level thicknesses

      WRITE(6,*)                                                        &
     & 'Partitioning soil moisture in unfrozen and frozen fractions'
! DEPENDS ON: freeze_soil
      CALL FREEZE_SOIL(LAND_FIELD,SM_LEVELS,                            &
     &                 D1(JCLAPP_HORN),                                 &
     &                 DZ_SOIL,                                         &
     &                 D1(JSAT_SOILW_SUCTION),                          &
     &                 D1(JSMCL(1)),                                    &
     &                 D1(J_DEEP_SOIL_TEMP(1)),                         &
     &                 D1(JVOL_SMC_SAT),                                &
     &                 D1(JSTHU(1)),                                    &
     &                 D1(JSTHF(1)))


      Return
      END SUBROUTINE UPDATE_SMC
#endif
