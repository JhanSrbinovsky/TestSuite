#if defined(CONTROL)
#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL--------------- SUBROUTINE IN_ACCTL --------------------------------
!LL
!LL Purpose : Calls the initialisation program for the data assimilation
!LL assimilation section (P3)
!LL
!LL  Control routine for CRAY YMP
!LL
!LL M.Bell      <- programmer of some or all of previous code or changes
!LL
!LL Programming standard; U M Documentation Paper No. 3 version 1
!LL dated 15/01/90
!LL
!LL System components covered P0
!LL
!LL Documentation : UM documentation paper  no P0
!LL
!LL END
!----------------------------------------------------------------------

      SUBROUTINE IN_ACCTL(                                              &
#include "argduma.h"
#include "argptra.h"
#include "argppx.h"
     &                  ICODE,CMESSAGE)

      IMPLICIT NONE

!
      INTEGER       ICODE             ! OUT: Error return code
      CHARACTER*256 CMESSAGE          ! OUT: Error return message
!
#include "parparm.h"

#include "typsize.h"
#include "typduma.h"
#include "typptra.h"
#include "cocnindx.h"

#include "cmaxsize.h"
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "ctime.h"
#include "chsunits.h"
#include "ccontrol.h"

#if defined(A18_2A)
#include "csizeobs.h"
#include "chistory.h"
#endif

#if defined(A18_2A)

      IF(L_AC) THEN

! DEPENDS ON: ac_init
      CALL AC_INIT (MODEL_LEVELS, WET_LEVELS, BL_LEVELS, TR_LEVELS,     &
     &              ROWS, N_ROWS, ROW_LENGTH,                           &
     &              A_MAX_OBS_SIZE, A_MAX_NO_OBS,                       &
     &              SECS_PER_STEPim(atmos_im),                          &
     &              MODEL_BASIS_TIME(1), MODEL_BASIS_TIME(2),           &
     &              MODEL_BASIS_TIME(3), MODEL_BASIS_TIME(4),           &
     &              MODEL_BASIS_TIME(5),                                &
     &              A_REALHD(1), A_REALHD(2), A_REALHD(3),              &
     &              A_REALHD(4), A_REALHD(5), A_REALHD(6),              &
!  Interim substitution to allow compilation at 5.0. If this needs to
!  be replaced by eta_theta_levels/eta_rho_levels, these can be
!  accessed by adding </argcona.h> to argument list [will also need
!  </parvars.h>,<typcona.h>:
!    &              A_LEVDEPC(JAK), A_LEVDEPC(JBK),
     &              A_LEVDEPC(1), A_LEVDEPC(1),                         &
#include "argppx.h"
     &              ICODE,CMESSAGE)

      END IF

#endif
!
      RETURN
      END SUBROUTINE IN_ACCTL

#endif
#endif
