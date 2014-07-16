#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculate diagnostic quantities from the initial atmosphere dump
!
! Subroutine Interface:
      SUBROUTINE InitDiag(                                              &
#include "argd1.h"
#include "argduma.h"
#include "argsts.h"
#include "arg_atm_fields.h"
#include "argcona.h"
#include "arglndm.h"
#include "argppx.h"
     & Dummy)

      IMPLICIT NONE
!
! Description:
!   InitDiag processes diagnostic requests from the initial atmosphere
!   dump, including both prognostic variables resident in D1 - UM
!   STASH section 0 - and fields derived from physics and dynamics
!   variables, as calculated in UM STASH sections 15 and 16.
!
! Method:
!   1. Call STASH to process diagnostics requests for section 0.
!   2. Process dynamics-derived diagnostics (section 15):
!    2a. Check validity of diagnostic requests (avoids re-checking
!        within every later call of St_diag1).
!    2b. Call St_diag1 as interface to: Dyn_diag calculations of
!        dynamics derived variables and STASH (section 15).
!   3. Process physics-derived diagnostics (section 16):
!    3a. Check validity of diagnostic requests (avoids re-checking
!        within every later call of St_diag2).
!    3b. Call St_diag2 as interface to: Phy_diag calculations of
!        physics derived variables and STASH (section 16).
!
! Current Code Owner: R Rawlins
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 5.0    22/06/99  Original code based on a substantial revision to
!                  INITDIA1 deck, for the 'C-P C dynamics upgrade'
!                  project. R Rawlins.
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):

#include "parvars.h"
#include "cmaxsize.h"
#include "csubmodl.h"

#include "typsize.h"
#include "typd1.h"
#include "typduma.h"
#include "typsts.h"
#include "ctracera.h"
#include "typ_atm_fields.h"
#include "typcona.h"
#include "typlndm.h"
#include "ppxlook.h"

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER                                                           &
     & Dummy            ! Not used, needed to end arg list
!   Array  arguments with intent(in):

!   Scalar arguments with intent(InOut):

!   Array  arguments with intent(InOut):

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):

! Local parameters:
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='InitDiag')

! Local scalars:
!   ErrorStatus
      INTEGER      ErrorStatus          ! Error flag (0 = OK)
      CHARACTER*256                                                     &
     & CMessage         ! Error message if return code >0

      INTEGER                                                           &
     & im_index      !  Internal Model Index for stash arrays
! Local dynamic arrays:

! Function & Subroutine calls:
      External STASH,St_diag1,St_diag2,Ereport

!- End of header

!   0. Initialisation

      ErrorStatus = 0
      Cmessage=''
      im_index = internal_model_index(atmos_im)

!----------------------------------------------------------------------
!   1. Call STASH to process diagnostics requests for section 0.

! DEPENDS ON: stash
         CALL STASH(a_sm,a_im,0,D1,                                     &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &    ErrorStatus,Cmessage)

!----------------------------------------------------------------------
!   2. Process dynamics-derived diagnostics (section 15):
!    2a. Check validity of diagnostic requests (avoids re-checking
!        within every later call of St_diag1).
!    2b. Call St_diag1 as interface to: Dyn_diag calculations of
!        dynamics derived variables and STASH (section 15).

      IF(      SF(0,15)                                                 &
                                ! Diagnostics required for this section
     &   .AND. ErrorStatus == 0) THEN

! DEPENDS ON: st_diag1
      CALL St_diag1( STASH_MAXLEN(15,im_index),                         &
#include "argd1.h"
#include "argduma.h"
#include "argsts.h"
#include "arg_atm_fields.h"
#include "argcona.h"
#include "arglndm.h"
#include "argppx.h"
     &    ErrorStatus,Cmessage)

      ENDIF      ! Diagnostics required for this section

!----------------------------------------------------------------------
!   3. Process physics-derived diagnostics (section 16):
!    3a. Check validity of diagnostic requests (avoids re-checking
!        within every later call of St_diag2).
!    3b. Call St_diag2 as interface to: Phy_diag calculations of
!        physics derived variables and STASH (section 16).

      IF(      SF(0,16)                                                 &
                                ! Diagnostics required for this section
     &   .AND. ErrorStatus == 0) THEN

! DEPENDS ON: st_diag2
      CALL St_diag2( STASH_MAXLEN(16,im_index),                         &
#include "argd1.h"
#include "argduma.h"
#include "argsts.h"
#include "arg_atm_fields.h"
#include "argcona.h"
#include "arglndm.h"
#include "argppx.h"
     &    ErrorStatus,Cmessage)

      ENDIF      ! Diagnostics required for this section

!----------------------------------------------------------------------

! Check error condition
      IF(ErrorStatus >  0) THEN
! DEPENDS ON: ereport
         CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

      RETURN
      END SUBROUTINE InitDiag
#endif
