#if defined(CONTROL)
#if defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine INITMOS ------------------------------------------------
!LL
!LL Purpose: Control routine for initialisation of MOS grid information.
!LL          Called from model initialisation routine INITIAL.
!LL
!LL Version 1 Author:   T.Johns
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL  3.2  13/04/93  Dynamic allocation of main arrays. R T H Barnes.
!    4.2  10/01/97  For MPP runs, skip call to MOSGRID and initialise
!                   MOS_OUTPUT_LENGTH. Problem to be addressed at 4.3.
!                   D.Robinson
!LL  4.3  25/02/97  Modset GDR1F402 removed, MPP MOS code now operative
!LL                 Supply global sizes to MOS_GRID          P.Burton
!LL  5.0  21/06/99  Corrections to allow compilation under cpp. Note
!LL                 that MOS functionality is no longer required, but
!LL                 needs to be retained for later integration with
!LL                 FSSSI capability. R Rawlins
!  5.4     28/08/02    Bug Fix (Bi-cyclic LAM)           Carol Roadnight
!    5.4  21/03/02  Remove comment on same line as #include
!                                                S. Carroll
!LL
!LL Programming standard : UM Doc Paper no 3
!LL
!LL Logical components covered : C4
!LL
!LL Project task : C4
!LL
!LL External documentation : UMDP no C4
!LL
!L* Interface and arguments --------------------------------------------
!
      SUBROUTINE INITMOS(                                               &
#include "argduma.h"
#include "argsts.h"
     &                   ICODE,CMESSAGE)
!
      IMPLICIT NONE
!
#include "cmaxsize.h"
#include "csubmodl.h"
#include "parvars.h"
#include "typsize.h"
#include "typduma.h"
#include "typsts.h"

      INTEGER                                                           &
     &    ICODE          !OUT     RETURN CODE
      CHARACTER*80                                                      &
     &    CMESSAGE       !OUT     ANY ERROR MESSAGE PASSED BACK
!*---------------------------------------------------------------------
! Common blocks and PARAMETERs
!
! For Model_domain
#include "cntlatm.h"
! Print status information
#include "cprintst.h"
!
! Subroutines called
!
      EXTERNAL MOSGRID
!
! Local variables.
!
      LOGICAL                                                           &
     & Global_domain                                                    &
                       ! T if global model
     &,LAM_domain      ! T if limited area model

       Global_domain = (Model_domain == mt_global)
       LAM_domain    = (Model_domain == mt_lam .OR.                     &
     &                  Model_domain == mt_cyclic_lam .or.              &
     &                  Model_domain == mt_bi_cyclic_LAM)

!L
!L  Set up masking array for MOS OUTPUT Ie MOS_OUTPUT
!L
! DEPENDS ON: mosgrid
      CALL MOSGRID(Global_domain,LAM_domain,                            &
     &             A_REALHD(rh_rotlat), A_REALHD(rh_rotlong),           &
     &             A_REALHD(rh_deltaNS),A_REALHD(rh_deltaEW),           &
     &             A_REALHD(rh_baselat),A_REALHD(rh_baselong),          &
     &             glsize(1,fld_type_p),glsize(2,fld_type_p),           &
     &             MOS_MASK,MOS_OUTPUT_LENGTH,ICODE,CMESSAGE)

       IF(PrintStatus >= PrStatus_Oper) THEN
      WRITE(6,*)'INITMOS : MOS_OUTPUT_LENGTH = ',MOS_OUTPUT_LENGTH
       ENDIF  ! PrintStatus test

      RETURN
      END SUBROUTINE INITMOS
#endif
#endif
