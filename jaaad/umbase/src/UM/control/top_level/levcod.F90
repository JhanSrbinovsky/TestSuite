#if defined(CONTROL) || defined(MAKEBC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Decode the STASH level code

!******************************************************************
! WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
! Any Changes to this routine must be accompanied with equivalent
! changes to the deck RCF_LEVELCODE
!******************************************************************

! Subroutine Interface:
      SUBROUTINE LEVCOD(ILIN,ILOUT,ErrorStatus,CMESSAGE)
      IMPLICIT NONE
! Description:
!   Sets ILOUT to an appropriate level size according to the value of IL
!   Level sizes are parametrised in comdeck MODEL.
!
! Current code owner:  S.J.Swarbrick
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
! Global variables:

#include "csubmodl.h"
#include "version.h"
#include "parparm.h"
#include "typsize.h"
#include "model.h"
#include "cntlatm.h"

! Subroutine arguments:
!   Scalar arguments with intent(in):
      INTEGER ILIN    ! Model level code

!   Scalar arguments with intent(out):
      INTEGER ILOUT   ! An actual level
      CHARACTER*80 CMESSAGE

! Local scalars:
      INTEGER I
      INTEGER J

! ErrorStatus:
      INTEGER ErrorStatus

!- End of Header ---------------------------------------------------

      IF(ILIN == 1) THEN
! First atmos level
        ILOUT=1
      ELSE IF(ILIN == 2) THEN
! Top atmos level
        ILOUT=MODEL_LEVELS
      ELSE IF(ILIN == 3) THEN
! Top wet level
        ILOUT=WET_LEVELS
      ELSE IF(ILIN == 4) THEN
        ILOUT=MODEL_LEVELS-1
      ELSE IF(ILIN == 5) THEN
! First boundary layer level
        ILOUT=MIN(1,BL_LEVELS)
      ELSE IF(ILIN == 6) THEN
! Last boundary layer level
        ILOUT=BL_LEVELS
      ELSE IF(ILIN == 7) THEN
        ILOUT=BL_LEVELS+1
      ELSE IF(ILIN == 8) THEN
! First soil level
        ILOUT=MIN(1,ST_LEVELS)
      ELSE IF(ILIN == 9) THEN
! Last soil level
        ILOUT=ST_LEVELS
      ELSE IF(ILIN == 10) THEN
! First tracer level
        ILOUT=MODEL_LEVELS-TR_LEVELS+1
      ELSE IF(ILIN == 11) THEN
! Last tracer level
        ILOUT=MODEL_LEVELS
      ELSE IF(ILIN == 12) THEN
        ILOUT=MODEL_LEVELS+1
      ELSE IF(ILIN == 13) THEN
! First gravity wave drag level
        ILOUT=StLevGWdrag
      ELSE IF(ILIN == 14) THEN
! Last gravity wave drag level
        ILOUT=MODEL_LEVELS
      ELSE IF(ILIN == 15) THEN
        ILOUT=BotVDiffLev
      ELSE IF(ILIN == 16) THEN
        ILOUT=TopVDiffLev-1
      ELSE IF(ILIN == 17) THEN
        ILOUT=TopVDiffLev
      ELSE IF(ILIN == 18) THEN
        ILOUT=BL_LEVELS-1
      ELSE IF(ILIN == 19) THEN
        ILOUT=MODEL_LEVELS+1
      ELSE IF(ILIN == 20) THEN
        ILOUT=MIN(2,ST_LEVELS)
      ELSE IF(ILIN == 21) THEN
        ILOUT=1
      ELSE IF(ILIN == 22) THEN
        ILOUT=NLEVSO
      ELSE IF(ILIN == 23) THEN
        ILOUT=OZONE_LEVELS
      ELSE IF(ILIN == 24) THEN
        ILOUT=MODEL_LEVELS*H_SWBANDS
      ELSE IF(ILIN == 25) THEN
        ILOUT=(MODEL_LEVELS+1)*H_SWBANDS
      ELSE IF(ILIN == 26) THEN
        ILOUT=WET_LEVELS*H_SWBANDS
      ELSE IF(ILIN == 27) THEN
        ILOUT=MODEL_LEVELS*H_LWBANDS
      ELSE IF(ILIN == 28) THEN
        ILOUT=(MODEL_LEVELS+1)*H_LWBANDS
      ELSE IF(ILIN == 29) THEN
        ILOUT=WET_LEVELS*H_LWBANDS
      ELSE IF(ILIN == 30) THEN
        ILOUT=2
      ELSE IF(ILIN == 32) THEN
        ILOUT=H_SWBANDS
      ELSE IF(ILIN == 33) THEN
        ILOUT=H_LWBANDS
      ELSE IF(ILIN == 34) THEN
        ILOUT=SM_LEVELS
      ELSE IF(ILIN == 35) THEN
        ILOUT=CLOUD_LEVELS
      ELSE IF(ILIN == 38) THEN
! Surface level on Charney-Phillips theta grid
        ILOUT=0
      ELSE IF(ILIN == 39) THEN
! Number of ISCCP simulator levels
        ILOUT=7
      ELSE IF(ILIN == 41) THEN
        ILOUT=OASLEV(1)
!       Allow room for expansion of ocean assimilation groups.
      ELSE IF(ILIN == 42) THEN
        ILOUT=OASLEV(2)
!       Allow room for expansion of ocean assimilation groups.
      ELSE IF(ILIN == 43) THEN
        ILOUT=OASLEV(3)
!       Allow room for expansion of ocean assimilation groups.
      ELSE IF(ILIN == 44) THEN
        ILOUT=OASLEV(4)
!       Allow room for expansion of ocean assimilation groups.
      ELSE IF(ILIN == 45) THEN
        ILOUT=OASLEV(5)
!       Allow room for expansion of ocean assimilation groups.
      ELSE IF(ILIN == 46) THEN
        ILOUT=OASLEV(6)
!       Allow room for expansion of ocean assimilation groups.
      ELSE
        WRITE(6,*)'LEVCOD: IMPOSSIBLE LEVEL CODE FOUND ',ILIN
        ErrorStatus=1
      END IF

      RETURN
      END SUBROUTINE LEVCOD
#endif
