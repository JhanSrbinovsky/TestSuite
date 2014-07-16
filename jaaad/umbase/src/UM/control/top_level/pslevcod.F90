#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Set the STASH addresses for D1
! Subroutine Interface:

!- End of subroutine code -------------------------------------------


!+Compute data lengths and addresses for primary fields
! Subroutine Interface:

!+Test whether level type is discrete (model) or continuous (non-model)
! Function Interface:
!- End of Function code --------------------------------------------

!+Decode the STASH pseudo level code
! Subroutine Interface:
      SUBROUTINE PSLEVCOD(ILIN,ILOUT,SWTCH,ErrorStatus,CMESSAGE)
      USE CSENARIO_MOD
      IMPLICIT NONE
! Description:
!   Sets ILOUT to an appropriate pseudo level size according
!    to the value of IL
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
#include "model.h"
#include "parparm.h"
#include "typsize.h"
#include "cntlatm.h"
#include "nstypes.h"

! Subroutine arguments:
!   Scalar arguments with intent(in):
      INTEGER ILIN    ! Model pseudo level code
      CHARACTER*1 SWTCH

!   Scalar arguments with intent(out):
      INTEGER ILOUT   ! An actual pseudo level
      CHARACTER*80 CMESSAGE

! Local scalars:
      INTEGER I
      INTEGER J

! Error Status:
      INTEGER ErrorStatus

!- End of Header --------------------------------------------------

      IF (SWTCH == 'F') THEN
        IF(ILIN == 1) THEN
          ILOUT=1
! Ocean assimilation groups
        ELSE IF(ILIN == 41) THEN
          ILOUT=OASLEV(1)
        ELSE IF(ILIN == 42) THEN
          ILOUT=OASLEV(2)
        ELSE IF(ILIN == 43) THEN
          ILOUT=OASLEV(3)
        ELSE IF(ILIN == 44) THEN
          ILOUT=OASLEV(4)
        ELSE IF(ILIN == 45) THEN
          ILOUT=OASLEV(5)
        ELSE IF(ILIN == 46) THEN
          ILOUT=OASLEV(6)
        ELSE
          WRITE(6,*)                                                    &
     &   'MSG FROM PSLEVCOD: ',                                         &
     &   'INAPPROPRIATE FIRST PSEUDO LEVEL CODE FOUND ',ILIN
          ErrorStatus=2
        END IF
      ELSE IF (SWTCH == 'L') THEN
        IF(ILIN == 1) THEN
          ILOUT=H_SWBANDS
        ELSE IF(ILIN == 2) THEN
          ILOUT=H_LWBANDS
        ELSEIF ( ILIN  ==  6 ) THEN
! Last index for HadCM2 sulphate loading patterns.
          ILOUT = NSULPAT
        ELSEIF ( ILIN  ==  7 ) THEN
! All surface types
          ILOUT = NTYPE
        ELSEIF ( ILIN  ==  8 ) THEN
! Plant functional types only
          ILOUT = NPFT
        ELSEIF ( ILIN  ==  9 ) THEN
! All tiles
          ILOUT = NTILES
        ELSEIF ( ILIN  ==  10 )THEN
! All sea ice catagories
          ILOUT = NICE
        ELSE
          WRITE(6,*)                                                    &
     &   'MSG FROM PSLEVCOD: ',                                         &
     &   'INAPPROPRIATE LAST PSEUDO LEVEL CODE FOUND ',ILIN
          ErrorStatus=2
        END IF

      END IF

      RETURN
      END SUBROUTINE PSLEVCOD
#endif
