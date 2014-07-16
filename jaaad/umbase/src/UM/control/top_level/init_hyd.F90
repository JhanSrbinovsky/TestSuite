#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE INIT_HYD ----------------------------------------------
!LL
!LL CW,D.GREGORY<-PROGRAMMERS OF SOME OR ALL OF PREVIOUS CODE OR CHANGES
!LL
!LL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
!LL VERSION  DATE
!LL  3.2  13/04/93  DYNAMIC ALLOCATION OF MAIN ARRAYS. R T H BARNES.
!LL  4.4  14/10/97  Deal with tiled canopy water content.  Richard Betts
!LL  4.5  19/01/98  Replace JVEG_FLDS(6) with JSURF_CAP. D. Robinson.
!LL  5.1  22/02/00  Add PARPARM for TYPSIZE                 P.Burton
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
!LL  VERSION NO. 1
!LL
!LL  SYSTEMS TASK ##
!LL
!LL  PURPOSE : TO SET CANOPY WATER TO ZERO WHEN CANOPY CAPACITY IS ZERO
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P###
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE INIT_HYD(                                              &
#include "argd1.h"
#include "argptra.h"
     &                    ICODE,CMESSAGE)
!
      IMPLICIT NONE
#include "cmaxsize.h"
#include "parparm.h"
#include "typsize.h"
#include "nstypes.h"
#include "typd1.h"
#include "typptra.h"

      LOGICAL L_VEG_FRACS  ! TRUE if tiled land surface scheme in use.
      INTEGER                                                           &
     &       ICODE          ! RETURN CODE : 0 NORMAL EXIT, >0 ERROR
      CHARACTER*(*)                                                     &
     &       CMESSAGE       ! ERROR MESSAGE IF ICODE >0

!
!----------------------------------------------------------------------
! INTERNAL LOOP COUNTERS
!----------------------------------------------------------------------
!
      INTEGER                                                           &
     & I                   ! Loop counter for land points
!
!*---------------------------------------------------------------------
!
      DO I=1,LAND_FIELD
       IF ( D1(JCANOPY_WATER+I-1)  >   D1(JSURF_CAP+I-1) )              &
     &                    D1(JCANOPY_WATER+I-1) = D1(JSURF_CAP+I-1)
      END DO
!
      RETURN
      END SUBROUTINE INIT_HYD
#endif
