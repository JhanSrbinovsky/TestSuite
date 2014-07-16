#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE HMRTORH ------------------------------------------------
!LL
!LL  Purpose : Convert Humidity Mixing Ratio to Relative Humidity
!LL            and vice versa.  Option to re-initialise cloud water
!LL            from rh and temp.
!LL
!LL            KRMODE = 1 - For HMR to RH
!LL            KRMODE = 2 - For RH  to HMR
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL
!LL   3.1  10/02/93    Hardwired calls to TIMER removed
!LL                    Author: A. Dickinson    Reviewer: F. Rawlins
!LL   3.4  07/09/94    Remove cloud water/ice and KRMODE=3 option
!LL                    Remove RHCRIT from arg list & ref to REINITQC
!LL                                                      Bruce M
!LL   3.4  19/09/94    Make available to A18_2A    S Bell
!LL   5.2  30/11/00    revise for new dynamics    B Macpherson
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LLEND------------------------------------------------------------------
!*L   Arguments
      SUBROUTINE HMRTORH (KRMODE, EXNER, PRESSURE, THETA, RH,           &
     &                    P_FIELD, P_LEVELS, Q_LEVELS, ICODE, CMESSAGE)
      IMPLICIT NONE
#include "acparm.h"
#include "comacp.h"
!-----------------------------------------------------------------------
      INTEGER KRMODE
      INTEGER P_FIELD,P_LEVELS,Q_LEVELS
      REAL EXNER (P_FIELD,P_LEVELS)
      REAL PRESSURE (P_FIELD,P_LEVELS)
      REAL THETA (P_FIELD,P_LEVELS)
      REAL RH    (P_FIELD,Q_LEVELS)

      INTEGER ICODE
      CHARACTER*256 CMESSAGE
!=======================================================================
!     UM Constant comdeck
#include "c_r_cp.h"
!-----------------------------------------------------------------------
!     Dynamic allocation
      REAL TEMP     (P_FIELD)
      REAL SMR      (P_FIELD)
!=======================================================================
!     Local variables
      INTEGER JLEV,J

!=======================================================================
      EXTERNAL QSAT,TIMER
!=======================================================================


      DO JLEV=1,Q_LEVELS

!       Convert Theta to Temperature
!       ----------------------------
        DO J=1,P_FIELD
          TEMP(J) = THETA(J,JLEV) * EXNER(J,JLEV)
        ENDDO


!       Obtain Saturated Mixing Ratio
!       -----------------------------
! DEPENDS ON: qsat
        CALL QSAT (SMR,TEMP,PRESSURE(1,JLEV),P_FIELD)


        IF (KRMODE == 1) THEN

!         Convert Mixing Ratio to Relative Humidity
!         -----------------------------------------
          DO J=1,P_FIELD
            RH(J,JLEV) = RH(J,JLEV)/SMR(J)
            RH(J,JLEV) = RH(J,JLEV)*100.0
          ENDDO

        ELSEIF (KRMODE == 2) THEN

!         Convert Relative Humidity to Mixing Ratio
!         -----------------------------------------
          DO J=1,P_FIELD
            RH(J,JLEV) = RH(J,JLEV)*SMR(J)
            RH(J,JLEV) = RH(J,JLEV)*0.01
          ENDDO


        ENDIF

      ENDDO

      RETURN
      END SUBROUTINE HMRTORH
#endif
