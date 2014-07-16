#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE HINTMO -------------------------------------------------
!LL
!LL  Purpose :
!LL
!LL     Performs bilinear horizontal interpolation of a field
!LL     on the model grid to observation locations.
!LL
!LL  For use on Cray Y-MP
!LL  For Cray - Global  ; Enable defs GLOBAL
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LLEND------------------------------------------------------------------
!
!*L  Arguments:---------------------------------------------------
      SUBROUTINE HINTMO(FIELD,CF1PT,CF2PT,CF3PT,CF4PT,                  &
     &                  NP1PT,NP2PT,NP3PT,NP4PT,                        &
     &                  LENFLD,NLEV,LENOBT,BACK,                        &
     &                  ICODE,CMESSAGE)

      IMPLICIT NONE
!-----------------------------------------------------------------
!  Analysis Correction comdecks
!-----------------------------------------------------------------
#include "acparm.h"
#include "comacp.h"
!-----------------------------------------------------------------------
      INTEGER LENFLD,NLEV,LENOBT
      REAL    FIELD(LENFLD*NLEV),CF1PT(LENOBT),CF2PT(LENOBT),           &
     &        CF3PT(LENOBT),CF4PT(LENOBT),BACK(LENOBT)
      INTEGER NP1PT(LENOBT),NP2PT(LENOBT),NP3PT(LENOBT),NP4PT(LENOBT)
      INTEGER ICODE
      CHARACTER*256 CMESSAGE
!
!-INTENT=IN--------------------------------------------------------
!
!     FIELD        - model field to be interpolated
!     LENFLD       - length of one level of field
!     NLEV         - no of levels in field
!     LENOBT       - no of observations
!     CF1,2,3,4PT  - interpolation coefficients for 4 field points
!     NP1,2,3,4PT  - pointers to 4 field points for interpolation
!
!-INTENT=OUT-----------------------------------------------------
!
!     BACK            -  values of field interpolated to obs locations
!     ICODE,CMESSAGE  -  error code and message
!*
!----------------------------------------------------------------------
!*L   Workspace usage
!-----------------------------------------------------------------------
!     NONE
!*
!----------------------------------------------------------------------
!*L   External subroutine calls
!-----------------------------------------------------------------------
      EXTERNAL TIMER
!*
!----------------------------------------------------------------------
!     Define local variables
!----------------------------------------------------------------------
      INTEGER JOB
!
!     JOB   -  Loop counter in loop over obs
!-----------------------------------------------------------------------
!
!L--- 1.     HORIZONTAL INTERPOLATION OF MODEL
!
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('HINTMO  ',3)
!
!     LOOP OVER OBSERVATIONS
!     ORDER OF POINTS IS NEAREST POINT (I,J), (I+/-1,J)
!                                      (I,J+/-1), (I+/-1,J+/-1)
      DO 1000 JOB=1,LENOBT
!
      BACK(JOB) = CF1PT(JOB) *  FIELD ( NP1PT(JOB) ) +                  &
     &            CF2PT(JOB) *  FIELD ( NP2PT(JOB) ) +                  &
     &            CF3PT(JOB) *  FIELD ( NP3PT(JOB) ) +                  &
     &            CF4PT(JOB) *  FIELD ( NP4PT(JOB) )
!
1000  CONTINUE
!
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('HINTMO  ',4)
!
      RETURN
      END SUBROUTINE HINTMO
#endif
