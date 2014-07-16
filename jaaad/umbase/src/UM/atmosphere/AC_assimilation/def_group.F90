#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE DEF_GROUP ----------------------------------------------
!LL
!LL  Purpose : Initialise Group Dependent arrays
!LL            in comdeck COMACP and COMAG.
!LL
!LL  For Global runs : Enable defs GLOBAL
!LL
!LL  Written 30/3/92 by Dave Robinson.
!LL
!LL  Logical components covered:
!LL
!LL  Project Task : P3
!LL
!LL  Documentation:
!LL
!LLEND------------------------------------------------------------------
!
!*L  ARGUMENTS :--------------------------------------------------------
      SUBROUTINE DEF_GROUP (P_LEVELS,Q_LEVELS,BL_LEVELS,TR_LEVELS,      &
     &                      ICODE,CMESSAGE)

      IMPLICIT NONE

      INTEGER                                                           &
     &   P_LEVELS                                                       &
                       ! IN - No of Model Levels.
     & , Q_LEVELS                                                       &
                       ! IN - No of Model Wet Levels.
     & , BL_LEVELS                                                      &
                       ! IN - No of levels in boundary layer.
     & , TR_LEVELS                                                      &
                       ! IN - No of Tracer levels.
     & , ICODE         ! OUT - Return Code
      CHARACTER*256 CMESSAGE  ! OUT - Error Message

!     AC Comdecks
#include "acparm.h"
#include "comacp.h"
#include "comag.h"
!     UM Comdecks
#include "c_mdi.h"

!*L  Workspace Usage ---------------------------------------------------
!     Local array.
      INTEGER AC_GROUPS (NOBTYPMX)  ! Stores default groups of AC Types.
!     Local variables.
      INTEGER JOBT  !  Loop counter.
      INTEGER I     !  group identifier

!     Order of processing and grouping of Observation Types
!     =====================================================

!     Numbers = (Group No*1000) + Obs Type

!L Default groupings are:
!L 1:pstar
!L 2:upper level temperatures
!L 3:surface temperatures
!L 4:upper level winds
!L 5:surface winds
!L 6:upper level RH
!L 7:surface RH
!L 8:MOPS RH
!L 9:Cloud histograms
!L 10:MOPS precip rate/phase
!L 11:Tracers (NB. when the assimilation is actually run, each tracer
!L    must be in a separate group; they are only put together here for
!L    convenience, to define common assimilation parameters.
!L    Normally only a few would be used at once.)

      DATA AC_GROUPS /                                                  &
     &  1101,                                                           &
     &  2201,2203,2205,2206,2207,2208,2209,2211,                        &
     &  3202,3204,                                                      &
     &  4301,4303,4311,                                                 &
     &  5302,5304,5305,5306,                                            &
     &  6401,6403,6405,                                                 &
     &  7402,7404,                                                      &
     &  8406,                                                           &
     &  9407,                                                           &
     &  10506,                                                          &
     &  11601, 11602, 11603, 11604, 11605,                              &
     &  11606, 11607, 11608, 11609, 11610,                              &
     &  11611, 11612, 11613, 11614, 11615,                              &
     &  11616, 11617, 11618, 11619, 11620,                              &
     &  11621, 11622, 11623, 11624, 11625,                              &
     &  11626, 11627, 11628, 11629,                                     &
     &  12901, 126*0/
!     NB. the number of tracers (group 11) should correspond to
!     A_MAX_TRVARS (in CTRACERA), currently 29

!     Copy above list into COMACP array DEF_AC_ORDER
      DO JOBT=1,NOBTYPMX
        DEF_AC_ORDER(JOBT) = AC_GROUPS(JOBT)
      ENDDO

!     Group Dependent Variables
!     =========================

!     DEF_NO_ANAL_LEVS  No of analysis levels
!     DEF_NO_WT_LEVS    No of weight levels
!     DEF_NO_ITERATIONS No of iterations
!     DEF_INTERVAL_ITER Interval in timesteps between iterations
!     DEF_AGRES_ROWS Ratio of No of rows in Model Grid to Analysis Grid
!     DEF_AGRES_PTS  Ratio of No of pts in Model Grid to Analysis Grid
!     DEF_MODE_HANAL Mode of Horizontal Analysis
!     DEF_FI_VAR_FACTOR group dep scaling factor in FI
#if defined(GLOBAL)
!     DEF_NUDGE_NH   Nudging Coefficients for NH
!     DEF_NUDGE_TR   Nudging Coefficients for TR
!     DEF_NUDGE_SH   Nudging Coefficients for SH
#else
!     DEF_NUDGE_LAM  Nudging Coefficients for LAM
#endif

#if defined(GLOBAL)
      I = 1                             ! pstar
      DEF_NO_ANAL_LEVS(I)  = 1
      DEF_NO_WT_LEVS(I)    = 1
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 4.0E-4
        DEF_NUDGE_TR(I) = 4.0E-4
        DEF_NUDGE_SH(I) = 4.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 5.0E-4
        DEF_NUDGE_TR(I) = 3.0E-4
        DEF_NUDGE_SH(I) = 3.8E-4
      ENDIF

      I = 2                             ! upper level temps
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NO_ANAL_LEVS(I)  = P_LEVELS
      DEF_NO_WT_LEVS(I)    = P_LEVELS
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 4.0E-4
        DEF_NUDGE_TR(I) = 4.0E-4
        DEF_NUDGE_SH(I) = 4.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 5.0E-4
        DEF_NUDGE_TR(I) = 3.0E-4
        DEF_NUDGE_SH(I) = 3.8E-4
      ENDIF

      I = 3                             ! surf temps
      DEF_NO_ANAL_LEVS(I)  = MAX(BL_LEVELS-2,3)
      DEF_NO_WT_LEVS(I)    = 1
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 4.0E-4
        DEF_NUDGE_TR(I) = 4.0E-4
        DEF_NUDGE_SH(I) = 4.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 5.0E-4
        DEF_NUDGE_TR(I) = 3.0E-4
        DEF_NUDGE_SH(I) = 3.8E-4
      ENDIF

      I = 4                             !upper level winds
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35 * 3.
      DEF_NO_ANAL_LEVS(I)  = P_LEVELS
      DEF_NO_WT_LEVS(I)    = P_LEVELS
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 6.0E-4
        DEF_NUDGE_TR(I) = 6.0E-4
        DEF_NUDGE_SH(I) = 6.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 6.6E-4
        DEF_NUDGE_TR(I) = 4.0E-4
        DEF_NUDGE_SH(I) = 4.3E-4
      ENDIF

      I = 5                             ! surf winds (inc scat)
      DEF_NO_ANAL_LEVS(I)  = BL_LEVELS
      DEF_NO_WT_LEVS(I)    = 1
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35 * 3.
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 6.0E-4
        DEF_NUDGE_TR(I) = 6.0E-4
        DEF_NUDGE_SH(I) = 6.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 6.3E-4
        DEF_NUDGE_TR(I) = 3.8E-4

        DEF_NUDGE_SH(I) = 4.0E-4
      ENDIF

      I = 6                             ! upper level RH
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NO_ANAL_LEVS(I)  = MIN(Q_LEVELS,P_LEVELS)
      DEF_NO_WT_LEVS(I)    = MIN(Q_LEVELS,P_LEVELS)
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 4.0E-4
        DEF_NUDGE_TR(I) = 4.0E-4
        DEF_NUDGE_SH(I) = 4.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 5.0E-4
        DEF_NUDGE_TR(I) = 3.0E-4
        DEF_NUDGE_SH(I) = 3.5E-4
      ENDIF

      I = 7                             ! surface RH
      DEF_NO_ANAL_LEVS(I)  = MAX(BL_LEVELS-2,3)
      DEF_NO_WT_LEVS(I)    = 1
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 4.0E-4
        DEF_NUDGE_TR(I) = 4.0E-4
        DEF_NUDGE_SH(I) = 4.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 5.0E-4
        DEF_NUDGE_TR(I) = 3.0E-4
        DEF_NUDGE_SH(I) = 3.5E-4
      ENDIF

      I = 8                             ! MOPS RH
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 2
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NO_ANAL_LEVS(I)  = MIN(Q_LEVELS,P_LEVELS)
      DEF_NO_WT_LEVS(I)    = MIN(Q_LEVELS,P_LEVELS)
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 4.0E-4
        DEF_NUDGE_TR(I) = 4.0E-4
        DEF_NUDGE_SH(I) = 4.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 5.0E-4
        DEF_NUDGE_TR(I) = 3.0E-4
        DEF_NUDGE_SH(I) = 3.5E-4
      ENDIF

      I = 9                             ! Cloud histograms
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 2
      DEF_FI_VAR_FACTOR(I) = 2.0
      DEF_NO_ANAL_LEVS(I)  = MIN(Q_LEVELS,P_LEVELS)
      DEF_NO_WT_LEVS(I)    = MIN(Q_LEVELS,P_LEVELS)
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 4.0E-4
        DEF_NUDGE_TR(I) = 4.0E-4
        DEF_NUDGE_SH(I) = 4.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 5.0E-4
        DEF_NUDGE_TR(I) = 3.0E-4
        DEF_NUDGE_SH(I) = 3.5E-4
      ENDIF

      I = 10                  ! MOPS precip rate/phase
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 2
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NO_ANAL_LEVS(I)  = 1
      DEF_NO_WT_LEVS(I)    = 1
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 4.0E-4
        DEF_NUDGE_TR(I) = 4.0E-4
        DEF_NUDGE_SH(I) = 4.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 1.0E6
        DEF_NUDGE_TR(I) = 1.0E6
        DEF_NUDGE_SH(I) = 1.0E6
      ENDIF

      I = 11                            ! tracers
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NO_ANAL_LEVS(I)  = TR_LEVELS
      DEF_NO_WT_LEVS(I)    = TR_LEVELS
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 4.0E-4
        DEF_NUDGE_TR(I) = 4.0E-4
        DEF_NUDGE_SH(I) = 4.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 5.0E-4
        DEF_NUDGE_TR(I) = 3.0E-4
        DEF_NUDGE_SH(I) = 3.8E-4
      ENDIF

      I = 12                            ! surface LOG Visibility
      DEF_NO_ANAL_LEVS(I)  = BL_LEVELS
      DEF_NO_WT_LEVS(I)    = 1
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 3.0E-4
        DEF_NUDGE_TR(I) = 3.0E-4
        DEF_NUDGE_SH(I) = 3.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 5.0E-4
        DEF_NUDGE_TR(I) = 3.0E-4
        DEF_NUDGE_SH(I) = 3.5E-4
      ENDIF

#else
      I = 1                             !pstar
      DEF_NO_ANAL_LEVS(I)  = 1
      DEF_NO_WT_LEVS(I)    = 1
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NUDGE_LAM(I)   = 5.0E-4

      I = 2                             !upper level temps
      DEF_NO_ANAL_LEVS(I)  = P_LEVELS
      DEF_NO_WT_LEVS(I)    = P_LEVELS
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NUDGE_LAM(I)   = 5.0E-4

      I = 3                             !surf temps
      DEF_NO_ANAL_LEVS(I)  = MAX(BL_LEVELS-2,3)
      DEF_NO_WT_LEVS(I)    = 1
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NUDGE_LAM(I)   = 5.0E-4
      IF (LAC_MES) THEN
        DEF_NO_ANAL_LEVS(I)  = 6
      ENDIF

      I = 4                             !upper level winds
      DEF_NO_ANAL_LEVS(I)  = P_LEVELS
      DEF_NO_WT_LEVS(I)    = P_LEVELS
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35 * 3.
      DEF_NUDGE_LAM(I)   = 6.6E-4

      I = 5                            !surf winds (inc scat)
      DEF_NO_ANAL_LEVS(I)  = BL_LEVELS
      DEF_NO_WT_LEVS(I)    = 1
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35 * 3.
      DEF_NUDGE_LAM(I)   = 6.3E-4
      IF (LAC_MES) THEN
        DEF_NO_ANAL_LEVS(I)  = 6
      ENDIF

      I = 6                            !upper level RH
      DEF_NO_ANAL_LEVS(I)  = MIN(Q_LEVELS,P_LEVELS)
      DEF_NO_WT_LEVS(I)    = MIN(Q_LEVELS,P_LEVELS)
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NUDGE_LAM(I)   = 5.0E-4

      I = 7                            !surface RH
      DEF_NO_ANAL_LEVS(I)  = MAX(BL_LEVELS-2,3)
      DEF_NO_WT_LEVS(I)    = 1
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NUDGE_LAM(I)   = 5.0E-4
      IF (LAC_MES) THEN
        DEF_NO_ANAL_LEVS(I)  = 6
      ENDIF

      I = 8                            !MOPS RH
      DEF_NO_ANAL_LEVS(I)  = MIN(Q_LEVELS,P_LEVELS)
      DEF_NO_WT_LEVS(I)    = MIN(Q_LEVELS,P_LEVELS)
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 2
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NUDGE_LAM(I)   = 1.0E-3

      I = 9                             ! Cloud histograms
      DEF_NO_ANAL_LEVS(I)  = MIN(Q_LEVELS,P_LEVELS)
      DEF_NO_WT_LEVS(I)    = MIN(Q_LEVELS,P_LEVELS)
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 2
      DEF_FI_VAR_FACTOR(I) = 2.0
      DEF_NUDGE_LAM(I)   = 5.0E-4

      I = 10                  ! MOPS precip rate/phase
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 2
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NO_ANAL_LEVS(I)  = 1
      DEF_NO_WT_LEVS(I)    = 1
      DEF_NUDGE_LAM(I)     = 1.0E6

      I = 11                            ! tracers
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NO_ANAL_LEVS(I)  = TR_LEVELS
      DEF_NO_WT_LEVS(I)    = TR_LEVELS
      DEF_NUDGE_LAM(I)  = 5.0E-4

      I = 12                           !surface LOG Visibility
      DEF_NO_ANAL_LEVS(I)  = BL_LEVELS
      DEF_NO_WT_LEVS(I)    = 1
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NUDGE_LAM(I)   = 5.0E-4
      IF (LAC_MES) THEN
        DEF_NO_ANAL_LEVS(I)  = 6
      ENDIF
#endif

      DO JOBT=I+1,NOBTYPMX
        DEF_NO_ANAL_LEVS(JOBT)  = IMDI
        DEF_NO_WT_LEVS(JOBT)    = IMDI
        DEF_NO_ITERATIONS(JOBT) = IMDI
        DEF_INTERVAL_ITER(JOBT) = IMDI
        DEF_AGRES_ROWS(JOBT)    = IMDI
        DEF_AGRES_PTS(JOBT)     = IMDI
        DEF_MODE_HANAL(JOBT)    = IMDI
        DEF_FI_VAR_FACTOR(JOBT) = RMDI
#if defined(GLOBAL)
        DEF_NUDGE_NH(JOBT) = RMDI
        DEF_NUDGE_TR(JOBT) = RMDI
        DEF_NUDGE_SH(JOBT) = RMDI

#else
        DEF_NUDGE_LAM(JOBT) = RMDI
#endif
      ENDDO


      RETURN
      END SUBROUTINE DEF_GROUP
#endif
