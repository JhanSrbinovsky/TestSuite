#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE DEF_TYPE ----------------------------------------------
!LL
!LL  Purpose : Initialise Observation Type Dependent Arrays
!LL            in comdeck COMACP
!LL
!LL  For use on Cray
!LL  For Global runs : Enable defs GLOBAL
!LL
!LL  Written 30/3/92 by Dave Robinson.
!LL
!LL SB, RS      <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL 3.1    18/2/93  : Set new defaults under LAC_MES. (B Macpherson)
!LL 3.2     8/7/93  : tidy up rh radinf values
!LL                 : add type 306 (drifter winds)
!LL                 : and add OBTHIN                 (S Bell)
!LL 3.3    1/12/93  : Revise mes defaults         B Macpherson
!LL 3.3    11/11/93 : Cater for type 407 (cloud histograms)  N Richards
!LL 3.3    1/12/93  : Cater for type 506 (MOPS rain)  B Macpherson
!LL 3.4    03/08/94 : Cater for tracers (type 6nn) R Swinbank
!LL                   (29 types, corresponding to STASH codes 61-89)
!LL 3.4    01/09/94 : Amend UARS defaults R Swinbank
!LL                   (correlation scale, sat120 thinning)
!LL 3.4    9/09/94  : Revise type 506 defaults        C D Jones
!LL                 : Revise mes 401 & 204 defaults   B Macpherson
!LL                 : Set defaults for log vis        P Clark
!LL 4.0    5/07/95  : Revise defaults for synops  SBell
!LL 4.0   30/06/95  : Revise mes defaults for 201,203,205,206,207,209,
!LL                 : 301,302,303,304,306,403   (Bruce M)
!LL 4.0   02/03/95  : Removed comment from nupdate directive
!LL                   line (Tracey Smith)
!LL
!LL  Programming Standard : UM Doc Paper No 3 ; Version 4 ; 5/2/92
!LL
!LL  Project Task : P3
!LL
!LLEND------------------------------------------------------------------
!
!*L  ARGUMENTS:---------------------------------------------------------
      SUBROUTINE DEF_TYPE (ICODE,CMESSAGE)

      IMPLICIT NONE

      INTEGER ICODE           !  Return Code
      CHARACTER*256 CMESSAGE  !  Error Message

!     AC Comdecks
#include "acparm.h"
#include "comacp.h"
#include "ctracera.h"

!*L  WORKSPACE USAGE:---------------------------------------------------
!     Local variables
      INTEGER I,J,JOBT  !  Array index, loop counters.

!     Arrays initialised for each observation type.
!     MASTER_AC_TYPES   - AC Observation Types known in AC Scheme
!     DEF_TIMEB         - Insertion period before observation time
!     DEF_TIMEA         - Insertion period after  observation time
!     DEF_TGETOBB       - Time Window before obs time to fetch obs
!     DEF_TGETOBA       - Time Window after  obs time to fetch obs
!     DEF_RADINF        - Maximum Normalised Influence Radius.
!     DEF_OBTHIN        - Observation thinning factor
!     DEF_CSCALE_START  - Correlation Scale at start of insertion period
!     DEF_CSCALE_OBTIME -                      observation time
!     DEF_CSCALE_END    -                      end of insertion period.

      I=0

#if defined(GLOBAL)

      I = I + 1              !  Defaults for Type 101 pstar
      MASTER_AC_TYPES(I)   = 101
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 201 sonde temp
      MASTER_AC_TYPES(I)   = 201
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 202 ship surf temp
      MASTER_AC_TYPES(I)   = 202
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 203 airep temp
      MASTER_AC_TYPES(I)   = 203
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 204 synop surf temp
      MASTER_AC_TYPES(I)   = 204
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 200.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0

      I = I + 1              !  Defaults for Type 205 LASS temp
      MASTER_AC_TYPES(I)   = 205
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 3
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 206 SATEM temp
      MASTER_AC_TYPES(I)   = 206
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 207 SAT120 temp
      MASTER_AC_TYPES(I)   = 207
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 3
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 208 constrained LASS
      MASTER_AC_TYPES(I)   = 208
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 3
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 209 BOGUS 1000-500
      MASTER_AC_TYPES(I)   = 209
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 211 UARS Temp
      MASTER_AC_TYPES(I)   = 211
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 301 Sonde winds
      MASTER_AC_TYPES(I)   = 301
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 302 Ship surf wind
      MASTER_AC_TYPES(I)   = 302
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 303 Airep wind
      MASTER_AC_TYPES(I)   = 303
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 304 Synop surf wind
      MASTER_AC_TYPES(I)   = 304
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 240.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0

      I = I + 1              !  Defaults for Type 305 Scatwind
      MASTER_AC_TYPES(I)   = 305
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 5
      DEF_CSCALE_START(I)  = 240.0
      DEF_CSCALE_OBTIME(I) = 120.0
      DEF_CSCALE_END(I)    = 120.0

      I = I + 1              !  Defaults for Type 306 Drifter winds
      MASTER_AC_TYPES(I)   = 306
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 311 UARS wind
      MASTER_AC_TYPES(I)   = 311
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 401 Sonde RH
      MASTER_AC_TYPES(I)   = 401
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0

      I = I + 1              !  Defaults for Type 402 Ship surf RH
      MASTER_AC_TYPES(I)   = 402
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0

      I = I + 1              !  Defaults for Type 403 Bogus RH
      MASTER_AC_TYPES(I)   = 403
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0

      I = I + 1              !  Defaults for Type 404 Synop surf RH
      MASTER_AC_TYPES(I)   = 404
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 200.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0

      I = I + 1              !  Defaults for Type 405 LASS RH
      MASTER_AC_TYPES(I)   = 405
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 3
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0

      I = I + 1              !  Defaults for Type 406 MOPS RH
      MASTER_AC_TYPES(I)   = 406
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 0.01
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 999.0
      DEF_CSCALE_OBTIME(I) = 999.0
      DEF_CSCALE_END(I)    = 999.0

      I = I + 1              !  Defaults for Type 407 Cloud histograms
      MASTER_AC_TYPES(I)   = 407
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 100.0
      DEF_CSCALE_OBTIME(I) = 100.0
      DEF_CSCALE_END(I)    = 100.0

      I = I + 1              !  Defaults for type 506 MOPS precip
      MASTER_AC_TYPES(I)   = 506
      DEF_TIMEB(I)         = 180.
      DEF_TIMEA(I)         = 180.
      DEF_TGETOBB(I)       = 180.
      DEF_TGETOBA(I)       = 180.
      DEF_RADINF(I)        = 0.01
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 999.0
      DEF_CSCALE_OBTIME(I) = 999.0
      DEF_CSCALE_END(I)    = 999.0

      !  Defaults for Type 6nn (tracers)
      !  Allow for A_MAX_TRVARS (currently 29) possible tracers
      DO J=1,A_MAX_TRVARS
        I = I + 1
        MASTER_AC_TYPES(I) = 600+J
        DEF_TIMEB(I)       = 240.0
        DEF_TIMEA(I)       = 60.0
        DEF_TGETOBB(I)     = 240.0
        DEF_TGETOBA(I)     = 60.0
        DEF_RADINF(I)      = 3.5
        DEF_OBTHIN(I)      = 1
        DEF_CSCALE_START(I) = 360.0
        DEF_CSCALE_OBTIME(I) = 200.0
        DEF_CSCALE_END(I)  = 200.0
      END DO

      I = I + 1              !  Defaults for Type 901 LOG Visibility
      MASTER_AC_TYPES(I)   = 901
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0

!  Defaults for Limited Area.
#else
      I = I + 1              !  Defaults for Type 101 pstar
      MASTER_AC_TYPES(I)   = 101
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 225.0
        DEF_CSCALE_OBTIME(I) = 150.0
        DEF_CSCALE_END(I)    = 165.0
        DEF_TIMEB(I)         = 120.0
        DEF_TIMEA(I)         = 24.0
        DEF_TGETOBB(I)       = 120.0
        DEF_TGETOBA(I)       = 24.0
        DEF_RADINF(I)        = 1.75
      ENDIF

      I = I + 1              !  Defaults for Type 201 sonde temp
      MASTER_AC_TYPES(I)   = 201
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 120.0
        DEF_CSCALE_OBTIME(I) = 100.0
        DEF_CSCALE_END(I)    = 100.0
      ENDIF

      I = I + 1              !  Defaults for Type 202 ship surf temp
      MASTER_AC_TYPES(I)   = 202
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 225.0
        DEF_CSCALE_OBTIME(I) = 150.0
        DEF_CSCALE_END(I)    = 165.0
      ENDIF

      I = I + 1              !  Defaults for Type 203 airep temp
      MASTER_AC_TYPES(I)   = 203
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 120.0
        DEF_CSCALE_OBTIME(I) = 100.0
        DEF_CSCALE_END(I)    = 100.0
      ENDIF

      I = I + 1              !  Defaults for Type 204 synop surf temp
      MASTER_AC_TYPES(I)   = 204
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 150.0
      DEF_CSCALE_OBTIME(I) = 120.0
      DEF_CSCALE_END(I)    = 120.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 105.0
        DEF_CSCALE_OBTIME(I) = 85.0
        DEF_CSCALE_END(I)    = 85.0
        DEF_TIMEB(I)         = 120.0
        DEF_TIMEA(I)         = 24.0
        DEF_TGETOBB(I)       = 120.0
        DEF_TGETOBA(I)       = 24.0
        DEF_RADINF(I)        = 1.75
      ENDIF

      I = I + 1              !  Defaults for Type 205 LASS temp
      MASTER_AC_TYPES(I)   = 205
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 120.0
        DEF_CSCALE_OBTIME(I) = 100.0
        DEF_CSCALE_END(I)    = 100.0
      ENDIF

      I = I + 1              !  Defaults for Type 206 SATEM temp
      MASTER_AC_TYPES(I)   = 206
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 120.0
        DEF_CSCALE_OBTIME(I) = 100.0
        DEF_CSCALE_END(I)    = 100.0
      ENDIF

      I = I + 1              !  Defaults for Type 207 SAT120 temp
      MASTER_AC_TYPES(I)   = 207
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 120.0
        DEF_CSCALE_OBTIME(I) = 100.0
        DEF_CSCALE_END(I)    = 100.0
      ENDIF

      I = I + 1              !  Defaults for Type 208 constrained LASS
      MASTER_AC_TYPES(I)   = 208
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0

      I = I + 1              !  Defaults for Type 209 BOGUS 1000-500
      MASTER_AC_TYPES(I)   = 209
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 120.0
        DEF_CSCALE_OBTIME(I) = 100.0
        DEF_CSCALE_END(I)    = 100.0
      ENDIF

      I = I + 1              !  Defaults for Type 211 UARS Temp
      MASTER_AC_TYPES(I)   = 211
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0

      I = I + 1              !  Defaults for Type 301 Sonde winds
      MASTER_AC_TYPES(I)   = 301
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 120.0
        DEF_CSCALE_OBTIME(I) = 100.0
        DEF_CSCALE_END(I)    = 100.0
      ENDIF

      I = I + 1              !  Defaults for Type 302 Ship surf wind
      MASTER_AC_TYPES(I)   = 302
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  =  95.0
        DEF_CSCALE_OBTIME(I) =  75.0
        DEF_CSCALE_END(I)    =  75.0
      ENDIF

      I = I + 1              !  Defaults for Type 303 Airep wind
      MASTER_AC_TYPES(I)   = 303
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 120.0
        DEF_CSCALE_OBTIME(I) = 100.0
        DEF_CSCALE_END(I)    = 100.0
      ENDIF

      I = I + 1              !  Defaults for Type 304 Synop surf wind
      MASTER_AC_TYPES(I)   = 304
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 240.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  =  95.0
        DEF_CSCALE_OBTIME(I) =  75.0
        DEF_CSCALE_END(I)    =  75.0
        DEF_TIMEB(I)         = 120.0
        DEF_TIMEA(I)         = 24.0
        DEF_TGETOBB(I)       = 120.0
        DEF_TGETOBA(I)       = 24.0
        DEF_RADINF(I)        = 1.75
      ENDIF

      I = I + 1              !  Defaults for Type 305 Scatwind
      MASTER_AC_TYPES(I)   = 305
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 240.0
      DEF_CSCALE_OBTIME(I) = 120.0
      DEF_CSCALE_END(I)    = 120.0

      I = I + 1              !  Defaults for Type 306 Ship surf wind
      MASTER_AC_TYPES(I)   = 306
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  =  95.0
        DEF_CSCALE_OBTIME(I) =  75.0
        DEF_CSCALE_END(I)    =  75.0
      ENDIF

      I = I + 1              !  Defaults for Type 311 UARS wind
      MASTER_AC_TYPES(I)   = 311
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0

      I = I + 1              !  Defaults for Type 401 Sonde RH
      MASTER_AC_TYPES(I)   = 401
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 115.0
        DEF_CSCALE_OBTIME(I) = 85.0
        DEF_CSCALE_END(I)    = 85.0
      ENDIF

      I = I + 1              !  Defaults for Type 402 Ship surf RH
      MASTER_AC_TYPES(I)   = 402
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 240.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 225.0
        DEF_CSCALE_OBTIME(I) = 150.0
        DEF_CSCALE_END(I)    = 165.0
      ENDIF

      I = I + 1              !  Defaults for Type 403 Bogus RH
      MASTER_AC_TYPES(I)   = 403
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 115.0
        DEF_CSCALE_OBTIME(I) =  85.0
        DEF_CSCALE_END(I)    =  85.0
      ENDIF

      I = I + 1              !  Defaults for Type 404 Synop surf RH
      MASTER_AC_TYPES(I)   = 404
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 150.0
      DEF_CSCALE_OBTIME(I) = 120.0
      DEF_CSCALE_END(I)    = 120.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 100.0
        DEF_CSCALE_OBTIME(I) =  75.0
        DEF_CSCALE_END(I)    =  75.0
        DEF_TIMEB(I)         = 120.0
        DEF_TIMEA(I)         = 24.0
        DEF_TGETOBB(I)       = 120.0
        DEF_TGETOBA(I)       = 24.0
        DEF_RADINF(I)        = 1.75
      ENDIF

      I = I + 1              !  Defaults for Type 405 LASS RH
      MASTER_AC_TYPES(I)   = 405
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0

      I = I + 1              !  Defaults for Type 406 MOPS RH
      MASTER_AC_TYPES(I)   = 406
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 0.01
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 999.0
      DEF_CSCALE_OBTIME(I) = 999.0
      DEF_CSCALE_END(I)    = 999.0

      I = I + 1              !  Defaults for Type 407 Cloud histograms
      MASTER_AC_TYPES(I)   = 407
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 50.0
      DEF_CSCALE_OBTIME(I) = 50.0
      DEF_CSCALE_END(I)    = 50.0

      I = I + 1              !  Defaults for type 506 MOPS precip
      MASTER_AC_TYPES(I)   = 506
      DEF_TIMEB(I)         = 180.
      DEF_TIMEA(I)         = 180.
      DEF_TGETOBB(I)       = 180.
      DEF_TGETOBA(I)       = 180.
      DEF_RADINF(I)        = 0.01
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 999.0
      DEF_CSCALE_OBTIME(I) = 999.0
      DEF_CSCALE_END(I)    = 999.0

      I = I + 1              !  Defaults for Type 901 LOG Visibility
      MASTER_AC_TYPES(I)   = 901
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 100.0
        DEF_CSCALE_OBTIME(I) =  75.0
        DEF_CSCALE_END(I)    =  75.0
        DEF_TIMEB(I)         = 120.0
        DEF_TIMEA(I)         = 24.0
        DEF_TGETOBB(I)       = 120.0
        DEF_TGETOBA(I)       = 24.0
        DEF_RADINF(I)        = 1.75
      ENDIF

      !  Defaults for Type 6nn (tracers)
      !  Allow for A_MAX_TRVARS (currently 29) possible tracers
      DO J=1,A_MAX_TRVARS
        I = I + 1
        MASTER_AC_TYPES(I) = 600+J
        DEF_TIMEB(I)       = 150.0
        DEF_TIMEA(I)       = 30.0
        DEF_TGETOBB(I)     = 150.0
        DEF_TGETOBA(I)     = 30.0
        DEF_RADINF(I)      = 3.5
        DEF_OBTHIN(I)      = 1
        DEF_CSCALE_START(I) = 300.0
        DEF_CSCALE_OBTIME(I) = 180.0
        DEF_CSCALE_END(I)  = 180.0
      END DO

#endif

!     Modify defaults for UARS assimilation.
      IF (LAC_UARS) THEN
        DO JOBT=1,I
          DEF_CSCALE_START(JOBT)  = 600.0
          DEF_CSCALE_OBTIME(JOBT) = 400.0
          DEF_CSCALE_END(JOBT)    = 400.0
          IF (MASTER_AC_TYPES(JOBT) == 207) THEN
            DEF_OBTHIN(JOBT)      = 2
          ELSE
            DEF_OBTHIN(JOBT)      = 1
          END IF
        ENDDO
      ENDIF

!     Initialise rest of arrays.
      IF (I <  NOBTYPMX) THEN
        DO JOBT = I+1,NOBTYPMX
          MASTER_AC_TYPES(JOBT)   = 0
          DEF_TIMEB(JOBT)         = 0.0
          DEF_TIMEA(JOBT)         = 0.0
          DEF_TGETOBB(JOBT)       = 0.0
          DEF_TGETOBA(JOBT)       = 0.0
          DEF_RADINF(JOBT)        = 0.0
          DEF_OBTHIN(JOBT)        = 0
          DEF_CSCALE_START(JOBT)  = 0.0
          DEF_CSCALE_OBTIME(JOBT) = 0.0
          DEF_CSCALE_END(JOBT)    = 0.0
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE DEF_TYPE
#endif
