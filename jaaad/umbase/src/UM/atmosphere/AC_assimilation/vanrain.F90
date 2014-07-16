#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE VANRAIN-------------------------------------------------
!LL
!LL  Purpose : Performs 'vertical' analysis for precip rate/phase data.
!LL
!LL     When IPASS=1,
!LL     model field is interpolated to ob locations and increments
!LL     calculated. Preliminary weight normalisation factors are also
!LL     calculated.
!LL
!LL     When IPASS=2,
!LL     data density is interpolated to ob locations and
!LL     final weight normalisation factors are calculated.
!LL
!LL  For use on Cray
!LL  For Cray - Global  ; Enable defs GLOBAL
!LL
!LL  Model            Modification history:
!LL version  date
!LL  3.3    Nov 1993  : adapted from VANPST by C D Jones
!LL  3.4    Aug 1994  : Only use incrs within given range of radars,
!LL                   : with an option to flag obs outside a given
!LL                   : lat/lon range in verification.  (C D Jones)
!LL  4.0  3/5/95 : Calculate weight normalisation factors as a
!LL              :   function of range from nearest radar (Chris Jones)
!LL   4.1   May 1996 : Remove any use of precipitation phase, as
!LL                  : Nimrod data does not include this. (C D Jones)
!LL   4.2   29/11/96  T3e mods Stuart Bell
!LL   5.2   12/12/00  change A to Earth_radius
!LL                   and remove daco-related code   B Macpherson
!     5.3   12/07/01  remove code to weight obs by distance
!                     to radar & lat/lon box diags
!                                          B Macpherson
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LLEND------------------------------------------------------------------
!
!*L  Arguments:---------------------------------------------------
       SUBROUTINE VANRAIN (KACT,IPASS,LS_RAIN,LS_SNOW,                  &
     &                     CONV_RAIN,CONV_SNOW,LENFLD,                  &
     &                     OBDATA,CF1PT,CF2PT,CF3PT,CF4PT,              &
     &                     NP1PT,NP2PT,NP3PT,NP4PT,                     &
     &                     PRINC,NORMF,OBS_LAT,OBS_LONG,                &
     &                     OBS_NO,LMISSD,NPTOBT,                        &
     &                     P_LEVELS,LENOBT,NDV,LENOB,                   &
     &                     NO_ANAL_LEVS,NO_ANAL_VAR,                    &
     &                     ICODE,CMESSAGE)
!
      IMPLICIT NONE
!-----------------------------------------------------------------
!  Analysis Correction comdecks
!
#include "acparm.h"
#include "comacp.h"
#include "comobs.h"
#include "comacdg.h"
#if !defined(GLOBAL)
#include "commg.h"
#endif
#include "c_a.h"
#include "c_pi.h"
!
! ----------------------------------------------------------------------
      INTEGER LENFLD,LENOBT,LENOB,NDV
      INTEGER P_LEVELS
      INTEGER KACT,IPASS,NPTOBT,NO_ANAL_LEVS,NO_ANAL_VAR
      REAL                                                              &
     &               LS_RAIN(LENFLD),LS_SNOW(LENFLD),                   &
     &               CONV_RAIN(LENFLD),CONV_SNOW(LENFLD),               &
     &               OBDATA(LENOBT,NDV),                                &
     &               PRINC(LENOB),NORMF(LENOB),                         &
     &               OBS_LAT(LENOBT),OBS_LONG(LENOBT),                  &
     &               CF1PT(LENOBT),CF2PT(LENOBT),CF3PT(LENOBT),         &
     &               CF4PT(LENOBT)
!
      LOGICAL LMISSD(LENOB+1,NO_ANAL_LEVS)
!
      INTEGER OBS_NO(LENOBT)
      INTEGER NP1PT(LENOBT),NP2PT(LENOBT),NP3PT(LENOBT)
      INTEGER NP4PT(LENOBT)
!
      INTEGER ICODE
      CHARACTER*256 CMESSAGE
!
!-INTENT=IN--------------------------------------------------------
!     KACT      - Observation type
!     IPASS     - Calculate incrs and preliminary norm factors for
!                 IPASS=1 and final norm factors for IPASS=2.
!     LENOBT    - no of obs in type
!     LENOB     - no of obs in group
!     NPTOBT    - pointer to present ob type within group
!     NDV       - no of data values in observation
!     NO_ANAL_LEVS - no of analysis levels
!     P_LEVELS  - no of model level
!     LS_RAIN   - }  model precipitation rates (IPASS=1) for
!     LS_SNOW   - }  large-scale and convective rain and snow
!     CONV_RAIN - }  or data density on model grid (IPASS=2)
!     CONV_SNOW - }
!     OBDATA    - observed values
!     OBS_LAT   - ob co-lats
!     OBS_LONG  - ob longitudes
!     OBS_NO    - observation numbers
!     NP1PT,NP2PT,NP3PT,NP4PT,
!               - pointers to nearest model grid points for each ob
!     CF1PT,CF2PT,CF3PT,CF4PT,
!               - weighting of model points in bilinear interpolation
!
!-INTENT=INOUT-----------------------------------------------------
!     PRINC        -  ob-model increments
!     NORMF        -  normalisation factors
!-INTENT=OUT-----------------------------------------------------
!     LMISSD          - logical to indicate missing data
!     ICODE,CMESSAGE  - error code and message
!*
!---------------------------------------------------------------------
!*L   Workspace usage
!-----------------------------------------------------------------------
!     DYNAMIC ALLOCATION ON CRAY
      REAL           INTPR(LENOBT),TOTPR(LENFLD)
      REAL           PHSINT(LENOBT),PHSFLD(LENFLD)
      REAL           WLAT(LENOBT),WLON(LENOBT)
!     INTPR   - model rate field interpolated to ob locations
!     TOTPR   - total rate field, =sum of LS/CONV_RAIN/SNOW
!     WLAT    - Equatorial latitude of obs in degrees
!     WLON    - Equatorial longitude of obs in degrees
!*
!----------------------------------------------------------------------
!*L   External subroutine calls
!-----------------------------------------------------------------------
      EXTERNAL TIMER,DIAGOPR,HINTMO
#if !defined(GLOBAL)
      EXTERNAL EQTOLL
#endif
!
!----------------------------------------------------------------------
!     Define local variables
!-----------------------------------------------------------------------
      INTEGER KTYPE,NLEV,NOBIPT,NEROPT,JOB,JPT
      INTEGER JRAD
      REAL  OB_TO_RADAR, OB_TO_RADAR_MIN(LENOBT)
      REAL  F2
!     KTYPE    -   Observation type
!     NLEV     -   No of observed levels
!     NOBIPT   -   Pointer to data values within OBDATA array
!     NEROPT   -   Pointer to obs errors within OBDATA array
!     JOB      -   Loop counter in loops over obs
!     JPT      -   Loop counter in loops over model points
!     JRAD     -   Loop counter over radars
!     OB_TO_RADAR      - Distance (km) calculated from ob to radar
!     OB_TO_RADAR_MIN  - Min value of OB_TO_RADAR
!----------------------------------------------------------------------
!
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('VANRAIN ',3)
!
      KTYPE = LACT(KACT)

      NLEV  = NOBLEV(KACT)
      IF (KTYPE == 506) THEN
        NOBIPT = 1
      ENDIF
      NEROPT = NERLEV1(KACT)-NDVHDR
!
!   COMBINE model precipitation rates to give TOTPR
!

      IF (IPASS  ==  1) THEN
         DO 1000 JPT=1,LENFLD
            TOTPR(JPT)=  LS_RAIN(JPT)   + LS_SNOW(JPT)                  &
     &                 + CONV_RAIN(JPT) + CONV_SNOW(JPT)
1000     CONTINUE
      ELSEIF (IPASS  ==  2) THEN   ! TOTPR set to obs density
         DO 1010 JPT=1,LENFLD
            TOTPR(JPT) = LS_RAIN(JPT)
1010     CONTINUE
      ENDIF
!
!
!L*** 1.   OBSERVATIONS GIVING DATA FOR PRECIP RATE/PHASE
!L         ----------------------------------------------
!
!L*** 1.1  HORIZONTAL INTERP OF RATE FIELD ON MODEL GRID
!
!
! DEPENDS ON: hintmo
        CALL HINTMO(TOTPR,CF1PT,CF2PT,CF3PT,CF4PT,                      &
     &              NP1PT,NP2PT,NP3PT,NP4PT,                            &
     &              LENFLD,1,LENOBT,INTPR,                              &
     &              ICODE,CMESSAGE)
        IF (ICODE >  0) GO TO 999

!
!L*** 1.2  VERTICAL INTERP/PROCESSING OF MODEL
!          NONE.
!L*** 1.3  SUBTRACTION, CALCULATION OF PRELIM NORM FACTOR (IPASS=1)
!L                      CALCULATION OF FINAL NORM FACTOR (IPASS=2)
!
      IF (IPASS == 1) THEN
!       PR Increment = PR Observation - Background PR
        DO 1300 JOB=1,LENOBT
          PRINC(NPTOBT+JOB) = OBDATA(JOB,NOBIPT) - INTPR(JOB)
1300    CONTINUE
      ENDIF

!
!     SET WEIGHT FACTOR TO ZERO WHERE OB. VALUE OR ERROR IS FLAGGED
      DO 1310 JOB=1,LENOBT
        LMISSD(NPTOBT+JOB,1) = OBDATA(JOB,NOBIPT) == MISSD .OR.         &
     &                         OBDATA(JOB,NEROPT) == MISSD
1310  CONTINUE
!
      IF (IPASS == 1) THEN
        IF (MDATADFN == 1) THEN
          DO 1320 JOB=1,LENOBT
            IF (LMISSD(NPTOBT+JOB,1)) THEN
              NORMF(NPTOBT+JOB) = 0.0
            ELSE
              NORMF(NPTOBT+JOB) = 1.0 /                                 &
     &                  SQRT(OBDATA(JOB,NEROPT)*OBDATA(JOB,NEROPT)+1.0)
            ENDIF
1320      CONTINUE
        ELSEIF (MDATADFN == 2) THEN
          DO 1330 JOB=1,LENOBT
            IF (LMISSD(NPTOBT+JOB,1)) THEN
              NORMF(NPTOBT+JOB) = 0.0
            ELSE
              NORMF(NPTOBT+JOB) = 1.0 /                                 &
     &                        ( OBDATA(JOB,NEROPT)*OBDATA(JOB,NEROPT) )
            ENDIF
1330      CONTINUE
        ENDIF
!
!
      ELSEIF (IPASS == 2) THEN
        IF (MDATADFN == 1) THEN
          DO 1340 JOB=1,LENOBT
            IF (LMISSD(NPTOBT+JOB,1)) THEN
              NORMF(NPTOBT+JOB) = 0.0
            ELSE
              NORMF(NPTOBT*JOB) = NORMF(NPTOBT+JOB)*NORMF(NPTOBT+JOB) / &
     &         ( INTPR(JOB)*NORMF(NPTOBT+JOB) +1.0 -                    &
     &           NORMF(NPTOBT+JOB)*NORMF(NPTOBT+JOB) )
            ENDIF
1340      CONTINUE
        ELSEIF (MDATADFN == 2) THEN
          DO 1350 JOB=1,LENOBT
            IF (LMISSD(NPTOBT+JOB,1)) THEN
              NORMF(NPTOBT+JOB) = 0.0
            ELSE
              NORMF(NPTOBT+JOB) = NORMF(NPTOBT+JOB)/( INTPR(JOB)+1.0 )
            ENDIF
1350      CONTINUE
        ENDIF
!
      ENDIF
!
!L*** 1.4  VERTICAL INTERP/PROCESSING OF INCS. CALC'N OF WEIGHT FACTORS.
!          NONE.
!
!     DIAGNOSTICS ON FIT OF MODEL PRECIP RATE/PHASE TO OBS
      IF (LDIAGAC .AND. (IPASS == 1) ) THEN

! DEPENDS ON: diagopr
        CALL DIAGOPR (lenobt,INTPR,OBDATA,LENOBT,LENOB,NDV,             &
     &                LMISSD,NOBIPT,NPTOBT,NO_ANAL_LEVS)
      ENDIF

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('VANRAIN ',4)

 999  CONTINUE
      RETURN
      END SUBROUTINE VANRAIN
#endif
