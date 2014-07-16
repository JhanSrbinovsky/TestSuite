#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE DIAGO --------------------------------------------------
!LL
!LL  Purpose : Provide Statistics on 'Observation-Model' Increments
!LL
!LL            Calculate and print out the following :
!LL            Mean Increment
!LL            Mean Square Increment
!LL            Extreme Increment and position (Lat/Long) of obs
!LL
!LL            Can be used on Individual or group of AC Obs Types
!LL
!LL  For Cray - Global  ; Enable defs GLOBAL
!LL
!LL  Written by Dave Robinson.
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  3.2     10/6/93  For LLBAND option areas made consistent with CF
!LL                   verif areas; cater for type 306;
!LL                   and Eliminate QA FORTRAN complaints    S Bell
!LL  3.3     13/9/93  Remove WRITE of error message          S Bell
!LL  3.3     10/12/93  Extend use of KCALL to handle call from DIAGOCC
!LL                                                       B Macpherson
!LL  4.1   31/05/96     The number of v points to be processed on a
!LL                     C grid differs from u by row_length. u,v
!LL                     dimensioned separately in call to WLLTOEQ.
!LL                     Requirement for VAR.
!LL                     Author I.Edmond       Reviewer D. Goddard
!    4.2 25/11/96: T3E mods. Stuart Bell
!LL
!LL  4.1  4/09/96:  Port to CRAY T3E  Deborah Salmond
!LL  5.2  30/11/00:  remove use of ak,bk for new dynamics
!LL                                        B Macpherson
!LL  6.0  17/10/03  Replace SHMEM with GCOM for port to SX6. Clive Jones
!LL  6.2  21/10/05  Replace GSYNC with SSYNC. P.Selwood
!LL  6.2  15/08/05  Free format fixes. P. Selwood
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Logical components covered:
!LL
!LL  Project Task : P3
!LL
!LL  Documentation:
!LL
!LLEND------------------------------------------------------------------
!*L  Arguments----------------------------------------------------------
      SUBROUTINE DIAGO (CSUBR,KTYPE,KCALL,                              &
     &           OBS_INCR,NORMF,OBS_LAT,OBS_LONG,                       &
     &           LMISSD,LENOB,LENOBT,NPTOBT,                            &
     &           NO_ANAL_LEVS,NO_ANAL_VAR)

!     IMPLICIT NONE  !FAILS AUTOTASKING (UNDECLARED/UNTYPED EXTERNALS)

      CHARACTER *(*) CSUBR   !  Subroutine where DIAGO is called from
      INTEGER KTYPE,KCALL,LENOB,LENOBT,NPTOBT,NO_ANAL_LEVS
      INTEGER NO_ANAL_VAR
      REAL                                                              &
     & OBS_INCR(LENOB+1,NO_ANAL_LEVS,NO_ANAL_VAR),                      &
     & NORMF (LENOB+1,NO_ANAL_LEVS), OBS_LAT (LENOBT+1),                &
     & OBS_LONG (LENOBT+1)
      LOGICAL LMISSD(LENOB+1,NO_ANAL_LEVS)

!     INTENT=IN---------------------------------------------------------
!     KTYPE    : AC Observation Type
!     KCALL    : Indicator to where DIAGO called from
!              : 1 = FROM AC
!                2 = FROM VAN## ; IPASS=1 ; BEFORE VERTICAL FILTERING
!                3 = FROM VAN## ; IPASS=1 ; AFTER  VERTICAL FILTERING
!                4 = FROM VAN## ; IPASS=2 ; BEFORE VERTICAL FILTERING
!                5 = FROM VAN## ; IPASS=2 ; AFTER  VERTICAL FILTERING
!                6 = from DIAGOCC ; IPASS=1 ; for multi-level cloud
!                7 = from DIAGOCC ; IPASS=1 ; for low cloud
!                8 = from DIAGOCC ; IPASS=1 ; for medium cloud
!                9 = from DIAGOCC ; IPASS=1 ; for high cloud
!               10 = from DIAGOCC ; IPASS=1 ; for total cloud
!                VAN## is Vertical Analysis routine calling DIAGO.
!     INTENT=INOUT------------------------------------------------------
!     INC      : Increments - P*, T, U, RH
!     VINC     : Increments - V
!     NORMF    : Normalisation factor
!     OBS_LAT  : Observation co-latitudes
!     OBS_LONG : Observation longitudes
!     LMISSD   : Array to control which obs are used on each level
!     LENOB    : No of obs in group of observation types
!     LENOBT   : No of obs for Observation Type KTYPE
!     NPTOBT   : Offset to first observation for type KTYPE
!     NO_ANAL_LEVS : No of analysis levels
!     INTENT=OUT--------------------------------------------------------
!*   ------------------------------------------------------------------

!-----AC common blocks
#include "acparm.h"
#if defined(MPP)
#include "parvars.h"
#include "mppac.h"
      real a_maxinc(P_LEVELS_MAX)
      integer iproc,istat
#endif
#include "comacp.h"
#include "comacdg.h"
#include "docacdg.h"
#if !defined(GLOBAL)
#include "commg.h"
#endif
#include "c_pi.h"
#include "c_r_cp.h"

!*L   Workspace Usage:--------------------------------------------------
!     -----------------------------------------------------------
!     Local work arrays
!     -----------------------------------------------------------
!**   Dynamic allocation
      LOGICAL OBSUSE(LENOBT+1)
#if !defined(GLOBAL)
      REAL WLT(LENOBT+1), WLN(LENOBT+1), WLN2(LENOBT+1)
      REAL COEFF1(LENOBT+1), COEFF2(LENOBT+1)
      REAL UWK(LENOBT+1,NO_ANAL_LEVS), VWK(LENOBT+1,NO_ANAL_LEVS)
#endif
#if !defined(GLOBAL)
!     WLT/WLN  : Real lat/long of wind obs
!     WLN2     : ELF  longitude of wind obs
!     COEFF1/2 : Coefficients to rotate wind components
!     UWK/VWK  : Rotated u/v components on real lat/long grid
#endif
!*    ------------------------------------------------------------------

!*L   External Subroutine Calls:----------------------------------------
      EXTERNAL TIMER
#if !defined(GLOBAL)
      EXTERNAL EQTOLL,W_COEFF,W_EQTOLL
#endif
!*    ------------------------------------------------------------------

!     ------------------------------------------------------------------
!     Local arrays and variables

      INTEGER JVAR,JBAND,JLEV,JOB,JOBT,JNANL
      INTEGER NBAND,INOBS,INOBSLV,IOBT,ANAL_VAR
      REAL STAT(0:NO_ANAL_LEVS,8),LAT,LONG,UI,VI,NF,SQINC
      REAL BANDLAT(4),MAXINC,SUMUI,SUMVI,SUMNF,SUMSQI
      CHARACTER NS*1,WE*1,TITLE1*9,TITLE2*25,LAB1*10,LAB2*10
      CHARACTER*3 BANDC(4)
      LOGICAL LWIND

!     ANAL_VAR : Analysis varibale (=KTYPE/100)
!     BANDLAT  : Boundaries of latitude bands
!     BANDC    : Labels to print out boundaries of latitude bands
!     INOBS    : Number of observations in latitude band
!     INOBSLV  : Number of observations on level
!     IOBT     : Index to list of obs types in MDIAGTPS
!     JBAND    : Loop variable for latitude bands
!     JLEV     : Loop variable for analysis levels
!     JNANL    : Number of analysis levels to be processed
!     JOB      : Loop variable for observations
!     JOBT     : Loop variable for observation types
!     JVAR     : Loop variable for statistic variable
!     LAB1/2   : Labels for print out
!     LAT      : Absolute value of latitude of maximum increment
!     LONG     : Absolute value of longitude of maximum increment
!     LWIND    : Indicator for working with wind obs type
!     MAXINC   : Maximum increment for level
!     NBAND    : Number of latitude bands
!     NF       : Normalisation Factors
!     NS/WE    : Labels for printing out lat/long of maximum increment
!     SQINC    : Square of increments
!     STAT     : Array storing statistics data
!                STAT(0,#) Overall
!                STAT(JLEV,#) For level JLEV
!                STAT(#,1) No of observations on  level JLEV
!                STAT(#,2) Mean Norm factor
!                STAT(#,3) Mean Increment - p*, t, u, rh
!                STAT(#,4) Mean Increment - v
!                STAT(#,5) Mean square Increment
!                STAT(#,6) Extreme increment
!                STAT(#,7) Lat of extreme increment
!                STAT(#,8) Long of extreme increment
!     SUMNF    : Sum of norm. factors
!     SUMUI    : Sum of increments - p*, t, u, rh
!     SUMVI    : Sum of increments - v
!     SUMSQI   : Sum of squared increments
!     TITLE1/2 : Titles for print out
!     UI       : Increments - p*, t, u, rh
!     VI       : Increments - v

!     ------------------------------------------------------------------
      DATA BANDLAT / 0.0, 1.0472, 2.0944, 3.1416 /
      DATA BANDC/'90N','30N','30S','90S'/
!     -----------------------------------------------------------
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('DIAGO   ',3)

!     Check that 'Obs-Model' statistics required
#if defined(MPP)
      IF (.NOT.LLDAC(2)) GO TO 999
#else
      IF (.NOT.LLDAC(2) .OR. LENOBT <= 0) GO TO 999
#endif

!     For verification purposes, only get statistics on first pass.
      IF (LVERIF.AND.(KCALL /= 2.AND.KCALL /= 3)) GO TO 999

!     Check that statistics required for calling routine
      IF (MDIAG >  0 .AND. KCALL /= MDIAG) GO TO 999

!     Check that statistics required for observation type KTYPE
      IF (MDIAGTPS(1) >  0) THEN
        IOBT=0
        DO 5 JOBT=1,NOBTYPMX
        IF (KTYPE == MDIAGTPS(JOBT)) IOBT=JOBT
   5    CONTINUE
        IF (IOBT == 0) GO TO 999
      ENDIF

!---- Set up titles

      LAB1='    '
      LAB2='    '
      ANAL_VAR = KTYPE/100
      IF (ANAL_VAR == 1) THEN
        LAB1 = 'P* OBS INC'
      ELSEIF (ANAL_VAR == 2) THEN
        IF (LTEMP) THEN
          LAB1 = ' T OBS INC'
        ELSE
          LAB1 = 'TH OBS INC'
        ENDIF
      ELSEIF (ANAL_VAR == 3) THEN
        LAB1 = ' U OBS INC'
        LAB2 = ' V OBS INC'
      ELSEIF (ANAL_VAR == 4) THEN
        IF ( KCALL >  5 ) THEN
          LAB1 = 'CC OBS INC'
        ELSE
          LAB1 = 'RH OBS INC'
        ENDIF
      ELSEIF (ANAL_VAR == 9) THEN
        LAB1 = 'LOG(VIS) OBS INC'
      ELSE
      GOTO 999
      ENDIF

      IF (KCALL == 2 .OR. KCALL == 3) THEN
        TITLE1 = 'IPASS=1'
      ELSEIF (KCALL == 4 .OR. KCALL == 5) THEN
        TITLE1 = 'IPASS=2'
      ELSE
        TITLE1 = '  '
      ENDIF

      TITLE2 = ' '
      IF (KTYPE == 201 .OR. KTYPE == 205 .OR. KTYPE == 206 .OR.         &
     &    KTYPE == 207 .OR. KTYPE == 208 .OR. KTYPE == 209 .OR.         &
     &    KTYPE == 301 .OR. KTYPE == 401 ) THEN
        IF (KCALL == 2 .OR. KCALL == 4) THEN
          TITLE2 = 'BEFORE VERTICAL FILTERING'
        ELSEIF (KCALL == 3 .OR. KCALL == 5) THEN
          TITLE2 = 'AFTER VERTICAL FILTERING'
        ELSE
          TITLE2 = ' '
        ENDIF
      ENDIF

      LWIND = ANAL_VAR == 3

      JNANL = NO_ANAL_LEVS

      IF (KTYPE == 202 .OR. KTYPE == 204 .OR.                           &
     &    KTYPE == 302 .OR. KTYPE == 304 .OR.                           &
     &    KTYPE == 305 .OR. KTYPE == 306 .OR.                           &
     &    KTYPE == 402 .OR. KTYPE == 404                                &
     &    .OR. KTYPE  ==  901) THEN
        JNANL=1
      ENDIF

#if defined(GLOBAL) && !defined(MPP)
      IF (LLBAND) THEN
        NBAND=3
      ELSE
        NBAND=1
      ENDIF
#else
!     ONLY ONE BAND OF STATISTICS ALLOWED FOR ELF GRID
      NBAND=1
#endif

!     Loop over latitude bands
      DO 100 JBAND=1,NBAND

!     Initialise statistics data
      DO JVAR=1,8
        DO JLEV=0,NO_ANAL_LEVS
          STAT(JLEV,JVAR)=0.0
        ENDDO
      ENDDO

!L    1.0 Determine which observations to be used
#if defined(MPP)
      IF (LENOBT <= 0) THEN
       do jlev=1,jnanl
        do jvar=0,8
         s_stat(jlev,jvar)=0.0
        enddo
       enddo
      ELSE
#endif

!     1.1 Initialise so that all observations to be used
      DO JOB=1,LENOBT
        OBSUSE(JOB)=.TRUE.
      ENDDO

#if defined(GLOBAL)
!     1.2 Skip observations outside latitude band
      IF (LLBAND) THEN
        DO JOB=1,LENOBT
          OBSUSE(JOB) = OBS_LAT(JOB) >= BANDLAT(JBAND)   .AND.          &
     &                  OBS_LAT(JOB) <  BANDLAT(JBAND+1)
        ENDDO
      ENDIF
#endif

!     1.3 Count number of observations to be used
      INOBS = 0
      DO JOB=1,LENOBT
        IF (OBSUSE(JOB)) INOBS = INOBS+1
      ENDDO

#if !defined(MPP)
      IF (INOBS >  0) THEN
#endif

#if !defined(GLOBAL)

!       2. Co-ordinate tranformation  (ELF version only)
!       ------------------------------------------------

!       Get real Lat/Long of ELF observations

!       Lat/Lon values in OBS_LAT/OBS_LONG are for ELF grid.
!       First convert from radians to degrees.
        DO JOB=1,LENOBT
          WLT (JOB) = 90.0-OBS_LAT(JOB)*RECIP_PI_OVER_180
          WLN (JOB) = OBS_LONG(JOB)*RECIP_PI_OVER_180
          WLN2(JOB) = WLN(JOB)
        ENDDO
!       Get real Lat/Lon of wind obs in WLT/WLN, keep ELF lon in WLN2
! DEPENDS ON: eqtoll
        CALL EQTOLL (WLT,WLN,WLT,WLN,ELFPLAT,ELFPLON,LENOBT)

        IF (LWIND) THEN

!         Rotate wind components on ELF grid to real lat/lon system.

!         Calculate coefficients of rotation.
! DEPENDS ON: w_coeff
          CALL W_COEFF(COEFF1,COEFF2,WLN,WLN2,ELFPLAT,ELFPLON,LENOBT)
!         Rotation of wind components in loop over levels.
!FPP$ CNCALL
          DO JLEV=1,JNANL
! DEPENDS ON: w_eqtoll
            CALL W_EQTOLL (COEFF1,COEFF2,                               &
     &                     OBS_INCR(NPTOBT+1,JLEV,1),                   &
     &                     OBS_INCR(NPTOBT+1,JLEV,2),                   &
     &                     UWK(1,JLEV),VWK(1,JLEV),LENOBT,LENOBT)
          ENDDO

        ENDIF
#endif

!       Autotasking Directives
!       ----------------------
!       Note that all CMIC$ calls have been commented out as it
!       causes problems when autotasking. Problem has not yet been
!       investigated.

!CMIC$  DO ALL AUTOSCOPE
!CMIC$1 PRIVATE(WU,WU2,WV,WNF,JLEV,JOB,INOBSLV,SUM,MAXINC)
!CMIC$2 SHARED(JNANL,LENOBT,OBSUSE,LMISSD,NPTOBT)
!CMIC$3 SHARED(NORMF,INC,VINC,STAT,LWIND,ANAL_VAR)
#if !defined(GLOBAL)
!CMIC$4 SHARED(UWK,VWK)
#endif

!       3.0 Calculate statistics required
!       ---------------------------------

!-----  Loop over levels
        DO 300 JLEV=1,JNANL

        INOBSLV=0
        SUMNF  = 0.0
        SUMUI  = 0.0
        SUMVI  = 0.0
        SUMSQI = 0.0
        MAXINC = 0.0

!-----  Loop over observations
        DO 305 JOB=1,LENOBT

        NF    = 0.0
        UI    = 0.0
        VI    = 0.0
        SQINC = 0.0

        IF (OBSUSE(JOB) .AND. .NOT.LMISSD(NPTOBT+JOB,JLEV)) THEN

          INOBSLV = INOBSLV+1

          IF (LNORMF) THEN
!           Use NF as weighting
            NF = NORMF (NPTOBT+JOB,JLEV)
          ELSE
!           Do not use NF as weighting (set to 1.0)
            NF = 1.0
          ENDIF
#if defined(GLOBAL)
          UI = OBS_INCR(NPTOBT+JOB,JLEV,1)
          IF (LWIND) VI = OBS_INCR(NPTOBT+JOB,JLEV,2)
#else
          IF (LWIND) THEN
            UI = UWK(JOB,JLEV)
            VI = VWK(JOB,JLEV)
          ELSE
            UI = OBS_INCR(NPTOBT+JOB,JLEV,1)
          ENDIF
#endif


!-----    3.1 Sum of Normalistion factors used as weights

!         For IPASS=1, the preliminary norm factor (Q1M) calculated in
!                      VAN### is used to weight the obs increments.
!         For IPASS=2, the final norm factor (Q2M) calculated in
!                      VAN### is used to weight the obs increments.

          SUMNF = SUMNF + NF

!-----    3.2 Weighted sum of increments - P*, T, U, RH
          SUMUI = SUMUI + NF*UI

!-----    3.3 Weighted sum of increments - V
          IF (LWIND) SUMVI = SUMVI + NF*VI

!-----    3.4 Weighted sum of squared increments
          SQINC = UI*UI
          IF (LWIND) SQINC = SQINC + VI*VI
          SUMSQI = SUMSQI + NF*SQINC

!-----    3.5 Determine maximum (squared) increment
          IF (SQINC  >   MAXINC) THEN
            MAXINC = SQINC
            IF (LWIND) THEN
!             For wind, extreme is maximum speed of vector increment
              STAT(JLEV,6) = SQRT(SQINC)
            ELSE
!             Extreme is increment
              STAT(JLEV,6) = UI
            ENDIF
!           Lat/Long of obs
#if defined(GLOBAL)
            STAT(JLEV,7) = OBS_LAT(JOB)
            STAT(JLEV,8) = OBS_LONG(JOB)
!-----      Convert to degrees
            STAT(JLEV,7) = 90.0 - STAT(JLEV,7) * RECIP_PI_OVER_180
            STAT(JLEV,8) = STAT(JLEV,8) * RECIP_PI_OVER_180
            IF (STAT(JLEV,8) >  180.0) STAT(JLEV,8)=STAT(JLEV,8)-360.0
#else
            STAT(JLEV,7) = WLT(JOB)
            STAT(JLEV,8) = WLN(JOB)
!-----      Already in degrees
#endif
          ENDIF

        ENDIF

  305   CONTINUE
#if defined(MPP)
        s_stat(jlev,0)=maxinc
        s_stat(jlev,1)=real(inobslv)
        s_stat(jlev,2)=sumnf
        s_stat(jlev,3)=sumui
        if(lwind)s_stat(jlev,4)=sumvi
        s_stat(jlev,5)=sumsqi
        s_stat(jlev,6)=STAT(JLEV,6)
        s_stat(jlev,7)=STAT(JLEV,7)
        s_stat(jlev,8)=STAT(JLEV,8)
  300   CONTINUE
      ENDIF   ! LENOBT <= 0
      CALL GC_SSYNC(NPROC,ISTAT)
      if(mype == 0)then
!
!   set stat to zero
!
        do jlev=1,jnanl
          stat(jlev,1)=0.0
          stat(jlev,2)=0.0
          stat(jlev,3)=0.0
          stat(jlev,4)=0.0
          stat(jlev,5)=0.0
          a_maxinc(jlev)=0.0
        enddo
      endif  ! mype=0

! Gather stats on PE0
      if(mype == 0) then ! stats off PE0
        do jlev=1,jnanl
          stat(jlev,1)=s_stat(jlev,1)
          stat(jlev,2)=s_stat(jlev,2)
          stat(jlev,3)=s_stat(jlev,3)
          if(lwind)stat(jlev,4)=s_stat(jlev,4)
          stat(jlev,5)=s_stat(jlev,5)
          if(s_stat(jlev,0) >  a_maxinc(jlev))then
            a_maxinc(jlev)=s_stat(jlev,0)
            stat(jlev,6)=s_stat(jlev,6)
            stat(jlev,7)=s_stat(jlev,7)
            stat(jlev,8)=s_stat(jlev,8)
          endif
        enddo
      endif

      CALL GC_SSYNC(NPROC,ISTAT)

      DO IPROC=1,NPROC-1
        IF(mype == 0) THEN  ! PE0 receives
          CALL GC_RRECV(IPROC,9*p_levels_max,IPROC,ISTAT,r_stat,s_stat)
          do jlev=1,jnanl
            stat(jlev,1)=stat(jlev,1)+r_stat(jlev,1)
            stat(jlev,2)=stat(jlev,2)+r_stat(jlev,2)
            stat(jlev,3)=stat(jlev,3)+r_stat(jlev,3)
            if(lwind)stat(jlev,4)=stat(jlev,4)+r_stat(jlev,4)
            stat(jlev,5)=stat(jlev,5)+r_stat(jlev,5)
            if(r_stat(jlev,0) >  a_maxinc(jlev))then
              a_maxinc(jlev)=r_stat(jlev,0)
              stat(jlev,6)=r_stat(jlev,6)
              stat(jlev,7)=r_stat(jlev,7)
              stat(jlev,8)=r_stat(jlev,8)
            endif
          enddo
        ELSEIF(mype == IPROC) THEN  ! other PEs send
          CALL GC_RSEND(IPROC,9*p_levels_max,0,ISTAT,r_stat,s_stat)
        ENDIF
      ENDDO

      CALL GC_SSYNC(NPROC,ISTAT)


      if(mype == 0) then
#endif

!-----  3.6 Get means for level JLEV and accumulate for overall means.
#if defined(MPP)
      do jlev=1,jnanl
      inobslv=int(stat(jlev,1))
      sumnf=stat(jlev,2)
      sumui=stat(jlev,3)
      if(lwind)sumvi=stat(jlev,4)
      sumsqi=stat(jlev,5)
#endif

        IF (INOBSLV >  0) THEN

          STAT(JLEV,1) = REAL( INOBSLV )
          STAT(0,1)    = STAT(0,1)+STAT(JLEV,1)
          IF (SUMNF >  0.0) THEN
!           Mean weight
            STAT(JLEV,2) = SUMNF/STAT(JLEV,1)
            STAT(0,2)    = STAT(0,2)+SUMNF
!           Mean observation increment - P*, T, U, RH
            STAT(JLEV,3) = SUMUI/SUMNF
            STAT(0,3)    = STAT(0,3)+SUMUI
            IF (LWIND) THEN
!             Mean observation increment - V
              STAT(JLEV,4) = SUMVI/SUMNF
              STAT(0,4)    = STAT(0,4)+SUMVI
            ENDIF
!           Mean square observation increment
            STAT(JLEV,5) = SUMSQI/SUMNF
            STAT(0,5)    = STAT(0,5)+SUMSQI
!           Maximum Increment
            MAXINC = ABS(STAT(JLEV,6))
            IF (MAXINC  >=  ABS(STAT(0,6))) THEN
              STAT(0,6) = STAT(JLEV,6)
              STAT(0,7) = STAT(JLEV,7)
              STAT(0,8) = STAT(JLEV,8)
            ENDIF
          ENDIF

        ENDIF

#if !defined(MPP)
 300    CONTINUE
#else
      enddo
#endif

!       4. Get means for all levels
!       ---------------------------
        IF (STAT(0,1) >  0.0 .AND. STAT(0,2) >  0.0) THEN

          SUMNF = STAT(0,2)
!         Mean weight
          STAT(0,2) = STAT(0,2)/STAT(0,1)
!         Mean observation increment - P*, T, U, RH
          STAT(0,3) = STAT(0,3)/SUMNF
!         Mean observation increment - V
          IF (LWIND) STAT(0,4) = STAT(0,4)/SUMNF
!         Mean square observation increment
          STAT(0,5) = STAT(0,5)/SUMNF

        ENDIF

!       If p*, convert from pascals to mb
        IF (ANAL_VAR == 1) THEN
          STAT(0,3) = STAT(0,3) * 0.01    ! Mean Increment
          STAT(0,5) = STAT(0,5) * 0.0001  ! Mean Square Increment
          STAT(0,6) = STAT(0,6) * 0.01    ! Maximum Increment
        ENDIF

!       Get RMS values if required
        IF (LRMS) THEN
          DO JLEV=0,JNANL
            IF (STAT(JLEV,5) >  0.0) THEN
              STAT(JLEV,5) = SQRT (STAT(JLEV,5))
            ENDIF
          ENDDO
        ENDIF

#if defined(MPP)
      endif  !mype=0
      CALL GC_SSYNC(NPROC,ISTAT)
#else
      ENDIF
#endif

!     5. Print out results
!     --------------------

#if defined(MPP)
      if(mype == 0)then
#endif
!---- 5.1 Print out titles

      IF (JBAND == 1) THEN
        IF (KCALL >  1 .AND. KCALL  <  6) THEN
          PRINT '(/,'' DIAGO called from '',A,T30,''Obs Type '',I5)',   &
     &    CSUBR,KTYPE
          PRINT '('' '',A9,2X,A25)', TITLE1,TITLE2
        ELSEIF (KCALL >= 6) THEN
          PRINT '(/,'' Fit to '',A,T22,''cloud (oktas)'')', CSUBR
        ELSE
          PRINT '(/,'' DIAGO called from '',A)', CSUBR
        ENDIF
      ENDIF
#if defined(GLOBAL)
      IF (LLBAND) THEN
        PRINT '(/,'' LATITUDES '',A3,'' to '',A3,5X,                    &
     &  ''No of obs '',I6/)', BANDC(JBAND),BANDC(JBAND+1),INOBS
      ELSE
        PRINT '(/,'' LATITUDES '',A3,'' to '',A3/)',BANDC(1),BANDC(4)
      ENDIF
#else
      PRINT '(/,'' WHOLE LIMITED AREA '',/)'
#endif

      IF (LRMS) THEN
        PRINT '(A,3X,A10,3X,A10,2X,A)',                                 &
     &  ' LEVEL     N  NORM FACTOR',LAB1,LAB2,                          &
     &  '    RMS       EXTREME INC & POSITION'
      ELSE
        PRINT '(A,3X,A10,3X,A10,2X,A)',                                 &
     &  ' LEVEL     N  NORM FACTOR',LAB1,LAB2,                          &
     &  'MEAN SQUARE   EXTREME INC & POSITION'
      ENDIF

!---- 5.2 Print overall means

      NS = ' '
      IF (STAT(0,7) >  0.0) THEN
        NS = 'N'
      ELSE
        NS = 'S'
      ENDIF
      WE = ' '
      IF (STAT(0,8) >  0.0) THEN
        WE = 'E'
      ELSE
        WE = 'W'
      ENDIF
      INOBSLV = INT(STAT(0,1))
      LAT     = ABS(STAT(0,7))
      LONG    = ABS(STAT(0,8))
      PRINT '(A,I6,5E13.5,F5.1,A1,F6.1,A1)','   ALL',INOBSLV,           &
     & (STAT(0,JVAR),JVAR=2,6),LAT,NS,LONG,WE

!---- 5.3 Print means for individual levels

#if defined(MPP)
       IF (JNANL >  1) THEN
#else
      IF (JNANL >  1 .AND. INOBS >  0) THEN
#endif
        DO JLEV=JNANL,1,-1
          INOBSLV = INT( STAT(JLEV,1) )
          IF (INOBSLV >  0) THEN
            NS = ' '
            IF (STAT(JLEV,7) >  0.0) THEN
              NS = 'N'
            ELSE
              NS = 'S'
            ENDIF
            WE = ' '
            IF (STAT(JLEV,8) >  0.0) THEN
              WE = 'E'
            ELSE
              WE = 'W'
            ENDIF
            LAT  = ABS(STAT(JLEV,7))
            LONG = ABS(STAT(JLEV,8))
            PRINT '(2I6,5E13.5,F5.1,A1,F6.1,A1)',                       &
     &      JLEV,INOBSLV,(STAT(JLEV,JVAR),JVAR=2,6),LAT,NS,LONG,WE
          ENDIF
        ENDDO
      ENDIF

#if defined(MPP)
      endif
      CALL GC_SSYNC(NPROC,ISTAT)
#endif
100   CONTINUE
      GOTO 1000

!   if jumping to 999 a problem has occurred but treat as non fatal
!   since this is a diagnostic routine
999   CONTINUE
1000  CONTINUE
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('DIAGO   ',4)
      RETURN
      END SUBROUTINE DIAGO
#endif
