#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE AC_INIT ------------------------------------------------
!LL
!LL  Purpose : Initialise and set up for assimilation.
!LL          : via calls to ACP_NAMEL,ACDIAG_NAMEL,NUM_OBS,SETTPS
!LL
!LL  For use on Cray Y-MP
!LL  For Cray - Global  ; Enable defs GLOBAL
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  2.8  6/11/92: portable I/O version of OPEN
!LL              : provided that IOBS_FORMAT ne 1 (R.Rawlins)
!LL  3.1 22/12/92: ROW_LENGTH added to ACP_NAMEL call.
!LL              : (Phil Andrews)
!LL
!LL   3.1  05/02/93    Portable Fortran unit no assigns
!LL                    Author: A. Dickinson    Reviewer: R. Stratton
!LL   3.2  8/7/93      Eliminate QA FORTRAN complaints    S Bell
!LL        4/8/93      and OPEN ACOBS for read only       S Bell
!LL   3.3 17/9/93      Revise LAM Inner boundary definition S Bell
!LL   3.3  23/11/93    Change NUM_OBS argument list to enable grid
!LL                    and level consistency checks G Bason
!LL   3.4  03/08/94    Pass tracer levels for setting defaults
!LL                                                        R Swinbank
!LL   4.0  20/07/95    Remove reference to old format acobs GBason
!LL   3.5  24/03/95    Changed OPEN to FILE_OPEN  P.Burton
!LL   4.1  9/5/96      Check values in ACPARM are large enough. S Bell
!    4.1  18/06/96: Changes to cope with changes in STASH addressing
!                   Author D.M. Goddard.
!    4.2 25/11/96: T3E mods Stuart Bell
!LL
!LL  4.1  4/09/96:  Port to CRAY T3E  Deborah Salmond
!LL  4.5  5/8/98:  Increase size of FILENAME.  Stuart Bell
!LL  5.2  22/3/01:  amend glsize dimensions   B Macpherson
!    5.3  09/07/01:  amend XLATN definition   B Macpherson
!                    remove ak,bk from call to acp_namel
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LLEND------------------------------------------------------------------
!*L  Arguments:---------------------------------------------------
      SUBROUTINE AC_INIT (P_LEVELS, Q_LEVELS, BL_LEVELS, TR_LEVELS,     &
     &                    P_ROWS, U_ROWS, ROW_LENGTH,                   &
     &                    TNDVMAX, NOBSMAX, TIMESTEP,                   &
     &                    BASIS_TIME_YY, BASIS_TIME_MM, BASIS_TIME_DD,  &
     &                    BASIS_TIME_HH, BASIS_TIME_MIN,                &
     &                    REALHD1, REALHD2, REALHD3,                    &
     &                    REALHD4, REALHD5, REALHD6,                    &
     &                    AK, BK,                                       &
#include "argppx.h"
     &                    ICODE, CMESSAGE)

!     Do not multi-task this routine
!FPP$ NOCONCUR R

      IMPLICIT NONE
#include "acparm.h"
#include "cenvir.h"
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"

      INTEGER P_LEVELS
      INTEGER Q_LEVELS
      INTEGER BL_LEVELS
      INTEGER TR_LEVELS
      INTEGER P_ROWS
      INTEGER U_ROWS
      INTEGER ROW_LENGTH
      INTEGER TNDVMAX
      INTEGER TNDVMAX_TOTAL
      INTEGER NOBSMAX
      INTEGER NOBSMAX_TOTAL
      REAL TIMESTEP
      INTEGER BASIS_TIME_YY
      INTEGER BASIS_TIME_MM
      INTEGER BASIS_TIME_DD
      INTEGER BASIS_TIME_HH
      INTEGER BASIS_TIME_MIN
      REAL REALHD1,REALHD2,REALHD3,REALHD4,REALHD5,REALHD6
      REAL AK(P_LEVELS),BK(P_LEVELS)
      CHARACTER*256 CMESSAGE
      INTEGER ICODE
#if !defined(MPP)
      INTEGER  mype
      PARAMETER (mype = 0 ) ! always zero in non-MPP code
#endif
!-INTENT=IN--------------------------------------------------------
!     P_LEVELS      : TOTAL NUMBER OF LEVELS
!     Q_LEVELS      : TOTAL NUMBER OF wet LEVELS
!     BL_LEVELS     : TOTAL NUMBER OF boundary layer LEVELS
!     TR_LEVELS     : TOTAL NUMBER OF tracer LEVELS
!     ROW_LENGTH    : NUMBER OF POINTS ON ROW
!     P_ROWS        : NUMBER OF ROWS (FOR PSTAR)
!     U_ROWS        : NUMBER OF ROWS (FOR wind)
!     TIMESTEP      : TIMESTEP IN SECONDS
!     BASIS_TIME_## : DEFINES DATA TIME
!     REALHD1-6     : DEFINES HORIZONTAL GRID
!     AK,BK         : DEFINES VERTICAL GRID
!-INTENT=OUT-----------------------------------------------------
!     TNDVMAX       : MAX SIZE OF OBS ARRAY
!     NOBSMAX       : MAX NUMBER OF OBS
!     ICODE         : NON ZERO FOR FAILURE
!     CMESSAGE      : REASON FOR FAILURE
!*---------------------------------------------------------------------
      EXTERNAL ACP_NAMEL,ACDIAG_NAMEL,NUM_OBS,SETTPS
      EXTERNAL GET_FILE
!L---------------------------------------------------------------------
#if defined(MPP)
#include "parvars.h"
#endif
!     AC Comdecks
#include "comacp.h"
#include "comacdg.h"
#include "comag.h"
#include "commg.h"
#include "comobs.h"
!L---------------------------------------------------------------------
!     UM Constant Comdecks
#include "c_pi.h"
!-----------------------------------------------------------------------
      INTEGER JFILE
      INTEGER ICCODE
      INTEGER IUNIT7
      CHARACTER*256 FILENAME
#if defined(MPP)
      INTEGER ROW_LENGTH_GLOBAL,P_ROWS_GLOBAL
      ROW_LENGTH_GLOBAL=glsize(1,fld_type_p)
      P_ROWS_GLOBAL=glsize(2,fld_type_p)
#endif
!-----------------------------------------------------------------------

      if(mype == 0) PRINT *, ' IN AC_INIT'

!L 1. Initialise Unit Numbers used in AC
! Unit No for first AC Observation File
      OBS_UNITNO = 70

! Unit No for CACHE FILE to store observations between timesteps
! Stores arrays OBS and OBS_FLAG
      IUNITNO = 15

!     Unit No for printed output used to inform operators of
!     data problems
      IUNIT7 = 7

!L 2. Set up COMMG Variables and check dimensions
      DLONG  = REALHD1
      DLAT   = REALHD2
      XLATN  = REALHD3 + (P_ROWS_GLOBAL-1)*DLAT
      XLONGW = REALHD4

      IF(P_LEVELS >   P_LEVELS_MAX)THEN
       ICODE=1
       CMESSAGE = ' ACINIT: P_LEVELS_MAX too small'
       if(mype == 0) PRINT *,  CMESSAGE,                                &
     &      ' P_LEVELS=',P_LEVELS,' P_LEVELS_MAX=',P_LEVELS_MAX
       GOTO 999
      ENDIF

      IF(Q_LEVELS >   Q_LEVELS_MAX)THEN
       ICODE=1
       CMESSAGE = ' ACINIT: Q_LEVELS_MAX too small'
       if(mype == 0) PRINT *,  CMESSAGE,                                &
     &      ' Q_LEVELS=',Q_LEVELS,' Q_LEVELS_MAX=',Q_LEVELS_MAX
       GOTO 999
      ENDIF

      IF(P_ROWS >   P_ROWS_MAX)THEN
       ICODE=1
       CMESSAGE = ' ACINIT: P_ROWS_MAX too small'
       if(mype == 0) PRINT *,  CMESSAGE,                                &
     &      ' P_ROWS=',P_ROWS,' P_ROWS_MAX=',P_ROWS_MAX
       GOTO 999
      ENDIF

      IF(ROW_LENGTH >   ROW_LENGTH_MAX)THEN
       ICODE=1
       CMESSAGE = ' ACINIT: ROW_LENGTH_MAX too small'
       if(mype == 0) PRINT *,  CMESSAGE,                                &
     &      ' ROW_LENGTH=',ROW_LENGTH,' ROW_LENGTH_MAX=',ROW_LENGTH_MAX
       GOTO 999
      ENDIF

! Make sure western boundary longitude in range 0 - 360 degrees
      IF (XLONGW <  0.0) XLONGW = XLONGW + 360.0

#if !defined(GLOBAL)
! Real lat/lon of pseudo N. pole in degrees
      ELFPLAT = REALHD5
      ELFPLON = REALHD6

#endif
      WRITE (6,'(/,A,(T25,10F10.6))') ' DLAT,DLONG,XLATN,XLONGW',       &
     &                                  DLAT,DLONG,XLATN,XLONGW

!L 3. ACP namelist. Set defaults, read in and process.
! DEPENDS ON: acp_namel
      CALL ACP_NAMEL (P_LEVELS, Q_LEVELS, BL_LEVELS, TR_LEVELS,         &
     &  P_ROWS, U_ROWS, ROW_LENGTH,                                     &
     &  TIMESTEP, ICODE, CMESSAGE)

      IF (ICODE >  0) GO TO 999

!L 4. ADIAG namelist. Set defaults, read in and process.
! DEPENDS ON: acdiag_namel
      CALL ACDIAG_NAMEL (ICODE,CMESSAGE)
      IF (ICODE >  0) GO TO 999

!L 5. Open Cache file (unit 15)

      if(mype == 0)then
      CALL GET_FILE(IUNITNO,FILENAME,256,ICODE)
        OPEN(IUNITNO,FILE=FILENAME,FORM='UNFORMATTED')
      CALL GET_FILE(IUNIT7,FILENAME,256,ICODE)
        OPEN(IUNIT7,FILE=FILENAME)
      endif

!L 6. Read in AC Obs Files and compute NOBSMAX and TNDVMAX

! DEPENDS ON: num_obs
      CALL NUM_OBS (NO_OBS_FILES,NOBSMAX,TNDVMAX,P_LEVELS,Q_LEVELS,     &
     &              P_ROWS,ROW_LENGTH,AK,BK,REALHD1,REALHD2,            &
     &              REALHD3,REALHD4,REALHD5,REALHD6,                    &
#include "argppx.h"
     &              ICODE,CMESSAGE)
      IF (ICODE >  0) GO TO 999

!L 7. Set up list of obs types and groups to be processed this run
! DEPENDS ON: settps
      CALL SETTPS (ICODE,CMESSAGE)
      IF (ICODE >  0) GO TO 999

!L 8. Set up Analysis Grid in COMAG
! Initialise COMAG variables which don't change during a run.

! Model Grid Spacings
      DLATMG  = DLAT *PI_OVER_180
      DLONGMG = DLONG*PI_OVER_180

! First row/point on p*/theta grid
      ROW1MGTH = ( 90.0 -XLATN )*PI_OVER_180
      PT1MGTH  = XLONGW*PI_OVER_180

! Convert Lat at which pts/row starts decreasing to co-lat/radians.
      AGLATDEC = (90.0-AGLATDEC)*PI_OVER_180

! Width of area according to model grid dimensions.
#if defined(MPP)
#if defined(GLOBAL)
      AGROWLEN = DLONGMG*ROW_LENGTH_GLOBAL
#else
      AGROWLEN = DLONGMG*(ROW_LENGTH_GLOBAL-1)
#endif
#else
#if defined(GLOBAL)
      AGROWLEN = DLONGMG*ROW_LENGTH
#else
      AGROWLEN = DLONGMG*(ROW_LENGTH-1)
#endif
#endif

! First row/point on wind grid.
! Model's wind grid is staggered down 1/2 row, across 1/2 pt.
      ROW1MGUV = ROW1MGTH + 0.5*DLATMG
      PT1MGUV  = PT1MGTH  + 0.5*DLONGMG

! The COMAG variables DEF_AGRES_ROWS and DEF_AGRES_PTS are
! initialised in DEF_GROUP. Use &ACP namelist arrays AGRES_ROWS
! and AGRES_PTS to change initialised values.

! The remaining COMAG variables are set in SETAG.

!L 9. Set up COMOBS Variables
! Time Interval (mins) between reading AC Obs files
#if defined(GLOBAL)
      TIMEINT = 360.0
#else
      TIMEINT = 180.0
#endif

! Reference Time/Date which is start of assimilation
      OBS_REF_YY  = BASIS_TIME_YY
      OBS_REF_MM  = BASIS_TIME_MM
      OBS_REF_DD  = BASIS_TIME_DD
      OBS_REF_HH  = BASIS_TIME_HH
      OBS_REF_MIN = BASIS_TIME_MIN

! Set up Latitudes and Longitudes of grid boundaries
! for observations to be used in assimilation.

#if defined(GLOBAL)
! For global assimilations, use all observations. For wind obs
! the First and Last row of the model wind grid are the N
! and S limits - controlled in RDOBS2.
      OBS_LAT_N  = 90.01
      OBS_LAT_S  =-90.01
      OBS_LONG_E = 180.01
      OBS_LONG_W =-180.01
#else
! For Limited area assimilations, reject observations within
! one grid length of boundary.
      OBS_LAT_N  = XLATN + 0.01 - DLAT
#if !defined(MPP)
      OBS_LAT_S  = XLATN - 0.01 - (P_ROWS-2)*DLAT
#else
      OBS_LAT_S  = XLATN - 0.01 - (P_ROWS_GLOBAL-2)*DLAT
#endif
      OBS_LONG_W = XLONGW - 0.01 + DLONG
#if !defined(MPP)
      OBS_LONG_E = XLONGW + 0.01 + (ROW_LENGTH-2)*DLONG
#else
      OBS_LONG_E = XLONGW + 0.01 + (ROW_LENGTH_GLOBAL-2)*DLONG
#endif

! After rotation of the Lat/Long of Obs to ELF co-ords,
! longitude values will be in range 0-360 degrees, so
! make sure that boundary values are consistent. Note that
! if this leaves ZLONMN > ZLONMX, the test for obs in the area
! will assume that the Limited Area grid straddles the Meridian
! between ZLONMN and ZLONMX.
      IF (OBS_LONG_W  <   0.0) THEN
        OBS_LONG_W = OBS_LONG_W + 360.0

      ELSEIF (OBS_LONG_W  >   360.0) THEN
        OBS_LONG_W = OBS_LONG_W - 360.0

      ENDIF

      IF (OBS_LONG_E  <   0.0) THEN
        OBS_LONG_E = OBS_LONG_E + 360.0

      ELSEIF (OBS_LONG_E  >   360.0) THEN
        OBS_LONG_E = OBS_LONG_E - 360.0

      ENDIF
#endif

! Initialise time for next read in of observation files.
! Forces read on first timestep.
      TIMENEXT = -1440.0

 999  CONTINUE
      RETURN
      END SUBROUTINE AC_INIT
#endif
