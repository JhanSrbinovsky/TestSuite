#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINES ACP_NAMEL,ACDIAG_NAMEL  -------------------------------
!LL
!LL  Purpose : Read in AC Namelist (&ACP) and process.
!LL
!LL ACDIAG_NAMEL:
!LL Set defaults for ACDIAG namelist variables.
!LL                  Read in namelist and process.
!LL
!LL  For use on Cray Y-MP
!LL  For Cray - Global  ; Enable defs GLOBAL
!LL
!LL S.Bell      <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  3.1 18/11/92:New variable MTHIN205,MVINT205
!LL              : (Stuart Bell)
!LL  3.1 9/12/92 : WB_CUT_OFF_LEV,WB_PSTAR_METHOD removed from namelist
!LL              : WB_THETA_INC added to choose P* and theta or just p*
!LL              : WB_LAT_CC added to store latitudinal correlation
!LL              : coefficient for WINDBAL derived increments.
!LL              : WB_LAND_FACTOR,WB_LAND_SCALE,WB_CC_LAM,
!LL              : WB_CC_N, WB_CC_S, WB_START_LEV, WB_END_LEV added
!LL              : to namelist to specify WINDBAL correlation coeffs.
!LL              : (Phil Andrews,Stuart Bell)
!LL  3.1 18/2/93 : New defaults set under LAC_MES switch. (B Macpherson)
!LL  3.2 28/6/93 : New defaults for WINDBAL control       (S Bell)
!LL  3.2 10/7/93 : Abort if read problems with namelist.  (S Bell)
!LL              : New variables OBTHIN (replacing MTHINxxx), MGLOSSFN
!LL              : and WB_THETA_SURF               S Bell
!LL  3.3 1/12/93 : Change some Mes defaults    Bruce Macpherson
!LL  3.3 5/11/93 : New variables for limited area WINDBAL
!LL              : & new subroutine to calc subset of limited area grid
!LL              : suitable for using multigrid. (PhilAndrews)
!LL  3.3 1/12/93 : THRESH-.. defaults specified   B Macpherson
!LL  3.4 03/08/94 : Pass TR_LEVELS for setting tracer defaults
!LL                                               R Swinbank
!LL  3.4 01/09/94 : Sharper vertical correlation defaults for UARS
!LL                                               R Swinbank
!LL  3.3+ 26/4/94: Update WINDBAL defaults        S Bell
!LL  3.4 9/09/94 : Remove L_REINITQC, add CSCALE_VERT_AERO
!LL              : Revise mes default for vert scaling of horiz scale
!LL              : Specify defaults for RADAR_RANGE,LRADAR,LHYDROL,
!LL              : NORTH/SOUTHLAT,WEST/EASTLON,L_MOPS_EQUALS_RH
!LL              :                     B Macpherson
!LL  4.0  21/2/95: Add LCHECK_GRID :control call to CHECK_OBS  G Bason
!LL  4.0  3/5/95 : Specify defaults for LHN_RANGE, L_LHN, L_LHN_SCALE,
!LL              : L_LHN_SEARCH, LHN_DIAG, RADAR_LAT, RADAR_LON,
!LL              : RADAR_RANGE_MAX, EPSILON_LHN, RELAX_CF_LHN,
!LL              : F1_506 , F2_506 , F3_506 , L_506_OBERR.
!LL              :                              Chris Jones.
!LL  4.0 5/07/95 : Revise CSCALE_VERT (a)GL(T1), (b)LA(T1) and (V1)
!LL              : and set all LAM extratropic=tropic S Bell
!LL  4.0 10/07/95 : Add variable L_OBS_CHECK to ACP G.Bason
!LL  4.0 20/07/95 : Set OBS_FORMAT default to be 2 and error trap#
!LL                 use of any other format G.Bason
!LL  4.0 21/3/95 : Add NON_DIV_COR_10M, DEF_CSCALE_VERT_MES.
!LL              : Revise defaults for LRADAR, L_MOPS_EQUALS_RH
!LL              :                                Bruce M
!LL  4.1  6/6/96 : Set defaults for new LHN variables: ALPHA_LHN,
!LL              :     LHN_LIMIT, FI_SCALE_LHN, NPASS_RF_LHN,
!LL              :     L_LHN_LIMIT, L_LHN_FACT, L_LHN_FILT (Chris Jones)
!LL  4.1  6/6/96 : Set use of all radars to TRUE (Chris Jones)
!    4.2 25/11/96: T3E mods Stuart Bell
!LL
!LL  4.1  4/09/96:  Port to CRAY T3E  Deborah Salmond
!LL  4.4  7/08/97:  Allow OBS_FORMAT=3  Adam Maycock
!LL  4.5  6/10/98:  change tropical correlation scale Stuart Bell
!LL  5.2  22/3/01:  amend glsize dimensions   B Macpherson
!LL  5.3  07/06/01:  Change print to write for namelist.  A Van der Wal
!    5.3  29/6/01:  remove redundant code B Macpherson
!    6.2  22/08/05:  Restructure defs for ACP namelist. P.Selwood
!    6.2  15/08/05: Fixes for free-format. P.Selwood
!    6.2  23/02/05:  Add new logical REMOVE_NEG_LH. Mark Dixon
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Logical system components covered:
!LL
!LL  Project Task : P3
!LL
!LL  External documentation:
!LL
!LLEND -------------------------------------------------------------
!LLEND
!*L  Arguments:
      SUBROUTINE ACP_NAMEL (P_LEVELS, Q_LEVELS, BL_LEVELS, TR_LEVELS,   &
     & P_ROWS, U_ROWS,ROW_LEN,                                          &
     & TIMESTEP, ICODE, CMESSAGE)

! Do not multi-task this routine.
!FPP$ NOCONCUR R

      IMPLICIT NONE
#include "acparm.h"

! Imported variables (INTENT=IN):L
      INTEGER     P_LEVELS            ! Total number of levels
      INTEGER     Q_LEVELS            ! Total number of wet levels
      INTEGER     BL_LEVELS           ! total number of boundary layer
                                      ! levels
      INTEGER     TR_LEVELS           ! total number of tracer levels
      INTEGER     P_ROWS              ! Number of rows (for pstar)
      INTEGER     U_ROWS              ! Number of rows (for wind)
      INTEGER     ROW_LEN             ! Length of each row.

      REAL        TIMESTEP            ! timestep in seconds

! Exported variables (INTENT=OUT)
      INTEGER ICODE ! Non zero for failure

      CHARACTER*256 CMESSAGE          ! Reason for failure
!*
!L---------------------------------------------------------------------
!L UM Comdecks
#include "c_mdi.h"
#include "c_a.h"
!L---------------------------------------------------------------------
!L AC Comdecks
#include "comacp.h"
#include "comag.h"
#if defined(MPP)
#include "parvars.h"
      INTEGER P_ROWS_GLOBAL,U_ROWS_GLOBAL
#endif
!L---------------------------------------------------------------------
#if defined(GLOBAL)
      REAL    GEOWT_NH
      REAL    GEOWT_SH
      REAL    NUDGE_NH(NOBTYPMX)
      REAL    NUDGE_TR(NOBTYPMX)
      REAL    NUDGE_SH(NOBTYPMX)
      REAL    SCFACT_NH
      REAL    SCFACT_TR
      REAL    SCFACT_SH
      REAL    WB_CC_N                     ! Northern hemisphere
                                          ! corln. coeff. for WINDBAL.
      REAL    WB_CC_S                     ! Southern hemisphere
                                          ! corln. coeff. for WINDBAL.

#else
      REAL    GEOWT_LAM
      REAL    NUDGE_LAM(NOBTYPMX)
      REAL    WB_CC_LAM                   ! WINDBAL correlation

#endif
      REAL    GEOWT_VERT(21)
      REAL    CSCALE_VERT(21,4)
      REAL    DEF_CSCALE_VERT_MES(21,4)
      REAL    SCFACT_VERT(P_LEVELS_MAX)
      REAL    FI_VAR_FACTOR(NOBTYPMX)
      REAL    TEMP_COORD
      REAL    CSCALE_VERT_AERO

      INTEGER AC_ORDER     (NOBTYPMX)
      INTEGER NO_ITERATIONS(NOBTYPMX)
      INTEGER INTERVAL_ITER(NOBTYPMX)
      INTEGER N_ANAL_LEVS  (NOBTYPMX)
      INTEGER N_WT_LEVS    (NOBTYPMX)
      INTEGER AGRES_ROWS   (NOBTYPMX)
      INTEGER AGRES_PTS    (NOBTYPMX)
      INTEGER MODE_HANAL   (NOBTYPMX)
      INTEGER TIMEB        (NOBTYPMX)
      INTEGER TIMEA        (NOBTYPMX)
      INTEGER TGETOBB      (NOBTYPMX)
      INTEGER TGETOBA      (NOBTYPMX)
      INTEGER RADINF       (NOBTYPMX)
      INTEGER OBTHIN       (NOBTYPMX)
      INTEGER CSCALE_START (NOBTYPMX)
      INTEGER CSCALE_OBTIME(NOBTYPMX)
      INTEGER CSCALE_END   (NOBTYPMX)
      INTEGER J,JJ,JLEV,JROW,JTYP,JOBT
      INTEGER NCOUNT
      INTEGER WB_START_LEV                 ! First level for non zero
                                           ! vertical correlation coeff.
      INTEGER WB_END_LEV                   ! Last  level for non zero
                                           ! vertical correlation coeff.

#if !defined(MPP)
      INTEGER  mype
      PARAMETER (mype = 0 ) ! always zero in non-MPP code
#endif
      EXTERNAL DEF_GROUP,DEF_TYPE,GROUP_DEP_VAR,TYPE_DEP_VAR

! Vertical correlation scale defaults; 21 values at 50hPa intervals
! for extratropics T,tropics T,for extratropics V,tropics V
#if defined(GLOBAL)
      DATA CSCALE_VERT/                                                 &
     &       12.0, 8.0, 8.0, 8.0, 7.0, 6.0, 5.0, 4.5, 4.0, 4.0,         &
     &        4.0, 4.5, 5.0, 5.0, 5.5, 6.0, 5.5, 4.0, 3.5, 2.0, 2.0,    &
     &       12.0, 8.0, 8.0, 8.0, 7.0, 6.0, 5.5, 5.0, 4.5, 4.0,         &
     &        3.5, 3.0, 3.0, 3.0, 3.0, 4.0, 3.5, 3.0, 2.5, 1.5, 1.5,    &
     &        6.0, 6.0, 5.0, 5.0, 5.0, 5.0, 4.0, 3.5, 3.5, 3.0,         &
     &        3.0, 2.5, 2.5, 2.5, 2.5, 3.0, 2.5, 2.0, 2.0, 1.5, 1.5,    &
     &        6.5, 6.5, 6.0, 6.0, 5.5, 5.0, 5.0, 4.0, 3.5, 3.5,         &
     &        3.5, 3.5, 3.0, 3.0, 3.0, 3.0, 3.0, 2.5, 2.0, 1.5, 1.5/
#else
! for LAM Tropics and Extra-tropics set the same
      DATA CSCALE_VERT/                                                 &
     &       12.0, 8.0, 8.0, 8.0, 7.0, 6.0, 5.0, 4.5, 4.0, 4.0,         &
     &        4.0, 4.5, 5.0, 5.0, 5.5, 6.0, 5.5, 4.0, 3.5, 2.0, 2.0,    &
     &       12.0, 8.0, 8.0, 8.0, 7.0, 6.0, 5.0, 4.5, 4.0, 4.0,         &
     &        4.0, 4.5, 5.0, 5.0, 5.5, 6.0, 5.5, 4.0, 3.5, 2.0, 2.0,    &
     &        9.0, 6.0, 5.0, 5.0, 5.0, 5.0, 4.0, 3.5, 3.5, 3.0,         &
     &        3.0, 2.5, 2.5, 2.5, 2.5, 3.0, 2.5, 2.0, 2.0, 1.5, 1.5,    &
     &        9.0, 6.0, 5.0, 5.0, 5.0, 5.0, 4.0, 3.5, 3.5, 3.0,         &
     &        3.0, 2.5, 2.5, 2.5, 2.5, 3.0, 2.5, 2.0, 2.0, 1.5, 1.5/
#endif

      DATA DEF_CSCALE_VERT_MES/                                         &
     &       18.0,18.0,18.0, 8.0, 7.0, 6.0, 5.0, 4.5, 4.0, 4.0,         &
     &        4.0, 4.5, 5.0, 5.0, 5.5, 6.0, 5.5, 4.0, 3.5, 2.0, 2.0,    &
     &       18.0,18.0,18.0, 8.0, 7.0, 6.0, 5.0, 4.5, 4.0, 4.0,         &
     &        4.0, 4.5, 5.0, 5.0, 5.5, 6.0, 5.5, 4.0, 3.5, 2.0, 2.0,    &
     &       27.0,27.0,27.0,12.0, 8.0, 6.0, 4.0, 3.5, 3.5, 3.0,         &
     &        3.0, 2.5, 2.5, 2.5, 2.5, 3.0, 2.5, 2.0, 2.0, 1.5, 1.5,    &
     &       27.0,27.0,27.0,12.0, 8.0, 6.0, 4.0, 3.5, 3.5, 3.0,         &
     &        3.0, 2.5, 2.5, 2.5, 2.5, 3.0, 2.5, 2.0, 2.0, 1.5, 1.5/
!     NAMELIST ACP
!     ------------
#if defined(GLOBAL)
      NAMELIST /ACP/                                                    &
     & AC_ORDER, AC_OBS_TYPES,                                          &
     & AGRES_ROWS, AGRES_PTS, MIN_AGPTS, AGLATDEC,                      &
     & TIMEA,TIMEB,TGETOBA,TGETOBB,                                     &
     & N_ANAL_LEVS,N_WT_LEVS,MODE_HANAL,FI_VAR_FACTOR,                  &
     & OBTIME_NOM, LTIMER_AC, LGEO, LHYDR, LSYN, LAC_UARS,              &
     & VERT_FILT,NO_ITERATIONS,INTERVAL_ITER,                           &
     & IOMITOBS,OBTHIN,SPEED_LIMIT305,MGLOSSFN,                         &
     & MHCORFN,MACDIAG,MWTFN,MDATADFN,MRAMPFN,                          &
     & FI_SCALE,FI_SCALE_FACTOR,DF_COEFF,DF_SCALE,DF_SCALE_LEV,         &
     & NPASS_RF,TIMEF_START,TIMEF_OBTIME,TIMEF_END,                     &
     & CSCALE_VERT,VERT_CUTOFF_SL,VERT_CUTOFF_BW,VERT_CUTOFF_BH,        &
     & CSCALE_START,CSCALE_OBTIME,CSCALE_END,RADINF,                    &
     & SCFACT_VERT,NO_SCFACT,                                           &
     & MGEOWT,GEOWT_VERT,MVINT205,                                      &
     & NON_DIV_COR, LAC_MES,                                            &
     & NON_DIV_COR_10M,                                                 &
     & OBS_FORMAT,NO_OBS_FILES,DIAG_RDOBS,NPROG,L_OBS_CHECK,            &
     & LWBAL_SF, LWBAL_UA, WB_THETA_UA, WB_THETA_SF,                    &
     & WB_START_LEV, WB_END_LEV, WB_LAND_SCALE, WB_LAND_FACTOR,         &
     & THRESH_DL, THRESH_LM, THRESH_MH, THRESH_RMSF,                    &
     & RADAR_RANGE  ,  LRADAR  ,  LHYDROL ,   L_LATLON_PRVER ,          &
     & NORTHLAT ,  SOUTHLAT ,  WESTLON ,  EASTLON ,                     &
     & LHN_RANGE ,  L_LHN , L_LHN_SCALE ,                               &
     & L_LHN_SEARCH , LHN_DIAG ,                                        &
     & RADAR_LAT , RADAR_LON , RADAR_RANGE_MAX ,                        &
     & EPSILON_LHN , ALPHA_LHN , RELAX_CF_LHN , LHN_LIMIT ,             &
     & F1_506 , F2_506 , F3_506 ,                                       &
     & L_506_OBERR , L_VERIF_RANGE , L_LHN_LIMIT , L_LHN_FACT ,         &
     & L_LHN_FILT , FI_SCALE_LHN , NPASS_RF_LHN ,                       &
     & L_MOPS_EQUALS_RH , CSCALE_VERT_AERO, LCHECK_GRID,                &
     & GEOWT_NH,GEOWT_SH,                                               &
     & NUDGE_NH,NUDGE_TR,NUDGE_SH,                                      &
     & SCFACT_NH,SCFACT_TR,SCFACT_SH,                                   &
     & TROPLAT,TROPINT,                                                 &
     & WB_CC_N, WB_CC_S, REMOVE_NEG_LH, USE_CONV_IN_MOPS
#else
      NAMELIST /ACP/                                                    &
     & AC_ORDER, AC_OBS_TYPES,                                          &
     & AGRES_ROWS, AGRES_PTS, MIN_AGPTS, AGLATDEC,                      &
     & TIMEA,TIMEB,TGETOBA,TGETOBB,                                     &
     & N_ANAL_LEVS,N_WT_LEVS,MODE_HANAL,FI_VAR_FACTOR,                  &
     & OBTIME_NOM, LTIMER_AC, LGEO, LHYDR, LSYN, LAC_UARS,              &
     & VERT_FILT,NO_ITERATIONS,INTERVAL_ITER,                           &
     & IOMITOBS,OBTHIN,SPEED_LIMIT305,MGLOSSFN,                         &
     & MHCORFN,MACDIAG,MWTFN,MDATADFN,MRAMPFN,                          &
     & FI_SCALE,FI_SCALE_FACTOR,DF_COEFF,DF_SCALE,DF_SCALE_LEV,         &
     & NPASS_RF,TIMEF_START,TIMEF_OBTIME,TIMEF_END,                     &
     & CSCALE_VERT,VERT_CUTOFF_SL,VERT_CUTOFF_BW,VERT_CUTOFF_BH,        &
     & CSCALE_START,CSCALE_OBTIME,CSCALE_END,RADINF,                    &
     & SCFACT_VERT,NO_SCFACT,                                           &
     & MGEOWT,GEOWT_VERT,MVINT205,                                      &
     & NON_DIV_COR, LAC_MES,                                            &
     & NON_DIV_COR_10M,                                                 &
     & OBS_FORMAT,NO_OBS_FILES,DIAG_RDOBS,NPROG,L_OBS_CHECK,            &
     & LWBAL_SF, LWBAL_UA, WB_THETA_UA, WB_THETA_SF,                    &
     & WB_START_LEV, WB_END_LEV, WB_LAND_SCALE, WB_LAND_FACTOR,         &
     & THRESH_DL, THRESH_LM, THRESH_MH, THRESH_RMSF,                    &
     & RADAR_RANGE  ,  LRADAR  ,  LHYDROL ,   L_LATLON_PRVER ,          &
     & NORTHLAT ,  SOUTHLAT ,  WESTLON ,  EASTLON ,                     &
     & LHN_RANGE ,  L_LHN , L_LHN_SCALE ,                               &
     & L_LHN_SEARCH , LHN_DIAG ,                                        &
     & RADAR_LAT , RADAR_LON , RADAR_RANGE_MAX ,                        &
     & EPSILON_LHN , ALPHA_LHN , RELAX_CF_LHN , LHN_LIMIT ,             &
     & F1_506 , F2_506 , F3_506 ,                                       &
     & L_506_OBERR , L_VERIF_RANGE , L_LHN_LIMIT , L_LHN_FACT ,         &
     & L_LHN_FILT , FI_SCALE_LHN , NPASS_RF_LHN ,                       &
     & L_MOPS_EQUALS_RH , CSCALE_VERT_AERO, LCHECK_GRID,                &
     & WB_LonOffset, WB_LonPts, WB_LatOffset, WB_LatPts,                &
     & GEOWT_LAM, NUDGE_LAM, WB_CC_LAM, REMOVE_NEG_LH, USE_CONV_IN_MOPS
#endif

!-----------------------------------------------------------------------
#if defined(MPP)
      P_ROWS_GLOBAL=glsize(2,fld_type_p)
      U_ROWS_GLOBAL=glsize(2,fld_type_p)-1
#endif
!L    1. ACP NAMELIST
!L    ---------------
!L 1.1 Set namelist defaults for variables defined by UI

      DO JOBT = 1, NOBTYPMX
! Group dependent arrays in namelist.
        AC_ORDER(JOBT)       = IMDI

! Obs Type dependent arrays in namelist.
        AC_OBS_TYPES(JOBT)   = 0
      END DO

! Switch  for UARS assimilation.
      LAC_UARS = .FALSE.

! Switch  for Mesoscale UM assimilation.
      LAC_MES = .FALSE.
!  Latent Heat Nudging
      L_LHN = .FALSE.

! Format of ac observation file
! 2 - Current ; 1 - Old Format now redundant
      OBS_FORMAT = 2

! Number of observation files to be used
      NO_OBS_FILES = 2

! Switch for 'non-oper no obs' test in AC1
      L_OBS_CHECK = .TRUE.

!L 1.2 Read in Namelist set by UI
      REWIND 5
      READ (5, ACP, END=120, ERR=120)

! QA Fortran recommends use of IOSTAT but it doesn't work with namelist
! Jump over error trap
      GOTO 121

! NAMELIST read error trap:
120   CONTINUE
        ICODE    = 1
        CMESSAGE = ' ACP_NAMEL : Error reading 1st ACP namelist'
        if(mype == 0) PRINT *, CMESSAGE


       ! Jump to end of routine:
        GOTO 999

! No NAMELIST read error so continue...
121   CONTINUE

!L 1.3 Set namelist defaults for variables defined by user
      DO JOBT = 1, NOBTYPMX
       ! Group dependent arrays in namelist.

        NO_ITERATIONS(JOBT)  = IMDI
        INTERVAL_ITER(JOBT)  = IMDI
        N_ANAL_LEVS(JOBT)    = IMDI
        N_WT_LEVS(JOBT)      = IMDI
        AGRES_ROWS(JOBT)     = IMDI
        AGRES_PTS (JOBT)     = IMDI
        MODE_HANAL(JOBT)     = IMDI
        FI_VAR_FACTOR(JOBT)  = RMDI

#if defined(GLOBAL)
        NUDGE_NH(JOBT)       = RMDI
        NUDGE_TR(JOBT)       = RMDI
        NUDGE_SH(JOBT)       = RMDI

#else
        NUDGE_LAM(JOBT)      = RMDI

#endif

! Obs Type dependent arrays in namelist.
        TIMEB(JOBT)          = 0
        TIMEA(JOBT)          = 0
        TGETOBB(JOBT)        = 0
        TGETOBA(JOBT)        = 0
        RADINF(JOBT)         = 0
        OBTHIN(JOBT)         = 0
        CSCALE_START(JOBT)   = 0
        CSCALE_OBTIME(JOBT)  = 0
        CSCALE_END(JOBT)     = 0
        NO_SCFACT(JOBT)      = 0

      END DO
      NO_SCFACT(1) = 305

#if defined(GLOBAL)

! Global version:
! Nominal ob time relative to start of assm for synoptic insertion
        OBTIME_NOM = 300.0

! Correlation Scale Factors (Horizontal)
      IF(LAC_UARS)THEN
        SCFACT_NH = 1.0
        SCFACT_TR = 1.0
        SCFACT_SH = 1.0

      ELSE
        SCFACT_NH = 1.0
        SCFACT_TR = 1.3
        SCFACT_SH = 1.5

      ENDIF

! Scaling for geostrophic increments
      IF(LAC_UARS)THEN
        GEOWT_NH  = 0.7
        GEOWT_SH  = 0.7

      ELSE
        GEOWT_NH  = 0.5
        GEOWT_SH  = 0.7

      END IF

! Latitude at which mid-latitude Analysis Parameters
! begin to change to tropical values.
        TROPLAT=30.0

! Interval (in degrees Lat) over which Relaxation Coefficients
! and Correlation Scale Factor change.
        TROPINT = 10.0

! Latitude at which model and Analysis Grid become identical.
#if defined(MPP)
! analysis grid=model grid in MPP mode
        AGLATDEC = 90.0
#else
        AGLATDEC = 50.0
#endif

#else
! Limited Area Version

! Nominal ob time relative to start of assm for synoptic insertion
        OBTIME_NOM = 150.0

! Scaling for geostrophic increments
        GEOWT_LAM = 0.5

! Force Model and Analysis Grid to be identical for Limited Area
        AGLATDEC = 90.0
#endif

! Parameters for time factor
      IF (LAC_UARS) THEN
        TIMEF_START  = 0.01
        TIMEF_END    = 0.01

      ELSE
        TIMEF_START  = 0.05
        TIMEF_END    = 0.05

      END IF

      TIMEF_OBTIME = 1.0

! Min no. of points on analysis grid row
      MIN_AGPTS = 14

! Geostrophic & hydrostatic incs.
      LGEO  = .TRUE.
      LHYDR = .TRUE.

! Hydrology incs.
      LHYDROL = .TRUE.

! Variable (rh or cloud fraction) used in MOPS obs
      L_MOPS_EQUALS_RH = .FALSE.

! Lat/lon area selection for precip verif
      L_LATLON_PRVER = .FALSE.
! Default lat/lon area when above logical set to .TRUE.
! (this area is larger than current FRONTIERS area)
      NORTHLAT = 60.0
      SOUTHLAT = 45.0
      WESTLON  = -20.0
      EASTLON  = 10.0

! Balanced surface pressure & pot. temp. increments for wind ?
#if defined(GLOBAL)
      LWBAL_SF = .TRUE.
      LWBAL_UA = .TRUE.

#else
      LWBAL_SF = .FALSE.
      LWBAL_UA = .FALSE.

#endif

! Enable calculation of balanced pot. temp. incs. for wind ?
! if true WINDBAL theta incrs are calculated from upper air winds.
      WB_THETA_UA = .TRUE.

! if true WINDBAL theta incrs are calculated from surface winds
! (using HYDRST).
      WB_THETA_SF = .TRUE.

! Initialize variables to define the correlation coefficients for
! WINDBAL.
#if defined(GLOBAL)
      WB_CC_N       = 0.70
      WB_CC_S       = 0.80
#else
      WB_CC_LAM     = 0.70

! Initialize variables to define a subset of limited area model on
! which a multigrid can be run efficiently (for use in WINDBAL).
! The LAM values are given (with mesoscale values after the comment).
! N.B. the LAM is 229 E/W by 132 N/S, and the mesoscale model is 92 by
! 92.

#endif
      WB_START_LEV   = BL_LEVELS - 1   !first non zero level
      WB_END_LEV     = P_LEVELS - 2    !last  non zero level
      WB_LAND_SCALE  = .TRUE.  ! rescale WINDBAL incs over land ?
      WB_LAND_FACTOR = 0.5     ! factor for rescaling   "   "

! Synoptic insertion mode
      LSYN = .FALSE.

! Vertical filtering parameter (used in vrtf)
      VERT_FILT = 0.0

! vertical correlation function coeff for aerosol
      CSCALE_VERT_AERO = 18.0

! Vertical correlation function coeff
! EXP(-(CSCALE_VERT*(LN(P1)-LN(P2)))**2)
! VERT_COR_SCALE is interpolated from  CSCALE_VERT
! which is defined for a regular set of 21 eta levels
! (1 through 0 step 0.05)
! 4 copies for extratropical non-wind (,1);tropical non-wind (,2)
!          for extratropical     wind (,3);tropical     wind (,4)
      IF (LAC_UARS) THEN
       ! Altered to decay more sharply than before
       ! (b in equation 5.1.3 is 9 rather than 3 (old default)
       ! or 16 in tropics - same at all levels)
        DO J = 1, 21
          CSCALE_VERT(J,1) = 3.0
          CSCALE_VERT(J,2) = 4.0
          CSCALE_VERT(J,3) = 3.0
          CSCALE_VERT(J,4) = 4.0
        END DO
      ELSE IF (LAC_MES) THEN
       ! use revised mes values
       DO JJ = 1, 4
          DO J = 1, 21
            CSCALE_VERT(J,JJ) = DEF_CSCALE_VERT_MES(J,JJ)

          END DO
        END DO
!     ELSE for LAM and GL (oper) use DATA table
      END IF

! Cut off in vertical correlation function
! single lvl obs(203,303,403) infl only vert_cutoff_sl scale hts
! same cut-off correl. applies to extrap'n of incomplete soundings.
      IF (LAC_UARS) THEN
        VERT_CUTOFF_SL = 1.0

      ELSE
        VERT_CUTOFF_SL = 2.0

      END IF

! Cut off in vertical correlation function for bogus wind data
      IF (LAC_UARS) THEN
        VERT_CUTOFF_BW = 1.0

      ELSE
        VERT_CUTOFF_BW = 1.0

      END IF

! Cut off in vertical correlation function for bogus humidity data
      VERT_CUTOFF_BH = 0.10

! Correlation function option
!     MHCORFN = 3 For auto-regressive function
!             = 4 with mod which tends to zero at edge on infl
!               +16 For time factor option in INITAC & AC
      IF (LAC_MES) THEN
        MHCORFN = 3

      ELSE
        MHCORFN = 4

      END IF

! Set ac diagnostic mode (each timestep)
! +8 FOR DIAG ONLY ITERATION AT END OF ITERATIONS
! +32 FOR DIAGS ON (OVERRULING LDIAGAC)
! +64 FOR DIAG ONLY ITERATION BETWEEN ITERATIONS
      DO J = 1, MODEACP
        MACDIAG(J) = 0

      END DO

      IF (LAC_MES) THEN
       ! Set MACDIAG to get diagnostics on timesteps 1, 12, 24, 36.
        MACDIAG(1)   = 32
        MACDIAG(12)  = 32
        MACDIAG(24)  = 32
        MACDIAG(36)  = 32

      ELSE  !  global and lam
       ! Set MACDIAG to get diagnostics on timesteps 1, 18 and 36.
        MACDIAG(1)  = 32
        MACDIAG(18) = 32
        MACDIAG(36) = 32

      END IF

! set mode for weights formula
      MWTFN = 2

! set mode for time ramp (1=linear,2=quadratic)
      IF (LAC_UARS) THEN
        MRAMPFN = 1

      ELSE
        MRAMPFN = 2

      END IF

! set mode for data density formula
      MDATADFN = 2

! set mode for GLOSS (type 208) processing option
! 1 as 205,2 use GLOSS error ratio,
! 3 do constrained retrieval, 4 as 3 without a background correction
      MGLOSSFN = 1

!  Min speed limit for ERS-1 winds above which direction is believed
      SPEED_LIMIT305 = 4.0

! set parameters to define LASS/GLOSS vertical interpolation
      MVINT205 = 2     ! 1 for V_INT_T, 2 for V_INT_TP

! set mode for latiude weighting of geostrophic increments
      MGEOWT = 3

!     VERTICAL SCALING FACTOR FOR GEOSTROPHIC INCREMENTS
!     this is defined for a regular set of 21 eta levels
!       (1 through 0 step 0.05)
      IF (LAC_UARS) THEN
! UARS defaults give full weights, except near jet level
        DO J=1,21
          GEOWT_VERT(J) = 1.0
        END DO

        GEOWT_VERT(13)=0.8
        GEOWT_VERT(14)=0.6
        GEOWT_VERT(15)=0.4
        GEOWT_VERT(16)=0.3
        GEOWT_VERT(17)=0.4
        GEOWT_VERT(18)=0.6
        GEOWT_VERT(19)=0.8
      ELSE
! the defaults here give full wt below 500mb and zero wt above 250mb
        DO J=1,10
          GEOWT_VERT(J) = 1.0
        END DO

        GEOWT_VERT(11) = 0.820
        GEOWT_VERT(12) = 0.720
        GEOWT_VERT(13) = 0.620
        GEOWT_VERT(14) = 0.500
        GEOWT_VERT(15) = 0.367
        GEOWT_VERT(16) = 0.200

        DO J = 17, 21
          GEOWT_VERT(J) = 0.0
        END DO
      END IF

! Vertical scaling for horizontal correlation scale
      IF (LAC_MES) THEN
        DO J = 1, BL_LEVELS
          SCFACT_VERT(J)   = 0.875
          NSLABS_SCFACT(J) = 0
          CSCFACT_V(J)     = 0.0
        END DO
        DO J = BL_LEVELS+1, P_LEVELS_MAX
          SCFACT_VERT(J)   = 1.0
          NSLABS_SCFACT(J) = 0
          CSCFACT_V(J)     = 0.0
        END DO

      ELSE
        DO J = 1, P_LEVELS_MAX
          SCFACT_VERT(J)   = 1.0
          NSLABS_SCFACT(J) = 0
          CSCFACT_V(J)     = 0.0
        END DO

      ENDIF

! Defaults for Filtered Increment Mode of Horizontal Analysis
! For testing, FI_SCALE_FACTOR=1.0 will match current HORINF
! Eventually, higher levels can be given a larger scale.
! DF_COEFF=1.0 is initial estimate. DF_COEFF=0.0 turns off DIVFILT
! NPASS_RF=2 matches SOAR and is default.

      NPASS_RF = 2
      IF (LAC_MES) THEN
!       For MOPS data, use FI method without filtering on model grid.
!       (Assumes that no other mes data groups require standard FI)
        NPASS_RF = 0
      END IF

      FI_SCALE = 400000.0
      DF_SCALE = 400000.0

      DO JLEV = 1, P_LEVELS_MAX
        FI_SCALE_FACTOR(JLEV) = 1.0
        DF_COEFF(JLEV)        = 1.0
        DF_SCALE_LEV(JLEV)    = DF_SCALE * FI_SCALE_FACTOR(JLEV)

      END DO

! Factor for non-divergence correction term in correlations
      NON_DIV_COR = 0.8
!  same factor for 10m wind data
      IF(LAC_MES) THEN
        NON_DIV_COR_10M = 0.0
      ELSE
        NON_DIV_COR_10M = 0.8
      ENDIF

! Control of diagnostic output from RDOBS/RDOBS2
! 0 - No listing ; 1 - Standard Listing ; 2 - More detailed listing.
      DIAG_RDOBS = 1

! Thresholds (mm/hr) for rainfall verification
      THRESH_DL   = 0.03
      THRESH_LM   = 0.125
      THRESH_MH   = 0.5
      THRESH_RMSF = 0.125

! Radar data reliability ranges.
! RADAR_RANGE is limit of high reliability (used in LHN & HCS)
! RADAR_RANGE_MAX is limit of usability (used in LHN)
      RADAR_RANGE = 100.0
!  Max Radar range
      RADAR_RANGE_MAX = 200.0
!  --- SPECIFY RADAR LOCATIONS ------
      RADAR_LAT(1)   =  58.21
      RADAR_LON(1)   =  353.81      !  BEACON HILL

      RADAR_LAT(2)   =  57.43
      RADAR_LON(2)   =  357.97      !  HILL OF DUDWICK

      RADAR_LAT(3)   =  55.69
      RADAR_LON(3)   =  355.77      !  CORSE HILL

      RADAR_LAT(4)   =  54.50
      RADAR_LON(4)   =  353.65      !  CASTOR BAY

      RADAR_LAT(5)   =  53.43
      RADAR_LON(5)   =  353.75      !  DUBLIN

      RADAR_LAT(6)   =  52.69
      RADAR_LON(6)   =  351.07      !  SHANNON

      RADAR_LAT(7)   =  53.75
      RADAR_LON(7)   =  357.72      !  HAMELDON

      RADAR_LAT(8)   =  53.33
      RADAR_LON(8)   =  359.45      !  INGHAM

      RADAR_LAT(9)   =  52.40
      RADAR_LON(9)   =  357.40      !  CLEE HILL

      RADAR_LAT(10)  =  51.69
      RADAR_LON(10)  =  359.47      !  CHENIES

      RADAR_LAT(11)  =  51.98
      RADAR_LON(11)  =  355.55      !  CRUGYGORLLWYN

      RADAR_LAT(12)  =  50.96
      RADAR_LON(12)  =  356.55      !  COBBACOMBE

      RADAR_LAT(13)  =  50.82
      RADAR_LON(13)  =  357.44      !  WARDON HILL

      RADAR_LAT(14)  =  50.00
      RADAR_LON(14)  =  354.77      !  PREDANNACK

      RADAR_LAT(15)  =  49.18
      RADAR_LON(15)  =  357.78      !  JERSEY

! List of which radars to use
! (.TRUE. implies MOPS precip data close to that radar will be used
! in rainfall verification and hydrology correction)
      LRADAR(1)   = .TRUE.     !  Beacon Hill
      LRADAR(2)   = .TRUE.     !  Hill of Dudwick
      LRADAR(3)   = .TRUE.     !  Corse Hill
      LRADAR(4)   = .TRUE.     !  Castor Bay
      LRADAR(5)   = .TRUE.     !  Dublin
      LRADAR(6)   = .TRUE.     !  Shannon
      LRADAR(7)   = .TRUE.     !  Hameldon
      LRADAR(8)   = .TRUE.     !  Ingham
      LRADAR(9)   = .TRUE.     !  Clee
      LRADAR(10)  = .TRUE.     !  Chenies
      LRADAR(11)  = .TRUE.     !  Crugygorllwyn
      LRADAR(12)  = .TRUE.     !  Cobbacombe
      LRADAR(13)  = .TRUE.     !  Wardon Hill
      LRADAR(14)  = .TRUE.     !  Predannack
      LRADAR(15)  = .TRUE.     !  Jersey


!  506 observation error calculations (WP171, eqn10)
      L_506_OBERR = .FALSE.
      F1_506 = 10.0
      F2_506 = 1.0
      F3_506 = 0.0
!  Scale points by 1/EPSILON, in LHN, if no near rain found
      L_LHN_SCALE = .TRUE.
!  Use search routine for nearby rain in LHN
      L_LHN_SEARCH = .TRUE.
!  Display detailed diagnostics for LHN routine
      LHN_DIAG = .TRUE.
!  Max range (in grid points) to search to in Latent Heat Nudging code
      LHN_RANGE = 4
!  Verify up to reliable range, or max range
      L_VERIF_RANGE = .TRUE.
!  Minimum acceptable ratio of model ppn to observed, for LHN
      EPSILON_LHN = 0.333
!  Relaxation coefficient for theta increments from the LHN scheme
      RELAX_CF_LHN = 1.0
!  Minimum factor by which model ppn can be scaled to fit obs ppn
      ALPHA_LHN = 0.333
!  Set default for switch to limit size of increments
      L_LHN_LIMIT = .FALSE.
!  Maximum size of increment to theta within LHN (K/day)
      LHN_LIMIT = 1.0
!  Set default for switch to use limit in ALPHA_LHN
      L_LHN_FACT = .TRUE.
!  Set default for switch to filter LHN theta incrs
      L_LHN_FILT = .TRUE.
!  Recursive filter scale in metres
      FI_SCALE_LHN = 17000.0
!  Number of passes through recursive filter, for LHN increments
      NPASS_RF_LHN = 2
! List of observations to be omitted from assimilation
! ** Specify model observation type to omit obs **
      DO JTYP = 1, NANALTYP
        IOMITOBS(JTYP) = 0

      END DO

! Diagnostics options for programmers
      NPROG = 0

! Set timing switch
      LTIMER_AC = .FALSE.

! Set control for calling CHECK_OBS
      LCHECK_GRID = .FALSE.

! Set REMOVE_NEG_LH to false
      REMOVE_NEG_LH =  .FALSE.

! Set USE_CONV_IN_MOPS to true
      USE_CONV_IN_MOPS = .TRUE.


!L 1.4 Read in Namelist set by user
      READ (5, ACP, END=140, ERR=140)

! QA Fortran recommends use of IOSTAT but it doesn't work with namelist
! Jump over NAMELIST read error trap:
      GOTO 141

! Trap NAMELIST read error
140   CONTINUE
        ICODE    = 1
        CMESSAGE = ' ACP_NAMEL : Error reading 2nd ACP namelist'
        if(mype == 0) PRINT *, CMESSAGE

       ! Jump to end of routine
        GOTO 999

! No NAMELIST read error so continue...
141   CONTINUE

!L 1.5 Print out Namelist for this run
      if(mype == 0)then
       PRINT *, ' '
       PRINT *, ' Parameters in Namelist ACP for this run'
       WRITE(6,ACP)
      endif

!L 2   Check Validity of Namelist Parameters
#if defined(MPP)
! Check AGLATDEC
      IF(AGLATDEC <  89.)THEN
          AGLATDEC = 90.
          if(mype == 0)                                                 &
     &     PRINT *, ' AGLATDEC reset to 90. when running MPP '
      END IF
#endif
! Check that AC_OBS_TYPES has any obs types
      NCOUNT = 0

      DO JOBT = 1, NOBTYPMX
        IF (AC_OBS_TYPES(JOBT)  >   0) THEN
          NCOUNT = NCOUNT + 1

        END IF
      END DO

      IF (NCOUNT  ==  0) THEN
        ICODE    = 1
        CMESSAGE = ' ACP_NAMEL : No obs types in AC_OBS_TYPES.'
        GOTO 999

      END IF

! Convert any old type numbers to new type numbers.
      DO JOBT = 1, NOBTYPMX
        IF (AC_OBS_TYPES(JOBT)  ==  501) THEN
          AC_OBS_TYPES(JOBT) = 302
          if(mype == 0)                                                 &
     &     PRINT *, 'Type 501 in AC_OBS_TYPES changed to Type 302'

        END IF

        IF (AC_OBS_TYPES(JOBT) == 502) THEN
          AC_OBS_TYPES(JOBT) = 305
          if(mype == 0)                                                 &
     &     PRINT *, 'Type 502 in AC_OBS_TYPES changed to Type 305'

        END IF
      END DO

#if defined(GLOBAL)
! Check that troplat and tropint are correctly specified
      IF (TROPLAT + TROPINT  >   90.0) THEN
        ICODE    = 1
        CMESSAGE =                                                      &
     &          ' ACP_NL : TROPLAT/TROPINT bug, is TROPLAT-TROPINT >= 0'
        GOTO 999

      END IF

! Check geowt_nh non zero if geostrophic increments are to be used
      IF (GEOWT_NH  ==  0.0 .AND. (LGEO)) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : GEOWT_NH = 0 not allowed with LGEO = T'
        GOTO 999

      END IF

#else
! Check geowt_lam non zero if geostrophic increments are to be used
      IF (GEOWT_LAM  ==  0.0 .AND. (LGEO)) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : GEOWT_LAM=0 not allowed with LGEO = T'
        GOTO 999

      END IF
#endif

! Check windbal variables
#if defined(GLOBAL)
      IF (WB_CC_N  <   0.0 .OR. WB_CC_N  >   1.0) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : invalid WB_CC_N'
        GOTO 999

      ENDIF
      IF (WB_CC_S  <   0.0 .OR. WB_CC_S  >   1.0) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : invalid WB_CC_S'
        GOTO 999

      ENDIF
#else
      IF (WB_CC_LAM  <   0.0 .OR. WB_CC_LAM  >   1.0) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : invalid WB_CC_LAM'
        GOTO 999

      ENDIF
#endif
      IF (WB_START_LEV  <   1 .OR. WB_START_LEV  >   P_LEVELS-5) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : invalid WB_START_LEV'
        GOTO 999

      ENDIF
      IF (WB_END_LEV    <   5 .OR. WB_END_LEV    >   P_LEVELS-1 .OR.    &
     &    WB_END_LEV-WB_START_LEV  <   4) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : invalid WB_END_LEV'
        GOTO 999

      ENDIF

! Check MRAMPFN
      IF (MRAMPFN  /=  1 .AND. MRAMPFN  /=  2) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : invalid MRAMPFN'
        GOTO 999

      ENDIF

! Check MDATADFN
      IF (MDATADFN  /=  1 .AND. MDATADFN  /=  2) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : invalid MDATADFN'
        GOTO 999

      ENDIF

! Check MGLOSSFN
      IF (MGLOSSFN  <   1 .OR. MGLOSSFN  >   4) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : invalid MGLOSSFN'
        GOTO 999

      ENDIF

! Check MWTFN
      IF (MWTFN /= 1 .AND. MWTFN /= 2) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : invalid MWTFN'
        GOTO 999

      ENDIF

! Check MHCORFN
      IF (MOD(MHCORFN,8) <  1 .OR. MOD(MHCORFN,8) >  4) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : invalid MHCORFN'
        GOTO 999

      END IF

! Check MVINT205
      IF (MVINT205  /=  1 .AND. MVINT205  /=  2) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : invalid MVINT205'
        GOTO 999

      ENDIF

! Check NO_OBS_FILES
      IF (NO_OBS_FILES >  10) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : invalid NO_OBS_FILES (max allowed is 10)'
        GOTO 999
      END IF

! Check OBS_FORMAT - only allow portable UM format
      IF (OBS_FORMAT /= 2 .AND. OBS_FORMAT /= 3) THEN
       ICODE    = 1
       CMESSAGE = 'ACP_NL : OBS_FORMAT not equal to 2 or 3'
       GOTO 999
      END IF

!L 3. Set Defaults for Group Dependent arrays
! DEPENDS ON: def_group
      CALL DEF_GROUP (P_LEVELS, Q_LEVELS, BL_LEVELS, TR_LEVELS,         &
     & ICODE, CMESSAGE)

      IF (ICODE >  0) GO TO 999

!L 4. Set Defaults for Obs Type Dependent arrays
! DEPENDS ON: def_type
      CALL DEF_TYPE (ICODE,CMESSAGE)

      IF (ICODE >  0) GO TO 999

!L 5. Check/Process Group dependent variables
! DEPENDS ON: group_dep_var
      CALL GROUP_DEP_VAR (AC_ORDER, NO_ITERATIONS, INTERVAL_ITER,       &
     &     N_ANAL_LEVS, N_WT_LEVS,                                      &

#if defined(GLOBAL)
     &     NUDGE_NH, NUDGE_TR, NUDGE_SH,                                &

#else
     &     NUDGE_LAM,                                                   &

#endif
     &     AGRES_ROWS, AGRES_PTS, MODE_HANAL, FI_VAR_FACTOR,            &
     &     ICODE,CMESSAGE)

      IF (ICODE >  0) GO TO 999

!L 6. Check/Process Obs Type dependent variables
! DEPENDS ON: type_dep_var
      CALL TYPE_DEP_VAR (TIMEB,TIMEA,TGETOBB,TGETOBA,RADINF,OBTHIN,     &
     &                   CSCALE_START,CSCALE_OBTIME,CSCALE_END,         &
     &                   ICODE,CMESSAGE)

      IF (ICODE >  0) GO TO 999


!L 11. Convert rainfall thresholds from mm/hr to kgm-2s-1
       THRESH_DL   = THRESH_DL   / 3600
       THRESH_LM   = THRESH_LM   / 3600
       THRESH_MH   = THRESH_MH   / 3600
       THRESH_RMSF = THRESH_RMSF / 3600

!L  12. Check Lat/Lon verification coordinates

! Need longitudes in range :  -180 < LON <= 180
      IF (EASTLON  >   180.0) THEN
        EASTLON = EASTLON - 360.0
        if(mype == 0) PRINT *,                                          &
     &    "EASTLON changed from ",EASTLON+360.0," to ",EASTLON
      ENDIF

      IF (WESTLON  >   180.0) THEN
        WESTLON = WESTLON - 360.0
        if(mype == 0) PRINT *,                                          &
     &    "WESTLON changed from ",WESTLON+360.0," to ",WESTLON
      ENDIF

      IF (NORTHLAT  <   SOUTHLAT) THEN
        if(mype == 0) PRINT *,                                          &
     &    "Swapping NORTHLAT and SOUTHLAT so that North>South"
        TEMP_COORD = NORTHLAT
        NORTHLAT   = SOUTHLAT
        SOUTHLAT   = TEMP_COORD
      ENDIF

      IF (EASTLON  <   WESTLON) THEN
        if(mype == 0) PRINT *,                                          &
     &   "Swapping EASTLON and WESTLON so that East > West"
        TEMP_COORD = EASTLON
        EASTLON    = WESTLON
        WESTLON    = TEMP_COORD
      ENDIF

 999  CONTINUE
      RETURN
      END SUBROUTINE ACP_NAMEL
!*L  Arguments:---------------------------------------------------
#endif
