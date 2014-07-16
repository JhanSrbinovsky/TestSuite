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
!    6.2  22/08/05  Restructure defs for ACP namelist. P.Selwood
!    6.2  15/08/05: Fixes for free-format. P.Selwood
!    6.2  23/02/05     Add new logical REMOVE_NEG_LH.
!                      Mark Dixon.
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
!*L  Arguments:---------------------------------------------------
      SUBROUTINE ACDIAG_NAMEL (ICODE,CMESSAGE)
      IMPLICIT NONE

!     Do not multi-task this routine.
!FPP$ NOCONCUR R

! Global parameters:
! UM Constant Comdecks:
#include "c_pi.h"

! AC Comdecks:
#include "acparm.h"
#include "comacp.h"
#include "comacdg.h"
#if !defined(GLOBAL)
#include "commg.h"
#endif
#include "parvars.h"
!*
! Exported variables (INTENT=OUT)
      INTEGER       ICODE              ! Non zero for failure

      CHARACTER*256 CMESSAGE           ! Reason for failure

#if !defined(GLOBAL)
      REAL LAT(4),LON(4)
      REAL XTEMP
#endif
      REAL ZLAT,ZLONG
      INTEGER J
#if !defined(MPP)
      INTEGER  mype
      PARAMETER (mype = 0 ) ! always zero in non-MPP code
#endif

! Namelist acdiag
      NAMELIST /ACDIAG/ LDIAGAC,LLDAC,MODACO,MDACO,                     &
     &          DLATN,DLATS,DLONGW,DLONGE,NDACPRT,                      &
     &          MDIAGTPS,LLBAND,MDIAG,                                  &
     &          DAGLAT,DAGLON,LLDAG0,LLDAG,                             &
     &          LNORMF,LRMS,LTEMP,LVERIF
! End of header.

!L 1.1 Set namelist defaults for variables defined by UI

      LDIAGAC=.FALSE.

!L 1.2 Read in Namelist set by UI

! rewind required to ensure first ACDIAG is read
      REWIND 5
      READ (5, ACDIAG, END=120, ERR=120)

! QA Fortran recommends use of IOSTAT but it doesn't work with namelist
      GOTO 121
120   CONTINUE
        ICODE    = 1
        CMESSAGE = ' ACP_NAMEL : Error reading 1st ACDIAG namelist'
        if(mype == 0) PRINT *, CMESSAGE
        GOTO 999
121   CONTINUE

!L 1.3 Set namelist defaults for variables defined by user
      LLDAC(1) = .FALSE.
      LLDAC(2) = .TRUE.
      LLDAC(3) = .FALSE.
      LLDAC(4) = .FALSE.
      LLDAC(5) = .FALSE.

      LLDAG0   = .FALSE.
      LLDAG(1) = .FALSE.
      LLDAG(2) = .FALSE.
      LLDAG(3) = .FALSE.
      LLDAG(4) = .FALSE.
      LLDAG(5) = .FALSE.

      LNORMF = .TRUE.
      LRMS   = .TRUE.
      LTEMP  = .TRUE.
      LVERIF = .FALSE.

      DO J = 1, NDACOP
        MDACO(J) = 0

      END DO
      MODACO = 1

      DO J = 1, NOBTYPMX
        MDIAGTPS(J) = 0

      END DO

      LLBAND  = .FALSE.
      MDIAG   = 0
      NDACPRT = 10

      DLATN   = 55.0
      DLATS   = 50.0
      DLONGW  =-10.0
      DLONGE  =- 5.0

      DAGLAT  = 57.0
      DAGLON  = 340.0

      NDGPRT  = 6

!L 1.4 Read in Namelist set by user
      READ (5, ACDIAG, END=140, ERR=140)

! QA Fortran recommends use of IOSTAT but it doesn't work with namelist
! Jump over NAMELIST read error trap.
      GOTO 141

! Namelist read error trap:
140   CONTINUE
        ICODE    = 1
        CMESSAGE = ' ACP_NAMEL: Error reading 2nd ACDIAG namelist'
        if(mype == 0) PRINT *, CMESSAGE
        GOTO 999

141   CONTINUE

!L 2. Print out Namelist for this run
      if(mype == 0)Then
       PRINT *, ' '
       PRINT *, ' Parameters in Namelist ACDIAG for this run'
       WRITE(6,ACDIAG)
      endif

!L 3. Process the namelist
! Process ACDIAG namelist parameters and set up variables in
! COMACDG for diagnostic control in the assimilation.

! Check MODACO
      IF (MODACO <  1 .OR. MODACO >  3) THEN
        LLDAC(1) = .FALSE.
        if(mype == 0)PRINT *, ' Invalid value for MODACO - ',MODACO,    &
     &                                           ' ; LLDAC(1) set to F'

      END IF

! If 'obs-model' statistics required for verification purposes
! then reset various parameters in ACP and ACDIAG namelist.
      IF (LVERIF) THEN
        LGEO    = .FALSE.   ! Do not calculate geostrophic increments
        LHYDR   = .FALSE.   ! Do not calculate hydrostatic increments
        LHYDROL = .FALSE.   ! Do not calculate hydrology   increments
      L_LHN = .FALSE.       ! Do not calculate LHN increments
        LNORMF  = .FALSE.   ! Set normalisation factors to 1.0
        LWBAL_SF= .FALSE.   ! Do not calculate P* & theta from windincs
        LWBAL_UA= .FALSE.   !   "        "        "         "     "
        LRMS    = .TRUE.    ! Print RMS values.
        LTEMP   = .TRUE.    ! Use Temp Increments

        DO J = 1, NOBTYPMX
          DEF_TGETOBB(J) = 180.0   ! Time Window of Obs - Before
          DEF_TGETOBA(J) = 181.0   ! Time Window of Obs - After
          DEF_OBTHIN(J) = 1        ! Dont reduce data vols by thinning

        END DO
      END IF

#if !defined(GLOBAL)
      LAT(1) = DLATN
      LON(1) = DLONGW
      LAT(2) = LAT(1)
      LON(2) = DLONGE
      LAT(3) = DLATS
      LON(3) = LON(2)
      LAT(4) = LAT(3)
      LON(4) = LON(1)

! Transform corners of diagnostic area to elf co-ordinates
! DEPENDS ON: lltoeq
      CALL LLTOEQ (LAT,LON,LAT,LON,ELFPLAT,ELFPLON,4)

! Adjust area to be rectangular in elf lat/lon
      DLATS  = 0.5* ( LAT(3) + LAT(4) )
      DLONGW = 0.5* ( LON(1) + LON(4) )
      DLATN  = 0.5* ( LAT(1) + LAT(2) )
      DLONGE = 0.5* ( LON(2) + LON(3) )

! Make sure dlatn >  dlats in new elf co-ords
      IF (DLATS  >   DLATN) THEN
        XTEMP = DLATS
        DLATS = DLATN
        DLATN = XTEMP

      END IF

! The search for obs permits dlongw >  dlonge by assuming
! that an area crossing the meridian is intended, so no test
! as for latitudes.

! Transform diagnostic point to elf co-ordinates
! DEPENDS ON: lltoeq
      CALL LLTOEQ (DAGLAT,DAGLON,DAGLAT,DAGLON,ELFPLAT,ELFPLON,1)
#endif

! Convert Latitudes into Co-Latitudes (0-180) (NP-SP)
      DLATS = 90.0 - DLATS
      DLATN = 90.0 - DLATN

! Convert Co-Latitudes to Radians
      DLATS = DLATS * PI_OVER_180
      DLATN = DLATN * PI_OVER_180

! Check Longitudes within Range (0-360 DEGREES)
      IF (DLONGW <  0.0) DLONGW = DLONGW+360.0
      IF (DLONGE <  0.0) DLONGE = DLONGE+360.0

! Convert Longitudes to Radians
      DLONGW = DLONGW * PI_OVER_180
      DLONGE = DLONGE * PI_OVER_180

999   CONTINUE
      RETURN
      END SUBROUTINE ACDIAG_NAMEL
#endif
