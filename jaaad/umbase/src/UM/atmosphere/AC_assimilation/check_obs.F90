#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE NUM_OBS   NUM_OBS2   CHECK_OBS  ------------------------
!LL
!LL  3 Subroutines in deck : NUM_OBS , NUM_OBS2  and  CHECK_OBS
!LL
!LL  Purpose : Calculate no of observations and no of observation
!LL            values in AC observation files. These values are
!LL            used to dimension any observation data arrays in the
!LL            assimilation code.
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL   3.2   25/05/93  Call SETPOS to get to start of obs files. For
!LL                   new format files only. D. Robinson
!LL           8/7/93  Eliminate QA FORTRAN complaints    S Bell
!LL
!LL   3.3   23/11/93  Include ACOBS/MODEL level and grid consistency
!LL                   checks.  G. Bason.
!LL   4.0   21/2/95  Add control for CHECK_OBS,        G. Bason.
!LL                  and relax check of real constants between file
!LL                  header and model values.          R. Rawlins.
!LL
!LL   4.0   20/07/95  Remove reference to old format acobs GBason
!LL
!LL   4.0   02/11/95  Remove redundant code. (N.Farnon)
!     4.1   18/06/96  Changes to cope with changes in STASH addressing
!                     Author D.M. Goddard.
!    4.2 25/11/96: T3E mods Stuart Bell
!LL
!LL  4.1  4/09/96:  Port to CRAY T3E  Deborah Salmond
!    4.4  3/7/97:  Set PER_FILE_TNDVMAX.   Deborah Salmond
!    4.5  28/7/98: In MPP runs do on 1pe &en broadcast. Deborah Salmond
!    4.5  5/8/98:  increase OBS_FILE_NAME,OBS_DIR_NAME size. Stuart Bell
!    5.2 22/3/01:  amend glsize. Bruce Macpherson
!    5.3 25/9/01:  Portability changes.     Z. Gardner
!    5.3 03/10/01:  comment out call to check_obs
!                   until revised for ND. B Macpherson
!    6.0 10/10/03  : Replace SHMEM with GCOM for SX6 Clive Jones
!    6.1 18/08/04: Declare function LOC as an INTEGER.       P.Dando
!    6.2 21/10/05: Replace GSYNC with SSYNC. P.Selwood
!
!LL  Programming Standard : UM Doc Paper No 3 ; Version 4 ; 5/2/92
!LL
!LL  Project Task : P3
!LL
!LL  Documentation : UM Documentation Paper No ??
!LLEND------------------------------------------------------------------
!LL  SUBROUTINE NUM_OBS2 -----------------------------------------------
!LL
!LL  Written 14/7/92 by Dave Robinson.
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL  Programming Standard : UM Doc Paper No 3 ; Version 4 ; 5/2/92
!LL
!LL  Project Task : P3
!LL
!LL  Purpose : Read in header section of observation file and get
!LL            number of observations and data values in this obs file.
!LL
!LL  Documentation : UM Documentation Paper No ??
!LLEND------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE CHECK_OBS (OBS_ROW_LENGTH,OBS_P_ROWS,OBS_P_LEVELS,     &
     &                      OBS_Q_LEVELS,OBS_AK,OBS_BK,OBS_LONG_RES,    &
     &                      OBS_LAT_RES,OBS_START_LAT,OBS_START_LONG,   &
     &                      OBS_LAT_PSEUDO_POLE,OBS_LONG_PSEUDO_POLE,   &
     &                      P_LEVELS,Q_LEVELS,P_ROWS,ROW_LENGTH,AK,BK,  &
     &                      REALHD1,REALHD2,REALHD3,REALHD4,            &
     &                      REALHD5,REALHD6,ICODE,CMESSAGE)
!-----------------------------------------------------------------------
!     ARGUMENTS
!-----------------------------------------------------------------------
      INTEGER OBS_ROW_LENGTH,OBS_P_ROWS,OBS_P_LEVELS,OBS_Q_LEVELS
      REAL OBS_LONG_RES,OBS_LAT_RES,OBS_START_LAT
      REAL OBS_START_LONG,OBS_LAT_PSEUDO_POLE,OBS_LONG_PSEUDO_POLE
      INTEGER P_LEVELS      !  No of model levels
      INTEGER Q_LEVELS      !  No of model wet levels
      INTEGER P_ROWS        !  No of model rows
      INTEGER ROW_LENGTH    !  No of points on row
      REAL    OBS_AK(P_LEVELS),OBS_BK(P_LEVELS)
      REAL    AK(P_LEVELS)  !  Vertical grid
      REAL    BK(P_LEVELS)
      REAL REALHD1,REALHD2  !  Horizontal grid
      REAL REALHD3,REALHD4
      REAL REALHD5,REALHD6
      INTEGER ICODE           !  Return code
      CHARACTER*256 CMESSAGE  !  Error message if ICODE > 0
#if defined(MPP)
#include "parvars.h"
      INTEGER ROW_LENGTH_GLOBAL,P_ROWS_GLOBAL
#endif
!-----------------------------------------------------------------------
!     LOCAL VARIABLES
!-----------------------------------------------------------------------
      INTEGER JLEV

      REAL P1,P2
      LOGICAL LNER
      LNER(P1,P2) = ((ABS(P1-P2))  >   (1.E-6*ABS(P1+P2)))
!-----------------------------------------------------------------------
#if defined(MPP)
      ROW_LENGTH_GLOBAL=glsize(1,fld_type_p)
      P_ROWS_GLOBAL=glsize(2,fld_type_p)
#endif
      PRINT *, ' '
      PRINT *, 'NUM_OBS : Checking consistency of level and grid'
      PRINT *, '          between ACOBS file and MODEL'
      PRINT *, ' '
      ICODE = 0
#if defined(MPP)
      IF (OBS_ROW_LENGTH /= ROW_LENGTH_GLOBAL) THEN
        PRINT *, 'Inconsistency found:'
        PRINT *, 'row_length_acobs      = ',obs_row_length
        PRINT *, 'row_length_global      = ',row_length_global
        ICODE = 1
      ENDIF
      IF (OBS_P_ROWS /= P_ROWS_GLOBAL) THEN
        PRINT *, 'Inconsistency found:'
        PRINT *, 'p_rows_acobs          = ',obs_p_rows
        PRINT *, 'p_rows_global          = ',p_rows_global
        ICODE = 1
      ENDIF
#else
      IF (OBS_ROW_LENGTH /= ROW_LENGTH) THEN
        PRINT *, 'Inconsistency found:'
        PRINT *, 'row_length_acobs      = ',obs_row_length
        PRINT *, 'row_length_model      = ',row_length
        ICODE = 1
      ENDIF
      IF (OBS_P_ROWS /= P_ROWS) THEN
        PRINT *, 'Inconsistency found:'
        PRINT *, 'p_rows_acobs          = ',obs_p_rows
        PRINT *, 'p_rows_model          = ',p_rows
        ICODE = 1
      ENDIF
#endif
      IF (OBS_P_LEVELS /= P_LEVELS) THEN
        PRINT *, 'Inconsistency found:'
        PRINT *, 'p_levels_acobs        = ',obs_p_levels
        PRINT *, 'p_levels_model        = ',p_levels
        ICODE = 1
      ENDIF
      IF (OBS_Q_LEVELS /= Q_LEVELS) THEN
        PRINT *, 'Inconsistency found:'
        PRINT *, 'q_levels_acobs        = ',obs_q_levels
        PRINT *, 'q_levels_model        = ',q_levels
        ICODE = 1
      ENDIF

      DO JLEV=1,P_LEVELS
        IF ((OBS_AK(JLEV) /= AK(JLEV)) .OR.                             &
     &     (OBS_BK(JLEV) /= BK(JLEV))) THEN
          PRINT *, 'Inconsistency found:'
          PRINT *, 'Level ',JLEV
          PRINT *, 'ak_acobs            = ',obs_ak(jlev)
          PRINT *, 'ak_model            = ',ak(jlev)
          PRINT *, 'bk_acobs            = ',obs_bk(jlev)
          PRINT *, 'bk_model            = ',bk(jlev)
          IF (LNER(OBS_AK(JLEV),AK(JLEV)) .OR.                          &
     &        LNER(OBS_BK(JLEV),BK(JLEV))) THEN
            ICODE=1
          ELSE
            PRINT *, 'But inconsistency is within 32 bit accuracy'
          ENDIF
        ENDIF
      ENDDO

      IF (OBS_LONG_RES /= REALHD1) THEN
        PRINT *, 'Inconsistency found:'
        PRINT *, 'long_res_acobs        = ',obs_long_res
        PRINT *, 'long_res_model        = ',realhd1
       IF (LNER(OBS_LONG_RES,REALHD1)) THEN
          ICODE=1
       ELSE
          PRINT *, 'But inconsistency is within 32 bit accuracy'
       ENDIF
      ENDIF
      IF (OBS_LAT_RES /= REALHD2) THEN
        PRINT *, 'Inconsistency found:'
        PRINT *, 'lat_res_acobs         = ',obs_lat_res
        PRINT *, 'lat_res_model         = ',realhd2
       IF (LNER(OBS_LAT_RES,REALHD2)) THEN
          ICODE=1
       ELSE
          PRINT *, 'But inconsistency is within 32 bit accuracy'
       ENDIF
      ENDIF
      IF (OBS_START_LAT /= REALHD3) THEN
        PRINT *, 'Inconsistency found:'
        PRINT *, 'start_lat_acobs       = ',obs_start_lat
        PRINT *, 'start_lat_model       = ',realhd3
       IF (LNER(OBS_START_LAT,REALHD3)) THEN
          ICODE=1
       ELSE
          PRINT *, 'But inconsistency is within 32 bit accuracy'
       ENDIF
      ENDIF
      IF (OBS_START_LONG /= REALHD4) THEN
        PRINT *, 'Inconsistency found:'
        PRINT *, 'start_long_acobs      = ',obs_start_long
        PRINT *, 'start_long_model      = ',realhd4
       IF (LNER(OBS_START_LONG,REALHD4)) THEN
          ICODE=1
       ELSE
          PRINT *, 'But inconsistency is within 32 bit accuracy'
       ENDIF
      ENDIF
      IF (OBS_LAT_PSEUDO_POLE /= REALHD5) THEN
        PRINT *, 'Inconsistency found:'
        PRINT *, 'lat_pseudo_pole_acobs = ',obs_lat_pseudo_pole
        PRINT *, 'lat_pseudo_pole_model = ',realhd5
      IF (LNER(OBS_LAT_PSEUDO_POLE,REALHD5)) THEN
          ICODE=1
       ELSE
          PRINT *, 'But inconsistency is within 32 bit accuracy'
       ENDIF
      ENDIF
      IF (OBS_LONG_PSEUDO_POLE /= REALHD6) THEN
        PRINT *, 'Inconsistency found:'
        PRINT *, 'long_pseudo_pole_acobs= ',obs_long_pseudo_pole
        PRINT *, 'long_pseudo_pole_model= ',realhd6
      IF (LNER(OBS_LONG_PSEUDO_POLE,REALHD6)) THEN
          ICODE=1
       ELSE
          PRINT *, 'But inconsistency is within 32 bit accuracy'
       ENDIF
      ENDIF

      IF (ICODE == 1) THEN
        CMESSAGE = 'NUM_OBS : Failure in NUM_OBS - mismatch between     &
     &model and acobs level/grid information'
      ELSE
        PRINT *, 'NUM_OBS : ACOBS and MODEL are consistent'
      ENDIF

      RETURN
      END SUBROUTINE CHECK_OBS
#endif
