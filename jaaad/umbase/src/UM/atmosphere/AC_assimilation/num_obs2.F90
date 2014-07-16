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
      SUBROUTINE NUM_OBS2 (UNIT_NO,P_LEVELS,Q_LEVELS,P_ROWS,ROW_LENGTH, &
     &                     AK,BK,REALHD1,REALHD2,REALHD3,REALHD4,       &
     &                     REALHD5,REALHD6,                             &
#include "dump_ar2.h"
     &                     LEN_DATA,NOBTYP,TNOBS,TNDV,                  &
#include "argppx.h"
     &                     ICODE,CMESSAGE)
!L----------------------------------------------------------------------
!
!     Do not autotask this routine
!FPP$ NOCONCUR R
!
!     IMPLICIT NONE

#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
!     AC Comdecks
#include "acparm.h"
#include "comacp.h"
!-----------------------------------------------------------------------
!     ARGUMENTS
!-----------------------------------------------------------------------
      INTEGER UNIT_NO      ! IN  : Unit no of Observation file
      INTEGER P_LEVELS     ! IN  : No of model levels
      INTEGER Q_LEVELS     ! IN  : No of model wet levels
      INTEGER P_ROWS       ! IN  : No of model rows
      INTEGER ROW_LENGTH   ! IN  : No of points on row
      REAL    AK(P_LEVELS) ! IN  : Vertical grid
      REAL    BK(P_LEVELS)
      REAL REALHD1,REALHD2 ! IN  : Horizontal grid
      REAL REALHD3,REALHD4
      REAL REALHD5,REALHD6
      INTEGER LEN_DATA     ! IN  : Dimension of data section
      INTEGER NOBTYP       ! OUT : No of observation types
      INTEGER TNOBS        ! OUT : Total no of observations
      INTEGER TNDV         ! OUT : Total no of data values
      INTEGER ICODE        ! OUT : Return code
      CHARACTER*256 CMESSAGE  !  OUT : Error message if ICODE > 0
!-----------------------------------------------------------------------
!     LEVEL/GRID VARIABLES
!-----------------------------------------------------------------------
      INTEGER OBS_ROW_LENGTH,OBS_P_ROWS,OBS_P_LEVELS,OBS_Q_LEVELS
      REAL OBS_AK(P_LEVELS),OBS_BK(P_LEVELS)
      REAL OBS_LONG_RES,OBS_LAT_RES,OBS_START_LAT
      REAL OBS_START_LONG,OBS_LAT_PSEUDO_POLE,OBS_LONG_PSEUDO_POLE
!-----------------------------------------------------------------------
!     LOCAL VARIABLES
!-----------------------------------------------------------------------
      INTEGER JLEV
!-----------------------------------------------------------------------
#include "dump_len.h"
!-----------------------------------------------------------------------
!     Dynamic allocated arrays
#include "dump_dim.h"
!-----------------------------------------------------------------------
      EXTERNAL READHEAD,SETPOS,FORT_GET_ENV
!-----------------------------------------------------------------------

!     Go to start of obs file
! DEPENDS ON: setpos
      CALL SETPOS (UNIT_NO,0,ICODE)

!     Read in headers from obs file
! DEPENDS ON: readhead
      CALL READHEAD (UNIT_NO,                                           &
#include "dump_ar1.h"
     &               LEN_DATA,                                          &
#include "argppx.h"
     &               START_BLOCK,ICODE,CMESSAGE)
      IF (ICODE >  0) RETURN

      OBS_ROW_LENGTH = INTHD(6)         !  No of points on row
      OBS_P_ROWS     = INTHD(7)         !  No of rows
      OBS_P_LEVELS   = INTHD(8)         !  No of model levels
      OBS_Q_LEVELS   = INTHD(9)         !  No of model wet levels

      TNOBS  = INTHD(28)                !  Total no of observations
      TNDV   = INTHD(29)                !  Total no of data values
      NOBTYP = INTHD(32)                !  No of observation types

      DO JLEV=1,P_LEVELS
        OBS_AK(JLEV) = LEVDEPC(JLEV+2,NOBTYP+1) !  Vertical grid
        OBS_BK(JLEV) = LEVDEPC(JLEV+2,NOBTYP+2) !
      ENDDO

      OBS_LONG_RES         = REALHD(1)  !  Horizontal grid
      OBS_LAT_RES          = REALHD(2)  !
      OBS_START_LAT        = REALHD(3)  !
      OBS_START_LONG       = REALHD(4)  !
      OBS_LAT_PSEUDO_POLE  = REALHD(5)  !
      OBS_LONG_PSEUDO_POLE = REALHD(6)  !

!     Check model and acobs file are consistent
! *** NB CHECK_OBS needs rewriting for ND, so commented out
!      IF(LCHECK_GRID)THEN
!      CALL CHECK_OBS (OBS_ROW_LENGTH,OBS_P_ROWS,OBS_P_LEVELS,
!     +                OBS_Q_LEVELS,OBS_AK,OBS_BK,OBS_LONG_RES,
!     +                OBS_LAT_RES,OBS_START_LAT,OBS_START_LONG,
!     +                OBS_LAT_PSEUDO_POLE,OBS_LONG_PSEUDO_POLE,
!     +                P_LEVELS,Q_LEVELS,P_ROWS,ROW_LENGTH,AK,BK,
!     +                REALHD1,REALHD2,REALHD3,REALHD4,REALHD5,REALHD6,
!     +                ICODE,CMESSAGE)
!      IF (ICODE >  0) RETURN
!      END IF

      RETURN
      END SUBROUTINE NUM_OBS2
!-----------------------------------------------------------------------
#endif
