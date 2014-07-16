#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  4 Subroutines in deck : RDOBS, RDOBS2, RDOBS3 and DAYS -----
!LL
!LL  Purpose : Read from ACOBS Files,reformat and place OBS header
!LL            details in COMOBS. The bulk of the required OBS data
!LL            is put into dynamic work array OBS for transmission via
!LL            argument list to GETOBS. OBS is written out to a cache
!LL            file for subsequent reading at later timesteps.
!LL            Thus reread of ACOBS files only required intermittently
!LL            (The routine DAYS does a dd/mm/yy to dayno)
!LL
!LL  For use on Cray
!LL  For Cray - Global  ; Enable defs GLOBAL
!LL
!LL S.Bell      <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL   3.2   25/05/93  Call SETPOS to get to start of new format obs
!LL                   files. Call READ_FLH to get FIXHD(162). Correct
!LL                   argument list and dimensions in RDOBS3. Tidy up
!LL                   print out of observation levels. D Robinson.
!LL           8/7/93  Eliminate QA FORTRAN complaints   S Bell
!LL
!LL   3.2   10/07/93  Cater for type 306    S Bell
!LL   3.3   18/04/94  Modify arglist to READDUMP.    T Johns
!LL   3.3   25/11/93  Correct setting of headers in RDOBS2 Greg Bason
!LL   4.0   20/07/95  Remove references to old format acobs GBason
!LL
!LL   4.0   02/11/95  Remove redundant code. (N.Farnon)
!     4.1   18/06/96  Changes to cope with changes in STASH addressing
!                     Author D.M. Goddard.

!LL  4.1   31/05/96     The number of v points to be processed on a
!LL                     C grid differs from u by row_length. u,v
!LL                     dimensioned separately in call to WLLTOEQ.
!LL                     Requirement for VAR.
!LL                     Author I.Edmond       Reviewer D. Goddard
!     4.2   25/11/96  T3E mods + adjust time window  Stuart Bell

!LL
!LL  4.1  4/09/96:  Port to CRAY T3E  Deborah Salmond
!LL  4.2   5/9/96     Fix overwriting triggered by 4.1 mod above. S.Bell
!LL  4.3   19/03/97   Removed MPP args from READDUMP     P.Burton
!    4.3  31/1/97: Bugfix for LAM on T3E. Stuart Bell
!    4.4  2/6/97: Fix to use less memory and allow more obs. SB/DS
!    4.4  17/10/97: Fix to ensure obs on 1 PE only. Deborah Salmond
!    5.2  12/11/00: change attop,base,right,left to elements of
!                   at_extremity              B Macpherson
!    5.3  09/07/01: amend for S->N ND grid order
!                     Bruce Macpherson
!    5.3  05/12/01:  Remove reference to the shmcomm & iovarsac include
!                    files, use a local dynamic array rather than
!                    common block for array 'work' in rdobs2.  S.Cusack
!    6.0 11/09/03:   Removed double ? for IBM cpp.             P.Dando
!    6.0  10/10/03:  Replace SHMEM with GCOM for SX6. Clive Jones
!    6.0  30/12/03:  Make argument list to READACOBS consistent with
!                    subroutine. Dave Robinson.
!    6.1  17/08/04:  Amend boundaries for correct allocation of
!                    observations to processors. Adam Maycock,
!                    (J.Bornemann lodged).
!    6.2  21/10/05:  Replace GSYNC with SSYNC. P.Selwood
!    6.2  15/08/05   Free format fixes. P.Selwood
!    6.2  24/01/06:  Replace large hard-wired MAX_SHMEM_SIZE and
!                    obs(num)dim with dynamic allocation. R Barnes
!
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Logical components covered:
!LL
!LL  Project Task : P3
!LL
!LL  External documentation:
!LL
!LLEND------------------------------------------------------------------
      SUBROUTINE RDOBS (NFILES,TIMESTEP_NO,TIMESTEP,OBS,OBS_FLAG,       &
     &                  TNDV,TNOBS,P_LEVELS,Q_LEVELS,TNDVMAX,NOBSMAX,   &
#include "argppx.h"
     &                  ICODE,CMESSAGE)
!L----------------------------------------------------------------------
!L   INTENT IN:
!L      NFILES   : NO OF AC OBSERVATION FILES TO BE READ
!                : SEE ACP NAMELIST PARAMETER NO_OBS_FILES
!L      TIMESTEP_NO : TIMESTEP NO
!L      TIMESTEP : TIMESTEP IN SECONDS
!L      P_LEVELS : NUMBER OF MODEL LEVELS
!L      Q_LEVELS : NUMBER OF WET MODEL LEVELS
!L      TNDVMAX  : MAX SIZE OF OBS ARRAY
!L      NOBSMAX  : MAX NO OF OBS
!L   INTENT OUT:
!L      TNDV     : ACTUAL SIZE OF OBS ARRAY
!L      TNOBS    : ACTUAL NO OF OBS
!L      OBS      : OBS array
!L      OBS_FLAG : do not use flags
!L      ICODE/CMESSAGE: for error processing
!L----------------------------------------------------------------------

!     Do not autotask this routine
!FPP$ NOCONCUR R
!
      IMPLICIT NONE

! EXTERNAL SUBROUTINE CALLS
      EXTERNAL RDOBS2,TIMER
#include "acparm.h"
#include "comobs.h"
#include "comacp.h"
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
!-----------------------------------------------------------------------
!     ARGUMENTS
      INTEGER TNDV,TNOBS
      INTEGER P_LEVELS,Q_LEVELS
      INTEGER TNDVMAX,NOBSMAX
      INTEGER NFILES
      INTEGER TIMESTEP_NO
      REAL    TIMESTEP
      INTEGER ICODE
      CHARACTER*256 CMESSAGE
!-----------------------------------------------------------------------
!     LOCAL VARIABLES
      INTEGER J
      REAL    TIMEREL
!-----------------------------------------------------------------------
!     Dynamic allocation
      INTEGER OBS_FLAG(NOBSMAX)
      REAL    OBS(TNDVMAX)
#if defined(MPP)
#include "parvars.h"
#include "mppac.h"
#endif
!-----------------------------------------------------------------------
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER ('RDOBS   ',3)

!-----------------------------------------------------------------------
!L           SECTION 1: COPY INPUT FILES TO WORK AREA, ETC.
!-----------------------------------------------------------------------

!     Find relative time for this timestep

      TIMEREL = NINT ((TIMESTEP_NO-1)*TIMESTEP) / 60.0
      IF (ABS(TIMENEXT-TIMEREL) <  0.016) TIMEREL = NINT (TIMEREL)

!     Determine whether AC Observation files need
!     to be read in on this timestep.

      IF (TIMEREL <  TIMENEXT) THEN

!       The AC Observation files do not need to be read in.
!       Use the observations in temporary store.

#if !defined(MPP)
          REWIND IUNITNO
          READ (IUNITNO) TNOBS,TNDV
          IF (TNOBS >  0 .AND. TNDV >  0) THEN

            READ (IUNITNO) (OBS_FLAG(J),J=1,TNOBS),                     &
     &                     (OBS(J),     J=1,TNDV)
            IF (DIAG_RDOBS >  1)                                        &
     &      WRITE (6,*) ' OBS_FLAG AND OBS READ IN'

          ENDIF
#endif

      ELSE

!       Read in the AC Observation Files

! DEPENDS ON: rdobs2
        CALL RDOBS2(NFILES,TIMESTEP_NO,OBS,OBS_FLAG,TNDV,               &
     &              TNOBS,TNDVMAX,NOBSMAX,P_LEVELS,Q_LEVELS,TIMEREL,    &
#include "argppx.h"
     &               ICODE,CMESSAGE)
        IF (ICODE >  0) GOTO 999

      ENDIF

999   CONTINUE
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER ('RDOBS   ',4)
      RETURN
      END SUBROUTINE RDOBS
#endif
