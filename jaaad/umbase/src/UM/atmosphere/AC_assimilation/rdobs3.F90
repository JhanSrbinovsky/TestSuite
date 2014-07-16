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
      SUBROUTINE RDOBS3 (UNIT_NO,NOBTYPMX,NOBTYP,OBSTYP,NOBS,NDATAV,    &
     &                   NOBLEV,OBLEVTYP,OBS_LEVELS,OBS_DATA,           &
#include "dump_ar2.h"
     &                   LEN_DATA,MAX_NDV,TNOBS,MISSD,P_LEVELS,         &
#include "argppx.h"
     &                   ICODE,CMESSAGE                                 &
#if defined(MPP)
     &                         ,IPT                                     &
#endif
     &                             )

      IMPLICIT NONE

!L----------------------------------------------------------------------
      INTEGER IPT
      INTEGER                                                           &
     &   UNIT_NO                                                        &
                          ! IN  Unit no of observation file
     &  ,NOBTYPMX                                                       &
                          ! IN  Max no of observation types
     &  ,NOBTYP                                                         &
                          ! OUT No of observation types
     &  ,OBSTYP(NOBTYPMX)                                               &
                          ! OUT Observation types
     &  ,NOBS(NOBTYPMX)                                                 &
                          ! OUT No of Observations
     &  ,NDATAV(NOBTYPMX)                                               &
                          ! OUT No of data values
     &  ,NOBLEV(NOBTYPMX)                                               &
                          ! OUT No of observation levels
     &  ,OBLEVTYP(NOBTYPMX)                                             &
                           !OUT Observation level type
     &  ,LEN_DATA                                                       &
                          ! IN  Dimension of data section
     &  ,TNOBS                                                          &
                          ! OUT Total no of observations
     &  ,P_LEVELS                                                       &
                          ! IN  No of model levels
     &  ,ICODE            ! OUT Return Code

      REAL                                                              &
     &   MISSD                                                          &
                                   ! OUT Real missing data indicator
     &  ,OBS_LEVELS(P_LEVELS+1,*)                                       &
                                   ! OUT Observation levels
     &  ,OBS_DATA(LEN_DATA)        ! OUT Observation data

!DR   REAL DATALEVS(LEN1_LEVDEPC-2,*)     ! OUT Obs levels

      CHARACTER*(*) CMESSAGE ! OUT Error message if ICODE > 0
!L----------------------------------------------------------------------
!L
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "dump_len.h"
!     Dynamic allocate arrays for READDUMP
#include "dump_dim.h"

!     Do not autotask this routine
!FPP$ NOCONCUR R
!
! EXTERNAL SUBROUTINE CALLS
      EXTERNAL READACOBS
!     LOCAL VARIABLES
      INTEGER JOBT,JLEV
      INTEGER MAX_NDV   !  Maximum no of data values (for obs type)
!-----------------------------------------------------------------------

!     Go to start of obs file
! DEPENDS ON: setpos
      CALL SETPOS (UNIT_NO,0,ICODE)

!     Read in the observation file
! DEPENDS ON: readacobs
      CALL READACOBS (UNIT_NO,                                          &
     &  FIXHD,    LEN_FIXHD,                                            &
     &  INTHD,    LEN_INTHD,                                            &
     &  REALHD,   LEN_REALHD,                                           &
     &  LEVDEPC,  LEN1_LEVDEPC,LEN2_LEVDEPC,                            &
     &  ROWDEPC,  LEN1_ROWDEPC,LEN2_ROWDEPC,                            &
     &  COLDEPC,  LEN1_COLDEPC,LEN2_COLDEPC,                            &
     &  FLDDEPC,  LEN1_FLDDEPC,LEN2_FLDDEPC,                            &
     &  EXTCNST,  LEN_EXTCNST,                                          &
     &  DUMPHIST, LEN_DUMPHIST,                                         &
     &  CFI1,     LEN_CFI1,                                             &
     &  CFI2,     LEN_CFI2,                                             &
     &  CFI3,     LEN_CFI3,                                             &
     &  LOOKUP,                                                         &
     &  LEN1_LOOKUP_OBS,LEN2_LOOKUP_OBS,                                &
     &               LEN_DATA,OBS_DATA,                                 &
#include "argppx.h"
     &               ICODE,CMESSAGE                                     &
     &                         ,IPT                                     &
     &                             )

      IF (ICODE >  0) GOTO 9999

      TNOBS  = INTHD(28)   !  Total no of observations
      NOBTYP = INTHD(32)   !  No of observation types
      MISSD  = REALHD(29)  !  Real MDI used in observations

      IF (NOBTYP >  0) THEN

        DO JOBT = 1,NOBTYP
          OBSTYP(JOBT) = LOOKUP(65,JOBT)   ! Observation Types
          NOBS  (JOBT) = LOOKUP(66,JOBT)   ! No of obs
          NDATAV(JOBT) = LOOKUP(67,JOBT)   ! No of data values
          NOBLEV(JOBT) = LEVDEPC(2,JOBT)   ! No of obs levels
          OBLEVTYP(JOBT) = LEVDEPC(1,JOBT) ! Level type
          DO JLEV =1,LEN1_LEVDEPC-2
            OBS_LEVELS(JLEV,JOBT) = LEVDEPC(2+JLEV,JOBT) ! Obs levels
          ENDDO
        ENDDO

!       MAXNLEV1 = LEN1_LEVDEPC-2  ! Max no of levels + 1

      ENDIF
9999  CONTINUE
      RETURN
      END SUBROUTINE RDOBS3
#endif
