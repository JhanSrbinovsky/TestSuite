#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE GETOBS,GETOB2,GETOB3----------------------------------
!LL
!LL  Purpose : Extract appropriate observations from the COMOBS common
!LL            block and the OBS array (passed from RDOBS via argument
!LL            list)
!LL            GETOBS called from AC gets list of obs in time window
!LL            GETOB2 called from AC2 gets lat, long, time and
!LL                   model obs type for an obs type.
!LL            GETOB3 called from VERTANL gets data specific to a type
!LL                   (eg data values and assoc error ratios)
!LL
!LL  For Cray - Global  ; Enable defs GLOBAL
!LL
!LL  Transfer from Cyber by Dave Robinson.
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL  3.1     18/11/92  : Option to thin GLOSS      S Bell
!LL  3.2     14/6/93   : use OBTHIN in place of THINxxx; cater for 306;
!LL                    : Eliminate QA FORTRAN complaints    S Bell
!LL  3.3     2/12/93   : New ob type 506/407   Bruce M/ Nigel R
!LL  3.4     9/09/94   : New ob type 901 Pete Clark
!LL  4.0     11/8/95   : deleting redundant comments relating to
!LL                    : OBS_FORMAT=1  SB
!LL  4.1   31/05/96     The number of v points to be processed on a
!LL                     C grid differs from u by row_length. u,v
!LL                     dimensioned separately in call to WLLTOEQ.
!LL                     Requirement for VAR.
!LL                     Author I.Edmond       Reviewer D. Goddard
!LL  4.2   25/11/96     Mods for T3E + adjust time window Stuart Bell
!LL  5.2   18/12/00     remove refs to daco-related diagnostic routines
!LL                                                B Macpherson
!LL  6.0   10/10/03     Replace SHMEM with GCOM for SX6. Clive JOnes
!LL  6.2   21/10/05     Replace GSYNC with SSYNC. P.Selwood.
!LL
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LLEND------------------------------------------------------------------
      SUBROUTINE GETOBS (KACT,ASSM_TIME,INOBS,OBS_NO,LENOBT,            &
     &                   OBS,OBS_FLAG,                                  &
     &                   TIMESTEP_NO,DG_ONLY,                           &
     &                   INOBSDIM,ICODE,CMESSAGE)
!L
!L    --------------------------------------------------------
!L  CALLED BY SUBROUTINE AC
!L  THIS OBTAINS ALL OBSERVATIONS TOGETHER WITH THEIR POSITION
!L  OF TYPE LACT(KACT) WHICH ARE WITHIN THE INSERTION PERIOD RELATIVE
!L  TO THE ASSIMILATION TIME
!L  (for scatwinds and sat120 a subset may be chosen)
!L
      IMPLICIT NONE

!     UM AC Comdecks
#include "acparm.h"
#include "comobs.h"
#include "comacp.h"
#include "comacdg.h"
#if defined(MPP)
#include "parvars.h"
#include "mppac.h"
      INTEGER Iproc
#endif
!-----------------------------------------------------------------------
      INTEGER KACT      ! IN Index for this obs type (pos in LACT etc)
      INTEGER INOBS           ! IN  Total no of obs
      INTEGER INOBSDIM        ! IN  Total no of obs (for dimensioning)
      INTEGER OBS_NO(INOBSDIM)   ! IN  Pointers to obs to be assimilated
      INTEGER OBS_FLAG(INOBSDIM) ! IN  Observation flags
      INTEGER TIMESTEP_NO     ! IN  Timestep number
      INTEGER LENOBT          ! OUT No of obs to be assimilated

      REAL    OBS(INOBSDIM,*) ! IN  Observation data for this type
      REAL    ASSM_TIME       ! IN  Assm time relative to start

      LOGICAL DG_ONLY         ! IN  Diagnostic only switch
      INTEGER ICODE           ! IN  Return code
      CHARACTER*256 CMESSAGE  ! IN  Error message
!-----------------------------------------------------------------------
!     Local variables
      INTEGER INB       !  No of obs before insertion period
      INTEGER INA       !  No of obs after insertion period
      INTEGER INF       !  No of obs that have been flagged
      INTEGER INT       !  No of obs thinned out (skipped)
      INTEGER JOB       !  Loop counter over no of obs
      INTEGER ISTART    !  First obs
      INTEGER ISKIP     !  Skip ISKIP obs in loop over obs
      INTEGER IP_TIME   !  Pointer to times in observation data
      INTEGER OBTHIN    !  Thinning factor for this ob type
      INTEGER ISTAT     !  for calls to GCOM routines
      REAL    TDIFF     !  Time difference between obs and assm time
      REAL    TGETOBB   !  Insertion period before obs time
      REAL    TGETOBA   !  Insertion period after obs time
!-----------------------------------------------------------------------
      EXTERNAL TIMER
!-----------------------------------------------------------------------
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('GETOBS  ',3)

      INB=0
      INA=0
      INF=0
      INT=0

! SET UP INDEX VECTOR FOR GATHERING OBS
!  ALLOWING FOR OPTION TO THIN OBS :likely use GLOSS/SAT120/ERS-1 DATA
!  BUT NOT ON DIAGNOSTIC ONLY ITERATION
!     Get OBTHIN for this observation type
      OBTHIN = DEF_OBTHIN( TYPE_INDEX(KACT) )
      IF (DG_ONLY) THEN
        ISTART = 1
        ISKIP  = 1
      ELSE
        IF (OBTHIN >  1 .AND. INOBS >  OBTHIN) THEN
          ISTART = 1 + MOD(TIMESTEP_NO,OBTHIN)
          ISKIP  = OBTHIN
        ELSE
          ISTART = 1
          ISKIP  = 1
        ENDIF
      ENDIF

!     Get TGETOBB and TGETOBA for this observation type
      TGETOBB = DEF_TGETOBB( TYPE_INDEX(KACT) )
      TGETOBA = DEF_TGETOBA( TYPE_INDEX(KACT) )

!     Pointer to observation time for this obs type
      IP_TIME = 3

      IF(INOBS >  0)THEN


      DO JOB = ISTART,INOBS,ISKIP

!     Get time difference between obs time and current assm time
      TDIFF = OBS(JOB,IP_TIME)-ASSM_TIME

!     Use obs if : time diff <  tgetobb-0.5
!         and if : time diff > -(tgetoba-0.5)
!         and if : it has not been flagged
!   -0.5 minutes helps make results same on different machines

      IF (TDIFF >= TGETOBB-0.5) THEN    !  Obs yet to be assimilated
        INB=INB+1
      ELSEIF (TDIFF <= -TGETOBA+0.5) THEN !  Obs has been assimilated
        INA=INA+1
      ELSEIF (OBS_FLAG(JOB) /= 0) THEN  !  Obs has been flagged
        INF=INF+1
      ELSE                              !  Assimilate this observation
        LENOBT=LENOBT+1
        OBS_NO(LENOBT) = OBS_NO_ST(KACT)+JOB
      END IF
!
      ENDDO

      ENDIF  !INOBS >  0

      INT=INOBS-(LENOBT+INB+INA+INF)
!
      IF (LDIAGAC) THEN
#if defined(MPP)
        CountA(1)=LENOBT
        CountA(2)=INB
        CountA(3)=INA
        CountA(4)=INF
        CountA(5)=INT
        If(mype == 0)Then
          Do JOB=1,5
           CountC(JOB)=CountA(JOB)
          EndDo
        Endif

         CALL GC_SSYNC(NPROC,ISTAT)

         Do iproc=1,nproc-1
          IF(mype == 0) THEN
! PE0 receives
            CALL GC_IRECV(IPROC,5,IPROC,ISTAT,COUNTB,COUNTA)
            Do JOB=1,5
              CountC(JOB)=CountC(JOB)+CountB(JOB)
             EndDo
           ELSEIF(mype == IPROC) THEN
! other PEs send
             CALL GC_ISEND(IPROC,5,0,ISTAT,COUNTB,COUNTA)
           ENDIF
         EndDo

         CALL GC_SSYNC(NPROC,ISTAT)

         PRINT '(A,I4,A,I6,A,I6,A,I6,A,I6,A,I6,A)',                     &
     &   ' TYPE',LACT(KACT),'  LENOBT',COUNTC(1),'  OMITTED WERE :',    &
     &   COUNTC(2),' (BEFORE) ',COUNTC(3),' (AFTER) ',                  &
     &   COUNTC(4),' (FLAGGED) ',COUNTC(5),' (THINNED) '

        CALL GC_SSYNC(NPROC,ISTAT)
#else
      PRINT '(A,I4,A,I6,A,I6,A,I6,A,I6,A,I6,A)',                        &
     & ' TYPE',LACT(KACT),'  LENOBT',LENOBT,                            &
     & '  OMITTED WERE :',INB,' (BEFORE) ',INA,' (AFTER) ',             &
     &                    INF,' (FLAGGED) ',INT,' (THINNED) '
#endif
      ENDIF
!
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('GETOBS  ',4)
      RETURN
      END SUBROUTINE GETOBS
!
#endif
