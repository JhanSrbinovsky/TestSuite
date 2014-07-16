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
!
      SUBROUTINE GETOB3 (KACT,OBS,OBDATA,                               &
#if !defined(GLOBAL)
     &                   OBS_LAT,OBS_LONG,                              &
#endif
     &                   LENOBT,INOBS,NDV,OBS_NO,ICODE,CMESSAGE)

      IMPLICIT NONE
#include "acparm.h"
#include "comobs.h"
#include "comacp.h"
#include "comacdg.h"
#if !defined(GLOBAL)
#include "commg.h"
#include "c_pi.h"
#endif
!-----------------------------------------------------------------------
      INTEGER KACT                   ! IN  Index to obs type
      INTEGER LENOBT                 ! IN  No of obs to be assimilated
      INTEGER INOBS                  ! IN  No of obs for this type
      INTEGER NDV                    ! IN  No of data values
!                                    !     excluding header section
      REAL OBS (INOBS,*)             ! IN  OBS array from RDOBS
      REAL OBDATA (LENOBT,NDV)       ! OUT ob values and error ratios
#if !defined(GLOBAL)
      REAL OBS_LAT (LENOBT)          ! IN  latitudes
      REAL OBS_LONG(LENOBT)          ! IN  longitudes
#endif
      INTEGER OBS_NO(LENOBT)         ! IN  pointers to obs.
      INTEGER ICODE                  ! OUT error code and message
      CHARACTER*256 CMESSAGE
!-----------------------------------------------------------------------
#if !defined(GLOBAL)
!L local work arrays
      REAL WLT    (LENOBT)
      REAL WLN    (LENOBT)
      REAL WLN2   (LENOBT)
      REAL COEFF1 (LENOBT)
      REAL COEFF2 (LENOBT)
      REAL UWK    (LENOBT)
      REAL VWK    (LENOBT)
!-----------------------------------------------------------------------
#endif
      REAL ZMINEP
      PARAMETER(ZMINEP=1.E-10)
!     ZMINEP IS THE MINIMUM ALLOWED OBSERVATIONAL ERROR TO AVOID /0.
      CHARACTER*10 LABEL
      INTEGER NEROPT,KTYPE,NLEV,JDV,JOB,KDV
!-----------------------------------------------------------------------
      EXTERNAL TIMER
#if !defined(GLOBAL)
      EXTERNAL EQTOLL,W_COEFF,W_EQTOLL
#endif
!-----------------------------------------------------------------------
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('GETOB3  ',3)

      KTYPE  = LACT (KACT)     !  AC Observation Type
      NLEV   = NOBLEV(KACT)    !  No of levels for this type

!     Get observation data and error ratios
      DO JDV=1,NDV
        DO JOB=1,LENOBT
          OBDATA(JOB,JDV) = OBS(OBS_NO(JOB)-OBS_NO_ST(KACT),NDVHDR+JDV)
        ENDDO
      ENDDO

      IF (LDIAGAC) THEN


        IF(NDAC >  0)THEN

#if !defined(GLOBAL)
!      DO LAT/LON TRANSFORMATIONS AND GET COEFF1,COEFF2 OUTSIDE JDV LOOP
#endif

         DO JDV=1,NDV
         IF (KTYPE == 101) THEN
                 IF (JDV == 1) LABEL = 'P* DATA'
                 IF (JDV == 2) LABEL = 'EO/EP'
         ELSEIF (KTYPE == 201) THEN
                 IF (JDV <= NLEV) THEN
                    WRITE(LABEL,21)'T  DATA',JDV
                 ELSE
                    KDV=JDV-NLEV
                    WRITE(LABEL,21)'EO/EP  ',KDV
                 ENDIF
         ELSEIF (KTYPE == 202 .OR. KTYPE == 204) THEN
               IF (JDV == 1) LABEL = 'T DATA'
               IF (JDV == 2) LABEL = 'EO/EP'
         ELSE IF(KTYPE == 203) THEN
               IF (JDV == 1) LABEL = 'OBS LEV'
               IF (JDV == 2) LABEL = 'T DATA'
               IF (JDV == 3) LABEL = 'EO/EP'
         ELSE IF(KTYPE == 205.OR.KTYPE == 206.OR.                       &
     &           KTYPE == 207.OR.KTYPE == 209.OR.                       &
     &           KTYPE == 211) THEN
                 IF (JDV <= NLEV) THEN
                    WRITE(LABEL,21)'T  DATA',JDV
                 ELSE
                    KDV=JDV-NLEV
                    WRITE(LABEL,21)'EO/EP  ',KDV
                 ENDIF
         ELSEIF (KTYPE == 208) THEN
                 IF (JDV <= NLEV) THEN
                    WRITE(LABEL,21)'T  DATA',JDV
                 ELSEIF (JDV >  NLEV.AND.JDV <= 2*NLEV) THEN
                    KDV=JDV-NLEV
                    WRITE(LABEL,21)'T F/G  ',KDV
                 ELSEIF (JDV >  2*NLEV.AND.JDV <= 3*NLEV) THEN
                    KDV=JDV-NLEV*2
                    WRITE(LABEL,21)'EO/EP  ',KDV
                 ELSE
                    LABEL='QUALIND'
                 ENDIF
         ELSEIF (KTYPE == 301 .OR. KTYPE == 311) THEN
                 IF (JDV <= NLEV) THEN
                    WRITE(LABEL,21)'U  DATA',JDV
                 ELSEIF (JDV <= NLEV*2)THEN
                    KDV=JDV-NLEV
                    WRITE(LABEL,21)'V  DATA',KDV
                 ELSE
                    KDV=JDV-NLEV*2
                    WRITE(LABEL,21)'EO/EP  ',KDV
                 ENDIF
         ELSEIF (KTYPE == 302 .OR. KTYPE == 304 .OR.                    &
     &           KTYPE == 305 .OR. KTYPE == 306 ) THEN
                 IF (JDV == 1) LABEL = 'U DATA'
                 IF (JDV == 2) LABEL = 'V DATA'
                 IF (JDV == 3) LABEL = 'EO/EP'
         ELSEIF (KTYPE == 303) THEN
                 IF (JDV == 1) LABEL = 'OBS LEV'
                 IF (JDV == 2) LABEL = 'U DATA'
                 IF (JDV == 3) LABEL = 'V DATA'
                 IF (JDV == 4) LABEL = 'EO/EP'
         ELSEIF (KTYPE == 401 .OR. KTYPE == 405 .OR. KTYPE == 406) THEN
                 IF (JDV <= NLEV) THEN
                    WRITE(LABEL,21)'RH DATA',JDV
                 ELSE
                    KDV=JDV-NLEV
                    WRITE(LABEL,21)'EO/EP  ',KDV
                 ENDIF
         ELSEIF (KTYPE == 402 .OR. KTYPE == 404) THEN
                 IF (JDV == 1) LABEL = 'RH DATA'
                 IF (JDV == 2) LABEL = 'EO/EP'
         ELSEIF (KTYPE == 403) THEN
                 IF (JDV == 1) LABEL = 'OBS LEV'
                 IF (JDV == 2) LABEL = 'RH DATA'
                 IF (JDV == 3) LABEL = 'EO/EP'
         ELSEIF (KTYPE == 407) THEN
                 WRITE(LABEL,21)'CTT CNT',JDV
         ELSEIF (KTYPE == 506) THEN
                 IF (JDV == 1) LABEL = 'PR DATA'
                 IF (JDV == 2) LABEL = 'EO/EP'
         ELSEIF (KTYPE == 901) THEN
                 IF (JDV == 1) LABEL = 'VIS DATA'
                 IF (JDV == 2) LABEL = 'EO/EP'
         ELSE
         ICODE=1
         CMESSAGE=' GETOB3 : ILLEGAL OB TYPE'
         GOTO 999
         ENDIF
!
   21 FORMAT(A7,I3)
!
       ENDDO   !  Loop over JDV
       END IF
      END IF


      NEROPT = NERLEV1(KACT)-NDVHDR

      DO JDV=NEROPT,NEROPT+NLEV-1

!       ENSURE ALL OB ERRORS ARE >= ZMINEP TO AVOID
!       DIVIDING BY ZERO IN VAN## Routines ?? is this necessary

        DO JOB=1,LENOBT
          IF (OBDATA(JOB,JDV) <  ZMINEP .AND. OBDATA(JOB,JDV) /= MISSD) &
     &        OBDATA(JOB,JDV) = ZMINEP
        ENDDO

      ENDDO

999   CONTINUE
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('GETOB3  ',4)
      RETURN
      END SUBROUTINE GETOB3
#endif
