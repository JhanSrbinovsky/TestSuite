#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE SETDAC -------------------------------------------------
!LL
!LL  Purpose : Controls diagnostic output from AC Scheme.
!LL
!LL            Diagnostic options are set through ADIAG namelist
!LL            which is read in INITAC.
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!    4.2 25/11/96: T3E mods Stuart Bell
!LL
!LL  For use on Cray
!LL  For Cray - Global  ; Enable defs GLOBAL
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LLEND------------------------------------------------------------------
      SUBROUTINE SETDAC(OBS,TNDV)
!L    -----------------

!L  CALLED BY SUBROUTINE RDOBS
!L  THE DIAGNOSTICS ARE CONTROLLED BY VARIABLES IN COMACDG, THESE ARE
!L  DOCUMENTED IN DOCACDG

!L  THIS SUBROUTINE CALLS NO OTHERS

!FPP$ NOCONCUR R

      IMPLICIT NONE
#include "acparm.h"
#include "comobs.h"
#include "comacdg.h"
#include "comacp.h"
#include "parvars.h"
!-----------------------------------------------------------------------
      INTEGER TNDV
      REAL OBS(TNDV)
      INTEGER J,JOB,JOBT,INOBS
      INTEGER IP_LAT,IP_LONG
      REAL LAT,LONG
!
!     3. SET UP FOR DIAGNOSTIC TYPE 1:- DETAILS ABOUT SELECTED OBS.
!     -------------------------------------------------------------
      NDACO = 0
      NDAC  = 0
      IF (LDIAGAC) THEN
        IF (LLDAC(1)) THEN

          IF (MODACO == 1) THEN

            IF (NDACPRT >  NDACP) NDACPRT=NDACP
#if defined(MPP)
          if(mype == 0)then
#endif
            PRINT *,                                                    &
     &      'MODACO=1 : Diagnostics on first ',NDACPRT,                 &
     &      ' Observations of each type.'
#if defined(MPP)
          endif
#endif

          ELSEIF (MODACO == 2) THEN

!           MDACO SET IN NAMELIST. COUNT THE NUMBER SET & CLEAR REST.
            DO J=1,NDACOP
              IF (MDACO(J) >  0) NDACO = NDACO+1
            ENDDO
            LLDAC(1) = NDACO >  0
            IF (NDACO <  NDACOP) THEN
              DO J=NDACO+1,NDACOP
                MDACO(J) = 0
              ENDDO
            ENDIF

          ELSEIF (MODACO == 3) THEN

!           MDACO SET UP FROM OBS IN LTD AREA DEFINED IN NAMELIST.
!           CHECK LATITUDE BAND (REMEMBER COLATITUDES ARE USED)
!           CHECK LONGITUDE BAND (ALLOW FOR AREAS CROSSING LONGITUDE 0)

            INOBS=0

            DO JOBT=1,NOBTYP

            IF (NOBS(JOBT) >  0) THEN

              IP_LAT  = MDISPOBT(JOBT)
              IP_LONG = MDISPOBT(JOBT) + NOBS(JOBT)

              IF (DLONGE >= DLONGW) THEN

                DO JOB=1,NOBS(JOBT)
                  LAT  = OBS(IP_LAT +JOB)
                  LONG = OBS(IP_LONG+JOB)
                  IF ( LAT  >= DLATN  .AND. LAT <= DLATS  .AND.         &
     &                 LONG >= DLONGW .AND. LONG <= DLONGE ) THEN
                    INOBS = INOBS+1
                  IF (INOBS <= NDACOP) MDACO(INOBS)=OBS_NO_ST(JOBT)+JOB
                  ENDIF
                ENDDO

              ELSE

!               AREA ASSUMED TO CROSS MERIDIAN (ELSE MISTAKE IN SET UP)
                DO JOB=1,NOBS(JOBT)
                  LAT  = OBS(IP_LAT +JOB)
                  LONG = OBS(IP_LONG+JOB)
                  IF ( LAT >= DLATN   .AND. LAT <= DLATS  .AND.         &
     &               (LONG >= DLONGW  .OR.  LONG <= DLONGE) ) THEN
                    INOBS = INOBS+1
                  IF (INOBS <= NDACOP) MDACO(INOBS)=OBS_NO_ST(JOBT)+JOB
                  ENDIF
                ENDDO

              ENDIF
            ENDIF
            ENDDO

            NDACO    = MIN(NDACOP,INOBS)
            LLDAC(1) = NDACO >  0
#if defined(MPP)
          if(mype == 0)then
#endif
            IF (NDACO <  INOBS) THEN
              PRINT *, 'MODACO=3 ',INOBS,                               &
     &        ' Observations found in limited area for diagnostics'
              PRINT *, '         ',NDACO,                               &
     &        ' of these listed and used.'
            ELSE
              PRINT *, 'MODACO=3 ',INOBS,                               &
     &        ' Observations found in limited area for diagnostics'
            ENDIF
#if defined(MPP)
          endif
#endif

          ENDIF
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE SETDAC
#endif
