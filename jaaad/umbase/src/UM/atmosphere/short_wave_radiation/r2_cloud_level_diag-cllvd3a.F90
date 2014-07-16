#if defined(A01_3A) || defined(A01_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Subroutine to calculate an observed effective radius.
!
! Purpose:
!   An effective radius as observed from above the cloud-top is
!   calculated.
!
! Method:
!   For each type of cloud containing water in any layer the effective
!   radius is weighted with the product of the area of the cloud and the
!   probability that light emitted from the cloud reaches the observing
!   instrument.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.4             24-01-97                Original Code
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Move initialization
!                                               of diagnostics to a
!                                               higher level.
!                                               (J. M. Edwards)
!       5.3             11-10-01                Convert output
!                                               diagnostics
!                                               to 2-D fields to match
!                                               structure.
!                                               (J. M. Edwards)
!       5.4             29-05-02                Allow for diagnosis
!                                               of effective radius
!                                               only from clouds which
!                                               are warmer than the
!                                               freezing point as done
!                                               in AVHRR retrievals.
!                                               (A. Jones)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_CLOUD_LEVEL_DIAG(IERR, N_PROFILE, N_LAYER, NCLDS    &
     &   , I_GATHER                                                     &
     &   , I_CLOUD, I_CLOUD_REPRESENTATION                              &
     &   , T, W_CLOUD, FRAC_CLOUD, L_ALL_TEMPS                          &
     &   , CONDENSED_MIX_RATIO, CONDENSED_RE                            &
     &   , L_OBSERVED_RE, WEIGHTED_RE, SUM_WEIGHT_RE                    &
     &   , col_list, row_list, row_length, rows                         &
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER                            &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED.
#include "error3a.h"
#include "stdio3a.h"
#include "dimfix3a.h"
#include "cldcmp3a.h"
#include "cldtyp3a.h"
#include "clrepp3a.h"
#include "clschm3a.h"
#include "c_0_dg_c.h"
!
!
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
!
!     DIMENSIONS OF ARRAYS:
      Integer, Intent(IN) :: row_length
!                              Number of grid-points in EW-direction
!                              in the local domain
      Integer, Intent(IN) :: rows
!                              Number of grid-points in NS-direction
!                              in the local domain
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_FIELD                                                    &
!             SIZE OF ARRAY OF ARRAYS PASSED FROM MAIN CODE
     &   , NPD_PROFILE                                                  &
!             SIZE OF ARRAY OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     ACTUAL SIZES USED:
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             Number of atmospheric layers in radiation
     &   , NCLDS                                                        &
!             NUMBER OF CLOUDY LEVELS
     &   , I_GATHER(NPD_FIELD)
!             LIST OF GATHERED POINTS
      Integer, Intent(IN) :: col_list(npd_field)
!                              EW indices of gathered points in the
!                              2-D domain
      Integer, Intent(IN) :: row_list(npd_field)
!                              NS indices of gathered points in the
!                              2-D domain
!
!     LOGICAL FLAGS FOR DIAGNOSTICS
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_OBSERVED_RE                                                &
!             FLAG TO ENABLE DIAGNOSIS OF EFFECTIVE RADIUS SEEN FROM
!             SPACE (N.B. THE ROUTINE IS AT PRESENT CALLED ONLY IF
!             THIS IS TRUE, BUT ITS PRESENCE HERE ALLOWS FOR POSSIBLE
!             FUTURE EXTENSION OF THE ROUTINE).
     &   , L_ALL_TEMPS
!             IF TRUE, THE ROUTINE HAS BEEN CALLED TO OBTAIN DIAGNOSTICS
!             FOR CLOUDS CONSISTING OF LIQUID WATER AT ANY TEMPERATURE
!             (AS DONE IN MODIS RETRIEVALS). IF FALSE, ONLY CLOUDS WITH
!             TEMPERATURES ABOVE FREEZING ARE TO BE DIAGNOSED (AS DONE
!             IN AVHRR RETRIEVALS).
!
!     REPRESENTATION OF CLOUDS
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_CLOUD_REPRESENTATION                                       &
!             REPRESENTATION OF CLOUDS
     &   , I_CLOUD
!             TREATMENT OF OVERLAPS
!
      REAL                                                              &
                !, INTENT(IN)
     &     W_CLOUD(NPD_PROFILE, NPD_LAYER)                              &
!             TOTAL AMOUNTS OF CLOUD
     &   , FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)           &
!             FRACTION OF TYPES OF CLOUD
     &   , CONDENSED_RE(NPD_PROFILE, 0: NPD_LAYER, NPD_CLOUD_COMPONENT) &
!             EFFECTIVE RADII OF CLOUDY COMPONENTS
     &   , CONDENSED_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER                &
     &      , NPD_CLOUD_COMPONENT)                                      &
!             MASS MIXING RATIOS OF CONDENSED COMPONENTS
     &   , T(NPD_PROFILE, 0: NPD_LAYER)
!             TEMPERATURE
!
      REAL                                                              &
                !, INTENT(OUT)
     &     WEIGHTED_RE(row_length, rows)                                &
!             WEIGHTED SUM OF EFFECTIVE RADIUS AND WEIGHTING FUNCTION
     &   , SUM_WEIGHT_RE(row_length, rows)
!             SUM OF WEIGHTS FOR EFFECTIVE RADIUS
!
!
!
!     LOCAL VARIABLES:
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , L                                                            &
!             LOOP VARIABLE
     &   , I_INV
!             INVERTED LOOP INDEX
      REAL                                                              &
     &     TRANS_OVERLYING_SPACE(NPD_PROFILE)                           &
!             PROBABILITY OF A PHOTON IN CLEAR AIR IN THE LEVEL ABOVE
!             THE CURRENT ONE REACHING SPACE
     &   , AREA_EXPOSED(NPD_PROFILE)                                    &
!             TOTAL AREA OF CLOUD IN THE CURRENT LAYER EXPOSED TO
!             CLEAR AIR IN THE LAYER ABOVE
     &   , AREA_EXPOSED_ST(NPD_PROFILE)                                 &
!             TOTAL AREA OF STRATIFORM CLOUD IN THE CURRENT LAYER
!             EXPOSED TO CLEAR AIR IN THE LAYER ABOVE
     &   , AREA_EXPOSED_CNV(NPD_PROFILE)                                &
!             TOTAL AREA OF CONVECTIVE CLOUD IN THE CURRENT LAYER
!             EXPOSED TO CLEAR AIR IN THE LAYER ABOVE
     &   , AREA_CLEAR_ABOVE(NPD_PROFILE)                                &
!             AREA OF THE CLEAR SKY REGION IN THE LAYER ABOVE
     &   , AREA_STRAT(NPD_PROFILE)                                      &
!             AREA OF STRATIFORM CLOUD IN THE CURRENT LAYER
     &   , AREA_STRAT_ABOVE(NPD_PROFILE)                                &
!             AREA OF STRATIFORM CLOUD IN THE LAYER ABOVE
     &   , AREA_CONV(NPD_PROFILE)                                       &
!             AREA OF CONVECTIVE CLOUD IN THE CURRENT LAYER
     &   , AREA_CONV_ABOVE(NPD_PROFILE)                                 &
!             AREA OF CONVECTIVE CLOUD IN THE LAYER ABOVE
     &   , AREA_CLEAR_CLEAR(NPD_PROFILE)                                &
!             AREA OF BOUNDARY WHERE CLEAR SKY OVERLIES CLEAR SKY
     &   , AREA_CLEAR(NPD_PROFILE)                                      &
!             AREA OF CLEAR SKY IN THE CURRENT LAYER
!             DOWN TO A LEVEL
     &   , AREA_UNCORRELATED(NPD_PROFILE)                               &
!             UNCORRELATED REGION ON THE INTERFACE
     &   , WEIGHTED_RE_G(NPD_PROFILE)                                   &
!             WEIGHTED SUM OF EFFECTIVE RADIUS AND WEIGHTING FUNCTION
     &   , SUM_WEIGHT_RE_G(NPD_PROFILE)
!             SUM OF WEIGHTS FOR EFFECTIVE RADIUS
!
!     VARIABLES FOR GATHERING
      INTEGER                                                           &
     &     N_LIST                                                       &
!             NUMBER OF POINTS IN LIST
     &   , L_LIST(NPD_PROFILE)
!             INDICES OF POINTS IN LIST
!
!     INDICATOR FUNCTION
      REAL                                                              &
     &     CHI_CNV(NPD_PROFILE)                                         &
!             CONVECTIVE INDICATOR FUNCTION
     &   , CHI_ST(NPD_PROFILE)
!             STRATIFORM INDICATOR FUNCTION
!
!
!
!
!     INITIALIZATION OF LOCAL FIELDS.
!
      IF (L_OBSERVED_RE) THEN
         DO L=1, N_PROFILE
            WEIGHTED_RE_G(L)=0.0E+00
            SUM_WEIGHT_RE_G(L)=0.0E+00
         ENDDO
      ENDIF
!
!     INITIALIZE THE TRANSMISION ABOVE CLOUDS.
      DO L=1, N_PROFILE
         TRANS_OVERLYING_SPACE(L)=1.0E+00
         AREA_CLEAR_ABOVE(L)=1.0E+00
      ENDDO
      IF (L_OBSERVED_RE.AND.(I_CLOUD == IP_CLOUD_TRIPLE)) THEN
         DO L=1, N_PROFILE
            AREA_STRAT_ABOVE(L)=0.0E+00
            AREA_CONV_ABOVE(L)=0.0E+00
         ENDDO
      ENDIF
!
!     STEP DOWN THROUGH THE ATMOSPHERE CALCULATING CONTRIBUTIONS TO
!     THE DIAGNOSTICS AND SUBSEQUENTLY ALLOWING FOR TRANSMISSION
!     THROUGH THE CURRENT LAYER.
!
      DO I=NCLDS, 1, -1
         I_INV=N_LAYER+1-I
!
         DO L=1, N_PROFILE
            AREA_CLEAR(L)=1.0E+00-W_CLOUD(L, I_INV)
         ENDDO
!
!        CALCULATE THE LOCAL AREA OF CLOUD RADIATING INTO CLEAR AIR.
         IF (I_CLOUD == IP_CLOUD_MIX_RANDOM) THEN
            DO L=1, N_PROFILE
               AREA_EXPOSED(L)=W_CLOUD(L, I_INV)                        &
     &            *AREA_CLEAR_ABOVE(L)
            ENDDO
         ELSE IF ( (I_CLOUD == IP_CLOUD_MIX_MAX).OR.                    &
     &             (I_CLOUD == IP_CLOUD_TRIPLE) ) THEN
            DO L=1, N_PROFILE
               AREA_EXPOSED(L)=MAX(0.0E+00, (W_CLOUD(L, I_INV)          &
     &            +AREA_CLEAR_ABOVE(L)-1.0E+00))
            ENDDO
         ENDIF
!
!
!
!
         IF (L_OBSERVED_RE) THEN
!
!
!
            IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CONV_STRAT) THEN
!
!
               IF ( (I_CLOUD == IP_CLOUD_MIX_MAX).OR.                   &
     &              (I_CLOUD == IP_CLOUD_MIX_RANDOM) ) THEN
!
!                 IF THE OVERLAP OF CONVECTIVE CLOUD IS NOT ASSUMED
!                 TO BE COHERENT THE OVERALL EXPOSED AREA MAY BE
!                 PARTITIONED ACCORDING TO THE FRACTIONAL
!                 CONTRIBUTIONS OF CLOUD IN THE CURRENT LAYER.
!
                  DO L=1, N_PROFILE
                     AREA_EXPOSED_ST(L)=AREA_EXPOSED(L)                 &
     &                  *FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_STRAT)
                     AREA_EXPOSED_CNV(L)=AREA_EXPOSED(L)                &
     &                  *FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_CONV)
                  ENDDO
!
               ELSE IF (I_CLOUD == IP_CLOUD_TRIPLE) THEN
!
!                 HERE, THE DIFFERENT TYPES OF CLOUDS OVERLAP
!                 COHERENTLY SO STRATIFORM CLOUD WILL BE EXPOSED
!                 ONLY IF THERE IS LESS STRATIFORM CLOUD IN THE
!                 LAYER ABOVE AND MORE CLEAR AIR IN THE LAYER ABOVE:
!                 UNDER THESE CONDITIONS THE NON-CORRELATED AREAS
!                 OVERLAP RANDOMLY.
!
                  DO L=1, N_PROFILE
                     AREA_STRAT(L)=W_CLOUD(L, I_INV)                    &
     &                  *FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_STRAT)
                     AREA_CONV(L)=W_CLOUD(L, I_INV)                     &
     &                 *FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_CONV)
                     AREA_UNCORRELATED(L)=1.0E+00                       &
     &                 -MIN(AREA_CLEAR(L), AREA_CLEAR_ABOVE(L))         &
     &                 -MIN(AREA_STRAT(L), AREA_STRAT_ABOVE(L))         &
     &                 -MIN(AREA_CONV(L), AREA_CONV_ABOVE(L))
!                    First find the area of uncorrelated
!                    stratiform cloud.
                     AREA_EXPOSED_ST(L)=MAX(0.0E+00                     &
     &                  , (AREA_STRAT(L)-AREA_STRAT_ABOVE(L)))
                     AREA_EXPOSED_ST(L)=MAX(0.0E+00, AREA_EXPOSED_ST(L) &
     &                  *(AREA_CLEAR_ABOVE(L)-AREA_CLEAR(L)))
!                    Now normalize within the uncorrelated region.
!                    If the uncorrelated area is 0 the exposed area
!                    must be 0, so no second branch of the IF-test
!                    is required.
                     IF (AREA_UNCORRELATED(L) >  0.0E+00)               &
     &                 AREA_EXPOSED_ST(L)                               &
     &                   =AREA_EXPOSED_ST(L)/AREA_UNCORRELATED(L)
                     AREA_EXPOSED_CNV(L)                                &
     &                  =AREA_EXPOSED(L)-AREA_EXPOSED_ST(L)
                  ENDDO
               ELSE
                  WRITE(IU_ERR, '(/A)')                                 &
     &               '*** ERROR: THE DIAGNOSTIC OF OBSERVED RE HAS NOT '&
     &               //'BEEN IMPLEMENTED WITH THIS OVERLAP OPTION.'
                  IERR=I_ERR_FATAL
                  RETURN
               ENDIF
!
!              THE INDICATOR FUNCTIONS FOR LIQUID WATER IN
!              CONVECTIVE OR STRAIFORM CLOUDS ARE SET TO 1
!              IF THERE IS ANY LIQUID WATER AND TO 0 OTHERWISE.
               DO L=1, N_PROFILE
                  IF (CONDENSED_MIX_RATIO(L, I_INV, IP_CLCMP_CNV_WATER) &
     &                >  0.0E+00) THEN
                     CHI_CNV(L)=1.0E+00
                  ELSE
                     CHI_CNV(L)=0.0E+00
                  ENDIF
                  IF (CONDENSED_MIX_RATIO(L, I_INV, IP_CLCMP_ST_WATER)  &
     &                >  0.0E+00) THEN
                     CHI_ST(L)=1.0E+00
                  ELSE
                     CHI_ST(L)=0.0E+00
                  ENDIF
               ENDDO
!
!              INCLUDE CONTRIBUTIONS FROM CONVECTIVE AND STRATIFORM
!              WATER CLOUDS.
               DO L=1, N_PROFILE
                  WEIGHTED_RE_G(L)=WEIGHTED_RE_G(L)                     &
     &               +TRANS_OVERLYING_SPACE(L)                          &
     &               *(AREA_EXPOSED_CNV(L)*CHI_CNV(L)                   &
     &               *CONDENSED_RE(L, I_INV, IP_CLCMP_CNV_WATER)        &
     &               +AREA_EXPOSED_ST(L)*CHI_ST(L)                      &
     &               *CONDENSED_RE(L, I_INV, IP_CLCMP_ST_WATER))
                  SUM_WEIGHT_RE_G(L)=SUM_WEIGHT_RE_G(L)                 &
     &               +TRANS_OVERLYING_SPACE(L)                          &
     &               *(AREA_EXPOSED_CNV(L)*CHI_CNV(L)                   &
     &               +AREA_EXPOSED_ST(L)*CHI_ST(L))
               ENDDO
!
            ELSE IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CSIW) THEN
!
               IF ( (I_CLOUD == IP_CLOUD_MIX_MAX).OR.                   &
     &              (I_CLOUD == IP_CLOUD_MIX_RANDOM) ) THEN
!
!                 IF THE OVERLAP OF CONVECTIVE CLOUD IS NOT ASSUMED
!                 TO BE COHERENT THE OVERALL EXPOSED AREA MAY BE
!                 PARTITIONED ACCORDING TO THE FRACTIONAL
!                 CONTRIBUTIONS OF CLOUD IN THE CURRENT LAYER.
!                 THE EXPOSED AREAS INCLUDE ONLY THE PARTS OF THE
!                 CLOUDS CONTAINING WATER DROPLETS.
!
                  DO L=1, N_PROFILE
                     AREA_EXPOSED_ST(L)=AREA_EXPOSED(L)                 &
     &                  *FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_SW)
                     AREA_EXPOSED_CNV(L)=AREA_EXPOSED(L)                &
     &                  *FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_CW)
                  ENDDO
!
               ELSE IF (I_CLOUD == IP_CLOUD_TRIPLE) THEN
!
!                 HERE, THE DIFFERENT TYPES OF CLOUDS OVERLAP
!                 COHERENTLY SO STRATIFORM CLOUD WILL BE EXPOSED
!                 ONLY IF THERE IS LESS STRATIFORM CLOUD IN THE
!                 LAYER ABOVE AND MORE CLEAR AIR IN THE LAYER ABOVE:
!                 UNDER THESE CONDITIONS THE NON-CORRELATED AREAS
!                 OVERLAP RANDOMLY.
!                 THE ACTUAL EXPOSED AREAS OF CONVECTIVE OR
!                 STRATIFORM CLOUD MUST THEN BE WEIGHTED BY FACTORS
!                 REPRESENTING THE LIQUID PORTION OF EACH CLOUD, SINCE
!                 NOTHING IS RETRIEVED OVER ICE. (THE HORIZONTAL
!                 ARRANGEMENT OF ICE AND WATER WITHIN EITHER TYPE OF
!                 CLOUD IS RANDOM).
!
                  DO L=1, N_PROFILE
!
                     AREA_STRAT(L)=W_CLOUD(L, I_INV)                    &
     &                  *(FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_SW)        &
     &                  +FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_SI))
                     AREA_CONV(L)=W_CLOUD(L, I_INV)                     &
     &                  *(FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_CW)        &
     &                  +FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_CI))
                     AREA_UNCORRELATED(L)=1.0E+00                       &
     &                  -MIN(AREA_CLEAR(L), AREA_CLEAR_ABOVE(L))        &
     &                  -MIN(AREA_STRAT(L), AREA_STRAT_ABOVE(L))        &
     &                  -MIN(AREA_CONV(L), AREA_CONV_ABOVE(L))
                     AREA_EXPOSED_ST(L)=MAX(0.0E+00                     &
     &                  , (AREA_STRAT(L)-AREA_STRAT_ABOVE(L)))
                     IF (AREA_UNCORRELATED(L) >  0.0E+00) THEN
                        AREA_EXPOSED_ST(L)                              &
     &                     =MAX(0.0E+00, AREA_EXPOSED_ST(L)             &
     &                     *(AREA_CLEAR_ABOVE(L)-AREA_CLEAR(L)))        &
     &                     /AREA_UNCORRELATED(L)
                     ELSE
                        AREA_EXPOSED_ST(L)=0.0E+00
                     ENDIF
                     AREA_EXPOSED_CNV(L)                                &
     &                  =AREA_EXPOSED(L)-AREA_EXPOSED_ST(L)
!
                     IF (FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_CW)         &
     &                   >  0.0E+00) THEN
                        AREA_EXPOSED_CNV(L)=AREA_EXPOSED_CNV(L)         &
     &                     /(1.0E+00                                    &
     &                     +FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_CI)      &
     &                     /FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_CW))
                     ELSE
                        AREA_EXPOSED_CNV(L)=0.0E+00
                     ENDIF
!
                     IF (FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_SW)         &
     &                   >  0.0E+00) THEN
                        AREA_EXPOSED_ST(L)=AREA_EXPOSED_ST(L)           &
     &                     /(1.0E+00                                    &
     &                     +FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_SI)      &
     &                     /FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_SW))
                     ELSE
                        AREA_EXPOSED_ST(L)=0.0E+00
                     ENDIF
!
                  ENDDO
               ELSE
                  WRITE(IU_ERR, '(/A)')                                 &
     &               '*** ERROR: THE DIAGNOSTIC OF OBSERVED RE HAS NOT '&
     &               //'BEEN IMPLEMENTED WITH THIS OVERLAP OPTION.'
                  IERR=I_ERR_FATAL
                  RETURN
               ENDIF
!
!
               DO L=1, N_PROFILE
!
                IF ((T(L, I_INV)  >   TM) .OR. L_ALL_TEMPS) THEN
                  WEIGHTED_RE_G(L)=WEIGHTED_RE_G(L)                     &
     &               +TRANS_OVERLYING_SPACE(L)                          &
     &               *(AREA_EXPOSED_CNV(L)                              &
     &               *CONDENSED_RE(L, I_INV, IP_CLCMP_CNV_WATER)        &
     &               +AREA_EXPOSED_ST(L)                                &
     &               *CONDENSED_RE(L, I_INV, IP_CLCMP_ST_WATER))
                  SUM_WEIGHT_RE_G(L)=SUM_WEIGHT_RE_G(L)                 &
     &               +TRANS_OVERLYING_SPACE(L)                          &
     &               *(AREA_EXPOSED_CNV(L)+AREA_EXPOSED_ST(L))
                ENDIF
               ENDDO
!
            ENDIF
!
!
         ENDIF
!
!
!
!        ADVANCE THE STORED QUANTITIES REFFERRING TO OVERLYING LAYERS.
!
!
!        THE TRANSMISSION TO SPACE CURRENTLY HOLDS THE PROBABILITY THAT
!        A PHOTON TRAVELLING UPWARDS IN THE CLEAR AIR IN THE LAYER ABOVE
!        WILL ESCAPE TO SPACE WITHOUT ENCOUNTERING A CLOUD. TO ADVANCE
!        THIS TO THE CURRENT LAYER IT MUST BE MULTIPLIED BY A FACTOR
!        REPRESENTING THE OVERLAP ASSUMPTION AT THE TOP OF THE PRESENT
!        LAYER.
!
         IF (I_CLOUD == IP_CLOUD_MIX_RANDOM) THEN
!
            DO L=1, N_PROFILE
               TRANS_OVERLYING_SPACE(L)=TRANS_OVERLYING_SPACE(L)        &
     &            *AREA_CLEAR_ABOVE(L)
            ENDDO
!
         ELSE IF ( (I_CLOUD == IP_CLOUD_MIX_MAX).OR.                    &
     &             (I_CLOUD == IP_CLOUD_TRIPLE) ) THEN
!
            DO L=1, N_PROFILE
               AREA_CLEAR_CLEAR(L)=MIN(AREA_CLEAR(L)                    &
     &            , AREA_CLEAR_ABOVE(L))
               IF (AREA_CLEAR(L) >  0.0E+00) THEN
                  TRANS_OVERLYING_SPACE(L)=TRANS_OVERLYING_SPACE(L)     &
     &               *AREA_CLEAR_CLEAR(L)/AREA_CLEAR(L)
               ELSE
                  TRANS_OVERLYING_SPACE(L)=0.0E+00
               ENDIF
            ENDDO
!
         ENDIF
!
!        ADVANCE THE AREAS OF CLOUD.
         DO L=1, N_PROFILE
            AREA_CLEAR_ABOVE(L)=AREA_CLEAR(L)
         ENDDO
         IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CONV_STRAT) THEN
            DO L=1, N_PROFILE
               AREA_STRAT_ABOVE(L)=W_CLOUD(L, I_INV)                    &
     &            *FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_STRAT)
            ENDDO
         ELSE IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CSIW) THEN
            DO L=1, N_PROFILE
               AREA_STRAT_ABOVE(L)=AREA_STRAT(L)
               AREA_CONV_ABOVE(L)=AREA_CONV(L)
            ENDDO
         ENDIF
!
      ENDDO
!
!
!
      IF (L_OBSERVED_RE) THEN
!        SCATTER THE DIAGNOSTICS BACK TO THE OUTPUT ARRAYS AND CONVERT
!        TO MICRONS (TO AVOID FIELDS BEING CORRUPTED BY PACKING).
         DO L=1, N_PROFILE
            WEIGHTED_RE(col_list(l), row_list(l))                       &
     &        =1.0E+06*WEIGHTED_RE_G(L)
            SUM_WEIGHT_RE(col_list(l), row_list(l))                     &
     &        =SUM_WEIGHT_RE_G(L)
         ENDDO
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE R2_CLOUD_LEVEL_DIAG
#endif
