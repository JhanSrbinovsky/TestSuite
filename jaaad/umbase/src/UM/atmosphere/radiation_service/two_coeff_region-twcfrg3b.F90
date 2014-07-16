#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate two-stream coefficients in the regions.
!
! Method:
!   The coefficients for each region are determined and
!       averaged.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.2             15-05-96                Original Code
!                                               (J. M. Edwards)
!       4.3             20-02-97                Calls to vector
!                                               searching routine
!                                               WHENFGT removed.
!                                               (J. M. Edwards)
!       4.5             04-02-99                optimisations added
!       4.5             11-06-98                Optimised version
!                                               (P. Burton)
!       5.3             04-10-01                Number of regions
!                                               passed explicitly.
!                                               (J. M. Edwards)
!  6.0  21/08/03  NEC SX-6 optimisation - NODEP directives.
!                 R Barnes & J-C Rioual.
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE TWO_COEFF_REGION(IERR                                  &
     &   , N_PROFILE, N_LAYER, N_CLOUD_TOP                              &
     &   , I_2STREAM, L_IR_SOURCE_QUAD, N_SOURCE_COEFF                  &
     &   , N_CLOUD_TYPE, FRAC_CLOUD                                     &
     &   , N_REGION, I_REGION_CLOUD, FRAC_REGION                        &
     &   , ASYMMETRY_FREE, OMEGA_FREE, TAU_FREE                         &
     &   , ASYMMETRY_CLOUD, OMEGA_CLOUD, TAU_CLOUD                      &
     &   , ISOLIR, SEC_0                                                &
     &   , TRANS, REFLECT, TRANS_0                                      &
     &   , SOURCE_COEFF                                                 &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     SIZES OF DUMMY ARRAYS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     INCLUDE COMDECKS.
#include "dimfix3a.h"
#include "spcrg3a.h"
#include "error3a.h"
#include "cldreg3a.h"
!
!
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             NUMBER OF LAYERS
     &   , N_CLOUD_TOP                                                  &
!             TOPMOST CLOUDY LAYER
     &   , ISOLIR                                                       &
!             SPECTRAL REGION
     &   , N_CLOUD_TYPE                                                 &
!             NUMBER OF TYPES OF CLOUDS
     &   , I_2STREAM                                                    &
!             TWO STREAM SCHEME
     &   , N_SOURCE_COEFF
!             NUMBER OF SOURCE COEFFICIENTS
!
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_REGION                                                     &
!             Number of cloudy regions
     &   , I_REGION_CLOUD(NPD_CLOUD_TYPE)
!             REGIONS IN WHICH TYPES OF CLOUDS FALL
!
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_IR_SOURCE_QUAD
!             USE A QUADRATIC SOURCE IN THE INFRA-RED
!
!     OPTICAL PROPERTIES OF LAYER:
      REAL                                                              &
                !, INTENT(IN)
     &     FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)           &
!             FRACTIONS OF DIFFERENT TYPES OF CLOUDS
     &   , FRAC_REGION(NPD_PROFILE, NPD_LAYER, NPD_REGION)              &
!             FRACTIONS OF TOTAL CLOUD OCCUPIED BY EACH REGION
     &   , ASYMMETRY_FREE(NPD_PROFILE, NPD_LAYER)                       &
!             CLEAR-SKY ASYMMETRY FACTOR
     &   , OMEGA_FREE(NPD_PROFILE, NPD_LAYER)                           &
!             CLEAR-SKY ALBEDO OF SINGLE SCATTERING
     &   , TAU_FREE(NPD_PROFILE, NPD_LAYER)                             &
!             CLEAR-SKY OPTICAL DEPTH
     &   , ASYMMETRY_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)      &
!             ASYMMETRY FACTOR
     &   , OMEGA_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)          &
!             ALBEDO OF SINGLE SCATTERING
     &   , TAU_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             OPTICAL DEPTH
!
!     SOLAR BEAM
      REAL                                                              &
                !, INTENT(IN)
     &     SEC_0(NPD_PROFILE)
!             SECANT OF ZENITH ANGLE
!
!
!     COEFFICIENTS IN THE TWO-STREAM EQUATIONS:
      REAL                                                              &
                !, INTENT(OUT)
     &     TRANS(NPD_PROFILE, NPD_LAYER, NPD_REGION)                    &
!             DIFFUSE TRANSMISSION COEFFICIENT
     &   , REFLECT(NPD_PROFILE, NPD_LAYER, NPD_REGION)                  &
!             DIFFUSE REFLECTION COEFFICIENT
     &   , TRANS_0(NPD_PROFILE, NPD_LAYER, NPD_REGION)                  &
!             DIRECT TRANSMISSION COEFFICIENT
     &   , SOURCE_COEFF(NPD_PROFILE, NPD_LAYER                          &
     &      , NPD_SOURCE_COEFF, NPD_REGION)
!             SOURCE COEFFICIENTS IN TWO-STREAM EQUATIONS
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , J                                                            &
!             LOOP VARIABLE
     &   , K                                                            &
!             LOOP VARIABLE
     &   , L                                                            &
!             LOOP VARIABLE
     &   , I_REGION
!             LOOP VARIABLE OVER REGIONS
!
!     COEFFICIENTS IN THE TWO-STREAM EQUATIONS:
      REAL                                                              &
                !, INTENT(OUT)
     &     TRANS_TEMP(NPD_PROFILE, NPD_LAYER)                           &
!             TEMPORARY DIFFUSE TRANSMISSION COEFFICIENT
     &   , REFLECT_TEMP(NPD_PROFILE, NPD_LAYER)                         &
!             TEMPORARY DIFFUSE REFLECTION COEFFICIENT
     &   , TRANS_0_TEMP(NPD_PROFILE, NPD_LAYER)                         &
!             TEMPORARY DIRECT TRANSMISSION COEFFICIENT
     &   , SOURCE_COEFF_TEMP(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF)
!             TEMPORARY SOURCE COEFFICIENTS IN TWO-STREAM EQUATIONS
!
!     VARIABLES FOR GATHERING:
      INTEGER                                                           &
     &     N_LIST                                                       &
!             NUMBER OF POINTS IN LIST
     &   , L_LIST(NPD_PROFILE)                                          &
!             LIST OF COLLECTED POINTS
     &   , LL
      REAL                                                              &
     &     TAU_GATHERED(NPD_PROFILE, NPD_LAYER)                         &
!             GATHERED OPTICAL DEPTH
     &   , OMEGA_GATHERED(NPD_PROFILE, NPD_LAYER)                       &
!             GATHERED ALEBDO OF SINGLE SCATTERING
     &   , ASYMMETRY_GATHERED(NPD_PROFILE, NPD_LAYER)                   &
!             GATHERED ASYMMETRY
     &   , SEC_0_GATHERED(NPD_PROFILE)                                  &
!             GATHERED ASYMMETRY
     &   , TMP_INV(NPD_PROFILE)
!             Temporary work array

      REAL,PARAMETER:: ONE=1.0
      REAL,PARAMETER:: SMALLP=TINY(ONE)

!
!     SUBROUTINES CALLED:
      EXTERNAL                                                          &
     &     TWO_COEFF
!
!     CRAY DIRECTIVES FOR THE WHOLE ROUTINE:
!     POINTS ARE NOT REPEATED IN THE INDEXING ARRAY, SO IT IS SAFE
!     TO VECTORIZE OVER INDIRECTLY ADDRESSED ARRAYS.
!fpp$ NODEPCHK R
!
!
!
!     DETERMINE THE OPTICAL PROPERTIES OF THE CLEAR-SKY REGIONS OF
!     THE LAYERS.
!
! DEPENDS ON: two_coeff
      CALL TWO_COEFF(IERR                                               &
     &   , N_PROFILE, 1, N_LAYER                                        &
     &   , I_2STREAM, L_IR_SOURCE_QUAD                                  &
     &   , ASYMMETRY_FREE, OMEGA_FREE, TAU_FREE                         &
     &   , ISOLIR, SEC_0                                                &
     &   , TRANS(1, 1, IP_REGION_CLEAR)                                 &
     &   , REFLECT(1, 1, IP_REGION_CLEAR)                               &
     &   , TRANS_0(1, 1, IP_REGION_CLEAR)                               &
     &   , SOURCE_COEFF(1, 1, 1, IP_REGION_CLEAR)                       &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
      IF (IERR /= I_NORMAL) RETURN
!
!
!     NOW DEAL WITH CLOUDS.
!
!     INITIALIZE THE FULL ARRAYS FOR CLOUDY REGIONS.
!
      DO I_REGION=1, N_REGION
         IF (I_REGION /= IP_REGION_CLEAR) THEN
            IF(ISOLIR /= IP_SOLAR) THEN
               DO I=N_CLOUD_TOP, N_LAYER
                  DO L=1, N_PROFILE
                     TRANS(L, I, I_REGION)=0.0E+00
                     REFLECT(L, I, I_REGION)=0.0E+00
                  ENDDO
               ENDDO
            ELSE IF(ISOLIR == IP_SOLAR) THEN
               DO I=N_CLOUD_TOP, N_LAYER
                  DO L=1, N_PROFILE
                     TRANS(L, I, I_REGION)=0.0E+00
                     REFLECT(L, I, I_REGION)=0.0E+00
                     TRANS_0(L, I, I_REGION)=0.0E+00
                  ENDDO
               ENDDO
            ENDIF
            DO J=1, N_SOURCE_COEFF
               DO I=N_CLOUD_TOP, N_LAYER
                  DO L=1, N_PROFILE
                     SOURCE_COEFF(L, I, J, I_REGION)=0.0E+00
                  ENDDO
               ENDDO
            ENDDO
!
!
         ENDIF
!
      ENDDO
!
!
!
!     CONSIDER EACH TYPE OF CLOUD IN TURN, CHECKING WHICH REGION IT
!     CONTRUBUTES TO AND FORM WEIGHTED SUMS OF CLOUD PROPERTIES.
!
      DO K=1, N_CLOUD_TYPE
!
!
!        SET THE REGION IN WHICH CLOUDS OF THIS TYPE ARE INCLUDED.
         I_REGION=I_REGION_CLOUD(K)
!
         DO I=N_CLOUD_TOP, N_LAYER
!
!           FORM A LIST OF POINTS WHERE CLOUD OF THIS TYPE EXISTS
!           ON THIS ROW FOR GATHERING.
            N_LIST=0
            DO L=1, N_PROFILE
               IF (FRAC_CLOUD(L, I, K).GT.0.0E+00) THEN     
                  N_LIST=N_LIST+1
                  L_LIST(N_LIST)=L
               ENDIF
            ENDDO
!
!
            IF (N_LIST >  0) THEN
!
!              GATHER THE OPTICAL PROPERTIES. THOUGH WE CONSIDER ONLY
!              ONE LAYER AT A TIME THE LOWER ROUTINES WILL OPERATE ON
!              ARRAYS WITH VERTICAL STRUCTURE, SO THE GATHERED ARRAYS
!              ARE TWO-DIMENSIONAL.
!
               DO L=1, N_LIST
                  TAU_GATHERED(L, I)                                    &
     &              =TAU_CLOUD(L_LIST(L), I, K)
                  OMEGA_GATHERED(L, I)                                  &
     &              =OMEGA_CLOUD(L_LIST(L), I, K)
                  ASYMMETRY_GATHERED(L, I)                              &
     &              =ASYMMETRY_CLOUD(L_LIST(L), I, K)
               ENDDO
               IF (ISOLIR == IP_SOLAR) THEN
                  DO L=1, N_LIST
                     SEC_0_GATHERED(L)=SEC_0(L_LIST(L))
                  ENDDO
               ENDIF
!
!
! DEPENDS ON: two_coeff
               CALL TWO_COEFF(IERR                                      &
     &            , N_LIST, I, I                                        &
     &            , I_2STREAM, L_IR_SOURCE_QUAD                         &
     &            , ASYMMETRY_GATHERED, OMEGA_GATHERED                  &
     &            , TAU_GATHERED                                        &
     &            , ISOLIR, SEC_0_GATHERED                              &
     &            , TRANS_TEMP, REFLECT_TEMP, TRANS_0_TEMP              &
     &            , SOURCE_COEFF_TEMP                                   &
     &            , NPD_PROFILE, NPD_LAYER                              &
     &            )
               IF (IERR /= I_NORMAL) RETURN
!
!
!CDIR NODEP
               DO L=1, N_LIST
                  LL=L_LIST(L)
                  TRANS(LL, I, I_REGION)=TRANS(LL, I, I_REGION)         &
     &               +FRAC_CLOUD(LL, I, K)*TRANS_TEMP(L, I)
                  REFLECT(LL, I, I_REGION)=REFLECT(LL, I, I_REGION)     &
     &               +FRAC_CLOUD(LL, I, K)*REFLECT_TEMP(L, I)
               ENDDO
               DO J=1, N_SOURCE_COEFF
!CDIR NODEP
                  DO L=1, N_LIST
                     LL=L_LIST(L)
                     SOURCE_COEFF(LL, I, J, I_REGION)                   &
     &                  =SOURCE_COEFF(LL, I, J, I_REGION)               &
     &                  +FRAC_CLOUD(LL, I, K)                           &
     &                  *SOURCE_COEFF_TEMP(L, I, J)
                  ENDDO
               ENDDO
               IF (ISOLIR == IP_SOLAR) THEN
!CDIR NODEP
                  DO L=1, N_LIST
                     LL=L_LIST(L)
                     TRANS_0(LL, I, I_REGION)=TRANS_0(LL, I, I_REGION)  &
     &                  +FRAC_CLOUD(LL, I, K)*TRANS_0_TEMP(L, I)
                  ENDDO
               ENDIF
!
            ENDIF
!
         ENDDO
      ENDDO
!
!
!     FINALLY, SCALE THE WEIGHTED SUMS BY THE CLOUD FRACTIONS.
      DO I_REGION=1, N_REGION
         IF (I_REGION /= IP_REGION_CLEAR) THEN
            DO I=N_CLOUD_TOP, N_LAYER
!
!              GATHER POINTS WITHIN THIS REGION.
               N_LIST=0
               DO L=1,N_PROFILE
                  IF (ABS(FRAC_REGION(L, I, I_REGION))                  &
     &                                     .GE.SMALLP) THEN
                     N_LIST=N_LIST+1
                     L_LIST(N_LIST)=L
                  ENDIF
               ENDDO
               IF(ISOLIR /= IP_SOLAR) THEN
!CDIR NODEP
                  DO L=1, N_LIST
                     LL=L_LIST(L)
                     TMP_INV(L)=1.0/FRAC_REGION(LL,I,I_REGION)
                     TRANS(LL, I, I_REGION)=TRANS(LL, I, I_REGION)      &
     &                    *TMP_INV(L)
                     REFLECT(LL, I, I_REGION)=REFLECT(LL, I, I_REGION)  &
     &                    *TMP_INV(L)
                  ENDDO
               ELSE IF(ISOLIR == IP_SOLAR) THEN
!CDIR NODEP
                  DO L=1, N_LIST
                     LL=L_LIST(L)
                     TMP_INV(L)=1.0/FRAC_REGION(LL,I,I_REGION)
                     TRANS(LL, I, I_REGION)=TRANS(LL, I, I_REGION)      &
     &                    *TMP_INV(L)
                     REFLECT(LL, I, I_REGION)=REFLECT(LL, I, I_REGION)  &
     &                    *TMP_INV(L)
                     TRANS_0(LL, I, I_REGION)=TRANS_0(LL, I, I_REGION)  &
     &                    *TMP_INV(L)
                  ENDDO
               ENDIF
               DO J=1, N_SOURCE_COEFF
!CDIR NODEP
                  DO L=1, N_LIST
                     LL=L_LIST(L)
                     SOURCE_COEFF(LL, I, J, I_REGION)                   &
     &                  =SOURCE_COEFF(LL, I, J, I_REGION)               &
     &                  *TMP_INV(L)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE TWO_COEFF_REGION
#endif
#endif
