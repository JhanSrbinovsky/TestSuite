#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate coefficients in the two-stream equations.
!
! Method:
!       The basic two-stream coefficients in the differential equations
!       are calculated. These are then used to determine the
!       transmission and reflection coefficients. Coefficients for
!       determining the solar or infra-red source terms are calculated.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.2             Nov. 96   T3E migration: CALL WHENFGT replaced
!                                  by portable fortran code.
!                                                S.J.Swarbrick
!       4.4             04-02-99     optimisation added
!       4.5             11-06-98                Optimised Code
!                                               (P. Burton)
!       5.1             04-04-00                Tolerances replaced
!                                               by F90 intrinsics.
!                                               (J. M. Edwards)
!LL  5.3  24/09/01  Portability changes.    Z. Gardner
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE TWO_COEFF(IERR                                         &
     &   , N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST                       &
     &   , I_2STREAM, L_IR_SOURCE_QUAD                                  &
     &   , ASYMMETRY, OMEGA, TAU                                        &
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
     &   , I_LAYER_FIRST                                                &
!             FIRST LAYER TO CONSIDER
     &   , I_LAYER_LAST                                                 &
!             LAST LAYER TO CONSIDER
     &   , ISOLIR                                                       &
!             SPECTRAL REGION
     &   , I_2STREAM
!             TWO STREAM SCHEME
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_IR_SOURCE_QUAD
!               USE A QUADRATIC SOURCE FUNCTION
!
!     OPTICAL PROPERTIES OF LAYER:
      REAL                                                              &
                !, INTENT(IN)
     &     ASYMMETRY(NPD_PROFILE, NPD_LAYER)                            &
!             ASYMMETRY FACTOR
     &   , OMEGA(NPD_PROFILE, NPD_LAYER)                                &
!             ALBEDO OF SINGLE SCATTERING
     &   , TAU(NPD_PROFILE, NPD_LAYER)
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
     &     TRANS(NPD_PROFILE, NPD_LAYER)                                &
!             DIFFUSE TRANSMISSION COEFFICIENT
     &   , REFLECT(NPD_PROFILE, NPD_LAYER)                              &
!             DIFFUSE REFLECTION COEFFICIENT
     &   , TRANS_0(NPD_PROFILE, NPD_LAYER)                              &
!             DIRECT TRANSMISSION COEFFICIENT
     &   , SOURCE_COEFF(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF)
!             SOURCE COEFFICIENTS IN TWO-STREAM EQUATIONS
!
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , L                                                            &
!             LOOP VARIABLE
     &   , K                                                            &
!             LOOP VARIABLE
     &   , N_INDEX                                                      &
!             NUMBER OF INDICES SATISFYING TEST
     &   , INDEX(NPD_PROFILE)
!             INDICES OF TESTED POINTS
!
!     COEFFICIENTS IN THE TWO-STREAM EQUATIONS:
      REAL                                                              &
     &     LAMBDA(NPD_PROFILE, NPD_LAYER)                               &
!             COEFFICIENTS IN TWO-STREAM EQUATIONS
     &   , SUM(NPD_PROFILE, NPD_LAYER)                                  &
!             SUM OF ALPHA_1 AND ALPHA_2
     &   , DIFF(NPD_PROFILE, NPD_LAYER)                                 &
!             DIFFERENCE OF ALPHA_1 AND ALPHA_2
     &   , GAMMA_UP(NPD_PROFILE, NPD_LAYER)                             &
!             BASIC SOLAR COEFFICIENT FOR UPWARD RADIATION
     &   , GAMMA_DOWN(NPD_PROFILE, NPD_LAYER)
!             BASIC SOLAR COEFFICIENT FOR DOWNWARD RADIATION
!
      REAL                                                              &
     &     TARGET
!             TARGET TO SEARCH FOR

      REAL                                                              &
     &     TEMP(NPD_PROFILE)
!
!
!     SUBROUTINES CALLED:
      EXTERNAL                                                          &
     &     TWO_COEFF_BASIC, SOLAR_COEFFICIENT_BASIC                     &
     &   , TRANS_SOURCE_COEFF
!
!     CRAY DIRECTIVES FOR THE WHOLE ROUTINE:
!     POINTS ARE NOT REPEATED IN THE INDEXING ARRAY, SO IT IS SAFE
!     TO VECTORIZE OVER INDIRECTLY ADDRESSED ARRAYS.
!fpp$ NODEPCHK R
!
!
!
!     PERTURB THE SINGLE SCATTERING ALBEDO AWAY FROM 1 TO AVOID
!     LATER DIVISION BY 0.
      TARGET=1.0E+00-16.0E+00*EPSILON(TARGET)
      DO I=I_LAYER_FIRST, I_LAYER_LAST
!
         N_INDEX=0
         DO L   =1,N_PROFILE
           IF (OMEGA(L,I) >  TARGET) THEN
             N_INDEX =N_INDEX+1
             INDEX(N_INDEX)=L
           END IF
         END DO
!
         DO K=1, N_INDEX
            OMEGA(INDEX(K), I)=TARGET
         ENDDO
      ENDDO
!
!     CALCULATE THE BASIC TWO-STREAM COEFFICIENTS.
! DEPENDS ON: two_coeff_basic
      CALL TWO_COEFF_BASIC(IERR                                         &
     &   , N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST                       &
     &   , I_2STREAM                                                    &
     &   , ASYMMETRY, OMEGA                                             &
     &   , SUM, DIFF                                                    &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
      IF (IERR /= I_NORMAL) THEN
         RETURN
      ENDIF
!
!     LAMBDA IS NOW CALCULATED.
#if defined(VECTLIB)
      DO I=I_LAYER_FIRST, I_LAYER_LAST
         DO L=1, N_PROFILE
            TEMP(L)=SUM(L,I)*DIFF(L,I)
         ENDDO
         CALL SQRT_V(N_PROFILE,TEMP,LAMBDA(1,I))
      ENDDO
#else
      DO I=I_LAYER_FIRST, I_LAYER_LAST
         DO L=1, N_PROFILE
            LAMBDA(L, I)=SQRT(SUM(L, I)*DIFF(L, I))
         ENDDO
      ENDDO
#endif
!
!
!     CALCULATE THE BASIC COEFFICIENTS FOR THE SOLAR SOURCE TERMS.
      IF (ISOLIR == IP_SOLAR) THEN
!        LAMBDA MAY BE PERTURBED BY THIS ROUTINE TO AVOID
!        ILL-CONDITIONING FOR THE SINGULAR ZENITH ANGLE.
! DEPENDS ON: solar_coefficient_basic
         CALL SOLAR_COEFFICIENT_BASIC(IERR                              &
     &      , N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST                    &
     &      , OMEGA, ASYMMETRY, SEC_0                                   &
     &      , I_2STREAM                                                 &
     &      , SUM, DIFF, LAMBDA                                         &
     &      , GAMMA_UP, GAMMA_DOWN                                      &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )
         IF (IERR /= I_NORMAL) RETURN
      ENDIF
!
!
!     DETERMINE THE TRANSMISSION AND REFLECTION COEFFICIENTS.
! DEPENDS ON: trans_source_coeff
      CALL TRANS_SOURCE_COEFF(N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST    &
     &   , ISOLIR, L_IR_SOURCE_QUAD                                     &
     &   , TAU, SUM, DIFF, LAMBDA, SEC_0                                &
     &   , GAMMA_UP, GAMMA_DOWN                                         &
     &   , TRANS, REFLECT, TRANS_0, SOURCE_COEFF                        &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
!
!
      RETURN
      END SUBROUTINE TWO_COEFF
#endif
#endif
