#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the basic coefficients for the solar beam.
!
! Method:
!       Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.2             Nov. 96   T3E migration: CALL WHENFLT replaced
!                                  by portable fortran code.
!                                                S.J.Swarbrick
!       4.2             08-08-96                Extra two-stream option
!                                               added.
!                                               (J. M. Edwards)
!LL  4.5  27/04/98  Add Fujitsu vectorization directive.
!LL                                           RBarnes@ecmwf.int
!       5.1             04-04-00                Tolerances replaced
!                                               by F90 intrinsics.
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
! Fujitsu directive to encourage vectorization for whole routine
!OCL NOVREC
      SUBROUTINE SOLAR_COEFFICIENT_BASIC(IERR                           &
     &   , N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST                       &
     &   , OMEGA, ASYMMETRY, SEC_0                                      &
     &   , I_2STREAM                                                    &
     &   , SUM, DIFF, LAMBDA                                            &
     &   , GAMMA_UP, GAMMA_DOWN                                         &
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
     &     NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     INCLUDE COMDECKS.
#include "twostr3a.h"
#include "stdio3a.h"
#include "error3a.h"
!
!
!     DUMMY VARIABLES.
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
!             FIRST LAYER TO CONSIDER
     &   , I_2STREAM
!             TWO-STREAM SCHEME
!
      REAL                                                              &
                !, INTENT(IN)
     &     OMEGA(NPD_PROFILE, NPD_LAYER)                                &
!             ALBEDO OF SINGLE SCATTERING
     &   , ASYMMETRY(NPD_PROFILE, NPD_LAYER)                            &
!             ASYMMETRY
     &   , SEC_0(NPD_PROFILE)                                           &
!             SECANT OF SOLAR ZENITH ANGLE
     &   , SUM(NPD_PROFILE, NPD_LAYER)                                  &
!             SUM OF TWO-STREAM COEFFICIENTS
     &   , DIFF(NPD_PROFILE, NPD_LAYER)                                 &
!             DIFFERENCE OF TWO-STREAM COEFFICIENTS
     &   , LAMBDA(NPD_PROFILE, NPD_LAYER)
!             LAMBDA
!
!     BASIC TWO-STREAM COEFFICIENTS:
      REAL                                                              &
                !, INTENT(OUT)
     &     GAMMA_UP(NPD_PROFILE, NPD_LAYER)                             &
!             COEFFICIENT FOR UPWARD RADIATION
     &   , GAMMA_DOWN(NPD_PROFILE, NPD_LAYER)
!             COEFFICIENT FOR DOWNWAD RADIATION
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
     &   , KK                                                           &
!             TEMPORARY VARIABLE
     &   , N_INDEX                                                      &
!             NUMBER OF INDICES SATISFYING TEST
     &   , INDEX(NPD_PROFILE)
!             INDICES OF TESTED POINTS
      REAL                                                              &
     &     KSI_0(NPD_PROFILE, NPD_LAYER)                                &
!             DIFFERENCE IN SOLAR SCATTERING FRACTIONS
     &   , TEST_ARRAY(NPD_PROFILE)                                      &
!             ARRAY TO TEST
     &   , FACTOR                                                       &
!             TEMPORARY VARIABLE
     &   , TOL_SEC0
!             TOLERANCE USED TO PERTURB THE DIFFUSIVITY FACTOR
!             AWAY FROM THE SECANT OF THE ZENITH ANGLE
      REAL                                                              &
     &     ROOT_3
!             SQUARE ROOT OF 3
      PARAMETER(                                                        &
     &     ROOT_3=1.7320508075688772E+00                                &
     &   )
!
!     SUBROUTINES CALLED:
!
!     CRAY DIRECTIVES FOR THE WHOLE ROUTINE:
!     POINTS ARE NOT REPEATED IN THE INDEXING ARRAY, SO IT IS SAFE
!     TO VECTORIZE OVER INDIRECTLY ADDRESSED ARRAYS.
!fpp$ NODEPCHK R
!
!
!
!     IF LAMBDA IS TOO CLOSE TO SEC_0, FROM THE POINT OF VIEW OF THE
!     EQUATIONS, THE DIRECT AND DIFFUSE BEAMS ARE NOT PROPERLY
!     DISTINGUISHED AND A SINGULARITY WILL BE ENCOUNTERED. THE
!     DIFFUSIVITY FACTOR IS PERTURBED TO ACCOUNT FOR THIS. A TOLERANCE
!     TO DETECT THIS CONDITION IS FIRST DEFINED.
      TOL_SEC0=64.0E+00*EPSILON(TOL_SEC0)
!
!
      DO I=I_LAYER_FIRST, I_LAYER_LAST
         DO L=1, N_PROFILE
            TEST_ARRAY(L)=ABS(LAMBDA(L, I)-SEC_0(L))
         ENDDO
!
         N_INDEX=0
         DO L=1,N_PROFILE
           IF (TEST_ARRAY(L) <  TOL_SEC0) THEN
             N_INDEX=N_INDEX+1
             INDEX(N_INDEX)=L
           END IF
         END DO
!
!DIR$ IVDEP
         DO K=1, N_INDEX
            KK=INDEX(K)
            SUM(KK, I)=(1.0E+00+TOL_SEC0)*SUM(KK, I)
            DIFF(KK, I)=(1.0E+00+TOL_SEC0)*DIFF(KK, I)
            LAMBDA(KK, I)=(1.0E+00+TOL_SEC0)*LAMBDA(KK, I)
         ENDDO
      ENDDO
!
      IF ( (I_2STREAM == IP_EDDINGTON).OR.                              &
     &     (I_2STREAM == IP_ELSASSER).OR.                               &
     &     (I_2STREAM == IP_PIFM85).OR.                                 &
     &     (I_2STREAM == IP_2S_TEST).OR.                                &
     &     (I_2STREAM == IP_HEMI_MEAN).OR.                              &
     &     (I_2STREAM == IP_PIFM80) ) THEN
!
          DO I=I_LAYER_FIRST, I_LAYER_LAST
             DO L=1, N_PROFILE
                KSI_0(L, I)=1.5E+00*ASYMMETRY(L, I)/SEC_0(L)
             ENDDO
          ENDDO
!
       ELSE IF (I_2STREAM == IP_DISCRETE_ORD) THEN
!
          DO I=I_LAYER_FIRST, I_LAYER_LAST
             DO L=1, N_PROFILE
                KSI_0(L, I)=ROOT_3*ASYMMETRY(L, I)/SEC_0(L)
             ENDDO
          ENDDO
!
       ELSE
!
         WRITE(IU_ERR, '(/A)')                                          &
     &      '*** ERROR: AN ILLEGAL SOLAR TWO-STREAM SCHEME HAS '        &
     &      //'BEEN SELECTED.'
         IERR=I_ERR_FATAL
         RETURN
!
      ENDIF
!
!
!     DETERMINE THE BASIC SOLAR COEFFICIENTS FOR THE
!     TWO-STREAM EQUATIONS.
!
      DO I=I_LAYER_FIRST, I_LAYER_LAST
         DO L=1, N_PROFILE
            FACTOR=0.5E+00*OMEGA(L, I)*SEC_0(L)                         &
     &         /((LAMBDA(L, I)-SEC_0(L))*(LAMBDA(L, I)+SEC_0(L)))
            GAMMA_UP(L, I)=FACTOR*(SUM(L, I)-SEC_0(L)                   &
     &         -KSI_0(L, I)*(DIFF(L, I)-SEC_0(L)))
            GAMMA_DOWN(L, I)=FACTOR*(SUM(L, I)+SEC_0(L)                 &
     &         +KSI_0(L, I)*(DIFF(L, I)+SEC_0(L)))
         ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE SOLAR_COEFFICIENT_BASIC
#endif
#endif
