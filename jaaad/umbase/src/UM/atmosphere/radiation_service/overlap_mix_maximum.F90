#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to find maximally overlapped energy transfer coefficients.
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
!       4.2             Nov. 96   T3E migration: CALL WHENFGT replaced
!                                  by portable fortran code.
!                                                S.J.Swarbrick
!       4.5             18-05-98                Reference to obsolete
!                                               solvers removed.
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
      SUBROUTINE OVERLAP_MIX_MAXIMUM(N_PROFILE, N_LAYER, N_CLOUD_TOP    &
     &   , ISOLIR, I_SOLVER                                             &
     &   , W_CLOUD, W_FREE                                              &
     &   , CLOUD_OVERLAP                                                &
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
#include "solver3a.h"
#include "clcfpt3a.h"
!
!     DUMMY ARGUMENTS.
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
     &   , I_SOLVER
!             SOLVER TO BE USED
      REAL                                                              &
                !, INTENT(IN)
     &     W_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUD AMOUNTS
!
      REAL                                                              &
                !, INTENT(OUT)
     &     W_FREE(NPD_PROFILE, NPD_LAYER)                               &
!             CLOUD-FREE AMOUNTS
     &   , CLOUD_OVERLAP(NPD_PROFILE, 0: NPD_LAYER, NPD_OVERLAP_COEFF)
!             COEFFICIENTS FOR TRANSFER OF ENERGY AT INTERFACE
!
!
!     LOCAL ARGUMENTS.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , L                                                            &
!             LOOP VARIABLE
     &   , K                                                            &
!             LOOP VARIABLE
     &   , KK                                                           &
!             LOOP VARIABLE
     &   , N_INDEX                                                      &
!             NUMBER OF INDICES IN WHEN TEST
     &   , INDEX(NPD_PROFILE)
!             INDICES OF WHEN TEST
      REAL                                                              &
     &     G_FF(NPD_PROFILE, 0: NPD_LAYER)                              &
!             BASIC TRANSFER COEFFICIENT
     &   , G_FC(NPD_PROFILE, 0: NPD_LAYER)                              &
!             BASIC TRANSFER COEFFICIENT
     &   , G_CF(NPD_PROFILE, 0: NPD_LAYER)                              &
!             BASIC TRANSFER COEFFICIENT
     &   , G_CC(NPD_PROFILE, 0: NPD_LAYER)                              &
!             BASIC TRANSFER COEFFICIENT
     &   , B_FF(NPD_PROFILE, 0: NPD_LAYER)                              &
!             BASIC TRANSFER COEFFICIENT
     &   , B_FC(NPD_PROFILE, 0: NPD_LAYER)                              &
!             BASIC TRANSFER COEFFICIENT
     &   , B_CF(NPD_PROFILE, 0: NPD_LAYER)                              &
!             BASIC TRANSFER COEFFICIENT
     &   , B_CC(NPD_PROFILE, 0: NPD_LAYER)
!             BASIC TRANSFER COEFFICIENT
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
!     SET THE FREE FRACTIONS IN EACH LAYER.
      DO I=1, N_CLOUD_TOP-1
         DO L=1, N_PROFILE
            W_FREE(L, I)=1.0E+00
         ENDDO
      ENDDO
      DO I=N_CLOUD_TOP, N_LAYER
         DO L=1, N_PROFILE
            W_FREE(L, I)=1.0E+00-W_CLOUD(L, I)
         ENDDO
      ENDDO
!
!     EVALUATE THE EXTENT OF OVERLAP BETWEEN LAYERS OF CLOUD
!     AT THE INTERFACE BETWEEN THE ITH AND (I+1)ST LAYER ON THE
!     ASSUMPTION OF MAXIMUM OVERLAP BETWEEN ADJACENT LAYERS.
!     THE TOP AND BOTTOM BOUNDARIES ARE EXCEPTIONAL.
!
!     HERE IS IS EASIEST TO FIND THE BASIC TRANSFER COEFFICIENTS
!     AND DEDUCE ALL OTHERS FROM THEM.
!
!     COEFFICIENTS ARE NOT REQUIRED IF THERE ARE NO CLOUDY LAYERS,
!     IN WHICH CASE N_CLOUD_TOP WILL BE SET TO N_LAYER+1.
      IF (N_CLOUD_TOP <= N_LAYER) THEN
         DO L=1, N_PROFILE
!
!           TOPMOST CLOUDY BOUNDARY:
            G_FF(L, N_CLOUD_TOP-1)=W_FREE(L, N_CLOUD_TOP)
            G_FC(L, N_CLOUD_TOP-1)=W_CLOUD(L, N_CLOUD_TOP)
            G_CF(L, N_CLOUD_TOP-1)=0.0E+00
            G_CC(L, N_CLOUD_TOP-1)=1.0E+00
            B_FF(L, N_CLOUD_TOP-1)=1.0E+00
            B_FC(L, N_CLOUD_TOP-1)=1.0E+00
            B_CF(L, N_CLOUD_TOP-1)=0.0E+00
            B_CC(L, N_CLOUD_TOP-1)=0.0E+00
!
!           BOTTOM BOUNDARY:
            G_FF(L, N_LAYER)=1.0E+00
            G_FC(L, N_LAYER)=0.0E+00
            G_CF(L, N_LAYER)=1.0E+00
            G_CC(L, N_LAYER)=0.0E+00
            B_FF(L, N_LAYER)=W_FREE(L, N_LAYER)
            B_FC(L, N_LAYER)=0.0E+00
            B_CF(L, N_LAYER)=W_CLOUD(L, N_LAYER)
            B_CC(L, N_LAYER)=1.0E+00
!
         ENDDO
      ENDIF
!
      DO I=N_CLOUD_TOP, N_LAYER-1
!
!        SET DEFAULT VALUES FOR COMPLETELY CLEAR OR CLOUDY LAYERS
!        AND OVERWRITE IN NORMAL CASES.
!
!        GAMMAS:
!
         DO L=1, N_PROFILE
            G_CC(L, I)=1.0E+00
            G_FF(L, I)=1.0E+00
         ENDDO
!
         N_INDEX=0
         DO L   =1,N_PROFILE
           IF (W_CLOUD(L,I) >  EPSILON(W_CLOUD)) THEN
             N_INDEX =N_INDEX+1
             INDEX(N_INDEX)=L
           END IF
         END DO
!
         DO K=1, N_INDEX
            KK=INDEX(K)
            G_CC(KK, I)=MIN(1.0E+00, W_CLOUD(KK, I+1)/W_CLOUD(KK, I))
         ENDDO
!
         N_INDEX=0
         DO L   =1,N_PROFILE
           IF (W_FREE(L,I) >  EPSILON(W_FREE)) THEN
             N_INDEX =N_INDEX+1
             INDEX(N_INDEX)=L
           END IF
         END DO
!
         DO K=1, N_INDEX
            KK=INDEX(K)
            G_FF(KK, I)=MIN(1.0E+00, W_FREE(KK, I+1)/W_FREE(KK, I))
         ENDDO
!        INFER THE OTHER VALUES FROM GENERIC RELATIONS.
         DO L=1, N_PROFILE
            G_FC(L, I)=1.0E+00-G_FF(L, I)
            G_CF(L, I)=1.0E+00-G_CC(L, I)
         ENDDO
!
!
!        BETAS:
!
         DO L=1, N_PROFILE
            B_CC(L, I)=1.0E+00
            B_FF(L, I)=1.0E+00
         ENDDO
!
         N_INDEX=0
         DO L   =1,N_PROFILE
           IF (W_CLOUD(L,I+1) >  EPSILON(W_CLOUD)) THEN
             N_INDEX =N_INDEX+1
             INDEX(N_INDEX)=L
           END IF
         END DO
!
         DO K=1, N_INDEX
            KK=INDEX(K)
            B_CC(KK, I)=MIN(1.0E+00, W_CLOUD(KK, I)/W_CLOUD(KK, I+1))
         ENDDO
!
         N_INDEX=0
         DO L   =1,N_PROFILE
           IF (W_FREE(L,I+1) >  EPSILON(W_FREE)) THEN
             N_INDEX =N_INDEX+1
             INDEX(N_INDEX)=L
           END IF
         END DO
!
         DO K=1, N_INDEX
            KK=INDEX(K)
            B_FF(KK, I)=MIN(1.0E+00, W_FREE(KK, I)/W_FREE(KK, I+1))
         ENDDO
!        INFER THE OTHER VALUSE FROM GENERIC RELATIONS.
         DO L=1, N_PROFILE
            B_FC(L, I)=1.0E+00-B_CC(L, I)
            B_CF(L, I)=1.0E+00-B_FF(L, I)
         ENDDO
!
      ENDDO
!
!     NOW CALCULATE THE OVERLAPPED ARRAYS USING THE BASIC COEFFICIENTS.
!
!     IN THE SOLAR REGION COEFFICIENTS FOR DOWNWARD COUPLING OF THE
!     FLUXES ARE REQUIRED. THESE COEFFICIENTS ARE ALSO NEEDED FOR
!     INFRA-RED CALCULATIONS WITH APPROXIMATE SCATTERING, OR
!     WHENEVER THE DIRECT SOLVER IS USED.
!
      IF ( (I_SOLVER == IP_SOLVER_MIX_DIRECT).OR.                       &
     &     (I_SOLVER == IP_SOLVER_MIX_DIRECT_HOGAN).OR.                 &
     &     (ISOLIR == IP_SOLAR).OR.                                     &
     &     ( (ISOLIR == IP_INFRA_RED).AND.                              &
     &     (I_SOLVER == IP_SOLVER_MIX_APP_SCAT) ) ) THEN
!
         DO I=N_CLOUD_TOP-1, N_LAYER
            DO L=1, N_PROFILE
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GFF)=G_FF(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GFC)=G_FC(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GCF)=G_CF(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GCC)=G_CC(L, I)
            ENDDO
         ENDDO
!
      ENDIF
!
!     WITH APPROXIMATE SCATTERING IN THE LONGWAVE, OR WHENEVER
!     A DIRECT SOLVER IS USED, THE CORRESPONDING
!     UPWARD COEFFICIENTS ARE NEEDED.
!
      IF ( (I_SOLVER == IP_SOLVER_MIX_DIRECT).OR.                       &
     &     (I_SOLVER == IP_SOLVER_MIX_DIRECT_HOGAN).OR.                 &
     &     ( (ISOLIR == IP_INFRA_RED).AND.                              &
     &     (I_SOLVER == IP_SOLVER_MIX_APP_SCAT) ) ) THEN
!
         DO I=N_CLOUD_TOP-1, N_LAYER
            DO L=1, N_PROFILE
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BFF)=B_FF(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BFC)=B_FC(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BCF)=B_CF(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BCC)=B_CC(L, I)
            ENDDO
         ENDDO
      ENDIF
!
!
      IF (I_SOLVER == IP_SOLVER_MIX_11) THEN
!
         DO I=N_CLOUD_TOP-1, N_LAYER
            DO L=1, N_PROFILE
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GM)                        &
     &            =G_FF(L, I)+G_CC(L, I)-G_CF(L, I)-G_FC(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GP)                        &
     &            =G_FF(L, I)-G_CC(L, I)+G_CF(L, I)-G_FC(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BM)                        &
     &            =B_FF(L, I)+B_CC(L, I)-B_CF(L, I)-B_FC(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BP)                        &
     &            =B_FF(L, I)-B_CC(L, I)-B_CF(L, I)+B_FC(L, I)
            ENDDO
         ENDDO
!
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE OVERLAP_MIX_MAXIMUM
#endif
#endif
