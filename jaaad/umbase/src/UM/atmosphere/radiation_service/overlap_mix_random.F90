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
!       4.1             17-05-96                Add code for new
!                                               solvers.
!       4.5             18-05-98                Reference to obsolete
!                                               solvers removed.
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE OVERLAP_MIX_RANDOM(N_PROFILE, N_LAYER, N_CLOUD_TOP     &
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
     &   , L
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
!     ASSUMPTION OF RANDOM OVERLAP BETWEEN ADJACENT LAYERS.
!     THE TOP AND BOTTOM BOUNDARIES ARE EXCEPTIONAL.
!
!     IN THE SOLAR REGION COEFFICIENTS FOR DOWNWARD COUPLING OF THE
!     FLUXES ARE REQUIRED. THESE COEFFICIENTS ARE ALSO NEEDED FOR
!     INFRA-RED CALCULATIONS WITH APPROXIMATE SCATTERING.
!
      IF ( (I_SOLVER == IP_SOLVER_MIX_DIRECT).OR.                       &
     &     (I_SOLVER == IP_SOLVER_MIX_DIRECT_HOGAN).OR.                 &
     &     (ISOLIR == IP_SOLAR).OR.                                     &
     &     ( (ISOLIR == IP_INFRA_RED).AND.                              &
     &       (I_SOLVER == IP_SOLVER_MIX_APP_SCAT) ) ) THEN
!
         DO I=N_CLOUD_TOP-1, N_LAYER-1
            DO L=1, N_PROFILE
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GFF)=W_FREE(L, I+1)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GFC)=W_CLOUD(L, I+1)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GCF)=W_FREE(L, I+1)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GCC)=W_CLOUD(L, I+1)
            ENDDO
         ENDDO
!
         DO L=1, N_PROFILE
!           BOTTOM BOUNDARY:
            CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_GFF)=1.0E+00
            CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_GFC)=0.0E+00
            CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_GCF)=1.0E+00
            CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_GCC)=0.0E+00
         ENDDO
!
      ENDIF
!
!     WITH APPROXIMATE SCATTERING IN THE LONGWAVE THE CORRESPONDING
!     UPWARD COEFFICIENTS ARE NEEDED.
!
      IF ( (I_SOLVER == IP_SOLVER_MIX_DIRECT).OR.                       &
     &     (I_SOLVER == IP_SOLVER_MIX_DIRECT_HOGAN).OR.                 &
     &     ( (ISOLIR == IP_INFRA_RED).AND.                              &
     &       (I_SOLVER == IP_SOLVER_MIX_APP_SCAT) ) ) THEN
!
         DO L=1, N_PROFILE
!           TOP CLOUDY BOUNDARY:
            CLOUD_OVERLAP(L, N_CLOUD_TOP-1, IP_CLOVLP_BFF)=1.0E+00
            CLOUD_OVERLAP(L, N_CLOUD_TOP-1, IP_CLOVLP_BFC)=1.0E+00
            CLOUD_OVERLAP(L, N_CLOUD_TOP-1, IP_CLOVLP_BCF)=0.0E+00
            CLOUD_OVERLAP(L, N_CLOUD_TOP-1, IP_CLOVLP_BCC)=0.0E+00
         ENDDO
!
         DO I=N_CLOUD_TOP, N_LAYER
            DO L=1, N_PROFILE
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BFF)=W_FREE(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BFC)=W_FREE(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BCF)=W_CLOUD(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BCC)=W_CLOUD(L, I)
            ENDDO
         ENDDO
      ENDIF
!
!
      IF (I_SOLVER == IP_SOLVER_MIX_11) THEN
!
         DO L=1, N_PROFILE
!
!           TOP BOUNDARY:
            CLOUD_OVERLAP(L, N_CLOUD_TOP-1, IP_CLOVLP_GM)=0.0E+00
            CLOUD_OVERLAP(L, N_CLOUD_TOP-1, IP_CLOVLP_GP)               &
     &         =2.0E+00-4.0E+00*W_CLOUD(L, N_CLOUD_TOP)
            CLOUD_OVERLAP(L, N_CLOUD_TOP-1, IP_CLOVLP_BM)=0.0E+00
            CLOUD_OVERLAP(L, N_CLOUD_TOP-1, IP_CLOVLP_BP)=2.0E+00
!
!           BOTTOM BOUNDARY:
            CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_GM)=0.0E+00
            CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_GP)=2.0E+00
            CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_BM)=0.0E+00
            CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_BP)                     &
     &         =2.0E+00-4.0E+00*W_CLOUD(L, N_LAYER)
!
         ENDDO
!
         DO I=N_CLOUD_TOP, N_LAYER-1
            DO L=1, N_PROFILE
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GM)=0.0E+00
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GP)                        &
     &            =2.0E+00-4.0E+00*W_CLOUD(L, I+1)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BM)=0.0E+00
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BP)                        &
     &            =2.0E+00-4.0E+00*W_CLOUD(L, I)
            ENDDO
         ENDDO
!
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE OVERLAP_MIX_RANDOM
#endif
#endif
