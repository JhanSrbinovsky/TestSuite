#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calcaulate IR source function for differential flux.
!
! Method:
!       The linear contribution to the source function is proportional
!       to the absorption divided by the optical depth. A tolerance is
!       added to the optical depth to allow for the depth's being 0.
!       Corrections may also be made for cwa quadratic variation in the
!       temperature across the layer and for the effects of edges.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE IR_SOURCE(N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST       &
     &   , SOURCE_COEFF, DEL_PLANCK, L_IR_SOURCE_QUAD, DIFF_PLANCK_2    &
     &   , L_2_STREAM_CORRECT, PLANCK_SOURCE                            &
     &   , GROUND_EMISSION, N_LAYER                                     &
     &   , TAU, TRANS                                                   &
     &   , S_DOWN, S_UP                                                 &
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
#include "scfpt3a.h"
!
!     DUMMY VARIABLES.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , I_LAYER_FIRST                                                &
!             FIRST LAYER TO CONSIDER
     &   , I_LAYER_LAST                                                 &
!             LAST LAYER TO CONSIDER
     &   , N_LAYER
!             NUMBER OF LAYERS
!
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_IR_SOURCE_QUAD                                             &
!             USE A QUADRATIC REPRESENTATION
     &   , L_2_STREAM_CORRECT
!             EDGE CORRECTION TO 2-STREAM
!
      REAL                                                              &
                !, INTENT(IN)
     &     SOURCE_COEFF(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF)       &
!             COEFFICIENTS FOR SOURCE TERMS
     &   , DEL_PLANCK(NPD_PROFILE, NPD_LAYER)                           &
!             DIFFERENCE IN PLANCK FUNCTION ACROSS THE LAYER
     &   , DIFF_PLANCK_2(NPD_PROFILE, NPD_LAYER)                        &
!             2x2ND DIFFERENCE OF PLANCKIAN
     &   , TAU(NPD_PROFILE, NPD_LAYER)                                  &
!             OPTCIAL DEPTH
     &   , TRANS(NPD_PROFILE, NPD_LAYER)                                &
!             TRANSMISSION COEFFICIENT
     &   , PLANCK_SOURCE(NPD_PROFILE, 0: NPD_LAYER)                     &
!             PLANCKIAN SOURCE FUNCTION
     &   , GROUND_EMISSION(NPD_PROFILE)
!             TOTAL FLUX EMITTED FROM GROUND
!
      REAL                                                              &
                !, INTENT(OUT)
     &     S_DOWN(NPD_PROFILE, NPD_LAYER)                               &
!             UPWARD SOURCE FUNCTION
     &   , S_UP(NPD_PROFILE, NPD_LAYER)
!             UPWARD SOURCE FUNCTION
!
!
!     LOCAL VARIABLES.
!
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
      REAL                                                              &
     &     TAUC(NPD_PROFILE, 0: NPD_LAYER)                              &
!             CUMULATIVE OPTICAL DEPTH
     &   , PLANCK_AVE(NPD_PROFILE, 0: NPD_LAYER)                        &
!             AVERAGE PLANCKIAN
     &   , DELTA_TAU_UP_TOP                                             &
!             OPTICAL DEPTH: SURF-TOP OF LAYER
     &   , DELTA_TAU_UP_BASE
!             OPTICAL DEPTH: SURF-BASE OF LAYER
!
!     FUNCTIONS CALLED:
      REAL                                                              &
     &     E3_ACC01
!             THIRD EXPONENTIAL INTEGRAL TO 1%
      EXTERNAL                                                          &
     &     E3_ACC01
!
!
!
!     MULTIPLY THE SOURCE COEFFICIENTS BY THE PLANCKIAN DIFFERENCES
!     TO THE ORDER REQUIRED.
!
      IF (L_IR_SOURCE_QUAD) THEN
!
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
               S_UP(L, I)=SOURCE_COEFF(L, I, IP_SCF_IR_1D)              &
     &            *DEL_PLANCK(L, I)                                     &
     &            +SOURCE_COEFF(L, I, IP_SCF_IR_2D)                     &
     &            *DIFF_PLANCK_2(L, I)
               S_DOWN(L, I)=-SOURCE_COEFF(L, I, IP_SCF_IR_1D)           &
     &            *DEL_PLANCK(L, I)                                     &
     &            +SOURCE_COEFF(L, I, IP_SCF_IR_2D)                     &
     &            *DIFF_PLANCK_2(L, I)
            ENDDO
!
         ENDDO
!
      ELSE
!
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
               S_UP(L, I)=SOURCE_COEFF(L, I, IP_SCF_IR_1D)              &
     &            *DEL_PLANCK(L, I)
               S_DOWN(L, I)=-S_UP(L, I)
            ENDDO
         ENDDO
!
      ENDIF
!
!
!     EDGE CORRECTIONS TO 2-STREAM EQUATIONS.
!
      IF (L_2_STREAM_CORRECT) THEN
!
         DO L=1, N_PROFILE
            TAUC(L, 0)=0.0E+00
         ENDDO
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               TAUC(L, I)=TAUC(L, I-1)+TAU(L, I)
               PLANCK_AVE(L, I)                                         &
     &            =0.5E+00*(PLANCK_SOURCE(L, I-1)+PLANCK_SOURCE(L, I))
            ENDDO
         ENDDO
!
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               DELTA_TAU_UP_TOP=TAUC(L, N_LAYER)-TAUC(L, I-1)
               DELTA_TAU_UP_BASE=TAUC(L, N_LAYER)-TAUC(L, I)
               S_UP(L, I)=S_UP(L, I)                                    &
     &            +2.0E+00*(GROUND_EMISSION(L)-PLANCK_AVE(L, I))        &
! DEPENDS ON: e3_acc01
     &            *(E3_ACC01(DELTA_TAU_UP_TOP)                          &
! DEPENDS ON: e3_acc01
     &            -TRANS(L, I)*E3_ACC01(DELTA_TAU_UP_BASE))
               S_DOWN(L, I)=S_DOWN(L, I)                                &
     &            +2.0E+00*PLANCK_AVE(L, I)                             &
! DEPENDS ON: e3_acc01
     &            *(TRANS(L, I)*E3_ACC01(TAUC(L, I-1))                  &
! DEPENDS ON: e3_acc01
     &            -E3_ACC01(TAUC(L, I)))
            ENDDO
         ENDDO

      ENDIF
!
!
!
      RETURN
      END SUBROUTINE IR_SOURCE
#endif
#endif
