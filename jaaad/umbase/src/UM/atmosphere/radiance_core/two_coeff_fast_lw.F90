#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate coefficients in the two-stream equations.
!
! Method:
!       The basic two-stream coefficients in the differential
!       equations are calculated on the assumption that scattering
!       can be ignored this routine is therefore only suitable
!       for use in the IR region. These coefficients are then used
!       to determine the transmission and reflection coefficients.
!       Coefficients for determining the solar or infra-red source
!       terms are also calculated.
!
! Current owner of code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE TWO_COEFF_FAST_LW(N_PROFILE                            &
     &  , I_LAYER_FIRST, I_LAYER_LAST                                   &
     &  , L_IR_SOURCE_QUAD, TAU                                         &
     &  , TRANS, SOURCE_COEFF                                           &
     &  , ND_PROFILE, ND_LAYER, ID_LT, ID_LB, ND_SOURCE_COEFF           &
     &  )
!
!
!
      IMPLICIT NONE
!
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for atmospheric profiles
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
     &  , ID_LT                                                         &
!           Topmost declared layer of optical depths
     &  , ID_LB                                                         &
!           Bottom declared layer of optical depths
     &  , ND_SOURCE_COEFF
!           Size allocated for source coefficients
!
!     Include header files.
#include "c_kinds.h"
#include "source_coeff_pointer_pcf3z.h"
!
!
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , I_LAYER_FIRST                                                 &
!           First layer to process
     &  , I_LAYER_LAST
!           Last layer to process
      LOGICAL, INTENT(IN) ::                                            &
     &    L_IR_SOURCE_QUAD
!           Use a quadratic source function
!
!     Optical properties of layer:
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TAU(ND_PROFILE, ID_LT: ID_LB)
!           Optical depths
!
!
!     Coefficients in the two-stream equations:
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    TRANS(ND_PROFILE, ND_LAYER)                                   &
!           Diffuse transmission coefficient
     &  , SOURCE_COEFF(ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF)
!           Source coefficients in two-stream equations
!
!
!     Local variables.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , L                                                             &
!           Loop variable
     &  , N_IN
!           Number of elements for vector exponential
!
!     Variables related to the treatment of ill-conditioning
      REAL  (Real64) ::                                                 &
     &    EPS_R                                                         &
!           The smallest real number such that 1.0-EPS_R is not 1
!           to the computer's precision
     &  , SQ_EPS_R
!           The square root of the above
!
!
!
!     Set the tolerances used in avoiding ill-conditioning, testing
!     on any variable.
      EPS_R=EPSILON(TAU(1, ID_LT))
      SQ_EPS_R=SQRT(EPS_R)
!
      DO I=I_LAYER_FIRST, I_LAYER_LAST
        DO L=1, N_PROFILE
          TRANS(L, I)=-1.66E+00_Real64*TAU(L, I)
        ENDDO
        DO L=N_PROFILE+1, ND_PROFILE
          TRANS(L, I)=0.0E+00_Real64
        ENDDO
      ENDDO
      N_IN=ND_PROFILE*(I_LAYER_LAST-I_LAYER_FIRST+1)
! DEPENDS ON: exp_v
      CALL EXP_V(N_IN                                                   &
     &  , TRANS(1, I_LAYER_FIRST), TRANS(1, I_LAYER_FIRST))
!
      DO I=I_LAYER_FIRST, I_LAYER_LAST
        DO L=1, N_PROFILE
          SOURCE_COEFF(L, I, IP_SCF_IR_1D)                              &
     &      =(1.0E+00_Real64-TRANS(L, I)+SQ_EPS_R)                      &
     &      /(1.66E+00_Real64*TAU(L, I)+SQ_EPS_R)
        ENDDO
      ENDDO
!
      IF (L_IR_SOURCE_QUAD) THEN
        DO I=I_LAYER_FIRST, I_LAYER_LAST
          DO L=1, N_PROFILE
            SOURCE_COEFF(L, I, IP_SCF_IR_2D)                            &
     &        =-(1.0E+00_Real64+TRANS(L, I)                             &
     &        -2.0E+00_Real64*SOURCE_COEFF(L, I, IP_SCF_IR_1D))         &
     &        /(1.66E+00_Real64*TAU(L, I)+SQ_EPS_R)
          ENDDO
        ENDDO
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE TWO_COEFF_FAST_LW
#endif
