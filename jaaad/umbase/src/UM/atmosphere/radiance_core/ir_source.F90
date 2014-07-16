#if defined(A70_1Z)
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
!       A correction may also be made for a quadratic variation in the
!       temperature across the layer.
!
! Current owner of code: James Manners
!
! Description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE IR_SOURCE(N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST       &
     &   , SOURCE_COEFF, DEL_PLANCK, L_IR_SOURCE_QUAD, DIFF_PLANCK_2    &
     &   , S_DOWN, S_UP                                                 &
     &   , ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF                        &
     &   )
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
     &  , ND_SOURCE_COEFF
!           Size allocated for source coefficients
!
!     Include header files.
#include "c_kinds.h"
#include "source_coeff_pointer_pcf3z.h"
!
!     Dummy variables.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , I_LAYER_FIRST                                                 &
!           First layer to consider
     &  , I_LAYER_LAST
!           Last layer to consider
!
      LOGICAL, INTENT(IN) ::                                            &
     &    L_IR_SOURCE_QUAD
!           Use a quadratic representation
!
      REAL  (Real64), INTENT(IN) ::                                     &
     &    SOURCE_COEFF(ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF)           &
!           Coefficients for source terms
     &  , DEL_PLANCK(ND_PROFILE, ND_LAYER)                              &
!           Difference in Planckian function across the layer
     &  , DIFF_PLANCK_2(ND_PROFILE, ND_LAYER)
!             2x2nd difference of Planckian
!
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    S_DOWN(ND_PROFILE, ND_LAYER)                                  &
!           Upward source function
     &  , S_UP(ND_PROFILE, ND_LAYER)
!           Upward source function
!
!
!     Local variables.
!
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , L
!           Loop variable
!
!
!
!     Multiply the source coefficients by the Planckian differences
!     to the order required.
!
      IF (L_IR_SOURCE_QUAD) THEN
!
        DO I=I_LAYER_FIRST, I_LAYER_LAST
          DO L=1, N_PROFILE
            S_UP(L, I)=SOURCE_COEFF(L, I, IP_SCF_IR_1D)                 &
     &        *DEL_PLANCK(L, I)                                         &
     &        +SOURCE_COEFF(L, I, IP_SCF_IR_2D)                         &
     &        *DIFF_PLANCK_2(L, I)
            S_DOWN(L, I)=-SOURCE_COEFF(L, I, IP_SCF_IR_1D)              &
     &        *DEL_PLANCK(L, I)                                         &
     &        +SOURCE_COEFF(L, I, IP_SCF_IR_2D)                         &
     &        *DIFF_PLANCK_2(L, I)
          ENDDO
!
        ENDDO
!
      ELSE
!
        DO I=I_LAYER_FIRST, I_LAYER_LAST
          DO L=1, N_PROFILE
            S_UP(L, I)=SOURCE_COEFF(L, I, IP_SCF_IR_1D)                 &
     &        *DEL_PLANCK(L, I)
            S_DOWN(L, I)=-S_UP(L, I)
          ENDDO
        ENDDO
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE IR_SOURCE
#endif
