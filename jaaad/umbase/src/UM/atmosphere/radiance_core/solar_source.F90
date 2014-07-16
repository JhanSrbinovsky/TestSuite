#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the solar flux and source terms.
!
! Method:
!       Straightforward.
!
! Current owner of code: James Manners
!
! Description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SOLAR_SOURCE(N_PROFILE, N_LAYER                        &
     &   , FLUX_INC_DIRECT                                              &
     &   , TRANS_0, SOURCE_COEFF                                        &
     &   , L_SCALE_SOLAR, ADJUST_SOLAR_KE                               &
     &   , FLUX_DIRECT                                                  &
     &   , S_DOWN, S_UP                                                 &
     &   , ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF                        &
     &   )
!
!
      USE solinc_data, ONLY: lg_orog_corr, L_orog
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
!           Size allocated for coefficients in the source terms
!
!     Include header files.
#include "c_kinds.h"
#include "source_coeff_pointer_pcf3z.h"
!
!     Dummy variables.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER
!           Number of layers
!
      LOGICAL, INTENT(IN) ::                                            &
     &    L_SCALE_SOLAR
!           Scaling applied to solar beam
!
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FLUX_INC_DIRECT(ND_PROFILE)                                   &
!           Incident solar flux
     &  , TRANS_0(ND_PROFILE, ND_LAYER)                                 &
!           Direct transmission coefficient
     &  , SOURCE_COEFF(ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF)           &
!           Reflection coefficient
     &  , ADJUST_SOLAR_KE(ND_PROFILE, ND_LAYER)
!           Adjustment to solar flux
!
!
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_DIRECT(ND_PROFILE, 0: ND_LAYER)                          &
!           Direct flux
     &  , S_DOWN(ND_PROFILE, ND_LAYER)                                  &
!           Downward source function
     &  , S_UP(ND_PROFILE, ND_LAYER)
!           Upward source function
!
!
!     Local variables.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , L
!           Loop variable
!
!
!
      DO L=1, N_PROFILE
        FLUX_DIRECT(L, 0)=FLUX_INC_DIRECT(L)
      ENDDO
!
!     The solar flux may be multiplied by a scaling factor if an
!     equivalent extinction is used.
!
      IF (L_SCALE_SOLAR) THEN
!
        DO I=1, N_LAYER
          DO L=1, N_PROFILE
            FLUX_DIRECT(L, I)                                           &
     &        =FLUX_DIRECT(L, I-1)*TRANS_0(L, I)                        &
     &        *ADJUST_SOLAR_KE(L, I)
            S_UP(L, I)=SOURCE_COEFF(L, I, IP_SCF_SOLAR_UP)              &
     &        *FLUX_DIRECT(L, I-1)
            S_DOWN(L, I)=(SOURCE_COEFF(L, I, IP_SCF_SOLAR_DOWN)         &
     &        -TRANS_0(L, I))*FLUX_DIRECT(L, I-1)                       &
     &        +FLUX_DIRECT(L, I)
          ENDDO
        ENDDO
!
      ELSE
!
        DO I=1, N_LAYER
          DO L=1, N_PROFILE
            FLUX_DIRECT(L, I)                                           &
     &        =FLUX_DIRECT(L, I-1)*TRANS_0(L, I)
            S_UP(L, I)=SOURCE_COEFF(L, I, IP_SCF_SOLAR_UP)              &
     &        *FLUX_DIRECT(L, I-1)
            S_DOWN(L, I)=SOURCE_COEFF(L, I, IP_SCF_SOLAR_DOWN)          &
     &        *FLUX_DIRECT(L, I-1)
          ENDDO
        ENDDO
!
      ENDIF
!
!
!     Correct the direct flux at the ground for sloping terrain

      IF (L_orog) THEN
         FLUX_DIRECT(1:N_PROFILE, N_LAYER) =                            &
     &      FLUX_DIRECT(1:N_PROFILE, N_LAYER) *                         &
     &      lg_orog_corr(1:N_PROFILE)

         S_DOWN(1:N_PROFILE, N_LAYER) =                                 &
     &         S_DOWN(1:N_PROFILE, N_LAYER) +                           &
     &         FLUX_DIRECT(1:N_PROFILE, N_LAYER) *                      &
     &         (lg_orog_corr(1:N_PROFILE) - 1.0)
      ENDIF

      RETURN
      END SUBROUTINE SOLAR_SOURCE
#endif
