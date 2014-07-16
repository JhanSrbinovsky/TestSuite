#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to rescale the phase function.
!
! Method:
!       The standard rescaling of the phase function is used.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE RESCALE_PHASE_FNC(N_PROFILE                            &
     &  , I_LAYER_FIRST, I_LAYER_LAST, N_DIRECTION, COS_SOL_VIEW        &
     &  , N_ORDER_PHASE, PHASE_FNC, FORWARD_SCATTER, FORWARD_SOLAR      &
     &  , L_RESCALE_SOLAR_PHF, N_ORDER_PHASE_SOLAR, PHASE_FNC_SOLAR     &
     &  , ND_PROFILE, ND_RADIANCE_PROFILE, ND_LAYER, ID_1               &
     &  , ND_DIRECTION, ND_MAX_ORDER                                    &
     &  )
!
!
      IMPLICIT NONE
!
! Include Header Files
#include "c_kinds.h"
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for atmospheric profiles
     &  , ND_RADIANCE_PROFILE                                           &
!           Size allocated for atmospheric profiles used specifically
!           for calculating radiances
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
     &  , ID_1                                                          &
!           Topmost declared layer for optical properties
     &  , ND_MAX_ORDER                                                  &
!           Size allocated for orders of spherical harmonics
     &  , ND_DIRECTION
!           Size allocated for viewing directions
!
!     Dummy arguments
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , I_LAYER_FIRST                                                 &
!           First layer to rescale
     &  , I_LAYER_LAST                                                  &
!           Last layer to rescale
     &  , N_DIRECTION
!           Number of directions
      LOGICAL, INTENT(IN) ::                                            &
     &    L_RESCALE_SOLAR_PHF
!           Flag to rescale the singly scattered solar phase function
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FORWARD_SCATTER(ND_PROFILE, ID_1: ND_LAYER)                   &
!           Forward scattering
     &  , FORWARD_SOLAR(ND_PROFILE, ID_1: ND_LAYER)                     &
!           Forward scattering for the solar beam
     &  , COS_SOL_VIEW(ND_RADIANCE_PROFILE, ND_DIRECTION)
!           Cosines of the angles between the solar direction
!           and the viewing directions
      INTEGER, INTENT(IN) ::                                            &
     &    N_ORDER_PHASE                                                 &
!           Order of terms in the phase function to be retained
     &  , N_ORDER_PHASE_SOLAR
!           Order of terms retained in treating singly scattered
!           solar radiation
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    PHASE_FNC(ND_PROFILE, ID_1: ND_LAYER, ND_MAX_ORDER)           &
!           Phase function
     &  , PHASE_FNC_SOLAR(ND_RADIANCE_PROFILE, ID_1: ND_LAYER           &
     &      , ND_DIRECTION)
!           The phase function for single scattered solar radiation
!
!     Local variables
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , K                                                             &
!           Loop variable
     &  , L                                                             &
!           Loop variable
     &  , ID                                                            &
!           Loop variable
     &  , LS
!           Loop variable
!
!     Legendre polynomials:
      REAL  (Real64) ::                                                 &
     &    CNST1                                                         &
!           Constant in recurrence for Legendre polynomials
     &  , P_LEGENDRE_LS(ND_RADIANCE_PROFILE)                            &
!           Legendre polynomial at the current order
     &  , P_LEGENDRE_LS_M1(ND_RADIANCE_PROFILE)                         &
!           Legendre polynomial at the previous order
     &  , P_LEGENDRE_TMP(ND_RADIANCE_PROFILE)
!           Temporary Legendre polynomial
!
      REAL  (Real64) ::                                                 &
     &    PEAK(ND_PROFILE)
!           Forward scattering peak
!
!
!
      DO K=1, N_ORDER_PHASE
        DO I=I_LAYER_FIRST, I_LAYER_LAST
          DO L=1, N_PROFILE
            PHASE_FNC(L, I, K)                                          &
     &        =(PHASE_FNC(L, I, K)-FORWARD_SCATTER(L, I))               &
     &        /(1.0E+00_Real64-FORWARD_SCATTER(L, I))
          ENDDO
        ENDDO
      ENDDO
!
!
      IF (L_RESCALE_SOLAR_PHF) THEN
!
        DO ID=1, N_DIRECTION
!
!         As usual we do not store Legendre polynomials:
          DO L=1, N_PROFILE
            P_LEGENDRE_LS_M1(L)=1.0E+00_Real64
            P_LEGENDRE_LS(L)=COS_SOL_VIEW(L, ID)
            PEAK(L)=1.0E+00_Real64+P_LEGENDRE_LS(L)*REAL(2*1+1, Real64)
          ENDDO
!
          DO LS=2, N_ORDER_PHASE_SOLAR
!           Calculate higher orders by recurrences.
            CNST1=1.0E+00_Real64-1.0E+00_Real64/REAL(LS, Real64)
            DO L=1, N_PROFILE
              P_LEGENDRE_TMP(L)=P_LEGENDRE_LS(L)
              P_LEGENDRE_LS(L)                                          &
     &          =(1.0E+00_Real64+CNST1)*P_LEGENDRE_LS(L)                &
     &          *COS_SOL_VIEW(L, ID)-CNST1*P_LEGENDRE_LS_M1(L)
              P_LEGENDRE_LS_M1(L)=P_LEGENDRE_TMP(L)
              PEAK(L)=PEAK(L)+REAL(2*LS+1, Real64)*P_LEGENDRE_LS(L)
            ENDDO
          ENDDO
!
!         This is not precisely a rescaling because we do not conserve
!         the forward peak, but what is calculated is what contributes
!         to scattered radiation outside the aureole.
          DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
              PHASE_FNC_SOLAR(L, I, ID)=(PHASE_FNC_SOLAR(L, I, ID)      &
     &          -FORWARD_SOLAR(L, I)*PEAK(L))                           &
     &          /(1.0E+00_Real64-FORWARD_SCATTER(L, I))
            ENDDO
          ENDDO
!
        ENDDO
!

      ENDIF
!
!
!
      RETURN
      END SUBROUTINE RESCALE_PHASE_FNC
#endif
