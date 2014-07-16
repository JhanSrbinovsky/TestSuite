#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the optical properties of aerosols.
!
! Method:
!       If the optical properties come from an observational
!       distribution a separate subroutine is called. Otherwise
!       appropriate mean quantities in the layer are calculated
!       as the parametrization requires and these values are
!       substituted into the parametrization to give the optical
!       properties. Aerosol properties may depend on the humidity.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE OPT_PROP_AEROSOL(IERR                                  &
     &  , N_PROFILE, FIRST_LAYER, LAST_LAYER                            &
     &  , N_ORDER_PHASE, L_RESCALE, N_ORDER_FORWARD                     &
     &  , L_HENYEY_GREENSTEIN_PF                                        &
     &  , N_AEROSOL, AEROSOL_MIX_RATIO                                  &
     &  , I_AEROSOL_PARAMETRIZATION                                     &
     &  , I_HUMIDITY_POINTER, HUMIDITIES, DELTA_HUMIDITY                &
     &  , MEAN_REL_HUMIDITY                                             &
     &  , AEROSOL_ABSORPTION, AEROSOL_SCATTERING, AEROSOL_PHASE_FNC     &
     &  , L_SOLAR_PHF, N_ORDER_PHASE_SOLAR, N_DIRECTION, COS_SOL_VIEW   &
#if defined(OBSERVED)
     &  , P, DENSITY                                                    &
     &  , N_OPT_LEVEL_AEROSOL_PRSC, AEROSOL_PRESSURE_PRSC               &
     &  , AEROSOL_ABSORPTION_PRSC, AEROSOL_SCATTERING_PRSC              &
     &  , AEROSOL_PHASE_FNC_PRSC                                        &
#endif
     &  , K_EXT_TOT, K_EXT_SCAT, PHASE_FNC, FORWARD_SCATTER             &
     &  , FORWARD_SOLAR, PHASE_FNC_SOLAR                                &
     &  , ND_PROFILE, ND_RADIANCE_PROFILE, ND_LAYER, ID_LT, ID_LB       &
     &  , ND_AEROSOL_SPECIES, ND_HUMIDITIES                             &
     &  , ND_PHASE_TERM, ND_MAX_ORDER, ND_DIRECTION                     &
#if defined(OBSERVED)
     &  , ND_PROFILE_PRSC, ND_OPT_LEVEL_PRSC                            &
#endif
     &  )
!
!
!
      IMPLICIT NONE
!
!
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for profiles
     &  , ND_RADIANCE_PROFILE                                           &
!           Size allocated for profiles of quantities only
!           required when radiances are wanted
     &  , ND_LAYER                                                      &
!           Size allocated for layers
     &  , ID_LT                                                         &
!           Topmost declared layer for output optical properties
     &  , ID_LB                                                         &
!           Bottom declared layer for output optical properties
     &  , ND_PHASE_TERM                                                 &
!           Size allocated for terms in phase function
     &  , ND_MAX_ORDER                                                  &
!           Size allocated for orders of sperical harmonics
     &  , ND_Direction                                                  &
!           Size allocated for viewing directions
     &  , ND_AEROSOL_SPECIES                                            &
!           Size allocated for aerosol species
     &  , ND_HUMIDITIES
!           Size allocated for humidities
#if defined(OBSERVED)
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE_PRSC                                               &
!           Size allowed for profiles of prescribed properties
     &  , ND_OPT_LEVEL_PRSC
!           Size allowed for levels of prescribed properties
#endif
!
!     Inclusion of header files.
#include "c_kinds.h"
#include "def_std_io_icf3z.h"
#include "aerosol_parametrization_pcf3z.h"
#include "error_pcf3z.h"
!
!     Dummy variables.
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , FIRST_LAYER                                                   &
!           First layer where propeties are required
     &  , LAST_LAYER                                                    &
!           Last layer where propeties are required
     &  , N_ORDER_PHASE                                                 &
!           Number of terms to retain in the phase function
     &  , N_ORDER_PHASE_SOLAR                                           &
!           Number of terms to retain in calculating the angular
!           scattering of solar radiation
     &  , N_DIRECTION                                                   &
!           Number of viewing directions
     &  , N_ORDER_FORWARD
!           Order used in forming the forward scattering parameter
!
      LOGICAL, INTENT(IN) ::                                            &
     &    L_RESCALE                                                     &
!           Flag for delta-rescaling
     &  , L_HENYEY_GREENSTEIN_PF                                        &
!           Flag to use a Henyey-Greenstein phase function
     &  , L_SOLAR_PHF
!           Flag to use calculate a separate solar phase function
!
      INTEGER, INTENT(IN) ::                                            &
     &    N_AEROSOL                                                     &
!           Number of aerosol species
     &  , I_AEROSOL_PARAMETRIZATION(ND_AEROSOL_SPECIES)                 &
!           Parametrizations of aerosols
     &  , I_HUMIDITY_POINTER(ND_PROFILE,  ND_LAYER)
!           Pointer to aerosol look-up table
!
      REAL  (Real64), INTENT(IN) ::                                     &
     &    COS_SOL_VIEW(ND_RADIANCE_PROFILE, ND_DIRECTION)
!           Cosines of the angles between the solar direction
!           and the viewing direction
!
      REAL  (Real64), INTENT(IN) ::                                     &
     &    AEROSOL_MIX_RATIO(ND_PROFILE, ND_LAYER                        &
     &      , ND_AEROSOL_SPECIES)                                       &
!           Number densty of aerosols
     &  , AEROSOL_ABSORPTION(ND_HUMIDITIES, ND_AEROSOL_SPECIES)         &
!           Aerosol absorption in band/mix ratio
     &  , AEROSOL_SCATTERING(ND_HUMIDITIES, ND_AEROSOL_SPECIES)         &
!           Aerosol scattering in band/mix ratio
     &  , AEROSOL_PHASE_FNC(ND_HUMIDITIES                               &
     &      , ND_PHASE_TERM, ND_AEROSOL_SPECIES)                        &
!           Aerosol phase function in band
     &  , HUMIDITIES(ND_HUMIDITIES, ND_AEROSOL_SPECIES)                 &
!           Array of humidities
     &  , DELTA_HUMIDITY                                                &
!           Increment in humidity
     &  , MEAN_REL_HUMIDITY(ND_PROFILE, ND_LAYER)
!           Mixing ratio of water vapour
!
#if defined(OBSERVED)
      INTEGER, INTENT(IN) ::                                            &
     &    N_OPT_LEVEL_AEROSOL_PRSC(ND_AEROSOL_SPECIES)
!           Number of aerosol data layers
!
      REAL  (Real64), INTENT(IN) ::                                     &
     &    P(ND_PROFILE, ND_LAYER)                                       &
!           Pressure
     &  , DENSITY(ND_PROFILE, ND_LAYER)                                 &
!           Atmospheric density
     &  , AEROSOL_PRESSURE_PRSC(ND_PROFILE_PRSC, ND_OPT_LEVEL_PRSC      &
     &      , ND_AEROSOL_SPECIES)                                       &
!           Pressures at levels of prescribed aerosol properties
     &  , AEROSOL_ABSORPTION_PRSC(ND_PROFILE_PRSC, ND_OPT_LEVEL_PRSC    &
     &      , ND_AEROSOL_SPECIES)                                       &
!           Prescribed aerosol absorption
     &  , AEROSOL_SCATTERING_PRSC(ND_PROFILE_PRSC, ND_OPT_LEVEL_PRSC    &
     &      , ND_AEROSOL_SPECIES)                                       &
!           Prescribed aerosol scattering
     &  , AEROSOL_PHASE_FNC_PRSC(ND_PROFILE_PRSC, ND_OPT_LEVEL_PRSC     &
     &      , ND_PHASE_TERM, ND_AEROSOL_SPECIES)
!           Prescribed aerosol phase function
#endif
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    K_EXT_SCAT(ND_PROFILE, ID_LT: ID_LB)                          &
!           Scattering extinction
     &  , K_EXT_TOT(ND_PROFILE, ID_LT: ID_LB)                           &
!           Total extinction
     &  , PHASE_FNC(ND_PROFILE, ID_LT: ID_LB, ND_MAX_ORDER)             &
!           Phase function
     &  , FORWARD_SCATTER(ND_PROFILE, ID_LT: ID_LB)                     &
!           Forward scattering
     &  , FORWARD_SOLAR(ND_PROFILE, ID_LT: ID_LB)                       &
!           Forward scattering for the solar beam
     &  , PHASE_FNC_SOLAR(ND_RADIANCE_PROFILE                           &
     &      , ID_LT: ID_LB, ND_DIRECTION)
!           Phase function relative to the solar beam
!
!     Local variables.
      LOGICAL                                                           &
     &    L_USE_HG_PHF                                                  &
!           Flag to use Henyey-Greenstein phase functions
     &  , L_INTERPOLATE_HUM
!           Flag to interpolate optical properties through a look-up
!           table of humidities
      INTEGER                                                           &
     &    L                                                             &
!           Loop variable
     &  , J                                                             &
!           Loop variable
     &  , I                                                             &
!           Loop variable
     &  , ID                                                            &
!           Loop variable
     &  , LS                                                            &
!           Loop variable
     &  , I_POINTER
!           Temporary pointer
      REAL  (Real64) ::                                                 &
     &    K_SCATTER(ND_PROFILE)                                         &
!           Scattering of current extinction of the current aerosol
     &  , ASYMMETRY(ND_PROFILE)                                         &
!           Asymmetry of the current aerosol
     &  , KS_PHF(ND_PROFILE)                                            &
!           Scattering coefficient multiplied by a coefficient in the
!           phase function
     &  , PHF_COEFF(ND_PROFILE)
!           Coefficient in the phase function of the current aerosol
!
      REAL  (Real64) ::                                                 &
     &    WEIGHT_UPPER(ND_PROFILE)
!           Weighting towards the upper end of an interval
!           in a look-up table
!
#if defined(OBSERVED)
!     Optical properties interpolated from prescribed properties
!     for each component
      REAL  (Real64) ::                                                 &
     &    K_EXT_SCAT_COMP(ND_PROFILE, ID_LT: ID_LB)                     &
!           Scattering extinction of component
     &  , K_EXT_TOT_COMP(ND_PROFILE, ID_LT: ID_LB)                      &
!           Total extinction of component
     &  , PHASE_FNC_COMP(ND_PROFILE, ID_LT: ID_LB, ND_MAX_ORDER)        &
!           Phase function of component
     &  , PHASE_FNC_SOLAR_COMP(ND_RADIANCE_PROFILE                      &
     &      , ID_LT: ID_LB, ND_DIRECTION)                               &
!           Solar phase function of component
     &  , FORWARD_SCATTER_COMP(ND_PROFILE, ID_LT: ID_LB)                &
!           Forward scattering of component
     &  , FORWARD_SOLAR_COMP(ND_PROFILE, ID_LT: ID_LB)
!           Forward scattering of the solar beam for the component
#endif
!
!     Legendre polynomials:
      REAL  (Real64) ::                                                 &
     &    CNST1                                                         &
!           Constant in recurrence for Legendre polynomials
     &  , P_LEGENDRE_LS(ND_RADIANCE_PROFILE)                            &
!           Legendre polynomial at the current order
     &  , P_LEGENDRE_LS_M1(ND_RADIANCE_PROFILE)                         &
!           Legendre polynomial at the previous order
     &  , P_LEGENDRE_TMP(ND_RADIANCE_PROFILE)                           &
!           Temporary Legendre polynomial
     &  , PHASE_FNC_SOLAR_TMP(ND_RADIANCE_PROFILE)
!           Current contribution to the solar phase function
!
#if defined(OBSERVED)
!     Subroutines called:
      EXTERNAL                                                          &
     &   PRSC_OPT_PROP
#endif
!
!
!
      DO J=1, N_AEROSOL
!
!       Use the Henyey-Greenstein phase function if specifically
!       requested, or if using an old parametrization which gives
!       only an asymmetry.
        L_USE_HG_PHF=L_HENYEY_GREENSTEIN_PF.OR.                         &
     &    (I_AEROSOL_PARAMETRIZATION(J) == IP_AEROSOL_PARAM_DRY).OR.    &
     &    (I_AEROSOL_PARAMETRIZATION(J) == IP_AEROSOL_PARAM_MOIST)
!
!       Interpolate from the look-up table if using moist properties.
        L_INTERPOLATE_HUM=                                              &
     &    (I_AEROSOL_PARAMETRIZATION(J) == IP_AEROSOL_PARAM_MOIST).OR.  &
     &    (I_AEROSOL_PARAMETRIZATION(J) == IP_AEROSOL_PARAM_PHF_MOIST)
!
        IF ( (I_AEROSOL_PARAMETRIZATION(J) ==                           &
     &        IP_AEROSOL_PARAM_DRY).OR.                                 &
     &       (I_AEROSOL_PARAMETRIZATION(J) ==                           &
     &        IP_AEROSOL_PARAM_PHF_DRY).OR.                             &
     &       (I_AEROSOL_PARAMETRIZATION(J) ==                           &
     &        IP_AEROSOL_PARAM_MOIST).OR.                               &
     &       (I_AEROSOL_PARAMETRIZATION(J) ==                           &
     &        IP_AEROSOL_PARAM_PHF_MOIST) ) THEN
!
!
          DO I=FIRST_LAYER, LAST_LAYER

            IF (L_INTERPOLATE_HUM) THEN
!
!             Calculate the required weights for interpolation
!             in this layer.
              DO L=1, N_PROFILE
                I_POINTER=I_HUMIDITY_POINTER(L, I)
                WEIGHT_UPPER(L)=(MEAN_REL_HUMIDITY(L, I)                &
     &            -HUMIDITIES(I_POINTER, J))                            &
     &            /DELTA_HUMIDITY
              ENDDO
!
!             Interpolate the absorption and scattering.
              DO L=1, N_PROFILE
                I_POINTER=I_HUMIDITY_POINTER(L, I)
                K_EXT_TOT(L, I)=K_EXT_TOT(L, I)                         &
     &            +AEROSOL_MIX_RATIO(L, I, J)                           &
     &            *(AEROSOL_ABSORPTION(I_POINTER, J)                    &
     &            +WEIGHT_UPPER(L)                                      &
     &            *(AEROSOL_ABSORPTION(I_POINTER+1, J)                  &
     &            -AEROSOL_ABSORPTION(I_POINTER, J)))
                K_SCATTER(L)=AEROSOL_MIX_RATIO(L, I, J)                 &
     &            *(AEROSOL_SCATTERING(I_POINTER, J)                    &
     &            +WEIGHT_UPPER(L)                                      &
     &            *(AEROSOL_SCATTERING(I_POINTER+1, J)                  &
     &            -AEROSOL_SCATTERING(I_POINTER, J)))
                K_EXT_SCAT(L, I)=K_EXT_SCAT(L, I)                       &
     &            +K_SCATTER(L)
              ENDDO
!
            ELSE
!
              DO L=1, N_PROFILE
!
!               Calculate volume extinctions directly from the
!               mass extinction coefficients.
                K_EXT_TOT(L, I)=K_EXT_TOT(L, I)                         &
     &            +AEROSOL_MIX_RATIO(L, I, J)                           &
     &            *AEROSOL_ABSORPTION(1, J)
                K_SCATTER(L)=AEROSOL_MIX_RATIO(L, I, J)                 &
     &            *AEROSOL_SCATTERING(1, J)
                K_EXT_SCAT(L, I)=K_EXT_SCAT(L, I)                       &
     &            +K_SCATTER(L)
!
              ENDDO
!
            ENDIF
!
!           The phase function:
!
            IF (L_USE_HG_PHF) THEN
!
!             Note that there is an ambiguity in the definition of a
!             Henyey-Greenstein phase function when humidity is included
!             since one could set up the lookup table externally with
!             all moments at the reference points set to powers of the
!             appropriate asymmetries, but then linear interpolation in
!             the humidity would not give a true Henyey-Greenstein
!             phase function at intermediate points. Here we adopt a
!             true Henyey-Greenstein approach, calculating just the
!             asymmetry.
!
!             Calculate the asymmetry:
              IF (L_INTERPOLATE_HUM) THEN
                DO L=1, N_PROFILE
                  I_POINTER=I_HUMIDITY_POINTER(L, I)
                  ASYMMETRY(L)                                          &
     &              =AEROSOL_PHASE_FNC(I_POINTER, 1, J)                 &
     &              +WEIGHT_UPPER(L)                                    &
     &              *(AEROSOL_PHASE_FNC(I_POINTER+1, 1, J)              &
     &              -AEROSOL_PHASE_FNC(I_POINTER, 1, J))
                ENDDO
              ELSE
                DO L=1, N_PROFILE
                  ASYMMETRY(L)=AEROSOL_PHASE_FNC(1, 1, J)
                ENDDO
              ENDIF
!
!             Set the lowest order in the phase function (required
!             for two-stream calculations and other quadratures).
              DO L=1, N_PROFILE
                PHASE_FNC(L, I, 1)=PHASE_FNC(L, I, 1)                   &
     &            +K_SCATTER(L)*ASYMMETRY(L)
              ENDDO
!
!             Initialize the product of the scattering and the
!             current moment of the phase function. This repeats
!             part of the preceeding loop, but separating it saves
!             an assignment in the case of two-stream calculations.
              IF (L_RESCALE.OR.(N_ORDER_PHASE >= 2)) THEN
                DO L=1, N_PROFILE
                  KS_PHF(L)=K_SCATTER(L)*ASYMMETRY(L)
                ENDDO
              ENDIF
!
!             Calculate weighted higher moments recursively.
              DO LS=2, N_ORDER_PHASE
                DO L=1, N_PROFILE
                  KS_PHF(L)=KS_PHF(L)*ASYMMETRY(L)
                  PHASE_FNC(L, I, LS)                                   &
     &              =PHASE_FNC(L, I, LS)+KS_PHF(L)
                ENDDO
              ENDDO
!
!             Usually, we will retain terms as far as the order of
!             truncation, but in the case of two-stream methods the
!             order of truncation will exceed the order of retention
!             by 1.
              IF (L_RESCALE) THEN
!
                IF (N_ORDER_FORWARD == N_ORDER_PHASE) THEN
                  DO L=1, N_PROFILE
                    FORWARD_SCATTER(L, I)                               &
     &                =FORWARD_SCATTER(L, I)+KS_PHF(L)
                  ENDDO
                ELSE IF (N_ORDER_FORWARD == N_ORDER_PHASE+1) THEN
                  DO L=1, N_PROFILE
                    FORWARD_SCATTER(L, I)                               &
     &                =FORWARD_SCATTER(L, I)+KS_PHF(L)*ASYMMETRY(L)
                  ENDDO
                ELSE
!                 This case probably shouldn't arise so we use
!                 inefficient explicit exponentiation.
                  DO L=1, N_PROFILE
                    FORWARD_SCATTER(L, I)                               &
     &                =FORWARD_SCATTER(L, I)                            &
     &                +K_SCATTER(L)*ASYMMETRY(L)**N_ORDER_FORWARD
                  ENDDO
                ENDIF
!
              ENDIF
!

!
            ELSE
!
!             Calculate the phase function generally. We don't
!             separate the first order here, because it is unlikely
!             that this block will be used in the case in a
!             two-stream calculation.
              DO LS=1, N_ORDER_PHASE
                IF (L_INTERPOLATE_HUM) THEN
                  DO L=1, N_PROFILE
                    I_POINTER=I_HUMIDITY_POINTER(L, I)
                    PHF_COEFF(L)=AEROSOL_PHASE_FNC(I_POINTER, LS, J)    &
     &                +WEIGHT_UPPER(L)                                  &
     &                *(AEROSOL_PHASE_FNC(I_POINTER+1, LS, J)           &
     &                -AEROSOL_PHASE_FNC(I_POINTER, LS, J))
                  ENDDO
                ELSE
                  DO L=1, N_PROFILE
                    PHF_COEFF(L)=AEROSOL_PHASE_FNC(1, LS, J)
                  ENDDO
                ENDIF
                DO L=1, N_PROFILE
                  PHASE_FNC(L, I, LS)=PHASE_FNC(L, I, LS)               &
     &              +K_SCATTER(L)*PHF_COEFF(L)
                ENDDO
              ENDDO
!
              IF (L_RESCALE) THEN
                IF (L_INTERPOLATE_HUM) THEN
                  DO L=1, N_PROFILE
                    I_POINTER=I_HUMIDITY_POINTER(L, I)
                    PHF_COEFF(L)                                        &
     &                =AEROSOL_PHASE_FNC(I_POINTER, N_ORDER_FORWARD, J) &
     &                +WEIGHT_UPPER(L)                                  &
     &                *(AEROSOL_PHASE_FNC(I_POINTER+1                   &
     &                , N_ORDER_FORWARD, J)                             &
     &                -AEROSOL_PHASE_FNC(I_POINTER                      &
     &                , N_ORDER_FORWARD, J))
                  ENDDO
                ELSE
                  DO L=1, N_PROFILE
                    PHF_COEFF(L)                                        &
     &                =AEROSOL_PHASE_FNC(1, N_ORDER_FORWARD, J)
                  ENDDO
                ENDIF
                DO L=1, N_PROFILE
                  FORWARD_SCATTER(L, I)                                 &
     &              =FORWARD_SCATTER(L, I)                              &
     &              +K_SCATTER(L)*PHF_COEFF(L)
                ENDDO
              ENDIF
            ENDIF
!
            IF (L_SOLAR_PHF) THEN
!             Calculate the solar phase function to higher accuracy.
              DO ID=1, N_DIRECTION
!               The Legendre polynomials are not stored so as to reduce
!               the requirement for memory at very high orders of solar
!               truncation.
                IF (L_INTERPOLATE_HUM) THEN
                  DO L=1, N_PROFILE
                    I_POINTER=I_HUMIDITY_POINTER(L, I)
                    PHF_COEFF(L)=AEROSOL_PHASE_FNC(I_POINTER, 1, J)     &
     &                +WEIGHT_UPPER(L)                                  &
     &                *(AEROSOL_PHASE_FNC(I_POINTER+1, 1, J)            &
     &                -AEROSOL_PHASE_FNC(I_POINTER, 1, J))
                  ENDDO
                ELSE
                  DO L=1, N_PROFILE
                    PHF_COEFF(L)=AEROSOL_PHASE_FNC(1, 1, J)
                  ENDDO
                ENDIF
                DO L=1, N_PROFILE
!                 Initialize the Legendre polynomials at the zeroth and
!                 first orders.
                  P_LEGENDRE_LS_M1(L)=1.0E+00_Real64
                  P_LEGENDRE_LS(L)=COS_SOL_VIEW(L, ID)
                  PHASE_FNC_SOLAR_TMP(L)=1.0E+00_Real64+PHF_COEFF(L)    &
     &              *P_LEGENDRE_LS(L)*REAL(2*1+1, Real64)
                ENDDO
!
                IF (L_USE_HG_PHF) THEN
!                 Store the asymmetry in this case.
                  DO L=1, N_PROFILE
                    ASYMMETRY(L)=PHF_COEFF(L)
                  ENDDO
                ENDIF
!
                DO LS=2, N_ORDER_PHASE_SOLAR
!                 Calculate higher orders by recurrences.
                  CNST1=1.0E+00_Real64-1.0E+00_Real64/REAL(LS, Real64)
                  DO L=1, N_PROFILE
                    P_LEGENDRE_TMP(L)=P_LEGENDRE_LS(L)
                    P_LEGENDRE_LS(L)                                    &
     &                =(1.0E+00_Real64+CNST1)*P_LEGENDRE_LS(L)          &
     &                *COS_SOL_VIEW(L, ID)-CNST1*P_LEGENDRE_LS_M1(L)
                    P_LEGENDRE_LS_M1(L)=P_LEGENDRE_TMP(L)
                  ENDDO
!
!                 Calculate the next moment of the phase function.
                  IF (L_USE_HG_PHF) THEN
                    DO L=1, N_PROFILE
                      PHF_COEFF(L)=PHF_COEFF(L)*ASYMMETRY(L)
                    ENDDO
                  ELSE
                    IF (L_INTERPOLATE_HUM) THEN
                      DO L=1, N_PROFILE
                        I_POINTER=I_HUMIDITY_POINTER(L, I)
                        PHF_COEFF(L)                                    &
     &                    =AEROSOL_PHASE_FNC(I_POINTER, LS, J)          &
     &                    +WEIGHT_UPPER(L)                              &
     &                    *(AEROSOL_PHASE_FNC(I_POINTER+1, LS, J)       &
     &                    -AEROSOL_PHASE_FNC(I_POINTER, LS, J))
                      ENDDO
                    ELSE
                      DO L=1, N_PROFILE
                        PHF_COEFF(L)=AEROSOL_PHASE_FNC(1, LS, J)
                      ENDDO
                    ENDIF
                  ENDIF
                  DO L=1, N_PROFILE
                    PHASE_FNC_SOLAR_TMP(L)= PHASE_FNC_SOLAR_TMP(L)      &
     &                 + PHF_COEFF(L)                                   &
!     &                =AEROSOL_PHASE_FNC(1, LS, J)
     &                *REAL(2*LS+1, Real64)*P_LEGENDRE_LS(L)
                  ENDDO
                ENDDO
!               Increment the stored phase function.
                DO L=1, N_PROFILE
                  PHASE_FNC_SOLAR(L, I, ID)                             &
     &              =PHASE_FNC_SOLAR(L, I, ID)                          &
     &            +K_SCATTER(L)*PHASE_FNC_SOLAR_TMP(L)
                ENDDO
              ENDDO
!
!             Continue to an extra order to find the rescaling
!             for the solar beam.
              IF (L_RESCALE) THEN
                IF (L_USE_HG_PHF) THEN
                  DO L=1, N_PROFILE
                    PHF_COEFF(L)=PHF_COEFF(L)*ASYMMETRY(L)
                  ENDDO
                ELSE
                  LS=N_ORDER_PHASE_SOLAR+1
                  IF (L_INTERPOLATE_HUM) THEN
                    DO L=1, N_PROFILE
                      I_POINTER=I_HUMIDITY_POINTER(L, I)
                      PHF_COEFF(L)                                      &
     &                  =AEROSOL_PHASE_FNC(I_POINTER, LS, J)            &
     &                  +WEIGHT_UPPER(L)                                &
     &                  *(AEROSOL_PHASE_FNC(I_POINTER+1, LS, J)         &
     &                  -AEROSOL_PHASE_FNC(I_POINTER, LS, J))
                    ENDDO
                  ELSE
                    DO L=1, N_PROFILE
                      PHF_COEFF(L)=AEROSOL_PHASE_FNC(1, LS, J)
                    ENDDO
                  ENDIF
                ENDIF
                DO L=1, N_PROFILE
                  FORWARD_SOLAR(L, I)=FORWARD_SOLAR(L, I)               &
     &              +K_SCATTER(L)*PHF_COEFF(L)
                ENDDO
              ENDIF
!
            ENDIF
!
          ENDDO
!
!
!
#if defined(OBSERVED)
        ELSE IF (I_AEROSOL_PARAMETRIZATION(J) ==                        &
     &           IP_AEROSOL_UNPARAMETRIZED) THEN
!
           CALL PRSC_OPT_PROP(IERR                                      &
     &       , N_PROFILE, FIRST_LAYER, LAST_LAYER                       &
     &       , L_RESCALE, N_ORDER_FORWARD                               &
     &       , L_HENYEY_GREENSTEIN_PF, N_ORDER_PHASE                    &
     &       , P, DENSITY                                               &
     &       , N_OPT_LEVEL_AEROSOL_PRSC(J)                              &
     &       , AEROSOL_PRESSURE_PRSC(1, 1, J)                           &
     &       , AEROSOL_ABSORPTION_PRSC(1, 1, J)                         &
     &       , AEROSOL_SCATTERING_PRSC(1, 1, J)                         &
     &       , AEROSOL_PHASE_FNC_PRSC(1, 1, 1, J)                       &
     &       , K_EXT_TOT_COMP, K_EXT_SCAT_COMP, PHASE_FNC_COMP          &
     &       , FORWARD_SCATTER_COMP, FORWARD_SOLAR_COMP                 &
     &       , L_SOLAR_PHF, N_ORDER_PHASE_SOLAR                         &
     &       , N_DIRECTION, COS_SOL_VIEW                                &
     &       , PHASE_FNC_SOLAR_COMP                                     &
     &       , ND_PROFILE, ND_RADIANCE_PROFILE, ND_LAYER                &
     &       , ID_LT, ID_LB                                             &
     &       , ND_DIRECTION                                             &
     &       , ND_PROFILE_PRSC, ND_OPT_LEVEL_PRSC                       &
     &       , ND_PHASE_TERM, ND_MAX_ORDER                              &
     &       )
!
          DO I=FIRST_LAYER, LAST_LAYER
            DO L=1, N_PROFILE
              K_EXT_TOT(L, I)=K_EXT_TOT(L, I)                           &
     &          +K_EXT_TOT_COMP(L, I)
              K_EXT_SCAT(L, I)=K_EXT_SCAT(L, I)                         &
     &          +K_EXT_SCAT_COMP(L, I)
            ENDDO
          ENDDO
          DO LS=1, N_ORDER_PHASE
            DO I=FIRST_LAYER, LAST_LAYER
              DO L=1, N_PROFILE
                PHASE_FNC(L, I, LS)=PHASE_FNC(L, I, LS)                 &
     &            +PHASE_FNC_COMP(L, I, LS)
              ENDDO
            ENDDO
          ENDDO
          IF (L_RESCALE) THEN
            DO I=FIRST_LAYER, LAST_LAYER
              DO L=1, N_PROFILE
                FORWARD_SCATTER(L, I)=FORWARD_SCATTER(L, I)             &
     &            +FORWARD_SCATTER_COMP(L, I)
              ENDDO
            ENDDO
          ENDIF
          IF (L_SOLAR_PHF) THEN
            IF (L_RESCALE) THEN
              DO I=FIRST_LAYER, LAST_LAYER
                DO L=1, N_PROFILE
                  FORWARD_SOLAR(L, I)=FORWARD_SOLAR(L, I)               &
     &              +FORWARD_SOLAR_COMP(L, I)
                ENDDO
              ENDDO
            ENDIF
            DO ID=1, N_DIRECTION
              DO I=FIRST_LAYER, LAST_LAYER
                DO L=1, N_PROFILE
                  PHASE_FNC_SOLAR(L, I, ID)                             &
     &              =PHASE_FNC_SOLAR(L, I, ID)                          &
     &              +PHASE_FNC_SOLAR_COMP(L, I, ID)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
!
#endif
        ELSE
!
          WRITE(IU_ERR, '(/A, I3, A)')                                  &
     &      '*** Error: I_AEROSOL_PARAMETRIZATION for species '         &
     &      , J, ' has been set to an illegal value.'
          IERR=I_ERR_FATAL
          RETURN
!
        ENDIF
!
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE OPT_PROP_AEROSOL
#endif
