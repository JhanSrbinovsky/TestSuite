#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate optical properties of water clouds.
!
! Method:
!       If the optical properties come from an observational
!       distribution a separate subroutine is called. Otherwise
!       appropriate mean quantities in the layer are calculated
!       as the parametrization requires and these values are
!       substituted into the parametrization to give the optical
!       properties.
!
!       Note that this routine produces optical propeties for a
!       single condensed component of the cloud.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE OPT_PROP_WATER_CLOUD(IERR                              &
     &  , N_PROFILE, N_LAYER, N_CLOUD_TOP                               &
     &  , N_CLOUD_PROFILE, I_CLOUD_PROFILE                              &
     &  , N_ORDER_PHASE, L_RESCALE, N_ORDER_FORWARD                     &
     &  , L_HENYEY_GREENSTEIN_PF, L_SOLAR_PHF, N_ORDER_PHASE_SOLAR      &
     &  , N_DIRECTION, COS_SOL_VIEW                                     &
     &  , I_PARAMETRIZATION_DROP, CLOUD_PARAMETER                       &
     &  , LIQ_WATER_MASS_FRAC, RADIUS_EFFECT                            &
#if defined(OBSERVED)
     &  , P, DENSITY                                                    &
     &  , N_OPT_LEVEL_CLOUD_PRSC                                        &
     &  , DROP_PRESSURE_PRSC, DROP_ABSORPTION_PRSC                      &
     &  , DROP_SCATTERING_PRSC, DROP_PHASE_FNC_PRSC                     &
#endif
     &  , K_EXT_TOT_CLOUD, K_EXT_SCAT_CLOUD                             &
     &  , PHASE_FNC_CLOUD, FORWARD_SCATTER_CLOUD                        &
     &  , FORWARD_SOLAR_CLOUD, PHASE_FNC_SOLAR_CLOUD                    &
     &  , ND_PROFILE, ND_RADIANCE_PROFILE, ND_LAYER, ID_CT              &
     &  , ND_DIRECTION                                                  &
     &  , ND_PHASE_TERM, ND_MAX_ORDER, ND_CLOUD_PARAMETER               &
#if defined(OBSERVED)
     &  , ND_PROFILE_PRSC, ND_OPT_LEVEL_PRSC                            &
#endif
     &  )
!
!
      IMPLICIT NONE
!
!
!     Sizes of arrays
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for profiles
     &  , ND_RADIANCE_PROFILE                                           &
!           Size allocated for profiles
     &  , ND_LAYER                                                      &
!           Size allocated for layers
     &  , ID_CT                                                         &
!           Topmost declared cloudy layer
     &  , ND_DIRECTION                                                  &
!           Size allocated for viewing directions
     &  , ND_PHASE_TERM                                                 &
!           Size allocated for terms in phase function
     &  , ND_MAX_ORDER                                                  &
!           Size allocated for orders of spherical harmonics
     &  , ND_CLOUD_PARAMETER
!            Size allocated for cloud parameters
#if defined(OBSERVED)
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE_PRSC                                               &
!           Size allowed for profiles of prescribed optical properties
     &  , ND_OPT_LEVEL_PRSC
!           Size allowed for levels of prescribed optical properties
#endif
!
!     Inclusion of header files.
#include "c_kinds.h"
#include "def_std_io_icf3z.h"
#include "cloud_parametrization_pcf3z.h"
#include "error_pcf3z.h"
!
!     Dummy variables.
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , N_CLOUD_TOP                                                   &
!           Topmost cloudy layer
     &  , N_ORDER_PHASE                                                 &
!           Number of terms to retain in the phase function
     &  , N_ORDER_PHASE_SOLAR                                           &
!           Number of terms to retain in single scattered solar
!           phase function
     &  , N_ORDER_FORWARD                                               &
!           Order used in forming the forward scattering parameter
     &  , I_PARAMETRIZATION_DROP                                        &
!           Treatment of droplets
#if defined(OBSERVED)
     &  , N_OPT_LEVEL_CLOUD_PRSC                                        &
!           Number of levels of prescribed optical properties
#endif
     &  , N_CLOUD_PROFILE(ID_CT: ND_LAYER)                              &
!           Number of cloudy profiles
     &  , I_CLOUD_PROFILE(ND_PROFILE, ID_CT: ND_LAYER)
!           Profiles containing clouds
      LOGICAL, INTENT(IN) ::                                            &
     &    L_RESCALE                                                     &
!           Flag for delta-rescaling
     &  , L_HENYEY_GREENSTEIN_PF                                        &
!           Flag to use a Henyey-Greenstein phase function
     &  , L_SOLAR_PHF
!           Flag to use an extended solar phase function in
!           single scattering
!
!     Viewing directions:
      INTEGER, INTENT(IN) ::                                            &
     &    N_DIRECTION
!           Number of viewing dierctions
      REAL  (Real64), INTENT(IN) ::                                     &
     &    COS_SOL_VIEW(ND_RADIANCE_PROFILE, ND_DIRECTION)
!           Cosines of the angles between the solar direction
!           and the viewing direction
!
      REAL  (Real64), INTENT(IN) ::                                     &
     &    CLOUD_PARAMETER(ND_CLOUD_PARAMETER)                           &
!           Cloud parameters
     &  , LIQ_WATER_MASS_FRAC(ND_PROFILE, ID_CT: ND_LAYER)              &
!           Liquid water content
     &  , RADIUS_EFFECT(ND_PROFILE, ID_CT: ND_LAYER)
!           Effective radius
#if defined(OBSERVED)
      REAL  (Real64), INTENT(IN) ::                                     &
     &    P(ND_PROFILE, ND_LAYER)                                       &
!           Pressure
     &  , DENSITY(ND_PROFILE, ND_LAYER)                                 &
!           Atmospheric density
     &  , DROP_PRESSURE_PRSC(ND_PROFILE_PRSC, ND_OPT_LEVEL_PRSC)        &
!           Pressure levels where optical data are prescribed
     &  , DROP_ABSORPTION_PRSC(ND_PROFILE_PRSC, ND_OPT_LEVEL_PRSC)      &
!           Prescribed absorption by droplets
     &  , DROP_SCATTERING_PRSC(ND_PROFILE_PRSC, ND_OPT_LEVEL_PRSC)      &
!           Prescribed scattering by droplets
     &  , DROP_PHASE_FNC_PRSC(ND_PROFILE_PRSC, ND_OPT_LEVEL_PRSC        &
     &      , ND_PHASE_TERM)
!           Prescribed phase function of droplets
#endif
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    K_EXT_SCAT_CLOUD(ND_PROFILE, ID_CT: ND_LAYER)                 &
!           Scattering extinction
     &  , K_EXT_TOT_CLOUD(ND_PROFILE, ID_CT: ND_LAYER)                  &
!           Total extinction
     &  , PHASE_FNC_CLOUD(ND_PROFILE, ID_CT: ND_LAYER, ND_MAX_ORDER)    &
!           Cloudy phase function
     &  , PHASE_FNC_SOLAR_CLOUD(ND_RADIANCE_PROFILE, ID_CT: ND_LAYER    &
     &      , ND_DIRECTION)                                             &
!           Cloudy phase function for singly scattered solar radiation
     &  , FORWARD_SCATTER_CLOUD(ND_PROFILE, ID_CT: ND_LAYER)            &
!           Cloudy forward scattering
     &  , FORWARD_SOLAR_CLOUD(ND_PROFILE, ID_CT: ND_LAYER)
!           Cloudy forward scattering for the solar beam
!
!     Local variables.
      INTEGER                                                           &
     &    L                                                             &
!           Loop variable
     &  , LL                                                            &
!           Loop variable
     &  , I                                                             &
!           Loop variable
     &  , ID                                                            &
!           Loop variable
     &  , LS
!           Loop variable
      REAL  (Real64) ::                                                 &
     &    ASYMMETRY_PROCESS(ND_PROFILE)                                 &
!           Asymmetry of current process.
     &  , PHF_TMP
!           Temporary Phase Function
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
     &  , KS_PHF(ND_RADIANCE_PROFILE)
!           Product of the scattering and the current moment of
!           the phase function
!
#if defined(OBSERVED)
!     Subroutines called:
      EXTERNAL                                                          &
     &    PRSC_OPT_PROP
#endif
!
!
!****************************************
! ADDITIONAL CODE ADDED BY JCT 30/09/05 *
!****************************************
!

      IF ( (I_PARAMETRIZATION_DROP == IP_SLINGO_SCHRECKER).OR.          &
     &     (I_PARAMETRIZATION_DROP == IP_ACKERMAN_STEPHENS).OR.         &
     &     (I_PARAMETRIZATION_DROP == IP_DROP_PADE_2) .OR.              &
     &     ( L_HENYEY_GREENSTEIN_PF .AND.                               &
     &       (I_PARAMETRIZATION_DROP == IP_SLINGO_SCHR_PHF) ) ) THEN
!
!       Optical properties are calculated from parametrized data.
!
        DO I=N_CLOUD_TOP, N_LAYER

!
!         To avoid the repetition of blocks of code or excessive
!         use of memory it is easiest to have an outer loop over
!         layers.
!
!
          SELECT CASE(I_PARAMETRIZATION_DROP)

          CASE(IP_SLINGO_SCHRECKER, IP_SLINGO_SCHR_PHF)
!
            DO LL=1, N_CLOUD_PROFILE(I)
              L=I_CLOUD_PROFILE(LL, I)
              K_EXT_TOT_CLOUD(L, I)                                     &
     &          =LIQ_WATER_MASS_FRAC(L, I)*(CLOUD_PARAMETER(1)          &
     &          +CLOUD_PARAMETER(2)/RADIUS_EFFECT(L, I))
              K_EXT_SCAT_CLOUD(L, I)=K_EXT_TOT_CLOUD(L, I)              &
     &          *(1.0E+00_Real64-CLOUD_PARAMETER(3)                     &
     &          -CLOUD_PARAMETER(4)*RADIUS_EFFECT(L, I))
              ASYMMETRY_PROCESS(L)=                                     &
     &          CLOUD_PARAMETER(5)+CLOUD_PARAMETER(6)                   &
     &          *RADIUS_EFFECT(L, I)
              PHASE_FNC_CLOUD(L, I, 1)=                                 &
     &          K_EXT_SCAT_CLOUD(L, I)*ASYMMETRY_PROCESS(L)
            ENDDO
!
!
          CASE(IP_ACKERMAN_STEPHENS)
!
!
            DO LL=1, N_CLOUD_PROFILE(I)
              L=I_CLOUD_PROFILE(LL, I)
              K_EXT_TOT_CLOUD(L, I)=LIQ_WATER_MASS_FRAC(L, I)           &
     &          *(CLOUD_PARAMETER(1)+CLOUD_PARAMETER(2)                 &
     &          *EXP(CLOUD_PARAMETER(3)*LOG(RADIUS_EFFECT(L, I))))
              K_EXT_SCAT_CLOUD(L, I)=K_EXT_TOT_CLOUD(L, I)              &
     &          *(1.0E+00_Real64-CLOUD_PARAMETER(4)                     &
     &          -CLOUD_PARAMETER(5)*EXP(CLOUD_PARAMETER(6)              &
     &          *LOG(RADIUS_EFFECT(L, I))))
              ASYMMETRY_PROCESS(L)                                      &
     &          =CLOUD_PARAMETER(7)+CLOUD_PARAMETER(8)                  &
     &          *EXP(CLOUD_PARAMETER(9)*LOG(RADIUS_EFFECT(L, I)))
              PHASE_FNC_CLOUD(L, I, 1)                                  &
     &          =K_EXT_SCAT_CLOUD(L, I)*ASYMMETRY_PROCESS(L)
            ENDDO
!
!
          CASE(IP_DROP_PADE_2)
!
!
            DO LL=1, N_CLOUD_PROFILE(I)
              L=I_CLOUD_PROFILE(LL, I)
              K_EXT_TOT_CLOUD(L, I)=LIQ_WATER_MASS_FRAC(L, I)           &
     &          *(CLOUD_PARAMETER(1)+RADIUS_EFFECT(L, I)                &
     &          *(CLOUD_PARAMETER(2)+RADIUS_EFFECT(L, I)                &
     &          *CLOUD_PARAMETER(3)))                                   &
     &          /(1.0E+00_Real64+RADIUS_EFFECT(L, I)                    &
     &          *(CLOUD_PARAMETER(4)+RADIUS_EFFECT(L, I)                &
     &          *(CLOUD_PARAMETER(5)+RADIUS_EFFECT(L, I)                &
     &          *CLOUD_PARAMETER(6))))
              K_EXT_SCAT_CLOUD(L, I)=K_EXT_TOT_CLOUD(L, I)              &
     &          *(1.0E+00_Real64                                        &
     &          -(CLOUD_PARAMETER(7)+RADIUS_EFFECT(L, I)                &
     &          *(CLOUD_PARAMETER(8)+RADIUS_EFFECT(L, I)                &
     &          *CLOUD_PARAMETER(9)))                                   &
     &          /(1.0E+00_Real64+RADIUS_EFFECT(L, I)                    &
     &          *(CLOUD_PARAMETER(10)+RADIUS_EFFECT(L, I)               &
     &          *CLOUD_PARAMETER(11))))
              ASYMMETRY_PROCESS(L)                                      &
     &          =(CLOUD_PARAMETER(12)+RADIUS_EFFECT(L, I)               &
     &          *(CLOUD_PARAMETER(13)+RADIUS_EFFECT(L, I)               &
     &          *CLOUD_PARAMETER(14)))                                  &
     &          /(1.0E+00_Real64+RADIUS_EFFECT(L, I)                    &
     &          *(CLOUD_PARAMETER(15)+RADIUS_EFFECT(L, I)               &
     &          *CLOUD_PARAMETER(16)))
              PHASE_FNC_CLOUD(L, I, 1)                                  &
     &          =K_EXT_SCAT_CLOUD(L, I)*ASYMMETRY_PROCESS(L)
            ENDDO
!
!
          END SELECT
!
!
!         Since these parametrizations include only the asymmetry,
!         it seems reasonable to extend them to higher
!         truncations using the Henyey-Greenstein phase function.
!
          DO LS=2, N_ORDER_PHASE
            DO LL=1, N_CLOUD_PROFILE(I)
              L=I_CLOUD_PROFILE(LL, I)
              PHASE_FNC_CLOUD(L, I, LS)                                 &
     &          =PHASE_FNC_CLOUD(L, I, LS-1)*ASYMMETRY_PROCESS(L)
            ENDDO
          ENDDO
!
          IF (L_RESCALE) THEN
            DO LL=1, N_CLOUD_PROFILE(I)
              L=I_CLOUD_PROFILE(LL, I)
                FORWARD_SCATTER_CLOUD(L, I)                             &
     &            =K_EXT_SCAT_CLOUD(L, I)                               &
     &            *ASYMMETRY_PROCESS(L)**N_ORDER_FORWARD
            ENDDO
          ENDIF
!
          IF (L_SOLAR_PHF) THEN
!           Calculate the solar phase function to higher accuracy.
            DO ID=1, N_DIRECTION
!             The Legendre polynomials are not stored so as to reduce
!             the requirement for memory at very high orders of solar
!             truncation.
              DO LL=1, N_CLOUD_PROFILE(I)
                L=I_CLOUD_PROFILE(LL, I)
!               Initialize the Legendre polynomials at the zeroth and
!               first orders.
                P_LEGENDRE_LS_M1(L)=1.0E+00_Real64
                P_LEGENDRE_LS(L)=COS_SOL_VIEW(L, ID)
                KS_PHF(L)=K_EXT_SCAT_CLOUD(L, I)*ASYMMETRY_PROCESS(L)
                PHASE_FNC_SOLAR_CLOUD(L, I, ID)=K_EXT_SCAT_CLOUD(L, I)  &
     &            +KS_PHF(L)*P_LEGENDRE_LS(L)*REAL(2*1+1, Real64)
              ENDDO

              DO LS=2, N_ORDER_PHASE_SOLAR
!               Calculate higher orders by recurrences.
                CNST1=1.0E+00_Real64-1.0E+00_Real64/REAL(LS, Real64)
                DO LL=1, N_CLOUD_PROFILE(I)
                  L=I_CLOUD_PROFILE(LL, I)
                  P_LEGENDRE_TMP(L)=P_LEGENDRE_LS(L)
                  P_LEGENDRE_LS(L)                                      &
     &              =(1.0E+00_Real64+CNST1)*P_LEGENDRE_LS(L)            &
     &              *COS_SOL_VIEW(L, ID)-CNST1*P_LEGENDRE_LS_M1(L)
                  P_LEGENDRE_LS_M1(L)=P_LEGENDRE_TMP(L)
                  KS_PHF(L)=KS_PHF(L)*ASYMMETRY_PROCESS(L)
                  PHASE_FNC_SOLAR_CLOUD(L, I, ID)                       &
     &              =PHASE_FNC_SOLAR_CLOUD(L, I, ID)                    &
     &              +KS_PHF(L)*P_LEGENDRE_LS(L)                         &
     &              *REAL(2*LS+1, Real64)

                ENDDO

              ENDDO
!
            ENDDO
!
!           Continue to an extra order to find the rescaling
!           for the solar beam.
            IF (L_RESCALE) THEN
              DO LL=1, N_CLOUD_PROFILE(I)
                L=I_CLOUD_PROFILE(LL, I)
                FORWARD_SOLAR_CLOUD(L, I)                               &
     &            =KS_PHF(L)*ASYMMETRY_PROCESS(L)
              ENDDO
            ENDIF
!

          ENDIF
!
        ENDDO
!

      ELSE IF (.NOT. L_HENYEY_GREENSTEIN_PF .AND.                       &
     &        (I_PARAMETRIZATION_DROP == IP_SLINGO_SCHR_PHF) ) THEN

        DO I=N_CLOUD_TOP, N_LAYER

! To avoid the repetition of blocks of code or excessive
! use of memory it is easiest to have an outer loop over
! layers

           DO LL=1, N_CLOUD_PROFILE(I)
              L=I_CLOUD_PROFILE(LL, I)
              K_EXT_TOT_CLOUD(L, I)                                     &
     &          =LIQ_WATER_MASS_FRAC(L, I)*(CLOUD_PARAMETER(1)          &
     &          +CLOUD_PARAMETER(2)/RADIUS_EFFECT(L, I))
              K_EXT_SCAT_CLOUD(L, I)=K_EXT_TOT_CLOUD(L, I)              &
     &          *(1.0E+00_Real64-CLOUD_PARAMETER(3)                     &
     &          -CLOUD_PARAMETER(4)*RADIUS_EFFECT(L, I))

           ENDDO

           DO LS=1, N_ORDER_PHASE
            DO LL=1, N_CLOUD_PROFILE(I)
             L=I_CLOUD_PROFILE(LL, I)
             PHASE_FNC_CLOUD(L, I, LS)                                  &
     &        =K_EXT_SCAT_CLOUD(L,I)*(CLOUD_PARAMETER(2*LS+3)           &
     &        +CLOUD_PARAMETER(2*LS+4)*RADIUS_EFFECT(L,I))

              ENDDO
           ENDDO

           LS=N_ORDER_FORWARD

           DO LL=1, N_CLOUD_PROFILE(I)
            L=I_CLOUD_PROFILE(LL, I)
            FORWARD_SCATTER_CLOUD(L,I)                                  &
     &        =K_EXT_SCAT_CLOUD(L,I)*(CLOUD_PARAMETER(2*LS+3)           &
     &        +CLOUD_PARAMETER(2*LS+4)*RADIUS_EFFECT(L,I))

           ENDDO

           IF (L_SOLAR_PHF) THEN


! Calculate the solar phase function to higher accuracy.
             DO ID=1, N_DIRECTION
! The Legendre polynomials are not stored so as to reduce
! the requirement for memory at very high orders of solar
! truncation.
              DO LL=1, N_CLOUD_PROFILE(I)
                L=I_CLOUD_PROFILE(LL, I)
! Initialize the Legendre polynomials at the zeroth and
! first orders.
                P_LEGENDRE_LS_M1(L)=1.0E+00_Real64
                P_LEGENDRE_LS(L)=COS_SOL_VIEW(L, ID)
                PHASE_FNC_SOLAR_CLOUD(L, I, ID)=K_EXT_SCAT_CLOUD(L, I)  &
     &             + PHASE_FNC_CLOUD(L, I, 1)                           &
     &             * P_LEGENDRE_LS(L)*REAL(2*1+1, Real64)

              ENDDO

              DO LS=2, N_ORDER_PHASE_SOLAR
! Calculate higher orders by recurrences.
                CNST1=1.0E+00_Real64-1.0E+00_Real64/REAL(LS, Real64)
                DO LL=1, N_CLOUD_PROFILE(I)
                  L=I_CLOUD_PROFILE(LL, I)
                  P_LEGENDRE_TMP(L)=P_LEGENDRE_LS(L)
                  P_LEGENDRE_LS(L)                                      &
     &              =(1.0E+00_Real64+CNST1)*P_LEGENDRE_LS(L)            &
     &              *COS_SOL_VIEW(L, ID)-CNST1*P_LEGENDRE_LS_M1(L)
                  P_LEGENDRE_LS_M1(L)=P_LEGENDRE_TMP(L)

                  PHF_TMP=CLOUD_PARAMETER(2*LS+3)                       &
     &                   + RADIUS_EFFECT(L,I)                           &
     &                   * CLOUD_PARAMETER(2*LS+4)


                 KS_PHF(L)=K_EXT_SCAT_CLOUD(L,I)*PHF_TMP
                 PHASE_FNC_SOLAR_CLOUD(L, I, ID)                        &
     &              = PHASE_FNC_SOLAR_CLOUD(L, I, ID)                   &
     &              + KS_PHF(L)*P_LEGENDRE_LS(L)                        &
     &              * REAL(2*LS+1, Real64)

                ENDDO
              ENDDO
             ENDDO
!
! Continue to an extra order to find the rescaling
! for the solar beam.
!
             IF (L_RESCALE) THEN
              LS=N_ORDER_PHASE_SOLAR+1
              DO LL=1, N_CLOUD_PROFILE(I)
                L=I_CLOUD_PROFILE(LL, I)

                PHF_TMP=CLOUD_PARAMETER(2*LS+3)                         &
     &                +RADIUS_EFFECT(L, I)                              &
     &                *CLOUD_PARAMETER(2*LS+4)
!
                FORWARD_SOLAR_CLOUD(L, I)                               &
     &            =K_EXT_SCAT_CLOUD(L, I)*PHF_TMP

              ENDDO
             ENDIF
           ENDIF


        ENDDO

      ENDIF
!
#if defined(OBSERVED)
!
      IF (I_PARAMETRIZATION_DROP == IP_DROP_UNPARAMETRIZED) THEN
!
        CALL PRSC_OPT_PROP(IERR                                         &
     &    , N_PROFILE, N_CLOUD_TOP, N_LAYER                             &
     &    , L_RESCALE, N_ORDER_FORWARD                                  &
     &    , L_HENYEY_GREENSTEIN_PF, N_ORDER_PHASE                       &
     &    , P, DENSITY                                                  &
     &    , N_OPT_LEVEL_CLOUD_PRSC                                      &
     &    , DROP_PRESSURE_PRSC, DROP_ABSORPTION_PRSC                    &
     &    , DROP_SCATTERING_PRSC, DROP_PHASE_FNC_PRSC                   &
     &    , K_EXT_TOT_CLOUD, K_EXT_SCAT_CLOUD, PHASE_FNC_CLOUD          &
     &    , FORWARD_SCATTER_CLOUD, FORWARD_SOLAR_CLOUD                  &
     &    , L_SOLAR_PHF, N_ORDER_PHASE_SOLAR, N_DIRECTION, COS_SOL_VIEW &
     &    , PHASE_FNC_SOLAR_CLOUD                                       &
     &    , ND_PROFILE, ND_RADIANCE_PROFILE, ND_LAYER, ID_CT, ND_LAYER  &
     &    , ND_DIRECTION                                                &
     &    , ND_PROFILE_PRSC, ND_OPT_LEVEL_PRSC                          &
     &    , ND_PHASE_TERM, ND_MAX_ORDER                                 &
     &    )
!
      ENDIF
#endif
!
      IF ( (I_PARAMETRIZATION_DROP /= IP_SLINGO_SCHRECKER).AND.         &
     &     (I_PARAMETRIZATION_DROP /= IP_ACKERMAN_STEPHENS).AND.        &
#if defined(OBSERVED)
     &     (I_PARAMETRIZATION_DROP /= IP_DROP_UNPARAMETRIZED).AND.      &
#endif
     &     (I_PARAMETRIZATION_DROP /= IP_DROP_PADE_2).AND.              &
     &     (I_PARAMETRIZATION_DROP /= IP_SLINGO_SCHR_PHF) ) THEN
        WRITE(IU_ERR, '(/A)') '*** Error: An invalid parametrization '  &
     &    //'of cloud droplets has been selected.'
        IERR=I_ERR_FATAL
        RETURN
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE OPT_PROP_WATER_CLOUD
#endif
