#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate optical properties of ice clouds.
!
! Method:
!       If the optical properties come from an observational
!       distribution a separate subroutine is called. Otherwise
!       appropriate mean quantities in the layer are calculated
!       as the parametrization requires and these values are
!       substituted into the parametrization to give the optical
!       properties.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE OPT_PROP_ICE_CLOUD(IERR                                &
     &  , N_PROFILE, N_LAYER, N_CLOUD_TOP                               &
     &  , N_CLOUD_PROFILE, I_CLOUD_PROFILE                              &
     &  , N_ORDER_PHASE, L_RESCALE, N_ORDER_FORWARD                     &
     &  , L_HENYEY_GREENSTEIN_PF, L_SOLAR_PHF, N_ORDER_PHASE_SOLAR      &
     &  , N_DIRECTION, COS_SOL_VIEW                                     &
     &  , I_PARAMETRIZATION_ICE, ICE_CLOUD_PARAMETER                    &
     &  , ICE_MASS_FRAC, DIM_CHAR_ICE                                   &
#if defined(OBSERVED)
     &  , P, T, DENSITY                                                 &
     &  , N_OPT_LEVEL_CLOUD_PRSC, ICE_PRESSURE_PRSC                     &
     &  , ICE_ABSORPTION_PRSC, ICE_SCATTERING_PRSC                      &
     &  , ICE_PHASE_FNC_PRSC                                            &
#endif
#if defined(UM)
     &  , T, DENSITY                                                    &
#endif
     &  , K_EXT_TOT_CLOUD, K_EXT_SCAT_CLOUD                             &
     &  , PHASE_FNC_CLOUD, FORWARD_SCATTER_CLOUD                        &
     &  , FORWARD_SOLAR_CLOUD, PHASE_FNC_SOLAR_CLOUD                    &
     &  , ND_PROFILE, ND_RADIANCE_PROFILE, ND_LAYER, ID_CT              &
     &  , ND_DIRECTION, ND_PHASE_TERM, ND_MAX_ORDER                     &
     &  , ND_CLOUD_PARAMETER                                            &
#if defined(OBSERVED)
     &  , ND_PROFILE_PRSC, ND_OPT_LEVEL_PRSC                            &
#endif
     &  )
!
!
      IMPLICIT NONE
!
!
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for profiles
     &  , ND_RADIANCE_PROFILE                                           &
!           Size allocated for points where radiances are calculated
     &  , ND_LAYER                                                      &
!           Size allocated for layers
     &  , ID_CT                                                         &
!           Topmost declared cloudy layer
     &  , ND_DIRECTION                                                  &
!           Size allocated for viewing directions
     &  , ND_CLOUD_PARAMETER                                            &
!           Size allocated for cloud parameters
     &  , ND_PHASE_TERM                                                 &
!           Size allocated for terms in the phase function
     &  , ND_MAX_ORDER
!           Size allocated for orders of spherical harmonics
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
#include "ice_cloud_param_pcf3z.h"
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
     &  , N_ORDER_PHASE                                                 &
!           Order of the phase function
     &  , N_ORDER_PHASE_SOLAR                                           &
!           Number of terms to retain in single scattered solar
!           phase function
     &  , N_ORDER_FORWARD                                               &
!           Order used in forming the forward scattering parameter
     &  , N_CLOUD_TOP                                                   &
!           Topmost cloudy layer
     &  , I_PARAMETRIZATION_ICE                                         &
!           Treatment of ice crystals
#if defined(OBSERVED)
     &  , N_OPT_LEVEL_CLOUD_PRSC                                        &
!           Number of levels of prescribed data
#endif
     &  , N_CLOUD_PROFILE(ID_CT: ND_LAYER)                              &
!           Number of cloudy profiles
     &  , I_CLOUD_PROFILE(ND_PROFILE, ID_CT: ND_LAYER)
!           Profiles containing clouds
      LOGICAL, INTENT(IN) ::                                            &
     &    L_RESCALE                                                     &
!           Delta-rescaling required
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
     &    ICE_CLOUD_PARAMETER(ND_CLOUD_PARAMETER)                       &
!           Ice cloud parameters
     &  , ICE_MASS_FRAC(ND_PROFILE, ID_CT: ND_LAYER)                    &
!           Ice mass fraction
     &  , DIM_CHAR_ICE(ND_PROFILE, ID_CT: ND_LAYER)                     &
!           Characteristic dimension for crystals
     &  , T(ND_PROFILE, ND_LAYER)                                       &
!           Temperature
     &  , DENSITY(ND_PROFILE, ND_LAYER)
!           Atmospheric density
#if defined(OBSERVED)
      REAL  (Real64), INTENT(IN) ::                                     &
     &    P(ND_PROFILE, ND_LAYER)                                       &
!           Pressure
     &  , ICE_PRESSURE_PRSC(ND_PROFILE_PRSC, ND_OPT_LEVEL_PRSC)         &
!           Pressure at which optical properties are prescribed
     &  , ICE_ABSORPTION_PRSC(ND_PROFILE_PRSC, ND_OPT_LEVEL_PRSC)       &
!           Prescribed absorption by ice crystals
     &  , ICE_SCATTERING_PRSC(ND_PROFILE_PRSC, ND_OPT_LEVEL_PRSC)       &
!           Prescribed absorption by ice crystals
     &  , ICE_PHASE_FNC_PRSC(ND_PROFILE_PRSC, ND_OPT_LEVEL_PRSC         &
     &      , ND_PHASE_TERM)
!           Prescribed phase functions of ice crystals
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
!           Asymmetry factor for current process
     &  , OMEGA                                                         &
!           Albedo of single scattering for the current process
     &  , X                                                             &
!           Temporary algebraic variable
     &  , Y                                                             &
!           Temporary algebraic variable
     &  , T_CELSIUS                                                     &
!           Temperature in celsius
     &  , TEMP_CORRECTION                                               &
!           Temperature correction
     &  , PHF_TMP
!           Temporary phase function
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
!
      IF ( (I_PARAMETRIZATION_ICE == IP_SLINGO_SCHRECKER_ICE).OR.       &
     &     (I_PARAMETRIZATION_ICE == IP_ICE_ADT).OR.                    &
     &     (I_PARAMETRIZATION_ICE == IP_ICE_ADT_10).OR.                 &
     &     (I_PARAMETRIZATION_ICE == IP_ICE_FU_SOLAR).OR.               &
     &     (I_PARAMETRIZATION_ICE == IP_ICE_FU_IR).OR.                  &
     &     (I_PARAMETRIZATION_ICE == IP_SUN_SHINE_VN2_IR).OR.           &
     &     (I_PARAMETRIZATION_ICE == IP_SUN_SHINE_VN2_IR).OR.           &
     &     ( L_HENYEY_GREENSTEIN_PF .AND.                               &
     &       (I_PARAMETRIZATION_ICE == IP_SLINGO_SCHR_ICE_PHF) ).OR.    &
     &     ( L_HENYEY_GREENSTEIN_PF .AND.                               &
     &       (I_PARAMETRIZATION_ICE == IP_ICE_FU_PHF) ) ) THEN
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
          SELECT CASE(I_PARAMETRIZATION_ICE)
!
          CASE(IP_SLINGO_SCHRECKER_ICE, IP_SLINGO_SCHR_ICE_PHF)
!
!
            DO LL=1, N_CLOUD_PROFILE(I)
              L=I_CLOUD_PROFILE(LL, I)
              K_EXT_TOT_CLOUD(L, I)                                     &
     &          =ICE_MASS_FRAC(L, I)*(ICE_CLOUD_PARAMETER(1)            &
     &          +ICE_CLOUD_PARAMETER(2)/DIM_CHAR_ICE(L, I))
              K_EXT_SCAT_CLOUD(L, I)=K_EXT_TOT_CLOUD(L, I)              &
     &          *(1.0_Real64-ICE_CLOUD_PARAMETER(3)                     &
     &        -ICE_CLOUD_PARAMETER(4)*DIM_CHAR_ICE(L, I))
              ASYMMETRY_PROCESS(L)                                      &
     &          =ICE_CLOUD_PARAMETER(5)+ICE_CLOUD_PARAMETER(6)          &
     &          *DIM_CHAR_ICE(L, I)
              PHASE_FNC_CLOUD(L, I, 1)                                  &
     &          =K_EXT_SCAT_CLOUD(L, I)*ASYMMETRY_PROCESS(L)
            ENDDO
!
!
          CASE (IP_ICE_ADT)
!
!
            DO LL=1, N_CLOUD_PROFILE(I)
              L=I_CLOUD_PROFILE(LL, I)
              X=LOG(DIM_CHAR_ICE(L, I)/ICE_CLOUD_PARAMETER(10))
              IF (X >  0.0_Real64) THEN
!               Large mode.
                K_EXT_TOT_CLOUD(L, I)=ICE_MASS_FRAC(L, I)               &
     &            *EXP(ICE_CLOUD_PARAMETER(1)                           &
     &            +X*(ICE_CLOUD_PARAMETER(2)                            &
     &            +X*(ICE_CLOUD_PARAMETER(3)                            &
     &            +X*(ICE_CLOUD_PARAMETER(4)                            &
     &            +X*ICE_CLOUD_PARAMETER(5)))))
              ELSE IF (X <= 0.0_Real64) THEN
!               Small mode.
                K_EXT_TOT_CLOUD(L, I)=ICE_MASS_FRAC(L, I)               &
     &            *EXP(ICE_CLOUD_PARAMETER(1)                           &
     &            +X*(ICE_CLOUD_PARAMETER(6)                            &
     &            +X*(ICE_CLOUD_PARAMETER(7)                            &
     &            +X*(ICE_CLOUD_PARAMETER(8)                            &
     &            +X*ICE_CLOUD_PARAMETER(9)))))
              ENDIF
              X=LOG(DIM_CHAR_ICE(L, I)/ICE_CLOUD_PARAMETER(20))
              IF (X >  0.0_Real64) THEN
!               Large mode.
                K_EXT_SCAT_CLOUD(L, I)=K_EXT_TOT_CLOUD(L, I)            &
     &            *(1.0_Real64-(ICE_CLOUD_PARAMETER(11)                 &
     &            +X*(ICE_CLOUD_PARAMETER(12)                           &
     &            +X*(ICE_CLOUD_PARAMETER(13)                           &
     &            +X*(ICE_CLOUD_PARAMETER(14)                           &
     &            +X*ICE_CLOUD_PARAMETER(15))))))
              ELSE IF (X <= 0.0_Real64) THEN
!               Small mode.
                K_EXT_SCAT_CLOUD(L, I)=K_EXT_TOT_CLOUD(L, I)            &
     &            *(1.0_Real64-(ICE_CLOUD_PARAMETER(11)                 &
     &            +X*(ICE_CLOUD_PARAMETER(16)                           &
     &            +X*(ICE_CLOUD_PARAMETER(17)                           &
     &            +X*(ICE_CLOUD_PARAMETER(18)                           &
     &            +X*ICE_CLOUD_PARAMETER(19))))))
              ENDIF
              X=LOG(DIM_CHAR_ICE(L, I)/ICE_CLOUD_PARAMETER(30))
              IF (X >  0.0_Real64) THEN
!               Large mode.
                ASYMMETRY_PROCESS(L)=ICE_CLOUD_PARAMETER(21)            &
     &            +X*(ICE_CLOUD_PARAMETER(22)                           &
     &            +X*(ICE_CLOUD_PARAMETER(23)                           &
     &            +X*(ICE_CLOUD_PARAMETER(24)                           &
     &            +X*ICE_CLOUD_PARAMETER(25))))
              ELSE IF (X <= 0.0_Real64) THEN
!               Small mode.
                ASYMMETRY_PROCESS(L)=ICE_CLOUD_PARAMETER(21)            &
     &            +X*(ICE_CLOUD_PARAMETER(26)                           &
     &            +X*(ICE_CLOUD_PARAMETER(27)                           &
     &            +X*(ICE_CLOUD_PARAMETER(28)                           &
     &            +X*ICE_CLOUD_PARAMETER(29))))
              ENDIF
              PHASE_FNC_CLOUD(L, I, 1)                                  &
     &          =K_EXT_SCAT_CLOUD(L, I)*ASYMMETRY_PROCESS(L)
            ENDDO
!
!
          CASE(IP_ICE_ADT_10)
!
!
            DO LL=1, N_CLOUD_PROFILE(I)
              L=I_CLOUD_PROFILE(LL, I)
              X=DIM_CHAR_ICE(L, I)/ICE_CLOUD_PARAMETER(12)
              Y=ICE_CLOUD_PARAMETER(6)                                  &
     &          +X*(ICE_CLOUD_PARAMETER(7)                              &
     &          +X*(ICE_CLOUD_PARAMETER(8)                              &
     &          +X*(ICE_CLOUD_PARAMETER(9)                              &
     &          +X*(ICE_CLOUD_PARAMETER(10)                             &
     &          +X*ICE_CLOUD_PARAMETER(11)))))
              K_EXT_TOT_CLOUD(L, I)=ICE_MASS_FRAC(L, I)                 &
     &          *EXP(ICE_CLOUD_PARAMETER(1)                             &
     &          +X*(ICE_CLOUD_PARAMETER(2)                              &
     &          +X*(ICE_CLOUD_PARAMETER(3)                              &
     &          +X*(ICE_CLOUD_PARAMETER(4)                              &
     &          +X*(ICE_CLOUD_PARAMETER(5)                              &
     &          +X*Y)))))
              X=DIM_CHAR_ICE(L, I)/ICE_CLOUD_PARAMETER(24)
              Y=ICE_CLOUD_PARAMETER(18)                                 &
     &          +X*(ICE_CLOUD_PARAMETER(19)                             &
     &          +X*(ICE_CLOUD_PARAMETER(20)                             &
     &          +X*(ICE_CLOUD_PARAMETER(21)                             &
     &          +X*(ICE_CLOUD_PARAMETER(22)                             &
     &          +X*ICE_CLOUD_PARAMETER(23)))))
              K_EXT_SCAT_CLOUD(L, I)=K_EXT_TOT_CLOUD(L, I)              &
     &          *(1.0_Real64-(ICE_CLOUD_PARAMETER(13)                   &
     &          +X*(ICE_CLOUD_PARAMETER(14)                             &
     &          +X*(ICE_CLOUD_PARAMETER(15)                             &
     &          +X*(ICE_CLOUD_PARAMETER(16)                             &
     &          +X*(ICE_CLOUD_PARAMETER(17)                             &
     &          +X*Y))))))
              X=DIM_CHAR_ICE(L, I)/ICE_CLOUD_PARAMETER(36)
              Y=ICE_CLOUD_PARAMETER(30)                                 &
     &          +X*(ICE_CLOUD_PARAMETER(31)                             &
     &          +X*(ICE_CLOUD_PARAMETER(32)                             &
     &          +X*(ICE_CLOUD_PARAMETER(33)                             &
     &          +X*(ICE_CLOUD_PARAMETER(34)                             &
     &          +X*ICE_CLOUD_PARAMETER(35)))))
              ASYMMETRY_PROCESS(L)=ICE_CLOUD_PARAMETER(25)              &
     &          +X*(ICE_CLOUD_PARAMETER(26)                             &
     &          +X*(ICE_CLOUD_PARAMETER(27)                             &
     &          +X*(ICE_CLOUD_PARAMETER(28)                             &
     &          +X*(ICE_CLOUD_PARAMETER(29)                             &
     &          +X*Y))))
              PHASE_FNC_CLOUD(L, I, 1)                                  &
     &          =K_EXT_SCAT_CLOUD(L, I)*ASYMMETRY_PROCESS(L)
            ENDDO
!
!
          CASE (IP_SUN_SHINE_VN2_VIS)
!
!
            DO LL=1, N_CLOUD_PROFILE(I)
              L=I_CLOUD_PROFILE(LL, I)
              T_CELSIUS=T(L, I)-2.7316E+02_Real64
              TEMP_CORRECTION=1.047_Real64                              &
     &          +T_CELSIUS*(-9.13E-05_Real64+T_CELSIUS                  &
     &          *(2.026E-04_Real64-1.056E-05_Real64*T_CELSIUS))
              K_EXT_TOT_CLOUD(L, I)=TEMP_CORRECTION*ICE_MASS_FRAC(L, I) &
     &          /(3.05548E-02_Real64                                    &
     &          +2.54802E+02_Real64*DENSITY(L, I)*ICE_MASS_FRAC(L, I))
              K_EXT_SCAT_CLOUD(L, I)=K_EXT_TOT_CLOUD(L, I)              &
     &          *(1.0_Real64-ICE_CLOUD_PARAMETER(1)                     &
     &          *EXP(ICE_CLOUD_PARAMETER(2)                             &
     &          *LOG(DENSITY(L, I)*ICE_MASS_FRAC(L, I)+1.0E-12_Real64)))&
     &          *(1.0_Real64+ICE_CLOUD_PARAMETER(5)                     &
     &          *(TEMP_CORRECTION-1.0_Real64)/TEMP_CORRECTION)
              ASYMMETRY_PROCESS(L)                                      &
     &          =ICE_CLOUD_PARAMETER(3)*EXP(ICE_CLOUD_PARAMETER(4)      &
     &          *LOG(DENSITY(L, I)*ICE_MASS_FRAC(L, I)+1.0E-12_Real64)) &
     &          *(1.0_Real64+ICE_CLOUD_PARAMETER(6)                     &
     &          *(TEMP_CORRECTION-1.0_Real64)/TEMP_CORRECTION)
              PHASE_FNC_CLOUD(L, I, 1)                                  &
     &          =K_EXT_SCAT_CLOUD(L, I)*ASYMMETRY_PROCESS(L)
            ENDDO
!
!
          CASE(IP_SUN_SHINE_VN2_IR)
!
!
            DO LL=1, N_CLOUD_PROFILE(I)
              L=I_CLOUD_PROFILE(LL, I)
              T_CELSIUS=T(L, I)-2.7316E+02_Real64
              TEMP_CORRECTION=1.047_Real64+T_CELSIUS                    &
     &          *(-9.13E-05_Real64+T_CELSIUS                            &
     &          *(2.026E-04_Real64-1.056E-05_Real64*T_CELSIUS))
              K_EXT_TOT_CLOUD(L, I)=TEMP_CORRECTION*ICE_MASS_FRAC(L, I) &
     &          /(6.30689E-02_Real64                                    &
     &          +2.65874E+02_Real64*DENSITY(L, I)*ICE_MASS_FRAC(L, I))
            ENDDO
!
!
          CASE(IP_ICE_FU_SOLAR)
!
!
            DO LL=1, N_CLOUD_PROFILE(I)
              L=I_CLOUD_PROFILE(LL, I)
              K_EXT_TOT_CLOUD(L, I)                                     &
     &          =ICE_MASS_FRAC(L, I)*(ICE_CLOUD_PARAMETER(1)            &
     &          +ICE_CLOUD_PARAMETER(2)/DIM_CHAR_ICE(L, I))
              OMEGA=1.0_Real64-(ICE_CLOUD_PARAMETER(3)                  &
     &          +DIM_CHAR_ICE(L, I)*(ICE_CLOUD_PARAMETER(4)             &
     &          +DIM_CHAR_ICE(L, I)*(ICE_CLOUD_PARAMETER(5)             &
     &          +DIM_CHAR_ICE(L, I)*ICE_CLOUD_PARAMETER(6))))
              K_EXT_SCAT_CLOUD(L, I)=K_EXT_TOT_CLOUD(L, I)*OMEGA
              ASYMMETRY_PROCESS(L)                                      &
     &          =ICE_CLOUD_PARAMETER(7)+DIM_CHAR_ICE(L, I)              &
     &          *(ICE_CLOUD_PARAMETER(8)+DIM_CHAR_ICE(L, I)             &
     &          *(ICE_CLOUD_PARAMETER(9)+DIM_CHAR_ICE(L, I)             &
     &          *(ICE_CLOUD_PARAMETER(10))))
              PHASE_FNC_CLOUD(L, I, 1)                                  &
     &          =K_EXT_SCAT_CLOUD(L, I)*ASYMMETRY_PROCESS(L)
!             The forward scattering will be limited later.
              FORWARD_SCATTER_CLOUD(L, I)=K_EXT_SCAT_CLOUD(L, I)        &
     &          *(1.0_Real64                                            &
     &          /MAX(1.0_Real64, 2.0_Real64*OMEGA)                      &
     &          +ICE_CLOUD_PARAMETER(11)+DIM_CHAR_ICE(L, I)             &
     &          *(ICE_CLOUD_PARAMETER(12)+DIM_CHAR_ICE(L, I)            &
     &          *(ICE_CLOUD_PARAMETER(13)+DIM_CHAR_ICE(L, I)            &
     &          *(ICE_CLOUD_PARAMETER(14)))))
            ENDDO
!
!
          CASE(IP_ICE_FU_IR)
!
!
            DO LL=1, N_CLOUD_PROFILE(I)
              L=I_CLOUD_PROFILE(LL, I)
              K_EXT_TOT_CLOUD(L, I)=ICE_MASS_FRAC(L, I)                 &
     &          *((ICE_CLOUD_PARAMETER(3)/DIM_CHAR_ICE(L, I)            &
     &          +ICE_CLOUD_PARAMETER(2))/DIM_CHAR_ICE(L, I)             &
     &          +ICE_CLOUD_PARAMETER(1))
              K_EXT_SCAT_CLOUD(L, I)=K_EXT_TOT_CLOUD(L, I)              &
     &          -(ICE_MASS_FRAC(L, I)/DIM_CHAR_ICE(L, I))               &
     &          *(ICE_CLOUD_PARAMETER(4)+DIM_CHAR_ICE(L, I)             &
     &          *(ICE_CLOUD_PARAMETER(5)+DIM_CHAR_ICE(L, I)             &
     &          *(ICE_CLOUD_PARAMETER(6)+DIM_CHAR_ICE(L, I)             &
     &          *ICE_CLOUD_PARAMETER(7))))
              ASYMMETRY_PROCESS(L)                                      &
     &          =ICE_CLOUD_PARAMETER(8)+DIM_CHAR_ICE(L, I)              &
     &          *(ICE_CLOUD_PARAMETER(9)+DIM_CHAR_ICE(L, I)             &
     &          *(ICE_CLOUD_PARAMETER(10)+DIM_CHAR_ICE(L, I)            &
     &          *(ICE_CLOUD_PARAMETER(11))))
              PHASE_FNC_CLOUD(L, I, 1)                                  &
     &          =K_EXT_SCAT_CLOUD(L, I)*ASYMMETRY_PROCESS(L)
            ENDDO
!
!
          CASE(IP_ICE_FU_PHF)
!
!
            DO LL=1, N_CLOUD_PROFILE(I)
              L=I_CLOUD_PROFILE(LL, I)
              X=DIM_CHAR_ICE(L, I)/ICE_CLOUD_PARAMETER(4)
              K_EXT_TOT_CLOUD(L, I)=ICE_MASS_FRAC(L, I)                 &
     &          *((ICE_CLOUD_PARAMETER(3)/X                             &
     &          +ICE_CLOUD_PARAMETER(2))/X                              &
     &          +ICE_CLOUD_PARAMETER(1))
              X=DIM_CHAR_ICE(L, I)/ICE_CLOUD_PARAMETER(9)
              K_EXT_SCAT_CLOUD(L, I)=K_EXT_TOT_CLOUD(L, I)              &
     &          *(1.0_Real64                                            &
     &          -(ICE_CLOUD_PARAMETER(5)+X                              &
     &          *(ICE_CLOUD_PARAMETER(6)+X                              &
     &          *(ICE_CLOUD_PARAMETER(7)+X                              &
     &          *ICE_CLOUD_PARAMETER(8)))))
              X=DIM_CHAR_ICE(L, I)/ICE_CLOUD_PARAMETER(14)
              ASYMMETRY_PROCESS(L)                                      &
     &          =ICE_CLOUD_PARAMETER(10)+DIM_CHAR_ICE(L, I)             &
     &          *(ICE_CLOUD_PARAMETER(11)+DIM_CHAR_ICE(L, I)            &
     &          *(ICE_CLOUD_PARAMETER(12)+DIM_CHAR_ICE(L, I)            &
     &          *(ICE_CLOUD_PARAMETER(13))))
              PHASE_FNC_CLOUD(L, I, 1)                                  &
     &          =K_EXT_SCAT_CLOUD(L, I)*ASYMMETRY_PROCESS(L)
            ENDDO
!
!
          END SELECT
!
!
!         Parametrizations which do not include explicit
!         representation of the higher moments are extended using the
!         Henyey-Greenstein phase function.
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
!
!           For most parameterizations the forward scattering
!           is determined from the asymmetry, but in the case of
!           Fu's parametrization it is defined specially, but must
!           be limited to avoid negative moments in the phase function.
            IF (I_PARAMETRIZATION_ICE == IP_ICE_FU_SOLAR) THEN
              DO LL=1, N_CLOUD_PROFILE(I)
                L=I_CLOUD_PROFILE(LL, I)
                FORWARD_SCATTER_CLOUD(L, I)                             &
     &            =MIN(FORWARD_SCATTER_CLOUD(L, I)                      &
     &            , K_EXT_SCAT_CLOUD(L, I)                              &
     &            *ASYMMETRY_PROCESS(L)**(N_ORDER_FORWARD-1))
              ENDDO
            ELSE
              DO LL=1, N_CLOUD_PROFILE(I)
                L=I_CLOUD_PROFILE(LL, I)
                FORWARD_SCATTER_CLOUD(L, I)                             &
     &            =K_EXT_SCAT_CLOUD(L, I)                               &
     &            *ASYMMETRY_PROCESS(L)**N_ORDER_FORWARD
              ENDDO
            ENDIF
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
                P_LEGENDRE_LS_M1(L)=1.0_Real64
                P_LEGENDRE_LS(L)=COS_SOL_VIEW(L, ID)
                KS_PHF(L)=K_EXT_SCAT_CLOUD(L, I)*ASYMMETRY_PROCESS(L)
                PHASE_FNC_SOLAR_CLOUD(L, I, ID)=K_EXT_SCAT_CLOUD(L, I)  &
     &            +KS_PHF(L)*P_LEGENDRE_LS(L)*REAL(2*1+1, Real64)
              ENDDO

              DO LS=2, N_ORDER_PHASE_SOLAR
!               Calculate higher orders by recurrences.
                CNST1=1.0_Real64-1.0_Real64/REAL(LS, Real64)
                DO LL=1, N_CLOUD_PROFILE(I)
                  L=I_CLOUD_PROFILE(LL, I)
                  P_LEGENDRE_TMP(L)=P_LEGENDRE_LS(L)
                  P_LEGENDRE_LS(L)                                      &
     &              =(1.0_Real64+CNST1)*P_LEGENDRE_LS(L)                &
     &              *COS_SOL_VIEW(L, ID)-CNST1*P_LEGENDRE_LS_M1(L)
                  P_LEGENDRE_LS_M1(L)=P_LEGENDRE_TMP(L)
                  KS_PHF(L)=KS_PHF(L)*ASYMMETRY_PROCESS(L)
                  PHASE_FNC_SOLAR_CLOUD(L, I, ID)                       &
     &              =PHASE_FNC_SOLAR_CLOUD(L, I, ID)                    &
     &              +KS_PHF(L)*P_LEGENDRE_LS(L)                         &
     &              *REAL(2*LS+1, Real64)
                ENDDO

              ENDDO
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
      ELSE IF ( .NOT.L_HENYEY_GREENSTEIN_PF .AND.                       &
     &     ( (I_PARAMETRIZATION_ICE == IP_SLINGO_SCHR_ICE_PHF).OR.      &
     &       (I_PARAMETRIZATION_ICE == IP_ICE_FU_PHF) ) ) THEN
!
        DO I=N_CLOUD_TOP, N_LAYER

!
!         To avoid the repetition of blocks of code or excessive
!         use of memory it is easiest to have an outer loop over
!         layers.
!
!
          IF (I_PARAMETRIZATION_ICE == IP_SLINGO_SCHR_ICE_PHF) THEN
!
            DO LL=1, N_CLOUD_PROFILE(I)
              L=I_CLOUD_PROFILE(LL, I)
              K_EXT_TOT_CLOUD(L, I)                                     &
     &          =ICE_MASS_FRAC(L, I)*(ICE_CLOUD_PARAMETER(1)            &
     &          +ICE_CLOUD_PARAMETER(2)/DIM_CHAR_ICE(L, I))
              K_EXT_SCAT_CLOUD(L, I)=K_EXT_TOT_CLOUD(L, I)              &
     &          *(1.0_Real64-ICE_CLOUD_PARAMETER(3)                     &
     &        -ICE_CLOUD_PARAMETER(4)*DIM_CHAR_ICE(L, I))
            ENDDO
            DO LS=1, N_ORDER_PHASE
              DO LL=1, N_CLOUD_PROFILE(I)
                L=I_CLOUD_PROFILE(LL, I)
                PHASE_FNC_CLOUD(L, I, LS)                               &
     &            =K_EXT_SCAT_CLOUD(L, I)*(ICE_CLOUD_PARAMETER(2*LS+3)  &
     &            +ICE_CLOUD_PARAMETER(2*LS+4)*DIM_CHAR_ICE(L, I))
              ENDDO
            ENDDO
            LS=N_ORDER_FORWARD
            DO LL=1, N_CLOUD_PROFILE(I)
              L=I_CLOUD_PROFILE(LL, I)
              FORWARD_SCATTER_CLOUD(L, I)                               &
     &          =K_EXT_SCAT_CLOUD(L, I)*(ICE_CLOUD_PARAMETER(2*LS+3)    &
     &          +ICE_CLOUD_PARAMETER(2*LS+4)*DIM_CHAR_ICE(L, I))
            ENDDO
!
          ELSE IF (I_PARAMETRIZATION_ICE == IP_ICE_FU_PHF) THEN
!
            DO LL=1, N_CLOUD_PROFILE(I)
              L=I_CLOUD_PROFILE(LL, I)
              X=DIM_CHAR_ICE(L, I)/ICE_CLOUD_PARAMETER(4)
              K_EXT_TOT_CLOUD(L, I)=ICE_MASS_FRAC(L, I)                 &
     &          *((ICE_CLOUD_PARAMETER(3)/X                             &
     &          +ICE_CLOUD_PARAMETER(2))/X                              &
     &          +ICE_CLOUD_PARAMETER(1))
              X=DIM_CHAR_ICE(L, I)/ICE_CLOUD_PARAMETER(9)
              K_EXT_SCAT_CLOUD(L, I)=K_EXT_TOT_CLOUD(L, I)              &
     &          *(1.0_Real64                                            &
     &          -(ICE_CLOUD_PARAMETER(5)+X                              &
     &          *(ICE_CLOUD_PARAMETER(6)+X                              &
     &          *(ICE_CLOUD_PARAMETER(7)+X                              &
     &          *ICE_CLOUD_PARAMETER(8)))))
            ENDDO
            DO LS=1, N_ORDER_PHASE
              DO LL=1, N_CLOUD_PROFILE(I)
                L=I_CLOUD_PROFILE(LL, I)
                X=DIM_CHAR_ICE(L, I)/ICE_CLOUD_PARAMETER(5*LS+9)
                PHASE_FNC_CLOUD(L, I, LS)                               &
     &            =K_EXT_SCAT_CLOUD(L, I)*(ICE_CLOUD_PARAMETER(5*LS+5)  &
     &            +X*(ICE_CLOUD_PARAMETER(5*LS+6)                       &
     &            +X*(ICE_CLOUD_PARAMETER(5*LS+7)                       &
     &            +X*ICE_CLOUD_PARAMETER(5*LS+8))))
              ENDDO
            ENDDO
            LS=N_ORDER_FORWARD
            DO LL=1, N_CLOUD_PROFILE(I)
              L=I_CLOUD_PROFILE(LL, I)
              X=DIM_CHAR_ICE(L, I)/ICE_CLOUD_PARAMETER(5*LS+9)
              FORWARD_SCATTER_CLOUD(L, I)                               &
     &          =K_EXT_SCAT_CLOUD(L, I)*(ICE_CLOUD_PARAMETER(5*LS+5)    &
     &          +X*(ICE_CLOUD_PARAMETER(5*LS+6)                         &
     &          +X*(ICE_CLOUD_PARAMETER(5*LS+7)                         &
     &          +X*ICE_CLOUD_PARAMETER(5*LS+8))))
            ENDDO
!
          ENDIF
!
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
                P_LEGENDRE_LS_M1(L)=1.0_Real64
                P_LEGENDRE_LS(L)=COS_SOL_VIEW(L, ID)
                PHASE_FNC_SOLAR_CLOUD(L, I, ID)=K_EXT_SCAT_CLOUD(L, I)  &
     &            +PHASE_FNC_CLOUD(L, I, 1)                             &
     &            *P_LEGENDRE_LS(L)*REAL(2*1+1, Real64)
              ENDDO

              DO LS=2, N_ORDER_PHASE_SOLAR
!               Calculate higher orders by recurrences. Moments of
!               the phase function cannot be taken from above as
!               we will typically require a much higher order here.
                CNST1=1.0_Real64-1.0_Real64/REAL(LS, Real64)
                DO LL=1, N_CLOUD_PROFILE(I)
                  L=I_CLOUD_PROFILE(LL, I)
                  P_LEGENDRE_TMP(L)=P_LEGENDRE_LS(L)
                  P_LEGENDRE_LS(L)                                      &
     &              =(1.0_Real64+CNST1)*P_LEGENDRE_LS(L)                &
     &              *COS_SOL_VIEW(L, ID)-CNST1*P_LEGENDRE_LS_M1(L)
                  P_LEGENDRE_LS_M1(L)=P_LEGENDRE_TMP(L)
!
                  SELECT CASE(I_PARAMETRIZATION_ICE)
!
                    CASE(IP_SLINGO_SCHR_ICE_PHF)
                      PHF_TMP=ICE_CLOUD_PARAMETER(2*LS+3)               &
     &                  +DIM_CHAR_ICE(L, I)                             &
     &                  *ICE_CLOUD_PARAMETER(2*LS+4)
!
                    CASE(IP_ICE_FU_PHF)
                      X=DIM_CHAR_ICE(L, I)/ICE_CLOUD_PARAMETER(5*LS+9)
                      PHF_TMP                                           &
     &                  =(ICE_CLOUD_PARAMETER(5*LS+5)                   &
     &                  +X*(ICE_CLOUD_PARAMETER(5*LS+6)                 &
     &                  +X*(ICE_CLOUD_PARAMETER(5*LS+7)                 &
     &                  +X*ICE_CLOUD_PARAMETER(5*LS+8))))
!
                  END SELECT
!
                  KS_PHF(L)=K_EXT_SCAT_CLOUD(L, I)*PHF_TMP
                  PHASE_FNC_SOLAR_CLOUD(L, I, ID)                       &
     &              =PHASE_FNC_SOLAR_CLOUD(L, I, ID)                    &
     &              +KS_PHF(L)*P_LEGENDRE_LS(L)                         &
     &              *REAL(2*LS+1, Real64)
                ENDDO


              ENDDO
            ENDDO
!
!           Continue to an extra order to find the rescaling
!           for the solar beam.
            IF (L_RESCALE) THEN
              LS=N_ORDER_PHASE_SOLAR+1
              DO LL=1, N_CLOUD_PROFILE(I)
                L=I_CLOUD_PROFILE(LL, I)
!
                SELECT CASE(I_PARAMETRIZATION_ICE)
!
                  CASE(IP_SLINGO_SCHR_ICE_PHF)
                    PHF_TMP=ICE_CLOUD_PARAMETER(2*LS+3)                 &
     &                +DIM_CHAR_ICE(L, I)                               &
     &                *ICE_CLOUD_PARAMETER(2*LS+4)
!
                  CASE(IP_ICE_FU_PHF)
                    X=DIM_CHAR_ICE(L, I)/ICE_CLOUD_PARAMETER(5*LS+9)
                    PHF_TMP                                             &
     &                =(ICE_CLOUD_PARAMETER(5*LS+5)                     &
     &                +X*(ICE_CLOUD_PARAMETER(5*LS+6)                   &
     &                +X*(ICE_CLOUD_PARAMETER(5*LS+7)                   &
     &                +X*ICE_CLOUD_PARAMETER(5*LS+8))))
!
                END SELECT
!
                FORWARD_SOLAR_CLOUD(L, I)                               &
     &            =K_EXT_SCAT_CLOUD(L, I)*PHF_TMP
              ENDDO
            ENDIF
!
          ENDIF
!
        ENDDO
!
      ENDIF
!
#if defined(OBSERVED)
!
      IF (I_PARAMETRIZATION_ICE == IP_ICE_UNPARAMETRIZED) THEN
!
        CALL PRSC_OPT_PROP(IERR                                         &
     &    , N_PROFILE, N_CLOUD_TOP, N_LAYER                             &
     &    , L_RESCALE, N_ORDER_FORWARD                                  &
     &    , L_HENYEY_GREENSTEIN_PF, N_ORDER_PHASE                       &
     &    , P, DENSITY                                                  &
     &    , N_OPT_LEVEL_CLOUD_PRSC                                      &
     &    , ICE_PRESSURE_PRSC, ICE_ABSORPTION_PRSC, ICE_SCATTERING_PRSC &
     &    , ICE_PHASE_FNC_PRSC                                          &
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
      IF ( (I_PARAMETRIZATION_ICE /= IP_SLINGO_SCHRECKER_ICE).AND.      &
     &     (I_PARAMETRIZATION_ICE /= IP_ICE_ADT).AND.                   &
     &     (I_PARAMETRIZATION_ICE /= IP_ICE_ADT_10).AND.                &
     &     (I_PARAMETRIZATION_ICE /= IP_SUN_SHINE_VN2_VIS).AND.         &
     &     (I_PARAMETRIZATION_ICE /= IP_SUN_SHINE_VN2_IR).AND.          &
     &     (I_PARAMETRIZATION_ICE /= IP_ICE_FU_SOLAR).AND.              &
     &     (I_PARAMETRIZATION_ICE /= IP_ICE_FU_IR).AND.                 &
     &     (I_PARAMETRIZATION_ICE /= IP_SLINGO_SCHR_ICE_PHF).AND.       &
     &     (I_PARAMETRIZATION_ICE /= IP_ICE_FU_PHF).AND.                &
     &     (I_PARAMETRIZATION_ICE /= IP_ICE_UNPARAMETRIZED) ) THEN
!
        WRITE(IU_ERR, '(/A)') '*** Error: An invalid parametrization '  &
     &    //'of ice crystals has been used.'
        IERR=I_ERR_FATAL
        RETURN
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE OPT_PROP_ICE_CLOUD
#endif
