#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate grey optical properties.
!
! Method:
!       For each activated optical process, excluding gaseous
!       absorption, increments are calculated for the total and
!       scattering extinctions, and the products of the asymmetry
!       factor and the forward scattering factor in clear and
!       cloudy regions. These increments are summed, and the grey
!       total and scattering extinctions and the asymmetry and forward
!       scattering factors are thus calculated.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!DD+ -------------------------------------------------------------------
! Compiler directives for specific computer systems:
!
! Indirect addressing will inhibit vectorization, but here it is safe
! because all points are independent.
!
! Fujistu VPP700:
!OCL NOVREC
!
! Cray vector machines:
!fpp$ NODEPCHK R
!
!DD- -------------------------------------------------------------------
      SUBROUTINE GREY_OPT_PROP(IERR                                     &
     &  , N_PROFILE, N_LAYER, P, T, DENSITY                             &
     &  , N_ORDER_PHASE, L_RESCALE, N_ORDER_FORWARD                     &
     &  , L_HENYEY_GREENSTEIN_PF, L_SOLAR_PHF, N_ORDER_PHASE_SOLAR      &
     &  , N_DIRECTION, COS_SOL_VIEW                                     &
     &  , L_RAYLEIGH, RAYLEIGH_COEFF                                    &
     &  , L_CONTINUUM, N_CONTINUUM, I_CONTINUUM_POINTER, K_CONTINUUM    &
     &  , AMOUNT_CONTINUUM                                              &
     &  , L_AEROSOL, N_AEROSOL, AEROSOL_MIX_RATIO                       &
     &  , I_AEROSOL_PARAMETRIZATION                                     &
     &  , I_HUMIDITY_POINTER, HUMIDITIES, DELTA_HUMIDITY                &
     &  , MEAN_REL_HUMIDITY                                             &
     &  , AEROSOL_ABSORPTION, AEROSOL_SCATTERING, AEROSOL_PHASE_FNC     &
#if defined(OBSERVED)
     &  , N_OPT_LEVEL_AEROSOL_PRSC, AEROSOL_PRESSURE_PRSC               &
     &  , AEROSOL_ABSORPTION_PRSC, AEROSOL_SCATTERING_PRSC              &
     &  , AEROSOL_PHASE_FNC_PRSC                                        &
#endif
     &  , L_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE, N_CLOUD_TOP        &
     &  , N_CONDENSED, L_CLOUD_CMP, I_PHASE_CMP                         &
     &  , I_CONDENSED_PARAM, CONDENSED_N_PHF, CONDENSED_PARAM_LIST      &
     &  , CONDENSED_MIX_RATIO, CONDENSED_DIM_CHAR                       &
     &  , N_CLOUD_TYPE, I_CLOUD_TYPE                                    &
#if defined(OBSERVED)
     &  , N_OPT_LEVEL_DROP_PRSC                                         &
     &  , DROP_PRESSURE_PRSC, DROP_ABSORPTION_PRSC                      &
     &  , DROP_SCATTERING_PRSC, DROP_PHASE_FNC_PRSC                     &
     &  , N_OPT_LEVEL_ICE_PRSC, ICE_PRESSURE_PRSC                       &
     &  , ICE_ABSORPTION_PRSC, ICE_SCATTERING_PRSC, ICE_PHASE_FNC_PRSC  &
#endif
     &  , K_EXT_TOT_CLR, K_EXT_SCAT_CLR, PHASE_FNC_CLR                  &
     &  , FORWARD_SCATTER_CLR, FORWARD_SOLAR_CLR, PHASE_FNC_SOLAR_CLR   &
     &  , K_EXT_TOT, K_EXT_SCAT, PHASE_FNC, FORWARD_SCATTER             &
     &  , FORWARD_SOLAR, PHASE_FNC_SOLAR                                &
     &  , ND_PROFILE, ND_RADIANCE_PROFILE, ND_LAYER                     &
     &  , ND_LAYER_CLR, ID_CT                                           &
     &  , ND_CONTINUUM, ND_AEROSOL_SPECIES, ND_HUMIDITIES               &
     &  , ND_CLOUD_PARAMETER, ND_CLOUD_COMPONENT, ND_CLOUD_TYPE         &
     &  , ND_PHASE_TERM, ND_MAX_ORDER, ND_DIRECTION                     &
#if defined(OBSERVED)
     &  , ND_PROFILE_AEROSOL_PRSC, ND_PROFILE_CLOUD_PRSC                &
     &  , ND_OPT_LEVEL_AEROSOL_PRSC, ND_OPT_LEVEL_CLOUD_PRSC            &
#endif
     &  )
!
!
!
      IMPLICIT NONE
!
! Include Header Files
#include "c_kinds.h"
!
!     Sizes of dummy arrays
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for profiles
     &  , ND_RADIANCE_PROFILE                                           &
!           Size allocated for profiles of quantities specifically
!           used in calulating radiances
     &  , ND_LAYER                                                      &
!           Size allocated for layers
     &  , ND_LAYER_CLR                                                  &
!           Size allocated for completely clear layers
     &  , ID_CT                                                         &
!           Topmost declared cloudy layer
     &  , ND_DIRECTION                                                  &
!           Size allocated for viewing directions
     &  , ND_AEROSOL_SPECIES                                            &
!           Size allocated for aerosols
     &  , ND_HUMIDITIES                                                 &
!           Size allocated for humidities
     &  , ND_CONTINUUM                                                  &
!           Size allocated for continua
     &  , ND_PHASE_TERM                                                 &
!           Size allocated for terms in the phase function
     &  , ND_MAX_ORDER                                                  &
!           Size allocated for the order of the calculation
     &  , ND_CLOUD_PARAMETER                                            &
!           Size allocated for cloud parameters
     &  , ND_CLOUD_COMPONENT                                            &
!           Size allocated for components of clouds
     &  , ND_CLOUD_TYPE
!           Size allocated for types of clouds
#if defined(OBSERVED)
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE_AEROSOL_PRSC                                       &
!           Size allocated for profiles of prescribed
!           cloudy optical properties
     &  , ND_PROFILE_CLOUD_PRSC                                         &
!           Size allocated for profiles of prescribed
!           aerosol optical properties
     &  , ND_OPT_LEVEL_AEROSOL_PRSC                                     &
!           Size allocated for levels of prescribed
!           cloudy optical properties
     &  , ND_OPT_LEVEL_CLOUD_PRSC
!           Size allocated for levels of prescribed
!           aerosol optical properties
#endif
!
!     Inclusion of header files.
#include "def_std_io_icf3z.h"
#include "aerosol_parametrization_pcf3z.h"
#include "cloud_scheme_pcf3z.h"
#include "phase_pcf3z.h"
#include "error_pcf3z.h"
!
!
!
!     Dummy arguments.
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
!
!
!     Basic atmospheric properties:
!
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER
!           Number of layers
!
      REAL  (Real64), INTENT(IN) ::                                     &
     &    P(ND_PROFILE, ND_LAYER)                                       &
!           Pressure
     &  , T(ND_PROFILE, ND_LAYER)                                       &
!           Temperature
     &  , DENSITY(ND_PROFILE, ND_LAYER)
!           Density at levels
!
!
!     Optical switches:
      LOGICAL, INTENT(IN) ::                                            &
     &    L_RESCALE                                                     &
!           Delta-rescaling required
     &  , L_HENYEY_GREENSTEIN_PF                                        &
!           Flag to use a Henyey-Greenstein phase function
     &  , L_SOLAR_PHF
!           Flag to use an extended phase function for solar
!           radiation
      INTEGER, INTENT(IN) ::                                            &
     &    N_ORDER_PHASE                                                 &
!           Order of terms in the phase function
     &  , N_ORDER_PHASE_SOLAR                                           &
!           Order of truncation of the solar beam
     &  , N_ORDER_FORWARD
!           Order used in forming the forward scattering parameter
!
!     Directional information
      INTEGER, INTENT(IN) ::                                            &
     &    N_DIRECTION
!           Number of viewing directions
      REAL  (Real64), INTENT(IN) ::                                     &
     &    COS_SOL_VIEW(ND_RADIANCE_PROFILE, ND_DIRECTION)
!           Cosines of the angles between the solar direction
!           and the viewing direction
!
!
!     Rayleigh scattering:
!
      LOGICAL, INTENT(IN) ::                                            &
     &    L_RAYLEIGH
!           Rayleigh scattering activated
!
      REAL  (Real64), INTENT(IN) ::                                     &
     &    RAYLEIGH_COEFF
!           Rayleigh coefficient
!
!
!     Continuum processes:
      LOGICAL, INTENT(IN) ::                                            &
     &    L_CONTINUUM
!           Continuum absorption activated
!
      INTEGER, INTENT(IN) ::                                            &
     &    N_CONTINUUM                                                   &
!           Number of continua
     &  , I_CONTINUUM_POINTER(ND_CONTINUUM)
!           Pointers to active continua
!
      REAL  (Real64), INTENT(IN) ::                                     &
     &    K_CONTINUUM(ND_CONTINUUM)                                     &
!           Continuum extinction
     &  , AMOUNT_CONTINUUM(ND_PROFILE, ND_LAYER, ND_CONTINUUM)
!           Amounts for continua
!
!
!     Properties of aerosols:
!
      LOGICAL, INTENT(IN) ::                                            &
     &    L_AEROSOL
!           Aerosols activated
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
     &    AEROSOL_MIX_RATIO(ND_PROFILE, ND_LAYER                        &
     &      , ND_AEROSOL_SPECIES)                                       &
!           Number densty of aerosols
     &  , AEROSOL_ABSORPTION(ND_HUMIDITIES, ND_AEROSOL_SPECIES)         &
!           Aerosol absorption in band for a mixing ratio of unity
     &  , AEROSOL_SCATTERING(ND_HUMIDITIES, ND_AEROSOL_SPECIES)         &
!           Aerosol scattering in band for a mixing ratio of unity
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
!     Observational properties of aerosols:
      INTEGER, INTENT(IN) ::                                            &
     &    N_OPT_LEVEL_AEROSOL_PRSC(ND_AEROSOL_SPECIES)
!           Number of levels of prescribed optical properties
!           of aerosols
!
      REAL  (Real64), INTENT(IN) ::                                     &
     &    AEROSOL_PRESSURE_PRSC(ND_PROFILE_AEROSOL_PRSC                 &
     &      , ND_OPT_LEVEL_AEROSOL_PRSC, ND_AEROSOL_SPECIES)            &
!           Pressures at which optical properties of aerosols
!           are prescribed
     &  , AEROSOL_ABSORPTION_PRSC(ND_PROFILE_AEROSOL_PRSC               &
     &      , ND_OPT_LEVEL_AEROSOL_PRSC, ND_AEROSOL_SPECIES)            &
!           Prescribed absorption by aerosols
     &  , AEROSOL_SCATTERING_PRSC(ND_PROFILE_AEROSOL_PRSC               &
     &      , ND_OPT_LEVEL_AEROSOL_PRSC, ND_AEROSOL_SPECIES)            &
!           Prescribed scattering by aerosols
     &  , AEROSOL_PHASE_FNC_PRSC(ND_PROFILE_AEROSOL_PRSC                &
     &      , ND_OPT_LEVEL_AEROSOL_PRSC                                 &
     &      , ND_PHASE_TERM, ND_AEROSOL_SPECIES)
!           Prescribed phase functions of aerosols
#endif
!
!
!     Properties of clouds:
!
      LOGICAL, INTENT(IN) ::                                            &
     &    L_CLOUD
!           Clouds activated
!
!     Geometry of clouds:
!
      INTEGER, INTENT(IN) ::                                            &
     &    N_CLOUD_TOP                                                   &
!           Topmost cloudy layer
     &  , N_CLOUD_TYPE                                                  &
!           Number of types of clouds
     &  , N_CLOUD_PROFILE(ID_CT: ND_LAYER)                              &
!           Number of cloudy profiles in each layer
     &  , I_CLOUD_PROFILE(ND_PROFILE, ID_CT: ND_LAYER)                  &
!           Profiles containing clouds
     &  , I_CLOUD_TYPE(ND_CLOUD_COMPONENT)
!           Types of cloud to which each component contributes
!
!     Microphysical quantities:
      INTEGER, INTENT(IN) ::                                            &
     &    N_CONDENSED                                                   &
!           Number of condensed components
     &  , I_PHASE_CMP(ND_CLOUD_COMPONENT)                               &
!           Phases of cloudy components
     &  , I_CONDENSED_PARAM(ND_CLOUD_COMPONENT)                         &
!           Parametrization schemes for cloudy components
     &  , CONDENSED_N_PHF(ND_CLOUD_COMPONENT)
!           Number of terms in the phase function for each
!           cloudy component
!
      LOGICAL, INTENT(IN) ::                                            &
     &    L_CLOUD_CMP(ND_CLOUD_COMPONENT)
!           Flags to activate cloudy components
!
      REAL  (Real64), INTENT(IN) ::                                     &
     &    CONDENSED_PARAM_LIST(ND_CLOUD_PARAMETER                       &
     &      , ND_CLOUD_COMPONENT)                                       &
!           Coefficients in parametrization schemes
     &  , CONDENSED_MIX_RATIO(ND_PROFILE, ID_CT: ND_LAYER               &
     &      , ND_CLOUD_COMPONENT)                                       &
!           Mixing ratios of cloudy components
     &  , CONDENSED_DIM_CHAR(ND_PROFILE, ID_CT: ND_LAYER                &
     &      , ND_CLOUD_COMPONENT)
!           Characteristic dimensions of cloudy components
!
#if defined(OBSERVED)
!     Prescribed cloudy optical properties:
      INTEGER, INTENT(IN) ::                                            &
     &    N_OPT_LEVEL_DROP_PRSC                                         &
!           Number of levels of prescribed
!           optical properties of droplets
     &  , N_OPT_LEVEL_ICE_PRSC
!           Number of levels of prescribed
!           optical properties of ice crystals
!
      REAL  (Real64), INTENT(IN) ::                                     &
     &    DROP_PRESSURE_PRSC(ND_PROFILE_CLOUD_PRSC                      &
     &      , ND_OPT_LEVEL_CLOUD_PRSC)                                  &
!           Pressures at which optical properties of
!           droplets are prescribed
     &  , DROP_ABSORPTION_PRSC(ND_PROFILE_CLOUD_PRSC                    &
     &      , ND_OPT_LEVEL_CLOUD_PRSC)                                  &
!           Prescribed absorption by droplets
     &  , DROP_SCATTERING_PRSC(ND_PROFILE_CLOUD_PRSC                    &
     &      , ND_OPT_LEVEL_CLOUD_PRSC)                                  &
!           Prescribed scattering by droplets
     &  , DROP_PHASE_FNC_PRSC(ND_PROFILE_CLOUD_PRSC                     &
     &      , ND_OPT_LEVEL_CLOUD_PRSC, ND_PHASE_TERM)                   &
!           Prescribed phase function of droplets
     &  , ICE_PRESSURE_PRSC(ND_PROFILE_CLOUD_PRSC                       &
     &      , ND_OPT_LEVEL_CLOUD_PRSC)                                  &
!           Pressures at which optical properties of
!           ice crystals are prescribed
     &  , ICE_ABSORPTION_PRSC(ND_PROFILE_CLOUD_PRSC                     &
     &      , ND_OPT_LEVEL_CLOUD_PRSC)                                  &
!           Prescribed absorption by ice crystals
     &  , ICE_SCATTERING_PRSC(ND_PROFILE_CLOUD_PRSC                     &
     &      , ND_OPT_LEVEL_CLOUD_PRSC)                                  &
!           Prescribed scattering by ice crystals
     &  , ICE_PHASE_FNC_PRSC(ND_PROFILE_CLOUD_PRSC                      &
     &      , ND_OPT_LEVEL_CLOUD_PRSC, ND_PHASE_TERM)
!           Prescribed phase functions of ice crystals
#endif
!
!
!     Calculated optical properties:
!
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    K_EXT_SCAT_CLR(ND_PROFILE, ND_LAYER_CLR)                      &
!           Clear_sky scattering extinction
     &  , K_EXT_TOT_CLR(ND_PROFILE, ND_LAYER_CLR)                       &
!           Total clear_sky extinction
     &  , PHASE_FNC_CLR(ND_PROFILE, ND_LAYER_CLR, ND_MAX_ORDER)         &
!           Clear-sky phase function
     &  , PHASE_FNC_SOLAR_CLR(ND_RADIANCE_PROFILE                       &
     &      , ND_LAYER_CLR, ND_DIRECTION)                               &
!           Clear-sky phase function for the solar beam
     &  , FORWARD_SCATTER_CLR(ND_PROFILE, ND_LAYER_CLR)                 &
!           Clear-sky forward scattering
     &  , FORWARD_SOLAR_CLR(ND_PROFILE, ND_LAYER_CLR)                   &
!           Clear-sky forward scattering for the solar beam
     &  , K_EXT_SCAT(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)     &
!           Scattering extinction
     &  , K_EXT_TOT(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)      &
!           Total extinction
     &  , PHASE_FNC(ND_PROFILE, ID_CT: ND_LAYER, ND_MAX_ORDER           &
     &      , 0: ND_CLOUD_TYPE)                                         &
!           Phase function
     &  , PHASE_FNC_SOLAR(ND_RADIANCE_PROFILE, ID_CT: ND_LAYER          &
     &      , ND_DIRECTION, 0: ND_CLOUD_TYPE)                           &
!           Phase function for the solar beam
     &  , FORWARD_SCATTER(ND_PROFILE, ID_CT: ND_LAYER                   &
     &      , 0: ND_CLOUD_TYPE)                                         &
!           Forward scattering
     &  , FORWARD_SOLAR(ND_PROFILE, ID_CT: ND_LAYER                     &
     &      , 0: ND_CLOUD_TYPE)
!           Forward scattering for the solar beam
!
!
!
!     Local variables.
      INTEGER                                                           &
     &    I_CONTINUUM                                                   &
!           Temporary continuum `index'
     &  , L                                                             &
!           Loop variable
     &  , LL                                                            &
!           Loop variable
     &  , I                                                             &
!           Loop variable
     &  , ID                                                            &
!           Loop variable
     &  , J                                                             &
!           Loop variable
     &  , K                                                             &
!           Loop variable
     &  , LS                                                            &
!           Loop variable
     &  , N_INDEX                                                       &
!           Number of indices satisfying the test
     &  , INDEX(ND_PROFILE)
!           Indices satifying the test
!
!     Temporary optical properties:
!
      REAL  (Real64) ::                                                 &
     &    K_EXT_SCAT_CLOUD_COMP(ND_PROFILE, ID_CT: ND_LAYER)            &
!           Scattering extinction of cloudy component
     &  , K_EXT_TOT_CLOUD_COMP(ND_PROFILE, ID_CT: ND_LAYER)             &
!           Total extinction of cloudy component
     &  , PHASE_FNC_CLOUD_COMP(ND_PROFILE, ID_CT: ND_LAYER              &
     &      , ND_MAX_ORDER)                                             &
!           Phase function of cloudy components
     &  , PHASE_FNC_SOLAR_CLOUD_COMP(ND_RADIANCE_PROFILE                &
     &      , ID_CT: ND_LAYER, ND_DIRECTION)                            &
!           Phase function of cloudy components for singly scattered
!           solar radiation
     &  , FORWARD_SCATTER_CLOUD_COMP(ND_PROFILE, ID_CT: ND_LAYER)       &
!           Forward scattering of cloudy component
     &  , FORWARD_SOLAR_CLOUD_COMP(ND_PROFILE, ID_CT: ND_LAYER)
!           Forward scattering for the solar beam
!           in the cloudy component
!
!     Subroutines called:
      EXTERNAL                                                          &
     &    OPT_PROP_AEROSOL, OPT_PROP_WATER_CLOUD, OPT_PROP_ICE_CLOUD
!
!
!
!     Initialize the extinction coefficients and the phase function.
      DO I=1, N_CLOUD_TOP-1
        DO L=1, N_PROFILE
          K_EXT_TOT_CLR(L, I)=0.0E+00_Real64
          K_EXT_SCAT_CLR(L, I)=0.0E+00_Real64
        ENDDO
      ENDDO
      DO I=N_CLOUD_TOP, N_LAYER
        DO L=1, N_PROFILE
          K_EXT_TOT(L, I, 0)=0.0E+00_Real64
          K_EXT_SCAT(L, I, 0)=0.0E+00_Real64
        ENDDO
      ENDDO
      DO LS=1, N_ORDER_PHASE
        DO I=1, N_CLOUD_TOP-1
          DO L=1, N_PROFILE
            PHASE_FNC_CLR(L, I, LS)=0.0E+00_Real64
          ENDDO
        ENDDO
        DO I=N_CLOUD_TOP, N_LAYER
          DO L=1, N_PROFILE
            PHASE_FNC(L, I, LS, 0)=0.0E+00_Real64
          ENDDO
        ENDDO
      ENDDO
!     Forward scattering is required only when delta-rescaling
!     is performed.
      IF (L_RESCALE) THEN
        DO I=1, N_CLOUD_TOP-1
          DO L=1, N_PROFILE
            FORWARD_SCATTER_CLR(L, I)=0.0E+00_Real64
          ENDDO
        ENDDO
        DO I=N_CLOUD_TOP, N_LAYER
          DO L=1, N_PROFILE
            FORWARD_SCATTER(L, I, 0)=0.0E+00_Real64
          ENDDO
        ENDDO
      ENDIF
!     If using a separate solar phase function that must be
!     initialized.
      IF (L_SOLAR_PHF) THEN
        DO ID=1, N_DIRECTION
          DO I=1, N_CLOUD_TOP-1
            DO L=1, N_PROFILE
              PHASE_FNC_SOLAR_CLR(L, I, ID)=0.0E+00_Real64
            ENDDO
          ENDDO
          DO I=N_CLOUD_TOP, N_LAYER
            DO L=1, N_PROFILE
              PHASE_FNC_SOLAR(L, I, ID, 0)=0.0E+00_Real64
            ENDDO
          ENDDO
        ENDDO
        DO I=1, N_CLOUD_TOP-1
          DO L=1, N_PROFILE
            FORWARD_SOLAR_CLR(L, I)=0.0E+00_Real64
          ENDDO
        ENDDO
        DO I=N_CLOUD_TOP, N_LAYER
          DO L=1, N_PROFILE
            FORWARD_SOLAR(L, I, 0)=0.0E+00_Real64
          ENDDO
        ENDDO
      ENDIF
!
!
!
!
!
!     Consider each optical process in turn.
!
!     Rayleigh scattering:
!
      IF (L_RAYLEIGH) THEN
        DO I=1, N_CLOUD_TOP-1
          DO L=1, N_PROFILE
            K_EXT_SCAT_CLR(L, I)                                        &
     &        =K_EXT_SCAT_CLR(L, I)+RAYLEIGH_COEFF
          ENDDO
        ENDDO
        DO I=N_CLOUD_TOP, N_LAYER
          DO L=1, N_PROFILE
            K_EXT_SCAT(L, I, 0)                                         &
     &        =K_EXT_SCAT(L, I, 0)+RAYLEIGH_COEFF
          ENDDO
        ENDDO
!
!       Only the second Lengendre polynomial contributes.
        IF (N_ORDER_PHASE >= 2) THEN
          DO I=1, N_CLOUD_TOP-1
            DO L=1, N_PROFILE
              PHASE_FNC_CLR(L, I, 2)=PHASE_FNC_CLR(L, I, 2)             &
     &          +RAYLEIGH_COEFF*1.0E-01_Real64
            ENDDO
          ENDDO
          DO I=N_CLOUD_TOP, N_LAYER
            DO L=1, N_PROFILE
              PHASE_FNC(L, I, 2, 0)=PHASE_FNC(L, I, 2, 0)               &
     &          +RAYLEIGH_COEFF*1.0E-01_Real64
            ENDDO
          ENDDO
        ENDIF
!
!       No formal rescaling is applied to the phase function for
!       Rayleigh scattering, as only g_2 is non-zero.
!
        IF (L_SOLAR_PHF) THEN
!
          DO ID=1, N_DIRECTION
            DO I=1, N_CLOUD_TOP-1
              DO L=1, N_PROFILE
                PHASE_FNC_SOLAR_CLR(L, I, ID)                           &
     &           =PHASE_FNC_SOLAR_CLR(L, I, ID)                         &
     &           +RAYLEIGH_COEFF                                        &
     &           *0.75E+00_Real64*(1.0E+00_Real64+COS_SOL_VIEW(L,ID)**2)
              ENDDO
            ENDDO
            DO I=N_CLOUD_TOP, N_LAYER
              DO L=1, N_PROFILE
                PHASE_FNC_SOLAR(L, I, ID, 0)                            &
     &           =PHASE_FNC_SOLAR(L, I, ID, 0)                          &
     &           +RAYLEIGH_COEFF                                        &
     &           *0.75E+00_Real64*(1.0E+00_Real64+COS_SOL_VIEW(L,ID)**2)
              ENDDO
            ENDDO
          ENDDO
!
        ENDIF
!
      ENDIF
!
      IF (L_AEROSOL) THEN
!       Include the effects of aerosol.
!       Above clouds.
! DEPENDS ON: opt_prop_aerosol
        CALL OPT_PROP_AEROSOL(IERR                                      &
     &    , N_PROFILE, 1, N_CLOUD_TOP-1                                 &
     &    , N_ORDER_PHASE, L_RESCALE, N_ORDER_FORWARD                   &
     &    , L_HENYEY_GREENSTEIN_PF                                      &
     &    , N_AEROSOL, AEROSOL_MIX_RATIO                                &
     &    , I_AEROSOL_PARAMETRIZATION                                   &
     &    , I_HUMIDITY_POINTER, HUMIDITIES, DELTA_HUMIDITY              &
     &    , MEAN_REL_HUMIDITY                                           &
     &    , AEROSOL_ABSORPTION, AEROSOL_SCATTERING, AEROSOL_PHASE_FNC   &
     &    , L_SOLAR_PHF, N_ORDER_PHASE_SOLAR, N_DIRECTION, COS_SOL_VIEW &
#if defined(OBSERVED)
     &    , P, DENSITY                                                  &
     &    , N_OPT_LEVEL_AEROSOL_PRSC, AEROSOL_PRESSURE_PRSC             &
     &    , AEROSOL_ABSORPTION_PRSC, AEROSOL_SCATTERING_PRSC            &
     &    , AEROSOL_PHASE_FNC_PRSC                                      &
#endif
     &    , K_EXT_TOT_CLR, K_EXT_SCAT_CLR                               &
     &    , PHASE_FNC_CLR, FORWARD_SCATTER_CLR                          &
     &    , FORWARD_SOLAR_CLR, PHASE_FNC_SOLAR_CLR                      &
     &    , ND_PROFILE, ND_RADIANCE_PROFILE, ND_LAYER                   &
     &    , 1, ND_LAYER_CLR                                             &
     &    , ND_AEROSOL_SPECIES, ND_HUMIDITIES                           &
     &    , ND_PHASE_TERM, ND_MAX_ORDER, ND_DIRECTION                   &
#if defined(OBSERVED)
     &    , ND_PROFILE_AEROSOL_PRSC, ND_OPT_LEVEL_AEROSOL_PRSC          &
#endif
     &    )
!       Within clouds:
! DEPENDS ON: opt_prop_aerosol
        CALL OPT_PROP_AEROSOL(IERR                                      &
     &    , N_PROFILE, N_CLOUD_TOP, N_LAYER                             &
     &    , N_ORDER_PHASE, L_RESCALE, N_ORDER_FORWARD                   &
     &    , L_HENYEY_GREENSTEIN_PF                                      &
     &    , N_AEROSOL, AEROSOL_MIX_RATIO                                &
     &    , I_AEROSOL_PARAMETRIZATION                                   &
     &    , I_HUMIDITY_POINTER, HUMIDITIES, DELTA_HUMIDITY              &
     &    , MEAN_REL_HUMIDITY                                           &
     &    , AEROSOL_ABSORPTION, AEROSOL_SCATTERING, AEROSOL_PHASE_FNC   &
     &    , L_SOLAR_PHF, N_ORDER_PHASE_SOLAR, N_DIRECTION, COS_SOL_VIEW &
#if defined(OBSERVED)
     &    , P, DENSITY                                                  &
     &    , N_OPT_LEVEL_AEROSOL_PRSC, AEROSOL_PRESSURE_PRSC             &
     &    , AEROSOL_ABSORPTION_PRSC, AEROSOL_SCATTERING_PRSC            &
     &    , AEROSOL_PHASE_FNC_PRSC                                      &
#endif
     &    , K_EXT_TOT(1, ID_CT, 0), K_EXT_SCAT(1, ID_CT, 0)             &
     &    , PHASE_FNC(1, ID_CT, 1, 0), FORWARD_SCATTER(1, ID_CT, 0)     &
     &    , FORWARD_SOLAR(1, ID_CT, 0), PHASE_FNC_SOLAR(1, ID_CT, 1, 0) &
     &    , ND_PROFILE, ND_RADIANCE_PROFILE, ND_LAYER                   &
     &    , ID_CT, ND_LAYER                                             &
     &    , ND_AEROSOL_SPECIES, ND_HUMIDITIES                           &
     &    , ND_PHASE_TERM, ND_MAX_ORDER, ND_DIRECTION                   &
#if defined(OBSERVED)
     &    , ND_PROFILE_AEROSOL_PRSC, ND_OPT_LEVEL_AEROSOL_PRSC          &
#endif
     &    )
      ENDIF
!
      IF (L_CONTINUUM) THEN
!       Include continuum absorption.
        DO J=1, N_CONTINUUM
          I_CONTINUUM=I_CONTINUUM_POINTER(J)
          DO I=1, N_CLOUD_TOP-1
            DO L=1, N_PROFILE
              K_EXT_TOT_CLR(L, I)=K_EXT_TOT_CLR(L, I)                   &
     &          +K_CONTINUUM(I_CONTINUUM)                               &
     &          *AMOUNT_CONTINUUM(L, I, I_CONTINUUM)
            ENDDO
          ENDDO
          DO I=N_CLOUD_TOP, N_LAYER
            DO L=1, N_PROFILE
              K_EXT_TOT(L, I, 0)=K_EXT_TOT(L, I, 0)                     &
     &          +K_CONTINUUM(I_CONTINUUM)                               &
     &          *AMOUNT_CONTINUUM(L, I, I_CONTINUUM)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!
!
!     Add the scattering on to the total extinction. The final clear-sky
!     phase function not calculated here since the product of the phase
!     function and scattering is also needed to calculate the cloudy
!     phase function.
      DO I=1, N_CLOUD_TOP-1
        DO L=1, N_PROFILE
          K_EXT_TOT_CLR(L, I)=K_EXT_TOT_CLR(L, I)+K_EXT_SCAT_CLR(L, I)
        ENDDO
      ENDDO
      DO I=N_CLOUD_TOP, N_LAYER
        DO L=1, N_PROFILE
          K_EXT_TOT(L, I, 0)=K_EXT_TOT(L, I, 0)+K_EXT_SCAT(L, I, 0)
        ENDDO
      ENDDO
!
!
!     If there are no clouds calculate the final optical properties
!     and return to the calling routine.
!
      IF (.NOT.L_CLOUD) THEN
!
        DO LS=1, N_ORDER_PHASE
          DO I=1, N_CLOUD_TOP-1
            DO L=1, N_PROFILE
              IF (K_EXT_SCAT_CLR(L, I) >  0.0E+00_Real64) THEN
                PHASE_FNC_CLR(L, I, LS)=PHASE_FNC_CLR(L, I, LS)         &
     &            /K_EXT_SCAT_CLR(L, I)
              ENDIF
            ENDDO
          ENDDO
          DO I=N_CLOUD_TOP, N_LAYER
            DO L=1, N_PROFILE
              IF (K_EXT_SCAT(L, I, 0) >  0.0E+00_Real64) THEN
                PHASE_FNC(L, I, LS, 0)=PHASE_FNC(L, I, LS, 0)           &
     &            /K_EXT_SCAT(L, I, 0)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
!
        IF (L_RESCALE) THEN
          DO I=1, N_CLOUD_TOP-1
            DO L=1, N_PROFILE
              IF (K_EXT_SCAT_CLR(L, I) >  0.0E+00_Real64) THEN
                FORWARD_SCATTER_CLR(L, I)                               &
     &            =FORWARD_SCATTER_CLR(L, I)                            &
     &            /K_EXT_SCAT_CLR(L, I)
              ENDIF
            ENDDO
          ENDDO
          DO I=N_CLOUD_TOP, N_LAYER
            DO L=1, N_PROFILE
              IF (K_EXT_SCAT(L, I, 0) >  0.0E+00_Real64) THEN
                FORWARD_SCATTER(L, I, 0)                                &
     &            =FORWARD_SCATTER(L, I, 0)                             &
     &            /K_EXT_SCAT(L, I, 0)
              ENDIF
            ENDDO
          ENDDO
!
        ENDIF
!
        IF (L_SOLAR_PHF) THEN
!
          DO I=1, N_CLOUD_TOP-1
            DO L=1, N_PROFILE
              IF (K_EXT_SCAT_CLR(L, I) >  0.0E+00_Real64) THEN
                FORWARD_SOLAR_CLR(L, I)                                 &
     &            =FORWARD_SOLAR_CLR(L, I)                              &
     &            /K_EXT_SCAT_CLR(L, I)
              ENDIF
            ENDDO
          ENDDO
          DO I=N_CLOUD_TOP, N_LAYER
            DO L=1, N_PROFILE
              IF (K_EXT_SCAT(L, I, 0) >  0.0E+00_Real64) THEN
                FORWARD_SOLAR(L, I, 0)                                  &
     &            =FORWARD_SOLAR(L, I, 0)                               &
     &            /K_EXT_SCAT(L, I, 0)
              ENDIF
            ENDDO
          ENDDO
!
          DO ID=1, N_DIRECTION
            DO I=1, N_CLOUD_TOP-1
              DO L=1, N_PROFILE
                IF (K_EXT_SCAT_CLR(L, I) >  0.0E+00_Real64)             &
     &            PHASE_FNC_SOLAR_CLR(L, I, ID)                         &
     &              =PHASE_FNC_SOLAR_CLR(L, I, ID)                      &
     &              /K_EXT_SCAT_CLR(L, I)
              ENDDO
            ENDDO
            DO I=N_CLOUD_TOP, N_LAYER
              DO L=1, N_PROFILE
                IF (K_EXT_SCAT(L, I, 0) >  0.0E+00_Real64)              &
     &            PHASE_FNC_SOLAR(L, I, ID, 0)                          &
     &              =PHASE_FNC_SOLAR(L, I, ID, 0)                       &
     &              /K_EXT_SCAT(L, I, 0)
              ENDDO
            ENDDO
          ENDDO
!
        ENDIF
!
        RETURN
!
      ENDIF
!
!
!
!
!     Addition of cloudy properties:
!
!
!     Add in background contibutions:
!
!
!     All the processes occurring outside clouds also occur
!     within them.
      DO K=1, N_CLOUD_TYPE
        DO I=N_CLOUD_TOP, N_LAYER
          DO L=1, N_PROFILE
            K_EXT_TOT(L, I, K)=K_EXT_TOT(L, I, 0)
            K_EXT_SCAT(L, I, K)=K_EXT_SCAT(L, I, 0)
            FORWARD_SCATTER(L, I, K)=FORWARD_SCATTER(L, I, 0)
            FORWARD_SOLAR(L, I, K)=FORWARD_SOLAR(L, I, 0)
          ENDDO
          DO LS=1, N_ORDER_PHASE
            DO L=1, N_PROFILE
              PHASE_FNC(L, I, LS, K)=PHASE_FNC(L, I, LS, 0)
            ENDDO
          ENDDO
        ENDDO
!       If using a separate solar phase function that must
!       be initialized.
        IF (L_SOLAR_PHF) THEN
          DO ID=1, N_DIRECTION
            DO I=N_CLOUD_TOP, N_LAYER
              DO L=1, N_PROFILE
                PHASE_FNC_SOLAR(L, I, ID, K)                            &
     &            =PHASE_FNC_SOLAR(L, I, ID, 0)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO
!
!
!
!     Add on the terms representing processes within clouds.
!
!     Loop over the condensed components, calculating their optical
!     properties and then assign them to the arrays for the types of
!     cloud.
!
      DO K=1, N_CONDENSED
!
!       Flags for dealing with components were set in the subroutine
!       set_cloud_pointer. we now determine whether the component is
!       to be included and calculate its optical properties according
!       to the phase of the component. these contributions are added
!       to the arrays for the selected type of cloud.
!
        IF (L_CLOUD_CMP(K)) THEN
!
          IF (I_PHASE_CMP(K) == IP_PHASE_WATER) THEN
!
!           Include scattering by water droplets.
!
! DEPENDS ON: opt_prop_water_cloud
            CALL OPT_PROP_WATER_CLOUD(IERR                              &
     &        , N_PROFILE, N_LAYER, N_CLOUD_TOP                         &
     &        , N_CLOUD_PROFILE, I_CLOUD_PROFILE                        &
     &        , N_ORDER_PHASE, L_RESCALE, N_ORDER_FORWARD               &
     &        , L_HENYEY_GREENSTEIN_PF, L_SOLAR_PHF                     &
     &        , N_ORDER_PHASE_SOLAR, N_DIRECTION, COS_SOL_VIEW          &
     &        , I_CONDENSED_PARAM(K)                                    &
     &        , CONDENSED_PARAM_LIST(1, K)                              &
     &        , CONDENSED_MIX_RATIO(1, ID_CT, K)                        &
     &        , CONDENSED_DIM_CHAR(1, ID_CT, K)                         &
#if defined(OBSERVED)
     &        , P, DENSITY                                              &
     &        , N_OPT_LEVEL_DROP_PRSC, DROP_PRESSURE_PRSC               &
     &        , DROP_ABSORPTION_PRSC, DROP_SCATTERING_PRSC              &
     &        , DROP_PHASE_FNC_PRSC                                     &
#endif
     &        , K_EXT_TOT_CLOUD_COMP, K_EXT_SCAT_CLOUD_COMP             &
     &        , PHASE_FNC_CLOUD_COMP, FORWARD_SCATTER_CLOUD_COMP        &
     &        , FORWARD_SOLAR_CLOUD_COMP, PHASE_FNC_SOLAR_CLOUD_COMP    &
     &        , ND_PROFILE, ND_RADIANCE_PROFILE, ND_LAYER, ID_CT        &
     &        , ND_DIRECTION, ND_PHASE_TERM, ND_MAX_ORDER               &
     &        , ND_CLOUD_PARAMETER                                      &
#if defined(OBSERVED)
     &        , ND_PROFILE_CLOUD_PRSC, ND_OPT_LEVEL_CLOUD_PRSC          &
#endif
     &        )
!
          ELSE IF (I_PHASE_CMP(K) == IP_PHASE_ICE) THEN
!
!           Include scattering by ice crystals.
!
! DEPENDS ON: opt_prop_ice_cloud
            CALL OPT_PROP_ICE_CLOUD(IERR                                &
     &        , N_PROFILE, N_LAYER, N_CLOUD_TOP                         &
     &        , N_CLOUD_PROFILE, I_CLOUD_PROFILE                        &
     &        , N_ORDER_PHASE, L_RESCALE, N_ORDER_FORWARD               &
     &        , L_HENYEY_GREENSTEIN_PF, L_SOLAR_PHF                     &
     &        , N_ORDER_PHASE_SOLAR, N_DIRECTION, COS_SOL_VIEW          &
     &        , I_CONDENSED_PARAM(K)                                    &
     &        , CONDENSED_PARAM_LIST(1, K)                              &
     &        , CONDENSED_MIX_RATIO(1, ID_CT, K)                        &
     &        , CONDENSED_DIM_CHAR(1, ID_CT, K)                         &
#if defined(OBSERVED)
     &        , P, T, DENSITY                                           &
     &        , N_OPT_LEVEL_ICE_PRSC, ICE_PRESSURE_PRSC                 &
     &        , ICE_ABSORPTION_PRSC, ICE_SCATTERING_PRSC                &
     &        , ICE_PHASE_FNC_PRSC                                      &
#endif
#if defined(UM)
     &        , T, DENSITY                                              &
#endif
     &        , K_EXT_TOT_CLOUD_COMP, K_EXT_SCAT_CLOUD_COMP             &
     &        , PHASE_FNC_CLOUD_COMP, FORWARD_SCATTER_CLOUD_COMP        &
     &        , FORWARD_SOLAR_CLOUD_COMP, PHASE_FNC_SOLAR_CLOUD_COMP    &
     &        , ND_PROFILE, ND_RADIANCE_PROFILE, ND_LAYER, ID_CT        &
     &        , ND_DIRECTION                                            &
     &        , ND_PHASE_TERM, ND_MAX_ORDER, ND_CLOUD_PARAMETER         &
#if defined(OBSERVED)
     &        , ND_PROFILE_CLOUD_PRSC, ND_OPT_LEVEL_CLOUD_PRSC          &
#endif
     &        )
!
          ENDIF
!
!
!
!         Increment the arrays of optical properties.
!
!
          DO I=N_CLOUD_TOP, N_LAYER
            DO LL=1, N_CLOUD_PROFILE(I)
              L=I_CLOUD_PROFILE(LL, I)
              K_EXT_TOT(L, I, I_CLOUD_TYPE(K))                          &
     &          =K_EXT_TOT(L, I, I_CLOUD_TYPE(K))                       &
     &          +K_EXT_TOT_CLOUD_COMP(L, I)
              K_EXT_SCAT(L, I, I_CLOUD_TYPE(K))                         &
     &          =K_EXT_SCAT(L, I, I_CLOUD_TYPE(K))                      &
     &          +K_EXT_SCAT_CLOUD_COMP(L, I)
            ENDDO
            DO LS=1, N_ORDER_PHASE
              DO LL=1, N_CLOUD_PROFILE(I)
                L=I_CLOUD_PROFILE(LL, I)
                PHASE_FNC(L, I, LS, I_CLOUD_TYPE(K))                    &
     &            =PHASE_FNC(L, I, LS, I_CLOUD_TYPE(K))                 &
     &            +PHASE_FNC_CLOUD_COMP(L, I, LS)
              ENDDO
            ENDDO
          ENDDO
          IF (L_RESCALE) THEN
            DO I=N_CLOUD_TOP, N_LAYER
              DO LL=1, N_CLOUD_PROFILE(I)
                L=I_CLOUD_PROFILE(LL, I)
                FORWARD_SCATTER(L, I, I_CLOUD_TYPE(K))                  &
     &            =FORWARD_SCATTER(L, I, I_CLOUD_TYPE(K))               &
     &            +FORWARD_SCATTER_CLOUD_COMP(L, I)
              ENDDO
            ENDDO
          ENDIF
          IF (L_SOLAR_PHF) THEN
            DO I=N_CLOUD_TOP, N_LAYER
              DO LL=1, N_CLOUD_PROFILE(I)
                L=I_CLOUD_PROFILE(LL, I)
                FORWARD_SOLAR(L, I, I_CLOUD_TYPE(K))                    &
     &            =FORWARD_SOLAR(L, I, I_CLOUD_TYPE(K))                 &
     &            +FORWARD_SOLAR_CLOUD_COMP(L, I)
              ENDDO
            ENDDO
            DO I=N_CLOUD_TOP, N_LAYER
              DO ID=1, N_DIRECTION
                DO LL=1, N_CLOUD_PROFILE(I)
                  L=I_CLOUD_PROFILE(LL, I)
                  PHASE_FNC_SOLAR(L, I, ID, I_CLOUD_TYPE(K))            &
     &              =PHASE_FNC_SOLAR(L, I, ID, I_CLOUD_TYPE(K))         &
     &              +PHASE_FNC_SOLAR_CLOUD_COMP(L, I, ID)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
!
        ENDIF
!
      ENDDO
!
!
!
!
!     Calculate the final optical properties.
!     The scattering was included in the free total extinction earlier,
!     but we have yet to divide the product of the phase function and
!     the scattering by the mean scattering.
!
      DO I=1, N_CLOUD_TOP-1
!
        N_INDEX=0
        DO L=1, N_PROFILE
          IF (K_EXT_SCAT_CLR(L, I) >  0.0E+00_Real64) THEN
            N_INDEX=N_INDEX+1
            INDEX(N_INDEX)=L
          ENDIF
        ENDDO
!
        DO LS=1, N_ORDER_PHASE
          DO K=1, N_INDEX
            PHASE_FNC_CLR(INDEX(K), I, LS)                              &
     &        =PHASE_FNC_CLR(INDEX(K), I, LS)                           &
     &        /K_EXT_SCAT_CLR(INDEX(K), I)
          ENDDO
        ENDDO
!
        IF (L_RESCALE) THEN
          DO K=1, N_INDEX
            FORWARD_SCATTER_CLR(INDEX(K), I)                            &
     &        =FORWARD_SCATTER_CLR(INDEX(K), I)                         &
     &        /K_EXT_SCAT_CLR(INDEX(K), I)
          ENDDO
        ENDIF
!
        IF (L_SOLAR_PHF) THEN
          DO K=1, N_INDEX
            FORWARD_SOLAR_CLR(INDEX(K), I)                              &
     &        =FORWARD_SOLAR_CLR(INDEX(K), I)                           &
     &        /K_EXT_SCAT_CLR(INDEX(K), I)
          ENDDO
          DO ID=1, N_DIRECTION
            DO K=1, N_INDEX
              PHASE_FNC_SOLAR_CLR(INDEX(K), I, ID)                      &
     &          =PHASE_FNC_SOLAR_CLR(INDEX(K), I, ID)                   &
     &          /K_EXT_SCAT_CLR(INDEX(K), I)
            ENDDO
          ENDDO
        ENDIF
!
      ENDDO
!
      DO I=N_CLOUD_TOP, N_LAYER
!
        N_INDEX=0
        DO L=1, N_PROFILE
          IF (K_EXT_SCAT(L, I, 0) >  0.0E+00_Real64) THEN
            N_INDEX=N_INDEX+1
            INDEX(N_INDEX)=L
          ENDIF
        ENDDO
!
        DO LS=1, N_ORDER_PHASE
          DO K=1, N_INDEX
            PHASE_FNC(INDEX(K), I, LS, 0)                               &
     &        =PHASE_FNC(INDEX(K), I, LS, 0)                            &
     &        /K_EXT_SCAT(INDEX(K), I, 0)
          ENDDO
        ENDDO
!
        IF (L_RESCALE) THEN
          DO K=1, N_INDEX
            FORWARD_SCATTER(INDEX(K), I, 0)                             &
     &        =FORWARD_SCATTER(INDEX(K), I, 0)                          &
     &        /K_EXT_SCAT(INDEX(K), I, 0)
          ENDDO
        ENDIF
!
        IF (L_SOLAR_PHF) THEN
          DO K=1, N_INDEX
            FORWARD_SOLAR(INDEX(K), I, 0)                               &
     &        =FORWARD_SOLAR(INDEX(K), I, 0)                            &
     &        /K_EXT_SCAT(INDEX(K), I, 0)
          ENDDO
          DO ID=1, N_DIRECTION
            DO K=1, N_INDEX
              PHASE_FNC_SOLAR(INDEX(K), I, ID, 0)                       &
     &          =PHASE_FNC_SOLAR(INDEX(K), I, ID, 0)                    &
     &          /K_EXT_SCAT(INDEX(K), I, 0)
            ENDDO
          ENDDO
        ENDIF
!
      ENDDO
!
!
!     Repeat for clouds.
      DO K=1, N_CLOUD_TYPE
        DO I=N_CLOUD_TOP, N_LAYER
          N_INDEX=0
          DO L=1, N_PROFILE
            IF (K_EXT_SCAT(L, I, K) >  0.0E+00_Real64) THEN
              N_INDEX=N_INDEX+1
              INDEX(N_INDEX)=L
            ENDIF
          ENDDO
!
          DO LS=1, N_ORDER_PHASE
            DO J=1, N_INDEX
              PHASE_FNC(INDEX(J), I, LS, K)                             &
     &          =PHASE_FNC(INDEX(J), I, LS, K)                          &
     &          /K_EXT_SCAT(INDEX(J), I, K)
            ENDDO
          ENDDO
          IF (L_RESCALE) THEN
            DO J=1, N_INDEX
              FORWARD_SCATTER(INDEX(J), I, K)                           &
     &          =FORWARD_SCATTER(INDEX(J), I, K)                        &
     &          /K_EXT_SCAT(INDEX(J), I, K)
            ENDDO
          ENDIF
!
          IF (L_SOLAR_PHF) THEN
            DO J=1, N_INDEX
              FORWARD_SOLAR(INDEX(J), I, K)                             &
     &          =FORWARD_SOLAR(INDEX(J), I, K)                          &
     &          /K_EXT_SCAT(INDEX(J), I, K)
            ENDDO
            DO ID=1, N_DIRECTION
              DO J=1, N_INDEX
                PHASE_FNC_SOLAR(INDEX(J), I, ID, K)                     &
     &            =PHASE_FNC_SOLAR(INDEX(J), I, ID, K)                  &
     &            /K_EXT_SCAT(INDEX(J), I, K)
              ENDDO
            ENDDO
          ENDIF
!
        ENDDO
      ENDDO

!
!
!
      RETURN
      END SUBROUTINE GREY_OPT_PROP
#endif
