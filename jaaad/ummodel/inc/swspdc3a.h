! SWSPDC3A contains declarations for reduced sw-spectral file in
! two-stream radiation code.
! (note: SWSPDC3A, SWSPCM3A and SWSARG3A must be consistent)

      ! dimensions for the sw spectrum
      INTEGER ::NPD_TYPE_SW            ! number of types of data
      INTEGER ::NPD_BAND_SW            ! number of spectral bands
      INTEGER ::NPD_EXCLUDE_SW         ! number of excluded bands
      INTEGER ::NPD_SPECIES_SW         ! number of gaseous species
      INTEGER ::NPD_ESFT_TERM_SW       ! number of esft terms
      INTEGER ::NPD_SCALE_FNC_SW       ! number of scaling functions
      INTEGER ::NPD_SCALE_VARIABLE_SW  ! number of scaling variables
      INTEGER ::NPD_SURFACE_SW         ! number of surface types
      INTEGER ::NPD_ALBEDO_PARM_SW     ! number of albedo parameters
      INTEGER ::NPD_CONTINUUM_SW       ! number of continua
      INTEGER ::NPD_DROP_TYPE_SW       ! number of drop types
      INTEGER ::NPD_ICE_TYPE_SW        ! number of ice crystal types
      INTEGER ::NPD_AEROSOL_SPECIES_SW ! number of aerosol species
      INTEGER ::NPD_CLOUD_PARAMETER_SW ! max number of cloud parameters
      INTEGER ::NPD_HUMIDITIES_SW      ! maximum number of humidities
      INTEGER ::NPD_THERMAL_COEFF_SW   ! number of thermal coefficients
      INTEGER ::NPD_AOD_WAVEL_SW       ! number of wavelengths for AOD

! general fields:
      ! flag for types of data present
      LOGICAL :: L_PRESENT_SW(0: NPD_TYPE_SW)

      ! properties of the spectral bands:
      INTEGER :: N_BAND_SW ! number of spectral bands

      ! shorter wavelength limits
      REAL :: WAVE_LENGTH_SHORT_SW(NPD_BAND_SW)

      ! longer wavelength limits
      REAL :: WAVE_LENGTH_LONG_SW(NPD_BAND_SW)

! exclusion of specific bands from parts of the spectrum:

      ! number of excluded bands within each spectral band
      INTEGER :: N_BAND_EXCLUDE_SW(NPD_BAND_SW)

      ! indices of excluded bands
      INTEGER :: INDEX_EXCLUDE_SW(NPD_EXCLUDE_SW, NPD_BAND_SW)

! fields for the solar flux:

      ! fraction of the incident solar flux in each band

      REAL :: SOLAR_FLUX_BAND_SW(NPD_BAND_SW)

! fields for Rayleigh scattering:

      REAL :: RAYLEIGH_COEFFICIENT_SW(NPD_BAND_SW) ! Rayleigh coeffs

! fields for gaseous absorption:

      INTEGER :: N_ABSORB_SW ! number of absorbers

      ! number of absorbers in each band
      INTEGER :: N_BAND_ABSORB_SW(NPD_BAND_SW)

      ! list of absorbers in each band

      INTEGER :: INDEX_ABSORB_SW(NPD_SPECIES_SW, NPD_BAND_SW)

      ! types of each gas in the spectral file

      INTEGER :: TYPE_ABSORB_SW(NPD_SPECIES_SW)

      ! number of esft terms in each band for each gas
      INTEGER :: I_BAND_ESFT_SW(NPD_BAND_SW, NPD_SPECIES_SW)

      ! type of esft scaling
      INTEGER :: I_SCALE_ESFT_SW(NPD_BAND_SW, NPD_SPECIES_SW)

      ! type of scaling function
      INTEGER :: I_SCALE_FNC_SW(NPD_BAND_SW, NPD_SPECIES_SW)

      ! esft exponents

      REAL :: K_ESFT_SW(NPD_ESFT_TERM_SW, NPD_BAND_SW, NPD_SPECIES_SW)

      ! esft weights

      REAL :: W_ESFT_SW(NPD_ESFT_TERM_SW, NPD_BAND_SW, NPD_SPECIES_SW)

      ! scaling parameters for each absorber and term

      REAL :: SCALE_VECTOR_SW(NPD_SCALE_VARIABLE_SW, NPD_ESFT_TERM_SW,  &
     &  NPD_BAND_SW, NPD_SPECIES_SW)

      ! reference pressure for scaling function
      REAL :: P_REFERENCE_SW(NPD_SPECIES_SW, NPD_BAND_SW)

      ! reference temperature for scaling function
      REAL :: T_REFERENCE_SW(NPD_SPECIES_SW, NPD_BAND_SW)

! representation of the planckian:

      INTEGER :: N_DEG_FIT_SW ! degree of thermal polynomial

      ! coefficients in polynomial fit to source function
      REAL :: THERMAL_COEFFICIENT_SW(0: NPD_THERMAL_COEFF_SW-1,         &
     &  NPD_BAND_SW)

      REAL :: T_REF_PLANCK_SW ! planckian reference temperature

! surface properties:

      ! method of specifying properties of surface
      INTEGER :: I_SPEC_SURFACE_SW(NPD_SURFACE_SW)

      ! number of parameters fitting the direct albedo
      INTEGER :: N_DIR_ALBEDO_FIT_SW(NPD_SURFACE_SW)

      LOGICAL :: L_SURFACE_SW(NPD_SURFACE_SW) ! surface types included

      ! surface albedos

      REAL :: SURFACE_ALBEDO_SW(NPD_BAND_SW, NPD_SURFACE_SW)

      ! coefficients for fitting direct albedo

      REAL :: DIRECT_ALBEDO_PARM_SW(0: NPD_ALBEDO_PARM_SW,NPD_BAND_SW,  &
     &  NPD_SURFACE_SW)

      ! surface emissivities
      REAL :: EMISSIVITY_GROUND_SW(NPD_BAND_SW, NPD_SURFACE_SW)

! fields for continua:

      ! number of continua in each band
      INTEGER :: N_BAND_CONTINUUM_SW(NPD_BAND_SW)

      ! list of continua in each band
      INTEGER :: INDEX_CONTINUUM_SW(NPD_BAND_SW, NPD_CONTINUUM_SW)

      INTEGER :: INDEX_WATER_SW ! index of water vapour

      ! type of scaling function for continuum
      INTEGER :: I_SCALE_FNC_CONT_SW(NPD_BAND_SW, NPD_CONTINUUM_SW)

      ! grey extinction coefficients for continuum
      REAL :: K_CONTINUUM_SW(NPD_BAND_SW, NPD_CONTINUUM_SW)

      ! scaling parameters for continuum
      REAL :: SCALE_CONTINUUM_SW(NPD_SCALE_VARIABLE_SW,NPD_BAND_SW,     &
     &  NPD_CONTINUUM_SW)

      ! reference pressure for scaling of continuum
      REAL :: P_REF_CONTINUUM_SW(NPD_CONTINUUM_SW, NPD_BAND_SW)

      ! reference temperature for scaling of continuum
      REAL :: T_REF_CONTINUUM_SW(NPD_CONTINUUM_SW, NPD_BAND_SW)

! fields for water droplets:

      ! parametrization type of droplets
      INTEGER :: I_DROP_PARAMETRIZATION_SW(NPD_DROP_TYPE_SW)

      ! types of droplet present
      LOGICAL :: L_DROP_TYPE_SW(NPD_DROP_TYPE_SW)

      ! parameters used to fit optical properties of clouds

      REAL :: DROP_PARAMETER_LIST_SW(NPD_CLOUD_PARAMETER_SW,NPD_BAND_SW,&
     &  NPD_DROP_TYPE_SW)

      ! minimum dimension permissible in the parametrization
      REAL :: DROP_PARM_MIN_DIM_SW(NPD_DROP_TYPE_SW)

      ! maximum dimension permissible in the parametrization
      REAL :: DROP_PARM_MAX_DIM_SW(NPD_DROP_TYPE_SW)

!     fields for aerosols:

      INTEGER :: N_AEROSOL_SW ! number of species of aerosol

      ! types of aerosols
      INTEGER :: TYPE_AEROSOL_SW(NPD_AEROSOL_SPECIES_SW)

      ! parametrization of aerosols
      INTEGER :: I_AEROSOL_PARAMETRIZATION_SW(NPD_AEROSOL_SPECIES_SW)

      ! numbers of humidities
      INTEGER :: NHUMIDITY_SW(NPD_AEROSOL_SPECIES_SW)

      ! aerosol species included
      LOGICAL :: L_AEROSOL_SPECIES_SW(NPD_AEROSOL_SPECIES_SW)

      ! absorption by aerosols
      REAL :: AEROSOL_ABSORPTION_SW(NPD_HUMIDITIES_SW,                  &
     &  NPD_AEROSOL_SPECIES_SW, NPD_BAND_SW)

      ! scattering by aerosols
      REAL :: AEROSOL_SCATTERING_SW(NPD_HUMIDITIES_SW,                  &
     &  NPD_AEROSOL_SPECIES_SW, NPD_BAND_SW)

      ! asymmetry of aerosols
      REAL :: AEROSOL_ASYMMETRY_SW(NPD_HUMIDITIES_SW,                   &
     &  NPD_AEROSOL_SPECIES_SW, NPD_BAND_SW)

      ! humidities for components
      REAL :: HUMIDITIES_SW(NPD_HUMIDITIES_SW, NPD_AEROSOL_SPECIES_SW)

! fields for ice crystals:

      ! types of parametrization of ice crystals

      INTEGER :: I_ICE_PARAMETRIZATION_SW(NPD_ICE_TYPE_SW)

      ! types of ice crystal present
      LOGICAL :: L_ICE_TYPE_SW(NPD_ICE_TYPE_SW)

      ! parameters used to fit single scattering of ice crystals
      REAL :: ICE_PARAMETER_LIST_SW(NPD_CLOUD_PARAMETER_SW, NPD_BAND_SW,&
     &  NPD_ICE_TYPE_SW)

      ! minimum dimension permissible in the parametrization
      REAL :: ICE_PARM_MIN_DIM_SW(NPD_ICE_TYPE_SW)

      ! maximum dimension permissible in the parametrization
      REAL :: ICE_PARM_MAX_DIM_SW(NPD_ICE_TYPE_SW)

! fields for doppler broadening:

      ! flag for doppler broadening for each species
      LOGICAL :: L_DOPPLER_PRESENT_SW(NPD_SPECIES_SW)

      !  offset to pressure to represent doppler broadening
      REAL :: DOPPLER_CORRECTION_SW(NPD_SPECIES_SW)

! fields for aerosol optical depth

      ! number of wavelengths
      INTEGER :: N_AOD_WAVEL_SW

      ! monochromatic specific absorption coefficient
      REAL :: AOD_ABSORPTION_SW(NPD_HUMIDITIES_SW,                      &
     &          NPD_AEROSOL_SPECIES_SW, NPD_AOD_WAVEL_SW)

      ! monochromatic specific scattering coefficient
      REAL :: AOD_SCATTERING_SW(NPD_HUMIDITIES_SW,                      &
     &          NPD_AEROSOL_SPECIES_SW, NPD_AOD_WAVEL_SW)

      ! relationship between aerosol component and type
      INTEGER :: I_AOD_TYPE_SW(NPD_AEROSOL_SPECIES_SW)

! SWSPDC3A end
