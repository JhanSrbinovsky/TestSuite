! SPDEC3A contains declarations for spectral file in two-stream
! radiation code.

! general fields:

      LOGICAL:: L_PRESENT(0:NPD_TYPE) ! flag for types of data present

! properties of the spectral bands:

      INTEGER :: N_BAND ! number of spectral bands

      REAL ::  WAVE_LENGTH_SHORT(NPD_BAND) ! shorter wavelength limits
      REAL ::  WAVE_LENGTH_LONG(NPD_BAND)  ! longer wavelength limits

      ! exclusion of specific bands from parts of the spectrum:

      ! number of excluded bands within each spectral band
      INTEGER :: N_BAND_EXCLUDE(NPD_BAND)

      ! indices of excluded bands
      INTEGER :: INDEX_EXCLUDE(NPD_EXCLUDE,NPD_BAND)

! fields for the solar flux:

      ! fraction of the incident solar flux in each band
      REAL :: SOLAR_FLUX_BAND(NPD_BAND)

      ! fields for Rayleigh scattering:

      REAL ::  RAYLEIGH_COEFFICIENT(NPD_BAND) ! Rayleigh coefficients

      ! fields for gaseous absorption:

      INTEGER :: N_ABSORB ! number of absorbers

      ! number of absorbers in each band
      INTEGER :: N_BAND_ABSORB(NPD_BAND)

      ! list of absorbers in each band
      INTEGER ::  INDEX_ABSORB(NPD_SPECIES, NPD_BAND)

      ! types of each gas in the spectral file
      INTEGER ::  TYPE_ABSORB(NPD_SPECIES)

      ! number of esft terms in band for each gas
      INTEGER ::  I_BAND_ESFT(NPD_BAND, NPD_SPECIES)

      ! type of esft scaling
      INTEGER ::  I_SCALE_ESFT(NPD_BAND, NPD_SPECIES)

      ! type of scaling function
      INTEGER :: I_SCALE_FNC(NPD_BAND, NPD_SPECIES)

      ! esft exponents
      REAL :: K_ESFT(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)

      ! esft weights
      REAL :: W_ESFT(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)

      ! scaling parameters for each absorber and term
      REAL :: SCALE_VECTOR(NPD_SCALE_VARIABLE, NPD_ESFT_TERM, NPD_BAND, &
     &  NPD_SPECIES)

      ! reference pressure for scaling function
      REAL :: P_REFERENCE(NPD_SPECIES, NPD_BAND)

      ! reference temperature for scaling function
      REAL :: T_REFERENCE(NPD_SPECIES, NPD_BAND)

      ! representation of the Planckian:

      INTEGER :: N_DEG_FIT ! degree of thermal polynomial

      ! coefficients in polynomial fit to source function

      REAL :: THERMAL_COEFFICIENT(0: NPD_THERMAL_COEFF-1, NPD_BAND)

      REAL :: T_REF_PLANCK ! planckian reference temperature

! surface properties:

      ! method of specifying properties of surface

      INTEGER :: I_SPEC_SURFACE(NPD_SURFACE)

      ! number of parameters fitting the direct albedo
      INTEGER :: N_DIR_ALBEDO_FIT(NPD_SURFACE)

      LOGICAL :: L_SURFACE(NPD_SURFACE) ! surface types included

      REAL :: SURFACE_ALBEDO(NPD_BAND, NPD_SURFACE) ! surface albedos

! coefficients for fitting direct albedo

      REAL :: DIRECT_ALBEDO_PARM(0:NPD_ALBEDO_PARM,NPD_BAND,NPD_SURFACE)

      ! surface emissivities
      REAL ::   EMISSIVITY_GROUND(NPD_BAND, NPD_SURFACE)


! fields for continua:

      ! number of continua in each band
      INTEGER :: N_BAND_CONTINUUM(NPD_BAND)

      ! list of continua continuua in each band
      INTEGER :: INDEX_CONTINUUM(NPD_BAND, NPD_CONTINUUM)

      INTEGER :: INDEX_WATER ! index of water vapour

      ! type of scaling function for continuum
      INTEGER :: I_SCALE_FNC_CONT(NPD_BAND, NPD_CONTINUUM)

      ! grey extinction coefficients for continuum
      REAL :: K_CONTINUUM(NPD_BAND, NPD_CONTINUUM)

      ! scaling parameters for continuum
      REAL :: SCALE_CONTINUUM(NPD_SCALE_VARIABLE,NPD_BAND,NPD_CONTINUUM)

      ! reference pressure for scaling of continuum
      REAL :: P_REF_CONTINUUM(NPD_CONTINUUM, NPD_BAND)

      ! reference temperature for scaling of continuum
      REAL :: T_REF_CONTINUUM(NPD_CONTINUUM, NPD_BAND)

! fields for water droplets:

      ! parametrization type of droplets

      INTEGER :: I_DROP_PARAMETRIZATION(NPD_DROP_TYPE)
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
      INTEGER :: N_DROP_PHF_TERM(NPD_DROP_TYPE)
!       Number of terms in the phase function for droplets
#endif

      ! types of droplet present

      LOGICAL :: L_DROP_TYPE(NPD_DROP_TYPE)

      ! parameters used to fit optical properties of clouds
      REAL :: DROP_PARAMETER_LIST(NPD_CLOUD_PARAMETER,NPD_BAND,         &
     &  NPD_DROP_TYPE)

      ! minimum dimension permissible in the parametrization
      REAL :: DROP_PARM_MIN_DIM(NPD_DROP_TYPE)

      ! maximum dimension permissible in the parametrization
      REAL :: DROP_PARM_MAX_DIM(NPD_DROP_TYPE)

! fields for aerosols:

      INTEGER :: N_AEROSOL ! number of species of aerosol

      ! types of aerosols
      INTEGER :: TYPE_AEROSOL(NPD_AEROSOL_SPECIES)

      ! parametrization of aerosols
      INTEGER :: I_AEROSOL_PARAMETRIZATION(NPD_AEROSOL_SPECIES)
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
      INTEGER :: N_AEROSOL_PHF_TERM(NPD_AEROSOL_SPECIES)
!       Number of terms in the phase function
#endif

      ! numbers of humidities
      INTEGER :: NHUMIDITY(NPD_AEROSOL_SPECIES)

      ! aerosol species included

      LOGICAL :: L_AEROSOL_SPECIES(NPD_AEROSOL_SPECIES)

      ! absorption by aerosols
      REAL :: AEROSOL_ABSORPTION(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES,   &
     &  NPD_BAND)

      ! scattering by aerosols
      REAL :: AEROSOL_SCATTERING(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES,   &
     &  NPD_BAND)

#if defined(A01_3A) || defined(A02_3A)
      REAL :: AEROSOL_ASYMMETRY(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES     &
     &     ,  NPD_BAND)
!       Asymmetry of aerosols
#endif
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
      REAL :: AEROSOL_PHASE_FNC(NPD_HUMIDITIES, NPD_PHASE_TERM          &
     &     ,  NPD_AEROSOL_SPECIES, NPD_BAND)
      REAL :: AEROSOL_ASYMMETRY(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES     &
     &     ,  NPD_BAND)
!       Asymmetry of aerosols
#endif

      ! humidities for components
      REAL :: HUMIDITIES(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES)

! fields for ice crystals:

      ! types of parametrization of ice crystals
      INTEGER :: I_ICE_PARAMETRIZATION(NPD_ICE_TYPE)
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
      INTEGER :: N_ICE_PHF_TERM(NPD_ICE_TYPE)
!       number of terms in the phase function for ice crystals
#endif

      ! types of ice crystal present
      LOGICAL :: L_ICE_TYPE(NPD_ICE_TYPE)

      ! parameters used to fit single scattering of ice crystals
      REAL :: ICE_PARAMETER_LIST(NPD_CLOUD_PARAMETER,NPD_BAND,          &
     &   NPD_ICE_TYPE)

      ! minimum dimension permissible in the parametrization
      REAL :: ICE_PARM_MIN_DIM(NPD_ICE_TYPE)

      ! maximum dimension permissible in the parametrization
      REAL :: ICE_PARM_MAX_DIM(NPD_ICE_TYPE)

! fields for doppler broadening:

      ! flag for doppler broadening for each species
      LOGICAL :: L_DOPPLER_PRESENT(NPD_SPECIES)

      ! doppler correction terms
      REAL :: DOPPLER_CORRECTION(NPD_SPECIES)

! fields for aerosol optical depth

      ! number of wavelengths at which the AOD is computed
      INTEGER :: N_AOD_WAVEL

      ! Monochromatic specific absorption and scattering
      ! coefficients for aerosols
      REAL :: AOD_ABSORPTION(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES,       &
     &  NPD_AOD_WAVEL)
      REAL :: AOD_SCATTERING(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES,       &
     &  NPD_AOD_WAVEL)

      ! Aerosol type for each aerosol component
      INTEGER :: I_AOD_TYPE(NPD_AEROSOL_SPECIES)

! SPDEC3A end
