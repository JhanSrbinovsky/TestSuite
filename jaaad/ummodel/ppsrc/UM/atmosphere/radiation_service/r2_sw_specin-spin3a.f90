

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read a shortwave spectral namelist.
!
! Purpose:
!   To read a shortwave namelist into a spectral array.
!
! Method:
!   The spectrum is read into the dynamically allocated array
!   and then reduced to a more manageable size.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             14-05-96                Set lower limits
!                                               for reduced dimensions
!                                               to ensure that they
!                                               may never be 0.
!                                               (J. M. Edwards)
!       4.4             02-09-97                Aerosol flags passed
!                                               in to the code to
!                                               enable only those
!                                               required to be
!                                               selected. Spectral
!                                               data are now longer
!                                               compressed into a
!                                               single array.
!                                               Actual IOS code put
!                                               into CMESSAGE.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Coding to allow
!                                               selection of gases
!                                               from the spectral
!                                               file.
!                                               (J. M. Edwards)
!
!       4.5        April 1998   Allow soot spectral data to be read.
!                                                     Luke Robinson.
!       5.2             15-11-00                Coding to allow
!                                               sea-salt aerosol
!                                               to be read from
!                                               the spectral file
!                                               if required.
!                                               (A. Jones)
!       5.3             04-04-01                Allow aerosol data to be
!                                               read in if required by
!                                               mesoscale model.
!                                               (S. Cusack)
!       5.5             05-02-03                Allow biomass aerosol
!                                               spectral data to be
!                                               read in, if required.
!                                               (P. Davison)
!       5.5             21-02-03                Allow mineral dust
!                                               spectral data to be
!                                               read in, if required.
!                                               (S Woodward)
!       6.2             16-11-05                Argument list of
!                                               R2_COMPRESS_SPECTRUM
!                                               modified to match
!                                               changes in the LW.
!                                               (N Bellouin)
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SW_SPECIN(IERR, CMESSAGE                            &
     &   , L_O2                                                         &
     &   , L_CLIMAT_AEROSOL                                             &
     &   , L_USE_DUST, L_USE_ARCLDUST                                   &
     &   , L_USE_SULPC_DIRECT, L_USE_ARCLSULP                           &
     &   , L_USE_SOOT_DIRECT, L_USE_ARCLBLCK                            &
     &   , L_USE_BMASS_DIRECT, L_USE_ARCLBIOM                           &
     &   , L_USE_SEASALT_DIRECT, L_USE_ARCLSSLT                         &
     &   , L_USE_OCFF_DIRECT, L_USE_ARCLOCFF                            &
     &   , L_USE_BIOGENIC, L_USE_ARCLDLTA                               &
     &   , L_MURK_RAD                                                   &
     &   )
!
!
      IMPLICIT NONE
!
!
!     ------------------------------------------------------------------
!     MODULE SETTING MAXIMUM DIMENSIONS OF ARRAYS IN THE RADIATION CODE.
!   4.5   Aug 1998     Increment by 2 the no. of aerosol species
!                      affecting the radiation.    Luke Robinson
!
!
!   5.2   Dec 2000     Increment by 2 the no. of aerosol species
!                      affecting the radiation.    Andy Jones
!   5.3   Oct 2001     Allow six more types of aerosol in preparation
!                      for runs with dust.      (J. M. Edwards)
!   5.5   Feb 2003     Increase aerosol species to include dust
!                      and biomass    (S. Woodward)
!   5.5   Feb 2003     Increase the number of types of ice crystal.
!                                                 John Edwards
!   6.2   Nov 2005     Increase number of types of data.
!                      Add number of wavelengths for aerosol optical
!                      depth.
!                                                 N. Bellouin
!   6.2   Jan 2006     Increase values of the cloud parameters and
!                      of the number of terms in the phase function
!                      under versions 3C and 3Z of the radiation
!                      code.                         (J.-C. Thelen)
!
      INTEGER,PARAMETER :: NPD_TYPE      = 16 ! number of types of data
      INTEGER,PARAMETER :: NPD_BAND      = 20 ! number of spectral bands
      INTEGER,PARAMETER :: NPD_EXCLUDE   = 2  ! nubmer of excluded bands
      INTEGER,PARAMETER :: NPD_SPECIES   = 11 ! number of gaseous species
      INTEGER,PARAMETER :: NPD_ESFT_TERM = 16 ! number of esft terms
      INTEGER,PARAMETER :: NPD_SCALE_FNC = 4  ! number of scaling funcs

      ! number of scaling variables
      INTEGER,PARAMETER :: NPD_SCALE_VARIABLE=4

      INTEGER,PARAMETER :: NPD_SURFACE     = 1 ! no of surface types
      INTEGER,PARAMETER :: NPD_ALBEDO_PARM = 4 ! no of albedo parameters
      INTEGER,PARAMETER :: NPD_CONTINUUM   = 2 ! no of continua
      INTEGER,PARAMETER :: NPD_DROP_TYPE   = 6 ! no of drop types
      INTEGER,PARAMETER :: NPD_ICE_TYPE    = 8 ! no of ice crystal types


      ! max no of cloud parameters
      INTEGER,PARAMETER :: NPD_CLOUD_PARAMETER=30

      INTEGER, PARAMETER :: NPD_AEROSOL_SPECIES=24

      ! maximum no of humidities
      INTEGER,PARAMETER :: NPD_HUMIDITIES=21

      ! no of thermal coefficients
      INTEGER,PARAMETER :: NPD_THERMAL_COEFF=9

      ! number of wavelengths for aerosol optical depth
      INTEGER,PARAMETER :: NPD_AOD_WAVEL=6
! MXSIZE3A end
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET ERROR FLAGS IN THE RADIATION CODE.
!
      INTEGER                                                           &
     &     I_NORMAL                                                     &
!             ERROR FREE CONDITION
     &   , I_ERR_FATAL                                                  &
!             FATAL ERROR: IMMEDIATE RETURN
     &   , I_ABORT_CALCULATION                                          &
!             CALCULATION ABORTED
     &   , I_MISSING_DATA                                               &
!             MISSING DATA ERROR: CONDITIONAL
     &   , I_ERR_IO                                                     &
!             I/O ERROR
     &   , I_ERR_RANGE                                                  &
!             INTERPOLATION RANGE ERROR
     &   , I_ERR_EXIST
!             EXISTENCE ERROR
!
      PARAMETER(                                                        &
     &     I_NORMAL=0                                                   &
     &   , I_ERR_FATAL=1                                                &
     &   , I_ABORT_CALCULATION=2                                        &
     &   , I_MISSING_DATA=3                                             &
     &   , I_ERR_IO=4                                                   &
     &   , I_ERR_RANGE=5                                                &
     &   , I_ERR_EXIST=6                                                &
     &   )
!
!     ------------------------------------------------------------------
! STDIO3A defines unit numbers for standard i/o in two-stream radiation
! code.
      INTEGER,PARAMETER:: IU_STDIN=5
      INTEGER,PARAMETER:: IU_STDOUT=6
      INTEGER,PARAMETER:: IU_ERR=6
! STDIO3A end
!
!
!     DUMMY ARGUMENTS
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_O2                                                         &
!             ABSORPTION BY OXYGEN IS TO BE INCLUDED.
     &   , L_CLIMAT_AEROSOL                                             &
!             CLIMATOLOGICAL AEROSOLS ARE TO BE INCLUDED
     &   , L_USE_DUST                                                   &
                      ! direct effects of mineral dust aerosol
     &   , L_USE_ARCLDUST                                               &
                 ! direct effect of mineral dust aerosol from NWP clim
     &   , L_USE_SULPC_DIRECT                                           &
!             THE DIRECT EFFECTS OF SULPHATE AEROSOLS ARE
!             TO BE INCLUDED
     &   , L_USE_ARCLSULP                                               &
                 ! direct effect of sulphate aerosol from NWP clim
     &   , L_USE_SOOT_DIRECT                                            &
!             USE THE DIRECT RAD EFFECTS OF SOOT IN THE SW
     &   , L_USE_ARCLBLCK                                               &
                 ! direct effect of black-carbon aerosol from NWP clim
     &   , L_USE_BMASS_DIRECT                                           &
!             USE THE DIRECT RAD EFFECTS OF BIOMASS SMOKE IN THE SW
     &   , L_USE_ARCLBIOM                                               &
                 ! direct effect of biomass aerosol from NWP clim
     &   , L_USE_SEASALT_DIRECT                                         &
!             USE THE DIRECT RAD EFFECTS OF SEASALT IN THE SW
     &   , L_USE_ARCLSSLT                                               &
                 ! direct effect of sea-salt aerosol from NWP clim
     &   , L_USE_BIOGENIC                                               &
!             USE THE BIOGENIC AEROSOL DIRECT EFFECT IN THE SW
     &   , L_USE_OCFF_DIRECT                                            &
                 ! direct effect of fossil-fuel org.carb. in the SW
     &   , L_USE_ARCLOCFF                                               &
                 ! direct effect of fossil-fuel org.carb. from NWP clim
     &   , L_USE_ARCLDLTA                                               &
                 ! direct effect of delta aerosol from NWP climatology
     &   , L_MURK_RAD
!             MESOSCALE MODEL AEROSOLS ARE TO BE INCLUDED
!
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      CHARACTER*80                                                      &
                        !, INTENT(OUT)
     &     CMESSAGE
!
!
!
!     LOCAL VARIABLES.
!
!
!     RADIATIVE VARIABLES FOR REDUCING THE SPECTRUM
!
! AERPRM3A aerosol parametrizations for two-stream radiation code.

      INTEGER,PARAMETER:: IP_AEROSOL_PARAM_DRY=1
      INTEGER,PARAMETER:: IP_AEROSOL_PARAM_MOIST=2
      INTEGER,PARAMETER:: IP_AEROSOL_UNPARAMETRIZED=3 ! Observational

! AERPRM3A end
! AERCMP3A start
!     ------------------------------------------------------------------
!     MODULE TO SET INDICES OF AEROSOL COMPONENTS.
!   4.5   Aug 1998     Set indices for two soot aerosol species.
!                                                 Luke Robinson
!
!
!   5.2   Dec 2000     Set indices for two sea-salt aerosol modes.
!                                                 Andy Jones
!
!   5.4   May 2002     Correct value of NPD_AEROSOL_COMPONENT
!                      and add warning comment.
!                                                 Andy Jones
!   5.5   Feb 2003     Set indices for two modes of biomass
!                      aerosol.                   P Davison
!   5.5   Feb 2003     Set indices for mineral dust aerosol.
!                                        S Woodward
!
      ! maximum number of aerosol components
      ! N.B: this must be at least as large as the
      !      largest value in the list below
      INTEGER,PARAMETER:: NPD_AEROSOL_COMPONENT=28

      INTEGER,PARAMETER:: IP_WATER_SOLUBLE=1
      INTEGER,PARAMETER:: IP_DUST_LIKE=2
      INTEGER,PARAMETER:: IP_OCEANIC=3
      INTEGER,PARAMETER:: IP_SOOT=4
      INTEGER,PARAMETER:: IP_ASH=5
      INTEGER,PARAMETER:: IP_SULPHURIC=6
      INTEGER,PARAMETER:: IP_ACCUM_SULPHATE=10
      INTEGER,PARAMETER:: IP_AITKEN_SULPHATE=11
      INTEGER,PARAMETER:: IP_FRESH_SOOT=12
      INTEGER,PARAMETER:: IP_AGED_SOOT=13
      INTEGER,PARAMETER:: IP_SEASALT_FILM=15
      INTEGER,PARAMETER:: IP_SEASALT_JET=16
      INTEGER,PARAMETER:: IP_DUST_1=17
      INTEGER,PARAMETER:: IP_DUST_2=18
      INTEGER,PARAMETER:: IP_DUST_3=19
      INTEGER,PARAMETER:: IP_DUST_4=20
      INTEGER,PARAMETER:: IP_DUST_5=21
      INTEGER,PARAMETER:: IP_DUST_6=22
      INTEGER,PARAMETER:: IP_BIOMASS_1=23
      INTEGER,PARAMETER:: IP_BIOMASS_2=24
      INTEGER,PARAMETER:: IP_BIOGENIC=25
      INTEGER,PARAMETER:: IP_OCFF_FRESH=26
      INTEGER,PARAMETER:: IP_OCFF_AGED=27
      INTEGER,PARAMETER:: IP_DELTA=28
! AERCMP3A end
! GASID3A defines indexing numbers of gaseous absorbing species for
! two-stream radiation code.
! the numbering 1-12 corresponds to lowtran 7.

      INTEGER,PARAMETER:: NPD_GASES=19 ! Number of indexed gases

      INTEGER,PARAMETER:: IP_H2O=1
      INTEGER,PARAMETER:: IP_CO2=2
      INTEGER,PARAMETER:: IP_O3=3
      INTEGER,PARAMETER:: IP_N2O=4
      INTEGER,PARAMETER:: IP_CO=5
      INTEGER,PARAMETER:: IP_CH4=6
      INTEGER,PARAMETER:: IP_O2=7
      INTEGER,PARAMETER:: IP_NO=8
      INTEGER,PARAMETER:: IP_SO2=9
      INTEGER,PARAMETER:: IP_NO2=10
      INTEGER,PARAMETER:: IP_NH3=11
      INTEGER,PARAMETER:: IP_HNO3=12
      INTEGER,PARAMETER:: IP_N2=13
      INTEGER,PARAMETER:: IP_CFC11=14
      INTEGER,PARAMETER:: IP_CFC12=15
      INTEGER,PARAMETER:: IP_CFC113=16
      INTEGER,PARAMETER:: IP_HCFC22=17
      INTEGER,PARAMETER:: IP_HFC125=18
      INTEGER,PARAMETER:: IP_HFC134A=19

! GASID3A end
!
      CHARACTER*256                                                     &
     &     SW_SPECTRAL_FILE
!             NAME OF FILE CONTAINING THE SPECTRAL DATA
      INTEGER                                                           &
     &     IERR_GET_FILE                                                &
!             ERROR FLAG RETURNED BY GET_FILE (NOT NECESSARILY
!             CONSISTENT WITH THE FLAGS IN ERROR3A).
     &   , IOS
!             STATUS OF I/O
!
      LOGICAL                                                           &
     &     L_RETAIN_ABSORB(NPD_SPECIES)                                 &
!             FLAG SET TO .TRUE. IF THE ABSORBER IS TO BE RETAINED
     &   , L_GAS_INCLUDED(NPD_GASES)
!             LOGICAL TO TEST FOR ACTUAL GASES INCLUDED
      INTEGER                                                           &
     &     N_ABSORB_RETAIN                                              &
!             NUMBER OF ABSORBERS TO RETAIN
     &   , INDEX_ABSORB_RETAIN(NPD_SPECIES)                             &
!             INDICES OF ABSORBERS TO BE RETAINED
     &   , COMPRESSED_INDEX(NPD_SPECIES)                                &
!             MAPPING FROM ORIGINAL TO COMPRESSED INDICES OF ABSORBERS
     &   , N_AEROSOL_RETAIN                                             &
!             NUMBER OF AEROSOLS IN THE SPECTRAL FILE TO BE RETAINED
!             FOR THE RADIATIVE CALCULATION
     &   , INDEX_AEROSOL_RETAIN(NPD_AEROSOL_SPECIES)                    &
!             INDEXING NUMBERS OF THE RETAINED AEROSOLS
     &   , N_AEROSOL_FOUND
!             NUMBER OF AEROSOLS FOR THE CURRENT GROUP OF PROCESSES
!             FOUND IN THE SPECTRAL FILE
!
!
!
!     DECLARE THE ELEMENTS OF THE INITIAL SPECTRUM FOR DYNAMIC
!     ALLOCATION AND SET UP AN APPROPRIATE NAMELIST.
!
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

      REAL :: AEROSOL_ASYMMETRY(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES     &
     &     ,  NPD_BAND)
!       Asymmetry of aerosols

      ! humidities for components
      REAL :: HUMIDITIES(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES)

! fields for ice crystals:

      ! types of parametrization of ice crystals
      INTEGER :: I_ICE_PARAMETRIZATION(NPD_ICE_TYPE)

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
! SWSP3A declares elements of a spectral file as a namelist.
!
! 6.2  24/10/05 Functionality for radiative forcing, timestepping
!               and radiances under versions 3C and 3Z of radiation
!               code added                 (J.-C. Thelen)
!

      NAMELIST/R2SWSP/                                                  &
!                       Blocks Present
     &  L_PRESENT,                                                      &
!                       Block 0
     &  N_BAND, N_ABSORB, N_AEROSOL, TYPE_ABSORB, TYPE_AEROSOL,         &
!                       Block 1
     &  WAVE_LENGTH_SHORT, WAVE_LENGTH_LONG,                            &
!                       Block 2
     &  SOLAR_FLUX_BAND,                                                &
!                       Block 3
     &  RAYLEIGH_COEFFICIENT,                                           &
!                       Block 4
     &  N_BAND_ABSORB, INDEX_ABSORB,                                    &
!                       Block 5
     &  I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC,                         &
     &  P_REFERENCE, T_REFERENCE, K_ESFT, W_ESFT, SCALE_VECTOR,         &
!                       Block 6
     &  N_DEG_FIT, T_REF_PLANCK, THERMAL_COEFFICIENT,                   &
!                       Block 7
     &  I_SPEC_SURFACE, N_DIR_ALBEDO_FIT, L_SURFACE,                    &
     &  SURFACE_ALBEDO, DIRECT_ALBEDO_PARM, EMISSIVITY_GROUND,          &
!                       Block 8
     &  N_BAND_CONTINUUM, INDEX_CONTINUUM, INDEX_WATER,                 &
!                       Block 9
     &  I_SCALE_FNC_CONT, P_REF_CONTINUUM, T_REF_CONTINUUM,             &
     &  K_CONTINUUM, SCALE_CONTINUUM,                                   &
!                       Block 10
     &  I_DROP_PARAMETRIZATION, L_DROP_TYPE, DROP_PARAMETER_LIST,       &
     &  DROP_PARM_MIN_DIM, DROP_PARM_MAX_DIM,                           &
!                       Block 11
     &  I_AEROSOL_PARAMETRIZATION, NHUMIDITY, L_AEROSOL_SPECIES,        &
     &  AEROSOL_ABSORPTION, AEROSOL_SCATTERING, AEROSOL_ASYMMETRY,      &
     &  HUMIDITIES,                                                     &
!                       Block 12
     &  I_ICE_PARAMETRIZATION, L_ICE_TYPE, ICE_PARAMETER_LIST,          &
     &  ICE_PARM_MIN_DIM, ICE_PARM_MAX_DIM,                             &
!                       Block 13
     &  L_DOPPLER_PRESENT, DOPPLER_CORRECTION,                          &
!                       Block 14
     &  N_BAND_EXCLUDE, INDEX_EXCLUDE,                                  &
!                       Block 15
     &  N_AOD_WAVEL, AOD_ABSORPTION, AOD_SCATTERING, I_AOD_TYPE

! SWSP3A end
!
!
!     DECLARE THE REDUCED SW SPECTRAL FILE AND ITS HOLDING COMMON BLOCK.
!
! SWSPDL3A contains declarations for reduced sw-spectral file in
! two-stream radiation code.

!     note: SWSPDC3A, SWSPCM3A and SWSARG3A must be consistent
!     note: since the arrays here will be passed in a common block
!     their sizes must be fixed, even though variable sizes are used
!     lower in the code. they are accordingly defined as 1-dimensional
!     arrays with fixed maximum sizes at this level.

! dimensions for the SW spectrum

      INTEGER::  NPD_TYPE_SW            ! number of types of data
      INTEGER::  NPD_BAND_SW            ! number of spectral bands
      INTEGER::  NPD_EXCLUDE_SW         ! number of excluded bands
      INTEGER::  NPD_SPECIES_SW         ! number of gaseous species
      INTEGER::  NPD_ESFT_TERM_SW       ! number of esft terms
      INTEGER::  NPD_SCALE_FNC_SW       ! number of scaling functions
      INTEGER::  NPD_SCALE_VARIABLE_SW  ! number of scaling variables
      INTEGER::  NPD_SURFACE_SW         ! number of surface types
      INTEGER::  NPD_ALBEDO_PARM_SW     ! number of albedo parameters
      INTEGER::  NPD_CONTINUUM_SW       ! number of continua
      INTEGER::  NPD_DROP_TYPE_SW       ! number of drop types
      INTEGER::  NPD_ICE_TYPE_SW        ! number of ice crystal types
      INTEGER::  NPD_AEROSOL_SPECIES_SW ! number of aerosol species
      INTEGER::  NPD_CLOUD_PARAMETER_SW ! max number of cloud parameters
      INTEGER::  NPD_HUMIDITIES_SW      ! maximum number of humidities
      INTEGER::  NPD_THERMAL_COEFF_SW   ! number of thermal coefficients
      INTEGER::  NPD_AOD_WAVEL_SW       ! number of wavelengths for AOD

! general fields:

      ! flag for types of data present

      LOGICAL :: L_PRESENT_SW(0: NPD_TYPE)

! properties of the spectral bands:

      INTEGER :: N_BAND_SW  ! number of spectral bands

      REAL :: WAVE_LENGTH_SHORT_SW(NPD_BAND) ! shorter wavelength limits
      REAL :: WAVE_LENGTH_LONG_SW(NPD_BAND)  ! longer wavelength limits

! exclusion of specific bands from parts of the spectrum:

      ! number of excluded bands within each spectral band
      INTEGER :: N_BAND_EXCLUDE_SW(NPD_BAND)

      ! indices of excluded bands

      INTEGER :: INDEX_EXCLUDE_SW(NPD_EXCLUDE, NPD_BAND)

! fields for the solar flux:

      ! fraction of the incident solar flux in each band
      REAL :: SOLAR_FLUX_BAND_SW(NPD_BAND)

! fields for Rayleigh scattering:

      REAL :: RAYLEIGH_COEFFICIENT_SW(NPD_BAND) ! Rayleigh coefficients

! fields for gaseous absorption:

      INTEGER :: N_ABSORB_SW ! number of absorbers

      ! number of absorbers in each band
      INTEGER :: N_BAND_ABSORB_SW(NPD_BAND)

      ! list of absorbers in each band
      INTEGER :: INDEX_ABSORB_SW(NPD_SPECIES, NPD_BAND)

      ! types of each gas in the spectral file
      INTEGER :: TYPE_ABSORB_SW(NPD_SPECIES)

      ! number of esft terms in each band for each gas
      INTEGER :: I_BAND_ESFT_SW(NPD_BAND, NPD_SPECIES)

      ! type of esft scaling
      INTEGER :: I_SCALE_ESFT_SW(NPD_BAND, NPD_SPECIES)

      ! type of scaling function

      INTEGER :: I_SCALE_FNC_SW(NPD_BAND, NPD_SPECIES)

      !  esft exponents
      REAL :: K_ESFT_SW(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)

      ! esft weights
      REAL :: W_ESFT_SW(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)

      ! scaling parameters for each absorber and term

      REAL :: SCALE_VECTOR_SW(NPD_SCALE_VARIABLE, NPD_ESFT_TERM,        &
     &  NPD_BAND, NPD_SPECIES)

      ! reference pressure for scaling function
      REAL :: P_REFERENCE_SW(NPD_SPECIES, NPD_BAND)

      ! reference temperature for scaling function
      REAL :: T_REFERENCE_SW(NPD_SPECIES, NPD_BAND)

! representation of the Planckian:

      INTEGER :: N_DEG_FIT_SW ! degree of thermal polynomial

      ! coefficients in polynomial fit to source function
      REAL :: THERMAL_COEFFICIENT_SW(0: NPD_THERMAL_COEFF-1,NPD_BAND)

      REAL :: T_REF_PLANCK_SW ! Planckian reference temperature

! surface properties:

      ! method of specifying properties of surface

      INTEGER :: I_SPEC_SURFACE_SW(NPD_SURFACE)

      ! number of parameters fitting the direct albedo
      INTEGER :: N_DIR_ALBEDO_FIT_SW(NPD_SURFACE)


      LOGICAL :: L_SURFACE_SW(NPD_SURFACE) ! surface types included

      ! surface albedos
      REAL :: SURFACE_ALBEDO_SW(NPD_BAND, NPD_SURFACE)

      ! coefficients for fitting direct albedo

      REAL :: DIRECT_ALBEDO_PARM_SW(0: NPD_ALBEDO_PARM, NPD_BAND,       &
     &  NPD_SURFACE)

      ! surface emissivities
      REAL :: EMISSIVITY_GROUND_SW(NPD_BAND, NPD_SURFACE)

! fields for continua:

      ! number of continua in each band
      INTEGER :: N_BAND_CONTINUUM_SW(NPD_BAND)

      ! list of continua in each band
      INTEGER :: INDEX_CONTINUUM_SW(NPD_BAND, NPD_CONTINUUM)

      INTEGER :: INDEX_WATER_SW ! index of water vapour

      ! type of scaling function for continuum
      INTEGER :: I_SCALE_FNC_CONT_SW(NPD_BAND, NPD_CONTINUUM)

      ! grey extinction coefficients for continuum
      REAL :: K_CONTINUUM_SW(NPD_BAND, NPD_CONTINUUM)

      ! scaling parameters for continuum
      REAL :: SCALE_CONTINUUM_SW(NPD_SCALE_VARIABLE,NPD_BAND,           &
     &  NPD_CONTINUUM)

      ! reference pressure for scaling of continuum
      REAL :: P_REF_CONTINUUM_SW(NPD_CONTINUUM, NPD_BAND)

      ! reference temperature for scaling of continuum
      REAL :: T_REF_CONTINUUM_SW(NPD_CONTINUUM, NPD_BAND)

! fields for water droplets:

      ! parametrization type of droplets
      INTEGER :: I_DROP_PARAMETRIZATION_SW(NPD_DROP_TYPE)

      ! types of droplet present
      LOGICAL :: L_DROP_TYPE_SW(NPD_DROP_TYPE)


      ! parameters used to fit optical properties of clouds
      REAL :: DROP_PARAMETER_LIST_SW(NPD_CLOUD_PARAMETER, NPD_BAND,     &
     &  NPD_DROP_TYPE)

      ! minimum dimension permissible in the parametrization
      REAL ::DROP_PARM_MIN_DIM_SW(NPD_DROP_TYPE)

      ! maximum dimension permissible in the parametrization
      REAL :: DROP_PARM_MAX_DIM_SW(NPD_DROP_TYPE)

! fields for aerosols:

      INTEGER :: N_AEROSOL_SW ! number of species of aerosol

      ! types of aerosols
      INTEGER :: TYPE_AEROSOL_SW(NPD_AEROSOL_SPECIES)

      ! parametrization of aerosols
      INTEGER :: I_AEROSOL_PARAMETRIZATION_SW(NPD_AEROSOL_SPECIES)

      ! numbers of humidities
      INTEGER :: NHUMIDITY_SW(NPD_AEROSOL_SPECIES)

      ! aerosol species included
      LOGICAL :: L_AEROSOL_SPECIES_SW(NPD_AEROSOL_SPECIES)

      ! absorption by aerosols
      REAL :: AEROSOL_ABSORPTION_SW(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES,&
     &  NPD_BAND)

      ! scattering by aerosols
      REAL :: AEROSOL_SCATTERING_SW(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES,&
     &  NPD_BAND)

      ! asymmetry of aerosols
      REAL :: AEROSOL_ASYMMETRY_SW(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES, &
     &  NPD_BAND)

      ! humidities for components
      REAL :: HUMIDITIES_SW(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES)

! fields for ice crystals:

      ! types of parametrization of ice crystals
      INTEGER :: I_ICE_PARAMETRIZATION_SW(NPD_ICE_TYPE)

      ! types of ice crystal present
      LOGICAL :: L_ICE_TYPE_SW(NPD_ICE_TYPE)

      ! parameters used to fit single scattering of ice crystals
      REAL :: ICE_PARAMETER_LIST_SW(NPD_CLOUD_PARAMETER, NPD_BAND,      &
     &  NPD_ICE_TYPE)

      ! minimum dimension permissible in the parametrization

      REAL :: ICE_PARM_MIN_DIM_SW(NPD_ICE_TYPE)

      ! maximum dimension permissible in the parametrization
      REAL :: ICE_PARM_MAX_DIM_SW(NPD_ICE_TYPE)

! fields for doppler broadening:

      ! flag for doppler broadening for each species
      LOGICAL :: L_DOPPLER_PRESENT_SW(NPD_SPECIES)

      ! offset to pressure to represent doppler broadening

      REAL :: DOPPLER_CORRECTION_SW(NPD_SPECIES)

! fields for aerosol optical depth

      ! number of wavelengths
      INTEGER :: N_AOD_WAVEL_SW

      ! monochromatic specific coefficient of absorption
      REAL                                                              &
     &     AOD_ABSORPTION_SW(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES,       &
     &                       NPD_AOD_WAVEL)

      ! monochromatic specific coefficient of scattering
      REAL                                                              &
     &     AOD_SCATTERING_SW(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES,       &
     &                       NPD_AOD_WAVEL)

      ! relationship between aerosol component and type
      INTEGER                                                           &
     &     I_AOD_TYPE_SW(NPD_AEROSOL_SPECIES)

! SWSPDL3A end
! SWSPCM3A declares common block containing the reduced sw spectral
! file.(note: SWSPDC3A, SWSPCM3A and SWSARG3A must be consistent)

      COMMON/R2SWSPCM/                                                  &
      ! dimensions of arrays:
     &  NPD_TYPE_SW, NPD_BAND_SW, NPD_EXCLUDE_SW,                       &
     &  NPD_SPECIES_SW, NPD_ESFT_TERM_SW, NPD_SCALE_FNC_SW,             &
     &  NPD_SCALE_VARIABLE_SW,                                          &
     &  NPD_THERMAL_COEFF_SW,                                           &
     &  NPD_SURFACE_SW, NPD_ALBEDO_PARM_SW,                             &
     &  NPD_CONTINUUM_SW,                                               &
     &  NPD_DROP_TYPE_SW, NPD_ICE_TYPE_SW, NPD_CLOUD_PARAMETER_SW,      &
     &  NPD_AEROSOL_SPECIES_SW, NPD_HUMIDITIES_SW,                      &
     &  NPD_AOD_WAVEL_SW,                                               &
      ! general arrays:
     &  L_PRESENT_SW,                                                   &
      ! properties of bands:
     &  N_BAND_SW, WAVE_LENGTH_SHORT_SW, WAVE_LENGTH_LONG_SW,           &
      ! exclusions from bands:
     &  N_BAND_EXCLUDE_SW, INDEX_EXCLUDE_SW,                            &
      ! solar fields:
     &  SOLAR_FLUX_BAND_SW, RAYLEIGH_COEFFICIENT_SW,                    &
      ! gaseous absorption:
     &  N_ABSORB_SW, N_BAND_ABSORB_SW, INDEX_ABSORB_SW,                 &
     &  TYPE_ABSORB_SW,                                                 &
     &  I_BAND_ESFT_SW, I_SCALE_ESFT_SW, I_SCALE_FNC_SW,                &
     &  K_ESFT_SW, W_ESFT_SW,                                           &
     &  SCALE_VECTOR_SW, P_REFERENCE_SW, T_REFERENCE_SW,                &
      ! thermal source function:
     &  N_DEG_FIT_SW, THERMAL_COEFFICIENT_SW, T_REF_PLANCK_SW,          &
      ! surface properties:
     &  I_SPEC_SURFACE_SW, N_DIR_ALBEDO_FIT_SW, L_SURFACE_SW,           &
     &  SURFACE_ALBEDO_SW, DIRECT_ALBEDO_PARM_SW ,                      &
     &  EMISSIVITY_GROUND_SW,                                           &
      ! continua:
     &  N_BAND_CONTINUUM_SW, INDEX_CONTINUUM_SW, INDEX_WATER_SW,        &
     &  I_SCALE_FNC_CONT_SW, K_CONTINUUM_SW,                            &
     &  SCALE_CONTINUUM_SW, P_REF_CONTINUUM_SW, T_REF_CONTINUUM_SW,     &
      ! water droplets:
     &  I_DROP_PARAMETRIZATION_SW, L_DROP_TYPE_SW,                      &
     &  DROP_PARAMETER_LIST_SW,                                         &
     &  DROP_PARM_MIN_DIM_SW, DROP_PARM_MAX_DIM_SW,                     &
      ! aerosols:
     &  N_AEROSOL_SW, TYPE_AEROSOL_SW, I_AEROSOL_PARAMETRIZATION_SW,    &
     &  NHUMIDITY_SW, HUMIDITIES_SW, L_AEROSOL_SPECIES_SW,              &
     &  AEROSOL_ABSORPTION_SW, AEROSOL_SCATTERING_SW,                   &
     &  AEROSOL_ASYMMETRY_SW,                                           &
      ! ice crystals:
     &  I_ICE_PARAMETRIZATION_SW, L_ICE_TYPE_SW,                        &
     &  ICE_PARAMETER_LIST_SW,                                          &
     &  ICE_PARM_MIN_DIM_SW, ICE_PARM_MAX_DIM_SW,                       &
      ! doppler broadening:
     &  L_DOPPLER_PRESENT_SW, DOPPLER_CORRECTION_SW,                    &
      ! aerosol optical depth
     &  N_AOD_WAVEL_SW, AOD_ABSORPTION_SW, AOD_SCATTERING_SW,           &
     &  I_AOD_TYPE_SW

! SWSPCM3A end
!
!
!
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , J
!             LOOP VARIABLE
!
      CHARACTER                                                         &
     &     CH_IOS*5
!             CHARACTER STRING FOR IOS ERROR
!
!
!     SUBROUTINES CALLED
      EXTERNAL                                                          &
     &     R2_COMPRESS_SPECTRUM
!
!
!     EACH BLOCK IS INITIALIZED AS MISSING:
      DATA L_PRESENT/.FALSE., NPD_TYPE*.FALSE./
!
!     INITIALIZE THE RANGE OF VALIDITY OF THE PARAMETRIZATIONS OF
!     DROPLETS AND ICE CRYSTALS. OLD SPECTRAL FILES WILL NOT CONTAIN
!     SUCH DATA, SO THE LIMITS FOR DROPLETS ARE INITIALIZED TO THOSE
!     FORMERLY SET IN THE MICROPHYSICAL SCHEME (MRF/UMIST
!     PARAMETRIZATION) TO ENSURE THAT THE RESULTS ARE BIT-REPRODUCIBLE.
!     VALUES FOR ICE COVER THE RANGE OF EFFECTIVE RADII USED IN
!     GENERATING THE DATA FOR THE ORIGINAL PARAMETRIZATION OF ICE
!     CRYSTALS.
!     AT SOME FUTURE RELEASE IT MAY BE DESIRABLE TO REMOVE DEFAULT
!     SETTINGS.
      DATA DROP_PARM_MIN_DIM/NPD_DROP_TYPE*3.5E-07/
      DATA DROP_PARM_MAX_DIM/NPD_DROP_TYPE*3.7E-05/
      DATA ICE_PARM_MIN_DIM/NPD_ICE_TYPE*3.75E-07/
      DATA ICE_PARM_MAX_DIM/NPD_ICE_TYPE*8.0E-05/
!
!
!
!     READ THE SHORTWAVE SPECTRUM AS A NAMELIST.
      CALL GET_FILE(57, SW_SPECTRAL_FILE, 256, IERR_GET_FILE)
      IF (IERR_GET_FILE /= 0) THEN
!        CONVERT THE ERROR FLAG FROM GET_FILE TO A FLAG RECOGNISED
!        BY THE RADIATION CODE.
         IERR=I_ERR_IO
         CMESSAGE='Error reading name of shortwave spectral file.'
         RETURN
      ENDIF
      OPEN(UNIT=57, FILE=SW_SPECTRAL_FILE, IOSTAT=IOS)
      IF (IOS /= 0) THEN
         IERR=I_ERR_IO
      WRITE(CH_IOS, '(I5)') IOS
         CMESSAGE='Error opening shortwave spectral file.'              &
     &      //' IOSTAT='//CH_IOS
         RETURN
      ENDIF
      READ(57, R2SWSP)
      CLOSE(57)
!
!     TEST FOR MINIMAL REQUISITE INFORMATION.
      IF ( .NOT.(L_PRESENT(0).AND.                                      &
     &           L_PRESENT(2) ) ) THEN
         CMESSAGE='Shortwave spectrum is deficient.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!
!
!     SET REDUCED DIMENSIONS, EITHER FROM THE SIZES OF THE FIXED ARRAYS
!     OR FROM THE ARRAYS READ IN.
!
      NPD_TYPE_SW=NPD_TYPE
      NPD_BAND_SW=MAX(N_BAND, 1)
      NPD_SPECIES_SW=MAX(N_ABSORB, 1)
      NPD_ALBEDO_PARM_SW=NPD_ALBEDO_PARM
      NPD_SCALE_FNC_SW=NPD_SCALE_FNC
      NPD_SCALE_VARIABLE_SW=NPD_SCALE_VARIABLE
      NPD_SURFACE_SW=NPD_SURFACE
      NPD_CONTINUUM_SW=NPD_CONTINUUM
      NPD_CLOUD_PARAMETER_SW=NPD_CLOUD_PARAMETER
      NPD_THERMAL_COEFF_SW=1
      NPD_AOD_WAVEL_SW=1
!
!
!     SEARCH THE SPECTRUM TO FIND MAXIMUM DIMENSIONS.
!
      NPD_EXCLUDE_SW=1
      IF (L_PRESENT(14)) THEN
         DO I=1, N_BAND
            NPD_EXCLUDE_SW=MAX(NPD_EXCLUDE_SW, N_BAND_EXCLUDE(I))
         ENDDO
      ENDIF
!
!     Search the spectrum to find those gases to be retained.
!     Water vapour, carbon dioxide and ozone are included
!     if present, but a warning is printed if they are
!     not included.
      DO I=1, NPD_GASES
         L_GAS_INCLUDED(I)=.FALSE.
      ENDDO
      N_ABSORB_RETAIN=0
!
      DO I=1, N_ABSORB
!
         L_RETAIN_ABSORB(I)=.FALSE.
         COMPRESSED_INDEX(I)=0
!
         IF ( (TYPE_ABSORB(I) == IP_H2O).OR.                            &
     &        (TYPE_ABSORB(I) == IP_CO2).OR.                            &
     &        (TYPE_ABSORB(I) == IP_O3).OR.                             &
     &        ( (TYPE_ABSORB(I) == IP_O2).AND.L_O2 ) ) THEN
            N_ABSORB_RETAIN=N_ABSORB_RETAIN+1
            INDEX_ABSORB_RETAIN(N_ABSORB_RETAIN)=I
            COMPRESSED_INDEX(I)=N_ABSORB_RETAIN
            L_RETAIN_ABSORB(I)=.TRUE.
            L_GAS_INCLUDED(TYPE_ABSORB(I))=.TRUE.
         ENDIF
!
      ENDDO
!
!
!     Print warning messages if those gases normally expected
!     are not present.
      IF (.NOT.L_GAS_INCLUDED(IP_H2O)) THEN
         WRITE(IU_ERR, '(/A, /A)')                                      &
     &      '*** WARNING: Water vapour is not included in the '         &
     &      , 'shortwave spectral file.'
      ENDIF
!
      IF (.NOT.L_GAS_INCLUDED(IP_CO2)) THEN
         WRITE(IU_ERR, '(/A, /A)')                                      &
     &      '*** WARNING: Carbon dioxide is not included in the '       &
     &      , 'shortwave spectral file.'
      ENDIF
!
      IF (.NOT.L_GAS_INCLUDED(IP_O3)) THEN
         WRITE(IU_ERR, '(/A, /A)')                                      &
     &      '*** WARNING: Ozone is not included in the '                &
     &      , 'shortwave spectral file.'
      ENDIF
!
      IF ((.NOT.L_GAS_INCLUDED(IP_O2)).AND.L_O2) THEN
         WRITE(IU_ERR, '(/A, /A)')                                      &
     &      '*** ERROR: Oxygen is not included in the shortwave '       &
     &      , 'spectral file, but was requested in the run.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!     Set an appropriate reduced dimension.
      NPD_SPECIES_SW=MAX(N_ABSORB_RETAIN, 1)
!
!
      NPD_ESFT_TERM_SW=1
      IF (L_PRESENT(5)) THEN
         DO I=1, N_BAND
            DO J=1, N_BAND_ABSORB(I)
               IF (L_RETAIN_ABSORB(INDEX_ABSORB(J, I)))                 &
     &            NPD_ESFT_TERM_SW=MAX(NPD_ESFT_TERM_SW                 &
     &            , I_BAND_ESFT(I, INDEX_ABSORB(J, I)))
            ENDDO
         ENDDO
      ENDIF
!
      NPD_DROP_TYPE_SW=1
      IF (L_PRESENT(10)) THEN
         DO I=1, NPD_DROP_TYPE
            IF (L_DROP_TYPE(I)) THEN
               NPD_DROP_TYPE_SW=MAX(NPD_DROP_TYPE_SW, I)
            ENDIF
         ENDDO
      ENDIF
!
      NPD_ICE_TYPE_SW=1
      IF (L_PRESENT(12)) THEN
         DO I=1, NPD_ICE_TYPE
            IF (L_ICE_TYPE(I)) THEN
               NPD_ICE_TYPE_SW=MAX(NPD_ICE_TYPE_SW, I)
            ENDIF
         ENDDO
      ENDIF
!
!
!     Aerosols must be treated carefully to allow for various
!     different combinations without requiring the spectral file
!     to be too constrained. Only those required will be retained.
!
!     Basic initialization to safe values.
      NPD_HUMIDITIES_SW=1
      N_AEROSOL_RETAIN=0
!
!     Check the spectral file for climatological aerosols
      IF (L_CLIMAT_AEROSOL.OR.L_MURK_RAD) THEN
!
         IF (L_PRESENT(11)) THEN
!
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( (TYPE_AEROSOL(I) == IP_WATER_SOLUBLE).OR.           &
     &              (TYPE_AEROSOL(I) == IP_DUST_LIKE).OR.               &
     &              (TYPE_AEROSOL(I) == IP_OCEANIC).OR.                 &
     &              (TYPE_AEROSOL(I) == IP_SOOT).OR.                    &
     &              (TYPE_AEROSOL(I) == IP_SULPHURIC) ) THEN
                  N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                  INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                  N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF

            ENDDO
!
            IF (N_AEROSOL_FOUND /= 5) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The SW Spectral file lacks some '              &
     &            //'climatological aerosols.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='SW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF

!     Check the spectral file for mineral dust classes.
!     (Only required for the direct effect).
!
      IF (L_USE_DUST .OR. L_USE_ARCLDUST) THEN
!
         IF (L_PRESENT(11)) THEN ! aerosol block present in spec file
!
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( (TYPE_AEROSOL(I) == IP_DUST_1).OR.                  &
     &              (TYPE_AEROSOL(I) == IP_DUST_2).OR.                  &
     &              (TYPE_AEROSOL(I) == IP_DUST_3).OR.                  &
     &              (TYPE_AEROSOL(I) == IP_DUST_4).OR.                  &
     &              (TYPE_AEROSOL(I) == IP_DUST_5).OR.                  &
     &              (TYPE_AEROSOL(I) == IP_DUST_6) ) THEN
                  N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                  INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                  N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF
!
            ENDDO
!
            IF (N_AEROSOL_FOUND /= 6) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The SW Spectral file lacks some '              &
     &            //'mineral dust aerosol.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='SW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!
!     Check the spectral file for sulphate aerosols. (These are
!     required only for the direct effect).
!
      IF (L_USE_SULPC_DIRECT .OR. L_USE_ARCLSULP) THEN
!
         IF (L_PRESENT(11)) THEN
!
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( (TYPE_AEROSOL(I) == IP_ACCUM_SULPHATE).OR.          &
     &              (TYPE_AEROSOL(I) == IP_AITKEN_SULPHATE) ) THEN
                  N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                  INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                  N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF

            ENDDO
!
            IF (N_AEROSOL_FOUND /= 2) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The SW Spectral file lacks some '              &
     &            //'sulphate aerosols.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='SW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!
!     Check the spectral file for soot aerosol modes. (Also only
!     required for the direct effect).
!
      IF (L_USE_SOOT_DIRECT .OR. L_USE_ARCLBLCK) THEN
!
         IF (L_PRESENT(11)) THEN ! aerosol block present in spec file
!
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( (TYPE_AEROSOL(I) == IP_FRESH_SOOT).OR.              &
     &              (TYPE_AEROSOL(I) == IP_AGED_SOOT) ) THEN
                  N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                  INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                  N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF
!
            ENDDO
!
            IF (N_AEROSOL_FOUND /= 2) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The SW Spectral file lacks some '              &
     &            //'soot aerosol.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='SW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!     Check the spectral file for biomass aerosol modes.
!     (Only required for the direct effect).
!
      IF (L_USE_BMASS_DIRECT .OR. L_USE_ARCLBIOM) THEN
!
         IF (L_PRESENT(11)) THEN ! aerosol block present in spec file
!
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( (TYPE_AEROSOL(I) == IP_BIOMASS_1).OR.               &
     &              (TYPE_AEROSOL(I) == IP_BIOMASS_2) ) THEN
                  N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                  INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                  N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF
!
            ENDDO
!
            IF (N_AEROSOL_FOUND /= 2) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The SW Spectral file lacks some '              &
     &            //'biomass aerosol.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='SW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!
!     Check the spectral file for sea-salt aerosol modes. (Also only
!     required for the direct effect).
!
      IF (L_USE_SEASALT_DIRECT .OR. L_USE_ARCLSSLT) THEN
!
         IF (L_PRESENT(11)) THEN ! aerosol block present in spec file
!
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( (TYPE_AEROSOL(I) == IP_SEASALT_FILM).OR.            &
     &              (TYPE_AEROSOL(I) == IP_SEASALT_JET) ) THEN
                  N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                  INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                  N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF
!
            ENDDO
!
            IF (N_AEROSOL_FOUND /= 2) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The SW Spectral file lacks some '              &
     &            //'sea-salt aerosol.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='SW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!     Check the spectral file for biogenic aerosol. (direct effect
!     only).
!
      IF (L_USE_BIOGENIC) THEN
!
         IF (L_PRESENT(11)) THEN ! aerosol block present in spec file
!
!           Search for the aerosol required for this scheme
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( TYPE_AEROSOL(I) == IP_BIOGENIC ) THEN
                 N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                 INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                 N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF
!
            ENDDO
!
            IF (N_AEROSOL_FOUND /= 1) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The SW Spectral file lacks some '              &
     &            //'biogenic aerosol.'
               RETURN
!
            ENDIF

          ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='SW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!     Check the spectral file for fossil-fuel organic carbon aerosol. 
!     (direct effect only).
!
      IF (L_USE_OCFF_DIRECT .OR. L_USE_ARCLOCFF) THEN
!
         IF (L_PRESENT(11)) THEN ! aerosol block present in spec file
!
!           Search for the aerosol required for this scheme
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF (  (TYPE_AEROSOL(I) == IP_OCFF_FRESH) .OR.            &
     &               (TYPE_AEROSOL(I) == IP_OCFF_AGED) ) THEN
                 N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                 INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                 N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF
!
            ENDDO
!
            IF (N_AEROSOL_FOUND /= 2) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The SW Spectral file lacks some '              &
     &            //'fossil fuel org.carb. aerosol.'
               RETURN
!
            ENDIF

          ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='SW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!     Check the spectral file for delta aerosol. (direct effect
!     only).
!
      IF (L_USE_ARCLDLTA) THEN
!
         IF (L_PRESENT(11)) THEN ! aerosol block present in spec file
!
!           Search for the aerosol required for this scheme
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( TYPE_AEROSOL(I) == IP_DELTA ) THEN
                 N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                 INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                 N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF
!
            ENDDO
!
            IF (N_AEROSOL_FOUND /= 1) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The SW Spectral file lacks some '              &
     &            //'delta aerosol.'
               RETURN
!
            ENDIF

          ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='SW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!     Set an appropriate reduced dimension.
      NPD_AEROSOL_SPECIES_SW=MAX(N_AEROSOL_RETAIN, 1)
!
!     Set the allowed number of humidities from the number of
!     retained aerosols.
!
      IF (L_PRESENT(11)) THEN
         DO I=1, N_AEROSOL_RETAIN
            IF (I_AEROSOL_PARAMETRIZATION(INDEX_AEROSOL_RETAIN(I)) ==   &
     &         IP_AEROSOL_PARAM_MOIST) THEN
               NPD_HUMIDITIES_SW=MAX(NPD_HUMIDITIES_SW                  &
     &            , NHUMIDITY(INDEX_AEROSOL_RETAIN(I)))
            ENDIF
         ENDDO
      ENDIF
!
!
!
!
!     TRANSFER THE LARGE NAMELIST TO THE REDUCED SPECTRUM.
!
!
! DEPENDS ON: r2_compress_spectrum
      CALL R2_COMPRESS_SPECTRUM(                                        &
!                       Spectral Array in Namelist
     &     L_PRESENT                                                    &
     &   , N_BAND, WAVE_LENGTH_SHORT , WAVE_LENGTH_LONG                 &
     &   , N_BAND_EXCLUDE, INDEX_EXCLUDE                                &
     &   , SOLAR_FLUX_BAND, RAYLEIGH_COEFFICIENT                        &
     &   , N_ABSORB, N_BAND_ABSORB, INDEX_ABSORB, TYPE_ABSORB           &
     &   , L_RETAIN_ABSORB, N_ABSORB_RETAIN, INDEX_ABSORB_RETAIN        &
     &   , COMPRESSED_INDEX, I_BAND_ESFT, K_ESFT, W_ESFT, I_SCALE_ESFT  &
     &   , I_SCALE_FNC, SCALE_VECTOR, P_REFERENCE, T_REFERENCE          &
     &   , N_DEG_FIT, THERMAL_COEFFICIENT, T_REF_PLANCK                 &
     &   , I_SPEC_SURFACE, L_SURFACE, SURFACE_ALBEDO                    &
     &   , N_DIR_ALBEDO_FIT, DIRECT_ALBEDO_PARM, EMISSIVITY_GROUND      &
     &   , N_BAND_CONTINUUM, INDEX_CONTINUUM, INDEX_WATER               &
     &   , K_CONTINUUM, I_SCALE_FNC_CONT, SCALE_CONTINUUM               &
     &   , P_REF_CONTINUUM, T_REF_CONTINUUM                             &
     &   , L_DROP_TYPE, I_DROP_PARAMETRIZATION, DROP_PARAMETER_LIST     &
     &   , DROP_PARM_MIN_DIM, DROP_PARM_MAX_DIM                         &
     &   , L_ICE_TYPE, I_ICE_PARAMETRIZATION, ICE_PARAMETER_LIST        &
     &   , ICE_PARM_MIN_DIM, ICE_PARM_MAX_DIM                           &
     &   , N_AEROSOL, TYPE_AEROSOL                                      &
     &   , N_AEROSOL_RETAIN, INDEX_AEROSOL_RETAIN                       &
     &   , L_AEROSOL_SPECIES, AEROSOL_ABSORPTION                        &
     &   , AEROSOL_SCATTERING, AEROSOL_ASYMMETRY                        &
     &   , NHUMIDITY, HUMIDITIES, I_AEROSOL_PARAMETRIZATION             &
     &   , L_DOPPLER_PRESENT, DOPPLER_CORRECTION, .FALSE.               &
     &   , N_AOD_WAVEL, AOD_ABSORPTION, AOD_SCATTERING, I_AOD_TYPE      &
!                       Reduced Spectral Array
!     ------------------------------------------------------------------
!     ARGUMENT LIST FOR THE REDUCED SW SPECTRAL FILE.
!     (NOTE: SWSPDC3A, SWSPCM3A AND SWSARG3A MUST BE CONSISTENT)
!
     &   , NPD_TYPE_SW, NPD_BAND_SW, NPD_EXCLUDE_SW                     &
     &   , NPD_SPECIES_SW, NPD_ESFT_TERM_SW, NPD_SCALE_FNC_SW           &
     &   , NPD_SCALE_VARIABLE_SW, NPD_THERMAL_COEFF_SW                  &
     &   , NPD_SURFACE_SW, NPD_ALBEDO_PARM_SW, NPD_CONTINUUM_SW         &
     &   , NPD_DROP_TYPE_SW, NPD_ICE_TYPE_SW, NPD_CLOUD_PARAMETER_SW    &
     &   , NPD_AEROSOL_SPECIES_SW, NPD_HUMIDITIES_SW, NPD_AOD_WAVEL_SW  &
     &   , L_PRESENT_SW, N_BAND_SW, WAVE_LENGTH_SHORT_SW                &
     &   , WAVE_LENGTH_LONG_SW, N_BAND_EXCLUDE_SW, INDEX_EXCLUDE_SW     &
     &   , SOLAR_FLUX_BAND_SW, RAYLEIGH_COEFFICIENT_SW                  &
     &   , N_ABSORB_SW, N_BAND_ABSORB_SW, INDEX_ABSORB_SW               &
     &   , TYPE_ABSORB_SW, I_BAND_ESFT_SW, I_SCALE_ESFT_SW              &
     &   , I_SCALE_FNC_SW, K_ESFT_SW, W_ESFT_SW                         &
     &   , SCALE_VECTOR_SW, P_REFERENCE_SW, T_REFERENCE_SW              &
     &   , N_DEG_FIT_SW, THERMAL_COEFFICIENT_SW, T_REF_PLANCK_SW        &
     &   , I_SPEC_SURFACE_SW, N_DIR_ALBEDO_FIT_SW, L_SURFACE_SW         &
     &   , SURFACE_ALBEDO_SW, DIRECT_ALBEDO_PARM_SW                     &
     &   , EMISSIVITY_GROUND_SW, N_BAND_CONTINUUM_SW, INDEX_CONTINUUM_SW&
     &   , INDEX_WATER_SW, I_SCALE_FNC_CONT_SW, K_CONTINUUM_SW          &
     &   , SCALE_CONTINUUM_SW, P_REF_CONTINUUM_SW, T_REF_CONTINUUM_SW   &
     &   , I_DROP_PARAMETRIZATION_SW, L_DROP_TYPE_SW                    &
     &   , DROP_PARAMETER_LIST_SW                                       &
     &   , DROP_PARM_MIN_DIM_SW, DROP_PARM_MAX_DIM_SW                   &
     &   , N_AEROSOL_SW, TYPE_AEROSOL_SW, I_AEROSOL_PARAMETRIZATION_SW  &
     &   , NHUMIDITY_SW, HUMIDITIES_SW, L_AEROSOL_SPECIES_SW            &
     &   , AEROSOL_ABSORPTION_SW, AEROSOL_SCATTERING_SW                 &
     &   , AEROSOL_ASYMMETRY_SW                                         &
     &   , I_ICE_PARAMETRIZATION_SW, L_ICE_TYPE_SW                      &
     &   , ICE_PARAMETER_LIST_SW                                        &
     &   , ICE_PARM_MIN_DIM_SW, ICE_PARM_MAX_DIM_SW                     &
     &   , L_DOPPLER_PRESENT_SW, DOPPLER_CORRECTION_SW                  &
     &   , N_AOD_WAVEL_SW, AOD_ABSORPTION_SW, AOD_SCATTERING_SW         &
     &   , I_AOD_TYPE_SW                                                &
!
!     ------------------------------------------------------------------
     &   )
!
!
!
      RETURN
      END SUBROUTINE R2_SW_SPECIN
