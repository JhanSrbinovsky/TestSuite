

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read a longwave spectral namelist.
!
! Purpose:
!   To read a longwave namelist into a spectral array.
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
!
!       4.4             02-09-97                Aerosol flags passed
!                                               in to the code to
!                                               enable only those
!                                               required to be
!                                               selected. Spectral
!                                               data are no longer
!                                               compressed into a
!                                               single array.
!                                               IOSTAT error code
!                                               returned as part of
!                                               CMESSAGE.
!                                               (J. M. Edwards)
!       4.5        April 1998   Allow soot spectral data to be read.
!                                                     Luke Robinson.
!       4.5             18-05-98                Coding to allow
!                                               selection of gases
!                                               from the spectral
!                                               file.
!                                               (J. M. Edwards)
!       5.3             04-04-01                Allow aerosol data to be
!                                               read in if required by
!                                               mesoscale model.
!                                               (S. Cusack)
!       6.2             15-11-05                Added support of
!                                               block 15 (aerosol
!                                               optical depth) in the
!                                               LW spectral file.
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_LW_SPECIN(IERR, CMESSAGE                            &
     &   , L_CH4, L_N2O, L_CFC11, L_CFC12                               &
     &   , L_CFC113, L_HCFC22, L_HFC125, L_HFC134A                      &
     &   , L_CLIMAT_AEROSOL                                             &
     &   , L_USE_DUST, L_USE_ARCLDUST                                   &
     &   , L_USE_SULPC_DIRECT, L_USE_ARCLSULP                           &
     &   , L_USE_SOOT_DIRECT, L_USE_ARCLBLCK                            &
     &   , L_USE_BMASS_DIRECT, L_USE_ARCLBIOM                           &
     &   , L_USE_SEASALT_DIRECT, L_USE_ARCLSSLT                         &
     &   , L_USE_OCFF_DIRECT, L_USE_ARCLOCFF                            &
     &   , L_USE_BIOGENIC, L_USE_ARCLDLTA                               &
     &   , L_MURK_RAD                                                   &
     &   , L_USE_AOD                                                    &
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
     &     L_CH4                                                        &
!             ABSORPTION BY METHANE IS INCLUDED
     &   , L_N2O                                                        &
!             ABSORPTION BY NITROUS OXIDE IS INCLUDED
     &   , L_CFC11                                                      &
!             ABSORPTION BY CFC11 IS INCLUDED
     &   , L_CFC12                                                      &
!             ABSORPTION BY CFC12 IS INCLUDED
     &   , L_CFC113                                                     &
!             ABSORPTION BY CFC113 IS INCLUDED
     &   , L_HCFC22                                                     &
!             ABSORPTION BY HCFC22 IS INCLUDED
     &   , L_HFC125                                                     &
!             ABSORPTION BY HFC125 IS INCLUDED
     &   , L_HFC134A                                                    &
!             ABSORPTION BY HFC134A IS INCLUDED
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
!             USE THE DIRECT RAD EFFECTS OF SOOT IN THE LW
     &   , L_USE_ARCLBLCK                                               &
                 ! direct effect of black-carbon aerosol from NWP clim
     &   , L_USE_BMASS_DIRECT                                           &
!             USE THE DIRECT RAD EFFECTS OF BIOMASS SMOKE IN THE LW
     &   , L_USE_ARCLBIOM                                               &
                 ! direct effect of biomass aerosol from NWP clim
     &   , L_USE_SEASALT_DIRECT                                         &
!             USE THE DIRECT RAD EFFECTS OF SEASALT IN THE LW
     &   , L_USE_ARCLSSLT                                               &
                 ! direct effect of sea-salt aerosol from NWP clim
     &   , L_USE_BIOGENIC                                               &
!             USE THE BIOGENIC AEROSOL DIRECT EFFECT IN THE LW
     &   , L_USE_OCFF_DIRECT                                            &
                 ! direct effect of fossil-fuel org.carb in the LW
     &   , L_USE_ARCLOCFF                                               &
                 ! direct effect of fossil-fuel org.carb. from NWP clim
     &   , L_USE_ARCLDLTA                                               &
                 ! direct effect of delta aerosol from NWP climatology
     &   , L_MURK_RAD                                                   &
!             MESOSCALE MODEL AEROSOLS ARE TO BE INCLUDED
     &   , L_USE_AOD
!             AT LEAST ONE OF THE AEROSOL OPTICAL DEPTH DIAGNOSTICS
!             IS REQUESTED
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
     &     LW_SPECTRAL_FILE
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
     &   , N_AEROSOL_RETAIN                                             &
!             NUMBER OF AEROSOLS IN THE SPECTRAL FILE TO BE RETAINED
!             FOR THE RADIATIVE CALCULATION
     &   , INDEX_AEROSOL_RETAIN(NPD_AEROSOL_SPECIES)                    &
!             INDEXING NUMBERS OF THE RETAINED AEROSOLS
     &   , COMPRESSED_INDEX(NPD_SPECIES)                                &
!             MAPPING FROM OLD TO NEW INDICES OF ABSORBERS
     &   , N_AEROSOL_FOUND
!             NUMBER OF AEROSOLS FOR THE CURRENT GROUP OF PROCESSES
!             FOUND IN THE SPECTRAL FILE
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
! LWSP3A declares elements of a spectral file as a namelist.
!
! 6.2  24/10/05 Functionality for radiative forcing, timestepping
!               and radiances under versions 3C and 3Z of radiation
!               code added                 (J.-C. Thelen)
!

      NAMELIST/R2LWSP/                                                  &
        ! Blocks Present
     &  L_PRESENT,                                                      &
        ! Block 0
     &  N_BAND, N_ABSORB, N_AEROSOL, TYPE_ABSORB, TYPE_AEROSOL,         &
        ! Block 1
     &  WAVE_LENGTH_SHORT, WAVE_LENGTH_LONG,                            &
        ! Block 2
     &  SOLAR_FLUX_BAND,                                                &
        ! Block 3
     &  RAYLEIGH_COEFFICIENT,                                           &
        ! Block 4
     &  N_BAND_ABSORB, INDEX_ABSORB,                                    &
        ! Block 5
     &  I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC,                         &
     &  P_REFERENCE, T_REFERENCE, K_ESFT, W_ESFT, SCALE_VECTOR,         &
        ! Block 6
     &  N_DEG_FIT, T_REF_PLANCK, THERMAL_COEFFICIENT,                   &
        ! Block 7
     &  I_SPEC_SURFACE, N_DIR_ALBEDO_FIT, L_SURFACE,                    &
     &  SURFACE_ALBEDO, DIRECT_ALBEDO_PARM, EMISSIVITY_GROUND,          &
        ! Block 8
     &  N_BAND_CONTINUUM, INDEX_CONTINUUM, INDEX_WATER,                 &
        ! Block 9
     &  I_SCALE_FNC_CONT, P_REF_CONTINUUM, T_REF_CONTINUUM,             &
     &  K_CONTINUUM, SCALE_CONTINUUM,                                   &
        ! Block 10
     &  I_DROP_PARAMETRIZATION, L_DROP_TYPE, DROP_PARAMETER_LIST,       &
     &  DROP_PARM_MIN_DIM, DROP_PARM_MAX_DIM,                           &
        ! Block 11
     &  I_AEROSOL_PARAMETRIZATION, NHUMIDITY, L_AEROSOL_SPECIES,        &
     &  AEROSOL_ABSORPTION, AEROSOL_SCATTERING, AEROSOL_ASYMMETRY,      &
     &  HUMIDITIES,                                                     &
        ! Block 12
     &  I_ICE_PARAMETRIZATION, L_ICE_TYPE, ICE_PARAMETER_LIST,          &
     &  ICE_PARM_MIN_DIM, ICE_PARM_MAX_DIM,                             &
        ! Block 13
     &  L_DOPPLER_PRESENT, DOPPLER_CORRECTION,                          &
        ! Block 14
     &  N_BAND_EXCLUDE, INDEX_EXCLUDE,                                  &
        ! Block 15
     &  N_AOD_WAVEL, AOD_ABSORPTION, AOD_SCATTERING, I_AOD_TYPE

! LWSP3A end
!
!
!     DECLARE THE REDUCED LW SPECTRAL FILE AND ITS HOLDING COMMON BLOCK.
!
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE CONTAINING DECLARATIONS FOR REDUCED LW-SPECTRAL FILE.
!     NOTE: LWSPDC3A, LWSPCM3A AND LWSARG3A MUST BE CONSISTENT.
!     NOTE: SINCE THE ARRAYS HERE WILL BE PASSED IN A COMMON BLOCK
!     THEIR SIZES MUST BE FIXED, EVEN THOUGH VARIABLE SIZES ARE USED
!     LOWER IN THE CODE. THEY ARE ACCORDINGLY DEFINED AS 1-DIMENSIONAL
!     ARRAYS WITH FIXED MAXIMUM SIZES AT THIS LEVEL.
!
!     ------------------------------------------------------------------
!
!
!     DIMENSIONS FOR THE REDUCED LW SPECTRAL FILE
!
      INTEGER                                                           &
     &     NPD_TYPE_LW                                                  &
!             NUMBER OF TYPES OF DATA IN LW SPECTRUM
     &   , NPD_BAND_LW                                                  &
!             NUMBER OF SPECTRAL BANDS IN LW SPECTRUM
     &   , NPD_EXCLUDE_LW                                               &
!             NUMBER OF EXCLUDED BANDS IN LW SPECTRUM
     &   , NPD_SPECIES_LW                                               &
!             NUMBER OF GASEOUS SPECIES IN LW SPECTRUM
     &   , NPD_ESFT_TERM_LW                                             &
!             NUMBER OF ESFT TERMS IN LW SPECTRUM
     &   , NPD_SCALE_FNC_LW                                             &
!             NUMBER OF SCALING FUNCTIONS IN LW SPECTRUM
     &   , NPD_SCALE_VARIABLE_LW                                        &
!             NUMBER OF SCALING VARIABLES IN LW SPECTRUM
     &   , NPD_SURFACE_LW                                               &
!             NUMBER OF SURFACE TYPES IN LW SPECTRUM
     &   , NPD_ALBEDO_PARM_LW                                           &
!             NUMBER OF ALBEDO PARAMETERS IN LW SPECTRUM
     &   , NPD_CONTINUUM_LW                                             &
!             NUMBER OF CONTINUA IN LW SPECTRUM
     &   , NPD_DROP_TYPE_LW                                             &
!             NUMBER OF DROP TYPES IN LW SPECTRUM
     &   , NPD_ICE_TYPE_LW                                              &
!             NUMBER OF ICE CRYSTAL TYPES IN LW SPECTRUM
     &   , NPD_AEROSOL_SPECIES_LW                                       &
!             NUMBER OF AEROSOL SPECIES IN LW SPECTRUM
     &   , NPD_CLOUD_PARAMETER_LW                                       &
!             MAX NUMBER OF CLOUD PARAMETERS IN LW SPECTRUM
     &   , NPD_HUMIDITIES_LW                                            &
!             MAXIMUM NUMBER OF HUMIDITIES IN LW SPECTRUM
     &   , NPD_THERMAL_COEFF_LW                                         &
!             NUMBER OF THERMAL COEFFICIENTS IN LW SPECTRUM
     &   , NPD_AOD_WAVEL_LW
!             NUMBER OF WAVELENGTHS FOR AOD IN LW SPECTRUM
!
!
!     GENERAL FIELDS:
!
      LOGICAL                                                           &
     &     L_PRESENT_LW(0: NPD_TYPE)
!             FLAG FOR TYPES OF DATA PRESENT
!
!
!
!     PROPERTIES OF THE SPECTRAL BANDS:
!
      INTEGER                                                           &
     &     N_BAND_LW
!             NUMBER OF SPECTRAL BANDS
!
      REAL                                                              &
     &     WAVE_LENGTH_SHORT_LW(NPD_BAND)                               &
!             SHORTER WAVELENGTH LIMITS
     &   , WAVE_LENGTH_LONG_LW(NPD_BAND)
!             LONGER WAVELENGTH LIMITS
!
!
!
!     EXCLUSION OF SPECIFIC BANDS FROM PARTS OF THE SPECTRUM:
!
      INTEGER                                                           &
     &     N_BAND_EXCLUDE_LW(NPD_BAND)                                  &
!             NUMBER OF EXCLUDED BANDS WITHIN EACH SPECTRAL BAND
     &   , INDEX_EXCLUDE_LW(NPD_EXCLUDE, NPD_BAND)
!             INDICES OF EXCLUDED BANDS
!
!
!
!     FIELDS FOR THE SOLAR FLUX:
!
      REAL                                                              &
     &     SOLAR_FLUX_BAND_LW(NPD_BAND)
!             FRACTION OF THE INCIDENT SOLAR FLUX IN EACH BAND
!
!
!
!     FIELDS FOR RAYLEIGH SCATTERING:
!
      REAL                                                              &
     &     RAYLEIGH_COEFFICIENT_LW(NPD_BAND)
!             RAYLEIGH COEFFICIENTS
!
!
!
!     FIELDS FOR GASEOUS ABSORPTION:
!
      INTEGER                                                           &
     &     N_ABSORB_LW                                                  &
!             NUMBER OF ABSORBERS
     &   , N_BAND_ABSORB_LW(NPD_BAND)                                   &
!             NUMBER OF ABSORBERS IN EACH BAND
     &   , INDEX_ABSORB_LW(NPD_SPECIES, NPD_BAND)                       &
!             LIST OF ABSORBERS IN EACH BAND
     &   , TYPE_ABSORB_LW(NPD_SPECIES)                                  &
!             TYPES OF EACH GAS IN THE SPECTRAL FILE
     &   , I_BAND_ESFT_LW(NPD_BAND, NPD_SPECIES)                        &
!             NUMBER OF ESFT TERMS IN EACH BAND FOR EACH GAS
     &   , I_SCALE_ESFT_LW(NPD_BAND, NPD_SPECIES)                       &
!             TYPE OF ESFT SCALING
     &   , I_SCALE_FNC_LW(NPD_BAND, NPD_SPECIES)
!             TYPE OF SCALING FUNCTION
!
      REAL                                                              &
     &     K_ESFT_LW(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)              &
!             ESFT EXPONENTS
     &   , W_ESFT_LW(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)              &
!             ESFT WEIGHTS
     &   , SCALE_VECTOR_LW(NPD_SCALE_VARIABLE, NPD_ESFT_TERM            &
     &        , NPD_BAND, NPD_SPECIES)                                  &
!             SCALING PARAMETERS FOR EACH ABSORBER AND TERM
     &   , P_REFERENCE_LW(NPD_SPECIES, NPD_BAND)                        &
!             REFERENCE PRESSURE FOR SCALING FUNCTION
     &   , T_REFERENCE_LW(NPD_SPECIES, NPD_BAND)
!             REFERENCE TEMPERATURE FOR SCALING FUNCTION
!
!
!
!     REPRESENTATION OF THE PLANCKIAN:
!
      INTEGER                                                           &
     &     N_DEG_FIT_LW
!             DEGREE OF THERMAL POLYNOMIAL
!
      REAL                                                              &
     &     THERMAL_COEFFICIENT_LW(0: NPD_THERMAL_COEFF-1                &
     &        , NPD_BAND)                                               &
!             COEFFICIENTS IN POLYNOMIAL FIT TO SOURCE FUNCTION
     &   , T_REF_PLANCK_LW
!             PLANCKIAN REFERENCE TEMPERATURE
!
!
!
!     SURFACE PROPERTIES:
!
      INTEGER                                                           &
     &     I_SPEC_SURFACE_LW(NPD_SURFACE)                               &
!             METHOD OF SPECIFYING PROPERTIES OF SURFACE
     &   , N_DIR_ALBEDO_FIT_LW(NPD_SURFACE)
!             NUMBER OF PARAMETERS FITTING THE DIRECT ALBEDO
!
      LOGICAL                                                           &
     &     L_SURFACE_LW(NPD_SURFACE)
!             SURFACE TYPES INCLUDED
!
      REAL                                                              &
     &     SURFACE_ALBEDO_LW(NPD_BAND, NPD_SURFACE)                     &
!             SURFACE ALBEDOS
     &   , DIRECT_ALBEDO_PARM_LW(0: NPD_ALBEDO_PARM                     &
     &      , NPD_BAND, NPD_SURFACE)                                    &
!             COEFFICIENTS FOR FITTING DIRECT ALBEDO
     &   , EMISSIVITY_GROUND_LW(NPD_BAND, NPD_SURFACE)
!             SURFACE EMISSIVITIES
!
!
!
!     FIELDS FOR CONTINUA:
!
      INTEGER                                                           &
     &     N_BAND_CONTINUUM_LW(NPD_BAND)                                &
!             NUMBER OF CONTINUA IN EACH BAND
     &   , INDEX_CONTINUUM_LW(NPD_BAND, NPD_CONTINUUM)                  &
!             LIST OF CONTINUA IN EACH BAND
     &   , INDEX_WATER_LW                                               &
!             INDEX OF WATER VAPOUR
     &   , I_SCALE_FNC_CONT_LW(NPD_BAND, NPD_CONTINUUM)
!             TYPE OF SCALING FUNCTION FOR CONTINUUM
!
      REAL                                                              &
     &     K_CONTINUUM_LW(NPD_BAND, NPD_CONTINUUM)                      &
!             GREY EXTINCTION COEFFICIENTS FOR CONTINUUM
     &   , SCALE_CONTINUUM_LW(NPD_SCALE_VARIABLE                        &
     &      , NPD_BAND, NPD_CONTINUUM)                                  &
!             SCALING PARAMETERS FOR CONTINUUM
     &   , P_REF_CONTINUUM_LW(NPD_CONTINUUM, NPD_BAND)                  &
!             REFERENCE PRESSURE FOR SCALING OF CONTINUUM
     &   , T_REF_CONTINUUM_LW(NPD_CONTINUUM, NPD_BAND)
!             REFERENCE TEMPERATURE FOR SCALING OF CONTINUUM
!
!
!
!     FIELDS FOR WATER DROPLETS:
!
      INTEGER                                                           &
     &     I_DROP_PARAMETRIZATION_LW(NPD_DROP_TYPE)
!             PARAMETRIZATION TYPE OF DROPLETS
!
      LOGICAL                                                           &
     &     L_DROP_TYPE_LW(NPD_DROP_TYPE)
!             TYPES OF DROPLET PRESENT
!
      REAL                                                              &
     &     DROP_PARAMETER_LIST_LW(NPD_CLOUD_PARAMETER                   &
     &        , NPD_BAND, NPD_DROP_TYPE)                                &
!             PARAMETERS USED TO FIT OPTICAL PROPERTIES OF CLOUDS
     &   , DROP_PARM_MIN_DIM_LW(NPD_DROP_TYPE)                          &
!             MINIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
     &   , DROP_PARM_MAX_DIM_LW(NPD_DROP_TYPE)
!             MAXIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR AEROSOLS:
!
      INTEGER                                                           &
     &     N_AEROSOL_LW                                                 &
!             NUMBER OF SPECIES OF AEROSOL
     &   , TYPE_AEROSOL_LW(NPD_AEROSOL_SPECIES)                         &
!             TYPES OF AEROSOLS
     &   , I_AEROSOL_PARAMETRIZATION_LW(NPD_AEROSOL_SPECIES)            &
!             PARAMETRIZATION OF AEROSOLS
     &   , NHUMIDITY_LW(NPD_AEROSOL_SPECIES)
!             NUMBERS OF HUMIDITIES
!
      LOGICAL                                                           &
     &     L_AEROSOL_SPECIES_LW(NPD_AEROSOL_SPECIES)
!             AEROSOL SPECIES INCLUDED
!
      REAL                                                              &
     &     AEROSOL_ABSORPTION_LW(NPD_HUMIDITIES                         &
     &        , NPD_AEROSOL_SPECIES, NPD_BAND)                          &
!             ABSORPTION BY AEROSOLS
     &   , AEROSOL_SCATTERING_LW(NPD_HUMIDITIES                         &
     &        , NPD_AEROSOL_SPECIES, NPD_BAND)                          &
!             SCATTERING BY AEROSOLS
     &   , AEROSOL_ASYMMETRY_LW(NPD_HUMIDITIES                          &
     &        , NPD_AEROSOL_SPECIES, NPD_BAND)                          &
!             ASYMMETRY OF AEROSOLS
     &   , HUMIDITIES_LW(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES)
!             HUMIDITIES FOR COMPONENTS
!
!
!
!     FIELDS FOR ICE CRYSTALS:
!
      INTEGER                                                           &
     &     I_ICE_PARAMETRIZATION_LW(NPD_ICE_TYPE)
!             TYPES OF PARAMETRIZATION OF ICE CRYSTALS
!
      LOGICAL                                                           &
     &     L_ICE_TYPE_LW(NPD_ICE_TYPE)
!             TYPES OF ICE CRYSTAL PRESENT
!
      REAL                                                              &
     &     ICE_PARAMETER_LIST_LW(NPD_CLOUD_PARAMETER                    &
     &        , NPD_BAND, NPD_ICE_TYPE)                                 &
!             PARAMETERS USED TO FIT SINGLE SCATTERING OF ICE CRYSTALS
     &   , ICE_PARM_MIN_DIM_LW(NPD_ICE_TYPE)                            &
!             MINIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
     &   , ICE_PARM_MAX_DIM_LW(NPD_ICE_TYPE)
!             MAXIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR DOPPLER BROADENING:
!
      LOGICAL                                                           &
     &     L_DOPPLER_PRESENT_LW(NPD_SPECIES)
!             FLAG FOR DOPPLER BROADENING FOR EACH SPECIES
!
      REAL                                                              &
     &     DOPPLER_CORRECTION_LW(NPD_SPECIES)
!             OFFSET TO PRESSURE TO REPRESENT DOPPLER BROADENING
!
!
!
!     FIELDS FOR AEROSOL OPTICAL DEPTH
!
      INTEGER                                                           &
     &     N_AOD_WAVEL_LW
!             NUMBER OF WAVELENGTHS
!
      REAL                                                              &
     &     AOD_ABSORPTION_LW(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES,       &
     &                       NPD_AOD_WAVEL)
!             MONOCHROMATIC SPECIFIC COEFFICIENT OF ABSORPTION
!
      REAL                                                              &
     &     AOD_SCATTERING_LW(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES,       &
     &                       NPD_AOD_WAVEL)
!             MONOCHROMATIC SPECIFIC COEFFICIENT OF SCATTERING
!
      INTEGER                                                           &
     &     I_AOD_TYPE_LW(NPD_AEROSOL_SPECIES)
!             RELATIONSHIP BETWEEN AEROSOL COMPONENT AND TYPE
!
!
!
!
!    ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     MODULE DECLARING COMMON BLOCK CONTAINING THE REDUCED LW SPECTRAL
!     FILE.
!     (NOTE: LWSPDC3A, LWSPCM3A AND LWSARG3A MUST BE CONSISTENT)
!
!
      COMMON/R2LWSPCM/                                                  &
!
!     DIMENSIONS OF ARRAYS:
     &     NPD_TYPE_LW, NPD_BAND_LW, NPD_EXCLUDE_LW                     &
     &   , NPD_SPECIES_LW, NPD_ESFT_TERM_LW, NPD_SCALE_FNC_LW           &
     &   , NPD_SCALE_VARIABLE_LW                                        &
     &   , NPD_THERMAL_COEFF_LW                                         &
     &   , NPD_SURFACE_LW, NPD_ALBEDO_PARM_LW                           &
     &   , NPD_CONTINUUM_LW                                             &
     &   , NPD_DROP_TYPE_LW, NPD_ICE_TYPE_LW, NPD_CLOUD_PARAMETER_LW    &
     &   , NPD_AEROSOL_SPECIES_LW, NPD_HUMIDITIES_LW, NPD_AOD_WAVEL_LW  &
!
!     GENERAL ARRAYS:
     &   , L_PRESENT_LW                                                 &
!
!     PROPERTIES OF BANDS:
     &   , N_BAND_LW, WAVE_LENGTH_SHORT_LW, WAVE_LENGTH_LONG_LW         &
!
!     EXCLUSIONS FROM BANDS:
     &   , N_BAND_EXCLUDE_LW, INDEX_EXCLUDE_LW                          &
!
!     SOLAR FIELDS:
     &   , SOLAR_FLUX_BAND_LW, RAYLEIGH_COEFFICIENT_LW                  &
!
!     GASEOUS ABSORPTION:
     &   , N_ABSORB_LW, N_BAND_ABSORB_LW, INDEX_ABSORB_LW               &
     &   , TYPE_ABSORB_LW                                               &
     &   , I_BAND_ESFT_LW, I_SCALE_ESFT_LW, I_SCALE_FNC_LW              &
     &   , K_ESFT_LW, W_ESFT_LW                                         &
     &   , SCALE_VECTOR_LW, P_REFERENCE_LW, T_REFERENCE_LW              &
!
!     THERMAL SOURCE FUNCTION:
     &   , N_DEG_FIT_LW, THERMAL_COEFFICIENT_LW, T_REF_PLANCK_LW        &
!
!     SURFACE PROPERTIES:
     &   , I_SPEC_SURFACE_LW, N_DIR_ALBEDO_FIT_LW, L_SURFACE_LW         &
     &   , SURFACE_ALBEDO_LW, DIRECT_ALBEDO_PARM_LW                     &
     &   , EMISSIVITY_GROUND_LW                                         &
!
!     CONTINUA:
     &   , N_BAND_CONTINUUM_LW, INDEX_CONTINUUM_LW, INDEX_WATER_LW      &
     &   , I_SCALE_FNC_CONT_LW, K_CONTINUUM_LW                          &
     &   , SCALE_CONTINUUM_LW, P_REF_CONTINUUM_LW, T_REF_CONTINUUM_LW   &
!
!     WATER DROPLETS:
     &   , I_DROP_PARAMETRIZATION_LW, L_DROP_TYPE_LW                    &
     &   , DROP_PARAMETER_LIST_LW                                       &
     &   , DROP_PARM_MIN_DIM_LW, DROP_PARM_MAX_DIM_LW                   &
!
!     AEROSOLS:
     &   , N_AEROSOL_LW, TYPE_AEROSOL_LW, I_AEROSOL_PARAMETRIZATION_LW  &
     &   , NHUMIDITY_LW, HUMIDITIES_LW, L_AEROSOL_SPECIES_LW            &
     &   , AEROSOL_ABSORPTION_LW, AEROSOL_SCATTERING_LW                 &
     &   , AEROSOL_ASYMMETRY_LW                                         &
!
!     ICE CRYSTALS:
     &   , I_ICE_PARAMETRIZATION_LW, L_ICE_TYPE_LW                      &
     &   , ICE_PARAMETER_LIST_LW                                        &
     &   , ICE_PARM_MIN_DIM_LW, ICE_PARM_MAX_DIM_LW                     &
!
!     DOPPLER BROADENING:
     &   , L_DOPPLER_PRESENT_LW, DOPPLER_CORRECTION_LW                  &
!
!     AEROSOL OPTICAL DEPTH
     &   , N_AOD_WAVEL_LW, AOD_ABSORPTION_LW, AOD_SCATTERING_LW         &
     &   , I_AOD_TYPE_LW
!
!     ------------------------------------------------------------------
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
!             CHARACTER STRING FOR IOSTAT ERROR
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
!     READ THE LONGWAVE SPECTRUM AS A NAMELIST.
      CALL GET_FILE(80, LW_SPECTRAL_FILE, 256, IERR_GET_FILE)
      IF (IERR_GET_FILE /= 0) THEN
!        CONVERT THE ERROR FLAG FROM GET_FILE TO A FLAG RECOGNISED
!        BY THE RADIATION CODE.
         IERR=I_ERR_IO
         CMESSAGE='Error reading name of longwave spectral file.'
         RETURN
      ENDIF
      OPEN(UNIT=80, FILE=LW_SPECTRAL_FILE, IOSTAT=IOS)
      IF (IOS /= 0) THEN
         IERR=I_ERR_IO
      WRITE(CH_IOS, '(I5)') IOS
        CMESSAGE='Error opening longwave spectral file.'                &
     &      //' IOSTAT='//CH_IOS
         RETURN
      ENDIF
      READ(80, R2LWSP)
      CLOSE(80)
!
!     TEST FOR MINIMAL REQUISITE INFORMATION.
      IF ( .NOT.(L_PRESENT(0).AND.                                      &
     &           L_PRESENT(6) ) ) THEN
         CMESSAGE='Longwave spectrum is deficient.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!
!
!     SET REDUCED DIMENSIONS, EITHER FROM THE SIZES OF THE FIXED ARRAYS
!     OR FROM THE ARRAYS READ IN.
!
      NPD_TYPE_LW=NPD_TYPE
      NPD_BAND_LW=MAX(N_BAND, 1)
      NPD_SPECIES_LW=MAX(N_ABSORB, 1)
      NPD_ALBEDO_PARM_LW=NPD_ALBEDO_PARM
      NPD_SCALE_FNC_LW=NPD_SCALE_FNC
      NPD_SCALE_VARIABLE_LW=NPD_SCALE_VARIABLE
      NPD_SURFACE_LW=NPD_SURFACE
      NPD_CONTINUUM_LW=NPD_CONTINUUM
      NPD_THERMAL_COEFF_LW=N_DEG_FIT+1
      NPD_CLOUD_PARAMETER_LW=NPD_CLOUD_PARAMETER
      NPD_AOD_WAVEL_LW=1
!
!
!     SEARCH THE SPECTRUM TO FIND MAXIMUM DIMENSIONS.
!
      NPD_EXCLUDE_LW=1
      IF (L_PRESENT(14)) THEN
         DO I=1, N_BAND
            NPD_EXCLUDE_LW=MAX(NPD_EXCLUDE_LW, N_BAND_EXCLUDE(I))
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
     &        ( (TYPE_ABSORB(I) == IP_CH4).AND.L_CH4 ).OR.              &
     &        ( (TYPE_ABSORB(I) == IP_N2O).AND.L_N2O ).OR.              &
     &        ( (TYPE_ABSORB(I) == IP_CFC11).AND.L_CFC11 ).OR.          &
     &        ( (TYPE_ABSORB(I) == IP_CFC12).AND.L_CFC12 ).OR.          &
     &        ( (TYPE_ABSORB(I) == IP_CFC113).AND.L_CFC113 ).OR.        &
     &        ( (TYPE_ABSORB(I) == IP_HCFC22).AND.L_HCFC22 ).OR.        &
     &        ( (TYPE_ABSORB(I) == IP_HFC125).AND.L_HFC125 ).OR.        &
     &        ( (TYPE_ABSORB(I) == IP_HFC134A).AND.L_HFC134A ) ) THEN
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
     &      , 'longwave spectral file.'
      ENDIF
!
      IF (.NOT.L_GAS_INCLUDED(IP_CO2)) THEN
         WRITE(IU_ERR, '(/A, /A)')                                      &
     &      '*** WARNING: Carbon dioxide is not included in the '       &
     &      , 'longwave spectral file.'
      ENDIF
!
      IF (.NOT.L_GAS_INCLUDED(IP_O3)) THEN
         WRITE(IU_ERR, '(/A, /A)')                                      &
     &      '*** WARNING: Ozone is not included in the '                &
     &      , 'longwave spectral file.'
      ENDIF
!
      IF ((.NOT.L_GAS_INCLUDED(IP_CH4)).AND.L_CH4) THEN
         WRITE(IU_ERR, '(/A, /A)')                                      &
     &      '*** ERROR: Methane is not included in the longwave '       &
     &      , 'spectral file, but was requested in the run.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      IF ((.NOT.L_GAS_INCLUDED(IP_N2O)).AND.L_N2O) THEN
         WRITE(IU_ERR, '(/A, /A)')                                      &
     &      '*** ERROR: Nitrous oxide is not included in the longwave ' &
     &      , 'spectral file, but was requested in the run.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      IF ((.NOT.L_GAS_INCLUDED(IP_CFC11)).AND.L_CFC11) THEN
         WRITE(IU_ERR, '(/A, /A)')                                      &
     &      '*** ERROR: CFC11 is not included in the longwave '         &
     &      , 'spectral file, but was requested in the run.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      IF ((.NOT.L_GAS_INCLUDED(IP_CFC12)).AND.L_CFC12) THEN
         WRITE(IU_ERR, '(/A, /A)')                                      &
     &      '*** ERROR: CFC12 is not included in the longwave '         &
     &      , 'spectral file, but was requested in the run.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      IF ((.NOT.L_GAS_INCLUDED(IP_CFC113)).AND.L_CFC113) THEN
         WRITE(IU_ERR, '(/A, /A)')                                      &
     &      '*** ERROR: CFC113 is not included in the longwave '        &
     &      , 'spectral file, but was requested in the run.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      IF ((.NOT.L_GAS_INCLUDED(IP_HCFC22)).AND.L_HCFC22) THEN
         WRITE(IU_ERR, '(/A, /A)')                                      &
     &      '*** ERROR: HCFC22 is not included in the longwave '        &
     &      , 'spectral file, but was requested in the run.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      IF ((.NOT.L_GAS_INCLUDED(IP_HFC125)).AND.L_HFC125) THEN
         WRITE(IU_ERR, '(/A, /A)')                                      &
     &      '*** ERROR: HFC125 is not included in the longwave '        &
     &      , 'spectral file, but was requested in the run.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      IF ((.NOT.L_GAS_INCLUDED(IP_HFC134A)).AND.L_HFC134A) THEN
         WRITE(IU_ERR, '(/A, /A)')                                      &
     &      '*** ERROR: HFC134A is not included in the longwave '       &
     &      , 'spectral file, but was requested in the run.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!     Set an appropriate reduced dimension.
      NPD_SPECIES_LW=MAX(N_ABSORB_RETAIN, 1)
!
      NPD_ESFT_TERM_LW=1
      IF (L_PRESENT(5)) THEN
         DO I=1, N_BAND
            DO J=1, N_BAND_ABSORB(I)
               IF (L_RETAIN_ABSORB(INDEX_ABSORB(J, I)))                 &
     &            NPD_ESFT_TERM_LW=MAX(NPD_ESFT_TERM_LW                 &
     &            , I_BAND_ESFT(I, INDEX_ABSORB(J, I)))
            ENDDO
         ENDDO
      ENDIF
!
      NPD_DROP_TYPE_LW=1
      IF (L_PRESENT(10)) THEN
         DO I=1, NPD_DROP_TYPE
            IF (L_DROP_TYPE(I)) THEN
               NPD_DROP_TYPE_LW=MAX(NPD_DROP_TYPE_LW, I)
            ENDIF
         ENDDO
      ENDIF
!
      NPD_ICE_TYPE_LW=1
      IF (L_PRESENT(12)) THEN
         DO I=1, NPD_ICE_TYPE
            IF (L_ICE_TYPE(I)) THEN
               NPD_ICE_TYPE_LW=MAX(NPD_ICE_TYPE_LW, I)
            ENDIF
         ENDDO
      ENDIF
!
!
!
!     Aerosols must be treated carefully to allow for various
!     different combinations without requiring the spectral file
!     to be too constrained. Only those required will be retained.
!
!     Basic initialization to safe values.
      NPD_HUMIDITIES_LW=1
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
               CMESSAGE='The LW Spectral file lacks some '              &
     &            //'climatological aerosols.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='LW Spectral file contains no aerosol data.'
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
               CMESSAGE='The LW Spectral file lacks some '              &
     &            //'mineral dust aerosol.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='LW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!



!     Check the spectral file for soot aerosols.
      IF (L_USE_SOOT_DIRECT .OR. L_USE_ARCLBLCK) THEN
         IF (L_PRESENT(11)) THEN
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
               IF ((TYPE_AEROSOL(I) == IP_FRESH_SOOT) .OR.              &
     &             (TYPE_AEROSOL(I) == IP_AGED_SOOT)) THEN
                  N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                  INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                  N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF
            ENDDO
!
            IF (N_AEROSOL_FOUND /= 2) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The LW Spectral file lacks some '              &
     &            //'soot aerosol data.'
               RETURN
!
            ENDIF
!
         ELSE
!
!
            IERR=I_ERR_FATAL
            CMESSAGE='LW Spectral file contains no soot data.'
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
               CMESSAGE='The LW Spectral file lacks some '              &
     &            //'biomass aerosol.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='LW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
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
               CMESSAGE='The LW Spectral file lacks some '              &
     &            //'sulphate aerosols.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='LW Spectral file contains no aerosol data.'
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
               CMESSAGE='The LW Spectral file lacks some '              &
     &            //'sea-salt aerosol.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='LW Spectral file contains no aerosol data.'
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
               CMESSAGE='The LW Spectral file lacks some '              &
     &            //'biogenic aerosol.'
               RETURN
!
            ENDIF

          ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='LW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!     Check the spectral file for fossil-fuel organic carbon aerosol.
!     (Only required for the direct effect).
!
      IF (L_USE_OCFF_DIRECT .OR. L_USE_ARCLOCFF) THEN
!
         IF (L_PRESENT(11)) THEN ! aerosol block present in spec file
!
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( (TYPE_AEROSOL(I) == IP_OCFF_FRESH).OR.              &
     &              (TYPE_AEROSOL(I) == IP_OCFF_AGED) ) THEN
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
               CMESSAGE='The LW Spectral file lacks some '              &
     &            //'fossil fuel org.carb. aerosol.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='LW Spectral file contains no aerosol data.'
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
               CMESSAGE='The LW Spectral file lacks some '              &
     &            //'delta aerosol.'
               RETURN
!
            ENDIF

          ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='LW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!     Set an appropriate reduced dimension.
      NPD_AEROSOL_SPECIES_LW=MAX(N_AEROSOL_RETAIN, 1)
!
!     Set the allowed number of humidities from the number of
!     retained aerosols.
!
      IF (L_PRESENT(11)) THEN
         DO I=1, N_AEROSOL_RETAIN
            IF (I_AEROSOL_PARAMETRIZATION(INDEX_AEROSOL_RETAIN(I)) ==   &
     &         IP_AEROSOL_PARAM_MOIST) THEN
               NPD_HUMIDITIES_LW=MAX(NPD_HUMIDITIES_LW                  &
     &            , NHUMIDITY(INDEX_AEROSOL_RETAIN(I)))
            ENDIF
         ENDDO
      ENDIF
!
!     Check that block 15 is present if the aerosol optical depth
!     was requested.
      IF (L_USE_AOD) THEN
        IF(.NOT. L_PRESENT(15)) THEN
          IERR = I_ERR_FATAL
          CMESSAGE='Block 15 needed in the LW spectral file.'
          RETURN
        ENDIF
       ! Check that the number of wavelengths in the spectral
       ! file is not larger than the number set in the include
       ! file MXSIZE3A.
        IF(N_AOD_WAVEL  >   NPD_AOD_WAVEL) THEN
          IERR = I_ERR_FATAL
          CMESSAGE='Increase NPD_AOD_WAVEL in MXSIZE3A.'
          RETURN
        ENDIF
        NPD_AOD_WAVEL_LW = NPD_AOD_WAVEL
      ENDIF
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
     &   , L_DOPPLER_PRESENT, DOPPLER_CORRECTION, L_USE_AOD             &
     &   , N_AOD_WAVEL, AOD_ABSORPTION, AOD_SCATTERING, I_AOD_TYPE      &
!                       Reduced Spectral Array
!     ------------------------------------------------------------------
!     ARGUMENT LIST OF LW SPECTRAL DATA.
!     (NOTE: LWSPDC3A, LWSPCM3A AND LWSARG3A MUST BE CONSISTENT)
!
     &   , NPD_TYPE_LW, NPD_BAND_LW, NPD_EXCLUDE_LW                     &
     &   , NPD_SPECIES_LW, NPD_ESFT_TERM_LW, NPD_SCALE_FNC_LW           &
     &   , NPD_SCALE_VARIABLE_LW, NPD_THERMAL_COEFF_LW                  &
     &   , NPD_SURFACE_LW, NPD_ALBEDO_PARM_LW, NPD_CONTINUUM_LW         &
     &   , NPD_DROP_TYPE_LW, NPD_ICE_TYPE_LW, NPD_CLOUD_PARAMETER_LW    &
     &   , NPD_AEROSOL_SPECIES_LW, NPD_HUMIDITIES_LW, NPD_AOD_WAVEL_LW  &
     &   , L_PRESENT_LW, N_BAND_LW, WAVE_LENGTH_SHORT_LW                &
     &   , WAVE_LENGTH_LONG_LW, N_BAND_EXCLUDE_LW, INDEX_EXCLUDE_LW     &
     &   , SOLAR_FLUX_BAND_LW, RAYLEIGH_COEFFICIENT_LW, N_ABSORB_LW     &
     &   , N_BAND_ABSORB_LW, INDEX_ABSORB_LW, TYPE_ABSORB_LW            &
     &   , I_BAND_ESFT_LW, I_SCALE_ESFT_LW, I_SCALE_FNC_LW, K_ESFT_LW   &
     &   , W_ESFT_LW, SCALE_VECTOR_LW, P_REFERENCE_LW, T_REFERENCE_LW   &
     &   , N_DEG_FIT_LW, THERMAL_COEFFICIENT_LW, T_REF_PLANCK_LW        &
     &   , I_SPEC_SURFACE_LW, N_DIR_ALBEDO_FIT_LW, L_SURFACE_LW         &
     &   , SURFACE_ALBEDO_LW, DIRECT_ALBEDO_PARM_LW                     &
     &   , EMISSIVITY_GROUND_LW, N_BAND_CONTINUUM_LW                    &
     &   , INDEX_CONTINUUM_LW, INDEX_WATER_LW, I_SCALE_FNC_CONT_LW      &
     &   , K_CONTINUUM_LW, SCALE_CONTINUUM_LW, P_REF_CONTINUUM_LW       &
     &   , T_REF_CONTINUUM_LW, I_DROP_PARAMETRIZATION_LW                &
     &   , L_DROP_TYPE_LW, DROP_PARAMETER_LIST_LW                       &
     &   , DROP_PARM_MIN_DIM_LW, DROP_PARM_MAX_DIM_LW                   &
     &   , N_AEROSOL_LW, TYPE_AEROSOL_LW, I_AEROSOL_PARAMETRIZATION_LW  &
     &   , NHUMIDITY_LW, HUMIDITIES_LW, L_AEROSOL_SPECIES_LW            &
     &   , AEROSOL_ABSORPTION_LW, AEROSOL_SCATTERING_LW                 &
     &   , AEROSOL_ASYMMETRY_LW, I_ICE_PARAMETRIZATION_LW               &
     &   , L_ICE_TYPE_LW, ICE_PARAMETER_LIST_LW                         &
     &   , ICE_PARM_MIN_DIM_LW, ICE_PARM_MAX_DIM_LW                     &
     &   , L_DOPPLER_PRESENT_LW, DOPPLER_CORRECTION_LW                  &
     &   , N_AOD_WAVEL_LW, AOD_ABSORPTION_LW, AOD_SCATTERING_LW         &
     &   , I_AOD_TYPE_LW                                                &
!
!     ------------------------------------------------------------------
     &   )
!
!
!
      RETURN
      END SUBROUTINE R2_LW_SPECIN
