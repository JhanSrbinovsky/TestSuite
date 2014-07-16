#if defined(A70_1A) || defined(A70_1B)
#if defined(A01_3A) || defined(A02_3A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to transfer spectrum to reduced array.
!
! Purpose:
!       Spectral data from the large dynamically allocated array
!       are transferred to the reduced array.
!
! Method:
!       Elements are copied across.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.4             03-09-97                Coding changes
!                                               associated with the
!                                               removal of pointers
!                                               into the spectral data.
!                                               Capability to select
!                                               aerosols added.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Coding to allow
!                                               selection of gases
!                                               from the spectral
!                                               file.
!                                               (J. M. Edwards)
!       5.1             04-04-00                Correct the limits
!                                               of a loop to
!                                               initialize the
!                                               full array.
!                                               (J. M. Edwards)
!       5.2             21-03-01                Provide default
!                                               initilizations for
!                                               logical variables to
!                                               deal with cases when
!                                               particular blocks are
!                                               absent from the
!                                               spectral file.
!                                               (J. M. Edwards)
!       5.5             24-02-03                Code for aggregate
!                                               parametrization of
!                                               ice crystals included.
!                                               (J. M. Edwards)
!       6.2             15-11-05                Added code for aerosol
!                                               optical depth (block
!                                               15 of the LW spectral
!                                               file).
!                                               (N. Bellouin)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_COMPRESS_SPECTRUM(                                  &
!                       Original Spectrum
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
     &   , NPDR_TYPE, NPDR_BAND, NPDR_EXCLUDE                           &
     &   , NPDR_SPECIES, NPDR_ESFT_TERM, NPDR_SCALE_FNC                 &
     &   , NPDR_SCALE_VARIABLE, NPDR_THERMAL_COEFF                      &
     &   , NPDR_SURFACE, NPDR_ALBEDO_PARM                               &
     &   , NPDR_CONTINUUM, NPDR_DROP_TYPE, NPDR_ICE_TYPE                &
     &   , NPDR_CLOUD_PARAMETER, NPDR_AEROSOL_SPECIES                   &
     &   , NPDR_HUMIDITIES, NPDR_AOD_WAVEL                              &
     &   , L_PRESENT_RD                                                 &
     &   , N_BAND_RD, WAVE_LENGTH_SHORT_RD , WAVE_LENGTH_LONG_RD        &
     &   , N_BAND_EXCLUDE_RD, INDEX_EXCLUDE_RD                          &
     &   , SOLAR_FLUX_BAND_RD, RAYLEIGH_COEFFICIENT_RD                  &
     &   , N_ABSORB_RD, N_BAND_ABSORB_RD, INDEX_ABSORB_RD               &
     &   , TYPE_ABSORB_RD                                               &
     &   , I_BAND_ESFT_RD, I_SCALE_ESFT_RD, I_SCALE_FNC_RD              &
     &   , K_ESFT_RD, W_ESFT_RD, SCALE_VECTOR_RD                        &
     &   , P_REFERENCE_RD, T_REFERENCE_RD                               &
     &   , N_DEG_FIT_RD, THERMAL_COEFFICIENT_RD, T_REF_PLANCK_RD        &
     &   , I_SPEC_SURFACE_RD, N_DIR_ALBEDO_FIT_RD                       &
     &   , L_SURFACE_RD, SURFACE_ALBEDO_RD, DIRECT_ALBEDO_PARM_RD       &
     &   , EMISSIVITY_GROUND_RD                                         &
     &   , N_BAND_CONTINUUM_RD, INDEX_CONTINUUM_RD, INDEX_WATER_RD      &
     &   , I_SCALE_FNC_CONT_RD, K_CONTINUUM_RD, SCALE_CONTINUUM_RD      &
     &   , P_REF_CONTINUUM_RD, T_REF_CONTINUUM_RD                       &
     &   , I_DROP_PARAMETRIZATION_RD, L_DROP_TYPE_RD                    &
     &   , DROP_PARAMETER_LIST_RD                                       &
     &   , DROP_PARM_MIN_DIM_RD, DROP_PARM_MAX_DIM_RD                   &
     &   , N_AEROSOL_RD, TYPE_AEROSOL_RD, I_AEROSOL_PARAMETRIZATION_RD  &
     &   , NHUMIDITY_RD, HUMIDITIES_RD                                  &
     &   , L_AEROSOL_SPECIES_RD, AEROSOL_ABSORPTION_RD                  &
     &   , AEROSOL_SCATTERING_RD, AEROSOL_ASYMMETRY_RD                  &
     &   , I_ICE_PARAMETRIZATION_RD, L_ICE_TYPE_RD                      &
     &   , ICE_PARAMETER_LIST_RD                                        &
     &   , ICE_PARM_MIN_DIM_RD, ICE_PARM_MAX_DIM_RD                     &
     &   , L_DOPPLER_PRESENT_RD, DOPPLER_CORRECTION_RD                  &
     &   , N_AOD_WAVEL_RD, AOD_ABSORPTION_RD, AOD_SCATTERING_RD         &
     &   , I_AOD_TYPE_RD                                                &
     &   )
!
!
      IMPLICIT NONE
!
!
!
!     ------------------------------------------------------------------
!     DECLARATION OF INITIAL SPECTRUM.
!     ------------------------------------------------------------------
#include "mxsize3a.h"
#include "spdec3a.h"
!
!     AUXILIARY VARIABLES USED TO SELECT PARTS OF THE INITIAL SPECTRUM
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_RETAIN_ABSORB(NPD_SPECIES)
!             FLAGS FOR THE RETENTION OF GASES IN THE SPECTRAL FILE
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_ABSORB_RETAIN                                              &
!             NUMBER OF ABSORBERS TO BE RETAINED
     &   , INDEX_ABSORB_RETAIN(NPD_SPECIES)                             &
!             INDICES OF ABSORBERS TO BE RETAINED
     &   , COMPRESSED_INDEX(NPD_SPECIES)                                &
!             MAPPING FROM OLD TO NEW INDICES OF ABSORBERS
     &   , N_AEROSOL_RETAIN                                             &
!             NUMBER OF AEROSOLS IN THE INITIAL SPECTRUM TO BE USED
!             IN THE CALCULATION
     &   , INDEX_AEROSOL_RETAIN(NPD_AEROSOL_SPECIES)
!             INDICES OF THE RETAINED AEROSOLS
!
!
!     ------------------------------------------------------------------
!     DECLARATION OF REDUCED SPECTRUM.
!     ------------------------------------------------------------------
!
!     DIMENSIONS OF REDUCED ARRAY:
!
      INTEGER                                                           &
     &     NPDR_BAND                                                    &
!             NUMBER OF SPECTRAL BANDS
     &   , NPDR_EXCLUDE                                                 &
!             NUMER OF EXCLUDED BANDS
     &   , NPDR_ESFT_TERM                                               &
!             NUMBER OF ESFT TERMS
     &   , NPDR_TYPE                                                    &
!             NUMBER OF DATA TYPES
     &   , NPDR_SPECIES                                                 &
!             NUMBER OF GASEOUS SPECIES
     &   , NPDR_SCALE_FNC                                               &
!             NUMBER OF SCALING FUNCTIONS
     &   , NPDR_SCALE_VARIABLE                                          &
!             NUMBER OF SCALING VARIABLES
     &   , NPDR_SURFACE                                                 &
!             NUMBER OF SURFACE TYPES
     &   , NPDR_ALBEDO_PARM                                             &
!             NUMBER OF ALBEDO PARAMETERS
     &   , NPDR_CONTINUUM                                               &
!             NUMBER OF CONTINUA
     &   , NPDR_DROP_TYPE                                               &
!             NUMBER OF DROP TYPES
     &   , NPDR_ICE_TYPE                                                &
!             NUMBER OF ICE CRYSTAL TYPES
     &   , NPDR_AEROSOL_SPECIES                                         &
!             NUMBER OF AEROSOL SPECIES
     &   , NPDR_THERMAL_COEFF                                           &
!             NUMBER OF THERMAL COEFFICIENTS
     &   , NPDR_CLOUD_PARAMETER                                         &
!             MAX NUMBER OF CLOUD PARAMETERS
     &   , NPDR_HUMIDITIES                                              &
!             MAXIMUM NUMBER OF HUMIDITIES
     &   , NPDR_AOD_WAVEL
!             NUMBER OF WAVELENGTHS FOR AEROSOL OPTICAL DEPTH
!
!
!     GENERAL FIELDS:
!
      LOGICAL                                                           &
     &     L_PRESENT_RD(0: NPDR_TYPE)
!             FLAG FOR TYPES OF DATA PRESENT
!
!
!
!     PROPERTIES OF THE SPECTRAL BANDS:
!
      INTEGER                                                           &
     &     N_BAND_RD
!             NUMBER OF SPECTRAL BANDS
!
      REAL                                                              &
     &     WAVE_LENGTH_SHORT_RD(NPDR_BAND)                              &
!             SHORTER WAVELENGTH LIMITS
     &   , WAVE_LENGTH_LONG_RD(NPDR_BAND)
!             LONGER WAVELENGTH LIMITS
!
!
!
!     EXCLUSION OF SPECIFIC BANDS FROM PARTS OF THE SPECTRUM:
      INTEGER                                                           &
     &     N_BAND_EXCLUDE_RD(NPDR_BAND)                                 &
!             NUMBER OF EXCLUDED BANDS WITHIN EACH SPECTRAL BAND
     &   , INDEX_EXCLUDE_RD(NPDR_EXCLUDE, NPDR_BAND)
!             INDICES OF EXCLUDED BANDS
!
!
!
!     FIELDS FOR THE SOLAR FLUX:
!
      REAL                                                              &
     &     SOLAR_FLUX_BAND_RD(NPDR_BAND)
!             FRACTION OF THE INCIDENT SOLAR FLUX IN EACH BAND
!
!
!
!     FIELDS FOR RAYLEIGH SCATTERING:
!
      REAL                                                              &
     &     RAYLEIGH_COEFFICIENT_RD(NPDR_BAND)
!             RAYLEIGH COEFFICIENTS
!
!
!
!     FIELDS FOR GASEOUS ABSORPTION:
!
      INTEGER                                                           &
     &     N_ABSORB_RD                                                  &
!             NUMBER OF ABSORBERS
     &   , N_BAND_ABSORB_RD(NPDR_BAND)                                  &
!             NUMBER OF ABSORBERS IN EACH BAND
     &   , INDEX_ABSORB_RD(NPDR_SPECIES, NPDR_BAND)                     &
!             LIST OF ABSORBERS IN EACH BAND
     &   , TYPE_ABSORB_RD(NPDR_SPECIES)                                 &
!             TYPES OF EACH GAS IN THE SPECTRAL FILE
     &   , I_BAND_ESFT_RD(NPDR_BAND, NPDR_SPECIES)                      &
!             NUMBER OF ESFT TERMS IN BAND FOR EACH GAS
     &   , I_SCALE_ESFT_RD(NPDR_BAND, NPDR_SPECIES)                     &
!             TYPE OF ESFT SCALING
     &   , I_SCALE_FNC_RD(NPDR_BAND, NPDR_SPECIES)
!             TYPE OF SCALING FUNCTION
!
      REAL                                                              &
     &     K_ESFT_RD(NPDR_ESFT_TERM, NPDR_BAND, NPDR_SPECIES)           &
!             ESFT EXPONENTS
     &   , W_ESFT_RD(NPDR_ESFT_TERM, NPDR_BAND, NPDR_SPECIES)           &
!             ESFT WEIGHTS
     &   , SCALE_VECTOR_RD(NPDR_SCALE_VARIABLE, NPDR_ESFT_TERM          &
     &        , NPDR_BAND, NPDR_SPECIES)                                &
!             SCALING PARAMETERS FOR EACH ABSORBER AND TERM
     &   , P_REFERENCE_RD(NPDR_SPECIES, NPDR_BAND)                      &
!             REFERENCE PRESSURE FOR SCALING FUNCTION
     &   , T_REFERENCE_RD(NPDR_SPECIES, NPDR_BAND)
!             REFERENCE TEMPERATURE FOR SCALING FUNCTION
!
!
!
!     REPRESENTATION OF THE PLANCKIAN:
!
      INTEGER                                                           &
     &     N_DEG_FIT_RD
!             DEGREE OF THERMAL POLYNOMIAL
!
      REAL                                                              &
     &     THERMAL_COEFFICIENT_RD(0: NPDR_THERMAL_COEFF-1, NPDR_BAND)   &
!             COEFFICIENTS IN POLYNOMIAL FIT TO SOURCE FUNCTION
     &   , T_REF_PLANCK_RD
!             PLANCKIAN REFERENCE TEMPERATURE
!
!
!
!     SURFACE PROPERTIES:
!
      INTEGER                                                           &
     &     I_SPEC_SURFACE_RD(NPDR_SURFACE)                              &
!             METHOD OF SPECIFYING PROPERTIES OF SURFACE
     &   , N_DIR_ALBEDO_FIT_RD(NPDR_SURFACE)
!             NUMBER OF PARAMETERS FITTING THE DIRECT ALBEDO
!
      LOGICAL                                                           &
     &     L_SURFACE_RD(NPDR_SURFACE)
!             SURFACE TYPES INCLUDED
!
      REAL                                                              &
     &     SURFACE_ALBEDO_RD(NPDR_BAND, NPDR_SURFACE)                   &
!             SURFACE ALBEDOS
     &   , DIRECT_ALBEDO_PARM_RD(0: NPDR_ALBEDO_PARM                    &
     &        , NPD_BAND, NPD_SURFACE)                                  &
!             COEFFICIENTS FOR FITTING DIRECT ALBEDO
     &   , EMISSIVITY_GROUND_RD(NPDR_BAND, NPDR_SURFACE)
!             SURFACE EMISSIVITIES
!
!
!
!     FIELDS FOR CONTINUA:
!
      INTEGER                                                           &
     &     N_BAND_CONTINUUM_RD(NPDR_BAND)                               &
!             NUMBER OF CONTINUA IN EACH BAND
     &   , INDEX_CONTINUUM_RD(NPDR_BAND, NPDR_CONTINUUM)                &
!             LIST OF CONTINUA IN EACH BAND
     &   , INDEX_WATER_RD                                               &
!             INDEX OF WATER VAPOUR
     &   , I_SCALE_FNC_CONT_RD(NPDR_BAND, NPDR_CONTINUUM)
!             TYPE OF SCALING FUNCTION FOR CONTINUUM
!
      REAL                                                              &
     &     K_CONTINUUM_RD(NPDR_BAND, NPDR_CONTINUUM)                    &
!             GREY EXTINCTION COEFFICIENTS FOR CONTINUUM
     &   , SCALE_CONTINUUM_RD(NPDR_SCALE_VARIABLE                       &
     &        , NPDR_BAND, NPDR_CONTINUUM)                              &
!             SCALING PARAMETERS FOR CONTINUUM
     &   , P_REF_CONTINUUM_RD(NPDR_CONTINUUM, NPDR_BAND)                &
!             REFERENCE PRESSURE FOR SCALING OF CONTINUUM
     &   , T_REF_CONTINUUM_RD(NPDR_CONTINUUM, NPDR_BAND)
!             REFERENCE TEMPERATURE FOR SCALING OF CONTINUUM
!
!
!
!     FIELDS FOR WATER DROPLETS:
!
      INTEGER                                                           &
     &     I_DROP_PARAMETRIZATION_RD(NPDR_DROP_TYPE)
!             PARAMETRIZATION TYPE OF DROPLETS
!
      LOGICAL                                                           &
     &     L_DROP_TYPE_RD(NPDR_DROP_TYPE)
!             TYPES OF DROPLET PRESENT
!
      REAL                                                              &
     &     DROP_PARAMETER_LIST_RD(NPDR_CLOUD_PARAMETER                  &
     &        , NPDR_BAND, NPDR_DROP_TYPE)                              &
!             PARAMETERS USED TO FIT OPTICAL PROPERTIES OF CLOUDS
     &   , DROP_PARM_MIN_DIM_RD(NPDR_DROP_TYPE)                         &
!             MINIMUM SIZE OF DROPLET PERMITTED IN THE PARAMETRIZATION
     &   , DROP_PARM_MAX_DIM_RD(NPDR_DROP_TYPE)
!             MAXIMUM SIZE OF DROPLET PERMITTED IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR AEROSOLS:
!
      INTEGER                                                           &
     &     N_AEROSOL_RD                                                 &
!             NUMBER OF SPECIES OF AEROSOL
     &   , TYPE_AEROSOL_RD(NPDR_AEROSOL_SPECIES)                        &
!             TYPES OF AEROSOLS
     &   , I_AEROSOL_PARAMETRIZATION_RD(NPDR_AEROSOL_SPECIES)           &
!             PARAMETRIZATION OF AEROSOLS
     &   , NHUMIDITY_RD(NPDR_AEROSOL_SPECIES)
!             NUMBERS OF HUMIDITIES
!
      LOGICAL                                                           &
     &     L_AEROSOL_SPECIES_RD(NPDR_AEROSOL_SPECIES)
!             AEROSOL SPECIES INCLUDED
!
      REAL                                                              &
     &     AEROSOL_ABSORPTION_RD(NPDR_HUMIDITIES, NPDR_AEROSOL_SPECIES  &
     &        , NPDR_BAND)                                              &
!             ABSORPTION BY AEROSOLS
     &   , AEROSOL_SCATTERING_RD(NPDR_HUMIDITIES, NPDR_AEROSOL_SPECIES  &
     &        , NPDR_BAND)                                              &
!             SCATTERING BY AEROSOLS
     &   , AEROSOL_ASYMMETRY_RD(NPDR_HUMIDITIES, NPDR_AEROSOL_SPECIES   &
     &        , NPDR_BAND)                                              &
!             ASYMMETRY OF AEROSOLS
     &   , HUMIDITIES_RD(NPDR_HUMIDITIES, NPDR_AEROSOL_SPECIES)
!             HUMIDITIES FOR COMPONENTS
!
!
!
!     FIELDS FOR ICE CRYSTALS:
!
      INTEGER                                                           &
     &     I_ICE_PARAMETRIZATION_RD(NPDR_ICE_TYPE)
!             TYPES OF PARAMETRIZATION OF ICE CRYSTALS
!
      LOGICAL                                                           &
     &     L_ICE_TYPE_RD(NPDR_ICE_TYPE)
!             TYPES OF ICE CRYSTAL PRESENT
!
      REAL                                                              &
     &     ICE_PARAMETER_LIST_RD(NPDR_CLOUD_PARAMETER                   &
     &        , NPDR_BAND, NPDR_ICE_TYPE)                               &
!             PARAMETERS USED TO FIT SINGLE SCATTERING OF ICE CRYSTALS
     &   , ICE_PARM_MIN_DIM_RD(NPDR_ICE_TYPE)                           &
!             MINIMUM SIZE OF ICE CRYSTAL PERMITTED
!             IN THE PARAMETRIZATION
     &   , ICE_PARM_MAX_DIM_RD(NPDR_ICE_TYPE)
!             MAXIMUM SIZE OF ICE CRYSTAL PERMITTED
!             IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR DOPPLER BROADENING:
!
      LOGICAL                                                           &
     &     L_DOPPLER_PRESENT_RD(NPDR_SPECIES)
!             FLAG FOR DOPPLER BROADENING FOR EACH SPECIES
!
      REAL                                                              &
     &     DOPPLER_CORRECTION_RD(NPDR_SPECIES)
!             DOPPLER CORRECTION TERMS
!
!     FIELDS FOR THE AEROSOL OPTICAL DEPTH
!
      LOGICAL L_USE_AOD
!             At least one of the aerosol optical depth diagnostics
!             is requested. Otherwise, block 15 of the LW spectral
!             file does not need to be read in.

      INTEGER N_AOD_WAVEL_RD
      REAL AOD_ABSORPTION_RD(NPDR_HUMIDITIES,                           &
     &                       NPDR_AEROSOL_SPECIES,                      &
     &                       NPDR_AOD_WAVEL),                           &
     &     AOD_SCATTERING_RD(NPDR_HUMIDITIES,                           &
     &                       NPDR_AEROSOL_SPECIES,                      &
     &                       NPDR_AOD_WAVEL)
      INTEGER I_AOD_TYPE_RD(NPDR_AEROSOL_SPECIES)
!
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , J                                                            &
!             LOOP VARIABLE
     &   , K                                                            &
!             LOOP VARIABLE
     &   , L                                                            &
!             LOOP VARIABLE
     &   , N_PARAMETER                                                  &
!             NUMBER OF PARAMETERS IN SCHEME.
     &   , I_SPECIES                                                    &
!             SPECIES OF GAS
     &   , I_CONTINUUM                                                  &
!             TYPE OF CONTINUUM
     &   , I_INITIAL
!             INDEXING NUMBER IN INITIAL SPECTRAL FILE
!
!
#include "sclfnc3a.h"
#include "wclprm3a.h"
#include "iclprm3a.h"
#include "aerprm3a.h"
#include "sclfnd3a.h"
!
!
!
!
!
!     Initialize all blocks of the compressed spectrum to .FALSE.
!     and do the same for all flags within the spectrum.
      DO I=0, NPDR_TYPE
         L_PRESENT_RD(I)=.FALSE.
      ENDDO
      DO I=1, NPDR_DROP_TYPE
        L_DROP_TYPE_RD(I)=.FALSE.
      ENDDO
      DO I=1, NPDR_AEROSOL_SPECIES
        L_AEROSOL_SPECIES_RD(I)=.FALSE.
      ENDDO
      DO I=1, NPDR_ICE_TYPE
        L_ICE_TYPE_RD(I)=.FALSE.
      ENDDO
      DO I=1, NPDR_SPECIES
        L_DOPPLER_PRESENT_RD(I)=.FALSE.
      ENDDO
!
!
!     PROCEED THROUGH EACH BLOCK OF THE SPECTRAL FILE TRANSFERRING
!     THE DATA FROM THE INPUT ARRAY TO THE REDUCED ARRAY.
!
!
!     BLOCK 0:
!
      IF (L_PRESENT(0)) THEN
         L_PRESENT_RD(0)=.TRUE.
         N_BAND_RD=N_BAND
         N_ABSORB_RD=N_ABSORB_RETAIN
         N_AEROSOL_RD=N_AEROSOL_RETAIN
         DO I=1, N_ABSORB_RETAIN
            TYPE_ABSORB_RD(I)=TYPE_ABSORB(INDEX_ABSORB_RETAIN(I))
         ENDDO
         DO I=1, N_AEROSOL_RETAIN
            TYPE_AEROSOL_RD(I)=TYPE_AEROSOL(INDEX_AEROSOL_RETAIN(I))
         ENDDO
      ENDIF
!
!     BLOCK 1:
      IF (L_PRESENT(1)) THEN
         L_PRESENT_RD(1)=.TRUE.
         DO I=1, N_BAND
            WAVE_LENGTH_SHORT_RD(I)=WAVE_LENGTH_SHORT(I)
            WAVE_LENGTH_LONG_RD(I)=WAVE_LENGTH_LONG(I)
         ENDDO
      ENDIF
!
!     BLOCK 2:
      IF (L_PRESENT(2)) THEN
         L_PRESENT_RD(2)=.TRUE.
         DO I=1, N_BAND
            SOLAR_FLUX_BAND_RD(I)=SOLAR_FLUX_BAND(I)
         ENDDO
      ENDIF
!
!     BLOCK 3:
      IF (L_PRESENT(3)) THEN
         L_PRESENT_RD(3)=.TRUE.
         DO I=1, N_BAND
            RAYLEIGH_COEFFICIENT_RD(I)=RAYLEIGH_COEFFICIENT(I)
         ENDDO
      ENDIF
!
!     BLOCK 4:
      IF (L_PRESENT(4)) THEN
         L_PRESENT_RD(4)=.TRUE.
         DO I=1, N_BAND
            N_BAND_ABSORB_RD(I)=0
            DO J=1, N_BAND_ABSORB(I)
               IF (L_RETAIN_ABSORB(INDEX_ABSORB(J, I))) THEN
                  N_BAND_ABSORB_RD(I)=N_BAND_ABSORB_RD(I)+1
                  INDEX_ABSORB_RD(N_BAND_ABSORB_RD(I), I)               &
     &               =COMPRESSED_INDEX(INDEX_ABSORB(J, I))
               ENDIF
            ENDDO
         ENDDO
      ENDIF
!
!     BLOCK 5:
      IF (L_PRESENT(5)) THEN
         L_PRESENT_RD(5)=.TRUE.
         DO I=1, N_BAND
            DO J=1, N_BAND_ABSORB_RD(I)
               I_SPECIES=INDEX_ABSORB_RD(J, I)
               I_INITIAL=INDEX_ABSORB_RETAIN(I_SPECIES)
               I_BAND_ESFT_RD(I, I_SPECIES)=I_BAND_ESFT(I, I_INITIAL)
               I_SCALE_ESFT_RD(I, I_SPECIES)=I_SCALE_ESFT(I, I_INITIAL)
               I_SCALE_FNC_RD(I, I_SPECIES)=I_SCALE_FNC(I, I_INITIAL)
               P_REFERENCE_RD(I_SPECIES, I)=P_REFERENCE(I_INITIAL, I)
               T_REFERENCE_RD(I_SPECIES, I)=T_REFERENCE(I_INITIAL, I)
               DO K=1, I_BAND_ESFT(I, I_INITIAL)
                  K_ESFT_RD(K, I, I_SPECIES)=K_ESFT(K, I, I_INITIAL)
                  W_ESFT_RD(K, I, I_SPECIES)=W_ESFT(K, I, I_INITIAL)
                  DO L=1, N_SCALE_VARIABLE(I_SCALE_FNC(I, I_INITIAL))
                     SCALE_VECTOR_RD(L, K, I, I_SPECIES)                &
     &                  =SCALE_VECTOR(L, K, I, I_INITIAL)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF
!
!     BLOCK 6:
      IF (L_PRESENT(6)) THEN
         L_PRESENT_RD(6)=.TRUE.
         N_DEG_FIT_RD=N_DEG_FIT
         T_REF_PLANCK_RD=T_REF_PLANCK
         DO I=1, N_BAND
            DO J=0, N_DEG_FIT
               THERMAL_COEFFICIENT_RD(J, I)=THERMAL_COEFFICIENT(J, I)
            ENDDO
         ENDDO
      ENDIF
!
!     BLOCK 7:
!
!     OMITTED SINCE SURFACE ALBEDOS ARE PROVIDED BY THE MODEL.
!
!     BLOCK 8:
      IF (L_PRESENT(8)) THEN
         L_PRESENT_RD(8)=.TRUE.
         DO I=1, N_BAND
            N_BAND_CONTINUUM_RD(I)=N_BAND_CONTINUUM(I)
            DO J=1, N_BAND_CONTINUUM(I)
               INDEX_CONTINUUM_RD(I, J)=INDEX_CONTINUUM(I, J)
            ENDDO
         ENDDO
!
         INDEX_WATER_RD=0
         DO I=1, N_ABSORB_RETAIN
            IF (INDEX_ABSORB_RETAIN(I) == INDEX_WATER) THEN
               INDEX_WATER_RD=I
            ENDIF
         ENDDO
!
      ENDIF
!
!     BLOCK 9:
      IF (L_PRESENT(9)) THEN
         L_PRESENT_RD(9)=.TRUE.
         DO I=1, N_BAND
            DO J=1, N_BAND_CONTINUUM(I)
               I_CONTINUUM=INDEX_CONTINUUM(I, J)
               I_SCALE_FNC_CONT_RD(I, I_CONTINUUM)                      &
     &            =I_SCALE_FNC_CONT(I, I_CONTINUUM)
               P_REF_CONTINUUM_RD(I_CONTINUUM, I)                       &
     &            =P_REF_CONTINUUM(I_CONTINUUM, I)
               T_REF_CONTINUUM_RD(I_CONTINUUM, I)                       &
     &            =T_REF_CONTINUUM(I_CONTINUUM, I)
               K_CONTINUUM_RD(I, I_CONTINUUM)                           &
     &            =K_CONTINUUM(I, I_CONTINUUM)
               DO L=1, N_SCALE_VARIABLE(I_SCALE_FNC_CONT                &
     &               (I, I_CONTINUUM))
                  SCALE_CONTINUUM_RD(L, I, I_CONTINUUM)                 &
     &               =SCALE_CONTINUUM(L, I, I_CONTINUUM)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
!
!     BLOCK 10:
      IF (L_PRESENT(10)) THEN
         L_PRESENT_RD(10)=.TRUE.
         DO I=1, NPDR_DROP_TYPE
            IF (L_DROP_TYPE(I)) THEN
               L_DROP_TYPE_RD(I)=.TRUE.
               I_DROP_PARAMETRIZATION_RD(I)=I_DROP_PARAMETRIZATION(I)
               DROP_PARM_MIN_DIM_RD(I)=DROP_PARM_MIN_DIM(I)
               DROP_PARM_MAX_DIM_RD(I)=DROP_PARM_MAX_DIM(I)
               IF (I_DROP_PARAMETRIZATION(I)                            &
     &             == IP_SLINGO_SCHRECKER) THEN
                  N_PARAMETER=6
               ELSE IF (I_DROP_PARAMETRIZATION(I)                       &
     &             == IP_ACKERMAN_STEPHENS) THEN
                  N_PARAMETER=9
               ELSE IF (I_DROP_PARAMETRIZATION(I)                       &
     &             == IP_DROP_PADE_2) THEN
                  N_PARAMETER=16
               ENDIF
!
               DO J=1, N_PARAMETER
                  DO K=1, N_BAND
                     DROP_PARAMETER_LIST_RD(J, K, I)                    &
     &                  =DROP_PARAMETER_LIST(J, K, I)
                  ENDDO
               ENDDO
            ELSE
               L_DROP_TYPE_RD(I)=.FALSE.
            ENDIF
         ENDDO
      ENDIF
!
!     BLOCK 11:
      IF (L_PRESENT(11)) THEN
         L_PRESENT_RD(11)=.TRUE.
         DO I=1, N_AEROSOL_RETAIN
            I_INITIAL=INDEX_AEROSOL_RETAIN(I)
            IF (L_AEROSOL_SPECIES(I_INITIAL)) THEN
               L_AEROSOL_SPECIES_RD(I)=.TRUE.
               I_AEROSOL_PARAMETRIZATION_RD(I)                          &
     &            =I_AEROSOL_PARAMETRIZATION(I_INITIAL)
               IF (I_AEROSOL_PARAMETRIZATION(I_INITIAL)                 &
     &             == IP_AEROSOL_PARAM_DRY) THEN
                  NHUMIDITY_RD(I)=0
                  DO K=1, N_BAND
                     AEROSOL_ABSORPTION_RD(1, I, K)                     &
     &                  =AEROSOL_ABSORPTION(1, I_INITIAL, K)
                     AEROSOL_SCATTERING_RD(1, I, K)                     &
     &                  =AEROSOL_SCATTERING(1, I_INITIAL, K)
                     AEROSOL_ASYMMETRY_RD(1, I, K)                      &
     &                  =AEROSOL_ASYMMETRY(1, I_INITIAL, K)
                  ENDDO
               ELSE IF (I_AEROSOL_PARAMETRIZATION(I_INITIAL)            &
     &             == IP_AEROSOL_PARAM_MOIST) THEN
                  INDEX_WATER_RD=INDEX_WATER
                  NHUMIDITY_RD(I)=NHUMIDITY(I_INITIAL)
                  DO J=1, NHUMIDITY(I_INITIAL)
                     HUMIDITIES_RD(J, I)=HUMIDITIES(J, I_INITIAL)
                     DO K=1, N_BAND
                        AEROSOL_ABSORPTION_RD(J, I, K)                  &
     &                     =AEROSOL_ABSORPTION(J, I_INITIAL, K)
                        AEROSOL_SCATTERING_RD(J, I, K)                  &
     &                     =AEROSOL_SCATTERING(J, I_INITIAL, K)
                        AEROSOL_ASYMMETRY_RD(J, I, K)                   &
     &                     =AEROSOL_ASYMMETRY(J, I_INITIAL, K)
                     ENDDO
                  ENDDO
               ENDIF
!
            ELSE
               L_AEROSOL_SPECIES_RD(I)=.FALSE.
            ENDIF
         ENDDO
      ENDIF
!
!     BLOCK 12:
      IF (L_PRESENT(12)) THEN
         L_PRESENT_RD(12)=.TRUE.
         DO I=1, NPDR_ICE_TYPE
            IF (L_ICE_TYPE(I)) THEN
               L_ICE_TYPE_RD(I)=.TRUE.
               ICE_PARM_MIN_DIM_RD(I)=ICE_PARM_MIN_DIM(I)
               ICE_PARM_MAX_DIM_RD(I)=ICE_PARM_MAX_DIM(I)
!
               I_ICE_PARAMETRIZATION_RD(I)=I_ICE_PARAMETRIZATION(I)
               IF (I_ICE_PARAMETRIZATION(I)                             &
     &             == IP_SLINGO_SCHRECKER_ICE) THEN
                  N_PARAMETER=6
               ELSE IF (I_ICE_PARAMETRIZATION(I)                        &
     &             == IP_SUN_SHINE_VN2_VIS) THEN
                  N_PARAMETER=6
               ELSE IF (I_ICE_PARAMETRIZATION(I)                        &
     &             == IP_SUN_SHINE_VN2_IR) THEN
                  N_PARAMETER=0
               ELSE IF (I_ICE_PARAMETRIZATION(I)                        &
     &             == IP_ICE_ADT) THEN
                  N_PARAMETER=30
               ELSE IF (I_ICE_PARAMETRIZATION(I)                        &
     &             == IP_ICE_AGG_DE) THEN
                  N_PARAMETER=14
               ENDIF
!
               DO J=1, N_PARAMETER
                  DO K=1, N_BAND
                     ICE_PARAMETER_LIST_RD(J, K, I)                     &
     &                  =ICE_PARAMETER_LIST(J, K, I)
                  ENDDO
               ENDDO
            ELSE
               L_ICE_TYPE_RD(I)=.FALSE.
            ENDIF
         ENDDO
      ENDIF
!
!     BLOCK 13:
      IF (L_PRESENT(13)) THEN
         L_PRESENT_RD(13)=.TRUE.
         DO I=1, N_ABSORB
            IF (L_RETAIN_ABSORB(I)) THEN
               L_DOPPLER_PRESENT_RD(COMPRESSED_INDEX(I))                &
     &            =L_DOPPLER_PRESENT(I)
               IF (L_DOPPLER_PRESENT(I))                                &
     &            DOPPLER_CORRECTION_RD(COMPRESSED_INDEX(I))            &
     &               =DOPPLER_CORRECTION(I)
            ENDIF
         ENDDO
      ENDIF
!
!
!     BLOCK 14:
      IF (L_PRESENT(14)) THEN
         L_PRESENT_RD(14)=.TRUE.
         DO I=1, N_BAND
             N_BAND_EXCLUDE_RD(I)=N_BAND_EXCLUDE(I)
             DO J=1, N_BAND_EXCLUDE(I)
                INDEX_EXCLUDE_RD(J, I)=INDEX_EXCLUDE(J, I)
             ENDDO
         ENDDO
      ENDIF
!
!
!     BLOCK 15 (we rely on the work done for block 11)
      IF (L_PRESENT(15).AND.L_USE_AOD) THEN
         L_PRESENT_RD(15)=.TRUE.
         N_AOD_WAVEL_RD = N_AOD_WAVEL

         DO I=1, N_AEROSOL_RETAIN
            I_INITIAL=INDEX_AEROSOL_RETAIN(I)
            IF (L_AEROSOL_SPECIES(I_INITIAL)) THEN
               I_AOD_TYPE_RD(I) = I_AOD_TYPE(I_INITIAL)
               IF (I_AEROSOL_PARAMETRIZATION(I_INITIAL)                 &
     &              ==  IP_AEROSOL_PARAM_DRY) THEN
                  DO K=1, N_AOD_WAVEL
                    AOD_ABSORPTION_RD(1, I, K) =                        &
     &               AOD_ABSORPTION(1, I_INITIAL, K)
                    AOD_SCATTERING_RD(1, I, K) =                        &
     &               AOD_SCATTERING(1, I_INITIAL, K)
                  ENDDO
               ELSE IF (I_AEROSOL_PARAMETRIZATION(I_INITIAL)            &
     &                   ==  IP_AEROSOL_PARAM_MOIST) THEN
                  DO J=1, NHUMIDITY(I_INITIAL)
                    DO K=1, N_AOD_WAVEL
                      AOD_ABSORPTION_RD(J, I, K) =                      &
     &                 AOD_ABSORPTION(J, I_INITIAL, K)
                      AOD_SCATTERING_RD(J, I, K) =                      &
     &                 AOD_SCATTERING(J, I_INITIAL, K)
                    ENDDO
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
      ENDIF
!
      RETURN
      END SUBROUTINE R2_COMPRESS_SPECTRUM
#endif
#endif
