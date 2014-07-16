#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set moist aerosol properties independent of bands.
!
! Method:
!       The mean relative humidities are calculated and pointers to
!       the lookup tables are set.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SET_MOIST_AEROSOL_PROPERTIES(IERR                      &
     &  , N_PROFILE, N_LAYER                                            &
     &  , N_AEROSOL, I_AEROSOL_PARAMETRIZATION, NHUMIDITY               &
     &  , WATER_MIX_RATIO, T, P, DELTA_HUMIDITY                         &
     &  , MEAN_REL_HUMIDITY, I_HUMIDITY_POINTER                         &
     &  , ND_PROFILE, ND_LAYER, ND_AEROSOL_SPECIES                      &
     &  )
!
!
      IMPLICIT NONE
!
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Maximum number of profiles
     &  , ND_LAYER                                                      &
!           Maximum number of layers
     &  , ND_AEROSOL_SPECIES
!           Maximum number of aerosols
!
!     Include header files.
#include "c_kinds.h"
#include "aerosol_parametrization_pcf3z.h"
#include "def_std_io_icf3z.h"
#include "error_pcf3z.h"
!
!     Dummy arguments.
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
!
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , N_AEROSOL                                                     &
!           Number of aerosol species
     &  , I_AEROSOL_PARAMETRIZATION(ND_AEROSOL_SPECIES)                 &
!           Parametrizations of aerosol species
     &  , NHUMIDITY(ND_AEROSOL_SPECIES)
!           Number of humidity values
      INTEGER, INTENT(OUT) ::                                           &
     &    I_HUMIDITY_POINTER(ND_PROFILE, ND_LAYER)
!           Pointers to look-up tables
      REAL  (Real64), INTENT(IN) ::                                     &
     &    WATER_MIX_RATIO(ND_PROFILE, ND_LAYER)                         &
!           Mixing ratio of water vapour
     &  , T(ND_PROFILE, ND_LAYER)                                       &
!           Temperatures
     &  , P(ND_PROFILE, ND_LAYER)
!           Pressures
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    MEAN_REL_HUMIDITY(ND_PROFILE, ND_LAYER)                       &
!           Mean humidities of layers
     &  , DELTA_HUMIDITY
!           Increment in humidity
!
!
!     local variables.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , J                                                             &
!           Loop variable
     &  , L                                                             &
!           Loop variable
     &  , NHUMIDITY_COMMON
!           Common number of humidities for moist aerosols
      REAL  (Real64) ::                                                 &
     &    MIX_RATIO_SAT(ND_PROFILE, ND_LAYER)
!           Saturated humidity mixing ratio
!
!     Subroutines called:
#if defined(STANDARD)
      EXTERNAL                                                          &
     &    QSAT_GILL
#endif
#if defined(UM)
      EXTERNAL                                                          &
     &    QSAT_WAT
#endif
!
!
!
!     Set up array of pointers to `include' the effects of humidity.
!     Calculate the saturated mixing ratio.
#if defined(STANDARD)
! DEPENDS ON: qsat_gill
      CALL QSAT_GILL(MIX_RATIO_SAT, T, P                                &
     &  , N_PROFILE, N_LAYER                                            &
     &  , ND_PROFILE, ND_LAYER)
#endif
#if defined(UM)
      DO I=1, N_LAYER
! DEPENDS ON: qsat_wat
        CALL QSAT_WAT(MIX_RATIO_SAT(1, I), T(1, I), P(1, I)             &
     &    , N_PROFILE)
      ENDDO
#endif
!
!     Determine the number of humidities to be used for moist
!     aerosols. This must be the same for all moist aerosols
!     in the current version of the code.
      NHUMIDITY_COMMON=0
      DO J=1, N_AEROSOL
        IF (I_AEROSOL_PARAMETRIZATION(J) == IP_AEROSOL_PARAM_MOIST)     &
     &       THEN
          IF (NHUMIDITY(J) >  0) THEN
!           Set the actual common value.
            IF (NHUMIDITY_COMMON == 0) THEN
              NHUMIDITY_COMMON=NHUMIDITY(J)
            ELSE IF (NHUMIDITY(J) /= NHUMIDITY_COMMON) THEN
!             There is an inconsistency.
              WRITE(IU_ERR, '(/A)')                                     &
     &          '***Error: The look-up tables for moist aerosols '      &
     &          , 'are of different sizes. Tgis is not permitted.'
              IERR=I_ERR_FATAL
              RETURN
            ENDIF
          ENDIF
        ENDIF
      ENDDO
!     The look-up table is assumed to be uniform in humidity.
      DELTA_HUMIDITY=1.0E+00_Real64                                     &
     &  /(REAL(NHUMIDITY_COMMON, Real64)-1.0E+00_Real64)
      DO I=1, N_LAYER
        DO L=1, N_PROFILE
          MEAN_REL_HUMIDITY(L, I)                                       &
     &      =WATER_MIX_RATIO(L,I)*(1.0E+00_Real64-MIX_RATIO_SAT(L,I))   &
     &      /((1.0E+00_Real64-WATER_MIX_RATIO(L,I))*MIX_RATIO_SAT(L,I))
!         Check that the mean relative humidity
!         does not exceed or equal 1.0.
          MEAN_REL_HUMIDITY(L, I)=MIN(MEAN_REL_HUMIDITY(L, I)           &
     &      , 9.9999E-01_Real64)
          I_HUMIDITY_POINTER(L, I)=1                                    &
     &      +INT(MEAN_REL_HUMIDITY(L, I)*(NHUMIDITY_COMMON-1))
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE SET_MOIST_AEROSOL_PROPERTIES
#endif
