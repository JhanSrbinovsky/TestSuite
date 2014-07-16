#if defined(A70_1A) || defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to compute the aerosol optical depth (AOD).
!
! Method:
!       The optical depth is the vertical integral of the aerosol
!       specific extinction weighted by the aerosol mass.
!
! Current Owner of Code: N. Bellouin
!
! History:
!
!       6.2             15/11/05                New deck.
!                                                          N. Bellouin
!
! Description of Code:
!   FORTRAN 77
!
!-----------------------------------------------------------------------
      SUBROUTINE COMPUTE_AOD (                                          &
! actual and fixed array dimensions
     &   N_AEROSOL, NPD_AEROSOL,                                        &
     &   N_AOD_WAVEL, NPD_AOD_WAVEL,                                    &
     &   N_PROFILE, NPD_PROFILE,                                        &
     &   N_LAYER, NPD_LAYER,                                            &
     &   N_HUMIDITIES, NPD_HUMIDITIES,                                  &
! variables with intent in
     &   TYPE_AEROSOL,                                                  &
     &   I_AOD_TYPE, I_TYPE_WANTED,                                     &
     &   AEROSOL_MIX_RATIO, D_MASS,                                     &
     &   I_AEROSOL_PARAMETRIZATION,                                     &
     &   AOD_ABSORPTION, AOD_SCATTERING,                                &
     &   I_HUMIDITY_POINTER, MEAN_RH,                                   &
     &   HUMIDITIES, DELTA_HUMIDITY,                                    &
! variable with intent out
     &   AOD                                                            &
     & )
      IMPLICIT NONE

! declaration of arguments
!
! array dimensions (actual and fixed)
!
!  number of aerosol components
      INTEGER N_AEROSOL
      INTEGER NPD_AEROSOL
!  number of wavelengths at which the AOD is to be computed
      INTEGER N_AOD_WAVEL
      INTEGER NPD_AOD_WAVEL
!  number of grid-boxes
      INTEGER N_PROFILE
      INTEGER NPD_PROFILE
!  number of vertical layers
      INTEGER N_LAYER
      INTEGER NPD_LAYER
!  number of humidities used in moist aerosol parameterisation
      INTEGER N_HUMIDITIES
      INTEGER NPD_HUMIDITIES
!
! arguments with intent in
!
!  array giving the type of each aerosol component (see AERCMP3A)
      INTEGER TYPE_AEROSOL(NPD_AEROSOL)
!  array connecting an aerosol component to an aerosol type
      INTEGER I_AOD_TYPE(NPD_AEROSOL)
!  aerosol type to be considered in this call
      INTEGER I_TYPE_WANTED
!  aerosol component mass mixing ratio
      REAL AEROSOL_MIX_RATIO(NPD_PROFILE, 0:NPD_LAYER, NPD_AEROSOL)
!  mass thickness of vertical layers
      REAL D_MASS(NPD_PROFILE, NPD_LAYER)
!  aerosol parameterisation (dry or moist)
      INTEGER I_AEROSOL_PARAMETRIZATION(NPD_AEROSOL)
!  aerosol specific coefficients for absorption and scattering
!  (monochromatic)
      REAL AOD_ABSORPTION(NPD_HUMIDITIES, NPD_AEROSOL, NPD_AOD_WAVEL)
      REAL AOD_SCATTERING(NPD_HUMIDITIES, NPD_AEROSOL, NPD_AOD_WAVEL)
!  mean relative humidities, and indices to the look-up tables
!  it may be the grid-box mean or clear-sky mean relative humidity,
!  depending on calculations made in FLUX_CALC
      INTEGER I_HUMIDITY_POINTER(NPD_PROFILE, NPD_LAYER)
      REAL MEAN_RH(NPD_PROFILE, NPD_LAYER)
      REAL HUMIDITIES(NPD_HUMIDITIES, NPD_AEROSOL)
      REAL DELTA_HUMIDITY
!
! arguments with intent out
!
!  computed aerosol optical depth
      REAL AOD(NPD_PROFILE, NPD_AOD_WAVEL)

!
! local variables
!
!  loop indices
      INTEGER I, J, JJ, K, L
!  number of aerosol components included in the current aerosol type
      INTEGER N_AER_IN_TYPE
!  indices of those components (at most, all components in the same
!                               type)
      INTEGER N_AER_IN_TYPE_INDEX(NPD_AEROSOL)
!  aerosol extinction (absorption + scattering)
      REAL EXTINCTION
!  variables needed for the interpolation on humidities
      INTEGER I_POINTER
      REAL WEIGHT_UPPER
      REAL WEIGHT_LOWER
!  parameters coding the different aerosol parameterisations
#include "aerprm3a.h"

!

      DO K = 1, N_AOD_WAVEL
        DO L = 1, N_PROFILE
          AOD(L, K) = 0.0
        END DO
      END DO

      ! FIND WHICH AEROSOL COMPONENT BELONGS TO THE CURRENT TYPE
      ! DO NOT CONSIDER CLIMATOLOGICAL AEROSOLS (TYPE_AEROSOL IS
      ! SMALLER THAN 10)
      N_AER_IN_TYPE = 0
      DO J = 1, N_AEROSOL
        IF(TYPE_AEROSOL(J) >= 10 .AND.                                  &
     &     I_AOD_TYPE(J) == I_TYPE_WANTED) THEN
          N_AER_IN_TYPE = N_AER_IN_TYPE + 1
          N_AER_IN_TYPE_INDEX(N_AER_IN_TYPE) = J
        END IF
      END DO

      ! COMPUTE THE AOD IF AT LEAST ONE COMPONENT MATCHES
      ! (IF NO MATCH, THEN THE AOD WILL REMAIN AT ZERO)
      IF (N_AER_IN_TYPE  >   0) THEN

        DO JJ = 1, N_AER_IN_TYPE
          J = N_AER_IN_TYPE_INDEX(JJ)

          IF(I_AEROSOL_PARAMETRIZATION(J)  ==                           &
     &       IP_AEROSOL_PARAM_DRY) THEN

            ! Non-hygroscopic aerosol

            DO K = 1, N_AOD_WAVEL
              EXTINCTION = AOD_ABSORPTION(1, J, K) +                    &
     &                     AOD_SCATTERING(1, J, K)
              DO I = 1, N_LAYER
                DO L = 1, N_PROFILE
                  AOD(L, K) = AOD(L, K) +                               &
     &              AEROSOL_MIX_RATIO(L, I, J) *                        &
     &              D_MASS(L, I) * EXTINCTION
                END DO ! L
              END DO ! I
            END DO ! K

          ELSE IF(I_AEROSOL_PARAMETRIZATION(J)  ==                      &
     &            IP_AEROSOL_PARAM_MOIST) THEN

            ! Hygroscopic aerosol
            ! interpolation on the mean relative humidity

            DO K = 1, N_AOD_WAVEL
              DO I = 1, N_LAYER
                DO L = 1, N_PROFILE
                  I_POINTER = I_HUMIDITY_POINTER(L, I)
                  WEIGHT_UPPER = ( MEAN_RH(L, I)                        &
     &                   - HUMIDITIES(I_POINTER, J))                    &
     &                   / DELTA_HUMIDITY
                  WEIGHT_LOWER = 1.00E+00 - WEIGHT_UPPER

                  EXTINCTION =                                          &
     &              (AOD_ABSORPTION(I_POINTER, J, K) +                  &
     &               AOD_SCATTERING(I_POINTER, J, K))                   &
     &              * WEIGHT_LOWER + WEIGHT_UPPER *                     &
     &              (AOD_ABSORPTION(I_POINTER+1,J,K) +                  &
     &               AOD_SCATTERING(I_POINTER+1,J,K))

                  AOD(L, K) = AOD(L, K) +                               &
     &              AEROSOL_MIX_RATIO(L, I, J) *                        &
     &              D_MASS(L, I) * EXTINCTION
                END DO ! L
              END DO ! I
            END DO ! K
          END IF

        END DO ! JJ

      END IF ! N_AER_IN_TYPE > 0

      RETURN
      END SUBROUTINE COMPUTE_AOD
#endif
#endif
