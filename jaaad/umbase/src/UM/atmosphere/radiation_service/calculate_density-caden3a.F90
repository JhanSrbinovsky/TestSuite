#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate densities.
!
! Method:
!       This routine calculates the density of air and the molar
!       densities of the broadening species for the self and foreign-
!       broadened continua using the gas law including the effect
!       of water vapour.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE CALCULATE_DENSITY(N_PROFILE, N_LAYER, L_CONTINUUM      &
     &   , WATER_FRAC, P, T, I_TOP                                      &
     &   , DENSITY, MOLAR_DENSITY_WATER, MOLAR_DENSITY_FRN              &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
!
      IMPLICIT NONE
!
!
!     SIZES OF DUMMY ARRAYS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     INCLUDE COMDECKS
#include "c_r_cp.h"
#include "c_epslon.h"
#include "phycn03a.h"
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             NUMBER OF LAYERS
     &   , I_TOP
!             TOP VERTICAL INDEX
      LOGICAL                                                           &
     &     L_CONTINUUM
!             CONTINUUM FLAG
      REAL                                                              &
                !, INTENT(IN)
     &     WATER_FRAC(NPD_PROFILE, 0: NPD_LAYER)                        &
!             MASS FRACTION OF WATER
     &   , P(NPD_PROFILE, 0: NPD_LAYER)                                 &
!             PRESSURE
     &   , T(NPD_PROFILE, 0: NPD_LAYER)
!             TEMPERATURE
      REAL                                                              &
                !, INTENT(OUT)
     &     DENSITY(NPD_PROFILE, 0: NPD_LAYER)                           &
!             AIR DENSITY
     &   , MOLAR_DENSITY_WATER(NPD_PROFILE, 0: NPD_LAYER)               &
!             MOLAR DENSITY OF WATER
     &   , MOLAR_DENSITY_FRN(NPD_PROFILE, 0: NPD_LAYER)
!             MOLAR DENSITY OF FOREIGN SPECIES
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     L                                                            &
!             LOOP VARIABLE
     &   , I
!             LOOP VARIABLE
!
!
!     FIND THE AIR DENSITY FIRST.
      DO I=I_TOP, N_LAYER
         DO L=1, N_PROFILE
            DENSITY(L, I)=P(L, I)/(R*T(L, I)                            &
     &         *(1.0E+00+C_VIRTUAL*WATER_FRAC(L, I)))
         ENDDO
      ENDDO
!
      IF (L_CONTINUUM) THEN
         DO I=I_TOP, N_LAYER
            DO L=1, N_PROFILE
               MOLAR_DENSITY_FRN(L, I)=DENSITY(L, I)                    &
     &            *(1.0E+00-WATER_FRAC(L, I))/MOL_WEIGHT_AIR
               MOLAR_DENSITY_WATER(L, I)=DENSITY(L, I)                  &
     &            *WATER_FRAC(L, I)/(EPSILON*MOL_WEIGHT_AIR)
            ENDDO
         ENDDO
      ENDIF
!
!
      RETURN
      END SUBROUTINE CALCULATE_DENSITY
#endif
#endif
