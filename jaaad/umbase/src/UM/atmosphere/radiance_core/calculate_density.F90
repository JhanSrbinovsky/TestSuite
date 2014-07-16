#if defined(A70_1Z)
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
! Current owner of code: James Manners
!
! Description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE CALCULATE_DENSITY(N_PROFILE, N_LAYER, L_CONTINUUM      &
     &   , WATER_FRAC, P, T, I_TOP                                      &
     &   , DENSITY, MOLAR_DENSITY_WATER, MOLAR_DENSITY_FRN              &
     &   , ND_PROFILE, ND_LAYER                                         &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Maximum number of profiles
     &  , ND_LAYER
!           Maximum number of layers
!
!     Include header files.
#include "c_kinds.h"
#if defined(UM)
#include "c_r_cp.h"
#include "c_epslon.h"
#endif
#include "physical_constants_0_ccf3z.h"
#if defined(STANDARD)
#include "physical_constants_1_ccf3z.h"
#endif
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , I_TOP
!           Top vertical `index'
      LOGICAL                                                           &
     &    L_CONTINUUM
!           Continuum flag
      REAL  (Real64), INTENT(IN) ::                                     &
     &    WATER_FRAC(ND_PROFILE, ND_LAYER)                              &
!           Mass fraction of water
     &  , P(ND_PROFILE, ND_LAYER)                                       &
!           Pressure
     &  , T(ND_PROFILE, ND_LAYER)
!           Temperature
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    DENSITY(ND_PROFILE, ND_LAYER)                                 &
!           Air density
     &  , MOLAR_DENSITY_WATER(ND_PROFILE, ND_LAYER)                     &
!           Molar density of water
     &  , MOLAR_DENSITY_FRN(ND_PROFILE, ND_LAYER)
!           Molar density of foreign species
!
!     Local variables.
      INTEGER                                                           &
     &    L                                                             &
!           Loop variable
     &  , I
!           Loop variable
!
!
!
!     find the air density first.
      DO I=I_TOP, N_LAYER
        DO L=1, N_PROFILE
#if defined(STANDARD)
           DENSITY(L, I)=P(L, I)/(R_GAS_DRY*T(L, I)                     &
     &       *(1.0E+00_Real64+(RATIO_MOLAR_WEIGHT-1.00E+00_Real64)      &
     &       *WATER_FRAC(L, I)))
#endif
#if defined(UM)
           DENSITY(L, I)=P(L, I)/(R*T(L, I)                             &
     &       *(1.0E+00_Real64+C_VIRTUAL*WATER_FRAC(L, I)))
#endif
        ENDDO
      ENDDO
!
      IF (L_CONTINUUM) THEN
        DO I=I_TOP, N_LAYER
          DO L=1, N_PROFILE
            MOLAR_DENSITY_FRN(L, I)=DENSITY(L, I)                       &
     &        *(1.0E+00_Real64-WATER_FRAC(L, I))/MOL_WEIGHT_AIR
#if defined(STANDARD)
            MOLAR_DENSITY_WATER(L, I)=DENSITY(L, I)                     &
     &        *WATER_FRAC(L, I)*(RATIO_MOLAR_WEIGHT/MOL_WEIGHT_AIR)
#endif
#if defined(UM)
            MOLAR_DENSITY_WATER(L, I)=DENSITY(L, I)                     &
     &        *WATER_FRAC(L, I)/(EPSILON*MOL_WEIGHT_AIR)
#endif
          ENDDO
        ENDDO
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE CALCULATE_DENSITY
#endif
