#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to apply a path-length scaling to the continuum.
!
! Method:
!       The scaling function is calculated. This is multpiled by a
!       suitable "amount" of continuum incorporating a broadening
!       density.
!
! Current owner of code: James Manners
!
! description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE RESCALE_CONTINUUM(N_PROFILE, N_LAYER, I_CONTINUUM      &
     &   , P, T, I_TOP                                                  &
     &   , DENSITY, MOLAR_DENSITY_WATER, MOLAR_DENSITY_FRN              &
     &   , WATER_FRAC                                                   &
     &   , AMOUNT_CONTINUUM                                             &
     &   , I_FNC                                                        &
     &   , P_REFERENCE, T_REFERENCE, SCALE_PARAMETER                    &
     &   , ND_PROFILE, ND_LAYER                                         &
     &   , ND_SCALE_VARIABLE                                            &
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
!           Size allocated for atmospheric profiles
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
     &  , ND_SCALE_VARIABLE
!           Size allocated for scaling variables
!
!     Include header files.
#include "c_kinds.h"
#include "physical_constants_0_ccf3z.h"
#include "continuum_pcf3z.h"
#include "scale_fnc_pcf3z.h"
!
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , I_CONTINUUM                                                   &
!           Continuum type
     &  , I_FNC                                                         &
!           Scaling function
     &  , I_TOP
!           Top `index' of arrays
      REAL  (Real64), INTENT(IN) ::                                     &
     &    WATER_FRAC(ND_PROFILE, ND_LAYER)                              &
!           Mass fraction of water
     &  , P(ND_PROFILE, ND_LAYER)                                       &
!           Pressure
     &  , T(ND_PROFILE, ND_LAYER)                                       &
!           Temperature
     &  , DENSITY(ND_PROFILE, ND_LAYER)                                 &
!           Overall density
     &  , MOLAR_DENSITY_WATER(ND_PROFILE, ND_LAYER)                     &
!           Molar density of water vapour
     &  , MOLAR_DENSITY_FRN(ND_PROFILE, ND_LAYER)                       &
!           Molar density of foreign species
     &  , P_REFERENCE                                                   &
!           Reference pressure
     &  , T_REFERENCE                                                   &
!           Reference pressure
     &  , SCALE_PARAMETER(ND_SCALE_VARIABLE)
!           Scaling paramters
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    AMOUNT_CONTINUUM(ND_PROFILE, ND_LAYER)
!           Amount of continuum
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
      IF (I_FNC == IP_SCALE_POWER_LAW) THEN
        DO I=I_TOP, N_LAYER
          DO L=1, N_PROFILE
            AMOUNT_CONTINUUM(L, I)                                      &
     &        =EXP(SCALE_PARAMETER(1)*LOG(P(L, I)/P_REFERENCE)          &
     &        +SCALE_PARAMETER(2)*LOG(T(L, I)/T_REFERENCE))
          ENDDO
        ENDDO
      ELSE IF(I_FNC == IP_SCALE_POWER_QUAD) THEN
        DO I=I_TOP, N_LAYER
          DO L=1, N_PROFILE
            AMOUNT_CONTINUUM(L, I)                                      &
     &        =EXP(SCALE_PARAMETER(1)*LOG(P(L, I)/P_REFERENCE))         &
     &        *(1.0E+00_Real64+SCALE_PARAMETER(2)*(T(L, I)              &
     &        /T_REFERENCE-1.0E+00_Real64)                              &
     &        +SCALE_PARAMETER(3)*(T(L, I)                              &
     &        /T_REFERENCE-1.0E+00_Real64)**2)
          ENDDO
        ENDDO
      ENDIF
!
      IF (I_CONTINUUM == IP_SELF_CONTINUUM) THEN
        DO I=1, N_LAYER
          DO L=1, N_PROFILE
            AMOUNT_CONTINUUM(L, I)=AMOUNT_CONTINUUM(L, I)               &
     &        *MOLAR_DENSITY_WATER(L, I)*WATER_FRAC(L, I)
          ENDDO
        ENDDO
      ELSE IF (I_CONTINUUM == IP_FRN_CONTINUUM) THEN
        DO I=1, N_LAYER
          DO L=1, N_PROFILE
            AMOUNT_CONTINUUM(L, I)=AMOUNT_CONTINUUM(L, I)               &
     &        *MOLAR_DENSITY_FRN(L, I)*WATER_FRAC(L, I)
          ENDDO
        ENDDO
      ELSE IF (I_CONTINUUM == IP_N2_CONTINUUM) THEN
        DO I=1, N_LAYER
          DO L=1, N_PROFILE
            AMOUNT_CONTINUUM(L, I)=AMOUNT_CONTINUUM(L, I)               &
     &        *N2_MASS_FRAC*DENSITY(L, I)
          ENDDO
        ENDDO
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE RESCALE_CONTINUUM
#endif
