#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to scale amounts of absorbers.
!
! Method:
!       The mixing ratio is multiplied by a factor determined
!       by the type of scaling selected.
!
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SCALE_ABSORB(IERR, N_PROFILE, N_LAYER                  &
     &  , GAS_MIX_RATIO, P, T, I_TOP                                    &
     &  , GAS_FRAC_RESCALED                                             &
     &  , I_FNC, P_REFERENCE, T_REFERENCE, SCALE_PARAMETER              &
     &  , L_DOPPLER, DOPPLER_CORRECTION                                 &
     &  , ND_PROFILE, ND_LAYER                                          &
     &  , ND_SCALE_VARIABLE                                             &
     &  )
!
!
      IMPLICIT NONE
!
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for profiles
     &  , ND_LAYER                                                      &
!           Size allocated for layers
     &  , ND_SCALE_VARIABLE
!           Size allocated for of scaling variables
!
!     Include header files.
#include "c_kinds.h"
#include "def_std_io_icf3z.h"
#include "scale_fnc_pcf3z.h"
#include "error_pcf3z.h"
!
!     Dummy arguments.
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , I_FNC                                                         &
!           Type of scaling function
     &  , I_TOP
!           Uppermost `index' for scaling (this will be 1 for fields
!           Given in layers, as in the unified model, or 0 for
!           Fields given at the boundaries of layers)
      LOGICAL, INTENT(IN) ::                                            &
     &    L_DOPPLER
!           Flag for Doppler term
      REAL  (Real64), INTENT(IN) ::                                     &
     &    GAS_MIX_RATIO(ND_PROFILE, ND_LAYER)                           &
!           Mass mixing ratio of gas
     &  , P(ND_PROFILE, ND_LAYER)                                       &
!           Pressure
     &  , T(ND_PROFILE, ND_LAYER)                                       &
!           Temperature
     &  , P_REFERENCE                                                   &
!           Reference pressure
     &  , T_REFERENCE                                                   &
!           Reference temperature
     &  , SCALE_PARAMETER(ND_SCALE_VARIABLE)                            &
!           Scaling paramters
     &  , DOPPLER_CORRECTION
!           Doppler-broadening correction
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    GAS_FRAC_RESCALED(ND_PROFILE, ND_LAYER)
!           Mass fraction of gas
!
!     Local variables.
      INTEGER                                                           &
     &    L                                                             &
!           Loop variable
     &  , I
!           Loop variable
      REAL  (Real64) ::                                                 &
     &    PRESSURE_OFFSET
!           Offset to pressure
!
!
!
!     Set the offset to the pressure for the Doppler correction.
      IF (L_DOPPLER) THEN
        PRESSURE_OFFSET=DOPPLER_CORRECTION
      ELSE
        PRESSURE_OFFSET=0.0E+00_Real64
      ENDIF
!
!     The array gas_frac_rescaled is used initially to hold only the
!     scaling functions, and only later is it multiplied by the
!     mixing ratios
!
      IF (I_FNC == IP_SCALE_POWER_LAW) THEN
        DO I=I_TOP, N_LAYER
          DO L=1, N_PROFILE
            GAS_FRAC_RESCALED(L, I)=                                    &
     &        EXP(SCALE_PARAMETER(1)*LOG((P(L, I)                       &
     &        +PRESSURE_OFFSET)                                         &
     &        /(P_REFERENCE+PRESSURE_OFFSET))                           &
     &        +SCALE_PARAMETER(2)*LOG(T(L, I)/T_REFERENCE))
          ENDDO
        ENDDO
      ELSE IF (I_FNC == IP_SCALE_FNC_NULL) THEN
        RETURN
      ELSE IF (I_FNC == IP_SCALE_POWER_QUAD) THEN
        DO I=I_TOP, N_LAYER
          DO L=1, N_PROFILE
            GAS_FRAC_RESCALED(L, I)=                                    &
     &        EXP(SCALE_PARAMETER(1)                                    &
     &        *LOG((P(L, I)+PRESSURE_OFFSET)                            &
     &        /(P_REFERENCE+PRESSURE_OFFSET)))                          &
     &        *(1.0E+00_Real64                                          &
     &        +SCALE_PARAMETER(2)*(T(L, I)/T_REFERENCE-1.0E+00_Real64)  &
     &        +SCALE_PARAMETER(3)*(T(L, I)                              &
     &        /T_REFERENCE-1.0E+00_Real64)**2)
          ENDDO
        ENDDO
      ELSE IF (I_FNC == IP_SCALE_DOPPLER_QUAD) THEN
!       There is no Doppler term here since it is implicitly included
!       in the scaling.
        DO I=I_TOP, N_LAYER
          DO L=1, N_PROFILE
            GAS_FRAC_RESCALED(L, I)=                                    &
     &        EXP(SCALE_PARAMETER(1)                                    &
     &        *LOG((P(L, I)+SCALE_PARAMETER(2))                         &
     &        /(P_REFERENCE+SCALE_PARAMETER(2))))                       &
     &        *(1.0E+00_Real64                                          &
     &        +SCALE_PARAMETER(3)*(T(L, I)/T_REFERENCE-1.0E+00_Real64)  &
     &        +SCALE_PARAMETER(4)*(T(L, I)                              &
     &        /T_REFERENCE-1.0E+00_Real64)**2)
          ENDDO
        ENDDO
      ELSE
        WRITE(IU_ERR, '(/A)')                                           &
     &    '*** Error: an illegal type of scaling has been given.'
        IERR=I_ERR_FATAL
        RETURN
      ENDIF
!
!     Multiply by the mixing ratio and limit negative scalings.
      DO I=N_LAYER, 1, -1
        DO L=1, N_PROFILE
          GAS_FRAC_RESCALED(L, I)=MAX(0.0E+00_Real64                    &
     &      , GAS_FRAC_RESCALED(L, I)*GAS_MIX_RATIO(L, I))
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE SCALE_ABSORB
#endif
