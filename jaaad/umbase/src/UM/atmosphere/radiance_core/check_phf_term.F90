#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to check the number of terms in the phase function.
!
! Purpose:
!   This subroutine checks the prescription of the phase function
!   against the specified options to ensure that information is
!   present to define all required moments.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE CHECK_PHF_TERM(IERR                                    &
     &  , L_AEROSOL, N_AEROSOL, I_AEROSOL_PARAMETRIZATION               &
     &  , N_AEROSOL_PHF_TERM                                            &
#if defined(OBSERVED)
     &  , N_PHASE_TERM_AEROSOL_PRSC                                     &
#endif
     &  , L_CLOUD, N_CONDENSED, I_CONDENSED_PARAM, I_PHASE_CMP          &
     &  , CONDENSED_N_PHF                                               &
#if defined(OBSERVED)
     &  , N_PHASE_TERM_DROP_PRSC, N_PHASE_TERM_ICE_PRSC                 &
#endif
     &  , N_ORDER_PHASE, L_HENYEY_GREENSTEIN_PF                         &
     &  , L_RESCALE, N_ORDER_FORWARD, L_SOLAR_PHF, N_ORDER_PHASE_SOLAR  &
     &  , ND_AEROSOL_SPECIES, ND_CONDENSED                              &
     &  )
!
!
      IMPLICIT NONE
!
!
!     Include header files
#include "c_kinds.h"
#include "def_std_io_icf3z.h"
#include "error_pcf3z.h"
#include "aerosol_parametrization_pcf3z.h"
#include "cloud_parametrization_pcf3z.h"
#include "ice_cloud_param_pcf3z.h"
#include "phase_pcf3z.h"
!
!
!     Dummy arguments:
      INTEGER                                                           &
     &    IERR
!           Error flag
!
!     Dimensions of arrays:
      INTEGER, INTENT(IN) ::                                            &
     &    ND_AEROSOL_SPECIES                                            &
!           Size allocated for species of aerosols
     &  , ND_CONDENSED
!           Size allocated for condensed components
!
!
!     Generic variables:
      LOGICAL, INTENT(IN) ::                                            &
     &    L_HENYEY_GREENSTEIN_PF                                        &
!           Flag for Henyey-Greenstein phase functions
     &  , L_RESCALE                                                     &
!           Flag for rescaling of the phase functions
     &  , L_SOLAR_PHF
!           Flag to use a separate treatment of the solar beam
      INTEGER, INTENT(IN) ::                                            &
     &    N_ORDER_PHASE                                                 &
!           Order of terms required in the phase function
     &  , N_ORDER_FORWARD                                               &
!           Order of the term in the phase function used for rescaling
     &  , N_ORDER_PHASE_SOLAR
!           Order of the phase function used in solar calculations
!
!     Aerosol Fields
      LOGICAL, INTENT(IN) ::                                            &
     &    L_AEROSOL
!           Flag to use aerosols
      INTEGER, INTENT(IN) ::                                            &
     &    N_AEROSOL                                                     &
!           Number of aerosols
     &  , I_AEROSOL_PARAMETRIZATION(ND_AEROSOL_SPECIES)                 &
!           Parametrizations adopted for aerosols
     &  , N_AEROSOL_PHF_TERM(ND_AEROSOL_SPECIES)
!           Number of terms in the phase function
#if defined(OBSERVED)
      INTEGER, INTENT(IN) ::                                            &
     &    N_PHASE_TERM_AEROSOL_PRSC(ND_AEROSOL_SPECIES)
!           Number of terms in the prescribed phase functions
!           for each species
#endif
!
!     Cloudy Fields
      LOGICAL, INTENT(IN) ::                                            &
     &    L_CLOUD
!           Flag to include clouds
      INTEGER, INTENT(IN) ::                                            &
     &    N_CONDENSED                                                   &
!           Number of condensed components
     &  , I_CONDENSED_PARAM(ND_CONDENSED)                               &
!           Parametrizations adopted for condensed components
     &  , I_PHASE_CMP(ND_CONDENSED)                                     &
!           Phases of the condensed components
     &  , CONDENSED_N_PHF(ND_CONDENSED)
!           Number of terms in the phase function
#if defined(OBSERVED)
      INTEGER, INTENT(IN) ::                                            &
     &    N_PHASE_TERM_DROP_PRSC                                        &
!           Number of terms in the prescribed phase functions
!           for each droplets
     &  , N_PHASE_TERM_ICE_PRSC
!           Number of terms in the prescribed phase functions
!           for each ice crystals
#endif
!
!
!     Local variables:
      LOGICAL                                                           &
     &    L_INADEQUATE
!           Flag for inadequate information
      INTEGER                                                           &
     &    J                                                             &
!           Loop variable
     &  , N_ORDER_REQUIRED
!           Order of terms which are required in the phase function
!
!
!
!     Determine the order of terms for which information in the
!     phase function will be required.
      IF (L_HENYEY_GREENSTEIN_PF) THEN
        N_ORDER_REQUIRED=1
      ELSE
        IF (L_RESCALE) THEN
          N_ORDER_REQUIRED=MAX(N_ORDER_PHASE, N_ORDER_FORWARD)
        ELSE
          N_ORDER_REQUIRED=N_ORDER_PHASE
        ENDIF
!       If the solar beam is to be treated separately more terms
!       may be required.
        IF (L_SOLAR_PHF) THEN
          N_ORDER_REQUIRED=                                             &
     &             MAX(N_ORDER_PHASE, N_ORDER_PHASE_SOLAR)
          IF (L_RESCALE) N_ORDER_REQUIRED=N_ORDER_REQUIRED+1
        ENDIF
      ENDIF
!
!     If aerosols are included carry out the required checks.
      IF (L_AEROSOL) THEN
        L_INADEQUATE=.FALSE.
        DO J=1, N_AEROSOL
          IF ( (I_AEROSOL_PARAMETRIZATION(J) ==                         &
     &          IP_AEROSOL_PARAM_DRY).OR.                               &
     &         (I_AEROSOL_PARAMETRIZATION(J) ==                         &
     &          IP_AEROSOL_PARAM_MOIST) ) THEN
!           In this case information will be extended as a
!           Henyey-Greenstein phase function; and the available
!           information will include the asymmetry.
            CONTINUE
          ELSE IF ( (I_AEROSOL_PARAMETRIZATION(J) ==                    &
     &               IP_AEROSOL_PARAM_PHF_DRY).OR.                      &
     &              (I_AEROSOL_PARAMETRIZATION(J) ==                    &
     &               IP_AEROSOL_PARAM_PHF_MOIST) ) THEN
            L_INADEQUATE=(N_ORDER_REQUIRED >  N_AEROSOL_PHF_TERM(J))
#if defined(OBSERVED)
          ELSE IF (I_AEROSOL_PARAMETRIZATION(J) ==                      &
     &           IP_AEROSOL_UNPARAMETRIZED) THEN
            L_INADEQUATE=(N_ORDER_REQUIRED >                            &
     &                    N_PHASE_TERM_AEROSOL_PRSC(J))
#endif
          ENDIF
!
          IF (L_INADEQUATE) THEN
            WRITE(IU_ERR, '(/A, /A, I3, A)')                            &
     &        '*** Error: There is not enough information to define'    &
     &        , 'the phase function for aerosol ', J                    &
     &        , ' to the desired order.'
            IERR=I_ERR_FATAL
            RETURN
          ENDIF
        ENDDO
      ENDIF
!
      IF (L_CLOUD) THEN
        L_INADEQUATE=.FALSE.
        DO J=1, N_CONDENSED
          IF (I_PHASE_CMP(J) == IP_PHASE_WATER) THEN
            IF ( (I_CONDENSED_PARAM(J) == IP_SLINGO_SCHRECKER).OR.      &
     &           (I_CONDENSED_PARAM(J) == IP_ACKERMAN_STEPHENS).OR.     &
     &           (I_CONDENSED_PARAM(J) == IP_DROP_PADE_2) ) THEN
!             The phase function will be extended as a
!             Henyey-Greenstein phase function from information
!             already present.
              CONTINUE
            ELSE IF ( (I_CONDENSED_PARAM(J) ==                          &
     &                   IP_SLINGO_SCHR_PHF) ) THEN
              L_INADEQUATE=(N_ORDER_REQUIRED >  CONDENSED_N_PHF(J))
#if defined(OBSERVED)
            ELSE IF (I_CONDENSED_PARAM(J) ==                            &
     &               IP_DROP_UNPARAMETRIZED) THEN
              L_INADEQUATE=(N_ORDER_REQUIRED >  N_PHASE_TERM_DROP_PRSC)
#endif
            ENDIF
            IF (L_INADEQUATE) THEN
              WRITE(IU_ERR, '(/A, /A, I3, A, /A)')                      &
     &          '*** Error: There is not enough information to define'  &
     &          , 'the phase function for condensed species ', J        &
     &          , ' (water droplets) ', 'to the desired order.'
              IERR=I_ERR_FATAL
              RETURN
            ENDIF
          ELSE IF (I_PHASE_CMP(J) == IP_PHASE_ICE) THEN
            IF ( (I_CONDENSED_PARAM(J) == IP_SLINGO_SCHRECKER_ICE).OR.  &
     &           (I_CONDENSED_PARAM(J) == IP_ICE_ADT).OR.               &
     &           (I_CONDENSED_PARAM(J) == IP_ICE_FU_SOLAR).OR.          &
     &           (I_CONDENSED_PARAM(J) == IP_ICE_FU_IR).OR.             &
     &           (I_CONDENSED_PARAM(J) == IP_ICE_ADT_10) ) THEN
!             The phase function will be extended as a
!             Henyey-Greenstein phase function from information
!             already present.
              CONTINUE
            ELSE IF ( (I_CONDENSED_PARAM(J) ==                          &
     &                   IP_SLINGO_SCHR_ICE_PHF).OR.                    &
     &                (I_CONDENSED_PARAM(J) == IP_ICE_FU_PHF) ) THEN
              L_INADEQUATE=(N_ORDER_REQUIRED >  CONDENSED_N_PHF(J))
#if defined(OBSERVED)
            ELSE IF (I_CONDENSED_PARAM(J) ==                            &
     &               IP_ICE_UNPARAMETRIZED) THEN
              L_INADEQUATE=(N_ORDER_REQUIRED >  N_PHASE_TERM_ICE_PRSC)
#endif
            ENDIF
            IF (L_INADEQUATE) THEN
              WRITE(IU_ERR, '(/A, /A, I3, A)')                          &
     &          '*** Error: There is not enough information to define'  &
     &          , 'the phase function for condensed species ', J        &
     &          , ' (ice crystals) to the desired order.'
              IERR=I_ERR_FATAL
              RETURN
            ENDIF
          ENDIF
!
        ENDDO
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE CHECK_PHF_TERM
#endif
