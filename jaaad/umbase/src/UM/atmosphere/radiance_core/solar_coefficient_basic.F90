#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the basic coefficients for the solar beam.
!
! Method:
!       Straightforward.
!
! Current owner of code: James Manners
!
! Description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SOLAR_COEFFICIENT_BASIC(IERR                           &
     &  , N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST                        &
     &  , OMEGA, ASYMMETRY, SEC_0                                       &
     &  , I_2STREAM                                                     &
     &  , SUM, DIFF, LAMBDA                                             &
     &  , GAMMA_UP, GAMMA_DOWN                                          &
     &  , ND_PROFILE, ID_LT, ID_LB                                      &
     &  )
!
!
!
!
      IMPLICIT NONE
!
!
!     Sizes of dummy arrays.
      INTEGER                                                           &
     &    ND_PROFILE                                                    &
!           Size allocated for atmospheric profiles
     &  , ID_LT                                                         &
!           Topmost declared layer
     &  , ID_LB
!           Bottom declared layer
!
!     Include header files.
#include "c_kinds.h"
#include "two_stream_scheme_pcf3z.h"
#include "def_std_io_icf3z.h"
#include "error_pcf3z.h"
!
!
!     Dummy variables.
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , I_LAYER_FIRST                                                 &
!           First layer to consider
     &  , I_LAYER_LAST                                                  &
!           First layer to consider
     &  , I_2STREAM
!           Two-stream scheme
!
      REAL  (Real64), INTENT(IN) ::                                     &
     &    OMEGA(ND_PROFILE, ID_LT: ID_LB)                               &
!           Albedo of single scattering
     &  , ASYMMETRY(ND_PROFILE, ID_LT: ID_LB)                           &
!           Asymmetry
     &  , SEC_0(ND_PROFILE)
!           Secant of solar zenith angle
!
!     Basic two-stream coefficients:
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    SUM(ND_PROFILE, ID_LT: ID_LB)                                 &
!           Sum of two-stream coefficients
     &  , DIFF(ND_PROFILE, ID_LT: ID_LB)                                &
!           Difference of two-stream coefficients
     &  , LAMBDA(ND_PROFILE, ID_LT: ID_LB)
!           Lambda
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    GAMMA_UP(ND_PROFILE, ID_LT: ID_LB)                            &
!           Coefficient for upward radiation
     &  , GAMMA_DOWN(ND_PROFILE, ID_LT: ID_LB)
!           Coefficient for downwad radiation
!
!
!     Local variables.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , L
!           Loop variable
      REAL  (Real64) ::                                                 &
     &    KSI_0(ND_PROFILE, ID_LT: ID_LB)                               &
!           Difference in solar scattering fractions
     &  , FACTOR
!           Temporary variable
      REAL  (Real64) ::                                                 &
     &    ROOT_3
!           Square root of 3
      PARAMETER(                                                        &
     &    ROOT_3=1.7320508075688772E+00_Real64                          &
     &  )
!
!     Variables related to the treatment of ill-conditioning
      REAL  (Real64) ::                                                 &
     &    TOL_PERTURB
!           The tolerance used to judge where the two-stream
!           expressions for the solar source become ill-conditioned
!
!
!
!     Set the tolerances used in avoiding ill-conditioning, testing
!     on any variable.
      TOL_PERTURB=3.2E+01_Real64*EPSILON(SEC_0(1))
!
!     If LAMBDA is too close to SEC_0 it must be perturbed.
      DO I=I_LAYER_FIRST, I_LAYER_LAST
        DO L=1, N_PROFILE
          IF ((ABS(LAMBDA(L, I)-SEC_0(L))) <  TOL_PERTURB) THEN
            SUM(L, I)=(1.0E+00_Real64+TOL_PERTURB)*SUM(L, I)
            DIFF(L, I)=(1.0E+00_Real64+TOL_PERTURB)*DIFF(L, I)
            LAMBDA(L, I)=(1.0E+00_Real64+TOL_PERTURB)*LAMBDA(L, I)
          ENDIF
        ENDDO
      ENDDO
!
      IF ( (I_2STREAM == IP_EDDINGTON).OR.                              &
     &     (I_2STREAM == IP_ELSASSER).OR.                               &
     &     (I_2STREAM == IP_PIFM85).OR.                                 &
     &     (I_2STREAM == IP_2S_TEST).OR.                                &
     &     (I_2STREAM == IP_HEMI_MEAN).OR.                              &
     &     (I_2STREAM == IP_PIFM80) ) THEN
!
        DO I=I_LAYER_FIRST, I_LAYER_LAST
          DO L=1, N_PROFILE
            KSI_0(L, I)=1.5E+00_Real64*ASYMMETRY(L, I)/SEC_0(L)
          ENDDO
        ENDDO
!
      ELSE IF (I_2STREAM == IP_DISCRETE_ORD) THEN
!
        DO I=I_LAYER_FIRST, I_LAYER_LAST
          DO L=1, N_PROFILE
            KSI_0(L, I)=ROOT_3*ASYMMETRY(L, I)/SEC_0(L)
          ENDDO
        ENDDO
!
      ELSE
!
        WRITE(IU_ERR, '(/A)')                                           &
     &    '*** Error: An illegal solar two-stream scheme has '          &
     &    //'been selected.'
        IERR=I_ERR_FATAL
        RETURN
!
      ENDIF
!
!
!     Determine the basic solar coefficients for the
!     two-stream equations.
!
      DO I=I_LAYER_FIRST, I_LAYER_LAST
        DO L=1, N_PROFILE
          FACTOR=0.5E+00_Real64*OMEGA(L, I)*SEC_0(L)                    &
     &      /((LAMBDA(L, I)-SEC_0(L))*(LAMBDA(L, I)+SEC_0(L)))
          GAMMA_UP(L, I)=FACTOR*(SUM(L, I)-SEC_0(L)                     &
     &      -KSI_0(L, I)*(DIFF(L, I)-SEC_0(L)))
          GAMMA_DOWN(L, I)=FACTOR*(SUM(L, I)+SEC_0(L)                   &
     &      +KSI_0(L, I)*(DIFF(L, I)+SEC_0(L)))
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE SOLAR_COEFFICIENT_BASIC
#endif
