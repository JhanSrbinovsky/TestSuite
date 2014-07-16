#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate cloudy two-stream coefficients.
!
! Method:
!       The coeffients for each type of cloud are determined and
!       averaged.
!
! Current owner of code: James Manners
!
! Description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!DD+ -------------------------------------------------------------------
! Compiler directives for specific computer systems:
!
! Indirect addressing will inhibit vectorization, but here it is safe
! because all points are independent.
!
! Fujistu VPP700:
!OCL NOVREC
!
! Cray vector machines:
!fpp$ NODEPCHK R
!
!DD- -------------------------------------------------------------------
      SUBROUTINE TWO_COEFF_CLOUD(IERR                                   &
     &   , N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST                       &
     &   , I_2STREAM, L_IR_SOURCE_QUAD, N_SOURCE_COEFF                  &
     &   , N_CLOUD_TYPE, FRAC_CLOUD                                     &
     &   , PHASE_FNC_CLOUD, OMEGA_CLOUD, TAU_CLOUD                      &
     &   , ISOLIR, SEC_0                                                &
     &   , TRANS_CLOUD, REFLECT_CLOUD, TRANS_0_CLOUD                    &
     &   , SOURCE_COEFF_CLOUD                                           &
     &   , ND_PROFILE, ND_LAYER, ID_CT, ND_MAX_ORDER                    &
     &   , ND_SOURCE_COEFF, ND_CLOUD_TYPE                               &
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
     &  , ND_LAYER                                                      &
!           Maximum number of layers
     &  , ID_CT                                                         &
!           Topmost declared potentially cloudy layer
     &  , ND_MAX_ORDER                                                  &
!           Size allocated for orders of spherical harmonics
     &  , ND_SOURCE_COEFF                                               &
!           Size allocated for source coefficients
     &  , ND_CLOUD_TYPE
!           Maximum number of types of cloud
!
!     Include header files.
#include "c_kinds.h"
#include "spectral_region_pcf3z.h"
#include "error_pcf3z.h"
!
!     Dummy arguments.
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , I_LAYER_FIRST                                                 &
!           First layer to consider
     &  , I_LAYER_LAST                                                  &
!           Last layer to consider
     &  , ISOLIR                                                        &
!           Spectral region
     &  , N_CLOUD_TYPE                                                  &
!           Number of types of clouds
     &  , I_2STREAM                                                     &
!           Two stream scheme
     &  , N_SOURCE_COEFF
!           Number of source coefficients
!
      LOGICAL, INTENT(IN) ::                                            &
     &    L_IR_SOURCE_QUAD
!           Use a quadratic source in the infra-red
!
!     Optical properties of layer:
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FRAC_CLOUD(ND_PROFILE, ID_CT: ND_LAYER, ND_CLOUD_TYPE)        &
!           Fractions of different types of clouds
     &  , PHASE_FNC_CLOUD(ND_PROFILE, ID_CT: ND_LAYER                   &
     &      , ND_MAX_ORDER, ND_CLOUD_TYPE)                              &
!           Phase functions
     &  , OMEGA_CLOUD(ND_PROFILE, ID_CT: ND_LAYER, ND_CLOUD_TYPE)       &
!           Albedo of single scattering
     &  , TAU_CLOUD(ND_PROFILE, ID_CT: ND_LAYER, ND_CLOUD_TYPE)
!           Optical depth
!
!     Solar beam
      REAL  (Real64), INTENT(IN) ::                                     &
     &    SEC_0(ND_PROFILE)
!           Secant of zenith angle
!
!
!     Coefficients in the two-stream equations:
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    TRANS_CLOUD(ND_PROFILE, ND_LAYER)                             &
!           Mean diffuse transmission coefficient
     &  , REFLECT_CLOUD(ND_PROFILE, ND_LAYER)                           &
!           Mean diffuse reflection coefficient
     &  , TRANS_0_CLOUD(ND_PROFILE, ND_LAYER)                           &
!           Mean direct transmission coefficient
     &  , SOURCE_COEFF_CLOUD(ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF)
!           Mean source coefficients in two-stream equations
!
!
!     Local variables.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , J                                                             &
!           Loop variable
     &  , K                                                             &
!           Loop variable
     &  , L
!           Loop variable
!
!     Coefficients in the two-stream equations:
      REAL  (Real64) ::                                                 &
     &    TRANS_TEMP(ND_PROFILE, 1)                                     &
!           Temporary diffuse transmission coefficient
     &  , REFLECT_TEMP(ND_PROFILE, 1)                                   &
!           Temporary diffuse reflection coefficient
     &  , TRANS_0_TEMP(ND_PROFILE, 1)                                   &
!           Temporary direct transmission coefficient
     &  , SOURCE_COEFF_TEMP(ND_PROFILE, 1, ND_SOURCE_COEFF)
!           Temporary source coefficients in two-stream equations
!
!     Variables for gathering:
      INTEGER                                                           &
     &    N_LIST                                                        &
!           Number of points in list
     &  , L_LIST(ND_PROFILE)                                            &
!           List of collected points
     &  , LL
!           Loop variable
      REAL  (Real64) ::                                                 &
     &    TAU_GATHERED(ND_PROFILE, 1)                                   &
!           Gathered optical depth
     &  , OMEGA_GATHERED(ND_PROFILE, 1)                                 &
!           Gathered alebdo of single scattering
     &  , ASYMMETRY_GATHERED(ND_PROFILE, 1)                             &
!           Gathered asymmetry
     &  , SEC_0_GATHERED(ND_PROFILE)
!           Gathered asymmetry
!
!     Subroutines called:
      EXTERNAL                                                          &
     &     TWO_COEFF
!
!
!
!     Initialize the full arrays.
!
      DO I=I_LAYER_FIRST, I_LAYER_LAST
        DO L=1, N_PROFILE
          TRANS_CLOUD(L, I)=0.0E+00_Real64
          REFLECT_CLOUD(L, I)=0.0E+00_Real64
        ENDDO
      ENDDO
      DO J=1, N_SOURCE_COEFF
        DO I=I_LAYER_FIRST, I_LAYER_LAST
          DO L=1, N_PROFILE
            SOURCE_COEFF_CLOUD(L, I, J)=0.0E+00_Real64
          ENDDO
        ENDDO
      ENDDO
!
      IF (ISOLIR == IP_SOLAR) THEN
        DO I=I_LAYER_FIRST, I_LAYER_LAST
          DO L=1, N_PROFILE
            TRANS_0_CLOUD(L, I)=0.0E+00_Real64
          ENDDO
        ENDDO
      ENDIF
!
!
!     Calculate the transmission and reflection coefficients for
!     each type of cloud and increment the totals, weighting with
!     the cloud fraction.
!
      DO K=1, N_CLOUD_TYPE
!
        DO I=I_LAYER_FIRST, I_LAYER_LAST
!
!         Determine where cloud of the current type exists
!         in this row and gather the points.
          N_LIST=0
          DO L=1, N_PROFILE
            IF (FRAC_CLOUD(L, I, K) >  0.0E+00_Real64) THEN
              N_LIST=N_LIST+1
              L_LIST(N_LIST)=L
            ENDIF
          ENDDO
!
!
          IF (N_LIST >  0) THEN
!
!           Gather the optical properties.
!           Here we must consider one layer at a time. To reduce
!           storage the temporary arrays are only one layer thick,
!           but they will be past to the subroutine where they
!           will be declared as running from the Ith to the Ith layer
!           to make the code more readable at the lower level.
!
            DO L=1, N_LIST
              TAU_GATHERED(L, 1)                                        &
     &          =TAU_CLOUD(L_LIST(L), I, K)
              OMEGA_GATHERED(L, 1)                                      &
     &          =OMEGA_CLOUD(L_LIST(L), I, K)
              ASYMMETRY_GATHERED(L, 1)                                  &
     &          =PHASE_FNC_CLOUD(L_LIST(L), I, 1, K)
            ENDDO
            IF (ISOLIR == IP_SOLAR) THEN
              DO L=1, N_LIST
                SEC_0_GATHERED(L)=SEC_0(L_LIST(L))
              ENDDO
            ENDIF
!
!
! DEPENDS ON: two_coeff
            CALL TWO_COEFF(IERR                                         &
     &        , N_LIST, I, I                                            &
     &        , I_2STREAM, L_IR_SOURCE_QUAD                             &
     &        , ASYMMETRY_GATHERED, OMEGA_GATHERED                      &
     &        , TAU_GATHERED                                            &
     &        , ISOLIR, SEC_0_GATHERED                                  &
     &        , TRANS_TEMP, REFLECT_TEMP, TRANS_0_TEMP                  &
     &        , SOURCE_COEFF_TEMP                                       &
     &        , ND_PROFILE, I, I, I, I, ND_SOURCE_COEFF                 &
     &        )
            IF (IERR /= I_NORMAL) RETURN
!
            DO L=1, N_LIST
              LL=L_LIST(L)
              TRANS_CLOUD(LL, I)=TRANS_CLOUD(LL, I)                     &
     &          +FRAC_CLOUD(LL, I, K)*TRANS_TEMP(L, 1)
              REFLECT_CLOUD(LL, I)=REFLECT_CLOUD(LL, I)                 &
     &          +FRAC_CLOUD(LL, I, K)*REFLECT_TEMP(L, 1)
            ENDDO
            DO J=1, N_SOURCE_COEFF
              DO L=1, N_LIST
                LL=L_LIST(L)
                SOURCE_COEFF_CLOUD(LL, I, J)                            &
     &            =SOURCE_COEFF_CLOUD(LL, I, J)                         &
     &            +FRAC_CLOUD(LL, I, K)*SOURCE_COEFF_TEMP(L, 1, J)
              ENDDO
            ENDDO
            IF (ISOLIR == IP_SOLAR) THEN
              DO L=1, N_LIST
                LL=L_LIST(L)
                TRANS_0_CLOUD(LL, I)=TRANS_0_CLOUD(LL, I)               &
     &             +FRAC_CLOUD(LL, I, K)*TRANS_0_TEMP(L, 1)
              ENDDO
            ENDIF
          ENDIF
!
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE TWO_COEFF_CLOUD
#endif
