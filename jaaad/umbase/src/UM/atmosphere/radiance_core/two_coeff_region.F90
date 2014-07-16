#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate two-stream coefficients in the regions.
!
! Method:
!       The coeffients for each region are determined and
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
      SUBROUTINE TWO_COEFF_REGION(IERR                                  &
     &   , N_PROFILE, N_LAYER, N_CLOUD_TOP                              &
     &   , I_2STREAM, L_IR_SOURCE_QUAD, N_SOURCE_COEFF                  &
     &   , N_CLOUD_TYPE, FRAC_CLOUD                                     &
     &   , N_REGION, I_REGION_CLOUD, FRAC_REGION                        &
     &   , PHASE_FNC_CLR, OMEGA_CLR, TAU_CLR                            &
     &   , PHASE_FNC, OMEGA, TAU                                        &
     &   , ISOLIR, SEC_0                                                &
     &   , TRANS, REFLECT, TRANS_0                                      &
     &   , SOURCE_COEFF                                                 &
     &   , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT                    &
     &   , ND_MAX_ORDER, ND_SOURCE_COEFF                                &
     &   , ND_CLOUD_TYPE, ND_REGION                                     &
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
     &  , ND_LAYER_CLR                                                  &
!           Size allocated for completely clear atmospheric layers
     &  , ID_CT                                                         &
!           Topmost declared cloudy layer
     &  , ND_MAX_ORDER                                                  &
!           Size allocated for orders of spherical harmonics
     &  , ND_SOURCE_COEFF                                               &
!           Size allocated for source coefficients
     &  , ND_CLOUD_TYPE                                                 &
!           Maximum number of types of cloud
     &  , ND_REGION
!           Maximum number of cloudy regions
!
!     Include header files.
#include "c_kinds.h"
#include "spectral_region_pcf3z.h"
#include "error_pcf3z.h"
#include "cloud_region_pcf3z.h"
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
     &  , N_CLOUD_TOP                                                   &
!           Topmost cloudy layer
     &  , ISOLIR                                                        &
!           Spectral region
     &  , N_CLOUD_TYPE                                                  &
!           Number of types of clouds
     &  , I_2STREAM                                                     &
!           Two stream scheme
     &  , N_SOURCE_COEFF
!           Number of source coefficients
!
      INTEGER, INTENT(IN) ::                                            &
     &    N_REGION                                                      &
!           Number of cloudy regions
     &  , I_REGION_CLOUD(ND_CLOUD_TYPE)
!           Regions in which types of clouds fall
!
      LOGICAL, INTENT(IN) ::                                            &
     &    L_IR_SOURCE_QUAD
!           Use a quadratic source in the infra-red
!
!     Optical properties of layer:
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FRAC_CLOUD(ND_PROFILE, ID_CT: ND_LAYER, ND_CLOUD_TYPE)        &
!           Fractions of different types of clouds
     &  , FRAC_REGION(ND_PROFILE, ID_CT: ND_LAYER, ND_REGION)           &
!           Fractions of total cloud occupied by each region
     &  , PHASE_FNC_CLR(ND_PROFILE, ND_LAYER_CLR, ND_MAX_ORDER)         &
!           Phase function in clear-sky
     &  , OMEGA_CLR(ND_PROFILE, ND_LAYER_CLR)                           &
!           Clear-sky albedo of single scattering
     &  , TAU_CLR(ND_PROFILE, ND_LAYER_CLR)                             &
!           Clear-sky optical depth
     &  , TAU(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)            &
!           Optical depth
     &  , OMEGA(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)          &
!           Albedo of single scattering
     &  , PHASE_FNC(ND_PROFILE, ID_CT: ND_LAYER                         &
     &      , ND_MAX_ORDER, 0: ND_CLOUD_TYPE)
!           Phase function
!
!     Solar beam
      REAL  (Real64), INTENT(IN) ::                                     &
     &    SEC_0(ND_PROFILE)
!           Secant of zenith angle
!
!
!     Coefficients in the two-stream equations:
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    TRANS(ND_PROFILE, ND_LAYER, ND_REGION)                        &
!           Diffuse transmission coefficient
     &  , REFLECT(ND_PROFILE, ND_LAYER, ND_REGION)                      &
!           Diffuse reflection coefficient
     &  , TRANS_0(ND_PROFILE, ND_LAYER, ND_REGION)                      &
!           Direct transmission coefficient
     &  , SOURCE_COEFF(ND_PROFILE, ND_LAYER                             &
     &    , ND_SOURCE_COEFF, ND_REGION)
!           Source coefficients in two-stream equations
!
!     Local variables.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , J                                                             &
!           Loop variable
     &  , K                                                             &
!           Loop variable
     &  , L                                                             &
!           Loop variable
     &  , I_REGION
!           Loop variable over regions
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
!     Determine the optical properties of the clear-sky regions of
!     the layers.
!
! DEPENDS ON: two_coeff
      CALL TWO_COEFF(IERR                                               &
     &  , N_PROFILE, 1, N_CLOUD_TOP-1                                   &
     &  , I_2STREAM, L_IR_SOURCE_QUAD                                   &
     &  , PHASE_FNC_CLR(1, 1, 1), OMEGA_CLR, TAU_CLR                    &
     &  , ISOLIR, SEC_0                                                 &
     &  , TRANS(1, 1, IP_REGION_CLEAR)                                  &
     &  , REFLECT(1, 1, IP_REGION_CLEAR)                                &
     &  , TRANS_0(1, 1, IP_REGION_CLEAR)                                &
     &  , SOURCE_COEFF(1, 1, 1, IP_REGION_CLEAR)                        &
     &  , ND_PROFILE, 1, ND_LAYER_CLR, 1, ND_LAYER, ND_SOURCE_COEFF     &
     &  )
! DEPENDS ON: two_coeff
      CALL TWO_COEFF(IERR                                               &
     &  , N_PROFILE, N_CLOUD_TOP, N_LAYER                               &
     &  , I_2STREAM, L_IR_SOURCE_QUAD                                   &
     &  , PHASE_FNC(1, ID_CT, 1, 0)                                     &
     &  , OMEGA(1, ID_CT, 0), TAU(1, ID_CT, 0)                          &
     &  , ISOLIR, SEC_0                                                 &
     &  , TRANS(1, 1, IP_REGION_CLEAR)                                  &
     &  , REFLECT(1, 1, IP_REGION_CLEAR)                                &
     &  , TRANS_0(1, 1, IP_REGION_CLEAR)                                &
     &  , SOURCE_COEFF(1, 1, 1, IP_REGION_CLEAR)                        &
     &  , ND_PROFILE, ID_CT, ND_LAYER, 1, ND_LAYER                      &
     &  , ND_SOURCE_COEFF                                               &
     &  )
      IF (IERR /= I_NORMAL) RETURN
!
!
!     Now deal with clouds.
!
!     Initialize the full arrays for cloudy regions.
!
      DO I_REGION=1, N_REGION
        IF (I_REGION /= IP_REGION_CLEAR) THEN
          DO I=N_CLOUD_TOP, N_LAYER
            DO L=1, N_PROFILE
              TRANS(L, I, I_REGION)=0.0E+00_Real64
              REFLECT(L, I, I_REGION)=0.0E+00_Real64
            ENDDO
          ENDDO
          DO J=1, N_SOURCE_COEFF
            DO I=N_CLOUD_TOP, N_LAYER
              DO L=1, N_PROFILE
                SOURCE_COEFF(L, I, J, I_REGION)=0.0E+00_Real64
              ENDDO
            ENDDO
          ENDDO
!
          IF (ISOLIR == IP_SOLAR) THEN
            DO I=N_CLOUD_TOP, N_LAYER
              DO L=1, N_PROFILE
                TRANS_0(L, I, I_REGION)=0.0E+00_Real64
              ENDDO
            ENDDO
          ENDIF
!
        ENDIF
!
      ENDDO
!
!
!
!     Consider each type of cloud in turn, checking which region it
!     contrubutes to and form weighted sums of cloud properties.
!
      DO K=1, N_CLOUD_TYPE
!
!
!       Set the region in which clouds of this type are included.
        I_REGION=I_REGION_CLOUD(K)
!
        DO I=N_CLOUD_TOP, N_LAYER
!
!         Form a list of points where cloud of this type exists
!         on this row for gathering.
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
!           Gather the optical properties. Though we consider only
!           one layer at a time the lower routines will operate on
!           arrays with vertical structure, so the gathered arrays
!           are two-dimensional, however, it is only necessary to
!           have one layer in the temporary arrays.
!
            DO L=1, N_LIST
              TAU_GATHERED(L, 1)                                        &
     &          =TAU(L_LIST(L), I, K)
              OMEGA_GATHERED(L, 1)                                      &
     &          =OMEGA(L_LIST(L), I, K)
              ASYMMETRY_GATHERED(L, 1)                                  &
     &          =PHASE_FNC(L_LIST(L), I, 1, K)
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
!
            DO L=1, N_LIST
              LL=L_LIST(L)
              TRANS(LL, I, I_REGION)=TRANS(LL, I, I_REGION)             &
     &          +FRAC_CLOUD(LL, I, K)*TRANS_TEMP(L, 1)
              REFLECT(LL, I, I_REGION)=REFLECT(LL, I, I_REGION)         &
     &          +FRAC_CLOUD(LL, I, K)*REFLECT_TEMP(L, 1)
            ENDDO
            DO J=1, N_SOURCE_COEFF
              DO L=1, N_LIST
                LL=L_LIST(L)
                SOURCE_COEFF(LL, I, J, I_REGION)                        &
     &            =SOURCE_COEFF(LL, I, J, I_REGION)                     &
     &            +FRAC_CLOUD(LL, I, K)                                 &
     &            *SOURCE_COEFF_TEMP(L, 1, J)
              ENDDO
            ENDDO
            IF (ISOLIR == IP_SOLAR) THEN
              DO L=1, N_LIST
                LL=L_LIST(L)
                TRANS_0(LL, I, I_REGION)=TRANS_0(LL, I, I_REGION)       &
     &            +FRAC_CLOUD(LL, I, K)*TRANS_0_TEMP(L, 1)
              ENDDO
            ENDIF
!
          ENDIF
!
        ENDDO
      ENDDO
!
!
!     Finally, scale the weighted sums by the cloud fractions.
      DO I_REGION=1, N_REGION
        IF (I_REGION /= IP_REGION_CLEAR) THEN
          DO I=N_CLOUD_TOP, N_LAYER
!
!           Gather points within this region.
            N_LIST=0
            DO L=1,N_PROFILE
              IF (FRAC_REGION(L, I, I_REGION) >  0.0E+00_Real64) THEN
                N_LIST=N_LIST+1
                L_LIST(N_LIST)=L
              ENDIF
            ENDDO
            DO L=1, N_LIST
              LL=L_LIST(L)
              TRANS(LL, I, I_REGION)=TRANS(LL, I, I_REGION)             &
     &          /FRAC_REGION(LL, I, I_REGION)
              REFLECT(LL, I, I_REGION)=REFLECT(LL, I, I_REGION)         &
     &          /FRAC_REGION(LL, I, I_REGION)
            ENDDO
            DO J=1, N_SOURCE_COEFF
              DO L=1, N_LIST
                LL=L_LIST(L)
                SOURCE_COEFF(LL, I, J, I_REGION)                        &
     &            =SOURCE_COEFF(LL, I, J, I_REGION)                     &
     &            /FRAC_REGION(LL, I, I_REGION)
              ENDDO
            ENDDO
            IF (ISOLIR == IP_SOLAR) THEN
              DO L=1, N_LIST
                LL=L_LIST(L)
                TRANS_0(LL, I, I_REGION)=TRANS_0(LL, I, I_REGION)       &
     &            /FRAC_REGION(LL, I, I_REGION)
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE TWO_COEFF_REGION
#endif
