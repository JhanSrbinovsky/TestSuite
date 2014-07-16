#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate two-stream coefficients in cloudy regions.
!
! Method:
!       The coefficients for each region are determined and
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
      SUBROUTINE TWO_COEFF_REGION_FAST_LW(IERR                          &
     &  , N_PROFILE, N_LAYER, N_CLOUD_TOP                               &
     &  , L_IR_SOURCE_QUAD, N_SOURCE_COEFF                              &
     &  , N_CLOUD_TYPE, FRAC_CLOUD                                      &
     &  , N_REGION, I_REGION_CLOUD, FRAC_REGION                         &
     &  , TAU_CLR, TAU                                                  &
     &  , ISOLIR                                                        &
     &  , TRANS, REFLECT, SOURCE_COEFF                                  &
     &  , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, ND_SOURCE_COEFF    &
     &  , ND_CLOUD_TYPE, ND_REGION                                      &
     &  )
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
     &  , ND_LAYER_CLR                                                  &
!           Maximum number of completely clear layers
     &  , ID_CT                                                         &
!           Topmost declared cloudy layer
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
#include "def_std_io_icf3z.h"
!
!
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
!     optical properties of layer:
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FRAC_CLOUD(ND_PROFILE, ND_LAYER, ND_CLOUD_TYPE)               &
!           Fractions of different types of clouds
     &  , FRAC_REGION(ND_PROFILE, ND_LAYER, ND_REGION)                  &
!           Fractions of total cloud occupied by each region
     &  , TAU_CLR(ND_PROFILE, ND_LAYER_CLR)                             &
!           Clear-sky optical depth
     &  , TAU(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)
!           Optical depth
!
!     Coefficients in the two-stream equations:
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    TRANS(ND_PROFILE, ND_LAYER, ND_REGION)                        &
!           Diffuse transmission coefficient
     &  , REFLECT(ND_PROFILE, ND_LAYER, ND_REGION)                      &
!           Diffuse reflection coefficient
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
     &    TRANS_TEMP(ND_PROFILE, ND_LAYER)                              &
!           Temporary diffuse transmission coefficient
     &  , SOURCE_COEFF_TEMP(ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF)
!           Temporary source coefficients in two-stream equations
!
!     Variables for gathering:
      INTEGER                                                           &
     &    N_LIST                                                        &
!           Number of points in list
     &  , L_LIST(ND_PROFILE)                                            &
!           List of collected points
     &  , LL
      REAL  (Real64) ::                                                 &
     &    TAU_GATHERED(ND_PROFILE, ND_LAYER)                            &
!           Gathered optical depth
     &  , TMP_INV(ND_PROFILE)
!           Temporary work array
!
!     Subroutines called:
      EXTERNAL                                                          &
     &    TWO_COEFF_FAST_LW
!
!
!
!     This routine should not be used outside the IR.
      IF (ISOLIR /= IP_INFRA_RED) THEN
        WRITE(IU_ERR, '(/A)')                                           &
     &    '*** Erroneous use of non-scattering code.'
        IERR=I_ERR_FATAL
        RETURN
      ENDIF
!
!     Determine the optical properties of the clear-sky regions of
!     the layers.
!
! DEPENDS ON: two_coeff_fast_lw
      CALL TWO_COEFF_FAST_LW(N_PROFILE, 1, N_CLOUD_TOP-1                &
     &  , L_IR_SOURCE_QUAD, TAU_CLR                                     &
     &  , TRANS(1, 1, IP_REGION_CLEAR)                                  &
     &  , SOURCE_COEFF(1, 1, 1, IP_REGION_CLEAR)                        &
     &  , ND_PROFILE, ND_LAYER, 1, ND_LAYER_CLR, ND_SOURCE_COEFF        &
     &  )
! DEPENDS ON: two_coeff_fast_lw
      CALL TWO_COEFF_FAST_LW(N_PROFILE, N_CLOUD_TOP, N_LAYER            &
     &  , L_IR_SOURCE_QUAD, TAU_CLR                                     &
     &  , TRANS(1, 1, IP_REGION_CLEAR)                                  &
     &  , SOURCE_COEFF(1, 1, 1, IP_REGION_CLEAR)                        &
     &  , ND_PROFILE, ND_LAYER, ID_CT, ND_LAYER, ND_SOURCE_COEFF        &
     &  )
      IF (IERR /= I_NORMAL) RETURN
      DO I=1, N_LAYER
        DO L=1, N_PROFILE
          REFLECT(L, I, IP_REGION_CLEAR)=0.0E+00_Real64
        ENDDO
      ENDDO
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
!           are two-dimensional.
!
            DO L=1, N_LIST
              TAU_GATHERED(L, I)=TAU(L_LIST(L), I, K)
            ENDDO
!
! DEPENDS ON: two_coeff_fast_lw
            CALL TWO_COEFF_FAST_LW(N_LIST, I, I                         &
     &        , L_IR_SOURCE_QUAD, TAU_GATHERED                          &
     &        , TRANS_TEMP                                              &
     &        , SOURCE_COEFF_TEMP                                       &
     &        , ND_PROFILE, ND_LAYER, ID_CT, ND_LAYER, ND_SOURCE_COEFF  &
     &        )
            IF (IERR /= I_NORMAL) RETURN
!
            DO LL=1, N_LIST
              L=L_LIST(LL)
              TRANS(L, I, I_REGION)=TRANS(L, I, I_REGION)               &
     &          +FRAC_CLOUD(L, I, K)*TRANS_TEMP(LL, I)
            ENDDO
            DO J=1, N_SOURCE_COEFF
              DO LL=1, N_LIST
                L=L_LIST(LL)
                SOURCE_COEFF(L, I, J, I_REGION)                         &
     &            =SOURCE_COEFF(L, I, J, I_REGION)                      &
     &            +FRAC_CLOUD(L, I, K)                                  &
     &            *SOURCE_COEFF_TEMP(LL, I, J)
              ENDDO
            ENDDO
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
            DO LL=1, N_LIST
              L=L_LIST(LL)
              TMP_INV(LL)=1.0E+00_Real64/FRAC_REGION(L, I, I_REGION)
              TRANS(L, I, I_REGION)=TRANS(L, I, I_REGION)               &
     &          *TMP_INV(LL)
            ENDDO
            DO J=1, N_SOURCE_COEFF
              DO LL=1, N_LIST
                L=L_LIST(LL)
                SOURCE_COEFF(L, I, J, I_REGION)                         &
     &            =SOURCE_COEFF(L, I, J, I_REGION)                      &
     &            *TMP_INV(LL)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE TWO_COEFF_REGION_FAST_LW
#endif
