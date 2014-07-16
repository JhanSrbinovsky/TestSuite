#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to find single scattering properties of all regions.
!
! Method:
!       Straightforward.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SINGLE_SCATTERING_ALL(I_SCATTER_METHOD_BAND            &
!                     Atmospheric Properties
     &  , N_PROFILE, N_LAYER, D_MASS                                    &
!                     Cloudy Properties
     &  , L_CLOUD, N_CLOUD_TOP, N_CLOUD_TYPE                            &
!                     Optical Properties
     &  , K_GREY_TOT_CLR, K_EXT_SCAT_CLR                                &
     &  , K_GREY_TOT, K_EXT_SCAT                                        &
     &  , K_GAS_ABS                                                     &
!                     Single Scattering Properties
     &  , TAU_CLR, OMEGA_CLR                                            &
     &  , TAU, OMEGA                                                    &
!                     Dimensions of Arrays
     &  , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, ND_CLOUD_TYPE      &
     &  )
!
!
      IMPLICIT NONE
!
! Include Header Files
#include "c_kinds.h"
!
!     Sizes of arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for profiles
     &  , ND_LAYER                                                      &
!           Size allocated for layers
     &  , ND_LAYER_CLR                                                  &
!           Size allocated for completely clear layers
     &  , ID_CT                                                         &
!           Topmost declared cloudy layer
     &  , ND_CLOUD_TYPE
!           Size allocated for types of clouds
!
!
!     Dummy variables.
      INTEGER, INTENT(IN) ::                                            &
     &    I_SCATTER_METHOD_BAND
!           Treatment of scattering in the band
!
!                     Atmospheric properties
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER
!           Number of layers
      REAL  (Real64), INTENT(IN) ::                                     &
     &    D_MASS(ND_PROFILE, ND_LAYER)
!           Mass thickness of each layer
!
!                     Cloudy properties
      LOGICAL, INTENT(IN) ::                                            &
     &    L_CLOUD
!           Flag for clouds
      INTEGER, INTENT(IN) ::                                            &
     &    N_CLOUD_TOP                                                   &
!           Topmost cloudy layer
     &  , N_CLOUD_TYPE
!           Number of types of clouds
!
!                     Optical properties
      REAL  (Real64), INTENT(IN) ::                                     &
     &    K_GREY_TOT_CLR(ND_PROFILE, ND_LAYER_CLR)                      &
!           Clear-sky absorptive extinction
     &  , K_EXT_SCAT_CLR(ND_PROFILE, ND_LAYER_CLR)                      &
!           Clear-sky scattering extinction
     &  , K_GREY_TOT(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)     &
!           Absorptive extinction
     &  , K_EXT_SCAT(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)     &
!           Scattering extinction
     &  , K_GAS_ABS(ND_PROFILE, ND_LAYER)
!           Gaseous extinction
!
!                     Single scattering properties
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    TAU_CLR(ND_PROFILE, ND_LAYER_CLR)                             &
!           Clear-sky optical depth
     &  , OMEGA_CLR(ND_PROFILE, ND_LAYER_CLR)                           &
!           Clear-sky single scattering albedo
     &  , TAU(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)            &
!           Optical depth
     &  , OMEGA(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)
!           Single scattering albedo
!
!
!     Local variables.
      INTEGER                                                           &
     &    K
!           Loop variable
!
!     Subroutines called:
      EXTERNAL                                                          &
     &    SINGLE_SCATTERING
!
!
!
!     Clear-sky properties:
!
!     In the following call K_GAS_ABS can be used as if it had the
!     smaller dimension ND_LAYER_CLR as long as the last dimension
!     is over atmospheric layers.
! DEPENDS ON: single_scattering
      CALL SINGLE_SCATTERING(I_SCATTER_METHOD_BAND                      &
     &  , N_PROFILE, 1, N_CLOUD_TOP-1                                   &
     &  , D_MASS                                                        &
     &  , K_GREY_TOT_CLR, K_EXT_SCAT_CLR, K_GAS_ABS                     &
     &  , TAU_CLR, OMEGA_CLR                                            &
     &  , ND_PROFILE, ND_LAYER, 1, ND_LAYER_CLR                         &
     &  )
! DEPENDS ON: single_scattering
      CALL SINGLE_SCATTERING(I_SCATTER_METHOD_BAND                      &
     &  , N_PROFILE, N_CLOUD_TOP, N_LAYER                               &
     &  , D_MASS                                                        &
     &  , K_GREY_TOT(1, ID_CT, 0)                                       &
     &  , K_EXT_SCAT(1, ID_CT, 0)                                       &
     &  , K_GAS_ABS                                                     &
     &  , TAU(1, ID_CT, 0), OMEGA(1, ID_CT, 0)                          &
     &  , ND_PROFILE, ND_LAYER, ID_CT, ND_LAYER                         &
     &  )
!
      IF (L_CLOUD) THEN
        DO K=1, N_CLOUD_TYPE
! DEPENDS ON: single_scattering
          CALL SINGLE_SCATTERING(I_SCATTER_METHOD_BAND                  &
     &      , N_PROFILE, N_CLOUD_TOP, N_LAYER                           &
     &      , D_MASS                                                    &
     &      , K_GREY_TOT(1, ID_CT, K)                                   &
     &      , K_EXT_SCAT(1, ID_CT, K)                                   &
     &      , K_GAS_ABS                                                 &
     &      , TAU(1, ID_CT, K), OMEGA(1, ID_CT, K)                      &
     &      , ND_PROFILE, ND_LAYER, ID_CT, ND_LAYER                     &
     &      )
        ENDDO
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE SINGLE_SCATTERING_ALL
#endif
