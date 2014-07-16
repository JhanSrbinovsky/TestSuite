


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to find single scattering propeties of all regions.
!
! Method:
!       Straightforward.
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
      SUBROUTINE SINGLE_SCATTERING_ALL(I_SCATTER_METHOD_BAND            &
!                       Atmospheric Propeties
     &   , N_PROFILE, N_LAYER, D_MASS                                   &
!                       Cloudy Propeties
     &   , L_CLOUD, N_CLOUD_TOP, N_CLOUD_TYPE                           &
!                       Optical Propeties
     &   , K_GREY_TOT_FREE, K_EXT_SCAT_FREE                             &
     &   , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD                           &
     &   , K_GAS_ABS                                                    &
!                       Single Scattering Propeties
     &   , TAU_FREE, OMEGA_FREE                                         &
     &   , TAU_CLOUD, OMEGA_CLOUD                                       &
!                       Dimensions of Arrays
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
!
      IMPLICIT NONE
!
!
!     SIZES OF ARRAYS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     INCLUDE COMDECKS.
! DIMFIX3A defines internal dimensions tied to algorithms for
! two-stream radiation code, mostly for clouds

      ! number of components of clouds
      INTEGER,PARAMETER:: NPD_CLOUD_COMPONENT=4

      ! number of permitted types of clouds.
      INTEGER,PARAMETER:: NPD_CLOUD_TYPE=4

      ! number of permitted representations of clouds.
      INTEGER,PARAMETER:: NPD_CLOUD_REPRESENTATION=4

      ! number of overlap coefficients for clouds
      INTEGER,PARAMETER:: NPD_OVERLAP_COEFF=18

      ! number of coefficients for two-stream sources
      INTEGER,PARAMETER:: NPD_SOURCE_COEFF=2

      ! number of regions in a layer
      INTEGER,PARAMETER:: NPD_REGION=3

! DIMFIX3A end
!
!
!     DUMMY VARIABLES.
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_SCATTER_METHOD_BAND
!             TREATMENT OF SCATTERING IN THE BAND
!
!                       Atmospheric Properties
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
      REAL                                                              &
                !, INTENT(IN)
     &     D_MASS(NPD_PROFILE, NPD_LAYER)
!             MASS THICKNESS OF EACH LAYER
!
!                       Cldoudy Propeties
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_CLOUD
!             FLAG FOR CLOUDS
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_CLOUD_TOP                                                  &
!             TOPMOST CLOUDY LAYER
     &   , N_CLOUD_TYPE
!             NUMBER OF TYPES OF CLOUDS
!
!                       Optical Properties
      REAL                                                              &
                !, INTENT(IN)
     &     K_GREY_TOT_FREE(NPD_PROFILE, NPD_LAYER)                      &
!             FREE ABSORPTIVE EXTINCTION
     &   , K_EXT_SCAT_FREE(NPD_PROFILE, NPD_LAYER)                      &
!             FREE SCATTERING EXTINCTION
     &   , K_GREY_TOT_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)     &
!             STRATIFORM ABSORPTIVE EXTINCTION
     &   , K_EXT_SCAT_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)     &
!             STRATIFORM SCATTERING EXTINCTION
     &   , K_GAS_ABS(NPD_PROFILE, NPD_LAYER)
!             GASEOUS EXTINCTION
!
!                       Single Scattering Properties
      REAL                                                              &
                !, INTENT(OUT)
     &     TAU_FREE(NPD_PROFILE, NPD_LAYER)                             &
!             FREE OPTICAL DEPTH
     &   , OMEGA_FREE(NPD_PROFILE, NPD_LAYER)                           &
!             FREE SINGLE SCATTERING ALBEDO
     &   , TAU_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)            &
!             CLOUDY OPTICAL DEPTH
     &   , OMEGA_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY SINGLE SCATTERING ALBEDO
!
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     K
!             LOOP VARIABLE
!
!     SUBROUTINES CALLED:
      EXTERNAL                                                          &
     &     SINGLE_SCATTERING
!
!
!
!     CLEAR-SKY PROPERTIES:
!
! DEPENDS ON: single_scattering
      CALL SINGLE_SCATTERING(I_SCATTER_METHOD_BAND                      &
     &   , N_PROFILE, N_LAYER, 1                                        &
     &   , D_MASS                                                       &
     &   , K_GREY_TOT_FREE, K_EXT_SCAT_FREE, K_GAS_ABS                  &
     &   , TAU_FREE, OMEGA_FREE                                         &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
      IF (L_CLOUD) THEN
         DO K=1, N_CLOUD_TYPE
! DEPENDS ON: single_scattering
            CALL SINGLE_SCATTERING(I_SCATTER_METHOD_BAND                &
     &         , N_PROFILE, N_LAYER, N_CLOUD_TOP                        &
     &         , D_MASS                                                 &
     &         , K_GREY_TOT_CLOUD(1, 1, K)                              &
     &         , K_EXT_SCAT_CLOUD(1, 1, K)                              &
     &         , K_GAS_ABS                                              &
     &         , TAU_CLOUD(1, 1, K), OMEGA_CLOUD(1, 1, K)               &
     &         , NPD_PROFILE, NPD_LAYER                                 &
     &         )
         ENDDO
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE SINGLE_SCATTERING_ALL
