#if defined(A70_1B) || defined(A70_1C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to initialize diagnostics and coupling arrays.
!
! Purpose:
!   The coupling and diagnostic arrays are zeroed.
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
      SUBROUTINE R2_INIT_COUPLE_DIAG(N_PROFILE                          &
        , SEA_FLUX                                                      &
        , L_SURFACE_DOWN_FLUX, SURFACE_DOWN_FLUX                        &
        , L_SURF_DOWN_CLR, SURF_DOWN_CLR                                &
        , L_SURF_UP_CLR, SURF_UP_CLR                                    &
        , L_FLUX_DIFFUSE, FLUX_DIFFUSE                                  &
        , L_UVFLUX_DIRECT, UV_FLUX_DIRECT                               &
        , L_UVFLUX_UP, UV_FLUX_UP                                       &
        , L_UVFLUX_DOWN, UV_FLUX_DOWN                                   &
        , L_FLUX_BELOW_690NM_SURF                                       &
        , FLUX_BELOW_690NM_SURF, FL_SEA_BELOW_690NM_SURF                &
        , L_MOSES_II                                                    &
        , SURF_VIS_DIR, SURF_VIS_DIF, SURF_NIR_DIR, SURF_NIR_DIF        &
        , L_CLOUD_EXTINCTION, CLOUD_EXTINCTION                          &
        , CLOUD_WEIGHT_EXTINCTION                                       &
        , L_LS_CLOUD_EXTINCTION, LS_CLOUD_EXTINCTION                    &
        , LS_CLOUD_WEIGHT_EXTINCTION                                    &
        , L_CNV_CLOUD_EXTINCTION, CNV_CLOUD_EXTINCTION                  &
        , CNV_CLOUD_WEIGHT_EXTINCTION                                   &
        , L_CLOUD_ABSORPTIVITY, CLOUD_ABSORPTIVITY                      &
        , CLOUD_WEIGHT_ABSORPTIVITY                                     &
        , L_LS_CLOUD_ABSORPTIVITY, LS_CLOUD_ABSORPTIVITY                &
        , LS_CLOUD_WEIGHT_ABSORPTIVITY                                  &
        , L_CNV_CLOUD_ABSORPTIVITY, CNV_CLOUD_ABSORPTIVITY              &
        , CNV_CLOUD_WEIGHT_ABSORPTIVITY                                 &
        , NPD_PROFILE, NPD_LAYER                                        &
        )
!
!
!
      IMPLICIT NONE
!
!
!     Dummy arguments
!
!     Dimensions of arrays
      INTEGER                                                           &
                !, INTENT(IN)
          NPD_PROFILE                                                   &
!           Maximum number of atmospheric profiles
        , NPD_LAYER
!             Maximum mumber of layers
!
      INTEGER                                                           &
                !, INTENT(IN)
          N_PROFILE
!           Number of atmospheric profiles
!
!     Switches for diagnostics:
      LOGICAL                                                           &
                !, INTENT(IN)
          L_FLUX_BELOW_690NM_SURF                                       &
!           Flux below 690nm at surface to be calculated
        , L_MOSES_II                                                    &
!           Surface SW fluxes required for MOSES II
        , L_SURFACE_DOWN_FLUX                                           &
!           Downward surface flux required
        , L_SURF_DOWN_CLR                                               &
!           Calculate downward clear flux
        , L_SURF_UP_CLR                                                 &
!           Calculate upward clear flux
        , L_FLUX_DIFFUSE                                                &  
!           Calculate the diffuse downward flux
        , L_UVFLUX_DIRECT                                               &
!           Calculate the direct UV flux
        , L_UVFLUX_UP                                                   &
!           Calculate the upward UV flux
        , L_UVFLUX_DOWN
!           Calculate the downward UV Flux           
!
      LOGICAL                                                           &
                !, INTENT(IN)
          L_CLOUD_EXTINCTION                                            &
!           Flag to calculate cloudy extinction
        , L_CLOUD_ABSORPTIVITY                                          &
!           Flag to calculate cloudy absorptivity
        , L_LS_CLOUD_EXTINCTION                                         &
!           Flag to calculate cloudy extinction
        , L_LS_CLOUD_ABSORPTIVITY                                       &
!           Flag to calculate cloudy absorptivity
        , L_CNV_CLOUD_EXTINCTION                                        &
!           Flag to calculate cloudy extinction
        , L_CNV_CLOUD_ABSORPTIVITY
!           Flag to calculate cloudy absorptivity
!
!     Surface fluxes for coupling or diagnostic use
      REAL                                                              &
                !, INTENT(OUT)
          SEA_FLUX(NPD_PROFILE)                                         &
!           Net downward flux into sea
        , SURFACE_DOWN_FLUX(NPD_PROFILE)                                &
!           Downward flux at surface
        , SURF_DOWN_CLR(NPD_PROFILE)                                    &
!           Clear-sky downward flux at surface
        , SURF_UP_CLR(NPD_PROFILE)                                      &
!           Clear-sky upward flux at surface
        , FLUX_DIFFUSE(NPD_PROFILE, 0:NPD_LAYER)                        &
!           Diffuse downward flux.   
        , UV_FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)                     &
!             Direct UV Flux
        , UV_FLUX_UP(NPD_PROFILE, 0: NPD_LAYER)                         &
!             Upward UV Flux
        , UV_FLUX_DOWN(NPD_PROFILE, 0: NPD_LAYER)                       &
!             Downward UV Flux  
        , FLUX_BELOW_690NM_SURF(NPD_PROFILE)                            &
!          GRID BOX MEAN NET SURFACE FLUX BELOW 690NM
        , FL_SEA_BELOW_690NM_SURF(NPD_PROFILE)                          &
!          OPEN SEA NET SURFACE FLUX BELOW 690NM
        , SURF_VIS_DIR(NPD_PROFILE)                                     &
!           Downward surface direct beam visible flux
        , SURF_VIS_DIF(NPD_PROFILE)                                     &
!           Downward surface diffuse visible flux
        , SURF_NIR_DIR(NPD_PROFILE)                                     &
!           Downward surface direct beam near-infrared flux
        , SURF_NIR_DIF(NPD_PROFILE)
!           Downward surface diffuse near-infrared flux
!
!
      REAL                                                              &
                !, INTENT(OUT)
          CLOUD_EXTINCTION(NPD_PROFILE, NPD_LAYER)                      &
!             overall extinction of clouds
        , CLOUD_WEIGHT_EXTINCTION(NPD_PROFILE, NPD_LAYER)               &
!             weighting for cloud extinction
        , LS_CLOUD_EXTINCTION(NPD_PROFILE, NPD_LAYER)                   &
!             extinction of LS clouds
        , LS_CLOUD_WEIGHT_EXTINCTION(NPD_PROFILE, NPD_LAYER)            &
!             weighting for LS cloud extinction
        , CNV_CLOUD_EXTINCTION(NPD_PROFILE, NPD_LAYER)                  &
!             extinction of conv clouds
        , CNV_CLOUD_WEIGHT_EXTINCTION(NPD_PROFILE, NPD_LAYER)           &
!             weighting for conv cloud extinction
        , CLOUD_ABSORPTIVITY(NPD_PROFILE, NPD_LAYER)                    &
!             overall absorptivity of clouds
        , CLOUD_WEIGHT_ABSORPTIVITY(NPD_PROFILE, NPD_LAYER)             &
!             weighting for cloud absorptivity
        , LS_CLOUD_ABSORPTIVITY(NPD_PROFILE, NPD_LAYER)                 &
!             absorptivity of LS clouds
        , LS_CLOUD_WEIGHT_ABSORPTIVITY(NPD_PROFILE, NPD_LAYER)          &
!             weighting for LS cloud absorptivity
        , CNV_CLOUD_ABSORPTIVITY(NPD_PROFILE, NPD_LAYER)                &
!             absorptivity of conv clouds
        , CNV_CLOUD_WEIGHT_ABSORPTIVITY(NPD_PROFILE, NPD_LAYER)
!             weighting for conv cloud absorptivity

      SEA_FLUX = 0.0
!
      IF (L_SURFACE_DOWN_FLUX) SURFACE_DOWN_FLUX = 0.0
      IF (L_SURF_DOWN_CLR)     SURF_DOWN_CLR     = 0.0
      IF (L_SURF_UP_CLR)       SURF_UP_CLR       = 0.0
      IF (L_FLUX_DIFFUSE)      FLUX_DIFFUSE      = 0.0
      IF (L_UVFLUX_DIRECT)     UV_FLUX_DIRECT    = 0.0
      IF (L_UVFLUX_UP)         UV_FLUX_UP        = 0.0
      IF (L_UVFLUX_DOWN)       UV_FLUX_DOWN      = 0.0

!
      IF (L_FLUX_BELOW_690NM_SURF) THEN
         FLUX_BELOW_690NM_SURF   = 0.0
         FL_SEA_BELOW_690NM_SURF = 0.0
      ENDIF
!
      IF (L_MOSES_II) THEN
         SURF_VIS_DIR=0.0
         SURF_VIS_DIF=0.0
         SURF_NIR_DIR=0.0
         SURF_NIR_DIF=0.0
      ENDIF
!
      IF (L_CLOUD_EXTINCTION) THEN
         CLOUD_EXTINCTION        = 0.0
         CLOUD_WEIGHT_EXTINCTION = 0.0
      ENDIF
!
      IF (L_LS_CLOUD_EXTINCTION) THEN
         LS_CLOUD_EXTINCTION        = 0.0
         LS_CLOUD_WEIGHT_EXTINCTION = 0.0
      ENDIF
!
      IF (L_CNV_CLOUD_EXTINCTION) THEN 
         CNV_CLOUD_EXTINCTION        = 0.0
         CNV_CLOUD_WEIGHT_EXTINCTION = 0.0
      ENDIF
!
      IF (L_CLOUD_ABSORPTIVITY) THEN
         CLOUD_ABSORPTIVITY        = 0.0
         CLOUD_WEIGHT_ABSORPTIVITY = 0.0
      ENDIF
!
      IF (L_LS_CLOUD_ABSORPTIVITY) THEN
         LS_CLOUD_ABSORPTIVITY        = 0.0
         LS_CLOUD_WEIGHT_ABSORPTIVITY = 0.0
      ENDIF
!
      IF (L_CNV_CLOUD_ABSORPTIVITY) THEN
         CNV_CLOUD_ABSORPTIVITY        = 0.0
         CNV_CLOUD_WEIGHT_ABSORPTIVITY = 0.0
      ENDIF
!
      RETURN
      END SUBROUTINE R2_INIT_COUPLE_DIAG
#endif
