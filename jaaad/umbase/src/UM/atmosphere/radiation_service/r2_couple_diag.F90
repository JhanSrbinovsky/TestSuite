#if defined(A70_1B) || defined(A70_1C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate spectral diagnostics and coupling arrays.
!
! Purpose:
!   The coupling and diagnostic arrays are calculated.
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
      SUBROUTINE R2_COUPLE_DIAG(N_PROFILE, N_LAYER, ISOLIR              &
        , ALBEDO_FIELD_DIFF, ALBEDO_FIELD_DIR                           &
        , ALBEDO_SEA_DIFF, ALBEDO_SEA_DIR                               &
        , FLANDG, ICE_FRACTION                                          &
        , PLANCK_FREEZE_SEA, PLANCK_LEADS_SEA                           &
        , PLANCK_AIR_SURFACE, THERMAL_SOURCE_GROUND                     &
        , FLUX_DOWN, FLUX_UP, FLUX_DIRECT                               &
        , FLUX_DOWN_CLEAR, FLUX_UP_CLEAR, FLUX_DIRECT_CLEAR             &
        , WEIGHT_690NM, WEIGHT_UV                                       &
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
        , L_MOSES_II, L_CTILE                                           &
        , SURF_VIS_DIR, SURF_VIS_DIF, SURF_NIR_DIR, SURF_NIR_DIF        &
        , L_CLOUD_EXTINCTION, CLOUD_EXTINCTION_BAND                     &
        , L_LS_CLOUD_EXTINCTION, LS_CLOUD_EXTINCTION_BAND               &
        , L_CNV_CLOUD_EXTINCTION, CNV_CLOUD_EXTINCTION_BAND             &
        , FLUX_DIRECT_BAND, FLUX_DIRECT_CLEAR_BAND                      &
        , N_CLOUD_TOP, N_CLOUD_PROFILE, I_CLOUD_PROFILE                 &
        , W_CLOUD, FRAC_CLOUD                                           &
        , CLOUD_EXTINCTION, CLOUD_WEIGHT_EXTINCTION                     &
        , LS_CLOUD_EXTINCTION, LS_CLOUD_WEIGHT_EXTINCTION               &
        , CNV_CLOUD_EXTINCTION, CNV_CLOUD_WEIGHT_EXTINCTION             &
        , L_CLOUD_ABSORPTIVITY, CLOUD_ABSORPTIVITY_BAND                 &
        , L_LS_CLOUD_ABSORPTIVITY, LS_CLOUD_ABSORPTIVITY_BAND           &
        , L_CNV_CLOUD_ABSORPTIVITY, CNV_CLOUD_ABSORPTIVITY_BAND         &
        , FLUX_TOTAL_BAND, FLUX_TOTAL_CLEAR_BAND                        &
        , CLOUD_ABSORPTIVITY, CLOUD_WEIGHT_ABSORPTIVITY                 &
        , LS_CLOUD_ABSORPTIVITY, LS_CLOUD_WEIGHT_ABSORPTIVITY           &
        , CNV_CLOUD_ABSORPTIVITY, CNV_CLOUD_WEIGHT_ABSORPTIVITY         &
        , NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE                        &
        )
!
!
!
      IMPLICIT NONE
!
!
!     Comdecks included
!     Spectral regions
#include "spcrg3a.h"
!
!     Dummy Arguments
!
!     Dimensions of arrays
      INTEGER                                                           &
                !, INTENT(IN)
          NPD_PROFILE                                                   &
!           Maximum number of atmospheric profiles
        , NPD_LAYER                                                     &
!           Maximum number of atmospheric layers
        , NPD_CLOUD_TYPE
!           Maximum number of cloud types
!
      INTEGER                                                           &
                !, INTENT(IN)
          N_PROFILE                                                     &
!           Number of atmospheric profiles
        , N_LAYER                                                       &
!           Number of layers in the atmosphere
        , ISOLIR
!           Spectral region
!
!     Logical switches for the code
      LOGICAL                                                           &
                !, INTENT(IN)
          L_NET
!           Flag for net fluxes
!
!     Switches for diagnostics:
      LOGICAL                                                           &
                !, INTENT(IN)
          L_FLUX_BELOW_690NM_SURF                                       &
!           Flux below 690nm at surface to be calculated
        , L_MOSES_II                                                    &
!           Surface SW required for MOSES II
        , L_CTILE                                                       &
!           Switch for coastal tiling
        , L_SURFACE_DOWN_FLUX                                           &
!           Downward surface flux required
        , L_SURF_DOWN_CLR                                               &
!           Calculate downward clear flux
        , L_SURF_UP_CLR                                                 &
!           Calculate upward clear flux
        , L_FLUX_DIFFUSE                                                &
!           Calculate the diffuse downward flux 
        , L_UVFLUX_DIRECT                                               &
!           Calculate the Direct UV Flux 
        , L_UVFLUX_UP                                                   &
!           Calculate the Upward UV Flux
        , L_UVFLUX_DOWN                                                 &
!           Calculate the Downward UV Flux
        , L_CLOUD_EXTINCTION                                            &
!           Calculate cloud extinction
        , L_CLOUD_ABSORPTIVITY                                          &
!           Calculate cloud absorptivity
        , L_LS_CLOUD_EXTINCTION                                         &
!           Calculate cloud extinction
        , L_LS_CLOUD_ABSORPTIVITY                                       &
!           Calculate cloud absorptivity
        , L_CNV_CLOUD_EXTINCTION                                        &
!           Calculate cloud extinction
        , L_CNV_CLOUD_ABSORPTIVITY
!           Calculate cloud absorptivity
!
!     Albedos
      REAL                                                              &
                !, INTENT(IN)
          ALBEDO_FIELD_DIFF(NPD_PROFILE)                                &
!           Diffuse albedo meaned over grid box
        , ALBEDO_FIELD_DIR(NPD_PROFILE)                                 &
!           Direct albedo meaned over grid box
        , ALBEDO_SEA_DIFF(NPD_PROFILE)                                  &
!           Diffuse albedo of open sea
        , ALBEDO_SEA_DIR(NPD_PROFILE)
!           Direct albedo meaned of open sea
!
      REAL                                                              &
                !, INTENT(IN)
          THERMAL_SOURCE_GROUND(NPD_PROFILE)                            &
!           Thermal source at ground
        , PLANCK_AIR_SURFACE(NPD_PROFILE)
!           Planck function at near-surface air temperature in band
!
!     Arguments relating to sea ice.
      REAL                                                              &
            !, INTENT(IN)
          PLANCK_FREEZE_SEA                                             &
!           Planck function over freezing sea
         , PLANCK_LEADS_SEA(NPD_PROFILE)                                &
!           Planck function over sea leads
         , FLANDG(NPD_PROFILE)                                          &
!            Land fraction
         , ICE_FRACTION(NPD_PROFILE)
!            FRACTION OF SEA-ICE IN SEA PORTION OF GRID BOX!
!
      REAL                                                              &
                !, INTENT(IN)
          WEIGHT_690NM                                                  &
!           Weighting applied to band for region below 690 nm
        , WEIGHT_UV
!           Weighting for UV Fluxes.
!
!     Calculated fluxes
      REAL                                                              &
                !, INTENT(IN)
          FLUX_DOWN(NPD_PROFILE)                                        &
!           Total downward or net flux at surface
        , FLUX_DIRECT(NPD_PROFILE)                                      &
!           Direct solar flux at surface
        , FLUX_UP(NPD_PROFILE)                                          &
!           Upward flux at surface
        , FLUX_DOWN_CLEAR(NPD_PROFILE)                                  &
!           Total clear-sky downward or net flux at surface
        , FLUX_UP_CLEAR(NPD_PROFILE)                                    &
!           Clear-sky upward flux at surface
        , FLUX_DIFFUSE(NPD_PROFILE, 0:NPD_LAYER)                        &
!           The diffuse downward flux 
        , FLUX_DIRECT_CLEAR(NPD_PROFILE)
!           Clear-sky direct solar flux at surface
!
! UV Fluxes:
!
       REAL                                                             &
            !, INTENT(OUT)
           UV_FLUX_DIRECT(NPD_PROFILE, 0:NPD_LAYER)                     &
!              DIRECT UV_FLUX
         , UV_FLUX_UP(NPD_PROFILE, 0:NPD_LAYER)                         &
!              UPWARD UV-FLUX
         , UV_FLUX_DOWN(NPD_PROFILE, 0:NPD_LAYER)
!              DOWNWARD UV-FLUX
!
!     Properties of clouds:
!
      INTEGER                                                           &
                !, INTENT(IN)
          N_CLOUD_TOP                                                   &
!             Topmost cloudy layer
        , N_CLOUD_PROFILE(NPD_LAYER)                                    &
!             Number of cloudy profiles in each layer
        , I_CLOUD_PROFILE(NPD_PROFILE, NPD_LAYER)
!             Profiles containing clouds
!
      REAL                                                              &
                !, INTENT(IN)
          CLOUD_EXTINCTION_BAND(NPD_PROFILE, NPD_LAYER)                 &
!             Extinction of cloud in the current band
        , CLOUD_ABSORPTIVITY_BAND(NPD_PROFILE, NPD_LAYER)               &
!             Absorptivity of cloud in the current band
        , LS_CLOUD_EXTINCTION_BAND(NPD_PROFILE, NPD_LAYER)              &
!             Extinction of cloud in the current band
        , LS_CLOUD_ABSORPTIVITY_BAND(NPD_PROFILE, NPD_LAYER)            &
!             ABSORPTIVITY OF CLOUD IN THE CURRENT BAND
        , CNV_CLOUD_EXTINCTION_BAND(NPD_PROFILE, NPD_LAYER)             &
!             Extinction of cloud in the current band
        , CNV_CLOUD_ABSORPTIVITY_BAND(NPD_PROFILE, NPD_LAYER)           &
!             Absorptivity of cloud in the current band
        , W_CLOUD(NPD_PROFILE, NPD_LAYER)                               &
!             Total cloud fraction
        , FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             Fractions of different types of cloud
!
!     Fluxes for weighting cloudy diagnostics
      REAL                                                              &     
          FLUX_DIRECT_CLEAR_BAND(NPD_PROFILE, 0: NPD_LAYER)             &
!             Direct clear-sky flux
        , FLUX_DIRECT_BAND(NPD_PROFILE, 0: NPD_LAYER)                   &
!             Direct flux     
        , FLUX_TOTAL_CLEAR_BAND(NPD_PROFILE, 2*NPD_LAYER+2)             &
!             Total (differential) clear-sky flux
        , FLUX_TOTAL_BAND(NPD_PROFILE, 2*NPD_LAYER+2)
!             Total (differential) clear-sky flux
!
!
!
!     Surface fluxes for coupling or diagnostic use
      REAL                                                              &
                !, INTENT(INOUT)
          SEA_FLUX(NPD_PROFILE)                                         &
!           Net downward flux into sea
        , SURFACE_DOWN_FLUX(NPD_PROFILE)                                &
!           Downward flux at surface
        , SURF_DOWN_CLR(NPD_PROFILE)                                    &
!           Clear-sky downward flux at surface
        , SURF_UP_CLR(NPD_PROFILE)                                      &
!           Clear-sky upward flux at surface
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
!     Incremented cloudy fields:
      REAL                                                              &
                !, INTENT(INOUT)
          CLOUD_EXTINCTION(NPD_PROFILE, NPD_LAYER)                      &
!             Overall extinction of clouds weighted with cloud amount
!             and clear-sky direct flux.
        , CLOUD_WEIGHT_EXTINCTION(NPD_PROFILE, NPD_LAYER)               &
!             Weighting for cloud extinction
        , CLOUD_ABSORPTIVITY(NPD_PROFILE, NPD_LAYER)                    &
!             Overall absorptivity of clouds weighted with cloud amount
!             and clear-sky flux.
        , CLOUD_WEIGHT_ABSORPTIVITY(NPD_PROFILE, NPD_LAYER)             &
!             Weighting for cloud absorptivity
        , LS_CLOUD_EXTINCTION(NPD_PROFILE, NPD_LAYER)                   &
!             Extinction of LS clouds weighted with cloud amount
!             and clear-sky direct flux.
        , LS_CLOUD_WEIGHT_EXTINCTION(NPD_PROFILE, NPD_LAYER)            &
!             Weighting for LS cloud extinction
        , LS_CLOUD_ABSORPTIVITY(NPD_PROFILE, NPD_LAYER)                 &
!             Absorptivity of LS clouds weighted with cloud amount
!             and clear-sky flux.
        , LS_CLOUD_WEIGHT_ABSORPTIVITY(NPD_PROFILE, NPD_LAYER)          &
!             Weighting for LS cloud absorptivity
        , CNV_CLOUD_EXTINCTION(NPD_PROFILE, NPD_LAYER)                  &
!             Extinction of conv clouds weighted with cloud amount
!             and clear-sky direct flux.
        , CNV_CLOUD_WEIGHT_EXTINCTION(NPD_PROFILE, NPD_LAYER)           &
!             Weighting for conv cloud extinction
        , CNV_CLOUD_ABSORPTIVITY(NPD_PROFILE, NPD_LAYER)                &
!             Absorptivity of conv clouds weighted with cloud amount
!             and clear-sky flux.
        , CNV_CLOUD_WEIGHT_ABSORPTIVITY(NPD_PROFILE, NPD_LAYER)
!             Weighting for conv cloud absorptivity
!
!
!     Local variables
      INTEGER                                                           &
          L                                                             &
!           Loop variable
        , I                                                             &
!           Loop variable
        , LL
!           Gathered loop variable
!
!
!
!     This is the flux into the sea over the ice-free parts of the
!     grid-box. The model is that the downward fluxes are uniform
!     across the grid-box, but the upward fluxes are not. At this
!     stage no weighting by the actual ice-free area is carried out:
!     that will be done later.
      IF (ISOLIR == IP_SOLAR) THEN
        DO L=1, N_PROFILE
          IF (FLANDG(L) <  1.0.AND.ICE_FRACTION(L) <  1.0) THEN
            SEA_FLUX(L)=SEA_FLUX(L)                                     &
              +FLUX_DOWN(L)*(1.0E+00-ALBEDO_SEA_DIFF(L))                &
              +FLUX_DIRECT(L)*(ALBEDO_SEA_DIFF(L)-ALBEDO_SEA_DIR(L))
          ENDIF
        ENDDO
      ELSE IF (ISOLIR == IP_INFRA_RED) THEN
        IF(L_CTILE)THEN
        DO L=1, N_PROFILE
          IF (FLANDG(L) <  1.0.AND.ICE_FRACTION(L) <  1.0) THEN
            SEA_FLUX(L)=SEA_FLUX(L)                                     &
              +(1.0E+00-ALBEDO_SEA_DIFF(L))                             &
              *(FLUX_DOWN(L)+PLANCK_AIR_SURFACE(L)                      &
              -PLANCK_LEADS_SEA(L))
          ENDIF
        ENDDO
        ELSE
        DO L=1, N_PROFILE
          IF (FLANDG(L) <  1.0) THEN
            SEA_FLUX(L)=SEA_FLUX(L)                                     &
              +(1.0E+00-ALBEDO_SEA_DIFF(L))                             &
              *(FLUX_DOWN(L)+PLANCK_AIR_SURFACE(L)                      &
              -PLANCK_FREEZE_SEA)
          ENDIF
        ENDDO
        ENDIF
      ENDIF
!
      IF (L_SURFACE_DOWN_FLUX) THEN
        IF (ISOLIR == IP_SOLAR) THEN
          DO L=1, N_PROFILE
            SURFACE_DOWN_FLUX(L)=SURFACE_DOWN_FLUX(L)                   &
              +FLUX_DOWN(L)
          ENDDO
        ELSE IF (ISOLIR == IP_INFRA_RED) THEN
          DO L=1, N_PROFILE
            SURFACE_DOWN_FLUX(L)=SURFACE_DOWN_FLUX(L)                   &
              +FLUX_DOWN(L)+PLANCK_AIR_SURFACE(L)
          ENDDO
        ENDIF
      ENDIF
      
      
      IF (ISOLIR == IP_SOLAR) THEN 
      
!
! CALCULATE THE DIFFUSE DOWNWARD FLUX
!

        IF (L_FLUX_DIFFUSE) THEN        
          DO I=0, N_LAYER
            DO L=1, N_PROFILE
               FLUX_DIFFUSE(L,I) = FLUX_DIFFUSE(L,I)                    &
                  + (FLUX_TOTAL_BAND(L,2*I+2)-FLUX_DIRECT_BAND(L,I))
            ENDDO 
          ENDDO
        ENDIF
     
!
! CALCULATE THE UV_FLUXES
!
      
        IF (L_UVFLUX_DIRECT) THEN
          DO I=0, N_LAYER
             DO L=1, N_PROFILE
                UV_FLUX_DIRECT(L,I)=UV_FLUX_DIRECT(L,I)                 &
                      + WEIGHT_UV*FLUX_DIRECT_BAND(L,I)
             ENDDO
          ENDDO
        ENDIF
        IF (L_UVFLUX_UP) THEN
          DO I=0, N_LAYER
            DO L=1, N_PROFILE
               UV_FLUX_UP(L,I)=UV_FLUX_UP(L,I)                         &
                      + WEIGHT_UV*FLUX_TOTAL_BAND(L,2*I+1)          
            ENDDO
          ENDDO  
        ENDIF
        IF (L_UVFLUX_DOWN) THEN
          DO I=0, N_LAYER
            DO L=1, N_PROFILE
                UV_FLUX_DOWN(L,I)=UV_FLUX_DOWN(L,I)                     &
                      + WEIGHT_UV*FLUX_TOTAL_BAND(L,2*I+2)          
            ENDDO
          ENDDO           
        ENDIF     
         
      ENDIF
         
      IF (L_SURF_DOWN_CLR) THEN
        IF (ISOLIR == IP_SOLAR) THEN
          DO L=1, N_PROFILE
            SURF_DOWN_CLR(L)=SURF_DOWN_CLR(L)                           &
              +FLUX_DOWN_CLEAR(L)
          ENDDO
        ELSE IF (ISOLIR == IP_INFRA_RED) THEN
          DO L=1, N_PROFILE
            SURF_DOWN_CLR(L)=SURF_DOWN_CLR(L)                           &
              +FLUX_DOWN_CLEAR(L)+PLANCK_AIR_SURFACE(L)
          ENDDO
        ENDIF
      ENDIF
!
      IF (L_SURF_UP_CLR) THEN
        IF (ISOLIR == IP_SOLAR) THEN
          DO L=1, N_PROFILE
            SURF_UP_CLR(L)=SURF_UP_CLR(L)                               &
              +FLUX_UP_CLEAR(L)
          ENDDO
        ELSE IF (ISOLIR == IP_INFRA_RED) THEN
          DO L=1, N_PROFILE
            SURF_UP_CLR(L)=SURF_UP_CLR(L)                               &
              +FLUX_UP_CLEAR(L)+PLANCK_AIR_SURFACE(L)
          ENDDO
        ENDIF
      ENDIF
!
!     This diagnostic is available only in the solar region. Over
!     sea-points it refers only to the flux over the open sea
!     (see the comments about SEA_FLUX above). Over land, both the
!     upward and downward fluxes are taken as uniform.
      IF (L_FLUX_BELOW_690NM_SURF) THEN
        IF (ISOLIR == IP_SOLAR) THEN
         IF (L_CTILE) THEN
          DO L=1, N_PROFILE
            FLUX_BELOW_690NM_SURF(L)=FLUX_BELOW_690NM_SURF(L)           &
                +WEIGHT_690NM*(FLUX_DOWN(L)-FLUX_UP(L))
            IF(FLANDG(L) <  1.0.AND.ICE_FRACTION(L) <  1.0) THEN
              FL_SEA_BELOW_690NM_SURF(L)                                &
                =FL_SEA_BELOW_690NM_SURF(L)                             &
                +WEIGHT_690NM                                           &
                *(FLUX_DOWN(L)*(1.0E+00-ALBEDO_SEA_DIFF(L))             &
                +FLUX_DIRECT(L)*(ALBEDO_SEA_DIFF(L)                     &
                -ALBEDO_SEA_DIR(L)))
            ENDIF
          ENDDO
         ELSE
          DO L=1, N_PROFILE
            IF (FLANDG(L) >  0.0) THEN
              FLUX_BELOW_690NM_SURF(L)=FLUX_BELOW_690NM_SURF(L)         &
                +WEIGHT_690NM*(FLUX_DOWN(L)-FLUX_UP(L))
            ELSE
              FLUX_BELOW_690NM_SURF(L)=FLUX_BELOW_690NM_SURF(L)         &
                +WEIGHT_690NM                                           &
                *(FLUX_DOWN(L)*(1.0E+00-ALBEDO_SEA_DIFF(L))             &
                +FLUX_DIRECT(L)*(ALBEDO_SEA_DIFF(L)-ALBEDO_SEA_DIR(L)))
            ENDIF
          ENDDO
         ENDIF
        ENDIF
      ENDIF
!
!     Surface shortwave diagnostics required for MOSES II
      IF (L_MOSES_II) THEN
        IF (ISOLIR == IP_SOLAR) THEN
          DO L=1, N_PROFILE
            SURF_VIS_DIR(L) = SURF_VIS_DIR(L) +                         &
                              WEIGHT_690NM*FLUX_DIRECT(L)
            SURF_NIR_DIR(L) = SURF_NIR_DIR(L) +                         &
                              (1. - WEIGHT_690NM)*FLUX_DIRECT(L)
            SURF_VIS_DIF(L) = SURF_VIS_DIF(L) +                         &
                              WEIGHT_690NM*(FLUX_DOWN(L)-FLUX_DIRECT(L))
            SURF_NIR_DIF(L) = SURF_NIR_DIF(L) +                         &
                       (1. - WEIGHT_690NM)*(FLUX_DOWN(L)-FLUX_DIRECT(L))
          ENDDO
        ENDIF
      ENDIF

!
!     Spectral diagnostics for clouds:
!
      IF (L_CLOUD_EXTINCTION) THEN
!
!        Increment the arrays of diagnostics. The extinction
!        calculated in this band is a mean value weighted with the
!        fractions of individual types of cloud (which sum to 1).
!        Here it is weighted with the clear-sky direct solar
!        flux in the band at the top of the current
!        layer and the total amount of cloud in the grid-box.
!        This definition has the advantage of convenience, but there
!        appears to be no optimal definition of an average extinction.
!
         DO I=N_CLOUD_TOP, N_LAYER
!CDIR NODEP
            DO LL=1, N_CLOUD_PROFILE(I)
               L=I_CLOUD_PROFILE(LL, I)
               CLOUD_WEIGHT_EXTINCTION(L, I)                            &
                  =CLOUD_WEIGHT_EXTINCTION(L, I)                        &
                  +W_CLOUD(L, I)*FLUX_DIRECT_CLEAR_BAND(L, I-1)
               CLOUD_EXTINCTION(L, I)                                   &
                  =CLOUD_EXTINCTION(L, I)                               &
                  +W_CLOUD(L, I)*FLUX_DIRECT_CLEAR_BAND(L, I-1)         &
                  *CLOUD_EXTINCTION_BAND(L, I)
            ENDDO
         ENDDO
!
      ENDIF
!
      IF (L_CLOUD_ABSORPTIVITY) THEN
!
!        Increment the arrays of diagnostics. The absorptivity
!        calculated in this band is a mean value weighted with the
!        fractions of individual types of cloud (which sum to 1).
!        Here it is weighted with the modulus of the clear_sky
!        differential flux in the band at the top of the current
!        layer and the total amount of cloud in the grid-box, as the
!        diagnostic is a measure of the effect of introducing an
!        infinitesimal layer of layer at the top of the current
!        layer into a clear atmosphere on the upward flux at the
!        top of the cloud.
!
         DO I=N_CLOUD_TOP, N_LAYER
!CDIR NODEP
            DO LL=1, N_CLOUD_PROFILE(I)
               L=I_CLOUD_PROFILE(LL, I)
               CLOUD_WEIGHT_ABSORPTIVITY(L, I)                          &
                  =CLOUD_WEIGHT_ABSORPTIVITY(L, I)                      &
                  +W_CLOUD(L, I)*ABS(FLUX_TOTAL_CLEAR_BAND(L, 2*I-1))
               CLOUD_ABSORPTIVITY(L, I)                                 &
                  =CLOUD_ABSORPTIVITY(L, I)                             &
                  +W_CLOUD(L, I)*ABS(FLUX_TOTAL_CLEAR_BAND(L, 2*I-1))   &
                  *CLOUD_ABSORPTIVITY_BAND(L, I)
            ENDDO
         ENDDO
!
      ENDIF
      IF (L_LS_CLOUD_EXTINCTION) THEN
!
!        Increment the arrays of diagnostics. The extinction
!        calculated in this band is a mean value weighted with the
!        fractions of individual types of cloud (which sum to 1).
!        Here it is weighted with the clear-sky direct solar
!        flux in the band at the top of the current
!        layer and the total amount of cloud in the grid-box.
!        This definition has the advantage of convenience, but there
!        appears to be no optimal definition of an average extinction.
!
         DO I=N_CLOUD_TOP, N_LAYER
!CDIR NODEP
            DO LL=1, N_CLOUD_PROFILE(I)
               L=I_CLOUD_PROFILE(LL, I)
               LS_CLOUD_WEIGHT_EXTINCTION(L, I)                         &
                  =LS_CLOUD_WEIGHT_EXTINCTION(L, I)                     &
                  +W_CLOUD(L, I)*(FRAC_CLOUD(L,I,1)+FRAC_CLOUD(L,I,2))  &
                   *FLUX_DIRECT_CLEAR_BAND(L, I-1)
               LS_CLOUD_EXTINCTION(L, I)                                &
                  =LS_CLOUD_EXTINCTION(L, I)                            &
                  +W_CLOUD(L, I)*(FRAC_CLOUD(L,I,1)+FRAC_CLOUD(L,I,2))  &
                  *FLUX_DIRECT_CLEAR_BAND(L, I-1)                       &
                  *LS_CLOUD_EXTINCTION_BAND(L, I)
            ENDDO
         ENDDO
!
      ENDIF
!
      IF (L_LS_CLOUD_ABSORPTIVITY) THEN
!
!        Increment the arrays of diagnostics. The absorptivity
!        calculated in this band is a mean value weighted with the
!        fractions of individual types of cloud (which sum to 1).
!        Here it is weighted with the modulus of the clear_sky
!        differential flux in the band at the top of the current
!        layer and the total amount of cloud in the grid-box as the
!        diagnostic is a measure of the effect of introducing an
!        infinitesimal layer of layer at the top of the current
!        layer into a clear atmosphere on the upward flux at the
!        top of the cloud.
!
         DO I=N_CLOUD_TOP, N_LAYER
!CDIR NODEP
            DO LL=1, N_CLOUD_PROFILE(I)
               L=I_CLOUD_PROFILE(LL, I)
               LS_CLOUD_WEIGHT_ABSORPTIVITY(L, I)                       &
                  =LS_CLOUD_WEIGHT_ABSORPTIVITY(L, I)                   &
                  +W_CLOUD(L, I)*(FRAC_CLOUD(L,I,1)+FRAC_CLOUD(L,I,2))  &
                   *ABS(FLUX_TOTAL_CLEAR_BAND(L, 2*I-1))
               LS_CLOUD_ABSORPTIVITY(L, I)                              &
                  =LS_CLOUD_ABSORPTIVITY(L, I)                          &
                  +W_CLOUD(L, I)*(FRAC_CLOUD(L,I,1)+FRAC_CLOUD(L,I,2))  &
                  *ABS(FLUX_TOTAL_CLEAR_BAND(L, 2*I-1))                 &
                  *LS_CLOUD_ABSORPTIVITY_BAND(L, I)
            ENDDO
         ENDDO
!
      ENDIF
      IF (L_CNV_CLOUD_EXTINCTION) THEN
!
!        Increment the arrays of diagnostics. The extinction
!        calculated in this band is a mean value weighted with the
!        fractions of individual types of cloud (which sum to 1).
!        Here it is weighted with the clear-sky direct solar
!        flux in the band at the top of the current
!        layer and the total amount of cloud in the grid-box.
!        This definition has the advantage of convenience, but there
!        appears to be no optimal definition of an average extinction.
!
         DO I=N_CLOUD_TOP, N_LAYER
!CDIR NODEP
            DO LL=1, N_CLOUD_PROFILE(I)
               L=I_CLOUD_PROFILE(LL, I)
               CNV_CLOUD_WEIGHT_EXTINCTION(L, I)                        &
                  =CNV_CLOUD_WEIGHT_EXTINCTION(L, I)                    &
                  +W_CLOUD(L, I)*(FRAC_CLOUD(L,I,3)+FRAC_CLOUD(L,I,4))  &
                   *FLUX_DIRECT_CLEAR_BAND(L, I-1)
               CNV_CLOUD_EXTINCTION(L, I)                               &
                  =CNV_CLOUD_EXTINCTION(L, I)                           &
                  +W_CLOUD(L, I)*(FRAC_CLOUD(L,I,3)+FRAC_CLOUD(L,I,4))  &
                  *FLUX_DIRECT_CLEAR_BAND(L, I-1)                       &
                  *CNV_CLOUD_EXTINCTION_BAND(L, I)
            ENDDO
         ENDDO
!
      ENDIF
!
      IF (L_CNV_CLOUD_ABSORPTIVITY) THEN
!
!        Increment the arrays of diagnostics. The absorptivity
!        calculated in this band is a mean value weighted with the
!        fractions of individual types of cloud (which sum to 1).
!        Here it is weighted with the modulus of the clear_sky
!        differential flux in the band at the top of the current
!        layer and the total amount of cloud in the grid-box as the
!        diagnostic is a measure of the effect of introducing an
!        infinitesimal layer of layer at the top of the current
!        layer into a clear atmosphere on the upward flux at the
!        top of the cloud.
!
         DO I=N_CLOUD_TOP, N_LAYER
!CDIR NODEP
            DO LL=1, N_CLOUD_PROFILE(I)
               L=I_CLOUD_PROFILE(LL, I)
               CNV_CLOUD_WEIGHT_ABSORPTIVITY(L, I)                      &
                  =CNV_CLOUD_WEIGHT_ABSORPTIVITY(L, I)                  &
                  +W_CLOUD(L, I)*(FRAC_CLOUD(L,I,3)+FRAC_CLOUD(L,I,4))  &
                   *ABS(FLUX_TOTAL_CLEAR_BAND(L, 2*I-1))
               CNV_CLOUD_ABSORPTIVITY(L, I)                             &
                  =CNV_CLOUD_ABSORPTIVITY(L, I)                         &
                  +W_CLOUD(L, I)*(FRAC_CLOUD(L,I,3)+FRAC_CLOUD(L,I,4))  &
                  *ABS(FLUX_TOTAL_CLEAR_BAND(L, 2*I-1))                 &
                  *CNV_CLOUD_ABSORPTIVITY_BAND(L, I)
            ENDDO
         ENDDO
!
      ENDIF
!
      RETURN
      END SUBROUTINE R2_COUPLE_DIAG
#endif
