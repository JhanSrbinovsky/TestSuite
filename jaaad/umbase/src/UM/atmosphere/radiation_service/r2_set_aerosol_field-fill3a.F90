#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set fields of aerosols.
!
! Purpose:
!   The mixing ratios of aerosols are transferred to the large array.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: James Manners
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SET_AEROSOL_FIELD(IERR                              &
     &   , N_PROFILE, NLEVS, N_LAYER, N_AEROSOL, TYPE_AEROSOL           &
     &   , I_GATHER, L_EXTRA_TOP                                        &
     &   , L_CLIMAT_AEROSOL, L_CLIM_AERO_HGT, L_HadGEM1_Clim_Aero       &
     &   , BL_DEPTH, T, N_LEVELS_BL, L_MURK_RAD, AERO_MESO              &
     &   , L_USE_DUST, DUST_DIM1, DUST_DIM2                             &
     &   , DUST_1, DUST_2, DUST_3, DUST_4, DUST_5, DUST_6               &
     &   , L_USE_BIOGENIC, BIOGENIC_DIM1, BIOGENIC_DIM2, BIOGENIC       &
     &   , L_USE_SULPC_DIRECT                                           &
     &   , SULP_DIM1, SULP_DIM2                                         &
     &   , ACCUM_SULPHATE, AITKEN_SULPHATE                              &
     &   , L_VOLCTS, VOLCMASS                                           &
     &   , L_USE_SEASALT_DIRECT, SALT_DIM_A, SALT_DIM_B                 &
     &   , SEA_SALT_FILM, SEA_SALT_JET, P                               &
     &,L_USE_SOOT_DIRECT, SOOT_DIM1, SOOT_DIM2, FRESH_SOOT, AGED_SOOT   &
     &   , L_USE_BMASS_DIRECT, BMASS_DIM1, BMASS_DIM2                   &
     &   , FRESH_BMASS, AGED_BMASS                                      &
     &   , L_USE_OCFF_DIRECT, OCFF_DIM1, OCFF_DIM2                      &
     &   , FRESH_OCFF, AGED_OCFF                                        &
     &   , N_ARCL_SPECIES, N_ARCL_COMPNTS, I_ARCL_COMPNTS               &
     &   , L_USE_ARCL, ARCL_DIM1, ARCL_DIM2, ARCL                       &
     &   , LAND, LYING_SNOW, PSTAR                                      &
     &   , P_LAYER_BOUNDARIES                                           &
     &   , TRINDX                                                       &
     &   , AEROSOL_MIX_RATIO                                            &
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_AEROSOL_SPECIES       &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED.
#include "c_g.h"
#include "stdio3a.h"
#include "error3a.h"
#include "aercmp3a.h"
#include "arcl_ids.h"
#include "c_r_cp.h"
#include "c_pi.h"
!
!     DUMMY ARGUMENTS.
!
!     SIZES OF ARRAYS:
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_FIELD                                                    &
!             FIELD SIZE IN CALLING PROGRAM
     &   , NPD_PROFILE                                                  &
!             SIZE OF ARRAY OF PROFILES
     &   , NPD_LAYER                                                    &
!             MAXIMUM NUMBER OF LAYERS
     &   , NPD_AEROSOL_SPECIES
!             MAXIMUM NUMBER OF AEROSOL SPECIES
!
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
!
!     ACTUAL SIZES USED:
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , NLEVS                                                        &
!             Number of layers used outside the radiation scheme
     &   , N_LAYER                                                      &
!             Number of layers seen by radiation
     &   , N_LEVELS_BL                                                  &
!             Number of layers occupied by boundary-layer aerosol
!             if L_CLIM_AERO_HGT is false.
     &   , N_AEROSOL                                                    &
!             NUMBER OF AEROSOLS IN SPECTRAL FILE
     &   , TYPE_AEROSOL(NPD_AEROSOL_SPECIES)
!             ACTUAL TYPES OF AEROSOLS
!
!     GATHERING ARRAY:
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_GATHER(NPD_FIELD)
!             LIST OF POINTS TO GATHER
!
!     FLAG FOR THE CLIMATOLOGICAL AEROSOL DISTRIBUTION.
      LOGICAL                                                           &
                   !, INTENT(IN)
     &     L_CLIMAT_AEROSOL                                             &
!             FLAG FOR CLIMATOLOGICAL AEROSOL DISTRIBUTION
     &   , L_CLIM_AERO_HGT                                              &
!             Flag to use the boundary layer depth in setting the
!             climatological aerosol
     &   , L_HadGEM1_Clim_Aero                                          &
!             Flag to use HadGEM1 setting for climatological aerosols
     &   , L_MURK_RAD                                                   &
!             FLAG FOR MESOSCALE MODEL AEROSOL
     &   , L_EXTRA_TOP
!             Flag to include an extra top layer in the radiation
!             scheme
!
! Declare mineral dust aerosol variables:
      LOGICAL L_USE_DUST !Use direct radiative effect of mineral dust
      INTEGER DUST_DIM1,DUST_DIM2
!        Dimensions for mineral dust arrays (P_FIELD,P_LEVELS or 1,1)
      REAL DUST_1(DUST_DIM1, DUST_DIM2)                                 &
                                        !Mass mixing ratio of div1 dust
     &   , DUST_2(DUST_DIM1, DUST_DIM2)                                 &
                                        !Mass mixing ratio of div2 dust
     &   , DUST_3(DUST_DIM1, DUST_DIM2)                                 &
                                        !Mass mixing ratio of div3 dust
     &   , DUST_4(DUST_DIM1, DUST_DIM2)                                 &
                                        !Mass mixing ratio of div4 dust
     &   , DUST_5(DUST_DIM1, DUST_DIM2)                                 &
                                        !Mass mixing ratio of div5 dust
     &   , DUST_6(DUST_DIM1, DUST_DIM2) !Mass mixing ratio of div6 dust
!     VARIABLES FOR THE BIOGENIC AEROSOL
      LOGICAL L_USE_BIOGENIC
      INTEGER BIOGENIC_DIM1, BIOGENIC_DIM2
      REAL BIOGENIC(BIOGENIC_DIM1, BIOGENIC_DIM2)

!     VARIABLES FOR THE SULPHUR CYCLE:
      LOGICAL                                                           &
                   !, INTENT(IN)
     &     L_USE_SULPC_DIRECT
!             FLAG TO USE SULPHUR CYCLE FOR DIRECT EFFECT
      INTEGER                                                           &
                   !, INTENT(IN)
     &     SULP_DIM1,SULP_DIM2
!             DIMENSIONS FOR _SULPHATE ARRAYS, (P_FIELD,P_LEVELS or 1,1)
      REAL                                                              &
                !, INTENT(IN)
     &     ACCUM_SULPHATE(SULP_DIM1, SULP_DIM2)                         &
!             MASS MIXING RATIOS OF ACCUMULATION MODE AEROSOL
     &   , AITKEN_SULPHATE(SULP_DIM1, SULP_DIM2)
!             MASS MIXING RATIOS OF AITKEN MODE AEROSOL
!
! Declare soot variables:
      LOGICAL L_USE_SOOT_DIRECT !USE DIRECT RAD. EFFECT OF SOOT AEROSOL
      INTEGER SOOT_DIM1,SOOT_DIM2
                !DIMENSIONS FOR SOOT ARRAYS, (P_FIELD,P_LEVELS or 1,1)
      REAL FRESH_SOOT(SOOT_DIM1, SOOT_DIM2)                             &
                                                 ! MMR OF FRESH SOOT
     &   , AGED_SOOT(SOOT_DIM1, SOOT_DIM2)       ! MMR OF AGED SOOT
!
! Declare biomass smoke aerosol variables:
      LOGICAL L_USE_BMASS_DIRECT
!              Use direct radiative effect of biomass smoke
      INTEGER BMASS_DIM1,BMASS_DIM2
!              Dimensions for biomass arrays (P_FIELD,P_LEVELS or 1,1)
      REAL FRESH_BMASS(BMASS_DIM1, BMASS_DIM2)                          &
                                                    ! MMR OF FRESH SMOKE
     &   , AGED_BMASS(BMASS_DIM1, BMASS_DIM2)       ! MMR OF AGED SMOKE
!
! Declare fossil-fuel organic carbon aerosol variables:
      LOGICAL L_USE_OCFF_DIRECT
!              Use direct radiative effect of OCFF aerosol
      INTEGER OCFF_DIM1, OCFF_DIM2
!              Dimensions for OCFF arrays (P_FIELD,P_LEVELS or 1,1)
      REAL FRESH_OCFF(OCFF_DIM1, OCFF_DIM2)                             &
                                                    ! MMR OF FRESH OCFF
     &   , AGED_OCFF(OCFF_DIM1, OCFF_DIM2)          ! MMR OF AGED OCFF
!
      LOGICAL                                                           &
     &     L_USE_SEASALT_DIRECT
!             Flag for direct effect of interactive sea-salt aerosol
!
      INTEGER                                                           &
     &     SALT_DIM_A, SALT_DIM_B
!             Array sizes of sea-salt aerosols
!                (either P_FIELD,P_LEVELS or 1,1)
!
      REAL                                                              &
     &     SEA_SALT_FILM(SALT_DIM_A, SALT_DIM_B)                        &
!             On input, number concentration (m-3) of film-mode
!             sea-salt aerosols; converted to mass mixing ratio.
     &   , SEA_SALT_JET(SALT_DIM_A, SALT_DIM_B)
!             On input, number concentration (m-3) of jet-mode
!             sea-salt aerosols; converted to mass mixing ratio.
!
! Aerosol climatology for NWP:
!
#include "arcl_dim.h"

      ! Number of requested species within the climatology
      Integer N_ARCL_SPECIES
      
      ! Corresponding number of requested components
      Integer N_ARCL_COMPNTS
      
      ! Model switches for each species
      Logical, dimension(NPD_ARCL_SPECIES) :: L_USE_ARCL
      
      ! Array index of components
      Integer, dimension(NPD_ARCL_COMPNTS) :: I_ARCL_COMPNTS
      
      ! Array dimensions
      Integer ARCL_DIM1, ARCL_DIM2
      
      ! Mass-mixing ratios
      Real                                                              &
     &    ARCL(ARCL_DIM1, ARCL_DIM2, N_ARCL_COMPNTS)
      
      REAL                                                              &
     &     AERO_MESO(NPD_FIELD, NLEVS)
!             MIXING RATIO OF 'URBAN' AEROSOL OF MESOSCALE MODEL
!     GENERAL ATMOSPHERIC PROPERTIES:
      INTEGER                                                           &
                !, INTENT(IN)
     &     TRINDX(NPD_FIELD)
!             LAYER BOUNDARY OF TROPOPAUSE
      REAL                                                              &
                !, INTENT(IN)
     &     PSTAR(NPD_FIELD)                                             &
!             SURFACE PRESSURES
     &   , VOLCMASS(NPD_FIELD)						&
!             Mass of stratospheric volcanic aerosol at each point 
     &,    P_LAYER_BOUNDARIES(NPD_FIELD,0:NLEVS)
!             PRESSURE AT BOUNDARIES OF LAYERS
!
      LOGICAL  L_VOLCTS

!     SURFACE FIELDS
      LOGICAL                                                           &
                !, INTENT(IN)
     &     LAND(NPD_FIELD)
!             LAND SEA MASK
      REAL                                                              &
                !, INTENT(IN)
     &     LYING_SNOW(NPD_FIELD)                                        &
!             DEPTH OF LYING SNOW
     &   , BL_DEPTH(NPD_FIELD)                                          &
!             Depth of the boundary layer
     &   , T(NPD_PROFILE, 0:NPD_LAYER)                                  &
!             Temperatures of atmospheric layers
     &   , P(NPD_PROFILE, 0:NPD_LAYER)
!             Pressures at layer centres
!
      REAL                                                              &
                !, INTENT(OUT)
     &     AEROSOL_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER                  &
     &        , NPD_AEROSOL_SPECIES)
!             MIXING RATIOS OF AEROSOLS
!
!
!
!     LOCAL VARIABLES:
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , J                                                            &
!             LOOP VARIABLE
     &   , L                                                            &
!             LOOP VARIABLE
     &   , LG                                                           &
!             INDEX FOR GATHERING
     &   , BLTOP                                                        &
!             INDEX OF UPPER BOUNDARY OF PLANETARY BOUNDARY LAYER
     &   , I_AEROSOL                                                    &
!             ACTUAL TYPE OF AEROSOL BEING CONSIDERED
     &   , I_TOP_COPY
!             Topmost radiative layer to be set directly from the
!             profile supplied
!
!
!     ARRAYS FOR THE CLIMATOLOGICAL AEROSOL MODEL
      Integer :: ii
!       Variable used in initialization of arrays.
      Logical, Dimension(NPD_AEROSOL_COMPONENT) :: L_IN_CLIMAT =        &
     &  (/ (.TRUE., ii=1, 4), .FALSE., .TRUE.,                          &
     &     (.FALSE., ii=7, NPD_AEROSOL_COMPONENT) /)
!       Flags to indicate which aerosols are included in the
!       climatology: this may be used to enable various components
!       to be replaced by fully prognostic schemes.
      Integer, Dimension(NPD_AEROSOL_COMPONENT) :: I_CLIM_POINTER =     &
     &    (/ 1, 2, 3, 4, 0, 5,                                          &
     &    (0, ii=7, NPD_AEROSOL_COMPONENT) /)

!       Pointers to the indices of the original climatological
!       aerosol model.
      REAL                                                              &
     &     AEROSOL_MIX_RATIO_CLIM(NPD_PROFILE, 0: NPD_LAYER, 5)
!             MIXING RATIOS OF THE CLIMATOLOGICAL AEROSOLS
!
      REAL                                                              &
     &     MODE_RADIUS_SS_FILM                                          &
!            Mode radius of film-mode sea-salt aerosol (m)
     &   , MODE_RADIUS_SS_JET                                           &
!            Mode radius of jet-mode sea-salt aerosol (m)
     &   , SIGMA_SS_FILM                                                &
!            Geometric standard deviation of film-mode sea-salt aerosol
     &   , SIGMA_SS_JET                                                 &
!            Geometric standard deviation of jet-mode sea-salt aerosol
     &   , DENSITY_SEA_SALT
!            Bulk density of sodium chloride, taken as representative
!            of sea-salt aerosol (kg m-3); taken from Pruppacher & Klett
!            2nd edition (1997), p. 244.
!
      Integer im
      Real, Dimension(NPD_AEROSOL_COMPONENT) :: Meso_frac =             &
     &  (/ 0.61, 0.17, 0.0, 0.22, (0.0, im=5,NPD_AEROSOL_COMPONENT) /)
      REAL INV17

      PARAMETER(                                                        &
     &     MODE_RADIUS_SS_FILM=0.1E-06                                  &
     &   , MODE_RADIUS_SS_JET=1.0E-06                                   &
     &   , SIGMA_SS_FILM=1.9                                            &
     &   , SIGMA_SS_JET=2.0                                             &
     &   , DENSITY_SEA_SALT=2165.0                                      &
     &   )
!
!
!
!     SUBROUTINES CALLED:
      EXTERNAL                                                          &
     &     R2_SET_AERO_CLIM_HADCM3
!
!
!
!     If an extra layer is used in the radiation scheme, any
!     aerosol profiles supplied will be copied into the profiles
!     seen by the radiation scheme starting with the second layer
!     down, rather than the first.
      I_TOP_COPY=1
      IF (L_EXTRA_TOP) I_TOP_COPY=2
!
!  Use HadGEM1 settings to remove boundary layer climatological aerosols.
      IF (L_HadGEM1_Clim_Aero) THEN
        DO ii=1,N_AEROSOL
          IF ((TYPE_AEROSOL(ii) == IP_WATER_SOLUBLE).OR.               &
     &        (TYPE_AEROSOL(ii) == IP_DUST_LIKE) .OR.                  &
     &        (TYPE_AEROSOL(ii) == IP_OCEANIC).OR.                     &
     &        (TYPE_AEROSOL(ii) == IP_SOOT)) THEN
            L_IN_CLIMAT(ii) = .FALSE.
          END IF
        END DO
      END IF
!
!  Use climatological soot if climatological aerosols are on and not
!  using interactive soot.
       L_IN_CLIMAT(IP_SOOT) = L_IN_CLIMAT(IP_SOOT)                      &
     &                     .AND.(.NOT.L_USE_SOOT_DIRECT)
!
!  Use climatological "oceanic" aerosol if climatological aerosols have
!  been selected and not using sea-salt parametrization.
      L_IN_CLIMAT(IP_OCEANIC) = L_IN_CLIMAT(IP_OCEANIC)                 &
     &                          .AND.(.NOT. L_USE_SEASALT_DIRECT)
!
!  The climatological water-soluble aerosol should not be used
!  if the direct effects of sulphate aerosols are included.
!  (Note that this differs the situation applying in earlier
!  versions of the model, such as HadAM3).
      L_IN_CLIMAT(IP_WATER_SOLUBLE) = L_IN_CLIMAT(IP_WATER_SOLUBLE)     &
     &                     .AND.(.NOT.L_USE_SULPC_DIRECT)
!
!
!
      IF (L_CLIMAT_AEROSOL) THEN
!
!        SET THE MIXING RATIOS OF THE CLIMATOLOGICAL AEROSOLS
!        USED IN THE CLIMATOLOGY OF HADCM3. A SEPARATE SUBROUTINE
!        IS USED TO ENSURE BIT-REPRODUCIBLE RESULTS BY USING
!        EARLIER CODE. THIS COULD BE ALTERED IF A NEW CLIMATOLOGY WERE
!        USED.
!
! DEPENDS ON: r2_set_aero_clim_hadcm3
         CALL R2_SET_AERO_CLIM_HADCM3(N_PROFILE, NLEVS, N_LAYER         &
     &      , I_GATHER, L_EXTRA_TOP                                     &
     &      , L_CLIM_AERO_HGT, BL_DEPTH, T, N_LEVELS_BL                 &
     &      , LAND, LYING_SNOW, PSTAR, P_LAYER_BOUNDARIES, TRINDX       &
     &      , L_VOLCTS, VOLCMASS                                        &
     &      , AEROSOL_MIX_RATIO_CLIM                                    &
     &      , NPD_FIELD, NPD_PROFILE, NPD_LAYER                         &
     &      )
!
      ENDIF
!
!
!     THE AEROSOLS REQUIRED BY FOR THE CALCULATION SHOULD HAVE BEEN
!     SELECTED WHEN THE SPECTRAL FILE WAS READ IN. EACH TYPE SHOULD
!     BE SET APPROPRIATELY.
!
      DO J=1, N_AEROSOL
!
         I_AEROSOL=TYPE_AEROSOL(J)
!
         IF (L_CLIMAT_AEROSOL.AND.L_IN_CLIMAT(I_AEROSOL)) THEN
!
!           Here mxing ratios in all layers are set because the
!           possible presence of an extra layer has been allowed
!           for in the subroutine.
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  AEROSOL_MIX_RATIO(L, I, J)                            &
     &               =AEROSOL_MIX_RATIO_CLIM(L, I                       &
     &               , I_CLIM_POINTER(I_AEROSOL))
               ENDDO
            ENDDO
!
! For the other aerosol species, they may be supplied from dedicated
! schemes (e.g. climate model) or climatologies (e.g. NWP model).
! The former are indicated by L_USE_<SPECIES>, the latter by
! L_USE_ARCL(IP_ARCL_<SPECIES>). Should both logical be set to true
! for a given aerosol species (i.e. an aerosol is provided by a 
! dedicated scheme AND a climatology), the climatology wins and will 
! be seen by radiation.
!
!
!
!     Mineral dust aerosols:
!
         ELSE IF (I_AEROSOL == IP_DUST_1) THEN

           IF (L_USE_ARCL(IP_ARCL_DUST)) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J)                            &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_DUST_B1))
                END DO
              END DO           
           ELSE IF (L_USE_DUST) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J)                            &
     &               =DUST_1(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF
!
         ELSE IF (I_AEROSOL == IP_DUST_2) THEN

           IF (L_USE_ARCL(IP_ARCL_DUST)) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J)                            &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_DUST_B2))
                END DO
              END DO           
           ELSE IF (L_USE_DUST) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J)                            &
     &               =DUST_2(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF
!
         ELSE IF (I_AEROSOL == IP_DUST_3) THEN

           IF (L_USE_ARCL(IP_ARCL_DUST)) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J)                            &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_DUST_B3))
                END DO
              END DO           
           ELSE IF (L_USE_DUST) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J)                            &
     &               =DUST_3(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF
!
         ELSE IF (I_AEROSOL == IP_DUST_4) THEN

           IF (L_USE_ARCL(IP_ARCL_DUST)) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J)                            &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_DUST_B4))
                END DO
              END DO           
           ELSE IF (L_USE_DUST) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J)                            &
     &               =DUST_4(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF
!
         ELSE IF (I_AEROSOL == IP_DUST_5) THEN

           IF (L_USE_ARCL(IP_ARCL_DUST)) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J)                            &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_DUST_B5))
                END DO
              END DO           
           ELSE IF (L_USE_DUST) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J)                            &
     &               =DUST_5(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF
!
         ELSE IF (I_AEROSOL == IP_DUST_6) THEN

           IF (L_USE_ARCL(IP_ARCL_DUST)) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J)                            &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_DUST_B6))
                END DO
              END DO           
           ELSE IF (L_USE_DUST) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J)                            &
     &               =DUST_6(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF
!
         ELSE IF (I_AEROSOL == IP_ACCUM_SULPHATE) THEN
!
!           Aerosols related to the sulphur cycle (note that dissolved
!           sulphate does not contribute to the direct effect):
!
           IF (L_USE_ARCL(IP_ARCL_SULP)) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J)                           &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_SULP_AC))
                END DO
              END DO
           ELSE IF (L_USE_SULPC_DIRECT) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J)                           &            
     &               =ACCUM_SULPHATE(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF
!
         ELSE IF (I_AEROSOL == IP_AITKEN_SULPHATE) THEN
           
           IF (L_USE_ARCL(IP_ARCL_SULP)) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J)                           &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_SULP_AK))
                END DO
              END DO
           ELSE IF (L_USE_SULPC_DIRECT) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J)                           &            
     &               =AITKEN_SULPHATE(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF
!
         ELSE IF (I_AEROSOL == IP_FRESH_SOOT) THEN

           IF (L_USE_ARCL(IP_ARCL_BLCK)) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J)                           &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_BLCK_FR))
                END DO
              END DO
           ELSE IF (L_USE_SOOT_DIRECT) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J)                           &            
     &               =FRESH_SOOT(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF
!
         ELSE IF (I_AEROSOL == IP_AGED_SOOT) THEN

           IF (L_USE_ARCL(IP_ARCL_BLCK)) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J)                           &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_BLCK_AG))
                END DO
              END DO
           ELSE IF (L_USE_SOOT_DIRECT) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J)                           &            
     &               =AGED_SOOT(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF
!
         ELSE IF (I_AEROSOL == IP_BIOMASS_1) THEN

           IF (L_USE_ARCL(IP_ARCL_BIOM)) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J)                           &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_BIOM_FR))
                END DO
              END DO
           ELSE IF (L_USE_BMASS_DIRECT) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J)                           &            
     &               =FRESH_BMASS(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF
!
         ELSE IF (I_AEROSOL == IP_BIOMASS_2) THEN

           IF (L_USE_ARCL(IP_ARCL_BIOM)) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J)                           &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_BIOM_AG))
                END DO
              END DO
           ELSE IF (L_USE_BMASS_DIRECT) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J)                           &            
     &               =AGED_BMASS(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF
!
         ELSE IF (I_AEROSOL == IP_SEASALT_FILM) THEN
     
           IF (L_USE_ARCL(IP_ARCL_SSLT)) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J)=                          &
     &                    (ARCL(LG, N_LAYER+1-I,                        &
     &                     I_ARCL_COMPNTS(IP_ARCL_SSLT_FI))*4.0*PI      &
     &                   *DENSITY_SEA_SALT*(MODE_RADIUS_SS_FILM**3.0)   &
     &                   *EXP(4.5*(ALOG(SIGMA_SS_FILM))**2.0))          &
     &                   /(3.0*P(L, I)/(R*T(L, I)))
                END DO
              END DO
           ELSE IF (L_USE_SEASALT_DIRECT) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J)=                          &
     &                    (SEA_SALT_FILM(LG, N_LAYER+1-I)*4.0*PI        &
     &                   *DENSITY_SEA_SALT*(MODE_RADIUS_SS_FILM**3.0)   &
     &                   *EXP(4.5*(ALOG(SIGMA_SS_FILM))**2.0))          &
     &                   /(3.0*P(L, I)/(R*T(L, I)))
                END DO
              END DO
           END IF
!
         ELSE IF (I_AEROSOL == IP_SEASALT_JET) THEN

           IF (L_USE_ARCL(IP_ARCL_SSLT)) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J)=                          &
     &                    (ARCL(LG, N_LAYER+1-I,                        &
     &                     I_ARCL_COMPNTS(IP_ARCL_SSLT_JT))*4.0*PI      &
     &                   *DENSITY_SEA_SALT*(MODE_RADIUS_SS_JET**3.0)    &
     &                   *EXP(4.5*(ALOG(SIGMA_SS_JET))**2.0))           &
     &                   /(3.0*P(L, I)/(R*T(L, I)))
                END DO
              END DO
           ELSE IF (L_USE_SEASALT_DIRECT) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J)=                          &
     &                    (SEA_SALT_JET(LG, N_LAYER+1-I)*4.0*PI         &
     &                   *DENSITY_SEA_SALT*(MODE_RADIUS_SS_JET**3.0)    &
     &                   *EXP(4.5*(ALOG(SIGMA_SS_JET))**2.0))           &
     &                   /(3.0*P(L, I)/(R*T(L, I)))
                END DO
              END DO
           END IF
!
         ELSE IF ((I_AEROSOL == IP_BIOGENIC) .AND. L_USE_BIOGENIC) THEN
     
           DO I=I_TOP_COPY, N_LAYER
             DO L=1, N_PROFILE
                LG=I_GATHER(L)
                AEROSOL_MIX_RATIO(L, I, J)=                             &
     &                   BIOGENIC(LG, N_LAYER+1-I)
             END DO
           END DO
!
         ELSE IF (I_AEROSOL == IP_OCFF_FRESH) THEN
         
           IF (L_USE_ARCL(IP_ARCL_OCFF)) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J)=                           &
                         ARCL(LG, N_LAYER+1-I,                          &
     &                     I_ARCL_COMPNTS(IP_ARCL_OCFF_FR))
                END DO
              END DO
           ELSE IF (L_USE_OCFF_DIRECT) THEN
             DO I=I_TOP_COPY, N_LAYER
               DO L=1, N_PROFILE
                 LG=I_GATHER(L)
                 AEROSOL_MIX_RATIO(L, I, J)=                            &
     &                   FRESH_OCFF(LG, N_LAYER+1-I)
               END DO
             END DO
           END IF
!
         ELSE IF (I_AEROSOL == IP_OCFF_AGED) THEN
         
           IF (L_USE_ARCL(IP_ARCL_OCFF)) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J)=                           &
                         ARCL(LG, N_LAYER+1-I,                          &
     &                     I_ARCL_COMPNTS(IP_ARCL_OCFF_AG))
                END DO
              END DO
           ELSE IF (L_USE_OCFF_DIRECT) THEN
             DO I=I_TOP_COPY, N_LAYER
               DO L=1, N_PROFILE
                 LG=I_GATHER(L)
                 AEROSOL_MIX_RATIO(L, I, J)=                            &
     &                   AGED_OCFF(LG, N_LAYER+1-I)
               END DO
             END DO
           END IF
!
         ELSE IF (I_AEROSOL == IP_DELTA) THEN
         
           IF (L_USE_ARCL(IP_ARCL_DLTA)) THEN
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J)=                           &
                         ARCL(LG, N_LAYER+1-I,                          &
      &                    I_ARCL_COMPNTS(IP_ARCL_DLTA_DL))
                END DO
              END DO
           
           END IF
!         
         ELSE
!
!           The options to the radiation code do not require this
!           aerosol to be considered: its mixing ratio is set to 0.
!           This block of code should not normally be executed,
!           but may be required for ease of including modifications.
!
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J)=0.0E+00
               ENDDO
            ENDDO
!
!
         ENDIF
!
!        If using an extra top layer extrapolate the mixing ratio
!        from the adjacent layer, except in the case of climatological
!        aerosols which are specifically set.
         IF (L_EXTRA_TOP) THEN
           IF ( ( (I_AEROSOL == IP_ACCUM_SULPHATE).AND.                 &
     &            (L_USE_SULPC_DIRECT.OR.L_USE_ARCL(IP_ARCL_SULP)) )    &
     &          .OR.                                                    &
     &          ( (I_AEROSOL == IP_AITKEN_SULPHATE).AND.                &
     &            (L_USE_SULPC_DIRECT.OR.L_USE_ARCL(IP_ARCL_SULP)) )    &
     &          .OR.                                                    &
     &          ( (I_AEROSOL == IP_DUST_1).AND.                         &
     &            (L_USE_DUST.OR.L_USE_ARCL(IP_ARCL_DUST)) )            &
     &          .OR.                                                    &
     &          ( (I_AEROSOL == IP_DUST_2).AND.                         &
     &            (L_USE_DUST.OR.L_USE_ARCL(IP_ARCL_DUST)) )            &
     &          .OR.                                                    &
     &          ( (I_AEROSOL == IP_DUST_3).AND.                         &
     &            (L_USE_DUST.OR.L_USE_ARCL(IP_ARCL_DUST)) )            &
     &          .OR.                                                    &
     &          ( (I_AEROSOL == IP_DUST_4).AND.                         &
     &            (L_USE_DUST.OR.L_USE_ARCL(IP_ARCL_DUST)) )            &
     &          .OR.                                                    &
     &          ( (I_AEROSOL == IP_DUST_5).AND.                         &
     &            (L_USE_DUST.OR.L_USE_ARCL(IP_ARCL_DUST)) )            &
     &          .OR.                                                    &
     &          ( (I_AEROSOL == IP_DUST_6).AND.                         &
     &            (L_USE_DUST.OR.L_USE_ARCL(IP_ARCL_DUST)) )            &
     &          .OR.                                                    &
     &          ( (I_AEROSOL == IP_SEASALT_FILM).AND.                   &
     &            (L_USE_SEASALT_DIRECT.OR.L_USE_ARCL(IP_ARCL_SSLT)))   &
     &          .OR.                                                    &
     &          ( (I_AEROSOL == IP_SEASALT_JET).AND.                    &
     &            (L_USE_SEASALT_DIRECT.OR.L_USE_ARCL(IP_ARCL_SSLT)))   &
     &          .OR.                                                    &
     &          ( (I_AEROSOL == IP_FRESH_SOOT).AND.                     &
     &            (L_USE_SOOT_DIRECT.OR.L_USE_ARCL(IP_ARCL_BLCK)) )     &
     &          .OR.                                                    &
     &          ( (I_AEROSOL == IP_AGED_SOOT).AND.                      &
     &            (L_USE_SOOT_DIRECT.OR.L_USE_ARCL(IP_ARCL_BLCK)) )     &
     &          .OR.                                                    &
     &          ( (I_AEROSOL == IP_BIOMASS_1).AND.                      &
     &            (L_USE_BMASS_DIRECT.OR.L_USE_ARCL(IP_ARCL_BIOM)) )    &
     &          .OR.                                                    &
     &          ( (I_AEROSOL == IP_BIOMASS_2).AND.                      &
     &            (L_USE_BMASS_DIRECT.OR.L_USE_ARCL(IP_ARCL_BIOM)) )    &
     &          .OR.                                                    &
     &          ( (I_AEROSOL == IP_BIOGENIC).AND. L_USE_BIOGENIC )      &
     &          .OR.                                                    &
     &          ( (I_AEROSOL == IP_OCFF_FRESH).AND.                     &
     &            (L_USE_OCFF_DIRECT.OR.L_USE_ARCL(IP_ARCL_OCFF)) )     &
     &          .OR.                                                    &
     &          ( (I_AEROSOL == IP_OCFF_AGED).AND.                      &
     &            (L_USE_OCFF_DIRECT.OR.L_USE_ARCL(IP_ARCL_OCFF)) )     &
     &          .OR.                                                    &
     &          ( (I_AEROSOL == IP_DELTA) )                             &
     &       ) THEN
             DO L=1, N_PROFILE
               AEROSOL_MIX_RATIO(L, 1, J)                               &
     &           =AEROSOL_MIX_RATIO(L, 2, J)
             ENDDO
           ENDIF
         ENDIF
!
      ENDDO
!
! Only use aerosols in the boundary layer for the mesoscale model
      IF (L_MURK_RAD) THEN
        INV17=1./1.7
        DO J=1, N_AEROSOL
          DO I=(NLEVS-N_LEVELS_BL+1),NLEVS
            DO L=1, N_PROFILE
              LG=I_GATHER(L)
! Note that aerosol mass mixing ratios are in units of 1E9 kg/kg
              AEROSOL_MIX_RATIO(L,I,J) = AEROSOL_MIX_RATIO(L,I,J)       &
     &          + AERO_MESO(LG,NLEVS+1-I)*1E-9*MESO_FRAC(J)*INV17
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!
!
      RETURN
      END SUBROUTINE R2_SET_AEROSOL_FIELD
#endif
