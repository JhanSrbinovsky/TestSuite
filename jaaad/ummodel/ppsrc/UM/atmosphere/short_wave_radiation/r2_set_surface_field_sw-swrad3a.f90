
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set surface fields.
!
! Purpose:
!   The albedos and emissivity of the surface are set.
!
! Method:
!   Straightforward. Though the arrays passed to the code may depend
!   on the spectral band, the input arrays have no spectral dependence.
!
! Current Owner of Code: James Manners
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SET_SURFACE_FIELD_SW(                               &
     &     N_BAND                                                       &
     &   , NLIT, LIST                                                   &
     &   , I_SURFACE, I_SPEC_SURFACE, L_SURFACE                         &
     &   , L_MICROPHYSICS, L_MOSES_II, l_cable, L_CTILE                 &
     &   , L_USE_SPEC_SEA                                               &
     &   , LAND, LAND0P5, OPEN_SEA_ALBEDO                               &
     &   , LAND_ALB, SICE_ALB                                           &
     &   , FLANDG, ICE_FRACTION                                         &
     &   , LAND_ALBEDO, WEIGHT_690NM                                    &
     &   , EMISSIVITY_FIELD, ALBEDO_FIELD_DIR, ALBEDO_FIELD_DIFF        &
     &   , LAND_G, LAND0P5_G, FLANDG_G, ICE_FRACTION_G                  &
     &   , ALBEDO_SEA_DIFF, ALBEDO_SEA_DIR                              &
     &   , NPD_FIELD, NPD_PROFILE, NPD_BAND_SW, NPD_SURFACE_SW          &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED
! SRFSP3A defines permitted methods of specifying the  surface albedo and
! emissivity for two-stream radiation code.

      ! properties specified by surface type
      INTEGER,PARAMETER:: IP_SURFACE_SPECIFIED=1

      ! properties passed into code
      INTEGER,PARAMETER:: IP_SURFACE_INTERNAL=2

      ! direct albedo fitted as polynomial
      INTEGER,PARAMETER:: IP_SURFACE_POLYNOMIAL=3

      ! fit in the functional form used by payne
      INTEGER,PARAMETER:: IP_SURFACE_PAYNE=4
! SRFSP3A end
!
!     DUMMY VARIABLES:
!
!     DIMENSIONS OF ARRAYS:
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_FIELD                                                    &
!             SIZE OF INPUT FIELDS
     &   , NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF ATMOSPHERIC PROFILES
     &   , NPD_BAND_SW                                                  &
!             MAXIMUM NUMBER OF SPECTRAL BANDS
     &   , NPD_SURFACE_SW
!             MAXIMUM NUMBER OF SURFACES
!
!     ACTUAL SIZES USED:
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_BAND
!             NUMBER OF SPECTRAL BANDS
!
!     LIT POINTS:
      INTEGER                                                           &
                !, INTENT(IN)
     &     NLIT                                                         &
!             NUMBER OF LIT POINTS
     &   , LIST(NPD_FIELD)
!             LIST OF SUNLIT POINTS
!
!     PROPERTIES OF SURFACES
      INTEGER                                                           &
                !, INTENT(OUT)
     &     I_SURFACE(NPD_PROFILE)                                       &
!             TYPES OF SURFACES
     &   , I_SPEC_SURFACE(NPD_SURFACE_SW)
      LOGICAL                                                           &
                !, INTENT(OUT)
     &     L_SURFACE(NPD_SURFACE_SW)
!             FLAGS FOR TYPES OF SURFACES
!
!     PHYSICAL PROPERTIES OF SURFACES:
      LOGICAL                                                           &
                !, INTENT(IN)
     &     LAND(NPD_FIELD)                                              &
!             LAND MASK
     &   , LAND0P5(NPD_FIELD)
!             LAND MASK (TRUE if land fraction >0.5)
      REAL                                                              &
                !, INTENT(IN)
     &     OPEN_SEA_ALBEDO(NPD_FIELD, 2)                                &
!             DIFFUSE ALBEDO FIELD
     &   , FLANDG(NPD_FIELD)                                            &
!             LAND FRACTION
     &   , LAND_ALB(NPD_FIELD)                                          &
     &   , SICE_ALB(NPD_FIELD)                                          &
     &   , LAND_ALBEDO(NPD_FIELD,4)                                     &
!             MOSES II LAND SURFACE ALBEDO FIELDS
     &   , WEIGHT_690NM(NPD_BAND_SW)                                    &
!             WEIGHTS FOR EACH BAND FOR REGION BELOW 690 NM
     &   , ICE_FRACTION(NPD_FIELD)
!             FRACTION OF SEA ICE IN SEA PORTION OF GRID BOX
!
!     MISCELLANEOUS INPUTS
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_MICROPHYSICS                                               &
!             FLAG TO CALCULATE MICROPHYSICS
     &   , L_MOSES_II                                                   &
!             FLAG FOR MOSES II LAND SURFACE
     &   , L_cable                                                      &
!             FLAG FOR cable LAND SURFACE
     &   , L_CTILE                                                      &
!             FLAG FOR COASTAL TILING
     &   , L_USE_SPEC_SEA
!             FLAG FOR SPECTRALLY DEPENDENT SEA ALBEDOS
!
!
!     SURFACE PROPERTIES SET.
      REAL                                                              &
                !, INTENT(OUT)
     &     EMISSIVITY_FIELD(NPD_PROFILE, NPD_BAND_SW)                   &
!             EMISSIVITIES OF SURFACES
     &   , ALBEDO_FIELD_DIFF(NPD_PROFILE, NPD_BAND_SW)                  &
!             DIFFUSE ALBEDO OF SURFACE
     &   , ALBEDO_FIELD_DIR(NPD_PROFILE, NPD_BAND_SW)
!             DIRECT ALBEDO OF SURFACE
!
!     GATHERED SURFACE FIELDS
      LOGICAL                                                           &
                !, INTENT(OUT)
     &     LAND_G(NPD_PROFILE)                                          &
!             GATHERED LAND MASK
     &   , LAND0P5_G(NPD_PROFILE)
!             GATHERED LAND MASK (TRUE if land fraction >0.5)
      REAL                                                              &
                !, INTENT(OUT)
     &     ALBEDO_SEA_DIFF(NPD_PROFILE, NPD_BAND_SW)                    &
!             DIFFUSE ALBEDO OF OPEN SEA
     &   , ALBEDO_SEA_DIR(NPD_PROFILE, NPD_BAND_SW)                     &
!             DIRECT ALBEDO OF OPEN SEA
     &   , FLANDG_G(NPD_PROFILE)                                        &
!             GATHERED LAND FRACTION
     &   , ICE_FRACTION_G(NPD_PROFILE)
!             GATHERED SEA-ICE FRACTION IN SEA PORTION OF GRID BOX
!
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
      REAL                                                              &
     &     SpectralSea(N_BAND)
!             Spectrally dependent sea albedos
!
      SpectralSea = 1.0
!
!
      IF (L_USE_SPEC_SEA .AND. N_BAND == 6) THEN
         SpectralSea = (/ 0.000000, 1.444205, 1.799420,                 &
     &                   0.5470291, 0.000000, 0.000000 /)
      ENDIF
!
!
!     OVERRIDE ANY SURFACE PROERTIES READ IN FROM THE SPECTRAL FILE.
      DO L=1, NLIT
         I_SURFACE(L)=1
      ENDDO
      L_SURFACE(1)=.TRUE.
      I_SPEC_SURFACE(1)=IP_SURFACE_INTERNAL
!
!
!        GATHER THE ARRAY OF SURFACE FLAGS
!        (required unconditionally from vn5.1 since flags used in
!         calculating flux diagnostics)
         DO L=1, NLIT
            LAND_G(L)=LAND(LIST(L))
            LAND0P5_G(L)=LAND0P5(LIST(L))
            FLANDG_G(L)=FLANDG(LIST(L))
            ICE_FRACTION_G(L)=ICE_FRACTION(LIST(L))

         ENDDO
!
!
!     SET THE ALBEDO FIELDS: AN AVERAGE ALBEDO IS REQUIRED OVER WHERE
!     THERE IS A COMBINATION OF OPEN SEA, SEA-ICE AND LAND. SEPARATE
!     ALBEDOS ARE PROVIDED FOR FOR OPEN SEA. BAND-DEPENDENT COPIES
!     OF THE ALBEDOS MUST BE MADE FOR CALCULATING COUPLING FLUXES.
!
      DO I=1, N_BAND
         DO L=1, NLIT
!
            EMISSIVITY_FIELD(L, I)=0.0E+00
!
            ALBEDO_SEA_DIR(L, I)=0.0E+00
            ALBEDO_SEA_DIFF(L, I)=0.0E+00
            ALBEDO_FIELD_DIFF(L, I)=0.0E+00
            ALBEDO_FIELD_DIR(L, I)=0.0E+00

            IF (FLANDG(LIST(L)) <  1.0) THEN
               ALBEDO_FIELD_DIFF(L, I)                                  &
     &            =SICE_ALB(LIST(L))*ICE_FRACTION(LIST(L))              &
     &            +OPEN_SEA_ALBEDO(LIST(L), 2)                          &
     &            *SpectralSea(I)                                       &
     &            *(1.0E+00-ICE_FRACTION(LIST(L)))
               ALBEDO_FIELD_DIR(L, I)                                   &
     &            =SICE_ALB(LIST(L))*ICE_FRACTION(LIST(L))              &
     &            +OPEN_SEA_ALBEDO(LIST(L), 1)                          &
     &            *SpectralSea(I)                                       &
     &            *(1.0E+00-ICE_FRACTION(LIST(L)))
               ALBEDO_SEA_DIR(L, I)=OPEN_SEA_ALBEDO(LIST(L), 1)         &
     &            *SpectralSea(I)
               ALBEDO_SEA_DIFF(L, I)=OPEN_SEA_ALBEDO(LIST(L), 2)        &
     &            *SpectralSea(I)
            ENDIF
            IF (FLANDG(LIST(L)) >  0.0) THEN
               IF ( L_MOSES_II .or. l_cable ) THEN
               IF ( L_CTILE ) THEN
                 ALBEDO_FIELD_DIFF(L,I) =                               &
     &             (1.-FLANDG_G(L))*ALBEDO_FIELD_DIFF(L, I) +           &
     &               FLANDG_G(L)*(WEIGHT_690NM(I)*LAND_ALBEDO(LIST(L),2)&
     &                  + (1. - WEIGHT_690NM(I))*LAND_ALBEDO(LIST(L),4))
                 ALBEDO_FIELD_DIR(L,I) =                                &
     &             (1.-FLANDG_G(L))*ALBEDO_FIELD_DIR(L, I) +            &
     &               FLANDG_G(L)*(WEIGHT_690NM(I)*LAND_ALBEDO(LIST(L),1)&
     &                 + (1. - WEIGHT_690NM(I))*LAND_ALBEDO(LIST(L),3))
               ELSE

                 ALBEDO_FIELD_DIFF(L,I) =                               &
     &                            WEIGHT_690NM(I)*LAND_ALBEDO(LIST(L),2)&
     &                   + (1. - WEIGHT_690NM(I))*LAND_ALBEDO(LIST(L),4)
                 ALBEDO_FIELD_DIR(L,I) =                                &
     &                            WEIGHT_690NM(I)*LAND_ALBEDO(LIST(L),1)&
     &                   + (1. - WEIGHT_690NM(I))*LAND_ALBEDO(LIST(L),3)
               ENDIF
               ELSE
! For non MOSES_II cannot have coastal tiling, therefore
! must be completely land:
               ALBEDO_FIELD_DIFF(L, I)=LAND_ALB(LIST(L))
               ALBEDO_FIELD_DIR(L, I)=LAND_ALB(LIST(L))
               ENDIF
            ENDIF
!
         ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE R2_SET_SURFACE_FIELD_SW
