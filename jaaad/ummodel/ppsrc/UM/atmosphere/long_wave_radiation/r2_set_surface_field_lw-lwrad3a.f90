
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
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!  5.1    04/04/00  Tolerances replaced by F90 intrinsics.
!                                   (J. M. Edwards)
!
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SET_SURFACE_FIELD_LW(                               &
     &     N_PROFILE, i_gather, N_BAND                                  &
     &   , I_SURFACE, I_SPEC_SURFACE, L_SURFACE                         &
     &   , EMISSIVITY_FIELD, ALBEDO_FIELD_DIR, ALBEDO_FIELD_DIFF        &
     &   , ALBEDO_SEA_DIR, ALBEDO_SEA_DIFF                              &
     &   , FLANDG                                                       &
     &   , N_FRAC_SOL_POINT, I_FRAC_SOL_POINT, ICE_FRACTION             &
     &   , ice_fraction_g                                               &
     &   , npd_field, NPD_PROFILE, NPD_BAND_LW, NPD_SURFACE_LW          &
     &   )
!
!
!
      USE rad_switches_mod, ONLY: lrad_ctile_fix, LRAD_EMIS_LAND_GEN    &
     & ,RAD_EMIS_LAND_GEN
     
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
     &     npd_field                                                    &
!             Size allocated for input fields
     &   , NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF ATMOSPHERIC PROFILES
     &   , NPD_BAND_LW                                                  &
!             MAXIMUM NUMBER OF SPECTRAL BANDS
     &   , NPD_SURFACE_LW
!             MAXIMUM NUMBER OF SURFACES
!
!     ACTUAL SIZES USED:
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF ATMOSPHERIC PROFILES
     &   , N_BAND
!             NUMBER OF SPECTRAL BANDS
!
      Integer, Intent(IN) :: i_gather(npd_field)
!                              List of points to gather
!
!     Characteristics of points:
      Real, Intent(IN)    :: Ice_fraction(npd_field)
!                              Fraction of sea-ice in sea part
!                              of grid-box.
      Real, Intent(IN)    :: Flandg(npd_field)
!                              Fraction of land
!
!
!     PROPERTIES OF SURFACES
      INTEGER                                                           &
                !, INTENT(OUT)
     &     I_SURFACE(NPD_PROFILE)                                       &
!             TYPES OF SURFACES
     &   , I_SPEC_SURFACE(NPD_SURFACE_LW)
      LOGICAL                                                           &
                !, INTENT(OUT)
     &     L_SURFACE(NPD_SURFACE_LW)
!             FLAGS FOR TYPES OF SURFACES
!
!     Gathered surface fields
      Real, Intent(OUT)    :: Ice_fraction_g(npd_profile)
!                               Gathered fraction of sea-ice in
!                               sea part of grid-box.
!
!     SURFACE PROPERTIES SET.
      REAL                                                              &
                !, INTENT(OUT)
     &     EMISSIVITY_FIELD(NPD_PROFILE, NPD_BAND_LW)                   &
!             EMISSIVITIES OF SURFACES
     &   , ALBEDO_FIELD_DIFF(NPD_PROFILE, NPD_BAND_LW)                  &
!             DIFFUSE ALBEDO OF SURFACE
     &   , ALBEDO_FIELD_DIR(NPD_PROFILE, NPD_BAND_LW)                   &
!             DIRECT ALBEDO OF SURFACE
     &   , ALBEDO_SEA_DIFF(NPD_PROFILE, NPD_BAND_LW)                    &
!             DIFFUSE ALBEDO OF OPEN SEA
     &   , ALBEDO_SEA_DIR(NPD_PROFILE, NPD_BAND_LW)
!             DIRECT ALBEDO OF OPEN SEA
!
!     VARIABLES CONCERNED WITH FRACTIONAL SEA ICE
      INTEGER                                                           &
                !, INTENT(OUT)
     &     N_FRAC_SOL_POINT                                             &
!             NUMBER OF POINTS WITH FRACTIONAL ICE/LAND COVER
     &   , I_FRAC_SOL_POINT(NPD_PROFILE)
!             INDICES OF POINTS WITH FRACTIONAL ICE/LAND COVER
!
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
      Integer :: lg
!                  Gathered index of grid-box
!
      Real :: sea_fraction
!               Fraction of open sea in the grid-box
!
!
!
!     OVERRIDE ANY SURFACE PROERTIES READ IN FROM THE SPECTRAL FILE.
      DO L=1, N_PROFILE
         I_SURFACE(L)=1
      ENDDO
      L_SURFACE(1)=.TRUE.
      I_SPEC_SURFACE(1)=IP_SURFACE_INTERNAL
!
!     Set the radiative characteristics of the surface. Aggregate
!     land surface emissivity values are now allowed to be non-unity.
!
      If (LRAD_EMIS_LAND_GEN) Then
      DO I=1, N_BAND
         DO L=1, N_PROFILE
            EMISSIVITY_FIELD(L, I)=flandg(i_gather(l))                  &
     &      *RAD_EMIS_LAND_GEN + (1.0 - flandg(i_gather(l))) * 1.0 
            ALBEDO_FIELD_DIFF(L, I)=1.0 - emissivity_field(l, i)
            ALBEDO_FIELD_DIR(L, I)=0.0E+00
            ALBEDO_SEA_DIFF(L, I)=0.0E+00
            ALBEDO_SEA_DIR(L, I)=0.0E+00
         ENDDO
      ENDDO
      Else
      DO I=1, N_BAND
         DO L=1, N_PROFILE
            EMISSIVITY_FIELD(L, I)=1.0E+00
            ALBEDO_FIELD_DIFF(L, I)=0.0E+00
            ALBEDO_FIELD_DIR(L, I)=0.0E+00
            ALBEDO_SEA_DIFF(L, I)=0.0E+00
            ALBEDO_SEA_DIR(L, I)=0.0E+00
         ENDDO
      ENDDO
      Endif
!
!     SET THE FRACTIONAL OPEN SEA COVERAGE. POINTS ARE REQUIRED WHERE
!     THIS IS NEITHER 0 NOR 1.
      n_frac_sol_point=0
      If (lrad_ctile_fix) Then
        Do l=1, n_profile
          lg=i_gather(l)
          Ice_fraction_g(l)=Ice_fraction(lg)
          sea_fraction = ( 1.0 - flandg(lg) )                           &
     &                 * ( 1.0 - Ice_fraction(lg) )
          If ( sea_fraction * (1.0 - sea_fraction) > 0.0 ) Then
            n_frac_sol_point=n_frac_sol_point+1
            i_frac_sol_point(n_frac_sol_point)=l
          Endif
        Enddo
      Else
        Do l=1, n_profile
          lg=i_gather(l)
          Ice_fraction_g(l)=Ice_fraction(lg)
          sea_fraction = ( 1.0 - flandg(lg) )                           &
     &                 * ( 1.0 - Ice_fraction(lg) )
          If ( sea_fraction * (1.0 - sea_fraction) >                    &
     &         epsilon(sea_fraction) ) Then
            n_frac_sol_point=n_frac_sol_point+1
            i_frac_sol_point(n_frac_sol_point)=l
          Endif
        Enddo
      Endif
!
!
!
      RETURN
      END SUBROUTINE R2_SET_SURFACE_FIELD_LW
