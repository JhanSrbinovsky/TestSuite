#if defined(A70_1C) || defined(A70_1Z)
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read a shortwave spectral namelist.
!
! Purpose:
!   To read a shortwave namelist into a spectral array.
!
! Method:
!   The spectrum is read into the dynamically allocated array
!   and then reduced to a more manageable size.
!
! Current Owner of Code: J. M. Edwards
!
! History:
! Version     Date          Comment
!   6.2     24-10-05   New deck based on the 3A version of SPIN3A
!                      with added functionality for radiative
!                      timestepping, radiative forcing, orographic
!                      effects and radiance calculations.
!                                                 (J.-C. Thelen)
!
!   6.2     16-11-05   Argument list of R2_COMPRESS_SPECTRUM
!                      modified to match changes in the LW.
!                                                 (N Bellouin)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to read a longwave spectral namelist.
!
! Purpose:
!   To read a longwave namelist into a spectral array.
!
! Method:
!   The spectrum is read into the dynamically allocated array
!   and then reduced to a more manageable size.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!   6.2     24-10-05   New deck based on the 3A version of SPIN3A
!                      with added functionality for radiative
!                      timestepping, radiative forcing, orographic
!                      effects and radiance calculations.
!                                                 (J.-C. Thelen)
!   6.2     15-11-05   Added support of block 15 (aerosol
!                      optical depth) in the LW spectral file.
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to transfer spectrum to reduced array.
!
! Purpose:
!       Spectral data from the large dynamically allocated array
!       are transferred to the reduced array.
!
! Method:
!       Elements are copied across.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!   6.2     24-10-05   New deck based on the 3A version of SPIN3A
!                      with added functionality for radiative
!                      timestepping, radiative forcing, orographic
!                      effects and radiance calculations.
!                                                 (J.-C. Thelen)
!   6.2     15-11-05   Added code for aerosol optical depth (block
!                      15 of the LW spectral file).
!                                               (N. Bellouin)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------

!+ Subroutine to allocate a reduced spectrum.
!
! Purpose:
!   To allocate space for the arrays of the reduced spectrum.
!
! Method:
!   The sizes of the arrays have already been selected. Here we
!   set aside the necessary space.
!
! Current owner of code: J.-C. Thelen
!
! History:
!       Version         Date                    Comment
!       6.2             25/10/05                Original code
!                                               (J.-C. Thelen)
! Description of code:
!   FORTRAN 90
!
!- ---------------------------------------------------------------------
      subroutine r2_allocate_spectrum(spectrum_a)
!
!
!
!
!     Modules used:
      use dec_spec
!
      implicit none
!
!
!     Define the spectrum.
      type (spectrum) spectrum_a
!
!
!
!     The dimensions for the arrays have been set in the calling
!     routine. Space for each array is now allocated and initialised
!     to zero.
!
      allocate(spectrum_a%l_present(                                    &
     &    0: spectrum_a%npd_type                                        &
     &  ))
      spectrum_a%l_present(:)=.false.
!
      allocate(spectrum_a%wave_length_short(                            &
     &    spectrum_a%npd_band                                           &
     &  ))
      spectrum_a%wave_length_short(:)=0.0

      allocate(spectrum_a%wave_length_long(                             &
     &    spectrum_a%npd_band                                           &
     &  ))
      spectrum_a%wave_length_long(:)=0.0

      allocate(spectrum_a%n_band_exclude(                               &
     &    spectrum_a%npd_band                                           &
     &  ))
      spectrum_a%n_band_exclude(:)=0.0

      allocate(spectrum_a%index_exclude(                                &
     &    spectrum_a%npd_exclude                                        &
     &  , spectrum_a%npd_band                                           &
     &  ))
      spectrum_a%index_exclude(:,:)=0.0
!
      allocate(spectrum_a%solar_flux_band(                              &
     &    spectrum_a%npd_band                                           &
     &  ))
      spectrum_a%solar_flux_band(:)=0.0

      allocate(spectrum_a%rayleigh_coefficient(                         &
     &    spectrum_a%npd_band                                           &
     &  ))
      spectrum_a%rayleigh_coefficient(:)=0.0
!
      allocate(spectrum_a%n_band_absorb(                                &
     &    spectrum_a%npd_band                                           &
     &  ))
      spectrum_a%n_band_absorb(:)=0.0

      allocate(spectrum_a%index_absorb(                                 &
     &    spectrum_a%npd_species                                        &
     &  , spectrum_a%npd_band                                           &
     &  ))
      spectrum_a%index_absorb(:,:)=0.0

      allocate(spectrum_a%type_absorb(                                  &
     &    spectrum_a%npd_species                                        &
     &  ))
      spectrum_a%type_absorb(:)=0.0

      allocate(spectrum_a%i_band_esft(                                  &
     &    spectrum_a%npd_band                                           &
     &  , spectrum_a%npd_species                                        &
     &  ))
      spectrum_a%i_band_esft(:,:)=0.0

      allocate(spectrum_a%i_scale_esft(                                 &
     &    spectrum_a%npd_band                                           &
     &  , spectrum_a%npd_species                                        &
     &  ))
      spectrum_a%i_scale_esft(:,:)=0.0

      allocate(spectrum_a%i_scale_fnc(                                  &
     &    spectrum_a%npd_band                                           &
     &  , spectrum_a%npd_species                                        &
     &  ))
      spectrum_a%i_scale_fnc(:,:)=0.0

      allocate(spectrum_a%k_esft(                                       &
     &    spectrum_a%npd_esft_term                                      &
     &  , spectrum_a%npd_band                                           &
     &  , spectrum_a%npd_species                                        &
     &  ))
      spectrum_a%k_esft(:,:,:)=0.0

      allocate(spectrum_a%w_esft(                                       &
     &    spectrum_a%npd_esft_term                                      &
     &  , spectrum_a%npd_band                                           &
     &  , spectrum_a%npd_species                                        &
     &  ))
      spectrum_a%w_esft(:,:,:)=0.0

      allocate(spectrum_a%scale_vector(                                 &
     &    spectrum_a%npd_scale_variable                                 &
     &  , spectrum_a%npd_esft_term                                      &
     &  , spectrum_a%npd_band                                           &
     &  , spectrum_a%npd_species                                        &
     &  ))
      spectrum_a%scale_vector(:,:,:,:)=0.0

      allocate(spectrum_a%p_reference(                                  &
     &    spectrum_a%npd_species                                        &
     &  , spectrum_a%npd_band                                           &
     &  ))
      spectrum_a%p_reference(:,:)=0.0

      allocate(spectrum_a%t_reference(                                  &
     &    spectrum_a%npd_species                                        &
     &  , spectrum_a%npd_band                                           &
     &  ))
      spectrum_a%t_reference(:,:)=0.0
!
      allocate(spectrum_a%thermal_coefficient(                          &
     &    0: spectrum_a%npd_thermal_coeff-1                             &
     &  , spectrum_a%npd_band                                           &
     &  ))
      spectrum_a%thermal_coefficient(:,:)=0.0
!
      allocate( spectrum_a%i_spec_surface(                              &
     &    spectrum_a%npd_surface                                        &
     &  ))
      spectrum_a%i_spec_surface(:)=0.0

      allocate( spectrum_a%n_dir_albedo_fit(                            &
     &    spectrum_a%npd_surface                                        &
     &  ))
      spectrum_a%n_dir_albedo_fit(:)=0.0

      allocate( spectrum_a%l_surface(                                   &
     &    spectrum_a%npd_surface                                        &
     &  ))
      spectrum_a%l_surface(:)=.false.

      allocate( spectrum_a%surface_albedo(                              &
     &    spectrum_a%npd_band                                           &
     &  , spectrum_a%npd_surface                                        &
     &  ))
      spectrum_a%surface_albedo(:,:)=0.0

      allocate( spectrum_a%direct_albedo_parm(                          &
     &    0:spectrum_a%npd_albedo_parm                                  &
     &  , spectrum_a%npd_band                                           &
     &  , spectrum_a%npd_surface                                        &
     & ))
      spectrum_a%direct_albedo_parm(:,:,:)=0.0

      allocate( spectrum_a%emissivity_ground(                           &
     &    spectrum_a%npd_band                                           &
     &  , spectrum_a%npd_surface                                        &
     & ))
      spectrum_a%emissivity_ground(:,:)=0.0
!
      allocate(spectrum_a%n_band_continuum(                             &
     &    spectrum_a%npd_band                                           &
     &  ))
      spectrum_a%n_band_continuum(:)=0.0

      allocate(spectrum_a%index_continuum(                              &
     &    spectrum_a%npd_band                                           &
     &  , spectrum_a%npd_continuum                                      &
     &  ))
      spectrum_a%index_continuum(:,:)=0.0

      allocate(spectrum_a%i_scale_fnc_cont(                             &
     &    spectrum_a%npd_band                                           &
     &  , spectrum_a%npd_continuum                                      &
     &  ))
      spectrum_a%i_scale_fnc_cont(:,:)=0.0

      allocate(spectrum_a%k_continuum(                                  &
     &    spectrum_a%npd_band                                           &
     &  , spectrum_a%npd_continuum                                      &
     &  ))
      spectrum_a%k_continuum(:,:)=0.0

      allocate(spectrum_a%scale_continuum(                              &
     &    spectrum_a%npd_scale_variable                                 &
     &  , spectrum_a%npd_band                                           &
     &  , spectrum_a%npd_continuum                                      &
     &  ))
      spectrum_a%scale_continuum(:,:,:)=0.0

      allocate(spectrum_a%p_ref_continuum(                              &
     &    spectrum_a%npd_continuum                                      &
     &  , spectrum_a%npd_band                                           &
     &  ))
      spectrum_a%p_ref_continuum(:,:)=0.0

      allocate(spectrum_a%t_ref_continuum(                              &
     &    spectrum_a%npd_continuum                                      &
     &  , spectrum_a%npd_band                                           &
     &  ))
      spectrum_a%t_ref_continuum(:,:)=0.0
!
      allocate(spectrum_a%i_drop_parametrization(                       &
     &    spectrum_a%npd_drop_type                                      &
     &  ))
      spectrum_a%i_drop_parametrization(:)=0.0

      allocate(spectrum_a%l_drop_type(                                  &
     &    spectrum_a%npd_drop_type                                      &
     &  ))
      spectrum_a%l_drop_type(:)=.false.

      allocate(spectrum_a%n_drop_phf_term(                              &
     &    spectrum_a%npd_drop_type                                      &
     &  ))
      spectrum_a%n_drop_phf_term(:)=0.0

      allocate(spectrum_a%drop_parameter_list(                          &
     &    spectrum_a%npd_cloud_parameter                                &
     &  , spectrum_a%npd_band                                           &
     &  , spectrum_a%npd_drop_type                                      &
     &  ))
      spectrum_a%drop_parameter_list(:,:,:)=0.0

      allocate(spectrum_a%drop_parm_min_dim(                            &
     &    spectrum_a%npd_drop_type                                      &
     &  ))
      spectrum_a%drop_parm_min_dim(:)=0.0

      allocate(spectrum_a%drop_parm_max_dim(                            &
     &    spectrum_a%npd_drop_type                                      &
     &  ))
      spectrum_a%drop_parm_max_dim(:)=0.0
!
      allocate(spectrum_a%type_aerosol(                                 &
     &    spectrum_a%npd_aerosol_species                                &
     &  ))
      spectrum_a%type_aerosol(:)=0.0

      allocate(spectrum_a%i_aerosol_parametrization(                    &
     &    spectrum_a%npd_aerosol_species                                &
     &  ))
      spectrum_a%i_aerosol_parametrization(:)=0.0

      allocate(spectrum_a%n_aerosol_phf_term(                           &
     &    spectrum_a%npd_aerosol_species                                &
     &  ))
      spectrum_a%n_aerosol_phf_term(:)=0.0

      allocate(spectrum_a%nhumidity(                                    &
     &    spectrum_a%npd_aerosol_species                                &
     &  ))
      spectrum_a%nhumidity(:)=0.0

      allocate(spectrum_a%l_aerosol_species(                            &
     &    spectrum_a%npd_aerosol_species                                &
     &  ))
      spectrum_a%l_aerosol_species(:)=.false.

      allocate(spectrum_a%aerosol_absorption(                           &
     &    spectrum_a%npd_humidities                                     &
     &  , spectrum_a%npd_aerosol_species                                &
     &  , spectrum_a%npd_band                                           &
     &  ))
      spectrum_a%aerosol_absorption(:,:,:)=0.0

      allocate(spectrum_a%aerosol_scattering(                           &
     &    spectrum_a%npd_humidities                                     &
     &  , spectrum_a%npd_aerosol_species                                &
     &  , spectrum_a%npd_band                                           &
     &  ))
      spectrum_a%aerosol_scattering(:,:,:)=0.0

      allocate(spectrum_a%aerosol_phase_fnc(                            &
     &    spectrum_a%npd_humidities                                     &
     &  , spectrum_a%npd_phase_term                                     &
     &  , spectrum_a%npd_aerosol_species                                &
     &  , spectrum_a%npd_band                                           &
     &  ))
      spectrum_a%aerosol_phase_fnc(:,:,:,:)=0.0

      allocate(spectrum_a%humidities(                                   &
     &    spectrum_a%npd_humidities                                     &
     &  , spectrum_a%npd_aerosol_species                                &
     &  ))
      spectrum_a%humidities(:,:)=0.0
!
      allocate(spectrum_a%i_ice_parametrization(                        &
     &    spectrum_a%npd_ice_type                                       &
     &  ))
      spectrum_a%i_ice_parametrization(:)=0.0

      allocate(spectrum_a%l_ice_type(                                   &
     &    spectrum_a%npd_ice_type                                       &
     &  ))
      spectrum_a%l_ice_type(:)=.false.

      allocate(spectrum_a%n_ice_phf_term(                               &
     &    spectrum_a%npd_ice_type                                       &
     &  ))
      spectrum_a%n_ice_phf_term(:)=0.0

      allocate(spectrum_a%ice_parameter_list(                           &
     &    spectrum_a%npd_cloud_parameter                                &
     &  , spectrum_a%npd_band                                           &
     &  , spectrum_a%npd_ice_type                                       &
     &  ))
      spectrum_a%ice_parameter_list(:,:,:)=0.0

      allocate(spectrum_a%ice_parm_min_dim(                             &
     &    spectrum_a%npd_ice_type                                       &
     &  ))
      spectrum_a%ice_parm_min_dim(:)=0.0

      allocate(spectrum_a%ice_parm_max_dim(                             &
     &    spectrum_a%npd_ice_type                                       &
     &  ))
      spectrum_a%ice_parm_max_dim(:)=0.0
!
      allocate(spectrum_a%l_doppler_present(                            &
     &    spectrum_a%npd_species                                        &
     &  ))
      spectrum_a%l_doppler_present(:)=.false.

      allocate(spectrum_a%doppler_correction(                           &
     &    spectrum_a%npd_species                                        &
     &  ))
      spectrum_a%doppler_correction(:)=0.0
!
      allocate(spectrum_a%i_aod_type(                                   &
     &    spectrum_a%npd_aerosol_species                                &
     &  ))
      spectrum_a%i_aod_type(:)=0.0

      allocate(spectrum_a%aod_absorption(                               &
     &    spectrum_a%npd_humidities                                     &
     &  , spectrum_a%npd_aerosol_species                                &
     &  , spectrum_a%npd_aod_wavel                                      &
     &  ))
      spectrum_a%aod_absorption(:,:,:)=0.0

      allocate(spectrum_a%aod_scattering(                               &
     &    spectrum_a%npd_humidities                                     &
     &  , spectrum_a%npd_aerosol_species                                &
     &  , spectrum_a%npd_aod_wavel                                      &
     &  ))
      spectrum_a%aod_scattering(:,:,:)=0.0
!
!
      return
      END SUBROUTINE r2_allocate_spectrum
#endif
#endif
