#if defined(A70_1C) || defined(A70_1Z) \
 || defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   This module defines the elements of the structure defining
!   algorithmic control of the radiation code.
!
! Current Code Owner: J.-C. Thelen
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   6.2   20/08/04   Original code.
!                  J.-C. Thelen
!
! Code Description:
!   Language: FORTRAN 90
!
!- End of header
!
MODULE dec_spec
!
  TYPE spectrum
!
!
!   Sizes of arrays:
!
    INTEGER :: npd_band
!     Number of spectral bands
    INTEGER :: npd_exclude
!     Numer of excluded bands
    INTEGER :: npd_esft_term
!     Number of esft terms
    INTEGER :: npd_type
!     Number of data types
    INTEGER :: npd_species
!     Number of gaseous species
    INTEGER :: npd_scale_fnc
!     Number of scaling functions
    INTEGER :: npd_scale_variable
!     Number of scaling variables
    INTEGER :: npd_surface
!     Number of surface types in LW spectrum
    INTEGER :: npd_albedo_parm
!     Number of albedo parameters in LW spectrum
    INTEGER :: npd_continuum
!     Number of continua
    INTEGER :: npd_drop_type
!     Number of drop types
    INTEGER :: npd_ice_type
!     Number of ice crystal types
    INTEGER :: npd_aerosol_species
!     Number of aerosol species
    INTEGER :: npd_thermal_coeff
!     Number of thermal coefficients
    INTEGER :: npd_cloud_parameter
!     Number of cloud parameters
    INTEGER :: npd_humidities
!     Number of humidities
    INTEGER :: npd_aod_wavel
!     Number of wavelengths for aod in LW spectrum
    INTEGER :: npd_phase_term
!     Number of terms in the phase function
!
! General fields:
!
    LOGICAL, Pointer ::    l_present(:)
!     Flag for types of data present
!
! Properties of the spectral bands:
!
    INTEGER :: n_band
!     Number of spectral bands
!
    REAL, Pointer :: wave_length_short(:)
!     Shorter wavelength limits
    REAL, Pointer :: wave_length_long(:)
!     Longer wavelength limits
!
!
!
! Exclusion of specific bands from parts of the spectrum:
!
    INTEGER, Pointer :: n_band_exclude(:)
!     Number of excluded bands within each spectral band
    INTEGER, Pointer :: index_exclude(:, :)
!     Indices of excluded bands
!
!
!
! Fields for the solar flux:
!
    REAL, Pointer :: solar_flux_band(:)
!     Fraction of the incident solar flux in each band
!
!
!
! Fields for rayleigh scattering:
!
    REAL, Pointer :: rayleigh_coefficient(:)
!     Rayleigh coefficients
!
!
!
! Fields for gaseous absorption:
!
    INTEGER :: n_absorb
!     Number of absorbers
    INTEGER, Pointer :: n_band_absorb(:)
!     Number of absorbers in each band
    INTEGER, Pointer :: index_absorb(:, :)
!     List of absorbers in each band
    INTEGER, Pointer :: type_absorb(:)
!     Types of each gas in the spectral file
    INTEGER, Pointer :: i_band_esft(:, :)
!     Number of esft terms in band for each gas
    INTEGER, Pointer :: i_scale_esft(:, :)
!     Type of esft scaling
    INTEGER, Pointer :: i_scale_fnc(:, :)
!     Type of scaling function
!
    REAL, Pointer :: k_esft(:, :, :)
!     ESFT exponents
    REAL, Pointer :: w_esft(:, :, :)
!     ESFT weights
    REAL, Pointer :: scale_vector(:, :, :, :)
!     Scaling parameters for each absorber and term
    REAL, Pointer :: p_reference(:, :)
!     Reference pressure for scaling function
    REAL, Pointer :: t_reference(:, :)
!     Reference temperature for scaling function
!
!
!
! Representation of the Planckian:
!
    INTEGER :: n_deg_fit
!     Degree of thermal polynomial
!
    REAL, Pointer :: thermal_coefficient(:, :)
!     Coefficients in polynomial fit to source function
    REAL :: t_ref_planck
!     Planckian reference temperature
!
!
! Surface properties
!
    INTEGER, pointer :: i_spec_surface(:)
!     Method of specifying properties of surface
    INTEGER, pointer :: n_dir_albedo_fit(:)
!     Number of parameters fitting the direct albedo

    LOGICAL, pointer :: l_surface(:)
!     Surface types included

    REAL, pointer :: surface_albedo(:,:)
!     Surface albedos
    REAL, pointer :: direct_albedo_parm(:,:,:)
!     Coefficients for fitting the direct albedo
    REAL, pointer :: emissivity_ground(:,:)
!     Surface emissivities

!
! Fields for continua:
!
    INTEGER, Pointer :: n_band_continuum(:)
!     Number of continua in each band
    INTEGER, Pointer :: index_continuum(:, :)
!     List of continua continuua in each band
    INTEGER :: index_water
!     Index of water vapour
    INTEGER, Pointer :: i_scale_fnc_cont(:, :)
!     Type of scaling function for continuum
!
    REAL, Pointer :: k_continuum(:, :)
!     Grey extinction coefficients for continuum
    REAL, Pointer :: scale_continuum(:, :, :)
!     Scaling parameters for continuum
    REAL, Pointer :: p_ref_continuum(:, :)
!     Reference pressure for scaling of continuum
    REAL, Pointer :: t_ref_continuum(:, :)
!     Reference temperature for scaling of continuum
!
!
!
! Fields for water droplets:
!
    INTEGER, Pointer :: i_drop_parametrization(:)
!     Parametrization type of droplets
!
    LOGICAL, Pointer :: l_drop_type(:)
!     Types of droplet present
!
    INTEGER, Pointer :: n_drop_phf_term(:)
!     Number of moments in the phase function for droplets
    REAL, Pointer :: drop_parameter_list(:, :, :)
!     Parameters used to fit optical properties of clouds
    REAL, Pointer :: drop_parm_min_dim(:)
!     Minimum dimension permissible in the parametrization
    REAL, Pointer :: drop_parm_max_dim(:)
!     Maximum dimension permissible in the parametrization
!
!
!
! Fields for aerosols:
!
    INTEGER :: n_aerosol
!     Number of species of aerosol
    INTEGER, Pointer :: type_aerosol(:)
!     Types of aerosols
    INTEGER, Pointer :: i_aerosol_parametrization(:)
!     Parametrization of aerosols
    INTEGER, Pointer :: n_aerosol_phf_term(:)
!     Number of terms in the phase function
    INTEGER, Pointer :: nhumidity(:)
!     Numbers of humidities
!
    LOGICAL, Pointer :: l_aerosol_species(:)
!     Aerosol species included
!
    REAL, Pointer :: aerosol_absorption(:, :, :)
!     Absorption by aerosols
    REAL, Pointer :: aerosol_scattering(:, :, :)
!     Scattering by aerosols
    REAL, Pointer :: aerosol_phase_fnc(:, :, :, :)
!     Phase function of aerosols
    REAL, Pointer :: humidities(:, :)
!     Humidities for components
!
!
!
! Fields for ice crystals:
!
    INTEGER, Pointer :: i_ice_parametrization(:)
!     Types of parametrization of ice crystals
!
    LOGICAL, Pointer :: l_ice_type(:)
!     Types of ice crystal present
!
    INTEGER, Pointer :: n_ice_phf_term(:)
!     Number of moments in the phase function for ice crystals
    REAL, Pointer :: ice_parameter_list(:, :, :)
!     Parameters used to fit single scattering of ice crystals
    REAL, Pointer :: ice_parm_min_dim(:)
!     Minimum dimension permissible in the parametrization
    REAL, Pointer :: ice_parm_max_dim(:)
!     Maximum dimension permissible in the parametrization
!
!
!
! Fields for Doppler broadening:
!
    LOGICAL, Pointer :: l_doppler_present(:)
!     Flag for Doppler broadening for each species
!
    REAL, Pointer :: doppler_correction(:)
!     Doppler correction terms
!
!
!
! Fields for aerosol optical depth:

    INTEGER :: n_aod_wavel
!     number of wavelengths
    INTEGER, Pointer :: i_aod_type(:)
!     relationship between aerosol component and type
    REAL, Pointer :: aod_absorption(:, :, :)
!     monochromatic specific absorption coefficient
    REAL, Pointer :: aod_scattering(:, :, :)
!     monochromatic specific scattering coefficient

  END TYPE spectrum
!
!
!
END MODULE dec_spec
#endif
