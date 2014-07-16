#if defined(A70_1C) || defined(A70_1Z)
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
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
      subroutine r2_compress_spectrum(                                  &
!                       Original spectrum
     &    l_present                                                     &
     &  , n_band, wave_length_short , wave_length_long                  &
     &  , n_band_exclude, index_exclude                                 &
     &  , solar_flux_band, rayleigh_coefficient                         &
     &  , n_absorb, n_band_absorb, index_absorb, type_absorb            &
     &  , l_retain_absorb, n_absorb_retain, index_absorb_retain         &
     &  , compressed_index, i_band_esft, k_esft, w_esft, i_scale_esft   &
     &  , i_scale_fnc, scale_vector, p_reference, t_reference           &
     &  , n_deg_fit, thermal_coefficient, t_ref_planck                  &
     &  , i_spec_surface, l_surface, surface_albedo                     &
     &  , n_dir_albedo_fit, direct_albedo_parm, emissivity_ground       &
     &  , n_band_continuum, index_continuum, index_water                &
     &  , k_continuum, i_scale_fnc_cont, scale_continuum                &
     &  , p_ref_continuum, t_ref_continuum                              &
     &  , l_drop_type, i_drop_parametrization                           &
     &  , n_drop_phf_term, drop_parameter_list                          &
     &  , drop_parm_min_dim, drop_parm_max_dim                          &
     &  , l_ice_type, i_ice_parametrization                             &
     &  , n_ice_phf_term, ice_parameter_list                            &
     &  , ice_parm_min_dim, ice_parm_max_dim                            &
     &  , n_aerosol, type_aerosol                                       &
     &  , n_aerosol_retain, index_aerosol_retain                        &
     &  , l_aerosol_species, aerosol_absorption                         &
     &  , aerosol_scattering, n_aerosol_phf_term, aerosol_phase_fnc     &
     &  , nhumidity, humidities, i_aerosol_parametrization              &
     &  , l_doppler_present, doppler_correction, l_use_aod              &
     &  , N_AOD_WAVEL, AOD_ABSORPTION, AOD_SCATTERING, I_AOD_TYPE       &
!                       Reduced spectral array
     &  , spectrum_rd                                                   &
     &  )
!
!
!
!     Modules used:
      use dec_spec
!
      implicit none
!
!
!     ------------------------------------------------------------------
!     Declaration of initial spectrum.
!     ------------------------------------------------------------------
#include "mxsize3a.h"
#include "spdec3a.h"
!
!     ------------------------------------------------------------------
!     Auxiliary variables used to select parts of the initial spectrum
!     ------------------------------------------------------------------
!
      LOGICAL, Intent(IN) :: l_retain_absorb(npd_species)
!           Flags for the retention of gases in the spectral file
      LOGICAL, Intent(IN) :: l_use_aod
!           Aerosol optical depth diagnostics are requested.
      INTEGER, Intent(IN) :: n_absorb_retain
!           Number of absorbers to be retained
      INTEGER, Intent(IN) :: index_absorb_retain(npd_species)
!           Indices of absorbers to be retained
      INTEGER, Intent(IN) :: compressed_index(npd_species)
!           Mapping from old to new indices of absorbers
      INTEGER, Intent(IN) :: n_aerosol_retain
!           Number of aerosols in the initial spectrum to be used
!           In the calculation
      INTEGER, Intent(IN) :: index_aerosol_retain(npd_aerosol_species)
!           Indices of the retained aerosols
!
!
!     ------------------------------------------------------------------
!     Declaration of reduced spectrum.
!     ------------------------------------------------------------------
!
      TYPE (spectrum) :: spectrum_rd
!
!
!
!     Local variables
      integer                                                           &
     &    i                                                             &
!           Loop variable
     &  , j                                                             &
!           Loop variable
     &  , k                                                             &
!           Loop variable
     &  , l                                                             &
!           Loop variable
     &  , ls                                                            &
!           Loop variable
     &  , n_parameter                                                   &
!           Number of parameters in scheme
     &  , i_initial                                                     &
!           Indexing number in initial spectrum
     &  , i_species                                                     &
!           Species of gas
     &  , i_continuum
!           Type of continuum
!
!
#include "sclfnc3a.h"
#include "wclprm3a.h"
#include "iclprm3a.h"
#include "aerprm3a.h"
#include "sclfnd3a.h"
!
!
!
!
!
!     Initailize all blocks of the compressed spectrum to .FALSE.
      do i=0, spectrum_rd%npd_type
        spectrum_rd%l_present(i)=.false.
      enddo
!
!
!     Proceed through each block of the spectral file transferring
!     the data from the input array to the reduced array.
!
!
!     Block 0:
!
      if (l_present(0)) then
        spectrum_rd%l_present(0)=.true.
        spectrum_rd%n_band=n_band
        spectrum_rd%n_absorb=n_absorb_retain
        spectrum_rd%n_aerosol=n_aerosol_retain
        do i=1, n_absorb_retain
          spectrum_rd%type_absorb(i)                                    &
     &      =type_absorb(index_absorb_retain(i))
        enddo
        do i=1, n_aerosol_retain
          spectrum_rd%type_aerosol(i)                                   &
     &      =type_aerosol(index_aerosol_retain(i))
        enddo
      endif
!
!     Block 1:
      if (l_present(1)) then
        spectrum_rd%l_present(1)=.true.
        do i=1, n_band
          spectrum_rd%wave_length_short(i)=wave_length_short(i)
          spectrum_rd%wave_length_long(i)=wave_length_long(i)
        enddo
      endif
!
!     Block 2:
      if (l_present(2)) then
        spectrum_rd%l_present(2)=.true.
        do i=1, n_band
          spectrum_rd%solar_flux_band(i)=solar_flux_band(i)
        enddo
      endif
!
!     Block 3:
      if (l_present(3)) then
        spectrum_rd%l_present(3)=.true.
        do i=1, n_band
          spectrum_rd%rayleigh_coefficient(i)=rayleigh_coefficient(i)
        enddo
      endif
!
!     Block 4:
      if (l_present(4)) then
        spectrum_rd%l_present(4)=.true.
        do i=1, n_band
          spectrum_rd%n_band_absorb(i)=0
          do j=1, n_band_absorb(i)
            if (l_retain_absorb(index_absorb(j, i))) then
              spectrum_rd%n_band_absorb(i)                              &
     &          =spectrum_rd%n_band_absorb(i)+1
              spectrum_rd%index_absorb(spectrum_rd%n_band_absorb(i), i) &
     &          =compressed_index(index_absorb(j, i))
            endif
          enddo
        enddo
      endif
!
!     Block 5:
      if (l_present(5)) then
        spectrum_rd%l_present(5)=.true.
!
        do i=1, n_band
          do j=1, spectrum_rd%n_band_absorb(i)
!
            i_species=spectrum_rd%index_absorb(j, i)
            i_initial=index_absorb_retain(i_species)
!
            spectrum_rd%i_band_esft(i, i_species)                       &
     &        =i_band_esft(i, i_initial)
            spectrum_rd%i_scale_esft(i, i_species)                      &
     &        =i_scale_esft(i, i_initial)
            spectrum_rd%i_scale_fnc(i, i_species)                       &
     &        =i_scale_fnc(i, i_initial)
            spectrum_rd%p_reference(i_species, i)                       &
     &        =p_reference(i_initial, i)
            spectrum_rd%t_reference(i_species, i)                       &
     &        =t_reference(i_initial, i)

            if (i_scale_fnc(i, i_initial) == ip_scale_fnc_null) then
              n_parameter=0
            else if (i_scale_fnc(i, i_initial) ==                       &
     &               ip_scale_power_law) then
              n_parameter=2
            else if (i_scale_fnc(i, i_initial) ==                       &
     &               ip_scale_power_quad) then
              n_parameter=3
            else if (i_scale_fnc(i, i_initial) ==                       &
     &               ip_scale_doppler_quad) then
              n_parameter=4
            endif

            do k=1, i_band_esft(i, i_initial)
              spectrum_rd%k_esft(k, i, i_species)                       &
     &          =k_esft(k, i, i_initial)
              spectrum_rd%w_esft(k, i, i_species)                       &
     &          =w_esft(k, i, i_initial)
              do l=1, n_scale_variable(i_scale_fnc(i, i_initial))
                spectrum_rd%scale_vector(l, k, i, i_species)            &
     &            =scale_vector(l, k, i, i_initial)
              enddo
            enddo
!
          enddo
        enddo
      endif
!
!     Block 6:
      if (l_present(6)) then
        spectrum_rd%l_present(6)=.true.
        spectrum_rd%n_deg_fit=n_deg_fit
        spectrum_rd%t_ref_planck=t_ref_planck
        do i=1, n_band
          do j=0, n_deg_fit
            spectrum_rd%thermal_coefficient(j, i)                       &
     &        =thermal_coefficient(j, i)
          enddo
        enddo
      endif
!
!     Block 8:
      if (l_present(8)) then
        spectrum_rd%l_present(8)=.true.
        do i=1, n_band
          spectrum_rd%n_band_continuum(i)=n_band_continuum(i)
          do j=1, n_band_continuum(i)
            spectrum_rd%index_continuum(i, j)=index_continuum(i, j)
          enddo
        enddo
!
        spectrum_rd%index_water=0
        do i=1, n_absorb_retain
          if (index_absorb_retain(i) == index_water) then
            spectrum_rd%index_water=i
          endif
        enddo
!
      endif
!
!     Block 9:
      if (l_present(9)) then
        spectrum_rd%l_present(9)=.true.
        do i=1, n_band
          do j=1, n_band_continuum(i)
            i_continuum=index_continuum(i, j)
            spectrum_rd%i_scale_fnc_cont(i, i_continuum)                &
     &        =i_scale_fnc_cont(i, i_continuum)
            spectrum_rd%p_ref_continuum(i_continuum, i)                 &
     &        =p_ref_continuum(i_continuum, i)
            spectrum_rd%t_ref_continuum(i_continuum, i)                 &
     &        =t_ref_continuum(i_continuum, i)

             if (i_scale_fnc(i, i_continuum) == ip_scale_fnc_null) then
               n_parameter=0
             else if (i_scale_fnc(i, i_continuum) ==                    &
     &               ip_scale_power_law) then
               n_parameter=2
             else if (i_scale_fnc(i, i_continuum) ==                    &
     &               ip_scale_power_quad) then
               n_parameter=3
             else if (i_scale_fnc(i, i_continuum) ==                    &
     &               ip_scale_doppler_quad) then
               n_parameter=4
             endif

            spectrum_rd%k_continuum(i, i_continuum)                     &
     &        =k_continuum(i, i_continuum)
            do l=1, n_scale_variable(i_scale_fnc_cont                   &
     &              (i, i_continuum))
              spectrum_rd%scale_continuum(l, i, i_continuum)            &
     &          =scale_continuum(l, i, i_continuum)
            enddo
          enddo
        enddo
      endif
!
!     Block 10:
      if (l_present(10)) then
        spectrum_rd%l_present(10)=.true.
        do i=1, spectrum_rd%npd_drop_type
          if (l_drop_type(i)) then
            spectrum_rd%l_drop_type(i)=.true.
            spectrum_rd%i_drop_parametrization(i)                       &
     &        =i_drop_parametrization(i)
            spectrum_rd%n_drop_phf_term(i)                              &
     &        =n_drop_phf_term(i)
            spectrum_rd%drop_parm_min_dim(i)=drop_parm_min_dim(i)
            spectrum_rd%drop_parm_max_dim(i)=drop_parm_max_dim(i)
            if (i_drop_parametrization(i) == ip_slingo_schrecker) then
              n_parameter=6
            else if (i_drop_parametrization(i) ==                       &
     &                ip_slingo_schr_phf) then
              n_parameter=2*n_drop_phf_term(i)+4
            else if (i_drop_parametrization(i) ==                       &
     &               ip_ackerman_stephens) then
              n_parameter=9
            else if (i_drop_parametrization(i) == ip_drop_pade_2) then
              n_parameter=16
            endif
!
            do j=1, n_parameter
              do k=1, n_band
                spectrum_rd%drop_parameter_list(j, k, i)                &
     &            =drop_parameter_list(j, k, i)
              enddo
            enddo
          else
            spectrum_rd%l_drop_type(i)=.false.
          endif
        enddo
      endif
!
!     Block 11:
      if (l_present(11)) then
        spectrum_rd%l_present(11)=.true.
        do i=1, n_aerosol_retain
          i_initial=index_aerosol_retain(i)
          if (l_aerosol_species(i_initial)) then
            spectrum_rd%l_aerosol_species(i)=.true.
            spectrum_rd%i_aerosol_parametrization(i)                    &
     &        =i_aerosol_parametrization(i_initial)
!
!**** Changed by JCT 26/05/04 ****
!
! Consider Full Phase Function:
!
            if ( (i_aerosol_parametrization(i_initial) ==               &
     &            ip_aerosol_param_phf_dry).or.                         &
     &           (i_aerosol_parametrization(i_initial) ==               &
     &          ip_aerosol_param_phf_moist) ) then
              spectrum_rd%n_aerosol_phf_term(i)                         &
     &          =n_aerosol_phf_term(i_initial)
!
! Consider only the asymmetry
!
            else if ( (i_aerosol_parametrization(i_initial) ==          &
     &            ip_aerosol_param_dry).or.                             &
     &           (i_aerosol_parametrization(i_initial) ==               &
     &          ip_aerosol_param_moist) ) then
              n_aerosol_phf_term(i_initial)=1
              spectrum_rd%n_aerosol_phf_term(i)=1
            endif
!
! Read in Asymmetry or phase function for dry aerosols
!
            if ( (i_aerosol_parametrization(i_initial) ==               &
     &            ip_aerosol_param_dry).or.                             &
     &           (i_aerosol_parametrization(i_initial) ==               &
     &          ip_aerosol_param_phf_dry) ) then
              spectrum_rd%nhumidity(i)=0
              do k=1, n_band
                spectrum_rd%aerosol_absorption(1, i, k)                 &
     &            =aerosol_absorption(1, i_initial, k)
                spectrum_rd%aerosol_scattering(1, i, k)                 &
     &            =aerosol_scattering(1, i_initial, k)
                do ls=1, n_aerosol_phf_term(i_initial)
                  spectrum_rd%aerosol_phase_fnc(1, ls, i, k)            &
     &              =aerosol_phase_fnc(1, ls, i_initial, k)
                enddo
              enddo
!
! Read in Asymmetry or phase function for wet aerosols
!
            else if ( (i_aerosol_parametrization(i_initial) ==          &
     &               ip_aerosol_param_moist).or.                        &
     &                (i_aerosol_parametrization(i_initial) ==          &
     &               ip_aerosol_param_phf_moist) ) then
              spectrum_rd%index_water=index_water
              spectrum_rd%nhumidity(i)=nhumidity(i_initial)
              do j=1, nhumidity(i_initial)
                spectrum_rd%humidities(j, i)=humidities(j, i_initial)
                do k=1, n_band
                  spectrum_rd%aerosol_absorption(j, i, k)               &
     &              =aerosol_absorption(j, i_initial, k)
                  spectrum_rd%aerosol_scattering(j, i, k)               &
     &              =aerosol_scattering(j, i_initial, k)
                  do ls=1, n_aerosol_phf_term(i_initial)
                    spectrum_rd%aerosol_phase_fnc(j, ls, i, k)          &
     &                =aerosol_phase_fnc(j, ls, i_initial, k)
                  enddo
                enddo
              enddo
            endif
!
          else
            spectrum_rd%l_aerosol_species(i)=.false.
          endif
        enddo
      endif
!
!*********************************
!
!     Block 12:
      if (l_present(12)) then
        spectrum_rd%l_present(12)=.true.
        do i=1, spectrum_rd%npd_ice_type
          if (l_ice_type(i)) then
            spectrum_rd%l_ice_type(i)=.true.
            spectrum_rd%ice_parm_min_dim(i)=ice_parm_min_dim(i)
            spectrum_rd%ice_parm_max_dim(i)=ice_parm_max_dim(i)
!
            spectrum_rd%i_ice_parametrization(i)                        &
     &        =i_ice_parametrization(i)
            spectrum_rd%n_ice_phf_term(i)                               &
     &        =n_ice_phf_term(i)
            if (i_ice_parametrization(i) ==                             &
     &          ip_slingo_schrecker_ice) then
              n_parameter=6
            else if (i_ice_parametrization(i) ==                        &
     &               ip_sun_shine_vn2_vis) then
              n_parameter=6
            else if (i_ice_parametrization(i) ==                        &
     &               ip_sun_shine_vn2_ir) then
              n_parameter=0
            else if (i_ice_parametrization(i) == ip_ice_adt) then
              n_parameter=30
            else if (i_ice_parametrization(i) == ip_ice_adt_10) then
              n_parameter=36
            else if (i_ice_parametrization(i) == ip_ice_fu_phf) then
              n_parameter=5*n_ice_phf_term(i)+9
            else if (i_ice_parametrization(i) == ip_ice_agg_de) then
              n_parameter=14
            endif
!
            do j=1, n_parameter
              do k=1, n_band
                spectrum_rd%ice_parameter_list(j, k, i)                 &
     &            =ice_parameter_list(j, k, i)
              enddo
            enddo
!
          else
            spectrum_rd%l_ice_type(i)=.false.
          endif
        enddo
      endif
!
!     Block 13:
      if (l_present(13)) then
        spectrum_rd%l_present(13)=.true.
        do i=1, n_absorb
          if (l_retain_absorb(i)) then
            spectrum_rd%l_doppler_present(compressed_index(i))          &
     &        =l_doppler_present(i)
            if (l_doppler_present(i))                                   &
     &        spectrum_rd%doppler_correction(compressed_index(i))       &
     &          =doppler_correction(i)
          endif
        enddo
      else
        do i=1, n_absorb_retain
          spectrum_rd%l_doppler_present(i)=.false.
        enddo
      endif
!
!
!     Block 14:
      if (l_present(14)) then
        spectrum_rd%l_present(14)=.true.
        do i=1, n_band
          spectrum_rd%n_band_exclude(i)=n_band_exclude(i)
          do j=1, n_band_exclude(i)
            spectrum_rd%index_exclude(j, i)=index_exclude(j, i)
          enddo
        enddo
      endif
!
!
!     BLOCK 15 (we rely on the work done for block 11)
      IF (L_PRESENT(15).and.l_use_aod) THEN
         spectrum_rd%l_present(15)=.TRUE.
         spectrum_rd%n_aod_wavel = N_AOD_WAVEL

         DO I=1, N_AEROSOL_RETAIN
            I_INITIAL=INDEX_AEROSOL_RETAIN(I)
            IF (L_AEROSOL_SPECIES(I_INITIAL)) THEN
               spectrum_rd%i_aod_type(I) = I_AOD_TYPE(I_INITIAL)
               IF (I_AEROSOL_PARAMETRIZATION(I_INITIAL)                 &
     &              ==  IP_AEROSOL_PARAM_DRY) THEN
                  DO K=1, N_AOD_WAVEL
                    spectrum_rd%aod_absorption(1, I, K) =               &
     &               AOD_ABSORPTION(1, I_INITIAL, K)
                    spectrum_rd%aod_scattering(1, I, K) =               &
     &               AOD_SCATTERING(1, I_INITIAL, K)
                  ENDDO
               ELSE IF (I_AEROSOL_PARAMETRIZATION(I_INITIAL)            &
     &                   ==  IP_AEROSOL_PARAM_MOIST) THEN
                  DO J=1, NHUMIDITY(I_INITIAL)
                    DO K=1, N_AOD_WAVEL
                      spectrum_rd%aod_absorption(J, I, K) =             &
     &                 AOD_ABSORPTION(J, I_INITIAL, K)
                      spectrum_rd%aod_scattering(J, I, K) =             &
     &                 AOD_SCATTERING(J, I_INITIAL, K)
                    ENDDO
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
      ENDIF
!
!
      return
      END SUBROUTINE r2_compress_spectrum
#endif
#endif
