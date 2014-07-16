#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read a longwave spectral namelist.
!
! Purpose:
!   To read a longwave namelist into a spectral array.
!
! Method:
!   The spectrum is read into the dynamically allocated array
!   and then reduced to a more manageable size.
!
! Current Owner of Code: James Manners
!
!- ---------------------------------------------------------------------
      subroutine r2_lw_specin(ierr, cmessage                            &
     &  , lw_spectral_file                                              &
     &  , l_ch4, l_n2o, l_cfc11, l_cfc12                                &
     &  , l_cfc113, l_hcfc22, l_hfc125, l_hfc134a                       &
     &  , l_climat_aerosol                                              &
     &  , l_use_dust, l_use_arcldust                                    &
     &  , l_use_sulpc_direct, l_use_arclsulp                            &
     &  , l_use_soot_direct, l_use_arclblck                             &
     &  , l_use_bmass_direct, l_use_arclbiom                            &
     &  , l_use_seasalt_direct, l_use_arclsslt                          &
     &  , l_use_ocff_direct, l_use_arclocff                             &
     &  , l_use_biogenic, l_use_arcldlta                                &
     &  , l_murk_rad, l_use_aod                                         &
     &  , l_gas, l_continuum, l_drop, l_aerosol, l_ice                  &
     &  , lw_spectrum                                                   &
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
#include "error3a.h"
#include "stdio3a.h"
#include "c_mdi.h"
!
!     Dummy arguments
      CHARACTER  (LEN=80), Intent(IN) :: lw_spectral_file
!           Name of file containing the spectral data
      LOGICAL, Intent(IN) :: l_ch4
!           Absorption by methane is included
      LOGICAL, Intent(IN) :: l_n2o
!           Absorption by nitrous oxide is included
      LOGICAL, Intent(IN) :: l_cfc11
!           Absorption by cfc11 is included
      LOGICAL, Intent(IN) :: l_cfc12
!           Absorption by cfc12 is included
      LOGICAL, Intent(IN) :: l_cfc113
!           Absorption by cfc113 is included
      LOGICAL, Intent(IN) :: l_hcfc22
!           Absorption by hcfc22 is included
      LOGICAL, Intent(IN) :: l_hfc125
!           Absorption by hfc125 is included
      LOGICAL, Intent(IN) :: l_hfc134a
!           Absorption by hfc134a is included
      LOGICAL, Intent(IN) :: l_climat_aerosol
!           Climatological aerosols are to be included
      LOGICAL, Intent(IN) :: l_use_sulpc_direct
!           The direct effects of sulphate aerosols are
!           to be included
      LOGICAL, Intent(IN) :: l_use_arclsulp
!           The direct effects of sulphate aerosols from NWP climatology
!           are to be included
      LOGICAL, Intent(IN) :: l_use_dust
!           Use the direct radiative effects of dust in the LW
      LOGICAL, Intent(IN) :: l_use_arcldust
!           Use the direct radiative effects of mineral dust from
!           NWP climatology in the LW
      LOGICAL, Intent(IN) :: l_use_bmass_direct
!           Use the direct radiative effects of biomass in the LW
      LOGICAL, Intent(IN) :: l_use_arclbiom
!           Use the direct radiative effects of biomass from
!           NWP climatology in the LW
      LOGICAL, Intent(IN) :: l_use_soot_direct
!           Use the direct radiative effects of soot in the LW
      LOGICAL, Intent(IN) :: l_use_arclblck
!           Use the direct radiative effects of black-carbon from NWP
!           climatology in the LW
      LOGICAL, Intent(IN) :: l_use_seasalt_direct
!           Include the direct radiative effect of seasalt in the LW
      LOGICAL, Intent(IN) :: l_use_arclsslt
!           Include the direct radiative effect of sea-salt from NWP
!           climatology in the LW
      LOGICAL, Intent(IN) :: l_use_biogenic
!           Use the biogenic aerosol direct effect in the LW
      LOGICAL, Intent(IN) :: l_use_ocff_direct
!           Include the direct radiative effect of fossil-fuel OC in LW
      LOGICAL, Intent(IN) :: l_use_arclocff
!           Include the direct radiative effect of fossil-fuel organic
!           carbon aerosol from NWP climatology in the LW
      LOGICAL, Intent(IN) :: l_use_arcldlta
!           Include the direct radiative effect of delta aerosol
!           from NWP climatology in the LW
      LOGICAL, Intent(IN) :: l_murk_rad
!           Mesoscale aerosols are to be included
      LOGICAL, Intent(IN) :: l_use_aod
!           At least one of the aerosol optical depth diagnostics
!           is requested
!
!     Flags for physical processes
      LOGICAL, Intent(INOUT) :: l_gas
!           Flag to include gaseous absorption
      LOGICAL, Intent(INOUT) :: l_continuum
!           Flag to include continuum absorption
      LOGICAL, Intent(INOUT) :: l_drop
!           Flag to include scattering by droplets
      LOGICAL, Intent(INOUT) :: l_aerosol
!           Flag to include scattering by aerosols
      LOGICAL, Intent(INOUT) :: l_ice
!           Flag to include scattering by ice crystals
!
      INTEGER, Intent(OUT) :: ierr
!           Error flag
      CHARACTER(LEN=1024), Intent(OUT)   :: cmessage
!
!     Define the reduced LW spectrum.
      TYPE (spectrum) :: lw_spectrum
!
!
!     Local variables.
!
!     Declare the initial spectrum for the namelist.
!
#include "mxsize3a.h"
#include "spdec3a.h"
#include "lwsp3a.h"
!
!
!     Radiative variables for reducing the spectrum
!
#include "aerprm3a.h"
#include "aercmp3a.h"
#include "gasid3a.h"
!
      integer                                                           &
     &    ios
!           Status of I/O
!
      Logical                                                           &
     &    l_retain_absorb(npd_species)                                  &
!           Flag set to .TRUE. if the absorber is to be retained
     &  , l_gas_included(npd_gases)
!           Logical to test for actual gases included
      integer                                                           &
     &    n_absorb_retain                                               &
!           Number of absorbers to retain
     &  , index_absorb_retain(npd_species)                              &
!           Indices of absorbers to be retained
     &  , n_aerosol_retain                                              &
!           Number of aerosols in the spectral file to be retained
!           For the radiative calculation
     &  , index_aerosol_retain(npd_aerosol_species)                     &
!           Indexing numbers of the retained aerosols
     &  , compressed_index(npd_species)                                 &
!           Mapping from old to new indices of absorbers
     &  , n_aerosol_found
!           Number of aerosols for the current group of processes
!           Found in the spectral file
!
!     Other local variables.
!
      integer                                                           &
     &    i, j, k
!           Loop variable
!
      character                                                         &
     &     ch_ios*5
!           Character string for iostat error
!
!     Subroutines called
      external                                                          &
     &     r2_allocate_spectrum, r2_compress_spectrum
!
!
!     Each block is initialized as missing:
      l_present=.false.
!
!
! Some spectral files specify the aerosol_asymmetry and others specify
! the aerosol phase-function. We want to read in both types of spectral
! file and hence we need to find out whether the namelists contains
! one or the other. In order to do this we initialise aerosol_asymmetry
! and aerosol_phase_fnc to the missing data flag and check which one
! has been allocate proper values. If aerosol_phase_fnc is specified
! nothing else needs to be done, otherwise we copy aerosol_asymmetry
! into aerosol_phase_fnc.

      aerosol_asymmetry = rmdi
      nhumidity=1


!     Read the longwave spectrum as a namelist.
      open(unit=80, file=lw_spectral_file, iostat=ios)
      if (ios /= 0) then
        ierr=i_err_io
        write(ch_ios, '(i5)') ios
        cmessage='error opening longwave spectral file.'                &
     &    //' iostat='//ch_ios
        return
      endif
      read(80, r2lwsp)
      close(80)
!
!
!     Test for minimal requisite information.
      if ( .not.(l_present(0).and.                                      &
     &           l_present(6) ) ) then
        cmessage='longwave spectrum is deficient.'
        ierr=i_err_fatal
        return
      endif
!
      if (l_gas) then
        if (.not.l_present(5)) then
          write(iu_err, '(/a)')                                         &
     &      '*** warning: the lw spectrum contains no '                 &
     &      //'gaseous absorption data.'
          l_gas=.false.
        endif
      endif
!
      if (l_continuum) then
        if (.not.l_present(9)) then
          write(iu_err, '(/a)')                                         &
     &      '*** warning: the lw spectrum contains no '                 &
     &      //'continuum absorption data.'
          l_continuum=.false.
        endif
      endif
!
      if (l_drop) then
        if (.not.l_present(10)) then
          write(iu_err, '(/a)')                                         &
     &      '*** warning: the lw spectrum contains no '                 &
     &      //'data for water droplets.'
          l_drop=.false.
        endif
      endif
!
      if (l_aerosol) then
        if (.not.l_present(11)) then
          write(iu_err, '(/a)')                                         &
     &      '*** warning: the lw spectrum contains no '                 &
     &      //'data for aerosols.'
          l_aerosol=.false.
        endif
      endif
!
      if (l_ice) then
        if (.not.l_present(12)) then
          write(iu_err, '(/a)')                                         &
     &      '*** warning: the lw spectrum contains no '                 &
     &      //'data for ice crystals.'
          l_ice=.false.
        endif
      endif
!
!     Set reduced dimensions, either from the sizes of the fixed arrays
!     or from the arrays read in.
!
      lw_spectrum%npd_type=npd_type
      lw_spectrum%npd_band=max(n_band, 1)
      lw_spectrum%npd_species=max(n_absorb, 1)
      lw_spectrum%npd_scale_variable=npd_scale_variable
      lw_spectrum%npd_surface=npd_surface
      lw_spectrum%npd_continuum=npd_continuum
      lw_spectrum%npd_thermal_coeff=n_deg_fit+1
      lw_spectrum%npd_cloud_parameter=npd_cloud_parameter
      lw_spectrum%npd_aod_wavel=1
      lw_spectrum%npd_phase_term=1
!
!
!     Search the spectrum to find maximum dimensions.
!
      lw_spectrum%npd_exclude=1
      if (l_present(14)) then
        do i=1, n_band
          lw_spectrum%npd_exclude                                       &
     &      =max(lw_spectrum%npd_exclude, n_band_exclude(i))
        enddo
      endif
!
!     Search the spectrum to find those gases to be retained.
!     Water vapour, carbon dioxide and ozone are included
!     if present, but a warning is printed if they are
!     not included.
      do i=1, npd_gases
        l_gas_included(i)=.false.
      enddo
      n_absorb_retain=0
!
      do i=1, n_absorb
!
        l_retain_absorb(i)=.false.
        compressed_index(i)=0
!
        if ( (type_absorb(i) == ip_h2o).or.                             &
     &       (type_absorb(i) == ip_co2).or.                             &
     &       (type_absorb(i) == ip_o3).or.                              &
     &       ( (type_absorb(i) == ip_ch4).and.l_ch4 ).or.               &
     &       ( (type_absorb(i) == ip_n2o).and.l_n2o ).or.               &
     &       ( (type_absorb(i) == ip_cfc11).and.l_cfc11 ).or.           &
     &       ( (type_absorb(i) == ip_cfc12).and.l_cfc12 ).or.           &
     &       ( (type_absorb(i) == ip_cfc113).and.l_cfc113 ).or.         &
     &       ( (type_absorb(i) == ip_hcfc22).and.l_hcfc22 ).or.         &
     &       ( (type_absorb(i) == ip_hfc125).and.l_hfc125 ).or.         &
     &       ( (type_absorb(i) == ip_hfc134a).and.l_hfc134a ) ) then
          n_absorb_retain=n_absorb_retain+1
          index_absorb_retain(n_absorb_retain)=i
          compressed_index(i)=n_absorb_retain
          l_retain_absorb(i)=.true.
          l_gas_included(type_absorb(i))=.true.
        endif
!
      enddo
!
!
!     Print warning messages if those gases normally expected
!     are not present. This is done only for the call advancing
!     the integration. Subsequent calls are diagnostic and may
!     be made only for a limited spectral range when it is
!     not appropriate to include all gases.
!
        if (.not.l_gas_included(ip_h2o)) then
          write(iu_err, '(/a, /a)')                                     &
     &      '*** warning: water vapour is not included in the '         &
     &      , 'longwave spectral file.'
        endif
!
        if (.not.l_gas_included(ip_co2)) then
          write(iu_err, '(/a, /a)')                                     &
     &      '*** warning: carbon dioxide is not included in the '       &
     &      , 'longwave spectral file.'
        endif
!
        if (.not.l_gas_included(ip_o3)) then
          write(iu_err, '(/a, /a)')                                     &
     &      '*** warning: ozone is not included in the '                &
     &      , 'longwave spectral file.'
        endif
!
        if ((.not.l_gas_included(ip_ch4)).and.l_ch4) then
          write(iu_err, '(/a, /a)')                                     &
     &      '*** error: methane is not included in the longwave '       &
     &      , 'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          return
        endif
!
        if ((.not.l_gas_included(ip_n2o)).and.l_n2o) then
          write(iu_err, '(/a, /a)')                                     &
     &      '*** error: nitrous oxide is not included in the longwave ' &
     &      , 'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          return
        endif
!
        if ((.not.l_gas_included(ip_cfc11)).and.l_cfc11) then
          write(iu_err, '(/a, /a)')                                     &
     &      '*** error: cfc11 is not included in the longwave '         &
     &      , 'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          return
        endif
!
        if ((.not.l_gas_included(ip_cfc12)).and.l_cfc12) then
          write(iu_err, '(/a, /a)')                                     &
     &      '*** error: cfc12 is not included in the longwave '         &
     &      , 'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          return
        endif
!
        if ((.not.l_gas_included(ip_cfc113)).and.l_cfc113) then
          write(iu_err, '(/a, /a)')                                     &
     &      '*** error: cfc113 is not included in the longwave '        &
     &      , 'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          return
        endif
!
        if ((.not.l_gas_included(ip_hcfc22)).and.l_hcfc22) then
          write(iu_err, '(/a, /a)')                                     &
     &      '*** error: hcfc22 is not included in the longwave '        &
     &      , 'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          return
        endif
!
        if ((.not.l_gas_included(ip_hfc125)).and.l_hfc125) then
          write(iu_err, '(/a, /a)')                                     &
     &      '*** error: hfc125 is not included in the longwave '        &
     &      , 'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          return
        endif
!
        if ((.not.l_gas_included(ip_hfc134a)).and.l_hfc134a) then
          write(iu_err, '(/a, /a)')                                     &
     &      '*** error: hfc134a is not included in the longwave '       &
     &      , 'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          return
        endif
!
!
!     Set an appropriate reduced dimension.
      lw_spectrum%npd_species=max(n_absorb_retain, 1)
!
      lw_spectrum%npd_esft_term=1
      if (l_present(5)) then
        do i=1, n_band
          do j=1, n_band_absorb(i)
            if (l_retain_absorb(index_absorb(j, i)))                    &
     &        lw_spectrum%npd_esft_term                                 &
     &          =max(lw_spectrum%npd_esft_term                          &
     &          , i_band_esft(i, index_absorb(j, i)))
          enddo
        enddo
      endif
!
      lw_spectrum%npd_drop_type=1
      if (l_present(10)) then
        do i=1, npd_drop_type
          if (l_drop_type(i)) then
            lw_spectrum%npd_drop_type                                   &
     &        =max(lw_spectrum%npd_drop_type, i)
          endif
        enddo
      endif
!
      lw_spectrum%npd_ice_type=1
      if (l_present(12)) then
        do i=1, npd_ice_type
          if (l_ice_type(i)) then
            lw_spectrum%npd_ice_type                                    &
     &        =max(lw_spectrum%npd_ice_type, i)
          endif
        enddo
      endif
!
!
!
!     Aerosols must be treated carefully to allow for various
!     different combinations without requiring the spectral file
!     to be too constrained. Only those required will be retained.
!
!     Basic initialization to safe values.
      lw_spectrum%npd_humidities=1
      n_aerosol_retain=0
!
!     Check the spectral file for climatological aerosols
      if (l_climat_aerosol .OR. l_murk_rad) then
!
        if (l_present(11)) then
!
!         Search for the aerosols required for this scheme.
          n_aerosol_found=0
          do i=1, n_aerosol
!
            if ( (type_aerosol(i) == ip_water_soluble).or.              &
     &           (type_aerosol(i) == ip_dust_like).or.                  &
     &           (type_aerosol(i) == ip_oceanic).or.                    &
     &           (type_aerosol(i) == ip_soot).or.                       &
     &           (type_aerosol(i) == ip_sulphuric) ) then
              n_aerosol_retain=n_aerosol_retain+1
              index_aerosol_retain(n_aerosol_retain)=i
              n_aerosol_found=n_aerosol_found+1
            endif
          enddo
!
          if (n_aerosol_found /= 5) then
!
            ierr=i_err_fatal
            cmessage='the lw spectral file lacks some '                 &
     &        //'climatological aerosols.'
            return
!
          endif
        else
!
          ierr=i_err_fatal
          cmessage='lw spectral file contains no aerosol data.'
          return
!
        endif
!
      endif
!
!
      IF (L_USE_DUST .OR. L_USE_ARCLDUST) THEN
!
         IF (L_PRESENT(11)) THEN ! aerosol block present in spec file
!
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( (TYPE_AEROSOL(I) == IP_DUST_1).OR.                  &
     &              (TYPE_AEROSOL(I) == IP_DUST_2).OR.                  &
     &              (TYPE_AEROSOL(I) == IP_DUST_3).OR.                  &
     &              (TYPE_AEROSOL(I) == IP_DUST_4).OR.                  &
     &              (TYPE_AEROSOL(I) == IP_DUST_5).OR.                  &
     &              (TYPE_AEROSOL(I) == IP_DUST_6) ) THEN
                  N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                  INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                  N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF
!
            ENDDO
!
            IF (N_AEROSOL_FOUND /= 6) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The LW Spectral file lacks some '              &
     &            //'mineral dust aerosol.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='LW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!     Check the spectral file for soot aerosols.
      if (l_use_soot_direct .or. l_use_arclblck) then
        if (l_present(11)) then
!         Search for the aerosols required for this scheme.
          n_aerosol_found=0
          do i=1, n_aerosol
            if ( (type_aerosol(i) == ip_fresh_soot).or.                 &
     &           (type_aerosol(i) == ip_aged_soot)) then
              n_aerosol_retain=n_aerosol_retain+1
              index_aerosol_retain(n_aerosol_retain)=i
              n_aerosol_found=n_aerosol_found+1
            endif
          enddo
!
          if (n_aerosol_found /= 2) then
!
            ierr=i_err_fatal
            cmessage='the lw spectral file lacks some '                 &
     &        //'soot aerosol data.'
            return
!
          endif
!
        else
!
!
          ierr=i_err_fatal
          cmessage='lw spectral file contains no soot data.'
          return
!
        endif
!
      endif
!
!     Check the spectral file for biomass aerosol modes.
!     (Only required for the direct effect).
!
      IF (L_USE_BMASS_DIRECT .OR. L_USE_ARCLBIOM) THEN
!
         IF (L_PRESENT(11)) THEN ! aerosol block present in spec file
!
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( (TYPE_AEROSOL(I) == IP_BIOMASS_1).OR.               &
     &              (TYPE_AEROSOL(I) == IP_BIOMASS_2) ) THEN
                  N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                  INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                  N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF
!
            ENDDO
!
            IF (N_AEROSOL_FOUND /= 2) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The LW Spectral file lacks some '              &
     &            //'biomass aerosol.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='LW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!     Check the spectral file for sulphate aerosols. (These are
!     required only for the direct effect).
!
      if (l_use_sulpc_direct .or. l_use_arclsulp) then
!
        if (l_present(11)) then
!
!         Search for the aerosols required for this scheme.
          n_aerosol_found=0
          do i=1, n_aerosol
!
            if ( (type_aerosol(i) == ip_accum_sulphate).or.             &
     &           (type_aerosol(i) == ip_aitken_sulphate) ) then
              n_aerosol_retain=n_aerosol_retain+1
              index_aerosol_retain(n_aerosol_retain)=i
              n_aerosol_found=n_aerosol_found+1
            endif
          enddo
!
          if (n_aerosol_found /= 2) then
!
            ierr=i_err_fatal
            cmessage='the lw spectral file lacks some '                 &
     &        //'sulphate aerosols.'
            return
!
          endif
        else
!
          ierr=i_err_fatal
          cmessage='lw spectral file contains no aerosol data.'
          return
!
        endif
!
      endif
!
!
!     Check the spectral file for sea-salt aerosols. (These are
!     required only for the direct effect).
!
      if (l_use_seasalt_direct .or. l_use_arclsslt) then
!
        if (l_present(11)) then
!
!         Search for the aerosols required for this scheme.
          n_aerosol_found=0
          do i=1, n_aerosol
!
            if ( (type_aerosol(i) == ip_seasalt_film).or.               &
     &           (type_aerosol(i) == ip_seasalt_jet) ) then
              n_aerosol_retain=n_aerosol_retain+1
              index_aerosol_retain(n_aerosol_retain)=i
              n_aerosol_found=n_aerosol_found+1
            endif
          enddo
!
          if (n_aerosol_found /= 2) then
!
            ierr=i_err_fatal
            cmessage='the lw spectral file lacks some '                 &
     &        //'seasalt aerosols.'
            return
!
          endif
        else
!
          ierr=i_err_fatal
          cmessage='lw spectral file contains no aerosol data.'
          return
!
        endif
!
      endif
!
!     Check the spectral file for biogenic aerosol. (direct effect
!     only).
!
      if (l_use_biogenic) then
!
         if (l_present(11)) then ! aerosol block present in spec file
!
!           Search for the aerosol required for this scheme
            n_aerosol_found=0
            do i=1, n_aerosol
!
               if ( type_aerosol(i) == ip_biogenic ) then
                 n_aerosol_retain=n_aerosol_retain+1
                 index_aerosol_retain(n_aerosol_retain)=i
                 n_aerosol_found=n_aerosol_found+1
               endif
!
            enddo
!
            if (n_aerosol_found /= 1) then
!
               ierr=i_err_fatal
               cmessage='The LW Spectral file lacks some '              &
     &            //'biogenic aerosol.'
               return
!
            endif

          else
!
            ierr=i_err_fatal
            cmessage='LW Spectral file contains no aerosol data.'
            return
!
         endif
!
      endif
!
!     Check the spectral file for fossil-fuel organic carbon aerosol.
!     (Only required for the direct effect).
!
      IF (L_USE_OCFF_DIRECT .OR. L_USE_ARCLOCFF) THEN
!
         IF (L_PRESENT(11)) THEN ! aerosol block present in spec file
!
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( (TYPE_AEROSOL(I) == IP_OCFF_FRESH).OR.              &
     &              (TYPE_AEROSOL(I) == IP_OCFF_AGED) ) THEN
                  N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                  INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                  N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF
!
            ENDDO
!
            IF (N_AEROSOL_FOUND /= 2) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The LW Spectral file lacks some '              &
     &            //'fossil fuel org.carb. aerosol.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='LW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!     Check the spectral file for delta aerosol. (direct effect
!     only).
!
      if (l_use_arcldlta) then
!
         if (l_present(11)) then ! aerosol block present in spec file
!
!           Search for the aerosol required for this scheme
            n_aerosol_found=0
            do i=1, n_aerosol
!
               if ( type_aerosol(i) == ip_delta ) then
                 n_aerosol_retain=n_aerosol_retain+1
                 index_aerosol_retain(n_aerosol_retain)=i
                 n_aerosol_found=n_aerosol_found+1
               endif
!
            enddo
!
            if (n_aerosol_found /= 1) then
!
               ierr=i_err_fatal
               cmessage='The LW Spectral file lacks some '              &
     &            //'delta aerosol.'
               return
!
            endif

          else
!
            ierr=i_err_fatal
            cmessage='LW Spectral file contains no aerosol data.'
            return
!
         endif
!
      endif

!
!     Set an appropriate reduced dimension.
      lw_spectrum%npd_aerosol_species=max(n_aerosol_retain, 1)
!
!     Set the allowed number of humidities from the number of
!     retained aerosols.
!
      if (l_present(11)) then
        do i=1, n_aerosol_retain
          if ( (i_aerosol_parametrization(index_aerosol_retain(i)) ==   &
     &          ip_aerosol_param_phf_dry).or.                           &
     &         (i_aerosol_parametrization(index_aerosol_retain(i)) ==   &
     &          ip_aerosol_param_phf_moist) ) then
            lw_spectrum%npd_phase_term=max(lw_spectrum%npd_phase_term   &
     &        , n_aerosol_phf_term(index_aerosol_retain(i)))
          endif
          if ( (i_aerosol_parametrization(index_aerosol_retain(i)) ==   &
     &          ip_aerosol_param_moist).or.                             &
     &         (i_aerosol_parametrization(index_aerosol_retain(i)) ==   &
     &          ip_aerosol_param_phf_moist) ) then
            lw_spectrum%npd_humidities=max(lw_spectrum%npd_humidities   &
     &        , nhumidity(index_aerosol_retain(i)))
          endif
        enddo

! Check whether aerosol_phase_fnc or aerosol_asymmetry is specified.

        if (ANY(aerosol_asymmetry /= rmdi)) then
           aerosol_phase_fnc(:,1,:,:)=aerosol_asymmetry
        endif

      endif
!
!
!     Check that block 15 is present if the aerosol optical depth
!     was requested.
      IF (L_USE_AOD) THEN
        IF(.NOT. L_PRESENT(15)) THEN
          IERR = I_ERR_FATAL
          CMESSAGE='Block 15 needed in the LW spectral file.'
          RETURN
        ENDIF
       ! Check that the number of wavelengths in the spectral
       ! file is not larger than the number set in the include
       ! file MXSIZE3A.
        IF(N_AOD_WAVEL  >   NPD_AOD_WAVEL) THEN
          IERR = I_ERR_FATAL
          CMESSAGE='Increase NPD_AOD_WAVEL in MXSIZE3A.'
          RETURN
        ENDIF
        lw_spectrum%npd_aod_wavel = NPD_AOD_WAVEL
      ENDIF
!
!
!     Allocate space for the LW spectrum.
! DEPENDS ON: r2_allocate_spectrum
      call r2_allocate_spectrum(lw_spectrum)
!
!     Transfer the large namelist to the reduced spectrum.
!
!
! DEPENDS ON: r2_compress_spectrum
      call r2_compress_spectrum(                                        &
!                       Spectral array in namelist
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
     &  , lw_spectrum                                                   &
     &  )
!
!
!
      return
      END SUBROUTINE r2_lw_specin
#endif
