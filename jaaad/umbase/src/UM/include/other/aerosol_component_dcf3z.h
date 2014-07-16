!     ------------------------------------------------------------------
!     Extinctions at 550 nm. These are set only for aerosols in the SRA
!     profiles. Other values are set to 0.
!
      DATA                                                              &
     &    EXT_550NM_COMPONENT(IP_WATER_SOLUBLE)/4.721E+06_real64/       &
     &  , EXT_550NM_COMPONENT(IP_DUST_LIKE)/1.629E+05_real64/           &
     &  , EXT_550NM_COMPONENT(IP_OCEANIC)/7.210E+05_real64/             &
     &  , EXT_550NM_COMPONENT(IP_SOOT)/9.276E+06_real64/                &
     &  , EXT_550NM_COMPONENT(IP_ASH)/5.970E+06_real64/                 &
     &  , EXT_550NM_COMPONENT(IP_SULPHURIC)/6.382E+06_real64/           &
     &  , EXT_550NM_COMPONENT(IP_AMMONIUM_SULPHATE)/0.0E+00_real64/     &
     &  , EXT_550NM_COMPONENT(IP_SAHARAN_DUST)/0.0E+00_real64/          &
     &  , EXT_550NM_COMPONENT(IP_ACCUM_SULPHATE)/0.0E+00_real64/        &
     &  , EXT_550NM_COMPONENT(IP_AITKEN_SULPHATE)/0.0E+00_real64/       &
     &  , EXT_550NM_COMPONENT(IP_FRESH_SOOT)/0.0E+00_real64/            &
     &  , EXT_550NM_COMPONENT(IP_AGED_SOOT)/0.0E+00_real64/             &
     &  , EXT_550NM_COMPONENT(IP_SODIUM_CHLORIDE)/0.0E+00_real64/       &
     &  , EXT_550NM_COMPONENT(IP_SEASALT_FILM)/0.0E+00_real64/          &
     &  , EXT_550NM_COMPONENT(IP_SEASALT_JET)/0.0E+00_real64/           &
     &  , EXT_550NM_COMPONENT(IP_DUST_DIV1)/0.0E+00_real64/             &
     &  , EXT_550NM_COMPONENT(IP_DUST_DIV2)/0.0E+00_real64/             &
     &  , EXT_550NM_COMPONENT(IP_DUST_DIV3)/0.0E+00_real64/             &
     &  , EXT_550NM_COMPONENT(IP_DUST_DIV4)/0.0E+00_real64/             &
     &  , EXT_550NM_COMPONENT(IP_DUST_DIV5)/0.0E+00_real64/             &
     &  , EXT_550NM_COMPONENT(IP_DUST_DIV6)/0.0E+00_real64/             &
     &  , EXT_550NM_COMPONENT(IP_BIOMASS_1)/0.0E+00_real64/             &
     &  , EXT_550NM_COMPONENT(IP_BIOMASS_2)/0.0E+00_real64/
!
!
!     Densities for the sra aerosols. These are not actual densities,
!     which are not known, but rather conventional values. Extreme
!     care is required when setting mass mixing ratios in the radiation
!     scheme to ensure that the values are consistent with the
!     spectral data.
!
      DATA                                                              &
     &    DENSITY_COMPONENT(IP_WATER_SOLUBLE)/1.00E+03_real64/          &
     &  , DENSITY_COMPONENT(IP_DUST_LIKE)/1.00E+03_real64/              &
     &  , DENSITY_COMPONENT(IP_OCEANIC)/1.00E+03_real64/                &
     &  , DENSITY_COMPONENT(IP_SOOT)/1.00E+03_real64/                   &
     &  , DENSITY_COMPONENT(IP_ASH)/1.00E+03_real64/                    &
     &  , DENSITY_COMPONENT(IP_SULPHURIC)/1.00E+03_real64/
!
!     Density for Saharan dust, an extension to the SRA set. Again,
!     this is a conventional value to be used consistently.
      DATA                                                              &
     &    DENSITY_COMPONENT(IP_SAHARAN_DUST)/1.000E+03_real64/
!
!     Densities for other aerosols: again care is required over the
!     question of consistency, but these values are fairly widely
!     accepted, except perhaps in the case of fresh and aged soot,
!     where a value of 10^3 kgm-3 is often used, but a precise value
!     seems to be lacking. (Note 20/11/01: A value of 1900 kgm-3 is
!     now used).
!
      DATA                                                              &
     &    DENSITY_COMPONENT(IP_AMMONIUM_SULPHATE)/1.769E+03_real64/     &
     &  , DENSITY_COMPONENT(IP_ACCUM_SULPHATE)/1.769E+03_real64/        &
     &  , DENSITY_COMPONENT(IP_AITKEN_SULPHATE)/1.769E+03_real64/       &
     &  , DENSITY_COMPONENT(IP_FRESH_SOOT)/1.900E+03_real64/            &
     &  , DENSITY_COMPONENT(IP_AGED_SOOT)/1.900E+03_real64/             &
     &  , DENSITY_COMPONENT(IP_SODIUM_CHLORIDE)/2.165E+03_real64/       &
     &  , DENSITY_COMPONENT(IP_SEASALT_FILM)/2.165E+03_real64/          &
     &  , DENSITY_COMPONENT(IP_SEASALT_JET)/2.165E+03_real64/           &
     &  , DENSITY_COMPONENT(IP_DUST_DIV1)/2.65E+03_real64/              &
     &  , DENSITY_COMPONENT(IP_DUST_DIV2)/2.65E+03_real64/              &
     &  , DENSITY_COMPONENT(IP_DUST_DIV3)/2.65E+03_real64/              &
     &  , DENSITY_COMPONENT(IP_DUST_DIV4)/2.65E+03_real64/              &
     &  , DENSITY_COMPONENT(IP_DUST_DIV5)/2.65E+03_real64/              &
     &  , DENSITY_COMPONENT(IP_DUST_DIV6)/2.65E+03_real64/              &
     &  , DENSITY_COMPONENT(IP_BIOMASS_1)/1.50E+03_real64/              &
     &  , DENSITY_COMPONENT(IP_BIOMASS_2)/1.50E+03_real64/
!
!     ------------------------------------------------------------------
