!     Module setting physical constants.
!
      REAL  (real64) ::                                                 &
     &    GRAV_ACC                                                      &
!           Acceleration due to gravity
     &  , R_GAS                                                         &
!           Universal gas constant
     &  , R_GAS_DRY                                                     &
!           Gas constant for dry air
     &  , CP_AIR_DRY                                                    &
!           Specific heat of dry air
     &  , RATIO_MOLAR_WEIGHT
!           Molecular weight of dry air/ molecular weight of water
!
      PARAMETER(                                                        &
     &    GRAV_ACC=9.80866E+00_real64                                   &
     &  , R_GAS=8.3143E+00_real64                                       &
     &  , R_GAS_DRY=287.026E+00_real64                                  &
     &  , CP_AIR_DRY=1.005E+03_real64                                   &
     &  , RATIO_MOLAR_WEIGHT=28.966E+00_real64/18.0153E+00_real64       &
     &  )
!
!     ------------------------------------------------------------------
