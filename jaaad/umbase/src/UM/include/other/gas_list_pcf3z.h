!     ------------------------------------------------------------------
!     Module to set indexing numbers of gaseous absorbing species.
!     The numbering 1-12 corresponds to LOWTRAN 7.
!
      INTEGER                                                           &
     &    NPD_GASES
!           Number of indexed gases
      PARAMETER(NPD_GASES=19)
!
      INTEGER                                                           &
     &    IP_H2O                                                        &
!           Index number of water vapour
     &  , IP_CO2                                                        &
!           Index number of carbon dioxide
     &  , IP_O3                                                         &
!           Index number of ozone
     &  , IP_N2O                                                        &
!           Index number of dinitrogen oxide
     &  , IP_CO                                                         &
!           Index number of carbon monoxide
     &  , IP_CH4                                                        &
!           Index number of methane
     &  , IP_O2                                                         &
!           Index number of oxygen
     &  , IP_NO                                                         &
!           Index number of nitrogen monoxide
     &  , IP_SO2                                                        &
!           Index number of sulphur dioxide
     &  , IP_NO2                                                        &
!           Index number of nitrogen dioxide
     &  , IP_NH3                                                        &
!           Index number of ammonia
     &  , IP_HNO3                                                       &
!           Index number of nitric acid
     &  , IP_N2                                                         &
!           Index number of nitrogen
     &  , IP_CFC11                                                      &
!           Index number of CFC11 (CFCl3)
     &  , IP_CFC12                                                      &
!           Index number of CFC12 (CF2Cl2)
     &  , IP_CFC113                                                     &
!           Index number of CFC113 (CF2ClCFCl2)
     &  , IP_HCFC22                                                     &
!           Index number of HCFC22 (CHF2Cl)
     &  , IP_HFC125                                                     &
!           Index number of HFC125 (C2HF5)
     &  , IP_HFC134A
!           Index number of HFC134A (CF3CFH2)
!
      PARAMETER(                                                        &
     &    IP_H2O=1                                                      &
     &  , IP_CO2=2                                                      &
     &  , IP_O3=3                                                       &
     &  , IP_N2O=4                                                      &
     &  , IP_CO=5                                                       &
     &  , IP_CH4=6                                                      &
     &  , IP_O2=7                                                       &
     &  , IP_NO=8                                                       &
     &  , IP_SO2=9                                                      &
     &  , IP_NO2=10                                                     &
     &  , IP_NH3=11                                                     &
     &  , IP_HNO3=12                                                    &
     &  , IP_N2=13                                                      &
     &  , IP_CFC11=14                                                   &
     &  , IP_CFC12=15                                                   &
     &  , IP_CFC113=16                                                  &
     &  , IP_HCFC22=17                                                  &
     &  , IP_HFC125=18                                                  &
     &  , IP_HFC134A=19                                                 &
     &  )
!
!     ------------------------------------------------------------------
