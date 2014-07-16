! CANCLSTA start
!
! Purpose : Cross-Reference List of Ancillary Fields
!           Atmosphere. Ancillary Reference numbers,
!           Stash Codes, Model Codes and Logical file Numbers
!
! User changes:
!   05/10 kdcorbin - added tracer fluxes
!
! -------------------------------------------------------------------
!   Column A  : Ancillary Reference Number
!          B  : Internal Model Number
!          C  : Stash Item code
!          D  : Logical file number
!          E  : Field Description
!
!   A  B    C   D  E
!   ----------------------------------------------------------------
!   1  1   30   9  Land Sea Mask
!   2  1   33  10  Orography
!   3  1   34  10  Orographic Variance
!   4  1   35  10  Orographic gradient XX
!   5  1   36  10  Orographic gradient XY
!   6  1   37  10  Orographic gradient YY
!   7  1   60   1  Ozone
!   8              Not used
!   9  1   23   2  Snow Depth
!  10  1   20   3  Deep Soil Temperature
!  11  1   40   4  Vol SMC at Wilting
!  12  1   41   4  Vol SMC at Critical Point
!  13              Not used
!  14  1   43   4  Vol SMC at Saturation
!  15  1   44   4  Saturated Soil Conductivity
!  16              Not used
!  17  1   46   4  Thermal Capacity
!  18  1   47   4  Thermal Conductivity
!  19  1   50   5  Vegetation Fraction
!  20  1   51   5  Root Depth
!  21  1   52   5  Snow Free Surface Albedo
!  22  1   53   5  Deep Snow Surface Albedo
!  23  1   54   5  Surface Resistance to Evaporation
!  24  1   55   5  Surface Capacity
!  25  1   56   5  Infiltration Factor
!  26  1   26   5  Surface Roughness (vegetation)
!  27  1   31   7  Sea-ice Fraction
!  28  1   24   6  Sea-surface Temperature
!  29  1   32   7  Sea-ice Thickness
!  30  1   28   8  Surface Currents : u-component
!  31  1   29   8  Surface Currents : v-component
!  32  1   93   9  Runoff coastal outflow point
!  33              Not used (slab model)
!  34              Not used
!  35  1   48   4  Saturated soil water suction
!  36  1    9   2  Soil moisture in layers
!  37              Not used (SLAB) 
!  38              Not used (SLAB) 
!  39  1   58  12  Sulphur dioxide emission
!  40  1   59  12  Dimethyl sulphur emission
!  41  1   88  13  Sulphate aerosol mass mixing ratio
!  42  1   87  13  Sulphuric acid aerosol mass mixing ratio
!  43  1   85  13  Soot aerosol mass mixing ratio
!  44  1   57  14  Multi-level murk source term emission
!  45  1   90  14  Multi-level murk concentration
!  46  1   17  10  Silhouette area of orography (orog. roughness)
!  47  1   18  10  Peak to trough height (for orog. roughness scheme)
!  48  1  301  15  User ancillary field 1
!  49  1  302  15  User ancillary field 2
!  50  1  303  15  User ancillary field 3
!  51  1  304  15  User ancillary field 4
!  52  1  305  15  User ancillary field 5
!  53  1  306  15  User ancillary field 6
!  54  1  307  15  User ancillary field 7
!  55  1  308  15  User ancillary field 8
!  56  1  309  15  User ancillary field 9
!  57  1  310  15  User ancillary field 10
!  58  1  311  15  User ancillary field 11
!  59  1  312  15  User ancillary field 12
!  60  1  313  15  User ancillary field 13
!  61  1  314  15  User ancillary field 14
!  62  1  315  15  User ancillary field 15
!  63  1  316  15  User ancillary field 16
!  64  1  317  15  User ancillary field 17
!  65  1  318  15  User ancillary field 18
!  66  1  319  15  User ancillary field 19
!  67  1  320  15  User ancillary field 20
!  68  1  127  12  NH3 (ammonia) aerosol emission
!  69  1  128  23  Surface fresh soot aerosol emission
!  70  1  129  23  Elevated fresh soot aerosol emission
!  71              Not used
!  72  1  121  17  Natural Sulphur dioxide emissions
!  73  1  122  18  OH concentrations
!  74  1  123  18  HO2 concentrations
!  75  1  124  18  H2O2 concentrations
!  76  1  125  18  Ozone (CHEM) concentrations
!  77  1  126  12  Sulphur dioxide high level emission
!  78  1  251  24  Surface CO2 emissions
!  79  1  207   4  Clapp-Hornberger parameter
!  80  1  208   5  Leaf area index of vegetated fraction
!  81  1  209   5  Canopy height of vegetated fraction
!  82  1  160  19  Aerosol data for radiative forcing of climate change
!  83  1  216  20  Initial fractions of surface types
!  84  1  217  21  Initial leaf area index of plant functional types
!  85  1  218  21  Initial canopy height of plant functional types
!  86  1  213  21  Initial gridbox mean canopy conductance
!  87  1  219  22  Fraction of vegetation subject to disturbance
!  88  1  220   4  Snow free albedo of bare soil
!  89  1  223   4  Soil carbon content
!  90  1  321  16  User ancillary multi 1
!  91  1  322  16  User ancillary multi 2
!  92  1  323  16  User ancillary multi 3
!  93  1  324  16  User ancillary multi 4
!  94  1  325  16  User ancillary multi 5
!  95  1  326  16  User ancillary multi 6
!  96  1  327  16  User ancillary multi 7
!  97  1  328  16  User ancillary multi 8
!  98  1  329  16  User ancillary multi 9
!  99  1  330  16  User ancillary multi 10
! 100  1  331  16  User ancillary multi 11
! 101  1  332  16  User ancillary multi 12
! 102  1  333  16  User ancillary multi 13
! 103  1  334  16  User ancillary multi 14
! 104  1  335  16  User ancillary multi 15
! 105  1  336  16  User ancillary multi 16
! 106  1  337  16  User ancillary multi 17
! 107  1  338  16  User ancillary multi 18
! 108  1  339  16  User ancillary multi 19
! 109  1  340  16  User ancillary multi 20
! 110  1  341  25  Tropopause-based Ozone
! 111  1  505  26  Land Fraction
! 112  1  418  27  Dust parent soil clay fraction
! 113  1  419  27  Dust parent soil silt fraction
! 114  1  420  27  Dust parent soil sand fraction
! 115  1  421  27  Dust soil mass fraction div 1
! 116  1  422  27  Dust soil mass fraction div 2
! 117  1  423  27  Dust soil mass fraction div 3
! 118  1  424  27  Dust soil mass fraction div 4
! 119  1  425  27  Dust soil mass fraction div 5
! 120  1  426  27  Dust soil mass fraction div 6
! 121  1  130  28  Biomass surface emissions
! 122  1  131  28  Biomass elevated emissions
! 123  1  132  29  DMS concentration in seawater
! 124  1  153  30  River Water Storage
! 125  1  151  31  Riv Channel Sequence
! 126  1  152  31  Riv Channel Direction
!
! 127-154 not used in UM
!
! 155  1    5  10  Orographic X gradient
! 156  1    6  10  Orographic Y gradient 
! 157  1  351  38  Climatalogical biogenic aerosol mmr
! 158  1  352  39  Clim Biomass-burning (fresh) mmr
! 159  1  353  39  Clim Biomass-burning (aged) mmr
! 160  1  354  39  Clim Biomass-burning (in-cloud) mmr
! 161  1  355  40  Clim Black Carbon (fresh) mmr
! 162  1  356  40  Clim Black Carbon (aged) mmr
! 163  1  357  41  Clim Sea-salt (film mode) npm3
! 164  1  358  41  Clim Sea-salt (jet mode) npm3
! 165  1  359  42  Clim Sulphate (accumulation mode) mmr
! 166  1  360  42  Clim Sulphate (Aitken mode) mmr
! 167  1  361  42  Clim Sulphate (dissolved) mmr
! 168  1  362  43  Clim Dust size division 1 mmr
! 169  1  363  43  Clim Dust size division 2 mmr
! 170  1  364  43  Clim Dust size division 3 mmr
! 171  1  365  43  Clim Dust size division 4 mmr
! 172  1  366  43  Clim Dust size division 5 mmr
! 173  1  367  43  Clim Dust size division 6 mmr
! 174  1  368  44  Clim Organic Carbon from Fossil Fuels (fresh) mmr
! 175  1  369  44  Clim Organic Carbon from Fossil Fuels (aged) mmr
! 176  1  370  44  Clim Organic Carbon from Fossil Fuels (in-cloud) mmr
! 177  1  371  45  Clim Delta Aerosol mmr
! 178  1  480  46  Prognostic Ozone Tracer Cariolle Scheme
! 179  1  481  46  Cariolle Ozone Production - Loss (P-L)
! 180  1  482  46  Cariolle Ozone P-L wrt Ozone Mix Ratio
! 181  1  483  46  Cariolle Ozone Volume Mixing Ratio
! 182  1  484  46  Cariolle Ozone P-L wrt Temperature
! 183  1  485  46  Cariolle Ozone Clim Temp
! 184  1  486  46  Cariolle Ozone P-L wrt Ozone above a point
! 185  1  487  46  Cariolle Ozone Column above a point
! 186  1  134  47  Surface fresh fossil-fuel OC aerosol emissions
! 187  1  135  47  Elevated fresh fossil-fuel OC aerosol emissions
! 188  1  100  48  Flux of Tracer 1
! 189  1  100  49  Flux of Tracer 2
! 190  1  100  50  Flux of Tracer 3
! 191  1  100  51  Flux of Tracer 4
! 192  1  100  52  Flux of Tracer 5
! 193  1  100  53  Flux of Tracer 6
! 194  1  100  54  Flux of Tracer 7
! 195  1  100  55  Flux of Tracer 8
! 196  1  100  56  Flux of Tracer 9
! 197  1  100  57  Flux of Tracer 10
! 198  1  100  58  Flux of Tracer 11
! 199  1  100  59  Flux of Tracer 12
! 200  1  100  60  Flux of Tracer 13
! 201  1  100  61  Flux of Tracer 14
! 202  1  100  62  Flux of Tracer 15
! 203  1  100  63  Flux of Tracer 16
! 204  1  100  64  Flux of Tracer 17
! 205  1  100  65  Flux of Tracer 18
! 206  1  100  66  Flux of Tracer 19
! 207  1  100  67  Flux of Tracer 20
!  ------------------------------------------------------------------
! CANCLSTA end
