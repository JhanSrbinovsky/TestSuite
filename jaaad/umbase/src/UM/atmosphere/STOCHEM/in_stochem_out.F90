#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      MODULE IN_STOCHEM_OUT
! ----------------------------------------------------------------------
! Purpose:
! To set output and other control options for STOCHEM.
!
! Method:
!
! Original Programmer: Colin Johnson
!
! Current code owner: Colin Johnson
!
! History:
! Date        Version     Comment
! -------     -------     -------
! 24/5/01     1.0         Original              Colin Johnson
! 17/09/03    6.0         Correct non-standard F90 continuation lines.
!                                                              P.Dando
! 20/08/04    6.1         Changes to allow vectorisation on SX6.
!                                                   M. Sanderson.
!  6.2   15/11/05  Extra variables added. M.G. Sanderson
!  6.2   01/03/06  All chemical fluxes have explicit STASH codes.
!                                                 M.G. Sanderson
!
!VVV  Vn5.2.1  IN_STOCHEM_OUT 2/8/01 - Original version
! ----------------------------------------------------------------------

! 3-D Output Species: these can be specified in any order.
      INTEGER, PARAMETER :: numchem=28
      CHARACTER(LEN=8), DIMENSION(numchem), PARAMETER ::                &
     & outnames=(/                                                      &
     &  'OH      ','NO      ','NO2     ','CO      ','CH4     ',         &
     &  'HCHO    ','O3      ','H2      ','HNO3    ','H2O2    ',         &
     &  'HO2     ','C2H6    ','PAN     ','C2H4    ','C3H8    ',         &
     &  'SO2     ','DMS     ','SA      ','AMMSUL  ','NAER    ',         &
     &  'O3_STRAT','H2O     ','Be-7    ','Be-10   ','Rn-222  ',         &
     &  'Pb-210  ','STRAT   ','TROP    '                                &
     & /)

! Fixed length header size
      INTEGER, PARAMETER :: len_flhead=256

! Set length of flux name arrays
      INTEGER, PARAMETER :: len_flx_str=54

! Set up output reaction list array sizes.
      INTEGER, PARAMETER :: nphotnames    = 23 ! No. photolysis rates
      INTEGER, PARAMETER :: nfluxnamesa   = 85 ! No. chemical reactions
      INTEGER, PARAMETER :: nfluxnamesb   = 54 ! No. chemical reactions
      INTEGER, PARAMETER :: nemissnames   = 30 ! No. emissions
      INTEGER, PARAMETER :: nddepnames    = 24 ! No. dry deposition rate
      INTEGER, PARAMETER :: nwdepnames    = 21 ! No. wet deposition rate
      INTEGER, PARAMETER :: nozonenames   = 14 ! No. ozone prod/dest ter
      INTEGER, PARAMETER :: nadvfluxnames = 60 ! No. fluxes across bound
      INTEGER, PARAMETER :: nextraflux    = 25 ! No. miscellaneous terms

! Calculate total number of output fluxes
      INTEGER, PARAMETER :: numflux = nphotnames                        &
     &                              + nfluxnamesa                       &
     &                              + nfluxnamesb                       &
     &                              + nemissnames                       &
     &                              + nddepnames                        &
     &                              + nwdepnames                        &
     &                              + nadvfluxnames                     &
     &                              + nozonenames                       &
     &                              + nextraflux

! Select true if a restart dump for STOCHEM is needed in an NRUN,
! if false an internal initialistaion is used.
      LOGICAL            :: stochstartdump = .FALSE.

! Set these logicals to TRUE when trace gases from STOCHEM are to be
! coupled to UM radiation scheme
!     LOGICAL            :: l_coupled_ch4  = .FALSE.
!     LOGICAL            :: l_coupled_o3   = .FALSE.

! These are used if the restart dump is incompatible with the fluxes
! specified in this run (see below).
! done automatically now
!     INTEGER, PARAMETER :: numflux_r=233,num3dflux_r=11,numchem_r=17
      INTEGER, PARAMETER :: num3dflux_r=11 ! for inicon_r

! num3dflux, num3dflux_dim are now set in Initialise_Stochem.
! num3dflux_dim used for array size so that it can be made an odd number
! to prevent memory bank conflicts.
! bank conflicts.
      INTEGER :: num3dflux
      INTEGER :: num3dflux_dim
      INTEGER, PARAMETER :: max3dflux=100

! Output flux names, and 3-D output requests.
! Turn final digit to 1 for 3-D output
! The fluxnames are given in several arrays because of the continuation
! card limit (t3e=99)and recomposed into the FLUXNAMES array.
      CHARACTER(LEN=len_flx_str), DIMENSION(nphotnames), PARAMETER ::   &
     &  photnames=(/                                                    &
     &   'O3+hv -> O3P                             301  24151  0',      &
     &   'O3+hv -> O1D                             302  24152  0',      &
     &   'NO2+hv -> NO+O3P                         303  24153  0',      &
     &   'H2O2+hv -> 2OH                           304  24154  0',      &
     &   'HNO3+hv -> OH+NO2                        305  24155  0',      &
     &   'HCHO+hv -> CO+2HO2                       306  24156  0',      &
     &   'HCHO+hv -> H2+CO                         307  24157  0',      &
     &   'CH3CHO+hv -> CH3O2+HO2+CO                308  24158  0',      &
     &   'CH3COE+hv -> C2H5O2+CH3COO2              309  24159  0',      &
     &   'HO2NO2+hv -> HO2+NO2                     310  24160  0',      &
     &   'MGLYOX+hv -> CH3COO2+HO2+CO              311  24161  0',      &
     &   'GLYOX+hv  -> HCHO+CO                     312  24162  0',      &
     &   'NO3+hv -> NO+O2                          313  24163  0',      &
     &   'NO3+hv -> NO2+OP                         314  24164  0',      &
     &   'N2O5+hv -> NO2+NO3                       315  24165  0',      &
     &   'CH3OOH+hv -> HCHO+HO2+OH                 316  24166  0',      &
     &   'PAN+hv -> CH3COO2+NO2                    317  24167  0',      &
     &   'ACETONE+hv -> CH3COO2+CH3O2              318  24168  0',      &
     &   'C2H5OOH+hv -> CH3CHO+HO2+OH              319  24169  0',      &
     &   'C3H7OOH+hv -> ACETONE+HO2+OH             320  24170  0',      &
     &   'C4H9OOH+hv -> CH3COE+HO2+OH              321  24171  0',      &
     &   'ISOPOOH+hv -> MVK+HCHO+HO2               322  24172  0',      &
     &   'MVKOOH+hv -> MGLYOX+HCHO+HO2             323  24173  0'/)

      CHARACTER(LEN=len_flx_str), DIMENSION(nfluxnamesa), PARAMETER ::  &
     &  fluxnamesa=(/                                                   &
     &   'O3P+O2 -> O3                               1  24001  0',      &
     &   'O3P+NO -> NO2                              5  24005  0',      &
     &   'O1D+M -> O3P                               7  24007  0',      &
     &   'H2O+O1D -> 2OH                             8  24008  0',      &
     &   'NO+O3 -> NO2+O2                           11  24011  0',      &
     &   'NO2+O3 -> NO3+O2                          12  24012  0',      &
     &   'OH+O3 -> HO2+O2                           13  24013  0',      &
     &   'HO2+O3 -> OH+2O2                          14  24014  0',      &
     &   'NO+NO3 -> 2NO2                            15  24015  0',      &
     &   'NO2+O3P -> NO+O2                          16  24016  0',      &
     &   'HO2+NO -> NO2+OH                          17  24017  0',      &
     &   'NO2+NO3 -> NO+NO2                         19  24019  0',      &
     &   'NO3+NO2 -> N2O5                           20  24020  0',      &
     &   'NO2+OH -> HNO3                            21  24021  0',      &
     &   'NO2+HO2 -> HO2NO2                         22  24022  0',      &
     &   'HO2NO2+M -> HO2+NO2                       23  24023  0',      &
     &   'OH+HO2NO2 -> H2O+NO2                      24  24024  0',      &
     &   '2NO3 -> 2NO2+O2                           27  24027  0',      &
     &   'N2O5 -> NO2+NO3                           29  24029  0',      &
     &   'HO2+OH -> H2O+O2                          30  24030  0',      &
     &   'H2O2+OH -> HO2                            31  24031  0',      &
     &   'NO3+HO2 -> HNO3                           32  24032  0',      &
     &   'H2+OH -> HO2                              33  24033  0',      &
     &   'NO3+HO2 -> OH+NO2                         34  24034  0',      &
     &   'HNO3+OH -> NO3                            35  24035  0',      &
     &   '2HO2 -> H2O2+O2                           36  24036  0',      &
     &   'SO2+OH -> HO2+SA                          39  24039  0',      &
     &   'N2O5 -> 2NAER                             41  24041  0',      &
     &   'HNO3 -> NAER                              42  24042  0',      &
     &   'ORGNIT -> NAER                            43  24043  0',      &
     &   'CH3OOH+OH -> CH3O2+H2O                    44  24044  0',      &
     &   'CH3OOH+OH -> HCHO+OH+H2O                  45  24045  0',      &
     &   'C2H5OOH+OH -> C2H5O2+H2O                  46  24046  0',      &
     &   'C3H7OOH+OH -> C3H7O2+H2O                  47  24047  0',      &
     &   'C4H9OOH+OH -> SC4H9O2+H2O                 48  24048  0',      &
     &   'CH4+OH -> CH3O2+H2O                       59  24059  1',      &
     &   'CH3O2+NO -> HCHO+NO2+HO2                  60  24060  0',      &
     &   '2CH3O2 -> 2HCHO+2HO2                      61  24061  0',      &
     &   '2CH3O2 -> HCHO+CH3OH                      62  24062  0',      &
     &   'CH3OH+OH -> HCHO+HO2                      63  24063  0',      &
     &   'HO2+CH3O2 -> CH3OOH                       65  24065  0',      &
     &   'HCHO+OH -> HO2+CO                         66  24066  0',      &
     &   'HCHO+NO3 -> HNO3+HO2+CO                   67  24067  0',      &
     &   'CO+OH -> CO2+HO2                          70  24070  0',      &
     &   'C2H6+OH -> C2H5O2+H2O                     71  24071  0',      &
     &   'C2H5O2+NO -> CH3CHO+NO2+HO2               72  24072  0',      &
     &   'C2H5O2+CH3O2 -> CH3CHO+HCHO+HO2           73  24073  0',      &
     &   'CH3COO2+CH3O2 -> 2HCHO                    74  24074  0',      &
     &   'CH3CHO+OH -> CH3COO2                      75  24075  0',      &
     &   'CH3COO2+NO2 -> PAN                        77  24077  0',      &
     &   'PAN -> CH3COO2+NO2                        78  24078  0',      &
     &   'NO+CH3COO2 -> CH3O2+NO2                   79  24079  0',      &
     &   'CH3O2+CH3COO2 -> HCHO+CH3O2+HO2           80  24080  0',      &
     &   'nC4H10+OH -> SC4H9O2                      81  24081  0',      &
     &   'SC4H9O2+NO -> CH3COE+NO2+HO2              83  24083  0',      &
     &   'SC4H9O2+CH3O2 -> CH3COE+HCHO+HO2          84  24084  0',      &
     &   'CH3COE+OH -> CH3COX                       86  24086  0',      &
     &   '2C2H5O2 -> CH3CHO+HO2                     90  24090  0',      &
     &   '2CH3COO2 -> 2CH3O2                        91  24091  0',      &
     &   'OH+C3H8 -> C3H7O2                         92  24092  0',      &
     &   'C3H7O2+NO -> ACETONE+NO2+HO2              93  24093  0',      &
     &   'OH+ACETONE -> ACETO2                      94  24094  0',      &
     &   'ACETO2+NO -> CH3COO2+HCHO+NO2             95  24095  0',      &
     &   'ACETO2+CH3O2 -> CH3COO2+HCHO+HO2          96  24096  0',      &
     &   'CH3O2+C3H7O2 -> ACETONE+HCHO+HO2          97  24097  0',      &
     &   'OH+PAN -> HCHO+NO3                        98  24098  0',      &
     &   'HO2+C2H5O2 -> C2H5OOH                     99  24099  0',      &
     &   'C2H5OOH+OH -> CH3CHO+OH                  100  24100  0',      &
     &   'HO2+C3H7O2 -> C3H7OOH                    101  24101  0',      &
     &   'C3H7OOH+OH -> ACETONE+OH                 102  24102  0',      &
     &   'HO2+C4H9O2 -> C4H9OOH                    103  24103  0',      &
     &   'C4H9OOH+OH -> CH3COE+OH                  104  24104  0',      &
     &   'NO+CH3COX -> CH3COO2+CH3CHO+NO2          105  24105  0',      &
     &   'CH3O2+CH3COX -> CH3COO2+CH3CHO+HCHO+HO2  106  24106  0',      &
     &   'C2H4+OH -> CH2O2C                        109  24109  0',      &
     &   'NO+CH2O2C -> HCHO+NO2+HO2                110  24110  0',      &
     &   'CH3O2+CH2O2C -> HCHO+HO2                 111  24111  0',      &
     &   'O3+C2H4 -> 1.47HCHO+.13H2+.31CO+0.2HO2   112  24112  0',      &
     &   'O3+C3H6 -> HCHO+.3CH4+.4CO+.28OH+...     123  24123  0',      &
     &   'O3+C3H6 -> CH3CHO+.24H2+.58CO+.18HO2     124  24124  0',      &
     &   'C3H6+OH -> CH3CHX                        125  24125  0',      &
     &   'NO+CH3CHX -> CH3CHO+HCHO+NO2+HO2         126  24126  0',      &
     &   'CH3O2+CH3CHX -> CH3CHO+2HCHO+2HO2        127  24127  0',      &
     &   'O3+C5H8 -> MVK+.22HCHO+.78CO+.27HO2+...  128  24128  0',      &
     &   'O3+MVK -> MGLYOX+.24HCHO+.76CO+.36HO2+.  129  24129  0'/)

      CHARACTER(LEN=len_flx_str), DIMENSION(nfluxnamesb), PARAMETER ::  &
     &  fluxnamesb=(/                                                   &
     &   'C2H6+NO3 -> C2H5O2+HNO3                  201  24201  0',      &
     &   'nC4H10+NO3 -> SC4H9O2+HNO3               202  24202  0',      &
     &   'C2H4+NO3 -> RNC2H4+HO2                   203  24203  0',      &
     &   'OH+RNC2H4 -> HCHO+NO2                    204  24204  0',      &
     &   'C3H6+NO3 -> RNC3H6+HO2                   205  24205  0',      &
     &   'OH+RNC3H6 -> CH3CHO+NO2                  206  24206  0',      &
     &   'CH3CHO+NO3 -> CH3COO2+HNO3               208  24208  0',      &
     &   'C5H8+NO3 -> RNC5H8+HO2                   209  24209  0',      &
     &   'OH+RNC5H8 -> MVK+NO2                     210  24210  0',      &
     &   'HO2+TOLP1 -> MEMALD+GLYOX                212  24212  0',      &
     &   'NO2+TOLP1 -> ORGNIT                      213  24213  0',      &
     &   'OH+DMS -> CH3SO+HCHO                     214  24214  0',      &
     &   'OH+DMS -> DMSO+HO2                       215  24215  0',      &
     &   'CH3SO+O3 -> CH3SO2                       216  24216  0',      &
     &   'CH3SO+NO2 -> CH3SO2+NO                   217  24217  0',      &
     &   'CH3SO2+O3 -> CH3SO3                      218  24218  0',      &
     &   'CH3SO2+NO2 -> CH3SO3+NO                  219  24219  0',      &
     &   'CH3SO2+O2 -> CH3O2+SO2                   220  24220  0',      &
     &   'CH3SO3+HO2 -> MSA                        221  24221  0',      &
     &   'CH3SO3+O2 -> CH3O2+SA                    222  24222  0',      &
     &   'CH3SO3+HCHO -> MSA+CO+HO2                223  24223  0',      &
     &   'OH+DMSO -> DMSO2+HO2                     224  24224  0',      &
     &   'OH+DMSO2 -> DMSP                         225  24225  0',      &
     &   'DMSP+NO -> CH3SO2+HCHO+NO2               226  24226  0',      &
     &   'DMSP+CH3O2 -> CH3SO2+HCHO+HO2            227  24227  0',      &
     &   'NO3+DMS -> CH3SO+HCHO+HNO3               228  24228  0',      &
     &   'OXYL+OH -> MEMALD+MGLYOX+HO2             230  24230  0',      &
     &   'NO2+OXYL1 -> ORGNIT                      231  24231  0',      &
     &   'OH+MEMALD -> MEMALD1                     232  24232  0',      &
     &   'NO+MEMALD1 -> MGLYOX+GLYOX+NO2+HO2       233  24233  0',      &
     &   'TOLUENE+OH -> MEMALD+GLYOX               234  24234  0',      &
     &   'TOLUENE+OH -> TOLP1                      235  24235  0',      &
     &   'OXYL+OH -> OXYL1                         236  24236  0',      &
     &   'OH+ORGNIT -> MEMALD+GLYOX+NO2            237  24237  0',      &
     &   'NO3+ORGNIT -> MEMALD+GLYOX+2NO2          238  24238  0',      &
     &   'OXYL1+HO2 -> MGLYOX+MEMALD               240  24240  0',      &
     &   'MEMALD1+CH3O2 -> MGLYOX+GLYOX+HCHO+2HO2  241  24241  0',      &
     &   'RO2IP1+HO2 -> ISOPOOH                    244  24244  0',      &
     &   'ISOPOOH+OH -> MVK+HCHO+OH                245  24245  0',      &
     &   'RO2IP2+HO2 -> MVKOOH                     246  24246  0',      &
     &   'MVKOOH+OH -> MGLYOX+HCHO+OH              247  24247  0',      &
     &   'MGLYOX+OH -> CH3COO2+CO                  248  24248  0',      &
     &   'C5H8+OH -> RO2IP1                        251  24251  0',      &
     &   'NO+RO2IP1 -> MVK+HCHO+NO2+HO2            252  24252  0',      &
     &   'OH+MVK -> RO2IP2                         253  24253  0',      &
     &   'NO+RO2IP2 -> MGLYOX+HCHO+NO2+HO2         254  24254  0',      &
     &   'EXTRA+OH -> CH3O2+H2O                    255  24255  0',      &
     &   'HSO3(aq)+H2O2(aq) -> SO4(aq)             260  24260  0',      &
     &   'HSO3(aq)+O3(aq) -> SO4(aq)               261  24261  0',      &
     &   'SO3(aq)+O3(aq) -> SO4(aq)                262  24262  0',      &
     &   'SO2(aq) -> SO4(aq)                       263  24263  0',      &
     &   '2NH4(aq)+SO4(aq) -> (NH4)2SO4(aq)        264  24264  0',      &
     &   'Be-7 ->                                  290  24290  0',      &
     &   'Rn-222 ->                                291  24291  0'/)

      CHARACTER(LEN=len_flx_str), DIMENSION(nemissnames), PARAMETER ::  &
     &  emissnames=(/                                                   &
     &   'NO emission                              404  25104  1',      &
     &   'NO2 Stratospheric input                  405  25105  0',      &
     &   'CO emission                              408  25108  0',      &
     &   'CH4 emission                             409  25109  0',      &
     &   'HCHO emission                            410  25110  0',      &
     &   'O3 Stratospheric input                   411  25111  0',      &
     &   'H2 emission                              412  25112  0',      &
     &   'HNO3 Stratospheric input                 413  25113  0',      &
     &   'C2H6 emission                            417  25117  0',      &
     &   'CH3CHO emission                          419  25119  0',      &
     &   'nC4H10 emission                          423  25123  0',      &
     &   'SO2 emission                             426  25126  0',      &
     &   'C2H4 emission                            427  25127  0',      &
     &   'C3H6 emission                            428  25128  0',      &
     &   'oXYL emission                            429  25129  0',      &
     &   'C3H8 emission                            431  25131  0',      &
     &   'CH3OH emission                           436  25136  0',      &
     &   'ACETONE emission                         437  25137  0',      &
     &   'C5H8 emission                            448  25148  1',      &
     &   'TOLUENE emission                         454  25154  0',      &
     &   'DMS emission                             460  25160  0',      &
     &   'NH3 emission                             470  25170  0',      &
     &   'EXTRA emission                           452  25152  0',      &
     &   'Be-7 emission                            474  25174  0',      &
     &   'Be-10 emission                           475  25175  0',      &
     &   'Rn-222 emission                          476  25176  0',      &
     &   'STRAT emission                           478  25178  0',      &
     &   'TROP emission                            479  25179  0',      &
     &   'Lightning NOx emission                   480  25180  1',      &
     &   'Aircraft NOx emission                    481  25181  0'/)

      CHARACTER(LEN=len_flx_str), DIMENSION(nddepnames), PARAMETER ::   &
     &  ddepnames=(/                                                    &
     &   'NO2 dry dep                              505  25205  1',      &
     &   'CO dry dep                               508  25208  1',      &
     &   'CH4 dry dep                              509  25209  1',      &
     &   'O3 dry dep                               511  25211  1',      &
     &   'H2 dry dep                               512  25212  0',      &
     &   'HNO3 dry dep                             513  25213  1',      &
     &   'H2O2 dry dep                             514  25214  0',      &
     &   'PAN dry dep                              521  25221  0',      &
     &   'CH3OOH dry dep                           522  25222  0',      &
     &   'SO2 dry dep                              526  25226  1',      &
     &   'SA dry dep                               530  25230  1',      &
     &   'C3H7OOH dry dep                          533  25233  0',      &
     &   'C2H5OOH dry dep                          534  25234  0',      &
     &   'C4H9OOH dry dep                          535  25235  0',      &
     &   'NH42SO4 dry dep                          539  25239  1',      &
     &   'ISOPOOH dry dep                          552  25252  0',      &
     &   'MVKOOH dry dep                           553  25253  0',      &
     &   'NAER dry dep                             558  25258  0',      &
     &   'MSA dry dep                              568  25268  0',      &
     &   'ORGNIT dry dep                           569  25269  0',      &
     &   'NH3 dry dep                              570  25270  1',      &
     &   'Be-7 dry dep                             574  25274  0',      &
     &   'Be-10 dry dep                            575  25275  0',      &
     &   'Pb-210 dry dep                           577  25277  0'       &
     &/)

      CHARACTER(LEN=len_flx_str), DIMENSION(nwdepnames), PARAMETER ::   &
     &  wdepnames=(/                                                    &
     &   'N2O5 wet dep                             607  25307  0',      &
     &   'HCHO wet dep                             610  25310  0',      &
     &   'HNO3 wet dep                             613  25313  0',      &
     &   'H2O2 wet dep                             614  25314  0',      &
     &   'CH3OOH wet dep                           622  25322  0',      &
     &   'SO2 wet dep                              626  25326  0',      &
     &   'SA wet dep                               630  25330  0',      &
     &   'C3H7OOH wet dep                          633  25333  0',      &
     &   'C2H5OOH wet dep                          634  25334  0',      &
     &   'C4H9OOH wet dep                          635  25335  0',      &
     &   'CH3OH wet dep                            636  25336  0',      &
     &   'NH42SO4 wet dep                          639  25339  0',      &
     &   'ISOPOOH wet dep                          652  25352  0',      &
     &   'MVKOOH wet dep                           653  25353  0',      &
     &   'NAER wet dep                             658  25358  0',      &
     &   'MSA wet dep                              668  25368  0',      &
     &   'ORGNIT wet dep                           669  25369  0',      &
     &   'NH3 wet dep                              670  25370  0',      &
     &   'Be-7 wet dep                             674  25374  0',      &
     &   'Be-10 wet dep                            675  25375  0',      &
     &   'Pb-210 wet dep                           677  25377  0'       &
     &/)

      CHARACTER(LEN=len_flx_str),DIMENSION(nadvfluxnames), PARAMETER :: &
     &  advfluxnames=(/                                                 &
     &   'O3 across tropopause                     711  24301  1',      &
     &   'O3(s) across tropopause                  771  24302  0',      &
     &   'STRAT across tropopause                  778  24305  0',      &
     &   'TROP across tropopause                   779  24306  0',      &
     &   'Be-7 across tropopause                   774  24303  0',      &
     &   'Be-10 across tropopause                  775  24304  0',      &
     &   'O3 across 100 hPa                        811  24311  0',      &
     &   'O3(s) across 100 hPa                     871  24312  0',      &
     &   'STRAT across 100 hPa                     878  24315  0',      &
     &   'TROP across 100 hPa                      879  24316  0',      &
     &   'Be-7 across 100 hPa                      874  24313  0',      &
     &   'Be-10 across 100 hPa                     875  24314  0',      &
     &   'O3 across 200 hPa                        911  24321  1',      &
     &   'O3(s) across 200 hPa                     971  24322  1',      &
     &   'STRAT across 200 hPa                     978  24325  0',      &
     &   'TROP across 200 hPa                      979  24326  0',      &
     &   'Be-7 across 200 hPa                      974  24323  0',      &
     &   'Be-10 across 200 hPa                     975  24324  0',      &
     &   'O3 across 300 hPa                       1011  24331  0',      &
     &   'O3(s) across 300 hPa                    1071  24332  0',      &
     &   'STRAT across 300 hPa                    1078  24335  0',      &
     &   'TROP across 300 hPa                     1079  24336  0',      &
     &   'Be-7 across 300 hPa                     1074  24333  0',      &
     &   'Be-10 across 300 hPa                    1075  24334  0',      &
     &   'O3 across 400 hPa                       1111  24341  0',      &
     &   'O3(s) across 400 hPa                    1171  24342  0',      &
     &   'STRAT across 400 hPa                    1178  24345  0',      &
     &   'TROP across 400 hPa                     1179  24346  0',      &
     &   'Be-7 across 400 hPa                     1174  24343  0',      &
     &   'Be-10 across 400 hPa                    1175  24344  0',      &
     &   'O3 across 500 hPa                       1211  24351  0',      &
     &   'O3(s) across 500 hPa                    1271  24352  0',      &
     &   'STRAT across 500 hPa                    1278  24355  0',      &
     &   'TROP across 500 hPa                     1279  24356  0',      &
     &   'Be-7 across 500 hPa                     1274  24353  0',      &
     &   'Be-10 across 500 hPa                    1275  24354  0',      &
     &   'O3 across 600 hPa                       1311  24361  0',      &
     &   'O3(s) across 600 hPa                    1371  24362  0',      &
     &   'STRAT across 600 hPa                    1378  24365  0',      &
     &   'TROP across 600 hPa                     1379  24366  0',      &
     &   'Be-7 across 600 hPa                     1374  24363  0',      &
     &   'Be-10 across 600 hPa                    1375  24364  0',      &
     &   'O3 across 700 hPa                       1411  24371  0',      &
     &   'O3(s) across 700 hPa                    1471  24372  0',      &
     &   'STRAT across 700 hPa                    1478  24375  0',      &
     &   'TROP across 700 hPa                     1479  24376  0',      &
     &   'Be-7 across 700 hPa                     1474  24373  0',      &
     &   'Be-10 across 700 hPa                    1475  24374  0',      &
     &   'O3 across 800 hPa                       1511  24381  0',      &
     &   'O3(s) across 800 hPa                    1571  24382  0',      &
     &   'STRAT across 800 hPa                    1578  24385  0',      &
     &   'TROP across 800 hPa                     1579  24386  0',      &
     &   'Be-7 across 800 hPa                     1574  24383  0',      &
     &   'Be-10 across 800 hPa                    1575  24384  0',      &
     &   'O3 across 900 hPa                       1611  24391  0',      &
     &   'O3(s) across 900 hPa                    1671  24392  0',      &
     &   'STRAT across 900 hPa                    1678  24395  0',      &
     &   'TROP across 900 hPa                     1679  24396  0',      &
     &   'Be-7 across 900 hPa                     1674  24393  0',      &
     &   'Be-10 across 900 hPa                    1675  24394  0'/)

      CHARACTER(LEN=len_flx_str), DIMENSION(nozonenames), PARAMETER ::  &
     &  ozonenames=(/                                                   &
     &   'O3 Production:   HO2+NO                  350  24402  0',      &
     &   'O3 Production:   CH3O2+NO                351  24403  0',      &
     &   'O3 Production:   RO2+NO                  352  24404  0',      &
     &   'O3 Production:   Stratospheric input     353  24405  0',      &
     &   'O3 Production:   Total Production        354  24406  0',      &
     &   'O3 Loss:         O1D+H2O                 355  24407  0',      &
     &   'O3 Loss:         O3+OH                   356  24408  0',      &
     &   'O3 Loss:         O3+HO2                  357  24409  0',      &
     &   'O3 Loss:         O3+HCs                  358  24410  0',      &
     &   'O3 Loss:         Other                   359  24411  0',      &
     &   'O3 Loss:         O3 Dry Deposition       360  24412  0',      &
     &   'O3 Loss:         Total Loss              361  24413  0',      &
     &   'O3 Budget:       Net Production          362  24414  0',      &
     &   'O3 Budget:       Net Chemical Prodn.     363  24415  0'       &
     & /)

      CHARACTER(LEN=len_flx_str), DIMENSION(nextraflux), PARAMETER ::   &
     &  extranames=(/                                                   &
     &   'NO2 surface dry dep velocity            1705  25405  0',      &
     &   'O3 surface dry dep velocity             1711  25411  0',      &
     &   'PAN surface dry dep velocity            1721  25421  0',      &
     &   'SO2 surface dry dep velocity            1726  25426  0',      &
     &   'NH3 surface dry dep velocity            1770  25470  0',      &
     &   'CO surface dry dep velocity             1708  25408  0',      &
     &   'CH4 surface dry dep velocity            1709  25409  0',      &
     &   'H2 surface dry dep velocity             1712  25412  0',      &
     &   'HNO3 surface dry dep velocity           1713  25413  0',      &
     &   'H2O2 surface dry dep velocity           1714  25414  0',      &
     &   'CH3OOH surface dry dep velocity         1722  25422  0',      &
     &   'C2H5OOH surface dry dep velocity        1734  25434  0',      &
     &   'C3H7OOH surface dry dep velocity        1733  25433  0',      &
     &   'C4H9OOH surface dry dep velocity        1735  25435  0',      &
     &   'ISOPOOH surface dry dep velocity        1752  25452  0',      &
     &   'MVKOOH surface dry dep velocity         1753  25453  0',      &
     &   'SA surface dry dep velocity             1730  25430  0',      &
     &   'NH42SO4 surface dry dep velocity        1739  25439  0',      &
     &   'NAER surface dry dep velocity           1758  25458  0',      &
     &   'MSA surface dry dep velocity            1768  25468  0',      &
     &   'ORGNIT surface dry dep velocity         1769  25469  0',      &
     &   'Ozone deposition flux via stomata       1790  24421  0',      &
     &   'AOT40 ppb-h using surface layer only    1799  24422  0',      &
     &   'No. of cloud-to-ground lightning flashes1796  24431  0',      &
     &   'No. of intercloud lightning flashes     1797  24432  0'       &
     & /)

      CHARACTER(LEN=len_flx_str), DIMENSION(numflux), PARAMETER ::      &
     &  fluxnames=(/                                                    &
     &    photnames,                                                    &
     &    fluxnamesa,                                                   &
     &    fluxnamesb,                                                   &
     &    emissnames,                                                   &
     &    ddepnames,                                                    &
     &    wdepnames,                                                    &
     &    advfluxnames,                                                 &
     &    ozonenames,                                                   &
     &    extranames                                                    &
     &  /)

! Logical to turn on/off station data output
      LOGICAL, PARAMETER :: LSTAT=.FALSE.

! Number of stations to output concentrations for
      INTEGER, PARAMETER :: numstat = 10
      CHARACTER(LEN=40),DIMENSION(numstat),PARAMETER :: station_names = &
     &(/'Hohenpeissenberg      11.0     47.8',                          &
     &  'Garmisch              11.1     47.5',                          &
     &  'Payerne                7.0     46.8',                          &
     &  'S.P. Capufiume        11.6     44.6',                          &
     &  'Thessaloniki          23.0     40.5',                          &
     &  'Leibadi               23.3     40.5',                          &
     &  'Jungfraujoch           6.0     46.3',                          &
     &  'Mte. Cimone           10.4     44.1',                          &
     &  'Sonnblick             13.0     47.0',                          &
     &  'Zugspitze             11.0     47.4'/)

! Logical to turn on/off cell following data
      LOGICAL, PARAMETER :: lfol=.FALSE.

! Number of cells to output budgets for
      INTEGER, PARAMETER :: numfol=100
      INTEGER, DIMENSION(numfol), PARAMETER :: cell_follow=(/           &
     &      1,     1001,     2001,     3001,     4001,                  &
     &   5001,     6001,     7001,     8001,     9001,                  &
     &  10001,    11001,    12001,    13001,    14001,                  &
     &  15001,    16001,    17001,    18001,    19001,                  &
     &  20001,    21001,    22001,    23001,    24001,                  &
     &  25001,    26001,    27001,    28001,    29001,                  &
     &  30001,    31001,    32001,    33001,    34001,                  &
     &  35001,    36001,    37001,    38001,    39001,                  &
     &  40001,    41001,    42001,    43001,    44001,                  &
     &  45001,    46001,    47001,    48001,    49001,                  &
     &  50001,    51001,    52001,    53001,    54001,                  &
     &  55001,    56001,    57001,    58001,    59001,                  &
     &  60001,    61001,    62001,    63001,    64001,                  &
     &  65001,    66001,    67001,    68001,    69001,                  &
     &  70001,    71001,    72001,    73001,    74001,                  &
     &  75001,    76001,    77001,    78001,    79001,                  &
     &  80001,    81001,    82001,    83001,    84001,                  &
     &  85001,    86001,    87001,    88001,    89001,                  &
     &  90001,    91001,    92001,    93001,    94001,                  &
     &  95001,    96001,    97001,    98001,    99001/)

! Make HOURLY true for 1 hour pp output, false for monthly o/p.
      LOGICAL, PARAMETER :: lhourly=.FALSE.
      CHARACTER :: filetype2

      END MODULE IN_STOCHEM_OUT
#endif
