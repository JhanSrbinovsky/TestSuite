! TYPPTRA
! include TYPSIZE first
! History:
! Version  Date    Comment
!  3.4   18/05/94  Add pointers to new slab prognostics. J F Thomson.
!  3.5   19/05/95  Remove pointers JK1,JK2,JEXPK1,JEXPK2,JKDA,JKDF
!                  and JRHCRIT. D. Robinson
!  4.1   13/10/95  Add pointers for new soil moisture fraction
!                  prognostics,canopy conductance and
!                  vegetation J.Smith
!  4.1   26/04/96  Add pointers for Sulphur Cycle (12)   MJWoodage
!  4.3    18/3/97  Add pointers for HadCM2 sulphate loading patterns
!                                                   William Ingram
!  4.4   05/09/97  Add pointer for net energy flux prognostic
!                  S.D.Mullerworth
!  4.4   05/08/97  Add pointer for conv. cld amt on model levs. JMG
!  4.4  10/09/97   Added pointers for snow grain size and snow soot
!                  content used in prognostic snow albedo scheme
!                                                        R. Essery
!  4.4    16/9/97  Add pointers for new vegetation and land surface
!                  prognostics.                       Richard Betts
!  4.5    1/07/98  Add pointer for ocean CO2 flux. C.D.Jones
!  4.5    19/01/98 Replace JVEG_FLDS and JSOIL_FLDS with
!                  individual pointers. D. Robinson
!  4.5    04/03/98 Add 2 pointers for NH3 in S Cycle     M. Woodage
!                  Add 5 pointers for soot               M. Woodage
!                  Add pointer SO2_HILEM to CARGPT_ATMOS M. Woodage
!  4.5    08/05/98 Add 16 pointers for User Anc.         D. Goddard
!  4.5    15/07/98 Add pointers for new 3D CO2 array.    C.D.Jones
!  5.0    20/05/99 Extensive revision of pointers for new dynamics data
!                  variables. D.M. Goddard
!  5.1    25/02/00 Add new LBC pointers and labelled old JRIM
!                  pointers ready for removal.            P.Burton
!  5.1    28/04/00 Add extra level to p/exner on rho levels
!                                                 Andy Malcolm
!  5.1    29/02/00 Added JNET_MFLUX pointer for net moisture flux
!  5.1    03/05/00 Add JTSTAR_ANOM pointer to CARGPT_ATMOS. D Robinson
!  5.2    18/10/00 Inserted lines for slab model fields   K.D.Williams
!  5.1    31/08/00 Added JOROG_LBC pointer
!                  Move JEXNER_LBC to top of list
!                  Remove levels from the LBC J_pointers    P.Burton
!  5.2    25/09/00 Clear out diagnostic RHcrit legacy code. A.C.Bushell
!
!  5.2    13/09/00 Remove JSOOT, JSO4, JH2SO4 pointers. P.Selwood.
!  5.2    15/11/00 Re-Introduce pointers for MOSES 2 and Triffid
!                  M. Best.
!  5.2    09/03/01 Introduce extra pointers in level dependent constants
!                  for height definitions. R Rawlins
!  5.3    04/10/01 Added new slab prognostics             K.Williams
!  5.3       06/01 Introduce extra pointers needed for coastal
!                  tiling.                              Nic Gedney
!  5.3    19/06/01 Added JTPPSOZONE (tropopauase-based ozone) Dave Tan
!  5.3    01/11/01 Add J_CO2_EMITS and J_CO2FLUX to common block
!  5.4    02/09/02 Added new pointer JSNODEP_SEA        K.Williams
!  5.4    28/08/02 Add JCATCH_SNOW and JSNOW_GRND to common block
!  5.5    05/02/03 Add pointers for biomass aerosol scheme. P Davison
!  5.5    05/11/02 Large-scale hydrology scheme: Add pointers for
!                  gamma function, water table depth,
!                  surface saturated fraction and wetland fraction.
!                                                        N. Gedney
!  5.5    24/02/03 Add JICE_FRACT_CAT, JICE_THICK_CAT, JTI_CAT &
!                  JSNODEP_SEA_CAT to common block      J.Ridley
!  5.5    03/02/03 Add JQCF2,JQRAIN,JQGRAUP and LBC variables. R.Forbes
!  5.5    26/02/03 River routing support. P.Selwood.
!  6.1    07/09/04 Add pointers for STOCHEM fields.     C. Johnson
!  6.0    30/07/03 Add cloud fraction lbc variables. Damian Wilson
!  6.0    12/08/03 removes JRIV_GRIDBOX. C.Bunton.
!  6.0    12/09/03 Extra pointers for river routing 2A. V.A.Bell
!  6.1    07/04/04 Add pointer for DMS concentration.   A. Jones
!  6.2    24/02/06 Add pointer for DMS flux from ocean model    J.Gunson
!   6.2   21/2/06  Declare new pointer to
!                  re-route outflow from inland basins to soil moisture
!                  P. Falloon
!  6.2    01/03/06 Add pointer for direct component of total downward
!                  surface PAR flux.  M.G. Sanderson.
!  6.2    19/01/06 Add variables required for surface radiative
!                  forcing. Needed under versions 3C and 3Z of
!                  the radiation code.            (J.-C.Thelen)
!  6.2    01/10/05 Add JMURK_LBC variables.   R.M.Forbes
!  6.2    14/11/05 Add JTR_UKCA pointers for UKCA tracers. R Barnes
!  6.2    01/03/06 Add pointer for RothC pools and fluxes.   C.D. Jones
#if defined(ATMOS) || defined (A34_1A)
! Pointers for ATMOSPHERE model variables. Configuration dependent.
! Addresses in D1 array of model variables held in primary and
! secondary space.

      ! 1: Array Variables (dimensions are resolution dependent.)

      ! 1.1: Data variables stored in primary space.
      INTEGER :: JU(MODEL_LEVELS)      ! u component of wind
      INTEGER :: JV(MODEL_LEVELS)      ! v component of wind
      INTEGER :: JW(0:MODEL_LEVELS)    ! w component of wind
      INTEGER :: JRHO(MODEL_LEVELS)    ! Density
      INTEGER :: JTHETA(MODEL_LEVELS)  ! Potential temperature
      INTEGER :: JQ(WET_LEVELS)        ! Specific humidity
      INTEGER :: JQCL(WET_LEVELS)      ! qcl
      INTEGER :: JQCF(WET_LEVELS)      ! qcf
      INTEGER :: JQCF2(WET_LEVELS)     ! second ice
      INTEGER :: JQRAIN(WET_LEVELS)    ! rain
      INTEGER :: JQGRAUP(WET_LEVELS)   ! graupel

      ! Exner pressure on rho levels
      INTEGER :: JEXNER_RHO_LEVELS(MODEL_LEVELS+1)

      INTEGER :: JU_ADV(MODEL_LEVELS)   ! Advective u component of wind
      INTEGER :: JV_ADV(MODEL_LEVELS)   ! Advective v component of wind
      INTEGER :: JW_ADV(0:MODEL_LEVELS) ! Advective w component of wind

      ! 1.2: Data variables stored in secondary space.
      INTEGER :: JP(MODEL_LEVELS+1)                  ! Pressure on rho le

      ! Pressure on theta levels
      INTEGER :: JP_THETA_LEVELS(MODEL_LEVELS)

      ! Exner pressure on theta levels
      INTEGER :: JEXNER_THETA_LEVELS(MODEL_LEVELS)

      ! 1.3: Cloud Fields
      INTEGER :: JCCW_RAD(WET_LEVELS)            ! CCW profile to radiation
      INTEGER :: JCCA(N_CCA_LEV)                 ! Convective cloud amount
      INTEGER :: JCF_AREA(WET_LEVELS)            ! Area Cloud Fraction
      INTEGER :: JCF_BULK(WET_LEVELS)            ! Bulk Cloud Fraction
      INTEGER :: JCF_LIQUID(WET_LEVELS)          ! Liquid cloud fraction
      INTEGER :: JCF_FROZEN(WET_LEVELS)          ! Frozen cloud fraction

      ! 1.4: Soil Ancillary fields
      INTEGER :: J_DEEP_SOIL_TEMP(ST_LEVELS)   ! Deep soil temperature

      INTEGER :: JSMCL(SM_LEVELS)   ! Soil moisture content in layers
      INTEGER :: JSTHU(SM_LEVELS)   ! Unfrozen soil moisture fraction
      INTEGER :: JSTHF(SM_LEVELS)   ! Frozen soil moisture fraction

      ! 1.5: Radiation Increments
      INTEGER :: JSW_INCS(0:MODEL_LEVELS+1)    ! SW radiation increments
      INTEGER :: JLW_INCS(0:MODEL_LEVELS)      ! LW radiation increments
! PAR radiation increment
      INTEGER :: JDIRPAR

      ! 1.6: Ozone
      INTEGER :: JOZONE(OZONE_LEVELS)          ! Ozone
!  tropopause-based ozone
      INTEGER :: JTPPSOZONE(tpps_ozone_levels)
      ! 1.7: Tracer and aerosol fields
      INTEGER :: JTRACER(TR_LEVELS,TR_VARS+1)  ! Tracers
      INTEGER :: JTR_UKCA(TR_LEVELS,TR_UKCA+1)  ! UKCA Tracers
      INTEGER :: JMURK_SOURCE(MODEL_LEVELS)    ! multi-level murk source
      INTEGER :: JMURK(MODEL_LEVELS)           ! multi-level murk concent
      INTEGER :: JDUST_DIV1(MODEL_LEVELS)      ! dust mmr, division 1
      INTEGER :: JDUST_DIV2(MODEL_LEVELS)      ! dust mmr, division 2
      INTEGER :: JDUST_DIV3(MODEL_LEVELS)      ! dust mmr, division 3
      INTEGER :: JDUST_DIV4(MODEL_LEVELS)      ! dust mmr, division 4
      INTEGER :: JDUST_DIV5(MODEL_LEVELS)      ! dust mmr, division 5
      INTEGER :: JDUST_DIV6(MODEL_LEVELS)      ! dust mmr, division 6
      INTEGER :: JSO2(MODEL_LEVELS)            ! sulphur dioxide gas
      INTEGER :: JDMS(MODEL_LEVELS)            ! dimethyl sulphide gas
      INTEGER :: JSO4_AITKEN(MODEL_LEVELS)     ! Aitken mode sulphate aer
      INTEGER :: JSO4_ACCU(MODEL_LEVELS)       ! accumulation mode sulpha
      INTEGER :: JSO4_DISS(MODEL_LEVELS)       ! dissloved  sulphate aero
      INTEGER :: JH2O2(MODEL_LEVELS)           ! hydrogen peroxide mmr
      INTEGER :: JNH3(MODEL_LEVELS)            ! ammonia gas mmr
      INTEGER :: JSOOT_NEW(MODEL_LEVELS)       ! fresh soot mmr
      INTEGER :: JSOOT_AGD(MODEL_LEVELS)       ! aged soot mmr
      INTEGER :: JSOOT_CLD(MODEL_LEVELS)       ! soot in cloud mmr
      INTEGER :: JBMASS_NEW(MODEL_LEVELS)      ! fresh biomass mmr
      INTEGER :: JBMASS_AGD(MODEL_LEVELS)      ! aged biomass mmr
      INTEGER :: JBMASS_CLD(MODEL_LEVELS)      ! cloud biomass mmr
      INTEGER :: JOCFF_NEW(MODEL_LEVELS)       ! fresh OCFF mmr
      INTEGER :: JOCFF_AGD(MODEL_LEVELS)       ! aged OCFF mmr
      INTEGER :: JOCFF_CLD(MODEL_LEVELS)       ! OCFF in cloud mmr
      INTEGER :: JSO2_NATEM(MODEL_LEVELS)      ! natural SO2 emissions
      INTEGER :: JOH(MODEL_LEVELS)             ! hydroxyl radical ancilla
      INTEGER :: JHO2(MODEL_LEVELS)            ! hydrogen dioxide ancilla
      INTEGER :: JH2O2_LIMIT(MODEL_LEVELS)     ! limiting H2O2 ancillary
      INTEGER :: JO3_CHEM(MODEL_LEVELS)        ! ozone for chemistry anci
      INTEGER :: JHadCM2_SO4(2)                ! HadCM2 sulphate loading
      ! JHadCM2_SO4: Should really be NSULPHAT (= 2) but to use NSULPAT,
      !              must include CSENARIO before every include of TYPPTR
      INTEGER :: JCO2(MODEL_LEVELS)            ! 3D CO2 FIELD
      INTEGER :: JCH4_STOCH(MODEL_LEVELS)      ! STOCHEM CH4
      INTEGER :: JO3_STOCH(MODEL_LEVELS)       ! STOCHEM O3
      INTEGER :: JOH_UKCA(MODEL_LEVELS)        ! OH MMR from UKCA
      INTEGER :: JHO2_UKCA(MODEL_LEVELS)       ! HO2 MMR from UKCA
      INTEGER :: JH2O2_UKCA(MODEL_LEVELS)      ! H2O2 MMR from UKCA
      INTEGER :: JO3_UKCA(MODEL_LEVELS)        ! O3 MMR from UKCA
      INTEGER :: JOZONE_TRACER(MODEL_LEVELS)   ! Prognostic O3 Tracer(Cariol)
      INTEGER :: JO3_PROD_LOSS(MODEL_LEVELS)   ! Cariol O3 Prod-Loss (P-L)
      INTEGER :: JO3_P_L_VMR(MODEL_LEVELS)     ! Cariol O3 P-L wrt VMR 
      INTEGER :: JO3_VMR(MODEL_LEVELS)         ! Cariol O3 Vol Mix Ratio-VMR
      INTEGER :: JO3_P_L_TEMP(MODEL_LEVELS)    ! Cariol O3 P-L wrt temp 
      INTEGER :: JO3_TEMP(MODEL_LEVELS)        ! Cariol O3 temp  
      INTEGER :: JO3_P_L_COLO3(MODEL_LEVELS)   ! Cariol O3 P-L wrt colO3 
      INTEGER :: JO3_COLO3(MODEL_LEVELS)       ! Cariol O3 column (colO3) 
      INTEGER :: JARCLBIOG_BG(MODEL_LEVELS)    ! Biogenic aerosol climatology 
      INTEGER :: JARCLBIOM_FR(MODEL_LEVELS)    ! Biomass burning (fresh) aerosol clim
      INTEGER :: JARCLBIOM_AG(MODEL_LEVELS)    ! Biomass burning (aged) aerosol clim
      INTEGER :: JARCLBIOM_IC(MODEL_LEVELS)    ! Biomass burning (in-cloud) aerosol clim
      INTEGER :: JARCLBLCK_FR(MODEL_LEVELS)    ! Black carbon (fresh) aerosol clim
      INTEGER :: JARCLBLCK_AG(MODEL_LEVELS)    ! Black carbon (aged) aerosol clim
      INTEGER :: JARCLSSLT_FI(MODEL_LEVELS)    ! Sea salt (film mode) aerosol clim 
      INTEGER :: JARCLSSLT_JT(MODEL_LEVELS)    ! Sea salt (jet mode) aerosol clim
      INTEGER :: JARCLSULP_AC(MODEL_LEVELS)    ! Sulphate (accumulation mode) aero clim
      INTEGER :: JARCLSULP_AK(MODEL_LEVELS)    ! Sulphate (Aitken mode) aerosol clim 
      INTEGER :: JARCLSULP_DI(MODEL_LEVELS)    ! Sulphate (dissolved) aerosol clim
      INTEGER :: JARCLDUST_B1(MODEL_LEVELS)    ! Dust (bin 1) aerosol climatology 
      INTEGER :: JARCLDUST_B2(MODEL_LEVELS)    ! Dust (bin 2) aerosol climatology 
      INTEGER :: JARCLDUST_B3(MODEL_LEVELS)    ! Dust (bin 3) aerosol climatology 
      INTEGER :: JARCLDUST_B4(MODEL_LEVELS)    ! Dust (bin 4) aerosol climatology 
      INTEGER :: JARCLDUST_B5(MODEL_LEVELS)    ! Dust (bin 5) aerosol climatology 
      INTEGER :: JARCLDUST_B6(MODEL_LEVELS)    ! Dust (bin 6) aerosol climatology 
      INTEGER :: JARCLOCFF_FR(MODEL_LEVELS)    ! Org carbon fossil fuel (fresh) aero clim
      INTEGER :: JARCLOCFF_AG(MODEL_LEVELS)    ! Org carbon fossil fuel (aged) aero clim
      INTEGER :: JARCLOCFF_IC(MODEL_LEVELS)    ! Org carbon fossil fuel (in-cloud) aero clim
      INTEGER :: JARCLDLTA_DL(MODEL_LEVELS)    ! Delta aerosol climatology
!

! 1.8: Multi-level user ancillary fields
      INTEGER :: JUSER_MULT1(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT2(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT3(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT4(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT5(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT6(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT7(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT8(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT9(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT10(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT11(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT12(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT13(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT14(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT15(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT16(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT17(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT18(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT19(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT20(MODEL_LEVELS)    ! multi-level user ancilla

      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
      INTEGER :: JOROG_LBC                       ! Orography LBC
      INTEGER :: JU_LBC                          ! U LBC
      INTEGER :: JV_LBC                          ! V LBC
      INTEGER :: JW_LBC                          ! W LBC
      INTEGER :: JRHO_LBC                        ! RHO LBC
      INTEGER :: JTHETA_LBC                      ! Theta LBC
      INTEGER :: JQ_LBC                          ! Q LBC
      INTEGER :: JQCL_LBC                        ! QCL LBC
      INTEGER :: JQCF_LBC                        ! QCF LBC
      INTEGER :: JQCF2_LBC                       ! 2nd Ice LBC
      INTEGER :: JQRAIN_LBC                      ! Rain LBC
      INTEGER :: JQGRAUP_LBC                     ! Graupel LBC
      INTEGER :: JCF_BULK_LBC                    ! CF_BULK LBC
      INTEGER :: JCF_LIQUID_LBC                  ! CF_LIQUID_LBC
      INTEGER :: JCF_FROZEN_LBC                  ! CF_FROZEN_LBC
      INTEGER :: JEXNER_LBC                      ! Exner LBC
      INTEGER :: JU_ADV_LBC                      ! U_ADV LBC
      INTEGER :: JV_ADV_LBC                      ! V_ADV LBC
      INTEGER :: JW_ADV_LBC                      ! W_ADV LBC
      INTEGER :: JMURK_LBC                       ! Murk aerosol LBC
      INTEGER :: JTRACER_LBC(TR_VARS+1)          ! Tracer LBCs
      ! Lateral Boundary Condition tendencies
      INTEGER :: JU_LBC_TEND                     ! U LBC  tendencies
      INTEGER :: JV_LBC_TEND                     ! V LBC tendencies
      INTEGER :: JW_LBC_TEND                     ! W LBC tendencies
      INTEGER :: JRHO_LBC_TEND                   ! RHO LBC tendencies
      INTEGER :: JTHETA_LBC_TEND                 ! Theta LBC tendencies
      INTEGER :: JQ_LBC_TEND                     ! Q LBC tendencies
      INTEGER :: JQCL_LBC_TEND                   ! QCL LBC tendencies
      INTEGER :: JQCF_LBC_TEND                   ! QCF LBC tendencies
      INTEGER :: JQCF2_LBC_TEND                  ! 2nd Ice
      INTEGER :: JQRAIN_LBC_TEND                 ! Rain LBC tendencies
      INTEGER :: JQGRAUP_LBC_TEND                ! Graupel
      INTEGER :: JCF_BULK_LBC_TEND               ! CF_BULK LBC tend'cies
      INTEGER :: JCF_LIQUID_LBC_TEND             ! CF_LIQUID_LBC t'cies
      INTEGER :: JCF_FROZEN_LBC_TEND             ! CF_FROZEN_LBC t'cies
      INTEGER :: JEXNER_LBC_TEND                 ! Exner LBC tendencies
      INTEGER :: JU_ADV_LBC_TEND                 ! U_ADV LBC tendencies
      INTEGER :: JV_ADV_LBC_TEND                 ! V_ADV LBC tendencies
      INTEGER :: JW_ADV_LBC_TEND                 ! W_ADV LBCtendencies
      INTEGER :: JMURK_LBC_TEND                  ! Murk aerosol LBC tend
      INTEGER :: JTRACER_LBC_TEND(TR_VARS+1)     ! Tracer LBC tendencies

      ! 2: Scalar Variables

      ! 2.1: Data variables stored in primary space.
      INTEGER :: JTSTAR         ! Surface temperature
      INTEGER :: JLAND          ! Land sea mask
      INTEGER :: JTSTAR_ANOM    ! Surface temperature anolomy
!   2.15: Fields for coastal tiling
      INTEGER :: JFRAC_LAND  ! Land fraction in grid box
      INTEGER :: JTSTAR_LAND ! Land surface temperature
      INTEGER :: JTSTAR_SEA  ! Sea surface temperature
      INTEGER :: JTSTAR_SICE ! Sea-ice surface temperature
! Set pointers for sea-ice and land albedos
      INTEGER :: JSICE_ALB                   ! Sea-ice albedo
      INTEGER :: JLAND_ALB                   ! Mean land albedo

      ! 2.2: Data variables stored in secondary space.

      INTEGER :: JPSTAR          ! Surface pressure

      ! 2.3: Cloud fields
      INTEGER :: JLCBASE         ! Lowest Convective cloud base
      INTEGER :: JCCB            ! Convective cloud base
      INTEGER :: JCCT            ! Convective cloud top

      INTEGER :: JCCLWP          ! Convective cloud liquid water path

      ! 2.4: Boundary layer fields

      INTEGER :: JZH                         ! Boundary layer depth

      ! Standard deviation of turbulent fluctuations of layer 1
                                    ! temperature
      INTEGER :: JT1_SD

      ! Standard deviation of turbulent fluctuations of layer 1 humidity
      INTEGER :: JQ1_SD

      ! Number of model levels in the  turbulently mixed layer
      INTEGER :: JNTML

      ! Top level for turb mixing in any decoupled Sc layer
      INTEGER :: JNTDSC

      ! Bottom level for turb mixing in any decoupled Sc layer
      INTEGER :: JNBDSC

      INTEGER :: JCUMULUS      ! Boundary layer convection flag

      ! 2.4: Soil Ancillary fields

      INTEGER :: JSAT_SOILW_SUCTION          ! Saturated soil water sucti
      INTEGER :: JTHERM_CAP     ! Thermal capacity
      INTEGER :: JTHERM_COND    ! Thermal conductivity
      INTEGER :: JVOL_SMC_CRIT  ! Vol smc at critical point
      INTEGER :: JVOL_SMC_WILT  ! Vol smc at wilting point
      INTEGER :: JVOL_SMC_SAT   ! Vol smc at saturation
      INTEGER :: JSAT_SOIL_COND ! Saturated soil conductivity
      INTEGER :: JCLAPP_HORN    ! Clapp-Hornberger B coefficient

      ! 2.5: Vegetation Ancillary fields

      INTEGER :: JCANOPY_WATER  ! Canopy Water
      INTEGER :: JSURF_CAP      ! Surface Capacity
      INTEGER :: JSURF_RESIST   ! Surface Resistance
      INTEGER :: JROOT_DEPTH    ! Root depth
      INTEGER :: JINFILT        ! Infiltration factor
      INTEGER :: JVEG_FRAC      ! Vegetation fraction
      INTEGER :: JLAI           ! Leaf area index
      INTEGER :: JCANHT         ! Canopy height
      INTEGER :: JZ0            ! Vegetative Roughness length
      INTEGER :: JSFA           ! Snow free albedo
      INTEGER :: JMDSA          ! Cold deep snow albedo
      INTEGER :: JGS            ! Gridbox mean canopy conductance

      ! 2.6: Orographic Ancillary fields

      INTEGER :: JOROG          ! Orographic height
      INTEGER :: JOROG_SD       ! Standard Deviation of orography
      INTEGER :: JOROG_SIL      ! Silhouette area of orography
      INTEGER :: JOROG_HO2      ! Peak to trough height/(2*sqrt2)
      INTEGER :: JOROG_GRAD_X   ! Orographic gradient x
      INTEGER :: JOROG_GRAD_Y   ! Orographic gradient y
      INTEGER :: JOROG_GRAD_XX  ! Orographic gradient xx
      INTEGER :: JOROG_GRAD_XY  ! Orographic gradient xy
      INTEGER :: JOROG_GRAD_YY  ! Orographic gradient yy

      ! 2.7: Sea/Sea Ice

      INTEGER :: JU_SEA         ! Surface current (u component)
      INTEGER :: JV_SEA         ! Surface current (v component)
      INTEGER :: JU_0_P         ! Surace  current (u) on p-grid
      INTEGER :: JV_0_P         ! Surface current (v) on p-grid
      INTEGER :: JICE_FRACTION  ! Sea ice fraction
      INTEGER :: JICE_THICKNESS ! Sea ice depth
      INTEGER :: JTI            ! Sea ice temperature
      INTEGER :: JICE_FRACT_CAT ! Sea ice fraction on catagories
      INTEGER :: JICE_THICK_CAT ! Sea ice thickness on catagories
      INTEGER :: JTI_CAT        ! Sea ice temperature on catagories

      ! 2.8: Snow

      INTEGER :: JSNODEP        ! Snow depth on land
      INTEGER :: JSNODEP_SEA    ! Snow depth on sea ice
      INTEGER :: JSNODEP_SEA_CAT ! Snow depth on sea ice catagories
      INTEGER :: JCATCH_SNOW    ! Coniferous canopy
!                               ! snow capacity
      INTEGER :: JSNOW_GRND     ! Snow below canopy
      INTEGER :: JSNSOOT        ! Snow soot content

! 2.9: aerosol emission fields, including mineral dust parent soil props

      INTEGER :: JSOIL_CLAY                    ! soil clay fraction
      INTEGER :: JSOIL_SILT                    ! soil silt fraction
      INTEGER :: JSOIL_SAND                    ! soil sand fraction
      INTEGER :: JDUST_MREL1                   ! soil rel mass, div 1
      INTEGER :: JDUST_MREL2                   ! soil rel mass, div 2
      INTEGER :: JDUST_MREL3                   ! soil rel mass, div 3
      INTEGER :: JDUST_MREL4                   ! soil rel mass, div 4
      INTEGER :: JDUST_MREL5                   ! soil rel mass, div 5
      INTEGER :: JDUST_MREL6                   ! soil rel mass, div 6


      INTEGER :: JSO2_EM        ! sulphur dioxide emission
      INTEGER :: JDMS_EM        ! dimethyl sulphide emission
      INTEGER :: JSO2_HILEM     ! high level SO2 emissions
      INTEGER :: JNH3_EM        ! ammonia gas surface emiss
      INTEGER :: JSOOT_EM       ! fresh soot surface emissions
      INTEGER :: JSOOT_HILEM    ! fresh soot high lev emissions
      INTEGER :: JBMASS_EM       ! fresh bmass surface emissions
      INTEGER :: JBMASS_HILEM    ! fresh bmass high lev emissions
      INTEGER :: JOCFF_EM        ! fresh OCFF surface emissions
      INTEGER :: JOCFF_HILEM     ! fresh OCFF high-level emissions
      INTEGER :: JDMS_CONC       ! seawater dimethyl sulphide conc.
      INTEGER :: JDMS_OFLUX      ! DMS flux from ocean model

      ! 2.10: User ancillary fields
      INTEGER :: JUSER_ANC1         ! user ancillary field 1
      INTEGER :: JUSER_ANC2                  ! user ancillary field 2
      INTEGER :: JUSER_ANC3                  ! user ancillary field 3
      INTEGER :: JUSER_ANC4                  ! user ancillary field 4
      INTEGER :: JUSER_ANC5                  ! user ancillary field 5
      INTEGER :: JUSER_ANC6                  ! user ancillary field 6
      INTEGER :: JUSER_ANC7                  ! user ancillary field 7
      INTEGER :: JUSER_ANC8                  ! user ancillary field 8
      INTEGER :: JUSER_ANC9                  ! user ancillary field 9
      INTEGER :: JUSER_ANC10                 !user ancillary field 10
      INTEGER :: JUSER_ANC11                 ! user ancillary field 11
      INTEGER :: JUSER_ANC12                 ! user ancillary field 12
      INTEGER :: JUSER_ANC13                 ! user ancillary field 13
      INTEGER :: JUSER_ANC14                 ! user ancillary field 14
      INTEGER :: JUSER_ANC15                 ! user ancillary field 15
      INTEGER :: JUSER_ANC16                 ! user ancillary field 16
      INTEGER :: JUSER_ANC17                 ! user ancillary field 17
      INTEGER :: JUSER_ANC18                 ! user ancillary field 18
      INTEGER :: JUSER_ANC19                 ! user ancillary field 19
      INTEGER :: JUSER_ANC20                 ! user ancillary field 20

      ! Tracer Fluxes - kdcorbin, 05/10
      INTEGER :: JTRACER_FLUX1
      INTEGER :: JTRACER_FLUX2
      INTEGER :: JTRACER_FLUX3
      INTEGER :: JTRACER_FLUX4
      INTEGER :: JTRACER_FLUX5
      INTEGER :: JTRACER_FLUX6
      INTEGER :: JTRACER_FLUX7
      INTEGER :: JTRACER_FLUX8
      INTEGER :: JTRACER_FLUX9
      INTEGER :: JTRACER_FLUX10
      INTEGER :: JTRACER_FLUX11
      INTEGER :: JTRACER_FLUX12
      INTEGER :: JTRACER_FLUX13
      INTEGER :: JTRACER_FLUX14
      INTEGER :: JTRACER_FLUX15
      INTEGER :: JTRACER_FLUX16
      INTEGER :: JTRACER_FLUX17
      INTEGER :: JTRACER_FLUX18
      INTEGER :: JTRACER_FLUX19
      INTEGER :: JTRACER_FLUX20  

      !   2.11: Store arrays for energy correction calculation
      INTEGER :: JNET_FLUX                   ! Net energy flux
      INTEGER :: JNET_MFLUX                  ! Net moisture flux

      !   2.12: Tiled Vegetation and Triffid fields
      INTEGER :: JFRAC_TYP        ! Fractions of surface type
      INTEGER :: JFRAC_CON1       ! Fractions of surface type
      INTEGER :: JFRAC_CON2       ! Fractions of surface type
      INTEGER :: JFRAC_CON3       ! Fractions of surface type
      INTEGER :: JFRAC_CON4       ! Fractions of surface type
      INTEGER :: JFRAC_CON5       ! Fractions of surface type
      INTEGER :: JFRAC_CON6       ! Fractions of surface type
      INTEGER :: JFRAC_CON7       ! Fractions of surface type
      INTEGER :: JFRAC_CON8       ! Fractions of surface type
      INTEGER :: JFRAC_CON9       ! Fractions of surface type
      INTEGER :: JLAI_PFT         ! LAI of plant functional types
      INTEGER :: JCANHT_PFT       ! Canopy hght of plant func types
      INTEGER :: JDISTURB         ! Disturbed fraction of vegetation
      INTEGER :: JSOIL_ALB        ! Snow-free albedo of bare soil
      INTEGER :: JSOIL_CARB       ! Soil carbon content
      INTEGER :: JSOIL_CARB1      ! Soil carbon content DPM
      INTEGER :: JSOIL_CARB2      ! Soil carbon content RPM
      INTEGER :: JSOIL_CARB3      ! Soil carbon content BIO
      INTEGER :: JSOIL_CARB4      ! Soil carbon content HUM
      INTEGER :: JNPP_PFT_ACC     ! Accumulated NPP on PFTs
      INTEGER :: JG_LF_PFT_ACC    ! Accum. leaf turnover rate PFTs
      INTEGER :: JG_PHLF_PFT_ACC  ! Accumulated phenological leaf
                                    ! turnover rate on PFTs
      INTEGER :: JRSP_W_PFT_ACC   ! Accum. wood respiration on PFTs
      INTEGER :: JRSP_S_ACC       ! Accumulated soil respiration
      INTEGER :: JRSP_S_ACC1      ! Accumulated soil respiration DPM
      INTEGER :: JRSP_S_ACC2      ! Accumulated soil respiration RPM
      INTEGER :: JRSP_S_ACC3      ! Accumulated soil respiration BIO
      INTEGER :: JRSP_S_ACC4      ! Accumulated soil respiration HUM
      INTEGER :: JCAN_WATER_TILE  ! Canopy water content on tiles
      INTEGER :: JCATCH_TILE      ! Canopy capacity on tiles
      INTEGER :: JINFIL_TILE      ! Max infiltration rate on tiles
      INTEGER :: JRGRAIN_TILE     ! Snow grain size on tiles
      INTEGER :: JSNODEP_TILE     ! Snow depth on tiles
      INTEGER :: JTSTAR_TILE      ! Surface temperature on tiles
      INTEGER :: JZ0_TILE         ! Surface roughness on tiles
      INTEGER :: JDOLR            ! TOA - surface upward LW at
                                  ! radiation timestep
      INTEGER :: JLW_DOWN         ! Surface downward LW at
                                  ! radiation timestep
      INTEGER :: JSW_TILE         ! Surface net SW on land tiles at
                                  ! radiation timestep
      !   2.13: Slab Model
      INTEGER :: JTSLAB           ! Temperature of slab ocean.
      INTEGER :: JTCLIM           ! Climatological SST's
      INTEGER :: JHCLIM           ! Climatological ice depth
      INTEGER :: JCHEAT     ! Caryheat (from ice model to ocn)
      INTEGER :: JOIFLX     ! Ocean-Ice heat flux
      INTEGER :: JUICE      ! Zonal comp of ice velocity
      INTEGER :: JVICE      ! Meridional comp of ice velocity
      INTEGER :: JSIG11NE   !  Internal stresses for
      INTEGER :: JSIG11SE   !  EVP ice dynamics
      INTEGER :: JSIG11SW
      INTEGER :: JSIG11NW
      INTEGER :: JSIG12NE
      INTEGER :: JSIG12SE
      INTEGER :: JSIG12SW
      INTEGER :: JSIG12NW
      INTEGER :: JSIG22NE
      INTEGER :: JSIG22SE
      INTEGER :: JSIG22SW
      INTEGER :: JSIG22NW

!   2.14: Carbon cycle fields
      INTEGER :: J_CO2FLUX   ! Ocean CO2 flux (Kg CO2/m2/s1)
      INTEGER :: J_CO2_EMITS ! Surface CO2 emissions (Kg CO2/m2/s1)

!   2.15: Fields carried forward from previous version.
!         May not be required
      INTEGER :: JSURF_RESIST_NIT  ! Surface resistance on
                                    ! non-ice tiles
      INTEGER :: JROOT_DEPTH_NIT   ! Root depth on non-ice tiles
      INTEGER :: JZ0V_TYP          ! Vegetative Roughness length on
                                    ! tiles
      INTEGER :: JTSNOW            ! Snow surface layer temperature
      INTEGER :: JICE_EDGE
      INTEGER :: JOROG_TENDENCY    ! Orographic tendencies
      INTEGER :: JOROG_SD_TENDENCY ! Orographic variable tendency

      ! Pointers for ATMOSPHERE model constants. Scalars only.
      ! Addresses in level dependent constants array.
      INTEGER :: JETATHETA
      INTEGER :: JETARHO
      INTEGER :: JRHCRIT
      INTEGER :: JSOIL_THICKNESS
      ! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
      INTEGER :: Jzseak_theta ! zsea(k) on theta levels
      INTEGER :: JCk_theta    ! C(k)    on theta levels
      INTEGER :: Jzseak_rho   ! zsea(k) on rho levels
      INTEGER :: JCk_rho      ! C(k)    on rho levels
      ! Addresses in Row and Col  dependent constants array.
      INTEGER :: JLAMBDA_INPUT_P
      INTEGER :: JLAMBDA_INPUT_U
      INTEGER :: JPHI_INPUT_P
      INTEGER :: JPHI_INPUT_V
!   2.16: Fields for large-scale hydrology scheme.

      INTEGER :: JTI_MEAN          !Mean topographic index
      INTEGER :: JTI_SIG           !Standard dev. in topographic index
      INTEGER :: JFEXP             !Exponential decay in soil
!                                  ! saturated conductivity
      INTEGER :: JGAMMA_INT        !Integrated gamma function
      INTEGER :: JWATER_TABLE      !Water table depth
      INTEGER :: JFSFC_SAT         !Surface saturation fraction
      INTEGER :: JF_WETLAND        !Wetland fraction
      INTEGER :: JSTHZW           !Soil moist fract. in deep-zw layer.
      INTEGER :: JA_FSAT          !Fitting parameter for Fsat in LSH.
      INTEGER :: JC_FSAT          !Fitting parameter for Fsat in LSH.
      INTEGER :: JA_FWET          !Fitting parameter for Fwet in LSH.
      INTEGER :: JC_FWET          !Fitting parameter for Fwet in LSH.


!   2.17: Fields for River routing.
      INTEGER :: JRIV_SEQUENCE   ! River sequence
      INTEGER :: JRIV_DIRECTION  ! River direction
      INTEGER :: JRIV_STORAGE    ! River water storage
      INTEGER :: JTOT_SURFROFF   ! Accumulated surface runoff
      INTEGER :: JTOT_SUBROFF    !     "       sub-surface runoff
      INTEGER :: JRIV_INLANDATM       ! inland basin outflow
! Fields for grid-to-grid river routing (river routing 2A)
      INTEGER :: JRIV_IAREA      ! Drainage area
      INTEGER :: JRIV_SLOPE      ! Grid-cell slope
      INTEGER :: JRIV_FLOWOBS1   ! Initial values of flow
      INTEGER :: JRIV_INEXT      ! Flow direction (x)
      INTEGER :: JRIV_JNEXT      ! Flow direction (y)
      INTEGER :: JRIV_LAND       ! Land-type (land/river/sea)
      INTEGER :: JRIV_SUBSTORE   ! Subsurface storage
      INTEGER :: JRIV_SURFSTORE  ! Surface storage
      INTEGER :: JRIV_FLOWIN     ! Surface inflow
      INTEGER :: JRIV_BFLOWIN    ! Subsurface inflow

      INTEGER :: JC_SOLAR 
      INTEGER :: JC_BLUE 
      INTEGER :: JC_DOWN 
      INTEGER :: JC_LONGWAVE 
      INTEGER :: JC_TAUX 
      INTEGER :: JC_TAUY 
      INTEGER :: JC_WINDMIX 
      INTEGER :: JC_SENSIBLE 
      INTEGER :: JC_SUBLIM 
      INTEGER :: JC_EVAP 
      INTEGER :: JC_BOTMELTN 
      INTEGER :: JC_TOPMELTN 
      INTEGER :: JC_LSRAIN 
      INTEGER :: JC_LSSNOW 
      INTEGER :: JC_CVRAIN 
      INTEGER :: JC_CVSNOW 
      INTEGER :: JC_RIVEROUT 

      INTEGER :: JC_PRESS

!    Fields for CABLE
      INTEGER :: JTSOIL_TILE(SM_LEVELS)  ! Tiled soil temperature
      INTEGER :: JSMCL_TILE(SM_LEVELS)   ! Tiled soil moisture content in layers
! Not needed MRD
!!!      INTEGER :: JSTHU_TILE(SM_LEVELS)   ! Tiled unfrozen soil moisture fraction
      INTEGER :: JSTHF_TILE(SM_LEVELS)   ! Tiled frozen soil moisture fraction
      INTEGER :: JSNOW_DEPTH3L(3)        ! Tiled snow depth
      INTEGER :: JSNOW_MASS3L(3)         ! Tiled snow mass
      INTEGER :: JSNOW_TMP3L(3)          ! Tiled snow temperature
      INTEGER :: JSNOW_RHO3L(3)          ! Tiled snow density
      INTEGER :: JSNOW_RHO1L             ! Tiled mean snow density
      INTEGER :: JSNOW_AGE               ! Tiled snow age
      INTEGER :: JSNOW_FLG3L             ! Flag for use of 3 level snow scheme
! Lestevens Sept 2012 - adding new progs for CASA-CNP
      INTEGER :: JCPOOL_TILE(10)
      INTEGER :: JNPOOL_TILE(10)
      INTEGER :: JPPOOL_TILE(12)
      INTEGER :: JSOIL_ORDER
      INTEGER :: JNIDEP
      INTEGER :: JNIFIX
      INTEGER :: JPWEA
      INTEGER :: JPDUST
      INTEGER :: JGLAI
      INTEGER :: JPHENPHASE
      INTEGER :: JCABLE_LAI(1)   

! Addresses in D1 array of primary variables: scalars
      COMMON/CARGPT_ATMOS/                                              &
! Data variables
     &  JTSTAR, JLAND, JPSTAR, JTSTAR_ANOM,                             &
     &  JFRAC_LAND, JTSTAR_LAND, JTSTAR_SEA, JTSTAR_SICE,               &
     &  JSICE_ALB, JLAND_ALB,                                           &
      ! Cloud fields
     &  JCCB, JCCT, JCCLWP,                                             &
      ! Boundary layer fields
     &  JZH, JT1_SD, JQ1_SD, JNTML, JNTDSC, JNBDSC, JCUMULUS,           &
      ! Soil Ancillary fields
     &  JSAT_SOILW_SUCTION, JTHERM_CAP, JTHERM_COND,                    &
     &  JVOL_SMC_CRIT, JVOL_SMC_WILT, JVOL_SMC_SAT,                     &
     &  JSAT_SOIL_COND, JCLAPP_HORN,                                    &
      ! Vegetation Ancillary fields
     &  JCANOPY_WATER, JSURF_CAP, JSURF_RESIST, JROOT_DEPTH,            &
     &  JINFILT, JVEG_FRAC, JLAI, JCANHT,                               &
     &  JZ0, JSFA, JMDSA, JGS,                                          &
      ! Orographic Ancillary fields
     &  JOROG, JOROG_SD, JOROG_SIL, JOROG_HO2,                          &
     &  JOROG_GRAD_X, JOROG_GRAD_Y,                                     &
     &  JOROG_GRAD_XX, JOROG_GRAD_XY, JOROG_GRAD_YY,                    &
      ! Sea/Sea Ice
     &  JU_SEA, JV_SEA, JU_0_P, JV_0_P,                                 &
     &  JICE_FRACTION, JICE_THICKNESS, JTI,                             &
     &  JICE_FRACT_CAT, JICE_THICK_CAT, JTI_CAT,                        &
      ! Snow
     &  JSNODEP, JSNODEP_SEA, JSNSOOT, JCATCH_SNOW, JSNOW_GRND,         &
     &  JSNODEP_SEA_CAT,                                                &

      ! Aerosol emission fields,including mineral dust parent soil props
     &  JSOIL_CLAY, JSOIL_SILT, JSOIL_SAND,JDUST_MREL1,JDUST_MREL2,     &
     &  JDUST_MREL3,JDUST_MREL4,JDUST_MREL5,JDUST_MREL6,                &
     &  JSO2_EM, JDMS_EM, JSO2_HILEM, JNH3_EM, JSOOT_EM, JSOOT_HILEM,   &
     &  JBMASS_EM, JBMASS_HILEM, JOCFF_EM, JOCFF_HILEM, JDMS_CONC,      &
     &  JDMS_OFLUX,                                                     &
      ! User ancillary fields
     &  JUSER_ANC1,  JUSER_ANC2, JUSER_ANC3, JUSER_ANC4,                &
     &  JUSER_ANC5,  JUSER_ANC6, JUSER_ANC7, JUSER_ANC8,                &
     &  JUSER_ANC9,  JUSER_ANC10, JUSER_ANC11, JUSER_ANC12,             &
     &  JUSER_ANC13,  JUSER_ANC14, JUSER_ANC15, JUSER_ANC16,            &
     &  JUSER_ANC17,  JUSER_ANC18, JUSER_ANC19, JUSER_ANC20,            &
      ! Tracer fluxes - kdcorbin, 04/10
     &  JTRACER_FLUX1, JTRACER_FLUX2, JTRACER_FLUX3, JTRACER_FLUX4,     &
     &  JTRACER_FLUX5, JTRACER_FLUX6, JTRACER_FLUX7, JTRACER_FLUX8,     &
     &  JTRACER_FLUX9, JTRACER_FLUX10,JTRACER_FLUX11,JTRACER_FLUX12,    &
     &  JTRACER_FLUX13,JTRACER_FLUX14,JTRACER_FLUX15,JTRACER_FLUX16,    &
     &  JTRACER_FLUX17,JTRACER_FLUX18,JTRACER_FLUX19,JTRACER_FLUX20,    &
     &  JTSLAB,JTCLIM,JHCLIM,JCHEAT,JOIFLX,JUICE,JVICE,                 &
     &  JSIG11NE,JSIG11SE,JSIG11SW,JSIG11NW,                            &
     &  JSIG12NE,JSIG12SE,JSIG12SW,JSIG12NW,                            &
     &  JSIG22NE,JSIG22SE,JSIG22SW,JSIG22NW,                            &
      ! Pointers for ATMOSPHERE model constants. Scalars only.
     &  JETATHETA, JETARHO, JRHCRIT, JSOIL_THICKNESS,                   &
     &  Jzseak_theta,Jck_theta,Jzseak_rho,Jck_rho,                      &
      ! pointers for input variable grid info.
     &  JLAMBDA_INPUT_P,  JLAMBDA_INPUT_U,                              &
     &  JPHI_INPUT_P, JPHI_INPUT_V,                                     &
      ! pointers for tiled vegetation and triffid
     &  JFRAC_TYP, JLAI_PFT, JCANHT_PFT, JDISTURB, JSOIL_ALB,           &
     &  JFRAC_CON1, JFRAC_CON2, JFRAC_CON3, JFRAC_CON4, JFRAC_CON5,     &
     &  JFRAC_CON6, JFRAC_CON7, JFRAC_CON8, JFRAC_CON9,                 &
     &  JSOIL_CARB,                                                     &
     &  JSOIL_CARB1, JSOIL_CARB2, JSOIL_CARB3, JSOIL_CARB4,             &
     &  JNPP_PFT_ACC, JG_LF_PFT_ACC, JG_PHLF_PFT_ACC,                   &
     &  JRSP_W_PFT_ACC, JRSP_S_ACC,                                     &
     &  JRSP_S_ACC1, JRSP_S_ACC2, JRSP_S_ACC3,                          &
     &  JRSP_S_ACC4, JCAN_WATER_TILE, JCATCH_TILE,                      &
     &  JINFIL_TILE, JRGRAIN_TILE, JSNODEP_TILE, JTSTAR_TILE, JZ0_TILE, &
     &  JDOLR, JLW_DOWN, JSW_TILE,JNET_FLUX,JNET_MFLUX,                 &
      ! pointers for carbon cycle
     &  J_CO2FLUX,J_CO2_EMITS,                                          &
      ! pointers for large-scale hydrology
     &  JFEXP, JTI_MEAN, JTI_SIG, JGAMMA_INT,                           &
     &  JWATER_TABLE, JFSFC_SAT, JF_WETLAND, JSTHZW,                    &
     &  JA_FSAT,      JC_FSAT,   JA_FWET,    JC_FWET,                   &
      ! pointers for river routing
     &  JRIV_SEQUENCE, JRIV_DIRECTION, JRIV_STORAGE,                    &
     &  JTOT_SURFROFF, JTOT_SUBROFF, JRIV_INLANDATM                     &
      ! pointers for coupling fields
     & , JC_SOLAR, JC_BLUE, JC_DOWN, JC_LONGWAVE, JC_TAUX               &
     & , JC_TAUY, JC_WINDMIX, JC_SENSIBLE, JC_SUBLIM                    &
     & , JC_EVAP, JC_BOTMELTN, JC_TOPMELTN, JC_LSRAIN                   & 
     & , JC_LSSNOW, JC_CVRAIN, JC_CVSNOW, JC_RIVEROUT,JC_PRESS

! TYPPTRA end
#endif
