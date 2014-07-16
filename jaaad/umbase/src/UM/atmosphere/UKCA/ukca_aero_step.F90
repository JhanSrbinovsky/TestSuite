#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject 
! to the terms and conditions set out therein.
! [Met Office Ref SC138] 
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!  Does one chemistry time step of the UKCA-MODE aerosol model.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
!  Current Code Owner:       Graham Mann
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
! Subroutine Interface:
      SUBROUTINE UKCA_AERO_STEP(NBOX,                                   &
       ND,MDT,MD,MDWAT,S0G,DRYDP,WETDP,RHOPAR,DVOL,WVOL,SM,             &
       AIRD,AIRDM3,RHOA,MFPA,DVISC,T,TSQRT,RH,S,PMID,PUPPER,PLOWER,     &
       EMC,EMCBM,ZO3,ZHO2,ZH2O2,USTR,US10M,ZNOT,                        &
       SURTP,LAND_FRAC,SURF,SEAICE,                                     &
       CRAIN,DRAIN,CRAIN_UP,DRAIN_UP,FCONV_CONV,LOWCLOUD,VFAC,          &
       EMANSO2,EMVOLCONSO2,NEMVOLCONSO2,EMVOLEXPSO2,NEMVOLEXPSO2,       &
       EMBIOMSO2,ISO2EMS,                                               &
       DTC,DTM,DTZ,NMTS,NZTS,LDAY,ACT,BUD_AER_MAS,                      &
       PRIMSU_ON,PRIMBCOC_ON,PRIMSS_ON,PRIMDU_ON,RAINOUT_ON,            &
       IMSCAV_ON,WETOX_ON,DDEPAER_ON,SEDI_ON,                           &
       DRYOX_IN_AER,WETOX_IN_AER,DELSO2,DELSO2_2,                       &
       COND_ON,NUCL_ON,COAG_ON,ICOAG,IMERGE,IFUCHS,IWVOLMETHOD,         &
       IACTMETHOD,IDDEPAER,INUCSCAV,VERBOSE,CHECKMD_ND,JLAT,INTRAOFF,   &
       INTEROFF,IDUSTEMS,S0G_DOT_CONDENSABLE,LWC,JLABOVE,N_MERGE_1D)

!-----------------------------------------------------------------------
!  Inputs
!  ------
!  NBOX        : Number of grid boxes
!  ND          : Aerosol ptcl number density for mode (cm^-3)
!  MDT         : Avg tot mass of aerosol ptcl in mode (particle^-1)
!  MD          : Avg cpt mass of aerosol ptcl in mode (particle^-1)
!  MDWAT       : Molecular concentration of water (molecules per ptcl)
!  S0G         : Partial masses of gas phase species (kg per gridbox)
!  DRYDP       : Geometric mean dry diameter for each mode (m)
!  WETDP       : Geometric mean wet diameter for each mode (m)
!  RHOPAR      : Particle density [incl. H2O & insoluble cpts] (kgm^-3)
!  DVOL        : Geometric mean dry volume for each mode (m^3)
!  WVOL        : Geometric mean wet volume for each mode (m^3)
!  SM          : Grid box mass of air (kg)
!  AIRD        : Number density of air (per cm3)
!  AIRDM3      : Number density of air (per m3)
!  RHOA        : Air density (kg/m3)
!  MFPA        : Mean free path of air (m)
!  DVISC       : Dynamic viscosity of air (kg m^-1 s^-1)
!  T           : Centre level temperature (K)
!  TSQRT       : Square-root of centre level temperature (K)
!  RH          : Relative humidity (dimensionless 0-1)
!  S           : Specific humidity (kg/kg)
!  PMID        : Centre level pressure (Pa)
!  PUPPER      : Upper interface pressure (Pa)
!  PLOWER      : Lower interface pressure (Pa)
!  EMC         : BC/OC ems rates from bio- & fossil-fuels (kgC/box/s)
!  EMCBM       : BC/OC ems rates from biomass burning (kgC/box/s)
!  ZO3         : Backgrnd vmr of O3
!  ZHO2        : Backgrnd conc. of HO2 (molecules per cc)
!  ZH2O2       : Backgrnd conc. of H2O2 (molecules per cc)
!  USTR        : Surface friction velocity (m/s)
!  US10M       : Scalar wind at 10m (ms-1)
!  ZNOT        : Roughness length (m)
!  SURTP       : Surface type: 0=seasurf,1=landsurf,2=oversea,3=overland
!  LAND_FRAC   : Fraction of horizontal gridbox area covered by land
!  SURF        : Surface area of box (horizontal) (m^2)
!  SEAICE      : Fraction of horizontal gridbox area containing seaice
!  CRAIN       : Rain rate for conv precip. in box (kgm^-2s^-1)
!  DRAIN       : Rain rate for dyn. precip. in box (kgm^-2s^-1)
!  CRAIN_UP    : Rain rate for conv precip. in box above (kgm^-2s^-1)
!  DRAIN_UP    : Rain rate for dyn. precip. in box above (kgm^-2s^-1)
!  FCONV_CONV  : Fraction of box condensate --> rain in 6 hours (conv)
!  LOWCLOUD    : Horizontal low cloud fraction
!  VFAC        : Vertical low cloud fraction
!  EMANSO2     : Anthrop. SO2 ems rates, low sources (kgSO2/box/s)
!  EMVOLCONSO2(NBOX,10) Volcanic SO2 ems rates (cont. src) (kgSO2/box/s)
!  NEMVOLCONSO2(NBOX)   # of 1x1 cont. volc. SO2 ems sources in gridbox
!  EMVOLEXPSO2(NBOX,10) Volcanic SO2 ems rates (expl. src) (kgSO2/box/s)
!  NEMVOLEXPSO2(NBOX)   # of 1x1 expl. volc. SO2 ems sources in gridbox
!  EMBIOMSO2(NBOX,6) Biomass SO2 ems rates (kgSO2/box/s) [6 alt ranges]
!  ISO2EMS     : Switch for scheme for primary H2SO4 aerosol emissions
!  DTC         : Chemistry time step (s)
!  DTM         : Microphysics time step (s)
!  DTZ         : Competition (cond/nucl) time step (s)
!  NMTS        : Number of microphysics timesteps per DTC
!  NZTS        : Number of competition timesteps per DTM
!  LDAY        : 1-Day, 0-Night
!  MODESOL     : Specifies whether mode is soluble/insoluble (=1/0)
!  ACT         : Particle dry radius above which activation is assumed
!  PRIMSU_ON   : Switch : primary SO4 ptcl emissions are on/off (1/0)
!  PRIMBCOC_ON : Switch : primary BC/OC ptcl emissions are on/off (1/0)
!  PRIMSS_ON   : Switch : primary seasalt ptcl ems are on/off (1/0)
!  PRIMDU_ON   : Switch : primary dust particle ems are on/off (1/0)
!  RAINOUT_ON  : Switch : rainout (nucl. scav.) is on/off (1/0)
!  IMSCAV_ON   : Switch : impaction scavenging is on/off (1/0)
!  WETOX_ON    : Switch : wet oxidation [+ cloudproc] is on/off (1/0)
!  DDEPAER_ON  : Switch : aerosol dry deposition is on/off (1/0)
!  SEDI_ON     : Switch : aerosol sedimentation is on/off (1/0)
!  DRYOX_IN_AER: Switch : update gas phase condensible concentrations
!                         on competition tstep in aerosol code? (1/0)
!  WETOX_IN_AER: Switch : calc. wet ox. of SO2 in aerosol code? (1/0)
!  COND_ON     : Switch : vapour condensation is  on/off (1/0)
!  NUCL_ON     : Switch : binary nucleation is on/off (1/0)
!  COAG_ON     : Switch : coagulation is on/off (1/0)
!  ICOAG       : KIJ method (1:GLOMAP, 2: M7, 3: UMorig, 4:UMorig MFPP)
!  IMERGE      : Switch to use mid-pts (=1), edges (2) or dynamic (=3)
!  IFUCHS      : Switch : Fuchs (1964) or Fuchs-Sutugin (1971) for CC
!  IWVOLMETHOD : Switch : wet volume method (1=as-H2SO4,2=multi-cpt)
!  IACTMETHOD  : Switch : activation method (0=off,1=fixed ract,2=NSO3)
!  IDDEPAER    : Switch : dry dep method (1=as in Spr05, 2=incl. sedi)
!  INUCSCAV    : Switch : scheme for removal by nucl scav
!                (1=as GLOMAP [accsol,corsol], 2=scav params as Stier05]
!  VERBOSE     : If =1 prints min/max (ND,MDT etc) after each process
!                for 1st grid box (for box model tests)
!  INTRAOFF    : Switch to turn off intra-modal coagulation
!  INTEROFF    : Switch to turn off intra-modal coagulation
!  CHECKMD_ND  : Switch : check for values of MD, ND out of range? (1/0)
!  JLAT        : Index of latitude slice currently being processed
!  IDUSTEMS    : Switch : dust emissions scheme
!              :          (1=Pringle scheme, 2=AEROCOM00 daily)
!  S0G_DOT_CONDENSABLE : Gas phase chemistry tendencies for change in
!                        condensable gas phase species (vmr per s)
!                        due to chem, dry and wet deposition.
!  LWC         : Cloud liquid water content [kg/m3]
!  JLABOVE     : Index of box directly above this grid box
!  N_MERGE_1D  : Count number of mode-merges in each box, each mode
!  DELSO2      : S(IV) --> S(VI) by H2O2 (molecules per cc)
!                [input from gas phase chem. scheme if WETOX_IN_AER=0]
!  DELSO2_2    : S(IV) --> S(VI) by O3u  (molecules per cc)
!                [input from gas phase chem. scheme if WETOX_IN_AER=0]
!
!  Outputs
!  -------
!  ND          : Aerosol ptcl number density for mode (cm^-3)
!  MDT         : Avg tot mass of aerosol ptcl in mode (particle^-1)
!  MD          : Avg cpt mass of aerosol ptcl in mode (particle^-1)
!  MDWAT       : molecular concentration of water (molecules per ptcl)
!  S0G         : Partial masses of gas phase species (kg per gridbox)
!  BUD_AER_MAS: Aerosol budget terms for mass
!  DELSO2      : S(IV) oxidised to S(VI) by H2O2 (molecules per cc)
!                [input as zero if WETOX_IN_AER=1]
!  DELSO2_2    : S(IV) oxidised to S(VI) by O3 (molecules per cc)
!                [input as zero if WETOX_IN_AER=1]
!
!  Local Variables
!  ---------------
!  SO2         : Sulfur dioxide conc (molecules per cc)
!  SO2_VMR     : Sulfur dioxode vmr
!  H2SO4       : Sulfuric acid vapour conc (moleculesH2SO4/cc)
!  H2O2        : Hygrogen peroxide conc (per cc)
!  GC          : Condensable vapour conc (molecules cpt per cc)
!  GCOLD       : Condensable vapour conc b4 process (molecules cpt/cm3)
!  DELH2O2     : Change in hygrogen peroxide conc due to process (/cm3)
!  FRAC_AQ_ACC : Fraction of wet oxidised SO2 -> soluble accum. mode
!  FRAC_AQ_COR : Fraction of wet oxidised SO2 -> soluble coarse mode
!  TOTWETOX    : Total wet oxidation rate of SO2 (molecules/cm3)
!  DELTAGC_COND: Change in vapour conc due to cond (molecules cpt/cm3)
!  DELTAGC_NUCL: Change in vapour conc due to nucl (molecules cpt/cm3)
!  DELTAH2SO4_NUCL: Change in H2SO4 conc due to nucl (molecules/cm3)
!  DELTAS0G    : Overall change in partial masses of condensable gases
!                after all NZTS competition steps (kg/box).
!  S0G_TO_GC   : Molar mass ratio : condensing gas to aerosol phase cpt
!                e.g. MM_GAS(MSEC_ORG)=0.150 kg/mol but it may condense
!                into organic carbon component with MM(CP_OC)=0.0168.
!                Need to apply conversion when going from gas-->aerosol
!  JRATE       : H2SO4 depletion rate by BH nucleation (molecules/cm3/s)
!  AGETERM1    : Depletion rate of each component (molecules cpt/cc/DTZ)
!                from condensation onto the 3 insoluble modes.
!                Used for calculation of ageing rate in UKCA_AGEING.
!  AGETERM2    : Rate of accomodation of material to each insoluble mode
!                as a result of coagulation with smaller soluble modes
!                (in molecules cpt /cm3/DTZ)
!                Used for calculation of ageing rate in UKCA_AGEING.
!  PROCESS     : Character string to store process currently being done
!                for error checking in box model
!  KII_ARR     : Coag coeff for intra-modal coag (IMODE-IMODE) (cm^3/s)
!  KIJ_ARR     : Coag coeff for inter-modal coag (IMODE-JMODE) (cm^3/s)
!
!  Inputted by module UKCA_CONSTANTS
!  ---------------------------------
!  PPI         : 3.1415927...........
!  AVC         : Avogadros constant (mol-1)
!  ZBOLTZ      : Boltzmann's constant (kg m2 s-2 K-1 molec-1)
!  VKARMN      : Von Karman's constant = 0.4
!  RA          : Dry air gas constant = 287.05 Jkg^-1 K^-1
!  RR          : Universal gas constant = 8.314 J/mol/K
!  GG          : Gravitational acceleration = 9.80665 ms^-2
!  NMOL        : Number of molecules per particle at nucleation
!  MM_DA       : Molar mass of dry air (kg/mol)
!  CONC_EPS    : Threshold for molecular conc. (molecules per cc)
!  EMS_EPS     : Threshold for emissions fluxes (kg/gridbox/s)
!  DN_EPS      : Value of DELN below which do not carry out process
!
!  Inputted by module UKCA_MODE_SETUP
!  ----------------------------------
!  NMODES      : Number of possible aerosol modes
!  NCP         : Number of possible aerosol components
!  MODE        : Logical variable defining which modes are set.
!  COMPONENT   : Logical variable defining which cpt are in which dsts
!  SOLUBLE     : Logical variable defining which cpts are soluble
!  MM          : Molar masses of components (kg/mole)
!  RHOCOMP     : Densities (dry) of each component (kg/m^3)
!  NO_IONS     : Number of dissociating ions in soluble components
!  DDPLIM0     : Lower limit for dry diameter in mode (m)
!  DDPLIM1     : Upper limit for dry diameter in mode (m)
!  DDPMID      : Mid-point of size mode = exp(0.5*(lndp0+lndp1)) (m)
!  MFRAC_0     : Initial mass fraction to set when no particles.
!  SIGMA       : Geometric standard deviation for each mode
!  X           : EXP((9/2)*LOG^2(SIGMA_G))
!  NUM_EPS     : Value of NEWN below which do not carry out process
!  MMID        : Mass of particle with dry diameter
!                dp=dpmed_g=exp(0.5*(lndp0+lndp1)) (ptcl^-1)
!  COAG_MODE   : Switch to defines which mode an IMODE-JMODE
!                coagulation goes into
!  COLLEFF4    : Array of aerosol-raindrop collision efficiencies
!                for impaction scavenging (l-u table) [NCOLL,NROW]
!  RADDROP     : Raindrop radii (microns) at mid-points in raindrop grid
!                of dimension NROW -- for impaction scavenging
!  NCOLL       : # of columns in COLLEFF4 (=20) corresp to aerosol grid
!  NROW        : # of rows in COLLEFF4 (=19) corresp to raindrop grid
!  CP_SU       : Component where sulfate is stored
!  CP_BC       : Component where black carbon is stored
!  CP_OC       : Component where organic carbon is stored
!  CP_CL       : Component where NaCl is stored
!  CP_DU       : Component where dust is stored
!  CP_SO       : Component where condensible organic species is stored
!
!  Inputted by module UKCA_SETUP_INDICES
!  -------------------------------------
!  ICHEM       : Switch for gas phase chemistry scheme (0=off,1=on)
!  NADVG       : # of advected gas phase tracers
!  NCHEMG      : # of gas phase tracers for gas phase chemistry scheme
!  MM_GAS      : Molar masses for gas phase species (kg/mol)
!  MSOTWO      : Index of MM_GAS, WTRATC and S0G for SO2
!  MH2SO4      : Index of MM_GAS, WTRATC and S0G for H2SO4
!  MH2O2       : Index of MM_GAS, WTRATC and S0G for H2O2
!  MSEC_ORG    : Index of MM_GAS, WTRATC and S0G for SEC_ORG
!  CONDENSABLE : Logical variable defining which cpts are condensable
!  DIMEN       : Molecular diamters of condensable components (m)
!  Various indices for budget terms in BUD_AER_MAS
!
!-----------------------------------------------------------------------
      USE UKCA_CONSTANTS
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES

      IMPLICIT NONE
!
!     Inputs
      INTEGER :: NBOX
      INTEGER :: IDUSTEMS
      INTEGER :: NEMVOLCONSO2(NBOX)
      INTEGER :: NEMVOLEXPSO2(NBOX)
      INTEGER :: ISO2EMS
      INTEGER :: NMTS
      INTEGER :: NZTS
      INTEGER :: LDAY(NBOX)
      INTEGER :: PRIMSU_ON
      INTEGER :: PRIMBCOC_ON
      INTEGER :: PRIMSS_ON
      INTEGER :: PRIMDU_ON
      INTEGER :: RAINOUT_ON
      INTEGER :: IMSCAV_ON
      INTEGER :: WETOX_ON
      INTEGER :: DDEPAER_ON
      INTEGER :: SEDI_ON
      INTEGER :: DRYOX_IN_AER
      INTEGER :: WETOX_IN_AER
      INTEGER :: COND_ON
      INTEGER :: NUCL_ON
      INTEGER :: COAG_ON
      INTEGER :: ICOAG
      INTEGER :: IMERGE
      INTEGER :: IFUCHS
      INTEGER :: IWVOLMETHOD
      INTEGER :: IACTMETHOD
      INTEGER :: IDDEPAER
      INTEGER :: INUCSCAV
      INTEGER :: VERBOSE
      INTEGER :: CHECKMD_ND
      INTEGER :: JLAT
      INTEGER :: INTRAOFF
      INTEGER :: INTEROFF
      INTEGER :: JLABOVE(NBOX)
      INTEGER :: N_MERGE_1D(NBOX,NMODES)

      REAL :: ND(NBOX,NMODES)
      REAL :: MDT(NBOX,NMODES)
      REAL :: MD(NBOX,NMODES,NCP)
      REAL :: S0G(NBOX,NADVG)
      REAL :: MDWAT(NBOX,NMODES)
      REAL :: DRYDP(NBOX,NMODES)
      REAL :: WETDP(NBOX,NMODES)
      REAL :: RHOPAR(NBOX,NMODES)
      REAL :: DVOL(NBOX,NMODES)
      REAL :: WVOL(NBOX,NMODES)
      REAL :: SM(NBOX)
      REAL :: AIRD(NBOX)
      REAL :: AIRDM3(NBOX)
      REAL :: RHOA(NBOX)
      REAL :: MFPA(NBOX)
      REAL :: DVISC(NBOX)
      REAL :: T(NBOX)
      REAL :: TSQRT(NBOX)
      REAL :: RH(NBOX)
      REAL :: S(NBOX)
      REAL :: PMID(NBOX)
      REAL :: PUPPER(NBOX)
      REAL :: PLOWER(NBOX)
      REAL :: ZO3(NBOX)
      REAL :: ZHO2(NBOX)
      REAL :: ZH2O2(NBOX)
      REAL :: USTR(NBOX)
      REAL :: US10M(NBOX)
      REAL :: ZNOT(NBOX)
      REAL :: SURTP(NBOX)
      REAL :: LAND_FRAC(NBOX)
      REAL :: SURF(NBOX)
      REAL :: SEAICE(NBOX)
      REAL :: CRAIN(NBOX)
      REAL :: DRAIN(NBOX)
      REAL :: CRAIN_UP(NBOX)
      REAL :: DRAIN_UP(NBOX)
      REAL :: FCONV_CONV(NBOX)
      REAL :: LOWCLOUD(NBOX)
      REAL :: VFAC(NBOX)
      REAL :: EMANSO2(NBOX,6)
      REAL :: EMVOLCONSO2(NBOX,10)
      REAL :: EMVOLEXPSO2(NBOX,10)
      REAL :: EMBIOMSO2(NBOX,6)
      REAL :: EMC(NBOX,4)
      REAL :: EMCBM(NBOX,2,6)
      REAL :: DTC
      REAL :: DTM
      REAL :: DTZ
      REAL :: ACT
      REAL :: S0G_DOT_CONDENSABLE(NBOX,NCHEMG)
      REAL :: DELSO2(NBOX)
      REAL :: DELSO2_2(NBOX)
      REAL :: LWC(NBOX)
!
!     Outputs for budget calculations
      REAL :: BUD_AER_MAS(NBOX,NBUDAER)
!
!     Local variables
      INTEGER :: JV
      INTEGER :: ICP
      INTEGER :: IMTS
      INTEGER :: IZTS
      REAL :: SO2(NBOX)
      REAL :: H2O2(NBOX)
      REAL :: H2SO4(NBOX)
      REAL :: DELH2O2(NBOX)
      REAL :: SO2_VMR(NBOX)
      REAL :: GC(NBOX,NCHEMG)
      REAL :: GCOLD(NBOX,NCHEMG)
      REAL :: DELTAGC_COND(NBOX,NCHEMG)
      REAL :: DELTAGC_NUCL(NBOX,NCHEMG)
      REAL :: DELTAS0G(NBOX)
      REAL :: DELTAH2SO4_NUCL(NBOX)
      REAL :: JRATE(NBOX)
      REAL :: AGETERM1(NBOX,3,NCHEMG)
      REAL :: AGETERM2(NBOX,4,3,NCP)
      REAL :: S0G_TO_GC
      REAL :: FRAC_AQ_ACC(NBOX)
      REAL :: FRAC_AQ_COR(NBOX)
      REAL :: TOTWETOX(NBOX)
      REAL :: KII_ARR(NBOX,NMODES)
      REAL :: KIJ_ARR(NBOX,NMODES,NMODES)
      CHARACTER*30 :: PROCESS
      CHARACTER*2  :: STRIZTS
      LOGICAL :: MASK(NBOX)
      LOGICAL :: MASK1(NBOX),MASK2(NBOX),MASK3(NBOX)
!
!     CHECKMD_ND set to 1/0 for whether/not to check for
!                values of ND, MD and MDT out of range.
      IF(CHECKMD_ND == 1) THEN
        PROCESS='At start of UKCA_AERO_STEP    '
! DEPENDS ON: ukca_check_md_nd
        CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,JLAT,ND,MD,MDT)
      ENDIF
!
!     VERBOSE included so that max,min,mean of ND/MD can be
!             written out and monitored after each process
!             when error checking in box model.
      IF(VERBOSE >= 2) THEN
       write(6,*) 'At start of UKCA_AERO_STEP :   ND,MD,MDT'
! DEPENDS ON: ukca_calcminmaxndmdt
       CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT)
      ENDIF
!
!     Calculate primary sulfate aerosol emissions
      IF(ICHEM == 1) THEN
       IF(PRIMSU_ON == 1) THEN
!
! DEPENDS ON: ukca_prim_su
        CALL UKCA_PRIM_SU(NBOX,ND,MDT,MD,                               &
       EMANSO2,EMVOLCONSO2,NEMVOLCONSO2,                                &
       EMVOLEXPSO2,NEMVOLEXPSO2,EMBIOMSO2,                              &
       DTC,SM,ISO2EMS,AIRD,BUD_AER_MAS)
!
        IF(CHECKMD_ND == 1) THEN
         PROCESS='Done primary sulfate emissions'
! DEPENDS ON: ukca_check_md_nd
         CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,JLAT,ND,MD,MDT)
        ENDIF
!
        IF(VERBOSE >= 2) THEN
         write(6,*) 'After UKCA_PRIM_SU has updated ND,MD,MDT'
! DEPENDS ON: ukca_calcminmaxndmdt
         CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT)
        ENDIF
!
       ENDIF
      ENDIF
!
!     Calculate primary carbonaceous aerosol emissions
      IF(PRIMBCOC_ON == 1) THEN
!
! DEPENDS ON: ukca_prim_car
        CALL UKCA_PRIM_CAR(NBOX,ND,MDT,MD,                              &
       EMC,EMCBM,DTC,SM,AIRD,SURTP,ISO2EMS,BUD_AER_MAS)
!
        IF(CHECKMD_ND == 1) THEN
          PROCESS='Done primary bc/oc emissions  '
! DEPENDS ON: ukca_check_md_nd
          CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,JLAT,ND,MD,MDT)
        ENDIF
!
        IF(VERBOSE >= 2) THEN
         write(6,*) 'After UKCA_PRIM_CAR has updated ND,MD,MDT'
! DEPENDS ON: ukca_calcminmaxndmdt
         CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT)
        ENDIF
!
      ENDIF
!
!     Calculate primary sea salt aerosol emissions
      IF(PRIMSS_ON == 1) THEN
!
! DEPENDS ON: ukca_prim_ss
       CALL UKCA_PRIM_SS(NBOX,ND,MDT,MD,                                &
         SURTP,LAND_FRAC,SURF,SEAICE,DTC,SM,US10M,ZNOT,AIRD,            &
         BUD_AER_MAS,VERBOSE)
!
       IF(CHECKMD_ND == 1) THEN
        PROCESS='Done primary seasalt emissions'
! DEPENDS ON: ukca_check_md_nd
        CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,JLAT,ND,MD,MDT)
       ENDIF
!
       IF(VERBOSE >= 2) THEN
        write(6,*) 'After UKCA_PRIM_SS has updated ND,MD,MDT'
! DEPENDS ON: ukca_calcminmaxndmdt
        CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT)
       ENDIF
!
      ENDIF
!
!     Calculate primary dust aerosol emissions
!      IF(PRIMDU_ON == 1) THEN
!
!! DEPENDS ON: ukca_prim_du
!        CALL UKCA_PRIM_DU
!
!        IF(CHECKMD_ND == 1) THEN
!          PROCESS='Done primary dust emissions   '
!! DEPENDS ON: ukca_check_md_nd
!          CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,JLAT,ND,MD,MDT)
!        ENDIF
!
!        IF(VERBOSE >= 2) THEN
!         write(6,*) 'After UKCA_PRIM_DU has updated ND,MD,MDT'
!! DEPENDS ON: ukca_calcminmaxndmdt
!         CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT)
!        ENDIF
!
!      ELSE
!       MASPRIMDUACCINS(:)=0.0
!       MASPRIMDUCORINS(:)=0.0
!      ENDIF
!
!     Recalculate dry and wet diameter and volume
! DEPENDS ON: ukca_calc_drydiam
      CALL UKCA_CALC_DRYDIAM(NBOX,ND,MD,MDT,DRYDP,DVOL,VERBOSE)
!
! DEPENDS ON: ukca_volume_mode
      CALL UKCA_VOLUME_MODE(NBOX,ND,MD,MDT,                             &
        RH,PMID,T,WVOL,WETDP,RHOPAR,IWVOLMETHOD,                        &
        DVOL,DRYDP,MDWAT,VERBOSE)
!
!     Calculate impaction scavenging of aerosol (washout)
      IF(IMSCAV_ON == 1) THEN
!
! DEPENDS ON: ukca_impc_scav
        CALL UKCA_IMPC_SCAV(NBOX,ND,MD,                                 &
       CRAIN,DRAIN,WETDP,DTC,BUD_AER_MAS)
!
        IF(CHECKMD_ND == 1) THEN
         PROCESS='Done impaction scavenging     '
! DEPENDS ON: ukca_check_md_nd
         CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,JLAT,ND,MD,MDT)
        ENDIF
!
        IF(VERBOSE >= 2) THEN
         write(6,*) 'After UKCA_IMPC_SCAV hasupdatedND'
! DEPENDS ON: ukca_calcminmaxndmdt
         CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT)
        ENDIF
!
      ENDIF
!
!     Calculate in cloud, Nucleation Scavenging
      IF(RAINOUT_ON == 1) THEN
!
! DEPENDS ON: ukca_rainout
         CALL UKCA_RAINOUT(NBOX,ND,MD,FCONV_CONV,                       &
       CRAIN,DRAIN,CRAIN_UP,DRAIN_UP,T,DTC,                             &
       BUD_AER_MAS,INUCSCAV)
!
        IF(CHECKMD_ND == 1) THEN
          PROCESS='Done nucleation scavenging    '
! DEPENDS ON: ukca_check_md_nd
          CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,JLAT,ND,MD,MDT)
        ENDIF
!
        IF(VERBOSE >= 2) THEN
         write(6,*) 'After UKCA_RAINOUT has updated ND'
! DEPENDS ON: ukca_calcminmaxndmdt
         CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT)
        ENDIF
!
      ENDIF
!
      IF(ICHEM == 1) THEN
!
!      Calculate aqueous chemistry
       IF(WETOX_ON == 1) THEN
!
!       Initialise variables for aqueous chemistry
        SO2     (:)=0.0
        H2O2    (:)=0.0
        DELH2O2 (:)=0.0
!
        IF(WETOX_IN_AER == 1) THEN
!
         MASK(:)=(S0G(:,MSOTWO) > 0.0)
         WHERE(MASK(:))
          SO2(:)=S0G(:,MSOTWO)*AIRD(:)/SM(:)
         ENDWHERE
         SO2_VMR(:)=SO2(:)/AIRD(:)
!
         MASK(:)=(S0G(:,MH2O2) > 0.0)
         WHERE(MASK(:))
          H2O2(:)=S0G(:,MH2O2)*AIRD(:)/SM(:)
         ENDWHERE
!
! .. below calculates DELSO2, DELSO2_2, DELH2O2 in the
! .. aerosol module (WETOX_IN_AER=1)
!
! DEPENDS ON: ukca_wetox
         CALL UKCA_WETOX(NBOX,ND,DELSO2,DELSO2_2,                       &
           DELH2O2,LOWCLOUD,VFAC,SO2_VMR,H2O2,ZH2O2,ZO3,                &
           ZHO2,PMID,T,AIRD,S,DTC,LDAY,LWC)

!        Update SO2 partial mass after reaction with H2O2.
         S0G(:,MSOTWO)=S0G(:,MSOTWO)-DELSO2(:)*SM(:)/AIRD(:)
!        Update SO2 partial mass after reaction with O3.
         S0G(:,MSOTWO)=S0G(:,MSOTWO)-DELSO2_2(:)*SM(:)/AIRD(:)
         MASK(:)=(S0G(:,MSOTWO) < 0.0)
         WHERE(MASK(:))
          S0G(:,MSOTWO)=0.0
         ENDWHERE

!        Update H2O2 partial mass after SO2 wetox & replenishment
!                                      (both included in DELH2O2)
         S0G(:,MH2O2)=S0G(:,MH2O2)+DELH2O2(:)*SM(:)/AIRD(:)
!
        ENDIF ! IF WETOX_IN_AER=1
! n.b. if WETOX_IN_AER.NE.1, then DELSO2,DELSO2_2 passed in to AERO_STEP
!
       ELSE
        DELSO2  (:)=0.0
        DELSO2_2(:)=0.0
       ENDIF ! if WETOX_ON=1
!
       TOTWETOX(:)=DELSO2(:)+DELSO2_2(:) ! this is zero if WETOX_ON=0
!
       IF(IACTMETHOD > 0) THEN ! if cloud processing on
!
! .. below cloud-processes those aerosol in the Aitken soluble
! .. mode which are larger than ACT to accumulation mode
! .. representing the fact that they have activated so
! .. making minimum between Aitsol and accsol respond to activation
!
! DEPENDS ON: ukca_cloudproc
         CALL UKCA_CLOUDPROC(NBOX,ND,MD,MDT,DRYDP,                      &
       LOWCLOUD,VFAC,ACT,VERBOSE,IACTMETHOD,BUD_AER_MAS)
!
       ENDIF ! IF IACTMETHOD > 0
!
       MASK1(:)=((ND(:,3)+ND(:,4)) > NUM_EPS(4))
       MASK2(:)=(MASK1(:).AND.(ND(:,3) > NUM_EPS(3)))
       MASK3(:)=(MASK1(:).AND.(ND(:,4) > NUM_EPS(4)))
!
       WHERE(MASK1(:))
!
!       Calculate number fractions to partition between acc and cor
        FRAC_AQ_ACC(:)=ND(:,3)/(ND(:,3)+ND(:,4))
        FRAC_AQ_COR(:)=ND(:,4)/(ND(:,3)+ND(:,4))
!
       ENDWHERE
!
       WHERE(MASK2(:))
!       Update soluble accum mode H2SO4 mass due to aqueous chem.
        MD(:,3,CP_SU)=(MD(:,3,CP_SU)*ND(:,3)+                           &
                       FRAC_AQ_ACC(:)*TOTWETOX(:))/ND(:,3)
        MDT(:,3)=(MDT(:,3)*ND(:,3)+FRAC_AQ_ACC(:)*TOTWETOX(:))/ND(:,3)
       ENDWHERE ! if some particles in soluble accum mode
!
       IF(NMASCLPRSUACCSOL1 > 0) THEN
        WHERE(MASK2(:))
         BUD_AER_MAS(:,NMASCLPRSUACCSOL1)=                              &
           FRAC_AQ_ACC(:)*DELSO2  (:)*SM(:)/AIRD(:)
        ENDWHERE
       ENDIF
!
       IF(NMASCLPRSUACCSOL2 > 0) THEN
        WHERE(MASK2(:))
         BUD_AER_MAS(:,NMASCLPRSUACCSOL2)=                              &
           FRAC_AQ_ACC(:)*DELSO2_2(:)*SM(:)/AIRD(:)
        ENDWHERE
       ENDIF
!
       WHERE(MASK3(:))
!       Update soluble coarse mode H2SO4 mass due to aqueous chem.
        MD(:,4,CP_SU)=(MD(:,4,CP_SU)*ND(:,4)+                           &
                       FRAC_AQ_COR(:)*TOTWETOX(:))/ND(:,4)
        MDT(:,4)=(MDT(:,4)*ND(:,4)+FRAC_AQ_COR(:)*TOTWETOX(:))/ND(:,4)
       ENDWHERE ! if some particles in soluble coarse mode
!
       IF(NMASCLPRSUCORSOL1 > 0) THEN
        WHERE(MASK3(:))
         BUD_AER_MAS(:,NMASCLPRSUCORSOL1)=                              &
           FRAC_AQ_COR(:)*DELSO2  (:)*SM(:)/AIRD(:)
        ENDWHERE
       ENDIF
!
       IF(NMASCLPRSUCORSOL2 > 0) THEN
        WHERE(MASK3(:))
         BUD_AER_MAS(:,NMASCLPRSUCORSOL2)=                              &
           FRAC_AQ_COR(:)*DELSO2_2(:)*SM(:)/AIRD(:)
        ENDWHERE
       ENDIF
!
       IF(CHECKMD_ND == 1) THEN
        PROCESS='Done aqueous phase chemistry  '
! DEPENDS ON: ukca_check_md_nd
        CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,JLAT,ND,MD,MDT)
       ENDIF
!
       IF(VERBOSE >= 2) THEN
        write(6,*) 'After UKCA_WETOX has updated MD,MDT'
! DEPENDS ON: ukca_calcminmaxndmdt
        CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT)
       ENDIF
!
      ENDIF ! ICHEM=1
!
!     Calculate aerosol dry deposition
      IF(DDEPAER_ON == 1) THEN
       IF(IDDEPAER == 1) THEN
! DEPENDS ON: ukca_ddepaer
        CALL UKCA_DDEPAER(NBOX,ND,MD,MDT,                               &
          RHOPAR,ZNOT,SEAICE,                                           &
          DTC,WETDP,USTR,PMID,PUPPER,PLOWER,T,SURTP,                    &
          RHOA,MFPA,DVISC,BUD_AER_MAS)
       ENDIF
       IF(IDDEPAER == 2) THEN
! DEPENDS ON: ukca_ddepaer_incl_sedi
        CALL UKCA_DDEPAER_INCL_SEDI(NBOX,ND,MD,MDT,RHOPAR,ZNOT,         &
          DTC,WETDP,USTR,PMID,PUPPER,PLOWER,T,SURTP,SEAICE,             &
          RHOA,MFPA,DVISC,BUD_AER_MAS,JLABOVE,SEDI_ON,SM)
       ENDIF
!
       IF(CHECKMD_ND == 1) THEN
        PROCESS='Done aerosol dry deposition   '
! DEPENDS ON: ukca_check_md_nd
        CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,JLAT,ND,MD,MDT)
       ENDIF
!
       IF(VERBOSE >= 2) THEN
        write(6,*) 'After DDEPAER has updated ND'
! DEPENDS ON: ukca_calcminmaxndmdt
        CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT)
       ENDIF
!
      ENDIF ! DDEPAER_ON=1
!
!     Loop over number of microphysics timesteps NMTS
      DO IMTS=1,NMTS
!
!     Calculate rate of condensation of gases onto existing particles
!     and nucleation rate
!
       DO JV=1,NCHEMG ! loop over gas phase components
        IF(CONDENSABLE(JV)) THEN
!
!        Set index of cpt into which condensable gas to be stored
         ICP=CONDENSABLE_CHOICE(JV)
!
!        Calculate ratio of gas phase to aerosol component molar masses
         S0G_TO_GC=MM_GAS(JV)/MM(ICP)
!
         MASK(:)=(S0G(:,JV) > 0.0)
!
         WHERE(MASK(:))
!
!         GC is gas phase molecular concentration of condensables
!         but molecules refer to molar mass of aerosol cpt
          GC(:,JV)=S0G_TO_GC*S0G(:,JV)*AIRD(:)/SM(:)
!
         ENDWHERE
!
         WHERE(.NOT.MASK(:))
!
          GC(:,JV)=0.0
!
         ENDWHERE
!
         GCOLD(:,JV)=GC(:,JV)
!
        ENDIF
       ENDDO
!
!      Split DTM substep into NZTS subsubsteps to allow for competition
!      between nucleation and condensation (and gas phase production
!      if DRYOX_IN_AER=1) to compete on short timesteps.
!
! DEPENDS ON: ukca_calc_coag_kernel
       CALL UKCA_CALC_COAG_KERNEL(NBOX,KII_ARR,KIJ_ARR,                 &
        DRYDP,DVOL,WETDP,WVOL,RHOPAR,MFPA,DVISC,T,PMID,                 &
        COAG_ON,ICOAG)
!
       DO IZTS=1,NZTS
!
        IF(VERBOSE >= 2) THEN
         WRITE(STRIZTS,'(1i2.2)') IZTS
         PROCESS='Start of loop, IZTS='//STRIZTS//' (GC)   '
! DEPENDS ON: ukca_calcminmaxgc
         CALL UKCA_CALCMINMAXGC(PROCESS,NBOX,GC   )
         PROCESS='Start of loop, IZTS='//STRIZTS//' (GCOLD)'
! DEPENDS ON: ukca_calcminmaxgc
         CALL UKCA_CALCMINMAXGC(PROCESS,NBOX,GCOLD)
        ENDIF
!
        IF(ICHEM == 1) THEN
!
         IF(DRYOX_IN_AER > 0) THEN
!
          DO JV=1,NCHEMG
           IF(CONDENSABLE(JV)) THEN
!
!           Set index of cpt into which condensable gas will be stored
            ICP=CONDENSABLE_CHOICE(JV)
!
!           Calculate ratio of gas phase to aerosol cpt molar masses
            S0G_TO_GC=MM_GAS(JV)/MM(ICP)
!
            GC(:,JV)=GC(:,JV)+DTZ*S0G_DOT_CONDENSABLE(:,JV)             &
                             *AIRD(:)*S0G_TO_GC
!
! .. update condensable gas phase species by chemical tendencies on
! .. competition tstep
!
! .. n.b. S0G_DOT_CONDENSABLE is in vmr of S0G(:,JV) per s,
! .. so need to * by S0G_TO_GC
!
           ENDIF
          ENDDO
!
          IF(VERBOSE >= 2) THEN
           WRITE(STRIZTS,'(1i2.2)') IZTS
           PROCESS='Done chem upd, IZTS='//STRIZTS//' (GC)   '
! DEPENDS ON: ukca_calcminmaxgc
           CALL UKCA_CALCMINMAXGC(PROCESS,NBOX,GC   )
           PROCESS='Done chem upd, IZTS='//STRIZTS//' (GCOLD)'
! DEPENDS ON: ukca_calcminmaxgc
           CALL UKCA_CALCMINMAXGC(PROCESS,NBOX,GCOLD)
          ENDIF
!
         ENDIF ! DRYOX_IN_AER
!
!        Carry out uptake of condensable gas phase species onto aerosol
         IF(COND_ON == 1) THEN
! DEPENDS ON: ukca_conden
          CALL UKCA_CONDEN(NBOX,GC,ND,MD,MDT,                           &
       DTZ,DRYDP,WETDP,TSQRT,RHOA,AIRDM3,DELTAGC_COND,                  &
       IFUCHS,AGETERM1,BUD_AER_MAS)
!
          IF(VERBOSE >= 2) THEN
           WRITE(STRIZTS,'(1i2.2)') IZTS
           PROCESS='Done cond upd, IZTS='//STRIZTS//' (GC)   '
! DEPENDS ON: ukca_calcminmaxgc
           CALL UKCA_CALCMINMAXGC(PROCESS,NBOX,GC   )
           PROCESS='Done cond upd, IZTS='//STRIZTS//' (GCOLD)'
! DEPENDS ON: ukca_calcminmaxgc
           CALL UKCA_CALCMINMAXGC(PROCESS,NBOX,GCOLD)
          ENDIF
!
          IF(CHECKMD_ND == 1) THEN
           PROCESS='Done conden of H2SO4 & SEC_ORG'
! DEPENDS ON: ukca_check_md_nd
           CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,JLAT,ND,MD,MDT)
          ENDIF
!
          IF(VERBOSE >= 2) THEN
           write(6,*) 'After UKCA_CONDEN has updated MD,MDT'
! DEPENDS ON: ukca_calcminmaxndmdt
           CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT)
          ENDIF
!
         ENDIF
!
!        Calculate rate of binary H2SO4-H2O nucleation
         IF(NUCL_ON == 1) THEN
!
!         Set index of aerosol cpt which condensed H2SO4 to be stored
          ICP=CONDENSABLE_CHOICE(MH2SO4)

!         Calculate ratio of gas phase to aerosol cpt molar mass (H2SO4)
          S0G_TO_GC=MM_GAS(MH2SO4)/MM(ICP)
!
          MASK(:)=(GC(:,MH2SO4) > 0.0)
!
          WHERE(MASK(:))
          H2SO4(:)=GC(:,MH2SO4)/S0G_TO_GC
         ENDWHERE
! .. Set (competition-step-updated) H2SO4 molecular concentration
! .. (using MM_GAS) for calculation of BHN rate in UKCA_CALCNUCRATE
!
          WHERE(.NOT.MASK(:))
          H2SO4(:)=0.0
         ENDWHERE
!
! DEPENDS ON: ukca_calcnucrate
          CALL UKCA_CALCNUCRATE(NBOX,DTZ,T,S,RH,AIRD,H2SO4,             &
                DELTAH2SO4_NUCL,JRATE)
!
          DELTAGC_NUCL(:,MH2SO4)=DELTAH2SO4_NUCL(:)*S0G_TO_GC
!
          IF(NMASNUCLSUNUCSOL > 0)                                      &
           BUD_AER_MAS(:,NMASNUCLSUNUCSOL)=                             &
           BUD_AER_MAS(:,NMASNUCLSUNUCSOL)+DELTAGC_NUCL(:,MH2SO4)
!
          IF(VERBOSE >= 2) THEN
           WRITE(STRIZTS,'(1i2.2)') IZTS
           PROCESS='Done nucl upd, IZTS='//STRIZTS//' (GC)   '
! DEPENDS ON: ukca_calcminmaxgc
           CALL UKCA_CALCMINMAXGC(PROCESS,NBOX,GC   )
           PROCESS='Done nucl upd, IZTS='//STRIZTS//' (GCOLD)'
! DEPENDS ON: ukca_calcminmaxgc
           CALL UKCA_CALCMINMAXGC(PROCESS,NBOX,GCOLD)
          ENDIF
!
         ELSE
          DELTAGC_NUCL(:,MH2SO4)=0.0
          JRATE(:)=0.0
         ENDIF
        ELSE
         DELTAGC_NUCL(:,MH2SO4)=0.0
         JRATE(:)=0.0
        ENDIF
!
        IF(VERBOSE >= 2) THEN
         write(6,*) 'After DRYDIAM & VOLhas updated ND'
! DEPENDS ON: ukca_calcminmaxndmdt
         CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT)
        ENDIF
!
        IF((COAG_ON == 1).OR.(NUCL_ON == 1)) THEN
         IF(VERBOSE >= 2) THEN
          write(6,*) 'About to call UKCA_COAGWITHNUCL  '
! DEPENDS ON: ukca_calcminmaxndmdt
          CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT)
         ENDIF
!
!        Update ND & MD due to combined coagulation-nucleation
! DEPENDS ON: ukca_coagwithnucl
         CALL UKCA_COAGWITHNUCL(NBOX,ND,MD,MDT,DELTAGC_NUCL,DTZ,JRATE,  &
       AGETERM2,INTRAOFF,INTEROFF,VERBOSE,BUD_AER_MAS,KII_ARR,KIJ_ARR)
!
         IF(CHECKMD_ND == 1) THEN
          PROCESS='Done combined coag./nucleation'
! DEPENDS ON: ukca_check_md_nd
          CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,JLAT,ND,MD,MDT)
         ENDIF
!
         IF(VERBOSE >= 2) THEN
          IF((COAG_ON == 1).AND.(NUCL_ON == 1)) THEN
           write(6,*) 'After COAG & NUCL have updated ND,MD,MDT'
          ENDIF
          IF((COAG_ON == 0).AND.(NUCL_ON == 1)) THEN
           write(6,*) 'After NUCL has updated ND,MD,MDT'
          ENDIF
          IF((COAG_ON == 1).AND.(NUCL_ON == 0)) THEN
           write(6,*) 'After COAG has updated ND,MD,MDT'
          ENDIF
! DEPENDS ON: ukca_calcminmaxndmdt
          CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT)
         ENDIF
        ENDIF ! if COAG_ON = 1 or NUCL_ON = 1
!
!       Apply ageing -- transfer ND,MD from insol. to sol. modes
! DEPENDS ON: ukca_ageing
        CALL UKCA_AGEING(NBOX,ND,MD,MDT,                                &
          AGETERM1,AGETERM2,WETDP,VERBOSE,BUD_AER_MAS)
!
        IF(CHECKMD_ND == 1) THEN
         PROCESS='Done ageing of insoluble modes'
! DEPENDS ON: ukca_check_md_nd
         CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,JLAT,ND,MD,MDT)
        ENDIF
!
        IF(VERBOSE >= 2) THEN
         write(6,*) 'After UKCA_AGEING has updated ND,MD,MDT'
! DEPENDS ON: ukca_calcminmaxndmdt
         CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT)
        ENDIF
!
       ENDDO ! end loop over competition subsubtimesteps
!
!      Recalculate dry and wet diameter and volume
! DEPENDS ON: ukca_calc_drydiam
       CALL UKCA_CALC_DRYDIAM(NBOX,ND,MD,MDT,DRYDP,DVOL,VERBOSE)
!
! DEPENDS ON: ukca_volume_mode
       CALL UKCA_VOLUME_MODE(NBOX,ND,MD,MDT,                            &
         RH,PMID,T,WVOL,WETDP,RHOPAR,IWVOLMETHOD,                       &
         DVOL,DRYDP,MDWAT,VERBOSE)
!
!      Apply mode-merging where necessary
! DEPENDS ON: ukca_remode
       CALL UKCA_REMODE(NBOX,ND,MD,MDT,DRYDP,WETDP,VERBOSE,             &
         IMERGE,BUD_AER_MAS,N_MERGE_1D)
!
       IF(CHECKMD_ND == 1) THEN
        PROCESS='Done UKCA_REMODE             '
! DEPENDS ON: ukca_check_md_nd
        CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,JLAT,ND,MD,MDT)
       ENDIF
!
       IF(VERBOSE >= 2) THEN
        write(6,*) 'After UKCA_REMODE has updated ND,MD,MDT'
! DEPENDS ON: ukca_calcminmaxndmdt
        CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT)
       ENDIF
!
       IF(ICHEM == 1) THEN
!
        DO JV=1,NCHEMG
         IF(CONDENSABLE(JV)) THEN
!
          ICP=CONDENSABLE_CHOICE(JV)
!
!         Calculate ratio of gas phase to aerosol cpt molar masses
          S0G_TO_GC=MM_GAS(JV)/MM(ICP)
!
          DELTAS0G(:)=(GC(:,JV)-GCOLD(:,JV))*(SM(:)/AIRD(:))/S0G_TO_GC
!
          MASK(:)=(DELTAS0G(:) < -S0G(:,JV))
! above limits deltaS0G to be -S0G (stop -ves)
          WHERE(MASK(:))
          DELTAS0G(:)=-S0G(:,JV)
         ENDWHERE
!
          S0G(:,JV)=S0G(:,JV)+DELTAS0G(:)
!
         ENDIF
        ENDDO
!
       ENDIF ! if ICHEM = 1
!
      ENDDO ! end loop over NMTS

      RETURN
      END SUBROUTINE UKCA_AERO_STEP
#endif
