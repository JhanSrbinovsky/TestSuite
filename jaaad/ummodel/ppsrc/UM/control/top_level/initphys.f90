
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL--------------- SUBROUTINE INITPHYS --------------------------------
!LL
!LL Purpose : Calls the initialisation program for the data assimilation
!LL assimilation section (P3)
!LL  Control routine for CRAY YMP
!LL
!LL T.Johns     <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL 3.2    27/03/93 Change INITPHYS for dynamic allocation. Remove
!LL                 size parameters no longer needed.       R.Rawlins.
!   4.0    06/07/95 CNTLATM called to allow access to version numbers
!                   to pick up new routines for reading spectral
!                   namelists for version 3A of the SW or LW radiation.
!                                             J. M. Edwards
!   4.4    02/09/97 Error checking for different sections improved.
!                   Logical flags for aerosols passed down to allow
!                   earlier checking of the spectral data.
!                                             J. M. Edwards
!   4.5    18/05/98 Logicals for including gases are now passed down
!                   to lower levels to select gases from the
!                   spectral file.
!                                             J. M. Edwards
!  4.5  April 1998  Pass soot logical to R2_SW_SPECIN. Luke Robinson.
!LL 5.0    30/06/99 Add call to CRUNTIMC to get L_CLIMAT_AEROSOL, which
!LL                 used to be in cntlatm at 4.5. M.L.Gallani
!   5.2    15/11/00 Add flags for direct effect of sea-salt aerosol
!                   to calls to R2_SW_SPECIN and R2_LW_SPECIN.
!                                             A. Jones
!   5.3    04/04/01 Pass mesoscale aerosol logical L_MURK_RAD
!                   to lower levels.                      S. Cusack
!   5.5    05/02/03 Add logical for direct effect of biomass aerosol
!                   to calls to R2_SW_SPECIN and R2_LW_SPECIN.
!                                                P Davison
!   5.5    21/02/03 Add logical for direct effect of mineral dust
!                   to calls to R2_SW_SPECIN and R2_LW_SPECIN.
!                                                S Woodward
!   6.2    15/11/05 Add logical saying whether the aerosol optical
!                   depth block is needed to call to R2_LW_SPECIN.
!                                                N Bellouin
! 6.2  25/11/05 Functionality for improved time stepping, radiative
!               forcing and radiance code added for versions 3C
!               and 3Z of radiation code             (J.-C. Thelen)
!   6.2    21/03/06 Included nstypes.h J Ridley.
!LL
!LL Programming standard; U M Documentation Paper No. 3 version 1
!LL dated 15/01/90
!LL
!LL System components covered P0
!LL
!LL Documentation : UM documentation paper  no P0
!LL
!LL END


      SUBROUTINE INITPHYS(ICODE,CMESSAGE)

!
! If new version of radiation code is required then we
! need to use the modules for improved time-stepping
!
!
      IMPLICIT NONE

      INTEGER    ICODE         ! Return code : 0 Normal exit
!                              !             : >0 Error
      CHARACTER*80 CMESSAGE    ! Error message if ICODE > 0

!*
!*L Subroutines called
      EXTERNAL R2_SW_SPECIN, R2_LW_SPECIN
!*
!------------------------ nstypes.h ----------------------------------
!jhan:further renovation of ths file may be necessary params are dependent on dataset
!jhan: ALSO nstypes_cable.h should be unecessary nsoil/soil is only difference
      !--- Number of non-vegetation surface types
      Integer, Parameter :: NNVG  = 4

      !--- Number of plant functional types.
      Integer, Parameter :: NPFT  = 13
      
      !--- Number of surface types.
      Integer, Parameter :: NTYPE =17 
      
      !--- Index of the surface type 'Soil'
      !Integer, Parameter :: SOIL  = 16 
      !dhb599, 20110615: change made as per Peter Vohralik, item 1:
      Integer, Parameter :: SOIL  = 14

!--- Land surface types :
!--- original veg. tiles 
!     1 - Broadleaf Tree
!     2 - Needleleaf Tree
!     3 - C3 Grass
!     4 - C4 Grass
!     5 - Shrub
!--- for testing these tiles are set = 1:5 
!     6 - Broadleaf Tree
!     7 - Needleleaf Tree
!     8 - C3 Grass
!     9 - C4 Grass
!    10 - Shrub
!--- for testing these tiles are set = 0
!    11 - 0 
!    11 - 0
!    11 - 0
!--- original non-veg tiles moved to these indices
!     14 - Urban
!     15 - Water
!     16 - Soil
!     17 - Ice


! CMAXSIZE maximum sizes for dimensioning arrays
! of model constants whose sizes are configuration dependent. This
! allows constants to be read in from a NAMELIST file and maintain
! the flexibility of dynamic allocation for primary variables. The
! maximum sizes should agree with the maximum sizes implicit in the
! front-end User Interface.

!
!  Model            Modification history:
! version  Date
! 3.2  26/03/93  New COMDECK. Author R.Rawlins
! 3.4  06/08/94: Parameter MAX_NO_OF_SEGS used to dimension addresses
!                in macro-tasked calls to SWRAD, LWRAD & CONVECT.
!                Authors: A.Dickinson, D.Salmond, Reviewer: R.Barnes
! 3.5  22/05/95  Add MAX_N_INTF. D. Robinson
! 4.5  29/07/98  Increase MAX_N_INTF/MAX_N_INTF_A to 8. D. Robinson.
! 5.0  20/04/99  Changes for conversion to C-P C dynamics grid.
!                R. Rawlins
!  6.1   04/08/04  Add diffusion variable max_power     Terry Davies
! 6.2  25/12/05  Add max_updiff_levels/max_sponge_width   Terry Davies

      INTEGER,PARAMETER::max_model_levels = 100 ! Maximum no. of levels

      ! Max levels in boundary layer
      INTEGER,PARAMETER:: max_bl_levels = max_model_levels

      ! Max size of alpha_Cd
      INTEGER,PARAMETER :: max_number_alpha_cds = max_bl_levels

      ! Max no. of levels for pvort output
      INTEGER,PARAMETER :: MAX_REQ_THPV_LEVS = max_model_levels

      ! Max no. 1-2-1 rows in polar filter
      INTEGER,PARAMETER ::  max_121_rows =  8
      ! 0 is used for horizontal diffusion pointer

      ! Max no. of levels (from top) to apply upper level diffusion
      INTEGER,PARAMETER ::  max_updiff_levels = 10

      ! Max size of any sponge zones
      INTEGER,PARAMETER ::  max_sponge_width = 10

      ! Max size of look-up tables for searches
      INTEGER,PARAMETER ::  max_look = 2048

      ! Max no. of atmos interface areas
      INTEGER,PARAMETER :: MAX_N_INTF_A =  8

      ! Max no. of points in LBC      
      INTEGER,PARAMETER :: max_intf_lbcrow_length = 1000
      INTEGER,PARAMETER :: max_intf_lbcrows = 1000
        
      ! Max no. of atmos interface levels
      INTEGER,PARAMETER :: MAX_INTF_LEVELS = max_model_levels

      ! Maximum number of physics segments
      INTEGER,PARAMETER :: MAX_NO_OF_SEGS = 200
      ! MAX_N_INTF/MAX_N_INTF_A to be sorted out in later version
      ! Max no. of interface areas
      INTEGER, PARAMETER :: MAX_N_INTF =  8
! CMAXSIZE end
! ----------------------- Header file CRUNTIMC  -----------------------
! Description: Run-time constants for the Atmosphere model (read only).
!              Contains variables that define parametrization values
!              chosen for atmosphere physics and dynamics schemes.
!              [Note that CNTLATM holds accompanying control switches
!              needed for addressing.]
!
! Current Code Owner: R. Rawlins
!
!------------------   Physics:   --------------------------------------
! Generalised physics switches:


! Boundary layer:

      LOGICAL :: l_use_bl_diag_term
      LOGICAL :: L_SBLeq
      LOGICAL :: L_SBLco
      LOGICAL :: L_LAMBDAM2    ! LambdaM=2*LambdaH (operational setting)
      LOGICAL :: L_FULL_LAMBDAS! Lambdas NOT reduced above NTML_LOCAL+1
      LOGICAL :: L_us_blsol    ! Switch for stable and non-oscillatory
                               !   BL vertical diffusion scheme
      REAL :: alpha_Cd (max_number_alpha_cds)
      REAL :: Muw_SBL
      REAL :: Mwt_SBL
      REAL :: Charnock
      REAL :: Z0HSEA           ! Roughness lengths for heat and moisture
      REAL :: Z0MIZ            !   over sea, marginal ice zone and
      REAL :: Z0SICE           !   sea ice (m)
      REAL :: ALBSNC_NVG(NNVG) ! Snow-covered albedo.
      REAL :: ALBSNF_NVG(NNVG) ! Snow-free albedo.
      REAL :: CATCH_NVG(NNVG)  ! Canopy capacity (kg/m2).
      REAL :: GS_NVG(NNVG)     ! Surface conductance (m/s).
      REAL :: INFIL_NVG(NNVG)  ! Infiltration enhancement factor.
      REAL :: ROOTD_NVG(NNVG)  ! Rootdepth (m).
      REAL :: Z0_NVG(NNVG)     ! Roughness length (m).
      REAL :: CH_NVG(NNVG)     ! "Canopy" heat capacity (J/K/m2)
      REAL :: VF_NVG(NNVG)     ! Fractional "canopy" coverage
      REAL :: Puns, Pstb ! parameters for uncond stable numerical solver
                         ! Puns : used in an unstable BL column
                         ! Pstb : used in an stable BL column
      REAL :: OROG_DRAG_PARAM  ! Drag coefficient for orographic form drag
      LOGICAL :: L_EMIS_LAND_GEN ! Aggregate land surface emissivity
      REAL :: EMIS_LAND_GEN      ! Aggregate land surface emissivity
      REAL :: WeightLouisToLong_in
!                        ! Weighting of Louis tails towards long
!                        ! tails: this should be removed in favour
!                        ! of placing RUN_BL in the module bl_option_mod 
!                        ! at a later release.

      INTEGER, DIMENSION(20) :: BL_OPTIONS
      INTEGER :: DECFIX
      INTEGER :: STOPWE_SBL
      INTEGER :: TRWEIGHTS1
      INTEGER :: FLUX_GRAD
      INTEGER :: NG_STRESS
      INTEGER :: SBL_OP
      INTEGER :: ISHEAR_BL
      INTEGER :: FORMDRAG, COR_UST, COR_MO_ITER
      INTEGER :: NON_LOCAL_BL
      INTEGER :: LOCAL_FA
      INTEGER :: PRANDTL
      INTEGER :: I_SCRN_T_DIAG
      INTEGER :: Buddy_sea
      INTEGER :: FD_stab_dep
      INTEGER :: Keep_Ri_FA
      INTEGER :: NL_BL_LEVELS

      INTEGER :: CAN_RAD_MOD              ! For radiation model (Pete Falloon)
      INTEGER :: ILAYERS                  ! For radiation model (Pete Falloon)

      COMMON  /RUN_BLVEG/ALBSNC_NVG,ALBSNF_NVG,CATCH_NVG,GS_NVG,        &
     &  INFIL_NVG,ROOTD_NVG,Z0_NVG,CH_NVG,VF_NVG,                       &
     &  CAN_RAD_MOD,ILAYERS
       NAMELIST/RUN_BLVEG/ALBSNC_NVG,ALBSNF_NVG,CATCH_NVG,GS_NVG,       &
     &  INFIL_NVG,ROOTD_NVG,Z0_NVG,CH_NVG,VF_NVG,                       &
     &  CAN_RAD_MOD,ILAYERS

      COMMON  /RUN_BLICE/Z0HSEA,Z0MIZ,Z0SICE
      NAMELIST /RUN_BLICE/Z0HSEA,Z0MIZ,Z0SICE
      REAL :: SeaSalinityFactor
!       Scaling of qsat allowing for salinity of sea water
      INTEGER :: ISeaZ0T
!       Option for specifying the thermal roughness length
!       at sea points
      INTEGER :: ISeaDynDiag
!       Option for dynamic diagnosis of boundary layer types
       ! Maximum and minimum values for the STPH_RP scheme
       ! Boundary Layer
       REAL         :: par_mezcla_max,par_mezcla,par_mezcla_min ! Max, 
                      ! mean and min value for the neutral mixing length
       REAL         :: G0_max,G0_RP,G0_min ! Max,mean and min values
                      ! for the flux profile parameter
       REAL         :: Charnock_max,Charnock_min ! Max and min values for
                      ! the charnock parameter

      COMMON  /RUN_BL/ FORMDRAG, OROG_DRAG_PARAM, l_use_bl_diag_term         &
     &  ,L_SBLeq, L_SBLco, alpha_Cd, Muw_SBL, Mwt_SBL, BL_OPTIONS            &
     &  ,par_mezcla_max, par_mezcla, par_mezcla_min                          &
     &  ,G0_max, G0_RP, G0_min, Charnock_max, Charnock, Charnock_min         &
     &  ,L_LAMBDAM2, L_FULL_LAMBDAS, SeaSalinityFactor, ISeaZ0T              &
     &  ,ISeaDynDiag ,L_us_blsol, Puns, Pstb, L_EMIS_LAND_GEN, EMIS_LAND_GEN &
     &  ,WeightLouisToLong_in

      NAMELIST/RUN_BL/ FORMDRAG, OROG_DRAG_PARAM                             &
     &  ,l_use_bl_diag_term, SBL_OP, COR_UST,  COR_MO_ITER                   &
     &  ,L_SBLeq, L_SBLco, alpha_Cd, Muw_SBL, Mwt_SBL, ISHEAR_BL, NG_STRESS  &
     &  ,DECFIX, STOPWE_SBL, TRWEIGHTS1, FLUX_GRAD                           &
     &  ,NON_LOCAL_BL, NL_BL_LEVELS                                          &
     &  ,par_mezcla_max, par_mezcla, par_mezcla_min                          &
     &  ,G0_max, G0_RP, G0_min, Charnock_max, Charnock, Charnock_min         &
     &  ,L_LAMBDAM2, L_FULL_LAMBDAS, LOCAL_FA, Prandtl, SeaSalinityFactor    &
     &  ,ISeaZ0T, ISeaDynDiag, Buddy_sea, L_us_blsol, Puns, Pstb             &
     &  ,FD_stab_dep, Keep_Ri_FA                                             &
     &  ,L_EMIS_LAND_GEN, EMIS_LAND_GEN, I_SCRN_T_DIAG                       &
     &  ,WeightLouisToLong_in

! Surface parameters for plant functional types

      REAL    :: ALBSNC_MAX(NPFT)
      REAL    :: ALBSNC_MIN(NPFT)
      REAL    :: ALBSNF_MAX(NPFT)
      REAL    :: DZ0V_DH(NPFT)
      REAL    :: CATCH0(NPFT)
      REAL    :: DCATCH_DLAI(NPFT)
      REAL    :: INFIL_F(NPFT)
      REAL    :: KEXT(NPFT)
      REAL    :: ROOTD_FT(NPFT)

      COMMON  /RUN_PFT/ALBSNC_MAX,ALBSNC_MIN,ALBSNF_MAX,DZ0V_DH,        &
     &  CATCH0,DCATCH_DLAI,INFIL_F,KEXT,ROOTD_FT
      NAMELIST/RUN_PFT/ALBSNC_MAX,ALBSNC_MIN,ALBSNF_MAX,DZ0V_DH,        &
     &  CATCH0,DCATCH_DLAI,INFIL_F,KEXT,ROOTD_FT


      ! Large scale precipitation:

      ! Threshold cloud liquid water content over land/sea for
      ! conversion to ppn (kg water per m**3)

      ! Defunct: not supported
      REAL,PARAMETER:: CW_LAND = 8.0E-4
      REAL,PARAMETER:: CW_SEA  = 2.0E-4


      LOGICAL :: L_cry_agg_dep ! Limit deposition to moisture available
      LOGICAL :: L_it_melting  ! Use Iterative melting
      LOGICAL :: L_seq_mcr     ! Use sequential microphysics updating
      LOGICAL :: L_autoc_3b    ! Use 3B autoconversion method
      LOGICAL :: L_autolim_3b  ! Use 3B autoconversion limit
      LOGICAL :: L_autoconv_murk  ! Use murk aerosol to calculate
                                  ! droplet number
      LOGICAL :: L_droplet_settle ! Allow cloud droplets to settle
      LOGICAL :: L_psd         ! Use generic ice particle size distn.

      REAL :: ec_auto          ! Collision / collection efficiency
      REAL :: N_drop_land      ! Droplet concentration over land
      REAL :: N_drop_sea       ! Droplet concentration over sea
      REAL :: Ntot_land        ! Droplet concentration over land
      REAL :: Ntot_sea         ! Droplet concentration over sea
      REAL :: N_drop_land_cr   ! N_drop_land ^ (-1/3)
      REAL :: N_drop_sea_cr    ! N_drop_sea ^ (-1/3)
      REAL :: X1R              ! Intercept parameter for raindrop
                               ! size distribution
      REAL :: X2R              ! Scaling parameter for rain DSD.
      REAL :: X4R              ! Shape factor for rain DSD.
      REAL :: X1I              ! Intercept parameter for aggregate
                               ! size distribution
      REAL :: X1IC             ! Intercept parameter for crystal
                               ! size distribution
      REAL :: AI,  BI          ! Aggregate mass-size params m(D)=AI D^BI
      REAL :: AIC, BIC         ! Crystal mass-size params m(D)=AIC D^BIC
      REAL :: LSP_EI, LSP_FI   ! Aggregate Best-Reynolds relationship
                               ! Re = LSP_EI BE^LSP_FI
      REAL :: LSP_EIC, LSP_FIC ! Crystal Best-Reynolds relationship
                               ! Re = LSP_EIC BE^LSP_FIC

      ! Critical RH for layer cloud formation
      REAL :: RHCRIT(max_model_levels)
       ! Maximum and minimum values for the STPH_RP scheme
       ! Large Scale Precipitation
      REAL         :: RHCRIT_max,RHCRIT_min ! Max and min values
                       ! for the critical relative humidity
      REAL         :: CI_max,CI_min ! Max and min values for the
                       ! CI parameter, which controls the ice fall speed
      REAL         :: M_CI_max,M_CI,M_CI_min ! Max, mean, min values for
                       ! the multiplication factor for CI and CIC in 3C
                       ! Microphysics
 
      COMMON  /RUN_Precip/RHCRIT                                        &
     & ,l_cry_agg_dep,L_it_melting,L_seq_mcr,L_autoc_3b,L_autolim_3b,   &
     &  l_psd,                                                          &
     &  l_autoconv_murk,ec_auto,N_drop_land,N_drop_sea,Ntot_land,       &
     &  Ntot_sea,N_drop_land_cr,N_drop_sea_cr,X1R,X2R,X4R,X1I,X1IC      &
     &  ,AI,BI,AIC,BIC,LSP_EI,LSP_FI,LSP_EIC,LSP_FIC                    &
     &  ,RHCRIT_max,RHCRIT_min                                          &
     &  ,CI_max,CI_min,M_CI_max,M_CI,M_CI_min,L_droplet_settle
      NAMELIST/RUN_Precip/RHCRIT                                        &
     & ,l_cry_agg_dep,L_it_melting,L_seq_mcr,L_autoc_3b,L_autolim_3b,   &
     &  l_psd,                                                          &
     &  l_autoconv_murk,ec_auto,N_drop_land,N_drop_sea,Ntot_land,       &
     &  Ntot_sea,N_drop_land_cr,N_drop_sea_cr,X1R,X2R,X4R,X1I,X1IC      &
     &  ,AI,BI,AIC,BIC,LSP_EI,LSP_FI,LSP_EIC,LSP_FIC                    &
     &  ,RHCRIT_max,RHCRIT_min                                          &
     &  ,CI_max,CI_min,M_CI_max,M_CI,M_CI_min,L_droplet_settle

      ! Large scale cloud:

      LOGICAL :: L_eacf                ! Use empirically adjusted
                                       ! cloud fraction

      INTEGER :: cloud_fraction_method ! Selects total cloud fraction
                                       ! calculation method
      INTEGER :: ice_fraction_method   ! Selects ice cloud fraction
                                       ! calculation method

      REAL    :: overlap_ice_liquid    ! Generic overlap parameter
                                       ! between ice and liquid phases
      REAL    :: ctt_weight            ! Cloud top temperature weight
      REAL    :: t_weight              ! Local temperature weight
      REAL    :: qsat_fixed            ! Prescribed qsat value
      REAL    :: sub_cld               ! Scaling factor
      REAL    :: dbsdtbs_turb_0        ! PC2 erosion rate / s-1

      COMMON  /RUN_Cloud/L_eacf,cloud_fraction_method,                  &
     &  overlap_ice_liquid,ice_fraction_method,ctt_weight,t_weight,     &
     &  qsat_fixed,sub_cld,dbsdtbs_turb_0

      NAMELIST/RUN_Cloud/L_eacf,cloud_fraction_method,                  &
     &  overlap_ice_liquid,ice_fraction_method,ctt_weight,t_weight,     &
     &  qsat_fixed,sub_cld,dbsdtbs_turb_0


      ! Radiation:

      ! True if climatological aerosol is included.
      LOGICAL :: L_climat_aerosol


      ! True to use real boundary layer heights to specify the boundary
      ! layer aerosol.
      LOGICAL :: L_clim_aero_hgt
      LOGICAL :: L_HadGEM1_Clim_Aero
!             Flag to use HadGEM1 setting for climatological aerosols

      LOGICAL :: L_SEC_VAR ! true if using time varying astronomy
      LOGICAL :: L_EqT     ! True if including the equation of time

      INTEGER :: A_SW_SEGMENTS ! No of batches used in shortwave code
      INTEGER :: A_SW_SEG_SIZE ! Size of sw batches. 
      INTEGER :: A_LW_SEGMENTS ! No of batches used in longwave code
      INTEGER :: A_LW_SEG_SIZE ! Size of lw batches.
      INTEGER :: AERO_BL_LEVELS !Common number of layers taken to be 
!                                occupied by the boundary-layer
!                                aerosol if the boundary layer
!                                depth is not used to determine the 
!                                number separately at each grid-point
!                                In previous versions of the code,
!                                this was taken to be BL_LEVELS 

      Logical :: l_ovrlap    ! Requires l_ccrad=.TRUE.
                             ! Allows Convective and LS Cloud to overlap 
                             ! for radiative impacts.
                             ! (Experimental, defaulted to FALSE,
                             ! requires a hand-edit to change)

      Logical :: l_ccw_scav  ! Requires l_ccrad=.TRUE. .AND. l_ovrlap=.TRUE.
                             ! Allows Convective Cloud Water (CCW) to
                             ! compensate for LS Cloud water in overlapping
                             ! LS/CCA fractions.
                             ! (Experimental, defaulted to FALSE, requires a
                             ! hand-edit to change)

      Logical :: l_cldtop_t_fix
                              ! .TRUE. Corrects Cloud Top Temperature in
                              !        radiation

      ! Valid options are specified in the include file o3intp.h.

      ! radiation block is already a separate module for some variables:
!
! Description:
!   Common Block for radiation
!   Contains:
!   Parameters used to increase the true surface albedo linearly from
!   a minimum value, alphaM at melting point, to a "cold-ice" value,
!   alphaC, over a temperature range, dTice
!   Mass Mixing Ratios of minor Gases N2O,CH4,CFC11,CFC12,O2
!   CFC113,HCFC22,HFC125, and HFC134A
!
! Current Code Owner: S Woodward
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 3.4      09/09/94 Original code. C Hewitt
! 3.4      02/10/94 Taken over, and MMRs added  S Woodward
! 4.4      18/09/97 Mass Mixing Ratio of Oxygen added
! 4.5      21/08/98 ALPHAB introduced and ALPHAM,ALPHAC reinterpreted
!          if l_ssice_albedo i.e. if snow on sea-ice affects albedo.
!          Jonathan Gregory
! 4.5      18/05/98 Mass Mixing Ratios of further halocarbons added.
!                                    J. M. Edwards
! 5.5      11/02/03 Added IS_NCOL    K.D.Williams
! 6.2      06/01/06 sea ice albedo parameters added  J.Ridley
!
! Declarations:
! Imported global variables (*CALLed COMDECKs etc...)

! Global parameters:
! Albedo of sea-ice at melting point (TM) if .not.l_ssice_albedo, or
! Albedo of snow on sea-ice at melting point (TM) if l_ssice_albedo
      REAL ALPHAM  ! "M" for "melting"
! Albedo of sea-ice at and below TM-DTICE if .not.l_ssice_albedo, or
! Albedo of snow on sea-ice at and below TM-DTICE if l_ssice_albedo
      REAL ALPHAC  ! "C" for "cold"
! Albedo of snow-free sea-ice if l_ssice_albedo
      REAL ALPHAB  ! "B" for "bare"
! Temperature range in which albedo of sea-ice, if .not.l_ssice_albedo,
! or of snow on sea-ice, if l_ssice_albedo, varies between its limits
      REAL DTICE
! Temperature range below TM over which meltponds form if
! l_sice_meltponds and l_ssice_albedo
      REAL DT_BARE
! Increment to albedo for each degree temperature rises above
! TM-DT_BARE. Only used if l_sice_meltponds and l_ssice_albedo
      REAL DALB_BARE_WET
! Fraction of SW radiation that penetrates seaice and scatters
! back causing an addition to the albedo. Only active if l_ssice_albedo
! and l_sice_scattering
      REAL PEN_RAD_FRAC
! attenutation parameter for SW in seaice which controls the
! additional albedo due to internal scattering
      REAL SW_BETA
      REAL N2OMMR  ! N2O mmr (Mass Mixing Ratio)
      REAL CH4MMR  ! CH4 mmr
      REAL C11MMR  ! CFC11 mmr
      REAL C12MMR  ! CFC12 mmr
      REAL O2MMR   ! O2 mmr
      REAL C113MMR      ! CFC113 mmr
      REAL HCFC22MMR    ! HCFC22 mmr
      REAL HFC125MMR    ! HFC125 mmr
      REAL HFC134AMMR   ! HFC134A mmr
! Number of columns used for internal sampling by the ISCCP simulator
      INTEGER IS_NCOL

! Global scalars:

! Global dynamic arrays:

! COMMON blocks:
      COMMON /RAD_COM/ ALPHAM, ALPHAC, ALPHAB, DTICE, N2OMMR, CH4MMR,   &
     & C11MMR, C12MMR, O2MMR, C113MMR, HCFC22MMR, HFC125MMR, HFC134AMMR,&
     & IS_NCOL, DT_BARE, DALB_BARE_WET, PEN_RAD_FRAC, SW_BETA

!- End of COMDECK declaration

      REAL :: CO2_MMR                ! CO2 concentration (if constant)

! DIMFIX3A defines internal dimensions tied to algorithms for
! two-stream radiation code, mostly for clouds

      ! number of components of clouds
      INTEGER,PARAMETER:: NPD_CLOUD_COMPONENT=4

      ! number of permitted types of clouds.
      INTEGER,PARAMETER:: NPD_CLOUD_TYPE=4

      ! number of permitted representations of clouds.
      INTEGER,PARAMETER:: NPD_CLOUD_REPRESENTATION=4

      ! number of overlap coefficients for clouds
      INTEGER,PARAMETER:: NPD_OVERLAP_COEFF=18

      ! number of coefficients for two-stream sources
      INTEGER,PARAMETER:: NPD_SOURCE_COEFF=2

      ! number of regions in a layer
      INTEGER,PARAMETER:: NPD_REGION=3

! DIMFIX3A end
      ! Scaling factors to simulate inhomogeneous cloud.
      REAL, DIMENSION(NPD_CLOUD_COMPONENT) :: INHOM_CLOUD_SW
      REAL, DIMENSION(NPD_CLOUD_COMPONENT) :: INHOM_CLOUD_LW


! sza: add fraction standard diviation of inhomogeneous cloud
!      parameter for triple_clouds
!
      Real FW_STD

      ! Decorrelation pressure scale for large scale cloud
      REAL :: DP_CORR_STRAT
      ! Decorrelation pressure scale for convective cloud
      REAL :: DP_CORR_CONV

      COMMON  /RUN_Radiation/L_climat_aerosol, L_clim_aero_hgt,         &
     &  L_HadGEM1_Clim_Aero, FW_STD,                                    &
     &  A_SW_SEGMENTS,A_SW_SEG_SIZE,A_LW_SEGMENTS,A_LW_SEG_SIZE,CO2_MMR,&
     &  L_SEC_VAR,L_EqT,INHOM_CLOUD_SW,INHOM_CLOUD_LW,DP_CORR_STRAT,    &
     &  DP_CORR_CONV,AERO_BL_LEVELS, l_ovrlap, l_ccw_scav, l_cldtop_t_fix

      NAMELIST/RUN_Radiation/L_climat_aerosol, L_clim_aero_hgt,         &
     &  L_HadGEM1_Clim_Aero, FW_STD,                                    &
     &  A_SW_SEGMENTS,A_SW_SEG_SIZE,A_LW_SEGMENTS,A_LW_SEG_SIZE,CO2_MMR,&
     &  L_SEC_VAR,L_EqT,INHOM_CLOUD_SW,INHOM_CLOUD_LW,DP_CORR_STRAT,    &
     &  DP_CORR_CONV, AERO_BL_LEVELS,l_ovrlap,l_ccw_scav,l_cldtop_t_fix,&
      ! list of variables declared in <rad_com/rad_com.h> to be added:
     &  ALPHAM, ALPHAC, ALPHAB, DTICE,                                  &
     &  DT_BARE,DALB_BARE_WET,PEN_RAD_FRAC,SW_BETA,                     &
     &  N2OMMR, CH4MMR, C11MMR, C12MMR,                                 &
     &  O2MMR, C113MMR, HCFC22MMR, HFC125MMR, HFC134AMMR, IS_NCOL

      ! Gravity wave drag:

      REAL :: KAY_GWAVE       ! Surface stress constant for GW drag

      REAL :: GWD_FRC         ! Critical Froude Number for 4A Scheme
       ! Maximum and minimum values for the STPH_RP scheme
       ! Gravity Wave drag
       REAL         :: GWD_FRC_max, GWD_FRC_min ! Max and min values
                       ! for the critical froud number
       REAL         :: KAY_GWAVE_max, KAY_GWAVE_min ! Max and min values
                       ! for the gravity wave constant

      REAL :: GWD_FSAT       ! Critical Froude number for wave breaking
                             ! used in the amplitude based saturation
                             ! test

      INTEGER :: SAT_SCHEME  ! Switch to determine whether to use
                             ! stress based    (sat_scheme=0) or
                             ! amplitude based (sat_scheme=1)
                             ! saturation test

      LOGICAL :: L_TAUS_SCALE  ! FROUDE NO. DEPENDENT SURFACE STRESS
      LOGICAL :: L_FIX_GWSATN  ! BUG FIXES TO GWSATN
      LOGICAL :: L_GWD_40km    ! Turn off orographic GWD above 40km
      LOGICAL :: L_USSP_OPAQUE ! OPAQUE LID FOR USSP SCHEME

      COMMON  /RUN_GWD/KAY_GWAVE, GWD_FRC, GWD_FSAT                     &
                      ,L_TAUS_SCALE, L_FIX_GWSATN, l_gwd_40km           &
                      ,L_USSP_OPAQUE, SAT_SCHEME                        &
                      ,GWD_FRC_max, GWD_FRC_min                         &
                      ,KAY_GWAVE_max, KAY_GWAVE_min

      NAMELIST/RUN_GWD/KAY_GWAVE, GWD_FRC, GWD_FSAT                     &
                      ,L_TAUS_SCALE, L_FIX_GWSATN, l_gwd_40km           &
                      ,L_USSP_OPAQUE, SAT_SCHEME                        &
                      ,GWD_FRC_max, GWD_FRC_min                         &
                      ,KAY_GWAVE_max, KAY_GWAVE_min

!------------------   End of Physics   ---------------------------------

      ! River Routing
      REAL :: RIVER_VEL       ! River velocity
      REAL :: RIVER_MCOEF     ! Meandering coefficient

      COMMON /RUN_RIVERS/ RIVER_VEL, RIVER_MCOEF
      NAMELIST /RUN_RIVERS/ RIVER_VEL, RIVER_MCOEF


      ! Aerosol Modelling

      INTEGER :: SO2_HIGH_LEVEL   ! Model level for chimney SO2 emiss
      INTEGER :: SOOT_HIGH_LEVEL  ! Model level for chimney soot emiss
      INTEGER :: BMASS_HIGH_LEVEL_1  ! Lowest and highest model levels
      INTEGER :: BMASS_HIGH_LEVEL_2  ! with biomass emissions.
      INTEGER :: OCFF_HIGH_LEVEL  ! Model level for chimney OCFF emiss

      ! improved nonhydrostatic  weights in tracer1 calculation
      LOGICAL :: L_tracer1_non_hydro

      COMMON  /RUN_Aerosol/                                             &
     & SO2_HIGH_LEVEL, SOOT_HIGH_LEVEL, BMASS_HIGH_LEVEL_1,             &
     & BMASS_HIGH_LEVEL_2, OCFF_HIGH_LEVEL, L_tracer1_non_hydro
      NAMELIST/RUN_Aerosol/                                             &
     & SO2_HIGH_LEVEL, SOOT_HIGH_LEVEL, BMASS_HIGH_LEVEL_1,             &
     & BMASS_HIGH_LEVEL_2, OCFF_HIGH_LEVEL, L_tracer1_non_hydro

!     Declarations for UKCA sub-model

      LOGICAL :: L_ukca_chem      ! True when UKCA chemistry is on
      LOGICAL :: L_ukca_family    ! True when using family chemistry
      LOGICAL :: L_ukca_advh2o    ! True when advecting H2O tracer
      LOGICAL :: L_ukca_phot2d    ! True when using 2D photolysis
      LOGICAL :: L_ukca_fastj     ! True when using Fastj photolysis
      LOGICAL :: L_ukca_trop      ! True for tropospheric chemistry
      LOGICAL :: L_ukca_tropisop  ! True for trop chemistry + isoprene
      LOGICAL :: L_ukca_strat     ! True for strat+reduced trop chemistry
      LOGICAL :: L_ukca_strattrop ! True for std strat+trop chemistry
      LOGICAL :: L_ukca_stratcfc  ! True for extended strat chemistry
      LOGICAL :: L_ukca_wachem    ! True for whole atm chemistry
      LOGICAL :: L_ukca_aerchem   ! True for trop+aerosol chemistry
      LOGICAL :: L_ukca_user      ! True for user-defined chemistry
      LOGICAL :: L_ukca_mode      ! True for UKCA-MODE aerosol scheme
      LOGICAL :: L_ukca_dust      ! True for UKCA dust scheme
      LOGICAL :: L_ukca_rnpb      ! True for Radon/Lead scheme
      LOGICAL :: L_ukca_o3budget  ! True for UKCA o3 budget
      LOGICAL :: L_ukca_qch4inter ! True for interact wetland CH4 ems
      LOGICAL :: L_ukca_isopinter ! True for interact isoprene ems
      LOGICAL :: L_ukca_terpinter ! True for interact terpene ems
      LOGICAL :: L_ukca_budget2   ! True for UKCA budget2 calculations
      LOGICAL :: L_ukca_qf11f12mbr  ! True for CFC emissions
      LOGICAL :: L_ukca_useumuivals ! True when using UMUI CFC values
      LOGICAL :: L_ukca_nat_sedi  ! True for NAT sedimentation
      LOGICAL :: L_ukca_het_psc   ! True for Het/PSC chemistry
      LOGICAL :: L_ukca_h2o_feedback ! True for H2O feedback from chem
      LOGICAL :: L_ukca_clbrcons  ! True for Cl/Br conservation
      LOGICAL :: L_ukca_usero3    ! True when using RO3 in photolysis
      LOGICAL :: L_ukca_useco3    ! True when using clim O3 in photolysis
      LOGICAL :: L_ukca_userelaxo3 ! True when using clim O3
!                                    above 0.3hPa in photolysis
      LOGICAL :: L_ukca_rado3   ! T when using UKCA O3 in radiation
      LOGICAL :: L_ukca_radch4  ! T when using UKCA CH4 in radiation
      LOGICAL :: L_ukca_radn2o  ! T when using UKCA N2O in radiation
      LOGICAL :: L_ukca_radf11  ! T when using UKCA CFC-11 in radn
      LOGICAL :: L_ukca_radf12  ! T when using UKCA CFC-12 in radn
      LOGICAL :: L_ukca_radf113 ! T when using UKCA CFC-113 in radn
      LOGICAL :: L_ukca_radf22  ! T when using UKCA HCFC-22 in radn
      LOGICAL :: L_ukca_radch2o ! T when using UKCA adv H2O in radn
      LOGICAL :: L_ukca_useoxid ! T when using UKCA oxidants in aer chem
      LOGICAL :: L_ukca_intdd   ! T when using interact dry deposition

! Tracers and chemistry integers:
      Integer :: jpctr          ! No. of chemical tracers
      Integer :: jpspec         ! No. of chemical species
      Integer :: jpbk           ! No. of bimolecular reactions
      Integer :: jptk           ! No. of termolecular reactions
      Integer :: jppj           ! No. of photolytic reactions
      Integer :: jphk           ! No. of heterogonous reactions
      Integer :: jpnr           ! jpbk + jptk + jppj + jphk
      Integer :: jpdd           ! No. of dry deposited species
      Integer :: jpdw           ! No. of wet deposited species
      Integer :: jpeq           ! No. of species in aqueous phase

! UKCA_MODE control features:
      LOGICAL :: L_ukca_primsu   ! T for primary sulphate aerosol emissions
      LOGICAL :: L_ukca_primss   ! T for primary sea-salt aerosol emissions
      LOGICAL :: L_ukca_primbcoc ! T for primary BC and OC aerosol emissions
      LOGICAL :: L_ukca_sedi     ! T for aerosol sedimentation
      LOGICAL :: L_ukca_nucl     ! T for aerosol sedimentation

      INTEGER :: i_mode_setup    ! Defines MODE aerosol scheme
      INTEGER :: i_mode_sizeprim ! Defines MODE size parameters
      INTEGER :: i_mode_nucscav  ! Defines nucleation scavenging method
      INTEGER :: i_mode_ddepaer  ! Defines dry deposition method
      INTEGER :: i_mode_nzts     ! No. of substeps for nucleation/sedimentation

      REAL :: mode_actdryr       ! Activation dry radius for MODE

      COMMON  /RUN_UKCA/ L_ukca_chem, L_ukca_family, L_ukca_advh2o,     &
     &         L_ukca_phot2d, L_ukca_fastj, L_ukca_trop,                &
     &         L_ukca_tropisop, L_ukca_strat, L_ukca_strattrop,         &
     &         L_ukca_stratcfc, L_ukca_wachem, L_ukca_aerchem,          &
     &         L_ukca_user, L_ukca_mode, L_ukca_dust,                   &
     &         L_ukca_rnpb, L_ukca_o3budget, L_ukca_qch4inter,          &
     &         L_ukca_isopinter, L_ukca_terpinter, L_ukca_budget2,      &
     &         L_ukca_qf11f12mbr, L_ukca_useumuivals, L_ukca_nat_sedi,  &
     &         L_ukca_het_psc, L_ukca_h2o_feedback, L_ukca_clbrcons,    &
     &         L_ukca_usero3, L_ukca_useco3, L_ukca_userelaxo3,         &
     &         L_ukca_rado3, L_ukca_radch4, L_ukca_radn2o,              &
     &         L_ukca_radf11, L_ukca_radf12, L_ukca_radf113,            &
     &         L_ukca_radf22, L_ukca_radch2o, L_ukca_useoxid,           &
     &         L_ukca_intdd, L_ukca_primsu, L_ukca_primss,              &
     &         L_ukca_primbcoc, L_ukca_sedi, L_ukca_nucl,               &
     &         jpctr, jpspec, jpbk, jptk, jppj, jphk, jpnr,             &
     &         jpdd, jpdw, jpeq,                                        &
     &         i_mode_setup, i_mode_sizeprim, i_mode_nucscav,           &
     &         i_mode_ddepaer, i_mode_nzts, mode_actdryr

      NAMELIST/RUN_UKCA/ L_ukca_chem, L_ukca_family, L_ukca_advh2o,     &
     &         L_ukca_phot2d, L_ukca_fastj, L_ukca_trop,                &
     &         L_ukca_tropisop, L_ukca_strat, L_ukca_strattrop,         &
     &         L_ukca_stratcfc, L_ukca_wachem, L_ukca_aerchem,          &
     &         L_ukca_user, L_ukca_mode, L_ukca_dust,                   &
     &         L_ukca_rnpb, L_ukca_o3budget, L_ukca_qch4inter,          &
     &         L_ukca_isopinter, L_ukca_terpinter, L_ukca_budget2,      &
     &         L_ukca_qf11f12mbr, L_ukca_useumuivals, L_ukca_nat_sedi,  &
     &         L_ukca_het_psc, L_ukca_h2o_feedback, L_ukca_clbrcons,    &
     &         L_ukca_usero3, L_ukca_useco3, L_ukca_userelaxo3,         &
     &         L_ukca_rado3, L_ukca_radch4, L_ukca_radn2o,              &
     &         L_ukca_radf11, L_ukca_radf12, L_ukca_radf113,            &
     &         L_ukca_radf22, L_ukca_radch2o, L_ukca_useoxid,           &
     &         L_ukca_intdd, L_ukca_primsu, L_ukca_primss,              &
     &         L_ukca_primbcoc, L_ukca_sedi, L_ukca_nucl,               &
     &         jpctr, jpspec, jpbk, jptk, jppj, jphk, jpnr,             &
     &         jpdd, jpdw, jpeq,                                        &
     &         i_mode_setup, i_mode_sizeprim, i_mode_nucscav,           &
     &         i_mode_ddepaer, i_mode_nzts, mode_actdryr

      INTEGER :: Instability_diagnostics  ! >0 if wanted, 0 otherwise
      INTEGER :: print_step    ! To control diagnostic printing interval
      INTEGER :: diag_interval ! diagnostic printing sampling frequency
      INTEGER :: norm_lev_start ! start level for norm diagnostics
      INTEGER :: norm_lev_end   ! end level for norm diagnostics
      INTEGER :: first_norm_print ! first timestep for norm printing
      INTEGER :: rpemax ! array size needed for diagnostic printing
      INTEGER :: rpemin ! array size needed for diagnostic printing
      INTEGER :: rpesum ! array size needed for diagnostic printing
      INTEGER :: ipesum ! array size needed for diagnostic printing
      INTEGER :: time_theta1_min      ! Timestep of min level 1 theta
      INTEGER :: time_w_max(max_model_levels) ! Timestep of max w
      INTEGER :: time_div_max(max_model_levels) ! Timestep of max div
      INTEGER :: time_div_min(max_model_levels) ! Timestep of min div
      INTEGER :: time_lapse_min(max_model_levels) ! Timestep of min
      INTEGER :: time_max_shear(max_model_levels) !Timestep max shear
      INTEGER :: time_max_wind(max_model_levels) ! Timestep of max wind
      INTEGER :: time_KE_max(max_model_levels) ! Timestep of max KE
      INTEGER :: time_KE_min(max_model_levels) ! Timestep of min KE
      INTEGER :: time_noise_max(max_model_levels) ! Timestep of max

      REAL:: frictional_timescale(max_model_levels) ! For idealised case
      REAL :: w_print_limit ! w Threshold for diagnostic printing
      REAL :: w_conv_limit  ! w Threshold for limiting convection
      REAL :: tropics_deg  ! define latitude for tropics
      REAL :: min_theta1_run                  ! Min theta level 1
      REAL :: dtheta1_run   ! Largest -ve delta theta at min theta1
      REAL :: max_w_run(max_model_levels)  ! Max w at a level
      REAL :: max_div_run(max_model_levels) ! Max divergence at a level
      REAL :: min_div_run(max_model_levels) ! Min divergence at a level
      REAL :: min_lapse_run(max_model_levels) ! Min dtheta/dz at a level
      REAL :: max_shear_run(max_model_levels) ! Max shear at a level
      REAL :: max_wind_run(max_model_levels) ! Max wind at a level
      REAL :: max_KE_run(max_model_levels)   ! Max KE at a level
      REAL :: min_KE_run(max_model_levels)   ! Min KE at a level
      REAL :: max_noise_run(max_model_levels) ! Max noise at a level

!     Problem_number not set here  ! Now controlled by namelist input
!     Instability_diagnostics      ! Now controlled by namelist input
!     frictional_timescale         ! Now intitialised in SETCONA

      LOGICAL :: L_idealised_data
      LOGICAL :: L_simple_friction
      LOGICAL :: L_print_w
      LOGICAL :: L_print_div
      LOGICAL :: L_diag_print_ops     ! diagnostic prints for ops
      LOGICAL :: L_print_pe     ! print diagnostics on all pe's if true
      LOGICAL :: L_print_shear        ! wind diagnostic prints
      LOGICAL :: L_print_max_wind     ! wind diagnostic prints
      LOGICAL :: L_diag_L2norms  ! l2norm diagnostic prints
      LOGICAL :: L_diag_L2helm   ! l2norm diagnostic prints from solver
      LOGICAL :: L_diag_wind     ! wind diagnostic prints
      LOGICAL :: L_diag_noise    ! w diagnostic prints
      LOGICAL :: L_flush6        ! if T then flush buffers on failure
      LOGICAL :: L_inviscid  ! Use prognostic inviscid theta at surface
      LOGICAL :: L_deep        ! Deep/shallow switch
      LOGICAL :: L_hydrostatic ! Hydrostatic option
      LOGICAL :: L_geostrophic ! Geostrophic balance at start
      LOGICAL :: L_solver      ! Solver active if true
      LOGICAL :: L_trap_errors ! Error trapping
      LOGICAL :: L_uv_zero     ! Set winds to zero
      LOGICAL :: L_diag_print  ! Print diagnostics
      LOGICAL :: L_print_lapse ! Print lapse_rate diagnostics
      LOGICAL :: L_print_wmax ! Print max w diagnostic
      LOGICAL :: L_print_theta1 ! Print level 1 theta diagnostic

!     L_idealised_data             ! Now controlled by namelist input

!------------------   Diagnostics:   --------------------------------

      COMMON  /RUN_Diagnostics/                                         &
     &  L_diag_print, L_diag_print_ops, L_flush6, print_step,           &
     &  diag_interval,w_print_limit, rpemax, rpemin, ipesum, rpesum,    &
     &  L_print_w, L_print_wmax,L_print_lapse, L_print_theta1,          &
     &  L_print_div, L_diag_wind, L_print_shear, L_print_max_wind,      &
     &  L_print_pe, L_diag_noise, L_diag_L2norms, L_diag_L2helm,        &
     &  norm_lev_start, norm_lev_end, first_norm_print,                 &
     &  max_w_run, min_theta1_run, dtheta1_run,                         &
     &  max_div_run, min_div_run, min_lapse_run,                        &
     &  max_shear_run, max_wind_run,                                    &
     &  max_noise_run, max_KE_run, min_KE_run,                          &
     &  time_KE_max, time_KE_min,                                       &
     &  time_w_max, time_div_max, time_div_min, time_lapse_min,         &
     &  time_max_shear, time_max_wind,                                  &
     &  time_theta1_min, time_noise_max

!------------------   Dynamics:   --------------------------------------

      ! Generalised dynamics:

      LOGICAL :: L_Physics   ! T: physics to be included
      LOGICAL :: L_Run_With_Physics2   ! T: physics2 to be included
      LOGICAL :: L_perturb_IC_theta    ! T: perturb theta on ts1   
      LOGICAL :: L_Backwards ! F: Integrate backwards without physics.
      LOGICAL :: L_Primitive ! T: Primitive Equation model
                   ! F: Semi-Geostrophic model (no dynamics)
      LOGICAL :: L_dry         ! T: run with no moisture
      LOGICAL :: L_adjust_wet  ! T: perform simple dry&moist adjustment
      LOGICAL :: L_mix_ratio

      LOGICAL :: L_LBC_balance  ! T: impose vertically balanced Exners
                                !    and rho in LBCs (set w=0)
                                ! F: leave Exners alone

      ! Switch for new (correct) lbcs developed in 2008
      LOGICAL :: L_lbc_new       ! .true. for new lbc switch

      LOGICAL :: L_free_slip    ! .true. for free-slip lower boundary

      ! T: reset polar values to mean value every polar_reset_timesteps
      LOGICAL, PARAMETER :: L_polar_reset = .false.

      ! interval for resetting polar to mean
      INTEGER, PARAMETER :: polar_reset_timesteps = 1

      INTEGER :: IntRand_seed

      ! Time weights for integration scheme:
      REAL :: alpha_1
      REAL :: alpha_2
      REAL :: alpha_3
      REAL :: alpha_4

!  variable resolution control
      LOGICAL :: L_regular  ! true if NOT variable resolution
      LOGICAL :: L_qwaterload  ! true if using adding waterloading terms

      INTEGER :: lam_var    ! number of variable res. lambda intervals
      INTEGER :: phi_var    ! number of variable res. phi intervals

      REAL :: var_ratio  ! grid-stretcing ratio for var grid
      REAL :: lam_ratio  ! scaling of original grid to high res grid
      REAL :: phi_ratio  ! scaling of original grid to high res grid
      REAL :: lam_frac  ! proportion of reg. points in West of domain
      REAL :: phi_frac  ! proportion of reg. points in South of domain
      REAL :: lambda_p_end
      REAL :: phi_p_end
      REAL :: lambda_u_end
      REAL :: phi_v_end
      REAL :: dlambda_p_end
      REAL :: dphi_p_end
      REAL :: dlambda_u_end
      REAL :: dphi_v_end

      ! GCR scheme= Generalized conjugate residual elliptic equation
      ! solver:
      LOGICAL :: GCR_use_tol_abs
      LOGICAL :: GCR_zero_init_guess
      LOGICAL :: GCR_use_residual_Tol

      ! T then use full equation on RHS on second and subsequent ADI
      ! timesteps
      LOGICAL :: GCR_adi_add_full_soln

      LOGICAL :: L_gcr_fast_x ! use faster non reproducible code
! Speed up solver convergence when dynamics-physics cycling
! is used.
       Logical :: L_GCR_cycle_opt

      LOGICAL :: L_interp_depart ! interpolate to find u,v departure pts

! ----------------------- include file: GCR_ITER_DIM -------------------
! Description: Introduce integer parameters to:
!              (1) declare the maximum allowed number of dynamics 
!                  cycles (timestep iterations)
!              (2) declare the number of GCR Iterations analysis steps
!
! Code Author: M. Diamantakis
! Current Code Owner: A. Malcom
!
      INTEGER, PARAMETER :: MAX_NUMCYCLES = 10 ! Max number of 
                                               ! dynamics cycles
      INTEGER, PARAMETER :: GCR_ANAL_STEPS = 3 ! Number of GCR Iterations
                                               ! analysis step

      INTEGER :: GCR_max_iterations
      INTEGER :: GCR_diagnostics    !RR replace by STASHflag?
      INTEGER :: GCR_its_avg_step(GCR_anal_steps)   
                                             ! GCR Iterations analysis 
                                             ! steps
      INTEGER, DIMENSION(max_numcycles) :: GCR_its_switch
                                             ! GCR Iterations analysis 
                                             ! switch
      INTEGER, DIMENSION(max_numcycles) :: GCR_max_its  
                                             ! GCR Max iterations
      INTEGER, DIMENSION(max_numcycles) :: GCR_min_its  
                                             ! GCR Min iterations
      INTEGER, DIMENSION(max_numcycles) :: GCR_max_time 
                                             ! Timestep number
                                             ! GCR Max iterations
      INTEGER, DIMENSION(max_numcycles) :: GCR_min_time 
                                             ! Timestep number
                                             ! GCR Min iterations
      INTEGER, DIMENSION(max_numcycles) :: GCR_sum_its  
                                             ! Sum iterations
                                             ! over test period
      INTEGER :: GCR_precon_option
      INTEGER :: GCR_n_ADI_pseudo_timesteps

      ! After how many iterations do we restart
      INTEGER :: GCR_Restart_value

      REAL :: GCR_tol_res
      REAL :: GCR_tol_abs
      REAL :: GCR_ADI_pseudo_timestep
      REAL :: G_term_tol
!      G_term_tol             ! Now controlled by namelist input
!
!  Parameters for dynamics-physics cycling
       INTEGER :: NumCycles
       LOGICAL :: L_new_tdisc
       REAL    :: extrp_weight
       REAL    :: GCR_tol_abs2
       REAL    :: GCR_tol_res2
       REAL    :: alpha_1_2, alpha_2_2, alpha_3_2, alpha_4_2
! Parameter for fully-interpolating theta option
       LOGICAL :: L_fint_theta
! Number of physics2 substeps (default is 1)
       INTEGER :: Num_Substeps
! Logical switch to enable phys2 substepping
       Logical :: L_phys2_substep
      COMMON  /RUN_Dyn/                                                 &
     &  L_idealised_data,                                               &
     &  L_Physics, L_phys2_substep, L_Backwards, L_Primitive,           &
     &  L_GCR_cycle_opt, L_Run_With_Physics2, L_perturb_IC_theta,       &
     &  L_mix_ratio, L_new_tdisc, L_fint_theta, IntRand_seed,           &
     &  L_free_slip,                                                    &
     &  L_simple_friction, L_inviscid, L_deep, L_hydrostatic,           &
     &  L_dry, L_adjust_wet, L_LBC_balance, L_lbc_new,                  &
     &  L_geostrophic, L_solver, L_trap_errors, L_uv_zero,              &
     &  L_regular, L_qwaterload, n_rims_to_do, lam_var, phi_var,        &
     &  var_ratio, lam_ratio, phi_ratio, lam_frac, phi_frac,            &
     &  lambda_p_end,  phi_p_end, lambda_u_end,  phi_v_end,             &
     &  dlambda_p_end,  dphi_p_end, dlambda_u_end,  dphi_v_end,         &
     &  alpha_1,alpha_2,alpha_3,alpha_4,                                &
     &  alpha_1_2, alpha_2_2, alpha_3_2, alpha_4_2,                     &
     &  Num_Substeps, NumCycles,                                        &
     &  extrp_weight,                                                   &
     &  GCR_diagnostics, GCR_max_its, GCR_min_its,                      &
     &  GCR_max_time, GCR_min_time,                                     &
     &  GCR_sum_its, GCR_its_switch, GCR_its_avg_step,                  &
     &  GCR_use_tol_abs,GCR_zero_init_guess,GCR_use_residual_Tol,       &
     &  GCR_adi_add_full_soln,L_gcr_fast_x,                             &
     &  L_interp_depart,                                                &
     &  GCR_max_iterations, GCR_precon_option,                          &
     &  GCR_n_ADI_pseudo_timesteps,GCR_Restart_value,GCR_tol_res,       &
     &  GCR_tol_res2,GCR_tol_abs,GCR_tol_abs2,GCR_ADI_pseudo_timestep,  &
     &  G_term_tol
      NAMELIST/RUN_Dyn/                                                 &
     &  L_Physics,L_phys2_substep,L_Backwards,L_Primitive,              &
     &  L_GCR_cycle_opt, L_Run_With_Physics2, L_perturb_IC_theta,       &
     &  L_mix_ratio, L_new_tdisc, L_fint_theta, IntRand_seed,           &
     &  L_free_slip, L_dry, L_adjust_wet,                               &
     &  L_LBC_balance, L_lbc_new,                                       &
     &  L_regular, L_qwaterload, n_rims_to_do, lam_var, phi_var,        &
     &  var_ratio, lam_ratio, phi_ratio, lam_frac, phi_frac,            &
     &  lambda_p_end,  phi_p_end, lambda_u_end,  phi_v_end,             &
     &  dlambda_p_end,  dphi_p_end, dlambda_u_end,  dphi_v_end,         &
     &  alpha_1,alpha_2,alpha_3,alpha_4,                                &
     &  alpha_1_2, alpha_2_2, alpha_3_2, alpha_4_2,                     &
     &  Num_Substeps, NumCycles,                                        &
     &  extrp_weight,                                                   &
     &  GCR_use_tol_abs,GCR_zero_init_guess,GCR_use_residual_Tol,       &
     &  GCR_adi_add_full_soln,L_gcr_fast_x,                             &
     &  L_interp_depart,                                                &
     &  GCR_max_iterations,GCR_diagnostics,GCR_precon_option,           &
     &  GCR_n_ADI_pseudo_timesteps,GCR_Restart_value,GCR_tol_res,       &
     &  GCR_tol_res2,GCR_tol_abs,GCR_tol_abs2,GCR_ADI_pseudo_timestep,  &

     &  G_term_tol, GCR_its_avg_step

! Semi-Lagrangian advection control options:


      INTEGER, PARAMETER :: Number_SL_choices = 4
      INTEGER, PARAMETER :: Theta_SL = 1
      INTEGER, PARAMETER :: moist_SL = 2
      INTEGER, PARAMETER :: Wind_SL = 3
      INTEGER, PARAMETER :: rho_SL = 4

      ! T if conservation to be enforced.
      LOGICAL:: L_conserv(Number_SL_choices)

      ! T if monotonicity to be enforced.
      LOGICAL:: L_mono(Number_SL_choices)

      ! T if high order interpolation scheme to be used.
      LOGICAL:: L_high(Number_SL_choices)

      ! T if high order scheme to be used in Ritchie routine.
      LOGICAL:: L_Ritchie_high

      ! T if monotone scheme to be used in Ritchie routine.
      LOGICAL:: L_Ritchie_mono

      ! T then only perform vector co-ordinate geometry in 2d.
      LOGICAL:: L_2d_sl_geometry


      ! T: sl code bit repoducible with any sensible halo size
      LOGICAL:: L_sl_halo_reprod

      ! improved nonhydrostatic conservation  of q, qcl and qcf
      LOGICAL:: L_moist_nonhydro_conserve

      ! Run tracer advection code with conservation on          
      LOGICAL :: L_conserve_tracers

      ! a code saying which high monotone scheme to use.
      ! 1 = tensor tri-cubic lagrange (j,i,k)
      ! no other options available at present.
      INTEGER :: high_order_scheme(Number_SL_choices)

      ! a code saying which  monotone scheme to use.
      ! 1 = tri-linear order (j,i,k)
      ! no other options available at present.
      INTEGER :: monotone_scheme(Number_SL_choices)
      REAL :: thmono_height       ! top for monotone fully-interp theta
      INTEGER :: thmono_levels ! levels for monotone fully-interp theta
!                 increments used to limit non-interp increments

      ! a code saying which high  order scheme to use in Ritchie routine
      INTEGER :: Ritchie_high_order_scheme

      ! a code saying which monotone scheme to use in Ritchie routine.
      INTEGER :: Ritchie_monotone_scheme

      INTEGER :: Depart_scheme     ! which departure point scheme to use

      ! for the chosen departure point scheme how many iterations/terms
      ! to use
      INTEGER :: Depart_order

      ! used in interpolation code:
      ! number of levels either side of default level to search.
      INTEGER :: interp_vertical_search_tol

      ! Set of look-up tables for searching on variable grids
      INTEGER :: look_lam(max_look)  ! lambda p search
      INTEGER :: look_phi(max_look)  ! phi p search

      ! Set of minumum search grid lengths for variable grids
      REAL ::  recip_dlam   ! smallest delta_lambda p
      REAL ::  recip_dphi   ! smallest delta_phi p

      INTEGER :: halo_lam  ! halo for lamp look-up table
      INTEGER :: halo_phi  ! halo for phip look-up table


      COMMON  /RUN_SL/L_conserv,L_mono,L_high,L_Ritchie_high,           &
     &  L_Ritchie_mono, thmono_height, thmono_levels, L_2d_sl_geometry, &
     &  L_sl_halo_reprod,high_order_scheme,monotone_scheme,             &
     &  Ritchie_high_order_scheme,Ritchie_monotone_scheme,              &
     &  Depart_scheme,Depart_order,interp_vertical_search_tol,          &
     &  recip_dlam, recip_dphi, halo_lam, halo_phi,                     &
     &  look_lam, look_phi,                                             &
     &  L_moist_nonhydro_conserve, L_conserve_tracers,                  &
     &  Instability_diagnostics

      NAMELIST/RUN_SL/L_conserv,L_mono,L_high,L_Ritchie_high,           &
     &  L_Ritchie_mono, thmono_height, L_2d_sl_geometry,                &
     &  L_sl_halo_reprod,high_order_scheme,monotone_scheme,             &
     &  Ritchie_high_order_scheme,Ritchie_monotone_scheme,              &
     &  Depart_scheme,Depart_order,interp_vertical_search_tol,          &
     &  L_moist_nonhydro_conserve, L_conserve_tracers

      ! Diffusion control options:
      LOGICAL :: L_diffusion
      LOGICAL :: L_upper_ramp   ! ramp upper-level diffusion
      LOGICAL :: L_adjust_theta ! activate convective adjustment
      LOGICAL :: L_vdiff_uv    ! activate targeted diffusion uv
      LOGICAL :: L_filter    ! activate polar filter or diffusion
      LOGICAL :: L_filter_incs    ! activate polar filter of incs
      LOGICAL :: L_pftheta   ! activate polar filter of theta
      LOGICAL :: L_pfuv      ! activate polar filter of u,v
      LOGICAL :: L_pfw       ! activate polar filter of w
      LOGICAL :: L_pofil_new  ! activate new polar filter
      LOGICAL :: L_pfexner  ! activate polar filter of Exner pressure
      LOGICAL :: L_pofil_hadgem2 ! run with HadGEM2 polar filter setting   
      LOGICAL :: L_diff_exner  ! activate diffusion of Exner pressure
      LOGICAL :: L_pfcomb   ! combined polar filter/diffusion active
      LOGICAL :: L_sponge   ! activate lateral boundariessponge zones
      LOGICAL :: L_pfincs    ! activate polar filter of incs
      LOGICAL :: L_diff_thermo ! horiz. diffusion of theta
      LOGICAL :: L_diff_wind   ! horiz. diffusion of u, v
      LOGICAL :: L_diff_w      ! horiz. diffusion of w
      LOGICAL :: L_diff_incs   ! horiz. diffusion of increments
      LOGICAL :: L_diff_auto   ! UM calculates diffusion parameters
      LOGICAL :: L_tardiff_q     ! activate targeted diffusion q
      LOGICAL :: L_diff_ctl      ! general diffusion control
! L_diff_ctl == L_diffusion .or. L_cdiffusion .or. L_vertical_diffusion
!               .or. L_divdamp .or.  L_tardiff_q .or. L_diag_print
      LOGICAL :: L_cdiffusion
      LOGICAL :: L_subfilter_horiz  ! subgrid turbulence in horizontal
      LOGICAL :: L_subfilter_vert   ! subgrid turbulence in vertical  
      LOGICAL :: L_subfilter_blend  ! blend diffusion coeffs

      ! value * del^2: order=2
      INTEGER  :: diffusion_order_thermo(max_model_levels)

      ! gives del^4 diffusion
      INTEGER  :: diffusion_order_wind(max_model_levels)
      INTEGER  :: diffusion_order_w(max_model_levels)
      INTEGER  :: diffusion_order_q(max_model_levels)
      INTEGER :: u_begin(0:max_121_rows) ! Sweep control on 121 filter
      INTEGER :: u_end(0:max_121_rows) ! Sweep control on 121 filter
      INTEGER :: v_begin(0:max_121_rows) ! Sweep control on 121 filter
      INTEGER :: v_end(0:max_121_rows) ! Sweep control on 121 filter
      INTEGER :: u_sweeps(max_121_rows) ! Sweep control on 121 filter
      INTEGER :: v_sweeps(max_121_rows) ! Sweep control on 121 filter
      INTEGER :: max_sweeps ! Max sweeps wanted for 121 filter
      INTEGER :: global_u_filter ! Sweep control on 121 filter
      INTEGER :: global_v_filter ! Sweep control on 121 filter
      INTEGER :: diff_order_thermo   ! diffusion order for theta
      INTEGER :: diff_order_wind ! diffusion order for winds
      INTEGER :: diff_timescale_thermo  ! diffusion timescale for theta
      INTEGER :: diff_timescale_wind    ! diffusion timescale for wind
      INTEGER :: vdiffuv_timescale ! diffusion e-folding timesteps
      INTEGER :: vdiffuv_start ! start level  targeted diffusion
      INTEGER :: vdiffuv_end ! end level targeted diffusion
      INTEGER :: adjust_theta_start ! start level convective adjustment
      INTEGER :: adjust_theta_end   ! end levelconvective adjustment
      INTEGER :: top_filt_start ! start level upper-level diffusion
      INTEGER :: top_filt_end   ! end level upper-level diffusion
      INTEGER :: sponge_ew   ! left/right boundaries sponge zone width
      INTEGER :: sponge_ns   ! north/south boundaries sponge zone width
      INTEGER :: sponge_power ! sponge zone weighting order
      INTEGER :: tardiffq_test ! test level test w targetted diffusion
      INTEGER :: tardiffq_start ! start level test w targetted diffusion
      INTEGER :: tardiffq_end ! end level test w targetted diffusion
      INTEGER :: turb_startlev_horiz   ! 1st lev for horiz subgrid turb
      INTEGER :: turb_endlev_horiz    ! last lev for horiz subgrid turb
      INTEGER :: turb_startlev_vert   ! 1st lev for vert subgrid turb
      INTEGER :: turb_endlev_vert     ! last lev for vert subgrid turb

      !  level - assume surfaces are horizontal
      INTEGER  :: horizontal_level
      INTEGER  ::  tar_horizontal ! steep slope test targeted diffusion

      REAL :: diffusion_coefficient_thermo(max_model_levels)
      REAL :: diffusion_coefficient_wind(max_model_levels)
      REAL :: diffusion_coefficient_w(max_model_levels)
      REAL :: diffusion_coefficient_q(max_model_levels)
      REAL :: vdiffuv_test ! test to activate shear diffusion of u,v
      REAL :: vdiffuv_factor ! vertical diffusion coeff for shear diff
      REAL :: scale_ratio ! Pass control on 121 filter
      REAL :: ref_lat_deg ! Reference lat for auto diffusion
      REAL :: top_diff    !  upper-level diffusion coefficient
      REAL :: up_diff_scale    !  upper-level diffusion ramping factor

      REAL :: up_diff(max_updiff_levels) ! upper-level diffusion coeff
      REAL :: sponge_wts_ew(max_sponge_width) ! sponge weights
      REAL :: sponge_wts_ns(max_sponge_width) ! sponge weights

      REAL :: diff_coeff_ref ! EW diffusion coefficient at polar cap
      REAL :: diff_coeff_thermo  ! NS theta diffusion coeff
      REAL :: diff_coeff_wind    ! NS u,v diffusion coeff
      REAL :: diff_coeff_phi    ! North-South diffusion coefficient
!   reference latitudes for filtering and diffusion
      REAL :: polar_cap  ! Apply 1-2-1 filter polewards
      REAL :: tardiffq_factor ! targeted diffusion coefficient
      REAL :: diff_factor
      REAL :: mix_factor

      ! Divergence damping control options:
      LOGICAL :: L_divdamp

      REAL :: div_damp_coefficient(max_model_levels)

      ! Polar filter control options:
      LOGICAL :: L_polar_filter
      ! T: use polar filter to filter increment
      LOGICAL :: L_polar_filter_incs

      REAL :: polar_filter_north_lat_limit
      REAL :: polar_filter_south_lat_limit
      REAL :: polar_filter_coefficient

      ! amount in radians to increment start latitude by per sweep
      REAL :: polar_filter_step_per_sweep

      ! max latitude at which filter can star
      REAL :: polar_filter_lat_limit

      ! number of sweeps of filter to do
      INTEGER :: polar_filter_n_sweeps

      ! Vertical Diffusion control options:
      LOGICAL :: L_vertical_diffusion
      LOGICAL :: L_ramp

      INTEGER :: level_start_wind
      INTEGER :: level_stop_wind
      INTEGER :: level_start_theta
      INTEGER :: level_stop_theta
      INTEGER :: level_start_q
      INTEGER :: level_stop_q

      REAL :: vert_diffusion_coeff_wind
      REAL :: vert_diffusion_coeff_theta
      REAL :: vert_diffusion_coeff_q
      REAL :: ramp_lat_radians

      ! Moisture resetting control options:
      LOGICAL :: l_qpos         ! logical to run qpos code

      ! true to do Method 2 false to do Method 1
      LOGICAL :: l_q_pos_local


      REAL :: qlimit         !    lowest allowed value of q

      COMMON  /RUN_Diffusion/L_diffusion,diffusion_order_thermo,        &
     &  diffusion_order_wind,diffusion_order_w,diffusion_order_q,       &
     &  diffusion_coefficient_thermo,diffusion_coefficient_wind,        &
     &  diffusion_coefficient_w,diffusion_coefficient_q,                &
     &  diff_order_thermo, diff_timescale_thermo, diff_coeff_thermo,    &
     &  diff_order_wind, diff_timescale_wind, diff_coeff_wind,          &
     &  horizontal_level, tar_horizontal, tropics_deg,                  &
     &  L_cdiffusion,L_ramp, ramp_lat_radians,                          &
     &  L_divdamp,div_damp_coefficient,                                 &
     &  L_polar_filter,L_polar_filter_incs,                             &
     &  polar_filter_north_lat_limit,polar_filter_south_lat_limit,      &
     &  polar_filter_coefficient,polar_filter_step_per_sweep,           &
     &  polar_filter_lat_limit,polar_filter_n_sweeps,                   &
     &  l_qpos,  l_q_pos_local, qlimit,                                 &
     &  L_diff_ctl, L_tardiff_q, w_conv_limit,                          &
     &  tardiffq_factor, tardiffq_test, tardiffq_start, tardiffq_end,   &
     &  L_subfilter_horiz, L_subfilter_vert, L_subfilter_blend,         &
     &  diff_factor, mix_factor, turb_startlev_horiz, turb_endlev_horiz,&
     &  turb_startlev_vert, turb_endlev_vert,                           &
     &  L_vertical_diffusion, L_pofil_hadgem2,                          &
     &  L_pofil_new, L_pfcomb, L_pftheta, L_pfuv, L_pfw, L_pfincs,      &
     &  L_pfexner, L_filter, L_filter_incs, L_diff_exner, L_diff_auto,  &
     &  L_diff_incs, L_diff_thermo, L_diff_wind, L_diff_w,              &
     &  L_vdiff_uv, L_adjust_theta,                                     &
     &  level_start_wind, level_stop_wind,                              &
     &  level_start_theta, level_stop_theta,                            &
     &  level_start_q, level_stop_q,                                    &
     &  polar_cap, scale_ratio, ref_lat_deg, max_sweeps,                &
     &  L_upper_ramp, top_filt_start, top_filt_end, top_diff,           &
     &  up_diff, up_diff_scale, sponge_wts_ew, sponge_wts_ns,           &
     &  u_begin, u_end, v_begin, v_end, u_sweeps, v_sweeps,             &
     &  global_u_filter, global_v_filter,                               &
     &  diff_coeff_ref, diff_coeff_phi, vdiffuv_timescale,              &
     &  vdiffuv_test, vdiffuv_factor, vdiffuv_start, vdiffuv_end,       &
     &  adjust_theta_start, adjust_theta_end,                           &
     &  L_sponge, sponge_power, sponge_ew, sponge_ns,                   &
     &  vert_diffusion_coeff_wind, vert_diffusion_coeff_theta,          &
     &  vert_diffusion_coeff_q

      NAMELIST/RUN_Diffusion/L_diffusion,diffusion_order_thermo,        &
     &  diffusion_order_wind,diffusion_order_w,diffusion_order_q,       &
     &  diffusion_coefficient_thermo,diffusion_coefficient_wind,        &
     &  diffusion_coefficient_w,diffusion_coefficient_q,                &
     &  diff_order_thermo, diff_timescale_thermo,                       &
     &  diff_order_wind, diff_timescale_wind,                           &
     &  horizontal_level, tar_horizontal,                               &
     &  L_cdiffusion,L_ramp, ramp_lat_radians,                          &
     &  L_divdamp,div_damp_coefficient,                                 &
     &  L_polar_filter,L_polar_filter_incs,                             &
     &  polar_filter_north_lat_limit,polar_filter_south_lat_limit,      &
     &  polar_filter_coefficient,polar_filter_step_per_sweep,           &
     &  polar_filter_lat_limit,polar_filter_n_sweeps,                   &
     &  l_qpos,  l_q_pos_local, qlimit,                                 &
     &  L_diff_ctl, L_tardiff_q, w_conv_limit,                          &
     &  tardiffq_factor, tardiffq_test, tardiffq_start, tardiffq_end,   &
     &  L_subfilter_horiz, L_subfilter_vert, L_subfilter_blend,         &
     &  diff_factor, mix_factor, turb_startlev_horiz, turb_endlev_horiz,&
     &  turb_startlev_vert, turb_endlev_vert,                           &
     &  L_vertical_diffusion,                                           &
     &  L_pftheta, L_pfuv, L_pfw, L_pfincs, L_pfexner, L_pofil_hadgem2, &
     &  L_diff_incs, L_diff_thermo, L_diff_wind, L_diff_w,              &
     &  level_start_wind, level_stop_wind,                              &
     &  level_start_theta, level_stop_theta,                            &
     &  level_start_q, level_stop_q,                                    &
     &  L_pofil_new, L_diff_auto,                                       &
     &  diff_coeff_ref, polar_cap, scale_ratio, ref_lat_deg,            &
     &  max_sweeps, L_upper_ramp, up_diff_scale, top_diff,              &
     &  top_filt_start, top_filt_end, L_vdiff_uv, vdiffuv_timescale,    &
     &  vdiffuv_test, vdiffuv_factor, vdiffuv_start, vdiffuv_end,       &
     &  L_adjust_theta, adjust_theta_start, adjust_theta_end,           &
     &  L_sponge, sponge_power, sponge_ew, sponge_ns,                   &
     &  vert_diffusion_coeff_wind, vert_diffusion_coeff_theta,          &
     &  vert_diffusion_coeff_q,                                         &
     &  L_diag_print, L_diag_print_ops, L_print_pe,                     &
     &  L_print_w, L_print_wmax, L_print_lapse, L_print_theta1,         &
     &  L_print_div, L_diag_wind, L_print_shear, L_print_max_wind,      &
     &  L_diag_noise, L_diag_L2norms, L_diag_L2helm,                    &
     &  norm_lev_start, norm_lev_end, first_norm_print,                 &
     &  print_step, diag_interval, w_print_limit, L_Flush6

!------------------   Dynamical core   -------------------------------
! Suarez-Held variables
      REAL :: SuHe_newtonian_timescale_ka
      REAL :: SuHe_newtonian_timescale_ks
      REAL :: SuHe_pole_equ_deltaT
      REAL :: SuHe_static_stab
      REAL :: base_frictional_timescale
      REAL :: SuHe_sigma_cutoff
      REAL :: SuHe_level_weight(max_model_levels)
      REAL :: friction_level(max_model_levels)

      INTEGER :: SuHe_relax
      INTEGER :: SuHe_fric

      LOGICAL :: L_SH_Williamson

      COMMON/Run_Dyncore/                                               &
     &  SuHe_newtonian_timescale_ka, SuHe_newtonian_timescale_ks,       &
     &  SuHe_pole_equ_deltaT, SuHe_static_stab,                         &
     &  base_frictional_timescale, SuHe_sigma_cutoff,                   &
     &  L_SH_Williamson, SuHe_relax, SuHe_fric,                         &
     &  SuHe_level_weight, frictional_timescale, friction_level

!------------------  Idealised model   ----------------------------

      INTEGER,PARAMETER:: max_num_profile_data = 100
      INTEGER,PARAMETER:: max_num_force_times = 100
      INTEGER,PARAMETER:: idl_max_num_bubbles = 3

! Idealised  variables
      REAL :: h_o
      REAL :: h_o_actual  ! height of growing mountain
      REAL :: h_o_per_step ! height change per step of growing mountain
      REAL :: lambda_fraction
      REAL :: phi_fraction
      REAL :: half_width_x
      REAL :: half_width_y
      REAL :: Witch_power
      REAL :: plat_size_x
      REAL :: plat_size_y
      REAL :: height_domain
      REAL :: delta_x
      REAL :: delta_y
      REAL :: big_factor
      REAL :: mag
      REAL :: first_theta_height
      REAL :: thin_theta_height
      REAL :: p_surface
      REAL :: theta_surface
      REAL :: dtheta_dz1(3)
      REAL :: height_dz1(3)
      REAL :: Brunt_Vaisala
      REAL :: u_in(4)
      REAL :: v_in(4)
      REAL :: height_u_in(3)
      REAL :: u_ramp_start
      REAL :: u_ramp_end
      REAL :: ujet_lat
      REAL :: ujet_width
      REAL :: t_horizfn_data(10)
      REAL :: q1
      REAL :: theta_ref(max_model_levels)
      REAL :: rho_ref(max_model_levels)
      REAL :: exner_ref(max_model_levels + 1)
      REAL :: q_ref(max_model_levels)
      REAL :: u_ref(max_model_levels)
      REAL :: v_ref(max_model_levels)
      REAL :: z_orog_print(0:max_model_levels)
      REAL :: f_plane
      REAL :: ff_plane
      REAL :: r_plane
      REAL :: zprofile_data(max_num_profile_data)
      REAL :: tprofile_data(max_num_profile_data)
      REAL :: qprofile_data(max_num_profile_data)
      REAL :: z_uvprofile_data(max_num_profile_data)
      REAL :: uprofile_data(max_num_profile_data)
      REAL :: vprofile_data(max_num_profile_data)
      REAL :: tforce_time_interval
      REAL :: qforce_time_interval
      REAL :: uvforce_time_interval
      REAL :: newtonian_timescale
      REAL :: z_tforce_data(max_num_profile_data)
      REAL :: z_qforce_data(max_num_profile_data)
      REAL :: z_uvforce_data(max_num_profile_data)
      REAL :: tforce_data(max_num_profile_data, max_num_force_times)
      REAL :: qforce_data(max_num_profile_data, max_num_force_times)
      REAL :: uforce_data(max_num_profile_data, max_num_force_times)
      REAL :: vforce_data(max_num_profile_data, max_num_force_times)
      REAL :: tforce_data_modlev(max_model_levels, max_num_force_times)
      REAL :: qforce_data_modlev(max_model_levels, max_num_force_times)
      REAL :: uforce_data_modlev(max_model_levels, max_num_force_times)
      REAL :: vforce_data_modlev(max_model_levels, max_num_force_times)
      REAL :: pforce_time_interval
      REAL :: p_surface_data(max_num_force_times)
      REAL :: perturb_factor
      REAL :: perturb_magnitude_t
      REAL :: perturb_magnitude_q
      REAL :: perturb_height(2)
      REAL :: orog_hgt_lbc
      REAL :: zprofile_orog
      REAL :: hf
      REAL :: cool_rate
      REAL :: IdlSurfFluxSeaParams(10) ! Idealised surface flux params
      REAL :: roughlen_z0m   
      REAL :: roughlen_z0h
      ! Idealised bubbles
      REAL :: idl_bubble_max(idl_max_num_bubbles) ! Bubble max amplitude
      REAL :: idl_bubble_height(idl_max_num_bubbles)  ! Bubble height
      REAL :: idl_bubble_width(idl_max_num_bubbles)   ! Bubble width
      REAL :: idl_bubble_depth(idl_max_num_bubbles)   ! Bubble depth
      ! Bubble x-offset, y-offset in normalised units (0:1)
      ! (0.5=domain centre)
      REAL :: idl_bubble_xoffset(idl_max_num_bubbles)
      REAL :: idl_bubble_yoffset(idl_max_num_bubbles)
      REAL :: DMPTIM, HDMP, ZDMP   ! Damping layer values
      REAL :: u_geo, v_geo         ! Geostrophic wind


      INTEGER :: n_rims_to_do   ! rim size for LAM
      INTEGER :: surface_type
      INTEGER :: grow_steps
      INTEGER :: grid_number
      INTEGER :: grid_flat
      INTEGER :: tprofile_number
      INTEGER :: qprofile_number
      INTEGER :: uvprofile_number
      INTEGER :: num_profile_data
      INTEGER :: num_uvprofile_data
      INTEGER :: t_horizfn_number
      INTEGER :: uv_horizfn_number
      INTEGER :: pforce_option
      INTEGER :: num_pforce_times
      INTEGER :: tforce_option
      INTEGER :: qforce_option
      INTEGER :: uvforce_option
      INTEGER :: num_tforce_levels
      INTEGER :: num_tforce_times
      INTEGER :: num_qforce_levels
      INTEGER :: num_qforce_times
      INTEGER :: num_uvforce_levels
      INTEGER :: num_uvforce_times
      INTEGER :: IdlSurfFluxSeaOption  ! Idealised surface flux option
      INTEGER :: first_constant_r_rho_level_new
      INTEGER :: big_layers
      INTEGER :: transit_layers
      INTEGER :: mod_layers
      INTEGER :: idl_bubble_option(idl_max_num_bubbles) ! Bubble option
      INTEGER :: idl_interp_option  ! Profile interpolation option
      INTEGER :: perturb_type

      LOGICAL :: L_initialise_data
      LOGICAL :: L_constant_dz
      LOGICAL :: L_trivial_trigs !.false. for Cartesian coords (lat=0.0)
      LOGICAL :: L_idl_bubble_saturate(idl_max_num_bubbles)
      LOGICAL :: L_fixed_lbcs
      LOGICAL :: L_fix_orog_hgt_lbc
      LOGICAL :: L_pressure_balance
      LOGICAL :: L_wind_balance
      LOGICAL :: L_rotate_winds
      LOGICAL :: L_polar_wind_zero
      LOGICAL :: L_vert_Coriolis
      LOGICAL :: L_rotating     ! .true. for Earth's rotation
      LOGICAL :: L_perturb      ! add random perturb. to surface theta
      LOGICAL :: L_code_test    ! User switch for testing code
      LOGICAL :: L_pforce
      LOGICAL :: L_baroclinic
      LOGICAL :: L_cyclone
      LOGICAL :: L_force
      LOGICAL :: L_force_lbc
      LOGICAL :: L_perturb_t
      LOGICAL :: L_perturb_q
      LOGICAL :: L_perturb_correlate_tq
      LOGICAL :: L_perturb_correlate_vert
      LOGICAL :: L_perturb_correlate_time
      LOGICAL :: L_damp      ! Logical for damping layer
      LOGICAL :: L_geo_for ! Logical for geostrophic wind forcing
      LOGICAL :: L_bomex     ! Logical for BOMEX set up
      LOGICAL :: L_spec_z0   ! specification of roughness length    

      COMMON  /RUN_Ideal/                                               &
     &  h_o, h_o_actual, h_o_per_step,                                  &
     &  lambda_fraction, phi_fraction, half_width_x, half_width_y,      &
     &  Witch_power, plat_size_x, plat_size_y,                          &
     &  height_domain, delta_x, delta_y, big_factor, mag,               &
     &  first_theta_height, thin_theta_height, p_surface,               &
     &  theta_surface, dtheta_dz1, height_dz1, Brunt_Vaisala,           &
     &  u_in, v_in, height_u_in, u_ramp_start, u_ramp_end, q1,          &
     &  ujet_lat, ujet_width,                                           &
     &  t_horizfn_number, t_horizfn_data, uv_horizfn_number,            &
     &  u_ref, v_ref, theta_ref, exner_ref, rho_ref, q_ref,             &
     &  z_orog_print, grow_steps,                                       &
     &  surface_type, grid_number, grid_flat,                           &
     &  tprofile_number, qprofile_number, uvprofile_number,             &
     &  num_profile_data, num_uvprofile_data,                           &
     &  tforce_option, qforce_option, uvforce_option,                   &
     &  num_tforce_levels, num_tforce_times,                            &
     &  num_qforce_levels, num_qforce_times,                            &
     &  num_uvforce_levels, num_uvforce_times,                          &
     &  L_pforce, pforce_option, num_pforce_times,                      &
     &  first_constant_r_rho_level_new,                                 &
     &  big_layers, transit_layers, mod_layers,                         &
     &  zprofile_data, tprofile_data, qprofile_data,                    &
     &  z_uvprofile_data, uprofile_data, vprofile_data,                 &
     &  tforce_time_interval, qforce_time_interval,                     &
     &  uvforce_time_interval, newtonian_timescale,                     &
     &  z_tforce_data, z_qforce_data, z_uvforce_data,                   &
     &  tforce_data, qforce_data, uforce_data, vforce_data,             &
     &  tforce_data_modlev, qforce_data_modlev,                         &
     &  uforce_data_modlev, vforce_data_modlev,                         &
     &  pforce_time_interval, p_surface_data,                           &
     &  L_initialise_data,                                              &
     &  L_perturb_t, perturb_magnitude_t,                               &
     &  L_perturb_q, perturb_magnitude_q,                               &
     &  L_perturb_correlate_tq,                                         &
     &  L_perturb_correlate_vert,                                       &
     &  L_perturb_correlate_time,                                       &
     &  perturb_type, perturb_height,                                   &
     &  L_constant_dz, L_polar_wind_zero,                               &
     &  L_wind_balance, L_rotate_winds,                                 &
     &  IdlSurfFluxSeaOption, IdlSurfFluxSeaParams,                     &
     &  L_spec_z0, roughlen_z0m, roughlen_z0h,                          &
     &  L_pressure_balance, L_vert_Coriolis,                            &
     &  cool_rate, L_force, L_force_lbc,                                &
     &  zprofile_orog, idl_interp_option, hf,                           &
     &  L_fix_orog_hgt_lbc, orog_hgt_lbc,                               &
     &  L_trivial_trigs, f_plane, ff_plane, r_plane,                    &
     &  idl_bubble_option, idl_bubble_max                               &
     &, idl_bubble_height, idl_bubble_width, idl_bubble_depth           &
     &, idl_bubble_xoffset,idl_bubble_yoffset                           &
     &, L_idl_bubble_saturate,                                          &
     &  L_rotating, L_fixed_lbcs, L_code_test,                          &
     &  L_baroclinic, L_cyclone,                                        &
     &  L_damp, L_geo_for, L_bomex,                                     &
     &  DMPTIM, HDMP, ZDMP,                                             &
     &  u_geo, v_geo
! CRUNTIMC end
! CCONSTS start
! Description:
!    This file contains declarations for derived constants within
!   the atmospheric model. Where necessary PARAMETERS are defined to
!   dimension these constants. All constants are placed in the common
!   block CDERIVED, except hardwired constants, e.g. ETA_SPLIT and LENs.
!   file CMAXSIZE must be included first
!
!   The derived constants are calculated in the routine SETCONA1.
!
! PA, WI      <- programmer of some or all of previous code or changes
!
!  Model            Modification history from model version 3.0:
! version  Date
! 3.2   26/03/93  Remove resolution dependent variables for dynamic
!                 allocation. R.Rawlins
!   4.0   20/07/95  Sizes of tables expanded for version 3A
!                   of the radiation in sections 1 and 2. J.M. Edwards
!   5.0   21/06/99  Remove obsolete constants, for C-P C-grid dynamics.
!                   M L Gallani
!   5.1   25/02/00  Replace Data ETA_SPLIT by hard-wire h_split in
!                   Readlsta. Also remove obsolete references to
!                   radiation tables. R Rawlins
! Logical component: F011

      ! No of cloud types ie low/med/high
      INTEGER, PARAMETER :: NUM_CLOUD_TYPES = 3

      ! derived constants:
      INTEGER :: LOW_BOT_LEVEL      ! Bottom level of lowest cloud type
      INTEGER :: LOW_TOP_LEVEL      ! Top      "    "   "       "    "
      INTEGER :: MED_BOT_LEVEL      ! Bottom   "    "  med      "    "
      INTEGER :: MED_TOP_LEVEL      ! Top      "    "   "       "    "
      INTEGER :: HIGH_BOT_LEVEL     ! Bottom   "    "  top      "    "
      INTEGER :: HIGH_TOP_LEVEL     ! Top      "    "   "       "    "

      ! height values to split model levels into l/m/h cloud
      REAL ::    h_split(NUM_CLOUD_TYPES+1)

      LOGICAL :: ELF                ! T if atmosphere model on LAM grid

      ! Constants for dynamics output independent of resolution but
      ! dependent on choice of levels for output.
      REAL :: REQ_THETA_PV_LEVS(MAX_REQ_THPV_LEVS)

      COMMON /CDERIVED/                                                 &
     &  h_split,LOW_BOT_LEVEL,LOW_TOP_LEVEL,MED_BOT_LEVEL,MED_TOP_LEVEL,&
     &  HIGH_BOT_LEVEL, HIGH_TOP_LEVEL,ELF,REQ_THETA_PV_LEVS
! CCONSTS end
! ----------------------- header file: CNTLATM  -----------------------
! Description: Control variables for the Atmosphere model (read only).
!              Contains logical switches for science options needed
!              for addressing calculations and intermediate control.
!              Predominantly used for holding logical flags set up by
!              the UMUI - read in by a namelist, but can also hold
!              derived control flags set by the model.
!              [Note that CRUNTIMC holds accompanying run-time
!              constants.]
!
! Author : R. Rawlins
!
! History:
! Version  Date      Comment.
!  5.0 20/04/99  Extensive revision to earlier comdeck for conversion
!                to C-P 'C' dynamics grid. R. Rawlins
!  5.1  4/02/00  Restore energy correction switches removed at UM5.0
!                Also added additional switches for mass and moisture
!                R A Stratton.
!  5.1 13/04/00  IAU control moved to CTFilt. Adam Clayton
!  5.2 22/08/00  Reinstate Murk and Tracer switches. P.Selwood.
!  5.2 15/11/00  Reintroduce logicals for MOSES 2 and triffid. M.Best.
!  5.3   12/10/01   Remove hard-wired L_trivial_trigs.
!                   put problem_number in CNTLATM.    Terry Davies
!  5.3 27/07/01  Add logical switch L_MURK_RAD   S. Cusack
!  5.3    06/01  Introduce and hardwire logical for setting of
!                leads temperature.                         Nic Gedney
!  5.3 15/10/01  Added L_USE_METHOX switch. David Jackson
!  5.3 29/08/01  Sulphate and sea-salt logicals split into separate
!                versions for different processes.  A. Jones
!  5.3 09/04/01  Add logical switch for spatial degradation of
!                E-S radiation calculations.             S. Cusack
!
!  5.3 19/06/01   Stuff to handle tropopause-based ozone added
!                 -- see Radiation block              Dave Tan
!  5.4 15/08/02   Reconfiguration use changes. P.Selwood.
!  5.4  2/08/02  Add logical switch for PC2 cloud and condensation
!                scheme.                              Damian Wilson
!  5.4 24/10/02  Moved L_GWD and L_USE_USSP from CRUNTIMC. S. Webster
!  5.4 11/03/02  Remove comment lines from after #include
!                                                 S. Carroll
!  5.5 06/02/03  River routing support. P.Selwood.
!  5.5 05/02/03  Add logicals for biomass aerosol scheme     P Davison

!  5.5 17/02/03  Add two large-scale hydrology logicals.
!                L_TOP and L_PDM.                  Nic Gedney
!  5.5 08/01/03  Remove L_3DVAR, L_4DVAR, L_3DVAR_BG, LSINGLE_HYDROL,
!                L_H2_SULPH and L_LSPICE_BDY. D Robinson
!  5.5 21/01/03  Move L_USE_TPPS_OZONE to NLSTCATM namelist. D Robinson
!  5.5 13/03/03  Move I_OZONE_INT here from CRUNTIMC.
!                                                  Jean-Claude Thelen
!  5.5 03/02/03  Include L_mcr logicals.    R.M.Forbes
!  5.5 19/02/03  Remove redundant L_BL_LSPICE and L_RMBL   A.Lock
!  6.0 30/07/03  Include l_pc2_lbc for cloud frac. lbcs. Damian Wilson
!  6.0 19/08/03  Add runtime controls for 4 physics sections
!  6.1  02/08/04  Add logicals for stochem coupling. R Barnes
!  6.1 13/05/04  Add super_array_size variable                 A.Malcolm
!  6.1 07/04/04  Add logical for autoconversion de-biasing.   A. Jones
!  6.1 07/04/04  Add logicals for interactive DMS emissions.  A. Jones
!  6.2 25/01/06  Add iteration count logical test variable    T.Edwards
!  6.2 15/11/05  Add logical for aerosol optical depth     N. Bellouin
!  6.2 23/11/05  Add logical for RH and hygroscopic aerosols N.Bellouin
!  6.2 01/03/06  Add L_UPPER_STRATQ switch - David Jackson
!  6.2 06/01/06  Add logicals for seaice albedo options. J.Ridley
!  6.2 23/02/06  Add logicals for UKCA sub-model.  F.O'Connor
!  6.1 07/04/04  Add logical for RothC temperature function.  C.D. Jones
!  6.2 21/07/05  Add moisture_array_size variable              A.Malcolm
!  6.2 01/10/05  Include L_murk_lbc logical.    R.M.Forbes
!  6.2 07/11/05  Add L_MOD_BARKER_ALBEDO and L_USE_SPEC_SEA
!                to NLSTCATM.        James Manners
!  6.2 25/01/06  Add L_INHOM_CLOUD to NLSTCATM.    James Manners
!  6.2 09/03/06  Add logicals for inland basin rerouting. P.Falloon
!  6.2 24/02/06  Add logical for 10m windspeed pass A2O  J.Gunson

!  6.2 24/10/05 Functionality for radiative forcing, timestepping
!               and radiances under versions 3C and 3Z of radiation
!               code added                 (J.-C. Thelen)
!  6.2 07/11/05 Add L_use_orog_corr to NLSTCATM.     James Manners
!  6.2 24/02/06  Add logical to allow DMS flux from ocean model J.Gunson
!  6.4 19/01/07 Removed a comment relating to A05_3c scheme. R A Stratton
!  6.4 08/01/06 Include Brooks cloud fraction logicals. Damian Wilson
!  6.4 10/01/07 Include mixing ratio control logical. Damian Wilson
!  6.6.2 10/06/09 Logicals for reading solar/volcanic forcing. Chris Jones
!------------------   General:  -------------------------------------
      INTEGER Model_domain        ! Domain of atmosphere model:
!                                   global,LAM,cyclic LAM,single column
! Model_domain meaningful names
! DOMTYP contains different model domain types
!
! Author : P.Burton
! History:
! Version  Date      Comment.
! 5.0      15/04/99  New comdeck
! 5.2      15/11/00  add bi_cyclic_lam domain   A. Malcolm

      INTEGER,PARAMETER:: mt_global        = 1
      INTEGER,PARAMETER:: mt_lam           = 2
      INTEGER,PARAMETER:: mt_cyclic_lam    = 3
      INTEGER,PARAMETER:: mt_bi_cyclic_lam = 4
      INTEGER,PARAMETER:: mt_single_column = 5
! DOMTYP end

      INTEGER :: problem_number      ! type of problem to be solved

      INTEGER MAXSECTS            ! Max. no. of code sections
      PARAMETER (MAXSECTS=99)
      CHARACTER*3 H_SECT(0:MAXSECTS) ! Array of code section-versions

      ! Physics:   ------------------------------------
      LOGICAL :: l_ssice_albedo     ! Sea-ice albedo affected by snow

      LOGICAL :: l_sice_heatflux    ! Semi-impl update of sea-ice temp

      LOGICAL :: l_sice_meltponds ! Sea-ice albedo affected by meltponds

      LOGICAL :: l_sice_scattering  ! Sea-ice albedo affected scattering

      LOGICAL :: l_sice_hadgem1a ! HadGEM1 sea-ice albedo bug corrected

      LOGICAL :: L_NEG_TSTAR        ! Test for -ve surface temperature.
!
      ! Use sulphate aerosol in microphysics
      LOGICAL :: L_USE_SULPHATE_AUTOCONV

      ! Use sea-salt aerosol in microphysics
      LOGICAL :: L_USE_SEASALT_AUTOCONV

      ! Use soot aerosol in microphysics
      LOGICAL :: L_USE_SOOT_AUTOCONV

      ! Use biomass aerosol in microphysics
      LOGICAL :: L_USE_BMASS_AUTOCONV
      
      ! Use fossil-fuel organic carbon in microphysics
      LOGICAL :: L_USE_OCFF_AUTOCONV
      
      ! Use autoconversion de-biasing scheme in microphysics
      LOGICAL :: L_AUTO_DEBIAS
      ! Use sulphate aerosol no. in S-cycle
      LOGICAL :: L_USE_SULPHATE_SULPC

      ! Use sea-salt aerosol no. in S-cycle
      LOGICAL :: L_USE_SEASALT_SULPC

      ! Use soot aerosol no. in S-cycle
      LOGICAL :: L_USE_SOOT_SULPC

      ! Use biomass aerosol no. in S-cycle
      LOGICAL :: L_USE_BMASS_SULPC
      
      ! Use fossil-organic carbon aerosol no. in S-cycle
      LOGICAL :: L_USE_OCFF_SULPC
      
      ! Energy correction:
      LOGICAL :: L_emcorr    ! T: turns on energy correction code
      LOGICAL :: LMASS_corr  ! T: Apply mass correction
      LOGICAL :: LQT_corr    ! T: Apply total moisture correction
      LOGICAL :: LEMQ_print  ! T: Print additional info from em code
      LOGICAL :: LENERGY     ! T: if timestep to cal energy correction
      LOGICAL :: LFLUX_RESET ! T: if timestep to reset flux array in D1

      ! number of model timesteps per energy correction period.
      INTEGER :: A_ENERGYSTEPS

      ! Radiation:

      LOGICAL :: L_radiation     !  F: Turns off radiation code
      LOGICAL :: L_MICROPHY           !  Microphysics in sw rad scheme

!     Use mineral dust in radiation calculations
      LOGICAL :: L_USE_DUST

!     Use biogenic aerosol in radiation code
      LOGICAL :: L_USE_BIOGENIC

      ! Use SO4 aerosol from sulphur cycle for direct/indirect effect
      ! in radiation, the latter for both SW and LW.
      LOGICAL :: L_USE_SULPC_DIRECT
      LOGICAL :: L_USE_SULPC_INDIRECT_SW
      LOGICAL :: L_USE_SULPC_INDIRECT_LW
      ! Indirect radiative effect of sea-salt
      LOGICAL :: L_USE_SEASALT_INDIRECT

      ! Direct radiative effect of sea-salt
      LOGICAL :: L_USE_SEASALT_DIRECT
      LOGICAL :: L_USE_SOOT_DIRECT  ! direct radiative effects of soot
      LOGICAL :: L_USE_SOOT_INDIRECT  ! indirect effects of soot
      ! Use biomass aerosol for direct/indirect effect in radiation.
      LOGICAL :: L_USE_BMASS_DIRECT
      LOGICAL :: L_USE_BMASS_INDIRECT
      
      ! Use fossil-fuel organic carbon aerosol for direct/indirect
      ! effect in radiation
      LOGICAL :: L_USE_OCFF_DIRECT
      LOGICAL :: L_USE_OCFF_INDIRECT
      
!     Use aerosol climatologies in radiation instead of prognostic variables
!     Set on a species by species basis
      LOGICAL :: L_USE_ARCLBIOM   ! biomass burning aerosol
      LOGICAL :: L_USE_ARCLBLCK   ! black carbon
      LOGICAL :: L_USE_ARCLSSLT   ! sea salt
      LOGICAL :: L_USE_ARCLSULP   ! sulpahtes
      LOGICAL :: L_USE_ARCLDUST   ! mineral dust
      LOGICAL :: L_USE_ARCLOCFF   ! organic carbon (fossil fuel)
      LOGICAL :: L_USE_ARCLDLTA   ! delta aerosol

      ! Aerosol optical depth diagnostic was requested
      LOGICAL :: L_USE_AOD

      ! Clear-sky relative humidity is to be used instead of
      ! grid-box mean RH, for hygroscopic aerosols
      LOGICAL :: L_USE_CLEARRH

      ! controls the use of spatial degradation of radiation calc.
      LOGICAL :: L_rad_deg
      LOGICAL :: L_CTILE       ! Switch for coastal tiling.
!                              ! If TRUE then land and sea can
!                              ! coexist in the same gridbox.
!                              ! If FALSE the land fraction
!                              ! must be equal to 0 to 1
      LOGICAL :: L_TOP         ! If TRUE then TOPMODEL-based
!                              ! hydrology scheme.
      LOGICAL :: L_PDM         ! If TRUE then PDM-based
!                              ! hydrology scheme.
      LOGICAL :: L_SOIL_SAT_DOWN ! If TRUE then super-saturated soil
!                                ! moisture moves downward, else upward

      INTEGER :: H_SWBANDS    ! Number of shortwave radiation bands
      INTEGER :: H_LWBANDS    ! Number of longwave radiation bands
      INTEGER :: A_SW_RADSTEP ! Number of advection steps per SW step
      INTEGER :: A_LW_RADSTEP ! Number of advection steps per LW step
! Number of advection steps per prognostic/diagnostic SW and LW step.
!'Prognostic' and 'Diagnostic' refer to the frequency of the calls
! to radiation code in the Unified Model.
! In the case of time stepping prognostic and diagnostic refer to the
! slow and fast radiative timestep respectively. In the case of radiative
! forcing they refer to the prognostic and diagnostic calls to radiation.
      INTEGER :: A_SW_RADSTEP_DIAG
! Number of advection steps per 'fast' SW step (3C)
      INTEGER :: A_LW_RADSTEP_DIAG
! Number of advection steps per 'fast' LW step (3C)
      INTEGER :: A_SW_RADSTEP_PROG
! Number of advection steps per 'slow' LW step (3C)
      INTEGER :: A_LW_RADSTEP_PROG
! Number of advection steps per 'slow' LW step (3C)

      INTEGER :: i_ozone_int  ! Option for interpolation of ozone

      ! Cloud:

      LOGICAL:: L_CLD_AREA           ! controls cloud area parametriz.
      LOGICAL:: L_ACF_CUSACK         ! ... whether to have Cusack
      LOGICAL:: L_ACF_BROOKS         ! ... or Brooks
      LOGICAL:: L_PC2                ! controls PC2 cloud scheme
      LOGICAL:: L_PC2_RESET          ! run PC2 scheme diagnostically
      LOGICAL:: L_PC2_LBC            ! LBC's contain cloud fractions
      LOGICAL:: L_PC2_DIAG_SH        ! Use diagnostic convective shallow cloud
                                     ! in PC2.

      ! Assimilation:

      ! Switches for assm mode
      LOGICAL:: L_AC

      ! T: Use V_INT_TP to output Temp on model levels.
      LOGICAL:: L_VINT_TP

      ! UM6.5 - MODEL_ANALYSIS_HRS replaced by MODEL_ANALYSIS_MINS - 
      !         change A_ASSIM_START_HR and A_ASSIM_END_HR to 
      !         A_ASSIM_START_MIN, A_ASSIM_END_MIN
      !         so that all three variables  have the same flexibility
      ! Time at which data assimilation starts (Hours after Basis Time)
      INTEGER :: A_ASSIM_START_MIN
      ! Time at which data assimilation  ends
      INTEGER :: A_ASSIM_END_MIN

      ! Switch for assimilation mode
      CHARACTER(LEN=5) :: A_ASSIM_MODE
      
      !---  PMSL diagnostic  ---
      ! Orographic height threshold for new pmsl calculation
      ! from ATMOS_STASH_Misc in UMUI for Calc_NPMSL routine
      REAL :: NPMSL_HEIGHT

      ! Switch for interpolated winds in lbcs
      LOGICAL :: L_int_uvw_lbc  ! .true. for advecting winds
                                ! interpolated in boundary zone
      
      !---  Tracers ---
      ! Aerosol

      LOGICAL:: L_MURK          !      :Total aerosol field
      LOGICAL:: L_MURK_ADVECT   !      :Aerosol advection
      LOGICAL:: L_MURK_SOURCE   !Bndry :Aerosol source & sink terms
      LOGICAL:: L_MURK_BDRY     !Layer :UK Mes bndry model
      LOGICAL:: L_BL_TRACER_MIX !model :Bndry layer tracer mixing
      LOGICAL :: L_MURK_RAD
      LOGICAL :: L_murk_lbc    !  Murk aerosol lbcs active

      ! For Aero_Ctl (Sulphur cycle or Soot)

      INTEGER CALL_CHEM_FREQ     !No. times chem called per atm tstep

!     Mineral dust aerosol
      LOGICAL :: L_DUST
!     Use old version of dust_uplift scheme used in CAM NWP models
      LOGICAL :: L_CAM_DUST

      ! Sulphur cycle

      LOGICAL :: L_SULPC_SO2   ! S Cycle: SO2 MMR included
      LOGICAL :: L_SULPC_DMS   ! S Cycle: DMS MMR included
      LOGICAL :: L_SULPC_OZONE ! S Cycle: Ozone included for oxidation 
                               !          of DMS and SO2
      LOGICAL :: L_SULPC_SO2_O3_NONBUFFERED ! S Cycle: SO2+O3 reaction NOT
                                            ! buffered by NH3.
      LOGICAL :: L_SULPC_ONLINE_OXIDANTS ! Sulphur Cycle : Use online
                                         ! oxidants from UKCA
      LOGICAL :: L_SO2_SURFEM  ! SO2 Surface Emissions
      LOGICAL :: L_SO2_HILEM   ! SO2 High Level Emissions
      LOGICAL :: L_SO2_NATEM   ! SO2 Natural Emissions
      LOGICAL :: L_DMS_EM      ! DMS Emissions
      LOGICAL :: L_DMS_EM_INTER      ! Interactive DMS Emissions
      LOGICAL :: L_DMS_Ointer        ! DMS emissions from ocean model
      LOGICAL :: L_DMS_Liss_Merlivat ! Switches to determine which
      LOGICAL :: L_DMS_Wanninkhof    !    scheme to use for interactive
      LOGICAL :: L_DMS_Nightingale   !    sea-air exchange of DMS
      LOGICAL :: L_SULPC_NH3   ! S Cycle: NH3 tracer included
      LOGICAL :: L_NH3_EM      ! S Cycle: NH3 emiss included

      ! Soot cycle

      LOGICAL :: L_SOOT                ! Soot included
      LOGICAL :: L_SOOT_SUREM          ! surface Soot emiss included
      LOGICAL :: L_SOOT_HILEM          ! elevated Soot emiss included

      ! Biomass aerosol

      LOGICAL :: L_BIOMASS             ! Biomass aerosol included
      LOGICAL :: L_BMASS_SUREM         ! Sfc biomass emiss included
      LOGICAL :: L_BMASS_HILEM         ! Elevated bmass emiss included
      
      ! Fossil-fuel organic carbon aerosol
      
      LOGICAL :: L_OCFF                ! OCFF aerosol included
      LOGICAL :: L_OCFF_SUREM          ! Surface OCFF emissions included
      LOGICAL :: L_OCFF_HILEM          ! Elevated OCFF emiss included
      
      ! Carbon cycle

      ! Interactive 3D CO2 field for use with carbon cycle model
      LOGICAL :: L_CO2_INTERACTIVE
      ! Switch for Radiation Interaction with CO2 - kdcorbin, 06/10
      LOGICAL :: L_CO2_RADIATION   
      ! Switch for CABLE - kdcorbin, 03/10
      LOGICAL :: l_cable
      
      ! Switch for calculating CO2/tracer atmospheric mass - kdcorbin, 05/10
      LOGICAL :: L_TRACER_MASS, L_CO2_MASS
      INTEGER :: I_TRACERMASS_START
      ! Switch for running passive tracers using CO2 fluxes - rml, 1/7/13
      LOGICAL :: L_CO2_TRACER

      ! Switch for calculating methane atmospheric loss - kdcorbin, 05/10
      LOGICAL :: L_METHANE_LOSS
      INTEGER :: I_METHANE_TRACERS

      ! Switch for calculating mcf atmospheric/ocean loss - kdcorbin, 05/10
      LOGICAL :: L_MCF_LOSS
      INTEGER :: I_MCF_TRACERNUMBER

      ! Switch for calculating radon decay - kdcorbin, 05/10
      LOGICAL :: L_RADON_DECAY
      INTEGER :: I_RADON_TRACERNUMBER

      ! CO2 Mass - kdcorbin, 05/10
      REAL :: CO2EMITMASS

     ! Tracer Mass - kdcorbin, 05/10
      REAL :: TMASS(21)

      LOGICAL :: L_CO2_EMITS          ! Include surface emissions

      ! 10m windspeed for air/sea gas flux calculations
      LOGICAL :: L_PSSWND          ! pass 10m windspeed to the ocean

      ! Dust deposition for ocean biology
      LOGICAL :: L_DUST2OCN        ! pass dust dep fields to the ocean

      LOGICAL :: L_Q10                  ! control T fn for soil resp

      ! Switch for turning off boundary layer code

      LOGICAL :: L_bl            !  F: Turns off boundary layer code
      ! MOSES II and Triffid logicals--------------------

      LOGICAL :: L_VEG_FRACS          ! Switch for vegetation fractions
      LOGICAL :: L_SNOW_ALBEDO        ! Prognostic snow albedo
      LOGICAL :: L_TRIFFID            ! Switch for interactive veg model
      LOGICAL :: L_PHENOL             ! Switch for leaf phenology

      ! Switch for running TRIFFID in equilibrium mode
      LOGICAL :: L_TRIF_EQ

      ! Switch for starting NRUN mid-way through a TRIFFID calling
      ! period
      LOGICAL :: L_NRUN_MID_TRIF
      LOGICAL :: L_DISTURB      ! Switch for anthropogenic disturbance

      INTEGER :: CAN_MODEL ! Switch for thermal vegetation canopy

      ! Vegetation:

      ! Update frequency for leaf phenology (days)
      INTEGER :: PHENOL_PERIOD

      INTEGER :: TRIFFID_PERIOD ! Update frequency for TRIFFID (days)

      ! Switch for anthropogenic heat source 
      LOGICAL :: l_anthrop_heat_src 

      ! Hardwire until re-assessment of whether these need to be
      ! re-introduced as UMUI-set switches.
! RR old switches needed for addressing but should be defunct?
! RR - can be set with parameters if needed in the interim. ->

      ! Large scale precipitation:

      LOGICAL :: L_rain          !  F: Turns off precipitation code

      ! 'New' cloud/precip microphysics, Defunct, only mixed phase
      ! phys supported
      LOGICAL, PARAMETER :: L_LSPICE    =.true.

      ! Microphysics complexity
      LOGICAL :: L_mcr_qcf2    !  Include second ice variable
      LOGICAL :: L_mcr_qrain   !  Include prognostic rain
      LOGICAL :: L_mcr_qgraup  !  Include prognosic graupel
      LOGICAL :: L_mcr_qcf2_lbc    !  Second ice variable lbcs active
      LOGICAL :: L_mcr_qrain_lbc   !  Prognostic rain lbcs active
      LOGICAL :: L_mcr_qgraup_lbc  !  Prognosic graupel lbcs active

      ! Controls the use of new RHcrit parametrization, option in Sec 9
      ! vn 2A.
      LOGICAL :: L_RHCPT

! Logicals for different radiation packages
      LOGICAL :: L_FORCING
! Calculate radiative forcings (3C)
      LOGICAL :: L_RADIANCE
! Calculate radiances          (3C)
      LOGICAL :: L_TIMESTEP
! Use new timestepping scheme  (3C)
      LOGICAL :: L_WENYI       
! Include Wenyi's pressure & temperature scaling (3A/3C)

      ! Convection:

      LOGICAL :: L_3D_CCA             ! Use 3D conv cloud amount
      LOGICAL :: L_PHASE_LIM          ! Limits phase change of precip
                                      ! in convective downdraught
      LOGICAL :: L_CCRAD              ! Main logical, will remove code
                                      ! connected with CCRad
                                      ! (including bugfixes)
      LOGICAL :: L_3D_CCW             ! Requires l_ccrad=.TRUE.
                                      ! .TRUE. : Radiation to use 3d ccw
                                      ! profile passed to it from
                                      ! convection.
                                      ! .FALSE.: Radiation constructs
                                      ! mean CCW profile from cclwp,ccb
                                      ! and cct as in original.

      ! Timestep frequency for calling convection
      ! Hardwired to calling every timestep
      INTEGER,PARAMETER :: A_CONV_STEP = 1

      ! GWD scheme:
      LOGICAL :: L_GWD        ! Use SSO drag scheme
      LOGICAL :: L_USE_USSP   ! Use spectral GWD scheme

      ! Radiation:

      ! Changes to open sea albedo for HadGEM1
      LOGICAL :: L_MOD_BARKER_ALBEDO ! Modified Barker albedo
      LOGICAL :: L_USE_SPEC_SEA      ! Spectr. dep. sea albedos

      ! Use modulus of fluxes to remove negative effective extinctions
      LOGICAL :: L_MOD_K_FLUX

      ! Fix the selection of fractional sea points in LW radiation
      LOGICAL :: L_CTILE_FIX

      ! Fix instability in quadratic correction to LW source term
      LOGICAL :: L_QUAD_SRC_FIX

      ! Scale the condensed water content to simulate
      ! inhomogeneous clouds
      LOGICAL :: L_INHOM_CLOUD

! Orography correction to SW radiation
      LOGICAL :: L_use_orog_corr    !  Find gradients from mean orog
      LOGICAL :: L_use_grad_corr    !  Use ancillary X & Y gradients

      ! Tropopause-based Ozone Scheme
      LOGICAL :: L_use_tpps_ozone   !  Use TPPS ozone scheme

! Methane oxidation
      REAL    :: Z_TOP
      LOGICAL :: L_USE_METHOX

! STOCHEM coupling to radiation
      LOGICAL :: L_USE_STOCHEM_CH4   ! for methane
      LOGICAL :: L_USE_STOCHEM_O3    ! for ozone
      ! River Routing
      LOGICAL :: L_RIVERS
      LOGICAL :: L_INLAND   ! control rerouting of inland basin water
      REAL    :: RIVER_STEP

      ! Hydrology:

      LOGICAL :: L_hydrology     !  F: Turns off hydrology code

! Max humidity in STRATQ
      LOGICAL :: L_UPPER_STRATQ

      LOGICAL, PARAMETER :: LMOSES        =.TRUE.  ! MOSES hydrology
      LOGICAL :: L_ICOUNT       !  T: Output iteration counts

      ! Mixing ratios:

      Logical :: l_mr_physics1            ! Use mixing ratio in
                                          ! atmos_physics1
      Logical :: l_mr_physics2            ! Use mixing ratio in
                                          ! atmos_physics2
! Stochastic Physics Random Parameters      
      LOGICAL :: L_RPSEED_READ  !  T: Read in previously specified seed
      LOGICAL :: L_RPSEED_WRITE !  T: WRITE out seed


! Ozone tracer as input to radiation scheme      
      LOGICAL :: L_USE_CARIOLLE
      LOGICAL :: L_USE_OZONEINRAD

! RR old switches---------------------------------------- <-

! OASIS coupling
      LOGICAL :: L_OASIS   ! OASIS coupling switch
      LOGICAL :: L_COUPLE_MASTER    ! Couple through master PE
      INTEGER :: OASIS_COUPLE_FREQ  ! Coupling frequency in
                                    ! number of timesteps. 

!     Logicals for UK Chemistry and Aerosols (UKCA) Model

      LOGICAL :: L_ukca           ! True when UKCA is switched on

! Natural climate forcing
      LOGICAL :: L_SCVARY            ! time varying solar forcing
      LOGICAL :: L_VOLCTS            ! time varying volcanic forcing

      COMMON/ CNTLCATM/                                                 &
     &  Model_domain,L_emcorr,                                          &
     &  L_OASIS, OASIS_COUPLE_FREQ, L_COUPLE_MASTER,                    &
     &  L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,                          &
     &  L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc,              &
     &  problem_number,                                                 &
     &  L_SNOW_ALBEDO,l_ssice_albedo,L_MICROPHY,H_SWBANDS,H_LWBANDS,    &
     &  A_SW_RADSTEP,A_LW_RADSTEP,                                      &
     &  A_SW_RADSTEP_DIAG,A_LW_RADSTEP_DIAG,                            &
     &  A_SW_RADSTEP_PROG,A_LW_RADSTEP_PROG,                            &
     &  L_CLD_AREA,L_ACF_CUSACK,L_ACF_BROOKS,L_PC2,L_PC2_RESET,         &
     &  L_PC2_LBC,L_PC2_diag_sh,L_rad_deg, L_AC ,L_VINT_TP,             &
     &  A_ASSIM_START_MIN, A_ASSIM_END_MIN,                             &
     &  NPMSL_HEIGHT,                                                   &
     &  LENERGY, LFLUX_RESET, LMASS_CORR , LQT_CORR, LEMQ_PRINT,        &
     &  A_ENERGYSTEPS,L_GWD,L_USE_USSP,                                 &
     &  L_Murk, L_MURK_ADVECT, L_MURK_SOURCE, L_MURK_BDRY,              &
     &  L_MURK_RAD, L_murk_lbc, L_BL_TRACER_MIX, L_int_uvw_lbc,         &
     &  L_DUST,L_CAM_DUST,L_SULPC_SO2,L_SULPC_DMS,L_SULPC_OZONE,        &
     &  L_SULPC_SO2_O3_NONBUFFERED,L_SO2_SURFEM, L_SO2_HILEM,           &
     &  L_SO2_NATEM,L_DMS_EM, L_DMS_EM_INTER,                           &
     &  L_SULPC_ONLINE_OXIDANTS,                                        &
     &  L_DMS_Ointer,                                                   &
     &  L_DMS_Liss_Merlivat, L_DMS_Wanninkhof, L_DMS_Nightingale,       &
     &  L_SULPC_NH3, L_NH3_EM, L_SOOT, L_SOOT_SUREM, L_SOOT_HILEM,      &
     &  L_BIOMASS, L_BMASS_SUREM, L_BMASS_HILEM,                        &
     &  L_PSSWND, L_DUST2OCN,                                           &
     &  L_RHCPT, L_CCRAD, L_3D_CCA, L_3D_CCW, L_PHASE_LIM,              &
     &  L_CO2_INTERACTIVE,L_CO2_EMITS,l_cable,                          &
! rml 1/7/13
     &  L_CO2_TRACER,                                                   &
!kdcorbin, 08/10
     &  L_CO2_RADIATION,L_TRACER_MASS,L_CO2_MASS,I_TRACERMASS_START,    &
     &  L_METHANE_LOSS,I_METHANE_TRACERS,L_MCF_LOSS,I_MCF_TRACERNUMBER, &
     &  L_RADON_DECAY,I_RADON_TRACERNUMBER,                             &
     &  L_Q10, L_NEG_TSTAR, L_VEG_FRACS, L_TRIFFID, L_PHENOL,           &
     &  L_TRIF_EQ, L_NRUN_MID_TRIF, L_DISTURB,                          &
     &  CAN_MODEL, PHENOL_PERIOD, TRIFFID_PERIOD,                       &
     &  L_USE_SEASALT_INDIRECT, L_USE_BIOGENIC,                         &
     &  L_USE_SEASALT_DIRECT, L_USE_DUST, L_USE_SULPC_INDIRECT_SW,      &
     &  L_USE_SULPC_INDIRECT_LW, L_USE_SULPC_DIRECT,                    &
     &  L_USE_SULPHATE_AUTOCONV, L_USE_SEASALT_AUTOCONV, L_AUTO_DEBIAS, &
     &  L_USE_SULPHATE_SULPC, L_USE_SEASALT_SULPC,                      &
     &  L_OCFF, L_OCFF_SUREM, L_OCFF_HILEM, L_USE_OCFF_AUTOCONV,        &
     &  L_USE_OCFF_SULPC, L_USE_OCFF_DIRECT, L_USE_OCFF_INDIRECT,       &
     &  L_USE_STOCHEM_CH4, L_USE_STOCHEM_O3,                            &
     &  L_MOD_BARKER_ALBEDO, L_USE_SPEC_SEA, L_MOD_K_FLUX, L_CTILE_FIX, &
     &  L_QUAD_SRC_FIX, L_USE_TPPS_OZONE, I_OZONE_INT,                  &
     &  L_ICOUNT,                                                       &
     &  L_use_orog_corr, L_use_grad_corr,                               &
     &  CALL_CHEM_FREQ, L_USE_SOOT_DIRECT, L_USE_SOOT_INDIRECT,         &
     &  L_USE_SOOT_AUTOCONV, L_USE_SOOT_SULPC, L_USE_BMASS_DIRECT,      &
     &  L_USE_BMASS_INDIRECT, L_USE_BMASS_AUTOCONV, L_USE_BMASS_SULPC,  &
     &  L_USE_ARCLBIOM, L_USE_ARCLBLCK,  L_USE_ARCLSSLT,                &
     &  L_USE_ARCLSULP, L_USE_ARCLDUST,  L_USE_ARCLOCFF, L_USE_ARCLDLTA,&
     &  L_ukca,                                                         &
     &  L_CTILE, L_RIVERS,L_INLAND, RIVER_STEP,                         &
     &  l_sice_meltponds, l_sice_scattering, l_sice_hadgem1a,           &
     &  L_USE_METHOX,Z_TOP,l_mr_physics1,l_mr_physics2,                 &
     &  L_INHOM_CLOUD,                                                  &
     &  L_radiation,L_FORCING,L_TIMESTEP,                               &
     &  L_RADIANCE, L_WENYI, L_bl, L_rain, L_hydrology,                 &
     &  L_TOP,L_PDM,L_USE_AOD,L_USE_CLEARRH,L_UPPER_STRATQ,             &
     &  L_sice_heatflux,L_SOIL_SAT_DOWN,                                &
     &  L_SCVARY,L_VOLCTS,                                              &
     &  l_anthrop_heat_src,                                             &
     &  L_USE_CARIOLLE,L_USE_OZONEINRAD,                                & 
     &  L_RPSEED_READ, L_RPSEED_WRITE,                                  &
     ! Character variables need to be at the end.
     &  H_SECT, A_ASSIM_MODE

      NAMELIST/NLSTCATM/                                                &
     &  Model_domain,L_emcorr,                                          &
     &  L_OASIS, OASIS_COUPLE_FREQ, L_COUPLE_MASTER,                    &
     &  L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,                          &
     &  L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc,              &
     &  problem_number,                                                 &
     &  L_SNOW_ALBEDO,l_ssice_albedo,L_MICROPHY,H_SWBANDS,H_LWBANDS,    &
     &  A_SW_RADSTEP,A_LW_RADSTEP,                                      &
     &  A_SW_RADSTEP_DIAG,A_LW_RADSTEP_DIAG,                            &
     &  A_SW_RADSTEP_PROG,A_LW_RADSTEP_PROG,                            &
     &  L_CLD_AREA,L_ACF_CUSACK,L_ACF_BROOKS,L_PC2,L_PC2_RESET,         &
     &  L_PC2_LBC,L_PC2_diag_sh, L_rad_deg, L_AC, L_VINT_TP,            &
     &  A_ASSIM_START_MIN, A_ASSIM_END_MIN,                             &
     &  NPMSL_HEIGHT,                                                   &
     &  LMASS_CORR , LQT_CORR, LEMQ_PRINT ,A_ENERGYSTEPS,               &
     &  L_GWD,L_USE_USSP,                                               &
     &  L_Murk, L_MURK_ADVECT, L_MURK_SOURCE, L_MURK_BDRY,              &
     &  L_MURK_RAD, L_murk_lbc, L_BL_TRACER_MIX, L_int_uvw_lbc,         &
     &  L_DUST,L_CAM_DUST,L_SULPC_SO2,L_SULPC_DMS,L_SULPC_OZONE,        &
     &  L_SULPC_SO2_O3_NONBUFFERED, L_SO2_SURFEM, L_SO2_HILEM,          &
     &  L_SO2_NATEM,                                                    &
     &  L_SULPC_ONLINE_OXIDANTS,                                        &
     &  L_DMS_EM, L_DMS_EM_INTER,                                       &
     &  L_DMS_Ointer,                                                   &
     &  L_DMS_Liss_Merlivat, L_DMS_Wanninkhof, L_DMS_Nightingale,       &
     &  L_SULPC_NH3, L_NH3_EM, L_SOOT, L_SOOT_SUREM, L_SOOT_HILEM,      &
     &  L_BIOMASS, L_BMASS_SUREM, L_BMASS_HILEM,                        &
     &  L_PSSWND, L_DUST2OCN,                                           &
     &  L_RHCPT, L_CCRAD, L_3D_CCA, L_3D_CCW, L_PHASE_LIM,              &
     &  L_CO2_INTERACTIVE,L_CO2_EMITS,l_cable,                          &
! rml 1/7/13
     &  L_CO2_TRACER,                                                   &
!kdcorbin, 08/10
     &  L_CO2_RADIATION,L_TRACER_MASS,L_CO2_MASS,I_TRACERMASS_START,    &
     &  L_METHANE_LOSS,I_METHANE_TRACERS,L_MCF_LOSS,I_MCF_TRACERNUMBER, &
     &  L_RADON_DECAY,I_RADON_TRACERNUMBER,                             &  
     &  L_Q10, L_NEG_TSTAR, L_VEG_FRACS, L_TRIFFID, L_PHENOL,           &
     &  L_TRIF_EQ, L_NRUN_MID_TRIF, L_DISTURB,                          &
     &  CAN_MODEL, PHENOL_PERIOD, TRIFFID_PERIOD,                       &
     &  L_USE_SEASALT_INDIRECT, L_USE_BIOGENIC,                         &
     &  L_USE_SEASALT_DIRECT, L_USE_DUST, L_USE_SULPC_INDIRECT_SW,      &
     &  L_USE_SULPC_INDIRECT_LW, L_USE_SULPC_DIRECT,                    &
     &  L_USE_SULPHATE_AUTOCONV, L_USE_SEASALT_AUTOCONV, L_AUTO_DEBIAS, &
     &  L_USE_SULPHATE_SULPC, L_USE_SEASALT_SULPC,                      &
     &  L_OCFF, L_OCFF_SUREM, L_OCFF_HILEM, L_USE_OCFF_AUTOCONV,        &
     &  L_USE_OCFF_SULPC, L_USE_OCFF_DIRECT, L_USE_OCFF_INDIRECT,       &
     &  L_USE_STOCHEM_CH4, L_USE_STOCHEM_O3,                            &
     &  L_MOD_BARKER_ALBEDO, L_USE_SPEC_SEA, L_MOD_K_FLUX, L_CTILE_FIX, &
     &  L_QUAD_SRC_FIX, L_USE_TPPS_OZONE, I_OZONE_INT,                  &
     &  L_ICOUNT,                                                       &
     &  L_use_orog_corr, L_use_grad_corr,                               &
     &  CALL_CHEM_FREQ, L_USE_SOOT_DIRECT, L_USE_SOOT_INDIRECT,         &
     &  L_USE_SOOT_AUTOCONV, L_USE_SOOT_SULPC, L_USE_BMASS_DIRECT,      &
     &  L_USE_BMASS_INDIRECT, L_USE_BMASS_AUTOCONV, L_USE_BMASS_SULPC,  &
     &  L_USE_ARCLBIOM, L_USE_ARCLBLCK,  L_USE_ARCLSSLT,                &
     &  L_USE_ARCLSULP, L_USE_ARCLDUST,  L_USE_ARCLOCFF, L_USE_ARCLDLTA,&
     &  L_ukca,                                                         &
     &  L_CTILE, L_RIVERS,L_INLAND, RIVER_STEP,                         &
     &  l_sice_meltponds, l_sice_scattering, l_sice_hadgem1a,           &
     &  L_USE_METHOX,Z_TOP,l_mr_physics1,l_mr_physics2,                 &
     &  L_INHOM_CLOUD,                                                  &
     &  L_radiation,L_FORCING,L_TIMESTEP,                               &
     &  L_RADIANCE, L_WENYI, L_bl, L_rain, L_hydrology,                 &
     &  L_TOP,L_PDM,L_USE_AOD,L_USE_CLEARRH,L_UPPER_STRATQ,             &
     &  L_SCVARY,L_VOLCTS,                                              &
     &  L_sice_heatflux, L_SOIL_SAT_DOWN,                               &
     &  l_anthrop_heat_src,                                             &
     &  L_USE_CARIOLLE,L_USE_OZONEINRAD,                                &
     &  L_RPSEED_READ, L_RPSEED_WRITE,                                  &
     &  H_SECT, A_ASSIM_MODE

      ! Control switches derived and set within the model, hence not
      ! passed in from the UMUI via namelist:

      ! Radiation
      LOGICAL :: L_SW_RADIATE  ! Activate SW radiation this timestep
      LOGICAL :: L_LW_RADIATE  ! Activate LW radiation this timestep
      LOGICAL :: L_SW_RADIATE_DIAG
! Activate fast SW radiation this timestep (3C)
      LOGICAL :: L_LW_RADIATE_DIAG
! Activate fast LW radiation this timestep (3C)
      LOGICAL :: L_SW_RADIATE_PROG
! Activate slow SW radiation this timestep (3C)
      LOGICAL :: L_LW_RADIATE_PROG
! Activate slow LW radiation this timestep (3C)
      LOGICAL :: Lexpand_ozone ! Convert zonal mean ozone to field

      ! convert zonal mean tpps ozone to field
      LOGICAL :: Lexpand_tpps_ozone

      INTEGER :: I_tpps_ozone_opts ! options for tropopause-based ozone

      ! size of super array holding all tracers
      Integer ::   super_array_size
      Integer ::   moisture_array_size

      COMMON/CNTLCATM2/                                                 &
       !  super tracer array size
     &  super_array_size, moisture_array_size,                          &
     &  L_SW_RADIATE,L_LW_RADIATE,                                      &
     &  L_SW_RADIATE_DIAG,L_LW_RADIATE_DIAG,                            &
     &  L_SW_RADIATE_PROG,L_LW_RADIATE_PROG,                            &
     &  Lexpand_ozone,                                                  &
     & Lexpand_tpps_ozone, I_tpps_ozone_opts
! options for tropopause-based ozone



      LOGICAL :: LTLEADS ! Switch for Leads temperature.
!                      ! If FALSE, they are assumed to be TFS
!                      ! Else they are prognostic.
! Default setting: leads temperatures are set to TFS
! HARDWIRE LTLEADS
      PARAMETER(LTLEADS =.FALSE.)
! SWOPT3A defines SW options for the radiation code.
!
!     NOTE: SWOPT3A and SWCAVR3A must be consistent.
!
!     Current Owner of Code: J. M. Edwards
!
!     History:
!     Version  Date      Comment.
!     5.2      05/12/00  Logical for extra top level added.
!                                J. M. Edwards
! 6.2  24/10/05 Functionality for radiative forcing, timestepping
!               and radiances under versions 3C and 3Z of radiation
!               code added                 (J.-C. Thelen)
!
!
      INTEGER :: I_2STREAM_SW     ! two_stream scheme selected
      INTEGER :: I_GAS_OVERLAP_SW ! treatment of gaseous overlaps
      INTEGER :: I_CLOUD_SW       ! treatment of cloud overlaps
      INTEGER :: I_CLOUD_REPRESENTATION_SW  ! representation of clouds
      INTEGER :: I_SOLVER_SW                ! solver selected

      LOGICAL :: L_O2_SW ! flag for oxygen

      ! types of droplets or ice crystals used for parametrizations
      INTEGER :: I_ST_WATER_SW  ! type for stratiform water
      INTEGER :: I_CNV_WATER_SW ! type for convective water
      INTEGER :: I_ST_ICE_SW    ! type for stratiform ice
      INTEGER :: I_CNV_ICE_SW   ! type for convective ice

      ! flag to partition convective clouds using the local temperature
      LOGICAL :: L_LOCAL_CNV_PARTITION_SW


      ! Flag to include an extra layer at the top of the atmosphere in
      ! radiative calculations
      LOGICAL :: L_EXTRA_TOP_SW

! SWOPT3A end
!     ------------------------------------------------------------------
!     COMMON BLOCK CONTAINING OPTIONS FOR 3A-RADIATION CODE.
!     VERSION FOR SHORTWAVE CALCULATIONS.
!
      COMMON/R2SWOPT/                                                   &
     &  I_2STREAM_SW, I_GAS_OVERLAP_SW, I_CLOUD_SW,                     &
     &  I_CLOUD_REPRESENTATION_SW, I_SOLVER_SW,                         &
     &  L_O2_SW,                                                        &
     &  I_ST_WATER_SW, I_CNV_WATER_SW, I_ST_ICE_SW, I_CNV_ICE_SW,       &
     &  L_LOCAL_CNV_PARTITION_SW, L_EXTRA_TOP_SW
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     Module defining LW options for the radiation code.
!
!     NOTE: LWOPT3A and LWCAVR3A must be consistent.
!
!     Current Owner of Code: J. M. Edwards
!
!     History:
!     Version  Date      Comment.
!     5.2      05/12/00  Logical for extra top level added.
!                                J. M. Edwards
!     5.3      04/10/01  Option for LW scattering added.
!                                J. M. Edwards
! 6.2  24/10/05 Functionality for radiative forcing, timestepping
!               and radiances under versions 3C and 3Z of radiation
!               code added                 (J.-C. Thelen)
!
      INTEGER                                                           &
     &     I_2STREAM_LW                                                 &
!             TWO-STREAM SCHEME
     &   , I_GAS_OVERLAP_LW                                             &
!             TREATMENT OF GASEOUS OVERLAPS
     &   , I_CLOUD_LW                                                   &
!             TREATMENT OF CLOUDY OVERLAPS
     &   , I_CLOUD_REPRESENTATION_LW                                    &
!             REPRESENTATION OF CLOUDS
     &   , I_SOLVER_LW                                                  &
!             SOLVER SELECTED
     &   , I_SCATTER_METHOD_LW
!             TREATMENT OF SCATTERING (THIS CONTROLS THE CALCULATION
!             OF OPTICAL PROPERTIES)
      LOGICAL                                                           &
     &     L_IR_SOURCE_QUAD_LW                                          &
!             REPRESENTATION OF THE IR-SOURCE TERM
     &   , L_MICROPHYSICS_LW                                            &
!             FLAG FOR MICROPHYSICS IN LONG WAVE
     &   , L_LOCAL_CNV_PARTITION_LW                                     &
!             FLAG TO PARTITION CONVECTIVE CLOUD USING THE
!             LOCAL TEMPERATURE.
     &   , L_EXTRA_TOP_LW
!             Flag to include an extra layer at the top of the
!             atmosphere in radiative calculations
!     OPTIONS FOR TRACE GASES:
      LOGICAL                                                           &
     &     L_N2O_LW                                                     &
!             FLAG FOR NITROUS OXIDE
     &   , L_CH4_LW                                                     &
!             FLAG FOR METHANE
     &   , L_CFC11_LW                                                   &
!             FLAG FOR CFC11
     &   , L_CFC12_LW                                                   &
!             FLAG FOR CFC12
     &   , L_CFC113_LW                                                  &
!             FLAG FOR CFC113
     &   , L_HCFC22_LW                                                  &
!             FLAG FOR HCFC22
     &   , L_HFC125_LW                                                  &
!             FLAG FOR HFC125
     &   , L_HFC134A_LW
!             FLAG FOR HFC134A
!
!     TYPES OF DROPLETS OR ICE CRYSTALS USED FOR PARAMETRIZATIONS
      INTEGER                                                           &
     &     I_ST_WATER_LW                                                &
!             TYPE FOR STRATIFORM WATER
     &   , I_CNV_WATER_LW                                               &
!             TYPE FOR CONVECTIVE WATER
     &   , I_ST_ICE_LW                                                  &
!             TYPE FOR STRATIFORM ICE
     &   , I_CNV_ICE_LW
!             TYPE FOR CONVECTIVE ICE
!
!     ------------------------------------------------------------------
! LWCOPT3A contains options for 3a-radiation code version for
! longwave calculations.

      COMMON/R2LWOPT/                                                   &
     &  I_2STREAM_LW, L_IR_SOURCE_QUAD_LW, I_GAS_OVERLAP_LW,            &
     &  I_CLOUD_LW, I_CLOUD_REPRESENTATION_LW, I_SOLVER_LW,             &
     &  I_SCATTER_METHOD_LW,                                            &
     &  L_N2O_LW, L_CH4_LW, L_CFC11_LW , L_CFC12_LW, L_CFC113_LW,       &
     &  L_HCFC22_LW, L_HFC125_LW, L_HFC134A_LW,                         &
     &  I_ST_WATER_LW, I_CNV_WATER_LW, I_ST_ICE_LW,                     &
     &  I_CNV_ICE_LW, L_MICROPHYSICS_LW,                                &
     &  L_LOCAL_CNV_PARTITION_LW, L_EXTRA_TOP_LW
! LWCOPT3A end

!
!
!     ------------- Shortwave Radiation -------------------------
!
      IF (H_SECT(1) == '03A') THEN
!
! DEPENDS ON: r2_sw_specin
        CALL R2_SW_SPECIN(ICODE, CMESSAGE                               &
     &    , L_O2_SW                                                     &
     &    , L_CLIMAT_AEROSOL                                            &
     &    , L_USE_DUST, L_USE_ARCLDUST                                  &
     &    , L_USE_SULPC_DIRECT, L_USE_ARCLSULP                          &
     &    , L_USE_SOOT_DIRECT, L_USE_ARCLBLCK                           &
     &    , L_USE_BMASS_DIRECT, L_USE_ARCLBIOM                          &
     &    , L_USE_SEASALT_DIRECT, L_USE_ARCLSSLT                        &
     &    , L_USE_OCFF_DIRECT, L_USE_ARCLOCFF                           &
     &    , L_USE_BIOGENIC, L_USE_ARCLDLTA                              &
     &    , L_MURK_RAD                                                  &
     &   )
        IF (ICODE /= 0) RETURN
!
      ELSE
!
        ICODE=1
        CMESSAGE='Unknown version of SW radiation encountered.'
        RETURN
!
      ENDIF
!
!
!     ------------- Longwave Radiation -------------------------
!
      IF (H_SECT(2) == '03A') THEN
!
! DEPENDS ON: r2_lw_specin
        CALL R2_LW_SPECIN(ICODE, CMESSAGE                               &
     &    , L_CH4_LW, L_N2O_LW, L_CFC11_LW, L_CFC12_LW                  &
     &    , L_CFC113_LW, L_HCFC22_LW, L_HFC125_LW, L_HFC134A_LW         &
     &    , L_CLIMAT_AEROSOL                                            &
     &    , L_USE_DUST, L_USE_ARCLDUST                                  &
     &    , L_USE_SULPC_DIRECT, L_USE_ARCLSULP                          &
     &    , L_USE_SOOT_DIRECT, L_USE_ARCLBLCK                           &
     &    , L_USE_BMASS_DIRECT, L_USE_ARCLBIOM                          &
     &    , L_USE_SEASALT_DIRECT, L_USE_ARCLSSLT                        &
     &    , L_USE_OCFF_DIRECT, L_USE_ARCLOCFF                           &
     &    , L_USE_BIOGENIC, L_USE_ARCLDLTA                              &
     &    , L_MURK_RAD                                                  &
     &    , L_USE_AOD                                                   &
     &   )
!
      ELSE
!
        ICODE=1
        CMESSAGE='Unknown version of LW radiation encountered.'
        RETURN
!
      ENDIF
!
!
!
!
      RETURN
      END SUBROUTINE INITPHYS
