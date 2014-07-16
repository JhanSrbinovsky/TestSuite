
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!*LL  SUBROUTINES LS_PPN and LS_PPNC------------------------------------
!LL  Purpose:
!LL          LS_PPN and LS_PPNC:
!LL           Calculate large-scale (dynamical) precipitation.
!LL           LS_PPNC is the gather/scatter routine which then
!LL           calls LSP_ICE.
!LL  Note: in all cases, level counters (incl subscripts) run from 1
!LL        (lowest model layer) to Q_LEVELS (topmost "wet" model
!LL        layer) - it is assumed that the bottom Q_LEVELS layers are
!LL        the "wet" layers.
!LL
!LL  Put through fpp on Cray.  Activate *IF definition CRAY if running
!LL  on the Cray.
!LL
!LL  Modification History from Version 4.4
!LL     Version    Date
!LL       4.5      March 98          New Deck        Damian Wilson
!
!         5.1      16-12-99   Change to allow 3D RHcrit diagnostic
!                             in place of 1D parameter. AC Bushell
!         5.2      25/10/00   Reintroduction of tracers. P.Selwood
!         5.2      21-11-00   Allow interfacting with 3C scheme.
!                             Return diagnostics. Damian Wilson
!LL
!         5.2      28-11-00   Pass down sea-salt arrays.
!                                                       A. Jones
!         5.3      29-08-01   Introduce logical to control effect of
!                             sulphate aerosol on autoconversion.
!                                                       A. Jones
!         5.3      24-09-01   Tidy up code (mostly redundant) relating
!                             to the sulphur cycle.     A. Jones
!         5.4      30-08-02   Remove lscav_agedsoot, which is now
!                             calculated in microphys_ctl.  P. Davison
!         5.4      27-07-02   Change cloud fractions to prognostics
!                             for PC2                   D. Wilson
!
!         5.4      25-04-02   Pass land fraction down to LSP_ICE.
!                                                       A. Jones
!         5.5      03-02-03   Pass extra microphysics variables
!                             down to LSP_ICE.           R.M.Forbes
!         5.5      03-02-03   Include extra microphysics variables
!                             (qcf2,qrain,qgraup)          R.M.Forbes
!         6.1      01-08-04   Include variables for prognostic rain
!                             and graupel.             R.M.Forbes
!         6.1      07-04-04   Add biomass smoke to call to LSP_ICE.
!                                                          A. Jones
!         6.1      07-04-04   Pass switch for autoconversion de-biasing
!                             down to LSP_ICE.              A. Jones
!         6.2      22-08-05   Remove commented out code. P.Selwood.
!         6.2      23-11-05   Pass through precip variables from
!                             UMUI. Damian Wilson
!         6.2      11-01-06   Include non-hydrostatic calculation
!                             of air density and model thickness.
!                                                       D. Wilson
!         6.2      31-01-06   Pass process rate diags through. R.Forbes
!         6.4      10-01-07   Provide mixing ratio logical control
!                                                       D. Wilson
!LL  Programming standard: Unified Model Documentation Paper No 3
!LL
!LL  Documentation: UM Documentation Paper 26.
!LL
!*L  Arguments:---------------------------------------------------------
      SUBROUTINE LS_PPN(                                                &
     &halo_i, halo_j, off_x, off_y,                                     &
     &p_layer_boundaries,p_theta_levels,TIMESTEP,BLAND,                 &
     &CW_SEA,CW_LAND,                                                   &
     &L_seq_mcr,L_autoc_3b,L_autolim_3b,L_autoconv_murk,                &
     &ec_auto,N_drop_land,                                              &
     &N_drop_sea,N_drop_land_cr,N_drop_sea_cr,Ntot_land, Ntot_sea,      &
     &x1i,x1ic,x1r,x2r,x4r,l_psd,ai,bi,aic,bic,                         &
     &lsp_ei,lsp_fi,lsp_eic,lsp_fic,                                    &
     &CF,CFL,CFF,                                                       &
     &RHCRIT,                                                           &
     &Q_LEVELS, model_levels, bl_levels                                 &
     &,lspice_dim1,lspice_dim2,lspice_dim3                              &
     &,rho_r2,r_rho_levels, r_theta_levels, q, qcf, qcl, T,             &
     & qcf2, qrain, qgraup,                                             &
     & L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_it_melting,             &
     &L_USE_SULPHATE_AUTOCONV,                                          &
     &L_AUTO_DEBIAS, l_mixing_ratio, l_cry_agg_dep, l_droplet_settle,   &
     &SEA_SALT_FILM, SEA_SALT_JET,                                      &
     &L_SEASALT_CCN, salt_dim1, salt_dim2, salt_dim3,                   &
     &L_use_biogenic, biogenic,                                         &
     &SNOW_DEPTH, LAND_FRACT,                                           &
     &SO2,L_SULPC_SO2,                                                  &
     &NH3,L_SULPC_NH3,                                                  &
     &SO4_AIT,SO4_ACC,SO4_DIS,                                          &
     &AGED_SOOT, L_SOOT,                                                &
     &BMASS_AGD,BMASS_CLD,L_BIOMASS_CCN,                                &
     &OCFF_AGD,OCFF_CLD,L_OCFF_CCN,                                     &
     &AEROSOL,L_MURK, L_pc2,                                            &
     &LSRAIN,LSSNOW,                                                    &
     &LSCAV_SO2,LSCAV_NH3,LSCAV_SO4AIT,LSCAV_SO4ACC,LSCAV_SO4DIS,       &
     &LSRAIN3D,LSSNOW3D,RAINFRAC3D,                                     &
     &row_length,rows,                                                  &
     &rhc_row_length, rhc_rows,                                         &
      ! Process rate diagnostics
     &  PSDEP,PSAUT,PSACW,PSACR,PSACI,PSMLT,PSMLTEVP                    &
     &, PRAUT,PRACW,PREVP                                               &
     &, PGAUT,PGACW,PGACS,PGMLT                                         &
     &, PIFRW,PIPRM,PIDEP,PIACW,PIACR,PIMLT,PIMLTEVP                    &
     &, PIFALL,PSFALL,PRFALL,PGFALL,PLSET,PLEVPSET                      &
      ! Process rate diagnostic switches
     &, L_PSDEP_diag,L_PSAUT_diag,L_PSACW_diag,L_PSACR_diag             &
     &, L_PSACI_diag,L_PSMLT_diag,L_PSMLTEVP_diag                       &
     &, L_PRAUT_diag,L_PRACW_diag,L_PREVP_diag                          &
     &, L_PGAUT_diag,L_PGACW_diag,L_PGACS_diag,L_PGMLT_diag             &
     &, L_PIFRW_diag,L_PIPRM_diag,L_PIDEP_diag,L_PIACW_diag             &
     &, L_PIACR_diag,L_PIMLT_diag,L_PIMLTEVP_diag                       &
     &, L_PIFALL_diag,L_PSFALL_diag,L_PRFALL_diag,L_PGFALL_diag         &
     &, L_PLSET_diag,L_PLEVPSET_diag                                    &
      ! Variables for stochastic physics random parameters2
     &, M_CI,                                                           &
     &maxsects, h_sect, ERROR                                           &
     &)

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::                                            &
     &  halo_i                                                          &
                 ! Halo size in x direction
     &, halo_j                                                          &
                 ! Halo size in y direction
     &, off_x                                                           &
                 ! Small halo size in x direction
     &, off_y                                                           &
                 ! Small halo size in y direction
     &, Q_LEVELS                                                        &
                 ! Number of "wet" levels in the model.
     &, model_levels                                                    &
                     ! Number of model levels
     &, bl_levels                                                       &
                     ! Number of boundary layer levels
     &, row_length,rows                                                 &
     &, rhc_row_length,rhc_rows                                         &
     &, lspice_dim1,lspice_dim2,lspice_dim3                             &
! Dimensions for 3D diagnostic arrays
     &, salt_dim1                                                       &
                     ! Array dimensions for sea-salt arrays (equal
     &, salt_dim2                                                       &
                     ! either to row_length, rows and Q_LEVELS, or
     &, salt_dim3    ! else 1,1,1, depending on L_SEASALT_CCN).

       INTEGER                                                          &
     &  salt_dim_ice                                                    &
                     ! Array dimension for passing to LSP_ICE.
     &, sea_salt_ptr                                                    &
                     ! Pointer for sea-salt arrays.
     &, biog_dim_ice                                                    &
                     ! Array dimension for passing to LSP_ICE.
     &, biogenic_ptr ! Pointer for biogenic aerosol array.
      REAL                                                              &
     & CF(row_length,rows,Q_LEVELS)                                     &
                                    ! IN Cloud fraction.
     &,p_theta_levels(row_length,rows,Q_LEVELS)                         &
     &,p_layer_boundaries(row_length, rows, 0:Q_LEVELS)                 &
     &,RHCRIT(rhc_row_length,rhc_rows,Q_LEVELS)                         &
                                                 ! IN Critical humidity
!                                                  for cloud formation.
     &,CFL(row_length, rows,Q_LEVELS)                                   &
                                        !IN Cloud liquid fraction.
     &,CFF(row_length, rows,Q_LEVELS)                                   &
                                        !IN Cloud ice fraction.
     &,rho_r2(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels) &
                                ! IN Air density * earth radius**2
     &,r_rho_levels(1-halo_i:row_length+halo_i,                         &
     &              1-halo_j:rows+halo_j,model_levels)                  &
                                ! IN Earths radius at each rho level
     &,r_theta_levels(1-halo_i:row_length+halo_i,                       &
     &              1-halo_j:rows+halo_j,0:model_levels)
                                ! IN Earths radius at each theta level

      REAL TIMESTEP                                                     &
                           ! IN Timestep (sec).
     &    ,CW_SEA                                                       &
                           ! IN threshold cloud liquid water content
!                               over sea for conversion to ppn
!                               (kg water per m**3)
     &    ,CW_LAND         ! IN threshold cloud liquid water content
!                               over land for conversion to ppn
!                               (kg water per m**3)
      LOGICAL BLAND(row_length,rows)                                    &
                                     ! IN Land/sea mask
     &,       L_MURK        ! IN Aerosol needs scavenging.

      LOGICAL L_USE_SULPHATE_AUTOCONV                                   &
                                       ! IN   True if sulphate aerosol
!                                      !   used in autoconversion (i.e.
!                                      !   second indirect effect)
     &,       L_SEASALT_CCN                                             &
                                       ! IN   True if sea-salt required
!                                      !   for second indirect effect.
     &,       L_BIOMASS_CCN                                             &
                                       ! IN   True if biomass smoke
!                                      !   used for 2nd indirect effect.
     &,       L_OCFF_CCN                                                &
                                       ! IN   True if fossil-fuel
                                       !  organic carbon aerosol is used
                                       !  for 2nd indirect effect
     &,       L_use_biogenic                                            &
                                       ! IN   True if biogenic aerosol
!                                      !   used for 2nd indirect effect.
     &,       L_AUTO_DEBIAS            ! IN   True if autoconversion
!                                      !   de-biasing scheme selected.

      LOGICAL L_SULPC_SO2  !IN Sulphur Cycle on, tracers scavenged if T
      LOGICAL L_SULPC_NH3  !IN Sul. Cycle NH3 on, tracers scav'd if T
      LOGICAL L_SOOT       !IN Soot Cycle on, tracers scavenged if T
      LOGICAL L_mcr_qcf2   !IN Use second prognostic ice if T
      LOGICAL L_mcr_qrain  !IN Use prognostic rain if T
      LOGICAL L_mcr_qgraup !IN Use prognostic graupel if T
      LOGICAL L_it_melting !IN Use iterative melting
      LOGICAL L_pc2        !IN Use the PC2 cloud and condenstn scheme
      Logical L_psd        !In Use generic ice particle size distributn

      REAL, INTENT(INOUT) ::                                            &
     & Q(row_length,rows,Q_LEVELS)                                      &
                                         ! Specific humidity (kg water
     &,QCF(row_length,rows,Q_LEVELS)                                    &
                                         ! Cloud ice (kg per kg air).
     &,QCL(row_length,rows,Q_LEVELS)                                    &
                                         ! Cloud liquid water (kg per
     &,qcf2(row_length,rows,q_levels)                                   &
                                         ! Ice (kg per kg air)
     &,qrain(row_length,rows,q_levels)                                  &
                                         ! Rain (kg per kg air)
     &,qgraup(row_length,rows,q_levels)                                 &
                                         ! Graupel (kg per kg air)
     &,T(row_length,rows,Q_LEVELS)                                      &
                                         ! Temperature (K).
     &,AEROSOL(row_length,rows,Q_LEVELS) ! Aerosol (K).

      REAL,INTENT(INOUT) ::                                             &
                                     !Sulphur Cycle tracers (mmr kg/kg)
     &    SO2(row_length,rows,Q_LEVELS)                                 &
     &   ,NH3(row_length,rows,Q_LEVELS)                                 &
     &   ,SO4_AIT(row_length,rows,Q_LEVELS)                             &
     &   ,SO4_ACC(row_length,rows,Q_LEVELS)                             &
     &   ,SO4_DIS(row_length,rows,Q_LEVELS)

      REAL,INTENT(INOUT) ::                                             &
                                     !Biomass smoke tracers
     &    BMASS_AGD(row_length, rows, Q_LEVELS)                         &
     &   ,BMASS_CLD(row_length, rows, Q_LEVELS)
      
      REAL,INTENT(INOUT) ::                                             &
                                     !Fossil-fuel organic carbon tracers
     &    OCFF_AGD(row_length, rows, Q_LEVELS)                          &
     &   ,OCFF_CLD(row_length, rows, Q_LEVELS)
      
      REAL,INTENT(INOUT) ::                                             &
                                     !Soot cycle tracers
     &    AGED_SOOT(row_length, rows, Q_LEVELS)
!
      REAL                                                              &
     &   SNOW_DEPTH(row_length,rows)                                    &
                                     ! IN Snow depth (m)
     &  ,LAND_FRACT(row_length,rows) ! IN Land fraction
!
      REAL                                                              &
     &    SEA_SALT_FILM(salt_dim1, salt_dim2, salt_dim3)                &
                                                         ! (m-3)
     &   ,SEA_SALT_JET(salt_dim1, salt_dim2, salt_dim3)  ! (m-3)
!
      REAL                                                              &
     &    biogenic(row_length, rows, Q_LEVELS)           ! (m.m.r.)
!
      Logical                                                           &
                   !, Intent(IN)
     &    L_seq_mcr                                                     &
                             ! Use sequential updating of mphys
     &,   L_autoc_3b                                                    &
                             ! Use 3B autoconversion method
     &,   L_autolim_3b                                                  &
                             ! Use fixed 3B values for the
                             ! autoconversion limit 
     &,   L_mixing_ratio                                                &
                             ! Use mixing ratios
     &,   L_cry_agg_dep                                                 &
                             ! Limit the supersaturation that can be
                             ! removed by deposition depending on
                             ! amount of ice in each category
     &,   l_droplet_settle                                              &
                             ! Allow cloud droplets to settle
     &,   L_autoconv_murk    ! Use murk aerosol to calc. drop number

      Real                                                              &
                   !, Intent(IN)
     &    ec_auto                                                       &
                             ! Collision coalescence efficiency
     &   ,N_drop_land                                                   &
                             ! Number of droplets over land / m-3
     &   ,N_drop_sea                                                    &
                             ! Number of droplets over sea / m-3
     &   ,N_drop_land_cr                                                &
                             ! N_drop_land ^ (-1/3) / m
     &   ,N_drop_sea_cr                                                 &
                             ! N_drop_sea ^ (-1/3) / m
     &   ,Ntot_land                                                     &
                             ! Number of droplets over land / m-3
     &   ,Ntot_sea           ! Number of droplets over sea / m-3

      Real                                                              &
                  !, Intent(IN)
     &    x1i                                                           &
                    ! Intercept of aggregate size distribution
     &   ,x1ic                                                          &
                    ! Intercept of crystal size distribution
     &   ,x1r                                                           &
                    ! Intercept of raindrop size distribution
     &   ,x2r                                                           &
                    ! Scaling parameter of raindrop size distribn
     &   ,x4r                                                           &
                    ! Shape parameter of raindrop size distribution
     &   ,ai, bi                                                        &
                    ! Ice aggregate mass size relationship m(D)=ai D^bi
     &   ,aic, bic                                                      &
                    ! Ice crystal mass size relationship m(D)=aic D^bic

     &   ,lsp_ei,  lsp_fi                                               &
                    ! Ice aggregate Best number and Reynolds number 
                    ! relationship: Re(D) =LSP_EI Be^LSP_FI
     &   ,lsp_eic, lsp_fic                                                      
                    ! Ice crystal Best number and Reynolds number 
                    ! relationship: Re(D) =LSP_EI Be^LSP_FI


! Microphysical process rate diagnostics (2D arrays on one level)
! Note: These arrays will only increase memory usage and are
!       only referenced if the particular diagnostic is active

      REAL, Intent(InOut) ::                                            &
     &  PSDEP(row_length,rows,Q_LEVELS)                                 &
                                          ! Deposition of vapour to snow
     &, PSAUT(row_length,rows,Q_LEVELS)                                 &
                                          ! Autoconversion of snow
     &, PSACW(row_length,rows,Q_LEVELS)                                 &
                                          ! Accretion of liq. by snow
     &, PSACR(row_length,rows,Q_LEVELS)                                 &
                                          ! Collection of rain by snow
     &, PSACI(row_length,rows,Q_LEVELS)                                 &
                                          ! Collection of ice crystals
     &, PSMLT(row_length,rows,Q_LEVELS)                                 &
                                          ! Melting of snow aggregates
     &, PSMLTEVP(row_length,rows,Q_LEVELS)! Evap. of melting aggregates
      REAL, Intent(InOut) ::                                            &
     &  PRAUT(row_length,rows,Q_LEVELS)                                 &
                                          ! Autoconversion of cloud
     &, PRACW(row_length,rows,Q_LEVELS)                                 &
                                          ! Accretion of liq. by rain
     &, PREVP(row_length,rows,Q_LEVELS)   ! Evaporation of rain
      REAL, Intent(InOut) ::                                            &
     &  PGAUT(row_length,rows,Q_LEVELS)                                 &
                                          ! Autoconversion of graupel
     &, PGACW(row_length,rows,Q_LEVELS)                                 &
                                          ! Accretion of liq. by graup
     &, PGACS(row_length,rows,Q_LEVELS)                                 &
                                          ! Collection of snow by graup
     &, PGMLT(row_length,rows,Q_LEVELS)   ! Melting of graupel
      REAL, Intent(InOut) ::                                            &
     &  PIFRW(row_length,rows,Q_LEVELS)                                 &
                                          ! Homogeneous freezing nucl.
     &, PIPRM(row_length,rows,Q_LEVELS)                                 &
                                          ! Heterogeneous nucl.
     &, PIDEP(row_length,rows,Q_LEVELS)                                 &
                                          ! Deposition of vapour to ice
     &, PIACW(row_length,rows,Q_LEVELS)                                 &
                                          ! Accretion of liq. by ice
     &, PIACR(row_length,rows,Q_LEVELS)                                 &
                                          ! Collection of rain by ice
     &, PIMLT(row_length,rows,Q_LEVELS)                                 &
                                          ! Melting of ice crystals
     &, PIMLTEVP(row_length,rows,Q_LEVELS)! Evaporation of melting ice
      REAL, Intent(InOut) ::                                            &
     &  PIFALL(row_length,rows,Q_LEVELS)                                &
                                          ! Sedimentation of ice crystal
     &, PSFALL(row_length,rows,Q_LEVELS)                                &
                                          ! Sedimentation of aggregates
     &, PRFALL(row_length,rows,Q_LEVELS)                                &
                                          ! Sedimentation of rain
     &, PGFALL(row_length,rows,Q_LEVELS)  ! Sedimentation of graupel
      REAL, Intent(InOut) ::                                            &
     &  PLSET(row_length,rows,Q_LEVELS)                                 &
                                          ! Droplet settling of liquid
     &, PLEVPSET(row_length,rows,Q_LEVELS)! Evaporated settled droplets

! Microphysical process rate diagnostic logical switches
      LOGICAL, Intent(In) ::                                            &
     &  L_PSDEP_diag                                                    &
                       ! Deposition of vapour to snow agg.
     &, L_PSAUT_diag                                                    &
                       ! Autoconversion of aggregates from cry
     &, L_PSACW_diag                                                    &
                       ! Accretion of liq. water by snow agg.
     &, L_PSACR_diag                                                    &
                       ! Collection of rain by snow aggregates
     &, L_PSACI_diag                                                    &
                       ! Collection of ice crystals by agg.
     &, L_PSMLT_diag                                                    &
                       ! Melting of snow aggregates
     &, L_PSMLTEVP_diag! Evaporation of melting aggregates
      LOGICAL, Intent(In) ::                                            &
     &  L_PRAUT_diag                                                    &
                       ! Autoconversion of cloud drops to rain
     &, L_PRACW_diag                                                    &
                       ! Accretion of liq. water by rain
     &, L_PREVP_diag   ! Evaporation of rain
      LOGICAL, Intent(In) ::                                            &
     &  L_PGAUT_diag                                                    &
                       ! Autoconversion of graupel from agg.
     &, L_PGACW_diag                                                    &
                       ! Accretion of liq. water by graupel
     &, L_PGACS_diag                                                    &
                       ! Collection of snow agg. by graupel
     &, L_PGMLT_diag   ! Melting of graupel
      LOGICAL, Intent(In) ::                                            &
     &  L_PIFRW_diag                                                    &
                       ! Homogeneous freezing nucleation
     &, L_PIPRM_diag                                                    &
                       ! Heterogeneous (primary) nucleation
     &, L_PIDEP_diag                                                    &
                       ! Deposition of vapour to ice crystals
     &, L_PIACW_diag                                                    &
                       ! Accretion of liq. water by ice cry.
     &, L_PIACR_diag                                                    &
                       ! Collection of rain by ice crystals
     &, L_PIMLT_diag                                                    &
                       ! Melting of ice crystals
     &, L_PIMLTEVP_diag! Evaporation of melting ice crystals
      LOGICAL, Intent(In) ::                                            &
     &  L_PIFALL_diag                                                   &
                       ! Sedimentation of ice crystals
     &, L_PSFALL_diag                                                   &
                       ! Sedimentation of aggregates
     &, L_PRFALL_diag                                                   &
                       ! Sedimentation of rain
     &, L_PGFALL_diag  ! Sedimentation of graupel
      LOGICAL, Intent(In) ::                                            &
     &  L_PLSET_diag                                                    &
                       ! Droplet settling of liquid water
     &, L_PLEVPSET_diag! Evaporated settled droplets

      REAL                                                              &
     & LSRAIN(row_length,rows)                                          &
                               ! OUT Surface rainfall rate (kg / sq m /
     &,LSSNOW(row_length,rows)                                          &
                               ! OUT Surface snowfall rate (kg / sq m /
     &,LSGRAUP(row_length,rows) ! Graupel fall rate (kg/m2/s)

      REAL                                                              &
                         ! OUT column totals of S Cycle tracers scav
     &    LSCAV_SO2(row_length,rows)                                    &
     &   ,LSCAV_NH3(row_length,rows)                                    &
     &   ,LSCAV_SO4AIT(row_length,rows)                                 &
     &   ,LSCAV_SO4ACC(row_length,rows)                                 &
     &   ,LSCAV_SO4DIS(row_length,rows)

      REAL                                                              &
     &    LSRAIN3D(lspice_dim1,lspice_dim2,lspice_dim3)                 &
                                                        ! OUT
!                           Rain rate out of each model layer
     &   ,LSSNOW3D(lspice_dim1,lspice_dim2,lspice_dim3)                 &
                                                        ! OUT
!                           Snow rate out of each model layer
     &   ,RAINFRAC3D(lspice_dim1,lspice_dim2,lspice_dim3) ! OUT
!                           rain fraction out of each model layer
! Variables for stochastic physics random parameters
      REAL,    INTENT(IN) :: M_CI   ! used to modify ice fall speed. 
!
      INTEGER                                                           &
     & ERROR          ! OUT Return code - 0 if OK,
!                                         1 if bad arguments.
      Integer                                                           &
     &   maxsects

      Character*3                                                       &
     &   h_sect(0:maxsects)

!*L  Workspace usage ---------------------------------------------------
!
      LOGICAL L_non_hydrostatic ! Use non-hydrostatic formulation of
                                ! layer thicknesses

!ajm      LOGICAL
!ajm     & H(PFIELD)      ! Used as "logical" in compression.
      INTEGER                                                           &
     & IX(row_length*rows,2)                                            &
                                 ! Index for compress/expand.
     &, n_iterations             ! Number of iterations

      REAL VFALL(row_length,rows)        ! snow fall velocity (m per s).
      REAL VFALL2(row_length,rows)       ! fall velocity for qcf2 (m/s)
      REAL LSSNOW2(row_length,rows)      ! snowfall rate for qcf2
      REAL droplet_flux(row_length,rows) ! water drop flux / kg m-2 s-1
      REAL VFALL_RAIN(row_length,rows)   ! fall velocity for rain (m/s)
      REAL VFALL_GRAUP(row_length,rows)  ! fall vel. for graupel (m/s)
      REAL CTTEMP(row_length,rows)
      REAL RAINFRAC(row_length,rows)
      REAL FRAC_ICE_ABOVE(row_length,rows) ! Cloud ice fraction
!                                            in layer above
      REAL layer_thickness(row_length,rows)
      REAL rho1(row_length,rows)
      REAL rho2(row_length,rows)
      REAL deltaz(row_length,rows,q_levels)
      REAL rhodz_dry(row_length,rows,q_levels)
      REAL rhodz_moist(row_length,rows,q_levels)
      REAL q_total(row_length,rows)


! Allocate CX and CONSTP arrays
! Start C_LSPSIZ
! Description: Include file containing idealised forcing options
! Author:      R. Forbes
!
! History:
! Version  Date      Comment
! -------  ----      -------
!   6.1    01/08/04  Increase dimension for rain/graupel.  R.Forbes
!   6.2    22/08/05  Include the step size between ice categories.
!                                                   Damian Wilson

! Sets up the size of arrays for CX and CONSTP
      REAL CX(100),CONSTP(100)
      INTEGER,PARAMETER:: ice_type_offset=20

! End C_LSPSIZ
! --------------------------COMDECK C_LSPMIC----------------------------
! SPECIFIES MICROPHYSICAL PARAMETERS FOR AUTOCONVERSION, HALLETT MOSSOP
! PROCESS, ICE NUCLEATION. ALSO SPECIFIES NUMBER OF ITERATIONS OF
! THE MICROPHYSICS AND ICE CLOUD FRACTION METHOD
! ----------------------------------------------------------------------
!
! History:
!
! Version    Date     Comment
! -------    ----     -------
!   5.4    16/08/02   Correct comment line, add PC2 parameters and
!                     move THOMO to c_micro     Damian Wilson
!   6.0    11/08/03   Correct value of wind_shear_factor for PC2
!                                                          Damian Wilson
!   6.2    17/11/05   Remove variables that are now in UMUI. D. Wilson
!   6.2    03/02/06   Include droplet settling logical. Damian Wilson
!
! ----------------------------------------------------------------------
!      AUTOCONVERSION TERMS
! ----------------------------------------------------------------------
!
!     LOGICAL, PARAMETER :: L_AUTOCONV_MURK is set in UMUI
! Set to .TRUE. to calculate droplet concentration from MURK aerosol,
! which will override L_USE_SULPHATE_AUTOCONV (second indirect effect
! of sulphate aerosol). If both are .FALSE., droplet concentrations
! from comdeck C_MICRO are used to be consistent with the values
      ! used in the radiation scheme.

      ! This next set of parameters is to allow the 3B scheme to
      ! be replicated at 3C/3D
        ! Inhomogeneity factor for autoconversion rate
        REAL,PARAMETER:: INHOMOG_RATE=1.0

        ! Inhomogeneity factor for autoconversion limit
        REAL,PARAMETER:: INHOMOG_LIM=1.0

        ! Threshold droplet radius for autoconversion
        REAL,PARAMETER:: R_THRESH=7.0E-6
      ! End of 3B repeated code

      !Do not alter R_AUTO and N_AUTO since these values are effectively
      ! hard wired into a numerical approximation in the autoconversion
      ! code. EC_AUTO will be multiplied by CONSTS_AUTO

      ! Threshold radius for autoconversion
      REAL, PARAMETER :: R_AUTO=20.0E-6

      ! Critical droplet number for autoconversion
      REAL, PARAMETER :: N_AUTO=1000.0

      ! Collision coalesence efficiency for autoconversion
!      REAL, PARAMETER :: EC_AUTO is set in UMUI

      ! The autoconversion powers define the variation of the rate with
      ! liquid water content and droplet concentration. The following are
      ! from Tripoli and Cotton

      !  Dependency of autoconversion rate on droplet concentration
      REAL, PARAMETER :: POWER_DROPLET_AUTO=-0.33333

      ! Dependency of autoconversion rate on water content
      REAL, PARAMETER :: POWER_QCL_AUTO=2.33333

      ! Dependency of autoconversion rate on air density
      REAL, PARAMETER :: power_rho_auto=1.33333

      ! CONSTS_AUTO = (4 pi)/( 18 (4 pi/3)^(4/3)) g /  mu (rho_w)^(1/3)
      ! See UM documentation paper 26, equation P26.132

      ! Combination of physical constants
      REAL, PARAMETER :: CONSTS_AUTO=5907.24

      ! Quantites for calculation of drop number by aerosols.
      ! Need only set if L_AUTOCONV_MURK=.TRUE.  See file C_VISBTY

      ! Scaling concentration (m-3) in droplet number concentration
      REAL, PARAMETER :: N0_MURK=500.0E6

      ! Scaling mass (kg/kg) in droplet number calculation from aerosols
      REAL, PARAMETER :: M0_MURK=1.458E-8

      ! Power in droplet number calculation from aerosols
      REAL, PARAMETER :: POWER_MURK=0.5

      ! Ice water content threshold for graupel autoconversion (kg/m^3)
      REAL, PARAMETER :: AUTO_GRAUP_QCF_THRESH = 3.E-4

      ! Temperature threshold for graupel autoconversion (degC)
      REAL, PARAMETER :: AUTO_GRAUP_T_THRESH = -4.0

      ! Temperature threshold for graupel autoconversion
      REAL, PARAMETER :: AUTO_GRAUP_COEFF = 0.5

      !-----------------------------------------------------------------
      ! Iterations of microphysics
      !-----------------------------------------------------------------

      ! Number of iterations in microphysics.
      INTEGER,PARAMETER :: LSITER=1
      ! Advise 1 iteration for every 10 minutes or less of timestep.

      !-----------------------------------------------------------------
      ! Nucleation of ice
      !-----------------------------------------------------------------

      ! Note that the assimilation scheme uses temperature thresholds
      ! in its calculation of qsat.

      ! Nucleation mass
      REAL, PARAMETER :: M0=1.0E-12

      ! Maximum Temp for ice nuclei nucleation (deg C)
      REAL, PARAMETER :: TNUC=-10.0

      ! Maximum temperature for homogenous nucleation is now in c_micro
      ! so that it is available to code outside of section A04.

      !  1.0/Scaling quantity for ice in crystals
      REAL, PARAMETER :: QCF0=1.0E4       ! This is an inverse quantity

      ! Minimum allowed QCF after microphysics
      REAL,PARAMETER:: QCFMIN=1.0E-8

      ! 1/scaling temperature in aggregate fraction calculation
      REAL, PARAMETER :: T_SCALING=0.0384

      !  Minimum temperature limit in calculation  of N0 for ice (deg C)
      REAL, PARAMETER :: T_AGG_MIN=-45.0

      !-----------------------------------------------------------------
      ! Hallett Mossop process
      !-----------------------------------------------------------------

      ! Switch off Hallett Mossop in this version but allow
      ! functionality

      ! Min temp for production of Hallett Mossop splinters (deg C)
      REAL, PARAMETER :: HM_T_MIN=-8.0

      ! Max temp for production of Hallett Mossop splinters (deg C)
      REAL, PARAMETER :: HM_T_MAX=-273.0
      ! REAL, PARAMETER :: HM_T_MAX=-3.0

      !  Residence distance for Hallett Mossop splinters (1/deg C)
      REAL, PARAMETER :: HM_DECAY=1.0/7.0

      ! Reciprocal of scaling liquid water content for HM process
      REAL, PARAMETER :: HM_RQCL=1.0/0.1E-3

      !-----------------------------------------------------------------
      ! PC2 Cloud Scheme Terms
      !-----------------------------------------------------------------

      ! Specifies the ice content (in terms of a fraction of qsat_liq)
      ! that corresponds to a factor of two reduction in the width of
      ! the vapour distribution in the liquid-free part of the gridbox.
      REAL, PARAMETER :: ICE_WIDTH=0.04

      ! Parameter that governs the rate of spread of ice cloud fraction
      ! due to windshear
      REAL, PARAMETER :: WIND_SHEAR_FACTOR = 1.5E-4

!  External subroutines called -----------------------------------------
      EXTERNAL LS_PPNC,LSPCON
!*----------------------------------------------------------------------
!  Physical constants -------------------------------------------------
      REAL CFMIN
      PARAMETER (                                                       &
     & CFMIN=1.0E-3                                                     &
                           ! Used for LS_PPNC  compress.
     &)
!  Define local variables ----------------------------------------------
      INTEGER I,K,j                                                     &
                        ! Loop counters: I - horizontal field index;
!                                      K - vertical level index.
     &,N              ! "nval" for WHEN routine.
!
      REAL work  ! work variable

!-----------------------------------------------------------------------
      ERROR=0

! Define l_non_hydrostatic and l_mixing_ratio. These would be best
! passed from the UMUI, but we will keep them here for the moment.
      if (l_mixing_ratio .or.                                           &
     &   (h_sect(4)  /=  '03B' .and. h_sect(4)  /=  '03C') ) then
        l_non_hydrostatic = .true.   ! This uses most physically accurate
                                     ! code.
      else
        l_non_hydrostatic = .false.  ! This ensures bit reproducibility
      end if


! Define CX and CONSTP values
! DEPENDS ON: lspcon
    CALL LSPCON(CX,CONSTP,x1i,x1ic,x1r,x2r,x4r,m_ci,l_psd,ai,bi,aic,bic,&
     &          lsp_ei,lsp_fi,lsp_eic,lsp_fic                       )
!-----------------------------------------------------------------------
!L Internal structure.
!L 1. Initialise rain and snow to zero.
!   Initialise scavenged amounts of S Cycle tracers to 0 for full field
!-----------------------------------------------------------------------
      do j=1,rows
        DO I=1,row_length
        LSRAIN(I,j)=0.0
        LSSNOW(I,j)=0.0
        LSSNOW2(i,j)=0.0
        LSGRAUP(i,j)=0.0
        droplet_flux(i,j)=0.0
        CTTEMP(I,J)=0.0
        RAINFRAC(I,J)=0.0
        FRAC_ICE_ABOVE(I,J)=0.0
        VFALL(I,j)=0.0
        VFALL2(i,j)=0.0
        VFALL_RAIN(i,j)=0.0
        VFALL_GRAUP(i,j)=0.0
      END DO ! Loop over points,i
      END DO ! Loop over points,j

      LSCAV_SO2(:,:)=0.0
      LSCAV_NH3(:,:)=0.0
      LSCAV_SO4AIT(:,:)=0.0
      LSCAV_SO4ACC(:,:)=0.0
      LSCAV_SO4DIS(:,:)=0.0

      If (l_non_hydrostatic) Then
! ----------------------------------------------------------------------
! Calculate the (non-hydrostatic) layer thicknesses (deltaz) and air
! densities multiplied by deltaz (rhodz_moist and rhodz_dry).
! ----------------------------------------------------------------------
! We should note that this formulation, although better than the
! hydrostatic formulation, is still not entirely conservative. To ensure
! conservation we would need to rewrite the large-scale precipitation
! scheme to consider masses in terms of rho<q>, and
! not the current <rho>q formulation.

        ! We only need to calculate averages for the moist levels
        Do k = 1, q_levels
          Do j = 1, rows
            Do i = 1, row_length

              ! Calculate densities at the boundaries of the layer
              ! by removing the r**2 term from rho_r2.
              ! Rho1 is the density at the lower boundary.
              rho1(i,j)= rho_r2(i,j,k)/( r_rho_levels(i,j,k) *          &
     &                             r_rho_levels(i,j,k) )

              ! Check whether there is a rho level above the current
              ! moist level.
              If (k  <   model_levels) Then
                ! Rho2 is the density at the upper boundary.
                rho2(i,j)= rho_r2(i,j,k+1)/( r_rho_levels(i,j,k+1) *    &
     &                               r_rho_levels(i,j,k+1) )

                ! Calculate the average value of rho across the layer
                ! multiplied by the layer thickness and the layer
                ! thickness.
                rhodz_moist(i,j,k) =                                    &
     &                      rho2(i,j) * ( r_theta_levels(i,j,k) -       &
     &                                         r_rho_levels(i,j,k) )    &
     &                   +  rho1(i,j) * ( r_rho_levels(i,j,k+1) -       &
     &                                    r_theta_levels(i,j,k) )
                deltaz(i,j,k) = r_rho_levels(i,j,k+1)                   &
     &                          - r_rho_levels(i,j,k)

                If (k  ==  1) Then
                  ! For the lowest layer we need to extend the lower
                  ! boundary from the first rho level to the surface.
                  ! The surface is the 0'th theta level.
                  deltaz(i,j,1) = r_rho_levels(i,j,2)                   &
     &                            - r_theta_levels(i,j,0)
                  rhodz_moist(i,j,1) = rhodz_moist(i,j,1)*deltaz(i,j,1) &
     &                   / (r_rho_levels(i,j,2)-r_rho_levels(i,j,1))
                End if  ! k  ==  1

              Else
                ! For a top layer higher than the highest rho level
                ! we can calculate a pseudo rho level. We will assume
                ! it has a similar density to the rho level below
                ! and that the intervening theta level is in the centre
                ! of the layer.
                deltaz(i,j,k) = 2.0*(r_theta_levels(i,j,k)              &
     &                              -r_rho_levels(i,j,k))
                rhodz_moist(i,j,k) = rho1(i,j) * deltaz(i,j,k)

              End if  ! k  <   model_levels

              ! Calculate total moisture
              q_total(i,j) = q(i,j,k) + qcl(i,j,k) + qcf(i,j,k)

              If (l_mcr_qcf2) Then
               q_total(i,j) = q_total(i,j) + qcf2(i,j,k)
              End if  ! l_mcr_qcf2

              If (l_mcr_qrain) Then
               q_total(i,j) = q_total(i,j) + qrain(i,j,k)
              End if  ! l_mcr_qrain

              If (l_mcr_qgraup) Then
               q_total(i,j) = q_total(i,j) + qgraup(i,j,k)
              End if  ! l_mcr_qgraup


              ! Rho_r2 uses the moist density of air. If the mixing
              ! ratio framework is in place then we need to also know
              ! the dry density of air.
              If (l_mixing_ratio) Then
                rhodz_dry(i,j,k) = rhodz_moist(i,j,k)                   &
     &                             / (1.0 + q_total(i,j))
              Else
                rhodz_dry(i,j,k) = rhodz_moist(i,j,k)                   &
     &                             * (1.0 - q_total(i,j))
              End if  ! l_mixing_ratio

            End Do  ! i
          End Do  ! j
        End Do  ! k

      End if  ! L_non_hydrostatic
!
!-----------------------------------------------------------------------
!L 2. Loop round levels from top down (counting bottom level as level 1,
!L    as is standard in the Unified model).
!-----------------------------------------------------------------------
!
      DO K=Q_LEVELS,1,-1

        If (l_it_melting) Then
!-----------------------------------------------------------------------
!      Calculate the number of fall and melting iterations
!-----------------------------------------------------------------------
          n_iterations = min( bl_levels / k + 1 , 5)

        Else
          ! Do not use iterative melting
          n_iterations = 1

        End if  ! l_it_melting

!-----------------------------------------------------------------------
!L 2.5 Form INDEX IX to gather/scatter variables in LS_PPNC
!-----------------------------------------------------------------------
!
!  Set index where cloud fraction > CFMIN or where non-zero pptn
!  Note: whenimd is functionally equivalent to WHENILE (but autotasks).
!
!
        N=0
        Do j = 1,rows
          Do i = 1,row_length
            layer_thickness(i,j) =                                      &
     &       p_layer_boundaries(i,j,k) - p_layer_boundaries(i,j,k-1)

            ! Set up IF statement to determine whether to call the
            ! microphysics code for this grid box (i.e. if there is
            ! already condensate in the grid box or there is
            ! precipitation about to fall into the grid box)
            work = QCF(i,j,k)

            ! Include extra microphysics variables if in use
            If (L_mcr_qcf2 )  work = work + QCF2(i,j,k) + LSSNOW2(i,j)
            If (L_mcr_qrain)  work = work + QRAIN(i,j,k)
            If (L_mcr_qgraup) work = work + QGRAUP(i,j,k)+LSGRAUP(i,j)
            If (L_droplet_settle) work = work + droplet_flux(i,j)

            If (CFL(i,j,k) > CFMIN .OR.                                 &
     &         (LSRAIN(i,j)+LSSNOW(i,j)) > 0.0 .OR. work > 0.0) Then
              ! include this grid box.
              ! Strictly speaking the CFL > CFMIN clause is too
              ! restrictive since ice nucleation does not require
              ! liquid water, but the code would be very messy.
              N = N + 1
              IX(N,1) = i
              IX(N,2) = j
            End If
          End Do ! Loop over points,i
        End Do ! Loop over points,j
!
        IF(N >  0)THEN

          IF (L_SEASALT_CCN) THEN
            sea_salt_ptr=K
            salt_dim_ice=N
          ELSE
            sea_salt_ptr=1
            salt_dim_ice=1
          ENDIF

          IF (L_use_biogenic) THEN
            biogenic_ptr=K
            biog_dim_ice=N
          ELSE
            biogenic_ptr=1
            biog_dim_ice=1
          ENDIF

! DEPENDS ON: ls_ppnc
          CALL LS_PPNC(IX,N,TIMESTEP,n_iterations,                      &
     &                 LSRAIN,LSSNOW,LSSNOW2,LSGRAUP,droplet_flux,      &
     &                 CF(1,1,K),CFL(1,1,K),CFF(1,1,K),                 &
     &                 QCF(1,1,K),QCL(1,1,K),T(1,1,K),                  &
     &                 QCF2(1,1,K),QRAIN(1,1,K),QGRAUP(1,1,K),          &
      ! Process rate diagnostics
     &  PSDEP(1,1,K),PSAUT(1,1,K),PSACW(1,1,K),PSACR(1,1,K)             &
     &, PSACI(1,1,K),PSMLT(1,1,K),PSMLTEVP(1,1,K)                       &
     &, PRAUT(1,1,K),PRACW(1,1,K),PREVP(1,1,K)                          &
     &, PGAUT(1,1,K),PGACW(1,1,K),PGACS(1,1,K),PGMLT(1,1,K)             &
     &, PIFRW(1,1,K),PIPRM(1,1,K),PIDEP(1,1,K),PIACW(1,1,K)             &
     &, PIACR(1,1,K),PIMLT(1,1,K),PIMLTEVP(1,1,K)                       &
     &, PIFALL(1,1,K),PSFALL(1,1,K),PRFALL(1,1,K),PGFALL(1,1,K)         &
     &, PLSET(1,1,K),PLEVPSET(1,1,K)                                    &
      ! Process rate diagnostic switches
     &, L_PSDEP_diag,L_PSAUT_diag,L_PSACW_diag,L_PSACR_diag             &
     &, L_PSACI_diag,L_PSMLT_diag,L_PSMLTEVP_diag                       &
     &, L_PRAUT_diag,L_PRACW_diag,L_PREVP_diag                          &
     &, L_PGAUT_diag,L_PGACW_diag,L_PGACS_diag,L_PGMLT_diag             &
     &, L_PIFRW_diag,L_PIPRM_diag,L_PIDEP_diag,L_PIACW_diag             &
     &, L_PIACR_diag,L_PIMLT_diag,L_PIMLTEVP_diag                       &
     &, L_PIFALL_diag,L_PSFALL_diag,L_PRFALL_diag,L_PGFALL_diag         &
     &, L_PLSET_diag,L_PLEVPSET_diag,                                   &
     &                 L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_psd,    &
     &                 L_it_melting,                                    &
     &                 l_non_hydrostatic, l_mixing_ratio, l_cry_agg_dep,&
     &                 l_droplet_settle,                                &
     &                 SO2(1,1,K),L_SULPC_SO2,                          &
     &                 NH3(1,1,K),L_SULPC_NH3,                          &
     &                 SO4_AIT(1,1,K),SO4_ACC(1,1,K),SO4_DIS(1,1,K),    &
     &                 AGED_SOOT(1,1,K),L_SOOT,                         &
     &                 BMASS_AGD(1,1,K),BMASS_CLD(1,1,K),L_BIOMASS_CCN, &
     &                 OCFF_AGD(1,1,K),OCFF_CLD(1,1,K),L_OCFF_CCN,      &
     &                 AEROSOL(1,1,K),L_MURK, l_pc2,                    &
     &                 LSCAV_SO2, LSCAV_NH3, LSCAV_SO4AIT,              &
     &                 LSCAV_SO4ACC,LSCAV_SO4DIS,                       &
     &                 SNOW_DEPTH, LAND_FRACT,                          &
     &                 L_USE_SULPHATE_AUTOCONV,                         &
     &                 L_AUTO_DEBIAS,                                   &
     &                 SEA_SALT_FILM(1,1,sea_salt_ptr),                 &
     &                 SEA_SALT_JET(1,1,sea_salt_ptr), L_SEASALT_CCN,   &
     &                 salt_dim1, salt_dim2, salt_dim_ice,              &
     &                 L_use_biogenic, biogenic(1, 1, biogenic_ptr),    &
     &                 biog_dim_ice,                                    &
     &                 Q(1,1,K),                                        &
     &                 p_theta_levels(1,1,K),layer_thickness,           &
     &                 deltaz(1,1,k), rhodz_dry(1,1,k),                 &
     &                 rhodz_moist(1,1,k),                              &
     &                 row_length,rows,rhc_row_length,rhc_rows,         &
     &                 BLAND,CW_SEA,                                    &
     &                 CW_LAND,                                         &
     &                 RHCRIT(1,1,K),                                   &
     &                 VFALL,VFALL2,VFALL_RAIN,VFALL_GRAUP,             &
     &                 FRAC_ICE_ABOVE,                                  &
     &                 CTTEMP,RAINFRAC,CX,CONSTP,                       &
     &                 L_seq_mcr,L_autoc_3b,L_autolim_3b,               &
     &                 L_autoconv_murk,ec_auto,                         &
     &                 N_drop_land,N_drop_sea,N_drop_land_cr,           &
     &                 N_drop_sea_cr,Ntot_land, Ntot_sea,               &
     &                 ai, bi, aic, bic                                 &
     &                 )
        ENDIF
!
! Copy rainfall and snowfall rates to 3D fields for diagnostic output
!
        IF (lspice_dim1  ==  row_length .AND. lspice_dim2  ==  rows     &
     &      .AND. lspice_dim3  ==  Q_LEVELS) THEN
! Only copy rain and snow to 3D fields if arrays are dimensionalized.
          do j=1,rows
            do i=1,row_length
              LSRAIN3D(I,J,K)=LSRAIN(I,J) + droplet_flux(i,j)
              LSSNOW3D(I,J,K)=LSSNOW(I,J)+LSSNOW2(I,J)+LSGRAUP(I,J)
              RAINFRAC3D(I,J,K)=RAINFRAC(I,J)
            end do
          end do
        ENDIF
!
      END DO ! Loop over K

! Add together ice crystals, snow aggregates and graupel
! for surface snow rate (kg/m2/s)

      If (L_mcr_qcf2) Then
        Do j = 1,rows
          Do i = 1,row_length
            LSSNOW(i,j) = LSSNOW(i,j) + LSSNOW2(i,j) + LSGRAUP(i,j)
          End Do
        End Do
      End If

      If (L_droplet_settle) Then
! Add droplet settling to the large-scale precip rate in order
! to be able to have a full water budget with the precipitation
! diagnostics.
        Do j = 1,rows
          Do i = 1,row_length
            LSRAIN(i,j) = LSRAIN(i,j) + droplet_flux(i,j)
          End Do
        End Do
      End If

 20   CONTINUE                  ! Branch for error exit
      RETURN
      END SUBROUTINE LS_PPN
!*LL  SUBROUTINE LS_PPNC------------------------------------------------
!*L  Arguments:---------------------------------------------------------
