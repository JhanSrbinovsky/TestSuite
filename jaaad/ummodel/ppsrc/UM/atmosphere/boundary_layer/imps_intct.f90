
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!! Subroutine IMPS_INTCT -------------------------------------------
!!!
!!! Purpose : Intermediate control level to call requested version of
!!!           IMP_SOLVER with the appropriate arguments.
!!!
!!!  Model            Modification history:
!!! version  Date
!!!  5.2   15/11/00   New Deck         M. Best
!    5.3  25/04/01  Add coastal tiling. Nic Gedney
!!!  5.4   15/05/02   Pass through heights arrays for implicit
!!!                   solver on a sphere   Adrian Lock
!!!  5.5   17/02/03   Pass through ice catagory arrays. J. Ridley
!!!  5.5   19/02/03   Remove redundent RML and L_BL_LSPICE code.
!!!                                        Adrian Lock
!!!  6.1  01/09/04  Pass potential evaporation related variables.
!                                                          Nic Gedney
!!!  6.1  17/05/04  Changes to the BL solver to enable phys2
!!!                 substepping.                       M. Diamantakis
! 6.2      21/02/06    Switch for mixing ratios   A.P.Lock
!!!  6.2  11/01/06  Remove MOSES I code                A.P.Lock
!!!  6.2  02/02/06  Passes L_flux_bc through argument list to allow
!!!                 settings for prescribed surface flux forcing
!!!                                                           R. Wong
!    6.4   10/01/07  Introduce unconditionally stable and non-oscillatory
!                    BL numerical solver               M. Diamantakis 
!!!
!!! Programming standard : unified model documentation paper No 3
!!!
!!! System components covered : P24
!!!
!!! System task : P0
!!!
!!!END -----------------------------------------------------------------

!    Arguments :-
      SUBROUTINE IMPS_INTCT (                                           &

! IN mpp variables
     & halo_i, halo_j, off_x, off_y, row_length, rows, n_rows,          &
     & global_row_length, proc_row_group, at_extremity,                 &
     & n_proc, n_procx, n_procy, neighbour,                             &

! IN values defining field dimensions and subset to be processed :
     & NTILES, land_pts, NICE,                                          &

! IN values defining vertical grid of model atmosphere :
     & MODEL_DOMAIN,                                                    &
     & BL_LEVELS,                                                       &
     & r_rho_levels, r_theta_levels,                                    &
     & GAMMA,                                                           &

! IN Substepping info 
     & Substep_Number, Num_Substeps, L_phys2_substep,                   &

! IN New BL solver parameters
     & L_us_blsol, Puns, Pstb,                                          &

! IN U and V momentum fields.
     & U, V, U_conv, V_conv,                                            &
! IN Non turbulent increments to momentum
!  (New dynamics only).
     & R_u,R_v,                                                         &

! IN soil/vegetation/land surface data :
     & LAND_MASK,LAND_INDEX,                                            &
     & ST_LEVELS,SM_LEVELS,                                             &

! IN sea/sea-ice data :
     & DI,ICE_FRACT,DI_NCAT,ICE_FRACT_NCAT,U_0,V_0,                     &

! IN cloud data :
     & Q,QCF,QCL,Q_conv,QCF_conv,QCL_conv,QCF_latest,QCL_latest,        &
     & T,T_conv,                                                        &

! IN everything not covered so far :
     & PSTAR,SURF_RADFLUX,TIMESTEP,L_SICE_HEATFLUX,                     &

! IN variables from BDY_LAYR (that used to be local arrays)
     & ALPHA1_SICE,ASHTF,DTRDZ_CHARNEY_GRID,RDZ_CHARNEY_GRID,           &
     & DTRDZ_U,DTRDZ_V,RDZ_U,RDZ_V,                                     &
     & CDR10M_U,CDR10M_V,Z1_TQ,                                         &
!ajm extra variable added
     & RHOKM_u,RHOKM_v,                                                 &

! IN variables for new BL solver
     & bl_type_1,bl_type_2,bl_type_3,bl_type_4,                         &
     & bl_type_5,bl_type_6,bl_type_7,                                   &

! IN additional variables for MOSES II
     & TILE_PTS,TILE_INDEX,TILE_FRAC,CANOPY,                            &
     & L_NEG_TSTAR,ALPHA1,FRACA,RHOKH_TILE,SMC,CHR1P5M,                 &
     & RESFS,Z0HSSI,Z0MSSI,CANHC_TILE,FLAKE,                            &
     & WT_EXT_TILE,LW_DOWN,SW_TILE,ASHTF_TILE,                          &
     & FQW_ICE,FTL_ICE,RESFT,RHOKH_SICE,RHOKPM,RHOKPM_POT,RHOKPM_SICE,  &
     & Z0H_TILE,Z0M_TILE,CHR1P5M_SICE,                                  &
     & FLAND,FLANDG,FLANDG_U,FLANDG_V,TSTAR_SEA,ANTHROP_HEAT,           &

!    EAK
!    IN
     &  l_cable                                                         &
     &, surf_down_sw,alb_tile,cos_zenith_angle,l_tile_pts               &
!     &, surf_down_sw,alb_tile,cos_zenith_angle               &
     &, lat,long,day,time_sec                                           &
     &, ls_rain, ls_snow, conv_rain, conv_snow                          &
     &, SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L                   &
     &, SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE             &
     &, TSOIL_TILE,T_SURF_TILE,HCONS,clapp                              &
     &, SATHH,SATCON,HCAP,HCON,Z1_UV                                    &
     &, smvccl, smvcwt, smvcst, sthf, sthu                              &
     &, slm                                                             &
     &, SOIL_TYPE,VEG_TYPE                                              &
     &, ISNOW_FLG3L,total_nsteps                                        &
     &, FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB                         &
     &, TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB                     &
     &, USTAR_CAB,SURF_HTF_CAB                                          &

!sxy
     &, TOT_ALB                                                         &
     &, SNAGE_TILE,RTSOIL_TILE                                          &
     &, GFLUX_TILE,SGFLUX_TILE                                          &
     &, GC,GS,CANOPY_GB                                                 & ! Added GC, pfv 25oct13a
! Lestevens 18 March 2010
     &, DIM_CS1, DIM_CS2    &
     &, NPP, NPP_FT         &
     &, GPP, GPP_FT         &
     &, RESP_S, RESP_S_TOT  &
     &, RESP_S_TILE         & !kdcorbin, 10/10
     &, RESP_P, RESP_P_FT,  &
     &  G_LEAF,             & !kdcorbin, 10/10
     &  TRANSP_TILE,        & !Lestevens 3Nov11
! Lestevens Sept2012: CasaCNP variables 
     &  CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,                               &
     &  SOIL_ORDER,GLAI,PHENPHASE,                                      &
     ! Lestevens 23apr13
     &  NPP_FT_ACC,RESP_W_FT_ACC,                                       &

! INOUT data :
     & T_SOIL, TI, TI_GB, TSTAR, T_latest,Q_latest,                     &

! INOUT Diagnostics started in BDY_LAYR not requiring STASH flags :
     & E_SEA,FQW,FTL,H_SEA,RHOKH,TAUX,TAUY,                             &

! INOUT additional variables for MOSES II
     & TSTAR_TILE,FQW_TILE,EPOT_TILE,FTL_TILE,                          &
     & SNOW_TILE,LE_TILE,RADNET_SICE,RADNET_TILE,OLR,                   &
     & TSTAR_SICE,TSTAR_SSI,                                            &
     & TAUX_LAND,TAUX_SSI,TAUY_LAND,TAUY_SSI,                           &

! OUT Increments to U and V momentum fields and Tl qw
!  (New dynamics only).
     & DU,DV,                                                           &

! OUT Diagnostic not requiring STASH flags :
     & RHOKH_mix,SEA_ICE_HTF,SURF_HT_FLUX,                              &
     & SURF_HT_FLUX_LAND,SURF_HT_FLUX_SICE,                             &

! OUT diagnostics requiring STASH flags :
     & SICE_MLT_HTF,SNOMLT_SURF_HTF,LATENT_HEAT,                        &
     & Q1P5M,T1P5M,U10M,V10M,                                           &

! (IN) STASH flags :-
     & SIMLT,SMLT,SLH,SQ1P5,ST1P5,SU10,SV10,l_ftl,l_fqw,l_taux,l_tauy,  &

! OUT additional variables for MOSES II
     & ESOIL_TILE,ES,EI_TILE,                                           &
     & Q1P5M_TILE,T1P5M_TILE,ECAN_TILE,MELT_TILE,                       &
     & TSTAR_LAND,E_SSI,EI_SICE,FTL_SSI,                                &
     & ERROR,                                                           &

! OUT data required elsewhere in UM system :
     & ECAN,EI,EXT,SNOWMELT,                                            &
     & lq_mix_bl,                                                       &

! SCM namelist logical
     & L_FLUX_BC,                                                       &
     & LTIMER                                                           &
     & )

      IMPLICIT NONE

!  Inputs :-

! (a) Defining horizontal grid and subset thereof to be processed.

      INTEGER                                                           &
! cjj additions - mpp variables.
     &  row_length                                                      &
                   ! Local number of points on a row
     &, rows                                                            &
                   ! Local number of rows in a theta field
     &, n_rows                                                          &
                   ! Local number of rows in a v field
     &, halo_i                                                          &
                   ! Size of halo in i direction.
     &, halo_j                                                          &
                   ! Size of halo in j direction.
     &, off_x                                                           &
                   ! Size of small halo in i
     &, off_y                                                           &
                   ! Size of small halo in j.
     &, global_row_length                                               &
                           ! number of points on a row
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &, n_proc                                                          &
                   ! Total number of processors
     &, n_procx                                                         &
                   ! Number of processors in longitude
     &, n_procy                                                         &
                   ! Number of processors in latitude
     &, neighbour(4)   ! Array with the Ids of the four neighbours in
                       !   the horizontal plane

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         !   south, east or west of the processor grid

! Parameters
      Integer                                                           &
     &   PNorth,                                                        &
                      ! North processor address in the neighbor array
     &   PEast,                                                         &
                      ! East processor address in the neighbor array
     &   PSouth,                                                        &
                      ! South processor address in the neighbor array
     &   PWest,                                                         &
                      ! West processor address in the neighbor array
     &   NoDomain     ! Value in neighbor array if the domain has
                      !  no neighbor in this direction. Otherwise
                      !  the value will be the tid of the neighbor
      Parameter (                                                       &
     &   PNorth   = 1,                                                  &
     &   PEast    = 2,                                                  &
     &   PSouth   = 3,                                                  &
     &   PWest    = 4,                                                  &
     &   NoDomain = -1)

      Integer                                                           &
     & NTILES                                                           &
                                 ! IN number of land tiles
     &,land_pts                                                         &
                                 ! IN No.of land points in whole grid.
     &,NICE                      ! IN No. of sea ice catagories

! (b) Defining vertical grid of model atmosphere.

      INTEGER                                                           &
     & BL_LEVELS                                                        &
                                   ! IN Max. no. of "boundary" levels
!                                     allowed.Assumed <= 30 for dim-
!                                     sioning of GAMMA in common deck
!                                     C_GAMMA used in SF_EXCH and KMKH
     &,MODEL_DOMAIN                                                     &

     &,Substep_Number                                                   &
     &,Num_Substeps

! Switch for calculating exchange coeffs from latest values.
      Logical                                                           &
     & L_phys2_substep                                                  &
! SCM logical for prescribing surface flux forcing
     &, L_flux_bc

      Logical :: L_us_blsol   ! switch for BL numerical solver

      Real :: Puns, Pstb ! parameters for uncond stable numerical solver
                         ! Puns : used in an unstable BL column
                         ! Pstb : used in an stable BL column
      Real                                                              &
     &  bl_type_1(row_length,rows)                                      &
                                  ! IN Indicator set to 1.0 if stable
!                                 !     b.l. diagnosed, 0.0 otherwise.
     &, bl_type_2(row_length,rows)                                      &
                                  ! IN Indicator set to 1.0 if Sc over
!                                 !     stable surface layer diagnosed,
!                                 !     0.0 otherwise.
     &, bl_type_3(row_length,rows)                                      &
                                  ! IN Indicator set to 1.0 if well
!                                 !     mixed b.l. diagnosed,
!                                 !     0.0 otherwise.
     &, bl_type_4(row_length,rows)                                      &
                                  ! IN Indicator set to 1.0 if
!                                 !     decoupled Sc layer (not over
!                                 !     cumulus) diagnosed,
!                                 !     0.0 otherwise.
     &, bl_type_5(row_length,rows)                                      &
                                  ! IN Indicator set to 1.0 if
!                                 !     decoupled Sc layer over cumulus
!                                 !     diagnosed, 0.0 otherwise.
     &, bl_type_6(row_length,rows)                                      &
                                  ! IN Indicator set to 1.0 if a
!                                 !     cumulus capped b.l. diagnosed,
!                                 !     0.0 otherwise.
     &, bl_type_7(row_length,rows)                                      
                                  ! IN Indicator set to 1.0 if a
!                                 !     shear-dominated b.l. diagnosed,
!                                 !     0.0 otherwise.

      LOGICAL                                                           &
     & lq_mix_bl              ! TRUE if mixing ratios used in
!                             ! boundary layer code

! Co-ordinate arrays
      Real                                                              &
           ! local vertical co-ordinate information
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                   1-halo_j:rows+halo_j,0:bl_levels)              &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &                 1-halo_j:rows+halo_j, bl_levels)

      REAL                                                              &
     &  GAMMA(bl_levels)                                                &
     &, R_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &        BL_LEVELS)                                                &
                            ! non-turbulent increment to u wind field
     &, R_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,             &
     &        BL_LEVELS)                                                &
                            ! non-turbulen increment to v wind field
     &, U(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      BL_LEVELS)                                                  &
     &, V(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      BL_LEVELS)

! (c) Soil/vegetation/land surface parameters (mostly constant).

      LOGICAL                                                           &
     & LAND_MASK(row_length,rows)  ! IN T if land, F elsewhere.

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




      INTEGER                                                           &
     & LAND_INDEX(land_pts)      ! IN LAND_INDEX(I)=J => the Jth
!                                     point in P_FIELD is the Ith
!                                     land point.

      INTEGER                                                           &
     & ST_LEVELS                                                        &
                                   ! IN No. of deep soil temp. levels
     &,SM_LEVELS                   ! IN No. of soil moisture levels


! (d) Sea/sea-ice data.

      REAL                                                              &
     & DI(row_length,rows)                                              &
                                   ! IN "Equivalent thickness" of
!                                     sea-ice(m).
     &,ICE_FRACT(row_length,rows)                                       &
                                   ! IN Fraction of gridbox covered by
!                                     sea-ice (decimal fraction).
     &,DI_NCAT(row_length,rows,nice)                                    &
                                            !IN "Equivalent thickness"
!                                     of ice catagory in grid box (m).
     &,ICE_FRACT_NCAT(row_length,rows,nice)                             &
                                            !IN Fraction of ice catagory
!                                     in gridbox (decimal fraction).
     &,U_0(row_length,rows)                                             & 
                                   ! IN W'ly component of surface
!                                     current (m/s).
     &,V_0(row_length,n_rows)      ! IN S'ly component of surface
!                                     current (m/s).

! (e) Cloud data.

      REAL                                                              &
     & QCF(1-halo_i:row_length+halo_i,                                  &
     &     1-halo_j:rows+halo_j,BL_LEVELS)                              &
                                          ! IN Cloud ice (kg per kg air)
     &,QCL(1-halo_i:row_length+halo_i,                                  &
     &     1-halo_j:rows+halo_j,BL_LEVELS)                              &
                                          ! IN Cloud liquid water
     &,Q(1-halo_i:row_length+halo_i,                                    &
     &   1-halo_j:rows+halo_j,BL_LEVELS)                                &
                                          ! IN specific humidity
     &,T(row_length,rows,BL_LEVELS)                                     &
                                          ! IN temperature
                          ! Latest estimates to time level n+1 values
     &,QCF_latest(row_length,rows,BL_LEVELS)                            &
                                ! IN Cloud ice (kg per kg air)
     &,QCL_latest(row_length,rows,BL_LEVELS)                            &
                                             ! IN Cloud liquid water
     &,T_conv(row_length, rows, bl_levels)                              &
     &,q_conv(row_length, rows, bl_levels)                              &
     &,qcl_conv(row_length, rows, bl_levels)                            &
     &,qcf_conv(row_length, rows, bl_levels)                            &
     &,u_conv(1-off_x:row_length+off_x,1-off_y:rows+off_y,              &
     &              bl_levels)                                          &
     &,v_conv(1-off_x:row_length+off_x,1-off_y:n_rows+off_y,            &
     &              bl_levels)

! (f) Atmospheric + any other data not covered so far, incl control.

      REAL                                                              &
     & PSTAR(row_length,rows)                                           &
                                   ! IN Surface pressure (Pascals).
     &,SURF_RADFLUX(row_length,rows)                                    &
                                    ! IN Surface net radiation
!                                    (W/sq m, positive downwards).
     &,TIMESTEP                    ! IN Timestep (seconds).

      LOGICAL                                                           &
     & L_SICE_HEATFLUX             ! IN T: semi-implicit sea-ice temp

! IN variables from BDY_LAYR (that used to be local arrays)
      REAL                                                              &
     & ALPHA1_SICE(ROW_LENGTH,ROWS)                                     &
!                                 ! IN ALPHA1 for sea-ice.
!                                 !    gridbox ALPHA1 for MOSES I
     &,ASHTF(row_length,rows)                                           &
                                  ! IN Coefficient to calculate surface
!                                 heat flux into soil or sea-ice.
     &,DTRDZ_CHARNEY_GRID(row_length,rows,BL_LEVELS)                    &
                                ! IN -g.dt/dp for model layers.
     &,RDZ_CHARNEY_GRID(row_length,rows,BL_LEVELS)                      &
                                ! IN RDZ(,1) is the reciprocal of the
!                                 height of level 1, i.e. of the
!                                 middle of layer 1.  For K > 1,
!                                 RDZ(,K) is the reciprocal
!                                 of the vertical distance
!                                 from level K-1 to level K.
     &,DTRDZ_U(row_length,rows,BL_LEVELS)                               &
                                          ! IN
     &,DTRDZ_V(row_length,n_rows,BL_LEVELS)                             &
                                           ! IN
!                                 -g.dt/dp for model wind layers.
     &,RDZ_U(row_length,rows,2:BL_LEVELS)                               &
                                          !  IN RDZ (K > 1) on UV-grid.
     &,RDZ_V(row_length,n_rows,2:BL_LEVELS)                             &
                                !  IN RDZ (K > 1) on UV-grid.
     &,CDR10M_U(row_length,rows)                                        &
                                ! IN Ratio of CD's reqd for calculation
     &,CDR10M_V(row_length,n_rows)                                      &
                                ! IN Ratio of CD's reqd for calculation
     &,Z1_TQ(row_length,rows)   ! IN Height of lowest theta level.
                                ! ajm  extra variable added

      LOGICAL LTIMER               ! Logical switch for TIMER diags

! IN additional variables for MOSES II

      LOGICAL                                                           &
     & L_NEG_TSTAR                 ! IN Switch for -ve TSTAR error check

      INTEGER                                                           &
     & TILE_PTS(NTYPE)                                                  &
                                 ! IN Number of tile points.
     &,TILE_INDEX(LAND_PTS,NTYPE)&
!                                ! IN Index of tile points.
! Lestevens March 2010
     &, DIM_CS1 &
     &, DIM_CS2

      REAL                                                              &
     & TILE_FRAC(land_pts,ntiles)                                       &
                                ! IN fractional coverage for each
                                !    surface tile
     &,CANOPY(land_pts,ntiles)                                          &
                                ! IN Surface/canopy water (kg/m2)
     &,ALPHA1(land_pts,ntiles)                                          &
                               ! IN Mean gradient of saturated
!                                 specific humidity with
!                                 respect to temperature between
!                                 the bottom model layer and the
!                                 tile surfaces.
     &,FRACA(land_pts,ntiles)                                           &
                                ! IN Fraction of surface
                                !            moisture flux with only
                                !            aerodynamic resistance.
     &,RHOKH_TILE(land_pts,ntiles)                                      &
                                   ! IN
!                                 Tile surface exchange coefficients
!                                 for heat
     &,SMC(land_pts,NTILES)                                             &
                               ! IN Soil moisture content in root depth
!                                  (kg/m2).
     &,CHR1P5M(land_pts,ntiles)                                         &
                                ! IN Ratio of coefficients reqd for
!                                 calculation of 1.5 m T.
     &,RESFS(land_pts,ntiles)                                           &
                              ! IN Combined soil, stomatal
!                                 and aerodynamicresistance
!                                 factor = PSIS/(1+RS/RA) for
!                                 fraction (1-FRACA)
     &,Z0HSSI(row_length,rows)                                          &
     &,Z0MSSI(row_length,rows)                                          &
                                ! IN Roughness lengths over sea
     &,CANHC_TILE(LAND_PTS,NTILES)                                      &
!                               ! IN Areal heat capacity of canopy
!                               !    for land tiles (J/K/m2).
     &,FLAKE(LAND_PTS,NTILES)                                           &
                                ! IN Lake fraction.
     &,WT_EXT_TILE(LAND_PTS,SM_LEVELS,NTILES)                           &
!                               ! IN Fraction of evapotranspiration
!                               !    which is extracted from each
!                               !    soil layer by each tile.
     &,LW_DOWN(ROW_LENGTH,ROWS)                                         &
                                ! IN Surface downward LW radiation
!                               !    (W/m2).
     &,SW_TILE(LAND_PTS,NTILES)                                         &
                                ! IN Surface net SW radiation on land
!                               !    tiles (W/m2).
     &,ASHTF_TILE(LAND_PTS,NTILES)                                      &
!                               ! IN Coefficient to calculate
!                               !    surface heat flux into land
!                               !    tiles.
     &,FQW_ICE(ROW_LENGTH,ROWS)                                         &
                                ! IN Surface FQW for sea-ice
     &,FTL_ICE(ROW_LENGTH,ROWS)                                         &
                                ! IN Surface FTL for sea-ice
     &,RESFT(LAND_PTS,NTILES)                                           &
                                ! IN Total resistance factor.
!                               !    FRACA+(1-FRACA)*RESFS for
!                               !    snow-free land, 1 for snow.
     &,RHOKH_SICE(ROW_LENGTH,ROWS)                                      &
!                               ! IN Surface exchange coefficients
!                               !    for sea and sea-ice
     &,RHOKPM(LAND_PTS,NTILES)                                          &
                                ! IN Land surface exchange coeff.
     &,RHOKPM_POT(LAND_PTS,NTILES)                                      &
!                               ! IN Land surface exchange coeff.
!                                    for potential evaporation.
     &,RHOKPM_SICE(ROW_LENGTH,ROWS)                                     &
!                               ! IN Sea-ice surface exchange coeff.
     &,Z0H_TILE(LAND_PTS,NTILES)                                        &
                                ! IN Tile roughness lengths for heat
!                               !    and moisture (m).
     &,Z0M_TILE(LAND_PTS,NTILES)                                        &
                                ! IN Tile roughness lengths for
!                               !    momentum.
     &,CHR1P5M_SICE(ROW_LENGTH,ROWS)                                    &
!                               ! IN CHR1P5M for sea and sea-ice
!                               !    (leads ignored).
     &,FLAND(LAND_PTS)                                                  &
                                ! IN Land fraction on land points.
     &,FLANDG(ROW_LENGTH,ROWS)                                          &
                                ! IN Land fraction on all points.
     &,FLANDG_U(ROW_LENGTH,ROWS)                                        &
!                               ! IN Land frac (on U-grid, with 1st
!                               !    and last rows undefined or, at
!                               !    present, set to "missing data")
     &,FLANDG_V(ROW_LENGTH,N_ROWS)                                      &
!                               ! IN Land frac (on V-grid, with 1st
!                               !    and last rows undefined or, at
!                               !    present, set to "missing data")
     &,TSTAR_SEA(ROW_LENGTH,ROWS)                                       &
!                               ! IN Open sea sfc temperature (K).
     &,ANTHROP_HEAT(NTILES)
!                               ! IN Additional heat source on tiles
!                               !    for anthropogenic urban heat
!                               !    source (W/m2)

!  STASH flags :-

      LOGICAL                                                           &
     & SIMLT                                                            &
               ! IN Flag for SICE_MLT_HTF (q.v.)
     &,SMLT                                                             &
               ! IN Flag for SNOMLT_SURF_HTF (q.v.)
     &,SLH                                                              &
               ! IN Flag for LATENT_HEAT (q.v.)
     &,SQ1P5                                                            &
               ! IN Flag for Q1P5M (q.v.)
     &,ST1P5                                                            &
               ! IN Flag for T1P5M (q.v.)
     &,SU10                                                             &
               ! IN Flag for U10M (q.v.)
     &,SV10                                                             &
               ! IN Flag for V10M (q.v.)

     &,l_ftl                                                            &
     &,l_fqw                                                            &
     &,l_taux                                                           &
     &,l_tauy
!  In/outs :-

      REAL                                                              &
     & Q_latest(row_length,rows,BL_LEVELS)                              &
                                            ! IN specific humidity
     &,T_latest(row_length,rows,BL_LEVELS)                              &
                                            ! IN temperature
     &,T_SOIL(land_pts,SM_LEVELS)                                       &
                                  ! INOUT Soil temperatures (K).
     &,TI(row_length,rows,nice)                                         &
                                ! INOUT Sea-ice surface layer
!                                      temperature on catagories (K).
     &,TI_GB(row_length,rows)                                           &
                                ! INOUT Sea-ice surface temp (gridbox).
     &,TSTAR(row_length,rows)   ! INOUT Surface temperature (K).


! INOUT additional variables for MOSES II
      REAL                                                              &
     & TSTAR_TILE(land_pts,ntiles)                                      &
                                ! INOUT Surface tile temperature
     &,FQW_TILE(land_pts,ntiles)                                        &
                                ! INOUT surface tile moisture flux
     &,EPOT_TILE(land_pts,ntiles)                                       &
!                               ! INOUT surface tile potential
!                               !       evaporation
     &,FTL_TILE(land_pts,ntiles)                                        &
!                               ! INOUT surface tile heat flux
     &,SNOW_TILE(LAND_PTS,NTILES)                                       &
!                               ! INOUT Snow on tiles (kg/m2).
     &,LE_TILE(LAND_PTS,NTILES)                                         &
                                ! INOUT Surface latent heat flux for
!                               !       land tiles (W/m2).
     &,RADNET_SICE(ROW_LENGTH,ROWS)                                     &
!                               ! INOUT Sea-ice surface net radiation.
     &,RADNET_TILE(LAND_PTS,NTILES)                                     &
!                               ! INOUT Tile surface net radiation.
     &,OLR(ROW_LENGTH,ROWS)                                             &
                                ! IN    TOA - surface upward LW on
!                               !       last radiation timestep
!                               ! OUT   Corrected TOA outward LW
     &,TSTAR_SICE(ROW_LENGTH,ROWS)                                      &
!                               ! INOUT Sea-ice sfc temperature (K).
     &,TSTAR_SSI(ROW_LENGTH,ROWS)
!                               ! INOUT Sea mean sfc temperature (K).

!     EAK
      REAL                                                              &
     &  alb_tile(land_pts,ntiles,4)                                  &
     &, surf_down_sw(row_length,rows,4)                                 &
     &, cos_zenith_angle(row_length,rows)                               &
     &, lat(row_length, rows)                                           &
                                          ! Lat. of gridpoint chosen
     &, long(row_length, rows)                                          &
                                          ! Long. of gridpoint chosen
     &, ls_rain(row_length, rows)                                       &
     &, ls_snow(row_length, rows)                                       &
     &, conv_rain(row_length, rows)                                     &
     &, conv_snow(row_length, rows)                                     &
     &, smvccl (land_pts)                                            &
                             ! soil/qrparm.soil.crit
     &, smvcwt (land_pts)                                            &
                             ! soil/qrparm.soil.wilt
     &, smvcst (land_pts)                                            &
                             ! soil/qrparm.soil.satn
     &, sthf(land_pts,sm_levels)                                    &
                                ! IN Frozen soil moisture content of
                                !     each layer as a fraction of
                                !     saturation.
     &, sthu(land_pts,sm_levels)                                    &
                                ! IN Unfrozen soil moisture content
                                !    of each layer as a fraction of
                                !    saturation.
     &, time_sec                                                        &
                                ! actual time of day in secs.
     &, SNOW_DEPTH3L(land_pts,NTILES,3)                          &
     &, SNOW_MASS3L(land_pts,NTILES,3)                           &
     &, SNOW_COND(land_pts,NTILES,3)                             &
     &, SNOW_TMP3L(land_pts,NTILES,3)                            &
     &, SNOW_RHO3L(land_pts,NTILES,3)                            &
     &, SNOW_RHO1L(land_pts,NTILES)                              &
     &, SNAGE_TILE(land_pts,NTILES)                              &
     &, SMCL_TILE(land_pts,NTILES,sm_levels)                     &
     &, STHU_TILE(land_pts,NTILES,sm_levels)                     &
     &, STHF_TILE(land_pts,NTILES,sm_levels)                     &
     &, TSOIL_TILE(land_pts,NTILES,sm_levels)                    &
     &, T_SURF_TILE(land_pts,NTILES)                             &
     &, RTSOIL_TILE(land_pts,NTILES)                             &
     &, GFLUX_TILE(land_pts,NTILES)                              &
     &, SGFLUX_TILE(land_pts,NTILES)                             &
     &, Z1_UV(row_length,rows)                                   &
                                ! Height of lowest u,v level.
     &, HCONS(land_pts)                                          &
! Lestevens Sept2012: CasaCNP variables 
     &, CPOOL_TILE(land_pts,NTILES,10)                           &
     &, NPOOL_TILE(land_pts,NTILES,10)                           &
     &, PPOOL_TILE(land_pts,NTILES,12)                           &
     &, SOIL_ORDER(land_pts)                                     &
     &, GLAI(land_pts,NTILES)                                    &
     &, PHENPHASE(land_pts,NTILES)                               &
     ! Lestevens 23apr13
     &, NPP_FT_ACC(land_pts,NTILES)                              &
     &, RESP_W_FT_ACC(land_pts,NTILES)

       Real                                                              &
     &  clapp(land_pts)                                              &
                             !  qrparm.soil.bwag ?
!                               Clapp-Hornberger exponent.
     &, satcon(land_pts)                                             &
                                 !  qrparm.soil.satcon
     &, sathh(land_pts)                                              &
                             !  soil water suction
     &, hcon (land_pts)                                              &
                             ! soil/qrparm.soil.hcond
     &, hcap (land_pts)                                              &
                             ! soil/qrparm.soil.hcap
     &, slm(land_pts, sm_levels)
                             !  qrclim.smc_pm(lev).(month)

!  EAK    diagnostic variables for CABLE output
      Real                                                           &
     & FTL_TILE_CAB(LAND_PTS,NTILES)                                 &
     &,FTL_CAB(LAND_PTS)                                             &
     &,LE_TILE_CAB(LAND_PTS,NTILES)                                  &
     &,LE_CAB(LAND_PTS)                                              &
     &,TSTAR_TILE_CAB(LAND_PTS,NTILES)                               &
     &,TSTAR_CAB(LAND_PTS)                                           &
     &,SMCL_CAB(LAND_PTS,SM_LEVELS)                                  &
     &,TSOIL_CAB(LAND_PTS,SM_LEVELS)                                 &
     &,USTAR_CAB(LAND_PTS)                                           &
     &,SURF_HTF_CAB(LAND_PTS)                                        &
     &,TOT_ALB(LAND_PTS,NTILES)                                      & 
     &,CANOPY_GB(LAND_PTS)                                           &
     &,GC(LAND_PTS,NTILES)                                           & ! Added GC, pfv 25oct13a
     &,GS(LAND_PTS)                                                  &
! Lestevens March 2010 - passing CO2 fluxes
     !kdcorbin, 11/10 - changed from NPFT
     &,NPP_FT(LAND_PTS,NTILES)                                       &
     &,NPP(LAND_PTS)                                                 &
     !kdcorbin, 11/10 - changed from NPFT
     &,GPP_FT(LAND_PTS,NTILES)                                       &
     &,GPP(LAND_PTS)                                                 &
     &,RESP_S(LAND_PTS,DIM_CS1)                                      &
     &,RESP_S_TOT(DIM_CS2)                                           &
     &,RESP_S_TILE(LAND_PTS,NTILES)                  & !kdcorbin, 10/10
     &,RESP_P(LAND_PTS)                                              &
     !kdcorbin, 11/10 - changed from NPFT
     &,RESP_P_FT(LAND_PTS,NTILES)                                      &
     &,G_LEAF(LAND_PTS,NTILES)      & !kdcorbin, 10/10
     &,TRANSP_TILE(LAND_PTS,NTILES)   !Lestevens 3Nov11

      Integer                                                           &
     &  day                                                             &
     &, total_nsteps                                                    &
                                ! Total number of steps in run
     &, SOIL_TYPE(row_length,rows)                                      &
     &, VEG_TYPE(row_length,rows)                                       &
     &, ISNOW_FLG3L(land_pts,NTILES)

      Logical                                                           &
     &  l_tile_pts(land_pts,ntiles)

      LOGICAL                                                           &
     & l_cable

!  Outputs :-
!-1 Diagnostic (or effectively so - includes coupled model requisites):-

!  (a) Calculated anyway (use STASH space from higher level) :-
!
      REAL                                                              &
     & E_SEA(row_length,rows)                                           &
                                ! OUT Evaporation from sea times
!                                     leads fraction. Zero over land.
!                                     (kg per square metre per sec).
     &,FQW(row_length,rows,BL_LEVELS)                                   &
                                ! OUT Moisture flux between layers
!                                     (kg per square metre per sec).
!                                     FQW(,1) is total water flux
!                                     from surface, 'E'.
     &,FTL(row_length,rows,BL_LEVELS)                                   &
                                ! OUT FTL(,K) contains net turbulent
!                                     sensible heat flux into layer K
!                                     from below; so FTL(,1) is the
!                                     surface sensible heat, H. (W/m2)
     &,H_SEA(row_length,rows)                                           &
                                ! OUT Surface sensible heat flux over
!                                     sea times leads fraction. (W/m2)
     &,RHOKH(row_length,rows,BL_LEVELS)                                 &
                                ! OUT Exchange coeffs for moisture.
     &,RHOKH_mix(row_length,rows,BL_LEVELS)                             &
                                ! OUT Exchange coeffs for moisture.
! for use in tracer mixing routines
     &,RHOKM_U(row_length,rows,BL_LEVELS)                               &
                                ! OUT Exchange coefficients for u
     &     ,RHOKM_V(row_length,n_rows,BL_LEVELS)                        &
                                ! OUT Exchange coefficients for v
     &,SEA_ICE_HTF(row_length,rows,nice)                                &
                                ! OUT Heat flux through sea-ice
!                                     (W/m2, positive downwards).
     &,SURF_HT_FLUX(row_length,rows)                                    &
                                ! OUT Net downward heat flux at
!                                     surface over land or sea-ice
!                                     fraction of gridbox (W/m2).
     &,SURF_HT_FLUX_LAND(ROW_LENGTH,ROWS)                               &
!                               ! OUT Net downward heat flux at
!                               !     surface over land
!                               !     fraction of gridbox (W/m2).
     &,SURF_HT_FLUX_SICE(ROW_LENGTH,ROWS)                               &
!                               ! OUT Net downward heat flux at
!                               !     surface over sea-ice
!                               !     fraction of gridbox (W/m2).
     &,TAUX(row_length,rows,BL_LEVELS)                                  &
                                ! OUT W'ly component of surface wind
!                                     stress (N/sq m).(On UV-grid with
!                                     first and last rows undefined or
!                                     at present, set to missing data
     &,TAUY(row_length,n_rows,BL_LEVELS)
                                ! OUT S'ly component of surface wind
!                                     stress (N/sq m).  On UV-grid;
!                                     comments as per TAUX.

!  (b) Not passed between lower-level routines (not in workspace at this
!      level) :-

      REAL                                                              &
     & SICE_MLT_HTF(row_length,rows,nice)                               &
                                ! OUT Heat flux due to melting of sea-
!                                   ice (Watts per sq metre).
     &,SNOMLT_SURF_HTF(row_length,rows)                                 &
                                ! OUT Heat flux required for surface
!                                   melting of snow (W/m2).
     &,LATENT_HEAT(row_length,rows)                                     &
                                    ! OUT Surface latent heat flux, +ve
!                                   upwards (Watts per sq m).
     &,Q1P5M(row_length,rows)                                           &
                                ! OUT Q at 1.5 m (kg water per kg air).
     &,T1P5M(row_length,rows)                                           &
                                ! OUT T at 1.5 m (K).
     &,U10M(row_length,rows)                                            &
                                ! OUT U at 10 m (m per s).
     &,V10M(row_length,n_rows)  ! OUT V at 10 m (m per s).

!-2 Genuinely output, needed by other atmospheric routines :-

      REAL                                                              &
     & EI(row_length,rows)                                              &
                                ! OUT Sublimation from lying snow or
!                                   sea-ice (kg/m2/s).
     &,ECAN(row_length,rows)                                            &
                                ! OUT Gridbox mean evaporation from
!                                   canopy/surface store (kg/m2/s).
!                                   Zero over sea.
     &,EXT(land_pts,SM_LEVELS)                                          &
                                ! OUT Extraction of water from each
!                               !     soil layer (kg/m2/s).
     &,SNOWMELT(row_length,rows)                                        &
                                ! OUT Snowmelt (kg/m2/s).
     &,DU(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      BL_LEVELS)                                                  &
                                ! OUT BL increment to u wind field
     &,DV(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      BL_LEVELS)          ! OUT BL increment to u wind field


! OUT additional variables for MOSES II
      REAL                                                              &
     & ESOIL_TILE(land_pts,ntiles)                                      &
                                ! OUT Evaporation from bare soil (kg/m2)
     &,ES(row_length,rows)                                              &
                                ! OUT Surface evapotranspiration from
!                               !     soil moisture store (kg/m2/s).
     &,EI_TILE(LAND_PTS,NTILES)                                         &
                                ! OUT EI for land tiles
     &,Q1P5M_TILE(LAND_PTS,NTILES)                                      &
!                               ! OUT Q1P5M over land tiles.
     &,T1P5M_TILE(LAND_PTS,NTILES)                                      &
!                               ! OUT T1P5M over land tiles.
     &,ECAN_TILE(LAND_PTS,NTILES)                                       &
                                ! OUT ECAN for land tiles
     &,MELT_TILE(LAND_PTS,NTILES)                                       &
!                               ! OUT Snowmelt on tiles (kg/m2/s).
     &,TSTAR_LAND(ROW_LENGTH,ROWS)                                      &
!                               ! OUT Land mean sfc temperature (K)
     &,E_SSI(ROW_LENGTH,ROWS)                                           &
                                ! OUT Mean sea moisture heat flux.
     &,EI_SICE(ROW_LENGTH,ROWS)                                         &
                                ! OUT Sea-ice sublimation rate
!                               !     (sea mean).
     &,FTL_SSI(ROW_LENGTH,ROWS)                                         &
                                ! OUT Mean sea surface heat flux.
     &,TAUX_LAND(ROW_LENGTH,ROWS)                                       &
                                ! OUT W'ly component of land sfc wind
!                               !     stress (N/sq m). (On U-grid
!                               !     with first and last rows
!                               !     undefined or, at present,
!                               !     set to missing data
     &,TAUX_SSI(ROW_LENGTH,ROWS)                                        &
!                               ! OUT W'ly compt of mean sea sfc wind
!                               !     stress (N/sq m). (On U-grid
!                               !     with first and last rows
!                               !     undefined or, at present,
!                               !     set to missing data
     &,TAUY_LAND(ROW_LENGTH,N_ROWS)                                     &
!                               ! OUT S'ly componentt of land sfc wind
!                               !     stress (N/sq m).  On V-grid;
!                               !     comments as per TAUX.
     &,TAUY_SSI(ROW_LENGTH,N_ROWS)
!                               ! OUT S'ly compt of mean sea sfc wind
!                               !     stress (N/sq m).  On V-grid;
!                               !     comments as per TAUX.

      INTEGER                                                           &
     & ERROR            ! OUT 0 - AOK;
!                       !     1 to 7  - bad grid definition detected;


!---------------------------------------------------------------------
!  External routines called :-

      EXTERNAL BDY_IMPL1,SF_IMPL,BDY_IMPL2
      EXTERNAL TIMER

!-----------------------------------------------------------------------
!   Symbolic constants (parameters) reqd in top-level routine :-

!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
!---Soil layer thicknesses (m)
!---6 layers => CABLE else revert to MOSES
   REAL,PARAMETER:: DZSOIL(6) =(/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/)

! C_VKMAN start
      REAL,PARAMETER:: VKMAN=0.4 ! Von Karman's constant
! C_VKMAN end

! Derived local parameters.

      REAL LCRCP,LS,LSRCP

      PARAMETER (                                                       &
     & LCRCP=LC/CP                                                      &
                             ! Evaporation-to-dT conversion factor.
     &,LS=LF+LC                                                         &
                             ! Latent heat of sublimation.
     &,LSRCP=LS/CP                                                      &
                             ! Sublimation-to-dT conversion factor.
     &  )

      INTEGER I, J, K

!  Workspace :-

      REAL                                                              &
     &  DU_NT(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &        BL_LEVELS)                                                &
                            ! non-turbulent increment to u wind field
     &, DV_NT(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,           &
     &        BL_LEVELS)    ! non-turbulen increment to v wind field

!-----------------------------------------------------------------------


      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('IMPS_INTCT ',3)
      ENDIF

! Default solver without substepping
      IF ( .NOT. L_PHYS2_SUBSTEP ) THEN

!      print *,'bef IMP_SOLVER',T_SOIL

! DEPENDS ON: imp_solver
      CALL IMP_SOLVER (                                                 &

! IN mpp variables
     & halo_i, halo_j, off_x, off_y, row_length, rows, n_rows,          &
     & global_row_length, proc_row_group, at_extremity,                 &
     & n_proc, n_procx, n_procy, neighbour,                             &

! IN values defining field dimensions and subset to be processed :
     & NTILES,land_pts,NICE,                                            &

! IN values defining vertical grid of model atmosphere :
     & MODEL_DOMAIN,                                                    &
     & BL_LEVELS,                                                       &
     & r_rho_levels, r_theta_levels,                                    &
     & GAMMA,                                                           &

! IN Substepping information 
     & Substep_Number, Num_Substeps, L_phys2_substep,                   &

! IN New BL solver parameters
     & L_us_blsol, Puns, Pstb,                                          &

! IN U and V momentum fields.
     & U, V,                                                            &
! IN Non turbulent increments to momentum
!  (New dynamics only).
     & R_u, R_v,                                                        &

! IN soil/vegetation/land surface data :
     & LAND_MASK,LAND_INDEX,                                            &
     & ST_LEVELS,SM_LEVELS,TILE_FRAC,CANOPY,                            &
     & FLAND,FLANDG,                                                    &

! IN sea/sea-ice data :
     & DI,ICE_FRACT,DI_NCAT,ICE_FRACT_NCAT,U_0,V_0,                     &

! IN cloud data :
     & Q(1:row_length,1:rows,1:BL_LEVELS),                              &
     & QCF(1:row_length,1:rows,1:BL_LEVELS),                            &
     & QCL(1:row_length,1:rows,1:BL_LEVELS),                            &
     & QCF_latest,QCL_latest, T,                                        &

! IN everything not covered so far :
     & PSTAR,SURF_RADFLUX,TIMESTEP,L_SICE_HEATFLUX,                     &

! IN variables from BDY_LAYR (that used to be local arrays)
     & ALPHA1,ASHTF,DTRDZ_CHARNEY_GRID,RDZ_CHARNEY_GRID,                &
     & DTRDZ_U,DTRDZ_V,RDZ_U,RDZ_V,                                     &
     & FRACA,RHOKH_TILE,SMC,                                            &
     & CHR1P5M,RESFS,Z0HSSI,Z0MSSI,CDR10M_U,CDR10M_V,Z1_TQ,             &
!ajm extra variable added
     & RHOKM_u,RHOKM_v,                                                 &
! needed for new BL numerical solver
     & bl_type_1,bl_type_2,bl_type_3,bl_type_4,                         &
     & bl_type_5,bl_type_6,bl_type_7,                                   &

! IN additional variables for MOSES II
     & TILE_PTS,TILE_INDEX,                                             &
     & L_NEG_TSTAR,                                                     &
     & CANHC_TILE,FLAKE,WT_EXT_TILE,                                    &
     & LW_DOWN,SW_TILE,ALPHA1_SICE,ASHTF_TILE,                          &
     & FQW_ICE,FTL_ICE,RESFT,RHOKH_SICE,RHOKPM,RHOKPM_POT,RHOKPM_SICE,  &
     & Z0H_TILE,Z0M_TILE,CHR1P5M_SICE,                                  &
     & FLANDG_U,FLANDG_V,ANTHROP_HEAT,                                  &

!    EAK
!    IN
     &  l_cable                                                         &
     &, surf_down_sw,alb_tile,cos_zenith_angle,l_tile_pts               &
!     &, surf_down_sw,alb_tile,cos_zenith_angle               &
     &, lat,long,day,time_sec                                           &
     &, ls_rain, ls_snow, conv_rain, conv_snow                          &
     &, SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L                   &
     &, SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE             &
     &, TSOIL_TILE,T_SURF_TILE,HCONS,clapp                              &
     &, SATHH,SATCON,HCAP,HCON,Z1_UV                                    &
     &, smvccl, smvcwt, smvcst, sthf, sthu                              &
     &, slm                                                             &
     &, SOIL_TYPE,VEG_TYPE                                              &
     &, ISNOW_FLG3L,total_nsteps                                        &
     &, FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB                         &
     &, TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB                     &
     &, USTAR_CAB,SURF_HTF_CAB                                          &

     &, TOT_ALB                                                         &
     &, SNAGE_TILE,RTSOIL_TILE                                          &
     &, GFLUX_TILE,SGFLUX_TILE                                          &
     &, GC,GS,CANOPY_GB                                                 & ! Added GC, pfv 25oct13a
! Lestevens 18 March 2010
     &, DIM_CS1, DIM_CS2    &
     &, NPP, NPP_FT         &
     &, GPP, GPP_FT         &
     &, RESP_S, RESP_S_TOT  &
     &, RESP_S_TILE         &  !kdcorbin, 10/10
     &, RESP_P, RESP_P_FT   &
     &, G_LEAF              &  !kdcorbin, 10/10
     &, TRANSP_TILE,        &  !Lestevens 3Nov11
! Lestevens Sept2012: CasaCNP variables 
     &  CPOOL_TILE,NPOOL_TILE,PPOOL_TILE                                &
     &, SOIL_ORDER,GLAI,PHENPHASE                                       &
     ! Lestevens 23apr13
     &, NPP_FT_ACC,RESP_W_FT_ACC,                                       &
 
! INOUT data :
     & T_SOIL,TI,TI_GB,TSTAR,                                           &
     & TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,TSTAR_SSI,                       &
     & TSTAR_TILE, T_latest,Q_latest,                                   &

! INOUT Diagnostics started in BDY_LAYR not requiring STASH flags :
     & E_SEA,FQW,FQW_TILE,EPOT_TILE,FTL,FTL_TILE,H_SEA,RHOKH,           &
     & TAUX,TAUY,                                                       &
     & TAUX_LAND,TAUX_SSI,TAUY_LAND,TAUY_SSI,                           &

! INOUT additional variables for MOSES II
     & SNOW_TILE,LE_TILE,RADNET_SICE,RADNET_TILE,OLR,                   &

! OUT Increments to U and V momentum fields and Tl qw
!  (New dynamics only).
     & DU,DV,                                                           &

! OUT Diagnostic not requiring STASH flags :
     & RHOKH_mix,SEA_ICE_HTF,                                           &

! OUT diagnostics requiring STASH flags :
     & SICE_MLT_HTF,SNOMLT_SURF_HTF,LATENT_HEAT,                        &
     & Q1P5M,T1P5M,U10M,V10M,                                           &

! (IN) STASH flags :-
     & SIMLT,SMLT,SLH,SQ1P5,ST1P5,SU10,SV10,                            &
     & l_ftl, l_fqw, l_taux, l_tauy,                                    &

! OUT additional variables for MOSES II
     & ESOIL_TILE,SURF_HT_FLUX,SURF_HT_FLUX_LAND,SURF_HT_FLUX_SICE,     &
     & EI_TILE,                                                         &
     & Q1P5M_TILE,T1P5M_TILE,ECAN_TILE,MELT_TILE,                       &
     & E_SSI,EI_SICE,FTL_SSI,                                           &
     & ERROR,                                                           &

! OUT data required elsewhere in UM system :
     & ECAN,EI,ES,EXT,SNOWMELT,                                         &
     & lq_mix_bl,                                                       &
! IN SCM namelist data
     & L_FLUX_BC,                                                       &
     & LTIMER                                                           &
     & )

! Call BL solver when substepping is active
      ELSE

! Form non-turbulent term for u, v component. These are basically
! equal to the latest increment due to the convection scheme.

      Do k=1, BL_LEVELS
        Do j=1, ROWS
          Do i=1, ROW_LENGTH
            DU_NT(i,j,k) = (U(i,j,k)+R_u(i,j,k))-U_conv(i,j,k)
          Enddo
        Enddo
        Do j=1, N_ROWS
          Do i=1, ROW_LENGTH
            DV_NT(i,j,k) = (V(i,j,k)+R_v(i,j,k))-V_conv(i,j,k)
          Enddo
        Enddo
      Enddo

!      print *,'impsintc FQW',fqw(:,:,1),fqw_tile
!      print *,'impsintc FQWt',ftl(:,:,1),ftl_tile

! DEPENDS ON: imp_solver
      CALL IMP_SOLVER (                                                 &

! IN mpp variables
     & halo_i, halo_j, off_x, off_y, row_length, rows, n_rows,          &
     & global_row_length, proc_row_group, at_extremity,                 &
     & n_proc, n_procx, n_procy, neighbour,                             &

! IN values defining field dimensions and subset to be processed :
     & NTILES,land_pts,NICE,                                            &

! IN values defining vertical grid of model atmosphere :
     & MODEL_DOMAIN,                                                    &
     & BL_LEVELS,                                                       &
     & r_rho_levels, r_theta_levels,                                    &
     & GAMMA,                                                           &

! IN Substepping info 
     & Substep_Number, Num_Substeps, L_phys2_substep,                   &

! IN New BL solver parameters
     & L_us_blsol, Puns, Pstb,                                          &

! IN U and V momentum fields.
     & U_conv, V_conv,                                                  &
! IN Non turbulent increments to momentum
!  (New dynamics only).
     & DU_NT,DV_NT,                                                     &

! IN soil/vegetation/land surface data :
     & LAND_MASK,LAND_INDEX,                                            &
     & ST_LEVELS,SM_LEVELS,TILE_FRAC,CANOPY,                            &
     & FLAND,FLANDG,                                                    &

! IN sea/sea-ice data :
     & DI,ICE_FRACT,DI_NCAT,ICE_FRACT_NCAT,U_0,V_0,                     &

! IN cloud data :
     & Q_conv,QCF_conv,QCL_conv,QCF_latest,QCL_latest,T_conv,           &

! IN everything not covered so far :
     & PSTAR,SURF_RADFLUX,TIMESTEP,L_SICE_HEATFLUX,                     &

! IN variables from BDY_LAYR (that used to be local arrays)
     & ALPHA1,ASHTF,DTRDZ_CHARNEY_GRID,RDZ_CHARNEY_GRID,                &
     & DTRDZ_U,DTRDZ_V,RDZ_U,RDZ_V,                                     &
     & FRACA,RHOKH_TILE,SMC,                                            &
     & CHR1P5M,RESFS,Z0HSSI,Z0MSSI,CDR10M_U,CDR10M_V,Z1_TQ,             &
!ajm extra variable added
     & RHOKM_u,RHOKM_v,                                                 &
! needed for new BL numerical solver
     & bl_type_1,bl_type_2,bl_type_3,bl_type_4,                         &
     & bl_type_5,bl_type_6,bl_type_7,                                   &

! IN additional variables for MOSES II
     & TILE_PTS,TILE_INDEX,                                             &
     & L_NEG_TSTAR,                                                     &
     & CANHC_TILE,FLAKE,WT_EXT_TILE,                                    &
     & LW_DOWN,SW_TILE,ALPHA1_SICE,ASHTF_TILE,                          &
     & FQW_ICE,FTL_ICE,RESFT,RHOKH_SICE,RHOKPM,RHOKPM_POT,RHOKPM_SICE,  &
     & Z0H_TILE,Z0M_TILE,CHR1P5M_SICE,                                  &
     & FLANDG_U,FLANDG_V,ANTHROP_HEAT,                                  &
!    EAK
!    IN
     &  l_cable                                                         &
     &, surf_down_sw,alb_tile,cos_zenith_angle,l_tile_pts               &
!     &, surf_down_sw,alb_tile,cos_zenith_angle               &
     &, lat,long,day,time_sec                                           &
     &, ls_rain, ls_snow, conv_rain, conv_snow                          &
     &, SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L                   &
     &, SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE             &
     &, TSOIL_TILE,T_SURF_TILE,HCONS,clapp                              &
     &, SATHH,SATCON,HCAP,HCON,Z1_UV                                    &
     &, smvccl, smvcwt, smvcst, sthf, sthu                              &
     &, slm                                                             &
     &, SOIL_TYPE,VEG_TYPE                                              &
     &, ISNOW_FLG3L,total_nsteps                                        &
     &, FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB                         &
     &, TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB                     &
     &, USTAR_CAB,SURF_HTF_CAB                                          &

!
     &, TOT_ALB                                                         &
     &, SNAGE_TILE,RTSOIL_TILE                                          &
     &, GFLUX_TILE,SGFLUX_TILE                                          &
     &, GC,GS,CANOPY_GB                                                 & ! Added GC, pfv 25oct13a
! Lestevens March 2010
! Lest. Jan 2013 - missing dim_* ?
     &, DIM_CS1, DIM_CS2    &
     &, NPP, NPP_FT         &
     &, GPP, GPP_FT         &
     &, RESP_S, RESP_S_TOT  &
     &, RESP_S_TILE         & !kdcorbin, 10/10
     &, RESP_P, RESP_P_FT   &
     &, G_LEAF              & !kdcorbin, 10/10
     &, TRANSP_TILE,        & !Lestevens 3Nov11
! Lestevens Sept2012: CasaCNP variables 
     &  CPOOL_TILE,NPOOL_TILE,PPOOL_TILE                                & 
     &, SOIL_ORDER,GLAI,PHENPHASE                                       &
     ! Lestevens 23apr13
     &, NPP_FT_ACC,RESP_W_FT_ACC,                                       &

! INOUT data :
     & T_SOIL,TI,TI_GB,TSTAR,                                           &
     & TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,TSTAR_SSI,                       &
     & TSTAR_TILE, T_latest,Q_latest,                                   &

! INOUT Diagnostics started in BDY_LAYR not requiring STASH flags :
     & E_SEA,FQW,FQW_TILE,EPOT_TILE,FTL,FTL_TILE,H_SEA,RHOKH,           &
     & TAUX,TAUY,                                                       &
     & TAUX_LAND,TAUX_SSI,TAUY_LAND,TAUY_SSI,                           &

! INOUT additional variables for MOSES II
     & SNOW_TILE,LE_TILE,RADNET_SICE,RADNET_TILE,OLR,                   &

! OUT Increments to U and V momentum fields and Tl qw
!  (New dynamics only).
     & DU,DV,                                                           &

! OUT Diagnostic not requiring STASH flags :
     & RHOKH_mix,SEA_ICE_HTF,                                           &

! OUT diagnostics requiring STASH flags :
     & SICE_MLT_HTF,SNOMLT_SURF_HTF,LATENT_HEAT,                        &
     & Q1P5M,T1P5M,U10M,V10M,                                           &

! (IN) STASH flags :-
     & SIMLT,SMLT,SLH,SQ1P5,ST1P5,SU10,SV10,                            &
     & l_ftl, l_fqw, l_taux, l_tauy,                                    &

! OUT additional variables for MOSES II
     & ESOIL_TILE,SURF_HT_FLUX,SURF_HT_FLUX_LAND,SURF_HT_FLUX_SICE,     &
     & EI_TILE,                                                         &
     & Q1P5M_TILE,T1P5M_TILE,ECAN_TILE,MELT_TILE,                       &
     & E_SSI,EI_SICE,FTL_SSI,                                           &
     & ERROR,                                                           &


! OUT data required elsewhere in UM system :
     & ECAN,EI,ES,EXT,SNOWMELT,                                         &
     & lq_mix_bl,                                                       &
! IN SCM namelist data
     & L_FLUX_BC,                                                       &
     & LTIMER                                                           &
     & )

! Adjust U, V increment to be from Un, Vn as required by the model.
        Do k=1, BL_LEVELS
          Do j=1, ROWS
            Do i=1, ROW_LENGTH
              DU(i,j,k) = DU(i,j,k)+(U_conv(i,j,k)-U(i,j,k))
            Enddo
          Enddo
          Do j=1, N_ROWS
            Do i=1, ROW_LENGTH
              DV(i,j,k) = DV(i,j,k)+(V_conv(i,j,k)-V(i,j,k))
            Enddo
          Enddo
        Enddo

      ENDIF

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('IMPS_INTCT ',4)
      ENDIF

      RETURN
      END SUBROUTINE IMPS_INTCT
