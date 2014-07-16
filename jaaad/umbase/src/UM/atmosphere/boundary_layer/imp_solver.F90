#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE IMP_SOLVER---------------------------------------------
!!!
!!!  Purpose: implicit solver for diffusion equation
!!!           split from bdy_layr routine
!!!
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!  5.2   15/11/00   New Deck         M. Best
!    5.3  25/04/01  Add coastal tiling. Nic Gedney
!!!  5.4  15/05/02  Pass through heights arrays for implicit
!!!                   solver on a sphere   Adrian Lock
!!!  5.5  12/02/03  Pass through seaice catagories. J.Ridley
!!!  5.5  19/02/03  Remove redundent L_BL_LSPICE. Adrian Lock
!!!  6.1  01/09/04  Pass potential evaporation related variables.
!                                                          Nic Gedney
!  6.1  17/05/04  Changes to the BL solver to enable phys2 substepping.
!                                                       M. Diamantakis
!!!  6.2  21/03/05  Pass through implicit scheme weights.
!!!                                       M. Diamantakis
!!!
! 6.2      21/02/06    Switch for mixing ratios   A.P.Lock
!!!  6.2  02/02/06  Passes L_flux_bc through argument list to allow
!!!                 settings for prescribed surface flux forcing
!!!                                                           R. Wong
!    6.4  10/01/07  Introduce unconditionally stable and non-oscillatory
!                    BL numerical solver    M. Diamantakis 
!!!
!!! Programming standard : unified model documentation paper No 3
!!!
!!!  Documentation: UMDP 24.
!!!
!!!---------------------------------------------------------------------

!    Arguments :-
      SUBROUTINE IMP_SOLVER (                                           &

! IN MPP variables
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
     & U, V,                                                            &
! IN Non turbulent increments to momentum
!  (New dynamics only).
     & DU_NT,DV_NT,                                                     &

! IN soil/vegetation/land surface data :
     & LAND_MASK,LAND_INDEX,                                            &
     & ST_LEVELS,SM_LEVELS,TILE_FRAC,CANOPY,                            &
     & FLAND,FLANDG,                                                    &

! IN sea/sea-ice data :
     & DI, ICE_FRACT, DI_NCAT, ICE_FRACT_NCAT, U_0, V_0,                &

! IN cloud data :
     & Q,QCF,QCL,QCF_latest,QCL_latest, T,                              &

! IN everything not covered so far :
     & PSTAR,RAD_SICE,TIMESTEP,L_SICE_HEATFLUX,                         &

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
! IN SCM LOGICAL
     & L_FLUX_BC,                                                       &
     & LTIMER                                                           &
     & )

      IMPLICIT NONE

!  Inputs :-

! (a) Defining horizontal grid and subset thereof to be processed.

      INTEGER                                                           &
! cjj additions - MPP variables.
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
     &,NICE                      ! IN No.of sea ice catagories

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
     & L_phys2_substep

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
                                 ! IN vertical co-ordinates
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j,0:bl_levels)                &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j, bl_levels)

      REAL                                                              &
     &  GAMMA(bl_levels)                                                &
     &, DU_NT(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &        BL_LEVELS)                                                &
                            ! non-turbulent increment to u wind field
     &, DV_NT(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,           &
     &        BL_LEVELS)                                                &
                            ! non-turbulen increment to v wind field
     &, U(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      BL_LEVELS)                                                  &
     &, V(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      BL_LEVELS)


! (c) Soil/vegetation/land surface parameters (mostly constant).

      LOGICAL                                                           &
     & LAND_MASK(row_length,rows)  ! IN T if land, F elsewhere.

#include "nstypes.h"


      INTEGER                                                           &
     & LAND_INDEX(land_pts)      ! IN LAND_INDEX(I)=J => the Jth
!                                     point in P_FIELD is the Ith
!                                     land point.

      INTEGER                                                           &
     & ST_LEVELS                                                        &
                                   ! IN No. of deep soil temp. levels
     &,SM_LEVELS             &     ! IN No. of soil moisture levels
! Lestevens March 2010
     &, DIM_CS1 &
     &, DIM_CS2

      REAL                                                              &
     & TILE_FRAC(land_pts,ntiles)                                       &
                                ! IN fractional coverage for each
                                !    surface tile
     &,CANOPY(land_pts,ntiles)                                          &
                                ! IN Surface/canopy water (kg/m2)
     &,FLAND(LAND_PTS)                                                  &
                                ! IN Land fraction on land tiles.
     &,FLANDG(ROW_LENGTH,ROWS)
!                               ! IN Land fraction on all tiles.

! (d) Sea/sea-ice data.

      REAL                                                              &
     & DI(row_length,rows)                                              &
                                   ! IN "Equivalent thickness" of
!                                     sea-ice(m).
     &,ICE_FRACT(row_length,rows)                                       &
                                   ! IN Fraction of gridbox covered by
!                                     sea-ice (decimal fraction).
     &,DI_NCAT(row_length,rows,nice)                                    &
                                      ! IN "Equivalent thickness" of
!                                     sea-ice on catagories(m).
     &,ICE_FRACT_NCAT(row_length,rows,nice)                             &
                                             ! IN Fraction of gridbox
!                                    covered by sea-ice catagory.
     &,U_0(row_length,rows)                                             & 
                                   ! IN W'ly component of surface
!                                     current (m/s).
     &,V_0(row_length,n_rows)      ! IN S'ly component of surface
!                                     current (m/s).

! (e) Cloud data.

      REAL                                                              &
     & QCF(row_length,rows,BL_LEVELS)                                   &
                                       ! IN Cloud ice (kg per kg air)
     &,QCL(row_length,rows,BL_LEVELS)                                   &
                                       ! IN Cloud liquid water
     &,Q(1,row_length,rows,BL_LEVELS)                                   &
                                       ! IN specific humidity
     &,T(row_length,rows,BL_LEVELS)                                     &
                                       ! IN temperature
                          ! Latest estimates to time level n+1 values
     &,QCF_latest(row_length,rows,BL_LEVELS)                            &
                                ! IN Cloud ice (kg per kg air)
     &,QCL_latest(row_length,rows,BL_LEVELS) ! IN Cloud liquid water

! (f) Atmospheric + any other data not covered so far, incl control.

      REAL                                                              &
     & PSTAR(row_length,rows)                                           &
                                   ! IN Surface pressure (Pascals).
     &,RAD_SICE(ROW_LENGTH,ROWS)                                        &
                                   ! IN Surface net SW and downward LW
!                                  !    radiation on sea-ice fraction
!                                  !    (W/sq m, positive downwards).
     &,TIMESTEP                    ! IN Timestep (seconds).

      LOGICAL                                                           &
     & L_SICE_HEATFLUX             ! IN T: semi-implicit sea-ice temp

! IN variables from BDY_LAYR (that used to be local arrays)
      REAL                                                              &
     & ALPHA1(land_pts,ntiles)                                          &
                               ! IN Mean gradient of saturated
!                                 specific humidity with
!                                 respect to temperature between
!                                 the bottom model layer and the
!                                 tile surfaces.
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
!                               ! IN 1/(Z_U(K)-Z_U(K-1)) m^{-1}
     &,RDZ_V(row_length,n_rows,2:BL_LEVELS)                             &
                                ! IN 1/(Z_V(K)-Z_V(K-1)) m^{-1}
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
     &,CDR10M_U(row_length,rows)                                        &
                                ! IN Ratio of CD's reqd for calculation
     &,CDR10M_V(row_length,n_rows)                                      &
                                ! IN Ratio of CD's reqd for calculation
     &,Z1_TQ(row_length,rows)                                           &
                                ! IN Height of lowest theta level.
!ajm  extra variable added
     &,RHOKM_U(row_length,rows,BL_LEVELS)                               &
                                ! IN Exchange coefficients for u
     &,RHOKM_V(row_length,n_rows,BL_LEVELS)                             &
                                ! IN Exchange coefficients for v
     &,ANTHROP_HEAT(NTILES)
                                ! IN Additional heat source on tiles
                                !    for anthropogenic urban heat
                                !    source (W/m2)

! IN additional variables for MOSES II
      LOGICAL                                                           &
     & L_NEG_TSTAR                 ! IN Switch for -ve TSTAR error check

      INTEGER                                                           &
     & TILE_PTS(NTYPE)                                                  &
                                 ! IN Number of tile points.
     &,TILE_INDEX(LAND_PTS,NTYPE)
!                                ! IN Index of tile points.

      REAL                                                              &
     & CANHC_TILE(LAND_PTS,NTILES)                                      &
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
     &,ALPHA1_SICE(ROW_LENGTH,ROWS)                                     &
!                               ! IN ALPHA1 for sea-ice.
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
     &,FLANDG_U(ROW_LENGTH,ROWS)                                        &
!                               ! IN Land frac (on U-grid, with 1st
!                               !    and last rows undefined or, at
!                               !    present, set to "missing data")
     &,FLANDG_V(ROW_LENGTH,N_ROWS)
!                               ! IN Land frac (on V-grid, with 1st
!                               !    and last rows undefined or, at
!                               !    present, set to "missing data")


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
!                                      temperature in catagory (K).
     &,TI_GB(row_length,rows)                                           &
                                  ! INOUT Ice temperature GBM (K).
     &,TSTAR(row_length,rows)                                           &
                                ! INOUT Surface temperature (K).
     &,TSTAR_LAND(ROW_LENGTH,ROWS)                                      &
                                   ! OUT   Land mean sfc temperature (K)
     &,TSTAR_SEA(ROW_LENGTH,ROWS)                                       &
                                   ! IN    Open sea sfc temperature (K).
     &,TSTAR_SICE(ROW_LENGTH,ROWS)                                      &
                                   ! INOUT Sea-ice sfc temperature (K).
     &,TSTAR_SSI(ROW_LENGTH,ROWS)                                       &
                                   ! INOUT Sea mean sfc temperature (K).
     &,TSTAR_TILE(land_pts,ntiles)
                                ! INOUT Surface tile temperature
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
     &, Z1_UV(row_length,rows)                                          &
                                ! Height of lowest u,v level.
     &, HCONS(land_pts)
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
!
     &,TOT_ALB(LAND_PTS,NTILES)                                      &
     &,CANOPY_GB(LAND_PTS)                                           &  
     &,GC(LAND_PTS,NTILES)                                           & ! Added GC, pfv 25oct13a
     &,GS(LAND_PTS)                                                  &  
! Lestevens 17 Feb 2010 - passing CO2 fluxes
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
     &,RESP_P_FT(LAND_PTS,NTILES)                                    &
     &,G_LEAF(LAND_PTS,NTILES)  & !kdcorbin, 10/10
     &,TRANSP_TILE(LAND_PTS,NTILES)  & !Lestevens 3Nov11
! Lestevens Sept2012: CasaCNP variables 
     &,CPOOL_TILE(LAND_PTS,NTILES,10) &
     &,NPOOL_TILE(LAND_PTS,NTILES,10) &
     &,PPOOL_TILE(LAND_PTS,NTILES,12) &
     &,SOIL_ORDER(LAND_PTS)           &
     &,GLAI(LAND_PTS,NTILES)          &
     &,PHENPHASE(LAND_PTS,NTILES)     &
     ! Lestevens 23apr13
     &,NPP_FT_ACC(land_pts,NTILES)    &
     &,RESP_W_FT_ACC(land_pts,NTILES)

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


!  In/Outputs :-
!-1 Diagnostic (or effectively so - includes coupled model requisites):-

!  (a) Calculated anyway (use STASH space from higher level) :-
!
      REAL                                                              &
     & E_SEA(row_length,rows)                                           &
                                ! INOUT Evaporation from sea times
!                                     leads fraction. Zero over land.
!                                     (kg per square metre per sec).
     &,FQW(row_length,rows,BL_LEVELS)                                   &
                                ! INOUT Moisture flux between layers
!                                     (kg per square metre per sec).
!                                     FQW(,1) is total water flux
!                                     from surface, 'E'.
     &,FQW_TILE(land_pts,ntiles)                                        &
                                ! INOUT surface tile moisture flux
     &,EPOT_TILE(land_pts,ntiles)                                       &
!                               ! INOUT surface tile potential
!                               !       evaporation
     &,E_SSI(ROW_LENGTH,ROWS)                                           &
                                ! OUT   Surface FQW for mean sea.
     &,EI_SICE(ROW_LENGTH,ROWS)                                         &
                                ! OUT   Sea-ice sumblimation
!                               !       (sea mean).
     &,FTL(row_length,rows,BL_LEVELS)                                   &
                                ! INOUT FTL(,K) contains net turbulent
!                                     sensible heat flux into layer K
!                                     from below; so FTL(,1) is the
!                                     surface sensible heat, H. (W/m2)
     &,FTL_TILE(land_pts,ntiles)                                        &
                                ! INOUT surface tile heat flux
     &,FTL_SSI(row_length,rows)                                         &
                                ! OUT sea mean surface heat flux
     &,H_SEA(row_length,rows)                                           &
                                ! INOUT Surface sensible heat flux over
!                                     sea times leads fraction. (W/m2)
     &,RHOKH(row_length,rows,BL_LEVELS)                                 &
                                ! INOUT Exchange coeffs for moisture.
     &,TAUX(row_length,rows,BL_LEVELS)                                  &
                                ! INOUT W'ly component of surface wind
!                                     stress (N/sq m).(On UV-grid with
!                                     first and last rows undefined or
!                                     at present, set to missing data
     &,TAUX_LAND(ROW_LENGTH,ROWS)                                       &
!                                ! INOUT W'ly component of land sfc wind
!                                !     stress (N/sq m). (On U-grid
!                                !     with first and last rows
!                                !     undefined or, at present,
!                                !     set to missing data
     &,TAUX_SSI(ROW_LENGTH,ROWS)                                        &
!                                ! INOUT W'ly compt of mean sea sfc wind
!                                !     stress (N/sq m). (On U-grid
!                                !     with first and last rows
!                                !     undefined or, at present,
!                                !     set to missing data
     &,TAUY(row_length,n_rows,BL_LEVELS)                                &
                                ! INOUT S'ly component of surface wind
!                                     stress (N/sq m).  On UV-grid;
!                                     comments as per TAUX.
     &,TAUY_LAND(ROW_LENGTH,N_ROWS)                                     &
!                                ! INOUT S'ly component of land sfc wind
!                                !     stress (N/sq m).  On V-grid;
!                                !     comments as per TAUX.
     &,TAUY_SSI(ROW_LENGTH,N_ROWS)
!                                ! INOUT S'ly compt of mean sea sfc wind
!                                !     stress (N/sq m).  On V-grid;
!                                !     comments as per TAUX.


! INOUT additional variables for MOSES II
      REAL                                                              &
     & SNOW_TILE(LAND_PTS,NTILES)                                       &
!                               ! INOUT Snow on tiles (kg/m2).
     &,LE_TILE(LAND_PTS,NTILES)                                         &
                                ! INOUT Surface latent heat flux for
!                               !       land tiles (W/m2).
     &,RADNET_SICE(ROW_LENGTH,ROWS)                                     &
!                               ! INOUT Sea-ice surface net radiation.
     &,RADNET_TILE(LAND_PTS,NTILES)                                     &
!                               ! INOUT Tile surface net radiation.
     &,OLR(ROW_LENGTH,ROWS)     ! IN    TOA - surface upward LW on
!                               !       last radiation timestep
!                               ! OUT   Corrected TOA outward LW


!  Outputs :-

! (New dynamics only)
      REAL                                                              &
     & DU(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      BL_LEVELS)                                                  &
                                ! OUT BL increment to u wind field
     &,DV(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      BL_LEVELS)          ! OUT BL increment to u wind field


!-1 Diagnostic (or effectively so - includes coupled model requisites):-

!  (a) Calculated anyway (use STASH space from higher level) :-
!
      REAL                                                              &
     & RHOKH_mix(row_length,rows,BL_LEVELS)                             &
                                ! OUT Exchange coeffs for moisture.
! for use in tracer mixing routines
     &,SEA_ICE_HTF(row_length,rows,nice)
                                ! OUT Heat flux through sea-ice
!                                     (W/m2, positive downwards).

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
     &,ES(row_length,rows)                                              &
                                ! OUT Surface evapotranspiration from
!                               !     soil moisture store (kg/m2/s).
     &,EXT(land_pts,SM_LEVELS)                                          &
                                ! OUT Extraction of water from each
!                               !     soil layer (kg/m2/s).
     &,SNOWMELT(row_length,rows)! OUT Snowmelt (kg/m2/s).


! OUT additional variables for MOSES II
      REAL                                                              &
     & ESOIL_TILE(land_pts,ntiles)                                      &
                                ! OUT Evaporation from bare soil (kg/m2)
     &,SURF_HT_FLUX(row_length,rows)                                    &
!                               ! OUT Net downward heat flux at surface
!                                     over land and sea-ice fraction of
!                                     gridbox (W/m2).
     &,SURF_HT_FLUX_LAND(ROW_LENGTH,ROWS)                               &
!                               ! OUT Net downward heat flux at
!                               !     surface over land
!                               !     fraction of gridbox (W/m2).
     &,SURF_HT_FLUX_SICE(ROW_LENGTH,ROWS)                               &
!                               ! OUT Net downward heat flux at
!                               !     surface over sea-ice
!                               !     fraction of gridbox (W/m2).
     &,EI_TILE(LAND_PTS,NTILES)                                         &
                                ! OUT EI for land tiles
     &,Q1P5M_TILE(LAND_PTS,NTILES)                                      &
!                               ! OUT Q1P5M over land tiles.
     &,T1P5M_TILE(LAND_PTS,NTILES)                                      &
!                               ! OUT T1P5M over land tiles.
     &,ECAN_TILE(LAND_PTS,NTILES)                                       &
                                ! OUT ECAN for land tiles
     &,MELT_TILE(LAND_PTS,NTILES)
!                               ! OUT Snowmelt on tiles (kg/m2/s).

      INTEGER                                                           &
     & ERROR                    ! OUT 0 - AOK;
!                               ! 1 to 7  - bad grid definition detected;


      LOGICAL LTIMER               ! Logical switch for TIMER diags
      LOGICAL L_FLUX_BC            ! Logical for prescribing surface
                                   ! flux forcing in the SCM

!---------------------------------------------------------------------
!  External routines called :-

      EXTERNAL BDY_IMPL1,SF_IMPL,BDY_IMPL2
      EXTERNAL TIMER

!-----------------------------------------------------------------------
!   Symbolic constants (parameters) reqd in top-level routine :-

#include "c_r_cp.h"
#include "c_g.h"
#include "c_lheat.h"
#include "soil_thick.h"
#include "c_vkman.h"

! Derived local parameters.

      REAL LCRCP,LS,LSRCP
      REAL Pnonl,P1,P2,I1,E1,E2 ! parameters for new BL solver 
      REAL SQRT2                ! SQRT(2.)
      
      PARAMETER (                                                       &
     & LCRCP=LC/CP                                                      &
                             ! Evaporation-to-dT conversion factor.
     &,LS=LF+LC                                                         &
                             ! Latent heat of sublimation.
     &,LSRCP=LS/CP                                                      &
                             ! Sublimation-to-dT conversion factor.
     &  )

      REAL                                                              &
     &  GAMMA1(row_length,rows)                                         &
     & ,GAMMA2(row_length,rows)

!-----------------------------------------------------------------------

!  Workspace :-

      REAL                                                              &
     & QW(row_length, rows, BL_LEVELS)                                  &
                                        ! LOCAL total water
     &,TL(row_length, rows, BL_LEVELS)                                  &
                                        ! LOCAL liquid water temperature
     &,DQW(ROW_LENGTH,ROWS,BL_LEVELS)                                   &
                                        ! LOCAL BL increment to q field
     &,DTL(ROW_LENGTH,ROWS,BL_LEVELS)                                   &
                                        ! LOCAL BL increment to T field
! DU_STAR, DV_STAR: 1st stage Temporary BL incr to u, v wind 
! components from new (stable) BL solver.
     &,DU_STAR(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         BL_LEVELS)                                               &
     &,DV_STAR(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,          &
     &         BL_LEVELS)                                               & 
     &,DQW_NT(ROW_LENGTH,ROWS,BL_LEVELS)                                &
                                        ! OUT NT incr to qw
     &,DTL_NT(ROW_LENGTH,ROWS,BL_LEVELS)                                &
                                        ! OUT NT incr to TL
     &,CT_CTQ(ROW_LENGTH,ROWS,BL_LEVELS)                                &
                                        ! LOCAL Coefficient in T and q
!                                       !       tri-diagonal implicit
!                                       !       matrix
     &,CQ_CM_U(ROW_LENGTH,ROWS,BL_LEVELS)                               &
!                                       ! LOCAL Coefficient in U
!                                       !       tri-diagonal implicit
!                                       !       matrix
     &,CQ_CM_V(ROW_LENGTH,ROWS,BL_LEVELS)                               &
!                                       ! LOCAL Coefficient in V
!                                       !       tri-diagonal implicit
!                                       !       matrix
!
     &,E_LAND(ROW_LENGTH,ROWS)                                          &
                                        ! LOCAL FQW over mean land
     &,EI_LAND(ROW_LENGTH,ROWS)                                         &
                                        ! LOCAL EI over mean land
     &,FTL_LAND(ROW_LENGTH,ROWS)                                        
                                        ! LOCAL FTL over mean land

! The three set of arrays below are needed by the uncond stable 
! BL numerical solver

      Real, DIMENSION(:,:,:), ALLOCATABLE ::                            &
     &  DQW1, DTL1, CTCTQ1 
!
! Following arrays: stable & non-oscillatory solver fluxes for 1st
!                   stage (predictor)
!
      Real, DIMENSION(:,:,:), ALLOCATABLE ::                            &
     &  FQW_star, FTL_star, TAUX_star, TAUY_star 

!
! Automatic arrays for land and sea surface stress diagnostics (used 
! as inputs to the ocean when the coupled model is selected)
!
      Real, DIMENSION(ROW_LENGTH,ROWS) ::                               &
        TAUX_LAND_star, TAUX_SSI_star

      Real, DIMENSION(ROW_LENGTH,N_ROWS) ::                             &
     &  TAUY_LAND_star, TAUY_SSI_star


      INTEGER I,J,K,L,N

      LOGICAL L_CORRECT 

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('IMPSOLV ',3)
      ENDIF

! Use standard BL solver
      IF ( .NOT. L_us_blsol ) Then 

!      print *,'impsolver FQW',T_SOIL,fqw(:,:,1),fqw_tile
!      print *,'impsolver FQWt',ftl(:,:,1),ftl_tile

! DEPENDS ON: bdy_impl1
      CALL BDY_IMPL1 (                                                  &
     & HALO_I,HALO_J,OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,BL_LEVELS,      &
     &   L_phys2_substep,                                               &
     & r_rho_levels,r_theta_levels,                                     &
     & Q,QCL,QCF,Q_LATEST,QCL_LATEST,QCF_LATEST,                        &
     & T,T_LATEST,                                                      &
     & DTRDZ_CHARNEY_GRID,DTRDZ_U,DTRDZ_V,                              &
     & RHOKH,RHOKM_U,RHOKM_V,                                           &
     & RDZ_CHARNEY_GRID,RDZ_U,RDZ_V,GAMMA,                              &
     & DU_NT,DV_NT,                                                     &
     & FQW,FTL,TAUX,TAUY,                                               &
     & QW,TL,                                                           &
     & CT_CTQ,DQW,DTL,CQ_CM_U,CQ_CM_V,DU,DV,                            &
     & LTIMER                                                           &
     &  )

!      print *,'impsolver FQW 1',fqw(:,:,1),fqw_tile
!      print *,'bef sf_impl GAMMA(1)', GAMMA(1)

! DEPENDS ON: sf_impl
      CALL SF_IMPL (                                                    &

! IN values defining field dimensions and subset to be processed :
     & OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,LAND_PTS,                     &

! IN soil/vegetation/land surface data :
     & LAND_INDEX,LAND_MASK,NICE,                                       &
     & NTILES,TILE_INDEX,TILE_PTS,SM_LEVELS,                            &
     & CANHC_TILE,CANOPY,FLAKE,SMC,                                     &
     & TILE_FRAC,WT_EXT_TILE,                                           &
     & FLAND,FLANDG,                                                    &

! IN sea/sea-ice data :
     & DI,ICE_FRACT,DI_NCAT,ICE_FRACT_NCAT,U_0,V_0,                     &

! IN everything not covered so far :
     & PSTAR,LW_DOWN,RAD_SICE,SW_TILE,TIMESTEP,                         &
     & T_SOIL,QW,TL,U,V,                                                &
     & RHOKM_U,RHOKM_V,GAMMA(1),                                        &
     & ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,                             &
     & DTRDZ_CHARNEY_GRID,DU,DV,                                        &
     & FQW_TILE,EPOT_TILE,FQW_ICE,FTL_ICE,                              &
     & FRACA,RESFS,RESFT,RHOKH,RHOKH_TILE,RHOKH_SICE,RHOKPM,RHOKPM_POT, &
     & RHOKPM_SICE,Z1_TQ,Z0HSSI,Z0MSSI,Z0H_TILE,Z0M_TILE,               &
     & CDR10M_U,CDR10M_V,CHR1P5M,CHR1P5M_SICE,CT_CTQ,DQW,DTL,           &
     & CQ_CM_U,CQ_CM_V,                                                 &
     & L_NEG_TSTAR,                                                     &
     & FLANDG_U,FLANDG_V,ANTHROP_HEAT,L_SICE_HEATFLUX,                  &

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

! IN STASH flags :-
     & SIMLT,SMLT,SLH,SQ1P5,ST1P5,SU10,SV10,                            &

! INOUT data :
     & TI,TI_GB,TSTAR,                                                  &
     & TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,TSTAR_SSI,                       &
     & TSTAR_TILE,SNOW_TILE,                                            &
     & LE_TILE,RADNET_SICE,RADNET_TILE,                                 &
     & E_SEA,FQW,FTL,FTL_TILE,H_SEA,OLR,                                &
     & TAUX,TAUY,                                                       &
     & TAUX_LAND,TAUX_SSI,TAUY_LAND,TAUY_SSI,                           &

! OUT Diagnostic not requiring STASH flags :
     & ECAN,EI_TILE,ESOIL_TILE,                                         &
     & SEA_ICE_HTF,SURF_HT_FLUX,SURF_HT_FLUX_LAND,SURF_HT_FLUX_SICE,    &

! OUT diagnostic requiring STASH flags :
     & SICE_MLT_HTF,SNOMLT_SURF_HTF,LATENT_HEAT,                        &
     & Q1P5M,Q1P5M_TILE,T1P5M,T1P5M_TILE,U10M,V10M,                     &

! OUT data required elsewhere in UM system :
     & ECAN_TILE,EI,ES,EXT,SNOWMELT,MELT_TILE,RHOKH_MIX,                &
     & ERROR,                                                           &

! LOGICAL LTIMER
     & lq_mix_bl,                                                       &
     & L_FLUX_BC,                                                       &
     & LTIMER                                                           &
     & )


! DEPENDS ON: bdy_impl2
      CALL BDY_IMPL2 (                                                  &

! IN values defining field dimensions and subset to be processed :
     & OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,                              &

! IN values defining vertical grid of model atmosphere :
     & BL_LEVELS,                                                       &

! IN Diagnostic switches
     &   l_ftl, l_fqw, l_taux, l_tauy,                                  &
! IN data :
     & GAMMA,                                                           &
     & RHOKH,RHOKM_U,RHOKM_V,RDZ_CHARNEY_GRID,RDZ_U,RDZ_V,              &

! INOUT data :
     & QW,TL,FQW,FTL,TAUX,TAUY,                                         &
     & DU,DV,CT_CTQ,DQW,DTL,CQ_CM_U,CQ_CM_V,                            &

! OUT data :
     & T_LATEST,Q_LATEST,RHOKH_MIX,                                     &

! LOGICAL LTIMER
     & LTIMER                                                           &
     &  )

      ELSE
!
! Use new unconditionally stable and non-oscillatory BL solver
!
      Allocate (DQW1(ROW_LENGTH,ROWS,BL_LEVELS))
      Allocate (DTL1(ROW_LENGTH,ROWS,BL_LEVELS))   
      Allocate (CTCTQ1(ROW_LENGTH,ROWS,BL_LEVELS))
!
! Following arrays: stable & non-oscillatory solver fluxes for 1st
!                   stage (predictor)
!
      Allocate (FQW_star(row_length,rows,BL_LEVELS))
      Allocate (FTL_star(row_length,rows,BL_LEVELS))
      Allocate (TAUX_star(row_length,rows,BL_LEVELS))
      Allocate (TAUY_star(row_length,n_rows,BL_LEVELS))
!
!----------------------------------------------------------------------
! Compute 1st stage solution (predictor).
!----------------------------------------------------------------------
!
! First compute the scheme coefficients for the 1st stage. Make 
! coefficients dependent on the BL type for achieving better balance 
! between stability-accuracy: stable BL can be strongly nonlinear and 
! stiff and thus numerically unstable, so choose a large P value. 
! Unstable BL are weakly nonlinear so the solver should be able to cope 
! with small P. 
!
!----------------------------------------------------------------------
      SQRT2 = SQRT(2.)

      DO J=1, ROWS
        DO I=1, ROW_LENGTH
          P1=bl_type_1(i,j)*Pstb+(1.-bl_type_1(i,j))*Puns
          P2=bl_type_2(i,j)*Pstb+(1.-bl_type_2(i,j))*Puns
          Pnonl=max(P1,P2)
          I1 = (1.+1./SQRT2)*(1.+Pnonl)
          E1 = (1.+1./SQRT2)*( Pnonl + (1./SQRT2) +                     &
     &                        SQRT(Pnonl*(SQRT2-1.)+0.5) )
          GAMMA1(I,J) = I1
          GAMMA2(I,J) = I1 - E1
        END DO
      END DO

      L_CORRECT = .FALSE.

! DEPENDS ON: BDY_IMPL3
      CALL BDY_IMPL3 (                                                  &
     & HALO_I,HALO_J,OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,BL_LEVELS,      &
     & L_CORRECT,r_rho_levels,r_theta_levels,                           &
     & Q,QCL,QCF,Q_LATEST,QCL_LATEST,QCF_LATEST,T,T_LATEST,             &
     & DTRDZ_CHARNEY_GRID,DTRDZ_U,DTRDZ_V,                              &
     & RHOKH(1,1,2),RHOKM_U(1,1,2),RHOKM_V(1,1,2),                      &
     & RDZ_CHARNEY_GRID,RDZ_U,RDZ_V,GAMMA1,GAMMA2,GAMMA,                &
     & DU_NT,DV_NT,DQW_NT,DTL_NT,FQW,FTL,TAUX,TAUY,                     &
     & QW,TL,DQW1,DTL1,CT_CTQ,CTCTQ1,DQW,DTL,CQ_CM_U,CQ_CM_V,DU,DV,     &
     & LTIMER                                                           &
     &  )

! DEPENDS ON: SF_IMPL2 
      CALL SF_IMPL2 (                                                   &

! IN values defining field dimensions and subset to be processed :
     & OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,LAND_PTS,                     &

! IN soil/vegetation/land surface data :
     & LAND_INDEX,LAND_MASK,NICE,NTILES,TILE_INDEX,TILE_PTS,SM_LEVELS,  &
     & CANHC_TILE,CANOPY,FLAKE,SMC,TILE_FRAC,WT_EXT_TILE,FLAND,FLANDG,  &

! IN sea/sea-ice data :
     & DI,ICE_FRACT,DI_NCAT,ICE_FRACT_NCAT,U_0,V_0,                     &

! IN everything not covered so far :
     & PSTAR,LW_DOWN,RAD_SICE,SW_TILE,TIMESTEP,T_SOIL,QW,TL,U,V,        &
     & RHOKM_U,RHOKM_V,GAMMA(1),GAMMA1,GAMMA2,ALPHA1,ALPHA1_SICE,       &
     & ASHTF,ASHTF_TILE,DTRDZ_CHARNEY_GRID,DU,DV,FQW_TILE,EPOT_TILE,    &
     & FQW_ICE,FTL_ICE,FRACA,RESFS,RESFT,RHOKH,RHOKH_TILE,RHOKH_SICE,   & 
     & RHOKPM,RHOKPM_POT,RHOKPM_SICE,Z1_TQ,                             &
     & Z0HSSI,Z0MSSI,Z0H_TILE,Z0M_TILE,                                 &
     & CDR10M_U,CDR10M_V,CHR1P5M,CHR1P5M_SICE,CT_CTQ,CTCTQ1,DQW,DTL,    &
     & DQW1,DTL1,DU_STAR,DV_STAR,CQ_CM_U,CQ_CM_V,                       &
     & L_NEG_TSTAR,L_CORRECT,FLANDG_U,FLANDG_V,ANTHROP_HEAT,            &
     & L_SICE_HEATFLUX,                                                 &

! IN STASH flags :-
     & SIMLT,SMLT,SLH,SQ1P5,ST1P5,SU10,SV10,                            &

! INOUT data :
     & TI,TI_GB,TSTAR,TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,TSTAR_SSI,        &
     & TSTAR_TILE,SNOW_TILE,LE_TILE,RADNET_SICE,RADNET_TILE,E_SEA,      &
     & FQW,FTL,FTL_TILE,H_SEA,OLR,TAUX,TAUY,                            &
     & TAUX_LAND,TAUX_LAND_star,TAUX_SSI,TAUX_SSI_star,TAUY_LAND,       &
     & TAUY_LAND_star,TAUY_SSI,TAUY_SSI_star,                           &

! OUT Diagnostic not requiring STASH flags :
     & ECAN,EI_TILE,ESOIL_TILE,                                         &
     & SEA_ICE_HTF,SURF_HT_FLUX,SURF_HT_FLUX_LAND,SURF_HT_FLUX_SICE,    &

! OUT diagnostic requiring STASH flags :
     & SICE_MLT_HTF,SNOMLT_SURF_HTF,LATENT_HEAT,                        &
     & Q1P5M,Q1P5M_TILE,T1P5M,T1P5M_TILE,U10M,V10M,                     &

! OUT data required elsewhere in UM system :
     & ECAN_TILE,EI,ES,EXT,SNOWMELT,MELT_TILE,RHOKH_MIX,ERROR,          &
     & lq_mix_bl, L_FLUX_BC, LTIMER                                     &
     & )

! DEPENDS ON: BDY_IMPL4
      CALL BDY_IMPL4 (                                                  &

! IN values defining field dimensions and subset to be processed :
     & HALO_I,HALO_J,OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,BL_LEVELS,      &

! IN Diagnostic switches
     & L_CORRECT, l_ftl, l_fqw, l_taux, l_tauy,                         &

! IN data :
     & GAMMA1,GAMMA2,R_RHO_LEVELS,RHOKH(1,1,2),RHOKM_U(1,1,2),          &
     & RHOKM_V(1,1,2),RDZ_CHARNEY_GRID,DTRDZ_CHARNEY_GRID,RDZ_U,RDZ_V,  &

! INOUT data :
     & QW,TL,FQW,FTL,TAUX,TAUY,FQW_star,FTL_star,TAUX_star,TAUY_star,   &
     & DU,DV,DU_STAR,DV_STAR,CT_CTQ,CTCTQ1,DQW,DTL,                     &
     & DQW_NT,DTL_NT,CQ_CM_U,CQ_CM_V,                                   &

! OUT data:
     & T_LATEST,Q_LATEST,RHOKH_MIX,                                     &

! LOGICAL LTIMER
     & LTIMER                                                           &
     &  )
!
!----------------------------------------------------------------------
! Compute 2nd stage (final) solution (corrector).
!----------------------------------------------------------------------
!
! First compute the scheme coefficients for the 2nd stage. Make 
! coefficients dependent on the BL type for achieving better balance 
! between stability-accuracy: stable BL can be strongly nonlinear and 
! stiff and thus numerically unstable, so choose a large P value. 
! Unstable BL are weakly nonlinear so the solver should be able to cope 
! with small P. 
!
!----------------------------------------------------------------------
      DO j=1, ROWS
        DO i=1, ROW_LENGTH
          P1=bl_type_1(i,j)*Pstb+(1.-bl_type_1(i,j))*Puns
          P2=bl_type_2(i,j)*Pstb+(1.-bl_type_2(i,j))*Puns
          Pnonl=max(P1,P2)
          I1=(1.+1./SQRT2)*(1.+Pnonl)
          E2=(1.+1./SQRT2)*( Pnonl+(1./SQRT2) -                         &
     &                      SQRT(Pnonl*(SQRT2-1.)+0.5))
          GAMMA1(I,J) = I1
          GAMMA2(I,J) = I1 - E2
        END DO
      END DO
!
      L_CORRECT = .TRUE.
! DEPENDS ON: BDY_IMPL3
      CALL BDY_IMPL3 (                                                  &
     & HALO_I,HALO_J,OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,BL_LEVELS,      &
     & L_CORRECT,r_rho_levels,r_theta_levels,                           &
     & Q,QCL,QCF,Q_LATEST,QCL_LATEST,QCF_LATEST,T,T_LATEST,             &
     & DTRDZ_CHARNEY_GRID,DTRDZ_U,DTRDZ_V,                              &
     & RHOKH(1,1,2),RHOKM_U(1,1,2),RHOKM_V(1,1,2),                      &
     & RDZ_CHARNEY_GRID,RDZ_U,RDZ_V,GAMMA1,GAMMA2,GAMMA,                &
     & DU_NT,DV_NT,DQW_NT,DTL_NT,FQW,FTL,TAUX,TAUY,                     &
     & QW,TL,DQW1,DTL1,CT_CTQ,CTCTQ1,DQW,DTL,CQ_CM_U,CQ_CM_V,DU,DV,     &
     & LTIMER                                                           &
     &  )

! DEPENDS ON: SF_IMPL2 
      CALL SF_IMPL2 (                                                   &

! IN values defining field dimensions and subset to be processed :
     & OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,LAND_PTS,                     &

! IN soil/vegetation/land surface data :
     & LAND_INDEX,LAND_MASK,NICE,NTILES,TILE_INDEX,TILE_PTS,SM_LEVELS,  &
     & CANHC_TILE,CANOPY,FLAKE,SMC,TILE_FRAC,WT_EXT_TILE,FLAND,FLANDG,  &

! IN sea/sea-ice data :
     & DI,ICE_FRACT,DI_NCAT,ICE_FRACT_NCAT,U_0,V_0,                     &

! IN everything not covered so far :
     & PSTAR,LW_DOWN,RAD_SICE,SW_TILE,TIMESTEP,T_SOIL,QW,TL,U,V,        &
     & RHOKM_U,RHOKM_V,GAMMA(1),GAMMA1,GAMMA2,ALPHA1,ALPHA1_SICE,       &
     & ASHTF,ASHTF_TILE,DTRDZ_CHARNEY_GRID,DU,DV,FQW_TILE,EPOT_TILE,    &
     & FQW_ICE,FTL_ICE,FRACA,RESFS,RESFT,RHOKH,RHOKH_TILE,RHOKH_SICE,   & 
     & RHOKPM,RHOKPM_POT,RHOKPM_SICE,Z1_TQ,                             &
     & Z0HSSI,Z0MSSI,Z0H_TILE,Z0M_TILE,                                 &
     & CDR10M_U,CDR10M_V,CHR1P5M,CHR1P5M_SICE,CT_CTQ,CTCTQ1,DQW,DTL,    &
     & DQW1,DTL1,DU_STAR,DV_STAR,CQ_CM_U,CQ_CM_V,                       &
     & L_NEG_TSTAR,L_CORRECT,FLANDG_U,FLANDG_V,ANTHROP_HEAT,            &
     & L_SICE_HEATFLUX,                                                 &

! IN STASH flags :-
     & SIMLT,SMLT,SLH,SQ1P5,ST1P5,SU10,SV10,                            &

! INOUT data :
     & TI,TI_GB,TSTAR,TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,TSTAR_SSI,        &
     & TSTAR_TILE,SNOW_TILE,LE_TILE,RADNET_SICE,RADNET_TILE,E_SEA,      &
     & FQW,FTL,FTL_TILE,H_SEA,OLR,TAUX,TAUY,                            &
     & TAUX_LAND,TAUX_LAND_star,TAUX_SSI,TAUX_SSI_star,TAUY_LAND,       &
     & TAUY_LAND_star,TAUY_SSI,TAUY_SSI_star,                           &

! OUT Diagnostic not requiring STASH flags :
     & ECAN,EI_TILE,ESOIL_TILE,                                         &
     & SEA_ICE_HTF,SURF_HT_FLUX,SURF_HT_FLUX_LAND,SURF_HT_FLUX_SICE,    &

! OUT diagnostic requiring STASH flags :
     & SICE_MLT_HTF,SNOMLT_SURF_HTF,LATENT_HEAT,                        &
     & Q1P5M,Q1P5M_TILE,T1P5M,T1P5M_TILE,U10M,V10M,                     &

! OUT data required elsewhere in UM system :
     & ECAN_TILE,EI,ES,EXT,SNOWMELT,MELT_TILE,RHOKH_MIX,ERROR,          &
     & lq_mix_bl, L_FLUX_BC, LTIMER                                     &
     & )

! DEPENDS ON: BDY_IMPL4
      CALL BDY_IMPL4 (                                                  &

! IN values defining field dimensions and subset to be processed :
     & HALO_I,HALO_J,OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,BL_LEVELS,      &

! IN Diagnostic switches
     & L_CORRECT, l_ftl, l_fqw, l_taux, l_tauy,                         &

! IN data :
     & GAMMA1,GAMMA2,R_RHO_LEVELS,RHOKH(1,1,2),RHOKM_U(1,1,2),          &
     & RHOKM_V(1,1,2),RDZ_CHARNEY_GRID,DTRDZ_CHARNEY_GRID,RDZ_U,RDZ_V,  &

! INOUT data :
     & QW,TL,FQW,FTL,TAUX,TAUY,FQW_star,FTL_star,TAUX_star,TAUY_star,   &
     & DU,DV,DU_STAR,DV_STAR,CT_CTQ,CTCTQ1,DQW,DTL,                     &
     & DQW_NT,DTL_NT,CQ_CM_U,CQ_CM_V,                                   &

! OUT data:
     & T_LATEST,Q_LATEST,RHOKH_MIX,                                     &

! LOGICAL LTIMER
     & LTIMER                                                           &
     &  )

      Deallocate (DQW1)
      Deallocate (DTL1)
      Deallocate (CTCTQ1)
!
! Following arrays: stable & non-oscillatory solver fluxes for 1st
!                   stage (predictor)
!
      Deallocate (FQW_star)
      Deallocate (FTL_star)
      Deallocate (TAUX_star)
      Deallocate (TAUY_star)

      ENDIF ! IF L_us_blsol 

      IF ( Substep_Number  ==  Num_Substeps ) Then

      DO J=1,ROWS
        DO I=1,ROW_LENGTH
          FTL_LAND(I,J)=0.0
          FTL_SSI(I,J)=0.0
          E_LAND(I,J)=0.0
          E_SSI(I,J)=0.0
          EI_LAND(I,J)=0.0
          EI_SICE(I,J)=0.0
        ENDDO
      ENDDO

      DO N=1,NTILES
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          FTL_LAND(I,J)=FTL_LAND(I,J) +                                 &
     &      FTL_TILE(L,N)*TILE_FRAC(L,N)
          E_LAND(I,J)=E_LAND(I,J) +                                     &
     &      FQW_TILE(L,N)*TILE_FRAC(L,N)
          EI_LAND(I,J)=EI_LAND(I,J) +                                   &
     &      EI_TILE(L,N)*TILE_FRAC(L,N)
        ENDDO
      ENDDO

      DO J=1,ROWS
        DO I=1,ROW_LENGTH
          IF(FLANDG(I,J) <  1.0)THEN
            FTL_SSI(I,J)=(FTL(I,J,1)-FTL_LAND(I,J)*FLANDG(I,J))         &
     &        /(1.0-FLANDG(I,J))
            E_SSI(I,J)=(FQW(I,J,1)-E_LAND(I,J)*FLANDG(I,J))             &
     &        /(1.0-FLANDG(I,J))
            EI_SICE(I,J)=(EI(I,J)-EI_LAND(I,J)*FLANDG(I,J))             &
     &        /(1.0-FLANDG(I,J))
          ENDIF
        ENDDO
      ENDDO


      Endif

      IF (SQ1P5 .AND. lq_mix_bl) THEN 
!----------------------------------------------------------------------
! Convert 1.5m mixing ratio to specific humidity
!-----------------------------------------------
! Not immediately clear what to do about QCL and QCF at 1.5m so assume 
! the same as level 1.  An alternative would be to assume zero.
!----------------------------------------------------------------------
         DO J=1,ROWS
         DO I=1,ROW_LENGTH
           Q1P5M(I,J)=Q1P5M(I,J)/(1.0+Q1P5M(I,J)+QCL(I,J,1)+QCF(I,J,1))
         ENDDO
         ENDDO

         DO N=1,NTILES
         DO K=1,TILE_PTS(N)
           L = TILE_INDEX(K,N)
           J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
           I = LAND_INDEX(L) - (J-1)*ROW_LENGTH

           Q1P5M_TILE(L,N) = Q1P5M_TILE(L,N)/ (1.0 + Q1P5M_TILE(L,N)    &
     &                                          +QCL(I,J,1)+QCF(I,J,1) )
         ENDDO
         ENDDO
      ENDIF

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('IMPSOLV ',4)
      ENDIF

      RETURN
      END SUBROUTINE IMP_SOLVER
#endif
