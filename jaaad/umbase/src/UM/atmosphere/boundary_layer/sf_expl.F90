#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE SF_EXPL------------------------------------------------
!!!
!!!  Purpose: Calculate explicit surface fluxes of heat, moisture and
!!!           momentum. Also calculates surface exchange coefficients
!!!           required for implicit update of surface fluxes and surface
!!!           information required by the explicit boundary layer
!!!           routine
!!!
!!!
!!!
!!!  Programming standard: Unified Model Documentation Paper No 4,
!!!                        Version ?, dated ?.
!!!
!!!  System component covered: P24.
!!!
!!!  Project task:
!!!
!!!  Documentation: UMDP 24.
!!!
!!!---------------------------------------------------------------------

!    Arguments :-
      SUBROUTINE SF_EXPL (                                              &

! IN values defining field dimensions and subset to be processed :
     & HALO_I,HALO_J,OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,                &
     & LAND_PTS,LAND_PTS_TRIF,NPFT_TRIF,                                &
     & DIM_CS1, DIM_CS2,                                                &

! IN  parameters for iterative SISL scheme
     & NumCycles, CycleNo,                                              &
! IN parameters required from boundary-layer scheme :
     & BQ_1,BT_1,Z1_UV,Z1_TQ,QW_1,TL_1,                                 &

! IN soil/vegetation/land surface data :
     & LAND_INDEX,LAND_MASK,FORMDRAG,FD_stab_dep,OROG_DRAG_PARAM,       &
     & NTILES,SM_LEVELS,                                                &
     & CANOPY,CATCH,CATCH_SNOW,HCON,HO2R2_OROG,                         &

     & FLAND,FLANDG,                                                    &
     & SNOW_TILE,SIL_OROG_LAND,SMVCCL,SMVCST,SMVCWT,                    &
     & STHF,STHU,Z0_TILE,                                               &

! IN sea/sea-ice data :
     & ICE_FRACT,U_0,V_0,U_0_P,V_0_P,Charnock,SeaSalinityFactor,        &

! IN everything not covered so far :
     & PSTAR,LW_DOWN,RAD_SICE,SW_TILE,TIMESTEP,ZH,                      &
     & CO2_MMR,CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,L_CO2_INTERACTIVE,        &
     & L_PHENOL,L_TRIFFID,L_Q10,ASTEPS_SINCE_TRIFFID,CAN_MODEL,         &
     & CS,FRAC,CANHT_FT,PHOTOSYNTH_ACT_RAD,LAI_FT,lq_mix_bl,            &
     & T_SOIL,TI,TSTAR,                                                 &
     & TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,TSTAR_SSI,                       &
     & TSTAR_TILE,Z_LAND,L_CTILE,COR_UST,COR_MO_ITER,I_SCRN_T_DIAG,     &
     & ALBSOIL,COS_ZENITH_ANGLE,                                        &
     & CAN_RAD_MOD, ILAYERS,                                            &
     & U_1,V_1,U_1_P,V_1_P,                                             &
     & L_DUST,ISeaZ0T,Buddy_Sea,ANTHROP_HEAT,                           &

! IN STASH flags :-
     & SFME,SQ1P5,ST1P5,SU10,SV10,SZ0HEFF,                              &

! INOUT data :
     & Z0MSEA,L_spec_z0,Z0M_SCM,Z0H_SCM,GS,                             &
     & G_LEAF_ACC,NPP_FT_ACC,RESP_W_FT_ACC,RESP_S_ACC,                  &

! OUT Diagnostic not requiring STASH flags :
     & CD,CH,RECIP_L_MO_SEA,E_SEA,FQW_1,                                &
     & FTL_1,FTL_TILE,LE_TILE,H_SEA,RADNET_SICE,RADNET_TILE,            &
     & RHOKM_1,RHOKM_U_1,RHOKM_V_1,RIB,RIB_TILE,TAUX_1,TAUY_1,          &
     & TAUX_LAND,TAUX_SSI,TAUY_LAND,TAUY_SSI,                           &

! OUT diagnostic requiring STASH flags :
     & FME,                                                             &

! OUT diagnostics required for soil moisture nudging scheme :
     & WT_EXT,RA,                                                       &

! OUT data required for tracer mixing :
     & RHO_ARESIST,ARESIST,RESIST_B,                                    &
     & RHO_ARESIST_TILE,ARESIST_TILE,RESIST_B_TILE,                     &

!OUT data required for mineral dust scheme
     & R_B_DUST,CD_STD_DUST,U_S_STD_TILE,                               &


! OUT data required for 4D-VAR :
     & RHO_CD_MODV1,                                                    &

! OUT data required elsewhere in UM system :
     & FB_SURF,U_S,T1_SD,Q1_SD,                                         &

! OUT data required elsewhere in boundary layer or surface code
     & ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,FQW_TILE,EPOT_TILE,          &
     & FQW_ICE,FTL_ICE,FRACA,RHOSTAR,RESFS,RESFT,                       &
     & RHOKH,RHOKH_TILE,RHOKH_SICE,RHOKPM,RHOKPM_POT,RHOKPM_SICE,       &
     & H_BLEND_OROG,Z0HSSI,Z0H_TILE,Z0H_EFF,Z0M_GB,Z0MSSI,Z0M_TILE,     &
     & Z0M_EFF,CDR10M_U,CDR10M_V,CHR1P5M,CHR1P5M_SICE,SMC,              &
     & VSHR,VSHR_LAND,VSHR_SSI,                                         &
     & GPP,NPP,RESP_P,G_LEAF,GPP_FT,NPP_FT,                             &
     & RESP_P_FT,RESP_S,RESP_S_TOT,CLAY_LAND,RESP_W_FT,                 &
     & GC,CANHC_TILE,WT_EXT_TILE,FLAKE,                                 &
     & TILE_INDEX,TILE_PTS,TILE_FRAC,FSMC,                              &
     & FLANDG_U,FLANDG_V,                                               &

! LOGICAL LTIMER
     & LTIMER,BL_diag                                                   &
     &, L_ukca,                                                         &
! EAK
     & l_cable,                                                         &
     & surf_down_sw,alb_tile,l_tile_pts,                                &
     & ls_rain,ls_snow,SW_DOWN,                                         &
     & lat,long,day,time_sec,                                           &
     & SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L,                   &
     & SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE,             &
     & TSOIL_TILE,T_SURF_TILE,HCONS,BEXP,                               &
     & SATHH,SATCON,HCAP,                                               &
!     & SOIL_TYPE,VEG_TYPE,                                              &
     & SOIL_TYPE,VEG_TYPE,                                              &

     & ISNOW_FLG3L,total_nsteps,                                        &
     &            FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB,              &
     &            TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB,          &
     &            USTAR_CAB,SURF_HTF_CAB,                               &
     &            U_S_CAB,CH_CAB,CD_CAB,                                &
     & SOIL_LAYER_MOISTURE,                                             &
!sxy
     & SNAGE_TILE,RTSOIL_TILE,                                          &
     & GFLUX_TILE,SGFLUX_TILE,                                          &
     & F_ROOT, sin_theta_latitude,                                      &
     & CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,SOIL_ORDER,                     &
     & NIDEP,NIFIX,PWEA,PDUST,GLAI,PHENPHASE,                           &
      ! end step of experiment, this step, step width, processor num
     & endstep, timestep_number, mype                                   &
     & )
 
      USE rad_switches_mod, ONLY: LRAD_EMIS_LAND_GEN,RAD_EMIS_LAND_GEN

      Use bl_diags_mod, Only:                                           &
          strnewbldiag

      IMPLICIT NONE


#include "c_dust_ndiv.h"


!  Inputs :-

! (a) Defining horizontal grid and subset thereof to be processed.
!    Checked for consistency.

      INTEGER                                                           &
     & ROW_LENGTH                                                       &
                  ! Local number of points on a row
     &,ROWS                                                             &
                  ! Local number of rows in a theta field
     &,N_ROWS                                                           &
                  ! Local number of rows in a v field
     &,OFF_X                                                            &
                  ! Size of small halo in i.
     &,OFF_Y                                                            &
                  ! Size of small halo in j.
     &,HALO_I                                                           &
                  ! Size of halo in i direction.
     &,HALO_J                                                           &
                  ! Size of halo in j direction.
     &,LAND_PTS                                                         &
                  ! IN No of land points being processed.
     &,LAND_PTS_TRIF                                                    &
!                 ! IN For dimensioning land fields
     &,NPFT_TRIF  ! IN For dimensioning PFT fields available only
!                 !    with TRIFFID. Set to NPFT when TRIFFID on,
!                 !    set to 1 when TRIFFID off.
      INTEGER                                                           &
     & CAN_RAD_MOD                                                      &
!                                  !Switch for canopy radiation model
     & ,ILAYERS                                                         &
!                                  !No of layers in canopy radiation model
     &,NumCycles                                                        &
                  ! Number of cycles (iterations) for iterative SISL.
     &,CycleNo                                                          &
                  ! Iteration no
     &,DIM_CS1, DIM_CS2      ! soil carbon dimensions

      INTEGER, Intent(IN) :: I_SCRN_T_DIAG
!       Method of diagnosing the screen temperature

! Defining vertical grid of model atmosphere.
      REAL                                                              &
     & BQ_1(ROW_LENGTH,ROWS)                                            &
                                   ! IN A buoyancy parameter
!                                  !    (beta q tilde).
     &,BT_1(ROW_LENGTH,ROWS)                                            &
                                   ! IN A buoyancy parameter
!                                  !    (beta T tilde).
     &,Z1_UV(ROW_LENGTH,ROWS)                                           &
                                   ! IN Height of lowest uv level (m).
     &,Z1_TQ(ROW_LENGTH,ROWS)                                           &
                                   ! IN Height of lowest tq level (m).
!                                  !    Note, if the grid used is
!                                  !    staggered in the vertical,
!                                  !    Z1_UV and Z1_TQ can be
!                                  !    different.
     &,QW_1(ROW_LENGTH,ROWS)                                            &
                                   ! IN Total water content
     &,TL_1(ROW_LENGTH,ROWS)       ! IN Ice/liquid water temperature


! (c) Soil/vegetation/land surface parameters (mostly constant).

      LOGICAL                                                           &
     & LAND_MASK(ROW_LENGTH,ROWS)                                       &
                                   ! IN T if land, F elsewhere.
     &,L_CO2_INTERACTIVE                                                &
                                   ! IN Switch for 3D CO2 field
     &,L_PHENOL                                                         &
                                   ! IN Indicates whether phenology
!                                  !    in use
     &,L_TRIFFID                                                        &
                                   ! IN Indicates whether TRIFFID
!                                  !    in use.
     &,L_CTILE                                                          &
                                   ! IN True if coastal tiling
     &,L_spec_z0                                                        &
                                   ! IN T if using prescribed
!                                  !    sea surface roughness lengths
     &,L_Q10                       ! IN True if using Q10 for soil resp

      INTEGER                                                           &
     & LAND_INDEX(LAND_PTS)        ! IN LAND_INDEX(I)=J => the Jth
!                                  !    point in ROW_LENGTH,ROWS is the
!                                  !    land point.

      INTEGER                                                           &
     & SM_LEVELS                                                        &
                                   ! IN No. of soil moisture levels
     &,NTILES                                                           &
                                   ! IN No. of land-surface tiles
     &,CO2_DIM_LEN                                                      &
                                   ! IN Length of a CO2 field row.
     &,CO2_DIM_ROW                                                      &
                                   ! IN Number of CO2 field rows.
     &,ASTEPS_SINCE_TRIFFID                                             &
                                   ! IN Number of atmospheric
!                                  !    timesteps since last call
!                                  !    to TRIFFID.
     &,CAN_MODEL                                                        &
                                   ! IN Swith for thermal vegetation
!                                  !    canopy
     &,FORMDRAG                                                         &
                                   ! IN Switch for orographic drag
     &,FD_stab_dep                 ! IN Switch to implement stability
!                                  !    dependence of orog form drag

      REAL                                                              &
     & CANOPY(LAND_PTS,NTILES)                                          &
                                   ! IN Surface/canopy water for
!                                  !    snow-free land tiles (kg/m2)
     &,CATCH(LAND_PTS,NTILES)                                           &
                                   ! IN Surface/canopy water capacity
!                                  !    of snow-free land tiles (kg/m2).
     &,CATCH_SNOW(LAND_PTS)                                             &
                                   ! IN Snow interception capacity of
!                                  !    NLT tile (kg/m2).
     &,HCON(LAND_PTS)                                                   &
                                   ! IN Soil thermal conductivity
!                                  !    (W/m/K).
     &,SNOW_TILE(LAND_PTS,NTILES)                                       &
                                   ! IN Lying snow on tiles (kg/m2)
     &,SMVCCL(LAND_PTS)                                                 &
                                   ! IN Critical volumetric SMC
!                                  !    (cubic m per cubic m of soil).
     &,SMVCST(LAND_PTS)                                                 &
                                   ! IN Volumetric saturation point
!                                  !    (m3/m3 of soil).
     &,SMVCWT(LAND_PTS)                                                 &
                                   ! IN Volumetric wilting point
!                                  !    (cubic m per cubic m of soil).
     &,STHF(LAND_PTS,SM_LEVELS)                                         &
                                   ! IN Frozen soil moisture content of
!                                  !    each layer as a fraction of
!                                  !    saturation.
     &,STHU(LAND_PTS,SM_LEVELS)                                         &
                                   ! IN Unfrozen soil moisture content
!                                  !    of each layer as a fraction of
!                                  !    saturation.
     &,Z0_TILE(LAND_PTS,NTILES)                                         &
                                   ! IN Tile roughness lengths (m).
     &,SIL_OROG_LAND(LAND_PTS)                                          &
                                   ! IN Silhouette area of unresolved
!                                  !    orography per unit horizontal
!                                  !    area on land points only.
     &,HO2R2_OROG(LAND_PTS)                                             &
                                   ! IN Standard Deviation of orography.
!                                  !    equivilent to peak to trough
!                                  !    height of unresolved orography
     &,OROG_DRAG_PARAM                                                  &
!                                  ! IN Orographic form drag coefficient
     &,FLAND(LAND_PTS)                                                  &
                                   ! IN Land fraction on land tiles.
     &,FLANDG(ROW_LENGTH,ROWS)
!                                  ! IN Land fraction on all tiles.
!                                  !    divided by 2SQRT(2) on land
!                                  !    points only (m)

! (d) Sea/sea-ice data.

      REAL                                                              &
     & ICE_FRACT(ROW_LENGTH,ROWS)                                       &
                                   ! IN Fraction of gridbox covered by
!                                  !    sea-ice (decimal fraction).
     &,U_0(ROW_LENGTH,ROWS)                                             &
                                   ! IN W'ly component of surface
!                                  !    current (m/s).
     &,V_0(ROW_LENGTH,N_ROWS)                                           &
                                   ! IN S'ly component of surface
!                                  !    current (m/s).
     &,U_0_P(ROW_LENGTH,ROWS)                                           &
                                   ! IN W'ly component of surface
!                                       current (m/s). P grid
     &,V_0_P(ROW_LENGTH,ROWS)                                           &
                                   ! IN S'ly component of surface
!                                       current (m/s). P grid
     & ,Charnock   ! Charnock parameter for sea surface


! (f) Atmospheric + any other data not covered so far, incl control.

#include "nstypes.h"

      REAL                                                              &
     & PSTAR(ROW_LENGTH,ROWS)                                           &
                                   ! IN Surface pressure (Pascals).
     &,LW_DOWN(ROW_LENGTH,ROWS)                                         &
                                   ! IN Surface downward LW radiation
!                                  !    (W/m2).
     &,SW_DOWN(ROW_LENGTH,ROWS)                                         & 
                                   ! IN Surface downward SW radiation (W/m2).
     &,RAD_SICE(ROW_LENGTH,ROWS)                                        &
                                   ! IN Surface net shortwave and
!                                  !    downward LWradiation for
!                                  !    sea-ice (W/sq m).
     &,SW_TILE(LAND_PTS,NTILES)                                         &
                                   ! IN Surface net SW radiation on
!                                  !    land tiles (W/m2).
     &,TIMESTEP                                                         &
                                   ! IN Timestep (seconds).
     &,ZH(ROW_LENGTH,ROWS)                                              &
                                   ! IN Height above surface of top of
!                                  !    boundary layer (metres).
     &,CO2_MMR                                                          &
                                   ! IN CO2 Mass Mixing Ratio
     &,CO2_3D(CO2_DIM_LEN,CO2_DIM_ROW)                                  &
!                                  ! IN 3D CO2 field if required.
     &,CS(LAND_PTS,DIM_CS1)                                             &
                              ! IN Soil carbon (kg C/m2).
     &,FRAC(LAND_PTS,NTYPE)                                             &
                                   ! IN Fractions of surface types.
     &,CANHT_FT(LAND_PTS,NPFT)                                          &
                                   ! IN Canopy height (m)
     &,PHOTOSYNTH_ACT_RAD(ROW_LENGTH,ROWS)                              &
!                                  ! IN Net downward shortwave radiation
!                                  !    in band 1 (w/m2).
     &,LAI_FT(LAND_PTS,NPFT)                                            &
                                   ! IN Leaf area index
     &,T_SOIL(LAND_PTS,SM_LEVELS)                                       &
                                   ! IN Soil temperatures (K).
     &,TI(ROW_LENGTH,ROWS)                                              &
                                   ! IN Sea-ice surface layer
     &,TSTAR_LAND(ROW_LENGTH,ROWS)                                      &
                                   ! IN Land mean surface temperature (K
     &,TSTAR_SEA(ROW_LENGTH,ROWS)                                       &
                                   ! IN Open sea surface temperature (K)
     &,TSTAR_SICE(ROW_LENGTH,ROWS)                                      &
                                   ! IN Sea-ice surface temperature (K).
     &,TSTAR_SSI(ROW_LENGTH,ROWS)                                       &
                                   ! IN mean sea surface temperature (K)
!                                  !    temperature (K).
     &,TSTAR(ROW_LENGTH,ROWS)                                           &
                                   ! IN GBM surface temperature (K).
     &,TSTAR_TILE(LAND_PTS,NTILES)                                      &
                                   ! IN Surface tile temperatures
     &,Z_LAND(ROW_LENGTH,ROWS)                                          &
                                   ! IN Land height (m).
     &,ALBSOIL(LAND_PTS)                                                &
!                                  ! Soil albedo.
     &, COS_ZENITH_ANGLE(row_length, rows)                              &
!                                  ! Cosine of the zenith angle
     &, sin_theta_latitude(row_length, rows)                            &
!                                  ! 
     &,U_1(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y)                 &
!                                  ! IN W'ly wind component (m/s)
     &,V_1(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:N_ROWS+OFF_Y)               &
!                                  ! IN S'ly wind component (m/s)
     &,U_1_P(ROW_LENGTH,ROWS)                                           &
                                   ! IN U_1 on P-grid.
     &,V_1_P(ROW_LENGTH,ROWS)                                           &
                                   ! IN V_1 on P-grid.
     &,ANTHROP_HEAT(NTILES)
!                                  ! IN Additional heat source on tiles 
!                                  !    used for anthropgenic urban 
!                                  !    heat source (W/m2)
!
!------------------------------------------------------
!     EAK
      Real                                                              &
     &  alb_tile(land_pts,ntiles,4)                                     &
     &, SNOW_DEPTH3L(LAND_PTS,NTILES,3)                          &
     &, SNOW_MASS3L(LAND_PTS,NTILES,3)                           &
     &, SNOW_COND(LAND_PTS,NTILES,3)                             &
     &, SNOW_TMP3L(LAND_PTS,NTILES,3)                            &
     &, SNOW_RHO3L(LAND_PTS,NTILES,3)                            &
     &, SNOW_RHO1L(LAND_PTS,NTILES)                              &
     &, SNAGE_TILE(LAND_PTS,NTILES)                              &
     &, SMCL_TILE(LAND_PTS,NTILES,sm_levels)                     &
     &, STHU_TILE(LAND_PTS,NTILES,sm_levels)                     &
     &, STHF_TILE(LAND_PTS,NTILES,sm_levels)                     &
     &, TSOIL_TILE(LAND_PTS,NTILES,sm_levels)                    &
     &, T_SURF_TILE(LAND_PTS,NTILES)                             &
     &, RTSOIL_TILE(LAND_PTS,NTILES)                             &
     &, GFLUX_TILE(LAND_PTS,NTILES)                              &
     &, SGFLUX_TILE(LAND_PTS,NTILES)                             &
     &, BEXP(LAND_PTS)                                           &
     &, SATCON(LAND_PTS)                                         &
     &, SATHH(LAND_PTS)                                          &
     &, HCAP(LAND_PTS)                                           &
     &, HCONS(LAND_PTS)                                          &
     &, SOIL_LAYER_MOISTURE(LAND_PTS,SM_LEVELS)                  &
                                               !IN soil moisture
!                                              ! per layer (kg m-2)
     &, CPOOL_TILE(LAND_PTS,NTILES,10)                           &
     &, NPOOL_TILE(LAND_PTS,NTILES,10)                           &
     &, PPOOL_TILE(LAND_PTS,NTILES,12)                           &
     &, SOIL_ORDER(LAND_PTS)                                     &
     &, NIDEP(land_pts)                                          &
     &, NIFIX(land_pts)                                          &
     &, PWEA(land_pts)                                           &
     &, PDUST(land_pts)                                          &
     &, GLAI(land_pts,NTILES)                                    &
     &, PHENPHASE(land_pts,NTILES) &

     &, surf_down_sw(row_length,rows,4)                                 &
     &, ls_rain(row_length, rows)                                       &
     &, ls_snow(row_length, rows)                                       &
     &, conv_rain(row_length, rows)                                     &
     &, conv_snow(row_length, rows)                                     &
     &, LAT(ROW_LENGTH,ROWS)                                            &
     &, LONG(ROW_LENGTH,ROWS)                                           &
     &, F_ROOT(sm_levels)                                               &
     &, time_sec                                                        

      Real                                                           &
     & FTL_TILE_CAB(LAND_PTS,NTILES)                                 &
     &,FTL_CAB(LAND_PTS)                                             &
     &,LE_TILE_CAB(LAND_PTS,NTILES)                                  &
     &,LE_CAB(LAND_PTS)                                              &
     &,TSTAR_TILE_CAB(LAND_PTS,NTILES)                               &
     &,TSTAR_CAB(LAND_PTS)                                           &
     &,SMCL_CAB(LAND_PTS,SM_LEVELS)                          &
     &,TSOIL_CAB(LAND_PTS,SM_LEVELS)                         &
     &,USTAR_CAB(LAND_PTS)                                           &
     &,SURF_HTF_CAB(LAND_PTS)                                        &
!sxy
     &,CH_CAB(LAND_PTS)                                              &
     &,CD_CAB(LAND_PTS)                                              &
     &,U_S_CAB(LAND_PTS)

      Logical                                                           &
     & l_cable
      LOGICAL, DIMENSION(LAND_PTS,NTILES) :: L_TILE_PTS

      Integer                                                           &
     &  day                                                             &
     &, total_nsteps                                                    &
                                ! Total number of steps in run
     &, SOIL_TYPE(row_length,rows)                                      &
     &, VEG_TYPE(row_length,rows)                                       &
     &, ISNOW_FLG3L(LAND_PTS,NTILES)
!
!--------------------------------------------------------------------
      LOGICAL                                                           &
     & LTIMER                                                           &
                                   ! IN Logical switch for TIMER diags
     &,L_DUST                      !IN switch for mineral dust
!
      REAL, Intent(IN) :: SeaSalinityFactor
!                                  ! Factor allowing for the effect
!                                  ! of the salinity of sea water
!                                  ! on the evaporative flux.
      INTEGER, Intent(IN) :: ISeaZ0T
!                                  ! Switch for the definition of
!                                  ! the thermal roughness
!                                  ! length over the sea.
      INTEGER, Intent(IN) :: COR_UST 
!                                  ! Switch to use correct 
!                                  ! friction velocity
      INTEGER, Intent(IN) :: COR_MO_ITER
!                                  ! Switch for MO iteration correction
      INTEGER, Intent(IN) :: Buddy_sea
!                                  ! Switch to use the wind speed from 
!                                  ! adjacent sea points for the 
!                                  ! sea part of coastal grid points
!
      LOGICAL                                                           &
     & lq_mix_bl
      Logical :: L_ukca   ! switch for UKCA scheme
!  STASH flags :-

      LOGICAL                                                           &
     & SFME                                                             &
               ! IN Flag for FME (q.v.).
     &,SZ0HEFF                                                          &
               ! IN Flag for Z0H_EFF
     &,SQ1P5                                                            &
               ! IN Flag for Q1P5M (q.v.)
     &,ST1P5                                                            &
               ! IN Flag for T1P5M (q.v.)
     &,SU10                                                             &
               ! IN Flag for U10M (q.v.)
     &,SV10    ! IN Flag for V10M (q.v.)

!  In/outs :-

      REAL                                                              &
     & Z0MSEA(ROW_LENGTH,ROWS)                                          &
                                   ! INOUT Sea-surface roughness
!                                  !       length for momentum (m).
     &,Z0M_SCM(ROW_LENGTH,ROWS)                                         &
                                   ! IN Fixed Sea-surface roughness
!                                  !    length for momentum (m).(SCM)
     &,Z0H_SCM(ROW_LENGTH,ROWS)                                         &
                                   ! IN Fixed Sea-surface roughness
!                                  !    length for heat (m). (SCM)
     &,GS(LAND_PTS)                                                     &
                                   ! INOUT "Stomatal" conductance to
!                                  !        evaporation (m/s).
     &,G_LEAF_ACC(LAND_PTS,NPFT)                                        &
                                   ! INOUT Accumulated G_LEAF
! Lestevens 17may13 - change npft to ntiles
     &,NPP_FT_ACC(LAND_PTS_TRIF,NTILES)                                 &
!    &,NPP_FT_ACC(LAND_PTS_TRIF,NPFT_TRIF)                              &
!                                  ! INOUT Accumulated NPP_FT
     &,RESP_W_FT_ACC(LAND_PTS_TRIF,NTILES)                              &
!    &,RESP_W_FT_ACC(LAND_PTS_TRIF,NPFT_TRIF)                           &
!                                  ! INOUT Accum RESP_W_FT
     &,RESP_S_ACC(LAND_PTS_TRIF,DIM_CS1) ! INOUT Accumulated RESP_S

!  Outputs :-
!-1 Diagnostic (or effectively so - includes coupled model requisites):-

      INTEGER                                                           &
     & TILE_INDEX(LAND_PTS,NTYPE)                                       &
                                   ! OUT Index of tile points
     &,TILE_PTS(NTYPE)             ! OUT Number of tile points


!  (a) Calculated anyway (use STASH space from higher level) :-
!
      REAL                                                              &
     & CD(ROW_LENGTH,ROWS)                                              &
                                   ! OUT Turbulent surface exchange
!                                  !     (bulk transfer) coefficient for
!                                  !     momentum.
     &,CH(ROW_LENGTH,ROWS)                                              &
                                   ! OUT Turbulent surface exchange
!                                  !     (bulk transfer) coefficient for
!                                  !     heat and/or moisture.
     &,RECIP_L_MO_SEA(row_length,rows)                                  &
!                                  ! OUT Reciprocal of the surface
!                                  !     Obukhov  length at sea
!                                  !     points. (m-1).
     &,E_SEA(ROW_LENGTH,ROWS)                                           &
                                   ! OUT Evaporation from sea times
!                                  !     leads fraction. Zero over land.
!                                  !     (kg per square metre per sec).
     &,FQW_1(ROW_LENGTH,ROWS)                                           &
                                   ! OUT Moisture flux between layers
!                                  !     (kg per square metre per sec).
!                                  !     FQW(,1) is total water flux
!                                  !     from surface, 'E'.
     &,FTL_1(ROW_LENGTH,ROWS)                                           &
                                   ! OUT FTL(,K) contains net turbulent
!                                  !     sensible heat flux into layer K
!                                  !     from below; so FTL(,1) is the
!                                  !     surface sensible heat, H.(W/m2)
     &,FTL_TILE(LAND_PTS,NTILES)                                        &
                                   ! OUT Surface FTL for land tiles
     &,LE_TILE(LAND_PTS,NTILES)                                         &
                                   ! OUT Surface latent heat flux for
!                                  !     land tiles
     &,H_SEA(ROW_LENGTH,ROWS)                                           &
                                   ! OUT Surface sensible heat flux over
!                                  !     sea times leads fraction (W/m2)
     &,RADNET_SICE(ROW_LENGTH,ROWS)                                     &
                                   ! OUT Surface net radiation on
!                                  !     sea-ice (W/m2)
     &,RADNET_TILE(LAND_PTS,NTILES)                                     &
                                   ! OUT Surface net radiation on
     &,RHOKM_LAND(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y)          &
!                                  ! OUT Exchange coefficients for
!                                  !     momentum on P-grid
     &,RHOKM_SSI(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y)           &
!                                  ! OUT Exchange coefficients for
!                                  !     momentum on P-grid
     &,FLANDG_U(ROW_LENGTH,ROWS)                                        &
                                   ! OUT Land frac (on U-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
!                                  !     land tiles (W/m2)
     &,RHOKM_1(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y)             &
!                                  ! OUT Exchange coefficients for
!                                  !     momentum on P-grid
     &,RHOKM_U_LAND(ROW_LENGTH,ROWS)                                    &
                                      ! OUT Exchange coefficients for
!                                  !     momentum (on U-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
     &,RHOKM_U_SSI(ROW_LENGTH,ROWS)                                     &
                                     ! OUT Exchange coefficients for
!                                  !     momentum (on U-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
     &,FLANDG_V(ROW_LENGTH,N_ROWS)                                      &
                                   ! OUT Land frac (on V-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
     &,RHOKM_U_1(ROW_LENGTH,ROWS)                                       &
                                   ! OUT Exchange coefficients for
!                                  !     momentum (on U-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
     &,RHOKM_V_LAND(ROW_LENGTH,N_ROWS)                                  &
!                                  ! OUT Exchange coefficients for
!                                  !     momentum (on V-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
     &,RHOKM_V_SSI(ROW_LENGTH,N_ROWS)                                   &
!                                  ! OUT Exchange coefficients for
!                                  !     momentum (on V-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
     &,RHOKM_V_1(ROW_LENGTH,N_ROWS)                                     &
                                   ! OUT Exchange coefficients for
!                                  !     momentum (on V-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
     &,RIB(ROW_LENGTH,ROWS)                                             &
                                   ! OUT Mean bulk Richardson number for
!                                  !     lowest layer.
     &,RIB_TILE(LAND_PTS,NTILES)                                        &
                                   ! OUT RIB for land tiles.
     &,TAUX_1(ROW_LENGTH,ROWS)                                          &
                                   ! OUT W'ly component of surface wind
     &,TAUX_LAND(ROW_LENGTH,ROWS)                                       &
                                   ! OUT W'ly component of sfc wind
!                                  !     stress (N/sq m). (On U-grid
!                                  !     with first and last rows
!                                  !     undefined or, at present,
!                                  !     set to missing data
     &,TAUX_SSI(ROW_LENGTH,ROWS)                                        &
                                   ! OUT W'ly component of sfc wind
!                                  !     stress (N/sq m). (On U-grid
!                                  !     with first and last rows
!                                  !     undefined or, at present,
!                                  !     set to missing data
!                                  !     stress (N/sq m). (On U-grid
!                                  !     with first and last rows
!                                  !     undefined or, at present,
     &,TAUY_LAND(ROW_LENGTH,N_ROWS)                                     &
                                   ! OUT S'ly component of sfc wind
!                                  !     stress (N/sq m).  On V-grid;
!                                  !     comments as per TAUX.
     &,TAUY_SSI(ROW_LENGTH,N_ROWS)                                      &
                                   ! OUT S'ly component of sfc wind
!                                  !     stress (N/sq m).  On V-grid;
!                                  !     comments as per TAUX.
!                                  !     set to missing data
     &,TAUY_1(ROW_LENGTH,N_ROWS)                                        &
                                   ! OUT S'ly component of surface wind
!                                  !     stress (N/sq m).  On V-grid;
!                                  !     comments as per TAUX.
     &,RHO_CD_MODV1(ROW_LENGTH,ROWS)                                    &
!                                  ! OUT Surface air density * drag coef
!                                  !     *mod(v1 - v0) before interp
     &,RHO_ARESIST(ROW_LENGTH,ROWS)                                     &
                                   ! OUT RHOSTAR*CD_STD*VSHR for Sulphur
!                                  !     cycle
     &,ARESIST(ROW_LENGTH,ROWS)                                         &
                                   ! OUT 1/(CD_STD*VSHR) for Sulphur
!                                  !     cycle
     &,RESIST_B(ROW_LENGTH,ROWS)                                        &
                                   ! OUT (1/CH-1/(CD_STD)/VSHR for
!                                  !     Sulphur cycle
     &,RHO_ARESIST_TILE(LAND_PTS,NTILES)                                &
!                                  ! OUT RHOSTAR*CD_STD*VSHR on land
!                                  !     tiles
     &,ARESIST_TILE(LAND_PTS,NTILES)                                    &
!                                  ! OUT 1/(CD_STD*VSHR) on land tiles
     &,RESIST_B_TILE(LAND_PTS,NTILES)                                   &
!                                  ! OUT (1/CH-1/CD_STD)/VSHR on land
!                                  !     tiles
     &, R_B_DUST(ROW_LENGTH,ROWS,NDIV)                                  &
                                      !OUT surf layer res for dust
     &, CD_STD_DUST(ROW_LENGTH,ROWS)                                    &
                                    ! OUT Bulk transfer coef. for
!                             ! momentum, excluding orographic effects
     &, U_S_STD_TILE(LAND_PTS,NTILES)                                   &
                                      ! OUT Surface friction velocity
!                             !     (standard value)
     &,WT_EXT(LAND_PTS,SM_LEVELS)                                       &
                                 !OUT cumulative fraction of transp'n
     &,RA(LAND_PTS)               ! OUT Aerodynamic resistance (s/m).


!  (b) Not passed between lower-level routines (not in workspace at this
!      level) :-

      REAL                                                              &
     & FME(ROW_LENGTH,ROWS)        ! OUT Wind mixing "power" (W/m2).

!-2 Genuinely output, needed by other atmospheric routines :-

      REAL                                                              &
     & FB_SURF(ROW_LENGTH,ROWS)                                         &
                                   ! OUT Surface flux buoyancy over
!                                  !     density (m^2/s^3)
     &,U_S(ROW_LENGTH,ROWS)                                             &
                                   ! OUT Surface friction velocity (m/s)
     &,T1_SD(ROW_LENGTH,ROWS)                                           &
                                   ! OUT Standard deviation of turbulent
!                                  !     fluctuations of layer 1 temp;
!                                  !     used in initiating convection.
     &,Q1_SD(ROW_LENGTH,ROWS)      ! OUT Standard deviation of turbulent
!                                  !     flucs of layer 1 humidity;
!                                  !     used in initiating convection.

      REAL                                                              &
     & ALPHA1(LAND_PTS,NTILES)                                          &
                                   ! OUT Mean gradient of saturated
!                                  !     specific humidity with respect
!                                  !     to temperature between the
!                                  !     bottom model layer and tile
!                                  !     surfaces
     &,ALPHA1_SICE(ROW_LENGTH,ROWS)                                     &
                                   ! OUT ALPHA1 for sea-ice.
     &,ASHTF(ROW_LENGTH,ROWS)                                           &
                                   ! OUT Coefficient to calculate
!                                  !     surface heat flux into soil or
!                                  !     sea-ice.
     &,ASHTF_TILE(LAND_PTS,NTILES)                                      &
                                   ! OUT Coefficient to calculate
!                                  !     surface heat flux into land
!                                  !     tiles.
     &,FQW_TILE(LAND_PTS,NTILES)                                        &
                                   ! OUT Surface FQW for land tiles
     &,EPOT_TILE(LAND_PTS,NTILES)                                       &
                                   ! OUT Local EPOT for land tiles.
!*I SFEXPL8A.407
     &,RHOKPM_POT(LAND_PTS,NTILES)                                      &
!                                  ! OUT Potential evaporation
!                                  !     exchange coeff.
     &,FQW_ICE(ROW_LENGTH,ROWS)                                         &
                                   ! OUT Surface FQW for sea-ice
     &,FTL_ICE(ROW_LENGTH,ROWS)                                         &
                                   ! OUT Surface FTL for sea-ice
     &,FRACA(LAND_PTS,NTILES)                                           &
                                   ! OUT Fraction of surface moisture
!                                  !     flux with only aerodynamic
!                                  !     resistance for snow-free land
!                                  !     tiles.
     &,RHOSTAR(ROW_LENGTH,ROWS)                                         &
                                   ! OUT Surface air density
     &,RESFS(LAND_PTS,NTILES)                                           &
                                   ! OUT Combined soil, stomatal
!                                  !     and aerodynamic resistance
!                                  !     factor for fraction (1-FRACA)
!                                  !     of snow-free land tiles.
     &,RESFT(LAND_PTS,NTILES)                                           &
                                   ! OUT Total resistance factor.
!                                  !     FRACA+(1-FRACA)*RESFS for
!                                  !     snow-free land, 1 for snow.
     &,RHOKH(ROW_LENGTH,ROWS)                                           &
                                   ! OUT Grid-box surface exchange
!                                  !     coefficients
     &,RHOKH_TILE(LAND_PTS,NTILES)                                      &
                                   ! OUT Surface exchange coefficients
!                                  !     for land tiles
     &,RHOKH_SICE(ROW_LENGTH,ROWS)                                      &
                                   ! OUT Surface exchange coefficients
!                                  !     for sea and sea-ice
     &,RHOKPM(LAND_PTS,NTILES)                                          &
                                   ! OUT Land surface exchange coeff.
     &,RHOKPM_SICE(ROW_LENGTH,ROWS)                                     &
                                   ! OUT Sea-ice surface exchange coeff.
     &,H_BLEND_OROG(ROW_LENGTH,ROWS)                                    &
!                                  ! OUT Blending height used as part of
!                                  !     effective roughness scheme
     &,Z0HSSI(ROW_LENGTH,ROWS)                                          &
                                   ! OUT Roughness length for heat and
!                                  !     moisture over sea (m).
     &,Z0MSSI(ROW_LENGTH,ROWS)                                          &
                                   ! OUT Roughness length for momentum
!                                  !     over sea (m).
     &,Z0H_TILE(LAND_PTS,NTILES)                                        &
                                   ! OUT Tile roughness lengths for heat
!                                  !     and moisture (m).
     &,Z0H_EFF(ROW_LENGTH,ROWS)                                         &
                                   ! OUT Effective grid-box roughness
!                                  !     length for heat, moisture (m)
     &,Z0M_GB(ROW_LENGTH,ROWS)                                          &
                                   ! OUT Gridbox mean roughness length
!                                  !     for momentum (m).
     &,Z0M_TILE(LAND_PTS,NTILES)                                        &
                                   ! OUT Tile roughness lengths for
!                                  !     momentum.
     &,Z0M_EFF(ROW_LENGTH,ROWS)                                         &
                                   ! OUT Effective grid-box roughness
!                                  !     length for momentum
     &,CDR10M_U(ROW_LENGTH,ROWS)                                        &
                                   ! OUT Ratio of CD's reqd for
!                                  !     calculation of 10 m wind. On
!                                  !     U-grid; comments as per RHOKM.
     &,CDR10M_V(ROW_LENGTH,N_ROWS)                                      &
                                   ! OUT Ratio of CD's reqd for
!                                  !     calculation of 10 m wind. On
!                                  !     V-grid; comments as per RHOKM.
     &,CHR1P5M(LAND_PTS,NTILES)                                         &
                                   ! OUT Ratio of coefffs for
!                                  !     calculation of 1.5m temp for
!                                  !     land tiles.
     &,CHR1P5M_SICE(ROW_LENGTH,ROWS)                                    &
!                                  ! OUT CHR1P5M for sea and sea-ice
!                                  !     (leads ignored).
     &,SMC(LAND_PTS)                                                    &
                                   ! OUT Available moisture in the
!                                  !     soil profile (mm).
     &,VSHR(ROW_LENGTH,ROWS)                                            &
                                   ! OUT Magnitude of surface-to-lowest
!                                  !     atm level wind shear (m per s).
     &,VSHR_LAND(ROW_LENGTH,ROWS)                                       &
                                   ! OUT Magnitude of surface-to-lowest
!                                  !     atm level wind shear (m per s).
     &,VSHR_SSI(ROW_LENGTH,ROWS)                                        &
                                   ! OUT Magnitude of surface-to-lowest
!                                  !     atm level wind shear (m per s).
     &,GPP(LAND_PTS)                                                    &
                                   ! OUT Gross primary productivity
!                                  !     (kg C/m2/s).
     &,NPP(LAND_PTS)                                                    &
                                   ! OUT Net primary productivity
!                                  !     (kg C/m2/s).
     &,RESP_P(LAND_PTS)                                                 &
                                   ! OUT Plant respiration (kg C/m2/s).
     !kdcorbin, 11/10 - changed from NPFT
     &,G_LEAF(LAND_PTS,NTILES)                                          &
                                   ! OUT Leaf turnover rate (/360days).
     !kdcorbin, 11/10 - changed from NPFT
     &,GPP_FT(LAND_PTS,NTILES)                                          &
                                   ! OUT Gross primary productivity
                                   !     on PFTs (kg C/m2/s).
     !kdcorbin, 11/10 - changed from NPFT
     &,NPP_FT(LAND_PTS,NTILES)                                          &
                                   ! OUT Net primary productivity
                                   !     (kg C/m2/s).
     !kdcorbin, 11/10 - changed from NPFT
     &,RESP_P_FT(LAND_PTS,NTILES)                                       &
                                   ! OUT Plant respiration on PFTs
!                                  !     (kg C/m2/s).
     &,RESP_S(LAND_PTS,DIM_CS1)                                         &
                                ! OUT Soil respiration (kg C/m2/s).
     &,RESP_S_TOT(DIM_CS2)                                              &
                                  ! OUT Total soil respiration
                                  ! (kg C/m2/s).
     &,CLAY_LAND(DIM_CS2)                                               &
                                  ! IN clay fraction on land points
     &,RESP_FRAC(DIM_CS2)                                               &
                                  ! respired fraction of RESP_S
     &,WORK_CLAY(DIM_CS2)                                               &
                                  ! working variable
     &,RESP_W_FT(LAND_PTS,NPFT)                                         &
                                   ! OUT Wood maintenance respiration
!                                  !     (kg C/m2/s).
     &,GC(LAND_PTS,NTILES)                                              &
                                   ! OUT "Stomatal" conductance to
!                                  !      evaporation for land tiles
!                                  !      (m/s).
     &,CANHC_TILE(LAND_PTS,NTILES)                                      &
                                   ! IN Areal heat capacity of canopy
!                                  !    for land tiles (J/K/m2).
     &,WT_EXT_TILE(LAND_PTS,SM_LEVELS,NTILES)                           &
!                                  ! IN Fraction of evapotranspiration
!                                  !    which is extracted from each
!                                  !    soil layer by each tile.
     &,FLAKE(LAND_PTS,NTILES)                                           &
                                   ! IN Lake fraction.
     &,TILE_FRAC(LAND_PTS,NTILES)                                       &
                                   ! OUT Tile fractions including
!                                  !     snow cover in the ice tile.
     &,FSMC(LAND_PTS,NPFT)         ! OUT Moisture availability factor.


!     Declaration of new BL diagnostics.
      Type (Strnewbldiag) :: BL_diag

!---------------------------------------------------------------------
!  External routines called :-

      EXTERNAL TILEPTS,PHYSIOL,HEAT_CON,SF_EXCH,                        &
     & Swap_Bounds,FILL_EXTERNAL_HALOS,P_TO_U,P_TO_V
      EXTERNAL TIMER
#include "c_0_dg_c.h"

!-----------------------------------------------------------------------
!   Symbolic constants (parameters) reqd in top-level routine :-

#include "c_r_cp.h"
#include "c_g.h"
#include "c_lheat.h"
#include "csigma.h"
#include "soil_thick.h"
#include "fldtype.h"
#include "blopt8a.h"

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

!-----------------------------------------------------------------------

!  Workspace :-

      REAL                                                              &
     & VFRAC_TILE(LAND_PTS,NTILES)                                      &
!                               ! Fractional canopy coverage for
!                               ! land tiles.
!     &,HCONS(LAND_PTS)                                                  &
                                ! Soil thermal conductivity including
     &,FLANDG_X(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y)            &
!                               ! Land fraction on P-grid
!                               ! the effects of water and ice (W/m2)
     &,CDR10M(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y)
!                               ! Ratio of CD's reqd for calculation
!                               ! of 10 m wind. On P-grid

!  Expanded wind arrays - only needed for Buddy_sea coastal tiling
      REAL                                                              &
     & U_1_PX(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y)              &
     &,V_1_PX(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y)              &
     &,U_0_PX(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y)              &
     &,V_0_PX(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y)

!     Factors used to weight windspeeds over coastal points
      REAL                                                              &
     &  FLANDFAC(ROW_LENGTH,ROWS)                                       &
     & ,FLANDFAC_U(ROW_LENGTH,ROWS)                                     &
     & ,FLANDFAC_V(ROW_LENGTH,N_ROWS)                                   &
     & ,FLANDFAC_X(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y)         &
     & ,FSEAFAC(ROW_LENGTH,ROWS)                                        &
     & ,FSEAFAC_U(ROW_LENGTH,ROWS)                                      &
     & ,FSEAFAC_V(ROW_LENGTH,N_ROWS)                                    &
     & ,FSEAFAC_X(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y)

      REAL                                                              &
     & SECS_PER_360DAYS         ! LOCAL Number of seconds in 360 days

      PARAMETER(SECS_PER_360DAYS=31104000.0)


!  Local scalars :-

      INTEGER                                                           &
     & I,J,K,L                                                          &
                  ! LOCAL Loop counter (horizontal field index).
     &,IS,JS                                                            &
                  ! Loop counter for coastal point stencil
     &,N                                                                &
                  ! LOCAL Loop counter (tile index).
     &,COUNT      ! Counter for average wind speed

      REAL                                                              &
     & USHEAR                                                           &
                    ! U-component of surface-to-lowest-level wind shear.
     &,VSHEAR                                                           &
                    ! V-component of surface-to-lowest-level wind shear.
     &,VSHR2        ! Square of magnitude of surface-to-lowest-level
!                   ! wind shear.

      REAL SEAWIND  ! average wind speed adjacent to coast
      REAL FSEAMAX  ! Maximum factor to apply to coast wind speed

     ! Minimum factor allowed to convert coastal wind speed to land part
      REAL FLANDMIN
      PARAMETER(FLANDMIN=0.2)

      ! end step of experiment,  step width, processor num
      integer :: endstep, timestep_number, mype

     IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SF_EXPL ',3)
      ENDIF

! Expand land fraction to global field:
      DO J=1,ROWS
       DO I=1,ROW_LENGTH
           FLANDG(I,J)=0.0
       ENDDO
      ENDDO
      DO L=1,LAND_PTS
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
        FLANDG(I,J)=FLAND(L)
      ENDDO

!-----------------------------------------------------------------------
! Call TILEPTS to calculate TILE_PTS and TILE_INDEX for surface types
!-----------------------------------------------------------------------
! DEPENDS ON: tilepts
      CALL TILEPTS(LAND_PTS,FRAC,TILE_PTS,TILE_INDEX)
       
!-----------------------------------------------------------------------
! Calculate wind shear between level 1 and the surface
!-----------------------------------------------------------------------

      DO J=1,ROWS
        DO I=1,ROW_LENGTH
        IF(FLANDG(I,J) <  1.0)THEN
          USHEAR = U_1_P(I,J) - U_0_P(I,J)
          VSHEAR = V_1_P(I,J) - V_0_P(I,J)
          VSHR2 = MAX (1.0E-6 , USHEAR*USHEAR + VSHEAR*VSHEAR)
          VSHR_SSI(I,J) = SQRT(VSHR2)
        ELSE
          VSHR_SSI(I,J) = 0.0
        ENDIF
!
        IF(FLANDG(I,J) >  0.0)THEN
        VSHR2 = MAX (1.0E-6 , U_1_P(I,J)*U_1_P(I,J)                     &
     &    + V_1_P(I,J)*V_1_P(I,J))
        VSHR_LAND(I,J) = SQRT(VSHR2)
        ELSE
          VSHR_LAND(I,J) = 0.0
        ENDIF
!
        VSHR(I,J)= FLANDG(I,J)*VSHR_LAND(I,J)                           &
     &    + (1.0 - FLANDG(I,J))*VSHR_SSI(I,J)
        ENDDO
      ENDDO

#if !defined(SCMA)

      DO J= 1, rows
      DO I= 1, row_length
        FLANDG_X(I,J)  = FLANDG(I,J)
      ENDDO
      ENDDO

! DEPENDS ON: swap_bounds
      Call Swap_Bounds(FLANDG_X, row_length, rows, 1,                   &
     &                         off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(FLANDG_X,ROW_LENGTH,ROWS,1,              &
     &                         off_x,off_y)

      If (L_CTILE .AND. Buddy_sea == ON) Then

        DO J= 1, rows
        DO I= 1, row_length
          U_1_PX(I,J) = U_1_P(I,J)
          V_1_PX(I,J) = V_1_P(I,J)
          U_0_PX(I,J) = U_0_P(I,J)
          V_0_PX(I,J) = V_0_P(I,J)
        ENDDO
        ENDDO

! DEPENDS ON: swap_bounds
        Call Swap_Bounds(U_1_PX, row_length, rows, 1,                   &
     &                         off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: swap_bounds
        Call Swap_Bounds(V_1_PX, row_length, rows, 1,                   &
     &                         off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: swap_bounds
        Call Swap_Bounds(U_0_PX, row_length, rows, 1,                   &
     &                         off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: swap_bounds
        Call Swap_Bounds(V_0_PX, row_length, rows, 1,                   &
     &                         off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(U_1_PX,ROW_LENGTH,ROWS,1,              &
     &                         off_x,off_y)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(V_1_PX,ROW_LENGTH,ROWS,1,              &
     &                         off_x,off_y)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(U_0_PX,ROW_LENGTH,ROWS,1,              &
     &                         off_x,off_y)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(V_0_PX,ROW_LENGTH,ROWS,1,              &
     &                         off_x,off_y)

        DO J=1,ROWS
        DO I=1,ROW_LENGTH
          FSEAFAC(I,J)  = 1.0
          FLANDFAC(I,J) = 1.0

          IF ( FLANDG(I,J) > 0.01 .AND. FLANDG(I,J) < 0.99 ) THEN
!           !-----------------------------------------------------
!           ! Calculate average windspeed over adjacent sea points
!           !-----------------------------------------------------
            SEAWIND=0.0
            COUNT = 0
            DO IS=I-1,I+1
            DO JS=J-1,J+1
              IF( FLANDG_X(IS,JS) < 0.001 )THEN
!               ! ie. this is basically a sea point
                USHEAR = U_1_PX(IS,JS) - U_0_PX(IS,JS)
                VSHEAR = V_1_PX(IS,JS) - V_0_PX(IS,JS)
                VSHR2 = MAX (1.0E-10 , USHEAR*USHEAR + VSHEAR*VSHEAR)
                SEAWIND = SEAWIND + SQRT( VSHR2 )
                COUNT = COUNT + 1
              ENDIF
            END DO
            END DO
!           !-----------------------------------------------------
!           ! Calculate multiplicative factor, FSEAFAC, to convert
!           ! from the GBM VSHR to an appropriate marine VSHR
!           !-----------------------------------------------------
            IF (COUNT > 0) THEN
              SEAWIND = SEAWIND/FLOAT(COUNT)
!             ! Restrict FSEAFAC so FLANDFAC>FLANDMIN
              FSEAMAX = MIN( 1./FLANDMIN,                               &
     &                      (1.-FLANDMIN*FLANDG(I,J))/(1.-FLANDG(I,J)) )
!             ! First limit is to keep fseamax sensible as FLANDG -> 1
!             ! Second limit is to keep fland > flandmin, remembering 
!             !   that the we want FLANDG-weighted sum of factors =1
!             !   to preserve the gridbox mean VSHR
              FSEAFAC(I,J) = MAX(1.0,                                   &
     &                       MIN( FSEAMAX, SEAWIND/VSHR(I,J) ))
            ENDIF

            VSHR_SSI(I,J) = VSHR(I,J) * FSEAFAC(i,j)

            FLANDFAC(i,j) = ( 1.0 - FSEAFAC(i,j)*(1.0-FLANDG(I,J)) )    &
     &                    / FLANDG(I,J)
            VSHR_LAND(I,J) = VSHR(I,J) * FLANDFAC(i,j)


            VSHR(I,J)= FLANDG(I,J)*VSHR_LAND(I,J)                       &
     &        + (1.0 - FLANDG(I,J))*VSHR_SSI(I,J)
            ENDIF
        ENDDO
        ENDDO

      End if  ! test on buddy_sea switch
#endif
!-----------------------------------------------------------------------
! Call MOSES II physiology routine to calculate surface conductances
! and carbon fluxes.
!-----------------------------------------------------------------------

! DEPENDS ON: physiol
       if ( .not. l_cable) CALL PHYSIOL (                               &
     & ROW_LENGTH,ROWS,LAND_PTS,LAND_INDEX,                             &
     & SM_LEVELS,NTILES,TILE_PTS,TILE_INDEX,                            &
     & DIM_CS1, DIM_CS2,                                                &
     & CO2_MMR,CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,L_CO2_INTERACTIVE,        &
     & L_TRIFFID, L_Q10,                                                &
     & CAN_MODEL,CS,FRAC,CANHT_FT,PHOTOSYNTH_ACT_RAD,                   &
     & LAI_FT,PSTAR,QW_1,STHU,TIMESTEP,T_SOIL,TSTAR_TILE,               &
     & SMVCCL,SMVCST,SMVCWT,VSHR,Z0_TILE,Z1_UV,                         &
     & CANHC_TILE,VFRAC_TILE,FLAKE,                                     &
     & G_LEAF,GS,GC,GPP,GPP_FT,NPP,NPP_FT,                              &
     & RESP_P,RESP_P_FT,RESP_S,RESP_W_FT,SMC,WT_EXT_TILE,FSMC,          &
     & WT_EXT,RA,F_ROOT,ALBSOIL,COS_ZENITH_ANGLE                        &
     &,CAN_RAD_MOD, ILAYERS)

! Lestevens 18 March 2010: Remove moses fluxes by setting to zero
      if (l_cable) then
          NPP=0.0
          GPP=0.0
          NPP_FT=0.0
          GPP_FT=0.0
          RESP_S=0.0
          RESP_S_TOT=0.0
          RESP_P=0.0
          RESP_P_FT=0.0
         !RESP_W=0.0
         !RESP_W_FT=0.0
          G_LEAF=0.0
!         Need to replace the gs_tile calculation in physiol?
!         Just use zero to have something defined.
!         ACCESS-1.3-CABLE-1.8: used GS=GC=0 here:
!         GS = 0. ! Used by ACCESS-1.3-CABLE-1.8
!         GC = 0. ! Used by ACCESS-1.3-CABLE-1.8
!         ACCESS-1.4-CABLE-2.1: GS and GC are returned to the UM from CABLE (pfv, 25oct13)
      endif

!----------------------------------------------------------------------
! If TRIFFID is being used apply any correction to the land-atmosphere
! fluxes on the first timestep after the last TRIFFID call. Such a
! correction will typically be associated with a total depletion of
! carbon or with maintanence of the seed fraction. The corrections
! are stored in the accumulation variables after the call to TRIFFID.
! The correction is added to the instantaneous land-atmosphere fluxes
! (so that the atmospheric carbon budget is corrected) but is not
! included in the accumulation variables which drive TRIFFID, since
! this has already been dealt with during the last TRIFFID call.
!----------------------------------------------------------------------
      IF (L_TRIFFID .AND.(ASTEPS_SINCE_TRIFFID==1)                      &
     &    .AND. CycleNo==NumCycles) THEN
        DO N=1,NPFT
          DO L=1,LAND_PTS
            NPP_FT(L,N)=NPP_FT(L,N)+NPP_FT_ACC(L,N)/TIMESTEP
            RESP_P_FT(L,N)=RESP_P_FT(L,N)-NPP_FT_ACC(L,N)/TIMESTEP
            NPP_FT_ACC(L,N)=-NPP_FT_ACC(L,N)
          ENDDO
        ENDDO
        DO L=1,LAND_PTS
           RESP_S(L,1)=RESP_S(L,1)+RESP_S_ACC(L,1)/TIMESTEP
           RESP_S(L,2)=RESP_S(L,2)+RESP_S_ACC(L,2)/TIMESTEP
           RESP_S(L,3)=RESP_S(L,3)+RESP_S_ACC(L,3)/TIMESTEP
           RESP_S(L,4)=RESP_S(L,4)+RESP_S_ACC(L,4)/TIMESTEP
           RESP_S_ACC(L,1)=-RESP_S_ACC(L,1)
           RESP_S_ACC(L,2)=-RESP_S_ACC(L,2)
           RESP_S_ACC(L,3)=-RESP_S_ACC(L,3)
           RESP_S_ACC(L,4)=-RESP_S_ACC(L,4)
        ENDDO
      ENDIF

!----------------------------------------------------------------------
! Increment accumulation of leaf turnover rate.
! This is required for leaf phenology and/or TRIFFID, either of
! which can be enabled independently of the other.
!----------------------------------------------------------------------
      IF ( CycleNo == NumCycles ) THEN
      IF (L_PHENOL.OR.L_TRIFFID) THEN
        DO N=1,NPFT
          DO L=1,LAND_PTS
            G_LEAF_ACC(L,N) = G_LEAF_ACC(L,N) +                         &
     &      G_LEAF(L,N)*(TIMESTEP/SECS_PER_360DAYS)
          ENDDO
        ENDDO
      ENDIF

!----------------------------------------------------------------------
! Increment accumulation prognostics for TRIFFID
!----------------------------------------------------------------------
      IF (L_TRIFFID) THEN
        DO N=1,NPFT
          DO L=1,LAND_PTS
            NPP_FT_ACC(L,N) = NPP_FT_ACC(L,N) + NPP_FT(L,N)*TIMESTEP
            RESP_W_FT_ACC(L,N) = RESP_W_FT_ACC(L,N)                     &
     &                                      + RESP_W_FT(L,N)*TIMESTEP
          ENDDO
        ENDDO
        DO L=1,LAND_PTS
          RESP_S_ACC(L,1)=RESP_S_ACC(L,1)+RESP_S(L,1)*TIMESTEP
          RESP_S_ACC(L,2)=RESP_S_ACC(L,2)+RESP_S(L,2)*TIMESTEP
          RESP_S_ACC(L,3)=RESP_S_ACC(L,3)+RESP_S(L,3)*TIMESTEP
          RESP_S_ACC(L,4)=RESP_S_ACC(L,4)+RESP_S(L,4)*TIMESTEP
        ENDDO
      ENDIF
      ENDIF ! CycleNo == NumCycles

!-----------------------------------------------------------------------
! calculate CO2:(BIO+HUM) ratio, dependent on soil clay content, and
! sum soil respiration components
! (RESP_FRAC here then contains the fraction of soil respiration which
! is respired to the atmos. the rest is re-partitioned into BIO+HUM)
!
! RESP_S_ACC contains the full amount, and this is carried forward to
! VEG_CTL for use in updating soil carbon pools. RESP_S_TOT calculated
! here is passed to BL_TRMIX as the fraction which is respired as CO2
! to the atmosphere. RESP_S_TOT, and RESP_S are also passed out for
! storage in diagnostics 3293, and 3467-470.
!
!-----------------------------------------------------------------------
      if (L_TRIFFID) then
        DO I=1,LAND_PTS
          WORK_CLAY(I) = exp(-0.0786 * 100.0*clay_land(i))
          RESP_FRAC(I) = (3.0895 + 2.672*WORK_CLAY(I)) /                &
     &                   (4.0895 + 2.672*WORK_CLAY(I))
          RESP_S(I,1)  = RESP_S(I,1) * RESP_FRAC(I)
          RESP_S(I,2)  = RESP_S(I,2) * RESP_FRAC(I)
          RESP_S(I,3)  = RESP_S(I,3) * RESP_FRAC(I)
          RESP_S(I,4)  = RESP_S(I,4) * RESP_FRAC(I)
          RESP_S_TOT(I) = RESP_S(I,1) + RESP_S(I,2) +                   &
     &                    RESP_S(I,3) + RESP_S(I,4)
        ENDDO
      ENDIF

!-----------------------------------------------------------------------
! Reset TILE_PTS and TILE_INDEX and set tile fractions to 1 if aggregate
! tiles are used (NTILES=1).
! Otherwise, set tile fractions to surface type fractions.
!-----------------------------------------------------------------------

      IF (NTILES == 1) THEN
        TILE_PTS(1) = LAND_PTS
        DO L=1,LAND_PTS
          TILE_FRAC(L,1) = 1.
          TILE_INDEX(L,1) = L
        ENDDO
      ELSE
        DO N=1,NTYPE
          DO L = 1, LAND_PTS
            TILE_FRAC(L,N) = FRAC(L,N)
            VFRAC_TILE(L,N) = FRAC(L,N)
            !kdcorbin, 08/10 - added CABLE test/initialization
            IF (L_CABLE) Then
               FLAKE(L,N) = 0.
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      IF (LAND_PTS >  0) THEN    ! Omit if no land points

!-----------------------------------------------------------------------
! Calculate the thermal conductivity of the top soil layer.
!-----------------------------------------------------------------------
! DEPENDS ON: heat_con
        CALL HEAT_CON (LAND_PTS,HCON,STHU,STHF,SMVCST,HCONS,LTIMER)

      ENDIF                     ! End test on land points

!-----------------------------------------------------------------------
!! Calculate net radiation on land tiles and sea-ice
!! The calculation on land tiles is repeated for LRAD_EMIS_LAND_GEN
!! equals false in order to achieve bit reproducibility.
!! It would be tidier to move the if statement within the do loop
!! but this affects vectorisation.
!-----------------------------------------------------------------------
 
          RADNET_TILE(:,:) = 0.
          LE_TILE(:,:) = 0.

      If (LRAD_EMIS_LAND_GEN) Then
        DO N=1,NTILES
          DO K=1,TILE_PTS(N)
            L = TILE_INDEX(K,N)
            J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
            I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
! EAK   ****************
!         should be IF( .NOT. l_cable ) RADNET ...
            RADNET_TILE(L,N) = SW_TILE(L,N) +                           &
     &          RAD_EMIS_LAND_GEN*(LW_DOWN(I,J) - SBCON*T_SOIL(L,1)**4)
          ENDDO
        ENDDO
      
      Else
      
        DO N=1,NTILES
          DO K=1,TILE_PTS(N)
            L = TILE_INDEX(K,N)
            J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
            I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
            RADNET_TILE(L,N) = SW_TILE(L,N) +                           &
     &             LW_DOWN(I,J) - SBCON*T_SOIL(L,1)**4
          ENDDO
        ENDDO
      Endif
         
      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        RADNET_SICE(I,J) = 0.
        IF (FLANDG(I,J) <  1.0 .AND. ICE_FRACT(I,J) >  0.)              &
     &    RADNET_SICE(I,J) = RAD_SICE(I,J) -                            &
     &                       ICE_FRACT(I,J)*SBCON*TI(I,J)**4
       ENDDO
      ENDDO

!-----------------------------------------------------------------------
!! 4.  Surface turbulent exchange coefficients and "explicit" fluxes
!!     (P243a, routine SF_EXCH).
!!     Wind mixing "power" and some values required for other, later,
!!     diagnostic calculations, are also evaluated if requested.
!-----------------------------------------------------------------------

! DEPENDS ON: sf_exch
      CALL SF_EXCH (                                                    &
     & ROW_LENGTH,ROWS,OFF_X,OFF_Y,HALO_I,HALO_J,                       &
     & LAND_PTS,NTILES,LAND_INDEX,                                      &
     & TILE_INDEX,TILE_PTS,FLAND,FLANDG,                                &
     & BQ_1,BT_1,CANHC_TILE,CANOPY,CATCH,DZSOIL(1),FLAKE,GC,HCONS,      &
     & CAN_MODEL,CATCH_SNOW,lq_mix_bl,                                  &
     & HO2R2_OROG,ICE_FRACT,SNOW_TILE,PSTAR,QW_1,RADNET_SICE,           &
     & RADNET_TILE,SIL_OROG_LAND,SMVCST,TILE_FRAC,TIMESTEP,             &
     & TL_1,TI,T_SOIL(1,1), COR_UST, COR_MO_ITER,                       &
     & TSTAR_TILE,TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,TSTAR_SSI,Z_LAND,     &
     & L_CTILE,SeaSalinityFactor,ISeaZ0T,I_SCRN_T_DIAG,                 &
     & TSTAR,L_spec_z0,Z0M_SCM,Z0H_SCM,L_DUST,                          &
     & VFRAC_TILE,VSHR_LAND,VSHR_SSI,ZH,Z0_TILE,Z1_UV,Z1_TQ,LAND_MASK,  &
     & SU10,SV10,SQ1P5,ST1P5,SFME,SZ0HEFF,LTIMER,FORMDRAG,FD_stab_dep,  &
     & OROG_DRAG_PARAM,Z0MSEA,                                          &
     & ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,CD,CH,                       &
     & RECIP_L_MO_SEA,CDR10M,CHR1P5M,                                   &
     & CHR1P5M_SICE,E_SEA,FME,FQW_1,FQW_TILE,EPOT_TILE,FQW_ICE,         &
     & FTL_1,FTL_TILE,FTL_ICE,FRACA,H_BLEND_OROG,H_SEA,Charnock,        &
     & RHOSTAR,RESFS,RESFT,RIB,RIB_TILE,                                &
     & FB_SURF,U_S,Q1_SD,T1_SD,Z0HSSI,Z0H_TILE,Z0H_EFF,                 &
     & Z0M_GB,Z0MSSI,Z0M_TILE,Z0M_EFF,RHO_ARESIST,ARESIST,RESIST_B,     &
     & RHO_ARESIST_TILE,ARESIST_TILE,RESIST_B_TILE,                     &
     & R_B_DUST,CD_STD_DUST,U_S_STD_TILE,                               &
     & RHO_CD_MODV1,RHOKH_TILE,RHOKH_SICE,RHOKM_1,RHOKM_LAND,RHOKM_SSI, &
     & RHOKPM,RHOKPM_POT,RHOKPM_SICE,BL_diag,ANTHROP_HEAT               &
! EAK
     &,l_cable                                                          & 
     &,SM_LEVELS,N_ROWS,LAND_PTS_TRIF,NPFT_TRIF,DIM_CS1,DIM_CS2,        &
     & SMVCCL,SMVCWT,                                                   &
     & CO2_MMR,CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,L_CO2_INTERACTIVE,        &
     & SW_TILE,LW_DOWN,SW_DOWN,                                         &
     & U_1,V_1,U_1_P,V_1_P,                                             &
     & surf_down_sw,alb_tile,cos_zenith_angle,l_tile_pts,               &
!     & surf_down_sw,alb_tile,cos_zenith_angle,               &
     & ls_rain,ls_snow,                                                 &
     & lat,long,day,time_sec,                                           &
!     & )
     & SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L,                   &
     & SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE,             &
     & TSOIL_TILE,T_SURF_TILE,BEXP,                                     &
     & SATHH,SATCON,HCAP,HCON,                                          &
!     & SOIL_TYPE,VEG_TYPE,                                              &
     & SOIL_TYPE,VEG_TYPE,albsoil,                                      &
     & ISNOW_FLG3L,total_nsteps,                                        &
     &            FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB,              &
     &            TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB,          &
     &            USTAR_CAB,SURF_HTF_CAB,                               &
!     & )
     & SOIL_LAYER_MOISTURE,STHF,STHU,T_SOIL,                            &
     & CS,FRAC,CANHT_FT,PHOTOSYNTH_ACT_RAD,LAI_FT,                      &
     & G_LEAF,GS,                                                       &
     & LE_TILE,                                                         &
!     & RHOKM_U_1,RHOKM_V_1,                                            &
     & TAUX_LAND,TAUX_SSI,TAUY_LAND,TAUY_SSI,                           &
     & CDR10M_U,CDR10M_V,                                               &
     & GPP,NPP,RESP_P,GPP_FT,NPP_FT,                                    &
     & RESP_P_FT,RESP_S,RESP_S_TOT,RESP_W_FT,                           &
     ! Lestevens 23apr13
     & NPP_FT_ACC,RESP_W_FT_ACC,                                        &
     & WT_EXT_TILE,FSMC,                                                &
     & FLANDG_U,FLANDG_V,                                               &
     & U_S_CAB,CH_CAB,CD_CAB,                                           &
!sxy
     & SNAGE_TILE,RTSOIL_TILE,                                          &
     & GFLUX_TILE,SGFLUX_TILE,                                          &
     & F_ROOT, sin_theta_latitude,                                      &
     & CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,SOIL_ORDER,                     &
     & NIDEP,NIFIX,PWEA,PDUST,GLAI,PHENPHASE,                           &
      ! end step of experiment, this step, step width, processor num
     & endstep, timestep_number, mype                                   &
     & )
!
#if !defined(SCMA)

! DEPENDS ON: swap_bounds
      Call Swap_Bounds(RHOKM_1, row_length, rows,                       &
     &                 1 , off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: swap_bounds
      Call Swap_Bounds(RHOKM_LAND, row_length, rows,                    &
     &                 1 , off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: swap_bounds
      Call Swap_Bounds(RHOKM_SSI, row_length, rows,                     &
     &                 1 , off_x, off_y, fld_type_p, .false.)

! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(RHOKM_1,ROW_LENGTH,ROWS,1,               &
     &                         off_x,off_y)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(RHOKM_LAND,ROW_LENGTH,ROWS,1,            &
     &                         off_x,off_y)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(RHOKM_SSI,ROW_LENGTH,ROWS,1,             &
     &                         off_x,off_y)

! DEPENDS ON: p_to_u
      CALL P_TO_U(RHOKM_1,ROW_LENGTH,ROWS,1,                            &
     &            off_x, off_y, RHOKM_U_1)
! DEPENDS ON: p_to_u_land
      CALL P_TO_U_LAND(RHOKM_LAND,FLANDG_X,ROW_LENGTH,ROWS,1,           &
     &            off_x, off_y, RHOKM_U_LAND)
! DEPENDS ON: p_to_u_sea
      CALL P_TO_U_SEA(RHOKM_SSI,FLANDG_X,ROW_LENGTH,ROWS,1,             &
     &            off_x, off_y, RHOKM_U_SSI)

! DEPENDS ON: p_to_v
      CALL P_TO_V(RHOKM_1,ROW_LENGTH, ROWS, n_rows,                     &
     &            1, off_x, off_y, RHOKM_V_1)
! DEPENDS ON: p_to_v_land
      CALL P_TO_V_LAND(RHOKM_LAND,FLANDG_X,ROW_LENGTH, ROWS, n_rows,    &
     &            1, off_x, off_y, RHOKM_V_LAND)
! DEPENDS ON: p_to_v_sea
      CALL P_TO_V_SEA(RHOKM_SSI,FLANDG_X,ROW_LENGTH, ROWS, n_rows,      &
     &            1, off_x, off_y, RHOKM_V_SSI)

! DEPENDS ON: p_to_u
      CALL P_TO_U(FLANDG_X,ROW_LENGTH,ROWS,1,                           &
     &            off_x, off_y, FLANDG_U)
! DEPENDS ON: p_to_v
      CALL P_TO_V(FLANDG_X,ROW_LENGTH, ROWS, n_rows,                    &
     &            1, off_x, off_y, FLANDG_V)

      If (L_CTILE .AND. Buddy_sea == ON) Then
!       ! Interpolate wind speed factors to u and v columns

        DO J= 1, rows
        DO I= 1, row_length
          FLANDFAC_X(I,J)= FLANDFAC(I,J)
          FSEAFAC_X(I,J) = FSEAFAC(I,J)
        ENDDO
        ENDDO

! DEPENDS ON: swap_bounds
        Call Swap_Bounds(FLANDFAC_X, row_length, rows, 1,               &
     &                         off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(FLANDFAC_X,ROW_LENGTH,ROWS,1,          &
     &                         off_x,off_y)
! DEPENDS ON: p_to_u
        CALL P_TO_U(FLANDFAC_X,ROW_LENGTH,ROWS,1,                       &
     &                         off_x, off_y, FLANDFAC_U)
! DEPENDS ON: p_to_v
        CALL P_TO_V(FLANDFAC_X,ROW_LENGTH, ROWS, n_rows, 1,             &
     &                         off_x, off_y, FLANDFAC_V)
! DEPENDS ON: swap_bounds
        Call Swap_Bounds(FSEAFAC_X, row_length, rows, 1,                &
     &                         off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(FSEAFAC_X,ROW_LENGTH,ROWS,1,           &
     &                         off_x,off_y)
! DEPENDS ON: p_to_u
        CALL P_TO_U(FSEAFAC_X,ROW_LENGTH,ROWS,1,                        &
     &                         off_x, off_y, FSEAFAC_U)
! DEPENDS ON: p_to_v
        CALL P_TO_V(FSEAFAC_X,ROW_LENGTH, ROWS, n_rows, 1,              &
     &                         off_x, off_y, FSEAFAC_V)

     END IF  ! test on Buddy_sea

#else
      DO J= 1, rows
       DO I= 1, row_length
         RHOKM_U_1(I,J) = RHOKM_1(I,J)
         RHOKM_U_LAND(I,J) = RHOKM_LAND(I,J)
         RHOKM_V_LAND(I,J) = RHOKM_LAND(I,J)
         RHOKM_U_SSI(I,J) = RHOKM_SSI(I,J)
         RHOKM_V_SSI(I,J) = RHOKM_SSI(I,J)
         RHOKM_V_1(I,J) = RHOKM_1(I,J)
         FLANDG_U(I,J) = FLANDG(I,J)
         FLANDG_V(I,J) = FLANDG(I,J)
         FLANDFAC_U(I,J) = 1.0
         FLANDFAC_V(I,J) = 1.0
         FSEAFAC_U(I,J)  = 1.0
         FSEAFAC_V(I,J)  = 1.0
       END DO
      END DO
#endif

      If (L_CTILE .AND. Buddy_sea == ON) Then

        DO J=1,ROWS
        DO I=1,ROW_LENGTH
          TAUX_LAND(I,J) = RHOKM_U_LAND(I,J)* U_1(I,J) * FLANDFAC_U(I,J)
          TAUX_SSI(I,J)  = RHOKM_U_SSI(I,J) * ( U_1(I,J) - U_0(I,J) )   &
     &                                      * FSEAFAC_U(I,J)
          TAUX_1(I,J) = FLANDG_U(I,J)*TAUX_LAND(I,J)                    &
     &                  + (1.-FLANDG_U(I,J))*TAUX_SSI(I,J)
        END DO
        END DO

        DO J=1,N_ROWS
        DO I=1,ROW_LENGTH
          TAUY_LAND(I,J) = RHOKM_V_LAND(I,J)* V_1(I,J) * FLANDFAC_V(I,J)
          TAUY_SSI(I,J)  = RHOKM_V_SSI(I,J) * ( V_1(I,J) - V_0(I,J) )   &
     &                                      * FSEAFAC_V(I,J)
          TAUY_1(I,J) = FLANDG_V(I,J)*TAUY_LAND(I,J)                    &
     &                  + (1.-FLANDG_V(I,J))*TAUY_SSI(I,J)
        END DO
        END DO

      Else   ! Standard code

        DO J=1,ROWS
        DO I=1,ROW_LENGTH
          TAUX_LAND(I,J) = RHOKM_U_LAND(I,J) * U_1(I,J)
          TAUX_SSI(I,J) = RHOKM_U_SSI(I,J) * ( U_1(I,J) - U_0(I,J) )

          TAUX_1(I,J) = FLANDG_U(I,J)*TAUX_LAND(I,J)                    &
     &                  + (1.-FLANDG_U(I,J))*TAUX_SSI(I,J)
        END DO
        END DO

        DO J=1,N_ROWS
        DO I=1,ROW_LENGTH
         TAUY_LAND(I,J) = RHOKM_V_LAND(I,J) * V_1(I,J)
         TAUY_SSI(I,J) = RHOKM_V_SSI(I,J) * ( V_1(I,J) - V_0(I,J) )

         TAUY_1(I,J) = FLANDG_V(I,J)*TAUY_LAND(I,J)                     &
     &                 + (1.-FLANDG_V(I,J))*TAUY_SSI(I,J)
        END DO
        END DO

      End If


#if !defined(SCMA)
      If (SU10 .or. SV10) Then
! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &                 CDR10M, row_length, rows,                        &
     &                 1, off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
      Call FILL_EXTERNAL_HALOS(CDR10M,ROW_LENGTH,ROWS,1,off_x,off_y)
      End If

      IF (SU10)THEN
! DEPENDS ON: p_to_u
        CALL P_TO_U(CDR10M,ROW_LENGTH,ROWS,1,                           &
     &              off_x,off_y,                                        &
     &              CDR10M_U)

      ENDIF

      IF (SV10)THEN
! DEPENDS ON: p_to_v
        CALL P_TO_V(CDR10M,ROW_LENGTH,ROWS,N_ROWS,1,                    &
     &            off_x,off_y, CDR10M_V)
      ENDIF
#else
       DO J= 1, rows
        DO I= 1, row_length
          CDR10M_U(i,j) = CDR10M(i,j)
          CDR10M_V(i,j) = CDR10M(i,j)
        enddo
       enddo
#endif


!-----------------------------------------------------------------------
!! Set grid-box surface exchange coefficients
!-----------------------------------------------------------------------
      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        IF( FLANDG(I,J)  <   1.0) THEN
          RHOKH(I,J) = (1.0 - FLANDG(I,J))*RHOKH_SICE(I,J)
        ELSE
          RHOKH(I,J) = 0.0
        ENDIF
       ENDDO
      ENDDO

      DO N=1,NTILES
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          RHOKH(I,J) = RHOKH(I,J)                                       &
     &      + FLANDG(I,J)*TILE_FRAC(L,N)*RHOKH_TILE(L,N)
        ENDDO
      ENDDO


      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SF_EXPL ',4)
      ENDIF

      RETURN
      END SUBROUTINE SF_EXPL
#endif
