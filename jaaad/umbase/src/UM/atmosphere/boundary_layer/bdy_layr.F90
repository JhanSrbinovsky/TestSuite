#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE BDY_LAYR-----------------------------------------------
!!!
!!!  Purpose: Calculate turbulent fluxes of heat, moisture and momentum
!!!           between (a) surface and atmosphere, (b) atmospheric levels
!!!           within the boundary layer, and/or the effects of these
!!!           fluxes on the primary model variables.  The flux of heat
!!!           into and through the soil is also modelled.  Numerous
!!!           related diagnostics are also calculated.
!!!
!!!  Model            
!!! Programming standard : unified model documentation paper No 3
!!!
!!!  Documentation: UMDP 24.
!!!
!!!---------------------------------------------------------------------

!    Arguments :-
      SUBROUTINE BDY_LAYR (                                             &

! IN MPP variables
     & halo_i, halo_j, off_x, off_y, row_length, rows, n_rows,          &
     & global_row_length, proc_row_group, at_extremity,                 &
     & n_proc, n_procx, n_procy, neighbour,                             &
     & DIM_CS1, DIM_CS2,                                                &

! IN  parameters for iterative SISL scheme
     & NumCycles, CycleNo,                                              &

! IN Substep number for ATMPHYS2
     & Substep_Number, Num_Substeps, L_phys2_substep,                   &

! IN values defining field dimensions and subset to be processed :
     & NTILES,land_pts,                                                 &

! IN values defining vertical grid of model atmosphere :
     & MODEL_DOMAIN,                                                    &
     & BL_LEVELS,                                                       &
     & r_rho_levels,r_theta_levels,eta_theta_levels,                    &
     & P,P_theta_levels, rho_rsq, rho_wet, rho_dry,                     &
     & GAMMA,SIN_THETA_LONGITUDE,COS_THETA_LONGITUDE,                   &
     & sin_theta_latitude,                                              &

! IN U, V and W momentum fields.
     & U, V, W, ETADOT,                                                 &
! IN Non turbulent increments to momentum
!  (New dynamics only).
     & DU_NT,DV_NT,                                                     &
! variables for subgrid turbulence scheme
     & visc_bl_m,delta_lambda, delta_phi,FM_3D, FH_3D,L_subfilter_vert, &
     & L_subfilter_horiz, L_subfilter_blend,max_diff,turb_startlev_vert,&
     & turb_endlev_vert,BL_COEF_KM, BL_COEF_KH,                         &

! IN soil/vegetation/land surface data :
     & LAND_MASK,LAND_INDEX,ST_LEVELS,SM_LEVELS,                        &
     & CANOPY,CATCH,CATCH_SNOW,HCON,SMVCCL,SMVCST,SMVCWT,STHF,STHU,     &
     & SIL_OROG_LAND,FORMDRAG,OROG_DRAG_PARAM,HO2R2_OROG,               &
     & SOIL_LAYER_MOISTURE,                                             &

! IN sea/sea-ice data :
     & ICE_FRACT,U_0,V_0,U_0_P,V_0_P,Charnock,SeaSalinityFactor,        &

! IN cloud data :
     & CF,Q,QCF,QCL,CCA,CCB,CCT, T,                                     &

! IN everything not covered so far :
     & L_LAMBDAM2,L_FULL_LAMBDAS,                                       &
     & CO2_MMR,PHOTOSYNTH_ACT_RAD,PSTAR,RAD_SICE,                       &
     & RAD_HR,MICRO_TENDS,TIMESTEP,lq_mix_bl,ZH_PREV,                   &
     & BL_OPTIONS,L_SBLeq,L_SBLco,Muw_SBL,Mwt_SBL,                      &
! IN SCM variables
     & flux_e, flux_h, L_flux_bc, L_spec_z0, Z0M_SCM, Z0H_SCM,          &

! IN variables required for mineral dust scheme
     & L_DUST,L_CAM_DUST,                                               &
     & SOIL_CLAY,SOIL_SILT,SOIL_SAND,                                   &
     & DUST_MREL1,DUST_MREL2,DUST_MREL3,                                &
     & DUST_MREL4,DUST_MREL5,DUST_MREL6,                                &

! IN additional variables for MOSES II
     & LAND_PTS_TRIF,NPFT_TRIF,                                         &
     & SNOW_TILE,Z0_TILE,LW_DOWN,SW_TILE,                               &
     & CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,L_CO2_INTERACTIVE,                &
     & L_PHENOL,L_TRIFFID,L_Q10,ASTEPS_SINCE_TRIFFID,CAN_MODEL,         &
     & CS,FRAC,CANHT_FT,LAI_FT,                                         &
     & FLAND,FLANDG,TSTAR_SEA,Z_LAND,L_CTILE,ANTHROP_HEAT,              &
     & ALBSOIL,COS_ZENITH_ANGLE,CAN_RAD_MOD,ILAYERS,                    &

!    EAK      
!    IN
     & surf_down_sw,alb_tile,l_tile_pts,                                &
!     & surf_down_sw,alb_tile,cos_zenith_angle,               &
     & ls_rain,ls_snow,SW_DOWN,                                         &
     & lat,long,day,time_sec,                                           &
     & SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L,                   &
     & SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE,             &
     & TSOIL_TILE,T_SURF_TILE,HCONS,BEXP,                               &
     &  SATHH,SATCON,HCAP,                                              &
!     & SOIL_TYPE,VEG_TYPE,                                             &
     & SOIL_TYPE,VEG_TYPE,                                              &
     & ISNOW_FLG3L,total_nsteps,                                        &
     &            FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB,              &
     &            TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB,          &
     &            USTAR_CAB,SURF_HTF_CAB,                               &
     &             U_S_CAB,CH_CAB,CD_CAB,                               &
     &            l_cable,                                              &
!sxy
     &           SNAGE_TILE,RTSOIL_TILE,                                &
     &           GFLUX_TILE,SGFLUX_TILE,                                &
! Lestevens Sept2012: CasaCNP variables 
     &           CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,                      &
     &           SOIL_ORDER,NIDEP,NIFIX,PWEA,PDUST,                     &
     &           GLAI,PHENPHASE,                                        &

! INOUT data :
     & GS,T_SOIL,TI,TSTAR,TSTAR_TILE,Z0MSEA,                            &

! INOUT additional variables for MOSES II
     & G_LEAF_ACC,NPP_FT_ACC,RESP_W_FT_ACC,RESP_S_ACC,                  &
     & TSTAR_LAND,TSTAR_SICE,TSTAR_SSI,                                 &

! INOUT data for STPH_RP
     & G0_RP,par_mezcla,                                                &
! OUT Diagnostic not requiring STASH flags :
     & CD,CH,E_SEA,FQW,FQW_TILE,EPOT_TILE,                              &
     & FTL,FTL_TILE,H_SEA,RHOKH,RHOKM,                                  &
     & RHOKM_u,RHOKM_v, RIB_GB,TAUX,TAUY,VSHR,zht,                      &
     & SHALLOWC,CU_OVER_OROG,                                           &
     & BL_TYPE_1,BL_TYPE_2,BL_TYPE_3,BL_TYPE_4,BL_TYPE_5,BL_TYPE_6,     &
     & BL_TYPE_7,                                                       &
     & Z0M_EFF_GB,Z0H_EFF_GB,                                           &

! OUT diagnostic requiring STASH flags :
     & FME,                                                             &

! OUT diagnostics required for soil moisture nudging scheme :
     & WT_EXT,RA,                                                       &

! (IN) STASH flags :-
     & SFME,SIMLT,SMLT,SLH,SQ1P5,ST1P5,SU10,SV10,SZ0HEFF,               &
!
! SCM Diagnostics (dummy values in full UM)
     & nSCMDpkgs,L_SCMDiags,                                            &

! OUT data required for tracer mixing :
     & RHO_ARESIST,ARESIST,RESIST_B,                                    &
     & NTML,                                                            &
     & KENT, WE_LIM, T_FRAC, ZRZI,                                      &
     & KENT_DSC, WE_LIM_DSC, T_FRAC_DSC, ZRZI_DSC,                      &
     & ZHSC, Z_HALF,                                                    &

!OUT variables required for mineral dust scheme
     & R_B_DUST,DUST_FLUX,                                              &
     & U_S_T_TILE,U_S_T_DRY_TILE,U_S_STD_TILE,                          &


! OUT data required for 4D-VAR :
     & RHO_CD_MODV1,RHO_KM,                                             &

! OUT variables required in IMP_SOLVER
     & ALPHA1,ASHTF,                                                    &
     & DTRDZ_CHARNEY_GRID,RDZ_CHARNEY_GRID,                             &
     & DTRDZ_U,DTRDZ_V,RDZ_U,RDZ_V,                                     &
     & FRACA,RHOKH_TILE,SMC,CHR1P5M,RESFS,Z0HSSI,Z0MSSI,                &
     & CDR10M_U,CDR10M_V,Z1_TQ,                                         &

! OUT variables which need to be maintained between substeps
     & RHO_UV,RHO_TQ,DZL_CHARNEY,RDZ,                                   &
     & Z1_UV,Z_FULL,Z_UV,Z_TQ,P_HALF,DELTAP,                            &

! OUT additional variables for MOSES II
     & LE_TILE,RADNET_SICE,                                             &
     & RADNET_TILE,RIB_TILE,RHO_ARESIST_TILE,ARESIST_TILE,              &
     & RESIST_B_TILE,ALPHA1_SICE,ASHTF_TILE,FQW_ICE,                    &
     & FTL_ICE,RESFT,RHOKH_SICE,RHOKPM,RHOKPM_POT,                      &
     & RHOKPM_SICE,Z0H_TILE,Z0M_GB,Z0M_TILE,CHR1P5M_SICE,               &
     & G_LEAF,GPP_FT,NPP_FT,                                            &
     & RESP_P_FT,RESP_S,RESP_S_TOT,RESP_W_FT,                           &
     & GC,CANHC_TILE,WT_EXT_TILE,FLAKE,                                 &
     & TILE_INDEX,TILE_PTS,TILE_FRAC,FSMC,                              &
     & RIB_SSI,TAUX_LAND,TAUX_SSI,TAUY_LAND,TAUY_SSI,VSHR_LAND,VSHR_SSI,&
     & FLANDG_U,FLANDG_V,                                               &

! OUT data required elsewhere in UM system :
     & ZH,GPP,NPP,RESP_P,T1_SD,Q1_SD,                                   &
     & NTPAR,NLCL,ZHPAR,Z_LCL,L_SHALLOW,                                &
     & NTDSC,NBDSC,CUMULUS,WSTAR,WTHVS,DELTHVU,                         &
     & UW0,VW0,                                                         &
     & ERROR,LTIMER,BL_diag                                             &
     &, L_ukca,                                                         &
      ! end step of experiment, this step, step width, processor num
     & endstep, timestep_number, mype                                   &
     & )

      Use bl_diags_mod, Only :                                          &
          strnewbldiag

      IMPLICIT NONE

#include "c_dust_ndiv.h"
#include "c_r_cp.h"
#include "c_g.h"
#include "c_lheat.h"
#include "blopt8a.h"

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
     &, neighbour(4)                                                    &
                      ! Array with the Ids of the four neighbours in
                   ! in the horizontal plane
     &,NumCycles                                                        &
                  ! Number of cycles (iterations) for iterative SISL.
     &,CycleNo                                                          &
                  ! Iteration no
     &, DIM_CS1, DIM_CS2      ! soil carbon dimensions

! IN Substep number for ATMPHYS2
      Integer                                                           &
     & Substep_Number                                                   &
     &,Num_Substeps

! Switch for calculating exchange coeffs from latest values.
      Logical                                                           &
     & L_phys2_substep
      LOGICAL                                                           &
     & lq_mix_bl

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         !  south,east or west of the processor grid

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
                                   ! IN No. of land-surface tiles
!                                  !    (MOSES II)
     &,land_pts                    ! IN No.of land points in whole grid.
! (b) Defining vertical grid of model atmosphere.

      INTEGER                                                           &
     & BL_LEVELS                                                        &
                                   ! IN Max. no. of "boundary" levels
!                                     allowed.Assumed <= 30 for dim-
!                                     sioning of GAMMA in common deck
!                                     C_GAMMA used in SF_EXCH and KMKH
     &,MODEL_DOMAIN                                                     &
     &,I,J,K,L,N

!
! Boundary Layer
      REAL, Intent(IN) :: SeaSalinityFactor
!                               ! Factor allowing for the effect
!                               ! of the salinity of sea water on
!                               ! the evaporative flux.
!
      REAL, Intent(IN) :: orog_drag_param
!                               ! Drag coefficient for orographic
!                               ! form drag
!
      INTEGER, DIMENSION(20) :: BL_OPTIONS   ! IN BL switches

      INTEGER                                                           &
     & FORMDRAG                    ! IN switch for orographic form drag
      REAL                                                              &
     &  Gamma(bl_levels)                                                &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j,bl_levels+1)                  &
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j,0:bl_levels+1)              &
     &, eta_theta_levels(0:bl_levels)                                   &
     &, COS_THETA_LONGITUDE (row_length,rows)                           &
     &, SIN_THETA_LONGITUDE (row_length,rows)                           &
     &, sin_theta_latitude(row_length,rows)                             &
     &, P(1-off_x:row_length+off_x,                                     &
     &    1-off_y:rows+off_y, BL_LEVELS+1)                              &
     &, P_theta_levels(row_length, rows, BL_LEVELS+1)                   &
     &, RHO_RSQ(1-off_x:row_length+off_x,                               &
     &          1-off_y:rows+off_y,BL_LEVELS+1)                         &
                                               ! IN Density * R**2
     &, rho_wet(row_length, rows, bl_levels+1)                          &
!                       ! wet density on rho levels (kg/m3)
     &, rho_dry(row_length, rows, bl_levels+1)                          &
!                       ! dry density on rho levels (kg/m3)
     &, DU_NT(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &        BL_LEVELS)                                                &
                            ! non-turbulent increment to u wind field
     &, DV_NT(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,           &
     &        BL_LEVELS)                                                &
                            ! non-turbulen increment to v wind field
     &, U(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      BL_LEVELS)                                                  &
     &, V(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      BL_LEVELS)                                                  &
     &, W(row_length, rows, 0:BL_LEVELS)                                &
     &, ETADOT(row_length, rows, 0:BL_LEVELS)                           &
     &, visc_BL_m(1:row_length, 1:rows, bl_levels)
!            ! visc_m only on BL levels
      Real                                                              &
     &  delta_lambda                                                    &
     &, delta_phi                                                       &
     &, FM_3D(row_length,rows,BL_LEVELS)                                &
!            ! stability function for momentum transport.
!            ! level 1 value is dummy for use in diagnostics
     &, FH_3D(row_length,rows,BL_LEVELS)                                &
!             ! stability function for heat and moisture.
!             ! level 1 value is dummy for use in diagnostics
     &, max_diff                                                        &
                  ! max diffusion coeff for run
     &, BL_COEF_KM(1:row_length, 1:rows, bl_levels-1)                   &
!             ! RHOKM from BL scheme
     &, BL_COEF_KH(1:row_length, 1:rows, bl_levels-1)
!             ! RHOKH from BL scheme

      logical                                                           &
     &  L_subfilter_vert                                                &
                            ! subgrid turbulence scheme in vertical     
     &, L_subfilter_horiz                                               &
                            ! subgrid turbulence scheme in horizontal   
     &, L_subfilter_blend   
                            ! blending BL and Smag coefficients    

      Integer                                                           &
     & turb_startlev_vert                                               &
                            ! start and end vertical levels for         
     &,turb_endlev_vert     
                            ! 3D turbulence scheme

! (c) Soil/vegetation/land surface parameters (mostly constant).

      LOGICAL                                                           &
     & LAND_MASK(row_length,rows)                                       &
                                   ! IN T if land, F elsewhere.
     &,L_CTILE                                                          &
                                  ! IN true if coastal tiling
     &,L_LAMBDAM2                                                       &
                   ! IN LambdaM=2*LambdaH (operational setting).
     &,L_FULL_LAMBDAS ! Lambdas NOT reduced above NTML_LOCAL+1

      INTEGER                                                           &
     & LAND_INDEX(land_pts)      ! IN LAND_INDEX(I)=J => the Jth
!                                     point in P_FIELD is the Ith
!                                     land point.

      INTEGER                                                           &
     & ST_LEVELS                                                        &
                                   ! IN No. of deep soil temp. levels
     &,SM_LEVELS                   ! IN No. of soil moisture levels

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
     &,HCON(land_pts)                                                   &
                                 ! IN Soil thermal conductivity
!                                     (W/m/K).
     &,SMVCCL(land_pts)                                                 &
                                 ! IN Critical volumetric SMC (m3/m3
!                                     of soil).
     &,SMVCST(land_pts)                                                 &
                                 ! IN Volumetric saturation point
!                                     (m3/m3 of soil).
     &,SMVCWT(land_pts)                                                 &
                                 ! IN Volumetric wilting point (m3/m3
!                                     of soil).
     &,STHF(land_pts,SM_LEVELS)                                         &
                                 ! IN Frozen soil moisture content of
!                                     each layer as a fraction of
!                                     saturation.
     &,STHU(land_pts,SM_LEVELS)                                         &
                                 ! IN Unfrozen soil moisture content
!                                     of each layer as a fraction of
!                                     saturation.
     &,SIL_OROG_LAND(land_pts)                                          &
                                 ! IN Silhouette area of unresolved
!                                     orography per unit horizontal area
!                                     on land points only.
     &,HO2R2_OROG(land_pts)                                             &
                                 ! IN Standard Deviation of orography.
!                                     equivilent to peak to trough
!                                     height of unresolved orography
!                                     devided by 2SQRT(2) on land
!                                     points only (m)
     &, SOIL_LAYER_MOISTURE(LAND_PTS,SM_LEVELS)                         &
                                               !IN soil moisture
!                 ! per layer (kg m-2)
!                 ! per layer (kg m-2)
     &,ZH_PREV(row_length, rows) ! IN boundary layer height from
!                                !    previous timestep

! (d) Sea/sea-ice data.

      REAL                                                              &
     & ICE_FRACT(row_length,rows)                                       &
                                   ! IN Fraction of gridbox covered by
!                                       sea-ice (decimal fraction).
     &,U_0(row_length,rows)                                             &
                                   ! IN W'ly component of surface
!                                       current (m/s).
     &,V_0(row_length,n_rows)                                           &
                                   ! IN S'ly component of surface
!                                       current (m/s).
     &,U_0_P(row_length,rows)                                           &
                                   ! IN W'ly component of surface
!                                       current (m/s). P grid
     &,V_0_P(row_length,rows)      ! IN S'ly component of surface
!                                       current (m/s). P grid

! (e) Cloud data.

      REAL                                                              &
     & CF(row_length, rows, BL_LEVELS)                                  &
                                          ! IN Cloud fraction (decimal).
     &,QCF(row_length,rows,BL_LEVELS)                                   &
                                          ! IN Cloud ice (kg per kg air)
     &,QCL(row_length,rows,BL_LEVELS)                                   &
                                          ! IN Cloud liquid water
     &,Q(row_length,rows,BL_LEVELS)                                     &
                                          ! IN specific humidity
     &,T(row_length,rows,BL_LEVELS)                                     &
                                          ! IN temperature
     &,CCA(row_length,rows)         ! IN Convective Cloud Amount
!                                     (decimal)

      INTEGER                                                           &
     & CCB(row_length,rows)                                             &
                                    ! IN Convective Cloud Base
     &,CCT(row_length,rows)         ! IN Convective Cloud Top

! (f) Atmospheric + any other data not covered so far, incl control.

      REAL                                                              &
     & CO2_MMR                                                          &
                                   ! IN CO2 Mass Mixing Ratio
     &,PHOTOSYNTH_ACT_RAD(ROW_LENGTH,ROWS)                              &
!                                  ! IN Net downward shortwave radiation
!                                  !    in band 1 (w/m2).
     &,PSTAR(row_length,rows)                                           &
                                    ! IN Surface pressure (Pascals).
     &,RAD_SICE(ROW_LENGTH,ROWS)                                        &
                                   ! IN Surface net shortwave and
!                                  !    downward LWradiation for
!                                  !    sea-ice (W/sq m).
     &,RAD_HR(row_length,rows,BL_LEVELS,2)                              &
                                    ! IN (LW,SW) rad heating rate (K/s)
     &, MICRO_TENDS(row_length, rows, bl_levels, 2)                     &
!                          ! Tendencies from microphys within BL levels
!                          ! (TL, K/s; QW, kg/kg/s)
     &,Muw_SBL,Mwt_SBL                                                  &
                                    ! IN Powers to use in prescription
!                                   !    of equilibrium profiles of
!                                   !    stress and buoyancy flux in
!                                   !    Equilibrium SBL model
     & ,Charnock                                                        &
                   ! Charnock parameter for sea surface
     &, SOIL_CLAY ( ROW_LENGTH, ROWS )                                  &
     &, SOIL_SILT ( ROW_LENGTH, ROWS )                                  &
     &, SOIL_SAND ( ROW_LENGTH, ROWS )                                  &
     &, DUST_MREL1 ( ROW_LENGTH, ROWS )                                 &
     &, DUST_MREL2 ( ROW_LENGTH, ROWS )                                 &
     &, DUST_MREL3 ( ROW_LENGTH, ROWS )                                 &
     &, DUST_MREL4 ( ROW_LENGTH, ROWS )                                 &
     &, DUST_MREL5 ( ROW_LENGTH, ROWS )                                 &
     &, DUST_MREL6 ( ROW_LENGTH, ROWS )                                 &

     &,TIMESTEP                                                         &
                                    ! IN Timestep (seconds).
     &,flux_e(row_length,rows)                                          &
                                   ! IN Surf. lat. heat flux   (W/m^2)
     &,flux_h(row_length,rows)     ! IN Surf. sens. heat flux  (W/m^2)
      LOGICAL                                                           &
     & L_flux_bc                                                        &
                                   ! IN Switches for prescribed surface
     &,L_spec_z0                   !    fluxes and roughness lengths

      LOGICAL                                                           &
     & L_SBLeq                                                          &
                                 ! IN Switch for Equilibrium SBL model
     &,L_SBLco                                                          &
                                 ! IN Switch for coupled gradient
!                                !    method in Equilibrium SBL model
     &,L_DUST                                                           &  
                                 ! IN Switch for mineral dust
     &,L_CAM_DUST                ! IN Use old version of dust_uplift 
                                 ! scheme for use in CAM NWP models

! Additional MOSES II variables

#include "nstypes.h"

      INTEGER                                                           &
     & LAND_PTS_TRIF                                                    &
                                   ! IN For dimensioning land fields
     &,NPFT_TRIF                                                        &
                                   ! IN For dimensioning PFT fields
!                                  !    available only with TRIFFID.
!                                  !    Set to NPFT when TRIFFID on,
!                                  !    set to 1 when TRIFFID off.
     &,CO2_DIM_LEN                                                      &
                                   ! IN Length of a CO2 field row.
     &,CO2_DIM_ROW                                                      &
                                   ! IN Number of CO2 field rows.
     &,ASTEPS_SINCE_TRIFFID                                             &
                                   ! IN Number of atmospheric
!                                  !    timesteps since last call
!                                  !    to TRIFFID.
     &,CAN_MODEL                   ! IN Swith for thermal vegetation
!                                  !    canopy

      REAL                                                              &
     & SNOW_TILE(LAND_PTS,NTILES)                                       &
                                   ! IN Lying snow on tiles (kg/m2)
     &,Z0_TILE(LAND_PTS,NTILES)                                         &
                                   ! IN Tile roughness lengths (m).
     &,LW_DOWN(ROW_LENGTH,ROWS)                                         &
                                   ! IN Surface downward LW radiation
!                                  !    (W/m2).
     &,SW_DOWN(ROW_LENGTH,ROWS)                                         &
                                   ! Surface downward SW radiation (W/m2).
     &,SW_TILE(LAND_PTS,NTILES)                                         &
                                   ! IN Surface net SW radiation on
!                                  !    land tiles (W/m2).
     &,CO2_3D(CO2_DIM_LEN,CO2_DIM_ROW)                                  &
!                                  ! IN 3D CO2 field if required.
     &,CS(LAND_PTS,DIM_CS1)                                             &
                              ! IN Soil carbon (kg C/m2).
     &,FRAC(LAND_PTS,NTYPE)                                             &
                                   ! IN Fractions of surface types.
     &,FLAND(LAND_PTS)                                                  &
                                   ! IN Land fraction on land tiles.
     &,FLANDG(ROW_LENGTH,ROWS)                                          &
!                                  ! IN Land fraction on all tiles.
     &,CANHT_FT(LAND_PTS,NPFT)                                          &
                                   ! IN Canopy height (m)
     &,LAI_FT(LAND_PTS,NPFT)                                            &
                                   ! IN Leaf area index
     &,ANTHROP_HEAT(NTILES)        ! IN Additional heat source on tiles 
!                                  !   for anthropogenic urban heat 
!                                  !   source (W/m2)  

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
     &, CPOOL_TILE(LAND_PTS,NTILES,10)                           &
     &, NPOOL_TILE(LAND_PTS,NTILES,10)                           &
     &, PPOOL_TILE(LAND_PTS,NTILES,12)                           &
     &, SOIL_ORDER(LAND_PTS)                                     &
     &, NIDEP(land_pts)                                          &
     &, NIFIX(land_pts)                                          &
     &, PWEA(land_pts)                                           &
     &, PDUST(land_pts)                                          &
     &, GLAI(land_pts,NTILES)                                    &
     &, PHENPHASE(land_pts,NTILES)                               &

     &, surf_down_sw(row_length,rows,4)                                 &
     &, ls_rain(row_length, rows)                                       &
     &, ls_snow(row_length, rows)                                       &
     &, conv_rain(row_length, rows)                                     &
     &, conv_snow(row_length, rows)                                     &
     &, lat(row_length, rows)                                           &
     &, long(row_length, rows)                                          &
     &, F_ROOT(sm_levels)                                               &
     &, time_sec                                                        

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
!sxy
     &,CH_CAB(LAND_PTS)                                              &
     &,CD_CAB(LAND_PTS)                                              &
     &,U_S_CAB(LAND_PTS)

      Integer                                                           &
     &  day                                                             &
     &, total_nsteps                                                    &
                                ! Total number of steps in run          
     &, SOIL_TYPE(row_length,rows)                                      &
     &, VEG_TYPE(row_length,rows)                                       &
     &, ISNOW_FLG3L(LAND_PTS,NTILES)

      Logical l_cable
      LOGICAL, DIMENSION(LAND_PTS,NTILES) :: L_TILE_PTS

!-----------------------------------------------------------

      LOGICAL                                                           &
     & L_CO2_INTERACTIVE                                                &
                                   ! IN Switch for 3D CO2 field
     &,L_PHENOL                                                         &
                                   ! IN Indicates whether phenology
!                                  !    in use
     &,L_TRIFFID                                                        &
                                   ! IN Indicates whether TRIFFID
!                                  !    in use.
     &,L_Q10                       ! IN Indicates Q10 for soil resp'n


      LOGICAL LTIMER                ! Logical switch for TIMER diags
      Logical :: L_ukca   ! switch for UKCA scheme

!  STASH flags :-

      LOGICAL                                                           &
     & SFME                                                             &
               ! IN Flag for FME (q.v.).
     &,SZ0HEFF                                                          &
               ! IN Flag for Z0H_EFF
     &,SIMLT                                                            &
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
     &,SV10    ! IN Flag for V10M (q.v.)
!
! Additional variables for SCM diagnostics which are dummy in full UM
      INTEGER                                                           &
     & nSCMDpkgs             ! No of SCM diagnostics packages

      LOGICAL                                                           &
     & L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

!  In/outs :-

      REAL                                                              &
     & GS(LAND_PTS)                                                     &
                                   ! INOUT "Stomatal" conductance to
!                                  !        evaporation (m/s).
     &,T_SOIL(land_pts,SM_LEVELS)                                       &
                                   ! INOUT Soil temperatures (K).
     &,TI(row_length,rows)                                              &
                                   ! INOUT Sea-ice surface layer
     &,TSTAR_LAND(ROW_LENGTH,ROWS)                                      &
                                   ! OUT   Land mean sfc temperature (K)
     &,TSTAR_SEA(ROW_LENGTH,ROWS)                                       &
                                   ! IN    Open sea sfc temperature (K).
     &,TSTAR_SICE(ROW_LENGTH,ROWS)                                      &
                                   ! INOUT Sea-ice sfc temperature (K).
     &,TSTAR_SSI(ROW_LENGTH,ROWS)                                       &
                                   ! INOUT Sea mean sfc temperature (K).
!                                      temperature (K).
     &,TSTAR(row_length,rows)                                           &
                                   ! INOUT Surface temperature (K).
     &,TSTAR_TILE(LAND_PTS,NTILES)                                      &
                                   ! INOUT Surface tile temperatures
     &,Z_LAND(ROW_LENGTH,ROWS)                                          &
                                   ! IN    Land height (m).
     &,Z0MSEA(row_length,rows)                                          
                                   ! INOUT Sea-surface roughness
!                                      length for momentum (m).
!                                      NB: same storage is used
!                                      for Z0V, so the intent is
!                                      IN for land points.
      REAL                                                              &
     & ALBSOIL(LAND_PTS)                                                &
!                                  ! Soil albedo.
     &, COS_ZENITH_ANGLE(row_length, rows)
!                                  ! Cosine of the zenith angle

      INTEGER                                                           &
     & CAN_RAD_MOD                                                      &
                                   !Switch for canopy radiation model
     & ,ILAYERS                                                         &
                                   !No of layers in canopy radiation model

     &,Z0M_SCM(row_length,rows)                                         &
                                   ! Fixed Sea-srf. roughness length(m)
                                   ! for momentum (SCM namelist)
     &,Z0H_SCM(row_length,rows)    ! Fixed Sea-srf. roughness length(m)
                                   ! for heat (SCM namelist)

      LOGICAL CUMULUS(row_length,rows)                                  &
                                   ! INOUT Logical switch for trade Cu
     &,L_SHALLOW(row_length,rows)  ! INOUT Flag to indicate shallow
!                                  !     convection

      INTEGER                                                           &
     & NTML(row_length,rows)                                            &
                                 ! INOUT Number of model levels in the
!                                !     surface-based turbulently mixed
!                                !     layer.
     &,NTPAR(row_length,rows)                                           &
                                 ! INOUT Top level of initial parcel
!                                !  ascent. Used in convection scheme.
     &,NLCL(row_length,rows)     ! INOUT No of levels to LCL
!
      REAL                                                              &
     & DELTHVU(row_length,rows)                                         &
                                 ! INOUT Integral of undilute parcel
!                                !     buoyancy over convective cloud
!                                !     layer (for convection scheme)
     &,ZHPAR(row_length,rows)                                           &
                                 ! INOUT Height of top of initial
!                                !     parcel ascent
     &,Z_LCL(row_length,rows)    ! INOUT Height of LCL

! Additional MOSES II variables
      REAL                                                              &
     & G_LEAF_ACC(LAND_PTS,NPFT)                                        &
                                   ! INOUT Accumulated G_LEAF
! Lestevens 17may13 - change npft to ntiles
     &,NPP_FT_ACC(LAND_PTS_TRIF,NTILES)                                 &
!     &,NPP_FT_ACC(LAND_PTS_TRIF,NPFT_TRIF)                              &
!                                  ! INOUT Accumulated NPP_FT
     &,RESP_W_FT_ACC(LAND_PTS_TRIF,NTILES)                              &
!     &,RESP_W_FT_ACC(LAND_PTS_TRIF,NPFT_TRIF)                           &
!                                  ! INOUT Accum RESP_W_FT
     &,RESP_S_ACC(LAND_PTS_TRIF,DIM_CS1) ! INOUT Accumulated RESP_S


! Definition of variables for STPH_RP
      REAL,INTENT(InOut) :: par_mezcla ! Used to modify
                                       ! LAMBDAH,LAMBDAM in EXCOEF
      REAL,INTENT(InOut) :: G0_RP ! Stability function param in EXCOEF

!  Outputs :-
!-1 Diagnostic (or effectively so - includes coupled model requisites):-

!  (a) Calculated anyway (use STASH space from higher level) :-
!
      REAL                                                              &
     & CD(row_length,rows)                                              &
                                   ! OUT Turbulent surface exchange
!                                     (bulk transfer) coefficient for
!                                     momentum.
     &,CH(row_length,rows)                                              &
                                   ! OUT Turbulent surface exchange
!                                     (bulk transfer) coefficient for
!                                     heat and/or moisture.
     &,E_SEA(row_length,rows)                                           &
                                   ! OUT Evaporation from sea times
!                                     leads fraction. Zero over land.
!                                     (kg per square metre per sec).
     &,FQW(row_length,rows,BL_LEVELS)                                   &
                                   ! OUT Moisture flux between layers
!                                     (kg per square metre per sec).
!                                     FQW(,1) is total water flux
!                                     from surface, 'E'.
     &,FQW_TILE(LAND_PTS,NTILES)                                        &
                                   ! OUT Surface FQW for land tiles
     &,EPOT_TILE(LAND_PTS,NTILES)                                       &
                                   ! OUT Local EPOT for land tiles.
     &,FTL(row_length,rows,BL_LEVELS)                                   &
                                   ! OUT FTL(,K) contains net turbulent
!                                     sensible heat flux into layer K
!                                     from below; so FTL(,1) is the
!                                     surface sensible heat, H. (W/m2)
     &,FTL_TILE(LAND_PTS,NTILES)                                        &
                                   ! OUT Surface FTL for land tiles
     &,H_SEA(row_length,rows)                                           &
                                   ! OUT Surface sensible heat flux over
!                                     sea times leads fraction. (W/m2)
     &,RHOKH(row_length,rows,BL_LEVELS)                                 &
                                   ! OUT Exchange coeffs for moisture.
     &,RHOKM(1-off_x:row_length+off_x,1-off_y:rows+off_y,BL_LEVELS)     &
!                                  ! OUT Exchange coefficients for
!                                  !     momentum on P-grid
     &,RHOKM_U(row_length,rows,BL_LEVELS)                               &
                                   ! OUT Exchange coefficients for u
     &,RHOKM_V(row_length,n_rows,BL_LEVELS)                             &
                                   ! OUT Exchange coefficients for v
     &,RIB_GB(row_length,rows)                                          &
                                   ! OUT Mean bulk Richardson number for
!                                     lowest layer.
     &,RIB_LAND(row_length,rows)                                        &
                                   !     Land mean bulk Richardson no.
!                                        for lowest layer.
     &,RIB_SSI(row_length,rows)                                         &
                                   ! OUT Sea mean bulk Richardson no.
!                                        for lowest layer.
     &,TAUX(row_length,rows,BL_LEVELS)                                  &
                                   ! OUT W'ly component of surface wind
!                                     stress (N/sq m).(On UV-grid with
     &,TAUX_LAND(ROW_LENGTH,ROWS)                                       &
                                   ! OUT W'ly compt of land sfc wind
!                                  !     stress (N/sq m). (On U-grid
!                                  !     with first and last rows
!                                  !     undefined or, at present,
!                                  !     set to missing data
     &,TAUX_SSI(ROW_LENGTH,ROWS)                                        &
                                   ! OUT W'ly compt of mean sea sfc wind
!                                  !     stress (N/sq m). (On U-grid
!                                  !     with first and last rows
!                                  !     undefined or, at present,
!                                  !     set to missing data
!                                     first and last rows undefined or
!                                     at present, set to missing data
     &,TAUY(row_length,n_rows,BL_LEVELS)                                &
                                   ! OUT S'ly component of surface wind
     &,TAUY_LAND(ROW_LENGTH,N_ROWS)                                     &
                                   ! OUT S'ly compt of land sfc wind
!                                  !     stress (N/sq m).  On V-grid;
!                                  !     comments as per TAUX.
     &,TAUY_SSI(ROW_LENGTH,N_ROWS)                                      &
                                   ! OUT S'ly compt of mean sea sfc wind
!                                  !     stress (N/sq m).  On V-grid;
!                                  !     comments as per TAUX.
!                                     stress (N/sq m).  On UV-grid;
!                                     comments as per TAUX.
     &,VSHR(row_length,rows)                                            &
                                   ! OUT Magnitude of surface-to-lowest
!                                     atm level wind shear (m per s).
     &,VSHR_LAND(row_length,rows)                                       &
                                   ! OUT Magnitude of land sfc-to-lowest
!                                     atm level wind shear (m per s).
     &,VSHR_SSI(row_length,rows)                                        &
                                   ! OUT Mag. of mean sea sfc-to-lowest
!                                     atm level wind shear (m per s).
     &,ZHT(row_length, rows)                                            &
                                   ! OUT Max height of turb mixing
     &,RHO_CD_MODV1(row_length,rows)                                    &
                                   ! OUT Surface air density * drag coef
!                                  !  *mod(v1 - v0) before interpolation
     &,BL_TYPE_1(row_length,rows)                                       &
                                 ! OUT Indicator set to 1.0 if stable
!                                  !     b.l. diagnosed, 0.0 otherwise.
     &,BL_TYPE_2(row_length,rows)                                       &
                                 ! OUT Indicator set to 1.0 if Sc over
!                                  !     stable surface layer diagnosed,
!                                  !     0.0 otherwise.
     &,BL_TYPE_3(row_length,rows)                                       &
                                 ! OUT Indicator set to 1.0 if well
!                                  !     mixed b.l. diagnosed,
!                                  !     0.0 otherwise.
     &,BL_TYPE_4(row_length,rows)                                       &
                                 ! OUT Indicator set to 1.0 if
!                                  !     decoupled Sc layer (not over
!                                  !     cumulus) diagnosed,
!                                  !     0.0 otherwise.
     &,BL_TYPE_5(row_length,rows)                                       &
                                 ! OUT Indicator set to 1.0 if
!                                  !     decoupled Sc layer over cumulus
!                                  !     diagnosed, 0.0 otherwise.
     &,BL_TYPE_6(row_length,rows)                                       &
                                 ! OUT Indicator set to 1.0 if a
!                                  !     cumulus capped b.l. diagnosed,
!                                  !     0.0 otherwise.
     &,BL_TYPE_7(row_length,rows)                                       &
                                 ! OUT Indicator set to 1.0 if a
!                                  !     Shear-dominated unstable b.l.
!                                  !     diagnosed, 0.0 otherwise.
     &,RHO_KM(row_length,rows,2:BL_LEVELS)                              &
                                   ! OUT Air density * turbulent mixing
!                                     coefficient for momentum before
!                                     interpolation.
     &,RHO_ARESIST(row_length,rows)                                     &
                                   ! OUT RHOSTAR*CD_STD*VSHR for SULPHUR
!                                     cycle
     &,ARESIST(row_length,rows)                                         &
                                   ! OUT 1/(CD_STD*VSHR) for Sulphur
!                                     cycle
     &,RESIST_B(row_length,rows)                                        &
                                   ! OUT (1/CH-1/(CD_STD)/VSHR for
!                                     Sulphur cycle
     &, WE_LIM(row_length,rows,3)                                       &
                                    ! OUT rho*entrainment rate implied b
!                                   !     placing of subsidence
     &, ZRZI(row_length,rows,3)                                         &
                                    ! OUT (z-z_base)/(z_i-z_base)
     &, T_FRAC(row_length,rows,3)                                       &
                                    ! OUT a fraction of the timestep
     &, WE_LIM_DSC(row_length,rows,3)                                   &
!                                   ! OUT rho*entrainment rate implied b
!                                   !     placing of subsidence
     &, ZRZI_DSC(row_length,rows,3)                                     &
                                    ! OUT (z-z_base)/(z_i-z_base)
     &, T_FRAC_DSC(row_length,rows,3)                                   &
!                                   ! OUT a fraction of the timestep
     &, Z_HALF(row_length,rows,BL_LEVELS)                               &
                                    ! OUT Z_HALF(*,K) is height of half
!                                   ! level k-1/2.
     &, ZHSC(row_length,rows)                                           &
                                    ! OUT Top of decoupled layer
     &, R_B_DUST(ROW_LENGTH,ROWS,NDIV)                                  &
                                      !OUT surface layer resist for dust
     &, WSTAR(row_length, rows)                                         &
                                   ! OUT Convective velocity scale (m/s)
     &, WTHVS(row_length, rows)                                         &
                                   ! OUT surface flux of thv (Km/s)
     &, SHALLOWC(row_length,rows)                                       &
                                   ! OUT Shallow Cu diagnostic
!                                  !   Indicator set to 1.0 if shallow,
!                                  !   0.0 if not shallow or not cumulus
     &, CU_OVER_OROG(row_length,rows)                                   &
!                                  ! OUT Indicator for cumulus
!                                  !     over steep orography
!                                  !   Indicator set to 1.0 if true,
!                                  !   0.0 if false. Exclusive.
     &,WT_EXT(LAND_PTS,sm_levels)                                       &
                                  !OUT cumulative fraction of transp'n
     &,RA(LAND_PTS)               !OUT Aerodynamic resistance (s/m)

      INTEGER                                                           &
     & NTDSC(row_length,rows)                                           &
                                   ! OUT Top level for turb mixing in
!                                           any decoupled Sc layer
     &,NBDSC(row_length,rows)                                           &
                                   ! OUT Bottom level of any decoupled
!                                  !     turbulently-mixed Sc layer.
     &, KENT(row_length,rows)                                           &
                                    ! OUT grid-level of SML inversion
     &, KENT_DSC(row_length,rows)   ! OUT grid-level of DSC inversion


!  (b) Not passed between lower-level routines (not in workspace at this
!      level) :-

      REAL                                                              &
     & FME(row_length,rows)      ! OUT Wind mixing "power" (W per sq m).
!
!-2 Genuinely output, needed by other atmospheric routines :-

      REAL                                                              &
     & ZH(row_length,rows)                                              &
                                 ! INOUT Height above surface of top of
!                                   boundary layer (metres).
     &,GPP(LAND_PTS)                                                    &
                                 ! OUT Gross primary productivity
!                                !     (kg C/m2/s).
     &,NPP(LAND_PTS)                                                    &
                                 ! OUT Net primary productivity
!                                !     (kg C/m2/s).
     &,RESP_P(LAND_PTS)                                                 &
                                 ! OUT Plant respiration (kg C/m2/s).
     &,T1_SD(row_length,rows)                                           &
                                 ! OUT Standard deviation of turbulent
!                                   fluctuations of layer 1 temperature;
!                                   for use in initiating convection.
     &,Q1_SD(row_length,rows)                                           &
                                 ! OUT Standard deviation of turbulent
!                                   fluctuations of layer 1 humidity;
!                                   for use in initiating convection.
     &,z0m_eff_gb(row_length,rows)                                      &
     &,z0h_eff_gb(row_length,rows)
                                 ! OUT Effective grid-box roughness
!                                   lengths for momentum and for 
!                                   heat, moisture
      Real                                                              &
     &  uw0(row_length,rows)                                            &
                             ! OUT U-component of surface wind stress
!                            !     on P-grid
     &, vw0(row_length,rows) ! OUT V-component of surface wind stress
!                            !     on P-grid

! OUT variables for IMP_SOLVER (that used to be local arrays)
      REAL                                                              &
     & ALPHA1(LAND_PTS,NTILES)                                          &
                                  ! OUT Mean gradient of saturated
!                                 !     specific humidity with respect
!                                 !     to temperature between the
!                                 !     bottom model layer and tile
!                                 !     surfaces
     &,ASHTF(row_length,rows)                                           &
                                  ! OUT Coefficient to calculate surface
!                                 heat flux into soil or sea-ice.
     &,DTRDZ_CHARNEY_GRID(row_length,rows,BL_LEVELS)                    &
!                                 ! OUT dt/(rho*r*r*dz) for scalar
!                                 !     flux divergence
     &,RDZ_CHARNEY_GRID(row_length,rows,BL_LEVELS)                      &
!                               ! OUT RDZ(,1) is the reciprocal of the
!                                 height of level 1, i.e. of the
!                                 middle of layer 1.  For K > 1,
!                                 RDZ(,K) is the reciprocal
!                                 of the vertical distance
!                                 from level K-1 to level K.
     &,DTRDZ_U(row_length,rows,BL_LEVELS)                               &
                                           ! OUT dt/(rho*r*r*dz) for
     &,DTRDZ_V(row_length,n_rows,BL_LEVELS)                             &
                                           ! OUT U,V flux divergence
     &,RDZ_U(row_length,rows,2:BL_LEVELS)                               &
                                ! OUT RDZ (K > 1) on UV-grid.
     &,RDZ_V(row_length,n_rows,2:BL_LEVELS)                             &
                                ! OUT RDZ (K > 1) on UV-grid.
     &,FRACA(LAND_PTS,NTILES)                                           &
                                ! OUT Fraction of surface moisture
!                               !     flux with only aerodynamic
!                               !     resistance for snow-free land
!                               !     tiles.
     &,RHOKH_TILE(LAND_PTS,NTILES)                                      &
!                               ! OUT Surface exchange coefficients
!                               !     for land tiles
     &,SMC(LAND_PTS)                                                    &
                                ! OUT Available moisture in the
!                               !     soil profile (mm).
     &,CHR1P5M(LAND_PTS,NTILES)                                         &
                                ! OUT Ratio of coefffs for
!                               !     calculation of 1.5m temp for
!                               !     land tiles.
     &,RESFS(LAND_PTS,NTILES)                                           &
                                ! OUT Combined soil, stomatal
!                               !     and aerodynamic resistance
!                               !     factor for fraction (1-FRACA)
!                               !     of snow-free land tiles.
     &,Z0HSSI(ROW_LENGTH,ROWS)                                          &
     &,Z0MSSI(ROW_LENGTH,ROWS)                                          &
                                ! OUT Roughness lengths over sea (m).
     &,CDR10M_U(row_length,rows)                                        &
                                ! OUT Ratio of CD's reqd for calculation
     &,CDR10M_V(row_length,n_rows)                                      &
                                ! OUT Ratio of CD's reqd for calculation
     &,Z1_TQ(row_length,rows)   ! OUT Height of lowest theta level.


! Additional MOSES II variables
      INTEGER                                                           &
     & TILE_INDEX(LAND_PTS,NTYPE)                                       &
                                   ! OUT Index of tile points
     &,TILE_PTS(NTYPE)             ! OUT Number of tile points


      REAL                                                              &
     & LE_TILE(LAND_PTS,NTILES)                                         &
                                   ! OUT Surface latent heat flux for
!                                  !     land tiles
     &,RADNET_SICE(ROW_LENGTH,ROWS)                                     &
                                   ! OUT Surface net radiation on
!                                  !     sea-ice (W/m2)
     &,RADNET_TILE(LAND_PTS,NTILES)                                     &
                                   ! OUT Surface net radiation on
!                                  !     land tiles (W/m2)
     &,RIB_TILE(LAND_PTS,NTILES)                                        &
                                   ! OUT RIB for land tiles.
     &,RHO_ARESIST_TILE(LAND_PTS,NTILES)                                &
!                                  ! OUT RHOSTAR*CD_STD*VSHR on land
!                                  !     tiles
     &,ARESIST_TILE(LAND_PTS,NTILES)                                    &
!                                  ! OUT 1/(CD_STD*VSHR) on land tiles
     &,RESIST_B_TILE(LAND_PTS,NTILES)                                   &
!                                  ! OUT (1/CH-1/CD_STD)/VSHR on land
!                                  !     tiles
     &,ALPHA1_SICE(ROW_LENGTH,ROWS)                                     &
                                   ! OUT ALPHA1 for sea-ice.
     &,ASHTF_TILE(LAND_PTS,NTILES)                                      &
                                   !OUT Coefficient to calculate
!                                  !     surface heat flux into land
!                                  !     tiles.
     &,FQW_ICE(ROW_LENGTH,ROWS)                                         &
                                   ! OUT Surface FQW for sea-ice
     &,FTL_ICE(ROW_LENGTH,ROWS)                                         &
                                   ! OUT Surface FTL for sea-ice
     &,RESFT(LAND_PTS,NTILES)                                           &
                                   ! OUT Total resistance factor.
!                                  !     FRACA+(1-FRACA)*RESFS for
!                                  !     snow-free land, 1 for snow.
     &,RHOKH_SICE(ROW_LENGTH,ROWS)                                      &
                                   ! OUT Surface exchange coefficients
!                                  !     for sea and sea-ice
     &,RHOKPM(LAND_PTS,NTILES)                                          &
                                   ! OUT Surface exchange coefficient.
     &,RHOKPM_POT(LAND_PTS,NTILES)                                      &
                                   ! OUT Potential evaporation
!                                  !     exchange coefficient.
     &,RHOKPM_SICE(ROW_LENGTH,ROWS)                                     &
                                   ! OUT Sea-ice surface exchange coeff.
     &,Z0H_TILE(LAND_PTS,NTILES)                                        &
                                   ! OUT Tile roughness lengths for heat
!                                  !     and moisture (m).
     &,Z0M_GB(ROW_LENGTH,ROWS)                                          &
                                   ! OUT Gridbox mean Roughness length
!                                  !      for momentum (m).
     &,Z0M_TILE(LAND_PTS,NTILES)                                        &
                                   ! OUT Tile roughness lengths for
!                                  !     momentum.
     &,CHR1P5M_SICE(ROW_LENGTH,ROWS)                                    &
!                                  ! OUT CHR1P5M for sea and sea-ice
!                                  !     (leads ignored).
     !kdcorbin, 11/10 - changed from NPFT
     &,G_LEAF(LAND_PTS,NTILES)                                          &
                                   ! OUT Leaf turnover rate (/360days).
     !kdcorbin, 11/10 - changed from NPFT
     &,GPP_FT(LAND_PTS,NTILES)                                          &
                                   ! OUT Gross primary productivity
!                                  !     on PFTs (kg C/m2/s).
     !kdcorbin, 11/10 - changed from NPFT
     &,NPP_FT(LAND_PTS,NTILES)                                          &
                                   ! OUT Net primary productivity
!                                  !     (kg C/m2/s).
     !kdcorbin, 11/10 - changed from NPFT
     &,RESP_P_FT(LAND_PTS,NTILES)                                       &
                                   ! OUT Plant respiration on PFTs
!                                  !     (kg C/m2/s).
     &,RESP_S(LAND_PTS,DIM_CS1)                                         &
                               ! OUT Soil respiration (kg C/m2/s).
     &,RESP_S_TOT(DIM_CS2)                                              &
                              ! OUT Total soil respiration
                              ! (kg C/m2/s).
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
     &,FLANDG_U(ROW_LENGTH,ROWS)                                        &
                                   ! OUT Land frac (on U-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
     &,FLANDG_V(ROW_LENGTH,N_ROWS)                                      &
                                   ! OUT Land frac (on V-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
     &,TILE_FRAC(LAND_PTS,NTILES)                                       &
                                   ! OUT Tile fractions including
!                                  !     snow cover in the ice tile.
     &,FSMC(LAND_PTS,NPFT)         ! OUT Moisture availability factor.



      INTEGER                                                           &
     & ERROR                    ! OUT 0 - AOK;
!                               !     1 to 7  - bad grid definition
                                !     detected;

!---------------------------------------------------------------------
!  External routines called :-

      EXTERNAL BDY_EXPL1,SF_EXPL,BDY_EXPL2
      EXTERNAL DUST_SRCE, DUST_SRCE_CAM
      EXTERNAL TIMER

!-----------------------------------------------------------------------

!  Workspace :-

       INTEGER                                                          &
     & IDIV                                                             &
            ! loop counter, mineral dust divisions
     &,M !loop counter

      INTEGER ISeaZ0T              ! Switch for the definition of
!                                  ! the thermal roughness
!                                  ! length over the sea.
      INTEGER COR_UST              ! Switch to use correct 
!                                  ! friction velocity
      INTEGER COR_MO_ITER          ! Switch for MO iteration correction
      INTEGER Buddy_sea            ! Switch to use the wind speed from 
!                                  ! adjacent sea points for the 
!                                  ! sea part of coastal grid points
      INTEGER :: I_SCRN_T_DIAG
!                                  ! Option for the diagnosis of
!                                  ! screen temperature
      INTEGER FD_stab_dep          ! Switch to implement stability
!                                  ! dependence of orographic form drag

      REAL                                                              &
     & QW(row_length, rows, BL_LEVELS)                                  &
     &,TL(row_length, rows, BL_LEVELS)                                  &
     &,A_DQSDT(row_length,rows,BL_LEVELS)                               &
!                               ! Saturated lapse rate factor
!                               ! on p,T,q-levels (full levels).
     &,A_QS(row_length,rows,BL_LEVELS)                                  &
                                ! Saturated lapse rate factor
!                               ! on p,T,q-levels (full levels).
     &,BQ(row_length,rows,BL_LEVELS)                                    &
                                ! A buoyancy parameter for clear air
!                               ! on p,T,q-levels (full levels).
     &,BQ_CLD(row_length,rows,BL_LEVELS)                                &
                                ! A buoyancy parameter for cloudy air
!                               ! on p,T,q-levels (full levels).
     &,BQ_GB(row_length,rows,BL_LEVELS)                                 &
                                ! A grid-box mean buoyancy parameter
!                               ! on p,T,q-levels (full levels).
     &,BT(row_length,rows,BL_LEVELS)                                    &
                                ! A buoyancy parameter for clear air
!                               ! on p,T,q-levels (full levels).
     &,BT_CLD(row_length,rows,BL_LEVELS)                                &
!                               ! A buoyancy parameter for cloudy air
!                               ! on p,T,q-levels (full levels).
     &,BT_GB(row_length,rows,BL_LEVELS)                                 &
                                ! A grid-box mean buoyancy parameter
!                               ! on p,T,q-levels (full levels).
     &,DELTAP(row_length,rows,BL_LEVELS)                                &
                                ! Difference in pressure between levels
     &,DQSDT(row_length,rows,BL_LEVELS)                                 &
                                ! Derivative of q_SAT w.r.t. T
     &,DZL_charney(row_length,rows,BL_LEVELS)                           &
                                ! DZL(,K) is depth in m of theta level
!                                 K, i.e. distance from boundary
!                                 K-1/2 to boundary K+1/2.
     &,FB_SURF(row_length,rows)                                         &
                                ! Surface flux buoyancy over density
!                               ! (m^2/s^3)
     &,P_HALF(row_length,rows,BL_LEVELS)                                &
                                ! P_HALF(*,K) is pressure at half
!                               ! level k-1/2.
     &,Z_FULL(row_length,rows,BL_LEVELS)                                &
                                ! Z_FULL(*,K) is height of full level k.
     &,Z_UV(row_length,rows,BL_LEVELS)                                  &
                                ! Z_UV(*,K) is height of half level
!                               ! k-1/2.
     &,Z_TQ(row_length,rows,BL_LEVELS)                                  &
                                ! Z_TQ(*,K) is height of half level
!                               ! k+1/2.
     &,RDZ(row_length,rows,BL_LEVELS)                                   &
                                ! RDZ(,1) is the reciprocal of the
!                                 height of level 1, i.e. of the
!                                 middle of layer 1.  For K > 1,
!                                 RDZ(,K) is the reciprocal
!                                 of the vertical distance
!                                 from level K-1 to level K.
     &,RHO_UV(row_length,rows,BL_LEVELS+1)                              &
!                               ! OUT density on UV (ie. rho) levels;
!                               !    used in RHOKH so dry density if
!                               !    Lq_mix_bl is true
     &,RHO_TQ(row_length,rows,BL_LEVELS)                                &
!                               ! OUT density on TQ (ie. theta) levels;
!                               !    used in RHOKM so wet density
     &,RHO_DRY_TQ(row_length,rows,BL_LEVELS)                            &
!                               ! OUT density on TQ (ie. theta) levels;
!                               !    used in non-turb flux integration
!                               !    so dry density if Lq_mix_bl is true
     &,U_P(row_length,rows,BL_LEVELS)                                   &
                                        ! U on P-grid.
     &,U_S(row_length,rows)                                             &
                                ! Surface friction velocity (m/s)
     &,V_P(row_length,rows,BL_LEVELS)                                   &
                                        ! V on P-grid.
     &,RECIP_L_MO_SEA(row_length,rows)
!                                ! Reciprocal of the surface
!                                ! Obukhov length at sea points.
!                                ! (m-1).

       REAL                                                             &
     & TAUX_FD_U(row_length,rows,BL_LEVELS)                             &
!                                          ! X comp of orographic stress
!                                          ! interpolated to U points
     &,TAUY_FD_V(row_length,n_rows,BL_LEVELS),                          &
!                                          ! Y comp of orographic stress
!                                          ! interpolated to V points

     & Z1_UV(row_length,rows)                                           &
                                ! Height of lowest u,v level.
     &,H_BLEND_OROG(row_length,rows)                                    &
                                ! Blending height used as part of
!                                 effective roughness scheme
     &,RHOSTAR(ROW_LENGTH,ROWS)                                         &
                                ! Surface air density
     &,DUST_FLUX(ROW_LENGTH,ROWS,NDIV)                                  &
                                       !dust production flux(kg m-2 s-1)
     &,DUST_FLUX_TILE(LAND_PTS,NTILES,NDIV)                             &
                                            !production flux from tiles
     &,CD_STD_DUST(ROW_LENGTH,ROWS)                                     &
                                   !  Bulk transfer coef. for
!                             ! momentum, excluding orographic effects
     &,U_S_STD_TILE(LAND_PTS,NTILES)                                    &
                                     ! Surface friction velocity
     &,CLAY_LAND(LAND_PTS)                                              &
                           ! soil clay fraction on land pts
     &,SAND_LAND(LAND_PTS)                                              &
                           ! soil sand fraction on land pts
     &,PSTAR_LAND(LAND_PTS)                                             &
                            ! surface pressure on land pts
     &,RHOSTAR_LAND(LAND_PTS)                                           &
                              ! surface air density on land pts
     &,MREL_LAND(LAND_PTS,NDIV)                                         &
                                ! soil size fraction on land pts
     &,U_S_T_TILE(LAND_PTS,NTILES,NDIV)                                 &
                                        !threshold friction vel for dust
     &,U_S_T_DRY_TILE(LAND_PTS,NTILES,NDIV)                             &
                                            !threshold friction velocity
     &,U_S_T(ROW_LENGTH,ROWS,NDIV)                                      &
                                   !threshold friction vel. for dust
     &,U_S_T_DRY(ROW_LENGTH,ROWS,NDIV) !threshold friction velocity
!                             !for dust, excluding soil moisture effects

       REAL                                                             &
     & RHOLEM                                                           &
               ! surface density in LEM
     &,TV1_SD                                                           &
               ! virt T standard deviation (approx)
     &,W_M     ! velocity scale
!     Declaration of new BL diagnostics.
      Type (Strnewbldiag) :: BL_diag

      ! end step of experiment, step width, processor num
      integer :: endstep, timestep_number, mype

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('BDYLAYR ',3)
      ENDIF
      ERROR = 0

!-----------------------------------------------------------------------
!! 0. Verify grid/subset definitions.  Arakawa 'B' grid with P-rows at
!!    extremes is assumed.  Extreme-most P-rows are ignored; extreme-
!!    most UV-rows are used only for interpolation and are not updated.
!-----------------------------------------------------------------------

      IF ( BL_LEVELS <  1 .OR. ST_LEVELS <  1 .OR. SM_LEVELS <  1       &
     & .OR. ROWS <  1 ) THEN
        ERROR = 1
        GOTO 999
      ELSEIF ( land_pts >  rows*row_length ) THEN
        ERROR = 7
        GOTO 999
      ENDIF

      IF (L_flux_bc) THEN
!       ! For specified surface fluxes impose uniform TSTAR,
!       ! calculated in CONV_DIAG to be consistent with fluxes
        DO J=1, rows
        DO I=1, row_length
         TSTAR_LAND(I,J) = TSTAR(I,J)
         TSTAR_SEA(I,J)  = TSTAR(I,J)
         TSTAR_SICE(I,J) = TSTAR(I,J)
         TSTAR_SSI(I,J)  = TSTAR(I,J)
        ENDDO
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!! Extract local switches from BL_OPTIONS array
!-----------------------------------------------------------------------
      ISeaZ0T   = BL_OPTIONS(8)
      COR_UST   = BL_OPTIONS(10)
      Buddy_Sea = BL_OPTIONS(14)
      I_SCRN_T_DIAG = BL_OPTIONS(15)
      COR_MO_ITER = BL_OPTIONS(16)
      FD_stab_dep = BL_OPTIONS(19)

! DEPENDS ON: bdy_expl1
      CALL BDY_EXPL1 (                                                  &

! IN MPP variables
     & halo_i, halo_j, off_x, off_y, row_length, rows, n_rows,          &
     & global_row_length, proc_row_group, at_extremity,                 &

!IN  Substepping information
     & Substep_Number,                                                  &
! IN values defining vertical grid of model atmosphere :
     & MODEL_DOMAIN,                                                    &
     & BL_LEVELS,                                                       &
     & r_rho_levels,r_theta_levels,                                     &
     & P,P_theta_levels, rho_rsq, rho_wet, rho_dry,                     &
     & SIN_THETA_LONGITUDE,COS_THETA_LONGITUDE,                         &

! IN U and V momentum fields.
     & U, V,                                                            &

! IN cloud data :
     & CF,Q,QCF,QCL,T,                                                  &

! IN everything not covered so far :
     & PSTAR,TIMESTEP,LQ_MIX_BL,                                        &

! OUT
     & DTRDZ_CHARNEY_GRID,RDZ_CHARNEY_GRID,DTRDZ_U,DTRDZ_V,             &
     & RDZ_U,RDZ_V,RHO_UV,RHO_TQ,RHO_DRY_TQ,DZL_charney,RDZ,            &
     & Z1_TQ,Z1_UV,Z_FULL,Z_HALF,Z_UV,Z_TQ,                             &
     & P_HALF,DELTAP,U_P,V_P,QW,TL,                                     &
     & BT,BQ,BT_CLD,BQ_CLD,BT_GB,BQ_GB,A_QS,A_DQSDT,DQSDT,              &

     & LTIMER                                                           &
     &  )



         DO L = 1,LAND_PTS
           J = (LAND_INDEX(L)-1)/ROW_LENGTH + 1
           I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
           CLAY_LAND(L) = SOIL_CLAY(I,J)
! NOTE!! NASTY HACK FOR NOW UNTIL ANCILLARY FILE IS SORTED OUT...
           IF (CLAY_LAND(L)  <   -1e9) CLAY_LAND(L)=0.23
         ENDDO !LAND_PTS

! DEPENDS ON: sf_expl
      CALL SF_EXPL (                                                    &

! IN values defining field dimensions and subset to be processed :
     & HALO_I,HALO_J,OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,                &
     & LAND_PTS,LAND_PTS_TRIF,NPFT_TRIF,                                &
     & DIM_CS1, DIM_CS2,                                                &

! IN  parameters for iterative SISL scheme
     & NumCycles, CycleNo,                                              &

! IN parameters required from boundary-layer scheme :
     & BQ_GB,BT_GB,Z1_UV,Z1_TQ,QW,TL,                                   &

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
     & T_SOIL,TI,TSTAR,TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,TSTAR_SSI,       &
     & TSTAR_TILE,Z_LAND,L_CTILE,COR_UST,COR_MO_ITER,I_SCRN_T_DIAG,     &
     & ALBSOIL,COS_ZENITH_ANGLE,CAN_RAD_MOD,ILAYERS,                    &
     & U,V,U_P,V_P,                                                     &
     & L_DUST,ISeaZ0T,Buddy_Sea, ANTHROP_HEAT,                          &

! IN STASH flags :-
     & SFME,SQ1P5,ST1P5,SU10,SV10,SZ0HEFF,                              &

! INOUT data :
     & Z0MSEA,L_spec_z0,Z0M_SCM,Z0H_SCM,GS,                             &
     & G_LEAF_ACC,NPP_FT_ACC,RESP_W_FT_ACC,RESP_S_ACC,                  &

! OUT Diagnostic not requiring STASH flags :
     & CD,CH,RECIP_L_MO_SEA,E_SEA,FQW,                                  &
     & FTL,FTL_TILE,LE_TILE,H_SEA,RADNET_SICE,RADNET_TILE,              &
     & RHOKM,RHOKM_U,RHOKM_V,RIB_GB,RIB_TILE,                           &
     & TAUX,TAUY,                                                       &
     & TAUX_LAND,TAUX_SSI,TAUY_LAND,TAUY_SSI,                           &

! OUT diagnostic requiring STASH flags :
     & FME,                                                             &

! OUT diagnostics required for soil moisture nudging scheme :
     & WT_EXT,RA,                                                       &

! OUT data required for tracer mixing :
     & RHO_ARESIST,ARESIST,RESIST_B,                                    &
     & RHO_ARESIST_TILE,ARESIST_TILE,RESIST_B_TILE,                     &
     & R_B_DUST,CD_STD_DUST,U_S_STD_TILE,                               &

! OUT data required for 4D-VAR :
     & RHO_CD_MODV1,                                                    &

! OUT data required elsewhere in UM system :
     & FB_SURF,U_S,T1_SD,Q1_SD,                                         &

! OUT data required elsewhere in boundary layer or surface code
     & ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,FQW_TILE,EPOT_TILE,          &
     & FQW_ICE,FTL_ICE,FRACA,RHOSTAR,RESFS,RESFT,                       &
     & RHOKH,RHOKH_TILE,RHOKH_SICE,RHOKPM,RHOKPM_POT,RHOKPM_SICE,       &
     & H_BLEND_OROG,Z0HSSI,Z0H_TILE,Z0H_EFF_GB,Z0M_GB,Z0MSSI,Z0M_TILE,  &
     & Z0M_EFF_GB,CDR10M_U,CDR10M_V,CHR1P5M,CHR1P5M_SICE,SMC,           &
     & VSHR,VSHR_LAND,VSHR_SSI,                                         &
     & GPP,NPP,RESP_P,G_LEAF,GPP_FT,NPP_FT,                             &
     & RESP_P_FT,RESP_S,RESP_S_TOT,CLAY_LAND,RESP_W_FT,                 &
     & GC,CANHC_TILE,WT_EXT_TILE,FLAKE,                                 &
     & TILE_INDEX,TILE_PTS,TILE_FRAC,FSMC,                              &
     & FLANDG_U,FLANDG_V,                                               &

! LOGICAL LTIMER
     & LTIMER,BL_diag                                                   &
     &, L_ukca                                                          &
! EAK
!    IN
     & ,l_cable,                                                        &
     & surf_down_sw,alb_tile,l_tile_pts,                                &
     & ls_rain,ls_snow,SW_DOWN,                                         &
     & lat,long,day,time_sec,                                           &
     & SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L,                   &
     & SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE,             &
     & TSOIL_TILE,T_SURF_TILE,HCONS,BEXP,                               &
     &  SATHH,SATCON,HCAP,                                              &
     & SOIL_TYPE,VEG_TYPE,                                              &
     & ISNOW_FLG3L,total_nsteps,                                        &
     &            FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB,              &
     &            TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB,          &
     &            USTAR_CAB,SURF_HTF_CAB,                               &
     &            U_S_CAB,CH_CAB,CD_CAB,                                &
     & SOIL_LAYER_MOISTURE,                                             &
     & SNAGE_TILE,RTSOIL_TILE,                                          &
     & GFLUX_TILE,SGFLUX_TILE,                                          &
     & F_ROOT, sin_theta_latitude,                                      &
! Lestevens Sept2012: CasaCNP variables 
     & CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,SOIL_ORDER,                     &
     & NIDEP,NIFIX,PWEA,PDUST,GLAI,PHENPHASE,                           &
      ! end step of experiment, this step, step width, processor num
     & endstep, timestep_number, mype                                   &
     & )

!       print *,'bdy_layr 3 soil_type,veg_type froot',soil_type,veg_type 
!       print *,'bdy_layr FQW after sf_xepl',fqw,FQW_TILE,ftl,ftl_tile
!        print *,'bdy_layr after sf_xepl snow_tile', SNOW_TILE

      conv_rain = 0. ! convective rain and snow calculated after call to ni_bl_ctl
      conv_snow = 0.

      IF (L_flux_bc) THEN
        ! Surface fluxes calculated in SFEXPL are substituted with
        ! forcing values.  NOTE: Surface calculation also made
        ! explicit (time weight set to zero).
        DO I=1, row_length
        DO J=1, rows
         !..Converts Fluxes from W/m^2 to rho*K/s
         RHOLEM = RHOSTAR(I,J)

         !..If comparing against LES with rho ne rhostar then match
         ! w'theta' rather than rho*wtheta

         ! RHOLEM = 1.0

         FQW(I,J,1)   = (RHOSTAR(I,J)*flux_e(I,J))/(LC*RHOLEM)
         FTL(I,J,1)   = (RHOSTAR(I,J)*flux_h(I,J))/(CP*RHOLEM)

         FB_SURF(I,J) = G * ( BT_GB(I,J,1)*FTL(I,J,1) +                 &
     &                        BQ_GB(I,J,1)*FQW(I,J,1) ) /RHOSTAR(I,J)
         IF ( FB_SURF(I,J)  >   0.0) THEN
          W_M        = ( 0.25*ZH(I,J)*FB_SURF(I,J) +                    &
     &                   U_S(I,J)*U_S(I,J)*U_S(I,J) ) ** (1.0/3.0)
          T1_SD(I,J) = 1.93 * FTL(I,J,1) / (RHOSTAR(I,J) * W_M)
          Q1_SD(I,J) = 1.93 * FQW(I,J,1) / (RHOSTAR(I,J) * W_M)
          TV1_SD     = T(I,J,1) * ( BT_GB(I,J,1)*T1_SD(I,J) +           &
     &                              BQ_GB(I,J,1)*Q1_SD(I,J) )
          T1_SD(I,J) = MAX ( 0.0 , T1_SD(I,J) )
          Q1_SD(I,J) = MAX ( 0.0 , Q1_SD(I,J) )
           IF (TV1_SD  <=  0.0) THEN
            T1_SD(I,J) = 0.0
            Q1_SD(I,J) = 0.0
           ENDIF
         ELSE
           T1_SD(I,J) = 0.0
           Q1_SD(I,J) = 0.0
         ENDIF
        ENDDO    ! J
       ENDDO    ! I
       DO I=1,land_pts
        DO L=1,NTILES
         FQW_TILE(I,L) = FQW(1,1,1)
         FTL_TILE(I,L) = FTL(1,1,1)
        ENDDO ! L
        ENDDO ! I
      ENDIF



! DEPENDS ON: bdy_expl2
      CALL BDY_EXPL2 (                                                  &

! IN MPP variables
     & halo_i, halo_j, off_x, off_y, row_length, rows, n_rows,          &
     & n_proc,                                                          &

! IN values defining vertical grid of model atmosphere :
     & BL_LEVELS,                                                       &
     & r_rho_levels,r_theta_levels,eta_theta_levels,                    &
     & P,P_theta_levels,                                                &

! IN U, V and W momentum fields.
     & U, V, W, ETADOT,                                                 &
! variables for subgrid turbulence scheme
     & visc_bl_m,delta_lambda, delta_phi,FM_3D, FH_3D,L_subfilter_vert, &
     & L_subfilter_horiz, L_subfilter_blend,max_diff,turb_startlev_vert,&
     & turb_endlev_vert,BL_COEF_KM, BL_COEF_KH,                         &

! IN from other part of explicit boundary layer code
     & RHO_UV,RHO_TQ,RHO_DRY_TQ,DZL_charney,RDZ,                        &
     & Z_FULL,Z_UV,Z_TQ,P_HALF,DELTAP,RHOSTAR,                          &
     & U_P,V_P,BT,BQ,BT_CLD,BQ_CLD,BT_GB,BQ_GB,A_QS,A_DQSDT,DQSDT,      &
     & RECIP_L_MO_SEA,FLANDG,                                           &
     & FORMDRAG, FD_stab_dep, OROG_DRAG_PARAM, RIB_GB, SIL_OROG_LAND,   &

! IN cloud data :
     & CF,Q,QCF,QCL,CCA,CCB,CCT,T,                                      &

! IN everything not covered so far :
     & RAD_HR,MICRO_TENDS,FB_SURF,U_S,T1_SD,Q1_SD,H_BLEND_OROG,         &
     & TIMESTEP, lq_mix_bl, ZH_PREV, L_LAMBDAM2, L_FULL_LAMBDAS,        &
     & BL_OPTIONS,L_SBLeq,L_SBLco,Muw_SBL,Mwt_SBL,                      &

!
! SCM Diagnostics (dummy values in full UM)
     & nSCMDpkgs,L_SCMDiags,                                            &
!
! INOUT variables
     & RDZ_CHARNEY_GRID,RDZ_U,RDZ_V,QW,TL,ZH,                           &

! INOUT data for STPH_RP
     & G0_RP,par_mezcla,                                                &
! OUT Diagnostic not requiring STASH flags :
     & FQW,FTL,RHOKH,                                                   &
     & RHOKM,RHOKM_u,RHOKM_v,TAUX,TAUY,                                 &
     & zht,                                                             &
     & SHALLOWC,CU_OVER_OROG,                                           &
     & BL_TYPE_1,BL_TYPE_2,BL_TYPE_3,BL_TYPE_4,BL_TYPE_5,BL_TYPE_6,     &
     & BL_TYPE_7,                                                       &
     & Z0M_EFF_GB,                                                      &

! OUT data required for tracer mixing :
     & NTML,                                                            &
     & KENT, WE_LIM, T_FRAC, ZRZI,                                      &
     & KENT_DSC, WE_LIM_DSC, T_FRAC_DSC, ZRZI_DSC,                      &
     & ZHSC, Z_HALF,                                                    &

! OUT data required for 4D-VAR :
     & RHO_KM,                                                          &

! OUT Explicit orographic stresses
     & TAUX_FD_U,TAUY_FD_V,                                             &

! OUT data required elsewhere in UM system :
     & NTPAR,NLCL,ZHPAR,Z_LCL,L_SHALLOW,                                &
     & NTDSC,NBDSC,CUMULUS,WSTAR,WTHVS,DELTHVU,                         &
     & UW0,VW0,U_0_P,V_0_P,LAND_PTS,                                    &
     & LAND_INDEX,LAND_MASK,HO2R2_OROG,                                 &
     & LTIMER,BL_diag                                                   &

     &  )


! Add additional orographic stress to surface stress over land

      IF(FORMDRAG  ==  Explicit_stress) THEN

       DO J=1,ROWS
        DO I=1,ROW_LENGTH
          IF(FLANDG_U(I,J) >  0.0)THEN
            TAUX_LAND(I,J) = TAUX_LAND(I,J) + TAUX_FD_U(I,J,1)
          ENDIF
          IF(FLANDG_U(I,J) <  1.0)THEN
            TAUX_SSI(I,J) = TAUX_SSI(I,J) + TAUX_FD_U(I,J,1)
          ENDIF
        END DO
       END DO

       DO J=1,N_ROWS
        DO I=1,ROW_LENGTH

          IF(FLANDG_V(I,J) >  0.0)THEN
            TAUY_LAND(I,J) = TAUY_LAND(I,J) + TAUY_FD_V(I,J,1)
          ENDIF

          IF(FLANDG_V(I,J) <  1.0)THEN
            TAUY_SSI(I,J) = TAUY_SSI(I,J) + TAUY_FD_V(I,J,1)
          ENDIF

        END DO
       END DO

      ENDIF


      DO J=1,ROWS
        DO I=1,ROW_LENGTH
          RIB_LAND(I,J)=0.0
          RIB_SSI(I,J)=0.0
        ENDDO
      ENDDO

      DO N=1,NTILES
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          RIB_LAND(I,J)=RIB_LAND(I,J) +                                 &
     &      RIB_TILE(L,N)*TILE_FRAC(L,N)
        ENDDO
      ENDDO

      DO J=1,ROWS
        DO I=1,ROW_LENGTH
          IF(FLANDG(I,J) <  1.0)                                        &
     &      RIB_SSI(I,J)=(RIB_GB(I,J)-RIB_LAND(I,J)*FLANDG(I,J))        &
     &        /(1.0-FLANDG(I,J))
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
! Mineral dust production
!-----------------------------------------------------------------------
!
       IF (L_DUST) THEN
!
!initialisation
!
         DO IDIV=1,NDIV
           DO J=1,ROWS
             DO I=1,ROW_LENGTH
               DUST_FLUX(I,J,IDIV)=0.0
               U_S_T(I,J,IDIV)=0.0
               U_S_T_DRY(I,J,IDIV)=0.0
             ENDDO !ROW_LENGTH
           ENDDO !ROWS
           DO K = 1,NTILES
             DO L = 1,LAND_PTS
               U_S_T_TILE(L,K,IDIV) = 0.
               U_S_T_DRY_TILE(L,K,IDIV) = 0.
               DUST_FLUX_TILE(L,K,IDIV) = 0.
             ENDDO !LANDPTS
           ENDDO !NTILES
         ENDDO !NDIV

!
!put fields into land arrays
!
         DO L = 1,LAND_PTS
           J = (LAND_INDEX(L)-1)/ROW_LENGTH + 1
           I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
           PSTAR_LAND(L) = PSTAR(I,J)
           RHOSTAR_LAND(L) = RHOSTAR(I,J)
           SAND_LAND(L) = SOIL_SAND(I,J)
           CLAY_LAND(L) = SOIL_CLAY(I,J)
        ENDDO !LAND_PTS

        DO L = 1,LAND_PTS
          J = (LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          MREL_LAND(L,1) = DUST_MREL1(I,J)
          MREL_LAND(L,2) = DUST_MREL2(I,J)
          MREL_LAND(L,3) = DUST_MREL3(I,J)
          MREL_LAND(L,4) = DUST_MREL4(I,J)
          MREL_LAND(L,5) = DUST_MREL5(I,J)
          MREL_LAND(L,6) = DUST_MREL6(I,J)
        ENDDO !LAND_PTS
!
!calculate dust source flux
!
        IF (L_CAM_DUST) THEN
!
! Use NWP CAM version of the dust uplift scheme  
!
! DEPENDS ON: dust_srce_cam
           CALL DUST_SRCE_CAM(                                             &
     & LAND_PTS,NTILES,SM_LEVELS,TILE_PTS,TILE_INDEX,TILE_FRAC,FLAND,   &
     & PSTAR_LAND,TSTAR_TILE,RHOSTAR_LAND,SOIL_LAYER_MOISTURE,SNOW_TILE,&
     & U_S_STD_TILE,MREL_LAND,CLAY_LAND,                                &
     & HO2R2_OROG,                                                      &
     & DUST_FLUX_TILE,U_S_T_TILE,U_S_T_DRY_TILE,                        &
     & SAND_LAND                                                        &
     &   )
        ELSE
!
! Use HadGEM version of the dust uplift scheme  
!
! DEPENDS ON: dust_srce
           CALL DUST_SRCE(                                                 &
     & LAND_PTS,NTILES,SM_LEVELS,TILE_PTS,TILE_INDEX,TILE_FRAC,FLAND,   &
     & PSTAR_LAND,TSTAR_TILE,RHOSTAR_LAND,SOIL_LAYER_MOISTURE,SNOW_TILE,&
     & U_S_STD_TILE,MREL_LAND,CLAY_LAND,                                &
     & HO2R2_OROG,                                                      &
     & DUST_FLUX_TILE,U_S_T_TILE,U_S_T_DRY_TILE,                        &
     & SAND_LAND                                                        &
     &   )
        ENDIF !L_CAM_DUST
!
!calculate gridbox mean values
!
        DO IDIV = 1,NDIV
          DO M = 1,NTILES
            DO N = 1,TILE_PTS(M)
              L = TILE_INDEX(N,M)
              J = (LAND_INDEX(L)-1)/ROW_LENGTH + 1
              I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
              DUST_FLUX(I,J,IDIV) = DUST_FLUX(I,J,IDIV) +               &
     &        DUST_FLUX_TILE(L,M,IDIV)*TILE_FRAC(L,M)
               U_S_T(I,J,IDIV) = U_S_T(I,J,IDIV) +                      &
     &        U_S_T_TILE(L,M,IDIV)*TILE_FRAC(L,M)
               U_S_T_DRY(I,J,IDIV) = U_S_T_DRY(I,J,IDIV) +              &
     &        U_S_T_DRY_TILE(L,M,IDIV)*TILE_FRAC(L,M)
            ENDDO !TILE_PTS
          ENDDO !NTILES
        ENDDO !NDIV
!
!
      ENDIF !L_DUST
!

  999  CONTINUE  ! Branch for error exit.

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('BDYLAYR ',4)
      ENDIF

      RETURN
      END SUBROUTINE BDY_LAYR
#endif
