#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!! Subroutine BL_INTCT -------------------------------------------
!!!
!!! Purpose : Intermediate control level to call requested version of
!!!           BDY_LAYR with the appropriate arguments.
!!!
!!! Level 3 control routine
!!! version for CRAY YMP
!!!
!!! Programming standard : unified model documentation paper No 3
!!!
!!! System components covered : P24
!!!
!!! System task : P0
!!!
!!!END -----------------------------------------------------------------

!    Arguments :-
      SUBROUTINE BL_INTCT (                                             &

! IN MPP variables
     & halo_i, halo_j, off_x, off_y, row_length, rows, n_rows,          &
     & global_row_length, proc_row_group, at_extremity,                 &
     & n_proc, n_procx, n_procy, neighbour,                             &

! IN Substep number for ATMPHYS2
     & Substep_Number, Num_Substeps, L_phys2_substep,NumCycles,CycleNo, &


! IN values defining field dimensions and subset to be processed :
     & NTILES,land_pts,                                                 &

! IN values defining vertical grid of model atmosphere :
     & MODEL_DOMAIN,                                                    &
     & BL_LEVELS,                                                       &
     & r_rho_levels,r_theta_levels,eta_theta_levels,                    &
     & P,P_theta_levels, rho_rsq, rho_wet, rho_dry,                     &
     & GAMMA,SIN_THETA_LONGITUDE,COS_THETA_LONGITUDE,                   &
     & sin_theta_latitude,                                              &

! IN U and V momentum fields.
     & U, V, W, ETADOT, U_conv, V_conv,                                 &
! IN Non turbulent increments to momentum
!  (New dynamics only).
     & DU_NT,DV_NT,                                                     &

! variables for subgrid turbulence scheme
     & visc_BL_m,delta_lambda,delta_phi, FM_3D, FH_3D, L_subfilter_vert,&
     & L_subfilter_horiz,L_subfilter_blend, max_diff,turb_startlev_vert,&
     & turb_endlev_vert,BL_COEF_KM, BL_COEF_KH,                         &
! IN soil/vegetation/land surface data :
     & LAND_MASK,LAND_INDEX,ST_LEVELS,SM_LEVELS,                        &
     & DIM_CS1, DIM_CS2,                                                &
     & HCON,SMVCCL,SMVCST,SMVCWT,STHF,STHU,                             &
     & SIL_OROG_LAND,FORMDRAG,OROG_DRAG_PARAM,HO2R2_OROG,               &
     & SOIL_LAYER_MOISTURE,                                             &

! IN sea/sea-ice data :
     & ICE_FRACT,U_0,V_0,U_0_P,V_0_P,Charnock,SeaSalinityFactor,        &

! IN cloud data :
     & CF,Q,QCF,QCL,QCF_latest,QCL_latest,CCA,CCB,CCT, T,               &

! IN everything not covered so far :
     & L_LAMBDAM2,L_FULL_LAMBDAS,                                       &
     & CO2_MMR,PHOTOSYNTH_ACT_RAD,PSTAR,SURF_RADFLUX,                   &
     & RAD_HR,MICRO_TENDS,TIMESTEP,lq_mix_bl,ZH_PREV,                   &
     & L_ctile,BL_OPTIONS,L_SBLeq,L_SBLco,Muw_SBL,Mwt_SBL,              &


! IN Variables for: prescribed surface flux forcing
!                   prescribed sea-surface roughness lengths
     & flux_e, flux_h, L_flux_bc, L_spec_z0, Z0M_SCM, Z0H_SCM,          &

! IN variables required for mineral dust scheme
     & L_DUST, L_CAM_DUST,                                              &
     & SOIL_CLAY,SOIL_SILT,SOIL_SAND,                                   &
     & DUST_MREL1,DUST_MREL2,DUST_MREL3,                                &
     & DUST_MREL4,DUST_MREL5,DUST_MREL6,                                &

! IN additional variables for MOSES II
     & LAND_PTS_TRIF,NPFT_TRIF,                                         &
     & CANOPY,CATCH,CATCH_SNOW,                                         &
     & SNOW_TILE,Z0_TILE,LW_DOWN,                                       &
     & SW_TILE,TSTAR_TILE,                                              &
     & CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,L_CO2_INTERACTIVE,                &
     & L_PHENOL,L_TRIFFID,L_Q10,ASTEPS_SINCE_TRIFFID,CAN_MODEL,         &
     & CS,FRAC,CANHT_FT,LAI_FT,                                         &
     & FLAND,FLANDG,TSTAR_SEA,Z_LAND,ANTHROP_HEAT,                      &
     & ALBSOIL,COS_ZENITH_ANGLE,CAN_RAD_MOD,ILAYERS,                    &

    
!    EAK      
!    IN       
     &        surf_down_sw,alb_tile,l_tile_pts,                         &
!     &        surf_down_sw,alb_tile,cos_zenith_angle,        &
     &        ls_rain,ls_snow,SW_DOWN,                                  &
     &        lat,long,day,time_sec,                                    &
     &  SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L,                  &
     &  SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE,            &
     &  TSOIL_TILE,T_SURF_TILE,HCONS,BEXP,                              &
     &  SATHH,SATCON,HCAP,                                              &
     &  SOIL_TYPE,VEG_TYPE,                                             &
     &  ISNOW_FLG3L,total_nsteps,                                       &
     &            FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB,              &
     &            TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB,          &
     &            USTAR_CAB,SURF_HTF_CAB,                               &
     &            U_S_CAB,CH_CAB,CD_CAB,                                &
     &            l_cable,                                              &
!sxy
     &            SNAGE_TILE,RTSOIL_TILE,                               &
     &            GFLUX_TILE,SGFLUX_TILE,                               &
! Lestevens Sept2012: CasaCNP variables 
     &            CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,                     &
     &            SOIL_ORDER,NIDEP,NIFIX,PWEA,PDUST,                    &
     &            GLAI,PHENPHASE,                                       &

! INOUT data :
! Variables for STPH_RP (G0_RP,par_mezcla)
     & T_SOIL,TI,TSTAR,Z0MSEA, T_latest,Q_latest,G0_RP,par_mezcla,      &

! INOUT additional variables for MOSES II
     & GS,G_LEAF_ACC,NPP_FT_ACC,RESP_W_FT_ACC,RESP_S_ACC,               &
     & TSTAR_SICE,TSTAR_SSI,                                            &

! OUT Diagnostic not requiring STASH flags :
     & CD,CH,E_SEA,FQW,FTL,H_SEA,RHOKH,RHOKM,RHOKM_u,                   &
     & RHOKM_v, RIB_GB,TAUX,TAUY,VSHR,zht,shallowc,cu_over_orog,        &
     & BL_TYPE_1,BL_TYPE_2,BL_TYPE_3,BL_TYPE_4,BL_TYPE_5,BL_TYPE_6,     &
     & bl_type_7,                                                       &
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

!OUT variables required for mineral dust scheme
     & R_B_DUST,DUST_FLUX,                                              &
     & U_S_T_TILE,U_S_T_DRY_TILE,U_S_STD_TILE,                          &
     & NTML,                                                            &
     & KENT, WE_LIM, T_FRAC, ZRZI,                                      &
     & KENT_DSC, WE_LIM_DSC, T_FRAC_DSC, ZRZI_DSC,                      &
     & ZHSC, Z_HALF,                                                    &

! OUT data required for 4D-VAR :
     & RHO_CD_MODV1,RHO_KM,                                             &

! OUT variables required in IMP_SOLVER
     & ALPHA1_SICE,ASHTF,                                               &
     & DTRDZ_CHARNEY_GRID,RDZ_CHARNEY_GRID,                             &
     & DTRDZ_U,DTRDZ_V,RDZ_U,RDZ_V,                                     &
     & CDR10M_U,CDR10M_V,                                               &
     & Z1_TQ,                                                           &
! OUT variables which need to be maintained between substeps
     & RHO_UV,RHO_TQ,DZL_CHARNEY,RDZ,                                   &
     & Z1_UV,Z_FULL,Z_UV,Z_TQ,P_HALF,DELTAP,                            &

! OUT additional variables for MOSES II
     & FTL_TILE,LE_TILE,RADNET_SICE,                                    &
     & RADNET_TILE,RIB_TILE,RHO_ARESIST_TILE,ARESIST_TILE,              &
     & RESIST_B_TILE,ALPHA1,ASHTF_TILE,FQW_TILE,EPOT_TILE,FQW_ICE,      &
     & FTL_ICE,FRACA,RESFS,RESFT,RHOKH_TILE,RHOKH_SICE,                 &
     & RHOKPM,RHOKPM_POT,RHOKPM_SICE,                                   &
     & Z0HSSI,Z0H_TILE,Z0M_GB,Z0MSSI,Z0M_TILE,CHR1P5M,CHR1P5M_SICE,     &
     & SMC,GPP,NPP,RESP_P,G_LEAF,GPP_FT,NPP_FT,                         &
     & RESP_P_FT,RESP_S,RESP_S_TOT,RESP_W_FT,                           &
     & GC,CANHC_TILE,WT_EXT_TILE,FLAKE,                                 &
     & TILE_INDEX,TILE_PTS,TILE_FRAC,FSMC,                              &
     & TSTAR_LAND,RIB_SSI,TAUX_LAND,TAUX_SSI,TAUY_LAND,TAUY_SSI,        &
     & VSHR_LAND,VSHR_SSI,                                              &
     & FLANDG_U,FLANDG_V,                                               &

! OUT data required elsewhere in UM system :
     & ZH,T1_SD,Q1_SD,                                                  &
     & NTPAR,NLCL,ZHPAR,Z_LCL,L_SHALLOW,                                &
     & NTDSC,NBDSC,CUMULUS,WSTAR,WTHVS,DELTHVU,                         &
     & uw0,vw0,                                                         &
     & L_ukca,                                                          &
     & ERROR,LTIMER,BL_diag,                                            &
      ! end step of experiment, this step, step width, processor num
     & endstep, timestep_number, mype )

      Use bl_diags_mod, Only :                                          &
          strnewbldiag

      IMPLICIT NONE

#include "c_dust_ndiv.h"


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

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         !  south,east or west of the processor grid

! Substep nummber for ATMPHYS2
      Integer                                                           &
     & Substep_Number                                                   &
     &,Num_Substeps

! Switch for calculating exchange coeffs from latest values.
      Logical                                                           &
     & L_phys2_substep
      LOGICAL                                                           &
     & lq_mix_bl

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
      INTEGER                                                           &
     & FORMDRAG     ! IN switch for orographic form drag

! (b) Defining vertical grid of model atmosphere.

      INTEGER                                                           &
     & BL_LEVELS                                                        &
                                   ! IN Max. no. of "boundary" levels
!                                     allowed.Assumed <= 30 for dim-
!                                     sioning of GAMMA in common deck
!                                     C_GAMMA used in SF_EXCH and KMKH
     &,MODEL_DOMAIN

!
! Boundary Layer
      REAL, Intent(IN) :: SeaSalinityFactor
!                                  ! Factor allowing for the effect
!                                  ! of the salinity of sea water
!                                  ! on the evaporative flux.
!
      REAL, Intent(IN) :: orog_drag_param
!                                  ! Drag coefficient for orographic
!                                  ! form drag scheme
!
      INTEGER, DIMENSION(20) :: BL_OPTIONS   ! IN BL switches

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
     &, W(row_length,rows,0:BL_LEVELS)                                  &
     &, ETADOT(row_length,rows,0:BL_LEVELS)                             &
     &,  U_conv(1-off_x:row_length+off_x,1-off_y:rows+off_y,            &
     &              bl_levels)                                          &
     &,  V_conv(1-off_x:row_length+off_x,1-off_y:n_rows+off_y,          &
     &              bl_levels)                                          &
     &, visc_BL_m(1:row_length, 1:rows, bl_levels)                      &
!            ! visc_m only on BL levels
     &, delta_lambda                                                    &
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
!            ! RHOKM from BL scheme
     &, BL_COEF_KH(1:row_length, 1:rows, bl_levels-1)
!            ! RHOKH from BL scheme

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
                            ! 3D turbulence scheme.

! (c) Soil/vegetation/land surface parameters (mostly constant).

      LOGICAL                                                           &
     & LAND_MASK(row_length,rows)                                       &
                                   ! IN T if land, F elsewhere.
     &,L_CTILE                                                          &
                                   ! IN Switch for coastal tiling
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
     & HCON(land_pts)                                                   &
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
     &, SOIL_LAYER_MOISTURE(LAND_PTS,SM_LEVELS)!IN soil moisture
!                 ! per layer (kg m-2)

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
     &,QCF(1-halo_i:row_length+halo_i,                                  &
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
                                    ! IN Cloud ice (kg kg^-1 air)
     &,QCL_latest(row_length,rows,BL_LEVELS)                            &
                                    ! IN Cloud liquid water
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
     &,SURF_RADFLUX(row_length,rows)                                    &
                                    ! IN Surface net radiation (W/sq m,
!                                     positive downwards).
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
     &, Charnock                                                        &
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
     &,ZH_PREV(row_length, rows)                                        &
                                    ! IN boundary layer height from
!                                   !    previous timestep
     &,flux_e(row_length,rows)                                          &
     &,flux_h(row_length,rows)
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
     &,L_CAM_DUST !Old version of dust_uplift scheme used in CAM NWP models


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
     & CANOPY(LAND_PTS,NTILES)                                          &
                                   ! IN Surface/canopy water for
!                                  !    snow-free land tiles (kg/m2)
     &,CATCH(LAND_PTS,NTILES)                                           &
                                   ! IN Surface/canopy water capacity
!                                  !    of snow-free land tiles (kg/m2).
     &,CATCH_SNOW(LAND_PTS)                                             &
                                   ! IN Snow interception capacity of
!                                  !    NLT tile (kg/m2).
     &,SNOW_TILE(LAND_PTS,NTILES)                                       &
                                   ! IN Lying snow on tiles (kg/m2)
     &,Z0_TILE(LAND_PTS,NTILES)                                         &
                                   ! IN Tile roughness lengths (m).
     &,LW_DOWN(ROW_LENGTH,ROWS)                                         &
                                   ! IN Surface downward LW radiation
!                                  !    (W/m2).
     &,SW_DOWN(ROW_LENGTH,ROWS)                                         &
                                   ! Surface downward SW radiation (W/m2).
     &,SW_TILE(LAND_PTS,NTILES)
                                   ! IN Surface net SW radiation on
!                                  !    land tiles (W/m2).
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
! Lestevens Sept2012: CasaCNP variables 
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
     &, lat(row_length, rows)                                           &
     &, long(row_length, rows)                                          &
     &, time_sec                                                        

      Real                                                                  &
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
      Logical                                                           &
     &  l_tile_pts(land_pts,ntiles)
!------------------------------------------------------

      REAL                                                              &
     & TSTAR_TILE(LAND_PTS,NTILES)                                      &
                                   ! IN Surface tile temperatures
     &,CO2_3D(CO2_DIM_LEN,CO2_DIM_ROW)                                  &
!                                  ! IN 3D CO2 field if required.
     &,CS(LAND_PTS,DIM_CS1)                                             &
                              ! IN Soil carbon (kg C/m2).
     &,FRAC(LAND_PTS,NTYPE)                                             &
                                   ! IN Fractions of surface types.
     &,CANHT_FT(LAND_PTS,NPFT)                                          &
                                   ! IN Canopy height (m)
     &,LAI_FT(LAND_PTS,NPFT)                                            &
                                   ! IN Leaf area index
     &,FLAND(LAND_PTS)                                                  &
                                   ! IN Land fraction on land tiles.
     &,FLANDG(ROW_LENGTH,ROWS)                                          &
     &,TSTAR_SEA(ROW_LENGTH,ROWS)                                       &
                                   ! IN Open sea sfc temperature (K).
!                                  ! IN Land fraction on all tiles.
     &,Z_LAND(ROW_LENGTH,ROWS)                                          &
                                   ! IN Land height (m).
     &,ALBSOIL(LAND_PTS)                                                &
!                                  ! Soil albedo.
     &, COS_ZENITH_ANGLE(row_length, rows)
!                                  ! Cosine of the zenith angle
      INTEGER                                                           &
     & CAN_RAD_MOD                                                      &
                                   !Switch for canopy radiation model
     & ,ILAYERS                    !No of layers in canopy radiation model



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
!     Declaration of new BL diagnostics.
      Type (Strnewbldiag) :: BL_diag


! Additional variables for SCM diagnostics which are dummy in full UM
      INTEGER                                                           &
     & nSCMDpkgs             ! No of SCM diagnostics packages

      LOGICAL                                                           &
     & L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages
!
!  In/outs :-

! Definition of variables for STPH_RP
      REAL,INTENT(InOut) :: par_mezcla ! Used to modify LAMBDAH,LAMBDAM
                                       !(neutral mixing length) EXCOEF
      REAL,INTENT(InOut) :: G0_RP!Stability function parameter in EXCOEF
      REAL                                                              &
     & Q_latest(row_length,rows,BL_LEVELS)                              &
                                            ! IN specific humidity
     &,T_latest(row_length,rows,BL_LEVELS)                              &
                                            ! IN temperature
     &,T_SOIL(land_pts,SM_LEVELS)                                       &
                                   ! INOUT Soil temperatures (K).
     &,TI(row_length,rows)                                              &
                                   ! INOUT Sea-ice surface layer
!                                      temperature (K).
     &,TSTAR(row_length,rows)                                           &
                                   ! INOUT Surface temperature (K).
     &,Z0MSEA(row_length,rows)                                          &
                                   ! INOUT Sea-surface roughness
!                                      length for momentum (m).
!                                      NB: same storage is used
!                                      for Z0V, so the intent is
!                                      IN for land points.
     &,Z0M_SCM(row_length,rows)                                         &
                                   ! IN Fixed Sea-surface roughness
                                   ! length(m) for momentum (SCM)
     &,Z0H_SCM(row_length,rows)    ! IN Fixed Sea-surface roughness
                                   ! length(m) for heat (SCM)

! Additional MOSES II variables
      REAL                                                              &
     & GS(LAND_PTS)                                                     &
                                   ! INOUT "Stomatal" conductance to
!                                  !        evaporation (m/s).
     &,G_LEAF_ACC(LAND_PTS,NPFT)                                        &
                                   ! INOUT Accumulated G_LEAF
! Lestevens 17may13 - changed npft to ntiles
     &,NPP_FT_ACC(LAND_PTS_TRIF,NTILES)                                 &
!    &,NPP_FT_ACC(LAND_PTS_TRIF,NPFT_TRIF)                              &
!                                  ! INOUT Accumulated NPP_FT
     &,RESP_W_FT_ACC(LAND_PTS_TRIF,NTILES)                              &
!    &,RESP_W_FT_ACC(LAND_PTS_TRIF,NPFT_TRIF)                           &
!                                  ! INOUT Accum RESP_W_FT
     &,RESP_S_ACC(LAND_PTS_TRIF,DIM_CS1)                                &
                                         ! INOUT Accumulated RESP_S
     &,TSTAR_SICE(ROW_LENGTH,ROWS)                                      &
                                   ! INOUT Sea-ice sfc temperature (K).
     &,TSTAR_SSI(ROW_LENGTH,ROWS)  ! INOUT Sea mean sfc temperature (K).



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
     &,TAUX(row_length,rows,BL_LEVELS)                                  &
                                   ! OUT W'ly component of surface wind
!                                     stress (N/sq m).(On UV-grid with
!                                     first and last rows undefined or
!                                     at present, set to missing data
     &,TAUY(row_length,n_rows,BL_LEVELS)                                &
                                   ! OUT S'ly component of surface wind
!                                     stress (N/sq m).  On UV-grid;
!                                     comments as per TAUX.
     &,VSHR(row_length,rows)                                            &
                                   ! OUT Magnitude of surface-to-lowest
!                                     atm level wind shear (m per s).
     &,RHO_CD_MODV1(row_length,rows)                                    &
                                   ! OUT Surface air density * drag coef
!                                  !  *mod(v1 - v0) before interpolation
     &, shallowc(row_length,rows)                                       &
                                   !OUT Indicator set to 1.0 if shallow,
!                                  !   0.0 if not shallow or not cumulus
     &, cu_over_orog(row_length,rows)                                   &
!                                  ! OUT Indicator for cumulus
!                                  !     over steep orography
!                                  !   Indicator set to 1.0 if true,
!                                  !   0.0 if false. Exclusive.
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
     &, bl_type_7(row_length,rows)                                      &
                                   !  Indicator set to 1.0 if a
!                                  !     shear-dominated b.l.
!                                  !      diagnosed, 0.0 otherwise.
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
     &, R_B_DUST(ROW_LENGTH,ROWS,NDIV)                                  &
                                      !OUT surface layer resist for dust
     &, DUST_FLUX(ROW_LENGTH,ROWS,NDIV)                                 &
                                       !OUT dust emissions (kg m-2 s-1)
     &, U_S_T_TILE(LAND_PTS,NTILES,NDIV)                                &
                                        !OUT threshold frict. vel
     &, U_S_T_DRY_TILE(LAND_PTS,NTILES,NDIV)                            &
                                            !OUT dry soil value
     &, U_S_STD_TILE(LAND_PTS,NTILES)                                   &
                                     !OUT friction velocity

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
     &, WT_EXT(LAND_PTS,sm_levels)                                      &
                                  !OUT cumulative fraction of transp'n
     &, RA(LAND_PTS)               !OUT Aerodynamic resistance (s/m)

!*APL*DIAGS
      REAL                                                              &
     &  ZHT(row_length, rows)      ! OUT Max height of turb mixing

      INTEGER                                                           &
     & NTML(row_length,rows)                                            &
                                   ! OUT Number of model levels in the
!                                           turbulently mixed layer.
     &,NTPAR(row_length,rows)                                           &
                                 ! INOUT Top level of initial parcel
!                                !  ascent. Used in convection scheme.
     &,NLCL(row_length,rows)                                            &
                                 ! INOUT No of levels to LCL
!
     &,NTDSC(row_length,rows)                                           &
                                   ! OUT Top level for turb mixing in
!                                           any decoupled Sc layer
     &,NBDSC(row_length,rows)                                           &
                                   ! OUT Bottom level of any decoupled
!                                  !     turbulently-mixed Sc layer.
     &, KENT(row_length,rows)                                           &
                                    ! OUT grid-level of SML inversion
     &, KENT_DSC(row_length,rows)   ! OUT grid-level of DSC inversion

      LOGICAL CUMULUS(row_length,rows)                                  &
                                   ! OUT Logical switch for trade Cu
     &, l_shallow(row_length, rows)      ! Logical indicator of
!                                        ! shallow convection
      Real                                                              &
     &  delthvu(row_length, rows)                                       &
                                     ! buoyancy integral
     &, zhpar(row_length, rows)                                         &
                                     ! height of ntpar
     &, z_lcl(row_length, rows)                                         &
                                      ! height of nlcl
     &, wstar(row_length, rows)                                         &
                                     ! surface-based mixed layer
!                                    ! velocity scale
     &, wthvs(row_length, rows)      ! surface buoyancy flux

      Real                                                              &
     &  uw0(row_length, rows)                                           &
!                       ! U-component of surface wind stress (P-grid)
     &, vw0(row_length, rows)
!                       ! V-component of surface wind stress (P-grid)

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

! OUT variables for IMP_SOLVER (that used to be local arrays)
      REAL                                                              &
     & ALPHA1_SICE(ROW_LENGTH,ROWS)                                     &
                                   ! OUT ALPHA1 for sea-ice.
!                                  !     Grid-box average for MOSES I
     &,ASHTF(row_length,rows)                                           &
                                  ! OUT Coefficient to calculate surface
!                                 heat flux into soil or sea-ice.
     &,DTRDZ_CHARNEY_GRID(row_length,rows,BL_LEVELS)                    &
!                                       ! OUT -g.dt/dp for model layers.
     &,RDZ_CHARNEY_GRID(row_length,rows,BL_LEVELS)                      &
!                               ! OUT RDZ(,1) is the reciprocal of the
!                                 height of level 1, i.e. of the
!                                 middle of layer 1.  For K > 1,
!                                 RDZ(,K) is the reciprocal
!                                 of the vertical distance
!                                 from level K-1 to level K.
     &,DTRDZ_U(row_length,rows,BL_LEVELS)                               &
                                          ! OUT
     &,DTRDZ_V(row_length,n_rows,BL_LEVELS)                             &
                                           ! OUT
!                                 -g.dt/dp for model wind layers.
     &,RDZ_U(row_length,rows,2:BL_LEVELS)                               &
                                ! OUT RDZ (K > 1) on UV-grid.
     &,RDZ_V(row_length,n_rows,2:BL_LEVELS)                             &
                                ! OUT RDZ (K > 1) on UV-grid.
     &,CDR10M_U(row_length,rows)                                        &
                                ! OUT Ratio of CD's reqd for calculation
     &,CDR10M_V(row_length,n_rows)                                      &
                                ! OUT Ratio of CD's reqd for calculation
     &,Z1_TQ(row_length,rows)   ! OUT Height of lowest theta level.


! OUT variables which need to be maintained between substeps. They will
!     be used by BDY_LAYR()
      REAL                                                              &
     & RHO_UV(row_length,rows,BL_LEVELS+1)                              &
!                               ! RHO_UV(*,K) is the density at half
!                               ! level k-1/2.
     &,RHO_TQ(row_length,rows,BL_LEVELS)                                &
!                               ! RHO_TQ(*,K) is the density at half
!                               ! level k+1/2.
     &,DZL_charney(row_length,rows,BL_LEVELS)                           &
                                ! DZL(,K) is depth in m of layer
!                                 K, i.e. distance from boundary
!                                 K-1/2 to boundary K+1/2.
     &,RDZ(row_length,rows,BL_LEVELS)                                   &
                                ! RDZ(,1) is the reciprocal of the
!                                 height of level 1, i.e. of the
!                                 middle of layer 1.  For K > 1,
!                                 RDZ(,K) is the reciprocal
!                                 of the vertical distance
!                                 from level K-1 to level K.
     &,Z1_UV(row_length,rows)                                           &
                                ! Height of lowest u,v level.
     &,Z_FULL(row_length,rows,BL_LEVELS)                                &
                                ! Z_FULL(*,K) is height of full level k.
     &,Z_UV(row_length,rows,BL_LEVELS)                                  &
                                ! Z_UV(*,K) is height of half level
!                               ! k-1/2.
     &,Z_TQ(row_length,rows,BL_LEVELS)                                  &
                                ! Z_TQ(*,K) is height of half level
!                               ! k+1/2.
     &,P_HALF(row_length,rows,BL_LEVELS)                                &
                                ! P_HALF(*,K) is pressure at half
!                               ! level k-1/2.
     &,DELTAP(row_length,rows,BL_LEVELS)
                                ! Difference in pressure between levels

! Additional MOSES II variables
      INTEGER                                                           &
     & TILE_INDEX(LAND_PTS,NTYPE)                                       &
                                   ! OUT Index of tile points
     &,TILE_PTS(NTYPE)             ! OUT Number of tile points


      REAL                                                              &
     & FTL_TILE(LAND_PTS,NTILES)                                        &
                                   ! OUT Surface FTL for land tiles
     &,LE_TILE(LAND_PTS,NTILES)                                         &
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
     &,ALPHA1(LAND_PTS,NTILES)                                          &
                                   ! OUT Mean gradient of saturated
!                                  !     specific humidity with respect
!                                  !     to temperature between the
!                                  !     bottom model layer and tile
!                                  !     surfaces
     &,ASHTF_TILE(LAND_PTS,NTILES)                                      &
                                   !OUT Coefficient to calculate
!                                  !     surface heat flux into land
!                                  !     tiles.
     &,FQW_TILE(LAND_PTS,NTILES)                                        &
                                   ! OUT Surface FQW for land tiles
     &,EPOT_TILE(LAND_PTS,NTILES)                                       &
                                   ! OUT Local EPOT for land tiles.
     &,FQW_ICE(ROW_LENGTH,ROWS)                                         &
                                   ! OUT Surface FQW for sea-ice
     &,FTL_ICE(ROW_LENGTH,ROWS)                                         &
                                   ! OUT Surface FTL for sea-ice
     &,FRACA(LAND_PTS,NTILES)                                           &
                                   ! OUT Fraction of surface moisture
!                                  !     flux with only aerodynamic
!                                  !     resistance for snow-free land
!                                  !     tiles.
     &,RESFS(LAND_PTS,NTILES)                                           &
                                   ! OUT Combined soil, stomatal
!                                  !     and aerodynamic resistance
!                                  !     factor for fraction (1-FRACA)
!                                  !     of snow-free land tiles.
     &,RESFT(LAND_PTS,NTILES)                                           &
                                   ! OUT Total resistance factor.
!                                  !     FRACA+(1-FRACA)*RESFS for
!                                  !     snow-free land, 1 for snow.
     &,RHOKH_TILE(LAND_PTS,NTILES)                                      &
                                   ! OUT Surface exchange coefficients
!                                  !     for land tiles
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
     &,Z0HSSI(ROW_LENGTH,ROWS)                                          &
     &,Z0MSSI(ROW_LENGTH,ROWS)                                          &
                                   ! OUT Roughness lengths over sea (m).
     &,Z0H_TILE(LAND_PTS,NTILES)                                        &
                                   ! OUT Tile roughness lengths for heat
!                                  !     and moisture (m).
     &,Z0M_GB(ROW_LENGTH,ROWS)                                          &
                                   ! OUT Gridbox mean roughness length 
!                                  !     for momentum (m).
     &,Z0M_TILE(LAND_PTS,NTILES)                                        &
                                   ! OUT Tile roughness lengths for
!                                  !     momentum.
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
     &,GPP(LAND_PTS)                                                    &
                                   ! OUT Gross primary productivity
!                                  !     (kg C/m2/s).
     &,NPP(LAND_PTS)                                                    &
                                   ! OUT Net primary productivity
!                                  !     (kg C/m2/s).
     &,RESP_P(LAND_PTS)                                                 &
                                   ! OUT Plant respiration (kg C/m2/s).
     !kdcorbin, 11/10 - changed from NPFT
     &,G_LEAF(LAND_PTS,NTILES)                                           &
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
     &,FSMC(LAND_PTS,NPFT)                                              &
                                   ! OUT Moisture availability factor.
     &,TSTAR_LAND(ROW_LENGTH,ROWS)                                      &
                                   ! OUT Land mean sfc temperature (K)
     &,RIB_SSI(row_length,rows)                                         &
                                   ! OUT Sea mean bulk Richardson number
!                                        for lowest layer.
     &,TAUX_LAND(ROW_LENGTH,ROWS)                                       &
                                   ! OUT W'ly component of land sfc wind
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
     &,TAUY_LAND(ROW_LENGTH,N_ROWS)                                     &
                                   ! OUT S'ly component of land sfc wind
!                                  !     stress (N/sq m).  On V-grid;
!                                  !     comments as per TAUX.
     &,TAUY_SSI(ROW_LENGTH,N_ROWS)                                      &
                                   ! OUT S'ly compt of mean sea sfc wind
!                                  !     stress (N/sq m).  On V-grid;
!                                  !     comments as per TAUX.
     &,VSHR_LAND(row_length,rows)                                       &
                                   ! OUT Magnitude of land sfc-to-lowest
!                                        atm level wind shear (m per s).
     &,VSHR_SSI(row_length,rows)                                        &
                                   ! OUT Mag. of mean sea sfc-to-lowest
!                                        atm level wind shear (m per s).
     &,FLANDG_U(ROW_LENGTH,ROWS)                                        &
                                   ! OUT Land frac (on U-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
     &,FLANDG_V(ROW_LENGTH,N_ROWS)                                      &
                                   ! OUT Land frac (on V-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
     &,ANTHROP_HEAT(NTILES)        ! IN Additional heat source on urban
                                   !    tiles for anthropogenic heat
                                   !    source (W/m2) 


      INTEGER                                                           &
     & ERROR                    ! OUT 0 - AOK;
!                               !     1 to 7  - bad grid definition
                                !     detected;
      ! end step of experiment, this step, step width, processor num
      integer :: endstep, timestep_number,mype
      !timestep_width = timestep

!---------------------------------------------------------------------
!  External routines called :-

      EXTERNAL BDY_LAYR
      EXTERNAL TIMER

!-----------------------------------------------------------------------

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('BL_INTCT ',3)
      ENDIF

      IF ( .NOT. L_PHYS2_SUBSTEP ) THEN

! DEPENDS ON: bdy_layr
      CALL BDY_LAYR (                                                   &

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
     & CF,Q(1:row_length,1:rows,1:bl_levels),                           &
     & QCF(1:row_length,1:rows,1:bl_levels),                            &
     & QCL(1:row_length,1:rows,1:bl_levels),                            &
     & CCA,CCB,CCT, T,                                                  &

! IN everything not covered so far :
     & L_LAMBDAM2,L_FULL_LAMBDAS,                                       &
     & CO2_MMR,PHOTOSYNTH_ACT_RAD,PSTAR,SURF_RADFLUX,                   &
     & RAD_HR,MICRO_TENDS,TIMESTEP, lq_mix_bl, ZH_PREV,                 &
     & BL_OPTIONS,L_SBLeq,L_SBLco,Muw_SBL,Mwt_SBL,                      &
! IN Variables for: Prescribed surface flux forcing
!                   Prescribed sea-surface roughness lengths
     & flux_e, flux_h, L_flux_bc, L_spec_z0, Z0M_SCM, Z0H_SCM,          &


! IN variables required for mineral dust scheme
     & L_DUST, L_CAM_DUST,                                              &
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
     & ls_rain,ls_snow,SW_DOWN,                                         &
     & lat,long,day,time_sec,                                           &
     & SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L,                   &
     & SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE,             &
     & TSOIL_TILE,T_SURF_TILE,HCONS,BEXP,                               &
     & SATHH,SATCON,HCAP,                                               &
     & SOIL_TYPE,VEG_TYPE,                                              &
     & ISNOW_FLG3L,total_nsteps,                                        &
     &            FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB,              &
     &            TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB,          &
     &            USTAR_CAB,SURF_HTF_CAB,                               &
     &            U_S_CAB,CH_CAB,CD_CAB,                                &
     &            l_cable,                                              &
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
     & bl_type_7,                                                       &
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

! Call explicit BL when substepping is active.
      ELSE

! DEPENDS ON: bdy_layr
      CALL BDY_LAYR (                                                   &

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

! IN U and V momentum fields
     & U_conv, V_conv, W, ETADOT,                                       &
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
     & CF,Q_latest,QCF_latest,QCL_latest,CCA,CCB,CCT, T_latest,         &

! IN everything not covered so far :
     & L_LAMBDAM2,L_FULL_LAMBDAS,                                       &
     & CO2_MMR,PHOTOSYNTH_ACT_RAD,PSTAR,SURF_RADFLUX,                   &
     & RAD_HR,MICRO_TENDS,TIMESTEP,lq_mix_bl,ZH_PREV,                   &
     & BL_OPTIONS,L_SBLeq,L_SBLco,Muw_SBL,Mwt_SBL,                      &

! IN Variables for: Prescribed surface flux forcing
!                   Prescribed sea-surface roughness lengths
     & flux_e, flux_h, L_flux_bc, L_spec_z0, Z0M_SCM, Z0H_SCM,          &
! IN variables required for mineral dust scheme
     & L_DUST, L_CAM_DUST,                                              &
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
     & ls_rain,ls_snow,SW_DOWN,                                         &
     & lat,long,day,time_sec,                                           &
     & SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L,                   &
     & SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE,             &
     & TSOIL_TILE,T_SURF_TILE,HCONS,BEXP,                               &
     & SATHH,SATCON,HCAP,                                               &
     & SOIL_TYPE,VEG_TYPE,                                              &
     & ISNOW_FLG3L,total_nsteps,                                        &
     &            FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB,              &
     &            TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB,          &
     &            USTAR_CAB,SURF_HTF_CAB,                               &
     &            U_S_CAB,CH_CAB,CD_CAB,                                &
     &            l_cable,                                              &
     &            SNAGE_TILE,RTSOIL_TILE,                               &
     &            GFLUX_TILE,SGFLUX_TILE,                               &
! Lestevens Sept2012: CasaCNP variables 
     &            CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,                     &
     &            SOIL_ORDER,NIDEP,NIFIX,PWEA,PDUST,                    &
     &            GLAI,PHENPHASE,                                       &

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
     & bl_type_7,                                                       &
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

       ENDIF


      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('BL_INTCT ',4)
      ENDIF

      RETURN
      END SUBROUTINE BL_INTCT
#endif
