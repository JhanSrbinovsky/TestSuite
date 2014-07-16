#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE BDY_EXPL2----------------------------------------------
!!!
!!!  Purpose: Calculate the explicit turbulent fluxes of heat, moisture
!!!           and momentum between atmospheric levels
!!!           within the boundary layer, and/or the effects of these
!!!           fluxes on the primary model variables.
!!!
!!! Programming standard : unified model documentation paper No 3
!!!
!!!  Documentation: UMDP 24.
!!!
!!!---------------------------------------------------------------------

!    Arguments :-
      SUBROUTINE BDY_EXPL2 (                                            &

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
     & visc_bl_m,delta_lambda, delta_phi, FM_3D, FH_3D,L_subfilter_vert,&
     & L_subfilter_horiz,L_subfilter_blend,max_diff,turb_startlev_vert, &
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
     & )

      Use cv_run_mod, Only:                                             &
          l_conv4a, l_rp,l_rp2

      Use bl_diags_mod, Only:                                           &
          strnewbldiag

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
     &, n_proc                                                          &
                   ! Total number of processors
     &, land_pts   ! No.of land points in whole grid.


! (b) Defining vertical grid of model atmosphere.

      INTEGER                                                           &
     & BL_LEVELS                   ! IN Max. no. of "boundary" levels
!                                     allowed.Assumed <= 30 for dim-
!                                     sioning of GAMMA in common deck
!                                     C_GAMMA used in SF_EXCH and KMKH

      INTEGER, DIMENSION(20) :: BL_OPTIONS   ! IN BL switches


!     Declaration of new BL diagnostics.
      Type (Strnewbldiag) :: BL_diag

      REAL                                                              &
     &  r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j,bl_levels+1)                  &
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j,0:bl_levels+1)              &
     &, eta_theta_levels(0:bl_levels)                                   &
     &, P(1-off_x:row_length+off_x,                                     &
     &    1-off_y:rows+off_y, BL_LEVELS+1)                              &
     &, P_theta_levels(row_length, rows, BL_LEVELS+1)                   &
     &, U(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      BL_LEVELS)                                                  &
     &, V(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      BL_LEVELS)                                                  &
     &, W(row_length,rows, 0:BL_LEVELS)                                 &
     &, ETADOT(row_length,rows, 0:BL_LEVELS)                            &
     &,RHO_UV(row_length,rows,BL_LEVELS+1)                              &
!                                  ! IN density on UV (ie. rho) levels;
!                                  !    used in RHOKH so dry density if
!                                  !    Lq_mix_bl is true
     &,RHO_TQ(row_length,rows,BL_LEVELS)                                &
!                                  ! IN density on TQ (ie. theta) levels;
!                                  !    used in RHOKM so wet density
     &,RHO_DRY_TQ(row_length,rows,BL_LEVELS)                            &
!                                  ! IN density on TQ (ie. theta) levels;
!                                  !    used in non-turb flux integration
!                                  !    so dry density if Lq_mix_bl is true
     &,DZL_charney(row_length,rows,BL_LEVELS)                           &
                                   ! IN DZL(,K) is depth in m of theta
!                                  !    level K, i.e. distance from
!                                  !    boundary K-1/2 to boundary K+1/2
     &,RDZ(row_length,rows,BL_LEVELS)                                   &
                                   ! IN RDZ(,1) is the reciprocal of
!                                  !    the height of level 1, i.e. of
!                                  !    the middle of layer 1.  For
!                                  !    K > 1, RDZ(,K) is the
!                                  !    reciprocal of the vertical
!                                  !    distance from level K-1 to
!                                  !    level K.
     &,Z_FULL(row_length,rows,BL_LEVELS)                                &
                                   ! IN Z_FULL(*,K) is height of full
!                                  !    level k.
     &,Z_UV(row_length,rows,BL_LEVELS)                                  &
                                   ! IN Z_UV(*,K) is height of u level k
     &,Z_TQ(row_length,rows,BL_LEVELS)                                  &
                                   ! IN Z_TQ(*,K) is height of theta
!                                  !    level k.
     &,P_HALF(row_length,rows,BL_LEVELS)                                &
                                   ! IN P_HALF(*,K) is pressure at
!                                  !    half level k-1/2.
     &,DELTAP(row_length,rows,BL_LEVELS)                                &
                                   ! IN Difference in pressure between
!                                  !    levels
     &,RHOSTAR(ROW_LENGTH,ROWS)                                         &
                                   ! IN Surface air density
     &,U_P(row_length,rows,BL_LEVELS)                                   &
!                                  ! IN U on P-grid.
     &,V_P(row_length,rows,BL_LEVELS)                                   &
!                                  ! IN V on P-grid.
     &,BT(row_length,rows,BL_LEVELS)                                    &
                                   ! IN A buoyancy parameter for clear
!                                  !    air on p,T,q-levels
!                                  !    (full levels).
     &,BQ(row_length,rows,BL_LEVELS)                                    &
                                   ! IN A buoyancy parameter for clear
!                                  !    air on p,T,q-levels
!                                  !    (full levels).
     &,BT_CLD(row_length,rows,BL_LEVELS)                                &
!                                  ! IN A buoyancy parameter for cloudy
!                                  !    air on p,T,q-levels
!                                  !    (full levels).
     &,BQ_CLD(row_length,rows,BL_LEVELS)                                &
                                   ! IN A buoyancy parameter for cloudy
!                                  !    air on p,T,q-levels
!                                  !    (full levels).
     &,BT_GB(row_length,rows,BL_LEVELS)                                 &
                                   ! IN A grid-box mean buoyancy param
!                                  ! on p,T,q-levels (full levels).
     &,BQ_GB(row_length,rows,BL_LEVELS)                                 &
                                   ! IN A grid-box mean buoyancy param
!                                  ! on p,T,q-levels (full levels).
     &,A_QS(row_length,rows,BL_LEVELS)                                  &
                                   ! IN Saturated lapse rate factor
!                                  !    on p,T,q-levels (full levels).
     &,A_DQSDT(row_length,rows,BL_LEVELS)                               &
!                                  ! IN Saturated lapse rate factor
!                                  !    on p,T,q-levels (full levels).
     &,DQSDT(row_length,rows,BL_LEVELS)                                 &
                                   ! IN Derivative of q_SAT w.r.t. T
     &,visc_BL_m(1:row_length, 1:rows, bl_levels)
!         ! IN: lambda^2*S OUT: rho*lambda^2*S*FM (on BL_LEVELS)

      logical                                                           &
     & L_subfilter_vert                                                 &
!                         ! subgrid turbulence scheme in vertical        
     &,L_subfilter_horiz                                                &
!                         ! subgrid turbulence scheme in horizontal      
     &,L_subfilter_blend 
!                         ! blending BL and Smag coefficients

      Integer                                                           &
     & turb_startlev_vert                                               & 
!                         ! start and end vertical levels for          
     &,turb_endlev_vert  
!                         ! 3D turbulence scheme

      REAL, DIMENSION (:,:,:), ALLOCATABLE ::                           &
     & visc_BL_h                                                        &
                        ! OUT: lambda^2*S*FH (on BL_LEVELS)
     &,visc_BL_h_rho    ! visc_BL_h on rho levels

      Real                                                              &
     & delta_lambda                                                     &
     &,delta_phi                                                        &
     &,FM_3D(row_length,rows,BL_LEVELS)                                 &
!            ! stability function for momentum transport.
!            ! level 1 value is dummy for use in diagnostics
     &,FH_3D(row_length,rows,BL_LEVELS)                                 &
!             ! stability function for heat and moisture.
!             ! level 1 value is dummy for use in diagnostics
     &,max_diff                                                         &
                     ! maximum diffusion coefficient for this run.
     &,BL_COEF_KM(1:row_length, 1:rows, bl_levels-1)                    &
!            ! RHOKM from BL scheme
     &,BL_COEF_KH(1:row_length, 1:rows, bl_levels-1)                    &
!            ! RHOKH from BL scheme
     &,RECIP_L_MO_SEA(row_length,rows)                                  &
!                                  ! IN Reciprocal of the surface
!                                  !    Obukhov length at sea points.
!                                  !    (m-1).
     &,FLANDG(ROW_LENGTH,ROWS)
!                                  ! IN Land fraction on all tiles




! (f) Atmospheric + any other data not covered so far, incl control.

      REAL                                                              &
     & RAD_HR(row_length,rows,BL_LEVELS,2)                              &
                                    ! IN (LW,SW) rad heating rate (K/s)
     &, MICRO_TENDS(row_length, rows, bl_levels, 2)                     &
!                          ! Tendencies from microphys within BL levels
!                          ! (TL, K/s; QW, kg/kg/s)
     &,FB_SURF(row_length,rows)                                         &
                                    ! IN Surface flux buoyancy over
!                                   ! density (m^2/s^3)
!
     &,U_S(row_length,rows)                                             &
                                    ! IN Surface friction velocity
!                                   !    (m/s)
     &,T1_SD(row_length,rows)                                           &
                                    ! IN Standard deviation of
!                                   ! turbulent fluctuations of layer 1
!                                   ! temperature; for use in
!                                   ! initiating convection.
     &,Q1_SD(row_length,rows)                                           &
                                    ! IN Standard deviation of turbulent
!                                   !    fluctuations of layer 1
!                                   !    humidity; for use in initiating
!                                   !    convection.
     &,H_BLEND_OROG(row_length,rows)                                    &
                                    ! IN Blending height used as part
!                                   ! of effective roughness scheme
     &,Muw_SBL,Mwt_SBL                                                  &
                                    ! IN Powers to use in prescription
!                                   !    of equilibrium profiles of
!                                   !    stress and buoyancy flux in
!                                   !    Equilibrium SBL model
     &,TIMESTEP                                                         &
                                    ! IN Timestep (seconds).
     &,ZH_PREV(row_length, rows)                                        &
                                    ! IN boundary layer height from
!                                   !    previous timestep
     &,RIB_GB(row_length,rows)                                          &
                                 ! IN  Bulk Richardson number for lowest
!                                ! layer
     &,SIL_OROG_LAND(land_pts)   ! IN Silhouette area of unresolved
!                                ! orography per unit horizontal area

      INTEGER                                                           &
     & FORMDRAG                                                         &
                                 ! IN switch for orographic form drag
     &,FD_stab_dep               ! IN Switch to implement stability
!                                !    dependence of orog form drag

      REAL                                                              &
     & OROG_DRAG_PARAM                                                  &
!                                ! IN drag coefficient for form drag
     &,l_int                     ! Length scale for TKE diagnostic

      LOGICAL                                                           &
     & L_SBLeq                                                          &
                                 ! IN Switch for Equilibrium SBL model
     &,L_SBLco                                                          &
                                 ! IN Switch for coupled gradient
!                                !    method in Equilibrium SBL model
     &,L_LAMBDAM2                                                       &
                                 ! IN Operational settings for local
!                                !    mixing lengths: LambdaM=2*LambdaH
     &,L_FULL_LAMBDAS            !    and Lambdas NOT reduced above 
!                                !    NTML_LOCAL+1

! Additional variables for SCM diagnostics which are dummy in full UM
      INTEGER                                                           &
     & nSCMDpkgs             ! No of SCM diagnostics packages

      LOGICAL                                                           &
     & L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages
!
      LOGICAL                                                           &
     & lq_mix_bl
      REAL                                                              &
     & U_0_P(row_length,rows)                                           & 
                                   ! IN W'ly component of surface
!                                       current (m/s). P grid
     &,V_0_P(row_length,rows)                                           & 
                                   ! IN S'ly component of surface
!                                       current (m/s). P grid
     &,HO2R2_OROG(land_pts)        ! IN Standard Deviation of orography.
!                                     equivilent to peak to trough
!                                     height of unresolved orography
!                                     devided by 2SQRT(2) on land
!                                     points only (m)

      INTEGER                                                           &
     & LAND_INDEX(land_pts)        ! IN LAND_INDEX(I)=J => the Jth
!                                     point in P_FIELD is the Ith
!                                     land point.
      LOGICAL                                                           &
     & LAND_MASK(row_length,rows)  ! IN T if land, F elsewhere.

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

      LOGICAL LTIMER                ! Logical switch for TIMER diags


! INOUT variables
      REAL                                                              &
     & RDZ_CHARNEY_GRID(row_length,rows,BL_LEVELS)                      &
!                                  ! INOUT RDZ(,1) is the reciprocal of
!                                  !       the height of level 1,
!                                  !       i.e. of the middle of layer 1
!                                  !       For K > 1, RDZ(,K) is the
!                                  !       reciprocal of the vertical
!                                  !       distance from level K-1 to
!                                  !       level K.
     &,RDZ_U(row_length,rows,2:BL_LEVELS)                               &
                                   ! INOUT RDZ (K > 1) on UV-grid.
     &,RDZ_V(row_length,n_rows,2:BL_LEVELS)                             &
                                   ! INOUT RDZ (K > 1) on UV-grid.
     &,QW(row_length, rows, BL_LEVELS)                                  &
!                                  ! INOUT Total water content
     &,TL(row_length, rows, BL_LEVELS)                                  &
!                                  ! INOUT Ice/liquid water temperature
     &,ZH(row_length,rows)         ! INOUT Height above surface of top
!                                  !       of boundary layer (metres).

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

! Definition of variables for STPH_RP
      REAL,INTENT(InOut) :: par_mezcla ! Used to modify LAMBDAH,LAMBDAM
                                       !(neutral mixing length in EXCOEF
      REAL,INTENT(InOut) :: G0_RP ! Stability function parameter in EXCOEF


!  Outputs :-
!-1 Diagnostic (or effectively so - includes coupled model requisites):-

!  (a) Calculated anyway (use STASH space from higher level) :-
!
      REAL                                                              &
     & FQW(row_length,rows,BL_LEVELS)                                   &
                                   ! OUT Moisture flux between layers
!                                     (kg per square metre per sec).
!                                     FQW(,1) is total water flux
!                                     from surface, 'E'.
     &,FTL(row_length,rows,BL_LEVELS)                                   &
                                   ! OUT FTL(,K) contains net turbulent
!                                     sensible heat flux into layer K
!                                     from below; so FTL(,1) is the
!                                     surface sensible heat, H. (W/m2)
     &,RHOKH(row_length,rows,BL_LEVELS)                                 &
                                   ! OUT Exchange coeffs for moisture.
     &,RHOKM_U(row_length,rows,BL_LEVELS)                               &
                                   ! OUT Exchange coefficients for u
     &,RHOKM_V(row_length,n_rows,BL_LEVELS)                             &
                                   ! OUT Exchange coefficients for v
     &,TAUX(row_length,rows,BL_LEVELS)                                  &
                                   ! OUT W'ly component of surface wind
!                                     stress (N/sq m).(On UV-grid with
!                                     first and last rows undefined or
!                                     at present, set to missing data
     &,TAUY(row_length,n_rows,BL_LEVELS)                                &
                                   ! OUT S'ly component of surface wind
!                                     stress (N/sq m).  On UV-grid;
!                                     comments as per TAUX.
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
     &,RHO_KM(row_length,rows,2:BL_LEVELS)
                                   ! OUT Air density * turbulent mixing
!                                     coefficient for momentum before
!                                     interpolation.
      REAL                                                              &
     & TAUX_FD_U(row_length,rows,BL_LEVELS)                             &
!                                  ! OUT  X comp. of orographic stress
!                                  ! interpolated to U points
     &,TAUY_FD_V(row_length,n_rows,BL_LEVELS)
!                                  ! OUT  Y comp. of orographic stress
!                                  ! interpolated to V points

!*APL*DIAGS
      REAL                                                              &
     &  ZHT(row_length, rows)                                           &
                                   ! OUT Max height of turb mixing
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
     &, ZHSC(row_length,rows)       ! OUT Top of decoupled layer

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

!-2 Genuinely output, needed by other atmospheric routines :-

      REAL                                                              &
     & Z0M_EFF_GB(row_length,rows) ! INOUT Effective grid-box roughness
!                                 length for momentum

      Real                                                              &
     &  uw0(row_length,rows)                                            &
                             ! OUT U-component of surface wind stress
!                            !     on P-grid
     &, vw0(row_length,rows) ! OUT V-component of surface wind stress
!                            !     on P-grid

!---------------------------------------------------------------------
!  External routines called :-

      EXTERNAL BTQ_INT,KMKHZ,EX_COEF,KMKH,Swap_Bounds,                  &
     & FILL_EXTERNAL_HALOS,P_TO_U,P_TO_V,EX_FLUX_TQ,EX_FLUX_UV,         &
     & FM_DRAG
      EXTERNAL TIMER

!-----------------------------------------------------------------------
!   Symbolic constants (parameters) reqd in top-level routine :-

#include "c_r_cp.h"
#include "c_g.h"
#include "c_lheat.h"
#include "soil_thick.h"
#include "c_vkman.h"
#include "fldtype.h"
#include "domtyp.h"

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

#include "blopt8a.h"
! Parameters also passed to EX_COEF
      REAL       LAMBDA_MIN  ! Min value of length scale LAMBDA.
!!      No longer a parameter, in order to vary it!
!!      PARAMETER (LAMBDA_MIN=40.0)
! Layer interface K_LOG_LAYR-1/2 is the highest which requires log
! profile correction factors to the vertical finite differences.
! The value should be reassessed if the vertical resolution is changed.
! We could set K_LOG_LAYR = BL_LEVELS and thus apply the correction
! factors for all the interfaces treated by the boundary layer scheme;
! this would be desirable theoretically but expensive computationally
! because of the use of the log function.
      INTEGER    K_LOG_LAYR
      PARAMETER (K_LOG_LAYR=2)
!-----------------------------------------------------------------------

!  Workspace :-

      REAL                                                              &
     & A_DQSDTM(row_length,rows,BL_LEVELS)                              &
!                               ! Saturated lapse rate factor
!                               ! on intermediate levels (half levels).
     &,A_QSM(row_length,rows,BL_LEVELS)                                 &
!                               ! Saturated lapse rate factor
!                               ! on intermediate levels (half levels).
     &,BQM(row_length,rows,BL_LEVELS)                                   &
                                ! A buoyancy parameter for clear air
!                               ! on intermediate levels (half levels).
     &,BQM_CLD(row_length,rows,BL_LEVELS)                               &
!                               ! A buoyancy parameter for cloudy air
!                               ! on intermediate levels (half levels).
     &,BTM(row_length,rows,BL_LEVELS)                                   &
                                ! A buoyancy parameter for clear air
!                               ! on intermediate levels (half levels).
     &,BTM_CLD(row_length,rows,BL_LEVELS)                               &
!                               ! A buoyancy parameter for cloudy air
!                               ! on intermediate levels (half levels).
     &,CFM(row_length,rows,BL_LEVELS)                                   &
!                               ! Estimate of cloud fraction
!                               ! on intermediate levels (half levels).
     &,DBDZ(row_length,rows,2:BL_LEVELS)                                &
!                               ! Buoyancy gradient across layer
                                !  interface.
     &,DVDZM(row_length,rows,2:BL_LEVELS)                               &
                                         ! Modulus of wind shear.
     &,DELTAP_UV(row_length,rows,BL_LEVELS)                             &
!                                 Difference in pressure between levels
!                                 on UV points
     &,DQW_1(row_length,rows)                                           &
                                ! Increment for QW(,1).
     &,F_NGSTRESS_HALO(1-off_x:row_length+off_x,1-off_y:rows+off_y,     &
     &                                                  2:BL_LEVELS)    &
                                ! dimensionless function for
!                               ! non-gradient stresses
     &,F_NGSTRESS(row_length,rows,2:BL_LEVELS)                          &
                                ! dimensionless function for
!                               ! non-gradient stresses
     &,F_NGSTRESS_U(row_length,rows,2:BL_LEVELS)                        &
                                ! dimensionless function for
!                               ! non-gradient U stress
     &,F_NGSTRESS_V(row_length,n_rows,2:BL_LEVELS)                      &
                                ! dimensionless function for
!                               ! non-gradient V stress
     &,GRAD_Q_ADJ(row_length,rows)                                      &
                                ! Humidity gradient adjustment
!                                 for non-local mixing in unstable
!                                 turbulent boundary layer.
     &,GRAD_T_ADJ(row_length,rows)                                      &
                                   ! Temperature gradient adjustment
!                                 for non-local mixing in unstable
!                                 turbulent boundary layer.
     &,RHOKHZ(row_length,rows,2:BL_LEVELS)                              &
!                               ! Non-local turbulent mixing
!                                 coefficient for heat and moisture.
     &,RHOKH_TOP(row_length,rows,2:BL_LEVELS)                           &
!                               ! Non-local turbulent mixing coefficient
!                               ! for top-down mixing of heat and
!                               ! moisture.
     &,RHOKH_uv(row_length,rows,BL_LEVELS)                              &
     &,RHOKM(1-off_x:row_length+off_x,1-off_y:rows+off_y,BL_LEVELS)     &
!            Exchange coefficients for momentum on P-grid
     &,RHOKMZ(row_length,rows,2:BL_LEVELS)                              &
!                               ! Non-local turbulent mixing
!                                 coefficient for momentum.
     &,RHOKM_TOP(row_length,rows,2:BL_LEVELS)                           &
!                               ! Non-local turbulent mixing coefficient
!                               ! for top-down mixing of momentum.
     &,ZLB(row_length,rows,0:BL_LEVELS)                                 &
                                ! ZLB(,K) is the height of
!                                 theta level K ( = 0.0 for "K=0")
     &,SIGMA_H(row_length,rows)
                                ! Standard deviation of subgrid 
                                ! orography (m) [= 2root2 * ho2r2_orog]
!
!     ! Terms for non-gradient flux parametrization
!     !  (=0 unless using 8C code with FLUX_GRAD=LockWhelan2006)
!
      REAL, DIMENSION(row_length,rows, BL_LEVELS+1) ::                  &
     &  FT_NT                                                           &
                                ! Non-turbulent heat and moisture flux
     &, FQ_NT                   !  (on rho levels, surface flux(K=1)=0)
      REAL, DIMENSION(row_length,rows, 2:BL_LEVELS) ::                  &
     &  RHOF2                                                           &
                                ! f2 and fsc term shape profiles
     &, RHOFSC                  !

      REAL, DIMENSION(row_length,rows) ::                               &
     &  TOTHF_ZH                                                        &
                                ! Total heat fluxes at inversions
     &, TOTHF_ZHSC                                                      &
                                !
     &, TOTQF_ZH                                                        &
                                ! Total moisture fluxes at inversions
     &, TOTQF_ZHSC                                                      &
                                !
     &, FT_NT_DSCB                                                      &
                                ! Non-turbulent heat and moisture flux
     &, FQ_NT_DSCB              !    at the base of the DSC layer.

       REAL                                                             &
     & ZH_LOCAL(row_length,rows)                                        &
                                   ! Height above surface of top of
!                                  !  boundary layer (metres) as
!                                  !  determined from the local
!                                  !  Richardson number profile.
     &,ZC_OUT(row_length, rows)                                         &
                                ! Sc cloud-depth
     &,LWP(row_length, rows)                                            &
                                ! Integrated (QCL+QCF)*DZ
     &,DTLDZ(row_length,rows,2:BL_LEVELS)                               &
!                                  ! TL+gz/cp gradient between
!                                  ! levels K and K-1
     &,DQWDZ(row_length,rows,2:BL_LEVELS)                               &
!                                  ! QW gradient between
     &,TAU_FD_X(row_length,rows,BL_LEVELS)                              &
                                           ! X comp of orographic stress
     &,TAU_FD_Y(row_length,rows,BL_LEVELS)                              &
                                           ! Y comp of orographic stress
     &,TAUX_HAL(1-off_x:row_length+off_x,1-off_y:rows+off_y,BL_LEVELS)  &
!                                          ! X comp of orographic stress
     &,TAUY_HAL(1-off_x:row_length+off_x,1-off_y:rows+off_y,BL_LEVELS)
!                                          ! Y comp of orographic stress
!
      INTEGER                                                           &
     & NTML_LOCAL(row_length,rows)                                      &
                                   ! Number of model layers in the
!                                    turbulently mixed layer as
!                                    determined from the local
!                                    Richardson number profile.
     &,NTML_NL(row_length,rows)                                         &
                                   ! Number of model layers in the
!                                    turbulently mixed layer as
!                                    determined from the parcel ascent.
     &,SML_DISC_INV(row_length,rows)                                    &
                                    ! Flags for whether discontinuous
     &,DSC_DISC_INV(row_length,rows)! inversions are diagnosed

      LOGICAL                                                           &
     & UNSTABLE(row_length,rows)                                        &
                                 ! Logical switch for unstable
!                                !    surface layer.
     &,DSC(row_length,rows)                                             &
                                 ! Flag set if decoupled
!                                ! stratocumulus layer found
     &,COUPLED(row_length,rows)  ! Flag to indicate Sc layer weakly
!                                ! coupled to surface (ie weakly
                                 ! decoupled)

#if defined(SCMA)
! INOUT SCMop is declared in here
#include "s_scmop.h"
#endif
       REAL                                                             &
     & SL(row_length,rows,bl_levels)                                    &
                                                  ! Static energy
     &,TAU_GRAD_U(row_length,rows,bl_levels)                            &
                                                  ! K*du/dz
     &,TAU_NON_GRAD_U(row_length,rows,bl_levels)                        &
                                                  ! Non-grad stress
     &,TAU_GRAD_V(row_length,n_rows,bl_levels)                          &
                                                  ! K*du/dz
     &,TAU_NON_GRAD_V(row_length,n_rows,bl_levels)! Non-grad stress

!
!  Local scalars :-

      Character*(*), Parameter ::  RoutineName = 'bdy_expl2'

      REAL                                                              &
     & DZU                                                              &
              ! Westerly wind shear between levels K+1 and K.
     &,DZV                                                              &
              ! Southerly wind shear between levels K+1 and K.
     &,ELH                                                              &
              ! Mixing length for heat & moisture at lower layer bdy.
     &,LAMBDAH                                                          &
              ! Asymptotic mixing length for turbulent transport
              ! of heat/moisture.
     &,VKZ                                                              &
              ! Temporary in calculation of ELH.
     &,F_LOG                                                            &
              ! Temporary in calculation of logarithmic correction
     &,ZMAX_FOR_DSC
              ! Maximum height to look for DSC cloud-base

      REAL                                                              &
     &  WEIGHT1                                                         &
     &, WEIGHT2                                                         &
     &, WEIGHT3                                                         &
     &, Z_SCALE                                                         &
     &, GRCP                                                            &
               ! G/CP
     &, DTLDZM                                                          &
               ! TL+gz/cp gradient interpolated to Z_TQ
     &, DQWDZM ! QW gradient interpolated to Z_TQ

      INTEGER                                                           &
     & I,J                                                              &
                  ! LOCAL Loop counter (horizontal field index).
     &,K,IENT                                                           &
                  ! LOCAL Loop counter (vertical level index).
     &,KP,KM                                                            &
                  ! K+/-1,
     &,L          ! LOCAL Loop counter for land points                  

      INTEGER                                                           &
     & SBL_OP                   !  stable boundary layer option

!
! Local BL swtiches, passed down in BL_OPTIONS
!
      INTEGER                                                           &
     & NG_STRESS                                                        &
                   ! switch for non-gradient stress
     &,ISHEAR_BL                                                        &
                   ! switch for shear-dominated b.l.
     &,DECFIX                                                           &
                   ! switch for correction to decoupling diagnosis
     &,STOPWE_SBL                                                       &
                   ! switch for spurious entrainment in SBLs
     &,FLUX_GRAD                                                        &
                   ! switch for revised flux-gradient relationships
     &,NON_LOCAL_BL                                                     &
                   ! switch on the non-local scheme
     &,LOCAL_FA                                                         &
                   ! switch for free atmospheric mixing options
     &,Keep_Ri_FA                                                       &
                   ! switch to keep local mixing in the free atmosphere
     &,Prandtl                                                          &
                   ! switch for Prandtl number options
     &,NL_BL_LEVELS
                   ! number of levels for the non-local scheme

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('BDY_EXPL2 ',3)
      END IF

!-----------------------------------------------------------------------
! Convert BL options to local variables
!-----------------------------------------------------------------------
      ISHEAR_BL   = BL_OPTIONS(1)
      NG_STRESS   = BL_OPTIONS(2)
      DECFIX      = BL_OPTIONS(3)
      STOPWE_SBL  = BL_OPTIONS(4)
!     TRWEIGHTS1  = BL_OPTIONS(5)   ! only needed in IMP_CTL
      FLUX_GRAD   = BL_OPTIONS(6)
      NON_LOCAL_BL= BL_OPTIONS(11)
      LOCAL_FA    = BL_OPTIONS(12)
      Prandtl     = BL_OPTIONS(13)
      Keep_Ri_FA  = BL_OPTIONS(17)
      NL_BL_LEVELS= BL_OPTIONS(18)

      IF ( NL_BL_LEVELS == OFF .OR.                                     &
     &     NL_BL_LEVELS > BL_LEVELS ) NL_BL_LEVELS = BL_LEVELS
                                      ! if unset or too large 
!-----------------------------------------------------------------------

      IF (FORMDRAG ==  Explicit_stress) THEN

!------------------------------------------------------------------
!      Set stresses to zero
!------------------------------------------------------------------

       DO K=1,BL_LEVELS
        DO J = 1, ROWS
         DO I=1, ROW_LENGTH
           TAU_FD_X(I,J,K) = 0.0
           TAU_FD_Y(I,J,K) = 0.0
         END DO
        END DO
       END DO

!------------------------------------------------------------------
!      Calculate stress profiles
!------------------------------------------------------------------
! DEPENDS ON: fm_drag
       CALL FM_DRAG ( ROW_LENGTH, ROWS                                  &
     &,  LAND_PTS, LAND_INDEX, FD_stab_dep                              &
     &,  BL_LEVELS, U_P, V_P, RHO_TQ, Z_UV, Z_TQ, Z0M_EFF_GB            &
     &,  TAU_FD_X, TAU_FD_Y, ZH_PREV, RIB_GB, SIL_OROG_LAND             &
     &,  OROG_DRAG_PARAM, LTIMER,  BL_diag )

!------------------------------------------------------------------
!      Orographic stress diagnostics
!------------------------------------------------------------------
       DO K=1,BL_LEVELS
        DO J = 1, ROWS
         DO I=1, ROW_LENGTH
           IF (BL_diag%L_ostressx)                                      &
     &         BL_diag%ostressx(i,j,k)=TAU_FD_X(i,j,k)
           IF (BL_diag%L_ostressy)                                      &
     &         BL_diag%ostressy(i,j,k)=TAU_FD_Y(i,j,k)
         END DO
        END DO
       END DO

      END IF

!
!-----------------------------------------------------------------------
! Apply dynamic diagnosis of shear-driven layers.
!-----------------------------------------------------------------------
!
      IF (BL_OPTIONS(9) == DynDiag_ZL) THEN
!
!       In cases where the parcel ascent continues right through the
!       boundary layer, we diagnose cumulus only where the surface
!       buoyancy flux is sufficiently unstable, testing the ratio of
!       the depth of the inversion to the Obukhov length. A value of
!       1 -- 2 is reasonable for this test and 1.6 is selected, but
!       no great precision is attached to this value. Since this is
!       of importance mainly at sea-points, to avoid complications
!       with coastal tiling, the scheme operates only at points
!       where the land fraction is below 0.5.
!
        DO j = 1, rows
          DO I=1,row_length
!
            IF ( FLANDG(I,j) < 0.5 ) THEN
              IF ( -Z_HALF(I,j,MIN(NTPAR(I,j),BL_LEVELS-1)+1) *         &
     &              RECIP_L_MO_SEA(I,j) < 1.6 ) THEN
                CUMULUS(I,j)   = .FALSE.
                L_SHALLOW(I,j) = .FALSE.
                NTML(I,j)      = MIN(NTPAR(I,j), BL_LEVELS-1)
                ZH(I,j)        = Z_HALF(I,j,NTML(I,j)+1)
              END IF
!             Override the provisional cumulus diagnosis if the
!             actual surface buoyancy flux indicates stability.
              IF ( FB_SURF(I,j) < 0.0 ) THEN
                CUMULUS(I,j)   = .FALSE.
                L_SHALLOW(I,j) = .FALSE.
                NTML(I,j)      = 1
                ZH(I,j)        = Z_HALF(I,j,2)
              END IF
            END IF
!
          END DO
        END DO
!
      END IF
!
!-----------------------------------------------------------------------
! Set NTML_NL to NTML as passed in from initial diagnosis routine
!-----------------------------------------------------------------------

      DO j = 1, rows
      DO I=1,row_length
        NTML_NL(I,j) = NTML(I,j)
        ZMAX_FOR_DSC = 2500.0
        IF ( NL_BL_LEVELS < BL_LEVELS .AND.                             &
     &       NTML_NL(I,j) > NL_BL_LEVELS-1 ) THEN
           NTML_NL(I,j) = NL_BL_LEVELS-1
           ZH(I,j)      = Z_HALF(I,j,NTML_NL(I,j)+1)
           ZMAX_FOR_DSC = Z_HALF(I,j,NL_BL_LEVELS-1)
        END IF
      END DO
      END DO

!-----------------------------------------------------------------------
!! 5.  Turbulent exchange coefficients and "explicit" fluxes between
!!     model layers in the boundary layer (P243b, routine KMKH).
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!! 5.2  Interpolate BT and BQ to half levels.
!-----------------------------------------------------------------------

! DEPENDS ON: btq_int
      CALL BTQ_INT (                                                    &
     & row_length,rows,BL_LEVELS                                        &
     &,z_full,z_half,BQ,BT,BQ_CLD,BT_CLD,A_QS,A_DQSDT                   &
     &,BQM,BTM,BQM_CLD,BTM_CLD,A_QSM,A_DQSDTM                           &
     &,CF,CFM                                                           &
     &,LTIMER                                                           &
     &  )

!-----------------------------------------------------------------------
!! 5.3  Calculate the diffusion coefficients Km and Kh.
!-----------------------------------------------------------------------
!     ! Initialise non-local K and fluxes to zero; necessary for levels 
!     ! above NL_BL_LEVELS
      DO K =2, BL_LEVELS
        DO j = 1, rows
        DO I = 1, row_length 
          FTL(I,j,K) = 0.0
          FQW(I,j,K) = 0.0
          RHOKMZ(I,j,K) = 0.0
          RHOKHZ(I,j,K) = 0.0
          RHOKM_TOP(I,j,K) = 0.0
          RHOKH_TOP(I,j,K) = 0.0
          F_NGSTRESS(I,J,K) = 0.0
        END DO
        END DO
      END DO
!     ! Initialise Lock-Whelan non-gradient terms to zero
!     ! Only calculated in KMKHZ8C
      DO K=1,BL_LEVELS
        DO j=1, rows
        DO I=1, row_length
          FT_NT(I,J,K)  = 0.0
          FQ_NT(I,J,K)  = 0.0
        END DO
        END DO
      END DO
      DO K=2,BL_LEVELS
        DO j=1, rows
        DO I=1, row_length
          RHOF2(I,J,K)  = 0.0
          RHOFSC(I,J,K) = 0.0
        END DO
        END DO
      END DO
      DO j=1, rows
      DO I=1, row_length
        TOTHF_ZH(I,J)   = 0.0
        TOTHF_ZHSC(I,J) = 0.0
        TOTQF_ZH(I,J)   = 0.0
        TOTQF_ZHSC(I,J) = 0.0
        FT_NT_DSCB(I,J) = 0.0
        FQ_NT_DSCB(I,J) = 0.0
      END DO
      END DO

      IF (L_subfilter_vert .AND. .NOT. L_subfilter_blend) THEN
        NON_LOCAL_BL = OFF
      END IF

      IF (NON_LOCAL_BL == ON) THEN

! DEPENDS ON: kmkhz
        CALL KMKHZ (                                                    &
     &     row_length,rows,NL_BL_LEVELS,halo_i,halo_j,n_proc,           &
     &     DECFIX, STOPWE_SBL, NG_STRESS, FLUX_GRAD, lq_mix_bl,         &
     &     TIMESTEP,P_theta_levels,P_HALF,RHO_TQ,RHO_UV,RHO_DRY_TQ,     &
     &     T,Q,QCL,QCF,CF,QW,TL,                                        &
     &     DZL_charney,RDZ_charney_grid,Z_FULL,Z_HALF,Z_UV,Z_TQ,        &
     &     FTL,FQW,RHOKMZ,RHOKHZ,RHOKM_TOP,RHOKH_TOP,                   &
     &     W,ETADOT,ETA_THETA_LEVELS,RAD_HR,MICRO_TENDS,                &
     &     BT,BQ,BTM,BQM,DQSDT,BTM_CLD,BQM_CLD,A_QS,A_QSM,A_DQSDTM,CFM, &
     &     Z0M_EFF_GB,U_S,FB_SURF,T1_SD,Q1_SD,RHOSTAR,                  &
     &     UNSTABLE,DSC,COUPLED,SML_DISC_INV,DSC_DISC_INV,              &
     &     CUMULUS,L_SHALLOW,NTML_NL,NTPAR,NLCL,NTDSC,NBDSC,            &
     &     ZH,ZH_PREV,ZHPAR,Z_LCL,ZHSC,ZMAX_FOR_DSC,                    &
     &     F_NGSTRESS, GRAD_T_ADJ, GRAD_Q_ADJ,                          &
     &     RHOF2, RHOFSC, FT_NT, FQ_NT, FT_NT_DSCB, FQ_NT_DSCB,         &
     &     TOTHF_ZH, TOTHF_ZHSC, TOTQF_ZH, TOTQF_ZHSC,                  &
     &     KENT, WE_LIM, T_FRAC, ZRZI,                                  &
     &     KENT_DSC, WE_LIM_DSC, T_FRAC_DSC, ZRZI_DSC,                  &
     &     BL_diag, nSCMDpkgs,L_SCMDiags,                               &
     &     LTIMER                                                       &
     &     )

      ELSE   ! not NON_LOCAL_BL

!        !-------------------------------------------------------------
!        ! Set all variables from the non-local scheme to zero or "off"
!        !  - reset all fluxes and K's arising from the non-local scheme
!        !-------------------------------------------------------------

        DO j = 1, rows
        DO I = 1, row_length 
!         ! surface mixed layer
          CUMULUS(I,j) = .FALSE.
          L_SHALLOW(I,j) = .FALSE.
          SML_DISC_INV(I,j) = 0
          NTPAR(I,j)   = 0
          NTML_NL(I,j) = -1    ! to ensure correct diagnostics
          ZH(I,j)      = 0.0
          GRAD_T_ADJ(I,j) = 0.0
          GRAD_Q_ADJ(I,j) = 0.0
!         ! decoupled mixed layer
          DSC(I,j)     = .FALSE.
          DSC_DISC_INV(I,j) = 0
          NTDSC(I,j)   = 0
          NBDSC(I,j)   = 0
          ZHSC(I,j)    = 0.0
          COUPLED(I,j) = .FALSE.
!         ! entrainment variables for non-local tracer mixing
          KENT(I,j) = 2
          KENT_DSC(I,j) = 2
          DO IENT = 1,3
           T_FRAC(I,j,IENT) = 0.0
           ZRZI(I,j,IENT)   = 0.0
           WE_LIM(I,j,IENT) = 0.0
           T_FRAC_DSC(I,j,IENT) = 0.0
           ZRZI_DSC(I,j,IENT)   = 0.0
           WE_LIM_DSC(I,j,IENT) = 0.0
          END DO
        END DO
        END DO

      END IF  ! test on NON_LOCAL_BL
!-----------------------------------------------------------------------
!! Calculate lapse rates
!-----------------------------------------------------------------------
      GRCP = G/CP
      DO K =2, BL_LEVELS
        DO j = 1, rows
          DO I = 1, row_length
            DTLDZ(I,j,K) = ( TL(I,j,K) - TL(I,j,K-1) )                  &
     &                             * RDZ_charney_grid(I,J,K) + GRCP
            DQWDZ(I,j,K) = ( QW(I,j,K) - QW(I,j,K-1) )                  &
     &                             * RDZ_charney_grid(I,J,K)
          END DO
        END DO
      END DO
!
#if defined(SCMA)
!
!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
      IF (L_SCMDiags(SCMDiag_bl)) THEN
!
        DO K=1,bl_levels
          DO j=1, rows
            DO I=1, row_length
              SL(I,j,K) = TL(I,j,K) + GRCP*Z_FULL(I,j,K)
            END DO ! I
          END DO ! j
        END DO ! K
!
!       Output SL
!
! DEPENDS ON: scmoutput
        Call SCMoutput(SL,                                              &
             'SL','Liquid/frozen water static energy (IN)','K',         &
             t_avg, d_bl, default_streams, '',RoutineName)
!
!       Output QW
!
! DEPENDS ON: scmoutput
        Call SCMoutput(qw,                                              &
             'qw','Total water content (IN)','kg/kg',                   &
             t_avg, d_bl, default_streams, '',RoutineName)
!
      END IF ! L_SCMDiags(SCMDiag_bl)
#endif


! Local Ri-based calculation of RHOKM and RHOKH:
! Calculate `buoyancy' gradient, DBDZ, on theta-levels
! NOTE: DBDZ(K) is on theta-level K-1
      DO K =3, BL_LEVELS
        DO j = 1, rows
          DO I = 1, row_length
            WEIGHT1 = r_rho_levels(I,J,K) -                             &
     &                r_rho_levels(I,J,K-1)
            WEIGHT2 = r_theta_levels(I,J,K-1)-                          &
     &                r_rho_levels(I,J,K-1)
            WEIGHT3 = r_rho_levels(I,J,K) -                             &
     &                r_theta_levels(I,J,K-1)
            DTLDZM = WEIGHT2 * DTLDZ(I,J,K)                             &
     &             + WEIGHT3 * DTLDZ(I,J,K-1)
            DQWDZM = WEIGHT2 * DQWDZ(I,J,K)                             &
     &             + WEIGHT3 * DQWDZ(I,J,K-1)

            DBDZ(I,J,K) = G*( BT_GB(I,J,K-1)*DTLDZM +                   &
     &                        BQ_GB(I,J,K-1)*DQWDZM )/WEIGHT1

          END DO
        END DO
      END DO

      K = 2
      DO j = 1, rows
        DO I = 1, row_length
          DBDZ(I,J,K) = G*( BT_GB(I,J,K-1)*DTLDZ(I,J,K) +               &
     &                      BQ_GB(I,J,K-1)*DQWDZ(I,J,K) )
        END DO
      END DO
! calculate modulus of shear on theta-levels
      DO K =2, BL_LEVELS
        DO j = 1, rows
          DO I = 1, row_length
            DZU = U_P(I,j,K) - U_P(I,j,K-1)
            DZV = V_P(I,j,K) - V_P(I,j,K-1)
            DVDZM(I,J,K) = MAX ( 1.0E-12 ,                              &
     &               SQRT(DZU*DZU + DZV*DZV) * RDZ(I,j,K)  )
          END DO
        END DO
      END DO
!-----------------------------------------------------------------------
!!  Calculate standard deviation of subgrid orography.  Potentially 
!!  used to weight between stable stability functions in EX_COEF.
!!  Also set LAMBDA_MIN (larger under LOCAL_FA option 3)
!-----------------------------------------------------------------------
      LAMBDA_MIN=40.0
      IF (LOCAL_FA == 2 .OR. LOCAL_FA == 3) THEN
!       ! link scale height to subgrid orography
        do j=1,rows
        DO I=1,row_length
          SIGMA_H(I,J) = 0.0
        ENDDO
        ENDDO
        DO L=1,LAND_PTS
          J=(land_index(L)-1)/row_length + 1
          I=land_index(L) - (J-1)*row_length
!         ! Multiply by 2root2 to get back to sigma_h
          SIGMA_H(I,J) = 2.83*HO2R2_OROG(L)
        ENDDO
        IF (LOCAL_FA == 3) LAMBDA_MIN=150.0
      END IF
!-----------------------------------------------------------------------
! call local coeff calculation for levels 2 to bl_levels
!-----------------------------------------------------------------------
      SBL_OP=BL_OPTIONS(7)

! call local coeff calculation for levels 2 to bl_levels

! DEPENDS ON: ex_coef
      CALL EX_COEF (                                                    &
     & row_length,rows,off_x,off_y,BL_LEVELS                            &
     &,K_LOG_LAYR,LAMBDA_MIN,L_LAMBDAM2,L_FULL_LAMBDAS,lq_mix_bl        &
     &,LOCAL_FA,Prandtl,SIGMA_H,CCB,CCT,NTML_LOCAL,ISHEAR_BL            &
     &,L_SBLeq,L_SBLco,Muw_SBL,Mwt_SBL                                  &
     &,SBL_OP,FLANDG                                                    &
     &,CCA,DBDZ,DVDZM                                                   &
     &,RHO_TQ,ZH_LOCAL,Z_UV,Z_TQ,Z0M_EFF_GB,H_BLEND_OROG                &
     &,CUMULUS,NTPAR,NTML_NL                                            &
     &,U_P,V_P,U_S,FB_SURF,QW,TL                                        &
     &,FM_3D, FH_3D, L_subfilter_vert, L_subfilter_horiz                &
     &,RHOKM,RHOKH_uv                                                   &
! Variables for STPH_RP 
     &,G0_RP,par_mezcla                                                 &
     &,LTIMER,BL_diag                                                   &
     &,nSCMDpkgs,L_SCMDiags                                             &
     & )

! interpolate rhokh_uv to rhokh levels 2 to bl_levels
      DO K =2, BL_LEVELS
        DO j = 1, rows
          DO I = 1, row_length
! potential change in order to preserve inversion structure
!            IF ( K  ==  NTML_LOCAL(I,J) ) THEN
!              RHOKH(I,J,K) = RHOKH_uv(I,J,K+1)
!            ELSE
              WEIGHT1 = r_theta_levels(I,J,K) -                         &
     &                r_theta_levels(I,J, K-1)
              WEIGHT2 = r_theta_levels(I,J,K) -                         &
     &                r_rho_levels(I,J,K)
              WEIGHT3 = r_rho_levels(I,J,K) -                           &
     &                r_theta_levels(I,J,K-1)
              IF ( K  ==  BL_LEVELS ) THEN
!               ! assume RHOKH_uv(BL_LEVELS+1) is zero
                RHOKH(I,J,K) = ( WEIGHT2/WEIGHT1 ) * RHOKH_uv(I,J,K)
              ELSE
                RHOKH(I,J,K) =                                          &
     &              WEIGHT3/WEIGHT1 *                                   &
     &                      RHOKH_uv(I,J,K+1)                           &
     &             +WEIGHT2/WEIGHT1 *                                   &
     &                      RHOKH_uv(I,J,K)
              END IF
!            END IF
            IF ((.NOT.L_SBLeq).OR.(FB_SURF(I,j) >  0.0)) THEN
!-------------------------------------------------------------------
!  Code moved from EX_COEF to avoid interpolation:
!  Include mixing length, ELH, in RHOKH.
!  Also incorporate log profile correction to the vertical finite
!  differences.
!-------------------------------------------------------------------
             IF (L_RP .or. L_RP2) THEN
               LAMBDAH = MAX ( LAMBDA_MIN , par_mezcla*ZH_LOCAL(I,j) )
             ELSE
               LAMBDAH = MAX ( LAMBDA_MIN , 0.15*ZH_LOCAL(I,j) )
             END IF
!-----------------------------------------------------------------------
!            ! Optionally reduce mixing length above local BL top
!-----------------------------------------------------------------------
             IF (K >= NTML_LOCAL(I,j)+2 .and. .NOT.L_FULL_LAMBDAS) THEN
               LAMBDAH = LAMBDA_MIN
             END IF
             IF ( K >= NTML_LOCAL(I,j)+2 .and. L_FULL_LAMBDAS .and.     &
     &            LOCAL_FA == 1) THEN
!              ! Assuming only LOCAL_FA = 1 will have L_FULL_LAMBDAS
!              ! If other LOCAL_FA options are coded here then 
!              ! changes must be included in section 2.1 of ex_coef
               Z_SCALE = 1000.0
               WEIGHT1 = 0.5*( 1.0 -                                    & 
     &                     TANH(3.*((Z_UV(I,j,K)/Z_SCALE )-1.0) ) ) 
               LAMBDAH = LAMBDAH * WEIGHT1                              &
     &                      + LAMBDA_MIN*( 1.0 -  WEIGHT1) 

             END IF

             IF (K  <=  K_LOG_LAYR) THEN
               VKZ   = VKMAN * ( Z_TQ(I,j,K) - Z_TQ(I,j,K-1) )
               F_LOG = LOG( ( Z_TQ(I,j,K) + Z0M_EFF_GB(I,j)   ) /       &
     &                      ( Z_TQ(I,j,K-1) + Z0M_EFF_GB(I,j) ) )
               ELH   = VKZ / ( F_LOG + VKZ / LAMBDAH )
             ELSE
               VKZ = VKMAN * ( Z_UV(I,j,K) + Z0M_EFF_GB(I,j) )
               ELH = VKZ / (1.0 + VKZ/LAMBDAH )
             END IF

             RHOKH(I,J,K) = ELH * RHOKH(I,J,K)

            END IF   ! test on L_SBLeq

!           ! Finally multiply RHOKH by dry density
            If (Lq_mix_bl)                                              & 
     &         RHOKH(I,J,K) = RHO_UV(I,j,K) * RHOKH(I,J,K)

          END DO
        END DO
      END DO

      DO j=1,rows
      DO I=1,row_length
!       !----------------------------------------------------------
!       ! Use local NTML if significantly higher (to allow for
!       ! local averaging) than the non-local or if the non-local
!       ! is on the ground (=1)
!       !----------------------------------------------------------
        IF ( .NOT.CUMULUS(I,j) .AND.                                    &
     &            ( NTML_LOCAL(I,j)  >   NTML_NL(I,j)+1                 &
     &              .OR. NTML_NL(I,j)  ==  1 )            ) THEN
          NTML(I,j) = NTML_LOCAL(I,j)
          SML_DISC_INV(I,j) = 0   ! reset flag for subgrid inversion
        ELSE
          NTML(I,j) = NTML_NL(I,j)
        END IF
!       !----------------------------------------------------------
!       ! If local NTML is higher than NTDSC then ignore DSC layer
!       ! for diagnostics but keep mixing associated with it
!       !----------------------------------------------------------
        IF ( NTML_LOCAL(I,j)  >   NTDSC(I,j)+1 ) THEN
          DSC_DISC_INV(I,j) = 0
          NTDSC(I,j) = 0
          NBDSC(I,j) = 0
          ZHSC(I,j)  = 0.0
          DSC(I,j)   = .FALSE.
          COUPLED(I,j) = .FALSE.
        END IF
      END DO
      END DO

! Calculate max of two coeffs
! DEPENDS ON: kmkh
      CALL KMKH (                                                       &
     & row_length,rows,BL_LEVELS,off_x,off_y                            &
     &,Keep_Ri_FA,RHOKM,RHO_KM(1,1,2),RHOKH                             &
     &,RHOKMZ(1,1,2),RHOKHZ(1,1,2)                                      &
     &,NTML,CUMULUS,RHOKM_TOP(1,1,2),RHOKH_TOP(1,1,2),UNSTABLE          &
     &,NTDSC,DSC                                                        &
     &,SML_DISC_INV,DSC_DISC_INV                                        &
     &,nSCMDpkgs,L_SCMDiags                                             &
     &,LTIMER                                                           &
     & )

      IF (BL_diag%L_rhokm) THEN
        DO k = 1, bl_levels
          DO j = 1, rows
            DO i = 1, row_length
              BL_diag%rhokm(i,j,k)=RHOKM(i,j,k)
            END DO
          END DO
        END DO
      END IF

      IF (BL_diag%L_rhokh) THEN
        DO k = 1, bl_levels
          DO j = 1, rows
            DO i = 1, row_length
              BL_diag%rhokh(i,j,k)=RHOKH(i,j,k)
            END DO
          END DO
        END DO
      END IF
! Assignment of TKE diagnostic
      IF (BL_diag%L_tke) THEN
       DO k = 2, bl_levels
         DO j = 1, rows
            DO i = 1, row_length


             VKZ = VKMAN * ( Z_TQ(I,j,K-1) + Z0M_EFF_GB(I,j) )

             IF ( FB_SURF(I,j) > 0.0 ) THEN
!              ! convective bl length scale
               l_int=LAMBDA_MIN

               IF (Z_TQ(I,j,K-1) < ZH(i,j)) THEN
                 l_int= (ZH(i,j)-Z_TQ(I,j,K-1)) /                       &
     &           ( 1.0 +  (ZH(i,j)-Z_TQ(I,j,K-1))/VKZ )
               END IF      ! z < zh

             ELSE
! neutral or stable bl scale

               l_int= VKZ / (1.0 + VKZ/LAMBDA_MIN )
             END IF   ! FB_SURF(I,j) > 0.0

! TKE diagnostic
            BL_diag%tke(i,j,k)=                                         &
     &      (RHOKM(i,j,k)/(RHO_TQ(i,j,k-1)*l_int))**2



            END DO
         END DO
       END DO
! k=1 values for TKE diagnostic
          DO j = 1, rows
            DO i = 1, row_length


             VKZ = VKMAN *  Z0M_EFF_GB(I,j)

             IF ( FB_SURF(I,j) > 0.0 ) THEN
!              ! convective bl length scale
               l_int=  ZH(i,j) /                                        &
     &         ( 1.0 +  ZH(i,j)/VKZ )

             ELSE
! neutral or stable bl scale

               l_int= VKZ / (1.0 + VKZ/LAMBDA_MIN )
             END IF   ! FB_SURF(I,j) > 0.0

! TKE diagnostic
            BL_diag%tke(i,j,1)=                                         &
     &      (RHOKM(i,j,1)/(RHOSTAR(i,j)*l_int))**2



            END DO
          END DO
      END IF    ! BL_diag%L_tke

      If (L_subfilter_vert) then

        ALLOCATE (visc_BL_h(1:row_length, 1:rows, bl_levels))
        ALLOCATE (visc_BL_h_rho(1:row_length, 1:rows, bl_levels))

        Do k = 1, BL_LEVELS-1
          DO j = 1, rows
            DO i = 1, row_length
              BL_COEF_KM(i,j,k) = RHOKM(i,j,k+1)
              BL_COEF_KH(i,j,k) = RHOKH(i,j,k+1)
            End Do
          End Do
        End Do
!
! visc_bl_m is now lambda^2*S
!
        DO k = turb_startlev_vert, turb_endlev_vert
          DO j = 1, rows
            DO i = 1, row_length
              visc_bl_h(i,j,k) = visc_bl_m(i,j,k)*                      &
     &                           FH_3D(i,j,k)
              visc_bl_m(i,j,k) = visc_bl_m(i,j,k)*                      &
     &                           FM_3D(i,j,k)
              visc_bl_h(i,j,k) = min(visc_bl_h(i,j,k),max_diff)
              visc_bl_m(i,j,k) = min(visc_bl_m(i,j,k),max_diff)
              visc_bl_m(i,j,k) = visc_bl_m(i,j,k)*rho_tq(i,j,k-1)
!!  Don't multiply visc_bl_h by rho yet but after interpolation 
!!  to rho levels
            End Do
          End Do
        End Do
!
! visc_bl_m and visc_bl_h are now rho*lambda^2*S*FM and lambda^2*S*FH
!
        DO k = 2, bl_levels
          DO j = 1, rows
            DO i = 1, row_length
              WEIGHT1 = r_theta_levels(I,J,K) -                         &
     &                r_theta_levels(I,J, K-1)
              WEIGHT2 = r_theta_levels(I,J,K) -                         &
     &                r_rho_levels(I,J,K)
              WEIGHT3 = r_rho_levels(I,J,K) -                           &
     &                r_theta_levels(I,J,K-1)
              If ( K  ==  BL_LEVELS ) THEN
! assume visc_bl_h(BL_LEVELS+1) is zero
                visc_bl_h_rho(I,J,K) = ( WEIGHT2/WEIGHT1 )              &
     &                                * visc_bl_h(I,J,K)
              Else
                visc_bl_h_rho(I,J,K) =                                  &
     &              WEIGHT3/WEIGHT1 *                                   &
     &                      visc_bl_h(I,J,K+1)                          &
     &             +WEIGHT2/WEIGHT1 *                                   &
     &                      visc_bl_h(I,J,K)
              End If
            End Do
          End Do
        End Do
!
! Apply blending of coefficients from BL and 3D turb scheme
!
        If (L_subfilter_blend) then

! DEPENDS ON: blend_coeff
          Call Blend_coeff(visc_BL_m, visc_BL_h_rho, RHOKM, RHOKH       &
     &,                         row_length, rows, BL_LEVELS )

        End If

!
! Overwrite the diffusion coefficients from the local BL scheme
!(RHOKM and RHOKH) with those obtained from the Smagorinsky scheme.
!
        Do k = 2, BL_LEVELS
          If (k >= turb_startlev_vert .AND.                             &
     &        k <= turb_endlev_vert) Then
            Do j = 1, rows
              Do i = 1, row_length
                RHOKM(i,j,k) = visc_BL_m(i,j,k)
! Finally multiply visc_BL_h by density on uv levels
                visc_BL_h_rho(i,j,k)=visc_BL_h_rho(i,j,k)*rho_uv(i,j,k)
                RHOKH(i,j,k) = visc_BL_h_rho(i,j,k)
              End Do
            End Do
          Else
            Do j = 1, rows
              Do i = 1, row_length
                RHOKM(i,j,k) = 0.0
                RHOKH(i,j,k) = 0.0
              End Do
            End Do
          End If
        End Do

        DEALLOCATE (visc_BL_h)
        DEALLOCATE (visc_BL_h_rho)

      End If ! L_subfilter_vert
!
!-----------------------------------------------------------------------
! Diagnose boundary layer type.
!      Seven different types are considered:
!      1 - Stable b.l.
!      2 - Stratocumulus over a stable surface layer.
!      3 - Well mixed buoyancy-driven b.l. (possibly with stratocumulus)
!      4 - Decoupled stratocumulus (not over cumulus).
!      5 - Decoupled stratocumulus over cumulus.
!      6 - Cumulus capped b.l.
!      7 - Shear-dominated unstable b.l.
!-----------------------------------------------------------------------
!      First initialise the type variables and set the diagnostic ZHT.
!
      DO j = 1, rows
        DO I = 1, row_length
            BL_TYPE_1(I,j) = 0.0
            BL_TYPE_2(I,j) = 0.0
            BL_TYPE_3(I,j) = 0.0
            BL_TYPE_4(I,j) = 0.0
            BL_TYPE_5(I,j) = 0.0
            BL_TYPE_6(I,j) = 0.0
            BL_TYPE_7(I,j) = 0.0
          ZHT(I,j) = MAX( ZH(I,j) , ZHSC(I,j) )
          IF ( NTML(I,j)  >   NTML_NL(I,j) ) THEN
!           ! Higher local K allowed so reset ZH, ZHT diagnostics
            ZH(I,j) = MAX( ZH(I,j) , ZH_LOCAL(I,j) )
            IF ( .NOT.DSC(I,j) )                                        &
     &        ZHT(I,j) = MAX( ZHT(I,j) , ZH_LOCAL(I,j) )
          END IF
          LWP(I,j) = (QCL(I,j,1)+QCF(I,j,1) )*Z_UV(I,j,2)  ! ???
          DO K=2,BL_LEVELS
            LWP(I,j) = LWP(I,j)                                         &
     &                  + ( QCL(I,j,K)+QCF(I,j,K) )*DZL_charney(I,j,K)
          END DO
        END DO
      END DO
      DO j = 1, rows
      DO I = 1, row_length
        IF (.NOT.UNSTABLE(I,j) .AND. .NOT.DSC(I,j) .AND.                &
     &       .NOT.CUMULUS(I,j)) THEN
!         Stable b.l.
          BL_TYPE_1(I,j) = 1.0
        ELSEIF (.NOT.UNSTABLE(I,j) .AND. DSC(I,j) .AND.                 &
     &       .NOT.CUMULUS(I,j)) THEN
!         Stratocumulus over a stable surface layer
          BL_TYPE_2(I,j) = 1.0
        ELSEIF (UNSTABLE(I,j) .AND. .NOT.CUMULUS(I,j) .AND.             &
     &       .NOT.DSC(I,j) ) THEN
!         Well mixed b.l. (possibly with stratocumulus)
          IF ( NTML(I,j)  >   NTML_NL(I,j) ) THEN
!           ! shear-dominated - currently identified
!           ! by local NTML overriding non-local
            BL_TYPE_7(I,j) = 1.0
          ELSE
!           ! buoyancy-dominated
          BL_TYPE_3(I,j) = 1.0
          END IF
        ELSEIF (UNSTABLE(I,j) .AND. DSC(I,j) .AND.                      &
     &                                   .NOT.CUMULUS(I,j)) THEN
!         Decoupled stratocumulus (not over cumulus)
          BL_TYPE_4(I,j) = 1.0
        ELSEIF (DSC(I,j) .AND. CUMULUS(I,j)) THEN
!         Decoupled stratocumulus over cumulus
          BL_TYPE_5(I,j) = 1.0
        ELSEIF (.NOT.DSC(I,j) .AND. CUMULUS(I,j)) THEN
!         Cumulus capped b.l.
          BL_TYPE_6(I,j) = 1.0
        END IF
      END DO
      END DO
!
!-----------------------------------------------------------------------
!! 5.4 Interpolate RHOKM's to uv points ready for the
!!     calculation of the explcit fluxes TAU_X and TAU_Y at levels
!!     above the surface.
!-----------------------------------------------------------------------

#if !defined(SCMA)

      IF(FORMDRAG ==  Explicit_stress)THEN

       DO K=1,BL_LEVELS
        DO J = 1, ROWS
         DO I=1, ROW_LENGTH
           taux_hal(i,j,k) = tau_fd_x(i,j,k)
           tauy_hal(i,j,k) = tau_fd_y(i,j,k)
         END DO
        END DO
       END DO

! DEPENDS ON: swap_bounds
       CALL swap_bounds(taux_hal,row_length,rows,bl_levels,             &
     &      off_x, off_y, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
       CALL swap_bounds(tauy_hal,row_length,rows,bl_levels,             &
     &      off_x, off_y, fld_type_p, .false.)

! DEPENDS ON: fill_external_halos
       CALL fill_external_halos(taux_hal,row_length,rows,bl_levels,     &
     &                          off_x,off_y)
! DEPENDS ON: fill_external_halos
       CALL fill_external_halos(tauy_hal,row_length,rows,bl_levels,     &
     &                          off_x,off_y)

! DEPENDS ON: p_to_u
       CALL p_to_u(taux_hal,row_length,rows,bl_levels,                  &
     &             off_x, off_y,taux_fd_u)

! DEPENDS ON: p_to_v
       CALL p_to_v(tauy_hal,row_length, rows, n_rows,                   &
     &             bl_levels, off_x, off_y,tauy_fd_v)

      END IF

! DEPENDS ON: swap_bounds
      Call Swap_Bounds(                                                 &
     &                 RHOKM(1-off_x,1-off_y,2), row_length, rows,      &
     &                 bl_levels-1, off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(RHOKM(1-off_x,1-off_y,2),ROW_LENGTH,ROWS,&
     &                         BL_LEVELS-1,off_x,off_y)

! DEPENDS ON: p_to_u
      CALL P_TO_U(RHOKM(1-off_x,1-off_y,2), ROW_LENGTH, ROWS,           &
     &            BL_LEVELS-1, off_x, off_y, RHOKM_U(1,1,2))

! DEPENDS ON: p_to_v
      CALL P_TO_V(RHOKM(1-off_x,1-off_y,2),ROW_LENGTH, ROWS, n_rows,    &
     &            BL_LEVELS-1, off_x, off_y, RHOKM_V(1,1,2))

      DO K=2,BL_LEVELS
      DO J= 1, rows
      DO I= 1, row_length
        F_NGSTRESS_HALO(I,J,K) = F_NGSTRESS(I,J,K)
      END DO
      END DO
      END DO
! DEPENDS ON: swap_bounds
      Call Swap_Bounds(F_NGSTRESS_HALO, row_length, rows,               &
     &                 bl_levels-1, off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(F_NGSTRESS_HALO,ROW_LENGTH,ROWS,         &
     &                         BL_LEVELS-1,off_x,off_y)

! DEPENDS ON: p_to_u
      CALL P_TO_U(F_NGSTRESS_HALO,ROW_LENGTH,ROWS,BL_LEVELS-1,          &
     &            off_x, off_y, F_NGSTRESS_U)

! DEPENDS ON: p_to_v
      CALL P_TO_V(F_NGSTRESS_HALO,ROW_LENGTH, ROWS, n_rows,             &
     &            BL_LEVELS-1, off_x, off_y, F_NGSTRESS_V)

#else
      DO K=2,BL_LEVELS
      DO J= 1, rows
      DO I= 1, row_length
        RHOKM_V(I,J,K) = RHOKM(I,J,K)
        RHOKM_U(I,J,K) = RHOKM(I,J,K)
        IF(FORMDRAG ==  Explicit_stress)THEN
         taux_fd_u(i,j,k) = tau_fd_x(i,j,k)
         tauy_fd_v(i,j,k) = tau_fd_y(i,j,k)
        END IF

        F_NGSTRESS_V(I,J,K) = F_NGSTRESS(I,J,K)
        F_NGSTRESS_U(I,J,K) = F_NGSTRESS(I,J,K)
      END DO
      END DO
      END DO
      IF(FORMDRAG ==  Explicit_stress)THEN

       DO J= 1, rows
        DO I= 1, row_length
          taux_fd_u(i,j,1) = tau_fd_x(i,j,1)
          tauy_fd_v(i,j,1) = tau_fd_y(i,j,1)
        END DO
       END DO

      END IF

#endif

!-----------------------------------------------------------------------
!! 5.5 Calculation of explicit fluxes of T,Q
!-----------------------------------------------------------------------

! DEPENDS ON: ex_flux_tq
      CALL EX_FLUX_TQ (                                                 &
     &  row_length,rows,BL_LEVELS,FLUX_GRAD                             &
     &, TL,QW,RDZ_charney_grid,FTL,FQW                                  &
     &, RHOKH, RHOKHZ(1,1,2)                                            &
     &, GRAD_T_ADJ, GRAD_Q_ADJ, RHOF2, RHOFSC                           &
     &, FT_NT, FQ_NT, FT_NT_DSCB, FQ_NT_DSCB                            &
     &, TOTHF_ZH, TOTHF_ZHSC, TOTQF_ZH, TOTQF_ZHSC                      &
     &, NTML_NL, NTDSC, NBDSC                                           &
     &, nSCMDpkgs,L_SCMDiags                                            &
     &, LTIMER                                                          &
     &  )

!-----------------------------------------------------------------------
!! 5.6 Calculation of explicit fluxes of U and V.
!-----------------------------------------------------------------------

! DEPENDS ON: ex_flux_uv
      CALL EX_FLUX_UV (                                                 &
                        ! For U
     &  ROW_LENGTH, rows, off_x, off_y, BL_LEVELS                       &
     &, U,ZH,RDZ_U,RHOKM_U,NG_STRESS,F_NGSTRESS_U,TAUX                  &
     &, TAUX_FD_U,FORMDRAG                                              &
     &, TAU_GRAD_U,TAU_NON_GRAD_U                                       &
     &, LTIMER                                                          &
     &  )

! DEPENDS ON: ex_flux_uv
      CALL EX_FLUX_UV (                                                 &
                        ! For V
     &  ROW_LENGTH, n_rows, off_x, off_y, BL_LEVELS                     &
     &, V,ZH,RDZ_V,RHOKM_V,NG_STRESS,F_NGSTRESS_V,TAUY                  &
     &, TAUY_FD_V,FORMDRAG                                              &
     &, TAU_GRAD_V,TAU_NON_GRAD_V                                       &
     &, LTIMER                                                          &
     &  )

#if defined(SCMA)
!
!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
      IF (L_SCMDiags(SCMDiag_bl)) THEN
!
! DEPENDS ON: scmoutput
        Call SCMoutput(TAU_GRAD_U,                                      &
             'taux_grad','Gradient part of u stress','kg/m/s2',         &
             t_avg, d_bl, default_streams, '',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(TAU_NON_GRAD_U,                                  &
             'taux_nongrad','Non-gradient part of u stress',            &
             'kg/m/s2',                                                 &
             t_avg, d_bl, default_streams, '',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(TAU_GRAD_V,                                      &
             'tauy_grad','Gradient part of v stress','kg/m/s2',         &
             t_avg, d_bl, default_streams, '',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(TAU_NON_GRAD_V,                                  &
             'tauy_nongrad','Non-gradient part of v stress',            &
             'kg/m/s2',                                                 &
             t_avg, d_bl, default_streams, '',RoutineName)
!
      END IF ! L_SCMDiags(SCMDiag_bl)
#endif
!
!-----------------------------------------------------------------------
!! 5.6.1 Calculate explicit surface fluxes of U and V on
!!       P-grid for convection scheme
!-----------------------------------------------------------------------

      DO j=1,rows
       DO i=1,row_length
         uw0(I,j) = -RHOKM(I,j,1) *                                     &
     &                     ( U_P(I,j,1) - U_0_P(I,j) )
         vw0(I,j) = -RHOKM(I,j,1) *                                     &
     &                     ( V_P(I,j,1) - V_0_P(I,j) )
       END DO
      END DO

!-----------------------------------------------------------------------
!! 5.7 Set NTML to max number of turbulently mixed layers
!      Calculate quantities to pass to convection scheme.
!-----------------------------------------------------------------------
          DO J= 1, rows
            DO I= 1, row_length
             WSTAR(I,j) = 0.0
             WTHVS(I,j) = 0.0
             CU_OVER_OROG(I,J) = 0.0
             IF ( CUMULUS(I,J) ) THEN
               IF ( FB_SURF(I,j)  >   0.0 ) THEN
                 WSTAR(I,j) = ( ZH(I,j)*FB_SURF(I,j) )**(1.0/3.0)
                 WTHVS(I,j) = FB_SURF(I,j) / ( G * BT(I,j,1) )
               END IF
               WSTAR(I,j) = MAX( 0.1, WSTAR(I,j) )
               IF (.NOT. L_CONV4A) THEN
               NTML(I,j) = MAX( 2, NTML_NL(I,J) - 1 )
               END IF
             ELSE
               NTML(I,j) = MAX( NTML_NL(I,J) , NTDSC(I,J) )
             END IF
! Limit explicitly calculated surface stresses
! to a physically plausible level.
             IF ( UW0(I,j)  >=  5.0 ) THEN
               UW0(I,j) =  5.0
             ELSEIF ( UW0(I,j)  <=  -5.0 ) THEN
               UW0(I,j) = -5.0
             END IF
             IF ( VW0(I,j)  >=  5.0 ) THEN
               VW0(I,j) =  5.0
             ELSEIF ( VW0(I,j)  <=  -5.0 ) THEN
               VW0(I,j) = -5.0
             END IF
             IF (BL_diag%L_wstar .AND. (FB_SURF(I,j) >0.0))  THEN
               BL_diag%wstar(i,j)= (ZH(i,j)*FB_SURF(i,j))**(1.0/3.0)
             END IF
           END DO
          END DO

      IF (L_conv4a) THEN

! Check for CUMULUS having been diagnosed over steep orography.
! Reset to false but keep NTML at NLCL (though decrease by 2 so that
! coupling between BL and convection scheme can be maintained).
! Reset type diagnostics.

      DO L=1,LAND_PTS
        J=(land_index(L)-1)/row_length + 1
        I=land_index(L) - (J-1)*row_length
        IF (CUMULUS(I,J) .AND. HO2R2_OROG(L)  >   900.0) THEN
          CUMULUS(I,J) = .FALSE.
          L_SHALLOW(I,J) = .FALSE.
          BL_TYPE_5(I,J) = 0.0
          BL_TYPE_6(I,J) = 0.0
          CU_OVER_OROG(I,J) = 1.0
          IF (NTML(I,J)  >=  3) NTML(I,J) = NTML(I,J) - 2
        END IF
      END DO

! Check that CUMULUS and L_SHALLOW are still consistent

          DO J= 1, rows
            DO I= 1, row_length
              IF ( .NOT. CUMULUS(I,j) ) L_SHALLOW(I,j) = .FALSE.
            END DO
          END DO

      END IF    ! (L_conv4a)

!-----------------------------------------------------------------------
!     Set shallow convection diagnostic: 1.0 if L_SHALLOW (and CUMULUS)
!                                        0.0 if .NOT. CUMULUS
!-----------------------------------------------------------------------
          DO J= 1, rows
            DO I= 1, row_length
              IF ( CUMULUS(I,j) .AND. L_SHALLOW(I,j) ) THEN
                SHALLOWC(I,j) = 1.0
              ELSE
                SHALLOWC(I,j) = 0.0
              END IF
            END DO
          END DO


      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('BDY_EXPL2 ',4)
      END IF

      RETURN
      END SUBROUTINE BDY_EXPL2
#endif
