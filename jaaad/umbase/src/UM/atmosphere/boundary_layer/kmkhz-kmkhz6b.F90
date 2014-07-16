#if defined(A03_8B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE KMKHZ  --------------------------------------------------
!!!
!!!  Purpose: To calculate the non-local turbulent mixing
!!!           coefficients KM and KH
!!!
!!!  Programming standard:
!!!
!!!  Documentation: UMDP No.24
!!!
!!!---------------------------------------------------------------------

! Arguments :-

      SUBROUTINE KMKHZ (                                                &
     & row_length,rows,BL_LEVELS,halo_i,halo_j,n_proc,                  &
     & DECFIX, STOPWE_SBL, NG_STRESS, FLUX_GRAD, lq_mix_bl,             &
     & TIMESTEP,P,P_HALF,RHO_TQ,RHO_UV,RHO_DRY_TQ,                      &
     & T,Q,QCL,QCF,CF,QW,TL,                                            &
     & DZL,RDZ,Z_FULL,Z_HALF,Z_UV,Z_TQ,                                 &
     & FTL,FQW,RHOKM,RHOKH,RHOKM_TOP,RHOKH_TOP,                         &
     & W,ETADOT,ETA_THETA_LEVELS,RAD_HR,MICRO_TENDS,                    &
     & BT,BQ,BTM,BQM,DQSDT,BTM_CLD,BQM_CLD,A_QS,A_QSM,A_DQSDTM,CFM,     &
     & Z0M,V_S,FB_SURF,T1_SD,Q1_SD,RHOSTAR_GB,                          &
     & UNSTABLE,DSC,COUPLED,SML_DISC_INV,DSC_DISC_INV,                  &
     & CUMULUS,L_SHALLOW, NTML,NTPAR,NLCL,NTDSC,NBDSC,                  &
     & ZH,ZH_PREV,ZHPAR,Z_LCL,ZHSC,ZMAX_FOR_DSC,                        &
     & F_NGSTRESS, GRAD_T_ADJ, GRAD_Q_ADJ,                              &
     & RHOF2, RHOFSC, FT_NT, FQ_NT, FT_NT_DSCB, FQ_NT_DSCB,             &
     & TOTHF_ZH, TOTHF_ZHSC, TOTQF_ZH, TOTQF_ZHSC,                      &
     & KENT, WE_LIM, T_FRAC, ZRZI,                                      &
     & KENT_DSC, WE_LIM_DSC, T_FRAC_DSC, ZRZI_DSC,                      &
     & BL_diag,nSCMDpkgs,L_SCMDiags,                                    &
     & LTIMER                                                           &
     & )

      Use cv_run_mod, Only:                                             &
          l_conv4a

      Use bl_diags_mod, Only:                                           &
          strnewbldiag

      IMPLICIT NONE

      INTEGER                                                           &
     & row_length,rows,halo_i,halo_j,n_proc                             &
     &,BL_LEVELS                                                        &
                              ! IN No. of atmospheric levels for
!                             !    which boundary layer fluxes are
!                             !    calculated.
     &,NG_STRESS                                                        &
                              ! IN switch for non-gradient stress
     &,DECFIX                                                           &
                              ! IN correction to decoupling diagnosis
     &,STOPWE_SBL                                                       &
                              ! IN stop spurious entrainment in SBLs
     &,FLUX_GRAD              ! IN switch for flux-gradient formulation

      LOGICAL                                                           &
     & LTIMER                                                            
                              ! IN Flag for TIMER diagnostics
      LOGICAL                                                           &
     & lq_mix_bl
      INTEGER                                                           &
     & NTPAR(row_length,rows)                                           &
                              ! IN Top level of parcel ascent.
!                             !    Used in convection scheme.
!                             !    NOTE: CAN BE > BL_LEVELS-1
     &,NLCL(row_length,rows)  ! IN No. of model layers below the
!                             !    lifting condensation level.

      LOGICAL                                                           &
     & L_SHALLOW(row_length,rows)
!                             ! IN Flag to indicate shallow
!                             !    convection (only for A05_4A)
      REAL                                                              &
     & BQ(row_length,rows,BL_LEVELS)                                    &
!                             ! IN A buoyancy parameter for clear air
!                             !    on p,T,q-levels (full levels).
     &,BT(row_length,rows,BL_LEVELS)                                    &
!                             ! IN A buoyancy parameter for clear air
!                             !    on p,T,q-levels (full levels).
     &,BQM(row_length,rows,BL_LEVELS)                                   &
!                             ! IN A buoyancy parameter for clear air
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,BTM(row_length,rows,BL_LEVELS)                                   &
!                             ! IN A buoyancy parameter for clear air
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,BQM_CLD(row_length,rows,BL_LEVELS)                               &
!                             ! IN A buoyancy parameter for cloudy air
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,BTM_CLD(row_length,rows,BL_LEVELS)                               &
!                             ! IN A buoyancy parameter for cloudy air
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,A_QS(row_length,rows,BL_LEVELS)                                  &
                              ! IN Saturated lapse rate factor
!                             !    on p,T,q-levels (full levels).
     &,A_QSM(row_length,rows,BL_LEVELS)                                 &
!                             ! IN Saturated lapse rate factor
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,A_DQSDTM(row_length,rows,BL_LEVELS)                              &
!                             ! IN Saturated lapse rate factor
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,P(row_length,rows,BL_LEVELS)                                     &
                              ! IN P(*,K) is pressure at full level k.
     &,P_HALF(row_length,rows,BL_LEVELS)                                &
!                             ! IN P_HALF(*,K) is pressure at half
!                             !    level k-1/2.
     &,QW(row_length,rows,BL_LEVELS)                                    &
                              ! IN Total water content (kg per kg air).
     &,TL(row_length,rows,BL_LEVELS)                                    &
                              ! IN Liquid/frozen water temperature (K).
     &,CFM(row_length,rows,BL_LEVELS)                                   &
!                             ! IN Estimate of cloud fraction
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,T(row_length,rows,BL_LEVELS)                                     &
!                             ! IN Temperature (K).
     &,QCF(row_length,rows,BL_LEVELS)                                   &
                                     ! IN Cloud ice (kg per kg air)
     &,QCL(row_length,rows,BL_LEVELS)                                   &
                                     ! IN Cloud liquid water
     &,Q(row_length,rows,BL_LEVELS)                                     &
                                     ! IN specific humidity

!                             ! IN specific humidity
     &,CF(row_length,rows,BL_LEVELS)                                    &
                              ! IN Cloud fractions for boundary levs.
     &,DZL(row_length,rows,BL_LEVELS)                                   &
!                             ! IN Layer depths (m).  DZL(,K) is the
!                             !    distance from layer boundary K-1/2
!                             !    to layer boundary K+1/2.  For K=1
!                             !    the lower boundary is the surface.
     &,ZH_PREV(row_length, rows)                                        &
!                             ! IN boundary layer height (m) from
!                             !    previous timestep
     &,RDZ(row_length,rows,BL_LEVELS)                                   &
!                             ! IN Reciprocal of distance between
!                             !    full levels (m-1).  1/RDZ(,K) is
!                             !    the vertical distance from level
!                             !    K-1 to level K, except that for
!                             !    K=1 it is the height of the
!                             !    lowest atmospheric full level.
     &,Z0M(row_length,rows)                                             &
                              ! IN Roughness length for momentum (m).
     &,Z_FULL(row_length,rows,BL_LEVELS)                                &
!                             ! IN Z_FULL(*,K) is the height of the
!                             !    k-th full level above the surface.
     &,Z_HALF(row_length,rows,BL_LEVELS)                                &
!                             ! IN Z_HALF(*,K) is the height of level
!                             !       k-1/2 above the surface (m).
     &,RAD_HR(row_length,rows,BL_LEVELS,2)                              &
!                             ! IN (LW,SW) radiative heating rates
!                             !    (K/s)
     &,MICRO_TENDS(row_length, rows, bl_levels, 2)                      &
!                             ! IN Tendencies from microphysics
!                             !    (TL, K/s; QW, kg/kg/s)
     &,V_S(row_length,rows)                                             &
                              ! IN Surface friction velocity (m/s)
     &,FB_SURF(row_length,rows)                                         &
                               !IN Surface buoyancy flux over density
!                             !       (m^2/s^3).
     &,FQW(row_length,rows,BL_LEVELS)                                   &
!                             ! IN "Explicit" flux of QW (i.e.
!                             !        evaporation) from layer below
!                             !        on P-grid (kg per sq m per s).
     &,FTL(row_length,rows,BL_LEVELS)                                   &
!                             ! IN "Explicit" flux of TL = H/CP
!                             !       (sensible heat/CP) from layer
!                             !       below, on P-grid.
     &,T1_SD(row_length,rows)                                           &
                              ! IN Standard Deviation of level 1
!                             !    temperature (K).
     &,Q1_SD(row_length,rows)                                           &
                              ! IN Standard Deviation of level 1
!                             !    specific humidity (kg/kg).
     &,DQSDT(row_length,rows,BL_LEVELS)                                 &
!                             ! IN Partial derivative of QSAT w.r.t.
!                             !    temperature.
     &,Z_UV(row_length,rows,BL_LEVELS)                                  &
!                             ! IN For a vertically staggered grid
!                             !    with a u,v-level first above the
!                             !    surface, Z_UV(*,K) is the height
!                             !    of the k-th u,v-level (half level
!                             !    k-1/2) above the surface;
!                             !    for an unstaggered grid the
!                             !    heights of the half-levels
!                             !    0.5 to BL_LEVELS-0.5 should be
!                             !    input to elements 1 to BL_LEVELS.
!                             !    (1st value not used in either
!                             !     case.)
     &,Z_TQ(row_length,rows,BL_LEVELS)                                  &
!                             ! IN For a vertically staggered grid
!                             !    with a u,v-level first above the
!                             !    surface, Z_TQ(*,K) is the height
!                             !    of the k-th T,q-level (full level
!                             !    k) above the surface;
!                             !    for an unstaggered grid the
!                             !    heights of the half levels
!                             !    1.5 to BL_LEVELS+0.5 should be
!                             !    input to elements 1 to BL_LEVELS.
!                             !    (value for BL_LEVELS not used
!                             !    in either case.)
!
     &,RHO_UV(row_length,rows,BL_LEVELS)                                &
!                             ! IN density on UV (ie. rho) levels,
!                             !    used in RHOKH so dry density if
!                             !    Lq_mix_bl is true
     &,RHO_TQ(row_length,rows,BL_LEVELS)                                &
!                             ! IN density on TQ (ie. theta) levels,
!                             !    used in RHOKM so wet density
     &,RHO_DRY_TQ(row_length,rows,BL_LEVELS)                            &
!                             ! IN density on TQ (ie. theta) levels,
!                             !    used in non-turb flux integration
!                             !    so dry density if Lq_mix_bl is true
     &,RHOSTAR_GB(row_length,rows)                                      &
!                             ! IN Surface air density in kg per
!                             !    cubic metre.
     &,TIMESTEP                                                         &
!                             ! IN Model timestep (s)
     &,ZMAX_FOR_DSC                                                     &
!                             ! IN Max height to look for DSC cloud-base
     &,W(row_length, rows, 0:BL_LEVELS)                                 &
!                             ! IN Vertical velocity (m/s)
     &,ETADOT(row_length, rows, 0:BL_LEVELS)                            &
!                             ! IN d(ETA)/dt
     &,eta_theta_levels(0:bl_levels)                                    &
!                             ! IN eta coordinate
     &,Z_LCL(row_length,rows)                                           &
                              ! IN Height of lifting condensation
!                             !    level.
     &,ZHPAR(row_length,rows) ! IN Height of top of NTPAR
!                             !    NOTE: CAN BE ABOVE BL_LEVELS-1
      INTEGER                                                           &
     & NTML(row_length,rows)  ! INOUT Number of model levels in the
!                             !       turbulently mixed layer.

      LOGICAL                                                           &
     & CUMULUS(row_length,rows)
!                             ! INOUT Flag for Cu in the bl

      REAL                                                              &
     & ZH(row_length,rows)    ! INOUT Boundary layer height (m).

      INTEGER                                                           &
     & NTDSC(row_length,rows)                                           &
                              ! OUT Top level for turb mixing in
!                             !     cloud layer
     &,NBDSC(row_length,rows)                                           &
                              ! OUT Bottom level of any decoupled
!                             !     turbulently mixed Sc layer.
     &,SML_DISC_INV(row_length,rows)                                    &
!                             ! OUT Flags for whether discontinuous
     &,DSC_DISC_INV(row_length,rows)                                    &
!                             ! OUT inversions are diagnosed
     &,KENT(row_length,rows)                                            &
                              ! OUT grid-level of SML inversion
     &,KENT_DSC(row_length,rows)
!                             ! OUT grid-level of DSC inversion

      LOGICAL                                                           &
     & UNSTABLE(row_length,rows)                                        &
!                             ! OUT Flag to indicate an unstable
!                             !     surface layer.
     &,DSC(row_length,rows)                                             &
                              ! OUT Flag set if decoupled
!                             !    stratocumulus layer found
     &,COUPLED(row_length,rows)! OUT Flag to indicate Sc layer weakly
!                             !     coupled to surface (ie weakly
!                             !     decoupled)
      REAL                                                              &
     & GRAD_T_ADJ(row_length,rows)                                      &
!                             ! OUT Temperature gradient adjustment
!                             !     for non-local mixing in unstable
!                             !     turbulent boundary layer.
     &,GRAD_Q_ADJ(row_length,rows)                                      &
!                             ! OUT Humidity gradient adjustment
!                             !     for non-local mixing in unstable
!                             !     turbulent boundary layer.
     &,RHOKM(row_length,rows,2:BL_LEVELS)                               &
!                             ! OUT Non-local turbulent mixing
!                             !     coefficient for momentum.
     &,RHOKH(row_length,rows,2:BL_LEVELS)                               &
!                             ! OUT Non-local turbulent mixing
!                             !     coefficient for scalars.
     &,RHOKM_TOP(row_length,rows,2:BL_LEVELS)                           &
!                             ! OUT Top-down turbulent mixing
!                             !     coefficient for momentum.
     &,RHOKH_TOP(row_length,rows,2:BL_LEVELS)                           &
!                             ! OUT Top-down turbulent mixing
!                             !     coefficient for scalars.
     &,F_NGSTRESS(row_length,rows,2:BL_LEVELS)                          &
                              ! OUT dimensionless function for
!                             !     non-gradient stresses
     &,ZHSC(row_length,rows)                                            &
                              ! OUT Cloud layer height (m).
     &,WE_LIM(row_length,rows,3)                                        &
!                             ! OUT rho*entrainment rate implied by
!                             !     placing of subsidence
     &,ZRZI(row_length,rows,3)                                          &
                              ! OUT (z-z_base)/(z_i-z_base)
     &,T_FRAC(row_length,rows,3)                                        &
!                             ! OUT a fraction of the timestep
     &,WE_LIM_DSC(row_length,rows,3)                                    &
!                             ! OUT rho*entrainment rate implied by
!                             !     placing of subsidence
     &,ZRZI_DSC(row_length,rows,3)                                      &
!                             ! OUT (z-z_base)/(z_i-z_base)
     &,T_FRAC_DSC(row_length,rows,3)
!                             ! OUT a fraction of the timestep

!     Dummy variables (used in KMKHZ8C)
      REAL, DIMENSION(row_length,rows, BL_LEVELS+1) ::                  &
     &  FT_NT                                                           &
                              ! OUT Non-turbulent heat and moisture
     &, FQ_NT                 !       fluxes (rho*Km/s, rho*m/s)
      REAL, DIMENSION(row_length,rows, 2:BL_LEVELS) ::                  &
     &  RHOF2                                                           &
                              ! OUT f2 and fsc term shape profiles
     &, RHOFSC                !       multiplied by rho

      REAL, DIMENSION(row_length,rows) ::                               &
     &  TOTHF_ZH                                                        &
                              ! OUT Total heat fluxes at inversions
     &, TOTHF_ZHSC                                                      &
                              !     (rho*Km/s)
     &, TOTQF_ZH                                                        &
                              ! OUT Total moisture fluxes at
     &, TOTQF_ZHSC                                                      &
                              !      inversions (rho*m/s)
     &, FT_NT_DSCB                                                      &
                              ! OUT Non-turbulent heat and moisture
     &, FQ_NT_DSCB            !      flux at the base of the DSC layer

!!----------------------------------------------------------------------
!     Declaration of new BL diagnostics.
      Type (Strnewbldiag) :: BL_diag

! Additional variables for SCM diagnostics which are dummy in full UM
      INTEGER                                                           &
     &  nSCMDpkgs             ! No of SCM diagnostics packages

      LOGICAL                                                           &
     &  L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages
!!----------------------------------------------------------------------
!    External references :-
      EXTERNAL TIMER, QSAT_mix, EXCF_NL
!!----------------------------------------------------------------------
!    Local and other symbolic constants :-
#include "c_lheat.h"
#include "c_r_cp.h"
#include "c_g.h"
#include "c_0_dg_c.h"
#include "c_epslon.h"
#include "c_mdi.h"
#if defined(SCMA)
#include "s_scmop.h"
#endif
      Character*(*), Parameter ::  RoutineName = 'kmkhz'

      REAL ETAR,GRCP,LCRCP,LFRCP,LS,LSRCP
      PARAMETER (                                                       &
     & ETAR=1.0/(1.0-EPSILON)                                           &
                                ! Used in buoyancy parameter BETAC.
     &,GRCP=G/CP                                                        &
                                ! Adiabatic lapse rate.
     &,LCRCP=LC/CP                                                      &
                                ! Latent heat of condensation / CP.
     &,LFRCP=LF/CP                                                      &
                                ! Latent heat of fusion / CP.
     &,LS=LC+LF                                                         &
                                ! Latent heat of sublimation.
     &,LSRCP=LS/CP                                                      &
                                ! Latent heat of sublimation / CP.
     &)

      REAL A_PLUME,B_PLUME,A_GRAD_ADJ,MAX_T_GRAD,MAX_SVL_GRAD,SC_CFTOL, &
     & CT_RESID,SVL_COUP,DEC_SVL_GRAD,FGF
      PARAMETER (                                                       &
     & A_PLUME=0.2                                                      &
     &,B_PLUME=3.26                                                     &
     &,A_GRAD_ADJ=3.26                                                  &
     &,MAX_T_GRAD=1.0E-3                                                &
     &,MAX_SVL_GRAD=1.0E-3                                              &
                            ! maximum SVL gradient in a mixed layer
     &,DEC_SVL_GRAD=1.0E-3                                              &
                            ! SVL gradient required for weak decoupling
     &,SC_CFTOL=0.1                                                     &
                            ! CF required for a Sc layer to be diagnosed
     &,CT_RESID=200.                                                    &
                        ! Approx parcel cloud-top residence time (in s)
     &,SVL_COUP=0.5                                                     &
                        ! Parameter controlling positioning of
!                       ! surface-driven entrainment
     &,FGF=0.0)         ! Adiabatic gradient factor for ice
!
!  Define local storage.
!
!  (a) Workspace.
!
      LOGICAL                                                           &
     & CLOUD_BASE(row_length,rows)                                      &
                                   ! Flag set when cloud base
!                                  ! is reached.
     &,DSC_SAVE(row_length,rows)   ! Copy of DSC needed to indicate
!                                  ! decoupling diagnosed in EXCF_NL
      REAL                                                              &
     & QS(row_length,rows,BL_LEVELS)                                    &
                                    ! Saturated sp humidity at pressure
!                                  ! and temperature of sucessive
!                                  ! levels.
     &,QCL_IC_BOT(row_length,rows,BL_LEVELS)                            &
!                                  ! In-cloud liquid water content
!                                  ! at the bottom of the model layer
     &,QCF_IC_BOT(row_length,rows,BL_LEVELS)                            &
!                                  ! In-cloud frozen water content
!                                  ! at the bottom of the model layer
     &,QCL_IC_TOP(row_length,rows,BL_LEVELS)                            &
!                                  ! In-cloud liquid water content
!                                  ! at the top of the model layer.
     &,QCF_IC_TOP(row_length,rows,BL_LEVELS)                            &
!                                  ! In-cloud frozen water content
!                                  ! at the top of the model layer.
     &,CFL(row_length,rows,BL_LEVELS)                                   &
                                      ! Liquid cloud fraction.
     &,CFF(row_length,rows,BL_LEVELS)                                   &
                                      ! Frozen cloud fraction.
     &,DQCLDZ(row_length,rows,BL_LEVELS)                                &
!                                  ! Vertical gradient of
!                                  ! in-cloud liquid cloud water
!                                  ! in a well-mixed layer.
     &,DQCFDZ(row_length,rows,BL_LEVELS)                                &
!                                  ! Vertical gradient of in-cloud
!                                  ! frozen cloud water in a
!                                  ! well-mixed layer.
     &,CHI_S(row_length,rows,2:BL_LEVELS)                               &
!                                  ! Mixing fraction of just saturate
!                                  ! mixture.
     &,BFLUX_SURF(row_length,rows)                                      &
                                   ! Buoyancy flux at the surface.
     &,BFLUX_SURF_SAT(row_length,rows)                                  &
!                                  ! Saturated-air surface buoyancy flux
     &,DB_TOP(row_length,rows)                                          &
                                   ! Buoyancy jump at the top of the
!                                  ! boundary layer.
     &,DB_CLD(row_length,rows,2:BL_LEVELS)                              &
!                                  ! In-cloud buoyancy jump across
!                                  ! layer interface.
     &,DB(row_length,rows,2:BL_LEVELS)                                  &
                                       ! Buoyancy jump across layer
!                                      !  interface.
     &,DB_KSURF_DRY(row_length,rows,2:BL_LEVELS)                        &
                                                 ! Dry buoyancy jump in
!                                  ! flux integral calculation (m/s2)
     &,DB_KTOP_DRY(row_length,rows,2:BL_LEVELS)                         &
                                                ! Sat. buoyancy jump in
!                                  ! flux integral calculation (m/s2)
     &,DB_KSURF_CLD(row_length,rows,2:BL_LEVELS)                        &
                                                 ! Dry buoyancy jump in
!                                  ! flux integral calculation (m/s2)
     &,DB_KTOP_CLD(row_length,rows,2:BL_LEVELS)                         &
                                                ! Sat. buoyancy jump in
!                                  ! flux integral calculation (m/s2)
     &,DF_OVER_CP(row_length,rows,BL_LEVELS)                            &
!                                  ! Radiative flux change over layer
!                                  ! divided by c_P.
     &,DFLW_OVER_CP(row_length,rows,BL_LEVELS)                          &
!                                  ! LW radiative flux change over layer
!                                  ! divided by c_P.
     &,DFSW_OVER_CP(row_length,rows,BL_LEVELS)                          &
!                                  ! SW radiative flux change over layer
!                                  ! divided by c_P.
     &,TLS_INC(row_length,rows,BL_LEVELS)                               &
!                                  ! Static energy increment due to
!                                  !    large-scale vertical
!                                  !    advection (K s^-1)
     &,QLS_INC(row_length,rows,BL_LEVELS)                               &
!                                  ! Specific humidity increment
!                                  !    due to large-scale
!                                  !    vertical advection (s^-1)
     &,DF_TOP_OVER_CP(row_length,rows)                                  &
                                      ! Radiative flux change at cloud
!                                     ! top divided by c_P.
     &,SVL_PLUME(row_length,rows)                                       &
                                      ! Liquid/frozen water virtual
!                                     ! static energy over CP for a
!                                     ! plume rising without dilution
!                                     ! from level 1.
     &,SL_PLUME(row_length,rows)                                        &
                                      ! Liquid/frozen water
!                                     ! static energy over CP for a
!                                     ! plume rising without dilution
!                                     ! from level 1.
     &,QW_PLUME(row_length,rows)                                        &
                                      ! QW for a plume rising without
!                                     ! dilution from level 1.
     &,ENV_SVL_KM1(row_length,rows)                                     &
                                      ! Density potential temperature
!                                     ! for last layer considered
     &,Z_CLD(row_length,rows)                                           &
                                      ! Cloud fraction weighted
!                                     ! thickness of b.l. cloud.
     &,BT_TOP(row_length,rows)                                          &
                                      ! Buoyancy parameter at the top of
!                                     ! the b.l.
     &,BTT_TOP(row_length,rows)                                         &
                                      ! In-cloud buoyancy parameter at
!                                     ! the top of the b.l.
     &,BTC_TOP(row_length,rows)                                         &
                                      ! Cloud fraction weighted buoyancy
!                                     ! parameter at the top of the b.l.
     &,DB_TOP_CLD(row_length,rows)                                      &
                                      ! In-cloud buoyancy jump at the
!                                     ! top of the b.l.
     &,CLD_FACTOR(row_length,rows)                                      &
                                      ! Fraction of grid box potentially
!                                     ! giving evaporative entrainment.
     &,CHI_S_TOP(row_length,rows)                                       &
                                      ! Mixing fraction of just saturate
!                                     ! mixture at top of the b.l.
     &,ZETA_S(row_length,rows)                                          &
                                      ! Non-cloudy fraction of mixing
!                                     ! layer for surface forced
!                                     ! entrainment term.
     &,ZETA_R(row_length,rows)                                          &
                                      ! Non-cloudy fraction of mixing
!                                     ! layer for cloud top radiative
!                                     ! cooling entrainment term.
     &,ZC(row_length,rows)                                              &
                                      ! Cloud depth (not cloud fraction
!                                     ! weighted).
     &,DSCDEPTH(row_length,rows)                                        &
                                      ! Depth of cloud-layer (m)
     &,DB_DSCT(row_length,rows)                                         &
                                      ! Buoyancy jump at the top of the
!                                     ! surface mixed-layer.
     &,SVL(row_length,rows,BL_LEVELS)                                   &
                                      ! Liquid/frozen water virtual
!                                     ! temperature over CP.
     &,DF_DSCT_OVER_CP(row_length,rows)                                 &
                                       ! Radiative flux change at DSC
!                                     !  top divided by c_P.
     &,BT_DSCT(row_length,rows)                                         &
                                      ! Buoyancy parameter at the top of
!                                     ! the DSC
     &,BTT_DSCT(row_length,rows)                                        &
                                      ! In-cloud buoyancy parameter at
!                                     ! the top of the DSC
     &,BTC_DSCT(row_length,rows)                                        &
                                      ! Cloud fraction weighted buoyancy
!                                     ! parameter at the top of the DSC
     &,DB_DSCT_CLD(row_length,rows)                                     &
                                      ! In-cloud buoyancy jump at the
!                                     ! top of the DSC
     &,CHI_S_DSCT(row_length,rows)                                      &
                                      ! Mixing fraction of just saturate
!                                     ! mixture at top of the DSC
     &,ZETA_R_DSC(row_length,rows)                                      &
                                      ! Non-cloudy fraction of DSC
!                                     ! for cloud top radiative
!                                     ! cooling entrainment term.
     &,ZC_DSC(row_length,rows)                                          &
                                      ! Cloud depth (not cloud fraction
!                                     ! weighted).
     &,Z_CLD_DSC(row_length,rows)                                       &
                                      ! Cloud fraction weighted
!                                     ! thickness of DSC cloud
     &,CLD_FACTOR_DSC(row_length,rows)                                  &
                                      ! Fraction of grid box potentially
!                                     ! giving evaporative entrainment
     &,D_SIEMS(row_length,rows)                                         &
                                      ! Siems (1990) et al. cloud-top
!                                     !  entr.t instab. parm
     &,D_SIEMS_DSC(row_length,rows)                                     &
                                      ! Siems (1990) et al. cloud-top
!                                     ! entr.t instab. parm for DSC
                                      ! layer
     &,Z_TOP(row_length,rows,BL_LEVELS)                                 &
                                        ! Z_TOP(*,K) is the height of
!                                       ! level k+1/2 above the surface.
     &,TV1_SD(row_length,rows)                                          &
                                      ! Standard Deviation of level 1
!                                     ! virtual temperature (K).
     &,F_NET_SML(row_length,rows,BL_LEVELS)                             &
     &,F_NET_DSC(row_length,rows,BL_LEVELS)                             &
!                                     ! Net radiative fluxes relative
!                                     ! to SML and DSC base
     &,DF_NET_SML(row_length,rows)                                      &
                                      ! Net radiative divergences over
     &,DF_NET_DSC(row_length,rows)                                      &
                                      ! various layer depths
     &,CF_SML(row_length,rows)                                          &
                                      ! cloud fraction of SML
     &,CF_DSC(row_length,rows)                                          &
                                      ! cloud fraction of DSC layer
     &,Z_CF_BASE(row_length,rows)                                       &
                                      ! cloud base height from cld sch
     &,Z_CTOP(row_length,rows)                                          &
                                      ! cloud top height
     &,DQW_SML(row_length,rows)                                         &
                                      ! QW change across SML disc inv
     &,DTL_SML(row_length,rows)                                         &
                                      ! TL change across SML disc inv
     &,DQW_DSC(row_length,rows)                                         &
                                      ! QW change across DSC disc inv
     &,DTL_DSC(row_length,rows)                                         &
                                      ! TL change across DSC disc inv
     &,TLS_INV(row_length,rows)                                         &
                                      ! TL subs increment across inv
     &,RHOKH_SURF_ENT(row_length,rows)                                  &
                                      ! SML surf-driven entr. KH
     &,RHOKH_TOP_ENT(row_length,rows)                                   &
                                      ! SML top-driven entr. KH
     &,RHOKH_DSCT_ENT(row_length,rows)                                  &
                                      ! DSC top-driven entr. KH
     &,ZDSC_BASE(row_length,rows)                                       &
                                      ! Height of base of K_top in DSC
     &,WE_PARM(row_length,rows)                                         &
                                      ! param.d entrainment rates (m/s)
     &,WE_DSC_PARM(row_length,rows)                                     &
                                      ! for surf and DSC layers
     &,ZH_NP1(row_length,rows)                                          &
                                      ! estimate of ZH at end of timeste
     &,ZHSC_NP1(row_length,rows)      ! estimate of ZHSC at end of times

      REAL                                                              &
     & SCMOUT(row_length,rows)        ! For SCM diagnostics

      INTEGER                                                           &
     & K_CLOUD_TOP(row_length,rows)                                     &
                                      ! Level number of top of b.l.
!                                     ! cloud.
     &,K_CLOUD_DSCT(row_length,rows)                                    &
                                      ! Level number of top of dec.
!                                     ! cloud.
     &,NTML_SAVE(row_length,rows)                                       &
                                      ! Copy of NTML
     &,NTML_PREV(row_length,rows)                                       &
                                      ! NTML from previous timestep
     &,K_PLUME(row_length,rows)                                         &
                                      ! Start grid-level for
                                      ! surface-driven plume
     &,K_CBASE(row_length,rows)       ! grid-level above cloud-base

! NEC vectorization
      INTEGER                                                           &
     & K_LEVEL(row_length,rows)                                         &
                                      ! array to store level selection
     &,K_CFF(row_length,rows)         ! level counter for CFF

!  (b) Scalars.
!
      REAL                                                              &
     & VIRT_FACTOR                                                      &
                         ! Temporary in calculation of buoyancy
!                        ! parameters.
     &,DTLDZ                                                            &
                         ! Vertical gradient of TL in a well-mixed
!                        ! layer.
     &,DQW                                                              &
                         ! Total water content change across layer
!                        ! interface.
     &,DTL                                                              &
                         ! Liquid/ice static energy change across
!                        ! layer interface.
     &,DTL_GA                                                           &
                         ! As DTL but inc gradient adjustment
     &,DQW_GA                                                           &
                         ! As DQW but inc gradient adjustment
     &,DQCL                                                             &
                         ! Cloud liquid water change across layer
!                        ! interface.
     &,DQCF                                                             &
                         ! Cloud frozen water change across layer
!                        ! interface.
     &,Q_VAP_PARC                                                       &
                         ! Vapour content of parcel
     &,Q_LIQ_PARC                                                       &
                         ! Condensed water content of parcel
     &,Q_LIQ_ENV                                                        &
                         ! Condensed water content of environment
     &,T_PARC                                                           &
                         ! Temperature of parcel
     &,T_DENS_PARC                                                      &
                         ! Density potential temperature of parcel
     &,T_DENS_ENV                                                       &
                         ! Density potential temperature of
!                        ! environment
     &,DENV_BYDZ                                                        &
                         ! Gradient of density potential
!                        ! temperature in environment
     &,DPAR_BYDZ                                                        &
                         ! Gradient of density potential
!                        ! temperature of parcel
     &,RHO_DZ                                                           &
                         ! rho*dz
     &,SVL_LAPSE                                                        &
                      ! Lapse rate of SVL above inversion (K/m)
     &,SVL_LAPSE_BASE                                                   &
                      ! Lapse rate of SVL above inversion (K/m)
     &,SVL_TOP                                                          &
                   ! s_VL at half level above inversion (K)
     &,DSVL_TOP                                                         &
                   ! s_VL jump across inversion grid layer (K)
     &,TOTHF_INV                                                        &
                   ! total heat flux at discontinous inversion height
     &,DF_INV                                                           &
                   ! temporary in rad divergence calculation
     &,DZ_DISC_MIN                                                      &
                   ! smallest allowed DZ_DISC
     &,TL_K                                                             &
                   ! TL (full static energy) on level K
     &,TL_KP1                                                           &
                   ! TL (full static energy) on level K+1
     &,TL_KP2                                                           &
                   ! TL (full static energy) on level K+2
     &,DB_DISC                                                          &
                   ! Temporary disc inversion buoyancy jump
     &,W_S_ENT                                                          &
                   ! numerical (subsidence) entrainment rate
     &,W_LS                                                             &
                   ! large-scale (subs) velocity
     &,W_LS_DSC                                                         &
                   ! large-scale (subs) velocity
     &,DZ_DISC                                                          &
                   ! height of ZH below Z_HALF(NTML+2)
     &,Z_SURF                                                           &
                   ! approx height of top of surface layer
     &,QUAD_A                                                           & 
                   ! term `a' in quadratic solver for DZ_DISC
     &,QUAD_BM                                                          & 
                   ! term `-b'in quadratic solver for DZ_DISC
     &,QUAD_C                                                           & 
                   ! term `c' in quadratic solver for DZ_DISC
     &,W_M                                                              &
                   ! scaling velocity for standard deviations
     &,W_S_CUBED                                                        &
                   ! convective velocity scale
     &,Z_CBASE                                                          &
                   ! cloud base height (m)
     &,ZDSC_CBASE                                                       &
                   ! DSC cloud base height (m)
     &,Z_INT                                                            &
                   ! depth of wb integral
     &,Z_INT_TOP                                                        &
                   ! top of wb-integration depth
     &,CF_FOR_WB                                                        &
                   ! CF for use in wb calculation for decoupling
     &,DFSW_TOP    ! SW radiative flux change assoc with cloud-top
!
      INTEGER                                                           &
     & I                                                                &
               ! Loop counter (horizontal field index).
     &,J                                                                &
               ! Offset counter in certain I loops.
     &,K                                                                &
               ! Loop counter (vertical level index).
     &,KL                                                               &
               ! K
     &,KM1                                                              &
               ! K-1
     &,KP1                                                              &
               ! K+1
     &,KP2                                                              &
               ! K+2
     &,MBL                                                              &
               ! Maximum number of model layers allowed in the rapidly
!              ! mixing layer; set to BL_LEVELS-1.
     &,K_RAD_SMLT                                                       &
                  ! highest SML level for radiative
!                 !   divergence calculation
     &,K_RAD_LIM                                                        &
                  !
     &,IENT       ! loop counter

!-----------------------------------------------------------------------
!!  0.  Check that the scalars input to define the grid are consistent.
!!      See comments to routine SF_EXCH for details.
!-----------------------------------------------------------------------

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('KMKHZ   ',3)
      END IF

!  Set MBL, "maximum number of boundary levels" for the purposes of
!  boundary layer height calculation.

      MBL = BL_LEVELS - 1

!-----------------------------------------------------------------------
! Calculate Z_TOP (top of levels) and NTML from previous timestep
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
        NTML_PREV(I,j) = 1
      END DO
      END DO
      DO K = 1,BL_LEVELS-1
        do j=1,rows
        do i=1,row_length
          Z_TOP(I,j,K) = Z_HALF(I,j,K+1)
!         !------------------------------------------------------------
!         !find NTML from previous TS (for accurate gradient adjustment
!         !of profiles - also note that NTML LE BL_LEVELS-1)
!         !------------------------------------------------------------
          IF ( ZH_PREV(I,j)  >=  Z_HALF(I,j,K+1) ) NTML_PREV(I,j)=K
        END DO
        END DO
      END DO
      K = BL_LEVELS
      do j=1,rows
      do i=1,row_length
        Z_TOP(I,j,K) = Z_HALF(I,j,K) + DZL(I,j,K)
      END DO
      END DO
!-----------------------------------------------------------------------
! Calculate SVL: conserved variable used to test for well mixed layers
!-----------------------------------------------------------------------
      DO K = 1,BL_LEVELS
        do j=1,rows
        do i=1,row_length
         SVL(I,j,K) = ( TL(I,j,K) + GRCP * Z_FULL(I,j,K) )              &
     &              * ( 1.0 + C_VIRTUAL*QW(I,j,K) )
        END DO
        END DO
      END DO
!
      DO  K = 1,BL_LEVELS
! DEPENDS ON: qsat_mix
        CALL QSAT_mix(QS(1,1,K),T(1,1,K),P(1,1,K),row_length*rows       &
     & ,lq_mix_bl)
      END DO
!--------------------------------------------------------------------
!..Calculate subsidence increments:
!--------------------------------------------------------------------
      DO K = 1,BL_LEVELS
      do j=1,rows
      do i=1,row_length

        IF ( W(I,J,K)  <=  0.0 ) THEN
!         ! only needed in subsidence regions
          KP1 = MIN( BL_LEVELS, K+1 )
          TLS_INC(I,J,K) = - W(I,J,K) *RDZ(I,J,KP1)                     &
     &                  * ( TL(I,J,KP1) - TL(I,J,K) +                   &
     &                      G/( CP*RDZ(I,J,KP1) ) )
          QLS_INC(I,J,K) = - W(I,J,K) *RDZ(I,J,KP1)                     &
     &                                    * ( QW(I,J,KP1) - QW(I,J,K) )
        ELSE
          TLS_INC(I,J,K) = 0.0
          QLS_INC(I,J,K) = 0.0
        END IF

      END DO
      END DO
      END DO


      DO J=1,ROWS
      DO I=1,ROW_LENGTH
        UNSTABLE(i,j) = (FB_SURF(i,j) >  0.0)
        K_PLUME(i,j)  = -1
      END DO
      END DO
!         !------------------------------------------------------------
!         ! Find grid-level above top of surface layer, taken
!         ! to be at a height, z_surf, given by:
!         !       Z_SURF = 0.1*ZH_PREV
!         ! Use ZH_prev since that will have determined the shape
!         ! of the time-level n profiles.
!         !------------------------------------------------------------
      DO K=1,BL_LEVELS-1
        DO J=1,ROWS
        DO I=1,ROW_LENGTH

          IF ( UNSTABLE(I,j) ) THEN

            Z_SURF = 0.1 * ZH_PREV(I,j)
            IF ( Z_FULL(I,j,K) >= Z_SURF .AND. K_PLUME(i,j) == -1 ) THEN
                 !reached z_surf
              K_PLUME(i,j)=K
            END IF
            IF ( SVL(I,j,K+1) >= SVL(I,j,K)                             &
     &              .AND. K_PLUME(i,j) == -1 ) THEN
                 !reached inversion
              K_PLUME(i,j)=K
            END IF

          END IF
        END DO
        END DO
      END DO

      DO J=1,ROWS
      DO I=1,ROW_LENGTH
        IF (K_PLUME(i,j) == -1) K_PLUME(I,J)=1
      END DO
      END DO
!-----------------------------------------------------------------------
!! 0.3a IF NOT CUMULUS:
!L     Look for decoupled cloudy mixed-layer above b.l. top at ZH
!L     (starting from level 3 and below 2.5km):
!L     find cloud-base above SML inversion, ie. above NTML+1,
!L     then cloud-top (ie. CF < SC_CFTOL)
!L     and finally check that cloud is well-mixed.
!-----------------------------------------------------------------------
!      Initialise variables
!
      do j=1,rows
      do i=1,row_length
        CLOUD_BASE(I,j) = .FALSE.
        DSC(I,j) = .FALSE.
        COUPLED(I,j) = .FALSE.
        ZHSC(I,j)    = 0.0
        NTDSC(I,j)   = 0
      END DO
      END DO
!
      DO K=3,MBL
        do j=1,rows
        do i=1,row_length
!----------------------------------------------------------------------
!..Find cloud-base (where cloud here means CF > SC_CFTOL)
!----------------------------------------------------------------------
          IF ( .NOT. CUMULUS(I,j) .AND.                                 &
     &         Z_FULL(I,j,K) < ZMAX_FOR_DSC .AND.                       &
     &         K  >   NTML(I,j)+1 .AND. CF(I,j,K)  >   SC_CFTOL         &
     &                          .AND. .NOT.CLOUD_BASE(I,j)              &
!                                  not yet found cloud-base
     &                          .AND. .NOT.DSC(I,j) ) THEN
!                                  not yet found a Sc layer
            CLOUD_BASE(I,j) = .TRUE.
          END IF
          IF ( CLOUD_BASE(I,j) .AND. .NOT.DSC(I,j) .AND. (              &
!                  found cloud-base but not yet reached cloud-top
     &                            (CF(I,j,K+1) <  SC_CFTOL) )           &
!                  got to cloud-top
     &       ) THEN
            CLOUD_BASE(I,j) = .FALSE.         ! reset CLOUD_BASE
!           !-----------------------------------------------------------
!           ! Look to see if at least top of cloud is well mixed:
!           ! test SVL-gradient for top 2 pairs of levels, in case
!           ! cloud top extends into the inversion.
!           ! Parcel descent in Section 4.0 below will determine depth
!           ! of mixed layer.
!           !----------------------------------------------------------
           IF ( (SVL(I,j,K)-SVL(I,j,K-1))                               &
     &                  /(Z_FULL(I,j,K)-Z_FULL(I,j,K-1))                &
     &                                         <   MAX_SVL_GRAD ) THEN
            DSC(I,j) = .TRUE.
            NTDSC(I,j) = K
            ZHSC(I,j)  = Z_HALF(I,j,NTDSC(I,j)+1)
           ELSEIF ( (SVL(I,j,K-1)-SVL(I,j,K-2))                         &
     &                  /(Z_FULL(I,j,K-1)-Z_FULL(I,j,K-2))              &
     &                                      <   MAX_SVL_GRAD ) THEN
!             !---------------------------------------------------------
!             ! Well-mixed layer with top at k-1 or k.  Check whether
!             ! there is a buoyancy inversion between levels k-1 and k
!             ! in a manner similar to the surface-driven plume: compare
!             ! the buoyancy gradient between levels K-1 and K for an
!             ! undiluted parcel and the environment
!             !---------------------------------------------------------
              SL_PLUME(I,j) = TL(I,j,K-1) + GRCP * Z_FULL(I,j,K-1)
              QW_PLUME(I,j) = QW(I,j,K-1)
! -------------------------------------------------------------------
! calculate parcel water by linearising qsat about the environmental
! temperature.
! -------------------------------------------------------------------
              IF(T(I,j,K) >  TM) THEN
                Q_LIQ_PARC = MAX( 0.0, ( QW_PLUME(I,j) - QS(I,j,K) -    &
     &            DQSDT(I,j,K)*                                         &
     &            ( SL_PLUME(I,j)-GRCP*Z_FULL(I,j,K)-T(I,j,K) )         &
     &                                 ) *A_QS(I,j,K) )
                Q_LIQ_ENV = MAX( 0.0, ( QW(I,j,K) - QS(I,j,K)           &
     &                      -DQSDT(I,j,K)*( TL(I,j,K) - T(I,j,K) )      &
     &                                 ) *A_QS(I,j,K) )
! add on the difference in the environment's ql as calculated by the
! partial condensation scheme (using some RH_CRIT value) and what it
! would be if RH_CRIT=1. This then imitates partial condensation
! in the parcel.
                Q_LIQ_PARC = Q_LIQ_PARC + QCL(I,j,K) + QCF(I,j,K)       &
     &                         - Q_LIQ_ENV
                T_PARC = SL_PLUME(I,j) - GRCP * Z_FULL(I,j,K) +         &
     &                           LCRCP*Q_LIQ_PARC
              ELSE
                Q_LIQ_PARC = MAX( 0.0, ( QW_PLUME(I,j) - QS(I,j,K) -    &
     &            DQSDT(I,j,K)*                                         &
     &              ( SL_PLUME(I,j)-GRCP*Z_FULL(I,j,K)-T(I,j,K) )       &
     &                                 ) *A_QS(I,j,K) )
                Q_LIQ_ENV = MAX( 0.0, ( QW(I,j,K) - QS(I,j,K)           &
     &             -DQSDT(I,j,K)*( TL(I,j,K) - T(I,j,K) )               &
     &                                 ) *A_QS(I,j,K) )
! add on difference in environment's ql between RH_CRIT and RH_CRIT=1
                Q_LIQ_PARC = Q_LIQ_PARC + QCL(I,j,K) + QCF(I,j,K)       &
     &                         - Q_LIQ_ENV
                T_PARC = SL_PLUME(I,j) - GRCP * Z_FULL(I,j,K) +         &
     &                           LSRCP*Q_LIQ_PARC
              END IF
              Q_VAP_PARC=QW_PLUME(I,j)-Q_LIQ_PARC
!
              T_DENS_PARC=T_PARC*(1.0+C_VIRTUAL*Q_VAP_PARC-Q_LIQ_PARC)
              T_DENS_ENV=T(I,j,K)*                                      &
     &                   (1.0+C_VIRTUAL*Q(I,j,K)-QCL(I,j,K)-QCF(I,j,K))
! find vertical gradients in parcel and environment SVL (using values
! from level below (K-1))
              ENV_SVL_KM1(I,j) = T(I,j,K-1) * ( 1.0+C_VIRTUAL*Q(I,j,K-1)&
     &             -QCL(I,j,K-1)-QCF(I,j,K-1) ) + GRCP*Z_FULL(I,j,K-1)
              DPAR_BYDZ=(T_DENS_PARC+GRCP*Z_FULL(I,j,K)-                &
     &                    ENV_SVL_KM1(I,j)) /                           &
     &                (Z_FULL(I,j,K)-Z_FULL(I,j,K-1))
              DENV_BYDZ=(T_DENS_ENV+GRCP*Z_FULL(I,j,K)-                 &
     &                    ENV_SVL_KM1(I,j))/                            &
     &                (Z_FULL(I,j,K)-Z_FULL(I,j,K-1))
!
              IF ( DENV_BYDZ >  1.25*DPAR_BYDZ ) THEN
!               ! there is an inversion between levels K-1 and K
                IF ( K  >=  NTML(I,j)+3 ) THEN
!                 ! if NTDSC EQ NTML+1 then assume we're looking
!                 ! at the same inversion and so don't set DSC
                  NTDSC(I,j) = K-1
                  ZHSC(I,j)  = Z_HALF(I,j,NTDSC(I,j)+1)
                  DSC(I,j) = .TRUE.
                END IF
              ELSE
!               ! no inversion between levels K-1 and K, assume there
!               ! is an inversion between K and K+1 because of CF change
                NTDSC(I,j) = K
                ZHSC(I,j)  = Z_HALF(I,j,NTDSC(I,j)+1)
                DSC(I,j) = .TRUE.
              END IF
           END IF
          END IF
        END DO
        END DO
      END DO
!-----------------------------------------------------------------------
!! 0.4a If the layer to ZHPAR is a cumulus layer capped by cloud and
!!      an inversion, declare this layer a decoupled cloud layer and
!!      set ZHSC and NTDSC accordingly.
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
        IF ( (L_conv4a .AND.                                            &
     &          L_SHALLOW(I,j) .AND. NTPAR(I,j)  <   BL_LEVELS )        &
!               ! shallow cumulus layer within BL_LEVELS (A05_4A)
     &       .OR. (.NOT. L_conv4a .AND.                                 &
     &            CUMULUS(I,j) .AND. NTPAR(I,j)  <   BL_LEVELS ) ) THEN
!               ! cumulus layer and inversion found
          IF ( CF(I,j,NTPAR(I,j))  >   SC_CFTOL  .OR.                   &
     &         CF(I,j,NTPAR(I,j)+1)  >   SC_CFTOL ) THEN
!            ! cloudy
            DSC(I,j)  = .TRUE.
            ZHSC(I,j) = ZHPAR(I,j)
            NTDSC(I,j)= NTPAR(I,j)
          END IF
        END IF
      END DO
      END DO
!-----------------------------------------------------------------------
! 0.4b Calculate the radiative flux changes across cloud top for the
!      stratocumulus layer and thence a first guess for the top-down
!      mixing depth of this layer, DSCDEPTH.
!-----------------------------------------------------------------------
!     Initialise variables
!------------------------------
      do j=1,rows
      do i=1,row_length
        K_CLOUD_DSCT(I,j) = 0
        DF_DSCT_OVER_CP(I,j) = 0.0
      END DO
      END DO
!
      DO K = 1,BL_LEVELS
        do j=1,rows
        do i=1,row_length

          DF_OVER_CP(I,j,K) = - ( RAD_HR(I,j,K,1) + RAD_HR(I,j,K,2) )   &
     &                        * RHO_DRY_TQ(I,j,K) * DZL(I,j,K)
          DFLW_OVER_CP(I,j,K) = - RAD_HR(I,j,K,1)                       &
     &                        * RHO_DRY_TQ(I,j,K) * DZL(I,j,K)
          DFSW_OVER_CP(I,j,K) = - RAD_HR(I,j,K,2)                       &
     &                        * RHO_DRY_TQ(I,j,K) * DZL(I,j,K)
!         !-------------------------------------------------------------
!         ! Find the layer with the greatest LW radiative flux jump and
!         ! assume that this is the top DSC level.  Limit the
!         ! search to above the SML.
!         !-------------------------------------------------------------
          K_RAD_LIM = NTML(I,j)+2

          IF ( DSC(I,j) .AND. K  >=  K_RAD_LIM .AND. K  <=  NTDSC(I,j)+2&
     &        .AND. DFLW_OVER_CP(I,j,K)  >   DF_DSCT_OVER_CP(I,j) ) THEN
            K_CLOUD_DSCT(I,j) = K
            DF_DSCT_OVER_CP(I,j) = DFLW_OVER_CP(I,j,K)
          END IF

        END DO
        END DO
      END DO
!     !-----------------------------------------------------------------
!     !  In case the LW radiative divergence is spread over 2 or 3 levs,
!     !  add on any cooling either side of the level of maximum cooling.
!     !-----------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
          IF ( DSC(I,j) .AND. K_CLOUD_DSCT(I,j)  >   0 ) THEN
            K=K_CLOUD_DSCT(I,j)
            DFSW_TOP = DFSW_OVER_CP(I,j,K)
            IF ( K >  1 ) THEN
              DF_DSCT_OVER_CP(I,j) = DF_DSCT_OVER_CP(I,j)               &
     &                                + MAX(0.0, DFLW_OVER_CP(I,j,K-1) )
              DFSW_TOP = DFSW_TOP + MIN(0.0, DFSW_OVER_CP(I,j,K-1) )
            END IF
            IF ( K <  BL_LEVELS ) THEN
              DF_DSCT_OVER_CP(I,j) = DF_DSCT_OVER_CP(I,j)               &
     &                                + MAX(0.0, DFLW_OVER_CP(I,j,K+1) )
              DFSW_TOP = DFSW_TOP + MIN(0.0, DFSW_OVER_CP(I,j,K+1) )
            END IF
!           !-----------------------------------------------------------
!           ! Combine SW and LW cloud-top divergences into a net
!           ! divergence by estimating SW flux divergence at a given
!           ! LW divergence = DF_SW * (1-exp{-A*kappa_sw/kappa_lw})
!           ! Empirically (from LEM data) a reasonable fit is found
!           ! with A small and (1-exp{-A*kappa_sw/kappa_lw}) = 0.35
!           !-----------------------------------------------------------
              DF_DSCT_OVER_CP(I,j) = MAX( 0.0,                          &
     &                  DF_DSCT_OVER_CP(I,j) + 0.35 * DFSW_TOP )
          END IF
      END DO
      END DO
!-----------------------------------------------------------------------
!     Set NBDSC, the bottom level of the DSC layer.
!     Note that this will only be used to give an estimate of the layer
!     depth, DSCDEPTH, used to calculate the entrainment
!     rate (the dependence is only weak), and that a more accurate
!     algorithm is subsequently used to determine the depth over which
!     the top-down mixing profiles will be applied.  If DSC is FALSE,
!     DSCDEPTH = 0.  The plume descent here uses a radiative
!     perturbation to the cloud-layer SVL (use level NTDSC-1 in case
!     SVL is not yet well-mixed to NTDSC), based roughly
!     on a typical cloud-top residence time.  If the plume does not sink
!     and the cloud is decoupled from the surface (ie. above Alan
!     Grant's ZH), then it is assumed to be stable, ie. St rather than
!     Sc, and no mixing or entrainment is applied to it.
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
       NBDSC(I,j) = NTDSC(I,j)+1
       IF (DSC(I,j)) THEN
! The depth of the radiatively-cooled layer tends to be less than O(50m)
! and so RAD_HR will be an underestimate of the cooling tendency there.
! Compensate by multiplying by DZL/50. (~4)
! Recall that DF_OVER_CP(I,j,K) = RAD_HR * RHO_DRY_TQ * DZL
! Thus use cloud-top radiative forcing as follows:

        K = NTDSC(I,j)
        RHO_DZ = RHO_DRY_TQ(I,j,K) * DZL(I,j,K)
        SVL_PLUME(I,j)=SVL(I,j,K-1)                                     &
     &     - CT_RESID * DZL(I,j,K)*DF_DSCT_OVER_CP(I,j) / ( 50.*RHO_DZ )

       ELSE
        SVL_PLUME(I,j)=0.0
       END IF
      END DO
      END DO

      DO K=BL_LEVELS-1,1,-1
      do j=1,rows
      do i=1,row_length
       IF ( K <  NTDSC(I,j) .AND. SVL_PLUME(I,j)  <   SVL(I,j,K) ) THEN
         NBDSC(I,j) = K+1     ! marks lowest level within ML
       END IF
      END DO
      END DO
      END DO
!----------------------------------------------------------------------
!  0.4e Tidy up variables associated with decoupled layer
!       NOTE that NTDSC GE 3 if non-zero
!----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
! Note that ZHSC-Z_HALF(NTML+2) may = 0, so this test comes first!
       IF (CUMULUS(I,j) .AND. DSC(I,j))                                 &
     &                  NBDSC(I,j) = MAX( NBDSC(I,j), NTML(I,j)+2 )
       IF ( NTDSC(I,j)  >=  1 ) THEN
        IF ( NBDSC(I,j) <   NTDSC(I,j)+1 ) THEN
          DSCDEPTH(I,j) =                                               &
     &                Z_HALF(I,j,NTDSC(I,j)+1) - Z_HALF(I,j,NBDSC(I,j))
        ELSE
!        !----------------------------------------------------------
!        ! Indicates a layer of zero depth
!        !----------------------------------------------------------
         IF (NTDSC(I,j) == NTPAR(I,j)) THEN
!           !----------------------------------------------------------
!           ! Indicates a Sc layer at the top of Cu: force mixing
!           ! over single layer.
!           !----------------------------------------------------------
           DSCDEPTH(I,j) = DZL(I,j,NTDSC(I,j))
         ELSE
           DSC(I,j)=.FALSE.
           NTDSC(I,j)=0
           ZHSC(I,j)=0.0
           DF_DSCT_OVER_CP(I,j) = 0.0
           DSCDEPTH(I,j) = 0.0
         END IF
       END IF
       ELSE  ! NTDSC EQ 0, just to make sure!
        DSCDEPTH(I,j)=0.0
        DSC(I,j)=.FALSE.
        ZHSC(I,j)=0.0
        DF_DSCT_OVER_CP(I,j) = 0.0
       END IF
      END DO
      END DO
!----------------------------------------------------------------------
!..If decoupled cloud-layer found test to see if it is, in fact,
!..only weakly decoupled from the surface mixed-layer:
!..if SVL difference between NTML and NTDSC is less than SVL_COUP Kelvin
!..then assume there is still some coupling.  This will mean that
!..the surface-driven entrainment term will be applied at ZHSC, no
!..subgrid inversion or entrainment will be calculated for ZH and
!..ZHSC will be the length scale used in the entrainment inputs.
!..Note that for CUMULUS surface-driven entrainment will be done
!..by the convection scheme.
!----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
        COUPLED(I,j) = .FALSE.
        IF ( DSC(I,j) .AND. .NOT.CUMULUS(I,j) ) THEN
!         !------------------------------------------------------------
!         ! Note this IF test structure is required because if DSC is
!         ! false then NTDSC = 0 and cannot be used to index SVL.
!         !------------------------------------------------------------
          IF ( SVL(I,j,NTDSC(I,j)) - SVL(I,j,NTML(I,j))  <   SVL_COUP ) &
     &       COUPLED(I,j) = .TRUE.
        END IF
      END DO
      END DO
!-----------------------------------------------------------------------
! Assuming a discontinuous inversion structure.
!    - to this point in the code, ZH and ZHSC mark the half-level at
!      the base of the inversion
!    - now they will be interpolated into the level above assuming
!      SVL(NTML+1) is a volume average over a subgrid discontinuous
!      inversion structure
!    - discontinuous jumps of TL and QW (and thence buoyancy) can be
!      calculated and used to determine the entrainment rate
!    - parametrized grid-level fluxes at NTML,NTDSC can then be made
!      consistent with this assumed inversion structure
!-----------------------------------------------------------------------
!      ! If any `problems' are encountered with this interpolation of ZH
!      ! (such as ZH diagnosed at or below Z_HALF(NTML+1)), then NTML
!      ! is lowered a level and ZH is set fractionally below what has
!      ! become Z_HALF(NTML+2).  This distance is such that for a net
!      ! dZH/dt of 1.E-4 m/s, ZH will be diagnosed as spending at least
!      ! half the timestep in level NTML+2, leaving the growth only
!      ! marginally affected.  Conversely, it allows a subsiding
!      ! inversion to fall more readily.

      DZ_DISC_MIN = 0.5 * TIMESTEP * 1.E-4
         do j=1,rows
         do i=1,row_length

        SML_DISC_INV(I,j) = 0  ! initialise flags to indicate whether a
        DSC_DISC_INV(I,j) = 0  ! discontinuous inversion is diagnosed
!..First interpolate to find ZH

        K = NTML(I,j)
!..by default, keep ZH at the half-level where it was diagnosed
!..initially and use grid-level jumps

        DTL_SML(I,j) = TL(I,j,K+1) - TL(I,j,K) + GRCP/RDZ(I,j,K+1)
        DQW_SML(I,j) = QW(I,j,K+1) - QW(I,j,K)

        IF ( .NOT.CUMULUS(I,j) .AND. .NOT.COUPLED(I,j) .AND.            &
     &                K  >   1 .AND. K  <=  BL_LEVELS-2 ) THEN
          IF ( SVL(I,j,K+2)  >   SVL(I,j,K+1)                           &
     &                  .AND. SVL(I,j,K+1)  >   SVL(I,j,K) ) THEN

              TL_K = TL(I,j,K) + GRCP * Z_FULL(I,j,K)
              TL_KP1 = TL(I,j,K+1) + GRCP * Z_FULL(I,j,K+1)
              TL_KP2 = TL(I,j,K+2) + GRCP * Z_FULL(I,j,K+2)

              IF ( K  <=  BL_LEVELS-3 ) THEN
! need to test for K+1 to K+2 gradient in case profile is concave
! (would mess up the inversion diagnosis)
                SVL_LAPSE = MAX(0.0,                                    &
     &            MIN( ( SVL(I,j,K+2) - SVL(I,j,K+1) ) * RDZ(I,j,K+2),  &
     &                 ( SVL(I,j,K+3) - SVL(I,j,K+2) ) * RDZ(I,j,K+3)   &
     &               )         )
              ELSE
                SVL_LAPSE = 0.0
              END IF
              IF ( K  >=  K_PLUME(I,j)+2 ) THEN
!               ! Use mean mixed layer gradient (if resolved) to allow
!               ! for stablisation by gradient-adjustment
!               ! Ignore level K in case inversion is dropping
                SVL_LAPSE_BASE = ( SVL(I,j,K-1)-SVL(I,j,K_PLUME(I,j)) )/&
     &                        (Z_FULL(I,j,K-1)-Z_FULL(I,j,K_PLUME(I,j)))
                SVL_LAPSE_BASE = MAX( 0.0, SVL_LAPSE_BASE )
              ELSE
                SVL_LAPSE_BASE = 0.0
              END IF

              QUAD_A  = 0.5*( SVL_LAPSE - SVL_LAPSE_BASE )
              QUAD_BM = SVL(I,j,K+2) - SVL(I,j,K)                       &
     &            - SVL_LAPSE * ( Z_FULL(I,j,K+2)-Z_HALF(I,j,K+2) )     &
     &            - SVL_LAPSE_BASE * ( Z_HALF(I,j,K+1)-Z_FULL(I,j,K) +  &
     &                                                    DZL(I,j,K+1) )
              QUAD_C  = DZL(I,j,K+1)*( SVL(I,j,K+1) - SVL(I,j,K) -      &
     &            SVL_LAPSE_BASE * (                                    &
     &              Z_HALF(I,j,K+1)-Z_FULL(I,j,K) + 0.5*DZL(I,j,K+1) ) )

              IF ( QUAD_BM  >   0.0 ) THEN
                IF ( QUAD_C  <=  0.0) THEN
!                 ! SVL extrapolated from K to K+1 is greater than
!                 ! the level K+1 value - inversion needs to rise so
!                 ! place it as high as possible
                  DZ_DISC = DZ_DISC_MIN
                ELSEIF ( QUAD_BM*QUAD_BM  >=  4.0*QUAD_A*QUAD_C ) THEN
!                 ! solve equation for DZ_DISC...
                  IF ( QUAD_A  /=  0.0 ) THEN
!                   !   ...quadratic if QUAD_A NE 0
                    DZ_DISC = ( QUAD_BM - SQRT( QUAD_BM*QUAD_BM         &
     &                                       - 4.0*QUAD_A*QUAD_C )      &
     &                              ) / (2.0*QUAD_A)
                  ELSE
!                   !   ...linear if QUAD_A EQ 0
                    DZ_DISC = QUAD_C / QUAD_BM
                  END IF
                ELSE
                 DZ_DISC = 99999.9  ! large dummy value
                END IF

                IF ( DZ_DISC  >   0.9 * DZL(I,j,K+1) ) THEN
!               ! ZH diagnosed very close to or below Z_HALF(K+1):
                  IF ( SVL(I,j,K)-SVL(I,j,K-1)  >   0.0) THEN
!                   ! top of ML stably stratified so lower NTML but
!                   ! set ZH only fractionally (DZ_DISC_MIN)
!                   ! below the top of the inversion level.
                    NTML(I,j) = NTML(I,j) - 1
                    K=NTML(I,j)
                    TL_K = TL(I,j,K) + GRCP * Z_FULL(I,j,K)
                    TL_KP1 = TL(I,j,K+1) + GRCP * Z_FULL(I,j,K+1)
                    TL_KP2 = TL(I,j,K+2) + GRCP * Z_FULL(I,j,K+2)
                    DZ_DISC = DZ_DISC_MIN
                  ELSE
!                   ! top of ML well-mixed so don't lower the inversion
!                   ! level but set ZH just (DZ_DISC_MIN) above the
!                   ! half-level to allow the inversion to subside if
!                   ! necessary.
                    DZ_DISC = DZL(I,j,K+1) - DZ_DISC_MIN
                  END IF
                END IF

              ELSE
!.. ignoring lapse rates
                DSVL_TOP = SVL(I,j,K+2) - SVL(I,j,K)
                DZ_DISC = DZL(I,j,K+1) *                                &
     &                          (SVL(I,j,K+1)-SVL(I,j,K)) / DSVL_TOP
              END IF

            ZH(I,j) = Z_HALF(I,j,K+2) - DZ_DISC
            SML_DISC_INV(I,j) = 1 ! set flag to indicate disc inv found

!..Calculate discontinuous jumps of TL and QW:
!..only accurate enough for large enough DZ_DISC

            IF (DZ_DISC/DZL(I,j,K+1)  >   0.1) THEN
              DTL_SML(I,j) = (TL_KP1-TL_K) * DZL(I,j,K+1) / DZ_DISC

              IF ( TL_KP2  >   TL_KP1 .AND.                             &
     &                           TL_KP1  >   TL_K ) THEN
                DTL_SML(I,j) = MIN( TL_KP2-TL_K, DTL_SML(I,j) )
              ELSEIF ( TL_KP2  <   TL_KP1 .AND.                         &
     &                           TL_KP1  <   TL_K ) THEN
                DTL_SML(I,j) = MAX( TL_KP2-TL_K, DTL_SML(I,j) )
              ELSE  ! TL non-monotonic
                DTL_SML(I,j) = TL_KP1 - TL_K
              END IF

              DQW_SML(I,j) = (QW(I,j,K+1)-QW(I,j,K))                    &
     &                                     * DZL(I,j,K+1) / DZ_DISC

              IF ( QW(I,j,K+2)  >   QW(I,j,K+1) .AND.                   &
     &               QW(I,j,K+1)  >   QW(I,j,K) ) THEN
                DQW_SML(I,j) = MIN( QW(I,j,K+2)-QW(I,j,K), DQW_SML(I,j))
              ELSEIF ( QW(I,j,K+2)  <   QW(I,j,K+1) .AND.               &
     &               QW(I,j,K+1)  <   QW(I,j,K) ) THEN
                DQW_SML(I,j) = MAX( QW(I,j,K+2)-QW(I,j,K), DQW_SML(I,j))
              ELSE  ! QW non-monotonic
                DQW_SML(I,j) = QW(I,j,K+1)-QW(I,j,K)
              END IF

            ELSE
!..DZ_DISC small, ZH close to level above: use double grid-level jumps
              DTL_SML(I,j) = TL_KP2 - TL_K
              DQW_SML(I,j) = QW(I,j,K+2) - QW(I,j,K)
            END IF

          END IF  ! SVL increasing
         END IF ! not cumulus and not at top of bl_levels
!-----------------------------------------------------------------------
!..Second interpolate to find ZHSC
!-----------------------------------------------------------------------
         IF ( DSC(I,j) ) THEN
          K = NTDSC(I,j)
!..by default, keep ZHSC at the half-level where it was diagnosed
!..initially and use grid-level jumps
          DTL_DSC(I,j) = TL(I,j,K+1) - TL(I,j,K) + GRCP/RDZ(I,j,K+1)
          DQW_DSC(I,j) = QW(I,j,K+1) - QW(I,j,K)
          IF ( K  <=  BL_LEVELS-2 ) THEN
          IF ( SVL(I,j,K+2)  >   SVL(I,j,K+1)                           &
     &             .AND. SVL(I,j,K+1)  >   SVL(I,j,K) ) THEN

              TL_K = TL(I,j,K) + GRCP * Z_FULL(I,j,K)
              TL_KP1 = TL(I,j,K+1) + GRCP * Z_FULL(I,j,K+1)
              TL_KP2 = TL(I,j,K+2) + GRCP * Z_FULL(I,j,K+2)

              IF ( K  <=  BL_LEVELS-3 ) THEN
! need to test for K+1 to K+2 gradient in case profile is concave
! (would mess up the inversion diagnosis)
                SVL_LAPSE = MAX( 0.0,                                   &
     &            MIN( ( SVL(I,j,K+2) - SVL(I,j,K+1) ) * RDZ(I,j,K+2),  &
     &                 ( SVL(I,j,K+3) - SVL(I,j,K+2) ) * RDZ(I,j,K+3)   &
     &               )         )
              ELSE
                SVL_LAPSE = 0.0
              END IF
              IF ( K  >=  NBDSC(I,j)+2 ) THEN
!               ! Use mean mixed layer gradient (if resolved) to allow
!               ! for stablisation by gradient-adjustment
!               ! Ignore level K in case inversion is dropping
                SVL_LAPSE_BASE = ( SVL(I,j,K-1)-SVL(I,j,NBDSC(I,j)) )/  &
     &                        (Z_FULL(I,j,K-1)-Z_FULL(I,j,NBDSC(I,j)))
                SVL_LAPSE_BASE = MAX( 0.0, SVL_LAPSE_BASE )
              ELSE
                SVL_LAPSE_BASE = 0.0
              END IF

              QUAD_A  = 0.5*( SVL_LAPSE - SVL_LAPSE_BASE )
              QUAD_BM = SVL(I,j,K+2) - SVL(I,j,K)                       &
     &             - SVL_LAPSE * ( Z_FULL(I,j,K+2)-Z_HALF(I,j,K+2) )    &
     &             - SVL_LAPSE_BASE * ( Z_HALF(I,j,K+1)-Z_FULL(I,j,K) + &
     &                                                    DZL(I,j,K+1) )
              QUAD_C  = DZL(I,j,K+1)*( SVL(I,j,K+1) - SVL(I,j,K) -      &
     &             SVL_LAPSE_BASE * (                                   &
     &              Z_HALF(I,j,K+1)-Z_FULL(I,j,K) + 0.5*DZL(I,j,K+1) ) )

              IF ( QUAD_BM  >   0.0 ) THEN
                IF ( QUAD_C  <=  0.0) THEN
!                 ! SVL extrapolated from K to K+1 is greater than
!                 ! the level K+1 value - inversion needs to rise
                  DZ_DISC = DZ_DISC_MIN
                ELSEIF ( QUAD_BM*QUAD_BM  >=  4.0*QUAD_A*QUAD_C ) THEN
!                 ! solve equation for DZ_DISC...
                  IF ( QUAD_A  /=  0.0 ) THEN
!                   !   ...quadratic if QUAD_A NE 0
                    DZ_DISC = ( QUAD_BM - SQRT( QUAD_BM*QUAD_BM         &
     &                                       - 4.0*QUAD_A*QUAD_C )      &
     &                              ) / (2.0*QUAD_A)
                  ELSE
!                   !   ...linear if QUAD_A EQ 0
                    DZ_DISC = QUAD_C / QUAD_BM
                  END IF
                ELSE
                 DZ_DISC = 99999.9  ! large dummy value
                END IF

                IF ( DZ_DISC  >   0.9 * DZL(I,j,K+1) ) THEN
                 IF ( NTDSC(I,j) == 2 ) THEN
                  DZ_DISC = DZL(I,j,K+1)
                 ELSE
!               ! ZHSC diagnosed very close to or below Z_HALF(K+1):
                  IF ( SVL(I,j,K)-SVL(I,j,K-1)  >   0.0) THEN
!                   ! top of ML stably stratified so lower NTDSC but
!                   ! set ZHSC only fractionally (DZ_DISC_MIN)
!                   ! below the top of the inversion level.
                    NTDSC(I,j) = NTDSC(I,j) - 1
                    K=NTDSC(I,j)
                    TL_K = TL(I,j,K) + GRCP * Z_FULL(I,j,K)
                    TL_KP1 = TL(I,j,K+1) + GRCP * Z_FULL(I,j,K+1)
                    TL_KP2 = TL(I,j,K+2) + GRCP * Z_FULL(I,j,K+2)
                    DZ_DISC = DZ_DISC_MIN
                    DSCDEPTH(I,j) = DSCDEPTH(I,j) - DZL(I,j,K+1)
!                   ! Note that all but DZ_DISC_MIN of this layer will
!                   ! be added back on to DSCDEPTH a few lines below
                  ELSE
!                   ! top of ML well-mixed so don't lower the inversion
!                   ! level but set ZHSC just (DZ_DISC_MIN) above the
!                   ! half-level to allow the inversion to subside if
!                   ! necessary.
                    DZ_DISC = DZL(I,j,K+1) - DZ_DISC_MIN
                  END IF
                 END IF
                END IF

              ELSE  ! QUAD_BM le 0
!.. ignoring lapse rates
                DSVL_TOP = SVL(I,j,K+2) - SVL(I,j,K)
                DZ_DISC = DZL(I,j,K+1) *                                &
     &                          (SVL(I,j,K+1)-SVL(I,j,K)) / DSVL_TOP
              END IF

            ZHSC(I,j) = Z_HALF(I,j,K+2) - DZ_DISC
            DSCDEPTH(I,j) = DSCDEPTH(I,j) + ZHSC(I,j) - Z_HALF(I,j,K+1)
            DSC_DISC_INV(I,j) = 1  ! set flag to indicate disc inv found

!..Calculate discontinuous jumps of TL and QW:
!..only accurate enough for large enough DZ_DISC
            IF (DZ_DISC/DZL(I,j,K+1)  >   0.1) THEN

              DTL_DSC(I,j) = (TL_KP1-TL_K) * DZL(I,j,K+1) / DZ_DISC

              IF ( TL_KP2  >   TL_KP1 .AND.                             &
     &                           TL_KP1  >   TL_K ) THEN
                DTL_DSC(I,j) = MIN( TL_KP2-TL_K, DTL_DSC(I,j) )
              ELSEIF ( TL_KP2  <   TL_KP1 .AND.                         &
     &                           TL_KP1  <   TL_K ) THEN
                DTL_DSC(I,j) = MAX( TL_KP2-TL_K, DTL_DSC(I,j) )
              ELSE  ! TL non-monotonic
                DTL_DSC(I,j) = TL_KP1 - TL_K
              END IF

              DQW_DSC(I,j) = (QW(I,j,K+1)-QW(I,j,K))                    &
     &                                     * DZL(I,j,K+1) / DZ_DISC

              IF ( QW(I,j,K+2)  >   QW(I,j,K+1) .AND.                   &
     &               QW(I,j,K+1)  >   QW(I,j,K) ) THEN
                DQW_DSC(I,j) = MIN( QW(I,j,K+2)-QW(I,j,K), DQW_DSC(I,j))
              ELSEIF ( QW(I,j,K+2)  <   QW(I,j,K+1) .AND.               &
     &               QW(I,j,K+1)  <   QW(I,j,K) ) THEN
                DQW_DSC(I,j) = MAX( QW(I,j,K+2)-QW(I,j,K), DQW_DSC(I,j))
              ELSE  ! QW non-monotonic
                DQW_DSC(I,j) = QW(I,j,K+1)-QW(I,j,K)
              END IF

            ELSE
!..DZ_DISC small, ZHSC close to level above: use double jumps
              DTL_DSC(I,j) = TL_KP2 - TL_K
              DQW_DSC(I,j) = QW(I,j,K+2) - QW(I,j,K)
            END IF

          END IF  ! SVL increasing
          END IF ! test on K LT BL_LEVELS-2
         END IF ! test on DSC
      END DO
      END DO
!-----------------------------------------------------------------------
         do j=1,rows
         do i=1,row_length
!
!-----------------------------------------------------------------------
!! 0.5 Calculate the within-layer vertical gradients of cloud liquid
!!     and frozen water for the layer 1
!-----------------------------------------------------------------------
!
        VIRT_FACTOR = 1.0 + C_VIRTUAL*Q(I,j,1) - QCL(I,j,1) - QCF(I,j,1)
!
!       ! Here this is an estimate of the gradient adjustment applied
!       ! the previous timestep (assumes T1_SD has not changed much,
!       ! which in turn assumes RHOKH(1) has not)
        GRAD_T_ADJ(I,j) = MIN( MAX_T_GRAD ,                             &
     &                       A_GRAD_ADJ * T1_SD(I,j) / ZH_PREV(I,j) )
!        IF (T1_SD(I,j)  >   0.0) THEN
!          GRAD_Q_ADJ(I,j) = (Q1_SD(I,j) / T1_SD(I,j)) * GRAD_T_ADJ(I,j)
!        ELSE
          GRAD_Q_ADJ(I,j) = 0.0
!        END IF
        DTLDZ = -GRCP + GRAD_T_ADJ(I,j)
        DQCLDZ(I,j,1) = -( DTLDZ*DQSDT(I,j,1) +                         &
     &                   G*QS(I,j,1)/(R*T(I,j,1)*VIRT_FACTOR) )         &
     &                  / (1.0 + LCRCP*DQSDT(I,j,1))
        DQCFDZ(I,j,1) = -( DTLDZ*DQSDT(I,j,1) +                         &
     &                   G*QS(I,j,1)/(R*T(I,j,1)*VIRT_FACTOR) ) * FGF   &
     &                  / (1.0 + LSRCP*DQSDT(I,j,1))
!
!-----------------------------------------------------------------------
!! 0.6 Calculate the cloud liquid and frozen water contents at the
!!     top and bottom of layer 1
!-----------------------------------------------------------------------
!
!MHM limit calculation to greater than a small cloud fraction
        IF ( QCL(I,j,1) + QCF(I,j,1)  >   0.0                           &
     &       .and. CF(i,j,1)  >   1.e-3 ) THEN
          CFL(I,j,1) = CF(I,j,1) * QCL(I,j,1)/(QCL(I,j,1)+QCF(I,j,1))
          CFF(I,j,1) = CF(I,j,1) * QCF(I,j,1)/(QCL(I,j,1)+QCF(I,j,1))
        ELSE
          CFL(I,j,1) = 0.0
          CFF(I,j,1) = 0.0
        END IF
!
        IF (CFL(I,j,1)  >   0.0) THEN
          QCL_IC_TOP(I,j,1) = QCL(I,j,1) / CFL(I,j,1) +                 &
     &                       0.5*DZL(I,j,1)*DQCLDZ(I,j,1)
        ELSE
          QCL_IC_TOP(I,j,1) = 0.0
        END IF
!
        IF (CFF(I,j,1)  >   0.0) THEN
          QCF_IC_TOP(I,j,1) = QCF(I,j,1) / CFF(I,j,1) +                 &
     &                       0.5*DZL(I,j,1)*DQCFDZ(I,j,1)
        ELSE
          QCF_IC_TOP(I,j,1) = 0.0
        END IF
!
        QCL_IC_BOT(I,j,1) = 0.0
        QCF_IC_BOT(I,j,1) = 0.0
!
      END DO
      END DO
!
!-----------------------------------------------------------------------
!! 1.  First loop round boundary layer levels.
!-----------------------------------------------------------------------
!
      DO K=2,BL_LEVELS
!
         do j=1,rows
         do i=1,row_length
!
!-----------------------------------------------------------------------
!! 1.4 Calculate the within-layer vertical gradients of cloud liquid
!!     and frozen water for the current layer
!-----------------------------------------------------------------------
!
          IF (K  <=  NTML_PREV(I,j)) THEN
            DTLDZ = -GRCP + GRAD_T_ADJ(I,j)
          ELSE
            DTLDZ = -GRCP
          END IF
!
! RNBS correction  28/11/97
          VIRT_FACTOR = 1.0 + C_VIRTUAL*Q(I,j,k) - QCL(I,j,k) -         &
     &                        QCF(I,j,k)

          DQCLDZ(I,j,K) = -( DTLDZ*DQSDT(I,j,K)                         &
     &                   + G*QS(I,j,K)/(R*T(I,j,K)*VIRT_FACTOR) )       &
     &                    / ( 1.0 + LCRCP*DQSDT(I,j,K) )
          DQCFDZ(I,j,K) = -( DTLDZ*DQSDT(I,j,K)                         &
     &                   + G*QS(I,j,K)/(R*T(I,j,K)*VIRT_FACTOR) ) * FGF &
     &                    / ( 1.0 + LSRCP*DQSDT(I,j,K) )
!
!-----------------------------------------------------------------------
!! 1.5 Calculate the cloud liquid and frozen water contents at the
!!     top and bottom of the current layer
!-----------------------------------------------------------------------
!
!MHM limit calculation to greater than a small cloud fraction
        IF ( QCL(I,j,k) + QCF(I,j,k)  >   0.0                           &
     &       .and. CF(i,j,k)  >   1.e-3 ) THEN
       CFL(I,j,K) = CF(I,j,K) * QCL(I,j,K)/( QCL(I,j,K) + QCF(I,j,K) )
       CFF(I,j,K) = CF(I,j,K) * QCF(I,j,K)/( QCL(I,j,K) + QCF(I,j,K) )
          ELSE
            CFL(I,j,K) = 0.0
            CFF(I,j,K) = 0.0
          END IF
!
          IF (CFL(I,j,K)  >   0.0) THEN
            QCL_IC_TOP(I,j,K) = QCL(I,j,K) / CFL(I,j,K) +               &
     &                         0.5*DZL(I,j,K)*DQCLDZ(I,j,K)
            QCL_IC_BOT(I,j,K) = MAX( 0.0 , QCL(I,j,K) / CFL(I,j,K) -    &
     &                                    0.5*DZL(I,j,K)*DQCLDZ(I,j,K) )
          ELSE
            QCL_IC_TOP(I,j,K) = 0.0
            QCL_IC_BOT(I,j,K) = 0.0
          END IF
!
          IF (CFF(I,j,K)  >   0.0) THEN
            QCF_IC_TOP(I,j,K) = QCF(I,j,K) / CFF(I,j,K) +               &
     &                         0.5*DZL(I,j,K)*DQCFDZ(I,j,K)
            QCF_IC_BOT(I,j,K) = MAX( 0.0 , QCF(I,j,K) / CFF(I,j,K) -    &
     &                                    0.5*DZL(I,j,K)*DQCFDZ(I,j,K) )
          ELSE
            QCF_IC_TOP(I,j,K) = 0.0
            QCF_IC_BOT(I,j,K) = 0.0
          END IF
        END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!! 2.  Second loop round boundary layer levels.
!-----------------------------------------------------------------------
!
      DO K=2,BL_LEVELS
        KM1 = K-1
         do j=1,rows
         do i=1,row_length
!
!-----------------------------------------------------------------------
!! 2.1 Calculate the jumps of QW and TL across the layer interface
!!     at level k-1/2.
!-----------------------------------------------------------------------
!
          DQW = QW(I,j,K) - QW(I,j,KM1)              ! Used in P243.C2
!-----------------------------------------------------------------------
!         ! TL_K_BOT   = TL(I,j,K)   - 0.5*DZL(I,j,K)  *(-GRCP)
!         ! TL_KM1_TOP = TL(I,j,KM1) + 0.5*DZL(I,j,KM1)*(-GRCP)
!         ! DTL = TL_K_BOT - TL_KM1_TOP   so therefore
!-----------------------------------------------------------------------
          DTL = TL(I,j,K) - TL(I,j,KM1) + GRCP/RDZ(I,j,K)
                                ! Used in P243.C2
!
          DQCL = CFL(I,j,K)*QCL_IC_BOT(I,j,K) -                         &
     &             CFL(I,j,KM1)*QCL_IC_TOP(I,j,KM1)
          DQCF = CFF(I,j,K)*QCF_IC_BOT(I,j,K) -                         &
     &             CFF(I,j,KM1)*QCF_IC_TOP(I,j,KM1)
!
!-----------------------------------------------------------------------
!! 2.3 Calculate the buoyancy jumps across the interface between layers
!!     k and k-1
!-----------------------------------------------------------------------
!
          DB(I,j,K) = G * ( BTM(I,j,KM1)*DTL + BQM(I,j,KM1)*DQW +       &
     &               (LCRCP*BTM(I,j,KM1) - ETAR*BQM(I,j,KM1)) * DQCL +  &
     &                (LSRCP*BTM(I,j,KM1) - ETAR*BQM(I,j,KM1)) * DQCF )
!
          DB_CLD(I,j,K) =G * ( BTM_CLD(I,j,KM1)*DTL                     &
     &         + BQM_CLD(I,j,KM1)*DQW )
!
!ajm   + changed to - in chi_s calculation at v2p9
          CHI_S(I,j,K) = -QCL_IC_TOP(I,j,KM1) /                         &
     &                  (A_QSM(I,j,KM1)*DQW - A_DQSDTM(I,j,KM1)*DTL)
        END DO
        END DO
      END DO
!
!----------------------------------------------------------------------
!..For levels NTML+1 and NTDSC+1, replace grid-level jumps calculated
!..above with jumps calculated assuming a discontinuous subgrid
!..inversion structure.
!----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length

          IF ( SML_DISC_INV(I,j)  ==  1 )THEN
            K = NTML(I,j)+1
!
!..here _top means top of mixed layer, _bot means bottom of free atmos
            IF (CFL(I,j,K-1)  >   0.0) THEN
              QCL_IC_TOP(I,j,K) = QCL(I,j,K-1)/CFL(I,j,K-1) +           &
     &                   ( ZH(I,j)-Z_FULL(I,j,K-1) )*DQCLDZ(I,j,K-1)
            ELSE
              QCL_IC_TOP(I,j,K) = 0.0
            END IF

            IF (CFL(I,j,K+1)  >   0.0) THEN
              QCL_IC_BOT(I,j,K) = MAX( 0.0 , QCL(I,j,K+1)/CFL(I,j,K+1)  &
     &            - ( Z_FULL(I,j,K+1)-ZH(I,j) )*DQCLDZ(I,j,K+1) )
            ELSE
              QCL_IC_BOT(I,j,K) = 0.0
            END IF

            IF (CFF(I,j,K-1)  >   0.0) THEN
              QCF_IC_TOP(I,j,K) = QCF(I,j,K-1)/CFF(I,j,K-1) +           &
     &                  ( ZH(I,j)-Z_FULL(I,j,K-1) )*DQCFDZ(I,j,K-1)
            ELSE
              QCF_IC_TOP(I,j,K) = 0.0
            END IF

            IF (CFF(I,j,K+1)  >   0.0) THEN
              QCF_IC_BOT(I,j,K) = MAX( 0.0 , QCF(I,j,K+1)/CFF(I,j,K+1)  &
     &            - ( Z_FULL(I,j,K+1)-ZH(I,j) )*DQCFDZ(I,j,K+1) )
            ELSE
              QCF_IC_BOT(I,j,K) = 0.0
            END IF

            DQCL = CFL(I,j,K+1)*QCL_IC_BOT(I,j,K) -                     &
     &             CFL(I,j,K-1)*QCL_IC_TOP(I,j,K)
            DQCF = CFF(I,j,K+1)*QCF_IC_BOT(I,j,K) -                     &
     &             CFF(I,j,K-1)*QCF_IC_TOP(I,j,K)

            DB_DISC = G * ( BTM(I,j,K-1)*DTL_SML(I,j)                   &
     &                      + BQM(I,j,K-1)*DQW_SML(I,j) +               &
     &             (LCRCP*BTM(I,j,K-1) - ETAR*BQM(I,j,K-1)) * DQCL +    &
     &             (LSRCP*BTM(I,j,K-1) - ETAR*BQM(I,j,K-1)) * DQCF   )

            IF ( DB_DISC  <   0.0 ) THEN
!             ! Diagnosed inversion statically unstable:
!             ! use entrainment K (rather than fluxes) but also
!             ! ensure DB>0 so that entrainment is non-zero
              SML_DISC_INV(I,j) = 0
              ZH(I,j) = Z_HALF(I,j,NTML(I,j)+1)
            ELSE
              DB(I,j,K) = DB_DISC

              DB_CLD(I,j,K) = G * ( BTM_CLD(I,j,K-1)*DTL_SML(I,j) +     &
     &                            BQM_CLD(I,j,K-1)*DQW_SML(I,j) )

              CHI_S(I,j,K) = -QCL_IC_TOP(I,j,K) /                       &
     &            ( A_QSM(I,j,K-1)*DQW_SML(I,j)                         &
     &                  - A_DQSDTM(I,j,K-1)*DTL_SML(I,j) )
            END IF

          END IF  ! disc inversion diagnosed

          IF ( DSC_DISC_INV(I,j)  ==  1 ) THEN
            K = NTDSC(I,j)+1

!..here _top means top of mixed layer, _bot means bottom of free atmos
            IF (CFL(I,j,K-1)  >   0.0) THEN
              QCL_IC_TOP(I,j,K) = QCL(I,j,K-1)/CFL(I,j,K-1) +           &
     &            ( ZHSC(I,j)-Z_FULL(I,j,K-1) )*DQCLDZ(I,j,K-1)
            ELSE
              QCL_IC_TOP(I,j,K) = 0.0
            END IF

            IF (CFL(I,j,K+1)  >   0.0) THEN
              QCL_IC_BOT(I,j,K) = MAX( 0.0 , QCL(I,j,K+1)/CFL(I,j,K+1)  &
     &             - ( Z_FULL(I,j,K+1)-ZHSC(I,j) )*DQCLDZ(I,j,K+1) )
            ELSE
              QCL_IC_BOT(I,j,K) = 0.0
            END IF

            IF (CFF(I,j,K-1)  >   0.0) THEN
              QCF_IC_TOP(I,j,K) = QCF(I,j,K-1)/CFF(I,j,K-1) +           &
     &            ( ZHSC(I,j)-Z_FULL(I,j,K-1) )*DQCFDZ(I,j,K-1)
            ELSE
              QCF_IC_TOP(I,j,K) = 0.0
            END IF

            IF (CFF(I,j,K+1)  >   0.0) THEN
              QCF_IC_BOT(I,j,K) = MAX( 0.0 , QCF(I,j,K+1)/CFF(I,j,K+1)  &
     &             - ( Z_FULL(I,j,K+1)-ZHSC(I,j) )*DQCFDZ(I,j,K+1) )
            ELSE
              QCF_IC_BOT(I,j,K) = 0.0
            END IF

            DQCL = CFL(I,j,K+1)*QCL_IC_BOT(I,j,K) -                     &
     &             CFL(I,j,K-1)*QCL_IC_TOP(I,j,K)
            DQCF = CFF(I,j,K+1)*QCF_IC_BOT(I,j,K) -                     &
     &             CFF(I,j,K-1)*QCF_IC_TOP(I,j,K)

            DB_DISC = G * (                                             &
     &         BTM(I,j,K-1)*DTL_DSC(I,j) + BQM(I,j,K-1)*DQW_DSC(I,j) +  &
     &               (LCRCP*BTM(I,j,K-1) - ETAR*BQM(I,j,K-1)) * DQCL +  &
     &               (LSRCP*BTM(I,j,K-1) - ETAR*BQM(I,j,K-1)) * DQCF   )

            IF ( DB_DISC  <   0.0 ) THEN
!             ! Diagnosed inversion statically unstable:
!             ! use entrainment K (rather than fluxes) but also
!             ! ensure DB>0 so that entrainment is non-zero
              DSC_DISC_INV(I,j) = 0
              ZHSC(I,j) = Z_HALF(I,j,NTDSC(I,j)+1)
            ELSE
              DB(I,j,K) = DB_DISC

              DB_CLD(I,j,K) = G * ( BTM_CLD(I,j,K-1)*DTL_DSC(I,j)       &
     &                          + BQM_CLD(I,j,K-1)*DQW_DSC(I,j) )

              CHI_S(I,j,K) = -QCL_IC_TOP(I,j,K) /                       &
     &                 ( A_QSM(I,j,K-1)*DQW_DSC(I,j)                    &
     &                         - A_DQSDTM(I,j,K-1)*DTL_DSC(I,j) )
            END IF

          END IF  ! disc inversion diagnosed

      END DO
      END DO
!-----------------------------------------------------------------------
!L 3.0a Calculate surface buoyancy flux
!-----------------------------------------------------------------------
!
      do j=1,rows
      do i=1,row_length

! use mixed-layer average of buoyancy parameters
        BFLUX_SURF(I,j) = 0.5 * G * (                                   &
     &       (BTM(I,j,1)+BTM(I,j,NTML(I,j)))*FTL(I,j,1) +               &
     &       (BQM(I,j,1)+BQM(I,j,NTML(I,j)))*FQW(I,j,1) )

        IF ( BFLUX_SURF(I,j)   >   0.0 ) THEN
          BFLUX_SURF_SAT(I,j) = 0.5 * G * (                             &
     &       (BTM_CLD(I,j,1)+BTM_CLD(I,j,NTML(I,j)))*FTL(I,j,1) +       &
     &       (BQM_CLD(I,j,1)+BQM_CLD(I,j,NTML(I,j)))*FQW(I,j,1) )
          IF ( COUPLED(I,j) ) BFLUX_SURF_SAT(I,j) = 0.5 * G * (         &
     &       (BTM_CLD(I,j,1)+BTM_CLD(I,j,NTDSC(I,j)))*FTL(I,j,1) +      &
     &       (BQM_CLD(I,j,1)+BQM_CLD(I,j,NTDSC(I,j)))*FQW(I,j,1) )
        ELSE
          BFLUX_SURF_SAT(I,j) = 0.0
        END IF

      END DO
      END DO
!-----------------------------------------------------------------------
!! 3.0aa Calculate uniform mixed-layer cloud fraction and thence
!!       estimate Sc layer cloud depth (not cloud fraction weighted).
!!       (If DSC=.FALSE. then NTDSC=0 and ZC_DSC remains equal to 0.)
!-----------------------------------------------------------------------
! First the SML
!--------------
! First find cloud-base as seen by the cloud scheme, at grid-level
! K_LEVEL and height Z_CF_BASE, to use as first guess or lower limit
!-----------------------------------------------------------------------
      DO J=1,ROWS
      DO I=1,ROW_LENGTH
         K_LEVEL(I,J) = NTML(I,J)
         DO WHILE ( CF(I,J,K_LEVEL(I,J))  >   SC_CFTOL                  &
     &               .AND. K_LEVEL(I,J)  >=  2 )
           K_LEVEL(I,J) = K_LEVEL(I,J) - 1
         END DO
      END DO
      END DO
      do j=1,rows
      do i=1,row_length
        CLOUD_BASE(I,j)= .FALSE.
        ZC(I,j)        = 0.0
        K_CBASE(I,j)   = 0
        Z_CF_BASE(I,j) = ZH(I,j)
        Z_CTOP(I,j)    = ZH(I,j)
!       ! Use a single CF for whole mixed-layer (more realistic).
!       ! Include NTML+1 if a subgrid inversion has been diagnosed
        IF ( COUPLED(I,j) .OR. CUMULUS(I,j) ) THEN
          CF_SML(I,j)=0.0
        ELSE
          K = NTML(I,j)
          CF_SML(I,j) = MAX( CF(I,j,K), CF(I,j,K+SML_DISC_INV(I,j)) )
        END IF
      END DO
      END DO
!
! initialisation of z_cf_base and first guess for cloud depth, ZC
!
      do j=1,rows
      do i=1,row_length
        IF ( CF_SML(I,j)  >   SC_CFTOL ) THEN
          K = NTML(I,j)
          IF ( SML_DISC_INV(I,j)  ==  0 .AND.                           &
     &         CF(I,j,K+1)  >   SC_CFTOL) THEN
!           ! if no subgrid inversion and level NTML+1 is cloudy
            Z_CTOP(I,j) = Z_TOP(I,j,K+1)
          END IF
          IF ( K_LEVEL(I,J)  ==  1 .AND.                                &
     &         CF(I,J,K_LEVEL(I,J))  >   SC_CFTOL) THEN
            Z_CF_BASE(I,J) = 0.0
          ELSE
            Z_CF_BASE(I,J) = Z_HALF(I,J,K_LEVEL(I,J)+1)
          END IF
          ZC(I,j) = Z_CTOP(I,j) - Z_CF_BASE(I,j)
        END IF
      END DO
      END DO
!--------------------------------------------------
! Find lowest level within ML with max CF
!--------------------------------------------------
      DO K = BL_LEVELS,1,-1
        do j=1,rows
        do i=1,row_length
          IF ( .NOT.CLOUD_BASE(I,j) .AND. K  <=  NTML(I,j)+1 .AND.      &
     &         CF_SML(I,j)  >   SC_CFTOL ) THEN
!              ! within cloudy boundary layer
            IF ( K  ==  1) THEN
              CLOUD_BASE(I,j) = .TRUE.
            ELSE
             IF ( CF(I,j,K-1)  <   CF(I,j,K) ) CLOUD_BASE(I,j) = .TRUE.
            END IF
            K_CBASE(I,j) = K
          END IF
        end do
        end do
      end do

!Initialise K_CFF = lowest level with ice cloud
      DO J=1,ROWS
      DO I=1,ROW_LENGTH
        K_CFF(I,J) = K_CBASE(I,J)
        IF (K_CFF(I,J) > 1) THEN
          DO WHILE ( CFF(I,J,K_CFF(I,J))  >   SC_CFTOL                  &
     &             .AND. K_CFF(I,J)  >   1)
             K_CFF(I,J) = K_CFF(I,J) - 1
          END DO
        END IF
      END DO
      END DO

      do j=1,rows
      do i=1,row_length

!         !--------------------------------------------------
!         ! Use adiabatic qcl gradient to estimate cloud-base
!         ! from in-cloud qcl in level K_CBASE
!         ! If k_cbase = 0 then it hasn't been initialised
!         !--------------------------------------------------

          IF ( CLOUD_BASE(I,J) .AND. K_CBASE(I,J) /= 0 ) THEN
            Z_CBASE = Z_FULL(I,j,K_CBASE(I,j)) -                        &
     &                QCL(I,j,K_CBASE(I,j)) /                           &
     &                ( CF(I,j,K_CBASE(I,j))*DQCLDZ(I,j,K_CBASE(I,j)) )
            IF ( DQCFDZ(I,j,K_CBASE(I,j))  >   0.0 ) THEN
              Z_CBASE = MIN( Z_CBASE, Z_FULL(I,j,K_CBASE(I,j)) -        &
     &               QCF(I,j,K_CBASE(I,j)) /                            &
     &               ( CF(I,j,K_CBASE(I,j))*DQCFDZ(I,j,K_CBASE(I,j)) )  &
     &                     )
            ELSE
!             !---------------------------------------------------------
!             ! No adiabatic QCF gradient so find lowest level, K_CFF,
!             ! with CFF>SC_CFTOL and assume cloud-base within that leve
!             !---------------------------------------------------------
              IF ( CFF(I,j,K_CFF(i,j))  <=  SC_CFTOL .AND.              &
     &                      K_CFF(i,j)  <   K_CBASE(I,j) )              &
     &             K_CFF(i,j) = K_CFF(i,j) + 1
!                  ! will want to raise K_CFF back up one level unless
!                  ! level 1 is cloudy or no sig frozen cloud at all
              Z_CBASE = MIN( Z_CBASE, Z_TOP(I,j,K_CFF(i,j)) -           &
     &                  DZL(I,j,K_CFF(i,j))                             &
     &                * CFF(I,j,K_CFF(i,j))/CF(I,j,K_CFF(i,j)) )
            END IF
!           !------------------------------------------------------
!           ! use cloud-base as seen by cloud scheme as lower limit
!           ! and base of level NTML+1 as upper limit
!           !------------------------------------------------------
            Z_CBASE = MIN( Z_HALF(I,j,NTML(I,j)+1),                     &
     &                     MAX( Z_CF_BASE(I,j), Z_CBASE) )

            ZC(I,j) = Z_CTOP(I,j) - Z_CBASE
          END IF

      END DO
      END DO
!-----------------------------------------------------------------------
! Second DSC layer
!-----------------------------------------------------------------------
      DO J=1,ROWS
      DO I=1,ROW_LENGTH
        K_LEVEL(I,J) = NTDSC(I,J)
        IF ( K_LEVEL(I,J) >= 2) THEN
          DO WHILE ( CF(I,J,K_LEVEL(I,J))  >   SC_CFTOL                 &
     &             .AND. K_LEVEL(I,J)  >=  2 )
              K_LEVEL(I,J) = K_LEVEL(I,J) - 1
          END DO
        END IF
      END DO
      END DO

      DO J=1,ROWS
      DO I=1,ROW_LENGTH
        CLOUD_BASE(I,j) = .FALSE.
        ZC_DSC(I,j) = 0.0
        CF_DSC(I,j) = 0.0
        K_CBASE(I,j) = 0
        Z_CF_BASE(I,j) = ZHSC(I,j)
        Z_CTOP(I,j)    = ZHSC(I,j)

        IF ( DSC(I,j) ) THEN
          K = NTDSC(I,j)
          CF_DSC(I,j) = MAX( CF(I,j,K), CF(I,j,K+1) )
!         !-------------------------------------------------------------
!         ! Find cloud-base as seen by cloud scheme, Z_CF_BASE,
!         ! to use as first guess or lower limit and find cloud top.
!         !-------------------------------------------------------------
          IF ( CF_DSC(I,j)  >   SC_CFTOL ) THEN
            K = NTDSC(I,j)
            IF ( DSC_DISC_INV(I,j)  ==  0 .AND.                         &
     &           CF(I,j,K+1)  >   SC_CFTOL) THEN
!             ! if no subgrid inversion and level NTML+1 is cloudy
!             ! then include this layer in cloud-depth
              Z_CTOP(I,j) = Z_TOP(I,j,K+1)
            END IF
            IF ( K_LEVEL(i,j)  ==  1 .AND.                              &
     &            CF(I,j,K_LEVEL(i,j))  >   SC_CFTOL) THEN
              Z_CF_BASE(I,j) = 0.0
            ELSE
              Z_CF_BASE(I,j) = Z_HALF(I,j,K_LEVEL(i,j)+1)
            END IF
            ZC_DSC(I,j) = Z_CTOP(I,j) - Z_CF_BASE(I,j)   ! first guess
          END IF
        END IF
      END DO
      END DO
!--------------------------------------------------
! Find lowest level within ML with max CF
!--------------------------------------------------
      DO K = BL_LEVELS,1,-1
        do j=1,rows
        do i=1,row_length
          IF ( .NOT.CLOUD_BASE(I,j) .AND. K  <=  NTDSC(I,j)+1 .AND.     &
     &         CF_DSC(I,j)  >   SC_CFTOL ) THEN
!              ! within cloudy boundary layer
            IF ( K  ==  1) THEN
              CLOUD_BASE(I,j) = .TRUE.
            ELSE
              IF ( CF(I,j,K-1)  <   CF(I,j,K) ) CLOUD_BASE(I,j) = .TRUE.
            END IF
            K_CBASE(I,j) = K
          END IF
        END DO ! I
        END DO ! J
      END DO ! K

! Initialise K_CFF
      DO J=1,ROWS
      DO I=1,ROW_LENGTH
        K_CFF(I,J) = K_CBASE(I,J)
        IF (K_CFF(I,J) > 1) THEN
          DO WHILE ( CFF(I,J,K_CFF(I,J))  >   SC_CFTOL                  &
     &                  .AND. K_CFF(I,J)  >   1)
            K_CFF(I,J) = K_CFF(I,J) - 1
          END DO
        END IF
      END DO
      END DO

      DO J=1,ROWS
      DO I=1,ROW_LENGTH

!       !--------------------------------------------------
!       ! use adiabatic qcl gradient to estimate cloud-base
!       ! from in-cloud qcl in level K_CBASE
!       !--------------------------------------------------
        IF ( CLOUD_BASE(I,j) .AND. K_CBASE(i,j)  /= 0 ) THEN
          Z_CBASE = Z_FULL(I,j,K_CBASE(I,j)) -                          &
     &              QCL(I,j,K_CBASE(I,j)) /                             &
     &              ( CF(I,j,K_CBASE(I,j))*DQCLDZ(I,j,K_CBASE(I,j)) )
          IF ( DQCFDZ(I,j,K_CBASE(I,j))  >   0.0 ) THEN
            Z_CBASE = MIN( Z_CBASE, Z_FULL(I,j,K_CBASE(I,j)) -          &
     &            QCF(I,j,K_CBASE(I,j)) /                               &
     &            ( CF(I,j,K_CBASE(I,j))*DQCFDZ(I,j,K_CBASE(I,j)) )     &
     &                   )
          ELSE
!           !----------------------------------------------------------
!           ! No adiabatic QCF gradient so find lowest level, K_CFF,
!           ! with CFF>SC_CFTOL and assume cloud-base within that level
!           !----------------------------------------------------------
            IF ( CFF(I,j,K_CFF(i,j))  <=  SC_CFTOL .AND.                &
     &                    K_CFF(i,j)  <   K_CBASE(I,j) )                &
     &           K_CFF(i,j) = K_CFF(i,j) + 1
!                ! will want to raise K_CFF back up one level unless
!                ! level 1 is cloudy or no sig frozen cloud at all
            Z_CBASE = MIN( Z_CBASE, Z_TOP(I,j,K_CFF(i,j)) -             &
     &                DZL(I,j,K_CFF(i,j))                               &
     &               * CFF(I,j,K_CFF(i,j))/CF(I,j,K_CFF(i,j)) )
          END IF
!         !------------------------------------------------------
!         ! use cloud-base as seen by cloud scheme as lower limit
!         ! and base of level NTDSC+1 as upper limit
!         !------------------------------------------------------
          Z_CBASE = MIN( Z_HALF(I,j,NTDSC(I,j)+1),                      &
     &                   MAX( Z_CF_BASE(I,j) , Z_CBASE) )

          ZC_DSC(I,j) = Z_CTOP(I,j) - Z_CBASE
        END IF

      END DO !I
      END DO !J
!     !-----------------------------------------------------------------
!     !  Layer cloud depth cannot be > the layer depth itself.
!     !-----------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
        ZC_DSC(I,j) = MIN( ZC_DSC(I,j), DSCDEPTH(I,j) )
      END DO
      END DO
!-----------------------------------------------------------------------
!- Calculate buoyancy flux factor used in the diagnosis of decoupling:
!-----------------------------------------------------------------------
      DO K = 2,BL_LEVELS
        do j=1,rows
        do i=1,row_length
          KL = K

          Z_INT_TOP = Z_FULL(I,j,K)
          IF (.NOT.CUMULUS(I,j) .AND. .NOT.COUPLED(I,j) .AND.           &
     &        K  ==  NTML(I,j)+1 .AND. NTML(I,j)  >=  2) THEN
            KL = NTML(I,j)
            Z_INT_TOP = ZH(I,j)
          ELSEIF (K  ==  NTDSC(I,j)+1) THEN
            KL = NTDSC(I,j)
            Z_INT_TOP = ZHSC(I,j)
          END IF
          KM1 = KL-1
          Z_INT = Z_INT_TOP - Z_FULL(I,j,K-1)  ! integration depth

          DQW = QW(I,j,KL) - QW(I,j,KM1)
          DTL = TL(I,j,KL) - TL(I,j,KM1) + GRCP/RDZ(I,j,KL)
          IF (KL  <=  NTML_PREV(I,j)) THEN
            DTL_GA = DTL - GRAD_T_ADJ(I,j)/RDZ(I,j,KL)
            DQW_GA = DQW - GRAD_Q_ADJ(I,j)/RDZ(I,j,KL)
          ELSE
            DTL_GA = DTL
            DQW_GA = DQW
          END IF
!         !----------------------------------------------------------
!         ! CF_FOR_WB is uniform `bl' CF for use within cloud layers
!         !----------------------------------------------------------
          CF_FOR_WB = 0.0
          Z_CBASE = ZH(I,j)-ZC(I,j)
          ZDSC_CBASE = ZHSC(I,j)-ZC_DSC(I,j)
          IF ( K  <=  NTML(I,j)+1 .AND.                                 &
     &         Z_FULL(I,j,K)  >=  Z_CBASE) CF_FOR_WB = CF_SML(I,j)
          IF ( K  <=  NTDSC(I,j)+1 .AND.                                &
     &         Z_FULL(I,j,K)  >=  ZDSC_CBASE) CF_FOR_WB = CF_DSC(I,j)
!         !----------------------------------------------------------
!         ! WB = -K_SURF*(DB/DZ - gamma_buoy) - K_TOP*DB/DZ
!         ! This is integrated in EXCF_NL, iterating the K profiles.
!         ! Here the relevant integrated DB/DZ factors are calculated
!         !----------------------------------------------------------
          DB_KSURF_DRY(I,j,K) = - G*RDZ(I,j,KL)*Z_INT*                  &
     &             ( BTM(I,j,KM1)*DTL_GA + BQM(I,j,KM1)*DQW_GA )
          DB_KTOP_DRY(I,j,K)  = - G*RDZ(I,j,KL)*Z_INT*                  &
     &             ( BTM(I,j,KM1)*DTL + BQM(I,j,KM1)*DQW )
          DB_KSURF_CLD(I,j,K) = - G*RDZ(I,j,KL)*Z_INT*                  &
     &             ( BTM_CLD(I,j,KM1)*DTL_GA + BQM_CLD(I,j,KM1)*DQW_GA )
          DB_KTOP_CLD(I,j,K)  = - G*RDZ(I,j,KL)*Z_INT*                  &
     &             ( BTM_CLD(I,j,KM1)*DTL + BQM_CLD(I,j,KM1)*DQW )
!         !-------------------------------------------------------
!         ! Weight cloud layer factors with cloud fraction
!         !-------------------------------------------------------
          DB_KSURF_CLD(I,j,K) = DB_KSURF_DRY(I,j,K)*(1.0-CF_FOR_WB) +   &
     &                          DB_KSURF_CLD(I,j,K)*CF_FOR_WB
          DB_KTOP_CLD(I,j,K)  = DB_KTOP_DRY(I,j,K)*(1.0-CF_FOR_WB) +    &
     &                          DB_KTOP_CLD(I,j,K)*CF_FOR_WB

        END DO
      END DO
      END DO
!-----------------------------------------------------------------------
!L 3.1 Calculate inputs for the top of b.l. entrainment parametrization
!-----------------------------------------------------------------------
!
      do j=1,rows
      do i=1,row_length
       ZETA_R_DSC(I,j) = 0.0
       CHI_S_DSCT(I,j) = 0.0
       CLD_FACTOR_DSC(I,j) = 0.0
       BT_DSCT(I,j) = 0.0
       BTT_DSCT(I,j) = 0.0
       BTC_DSCT(I,j) = 0.0
       DB_DSCT(I,j) = 0.0
       DB_DSCT_CLD(I,j) = 0.0
       CHI_S_TOP(I,j) = 0.0
       CLD_FACTOR(I,j) = 0.0
       BT_TOP(I,j) = 0.0
       BTT_TOP(I,j) = 0.0
       BTC_TOP(I,j) = 0.0
       DB_TOP(I,j) = 0.0
       DB_TOP_CLD(I,j) = 0.0    ! default required if COUPLED
       Z_CLD(I,j) = 0.0
       Z_CLD_DSC(I,j) = 0.0
      END DO
      END DO
!
      DO K = 1,BL_LEVELS
        do j=1,rows
        do i=1,row_length
          IF ( K  <=  NTML(I,j)+1 ) THEN
!           !---------------------------------------------------
!           ! Calculation of cloud fraction weighted
!           ! thickness of cloud in the surface mixed layer
!           !---------------------------------------------------
            Z_CLD(I,j) = Z_CLD(I,j) +                                   &
     &                 CF(I,j,K) * 0.5 * DZL(I,j,K) +                   &
     &                  MIN( CFL(I,j,K) * 0.5 * DZL(I,j,K) ,            &
     &                          QCL(I,j,K) / DQCLDZ(I,j,K) )
            IF ( DQCFDZ(I,j,K)  >   0.0) THEN
              Z_CLD(I,j) = Z_CLD(I,j) +                                 &
     &                  MIN( CFF(I,j,K) * 0.5 * DZL(I,j,K) ,            &
     &                          QCF(I,j,K) / DQCFDZ(I,j,K) )
            ELSE
              Z_CLD(I,j) = Z_CLD(I,j) + CFF(I,j,K) * 0.5 * DZL(I,j,K)
            END IF
          END IF
!
          IF ( DSC(I,j) .AND. K <= NTDSC(I,j)+1 .AND.                   &
     &         ( COUPLED(I,j) .OR.                                      &
     &               Z_TOP(I,j,K) >= ZHSC(I,j)-ZC_DSC(I,j) ) ) THEN
!           !-----------------------------------------------------------
!           ! Calculation of cloud fraction weighted thickness of
!           ! cloud in the DSC layer (or to the surface if COUPLED)
!           !-----------------------------------------------------------
            Z_CLD_DSC(I,j) = Z_CLD_DSC(I,j) +                           &
     &                 CF(I,j,K) * 0.5 * DZL(I,j,K) +                   &
     &                  MIN( CFL(I,j,K) * 0.5 * DZL(I,j,K) ,            &
     &                          QCL(I,j,K) / DQCLDZ(I,j,K) )
            IF ( DQCFDZ(I,j,K)  >   0.0) THEN
              Z_CLD_DSC(I,j) = Z_CLD_DSC(I,j) +                         &
     &                  MIN( CFF(I,j,K) * 0.5 * DZL(I,j,K) ,            &
     &                          QCF(I,j,K) / DQCFDZ(I,j,K) )
            ELSE
              Z_CLD_DSC(I,j) = Z_CLD_DSC(I,j) +                         &
     &                                CFF(I,j,K) * 0.5 * DZL(I,j,K)
            END IF
          END IF
!
          IF (K  ==  NTML(I,j)  .AND. .NOT.COUPLED(I,j) ) THEN
!           !------------------------------------------------------
!           ! Calculation of SML inputs.  If COUPLED then these are
!           ! not used (as no entrainment is then applied at ZH)
!           !------------------------------------------------------
            CHI_S_TOP(I,j) = MAX( 0.0, MIN( CHI_S(I,j,K+1), 1.) )
            KP2=MIN(K+1+SML_DISC_INV(I,j),BL_LEVELS)
            CLD_FACTOR(I,j) = MAX( 0.0 , CF(I,j,K)-CF(I,j,KP2) )
            BT_TOP(I,j) = G * BTM(I,j,K)
            BTT_TOP(I,j) = G * BTM_CLD(I,j,K)
            BTC_TOP(I,j) = BTT_TOP(I,j)
            DB_TOP(I,j) = DB(I,j,K+1)
            DB_TOP_CLD(I,j) = DB_CLD(I,j,K+1)
            IF ( DB_TOP(I,j)  <   0.001 ) THEN
!             ! Diagnosed inversion statically unstable:
!             ! ensure DB>0 so that entrainment is non-zero and
!             ! instability can be removed.
              DB_TOP(I,j) = 0.001
              DB_TOP_CLD(I,j) = 0.0  ! set buoyancy reversal
              CHI_S_TOP(I,j) = 0.0   ! term to zero
          END IF
          END IF
          IF (K  ==  NTDSC(I,j)) THEN
!           !---------------------------------------------------
!           ! Calculation of DSC inputs
!           ! (if DSC=.FALSE. then K never equals NTDSC(=0))
!           !---------------------------------------------------
            CHI_S_DSCT(I,j) = MAX( 0.0, MIN( CHI_S(I,j,K+1), 1.) )
            KP2=MIN(K+1+DSC_DISC_INV(I,j),BL_LEVELS)
            CLD_FACTOR_DSC(I,j) = MAX( 0.0 , CF(I,j,K)-CF(I,j,KP2) )
            BT_DSCT(I,j) = G * BTM(I,j,K)
            BTT_DSCT(I,j) = G * BTM_CLD(I,j,K)
            BTC_DSCT(I,j) = BTT_DSCT(I,j)
            DB_DSCT(I,j) = DB(I,j,K+1)
            DB_DSCT_CLD(I,j) = DB_CLD(I,j,K+1)
            IF ( DB_DSCT(I,j)  <   0.001 ) THEN
!             ! Diagnosed inversion statically unstable:
!             ! ensure DB>0 so that entrainment is non-zero and
!             ! instability can be removed.
              DB_DSCT(I,j) = 0.001
              DB_DSCT_CLD(I,j) = 0.0  ! set buoyancy reversal
              CHI_S_DSCT(I,j) = 0.0   ! term to zero
          END IF
          END IF
        END DO
        END DO
      END DO
!-----------------------------------------------------------------------
! Next those terms which depend on the presence of buoyancy reversal
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
        Z_CLD(I,j) = MIN( Z_CLD(I,j), ZH(I,j) )
        Z_CLD_DSC(I,j) = MIN( Z_CLD_DSC(I,j), ZHSC(I,j) )
!       !---------------------------------------------------------------
!       ! First the surface mixed layer.
!       !---------------------------------------------------------------
        IF ( COUPLED(I,j) ) THEN
         ZETA_S(I,j) = 1.0 - Z_CLD_DSC(I,j) / ZHSC(I,j)
         ZETA_R(I,j) = 1.0 - ZC_DSC(I,j) / ZHSC(I,j)
        ELSE
         ZETA_S(I,j) = 1.0 - Z_CLD(I,j) / ZH(I,j)
         ZETA_R(I,j) = 1.0 - ZC(I,j) / ZH(I,j)
        END IF
!
        IF (DB_TOP_CLD(I,j)  >=  0.0) THEN
!         !--------------------------------------------------
!         ! i.e. no buoyancy reversal (or default if COUPLED)
!         !--------------------------------------------------
          DB_TOP_CLD(I,j) = 0.0
          D_SIEMS(I,j) = 0.0
        ELSE
!         !----------------------------
!         ! IF (DB_TOP_CLD(I,j)  <   0.0)
!         ! i.e. buoyancy reversal
!         !----------------------------
          DB_TOP_CLD(I,j) = -DB_TOP_CLD(I,j) * CLD_FACTOR(I,j)
          D_SIEMS(I,j) = MAX( 0.0,                                      &
     &         CHI_S_TOP(I,j) * DB_TOP_CLD(I,j) / (DB_TOP(I,j)+1.E-14) )
          ZETA_R(I,j) = MIN( ZETA_R(I,j)+10.0*(1.0-ZETA_R(I,j))         &
     &                                    *D_SIEMS(I,j), 1.0 )
        END IF
!       !---------------------------------------------------------------
!       ! Now the decoupled Sc layer (DSC).
!       !---------------------------------------------------------------
        IF (DSC(I,j)) THEN
          IF ( COUPLED(I,j) ) THEN
            ZETA_R_DSC(I,j) = 1.0 - ZC_DSC(I,j) / ZHSC(I,j)
          ELSE
            ZETA_R_DSC(I,j) = 1.0 - ZC_DSC(I,j) / DSCDEPTH(I,j)
          END IF

          IF (DB_DSCT_CLD(I,j)  >=  0.0) THEN
!           !----------------------------
!           ! i.e. no buoyancy reversal
!           !----------------------------
            DB_DSCT_CLD(I,j) = 0.0
            D_SIEMS_DSC(I,j) = 0.0
          ELSE
!           !----------------------------
!           ! IF (DB_DSCT_CLD(I,j)  <   0.0)
!           ! i.e. buoyancy reversal
!           !----------------------------
            DB_DSCT_CLD(I,j) = -DB_DSCT_CLD(I,j) * CLD_FACTOR_DSC(I,j)
            D_SIEMS_DSC(I,j) = MAX( 0.0, CHI_S_DSCT(I,j)                &
     &                      * DB_DSCT_CLD(I,j) / (DB_DSCT(I,j)+1.E-14) )
            ZETA_R_DSC(I,j) = MIN( ZETA_R_DSC(I,j) +                    &
     &              10.0*(1.0-ZETA_R_DSC(I,j)) * D_SIEMS_DSC(I,j), 1.0 )
          END IF
        END IF
      END DO
      END DO
!
!-----------------------------------------------------------------------
! 4. Calculate the radiative flux change across cloud top for mixed-
!    layer to ZH.  Restrict search for maximum divergence to below
!    NTML+2.  This may introduce errors if NTML changes a lot during
!    the radiative timestep but can't be helped.
!-----------------------------------------------------------------------
!     Initialise variables
!------------------------------
      do j=1,rows
      do i=1,row_length
        K_CLOUD_TOP(I,j) = 0
        DF_TOP_OVER_CP(I,j) = 0.0
      END DO
      END DO
!
      DO K = 1,BL_LEVELS
        do j=1,rows
        do i=1,row_length
! *APL*restrict search to `close' to ZH
          K_RAD_SMLT = NTML(I,j)+1
!         !-------------------------------------------------------------
!         ! Find the layer with the LW greatest rad flux jump below
!         ! K_RAD_SMLT and assume that this is the top of the SML.
!         !-------------------------------------------------------------
          IF (DFLW_OVER_CP(I,j,K)  >   DF_TOP_OVER_CP(I,j)              &
     &                .AND. K  <=  K_RAD_SMLT ) THEN
            K_CLOUD_TOP(I,j) = K
            DF_TOP_OVER_CP(I,j) = DFLW_OVER_CP(I,j,K)
          END IF

        END DO
        END DO
      END DO
!     !-----------------------------------------------------------------
!     !  In case the radiative divergence is spread over 2 or 3 levels,
!     !  add on any cooling either side of the level of maximum cooling.
!     !-----------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
          IF ( K_CLOUD_TOP(I,j)  >   0 ) THEN
            K=K_CLOUD_TOP(I,j)
            DFSW_TOP = DFSW_OVER_CP(I,j,K)
            IF ( K >  1 ) THEN
              DF_TOP_OVER_CP(I,j) = DF_TOP_OVER_CP(I,j)                 &
     &                              + MAX(0.0, DFLW_OVER_CP(I,j,K-1) )
              DFSW_TOP = DFSW_TOP + MIN(0.0, DFSW_OVER_CP(I,j,K-1) )
            END IF
            IF ( K <  BL_LEVELS ) THEN
              DF_TOP_OVER_CP(I,j) = DF_TOP_OVER_CP(I,j)                 &
     &                              + MAX(0.0, DFLW_OVER_CP(I,j,K+1) )
              DFSW_TOP = DFSW_TOP + MIN(0.0, DFSW_OVER_CP(I,j,K+1) )
            END IF
!           !-----------------------------------------------------------
!           ! Combine SW and LW cloud-top divergences into a net
!           ! divergence by estimating SW flux divergence at LW
!           ! extinction depth = DF_SW * (1-exp{-A*kappa_sw/kappa_lw})
!           ! Choose A=3 (to achieve 95% extinction of LW) and
!           ! kappa_lw = 3*kappa_sw
!           !-----------------------------------------------------------
            DF_TOP_OVER_CP(I,j) = MAX( 0.0,                             &
     &                  DF_TOP_OVER_CP(I,j) + 0.35 * DFSW_TOP )
          END IF
      END DO
      END DO
!-----------------------------------------------------------------------
!! 5.1 Subroutine EXCF_NL.
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
        NTML_SAVE(I,j) = NTML(I,j)  ! needed to identify changes
        DSC_SAVE(I,j) = DSC(I,j)    !      in excf_nl
      END DO
      END DO
!
! DEPENDS ON: excf_nl
      CALL EXCF_NL (                                                    &
     & row_length,rows,BL_LEVELS,                                       &
     & NTML,CF,                                                         &
     & RDZ,ZH,Z_UV,Z_TQ,RHO_UV,RHO_TQ,RHOSTAR_GB,                       &

     & Z0M,V_S,FB_SURF,DB_TOP,                                          &
     & BFLUX_SURF,BFLUX_SURF_SAT,ZETA_S,BT_TOP,BTT_TOP,                 &
     & DF_TOP_OVER_CP,ZETA_R,BTC_TOP,                                   &
     & DB_TOP_CLD,CHI_S_TOP,ZC,                                         &
     & RHOKM(1,1,2),RHOKH(1,1,2),RHOKM_TOP(1,1,2),RHOKH_TOP(1,1,2),     &
     & DECFIX,STOPWE_SBL,                                               &
     & NG_STRESS,F_NGSTRESS,ZHSC,DSCDEPTH,NTDSC,DB_DSCT,SVL,            &
     & BT_DSCT,BTT_DSCT,                                                &
     & DF_DSCT_OVER_CP,ZETA_R_DSC,BTC_DSCT,                             &
     & DB_DSCT_CLD,CHI_S_DSCT,ZC_DSC,COUPLED,                           &
     & D_SIEMS,D_SIEMS_DSC,NBDSC,                                       &
     & DB_KSURF_DRY,DB_KTOP_DRY,DB_KSURF_CLD,DB_KTOP_CLD,               &
     & DSC,CUMULUS,ZDSC_BASE,                                           &
     & RHOKH_TOP_ENT,RHOKH_DSCT_ENT,RHOKH_SURF_ENT,                     &
     & LTIMER                                                           &
     &)

!-----------------------------------------------------------------------
!-adjust SML/DSC properties depending on diagnoses in EXCF_NL
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
        IF ( DSC(I,j) .AND. .NOT.DSC_SAVE(I,j) ) THEN
!..decoupling diagnosed in EXCF_NL - change parameters around
          DTL_DSC(I,j) = DTL_SML(I,j)
          DQW_DSC(I,j) = DQW_SML(I,j)
          DSC_DISC_INV(I,j) = SML_DISC_INV(I,j)
          SML_DISC_INV(I,j) = 0
          DTL_SML(I,j) = 0.0
          DQW_SML(I,j) = 0.0
          K_CLOUD_DSCT(I,j) = K_CLOUD_TOP(I,j)
          K_CLOUD_TOP(I,j) = 0
        END IF
        IF ( .NOT.DSC(I,j) .AND. DSC_SAVE(I,j) ) THEN
!..decoupled layer removed in EXCF_NL; either...
          IF ( NTML_SAVE(I,j)  ==  NTML(I,j) ) THEN
!         !...had no turbulence forcing
            DTL_DSC(I,j) = 0.0
            DQW_DSC(I,j) = 0.0
            DSC_DISC_INV(I,j) = 0
            K_CLOUD_DSCT(I,j) = 0
          ELSE
!         !...recoupled with surface layer
            DTL_SML(I,j) = DTL_DSC(I,j)
            DQW_SML(I,j) = DQW_DSC(I,j)
            DTL_DSC(I,j) = 0.0
            DQW_DSC(I,j) = 0.0
            SML_DISC_INV(I,j) = DSC_DISC_INV(I,j)
            DSC_DISC_INV(I,j) = 0
            K_CLOUD_TOP(I,j) = K_CLOUD_DSCT(I,j)
            K_CLOUD_DSCT(I,j) = 0
          END IF
        END IF
      END DO
      END DO
!-----------------------------------------------------------------------
!L 6.  Calculate "explicit" entrainment fluxes of TL and QW.
!-----------------------------------------------------------------------
!..Integrate radiative divergence over layers:
!..  DF_NET (SML,DSC) are the total net radiative forcings of the
!..                   layers, used to calculate the total heat fluxes
!..                   at the discontinuous inversion levels
!..  F_NET (SML,DSC) are the net radiative divergences from the base
!..                  of the turbulently mixed layers, used to calculate
!..                  the appropriate grid-level turbulent fluxes
!..  All are in units of rho * Km/s
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
        DF_NET_SML(I,j) = 0.0
        DF_NET_DSC(I,j) = 0.0
        F_NET_SML(I,j,1) = 0.0
        F_NET_DSC(I,j,1) = 0.0
!       ! If no rad cooling in layer, need to set K_CLOUD_TOP to NTML
!       ! so that the correct entrainment heat flux is specified,
!       ! consistent with a linear total heat flux profile.
        IF ( K_CLOUD_TOP(I,j)  ==  0 ) K_CLOUD_TOP(I,j) = NTML(I,j)
        IF ( K_CLOUD_DSCT(I,j)  ==  0 ) K_CLOUD_DSCT(I,j) = NTDSC(I,j)
      END DO
      END DO
      DO K = 1,BL_LEVELS
      do j=1,rows
      do i=1,row_length
!..First the surface-forced mixed layer
          IF ( SML_DISC_INV(I,j)  ==  1 ) THEN
            IF ( K  <=  K_CLOUD_TOP(I,j) )                              &
     &        DF_NET_SML(I,j) = DF_NET_SML(I,j) + DF_OVER_CP(I,j,K)
          END IF
!..Second the decoupled mixed layer
          IF ( DSC_DISC_INV(I,j)  ==  1 ) THEN
            IF ( K  <=  K_CLOUD_DSCT(I,j) .AND. K  >=  NBDSC(I,j) )     &
     &         DF_NET_DSC(I,j) = DF_NET_DSC(I,j) + DF_OVER_CP(I,j,K)
          END IF
      END DO
      END DO
      END DO
!-----------------------------------------------------------------------
!..Assume DF_OVER_CP(K_cloud_top+2) is representative of clear-air rad
!..divergence and so extrapolate this down to discontinuous inversion
!..height and subtract this `clear-air' part from the grid-level
!..divergence.
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
!..First the surface-forced mixed layer
          IF ( SML_DISC_INV(I,j)  ==  1 ) THEN
            K = K_CLOUD_TOP(I,j)+1
            IF ( K  <   BL_LEVELS ) THEN
              DF_INV = DF_OVER_CP(I,j,K) - DF_OVER_CP(I,j,K+1) *        &
     &                       ( Z_HALF(I,j,NTML(I,j)+2) - ZH(I,j) )      &
     &                       / DZL(I,j,K+1)
            ELSEIF ( K  ==  BL_LEVELS ) THEN
              DF_INV = DF_OVER_CP(I,j,K)
            ELSE
              DF_INV = 0.0
            END IF
            DF_NET_SML(I,j) = DF_NET_SML(I,j) + MAX( DF_INV, 0.0 )
          END IF
!..Second the decoupled mixed layer
          IF ( DSC_DISC_INV(I,j)  ==  1 ) THEN
            K = K_CLOUD_DSCT(I,j)+1
            IF ( K  <   BL_LEVELS ) THEN
              DF_INV = DF_OVER_CP(I,j,K) - DF_OVER_CP(I,j,K+1) *        &
     &                       ( Z_HALF(I,j,NTDSC(I,j)+2) - ZHSC(I,j) )   &
     &                       / DZL(I,j,K+1)
            ELSEIF ( K  ==  BL_LEVELS ) THEN
              DF_INV = DF_OVER_CP(I,j,K)
            ELSE
              DF_INV = 0.0
            END IF
            DF_NET_DSC(I,j) = DF_NET_DSC(I,j) + MAX( DF_INV, 0.0 )
          END IF
      END DO
      END DO

!..Now the net radiative flux profiles from the layer base

      DO K = 2,BL_LEVELS
        do j=1,rows
        do i=1,row_length
!..First the surface-forced mixed layer
          IF ( SML_DISC_INV(I,j)  ==  1 ) THEN
             F_NET_SML(I,j,K) = F_NET_SML(I,j,K-1) + DF_OVER_CP(I,j,K-1)
          END IF
!..Second the decoupled mixed layer
          IF ( DSC_DISC_INV(I,j)  ==  1 ) THEN
            IF ( K  <   NBDSC(I,j) ) THEN
             F_NET_DSC(I,j,K) = 0.0
            ELSE
             F_NET_DSC(I,j,K) = F_NET_DSC(I,j,K-1) + DF_OVER_CP(I,j,K-1)
            END IF
          END IF
        END DO
        END DO
      END DO
!-----------------------------------------------------------------------
!..Specify entrainment fluxes at NTML+1 and NTDSC+1 directly through FTL
!..and FQW (and set the entrainment RHOKH to zero).
!..For QW, assume turbulent flux at subgrid ZH = - w_e * DQW_SML
!..The flux at the half-level below ZH can then be found by linear
!..interpolation between ZH and the surface (or base of the mixed-layer)
!..For TL, the procedure is basically the same, except that it is
!..the total (turb+rad) heat flux, TOTHF, that is linear.  Thus, TOTHF
!..at ZH is calculated, which is dependent on both the entrainment flux
!..and the net radiative divergence over the mixed layer, DF_NET_SML,
!..the linear profile of TOTHF is interpolated onto the half-level below
!..ZH and the radiative flux there, DF_NET_NTML, is sutracted off to
!..give the required turbulent component.  This has the happy
!..coincidence of ensuring that the `correct' entrainment flux is
!..specified for the given radiative flux (it even compensates for the
!..radiative flux having got a grid-level out, since radiation was last
!..called, for example).
!..For momentum, given the horizontal interpolation required, together
!..with the lack of accuracy in assuming a discontinuous inversion,
!..entrainment continues to be specified using the specified RHOKM.
!-----------------------------------------------------------------------
      DO K = 2,BL_LEVELS
      do j=1,rows
      do i=1,row_length
          FTL(I,j,K) = 0.0
          FQW(I,j,K) = 0.0
      END DO
      END DO
      END DO
      do j=1,rows
      do i=1,row_length

!..First the surface-based mixed layer (if entraining)
        K=NTML(I,j)+1
        KENT(I,j) = 2
        ZH_NP1(I,j) = 0.0
        DO IENT = 1,3
          T_FRAC(I,j,IENT) = 0.0
          ZRZI(I,j,IENT)   = 0.0
          WE_LIM(I,j,IENT) = 0.0
        END DO
        WE_PARM(I,j) = RDZ(I,j,K)*                                      &
     &                     ( RHOKH_TOP_ENT(I,j)+RHOKH_SURF_ENT(I,j) )   &
     &                                                / RHO_UV(I,j,K)

        IF ( SML_DISC_INV(I,j)  ==  1 .AND. .NOT.COUPLED(I,j) .AND.     &
     &       (RHOKH_TOP_ENT(I,j)+RHOKH_SURF_ENT(I,j))  >   0.0 ) THEN

          KENT(I,j) = K
!-----------------------------------------------------------------------
!..Calculate ZH at end of timestep, ZH_NP1
!-----------------------------------------------------------------------
!..linearly interpolate vertical velocity to ZH
          IF ( ZH(I,j)  >=  Z_FULL(I,j,K) ) THEN
            W_LS = W(I,j,K) + ( W(I,j,K+1) - W(I,j,K) )                 &
     &                  * (ZH(I,j)-Z_FULL(I,j,K)) * RDZ(I,j,K+1)
          ELSE
            W_LS = W(I,j,K) + ( W(I,j,K) - W(I,j,K-1) )                 &
     &                  * (ZH(I,j)-Z_FULL(I,j,K)) * RDZ(I,j,K)
          END IF
          W_LS = MIN ( W_LS, 0.0 )
!         ! only interested in subsidence

          ZH_NP1(I,j) = ZH(I,j) +                                       &
     &                    TIMESTEP * ( WE_PARM(I,j) + W_LS )
          ZH_NP1(I,j) = MAX( ZH_NP1(I,j), Z_HALF(I,j,K-1) )
          IF ( ZH_NP1(I,j)  >   Z_TOP(I,j,K+1) ) THEN
!           ! limit ZH and W_e (and therefore the entraiment fluxes)
!           ! because the inversion cannot rise more than one level
!           ! in a timestep.
            ZH_NP1(I,j) = Z_TOP(I,j,K+1)
            WE_PARM(I,j) =                                              &
     &              (Z_TOP(I,j,K+1) - ZH(I,j))/TIMESTEP - W_LS
          END IF

!..Linearly interpolate between the known total (turb+rad) TL fluxes at
!..the surface and the mean height of the discontinuous inversion,
!..taking care over where the inversion is during the timestep.

          IF ( ZH_NP1(I,j)  >   Z_HALF(I,j,K+1) ) THEN
!    ! ZH risen above level K+1 so specify appropriate flux at this leve

!           ! Adjust the entraiment rate to allow for that implied by
!           ! the subsidence increment applied below (combine increments
!           ! from levels NTML and NTML+1 to give the subsidence
!           ! associated with the inversion.
            TLS_INV(I,j) = TLS_INC(I,j,K) + TLS_INC(I,j,K-1)
            W_S_ENT = 0.0
            IF ( DTL_SML(I,j)  /=  0.0 ) W_S_ENT =                      &
     &             MIN( 0.0, -TLS_INV(I,j) * DZL(I,j,K) /DTL_SML(I,j) )
!           ! Only allow w_e to be reduced to zero!
            WE_LIM(I,j,3) = RHO_UV(I,j,K+1) *                           &
     &                      MAX( 0.0, WE_PARM(I,j) + W_S_ENT )

            TOTHF_INV = - WE_LIM(I,j,3)*DTL_SML(I,j) + DF_NET_SML(I,j)
!    ! T_FRAC is fraction of timestep inversion spends above Z_HALF(K+1)
            T_FRAC(I,j,3) = (ZH_NP1(I,j)-Z_HALF(I,j,K+1)) /             &
     &                      (ZH_NP1(I,j)-ZH(I,j))
!    ! Estimate average flux over timestep as flux for average ZH
            ZRZI(I,j,3)  = Z_HALF(I,j,K+1)*2.0/                         &
     &                     (Z_HALF(I,j,K+1)+ZH_NP1(I,j))

            FTL(I,j,K+1) = T_FRAC(I,j,3) * (                            &
     &         FTL(I,j,1) + ( TOTHF_INV - FTL(I,j,1) )*ZRZI(I,j,3)      &
     &                                          - F_NET_SML(I,j,K+1) )
            RHOKH_TOP(I,j,K+1) = 0.0   ! apply entrainment explicitly
            RHOKH(I,j,K+1) = 0.0       !      "

!    ! Layer should now be well-mixed to NTML+1 (by end of timestep)
!    ! so raise NTML by one (this means gradient-adjustment is also
!    ! applied at half-level old_NTML+1).
!    ! Note KH profiles should already be calculated at level NTML+1
!    ! because ZH is above this level.
            NTML(I,j) = NTML(I,j) + 1

          ELSE  ! ZH always below half-level K+1

           RHOKH_TOP(I,j,K) = 0.0   ! apply entrainment explicitly
           RHOKH(I,j,K) = 0.0       !      "

           IF ( ZH_NP1(I,j)  >=  Z_HALF(I,j,K) ) THEN
!..ZH always above level K so specify full entrainment flux here

!    ! Adjust the entraiment rate to allow for that impled by the
!    ! subsidence increment applied below.
            W_S_ENT = 0.0
            IF ( DTL_SML(I,j)  /=  0.0 ) W_S_ENT =                      &
     &        MIN( 0.0, -TLS_INC(I,j,K-1) * DZL(I,j,K-1) /DTL_SML(I,j) )
            WE_LIM(I,j,2) = RHO_UV(I,j,K) *                             &
     &                      MAX( 0.0, WE_PARM(I,j) + W_S_ENT )

            TOTHF_INV = - WE_LIM(I,j,2)*DTL_SML(I,j) + DF_NET_SML(I,j)
            T_FRAC(I,j,2) = 1.0
            ZRZI(I,j,2) = Z_HALF(I,j,K)*2.0/(ZH(I,j)+ZH_NP1(I,j))

            FTL(I,j,K) = FTL(I,j,1) + ( TOTHF_INV - FTL(I,j,1) )        &
     &                               *ZRZI(I,j,2) - F_NET_SML(I,j,K)
           ELSE
!..ZH dropped below level K so specify full entrainment flux at K-1...
            IF (K-1  >=  2) THEN     ! FTL(K=1) is surface flux!
             NTML(I,j) = NTML(I,j) - 1
             W_S_ENT = 0.0
             IF ( DTL_SML(I,j)  /=  0.0 ) W_S_ENT =                     &
     &        MIN( 0.0, -TLS_INC(I,j,K-2) * DZL(I,j,K-2) /DTL_SML(I,j) )
             WE_LIM(I,j,1) = RHO_UV(I,j,K-1) *                          &
     &                              MAX( 0.0, WE_PARM(I,j)+W_S_ENT )
             TOTHF_INV = - WE_LIM(I,j,1)*DTL_SML(I,j) + DF_NET_SML(I,j)
             T_FRAC(I,j,1) = 1.0
             ZRZI(I,j,1) = Z_HALF(I,j,K-1)*2.0/(ZH(I,j)+ZH_NP1(I,j))
             FTL(I,j,K-1) = FTL(I,j,1)                                  &
     &                     + (TOTHF_INV - FTL(I,j,1)) * ZRZI(I,j,1)     &
     &                                       - F_NET_SML(I,j,K-1)
             RHOKH_TOP(I,j,K-1) = 0.0   ! apply entrainment explicitly
             RHOKH(I,j,K-1) = 0.0      !      "
            END IF
!    ! ...and set specified flux at level K allowing for time spent
!    ! with ZH above this level

            W_S_ENT = 0.0
            IF ( DTL_SML(I,j)  /=  0.0 ) W_S_ENT =                      &
     &        MIN( 0.0, -TLS_INC(I,j,K-1) * DZL(I,j,K-1) /DTL_SML(I,j) )
            WE_LIM(I,j,2) = RHO_UV(I,j,K) *                             &
     &                       MAX( 0.0, WE_PARM(I,j) + W_S_ENT )
            TOTHF_INV = - WE_LIM(I,j,2)*DTL_SML(I,j) + DF_NET_SML(I,j)
!    ! T_FRAC is fraction of timestep inversion spends above Z_HALF(K)
            T_FRAC(I,j,2) = (ZH(I,j)-Z_HALF(I,j,K))                     &
     &                                 / (ZH(I,j)-ZH_NP1(I,j))
            ZRZI(I,j,2) = Z_HALF(I,j,K)*2.0/(Z_HALF(I,j,K)+ZH(I,j))
            FTL(I,j,K) = T_FRAC(I,j,2) * (                              &
     &         FTL(I,j,1) + (TOTHF_INV - FTL(I,j,1)) * ZRZI(I,j,2)      &
     &                                           - F_NET_SML(I,j,K) )
           END IF

          END IF  ! ZH gone above half-level K+1.

        ELSE   ! NOT specifying entrainment flux but KH

!         ! Add entrainment KH to K-profiles
!         ! (for COUPLED layers these will be zero)
          RHOKH_TOP(I,j,K) = RHOKH_TOP(I,j,K) + RHOKH_TOP_ENT(I,j)
          RHOKH(I,j,K) = RHOKH(I,j,K) + RHOKH_SURF_ENT(I,j)

        END IF  ! IF NOT CUMULUS OR COUPLED
!-------------------------------------------------
!..Second the decoupled mixed layer, if entraining
!-------------------------------------------------
        K=NTDSC(I,j)+1
        KENT_DSC(I,j) = 2
        ZHSC_NP1(I,j) = 0.0
        DO IENT = 1,3
          T_FRAC_DSC(I,j,IENT) = 0.0
          ZRZI_DSC(I,j,IENT)   = 0.0
          WE_LIM_DSC(I,j,IENT) = 0.0
        END DO
        WE_DSC_PARM(I,j) = RDZ(I,j,K)*RHOKH_DSCT_ENT(I,j)               &
     &                                             / RHO_UV(I,j,K)

        IF ( DSC_DISC_INV(I,j)  ==  1                                   &
     &                .AND. RHOKH_DSCT_ENT(I,j)  >   0.0 ) THEN

          KENT_DSC(I,j) = K
!-----------------------------------------------------------------------
!..Calculate ZHSC at end of timestep, ZHSC_NP1
!-----------------------------------------------------------------------
!..interpolate vertical velocity to ZH
          IF ( ZHSC(I,j)  >=  Z_FULL(I,j,K) ) THEN
            W_LS_DSC = W(I,j,K) + ( W(I,j,K+1) - W(I,j,K) ) *           &
     &                       (ZHSC(I,j)-Z_FULL(I,j,K)) * RDZ(I,j,K+1)
          ELSE
            W_LS_DSC = W(I,j,K) + ( W(I,j,K) - W(I,j,K-1) ) *           &
     &                       (ZHSC(I,j)-Z_FULL(I,j,K)) * RDZ(I,j,K)
          END IF
          W_LS_DSC = MIN ( W_LS_DSC, 0.0 )
!         ! only interested in subsidence

          ZHSC_NP1(I,j) = ZHSC(I,j) +                                   &
     &          TIMESTEP * ( WE_DSC_PARM(I,j) + W_LS_DSC )
          ZHSC_NP1(I,j) = MAX( ZHSC_NP1(I,j), Z_HALF(I,j,K-1) )
          IF ( ZHSC_NP1(I,j)  >   Z_TOP(I,j,K+1) ) THEN
!           ! limit ZHSC and W_e (and therefore the entrainment fluxes)
!           ! because the inversion cannot rise more than one level
!           ! in a timestep.
            ZHSC_NP1(I,j) = Z_TOP(I,j,K+1)
            WE_DSC_PARM(I,j) =                                          &
     &         (Z_TOP(I,j,K+1) - ZHSC(I,j))/TIMESTEP - W_LS_DSC
          END IF

!..linearly interpolate between the known total (turb+rad) TL fluxes at
!..the base (assumed zero) and the mean height of the discontinuous
!..inversion, taking care over where the inversion is during the
!..timestep.  Assume layer base is fixed during timestep.

          IF ( ZHSC_NP1(I,j)  >   Z_HALF(I,j,K+1) ) THEN
!    ! ZHSC risen above level K+1 so specify approp. flux at this level

!          ! Adjust the entraiment rate to allow for that implied by
!          ! the subsidence increment applied below (combine  increments
!          ! from levels NTDSC and NTDSC+1)
           TLS_INV(I,j) = TLS_INC(I,j,K) + TLS_INC(I,j,K-1)
           W_S_ENT = 0.0
           IF ( DTL_DSC(I,j)  /=  0.0 ) W_S_ENT =                       &
     &              MIN( 0.0, -TLS_INV(I,j) * DZL(I,j,K) /DTL_DSC(I,j) )
!          ! Only allow w_e to be reduced to zero!
           WE_LIM_DSC(I,j,3) = RHO_UV(I,j,K+1) *                        &
     &                     MAX( 0.0, WE_DSC_PARM(I,j) + W_S_ENT )

           TOTHF_INV = - WE_LIM_DSC(I,j,3)*DTL_DSC(I,j) +               &
     &                                                   DF_NET_DSC(I,j)
           ZRZI_DSC(I,j,3) = ( Z_HALF(I,j,K+1)-                         &
     &                                     (ZHSC(I,j)-DSCDEPTH(I,j)) )/ &
     &       ( 0.5*(ZHSC_NP1(I,j)+Z_HALF(I,j,K+1))-                     &
     &                                     (ZHSC(I,j)-DSCDEPTH(I,j)) )
           T_FRAC_DSC(I,j,3) = (ZHSC_NP1(I,j)-Z_HALF(I,j,K+1)) /        &
     &                                       (ZHSC_NP1(I,j)-ZHSC(I,j))
           FTL(I,j,K+1) = T_FRAC_DSC(I,j,3) * (                         &
     &                TOTHF_INV * ZRZI_DSC(I,j,3) - F_NET_DSC(I,j,K+1) )
           RHOKH_TOP(I,j,K+1) = 0.0   ! apply entrainment explicitly
           RHOKH(I,j,K+1) = 0.0       !      "

!    ! Layer should now be well-mixed to NTDSC+1 (by end of timestep)
!    ! so raise NTDSC by one:
           NTDSC(I,j) = NTDSC(I,j) + 1
!    ! Note KH profiles should already be calculated at level
!    ! NTDSC_old+1 because ZHSC is above this level.

          ELSE  ! ZHSC always below half-level K+1

           RHOKH(I,j,K) = 0.0       ! apply entrainment explicitly
           RHOKH_TOP(I,j,K) = 0.0   !      "
           IF ( ZHSC_NP1(I,j)  >=  Z_HALF(I,j,K) ) THEN
!..ZH always above level K so specify full entrainment flux here
            W_S_ENT = 0.0
            IF ( DTL_DSC(I,j)  /=  0.0 ) W_S_ENT =                      &
     &        MIN( 0.0, -TLS_INC(I,j,K-1) * DZL(I,j,K-1) /DTL_DSC(I,j) )
            WE_LIM_DSC(I,j,2) = RHO_UV(I,j,K) *                         &
     &                      MAX( 0.0, WE_DSC_PARM(I,j) + W_S_ENT )
            TOTHF_INV = - WE_LIM_DSC(I,j,2)*DTL_DSC(I,j) +              &
     &                                                 DF_NET_DSC(I,j)
            ZRZI_DSC(I,j,2) =( Z_HALF(I,j,K)-(ZHSC(I,j)-DSCDEPTH(I,j)) )&
     &                  /( DSCDEPTH(I,j)+0.5*(ZHSC_NP1(I,j)-ZHSC(I,j)) )
!           ! ZRZI_DSC = (z - z_base) / (zhsc - z_base)
            T_FRAC_DSC(I,j,2) = 1.0

            FTL(I,j,K) = TOTHF_INV * ZRZI_DSC(I,j,2) - F_NET_DSC(I,j,K)

           ELSE
!..ZH dropped below level K so specify full entrainment flux at K-1...

            NTDSC(I,j) = NTDSC(I,j) - 1  ! could reduce NTDSC to 1
            W_S_ENT = 0.0
            IF ( DTL_DSC(I,j)  /=  0.0 ) W_S_ENT = MIN( 0.0,            &
     &                 -TLS_INC(I,j,K-2) * DZL(I,j,K-2) /DTL_DSC(I,j) )
            WE_LIM_DSC(I,j,1) = RHO_UV(I,j,K-1) *                       &
     &                         MAX( 0.0, WE_DSC_PARM(I,j) + W_S_ENT )

            TOTHF_INV = - WE_LIM_DSC(I,j,1)*DTL_DSC(I,j) +              &
     &                                                   DF_NET_DSC(I,j)
            ZRZI_DSC(I,j,1) =                                           &
     &               ( Z_HALF(I,j,K-1)-(ZHSC(I,j)-DSCDEPTH(I,j)) )      &
     &                  /( DSCDEPTH(I,j)+0.5*(ZHSC_NP1(I,j)-ZHSC(I,j)) )
            T_FRAC_DSC(I,j,1) = 1.0
            FTL(I,j,K-1) = TOTHF_INV*ZRZI_DSC(I,j,1) -F_NET_DSC(I,j,K-1)

            RHOKH_TOP(I,j,K-1) = 0.0   ! apply entrainment explicitly
            RHOKH(I,j,K-1) = 0.0       !      "
!    ! ...and set specified flux at level K for time spent
!    ! with ZH above this level
            W_S_ENT = 0.0
            IF ( DTL_DSC(I,j)  /=  0.0 ) W_S_ENT = MIN( 0.0,            &
     &                  -TLS_INC(I,j,K-1) * DZL(I,j,K-1) /DTL_DSC(I,j) )
            WE_LIM_DSC(I,j,2) = RHO_UV(I,j,K) *                         &
     &                   MAX( 0.0, WE_DSC_PARM(I,j) + W_S_ENT )
            TOTHF_INV = - WE_LIM_DSC(I,j,2)*DTL_DSC(I,j) +              &
     &                                                   DF_NET_DSC(I,j)
            T_FRAC_DSC(I,j,2) = (ZHSC(I,j)-Z_HALF(I,j,K)) /             &
     &                                  (ZHSC(I,j)-ZHSC_NP1(I,j))
            ZRZI_DSC(I,j,2) =( Z_HALF(I,j,K)-(ZHSC(I,j)-DSCDEPTH(I,j)) )&
     &                  /( DSCDEPTH(I,j)+0.5*(ZHSC(I,j)-Z_HALF(I,j,K)) )
            FTL(I,j,K) = T_FRAC_DSC(I,j,2) * (                          &
     &                  TOTHF_INV * ZRZI_DSC(I,j,2) - F_NET_DSC(I,j,K) )
           END IF

          END IF  ! ZHSC gets above or stays below half-level K+1

        ELSEIF ( DSC(I,j) ) THEN

!         ! not specifying entrainment flux but KH
          RHOKH_TOP(I,j,K) = RHOKH_DSCT_ENT(I,j)

        END IF  ! IF NOT DSC

      END DO
      END DO
!-----------------------------------------------------------------------
! Specify QW entrainment fluxes
! (must be separate loops to avoid overwriting)
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
        DO IENT = 1,3
          IF ( KENT(I,j)-2+IENT  >   1 )                                &
     &      FQW(I,j,KENT(I,j)-2+IENT) = T_FRAC(I,j,IENT)*( FQW(I,j,1)   &
     &               - ( WE_LIM(I,j,IENT)*DQW_SML(I,j) + FQW(I,j,1) )   &
     &                                   *ZRZI(I,j,IENT) )

        END DO
        DO IENT = 1,3
          IF ( KENT_DSC(I,j) >= 3 .AND. T_FRAC_DSC(I,j,IENT) >  0.0 )   &
     &      FQW(I,j,KENT_DSC(I,j)-2+IENT) =  - T_FRAC_DSC(I,j,IENT) *   &
     &          WE_LIM_DSC(I,j,IENT) * DQW_DSC(I,j) * ZRZI_DSC(I,j,IENT)
        END DO
      END DO
      END DO
!-----------------------------------------------------------------------
!- Update standard deviations and gradient adjustment to use this
!- timestep's ZH (code from SF_EXCH)
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
        W_S_CUBED = 0.25 * ZH(I,j) * FB_SURF(I,j)
        IF (W_S_CUBED  >   0.0) THEN
          W_M  =                                                        &
     &   ( W_S_CUBED + V_S(I,j) * V_S(I,j) * V_S(I,j) ) ** (1.0/3.0)
          T1_SD(I,j) = 1.93 * FTL(I,j,1) / (RHOSTAR_GB(I,j) * W_M)
          Q1_SD(I,j) = 1.93 * FQW(I,j,1) / (RHOSTAR_GB(I,j) * W_M)
          TV1_SD(I,j) = T(I,j,1) *                                      &
     &        ( 1.0 + C_VIRTUAL*Q(I,j,1) - QCL(I,j,1) - QCF(I,j,1) ) *  &
     &        ( BT(I,j,1)*T1_SD(I,j) + BQ(I,j,1)*Q1_SD(I,j) )
          T1_SD(I,j) = MAX ( 0.0 , T1_SD(I,j) )
          Q1_SD(I,j) = MAX ( 0.0 , Q1_SD(I,j) )
          IF (TV1_SD(I,j)  <=  0.0) THEN
            TV1_SD(I,j) = 0.0
            T1_SD(I,j) = 0.0
            Q1_SD(I,j) = 0.0
          END IF
        END IF
        GRAD_T_ADJ(I,j) = MIN( MAX_T_GRAD ,                             &
     &                       A_GRAD_ADJ * T1_SD(I,j) / ZH(I,j) )
!        IF (T1_SD(I,j)  >   0.0) THEN
!          GRAD_Q_ADJ(I,j) = (Q1_SD(I,j) / T1_SD(I,j)) * GRAD_T_ADJ(I,j)
!        ELSE
          GRAD_Q_ADJ(I,j) = 0.0
!        END IF
      END DO
      END DO
!-----------------------------------------------------------------------
!- Save diagnostics
!-----------------------------------------------------------------------
        do j=1,rows
        do i=1,row_length
          IF ( DSC(I,j) ) THEN
            IF (BL_diag%L_dscbase) THEN
              BL_diag%dscbase(i,j)= ZHSC(I,j)-DSCDEPTH(I,j)
            END IF
            IF (BL_diag%L_cldbase) THEN
              BL_diag%cldbase(i,j)= ZHSC(I,j)-ZC_DSC(I,j)
            END IF
            IF (BL_diag%L_weparm_dsc) THEN
              BL_diag%weparm_dsc(i,j)= WE_DSC_PARM(I,j)
            END IF
          ELSE
            IF (BL_diag%L_dscbase) THEN
              BL_diag%dscbase(i,j)= RMDI
            END IF
            IF (BL_diag%L_cldbase) THEN
              BL_diag%cldbase(i,j)= ZH(I,j)-ZC(I,j)
            END IF
            IF (BL_diag%L_weparm_dsc) THEN
              BL_diag%weparm_dsc(i,j)= WE_PARM(I,j)
            END IF
          END IF
          IF (BL_diag%L_weparm) THEN
            BL_diag%weparm(i,j)= WE_PARM(I,j)
          END IF
        END DO
        END DO

#if defined(SCMA)
!
!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
      IF (L_SCMDiags(SCMDiag_bl)) THEN

!------------------------------------------------------
        DO J=1,rows
        DO I=1,row_length
          IF ( DSC(I,j) ) THEN
            SCMOUT(I,j)= ZHSC(I,j)-DSCDEPTH(I,j)
          ELSE
            SCMOUT(i,j)= ZH(I,j)
          END IF
        END DO
        END DO
! DEPENDS ON: scmoutput
        Call SCMoutput(SCMOUT,                                          &
             'DSCbase','Base of decoupled layer','m',                   &
             t_avg,d_sl,default_streams,'',RoutineName)
!------------------------------------------------------
        DO J=1,rows
        DO I=1,row_length
          IF ( DSC(I,j) ) THEN
            SCMOUT(I,j)= ZHSC(I,j)-ZC_DSC(I,j)
          ELSE
            SCMOUT(i,j)= ZH(I,j)-ZC(I,j)
          END IF
        END DO
        END DO
! DEPENDS ON: scmoutput
        Call SCMoutput(SCMOUT,                                          &
             'ScCldBase','Stratocumulus cloud base','m',                &
             t_avg,d_sl,default_streams,'',RoutineName)
!------------------------------------------------------
! DEPENDS ON: scmoutput
        Call SCMoutput(WE_PARM,                                         &
             'Entr_SML','SML-top entrainment rate','m/s',               &
             t_avg,d_sl,default_streams,'',RoutineName)
!------------------------------------------------------
        DO J=1,rows
        DO I=1,row_length
          IF ( DSC(I,j) ) THEN
            SCMOUT(I,j)= WE_DSC_PARM(I,j)
          ELSE
            SCMOUT(i,j)= WE_PARM(I,j)
          END IF
        END DO
        END DO
! DEPENDS ON: scmoutput
        Call SCMoutput(SCMOUT,                                          &
             'Entr_BL','BL-top entrainment rate','m/s',                 &
             t_avg,d_sl,default_streams,'',RoutineName)
!------------------------------------------------------

      END IF ! L_SCMDiags(SCMDiag_bl)
#endif


      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('KMKHZ   ',4)
      END IF

      RETURN
      END SUBROUTINE KMKHZ
#endif
