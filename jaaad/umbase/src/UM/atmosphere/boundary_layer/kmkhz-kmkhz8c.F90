#if defined(A03_8C)
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
     & DECFIX, STOPWE_SBL, NG_STRESS, FLUX_GRAD, Lq_mix_bl,             &
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
     & KENT, WE_LIM, T_FRAC_TR, ZRZI_TR,                                &
     & KENT_DSC, WE_LIM_DSC, T_FRAC_DSC_TR, ZRZI_DSC_TR,                &
     & BL_diag,nSCMDpkgs,L_SCMDiags,                                    &
     & LTIMER                                                           &
     & )

      Use cv_run_mod, Only:                                             &
          l_conv4a

      Use bl_diags_mod, Only :                                          &
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
     &,FLUX_GRAD                                                        &
                              ! IN switch for flux-gradient formulation
     &,DECFIX                                                           &
                              ! IN switch not used in 8C version
     &,STOPWE_SBL             !          "

      LOGICAL                                                           &
     & LTIMER                                                           &
                              ! IN Flag for TIMER diagnostics
     &,Lq_mix_bl              ! IN True if using mixing ratios

      REAL                                                              &
     & TIMESTEP                                                         &
                              ! IN Model timestep (s)
     &,ZMAX_FOR_DSC
!                             ! IN Max height to look for DSC cloud-base

      INTEGER                                                           &
     & NTPAR(row_length,rows)                                           &
                              ! IN Top level of parcel ascent.
!                             !    Used in convection scheme.
!                             !    NOTE: CAN BE > BL_LEVELS-1
     &,NLCL(row_length,rows)  ! IN No. of model layers below the
!                             !    lifting condensation level.

      LOGICAL                                                           &
     & L_SHALLOW(row_length,rows)! IN Flag to indicate shallow
!                                !    convection (only for A05_4A)

      REAL                                                              &
     & BQ(row_length,rows,BL_LEVELS)                                    &
                              ! IN A buoyancy parameter for clear air
!                             !    on p,T,q-levels (full levels).
     &,BT(row_length,rows,BL_LEVELS)                                    &
                              ! IN A buoyancy parameter for clear air
!                             !    on p,T,q-levels (full levels).
     &,BQM(row_length,rows,BL_LEVELS)                                   &
                              ! IN A buoyancy parameter for clear air
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,BTM(row_length,rows,BL_LEVELS)                                   &
                              ! IN A buoyancy parameter for clear air
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,BQM_CLD(row_length,rows,BL_LEVELS)                               &
                              ! IN A buoyancy parameter for cloudy air
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,BTM_CLD(row_length,rows,BL_LEVELS)                               &
                              ! IN A buoyancy parameter for cloudy air
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,A_QS(row_length,rows,BL_LEVELS)                                  &
                              ! IN Saturated lapse rate factor
!                             !    on p,T,q-levels (full levels).
     &,A_QSM(row_length,rows,BL_LEVELS)                                 &
                              ! IN Saturated lapse rate factor
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,A_DQSDTM(row_length,rows,BL_LEVELS)                              &
                              ! IN Saturated lapse rate factor
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,DQSDT(row_length,rows,BL_LEVELS)                                 &
                              ! IN Partial derivative of QSAT w.r.t.
!                             !    temperature.
     &,CFM (row_length,rows,BL_LEVELS)                                  &
                              ! IN Estimate of cloud fraction
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,P(row_length,rows,BL_LEVELS)                                     &
                              ! IN P(*,K) is pressure at full level k.
     &,P_HALF(row_length,rows,BL_LEVELS)                                &
                              ! IN P_HALF(*,K) is pressure at half
!                             !    level k-1/2.
     &,QW(row_length,rows,BL_LEVELS)                                    &
                              ! IN Total water content (kg per kg air).
     &,TL(row_length,rows,BL_LEVELS)                                    &
                              ! IN Liquid/frozen water temperature (K).
     &,T(row_length,rows,BL_LEVELS)                                     &
                              ! IN Temperature (K).
     &,QCF(row_length,rows,BL_LEVELS)                                   &
                              ! IN Cloud ice (kg per kg air)
     &,QCL(row_length,rows,BL_LEVELS)                                   &
                              ! IN Cloud liquid water
     &,Q(row_length,rows,BL_LEVELS)                                     &
                              ! IN specific humidity
     &,CF(row_length,rows,BL_LEVELS)                                    &
                              ! IN Cloud fractions for boundary levs.
     &,Z_FULL(row_length,rows,BL_LEVELS)                                &
                              ! IN Z_FULL(*,K) is the height of the
!                             !    k-th full level above the surface.
     &,Z_HALF(row_length,rows,BL_LEVELS)                                &
                              ! IN Z_HALF(*,K) is the height of level
!                             !       k-1/2 above the surface (m).
     &,Z_UV(row_length,rows,BL_LEVELS)                                  &
                              ! IN For a vertically staggered grid
!                             !    with a u,v-level first above the
!                             !    surface, Z_UV(*,K) is the height
!                             !    of the k-th u,v-level (half level
!                             !    k-1/2) above the surface
     &,Z_TQ(row_length,rows,BL_LEVELS)                                  &
                              ! IN For a vertically staggered grid
!                             !    with a u,v-level first above the
!                             !    surface, Z_TQ(*,K) is the height
!                             !    of the k-th T,q-level (full level
!                             !    k) above the surface
     &,DZL(row_length,rows,BL_LEVELS)                                   &
                              ! IN Layer depths (m).  DZL(,K) is the
!                             !    distance from layer boundary K-1/2
!                             !    to layer boundary K+1/2.  For K=1
!                             !    the lower boundary is the surface.
     &,RDZ(row_length,rows,BL_LEVELS)                                   &
                              ! IN Reciprocal of distance between
!                             !    full levels (m-1).  1/RDZ(,K) is
!                             !    the vertical distance from level
!                             !    K-1 to level K, except that for
!                             !    K=1 it is the height of the
!                             !    lowest atmospheric full level.
     &,RHO_UV(row_length,rows,BL_LEVELS)                                &
                              ! IN density on UV (ie. rho) levels,
!                             !    used in RHOKH so dry density if
!                             !    Lq_mix_bl is true
     &,RHO_TQ(row_length,rows,BL_LEVELS)                                &
                              ! IN density on TQ (ie. theta) levels,
!                             !    used in RHOKM so wet density
     &,RHO_DRY_TQ(row_length,rows,BL_LEVELS)
                              ! IN density on TQ (ie. theta) levels,
!                             !    used in non-turb flux integration
!                             !    so dry density if Lq_mix_bl is true

      REAL                                                              &
     & V_S(row_length,rows)                                             &
                              ! IN Surface friction velocity (m/s)
     &,FB_SURF(row_length,rows)                                         &
                              ! IN Surface buoyancy flux over density
!                             !       (m^2/s^3).
     &,RHOSTAR_GB(row_length,rows)                                      &
                              ! IN Surface air density in kg/m3
     &,T1_SD(row_length,rows)                                           &
                              ! IN Standard Deviation of level 1
!                             !    temperature (K).
     &,Q1_SD(row_length,rows)                                           &
                              ! IN Standard Deviation of level 1
!                             !    specific humidity (kg/kg).
     &,Z0M(row_length,rows)                                             &
                              ! IN Roughness length for momentum (m).
     &,Z_LCL(row_length,rows)                                           &
                              ! IN Height of lifting condensation
!                             !    level.
     &,ZHPAR(row_length,rows)                                           &
                              ! IN Height of top of NTPAR
!                             !    NOTE: CAN BE ABOVE BL_LEVELS-1
     &,ZH_PREV(row_length,rows)
!                             ! IN boundary layer height (m) from
!                             !    previous timestep

      REAL                                                              &
     & RAD_HR(row_length,rows,BL_LEVELS,2)                              &
                              ! IN (LW,SW) radiative heating rates (K/s)
     &,MICRO_TENDS(row_length,rows,BL_LEVELS,2)
                              ! IN Tendencies from microphysics
!                             !    (TL, K/s; QW, kg/kg/s)

      REAL                                                              &
     & W(row_length,rows,0:BL_LEVELS)                                   &
                              ! IN Vertical velocity (m/s)
     &,ETADOT(row_length,rows,0:BL_LEVELS)
                              ! IN d(ETA)/dt

      REAL                                                              &
     & eta_theta_levels(0:bl_levels)
                              ! IN eta coordinate

      INTEGER                                                           &
     & NTML(row_length,rows)  ! INOUT Number of model levels in the
!                             !       turbulently mixed layer.

      LOGICAL                                                           &
     & CUMULUS(row_length,rows)
!                             ! INOUT Flag for Cu in the bl

      REAL                                                              &
     & ZH(row_length,rows)                                              &
                              ! INOUT Boundary layer height (m).
     &,FTL(row_length,rows,BL_LEVELS)                                   &
     &,FQW(row_length,rows,BL_LEVELS)
!                             ! INOUT "Explicit" fluxes of TL and QW
!                             !       (rho*Km/s, rho*m/s)
!                             !       IN:  level 1 (surface flux)
!                             !       OUT: entrainment-level flux

      INTEGER                                                           &
     & NTDSC(row_length,rows)                                           &
                              ! OUT Top level for turb mixing in
!                             !       cloud layer
     &,NBDSC(row_length,rows)                                           &
                              ! OUT Bottom level of any decoupled
!                             !       turbulently mixed Sc layer
     &,SML_DISC_INV(row_length,rows)                                    &
                              ! OUT Flags for whether discontinuous
     &,DSC_DISC_INV(row_length,rows)                                    &
                              ! OUT   inversions are diagnosed
     &,KENT(row_length,rows)                                            &
                              ! OUT Grid-levels of SML and DSC
     &,KENT_DSC(row_length,rows)
                              ! OUT   inversions (for tracer mixing)

      LOGICAL                                                           &
     & UNSTABLE(row_length,rows)                                        &
                              ! OUT Flag to indicate an unstable
                              !     surface layer.
     &,DSC(row_length,rows)                                             &
                              ! OUT Flag set if decoupled
                              !     stratocumulus layer found
     &,COUPLED(row_length,rows)
                              ! OUT Flag to indicate Sc layer weakly
                              !     decoupled (implies mixing at SML
                              !     top is through K profiles rather
                              !     than entrainment parametrization)

      REAL                                                              &
     & RHOKM(row_length,rows, 2:BL_LEVELS)                              &
                              ! OUT Non-local turbulent mixing
!                             !     coefficient for momentum.
     &,RHOKH(row_length,rows, 2:BL_LEVELS)                              &
                              ! OUT Non-local turbulent mixing
!                             !     coefficient for scalars.
     &,RHOKM_TOP(row_length,rows, 2:BL_LEVELS)                          &
                              ! OUT Top-down turbulent mixing
!                             !     coefficient for momentum.
     &,RHOKH_TOP(row_length,rows, 2:BL_LEVELS)                          &
                              ! OUT Top-down turbulent mixing
!                             !     coefficient for scalars.
     &,F_NGSTRESS(row_length,rows, 2:BL_LEVELS)                         &
                              ! OUT dimensionless function for
!                             !     non-gradient stresses
     &,RHOF2(row_length,rows, 2:BL_LEVELS)                              &
     &,RHOFSC(row_length,rows, 2:BL_LEVELS)
                              ! OUT f2 and fsc term shape profiles
                              !       multiplied by rho

      REAL                                                              &
     & FT_NT(row_length,rows, BL_LEVELS+1)                              &
     &,FQ_NT(row_length,rows, BL_LEVELS+1)
                              ! OUT Non-turbulent heat and moisture
                              !       fluxes (rho*Km/s, rho*m/s)

      REAL                                                              &
     & TOTHF_ZH(row_length,rows)                                        &
                              ! OUT Total heat fluxes at inversions
     &,TOTHF_ZHSC(row_length,rows)                                      &
                              !     (rho*Km/s)
     &,TOTQF_ZH(row_length,rows)                                        &
                              ! OUT Total moisture fluxes at
     &,TOTQF_ZHSC(row_length,rows)                                      &
                              !      inversions (rho*m/s)
     &,FT_NT_DSCB(row_length,rows)                                      &
                              ! OUT Non-turbulent heat and moisture
     &,FQ_NT_DSCB(row_length,rows)                                      &
                              !      flux at the base of the DSC layer
     &,GRAD_T_ADJ(row_length,rows)                                      &
                              ! OUT Temperature gradient adjustment
!                             !     for non-local mixing in unstable
!                             !     turbulent boundary layer.
     &,GRAD_Q_ADJ(row_length,rows)                                      &
                              ! OUT Humidity gradient adjustment
!                             !     for non-local mixing in unstable
!                             !     turbulent boundary layer.
     &,ZHSC(row_length,rows)
                              ! OUT Cloud layer height (m).

!     ! The following are used in tracer mixing.
!     ! At 8B 3 elements were used - here just (i,j,2) is used.
      REAL                                                              &
     & WE_LIM(row_length,rows,3)                                        &
                              ! OUT rho*entrainment rate implied by
!                             !     placing of subsidence
     &,ZRZI_TR(row_length,rows,3)                                       &
                              ! OUT (z-z_base)/(z_i-z_base)
     &,T_FRAC_TR(row_length,rows,3)                                     &
                              ! OUT a fraction of the timestep
     &,WE_LIM_DSC(row_length,rows,3)                                    &
                              ! OUT rho*entrainment rate implied by
!                             !     placing of subsidence
     &,ZRZI_DSC_TR(row_length,rows,3)                                   &
                              ! OUT (z-z_base)/(z_i-z_base)
     &,T_FRAC_DSC_TR(row_length,rows,3)
                              ! OUT a fraction of the timestep

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
!
#if defined(SCMA)
#include "s_scmop.h"
#endif
#include "c_lheat.h"
#include "c_r_cp.h"
#include "c_g.h"
#include "c_0_dg_c.h"
#include "c_epslon.h"
#include "c_vkman.h"
#include "c_mdi.h"
#include "blopt8a.h"
! Options for flux gradient formulation, FLUX_GRAD:
!   Locketal2000   = Flux gradients as in Lock et al. (2000)
!   HoltBov1993    = Flux gradients as in Lock et al (2000) but using
!                    coefficients from Holtslag and Boville (1993)
!   LockWhelan2006 = Flux gradients as in Lock and Whelan (2006)

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

      REAL A_PLUME,B_PLUME,A_GRAD_ADJ,A_GA_HB93,A_GA_LW06,MAX_T_GRAD,   &
     &     MAX_SVL_GRAD,SC_CFTOL,CT_RESID,SVL_COUP,DEC_SVL_GRAD,FGF
      PARAMETER (                                                       &
     & A_PLUME=0.2                                                      &
     &,B_PLUME=3.26                                                     &
     &,A_GRAD_ADJ=3.26                                                  &
     &,A_GA_HB93=7.2                                                    &
     &,A_GA_LW06=10.0                                                   &
     &,MAX_T_GRAD=1.0E-3                                                &
!!     &,MAX_Q_GRAD=1.0E-5   ! Or 1.E-6????
     &,MAX_SVL_GRAD=1.0E-3                                              &
                            ! maximum SVL gradient in a mixed layer
     &,DEC_SVL_GRAD=1.0E-3                                              &
                            ! SVL gradient required for weak decoupling
     &,SC_CFTOL=0.1                                                     &
                            ! CF required for a Sc layer to be diagnosed
     &,CT_RESID=200.                                                    &
                            ! Parcel cloud-top residence time (in s)
     &,SVL_COUP=0.5                                                     &
                            ! Parameter controlling positioning of
!                           ! surface-driven entrainment
     &,FGF=0.0)             ! Adiabatic gradient factor for ice
!
!  Define local storage.
!
!  (a) Workspace.
!
      LOGICAL, DIMENSION(row_length,rows) ::                            &
     & CLOUD_BASE                                                       &
                          ! Flag set when cloud base is reached.
     &,DSC_SAVE           ! Copy of DSC needed to indicate
!                         !   decoupling diagnosed in EXCF_NL

      REAL, DIMENSION(row_length,rows,BL_LEVELS) ::                     &
     & QS                                                               &
                          ! Saturated sp humidity at pressure
!                         !   and temperature of sucessive levels.
     &,QCL_IC_BOT                                                       &
                          !
     &,QCF_IC_BOT                                                       &
                          ! In-cloud liquid and frozen water contents
     &,QCL_IC_TOP                                                       &
                          !   at the bottom and top of the model layer
     &,QCF_IC_TOP                                                       &
                          !
     &,CFL                                                              &
                          ! Liquid cloud fraction.
     &,CFF                                                              &
                          ! Frozen cloud fraction.
     &,DQCLDZ                                                           &
                          ! Vertical gradient of in-cloud liquid cloud
!                         !   water in a well-mixed layer.
     &,DQCFDZ                                                           &
                          ! Vertical gradient of in-cloud frozen cloud
!                         !   water in a well-mixed layer.
     &,SLS_INC                                                          &
                          ! SL and QW increments due to large-scale
     &,QLS_INC                                                          &
                          !    vertical advection (K s^-1, s^-1)
     &,DF_OVER_CP                                                       &
                          ! Radiative flux change over layer / c_P
     &,DFLW_OVER_CP                                                     &
                          ! LW radiative flux change over layer / c_P
     &,DFSW_OVER_CP                                                     &
                          ! SW radiative flux change over layer / c_P
     &,SVL                                                              &
                          ! Liquid/frozen water virtual temperature / CP
     &,SL                                                               &
                          ! TL + G*Z/CP (K)
     &,Z_TOP              ! Z_TOP(*,K) is the height of
!                         ! level k+1/2 above the surface.

       REAL, DIMENSION(row_length,rows,2:BL_LEVELS) ::                  &
     & CHI_S                                                            &
                          ! Mixing fraction of just saturate mixture.
     &,DB_CLD                                                           & 
                          ! In-cloud buoyancy jump across layer interf'e
     &,DB                                                               &
                          ! Buoyancy jump across layer interface.
     &,DB_GA_DRY                                                        &
                          ! Out-of-cloud (DRY) and in-cloud buoyancy
     &,DB_NOGA_DRY                                                      &
                          !   jumps used in flux integral calculation.
     &,DB_GA_CLD                                                        &
                          !   GA terms include gradient adjustment
     &,DB_NOGA_CLD        !   arising from non-gradient fluxes. (m/s2)

!-----------------------------------------------------------------------
! The following fluxes, flux changes are in units of rho*Km/s for heat
! and rho*m/s for humidity
!-----------------------------------------------------------------------
      REAL                                                              &
     & DFMIC (row_length,rows,BL_LEVELS,2)                              &
                                           ! Flux changes from microphys
     &,DFSUBS(row_length,rows,BL_LEVELS,2)                              &
                                           !   and subsidence
     &,FRAD (row_length,rows,BL_LEVELS+1)                               &
                                           ! Fluxes from net radiation,
     &,FRAD_LW(row_length,rows,BL_LEVELS+1)                             &
                                           !   LW,
     &,FRAD_SW(row_length,rows,BL_LEVELS+1)                             &
                                           !   SW,
     &,FMIC (row_length,rows,BL_LEVELS+1,2)                             &
                                           !   microphys and subsidence;
     &,FSUBS(row_length,rows,BL_LEVELS+1,2)!   for T, Q separately.

      REAL, DIMENSION(row_length,rows) ::                               &
     & BFLUX_SURF                                                       &
                          ! Buoyancy flux at the surface.
     &,BFLUX_SURF_SAT                                                   &
                          ! Saturated-air surface buoyancy flux
     &,DB_TOP                                                           &
                          ! Buoyancy jump at the top of the BL
     &,DB_DSCT                                                          &
                          ! Buoyancy jump at the DSC layer top
     &,DF_TOP_OVER_CP                                                   &
                          ! Radiative flux change at cloud top / c_P
     &,DF_DSCT_OVER_CP                                                  &
                          ! Radiative flux change at DSC top / CP
     &,SVL_PLUME                                                        &
                          ! SVL, SL and QW for a plume rising without
     &,SL_PLUME                                                         &
                          !   dilution from level 1.
     &,QW_PLUME                                                         &
                          !
     &,ENV_SVL_KM1                                                      &
                          ! Density potential temperature layer K-1
     &,BT_TOP                                                           &
                          ! Buoyancy parameter at the top of the SML
     &,BT_DSCT                                                          &
                          !   and DSC
     &,BTT_TOP                                                          &
                          ! In-cloud buoyancy param at the top of the BL
     &,BTT_DSCT                                                         &
                          !   and DSC
     &,BTC_TOP                                                          &
                          ! Cloud fraction weighted buoyancy parameter
     &,BTC_DSCT                                                         &
                          !   at the top of the SML and DSC
     &,DB_TOP_CLD                                                       &
                          ! In-cloud buoyancy jump at the top of the BL
     &,DB_DSCT_CLD                                                      &
                          !   and DSC
     &,CLD_FACTOR                                                       &
                          ! Fraction of grid box potentially giving
     &,CLD_FACTOR_DSC                                                   &
                          !   evaporative entrainment, for SML and DSC
     &,CHI_S_TOP                                                        &
                          ! Mixing fraction of just saturated mixture
     &,CHI_S_DSCT                                                       &
                          !   at top of the SML and DSC layer
     &,ZETA_S                                                           &
                          ! Non-cloudy fraction of mixing layer for
!                         !   surface forced entrainment term.
     &,ZETA_R                                                           &
                          ! Non-cloudy fraction of mixing layer for
!                         !   cloud-top radiative cooling entrainment
     &,ZETA_R_DSC                                                       &
                          !   term in SML and DSC layers
     &,ZC                                                               &
                          ! Cloud depth (not cloud fraction weighted).
     &,ZC_DSC                                                           &
                          !   for SML and DSC layer (m)
     &,Z_CLD                                                            &
                          ! Cloud fraction weighted depth of cloud.
     &,Z_CLD_DSC                                                        &
                          !   for SML and DSC layers (m)
     &,DSCDEPTH                                                         &
                          ! Depth of DSC layer (m)
     &,D_SIEMS                                                          &
                          ! Siems (1990) et al. cloud-top entrainment
     &,D_SIEMS_DSC                                                      &
                          ! instability parm for SML and DSC inversions
     &,TV1_SD                                                           &
                          ! Standard Deviation of level 1 Tv
     &,FT_NT_ZH                                                         &
                          ! FT_NT at ZH
     &,FT_NT_ZHSC                                                       &
                          ! FT_NT at ZHSC
     &,FQ_NT_ZH                                                         &
                          ! FQ_NT at ZH
     &,FQ_NT_ZHSC                                                       &
                          ! FQ_NT at ZHSC
     &,DF_INV_SML                                                       &
                          ! Radiative flux divergences
     &,DF_INV_DSC                                                       &
                          !   over inversion grid-level
     &,CF_SML                                                           &
                          ! cloud fraction of SML
     &,CF_DSC                                                           &
                          ! cloud fraction of DSC layer
     &,Z_CF_BASE                                                        &
                          ! cloud base height from cloud scheme
     &,Z_CTOP                                                           &
                          ! cloud top height
     &,DQW_SML                                                          &
                          ! QW and SL changes across SML disc inv
     &,DSL_SML                                                          &
                          !
     &,DQW_DSC                                                          &
                          ! QW and SL changes across DSC disc inv
     &,DSL_DSC                                                          &
                          !
     &,RHOKH_SURF_ENT                                                   &
                          ! SML surf-driven entrainment KH
     &,RHOKH_TOP_ENT                                                    &
                          ! SML top-driven entrainment KH
     &,RHOKH_DSCT_ENT                                                   &
                          ! DSC top-driven entrainment KH
     &,ZDSC_BASE                                                        &
                          ! Height of base of K_top in DSC
     &,WE_PARM                                                          &
                          ! Parametrised entrainment rates (m/s)
     &,WE_DSC_PARM                                                      &
                          !   for surf and DSC layers
     &,WE_RHO                                                           &
                          ! rho*entrainment rate
     &,WE_RHO_DSC                                                       &
                          ! rho*entrainment rate for DSC
     &,W_LS                                                             &
                          ! large-scale (subs) velocity
     &,W_LS_DSC                                                         &
                          !   at subgrid inversion heights
     &,ZH_NP1                                                           &
                          ! estimate of ZH at end of timestep
     &,ZHSC_NP1                                                         &
                          ! estimate of ZHSC at end of timestep
     &,ZH_FRAC                                                          &
                          ! (ZH-ZHALF)/DZ
     &,ZHSC_FRAC                                                        &
                          ! (ZHSC-ZHALF)/DZ
     &,ZRZI                                                             &
                          ! (z-z_base)/(z_i-z_base)
     &,ZRZI_DSC                                                         &
                          ! (z-z_base)/(z_i-z_base)
     &,T_FRAC                                                           &
                          ! Fraction of timestep inversion is above
     &,T_FRAC_DSC         !   entr.t flux-level for SML and DSC layers

      REAL                                                              &
     & SCMOUT(row_length,rows)   ! For SCM diagnostics

      INTEGER, DIMENSION(row_length,rows) ::                            &
     & K_CLOUD_TOP                                                      &
                          ! Level number of top of b.l. cloud.
     &,K_CLOUD_DSCT                                                     &
                          ! Level number of top of dec. cloud.
     &,NTML_SAVE                                                        &
                          ! Copy of NTML
     &,NTML_PREV                                                        &
                          ! NTML from previous timestep
     &,K_PLUME                                                          &
                          ! Start grid-level for surface-driven plume
     &,K_CBASE            ! grid-level above cloud-base

! NEC vectorization
      INTEGER, DIMENSION(row_length,rows) ::                            &
     & K_LEVEL                                                          &
                          ! array to store level selection
     &,K_CFF              ! level counter for CFF

!  (b) Scalars.
!
      REAL                                                              &
     & VIRT_FACTOR                                                      &
                         ! Temporary in calculation of buoyancy
!                        ! parameters.
     &,DSLDZ                                                            &
                         ! Vertical gradient of SL in a well-mixed
!                        ! layer.
     &,DQW                                                              &
                         ! Total water content change across layer
!                        ! interface.
     &,DSL                                                              &
                         ! Liquid/ice static energy change across
!                        ! layer interface.
     &,DSL_GA                                                           &
                         ! As DSL but inc gradient adjustment
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
     &,R_D_ETA                                                          &
                         ! 1/(eta(k+1)-eta(k))
     &,SVL_LAPSE                                                        &
                      ! Lapse rate of SVL above inversion (K/m)
     &,SL_LAPSE                                                         &
                      ! Lapse rate of SL above inversion (K/m)
     &,QW_LAPSE                                                         &
                      ! Lapse rate of QW above inversion (kg/kg/m)
     &,SVL_LAPSE_BASE                                                   &
                      ! Lapse rate of SVL above inversion (K/m)
     &,SVL_TOP                                                          &
                      ! s_VL at half level above inversion (K)
     &,DSVL_TOP                                                         &
                      ! s_VL jump across inversion grid layer (K)
     &,TOTHF_EFL                                                        &
                      ! total heat flux at entrainment flux grid-level
     &,TOTQF_EFL                                                        &
                      ! Total QW flux at entrainment flux grid-level
     &,ML_TEND                                                          &
                      ! mixed layer tendency (d/dt)
     &,FA_TEND                                                          &
                      ! free atmospheric tendency (d/dt)
     &,INV_TEND                                                         &
                      ! limit on inversion grid-level tendency (d/dt)
     &,DFLW_INV                                                         &
                      ! temporary in LW rad divergence calculation
     &,DFSW_INV                                                         &
                      ! temporary in SW rad divergence calculation
     &,DZ_DISC_MIN                                                      &
                      ! smallest allowed DZ_DISC
     &,DB_DISC                                                          &
                      ! Temporary disc inversion buoyancy jump
     &,W_S_ENT                                                          &
                      ! numerical (subsidence) entrainment rate
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
                      ! scaling velocity for momentum
     &,W_H                                                              &
                      ! scaling velocity for heat
     &,WSTAR3                                                           &
                      ! cube of convective velocity scale
     &,W_S_CUBED                                                        &
                      ! convective velocity scale
     &,Z_CBASE                                                          &
                      ! cloud base height (m)
     &,ZDSC_CBASE                                                       &
                      ! DSC cloud base height (m)
     &,CF_FOR_WB                                                        &
                      ! CF for use in wb calculation for decoupling
     &,DFSW_TOP                                                         &
                      ! SW radiative flux change assoc with cloud-top
     &,WB_TEST                                                          &
                      ! test wb (m2/s-3)
     &,RATIO                                                            &
                      ! temporary ratio
     &,C_WS                                                             &
                      ! Empirical constant multiplying Wstar
     &,PR_NEUT        ! Neutral Prandtl number
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
     &,KP,KM                                                            &
                   !
     &,K_RAD_LIM   ! limit on levels within which to search for
!                  !   the max LW radiative cooling

      LOGICAL                                                           &
     & MOISTEN
                   ! indicator of whether inversion grid-level should
!                  !   moisten this timestep (or dry)

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('KMKHZ   ',3)
      END IF

!-----------------------------------------------------------------------
! Index to subroutine KMKHZ8C
!
! 1. Set up local variables, etc
! 2. Look for decoupled cloudy mixed-layer above SML top
! 3. Diagnose a discontinuous inversion structure.
! 4. Calculate the within-layer vertical gradients of cloud liquid
!      and frozen water
! 5. Calculate uniform mixed-layer cloud fractions and cloud depths
! 6. Calculate buoyancy flux factor used in the diagnosis of decoupling
! 7. Calculate inputs for the top of b.l. entrainment parametrization
! 8. Calculate the radiative flux change across cloud top
! 9. Calculate the non-turbulent fluxes at the layer boundaries.
! 10.Call subroutine EXCF_NL
!      - calculates parametrized entrainment rate, K profiles and
!        non-gradient flux/stress functions
! 11.Calculate "explicit" entrainment fluxes of SL and QW.
!
!-----------------------------------------------------------------------
! 1. Set up local variables, etc
!-----------------------------------------------------------------------
! 1.1 Calculate Z_TOP (top of levels) and NTML from previous timestep
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
! 1.2 Calculate SVL: conserved buoyanvy-like variable
!-----------------------------------------------------------------------
      DO K = 1,BL_LEVELS
        do j=1,rows
        do i=1,row_length
         SL(I,j,K)  = TL(I,j,K) + GRCP * Z_FULL(I,j,K)
         SVL(I,j,K) = SL(I,j,K) * ( 1.0 + C_VIRTUAL*QW(I,j,K) )
        END DO
        END DO
      END DO
!
      DO  K = 1,BL_LEVELS
! DEPENDS ON: qsat_mix
        CALL QSAT_mix(QS(1,1,K),T(1,1,K),P(1,1,K),                      &
     &                row_length*rows,Lq_mix_bl)
      END DO
!--------------------------------------------------------------------
! 1.3 Integrate non-turbulent increments to give flux profiles:
!     FT_NT, FQ_NT  are the flux profiles from non-turbulent processes
!                  (consisting of radiative FRAD, subsidence FSUBS and
!                   microphysical FMIC fluxes)
!--------------------------------------------------------------------
! For heat, units of rho * Km/s
! For humidity, units of rho * m/s
!----------------------------------
      DO K = 1,BL_LEVELS
      do j=1,rows
      do i=1,row_length

        RHO_DZ = RHO_DRY_TQ(I,j,K) * DZL(I,j,K)

        DFLW_OVER_CP(I,j,K) = - RAD_HR(I,j,K,1) * RHO_DZ
        DFSW_OVER_CP(I,j,K) = - RAD_HR(I,j,K,2) * RHO_DZ
        DF_OVER_CP(I,j,K)   = DFLW_OVER_CP(I,j,K) + DFSW_OVER_CP(I,j,K)

        DFMIC(I,j,K,1)  = - MICRO_TENDS(I,j,K,1) * RHO_DZ
        DFMIC(I,j,K,2)  = - MICRO_TENDS(I,j,K,2) * RHO_DZ

        IF ( ETADOT(I,J,K)  <=  0.0 ) THEN
!         ! only needed in subsidence regions
          KP = MIN( BL_LEVELS, K+1 )
          KM = KP-1
          R_D_ETA = 1.0 /( eta_theta_levels(KP) - eta_theta_levels(KM) )
          SLS_INC(I,J,K) = - ETADOT(I,J,K) * R_D_ETA                    &
     &                                    * ( SL(I,J,KP) - SL(I,J,KM) )
          QLS_INC(I,J,K) = - ETADOT(I,J,K) * R_D_ETA                    &
     &                                    * ( QW(I,J,KP) - QW(I,J,KM) )
        ELSE
          SLS_INC(I,J,K) = 0.0
          QLS_INC(I,J,K) = 0.0
        END IF

        DFSUBS(I,j,K,1) = - SLS_INC(I,j,K)       * RHO_DZ
        DFSUBS(I,j,K,2) = - QLS_INC(I,j,K)       * RHO_DZ

      END DO
      END DO
      END DO
      do j=1,rows
      do i=1,row_length
!       ! Non-turbulent fluxes are defined relative to the surface
!       ! so set them to zero at the surface
        FRAD(I,j,1)  = 0.0
        FRAD_LW(I,j,1) = 0.0
        FRAD_SW(I,j,1) = 0.0
        FSUBS(I,j,1,1) = 0.0 ! for heat
        FSUBS(I,j,1,2) = 0.0 ! for humidity
        FMIC(I,j,1,1) = 0.0  ! for heat
        FMIC(I,j,1,2) = 0.0  ! for humidity
      END DO
      END DO
      DO K = 2,BL_LEVELS+1
        do j=1,rows
        do i=1,row_length
          FRAD(I,j,K)   = FRAD(I,j,K-1)    + DF_OVER_CP(I,j,K-1)
          FRAD_LW(I,j,K)= FRAD_LW(I,j,K-1) + DFLW_OVER_CP(I,j,K-1)
          FRAD_SW(I,j,K)= FRAD_SW(I,j,K-1) + DFSW_OVER_CP(I,j,K-1)
          FSUBS(I,j,K,1)= FSUBS(I,j,K-1,1) + DFSUBS(I,j,K-1,1)
          FSUBS(I,j,K,2)= FSUBS(I,j,K-1,2) + DFSUBS(I,j,K-1,2)
          FMIC(I,j,K,1) = FMIC(I,j,K-1,1)  + DFMIC(I,j,K-1,1)
          FMIC(I,j,K,2) = FMIC(I,j,K-1,2)  + DFMIC(I,j,K-1,2)
        END DO
        END DO
      END DO
      DO K = 1,BL_LEVELS+1
        do j=1,rows
        do i=1,row_length
          FT_NT(I,j,K) = FRAD(I,j,K) + FMIC(I,j,K,1) + FSUBS(I,j,K,1)
          FQ_NT(I,j,K) =               FMIC(I,j,K,2) + FSUBS(I,j,K,2)
        END DO
        END DO
      END DO
!-----------------------------------------------------------------------
! 1.4 Set UNSTABLE flag and find first level above surface layer
!-----------------------------------------------------------------------
      DO J=1,ROWS
      DO I=1,ROW_LENGTH
        UNSTABLE(i,j) = (FB_SURF(i,j) >  0.0)
        K_PLUME(i,j)  = -1
      END DO
      END DO
!     !------------------------------------------------------------
!     ! Find grid-level above top of surface layer, taken
!     ! to be at a height, z_surf, given by:
!     !       Z_SURF = 0.1*ZH_PREV
!     ! Use ZH_prev since that will have determined the shape
!     ! of the time-level n profiles.
!     !------------------------------------------------------------
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
! 2.  Look for decoupled cloudy mixed-layer above SML top
!-----------------------------------------------------------------------
! 2.1  (IF NOT CUMULUS: starting from level 3 and below 2.5km):
!      find cloud-base above SML inversion, ie. above NTML+1,
!      then cloud-top (ie. CF < SC_CFTOL)
!      and finally check that cloud is well-mixed.
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
      DO K=3,BL_LEVELS-1
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
! 2.2 If the layer to ZHPAR is a cumulus layer capped by cloud and
!       an inversion, declare this layer a decoupled cloud layer and
!       set ZHSC and NTDSC accordingly.
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
! 2.3 Calculate the radiative flux changes across cloud top for the
!      stratocumulus layer and thence a first guess for the top-down
!      mixing depth of this layer, DSCDEPTH.
!-----------------------------------------------------------------------
!     Initialise variables
!------------------------------
      do j=1,rows
      do i=1,row_length
        K_CLOUD_DSCT(I,j) = 0
        DF_DSCT_OVER_CP(I,j) = 0.0
        DF_INV_DSC(I,j) = 0.0
      END DO
      END DO
!
      DO K = 1,BL_LEVELS
        do j=1,rows
        do i=1,row_length
!         !-------------------------------------------------------------
!         ! Find the layer with the greatest LW radiative flux jump.
!         ! and assume that this marks the top of the DSC layer.
!         ! Necessary as radiation is not usually called every timestep.
!         !-------------------------------------------------------------
!         ! Limit the search to above the SML.
          K_RAD_LIM = NTML(I,j)+2

          IF ( DSC(I,j) .AND. K  >=  K_RAD_LIM .AND. K  <=  NTDSC(I,j)+2&
     &        .AND. DFLW_OVER_CP(I,j,K)  >   DF_DSCT_OVER_CP(I,j) ) THEN
            K_CLOUD_DSCT(I,j) = K
!           ! Set K_CLOUD_DSCT to the level below if its DF is greater
!           ! than half the maximum.  DF in level K_CLOUD_DSCT+1 is then
!           ! included as DF_INV_DSC below.
            IF (DFLW_OVER_CP(I,j,K-1)  >   0.5*DFLW_OVER_CP(I,j,K))     &
     &         K_CLOUD_DSCT(I,j) = K-1
            DF_DSCT_OVER_CP(I,j) = DFLW_OVER_CP(I,j,K)
          END IF

        END DO
        END DO
      END DO
!     !-----------------------------------------------------------------
!     !  Find bottom grid-level (K_LEVEL) for cloud-top radiative flux
!     !  divergence: higher of base of LW radiatively cooled layer,
!     !  ZH and 0.5*ZHSC, since cooling must be in upper part of layer
!     !  in order to generate turbulence.
!     !-----------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
        K_LEVEL(I,j) = K_CLOUD_DSCT(I,j)
        IF ( K_CLOUD_DSCT(I,j)  >   1 ) THEN
            K_RAD_LIM = NTML(I,j)+1
            K=K_CLOUD_DSCT(I,j)-1
            KL=MAX(1,K)  ! only to avoid out-of-bounds compiler warning
            DO WHILE ( K  >   K_RAD_LIM                                 &
     &                .AND. DFLW_OVER_CP(I,j,KL)  >   0.0               &
     &                .AND. Z_FULL(I,j,KL)  >   0.5*ZHSC(I,j) )
              K_LEVEL(I,j) = K
              K = K-1
              KL=MAX(1,K)
            END DO
        END IF
      END DO
      END DO
!     !-----------------------------------------------------------------
!     ! Calculate LW and SW flux divergences and combine into
!     ! cloud-top turbulence forcing.
!     ! Need to account for radiative divergence in cloud in inversion
!     ! grid-level, DF_INV_DSC. Assume DF_OVER_CP(K_cloud_dsct+2) is
!     ! representative of clear-air rad divergence and so subtract this
!     ! `clear-air' part from the grid-level divergence.
!     !-----------------------------------------------------------------
      do j=1,rows
      do i=1,row_length

        IF ( K_CLOUD_DSCT(I,j)  >   0 ) THEN
          DFLW_INV = 0.0
          DFSW_INV = 0.0
          IF ( K_CLOUD_DSCT(I,j) <  BL_LEVELS ) THEN
            K = K_CLOUD_DSCT(I,j)+1
            IF ( K  <   BL_LEVELS ) THEN
              DFLW_INV = DFLW_OVER_CP(I,j,K)                            &
     &                   - DFLW_OVER_CP(I,j,K+1)                        &
     &                          * DZL(I,j,K)/DZL(I,j,K+1)
              DFSW_INV = DFSW_OVER_CP(I,j,K)                            &
     &                   - DFSW_OVER_CP(I,j,K+1)                        &
     &                          * DZL(I,j,K)/DZL(I,j,K+1)
            ELSE
              DFLW_INV = DFLW_OVER_CP(I,j,K)
              DFSW_INV = DFSW_OVER_CP(I,j,K)
            END IF
            DFLW_INV = MAX( DFLW_INV, 0.0 )
            DFSW_INV = MIN( DFSW_INV, 0.0 )
          END IF
          DF_INV_DSC(I,j) = DFLW_INV + DFSW_INV

          DF_DSCT_OVER_CP(I,j) = FRAD_LW(I,j,K_CLOUD_DSCT(I,j)+1)       &
     &                         - FRAD_LW(I,j,K_LEVEL(I,j))              &
     &                         + DFLW_INV

          DFSW_TOP = FRAD_SW(I,j,K_CLOUD_DSCT(I,j)+1)                   &
     &             - FRAD_SW(I,j,K_LEVEL(I,j))                          &
     &             + DFSW_INV

!         !-----------------------------------------------------------
!         ! Combine SW and LW cloud-top divergences into a net
!         ! divergence by estimating SW flux divergence at a given
!         ! LW divergence = DF_SW * (1-exp{-A*kappa_sw/kappa_lw})
!         ! Empirically (from LEM data) a reasonable fit is found
!         ! with A small and (1-exp{-A*kappa_sw/kappa_lw}) = 0.35
!         !-----------------------------------------------------------
          DF_DSCT_OVER_CP(I,j) = MAX( 0.0,                              &
     &                  DF_DSCT_OVER_CP(I,j) + 0.35 * DFSW_TOP )
        END IF
      END DO
      END DO
!-----------------------------------------------------------------------
! 2.4 Set NBDSC, the bottom level of the DSC layer.
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
! 2.5 Tidy up variables associated with decoupled layer
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
           K_CLOUD_DSCT(I,j) = 0
           DF_INV_DSC(I,j)   = 0.0
           DSCDEPTH(I,j) = 0.0
         END IF
       END IF
       ELSE  ! NTDSC EQ 0, just to make sure!
        DSCDEPTH(I,j)=0.0
        DSC(I,j)=.FALSE.
        ZHSC(I,j)=0.0
        DF_DSCT_OVER_CP(I,j) = 0.0
        K_CLOUD_DSCT(I,j) = 0
        DF_INV_DSC(I,j)   = 0.0
       END IF
      END DO
      END DO
!----------------------------------------------------------------------
!2.6 If decoupled cloud-layer found test to see if it is, in fact,
!  only weakly decoupled from the surface mixed-layer:
!  if SVL difference between NTML and NTDSC is less than SVL_COUP (in K)
!  then assume there is still some coupling.  This will mean that
!  the surface-driven entrainment term will be applied at ZHSC, no
!  subgrid inversion or entrainment will be calculated for ZH and
!  ZHSC will be the length scale used in the entrainment inputs.
!  Note that for CUMULUS "surface-driven entrainment" will be done
!  by the convection scheme.
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
! 3. Diagnose a discontinuous inversion structure:
!    - to this point in the code, ZH and ZHSC mark the half-level at
!      the base of the inversion
!    - now they will be interpolated into the level above assuming
!      SVL(NTML+1) is a volume average over a subgrid discontinuous
!      inversion structure
!    - discontinuous jumps of SL and QW (and thence buoyancy) can be
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

        KP2 = MIN(K+2, BL_LEVELS)
        DSL_SML(I,j) = SL(I,j,KP2) - SL(I,j,K)
        DQW_SML(I,j) = QW(I,j,KP2) - QW(I,j,K)

        IF ( .NOT.CUMULUS(I,j) .AND. .NOT.COUPLED(I,j) .AND.            &
     &                K  >   1 .AND. K  <=  BL_LEVELS-2 ) THEN
          IF ( SVL(I,j,K+2)  >   SVL(I,j,K+1)                           &
     &                  .AND. SVL(I,j,K+1)  >   SVL(I,j,K) ) THEN

              IF ( K  <=  BL_LEVELS-3 ) THEN
! need to test for K+1 to K+2 gradient in case profile is concave
! (would mess up the inversion diagnosis so best just to ignore lapse)
                SVL_LAPSE = MAX(0.0,                                    &
     &                ( SVL(I,j,K+3) - SVL(I,j,K+2) ) * RDZ(I,j,K+3)  )
                IF ( SVL_LAPSE  >                                       &
     &                ( SVL(I,j,K+2) - SVL(I,j,K+1) ) * RDZ(I,j,K+2) )  &
     &                SVL_LAPSE = 0.0
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

!..Calculate discontinuous jumps of SL and QW:
!
!           ! Take double grid-level jumps as a starting point
            DSL_SML(I,j) = SL(I,j,K+2) - SL(I,j,K)
            DQW_SML(I,j) = QW(I,j,K+2) - QW(I,j,K)

!           ! Allow for lapse rate above inversion, if known
            IF ( K  <=  BL_LEVELS-3 ) THEN
              SL_LAPSE = MAX( 0.0,                                      &
     &           ( SL(I,j,K+3) - SL(I,j,K+2) )*RDZ(I,j,K+3) )
              QW_LAPSE = MIN( 0.0,                                      &
     &           ( QW(I,j,K+3) - QW(I,j,K+2) )*RDZ(I,j,K+3) )
              DZ_DISC = Z_FULL(I,j,K+2) - ZH(I,j)
!             ! Only reduce inversion jumps by at most half
              DSL_SML(I,j) = DSL_SML(I,j) -                             &
     &             MIN( 0.5*DSL_SML(I,j), SL_LAPSE*DZ_DISC )
              IF ( DQW_SML(I,j)  <   0.0 ) THEN
                DQW_SML(I,j) = DQW_SML(I,j) -                           &
     &             MAX( 0.5*DQW_SML(I,j), QW_LAPSE*DZ_DISC )
              END IF
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
          KP2 = MIN(K+2, BL_LEVELS)
          DSL_DSC(I,j) = SL(I,j,KP2) - SL(I,j,K)
          DQW_DSC(I,j) = QW(I,j,KP2) - QW(I,j,K)
          IF ( K  <=  BL_LEVELS-2 ) THEN
          IF ( SVL(I,j,K+2)  >   SVL(I,j,K+1)                           &
     &             .AND. SVL(I,j,K+1)  >   SVL(I,j,K) ) THEN

              IF ( K  <=  BL_LEVELS-3 ) THEN
! need to test for K+1 to K+2 gradient in case profile is concave
! (would mess up the inversion diagnosis so best just to ignore)
                SVL_LAPSE = MAX(0.0,                                    &
     &                ( SVL(I,j,K+3) - SVL(I,j,K+2) )*RDZ(I,j,K+3) )
                IF ( SVL_LAPSE  >                                       &
     &                ( SVL(I,j,K+2) - SVL(I,j,K+1) )*RDZ(I,j,K+2) )    &
     &                SVL_LAPSE = 0.0
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

!..Calculate discontinuous jumps of SL and QW:
!
!           ! Take double grid-level jumps as a starting point
            DSL_DSC(I,j) = SL(I,j,K+2) - SL(I,j,K)
            DQW_DSC(I,j) = QW(I,j,K+2) - QW(I,j,K)

!           ! Allow for lapse rate above inversion, if known
            IF ( K  <=  BL_LEVELS-3 ) THEN
              SL_LAPSE = MAX( 0.0,                                      &
     &           ( SL(I,j,K+3) - SL(I,j,K+2) )*RDZ(I,j,K+3) )
              QW_LAPSE = MIN( 0.0,                                      &
     &           ( QW(I,j,K+3) - QW(I,j,K+2) )*RDZ(I,j,K+3) )
              DZ_DISC = Z_FULL(I,j,K+2) - ZHSC(I,j)
!             ! Only reduce inversion jumps by at most half
              DSL_DSC(I,j) = DSL_DSC(I,j) -                             &
     &             MIN( 0.5*DSL_DSC(I,j), SL_LAPSE*DZ_DISC )
              IF ( DQW_DSC(I,j)  <   0.0 ) THEN
!               ! Only allow for QW lapse rate if both it and the
!               ! grid-level jump negative (expected sign)
                DQW_DSC(I,j) = DQW_DSC(I,j) -                           &
     &             MAX( 0.5*DQW_DSC(I,j), QW_LAPSE*DZ_DISC )
              END IF
            END IF

          END IF  ! SVL increasing
          END IF ! test on K LT BL_LEVELS-2
         END IF ! test on DSC
      END DO
      END DO
!-----------------------------------------------------------------------
! 4.  Calculate the within-layer vertical gradients of cloud liquid
!      and frozen water
! 4.1 First for level 1
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
!
        VIRT_FACTOR = 1.0 + C_VIRTUAL*Q(I,j,1) - QCL(I,j,1) - QCF(I,j,1)
        GRAD_T_ADJ(I,j) = 0.0
        GRAD_Q_ADJ(I,j) = 0.0
        IF ( UNSTABLE(I,j) ) THEN
!         ! Here this is an estimate of the gradient adjustment applied
!         ! the previous timestep (assumes T1_SD has not changed much,
!         ! which in turn assumes the surface fluxes have not)
          IF (FLUX_GRAD  ==  Locketal2000) THEN
            GRAD_T_ADJ(I,j) = MIN( MAX_T_GRAD ,                         &
     &                       A_GRAD_ADJ * T1_SD(I,j) / ZH_PREV(I,j) )
            GRAD_Q_ADJ(I,j) = 0.0
          ELSEIF (FLUX_GRAD  ==  HoltBov1993) THEN
!           ! Use constants from Holtslag and Boville (1993)
!           ! Conv limit GAMMA_TH = 10 *FTL1/(wstar*zh)
!           ! Neut limit GAMMA_TH = 7.2*wstar*FTL1/(ustar^2*zh)
            WSTAR3 = FB_SURF(I,j) * ZH_PREV(I,j)
            C_WS = 0.6
            W_M =( V_S(I,j)**3.0 + C_WS*WSTAR3 )**(1.0/3.0)

            GRAD_T_ADJ(I,j) = A_GA_HB93*(WSTAR3**(1.0/3.0))*FTL(I,j,1)  &
     &                        / ( RHOSTAR_GB(I,j)*W_M*W_M*ZH_PREV(I,j) )
!            GRAD_Q_ADJ(I,j) = A_GA_HB93*(WSTAR3**(1.0/3.0))*FQW(I,j,1)
!     &                       / ( RHOSTAR_GB(I,j)*W_M*W_M*ZH_PREV(I,j) )
            GRAD_Q_ADJ(I,j) = 0.0
          ELSEIF (FLUX_GRAD  ==  LockWhelan2006) THEN
!           ! Use constants from LockWhelan2006
!           ! Conv limit GAMMA_TH = 10 *FTL1/(wstar*zh)
!           ! Neut limit GAMMA_TH = 7.5*FTL1/(ustar*zh)
            WSTAR3  = FB_SURF(I,j) * ZH_PREV(I,j)
            C_WS    = 0.42   !  = 0.75^3
            PR_NEUT = 0.75
            W_H = ( ( V_S(I,j)**3.0+C_WS*WSTAR3 )**(1.0/3.0) )/ PR_NEUT

            GRAD_T_ADJ(I,j) = A_GA_LW06 * FTL(I,j,1)                    &
     &                         / ( RHOSTAR_GB(I,j)*W_H*ZH_PREV(I,j) )
            GRAD_Q_ADJ(I,j) = A_GA_LW06 * FQW(I,j,1)                    &
     &                         / ( RHOSTAR_GB(I,j)*W_H*ZH_PREV(I,j) )
          END IF
        END IF  ! test on UNSTABLE

        DSLDZ = -GRCP + GRAD_T_ADJ(I,j)
        DQCLDZ(I,j,1) = -( DSLDZ*DQSDT(I,j,1) +                         &
     &                   G*QS(I,j,1)/(R*T(I,j,1)*VIRT_FACTOR) )         &
     &                  / (1.0 + LCRCP*DQSDT(I,j,1))
        DQCFDZ(I,j,1) = -( DSLDZ*DQSDT(I,j,1) +                         &
     &                   G*QS(I,j,1)/(R*T(I,j,1)*VIRT_FACTOR) ) * FGF   &
     &                  / (1.0 + LSRCP*DQSDT(I,j,1))
!
!       ! limit calculation to greater than a small cloud fraction
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
!-----------------------------------------------------------------------
!! 4.2 Calculate the within-layer vertical gradients of cloud liquid
!!     and frozen water for the current layer
!-----------------------------------------------------------------------
      DO K=2,BL_LEVELS
!
         do j=1,rows
         do i=1,row_length
!
          IF (K  <=  NTML_PREV(I,j)) THEN
            DSLDZ = -GRCP + GRAD_T_ADJ(I,j)
          ELSE
            DSLDZ = -GRCP
          END IF
!
! RNBS correction  28/11/97
          VIRT_FACTOR = 1.0 + C_VIRTUAL*Q(I,j,k) - QCL(I,j,k) -         &
     &                        QCF(I,j,k)

          DQCLDZ(I,j,K) = -( DSLDZ*DQSDT(I,j,K)                         &
     &                   + G*QS(I,j,K)/(R*T(I,j,K)*VIRT_FACTOR) )       &
     &                    / ( 1.0 + LCRCP*DQSDT(I,j,K) )
          DQCFDZ(I,j,K) = -( DSLDZ*DQSDT(I,j,K)                         &
     &                   + G*QS(I,j,K)/(R*T(I,j,K)*VIRT_FACTOR) ) * FGF &
     &                    / ( 1.0 + LSRCP*DQSDT(I,j,K) )
!
!-----------------------------------------------------------------------
! 4.3  Calculate the cloud liquid and frozen water contents at the
!      top and bottom of the current layer
!-----------------------------------------------------------------------
!         ! limit calculation to greater than a small cloud fraction
          IF ( QCL(I,j,k) + QCF(I,j,k)  >   0.0                         &
     &         .and. CF(i,j,k)  >   1.e-3 ) THEN
            CFL(I,j,K) = CF(I,j,K) * QCL(I,j,K) /                       &
     &                   ( QCL(I,j,K) + QCF(I,j,K) )
            CFF(I,j,K) = CF(I,j,K) * QCF(I,j,K) /                       &
     &                   ( QCL(I,j,K) + QCF(I,j,K) )
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
!-----------------------------------------------------------------------
!! 4.4 Calculate the jumps of QW, SL and buoyancy across the layer
!!     interface at level k-1/2.
!-----------------------------------------------------------------------
      DO K=2,BL_LEVELS
        KM1 = K-1
         do j=1,rows
         do i=1,row_length
!
!
          DQW = QW(I,j,K) - QW(I,j,KM1)              ! Used in P243.C2
          DSL = SL(I,j,K) - SL(I,j,KM1)              ! Used in P243.C2
!
          DQCL = CFL(I,j,K)*QCL_IC_BOT(I,j,K) -                         &
     &             CFL(I,j,KM1)*QCL_IC_TOP(I,j,KM1)
          DQCF = CFF(I,j,K)*QCF_IC_BOT(I,j,K) -                         &
     &             CFF(I,j,KM1)*QCF_IC_TOP(I,j,KM1)
!
          DB(I,j,K) = G * ( BTM(I,j,KM1)*DSL + BQM(I,j,KM1)*DQW +       &
     &               (LCRCP*BTM(I,j,KM1) - ETAR*BQM(I,j,KM1)) * DQCL +  &
     &                (LSRCP*BTM(I,j,KM1) - ETAR*BQM(I,j,KM1)) * DQCF )
!
          DB_CLD(I,j,K) =G * ( BTM_CLD(I,j,KM1)*DSL                     &
     &         + BQM_CLD(I,j,KM1)*DQW )
!
          CHI_S(I,j,K) = -QCL_IC_TOP(I,j,KM1) /                         &
     &                  (A_QSM(I,j,KM1)*DQW - A_DQSDTM(I,j,KM1)*DSL)
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

            DB_DISC = G * ( BTM(I,j,K-1)*DSL_SML(I,j)                   &
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

              DB_CLD(I,j,K) = G * ( BTM_CLD(I,j,K-1)*DSL_SML(I,j) +     &
     &                            BQM_CLD(I,j,K-1)*DQW_SML(I,j) )

              CHI_S(I,j,K) = -QCL_IC_TOP(I,j,K) /                       &
     &            ( A_QSM(I,j,K-1)*DQW_SML(I,j)                         &
     &                  - A_DQSDTM(I,j,K-1)*DSL_SML(I,j) )
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
     &         BTM(I,j,K-1)*DSL_DSC(I,j) + BQM(I,j,K-1)*DQW_DSC(I,j) +  &
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

              DB_CLD(I,j,K) = G * ( BTM_CLD(I,j,K-1)*DSL_DSC(I,j)       &
     &                          + BQM_CLD(I,j,K-1)*DQW_DSC(I,j) )

              CHI_S(I,j,K) = -QCL_IC_TOP(I,j,K) /                       &
     &                 ( A_QSM(I,j,K-1)*DQW_DSC(I,j)                    &
     &                         - A_DQSDTM(I,j,K-1)*DSL_DSC(I,j) )
            END IF

          END IF  ! disc inversion diagnosed

      END DO
      END DO
!-----------------------------------------------------------------------
! 5.  Calculate uniform mixed-layer cloud fraction and thence
!        estimate Sc layer cloud depth (not cloud fraction weighted).
!        (If DSC=.FALSE. then NTDSC=0 and ZC_DSC remains equal to 0.)
!-----------------------------------------------------------------------
! First the SML
!---------------
      do j=1,rows
      do i=1,row_length
        CLOUD_BASE(I,j)= .FALSE.
        ZC(I,j)        = 0.0
        K_CBASE(I,j)   = 0
        Z_CF_BASE(I,j) = ZH(I,j)
        Z_CTOP(I,j)    = ZH(I,j)
!       ! Use a single CF for whole mixed-layer (more realistic).
!       ! Include NTML+1 if a subgrid inversion has been diagnosed
        IF ( COUPLED(I,j) .OR. CUMULUS(I,j) .OR. NTML(I,j) == 1 ) THEN
          CF_SML(I,j)=0.0
        ELSE
          K = NTML(I,j)
          CF_SML(I,j) = MAX( CF(I,j,K), CF(I,j,K-1) )
        END IF
      END DO
      END DO
!-----------------------------------------------------------------------
! First find cloud-base as seen by the cloud scheme, at grid-level
! K_LEVEL and height Z_CF_BASE, to use as first guess or lower limit
!-----------------------------------------------------------------------
      DO J=1,ROWS
      DO I=1,ROW_LENGTH
        K_LEVEL(I,J) = NTML(I,J)
        IF ( CF_SML(I,j)  >   SC_CFTOL ) THEN
          K_LEVEL(I,J) = NTML(I,J)-1
          DO WHILE ( CF(I,J,K_LEVEL(I,J))  >   SC_CFTOL                 &
     &               .AND. K_LEVEL(I,J)  >=  2 )
            K_LEVEL(I,J) = K_LEVEL(I,J) - 1
          END DO
        END IF
      END DO
      END DO
      do j=1,rows
      do i=1,row_length
        IF ( CF_SML(I,j)  >   SC_CFTOL ) THEN
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
          IF ( .NOT.CLOUD_BASE(I,j) .AND. K  <=  NTML(I,j) .AND.        &
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
     &                  .AND. K_CFF(I,J)  >   1 )
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
        CLOUD_BASE(I,j) = .FALSE.
        ZC_DSC(I,j) = 0.0
        CF_DSC(I,j) = 0.0
        K_CBASE(I,j) = 0
        Z_CF_BASE(I,j) = ZHSC(I,j)
        Z_CTOP(I,j)    = ZHSC(I,j)

        IF ( DSC(I,j) ) THEN
          K = NTDSC(I,j)
          CF_DSC(I,j) = MAX( CF(I,j,K), CF(I,j,K-1) )
        END IF
      END DO
      END DO
!-------------------------------------------------------------
! Find cloud-base as seen by cloud scheme, Z_CF_BASE,
! to use as first guess or lower limit and find cloud top.
!-------------------------------------------------------------
      DO J=1,ROWS
      DO I=1,ROW_LENGTH
        K_LEVEL(I,J) = NTDSC(I,J)
        IF ( CF_DSC(I,j)  >   SC_CFTOL ) THEN
!         ! assume level NTDSC is cloudy so start from NTDSC-1
          K_LEVEL(I,J) = MAX( 2, NTDSC(I,J) - 1 )
          DO WHILE ( CF(I,J,K_LEVEL(I,J))  >   SC_CFTOL                 &
     &             .AND. K_LEVEL(I,J)  >=  2 )
            K_LEVEL(I,J) = K_LEVEL(I,J) - 1
          END DO
        END IF
      END DO
      END DO

      DO J=1,ROWS
      DO I=1,ROW_LENGTH
        IF ( CF_DSC(I,j)  >   SC_CFTOL ) THEN
          IF ( K_LEVEL(i,j)  ==  1 .AND.                                &
     &         CF(I,j,K_LEVEL(i,j))  >   SC_CFTOL) THEN
            Z_CF_BASE(I,j) = 0.0
          ELSE
            Z_CF_BASE(I,j) = Z_HALF(I,j,K_LEVEL(i,j)+1)
          END IF
          ZC_DSC(I,j) = Z_CTOP(I,j) - Z_CF_BASE(I,j)   ! first guess
        END IF
      END DO
      END DO
!--------------------------------------------------
! Find lowest level within ML with max CF
!--------------------------------------------------
      DO K = BL_LEVELS,1,-1
        do j=1,rows
        do i=1,row_length
          IF ( .NOT.CLOUD_BASE(I,j) .AND. K  <=  NTDSC(I,j) .AND.       &
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
! 6. Calculate buoyancy flux factor used in the diagnosis of decoupling
!-----------------------------------------------------------------------
      DO K = 2,BL_LEVELS
        do j=1,rows
        do i=1,row_length
          DQW = QW(I,j,K) - QW(I,j,K-1)
          DSL = SL(I,j,K) - SL(I,j,K-1)
          IF (K  <=  NTML_PREV(I,j)) THEN
            DSL_GA = DSL - GRAD_T_ADJ(I,j)/RDZ(I,j,K)
            DQW_GA = DQW - GRAD_Q_ADJ(I,j)/RDZ(I,j,K)
          ELSE
            DSL_GA = DSL
            DQW_GA = DQW
          END IF
!         !----------------------------------------------------------
!         ! CF_FOR_WB is uniform `bl' CF for use within cloud layers
!         !----------------------------------------------------------
          CF_FOR_WB = 0.0
          Z_CBASE = ZH(I,j)-ZC(I,j)
          ZDSC_CBASE = ZHSC(I,j)-ZC_DSC(I,j)
          IF ( Z_FULL(I,j,K)  <=  ZH(I,j) .AND.                         &
     &         Z_FULL(I,j,K)  >=  Z_CBASE) CF_FOR_WB = CF_SML(I,j)
          IF ( Z_FULL(I,j,K)  <=  ZHSC(I,j) .AND.                       &
     &         Z_FULL(I,j,K)  >=  ZDSC_CBASE) CF_FOR_WB = CF_DSC(I,j)
!         !----------------------------------------------------------
!         ! WB = -K_SURF*(DB/DZ - gamma_buoy) - K_TOP*DB/DZ
!         ! This is integrated in EXCF_NL, iterating the K profiles.
!         ! Here the relevant integrated DB/DZ factors are calculated
!         !----------------------------------------------------------
          DB_GA_DRY(I,j,K) = - G *                                      &
     &             ( BTM(I,j,K-1)*DSL_GA + BQM(I,j,K-1)*DQW_GA )
          DB_NOGA_DRY(I,j,K)  = - G *                                   &
     &             ( BTM(I,j,K-1)*DSL + BQM(I,j,K-1)*DQW )
          DB_GA_CLD(I,j,K) = - G *                                      &
     &             ( BTM_CLD(I,j,K-1)*DSL_GA + BQM_CLD(I,j,K-1)*DQW_GA )
          DB_NOGA_CLD(I,j,K)  = - G *                                   &
     &             ( BTM_CLD(I,j,K-1)*DSL + BQM_CLD(I,j,K-1)*DQW )
!         !-------------------------------------------------------
!         ! Weight cloud layer factors with cloud fraction
!         !-------------------------------------------------------
          DB_GA_CLD(I,j,K) = DB_GA_DRY(I,j,K)*(1.0-CF_FOR_WB) +         &
     &                          DB_GA_CLD(I,j,K)*CF_FOR_WB
          DB_NOGA_CLD(I,j,K)  = DB_NOGA_DRY(I,j,K)*(1.0-CF_FOR_WB) +    &
     &                          DB_NOGA_CLD(I,j,K)*CF_FOR_WB
        END DO
      END DO
      END DO
!-----------------------------------------------------------------------
! 7. Calculate inputs for the top of b.l. entrainment parametrization
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
! 7.1 Calculate surface buoyancy flux
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length

!       ! use mixed-layer average of buoyancy parameters
        BFLUX_SURF(I,j) = 0.5 * G * (                                   &
     &       (BTM(I,j,1)+BTM(I,j,NTML(I,j)))*FTL(I,j,1) +               &
     &       (BQM(I,j,1)+BQM(I,j,NTML(I,j)))*FQW(I,j,1) )

        IF ( BFLUX_SURF(I,j)  >   0.0 ) THEN
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
! 7.2 Calculation of cloud fraction weighted thickness of
!     cloud in the SML and DSC layer (or to the surface if COUPLED)
!-----------------------------------------------------------------------
      DO K = 1,BL_LEVELS
        do j=1,rows
        do i=1,row_length
          IF ( K  <=  NTML(I,j)+1 ) THEN
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
        END DO
        END DO
      END DO
!-----------------------------------------------------------------------
! 7.3 Calculation of other SML and DSC inputs to entr param.
!     If COUPLED then SML are not used as no "entrainment" is then
!     applied at ZH.
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
!         !------------------------------------------------------
!         ! Calculation of SML inputs.
!         !------------------------------------------------------
          K = NTML(I,j)
          IF (.NOT.COUPLED(I,j) ) THEN
            CHI_S_TOP(I,j) = MAX( 0.0, MIN( CHI_S(I,j,K+1), 1.) )
            KP2=MIN(K+1+SML_DISC_INV(I,j),BL_LEVELS)
            CLD_FACTOR(I,j) = MAX( 0.0 , CF(I,j,K)-CF(I,j,KP2) )
            BT_TOP(I,j)  = G * BTM(I,j,K)
            BTT_TOP(I,j) = G * BTM_CLD(I,j,K)
            BTC_TOP(I,j) = BTT_TOP(I,j)
            DB_TOP(I,j)  = DB(I,j,K+1)
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
!         !---------------------------------------------------
!         ! Calculation of DSC inputs
!         !---------------------------------------------------
          IF (DSC(I,j)) THEN
            K = NTDSC(I,j)
            CHI_S_DSCT(I,j) = MAX( 0.0, MIN( CHI_S(I,j,K+1), 1.) )
            KP2=MIN(K+1+DSC_DISC_INV(I,j),BL_LEVELS)
            CLD_FACTOR_DSC(I,j) = MAX( 0.0 , CF(I,j,K)-CF(I,j,KP2) )
            BT_DSCT(I,j)  = G * BTM(I,j,K)
            BTT_DSCT(I,j) = G * BTM_CLD(I,j,K)
            BTC_DSCT(I,j) = BTT_DSCT(I,j)
            DB_DSCT(I,j)  = DB(I,j,K+1)
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
!-----------------------------------------------------------------------
! 7.4 Next those terms which depend on the presence of buoyancy reversal
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
! 8. Calculate the radiative flux change across cloud top for mixed-
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
        DF_INV_SML(I,j) = 0.0
      END DO
      END DO
!
      DO K = 1,BL_LEVELS
        do j=1,rows
        do i=1,row_length
!         ! restrict search to `close' to ZH
          K_RAD_LIM = NTML(I,j)+1
!         !-------------------------------------------------------------
!         ! Find the layer below K_RAD_LIM with the greatest LW
!         ! radiative flux jump and assume that this is
!         ! the top of the SML.
!         !-------------------------------------------------------------
          IF (DFLW_OVER_CP(I,j,K)  >   DF_TOP_OVER_CP(I,j)              &
     &                .AND. K  <=  K_RAD_LIM ) THEN
            K_CLOUD_TOP(I,j) = K
            IF ( K >  1 ) THEN
!             ! Set K_CLOUD_TOP to the level below if its DF is
!             ! greater than half the maximum.  DF in level
!             ! K_CLOUD_TOP+1 is then included as DF_INV_SML below.
              IF (DFLW_OVER_CP(I,j,K-1)  >   0.5*DFLW_OVER_CP(I,j,K))   &
     &          K_CLOUD_TOP(I,j) = K-1
            END IF
          END IF

        END DO
        END DO
      END DO
!     !-----------------------------------------------------------------
!     !  Find bottom grid-level (K_LEVEL) for cloud-top radiative fux
!     !  divergence: higher of base of LW radiatively cooled layer,
!     !  0.5*ZH, since cooling must be in upper part of layer
!     !-----------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
        K_LEVEL(I,j) = K_CLOUD_TOP(I,j)
        IF ( K_CLOUD_TOP(I,j)  >   1 ) THEN
            K_RAD_LIM = 1
            K=K_CLOUD_TOP(I,j)-1
            KL=MAX(1,K)  ! only to avoid out-of-bounds compiler warning
            DO WHILE ( K  >   K_RAD_LIM                                 &
     &                .AND. DFLW_OVER_CP(I,j,KL)  >   0.0               &
     &                .AND. Z_FULL(I,j,KL)  >   0.5*ZH(I,j) )
              K_LEVEL(I,j) = K
              K = K-1
              KL=MAX(1,K)
            END DO
        END IF
      END DO
      END DO
!     !-----------------------------------------------------------------
!     ! Calculate LW and SW flux divergences and combine into
!     ! cloud-top turbulence forcing.
!     ! Need to account for radiative divergences in cloud in inversion
!     ! grid-level, DF_INV_SML. Assume DF_OVER_CP(K_cloud_top+2) is
!     ! representative of clear-air rad divergence and so subtract this
!     ! `clear-air' part from the grid-level divergence.
!     !-----------------------------------------------------------------
      do j=1,rows
      do i=1,row_length

        IF ( K_CLOUD_TOP(I,j)  >   0 ) THEN
          DFLW_INV = 0.0
          DFSW_INV = 0.0
          IF ( K_CLOUD_TOP(I,j) <  BL_LEVELS ) THEN
            K = K_CLOUD_TOP(I,j)+1
            IF ( K  <   BL_LEVELS ) THEN
              DFLW_INV = DFLW_OVER_CP(I,j,K)                            &
     &                   - DFLW_OVER_CP(I,j,K+1)                        &
     &                          * DZL(I,j,K)/DZL(I,j,K+1)
              DFSW_INV = DFSW_OVER_CP(I,j,K)                            &
     &                   - DFSW_OVER_CP(I,j,K+1)                        &
     &                          * DZL(I,j,K)/DZL(I,j,K+1)
            ELSE
              DFLW_INV = DFLW_OVER_CP(I,j,K)
              DFSW_INV = DFSW_OVER_CP(I,j,K)
            END IF
            DFLW_INV = MAX( DFLW_INV, 0.0 )
            DFSW_INV = MIN( DFSW_INV, 0.0 )
          END IF
          DF_INV_SML(I,j) = DFLW_INV + DFSW_INV

          DF_TOP_OVER_CP(I,j) = FRAD_LW(I,j,K_CLOUD_TOP(I,j)+1)         &
     &                         - FRAD_LW(I,j,K_LEVEL(I,j))              &
     &                         + DFLW_INV

          DFSW_TOP = FRAD_SW(I,j,K_CLOUD_TOP(I,j)+1)                    &
     &             - FRAD_SW(I,j,K_LEVEL(I,j))                          &
     &             + DFSW_INV

!         !-----------------------------------------------------------
!         ! Combine SW and LW cloud-top divergences into a net
!         ! divergence by estimating SW flux divergence at a given
!         ! LW divergence = DF_SW * (1-exp{-A*kappa_sw/kappa_lw})
!         ! Empirically (from LEM data) a reasonable fit is found
!         ! with A small and (1-exp{-A*kappa_sw/kappa_lw}) = 0.35
!         !-----------------------------------------------------------
          DF_TOP_OVER_CP(I,j) = MAX( 0.0,                               &
     &                  DF_TOP_OVER_CP(I,j) + 0.35 * DFSW_TOP )
        END IF
      END DO
      END DO
! ------------------------------------------------------------------
! 9. Calculate the non-turbulent fluxes at the layer boundaries.
!  - the radiative flux at the inversion allows for an estimate
!    of the FA flux divergence within the inversion grid-level
!  - because the radiative time-step is usually longer the radiative
!    cloud-top grid-level (K_CLOUD_TOP) is allowed to differ from
!    the actual one (NTML)
!  - the subsidence flux at the inversion is taken from the
!    flux grid-level below it (assumes the divergence across
!    the inversion is physically above the BL)
!  - the microphysical flux at the inversion is taken from the
!    flux grid-level just above it (assumes the divergence across
!    the inversion grid-level is physically within the BL)
!  - if no rad cooling was identified in layer, need to set
!    K_CLOUD_TOP and K_CLOUD_DSCT to top level in mixed layer
! ------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
        IF ( K_CLOUD_TOP(I,j)  ==  0 ) K_CLOUD_TOP(I,j) = NTML(I,j)

        FT_NT_ZH(I,j)   = FRAD(I,j,K_CLOUD_TOP(I,j)+1)                  &
     &                       + DF_INV_SML(I,j)
        FT_NT_ZH(I,j)   = FT_NT_ZH(I,j)   + FMIC(I,j,NTML(I,j)+2,1)     &
     &                                    + FSUBS(I,j,NTML(I,j),1)
        FQ_NT_ZH(I,j)   = FMIC(I,j,NTML(I,j)+2,2)                       &
     &                  + FSUBS(I,j,NTML(I,j),2)

        FT_NT_ZHSC(I,j) = 0.0
        FT_NT_ZHSC(I,j) = 0.0
        FQ_NT_ZHSC(I,j) = 0.0
        IF ( DSC(I,j) ) THEN
          IF ( K_CLOUD_DSCT(I,j)  ==  0 ) K_CLOUD_DSCT(I,j) = NTDSC(I,j)
          FT_NT_ZHSC(I,j) = FRAD(I,j,K_CLOUD_DSCT(I,j)+1)               &
     &                         + DF_INV_DSC(I,j)
          FT_NT_ZHSC(I,j) = FT_NT_ZHSC(I,j) + FMIC(I,j,NTDSC(I,j)+2,1)  &
     &                                      + FSUBS(I,j,NTDSC(I,j),1)
          FQ_NT_ZHSC(I,j) = FMIC(I,j,NTDSC(I,j)+2,2)                    &
     &                    + FSUBS(I,j,NTDSC(I,j),2)
        END IF
      end do
      end do
!-----------------------------------------------------------------------
! 10.  Subroutine EXCF_NL.
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
        NTML_SAVE(I,j) = NTML(I,j)  ! needed to identify changes
        DSC_SAVE(I,j)  = DSC(I,j)   !      in excf_nl
      END DO
      END DO
!
! DEPENDS ON: excf_nl
      CALL EXCF_NL (                                                    &
     & row_length,rows,BL_LEVELS,n_proc,FLUX_GRAD,NG_STRESS,            &
     & DSC,CUMULUS,COUPLED,                                             &
     & CF,NTML,ZH,ZHSC,DSCDEPTH,ZDSC_BASE,NTDSC,NBDSC,                  &
     & RDZ,Z_UV,Z_TQ,RHO_UV,RHO_TQ,RHOSTAR_GB,                          &
     & Z0M,V_S,FB_SURF,                                                 &
     & BTM,BQM,BTM_CLD,BQM_CLD,CF_SML,CF_DSC,                           &
     & BFLUX_SURF,BFLUX_SURF_SAT,ZETA_S,                                &
     & DF_TOP_OVER_CP,ZETA_R,BT_TOP,BTT_TOP,BTC_TOP,                    &
     & DB_TOP,DB_TOP_CLD,CHI_S_TOP,ZC,D_SIEMS,                          &
     & DF_DSCT_OVER_CP,ZETA_R_DSC,BT_DSCT,BTT_DSCT,BTC_DSCT,            &
     & DB_DSCT,DB_DSCT_CLD,CHI_S_DSCT,ZC_DSC,D_SIEMS_DSC,               &
     & DB_GA_DRY,DB_NOGA_DRY,DB_GA_CLD,DB_NOGA_CLD,                     &
     & DSL_SML,DQW_SML,DSL_DSC,DQW_DSC,                                 &
     & RHOKM, RHOKH, RHOKM_TOP, RHOKH_TOP,                              &
     & RHOKH_TOP_ENT, RHOKH_DSCT_ENT, RHOKH_SURF_ENT,                   &
     & FT_NT, FT_NT_ZH, FT_NT_ZHSC, FQ_NT, FQ_NT_ZH, FQ_NT_ZHSC,        &
     & RHOF2,RHOFSC,F_NGSTRESS,                                         &
     & LTIMER                                                           &
     &)

!-----------------------------------------------------------------------
!-adjust SML/DSC properties depending on diagnoses in EXCF_NL
! Note that the non-turbulent fluxes at inversions will have been
! swapped in EXCF_NL (ie. FT/Q_NT_ZH/ZHSC)
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
        IF ( DSC(I,j) .AND. .NOT.DSC_SAVE(I,j) ) THEN
!..decoupling diagnosed in EXCF_NL - change parameters around
          DSL_DSC(I,j) = DSL_SML(I,j)
          DQW_DSC(I,j) = DQW_SML(I,j)
          DB_DSCT(I,j) = DB_TOP(I,j)  ! copy diagnostics across
          ZC_DSC(I,j)  = ZC(I,j)      !
          DSC_DISC_INV(I,j) = SML_DISC_INV(I,j)
          SML_DISC_INV(I,j) = 0
          DSL_SML(I,j) = 0.0
          DQW_SML(I,j) = 0.0
          DF_INV_DSC(I,j) = DF_INV_SML(I,j)
          DF_INV_SML(I,j) = 0.0
        END IF
        IF ( .NOT.DSC(I,j) .AND. DSC_SAVE(I,j) ) THEN
!..decoupled layer removed in EXCF_NL; either...
          IF ( NTML_SAVE(I,j)  ==  NTML(I,j) ) THEN
!         !...had no turbulence forcing
            DSL_DSC(I,j) = 0.0
            DQW_DSC(I,j) = 0.0
            DSC_DISC_INV(I,j) = 0
            DF_INV_DSC(I,j) = 0.0
          ELSE
!         !...recoupled with surface layer
            DSL_SML(I,j) = DSL_DSC(I,j)
            DQW_SML(I,j) = DQW_DSC(I,j)
            DSL_DSC(I,j) = 0.0
            DQW_DSC(I,j) = 0.0
            SML_DISC_INV(I,j) = DSC_DISC_INV(I,j)
            DSC_DISC_INV(I,j) = 0
            DF_INV_SML(I,j) = DF_INV_DSC(I,j)
            DF_INV_DSC(I,j) = 0.0
          END IF
        END IF
      END DO
      END DO
!-----------------------------------------------------------------------
! 11.  Calculate "explicit" entrainment fluxes of SL and QW.
!-----------------------------------------------------------------------
! Calculate the non-turbulent fluxes at the DSC base
! ------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
         FT_NT_DSCB(I,j) = FT_NT(I,j,1)
         IF ( NBDSC(I,j)  >   1 ) THEN
           K = NBDSC(I,j)  ! NBDSC marks the lowest flux-level
!                          !    within the DSC layer
!                          ! Interpolate non-turb flux to base
!                          !    of DSC layer:
           FT_NT_DSCB(I,j) = FT_NT(I,j,K-1) +                           &
     &               (FT_NT(I,j,K)-FT_NT(I,j,K-1))                      &
     &              *(ZDSC_BASE(I,j)-Z_HALF(I,j,K-1))/DZL(I,j,K-1)
         END IF
      END DO
      END DO
!-----------------------------------------------------------------------
!..Specify entrainment fluxes at NTML+1 and NTDSC+1 directly through FTL
!..and FQW (and set the entrainment RHOKH to zero).
!..The turbulent flux at the subgrid ZH is given by the entrainment
!..parametrization ( = -w_e*'jump across inversion').  Together with
!..the non-turbulent flux profile (rad+microphys+subs), this gives the
!..total flux at the subgrid ZH.  The linear total flux profile is then
!..interpolated onto the half-level below ZH (the entrainment flux
!..grid-level).  The total flux divergence across the inversion grid
!..level is then checked for consistency with the entrainment/subsidence
!..balance, eg. a falling inversion should warm the inversion grid-level
!..Finally, the non-turbulent flux is sutracted off to give the required
!..turbulent component at the entrainment flux grid-level.
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
!-----------------------------------------------------------------------
!..First the surface-based mixed layer (if entraining and a
!..discontinuous inversion structure was diagnosed -
!..ie. the inversion is well-defined)
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length

        ZH_NP1(I,j)   = 0.0
        T_FRAC(I,j)   = 0.0
        ZRZI(I,j)     = 0.0
        WE_RHO(I,j)   = 0.0
        TOTHF_ZH(I,j) = 0.0
        TOTQF_ZH(I,j) = 0.0
        K=NTML(I,j)+1

!       ! Only RHOKH_ENT is passed out of EXCFNL so recalculate WE:
        WE_PARM(I,j) = RDZ(I,j,K)*                                      &
     &                     ( RHOKH_TOP_ENT(I,j)+RHOKH_SURF_ENT(I,j) )   &
     &                                                / RHO_UV(I,j,K)

        IF ( SML_DISC_INV(I,j)  ==  1 .AND. .NOT.COUPLED(I,j) .AND.     &
     &       (RHOKH_TOP_ENT(I,j)+RHOKH_SURF_ENT(I,j))  >   0.0 ) THEN

!-----------------------------------------------------------------------
!..Calculate ZH at end of timestep, ZH_NP1
!-----------------------------------------------------------------------
!..linearly interpolate vertical velocity to ZH
          IF ( ZH(I,j)  >=  Z_FULL(I,j,K) ) THEN
            W_LS(I,j) = W(I,j,K) + ( W(I,j,K+1) - W(I,j,K) )            &
     &                  * (ZH(I,j)-Z_FULL(I,j,K)) * RDZ(I,j,K+1)
          ELSE
            W_LS(I,j) = W(I,j,K) + ( W(I,j,K) - W(I,j,K-1) )            &
     &                  * (ZH(I,j)-Z_FULL(I,j,K)) * RDZ(I,j,K)
          END IF
          W_LS(I,j) = MIN ( W_LS(I,j), 0.0 )
!         ! only interested in subsidence

          ZH_NP1(I,j) = ZH(I,j) +                                       &
     &                    TIMESTEP * ( WE_PARM(I,j) + W_LS(I,j) )
          ZH_NP1(I,j) = MAX( ZH_NP1(I,j), Z_HALF(I,j,K-1) )
          IF ( ZH_NP1(I,j)  >   Z_TOP(I,j,K+1) ) THEN
!           ! limit ZH and W_e (and therefore the entraiment fluxes)
!           ! because the inversion cannot rise more than one level
!           ! in a timestep.
            ZH_NP1(I,j) = Z_TOP(I,j,K+1)
            WE_PARM(I,j) =                                              &
     &              (Z_TOP(I,j,K+1) - ZH(I,j))/TIMESTEP - W_LS(I,j)
          END IF
!-----------------------------------------------------------------------
!..Decide on which grid-level to apply entrainment flux
!-----------------------------------------------------------------------
          IF ( ZH_NP1(I,j)  >   Z_HALF(I,j,NTML(I,j)+2) ) THEN
!           ! ZH risen above level K+1 so specify appropriate flux
!           ! at this level and raise NTML by one (this means
!           ! gradient-adjustment is also applied at half-level
!           ! old_NTML+1).  Note KH profiles should already be
!           ! calculated at level NTML+1 because ZH is above this level.
            NTML(I,j) = NTML(I,j) + 1
            K=NTML(I,j)+1

!           ! T_FRAC is fraction of timestep inversion is above
!           ! the entrainment flux grid-level (at Z_HALF(K))
            T_FRAC(I,j) = (ZH_NP1(I,j)-Z_HALF(I,j,K)) /                 &
     &                    (ZH_NP1(I,j)-ZH(I,j))
!           ! ZH_FRAC is the timestep-average fraction of mixed layer
!           ! air in the inversion grid-level, level NTML+1
            ZH_FRAC(I,j) = 0.5*T_FRAC(I,j)*(ZH_NP1(I,j)-Z_HALF(I,j,K) ) &
     &                     / DZL(I,j,K)

          ELSE IF ( ZH_NP1(I,j)  >=  Z_HALF(I,j,NTML(I,j)+1) ) THEN
!           ! ZH always between half-levels NTML+1 and NTML+2

            T_FRAC(I,j) = 1.0
            ZH_FRAC(I,j) = ( 0.5*(ZH(I,j)+ZH_NP1(I,j)) - Z_HALF(I,j,K) )&
     &                     / DZL(I,j,K)

          ELSE
!           ! ZH falls below half-level NTML+1
            IF (NTML(I,j)  >=  2) THEN     ! FTL(K=1) is surface flux!
              NTML(I,j) = NTML(I,j) - 1
              K=NTML(I,j)+1
              RHOKH_TOP(I,j,K+1) = 0.0   ! also need to remove diffusion
              RHOKH(I,j,K+1)     = 0.0   ! at old entrainment grid-level
!             ! First calculate fraction of time below Z_HALF(NTML+2)
              T_FRAC(I,j) = (Z_HALF(I,j,K+1)-ZH_NP1(I,j)) /             &
     &                      (ZH(I,j)-ZH_NP1(I,j))
!             ! ZH_FRAC is still the timestep-average fraction of ML
!             ! air in grid-level NTML+1 (so close to 1 here)
              ZH_FRAC(I,j) = 1.0 -                                      &
     &           0.5*T_FRAC(I,j)*(Z_HALF(I,j,K)-ZH_NP1(I,j))/ DZL(I,j,K)
!             ! Now set T_FRAC to fraction of time above Z_HALF(NTML+1)
              T_FRAC(I,j)  = 1.0
            ELSE
!             ! keep implicit (diffusive) entrainment
              T_FRAC(I,j) = 0.0
              ZH_FRAC(I,j) = 0.0
              SML_DISC_INV(I,j) = 0
            END IF

          END IF  ! test on where to apply entrainment flux

          WE_RHO(I,j) = RHO_UV(I,j,K) * WE_PARM(I,j)
          ZRZI(I,j)   = Z_HALF(I,j,K)*2.0/(ZH(I,j)+ZH_NP1(I,j))

        END IF   ! test on SML_DISC_INV, etc
      END DO
      END DO
!-----------------------------------------------------------------------
!..Linearly interpolate between the known total (turb+rad+subs+micro)
!..flux at the surface and the parametrized flux at the inversion
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length

!       ! Entrainment flux applied to level NTML+1 which is the
!       ! flux-level above the top of the SML
        K=NTML(I,j)+1

        IF ( T_FRAC(I,j)  >   0.0 ) THEN

          RHOKH_TOP(I,j,K) = 0.0   ! apply entrainment explicitly
          RHOKH(I,j,K) = 0.0       !      "

          TOTHF_ZH(I,j) = - WE_RHO(I,j)*DSL_SML(I,j) + FT_NT_ZH(I,j)
!         ! Linearly interpolate to entrainment flux grid-level
          TOTHF_EFL = FT_NT(I,j,1) + FTL(I,j,1) +                       &
     &               ( TOTHF_ZH(I,j)-FT_NT(I,j,1)-FTL(I,j,1) )*ZRZI(I,j)
!         ! Ensure total heat flux gradient in inversion grid-level is
!         ! consistent with inversion rising (ie. implies cooling in
!         ! level K relative to the mixed layer) or falling
!         ! (implies warming)

          ML_TEND = -( TOTHF_ZH(I,j)-FT_NT(I,j,1)-FTL(I,j,1) ) / ZH(I,j)
          FA_TEND = 0.0
          IF ( K+1  <=  BL_LEVELS )                                     &
     &        FA_TEND = - ( FT_NT(I,j,K+2) - FT_NT(I,j,K+1) )           &
     &                    / DZL(I,j,K+1)
          INV_TEND =       ZH_FRAC(I,j) * ML_TEND                       &
     &             + (1.0-ZH_FRAC(I,j)) * FA_TEND
          IF (WE_PARM(I,j)+W_LS(I,j)  >=  0.0) THEN
!           ! Ensure inversion level does cool relative to ML
            TOTHF_EFL = MIN( TOTHF_EFL,                                 &
     &                       FT_NT(I,j,K+1)+INV_TEND*DZL(I,j,K) )
!           ! Ensure inversion level won't end up colder than
!           ! NTML by end of timestep.
!           ! Set INV_TEND to max allowable cooling rate, also
!           ! allowing for change in ML_TEND arising from this change
!           ! to TOTHF_EFL:
            INV_TEND = (SL(I,j,K-1)-SL(I,j,K))/TIMESTEP                 &
     &                  + (FT_NT(I,j,1)+FTL(I,j,1))/Z_HALF(I,j,K)
            TOTHF_EFL = MAX( TOTHF_EFL,                                 &
     &                       (FT_NT(I,j,K+1)+INV_TEND*DZL(I,j,K))       &
     &                            /(1.+ DZL(I,j,K)/Z_HALF(I,j,K)) )
          ELSE  ! WE_PARM+W_LS < 0
!           ! Ensure inversion level does warm relative to ML
            TOTHF_EFL = MAX( TOTHF_EFL,                                 &
     &                       FT_NT(I,j,K+1)+INV_TEND*DZL(I,j,K) )
          END IF
!         ! Turbulent entrainment flux is then the residual of the total
!         ! flux and the net flux from other processes
          FTL(I,j,K) =  T_FRAC(I,j) * ( TOTHF_EFL - FT_NT(I,j,K) )

        ELSE   ! NOT specifying entrainment flux but KH
!         ! Add entrainment KH to K-profiles
!         ! (for COUPLED layers these will be zero)
          RHOKH_TOP(I,j,K) = RHOKH_TOP(I,j,K) + RHOKH_TOP_ENT(I,j)
          RHOKH(I,j,K) = RHOKH(I,j,K) + RHOKH_SURF_ENT(I,j)

        END IF  ! test on T_FRAC gt 0

      END DO
      END DO
!-------------------------------------------------
!..Second the decoupled mixed layer, if entraining
!-------------------------------------------------
      do j=1,rows
      do i=1,row_length

        ZHSC_NP1(I,j)   = 0.0
        T_FRAC_DSC(I,j) = 0.0
        ZRZI_DSC(I,j)   = 0.0
        WE_RHO_DSC(I,j) = 0.0
        TOTHF_ZHSC(I,j) = 0.0
        TOTQF_ZHSC(I,j) = 0.0

        K=NTDSC(I,j)+1
        WE_DSC_PARM(I,j) = RDZ(I,j,K)*RHOKH_DSCT_ENT(I,j)               &
     &                                             / RHO_UV(I,j,K)

        IF ( DSC_DISC_INV(I,j)  ==  1                                   &
     &                .AND. RHOKH_DSCT_ENT(I,j)  >   0.0 ) THEN

!-----------------------------------------------------------------------
!..Calculate ZHSC at end of timestep, ZHSC_NP1
!-----------------------------------------------------------------------
!..interpolate vertical velocity to ZH
          IF ( ZHSC(I,j)  >=  Z_FULL(I,j,K) ) THEN
            W_LS_DSC(I,j) = W(I,j,K) + ( W(I,j,K+1) - W(I,j,K) ) *      &
     &                       (ZHSC(I,j)-Z_FULL(I,j,K)) * RDZ(I,j,K+1)
          ELSE
            W_LS_DSC(I,j) = W(I,j,K) + ( W(I,j,K) - W(I,j,K-1) ) *      &
     &                       (ZHSC(I,j)-Z_FULL(I,j,K)) * RDZ(I,j,K)
          END IF
          W_LS_DSC(I,j) = MIN ( W_LS_DSC(I,j), 0.0 )
!         ! only interested in subsidence

          ZHSC_NP1(I,j) = ZHSC(I,j) +                                   &
     &          TIMESTEP * ( WE_DSC_PARM(I,j) + W_LS_DSC(I,j) )
          ZHSC_NP1(I,j) = MAX( ZHSC_NP1(I,j), Z_HALF(I,j,K-1) )
          IF ( ZHSC_NP1(I,j)  >   Z_TOP(I,j,K+1) ) THEN
!           ! limit ZHSC and W_e (and therefore the entrainment fluxes)
!           ! because the inversion cannot rise more than one level
!           ! in a timestep.
            ZHSC_NP1(I,j) = Z_TOP(I,j,K+1)
            WE_DSC_PARM(I,j) =                                          &
     &         (Z_TOP(I,j,K+1) - ZHSC(I,j))/TIMESTEP - W_LS_DSC(I,j)
          END IF
!-----------------------------------------------------------------------
!..Decide on which grid-level to apply entrainment flux
!-----------------------------------------------------------------------
          IF ( ZHSC_NP1(I,j)  >   Z_HALF(I,j,NTDSC(I,j)+2) ) THEN
!           ! ZHSC risen above level K+1 so specify appropriate
!           ! flux at this level and raise NTDSC by one

            NTDSC(I,j) = NTDSC(I,j) + 1
            K = NTDSC(I,j)+1
            T_FRAC_DSC(I,j) = (ZHSC_NP1(I,j)-Z_HALF(I,j,K)) /           &
     &                        (ZHSC_NP1(I,j)-ZHSC(I,j))

            ZHSC_FRAC(I,j) = 0.5*T_FRAC_DSC(I,j)*                       &
     &                       ( ZHSC_NP1(I,j)-Z_HALF(I,j,K) )/ DZL(I,j,K)

          ELSE IF ( ZHSC_NP1(I,j)  >   Z_HALF(I,j,NTDSC(I,j)+1) ) THEN
!           ! ZHSC always between half-levels NTDSC+1 and NTDSC+2

            T_FRAC_DSC(I,j) = 1.0
            ZHSC_FRAC(I,j) = ( 0.5*(ZHSC(I,j)+ZHSC_NP1(I,j))            &
     &                                     - Z_HALF(I,j,K) )/ DZL(I,j,K)

          ELSE
!           ! ZHSC falls below half-level NTDSC+1
            NTDSC(I,j) = NTDSC(I,j) - 1  ! could reduce NTDSC to 1
            K = NTDSC(I,j)+1
            RHOKH_TOP(I,j,K+1) = 0.0
            RHOKH(I,j,K+1)     = 0.0

!           ! First calculate fraction of time below Z_HALF(NTDSC+2)
            T_FRAC_DSC(I,j) = (Z_HALF(I,j,K+1)-ZHSC_NP1(I,j)) /         &
     &                        (ZHSC(I,j)-ZHSC_NP1(I,j))
!           ! ZHSC_FRAC is still the timestep-average fraction of ML
!           ! air in grid-level NTDSC+1 (so close to 1 here)
            ZHSC_FRAC(I,j) = 1.0 - 0.5*T_FRAC_DSC(I,j)*                 &
     &                      (Z_HALF(I,j,K)-ZHSC_NP1(I,j))/ DZL(I,j,K)
!           ! Now set T_FRAC to fraction of time above Z_HALF(NTDSC+1)
            T_FRAC_DSC(I,j)  = 1.0

          END IF  ! test on where to apply entrainment flux

          WE_RHO_DSC(I,j) = RHO_UV(I,j,K) * WE_DSC_PARM(I,j)
!         ! for z'/z_i' assume height of DSC base is fixed in time
          ZRZI_DSC(I,j) =( Z_HALF(I,j,K)-(ZHSC(I,j)-DSCDEPTH(I,j)) )    &
     &                  /( DSCDEPTH(I,j)+0.5*(ZHSC_NP1(I,j)-ZHSC(I,j)) )

        END IF   ! test on DSC_DISC_INV, etc
      END DO
      END DO
!-----------------------------------------------------------------------
!..Linearly interpolate between the known total (turb+rad+subs+micro)
!..flux at the DSC base and the parametrized flux at the inversion
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length

!       ! Entrainment flux applied to level NTDSC+1 which is the
!       ! flux-level above the top of the DSC layer
        K=NTDSC(I,j)+1

        IF ( T_FRAC_DSC(I,j)  >   0.0 ) THEN


          RHOKH_TOP(I,j,K) = 0.0   ! apply entrainment explicitly
          RHOKH(I,j,K)     = 0.0   !      "

          TOTHF_ZHSC(I,J) = - WE_RHO_DSC(I,j)*DSL_DSC(I,j)              &
     &                            + FT_NT_ZHSC(I,j)
          TOTHF_EFL = FT_NT_DSCB(I,j) +                                 &
     &                ( TOTHF_ZHSC(I,J)-FT_NT_DSCB(I,j) )*ZRZI_DSC(I,j)
!         ! Ensure total heat flux gradient in inversion grid-level is
!         ! consistent with inversion rising (implies cooling in
!         ! level K, relative to mixed layer) or falling
!         ! (implies warming)
          ML_TEND = - ( TOTHF_ZHSC(I,J)-FT_NT_DSCB(I,j) )/ DSCDEPTH(I,j)
          FA_TEND = 0.0
          IF ( K+1  <=  BL_LEVELS )                                     &
     &        FA_TEND = - ( FT_NT(I,j,K+2) - FT_NT(I,j,K+1) )           &
     &                    / DZL(I,j,K+1)
          INV_TEND =       ZHSC_FRAC(I,j) * ML_TEND                     &
     &             + (1.0-ZHSC_FRAC(I,j)) * FA_TEND
          IF (WE_DSC_PARM(I,j)+W_LS_DSC(I,j)  >=  0.0) THEN
!           ! Ensure inversion level does cool relative to ML
            TOTHF_EFL = MIN( TOTHF_EFL,                                 &
     &                       FT_NT(I,j,K+1)+INV_TEND*DZL(I,j,K) )
!           ! Ensure inversion level won't end up colder than
!           ! NTDSC by end of timestep.
            INV_TEND = (SL(I,j,K-1)-SL(I,j,K))/TIMESTEP                 &
     &                 + FT_NT_DSCB(I,j)/DSCDEPTH(I,j)
            TOTHF_EFL = MAX( TOTHF_EFL,                                 &
     &                    (FT_NT(I,j,K+1)+INV_TEND*DZL(I,j,K))          &
     &                     /(1.+ DZL(I,j,K)/DSCDEPTH(I,j))   )
          ELSE   ! WE_DSC_PARM+W_LS_DSC < 0
            TOTHF_EFL = MAX( TOTHF_EFL,                                 &
     &                       FT_NT(I,j,K+1)+INV_TEND*DZL(I,j,K) )
          END IF
!         ! Turbulent entrainment flux is then the residual of the total
!         ! flux and the net flux from other processes
          FTL(I,j,K) = T_FRAC_DSC(I,j) * ( TOTHF_EFL - FT_NT(I,j,K) )

        ELSEIF ( DSC(I,j) ) THEN

!         ! not specifying entrainment flux but KH
          RHOKH_TOP(I,j,K) = RHOKH_DSCT_ENT(I,j)

        END IF  ! IF NOT DSC

      END DO
      END DO
!-----------------------------------------------------------------------
! Specify QW entrainment fluxes
!-----------------------------------------------------------------------
! Calculate the non-turbulent fluxes at the layer boundaries
!  - the subsidence flux at the inversion is taken from the
!    flux grid-level below it (assumes the divergence across
!    the inversion is physically above the BL)
!  - the microphysical flux at the inversion is taken from the
!    flux grid-level just above it (assumes the divergence across
!    the inversion grid-level is physically within the BL)
! ------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
         FQ_NT_DSCB(I,j) = FQ_NT(I,j,1)
         IF ( NBDSC(I,j)  >   1 ) THEN
           K = NBDSC(I,j)  ! NBDSC marks the lowest flux-level
!                          !    within the DSC layer
!                          ! Interpolate non-turb flux to base
!                          !    of DSC layer:
           FQ_NT_DSCB(I,j) = FQ_NT(I,j,K-1) +                           &
     &               (FQ_NT(I,j,K)-FQ_NT(I,j,K-1))                      &
     &              *(ZDSC_BASE(I,j)-Z_HALF(I,j,K-1))/DZL(I,j,K-1)
         END IF
      END DO
      END DO
!-----------------------------------------------------------------------
! Calculate grid-level QW fluxes at inversion, ensuring the turbulent,
! microphysical and subsidence fluxes are correctly coupled.
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
        K = NTML(I,j)+1
        IF ( T_FRAC(I,j)  >   0.0 ) THEN
!         ! Calculate total (turb+micro+subs) QW flux at subgrid
!         ! inversion height
          TOTQF_ZH(I,J) = - WE_RHO(I,j)*DQW_SML(I,j) + FQ_NT_ZH(I,j)
!         ! Interpolate to entrainment flux-level below
          TOTQF_EFL = FQ_NT(I,j,1) + FQW(I,j,1) +  ZRZI(I,j) *          &
     &                     ( TOTQF_ZH(I,J) - FQ_NT(I,j,1) - FQW(I,j,1) )
!         ! Need to ensure the total QW flux gradient in inversion
!         ! grid-level is consistent with inversion rising or falling.
!         ! If QW(K) is drier than mixed layer then inversion rising
!         ! implies moistening in level K relative to mixed layer
!         ! while falling would imply relative drying of level K.
!         ! If QW(K) is moister than ML then want opposite tendencies.
          ML_TEND = - ( TOTQF_ZH(I,J)-FQ_NT(I,j,1)-FQW(I,j,1) ) /ZH(I,j)
          FA_TEND = 0.0
          IF ( K+1  <=  BL_LEVELS )                                     &
     &      FA_TEND = - ( FQ_NT(I,j,K+2)-FQ_NT(I,j,K+1) )               &
     &                  / DZL(I,j,K+1)
          INV_TEND =       ZH_FRAC(I,j) * ML_TEND                       &
     &             + (1.0-ZH_FRAC(I,j)) * FA_TEND

          MOISTEN = (WE_PARM(I,j)+W_LS(I,j)  >=  0.0)
          IF ( QW(I,j,K) >  QW(I,j,K-1) ) MOISTEN = .NOT.MOISTEN

          IF ( MOISTEN ) THEN
!           ! Ensure inversion level does moisten relative to ML
            TOTQF_EFL = MAX( TOTQF_EFL,                                 &
     &                      FQ_NT(I,j,K+1)+INV_TEND*DZL(I,j,K) )
            IF (WE_PARM(I,j)+W_LS(I,j)  >=  0.0) THEN
!             ! Ensure inversion level won't end up more moist than
!             ! NTML by end of timestep.
!             ! Set INV_TEND to max allowable moistening rate, also
!             ! allowing for change in ML_TEND arising from this change
!             ! to TOTQF_EFL:
              INV_TEND = (QW(I,j,K-1)-QW(I,j,K))/TIMESTEP               &
     &                    + (FQ_NT(I,j,1)+FQW(I,j,1))/Z_HALF(I,j,K)
              TOTQF_EFL = MIN( TOTQF_EFL,                               &
     &                      (FQ_NT(I,j,K+1)+INV_TEND*DZL(I,j,K))        &
     &                       /(1.+ DZL(I,j,K)/Z_HALF(I,j,K))   )
            END IF
          ELSE
            TOTQF_EFL = MIN( TOTQF_EFL,                                 &
     &                      FQ_NT(I,j,K+1)+INV_TEND*DZL(I,j,K) )
            IF (WE_PARM(I,j)+W_LS(I,j)  >=  0.0) THEN
!             ! Ensure inversion level won't end up drier than
!             ! NTML by end of timestep.
!             ! Set INV_TEND to max allowable drying rate:
              INV_TEND = (QW(I,j,K-1)-QW(I,j,K))/TIMESTEP               &
     &                    + (FQ_NT(I,j,1)+FQW(I,j,1))/Z_HALF(I,j,K)
              TOTQF_EFL = MAX( TOTQF_EFL,                               &
     &                      (FQ_NT(I,j,K+1)+INV_TEND*DZL(I,j,K))        &
     &                       /(1.+ DZL(I,j,K)/Z_HALF(I,j,K))   )
            END IF
          END IF
          FQW(I,j,K) = T_FRAC(I,j) *                                    &
     &                     ( TOTQF_EFL - FQ_NT(I,j,K) )
        END IF

      END DO
      END DO
!-----------------------------------------------------------------------
! Now decoupled layer
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length

        IF ( T_FRAC_DSC(I,j)  >   0.0 ) THEN

          K = NTDSC(I,j)+1

!         ! Calculate total (turb+micro) QW flux at subgrid inversion
          TOTQF_ZHSC(I,J) = - WE_RHO_DSC(I,j)*DQW_DSC(I,j)              &
     &                        + FQ_NT_ZHSC(I,j)
!         ! Interpolate to entrainment flux-level
          TOTQF_EFL = FQ_NT_DSCB(I,j) +                                 &
     &              ( TOTQF_ZHSC(I,J) - FQ_NT_DSCB(I,j) )*ZRZI_DSC(I,j)

          ML_TEND = - ( TOTQF_ZHSC(I,J)-FQ_NT_DSCB(I,j) )/DSCDEPTH(I,j)
          FA_TEND = 0.0
          IF ( K+1  <=  BL_LEVELS )                                     &
     &       FA_TEND = - ( FQ_NT(I,j,K+2)-FQ_NT(I,j,K+1) )              &
     &                   / DZL(I,j,K+1)
          INV_TEND =       ZHSC_FRAC(I,j) * ML_TEND                     &
     &             + (1.0-ZHSC_FRAC(I,j)) * FA_TEND

          MOISTEN = (WE_DSC_PARM(I,j)+W_LS_DSC(I,j)  >=  0.0)
          IF ( QW(I,j,K) >  QW(I,j,K-1) ) MOISTEN = .NOT.MOISTEN
          IF ( MOISTEN ) THEN
            TOTQF_EFL = MAX( TOTQF_EFL,                                 &
     &                        FQ_NT(I,j,K+1)+INV_TEND*DZL(I,j,K) )
            IF (WE_DSC_PARM(I,j)+W_LS_DSC(I,j)  >=  0.0) THEN
!             ! Ensure inversion level won't end up more moist than
!             ! NTDSC by end of timestep.
              INV_TEND = (QW(I,j,K-1)-QW(I,j,K))/TIMESTEP               &
     &                   + FQ_NT_DSCB(I,j)/DSCDEPTH(I,j)
              TOTQF_EFL = MIN( TOTQF_EFL,                               &
     &                      (FQ_NT(I,j,K+1)+INV_TEND*DZL(I,j,K))        &
     &                       /(1.+ DZL(I,j,K)/DSCDEPTH(I,j))   )
            END IF
          ELSE
            TOTQF_EFL = MIN( TOTQF_EFL,                                 &
     &                        FQ_NT(I,j,K+1)+INV_TEND*DZL(I,j,K) )
            IF (WE_DSC_PARM(I,j)+W_LS_DSC(I,j)  >=  0.0) THEN
!             ! Ensure inversion level won't end up drier than
!             ! NTDSC by end of timestep.
!             ! Set INV_TEND to max allowable drying rate:
              INV_TEND = (QW(I,j,K-1)-QW(I,j,K))/TIMESTEP               &
     &                   + FQ_NT_DSCB(I,j)/DSCDEPTH(I,j)
              TOTQF_EFL = MAX( TOTQF_EFL,                               &
     &                      (FQ_NT(I,j,K+1)+INV_TEND*DZL(I,j,K))        &
     &                       /(1.+ DZL(I,j,K)/DSCDEPTH(I,j))   )
            END IF
          END IF
          FQW(I,j,K) = T_FRAC_DSC(I,j) * ( TOTQF_EFL - FQ_NT(I,j,K) )
        END IF
      END DO
      END DO
!-----------------------------------------------------------------------
! Calculate effective entrainment (ie. reduced to allow for subsidence
! increments in the ML) for use in tracer mixing.  Take theta_l as a
! representative scalar field since jump should always be the same sign
! and code therefore simpler.
! In this version the inversion fluxes are only implemented at one
! grid-level so only one element of these 3D arrays is used.
!-----------------------------------------------------------------------
      do K=1,3
      do j=1,rows
      do i=1,row_length
        WE_LIM(I,j,K)    = 0.0
        T_FRAC_TR(I,j,K) = 0.0
        ZRZI_TR(I,j,K)   = 0.0
        WE_LIM_DSC(I,j,K)    = 0.0
        T_FRAC_DSC_TR(I,j,K) = 0.0
        ZRZI_DSC_TR(I,j,K)   = 0.0
      END DO
      END DO
      END DO
      do j=1,rows
      do i=1,row_length
        KENT(I,j)   = NTML(I,j)+1
        T_FRAC_TR(I,j,2) = T_FRAC(I,j)
        ZRZI_TR(I,j,2) = ZRZI(I,j)
        IF ( T_FRAC(I,j)  >   0.0 ) THEN
          W_S_ENT = 0.0
          K = NTML(I,j)
          IF ( DSL_SML(I,j)  /=  0.0 ) W_S_ENT =                        &
     &        MIN( 0.0, -SLS_INC(I,j,K) * DZL(I,j,K) /DSL_SML(I,j) )
!         ! Only allow w_e to be reduced to zero!
          WE_LIM(I,j,2) = RHO_UV(I,j,K+1) *                             &
     &                      MAX( 0.0, WE_PARM(I,j) + W_S_ENT )
        END IF
        KENT_DSC(I,j)   = NTDSC(I,j)+1
        T_FRAC_DSC_TR(I,j,2) = T_FRAC_DSC(I,j)
        ZRZI_DSC_TR(I,j,2) = ZRZI_DSC(I,j)
        IF ( T_FRAC_DSC(I,j)  >   0.0 ) THEN
          W_S_ENT = 0.0
          K = NTDSC(I,j)
          IF ( DSL_DSC(I,j)  /=  0.0 ) W_S_ENT =                        &
     &        MIN( 0.0, -SLS_INC(I,j,K) * DZL(I,j,K) /DSL_DSC(I,j) )
!         ! Only allow w_e to be reduced to zero!
          WE_LIM_DSC(I,j,2) = RHO_UV(I,j,K) *                           &
     &                      MAX( 0.0, WE_DSC_PARM(I,j) + W_S_ENT )
        END IF
      END DO
      END DO
!-----------------------------------------------------------------------
! 12. Update standard deviations and gradient adjustment to use this
!     timestep's ZH (code from SF_EXCH)
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
        IF ( UNSTABLE(I,j) ) THEN
          IF (FLUX_GRAD  ==  Locketal2000) THEN
            W_S_CUBED = 0.25 * ZH(I,j) * FB_SURF(I,j)
            IF (W_S_CUBED  >   0.0) THEN
              W_M  =                                                    &
     &       ( W_S_CUBED + V_S(I,j) * V_S(I,j) * V_S(I,j) ) ** (1.0/3.0)
              T1_SD(I,j) = 1.93 * FTL(I,j,1) / (RHOSTAR_GB(I,j) * W_M)
              Q1_SD(I,j) = 1.93 * FQW(I,j,1) / (RHOSTAR_GB(I,j) * W_M)
              TV1_SD(I,j) = T(I,j,1) *                                  &
     &          ( 1.0 + C_VIRTUAL*Q(I,j,1) - QCL(I,j,1) - QCF(I,j,1) ) *&
     &          ( BT(I,j,1)*T1_SD(I,j) + BQ(I,j,1)*Q1_SD(I,j) )
              T1_SD(I,j) = MAX ( 0.0 , T1_SD(I,j) )
              Q1_SD(I,j) = MAX ( 0.0 , Q1_SD(I,j) )
              IF (TV1_SD(I,j)  <=  0.0) THEN
                TV1_SD(I,j) = 0.0
                T1_SD(I,j) = 0.0
                Q1_SD(I,j) = 0.0
              END IF
            END IF
            GRAD_T_ADJ(I,j) = MIN( MAX_T_GRAD ,                         &
     &                       A_GRAD_ADJ * T1_SD(I,j) / ZH(I,j) )
            GRAD_Q_ADJ(I,j) = 0.0
          ELSEIF (FLUX_GRAD  ==  HoltBov1993) THEN
!           ! Use constants from Holtslag and Boville (1993)
!           ! Conv limit GAMMA_TH = 10 *FTL1/(wstar*zh)
!           ! Neut limit GAMMA_TH = 7.2*wstar*FTL1/(ustar^2*zh)
            WSTAR3 = FB_SURF(I,j) * ZH(I,j)
            W_M =( V_S(I,j)**3.0 + 0.6*WSTAR3 )**(1.0/3.0)

            GRAD_T_ADJ(I,j) = A_GA_HB93*(WSTAR3**(1.0/3.0))*FTL(I,j,1)  &
     &                        / ( RHOSTAR_GB(I,j)*W_M*W_M*ZH(I,j) )
!            GRAD_Q_ADJ(I,j) = A_GA_HB93*(WSTAR3**(1.0/3.0))*FQW(I,j,1)
!     &                       / ( RHOSTAR_GB(I,j)*W_M*W_M*ZH(I,j) )
!          ! Set q term to zero for same empirical reasons as Lock et al
            GRAD_Q_ADJ(I,j) = 0.0
          ELSEIF (FLUX_GRAD  ==  LockWhelan2006) THEN
!           ! Use constants LockWhelan2006
!           ! Conv limit GAMMA_TH = 10 *FTL1/(wstar*zh)
!           ! Neut limit GAMMA_TH = 7.5*FTL1/(ustar*zh)
            WSTAR3 = FB_SURF(I,j) * ZH(I,j)
            W_H =( ((4.0/3.0)*V_S(I,j))**3.0 + WSTAR3 )**(1.0/3.0)

            GRAD_T_ADJ(I,j) = A_GA_LW06 * FTL(I,j,1)                    &
     &                         / ( RHOSTAR_GB(I,j)*W_H*ZH(I,j) )
            GRAD_Q_ADJ(I,j) = A_GA_LW06 * FQW(I,j,1)                    &
     &                         / ( RHOSTAR_GB(I,j)*W_H*ZH(I,j) )
          END IF
        END IF  ! test on UNSTABLE
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
