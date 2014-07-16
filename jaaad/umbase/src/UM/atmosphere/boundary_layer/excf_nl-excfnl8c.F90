#if defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!
!!!  SUBROUTINE EXCF_NL ----------------------------------------------
!!!
!!!  Purpose: To calculate non-local exchange coefficients,
!!!           entrainment parametrization and non-gradient flux terms.
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!  6.2    Feb 2006  8C version of EXCF_NL.
!!!                   Includes a separate integration of the buoyancy
!!!                   flux close to the BL top and incorporates option
!!!                   for a revised flux-gradient relationship,
!!!                   generic to all scalars.              Adrian Lock
!!!  6.4    Nov 2006  Remove out-of-bounds errors.         Adrian Lock
!!!
!!!  Programming standard:
!!!
!!!  Documentation: UMDP No.24
!!!
!!!---------------------------------------------------------------------
!
!!   Arguments :-
      SUBROUTINE EXCF_NL (                                              &
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
!
      IMPLICIT NONE

      LOGICAL LTIMER

      INTEGER                                                           &
     & row_length,rows                                                  &
     &,n_proc                                                           &
                   ! IN Total number of processors
     &,BL_LEVELS                                                        &
                   ! IN maximum number of boundary layer levels
     &,NG_STRESS                                                        &
                   ! IN switch for non-gradient stress
     &,FLUX_GRAD   ! IN switch for flux-gradient formulation

      REAL, DIMENSION(row_length,rows,BL_LEVELS) ::                     &
     & CF                                                               &
                                ! IN Cloud fraction
     &,RDZ                                                              &
                                ! IN Reciprocal of distance between
                                !    T,q-levels (m^-1). 1/RDZ(,K) is
                                !    the vertical distance from level
                                !    K-1 to level K, except that for
                                !    K=1 it is just the height of the
                                !    lowest atmospheric level.
     &,Z_UV                                                             &
                                ! IN For a vertically staggered grid
                                !    with a u,v-level first above the
                                !    surface, Z_UV(*,K) is the height
                                !    of the k-th u,v-level (half level
                                !    k-1/2) above the surface
     &,Z_TQ                                                             &
                                ! IN For a vertically staggered grid
                                !    with a u,v-level first above the
                                !    surface, Z_TQ(*,K) is the height
                                !    of the k-th T,q-level (full level
                                !    k) above the surface
     &,RHO_UV                                                           &
                                ! IN For a vertically staggered grid
                                !    with a u,v-level first above the
                                !    surface, RHO_UV(*,K) is the
                                !    density at the k-th u,v-level
                                !    above the surface
     &,RHO_TQ                                                           &
                                ! IN For a vertically staggered grid
                                !    with a u,v-level first above the
                                !    surface, RHO_TQ(*,K) is the
                                !    density of the k-th T,q-level
                                !    above the surface
     &,BQM                                                              &
                                ! IN Buoyancy parameters for clear and
     &,BTM                                                              &
                                !    cloudy air on half levels
     &,BQM_CLD                                                          &
                                !    (*,K) elements are k+1/2 values
     &,BTM_CLD                  !

      REAL, DIMENSION(row_length,rows) ::                               &
     & RHOSTAR_GB                                                       &
                                ! IN Surface density (kg/m3)
     &,Z0M                                                              &
                                ! IN Roughness length for momentum (m).
     &,V_S                                                              &
                                ! IN Surface friction velocity (m/s).
     &,FB_SURF                                                          &
                                ! IN Buoyancy flux at the surface over
!                               !    density (m^2/s^3).
     &,BFLUX_SURF                                                       &
                                ! IN Surface buoyancy flux (kg/m/s^3).
     &,BFLUX_SURF_SAT                                                   &
                                ! IN Saturated-air surface buoyancy
                                !    flux.
     &,DB_TOP                                                           &
                                ! IN Buoyancy jump across the top of
                                !    the SML (m/s^2).
     &,DF_TOP_OVER_CP                                                   &
                                ! IN Radiative flux change at cloud top
                                !    divided by c_P (K.kg/m^2/s).
     &,BT_TOP                                                           &
                                ! IN Buoyancy parameter at the top of
                                !    the b.l. (m/s^2/K).
     &,BTT_TOP                                                          &
                                ! IN In-cloud buoyancy parameter at
                                !    the top of the b.l. (m/s^2/K).
     &,BTC_TOP                                                          &
                                ! IN Cloud fraction weighted buoyancy
                                !    parameter at the top of the b.l.
     &,DB_TOP_CLD                                                       &
                                ! IN In-cloud buoyancy jump at the
                                !    top of the b.l. (m/s^2).
     &,CHI_S_TOP                                                        &
                                ! IN Mixing fraction of just saturated
                                !    mixture at top of the b.l.
     &,ZETA_S                                                           &
                                ! IN Non-cloudy fraction of mixing
                                !    layer for surface forced
                                !    entrainment term.
     &,ZETA_R                                                           &
                                ! IN Non-cloudy fraction of mixing
                                !    layer for cloud top radiative
                                !    cooling entrainment term.
     &,ZC                                                               &
                                ! IN Cloud depth (not cloud fraction
                                !    weighted) (m).
     &,DB_DSCT                                                          &
                                ! IN Buoyancy jump across the top of
                                !    the DSC layer (m/s^2).
     &,DF_DSCT_OVER_CP                                                  &
                                ! IN Radiative flux change at DSC top
                                !    divided by c_P (K.kg/m^2/s).
     &,BT_DSCT                                                          &
                                ! IN Buoyancy parameters at the top of
     &,BTT_DSCT                                                         &
                                !    the DSC layer (m/s^2/K)
     &,BTC_DSCT                                                         &
                                !
     &,DB_DSCT_CLD                                                      &
                                ! IN In-cloud buoyancy jump at the
                                !    top of the DSC (m/s^2).
     &,CHI_S_DSCT                                                       &
                                ! IN Mixing fraction of just saturated
                                !    mixture at top of the DSC
     &,ZETA_R_DSC                                                       &
                                ! IN Non-cloudy fraction of DSC
                                !    for cloud top radiative
                                !    cooling entrainment term.
     &,ZC_DSC                                                           &
                                ! IN Cloud depth (not cloud fraction
                                !    weighted) (m).
     &,D_SIEMS                                                          &
                                ! IN Siems (1990) et al. cloud-top
     &,D_SIEMS_DSC                                                      &
                                !    entrainmentt instability params
     &,DQW_SML                                                          &
                                ! IN QW change across SML disc inv
     &,DSL_SML                                                          &
                                ! IN SL change across SML disc inv
     &,DQW_DSC                                                          &
                                ! IN QW change across DSC disc inv
     &,DSL_DSC                                                          &
                                ! IN SL change across DSC disc inv
     &,CF_SML                                                           &
                                ! IN cloud fraction of SML
     &,CF_DSC                   ! IN cloud fraction of DSC layer

      REAL, DIMENSION(row_length,rows, 2:BL_LEVELS) ::                  &
     & DB_GA_DRY                                                        &
                                ! IN Cloudy and cloud-free buoyancy
     &,DB_NOGA_DRY                                                      &
                                !    jumps for flux integral
     &,DB_GA_CLD                                                        &
                                !    calculation (m/s2):
     &,DB_NOGA_CLD
!
      REAL, DIMENSION(row_length,rows, BL_LEVELS+1) ::                  &
     & FT_NT                                                            &
                                ! IN Non-turbulent heat (rho*Km/s) and
     &,FQ_NT                    !      moisture (rho*m/s) fluxes

      LOGICAL, DIMENSION(row_length,rows) ::                            &
     & COUPLED                                                          &
                                ! INOUT Flag to indicate Sc layer
!                               !       weakly decoupled
     &,CUMULUS                                                          &
                                ! INOUT Flag for cumulus
     &,DSC                      ! INOUT Flag set if decoupled stratocu
!                               !       layer found.

      INTEGER, DIMENSION(row_length,rows) ::                            &
     & NTML                                                             &
                                ! INOUT  Number of turbulently mixed
!                               !        layers.
     &,NTDSC                    ! INOUT  Top level of any decoupled
!                               !        turbulently mixed Sc layer
!
      REAL, DIMENSION(row_length,rows) ::                               &
     & ZHSC                                                             &
                                ! INOUT Cloud-layer height (m)
     &,ZH                                                               &
                                ! INOUT Boundary layer height (m)
     &,DSCDEPTH                                                         &
                                ! INOUT Decoupled cloud-layer depth (m)
     &,FT_NT_ZH                                                         &
                                ! INOUT Non-turbulent heat (rho*Km/s)
     &,FT_NT_ZHSC                                                       &
                                !        and moisture (rho*m/s) fluxes
     &,FQ_NT_ZH                                                         &
                                !        evaluated at the SML and DSC
     &,FQ_NT_ZHSC               !        inversions

      INTEGER                                                           &
     & NBDSC(row_length,rows)   ! OUT Bottom level of any decoupled
!                               !     turbulently mixed Sc layer.

      REAL, DIMENSION(row_length,rows, 2:BL_LEVELS) ::                  &
     & RHOKM                                                            &
                                ! OUT Layer k-1 - to - layer k
                                !     turbulent mixing coefficient
                                !     for momentum (kg/m/s).
     &,RHOKH                                                            &
                                ! OUT Layer k-1 - to - layer k
                                !     turbulent mixing coefficient
                                !     for heat and moisture (kg/m/s).
     &,RHOKM_TOP                                                        &
                                ! OUT exchange coefficient for
                                !     momentum due to top-down mixing
     &,RHOKH_TOP                                                        &
                                ! OUT exchange coefficient for
                                !     heat and moisture due to top-down
                                !     mixing
     &,F_NGSTRESS                                                       &
                                ! OUT dimensionless function for
!                               !     non-gradient stresses
     &,RHOF2                                                            &
                                ! OUT f2 and fsc term shape profiles
     &,RHOFSC                   !       multiplied by rho

      REAL, DIMENSION(row_length,rows) ::                               &
     & ZDSC_BASE                                                        &
                                ! OUT Height of base of K_top in DSC
     &,RHOKH_SURF_ENT                                                   &
                                ! OUT SML surface-driven entrainment KH
     &,RHOKH_TOP_ENT                                                    &
                                ! OUT SML top-driven entrainment KH
     &,RHOKH_DSCT_ENT           ! OUT DSC top-driven entrainment KH

!*L---------------------------------------------------------------------
      EXTERNAL TIMER, excfnl_cci, excfnl_compin
!*L---------------------------------------------------------------------
!    Local and other symbolic constants :-
#include "c_lheat.h"
#include "c_r_cp.h"
#include "c_g.h"
#include "c_epslon.h"
#include "c_vkman.h"
#include "blopt8a.h"
! Options for flux gradient formulation, switch FLUX_GRAD:
!   Locketal2000   = Flux gradients as in Lock et al. (2000)
!   HoltBov1993    = Flux gradients as in Lock et al (2000) but using
!                    coefficients from Holtslag and Boville (1993)
!   LockWhelan2006 = Flux gradients as in Lock and Whelan (2006)

      REAL A_ENT_1,A_ENT_2,C_T,A_ENT_SHR,DEC_THRES_CLEAR,DEC_THRES_CLOUD
      INTEGER N_STEPS
      PARAMETER (                                                       &
     & A_ENT_1=0.23                                                     &
                                ! Entrainment parameter.
     &,A_ENT_2=0.056                                                    &
                                ! Entrainment parameter.
     &,A_ENT_SHR=5.                                                     &
                                ! Entrainment parameter.
     &,C_T=1.0                                                          &
                                ! Parameter in Zilitinkevich term.
     &,DEC_THRES_CLEAR=1.0                                              &
                                ! Decoupling thresholds for clear and
     &,DEC_THRES_CLOUD=0.05                                             &
                                ! cloudy layers (larger makes
!                               ! decoupling less likely)
     &,N_STEPS=3                                                        &
                                ! Number of steps through the mixed
!                               ! layer per sweep
     &)

      REAL S_M,A_NGS
      PARAMETER (                                                       &
     & S_M   = 1.0                                                      &
                                ! empirical parameters in
     &,A_NGS = 2.7                                                      &
                                ! non-gradient stresses
     &)
!*
!
!  Define local storage.
!
!  (a) Workspace.
!
      INTEGER                                                           &
     & KSURF(row_length,rows)   ! First Theta-level above surface layer
                                !  well-mixed SC layer
      LOGICAL, DIMENSION(row_length,rows) ::                            &
     & SCBASE                                                           &
                                ! Flag to signal base of CML reached
     &,TEST_WELL_MIXED                                                  &
                                ! Flag to test wb integration
!                               ! for a well-mixed layer
     &,KSURF_ITERATE                                                    &
                                ! Flag to perform iteration to
!                               ! find top of Ksurf
     &,KTOP_ITERATE             ! Flag to perform iteration to
!                               ! find base of Ktop
      REAL                                                              &
     & KH_SURF(row_length,rows,2:BL_LEVELS)                             &
                                ! Shape factor for non-local
!                               ! turbulent mixing coefficient
     &,WBMIX(row_length,rows,BL_LEVELS)                                 &
                                        ! WB*DZ if were diag as mixed
     &,WBEND(row_length,rows,BL_LEVELS) ! WB*DZ after dec diag

      REAL, DIMENSION(row_length,rows) ::                               &
     & W_M_TOP                                                          &
                      ! Turbulent velocity scale for momentum
!                     !   evaluated at the top of the b.l.
     &,W_H_TOP                                                          &
                      ! Turbulent velocity scale for scalars
!                     !   evaluated at the top of the b.l.
     &,PRANDTL_TOP                                                      &
                      ! Turbulent Prandtl number
                      !   evaluated at the top of the b.l.
     &,KH_TOP_FACTOR                                                    &
                      ! Factors to ensure K_H and K_M profiles are
     &,KM_TOP_FACTOR                                                    &
                      !   continuous at ZH and ZHSC
     &,KH_SCT_FACTOR                                                    &
                      !               "
     &,KM_SCT_FACTOR                                                    &
                      !               "
     &,KH_DSCT_FACTOR                                                   &
                      !               "
     &,KM_DSCT_FACTOR                                                   &
                      !               "
     &,V_TOP                                                            &
                      ! Velocity scale for top-down convection
     &,V_TOP_DSC                                                        &
                      !
     &,V_SUM                                                            &
                      ! total velocity scale
     &,V_SUM_DSC                                                        &
                      !
     &,ZSML_TOP                                                         &
                      ! Height of top of surf-driven K in SML
     &,ZSML_BASE                                                        &
                      ! Height of base of top-driven K in SML
     &,V_SURF                                                           &
                      ! Velocity scale for surface-up conve
     &,SCDEPTH                                                          &
                      ! Depth of top-driven mixing in SML
     &,Z_INV                                                            &
                      ! inversion height (top of K profiles)
     &,WB_SURF_INT                                                      &
                      ! Estimate of wb integrated over surface layer
     &,WB_DZRAD_INT                                                     &
                      ! Estimate of wb integrated over cloud-top region
     &,DZRAD                                                            &
                      ! Depth of cloud-top (radiatively cooled) region
     &,V_KTOP                                                           &
                      ! velocity scale for K_top profile
     &,V_KSUM                                                           &
                      ! total velocity scale
     &,Z_CBASE                                                          &
                      ! cloud base height
     &,WB_RATIO                                                         &
                      ! WBN_INT/WBP_INT
     &,DEC_THRES                                                        &
                      ! Local decoupling threshold
     &,WBP_INT                                                          &
                      ! Positive part of buoyancy flux integral
     &,WBN_INT                                                          &
                      ! Negative part of buoyancy flux integral
     &,ZINV_PR                                                          &
                      ! Height of layer top above surface
     &,KHTOP                                                            &
                      ! temporary KH_top in wb integration
     &,KHSURF                                                           &
                      ! temporary KH_surf in wb integration
     &,ZWB0                                                             &
                      ! height at which wb assumed to go to zero
     &,Z_TOP_LIM                                                        &
                      ! upper height limit on K profile
     &,Z_BOT_LIM                                                        &
                      ! lower height limit on K profile
     &,Z_INC                                                            &
                      ! Step size (m)
     &,RHO_WE                                                           &
                      ! rho*param.d entrainment rates...
     &,RHO_WE_SML                                                       &
                      !  ...for surf and DSC layers (kg/m2/s)
     &,RHO_WE_DSC                                                       &
                      !
     &,CF_ML                                                            &
                      ! Mixed layer cloud fraction
     &,DF_CTOP                                                          &
                      ! Cloud-top radiative flux divergence
     &,DQW                                                              &
                      ! QW jump across inversion
     &,DSL                                                              &
                      ! SL jump across inversion
     &,FT_NT_DSCB                                                       &
                      ! Non-turbulent heat flux at DSC base
     &,FQ_NT_DSCB                                                       &
                      ! Non-turbulent moisture flux at DSC base
     &,TOTHF_ZI                                                         &
                      ! Total heat and moisture fluxes at
     &,TOTQF_ZI       !   the inversion height, Zi
!
!  (b) Scalars.
!
      REAL                                                              &
     & PRANDTL                                                          &
                    ! Turbulent Prandtl number.
     &,PR_NEUT                                                          &
                    ! Neutral limit for Prandtl number
     &,PR_CONV                                                          &
                    ! Convective limit for Prandtl number
     &,ZK_UV                                                            &
                    ! Height above surface of u,v-level.
     &,ZK_TQ                                                            &
                    ! Height above surface of T,q-level.
     &,WSTAR3                                                           &
                    ! Cube of free-convective velocity scale
     &,C_WS                                                             &
                    ! Empirical constant multiplying Wstar
     &,W_S_CUBED_UV                                                     &
                    ! WSTAR for u,v-level
     &,W_S_CUBED_TQ                                                     &
                    !   and T,q-level
     &,W_M_UV                                                           &
                    ! Turbulent velocity scale for momentum: u,v-level
     &,W_M_TQ                                                           &
                    !   and T,q-level
     &,W_H_UV                                                           &
                    ! Turbulent velocity scale for scalars: u,v-level
     &,W_H_TQ                                                           &
                    !   and T,q-level
     &,W_M_HB_3                                                         &
                    ! Cube of W_M, as given by Holtslag and Boville, 93
     &,W_M_NEUT                                                         &
                    ! Neutral limit for W_M
     &,SF_TERM                                                          &
                    ! Surface flux term for entrainment parametrization.
     &,SF_SHEAR_TERM                                                    &
                    ! Surface shear term for entrainment paramn.
     &,IR_TERM                                                          &
                    ! Indirect radiative term for entrainment paramn.
     &,DR_TERM                                                          &
                    ! Direct radiative term for entrainment paramn.
     &,EVAP_TERM                                                        &
                    ! Evaporative term in entrainment parametrization.
     &,ZIL_CORR                                                         &
                    ! Zilitinkevich correction term in entrn. paramn.
     &,ZETA_S_FAC                                                       &
                    ! Factor involving ZETA_S.
     &,ZETA_R_SQ                                                        &
                    ! ZETA_R squared.
     &,ZR                                                               &
                    ! Ratio ZC/ZH.
     &,Z_PR                                                             &
                    ! Height above surface layer
     &,ZH_PR                                                            &
                    ! Height of layer top above surface
     &,Z_RATIO                                                          &
                    ! Ratio of heights
     &,ZCML_BASE                                                        &
                    ! Height of base of cloud mixed layer
     &,RHOKH_ENT                                                        &
                    ! entrainment eddy viscosity
     &,FRAC_TOP                                                         &
                    ! Fraction of turbulent mixing driven from the top
     &,FACTOR                                                           &
                    ! Temporary scalar
     &,ALPHA_T                                                          &
                    ! Parametrized fraction of cloud-top
!                   ! radiative cooling within the inversion
     &,DZ_INV                                                           &
                    ! Parametrizzed inversion thickness (m)
     &,L_RAD                                                            &
                    ! Estimate of e-folding radiative flux
!                   ! decay depth (assumed >= 25m)
     &,WB_CLD                                                           &
                     ! Cloud layer buoyancy flux
     &,WB_SCLD                                                          &
                     ! Sub-cloud layer buoyancy flux
     &,CLD_FRAC                                                         &
                     ! Vertical fraction of layer containing cloud
     &,ZB_KTOP                                                          &
                     ! height of base of K_top profile
     &,DB_RATIO                                                         &
                     ! Temporary in ZWB0 calculation
     &,GAMMA_WBS                                                        &
                     ! Surface layer wb gradient
     &,WSL_DZRAD_INT                                                    &
                     ! Estimate of wsl and wqw integrated over
     &,WQW_DZRAD_INT                                                    &
                     !   the cloud-top region
     &,WSLNG                                                            &
                     ! Non-gradient part of SL flux
     &,WQWNG                                                            &
                     ! Non-gradient part of QW flux
     &,F2, FSc       ! Shape functions for non-gradient fluxes
!
      INTEGER                                                           &
     & I,j                                                              &
                    ! Loop counter (horizontal field index).
     &,K                                                                &
                    ! Loop counter (vertical level index).
     &,N_SWEEP                                                          &
                    ! sweep counter
     &,NS           ! step counter
!
! 2D arrays for optimisation
!
      INTEGER, DIMENSION(row_length,rows) ::                            &
     & ntop                                                             &
                     ! top level of surf-driven K profile
     &,ntml_new                                                         &
                     ! temporary in NTML calculation
     &,kwb0          ! level at which wb assumed to go to zero

      INTEGER, DIMENSION(rows*row_length) ::                            &
     & up            ! indicator of upward/downward sweep

      LOGICAL, DIMENSION(row_length,rows) :: status_ntml

      INTEGER :: ksurf_min, ntop_max, ntml_max, ntdsc_max

! Array introduced to calculate kwb0

      logical, dimension(row_length,rows) :: kstatus
      integer                             :: max_ntml

! Variables for vector compression

      integer :: ij_len
      integer :: ic
      integer :: c_len
      logical,dimension(row_length*rows)  :: to_do
      integer,dimension(row_length*rows)  :: ind_todo
      integer :: c_len_i
      logical,dimension(row_length*rows)  :: todo_inner
      integer,dimension(row_length*rows)  :: ind_todo_i

      integer :: i1, j1, l
!
!     Status from call to GCOM
      INTEGER :: istat_gc

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('EXCF_NL  ',3)
      END IF
!
!-----------------------------------------------------------------------
! Index to subroutine EXCFNL8C
!
! 0. Calculate top-of-b.l. velocity scales and Prandtl number.
! 1. Calculate the top-of-b.l. entrainment parametrization
! 2. Estimate the depths of top-down and surface-up mixing.
!   2.1 First test for well-mixed boundary layer
!   2.2 Iterate to find top of surface-driven mixing, ZSML_TOP,
!   2.3 Iterate to find the base of the top-driven K profile, ZDSC_BASE.
! 3. Calculate factors required to ensure that the K profiles are
!    continuous at the inversion
! 4. Calculate height dependent turbulent transport coefficients
!    within the mixing layers.
!
!-----------------------------------------------------------------------
!! 0.  Calculate top-of-b.l. velocity scales and Prandtl number.
!-----------------------------------------------------------------------
!
!cdir collapse
       do j=1,rows
       do i=1,row_length

         RHOKH_SURF_ENT(I,j) = 0.0
         RHOKH_TOP_ENT(I,j) = 0.0
         RHOKH_DSCT_ENT(I,j) = 0.0
         V_TOP(I,j) = 0.0
         V_SURF(I,j)= 0.0
         V_SUM(I,j) = 0.0
         V_TOP_DSC(I,j) = 0.0
         V_SUM_DSC(I,j) = 0.0
         RHO_WE(I,j)     = 0.0
         RHO_WE_SML(I,j) = 0.0
         RHO_WE_DSC(I,j) = 0.0

         IF (FB_SURF(I,j)  >=  0.0) THEN
!
!         ! Free-convective velocity scale cubed
!
          IF (COUPLED(I,j)) THEN
            WSTAR3 = ZHSC(I,j) * FB_SURF(I,j)
          ELSE
            WSTAR3 =   ZH(I,j) * FB_SURF(I,j)
          END IF

          IF (FLUX_GRAD  ==  Locketal2000) THEN
!
!           ! Turbulent velocity scale for momentum
!
            C_WS = 0.25
            W_M_TOP(I,j) = (V_S(I,j)*V_S(I,j)*V_S(I,j) +                &
     &                      C_WS*WSTAR3)**(1.0/3.0)
!
!           ! Turbulent Prandtl number and velocity scale for scalars
!           ! gives 0.375<Pr<0.75 for convective to neutral conditions
            PRANDTL_TOP(I,j) = 0.75 *                                   &
     &                   ( V_S(I,j)*V_S(I,j)*V_S(I,j)*V_S(I,j) +        &
     &                       (1.0/25.0)*WSTAR3*W_M_TOP(I,j) ) /         &
     &                   ( V_S(I,j)*V_S(I,j)*V_S(I,j)*V_S(I,j) +        &
     &                       (2.0/25.0)*WSTAR3*W_M_TOP(I,j) )
            W_H_TOP(I,j) = W_M_TOP(I,j) / PRANDTL_TOP(I,j)
          ELSEIF (FLUX_GRAD  ==  HoltBov1993) THEN
            C_WS = 0.6
            W_M_TOP(I,j) = (V_S(I,j)*V_S(I,j)*V_S(I,j) +                &
     &                      C_WS*WSTAR3)**(1.0/3.0)
!           ! Using Lock et al interpolation but with
!           ! HB93 range of 0.6<Pr<1
            PR_NEUT = 1.0
            PR_CONV = 0.6
            PRANDTL_TOP(I,j) = PR_NEUT *                                &
     &            ( V_S(I,j)*V_S(I,j)*V_S(I,j)*V_S(I,j) +               &
     &                              (1.0/25.0)*WSTAR3*W_M_TOP(I,j) ) /  &
     &            ( V_S(I,j)*V_S(I,j)*V_S(I,j)*V_S(I,j) +               &
     &                (PR_NEUT/(25.0*PR_CONV))*WSTAR3*W_M_TOP(I,j) )
            W_H_TOP(I,j) = W_M_TOP(I,j) / PRANDTL_TOP(I,j)
          ELSEIF (FLUX_GRAD  ==  LockWhelan2006) THEN
            PR_NEUT = 0.75
            PR_CONV = 0.6
            C_WS    = 0.42   !  ~ 0.75^3
!           ! Slightly contrived notation since really we know W_H_TOP
!           ! and Prandtl range but this makes similarity to
!           ! Lock et al clearer (possibly!)
            W_M_NEUT = ( V_S(I,j)*V_S(I,j)*V_S(I,j) +                   &
     &                   C_WS*WSTAR3 )**(1.0/3.0)
            W_H_TOP(I,j) = W_M_NEUT / PR_NEUT

            PRANDTL_TOP(I,j) = PR_NEUT *                                &
     &            ( V_S(I,j)*V_S(I,j)*V_S(I,j)*V_S(I,j) +               &
     &                (1.0/25.0)*WSTAR3*W_M_NEUT ) /                    &
     &            ( V_S(I,j)*V_S(I,j)*V_S(I,j)*V_S(I,j) +               &
     &                (PR_NEUT/(25.0*PR_CONV))*WSTAR3*W_M_NEUT )

            W_M_TOP(I,j) = W_H_TOP(I,j) * PRANDTL_TOP(I,j)
          END IF

        ELSE
          W_M_TOP(I,j) = V_S(I,j)
          PRANDTL_TOP(I,j) = 0.75
          W_H_TOP(I,j) = W_M_TOP(I,j) / PRANDTL_TOP(I,j)
        END IF
       END DO
      END DO
!
!-----------------------------------------------------------------------
!! 1. Calculate the top-of-b.l. entrainment parametrization
!-----------------------------------------------------------------------
!! 1.1 Initialise 3D arrays
!-----------------------------------------------------------------------
      DO K=2,BL_LEVELS
!cdir collapse
        do j=1,rows
        do i=1,row_length
           RHOKH(I,j,K) = 0.0
           RHOKM(I,j,K) = 0.0
           RHOKH_TOP(I,j,K) = 0.0
           RHOKM_TOP(I,j,K) = 0.0
           RHOF2(I,j,K)  = 0.0
           RHOFSC(I,j,K) = 0.0
           F_NGSTRESS(I,j,K) = 0.0
           KH_SURF(I,j,K) = 0.0
        END DO
        END DO
      END DO
!cdir collapse
      do j=1,rows
      do i=1,row_length
!-----------------------------------------------------------------------
!! 1.2 Calculate top-of-b.l. entrainment mixing coefficients
!!     and store b.l. top quantities for later use.
!-----------------------------------------------------------------------
!      FIRST the top of the SML (if not coupled)
!-----------------------------------------------
        K = NTML(I,j)+1
        IF ( .NOT.COUPLED(I,j) .AND. FB_SURF(I,j)  >=  0.0 ) THEN
!               ! Correction (controlled by STOPWE_SBL in 6B scheme):
!               !   require FB_SURF>0 too
!           !-----------------------------------------------------------
!           ! Calculate the surface buoyancy flux term
!           !-----------------------------------------------------------
            ZETA_S_FAC = (1.0 - ZETA_S(I,j)) * (1.0 - ZETA_S(I,j))
            SF_TERM = A_ENT_1 * MAX ( 0.0 ,                             &
     &                          ( (1.0 - ZETA_S_FAC) * BFLUX_SURF(I,j)  &
     &                            + ZETA_S_FAC * BFLUX_SURF_SAT(I,j) ) )
!           !-----------------------------------------------------------
!           ! Calculate the surface shear term
!           !-----------------------------------------------------------
            SF_SHEAR_TERM =  A_ENT_SHR * V_S(I,j) * V_S(I,j) * V_S(I,j) &
     &                      * RHO_UV(I,j,K)  / ZH(I,j)
!           !-----------------------------------------------------------
!           ! Calculate the indirect radiative term
!           !-----------------------------------------------------------
            ZETA_R_SQ = ZETA_R(I,j)*ZETA_R(I,j)
            IR_TERM = ( BT_TOP(I,j)*ZETA_R_SQ +                         &
     &                  BTT_TOP(I,j)*(1.0-ZETA_R_SQ) )                  &
     &                * A_ENT_1 * DF_TOP_OVER_CP(I,j)
!           !-----------------------------------------------------------
!           ! Calculate the evaporative term
!           !-----------------------------------------------------------
            IF ( DB_TOP(I,j)  >   0.0) THEN
              ZR = SQRT( ZC(I,j) / ZH(I,j) )
              EVAP_TERM = A_ENT_2 * RHO_UV(I,j,K)                       &
     &                  * CHI_S_TOP(I,j) * CHI_S_TOP(I,j)               &
     &                  * ZR * ZR * ZR * DB_TOP_CLD(I,j)                &
     &                  * SQRT( ZH(I,j) * DB_TOP(I,j) )
            ELSE
              EVAP_TERM = 0.0
            END IF
!           !-----------------------------------------------------------
!           ! Combine forcing terms to calculate the representative
!           ! velocity scales
!           !-----------------------------------------------------------
            V_SUM(I,j) = ( (SF_TERM + SF_SHEAR_TERM +                   &
     &                      IR_TERM + EVAP_TERM)                        &
     &                   * ZH(I,j) /(A_ENT_1*RHO_UV(I,j,K)) )**(1.0/3.0)
            V_TOP(I,j) = ( (IR_TERM+EVAP_TERM) * ZH(I,j)                &
     &                           / (A_ENT_1*RHO_UV(I,j,K)) )**(1.0/3.0)
            V_SURF(I,j) = ( (SF_TERM) * ZH(I,j)                         &
     &                           / (A_ENT_1*RHO_UV(I,j,K)) )**(1.0/3.0)
!           !-----------------------------------------------------------
!           ! Calculate the direct radiative term
!           !  can only calculate for DB_TOP > 0
!           !-----------------------------------------------------------
            IF ( DB_TOP(I,j)  >   0.0) THEN
             DZ_INV  = MIN( V_SUM(I,j)*V_SUM(I,j) / DB_TOP(I,j) ,100.0 )
             L_RAD   = 15.0 * MAX( 1.0 , 200./(ZC(I,j)+1.0E-14) )
             ALPHA_T = 1.0 - EXP(-0.5*DZ_INV/L_RAD)
             IF (D_SIEMS(I,j)  >   0.0) ALPHA_T =                       &
     &             MIN( 1.0, ALPHA_T + 10.0*D_SIEMS(I,j)*(1.0-ALPHA_T) )
             DR_TERM = BTC_TOP(I,j) * ALPHA_T * DF_TOP_OVER_CP(I,j)
!            !----------------------------------------------------------
!            ! Combine terms to calculate the entrainment
!            ! mixing coefficients
!            !----------------------------------------------------------
             ZIL_CORR = C_T * ( (SF_TERM + SF_SHEAR_TERM +              &
     &                           IR_TERM + EVAP_TERM) /                 &
     &                   (RHO_UV(I,j,K) * SQRT(ZH(I,j))) )**(2.0/3.0)

             RHO_WE_SML(I,j) = (SF_TERM + SF_SHEAR_TERM                 &
     &                   +      IR_TERM + EVAP_TERM + DR_TERM)          &
     &                      / ( DB_TOP(I,j) + ZIL_CORR )

             RHOKH_ENT = RHO_WE_SML(I,j)/ RDZ(I,j,K)

             FRAC_TOP = V_TOP(I,j) / ( V_TOP(I,j)+W_H_TOP(I,j)+1.0E-14 )
             RHOKH_SURF_ENT(I,j) = RHOKH_ENT * ( 1.0 - FRAC_TOP )
             RHOKH_TOP_ENT(I,j) = RHOKH_ENT * FRAC_TOP

             RHOKM(I,j,K) = PRANDTL_TOP(I,j) * RHOKH_SURF_ENT(I,j)      &
     &                    * RDZ(I,j,K) * (Z_UV(I,j,K)-Z_UV(I,j,K-1))    &
     &                    * RHO_TQ(I,j,K-1) / RHO_UV(I,j,K)
             RHOKM_TOP(I,j,K) = 0.75 * RHOKH_TOP_ENT(I,j)               &
     &                    * RDZ(I,j,K) * (Z_UV(I,j,K)-Z_UV(I,j,K-1))    &
     &                    * RHO_TQ(I,j,K-1) / RHO_UV(I,j,K)
            END IF    ! test on DB_TOP GT 0
        END IF
!----------------------------------------------------------------
!      THEN the top of the DSC (if coupled use ZHSC length-scale)
!----------------------------------------------------------------
        IF ( NTDSC(I,j)  >   0 ) THEN
          K = NTDSC(I,j)+1
          IF (COUPLED(I,j)) THEN
!             !--------------------------------------------------------
!             ! Calculate the surface buoyancy flux term
!             !--------------------------------------------------------
              ZETA_S_FAC = (1.0 - ZETA_S(I,j)) * (1.0 - ZETA_S(I,j))
              SF_TERM = A_ENT_1 * MAX ( 0.0 ,                           &
     &                          ( (1.0 - ZETA_S_FAC) * BFLUX_SURF(I,j)  &
     &                            + ZETA_S_FAC * BFLUX_SURF_SAT(I,j) ) )
!             !--------------------------------------------------------
!             ! Calculate the surface shear term
!             !--------------------------------------------------------
              IF (FB_SURF(I,j)  >=  0.0) THEN
                SF_SHEAR_TERM = A_ENT_SHR * V_S(I,j)*V_S(I,j)*V_S(I,j)  &
     &                          * RHO_UV(I,j,K)  / ZHSC(I,j)
              ELSE
                SF_SHEAR_TERM = 0.0
              END IF
              V_SURF(I,j) = ( (SF_TERM) * ZHSC(I,j)                     &
     &                            / (A_ENT_1*RHO_UV(I,j,K)) )**(1.0/3.0)
          ELSE
              SF_TERM = 0.0
              SF_SHEAR_TERM = 0.0
          END IF
!         !-----------------------------------------------------------
!         ! Calculate the indirect radiative term
!         !-----------------------------------------------------------
          ZETA_R_SQ = ZETA_R_DSC(I,j)*ZETA_R_DSC(I,j)
          IR_TERM = ( BT_DSCT(I,j)*ZETA_R_SQ +                          &
     &                  BTT_DSCT(I,j)*(1.0-ZETA_R_SQ) )                 &
     &                  * A_ENT_1 * DF_DSCT_OVER_CP(I,j)
!         !-----------------------------------------------------------
!         ! Calculate the evaporative term
!         !-----------------------------------------------------------
          IF (DB_DSCT(I,j)  >   0.0) THEN
              ZR = SQRT( ZC_DSC(I,j) / DSCDEPTH(I,j) )
              EVAP_TERM = A_ENT_2 * RHO_UV(I,j,K)                       &
     &                  * CHI_S_DSCT(I,j) * CHI_S_DSCT(I,j)             &
     &                  * ZR * ZR * ZR * DB_DSCT_CLD(I,j)               &
     &                  * SQRT( DSCDEPTH(I,j) * DB_DSCT(I,j) )
          ELSE
              EVAP_TERM = 0.0
          END IF
!         !-----------------------------------------------------------
!         ! Combine forcing terms to calculate the representative
!         ! velocity scales
!         !-----------------------------------------------------------
          V_SUM_DSC(I,j) = ( (SF_TERM + SF_SHEAR_TERM +                 &
     &                          IR_TERM + EVAP_TERM)                    &
     &              * DSCDEPTH(I,j) / (A_ENT_1*RHO_UV(I,j,K)) )**(1./3.)
          V_TOP_DSC(I,j) =( (IR_TERM + EVAP_TERM) * DSCDEPTH(I,j) /     &
     &                         (A_ENT_1*RHO_UV(I,j,K)) )**(1.0/3.0)
!         !-----------------------------------------------------------
!         ! Calculate the direct radiative term
!         !-----------------------------------------------------------
          IF (DB_DSCT(I,j)  >   0.0) THEN
             DZ_INV  = MIN( V_SUM_DSC(I,j)*V_SUM_DSC(I,j)/DB_DSCT(I,j), &
     &                      100.0 )
             L_RAD   = 15.0 * MAX( 1.0 , 200./(ZC_DSC(I,j)+1.0) )
             ALPHA_T = 1.0 - EXP(-0.5*DZ_INV/L_RAD)
             IF (D_SIEMS_DSC(I,j)  >   0.0) ALPHA_T =                   &
     &         MIN( 1.0, ALPHA_T + 10.0*D_SIEMS_DSC(I,j)*(1.0-ALPHA_T) )
             DR_TERM = BTC_TOP(I,j) * ALPHA_T * DF_DSCT_OVER_CP(I,j)
!            !----------------------------------------------------------
!            ! Finally combine terms to calculate the entrainment
!            ! rate and mixing coefficients
!            !----------------------------------------------------------
             ZIL_CORR = C_T * ( (SF_TERM + SF_SHEAR_TERM +              &
     &                           IR_TERM + EVAP_TERM) /                 &
     &                (RHO_UV(I,j,K) * SQRT(DSCDEPTH(I,j))) )**(2.0/3.0)
             RHO_WE_DSC(I,j) = ( SF_TERM + SF_SHEAR_TERM                &
     &               +           IR_TERM + EVAP_TERM + DR_TERM )        &
     &                          / ( DB_DSCT(I,j) + ZIL_CORR )

             RHOKH_DSCT_ENT(I,j) = RHO_WE_DSC(I,j)/ RDZ(I,j,K)
             RHOKM_TOP(I,j,K) = 0.75 * RHOKH_DSCT_ENT(I,j)              &
     &                    * RDZ(I,j,K) * (Z_UV(I,j,K)-Z_UV(I,j,K-1))    &
     &                    * RHO_TQ(I,j,K-1) / RHO_UV(I,j,K)
          END IF   ! test on DB_DSCT gt 0
        END IF
      END DO
      END DO
!
!  If there is no turbulence generation in DSC layer, ignore it.
!
!cdir collapse
      DO j=1,rows
      DO I=1,row_length
        IF ( V_TOP_DSC(I,j)  <=  0.0 ) THEN
          DSC(I,j) = .FALSE.
          NTDSC(I,j) = 0
          ZHSC(I,j) = 0.0
          ZC_DSC(I,j) = 0.0
          DSCDEPTH(I,j) = 0.0
          COUPLED(I,j) = .FALSE.
        END IF
      END DO
      END DO
! ----------------------------------------------------------------------
! 2.0 Estimate the depths of top-down and surface-up mixing.
!     These amount to diagnoses of recoupling and decoupling.
!     The K_H profiles are applied over layers such that the ratio
!        WBN_INT/WBP_INT = DEC_THRES (parameter),
!     where WBN_INT and WBP_INT are the magnitudes of the integrals of
!     the negative and positive parts, respectively, of the resulting
!     buoyancy flux profile (given by - KH * DB_FOR_FLUX).
! ----------------------------------------------------------------------
! 2.1 First test for well-mixed boundary layer
!     (ie. both KH profiles extending from cloud-top to the surface).
!     If the parcel ascent diagnosed:
!        DSC    - test for well-mixed up to ZHSC = recoupling
!        no DSC - test for well-mixed up to ZH   = decoupling
! -----------------------------------------------------------
      DO K=1,BL_LEVELS
!cdir collapse
      DO j=1,rows
      DO I=1,row_length
!           ! potentially useful diagnostics
            WBMIX(I,j,K)=0.0  ! WB if were diag as well-mixed
            WBEND(I,j,K)=0.0  ! WB after dec diag
      END DO
      END DO
      END DO
!
! Default settings
!
!cdir collapse
      DO j=1,rows
      DO I=1,row_length
        ZSML_TOP(I,j)  = ZH(I,j)
        ZSML_BASE(I,j) = 0.1 * ZH(I,j)
        ZDSC_BASE(I,j) = 0.1 * ZHSC(I,j)
        DEC_THRES(I,j) = DEC_THRES_CLOUD  ! use cloudy by default
        Z_INV(i,j) = 0.0    ! inversion height (top of K profiles)
        ZWB0(i,j)  = 0.0    ! height at which WB goes to zero
        WBP_INT(i,j) = 0.0
        WBN_INT(i,j) = 0.0
        WB_SURF_INT(i,j) = 0.0
        WB_DZRAD_INT(i,j) = 0.0
        DZRAD(I,j) = 100.0
        TOTHF_ZI(I,j) = 0.0
        TOTQF_ZI(I,j) = 0.0
        kstatus(i,j)= .true.
        kwb0(i,j)  = 2
        NTOP(i,j)  = -1
        KSURF(i,j) = 1
      END DO
      END DO
!
! Find KSURF, the first theta-level above the surface layer
!
      DO K=2,BL_LEVELS
!cdir collapse
      DO j=1,rows
      DO I=1,row_length
          IF ( Z_TQ(I,j,K-1)  <   0.1*ZH(I,j) ) KSURF(I,j) = K
      END DO
      END DO
      END DO
!
! Find kwb0, level with lowest positive cloud-free buoyancy gradient
!
      max_ntml=maxval(ntml)
#if !defined(SCMA)
!     For bit-reproducibility make consistent 
!     for all PE decompositions.
      CALL gc_imax(1, n_proc, istat_gc, max_ntml)
#endif
      do k=2,max_ntml
        do j=1, rows
        do i=1,row_length
          if (kstatus(i,j)) then
            if ( (db_ga_dry(i,J,k) <=  0.0) .or.                        &
     &           (k >= ntml(i,j)) ) then
                       kstatus(i,j)=.false.
                       kwb0(i,j)=k
            end if
          end if
        end do
        end do
      end do
!
! Set flags for iterating wb integral to calculate depth of mixing,
! one each for KSURF and K_TOP.  Note these will be updated depending on
! what happpens on testing for a well-mixed layer in section 2.2.
!
!cdir collapse
      DO j=1,rows
      DO I=1,row_length
        TEST_WELL_MIXED(I,j) = .FALSE.
        KSURF_ITERATE(I,j)= .FALSE.
        KTOP_ITERATE(I,j) = .FALSE.
!       ! Following test not needed!!
!        IF ( BFLUX_SURF(I,j)  >   0.0 .AND.
!            ! Otherwise: surface mixing generated by stable scheme
!     &       .NOT.CUMULUS(I,j) .AND.
!            ! Rule out CUMULUS layers from iteration of ZH
!     &       NTDSC(I,j)  >   2
!            ! Otherwise: layer too shallow to resolve decoupling
!     &     ) THEN
!          KSURF_ITERATE(I,j)= .TRUE.
!        END IF
        IF ( NTDSC(I,j)  >   2 ) THEN
          KTOP_ITERATE(I,j) = .TRUE.
        END IF
      END DO ! I
      END DO ! J
!------------------------------------------------------------
! First test buoyancy flux integral for well-mixed layer
! -----------------------------------------------------------
!cdir collapse
      DO j=1,rows
      DO I=1,row_length
        IF ( BFLUX_SURF(I,j)  >   0.0) THEN
!         ! can only be coupled to an unstable SML
          IF ( NTDSC(I,j)  >   2 ) THEN
!           ! well-resolved DSC layer (6B DECFIX included)
!           ! ...test for recoupling with SML
            TEST_WELL_MIXED(I,j) = .TRUE.
            Z_INV(i,j)  = ZHSC(I,j)
            Z_CBASE(i,j)= Z_INV(i,j) - ZC_DSC(I,j)
            CF_ML(I,j)  = CF_DSC(I,j)
            V_KTOP(i,j) = V_TOP_DSC(I,j)
            V_KSUM(i,j) = V_SUM_DSC(I,j)
            DF_CTOP(I,j)= DF_DSCT_OVER_CP(I,j)
            RHO_WE(I,j) = RHO_WE_DSC(I,j)
            DSL(I,j)    = DSL_DSC(I,j)
            DQW(I,j)    = DQW_DSC(I,j)
            TOTHF_ZI(I,j)= - RHO_WE(I,j)*DSL(I,j) + FT_NT_ZHSC(I,j)
            TOTQF_ZI(I,j)= - RHO_WE(I,j)*DQW(I,j) + FQ_NT_ZHSC(I,j)
            NTOP(i,j)   = NTDSC(I,j) - 1
!           ! assuming wb goes to zero by the lowest of ZH
!           ! or cloud-base, but above surface layer
          ELSE IF ( .NOT.DSC(I,j) .AND. .NOT.CUMULUS(I,j) .AND.         &
     &              NTML(I,j)  >   2) THEN
!           ! well-resolved SML       (6B DECFIX included)
!           ! ...test for decoupling
!           ! Note: code can only deal with one DSC layer at a time so
!           ! can't decouple SML if a DSC layer already exists.
!           ! ---------------------------------------------------------
!           ! 6B DECFIX switch included:
!           ! If the BL layer is cloud-free then use a less restrictive
!           ! threshold - ideally, the parcel ascent would have
!           ! found the correct BL top in this case but this test is
!           ! kept to keep negative buoyancy fluxes under control
!           ! (ie. DEC_THRES_CLEAR=1 ensures wbn_int < |wbp_int|)
            IF (ZC(I,j)  ==  0.0) DEC_THRES(I,j) = DEC_THRES_CLEAR
!           ! ---------------------------------------------------------
            TEST_WELL_MIXED(I,j) = .TRUE.
            Z_INV(i,j)  = ZH(I,j)
            Z_CBASE(i,j)= Z_INV(i,j) - ZC(I,j)
            CF_ML(I,j)  = CF_SML(I,j)
            V_KTOP(i,j) = V_TOP(I,j)
            V_KSUM(i,j) = V_SUM(I,j)
            DF_CTOP(I,j)= DF_TOP_OVER_CP(I,j)
            RHO_WE(I,j) = RHO_WE_SML(I,j)
            DSL(I,j)    = DSL_SML(I,j)
            DQW(I,j)    = DQW_SML(I,j)
            TOTHF_ZI(I,j)= - RHO_WE(I,j)*DSL(I,j) + FT_NT_ZH(I,j)
            TOTQF_ZI(I,j)= - RHO_WE(I,j)*DQW(I,j) + FQ_NT_ZH(I,j)
            NTOP(i,j)   = NTML(I,j) - 1
          END IF
        END IF
      END DO ! I
      END DO ! J
!cdir collapse
      DO j=1,rows
      DO I=1,row_length
        IF ( TEST_WELL_MIXED(I,j) ) THEN
!         !-------------------------------------------------------------
!         ???  This code also invoked for cloud-free cases
!               - does this matter?
!               - only ignoring top grid level of integration which
!                 shouldn't be important/relevant for cloud-free layers.
!         ! Estimate wb integral over radiatively cooled cloud-top
!         ! region, from Z_INV to Z_INV-DZRAD.
!         ! DZRAD taken to be constant (100m) for simplicity, but also
!         ! integration depth taken to extend down to at least
!         ! Z_TQ(NTML-1) but without going below cloud-base.
!         !-------------------------------------------------------------
          DZRAD(I,j) = 100.0
          DO WHILE ( Z_TQ(I,j,NTOP(I,j))  >   Z_INV(I,j)-DZRAD(I,j)     &
     &         .AND. Z_TQ(I,j,NTOP(I,j)-1)  >   Z_INV(I,j)-Z_CBASE(I,j) &
     &         .AND. NTOP(I,j)  >   2 )
            NTOP(i,j) = NTOP(i,j) - 1
          END DO
          DZRAD(I,j) = Z_INV(I,j) - Z_TQ(I,j,NTOP(I,j))

          WSL_DZRAD_INT = DZRAD(I,j) *                                  &
     &                    ( 0.66*DF_CTOP(I,j) - RHO_WE(I,j)*DSL(I,j) )
          WQW_DZRAD_INT = - DZRAD(I,j) * RHO_WE(I,j) * DQW(I,j)

          WB_DZRAD_INT(I,j) = WSL_DZRAD_INT * (                         &
     &                          (1.0-CF_ML(I,j))*BTM(I,j,NTOP(i,j)+1) + &
     &                            CF_ML(I,j)*BTM_CLD(I,j,NTOP(i,j)+1) ) &
     &                      + WQW_DZRAD_INT * (                         &
     &                          (1.0-CF_ML(I,j))*BQM(I,j,NTOP(i,j)+1) + &
     &                            CF_ML(I,j)*BQM_CLD(I,j,NTOP(i,j)+1) )
          WB_DZRAD_INT(I,j) = G * WB_DZRAD_INT(I,j)
          WB_DZRAD_INT(I,j) = MAX( 0.0, WB_DZRAD_INT(I,j) )

!         ! Include WB_DZRAD_INT in WBP_INT as it set to be >0
          WBP_INT(i,j) = WBP_INT(i,j) + WB_DZRAD_INT(I,j)

        ELSE

          WB_DZRAD_INT(I,j) = -1.0  ! To identify not calculated

        END IF
      END DO ! I
      END DO ! J
!
! For WB diagnostics, convert integrated WB to uniform profile
!
      DO K = 2, BL_LEVELS
        DO j=1,rows
        DO I=1,row_length
          IF ( TEST_WELL_MIXED(I,j) .AND. K >=  NTOP(I,j)+1 ) THEN
            IF ( K <= NTDSC(I,j)+1 ) THEN
              WBEND(I,j,K) = WB_DZRAD_INT(I,j)/DZRAD(I,j)
              WBMIX(I,j,K) = WB_DZRAD_INT(I,j)/DZRAD(I,j)
            ELSE IF ( K <= NTML(I,j)+1 ) THEN
              WBEND(I,j,K) = WB_DZRAD_INT(I,j)/DZRAD(I,j)
              WBMIX(I,j,K) = WB_DZRAD_INT(I,j)/DZRAD(I,j)
            END IF
          END IF
        END DO ! I
        END DO ! J
      END DO ! K

!cdir collapse
      DO j=1,rows
      DO I=1,row_length
        IF ( TEST_WELL_MIXED(I,j) ) THEN
!         !----------------------------------------------------
!         ! Estimate wb integral over surface layer
!         ! (and up to next theta-level, namely Z_TQ(KSURF) )
!         ! assuming linear profile going to zero at ZWB0
!         !----------------------------------------------------
          IF ( Kwb0(i,j)  ==  NTML(I,j) ) THEN
            ZWB0(i,j) = ZH(I,j)
          ELSEIF ( Kwb0(i,j)  ==  2 ) THEN
            ZWB0(i,j) = Z_UV(I,j,2)
          ELSE
            K=Kwb0(i,j)
!           ! now DB_GA_DRY(K) LE 0 and DB_GA_DRY(K-1) GT 0
!           ! so interpolate:
            DB_RATIO = DB_GA_DRY(I,j,K-1)                               &
     &              / ( DB_GA_DRY(I,j,K-1) - DB_GA_DRY(I,j,K) )
            DB_RATIO = MAX( 0.0, DB_RATIO )  ! trap for rounding error
            ZWB0(i,j)=Z_UV(I,j,K-1) +                                   &
     &                DB_RATIO * (Z_UV(I,j,K)-Z_UV(I,j,K-1))
          END IF
          WB_SURF_INT(i,j) = BFLUX_SURF(I,j) * Z_TQ(I,j,KSURF(I,j)) *   &
     &                  ( 1.0 - Z_TQ(I,j,KSURF(I,j))/(2.0*ZWB0(i,j)))
          WB_SURF_INT(i,j) = MAX( 1.0E-14, WB_SURF_INT(i,j) )
        ELSE
!         ! only include surface layer contribution for unstable mixing
          WB_SURF_INT(i,j) = 1.0E-14
        END IF

        WBP_INT(i,j) = WBP_INT(i,j) + WB_SURF_INT(i,j) ! must be >0

!       ! Save surface and bl-top layer integral for diagnostics
        WBMIX(I,j,KSURF(I,j)) = WB_SURF_INT(i,j)
        WBEND(I,j,KSURF(I,j)) = WB_SURF_INT(i,j)

      END DO ! I
      END DO ! J

      ksurf_min=minval(KSURF)
#if !defined(SCMA)
!     For bit-reproducibility make consistent
!     for all PE decompositions.
      CALL gc_imin(1, n_proc, istat_gc, ksurf_min)
#endif
      ntop_max=maxval(NTOP)
#if !defined(SCMA)
!     For bit-reproducibility make consistent
!     for all PE decompositions.
      CALL gc_imax(1, n_proc, istat_gc, ntop_max)
#endif

      DO K=ksurf_min+1,ntop_max
        DO J=1, rows
        DO I=1, row_length

         IF ( TEST_WELL_MIXED(I,j) ) THEN
!         ! ----------------------------------------------
!         ! worth testing layer as well-mixed to cloud-top
!         ! ----------------------------------------------
          ZB_KTOP = 0.1*Z_INV(i,j)
          ZINV_PR(i,j) = Z_INV(i,j) - ZB_KTOP
!         ! DB(K)is the K to K-1 difference and already
!         ! integrated up to K_SURF, so start this loop at KSURF+1
          IF ( (K >= KSURF(I,j)+1) .AND. (K <= NTOP(i,j)) ) THEN

            KHTOP(i,j) = 0.0
            KHSURF(i,j)= 0.0
            F2         = 0.0
            FSc        = 0.0

            Z_PR = Z_UV(I,j,K) - ZB_KTOP
            IF (Z_PR  >   0.0 .AND. Z_PR  <   ZINV_PR(i,j)) THEN
              Z_RATIO = Z_PR/ZINV_PR(i,j)

              IF (FLUX_GRAD  ==  LockWhelan2006) THEN
                KHTOP(i,j) = 3.6 * VKMAN * RHO_UV(I,j,K) * V_KTOP(I,j)  &
     &                           * ZINV_PR(i,j) * (Z_RATIO**3.0)        &
     &                           * ( (1.0-Z_RATIO)*(1.0-Z_RATIO) )
                F2 = RHO_UV(I,j,K) * 0.5 * Z_RATIO                      &
     &                                   * 2.0**( Z_RATIO**4.0 )
                IF ( V_KSUM(i,j)  >   0.0 ) THEN
                  FSc = RHO_UV(I,j,K) * 3.5 * (V_KTOP(i,j)/V_KSUM(i,j)) &
     &                  * (Z_RATIO**3.0) * (1.0-Z_RATIO)
                END IF
              ELSE
                KHTOP(i,j) = 0.85 * VKMAN * RHO_UV(I,j,K) * V_KTOP(i,j) &
     &                            * (( 1.0 - Z_RATIO )**0.8)            &
     &                            * Z_PR * Z_RATIO
              END IF
            END IF

            Z_PR = Z_UV(I,j,K)
            IF ( Z_PR  <   Z_INV(i,j)) THEN
!             !--------------------------------
!             ! include surface-driven profile
!             !--------------------------------
              KHSURF(i,j) = VKMAN * RHO_UV(I,j,K) *                     &
     &                 W_H_TOP(I,j)*Z_PR*( 1.0 - Z_PR/Z_INV(i,j) )      &
     &                                  *( 1.0 - Z_PR/Z_INV(i,j) )
            END IF

            IF (FLUX_GRAD  ==  LockWhelan2006) THEN
              WSLNG = (F2+FSc)*TOTHF_ZI(i,j) - FT_NT(i,j,K)
              WQWNG = (F2+FSc)*TOTQF_ZI(i,j) - FQ_NT(i,j,K)
            END IF

            IF ( Z_TQ(I,j,K)  <=  Z_CBASE(i,j) ) THEN
!             ! Completely below cloud-base so use cloud-free formula
              WB_SCLD = KHSURF(i,j) * DB_GA_DRY(I,j,K) +                &
     &                  KHTOP(i,j) * DB_NOGA_DRY(I,j,K)
              IF (FLUX_GRAD  ==  LockWhelan2006) THEN
                WB_SCLD = WB_SCLD + ( G/RDZ(I,j,K) ) *                  &
     &            ( BTM(I,j,K-1)*WSLNG + BQM(I,j,K-1)*WQWNG )
              END IF
              WB_CLD  = 0.0
            ELSEIF (Z_TQ(I,j,K-1)  >=  Z_CBASE(i,j)) THEN
!             ! Completely above cloud-base so use cloudy formula
              WB_CLD = ( KHSURF(i,j) * DB_GA_CLD(I,j,K) +               &
     &                   KHTOP(i,j)  * DB_NOGA_CLD(I,j,K) )
              IF (FLUX_GRAD  ==  LockWhelan2006) THEN
                WB_CLD = WB_CLD + ( G/RDZ(I,j,K) ) * (                  &
     &                   ( BTM(I,j,K-1)*(1.0-CF_ML(I,j)) +              &
     &                     BTM_CLD(i,j,K-1)*CF_ML(i,j) )*WSLNG +        &
     &                   ( BQM(i,j,K-1)*(1.0-CF_ML(i,j)) +              &
     &                     BQM_CLD(i,j,K-1)*CF_ML(i,j) )*WQWNG )
              END IF
              WB_SCLD = 0.0
            ELSE
!             ! cloud-base within this integration range
!             ! so treat cloud and sub-cloud layer wb separately
              WB_SCLD = KHSURF(i,j) * DB_GA_DRY(I,j,K) +                &
     &                  KHTOP(i,j) * DB_NOGA_DRY(I,j,K)
              WB_CLD = ( KHSURF(i,j) * DB_GA_CLD(I,j,K) +               &
     &                   KHTOP(i,j)  * DB_NOGA_CLD(I,j,K) )
              IF (FLUX_GRAD  ==  LockWhelan2006) THEN
                WB_SCLD = WB_SCLD + ( G/RDZ(I,j,K) ) *                  &
     &            ( BTM(I,j,K-1)*WSLNG + BQM(I,j,K-1)*WQWNG )
                WB_CLD = WB_CLD + ( G/RDZ(I,j,K) ) * (                  &
     &                   ( BTM(I,j,K-1)*(1.0-CF_ML(I,j)) +              &
     &                     BTM_CLD(i,j,K-1)*CF_ML(i,j) )*WSLNG +        &
     &                   ( BQM(i,j,K-1)*(1.0-CF_ML(i,j)) +              &
     &                     BQM_CLD(i,j,K-1)*CF_ML(i,j) )*WQWNG )
              END IF
              CLD_FRAC = (Z_TQ(I,j,K)-Z_CBASE(i,j))                     &
     &                  /(Z_TQ(I,j,K)-Z_TQ(I,j,K-1))
              WB_CLD  = CLD_FRAC * WB_CLD
              WB_SCLD = (1.0-CLD_FRAC) * WB_SCLD
            END IF

            WBMIX(I,j,K) = WB_CLD+WB_SCLD
            WBEND(I,j,K) = WB_CLD+WB_SCLD

            IF (WB_CLD  >=  0.0) THEN
              WBP_INT(i,j) = WBP_INT(i,j) + WB_CLD
            ELSE
              WBN_INT(i,j) = WBN_INT(i,j) - WB_CLD
            END IF
            IF (WB_SCLD  >=  0.0) THEN
              WBP_INT(i,j) = WBP_INT(i,j) + WB_SCLD
            ELSE
              WBN_INT(i,j) = WBN_INT(i,j) - WB_SCLD
            END IF

          END IF ! K

          WB_RATIO(i,j) = WBN_INT(i,j)/WBP_INT(i,j)

         END IF ! TEST_WELL_MIXED
        END DO ! I
        END DO ! J
      END DO ! K
!
! Test WB_Ratio to see if well-mixed layer allowed
!
!cdir collapse
      DO J=1,rows
      DO I=1,row_length
        IF ( TEST_WELL_MIXED(I,j) ) THEN
          IF ( WB_RATIO(i,j)  <=  DEC_THRES(i,j) ) THEN
!           ! No need to test depth of mixing any further as
!           ! well-mixed layer buoyancy flux integral criteria.
!           ! SML will simply stay well-mixed (and so use defaults)
            KSURF_ITERATE(I,j)= .FALSE.
            KTOP_ITERATE(I,j) = .FALSE.
            IF ( DSC(I,j) ) THEN
!             ! Recouple DSC with SML:
!             ! move surface driven entrainment
!             ! RHOKH(z_i) = rho * w_e * DZL and w_e ~ 1/DB_TOP, so:
              IF ( DB_TOP(I,j) >  0.0 .AND. DB_DSCT(I,j) >  0.01 ) THEN
!                                         ! can't calc Zil. term
                RHOKH_SURF_ENT(I,j) = RHOKH_SURF_ENT(I,j) *             &
     &                 ( RHO_UV(I,j,NTDSC(I,j)+1) * DB_TOP(I,j) *       &
     &                   RDZ(I,j,NTML(I,j)+1) ) /                       &
     &                 ( RHO_UV(I,j,NTML(I,j)+1) * DB_DSCT(I,j) *       &
     &                                      RDZ(I,j,NTDSC(I,j)+1) )
                RHOKM(I,j,NTDSC(I,j)+1) = RHOKM(I,j,NTML(I,j)+1) *      &
     &                 ( RHO_TQ(I,j,NTDSC(I,j)) * DB_TOP(I,j) *         &
     &                   RDZ(I,j,NTML(I,j)+1) ) /                       &
     &                 ( RHO_TQ(I,j,NTML(I,j)) * DB_DSCT(I,j) *         &
     &                                    RDZ(I,j,NTDSC(I,j)+1) )
              END IF
!             ! redesignate top-driven entrainment at ZHSC
!             ! (ignore that calculated at ZH)
              RHOKH_TOP_ENT(I,j) = RHOKH_DSCT_ENT(I,j)
              ZH(I,j) = ZHSC(I,j)
              NTML(I,j) = NTDSC(I,j)
              V_TOP(I,j) = V_TOP_DSC(I,j)
              ZSML_BASE(I,j) = 0.1 * ZH(I,j)
              ZC(I,j) = ZC_DSC(I,j)
              ZHSC(I,j) = 0.0
              NTDSC(I,j) = 0
              V_TOP_DSC(I,j) = 0.0
              ZDSC_BASE(I,j) = 0.0
              ZC_DSC(I,j)    = 0.0
              FT_NT_ZH(I,j)   = FT_NT_ZHSC(I,j)
              FT_NT_ZHSC(I,j) = 0.0
              FQ_NT_ZH(I,j)   = FQ_NT_ZHSC(I,j)
              FQ_NT_ZHSC(I,j) = 0.0
              DSC(I,j) = .FALSE.
              CUMULUS(I,j) = .FALSE.
              COUPLED(I,j) = .FALSE.
            END IF  ! recoupled DSC layer
          ELSE   ! buoyancy flux threshold violated
!           !---------------------------------
!           ! Extent of mixing must be reduced
!           !---------------------------------
            IF ( .NOT.CUMULUS(I,j) ) KSURF_ITERATE(I,j) = .TRUE.
            KTOP_ITERATE(I,j)  = .TRUE.
            IF (.NOT.DSC(I,j)) THEN
!             ! Set up a `COUPLED' decoupled layer,
!             !   implies no explicit `entrainment' at ZH.
!             ! Note a new ZH (and thence NTML) will be calculated by
!             ! wb integral iteration.
              IF (CUMULUS(I,j)) ZK_UV=SQRT(ZH(I,j)-1000000.)
!                             ! APLTEST: shouldn't ever happen!
              DSC(I,j) = .TRUE.
              COUPLED(I,j) = .TRUE.
              NTDSC(I,j) = NTML(I,j)
              ZHSC(I,j) = ZH(I,j)
              ZC_DSC(I,j) = ZC(I,j)
              V_TOP_DSC(I,j) = V_TOP(I,j)
              V_SUM_DSC(I,j) = V_SUM(I,j)
              FT_NT_ZHSC(I,j) = FT_NT_ZH(I,j)
              FQ_NT_ZHSC(I,j) = FQ_NT_ZH(I,j)
!             ! put all entrainment into RHOKH_TOP
              RHOKH_DSCT_ENT(I,j) = RHOKH_TOP_ENT(I,j)                  &
     &                            + RHOKH_SURF_ENT(I,j)
              RHOKH_TOP_ENT(I,j) = 0.0
              RHOKH_SURF_ENT(I,j) = 0.0
              RHOKM_TOP(I,j,NTML(I,j)+1) = RHOKM_TOP(I,j,NTML(I,j)+1)   &
     &                                   + RHOKM(I,j,NTML(I,j)+1)
              RHOKM(I,j,NTML(I,j)+1) = 0.0
            END IF
          END IF   ! test on WB_RATIO LE DEC_THRES
        END IF   ! testing for well-mixed layer (TEST_WELL_MIXED)

      END DO ! I
      END DO ! J
! ----------------------------------------------------------------------
! 2.2 Start iteration to find top of surface-driven mixing, ZSML_TOP,
!     within predetermined maximum and minimum height limits.
!     The solution is the height that gives WB_RATIO = DEC_THRES.
!     Procedure used makes 3 sweeps (up, down and up again), using
!     progressively smaller increments (Z_INC), each time stopping when
!     the buoyancy flux threshold or the height limits are reached.
!--------------------------------------------------------------------
!     If boundary layer is stable then ignore surface driven mixing.
!--------------------------------------------------------------------
!cdir collapse
      DO J=1,rows
      DO I=1,row_length

        IF ( KSURF_ITERATE(I,j) ) THEN
!         !-----------------------------------------------------------
!         ! Mixing must extend just above surface layer
!         ! (not clear precisely how to define this here: for now use
!         !  K_SURF calculated from ZH)
!         !-----------------------------------------------------------
!         ! limit K-surf to below cloud-base
          Z_BOT_LIM(i,j)=Z_UV(I,j,KSURF(I,j)+1)                         &
     &           + 0.1 * (Z_UV(I,j,KSURF(I,j)+2)-Z_UV(I,j,KSURF(I,j)+1))
!         ! limit K-surf to below cloud-top radiatively cooled layer
          Z_TOP_LIM(i,j)=MAX( Z_BOT_LIM(i,j), ZHSC(I,j)- DZRAD(I,j) )

          Z_CBASE(i,j) = ZHSC(I,j) - ZC_DSC(I,j)
!         !-----------------------------------------------------
!         ! Initial increment to ZSML_TOP found by dividing
!         ! up depth of layer within which it is allowed:
!         ! Start with ZSML_TOP at lower limit and work upwards
!         !-----------------------------------------------------
          Z_INC(i,j)=(Z_TOP_LIM(i,j)-Z_BOT_LIM(i,j))                    &
     &                 / FLOAT(N_STEPS)
          ZSML_TOP(I,j) = Z_BOT_LIM(i,j)

          WB_RATIO(i,j) = DEC_THRES(i,j) - 1.0 ! to be < DEC_THRES

        END IF ! KSURF_ITERATE

      END DO
      END DO

      ij_len=row_length*rows
      do i=1,ij_len
        to_do(i)    = .false.
        ind_todo(i) = i
        UP(i)       = 1
      end do

      c_len=ij_len

      l=0
      do j=1, rows
      do i=1, row_length
        l=l+1
        if (ksurf_iterate(i,j)) then
          to_do(l)=.true.
        end if
      end do
      end do

      do n_sweep=1,3

!       ! Compress to_do and ind_todo (will have new length c_len)
! DEPENDS ON: excfnl_cci
        call excfnl_cci(row_length, rows, c_len, to_do, ind_todo)

!       ! Restart inner interation with the points of outer
        c_len_i = c_len
        todo_inner(1:c_len_i) = to_do(1:c_len_i)
        ind_todo_i(1:c_len_i) = ind_todo(1:c_len_i)

        do ns=1,n_steps

!         ! Calculate active elements and compress
! DEPENDS ON: excfnl_compin
          call excfnl_compin(row_length, rows, c_len_i,                 &
     &                       ind_todo_i, todo_inner,                    &
     &                       UP, WB_RATIO, dec_thres, 1)

!cdir nodep
          do ic=1, c_len_i
            j1=(ind_todo_i(ic)-1)/row_length+1
            i1=ind_todo_i(ic)-(j1-1)*row_length

            ZSML_TOP(i1,j1)=ZSML_TOP(i1,j1)+Z_INC(i1,j1)

!           ! assume wb goes to zero at ZSML_TOP
            WB_SURF_INT(i1,j1) =                                        &
     &           BFLUX_SURF(i1,j1) * Z_TQ(i1,j1,KSURF(i1,j1)) *         &
     &           ( 1.0 - Z_TQ(i1,j1,KSURF(i1,j1))/                      &
     &                                      (2.0*ZSML_TOP(i1,j1)) )
            WB_SURF_INT(i1,j1) = MAX(1.0E-14,WB_SURF_INT(i1,j1))
            WBEND(i1,j1,KSURF(i1,j1)) = WB_SURF_INT(i1,j1)
!           ! Note: WB_DZRAD_INT not included as K_SURF restricted
!           !       to below zi-dzrad
            WBP_INT(i1,j1) = WB_SURF_INT(i1,j1)  ! must be > 0
            WBN_INT(i1,j1) = 0.0

            Z_INV(i1,j1) = ZSML_TOP(i1,j1)

          end do ! ic c_len_i
!
!..Integrate buoyancy flux profile given this ZSML_TOP
!
          do k=ksurf_min+1,ntop_max
!cdir nodep
            do ic=1, c_len_i
              j1=(ind_todo_i(ic)-1)/row_length+1
              i1=ind_todo_i(ic)-(j1-1)*row_length

              IF ( K  >=  KSURF(i1,j1)+1 .AND.                          &
     &             K  <=  NTOP(i1,j1) ) THEN

                Z_PR = Z_UV(i1,j1,K)
                IF (Z_PR  <   Z_INV(i1,j1)) THEN
                  KH_SURF(i1,j1,K) = W_H_TOP(i1,j1)                     &
     &                                  *Z_PR*(1.0-Z_PR/Z_INV(i1,j1) )  &
     &                                       *(1.0-Z_PR/Z_INV(i1,j1) )
                ELSE
                  KH_SURF(i1,j1,K) = 0.0
                END IF
!               !-----------------------------------------------------
!               ! No F2 or FSc terms here because we're only
!               ! considering effects driven from the surface
!               !-----------------------------------------------------
                IF (Z_CBASE(i1,j1)  >   Z_TQ(i1,j1,K)) THEN
!                 ! cloud-base above this range so use dry WB
                  WB_SCLD= KH_SURF(i1,j1,K) * DB_GA_DRY(i1,j1,K)
                  WB_CLD = 0.0
                ELSEIF (Z_CBASE(i1,j1)  <   Z_TQ(i1,j1,K-1)) THEN
!                 ! cloud-base below this range so use cloudy WB
                  WB_CLD = KH_SURF(i1,j1,K) * DB_GA_CLD(i1,j1,K)
                  WB_SCLD=0.0
                ELSE
!                 ! cloud-base within this integration range
!                 ! so treat cloud and sub-cloud layer wb separately
                  CLD_FRAC = (Z_TQ(i1,j1,K)-Z_CBASE(i1,j1))             &
     &                      /(Z_TQ(i1,j1,K)-Z_TQ(i1,j1,K-1))
                  WB_CLD  = CLD_FRAC                                    &
     &                       * KH_SURF(i1,j1,K)*DB_GA_CLD(i1,j1,K)
                  WB_SCLD = (1.0-CLD_FRAC)                              &
     &                       * KH_SURF(i1,j1,K)*DB_GA_DRY(i1,j1,K)
                END IF

                WBEND(i1,j1,K) = WB_CLD+WB_SCLD

                IF (WB_CLD  >=  0.0) THEN
                  WBP_INT(i1,j1)= WBP_INT(i1,j1) + WB_CLD
                ELSE
                  WBN_INT(i1,j1)= WBN_INT(i1,j1) - WB_CLD
                END IF
                IF (WB_SCLD  >=  0.0) THEN
                  WBP_INT(i1,j1)= WBP_INT(i1,j1) + WB_SCLD
                ELSE
                  WBN_INT(i1,j1) = WBN_INT(i1,j1)- WB_SCLD
                END IF

              END IF ! K

            END DO ! ic c_len_i
          END DO ! K

          DO ic=1,c_len_i
            j1=(ind_todo_i(ic)-1)/row_length+1
            i1=ind_todo_i(ic)-(j1-1)*row_length
            WB_RATIO(i1,j1) = WBN_INT(i1,j1)/WBP_INT(i1,j1)
          END DO ! ic c_len_i


        END DO  ! loop stepping up through ML (N_steps)

!cdir nodep
        DO ic=1, c_len
          l=ind_todo(ic)
          j1=(l-1)/row_length+1
          i1=l-(j1-1)*row_length

!..sub-divide current Z_INC into one more part than there will be steps
!..as there is no need to calculate WB for ZSML_TOP at a current Z_INC
          Z_INC(i1,j1)= Z_INC(i1,j1)/FLOAT(N_STEPS+1)

          IF ((UP(l) == 1.AND.WB_RATIO(i1,j1) >= DEC_THRES(i1,j1)).OR.  &
!                ! hit thres while working up
     &        (UP(l) == 0.AND.WB_RATIO(i1,j1) <= DEC_THRES(i1,j1))) THEN
!                ! hit thres while working down
                  UP(l) = 1-UP(l)   ! change direction of sweep
                  Z_INC(i1,j1)= - Z_INC(i1,j1)
          ELSEIF (ZSML_TOP(i1,j1) >= Z_TOP_LIM(i1,j1)-1.0) THEN
!                ! hit upper height limit (give-or-take 1m) without
!                ! reaching threshold
                  to_do(ic)=.false.
                  ZSML_TOP(i1,j1) = Z_TOP_LIM(i1,j1)
          ELSEIF (ZSML_TOP(i1,j1)  <=  Z_BOT_LIM(i1,j1)+ 1.0) THEN
!                ! hit lower height limit (give-or-take 1m) without
!                ! reaching threshold
                  to_do(ic)=.false.
                  ZSML_TOP(i1,j1) = Z_BOT_LIM(i1,j1)
          END IF

!..Note that if the threshold has not been passed then the next sweep
!..continues in the same direction (but with reduced increment).

        END DO ! c_len

      END DO ! n_sweep

      DO J=1, rows
      DO I=1, row_length
        ntml_new(I,j) = 2
        status_ntml(I,j)=.true.
      END DO
      END DO

      DO K=2,BL_LEVELS-2
!cdir collapse
        DO J=1, rows
        DO I=1, row_length
          IF ( KSURF_ITERATE(I,j) .and. status_ntml(i,j) ) THEN
!           ! -------------
!           ! find new NTML
!           ! -------------
            IF  (Z_UV(I,j,K+1)  <   ZSML_TOP(I,j)) THEN
              NTML_NEW(i,j) = K+1
            ELSE
              status_ntml(i,j)=.false.
            END IF
!           ! --------------------------------------------------------
!           ! Rounding error previously found to give
!           !      ZSML_TOP > Z_TOP_LIM = ZHSC
!           ! Test on ZSML_TOP hitting thresholds consequently changed
!           ! but also include the following failsafe tests here.
!           ! --------------------------------------------------------
            NTML(I,j) = MIN( NTDSC(I,j), NTML_NEW(i,j)-1 )
            ZH(I,j)   = MIN(  ZHSC(I,j), ZSML_TOP(I,j) )

          END IF  ! KSURF_ITERATE true

        END DO
        END DO
      END DO
! ----------------------------------------------------------------------
! 2.3 Now repeat the above procedure to find the base of the
!     top-driven K profile, ZDSC_BASE.
! ----------------------------------------------------------------------
!cdir collapse
      DO J=1, rows
      DO I=1, row_length

        IF ( KTOP_ITERATE(I,j) ) THEN

!         ! Lower limit on base of DSC layer
          Z_BOT_LIM(i,j) = 0.1 * ZH(I,j)
!         ! Upper limit on base of DSC layer
          Z_TOP_LIM(i,j) = ZHSC(I,j) - DZRAD(I,j)
!         ! If cumulus limit base of top-driven mixing to above ZH
          IF ( CUMULUS(I,j) ) THEN
            Z_BOT_LIM(i,j) = MIN( ZH(I,j), Z_TOP_LIM(i,j) )
          END IF

          Z_CBASE(i,j) = ZHSC(I,j) - ZC_DSC(I,j)

!..Divide up depth of layer within which ZDSC_BASE is allowed
          Z_INC(i,j)=(Z_TOP_LIM(i,j)-Z_BOT_LIM(i,j))                    &
     &                   /FLOAT(N_STEPS)
          ZDSC_BASE(I,j) = Z_BOT_LIM(i,j)
!                          ! will start at Z_BOT_LIM+Z_INC

          WB_RATIO(i,j) = DEC_THRES(i,j) + 1.0 ! to be > DEC_THRES

        END IF ! KTOP_ITERATE

      END DO
      END DO
!cdir collapse
      DO j=1,rows
      DO I=1,row_length
        IF ( KTOP_ITERATE(I,j) .AND. WB_DZRAD_INT(I,j)  <   0.0 ) THEN
!         !-------------------------------------------------------------
!         ! Estimation of wb integral over radiatively cooled cloud-top
!         ! region not yet performed (ie. DSC over stable surface) so
!         ! do it now.
!         !-------------------------------------------------------------
          Z_INV(i,j)  = ZHSC(I,j)
          Z_CBASE(i,j)= Z_INV(i,j) - ZC_DSC(I,j)
          CF_ML(I,j)  = CF_DSC(I,j)
          DF_CTOP(I,j)= DF_DSCT_OVER_CP(I,j)
          RHO_WE(I,j) = RHO_WE_DSC(I,j)
          DSL(I,j)    = DSL_DSC(I,j)
          DQW(I,j)    = DQW_DSC(I,j)
          TOTHF_ZI(I,j)= - RHO_WE(I,j)*DSL(I,j) + FT_NT_ZHSC(I,j)
          TOTQF_ZI(I,j)= - RHO_WE(I,j)*DQW(I,j) + FQ_NT_ZHSC(I,j)
          NTOP(I,j)   = NTDSC(I,j) - 1
          DZRAD(I,j)  = 100.0
          DO WHILE ( Z_TQ(I,j,NTOP(I,j))  >   Z_INV(I,j)-DZRAD(I,j)     &
     &         .AND. Z_TQ(I,j,NTOP(I,j)-1)  >   Z_INV(I,j)-Z_CBASE(I,j) &
     &         .AND. NTOP(I,j)  >   1 )
            NTOP(i,j) = NTOP(i,j) - 1
          END DO
          DZRAD(I,j) = Z_INV(I,j) - Z_TQ(I,j,NTOP(I,j))

          WSL_DZRAD_INT = DZRAD(I,j) *                                  &
     &                    ( 0.66*DF_CTOP(I,j) - RHO_WE(I,j)*DSL(I,j) )
          WQW_DZRAD_INT = - DZRAD(I,j) * RHO_WE(I,j) * DQW(I,j)

          WB_DZRAD_INT(I,j) = WSL_DZRAD_INT * (                         &
     &                          (1.0-CF_ML(I,j))*BTM(I,j,NTOP(i,j)+1) + &
     &                            CF_ML(I,j)*BTM_CLD(I,j,NTOP(i,j)+1) ) &
     &                      + WQW_DZRAD_INT * (                         &
     &                          (1.0-CF_ML(I,j))*BQM(I,j,NTOP(i,j)+1) + &
     &                            CF_ML(I,j)*BQM_CLD(I,j,NTOP(i,j)+1) )
          WB_DZRAD_INT(I,j) = G * WB_DZRAD_INT(I,j)
          WB_DZRAD_INT(I,j) = MAX( 0.0, WB_DZRAD_INT(I,j) )

          DO K = NTOP(I,j)+1, NTDSC(I,j)+1
            WBEND(I,j,K) = WB_DZRAD_INT(I,j)/DZRAD(I,j)
          END DO

        END IF

        WB_DZRAD_INT(I,j) = MAX( 1.0E-14, WB_DZRAD_INT(I,j) )

      END DO ! I
      END DO ! J

      ntop_max  = maxval(NTOP)
#if !defined(SCMA)
!     For bit-reproducibility make consistent
!     for all PE decompositions.
      CALL gc_imax(1, n_proc, istat_gc, ntop_max)
#endif

      ij_len=row_length*rows
      do i=1,ij_len
        to_do(i)    = .false.
        ind_todo(i) = i
        UP(i)     = 1
      end do

      c_len=ij_len

      l=0
      do j=1, rows
      do i=1, row_length
        l=l+1
        if (ktop_iterate(i,j)) then
          to_do(l)=.true.
        end if
      end do
      end do

      do n_sweep=1,3

!       ! Compress to_do and ind_todo (will have new length c_len)
! DEPENDS ON: excfnl_cci
        call excfnl_cci(row_length, rows, c_len, to_do, ind_todo)

!       ! Restart inner interation with the points of outer
        c_len_i = c_len
        todo_inner(1:c_len_i) = to_do(1:c_len_i)
        ind_todo_i(1:c_len_i) = ind_todo(1:c_len_i)

        do ns=1,n_steps

!         ! Calculate active elements and compress
! DEPENDS ON: excfnl_compin
          call excfnl_compin(row_length, rows, c_len_i,                 &
     &                       ind_todo_i, todo_inner,                    &
     &                       UP, WB_RATIO, dec_thres, 2)

!cdir nodep
          do ic=1, c_len_i
            j1=(ind_todo_i(ic)-1)/row_length+1
            i1=ind_todo_i(ic)-(j1-1)*row_length

            ZDSC_BASE(i1,j1) = ZDSC_BASE(i1,j1)+Z_INC(i1,j1)
            SCBASE(i1,j1)    = .TRUE.
!           ! Flag for NBDSC found (only needed for LockWhelan2006)
            IF (FLUX_GRAD  ==  LockWhelan2006) SCBASE(i1,j1) = .FALSE.
            FT_NT_DSCB(i1,j1)= 0.0
            FQ_NT_DSCB(i1,j1)= 0.0

            WBN_INT(i1,j1) = 0.0
            WBP_INT(i1,j1) = WB_DZRAD_INT(i1,j1)
            WBEND(i1,j1,KSURF(i1,j1)) = 0.0
            IF ( KSURF_ITERATE(i1,j1) .AND.                             &
     &           ZDSC_BASE(i1,j1)  <   ZSML_TOP(i1,j1) ) THEN
!             ! only include surface flux if K_SURF is included
!             ! in the wb calculation and K profiles overlap
              WBP_INT(i1,j1) = WBP_INT(i1,j1) + WB_SURF_INT(i1,j1)
              WBN_INT(i1,j1) = 0.0
              WBEND(i1,j1,KSURF(i1,j1)) = WB_SURF_INT(i1,j1)
            END IF

            ZINV_PR(i1,j1) = ZHSC(i1,j1)-ZDSC_BASE(i1,j1)

          end do ! ic c_len_i
!
!..Integrate buoyancy flux profile given this ZDSC_BASE
!
          do k=ksurf_min+1,ntop_max
!cdir nodep
            do ic=1, c_len_i
              j1=(ind_todo_i(ic)-1)/row_length+1
              i1=ind_todo_i(ic)-(j1-1)*row_length

              IF ((K >= KSURF(i1,j1)+1).and.(K <= NTOP(i1,j1))) THEN

                KHTOP(i1,j1) = 0.0
                F2           = 0.0
                FSc          = 0.0
                IF (.NOT. SCBASE(i1,j1) ) THEN
                  FT_NT_DSCB(i1,j1) = FT_NT(i1,j1,K)
                  FQ_NT_DSCB(i1,j1) = FQ_NT(i1,j1,K)
                END IF
                Z_PR = Z_UV(i1,j1,K) - ZDSC_BASE(i1,j1)

                IF (Z_PR >   0.0 .AND.Z_PR <  ZINV_PR(i1,j1)) THEN

                  IF (.NOT. SCBASE(i1,j1) ) THEN
                    SCBASE(i1,j1) = .TRUE.
                    Z_RATIO = (ZDSC_BASE(i1,j1)-Z_UV(i1,j1,K-1))        &
     &                       /(Z_UV(i1,j1,K)-Z_UV(i1,j1,K-1))
                    FT_NT_DSCB(i1,j1) = FT_NT(i1,j1,K-1) +              &
     &                    (FT_NT(i1,j1,K)-FT_NT(i1,j1,K-1))*Z_RATIO
                    FQ_NT_DSCB(i1,j1) = FQ_NT(i1,j1,K-1) +              &
     &                    (FQ_NT(i1,j1,K)-FQ_NT(i1,j1,K-1))*Z_RATIO
                  END IF

                  Z_RATIO = Z_PR/ZINV_PR(i1,j1)

                  IF (FLUX_GRAD  ==  LockWhelan2006) THEN
                    KHTOP(i1,j1) = 3.6 * VKMAN * RHO_UV(i1,j1,K)        &
     &                           * V_TOP_DSC(i1,j1) * ZINV_PR(i1,j1)    &
     &                * (Z_RATIO**3.0) * ( (1.0-Z_RATIO)*(1.0-Z_RATIO) )
                    F2 = RHO_UV(i1,j1,K) * 0.5 * Z_RATIO                &
     &                                   * 2.0**( Z_RATIO**4.0 )
                    IF ( V_SUM_DSC(i1,j1)  >   0.0 ) THEN
                      FSc = 3.5 * RHO_UV(i1,j1,K)                       &
     &                         * (V_TOP_DSC(i1,j1)/V_SUM_DSC(i1,j1))    &
     &                         * (Z_RATIO**3.0) * (1.0-Z_RATIO)
                    END IF
                  ELSE
                    KHTOP(i1,j1) = 0.85 * VKMAN * RHO_UV(i1,j1,K)       &
     &                     * V_TOP_DSC(i1,j1) * (( 1.0 - Z_RATIO )**0.8)&
     &                     * Z_PR * Z_RATIO
                  END IF

                END IF

                KHSURF(i1,j1) = 0.0
                IF ( ZDSC_BASE(i1,j1)  <   ZSML_TOP(i1,j1) ) THEN
!                 ! only include K_surf if profiles overlap
!                 ! otherwise layers are independent
                  KHSURF(i1,j1) = KH_SURF(i1,j1,K)
                END IF

                IF (FLUX_GRAD  ==  LockWhelan2006) THEN
                  WSLNG = (F2+FSc)*(TOTHF_ZI(i1,j1)-FT_NT_DSCB(i1,j1))  &
     &                        - ( FT_NT(i1,j1,K)-FT_NT_DSCB(i1,j1) )
                  WQWNG = (F2+FSc)*(TOTQF_ZI(i1,j1)-FQ_NT_DSCB(i1,j1))  &
     &                        - ( FQ_NT(i1,j1,K)-FQ_NT_DSCB(i1,j1) )
                END IF

                IF ( Z_TQ(i1,j1,K)  <=  Z_CBASE(i1,j1) ) THEN
!                 ! Completely below cloud-base so use cloud-free form
                  WB_SCLD = KHSURF(i1,j1)* DB_GA_DRY(i1,j1,K) +         &
     &                      KHTOP(i1,j1) * DB_NOGA_DRY(i1,j1,K)
                  IF (FLUX_GRAD  ==  LockWhelan2006) THEN
                    WB_SCLD = WB_SCLD + ( G/RDZ(i1,j1,K) ) *            &
     &                 ( BTM(i1,j1,K-1)*WSLNG + BQM(i1,j1,K-1)*WQWNG )
                  END IF
                  WB_CLD  = 0.0
                ELSEIF (Z_TQ(i1,j1,K-1)  >=  Z_CBASE(i1,j1)) THEN
!                 ! Completely above cloud-base so use cloudy formula
                  WB_CLD = ( KHSURF(i1,j1) * DB_GA_CLD(i1,j1,K) +       &
     &                       KHTOP(i1,j1)  * DB_NOGA_CLD(i1,j1,K) )
                  IF (FLUX_GRAD  ==  LockWhelan2006) THEN
                    WB_CLD = WB_CLD + ( G/RDZ(i1,j1,K) ) * (            &
     &                   ( BTM(i1,j1,K-1)*(1.0-CF_ML(i1,j1)) +          &
     &                     BTM_CLD(i1,j1,K-1)*CF_ML(i1,j1) )*WSLNG +    &
     &                   ( BQM(i1,j1,K-1)*(1.0-CF_ML(i1,j1)) +          &
     &                     BQM_CLD(i1,j1,K-1)*CF_ML(i1,j1) )*WQWNG )
                  END IF
                  WB_SCLD = 0.0
                ELSE
!                 ! cloud-base within this integration range
!                 ! so treat cloud and sub-cloud layer wb separately
                  WB_SCLD = KHSURF(i1,j1) * DB_GA_DRY(i1,j1,K) +        &
     &                      KHTOP(i1,j1) * DB_NOGA_DRY(i1,j1,K)
                  WB_CLD = ( KHSURF(i1,j1) * DB_GA_CLD(i1,j1,K) +       &
     &                       KHTOP(i1,j1)  * DB_NOGA_CLD(i1,j1,K) )
                  IF (FLUX_GRAD  ==  LockWhelan2006) THEN
                    WB_SCLD = WB_SCLD + ( G/RDZ(i1,j1,K) ) *            &
     &                 ( BTM(i1,j1,K-1)*WSLNG + BQM(i1,j1,K-1)*WQWNG )
                    WB_CLD = WB_CLD + ( G/RDZ(i1,j1,K) ) * (            &
     &                   ( BTM(i1,j1,K-1)*(1.0-CF_ML(i1,j1)) +          &
     &                     BTM_CLD(i1,j1,K-1)*CF_ML(i1,j1) )*WSLNG +    &
     &                   ( BQM(i1,j1,K-1)*(1.0-CF_ML(i1,j1)) +          &
     &                     BQM_CLD(i1,j1,K-1)*CF_ML(i1,j1) )*WQWNG )
                  END IF
                  CLD_FRAC = (Z_TQ(i1,j1,K)-Z_CBASE(i1,j1))             &
     &                      /(Z_TQ(i1,j1,K)-Z_TQ(i1,j1,K-1))
                  WB_CLD  = CLD_FRAC * WB_CLD
                  WB_SCLD = (1.0-CLD_FRAC) * WB_SCLD
                END IF

                WBEND(i1,j1,K) = WB_CLD + WB_SCLD

                IF (WB_CLD  >=  0.0) THEN
                  WBP_INT(i1,j1) = WBP_INT(i1,j1)+WB_CLD
                ELSE
                  WBN_INT(i1,j1) = WBN_INT(i1,j1)-WB_CLD
                END IF
                IF (WB_SCLD  >=  0.0) THEN
                  WBP_INT(i1,j1) = WBP_INT(i1,j1)+WB_SCLD
                ELSE
                  WBN_INT(i1,j1) = WBN_INT(i1,j1)-WB_SCLD
                END IF

              END IF ! K
            END DO ! ic c_len_i
          END DO ! K

          DO ic=1,c_len_i
            j1=(ind_todo_i(ic)-1)/row_length+1
            i1=ind_todo_i(ic)-(j1-1)*row_length
            WB_RATIO(i1,j1)=WBN_INT(i1,j1)/WBP_INT(i1,j1)
          END DO ! ic c_len_i

        END DO  ! loop stepping up through ML

!cdir nodep
        DO ic=1, c_len
          l=ind_todo(ic)
          j1=(l-1)/row_length+1
          i1=l-(j1-1)*row_length

!..sub-divide current Z_INC into one more part than there will be steps
!..as there is no need to recalculate WB at a current Z_INC
          Z_INC(i1,j1)= Z_INC(i1,j1)/FLOAT(N_STEPS+1)

          IF (                                                          &
     &       (UP(l) == 1 .AND.WB_RATIO(i1,j1) <= DEC_THRES(i1,j1)).OR.  &
!               ! hit thres while working up
     &       (UP(l) == 0 .AND.WB_RATIO(i1,j1) >= DEC_THRES(i1,j1))) THEN
!               ! hit thres while working down
            UP(l) = 1-UP(l)   ! change direction of sweep
            Z_INC(i1,j1)=- Z_INC(i1,j1)
          ELSEIF ( ZDSC_BASE(i1,j1) >= Z_TOP_LIM(i1,j1)-1.0 .OR.        &
     &             ZDSC_BASE(i1,j1) <=  Z_BOT_LIM(i1,j1)+1.0 ) THEN
!           ! hit height limits (give-or-take 1m) without
!           ! reaching threshold
            to_do(ic)=.false.
          END IF

!..Note that if the threshold has not been passed then the next sweep
!..continues in the same direction (but with reduced increment).

        END DO ! c_len

      END DO  ! loop over sweeps
!
! Convert integrated WB to profiles of WB itself for diagnostics
!
      DO K=1,ntop_max
        DO j=1,rows
        DO I=1,row_length
!         ! convert to m2/s-3
          IF ( K <= KSURF(I,j) ) THEN
            GAMMA_WBS = ( (WBMIX(I,j,KSURF(I,j))/Z_TQ(I,j,KSURF(I,j)))  &
     &                    - BFLUX_SURF(I,j)  )*2.0/Z_TQ(I,j,KSURF(I,j))
            WBMIX(I,j,K) = BFLUX_SURF(I,j) + GAMMA_WBS*Z_UV(I,j,K)

            GAMMA_WBS = ( (WBEND(I,j,KSURF(I,j))/Z_TQ(I,j,KSURF(I,j)))  &
     &                    - BFLUX_SURF(I,j)  )*2.0/Z_TQ(I,j,KSURF(I,j))
            WBEND(I,j,K) =  BFLUX_SURF(I,j) + GAMMA_WBS*Z_UV(I,j,K)
          ELSE IF ( K <= NTOP(I,j) ) THEN
            WBMIX(I,j,K)=WBMIX(I,j,K)*RDZ(I,j,K)
            WBEND(I,j,K)=WBEND(I,j,K)*RDZ(I,j,K)
          END IF
        END DO
        END DO
      END DO

      ntml_max  = maxval(NTML)
#if !defined(SCMA)
!     For bit-reproducibility make consistent
!     for all PE decompositions.
      CALL gc_imax(1, n_proc, istat_gc, ntml_max)
#endif
      ntdsc_max = maxval(NTDSC)
#if !defined(SCMA)
!     For bit-reproducibility make consistent
!     for all PE decompositions.
      CALL gc_imax(1, n_proc, istat_gc, ntdsc_max)
#endif


      DO K = MAX(1,ntop_max+1), MAX(ntml_max,ntdsc_max)+1
        DO j=1,rows
        DO I=1,row_length
!       ! convert to m2/s-3
          IF ( K >=  NTOP(I,j)+1 .AND. K <= NTDSC(I,j)+1 ) THEN
            WBEND(I,j,K) = WB_DZRAD_INT(I,j)/DZRAD(I,j)
          END IF
          IF ( K >=  NTOP(I,j)+1 .AND. K <= NTML(I,j)+1 ) THEN
            WBMIX(I,j,K) = WB_DZRAD_INT(I,j)/DZRAD(I,j)
          END IF
        END DO
        END DO
      END DO

! ----------------------------------------------------------------------
! 2.4 Set depth of cloud-top driven mixing in SML when there is a DSC
!     layer above (eg. fog under Sc) to be the SML layer depth
! ----------------------------------------------------------------------
!cdir collapse
      DO j=1, rows
      DO i=1, row_length
        IF (CUMULUS(I,j) .OR. COUPLED(I,j) ) THEN
!         ! ignore SML `cloud-top' driven mixing
          ZSML_BASE(I,j) = ZH(I,j)
          V_TOP(I,j)     = 0.0
        ELSE
          ZSML_BASE(I,j) = 0.1*ZH(I,j)
        END IF
      END DO  ! loop over j
      END DO  ! loop over I
!-----------------------------------------------------------------------
! 3.  Calculate factors required to ensure that the non-local turbulent
!     mixing coefficient profiles are continuous as the entrainment
!     level is approached.
!-----------------------------------------------------------------------
!cdir collapse
      do j=1,rows
      do i=1,row_length
        K=NTML(I,j)+1
        KH_TOP_FACTOR(I,j) = MAX( 0.7 , 1.0 - SQRT(                     &
     &           RHOKH_SURF_ENT(I,j) /                                  &
     &                 ( RHO_UV(I,j,K)*W_H_TOP(I,j)*VKMAN*ZH(I,j) ) ) )
        KM_TOP_FACTOR(I,j) = MAX( 0.7 , 1.0 - SQRT( RHOKM(I,j,K) /      &
     &             ( RHO_TQ(I,j,K-1)*W_M_TOP(I,j)*VKMAN*ZH(I,j) ) ) )
        SCDEPTH(I,j) = ZH(I,j) - ZSML_BASE(I,j)
        FACTOR = 0.85 * RHO_UV(I,j,K) * V_TOP(I,j) *VKMAN *SCDEPTH(I,j)
        IF ( FACTOR  >   0.0) THEN
          KH_SCT_FACTOR(I,j) = 1.0 -                                    &
     &                         ( RHOKH_TOP_ENT(I,j) / FACTOR )**1.25
!                                                         ! 1.25=1/0.8
        ELSE
          KH_SCT_FACTOR(I,j) = 1.0
        END IF
        FACTOR = 0.85 * RHO_TQ(I,j,K-1) * V_TOP(I,j) *                  &
     &                  VKMAN * SCDEPTH(I,j) * 0.75
        IF ( FACTOR  >   0.0) THEN
          KM_SCT_FACTOR(I,j) = 1.0 -                                    &
     &                         ( RHOKM_TOP(I,j,K) / FACTOR )**1.25
!                                                         ! 1.25=1/0.8
        ELSE
          KM_SCT_FACTOR(I,j) = 1.0
        END IF
!
        IF (NTDSC(I,j)  >   0) THEN
!       !-------------------------------------------------------------
!       ! Set up factors to ensure K profile continuity at ZHSC;
!       ! no need to limit size of factor as precise shape of top-down
!       ! mixing profile not important.
!       ! Only calculate _DSCT_FACTORs when a decoupled stratocumulus
!       ! layer exists, i.e. NTDSC > 0.
!       !-------------------------------------------------------------
          K=NTDSC(I,j)+1
          DSCDEPTH(I,j) = ZHSC(I,j) - ZDSC_BASE(I,j)
          FACTOR = 0.85*RHO_UV(I,j,K)*V_TOP_DSC(I,j)*VKMAN*DSCDEPTH(I,j)
          IF ( FACTOR  >   0.0) THEN
            KH_DSCT_FACTOR(I,j) = 1.0 -                                 &
     &                          ( RHOKH_DSCT_ENT(I,j) / FACTOR )**1.25
!                                                         ! 1.25=1/0.8
          ELSE
            KH_DSCT_FACTOR(I,j) = 1.0
          END IF

          FACTOR = 0.75 * 0.85 * RHO_TQ(I,j,K-1) * V_TOP_DSC(I,j) *     &
     &                         VKMAN * DSCDEPTH(I,j)
          IF ( FACTOR  >   0.0) THEN
            KM_DSCT_FACTOR(I,j) = 1.0 -                                 &
     &                         ( RHOKM_TOP(I,j,K) / FACTOR )**1.25
!                                                         ! 1.25=1/0.8
          ELSE
            KM_DSCT_FACTOR(I,j) = 1.0
          END IF
        END IF
      END DO
      END DO
!
!-----------------------------------------------------------------------
!! 4.  Calculate height dependent turbulent
!!     transport coefficients within the mixing layer.
!-----------------------------------------------------------------------
!
! Reset identifiers of base of decoupled layer mixing
!
!cdir collapse
      do j=1,rows
      do i=1,row_length
        SCBASE(I,j) = .FALSE.
        NBDSC(I,j)  = 0
      END DO
      END DO
!-------------------------------------------------------------
! Calculate RHOK(H/M)_TOP, top-down turbulent mixing profiles
! for the surface mixed layer.
! This is a variation on an up-side-down version of the cubic
! surface-forced profiles below.  Implement between at least
! the top of the `surface layer' (at Z=0.1*ZH) and ZH.
! Note this may well include NTML+1: entrainment fluxes will
! be dealt with in KMKHZ.
!-------------------------------------------------------------
      DO K=2,BL_LEVELS
!cdir collapse
         do j=1,rows
         do i=1,row_length
!
!           Calculate the height of u,v-level above the surface
! *APL: z0m removed from z in K(z)
            ZK_UV = Z_UV(I,j,K)
!
!           Calculate the height of T,q-level above the surface
!
            ZK_TQ = Z_TQ(I,j,K-1)
!
          IF ( ZK_UV  <   ZH(I,j) .AND.                                 &
     &         ZK_UV  >   ZSML_BASE(I,j) ) THEN
            Z_PR  = ZK_UV - ZSML_BASE(I,j)
            ZH_PR = ZH(I,j) - ZSML_BASE(I,j)
            Z_RATIO = Z_PR/ZH_PR
            IF (FLUX_GRAD  ==  LockWhelan2006) THEN

              RHOKH_TOP(I,j,K) = 3.6*VKMAN * RHO_UV(I,j,K) * V_TOP(I,j) &
     &                           * ZH_PR * (Z_RATIO**3.0)               &
     &                           * (( 1.0 - Z_RATIO )**2.0)

              IF ( .NOT.COUPLED(I,j) ) THEN
                RHOF2(I,j,K)  = RHO_UV(i,j,K) * 0.5 * Z_RATIO           &
     &                                      * 2.0**( Z_RATIO**4.0 )
                IF ( V_SUM(i,j)  >   0.0 ) THEN
                  RHOFSC(I,j,K) = 3.5 * RHO_UV(i,j,K)                   &
     &                            * (V_TOP(i,j)/V_SUM(i,j))             &
     &                            * (Z_RATIO**3.0) * (1.0-Z_RATIO)
                END IF
              END IF

            ELSE  ! Not LockWhelan2006

              RHOKH_TOP(I,j,K) = RHO_UV(I,j,K) * V_TOP(I,j) * 0.85 *    &
     &          VKMAN * ( ( 1.0 - KH_SCT_FACTOR(I,j)*Z_RATIO )**0.8 )   &
     &                                         * Z_PR * Z_RATIO
            END IF

          END IF
!         !-------------------------------------------------------
!         !   For LockWhelan2006, KM_TOP could be changed to match
!         !   the shape of KH_TOP.  This has not been done on the
!         !   grounds that the change in shape arises with the
!         !   inclusion of the other non-gradient terms.
!         !-------------------------------------------------------
          IF ( ZK_TQ  <   ZH(I,j) .AND.                                 &
     &         ZK_TQ  >   ZSML_BASE(I,j) ) THEN
              Z_PR = ZK_TQ - ZSML_BASE(I,j)
              ZH_PR = ZH(I,j) - ZSML_BASE(I,j)
              RHOKM_TOP(I,j,K) = 0.75 * RHO_TQ(I,j,K-1) * V_TOP(I,j) *  &
     &              0.85 * VKMAN *                                      &
     &              ( ( 1.0 - KM_SCT_FACTOR(I,j)*Z_PR/ZH_PR )**0.8 )    &
     &                                         * Z_PR * Z_PR / ZH_PR
!                                                     ! PRANDTL=0.75
          END IF
!         !-------------------------------------------------------------
!         ! Add contribution to top-down mixing coefficient
!         ! profiles for decoupled stratocumulus layers when
!         ! one exists
!         !-------------------------------------------------------------
          IF ( ZK_UV  <   ZHSC(I,j) .AND.                               &
     &            ZK_UV  >   ZDSC_BASE(I,j) ) THEN
            IF (.NOT. SCBASE(I,j) ) THEN
              SCBASE(I,j) = .TRUE.
!             ! identifies lowest layer below which there is mixing
              NBDSC(I,j) = K
            END IF
!           !-----------------------------------------------------------
!           ! Calculate RHOK(H/M)_TOP, top-down turbulent mixing
!           ! profiles and add to any generated in the surface mixing
!           ! layer.
!           ! This is a variation on an up-side-down version of the
!           ! cubic surface-forced profiles above.  Implement between
!           ! at least the top of the `surface layer' (at Z=0.1*ZH) and
!           ! ZHSC.
!           !-----------------------------------------------------------
            Z_PR = ZK_UV - ZDSC_BASE(I,j)
            ZH_PR = ZHSC(I,j) - ZDSC_BASE(I,j)
            Z_RATIO = Z_PR/ZH_PR


            IF (FLUX_GRAD  ==  LockWhelan2006) THEN

              RHOKH_TOP(I,j,K) = 3.6*VKMAN * RHO_UV(I,j,K)              &
     &                         * V_TOP_DSC(I,j) * ZH_PR * (Z_RATIO**3.0)&
     &                         * (( 1.0 - Z_RATIO )**2.0)

              RHOF2(I,j,K)  = RHO_UV(i,j,K) * 0.5 * Z_RATIO             &
     &                                      * 2.0**( Z_RATIO**4.0 )
              IF ( V_SUM_DSC(i,j)  >   0.0 ) THEN
                RHOFSC(I,j,K) = 3.5 * RHO_UV(i,j,K)                     &
     &                            * (V_TOP_DSC(i,j)/V_SUM_DSC(i,j))     &
     &                            * (Z_RATIO**3.0) * (1.0-Z_RATIO)
              END IF

            ELSE  ! Not LockWhelan2006

              RHOKH_TOP(I,j,K) = RHOKH_TOP(I,j,K) +                     &
     &           RHO_UV(I,j,K)*V_TOP_DSC(I,j)*0.85*VKMAN*               &
     &              ( ( 1.0 - KH_DSCT_FACTOR(I,j)*Z_RATIO )**0.8 )      &
     &                                         * Z_PR * Z_RATIO
            END IF
          END IF
!         !-------------------------------------------------------------
!         ! Now momentum
!         !-------------------------------------------------------------
          IF ( ZK_TQ  <   ZHSC(I,j) .AND.                               &
     &         ZK_TQ  >   ZDSC_BASE(I,j) ) THEN
            Z_PR = ZK_TQ - ZDSC_BASE(I,j)
            ZH_PR = ZHSC(I,j) - ZDSC_BASE(I,j)
              RHOKM_TOP(I,j,K) = RHOKM_TOP(I,j,K) +                     &
     &           0.75*RHO_TQ(I,j,K-1)*V_TOP_DSC(I,j)*0.85*VKMAN*        &
     &              ( ( 1.0 - KM_DSCT_FACTOR(I,j)*Z_PR/ZH_PR )**0.8 )   &
     &                                      * Z_PR * Z_PR / ZH_PR
          END IF

        end do
        end do
      END DO
!----------------------------------------------------
! Now K_SURF profiles
!----------------------------------------------------
      IF (FLUX_GRAD  ==  LockWhelan2006) THEN
!----------------------------------------------------
! Lock and Whelan formulation
!----------------------------------------------------
      C_WS = 0.42     ! ~ PR_NEUT^3 by design
      PR_NEUT = 0.75
      PR_CONV = 0.6

      DO K=2,BL_LEVELS
!cdir collapse
        do j=1,rows
        do i=1,row_length
!
          ZK_UV = Z_UV(I,j,K)
          ZK_TQ = Z_TQ(I,j,K-1)
!
          IF (FB_SURF(I,j)  >=  0.0) THEN
!
!           Calculate the free-convective scaling velocity at z(k)
!
            IF (COUPLED(I,j)) THEN  !  coupled
              WSTAR3 = ZHSC(I,j) * FB_SURF(I,j)
            ELSE
              WSTAR3 = ZH(I,j) * FB_SURF(I,j)
            END IF

            IF (ZK_UV  <=  0.1*ZH(I,j)) THEN
!             Surface layer calculation
              W_S_CUBED_UV = 10.*C_WS * ZK_UV * FB_SURF(I,j)
            ELSE
!             Outer layer calculation
              W_S_CUBED_UV = C_WS * WSTAR3
            END IF

            IF (ZK_TQ  <=  0.1*ZH(I,j)) THEN
!             Surface layer calculation
              W_S_CUBED_TQ = 10.*C_WS * ZK_TQ * FB_SURF(I,j)
            ELSE
!             Outer layer calculation
              W_S_CUBED_TQ = C_WS * WSTAR3
            END IF
!
!           Turbulent velocity scale for scalars
!
            W_M_NEUT = ( V_S(I,j)*V_S(I,j)*V_S(I,j) + W_S_CUBED_UV )    &
     &                             **(1.0/3.0)
            W_H_UV = W_M_NEUT/PR_NEUT
!
!           ! Also calc on TQ levels for W_M_TQ
            W_M_NEUT = ( V_S(I,j)*V_S(I,j)*V_S(I,j) + W_S_CUBED_TQ )    &
     &                             **(1.0/3.0)
            W_H_TQ = W_M_NEUT/PR_NEUT
!
!           Turbulent Prandtl number and velocity scale for scalars
!
            PRANDTL = PR_NEUT*                                          &
     &      ( V_S(I,j)*V_S(I,j)*V_S(I,j)*V_S(I,j) +                     &
     &       (1.0/(C_WS*25.0))*W_S_CUBED_TQ*W_M_NEUT ) /                &
     &      ( V_S(I,j)*V_S(I,j)*V_S(I,j)*V_S(I,j) +                     &
     &       (1.0/(C_WS*25.0))*(PR_NEUT/PR_CONV)*W_S_CUBED_TQ*W_M_NEUT )

            W_M_TQ = PRANDTL * W_H_TQ
!
            IF ( ZK_UV  <   ZH(I,j) ) THEN
!             !---------------------------------------------------------
!             ! Calculate RHOKH(w_h,z/z_h)
!             !---------------------------------------------------------
!
              RHOKH(I,j,K) = RHO_UV(I,j,K) * W_H_UV * VKMAN * ZK_UV *   &
     &                              ( 1.0 - ( ZK_UV / ZH(I,j) ) ) *     &
     &                              ( 1.0 - ( ZK_UV / ZH(I,j) ) )

            END IF
            IF ( ZK_TQ  <   ZH(I,j) ) THEN
!             !---------------------------------------------------------
!             ! Calculate RHOKM(w_m,z/z_h)
!             !---------------------------------------------------------
!
              RHOKM(I,j,K) = RHO_TQ(I,j,K-1) * W_M_TQ * VKMAN * ZK_TQ * &
     &                              ( 1.0 - ( ZK_TQ / ZH(I,j) ) ) *     &
     &                              ( 1.0 - ( ZK_TQ / ZH(I,j) ) )

            END IF
          END IF
        END DO
        END DO
      END DO

      ELSE
!----------------------------------------------------
! Lock et al and Holtstalg and Boville formulations
!----------------------------------------------------
!     ! Default to Lock et al
      C_WS = 0.25
      PR_NEUT = 0.75
      PR_CONV = 0.375
      IF (FLUX_GRAD  ==  HoltBov1993) THEN
        C_WS = 0.6
        PR_NEUT = 1.0
        PR_CONV = 0.6
      END IF

      DO K=2,BL_LEVELS
!cdir collapse
        do j=1,rows
        do i=1,row_length
!
!         Calculate the height of u,v-level above the surface
! *APL: z0m removed from z in K(z)
          ZK_UV = Z_UV(I,j,K)
!
!         Calculate the height of T,q-level above the surface
!
          ZK_TQ = Z_TQ(I,j,K-1)
!
          IF (FB_SURF(I,j)  >=  0.0) THEN
!
!           Calculate the free-convective scaling velocity at z(k)
!
            IF (ZK_UV  <=  0.1*ZH(I,j)) THEN
!
!             Surface layer calculation
!
              W_S_CUBED_UV = 10.*C_WS * ZK_UV * FB_SURF(I,j)
            ELSE
!
!             Outer layer calculation
!
              IF (COUPLED(I,j)) THEN  !  coupled and cloudy
               W_S_CUBED_UV = C_WS * ZHSC(I,j) * FB_SURF(I,j)
              ELSE
               W_S_CUBED_UV = C_WS * ZH(I,j) * FB_SURF(I,j)
              END IF
            END IF

            IF (ZK_TQ  <=  0.1*ZH(I,j)) THEN
!
!             Surface layer calculation
!
              W_S_CUBED_TQ = 10.*C_WS * ZK_TQ * FB_SURF(I,j)
            ELSE
!
!             Outer layer calculation
!
              IF (COUPLED(I,j)) THEN  !  coupled and cloudy
                W_S_CUBED_TQ = C_WS * ZHSC(I,j) * FB_SURF(I,j)
              ELSE
                W_S_CUBED_TQ = C_WS * ZH(I,j) * FB_SURF(I,j)
              END IF
            END IF
!
!           Turbulent velocity scale for momentum
!
            W_M_UV = (V_S(I,j)*V_S(I,j)*V_S(I,j) + W_S_CUBED_UV)        &
     &                             **(1.0/3.0)
!
            W_M_TQ = (V_S(I,j)*V_S(I,j)*V_S(I,j) + W_S_CUBED_TQ)        &
     &                             **(1.0/3.0)
!
!           Turbulent Prandtl number and velocity scale for scalars
!
            PRANDTL = PR_NEUT*( V_S(I,j)*V_S(I,j)*V_S(I,j)*V_S(I,j) +   &
     &         (1.0/(C_WS*25.0))*W_S_CUBED_UV*W_M_UV ) /                &
     &                        ( V_S(I,j)*V_S(I,j)*V_S(I,j)*V_S(I,j) +   &
     &         (1.0/(C_WS*25.0))*(PR_NEUT/PR_CONV)*W_S_CUBED_UV*W_M_UV )
            W_H_UV = W_M_UV / PRANDTL
!
            IF ( ZK_UV  <   ZH(I,j) ) THEN
!             !---------------------------------------------------------
!             ! Calculate RHOKH(w_h,z/z_h)
!             !---------------------------------------------------------
!
              RHOKH(I,j,K) = RHO_UV(I,j,K) * W_H_UV * VKMAN * ZK_UV *   &
     &            ( 1.0 - KH_TOP_FACTOR(I,j) * ( ZK_UV / ZH(I,j) ) ) *  &
     &            ( 1.0 - KH_TOP_FACTOR(I,j) * ( ZK_UV / ZH(I,j) ) )

            END IF
            IF ( ZK_TQ  <   ZH(I,j) ) THEN
!             !---------------------------------------------------------
!             ! Calculate RHOKM(w_m,z/z_h)
!             !---------------------------------------------------------
!
              RHOKM(I,j,K) = RHO_TQ(I,j,K-1) * W_M_TQ * VKMAN * ZK_TQ * &
     &            ( 1.0 - KM_TOP_FACTOR(I,j) * ( ZK_TQ / ZH(I,j) ) ) *  &
     &            ( 1.0 - KM_TOP_FACTOR(I,j) * ( ZK_TQ / ZH(I,j) ) )

            END IF
          END IF
        END DO
        END DO
      END DO

      END IF  ! Test on Flux_grad

      IF ( NG_STRESS  ==  BrownGrant97 .OR.                             &
     &     NG_STRESS  ==  BrownGrant97_limited ) THEN
        DO K=2,BL_LEVELS
!cdir collapse
          do j=1,rows
          do i=1,row_length
            ZK_TQ = Z_TQ(I,j,K-1)   ! stresses are calc on theta-levs
            IF ( FB_SURF(I,j)  >   0.0 .AND. ZK_TQ  <   ZH(I,j) ) THEN
!             !---------------------------------------------------------
!             ! Calculate non-gradient stress function
!             ! (Brown and Grant 1997)
!             ! Shape function chosen such that non-gradient stress
!             ! goes to zero at 0.1*ZH and ZH
!             !---------------------------------------------------------
              IF ( ZK_TQ  >   0.1*ZH(I,j) ) THEN
                Z_PR = ZK_TQ - 0.1*ZH(I,j)
                ZH_PR = 0.9*ZH(I,j)
!               !
!               ! Outer layer calculation
!               !
                IF (COUPLED(I,j)) THEN  !  coupled and cloudy
                  WSTAR3 = ZHSC(I,j) * FB_SURF(I,j)
                ELSE
                  WSTAR3 =   ZH(I,j) * FB_SURF(I,j)
                END IF

!               ! Use the Holtslag and Boville velocity scale for
!               ! non-gradient stress stability dependence, as in BG97
                W_M_HB_3 = V_S(I,j)*V_S(I,j)*V_S(I,j) + 0.6*WSTAR3
                F_NGSTRESS(I,j,K) = ( RHO_TQ(I,j,K-1)/RHOSTAR_GB(I,j) ) &
     &            * S_M * ( A_NGS * WSTAR3 / W_M_HB_3 )                 &
     &               * ( Z_PR / ZH_PR ) * ( 1.0 -  ( Z_PR / ZH_PR ) ) * &
     &                                    ( 1.0 -  ( Z_PR / ZH_PR ) )
              END IF
            END IF
          END DO
          END DO
        END DO
      END IF

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('EXCF_NL  ',4)
      END IF
      RETURN
      END SUBROUTINE EXCF_NL
#endif
