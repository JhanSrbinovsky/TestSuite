#if defined(A03_8B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!!!  SUBROUTINE EXCF_NL ----------------------------------------------
!!!
!!!  Purpose: To calculate non-local exchange coefficients for
!!!           boundary layer subroutine KMKH.
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!  5.4    Aug 2002  New deck is a restructuring of the diagnosis
!!!                   of the depth of mixing in EXCFNL6A
!!!                   to improve clarity and allow correct treatment
!!!                   of gradient adjustment and more accurate
!!!                   treatment of cloud base in wb integration.
!!!                                               A.P.Lock
!!!  5.5    Feb 2003  Calculate non-gradient stress function A.P.Lock
!!!  6.1    Aug 2003  Correct out-of-bounds errors. Improvement to
!!!                   diagnosis of BL types. R. J. Beare
!!!  6.2    Jan 2006  Optimisation by J-C Rioual.  A P Lock
!!!  6.2    Jan 2006  Include new switches DECFIX, STOPWE_SBL
!!!                         Adrian Lock
!!!  6.2    Jan 2004  Switch STOPWE_SBL to stop radiatively driven
!!!                   entrainment in stable BLs.     A.P.Lock
!!!  6.2    Dec 2005  Convert WBMIX and WBEND to units of
!!!                   m2/s-3.
!!!   6.2   8/02/06   Correction to non-gradient stress to
!!!                   follow Brown and Grant 97 and include
!!!                   density properly.   Bob Beare
!!!   6.4   Nov 2006  Remove out-of-bounds errors    Adrian Lock
!!!
!!!  Programming standard:
!!!
!!!  Documentation: UMDP No.24
!!!
!!!---------------------------------------------------------------------
!
!!   Arguments :-
      SUBROUTINE EXCF_NL (                                              &
     & row_length,rows,BL_LEVELS,                                       &
     & NTML,CF,                                                         &
     & RDZ,ZH,Z_UV,Z_TQ,RHO_UV,RHO_TQ,RHOSTAR_GB,                       &
     & Z0M,V_S,FB_SURF,DB_TOP,                                          &
     & BFLUX_SURF,BFLUX_SURF_SAT,ZETA_S,BT_TOP,BTT_TOP,                 &
     & DF_TOP_OVER_CP,ZETA_R,BTC_TOP,                                   &
     & DB_TOP_CLD,CHI_S_TOP,ZC,                                         &
     & RHOKM,RHOKH,RHOKM_TOP,RHOKH_TOP,                                 &
     & DECFIX,STOPWE_SBL,NG_STRESS,F_NGSTRESS,                          &
     & ZHSC,DSCDEPTH,NTDSC,DB_DSCT,SVL,                                 &
     & BT_DSCT,BTT_DSCT,                                                &
     & DF_DSCT_OVER_CP,ZETA_R_DSC,BTC_DSCT,                             &
     & DB_DSCT_CLD,CHI_S_DSCT,ZC_DSC,COUPLED,                           &
     & D_SIEMS,D_SIEMS_DSC,NBDSC,                                       &
     & DB_KSURF_DRY,DB_KTOP_DRY,DB_KSURF_CLD,DB_KTOP_CLD,               &
     & DSC,CUMULUS,ZDSC_BASE,                                           &
     & RHOKH_TOP_ENT,RHOKH_DSCT_ENT,RHOKH_SURF_ENT,                     &
     & LTIMER                                                           &
     &)
!
      IMPLICIT NONE

      LOGICAL LTIMER

      INTEGER                                                           &
     & row_length,rows                                                  &
     &,BL_LEVELS                                                        &
                   ! IN maximum number of boundary layer levels
     &,NG_STRESS                                                        &
                   ! IN switch for non-gradient stress
     &,DECFIX                                                           &
                              ! IN correction to decoupling diagnosis
     &,STOPWE_SBL             ! IN stop spurious entrainment in SBLs

      LOGICAL                                                           &
     & COUPLED(row_length,rows)                                         &
                                ! IN  Flag to indicate Sc layer weakly
!                               !     coupled to surface (ie weakly
                                !     decoupled)
     &,CUMULUS(row_length,rows)                                         &
                                ! INOUT Flag for cumulus
     &,DSC(row_length,rows)     ! INOUT Flag set if decoupled stratocu
!                               !       layer found.
      INTEGER                                                           &
     & NTML(row_length,rows)                                            &
                                ! IN  Number of turbulently mixed
!                               !     layers.
     &,NTDSC(row_length,rows)   ! IN  Top level of any decoupled
!                               !     turbulently mixed Sc layer.
!
      REAL                                                              &
     & RDZ(row_length,rows,BL_LEVELS)                                   &
                                ! IN Reciprocal of distance between
                                !    T,q-levels (m^-1). 1/RDZ(,K) is
                                !    the vertical distance from level
                                !    K-1 to level K, except that for
                                !    K=1 it is just the height of the
                                !    lowest atmospheric level.
     &,ZH(row_length,rows)                                              &
                                ! IN Boundary layer height (m).
     &,Z_UV(row_length,rows,BL_LEVELS)                                  &
                                ! IN For a vertically staggered grid
                                !    with a u,v-level first above the
                                !    surface, Z_UV(*,K) is the height
                                !    of the k-th u,v-level (half level
                                !    k-1/2) above the surface;
                                !    for an unstaggered grid the
                                !    heights of the half-levels
                                !    0.5 to BL_LEVELS-0.5 should be
                                !    input to elements 1 to BL_LEVELS.
                                !    (1st value not used in either
                                !     case.)
     &,Z_TQ(row_length,rows,BL_LEVELS)                                  &
                                ! IN For a vertically staggered grid
                                !    with a u,v-level first above the
                                !    surface, Z_TQ(*,K) is the height
                                !    of the k-th T,q-level (full level
                                !    k) above the surface;
                                !    code no longer works for an
                                !    unstaggered grid as Z_TQ is used
!                               !    to calculate K_SURF
     &,RHO_UV(row_length,rows,BL_LEVELS)                                &
                                ! IN For a vertically staggered grid
                                !    with a u,v-level first above the
                                !    surface, RHO_UV(*,K) is the
                                !    density at the k-th u,v-level
                                !    above the surface;
                                !    for an unstaggered grid the
                                !    densities at the layer interfaces
                                !    (half-levels) 0.5 to BL_LEVELS-0.5
                                !    should be input to elements 1 to
                                !    BL_LEVELS.
                                !    (1st value not used in either
                                !    case.)
     &,RHO_TQ(row_length,rows,BL_LEVELS)                                &
                                ! IN For a vertically staggered grid
                                !    with a u,v-level first above the
                                !    surface, RHO_TQ(*,K) is the
                                !    density of the k-th T,q-level
                                !    above the surface;
                                !    for an unstaggered grid the
                                !    densities at the layer interfaces
                                !    (half-levels) 1.5 to BL_LEVELS+0.5
                                !    should be input to elements 1 to
                                !    BL_LEVELS.
                                !    (value for BL_LEVELS not used
                                !    in either case.)
     &,RHOSTAR_GB(row_length,rows)                                      &
                                ! IN Surface density (kg/m3)

     &,Z0M(row_length,rows)                                             &
                                ! IN Roughness length for momentum (m).
     &,V_S(row_length,rows)                                             &
                                ! IN Surface friction velocity (m/s).
     &,FB_SURF(row_length,rows)                                         &
                                ! IN Buoyancy flux at the surface over
!                               !    density (m^2/s^3).
     &,BFLUX_SURF(row_length,rows)                                      &
                                ! IN Surface buoyancy flux (kg/m/s^3).
     &,BFLUX_SURF_SAT(row_length,rows)                                  &
                                ! IN Saturated-air surface buoyancy
                                !    flux.
     &,DB_TOP(row_length,rows)                                          &
                                ! IN Buoyancy jump across top of b.l
!                               !    (m/s^2)
     &,DF_TOP_OVER_CP(row_length,rows)                                  &
                                ! IN Radiative flux change at cloud top
                                !    divided by c_P (K.kg/m^2/s).
     &,BT_TOP(row_length,rows)                                          &
                                ! IN Buoyancy parameter at the top of
                                !    the b.l. (m/s^2/K).
     &,BTT_TOP(row_length,rows)                                         &
                                ! IN In-cloud buoyancy parameter at
                                !    the top of the b.l. (m/s^2/K).
     &,BTC_TOP(row_length,rows)                                         &
                                ! IN Cloud fraction weighted buoyancy
                                !    parameter at the top of the b.l.
     &,DB_TOP_CLD(row_length,rows)                                      &
                                  ! IN In-cloud buoyancy jump at the
                                !    top of the b.l. (m/s^2).
     &,CHI_S_TOP(row_length,rows)                                       &
                                ! IN Mixing fraction of just saturated
                                !    mixture at top of the b.l.
     &,ZETA_S(row_length,rows)                                          &
                                ! IN Non-cloudy fraction of mixing
                                !    layer for surface forced
                                !    entrainment term.
     &,ZETA_R(row_length,rows)                                          &
                                ! IN Non-cloudy fraction of mixing
                                !    layer for cloud top radiative
                                !    cooling entrainment term.
     &,ZC(row_length,rows)                                              &
                                ! IN Cloud depth (not cloud fraction
                                !    weighted) (m).
     &,CF(row_length,rows,BL_LEVELS) ! IN Cloud fraction
!
      REAL                                                              &
     & ZHSC(row_length,rows)                                            &
                                ! IN Cloud-layer height (m).
     &,DSCDEPTH(row_length,rows)                                        &
                                ! IN Decoupled cloud-layer depth (m).
     &,DB_DSCT(row_length,rows)                                         &
                                ! IN Buoyancy parameter at the top of
                                !    the DSC layer (m/s^2/K).
     &,SVL(row_length,rows,BL_LEVELS)                                   &
                                ! IN s_VL (K) *APL*PD
     &,DF_DSCT_OVER_CP(row_length,rows)                                 &
                                ! IN Radiative flux change at DSC top
                                !    divided by c_P (K.kg/m^2/s).
     &,BT_DSCT(row_length,rows)                                         &
                                 ! IN Buoyancy parameter at the top of
                                 !    the DSC  (m/s^2/K).
     &,BTT_DSCT(row_length,rows)                                        &
                                 ! IN In-cloud buoyancy parameter at
                                 !    the top of the DSC (m/s^2/K).
     &,BTC_DSCT(row_length,rows)                                        &
                                 ! IN Cloud fraction weighted buoyancy
                                 !    parameter at the top of the DSC
     &,DB_DSCT_CLD(row_length,rows)                                     &
                                !     IN In-cloud buoyancy jump at the
                                !    top of the DSC (m/s^2).
     &,CHI_S_DSCT(row_length,rows)                                      &
                                ! IN Mixing fraction of just saturated
                                !    mixture at top of the DSC
     &,ZETA_R_DSC(row_length,rows)                                      &
                                ! IN Non-cloudy fraction of DSC
                                !    for cloud top radiative
                                !    cooling entrainment term.
     &,ZC_DSC(row_length,rows)                                          &
                                 ! IN Cloud depth (not cloud fraction
                                 !    weighted) for DSC (m).
     &,D_SIEMS(row_length,rows)                                         &
                                 ! IN Siems (1990) et al. cloud-top
                                 !    entr.t instab. parm
     &,D_SIEMS_DSC(row_length,rows)                                     &
                                 !IN Siems (1990) et al. cloud-top
                                 !   entr.t instab. parm for DSC layer
     &,DB_KSURF_DRY(row_length,rows,2:BL_LEVELS)                        &
                                                 ! IN Dry buoyancy jump
!                                ! flux integral calculation (m/s2)
     &,DB_KTOP_DRY(row_length,rows,2:BL_LEVELS)                         &
                                                ! IN Sat. buoyancy jump
!                                ! flux integral calculation (m/s2)
     &,DB_KSURF_CLD(row_length,rows,2:BL_LEVELS)                        &
                                                 ! IN Dry buoyancy jump
!                                ! flux integral calculation (m/s2)
     &,DB_KTOP_CLD(row_length,rows,2:BL_LEVELS) ! IN Sat. buoyancy jump
!                                ! flux integral calculation (m/s2)
      INTEGER                                                           &
     & NBDSC(row_length,rows)    ! OUT Bottom level of any decou
!                                !     turbulently mixed Sc layer.
      REAL                                                              &
     & RHOKM(row_length,rows,2:BL_LEVELS)                               &
                                 ! OUT Layer k-1 - to - layer k
                                 !     turbulent mixing coefficient
                                 !     for momentum (kg/m/s).
     &,RHOKH(row_length,rows,2:BL_LEVELS)                               &
                                         ! OUT Layer k-1 - to - layer k
                                !     turbulent mixing coefficient
                                !     for heat and moisture (kg/m/s).
     &,RHOKM_TOP(row_length,rows,2:BL_LEVELS)                           &
                                ! OUT exchange coefficient for
                                !     momentum due to top-down mixing
     &,RHOKH_TOP(row_length,rows,2:BL_LEVELS)                           &
                                ! OUT exchange coefficient for
                                !     heat and moisture due to top-down
                                !     mixing
     &,F_NGSTRESS(row_length,rows,2:BL_LEVELS)                          &
                                ! OUT dimensionless function for
!                               !     non-gradient stresses
     &,ZDSC_BASE(row_length,rows)                                       &
                                 ! OUT Height of base of K_top in DSC
     &,RHOKH_SURF_ENT(row_length,rows)                                  &
                                      ! OUT SML surf-driven entr. KH
     &,RHOKH_TOP_ENT(row_length,rows)                                   &
                                      ! OUT SML top-driven entr. KH
     &,RHOKH_DSCT_ENT(row_length,rows)! OUT DSC top-driven entr. KH
!*
!*L---------------------------------------------------------------------
      EXTERNAL TIMER
!*
!*L---------------------------------------------------------------------
!    Local and other symbolic constants :-
#include "blopt8a.h"
#include "c_lheat.h"
#include "c_r_cp.h"
#include "c_epslon.h"
#include "c_vkman.h"
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
      LOGICAL                                                           &
     & SCBASE(row_length,rows)                                          &
                                   ! Flag to signal base of CML reached
     &,KSURF_ITERATE(row_length,rows)                                   &
                                      ! Flag to perform iteration to
!                                     ! find top of Ksurf
     &,KTOP_ITERATE(row_length,rows)  ! Flag to perform iteration to
!                                     ! find base of Ktop

      REAL                                                              &
     & W_M_TOP(row_length,rows)                                         &
                                ! Turbulent velocity scale for momentum
!                               ! evaluated at the top of the b.l.
     &,W_H_TOP(row_length,rows)                                         &
                                ! Turbulent velocity scale for scalars
!                               ! evaluated at the top of the b.l.
     &,PRANDTL_TOP(row_length,rows)                                     &
                                ! Turbulent Prandtl number
                                ! evaluated at the top of the b.l.
     &,KH_TOP_FACTOR(row_length,rows)                                   &
                                ! Factor to ensure K_H profile is
!                               ! continuous at z_h.
     &,KM_TOP_FACTOR(row_length,rows)                                   &
                                ! Factor to ensure K_M profile is
!                               ! continuous at z_h.
     &,V_TOP(row_length,rows)                                           &
                                ! velocity scale for top-down convection
     &,V_TOP_DSC(row_length,rows)                                       &
                                ! velocity scale for top-down convection
     &,KH_SCT_FACTOR(row_length,rows)                                   &
                                ! Factor to ensure K_H profile is
!                               ! continuous at z_h.
     &,KM_SCT_FACTOR(row_length,rows)                                   &
                                ! Factor to ensure K_M profile is
!                               ! continuous at z_h.
     &,KH_DSCT_FACTOR(row_length,rows)                                  &
                                ! Factor to ensure K_H profile is
!                               ! continuous at z_h.
     &,KM_DSCT_FACTOR(row_length,rows)                                  &
                                ! Factor to ensure K_M profile is
!                               ! continuous at z_h.
     &,ZSML_TOP(row_length,rows)                                        &
                                ! Height of top of surf-driven K in S
     &,ZSML_BASE(row_length,rows)                                       &
                                ! Height of base of top-driven K in S
     &,V_SURF(row_length,rows)                                          &
                                ! Velocity scale for surface-up conve
     &,KH_SURF(row_length,rows,2:BL_LEVELS)                             &
                                ! Shape factor for non-local
!                               ! turbulent mixing coefficient
     &,SCDEPTH(row_length,rows)                                         &
                                ! Depth of top-driven mixing in
     &,WBMIX(row_length,rows,BL_LEVELS)                                 &
                                        ! WB*DZ if were diag as mixed
     &,WBEND(row_length,rows,BL_LEVELS) ! WB*DZ after dec diag

      INTEGER                                                           &
     & KSURF(row_length,rows)   ! First Theta-level above surface layer
                                !  well-mixed SC layer

!  (b) Scalars.
!
      REAL                                                              &
     & PRANDTL                                                          &
                    ! Turbulent Prandtl number.
     &,ZK_UV                                                            &
                    ! Height above surface of u,v-level.
     &,ZK_TQ                                                            &
                    ! Height above surface of T,q-level.
     &,W_S_CUBED_UV                                                     &
                    ! Cube of free-convective velocity scale (u,v-level)
     &,W_S_CUBED_TQ                                                     &
                    ! Cube of free-convective velocity scale (T,q-level)
     &,W_M_UV                                                           &
                    ! Turbulent velocity scale for momentum (u,v-level).
     &,W_M_TQ                                                           &
                    ! Turbulent velocity scale for momentum (T,q-level).
     &,W_H_UV                                                           &
                    ! Turbulent velocity scale for scalars (u,v-level).
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
     &,ZCML_BASE                                                        &
                    ! Height of base of cloud mixed layer
     &,RHOKH_ENT                                                        &
                    ! entrainment eddy viscosity
     &,FRAC_TOP                                                         &
                    ! Fraction of turbulent mixing driven from the top
     &,FACTOR                                                           &
                    ! Temporary scalar
     &,ENT_FACTOR                                                       &
                    ! Factor to weight entrainment by CF
     &,V_SUM                                                            &
                    ! generalised turbulent velocity scale (m/s)
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
     &,GAMMA_WBS     ! Surface layer wb gradient
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
      REAL, DIMENSION(row_length,rows) ::                               &
     & Z_INV                                                            &
                     ! inversion height (top of K profiles)
     &,WB_SURF_INT                                                      &
                     ! Estimate of wb intqegrated over surface layer
     &,V_KTOP                                                           &
                     ! velocity scale for K_top profile
     &,Z_CBASE                                                          &
                     ! cloud base height
     &,WB_RATIO                                                         &
                     ! WBN_INT/WBP_INT
     &,WBP_INT                                                          &
                     ! Positive part of buoyancy flux integral
     &,WBN_INT                                                          &
                     ! Negative part of buoyancy flux integral
     &,DEC_THRES                                                        &
                     ! Local decoupling threshold
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
     &,Z_INC         ! Step size (m)

      INTEGER, DIMENSION(row_length,rows) ::                            &
     & ntop                                                             &
                     ! top level of surf-driven K profile
     &,kwb0          ! level at which wb assumed to go to zero

      integer, dimension(rows*row_length) ::                            &
     & up            ! indicator of upward/downward sweep

      LOGICAL, DIMENSION(row_length,rows) :: status_ntop

      INTEGER :: ksurf_min, ntop_max, ntdsc_max

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

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('EXCF_NL  ',3)
      END IF
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
         V_SURF(I,j) = 0.0
         V_TOP_DSC(I,j) = 0.0

         IF (FB_SURF(I,j)  >=  0.0) THEN
!
!         By definition the top of the b.l. is in the 'outer layer' so
!         the free-convective velocity scale cubed is
!
          IF (COUPLED(I,j)) THEN
           W_S_CUBED_UV = 0.25 * ZHSC(I,j) * FB_SURF(I,j)
          ELSE
           W_S_CUBED_UV = 0.25 * ZH(I,j) * FB_SURF(I,j)
          END IF
!
!         Turbulent velocity scale for momentum
!
          W_M_TOP(I,j) = (V_S(I,j)*V_S(I,j)*V_S(I,j) +                  &
     &                    W_S_CUBED_UV)**(1.0/3.0)
!
!         Turbulent Prandtl number and velocity scale for scalars
!
          PRANDTL_TOP(I,j) = 0.75 *                                     &
     &                    ( V_S(I,j)*V_S(I,j)*V_S(I,j)*V_S(I,j) +       &
     &                        (4.0/25.0)*W_S_CUBED_UV*W_M_TOP(I,j) ) /  &
     &                         ( V_S(I,j)*V_S(I,j)*V_S(I,j)*V_S(I,j) +  &
     &                        (8.0/25.0)*W_S_CUBED_UV*W_M_TOP(I,j) )
          W_H_TOP(I,j) = W_M_TOP(I,j) / PRANDTL_TOP(I,j)
        ELSE
          W_M_TOP(I,j) = V_S(I,j)
          PRANDTL_TOP(I,j) = 0.75
          W_H_TOP(I,j) = W_M_TOP(I,j) / PRANDTL_TOP(I,j)
        END IF
       END DO
      END DO
!
!-----------------------------------------------------------------------
!! 1.  Loop round levels; calculate the top-of-b.l. entrainment
!!     mixing coefficients.
!-----------------------------------------------------------------------
!
      DO K=2,BL_LEVELS
!cdir collapse
           do j=1,rows
            do i=1,row_length
!
!-----------------------------------------------------------------------
!! 1.2 Calculate top-of-b.l. entrainment mixing coefficients
!!     and store b.l. top quantities for later use.
!-----------------------------------------------------------------------
!      FIRST the top of the SML (if not coupled)
!-----------------------------------------------
!..Initialise RHOKs: entrainment now added later for KH, in KMKHZ
           RHOKH(I,j,K) = 0.0
           RHOKM(I,j,K) = 0.0
           RHOKH_TOP(I,j,K) = 0.0
           RHOKM_TOP(I,j,K) = 0.0
           F_NGSTRESS(I,j,K) = 0.0
           KH_SURF(I,j,K) = 0.0
           IF ( K  ==  NTML(I,j)+1 .AND. .NOT.COUPLED(I,j) .AND.        &
     &          ( STOPWE_SBL  ==  0 .OR. FB_SURF(I,j)  >=  0.0 ) ) THEN
!            ! Correction controlled by STOPWE_SBL:
!            !   if ne 0, then requires FB_SURF>0 too

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
            IF (FB_SURF(I,j)  >=  0.0) THEN
            SF_SHEAR_TERM =  A_ENT_SHR * V_S(I,j) * V_S(I,j) * V_S(I,j) &
     &                      * RHO_UV(I,j,K)  / ZH(I,j)
            ELSE
              SF_SHEAR_TERM = 0.0
            END IF
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
            EVAP_TERM = A_ENT_2 * RHO_UV(I,j,K)                         &
     &                  * CHI_S_TOP(I,j) * CHI_S_TOP(I,j)               &
     &                  * ZR * ZR * ZR * DB_TOP_CLD(I,j)                &
     &                  * SQRT( ZH(I,j) * DB_TOP(I,j) )
            ELSE
              EVAP_TERM = 0.0
            END IF
!            IF (CF(I,j,K-1) >= 0.9) THEN
              ENT_FACTOR = 1.0
!            ELSE
!              ENT_FACTOR = EXP(-((0.90-CF(I,j,K-1))**3.0)/0.075)
!            END IF
!           !-----------------------------------------------------------
!           ! Combine forcing terms to calculate the representative
!           ! velocity scales
!           !-----------------------------------------------------------
            V_SUM   = ( (SF_TERM + SF_SHEAR_TERM +                      &
     &                   ENT_FACTOR*(IR_TERM + EVAP_TERM))              &
     &                   * ZH(I,j) /(A_ENT_1*RHO_UV(I,j,K)) )**(1.0/3.0)
            V_TOP(I,j) = ( ENT_FACTOR*(IR_TERM+EVAP_TERM) * ZH(I,j)     &
     &                           / (A_ENT_1*RHO_UV(I,j,K)) )**(1.0/3.0)
            V_SURF(I,j) = ( (SF_TERM) * ZH(I,j)                         &
     &                           / (A_ENT_1*RHO_UV(I,j,K)) )**(1.0/3.0)
!           !-----------------------------------------------------------
!           ! Calculate the direct radiative term
!           !  can only calculate for DB_TOP > 0
!           !-----------------------------------------------------------
            IF ( DB_TOP(I,j)  >   0.0) THEN
             DZ_INV  = MIN( V_SUM * V_SUM / DB_TOP(I,j) ,100.0 )
            L_RAD   = 15.0 * MAX( 1.0 , 200./(ZC(I,j)+1.0E-14) )
             ALPHA_T = 1.0 - EXP(-0.5*DZ_INV/L_RAD)
            IF (D_SIEMS(I,j)  >   0.0) ALPHA_T =                        &
     &             MIN( 1.0, ALPHA_T + 10.0*D_SIEMS(I,j)*(1.0-ALPHA_T) )
            DR_TERM = BTC_TOP(I,j) * ALPHA_T * DF_TOP_OVER_CP(I,j)
!           !-----------------------------------------------------------
!           ! Combine terms to calculate the entrainment
!           ! mixing coefficients
!           !-----------------------------------------------------------
            ZIL_CORR = C_T * ( (SF_TERM + SF_SHEAR_TERM + ENT_FACTOR *  &
     &               (IR_TERM + EVAP_TERM)) /                           &
     &                   (RHO_UV(I,j,K) * SQRT(ZH(I,j))) )**(2.0/3.0)

            RHOKH_ENT = (SF_TERM + SF_SHEAR_TERM                        &
     &                   + ENT_FACTOR*(IR_TERM + EVAP_TERM + DR_TERM))  &
     &                      / ((DB_TOP(I,j) + ZIL_CORR) * RDZ(I,j,K) )

             FRAC_TOP = V_TOP(I,j) / ( V_TOP(I,j)+W_H_TOP(I,j)+1.0E-14 )

             RHOKH_SURF_ENT(I,j) = RHOKH_ENT * ( 1.0 - FRAC_TOP )
             RHOKH_TOP_ENT(I,j) = RHOKH_ENT * FRAC_TOP

! APL change for C-P grid:
!             RHOKM(I,j,K) = PRANDTL_TOP(I,j) * RHOKH_SURF_ENT(I,j)
!     &                   * RHO_TQ(I,j,K-1) / RHO_UV(I,j,K)
             RHOKM(I,j,K) = PRANDTL_TOP(I,j) * RHOKH_SURF_ENT(I,j)      &
     &                    * RDZ(I,j,K) * (Z_UV(I,j,K)-Z_UV(I,j,K-1))    &
     &                    * RHO_TQ(I,j,K-1) / RHO_UV(I,j,K)
!             RHOKM_TOP(I,j,K) = 0.75 * RHOKH_TOP_ENT(I,j)
!     &                   * RHO_TQ(I,j,K-1) / RHO_UV(I,j,K)
             RHOKM_TOP(I,j,K) = 0.75 * RHOKH_TOP_ENT(I,j)               &
     &                    * RDZ(I,j,K) * (Z_UV(I,j,K)-Z_UV(I,j,K-1))    &
     &                    * RHO_TQ(I,j,K-1) / RHO_UV(I,j,K)
            END IF    ! test on DB_TOP GT 0
          END IF
!----------------------------------------------------------------
!      THEN the top of the DSC (if coupled use ZHSC length-scale)
!----------------------------------------------------------------
          IF ( (K  ==  NTDSC(I,j)+1) ) THEN
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
!           !-----------------------------------------------------------
!           ! Calculate the indirect radiative term
!           !-----------------------------------------------------------
            ZETA_R_SQ = ZETA_R_DSC(I,j)*ZETA_R_DSC(I,j)
            IR_TERM = ( BT_DSCT(I,j)*ZETA_R_SQ +                        &
     &                  BTT_DSCT(I,j)*(1.0-ZETA_R_SQ) )                 &
     &                  * A_ENT_1 * DF_DSCT_OVER_CP(I,j)
!           !-----------------------------------------------------------
!           ! Calculate the evaporative term
!           !-----------------------------------------------------------
            IF (DB_DSCT(I,j)  >   0.0) THEN
              ZR = SQRT( ZC_DSC(I,j) / DSCDEPTH(I,j) )
              EVAP_TERM = A_ENT_2 * RHO_UV(I,j,K)                       &
     &                  * CHI_S_DSCT(I,j) * CHI_S_DSCT(I,j)             &
     &                  * ZR * ZR * ZR * DB_DSCT_CLD(I,j)               &
     &                  * SQRT( DSCDEPTH(I,j) * DB_DSCT(I,j) )
            ELSE
              EVAP_TERM = 0.0
            END IF
!             IF (CF(I,j,K-1) >= 0.9) THEN
            ENT_FACTOR = 1.0
!             ELSE
!              ENT_FACTOR = EXP(-((0.90-CF(I,j,K-1))**3.0)/0.075)
!             END IF
!           !-----------------------------------------------------------
!           ! Combine forcing terms to calculate the representative
!           ! velocity scales
!           !-----------------------------------------------------------
            V_SUM   = ( (SF_TERM + SF_SHEAR_TERM +                      &
     &                   ENT_FACTOR*(IR_TERM + EVAP_TERM))              &
     &              * DSCDEPTH(I,j) / (A_ENT_1*RHO_UV(I,j,K)) )**(1./3.)
            V_TOP_DSC(I,j) =( ENT_FACTOR * (IR_TERM + EVAP_TERM) *      &
     &                      DSCDEPTH(I,j) /                             &
     &                         (A_ENT_1*RHO_UV(I,j,K)) )**(1.0/3.0)
!           !-----------------------------------------------------------
!           ! Calculate the direct radiative term
!           !-----------------------------------------------------------
            IF (DB_DSCT(I,j)  >   0.0) THEN
             DZ_INV  = MIN( V_SUM*V_SUM / DB_DSCT(I,j) ,100.0 )
             L_RAD   = 15.0 * MAX( 1.0 , 200./(ZC_DSC(I,j)+1.0) )
             ALPHA_T = 1.0 - EXP(-0.5*DZ_INV/L_RAD)
             IF (D_SIEMS_DSC(I,j)  >   0.0) ALPHA_T =                   &
     &         MIN( 1.0, ALPHA_T + 10.0*D_SIEMS_DSC(I,j)*(1.0-ALPHA_T) )
             DR_TERM = BTC_TOP(I,j) * ALPHA_T * DF_DSCT_OVER_CP(I,j)
!            !----------------------------------------------------------
!            ! Finally combine terms to calculate the entrainment
!            ! mixing coefficients
!            !----------------------------------------------------------
             ZIL_CORR = C_T * ( (SF_TERM + SF_SHEAR_TERM +              &
     &                 ENT_FACTOR*(IR_TERM + EVAP_TERM)) /              &
     &                (RHO_UV(I,j,K) * SQRT(DSCDEPTH(I,j))) )**(2.0/3.0)
             RHOKH_DSCT_ENT(I,j) = ( SF_TERM + SF_SHEAR_TERM            &
     &               + ENT_FACTOR * (IR_TERM + EVAP_TERM + DR_TERM) )   &
     &                / ((DB_DSCT(I,j) + ZIL_CORR) * RDZ(I,j,K) )

!             RHOKM_TOP(I,j,K) = 0.75 * RHOKH_DSCT_ENT(I,j)
!     &                   * RHO_TQ(I,j,K-1) / RHO_UV(I,j,K)
             RHOKM_TOP(I,j,K) = 0.75 * RHOKH_DSCT_ENT(I,j)              &
     &                    * RDZ(I,j,K) * (Z_UV(I,j,K)-Z_UV(I,j,K-1))    &
     &                    * RHO_TQ(I,j,K-1) / RHO_UV(I,j,K)
            END IF   ! test on DB_DSCT gt 0
          END IF
        END DO
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
!        WBN_INT/WBP_INT = DEC_THRES (decoupling parameter),
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
!cdir collapse
      DO j=1,rows
      DO I=1,row_length
        DO K=1,BL_LEVELS
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
        kstatus(i,j)= .true.
        kwb0(i,j)  = 2
        ntop(i,j)  = -1
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
      do k=2,max_ntml
        do j=1, rows
        do i=1,row_length
          if (kstatus(i,j)) then
            if ( (db_ksurf_dry(i,J,k) <=  0.0) .or.                     &
     &           (k >= ntml(i,j)) ) then
                       kstatus(i,j)=.false.
                       kwb0(i,j)=k
            end if
          end if
        end do
        end do
      end do
!
! Set flags for iteratiing wb integral to calculate depth of mixing,
! one each for KSURF and K_TOP.  Note these will be updated depending on
! what happpens on testing for a well-mixed layer in section 2.2.
!
!cdir collapse
      DO j=1,rows
      DO I=1,row_length
        KSURF_ITERATE(I,j)= .FALSE.
        KTOP_ITERATE(I,j) = .FALSE.
        IF ( BFLUX_SURF(I,j)  >   0.0 .AND.                             &
!            ! Otherwise: surface mixing generated by stable scheme
     &       .NOT.CUMULUS(I,j) .AND.                                    &
!            ! Rule out CUMULUS layers from iteration of ZH
     &       NTDSC(I,j)  >   2                                          &
!            ! Otherwise: layer too shallow to resolve decoupling
     &     ) THEN
          IF ( DECFIX  ==  ON .OR.                                      &
!              ! Revised version
     &         ( DECFIX  ==  OFF .AND.                                  &
!              ! Standard version more restrictive:
     &           COUPLED(I,j) .AND.                                     &
!                ! Otherwise: use parcel ascent diagnosis of ZH
     &           V_TOP_DSC(I,j)  >   0.5                                &
!                ! Otherwise: DSC layer not turbulently active enough,
!                !            particularly without diagnosis of shear
!                !            production (ZH would be capped at
!                !            cloud-base without sufficient K_top to
!                !            sustain mixing into the cloud-layer)
     &         ) ) THEN
             KSURF_ITERATE(I,j)= .TRUE.
          END IF
        END IF
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
!          ! well-resolved DSC layer
           IF ( DECFIX  ==  ON .OR.                                     &
!             ! Revised version
     &        ( DECFIX  ==  OFF .AND. V_TOP_DSC(I,j)  >   0.5 )         &
!             ! Standard version requires turbulently active DSC
     &        ) THEN
!           ! ...test for recoupling with SML
            Z_INV(i,j)  = ZHSC(I,j)
            Z_CBASE(i,j)= Z_INV(i,j) - ZC_DSC(I,j)
            V_KTOP(i,j) = V_TOP_DSC(I,j)
            NTOP(i,j)   = NTDSC(I,j)
           END IF
          ELSE IF ( .NOT.DSC(I,j) .AND. .NOT.CUMULUS(I,j) .AND.         &
     &              NTML(I,j)  >   2) THEN
!          ! well-resolved SML
           IF ( DECFIX  ==  ON .OR.                                     &
!             ! Revised version
     &        ( DECFIX  ==  OFF .AND. V_TOP(I,j)  >   0.5 )             &
!             ! Standard version requires turbulently active SML
     &        ) THEN
!           ! ...test for decoupling
!           ! Note: code can only deal with one DSC layer at a time so
!           ! can't decouple SML if a DSC layer already exists.
            IF ( DECFIX  ==  ON .AND. ZC(I,j)  ==  0.0) THEN
!             ! If the BL is cloud-free then use a less restrictive
!             ! threshold - ideally, the parcel ascent would have
!             ! found the correct BL top in this case but this test is
!             ! kept to keep negative buoyancy fluxes under control
!             ! (ie. DEC_THRES_CLEAR=1 ensures wbn_int < |wbp_int|)
              DEC_THRES(I,j) = DEC_THRES_CLEAR
            END IF
            Z_INV(i,j)  = ZH(I,j)
            Z_CBASE(i,j)= Z_INV(i,j) - ZC(I,j)
            V_KTOP(i,j) = V_TOP(I,j)
            NTOP(i,j)   = NTML(I,j)
           END IF
          END IF
!         !----------------------------------------------------
!         ! estimate wb integral over surface layer
!         ! (and up to next theta-level, namely Z_TQ(KSURF) )
!         ! assuming linear profile going to zero at ZWB0
!         !----------------------------------------------------
          IF ( Kwb0(i,j)  ==  NTML(I,j) ) THEN
            ZWB0(i,j) = ZH(I,j)
          ELSEIF ( Kwb0(i,j)  ==  2 ) THEN
            ZWB0(i,j) = Z_UV(I,j,2)
          ELSE
            K=Kwb0(i,j)
!           ! now DB_KSURF_DRY(K) LE 0 and DB_KSURF_DRY(K-1) GT 0
!           ! so interpolate:
            DB_RATIO = DB_KSURF_DRY(I,j,K-1)                            &
     &                  /(DB_KSURF_DRY(I,j,K-1)                         &
     &                  -DB_KSURF_DRY(I,j,K))
            DB_RATIO = MAX( 0.0, DB_RATIO )  ! trap for rounding error
            ZWB0(i,j)=Z_UV(I,j,K-1)+DB_RATIO*                           &
     &                  (Z_UV(I,j,K)-Z_UV(I,j,K-1))
          END IF
          WB_SURF_INT(i,j) = BFLUX_SURF(I,j) * Z_TQ(I,j,KSURF(I,j)) *   &
     &                  ( 1.0 - Z_TQ(I,j,KSURF(I,j))/(2.0*ZWB0(i,j)))
          WB_SURF_INT(i,j) = MAX( 1.0E-14, WB_SURF_INT(i,j) )
        ELSE
!         ! only include surface layer contribution for unstable mixing
          WB_SURF_INT(i,j) = 1.0E-14
        END IF

        WBP_INT(i,j) = WB_SURF_INT(i,j)  ! must be > 0
        WBN_INT(i,j) = 0.0

      END DO ! I
      END DO ! J

      ksurf_min=minval(KSURF)
      ntop_max=maxval(NTOP)

      DO K=ksurf_min+1,ntop_max+1
        DO J=1, rows
        DO I=1, row_length

         IF ( Z_INV(i,j)  >   0.0 ) THEN
!         ! ----------------------------------------------
!         ! worth testing layer as well-mixed to cloud-top
!         ! ----------------------------------------------
          WBMIX(I,j,KSURF(I,j)) = WB_SURF_INT(i,j)
          WBEND(I,j,KSURF(I,j)) = WB_SURF_INT(i,j)

          ZB_KTOP = 0.1*Z_INV(i,j)
          ZINV_PR(i,j) = Z_INV(i,j) - ZB_KTOP
!         ! DB(K)is the K to K-1 difference and already
!         ! integrated up to K_SURF, so start this loop at KSURF+1
          IF ( (K >= KSURF(I,j)+1) .and.  (K <= NTOP(i,j)+1) ) THEN

            KHTOP(i,j) = 0.0
            KHSURF(i,j) = 0.0
            Z_PR = Z_UV(I,j,K) - ZB_KTOP
            IF (Z_PR  >   0.0 .AND. Z_PR  <   ZINV_PR(i,j)) THEN
              KHTOP(i,j)=VKMAN*RHO_UV(I,j,K)*0.85*V_KTOP(i,j)*          &
     &                       (( 1.0 - Z_PR/ZINV_PR(i,j) )**0.8)         &
     &                             * Z_PR * Z_PR / ZINV_PR(i,j)
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

            IF (Z_CBASE(i,j)  >   Z_TQ(I,j,K)) THEN
!             ! cloud-base above this integration range so use dry WB
              WB_SCLD = ( KHSURF(i,j) * DB_KSURF_DRY(I,j,K) +           &
     &                    KHTOP(i,j) * DB_KTOP_DRY(I,j,K) )
              WB_CLD  = 0.0
            ELSEIF (Z_CBASE(i,j)  <   Z_TQ(I,j,K-1)) THEN
!             ! cloud-base below this integration range so use cloudy WB
              WB_CLD = ( KHSURF(i,j) * DB_KSURF_CLD(I,j,K) +            &
     &                   KHTOP(i,j)  * DB_KTOP_CLD(I,j,K) )
              WB_SCLD = 0.0
            ELSE
!             ! cloud-base within this integration range
!             ! so treat cloud and sub-cloud layer wb separately
              CLD_FRAC = (Z_TQ(I,j,K)-Z_CBASE(i,j))                     &
     &                  /(Z_TQ(I,j,K)-Z_TQ(I,j,K-1))
              WB_CLD = CLD_FRAC                                         &
     &                * ( KHSURF(i,j) * DB_KSURF_CLD(I,j,K) +           &
     &                    KHTOP(i,j) * DB_KTOP_CLD(I,j,K) )
              WB_SCLD = (1.0-CLD_FRAC) *                                &
     &                ( KHSURF(i,j) * DB_KSURF_DRY(I,j,K) +             &
     &                  KHTOP(i,j) * DB_KTOP_DRY(I,j,K) )
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

         END IF ! Z_INV GT 0
        END DO ! I
        END DO ! J
      END DO ! K
!
! Test WB_Ratio to see if layer should be well-mixed
!
!cdir collapse
      DO J=1,rows
      DO I=1,row_length
        IF ( Z_INV(i,j)  >   0.0 ) THEN
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
              ZC_DSC(I,j) = 0.0
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
!             ! Set up a `COUPLED' decoupled layer
!             ! Note a new ZH (and thence NTML) will be calculated by
!             ! wb integral iteration.
              IF (CUMULUS(I,j)) ZK_UV=SQRT(ZH(I,j)-1000000.)
!                             ! APLTEST: shouldn't ever happen!
              DSC(I,j) = .TRUE.
              COUPLED(I,j) = .TRUE.  ! as full entrainment applied at ZH
              NTDSC(I,j) = NTML(I,j)
              ZHSC(I,j) = ZH(I,j)
              ZC_DSC(I,j) = ZC(I,j)
              V_TOP_DSC(I,j) = V_TOP(I,j)
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
        END IF   ! testing for well-mixed layer (Z_INV GT 0)

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
          Z_BOT_LIM(i,j)=Z_UV(I,j,KSURF(I,j)+1)                         &
     &           + 0.1 * (Z_UV(I,j,KSURF(I,j)+2)-Z_UV(I,j,KSURF(I,j)+1))
!         ! limit K-surf to below cloud-base
          Z_TOP_LIM(i,j)=MAX(Z_BOT_LIM(i,j),ZHSC(I,j) )

          Z_CBASE(i,j) = ZHSC(I,j) - ZC_DSC(I,j)
!         !-----------------------------------------------------
!         ! Initial increment to ZSML_TOP found by dividing
!         ! up depth of layer within which it is allowed:
!         ! Start with ZSML_TOP at lower limit and work upwards
!         !-----------------------------------------------------
          Z_INC(i,j)=(Z_TOP_LIM(i,j)-Z_BOT_LIM(i,j))                    &
     &                 / FLOAT(N_STEPS)
          ZSML_TOP(I,j) = Z_BOT_LIM(i,j)

          N_SWEEP = 1
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

      ksurf_min=minval(KSURF)
      ntdsc_max=maxval(NTDSC)

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
            WBP_INT(i1,j1) = WB_SURF_INT(i1,j1)  ! must be > 0
            WBN_INT(i1,j1) = 0.0

            Z_INV(i1,j1) = ZSML_TOP(i1,j1)

          end do ! ic c_len_i
!
!..Integrate buoyancy flux profile given this ZSML_TOP
!
          do k=ksurf_min+1,ntdsc_max+1
!cdir nodep
            do ic=1, c_len_i
              j1=(ind_todo_i(ic)-1)/row_length+1
              i1=ind_todo_i(ic)-(j1-1)*row_length

              IF ( K  >=  KSURF(i1,j1)+1 .and.                          &
     &             K  <=  NTDSC(i1,j1)+1 ) then

                Z_PR = Z_UV(i1,j1,K)
                IF (Z_PR  <   Z_INV(i1,j1)) THEN
                  KH_SURF(i1,j1,K) = W_H_TOP(i1,j1)                     &
     &                                  *Z_PR*(1.0-Z_PR/Z_INV(i1,j1) )  &
     &                                       *(1.0-Z_PR/Z_INV(i1,j1) )
                ELSE
                  KH_SURF(i1,j1,K) = 0.0
                END IF

                IF (Z_CBASE(i1,j1)  >   Z_TQ(i1,j1,K)) THEN
!                 ! cloud-base above this range so use dry WB
                  WB_SCLD= KH_SURF(i1,j1,K) * DB_KSURF_DRY(i1,j1,K)
                  WB_CLD = 0.0
                ELSEIF (Z_CBASE(i1,j1)  <   Z_TQ(i1,j1,K-1)) THEN
!                 ! cloud-base below this range so use cloudy WB
                  WB_CLD = KH_SURF(i1,j1,K) * DB_KSURF_CLD(i1,j1,K)
                  WB_SCLD=0.0
                ELSE
!                 ! cloud-base within this integration range
!                 ! so treat cloud and sub-cloud layer wb separately
                  CLD_FRAC = (Z_TQ(i1,j1,K)-Z_CBASE(i1,j1))             &
     &                      /(Z_TQ(i1,j1,K)-Z_TQ(i1,j1,K-1))
                  WB_CLD  = CLD_FRAC                                    &
     &                       * KH_SURF(i1,j1,K)*DB_KSURF_CLD(i1,j1,K)
                  WB_SCLD = (1.0-CLD_FRAC)                              &
     &                       * KH_SURF(i1,j1,K)*DB_KSURF_DRY(i1,j1,K)
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
        ntop(I,j)=2
        status_ntop(I,j)=.true.
      END DO
      END DO

      DO K=2,BL_LEVELS-2
!cdir collapse
        DO J=1, rows
        DO I=1, row_length
          IF ( KSURF_ITERATE(I,j) .and. status_ntop(i,j) ) THEN
!           ! -------------
!           ! find new NTML
!           ! -------------
            IF  (Z_UV(I,j,K+1)  <   ZSML_TOP(I,j)) THEN
              NTOP(i,j) = K+1
            ELSE
              status_ntop(i,j)=.false.
            END IF
!           ! --------------------------------------------------------
!           ! Rounding error previously found to give
!           !      ZSML_TOP > Z_TOP_LIM = ZHSC
!           ! Test on ZSML_TOP hitting thresholds consequently changed
!           ! but also include the following failsafe tests here.
!           ! --------------------------------------------------------
            NTML(I,j) = MIN( NTDSC(I,j), NTOP(i,j)-1 )
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

            Z_BOT_LIM(i,j) = 0.1 * ZSML_TOP(I,j)
!           ! Ensure mixing extends slightly below base of level NTDSC
            Z_TOP_LIM(i,j) = Z_UV(I,j,NTDSC(I,j)) - 0.1 *               &
     &                  ( Z_UV(I,j,NTDSC(I,j))-Z_UV(I,j,NTDSC(I,j)-1) )
!           ! Limit base of top-driven mixing to above ZH if cumulus
            IF ( CUMULUS(I,j) ) THEN
              Z_BOT_LIM(i,j) = Z_UV(I,j,NTML(I,j)+1)
              IF (Z_TOP_LIM(i,j) <  Z_BOT_LIM(i,j) )                    &
     &             Z_BOT_LIM(i,j) = Z_TOP_LIM(i,j)
            END IF

            Z_CBASE(i,j) = ZHSC(I,j) - ZC_DSC(I,j)

!..Divide up depth of layer within which ZDSC_BASE is allowed
            Z_INC(i,j)=(Z_TOP_LIM(i,j)-Z_BOT_LIM(i,j))                  &
     &                   /FLOAT(N_STEPS)
            ZDSC_BASE(I,j) = Z_BOT_LIM(i,j)
!                            ! will start at Z_BOT_LIM+Z_INC

            N_SWEEP = 1
            WB_RATIO(i,j) = DEC_THRES(i,j) + 1.0 ! to be > DEC_THRES

        END IF ! KTOP_ITERATE

      END DO
      END DO

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
            WBN_INT(i1,j1) = 0.0
            WBP_INT(i1,j1) = 1.E-14
            WBEND(i1,j1,KSURF(i1,j1)) = 0.0
            IF ( KSURF_ITERATE(i1,j1) .AND.                             &
     &           ZDSC_BASE(i1,j1)  <   ZSML_TOP(i1,j1) ) THEN
!             ! only include surface flux if K_SURF is included
!             ! in the wb calculation and K profiles overlap
              WBP_INT(i1,j1) = WB_SURF_INT(i1,j1)
              WBN_INT(i1,j1) = 0.0
              WBEND(i1,j1,KSURF(i1,j1)) = WB_SURF_INT(i1,j1)
            END IF

            ZINV_PR(i1,j1) = ZHSC(i1,j1)-ZDSC_BASE(i1,j1)

          end do ! ic c_len_i
!
!..Integrate buoyancy flux profile given this ZDSC_BASE
!
          do k=ksurf_min+1,ntdsc_max+1
!cdir nodep
            do ic=1, c_len_i
              j1=(ind_todo_i(ic)-1)/row_length+1
              i1=ind_todo_i(ic)-(j1-1)*row_length

              IF ((K >= KSURF(i1,j1)+1).and.(K <= NTDSC(i1,j1)+1)) THEN
                Z_PR = Z_UV(i1,j1,K) - ZDSC_BASE(i1,j1)
                KHTOP(i1,j1) = 0.0
                IF (Z_PR >   0.0 .AND.Z_PR <  ZINV_PR(i1,j1)) THEN
                  KHTOP(i1,j1) = 0.85 * V_TOP_DSC(i1,j1) *              &
     &                           (( 1.0 - Z_PR/ZINV_PR(i1,j1) )**0.8)   &
     &                         * Z_PR * Z_PR / ZINV_PR(i1,j1)
                END IF
                KHSURF(i1,j1) = 0.0
                IF ( ZDSC_BASE(i1,j1)  <   ZSML_TOP(i1,j1) ) THEN
!                 ! only include K_surf if profiles overlap
!                 ! otherwise layers are independent
                  KHSURF(i1,j1) = KH_SURF(i1,j1,K)
                END IF

                IF (Z_CBASE(i1,j1)  >   Z_TQ(i1,j1,K)) THEN
!                 ! cloud-base above this range so use dry WB
                  WB_SCLD=( KHSURF(i1,j1)*DB_KSURF_DRY(i1,j1,K)+        &
     &                      KHTOP(i1,j1) *DB_KTOP_DRY(i1,j1,K) )
                  WB_CLD = 0.0
                ELSEIF (Z_CBASE(i1,j1)  <   Z_TQ(i1,j1,K-1)) THEN
!                 ! cloud-base below this range so use cloudy WB
                  WB_CLD = ( KHSURF(i1,j1)*DB_KSURF_CLD(i1,j1,K) +      &
     &                       KHTOP(i1,j1)*DB_KTOP_CLD(i1,j1,K) )
                  WB_SCLD = 0.0
                ELSE
!                 ! cloud-base within this integration range
!                 ! so treat cloud and sub-cloud layer wb separately
                  CLD_FRAC= (Z_TQ(i1,j1,K)-Z_CBASE(i1,j1))              &
     &                     /(Z_TQ(i1,j1,K)-Z_TQ(i1,j1,K-1))
                  WB_CLD  = CLD_FRAC *                                  &
     &                   ( KHSURF(i1,j1) * DB_KSURF_CLD(i1,j1,K) +      &
     &                     KHTOP(i1,j1)  * DB_KTOP_CLD(i1,j1,K) )
                  WB_SCLD = (1.0-CLD_FRAC) *                            &
     &                   ( KHSURF(i1,j1) * DB_KSURF_DRY(i1,j1,K) +      &
     &                     KHTOP(i1,j1) * DB_KTOP_DRY(i1,j1,K) )
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
!            ! hit thres while working up
     &       (UP(l) == 0 .AND.WB_RATIO(i1,j1) >= DEC_THRES(i1,j1))) THEN
!            ! hit thres while working down
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
!     Calculate factors required to ensure that the non-local turbulent
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
!! 2.  Loop around levels again calculating height dependent turbulent
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
!         !-------------------------------------------------------------
!         ! Calculate RHOK(H/M)_TOP, top-down turbulent mixing profiles
!         ! for the surface mixed layer.
!         ! This is a variation on an up-side-down version of the cubic
!         ! surface-forced profiles below.  Implement between at least
!         ! the top of the `surface layer' (at Z=0.1*ZH) and ZH.
!         ! Note this may well include NTML+1: entrainment fluxes will
!         ! be dealt with in KMKHZ.
!         !-------------------------------------------------------------
!
          IF ( ZK_UV  <   ZH(I,j) .AND.                                 &
     &         ZK_UV  >   ZSML_BASE(I,j) ) THEN
            Z_PR  = ZK_UV - ZSML_BASE(I,j)
            ZH_PR = ZH(I,j) - ZSML_BASE(I,j)
            RHOKH_TOP(I,j,K) = RHO_UV(I,j,K) * V_TOP(I,j) * 0.85 *      &
     &          VKMAN * ( ( 1.0 - KH_SCT_FACTOR(I,j)*Z_PR/ZH_PR )**0.8 )&
     &                                         * Z_PR * Z_PR / ZH_PR
          END IF
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
              RHOKH_TOP(I,j,K) = RHOKH_TOP(I,j,K) +                     &
     &           RHO_UV(I,j,K)*V_TOP_DSC(I,j)*0.85*VKMAN*               &
     &              ( ( 1.0 - KH_DSCT_FACTOR(I,j)*Z_PR/ZH_PR )**0.8 )   &
     &                                         * Z_PR * Z_PR / ZH_PR
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

          IF (FB_SURF(I,j)  >=  0.0) THEN
!
!           Calculate the free-convective scaling velocity at z(k)
!
            IF (ZK_UV  <=  0.1*ZH(I,j)) THEN
!
!             Surface layer calculation
!
              W_S_CUBED_UV = 2.5 * ZK_UV * FB_SURF(I,j)
            ELSE
!
!             Outer layer calculation
!
              IF (COUPLED(I,j)) THEN  !  coupled and cloudy
               W_S_CUBED_UV = 0.25 * ZHSC(I,j) * FB_SURF(I,j)
              ELSE
               W_S_CUBED_UV = 0.25 * ZH(I,j) * FB_SURF(I,j)
              END IF
            END IF

            IF (ZK_TQ  <=  0.1*ZH(I,j)) THEN
!
!             Surface layer calculation
!
              W_S_CUBED_TQ = 2.5 * ZK_TQ * FB_SURF(I,j)
            ELSE
!
!             Outer layer calculation
!
              IF (COUPLED(I,j)) THEN  !  coupled and cloudy
                W_S_CUBED_TQ = 0.25 * ZHSC(I,j) * FB_SURF(I,j)
              ELSE
              W_S_CUBED_TQ = 0.25 * ZH(I,j) * FB_SURF(I,j)
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
            PRANDTL = 0.75 * ( V_S(I,j)*V_S(I,j)*V_S(I,j)*V_S(I,j) +    &
     &                   (4.0/25.0)*W_S_CUBED_UV*W_M_UV ) /             &
     &                       ( V_S(I,j)*V_S(I,j)*V_S(I,j)*V_S(I,j) +    &
     &                   (8.0/25.0)*W_S_CUBED_UV*W_M_UV )
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
                  W_S_CUBED_TQ = 0.25 * ZHSC(I,j) * FB_SURF(I,j)
                ELSE
                  W_S_CUBED_TQ = 0.25 * ZH(I,j) * FB_SURF(I,j)
                END IF


! 4*W_S_CUBED_TQ = the convective boundary layer
! velocity scale cubed
! V_S = the friction velocity

                F_NGSTRESS(I,j,K) = ( RHO_TQ(I,j,K-1)/RHOSTAR_GB(I,j) ) &
     &            * S_M * ( A_NGS * 4.0*W_S_CUBED_TQ/                   &
     &             (V_S(I,j)*V_S(I,j)*V_S(I,j) + W_S_CUBED_TQ*4.*0.6 ) )&
     &               * ( Z_PR / ZH_PR ) * ( 1.0 -  ( Z_PR / ZH_PR ) ) * &
     &                                    ( 1.0 -  ( Z_PR / ZH_PR ) )

              END IF
            END IF
          END DO
          END DO
        END DO
      END IF

      DO K=1,BL_LEVELS
        DO j=1,rows
        DO I=1,row_length
!         ! convert to m2/s-3
          IF (K >  KSURF(I,j)) THEN
            WBMIX(I,j,K)=WBMIX(I,j,K)*RDZ(I,j,K)
            WBEND(I,j,K)=WBEND(I,j,K)*RDZ(I,j,K)
          ELSE
            GAMMA_WBS = ( (WBMIX(I,j,KSURF(I,j))/Z_TQ(I,j,KSURF(I,j)))  &
     &                    - BFLUX_SURF(I,j)  )*2.0/Z_TQ(I,j,KSURF(I,j))
            WBMIX(I,j,K) = BFLUX_SURF(I,j) + GAMMA_WBS*Z_UV(I,j,K)

            GAMMA_WBS = ( (WBEND(I,j,KSURF(I,j))/Z_TQ(I,j,KSURF(I,j)))  &
     &                    - BFLUX_SURF(I,j)  )*2.0/Z_TQ(I,j,KSURF(I,j))
            WBEND(I,j,K) =  BFLUX_SURF(I,j) + GAMMA_WBS*Z_UV(I,j,K)
          END IF
        END DO
        END DO
      END DO
      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('EXCF_NL  ',4)
      END IF
      RETURN
      END SUBROUTINE EXCF_NL
#endif
