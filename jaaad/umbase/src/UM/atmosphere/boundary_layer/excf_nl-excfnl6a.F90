#if defined(A03_8A)
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
!!!  4.4    Feb 1997  Written by R N B Smith
!!!  5.1    Apr 2000  Shear-driven entrainment removed if surface is
!!!                   stable.     A.P.Lock
!!!  5.2    Mar 2000  New diagnosis of decoupling/recoupling,
!!!                   APLHA_T parametrized,
!!!                   changes resulting from the introduction of
!!!                   subgrid ZH,ZHSC.    A.P.Lock
!!!  5.5    Feb 2003  Calculate non-gradient stress function A.P.Lock
!!!  6.2    Aug 2005  Add spaces to DO WHILE. P.Selwood
!!!   6.2   8/02/06   Correction to non-gradient stress to
!!!                   follow Brown and Grant 97 and include
!!!                   density properly.   Bob Beare

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
     & NG_STRESS,F_NGSTRESS,                                            &
     & ZHSC,DSCDEPTH,NTDSC,DB_DSCT,SVL,                                 &
     & BT_DSCT,BTT_DSCT,                                                &
     & DF_DSCT_OVER_CP,ZETA_R_DSC,BTC_DSCT,                             &
     & DB_DSCT_CLD,CHI_S_DSCT,ZC_DSC,COUPLED,                           &
     & D_SIEMS,D_SIEMS_DSC,NBDSC,                                       &
     & DB_FOR_FLUX,DSC,CUMULUS,ZDSC_BASE,                               &
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
     &,NG_STRESS   ! IN switch for non-gradient stress

      LOGICAL                                                           &
     & COUPLED(row_length,rows)                                         &
                                ! IN  Flag to indicate Sc layer weakly
!                               !     coupled to surface (ie weakly
                                !     decoupled)
     &,CUMULUS(row_length,rows)                                         &
                                ! INOUT Flag for cumulus
     &,DSC(row_length,rows)     ! INOUT Flag set if decoupled stra
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
                                !    for an unstaggered grid the
                                !    heights of the half levels
                                !    1.5 to BL_LEVELS+0.5 should be
                                !    input to elements 1 to BL_LEVELS.
                                !    (value for BL_LEVELS not used
                                !    in either case.)

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
!*APL*ALPHA
     &,D_SIEMS(row_length,rows)                                         &
                                 ! IN Siems (1990) et al. cloud-top
                                 !    entr.t instab. parm
     &,D_SIEMS_DSC(row_length,rows)                                     &
                                 !IN Siems (1990) et al. cloud-top
                                 !   entr.t instab. parm for DSC layer
     &,DB_FOR_FLUX(row_length,rows,2:BL_LEVELS) ! Buoyancy jump in flux
!                                     ! calculation (m/s2)
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
#include "c_lheat.h"
#include "c_r_cp.h"
#include "c_epslon.h"
#include "c_vkman.h"
      REAL A_ENT_1,A_ENT_2,C_T,A_ENT_SHR,DEC_THRES
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
     &,DEC_THRES=0.1                                                    &
                                ! Decoupling threshold (larger makes
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
     &,RECOUPLED(row_length,rows)  ! Flag to signal layer recoupled

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
     &,DB_SCALED(row_length,rows,2:BL_LEVELS)                           &
                                ! scaled buoyancy flux fac
     &,ZSML_TOP(row_length,rows)                                        &
                                ! Height of top of surf-driven K in S
     &,ZSML_BASE(row_length,rows)                                       &
                                ! Height of base of top-driven K in S
     &,V_SURF(row_length,rows)                                          &
                                ! Velocity scale for surface-up conve
     &,KH_SMLT_SHAPE(row_length,rows,2:BL_LEVELS)                       &
                                ! Shape factor for non
!                               ! turbulent mixing coefficient
     &,SCDEPTH(row_length,rows)
                                ! Depth of top-driven mixing in
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
     &,ZH_M                                                             &
                    ! Height taken as top of the turbulent mixing layer
!                   ! for momentum.
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
! *APL*U*
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
                    ! Height of cloud layer top above surface layer
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
!*APL*ALPHA
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
     &,RHOKH_SHAPE                                                      &
                     ! RHOKH shape factor
     &,WB                                                               &
                     ! Buoyancy flux
     &,WBN_INT                                                          &
                     ! Negative part of buoyancy flux integral
     &,WBP_INT                                                          &
                     ! Positive part of buoyancy flux integral
     &,WB_RATIO                                                         &
                     ! WBN_INT/WBP_INT
     &,Z_TOP_LIM                                                        &
                     ! upper height limit on K profile
     &,Z_BOT_LIM                                                        &
                     ! lower height limit on K profile
     &,Z_INC         ! Step size (m)
!
      INTEGER                                                           &
     & I,j                                                              &
                    ! Loop counter (horizontal field index).
     &,K                                                                &
                    ! Loop counter (vertical level index).
     &,UP                                                               &
                    ! indicator of upward/downward sweep
     &,N_SWEEP                                                          &
                    ! sweep counter
     &,NS                                                               &
                    ! step counter
     &,KTOP         ! top level of surf-driven K profile
!

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('EXCF_NL  ',3)
      ENDIF
!
!-----------------------------------------------------------------------
!! 0.  Calculate top-of-b.l. velocity scales and Prandtl number.
!-----------------------------------------------------------------------
!
       do j=1,rows
       do i=1,row_length

         RHOKH_SURF_ENT(I,j) = 0.0
         RHOKH_TOP_ENT(I,j) = 0.0
         RHOKH_DSCT_ENT(I,j) = 0.0
         V_TOP(I,j) = 0.0
         V_SURF(I,j) = 0.0
         V_TOP_DSC(I,j) = 0.0
         SCBASE(I,j)=.FALSE.

         IF (FB_SURF(I,j)  >=  0.0) THEN
!
!         By definition the top of the b.l. is in the 'outer layer' so
!         the free-convective velocity scale cubed is
!
          IF (COUPLED(I,j)) THEN
           W_S_CUBED_UV = 0.25 * ZHSC(I,j) * FB_SURF(I,j)
          ELSE
           W_S_CUBED_UV = 0.25 * ZH(I,j) * FB_SURF(I,j)
          ENDIF
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
        ENDIF
       ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!! 1.  Loop round levels; calculate the top-of-b.l. entrainment
!!     mixing coefficients.
!-----------------------------------------------------------------------
!
      DO K=2,BL_LEVELS
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
           KH_SMLT_SHAPE(I,j,K) = 0.0
           IF ( K  ==  NTML(I,j)+1 .AND. .NOT.COUPLED(I,j) ) THEN
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
            ENDIF
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
            ENDIF
!            IF (CF(I,j,K-1) >= 0.9) THEN
              ENT_FACTOR = 1.0
!            ELSE
!              ENT_FACTOR = EXP(-((0.90-CF(I,j,K-1))**3.0)/0.075)
!            ENDIF
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
            ENDIF    ! test on DB_TOP GT 0
          ENDIF
!----------------------------------------------------------------
!      THEN the top of the DSC (if coupled use ZHSC length-scale)
!----------------------------------------------------------------
          IF ( (K  ==  NTDSC(I,j)+1) ) THEN
             IF (COUPLED(I,j)) THEN
!              !--------------------------------------------------------
!              ! Calculate the surface buoyancy flux term
!              !--------------------------------------------------------
               ZETA_S_FAC = (1.0 - ZETA_S(I,j)) * (1.0 - ZETA_S(I,j))
               SF_TERM = A_ENT_1 * MAX ( 0.0 ,                          &
     &                          ( (1.0 - ZETA_S_FAC) * BFLUX_SURF(I,j)  &
     &                            + ZETA_S_FAC * BFLUX_SURF_SAT(I,j) ) )
!              !--------------------------------------------------------
!              ! Calculate the surface shear term
!              !--------------------------------------------------------
               IF (FB_SURF(I,j)  >=  0.0) THEN
               SF_SHEAR_TERM = A_ENT_SHR * V_S(I,j)*V_S(I,j)*V_S(I,j)   &
     &                          * RHO_UV(I,j,K)  / ZHSC(I,j)
             ELSE
                 SF_SHEAR_TERM = 0.0
               ENDIF
               V_SURF(I,j) = ( (SF_TERM) * ZHSC(I,j)                    &
     &                            / (A_ENT_1*RHO_UV(I,j,K)) )**(1.0/3.0)
             ELSE
               SF_TERM = 0.0
               SF_SHEAR_TERM = 0.0
             ENDIF
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
            EVAP_TERM = A_ENT_2 * RHO_UV(I,j,K)                         &
     &                  * CHI_S_DSCT(I,j) * CHI_S_DSCT(I,j)             &
     &                  * ZR * ZR * ZR * DB_DSCT_CLD(I,j)               &
     &                  * SQRT( DSCDEPTH(I,j) * DB_DSCT(I,j) )
            ELSE
             EVAP_TERM = 0.0
            ENDIF
!             IF (CF(I,j,K-1) >= 0.9) THEN
              ENT_FACTOR = 1.0
!             ELSE
!              ENT_FACTOR = EXP(-((0.90-CF(I,j,K-1))**3.0)/0.075)
!             ENDIF
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
            IF (D_SIEMS_DSC(I,j)  >   0.0) ALPHA_T =                    &
     &         MIN( 1.0, ALPHA_T + 10.0*D_SIEMS_DSC(I,j)*(1.0-ALPHA_T) )
             DR_TERM = BTC_TOP(I,j) * ALPHA_T * DF_DSCT_OVER_CP(I,j)
!           !-----------------------------------------------------------
!           ! Finally combine terms to calculate the entrainment
!           ! mixing coefficients
!           !-----------------------------------------------------------
            ZIL_CORR = C_T * ( (SF_TERM + SF_SHEAR_TERM +               &
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
            ENDIF   ! test on DB_DSCT gt 0
          ENDIF
        ENDDO
        ENDDO
      ENDDO
!
! ----------------------------------------------------------------------
! 2.  Estimate the depths of top-down and surface-up mixing.
!     These amount to diagnoses of recoupling and decoupling.
!     The K_H profiles are applied over layers such that the ratio
!        WBN_INT/WBP_INT = DEC_THRES (parameter set to ~0.1),
!     where WBN_INT and WBP_INT are the magnitudes of the integrals of
!     the negative and positive parts, respectively, of the resulting
!     buoyancy flux profile (given by - KH * DB_FOR_FLUX).
!     The only other conditions applied are that the base of any
!     decoupled layer overlying cumulus must be above Z_UV(NTML+2),
!     unless recoupling is diagnosed (both KH profiles are applied
!     from the surface to the top of the decoupled layer),
!     and that mixing must extend over at least two levels.
!     Surface-driven mixing is applied at least up to the base of
!     any top-down profile, but higher if the WB ratio allows.
! ----------------------------------------------------------------------
! 2.1 First consider decoupled layers and test for recoupling
!     (ie. both KH profiles extending from ZHSC to the surface) if
!     top-driven mixing is strong enough (V_TOP_DSC > 0.5, say).
!     Note that if the surface buoyancy flux is negative, V_TOP is zero
!     and the surface driven KH makes no contribution to wb.
! ----------------------------------------------------------
        do j=1,rows
        do i=1,row_length
        RECOUPLED(I,j) = .FALSE.    ! default
        IF ( V_TOP_DSC(I,j)  <=  0.0 ) THEN
          DSC(I,j) = .FALSE.
          NTDSC(I,j) = 0
          ZHSC(I,j) = 0.0
          DSCDEPTH(I,j) = 0.0
        ENDIF
      ENDDO
      ENDDO
      DO K=2,BL_LEVELS
      do j=1,rows
      do i=1,row_length
        IF ( DSC(I,j) ) THEN
            DB_SCALED(I,j,K) = - 0.85*VKMAN*RHO_UV(I,j,K)               &
     &                              *V_TOP_DSC(I,j)*DB_FOR_FLUX(I,j,K)
        ENDIF
      ENDDO
      ENDDO
      ENDDO

      do j=1,rows
      do i=1,row_length
        IF ( NTDSC(I,j)  ==  0 ) THEN
          ZDSC_BASE(I,j) = 0.0
        ELSEIF ( NTDSC(I,j)  <=  2 ) THEN
          ZDSC_BASE(I,j) = 0.1*ZH(I,j)
           ELSE
        WB = BFLUX_SURF(I,j) / RDZ(I,j,1)
        IF ( V_TOP_DSC(I,j)  >   0.5 .AND. V_SURF(I,j)  >   0.0 ) THEN
!       ! significant top and surface forcing so test for recoupling
          RECOUPLED(I,j) = .TRUE.    ! only indicates testing here
          Z_BOT_LIM = 0.1 * ZH(I,j)  ! in the end, KH_TOP must not
!                                  ! extend into the surface layer
          ZDSC_BASE(I,j) = 0.0
          NBDSC(I,j) = 2
          WBP_INT = WB
          WBN_INT = 0.0
            ELSE
!       ! Not safe to look for recoupling - start from maximum depth
          Z_BOT_LIM = 0.1 * ZH(I,j)
          NBDSC(I,j) = 2
!       ! Ensure depth of mixing is at least two layers:
          Z_TOP_LIM = Z_UV(I,j,NTDSC(I,j)-1)
!       ! Limit base of top-driven mixing to above ZH if cumulus-capped
          IF ( CUMULUS(I,j) ) THEN
            Z_BOT_LIM = Z_UV(I,j,NTML(I,j)+1)
            NBDSC(I,j) = NTML(I,j)+2
            IF ( Z_TOP_LIM  <   Z_BOT_LIM ) Z_TOP_LIM = Z_BOT_LIM
            ENDIF
          ZDSC_BASE(I,j) = Z_BOT_LIM
          WBN_INT = 0.0
          WBP_INT = 1.E-14
           ENDIF

        ZH_PR = ZHSC(I,j) - ZDSC_BASE(I,j)
        DO K = NBDSC(I,j), NTDSC(I,j)
         Z_PR = Z_UV(I,j,K) - ZDSC_BASE(I,j)
         IF (Z_PR  >   0.0 .AND. Z_PR  <   ZH_PR) THEN
           RHOKH_SHAPE = (( 1.0 - Z_PR/ZH_PR )**0.8)                    &
     &                                * Z_PR * Z_PR / ZH_PR
           IF ( RECOUPLED(I,j) )                                        &
!          !----------------------------------------
!          ! include surface-driven profile but
!          ! recall DB_SCALED is scaled on V_TOP_DSC
!          !----------------------------------------
     &       RHOKH_SHAPE = RHOKH_SHAPE + ( V_SURF(I,j)/V_TOP_DSC(I,j) ) &
     &                                    * Z_PR * ( 1.0 - Z_PR/ZH_PR ) &
     &                                           * ( 1.0 - Z_PR/ZH_PR )
           WB = RHOKH_SHAPE * DB_SCALED(I,j,K)
         ELSE
           WB = 0.0
          ENDIF

         IF (WB  >=  0.0) THEN
           WBP_INT = WBP_INT + WB
         ELSE
           WBN_INT = WBN_INT - WB
         ENDIF
        ENDDO

        WB_RATIO = WBN_INT/WBP_INT

        IF (WB_RATIO  <   DEC_THRES) THEN
!..No need to test depth of mixing any further - maximum permissable
!..extent of mixing satisfies buoyancy flux integral criteria

          IF ( RECOUPLED(I,j) ) THEN
!           ! move/adjust surface driven entrainment
!           ! RHOKH(z_i) = rho * w_e * DZL and w_e ~ 1/DB_TOP, so:
            IF ( DB_TOP(I,j) >  0.0 .AND. DB_DSCT(I,j) >  0.01 ) THEN
!                                         ! can't calc Zil. term
              RHOKH_SURF_ENT(I,j) = RHOKH_SURF_ENT(I,j) *               &
     &         ( RHO_UV(I,j,NTDSC(I,j)+1) * DB_TOP(I,j)                 &
     &                                        * RDZ(I,j,NTML(I,j)+1) ) /&
     &         ( RHO_UV(I,j,NTML(I,j)+1) * DB_DSCT(I,j)                 &
     &                                        * RDZ(I,j,NTDSC(I,j)+1) )

              RHOKM(I,j,NTDSC(I,j)+1) = RHOKM(I,j,NTML(I,j)+1) *        &
     &         (RHO_TQ(I,j,NTDSC(I,j)) * DB_TOP(I,j)                    &
     &               * (Z_UV(I,j,NTDSC(I,j)+1)-Z_UV(I,j,NTDSC(I,j))) ) /&
     &         (RHO_TQ(I,j,NTML(I,j)) * DB_DSCT(I,j)                    &
     &                 * (Z_UV(I,j,NTML(I,j)+1)-Z_UV(I,j,NTML(I,j))) )
          ENDIF
!           ! redesignate top-driven entrainment at ZHSC
!           ! (ignore that calculated at ZH)
            RHOKH_TOP_ENT(I,j) = RHOKH_DSCT_ENT(I,j)

            ZH(I,j) = ZHSC(I,j)
            NTML(I,j) = NTDSC(I,j)
            V_TOP(I,j) = V_TOP_DSC(I,j)
            ZSML_BASE(I,j) = 0.1 * ZH(I,j)
            ZHSC(I,j) = 0.0
            NTDSC(I,j) = 0
            V_TOP_DSC(I,j) = 0.0
            ZDSC_BASE(I,j) = 0.0
            DSC(I,j) = .FALSE.
            CUMULUS(I,j) = .FALSE.
            COUPLED(I,j) = .FALSE.
         ENDIF

        ELSE
!         !--------------------------------------------
!         ! Extent of top-driven mixing must be reduced
!         !--------------------------------------------
          RECOUPLED(I,j) = .FALSE.    ! recoupling not possible

!..Start iteration to find base of top-driven mixing, defined as the
!..height that gives WB_RATIO = DEC_THRES.
!..Procedure used makes 3 sweeps (up, down and back up again), using
!..progressively smaller increments (Z_INC), each time stopping when
!..the decoupling threshold or the height limits are passed.

!..Need to reset limits in case recoupling was attempted above
          SCBASE(I,j) = .FALSE. ! indicator of reaching height limits
          Z_TOP_LIM = Z_UV(I,j,NTDSC(I,j)-1)
          Z_BOT_LIM = 0.1 * ZH(I,j)
          NBDSC(I,j) = 2
!..Limit base of mixing to above ZH if cumulus-capped
          IF ( CUMULUS(I,j) ) THEN
            Z_BOT_LIM = Z_UV(I,j,NTML(I,j)+1)
            NBDSC(I,j) = NTML(I,j)+2
            IF ( Z_TOP_LIM  <   Z_BOT_LIM ) Z_TOP_LIM = Z_BOT_LIM
          ENDIF

!..Divide up depth of layer within which ZDSC_BASE is allowed
          Z_INC= ( Z_TOP_LIM - Z_BOT_LIM ) / FLOAT(N_STEPS-1)
          UP=1   ! indicates whether sweep is up or down

          ZDSC_BASE(I,j) = Z_BOT_LIM - Z_INC  ! to start at Z_BOT_LIM

          N_SWEEP = 1
          DO WHILE ( N_SWEEP  <=  3 .AND. .NOT. SCBASE(I,j) )

            NS = 1
            DO WHILE ( NS  <=  N_STEPS .AND. (                          &
     &                (UP  ==  1 .AND. WB_RATIO  >   DEC_THRES) .OR.    &
!                      ! work up while wb_ratio gt thres
     &                (UP  ==  0 .AND. WB_RATIO  <   DEC_THRES) ) )
!                      ! work down while wb_ratio lt thres

              WBP_INT = 1.0E-14
              WBN_INT = 0.0
              ZDSC_BASE(I,j) = ZDSC_BASE(I,j) + Z_INC
              ZH_PR = ZHSC(I,j) - ZDSC_BASE(I,j)

!..Calculate buoyancy flux profile given this ZDSC_BASE
              DO K = NBDSC(I,j), NTDSC(I,j)
                Z_PR = Z_UV(I,j,K) - ZDSC_BASE(I,j)
                IF (Z_PR  >   0.0 .AND. Z_PR  <   ZH_PR) THEN
                  RHOKH_SHAPE = (( 1.0 - Z_PR/ZH_PR )**0.8)             &
     &                                * Z_PR * Z_PR / ZH_PR
                  WB = RHOKH_SHAPE * DB_SCALED(I,j,K)
                ELSE
                  WB = 0.0
                ENDIF

                IF (WB  >=  0.0) THEN
                   WBP_INT = WBP_INT + WB
                ELSE
                   WBN_INT = WBN_INT - WB
                ENDIF
        ENDDO

              WB_RATIO = WBN_INT/WBP_INT
              NS = NS + 1

            ENDDO  ! loop stepping up through ML

!..sub-divide current Z_INC into one more part than there will be steps
!..as there is no need to recalculate WB at a current Z_INC
            Z_INC =  Z_INC / FLOAT(N_STEPS+1)

            IF ( (UP  ==  1 .AND. WB_RATIO  <=  DEC_THRES) .OR.         &
!                ! hit thres while working up
     &           (UP  ==  0 .AND. WB_RATIO  >=  DEC_THRES) ) THEN
!                ! hit thres while working down
              UP = 1-UP   ! change direction of sweep
              Z_INC = - Z_INC
            ELSEIF (ZDSC_BASE(I,j)  >=  Z_TOP_LIM - 1.0 .OR.            &
     &              ZDSC_BASE(I,j)  <=  Z_BOT_LIM + 1.0 ) THEN
!                ! hit height limits (give-or-take 1m) without
!                ! reaching threshold
              SCBASE(I,j) = .TRUE.
            ENDIF

!..Note that if the threshold has not been passed then the next sweep
!..continues in the same direction (but with reduced increment).

            N_SWEEP = N_SWEEP + 1
          ENDDO  ! loop over sweeps

        ENDIF  ! shallower mixing necessary

        ENDIF  ! NTDSC LE 2

      ENDDO  ! loop over I
      ENDDO  ! loop over j
! ----------------------------------------------------------------------
! 2.2 Second consider top-driven mixing in apparently surface-based MLs
!      = diagnose decoupling.  Start from the surface and work up.
!     Only do this if there isn't already a decoupled layer or Cu above;
!     if there is a cloud layer above, this is likely to insulate the
!     SML making any top-driven mixing only weak and, more importantly,
!     we can only deal with one decoupled layer at a time and so don't
!     want any more!
! ----------------------------------------------------------------------
      DO K=2,BL_LEVELS
        do j=1,rows
        do i=1,row_length
          DB_SCALED(I,j,K) = - 0.85*VKMAN*RHO_UV(I,j,K)*V_TOP(I,j)*     &
     &                                    DB_FOR_FLUX(I,j,K)
        ENDDO
        ENDDO
      ENDDO
        do j=1,rows
        do i=1,row_length
        IF ( NTML(I,j)  <=  2 .OR. CUMULUS(I,j) .OR. DSC(I,j) .OR.      &
     &         RECOUPLED(I,j) .OR. V_TOP(I,j)  <   0.1) THEN
!..no need to bother about decoupling or placing of top-driven mixing
          ZSML_BASE(I,j) = 0.1*ZH(I,j)
          ZSML_TOP(I,j) = ZH(I,j)
        ELSE
!..Start by testing wb integrals with mixing to the surface (or as close
!..as it is allowed to go without upsetting the surface layer
!..calculations).  As top-driven mixing is being considered here the
!..surface flux is ignored.
        ZSML_BASE(I,j) = 0.1*ZH(I,j)
        ZSML_TOP(I,j) = ZH(I,j)

        WBN_INT = 0.0
        WBP_INT = 1.0E-14
        ZH_PR = ZH(I,j) - ZSML_BASE(I,j)
        DO K=2,NTML(I,j)
         Z_PR = Z_UV(I,j,K) - ZSML_BASE(I,j)
         IF (Z_PR  >   0.0 .AND. Z_PR  <   ZH_PR) THEN
           RHOKH_SHAPE = (( 1.0 - Z_PR/ZH_PR )**0.8)                    &
     &                                * Z_PR * Z_PR / ZH_PR
           WB = RHOKH_SHAPE * DB_SCALED(I,j,K)
           ELSE
           WB = 0.0
         ENDIF

         IF (WB  >=  0.0) THEN
           WBP_INT = WBP_INT + WB
            ELSE
           WBN_INT = WBN_INT - WB
            ENDIF
        ENDDO

        WB_RATIO = WBN_INT/WBP_INT

        IF (WB_RATIO  >   DEC_THRES) THEN
!         !-----------------------------------------
!         ! Mixing down to the surface unsustainable
!         !-----------------------------------------

!..Start iteration to find base of top-driven mixing, defined as the
!..height that gives WB_RATIO = DEC_THRES.
!..Procedure used makes 3 sweeps (up, down and back up again), using
!..progressively smaller increments (Z_INC), each time stopping when
!..the decoupling threshold or the height limits are passed.
!..Limit depth of mixing to be at least two layers:
          SCBASE(I,j) = .FALSE. ! indicator of reaching height limits
          Z_TOP_LIM = Z_UV(I,j,NTML(I,j)-1)
          Z_BOT_LIM = 0.1 * ZH(I,j)

!..Divide up depth of layer within which ZSML_BASE is allowed
          Z_INC= ( Z_TOP_LIM - Z_BOT_LIM ) / FLOAT(N_STEPS)
          UP=1   ! indicates whether sweep is up or down

          N_SWEEP = 1
          DO WHILE ( N_SWEEP  <=  3 .AND. .NOT. SCBASE(I,j) )

            NS = 1
            DO WHILE ( NS  <=  N_STEPS .AND. (                          &
     &                (UP  ==  1 .AND. WB_RATIO  >   DEC_THRES) .OR.    &
!                      ! work up while wb_ratio gt thres
     &                (UP  ==  0 .AND. WB_RATIO  <   DEC_THRES) ) )
!                      ! work down while wb_ratio lt thres

              WBP_INT = 1.0E-14
              WBN_INT = 0.0
              ZSML_BASE(I,j) = ZSML_BASE(I,j) + Z_INC
              ZH_PR = ZH(I,j) - ZSML_BASE(I,j)

!..Calculate buoyancy flux profile given this ZSML_BASE
              DO K=2,NTML(I,j)
                Z_PR = Z_UV(I,j,K) - ZSML_BASE(I,j)
                IF (Z_PR  >   0.0 .AND. Z_PR  <   ZH_PR) THEN
                  KH_SMLT_SHAPE(I,j,K) = (( 1.0 - Z_PR/ZH_PR )**0.8)    &
     &                                * Z_PR * Z_PR / ZH_PR
                  WB = KH_SMLT_SHAPE(I,j,K) * DB_SCALED(I,j,K)
                ELSE
                  KH_SMLT_SHAPE(I,j,K) = 0.0
                  WB = 0.0
           ENDIF

                IF (WB  >=  0.0) THEN
                   WBP_INT = WBP_INT + WB
                ELSE
                   WBN_INT = WBN_INT - WB
          ENDIF
        ENDDO

              WB_RATIO = WBN_INT/WBP_INT

              NS = NS + 1

            ENDDO  ! loop stepping up through ML
!..sub-divide current Z_INC into one more part than there will be steps
!..as there is no need to WB for ZSML_BASE at a current Z_INC
            Z_INC =  Z_INC / FLOAT(N_STEPS+1)

            IF ( (UP  ==  1 .AND. WB_RATIO  <=  DEC_THRES) .OR.         &
!                ! hit thres while working up
     &           (UP  ==  0 .AND. WB_RATIO  >=  DEC_THRES) ) THEN
!                ! hit thres while working down
              UP = 1-UP   ! change direction of sweep
              Z_INC = - Z_INC
            ELSEIF (ZSML_BASE(I,j)  >=  Z_TOP_LIM - 1.0 .OR.            &
     &              ZSML_BASE(I,j)  <=  Z_BOT_LIM + 1.0 ) THEN
!                ! hit height limits (give-or-take 1m) without
!                ! reaching threshold
              SCBASE(I,j) = .TRUE.
            ENDIF

!..Note that if the threshold has not been passed then the next sweep
!..continues in the same direction (but with reduced increment).

            N_SWEEP = N_SWEEP + 1
          ENDDO  ! loop over sweeps

          IF ( SCBASE(I,j) ) THEN
!..hit height limits so mix up to ZH, in either case
            ZSML_TOP(I,j) = ZH(I,j)
          ELSE
!-----------------------------------------------------------------------
!..Decoupling apparent.  Repeat the above procedure to see where the top
!..of the surface-driven profile (at ZSML_TOP) should be (by requiring
!..the total buoyancy flux to satisfy the integral conditions).
!-----------------------------------------------------------------------
!..Start with ZSML_TOP at ZSML_BASE
            ZSML_TOP(I,j) = ZSML_BASE(I,j)

            WB = BFLUX_SURF(I,j) / RDZ(I,j,1)
            IF (WB  >   0.0) THEN
              WBP_INT = WB
              WBN_INT = 0.0
            ELSE
              WBP_INT = 1.E-14
              WBN_INT = - WB
          ENDIF

            ZH_PR = ZSML_TOP(I,j)
!..Calculate buoyancy flux profile given this ZSML_TOP
            DO K=2,NTML(I,j)
              Z_PR = Z_UV(I,j,K)
              IF (Z_PR  <   ZSML_TOP(I,j)) THEN
!               !  recall DB_SCALED is scaled on V_TOP
                RHOKH_SHAPE = ( V_SURF(I,j) / V_TOP(I,j) )              &
     &                          * Z_PR * ( 1.0 - Z_PR/ZH_PR )           &
     &                                 * ( 1.0 - Z_PR/ZH_PR )
              ELSE
                RHOKH_SHAPE = 0.0
         ENDIF
              WB = ( RHOKH_SHAPE + KH_SMLT_SHAPE(I,j,K) ) *             &
     &                                              DB_SCALED(I,j,K)

              IF (WB  >=  0.0) THEN
                WBP_INT = WBP_INT + WB
              ELSE
                WBN_INT = WBN_INT - WB
              ENDIF
        ENDDO

            WB_RATIO = WBN_INT/WBP_INT

            IF (WB_RATIO  <   DEC_THRES) THEN
!..WB_RATIO indicates surface profile can be extended higher
              SCBASE(I,j) = .FALSE. ! indicator of reaching height limit
              Z_TOP_LIM = ZH(I,j)
              Z_BOT_LIM = ZSML_BASE(I,j)

!..Divide up depth of layer within which ZSML_BASE is allowed
              Z_INC= ( Z_TOP_LIM - Z_BOT_LIM ) / FLOAT(N_STEPS)
              UP=1   ! indicates whether sweep is up or down

              N_SWEEP = 1
              DO WHILE ( N_SWEEP  <=  3 .AND. .NOT. SCBASE(I,j) )

                NS = 1
                DO WHILE ( NS  <=  N_STEPS .AND. (                      &
     &                (UP  ==  1 .AND. WB_RATIO  <   DEC_THRES) .OR.    &
!                      ! work up while wb_ratio lt thres
     &                (UP  ==  0 .AND. WB_RATIO  >   DEC_THRES) )       &
!                      ! work down while wb_ratio gt thres
     &                  .AND. .NOT.SCBASE(I,j) )
!                      ! stop if
                  WB = BFLUX_SURF(I,j) / RDZ(I,j,1)
                  IF (WB  >   0.0) THEN
                    WBP_INT = WB
                    WBN_INT = 0.0
                  ELSE
                    WBN_INT = - WB
                    WBP_INT = 1.E-14
                  ENDIF
                  ZSML_TOP(I,j) = ZSML_TOP(I,j) + Z_INC
                  ZH_PR = ZSML_TOP(I,j)

!..Calculate buoyancy flux profile given this ZSML_TOP
                  DO K=2,NTML(I,j)
                    Z_PR = Z_UV(I,j,K)
                    IF (Z_PR  <   ZSML_TOP(I,j)) THEN
                      RHOKH_SHAPE = ( V_SURF(I,j) / V_TOP(I,j) )        &
     &                          * Z_PR * ( 1.0 - Z_PR/ZH_PR )           &
     &                                 * ( 1.0 - Z_PR/ZH_PR )
!..recall DB_SCALED is scaled on V_TOP
                    ELSE
                      RHOKH_SHAPE = 0.0
                    ENDIF
                    WB = ( RHOKH_SHAPE+KH_SMLT_SHAPE(I,j,K) )           &
     &                                         * DB_SCALED(I,j,K)

                    IF (WB  >=  0.0) THEN
                       WBP_INT = WBP_INT + WB
                    ELSE
                       WBN_INT = WBN_INT - WB
                    ENDIF
        ENDDO

                  WB_RATIO = WBN_INT/WBP_INT
                  NS = NS + 1

                ENDDO  ! loop stepping up through ML

!..sub-divide current Z_INC into one more part than there will be steps
!..as there is no need to WB for ZSML_TOP at a current Z_INC
                Z_INC =  Z_INC / FLOAT(N_STEPS+1)

                IF ( (UP  ==  1 .AND. WB_RATIO  >=  DEC_THRES) .OR.     &
!                ! hit thres while working up
     &             (UP  ==  0 .AND. WB_RATIO  <=  DEC_THRES) ) THEN
!                ! hit thres while working down
                  UP = 1-UP   ! change direction of sweep
                  Z_INC = - Z_INC
                ELSEIF (ZSML_TOP(I,j)  >=  Z_TOP_LIM - 1.0 .OR.         &
     &                  ZSML_TOP(I,j)  <=  Z_BOT_LIM + 1.0) THEN
!                ! hit height limits (give-or-take 1m) without
!                ! reaching threshold
                  SCBASE(I,j) = .TRUE.
                ENDIF

!..Note that if the threshold has not been passed then the next sweep
!..continues in the same direction (but with reduced increment).

                N_SWEEP = N_SWEEP + 1

              ENDDO  ! loop over sweeps

            ENDIF  ! higher K_surf than K_top base possible

          ENDIF  ! test on whether to limit KH_SURF (test on SCBASE)
!---------------------------------------------------------------------
!..Have now calculated:
!       ZSML_TOP = top of surface-driven mixing
!       ZSML_BASE = base of top-driven mixing
!..If ZSML_TOP is less than ZH, put all entrainment into RHOKH_DSCT.
!..Call this a decoupled layer (as mixing does not extend over the
!..whole layer) and set flags accordingly.

          IF (ZSML_TOP(I,j)  <   ZH(I,j) ) THEN
            RHOKH_DSCT_ENT(I,j) = RHOKH_TOP_ENT(I,j)                    &
     &                             + RHOKH_SURF_ENT(I,j)
            RHOKH_TOP_ENT(I,j) = 0.0
            RHOKH_SURF_ENT(I,j) = 0.0
            RHOKM_TOP(I,j,NTML(I,j)+1) = RHOKM_TOP(I,j,NTML(I,j)+1)     &
     &                                 + RHOKM(I,j,NTML(I,j)+1)
            RHOKM(I,j,NTML(I,j)+1) = 0.0

            DSC(I,j) = .TRUE.
            NTDSC(I,j) = NTML(I,j)
            ZHSC(I,j) = ZH(I,j)
            ZDSC_BASE(I,j) = ZSML_BASE(I,j)
            V_TOP_DSC(I,j) = V_TOP(I,j)

            KTOP = 2
            DO WHILE ( Z_UV(I,j,KTOP+1)  <   ZSML_TOP(I,j) )
              KTOP = KTOP + 1
        ENDDO

            ZH(I,j) = ZSML_TOP(I,j)
            NTML(I,j) = KTOP - 1
!..switch off top-down mixing in SML because necessary V_TOP not known
            V_TOP(I,j) = 0.0
            ZSML_BASE(I,j) = ZH(I,j)

            COUPLED(I,j) = .TRUE.
!          ! because full entrainment applied at ZHSC
          ENDIF

        ENDIF  ! mixing down to the surface unsustainable

        ENDIF  ! NTML LE 2, CUMULUS, V_TOP < 0.1

      ENDDO  ! loop over I
      ENDDO  ! loop over j

!-----------------------------------------------------------------------
!     Calculate factors required to ensure that the non-local turbulent
!     mixing coefficient profiles are continuous as the entrainment
!     level is approached.
!-----------------------------------------------------------------------
!
      do j=1,rows
      do i=1,row_length
          K=NTML(I,j)+1
!       ! for cubic form of KH and KM:
            KH_TOP_FACTOR(I,j) = MAX( 0.7 , 1.0 - SQRT(                 &
     &           RHOKH_SURF_ENT(I,j) /                                  &
     &                 ( RHO_UV(I,j,K)*W_H_TOP(I,j)*VKMAN*ZH(I,j) ) ) )
            KM_TOP_FACTOR(I,j) = MAX( 0.7 , 1.0 - SQRT( RHOKM(I,j,K) /  &
     &             ( RHO_TQ(I,j,K-1)*W_M_TOP(I,j)*VKMAN*ZH(I,j) ) ) )
!
!       ! for quadratic form of KH and KM:
!           KH_TOP_FACTOR(I,j) = MAX( 0.9 , 1.0 - RHOKH_SURF_ENT(I,j) /
!    &                   ( RHO_UV(I,j,K)*W_H_TOP(I,j)*VKMAN*ZH(I,j) ) )
!           KM_TOP_FACTOR(I,j) = MAX( 0.9 , 1.0 - RHOKM(I,j,K) /
!    &              ( RHO_TQ(I,j,K-1)*W_M_TOP(I,j)*VKMAN*ZH(I,j) ) )
!
        SCDEPTH(I,j) = ZH(I,j) - ZSML_BASE(I,j)
        FACTOR = 0.85 * RHO_UV(I,j,K) * V_TOP(I,j) *VKMAN *SCDEPTH(I,j)
          IF ( FACTOR  >   0.0) THEN
            KH_SCT_FACTOR(I,j) = 1.0 -                                  &
     &     ( (RHOKH_TOP_ENT(I,j)*RHOKH_TOP_ENT(I,j)) / (FACTOR*FACTOR) )
          ELSE
            KH_SCT_FACTOR(I,j) = 1.0
          ENDIF
        FACTOR = 0.85 * RHO_TQ(I,j,K-1) * V_TOP(I,j) *                  &
     &                  VKMAN * SCDEPTH(I,j) * 0.75
          IF ( FACTOR  >   0.0) THEN
            KM_SCT_FACTOR(I,j) = 1.0 -                                  &
     &         ( (RHOKM_TOP(I,j,K)*RHOKM_TOP(I,j,K)) / (FACTOR*FACTOR) )
          ELSE
            KM_SCT_FACTOR(I,j) = 1.0
          ENDIF
!
!       !---------------------------------------------------------------
!       !  Set up factors to ensure K profile continuity at ZHSC;
!       !  no need to limit size of factor as precise shape of top-down
!       !  mixing profile not important.
!       !---------------------------------------------------------------
!
       IF (NTDSC(I,j)  >   0) THEN
!       !-------------------------------------------------------------
!       ! Only calculate _DSCT_FACTORs when a decoupled stratocumulus
!       ! layer exists, i.e. NTDSC > 0.
!       !-------------------------------------------------------------
          K=NTDSC(I,j)+1
        DSCDEPTH(I,j) = ZHSC(I,j) - ZDSC_BASE(I,j)
        FACTOR = 0.85*RHO_UV(I,j,K)*V_TOP_DSC(I,j)*VKMAN*DSCDEPTH(I,j)
          IF ( FACTOR  >   0.0) THEN
            KH_DSCT_FACTOR(I,j) = 1.0 -                                 &
     &        ( (RHOKH_DSCT_ENT(I,j)*RHOKH_DSCT_ENT(I,j)) /             &
     &                                                 (FACTOR*FACTOR) )
          ELSE
            KH_DSCT_FACTOR(I,j) = 1.0
          ENDIF

        FACTOR = 0.75 * 0.85 * RHO_TQ(I,j,K-1) * V_TOP_DSC(I,j) *       &
     &                         VKMAN * DSCDEPTH(I,j)
          IF ( FACTOR  >   0.0) THEN
            KM_DSCT_FACTOR(I,j) = 1.0 -                                 &
     &         ( (RHOKM_TOP(I,j,K)*RHOKM_TOP(I,j,K)) / (FACTOR*FACTOR) )
          ELSE
            KM_DSCT_FACTOR(I,j) = 1.0
          ENDIF
       ENDIF
        ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!! 2.  Loop around levels again calculating height dependent turbulent
!!     transport coefficients within the mixing layer.
!-----------------------------------------------------------------------
!
! Reset identifiers of base of decoupled layer mixing
!
      do j=1,rows
      do i=1,row_length
        SCBASE(I,j) = .FALSE.
        NBDSC(I,j)  = 0
      ENDDO
      ENDDO

      DO K=2,BL_LEVELS
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
          ENDIF
!
! Keep depth of mixing fixed for momentum
!
          ZCML_BASE = MAX( 0.1*Z_TQ(I,j,NTML(I,j)) ,                    &
     &                                Z_TQ(I,j,NTML(I,j))-SCDEPTH(I,j) )
          IF ( K  <=  NTML(I,j) .AND.                                   &
     &            ZK_TQ  >   ZCML_BASE ) THEN
              Z_PR = ZK_TQ - ZCML_BASE
              ZH_PR = Z_TQ(I,j,NTML(I,j)) - ZCML_BASE
              RHOKM_TOP(I,j,K) = 0.75 * RHO_TQ(I,j,K-1) * V_TOP(I,j) *  &
     &              0.85 * VKMAN *                                      &
     &              ( ( 1.0 - KM_SCT_FACTOR(I,j)*Z_PR/ZH_PR )**0.8 )    &
     &                                         * Z_PR * Z_PR / ZH_PR
!                                                     ! PRANDTL=0.75
          ENDIF
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
            ENDIF
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
          ENDIF
!         !-------------------------------------------------------------
!         ! Keep depth of mixing fixed for momentum.
!         !-------------------------------------------------------------
          IF ( K  <=  NTDSC(I,j) ) THEN
          ZCML_BASE = MAX( 0.1*Z_TQ(I,j,NTML(I,j)) ,                    &
     &                              Z_TQ(I,j,NTDSC(I,j))-DSCDEPTH(I,j) )
           IF ( ZK_TQ  >   ZCML_BASE ) THEN
              Z_PR = ZK_TQ - ZCML_BASE
              ZH_PR = Z_TQ(I,j,NTDSC(I,j)) - ZCML_BASE
              RHOKM_TOP(I,j,K) = RHOKM_TOP(I,j,K) +                     &
     &           0.75*RHO_TQ(I,j,K-1)*V_TOP_DSC(I,j)*0.85*VKMAN*        &
     &              ( ( 1.0 - KM_DSCT_FACTOR(I,j)*Z_PR/ZH_PR )**0.8 )   &
     &                                      * Z_PR * Z_PR / ZH_PR
          ENDIF
          ENDIF

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
              ENDIF
            ENDIF

            IF (ZK_TQ  <=  0.1*Z_TQ(I,j,NTML(I,j))) THEN
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
            ENDIF
            ENDIF
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
!             N.B. ZH(I,j) = Z_UV(I,j,NTML(I,j)+1)
!
!               with cubic form
!
              RHOKH(I,j,K) = RHO_UV(I,j,K) * W_H_UV * VKMAN * ZK_UV *   &
     &            ( 1.0 - KH_TOP_FACTOR(I,j) * ( ZK_UV / ZH(I,j) ) ) *  &
     &            ( 1.0 - KH_TOP_FACTOR(I,j) * ( ZK_UV / ZH(I,j) ) )
!
!               with quadratic form
!
!             RHOKH(I,j,K) = RHO_UV(I,j,K) * W_H_UV * VKMAN * ZK_UV *
!    &                ( 1.0 - KH_TOP_FACTOR(I,j) * ( ZK_UV / ZH(I,j) ) )
            ENDIF
            IF ( K  <=  NTML(I,j) ) THEN

!             !---------------------------------------------------------
!             ! Calculate RHOKM(w_m,z/z_h)
!             !---------------------------------------------------------
              ZH_M = Z_TQ(I,j,NTML(I,j))
!
!               with cubic form
!
              RHOKM(I,j,K) = RHO_TQ(I,j,K-1) * W_M_TQ * VKMAN * ZK_TQ * &
     &                ( 1.0 - KM_TOP_FACTOR(I,j) * ( ZK_TQ / ZH_M ) ) * &
     &                ( 1.0 - KM_TOP_FACTOR(I,j) * ( ZK_TQ / ZH_M ) )
!
!               with quadratic form
!
!             RHOKM(I,j,K) = RHO_TQ(I,j,K-1) * W_M_TQ * VKMAN * ZK_TQ *
!    &                ( 1.0 - KM_TOP_FACTOR(I,j) * ( ZK_TQ / ZH_M ) )

            ENDIF
          ENDIF
        ENDDO
        ENDDO
      ENDDO

      IF (NG_STRESS  ==  1) THEN
        DO K=2,BL_LEVELS
          do j=1,rows
          do i=1,row_length
            IF (FB_SURF(I,j)  >   0.0 .AND. K  <=  NTML(I,j) ) THEN
!             !---------------------------------------------------------
!             ! Calculate non-gradient stress function
!             ! (Brown and Grant 1997)
!             ! Shape function chosen such that non-gradient stress
!             ! goes to zero at 0.1*ZH and ZH
!             !---------------------------------------------------------
              ZK_TQ = Z_TQ(I,j,K-1)   ! stresses are calc on theta-levs
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
                ENDIF



! 4*W_S_CUBED_TQ = the convective boundary layer
! velocity scale cubed
! V_S = the friction velocity

                F_NGSTRESS(I,j,K) = ( RHO_TQ(I,j,K-1)/RHOSTAR_GB(I,j) ) &
     &            * S_M * ( A_NGS * 4.0*W_S_CUBED_TQ/                   &
     &             (V_S(I,j)*V_S(I,j)*V_S(I,j) + W_S_CUBED_TQ*4.*0.6 ) )&
     &               * ( Z_PR / ZH_PR ) * ( 1.0 -  ( Z_PR / ZH_PR ) ) * &
     &                                    ( 1.0 -  ( Z_PR / ZH_PR ) )


              ENDIF
            ENDIF
          ENDDO
          ENDDO
        ENDDO
      ENDIF

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('EXCF_NL  ',4)
      ENDIF
      RETURN
      END SUBROUTINE EXCF_NL
#endif
