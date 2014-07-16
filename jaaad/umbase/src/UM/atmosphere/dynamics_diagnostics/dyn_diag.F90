#if defined(A15_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculate various diagnostics related to dynamics variables.
!
! Subroutine Interface:
      SUBROUTINE Dyn_diag(                                              &
! Primary data: in
     & exner_rho_levels,rho,u,v,w                                       &
     &,exner_theta_levels                                               &
     &,theta                                                            &
! Grid sizes and definition: in
     &,rows,n_rows,row_length,model_levels,wet_levels,bl_levels         &
     &,global_rows,global_row_length                                    &
     &,theta_field_size,u_field_size,v_field_size                       &
     &,eta_theta_levels,eta_rho_levels                                  &
     &,Model_domain                                                     &
! Grid coordinates: in
     &,delta_lambda,delta_phi                                           &
     &,ew_space,ns_space,first_lat,first_long,phi_pole,lambda_pole      &
! Pre-calculated grid associated arrays: in
     &,r_at_u,r_at_v                                                    &
     &, r_theta_levels, r_rho_levels, sec_v_latitude                    &
     &, tan_v_latitude, sec_theta_latitude, f3_at_v                     &
     &,rot_coeff1,rot_coeff2                                            &
! Time information: in
     &,forecast_hrs                                                     &
! Theta levels for output arrays
     &, desired_theta                                                   &
! Pressure levels for output arrays: in
     &,ucomB_press,vcomB_press                                          &
     &,ucomp_press,vcomp_press,wcomp_press                              &
     &,testd_press                                                      &
     &,p_height,theta_height,rho_height                                 &
     &,w_height,u_height,v_height                                       &
     &,pv_press                                                         &
! Model levels for output arrays: in
     &,ucomB_model,vcomB_model                                          &
     &,testd_model                                                      &
     &,Htheta_model,Hrho_model                                          &
! Flags to request each diagnostic output field: in
! wind related diagnostics
     &,qucomB_m,qvcomB_m                                                &
     &,qucomB_p,qvcomB_p                                                &
     &,qucomp_p,qvcomp_p,qwcomp_p                                       &
     &,qu50mB_h,qv50mB_h                                                &
     &,qu50m_h,qv50m_h                                                  &
! PV related diagnostics
     &,qpotn_vort_theta,qtheta_potn_vort,qtheta_pv_points               &
     &,qpv_mod_levs,qpv_theta_mlev,qpotn_vort_press                     &
! test fields
     &,qdia1,qdia2,qdia3,qdia4                                          &
! flux diagnostics
     &,qrhow,qrhouw,qrhovw,qrhow_up                                     &
     &,qrhow_down,qrhowc_up,qrhowc_down                                 &
! height and height level diagnostics
     &,qHtheta_ml,qHrho_ml                                              &
     &,qpress_h,qtheta_h,qrho_h                                         &
     &,qu_h,qv_h,qw_h                                                   &
! other diagnostics
     &,spec_w,qtrue_density                                             &
! Flags for wind rotation (lam grid): in
     &,rot_uvcomB_p                                                     &
! Diagnostics lengths: in
     &,ucomB_m_levs,vcomB_m_levs                                        &
     &,ucomB_p_levs,vcomB_p_levs                                        &
     &,ucomp_p_levs,vcomp_p_levs,wcomp_p_levs                           &
     &,pv_theta_levs,pv_press_levs                                      &
     &,testd_p_levs,testd_m_levs                                        &
     &,Htheta_m_levs,Hrho_m_levs                                        &
     &,p_h_levs,theta_h_levs,rho_h_levs,w_h_levs,u_h_levs,v_h_levs      &
! Diagnostic arrays: out
! wind related diagnostics
     &,ucomB_m,vcomB_m                                                  &
     &,ucomB_p,vcomB_p                                                  &
     &,ucomp_p,vcomp_p,wcomp_p                                          &
     &,u50mB_h,v50mB_h                                                  &
     &,u50m_h,v50m_h                                                    &
! PV related diagnostics
     &,potn_vort_theta,theta_potn_vort,theta_pv_points                  &
     &,pv_mod_levs,pv_theta_mlev,potn_vort_press                        &
! test fields
     &,testdiag1,testdiag2,testdiag3,testdiag4                          &
! flux diagnostics
     &,rhow,rhouw,rhovw,rhow_up                                         &
     &,rhow_down,rhow_convup,rhow_convdown                              &
! height and height level diagnostics
     &,height_theta_ml,height_rho_ml                                    &
     &,press_h,theta_h,rho_h                                            &
     &,ucomp_h,vcomp_h,wcomp_h                                          &
! other diagnostics
     &,spec_3D,true_density)

#if defined(FLUME)
! FLUME-STASH 
USE MATMFlumeModel, only:flumeSendDiag
USE SharedID
USE flumerun
#endif

      IMPLICIT NONE
!
! Description:
!   Calculate dynamics-related diagnostics - held in STASH section 15 -
!   which may include interpolation onto pressure surfaces. Diagnostics
!   currently supported:
!   [All diagnostics on native 'C' grid unless specified.]
!   STASH item
!     2 u component of wind on model levels      on 'B' grid
!     3 v component of wind on model levels      on 'B' grid
!   201 u component of wind on pressure surfaces on 'B' grid
!   202 v component of wind on pressure surfaces on 'B' grid
!   243 u component of wind on pressure surfaces
!   244 v component of wind on pressure surfaces
!   242 w component of wind on pressure surfaces
!   245 u component of wind at 50m height
!   246 v component of wind at 50m height
!   212 u component of wind at 50m height on 'B' grid
!   213 v component of wind at 50m height on 'B' grid
!   214 potential vorticity on theta levels
!   229 potential vorticity on pressure levels
!   215 theta on potential vorticity = +/-2 surface
!   216 theta at potential vorticity points
!   217 potential vorticity on model levels
!   218 potential vorticity on model theta grid and theta levels
!   231 test analytic field on v grid - single level
!   232 test analytic field on p grid - single level
!   233 test analytic field on p grid - pressure levels
!   234 test analytic field on p grid - model levels
!   260 mass flux (rhow) on model levels
!   261 momentum flux (rhouw) on model levels
!   262 momentum flux (rhovw) on model levels
!   263 upward mass flux (rhow, w> 0m/s) on model levels
!   264 downward mass flux (rhow, w<0m/s) on model levels
!   265 upward convective mass flux (rhow, w>1m/s) on model levels
!   266 downward convective mass flux (rhow, w<-1m/s) on model levels
!   101 Height from sea level on theta model levels
!   102 Height from sea level on rho model levels
!   108 Pressure on height (from sea) levels
!   119 Potential temperature on height (from sea) levels
!   127 Rho (density) on height (from sea) levels
!   142 W component of wind on height (from sea) levels
!   143 U component of wind on height (from sea) levels
!   144 V component of wind on height (from sea) levels
!   270 real part of spectra (model levels)
!   271 true unscaled density on model levels
!
! Method:
!   Required level lists and logical switches are determined by the
!   calling routine from STASH requests and STASHflags.
!   Primary model data is input, and each diagnostic is calculated - if
!   its flag is set - in simple sequential order. Where the extraction
!   of the diagnostic quantity requires further calculation, a lower
!   level diagnostic-specific routine is called.
!
! STASH items 002,003     : u,v wind components on model levs 'B' grid:
!    Simple horizontal interpolation from u,v to uv 'B' position.
!
! STASH items 201,202     : u,v wind components on p surfaces 'B' grid:
! 1. Simple horizontal interpolation of p onto u,v staggering.
! 2. Perform vertical interpolation of u,v onto each o/p p surface.
! 3. Simple horizontal interpolation from u,v to uv 'B' position.
! 4. (Lam only) Rotate winds from native eq lat-long to standard grid
!
! STASH item 214: pv on theta levels:
! 1. Interpolate theta onto PV points.
! 2. Calculate PV on model levels.
! 3. Perform vertical interpolation of PV onto theta levels.

! STASH items 243,244,242 : u,v,w wind components on pressure surfaces:
! 1. Simple horizontal interpolation of p onto u,v staggering.
! 2. Perform vertical interpolation of u,v onto each o/p p surface.
!
! STASH items 212,213 : u,v wind components at 50m height 'B' grid:
! STASH items 244,245 : u,v wind components at 50m height:
! 1. Re-use p_at_u/v arrays to get height surface 50m above orography
! 2. Perform vertical interpolation of u,v onto height surface
! 3. Simple horizontal interpolation from u,v to uv 'B' position.
!
! STASH items 231,232,233,234: Test diagnostics 1-4:
! Call TestDiag routine. See UM Doc Paper D7.
!
! STASH item 215: theta on pv=+/-2 surface:
! 1. Interpolate theta onto PV points.
! 2. Calculate PV on model levels.
! 3. Set mod_PV = |PV|
! 4. Perform vertical interpolation of theta onto PV=2 surface
!
! STASH item 216: theta at pv points:
! 1. Interpolate theta onto PV points.
!
! STASH item 217: pv on model levels:
! 1. Calculate PV on model levels.
!
! STASH item 218: pv on model theta levels:
! 1. Calculate PV on model levels, using 'calc_pv_at_theta'.
!
! Current Code Owner: <Rick Rawlins>
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.0  02/06/99   Extensive revision of DYNDIA1A deck for 'C-P C
!                  dynamics upgrade' project. Rick Rawlins.
!  5.1  25/01/00   Change u,v diagnostics from 'C' to 'B' grid. Retain
!                  old diagnostics but with new STASHcodes.
!                  Add u,v on model levels on 'B' grid.
!                  Replace rmdi_pp by rmdi. R Rawlins
!  5.2  14/09/00   Add PV on theta levels.  Z. Gardner
!  5.2  19/03/01   Change to use exner pressures instead of p in
!                  vertical interpolation and change
!                  call to vert_interp_mdi2 to vert_interp2
!  5.2  31/07/00   (1) Correction for (15,002-003): add
!                  swapbounds to update halos.
!                  (2) Introduce rotation of winds for a subset of
!                  diagnostics in lam. R Rawlins
!  5.2  15/11/00   Set up stash variables levels required for pv on
!                  theta levels diagnostic.  Zoe Gardner
!  5.3  27/04/01   Correct non-reprod u&v on pressure levels
!                  caused by acw4f502 changes.   S. Cusack
!  5.3  06/09/01   Set up stash variables levels required for pv on
!                  pressure levels diagnostic.  D.M. Goddard
!  5.4  28/05/02   Add theta on pv=2 or -2 surface also theta at pv poin
!                  and pv on model levels. T.J. Hinton
!  5.4  16/07/02   Introduce new diagnostics: u,v,w,theta,rho,p on
!                  geometric height levels and height (in meters) of
!                  theta, rho model levels from sea level.
!                               M. Diamantakis
!  5.4     28/08/02    Bug Fix (Bi-cyclic LAM)           Carol Roadnight
!  5.4  15/05/02   Correction to get bit comparison of potential
!                  vorticity diagnostics over different processor
!                  configurations. D.M. Goddard
!  5.4  11/03/02   Remove comment on same line as #include
!                                               S. Carroll
!  5.5  28/02/03   Add new mass flux diagnostics       Carol Roadnight
!  5.5  26/02/03   Get averaged spectra of vertical velocity field
!                                                       Carol Roadnight
!  6.0  18/07/03   Merge theta on PV = +/- 2 fields T.J.Hinton
!  6.1  17/05/04   Calc_PV arguments corrected. Adam Clayton
!  6.2  15/08/05   Free format fixes. P.Selwood
!  6.4  16/11/06   Get out true unscaled density          Andy Malcolm
!  6.4  16/11/06   Remove unnecessary swap_bounds calls   Andy Malcolm
!  6.4  16/11/06   Move & combine swap_bounds calls logic Andy Malcolm
!  6.4  16/11/06   Change code so only 1 call to CALC_PV  Andy Malcolm
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Global variables (*CALLed COMDECKs etc...):
#include "parvars.h"
#include "cphyscon.h"
#include "interpor.h"
! Model_domain meaningful names
#include "domtyp.h"

! Subroutine arguments
!   Scalar arguments with intent(in):
! Grid sizes:
      INTEGER                                                           &
     & rows,n_rows,row_length,model_levels,wet_levels,bl_levels         &
     &, global_rows,global_row_length                                   &
     &,theta_field_size,u_field_size,v_field_size                       &
     &,Model_domain
! Grid coordinates: in
      REAL                                                              &
     & delta_lambda,delta_phi                                           &
     &,ew_space,ns_space,first_lat,first_long,phi_pole,lambda_pole
! Time information: in
      REAL                                                              &
     & forecast_hrs       ! T+forecast_hrs: hours after analysis time
                          ! UM6.5 - MODEL_ANALYSIS_HRS changed to REAL -
                          ! requires FORECAST_HRS changed to REAL also  
! Diagnostics lengths: in
      INTEGER                                                           &
     & ucomB_m_levs                                                     & 
                          ! No of levels for output u on m 'B' grid
     &,vcomB_m_levs                                                     & 
                          ! No of levels for output v on m 'B' grid
     &,ucomB_p_levs                                                     & 
                          ! No of levels for output u on p 'B' grid
     &,vcomB_p_levs                                                     & 
                          ! No of levels for output v on p 'B' grid
     &,ucomp_p_levs                                                     &
                          ! No of levels for output of u on p
     &,vcomp_p_levs                                                     &
                          ! No of levels for output of v on p
     &,wcomp_p_levs                                                     &
                          ! No of levels for output of w on p
     &, pv_theta_levs                                                   &
                         !No of levels for output of pv on theta
     &, pv_press_levs                                                   &
                            !No of levels for output of pv on pressure
     &,testd_p_levs                                                     &
                          ! No of levels for output of testdiag3
     &,testd_m_levs                                                     &
                          ! No of levels for output of testdiag4
     &,Htheta_m_levs                                                    &
                              ! Num of levs for output H on theta lev
     &,Hrho_m_levs                                                      &
                              ! Num of levs for output H on rho lev
     &,p_h_levs                                                         &
                              ! Num of levs for output p on H levs
     &,theta_h_levs                                                     &
                              ! Num of levs for output theta on H levs
     &,rho_h_levs                                                       &
                              ! Num of levs for output rho on H levs
     &,w_h_levs                                                         &
                              ! Num of levs for output w on H levs
     &,u_h_levs                                                         &
                              ! Num of levs for output u on H levs
     &,v_h_levs               ! Num of levs for output v on H levs
! Flags to request each diagnostic output field: IN
      LOGICAL                                                           &
     & qucomB_m                                                         & 
                     ! Flag for U wind on model    levels  'B' grid
     &,qvcomB_m                                                         & 
                     ! Flag for V wind on model    levels  'B' grid
     &,qucomB_p                                                         & 
                     ! Flag for U wind on pressure levels  'B' grid
     &,qvcomB_p                                                         & 
                     ! Flag for V wind on pressure levels  'B' grid
     &,qucomp_p                                                         &
                     ! Flag for U wind component on pressure levels
     &,qvcomp_p                                                         &
                     ! Flag for V wind component on pressure levels
     &,qwcomp_p                                                         &
                     ! Flag for W wind component on pressure levels
     &,qu50mB_h                                                         & 
                     ! Flag for U wind comp at 50m height  'B' grid
     &,qv50mB_h                                                         & 
                     ! Flag for V wind comp at 50m height  'B' grid
     &,qu50m_h                                                          &
                     ! Flag for U wind component at 50m height
     &,qv50m_h                                                          &
                     ! Flag for V wind component at 50m height
     &, qpotn_vort_theta                                                &
                            !Flag for pv on theta levels
     &, qpotn_vort_press                                                &
                            !Flag for pv on pressure levels
     &,qtheta_potn_vort                                                 &
                           !Flag for theta on pv=+/-2 surface
     &,qtheta_pv_points                                                 &
                           !Flag for theta at pv points
     &,qpv_mod_levs                                                     &
                           !Flag for pv on model levels
     &,qpv_theta_mlev                                                   &
                           !Flag for pv on model theta points and levels
     &,qdia1                                                            &
                     ! Flag for test diagnostic 1
     &,qdia2                                                            &
                     ! Flag for test diagnostic 2
     &,qdia3                                                            &
                     ! Flag for test diagnostic 3
     &,qdia4                                                            &
                     ! Flag for test diagnostic 4
     &,qrhow                                                            &
                     ! Flag for rhow
     &,qrhouw                                                           &
                     ! Flag for rhouw
     &,qrhovw                                                           &
                     ! Flag for rhovw
     &,qrhow_up                                                         &
                     ! Flag for rhow_up
     &,qrhow_down                                                       &
                     ! Flag for rhow_down
     &,qrhowc_up                                                        &
                     ! Flag for rhow_convup
     &,qrhowc_down                                                      &
                     ! Flag for rhow_convdown
     &,qHtheta_ml                                                       &
                                 ! Flag for height on theta levs
     &,qHrho_ml                                                         &
                                 ! Flag for height on rho levs
     &,qpress_H                                                         &
                               ! Flag for press on height levs
     &,qtheta_H                                                         &
                               ! Flag for theta on height levs
     &,qrho_H                                                           &
                               ! Flag for rho on height levs
     &,qw_H                                                             &
                               ! Flag for w on height levs
     &,qu_H                                                             &
                               ! Flag for u on height levs
     &,qv_H                                                             &
                               ! Flag for v on height levs
     &,spec_w                                                           &
                               ! Flag for real part of spectra          
     &,qtrue_density           ! flag for true unscaled density
! Flags for wind rotation (lam grid)  (rotate if .T. ) : in
      LOGICAL, INTENT(IN) ::                                            &
     & rot_uvcomB_p          ! u,v (B grid) on p levels
!   Array  arguments with intent(in):
! Primary data: IN
      REAL                                                              &
     & u  (1-offx:row_length+offx, 1-offy:rows+offy,    model_levels)   &
     &,v  (1-offx:row_length+offx, 1-offy:n_rows+offy,  model_levels)   &
     &,w  (1-offx:row_length+offx, 1-offy:rows+offy,  0:model_levels)   &
     &,rho(1-offx:row_length+offx, 1-offy:rows+offy,    model_levels)   &
     &,exner_rho_levels                                                 &
     &    (1-offx:row_length+offx, 1-offy:rows+offy, model_levels)      &
     &,exner_theta_levels                                               &
     &    (1-offx:row_length+offx, 1-offy:rows+offy,    model_levels)   &
     &, theta(1-offx:row_length+offx, 1-offy:rows+offy,model_levels)    &
! Vertical grid definition: IN
     &,eta_theta_levels(0:model_levels)                                 &
                                        ! vertical grid for theta vars
     &,eta_rho_levels    (model_levels)                                 &
                                        ! vertical grid for rho   vars
! Pre-calculated grid associated arrays: IN
     &,r_at_u (1-halo_i:row_length+halo_i, 1-halo_j:  rows+halo_j,      &
     &          model_levels)                                           &
     &,r_at_v (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,      &
     &          model_levels)                                           &
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                   1-halo_j:rows+halo_j,0:model_levels)           &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &                 1-halo_j:rows+halo_j, model_levels)              &
     &, sec_v_latitude (1-offx:row_length+offx,1-offy:n_rows+offy)      &
     &, tan_v_latitude(row_length,n_rows)                               &
     &, sec_theta_latitude (1-offx:row_length+offx,1-offy:rows+offy)    &
     &, f3_at_v (1-offx:row_length+offx, 1-offy: n_rows+offy)           &
     &,rot_coeff1(row_length,n_rows)                                    &
                                      ! for lam wind rotations
     &,rot_coeff2(row_length,n_rows)                                    &
                                      !     (on B grid)
! Pressure levels (units mb) for output arrays: IN
     &,ucomB_press(ucomB_p_levs)                                        & 
                                        ! for u wind  'B' grid  press
     &,vcomB_press(vcomB_p_levs)                                        & 
                                        ! for v wind  'B' grid  press
     &,ucomp_press(ucomp_p_levs)                                        &
                                        ! for u wind
     &,vcomp_press(vcomp_p_levs)                                        &
                                        ! for v wind
     &,wcomp_press(wcomp_p_levs)                                        &
                                        ! for w wind
     &,testd_press(testd_p_levs)                                        &
                                        ! for test diagnostic3
     &,p_height(p_h_levs)                                               &
                                        ! for diagnostic 108 (press)
     &,theta_height(theta_h_levs)                                       &
                                        ! for diagnostic 119 (theta)
     &,rho_height(rho_h_levs)                                           &
                                        ! for diagnostic 127 (rho)
     &,w_height(w_h_levs)                                               &
                                        ! for diagnostic 142 (w)
     &,u_height(u_h_levs)                                               &
                                        ! for diagnostic 143 (u)
     &,v_height(v_h_levs)                                               &
                                        ! for diagnostic 144 (v)
     &, desired_theta(pv_theta_levs)                                    &
                                      ! for potential vorticity
                                      !  on theta levels
     &, pv_press(pv_press_levs)                                         &
                                      ! for potential vorticity
                                      !  on pressure levels
! Model    levels for output arrays: IN
     &,testd_model(testd_m_levs)        ! for test diagnostic4

      INTEGER                                                           &
     & ucomB_model(ucomB_m_levs)                                        & 
                                        ! for u wind  'B' grid  model
     &,vcomB_model(vcomB_m_levs)                                        & 
                                        ! for v wind  'B' grid  model
     &,Htheta_model(Htheta_m_levs)                                      &
                                        ! for diagnostic 101
     &,Hrho_model(Hrho_m_levs)          ! for diagnostic 102
!   Scalar arguments with intent(InOut):

!   Array  arguments with intent(InOut):

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
! Diagnostic arrays: OUT
      REAL                                                              &
     & ucomB_m(row_length,n_rows,ucomB_m_levs)                          & 
                                               ! u 'B' grid at mod levs
     &,vcomB_m(row_length,n_rows,vcomB_m_levs)                          & 
                                               ! v 'B' grid at mod levs
     &,ucomB_p(row_length,n_rows,ucomB_p_levs)                          & 
                                               ! u 'B' grid at pressures
     &,vcomB_p(row_length,n_rows,vcomB_p_levs)                          & 
                                               ! v 'B' grid at pressures
     &,ucomp_p(row_length,  rows,ucomp_p_levs)                          &
                                               ! u at selected pressures
     &,vcomp_p(row_length,n_rows,ucomp_p_levs)                          &
                                               ! v at selected pressures
     &,wcomp_p(row_length,rows,  wcomp_p_levs)                          &
                                               ! w at selected pressures
     &,u50mB_h(row_length,n_rows)                                       & 
                                               ! u at 50m ht 'B' grid
     &,v50mB_h(row_length,n_rows)                                       & 
                                               ! v at 50m ht 'B' grid
     &,u50m_h(row_length,  rows)                                        &
                                               ! u at 50m height
     &,v50m_h(row_length,n_rows)                                        &
                                               ! v at 50m height
     &, potn_vort_theta(row_length, n_rows, pv_theta_levs)              &
!                                  pv on theta levels
     &, potn_vort_press(row_length, n_rows, pv_press_levs)              &
                                               ! pv on pressure levels
     &,theta_potn_vort(row_length, n_rows)                              &
                                               ! theta on pv+/-2
     &,theta_pv_points(row_length,n_rows,model_levels)                  &
                                                       ! theta at pv poi
     &,pv_mod_levs(row_length,n_rows,model_levels)                      &
                                                   ! pv on model levels
     &,pv_theta_mlev(row_length,rows,model_levels)                      &
                                              ! pv on model theta levels
     &,testdiag1(row_length,n_rows)                                     &
                                               ! testdiag1 -single level
     &,testdiag2(row_length,  rows)                                     &
                                               ! testdiag2 -single level
     &,testdiag3(row_length,rows,testd_p_levs)                          &
                                               ! testdiag3 -press levels
     &,testdiag4(row_length,rows,testd_m_levs)                          &
                                               ! testdiag4 -model levels
     &,rhow(row_length,rows,model_levels)                               &
                                   ! mass flux on model levels
     &,rhow_up(row_length,rows,model_levels)                            &
                                   ! up mass flux on model levels
     &,rhow_down(row_length,rows,model_levels)                          &
                                   ! down mass flux on model levels
     &,rhow_convup(row_length,rows,model_levels)                        &
                                   ! up conv mass flux on model levels
     &,rhow_convdown(row_length,rows,model_levels)                      &
                                   ! down conv mass flux on model levels
     &,rhouw(row_length,rows,model_levels)                              &
                                             ! momentum flux
     &,rhovw(row_length,rows,model_levels)                              &
                                             ! momentum flux
! stash work arrays for height on theta and rho model levels
     &,height_theta_ml(row_length,rows,Htheta_m_levs)                   &
     &,height_rho_ml(row_length,rows,Hrho_m_levs)                       &
! stash work arrays for height level diagnostics
     &,press_h(row_length,rows,p_h_levs)                                &
     &,theta_h(row_length,rows,theta_h_levs)                            &
     &,rho_h(row_length,rows,rho_h_levs)                                &
     &,wcomp_h(row_length,rows,w_h_levs)                                &
     &,ucomp_h(row_length,rows,u_h_levs)                                &
     &,vcomp_h(row_length,n_rows,v_h_levs)                              &
     &,true_density(row_length,rows,model_levels) ! unscaled density


! Local parameters:
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='Dyn_diag')
      REAL                                                              &
     & z_50m                    ! height=50m above earth's surface
      PARAMETER (                                                       &
     & z_50m = 50.                                                      &
     &)

! Local scalars:
!   ErrorStatus
      INTEGER      ErrorStatus                                          &
                                        ! Error flag (0 = OK)
     &,i,j,k                                                            &
                                        ! loop counters
     &,kk                                                               &
                                        ! k level index
     &,interp_order    !  order of vertical interpolation

#if defined(FLUME)
      ! Arguments for FlumeSendDiag
      INTEGER im ! model
      INTEGER is ! section
      INTEGER ie ! item
      INTEGER levels ! number of vertical levels
#endif
      CHARACTER*256                                                     &
     & CMessage         ! Error message if return code >0

      LOGICAL                                                           &
     & LAM              ! T: limited area model

      REAL                                                              &
     & dummy                                                            &
                                ! dummy argument - not referenced
     &,desired_potn_vort                                                &
                                ! value of pv surface for theta
     &,pressure_pa                                                      &
                                ! pressure in pascals
     &,pressure_ex                                                      &
                                ! exner pressure
     &,weight                                                           &
                           ! weight for calculating w_on_rho
     &,desired_r              ! height to interpolate

! Local dynamic arrays:
      REAL                                                              &
     & exner_at_u(row_length,rows,  model_levels)                       &
                                                  !pressure at u points
     &,exner_at_v(row_length,n_rows,model_levels)                       &
                                                  !pressure at v points
     &,exner_at_w(row_length,rows,  model_levels)                       &
                                                  !pressure at w points
     &,exner_at_pv(row_length, n_rows, model_levels)                    &
                                               !pressure at pv points
![Note: exner_at_u,exner_at_v are optionally re-used as height points.]
     &,theta_at_PV(row_length, n_rows, model_levels)                    &
     &,T (1-offx:row_length+offx, 1-offy:rows+offy, model_levels)       &
     &,PV(row_length, n_rows, model_levels)                             &
     &,work_1(row_length,rows)                                          &
                                              ! workspace
     &,work_2(row_length,rows)                                          &
                                              ! workspace
     &,mod_PV(row_length, n_rows, model_levels)                         &
                                                 ! for pv=-2 diag
     &,w_on_rho(row_length,rows,model_levels)                           &
                                          ! vertical velocity on rho pts
     &,u_on_rho(row_length,rows,model_levels)                           &
                                              ! u velocity on rho pts
     &,v_on_rho(row_length,rows,model_levels)                           &
                                              ! v velocity on rho pts
     &,true_rho(row_length,rows,model_levels) ! rho/(r_rho_levels^2)

! Spectral variables

      Real                                                              &
     & spec_3D(row_length,rows,model_levels)                            &
     &,spec_2D(1-offx:row_length+offx,1-offy:rows+offy)

      Integer                                                           &
     & klev                                                             &
     &,gath_proc

! Function & Subroutine calls:
      External Ereport,vert_interp2,vert_interp_mdi,TestDiag            &
     &        ,uC_to_uB,vC_to_vB, calc_pv

!- End of header


! ----------------------------------------------------------------------
! Section 0.  Initialisation.
! ----------------------------------------------------------------------

! Set Error code to zero
      ErrorStatus = 0

! Set order of vertical interpolation
      interp_order = interp_order_linear

! Determine whether limited area model
      IF (model_domain  ==  mt_lam .OR.                                 &
     &    model_domain  ==  mt_cyclic_lam .or.                          &
     &    model_domain  ==  mt_bi_cyclic_lam) THEN
        lam = .true.
      else
        lam = .false.
      ENDIF

!-----------------------------------------------
!    do all swapbounds calls needed
!-----------------------------------------------
      IF(qucomp_p .OR. qucomB_p .OR. qpotn_vort_press) THEN
! DEPENDS ON: swap_bounds
        call swap_bounds(exner_rho_levels,row_length,rows,              &
     &                 model_levels,offx,offy,fld_type_p,.false.)
      endif

      IF(qpotn_vort_theta .OR. qtheta_potn_vort .OR.                    &
     &   qtheta_pv_points .OR. qpotn_vort_press .OR. qpv_mod_levs .OR.  &
     &   qpv_theta_mlev) THEN  
! DEPENDS ON: swap_bounds
        CALL Swap_bounds(theta,row_length,rows,model_levels,            &
     &                   offx,offy,fld_type_p,.true.)
      endif

      IF(qucomB_m .OR. qpotn_vort_press .OR. qtheta_potn_vort .OR.      &
     &   qpv_mod_levs .OR. qrhouw .OR. qpotn_vort_theta .OR.            &
     &   qpv_theta_mlev) THEN  
! DEPENDS ON: swap_bounds
        CALL Swap_bounds(u,row_length,rows,model_levels,                &
     &                 offx,offy,fld_type_u,.true.)
      endif

      IF(qvcomB_m .OR. qpotn_vort_press .OR. qtheta_potn_vort .OR.      &
     &   qpv_mod_levs .OR. qrhovw .OR.qpotn_vort_theta .OR.             &
     &   qpv_theta_mlev) THEN  
! DEPENDS ON: swap_bounds
        CALL Swap_bounds(v,row_length,n_rows,model_levels,              &
     &                 offx,offy,fld_type_v,.true.)
      endif

      If (qpotn_vort_press .OR. qtheta_potn_vort .OR.                   &
     &    qpv_mod_levs .OR. qpotn_vort_theta .OR. qpv_theta_mlev) THEN
! DEPENDS ON: swap_bounds
        CALL Swap_bounds(rho,row_length,rows,model_levels,              &
     &                   offx,offy,fld_type_p,.true.)
      endif

      IF(qucomp_p.OR.qucomB_p) THEN
!  Calculate exner at u points. Store in exner_at_u
! Need call to swapbounds for exner_rho before interpolation

      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
             exner_at_u(i,j,k) = 0.5 *                                  &
     &             (exner_rho_levels(i,j,k) + exner_rho_levels(i+1,j,k))
           ENDDO ! i
        ENDDO ! j
      ENDDO ! k

      ENDIF ! on STASHflag

      IF(qvcomp_p.OR.qvcomB_p) THEN
!  Calculate exner at v points. Store in exner_at_v

      DO k = 1, model_levels
        DO j = 1, n_rows
          DO i = 1, row_length
              exner_at_v(i,j,k) = 0.5 *                                 &
     &             (exner_rho_levels(i,j,k) + exner_rho_levels(i,j+1,k))
           ENDDO ! i
        ENDDO ! j
      ENDDO ! k

      ENDIF ! on STASHflag

      IF(qwcomp_p) THEN
!  Calculate exner at w points. Store in exner_at_w

      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
             exner_at_w(i,j,k) = exner_theta_levels(i,j,k)
           ENDDO ! i
        ENDDO ! j
      ENDDO ! k

      ENDIF ! on STASHflag

      IF(qpotn_vort_theta.OR.qtheta_potn_vort.OR.                       &
     &   qtheta_pv_points) THEN

! Calculate theta at PV points. Store in theta_at_pv
! first interpolate theta to rho levels. Use linear interpolation.
! Store in T as work space.

        Do j = 1-offy, rows+offy
          Do i = 1-offx, row_length+offx
             T(i,j,1) = theta(i,j,1)
          End Do
        End Do

        Do k = 2, model_levels
          Do j = 1-offy, rows+offy
            Do i = 1-offx, row_length+offx
                    T(i,j,k) = (theta(i,j,k) *                          &
     &                        (r_rho_levels(i,j,k) -                    &
     &                         r_theta_levels(i,j,k-1) ) +              &
     &                        theta(i,j,k-1) *                          &
     &                        (r_theta_levels(i,j,k) -                  &
     &                         r_rho_levels(i,j,k) ) ) /                &
     &                        (r_theta_levels(i,j,k) -                  &
     &                         r_theta_levels(i,j,k-1) )
            End Do
          End Do
        End Do

        Do k = 1, model_levels
          Do j = 1, n_rows
            Do i = 1, row_length
               theta_at_pv(i,j,k) = .25*(T(i+1,j,k) + T(i,j,k) +        &
     &                                   T(i+1,j+1,k) + T(i,j+1,k) )
            End Do
          End Do
        End Do
      ENDIF ! on STASHflag

      IF(qpotn_vort_press) THEN

! Interpolate exner onto PV points.

        Do k = 1, model_levels
          Do j = 1, n_rows
            Do i = 1, row_length
               exner_at_pv(i,j,k) = .25*(exner_rho_levels(i+1,j,k) +    &
     &                                   exner_rho_levels(i,j,k) +      &
     &                                   exner_rho_levels(i+1,j+1,k) +  &
     &                                   exner_rho_levels(i,j+1,k) )
            End Do
          End Do
        End Do
      ENDIF ! on STASHflag

! STASH items 002,003     : u,v wind components on model levs 'B' grid
! ----------------------------------------------------------------------

      IF(qucomB_m) THEN

        DO  k=1,ucomB_m_levs
! Perform simple horizontal interpolation from 'C' to 'B' grid
! Halos already populated, so interpolate directly:

          kk = ucomB_model(k) ! selected model level

          DO j=1,n_rows
          DO i=1,row_length
            ucomB_m(i,j,k) = (u(i,j,kk)+u(i,j+1,kk)) * 0.5
          ENDDO ! i
          ENDDO ! j

        ENDDO  ! k model levels loop

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 2
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = ucomB_m_levs
          CALL flumeSendDiag (ucomB_m,im,is,ie,row_length,n_rows,levels)
        END IF
#endif

      ENDIF ! on STASHflag

      IF(qvcomB_m) THEN

        DO  k=1,vcomB_m_levs
! Perform simple horizontal interpolation from 'C' to 'B' grid
! Halos already populated, so interpolate directly:

          kk = vcomB_model(k) ! selected model level

          DO j=1,n_rows
          DO i=1,row_length
            vcomB_m(i,j,k) = (v(i,j,kk)+v(i+1,j,kk)) * 0.5
          ENDDO ! i
          ENDDO ! j

        ENDDO  ! k model levels loop

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 3
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = vcomB_m_levs
          CALL flumeSendDiag (vcomB_m,im,is,ie,row_length,n_rows,levels)
        END IF
#endif

      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH items 201,202     : u,v wind components on p surfaces 'B' grid
! ----------------------------------------------------------------------

      IF(qucomB_p) THEN

        DO  k=1,ucomB_p_levs

          pressure_pa = ucomB_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
! DEPENDS ON: vert_interp2
          CALL vert_interp2 (u, row_length, rows, model_levels          &
     &                          ,pressure_ex                            &
     &                          ,offx, offy, 0, 0                       &
     &                          ,exner_at_u, interp_order               &
     &                          ,work_1 )

! Perform simple horizontal interpolation from 'C' to 'B' grid

! DEPENDS ON: uc_to_ub
          CALL  uC_to_uB(work_1,                                        &
     &                   row_length,rows,n_rows,1,offx,offy,            &
     &                   ucomB_p(1,1,k))
        ENDDO  ! k pressure levels loop

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 201
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = ucomB_p_levs
          CALL flumeSendDiag (ucomB_p,im,is,ie,row_length,n_rows,levels)
        END IF
#endif

      ENDIF ! on STASHflag

      IF(qvcomB_p) THEN
        DO  k=1,vcomB_p_levs

          pressure_pa = vcomB_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
! DEPENDS ON: vert_interp2
          CALL vert_interp2 (v, row_length, n_rows, model_levels        &
     &                          ,pressure_ex                            &
     &                          ,offx, offy, 0, 0                       &
     &                          ,exner_at_v, interp_order               &
     &                          ,work_1 )

! Perform simple horizontal interpolation from 'C' to 'B' grid

! DEPENDS ON: vc_to_vb
          CALL  vC_to_vB(work_1,                                        &
     &                   row_length,n_rows,1,offx,offy,                 &
     &                   vcomB_p(1,1,k))

! Rotate winds from model to standard lat-long grid
          IF (lam .AND. rot_uvcomB_p) THEN
! First check valid requests: implicit assumption that u components
! and v components are requested for the same pressure levels
            IF(qucomB_p .AND. ucomB_press(k) == vcomB_press(k)) THEN

              DO j=1,n_rows
              DO i=1,row_length
                work_1(i,j) = ucomB_p(i,j,k)
                work_2(i,j) = vcomB_p(i,j,k)
              ENDDO
              ENDDO

! Rotation calculation on B grid
! DEPENDS ON: w_eqtoll
              CALL W_EqtoLL(rot_coeff1,rot_coeff2,work_1,work_2,        &
     &         ucomB_p(1,1,k),vcomB_p(1,1,k),v_field_size,v_field_size)

            ELSE

              ErrorStatus = -1        ! Warning
              Cmessage='wind diagnostics cannot be rotated: u and v '   &
     &           //'requested components on pressure levels must match'
! DEPENDS ON: ereport
              CALL Ereport(Routinename,ErrorStatus,Cmessage)
              exit                    ! jump out of levels loop

            ENDIF       ! Check valid request

          ENDIF  ! (lam wind rotation)

        ENDDO  ! k pressure levels loop

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 202
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = vcomB_p_levs
          CALL flumeSendDiag (vcomB_p,im,is,ie,row_length,n_rows,levels)
        END IF
#endif

      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! Calculate PV for use with STASH items 229,214,215,217
! ----------------------------------------------------------------------

      If (qpotn_vort_press .OR. qtheta_potn_vort .OR.                   &
     &    qpv_mod_levs .or. qpotn_vort_theta) THEN

! Calculate PV on model levels.

! DEPENDS ON: calc_pv
        Call Calc_PV                                                    &
     &            (u, v, theta, rho,                                    &
     &             r_theta_levels, r_rho_levels,                        &
     &             r_at_u, r_at_v,                                      &
     &             sec_v_latitude, tan_v_latitude,                      &
     &             sec_theta_latitude, f3_at_v,                         &
     &             delta_lambda, delta_phi,                             &
     &             row_length, rows, n_rows, model_levels,              &
     &             offx, offy, halo_i, halo_j,                          &
     &             at_extremity,                                        &
     &             PV)

      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH item 229 : potential vorticity on pressure levels
! ----------------------------------------------------------------------

      If (qpotn_vort_press) THEN

        Do k = 1, pv_press_levs

          pressure_pa = pv_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
! DEPENDS ON: vert_interp2
          Call vert_interp2 (pv, row_length, n_rows,                    &
     &                       model_levels,                              &
     &                       pressure_ex,                               &
     &                       0, 0, 0, 0,                                &
     &                       exner_at_pv, interp_order_cubic,           &
     &                       potn_vort_press(1,1,k) )

        ENDDO  ! k pressure levels loop

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 229
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = pv_press_levs
          CALL flumeSendDiag &
                   (potn_vort_press,im,is,ie,row_length,n_rows,levels)
        END IF
#endif

      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH item 214 : potential vorticity on theta levels
! ----------------------------------------------------------------------

      If (qpotn_vort_theta) THEN

        Do k = 1, pv_theta_levs

! DEPENDS ON: vert_interp_mdi
          Call vert_interp_mdi (PV, row_length, n_rows,                 &
     &                          model_levels,                           &
     &                          desired_theta(k),                       &
     &                          0, 0, 0, 0,                             &
     &                          theta_at_pv, interp_order_cubic,        &
     &                          rmdi, potn_vort_theta(1,1,k) )

        ENDDO  ! k theta levels loop

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 214
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = pv_theta_levs
          CALL flumeSendDiag&
            (potn_vort_theta,im,is,ie,row_length,n_rows,levels)
        END IF
#endif

      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH items 243,244,242 : u,v,w wind components on pressure surfaces
! ----------------------------------------------------------------------

      IF(qucomp_p) THEN

        DO  k=1,ucomp_p_levs

          pressure_pa = ucomp_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
! DEPENDS ON: vert_interp2
          CALL vert_interp2 (u, row_length, rows, model_levels          &
     &                          ,pressure_ex                            &
     &                          ,offx, offy, 0, 0                       &
     &                          ,exner_at_u, interp_order               &
     &                          ,ucomp_p(1,1,k) )

        ENDDO  ! k pressure levels loop

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 243
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = ucomp_p_levs
          CALL flumeSendDiag (ucomp_p,im,is,ie,row_length,n_rows,levels)
        END IF
#endif

      ENDIF ! on STASHflag

      IF(qvcomp_p) THEN
        DO  k=1,vcomp_p_levs

          pressure_pa = vcomp_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
! DEPENDS ON: vert_interp2
          CALL vert_interp2 (v, row_length, n_rows, model_levels        &
     &                          ,pressure_ex                            &
     &                          ,offx, offy, 0, 0                       &
     &                          ,exner_at_v, interp_order               &
     &                          ,vcomp_p(1,1,k) )

        ENDDO  ! k pressure levels loop

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 244
          levels = vcomp_p_levs
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          CALL flumeSendDiag (vcomp_p,im,is,ie,row_length,n_rows,levels)
        END IF
#endif

      ENDIF ! on STASHflag

      IF(qwcomp_p) THEN
        DO  k=1,wcomp_p_levs

          pressure_pa = wcomp_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
! DEPENDS ON: vert_interp2
          CALL vert_interp2 (w, row_length, rows, model_levels          &
     &                          ,pressure_ex                            &
     &                          ,offx, offy, 0, 0                       &
     &                          ,exner_at_w, interp_order               &
     &                          ,wcomp_p(1,1,k) )

        ENDDO  ! k pressure levels loop

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 242
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = wcomp_p_levs
          CALL flumeSendDiag (wcomp_p,im,is,ie,row_length,n_rows,levels)
        END IF
#endif

      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH items 212,213 : u,v wind components at 50m height 'B' grid
! STASH items 245,246 : u,v wind components at 50m height
! ----------------------------------------------------------------------

! Restrict interpolation to boundary levels, which should always
! greatly exceed levels in the vicinity of 50m.


      IF(qu50mB_h.OR.qu50m_h) THEN
! Generate height field above orography for lower (ie boundary) levels
! at u pts: (re-use array exner_at_u for workspace)
        DO k=1,bl_levels
          DO j=1,rows
            DO i=1,row_length
              exner_at_u(i,j,k)= r_at_u(i,j,k) - r_at_u(i,j,1)
            ENDDO ! i
          ENDDO ! j
        ENDDO ! k

! DEPENDS ON: vert_interp_mdi
         CALL vert_interp_mdi (u, row_length, rows,                     &
     &                         bl_levels, z_50m,                        &
     &                         offx, offy,                              &
     &                         0, 0,                                    &
     &                         exner_at_u, interp_order,                &
     &                         rmdi, work_1 )
         IF(qu50m_h) THEN
           DO j=1,rows
             DO i=1,row_length
               u50m_h(i,j) = work_1(i,j)
             ENDDO ! i
           ENDDO ! j

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 245
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          CALL flumeSendDiag (u50m_h,im,is,ie,row_length,n_rows)
        END IF
#endif

         ENDIF ! on STASHflag qu50m_h

       IF(qu50mB_h) THEN
! Perform simple horizontal interpolation from 'C' to 'B' grid

! DEPENDS ON: uc_to_ub
            CALL  uC_to_uB(work_1,                                      &
     &                     row_length,rows,n_rows,1,offx,offy,          &
     &                     u50mB_h)

#if defined(FLUME)     
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 212
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          CALL flumeSendDiag (u50mB_h,im,is,ie,row_length,n_rows)     
        END IF
#endif
     
       ENDIF ! on STASHflag qu50mB_h

      ENDIF ! on STASHflags  qu50mB_h.OR.qu50m_h

      IF(qv50mB_h.OR.qv50m_h) THEN

! Generate height field above orography for lower (ie boundary) levels
! at v pts: (re-use array exner_at_v for workspace)
        DO k=1,bl_levels
          DO j=1,n_rows
            DO i=1,row_length
              exner_at_v(i,j,k)= r_at_v(i,j,k) - r_at_v(i,j,1)
            ENDDO ! i
          ENDDO ! j
        ENDDO ! k

! DEPENDS ON: vert_interp_mdi
         CALL vert_interp_mdi (v, row_length, n_rows,                   &
     &                         bl_levels, z_50m,                        &
     &                         offx, offy,                              &
     &                         0, 0,                                    &
     &                         exner_at_v, interp_order,                &
     &                         rmdi, work_1)
         IF(qv50m_h) THEN
           DO j=1,n_rows
             DO i=1,row_length
               v50m_h(i,j) = work_1(i,j)
             ENDDO ! i
           ENDDO ! j

#if defined(FLUME)
!        FLUME-STASH 
         IF (Flume_run) THEN
           im = 1
           is = 15
           ie = 246
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
           model=im
           section=is
           item=ie
           CALL flumeSendDiag (v50m_h,im,is,ie,row_length,n_rows)
         ENDIF
#endif

         ENDIF ! on STASHflag qv50m_h

        IF(qv50mB_h) THEN
! Perform simple horizontal interpolation from 'C' to 'B' grid

! DEPENDS ON: vc_to_vb
            CALL  vC_to_vB(work_1,                                      &
     &                     row_length,n_rows,1,offx,offy,               &
     &                     v50mB_h)

#if defined(FLUME)
!        FLUME-STASH 
         IF (Flume_run) THEN
           im = 1
           is = 15
           ie = 213
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
           model=im
           section=is
           item=ie
           CALL flumeSendDiag (v50mB_h,im,is,ie,row_length,n_rows)
         ENDIF
#endif

         ENDIF ! on STASHflag qv50mB_h

      ENDIF ! on STASHflags qv50mB_h.OR.qv50m_h


! ----------------------------------------------------------------------
! STASH items 231,232,233,234: Test diagnostics 1-4
! ----------------------------------------------------------------------

      IF (qdia1.OR.qdia2.OR.qdia3.OR.qdia4) THEN

! DEPENDS ON: testdiag
        CALL TestDiag(                                                  &
     &  theta_field_size,v_field_size,rows,n_rows,row_length            &
     & ,ew_space,ns_space,first_lat,first_long,phi_pole,lambda_pole     &
     & ,lam                                                             &
     & ,testd_press,testd_p_levs                                        &
     & ,testd_model,testd_m_levs,forecast_hrs                           &
     & ,testdiag1,testdiag2,testdiag3,testdiag4                         &
     & ,qdia1,qdia2,qdia3,qdia4)

      ENDIF ! on STASHflag

! ---------------------------------------------------------------------
! STASH item 215 : theta on potential vorticity = +/-2 surface
! ---------------------------------------------------------------------

      IF (qtheta_potn_vort) THEN
        desired_potn_vort = 0.0000020    ! set pv surface to pv=2x10^-6

! Take the absolute value of pv so that interpolation is to pv = +/- 2

        Do k = 1, model_levels
          Do j = 1, n_rows
            Do i = 1, row_length
              mod_PV(i,j,k) = ABS(PV(i,j,k))
            End Do
          End Do
        End Do

! Interpolate theta onto pv=2 surface

! DEPENDS ON: vert_interp_mdi
        Call vert_interp_mdi (theta_at_pv, row_length, n_rows,          &
     &                         model_levels,                            &
     &                         desired_potn_vort,                       &
     &                         0, 0, 0, 0,                              &
     &                         mod_PV, interp_order_linear,             &
     &                         rmdi, theta_potn_vort(1,1) )

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 215
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          CALL flumeSendDiag (theta_potn_vort,im,is,ie,row_length,n_rows)
        ENDIF
#endif

      ENDIF ! on STASHflag

! ---------------------------------------------------------------------
! STASH item 216 : theta at potential vorticity points
! ---------------------------------------------------------------------

      IF (qtheta_pv_points) THEN

         Do k = 1, model_levels
            Do j = 1, n_rows
               Do i = 1, row_length
                  theta_pv_points(i,j,k) = theta_at_pv(i,j,k)
               End Do
            End Do
         End Do

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 216
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = model_levels
          CALL flumeSendDiag &
                  (theta_pv_points,im,is,ie,row_length,n_rows,levels)
        ENDIF
#endif

      ENDIF ! on STASHflag

! ---------------------------------------------------------------------
! STASH item 217 : potential vorticity on model levels
! ---------------------------------------------------------------------

      IF (qpv_mod_levs) THEN

         Do k = 1, model_levels
            Do j = 1, n_rows
               Do i = 1, row_length
                  pv_mod_levs(i,j,k) = PV(i,j,k)
               End Do
            End Do
         End Do

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 217
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = model_levels
          CALL flumeSendDiag&
             (pv_mod_levs,im,is,ie,row_length,n_rows,levels)
        ENDIF
#endif

      ENDIF ! on STASHflag

! ---------------------------------------------------------------------
! STASH item 218 : potential vorticity on model theta points and levels
! ---------------------------------------------------------------------

      IF (qpv_theta_mlev) THEN

! DEPENDS ON: calc_pv_at_theta
        Call Calc_PV_at_theta(u, v, theta, rho,                         &
     &                        r_theta_levels, r_rho_levels,             &
     &                        r_at_u, r_at_v,                           &
     &                        sec_v_latitude, tan_v_latitude,           &
     &                        sec_theta_latitude, f3_at_v,              &
     &                        delta_lambda, delta_phi,                  &
     &                        model_domain,                             &
     &                        pv_theta_mlev)

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 218
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = model_levels
          CALL flumeSendDiag&
             (pv_theta_mlev,im,is,ie,row_length,rows,levels)
        END IF
#endif

      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH items 260-266 : mass and momentum fluxes
! 260 = mass flux = rhow
! 261 = rhouw
! 262 = rhovw
! 263 = upward mass flux = rhow with w >0m/s
! 264 = downward mass flux = rhow with w <0m/s
! 265 = upward convective mass flux = rhow with w> 1m/s
! 266 = upward convective mass flux = rhow with w < -1m/s
! ----------------------------------------------------------------------
      IF (qrhow .OR. qrhouw .OR. qrhovw .OR. qrhow_up .OR.              &
     &     qrhow_down .OR. qrhowc_up .OR. qrhowc_down) THEN
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              weight = (r_theta_levels(i,j,k) - r_rho_levels(i,j,k))/   &
     &                 (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1))
              w_on_rho(i,j,k)= weight        * w(i,j,k-1) +             &
     &                        (1.0 - weight) * w(i,j,k)
              true_rho(i,j,k)=rho(i,j,k)/                               &
     &                        (r_rho_levels(i,j,k)*r_rho_levels(i,j,k))
              rhow(i,j,k) = true_rho(i,j,k)*w_on_rho(i,j,k)
              IF (qrhouw) THEN
                u_on_rho(i,j,k)=0.5*(u(i-1,j,k)+u(i,j,k))
                rhouw(i,j,k)=rhow(i,j,k)*u_on_rho(i,j,k)
              ENDIF
              IF (qrhovw) THEN
                v_on_rho(i,j,k)=0.5*(v(i,j-1,k)+v(i,j,k))
                rhovw(i,j,k)=rhow(i,j,k)*v_on_rho(i,j,k)
              ENDIF
              IF (qrhow_up) THEN
                IF (w_on_rho(i,j,k)  >   0.0) THEN
                  rhow_up(i,j,k)=rhow(i,j,k)
                ELSE
                  rhow_up(i,j,k)=0.0
                ENDIF
              ENDIF
              IF (qrhow_down) THEN
                IF (w_on_rho(i,j,k)  <   0.0) THEN
                  rhow_down(i,j,k) = rhow(i,j,k)
                ELSE
                  rhow_down(i,j,k)=0.0
                ENDIF
              ENDIF
              IF (qrhowc_up) THEN
                IF (w_on_rho(i,j,k)  >   1.0) THEN
                  rhow_convup(i,j,k)=rhow(i,j,k)
                ELSE
                  rhow_convup(i,j,k)=0.0
                ENDIF
              ENDIF
              IF (qrhowc_down) THEN
                IF (w_on_rho(i,j,k)  <   -1.0) THEN
                  rhow_convdown(i,j,k) = rhow(i,j,k)
                ELSE
                  rhow_convdown(i,j,k)=0.0
                ENDIF
              ENDIF
            End Do
          End Do
        End Do

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN     
          im = 1
          is = 15
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          levels = model_levels
          IF (qrhow) THEN
            ie = 260
            item=ie
            CALL flumeSendDiag (rhow,im,is,ie,row_length,n_rows,levels)
          END IF 
          IF (qrhouw) THEN
            ie = 261
            item=ie
            CALL flumeSendDiag (rhouw,im,is,ie,row_length,n_rows,levels)
          END IF 
          IF (qrhovw) THEN
            ie = 262
            item=ie
            CALL flumeSendDiag (rhovw,im,is,ie,row_length,n_rows,levels)
          END IF 
          IF (qrhow_up) THEN
            ie = 263
            item=ie
            CALL flumeSendDiag (rhow_up,im,is,ie,row_length,n_rows,levels)
          END IF 
          IF (qrhow_down) THEN
            ie = 264
            item=ie
            CALL flumeSendDiag (rhow_down,im,is,ie,row_length,n_rows,levels)
          END IF
          IF (qrhowc_up) THEN
            ie = 265
            item=ie
            CALL flumeSendDiag &
                     (rhow_convup,im,is,ie,row_length,n_rows,levels)
          END IF
          IF (qrhowc_down) THEN
            ie = 266
            item=ie
            CALL flumeSendDiag &
                     (rhow_convdown,im,is,ie,row_length,n_rows,levels)
          END IF        
        END IF ! Flume_run
#endif

      ENDIF ! on STASHflag
!-----------------------------------------------------------------------
! STASH items 101, 102, 108, 119, 127, 142, 143, 144
!-----------------------------------------------------------------------
! Height on theta-levels diagnostic

      If ( qHtheta_ml ) Then

         Do k=1, Htheta_m_levs
            Do j=1, rows
               Do i=1, row_length
                  kk = Htheta_model(k)
                  height_theta_ml(i,j,k)                                &
     &                        = r_theta_levels(i,j,kk) - Earth_Radius
               Enddo
            Enddo
         Enddo

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 101
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = Htheta_m_levs
          CALL flumeSendDiag&
             (height_theta_ml,im,is,ie,row_length,n_rows,levels)
        END IF
#endif
      Endif

! Height on rho-levels diagnostic

      If ( qHrho_ml ) Then

         Do k=1, Hrho_m_levs
            Do j=1, rows
               Do i=1, row_length
                  kk = Hrho_model(k)
                  height_rho_ml(i,j,k)                                  &
     &                 = r_rho_levels(i,j,kk) - Earth_Radius
               Enddo
            Enddo
         Enddo

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 102
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = Hrho_m_levs
          CALL flumeSendDiag &
              (height_rho_ml,im,is,ie,row_length,n_rows,levels)
        END IF
#endif
      End If

! pressure on height-levels diagnostic

      If ( qpress_H ) Then

         Do k=1, p_h_levs

            desired_r = p_height(k) + Earth_Radius

! DEPENDS ON: vert_interp_mdi
            Call vert_interp_mdi (exner_rho_levels,row_length,rows,     &
     &           model_levels,desired_r,offx,offy,halo_i,halo_j,        &
     &           r_rho_levels,interp_order,rmdi,press_h(:,:,k) )

! Convert to standard pressure

            Do j=1, rows
               Do i=1, row_length
                  If ( press_h(i,j,k)  /=  rmdi ) Then
                     press_h(i,j,k) = p_zero *                          &
     &                    press_h(i,j,k)**(1./kappa)
                  End If
               End Do
            End Do

         End do  ! end k-loop

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 108
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = p_h_levs
          ! rows=73, n_rows=72
          CALL flumeSendDiag (press_h,im,is,ie,row_length,rows,levels)
        END IF
#endif
      End if

! potential temperature on height-levels diagnostic

      If ( qtheta_H ) Then

         Do k=1,theta_h_levs

            desired_r = theta_height(k) + Earth_Radius

! DEPENDS ON: vert_interp_mdi
            Call vert_interp_mdi ( theta,row_length,rows,               &
     &           model_levels,desired_r,offx,offy,halo_i,halo_j,        &
     &           r_theta_levels(1-halo_i,1-halo_j,1),interp_order,      &
     &           rmdi,theta_h(:,:,k) )

         End do

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 119
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = theta_h_levs
          CALL flumeSendDiag (theta_h,im,is,ie,row_length,n_rows,levels)
        END IF
#endif
      End if

! density (rho) on height-levels diagnostic

      If ( qrho_H ) Then

         Do k=1,rho_h_levs

            desired_r = rho_height(k) + Earth_Radius

! DEPENDS ON: vert_interp_mdi
            Call vert_interp_mdi (rho,row_length,rows,model_levels,     &
     &           desired_r,offx,offy,halo_i,halo_j,r_rho_levels,        &
     &           interp_order,rmdi,rho_h(:,:,k))

! Convert to true density: divide by r_rho_levels**2

            Do j=1, rows
               Do i=1, row_length
                  If ( rho_h(i,j,k)  /=  rmdi ) Then
                     rho_h(i,j,k) = rho_h(i,j,k) /                      &
     &                    (r_rho_levels(i,j,k)*r_rho_levels(i,j,k))
                  End if
               End do
            End do

         End do

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 127
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = rho_h_levs
          CALL flumeSendDiag (rho_h,im,is,ie,row_length,n_rows,levels)
        END IF
#endif
      End if

! w on height-levels diagnostic

      If ( qw_H ) Then

         Do k=1, w_h_levs

            desired_r = w_height(k) + Earth_Radius

! DEPENDS ON: vert_interp_mdi
            Call vert_interp_mdi (w,row_length,rows,model_levels+1,     &
     &           desired_r,offx,offy,halo_i,halo_j,r_theta_levels,      &
     &           interp_order,rmdi,wcomp_h(:,:,k) )

         End do

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 142
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = w_h_levs
          CALL flumeSendDiag (wcomp_h,im,is,ie,row_length,n_rows,levels)
        END IF
#endif
      End if

! u on height-levels diagnostic

      If ( qu_H ) Then

         Do k=1, u_h_levs

            desired_r = u_height(k) + Earth_Radius

! DEPENDS ON: vert_interp_mdi
            Call vert_interp_mdi (u,row_length,rows,model_levels,       &
     &           desired_r,offx,offy,halo_i,halo_j,                     &
     &           r_at_u,interp_order,rmdi,ucomp_h(:,:,k))

         End do

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 143
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = u_h_levs
          CALL flumeSendDiag (ucomp_h,im,is,ie,row_length,n_rows,levels)
        END IF
#endif
      End if

! v on height-levels diagnostic

      If ( qv_H ) Then

         Do k=1, v_h_levs

            desired_r = v_height(k) + Earth_Radius

! DEPENDS ON: vert_interp_mdi
            Call vert_interp_mdi (v,row_length,n_rows,model_levels,     &
     &           desired_r,offx,offy,halo_i,halo_j,r_at_v,              &
     &           interp_order,rmdi,vcomp_h(:,:,k))

         End do

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 144
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = v_h_levs
          CALL flumeSendDiag (vcomp_h,im,is,ie,row_length,n_rows,levels)
        END IF
#endif
      End if

      IF (spec_w) then
        If (model_domain  ==  4) then
          gath_proc=0
          do klev=1,model_levels
! DEPENDS ON: calc_spectra
            CALL CALC_SPECTRA(w(1-offx,1-offy,klev),spec_2D,            &
     &             row_length+2*offx,rows+2*offy,                       &
     &             global_row_length,global_rows,                       &
     &             fld_type_p,halo_type_single,                         &
     &             gath_proc)
            do j = 1,rows
              do i = 1,row_length
                spec_3D(i,j,klev)=spec_2D(i,j)
              enddo
            enddo
          enddo
        Else
          write(6,*)'The spectra is not set up for this domain'
        Endif
      ENDIF

! -----------------------------------------------------
! stash item 271  : true density on model_levels
! -----------------------------------------------------
      if (qtrue_density) then
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              true_density(i,j,k)=rho(i,j,k)/                           &
     &                        (r_rho_levels(i,j,k)*r_rho_levels(i,j,k))
            end do
          end do
        end do

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 15
          ie = 271
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = model_levels
          CALL flumeSendDiag (true_density,im,is,ie,row_length,n_rows,levels)
        END IF
#endif
      endif

! ---------------------------------------------------------------------
! Check error condition
      IF(ErrorStatus >  0) THEN
! DEPENDS ON: ereport
         CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

      RETURN
      END SUBROUTINE Dyn_diag
#endif
