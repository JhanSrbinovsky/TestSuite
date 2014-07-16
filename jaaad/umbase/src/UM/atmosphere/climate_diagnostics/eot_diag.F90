#if defined(A30_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculates the required EOT climate diagnostics, and fills STAHSWORK.
!
! Subroutine Interface:
      SUBROUTINE Eot_diag(                                              &
! Primary data: in
     &  pstar,p,rho,u,v,w,theta,q,qcl,qcf                               &
     &  ,p_theta_levels                                                 &
     &  ,exner_rho_levels,exner_theta_levels                            &
     &  ,energy_corr_now                                                &
     &  ,inc_u, inc_v, inc_w, inc_t                                     &
     &  ,inc_q, inc_qcl, inc_qcf, inc_rho                               &
! Grid sizes and definition: in
     &  ,rows,n_rows,row_length,model_levels,wet_levels,bl_levels       &
     &  ,theta_field_size,u_field_size,v_field_size                     &
     &  ,eta_theta_levels,eta_rho_levels                                &
     &  ,r_theta_levels,r_rho_levels                                    &
     &  ,sin_theta_latitude,cos_theta_latitude                          &
     &  ,sin_theta_longitude,cos_theta_longitude                        &
     &  ,Model_domain                                                   &
     &  ,delta_lambda                                                   &
                                ! grid longitude spacing in radians
     &  ,delta_phi                                                      &
                                ! grid latitude  spacing in radians
     &, global_row_length                                               &
! Grid coordinates: in
     &  ,ew_space,ns_space,first_lat,first_long,phi_pole,lambda_pole    &
! Pre-calculated grid associated arrays: in
     &  ,r_at_u,r_at_v                                                  &
! Time information: in
     &  ,forecast_hrs                                                   &
! Pressure levels for output arrays: in
     &  ,u_press,v_press,w_press                                        &
     &  ,T_PRESS,Q_PRESS,RH_PRESS                                       &
     &  ,Z_PRESS,OM_PRESS,HEAVY_PRESS                                   &
! Flags to request each diagnostic output field: in
     &  ,qsf,  qv_m, qw_m,  qt_m, qq_m, qz_m, qke_m                     &
     &  ,qu_mm,qv_mm,qw_mm ,qt_mm,qq_mm,qz_mm,qke_mm                    &
     &  ,qteot_m,qwbig_m,qwbig2_m,qRH_m,qdry_mass_m                     &
     &  ,qt_inc,qq_inc,qqcl_inc                                         &
     &  ,qqcf_inc,qu_inc,qv_inc,qw_inc,qrho_inc                         &
     &  ,qt2_inc,qq2_inc,qqcl2_inc,qqcf2_inc,qu2_inc,qv2_inc            &
     &  ,qw2_inc,qrho2_inc                                              &
     &  ,qu_p,qv_p,qw_p                                                 &
     &  ,qt_p,qq_p,qrh_p                                                &
     &  ,qz_p,qom_p                                                     &
     &  ,qheavy_p,qtv_p,qtvom_p                                         &
     &  ,qtot_ke,qtot_ke_w,qtcm                                         &
     &  ,qtcmq,qtcqcl,qtcqcf                                            &
     &  ,qtmf_u,qtmf_v,qtmf_w,qmtorque                                  &
     &  ,qam_m1,qam_m2,qam_m3                                           &
     &  ,qam_w1,qam_w2,qam_w3                                           &
     &  ,qpstar,qpstar_uv                                               &
     &  ,qencorr,qtot_cvt,qtot_gr                                       &
     &  ,qcol_ugz, qcol_vgz, qcol_wgz, qcol_uT,  qcol_vt,  qcol_wt      &
     &  ,qcol_uq,  qcol_vq,  qcol_wq , qcol_uv,  qcol_uw,  qcol_vw      &
     &  ,qcol_uKe, qcol_vKe, qcol_wKe, qcol_u,   qcol_v,   qcol_w       &
     &  ,qm_spd_ew,qm_spd_ns                                            &
     &  ,qtrop_p,qtrop_t,qtrop_z,qtrop_icao                             &
! Diagnostics lengths: in
     &  ,u_m_levs,v_m_levs,w_m_levs,t_m_levs,q_m_levs,z_m_levs          &
     &  ,ke_m_levs,u_mm_levs,v_mm_levs,w_mm_levs,t_mm_levs,q_mm_levs    &
     &  ,z_mm_levs,ke_mm_levs                                           &
     &  ,teot_m_levs,wbig_m_levs,wbig2_m_levs,Rh_m_levs,dry_mass_m_levs &
     &  ,t_inc_levs,q_inc_levs,qcl_inc_levs                             &
     &  ,qcf_inc_levs,u_inc_levs,v_inc_levs                             &
     &  ,w_inc_levs,rho_inc_levs                                        &
     &  ,t2_inc_levs,q2_inc_levs,qcl2_inc_levs,qcf2_inc_levs            &
     &  ,u2_inc_levs,v2_inc_levs,w2_inc_levs,rho2_inc_levs              &
     &  ,u_p_levs,v_p_levs,w_p_levs                                     &
     &  ,t_p_levs,q_p_levs,rh_p_levs                                    &
     &  ,z_p_levs,om_p_levs,heavy_p_levs                                &
     &  ,tv_p_levs,tvom_p_levs                                          &
     &  ,prod_p_levs                                                    &
                     ! pressure levels required for each product
     &  ,npress_diags,nmodel_diags                                      &
! Tropopause index bounds
     &  ,min_trop_level,max_trop_level                                  &
! Model levels indicies: in
     &  ,u_m_list,v_m_list,w_m_list,t_m_list,q_m_list,z_m_list,ke_m_list&
     &  ,u_mm_list,v_mm_list,w_mm_list,t_mm_list,q_mm_list,z_mm_list    &
     &  ,ke_mm_list,teot_m_list,wbig_m_list,wbig2_m_list,RH_m_list      &
     &  ,dry_mass_m_list                                                &
     &  ,t_inc_list,q_inc_list,qcl_inc_list                             &
     &  ,qcf_inc_list,u_inc_list,v_inc_list                             &
     &  ,w_inc_list,rho_inc_list                                        &
     &  ,t2_inc_list,q2_inc_list,qcl2_inc_list,qcf2_inc_list            &
     &  ,u2_inc_list,v2_inc_list,w2_inc_list,rho2_inc_list              &
     &  ,prod_m_list                                                    &
                     ! index of model levels required for product
! Product levels: in
     &  ,prod_ind                                                       &
! Diagnostic arrays: out
     &  ,si,stashwork                                                   &
     &  ,u_m,v_m,w_m ,t_m,q_m,z_m,ke_m                                  &
     &  ,u_mm,v_mm,w_mm,t_mm,q_mm,z_mm,ke_mm                            &
     &  ,teot_m ,wbig_m,wbig2_m,RH_m,dry_mass_m                         &
     &  ,t_inc,q_inc,qcl_inc,qcf_inc,u_inc,v_inc,w_inc,rho_inc          &
     &  ,t2_inc,q2_inc,qcl2_inc,qcf2_inc,u2_inc,v2_inc,w2_inc,rho2_inc  &
     &  ,u_p,v_p,w_p,t_p,q_p,rh_p ,z_p,om_p,heavy_p                     &
     &  ,tv_p,tvom_p                                                    &
     &  ,tot_ke,tot_ke_w,tcm                                            &
     &  ,tcmq,tcqcl,tcqcf                                               &
     &  ,tmf_u,tmf_v,tmf_w,mtorque                                      &
     &  ,am_m1,am_m2,am_m3                                              &
     &  ,am_w1,am_w2,am_w3                                              &
     &  ,o_pstar,o_pstar_uv                                             &
     &  ,encorr,tot_cvt,tot_gr                                          &
     &  ,col_ugz, col_vgz, col_wgz, col_uT,  col_vt,  col_wt            &
     &  ,col_uq,  col_vq,  col_wq , col_uv,  col_uw,  col_vw            &
     &  ,col_uKe, col_vKe, col_wKe, col_u,   col_v,   col_w             &
     &  ,m_spd_ew,m_spd_ns                                              &
     &  ,trop_p,trop_t,trop_z,trop_icao                                 &
     &  )
      IMPLICIT NONE
!
! Description:
!   Calculate eot-related diagnostics - held in STASH section 30 -
!   The stash items are split into centuries, as follows:
!
!   000s   Standard Model level fields and products
!   100s   Other Model level fields
!   200s   Standard Pressure level fields and products
!   300s   Other Pressure level fields
!   400s   Surface fields and Column integrals
!
!   Diagnostics on model rho points
!    1  u component of wind
!    2  v component of wind
!    3  w component of wind
!    4  T temperature
!    5  q specific humidity
!    6  z height
!    7  KE = 0.5(u*u+v*v+w*w)
!   Mass weighted products on model rho levels
!   000+10*(1st code)+(2nd code) ie
!    11 u*u
!    12 u*v
!    17 u*ke
!   mass weighted fields on rho levels
!   where (field * rho *r*r*dr/a*a) a - radius of earth
!   101  u component of wind
!   102  v component of wind
!   103  w component of wind
!   104  T temperature
!   105  q specific humidity
!   106  z height
!   107  KE = 0.5(u*u+v*v+w*w)
!
!   Diagnostics on pressure Pressure surfaces.
!   STASH item
!   201 u component of wind on pressure surfaces
!   202 v component of wind on pressure surfaces
!   203 w component of wind on pressure surfaces
!   204 t temperature on pressure surfaces
!   205 q specific humidity on pressure surfaces
!   206 rh relative humidity on pressure surfaces
!   207 z geopotential on pressure surfaces
!   208 om omega on pressure surfaces
!
!   Products are also produced, the stashcodes are
!   200+10*(1st code)+(2nd code) ie
!   211 UU
!   212 UV
!   224 VT
!
!
!  Further user and technical documentation can be found under:
!  http://www-hc/~hadsu/
!
!
! Method:
!   Required level lists and logical switches are determined by the
!   calling routine from STASH requests and STASHflags.
!   Code 201-208 are processed sepertely, product codes
!   are calculated in a loop over single pressure surface fields
!
! Current Code Owner: <Simon Wilson>
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.1  02/02/2000   1st itteration
!  5.2  20/03/2001    Change to use exner pressures instead of p in
!                     vertical interpolation and change call to
!                     vert_interp_mdi2 to vert_interp2
!
!  5.2  15/08/2000   Pstar diagnostic added Simon Wilson
!
!  5.2  22/09/2000 Extra diagnostics and bug fixes. SS Wilson
!
!  5.2  19/03/2001   Linear interpolate in exner instead of height.
!                    (T on p levels). C. Wilson
!  5.2  19/03/2001   Linear interpolate in exner instead of height.
!                    C. Wilson
!LL  5.3  24/09/01  Portability changes.    Z. Gardner
!  5.3  27/04/2001   Correct non-reprod u&v on pressure levels
!                    caused by acw4f502 changes.  R A Stratton.
!
!  5.3  08/08/2001   Tropopause diagnostics added. DGH Tan
!  5.4  24/04/2002   Add RH on model levels, wbig, squared increment
!                    diagnostics.  R A Stratton.
!        4/07/2002   Add column integral of CvT and gr R A Stratton.
!  5.4     28/08/02    Bug Fix (Bi-cyclic LAM)           Carol Roadnight

!  5.4     01/10/02    correction to rho_dry                  A. Malcolm
!  5.4  11/03/2002   Remove comment on same line as #include
!                                                 S. Carroll
!  5.5  13/12/2002   Add wbig 0.1. Sort out products on model levels
!                    and add more. R A Stratton.
!       18/12/2002   Add more column integrals. R A Stratton
!
!  5.5  23/01/2003   Correct mountain torque diagnostic and add
!                    surface pressure drag diagnostics. S. Webster
!  6.0  09/02/2004   Call swap_bounds and fill_external_halos
!                    for exner_rho_levels and exner_theta_levels
!                    A. A. Dickinson
!  6.1  08/07/04   send P_theta_levels into vert_h_onto_p
!                                                    Michael Hughes
!  6.1   21/05/04    correct bit reproducibility problem     A. Malcolm
!  6.2   15/08/05    Free format fixes. P.Selwood
!  6.2   02/02/06 Modify argument list in call to tropin with SCM
!                 dummy variables to retain consistency in
!                 argument list.                       R. Wong
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
#include "typstsz.h"

! Subroutine arguments
!   Scalar arguments with intent(in):
! Grid sizes:
      INTEGER                                                           &
     &  rows,n_rows,row_length,model_levels,wet_levels,bl_levels        &
     &  ,theta_field_size,u_field_size,v_field_size                     &
     &  ,Model_domain                                                   &
     &  ,global_row_length ! global row length
! Grid coordinates: in
      REAL                                                              &
     & ew_space,ns_space,first_lat,first_long,phi_pole,lambda_pole
! Time information: in
      INTEGER                                                           &
     & forecast_hrs       ! T+forecast_hrs: hours after analysis time
! Diagnostics lengths: in
      INTEGER, intent(IN) ::                                            &
     & u_m_levs                                                         &
                      ! No of levels for output of U_M
     &,v_m_levs                                                         &
                      ! No of levels for output of V_M
     &,w_m_levs                                                         &
                      ! No of levels for output of W_M
     &,t_m_levs                                                         &
                      ! No of levels for output of T_M
     &,q_m_levs                                                         &
                      ! No of levels for output of Q_M
     &,z_m_levs                                                         &
                      ! No of levels for output of z_M
     &,ke_m_levs                                                        &
                       ! No of levels for output of ke_M
     &,wbig_m_levs                                                      &
                      ! No of levels for output of wbig_m
     &,wbig2_m_levs                                                     &
                       ! No of levels for output of wbig2_m
     &,Rh_m_levs                                                        &
                      ! No of levels for output of RH
     &,dry_mass_m_levs                                                  &
                       ! No of levels for output of dry mass weighting
     &,u_mm_levs                                                        &
                       ! No of levels for output of U_MM
     &,v_mm_levs                                                        &
                       ! No of levels for output of V_MM
     &,w_mm_levs                                                        &
                       ! No of levels for output of W_MM
     &,t_mm_levs                                                        &
                       ! No of levels for output of T_MM
     &,q_mm_levs                                                        &
                       ! No of levels for output of Q_MM
     &,z_mm_levs                                                        &
                       ! No of levels for output of z_MM
     &,ke_mm_levs                                                       &
                        ! No of levels for output of KE_MM
     &,teot_m_levs                                                      &
                         ! No of levels for output of TEOT_M
     &,u_inc_levs                                                       &
                        ! No of levels for output of U_INC
     &,v_inc_levs                                                       &
                        ! No of levels for output of V_INC
     &,w_inc_levs                                                       &
                        ! No of levels for output of W_INC
     &,t_inc_levs                                                       &
                        ! No of levels for output of T_INC
     &,q_inc_levs                                                       &
                        ! No of levels for output of Q_INC
     &,qcl_inc_levs                                                     &
                          ! No of levels for output of Q_INC
     &,qcf_inc_levs                                                     &
                          ! No of levels for output of Q_INC
     &,rho_inc_levs                                                     &
                          ! No of levels for output of Q_INC
     &,u2_inc_levs                                                      &
                         ! No of levels for output of U2_INC
     &,v2_inc_levs                                                      &
                         ! No of levels for output of V2_INC
     &,w2_inc_levs                                                      &
                         ! No of levels for output of W2_INC
     &,t2_inc_levs                                                      &
                         ! No of levels for output of T2_INC
     &,q2_inc_levs                                                      &
                         ! No of levels for output of Q2_INC
     &,qcl2_inc_levs                                                    &
                           ! No of levels for output of Qcl2_INC
     &,qcf2_inc_levs                                                    &
                           ! No of levels for output of Qcf2_INC
     &,rho2_inc_levs                                                    &
                           ! No of levels for output of rho2_INC
     &,u_p_levs                                                         &
                      ! No of levels for output of U_P
     &,v_p_levs                                                         &
                      ! No of levels for output of V_P
     &,w_p_levs                                                         &
                      ! No of levels for output of W_P
     &,t_p_levs                                                         &
                      ! No of levels for output of T_P
     &,q_p_levs                                                         &
                      ! No of levels for output of Q_P
     &,rh_p_levs                                                        &
                       ! No of levels for output of RH_P
     &,z_p_levs                                                         &
                      ! No of levels for output of Z_P
     &,om_p_levs                                                        &
                       ! No of levels for output of OM_P
     &,heavy_p_levs                                                     &
                       ! No of levels for output of HEAVY_P
     &,tv_p_levs                                                        &
                    ! No of levels for output of TV_P
     &,tvom_p_levs    ! No of levels for output of TVOM_P

! Flags to request each diagnostic output field: IN
      LOGICAL, intent(IN) ::                                            &
     &   qv_m                                                           &
               ! Flag for V wind component on model levels
     &  ,qw_m                                                           &
                 ! Flag for W wind component on model level
     &  ,qt_m                                                           &
                 ! Flag for T on model level
     &  ,qq_m                                                           &
                 ! Flag for Q on model level
     &  ,qz_m                                                           &
                 ! Flag for z on model level
     &  ,qke_m                                                          &
                  ! Flag for ke on model level
     &  ,qwbig_m                                                        &
                 ! Flag for wbig on model levels
     &  ,qwbig2_m                                                       &
                  ! Flag for wbig0.1 on model levels
     &  ,qRH_m                                                          &
                 ! Flag for RH on model levels
     &  ,qdry_mass_m                                                    &
                     ! Flag for dry mass weighting on model levels
     &  ,qu_mm                                                          &
                 ! Flag for U wind component on model levels
     &  ,qv_mm                                                          &
                ! Flag for V wind component on model levels
     &  ,qw_mm                                                          &
                  ! Flag for W wind component on model level
     &  ,qt_mm                                                          &
                  ! Flag for T on model level
     &  ,qq_mm                                                          &
                  ! Flag for Q on model level
     &  ,qz_mm                                                          &
                  ! Flag for z on model level
     &  ,qke_mm                                                         &
                   ! Flag for KE on model level
     &  ,qteot_m                                                        &
                    ! Flag for T at eot on model level
     &  ,qu_inc                                                         &
                  ! Flag for U wind component on model levels
     &  ,qv_inc                                                         &
                 ! Flag for V wind component on model levels
     &  ,qw_inc                                                         &
                   ! Flag for W wind component on model level
     &  ,qt_inc                                                         &
                   ! Flag for T on model level
     &  ,qq_inc                                                         &
                   ! Flag for Q on model level
     &  ,qqcl_inc                                                       &
                     ! Flag for QCL on model level
     &  ,qqcf_inc                                                       &
                     ! Flag for QCF on model level
     &  ,qrho_inc                                                       &
                     ! Flag for RHO on model level
     &  ,qu2_inc                                                        &
                    ! Flag for U2 wind component on model levels
     &  ,qv2_inc                                                        &
                    ! Flag for V2 wind component on model levels
     &  ,qw2_inc                                                        &
                    ! Flag for W2 wind component on model level
     &  ,qt2_inc                                                        &
                    ! Flag for T2 on model level
     &  ,qq2_inc                                                        &
                    ! Flag for Q2 on model level
     &  ,qqcl2_inc                                                      &
                      ! Flag for Qcl2 on model level
     &  ,qqcf2_inc                                                      &
                      ! Flag for Qcf2 on model level
     &  ,qrho2_inc                                                      &
                      ! Flag for RHO2 on model level
     &  ,qu_p                                                           &
                ! Flag for U wind component on pressure levels
     &  ,qv_p                                                           &
               ! Flag for V wind component on pressure levels
     &  ,qw_p                                                           &
                 ! Flag for W wind component on pressure level
     &  ,qt_p                                                           &
                 ! Flag for T on pressure level
     &  ,qq_p                                                           &
                 ! Flag for Q on pressure level
     &  ,qrh_p                                                          &
                 ! Flag for RH on pressure level
     &  ,qz_p                                                           &
                 ! Flag for Z on pressure level
     &  ,qom_p                                                          &
                 ! Flag for Omega on pressure level
     &  ,qheavy_p                                                       &
                   ! Flag for Heavyside fn on pressure level
     &  ,qtv_p                                                          &
                   ! Flag for virtual temp  on pressure level
     &  ,qtvom_p                                                        &
                   ! Flag for virtual temp*om on pressure level
     &  ,qsf(nitems)                                                    &
                     ! Flag for all diags
     &  ,qtot_ke                                                        &
     &  ,qtot_ke_w                                                      &
     &  ,qtot_cvt                                                       &
     &  ,qtot_gr                                                        &
     &  ,qtcm                                                           &
     &  ,qtcmq                                                          &
     &  ,qtcqcl                                                         &
     &  ,qtcqcf                                                         &
     &  ,qtmf_u                                                         &
     &  ,qtmf_v                                                         &
     &  ,qtmf_w                                                         &
     &  ,qmtorque                                                       &
     &  ,qm_spd_ew                                                      &
     &  ,qm_spd_ns                                                      &
     &  ,qam_m1                                                         &
     &  ,qam_m2                                                         &
     &  ,qam_m3                                                         &
     &  ,qam_w1                                                         &
     &  ,qam_w2                                                         &
     &  ,qam_w3                                                         &
     &  ,qpstar                                                         &
     &  ,qpstar_uv                                                      &
     &  ,qencorr                                                        &
     &  ,qtrop_p                                                        &
     &  ,qtrop_t                                                        &
     &  ,qtrop_z                                                        &
     &  ,qtrop_icao                                                     &
     &  ,qcol_ugz, qcol_vgz, qcol_wgz                                   &
     &  ,qcol_uT,  qcol_vt,  qcol_wt                                    &
     &  ,qcol_uq,  qcol_vq,  qcol_wq                                    &
     &  ,qcol_uv,  qcol_uw,  qcol_vw                                    &
     &  ,qcol_uKe, qcol_vKe, qcol_wKe                                   &
     &  ,qcol_u,   qcol_v,   qcol_w


!   Array  arguments with intent(in):
! Primary data: IN
      REAL                                                              &
     &  u(1-offx:row_length+offx, 1-offy:rows+offy,model_levels)        &
     &  ,v(1-offx:row_length+offx, 1-offy:n_rows+offy,model_levels)     &
     &  ,w(1-offx:row_length+offx, 1-offy:rows+offy,0:model_levels)     &
     &  ,theta(1-offx:row_length+offx, 1-offy:rows+offy,model_levels)   &
     &  ,q(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,wet_levels)  &
     &  ,qcl(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,wet_levels)&
     &  ,qcf(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,wet_levels)&
     &  ,rho(1-offx:row_length+offx, 1-offy:rows+offy,model_levels)     &
     &  ,p(1-offx:row_length+offx, 1-offy:rows+offy,model_levels)       &
     &  ,pstar(row_length, rows)                                        &
     &  ,p_theta_levels                                                 &
     &  (1-offx:row_length+offx, 1-offy:rows+offy,    model_levels)     &
     &  ,exner_rho_levels                                               &
     &  (1-offx:row_length+offx, 1-offy:rows+offy,    model_levels)     &
     &  ,exner_theta_levels                                             &
     &  (1-offx:row_length+offx, 1-offy:rows+offy,    model_levels)     &
     &  ,energy_corr_now                                                &
     &,  inc_u(1-offx:row_length+offx, 1-offy:rows+offy,                &
     &        model_levels)                                             &
     &, inc_v(1-offx:row_length+offx, 1-offy:n_rows+offy,               &
     &        model_levels)                                             &
     &, inc_w(row_length, rows, model_levels)                           &
     &, inc_t(1-offx:row_length+offx,                                   &
     &               1-offy:rows+offy, model_levels)                    &
     &, inc_q(1-halo_i:row_length+halo_i,                               &
     &               1-halo_j:rows+halo_j,wet_levels)                   &
     &, inc_qcl(1-halo_i:row_length+halo_i,                             &
     &               1-halo_j:rows+halo_j,wet_levels)                   &
     &, inc_qcf(1-halo_i:row_length+halo_i,                             &
     &               1-halo_j:rows+halo_j,wet_levels)                   &
     &, inc_rho(1-offx:row_length+offx,                                 &
     &             1-offy:rows+offy, model_levels)                      &


! Vertical grid definition: IN
     &  ,eta_theta_levels(0:model_levels)                               &
                                          !vertical grid for theta vars
     &  ,eta_rho_levels    (model_levels)                               &
                                          ! vertical grid for rho   vars
     &  ,r_theta_levels  (1-halo_i:row_length+halo_i,                   &
     &  1-halo_j:rows+halo_j,       0:model_levels)                     &
     &  ,r_rho_levels    (1-halo_i:row_length+halo_i,                   &
     &  1-halo_j:rows+halo_j,         model_levels)                     &
! Pre-calculated grid associated arrays: IN
     &  ,r_at_u (1-halo_i:row_length+halo_i, 1-halo_j:  rows+halo_j,    &
     &  model_levels)                                                   &
     &  ,r_at_v (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,    &
     &          model_levels)                                           &
     &  ,cos_theta_latitude (1-Offx:row_length+Offx,                    &
     &                      1-Offy:rows+Offy)                           &
                                                 ! cos(latitude)
     &  ,cos_theta_longitude (row_length, rows)                         &
     &  ,sin_theta_longitude (row_length, rows)                         &
     &  ,sin_theta_latitude (row_length, rows)                          &
! Pressure levels (units mb) for output arrays: IN
     &  ,u_press(u_p_levs)                                              &
                                ! for u wind
     &  ,v_press(v_p_levs)                                              &
                                ! for v wind
     &  ,w_press(w_p_levs)                                              &
                                ! for w wind
     &  ,t_press(t_p_levs)                                              &
                                ! for t
     &  ,q_press(q_p_levs)                                              &
                                ! for q
     &  ,rh_press(rh_p_levs)                                            &
                                ! for rh
     &  ,z_press(z_p_levs)                                              &
                                ! for z
     &  ,om_press(om_p_levs)                                            &
                                ! for om
     &  ,heavy_press(heavy_p_levs)                                      &
                                   ! for heavyside fn
     &  ,delta_lambda                                                   &
                                ! grid longitude spacing in radians
     &  ,delta_phi              ! grid latitude  spacing in radians


      integer, intent(in) ::                                            &
     &  npress_diags                                                    &
                                !# of diags on pressure levels
     &  ,nmodel_diags                                                   &
                                 !# of diags on model levels
     & ,prod_IND(NUM_STASH_LEVELS*2,npress_diags,npress_diags)          &
     & ,prod_p_levs(npress_diags,npress_diags)

! Tropopause index bounds.  Should be consistent with
! a statement of the form #include <typcona/typcona.h>.
! and also with Deck O3EXP1.84--87
      Integer, intent(IN) :: min_trop_level
!         Minimum permitted level of the tropopause
      Integer, intent(IN) :: max_trop_level
!         Maximum permitted level of the tropopause

      LOGICAL, Intent(IN) :: U_M_LIST(MODEL_LEVELS)                     &
     &  ,V_M_LIST(MODEL_LEVELS)                                         &
     &  ,W_M_LIST(MODEL_LEVELS),wbig_m_list(MODEL_LEVELS)               &
     &  ,wbig2_m_list(MODEL_LEVELS)                                     &
     &  ,T_M_LIST(MODEL_LEVELS),RH_m_list(MODEL_LEVELS)                 &
     &  ,Q_M_LIST(MODEL_LEVELS)                                         &
     &  ,Z_M_LIST(MODEL_LEVELS)                                         &
     &  ,KE_M_LIST(MODEL_LEVELS)                                        &
     &  ,U_MM_LIST(MODEL_LEVELS)                                        &
     &  ,V_MM_LIST(MODEL_LEVELS)                                        &
     &  ,W_MM_LIST(MODEL_LEVELS)                                        &
     &  ,T_MM_LIST(MODEL_LEVELS)                                        &
     &  ,Q_MM_LIST(MODEL_LEVELS)                                        &
     &  ,Z_MM_LIST(MODEL_LEVELS)                                        &
     &  ,KE_MM_LIST(MODEL_LEVELS)                                       &
     &  ,TEOT_M_LIST(MODEL_LEVELS)                                      &
     &  ,Dry_mass_m_list(MODEL_LEVELS)                                  &
     &  ,U_INC_LIST(MODEL_LEVELS)                                       &
     &  ,V_INC_LIST(MODEL_LEVELS)                                       &
     &  ,W_INC_LIST(MODEL_LEVELS)                                       &
     &  ,T_INC_LIST(MODEL_LEVELS)                                       &
     &  ,Q_INC_LIST(MODEL_LEVELS)                                       &
     &  ,QCL_INC_LIST(MODEL_LEVELS)                                     &
     &  ,QCF_INC_LIST(MODEL_LEVELS)                                     &
     &  ,RHO_INC_LIST(MODEL_LEVELS)                                     &
     &  ,U2_INC_LIST(MODEL_LEVELS)                                      &
     &  ,V2_INC_LIST(MODEL_LEVELS)                                      &
     &  ,W2_INC_LIST(MODEL_LEVELS)                                      &
     &  ,T2_INC_LIST(MODEL_LEVELS)                                      &
     &  ,Q2_INC_LIST(MODEL_LEVELS)                                      &
     &  ,QCL2_INC_LIST(MODEL_LEVELS)                                    &
     &  ,QCF2_INC_LIST(MODEL_LEVELS)                                    &
     &  ,RHO2_INC_LIST(MODEL_LEVELS)

      logical prod_m_list(model_levels,nmodel_diags,nmodel_diags)


!   Scalar arguments with intent(InOut):

!   Array  arguments with intent(InOut):
      integer si(nitems)
      real stashwork(*)

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
! Diagnostic arrays: OUT
      REAL, Intent(out) ::                                              &
     &  u_m(row_length,rows,  u_m_levs)                                 &
                                        ! u at selected model lev
     &  ,v_m(row_length,rows,  v_m_levs)                                &
                                         ! v at selected model lev
     &  ,w_m(row_length,rows,  w_m_levs)                                &
                                         ! w at selected model lev
     &  ,t_m(row_length,rows,  t_m_levs)                                &
                                         ! t at selected model lev
     &  ,q_m(row_length,rows,  q_m_levs)                                &
                                         ! q at selected model lev
     &  ,z_m(row_length,rows,  z_m_levs)                                &
                                         ! z at selected model lev
     &  ,ke_m(row_length,rows, ke_m_levs)                               &
                                          ! ke at selected model lev
     &  ,wbig_m(row_length,rows,wbig_m_levs)                            &
                                             !wbig at selected mod lev
     &  ,wbig2_m(row_length,rows,wbig2_m_levs)                          &
                                               !wbig at selected mod lev
     &  ,RH_m(row_length,rows,Rh_m_levs)                                &
                                           ! RH at selected model lev
     &  ,dry_mass_m(row_length,rows,dry_mass_m_levs)                    &
                                                     ! dry mass mod lev
     &  ,u_mm(row_length,rows,  u_mm_levs)                              &
                                           ! u at selected model lev
     &  ,v_mm(row_length,rows,  v_mm_levs)                              &
                                           ! v at selected model lev
     &  ,w_mm(row_length,rows,  w_mm_levs)                              &
                                           ! w at selected model lev
     &  ,t_mm(row_length,rows,  t_mm_levs)                              &
                                           ! t at selected model lev
     &  ,q_mm(row_length,rows,  q_mm_levs)                              &
                                           ! q at selected model lev
     &  ,z_mm(row_length,rows,  z_mm_levs)                              &
                                           ! z at selected model lev
     &  ,ke_mm(row_length,rows, ke_mm_levs)                             &
                                            ! ke at selected model lev
     &  ,u_inc(row_length,rows,  u_inc_levs)                            &
                                             ! u at selected model lev
     &  ,v_inc(row_length,n_rows,  v_inc_levs)                          &
                                               ! v at selected model lev
     &  ,w_inc(row_length,rows,  w_inc_levs)                            &
                                             ! w at selected model lev
     &  ,t_inc(row_length,rows,  t_inc_levs)                            &
                                             ! t at selected model lev
     &  ,q_inc(row_length,rows,  q_inc_levs)                            &
                                             ! q at selected model lev
     &  ,qcl_inc(row_length,rows,  qcl_inc_levs)                        &
                                                 ! qcl at  model lev
     &  ,qcf_inc(row_length,rows,  qcf_inc_levs)                        &
                                                 ! qcf at  model lev
     &  ,rho_inc(row_length,rows,  rho_inc_levs)                        &
                                                 ! rho at  model lev
     &  ,u2_inc(row_length,rows, u2_inc_levs)                           &
                                              ! u2 at selected mod lev
     &  ,v2_inc(row_length,n_rows, v2_inc_levs)                         &
                                                ! v2 at selected mod lev
     &  ,w2_inc(row_length,rows, w2_inc_levs)                           &
                                              ! w2 at selected mod lev
     &  ,t2_inc(row_length,rows, t2_inc_levs)                           &
                                              ! t2 at selected mod lev
     &  ,q2_inc(row_length,rows, q2_inc_levs)                           &
                                              ! q2 at selected mod lev
     &  ,qcl2_inc(row_length,rows, qcl2_inc_levs)                       &
                                                  ! qcl2 at selected lev
     &  ,qcf2_inc(row_length,rows, qcf2_inc_levs)                       &
                                                  ! qcf2 at selected lev
     &  ,rho2_inc(row_length,rows, rho2_inc_levs)                       &
                                                  ! rho2 at  model lev
     &  ,teot_m(row_length,rows,  teot_m_levs)                          &
                                               ! T at eot/selected
     &  ,u_p(row_length,n_rows,  u_p_levs)                              &
                                           ! u at selected pressures
     &  ,v_p(row_length,n_rows,  v_p_levs)                              &
                                           ! v at selected pressures
     &  ,w_p(row_length,n_rows,  w_p_levs)                              &
                                           ! w at selected pressures
     &  ,t_p(row_length,n_rows,  t_p_levs)                              &
                                           ! t at selected pressures
     &  ,q_p(row_length,n_rows,  q_p_levs)                              &
                                           ! q at selected pressures
     &  ,rh_p(row_length,n_rows,  rh_p_levs)                            &
                                             ! rh at selected pressures
     &  ,z_p(row_length,n_rows,  z_p_levs)                              &
                                           ! z at selected pressures
     &  ,om_p(row_length,n_rows,  om_p_levs)                            &
                                             ! om at selected pressures
     &  ,heavy_p(row_length,n_rows,  heavy_p_levs)                      &
                                                   ! heavy at pressures
     &  ,tv_p(row_length,n_rows,  tv_p_levs)                            &
                                             ! tv at pressures
     &  ,tvom_p(row_length,n_rows,  tvom_p_levs)                        &
                                                 ! tv at pressures
     &  ,tot_ke(row_length,rows)                                        &
     &  ,tot_ke_w(row_length,rows)                                      &
     &  ,tot_cvt(row_length,rows)                                       &
     &  ,tot_gr(row_length,rows)                                        &
     &  ,tcm(row_length,rows)                                           &
     &  ,tcmq(row_length,rows)                                          &
     &  ,tcqcl(row_length,rows)                                         &
     &  ,tcqcf(row_length,rows)                                         &
     &  ,tmf_u(row_length,rows)                                         &
     &  ,tmf_v(row_length,rows)                                         &
     &  ,tmf_w(row_length,rows)                                         &
     &  ,mtorque(row_length,rows)                                       &
     &  ,m_spd_ew(row_length,rows)                                      &
     &  ,m_spd_ns(row_length,n_rows)                                    &
     &  ,am_m1(row_length,rows)                                         &
     &  ,am_m2(row_length,rows)                                         &
     &  ,am_m3(row_length,rows)                                         &
     &  ,am_w1(row_length,rows)                                         &
     &  ,am_w2(row_length,rows)                                         &
     &  ,am_w3(row_length,rows)                                         &
     &  ,o_pstar(row_length,rows)                                       &
     &  ,o_pstar_uv(row_length,n_rows)                                  &
     &  ,encorr(row_length,rows)                                        &
     &  ,trop_p(row_length,rows)                                        &
     &  ,trop_t(row_length,rows)                                        &
     &  ,trop_z(row_length,rows)                                        &
     &  ,trop_icao(row_length,rows)                                     &
     &  ,col_ugz(row_length,rows)                                       &
     &  ,col_vgz(row_length,rows)                                       &
     &  ,col_wgz(row_length,rows)                                       &
     &  ,col_ut(row_length,rows)                                        &
     &  ,col_vt(row_length,rows)                                        &
     &  ,col_wt(row_length,rows)                                        &
     &  ,col_uq(row_length,rows)                                        &
     &  ,col_vq(row_length,rows)                                        &
     &  ,col_wq(row_length,rows)                                        &
     &  ,col_uv(row_length,rows)                                        &
     &  ,col_uw(row_length,rows)                                        &
     &  ,col_vw(row_length,rows)                                        &
     &  ,col_uke(row_length,rows)                                       &
     &  ,col_vke(row_length,rows)                                       &
     &  ,col_wke(row_length,rows)                                       &
     &  ,col_u(row_length,rows)                                         &
     &  ,col_v(row_length,rows)                                         &
     &  ,col_w(row_length,rows)


! SCM Dummy variables to keep call to tropin consistent.
       REAL                                                             &
     & scm_dummy_1d(1,1)                                                &
     &,scm_dummy_2d(1,1,0:model_levels)

! Local parameters:
      Logical ::                                                        &
     &   qu_m                                                           &
                ! Flag for U wind component on model levels
     &  ,qprod_p(nitems) ! Flag for prods on pressure level

      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='Eot_diag')

! Local scalars:
!   ErrorStatus
      INTEGER      ErrorStatus                                          &
                                        ! Error flag (0 = OK)
     &  ,i,j,k                                                          &
                                ! loop counters
     &  ,gi,gj                                                          &
     &  ,interp_order                                                   &
                                !  order of vertical interpolation
     &  ,iprod,jprod,code30,                                            &
     &  count

      CHARACTER*256                                                     &
     & CMessage         ! Error message if return code >0

      LOGICAL                                                           &
     & LAM              ! T: limited area model

      REAL                                                              &
     &  dummy                                                           &
                                ! dummy argument - not referenced
     &  ,pressure_pa                                                    &
                                ! pressure in pascals
     &  ,pressure_ex                                                    &
                                ! exner pressure
     &  ,factor                                                         &
                ! scaling factor
     &  ,grid_factor                                                    &
                     ! grid factor
     &  ,cos_long                                                       &
     &  ,sin_long                                                       &
     &  ,sin_lat                                                        &
     &  ,weight1,weight2,weight3,ww1,ww2,temp

! Local dynamic arrays:
      REAL ::                                                           &
     &   exner_at_u(row_length,rows, model_levels)                      &
                                                    ! exner at u points
     &  ,exner_at_v(row_length,n_rows,model_levels)                     &
                                                    ! exner at v points
     &  ,pstar_uv(row_length,n_rows)                                    &
                                     ! pstar at uv points
     &  ,work_1(row_length,rows)                                        &
                                 ! workspace
     &  ,work_2(row_length,rows,model_levels)                           &
                                              ! workspace
     &  ,work_rho(row_length,rows,model_levels)                         &
                                                ! workspace
     &  ,delr_rho(row_length,rows,model_levels)                         &
                                                ! difference
                                ! between rho levels
     &  ,pstar_halo(1-offx:row_length+offx, 1-offy:rows+offy)           &
     &  ,field_rho(row_length,rows,model_levels,nmodel_diags)           &
                                            ! Fields on rho grid points
     &  ,mass_we(row_length,rows,model_levels)    ! mass weighting

! Tropopause variables adapted from subroutine O3_to_3D
      Integer :: trindx(row_length, rows)
!         Points to the layer boundary below that containing the
!         tropopause
      Real :: T_filled_halo(1-offx:row_length+offx,                     &
     &                      1-offy:rows+offy,                           &
     &                      model_levels)



      REAL                                                              &
! prod1 at selected pressures
     &  field1_p(row_length,n_rows  ,u_p_levs)                          &
! prod2 at selected pressures
     &  ,field2_p(row_length,n_rows  ,u_p_levs)                         &
! prod1 model levels
     &  ,field1_m(row_length,rows  ,model_levels)                       &
! prod2 model levels
     &  ,field2_m(row_length,rows  ,model_levels)                       &
     &  ,longitude(row_length,rows)                                     &
     &  ,latitude(row_length,rows)                                      &
     &  ,cossq(row_length,rows)                                         &
     &  ,slcp(row_length,rows)                                          &
     &  ,clcp(row_length,rows)                                          &
     &  ,spslcp(row_length,rows)                                        &
     &  ,spclcp(row_length,rows)                                        &
     &  ,r3dr,rocos

      Real                                                              &
     &  mag_vector_np(model_levels)                                     &
                                     ! magnitude of the vector wind NP
     &, dir_vector_np(model_levels)                                     &
                                     ! direction of the vector wind NP
     &, mag_vector_sp(model_levels)                                     &
                                     ! magnitude of the vector wind SP
     &, dir_vector_sp(model_levels)  ! direction of the vector wind SP


! Local dynamic arrays:
      REAL                                                              &
     &  p_at_theta                                                      &
     &  (row_length,rows,model_levels)                                  &
                                       ! pressure at theta points
     &  ,exner_at_theta                                                 &
     &   (row_length,rows,model_levels)                                 &
                                          ! exner at theta points
     &  ,rh(row_length,rows,model_levels)                               &
                                          ! workspace for RH
     &  ,T(row_length,rows,model_levels)                                &
                                         ! T at theta levels
     &  , true_rho_halo(1-offx:row_length+offx, 1-offy:rows+offy        &
     &  , model_levels)                                                 &
     &  , tmp_data_halo(1-offx:row_length+offx, 1-offy:rows+offy        &
     &  , model_levels)


      logical :: qneed_u,qneed_v,qneed_w,qneed_t,qneed_q,qneed_ke       &
     &  , qneed_tropopause,qneed_z
!
! Tropopause variable
      logical QMODEL_HALF_HEIGHT ! Interim local definition
                                 ! Later pass as argument?

      INTEGER N_SUMS
      PARAMETER (N_SUMS=10)   ! number of global summations

      REAL                                                              &
     & vert_int_array(row_length,rows,n_sums)                           &
                                                 ! array to be summed
     &,flux_array(row_length,rows,3)       ! flux array

      Real                                                              &
     &  rho_dry(row_length,rows,model_levels)                           &
                                                    ! rho dry
     &, dry_to_wet(row_length,rows,model_levels)                        &
     &, wbig,wbig2                    ! restriction on w

      Integer                                                           &
     & n_proc                                                           &
                          ! Total number of processors       (dummy)
     &, n_procx                                                         &
                           ! Number of processors in longitude (dummy)
     &, n_procy            ! Number of processors in latitude  (dummy)


! pointers for vert_int_array
      Integer, parameter ::                                             &
     & ip_dry_mass = 1                                                  &
                            ! dry mass
     &,ip_wet_mass = 2                                                  &
                            ! wet mass
     &,ip_cvT      = 3                                                  &
                            ! cvT
     &,ip_gr       = 4                                                  &
                            ! gr
     &,ip_keu      = 5                                                  &
                            ! keu
     &,ip_kev      = 6                                                  &
                            ! kev
     &,ip_kew      = 7                                                  &
                            ! kew
     &,ip_q        = 8                                                  &
                            ! q
     &,ip_qcl      = 9                                                  &
                            ! qcl
     &,ip_qcf      =10                                                  &
                            ! qcf
     &,ip_qu   = 1                                                      &
                        ! qu
     &,ip_qv   = 2                                                      &
                        ! qv
     &,ip_qw   = 3      ! qw

! pointers for rho model level fields and products
      Integer, parameter ::                                             &
     &  irho_u = 1                                                      &
     &, irho_v = 2                                                      &
     &, irho_w = 3                                                      &
     &, irho_t = 4                                                      &
     &, irho_q = 5                                                      &
     &, irho_z = 6                                                      &
     &, irho_ke= 7

! Function & Subroutine calls:
      External Ereport,vert_interp2,vert_interp_mdi,TestDiag,           &
     &  uC_to_uB,vC_to_vB,pC_to_pB,T_vert_interp_to_p,                  &
     &  qsat,P_TO_T,VERT_ENG_MASSQ,swap_bounds
      External Fill_external_halos, tropin

!- End of header


! ----------------------------------------------------------------------
! Section 0.  Initialisation.
! ----------------------------------------------------------------------

! Set Error code to zero
      ErrorStatus = 0

!Set up variable obtained from sf(1,30)
      qprod_p=qsf
      qu_m=qsf(001)

! Factor to scale mass weighted fields by
! this is done because otherwise fields like u*z can take the same value
! as the missing data indicator which prevents meaning.

      factor=1.e-5

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

! Set SCM dummy values to zero
       scm_dummy_1d(:,:)   = 0.0
       scm_dummy_2d(:,:,:) = 0.0

!set up logicals for model fields
      qneed_ke=qke_m.or.qke_mm.or.qprod_p(17).or.qprod_p(27).or.        &
     &         qprod_p(37).or.qprod_p(47).or.qprod_p(57).or.            &
     &         qprod_p(67).or.qprod_p(77).or.                           &
     &         qcol_uke.or.qcol_vke.or.qcol_wke

      if (qneed_ke) then
        qneed_u=.true.
        qneed_v=.true.
        qneed_w=.true.
      else
        qneed_u=qprod_p(11).or.qprod_p(12).or.qprod_p(13).or.           &
     &          qprod_p(14).or.qprod_p(15).or.qprod_p(16).or.           &
     &          qu_m.or.qu_mm.or.qcol_u.or.qcol_ugz.or.                 &
     &          qcol_ut.or.qcol_uq.or.qcol_uv.or.qcol_uw

        qneed_v=qprod_p(12).or.qprod_p(22).or.qprod_p(23).or.           &
     &          qprod_p(24).or.qprod_p(25).or.qprod_p(26).or.           &
     &          qv_m.or.qv_mm.or.qcol_v.or.qcol_vgz.or.                 &
     &          qcol_vt.or.qcol_vq.or.qcol_vw.or.qcol_uv

        qneed_w=qprod_p(13).or.qprod_p(23).or.qprod_p(33).or.           &
     &          qprod_p(34).or.qprod_p(35).or.qprod_p(36).or.           &
     &          qw_m.or.qw_mm .or.qcol_w.or.qcol_wgz.or.                &
     &          qcol_wt.or.qcol_wq.or.qcol_vw.or.qcol_uw

      endif  ! test on qneed_ke

      qneed_tropopause = qtrop_p .or. qtrop_t                           &
     &  .or. qtrop_z  .or. qtrop_icao

      qneed_t=qprod_p(14).or.qprod_p(24).or.qprod_p(34).or.             &
     &        qprod_p(44).or.qprod_p(45).or.qprod_p(46).or.             &
     &        qt_m.or.qt_mm.or.qcol_wt.or.                              &
     &        qcol_ut.or.qcol_vt.or. qneed_tropopause

      qneed_q=qprod_p(15).or.qprod_p(25).or.qprod_p(35).or.             &
     &        qprod_p(45).or.qprod_p(55).or.qprod_p(56).or.             &
     &        qq_m.or.qq_mm .or.qcol_uq.or.qcol_vq.or.qcol_wq

      qneed_z=qprod_p(16).or.qprod_p(26).or.qprod_p(36).or.             &
     &        qprod_p(46).or.qprod_p(56).or.qprod_p(66).or.             &
     &        qz_m.or.qz_mm .or.qcol_ugz.or.qcol_vgz.or.qcol_wgz



!  Calculate exner at theta points
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            exner_at_theta(i,j,k) = exner_theta_levels(i,j,k)
          ENDDO                 ! i
        ENDDO                   ! j
      ENDDO                     ! k


      IF (qneed_u.or.qneed_v.or.qneed_w.or.                             &
     &  qneed_t.or.qneed_q.or.qdry_mass_m.or.                           &
     &  qam_m1.or.qam_m2.or.qam_m3.or.                                  &
     &  qam_w1.or.qam_w2.or.qam_w3) then

! DEPENDS ON: swap_bounds
        call swap_bounds(u,row_length,rows,model_levels,                &
     &    offx,offy,fld_type_u,.true.)

! DEPENDS ON: swap_bounds
        call swap_bounds(v,row_length,n_rows,model_levels,              &
     &    offx,offy,fld_type_v,.true.)

! DEPENDS ON: u_to_p
        CALL u_to_p(u,row_length,rows,model_levels,                     &
     &    offx,offy,model_domain,at_extremity,field_rho(1,1,1,irho_u))

! DEPENDS ON: v_to_p
        CALL v_to_p(v,row_length,rows,n_rows,model_levels,              &
     &    offx,offy,model_domain,at_extremity,field_rho(1,1,1,irho_v))

! problem of u & v at poles
        If (model_domain  ==  1 ) Then

! DEPENDS ON: polar_vector_wind_n
          Call Polar_vector_wind_n                                      &
     &      (v,                                                         &
     &      sin_theta_longitude,cos_theta_longitude,                    &
     &      row_length,n_rows, model_levels,                            &
     &      mag_vector_np,dir_vector_np,                                &
     &      mag_vector_sp,dir_vector_sp,                                &
     &      offx, offy, global_row_length,                              &
     &      gc_proc_row_group, at_extremity)

          If (at_extremity(PSouth) ) Then
            Do k = 1,model_levels
              Do i = 1, row_length
                field_rho(i,1,k,irho_u) = 0.0
              End Do
            End Do
            Do k = 1,model_levels
              Do i = 1, row_length
                field_rho(i,1,k,irho_v) = mag_vector_sp(k)
              End Do
            End Do
          endif

          If (at_extremity(PNorth) ) Then
            Do k = 1,model_levels
              Do i = 1, row_length
                field_rho(i,rows,k,irho_u) = 0.0
              End Do
            End Do
            Do k = 1,model_levels
              Do i = 1, row_length
                field_rho(i,rows,k,irho_v) = mag_vector_np(k)
              End Do
            End Do
          endif
        endif

!  Calculate delr_rho
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              delr_rho(i,j,k)=r_theta_levels(i,j,k)-                    &
     &          r_theta_levels(i,j,k-1)
            ENDDO               ! i
          ENDDO                 ! j
        ENDDO                   ! k
!----------------------------------------------------------------------
! convert rho to rho dry in the same way as done in dynamics - Flux_rho
!----------------------------------------------------------------------
        k = 1
        Do j = 1, rows
          Do i = 1, row_length
            temp = 1. - q(i,j,k)- qcl(i,j,k) - qcf(i,j,k)
            rho_dry(i,j,k) = rho(i,j,k) * temp
            dry_to_wet(i,j,k)= 1. / temp
          End Do
        End Do

        Do k = 2, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
              weight2 = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
              weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
            temp = ( weight2 *                                          &
     &             (1. - q(i,j,k)- qcl(i,j,k) - qcf(i,j,k) ) +          &
     &               weight1 *                                          &
     &             (1. - q(i,j,k-1)- qcl(i,j,k-1) - qcf(i,j,k-1) ) )    &
     &             / weight3
              dry_to_wet(i,j,k)= 1./temp
              rho_dry(i,j,k) = rho(i,j,k) * temp
            End Do
          End Do
        End Do

        Do k = wet_levels+1,model_levels
          Do j = 1, rows
            Do i = 1, row_length
              rho_dry(i,j,k) = rho(i,j,k)
              dry_to_wet(i,j,k)= 1.
            End Do
          End Do
        End Do
!-----------------------------------------------------------------------
! Mass weighting factor for each location ie rho_dry*r*r*dr/(a*a)
! Note currently using rho_dry not rho_wet
! Note rho_dry includes r*r factor
!-----------------------------------------------------------------------

       Do  k=1,model_levels
         Do j=1,rows
           Do i=1,row_length
              mass_we(i,j,k)=rho_dry(i,j,k)*delr_rho(i,j,k)             &
     &                                /(earth_radius*earth_radius)
           Enddo
         Enddo
       Enddo
      ENDIF     ! test on logicals


! ----------------------------------------------------------------------
! 115 dry mass weighting on model rho levels
! ----------------------------------------------------------------------
      IF(qdry_mass_m) THEN
        count=1
        DO  k=1,model_levels
          If (dry_mass_m_list(k)) then
            Do j=1,rows
              Do i=1,row_length
                dry_mass_m(i,j,count)=mass_we(i,j,k)
              Enddo
            Enddo
            count=count+1
          Endif
        Enddo
      Endif
!-----------------------------------------------------------------------
!  U and V on pressure levels
!-----------------------------------------------------------------------

      IF (qu_p.or.qv_p) THEN
! Need call to swapbounds for exner_rho before interpolation
! DEPENDS ON: swap_bounds
        call swap_bounds(exner_rho_levels,row_length,rows,              &
     &                   model_levels,offx,offy,fld_type_p,.false.)
      ENDIF

      IF(qu_p) THEN

!  Calculate exner at u points. Store in exner_at_u

        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              exner_at_u(i,j,k) = 0.5 *                                 &
     &             (exner_rho_levels(i,j,k) + exner_rho_levels(i+1,j,k))
            ENDDO               ! i
          ENDDO                 ! j
        ENDDO                   ! k

      ENDIF ! on STASHflag

      IF(qv_p) THEN
!  Calculate exner at v points. Store in exner_at_v

        DO k = 1, model_levels
          DO j = 1, n_rows
            DO i = 1, row_length
              exner_at_v(i,j,k) = 0.5 *                                 &
     &             (exner_rho_levels(i,j,k) + exner_rho_levels(i,j+1,k))
            ENDDO               ! i
          ENDDO                 ! j
        ENDDO                   ! k

      ENDIF ! on STASHflag

!   Calculate temperature at theta points
      IF(qt_p .OR. qRH_p.or.qneed_t.or.qteot_m .or. qRH_m) THEN

        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              T(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
              p_at_theta(i,j,k) = p_theta_levels(i,j,k)
              rh(i,j,k) = 0.
            ENDDO               ! i
          ENDDO                 ! j
        ENDDO                   ! k

      ENDIF ! on STASHflag

!  height of model rho levels
!   Remove radius of earth from height field
      IF(qz_p.or.qneed_z) THEN

        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              field_rho(i,j,k,irho_z) = r_rho_levels(i,j,k)             &
     &                                         - earth_radius
            ENDDO               ! i
          ENDDO                 ! j
        ENDDO                   ! k

      ENDIF ! on relevant STASHflags


! extra points done because of LAM domains
      Do j = 1, n_rows+1
        Do i = 1, row_length+1
          pstar_halo(i,j) = p(i,j,1) + g * rho(i,j,1)                   &
     &      * (r_rho_levels(i,j,1) -                                    &
     &      r_theta_levels(i,j,0) ) /                                   &
     &      (r_rho_levels(i,j,1) *                                      &
     &      r_rho_levels(i,j,1) )
        enddo
      enddo

! DEPENDS ON: swap_bounds
      call Swap_Bounds(pstar_halo,                                      &
     &  row_length, rows, 1,                                            &
     &  offx, offy, fld_type_p, .false.)

!CDIR NOUNROLL
      DO j=1,n_rows
        DO i=1,row_length
          pstar_uv(i,j)=(pstar_halo(i,j)+                               &
     &      pstar_halo(i+1,j)+                                          &
     &      pstar_halo(i,j+1)+                                          &
     &      pstar_halo(i+1,j+1)) * 0.25
        ENDDO                   ! i
      ENDDO                     ! j

! ----------------------------------------------------------------------
! STASH items 001,002 : u,v winds on model levels
! ----------------------------------------------------------------------


      IF(qu_m) THEN
        count=1
        DO  k=1,model_levels
          if (u_m_list(k)) then
            do j=1,rows
              DO i=1,row_length
                u_m(i,j,count)=field_rho(i,j,k,irho_u)
              enddo
            enddo
            count=count+1
          endif
        enddo
      Endif
      If(qu_mm) Then
        count=1
        DO  k=1,model_levels
          if (u_mm_list(k)) then
            do j=1,rows
              DO i=1,row_length
                u_mm(i,j,count)=field_rho(i,j,k,irho_u)*mass_we(i,j,k)
              enddo
            enddo
            count=count+1
          endif
        enddo

      ENDIF ! on STASHflag

      IF(qv_m) THEN
        count=1
        DO  k=1,model_levels
          if (v_m_list(k)) then
            do j=1,rows
              DO i=1,row_length
                v_m(i,j,count)=field_rho(i,j,k,irho_v)
              enddo
            enddo
            count=count+1
          endif
        enddo
      Endif
      If(qv_mm) then
        count=1
        DO  k=1,model_levels
          if (v_mm_list(k)) then
            do j=1,rows
              DO i=1,row_length
                v_mm(i,j,count)=field_rho(i,j,k,irho_v)*mass_we(i,j,k)
              enddo
            enddo
            count=count+1
          endif
        enddo
      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH items 003,004,005 : w wind, t and q on model levels
! ----------------------------------------------------------------------
      k=1       !  level 1
      If (qneed_t) then
        Do j = 1, rows
          DO I=1,row_length
              field_rho(i,j,k,irho_t)=t(i,j,k)
          END DO
        END DO
            endif
      If (qneed_q) then
        Do j = 1, rows
          DO I=1,row_length
              field_rho(i,j,k,irho_q)=q(i,j,k)
          END DO
        END DO
            endif
      If (qneed_w) then
        Do j = 1, rows
          DO I=1,row_length
            weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
            weight2 = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
            weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
            ww1 = weight1/weight3
            ww2 = weight2/weight3
            field_rho(i,j,k,irho_w)=ww2 * w(i,j,K)  + ww1 * w(i,j,K-1)
          END DO
        END DO
      endif

      if (qneed_w.or.qneed_t) then
        DO K=2,model_levels
          Do j = 1, rows
            DO I=1,row_length
              weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
              weight2 = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
              weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
              ww1 = weight1/weight3
              ww2 = weight2/weight3
              if (qneed_t) then
                field_rho(i,j,k,irho_t)=ww2 * t(i,j,K) + ww1*t(i,j,K-1)
              endif
              if (qneed_w) then
                field_rho(i,j,k,irho_w)=ww2 * w(i,j,K) + ww1*w(i,j,K-1)
              endif
            END DO
          END DO
        END DO
      endif

      if (qneed_q) then
        DO K=2,wet_levels
          Do j = 1, rows
            DO I=1,row_length
              weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
              weight2 = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
              weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
              ww1 = weight1/weight3
              ww2 = weight2/weight3
              field_rho(i,j,k,irho_q)=ww2 * q(i,j,K)+ ww1 * q(i,j,K-1)
            END DO
          END DO
        END DO
      endif

! ----------------------------------------------------------------------
! 003 w on model rho levels
! ----------------------------------------------------------------------
      IF(qw_m) THEN
        count=1
        DO  k=1,model_levels
          If (w_m_list(k)) then
            DO j=1,rows
              DO i=1,row_length
                w_m(i,j,count)=field_rho(i,j,k,irho_w)
              ENDDO
            ENDDO
            count=count+1
            endif
        ENDDO
      Endif
      If(qt_mm) Then
        count=1
        DO  k=1,model_levels
          If (w_mm_list(k)) then
            DO j=1,rows
              DO i=1,row_length
                w_mm(i,j,count)=field_rho(i,j,k,irho_w)*mass_we(i,j,k)
              ENDDO
            ENDDO
            count=count+1
          Endif
        ENDDO
      ENDIF ! on STASHflag
! ----------------------------------------------------------------------
! 004 T on model rho levels
! ----------------------------------------------------------------------
      IF(qt_m) THEN
        count=1
        DO  k=1,model_levels
          If (t_m_list(k)) then
            DO j=1,rows
              DO i=1,row_length
                t_m(i,j,count)=field_rho(i,j,k,irho_t)
              Enddo
            Enddo
            count=count+1
            endif
        Enddo
      Endif
      If(qt_mm) Then
            count=1
        DO  k=1,model_levels
          If (t_mm_list(k)) then
            DO j=1,rows
              DO i=1,row_length
                t_mm(i,j,count)=field_rho(i,j,k,irho_t)*mass_we(i,j,k)
              Enddo
            Enddo
            count=count+1
          Endif
        Enddo
      ENDIF ! on STASHflag
! ----------------------------------------------------------------------
! 005 q on model rho levels
! ----------------------------------------------------------------------
      IF(qq_m) THEN
        count=1
        DO  k=1,wet_levels
          If (q_m_list(k)) then
            DO j=1,rows
              DO i=1,row_length
                q_m(i,j,count)=field_rho(i,j,k,irho_q)
              Enddo
            Enddo
            count=count+1
          Endif
        Enddo
      Endif
      If(qq_mm) Then
        count=1
        DO  k=1,wet_levels
          If (q_mm_list(k)) then
            DO j=1,rows
              DO i=1,row_length
                q_mm(i,j,count)=field_rho(i,j,k,irho_q)*mass_we(i,j,k)
              Enddo
            Enddo
            count=count+1
          Endif
        Enddo
      ENDIF ! on STASHflag
! ----------------------------------------------------------------------
! 006 z on model rho levels
! ----------------------------------------------------------------------
      IF(qz_m) THEN
        count=1
        DO  k=1,model_levels
          If (z_m_list(k)) then
            DO j=1,rows
              DO i=1,row_length
                z_m(i,j,count)=field_rho(i,j,k,irho_z)
              Enddo
            Enddo
            count=count+1
          Endif
        Enddo
      Endif
      If(qz_mm) Then
        count=1
        DO  k=1,model_levels
          If (z_mm_list(k)) then
            DO j=1,rows
              DO i=1,row_length
                z_mm(i,j,count)=field_rho(i,j,k,irho_z)*mass_we(i,j,k)
              Enddo
            Enddo
            count=count+1
          Endif
        Enddo
      ENDIF ! on STASHflag
! ----------------------------------------------------------------------
! 007 KE on model rho levels
! ----------------------------------------------------------------------
      If (qneed_ke) then
        DO K=1,model_levels
                  Do j = 1, rows
            DO I=1,row_length
              field_rho(i,j,k,irho_ke) = 0.5*(                          &
     &        field_rho(i,j,k,irho_u)*field_rho(i,j,k,irho_u)+          &
     &        field_rho(i,j,k,irho_v)*field_rho(i,j,k,irho_v)+          &
     &        field_rho(i,j,k,irho_w)*field_rho(i,j,k,irho_w))
            END DO
          END DO
        END Do
      Endif

      IF(qke_m) THEN
        count=1
        DO  k=1,model_levels
          If (ke_m_list(k)) then
            DO j=1,rows
              DO i=1,row_length
                ke_m(i,j,count)=field_rho(i,j,k,irho_ke)
                    enddo
                  enddo
                  count=count+1
          Endif
              enddo
      Endif
      If(qke_mm) then
        count=1
        DO  k=1,model_levels
          If (ke_mm_list(k)) then
            DO j=1,rows
              DO i=1,row_length
               ke_mm(i,j,count)=field_rho(i,j,k,irho_ke)*mass_we(i,j,k)
              Enddo
            Enddo
            count=count+1
          Endif
        Enddo
      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH items :products on model levels
! At present all products mass weighted using dry density
! ----------------------------------------------------------------------

      Do iprod=1,nmodel_diags
        Do jprod=iprod,nmodel_diags
          code30=iprod*10+jprod
          IF (qprod_p(code30)) THEN
              count=1
              DO  k=1,model_levels !loop over levels
                if (prod_m_list(k,iprod,jprod)) then
                  Do j = 1, rows
                    Do i = 1, row_length
                    stashwork(si(code30)+ (i-1)+(j-1)*row_length+       &
     &                  (count-1)*rows*row_length)=                     &
     &                  field_rho(i,j,k,iprod)*field_rho(i,j,k,jprod)   &
     &                  *factor*mass_we(i,j,k)
                    enddo
                  enddo
                  count=count+1
                ENDIF
              enddo
            ENDIF
        enddo
      enddo


! ----------------------------------------------------------------------
! STASH item 111 : T at eot
! ----------------------------------------------------------------------
      if (qteot_m) then
        count=1
        DO  k=1,model_levels
          if (teot_m_list(k)) then
            do j=1,rows
              DO i=1,row_length
                teot_m(i,j,count)=t(i,j,k)
              enddo
            enddo
            count=count+1
          endif
        enddo
      endif

! ----------------------------------------------------------------------
! STASH item 181 : T increment
! ----------------------------------------------------------------------
      if (qt_inc) then
! DEPENDS ON: swap_bounds
        call Swap_Bounds(inc_t,                                         &
     &                   row_length, rows, model_levels,                &
     &                   offx, offy, fld_type_p, .true.)
        count=1
        DO  k=1,model_levels
          if (t_inc_list(k)) then
            do j=1,rows
              DO i=1,row_length
                t_inc(i,j,count)=inc_t(i,j,k)
              enddo
            enddo
            count=count+1
          endif
        enddo
      endif

! ----------------------------------------------------------------------
! STASH item 182 : Q increment
! ----------------------------------------------------------------------
      if (qq_inc) then
! DEPENDS ON: swap_bounds
        call Swap_Bounds(inc_q,                                         &
     &                   row_length, rows, wet_levels,                  &
     &                   halo_i, halo_j, fld_type_p, .false.)
        count=1
        DO  k=1,wet_levels
          if (q_inc_list(k)) then
            do j=1,rows
              DO i=1,row_length
                q_inc(i,j,count)=inc_q(i,j,k)
              enddo
            enddo
            count=count+1
          endif
        enddo
      endif

! ----------------------------------------------------------------------
! STASH item 183 : QCL increment
! ----------------------------------------------------------------------
      if (qqcl_inc) then
! DEPENDS ON: swap_bounds
        call Swap_Bounds(inc_qcl,                                       &
     &                   row_length, rows, wet_levels,                  &
     &                   halo_i, halo_j, fld_type_p, .false.)
        count=1
        DO  k=1,wet_levels
          if (qcl_inc_list(k)) then
            do j=1,rows
              DO i=1,row_length
                qcl_inc(i,j,count)=inc_qcl(i,j,k)
              enddo
            enddo
            count=count+1
          endif
        enddo
      endif

! ----------------------------------------------------------------------
! STASH item 184 : QCF increment
! ----------------------------------------------------------------------
      if (qqcf_inc) then
! DEPENDS ON: swap_bounds
        call Swap_Bounds(inc_qcf,                                       &
     &                   row_length, rows, wet_levels,                  &
     &                   halo_i, halo_j, fld_type_p, .false.)
        count=1
        DO  k=1,wet_levels
          if (qcf_inc_list(k)) then
            do j=1,rows
              DO i=1,row_length
                qcf_inc(i,j,count)=inc_qcf(i,j,k)
              enddo
            enddo
            count=count+1
          endif
        enddo
      endif

! ----------------------------------------------------------------------
! STASH item 185 : U increment
! ----------------------------------------------------------------------
      if (qu_inc) then
! DEPENDS ON: swap_bounds
        call Swap_Bounds(inc_u,                                         &
     &                   row_length, rows, model_levels,                &
     &                   offx, offy, fld_type_u, .true.)
        count=1
        DO  k=1,model_levels
          if (u_inc_list(k)) then
            do j=1,rows
              DO i=1,row_length
                u_inc(i,j,count)=inc_u(i,j,k)
              enddo
            enddo
            count=count+1
          endif
        enddo
      endif

! ----------------------------------------------------------------------
! STASH item 186 : V increment
! ----------------------------------------------------------------------
      if (qv_inc) then
! DEPENDS ON: swap_bounds
        call Swap_Bounds(inc_v,                                         &
     &                   row_length, n_rows, model_levels,              &
     &                   offx, offy, fld_type_v, .true.)
        count=1
        DO  k=1,model_levels
          if (v_inc_list(k)) then
            do j=1,n_rows
              DO i=1,row_length
                v_inc(i,j,count)=inc_v(i,j,k)
              enddo
            enddo
            count=count+1
          endif
        enddo
      endif

! ----------------------------------------------------------------------
! STASH item 187 : W increment
! ----------------------------------------------------------------------
      if (qw_inc) then
        count=1
        DO  k=1,model_levels
          if (w_inc_list(k)) then
            do j=1,rows
              DO i=1,row_length
                w_inc(i,j,count)=inc_w(i,j,k)
              enddo
            enddo
            count=count+1
          endif
        enddo
      endif

! ----------------------------------------------------------------------
! STASH item 188 : RHO increment
! ----------------------------------------------------------------------
      if (qrho_inc) then
! DEPENDS ON: swap_bounds
        call Swap_Bounds(inc_rho,                                       &
     &                   row_length, rows, model_levels,                &
     &                   offx, offy, fld_type_p, .true.)
        count=1
        DO  k=1,model_levels
          if (rho_inc_list(k)) then
            do j=1,rows
              DO i=1,row_length
                rho_inc(i,j,count)=inc_rho(i,j,k)/                      &
     &            (r_rho_levels(i,j,k)*r_rho_levels(i,j,k))
              enddo
            enddo
            count=count+1
          endif
        enddo
      endif

! ----------------------------------------------------------------------
! STASH item 171 : T increment**2
! ----------------------------------------------------------------------
      if (qt2_inc) then
! DEPENDS ON: swap_bounds
        call Swap_Bounds(inc_t,                                         &
     &                   row_length, rows, model_levels,                &
     &                   offx, offy, fld_type_p, .true.)
        count=1
        DO  k=1,model_levels
          if (t2_inc_list(k)) then
            do j=1,rows
              DO i=1,row_length
                t2_inc(i,j,count)=inc_t(i,j,k)*inc_t(i,j,k)
              enddo
            enddo
            count=count+1
          endif
        enddo
      endif

! ----------------------------------------------------------------------
! STASH item 172 : Q increment **2
! ----------------------------------------------------------------------
      if (qq2_inc) then
! DEPENDS ON: swap_bounds
        call Swap_Bounds(inc_q,                                         &
     &                   row_length, rows, wet_levels,                  &
     &                   halo_i, halo_j, fld_type_p, .false.)
        count=1
        DO  k=1,wet_levels
          if (q2_inc_list(k)) then
            do j=1,rows
              DO i=1,row_length
                q2_inc(i,j,count)=inc_q(i,j,k)*inc_q(i,j,k)
              enddo
            enddo
            count=count+1
          endif
        enddo
      endif

! ----------------------------------------------------------------------
! STASH item 173 : Qcl increment **2
! ----------------------------------------------------------------------
      if (qqcl2_inc) then
! DEPENDS ON: swap_bounds
        call Swap_Bounds(inc_qcl,                                       &
     &                   row_length, rows, wet_levels,                  &
     &                   halo_i, halo_j, fld_type_p, .false.)
        count=1
        DO  k=1,wet_levels

          if (qcl2_inc_list(k)) then
            do j=1,rows
              DO i=1,row_length
                qcl2_inc(i,j,count)=inc_qcl(i,j,k)*inc_qcl(i,j,k)
              enddo
            enddo
            count=count+1
          endif
        enddo
      endif

! ----------------------------------------------------------------------
! STASH item 174 : Qcf increment **2
! ----------------------------------------------------------------------
      if (qqcf2_inc) then
! DEPENDS ON: swap_bounds
        call Swap_Bounds(inc_qcf,                                       &
     &                   row_length, rows, wet_levels,                  &
     &                   halo_i, halo_j, fld_type_p, .false.)
        count=1
        DO  k=1,wet_levels
          if (qcf2_inc_list(k)) then
            do j=1,rows
              DO i=1,row_length
                qcf2_inc(i,j,count)=inc_qcf(i,j,k)*inc_qcf(i,j,k)
              enddo
            enddo
            count=count+1
          endif
        enddo
      endif
! ----------------------------------------------------------------------
! STASH item 175 : U increment **2
! ----------------------------------------------------------------------
      if (qu2_inc) then
! DEPENDS ON: swap_bounds
        call Swap_Bounds(inc_u,                                         &
     &                   row_length, rows, model_levels,                &
     &                   offx, offy, fld_type_u, .true.)
        count=1
        DO  k=1,model_levels
          if (u2_inc_list(k)) then
            do j=1,rows
              DO i=1,row_length
                u2_inc(i,j,count)=inc_u(i,j,k)*inc_u(i,j,k)
              enddo
            enddo
            count=count+1
          endif
        enddo
      endif

! ----------------------------------------------------------------------
! STASH item 176 : V increment **2
! ----------------------------------------------------------------------
      if (qv2_inc) then
! DEPENDS ON: swap_bounds
        call Swap_Bounds(inc_v,                                         &
     &                   row_length, n_rows, model_levels,              &
     &                   offx, offy, fld_type_v, .true.)
        count=1
        DO  k=1,model_levels
          if (v2_inc_list(k)) then
            do j=1,n_rows
              DO i=1,row_length
                v2_inc(i,j,count)=inc_v(i,j,k)*inc_v(i,j,k)
              enddo
            enddo
            count=count+1
          endif
        enddo
      endif

! ----------------------------------------------------------------------
! STASH item 177 : W increment **2
! ----------------------------------------------------------------------
      if (qw2_inc) then
        count=1
        DO  k=1,model_levels
          if (w2_inc_list(k)) then
            do j=1,rows
              DO i=1,row_length
                w2_inc(i,j,count)=inc_w(i,j,k)*inc_w(i,j,k)
              enddo
            enddo
            count=count+1
          endif
        enddo
      endif

! ----------------------------------------------------------------------
! STASH item 178 : RHO increment **2
! ----------------------------------------------------------------------
      if (qrho2_inc) then
! DEPENDS ON: swap_bounds
        call Swap_Bounds(inc_rho,                                       &
     &                   row_length, rows, model_levels,                &
     &                   offx, offy, fld_type_p, .true.)
        count=1
        DO  k=1,model_levels
          if (rho2_inc_list(k)) then
            do j=1,rows
              DO i=1,row_length
               rho2_inc(i,j,count)=(inc_rho(i,j,k)/                     &
     &            (r_rho_levels(i,j,k)*r_rho_levels(i,j,k)))**2

              enddo
            enddo
            count=count+1
          endif
        enddo
      endif

! ----------------------------------------------------------------------
! STASH items 201 : u wind components on pressure surfaces
! ----------------------------------------------------------------------
      IF(qu_p) THEN

        DO  k=1,u_p_levs

          pressure_pa = u_press(k)*100.0   ! convert to Pascals
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
     &      row_length,rows,n_rows,1,offx,offy,                         &
     &      u_p(1,1,k))

          DO j = 1, n_rows
            DO i = 1, row_length
              if (pstar_uv(i,j) <  u_press(k)*100.) then
                u_p(i,j,k)=0.0
              endif
            enddo
          enddo
        ENDDO  ! k pressure levels loop
      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH items 202 : v wind components on pressure surfaces
! ----------------------------------------------------------------------
      IF(qv_p) THEN
        DO  k=1,v_p_levs

          pressure_pa = v_press(k)*100.0   ! convert to Pascals
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
     &      row_length,n_rows,1,offx,offy,                              &
     &      v_p(1,1,k))
          DO j = 1, n_rows
            DO i = 1, row_length
              if (pstar_uv(i,j) <  v_press(k)*100.) then
                v_p(i,j,k)=0.0
              endif
            enddo
          enddo

        ENDDO  ! k pressure levels loop
      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH items 203 : w wind components on pressure surfaces
! ----------------------------------------------------------------------
      IF(qw_p) THEN
        DO  k=1,w_p_levs

          pressure_pa = w_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
! DEPENDS ON: vert_interp2
          Call vert_interp2 (w(1-offx,1-offy,1)                         &
     &      ,row_length, rows, model_levels                             &
     &      ,pressure_ex                                                &
     &      ,offx, offy, 0,0                                            &
     &      ,exner_at_theta, interp_order                               &
     &      ,work_1 )

! Perform simple horizontal interpolation from 'C' to 'B' grid

! DEPENDS ON: pc_to_pb
          CALL  pC_to_pB(work_1,                                        &
     &      row_length,rows,n_rows,1,offx,offy,                         &
     &      w_p(1,1,k))
          DO j = 1, n_rows
            DO i = 1, row_length
              if (pstar_uv(i,j) <  w_press(k)*100.) then
                w_p(i,j,k)=0.0
              endif
            enddo
          enddo

        ENDDO  ! k pressure levels loop
      ENDIF ! on STASHflag


! ----------------------------------------------------------------------
! STASH item 204: temperature              on pressure surfaces
! ----------------------------------------------------------------------
      IF(qt_p) THEN
        DO k = 1, t_p_levs

          pressure_pa = t_press(k)*100. ! convert to Pascals
! DEPENDS ON: t_vert_interp_to_p
          CALL T_vert_interp_to_p(                                      &
     &      T, theta, row_length, rows                                  &
     &      ,model_levels, pressure_pa, offx, offy ,halo_i,halo_j       &
     &      ,p_theta_levels, Lapse, R, g                                &
     &      ,bl_levels                                                  &
     &      ,exner_theta_levels                                         &
     &      ,r_theta_levels                                             &
     &      ,kappa, p_zero, work_1 )

! DEPENDS ON: pc_to_pb
          CALL  pC_to_pB(work_1,                                        &
     &      row_length,rows,n_rows,1,offx,offy,                         &
     &      t_p(1,1,k))
          DO j = 1, n_rows
            DO i = 1, row_length
              if (pstar_uv(i,j) <  t_press(k)*100.) then
                t_p(i,j,k)=0.0
              endif
            enddo
          enddo


        ENDDO ! over output STASH pressure levels
      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH item 205: Q              on pressure surfaces
! ----------------------------------------------------------------------
      IF(qq_p) THEN
        DO k = 1, q_p_levs

          pressure_pa = q_press(k)*100. ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
! DEPENDS ON: vert_interp2
          CALL vert_interp2(q, row_length, rows, model_levels           &
     &      ,pressure_ex                                                &
     &      ,halo_i, halo_j, 0,0                                        &
     &      ,exner_at_theta, interp_order                               &
     &      ,work_1 )

! DEPENDS ON: pc_to_pb
          CALL  pC_to_pB(work_1,                                        &
     &      row_length,rows,n_rows,1,offx,offy,                         &
     &      q_p(1,1,k))

          DO j = 1, n_rows
            DO i = 1, row_length
              if (pstar_uv(i,j) <  q_press(k)*100.) then
                q_p(i,j,k)=0.0
              endif
            enddo
          enddo

        ENDDO ! over output STASH pressure levels
      ENDIF ! on STASHflag
! ----------------------------------------------------------------------
! STASH item 206: RH              on pressure surfaces
! Calculation of RH also required by  item 113 RH on model levels
! ----------------------------------------------------------------------
      IF(qRH_p .or. qRH_m) THEN
        DO k = 1, wet_levels
!  Find humidity saturation at theta points - store in rh
! DEPENDS ON: qsat
          CALL QSAT(rh(1,1,k),T(1,1,k),p_at_theta(1,1,k),               &
     &      theta_field_size)
!  And convert to relative humidity
          DO j = 1, rows
            DO i = 1, row_length

              rh(i,j,k) = q(i,j,k)/rh(i,j,k)*100.
!  Supersaturation (>100%) can occur with mixed phase scheme but
!  negative humidity is removed from the diagnostic:
              IF(rh(i,j,k) <  0.0) THEN
                rh(i,j,k) = 0.
              ENDIF
            ENDDO               ! i
          ENDDO                 ! j
        ENDDO                   ! k wet_levels

!  Interpolate
        DO k = 1, rh_p_levs

          pressure_pa = RH_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
! DEPENDS ON: vert_interp2
          CALL vert_interp2(                                            &
     &      rh, row_length, rows, model_levels                          &
     &      ,pressure_ex                                                &
     &      ,0,0, 0,0                                                   &
     &      ,exner_at_theta, interp_order                               &
     &      ,work_1 )

! DEPENDS ON: pc_to_pb
          CALL  pC_to_pB(work_1,                                        &
     &      row_length,rows,n_rows,1,offx,offy,                         &
     &      rh_p(1,1,k))
          DO j = 1, n_rows
            DO i = 1, row_length
              if (pstar_uv(i,j) <  rh_press(k)*100.) then
                rh_p(i,j,k)=0.0
              endif
            enddo
          enddo

        ENDDO ! k over output STASH pressure levels

      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH item 207: Z              on pressure surfaces
! ----------------------------------------------------------------------
      IF(qz_p) THEN
        DO k = 1, z_p_levs

         pressure_pa = z_press(k)*100. ! convert to Pascals
! DEPENDS ON: vert_h_onto_p
         CALL vert_h_onto_p(                                            &
     &     field_rho(1,1,1,irho_z), row_length,rows, model_levels       &
     &     ,pressure_pa,r_rho_levels, r_theta_levels                    &
     &, p_theta_levels                                                  &
     &     ,theta, exner_theta_levels                                   &
     &       , exner_rho_levels                                         &
     &     ,R, g, Lapse,bl_levels                                       &
     &     ,offx, offy, halo_i, halo_j                                  &
     &     ,p, interp_order,kappa, p_zero,cp, work_1 )

! DEPENDS ON: pc_to_pb
         CALL  pC_to_pB(work_1,                                         &
     &     row_length,rows,n_rows,1,offx,offy,                          &
     &     z_p(1,1,k))
         DO j = 1, n_rows
           DO i = 1, row_length
             if (pstar_uv(i,j) <  z_press(k)*100.) then
               z_p(i,j,k)=0.0
             endif
           enddo
         enddo

       ENDDO                    ! over output STASH pressure levels
      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH item 208: Omega              on pressure surfaces
! ----------------------------------------------------------------------
      IF(qom_p) THEN
        Do k = 1, model_levels
          Do j = 1-offy, rows+offy
            Do i = 1-offx, row_length+offx
              true_rho_halo(i,j,k) = rho(i,j,k) /                       &
     &          (r_rho_levels(i,j,k) *                                  &
     &          r_rho_levels(i,j,k))
            End Do
          End Do
        End Do

! DEPENDS ON: p_to_t
        call P_TO_T (                                                   &
     &    row_length,rows, halo_i, halo_j,                              &
     &    offx,offy, model_levels-1                                     &
     &    , R_THETA_LEVELS, R_RHO_LEVELS                                &
     &    , true_rho_halo                                               &
     &    , tmp_data_halo                                               &
     &    )


        Do j = 1-offy, rows+offy
          Do i = 1-offx, row_length+offx
            tmp_data_halo(i,j,model_levels)=                            &
     &        true_rho_halo(i,j,model_levels)
          End Do
        End Do


        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              work_2(i,j,k)= -1*g*tmp_data_halo(i,j,k)*w(i,j,k)
            End Do
          End Do
        End Do

        DO  k=1,om_p_levs

          pressure_pa = om_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
! DEPENDS ON: vert_interp2
          CALL vert_interp2(work_2,                                     &
     &      row_length, rows, model_levels                              &
     &      ,pressure_ex                                                &
     &      ,0,0, 0,0                                                   &
     &      ,exner_at_theta, interp_order                               &
     &      ,work_1 )

! Perform simple horizontal interpolation from 'C' to 'B' grid

! DEPENDS ON: pc_to_pb
          CALL  pC_to_pB(work_1,                                        &
     &      row_length,rows,n_rows,1,offx,offy,                         &
     &      om_p(1,1,k))

          DO j = 1, n_rows
            DO i = 1, row_length
              if (pstar_uv(i,j) <  om_press(k)*100.) then
                om_p(i,j,k)=0.0
              endif
            enddo
          enddo
        ENDDO  ! k pressure levels loop
      ENDIF ! on STASHflag


! ----------------------------------------------------------------------
! STASH products items : all pressure surface products
! Loop over the pressure level diags and see if product is required
! If so process
! ----------------------------------------------------------------------
      do iprod=1,npress_diags ! code of 1st field
        do jprod=1,npress_diags ! code of 2nd field
          code30=200+iprod*10+jprod ! code of product
          IF(qprod_p(code30)) THEN   ! product required
            if (iprod  ==  1) then
              field1_p=u_p
            else if (iprod  ==  2) then
              field1_p=v_p
            else if (iprod  ==  3) then
              field1_p=w_p
            else if (iprod  ==  4) then
              field1_p=t_p
            else if (iprod  ==  5) then
              field1_p=q_p
            else if (iprod  ==  6) then
              field1_p=rh_p
            else if (iprod  ==  7) then
              field1_p=z_p
            else if (iprod  ==  8) then
              field1_p=om_p
            endif

            if (jprod  ==  1) then
              field2_p=u_p
            else if (jprod  ==  2) then
              field2_p=v_p
            else if (jprod  ==  3) then
              field2_p=w_p
            else if (jprod  ==  4) then
              field2_p=t_p
            else if (jprod  ==  5) then
              field2_p=q_p
            else if (jprod  ==  6) then
              field2_p=rh_p
            else if (jprod  ==  7) then
              field2_p=z_p
            else if (jprod  ==  8) then
              field2_p=om_p
            endif
            DO  k=1,prod_p_levs(iprod,jprod) !loop over levels
              Do j = 1, n_rows
                Do i = 1, row_length
                  stashwork(si(code30)+                                 &
     &              (i-1)+(j-1)*row_length+(k-1)*n_rows*row_length)=    &
     &              field1_p(i,j,prod_ind(k,iprod,jprod))*              &
     &              field2_p(i,j,prod_ind(k+prod_p_levs(iprod,jprod),   &
     &              iprod,jprod))
                enddo
              enddo
            enddo
          ENDIF                 ! on STASHflag
        enddo
      enddo

! ----------------------------------------------------------------------
! STASH item 301: Heavyside on pressure surfaces
! ----------------------------------------------------------------------
!L The Heavyside function is defined as 1.0 if the pressure level
!L  is above the surface (i.e. pstar) and 0.0 if below. A time mean of
!L  this will give information on the fraction of time a pressure
!L  level is above the land or sea surface.
      if (qheavy_p) then
        do k=1,heavy_p_levs
          Do j = 1, n_rows
            Do i = 1, row_length
              if (pstar_uv(i,j) <  heavy_press(k)*100.) then
                heavy_p(i,j,k)=0.0
              else
                heavy_p(i,j,k)=1.0
              endif
            enddo
          enddo
        enddo
      endif
! ----------------------------------------------------------------------
! STASH item 302: virtual temp on pressure surfaces
! ----------------------------------------------------------------------
      if (qt_p .and. qq_p .and. qtv_p) then
        DO  k=1,prod_p_levs(4,5) !loop over levels shared by t and q
          Do j = 1, n_rows
            Do i = 1, row_length
              tv_p(i,j,k)=t_p(i,j,prod_ind(k,4,5))*                     &
     &          (1+c_virtual*                                           &
     &          q_p(i,j,prod_ind(k+prod_p_levs(4,5),4,5)))
            enddo
          enddo
        enddo
! ----------------------------------------------------------------------
! STASH item 302: virtual temp* omega on pressure surfaces
! ----------------------------------------------------------------------
        if (qom_p .and. qtvom_p) then
          DO  k=1,tvom_p_levs
            Do j = 1, n_rows
              Do i = 1, row_length
                tvom_p(i,j,k)=tv_p(i,j,k)*om_p(i,j,k)
              enddo
            enddo
          enddo
        endif
      endif
! ----------------------------------------------------------------------
! STASH call routine for all column integral fields
! ----------------------------------------------------------------------
      if (qtot_ke.or.qtot_ke_w.or.qtcm.or.qtot_cvt.or.qtot_gr           &
     &  .or.qtcmq.or.qtcqcl.or.qtcqcf                                   &
     &  .or.qtmf_u.or.qtmf_v.or.qtmf_w                                  &
     &  .or.qencorr) then
        grid_factor=1/(earth_radius*earth_radius)
! DEPENDS ON: vert_eng_massq
        call vert_eng_massq(                                            &
     &    halo_i, halo_j, offx, offy                                    &
     &    ,global_row_length, gc_proc_row_group                         &
     &    ,at_extremity, n_proc, n_procx, n_procy                       &
     &    ,neighbour                                                    &
     &    ,mype                                                         &
     &    ,row_length, rows, n_rows                                     &
     &    ,model_domain                                                 &
     &    ,model_levels,wet_levels                                      &
     &    ,r_theta_levels,r_rho_levels                                  &
     &    ,cos_theta_longitude,sin_theta_longitude                      &
     &    ,theta ,u,v,w, rho, q,qcl,qcf                                 &
     &    ,exner_theta_levels                                           &
     &    ,.true.                                                       &
     &    ,rho_dry,dry_to_wet                                           &
     &    ,vert_int_array, flux_array )


      endif
! ----------------------------------------------------------------------
! STASH item 401: Total KE
! ----------------------------------------------------------------------

      if (qtot_ke) then
        Do j = 1, rows
          Do i = 1, row_length
            tot_ke(i,j)=(vert_int_array(i,j,ip_keu)+                    &
     &        vert_int_array(i,j,ip_kev))*grid_factor
          enddo
        enddo
      endif
! ----------------------------------------------------------------------
! STASH item 402: Total KE with W
! ----------------------------------------------------------------------
      if (qtot_ke_w) then
        Do j = 1, rows
          Do i = 1, row_length
            tot_ke_w(i,j)=(vert_int_array(i,j,ip_keu)+                  &
     &        vert_int_array(i,j,ip_kev)+                               &
     &        vert_int_array(i,j,ip_kew))*grid_factor
          enddo
        enddo
      endif
! ----------------------------------------------------------------------
! STASH item 403: Total column dry mass
! ----------------------------------------------------------------------
      if (qtcm) then
        Do j = 1, rows
          Do i = 1, row_length
            tcm(i,j)=vert_int_array(i,j,ip_dry_mass)*grid_factor
          enddo
        enddo
      endif
! ----------------------------------------------------------------------
! STASH item 404: Total column wet mass
! ----------------------------------------------------------------------
      if (qtcmq) then
        Do j = 1, rows
          Do i = 1, row_length
            tcmq(i,j)=vert_int_array(i,j,ip_wet_mass)*grid_factor
          enddo
        enddo
      endif
! ----------------------------------------------------------------------
! STASH item 405:Total column qcl
! ----------------------------------------------------------------------
      if (qtcqcl) then
        Do j = 1, rows
          Do i = 1, row_length
            tcqcl(i,j)=vert_int_array(i,j,ip_qcl)*grid_factor
          enddo
        enddo
      endif
! ----------------------------------------------------------------------
! STASH item 406: Total column qcf
! ----------------------------------------------------------------------
      if (qtcqcf) then
        Do j = 1, rows
          Do i = 1, row_length
            tcqcf(i,j)=vert_int_array(i,j,ip_qcf)*grid_factor
          enddo
        enddo
      endif
! ----------------------------------------------------------------------
! STASH item 407: Total moisture flux U
! ----------------------------------------------------------------------
      if (qtmf_u) then
        Do j = 1, rows
          Do i = 1, row_length
            tmf_u(i,j)=flux_array(i,j,ip_qu)*grid_factor
          enddo
        enddo
      endif
! ----------------------------------------------------------------------
! STASH item 408: Total moisture flux V
! ----------------------------------------------------------------------
      if (qtmf_v) then
        Do j = 1, rows
          Do i = 1, row_length
            tmf_v(i,j)=flux_array(i,j,ip_qv)*grid_factor
          enddo
        enddo
      endif
! ----------------------------------------------------------------------
! STASH item 409: Total moisture flux W
! ----------------------------------------------------------------------
      if (qtmf_w) then
        Do j = 1, rows
          Do i = 1, row_length
            tmf_w(i,j)=flux_array(i,j,ip_qw)*grid_factor
          enddo
        enddo
      endif
! ----------------------------------------------------------------------
! STASH item 410: Mountain torque
! ----------------------------------------------------------------------

      if (qmtorque) then
        Do j = 1, rows
          Do i = 1, row_length
            mtorque(i,j)= 0.5 * cos_theta_latitude(i,j)                 &
     &                        * ( pstar_halo(i,j) + pstar_halo(i+1,j)  )&
     &                         * (  r_theta_levels(i+1,j,0)             &
     &                            - r_theta_levels( i ,j,0)            )&
     &                         / ( delta_lambda                        )
          enddo
        enddo
      endif

! ----------------------------------------------------------------------
! STASH item 440 and 441: East-West and North-South components
!                         of the surface pressure drag
! ----------------------------------------------------------------------

      if (qm_spd_ew) then
        Do j = 1, rows
          Do i = 1, row_length
            m_spd_ew(i,j)= 0.5 * ( pstar_halo(i,j) + pstar_halo(i+1,j) )&
     &                         * (  r_theta_levels(i+1,j,0)             &
     &                            - r_theta_levels( i ,j,0)            )&
     &                         / ( delta_lambda * earth_radius         )
          enddo
        enddo
      endif

      if (qm_spd_ns) then
        Do j = 1, n_rows
          Do i = 1, row_length
            m_spd_ns(i,j)= 0.5 * ( pstar_halo(i,j) + pstar_halo(i,j+1) )&
     &                         * (  r_theta_levels(i,j+1,0)             &
     &                            - r_theta_levels(i, j ,0)            )&
     &                         / ( delta_phi    * earth_radius         )
          enddo
        enddo
      endif
! ----------------------------------------------------------------------
! STASH items 411-416:Angular momentum
! ----------------------------------------------------------------------
      if (qam_m1.or.qam_m2.or.qam_m3.or.                                &
     &  qam_w1.or.qam_w2.or.qam_w3) then
! constants
        factor=1.e-24

! calculate longitude  & latitude
        do j=1,rows
          gj = datastart(2) + j - 1
          do i=1,row_length
            gi = datastart(1) + i - 1
            cossq(i,j)=cos_theta_latitude(i,j)*cos_theta_latitude(i,j)
          enddo
        enddo
        if (qam_m1.or.qam_m2.or.qam_w1.or.qam_w2) then
          do j=1,rows
            do i=1,row_length
              spclcp(i,j)=sin_theta_latitude(i,j)*                      &
     &          cos_theta_longitude(i,j)*cos_theta_latitude(i,j)
              spslcp(i,j)=sin_theta_latitude(i,j)*                      &
     &          sin_theta_longitude(i,j)*cos_theta_latitude(i,j)
              clcp(i,j)=cos_theta_longitude(i,j)*cos_theta_latitude(i,j)
              slcp(i,j)=sin_theta_longitude(i,j)*cos_theta_latitude(i,j)
            enddo
          enddo
        endif
        do j=1,rows
          do i=1,row_length
            if (qam_m1) am_m1(i,j)=0.0
            if (qam_m2) am_m2(i,j)=0.0
            if (qam_m3) am_m3(i,j)=0.0
            if (qam_w1) am_w1(i,j)=0.0
            if (qam_w2) am_w2(i,j)=0.0
            if (qam_w3) am_w3(i,j)=0.0
          enddo
        enddo
        do k=1,model_levels
          do j=1,rows
            do i=1,row_length
              r3dr=rho_dry(i,j,k)*r_rho_levels(i,j,k)                   &
     &          *delr_rho(i,j,k)*factor
              rocos=omega*r_rho_levels(i,j,k)*cos_theta_latitude(i,j)

              if (qam_w1) am_w1(i,j)=am_w1(i,j)                         &
     &          +(field_rho(i,j,k,irho_u)*spclcp(i,j)                   &
     &           -field_rho(i,j,k,irho_v)*slcp(i,j))*r3dr
              if (qam_w2) am_w2(i,j)=am_w2(i,j)                         &
     &          +(field_rho(i,j,k,irho_u)*spslcp(i,j)                   &
     &           +field_rho(i,j,k,irho_v)*clcp(i,j))*r3dr
              if (qam_w3) am_w3(i,j)=                                   &
     &          am_w3(i,j)-field_rho(i,j,k,irho_u)*cossq(i,j)*r3dr

              if (qam_m1) am_m1(i,j)=am_m1(i,j)+rocos*spclcp(i,j)*r3dr
              if (qam_m2) am_m2(i,j)=am_m2(i,j)+rocos*spslcp(i,j)*r3dr
              if (qam_m3) am_m3(i,j)=am_m3(i,j)-rocos*cossq(i,j)*r3dr
            enddo
          enddo
        enddo
      endif

! ----------------------------------------------------------------------
! STASH items 417-418:Pstar
! ----------------------------------------------------------------------
      if (qpstar) then
        do j=1,rows
          do i=1,row_length
            o_pstar(i,j)=pstar(i,j)
          enddo
        enddo
      endif

      if (qpstar_uv) then
        do j=1,n_rows
          do i=1,row_length
            o_pstar_uv(i,j)=pstar_uv(i,j)
          enddo
        enddo
      endif
! ----------------------------------------------------------------------
! STASH item 419: Energy Correction
! ----------------------------------------------------------------------
      if (qencorr) then
        Do j = 1, rows
          Do i = 1, row_length
            encorr(i,j)=(cp-r)*energy_corr_now*                         &
     &        vert_int_array(i,j,ip_dry_mass)*grid_factor
          enddo
        enddo
      endif
! ----------------------------------------------------------------------
! STASH item 420:Total column cvT
! ----------------------------------------------------------------------
      if (qtot_cvt) then
        Do j = 1, rows
          Do i = 1, row_length
            tot_cvt(i,j)=vert_int_array(i,j,ip_cvt)*grid_factor
          enddo
        enddo

      endif
! ----------------------------------------------------------------------
! STASH item 421:Total column gr
! ----------------------------------------------------------------------
      if (qtot_gr) then
        Do j = 1, rows
          Do i = 1, row_length
            tot_gr(i,j)=vert_int_array(i,j,ip_gr)*grid_factor
          enddo
        enddo
      endif
! ----------------------------------------------------------------------
! STASH item 422: Column integral ugz dry mass weighted
! ----------------------------------------------------------------------
      If (qcol_ugz) then
        Do j = 1, rows
          Do i = 1, row_length
            col_ugz(i,j)= 0.0
          Enddo
        Enddo
        Do k=1,model_levels
          Do j = 1, rows
            Do i = 1, row_length
              col_ugz(i,j)= col_ugz(i,j)+field_rho(i,j,k,irho_u)*g*     &
     &                     field_rho(i,j,k,irho_z)*mass_we(i,j,k)
            Enddo
          Enddo
        Enddo
      Endif
! ----------------------------------------------------------------------
! STASH item 423: Column integral vgz dry mass weighted
! ----------------------------------------------------------------------
      If (qcol_vgz) then
        Do j = 1, rows
          Do i = 1, row_length
            col_vgz(i,j)= 0.0
          Enddo
        Enddo
        Do k=1,model_levels
          Do j = 1, rows
            Do i = 1, row_length
              col_vgz(i,j)= col_vgz(i,j)+field_rho(i,j,k,irho_v)*g*     &
     &                     field_rho(i,j,k,irho_z)*mass_we(i,j,k)
            Enddo
          Enddo
        Enddo
      Endif
! ----------------------------------------------------------------------
! STASH item 424: Column integral wgz dry mass weighted
! ----------------------------------------------------------------------
      If (qcol_wgz) then
        Do j = 1, rows
          Do i = 1, row_length
            col_wgz(i,j)= 0.0
          Enddo
        Enddo
        Do k=1,model_levels
          Do j = 1, rows
            Do i = 1, row_length
              col_wgz(i,j)= col_wgz(i,j)+field_rho(i,j,k,irho_w)*g*     &
     &                     field_rho(i,j,k,irho_z)*mass_we(i,j,k)
            Enddo
          Enddo
        Enddo
      Endif
! ----------------------------------------------------------------------
! STASH item 425: Column integral uT dry mass weighted
! ----------------------------------------------------------------------
      If (qcol_ut) then
        Do j = 1, rows
          Do i = 1, row_length
            col_ut(i,j)= 0.0
          Enddo
        Enddo
        Do k=1,model_levels
          Do j = 1, rows
            Do i = 1, row_length
              col_ut(i,j)=col_ut(i,j)+ field_rho(i,j,k,irho_u)*         &
     &                     field_rho(i,j,k,irho_t)*mass_we(i,j,k)
            Enddo
          Enddo
        Enddo
      Endif
! ----------------------------------------------------------------------
! STASH item 426: Column integral vt dry mass weighted
! ----------------------------------------------------------------------
      If (qcol_vt) then
        Do j = 1, rows
          Do i = 1, row_length
            col_vt(i,j)= 0.0
          Enddo
        Enddo
        Do k=1,model_levels
          Do j = 1, rows
            Do i = 1, row_length
              col_vt(i,j)=col_vt(i,j)+ field_rho(i,j,k,irho_v)*         &
     &                     field_rho(i,j,k,irho_t)*mass_we(i,j,k)
            Enddo
          Enddo
        Enddo
      Endif
! ----------------------------------------------------------------------
! STASH item 427: Column integral wt dry mass weighted
! ----------------------------------------------------------------------
      If (qcol_wt) then
        Do j = 1, rows
          Do i = 1, row_length
            col_wt(i,j)= 0.0
          Enddo
        Enddo
        Do k=1,model_levels
          Do j = 1, rows
            Do i = 1, row_length
              col_wt(i,j)=col_wt(i,j)+ field_rho(i,j,k,irho_w)*         &
     &                     field_rho(i,j,k,irho_t)*mass_we(i,j,k)
            Enddo
          Enddo
        Enddo
      Endif
! ----------------------------------------------------------------------
! STASH item 428: Column integral uq dry mass weighted
! ----------------------------------------------------------------------
      If (qcol_uq) then
        Do j = 1, rows
          Do i = 1, row_length
            col_uq(i,j)= 0.0
          Enddo
        Enddo
        Do k=1,model_levels
          Do j = 1, rows
            Do i = 1, row_length
              col_uq(i,j)= col_uq(i,j)+field_rho(i,j,k,irho_u)*         &
     &                     field_rho(i,j,k,irho_q)*mass_we(i,j,k)
            Enddo
          Enddo
        Enddo
      Endif
! ----------------------------------------------------------------------
! STASH item 429: Column integral vq dry mass weighted
! ----------------------------------------------------------------------
      If (qcol_vq) then
        Do j = 1, rows
          Do i = 1, row_length
            col_vq(i,j)= 0.0
          Enddo
        Enddo
        Do k=1,model_levels
          Do j = 1, rows
            Do i = 1, row_length
              col_vq(i,j)= col_vq(i,j)+field_rho(i,j,k,irho_v)*         &
     &                     field_rho(i,j,k,irho_q)*mass_we(i,j,k)
            Enddo
          Enddo
        Enddo
      Endif
! ----------------------------------------------------------------------
! STASH item 430: Column integral wq dry mass weighted
! ----------------------------------------------------------------------
      If (qcol_wq) then
        Do j = 1, rows
          Do i = 1, row_length
            col_wq(i,j)= 0.0
          Enddo
        Enddo
        Do k=1,model_levels
          Do j = 1, rows
            Do i = 1, row_length
              col_wq(i,j)= col_wq(i,j)+ field_rho(i,j,k,irho_w)*        &
     &                     field_rho(i,j,k,irho_q)*mass_we(i,j,k)
            Enddo
          Enddo
        Enddo
      Endif
! ----------------------------------------------------------------------
! STASH item 431: Column integral uv dry mass weighted
! ----------------------------------------------------------------------
      If (qcol_uv) then
        Do j = 1, rows
          Do i = 1, row_length
            col_uv(i,j)= 0.0
          Enddo
        Enddo
        Do k=1,model_levels
          Do j = 1, rows
            Do i = 1, row_length
              col_uv(i,j)= col_uv(i,j)+ field_rho(i,j,k,irho_u)*        &
     &                     field_rho(i,j,k,irho_v)*mass_we(i,j,k)
            Enddo
          Enddo
        Enddo
      Endif
! ----------------------------------------------------------------------
! STASH item 432: Column integral uw dry mass weighted
! ----------------------------------------------------------------------
      If (qcol_uw) then
        Do j = 1, rows
          Do i = 1, row_length
            col_uw(i,j)= 0.0
          Enddo
        Enddo
        Do k=1,model_levels
          Do j = 1, rows
            Do i = 1, row_length
              col_uw(i,j)= col_uw(i,j)+ field_rho(i,j,k,irho_u)*        &
     &                     field_rho(i,j,k,irho_w)*mass_we(i,j,k)
            Enddo
          Enddo
        Enddo
      Endif
! ----------------------------------------------------------------------
! STASH item 433: Column integral vw dry mass weighted
! ----------------------------------------------------------------------
      If (qcol_vw) then
        Do j = 1, rows
          Do i = 1, row_length
            col_vw(i,j)= 0.0
          Enddo
        Enddo
        Do k=1,model_levels
          Do j = 1, rows
            Do i = 1, row_length
              col_vw(i,j)=col_vw(i,j)+ field_rho(i,j,k,irho_v)*         &
     &                     field_rho(i,j,k,irho_w)*mass_we(i,j,k)
            Enddo
          Enddo
        Enddo
      Endif
! ----------------------------------------------------------------------
! STASH item 434: Column integral uke dry mass weighted
! ----------------------------------------------------------------------
      If (qcol_uke) then
        Do j = 1, rows
          Do i = 1, row_length
            col_uke(i,j)= 0.0
          Enddo
        Enddo
        Do k=1,model_levels
          Do j = 1, rows
            Do i = 1, row_length
              col_uke(i,j)=col_uke(i,j)+ field_rho(i,j,k,irho_u)*       &
     &                     field_rho(i,j,k,irho_ke)*mass_we(i,j,k)
            Enddo
          Enddo
        Enddo
      Endif
! ----------------------------------------------------------------------
! STASH item 435: Column integral vke dry mass weighted
! ----------------------------------------------------------------------
      If (qcol_vke) then
        Do j = 1, rows
          Do i = 1, row_length
            col_vke(i,j)= 0.0
          Enddo
        Enddo
        Do k=1,model_levels
          Do j = 1, rows
            Do i = 1, row_length
              col_vke(i,j)= col_vke(i,j)+field_rho(i,j,k,irho_v)*       &
     &                     field_rho(i,j,k,irho_ke)*mass_we(i,j,k)
            Enddo
          Enddo
        Enddo
      Endif
! ----------------------------------------------------------------------
! STASH item 436: Column integral wke dry mass weighted
! ----------------------------------------------------------------------
      If (qcol_wke) then
        Do j = 1, rows
          Do i = 1, row_length
            col_wke(i,j)= 0.0
          Enddo
        Enddo
        Do k=1,model_levels
          Do j = 1, rows
            Do i = 1, row_length
              col_wke(i,j)=col_wke(i,j)+field_rho(i,j,k,irho_w)*        &
     &                     field_rho(i,j,k,irho_ke)*mass_we(i,j,k)
            Enddo
          Enddo
        Enddo
      Endif
! ----------------------------------------------------------------------
! STASH item 437: Column integral u dry mass weighted
! ----------------------------------------------------------------------
      If (qcol_u) then
        Do j = 1, rows
          Do i = 1, row_length
            col_u(i,j)= 0.0
          Enddo
        Enddo
        Do k=1,model_levels
          Do j = 1, rows
            Do i = 1, row_length
              col_u(i,j)=col_u(i,j)+field_rho(i,j,k,irho_u)             &
     &                                         *mass_we(i,j,k)
            Enddo
          Enddo
        Enddo
      Endif
! ----------------------------------------------------------------------
! STASH item 438: Column integral v dry mass weighted
! ----------------------------------------------------------------------
      If (qcol_v) then
        Do j = 1, rows
          Do i = 1, row_length
            col_v(i,j)= 0.0
          Enddo
        Enddo
        Do k=1,model_levels
          Do j = 1, rows
            Do i = 1, row_length
              col_v(i,j)=col_v(i,j)+field_rho(i,j,k,irho_v)             &
     &                                         *mass_we(i,j,k)
            Enddo
          Enddo
        Enddo
      Endif
! ----------------------------------------------------------------------
! STASH item 439: Column integral w dry mass weighted
! ----------------------------------------------------------------------
      If (qcol_w) then
        Do j = 1, rows
          Do i = 1, row_length
            col_w(i,j)= 0.0
          Enddo
        Enddo
        Do k=1,model_levels
          Do j = 1, rows
            Do i = 1, row_length
              col_w(i,j)= col_w(i,j)+ field_rho(i,j,k,irho_w)           &
     &                                         *mass_we(i,j,k)
            Enddo
          Enddo
        Enddo
      Endif
! ----------------------------------------------------------------------
!  Stash items 112 : wbig w > wbig on model levels
! ----------------------------------------------------------------------
      IF(qwbig_m) THEN
        wbig=1.0     ! test value
        count=1
        DO k=1,model_levels     ! loop over levels
          IF (wbig_m_list(k)) THEN
             DO j = 1,rows
               DO i = 1,row_length
                 if (w(i,j,k) >  wbig) then
                    wbig_m(i,j,count) = 1.0
                 else
                    wbig_m(i,j,count) = 0.0
                 endif
               ENDDO
             ENDDO
             count=count+1
          ENDIF
        ENDDO
      ENDIF
! ----------------------------------------------------------------------
!  Stash items 114 : wbig2 w > wbig on model levels
! ----------------------------------------------------------------------
      IF(qwbig2_m) THEN
        wbig2=0.1     ! test value
        count=1
        DO k=1,model_levels     ! loop over levels
          IF (wbig2_m_list(k)) THEN
             DO j = 1,rows
               DO i = 1,row_length
                 if (w(i,j,k) >  wbig2) then
                    wbig2_m(i,j,count) = 1.0
                 else
                    wbig2_m(i,j,count) = 0.0
                 endif
               ENDDO
             ENDDO
             count=count+1
          ENDIF
        ENDDO
      ENDIF
! ----------------------------------------------------------------------
!  Stash items 113 : RH on model levels - see earlier calculation
! ----------------------------------------------------------------------
      IF(qRH_m) THEN
        count=1
        DO k=1,wet_levels     ! loop over levels
          IF (RH_m_list(k)) THEN
             DO j = 1,rows
               DO i = 1,row_length
                   RH_m(i,j,count) = RH(i,j,k)
               ENDDO
             ENDDO
             count=count+1
          ENDIF
        ENDDO
      ENDIF
! ----------------------------------------------------------------------
! STASH item 451,452,453,454: Tropopause diagnostics
! ----------------------------------------------------------------------
      IF (qneed_tropopause) THEN
!
!       The (thermal) tropopause is required.
!       This code should be maintained in parallel with
!       subroutines RAD_CTL and O3_to_3D.
!
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              T_filled_halo(i,j,k)=T(i,j,k)
            End Do
          End Do
        End Do
!
! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &    T_filled_halo, row_length, rows,                              &
     &    model_levels, offx, offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        Call Fill_external_halos(T_filled_halo,row_length,rows,         &
     &    model_levels,offx,offy)

! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &    exner_rho_levels, row_length, rows,                           &
     &    model_levels, offx, offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        Call Fill_external_halos(exner_rho_levels,row_length,rows,      &
     &    model_levels,offx,offy)
!
! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &    exner_theta_levels, row_length, rows,                         &
     &    model_levels, offx, offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        Call Fill_external_halos(exner_theta_levels,row_length,rows,    &
     &    model_levels,offx,offy)


!
! DEPENDS ON: tropin
        Call tropin (T_filled_halo,                                     &
     &                 exner_rho_levels, exner_theta_levels,            &
     &                 row_length, rows, model_levels, offx, offy,      &
     &                 at_extremity,scm_dummy_1d,scm_dummy_2d,          &
     &                 min_trop_level, max_trop_level, trindx )
!
!           Set the tropopause height at the grid-point.
!           To calculate at the next highest theta level, use
!              r_theta_levels(i, j, trindx(i, j)).
!           To calculate at the model layer boundary, use
!              r_rho_levels(i, j, trindx(i, j)).
!           Note that as TRINDX
!           is counted upward from the surface, while the first rho
!           level is omitted from the physics, no offset to the final
!           index of r_rho_levels is required.
!           (See Deck O3EXP1.301--304)
!           Values of r are measured from the centre of the Earth.
!           Hence subtraction of the surface r value must be performed.
        IF (QTROP_Z) THEN
          Do j=1, rows
            Do i=1, row_length
              trop_z(i,j) =   r_rho_levels(i, j, trindx(i, j))          &
     &                      - r_theta_levels(i, j, 0)
            Enddo
          Enddo
        ENDIF
!
!           Set the tropopause temperature at the grid-point.
!           To calculate at the next highest theta level, set
!              trop_t(i,j)=T_filled_halo(i, j, trindx(i, j))
!           Use linear interpolation in r
!              (consistent with Energy diagnostics).
!           Note that as TRINDX
!           is counted upward from the surface, while the first rho
!           level is omitted from the physics, no offset to the final
!           index of r_rho_levels is required.
!           (See Deck O3EXP1.301--304)
!
        IF (QTROP_T) THEN
          Do j=1, rows
            Do i=1, row_length
                k = trindx(i,j)
                weight1=r_theta_levels(i,j,k)-r_rho_levels(i,j,k)
                weight2=r_rho_levels(i,j,k)-r_theta_levels(i,j,k-1)
                weight3=r_theta_levels(i,j,k)-r_theta_levels(i,j,k-1)

                ww1 = weight1/weight3
                ww2 = weight2/weight3

                trop_t(i,j) =   ww2 * T_filled_halo(i,j,k)              &
     &                        + ww1 * T_filled_halo(i,j,k-1)

            Enddo
          Enddo
        ENDIF

!           Set the tropopause pressure at the grid-point.
!           To calculate at the next highest theta level, use
!              p_theta_levels(i, j, trindx(i, j)).
!           To calculate at the model layer boundary, use
!              p(i, j, trindx(i, j)).
!           Note that as TRINDX
!           is counted upward from the surface, while the first rho
!           level is omitted from the physics, no offset to the final
!           index of r_rho_levels is required.
!           (See Deck O3EXP1.301--304)
!
        IF (QTROP_P) THEN
          Do j=1, rows
            Do i=1, row_length
              trop_p(i,j)=p(i, j, trindx(i, j))
            Enddo
          Enddo
        ENDIF

!       Interim initialization of qmodel_half_height
        QMODEL_HALF_HEIGHT = .true.
!
        IF(QMODEL_HALF_HEIGHT) THEN
!-----------------------------------------------------------------------
!      Call tropopause program with quality control
!      Needs CALL to TROP here
!-----------------------------------------------------------------------
!         Calculates ICAO heights from Trop pressure field
!-----------------------------------------------------------------------
          IF(QTROP_ICAO) THEN
            Do j = 1, rows
              Do i = 1, row_length
                trop_icao(i,j)=13840. ! metres
              enddo
            enddo
! Interim values in the absence of subroutine ICAO_HT
          ENDIF
!
        ELSE
          WRITE(6,444)
          WRITE(6,444)
 444      FORMAT(' Subroutine TROP not called No MODEL_HALF_HEIGHTS')
        ENDIF
      ENDIF
! ---------------------------------------------------------------------
! Check error condition
      IF(ErrorStatus >  0) THEN
! DEPENDS ON: ereport
         CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

      RETURN
      END SUBROUTINE Eot_diag
#endif
