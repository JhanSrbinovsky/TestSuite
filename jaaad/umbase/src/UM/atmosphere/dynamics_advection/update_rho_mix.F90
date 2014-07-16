#if defined(A12_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  subroutine update_rho_mix
      subroutine update_rho_mix(                                        &
     &                      timestep, NumCycles, CycleNo, rows, n_rows, &
     &                      row_length, model_levels, wet_model_levels, &
     &                      model_domain, first_constant_r_rho_level,   &
     &                      alpha_1, alpha_2, rims_to_do,               &
     &                      halo_i, halo_j, nproc, gc_proc_row_group,   &
     &                      L_regular, at_extremity, global_row_length, &
     &                      off_x, off_y, halo_i_wind, halo_j_wind,     &
     &                      cos_v_latitude, delta_lambda, delta_phi,    &
     &                      lambda_p, phi_p, lambda_u, phi_v,           &
     &                      recip_dlambda_u, recip_dphi_p,              &
     &                      recip_dphi_v,                               &
     &                      wt_lambda_u, wt_phi_v,                      &
     &                      R_u, R_v, R_w,                              &
     &                      u, v, w, rho, rho_np1,                      &
     &                      mix_v, mix_cl, mix_cf,                      &
     &                      mix_cf2, mix_rain, mix_graup,               &
     &                      mix_v_star, mix_cl_star, mix_cf_star,       &
     &                      mix_cf2_star, mix_rain_star, mix_graup_star,&
     &                      mix_v_np1, mix_cl_np1, mix_cf_np1,          &
     &                      mix_cf2_np1, mix_rain_np1, mix_graup_np1,   &
     &                      L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,      &
     &                      r_theta_levels, r_rho_levels,               &
     &                      eta_theta_levels, eta_rho_levels,           &
     &                      FV_sec_theta_latitude,                      &
     &                      wet_to_dry_n, wet_to_dry_np1,               &
     &                      rho_n, inc_rho, L_do_increment, L_new_tdisc &
     &                      )

! Purpose:
!         subroutine to update rho, to tidy ATM_STEP
!
! Method:
!          Is described in ;
!
!
! Original Programmer: Andrew J. Malcolm
! Current code owner: Andrew J. Malcolm
!
! History:
! Version   Date       Comment
! ----     -------     -------
! 5.4     15/10/00     This deck introduced (based on update_rho)
!                                                         Andy Malcolm
! 5.5      14/02/03  improve tracer conservation properties A.Malcolm
!  5.5  26/03/03  put swapbounds calls before fill_halos     A.Malcolm
!   5.5     25/02/03  extra moisture variable changes    Andy Malcolm
! 6.2     25/11/05   Add code for cycling semi-Lagrangian scheme
!                                                        M. Diamantakis
!  6.2  25/12/05  Variable resolution changes            Yongming Tang
!  6.2  25/12/05  Change solver domain size              Terry  Davies
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
              ! model dimensions
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, n_rows                                                          &
                         ! number of v rows.
     &, model_levels                                                    &
                         ! number of model levels
     &, wet_model_levels                                                &
                         ! number of model levels where moisture
                         ! variables are held.
     &, halo_i                                                          &
                     ! Size of halo in i.
     &, halo_j                                                          &
                     ! Size of halo in j.
     &, off_x                                                           &
     &, off_y                                                           &
     &, halo_i_wind                                                     &
                     ! Size of halo in i for wind fields.
     &, halo_j_wind                                                     &
                     ! Size of halo in j for wind fields.
     &, nproc                                                           &
                    ! Total number of processors
     &, gc_proc_row_group                                               &
                          ! Group id for processors on the same row
     &, global_row_length                                               &
     &, rims_to_do                                                      &
                        ! rim size of lbc weights = 1
     &, NumCycles                                                       &
     &, CycleNo

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_regular   ! false if variable resolution


      Integer                                                           &
     &  first_constant_r_rho_level ! first rho level on which r
                                   ! is constant.

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

      Real                                                              &
     &  timestep                                                        &
     &, alpha_1                                                         &
     &, alpha_2

      Real                                                              &
     &  u(1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)   &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y, model_levels) &
     &, w(1-off_x:row_length+off_x, 1-off_y:rows+off_y, 0:model_levels) &
     &, rho_n (1-off_x:row_length+off_x,                                &
     &         1-off_y:rows+off_y, model_levels)                        &
     &, mix_v  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          wet_model_levels)                                       &
     &, mix_cl (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          wet_model_levels)                                       &
     &, mix_cf (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          wet_model_levels)                                       &
     &, mix_cf2 (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &       wet_model_levels)                                          &
     &, mix_rain (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,     &
     &       wet_model_levels)                                          &
     &, mix_graup (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,    &
     &       wet_model_levels)                                          &
     &, mix_v_star  (1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &               wet_model_levels)                                  &
     &, mix_cl_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &               wet_model_levels)                                  &
     &, mix_cf_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &               wet_model_levels)                                  &
     &, mix_cf2_star (1-off_x:row_length+off_x,                         &
     &               1-off_y:rows+off_y,wet_model_levels)               &
     &, mix_rain_star (1-off_x:row_length+off_x,                        &
     &               1-off_y:rows+off_y,wet_model_levels)               &
     &, mix_graup_star (1-off_x:row_length+off_x,                       &
     &               1-off_y:rows+off_y,wet_model_levels)               &
     &, mix_v_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,        &
     &               wet_model_levels)                                  &
     &, mix_cl_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &               wet_model_levels)                                  &
     &, mix_cf_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &               wet_model_levels)                                  &
     &, mix_cf2_np1 (1-off_x:row_length+off_x,                          &
     &               1-off_y:rows+off_y,wet_model_levels)               &
     &, mix_rain_np1 (1-off_x:row_length+off_x,                         &
     &               1-off_y:rows+off_y,wet_model_levels)               &
     &, mix_graup_np1 (1-off_x:row_length+off_x,                        &
     &               1-off_y:rows+off_y,wet_model_levels)

      Real, Intent(Out) ::                                              &
     &  inc_rho(1-off_x:row_length+off_x,                               &
     &             1-off_y:rows+off_y, model_levels)


      Real                                                              &
     &  R_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &        model_levels)                                             &
     &, R_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,             &
     &        model_levels)                                             &
     &, R_w(row_length, rows, model_levels)

      Real                                                              &
           ! trigonometric arrays.
     &  FV_sec_theta_latitude (1-off_x:row_length+off_x,                &
     &                         1-off_y:rows+off_y)                      &
     &, cos_v_latitude (1-off_x:row_length+off_x,                       &
     &                  1-off_y:n_rows+off_y)

      Real                                                              &
           ! vertical co-ordinate arrays.
     &  eta_theta_levels (0:model_levels)                               &
     &, eta_rho_levels (model_levels)                                   &
     &, r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)

      Real                                                              &
           ! horizontal co-ordinate spacing.
     &  delta_lambda                                                    &
     &, delta_phi

      Real                                                              &
           !VarRes horizontal co-ordinate information
     &  lambda_p(1-halo_i:row_length+halo_i)                            &
     &, phi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)         &
     &, lambda_u(1-halo_i:row_length+halo_i)                            &
     &, phi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)       &
     &, recip_dlambda_u(1-halo_i:row_length+halo_i)                     &
     &, recip_dphi_p(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j)   &
     &, recip_dphi_v(1-halo_i:row_length+halo_i,1-halo_j:n_rows+halo_j) &
     &, wt_lambda_u(1-halo_i:row_length+halo_i)                         &
     &, wt_phi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)
 

      Logical                                                           &
     &  L_do_increment                                                  &
     &, L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup                           &
     &, L_new_tdisc

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real                                                              &
     &  rho (1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
     &       model_levels)                                              &
     &, rho_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &       model_levels)

! Local Variables.

      Integer                                                           &
     &  i, j, k      ! Loop indices

! Local arrays

      Real                                                              &
     &  u1(1-halo_i_wind:row_length+halo_i_wind,                        &
     &     1-halo_j_wind:rows+halo_j_wind, model_levels)                &
     &, v1(1-halo_i_wind:row_length+halo_i_wind,                        &
     &     1-halo_j_wind:n_rows+halo_j_wind, model_levels)              &
     &, w1(1-halo_i_wind:row_length+halo_i_wind,                        &
     &     1-halo_j_wind:rows+halo_j_wind, 0:model_levels)              &
     &, moist     ( 1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,   &
     &          wet_model_levels)                                       &
     &, moist_star( 1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &          wet_model_levels)                                       &
     &, wet_to_dry_n (1-off_x:row_length+off_x, 1-off_y:rows+off_y,     &
     &          model_levels)                                           &
     &, wet_to_dry_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,   &
     &          model_levels)

      Real, Dimension(:,:,:), Allocatable::                             &
     &  moist_np1

#include "parparm.h"
#include "domtyp.h"

       moist = mix_v + mix_cl + mix_cf
       moist_star = mix_v_star + mix_cl_star + mix_cf_star

      If(L_mcr_qcf2)then
        moist      = moist + mix_cf2
        moist_star = moist_star + mix_cf2_star
      endif
      If(L_mcr_qrain)then
        moist      = moist + mix_rain
        moist_star = moist_star + mix_rain_star
      endif
      If(L_mcr_qgraup)then
        moist      = moist + mix_graup
        moist_star = moist_star + mix_graup_star
      endif

      If ( CycleNo > 1 .and. L_new_tdisc ) Then
        Allocate( moist_np1 (1-off_x:row_length+off_x,                  &
     &            1-off_y:rows+off_y, wet_model_levels) )
        moist_np1 = mix_v_np1 + mix_cl_np1 + mix_cf_np1
        If (L_mcr_qcf2)   moist_np1 = moist_np1 + mix_cf2_np1
        If (L_mcr_qrain)  moist_np1 = moist_np1 + mix_rain_np1
        If (L_mcr_qgraup) moist_np1 = moist_np1 + mix_graup_np1
      Else
        Allocate( moist_np1 (1,1,1) )
      End If

! set u1, v1 and w1.

        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
            u1(i,j,k)=u(i,j,k) + alpha_1 * R_u(i,j,k)
            End Do
          End Do
        End Do

        Do k = 1, model_levels
          Do j = 1, n_rows
            Do i = 1, row_length
            v1(i,j,k)=v(i,j,k) + alpha_1 * R_v(i,j,k)
            End Do
          End Do
        End Do

        Do k = 1, model_levels - 1
          Do j = 1, rows
            Do i = 1, row_length
              w1(i,j,k)=w(i,j,k) + alpha_2 * R_w(i,j,k)
            End Do
          End Do
        End Do
        Do k = 0, model_levels, model_levels
          Do j = 1, rows
            Do i =1, row_length
             w1(i,j,k)= 0.
            End Do
          End Do
        End Do

! Do halo swaps for u1 and v1 variables


! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                               &
     &                   u1, row_length, rows,                          &
     &                   model_levels,                                  &
     &                   halo_i_wind, halo_j_wind, fld_type_u, .true.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(u1,row_length,rows,                    &
     &                           model_levels,halo_i_wind,halo_j_wind)


! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                               &
     &                   v1, row_length, n_rows,                        &
     &                   model_levels,                                  &
     &                   halo_i_wind, halo_j_wind, fld_type_v, .true.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(v1,row_length,n_rows,                  &
     &                           model_levels,halo_i_wind,halo_j_wind)

! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                               &
     &                   w1, row_length, rows,                          &
     &                   model_levels+1,                                &
     &                   halo_i_wind, halo_j_wind, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(w1,row_length,rows,model_levels+1,     &
     &                           halo_i_wind,halo_j_wind)

! DEPENDS ON: flux_rho_mix
        Call Flux_Rho_mix(                                              &
     &                 u1, v1, w1, moist, moist_star, moist_np1,        &
     &                 r_theta_levels, r_rho_levels,                    &
     &                 eta_theta_levels, eta_rho_levels,                &
     &                 FV_sec_theta_latitude,                           &
     &                 cos_v_latitude, delta_lambda, delta_phi,         &
     &                 lambda_p, phi_p, lambda_u, phi_v,                &
     &                 wt_lambda_u,  wt_phi_v,                          &
     &                 recip_dlambda_u, recip_dphi_v,                   &
     &                 timestep, NumCycles, CycleNo,                    &
     &                 rows, n_rows, row_length,                        &
     &                 model_levels, wet_model_levels, model_domain,    &
     &                 first_constant_r_rho_level, rims_to_do,          &
     &                 halo_i, halo_j, nproc, gc_proc_row_group,        &
     &                 L_regular, at_extremity, global_row_length,      &
     &                 off_x,off_y,halo_i, halo_j,                      &
     &                 wet_to_dry_n, wet_to_dry_np1, alpha_1,           &
     &                 rho, rho_np1, L_new_tdisc )

      Deallocate ( moist_np1 )

      If (L_do_increment) then
! calculate increment to rho over timestep
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
            inc_rho(i,j,k) = rho(i,j,k) - rho_n(i,j,k)
            End Do
          End Do
        End Do
      endif

      return
      END SUBROUTINE update_rho_mix

#endif
