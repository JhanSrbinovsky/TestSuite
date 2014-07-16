
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Calc_PMSL

      Subroutine Calc_PMSL(                                             &
     &                   theta, exner_theta_levels, p,                  &
     &                   r_theta_levels, r_rho_levels, eta_theta_levels,&
     &                   g, R, Lapse_Rate, earth_radius,                &
     &                   row_length, rows, model_levels,                &
     &                   boundary_layer_levels,                         &
     &                   off_x, off_y, halo_i, halo_j,                  &
     &                   pmsl,p_star,                                   &
!
!-- New inputs to be added to call from diagnostics_end--------
!
     &                       me, n_proc, n_procx, n_procy,              &
     &                       g_datastart,                               &
     &                       neighbour, at_extremity,                   &
     &                       all_proc_group,model_domain,               &
     &                       delta_lambda, delta_phi,                   &
     &                       Cp, npmsl_height,                          &
     &                       sec_theta_latitude,                        &
     &                       global_row_length, global_rows)
!
!-------------------------------------------------------------
!

! Purpose:
!          Calculates Pressure at mean sea level.
!
! Method:
!          Modified version of that described in equations 3.9 and 3.11
!          in Unified Model Documentation
!          Paper Number S1. A. Dickinson
!
! Original Programmer: Mark H. Mawson
! Current code owner: Andrew J. Malcolm
!
! History:
!  Version     Date    Comment
! ----     -------     -------
! 6/6/01   5.3         Added portable alternative to T3E code
!                      Removed hardwired PI     P.Burton
!    5.3     07/11/01  Add Clive Wilson change to PMSL calc   A.Malcolm
!    5.4     01/05/02  Impovements for LAM calculations. P.Selwood
!                      (J. Bornemann's initial mod)
!    6.0     12/06/03  fix bug found in non-t3e code          A.Malcolm
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. i.e.: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                        ! number of points on a row
     &, rows                                                            &
                        ! number of rows of data
     &, model_levels                                                    &
                        ! number of levels of data
     &, boundary_layer_levels                                           &
                              ! number of boundary layer levels
     &, off_x                                                           &
     &, off_y                                                           &
     &, halo_i                                                          &
     &, halo_j

      Integer                                                           &
     &  me                                                              &
                            !IN. Processor number
     &, n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, g_datastart (3,0:n_proc-1)                                      &
     &, neighbour(4)                                                    &
                             ! Array with the Ids of the four neighbours
                             ! in the horizontal plane
     &, all_proc_group                                                  &
                       ! Group id for all processors
     &, model_domain                                                    &
                         ! indicator as to model type, ie global, lam
     &, global_row_length                                               &
                            !IN. NUMBER OF points on a global row
     &, global_rows         !IN. NUMBER OF global rows

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

      Real                                                              &
     & Cp

      Real                                                              &
     &  R                                                               &
     &, g                                                               &
     &, Lapse_Rate                                                      &
     &, earth_radius                                                    &
     &, delta_lambda                                                    &
     &, delta_phi                                                       &
     &, npmsl_height ! Orographic height above which relaxation occurs

      Real                                                              &
     &  p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &          model_levels)                                           &
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y, model_levels)            &
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j,0:model_levels)             &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j, model_levels)                &
     &, sec_theta_latitude (1-off_x:row_length+off_x,                   &
     &                      1-off_y:rows+off_y)                         &
     &, eta_theta_levels(0:model_levels) ! vertical grid for theta vars


! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
     &  pmsl (row_length, rows)                                         &
     &, p_star (row_length, rows)

! Local variables

      Integer                                                           &
     & i, j                                                             &
     &, k, LevelUpper

      Real                                                              &
     &  T_ref_level_1                                                   &
     &, T_at_mean_sea_level (row_length, rows)                          &
     &, power                                                           &
     &, phi_star (row_length, rows)                                     &
     &, cos_p_latitude (row_length, rows)                               &
     &, height_domain

      Logical l_pmsl_new

      REAL, PARAMETER :: UpperHeight = 2000.0
      ! Height of upper level above ground

! ----------------------------------------------------------------------
! Section 1. Calculate reference temperature at lowest pressure level
!            from value at level just above top of boundary layer and
!            assuming a fixed lapse rate with height.
!            Uses equation 3.9 from UMDP S1.
!            Calculate PMSL using equation 3.11 from UMDP S1 but with
!            level 1 values replacing * ones. the expression on the top
!            of the quotient in equation 3.11 is in fact simply the
!            temperature at mean sea level.
! ----------------------------------------------------------------------

! find upper level (replacing level just above top of boundary layer!)
! height_domain is the top model height and is the same everywhere
      height_domain = r_theta_levels(1,1,model_levels) - Earth_radius

      DO k = 1, model_levels
        LevelUpper = k
        IF ( height_domain * eta_theta_levels(k) > UpperHeight) EXIT
      END DO
      l_pmsl_new=.true.

      power = g / (R * Lapse_Rate)


      Do j = 1, rows
        Do i = 1, row_length

          T_ref_level_1 = theta(i,j,LevelUpper) *                       &
     &                 exner_theta_levels(i,j,LevelUpper)               &
     &                  + Lapse_Rate *                                  &
     &                  (r_theta_levels(i,j,LevelUpper) -               &
     &                   r_rho_levels(i,j,1) )

          T_at_mean_sea_level(i,j) = T_ref_level_1                      &
     &                          + Lapse_Rate *                          &
     &                           (r_rho_levels(i,j,1) -                 &
     &                            earth_radius )

          pmsl(i,j) = p(i,j,1) *                                        &
     &               (T_at_mean_sea_level(i,j) / T_ref_level_1 )        &
     &               ** power

        End Do
      End Do
!
!   New method for pmsl from Unified Model
!
      if (l_pmsl_new) then

      Do j = 1, rows
        Do i = 1, row_length
          phi_star(i,j)=g * (r_theta_levels(i,j,0)-earth_radius)
          cos_p_latitude(i,j)= 1./sec_theta_latitude(i,j)

        End Do
      End Do

! DEPENDS ON: calc_npmsl
      call  Calc_NPMSL(pmsl,p_star,                                     &
     &                 phi_star,theta,T_at_mean_sea_level,              &
     &                 cos_p_latitude,delta_lambda,delta_phi,           &
     &                 row_length,rows,                                 &
     &                 global_row_length,global_rows,                   &
     &                 g, R, Lapse_Rate, earth_radius,                  &
     &                 me, n_proc, n_procx, n_procy,                    &
     &                 off_x, off_y,                                    &
     &                 g_datastart(1,me),                               &
     &                 neighbour, at_extremity,                         &
     &                 all_proc_group,model_domain,                     &
     &                 Cp, npmsl_height)

      endif

      return
      END SUBROUTINE Calc_PMSL


