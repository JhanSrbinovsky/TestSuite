
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine update_fields
      subroutine update_fields(                                         &
     &                         NumCycles,CycleNo,L_new_tdisc,           &
     &                         L_mix_ratio, extrp_weight,               &
     &                         exner_rho_levels,                        &
     &                         exner_theta_levels,                      &
     &                         u, u_np1, v, v_np1, w, w_np1,            &
     &                         u_adv, v_adv, w_adv,                     &
     &                         theta, q, qcl, qcf,                      &
     &                         qcf2, qrain, qgraup,                     &
     &                         cf, cfl, cff,                            &
     &                         exner_prime, R_u, R_v, R_w,              &
     &                         theta_star, theta_np1, dtheta_dr_term,   &
     &                         q_star, q_np1, qcl_star, qcl_np1,        &
     &                         qcf_star, qcf_np1, qcf2_star, qcf2_np1,  &
     &                         qrain_star, qrain_np1,                   &
     &                         qgraup_star, qgraup_np1,                 &
     &                         cf_star, cfl_star, cff_star,             &
     &                         row_length, rows, n_rows, model_levels,  &
     &                         offx, offy, halo_i, halo_j,              &
     &                         sin_theta_longitude, cos_theta_longitude,&
     &                         mag_vector_np, dir_vector_np,            &
     &                         mag_vector_sp, dir_vector_sp,            &
     &                         global_row_length, gc_proc_row_group,    &
     &                         at_extremity, datastart, model_domain,   &
     &                         alpha_2 , timestep, wet_levels,          &
     &                         i_start, i_end, j_start, j_end,j0,j1,    &
     &                         delta_lambda, lambda_p, lambda_u,        &
     &                         inc_t, inc_u, inc_v, inc_w,              &
     &                         L_do_t, L_do_inc_vels, L_pc2,            &
     &                         L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,   &
     &                         L_regular )

! Purpose:
!          update fields in a subroutine to tidy ATM_STEP
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
!  5.2    15/10/00  This deck introduced                 Andy Malcolm
!  5.4    22/07/02  Update cloud fractions for PC2       D. Wilson
!  5.5    03/02/03   Update extra microphysics prognostics. R.M.Forbes
!  6.2    15/03/06  Enable iteration of phys-dyn and use of different
!    interpolation/extrapolation formulae.

!  6.1    15/09/04  Change two array names               A.Malcolm
!  6.2    04/10/05  Changes for cycling semi-Lagrangian scheme.
!                                                        M. Diamantakis
!  6.2    25/12/05  Variable resolution changes          Yongming Tang
!  6.4    11/12/06  1-point halo for q vars in iterative SL scheme 
!                                                        M. Diamantakis
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.
      INTEGER                                                           &
     &       NumCycles,                                                 &
     &       CycleNo,                                                   &
     &       ROW_LENGTH,                                                &
                                  ! IN: No of points per local row
     &       ROWS,                                                      &
                                  ! IN: No of local (theta) rows
     &       n_ROWS,                                                    &
                                  ! IN: No of local (v) rows
     &       MODEL_LEVELS,                                              &
                                  ! IN: No of model levels
     &       WET_LEVELS,                                                &
                                  ! IN: No of moist-levels
     &       Offx,                                                      &
                                  ! standard halo size in East-West
     &       Offy,                                                      &
                                  ! standard halo size in North-South
     &       halo_i,                                                    &
                                  ! extended halo size in East-West
     &       halo_j               ! extended halo size in North-South

      Real                                                              &
           !VarRes horizontal co-ordinate information
     &  lambda_p(1-halo_i:row_length+halo_i)                            &
     &, lambda_u(1-halo_i:row_length+halo_i)

      Real, intent(inout) ::                                            &
     &  u(1-offx:row_length+offx, 1-offy:rows+offy,                     &
     &      model_levels)                                               &
     &, u_np1(1-offx:row_length+offx, 1-offy:rows+offy,                 &
     &      model_levels)                                               &
     &, v(1-offx:row_length+offx, 1-offy:n_rows+offy,                   &
     &      model_levels)                                               &
     &, v_np1(1-offx:row_length+offx, 1-offy:n_rows+offy,               &
     &      model_levels)                                               &
     &, w(1-offx:row_length+offx, 1-offy:rows+offy,                     &
     &      0:model_levels)                                             &
     &, w_np1(1-offx:row_length+offx, 1-offy:rows+offy,                 &
     &      0:model_levels)                                             &
     &, u_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &          model_levels)                                           &
     &, v_adv(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,       &
     &          model_levels)                                           &
     &, w_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &          0:model_levels)                                         &
     &, theta(1-offx:row_length+offx, 1-offy:rows+offy,                 &
     &          model_levels)                                           &
     &, theta_np1(1-offx:row_length+offx,                               &
     &               1-offy:rows+offy, model_levels)                    &
     &, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &      wet_levels)                                                 &
     &, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_levels)                                               &
     &, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_levels)                                               &
     &, qcf2(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &        wet_levels)                                               &
     &, qrain(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &        wet_levels)                                               &
     &, qgraup(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &        wet_levels)                                               &
     &, cf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
     &      wet_levels)                                                 &
     &, cfl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_levels)                                               &
     &, cff(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_levels)                                               &
     &, exner_rho_levels(1-offx:row_length+offx,                        &
     &                   1-offy:rows+offy, model_levels+1)              &
     &, exner_theta_levels(1-offx:row_length+offx,                      &
     &                     1-offy:rows+offy, model_levels)

      Real, intent(in) ::                                               &
     &  theta_star(1-offx:row_length+offx,                              &
     &               1-offy:rows+offy, model_levels)                    &
     &, q_star(1-offx:row_length+offx,                                  &
     &           1-offy:rows+offy, wet_levels)                          &
     &, qcl_star(1-offx:row_length+offx,                                &
     &             1-offy:rows+offy, wet_levels)                        &
     &, qcf_star(1-offx:row_length+offx,                                &
     &             1-offy:rows+offy, wet_levels)                        &
     &, qcf2_star(1-offx:row_length+offx,                               &
     &             1-offy:rows+offy, wet_levels)                        &
     &, qrain_star(1-offx:row_length+offx,                              &
     &             1-offy:rows+offy, wet_levels)                        &
     &, qgraup_star(1-offx:row_length+offx,                             &
     &             1-offy:rows+offy, wet_levels)                        &
     &, cf_star(1-offx:row_length+offx,                                 &
     &           1-offy:rows+offy, wet_levels)                          &
     &, cfl_star(1-offx:row_length+offx,                                &
     &             1-offy:rows+offy, wet_levels)                        &
     &, cff_star(1-offx:row_length+offx,                                &
     &             1-offy:rows+offy, wet_levels)

      Real, intent(in) ::                                               &
     &  cos_theta_longitude (row_length, rows)                          &
     &, sin_theta_longitude (row_length, rows)

      Real, intent(in) ::                                               &
     &  R_u(1-offx:row_length+offx, 1-offy:rows+offy,                   &
     &        model_levels)                                             &
     &, R_v(1-offx:row_length+offx, 1-offy:n_rows+offy,                 &
     &        model_levels)                                             &
     &, R_w(row_length, rows, model_levels)

      Logical                                                           &
     &  L_do_t                                                          &
     &, L_do_inc_vels                                                   &
     &, L_regular                                                       &
     &, L_mcr_qcf2                                                      &
                      ! T => Use second ice prognostic
     &, L_mcr_qrain                                                     &
                       ! T => Use rain prognostic
     &, L_mcr_qgraup                                                    &
                       ! T => Use graupel prognostic
     &, L_pc2                                                           &
               ! Use PC2 cloud scheme
     &, L_new_tdisc                                                     &
     &, L_mix_ratio


      Real, intent(Out) ::                                              &
     &  inc_t(1-offx:row_length+offx,                                   &
     &        1-offy:rows+offy, model_levels)                           &
     &, inc_u(1-offx:row_length+offx, 1-offy:rows+offy,                 &
     &        model_levels)                                             &
     &, inc_v(1-offx:row_length+offx, 1-offy:n_rows+offy,               &
     &        model_levels)                                             &
     &, inc_w(row_length, rows, model_levels)                           &
     &, q_np1(1-offx:row_length+offx,                                   &
     &           1-offy:rows+offy, wet_levels)                          &
     &, qcl_np1(1-offx:row_length+offx,                                 &
     &           1-offy:rows+offy, wet_levels)                          &
     &, qcf_np1(1-offx:row_length+offx,                                 &
     &           1-offy:rows+offy, wet_levels)                          &
     &, qcf2_np1 (1-offx:row_length+offx, 1-offy:rows+offy,             &
     &          wet_levels)                                             &
     &, qrain_np1 (1-offx:row_length+offx, 1-offy:rows+offy,            &
     &          wet_levels)                                             &
     &, qgraup_np1 (1-offx:row_length+offx, 1-offy:rows+offy,           &
     &          wet_levels)

      Real, intent(in) ::                                               &
     &  exner_prime (1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)                                          &
     &, dtheta_dr_term (row_length, rows, model_levels)

! time interpolation interpolation weight for first SL iteration.
      Real extrp_weight

      Real                                                              &
     &  mag_vector_np (model_levels)                                    &
     &, dir_vector_np (model_levels)                                    &
     &, mag_vector_sp (model_levels)                                    &
     &, dir_vector_sp (model_levels)                                    &
     &, delta_lambda                                                    &
     &, alpha_2                                                         &
     &, timestep

        Integer                                                         &
     &  i_start, i_end                                                  &
     &, j_start, j_end  ,j0,j1                                          &
     &, model_domain , datastart(3)                                     &
     &, global_row_length                                               &
     &, gc_proc_row_group

      Integer                                                           &
     &  gi

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

!    local variables
      Integer                                                           &
     & i,j,k

!========================== COMDECK PARPARM ====================
!   Description:
!
!   This COMDECK contains PARAMETERs for the mpp-UM
!
!   Two sets of parameters are set up -
!     i)  for the mpp-UM itself.
!     ii) for the interface to the Message Passing Software.
!
      !=================================================================
      ! Parameters needed for the mpp-UM
      !=================================================================
      ! maximum number of spatial dimensions
      INTEGER,PARAMETER:: Ndim_max = 3 ! 3d data

      ! number of different halo types
      INTEGER,PARAMETER:: NHalo_max = 3 ! for N.D. atmos. model

      INTEGER,PARAMETER:: halo_type_single   = 1
      INTEGER,PARAMETER:: halo_type_extended = 2
      INTEGER,PARAMETER:: halo_type_no_halo  = 3

! FLDTYPE definitions for the different field types recognised on the
! decomposition
      INTEGER,PARAMETER:: Nfld_max=7 ! maximum number of field types
      INTEGER,PARAMETER:: fld_type_p=1       ! grid on P points
      INTEGER,PARAMETER:: fld_type_u=2       ! grid on U points
      INTEGER,PARAMETER:: fld_type_v=3       ! grid on V points
      INTEGER,PARAMETER:: fld_type_comp_wave  = 4
                              ! Compressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_full_wave  = 5
                              ! Uncompressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_rim_wave   = 6
                              ! Boundary data for WAM Wave Field
      INTEGER,PARAMETER:: fld_type_r=7       ! grid on river points
      INTEGER,PARAMETER:: fld_type_unknown=-1! non-standard grid
! FLDTYPE end

      ! Used in addressing to indicate if calculation is for a local or
      ! global (ie. disk dump) size

      INTEGER,PARAMETER:: local_data=1
      INTEGER,PARAMETER:: global_dump_data=2

      ! maximum permitted size of a halo
      INTEGER,PARAMETER:: Max_Halo_Size=10

      !=================================================================
      ! Parameters needed for the Message Passing Software
      !=================================================================
      INTEGER,PARAMETER:: Maxproc = 512 ! Max number of processors

      ! Processor addresses in the neighbour array
      INTEGER,PARAMETER:: PNorth   = 1
      INTEGER,PARAMETER:: PEast    = 2
      INTEGER,PARAMETER:: PSouth   = 3
      INTEGER,PARAMETER:: PWest    = 4

      ! Value in neighbour array if the domain has  no neighbour in this
      ! direction. Otherwise the value will be the tid of the neighbor
      INTEGER,PARAMETER:: NoDomain = -1

      INTEGER,PARAMETER:: BC_STATIC   = 1 ! Static boundary conditions
      INTEGER,PARAMETER:: BC_CYCLIC   = 2 ! Cyclic boundary conditions
      INTEGER,PARAMETER:: BC_OVERPOLE = 3 ! Transfer over pole
! PARPARM end
! DOMTYP contains different model domain types
!
! Author : P.Burton
! History:
! Version  Date      Comment.
! 5.0      15/04/99  New comdeck
! 5.2      15/11/00  add bi_cyclic_lam domain   A. Malcolm

      INTEGER,PARAMETER:: mt_global        = 1
      INTEGER,PARAMETER:: mt_lam           = 2
      INTEGER,PARAMETER:: mt_cyclic_lam    = 3
      INTEGER,PARAMETER:: mt_bi_cyclic_lam = 4
      INTEGER,PARAMETER:: mt_single_column = 5
! DOMTYP end

! Check if we are in the final SL advection cycle.
! If yes, then perform a full update_fields.

      If ( CycleNo == NumCycles ) Then

        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              exner_rho_levels(i,j,k) =                                 &
     &        exner_rho_levels(i,j,k) + exner_prime(i,j,k)
            End Do
          End Do
        End Do

!TD  Upper boundary condition d(exner_primed)/dr at top theta level
!TD  is assumed in Helmholtz equation. Exner_rho_primed(model_levels)
!TD  also needs to be added to Exner_rho at level above to satisfy ubc
        k =  model_levels
        Do j = 1, rows
          Do i = 1, row_length
            exner_rho_levels(i,j,k+1) =                                 &
     &      exner_rho_levels(i,j,k+1) + exner_prime(i,j,k)
          End Do
        End Do

        if( L_do_inc_vels ) then
! store EoT increments
        inc_u=u   ! store u at time level n (because of polar_vector_win
        inc_v=R_v   ! store v increment
        inc_w=R_w   ! store w increment
        endif

! Use wind field increments as produced by Helmholtz equation

! update winds first, then extrapolate them at tn+dt/2
! and store in _adv variables to use at next timestep.
        Do k = 1, model_levels
          Do j = j0, j1
            Do i = 1, row_length
              u_adv(i,j,k) = u(i,j,k) + extrp_weight*R_u(i,j,k)
              u(i,j,k)     = u(i,j,k) + R_u(i,j,k)
            End Do
          End Do
          Do j = 1, n_rows
            Do i = 1, row_length
              v_adv(i,j,k) = v(i,j,k) + extrp_weight*R_v(i,j,k)
              v(i,j,k)     = v(i,j,k) + R_v(i,j,k)
            End Do
          End Do
        End Do

        Do k = 1, model_levels - 1
          Do j = 1, rows
            Do i = 1, row_length
              w_adv(i,j,k) = w(i,j,k) + extrp_weight*R_w(i,j,k)
              w(i,j,k)     = w(i,j,k) + R_w(i,j,k)
            End Do
          End Do
        End Do

        If (model_domain == mt_global ) Then
! calculate u at the poles.

! DEPENDS ON: polar_vector_wind_n
          Call Polar_vector_wind_n(                                     &
     &                       v,                                         &
     &                       sin_theta_longitude,                       &
     &                       cos_theta_longitude, row_length,           &
     &                       n_rows, model_levels, mag_vector_np,       &
     &                       dir_vector_np, mag_vector_sp,              &
     &                       dir_vector_sp,                             &
     &                       offx, offy, global_row_length,             &
     &                       gc_proc_row_group, at_extremity )

          If ( at_extremity(PSouth) .and. L_regular ) Then
            j=1
            Do k = 1,model_levels
              Do i = 1, row_length
                gi = datastart(1) + i - 1
                u(i,j,k) = - mag_vector_sp(k) * sin ( (gi-.5)*          &
     &                       delta_lambda - dir_vector_sp(k))
              End Do
            End Do
          ElseIf (at_extremity(PSouth) ) Then
            j=1
            Do k = 1,model_levels
              Do i = 1, row_length
                u(i,j,k) = - mag_vector_sp(k) *                         &
     &                       sin ( lambda_u(i)- dir_vector_sp(k) )
              End Do
            End Do
          EndIf !  at_extremity(PSouth) .and. L_regular

          If ( at_extremity(PNorth) .and. L_regular ) Then
            j=rows
            Do k = 1,model_levels
              Do i = 1, row_length
                gi = datastart(1) + i - 1
                u(i,j,k) = mag_vector_np(k) *                           &
     &                     sin((gi-.5)*delta_lambda-dir_vector_np(k))
              End Do
            End Do
          ElseIf (at_extremity(PNorth) ) Then
            j=rows
            Do k = 1,model_levels
              Do i = 1, row_length
                u(i,j,k) = mag_vector_np(k) *                           &
     &                     sin ( lambda_u(i)- dir_vector_np(k) )
              End Do
            End Do
          End If  !  at_extremity(PNorth) .and. L_regular

        End If ! model_domain == mt_global

        if( L_do_inc_vels ) then
! calculate increment to u field
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              inc_u(i,j,k) = u(i,j,k) - inc_u(i,j,k)
            End Do
          End Do
        End Do
        endif

        If(L_do_t) then
! store T at time level n
        Do k = 1, model_levels
          Do j = j_start, j_end
            Do i = i_start, i_end
              inc_t(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
            End Do
          End Do

        End Do
        endif

! set theta at time level n+1
! do not update boundaries in LAM

        Do k = 1, model_levels - 1
          Do j = j_start, j_end
            Do i = i_start, i_end
              theta(i,j,k) = theta_star(i,j,k) - alpha_2                &
     &                     * timestep * R_w(i,j,k)                      &
     &                     * dtheta_dr_term(i,j,k)
            End Do
          End Do
        End Do

        k = model_levels
        Do j = j_start, j_end
          Do i = i_start, i_end
            theta(i,j,k) = theta_star(i,j,k)
          End Do
        End Do

        Do k = 1, wet_levels
          Do j = j_start, j_end
            Do i = i_start, i_end
              q(i,j,k)   = q_star(i,j,k)
              qcl(i,j,k) = qcl_star(i,j,k)
              qcf(i,j,k) = qcf_star(i,j,k)
              If (L_pc2) then
                cf(i,j,k)  = cf_star(i,j,k)
                cfl(i,j,k) = cfl_star(i,j,k)
                cff(i,j,k) = cff_star(i,j,k)
              End If  ! L_pc2
            End Do
          End Do
        End Do

        If (L_mcr_qcf2)                                                 &
     &      qcf2(i_start:i_end, j_start:j_end, :) =                     &
     &      qcf2_star(i_start:i_end, j_start:j_end, :)

        If (L_mcr_qrain)                                                &
     &      qrain(i_start:i_end, j_start:j_end, :) =                    &
     &      qrain_star(i_start:i_end, j_start:j_end, :)

        If (L_mcr_qgraup)                                               &
     &      qgraup(i_start:i_end, j_start:j_end, :) =                   &
     &      qgraup_star(i_start:i_end, j_start:j_end, :)

      Else

        If ( L_new_tdisc ) Then

! Update u_adv, v_adv, w_adv. These are winds interpolated at tn+dt/2.
! These values will be used to find the departure point at next cycle.
! Also compute u, v, w estimate at tn+1 (at current cycle).
          Do k = 1, model_levels
            Do j = j0, j1
              Do i = 1, row_length
                u_adv(i,j,k) = u(i,j,k) + 0.5*R_u(i,j,k)
                u_np1(i,j,k) = u(i,j,k) + R_u(i,j,k)
              End Do
            End Do
            Do j = 1, n_rows
              Do i = 1, row_length
                v_adv(i,j,k) = v(i,j,k) + 0.5*R_v(i,j,k)
                v_np1(i,j,k) = v(i,j,k) + R_v(i,j,k)
              End Do
            End Do
          End Do

          Do k = 1, model_levels - 1
            Do j = 1, rows
              Do i = 1, row_length
                w_adv(i,j,k) = w(i,j,k) + 0.5*R_w(i,j,k)
                w_np1(i,j,k) = w(i,j,k) + R_w(i,j,k)
              End Do
            End Do
          End Do

          If (model_domain == mt_global ) Then
! calculate u at the poles.

! DEPENDS ON: polar_vector_wind_n
            Call Polar_vector_wind_n(                                   &
     &                       v_np1,                                     &
     &                       sin_theta_longitude,                       &
     &                       cos_theta_longitude, row_length,           &
     &                       n_rows, model_levels, mag_vector_np,       &
     &                       dir_vector_np, mag_vector_sp,              &
     &                       dir_vector_sp,                             &
     &                       offx, offy, global_row_length,             &
     &                       gc_proc_row_group, at_extremity )

            If (at_extremity(PSouth) ) Then
              j=1
              Do k = 1,model_levels
                Do i = 1, row_length
                  gi = datastart(1) + i - 1
                  u_np1(i,j,k) = - mag_vector_sp(k) * sin ( (gi-.5)*    &
     &                             delta_lambda - dir_vector_sp(k))
                End Do
              End Do
            End If

            If (at_extremity(PNorth) ) Then
              j=rows
              Do k = 1,model_levels
                Do i = 1, row_length
                  gi = datastart(1) + i - 1
                  u_np1(i,j,k) = mag_vector_np(k) *                     &
     &                     sin((gi-.5)*delta_lambda-dir_vector_np(k))
                End Do
              End Do
            End If

          End If ! model_domain == mt_global

! set theta_np1 (theta at time level n+1 at end of iteration)
! always update boundaries (treated as theta_star, no-blending)
          Do k = 1, model_levels - 1
            Do j = 1, rows
              Do i = 1, row_length
                theta_np1(i,j,k) = theta_star(i,j,k) - alpha_2          &
     &                     * timestep*R_w(i,j,k)*dtheta_dr_term(i,j,k)
              End Do
            End Do
          End Do

          k = model_levels
          Do j = 1, rows
            Do i = 1, row_length
              theta_np1(i,j,k) = theta_star(i,j,k)
            End Do
          End Do

! set q_np1 etc (q at time level n+1 at end of iteration)
! always update boundaries (treated as q_star, no-blending)
          If ( .NOT. L_mix_ratio ) Then

            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  q_np1(i,j,k) = q_star(i,j,k)
                  qcl_np1(i,j,k) = qcl_star(i,j,k)
                  qcf_np1(i,j,k) = qcf_star(i,j,k)
                End Do
              End Do
            End Do

            If (L_mcr_qcf2)                                             &
     &        qcf2_np1(1:row_length,1:rows,:) =                         &
     &              qcf2_star(1:row_length,1:rows,:)

            If (L_mcr_qrain)                                            &
     &        qrain_np1(1:row_length,1:rows,:) =                        &
     &              qrain_star(1:row_length,1:rows,:)

            If (L_mcr_qgraup)                                           &
     &        qgraup_np1(1:row_length,1:rows,:) =                       &
     &              qgraup_star(1:row_length,1:rows,:)

          End If

        Else

! Update u_adv, v_adv, w_adv. These are winds interpolated at tn+dt/2.
! These values will be used to find the departure point at next cycle.
          Do k = 1, model_levels
            Do j = j0, j1
              Do i = 1, row_length
                u_adv(i,j,k) = u(i,j,k) + 0.5*R_u(i,j,k)
              End Do
            End Do
            Do j = 1, n_rows
              Do i = 1, row_length
                v_adv(i,j,k) = v(i,j,k) + 0.5*R_v(i,j,k)
              End Do
            End Do
          End Do

          Do k = 1, model_levels - 1
            Do j = 1, rows
              Do i = 1, row_length
                w_adv(i,j,k) = w(i,j,k) + 0.5*R_w(i,j,k)
              End Do
            End Do
          End Do

        End If ! L_new_tdisc

      End If ! CycleNo == NumCycles

! The following code needs to be executed regardless
! cycle no.

      If (model_domain == mt_global ) Then

! calculate vector wind at poles for u_adv.
! DEPENDS ON: polar_vector_wind_n
        Call Polar_vector_wind_n(                                       &
     &                       v_adv,                                     &
     &                       sin_theta_longitude,                       &
     &                       cos_theta_longitude, row_length,           &
     &                       n_rows, model_levels, mag_vector_np,       &
     &                       dir_vector_np, mag_vector_sp,              &
     &                       dir_vector_sp,                             &
     &                       halo_i, halo_j, global_row_length,         &
     &                       gc_proc_row_group, at_extremity)

        If ( at_extremity(PSouth) .and. L_regular ) Then
          j=1
          Do k = 1,model_levels
            Do i = 1, row_length
              gi = datastart(1) + i - 1
              u_adv(i,j,k) = - mag_vector_sp(k) *                       &
     &                sin ( (gi-.5)*delta_lambda - dir_vector_sp(k) )
            End Do
          End Do
        ElseIf ( at_extremity(PSouth) ) Then
          j=1
          Do k = 1,model_levels
            Do i = 1, row_length
              u_adv(i,j,k) = - mag_vector_sp(k) *                       &
     &                     sin ( lambda_u(i)- dir_vector_sp(k) )
            End Do
          End Do
        End If  !  at_extremity(PSouth) .and. L_regular

        If ( at_extremity(PNorth) .and. L_regular ) Then
          j=rows
          Do k = 1,model_levels
            Do i = 1, row_length
              gi = datastart(1) + i - 1
              u_adv(i,j,k) = mag_vector_np(k) *                         &
     &                sin ( (gi-.5)*delta_lambda - dir_vector_np(k) )
            End Do
          End Do
        ElseIf (at_extremity(PNorth) ) Then
          j=rows
          Do k = 1,model_levels
            Do i = 1, row_length
              u_adv(i,j,k) = mag_vector_np(k) *                         &
     &                     sin ( lambda_u(i)- dir_vector_np(k) )
            End Do
          End Do
        End If  !  at_extremity(PNorth).and. L_regular

      End If


      return
      END SUBROUTINE update_fields

