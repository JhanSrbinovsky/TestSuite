
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine filttarq
      subroutine filttarq(                                              &
     &                    field, w_adv,                                 &
     &                    row_length, rows, n_rows,                     &
     &                    levels, model_levels,                         &
     &                    first_constant_r_rho_level,                   &
     &                    r_theta_levels, r_rho_levels,                 &
     &                    sec_theta_latitude, cos_v_latitude,           &
     &                    off_x, off_y, halo_i, halo_j,                 &
     &                    halo_i_in, halo_j_in,                         &
     &                    at_extremity, proc_row_group,                 &
     &                    delta_lambda, delta_phi,                      &
     &                    model_domain, global_row_length,              &
     &                    w_limit, horizontal_level,                    &
     &                    factor, test_level, active_levels,            &
     &                    start_level, end_level, L_scale_w,            &
     &                    L_diag_w, w_local_mask)

! Purpose:
!          Applies conservative horizontal diffusion to a field
!          subject to w > w_limit, and has a steep slope test
!
! Method:
!          Is described in ;
!
!
! Original Programmer: Terry Davies
! Current code owner: Andrew J. Malcolm
!
! History:
! Version   Date      Comment
! ----     ------     -------
! 6.2   25/12/05  Based on tardiff_q_wss but called from ni_filter
!                                                       Terry Davies
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Logical, Intent(In) ::                                            &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_diag_w                                                        &
                         ! diagnostic control
     &, L_scale_w        ! scale diffusion coefficient by w_max/w_limit

      Integer, Intent(In) ::                                            &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, n_rows                                                          &
                         ! number of v rows.
     &, levels                                                          &
                         ! number of levels in field
     &, model_levels                                                    &
                         ! number of model levels
     &, first_constant_r_rho_level                                      &
     &, active_levels                                                   &
                         ! number of levels to apply diffusion
     &, off_x                                                           &
                         ! Size of small halo in i
     &, off_y                                                           &
                         ! Size of small halo in j.
     &, halo_i                                                          &
                         ! Size of halo in i direction.
     &, halo_j                                                          &
                         ! Size of halo in j direction.
     &, halo_i_in                                                       &
                         ! Halo size of field in i direction
     &, halo_j_in                                                       &
                         ! Halo size of field in j direction
     &, global_row_length                                               &
     &, proc_row_group

      Integer, Intent(In) ::                                            &
     &  model_domain                                                    &
                         ! holds integer code for model domain
     &, test_level                                                      &
     &, start_level                                                     &
     &, end_level                                                       &
     &, horizontal_level   !  upper test level for steep slopes

      Real, Intent(In) ::                                               &
     &  w_limit                                                         &
                 ! Vertical velocity test value
     &, factor   ! effective diffusion coefficient

      Real, Intent(In) ::                                               &
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)

      Real, Intent(In) ::                                               &
     &  sec_theta_latitude (1-off_x:row_length+off_x,                   &
     &                      1-off_y:rows+off_y)                         &
     &, cos_v_latitude (1-off_x:row_length+off_x,                       &
     &                  1-off_y:n_rows+off_y)

      Real, Intent(In) ::                                               &
     &  w_adv (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j         &
     &,        0:model_levels)

      Real, Intent(In) ::                                               &
     &  delta_lambda                                                    &
     &, delta_phi

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real, Intent(InOut) ::                                            &
     &  field (1-halo_i_in:row_length+halo_i_in,                        &
     &         1-halo_j_in:rows+halo_j_in, levels)                      &
     &, w_local_mask( row_length, rows)

! Local Variables.

      Integer                                                           &
     &  i, j, k, ka                                                     &
                                  ! Loop indices
     &, j_start, j_stop                                                 &
                              ! Loop indices
     &, info                                                            &
     &, level_flat

      Real                                                              &
     &  recip_delta_phi                                                 &
     &, recip_delta_phi2                                                &
     &, recip_delta_lam2                                                &
     &, pole_term

! Local arrays

      Real                                                              &
     &  lambda_term(row_length, rows, active_levels )                   &
     &, phi_term(row_length, rows, active_levels )                      &
     &, delta_z(1-off_x:row_length+off_x,                               &
     &          1-off_y:rows+off_y, active_levels )                     &
     &, temp(1-off_x:row_length+off_x, 1-off_y:rows+off_y)              &
     &, l_s_poles(row_length, active_levels )                           &
     &, l_n_poles(row_length, active_levels )                           &
     &, sum_pole(active_levels )                                        &
     &, diffc_i(1-off_x:row_length+off_x, 1-off_y:rows+off_y)           &
     &, diffc_j(1-off_x:row_length+off_x, 1-off_y:rows+off_y)           &
     &, mask_i(1-off_x:row_length, rows, active_levels)                 &
     &, mask_j(row_length, 1-off_y:rows, active_levels)                 &
     &, w_max(1-off_x:row_length+off_x, 1-off_y:rows+off_y)

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

! No External Routines:

      recip_delta_phi = 1.0/ delta_phi
      recip_delta_phi2 =  recip_delta_phi * recip_delta_phi
      recip_delta_lam2 = 1. / (delta_lambda * delta_lambda)

      j_start = 1
      j_stop = rows
      If (model_domain == mt_Global ) Then
        pole_term = 8. * recip_delta_phi / global_row_length
        If (at_extremity(PSouth)) j_start = 2
        If (at_extremity(PNorth)) j_stop = rows - 1
      End If

! initialise w_local_mask to zero

      If (L_diag_w) Then
        Do j = 1,rows
          Do i = 1, row_length
            w_local_mask(i,j) = 0.0
          End Do
        End Do
      EndIf !  L_diag_w

! ----------------------------------------------------------------------
! Section 1.   Calculate delta_z and diffusion coefficients
! ----------------------------------------------------------------------
! Originally calculate D(r)/D(eta) about q levels BUT
! instead calculate delta_z about q levels
! We can omit the d(eta) part since this is constant and can be
! brought out through the derivatives and cancelled
      Do k = 1, active_levels
        ka = k + start_level - 1
        If(k < model_levels) then
          Do j = 1-off_y, rows+off_y
            Do i = 1-off_x, row_length+off_x
               delta_z(i,j,k) = (r_rho_levels(i,j,ka+1) -               &
     &                           r_rho_levels(i,j,ka) )
            End Do
          End Do
        Else ! k = model_levels
! can use any constant value for delta_z since it will cancel
          Do j = 1-off_y, rows+off_y
            Do i = 1-off_x, row_length+off_x
              delta_z(i,j,k) = 1.0
            End Do
          End Do
        EndIf ! k < model_levels
      End Do   !  k = 1, active_levels

!  Initialise diffusion coefficients to zero
      Do j = 1-off_y, rows+off_y
        Do i = 1-off_x, row_length+off_x
          diffc_i(i,j) = 0.0
          diffc_j(i,j) = 0.0
          w_max(i,j) = 0.0
        Enddo  ! i = 1-off_x, row_length+off_x
      Enddo    ! j = 1-off_y, rows+off_y

!  Check to see if w > w_limit at each point
       Do k =  test_level, end_level
         Do j = 1-off_y, rows+off_y
           Do i = 1-off_x, row_length+off_x
!  Find if vertical velocity above threshold at this point
!  start at level 5 since delta_z small near surface
              if( w_max(i,j) < w_adv(i,j,k) ) then
                  w_max(i,j) = w_adv(i,j,k)
              endif ! w_max(i,j) < w_adv(i,j,k
            Enddo  ! i = 1-off_x, row_length+off_x
          Enddo    ! j = 1-off_y, rows+off_y
        EndDo  !   k =  test_level, end_level

         Do j = 1-off_y, rows+off_y
           Do i = 1-off_x, row_length+off_x
             if( w_max(i,j) > w_limit)then
               diffc_i(i,j) = factor
               diffc_j(i,j) = factor
             endif ! w_max(i,j) > w_limit
           Enddo  ! i = 1-off_x, row_length+off_x
         Enddo    ! j = 1-off_y, rows+off_y

        Do j = 1, rows+off_y
          Do i = 1, row_length+off_x
           if( w_max(i,j) > w_max(i-1,j))then
             diffc_i(i-1,j) = diffc_i(i,j)
           endif ! w_max(i,j) > w_max(i-1,j)
           if( w_max(i,j) > w_max(i,j-1))then
             diffc_j(i,j-1) = diffc_j(i,j)
           endif ! w_max(i,j) > w_max(i,j-1)
          Enddo  !i = 1, row_length+off_x
        Enddo    !j = 1, rows+off_y

      If (L_diag_w) Then
        Do j = 1,rows
          Do i = 1, row_length
            if( w_max(i,j) > w_limit)then
                w_local_mask(i,j) = 1.0
            endif ! w_max(i,j) > w_limit
          End Do
        End Do
      EndIf !  L_diag_w

! ----------------------------------------------------------------------
! Section 2.   Switch off theta diffusion at steep slopes
! ----------------------------------------------------------------------
!  Horizontal_level should be the first flat level and as such
!  does not need testing for a sloping surface.
!   i.e if level_flat < 2   then
!       diffusion behaves as original (along eta surfaces)
! To switch off slope test completetely then
! either hard-wire level_flat = 1 or ensure that the argument
!    horizontal level <= wet_model_levels - - model_levels + 1
!      level_flat = 1
       level_flat =  horizontal_level
       if (level_flat  <=  0) level_flat = 1

       Do k = 1, active_levels
        Do j = 1, rows
          Do i = 0, row_length
           mask_i(i,j,k) = diffc_i(i,j)
          Enddo
        Enddo
        Do j = 0, rows
          Do i = 1, row_length
           mask_j(i,j,k) = diffc_j(i,j)
          Enddo
        Enddo
      Enddo   !  k = 1, active_levels

      Do k = start_level, level_flat - 1
         ka = k - start_level + 1
        Do j = 1, rows
          Do i = 0, row_length
            if( r_theta_levels(i,j,k) < r_theta_levels(i+1,j,k-1) .or.  &
     &          r_theta_levels(i+1,j,k) < r_theta_levels(i,j,k-1)) then
              mask_i(i,j,ka) = 0.0
           endif
          Enddo
        Enddo
        Do j = 0, rows
          Do i = 1, row_length
            if( r_theta_levels(i,j,k) < r_theta_levels(i,j+1,k-1) .or.  &
     &          r_theta_levels(i,j+1,k) < r_theta_levels(i,j,k-1)) then
              mask_j(i,j,ka) = 0.0
           endif
          Enddo
        Enddo
      Enddo    ! ka = start_level, level_flat - 1

! ----------------------------------------------------------------------
! Section 2.1  Calculate lambda direction term.
! ----------------------------------------------------------------------

        Do k = 1, active_levels
           ka = k + start_level - 1
            Do j = 1, rows
              Do i = 1-off_x, row_length
                temp(i,j) = ( field(i+1,j,ka) * delta_z(i+1,j,k) -      &
     &                        field(i  ,j,ka) * delta_z(i ,j,k) ) *     &
     &                                              mask_i(i,j,k) *     &
     &             r_theta_levels(i,j,ka) * r_theta_levels(i,j,ka)
              End Do
              Do i = 1, row_length
                lambda_term(i,j,k) = temp(i,j) - temp(i-1,j)
              End Do
            End Do

! ----------------------------------------------------------------------
! Section 2.2  Calculate phi direction term.
! ----------------------------------------------------------------------

             Do j = j_start-1, j_stop
              Do i = 1, row_length
                temp(i,j) = ( field(i,j+1,ka) * delta_z(i,j+1,k) -      &
     &                        field(i,j,  ka) * delta_z(i,j,  k) ) *    &
     &                         mask_j(i,j,k) * cos_v_latitude(i,j) *    &
     &               r_theta_levels(i,j,ka) * r_theta_levels(i,j,ka)
              End Do
             End Do
             Do j = j_start, j_stop
              Do i = 1, row_length
                phi_term(i,j,k) = (temp(i,j) - temp(i,j-1)) *           &
     &                              sec_theta_latitude(i,j)
              End Do
             End Do

        If ( model_domain == mt_Global ) Then

              If(at_extremity(PSouth))then
                Do i = 1, row_length
                  l_s_poles(i,k) = ( field(i,2,ka) * delta_z(i,2,k) -   &
     &                               field(i,1,ka) * delta_z(i,1,k)) *  &
     &                           mask_j(i,1,k) * cos_v_latitude(i,1) *  &
     &               r_theta_levels(i,1,ka) * r_theta_levels(i,1,ka)
                End Do
              End If

              If(at_extremity(PNorth))then
                Do i = 1, row_length
                  l_n_poles(i,k) = (field(i,rows-1,ka) *                &
     &                                          delta_z(i,rows-1,k) -   &
     &                        field(i,rows,ka) * delta_z(i,rows,k)) *   &
     &                   mask_j(i,rows-1,k) * cos_v_latitude(i,rows-1) *&
     &            r_theta_levels(i,rows,ka) * r_theta_levels(i,rows,ka)
                End Do
              End If

            End If !model_domain == mt_Global

        EndDo ! k = 1, active_levels

        If ( model_domain == mt_Global ) Then

          If (at_extremity(PSouth)) Then
            Call gcg_rvecsumr(row_length, row_length, 1, active_levels, &
     &                        l_s_poles, proc_row_group, info, sum_pole)
            Do k = 1, active_levels
              sum_pole(k) = pole_term * sum_pole(k)
              Do i = 1, row_length
                phi_term(i,1,k) = sum_pole(k)
              End Do
            End Do   ! k = 1, active_levels
          End If  ! at_extremity(PSouth)

          If (at_extremity(PNorth)) Then
            Call gcg_rvecsumr(row_length, row_length, 1, active_levels, &
     &                        l_n_poles, proc_row_group, info, sum_pole)
            Do k = 1, active_levels
              sum_pole(k) = pole_term * sum_pole(k)
              Do i = 1, row_length
                phi_term(i,rows,k) = sum_pole(k)
              End Do
            End Do  ! k = 1, active_levels
          End If ! at_extremity(PNorth)

        End If !  model_domain == mt_Global

! ----------------------------------------------------------------------
! Section 3.   Diffusion for q
! ----------------------------------------------------------------------

       Do  k = 1, active_levels
         ka = k + start_level - 1
          Do j = 1, rows
            Do i = 1, row_length
              field(i,j,ka) = field(i,j,ka) +                           &
     &                        (lambda_term(i,j,k) * recip_delta_lam2 +  &
     &                            phi_term(i,j,k) * recip_delta_phi2) / &
     &               (r_theta_levels(i,j,ka) * r_theta_levels(i,j,ka) * &
     &                                                delta_z(i,j,k) )
            End Do
          End Do

      End Do   ! k = 1, active_levels

      return    ! End of routine
      END SUBROUTINE filttarq

