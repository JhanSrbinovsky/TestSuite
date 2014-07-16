
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine tardiff_q_w
      subroutine tardiff_q_w                                            &
     &                     (q, w_adv,                                   &
     &                      r_theta_levels, r_rho_levels,               &
     &                      sec_theta_latitude,                         &
     &                      cos_theta_latitude, cos_v_latitude,         &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      halo_i_star, halo_j_star,                   &
     &                      at_extremity, proc_row_group,               &
     &                      delta_lambda, delta_phi,                    &
     &                      timestep, rows, n_rows, row_length,         &
     &                      model_levels, wet_model_levels,             &
     &                      model_domain, global_row_length,            &
     &                      q_star, w_limit,                            &
     &                      factor, test_level,                         &
     &                      start_level, end_level,                     &
     &                      L_diag_w, w_local_mask)

! Purpose:
!          Calculates conservative horizontal diffusion increment to q
!          subject to w > w_limit
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
! 5.5   24/02/03  This deck introduced  based on h_cdiff_q
!                 Pass out information for diagnostic output
!                                                Terry Davies
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_diag_w         ! diagnostic control

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, n_rows                                                          &
                         ! number of v rows.
     &, model_levels                                                    &
                         ! number of model levels.
     &, wet_model_levels                                                &
                             ! number of model levels.
     &, halo_i_star                                                     &
                      ! Size of halo in i direction for star field.
     &, halo_j_star                                                     &
                      ! Size of halo in j direction for star field.
     &, halo_i                                                          &
                      ! Size of halo in i direction.
     &, halo_j                                                          &
                      ! Size of halo in j direction.
     &, off_x                                                           &
                      ! Size of small halo in i
     &, off_y                                                           &
                      ! Size of small halo in j.
     &, global_row_length                                               &
     &, proc_row_group

      Integer                                                           &
     &  model_domain                                                    &
                         ! holds integer code for model domain
     &, test_level                                                      &
     &, start_level                                                     &
     &, end_level

      Real                                                              &
     &  timestep                                                        &
     &, w_limit                                                         &
                 ! Vertical velocity test value
     &, factor   ! effective diffusion coefficient

      Real                                                              &
           ! vertical co-ordinate arrays.
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)

      Real                                                              &
           !  trigonometric functions
     &  sec_theta_latitude (1-off_x:row_length+off_x,                   &
     &                    1-off_y:rows+off_y)                           &
     &, cos_theta_latitude (1-off_x:row_length+off_x,                   &
     &                    1-off_y:rows+off_y)                           &
     &, cos_v_latitude (1-off_x:row_length+off_x,                       &
     &                    1-off_y:n_rows+off_y)

      Real                                                              &
     &  q (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j             &
     &,         wet_model_levels)                                       &
     &, w_adv (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j         &
     &,        0:model_levels)

      Real                                                              &
           ! horizontal co-ordinate spacing.
     &  delta_lambda                                                    &
     &, delta_phi

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real                                                              &
     &  q_star (1-halo_i_star:row_length+halo_i_star,                   &
     &          1-halo_j_star:rows+halo_j_star, wet_model_levels)       &
     &, w_local_mask(row_length,rows)

! Local Variables.

      Integer                                                           &
     &  i, j, k, ka                                                     &
                                  ! Loop indices
     &, j_start, j_stop                                                 &
                              ! Loop indices
     &, info                                                            &
     &, active_levels

      Real                                                              &
     &  recip_delta_phi                                                 &
     &, scalar1                                                         &
     &, scalar2

! Local arrays

      Real                                                              &
     &  lambda_term(row_length, rows, model_levels )                    &
     &, phi_term(row_length, rows, model_levels )                       &
     &, delta_z(1-off_x:row_length+off_x,                               &
     &          1-off_y:rows+off_y, model_levels )                      &
     &, temp(1-off_x:row_length+off_x, 1-off_y:rows+off_y)              &
     &, l_s_poles(row_length, model_levels )                            &
     &, l_n_poles(row_length, model_levels )                            &
     &, sum_s(model_levels )                                            &
     &, sum_n(model_levels )                                            &
     &, diffc_i(1-off_x:row_length+off_x, 1-off_y:rows+off_y)           &
     &, diffc_j(1-off_x:row_length+off_x, 1-off_y:rows+off_y)           &
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

      active_levels = end_level - start_level + 1
      scalar1 = 1. / (delta_lambda * delta_lambda)
      recip_delta_phi = 1.0/ delta_phi
      scalar2 = recip_delta_phi * recip_delta_phi

      j_start = 1
      j_stop = rows
      If (model_domain  ==  1) Then
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
        If(k < model_levels) then
          Do j = 1-off_y, rows+off_y
            Do i = 1-off_x, row_length+off_x
               delta_z(i,j,k) = (r_rho_levels(i,j,k+1) -                &
     &                            r_rho_levels(i,j,k) )
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
          diffc_j(i,j) =   0.0
          w_max(i,j) = 0.0
        Enddo  ! i = 1-off_x, row_length+off_x
      Enddo    ! j = 1-off_y, rows+off_y

!  Check to see if w > w_print_limit at each point
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
               diffc_i(i,j) = factor / scalar1
               diffc_j(i,j) = factor / scalar2
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
! Section 2.1  Calculate lambda direction term.
! ----------------------------------------------------------------------

        Do k = 1, active_levels
           ka = k + start_level - 1
            Do j = 1, rows
              Do i = 1-off_x, row_length
                temp(i,j) = (q(i+1,j,ka) * delta_z(i+1,j,k) -           &
     &                        q(i  ,j,ka) * delta_z(i ,j,k) ) *         &
     &                                           diffc_i(i,j) *         &
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
                temp(i,j) = (q(i,j+1,ka) * delta_z(i,j+1,k) -           &
     &                        q(i,j,  ka) * delta_z(i,j,  k) ) *        &
     &                      diffc_j(i,j) * cos_v_latitude(i,j) *        &
     &               r_theta_levels(i,j,ka) * r_theta_levels(i,j,ka)
              End Do
             End Do
             Do j = j_start, j_stop
              Do i = 1, row_length
                phi_term(i,j,k) = (temp(i,j) - temp(i,j-1)) *           &
     &                              sec_theta_latitude(i,j)
              End Do
             End Do

        If ( model_domain  ==  mt_Global ) Then

              If(at_extremity(PSouth))then
                Do i = 1, row_length
             l_s_poles(i,k) = (q(i,2,ka) * delta_z(i,2,k) -             &
     &                          q(i,1,ka) * delta_z(i,1,k)) *           &
     &                  diffc_j(i,1) * cos_v_latitude(i,1) *            &
     &             r_theta_levels(i,1,ka) * r_theta_levels(i,1,ka)
                End Do
              End If

              If(at_extremity(PNorth))then
                Do i = 1, row_length
            l_n_poles(i,k) = (q(i,rows-1,ka) * delta_z(i,rows-1,k) -    &
     &                         q(i,rows,ka) * delta_z(i,rows,k)) *      &
     &                   diffc_j(i,rows-1) * cos_v_latitude(i,rows-1) * &
     &         r_theta_levels(i,rows-1,ka) * r_theta_levels(i,rows-1,ka)
                End Do
              End If

            End If !model_domain  ==  mt_Global

        EndDo ! k = 1, active_levels

        If ( model_domain  ==  mt_Global ) Then

          If (at_extremity(PSouth)) Then
            Call gcg_rvecsumr(row_length, row_length, 1, active_levels, &
     &                        l_s_poles, proc_row_group, info, sum_s)
            Do k = 1, active_levels
              sum_s(k) = sum_s(k) * 8. * recip_delta_phi /              &
     &                                   global_row_length
              Do i = 1, row_length
                phi_term(i,1,k) = sum_s(k)
              End Do
            End Do   ! k = 1, active_levels
          End If  ! at_extremity(PSouth)

          If (at_extremity(PNorth)) Then
            Call gcg_rvecsumr(row_length, row_length, 1, active_levels, &
     &                        l_n_poles, proc_row_group, info, sum_n)
            Do k = 1, active_levels
              sum_n(k) = sum_n(k) * 8. * recip_delta_phi /              &
     &                                   global_row_length
              Do i = 1, row_length
                phi_term(i,rows,k) = sum_n(k)
              End Do
            End Do  ! k = 1, active_levels
          End If ! at_extremity(PNorth)

        End If !model_domain  ==  mt_Global

! ----------------------------------------------------------------------
! Section 3.   Diffusion for q
! ----------------------------------------------------------------------

      Do k = 1, active_levels
         ka = k + start_level - 1
          Do j = 1, rows
            Do i = 1, row_length
              q_star(i,j,ka) = q_star(i,j,ka) +                         &
     &                        ( lambda_term(i,j,k) * scalar1 +          &
     &                             phi_term(i,j,k) * scalar2 ) /        &
     &             (r_theta_levels(i,j,ka) * r_theta_levels(i,j,ka) *   &
     &                                              delta_z(i,j,k) )
            End Do
          End Do
      End Do   ! k = 1, active_levels

      return    ! End of routine
      END SUBROUTINE tardiff_q_w


