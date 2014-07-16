
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine h_diff_q
      subroutine h_diff_q                                               &
     &                     (q, r_theta_levels,                          &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      halo_i_star, halo_j_star,                   &
     &                      at_extremity, proc_row_group,               &
     &                      n_proc, n_procx, n_procy, neighbour,        &
     &                      delta_lambda, delta_phi,                    &
     &                      timestep, rows, row_length,                 &
     &                      model_levels, wet_model_levels,             &
     &                      model_domain, global_row_length,            &
     &                      diffusion_coefficient,                      &
     &                      diffusion_order,                            &
     &                      q_star,                                     &
     &                      horizontal_level)

! Purpose:
!          Calculates horizontal diffusion increment to q.
!
! Method:
!          Is described in ;
!
!
! Original Programmer: Mark H. Mawson
! Current code owner: Andrew J. Malcolm
!
! History:
! Version   Date      Comment
! ----     ------     -------
!LL 5.2   24/10/00    This deck introduced               Andy Malcolm
!   5.3   19/10/01    Use appropriate gcg routines.   S. Cusack
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
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
     &, proc_row_group                                                  &
     &, n_proc                                                          &
                   ! Total number of processors
     &, n_procx                                                         &
                   ! Number of processors in longitude
     &, n_procy                                                         &
                   ! Number of processors in latitude
     &, neighbour(4)         ! Array with the Ids of the four neighbours
                             ! in the horizontal plane

      Integer                                                           &
     &  model_domain                                                    &
                         ! holds integer code for model domain
     &, diffusion_order(wet_model_levels)

      Real                                                              &
     &  timestep

      Real                                                              &
     &  diffusion_coefficient(wet_model_levels)

      Real                                                              &
           ! vertical co-ordinate arrays.
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)

      Real                                                              &
     &  q (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
     &         wet_model_levels)

      Real                                                              &
           ! horizontal co-ordinate spacing.
     &  delta_lambda                                                    &
     &, delta_phi

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real                                                              &
     &  q_star (1-halo_i_star:row_length+halo_i_star,                   &
     &          1-halo_j_star:rows+halo_j_star, wet_model_levels)

      Integer                                                           &
     &  horizontal_level   ! level at which steep slope test no
!                          ! longer operates

! Local Variables.

      Integer                                                           &
     &  i, j, k, ka, order, j0, j1                                      &
                                        ! Loop indices
     &, info

      Real                                                              &
     &  scalar1                                                         &
     &, scalar2                                                         &
     &, sign

! Local arrays

      Real                                                              &
     &  lambda_term(row_length, rows, model_levels)                     &
     &, phi_term(row_length, rows, model_levels)                        &
     &, timestep_over_r_squared(row_length, rows, model_levels)         &
     &, field(1-off_x:row_length+off_x,                                 &
     &        1-off_y:rows+off_y, model_levels)                         &
     &, l_s_poles(row_length, model_levels)                             &
     &, l_n_poles(row_length, model_levels)                             &
     &, sum_s(model_levels)                                             &
     &, sum_n(model_levels)

      Integer                                                           &

     &  mask_i(1-off_x:row_length, rows, model_levels)                  &
     &, mask_j(row_length, 1-off_y:rows, model_levels)

      Integer                                                           &
     &  active_diffusion_order(model_levels)                            &
     &, n_active_levels                                                 &
     &, max_diffusion_order                                             &
     &, level_base                                                      &
     &, level_flat

      Real                                                              &
     &  active_diff_coeff(model_levels)

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

      scalar1 = 1. / (delta_lambda * delta_lambda)
      scalar2 = 1. / (delta_phi * delta_phi)
      j0 = 1
      j1 = rows
      If (model_domain  ==  1) Then
        If (at_extremity(PSouth)) j0 = 2
        If (at_extremity(PNorth)) j1 = rows-1
      End If

! ----------------------------------------------------------------------
! Section 1.   Copy q into field for those model levels with a
!              positive diffusion coefficient.
!              Copy coeffs and diffusion order for all active levels.
! ----------------------------------------------------------------------

      n_active_levels = 0
      max_diffusion_order = 0
      Do k = 1, wet_model_levels
        If (diffusion_coefficient(k)  >   0.0 .and.                     &
     &      diffusion_order(k)  >   0 )Then
          n_active_levels = n_active_levels + 1
          active_diffusion_order(n_active_levels) = diffusion_order(k)
          active_diff_coeff(n_active_levels) = diffusion_coefficient(k)
          Do j = 1-off_y, rows+off_y
            Do i = 1-off_x, row_length+off_x
              field(i,j,n_active_levels) = q(i,j,k)
            End Do
          End Do
          Do j = 1, rows
            Do i = 1, row_length
              timestep_over_r_squared(i,j,n_active_levels) =            &
     &                        timestep /                                &
     &                       ( r_theta_levels(i,j,k) *                  &
     &                         r_theta_levels(i,j,k) )
            End Do
          End Do
          If (diffusion_order(k)  >   max_diffusion_order )             &
     &        max_diffusion_order = diffusion_order(k)
        End If
      End Do

        level_base = wet_model_levels - n_active_levels
!  Last model-level upon which diffusion is inactive
! ----------------------------------------------------------------------
! Section 2.   Loop over the order to produce order * del squared
! ----------------------------------------------------------------------
!  Horizontal_level should be the first flat level and as such
!  does not need testing for a sloping surface.
!   i.e if level_flat < 2   then
!       diffusion behaves as original (along eta surfaces)
! To switch off slope test completetely then
! either hard-wire level_flat = 1 or ensure that the argument
!    horizontal level <= wet_model_levels - n_active_levels + 1
!      level_flat = 1
       level_flat = n_active_levels -                                   &
     &              wet_model_levels + horizontal_level
       if (level_flat  <=  0) level_flat = 1

       Do k = 1, n_active_levels
        Do j = 1, rows
          Do i = 0, row_length
           mask_i(i,j,k)=1
          Enddo
        Enddo
        Do j = 0, rows
          Do i = 1, row_length
           mask_j(i,j,k)=1
          Enddo
        Enddo
      Enddo

      Do ka = 1, level_flat - 1
        k = ka + level_base
        Do j = 1, rows
          Do i = 0, row_length
           if(r_theta_levels(i+1,j,k) <  r_theta_levels(i,j,k+1)        &
     &        .and.                                                     &
     &        r_theta_levels(i+1,j,k) >  r_theta_levels(i,j,k-1))then
           mask_i(i,j,ka)=1
           else
           mask_i(i,j,ka)=0
           endif
          Enddo
        Enddo
        Do j = 0, rows
          Do i = 1, row_length
           if(r_theta_levels(i,j+1,k) <  r_theta_levels(i,j,k+1)        &
     &        .and.                                                     &
     &        r_theta_levels(i,j+1,k) >  r_theta_levels(i,j,k-1))then
           mask_j(i,j,ka)=1
           else
           mask_j(i,j,ka)=0
           endif
          Enddo
        Enddo
      Enddo

      Do order = 1, max_diffusion_order

! ----------------------------------------------------------------------
! Section 2.1  Calculate lambda direction term.
! ----------------------------------------------------------------------

        Do k = 1, n_active_levels
! Only do calculations for levels whose diffusion order is less than
! the current value
          If (order  <=  active_diffusion_order(k)) Then
            Do j = 1, rows
              Do i = 1, row_length
                lambda_term(i,j,k) = (field(i+1,j,k) - field(i,j,k)) *  &
     &                            mask_i(i,j,k)*active_diff_coeff(k)    &
     &                         - (field(i,j,k) - field(i-1,j,k)) *      &
     &                            mask_i(i-1,j,k)*active_diff_coeff(k)
              End Do
            End Do

! ----------------------------------------------------------------------
! Section 2.2  Calculate phi direction term.
! ----------------------------------------------------------------------

            Do j = j0, j1
              Do i = 1, row_length
                phi_term(i,j,k) = (field(i,j+1,k) - field(i,j,k)) *     &
     &                         mask_j(i,j,k)*active_diff_coeff(k)       &
     &                      - (field(i,j,k) - field(i,j-1,k)) *         &
     &                         mask_j(i,j-1,k)*active_diff_coeff(k)
              End Do
            End Do

            If (model_domain  ==  1) Then
              If(at_extremity(PSouth))then
                Do i = 1, row_length
                  l_s_poles(i,k) = (field(i,2,k) - field(i,1,k)) *      &
     &                             active_diff_coeff(k)
                End Do
              End If
              If(at_extremity(PNorth))then
                Do i = 1, row_length
                  l_n_poles(i,k)=(field(i,rows-1,k) - field(i,rows,k))* &
     &                           active_diff_coeff(k)
                End Do
              End If
            End If

          End If ! On if (order  <=  active_diffusion_order(k))
        End Do ! end loop over model levels

        If (model_domain  ==  1) Then
          If (at_extremity(PSouth)) Then
            Call gcg_rvecsumr(row_length, row_length, 1,                &
     &                        n_active_levels,                          &
     &                        l_s_poles, proc_row_group, info, sum_s)
            Do k = 1, n_active_levels
              sum_s(k) = sum_s(k) * 2. / global_row_length
              Do i = 1, row_length
                phi_term(i,1,k) = sum_s(k)
              End Do
            End Do
          End If

          If (at_extremity(PNorth)) Then
            Call gcg_rvecsumr(row_length, row_length, 1,                &
     &                        n_active_levels,                          &
     &                        l_n_poles, proc_row_group, info, sum_n)
            Do k = 1, n_active_levels
              sum_n(k) = sum_n(k) * 2. / global_row_length
              Do i = 1, row_length
                phi_term(i,rows,k) = sum_n(k)
              End Do
            End Do
          End If
        End If

! ----------------------------------------------------------------------
! Section 2.3   Calculate new variable.
! ----------------------------------------------------------------------

        Do k = 1, n_active_levels
          If (order  <=  active_diffusion_order(k)) Then
            Do j = 1, rows
              Do i = 1, row_length
                field(i,j,k) = timestep_over_r_squared(i,j,k) *         &
     &                       (lambda_term(i,j,k) * scalar1 +            &
     &                        phi_term(i,j,k) * scalar2 )
              End Do
            End Do
          End If
        End Do

        If (order  /=  max_diffusion_order ) then
! DEPENDS ON: swap_bounds
          call Swap_Bounds(                                             &
     &                   field, row_length, rows, n_active_levels,      &
     &                   off_x, off_y, fld_type_p, .false.)
        End If

! End loop over order
      End Do

! ----------------------------------------------------------------------
! Section 3.   Diffusion increment is field * appropriate sign.
! ----------------------------------------------------------------------

! use ka to calculate which level of field maps to which level of
! input data.
      ka = 0
      Do k = 1, model_levels
        If (diffusion_coefficient(k)  >   0.0 .and.                     &
     &      diffusion_order(k)  >   0 )Then
          ka = ka + 1
          sign = (-1) ** (diffusion_order(k)-1)
          Do j = 1, rows
            Do i = 1, row_length
              q_star(i,j,k) = q_star(i,j,k) +                           &
     &                            field(i,j,ka) * sign
            End Do
          End Do
        End If
      End Do

! End of routine
      return
      END SUBROUTINE h_diff_q

