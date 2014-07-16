
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine turb_diff_v

      Subroutine turb_diff_v(                                           &
     &                       v, r_at_v,                                 &
     &                       r_theta_levels, r_rho_levels,              &
     &                       sec_v_latitude,                            &
     &                       cos_theta_latitude, sec_theta_latitude,    &
     &                       off_x, off_y, halo_i, halo_j,              &
     &                       delta_lambda, delta_phi, timestep,         &
     &                       rows, n_rows, row_length,                  &
     &                       model_levels, levels,                      &
     &                       diff_coeff_theta, diff_coeff_centre, R_v )

! Purpose:
!          Turbulent diffusion of a v-type field
!          Based on conservative diffusion operator.
!          3d diffusion coefficient
!          Active wherever diff_coeff > 0
!          If you need to switch off over steep slopes
!          then set diff_coeff = 0 at required points before entry.
!
! Method:
!          Is described in ;
!
!
! Original Programmer:   T Davies
! Current code owner: Andrew J. Malcolm
!
! History:
! Version    Date      Comment
! ----     -------     -------
! 6.2     04/10/05    New code based on conservative diffusion. T Davies
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, n_rows                                                          &
                         ! number of v-rows
     &, model_levels                                                    &
                         ! number of model levels.
     &, halo_i                                                          &
                      ! Size of halo in i direction.
     &, halo_j                                                          &
                      ! Size of halo in j direction.
     &, off_x                                                           &
                      ! Size of small halo in i
     &, off_y                                                           &
                      ! Size of small halo in j.
     &, levels   ! levels to process

      Real                                                              &
     &  delta_lambda                                                    &
     &, delta_phi                                                       &
     &, timestep

      Real                                                              &
           ! vertical co-ordinate arrays and input arrays
     &  r_at_v (1-halo_i:row_length+halo_i,                             &
     &          1-halo_j:n_rows+halo_j, model_levels)                   &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)               &
     &, r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, v (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,              &
     &     model_levels)

      Real                                                              &
           !  trigonometric functions
     &  sec_v_latitude (1-off_x:row_length+off_x,                       &
     &                  1-off_y:n_rows+off_y)                           &
     &, cos_theta_latitude (1-off_x:row_length+off_x,                   &
     &                      1-off_y:rows+off_y)                         &
     &, sec_theta_latitude (1-off_x:row_length+off_x,                   &
     &                      1-off_y:rows+off_y)                         &
!  diffusion coefficients common to both u and v
!  needed at theta points and at centre of grid cell
     &, diff_coeff_theta(1-off_x:row_length+off_x, 1-off_y:rows + off_y,&
     &                                                   levels)        &
     &, diff_coeff_centre(1-off_x:row_length, 1-off_y:rows, levels)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real                                                              &
     &  R_v (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,            &
     &       model_levels)

! Local Variables.

      Integer                                                           &
     &  i, j, k           ! Loop indices

      Real                                                              &
     &  recip_delta_lambda                                              &
     &, recip_delta_phi

! Local arrays

      Real                                                              &
     &  r_theta_at_v(1-off_x:row_length+off_x,                          &
     &               1-off_y:n_rows+off_y, 0:levels)                    &
     &, delta_z(1-off_x:row_length+off_x,                               &
     &          1-off_y:n_rows+off_y, levels )                          &
     &, recip_r_squared_delz(1-off_x:row_length+off_x,                  &
     &                       1-off_y:n_rows+off_y, levels )             &
     &, temp(1-off_x:row_length, n_rows+ 1)                             &
     &, diff_lambda(1-off_x:row_length, n_rows)                         &
     &, lambda_term(row_length, n_rows)                                 &
     &, phi_term(row_length, n_rows)

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

! ----------------------------------------------------------------------
! Section 1.   Set r field  and calculate delta_z
! ----------------------------------------------------------------------

      recip_delta_lambda = 1.0 / delta_lambda
      recip_delta_phi = 1.0 / delta_phi

      Do k = 0, levels
        Do j = 1-off_y, n_rows+off_y
          Do i = 1-off_x, row_length+off_x
            r_theta_at_v(i,j,k) = .5 * (r_theta_levels(i+1,j,k) +       &
     &                                  r_theta_levels(i  ,j,k) )
          End Do
        End Do
      End Do ! k = 0, levels

! call swap_bounds to set halo points

! DEPENDS ON: swap_bounds
      call Swap_Bounds                                                  &
     &                (r_theta_at_v,                                    &
     &                 row_length, n_rows, levels + 1,                  &
     &                 off_x, off_y, fld_type_v, .false.)

! calculate D(r)/D(eta) about theta levels
      Do k = 1, levels
        if ( k < model_levels ) then
          Do j = 1-off_y, n_rows+off_y
            Do i = 1-off_x, row_length+off_x
!
! "conservative diffusion" diffuses dzQ rather than Q 
! delta_z is set to 1.0 to avoid instabilities ocurring in the presence
! of orography due to volume differences between neighbouring gridboxes.
!
              delta_z(i,j,k) = 1.0
!              delta_z(i,j,k) = r_theta_at_v(i,j,k) -                   &
!     &                           r_theta_at_v(i,j,k-1)
              recip_r_squared_delz(i,j,k) = 1.0 / ( r_at_v(i,j,k) *     &
     &                                               r_at_v(i,j,k) *    &
     &                                               delta_z(i,j,k) )
            End Do
          End Do
        else ! ( k = model_levels
! By setting dr/deta=1 we effectively cancel out constant dr/deta
! at top level and then no need to take special action in later loops
          Do j = 1-off_y, n_rows+off_y
            Do i = 1-off_x, row_length+off_x
              delta_z(i,j,k) = 1.0
              recip_r_squared_delz(i,j,k) = 1.0 / ( r_at_v(i,j,k) *     &
     &                                               r_at_v(i,j,k) )
            End Do
          End Do
        endif ! k < model_levels
      End Do   !  k = 1, levels

! ----------------------------------------------------------------------
! Section 2.0 Diffusion
! ----------------------------------------------------------------------

      Do k = 1, levels

! ----------------------------------------------------------------------
! Section 2.1  Calculate lambda direction term.
! ----------------------------------------------------------------------

        Do j = 1, n_rows
          Do i = 1-off_x, row_length
            temp(i,j) = (v(i+1,j,k) * delta_z(i+1,j,k) -                &
     &                    v(i  ,j,k) * delta_z(i  ,j,k) ) *             &
     &                recip_delta_lambda * diff_coeff_centre(i,j,k) *   &
     &              sec_theta_latitude(i,j) * sec_theta_latitude(i,j)
          End Do
          Do i = 1, row_length
            lambda_term(i,j) = recip_delta_lambda *                     &
     &                           ( temp(i,j) - temp(i-1,j) )
          End Do
        End Do

! ----------------------------------------------------------------------
! Section 2.2  Calculate phi direction term.
! ----------------------------------------------------------------------

        Do j = 1, n_rows + 1
          Do i = 1, row_length
            temp(i,j) = (v(i,j,k) * delta_z(i,j,k) -                    &
     &                    v(i,j-1,k) * delta_z(i,j-1,k) ) *             &
     &                   recip_delta_phi * diff_coeff_theta(i,j,k) *    &
     &                                   cos_theta_latitude(i,j)
          End Do
        End Do
        Do j = 1, n_rows
          Do i = 1, row_length
            phi_term(i,j) = ( temp(i,j+1) - temp(i,j) ) *               &
     &                       recip_delta_phi * sec_v_latitude(i,j)
          End Do
        End Do

! ----------------------------------------------------------------------
! Section 2.3   Calculate new variable.
! ----------------------------------------------------------------------

        Do j = 1, n_rows
          Do i = 1, row_length
            R_v(i,j,k) = R_v(i,j,k) + timestep *                        &
     &                       recip_r_squared_delz(i,j,k) *              &
     &                     ( lambda_term(i,j) + phi_term(i,j) )
          End Do
        End Do

      End Do  ! k = 1, levels

! DEPENDS ON: swap_bounds
      call Swap_Bounds(                                                 &
     &                  R_v, row_length, n_rows, levels,                &
     &                  off_x, off_y, fld_type_v, .true.)

      return
      END SUBROUTINE turb_diff_v

