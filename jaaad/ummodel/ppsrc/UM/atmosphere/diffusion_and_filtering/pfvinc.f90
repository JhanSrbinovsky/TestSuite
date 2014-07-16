
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine pfvinc

      Subroutine pfvinc(                                                &
     &                  R_v, r_at_v, r_theta_levels, r_rho_levels,      &
     &                  off_x, off_y, halo_i, halo_j,                   &
     &                  sec_v_latitude, cos_theta_latitude,             &
     &                  me, n_procy, at_extremity,                      &
     &                  delta_lambda, delta_phi,                        &
     &                  rows, n_rows, row_length, model_levels,         &
     &                  max_filter_rows, v_begin, v_end,                &
     &                  v_sweeps, global_v_filter,                      &
     &                  horizontal_level, model_domain,                 &
     &                  diff_coeff_phi, diff_coeff_v, L_diff )

! Purpose:
!          Filter/diffusion of a v-type field
!          Based on conservative diffusion operator.
!          Applied as multiple EW 1-2-1 filter near poles.
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
! 6.2     14/02/06    New code based on pofil_v_incs + changes
!                                                         Terry Davies
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
     &, L_diff

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
     &, model_domain                                                    &
     &, n_procy                                                         &
                   ! Number of processors in latitude
     &, me        ! processor id

      Integer                                                           &
     &  max_filter_rows                                                 &
                             !  array dimension for v_begin etc.
     &, v_begin(0:max_filter_rows)                                      &
     &, v_end(0:max_filter_rows)                                        &
     &, v_sweeps(max_filter_rows)                                       &
     &, global_v_filter                                                 &
     &, horizontal_level   ! level at which steep slope test no
!                               ! longer operates

      Real                                                              &
     &  delta_lambda                                                    &
     &, delta_phi                                                       &
     &, diff_coeff_phi   ! NS diffusion coefficient

      Real                                                              &
           ! vertical co-ordinate arrays.
     &  r_at_v (1-halo_i:row_length+halo_i,                             &
     &          1-halo_j:n_rows+halo_j, model_levels)                   &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)               &
     &, r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)

      Real                                                              &
           !  trigonometric functions
     &  sec_v_latitude (1-off_x:row_length+off_x,                       &
     &                  1-off_y:n_rows+off_y)                           &
     &, cos_theta_latitude (1-off_x:row_length+off_x,                   &
     &                      1-off_y:rows+off_y)                         &
! EW diffusion coefficient
     &, diff_coeff_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real                                                              &
     &  R_v(1-off_x:row_length+off_x,                                   &
     &      1-off_y:n_rows+off_y, model_levels)

! Local Variables.

      Integer                                                           &
     &  i, j, k                                                         &
                          ! Loop indices
     &, j_begin, j_end                                                  &
                          ! Loop bounds
     &, j_store                                                         &
     &, level_flat                                                      &
     &, i_filter                                                        &
     &, i_sweep                                                         &
     &, num_pass

      Logical                                                           &
     &  L_cycle                                                         &
     &, L_combine

! Local arrays

      Real                                                              &
     &  r_theta_at_v(1-off_x:row_length+off_x,                          &
     &               1-off_y:n_rows+off_y, 0:model_levels)              &
     &, delta_z(1-off_x:row_length+off_x,                               &
     &          1-off_y:n_rows+off_y, model_levels )                    &
     &, recip_r_squared_delz(1-off_x:row_length+off_x,                  &
     &                       1-off_y:n_rows+off_y, model_levels )       &
     &, temp(1-off_x:row_length, 1-off_y:n_rows)                        &
     &, mask_i(1-off_x:row_length, n_rows, model_levels)                &
     &, mask_j(row_length, 1-off_y:n_rows, model_levels)

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
! Section 1.   Set up loop bounds and pointers
!            v_begin and v_end for each potential sweep are set SETCON
! ----------------------------------------------------------------------
! To switch off slope test completetely then
! either hard-wire level_flat = 1 or ensure that the argument
!    horizontal level <= wet_model_levels - - model_levels + 1
!      level_flat = 1
      level_flat =  horizontal_level
      if (level_flat  <=  0) level_flat = 1

! ----------------------------------------------------------------------
! Section 2.   Set r field  and calculate delta_z
! ----------------------------------------------------------------------

      Do k = 0, model_levels
        Do j = 1-off_y, n_rows+off_y
          Do i = 1-off_x, row_length+off_x
            r_theta_at_v(i,j,k) = .5 * (r_theta_levels(i,j+1,k) +       &
     &                                  r_theta_levels(i,j  ,k) )
          End Do
        End Do
      End Do ! k = 0, model_levels

! call swap_bounds to set halo points

! DEPENDS ON: swap_bounds
      call Swap_Bounds                                                  &
     &                (r_theta_at_v,                                    &
     &                 row_length, n_rows, model_levels + 1,            &
     &                 off_x, off_y, fld_type_v, .false.)

! calculate D(r)/D(eta) about theta levels
      Do k = 1, model_levels - 1
        Do j = 1-off_y, n_rows+off_y
          Do i = 1-off_x, row_length+off_x
            delta_z(i,j,k) = r_theta_at_v(i,j,k) -                      &
     &                         r_theta_at_v(i,j,k-1)
            recip_r_squared_delz(i,j,k) = 1.0 /                         &
     &                                 ( r_at_v(i,j,k) * r_at_v(i,j,k) *&
     &                                                  delta_z(i,j,k) )
            End Do
          End Do
      End Do   !  k = 1, model_levels - 1
! By setting dr/deta=1 we effectively cancel out constant dr/deta
! at top level and then no need to take special action in later loops
      k = model_levels
      Do j = 1-off_y, n_rows+off_y
        Do i = 1-off_x, row_length+off_x
          delta_z(i,j,k) = 1.0
          recip_r_squared_delz(i,j,k) = 1.0 /                           &
     &                                 ( r_at_v(i,j,k) * r_at_v(i,j,k) )
        End Do
      End Do

! ----------------------------------------------------------------------
! Section 3.   Switch off theta diffusion at steep slopes
! ----------------------------------------------------------------------
!  Horizontal_level should be the first flat level and as such
!  does not need testing for a sloping surface.
!   i.e if level_flat < 2   then
!       diffusion behaves as original (along eta surfaces)

      Do k = 1, model_levels
        Do j = 1, n_rows
          Do i = 1-off_x, row_length
            mask_i(i,j,k) = 1.0
          Enddo
        Enddo
        Do j = 0, n_rows
          Do i = 1, row_length
            mask_j(i,j,k) = 1.0
          Enddo
        Enddo
      Enddo  ! 1, model_levels

      Do j = 1, n_rows
        Do i = 1-off_x, row_length
          if( r_at_v(i,j,1) < r_theta_at_v(i+1,j,0) .or.                &
     &         r_at_v(i+1,j,1) < r_theta_at_v(i,j,0) )                  &
     &        mask_i(i,j,1) = 0.0
        Enddo
      Enddo
      Do j = 0, n_rows
        Do i = 1, row_length
          if( r_at_v(i,j,1) < r_theta_at_v(i,j+1,0) .or.                &
     &         r_at_v(i,j+1,1) < r_theta_at_v(i,j,0) )                  &
     &        mask_j(i,j,1) = 0.0
        Enddo
      Enddo
      Do k = 2, level_flat - 1
        Do j = 1, n_rows
          Do i = 1-off_x, row_length
            if( r_at_v(i,j,k) < r_at_v(i+1,j,k-1) .or.                  &
     &           r_at_v(i+1,j,k) < r_at_v(i,j,k-1) )                    &
     &          mask_i(i,j,k) = 0.0
          Enddo
        Enddo
        Do j = 0, n_rows
           Do i = 1, row_length
            if( r_at_v(i,j,k) < r_at_v(i,j+1,k-1) .or.                  &
     &           r_at_v(i,j+1,k) < r_at_v(i,j,k-1) )                    &
     &         mask_j(i,j,k) = 0.0
          Enddo
        Enddo
      Enddo  ! k = 2, level_flat - 1

! ----------------------------------------------------------------------
! Section 4.0  Polar Filtering and diffusion
! ----------------------------------------------------------------------

!        write(6,*)'pofil_v_incs global_v_filter = ', global_v_filter

      L_combine = .false.
      num_pass = global_v_filter
      j_begin = v_begin(1)
      j_end = v_end(1)
! Combine EW diffusion with first sweep of polar filter
      if ( L_diff ) then
        L_combine = .true.
        if( global_v_filter < 1 ) num_pass = 1
        if ( j_begin < 0 ) then
          j_begin = v_begin(0)
          j_end = v_end(0)
        elseif ( j_end < v_end(0) ) then
          j_end = v_end(0)
        elseif ( j_begin > v_begin(0) ) then
          j_begin = v_begin(0)
        endif !  j_begin < 0
      endif ! L_diff

      DO i_filter = 1, num_pass

        if( i_filter > 1 ) then
          L_combine = .false.
          j_begin = v_begin(i_filter)
          j_end = v_end(i_filter)
        endif !i_filter > 1

        L_cycle = .false.
!  If n_procy = 1, need to filter hemispheres separately
        if( n_procy == 1 ) L_cycle = .true.

        i_sweep = 1

        Do  ! Sweeping loop  i_sweep = 1, v_sweeps(i_filter)
!          write(6,*)'  i_filter = ', i_filter ,' i_sweep = ', i_sweep
!     &          ,' j_begin = ', j_begin,' j_end = ', j_end

! ----------------------------------------------------------------------
! Section 4.1   EW Filtering and diffusion
! ----------------------------------------------------------------------
          Do k = 1, model_levels

            Do j = j_begin, j_end
              Do i = 1-off_x, row_length
                temp(i,j) = (R_v(i+1,j,k) * delta_z(i+1,j,k) -          &
     &                         R_v(i  ,j,k) * delta_z(i  ,j,k) ) *      &
     &                           mask_i(i,j,k) * diff_coeff_v(i,j) *    &
     &                  0.25 * ( r_at_v(i,j,k) + r_at_v(i+1,j,k) ) *    &
     &                         ( r_at_v(i,j,k) + r_at_v(i+1,j,k) )
              End Do

              Do i = 1, row_length
                R_v(i,j,k) = R_v(i,j,k) + ( temp(i,j) - temp(i-1,j) ) * &
     &                                   recip_r_squared_delz(i,j,k)
              End Do
           End Do  ! j = j_begin, j_end

          End Do  ! k = 1, model_levels

! If n_procy=1, Northern hemisphere needs to be done after S. Hem
          If ( L_cycle ) then
            j_store = j_begin
            if ( i_sweep == 1 .and. L_combine ) then
              j_begin = j_end + 1
            else
              j_begin = n_rows - j_end + 1
            endif ! i_sweep == 1 .and. L_combine
            j_end = n_rows - j_store + 1
            L_cycle = .false.
        CYCLE
          endif ! L_cycle)

! Reset pointers for next sweep
! either because N Hem has just been done or
! 1st sweep was combined filter and diffusion
          j_begin = v_begin(i_filter)
          j_end = v_end(i_filter)
          if( n_procy == 1 ) L_cycle = .true.

! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                     R_v, row_length, n_rows, model_levels,       &
     &                     off_x, off_y, fld_type_v, .true.)

          i_sweep = i_sweep + 1
          if ( i_sweep > v_sweeps(i_filter) ) EXIT
        CYCLE
        EndDo ! sweeping  loop

      EndDO   !  i_filter = 1, num_pass

! ----------------------------------------------------------------------
! Section 4.2   NS diffusion
! ----------------------------------------------------------------------

      if( L_diff ) then
        j_begin = v_begin(0)
        j_end = v_end(0)
        if( n_procy == 1 ) j_end = n_rows - j_begin + 1
!          write(6,*)' NS diffusion  diff_coeff_phi = ', diff_coeff_phi
!     &          ,' j_begin = ', j_begin,' j_end = ', j_end

        Do k = 1, model_levels

          Do j = j_begin-off_y, j_end
            Do i = 1, row_length
              temp(i,j) = (R_v(i,j+1,k) * delta_z(i,j+1,k) -            &
     &                       R_v(i,j,  k) * delta_z(i,j, k) ) *         &
     &                               mask_j(i,j,k) * diff_coeff_phi *   &
     &                r_rho_levels(i,j+1,k) * r_rho_levels(i,j+1,k) *   &
     &                                  cos_theta_latitude(i,j+1)
            End Do
          End Do
          Do j = j_begin, j_end
            Do i = 1, row_length
              R_v(i,j,k) = R_v(i,j,k) + sec_v_latitude(i,j) *           &
     &                           ( temp(i,j) - temp(i,j-1) ) *          &
     &                         recip_r_squared_delz(i,j,k)
            End Do
          End Do

          If ( model_domain == mt_global ) Then
            If (at_extremity(PSouth) ) Then
              j = j_begin-1
              Do i = 1, row_length
                R_v(i,j,k) = R_v(i,j,k) + recip_r_squared_delz(i,j,k) * &
     &                               temp(i,j) * sec_v_latitude(i,j)
              End Do
            endIf !at_extremity(PSouth)
          If (at_extremity(PNorth) ) Then
            j = j_end + 1
              Do i = 1, row_length
                R_v(i,j,k) = R_v(i,j,k) - recip_r_squared_delz(i,j,k) * &
     &                             temp(i,j-1) * sec_v_latitude(i,j)
              End Do
            endIf !at_extremity(PNorth)
          endIf ! model_domain == mt_global )

        End Do  ! k = 1, model_levels

! DEPENDS ON: swap_bounds
        call Swap_Bounds(                                               &
     &                   R_v, row_length, n_rows, model_levels,         &
     &                   off_x, off_y, fld_type_v, .true.)
      endif ! L_diff


      return
      END SUBROUTINE pfvinc

