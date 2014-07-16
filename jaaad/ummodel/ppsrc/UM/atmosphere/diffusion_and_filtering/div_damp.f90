
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Div_damp

      Subroutine Div_damp(                                              &
     &                      u, v, rho, r_at_u, r_at_v,                  &
     &                      r_rho_levels,                               &
     &                      FV_sec_theta_latitude, cos_v_latitude,      &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      me, at_extremity, proc_row_group,           &
     &                      n_proc, n_procx, n_procy, neighbour,        &
     &                      delta_lambda, delta_phi,                    &
     &                      timestep, rows, n_rows, row_length,         &
     &                      model_levels, model_domain,                 &
     &                      g_row_length, div_damp_coefficient,         &
     &                      R_u, R_v)

! Purpose:
!          Performs divergence damping.
!
! Method:
!          Is described in ;
!
!
! Original Programmer: Mark H. Mawson
! Current code owner: Andrew J. Malcolm
!
! History:
! Date     Version     Comment
! ----     -------     -------
!LL   5.1   11/02/00  Use DOMTYP parameters                    P.Burton
!LL   5.2   7/11/00   Correct fld_type argument in SWAPBOUNDS call
!LL                                                       P.Burton
!     5.3  19/10/01   Use appropriate gcg routines.   S. Cusack
!   5.3     15/09/01  add mt_bi_cyclic_LAM code            A. Malcolm
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
     &, n_rows                                                          &
                         ! number of rows.
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
     &, proc_row_group                                                  &
     &, n_proc                                                          &
                   ! Total number of processors
     &, n_procx                                                         &
                   ! Number of processors in longitude
     &, n_procy                                                         &
                   ! Number of processors in latitude
     &, neighbour(4)                                                    &
                             ! Array with the Ids of the four neighbours
                             ! in the horizontal plane
     &, me                                                              &
     &, g_row_length

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

      Real                                                              &
     &  timestep

      Real                                                              &
     &  div_damp_coefficient(model_levels)

      Real                                                              &
           ! vertical co-ordinate arrays.
     &  r_at_u (1-halo_i:row_length+halo_i,                             &
     &          1-halo_j:rows+halo_j, model_levels)                     &
     &, r_at_v (1-halo_i:row_length+halo_i,                             &
     &          1-halo_j:n_rows+halo_j, model_levels)                   &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)

      Real                                                              &
     &  u (1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     &     model_levels)                                                &
     &, v (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,              &
     &     model_levels)                                                &
     &, rho (1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
     &       model_levels)

      Real                                                              &
           ! trigonometric arrays.
     &  FV_sec_theta_latitude (1-off_x:row_length+off_x,                &
     &                         1-off_y:rows+off_y)                      &
     &, cos_v_latitude (1-off_x:row_length+off_x,                       &
     &                  1-off_y:n_rows+off_y)

      Real                                                              &
           ! horizontal co-ordinate spacing.
     &  delta_lambda                                                    &
     &, delta_phi

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real                                                              &
     &  R_u (1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
     &       model_levels)                                              &
     &, R_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,             &
     &      model_levels)

! Local Variables.

      Integer                                                           &
     &  i, j, k, j0, j1                                                 &
                             ! Loop indices
     &, info

      Real                                                              &
     &  recip_delta_lambda                                              &
     &, recip_delta_phi                                                 &
     &, sum_n                                                           &
     &, sum_s

! Local arrays

      Real                                                              &
     &  divergence(1-off_x:row_length+off_x, 1-off_y:rows+off_y)        &
     &, l_n_poles(row_length)                                           &
     &, l_s_poles(row_length)

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

! External Routines:

! ----------------------------------------------------------------------
! Section 1. Calculate divergence
! ----------------------------------------------------------------------

      j0 = 1
      j1 = rows
      If (model_domain  /=  mt_bi_cyclic_lam) then
      If (at_extremity(PSouth)) j0 = 2
      If (at_extremity(PNorth)) j1 = rows-1
      Endif
      recip_delta_lambda = 1. / delta_lambda
      recip_delta_phi = 1. / delta_phi

      Do k = 1, model_levels
        If (div_damp_coefficient(k)  >   0.0 ) Then

! Calculate u derivative at rho points

          Do j = j0, j1
            Do i = 1, row_length
              divergence(i,j) = ( u(i,j,k) *                            &
     &                           (rho(i,j,k) + rho(i+1,j,k)) -          &
     &                           u(i-1,j,k) *                           &
     &                           (rho(i-1,j,k) + rho(i,j,k)) )*         &
     &                            recip_delta_lambda * 0.5*             &
     &                            FV_sec_theta_latitude(i,j)
            End Do
          End Do

! set values at northern and southern boundaries to zero.

          If (model_domain  /=  mt_bi_cyclic_lam) then
          If (at_extremity(PSouth)) Then
            Do i = 1, row_length
              divergence(i,1) = 0.
            End Do
          End If
          If (at_extremity(PNorth)) Then
            Do i = 1, row_length
              divergence(i,rows) = 0.
            End Do
          End If
          Endif

! Calculate v derivative at rho points

          Do j = j0, j1
            Do i = 1, row_length
              divergence(i,j) = divergence(i,j) +                       &
     &                         ( v(i,j,k) *                             &
     &                           (rho(i,j,k) + rho(i,j+1,k)) *          &
     &                           cos_v_latitude(i,j) -                  &
     &                           v(i,j-1,k) *                           &
     &                           (rho(i,j,k) + rho(i,j-1,k)) *          &
     &                           cos_v_latitude(i,j-1) ) * 0.5*         &
     &                          recip_delta_phi *                       &
     &                          FV_sec_theta_latitude(i,j)
            End Do
          End Do

          If (model_domain  ==  mt_Global) Then

! save values at poles for polar calculation
            If(at_extremity(PSouth))then
              Do i=1,row_length
                l_s_poles(i) = v(i,1,k) * .5 * ( rho(i,2,k)             &
     &                                          + rho(i,1,k))
              End Do
              Call gcg_rvecsumr(row_length, row_length, 1, 1,           &
     &                      l_s_poles, proc_row_group, info, sum_s)

              sum_s = sum_s * recip_delta_phi                           &
     &               * cos_v_latitude(1,1)                              &
     &               * FV_sec_theta_latitude(1,1) / g_row_length
              Do i = 1, row_length
                divergence(i,1) = sum_s
              End Do
            End If
            If(at_extremity(PNorth))then
              Do i=1,row_length
                l_n_poles(i) = v(i,n_rows,k) * .5 *                     &
     &                         ( rho(i,n_rows,k)+ rho(i,n_rows+1,k))
              End Do
              Call gcg_rvecsumr(row_length, row_length, 1, 1,           &
     &                         l_n_poles, proc_row_group, info, sum_n)
              sum_n = - sum_n * recip_delta_phi                         &
     &               * cos_v_latitude(1,n_rows)                         &
     &               * FV_sec_theta_latitude(1,rows) / g_row_length
              Do i = 1, row_length
                divergence(i,rows) = sum_n
              End Do
            End If
          End If

          Do j = 1, rows
            Do i = 1, row_length
              divergence(i,j) = divergence(i,j)                         &
     &                        / (rho(i,j,k) * r_rho_levels(i,j,k))
            End Do
          End Do

! swap divergence
! DEPENDS ON: swap_bounds
          call Swap_Bounds(                                             &
     &                   divergence, row_length, rows, 1,               &
     &                   off_x, off_y, fld_type_p, .false.)

! Calculate damping term applied to u
          Do j = j0, j1
            Do i = 1, row_length
              R_u(i,j,k) = R_u(i,j,k) +                                 &
     &                  div_damp_coefficient(k) * timestep *            &
     &                 (divergence(i+1,j) - divergence(i,j))            &
     &                 * recip_delta_lambda                             &
     &                 * FV_sec_theta_latitude(i,j)                     &
     &                 / r_at_u(i,j,k)
            End Do
          End Do

! Calculate damping term applied to v
          Do j = 1, n_rows
            Do i = 1, row_length
              R_v(i,j,k) = R_v(i,j,k) +                                 &
     &                  div_damp_coefficient(k) * timestep *            &
     &                 (divergence(i,j+1) - divergence(i,j))            &
     &                 * recip_delta_phi                                &
     &                 / r_at_v(i,j,k)
            End Do
          End Do

        End If ! on div damp coefficient
      End Do

! End of routine
      return
      END SUBROUTINE Div_damp

