
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine diffupper

      subroutine diffupper(                                             &
     &                     field, fld_type, off_v,                      &
     &                     levels, active_levels,                       &
     &                     in_rows, cos_rows, row_length,               &
     &                     off_x, off_y, halo_i, halo_j,                &
     &                     sec_latitude, cos_latitude,                  &
     &                     j_start, j_stop, pole_term,                  &
     &                     model_domain, at_extremity, proc_row_group,  &
     &                     start_level, up_diff, max_upd_levels )


! Purpose:
!          This routine used to evaluate diffusion at upper levels only.
!          Calculates conservative horizontal diffusion of a field
!          assuming that all model surfaces are at constant height.
!          This means we can cancel the delta_eta terms and the
!          r-terms can be absorbed into the diffusion coefficient.
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
! 6.2   24/12/05  This deck introduced based on tardiffq   Terry Davies
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Logical, Intent(In) ::                                            &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

      Integer, Intent(In) ::                                            &
     &  row_length                                                      &
                      ! number of point on a row.
     &, in_rows                                                         &
                      ! number of rows for field
     &, cos_rows                                                        &
                      ! number of rows for cos lat
     &, levels                                                          &
                      ! number of levels in field
     &, max_upd_levels                                                  &
                        ! max number of levels diffusion applied
     &, active_levels                                                   &
                      ! number of levels for diffusion
     &, start_level                                                     &
                      ! start level for diffusion
     &, j_start                                                         &
                      ! NS-Loop pointer
     &, j_stop                                                          &
                      ! NS-Loop pointer
     &, off_v                                                           &
                      ! Offset for cos pointer (v=1, u,theta=0)
     &, halo_i                                                          &
                      ! Size of halo in i direction.
     &, halo_j                                                          &
                      ! Size of halo in j direction.
     &, off_x                                                           &
                      ! Size of small halo in i
     &, off_y                                                           &
                      ! Size of small halo in j.
     &, model_domain                                                    &
                      ! holds integer code for model domain
     &, fld_type                                                        &
                      ! field type (v different from u, w, theta)
     &, proc_row_group

      Real, Intent(In) ::                                               &
     &  pole_term                                                       &
                     ! factor at pole calculated once
     &, up_diff(max_upd_levels)   ! effective diffusion coefficient

      Real, Intent(In) ::                                               &
     &  sec_latitude (1-off_x:row_length+off_x, 1-off_y:in_rows+off_y)  &
     &, cos_latitude (1-off_x:row_length+off_x, 1-off_y:cos_rows+off_y)


! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real, Intent(InOut) ::                                            &
     & field(1-halo_i:row_length+halo_i,1-halo_j:in_rows+halo_j,levels)

! Local Variables.

      Integer                                                           &
     &  i, j, k, ka                                                     &
                              ! Loop indices
     &, info

! Local arrays

      Real                                                              &
     &  lambda_term( row_length, in_rows, active_levels )               &
     &, phi_term( row_length, in_rows, active_levels )                  &
     &, temp(1-off_x:row_length+off_x,1-off_y:in_rows+off_y )           &
     &, l_s_poles(row_length, active_levels )                           &
     &, l_n_poles(row_length, active_levels )                           &
     &, sum_pole(active_levels )

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
! Section 1. Diffusion constant in wavenumber space
! ----------------------------------------------------------------------

      Do k = 1, active_levels
        ka = k + start_level - 1

! ----------------------------------------------------------------------
! Section 1.1  Calculate lambda direction term.
! ----------------------------------------------------------------------
        Do j = 1, in_rows
          Do i = 1-off_x, row_length
            temp(i,j) = ( field(i+1,j,ka) - field(i,j,ka) ) * up_diff(k)
          End Do
          Do i = 1, row_length
            lambda_term(i,j,k) = temp(i,j) - temp(i-1,j)
          End Do
        End Do

! ----------------------------------------------------------------------
! Section 1.2  Calculate phi direction term.
! ----------------------------------------------------------------------

        Do j = j_start-1, j_stop
          Do i = 1, row_length
            temp(i,j) = ( field(i,j+1,ka) - field(i,j,ka) ) *           &
     &                              cos_latitude(i,j+off_v) * up_diff(k)
          End Do
        End Do

        Do j = j_start, j_stop
          Do i = 1, row_length
            phi_term(i,j,k) = (temp(i,j) - temp(i,j-1)) *               &
     &                             sec_latitude(i,j)
          End Do
        End Do

! ----------------------------------------------------------------------
! Section 1.3  Polar term for theta, u points
! ----------------------------------------------------------------------
        If (model_domain == mt_Global .and. fld_type < fld_type_v ) Then

          If(at_extremity(PSouth))then
            Do i = 1, row_length
              l_s_poles(i,k) = (field(i,2,ka) - field(i,1,ka)) *        &
     &                          up_diff(k) * cos_latitude(i,1)
            End Do
          End If
          If(at_extremity(PNorth))then
            Do i = 1, row_length
              l_n_poles(i,k) = ( field(i, in_rows-1, ka) -              &
     &                           field(i, in_rows, ka) ) * up_diff(k) * &
     &                           cos_latitude(i,in_rows-1)
            End Do
          End If

        End If  ! model_domain == mt_Global .and. fld_type < fld_type_v

      EndDo ! k = 1, active_levels

      If (model_domain == mt_Global .and. fld_type < fld_type_v) Then

        If (at_extremity(PSouth)) Then
          Call gcg_rvecsumr(row_length, row_length, 1, active_levels,   &
     &                      l_s_poles, proc_row_group, info, sum_pole)
          Do k = 1, active_levels
            sum_pole(k) = pole_term * sum_pole(k)
            Do i = 1, row_length
              phi_term(i,1,k) = sum_pole(k)
            End Do
          End Do ! k = 1, active_levels
        End If  !  at_extremity(PSouth)

        If (at_extremity(PNorth)) Then
          Call gcg_rvecsumr(row_length, row_length, 1, active_levels,   &
     &                      l_n_poles, proc_row_group, info, sum_pole)
          Do k = 1, active_levels
            sum_pole(k) = pole_term * sum_pole(k)
            Do i = 1, row_length
              phi_term(i,in_rows,k) = sum_pole(k)
            End Do
          End Do ! k = 1, active_levels
        End If !  at_extremity(PNorth)

      End If !  model_domain == mt_Global .and. fld_type < fld_type_v

! ----------------------------------------------------------------------
! Section 2  Update fields
! ----------------------------------------------------------------------

      Do k = 1, active_levels
        ka = k + start_level - 1

        Do j = 1, in_rows
          Do i = 1, row_length
            field(i,j,ka) = field(i,j,ka) +                             &
     &                      lambda_term(i,j,k) + phi_term(i,j,k)
          End Do
        End Do

      EndDo ! k = 1, active_levels

      return    ! End subroutine diffupper
      END SUBROUTINE diffupper

