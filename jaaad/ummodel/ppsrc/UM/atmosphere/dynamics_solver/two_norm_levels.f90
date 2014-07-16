
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Two_Norm_levels

      Subroutine Two_Norm_levels(                                       &
     &                          field, row_length, rows, levels,        &
     &                          start_level, end_level,                 &
     &                          model_domain, halo_i, halo_j,           &
     &                          at_extremity, n_procx, n_procy,         &
     &                          global_row_length, global_rows,         &
     &                          gc_proc_col_group, gc_proc_row_group,   &
     &                          l_datastart, L_do_halos,                &
     &                          L_do_rims, rims_to_do, Two_Norm )

! Purpose:
!          Calculates the two norm over a range of levels of the
!          input field.
!
! Method:
!          Is described in ;
!
!          Documentation yet to be written
!
! Original Progammer: Terry Davies
! Current code owner: Andrew J. Malcolm
!
! History:
! Version   Date     Comment
! ----     -------   -------
!  6.2    20/10/05   New routine based on GCR_Two_Norm
!                                                    Terry Davies
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
     &, L_do_rims                                                       &
                       ! LAMs: include rims, Global do all polar points
     &, L_do_halos     ! include halos in calculation of norm
!                    NB this means that some points are counted twice

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, levels                                                          &
                         ! number of levels in field.
     &, start_level                                                     &
                         ! start_level for norm calculation
     &, end_level        ! end_level for norm calculation

      Integer                                                           &
     &  model_domain                                                    &
                         ! holds integer code for model domain
     &, halo_i                                                          &
     &, halo_j                                                          &
     &, rims_to_do                                                      &
     &, l_datastart(3)                                                  &
     &, global_row_length                                               &
     &, global_rows                                                     &
     &, gc_proc_col_group                                               &
     &, gc_proc_row_group                                               &
     &, n_procx                                                         &
     &, n_procy

      Real                                                              &
     &  field(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j, levels)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

! Arguments with Intent OUT. ie: variables Output only

      Real                                                              &
     &  Two_Norm

! Local Variables.

      Integer                                                           &
     &  i, j, k                                                         &
     &, i_start                                                         &
     &, i_stop                                                          &
     &, j_start                                                         &
     &, j_stop                                                          &
     &, rims2                                                           &
     &, ipoints                                                         &
     &, istat

      Real                                                              &
     &  points           ! number of points used in norm

! Local arrays for parallel code

      real                                                              &
     &  two_norm_component( row_length+2*halo_i, rows+2*halo_j )        &
     &, two_norm_rows( rows+2*halo_j )

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


!  External Routines:
      External                                                          &
     &  gcg_rvecsumr

! ----------------------------------------------------------------------
! Section 1.   Calculate Norm.
! ----------------------------------------------------------------------

      rims2 = 2* rims_to_do
      i_start = 1
      i_stop = row_length
      j_start = 1
      j_stop = rows

      if ( L_do_rims ) then

        ipoints = global_row_length * global_rows

        if( L_do_halos ) then
          i_start = 1 - halo_i
          i_stop = row_length + halo_i
          j_start = 1 - halo_j
          j_stop = rows + halo_j
          ipoints = (global_row_length + 2 * n_procx * halo_i ) *       &
     &              (global_rows + 2 * n_procy * halo_j )
        endif ! L_do_halos

      else  ! do not do halos or rims

        if( L_do_halos ) then
          write(6,*)' You must do rims if halos are included in L2norms'
        endif ! L_do_halos

        If ( model_domain == mt_global ) Then
          if(at_extremity(PSouth)) j_start = 2
          if(at_extremity(PNorth)) j_stop = rows - 1
          ipoints = global_row_length * (global_rows - 2) + 2
        Endif ! model_domain == mt_global

        If( model_domain == mt_lam ) then
          if(at_extremity(PSouth)) j_start = 1 + rims_to_do
          if(at_extremity(PNorth)) j_stop = rows - rims_to_do
          if(at_extremity(PEast)) i_stop = row_length - rims_to_do
          if(at_extremity(PWest)) i_start = 1 + rims_to_do
          ipoints = (global_row_length - rims2) *                       &
     &              (global_rows  - rims2)

        Elseif (model_domain == mt_cyclic_LAM) then
          if(at_extremity(PSouth)) j_start = 1 + rims_to_do
          if(at_extremity(PNorth)) j_stop = rows - rims_to_do
          ipoints = global_row_length * (global_rows - rims2)

        End If ! model_domain == mt_lam

      endif ! L_do_rims

      points = float(ipoints)

      Do j = 1, rows + 2 * halo_j
        Do i = 1, row_length + 2 * halo_i
          two_norm_component(i,j) = 0.0
        End Do
      End Do

      If ( model_domain == mt_Global .and. .not. L_do_rims) Then
! Global model only calculate norm for one of the polar points
!  if L_do_rims is .false.

        If(at_extremity(PSouth) .and. (l_datastart(1) == 1)) then
          Do k = start_level, end_level
            two_norm_component(1,1)= two_norm_component(1,1) +          &
     &                               field(1,1,k) * field(1,1,k)
          End Do
        End If
        If(at_extremity(PNorth) .and. (l_datastart(1) == 1)) then
          Do k = start_level, end_level
            two_norm_component(1,rows)= two_norm_component(1,rows) +    &
     &                             field(1,rows,k) * field(1,rows,k)
          End Do
        End If
      End If  !  model_domain == mt_Global .and. .not. L_do_rims

      Do k = start_level, end_level
        Do j = j_start , j_stop
          Do i = i_start, i_stop
            two_norm_component( i + halo_i, j + halo_j ) =              &
     &          two_norm_component( i + halo_i, j + halo_j ) +          &
     &                                 field(i,j,k) * field(i,j,k)
          End Do
        End Do
      End Do  ! k = start_level, end_level

      call gcg_rvecsumr( row_length + 2*halo_i, row_length + 2*halo_i,  &
     &                   1, rows + 2*halo_j, two_norm_component,        &
     &                   gc_proc_row_group, istat, two_norm_rows )
      call gcg_rvecsumr( rows + 2*halo_j, rows + 2*halo_j, 1, 1,        &
     &                   two_norm_rows, gc_proc_col_group,              &
     &                   istat, two_norm )

! Normalize by dividing by number of points

      Two_Norm = SQRT (Two_Norm / points)

      return   ! End of routine
      END SUBROUTINE Two_Norm_levels

