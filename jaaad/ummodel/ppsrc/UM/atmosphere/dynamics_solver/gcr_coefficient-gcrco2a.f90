
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine GCR_Coefficient

      Subroutine GCR_Coefficient(                                       &
     &                           first_term, second_term,               &
     &                           row_length, rows, model_levels,        &
     &                           model_domain, coefficient,             &
     &                           at_extremity, n_proc,                  &
     &                           gc_proc_col_group, gc_proc_row_group,  &
     &                           l_datastart                            &
     &                           )

! Purpose:
!          Calculates a coefficient given by
!
!                              < first_term, second_term >
!          inner_product = -   __________________________
!
!                              < second_term, second_term >
!
! Method:
!          Is described in ;
!
!          Documentation yet to be written
!
! Original Progammer: Mark H. Mawson
! Current code owner: Andrew J. Malcolm
!
! History:
! Date     Version     Comment
! ----     -------     -------
! 22/03/00 5.1        Changed j_start/stop for cyclic LAM  Andy Malcolm
!LL   5.1   10/02/00  Use DOMTYP parameters                    P.Burton
!   5.3     15/09/01  add mt_bi_cyclic_LAM code            A. Malcolm
!     5.3   19/10/01  Use appropriate gcg routines.   S. Cusack
!     6.2   21/10/05  Remove unused external.              P.Selwood
!     6.2   19/10/05  Unroll loop for better vectorisation.
!                     K. Ketelsen/P.Selwood
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
     &, model_levels     ! number of model levels.

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

      Real                                                              &
     &  first_term(row_length,rows,model_levels)                        &
                                                  ! first term in inner
                                                  ! product definition
     &, second_term(row_length,rows,model_levels)
                                                  ! second term in inner
                                                  ! product definition

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

! Arguments with Intent OUT. ie: variables Output only

      Real                                                              &
     &  coefficient

! Local Variables.

      Integer                                                           &
     &  i, j, k                                                         &
     &, i_start                                                         &
     &, i_stop                                                          &
     &, j_start                                                         &
     &, j_stop

      Integer                                                           &
     &  l_datastart(3)

      Integer                                                           &
     &  istat                                                           &
     &, gc_proc_col_group                                               &
     &, gc_proc_row_group                                               &
     &, n_proc

! parallel Local arrays
      Real                                                              &
     &  inner_product(2)                                                &
     &, term(row_length,rows,2)                                         &
     &, sum_temp(rows,2)

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

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



! ----------------------------------------------------------------------
! Section 1.   Calculate Error Norm.
! ----------------------------------------------------------------------

! top term is in term (i,j,1)
! bottom term is in term (i,j,2)
      Do i= 1, row_length
        Do j= 1, rows
          term(i,j,1)=0.0
          term(i,j,2)=0.0
        End Do
      End Do

      If (model_domain  ==  mt_global .or.                              &
     &    model_domain  ==  mt_bi_cyclic_LAM) Then
! Solve over full domain
        i_start = 1
        i_stop = row_length
        j_start = 1
        j_stop = rows
      If (model_domain  ==  mt_global) Then
        if(at_extremity(PSouth))j_start=2
        if(at_extremity(PNorth))j_stop=rows-1
      Endif
      Else If(model_domain  ==  mt_lam)then
! Solve over interior points only.
        i_start = 1
        i_stop = row_length
        j_start = 1
        j_stop = rows
        if(at_extremity(PSouth)) j_start = 2
        if(at_extremity(PNorth)) j_stop = rows-1
        if(at_extremity(PEast)) i_stop = row_length-1
        if(at_extremity(PWest)) i_start = 2
      Else
! Solve over interior points only periodic in x => i_start=1,
! i_stop=row_length.
        i_start = 1
        i_stop = row_length
        j_start = 1
        j_stop = rows
        if(at_extremity(PSouth)) j_start = 2
        if(at_extremity(PNorth)) j_stop = rows-1
      End If

      If (model_domain  ==  mt_global) Then
! Global model only calculate inner product for one of the polar points.

        If(at_extremity(PSouth).and.(l_datastart(1) == 1))then
          Do k = 1, model_levels
            term(1,1,1)= term(1,1,1) +                                  &
     &                   first_term(1,1,k) * second_term(1,1,k)
            term(1,1,2)= term(1,1,2) +                                  &
     &                   second_term(1,1,k) * second_term(1,1,k)
          End Do
        End If
        if(at_extremity(PNorth).and.(l_datastart(1) == 1))then
          Do k = 1, model_levels
            term(1,rows,1)= term(1,rows,1) + first_term(1,rows,k)       &
     &                      * second_term(1,rows,k)
            term(1,rows,2)= term(1,rows,2) + second_term(1,rows,k)      &
     &                      * second_term(1,rows,k)
          End Do
        End If
      End If

      Do k = 1, model_levels
        Do j = j_start , j_stop
          Do i = i_start, i_stop
            term(i,j,1) = term(i,j,1) + first_term(i,j,k)               &
     &                                * second_term(i,j,k)
            term(i,j,2) = term(i,j,2) + second_term(i,j,k)              &
     &                                * second_term(i,j,k)
          End Do
        End Do
      End Do

      Call gcg_rvecsumr(row_length, row_length, 1, 2*rows, term,        &
     &                  gc_proc_row_group, istat, sum_temp)
      Call gcg_rvecsumr(rows, rows, 1, 2, sum_temp,                     &
     &                  gc_proc_col_group, istat, inner_product)


! Calculate coefficient

      coefficient = - inner_product(1) / inner_product(2)

! End of routine
      return
      END SUBROUTINE GCR_Coefficient

