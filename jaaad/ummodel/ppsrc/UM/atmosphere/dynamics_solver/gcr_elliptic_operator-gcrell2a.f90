
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine GCR_Elliptic_Operator






      Subroutine GCR_Elliptic_Operator(                                 &
     &                                field, row_length, rows,          &
     &                                model_levels, model_domain,       &
     &                                first_constant_r_rho_level,       &
     &                                first_constant_r_rho_level_m1,    &
     &                                eta_theta_levels,                 &
     &                                eta_rho_levels,                   &
     &                                r_theta_levels, r_rho_levels,     &
     &                                FV_cos_theta_latitude,            &
     &                                HM_Cxx1, HM_Cxx2, HM_Cyy1,        &
     &                                HM_Cyy2, HM_Czz, HM_Cz,           &
     &                                HM_Cxz, HM_Cyz,  HM_Cxp, HM_Cyp,  &
     &                                HM_Cxy1, HM_Cxy2,                 &
     &                                HM_Cyx1, HM_Cyx2,                 &
     &                                HM_C2n,HM_C3, HM_C4, HM_C5,       &
     &                                weight_upper, weight_lower,       &
     &                                L_of_field,                       &
     &                                offx, offy, at_extremity,         &
     &                                n_rows, g_row_len, n_proc,        &
     &                                proc_row_group,                   &
     &                                halo_i,halo_j                     &
     &                                 )

! Purpose:
!          Applies Elliptic Operator L to input field.
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
! 22/03/00 5.1         Changed j bounds
!                      Changed logic for model_domain tests.
!                                                        Andy Malcolm
!LL   5.1   10/02/00  Use DOMTYP parameters                    P.Burton
!LL   5.2   27/09/00   bug fix (UMDP 15, page 8.6, aside)    A.Malcolm
!  5.3  17/07/01   add fujitsu compiler directives         Andy Malcolm
!  5.3  17/07/01   add fujitsu optimisations on switch     Andy Malcolm
!     5.3   19/10/01  Use appropriate gcg routines.   S. Cusack
!   5.3     15/09/01  add mt_bi_cyclic_LAM code            A. Malcolm
!  6.0  19/08/03  NEC SX-6 optimisation - loop unroll
!                 directives for NEC too.  R Barnes & J-C Rioual.
!  6.2  21/10/05  Remove commented out code P.Selwood
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
     &  offx,offy,n_rows,n_proc,halo_i,halo_j

      Integer                                                           &
     &  model_domain                                                    &
                         ! holds integer code for model domain
     &, first_constant_r_rho_level                                      &
                                   ! first rho level on which r
                                   ! is constant.
     &, first_constant_r_rho_level_m1 ! value used to dimension
                                      ! arrays, max of (1 and
                                      ! first_constant_r_rho_level)

      Real                                                              &
     &  field(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
                !input field

! interpolation weights for moving between theta levels and rho levels
      Real                                                              &
     &  weight_upper(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)                                          &
     &, weight_lower(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)

      Real                                                              &
            ! Coefficients of the elliptic operator
     &  HM_Cxx1 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, HM_Cxx2 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, HM_Cxy1 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, HM_Cxy2 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, HM_Cyy1 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, HM_Cyy2 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, HM_Cyx1 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, HM_Cyx2 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, HM_Czz (1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cz (1-offx:row_length+offx,1-offy:rows+offy,model_levels)    &
     &, HM_C2n (1-offx:row_length+offx,1-offy:rows+offy,                &
     &          first_constant_r_rho_level_m1)                          &
     &, HM_C3 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)    &
     &, HM_C4 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)    &
     &, HM_C5 (1-offx:row_length+offx,1-offy:rows+offy,                 &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cxz (1-offx:row_length+offx,1-offy:rows+offy,                &
     &          first_constant_r_rho_level_m1)                          &
     &, HM_Cyz (1-offx:row_length+offx,1-offy:rows+offy,                &
     &          first_constant_r_rho_level_m1)                          &
     &, HM_Cxp (1-offx:row_length+offx,1-offy:rows+offy,                &
     &          first_constant_r_rho_level_m1)                          &
     &, HM_Cyp (1-offx:row_length+offx,1-offy:rows+offy,                &
     &          first_constant_r_rho_level_m1)

      Real                                                              &
            ! Trigonometric functions
     &  FV_cos_theta_latitude (1-offx:row_length+offx,1-offy:rows+offy)
                              ! finite volume cosine array

      Real                                                              &
           ! vertical co-ordinate information
     &  r_theta_levels (1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,&
     &                  0:model_levels)                                 &
     &, r_rho_levels (1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,  &
     &                model_levels)                                     &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)

!   parallel variables Local
      Integer info

! parallel variables   IN
      Integer                                                           &
     &  proc_row_group                                                  &
     &, g_row_len

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid


! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

! Arguments with Intent OUT. ie: variables Output only

      Real                                                              &
     &  L_of_field(row_length,rows,model_levels)

! Local Variables.

      Integer                                                           &
     &  i, j, k                                                         &
                     ! Loop indices
     &, j_begin                                                         &
                     ! loop limits
     &, j_end                                                           &
                     ! loop limits
     &, i_start                                                         &
     &, i_stop                                                          &
     &, j_start                                                         &
     &, j_stop

      Real                                                              &
     &  term_1                                                          &
                   ! a term in the equation
     &, term_2                                                          &
                   ! another term in the equation
     &, factor_1                                                        &
                   ! another bit of the equation
     &, factor_2                                                        &
                   ! and another bit of the equation
     &, factor_3                                                        &
                   ! and yet another bit of the equation
     &, interp_upper                                                    &
                     ! weights for vertical averaging in eta
     &, interp_lower                                                    &
                     ! weights for vertical averaging in eta
     &, term_lower                                                      &
     &, recip_g_row_len

! Local arrays

      Real                                                              &
     &  first_deriv_x(1-offx:row_length+offx,1-offy:rows+offy)          &
     &, x_term(1-offx:row_length+offx,1-offy:rows+offy,model_levels)    &
     &, first_deriv_y(1-offx:row_length+offx,1-offy:n_rows+offy)        &
     &, y_term(1-offx:row_length+offx,1-offy:n_rows+offy,model_levels)  &
     &, interpx(1-offx:row_length+offx,1-offy:rows+offy)                &
     &, interpy(1-offx:row_length+offx,1-offy:n_rows+offy)              &
     &, term_upper(1-offx:row_length+offx,1-offy:rows+offy)             &
     &, l_s_poles(row_length,model_levels)                              &
     &, sum_s(model_levels)                                             &
     &, l_n_poles(row_length,model_levels)                              &
     &, sum_n(model_levels)                                             &
     &, factor_1_p(model_levels)

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
       External  swap_bounds

       integer i_start_no_halo,i_start_halo,i_end_halo,i_end_no_halo
       integer j_start_no_halo,j_start_halo,j_end_halo,j_end_no_halo
       integer j_start_all_rows,j_end_all_rows

       real tt,tp

! ----------------------------------------------------------------------
! Section 1.   Calculate Elliptic Operator applied to input field
! ----------------------------------------------------------------------

! Solve over full domain
        i_start = 1-offx
        i_stop = row_length+offx
        i_start_halo=1-offx
        i_start_no_halo=2-offx
        i_end_halo=row_length+offx
        i_end_no_halo=row_length+offx-1
        if(at_extremity(PSouth) .and.                                   &
     &     model_domain  /=  mt_bi_cyclic_LAM)then
          j_start_no_halo=2
          j_start_halo=2
          j_start_all_rows=1
          j_start = 1
        else
          j_start = 1-offy
          j_start_no_halo=2-offy
          j_start_halo=1-offy
          j_start_all_rows=1-offy
        endif
        if(at_extremity(PNorth) .and.                                   &
     &     model_domain  /=  mt_bi_cyclic_LAM)then
          j_stop = rows
          j_end_halo=rows-1
          j_end_no_halo=rows-1
          j_end_all_rows=rows
        else
          j_stop=rows+offy
          j_end_halo=rows+offy
          j_end_no_halo=rows+offy-1
          j_end_all_rows=rows+offy
        endif
! no lambda derivatives at poles.
        if(at_extremity(PSouth) .and.                                   &
     &     model_domain  /=  mt_bi_cyclic_LAM)then
          j_begin = 2
        else
          j_begin= 1-offy
        endif
        if(at_extremity(PNorth) .and.                                   &
     &     model_domain  /=  mt_bi_cyclic_LAM)then
          j_end = rows - 1
        else
          j_end = rows+offy
        endif
      If(model_domain  ==  mt_lam)then
! Lop off edge i points, j points are missing at edges automatically
! as the polar code is not executed
        If (at_extremity(PWest)) Then
          i_start = 1
          i_start_halo=1
          i_start_no_halo=2
        End If
        If (at_extremity(PEast)) Then
          i_stop = row_length
          i_end_halo=row_length
          i_end_no_halo=row_length-1
        End If
      Elseif (model_domain  ==  mt_cyclic_LAM) then
!  Solve over interior points only periodic in x => i_start=1,
! i_stop=row_length.
        i_start = 1
        i_stop = row_length
        if(at_extremity(PSouth))then
          j_start = 2
          j_begin = 2
        End If
        if(at_extremity(PNorth))then
          j_stop = rows - 1
          j_end = rows - 1
        End If
      Elseif(model_domain  ==  mt_bi_cyclic_LAM) then
        i_start = 1
        i_stop = row_length
        j_start = 1
        j_begin = 1
        j_stop = rows
        j_end = rows
      End If

      recip_g_row_len = 1. / g_row_len

      If (model_domain  ==  mt_lam) Then
! Initialise first derivatives to zero everywhere
        Do j = 1-offy,rows+offy
          Do i = 1-offx,row_length+offx
            first_deriv_x(i,j) = 0.0
          End Do
        End Do
        Do j = 1-offy,n_rows+offy
          Do i = 1-offx,row_length+offx
            first_deriv_y(i,j) = 0.0
          End Do
        End Do
      End If

! Loop over all levels

      Do k = 1, model_levels

! ----------------------------------------------------------------------
! Section 1.1   Calculate horizontal first derivatives.
! ----------------------------------------------------------------------

        If ( k  <   first_constant_r_rho_level ) Then

! first calculate first derivative with respect to x on constant r
! surface.

! calculate vertical derivative at all points.
! store in interpx to save re-calculation in y step.
          If (k  ==  1) Then
            do j=j_start_all_rows,j_end_all_rows
              do i=1-offx,row_length+offx
                term_upper(i,j) = HM_C2n(i,j,k) *                       &
     &                         (field(i,j,k+1) - field(i,j,k))
                interpx(i,j) = weight_upper(i,j,k) * term_upper(i,j)
              End Do
            End Do
          Else
!dir$ split
!dir$ unroll 4
            do j=j_start_all_rows,j_end_all_rows
!           Do i= i_start_halo,i_end_halo
             do i=1-offx,row_length+offx
                term_lower = term_upper(i,j)
                term_upper(i,j) = HM_C2n(i,j,k) *                       &
     &                         (field(i,j,k+1) - field(i,j,k))
                interpx(i,j) = weight_upper(i,j,k) * term_upper(i,j)    &
     &                       + weight_lower(i,j,k) * term_lower
              End Do
            End Do
          End If
! calculate first derivative on constant surface.
!!!   polar rows not included
          Do j = j_start_halo,j_end_halo
!dir$ unroll 8
            Do i = i_start_halo,i_end_no_halo
              first_deriv_x(i,j) = ( field(i+1,j,k) -                   &
     &                               field(i,j,k) )
            End Do
!dir$ unroll 4
            Do i = i_start_halo,i_end_no_halo
              first_deriv_x(i,j) = first_deriv_x(i,j)                   &
     &                               - HM_Cxp(i,j,k) * .5 *             &
     &                                ( interpx(i,j) +                  &
     &                                  interpx(i+1,j) )
            End Do
          End Do
! now calculate first derivative with respect to y on constant r
! surface.
!           do i=1-offx,row_length+offx
          Do j = j_start_all_rows,j_end_no_halo
!dir$ unroll 8
            Do i= i_start_halo,i_end_halo
              first_deriv_y(i,j) = ( field(i,j+1,k) -                   &
     &                               field(i,j,k) )
            End Do
!dir$ unroll 4
            Do i = i_start_halo,i_end_halo
              first_deriv_y(i,j) = first_deriv_y(i,j)                   &
     &                              - HM_Cyp(i,j,k) * .5 *              &
     &                                ( interpx(i,j) +                  &
     &                                  interpx(i,j+1) )
            End Do
          End Do

        Else
! surface is constant.
! first calculate first derivative with respect to x on constant r
! surface.
          Do j = j_start_halo,j_end_halo
!dir$ unroll 8
            Do i = i_start_halo,i_end_no_halo
              first_deriv_x(i,j) = field(i+1,j,k) -                     &
     &                             field(i,j,k)
            End Do
          End Do

! now calculate first derivative with respect to y on constant r
! surface.
          Do j = j_start_all_rows,j_end_no_halo
!dir$ unroll 8
            Do i= i_start_halo,i_end_halo
              first_deriv_y(i,j) = field(i,j+1,k) -                     &
     &                             field(i,j,k)
            End Do
          End Do

        End If

        If (model_domain  ==  mt_global) Then
          If(at_extremity(PSouth))then
            Do i=i_start_halo,i_end_no_halo
              first_deriv_x(i,1)=0.
            End Do
          End If
          If(at_extremity(PNorth))then
            Do i=i_start_halo,i_end_no_halo
              first_deriv_x(i,rows)=0.
            End Do
          End If
        Else If (model_domain  ==  mt_lam .or.                          &
     &           model_domain  ==  mt_cyclic_lam) Then
! zero derivatives on boundaries
          If(at_extremity(PSouth))then
            Do i=1-offx,row_length+offx
              first_deriv_x(i,1)=0.
              first_deriv_y(i,1)=0.
            End Do
          End If
          If(at_extremity(PNorth))then
            Do i=1-offx,row_length+offx
              first_deriv_x(i,rows)=0.
              first_deriv_y(i,n_rows)=0.
            End Do
          End If
        End if
        If (model_domain  ==  mt_lam) Then
          If(at_extremity(PWest))then
            Do j=1-offy,rows+offy
              first_deriv_x(1,j)=0.
            End Do
            Do j=1-offy,n_rows+offy
              first_deriv_y(1,j)=0.
            End Do
          End If
          If(at_extremity(PEast))then
            Do j=1-offy,rows+offy
              first_deriv_x(row_length,j)=0.
              first_deriv_x(row_length-1,j)=0.
            End Do
            Do j=1-offy,n_rows+offy
              first_deriv_y(row_length,j)=0.
            End Do
          End If
        End If

! ----------------------------------------------------------------------
! Section 1.2   Calculate the term to be differenced and averaged in x,
!               and the term to be differenced and averaged in y.
! ----------------------------------------------------------------------

! x_term.

        Do j=j_start_no_halo,j_end_no_halo

        i = i_start_halo-1

              tp   =        HM_Cxy2(i+1,j,k) * first_deriv_y(i+1,j) +   &
     &                      HM_Cxy2(i+1,j-1,k) * first_deriv_y(i+1,j-1)

!dir$ split
!dir$ unroll
          Do i=i_start_halo,i_end_no_halo

              tt   =        HM_Cxy2(i+1,j,k) * first_deriv_y(i+1,j) +   &
     &                      HM_Cxy2(i+1,j-1,k) * first_deriv_y(i+1,j-1)

            x_term(i,j,k) = HM_Cxx2(i,j,k) * first_deriv_x(i,j) +       &
     &                      .25 * HM_Cxy1(i,j,k) *                      &
     &                    ( tt + tp )
              tp=tt

          End Do
        End Do

! end rows
        If (model_domain  ==  mt_global) Then
          if(at_extremity(PSouth))then
            do i=i_start_halo,i_end_no_halo
              x_term(i,1,k)=0.
            end do
          endif
          if(at_extremity(PNorth))then
            do i=i_start_halo,i_end_no_halo
              x_term(i,rows,k)=0.
            end do
          endif
        End If

! y_term.

        Do j=j_start_all_rows,j_end_no_halo
        i=i_start_no_halo-1

                tp =        HM_Cyx2(i,j,k) * first_deriv_x(i,j) +       &
     &                      HM_Cyx2(i,j+1,k) * first_deriv_x(i,j+1)

!dir$ split
!dir$ unroll
          Do i=i_start_no_halo,i_end_no_halo

               tt =         HM_Cyx2(i,j,k) * first_deriv_x(i,j) +       &
     &                      HM_Cyx2(i,j+1,k) * first_deriv_x(i,j+1)

            y_term(i,j,k) = HM_Cyy2(i,j,k) * first_deriv_y(i,j) -       &
     &                      .25 * HM_Cyx1(i,j,k) *                      &
     &                    ( tt + tp )
               tp = tt

          End Do
        End Do

        If (model_domain  ==  mt_lam) Then
! Set zero derivatives on boundaries
          If (at_extremity(PWest)) Then
            do j=j_start_halo,j_end_halo
              x_term(1,j,k)=0.
            end do
          End If
          If (at_extremity(PEast)) Then
            do j=j_start_halo,j_end_halo
              x_term(row_length-1,j,k)=0.
              x_term(row_length,j,k)=0.
            end do
          End If
        End If
        If (model_domain  ==  mt_lam .or.                               &
     &      model_domain  ==  mt_cyclic_lam) Then
          If (at_extremity(PSouth)) Then
            do i=i_start_halo,i_end_halo
              y_term(i,1,k)=0.
            end do
          End If
          If (at_extremity(PNorth)) Then
            do i=i_start_halo,i_end_halo
              y_term(i,n_rows,k)=0.
            end do
          End If
        End If
! ----------------------------------------------------------------------
! Section 2.1   Calculate constant terms of operator applied to field.
! ----------------------------------------------------------------------

        Do j= 1, rows
!dir$ unroll 4
          Do i= 1, row_length
            L_of_field(i,j,k) = - HM_C4(i,j,k) * field(i,j,k)           &
     &                          * FV_cos_theta_latitude(i,j)
          End Do
        End Do

! ----------------------------------------------------------------------
! Section 2.2   Calculate lambda derivative terms of operator applied
!               to field, including d/dx(d/dy) terms.
! ----------------------------------------------------------------------

! Calculate Operator
        Do j=j_start_no_halo,j_end_no_halo
          Do i= i_start_no_halo, i_end_no_halo
            L_of_field(i,j,k) = L_of_field(i,j,k) +                     &
     &                       ( x_term(i,j,k) * HM_Cxx1(i,j,k) -         &
     &                         x_term(i-1,j,k) * HM_Cxx1(i-1,j,k))      &
!         End Do
!       End Do

! ----------------------------------------------------------------------
! Section 2.3   Calculate phi derivative terms of operator applied
!               to field, including d/dy(d/dx) terms..
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 2.3.1 Interior Points
! ----------------------------------------------------------------------

! Calculate Operator
!       Do j =j_start_no_halo,j_end_no_halo
!         Do i=i_start_no_halo,i_end_no_halo

!           L_of_field(i,j,k) = L_of_field(i,j,k) +
     &                                            +                     &
     &                    ( y_term(i,j,k) * HM_Cyy1(i,j,k) -            &
     &                      y_term(i,j-1,k) * HM_Cyy1(i,j-1,k) )
          End Do
        End Do

      End Do ! end loop over k levels

      If (model_domain  ==  mt_global) Then
! ----------------------------------------------------------------------
! Section 2.3.2 Poles in Global Model
! ----------------------------------------------------------------------

! South Pole
! Calculate Operator
        If(at_extremity(PSouth))then
          Do k=1, model_levels
            Do i=1,row_length
              l_s_poles(i,k) = y_term(i,1,k) * HM_Cyy1(i,1,k)
            End Do
          End Do
        End If

! North Pole
! Calculate Operator
        If(at_extremity(PNorth))then
          Do k=1, model_levels
            Do i=1,row_length
              l_n_poles(i,k) = y_term(i,rows-1,k)                       &
     &                         * HM_Cyy1(i,rows-1,k)
            End Do
          End Do
        End If


! average the value and add on constant term, note any other polar
! point will do as all values are the same.
        If (at_extremity(PSouth)) Then
          call gcg_rvecsumr(row_length,row_length,1,model_levels,       &
     &                      l_s_poles,proc_row_group,info,sum_s)
          Do k = 1, model_levels
            L_of_field(1,1,k) = L_of_field(1,1,k) +                     &
     &                          sum_s(k) * recip_g_row_len
! Copy answer at one point to all others
            Do i =2,row_length
              L_of_field(i,1,k) = L_of_field(1,1,k)
            End Do
          End Do
        End If
        If (at_extremity(PNorth)) Then
          call gcg_rvecsumr(row_length,row_length,1,model_levels,       &
     &                      l_n_poles,proc_row_group,info,sum_n)
          Do k = 1, model_levels
            L_of_field(1,rows,k) = L_of_field(1,rows,k)                 &
     &                             - sum_n(k) * recip_g_row_len
! Copy answer at one point to all others
            Do i =2,row_length
              L_of_field(i,rows,k) = L_of_field(1,rows,k)
            End Do
          End Do
        End If

      End If ! on model_domain type

      If (model_levels  >   1) Then
        Do k = 1, model_levels
! ----------------------------------------------------------------------
! Section 2.4   Calculate vertical derivative terms of operator applied
!               to field.
! ----------------------------------------------------------------------

          If ( k  ==  1) Then
            factor_1 = 1. / (eta_theta_levels(k) -                      &
     &                       eta_theta_levels(k-1))
            factor_2 = 1. / (eta_rho_levels(k+1) -                      &
     &                       eta_rho_levels(k))

            Do j= 1, rows
              Do i= 1, row_length
! Calculate upper derivative
                term_1 = ( field(i,j,k+1) - field(i,j,k) )              &
     &                    * factor_2

! Calculate lower derivative. This is zero by boundary condition.

! Calculate operator
                L_of_field(i,j,k) = L_of_field(i,j,k) +                 &
     &                              FV_cos_theta_latitude(i,j) *        &
     &                              ( factor_1 * term_1 *               &
     &                               HM_Czz(i,j,k) + HM_C3(i,j,k) *     &
     &                               term_1 * HM_Cz(i,j,k) )

              End Do
            End Do

          Else if (k  ==  model_levels) Then
            factor_1 = 1. / (eta_theta_levels(k) -                      &
     &                       eta_theta_levels(k-1))
            factor_3 = 1. / (eta_rho_levels(k) -                        &
     &                       eta_rho_levels(k-1))

            Do j= 1, rows
              Do i= 1, row_length

! Calculate upper derivative. This is zero by boundary condition.

! Calculate lower derivative.
                term_2 = ( field(i,j,k) - field(i,j,k-1) )              &
     &                    * factor_3

! Calculate operator
                L_of_field(i,j,k) = L_of_field(i,j,k) -                 &
     &                              FV_cos_theta_latitude(i,j) *        &
     &                              ( factor_1 *                        &
     &                                term_2 * HM_Czz(i,j,k-1)          &
     &                               - HM_C3(i,j,k) *                   &
     &                                term_2 * HM_Cz(i,j,k-1) *         &
     &                                weight_lower(i,j,k) )

              End Do
            End Do

          Else
            factor_1 = 1. / (eta_theta_levels(k) -                      &
     &                       eta_theta_levels(k-1))
            factor_2 = 1. / (eta_rho_levels(k+1) -                      &
     &                       eta_rho_levels(k))
            factor_3 = 1. / (eta_rho_levels(k) -                        &
     &                       eta_rho_levels(k-1))

            Do j= 1, rows
!dir$ split
!dir$ unroll 4
              Do i= 1, row_length

                L_of_field(i,j,k) = L_of_field(i,j,k) +                 &
     &                              FV_cos_theta_latitude(i,j) *        &
     &                    ( field(i,j,k+1) - field(i,j,k) )             &
     &                    * factor_2 * (                                &
     &                              factor_1 *                          &
     &                                 HM_Czz(i,j,k)                    &
     &                               + HM_C3(i,j,k) *                   &
     &                                 HM_Cz(i,j,k) *                   &
     &                                weight_upper(i,j,k) )
              Enddo
            Enddo
            Do j= 1, rows
!dir$ split
!dir$ unroll 4
              Do i= 1, row_length

                L_of_field(i,j,k) = L_of_field(i,j,k) -                 &
     &                              FV_cos_theta_latitude(i,j) *        &
     &                   ( field(i,j,k) - field(i,j,k-1) )              &
     &                    * factor_3 * (                                &
     &                              factor_1 *                          &
     &                                 HM_Czz(i,j,k-1)                  &
     &                               - HM_C3(i,j,k) *                   &
     &                                 HM_Cz(i,j,k-1) *                 &
     &                                weight_lower(i,j,k) )

              End Do
            End Do

          End If

! End loop over model levels
        End Do

! End if on there being more than 1 level
      End If

! ----------------------------------------------------------------------
! Section 2.5   Calculate mixed vertical derivative terms of operator
!               applied to field.
! ----------------------------------------------------------------------

      If (first_constant_r_rho_level  >   1 ) Then

        Do k = 1, first_constant_r_rho_level

          If (k  ==  1 ) Then

! Interpolate quantities to levels either side in vertical
           interp_upper = (eta_theta_levels(k) - eta_rho_levels(k)) /   &
     &                     (eta_rho_levels(k+1) - eta_rho_levels(k))
            interp_lower = 1.0 - interp_upper

            Do j=j_start_no_halo,j_end_no_halo
!dir$ unroll 4
              Do i=i_start_halo,i_end_no_halo
                interpx(i,j) = ( interp_upper * x_term(i,j,k+1)         &
     &                             + interp_lower * x_term(i,j,k) ) *   &
     &                           HM_Cxz(i,j,k)
              End Do
            End Do

            Do j=j_start_all_rows,j_end_no_halo
!dir$ unroll 4
              Do i=i_start_no_halo,i_end_no_halo
                interpy(i,j) = (interp_upper * y_term(i,j,k+1)          &
     &                            + interp_lower * y_term(i,j,k) ) *    &
     &                           HM_Cyz(i,j,k)
              End Do
            End Do
! Calculate 1. / delta eta.
            factor_1_p(k) = 1. / (eta_theta_levels(k) -                 &
     &                            eta_theta_levels(k-1))

            Do j= j_start_no_halo, j_end_no_halo
              Do i= i_start_no_halo,i_end_no_halo
                term_upper(i,j) = 0.5 * ( interpx(i,j)                  &
     &                                    + interpx(i-1,j)              &
     &                                    + interpy(i,j)                &
     &                                    + interpy(i,j-1) )            &
     &                          * HM_C5(i,j,k)

                L_of_field(i,j,k) = L_of_field(i,j,k) -                 &
     &                              term_upper(i,j)                     &
     &                              * factor_1_p(k)                     &
     &                              * FV_cos_theta_latitude(i,j)
              End Do
            End Do

! save terms to perform calculation at boundaries
! Boundaries assumed flat in limited area model so code omitted.
            If (model_domain  ==  mt_global) Then
! south pole
              If(at_extremity(PSouth))then
                Do i= 1, row_length
                  l_s_poles(i,k) = interpy(i,1)
                End Do
              End If
! north pole
              if(at_extremity(PNorth))then
                do i= 1, row_length
                  l_n_poles(i,k) = interpy(i,n_rows)
                End Do
              End If

            End If


          Else If (k  ==  first_constant_r_rho_level ) Then

! Calculate 1. / delta eta.
            factor_1_p(k) = 1. / (eta_theta_levels(k) -                 &
     &                            eta_theta_levels(k-1))

            Do j=j_start_no_halo,j_end_no_halo
              Do i=i_start_no_halo,i_end_no_halo

                L_of_field(i,j,k) = L_of_field(i,j,k) +                 &
     &                              term_upper(i,j)                     &
     &                           * factor_1_p(k)                        &
     &                           * FV_cos_theta_latitude(i,j)
              End Do
            End Do

          Else

! Interpolate quantities to levels either side in vertical
            interp_upper = (eta_theta_levels(k) - eta_rho_levels(k)) /  &
     &                     (eta_rho_levels(k+1) - eta_rho_levels(k))
            interp_lower = 1.0 - interp_upper

            Do j=j_start_no_halo,j_end_no_halo
              Do i=i_start_halo,i_end_no_halo

                interpx(i,j) = ( interp_upper * x_term(i,j,k+1)         &
     &                           + interp_lower * x_term(i,j,k) )       &
     &                           * HM_Cxz(i,j,k)
              End Do
            End Do

            Do j=j_start_all_rows,j_end_no_halo
              Do i=i_start_no_halo,i_end_no_halo

                interpy(i,j) = ( interp_upper * y_term(i,j,k+1)         &
     &                           + interp_lower * y_term(i,j,k) )       &
     &                           * HM_Cyz(i,j,k)
              End Do
            End Do

! Calculate 1. / delta eta.
            factor_1_p(k) = 1. / (eta_theta_levels(k) -                 &
     &                            eta_theta_levels(k-1))

            Do j=j_start_no_halo,j_end_no_halo
              Do i=i_start_no_halo,i_end_no_halo

                term_lower = term_upper(i,j)
                term_upper(i,j) = 0.5 * ( interpx(i,j)                  &
     &                                    + interpx(i-1,j)              &
     &                                    + interpy(i,j)                &
     &                                    + interpy(i,j-1) )            &
     &                          * HM_C5(i,j,k)

                L_of_field(i,j,k) = L_of_field(i,j,k) -                 &
     &                           (  term_upper(i,j)                     &
     &                            - term_lower )                        &
     &                           * factor_1_p(k)                        &
     &                           *  FV_cos_theta_latitude(i,j)
              End Do
            End Do

! save terms to perform calculation at boundaries
! Boundaries assumed flat in limited area model so code omitted.
            If (model_domain  ==  mt_global) Then
! south pole
              If(at_extremity(PSouth))then
                Do i= 1, row_length
                  l_s_poles(i,k) = interpy(i,1)
                End Do
              End If
! north pole
              if(at_extremity(PNorth))then
                do i= 1, row_length
                  l_n_poles(i,k) = interpy(i,n_rows)
                End Do
              End If

            End If

          End If

        End Do

! Now do calculations at boundary.

! average the value and add on constant term, note any other polar
! point will do as all values are the same.
        If (at_extremity(PSouth) .and.                                  &
     &      (model_domain  ==  mt_global)) Then
          call gcg_rvecsumr(row_length,row_length,1,                    &
     &                      first_constant_r_rho_level-1,               &
     &                      l_s_poles,proc_row_group,info,sum_s)

          Do k = 1, first_constant_r_rho_level

            If (k  ==  1) Then
             term_upper(1,1) = 2.0 * sum_s(k) * HM_C5(1,1,k) *          &
     &                          recip_g_row_len
              L_of_field(1,1,k) = L_of_field(1,1,k) - term_upper(1,1)   &
     &                            * factor_1_p(k)                       &
     &                            * FV_cos_theta_latitude(1,1)
              Do i = 2, row_length
                L_of_field(i,1,k) = L_of_field(1,1,k)
              End Do

            Else If (k  ==  first_constant_r_rho_level) Then
              L_of_field(1,1,k) = L_of_field(1,1,k) + term_upper(1,1)   &
     &                            * factor_1_p(k)                       &
     &                            * FV_cos_theta_latitude(1,1)
              Do i = 2, row_length
                L_of_field(i,1,k) = L_of_field(1,1,k)
              End Do

            Else
              term_lower = term_upper(1,1)
              term_upper(1,1) = 2.0 * sum_s(k) * HM_C5(1,1,k) *         &
     &                          recip_g_row_len
              L_of_field(1,1,k) = L_of_field(1,1,k) -                   &
     &                           ( term_upper(1,1) -                    &
     &                             term_lower )                         &
     &                              * factor_1_p(k)                     &
     &                              * FV_cos_theta_latitude(1,1)
              Do i = 2, row_length
                L_of_field(i,1,k) = L_of_field(1,1,k)
              End Do

            End If
          End Do
        End If

        If (at_extremity(PNorth) .and.                                  &
     &      (model_domain  ==  mt_global)) Then
          call gcg_rvecsumr(row_length,row_length,1,                    &
     &                      first_constant_r_rho_level-1,               &
     &                      l_n_poles,proc_row_group,info,sum_n)

          j = rows
          Do k = 1, first_constant_r_rho_level

            If (k  ==  1) Then
              term_upper(1,j) = 2.0 * sum_n(k) * HM_C5(1,j,k) *         &
     &                          recip_g_row_len
              L_of_field(1,j,k) = L_of_field(1,j,k) - term_upper(1,j)   &
     &                              * factor_1_p(k)                     &
     &                              * FV_cos_theta_latitude(1,j)
              Do i = 2, row_length
                L_of_field(i,j,k) = L_of_field(1,j,k)
              End Do

            Else If (k  ==  first_constant_r_rho_level) Then
              L_of_field(1,j,k) = L_of_field(1,j,k) + term_upper(1,j)   &
     &                              * factor_1_p(k)                     &
     &                              * FV_cos_theta_latitude(1,j)
              Do i = 2, row_length
                L_of_field(i,j,k) = L_of_field(1,j,k)
              End Do

            Else
              term_lower = term_upper(1,j)
              term_upper(1,j) = 2.0 * sum_n(k) * HM_C5(1,j,k) *         &
     &                          recip_g_row_len
              L_of_field(1,j,k) = L_of_field(1,j,k) -                   &
     &                              ( term_upper(1,j) -                 &
     &                                term_lower )                      &
     &                              * factor_1_p(k)                     &
     &                              * FV_cos_theta_latitude(1,j)
              Do i = 2, row_length
                L_of_field(i,j,k) = L_of_field(1,j,k)
              End Do

            End If
          End Do
        End If

      End If

! End of routine  GCR_Elliptic_Operator
      return
      END SUBROUTINE GCR_Elliptic_Operator

