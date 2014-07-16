
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine GCR_precon_ADI_setup

      Subroutine GCR_precon_ADI_setup(                                  &
     &                     HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,          &
     &                     HM_Czz, HM_Cz, HM_C3, HM_C4,                 &
     &                     eta_theta_levels, eta_rho_levels,            &
     &                     r_theta_levels, r_rho_levels,                &
     &                     row_length, rows, model_levels,              &
     &                     FV_sec_theta_latitude,                       &
     &                     model_domain, GCR_precon_option,             &
     &                     ADI_pseudo_timestep,                         &
     &                     rescale,                                     &
     &                     weight_upper, weight_lower,                  &
     &                     a0_z, a1_z, factor_z,                        &
     &                     a0_za, a1_za, factor_za,                     &
     &                     a0_y, a1_y, factor_y,                        &
     &                     a0_x, a1_x, factor_x,                        &
     &                     F_vector_x, G_vector_x,                      &
     &                     offx,offy,n_proc, global_row_length,         &
     &                     at_extremity, me,                            &
     &                     n_procx,n_procy,proc_row_group,              &
     &                     proc_col_group,halo_i,halo_j)

! Purpose:
!          Calculates ADI pre-conditioning matrix coefficients.
!          This code only for global model at the moment.
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
! Version  Date        Comment
! -------  ----        -------
!LL   5.1   11/02/00  Use DOMTYP parameters                    P.Burton
! 5.2      11/7/00     Modify calls to solver routines.  P.Burton
! 5.3      23/11/01    change magic_numbers.        A.Malcolm
! 5.3      19/10/01    Use appropriate gcg routines.   S. Cusack
!  6.0  18/08/03  NEC SX-6 optimisation - change dimension to avoid
!                 bank conflict.  R Barnes & J-C Rioual.
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit none

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, model_domain                                                    &
     &, GCR_precon_option  ! 0 = no preconditioning
                           ! 1 = Vertical block pre-conditioner
                           ! 2 = 1 iteration of 1 followed by 3D ADI
                           ! 3 = 3D ADI only
                           ! 4 = 1 iteration of 1 followed by xz ADI
                           ! 5 = 1 iteration of 1 followed by xz ADI

      Integer                                                           &
     &  offx,offy,n_proc, global_row_length,n_procx,n_procy             &
     &  ,halo_i,halo_j,me

      Integer                                                           &
     &  proc_row_group                                                  &
     &, proc_col_group

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid


! interpolation weights for moving between theta levels and rho levels
      Real                                                              &
     &  weight_upper(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)                                          &
     &, weight_lower(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)

      Real                                                              &
     &  ADI_pseudo_timestep

      Real                                                              &
     &  HM_Cxx1(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cxx2(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cyy1(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cyy2(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Czz(1-offx:row_length+offx,1-offy:rows+offy,model_levels)    &
     &, HM_Cz(1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, HM_C3(1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, HM_C4(1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, rescale (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, FV_sec_theta_latitude(1-offx:row_length+offx,1-offy:rows+offy)

      Real                                                              &
           ! vertical co-ordinate information
     &  r_theta_levels (1-halo_i:row_length+halo_i                      &
     &                 ,1-halo_j:rows+halo_j,0:model_levels)            &
     &, r_rho_levels (1-halo_i:row_length+halo_i                        &
     &                 ,1-halo_j:rows+halo_j,model_levels)              &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)

! Arguments with Intent OUT. ie: Output variables.
      Real                                                              &
     &  a0_za(row_length,rows,model_levels)                             &
     &, a1_za(row_length,rows,model_levels)                             &
     &, factor_za(row_length,rows,model_levels)                         &
     &, a0_z(row_length,rows,model_levels)                              &
     &, a1_z(row_length,rows,model_levels)                              &
     &, factor_z(row_length,rows,model_levels)                          &
     &, a0_x(row_length+1,rows,model_levels)                            &
     &, a1_x(row_length+1,rows,model_levels)                            &
     &, factor_x(row_length+1,rows,model_levels)                        &
     &, F_vector_x(row_length+1,rows,model_levels)                      &
     &, G_vector_x(2,rows,model_levels)                                 &
     &, a0_y(row_length,rows,model_levels)                              &
     &, a1_y(row_length,rows,model_levels)                              &
     &, factor_y(row_length,rows,model_levels)

! Local Variables.

      Integer                                                           &
     &  i,j,k, istat, ibase, len

      Real                                                              &
     &  factor_1                                                        &
     &, factor_2                                                        &
     &, factor_3                                                        &
     &, recip_ADI_pseudo_timestep

      Integer                                                           &
     &  i_start                                                         &
     &, i_end                                                           &
     &, j_start                                                         &
     &, j_end

! Local arrays.

      Real                                                              &
     &  sum_n(model_levels)                                             &
     &, sum_s(model_levels)                                             &
     &, sum_n_component(row_length,model_levels)                        &
     &, sum_s_component(row_length,model_levels)

      Real                                                              &
     &  a2_za(row_length,rows,model_levels)                             &
     &, a2_x(row_length,rows,model_levels)                              &
     &, a2_y(row_length,rows,model_levels)                              &
     &, f_vec_1(rows,model_levels)

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
! Description: COMDECK containing different preconditioner options
!
! Author : Andy Malcolm
! History:
! Version  Date      Comment.
! 5.3      23/10/01  New comdeck
!
      INTEGER                                                           &
     &     no_precon                                                    &
     &,    vert_precon                                                  &
     &,    vert_plus_xyz_ADI_precon                                     &
     &,    xyz_ADI_precon                                               &
     &,    vert_plus_xz_ADI_precon                                      &
     &,    xz_ADI_precon                                                &
     &,    Dufort_Frankel_precon

      PARAMETER(                                                        &
     &     no_precon               = 0                                  &
     &,    vert_precon             = 1                                  &
     &,    vert_plus_xyz_ADI_precon = 2                                 &
     &,    xyz_ADI_precon          = 3                                  &
     &,    vert_plus_xz_ADI_precon = 4                                  &
     &,    xz_ADI_precon           = 5                                  &
     &,    Dufort_Frankel_precon   = 6 )

!    External routines.
      External                                                          &
     &  GCG_rbcast, GCG_rvecsumr

!-----------------------------------------------------------------------
!     Section 1. Set initial guess to solution and initialise
!                variables.
!-----------------------------------------------------------------------

      recip_ADI_pseudo_timestep = 1./ ADI_pseudo_timestep

      j_start = 1
      j_end = rows
      If (at_extremity(PSouth) ) Then
        j_start = 2
      End If
      If (at_extremity(PNorth) ) Then
        j_end = rows-1
      End If

!-----------------------------------------------------------------------
!     Section 2. Set-up matrix coefficients in x direction.
!-----------------------------------------------------------------------

      If (model_domain  ==  mt_Global) Then
        i_start = 1
        i_end = row_length
        If (at_extremity(PWest) ) Then
          i_start = 2
        End If
        If (at_extremity(PEast) ) Then
          i_end = row_length-2
        End If

! Form matrix

        Do k = 1, model_levels

! zero f_vector
          Do j = 1, rows
            Do i = 1, row_length
              F_vector_x(i,j,k) = 0.
            End Do
          End Do

! Calculate coefficients at i=2, row_length-2
          Do j = j_start, j_end
            Do i = i_start, i_end
              a1_x(i,j,k) = HM_Cxx1(i,j,k)*HM_Cxx2(i,j,k) *             &
     &                      FV_sec_theta_latitude(i,j)
              a2_x(i,j,k) = HM_Cxx1(i-1,j,k)*HM_Cxx2(i-1,j,k) *         &
     &                      FV_sec_theta_latitude(i,j)
              a0_x(i,j,k) = - a2_x(i,j,k) - a1_x(i,j,k)                 &
     &                      - recip_ADI_pseudo_timestep                 &
     &                       * rescale(i,j,k)                           &
     &                      - HM_C4(i,j,k)
            End Do
          End Do

          If (at_extremity(PWest) ) Then
! Calculate coefficient at i=1
            i = 1
            Do j = j_start, j_end
              a1_x(i,j,k) = HM_Cxx1(i,j,k)*HM_Cxx2(i,j,k) *             &
     &                      FV_sec_theta_latitude(i,j)
              F_vector_x(i,j,k) = HM_Cxx1(i-1,j,k) *                    &
     &                          HM_Cxx2(i-1,j,k) *                      &
     &                      FV_sec_theta_latitude(i,j)
              a0_x(i,j,k) = - F_vector_x(i,j,k) - a1_x(i,j,k)           &
     &                      - recip_ADI_pseudo_timestep                 &
     &                       * rescale(i,j,k)                           &
     &                      - HM_C4(i,j,k)
            End Do
          End If

          If (at_extremity(PEast) ) Then
! Calculate coefficients at i= row_length-1
            i = row_length - 1
            Do j = j_start, j_end
              F_vector_x(i,j,k) = HM_Cxx1(i,j,k)*HM_Cxx2(i,j,k) *       &
     &                            FV_sec_theta_latitude(i,j)
              a2_x(i,j,k) = HM_Cxx1(i-1,j,k)*HM_Cxx2(i-1,j,k) *         &
     &                      FV_sec_theta_latitude(i,j)
              a0_x(i,j,k) = - a2_x(i,j,k) - F_vector_x(i,j,k)           &
     &                      - recip_ADI_pseudo_timestep                 &
     &                      * rescale(i,j,k)                            &
     &                      - HM_C4(i,j,k)
            End Do

! Calculate coefficients at i= row_length
            i = row_length
            Do j = j_start, j_end
              G_vector_x(1,j,k) = HM_Cxx1(i,j,k)*HM_Cxx2(i,j,k) *       &
     &                            FV_sec_theta_latitude(i,j)
              G_vector_x(2,j,k) = HM_Cxx1(i-1,j,k)*HM_Cxx2(i-1,j,k)*    &
     &                          FV_sec_theta_latitude(i,j)
              a0_x(i,j,k) = - G_vector_x(1,j,k) - G_vector_x(2,j,k)     &
     &                    - recip_ADI_pseudo_timestep                   &
     &                    * rescale(i,j,k)                              &
     &                      - HM_C4(i,j,k)
            End Do

          End If

        End Do ! end loop over levels

! solve cyclic tridiagonal system for solution
! reduce matrix to upper diagonal form.

        If (at_extremity(PWest) ) Then
          Do k = 1, model_levels
            Do j = j_start, j_end
              a0_x(1,j,k) = 1./a0_x(1,j,k)
            End Do
          End Do
        End If

!        Do i= 2, row_length - 1
!          Do k = 1, model_levels
!            Do j = j_start, j_end
!              factor_x(i,j,k) = a2_x(i,j,k) * a0_x(i-1,j,k)
!              a0_x(i,j,k) = 1./(a0_x(i,j,k) - factor_x(i,j,k)
!     &                          *a1_x(i-1,j,k))
!              F_vector_x(i,j,k) = F_vector_x(i,j,k) -
!     &                            factor_x(i,j,k)*F_vector_x(i-1,j,k)
!            End Do
!          End Do
!       End Do
! DEPENDS ON: solve_forward_x1
        Call solve_forward_x1(a0_x, a1_x, a2_x, factor_x, F_vector_x,   &
     &                      row_length, rows, model_levels,             &
     &                      j_start, j_end)

! back substitute to get solution

        If (at_extremity(PEast)) Then
          i = row_length
          Do k = 1, model_levels
            Do j = j_start, j_end
              F_vector_x(row_length-1,j,k) =                            &
     &                                   F_vector_x(row_length-1,j,k)   &
     &                                      * a0_x(row_length-1,j,k)
            End Do
          End Do
        End If

!        Do i= row_length - 2, 1, -1
!          Do k = 1, model_levels
!            Do j = j_start, j_end
!              F_vector_x(i,j,k) = a0_x(i,j,k) * (F_vector_x(i,j,k) -
!     &                                a1_x(i,j,k)*F_vector_x(i+1,j,k))
!            End Do
!          End Do
!        End Do

! DEPENDS ON: solve_backward_x1
        Call solve_backward_x1(F_vector_x, a0_x, a1_x,                  &
     &                       row_length, rows, model_levels,            &
     &                       j_start, j_end, 0, 0)

! Solve for last element of field
! First broadcast F_vector_1 to all row processors
        ibase = (me/n_procx) * n_procx
        If (at_extremity(PWest)) Then
          Do k = 1, model_levels
            Do j = j_start, j_end
              f_vec_1(j,k) = F_vector_x(1,j,k)
            End Do
          End Do
        End If

! Broadcast from western most (first on row ) processor on row to all
! others
        If (n_procx  >   1) Then
          len = model_levels * rows
          Call gcg_rbcast(201, len, ibase,                              &
     &                    proc_row_group, istat, f_vec_1)
        End If

        If (at_extremity(PEast)) Then
          i = row_length
          Do k = 1, model_levels
            Do j = j_start, j_end
              a0_x(i,j,k) = 1. / ( a0_x(i,j,k) - G_Vector_x(1,j,k)      &
     &                                        *F_vec_1(j,k) -           &
     &                                         G_vector_x(2,j,k) *      &
     &                                    F_vector_x(row_length-1,j,k))
            End Do
          End Do
        End If

!          Else
! Limited area not implemented so ignored.
      End If

!-----------------------------------------------------------------------
!     Section 3. Setup matrix coefficients in vertical direction.
!-----------------------------------------------------------------------

! Form matrix

      Do k = 1, model_levels
        If ( k  ==  1) Then
          factor_1 = 1. / (eta_theta_levels(k) -                        &
     &                     eta_theta_levels(k-1))
          factor_2 = 1. / (eta_rho_levels(k+1) -                        &
     &                       eta_rho_levels(k))
          Do j = 1, rows
            Do i = 1, row_length
              a1_za(i,j,k) = factor_2 * (factor_1 * HM_Czz(i,j,k) +     &
     &                      HM_C3(i,j,k) * HM_Cz(i,j,k) )
              a0_za(i,j,k) = - a1_za(i,j,k) - HM_C4(i,j,k)
            End Do
          End Do

        Else if (k  ==  model_levels) Then
          factor_1 = 1. / (eta_theta_levels(k) -                        &
     &                     eta_theta_levels(k-1))
          factor_3 = 1. / (eta_rho_levels(k) -                          &
     &                    eta_rho_levels(k-1))
          Do j = 1, rows
            Do i = 1, row_length
              a2_za(i,j,k) = factor_3 * (factor_1 * HM_Czz(i,j,k-1) -   &
     &                      HM_C3(i,j,k) * HM_Cz(i,j,k-1) *             &
     &                                weight_lower(i,j,k) )
              a0_za(i,j,k) = - a2_za(i,j,k) - HM_C4(i,j,k)
            End Do
          End Do
        Else
          factor_1 = 1. / (eta_theta_levels(k) -                        &
     &                     eta_theta_levels(k-1))
          factor_2 = 1. / (eta_rho_levels(k+1) -                        &
     &                     eta_rho_levels(k))
          factor_3 = 1. / (eta_rho_levels(k) -                          &
     &                     eta_rho_levels(k-1))

          Do j = 1, rows
            Do i = 1, row_length
              a1_za(i,j,k) = factor_2 * (factor_1 * HM_Czz(i,j,k) +     &
     &                       HM_C3(i,j,k) * HM_Cz(i,j,k) *              &
     &                                weight_upper(i,j,k) )
              a2_za(i,j,k) = factor_3 * (factor_1 * HM_Czz(i,j,k-1) -   &
     &                       HM_C3(i,j,k) * HM_Cz(i,j,k-1) *            &
     &                                weight_lower(i,j,k) )
              a0_za(i,j,k) = -a2_za(i,j,k)- a1_za(i,j,k)- HM_C4(i,j,k)
            End Do
          End Do
        End If
      End Do

        If (GCR_precon_option  ==  vert_plus_xyz_ADI_precon .or.        &
     &       GCR_precon_option  ==  vert_plus_xz_ADI_precon)            &
     &   Then
! need to set up coefficients for block vertical solve as well
! Copy a1_za to a1_z
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              a1_z(i,j,k) = a1_za(i,j,k)
            End Do
          End Do
        End Do

! Include horizontal second derivative terms in matrix.
! For model domain = 1 only

        Do k = 1, model_levels
          Do j = j_start, j_end
            Do i = 1, row_length
              a0_z(i,j,k) = a0_za(i,j,k)-FV_sec_theta_latitude(i,j) *   &
     &                  (HM_Cxx1(i-1,j,k)*HM_Cxx2(i-1,j,k) +            &
     &                   HM_Cxx1(i,j,k)*HM_Cxx2(i,j,k) +                &
     &                   HM_Cyy1(i,j,k)*HM_Cyy2(i,j,k) +                &
     &                   HM_Cyy1(i,j-1,k)*HM_Cyy2(i,j-1,k) )
            End Do
          End Do
        End Do

! polar boundaries
        If(at_extremity(PSouth))then
          Do k = 1, model_levels
            Do i = 1,row_length
              sum_s_component(i,k)= - HM_Cyy1(i,1,k)*HM_Cyy2(i,1,k)
            End Do
          End Do
        End If
        If(at_extremity(PNorth))then
          Do k = 1, model_levels
            Do i = 1,row_length
              sum_n_component(i,k)= - HM_Cyy1(i,rows-1,k)*              &
     &                                HM_Cyy2(i,rows-1,k)
            End Do
          End Do
        End If

        If(at_extremity(PSouth))then
          Call gcg_rvecsumr(row_length,row_length,1,                    &
     &                      model_levels,sum_s_component,               &
     &                      proc_row_group,istat,sum_s)
          Do k = 1, model_levels
            sum_s(k) = sum_s(k) * FV_sec_theta_latitude(1,1)            &
     &                   /global_row_length
            Do i = 1,row_length
              a0_z(i,1,k) = a0_za(i,1,k) + sum_s(k)
            End Do
          End Do
        End If

        If(at_extremity(PNorth))then
          Call gcg_rvecsumr(row_length,row_length,1,                    &
     &                      model_levels,sum_n_component,               &
     &                      proc_row_group,istat,sum_n)
          Do k = 1, model_levels
            sum_n(k) = sum_n(k) * FV_sec_theta_latitude(1,rows)         &
     &                   /global_row_length
            Do i = 1,row_length
              a0_z(i,rows,k) = a0_za(i,rows,k) + sum_n(k)
            End Do
          End Do
        End If

      End If ! on precon options

! For ADI coefficients add on relaxation pseudo_timestep
      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            a0_za(i,j,k) = a0_za(i,j,k)                                 &
     &                    - recip_ADI_pseudo_timestep                   &
     &                    * rescale(i,j,k)
          End Do
        End Do
      End Do

! Solve tridiagonal system.
! solution is in Soln
! reduce matrix to upper diagonal form.

      Do j = 1, rows
        Do i = 1, row_length
          a0_za(i,j,1) = 1./a0_za(i,j,1)
        End Do
      End Do
      Do k= 2, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            factor_za(i,j,k) = a2_za(i,j,k) * a0_za(i,j,k-1)
            a0_za(i,j,k) = 1./(a0_za(i,j,k) - factor_za(i,j,k)          &
     &                      * a1_za(i,j,k-1))
          End Do
        End Do
      End Do

        If (GCR_precon_option  ==  vert_plus_xyz_ADI_precon .or.        &
     &       GCR_precon_option  ==  vert_plus_xz_ADI_precon)            &
     &   Then
! Repeat for block vertical pre-conditioner
! Solve tridiagonal system.
! solution is in Soln
! reduce matrix to upper diagonal form.

        Do j = 1, rows
          Do i = 1, row_length
            a0_z(i,j,1) = 1./a0_z(i,j,1)
          End Do
        End Do
        Do k= 2, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              factor_z(i,j,k) = a2_za(i,j,k) * a0_z(i,j,k-1)
              a0_z(i,j,k) = 1./(a0_z(i,j,k) - factor_z(i,j,k)           &
     &                      * a1_z(i,j,k-1))
            End Do
          End Do
        End Do
      End If

!-----------------------------------------------------------------------
!     Section 4. Set-up matrix coefficients in y direction.
!-----------------------------------------------------------------------

        If (GCR_precon_option  ==  vert_plus_xyz_ADI_precon .or.        &
     &       GCR_precon_option  ==  xyz_ADI_precon) Then
! first do all non-polar points
        If (model_domain  ==  mt_Global) Then

          Do k = 1, model_levels
            Do j = j_start, j_end
              Do i = 1, row_length
                a1_y(i,j,k) = HM_Cyy1(i,j,k)*HM_Cyy2(i,j,k) *           &
     &                        FV_sec_theta_latitude(i,j)
                a2_y(i,j,k) = HM_Cyy1(i,j-1,k)*HM_Cyy2(i,j-1,k) *       &
     &                        FV_sec_theta_latitude(i,j)
                a0_y(i,j,k) = - a2_y(i,j,k) - a1_y(i,j,k)               &
     &                      - recip_ADI_pseudo_timestep                 &
     &                       * rescale(i,j,k)                           &
     &                      - HM_C4(i,j,k)
              End Do
            End Do

          End Do ! end loop over levels

! Solve tridiagonal system.
! solution is in Soln
! reduce matrix to upper diagonal form.

          If (at_extremity(PSouth)) Then
            Do k = 1, model_levels
              Do i = 1, row_length
                a0_y(i,2,k) = 1./a0_y(i,2,k)
              End Do
            End Do
          End If

!          Do k= 1, model_levels
!            Do j = 3, rows-1
!              Do i = 1, row_length
!                factor_y(i,j,k) = a2_y(i,j,k) * a0_y(i,j-1,k)
!                a0_y(i,j,k) = 1./(a0_y(i,j,k) - factor_y(i,j,k)
!     &                          *a1_y(i,j-1,k))
!              End Do
!            End Do
!          End Do

! DEPENDS ON: solve_forward_y1
           Call solve_forward_y1(a0_y, a1_y, a2_y, factor_y,            &
     &                         row_length, rows, model_levels)

        End If

      End If ! on precon option

!     end of routine GCR_precon_ADI_setup

      return
      END SUBROUTINE GCR_precon_ADI_setup

