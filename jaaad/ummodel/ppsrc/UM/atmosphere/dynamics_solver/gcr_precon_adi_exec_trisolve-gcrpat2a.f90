
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine GCR_precon_ADI_exec_trisolve

      Subroutine GCR_precon_ADI_exec_trisolve(j_start,j_end,            &
     &                     r, HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,       &
     &                     HM_Czz, HM_Cz,                               &
     &                     HM_Cxz, HM_Cyz, HM_Cxp, HM_Cyp,              &
     &                     HM_Cxy1, HM_Cxy2, HM_Cyx1, HM_Cyx2,          &
     &                     HM_C2n, HM_C3, HM_C4, HM_C5,                 &
     &                     row_length, rows, n_rows, model_levels,      &
     &                     first_constant_r_rho_level,                  &
     &                     first_constant_r_rho_level_m1,               &
     &                     eta_theta_levels, eta_rho_levels,            &
     &                     r_theta_levels, r_rho_levels,                &
     &                     FV_sec_theta_latitude,                       &
     &                     FV_cos_theta_latitude,                       &
     &                     model_domain, GCR_precon_option,             &
     &                     GCR_adi_add_full_soln,                       &
     &                     ADI_pseudo_timestep,                         &
     &                     GCR_n_ADI_pseudo_timesteps, rescale,         &
     &                     weight_upper, weight_lower,                  &
     &                     a0_z, a1_z, factor_z,                        &
     &                     a0_y, a1_y, factor_y,                        &

! new
     &                     a_plus_x, a_minus_x,                         &
     &                     factor_forward, factor_backward,             &
     &                     recip_a_central_x,                           &
     &                     bv_a_matrix_0, bv_a_matrix_np1,              &
     &                     recip_bv_a_matrix_diag, bv_a_matrix_sup,     &
     &                     bv_factor_forward, bv_factor_backward,       &
     &                     bv_soln_n_term,                              &
     &                     bv_soln_1_term1, bv_soln_1_term2,            &

     &                     offx,offy,n_proc, global_row_length,         &
     &                     at_extremity, neighbour, me,                 &
     &                     n_procx,n_procy,proc_row_group,              &
     &                     proc_col_group,halo_i,halo_j,                &
     &                     Soln)

! Purpose:
!          Calculates ADI pre-conditioning operator applied to field.
!
! Method:
!          Is described in ;
!
!          Documentation yet to be written
!
!
! Original Progammer: Mark H. Mawson
! Current code owner: Andrew J. Malcolm
!
! History:
! Version  Date        Comment
! -------  ----        -------
! 5.2      11/7/00     Modify calls to solver routines.  P.Burton
! 5.3      24/10/01    change magic_numbers. M.Diamantakis
! 5.3     19/10/01     Use appropriate gcg routines.   S. Cusack
!  6.0  18/08/03  NEC SX-6 optimisation - change dimension to avoid
!                 bank conflict and T3E optimised code replaced by
!                 straightforward loops.  R Barnes & J-C Rioual.
!  6.2  03/11/04  change dimensions to allow loop collapsing in
!                 Mpp_tri_solve_exec.
!                 Unroll loops by directive.      Klaus Ketelsen/
!                 Add jstrt and j_end to parameter list. R Barnes
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit none

! Arguments with Intent IN. ie: Input variables.

      integer, intent(in)                   :: j_start,j_end

      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, n_rows                                                          &
     &, model_levels                                                    &
     &, first_constant_r_rho_level                                      &
                                   ! first rho level on which r
                                   ! is constant.
     &, first_constant_r_rho_level_m1                                   &
                                      ! value used to dimension
                                      ! arrays, max of (1 and
                                      ! first_constant_r_rho_level)
     &, model_domain                                                    &
     &, GCR_n_ADI_pseudo_timesteps                                      &
     &, GCR_precon_option  ! 0 = no preconditioning
                           ! 1 = Vertical block pre-conditioner
                           ! 2 = 1 iteration of 1 followed by 3D ADI
                           ! 3 = 3D ADI only
                           ! 4 = 1 iteration of 1 followed by xz ADI
                                ! 5 = xz ADI only

      Logical                                                           &
     &  GCR_adi_add_full_soln ! true then use full equation on RHS
                            ! on second and subsequent ADI timesteps

      Integer                                                           &
     &  offx,offy,n_proc, global_row_length,n_procx,n_procy             &
     &  ,halo_i,halo_j, me                                              &
     &, neighbour(4)

      Integer                                                           &
     &  proc_row_group                                                  &
     &, proc_col_group

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Parameters
      Integer                                                           &
     &   PNorth,                                                        &
                      ! North processor address in the neighbor array
     &   PEast,                                                         &
                      ! East processor address in the neighbor array
     &   PSouth,                                                        &
                      ! South processor address in the neighbor array
     &   PWest,                                                         &
                      ! West processor address in the neighbor array
     &   NoDomain     ! Value in neighbor array if the domain has
                      !  no neighbor in this direction. Otherwise
                      !  the value will be the tid of the neighbor
      Parameter (                                                       &
     &   PNorth   = 1,                                                  &
     &   PEast    = 2,                                                  &
     &   PSouth   = 3,                                                  &
     &   PWest    = 4,                                                  &
     &   NoDomain = -1)

      Real                                                              &
     &  ADI_pseudo_timestep

! interpolation weights for moving between theta levels and rho levels
      Real                                                              &
     &  weight_upper(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)                                          &
     &, weight_lower(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)

      Real                                                              &
           ! vertical co-ordinate information
     &  r_theta_levels (1-halo_i:row_length+halo_i                      &
     &                 ,1-halo_j:rows+halo_j,0:model_levels)            &
     &, r_rho_levels (1-halo_i:row_length+halo_i                        &
     &                 ,1-halo_j:rows+halo_j,model_levels)              &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)

      Real                                                              &
     &  HM_Cxx1(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cxx2(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cxy1 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, HM_Cxy2 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, HM_Cyy1(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cyy2(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cyx1 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, HM_Cyx2 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, HM_Czz(1-offx:row_length+offx,1-offy:rows+offy,model_levels)    &
     &, HM_Cz(1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, HM_C2n (1-offx:row_length+offx,1-offy:rows+offy,                &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_C3(1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, HM_C4(1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, HM_C5 (1-offx:row_length+offx,1-offy:rows+offy,                 &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cxz (1-offx:row_length+offx,1-offy:rows+offy,                &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cyz (1-offx:row_length+offx,1-offy:rows+offy,                &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cxp (1-offx:row_length+offx,1-offy:rows+offy,                &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cyp (1-offx:row_length+offx,1-offy:rows+offy,                &
     &         first_constant_r_rho_level_m1)                           &
     &, r(row_length,rows,model_levels)                                 &
     &, rescale (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, FV_sec_theta_latitude(1-offx:row_length+offx,1-offy:rows+offy)  &
     &, FV_cos_theta_latitude(1-offx:row_length+offx,1-offy:rows+offy)

      Real                                                              &
     &  a_plus_x(row_length+1,rows,model_levels)                        &
     &, a_minus_x(row_length+1,rows,model_levels)

      Real                                                              &
     &  factor_forward(row_length+1,j_start:j_end,model_levels)         &
     &, factor_backward(row_length+1,j_start:j_end,model_levels)        &
     &, bv_a_matrix_0(2*n_procx,rows,model_levels)                      &
     &, bv_a_matrix_np1(2*n_procx,rows,model_levels)                    &
     &, bv_factor_forward(2,2*n_procx,rows,model_levels)                &
     &, bv_factor_backward(2*n_procx,rows,model_levels)                 &
     &, recip_bv_a_matrix_diag(2*n_procx,rows,model_levels)             &
     &, bv_a_matrix_sup(2*n_procx,rows,model_levels)                    &
     &, recip_a_central_x(row_length,rows,model_levels)

      Real                                                              &
     &  bv_soln_n_term(rows, model_levels)                              &
     &, bv_soln_1_term1(rows, model_levels)                             &
     &, bv_soln_1_term2(rows, model_levels)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.
      Real                                                              &
     &  Soln(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

! Local Variables.

      Real                                                              &
     &  term_1                                                          &
     &, term_2                                                          &
     &, factor_1                                                        &
     &, factor_2                                                        &
     &, factor_3

      Integer                                                           &
     &  i,j,k                                                           &
     &, pseudo_timestep_number

      Integer                 :: istat

      Real                                                              &
     &  recip_ADI_pseudo_timestep

! Local arrays.

!kk   1. dimension row_length+1 to avoid copy in Mpp_tri_solve_exec
      real,dimension(row_length+1,rows,model_levels)           :: RHS
!kk   To save copying, a secon array with j_start:j_end is required
      real,dimension(row_length+1,j_start:j_end,model_levels)  :: RHS_c

      Real                                                              &
     &  a0_z(row_length,rows,model_levels)                              &
     &, a1_z(row_length,rows,model_levels)                              &
     &, factor_z(row_length,rows,model_levels)                          &
     &, a0_y(row_length,rows,model_levels)                              &
     &, a1_y(row_length,rows,model_levels)                              &
     &, factor_y(row_length,rows,model_levels)                          &
     &, soln_prev(row_length,rows,model_levels)                         &
     &, L_of_soln(row_length,rows,model_levels)

      Real                                                              &
     &  sum_n(model_levels)                                             &
     &, sum_s(model_levels)                                             &
     &, sum_n_component(row_length,model_levels)                        &
     &, sum_s_component(row_length,model_levels)

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
     &  GCG_rvecsumr, GCG_rbcast,mpp_tri_solve_exec

!-----------------------------------------------------------------------
!     Section 1. Set initial guess to solution and initialise
!                variables.
!-----------------------------------------------------------------------


      recip_ADI_pseudo_timestep = 1./ ADI_pseudo_timestep

      Do pseudo_timestep_number = 1, GCR_n_ADI_pseudo_timesteps

! Form matrix
        If (pseudo_timestep_number  ==  1) Then
          Do k = 1, model_levels
            Do j = j_start,j_end
              Do i = 1, row_length
                RHS_c(i,j,k) = r(i,j,k) * FV_sec_theta_latitude(i,j)
              End Do
            End Do
!           If needed, set boundary rows directly in RHS
            Do j = 1,j_start-1
              Do i = 1, row_length
                RHS(i,j,k) = r(i,j,k) * FV_sec_theta_latitude(i,j)
              End Do
            End Do
            Do j = j_end+1, rows
              Do i = 1, row_length
                RHS(i,j,k) = r(i,j,k) * FV_sec_theta_latitude(i,j)
              End Do
            End Do
          End Do
        Else
! swop soln on all procesors
! DEPENDS ON: swap_bounds
        Call swap_bounds(soln,row_length,rows,model_levels,             &
     &                   offx,offy,fld_type_p,.FALSE.)

        If (GCR_adi_add_full_soln) Then
! DEPENDS ON: gcr_elliptic_operator
          Call GCR_Elliptic_Operator(                                   &
     &                           soln, row_length,                      &
     &                           rows, model_levels, model_domain,      &
     &                           first_constant_r_rho_level,            &
     &                           first_constant_r_rho_level_m1,         &
     &                           eta_theta_levels, eta_rho_levels,      &
     &                           r_theta_levels, r_rho_levels,          &
     &                           FV_cos_theta_latitude,                 &
     &                           HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,    &
     &                           HM_Czz, HM_Cz,                         &
     &                           HM_Cxz, HM_Cyz, HM_Cxp, HM_Cyp,        &
     &                           HM_Cxy1, HM_Cxy2, HM_Cyx1, HM_Cyx2,    &
     &                           HM_C2n, HM_C3, HM_C4, HM_C5,           &
     &                           weight_upper, weight_lower,            &
     &                           L_of_soln,                             &
     &                           offx, offy, at_extremity, n_rows,      &
     &                           global_row_length, n_proc,             &
     &                           proc_row_group, halo_i,halo_j          &
     &                           )

! add on to RHS and save solution from previous timestep.
            Do k = 1, model_levels
              Do j = 1, rows
                Do i = 1, row_length
                RHS(i,j,k) = r(i,j,k) * FV_sec_theta_latitude(i,j)      &
     &                       - L_of_soln(i,j,k)                         &
     &                      * FV_sec_theta_latitude(i,j)
                  soln_prev(i,j,k) = soln(i,j,k)
                End Do
              End Do
            End Do

          Else
! add on constant terms and save solution from previous timestep.
            Do k = 1, model_levels
              Do j = 1, rows
                Do i = 1, row_length
                RHS(i,j,k) = r(i,j,k) * FV_sec_theta_latitude(i,j)      &
     &                       + HM_C4(i,j,k) * soln(i,j,k)
                soln_prev(i,j,k) = soln(i,j,k)
                End Do
              End Do
            End Do

! Add on z, y terms
            Do k = 1, model_levels
! Include Vertical derivative term in right-hand-side
              If ( k  ==  1) Then
                factor_1 = 1. / (eta_theta_levels(k) -                  &
     &                         eta_theta_levels(k-1))
                factor_2 = 1. / (eta_rho_levels(k+1) -                  &
     &                         eta_rho_levels(k))

                Do j = 1, rows
                  Do i = 1, row_length

! Calculate upper derivative
                    term_1 = ( soln(i,j,k+1) - soln(i,j,k) )            &
     &                      * factor_2

! Calculate lower derivative. This is zero by boundary condition.

                  RHS(i,j,k) = RHS(i,j,k) -                             &
     &                              ( factor_1 *                        &
     &                               term_1 * HM_Czz(i,j,k)             &
     &                              + HM_C3(i,j,k) *                    &
     &                               term_1 * HM_Cz(i,j,k) )
                  End Do
                End Do

              Else if (k  ==  model_levels) Then
                factor_1 = 1. / (eta_theta_levels(k) -                  &
     &                         eta_theta_levels(k-1))
                factor_3 = 1. / (eta_rho_levels(k) -                    &
     &                         eta_rho_levels(k-1))

                Do j = 1, rows
                  Do i = 1, row_length

! Calculate upper derivative. This is zero by boundary condition.

! Calculate lower derivative.
                    term_2 = ( soln(i,j,k) - soln(i,j,k-1) )            &
     &                          * factor_3

                    RHS(i,j,k) = RHS(i,j,k) +                           &
     &                              ( factor_1 *                        &
     &                                term_2 * HM_Czz(i,j,k-1)          &
     &                               - HM_C3(i,j,k) *                   &
     &                                term_2 * HM_Cz(i,j,k-1) *         &
     &                                weight_lower(i,j,k) )

                  End Do
                End Do

              Else
                factor_1 = 1. / (eta_theta_levels(k) -                  &
     &                         eta_theta_levels(k-1))
                factor_2 = 1. / (eta_rho_levels(k+1) -                  &
     &                         eta_rho_levels(k))
                factor_3 = 1. / (eta_rho_levels(k) -                    &
     &                         eta_rho_levels(k-1))
                Do j = 1, rows
!dir$ split
!dir$ unroll 4
                  Do i = 1, row_length

                    RHS(i,j,k) = RHS(i,j,k) -                           &
     &                       ( soln(i,j,k+1) - soln(i,j,k) )            &
     &                        * factor_2 *                              &
     &                              ( factor_1 *                        &
     &                                 HM_Czz(i,j,k)                    &
     &                               + HM_C3(i,j,k) *                   &
     &                                 HM_Cz(i,j,k) *                   &
     &                                weight_upper(i,j,k) )
                  Enddo
                Enddo
                Do j = 1, rows
!dir$ split
!dir$ unroll 4
                  Do i = 1, row_length
                    RHS(i,j,k) = RHS(i,j,k) +                           &
     &                       ( soln(i,j,k) - soln(i,j,k-1) )            &
     &                        * factor_3 *                              &
     &                              ( factor_1 *                        &
     &                                  HM_Czz(i,j,k-1)                 &
     &                               - HM_C3(i,j,k) *                   &
     &                                 HM_Cz(i,j,k-1) *                 &
     &                                weight_lower(i,j,k)  )

                  End Do
                End Do

              End If

              Do j = j_start, j_end
!dir$ split
!dir$ unroll 4
                Do i = 1, row_length
                  RHS(i,j,k) = RHS(i,j,k) + (                           &
     &                       HM_Cyy1(i,j,k)*HM_Cyy2(i,j,k)              &
     &                       *(soln(i,j,k) - soln(i,j+1,k)) +           &
     &                       HM_Cyy1(i,j-1,k)*HM_Cyy2(i,j-1,k)          &
     &                       *(soln(i,j,k)- soln(i,j-1,k)))             &
     &                       *FV_sec_theta_latitude(i,j)
                End Do
              End Do

            End Do ! end loop over levels

            If(at_extremity(PSouth))then
              Do k = 1, model_levels
                Do i = 1,row_length
                sum_s_component(i,k)= HM_Cyy1(i,1,k)*HM_Cyy2(i,1,k)     &
     &                                *(soln(i,1,k)-soln(i,2,k))
                End Do
              End Do
            End If
            If(at_extremity(PNorth))then
              Do k = 1, model_levels
                Do i = 1,row_length
                sum_n_component(i,k)= - HM_Cyy1(i,rows-1,k)*            &
     &                                  HM_Cyy2(i,rows-1,k)             &
     &                          *(soln(i,rows-1,k)-soln(i,rows,k))
                End Do
              End Do
            End If

            If(at_extremity(PSouth))then
              Call gcg_rvecsumr(row_length,row_length,1,                &
     &                        model_levels,sum_s_component,             &
     &                        proc_row_group,istat,sum_s)
              Do k = 1, model_levels
                sum_s(k) = sum_s(k) * FV_sec_theta_latitude(1,1)        &
     &                     /global_row_length
                Do i = 1,row_length
                  RHS(i,1,k) = RHS(i,1,k) + sum_s(k)
                End Do
              End Do
            End If

            If(at_extremity(PNorth))then
              Call gcg_rvecsumr(row_length,row_length,1,                &
     &                        model_levels,sum_n_component,             &
     &                        proc_row_group,istat,sum_n)
              Do k = 1, model_levels
              sum_n(k) = sum_n(k) * FV_sec_theta_latitude(1,rows)       &
     &                     /global_row_length
                Do i = 1,row_length
                  RHS(i,rows,k) = RHS(i,rows,k) + sum_n(k)
                End Do
              End Do
            End If


           Do k = 1, model_levels
             Do j = j_start, j_end
               Do i = 1, row_length
                 RHS_c(i,j,k) = RHS(i,j,k) + (                          &
     &                     HM_Cxx1(i,j,k)*HM_Cxx2(i,j,k)                &
     &                     *(soln(i,j,k) - soln(i+1,j,k)) +             &
     &                     HM_Cxx1(i-1,j,k)*HM_Cxx2(i-1,j,k)            &
     &                     *(soln(i,j,k)- soln(i-1,j,k)))               &
     &                     *FV_sec_theta_latitude(i,j)
               End Do
             End Do
           End Do

          End If

        End If ! on first timestep or not

!-----------------------------------------------------------------------
!     Section 2. Perform ADI sweep in lambda direction.
!                Only implemented for cyclic domains.
!-----------------------------------------------------------------------

! copy rhs into L_of_Soln and use this to pass into subroutine.
! required since RHS passed into subroutine is changed on output and
! this routine requires an unchanged RHS later on.
!         Do k = 1, model_levels
!           Do j = j_start, j_end
!             Do i= 1, row_length
!               L_of_soln(i,j,k)= RHS(i,j,k)
!             End Do
!           End Do
!         End Do

! DEPENDS ON: mpp_tri_solve_exec
          Call mpp_tri_solve_exec(                                      &
     &                     row_length, rows, model_levels,              &
     &                     offx, offy, j_start, j_end,                  &
     &                     n_proc, n_procx, me,                         &
     &                     proc_row_group,                              &
     &                     a_plus_x, a_minus_x,                         &
     &                     factor_forward, factor_backward,             &
     &                     recip_a_central_x,                           &
     &                     bv_a_matrix_0, bv_a_matrix_np1,              &
     &                     recip_bv_a_matrix_diag, bv_a_matrix_sup,     &
     &                     bv_factor_forward, bv_factor_backward,       &
     &                     bv_soln_n_term,                              &
     &                     bv_soln_1_term1,                             &
     &                     bv_soln_1_term2,                             &
     &                     L_of_soln, RHS_c, soln)

! swop soln on all procesors
! DEPENDS ON: swap_bounds
        Call swap_bounds(                                               &
     &                   soln,row_length,rows,model_levels,             &
     &                   offx,offy,fld_type_p,.FALSE.)

! Modify RHS to contain x term correction

       Do k = 1, model_levels
         Do j = j_start, j_end
           Do i = 1, row_length
             RHS(i,j,k) = RHS_c(i,j,k) + (                              &
     &                       HM_Cxx1(i,j,k)*HM_Cxx2(i,j,k)              &
     &                       *(soln(i,j,k) - soln(i+1,j,k)) +           &
     &                       HM_Cxx1(i-1,j,k)*HM_Cxx2(i-1,j,k)          &
     &                       *(soln(i,j,k)- soln(i-1,j,k)))             &
     &                       *FV_sec_theta_latitude(i,j)
           End Do
         End Do
       End Do

!-----------------------------------------------------------------------
!     Section 3. Perform ADI sweep in phi direction.
!-----------------------------------------------------------------------

      If (GCR_precon_option  ==  vert_plus_xyz_ADI_precon .or.          &
     &    GCR_precon_option  ==  xyz_ADI_precon )                       &
     &    Then
          Do k = 1, model_levels
            Do i= 1, row_length
              Do j = j_start, j_end
                soln(i,j,k)= RHS(i,j,k)
              End Do
            End Do
          End Do

! Solve tridiagonal system.
! solution is in Soln
! reduce matrix to upper diagonal form.

!          Do k= 1, model_levels
!            Do j = 3, rows - 1
!              Do i = 1, row_length
!                soln(i,j,k)= soln(i,j,k)-
!     &                       factor_y(i,j,k)*soln(i,j-1,k)
!              End Do
!            End Do
!          End Do
! DEPENDS ON: solve_forward_y2
          Call solve_forward_y2(soln, factor_y,                         &
     &                        row_length, rows, model_levels,           &
     &                        offx, offy)

! Back substitute to get solution.

          If(at_extremity(PNorth)) Then
            Do k= 1, model_levels
              Do i = 1, row_length
                Soln(i,rows-1,k) = a0_y(i,rows-1,k) *                   &
     &                             soln(i,rows-1,k)
              End Do
            End Do
          End If

!          Do k= 1, model_levels
!            Do j = rows-2, 2, -1
!              Do i = 1, row_length
!                Soln(i,j,k) = a0_y(i,j,k) * ( soln(i,j,k) -
!     &                              a1_y(i,j,k) * Soln(i,j+1,k) )
!              End Do
!            End Do
!          End Do

! DEPENDS ON: solve_backward_y1
          Call solve_backward_y1(soln, a0_y, a1_y,                      &
     &                         row_length, rows, model_levels,          &
     &                         offx, offy)

! set soln at poles to zero
          If(at_extremity(PNorth)) Then
            Do k= 1, model_levels
              Do i = 1, row_length
                Soln(i,rows,k) = 0.
              End Do
            End Do
          End If
          If(at_extremity(PSouth)) Then
            Do k= 1, model_levels
              Do i = 1, row_length
                Soln(i,1,k) = 0.
              End Do
            End Do
          End If

! swop soln on all procesors
! DEPENDS ON: swap_bounds
          CALL swap_bounds(                                             &
     &                     soln,row_length,rows,model_levels,           &
     &                     offx,offy,fld_type_p,.FALSE.)

! add on y direction correction
          Do k = 1, model_levels
!CDIR unroll=4
            Do j = j_start, j_end
              Do i = 1, row_length
                RHS(i,j,k) = RHS(i,j,k) + (                             &
     &                     HM_Cyy1(i,j,k)*HM_Cyy2(i,j,k)                &
     &                       *(soln(i,j,k) - soln(i,j+1,k)) +           &
     &                       HM_Cyy1(i,j-1,k)*HM_Cyy2(i,j-1,k)          &
     &                       *(soln(i,j,k)- soln(i,j-1,k)))             &
     &                       *FV_sec_theta_latitude(i,j)
              End Do
            End Do
          End Do

          If(at_extremity(PSouth))then
            Do k = 1, model_levels
              Do i = 1,row_length
                sum_s_component(i,k)= HM_Cyy1(i,1,k)*HM_Cyy2(i,1,k)     &
     &                                *(soln(i,1,k)-soln(i,2,k))
              End Do
            End Do
          End If
          If(at_extremity(PNorth))then
            Do k = 1, model_levels
              Do i = 1,row_length
                sum_n_component(i,k)= - HM_Cyy1(i,rows-1,k)*            &
     &                                  HM_Cyy2(i,rows-1,k)             &
     &                           *(soln(i,rows-1,k)-soln(i,rows,k))
              End Do
            End Do
          End If

          If(at_extremity(PSouth))then
            Call gcg_rvecsumr(row_length,row_length,1,                  &
     &                        model_levels,sum_s_component,             &
     &                        proc_row_group,istat,sum_s)
            Do k = 1, model_levels
              sum_s(k) = sum_s(k) * FV_sec_theta_latitude(1,1)          &
     &                     /global_row_length
              Do i = 1,row_length
                RHS(i,1,k) = RHS(i,1,k) + sum_s(k)
              End Do
            End Do
          End If

          If(at_extremity(PNorth))then
            Call gcg_rvecsumr(row_length,row_length,1,                  &
     &                        model_levels,sum_n_component,             &
     &                        proc_row_group,istat,sum_n)
            Do k = 1, model_levels
              sum_n(k) = sum_n(k) * FV_sec_theta_latitude(1,rows)       &
     &                     /global_row_length
              Do i = 1,row_length
                RHS(i,rows,k) = RHS(i,rows,k) + sum_n(k)
              End Do
            End Do
          End If

        End If ! on precon option
!-----------------------------------------------------------------------
!     Section 4. Perform ADI sweep in vertical direction.
!-----------------------------------------------------------------------

! Solve tridiagonal system.
! solution is in Soln
! reduce matrix to upper diagonal form.

        Do j = 1, rows
!CDIR unroll=8
          Do k= 2, model_levels
            Do i = 1, row_length
              RHS(i,j,k) = RHS(i,j,k) - factor_z(i,j,k)*RHS(i,j,k-1)
            End Do
          End Do
        End Do

! Back substitute to get solution.

        Do j = 1, rows
          Do i = 1, row_length
            Soln(i,j,model_levels) = a0_z(i,j,model_levels) *           &
     &                               RHS(i,j,model_levels)
          End Do
        End Do
        Do j = 1, rows
!CDIR unroll=8
          Do k= model_levels-1, 1, -1
            Do i = 1, row_length
              Soln(i,j,k) = a0_z(i,j,k) * ( RHS(i,j,k) -                &
     &                            a1_z(i,j,k) * Soln(i,j,k+1) )
            End Do
          End Do
        End Do

! add on correction to previous solution
        If (pseudo_timestep_number  >   1) Then
          Do j = 1, rows
!CDIR unroll=8
            Do k= 1, model_levels
              Do i = 1, row_length
                Soln(i,j,k) = Soln(i,j,k) + Soln_prev(i,j,k)
              End Do
            End Do
          End Do
        End If

      End Do ! end loop over number of pseudo timesteps

!     end of routine GCR_precon_ADI_exec

      return
      END SUBROUTINE GCR_precon_ADI_exec_trisolve

