#if defined(A10_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine GCR_precon_ADI_exec

      Subroutine GCR_precon_ADI_exec(                                   &
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
     &                     a0_x, a1_x, factor_x,                        &
     &                     F_vector_x, G_vector_x,                      &
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
! 5.0                  This routine introduced           Andy Malcolm
! 5.2      11/7/00     Modify calls to solver routines.  P.Burton
! 5.3     23/10/01     Change from magic numbers         Andy Malcolm
! 5.3      19/10/01    Use appropriate gcg routines.   S. Cusack
!  6.0  18/08/03  NEC SX-6 optimisation - change dimension to avoid
!                 bank conflict.  R Barnes & J-C Rioual.
! 6.1     30/07/04     Optimisation to reduce message passing
!                                                           Andy Malcolm
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

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.
      Real                                                              &
     &  Soln(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

      Real                                                              &
     &solno(1-offx:row_length+offx+1,1-offy:rows+offy,model_levels)

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

      Integer                                                           &
     &  i_start                                                         &
     &, i_end                                                           &
     &, j_start                                                         &
     &, j_end                                                           &
     &, istat                                                           &
     &, ibase                                                           &
     &, len

      Real                                                              &
     &  recip_ADI_pseudo_timestep

! Local arrays.

      Real                                                              &
     &  RHS(row_length,rows,model_levels)                               &
     &, a0_z(row_length,rows,model_levels)                              &
     &, a1_z(row_length,rows,model_levels)                              &
     &, factor_z(row_length,rows,model_levels)                          &
     &, a0_y(row_length,rows,model_levels)                              &
     &, a1_y(row_length,rows,model_levels)                              &
     &, factor_y(row_length,rows,model_levels)                          &
     &, a0_x(row_length+1,rows,model_levels)                            &
     &, a1_x(row_length+1,rows,model_levels)                            &
     &, factor_x(row_length+1,rows,model_levels)                        &
     &, F_vector_x(row_length+1,rows,model_levels)                      &
     &, G_vector_x(2,rows,model_levels)                                 &
     &, soln_prev(row_length,rows,model_levels)                         &
     &, soln_end(rows,model_levels)                                     &
     &, L_of_soln(row_length,rows,model_levels)

      Real                                                              &
     &  sum_n(model_levels)                                             &
     &, sum_s(model_levels)                                             &
     &, sum_n_component(row_length,model_levels)                        &
     &, sum_s_component(row_length,model_levels)

#include "fldtype.h"
#include "precon.h"

!    External routines.
      External                                                          &
     &  GCG_rvecsumr, GCG_rbcast

!-----------------------------------------------------------------------
!     Section 1. Set initial guess to solution and initialise
!                variables.
!-----------------------------------------------------------------------

      j_start = 1
      j_end = rows
      If (at_extremity(PSouth) ) Then
        j_start = 2
      End If
      If (at_extremity(PNorth) ) Then
        j_end = rows-1
      End If

      recip_ADI_pseudo_timestep = 1./ ADI_pseudo_timestep

      Do pseudo_timestep_number = 1, GCR_n_ADI_pseudo_timesteps

! Form matrix
        If (pseudo_timestep_number  ==  1) Then
          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                RHS(i,j,k) = r(i,j,k) * FV_sec_theta_latitude(i,j)
              End Do
            End Do
          End Do

        Else
! swop soln on all procesors
! DEPENDS ON: swap_bounds
          Call swap_bounds(                                             &
     &                     soln,row_length,rows,model_levels,           &
     &                     offx,offy,fld_type_p,.FALSE.)

        If (GCR_adi_add_full_soln) Then
! DEPENDS ON: gcr_elliptic_operator
          Call GCR_Elliptic_Operator(                                   &
     &                          soln, row_length,                       &
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
     &                              ( factor_1 * term_1 *               &
     &                               HM_Czz(i,j,k) + HM_C3(i,j,k) *     &
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
                  Do i = 1, row_length

! Calculate upper derivative
                    term_1 = ( soln(i,j,k+1) - soln(i,j,k) )            &
     &                        * factor_2

! Calculate lower derivative
                    term_2 = ( soln(i,j,k) - soln(i,j,k-1) )            &
     &                        * factor_3

! Calculate operator
                    RHS(i,j,k) = RHS(i,j,k) -                           &
     &                              ( factor_1 *                        &
     &                              ( term_1 * HM_Czz(i,j,k) -          &
     &                                term_2 * HM_Czz(i,j,k-1) )        &
     &                               + HM_C3(i,j,k) *                   &
     &                              ( term_1 * HM_Cz(i,j,k) *           &
     &                                weight_upper(i,j,k) +             &
     &                                term_2 * HM_Cz(i,j,k-1) *         &
     &                                weight_lower(i,j,k) ) )

                  End Do
                End Do

              End If

              Do j = j_start, j_end
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
#if defined(REPROD)
              Call gcg_rvecsumr(row_length,row_length,1,                &
     &                        model_levels,sum_s_component,             &
     &                        proc_row_group,istat,sum_s)
#else
              Call gcg_rvecsumf(row_length,row_length,1,                &
     &                        model_levels,sum_s_component,             &
     &                        proc_row_group,istat,sum_s)
#endif
              Do k = 1, model_levels
                sum_s(k) = sum_s(k) * FV_sec_theta_latitude(1,1)        &
     &                     /global_row_length
                Do i = 1,row_length
                  RHS(i,1,k) = RHS(i,1,k) + sum_s(k)
                End Do
              End Do
            End If

            If(at_extremity(PNorth))then
#if defined(REPROD)
              Call gcg_rvecsumr(row_length,row_length,1,                &
     &                        model_levels,sum_n_component,             &
     &                        proc_row_group,istat,sum_n)
#else
              Call gcg_rvecsumf(row_length,row_length,1,                &
     &                        model_levels,sum_n_component,             &
     &                        proc_row_group,istat,sum_n)
#endif
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
                  RHS(i,j,k) = RHS(i,j,k) + (                           &
     &                       HM_Cxx1(i,j,k)*HM_Cxx2(i,j,k)              &
     &                       *(soln(i,j,k) - soln(i+1,j,k)) +           &
     &                       HM_Cxx1(i-1,j,k)*HM_Cxx2(i-1,j,k)          &
     &                       *(soln(i,j,k)- soln(i-1,j,k)))             &
     &                       *FV_sec_theta_latitude(i,j)
                End Do
              End Do
            End Do

          End If

        End If ! on first timestep or not

!-----------------------------------------------------------------------
!     Section 2. Perform ADI sweep in lambda direction.
!                Only implemented for cyclic domains.
!-----------------------------------------------------------------------

        Do k = 1, model_levels
          Do j = j_start, j_end
            Do i= 1, row_length
              soln(i,j,k)= RHS(i,j,k)
            End Do
          End Do
        End Do

! solve cyclic tridiagonal system for solution
! reduce matrix to upper diagonal form.

      IF (n_procx >  1) then
!        Do i= 2, row_length - 1
!          Do k = 1, model_levels
!            Do j = j_start, j_end
!              soln(i,j,k)= soln(i,j,k)-
!     &                     factor_x(i,j,k)*soln(i-1,j,k)
!            End Do
!          End Do
!        End Do

        do k=1,model_levels
        do j=1-offy,rows+offy
        do i=1-offx,row_length+offx
        solno(i,j,k)=soln(i,j,k)
        enddo
        enddo
        enddo

! DEPENDS ON: solve_forward_x2
        Call solve_forward_x2(solno, factor_x,                          &
     &                      row_length, rows, model_levels,             &
     &                      j_start, j_end, offx, offy)



! back substitute to get solution.

        If (at_extremity(PEast)) Then
          Do k = 1, model_levels
            Do j = j_start, j_end
              solno(row_length-1,j,k) = solno(row_length-1,j,k) *       &
     &                                 a0_x(row_length-1,j,k)
            End Do
          End Do
        End If

!        Do i= row_length - 2, 1, -1
!          Do k = 1, model_levels
!            Do j = j_start, j_end
!              soln(i,j,k) = a0_x(i,j,k)*(soln(i,j,k) -
!     &                                   a1_x(i,j,k)*soln(i+1,j,k))
!            End Do
!          End Do
!        End Do

! DEPENDS ON: solve_backward_x1
        Call solve_backward_x1(solno, a0_x, a1_x,                       &
     &                       row_length, rows, model_levels,            &
     &                       j_start, j_end, offx, offy)

        do k=1,model_levels
        do j=1-offy,rows+offy
        do i=1-offx,row_length+offx
        soln(i,j,k)=solno(i,j,k)
        enddo
        enddo
        enddo


! solve for last element of field.
! first broadcast soln(1,j,k) to all row processors from
! western most (first on row) processor
        ibase = (me/n_procx) * n_procx
        If (at_extremity(PWest)) Then
          Do k = 1, model_levels
            Do j = j_start, j_end
              soln_end(j,k) = soln(1,j,k)
            End Do
          End Do
        End If

        If (n_procx  >   1) Then
          len = model_levels * rows
          Call gcg_rbcast(201, len, ibase,                              &
     &                    proc_row_group, istat, soln_end)
        End If

        i_end = row_length
        If (at_extremity(PEast)) Then
          Do k = 1, model_levels
            Do j = j_start, j_end
              soln(row_length,j,k) = ( soln(row_length,j,k) -           &
     &                                G_vector_x(1,j,k)*soln_end(j,k)-  &
     &                                G_vector_x(2,j,k) *               &
     &                                soln(row_length-1,j,k) )          &
     &                                *a0_x(row_length,j,k)
            End Do
          End Do
          i_end = row_length - 1
        End If

! solve for all other elements of field.
! first broadcast soln(row_length,j,k) to all row processors from
! eastern most (last on row) processor
        ibase = (me/n_procx) * n_procx
        If (at_extremity(PEast)) Then
          Do k = 1, model_levels
            Do j = j_start, j_end
              soln_end(j,k) = soln(row_length,j,k)
            End Do
          End Do
        End If

        If (n_procx  >   1) Then
          len = model_levels * rows
          Call gcg_rbcast(201, len, ibase+n_procx-1,                    &
     &                    proc_row_group, istat, soln_end)
        End If

        Do k = 1, model_levels
          Do j = j_start, j_end
            Do i= 1, i_end
!              soln(i,j,k) = soln(i,j,k) - F_vector_x(i,j,k) *
!     &                                    soln(row_length,j,k)
              soln(i,j,k) = soln(i,j,k) - F_vector_x(i,j,k) *           &
     &                                    soln_end(j,k)
            End Do
          End Do
        End Do
      ELSE
! solve cyclic tridiagonal system for solution
! reduce matrix to upper diagonal form.

            Do i= 2, row_length - 1
              Do k = 1, model_levels
                Do j = j_start, j_end
             soln(i,j,k) = soln(i,j,k) - factor_x(i,j,k)*soln(i-1,j,k)
                End Do
              End Do
            End Do
! back substitute to get solution.

            Do k = 1, model_levels
              Do j = j_start, j_end
                soln(row_length-1,j,k) = soln(row_length-1,j,k) *       &
     &                                  a0_x(row_length-1,j,k)
              End Do
            End Do
            Do i= row_length - 2, 1, -1
              Do k = 1, model_levels
                Do j = j_start, j_end
                  soln(i,j,k) = a0_x(i,j,k) * (soln(i,j,k) -            &
     &                                    a1_x(i,j,k)*soln(i+1,j,k))
                End Do
              End Do
            End Do

! solve for last element of field.

            Do k = 1, model_levels
              Do j = j_start, j_end
                soln(row_length,j,k) = ( soln(row_length,j,k) -         &
     &                             G_vector_x(1,j,k) * soln(1,j,k) -    &
     &                                G_vector_x(2,j,k) *               &
     &                                soln(row_length-1,j,k) )          &
     &                             *a0_x(row_length,j,k)
              End Do
            End Do

! solve for all other elements of field.

            Do k = 1, model_levels
              Do j = j_start, j_end
                Do i= 1, row_length - 1
                  soln(i,j,k) = soln(i,j,k) - F_vector_x(i,j,k) *       &
     &                                       soln(row_length,j,k)
                End Do
              End Do
            End Do

        ENDIF   !n_procx

! swop soln on all procesors
! DEPENDS ON: swap_bounds
        CALL swap_bounds(                                               &
     &                   soln,row_length,rows,model_levels,             &
     &                   offx,offy,fld_type_p,.FALSE.)

! Modify RHS to contain x term correction

        Do k = 1, model_levels
          Do j = j_start, j_end
            Do i = 1, row_length
              RHS(i,j,k) = RHS(i,j,k) + (                               &
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

        If (GCR_precon_option  ==  vert_plus_xyz_ADI_precon .or.        &
     &      GCR_precon_option  ==  xyz_ADI_precon)                      &
     &    Then
          Do k = 1, model_levels
            Do i= 1, row_length
              Do j = j_start, j_end
                soln(i,j,k)= RHS(i,j,k)
              End Do
            End Do
          End Do

      IF (n_procy >  1) then
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

      ELSE
! Solve tridiagonal system.
! solution is in Soln
! reduce matrix to upper diagonal form.

            Do k= 1, model_levels
              Do j = 3, rows - 1
                Do i = 1, row_length
              soln(i,j,k) = soln(i,j,k) - factor_y(i,j,k)*soln(i,j-1,k)
                End Do
              End Do
            End Do

! Back substitute to get solution.

            Do k= 1, model_levels
              Do i = 1, row_length
                Soln(i,rows-1,k) = a0_y(i,rows-1,k) *                   &
     &                                   soln(i,rows-1,k)
              End Do
            End Do
            Do k= 1, model_levels
              Do j = rows-2, 2, -1
                Do i = 1, row_length
                  Soln(i,j,k) = a0_y(i,j,k) * ( soln(i,j,k) -           &
     &                               a1_y(i,j,k) * Soln(i,j+1,k) )
                End Do
              End Do
            End Do

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

      ENDIF !n_procy
! swop soln on all procesors
! DEPENDS ON: swap_bounds
          Call swap_bounds(                                             &
     &                     soln,row_length,rows,model_levels,           &
     &                     offx,offy,fld_type_p,.FALSE.)

! add on y direction correction
          Do k = 1, model_levels
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
#if defined(REPROD)
            Call gcg_rvecsumr(row_length,row_length,1,                  &
     &                        model_levels,sum_s_component,             &
     &                        proc_row_group,istat,sum_s)
#else
            Call gcg_rvecsumf(row_length,row_length,1,                  &
     &                        model_levels,sum_s_component,             &
     &                        proc_row_group,istat,sum_s)
#endif
            Do k = 1, model_levels
              sum_s(k) = sum_s(k) * FV_sec_theta_latitude(1,1)          &
     &                     /global_row_length
              Do i = 1,row_length
                RHS(i,1,k) = RHS(i,1,k) + sum_s(k)
              End Do
            End Do
          End If

          If(at_extremity(PNorth))then
#if defined(REPROD)
            Call gcg_rvecsumr(row_length,row_length,1,                  &
     &                        model_levels,sum_n_component,             &
     &                        proc_row_group,istat,sum_n)
#else
            Call gcg_rvecsumf(row_length,row_length,1,                  &
     &                        model_levels,sum_n_component,             &
     &                        proc_row_group,istat,sum_n)
#endif
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

        Do k= 2, model_levels
          Do j = 1, rows
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
        Do k= model_levels-1, 1, -1
          Do j = 1, rows
            Do i = 1, row_length
              Soln(i,j,k) = a0_z(i,j,k) * ( RHS(i,j,k) -                &
     &                            a1_z(i,j,k) * Soln(i,j,k+1) )
            End Do
          End Do
        End Do

! add on correction to previous solution
        If (pseudo_timestep_number  >   1) Then
          Do k= 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                Soln(i,j,k) = Soln(i,j,k) + Soln_prev(i,j,k)
              End Do
            End Do
          End Do
        End If

      End Do ! end loop over number of pseudo timesteps

!     end of routine GCR_precon_ADI_exec

      return
      END SUBROUTINE GCR_precon_ADI_exec

#endif
