
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Mono_Conserv_qcon.

      Subroutine Mono_Conserv_qcon(                                     &
     &                        Data_out_high, Data_out_mono,             &
     &                        Ext_Data, rho_dr_n, rho_dr_np1,           &
     &                        dim_i_in, dim_j_in, dim_k_in,             &
     &                        dim_i_out, dim_j_out, dim_k_out,          &
     &                        row_length, rows,                         &
     &                        lambda_p, phi_p, L_regular,               &
     &                        cos_latitude, off_x, off_y,               &
     &                        n_proc, halo_i, halo_j, model_domain,     &
     &                        me, proc_row_group, proc_col_group,       &
     &                        halo_data_out_i, halo_data_out_j,         &
     &                        Data_out, conserv_fail)

! Purpose:
!          Ensures scheme is conservative, if desired. This is only
!          possible if the input data points and the output data
!          points are the same.
!
! Method:
!          Is described in ;
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson
!
!          and is based on Priestley, 1993. The full reference
!          for which can be found in the above documentation.
!
! Original Programmer:  Nigel Wood
! Current code owner: Andrew J. Malcolm
!
! History:
! Version  Date     Comment
! -------  -----    -------
! 5.2     02/11/00  New routine. Improved conservation calculation
!                                                             A. Malcolm
! 5.3     19/10/01  Use appropriate gcg routines.   S. Cusack
! 5.5      14/02/03  improve tracer conservation properties A.Malcolm
! 6.0      03/11/03  add compiler directive for bit-comparison A.Malcolm
!  6.0  03/10/03  NEC SX-6 optimisation - R Barnes & J-C Rioual.
!  6.1  17/08/04  NEC SX-6 optimisation S.S.Wilson
! 6.1      13/05/04  remove unnecessary multiplications       A.Malcolm
! 6.2      21/10/05  Replace GSYNC with SSYNC. P.Selwood.
!  6.2  25/12/05  Variable resolution changes            Yongming Tang
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
     &  L_regular  ! false if variable resolution

      Integer                                                           &
     &  dim_i_in                                                        &
                    ! Dimension of input Data arrays in i direction.
     &, dim_j_in                                                        &
                    ! Dimension of input Data arrays in j direction.
     &, dim_k_in                                                        &
                    ! Dimension of input Data arrays in k direction.
     &, dim_i_out                                                       &
                    ! Dimension of output Data arrays in i direction.
     &, dim_j_out                                                       &
                    ! Dimension of output Data arrays in j direction.
     &, dim_k_out                                                       &
                    ! Dimension of output Data arrays in k direction.
     &, row_length                                                      &
                     ! number of points on a row for lambda_p, phi_p
     &, rows                                                            &
                     ! number of rows for phi_p
     &, n_proc                                                          &
                    ! Total number of processors
     &, me                                                              &
                    ! My processor number
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &, proc_col_group                                                  &
                       ! Group id for processors on the same column
     &, halo_i                                                          &
                    ! Size of halo in i direction.
     &, halo_j                                                          &
                    ! Size of halo in j direction.
     &, off_x                                                           &
                    ! Size of small halo in i direction.
     &, off_y                                                           &
                    ! Size of small halo in j direction.
     &, halo_data_out_i                                                 &
                        ! size of data out halo in i direction
     &, halo_data_out_j ! size of data out halo in j direction

      Integer :: model_domain      ! holds integer code for model domain   

      Real                                                              &
                  !  VarRes horizontal co-ordinate spacing.
     &  lambda_p (1 - halo_i : dim_i_in + halo_i)                       &
     &, phi_p (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)


      Real                                                              &
     &  Data_out_high (dim_i_out, dim_j_out, dim_k_out)                 &
                            ! data interpolated by high order scheme
     &, Data_out_mono (dim_i_out, dim_j_out, dim_k_out)                 &
                            ! data interpolated by monotone scheme
     &, rho_dr_n (dim_i_in, dim_j_in, dim_k_in)                         &
                             ! Weighted vertical layer thickness
                             ! of  data times density at time level n.
     &, rho_dr_np1 (dim_i_in, dim_j_in, dim_k_in)                       &
                             ! Weighted vertical layer thickness
                             ! of  data times density at time level n+1.
     &, Ext_Data (1-halo_i:dim_i_in+halo_i+1,                           &
     &            1-halo_j:dim_j_in+halo_j,-1:dim_k_in+2)
                              !data to be  interpolated

      Real                                                              &
     &  cos_latitude (1-off_x:dim_i_in+off_x,                           &
     &                1-off_y:dim_j_in+off_y)

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
           ! data at desired locations.
     &  Data_out (1-halo_data_out_i:dim_i_out+halo_data_out_i,          &
     &            1-halo_data_out_j:dim_j_out+halo_data_out_j,          &
     &            dim_k_out)

      Logical                                                           &
     &  conserv_fail ! true if conservation enforcement failed.

! Local Variables.

! scalars

      Integer                                                           &
     &  i, j, k                                                         &
                       ! Loop indices
     &, count

      Logical                                                           &
     &  L_conservation_enforced

      Real                                                              &
     &  C_star                                                          &
     &, Integral_beta                                                   &
     &, mass_n(dim_i_in,dim_j_in,dim_k_in)                              &
     &, mass_np1(dim_i_in,dim_j_in,dim_k_in)                            &
     &, surplus                                                         &
     &, average_alpha                                                   &
     &,   temp1(dim_i_out,dim_j_out)                                    &
     &,   temp2(dim_i_out,dim_j_out)                                    &
     &,   temp1_rows(dim_j_out)                                         &
     &,   temp2_rows(dim_j_out)                                         &
     &,   temp1_sum(1)                                                  &
     &,   temp2_sum(1)                                                  &
     &,   temp3_sum(1)

! arrays

      Logical                                                           &
     &  L_alpha_set (dim_i_out, dim_j_out, dim_k_out)

      Real                                                              &
     &  alpha_cons (dim_i_in, dim_j_in, dim_k_in)                       &
     &, beta (dim_i_in, dim_j_in, dim_k_in)                             &
     &, alpha_max (dim_i_out, dim_j_out, dim_k_out)                     &
     &, lamxphi (dim_i_in, dim_j_in)

      real fielda(dim_i_in, dim_j_in,2)
      real work(dim_i_out, dim_j_out, dim_k_out)
      real sum_rows(dim_j_in,2), sum_all(2)

      Integer                                                           &
     &  istat

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

! External Routines: None
! subroutines: None
! Functions: None

! ----------------------------------------------------------------------
!  Section 1.  Initial settings.
!              Set alpha_max to 1 as both input fields are monotone.
! ----------------------------------------------------------------------

      conserv_fail=.false.
      alpha_max(:,:,:) = 1.0

      if(L_regular) then
        Do k = 1, dim_k_in
          Do j = 1, dim_j_in
            Do i = 1, dim_i_in
              mass_n(i,j,k)   =   cos_latitude(i,j) * rho_dr_n(i,j,k)
              mass_np1(i,j,k) =   cos_latitude(i,j) * rho_dr_np1(i,j,k)
            End Do
          End Do
        End Do
      else ! variable resolution
        Do j = 1, dim_j_in
          Do i = 1, dim_i_in
            lamxphi(i,j) = (lambda_p(i+1) - lambda_p(i)) *              &
     &                        (phi_p(i,j+1) - phi_p(i,j))
          End Do
        End Do
       Do k = 1, dim_k_in
          Do j = 1, dim_j_in
            Do i = 1, dim_i_in
              mass_n(I,j,k)   = lamxphi(i,j) *  cos_latitude(i,j)       &
     &                          * rho_dr_n(i,j,k)
              mass_np1(i,j,k) = lamxphi(i,j) *  cos_latitude(i,j)       &
     &                          * rho_dr_np1(i,j,k)
            End Do
          End Do
        End Do
      endif       ! L_regular

! ----------------------------------------------------------------------
!  Section 2.   Enforce Conservation.
!               Section 2. in Priestley 1993.
!               Form output data from the new alpha values.
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
!  section 2.1  Calculate C_star, Beta, and integral of beta.
! ----------------------------------------------------------------------

        k = 1
          Do j = 1, dim_j_in
            Do i = 1, dim_i_in
              fielda(i,j,1) = Ext_Data (i,j,k) * mass_n(i,j,k) -        &
     &                        Data_out_mono(i,j,k) * mass_np1(i,j,k)
              beta (i,j,k) = (Data_out_high(i,j,k) -                    &
     &                        Data_out_mono(i,j,k)) * mass_np1(i,j,k)
              fielda(i,j,2) = beta(i,j,k)
            End Do
          End Do

        Do k = 2, dim_k_in
          Do j = 1, dim_j_in
            Do i = 1, dim_i_in
              work(i,j,k) =  Ext_Data (i,j,k) * mass_n(i,j,k) -         &
     &                       Data_out_mono(i,j,k) * mass_np1(i,j,k)
              beta (i,j,k) = (Data_out_high(i,j,k) -                    &
     &                        Data_out_mono(i,j,k)) * mass_np1(i,j,k)
            End Do
          End Do
        End Do

! now sum into field arrays
        Do k = 2, dim_k_in
!DIR$ SUPPRESS fielda
          Do j = 1, dim_j_in
            Do i = 1, dim_i_in
              fielda(i,j,1) = fielda(i,j,1) + work(i,j,k)
              fielda(i,j,2) = fielda(i,j,2) + beta(i,j,k)
            End Do
          End Do
        End Do

! calculate global sums using rvecsumr
      call gc_ssync(n_proc,istat)
      call gcg_rvecsumr(dim_i_in,dim_i_in,1,                            &
     &   dim_j_in*2,fielda, proc_row_group,istat,sum_rows)
      call gc_ssync(n_proc,istat)
      call gcg_rvecsumr(dim_j_in,dim_j_in,1,                            &
     &     2,sum_rows, proc_col_group,istat,sum_all)

      C_star = sum_all(1)
      Integral_beta = sum_all(2)

! ----------------------------------------------------------------------
!  section 2.2  Check to see if C_star = integral of beta, unlikely.
! ----------------------------------------------------------------------

        L_conservation_enforced = .false.

        If (C_star  ==  Integral_beta ) Then
          L_conservation_enforced = .true.

!cdir collapse
          Do k = 1, dim_k_in
            Do j = 1, dim_j_in
              Do i = 1, dim_i_in
                alpha_cons(i,j,k) = alpha_max(i,j,k)
                L_alpha_set(i,j,k) = .true.
              End Do
            End Do
          End Do

        Else If ( Integral_beta  <   C_star) Then

          C_star = - C_star

          beta(:,:,:) = - beta(:,:,:)

        End If

! ----------------------------------------------------------------------
!  section 2.3  set initial values of alpha if beta less than or equal
!               to zero, step 1 in Priestley's algorithm.
! ----------------------------------------------------------------------

        If (.not. L_conservation_enforced ) Then
!cdir collapse
          Do k = 1, dim_k_in
            Do j = 1, dim_j_in
              Do i = 1, dim_i_in
                If ( beta(i,j,k)  <=  0.) Then
                  alpha_cons(i,j,k) = alpha_max(i,j,k)
                  L_alpha_set(i,j,k) = .true.
                Else
                  alpha_cons(i,j,k) = 0.
                  L_alpha_set(i,j,k) = .false.
                End If
              End Do
            End Do
          End Do
        End If

! ----------------------------------------------------------------------
!  section 2.4  Loop over steps 2 to 5 in Priestley's algorithm.
! ----------------------------------------------------------------------

        Do while ( .not. L_conservation_enforced )

! calculate surplus and sum beta at points where alpha has not been
! set.

          k = 1
!cdir collapse
            Do j = 1, dim_j_in
              Do i = 1, dim_i_in
                If (L_alpha_set(i,j,k)) Then
                  fielda(i,j,1) = alpha_cons(i,j,k) * beta(i,j,k)
                  fielda(i,j,2) = 0.0
                Else
                  fielda(i,j,1) = 0.0
                  fielda(i,j,2) = beta(i,j,k)
                End If
              End Do
            End Do

          Do k = 2, dim_k_in
            Do j = 1, dim_j_in
              Do i = 1, dim_i_in
                If (L_alpha_set(i,j,k)) Then
                  fielda(i,j,1) = fielda(i,j,1) +                       &
     &                          alpha_cons(i,j,k) * beta(i,j,k)
                Else
                  fielda(i,j,2) = fielda(i,j,2) + beta(i,j,k)
                End If
              End Do
            End Do
          End Do

! Calculate global sum of surplus and integral_beta
! calculate global sums using rvecsumr

          call gc_ssync(n_proc,istat)
          call gcg_rvecsumr(dim_i_in,dim_i_in,1,                        &
     &          dim_j_in*2,fielda, proc_row_group,istat,sum_rows)


          call gc_ssync(n_proc,istat)

          call gcg_rvecsumr(dim_j_in,dim_j_in,1,                        &
     &          2,sum_rows, proc_col_group,istat,sum_all)





          surplus = C_star - sum_all(1)
          Integral_beta = sum_all(2)

! Check surplus is not negative and integral of beta is non-zero.
! If it is print warning message and
! terminate code returning montone solution.

          average_alpha = 0.0
          If (surplus  <   0. .or. Integral_beta  ==  0.) Then

! this terminates code
            L_conservation_enforced = .true.
            conserv_fail = .true.

! this causes code to use monotone alphas

!cdir collapse
            Do k = 1, dim_k_in
              Do j = 1, dim_j_in
                Do i = 1, dim_i_in
                  L_alpha_set(i,j,k) = .true.
                  alpha_cons(i,j,k) = alpha_max(i,j,k)
                End Do
              End Do
            End Do

! warning message

            If (me  ==  0) Then
            Print*, "****************** WARNING *******************"
            Print*, "** Conservation enforcement failed          **"
            Print*, "** Run continuing using best estimate       **"
            Print*, "****************** WARNING *******************"
            End If

          Else

! calculate average value of alpha.

            average_alpha = surplus / Integral_beta

! set points which have not been set.
! set values only where average alpha exceeds alpha max.

            count = 0

!cdir collapse
            Do k = 1, dim_k_in
              Do j = 1, dim_j_in
                Do i = 1, dim_i_in
                  If (.not. L_alpha_set(i,j,k) .and.                    &
     &                average_alpha  >   alpha_max(i,j,k) ) Then
                    alpha_cons (i,j,k) = alpha_max (i,j,k)
                    L_alpha_set(i,j,k) = .true.
                    count = count + 1
                  End If
                End Do
              End Do
            End Do

            call gc_isum(1, n_proc, istat, count)

! If average alpha is less than alpha max everywhere, then set all
! values to this average values.

            If (count  ==  0 ) Then
              L_conservation_enforced = .true.
            End If

! End if for negative surplus
          End If

! End do while loop
        End Do

! If average alpha is less than alpha max everywhere, then set all
! unset values to this average values.
!cdir collapse
        Do k = 1, dim_k_in
          Do j = 1, dim_j_in
            Do i = 1, dim_i_in
              If (.not. L_alpha_set(i,j,k) ) Then
                alpha_cons (i,j,k) = average_alpha
              End If
            End Do
          End Do
        End Do

! ----------------------------------------------------------------------
!  section 2.5  Set output data.
! ----------------------------------------------------------------------

        Do k = 1, dim_k_out
          Do j = 1, dim_j_out
            Do i = 1, dim_i_out
               Data_out (i,j,k) = (1.- alpha_cons(i,j,k)) *             &
     &                             Data_out_mono(i,j,k)                 &
     &                            + alpha_cons(i,j,k) *                 &
     &                             Data_out_high(i,j,k)
            End Do
          End Do
        End Do

      IF(conserv_fail  .and. (model_domain==mt_global) )then
        temp1(:,:)=0.0
        temp2(:,:)=0.0
        Do k = 1, dim_k_out
          Do j = 1, dim_j_out
            Do i = 1, dim_i_out
              temp1(i,j) = temp1(i,j)+ (Ext_Data (i,j,k)                &
     &           * mass_n(i,j,k))
              temp2(i,j) = temp2(i,j)+ (data_out (i,j,k)                &
     &           * mass_np1(i,j,k))
            End Do
          End Do
        End Do

        call gcg_rvecsumr(dim_i_out,dim_i_out,1,                        &
     &                dim_j_out,temp1,proc_row_group,istat,temp1_rows)
        call gcg_rvecsumr(dim_j_out,dim_j_out,1,                        &
     &                    1,temp1_rows,proc_col_group,istat,temp1_sum)
        call gcg_rvecsumr(dim_i_out,dim_i_out,1,                        &
     &                dim_j_out,temp2,proc_row_group,istat,temp2_rows)
        call gcg_rvecsumr(dim_j_out,dim_j_out,1,                        &
     &                    1,temp2_rows,proc_col_group,istat,temp2_sum)

!     fix to make mass conserve even if conservation fails
        if( temp2_sum(1) /= 0.0 ) then
          temp3_sum(1)=temp1_sum(1)/temp2_sum(1)
          data_out(:,:,:)=data_out(:,:,:)*temp3_sum(1)
        endif

      endif        ! conserv_fail

      return
      END SUBROUTINE Mono_Conserv_qcon

