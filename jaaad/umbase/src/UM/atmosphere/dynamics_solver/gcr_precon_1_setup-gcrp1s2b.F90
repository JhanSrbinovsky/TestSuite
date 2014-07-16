#if defined(A10_2B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine GCR_precon_1_setup

      Subroutine GCR_precon_1_setup(                                    &
     &                     HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,          &
     &                     HM_Czz, HM_Cz, HM_C3, HM_C4,                 &
     &                     eta_theta_levels, eta_rho_levels,            &
     &                     row_length, rows, n_rows, model_levels,      &
     &                     FV_sec_theta_latitude,                       &
     &                     dlambda_p, dphi_v,                           &
     &                     recip_dlamu, recip_dphiv,                    &
     &                     wt_lambda_u, wt_phi_v,                       &
     &                     model_domain,                                &
     &                     weight_upper, weight_lower,                  &
     &                     a0, a1, factor, Pi,                          &
     &                     offx, offy, halo_i, halo_j,                  &
     &                     global_row_length, proc_row_group,           &
     &                     at_extremity, L_regular,                     &
     &                     i_start, i_stop, j_start, j_stop,            &
     &                     j_begin, j_end                               &
     &                     )

! Purpose:
!          Calculates constant terms used inpre-conditioning operator
!          applied to field.
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
! Version   Date       Comment
! ----     -------     -------
!   6.2   25/12/05  Variable resolution changes    Yongming Tang
!   6.2   25/12/05  Pass loop bounds variables as arguments
!                                                     Terry Davies
!   6.2   21/10/05  Remove unused external.           P.Selwood.
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
     &, model_domain

      Real                                                              &
     &  Pi
      Integer                                                           &
     &  offx, offy                                                      &
     &, halo_i, halo_j                                                  &
     &, global_row_length                                               &
     &, i_start, i_stop                                                 &
                                     ! loop bounds set in PE_Helmholtz
     &, j_start, j_stop                                                 &
                                     ! loop bounds set in PE_Helmholtz
     &, j_begin, j_end               ! loop bounds set in PE_Helmholtz

! interpolation weights for moving between theta levels and rho levels
      Real                                                              &
     &  weight_upper(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)                                          &
     &, weight_lower(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)

      Real                                                              &
     &  HM_Cxx1(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cxx2(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cyy1(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cyy2(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Czz(1-offx:row_length+offx,1-offy:rows+offy,model_levels)    &
     &, HM_Cz(1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, HM_C3(1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, HM_C4(1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, FV_sec_theta_latitude(1-offx:row_length+offx,1-offy:rows+offy)

      Real                                                              &
           ! vertical co-ordinate information
     &  eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)

       Real                                                             &
                   !  VarRes horizontal co-ordinate spacing.
     &  dlambda_p(1-halo_i:row_length+halo_i)                           &
     &, dphi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)      &
     &, recip_dlamu(1-halo_i : row_length + halo_i)                     &
     &, recip_dphiv(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j) &
     &, wt_lambda_u(1-halo_i:row_length+halo_i)                         &
     &, wt_phi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)

! Arguments with Intent OUT. ie: Output variables.
      Real                                                              &
     &  a0(row_length,rows,model_levels)                                &
     &, a1(row_length,rows,model_levels)                                &
     &, factor(row_length,rows,model_levels)

! Local Variables.

      Integer                                                           &
     &  i,j,k                                                           &
     &, istat

      Real                                                              &
     &  factor_1                                                        &
     &, factor_2                                                        &
     &, factor_3

!  parallel variables
      Integer                                                           &
     &  proc_row_group

      Real                                                              &
     &  sum_n(model_levels)                                             &
     &, sum_s(model_levels)                                             &
     &, sum_n_component(row_length,model_levels)                        &
     &, sum_s_component(row_length,model_levels)

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_regular


! Local arrays.

      Real                                                              &
     &  a2(row_length,rows,model_levels)                                &
     &, work1(1-offx:row_length+offx, 1-offy:rows+offy)                 &
     &, work2(1-offx:row_length+offx, 1-offy:rows+offy)


#include "parparm.h"
#include "domtyp.h"
!   External routines.
      External                                                          &
     &   gcg_rvecsumr

!-----------------------------------------------------------------------
!     Section 1.
!-----------------------------------------------------------------------

      Do k = 1, model_levels

        If ( k  ==  1) Then
          factor_1 = 1. / (eta_theta_levels(k) -                        &
     &                     eta_theta_levels(k-1))
          factor_2 = 1. / (eta_rho_levels(k+1) -                        &
     &                     eta_rho_levels(k))
          Do j = j_start, j_stop
            Do i = i_start, i_stop
              a1(i,j,k) = factor_2 * (factor_1 * HM_Czz(i,j,k) +        &
     &                    HM_C3(i,j,k) * HM_Cz(i,j,k) )
              a0(i,j,k) = - a1(i,j,k) - HM_C4(i,j,k)
            End Do
          End Do

        Else if (k  ==  model_levels) Then
          factor_1 = 1. / (eta_theta_levels(k) -                        &
     &                     eta_theta_levels(k-1))
          factor_3 = 1. / (eta_rho_levels(k) -                          &
     &                     eta_rho_levels(k-1))

          Do j = j_start, j_stop
            Do i = i_start, i_stop
              a2(i,j,k) = factor_3 * (factor_1 * HM_Czz(i,j,k-1) -      &
     &                    HM_C3(i,j,k) * HM_Cz(i,j,k-1) *               &
     &                                weight_lower(i,j,k) )
              a0(i,j,k) = - a2(i,j,k) - HM_C4(i,j,k)
            End Do
          End Do

        Else  ! model_levels > k > 1

          factor_1 = 1. / (eta_theta_levels(k) -                        &
     &                     eta_theta_levels(k-1))
          factor_2 = 1. / (eta_rho_levels(k+1) -                        &
     &                     eta_rho_levels(k))
          factor_3 = 1. / (eta_rho_levels(k) -                          &
     &                     eta_rho_levels(k-1))

          Do j = j_start, j_stop
            Do i = i_start, i_stop
              a1(i,j,k) = factor_2 * (factor_1 * HM_Czz(i,j,k) +        &
     &                    HM_C3(i,j,k) * HM_Cz(i,j,k) *                 &
     &                                weight_upper(i,j,k) )
              a2(i,j,k) = factor_3 * (factor_1 * HM_Czz(i,j,k-1) -      &
     &                    HM_C3(i,j,k) * HM_Cz(i,j,k-1) *               &
     &                                weight_lower(i,j,k) )
              a0(i,j,k) = - a2(i,j,k) - a1(i,j,k)  - HM_C4(i,j,k)
            End Do
          End Do

        End If  !  k  ==  1

      End Do ! k = 1, model_levels

! Include horizontal second derivative terms in matrix.

! Interior points - Common to all model_domains.

      if ( L_regular ) then

        Do k = 1, model_levels
          Do j = j_begin, j_end
            Do i = i_start, i_stop
              a0(i,j,k) = a0(i,j,k) - FV_sec_theta_latitude(i,j) *      &
     &                        (HM_Cxx1(i-1,j,k)*HM_Cxx2(i-1,j,k) +      &
     &                             HM_Cxx1(i,j,k)*HM_Cxx2(i,j,k) +      &
     &                             HM_Cyy1(i,j,k)*HM_Cyy2(i,j,k) +      &
     &                         HM_Cyy1(i,j-1,k)*HM_Cyy2(i,j-1,k) )
            End Do
          End Do
        End Do ! k = 1, model_levels

        else  !  variable resolution

        Do k = 1, model_levels

          Do j = j_begin, j_end
            Do i = i_start-1, i_stop
              work1(i,j) = HM_Cxx1(i,j,k) * HM_Cxx2(i,j,k)
            End Do
          End Do
          Do j = j_begin-1, j_end
            Do i = i_start, i_stop
              work2(i,j) = HM_Cyy1(i,j,k) * HM_Cyy2(i,j,k)
            End Do
          End Do
          Do j = j_begin, j_end
            Do i = i_start, i_stop
              a0(i,j,k) = a0(i,j,k) - 2.0 * FV_sec_theta_latitude(i,j) *&
     &                           ( recip_dlamu(i-1) * recip_dlamu(i-1) *&
     &                                 ( work1(i-1,j) * wt_lambda_u(i) +&
     &                                     work1(i,j) *                 &
     &                                        (1.0 - wt_lambda_u(i)) ) +&
     &                         recip_dphiv(i,j-1) * recip_dphiv(i,j-1) *&
     &                                 ( work2(i,j) *                   &
     &                                           (1.0 - wt_phi_v(i,j)) +&
     &                                 work2(i,j-1) * wt_phi_v(i,j) ) )
            End Do
          End Do

      End Do ! k = 1, model_levels

      endif ! L_regular

! Treatment of North/South boundary.

      If (Model_domain  ==  mt_Global )Then

        If (at_extremity(PSouth) )then
          if ( L_regular ) then
            Do k = 1, model_levels
              Do i = 1,row_length
                sum_s_component(i,k) = - HM_Cyy1(i,1,k) * HM_Cyy2(i,1,k)
              End Do
            End Do  ! k = 1, model_levels
          else  !  variable resolution
            Do k = 1, model_levels
              Do i = 1,row_length
              sum_s_component(i,k) = - 2.0 * (1.0 - wt_phi_v(1,2)) *    &
     &                                 HM_Cyy1(i,1,k) * HM_Cyy2(i,1,k)  &
     &                                                * dlambda_p(i)
             End Do
            End Do  ! k = 1, model_levels
          endif ! L_regular
#if defined(REPROD)
          Call gcg_rvecsumr(row_length,row_length,1,                    &
     &                      model_levels,sum_s_component,               &
     &                      proc_row_group,istat,sum_s)
#else
          Call gcg_rvecsumf(row_length,row_length,1,                    &
     &                      model_levels,sum_s_component,               &
     &                      proc_row_group,istat,sum_s)
#endif
          if ( L_regular ) then
            Do k = 1, model_levels
              sum_s(k) = sum_s(k) * FV_sec_theta_latitude(1,1) /        &
     &                              global_row_length
              Do i = 1,row_length
                a0(i,1,k) = a0(i,1,k) + sum_s(k)
              End Do
            End Do  ! k = 1, model_levels
          else  !  variable resolution
            Do k = 1, model_levels
               sum_s(k) = sum_s(k) * FV_sec_theta_latitude(1,1) /       &
     &                  ( 2.0 * Pi * dphi_v(1,1) * dphi_v(1,1) )
              Do i = 1,row_length
                a0(i,1,k) = a0(i,1,k) + sum_s(k)
              End Do
            End Do  ! k = 1, model_levels
            endif ! L_regular
        End If  !  at_extremity(PSouth)

        If (at_extremity(PNorth)) then
          if ( L_regular ) then
            Do k = 1, model_levels
              Do i = 1,row_length
                sum_n_component(i,k) = - HM_Cyy1(i,rows-1,k) *          &
     &                                   HM_Cyy2(i,rows-1,k)
              End Do
            End Do  ! k = 1, model_levels
          else  !  variable resolution
            Do k = 1, model_levels
              Do i = 1,row_length
                sum_n_component(i,k) = - 2.0 * wt_phi_v(1,n_rows) *     &
     &                        HM_Cyy1(i,rows-1,k) * HM_Cyy2(i,rows-1,k) &
     &                                            * dlambda_p(i)
              End Do
            End Do  ! k = 1, model_levels
          endif ! L_regular
#if defined(REPROD)
          Call gcg_rvecsumr(row_length,row_length,1,                    &
     &                      model_levels,sum_n_component,               &
     &                      proc_row_group,istat,sum_n)
#else
          Call gcg_rvecsumf(row_length,row_length,1,                    &
     &                      model_levels,sum_n_component,               &
     &                      proc_row_group,istat,sum_n)
#endif
          if ( L_regular ) then
            Do k = 1, model_levels
              sum_n(k) = sum_n(k) * FV_sec_theta_latitude(1,rows) /     &
     &                              global_row_length
              Do i = 1,row_length
                a0(i,rows,k) = a0(i,rows,k) + sum_n(k)
              End Do
            End Do !  k = 1, model_levels
          else  !  variable resolution
            Do k = 1, model_levels
              sum_n(k) = sum_n(k) * FV_sec_theta_latitude(1,rows) /     &
     &            (2.0 * Pi * dphi_v(1,n_rows-1) * dphi_v(1,n_rows-1))
              Do i = 1,row_length
                a0(i,rows,k) = a0(i,rows,k) + sum_n(k)
              End Do
            End Do   !  k = 1, model_levels
          endif ! L_regular
        End If  !  at_extremity(PNorth)

      End If  !  Model_domain  ==  mt_Global

! Solve tridiagonal system.
! solution is in Soln
! reduce matrix to upper diagonal form.

! for global parallel code a0 only valid on interior points , NOT halos
      Do j = j_start, j_stop
        Do i = i_start, i_stop
          a0(i,j,1) = 1./a0(i,j,1)
        End Do
      End Do

      Do k= 2, model_levels
        Do j = j_start, j_stop
          Do i = i_start, i_stop
            factor(i,j,k) = a2(i,j,k) * a0(i,j,k-1)
            a0(i,j,k) = 1./(a0(i,j,k) - factor(i,j,k)*a1(i,j,k-1))
          End Do
        End Do
      End Do !  k= 2, model_levels

!     end of routine GCR_precon_1_setup

      return
      END SUBROUTINE GCR_precon_1_setup

#endif
