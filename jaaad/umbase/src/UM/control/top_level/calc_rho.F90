#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculates rho in the SCM.
!
! Subroutine Calc_rho
      Subroutine calc_rho( model_levels, rows,                          &
     &                     row_length, exner_rho_levels,                &
     &                     exner_theta_levels, r_rho_levels,            &
     &                     r_theta_levels, theta, p,                    &
! ARGUMENTS OUT
     &                     rho)


      IMPLICIT NONE

!
! Description:  To calculate rho for the Single Column Model
!
!
! Method:
!   It follows the code used in the main UM.
!
! Owner: UM System team
!
! History:
! Version  Date     Comment
! =======  ====     =======
! 5.3      30/04/01 Original code. (Zoe Gardner)
!
! Code description:
!   FORTRAN 77 + common extensions also in fortran 90.
!   This code is written to UM programming standards version 7.4

!     INCLUDED COMDECKS
#include "c_r_cp.h"

      !     Inputs

      INTEGER                                                           &
     &  model_levels                                                    &
     &, rows                                                            &
     &, row_length

      REAL                                                              &
     &  r_theta_levels(row_length, rows, 0:model_levels)                &
     &, r_rho_levels(row_length, rows, model_levels)                    &
     &, exner_rho_levels(row_length, rows, model_levels+1)              &
     &, exner_theta_levels(row_length, rows, model_levels)              &
     &, p(row_length, rows,model_levels+1)                              &
     &, theta(row_length, rows,model_levels)  ! Temperature(K)


!     Outputs

      REAL                                                              &
     &  rho(row_length, rows, model_levels)

!     Local variables

      INTEGER                                                           &
     & i,j,k

      Real                                                              &
     &  temp1

! 3-d work arrays
      Real                                                              &
     &  weight_upper(row_length, rows, model_levels)                    &
     &, weight_lower(row_length, rows, model_levels)

!----------------------------------------------------------------------
!      Calculate r**2 rho from equation of state
!----------------------------------------------------------------------

        k = 1
        Do j = 1, rows
          Do i = 1, row_length
            rho(i,j,k) = r_rho_levels(i,j,k) * r_rho_levels(i,j,k) *    &
     &                   p(i,j,k) /                                     &
     &                   (R * theta(i,j,k) * exner_rho_levels(i,j,k))
          End Do
        End Do

! set up vertical interpolation weights between theta levels and rho
! levels
      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            weight_upper(i,j,k) = (r_rho_levels(i,j,k)                  &
     &                             - r_theta_levels(i,j,k-1) )          &
     &                           / (r_theta_levels(i,j,k)               &
     &                              - r_theta_levels(i,j,k-1) )
            weight_lower(i,j,k) = 1.0 - weight_upper(i,j,k)
          End Do
        End Do
      End Do

      Do k = 2, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            temp1 = weight_upper(i,j,k) * theta(i,j,k) +                &
     &                 weight_lower(i,j,k) * theta(i,j,k-1)
            rho(i,j,k) = r_rho_levels(i,j,k)*r_rho_levels(i,j,k) *      &
     &             p(i,j,k) /(R * temp1 * exner_rho_levels(i,j,k))
          End Do
        End Do
      End Do


      Return
      END SUBROUTINE calc_rho

#endif
