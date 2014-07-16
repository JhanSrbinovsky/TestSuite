#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculates exner in the SCM.
!
!  Subroutine Calc_press
      SUBROUTINE Calc_press                                             &
! Input data
     &    (model_levels, wet_model_levels, rows, row_length, p, theta,  &
     &     q, r_theta_levels, r_rho_levels, l_calc_exner, l_calc_rho,   &
! In/Out
     &     rho,                                                         &
! Output data
     &     exner_theta_levels, exner_rho_levels, p_theta_levels,        &
     &     rp, rp_theta, p_star)

      IMPLICIT NONE

!
! Description:  To calculate exner, pressure on theta levels,
!               the reciprocol or pressure and rho (if required)
!               for the Single Column Model
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
#include "c_a.h"
#include "c_r_cp.h"
#include "c_g.h"

!     Inputs

      INTEGER                                                           &
     &  model_levels                                                    &
     &, wet_model_levels                                                &
     &, rows                                                            &
     &, row_length

      REAL                                                              &
     &  p(row_length, rows, model_levels+1)                             &
     &, q(row_length, rows,wet_model_levels)                            &
     &, theta(row_length, rows,model_levels)  ! Temperature(K)

      LOGICAL                                                           &
     &  l_calc_exner                                                    &
                      ! If true then exner is calculated from p,
                      ! if false then p is calculated from exner
     &, l_calc_rho    ! If true then rho is calculated

!     In/Out

      REAL                                                              &
     &  rho(row_length, rows, model_levels)

!     Outputs

      REAL                                                              &
     &  r_theta_levels(row_length, rows, 0:model_levels)                &
     &, r_rho_levels(row_length, rows, model_levels)                    &
     &, exner_rho_levels(row_length, rows, model_levels+1)              &
     &, exner_theta_levels(row_length, rows, model_levels)              &
     &, p_theta_levels(row_length, rows, model_levels)                  &
     &, rp_theta(row_length, rows, model_levels)                        &
                                                 ! reciprocol pressure
     &, rp(row_length, rows, model_levels+1)                            &
                                              ! reciprocol pressure
     &, p_star(row_length,rows) ! surface pressure

      CHARACTER*(80)                                                    &
     &       CMESSAGE              ! Error message if ICODE >0
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='CALCPRESS')

      INTEGER                                                           &
     &  icode

!     Local variables

      INTEGER                                                           &
     & i,j,k

      REAL                                                              &
     &  constant

!----------------------------------------------------------------------

! 0.1 Initialise secondary arrays.
! calculate p from exner_rho_levels if l_calc_exner is false
! calculate exner_rho_levels from p if l_calc_exner is true
! [halos required for diagnostic calculations at T+0 in INITDIAG.]
      constant = 1./ kappa
      If (l_calc_exner) then
        Do k = 1, model_levels+1
          Do j = 1, rows
            Do i = 1, row_length
              exner_rho_levels(i,j,k) =                                 &
     &                 (p(i,j,k)/p_zero)**(1.0/constant)
            End Do
          End Do
        End Do
      Else
        Do k = 1, model_levels+1
          Do j = 1, rows
            Do i = 1, row_length
              p(i,j,k) = (exner_rho_levels(i,j,k) ** constant) *        &
     &                            p_zero
            End Do
          End Do
        End Do
      End If

! calculate Exner at theta levels
! DEPENDS ON: calc_exner_at_theta
      Call Calc_Exner_at_theta (r_theta_levels, r_rho_levels,           &
     &                 exner_rho_levels,                                &
     &                 row_length, rows, model_levels,                  &
     &                 0,0,0,0,exner_theta_levels, .FALSE.)

! calculate p at theta_levels
! DEPENDS ON: calc_p_from_exner
      Call Calc_P_from_exner (p_theta_levels, kappa, p_zero,            &
     &                 row_length, rows, model_levels,                  &
     &                 0,0,                                             &
     &                 exner_theta_levels,.FALSE.)

! Calculate rho if required (ie if this is the start of the run).

      If (L_calc_rho) Then
! DEPENDS ON: calc_rho
        Call Calc_rho( model_levels, rows,                              &
     &                 row_length, exner_rho_levels,                    &
     &                 exner_theta_levels,                              &
     &                 r_rho_levels,                                    &
     &                 r_theta_levels, theta, p,                        &
     &                 rho)
      End If

! calculate p_star using rho and p on model levels
! DEPENDS ON: calc_p_star
      Call Calc_P_star (r_theta_levels, r_rho_levels, p, rho,           &
     &                  g, row_length, rows, model_levels,              &
     &                  0,0,0,0,                                        &
     &                  p_star)

! calculate reciprocol pressures

      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            rp(i,j,k) = 1.0/p(i,j,k)
            rp_theta(i,j,k) = 1.0/p_theta_levels(i,j,k)
          End Do
        End Do
      End Do

      Do j = 1, rows
        Do i = 1, row_length
          rp(i,j,model_levels+1) = 1.0/p(i,j,model_levels+1)
        End Do
      End Do

      Return
      END SUBROUTINE Calc_press

#endif
