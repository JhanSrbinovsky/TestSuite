#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculates the heights in the SCM.
!
! Subroutine Interface:
      SUBROUTINE Calc_levels                                            &
! Input data
     &    (eta_theta_levels, eta_rho_levels, z_top_of_model, orog,      &
     &     height_gen_method, boundary_layer_levels, model_levels,      &
     &     first_constant_r_rho_level, rows, row_length,                &
! Output data
     &     r_theta_levels, r_rho_levels)

      IMPLICIT NONE

!
! Description:  To set up r_theta_levels and r_rho_levels for the
!               Single Column Model
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

!     Inputs

      INTEGER                                                           &
     &  boundary_layer_levels                                           &
     &, model_levels                                                    &
     &, first_constant_r_rho_level                                      &
     &, rows                                                            &
     &, row_length                                                      &
     &, height_gen_method

      REAL                                                              &
     &  eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)                                    &
     &, z_top_of_model                                                  &
     &, orog(row_length, rows)

!     Outputs

      REAL                                                              &
     &  r_theta_levels(row_length, rows, 0:model_levels)                &
     &, r_rho_levels(row_length, rows, model_levels)

      CHARACTER*(80)                                                    &
     &       CMESSAGE              ! Error message if ICODE >0
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='CALCLEVS')

      INTEGER :: ErrorStatus

!     Local variables

      REAL                                                              &
     &  r_ref_theta(model_levels)                                       &
     &, r_ref_rho(model_levels)

      INTEGER                                                           &
     & i,j,k

      INTEGER, PARAMETER :: height_gen_original = 1 ! methods for height
      INTEGER, PARAMETER :: height_gen_smooth   = 2 ! generation

! External subroutines called:
      EXTERNAL                                                          &
     &  Ereport

!----------------------------------------------------------------------
!     Set up heights
!----------------------------------------------------------------------

! Set reference profile

      Do k = 1, model_levels
        r_ref_theta(k) = eta_theta_levels(k) * z_top_of_model
        r_ref_rho(k) = eta_rho_levels(k) * z_top_of_model
      End Do

! Set bottom level, ie orography
      Do j = 1, rows
        Do i = 1, row_length
          r_theta_levels(i,j,0) = orog(i,j) + Earth_radius
        End Do
      End Do
! For constant levels set r to be a constant on the level
      Do k = first_constant_r_rho_level, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            r_theta_levels(i,j,k) = Earth_radius + r_ref_theta(k)
            r_rho_levels(i,j,k) = Earth_radius + r_ref_rho(k)
          End Do
        End Do
      End Do

      Select Case( height_gen_method)
        Case( height_gen_original)
! The original version of height generation used in the SI dynamics
!
! For boundary layer levels set depth to be constant.
          Do k = 1, boundary_layer_levels
            Do j = 1, rows
              Do i = 1, row_length
                r_theta_levels(i,j,k) = r_theta_levels(i,j,0) +         &
     &                                     r_ref_theta(k)
                r_rho_levels(i,j,k) = r_theta_levels(i,j,0) +           &
     &                                     r_ref_rho(k)
              End Do
            End Do
          End Do
! For intemediate levels use linear relaxation to constant value.
! set orographic heights.
          Do k = boundary_layer_levels+1,                               &
     &                   first_constant_r_rho_level-1
            Do j = 1, rows
              Do i = 1, row_length
                r_rho_levels(i,j,k) =                                   &
     &            ( r_rho_levels(i,j,first_constant_r_rho_level) -      &
     &              r_theta_levels(i,j,boundary_layer_levels) ) *       &
     &            ( eta_rho_levels(k) -                                 &
     &              eta_theta_levels(boundary_layer_levels) )/          &
     &            (eta_rho_levels(first_constant_r_rho_level) -         &
     &             eta_theta_levels(boundary_layer_levels) )            &
     &            +  r_theta_levels(i,j,boundary_layer_levels)
                r_theta_levels(i,j,k) =                                 &
     &            ( r_rho_levels(i,j,first_constant_r_rho_level) -      &
     &              r_theta_levels(i,j,boundary_layer_levels) ) *       &
     &            ( eta_theta_levels(k) -                               &
     &              eta_theta_levels(boundary_layer_levels) ) /         &
     &            ( eta_rho_levels(first_constant_r_rho_level) -        &
     &              eta_theta_levels(boundary_layer_levels) )           &
     &            +  r_theta_levels(i,j,boundary_layer_levels)
              End Do
            End Do
          End Do

        Case( height_gen_smooth )
! A smooth quadratic height generation
          Do k = 1, first_constant_r_rho_level-1
            Do j = 1, rows
              Do i= 1, row_length
              r_rho_levels(i,j,k) = eta_rho_levels(k) * z_top_of_model +&
     &         Earth_radius + Orog(i,j) * (1.0 - eta_rho_levels(k)      &
     &              /eta_rho_levels(first_constant_r_rho_level))**2
              r_theta_levels(i,j,k) = eta_theta_levels(k) *             &
     &             z_top_of_model + Earth_radius + Orog(i,j) *          &
     &             (1.0 - eta_theta_levels(k) /                         &
     &              eta_rho_levels(first_constant_r_rho_level))**2
              End Do
            End Do
          End Do

        Case Default
          ErrorStatus = 10
          Write (Cmessage,*) 'Unrecognised height generation method - ',&
     &                       'Dump needs to be reconfigured'
! DEPENDS ON: ereport
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
      End Select

      Return
      END SUBROUTINE Calc_levels

#endif
