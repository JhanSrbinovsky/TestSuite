#if defined(A32_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to calculate ThetaL from Theta

      Subroutine Theta_to_ThetaL (                                      &
     &           Theta, QCF, QCL, Exner_P_Theta, ThetaL,                &
     &           row_length, rows, off_x, off_y, halo_x, halo_y,        &
     &           wet_levels, model_levels)

! Description:
!   This routine derives ThetaL from the Theta/QCL/QCF/Exner_P.
!   ThetaL is a prognostic for the UM prior to the New Dynamics.
!
! Method:
!   Formula used :
!   ThetaL = Theta - (LC*QCL + (LC+LF)*QCF) / (CP * Exner_P_Theta)
!
! Owner: Dave Robinson
!
! History:
! Version Date     Comment
! ------- ----     -------
! 6.1     18/08/04 Original code. Dave Robinson
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

      Implicit None

      Integer :: row_length
      Integer :: rows
      Integer :: off_x, off_y
      Integer :: halo_x, halo_y
      Integer :: wet_levels
      Integer :: model_levels
      Integer :: theta_field_size

      Real, Intent(IN) :: Theta (1-off_x:row_length+off_x,              &
     &                           1-off_y:rows+off_y,                    &
     &                           model_levels)
      Real, Intent(IN) :: QCF (1-halo_x:row_length+halo_x,              &
     &                         1-halo_y:rows+halo_y,                    &
     &                         wet_levels)
      Real, Intent(IN) :: QCL (1-halo_x:row_length+halo_x,              &
     &                         1-halo_y:rows+halo_y,                    &
     &                         wet_levels)
      Real, Intent(IN) :: Exner_P_Theta(1-off_x:row_length+off_x,       &
     &                                  1-off_y:rows+off_y,             &
     &                                  model_levels)
      Real, Intent(OUT) :: ThetaL(1-off_x:row_length+off_x,             &
     &                            1-off_y:rows+off_y,                   &
     &                            model_levels)

#include "c_lheat.h"
#include "c_r_cp.h"

      Integer :: i, j, k

      Do k = 1, wet_levels
        Do j = 1-off_y, rows+off_y
          Do i = 1-off_x, row_length+off_x
            ThetaL(i,j,k) =                                             &
     &      Theta (i,j,k) - ( LC*QCL(i,j,k) + (LC+LF)*QCF(i,j,k) )      &
     &                      / ( CP*Exner_P_Theta(i,j,k) )
          End Do
        End Do
      End Do

      If (wet_levels < model_levels) Then

        Do k = wet_levels+1, model_levels
          thetaL(:,:,k) = theta(:,:,k)
        End Do

      End If

      Return
      END SUBROUTINE Theta_to_ThetaL
#endif
