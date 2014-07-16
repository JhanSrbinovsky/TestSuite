#if defined(A12_2A) || defined(MAKEBC)
! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Calc_Exner_at_theta

      Subroutine Calc_Exner_at_theta(                                   &
     &                      r_theta_levels, r_rho_levels, exner_rho,    &
     &                      row_length, rows, model_levels,             &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      exner_theta_levels, L_include_halos)

! Purpose:
!          Calculates exner pressure at theta levels.
!
! Method:
!          Linear interpolation in height of exner (at rho levels)
!
!
!
! Original Programmer: Clive Wilson
! Current code owner: Andrew J. Malcolm
!
! History:
! Date     Version     Comment
! ----     -------     -------
! 27/07/01  5.3        Original code                        C Wilson
! 02/09/04  6.1        Added MAKEBC def                     R Sempers
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
     &, model_levels                                                    &
                         ! number of model levels.
     &, halo_i                                                          &
                   ! Size of halo in i direction.
     &, halo_j                                                          &
                   ! Size of halo in j direction.
     &, off_x                                                           &
                   ! Size of small halo in i
     &, off_y      ! Size of small halo in j.
      LOGICAL                                                           &
     &  L_include_halos  ! If .TRUE. then include halo regions
                         ! when performing the calculations

      Real                                                              &
     &  exner_rho(1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &      model_levels + 1)                                           &
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                   1-halo_j:rows+halo_j,0:model_levels)           &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &                 1-halo_j:rows+halo_j, model_levels)

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
     &  exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                 1-off_y:rows+off_y, model_levels)

! Local Variables.

      Integer                                                           &
     &  i, j, k                                                         &
                     ! Loop indices
     &,  i_start,i_end                                                  &
                         ! Loop bounds
     &,  j_start,j_end   ! Loop bounds

! No External Routines:

! ----------------------------------------------------------------------
! Section 1.   Calculate pressure at desired theta levels.
! ----------------------------------------------------------------------

      IF (L_include_halos) THEN
        i_start=1-off_x
        i_end=row_length+off_x
        j_start=1-Off_y
        j_end=rows+Off_y
      ELSE
        i_start=1
        i_end=row_length
        j_start=1
        j_end=rows
      ENDIF

      Do k = 1, model_levels - 1
        Do j = j_start, j_end
          Do i = i_start, i_end
            exner_theta_levels(i,j,k) = ( exner_rho(i,j,k) *            &
     &                                (r_rho_levels(i,j,k+1) -          &
     &                                 r_theta_levels(i,j,k) ) +        &
     &                                exner_rho(i,j,k+1) *              &
     &                                (r_theta_levels(i,j,k) -          &
     &                                 r_rho_levels(i,j,k) ) ) /        &
     &                                (r_rho_levels(i,j,k+1) -          &
     &                                 r_rho_levels(i,j,k) )
          End Do
        End Do
      End Do

        k = model_levels

!AM  extra pressure level above top theta level is same height above
!AM  as is the pressure below - hence weights are 0.5 and there is
!AM  no need to store the r_rho_level for the extra pressure
      Do j = j_start, j_end
        Do i = i_start, i_end
            exner_theta_levels(i,j,k) = 0.5 *                           &
     &                      ( exner_rho(i,j,k) + exner_rho(i,j,k+1) )
        End Do
      End Do

! End of routine
      return
      END SUBROUTINE Calc_Exner_at_theta
!
! Subroutine Calc_P_from_Exner
!


#endif
