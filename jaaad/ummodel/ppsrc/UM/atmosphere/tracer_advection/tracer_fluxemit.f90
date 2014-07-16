SUBROUTINE TRACER_FLUXEMIT(       &
   &        row_length, rows, mype, timestep, &
   &        flux, off_x, off_y, &
   &        FV_cos_theta_latitude, delta_lambda, delta_phi, &
   &        sumemit)

    IMPLICIT NONE
! constants required for calculation
!*L------------------COMDECK C_PI---------------------------------------
!LL
!LL 4.0 19/09/95  New value for PI. Old value incorrect
!LL               from 12th decimal place. D. Robinson
!LL 5.1 7/03/00   Fixed/Free format P.Selwood
!LL

      ! Pi
      Real, Parameter :: Pi                 = 3.14159265358979323846

      ! Conversion factor degrees to radians
      Real, Parameter :: Pi_Over_180        = Pi/180.0

      ! Conversion factor radians to degrees
      Real, Parameter :: Recip_Pi_Over_180  = 180.0/Pi

!*----------------------------------------------------------------------
!*L------------------COMDECK C_A----------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Replace variable A by more meaningful name for
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Convert to Fixed/Free format. P. Selwood

      ! Mean radius of Earth in m.
      Real, Parameter  :: Earth_Radius = 6371000.

!*----------------------------------------------------------------------

!---------------------------------
! EXTERNAL SUBROUTINE CALLS - 
!---------------------------------
     external do_sums
!---------------------------------

!---------------------------------
!INPUT variables
!--------------------------------

    INTEGER row_length, rows, mype, off_x, off_y
    REAL timestep
    REAL          &
     flux(row_length,rows)

! Input trigonometric functions (required by polar_wind)
      Real                                                              &
     & FV_cos_theta_latitude(1-off_x:row_length+off_x,                  &
     &                        1-off_y:rows+off_y)

! Input model grid spacing in radians

      REAL delta_lambda,delta_phi

!----------------------------
!LOCAL variables
!----------------------------

    integer, parameter :: n_sums=1
    real :: combo_results(row_length,rows,n_sums),sum_results(n_sums)
    real :: sumemit
    integer :: i,j,k

!-----------------------------
!Zero out summation arrays
!-----------------------------
do k=1,n_sums
   sum_results(k) = 0.0
enddo
do k=1,n_sums
   do j=1,rows
      do i=1,row_length
         combo_results(i,j,k) = 0.0
      enddo
   enddo
enddo

!----------------------------
!Sum over all latitudes and longitudes
   do i=1,row_length
     do j=1,rows  
        combo_results(i,j,1) = flux(i,j)*FV_cos_theta_latitude(i,j)
     enddo
   enddo

!---------------------------
!Combine values from all processors
! DEPENDS ON: do_sums
   call do_sums(combo_results,row_length,rows,0,0,1,3,n_sums,sum_results)

   sumemit = sum_results(1)*timestep*earth_radius*earth_radius  &
                 *delta_phi*delta_lambda

return

END SUBROUTINE TRACER_FLUXEMIT



