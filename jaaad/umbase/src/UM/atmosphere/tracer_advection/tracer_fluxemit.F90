SUBROUTINE TRACER_FLUXEMIT(       &
   &        row_length, rows, mype, timestep, &
   &        flux, off_x, off_y, &
   &        FV_cos_theta_latitude, delta_lambda, delta_phi, &
   &        sumemit)

    IMPLICIT NONE
! constants required for calculation
#include "c_pi.h"
#include "c_a.h"

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



