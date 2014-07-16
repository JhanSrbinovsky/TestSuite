
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!     subroutine P_TO_U  for calculating variables held at p points
!     at u points, on the new dynamics grid
!
!     This routine does interior points of array not halos,
!     but requires halo information to be set.
!
!     global code has E-W wrap around
!     LAM code has most east u point set to zero
! History:
! Version   Date     Comment
! -------   ----     -------
! 5.3  17/10/01 Changes required for Single Column Model
!                                             Z. Gardner
! 6.2  03/02/06 Move to c92_2a.  P.Selwood
!
      Subroutine P_TO_U(array_on_p_points,row_length,rows,levels,       &
     &                  halo_i,halo_j,                                  &
     &                  array_on_u_points)
!
      Integer                                                           &
               ! INTENT (IN)
     & row_length                                                       &
     &,rows                                                             &
     &,levels                                                           &
     &,halo_i                                                           &
     &,halo_j

      Real                                                              &
               ! INTENT (IN)
     & array_on_p_points(1-halo_i:row_length+halo_i,                    &
     &                   1-halo_j:rows+halo_j,levels)

      Real                                                              &
               ! INTENT (OUT)
     & array_on_u_points(row_length,rows,levels)

      Integer                                                           &
               ! local variables
     & i,j,k

      do k=1,levels
        do j=1,rows

          do i=1,row_length

            array_on_u_points(i,j,k)= 0.5*                              &
     &    ( array_on_p_points(i,j,k) + array_on_p_points(i+1,j,k) )




          end do

        end do
      end do

      return
      END SUBROUTINE P_TO_U
