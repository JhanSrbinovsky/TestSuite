#if defined(C92_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!     subroutine P_TO_V  for calculating variables held at p points
!     at v points, on the new dynamics grid
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
! 6.2  03/02/06 Move to C92_2A.  P.Selwood.
!
      Subroutine P_TO_V(array_on_p_points,row_length,rows,nrows,        &
     &                  levels,halo_i,halo_j,                           &
     &                  array_on_v_points)
!
      Integer                                                           &
               ! INTENT (IN)
     & row_length                                                       &
     &,rows                                                             &
     &,nrows                                                            &
     &,levels                                                           &
     &,halo_i                                                           &
     &,halo_j

      Real                                                              &
               ! INTENT (IN)
     & array_on_p_points(1-halo_i:row_length+halo_i,                    &
     &                   1-halo_j:rows+halo_j,levels)

      Real                                                              &
               ! INTENT (OUT)
     & array_on_v_points(row_length, nrows, levels)

      Integer                                                           &
               ! local variables
     & i,j,k

      do k=1,levels
        do j=1,nrows
          do i=1,row_length

#if !defined(SCMA)
            array_on_v_points(i,j,k)= 0.5 *                             &
     &    ( array_on_p_points(i,j,k) + array_on_p_points(i,j+1,k) )
#else
            array_on_v_points(i,j,k)= array_on_p_points(i,j,k)
#endif

          end do

        end do
      end do

      return
      END SUBROUTINE P_TO_V
#endif
