#if defined(C92_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Description:
!     subroutine P_TO_V_SEA  for calculating sea variables held at
!     p points at v points, on the new dynamics grid
!
!     This routine does interior points of array not halos,
!     but requires halo information to be set.
!
!     global code has E-W wrap around
!     LAM code has most east u point set to zero
!
!!!  Model            Modification history:
!!! version  Date
!!!  5.3   23/04/01   New Deck         N.Gedney
!    6.2   03/02/06   Move to C92_2A.  P.Selwood
!
!
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! Subroutine Interface:
      Subroutine P_TO_V_SEA(array_on_p_points,fland_on_p_points,        &
     &                  row_length,rows,nrows,                          &
     &                  levels,halo_i,halo_j,                           &
     &                  array_on_v_points)
!
      IMPLICIT NONE
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
     &                   1-halo_j:rows+halo_j,levels)                   &
     &,fland_on_p_points(1-halo_i:row_length+halo_i,                    &
     &                   1-halo_j:rows+halo_j)   !IN land fraction
!                                                !   on p points.

      Real                                                              &
               ! INTENT (OUT)
     & array_on_v_points(row_length, nrows, levels)

      Integer                                                           &
               ! local variables
     & i,j,k

      Real                                                              &
     & fsea_total    ! total land fraction on u grid point

      do k=1,levels
        do j=1,nrows
          do i=1,row_length
            fsea_total=2.0-(fland_on_p_points(i,j)                      &
     &         +fland_on_p_points(i,j+1))

            if(fsea_total >  0.0)then
             array_on_v_points(i,j,k)= 1./fsea_total*                   &
     &    ( (1.-fland_on_p_points(i,j))*array_on_p_points(i,j,k) +      &
     &      (1.-fland_on_p_points(i,j+1))*array_on_p_points(i,j+1,k) )
            else
             array_on_v_points(i,j,k)= 0.0
            endif

          end do
        end do
      end do

      return
      END SUBROUTINE P_TO_V_SEA
#endif
