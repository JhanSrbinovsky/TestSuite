
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Interpolates u-component from the u-grid to the p-grid.
!
! Subroutine Interface:

      Subroutine lbc_u_to_p ( u,                                        &
     &                        u_p,                                      &
     &                        row_length,                               &
     &                        rows,                                     &
     &                        offx,                                     &
     &                        offy )

      Implicit NONE
!
! Description:
!   Interpolates u-component from the u-grid to the p-grid.
!
! Method:
!   Linear Interpolation of u-comp from u-grid to p-grid.
!
! Original Author : Dave Robinson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.5    16/10/02  Original code. Dave Robinson
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Subroutine arguments

      Integer :: row_length    ! No of points
      Integer :: rows          ! No of rows
      Integer :: offx, offy    ! Halo sizes

!     u   : u-component on u-grid
!     u_p : u-component on p-grid

      Real    :: u  (1-offx:row_length+offx, 1-offy:rows+offy)
      Real    :: u_p(1-offx:row_length+offx, 1-offy:rows+offy)

! Local scalars:

      Integer :: i,j    !  Loop indices
!
!- End of header

!     ----------------------------------------------
!     Interpolate u from u-grid to p-grid ; u to u_p
!     ----------------------------------------------

      Do j = 1-offy, rows+offy

        Do i = 1-offx+1, row_length+offx
          u_p(i,j) = ( u(i-1,j) +  u(i,j) ) * 0.5
        End Do

        ! Copy for first column of P points.
        u_p(1-offx,j) = u(1-offx,j)

      End Do

      Return
      END SUBROUTINE lbc_u_to_p
