
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Interpolates u-component from the p-grid to the u-grid.
!
! Subroutine Interface:

      Subroutine lbc_p_to_u (u_p,                                       &
     &                       u,                                         &
     &                       row_length,                                &
     &                       rows,                                      &
     &                       offx,                                      &
     &                       offy )

      Implicit NONE
!
! Description:
!   Interpolates u-component from the p-grid to the u-grid.
!
! Method:
!   Linear Interpolation of u-comp from p-grid to u-grid.
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

!     u_p : u-component on p-grid
!     u   : u-component on u-grid

      Real    :: u_p(1-offx:row_length+offx, 1-offy:rows+offy)
      Real    :: u  (1-offx:row_length+offx, 1-offy:rows+offy)

! Local scalars:

      Integer :: i,j           ! Loop indices
!
!- End of header

!     ----------------------------------------------
!     Interpolate u from p-grid to u-grid ; u_p to u
!     ----------------------------------------------

      Do j = 1-offy, rows+offy

        Do i = 1-offx, row_length+offx-1
          u(i,j) = ( u_p(i,j) +  u_p(i+1,j) ) * 0.5
        End Do

      End Do

      Return
      END SUBROUTINE lbc_p_to_u
