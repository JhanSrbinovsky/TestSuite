
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Interpolates the v-component from the v-grid to the p-grid.
!
! Subroutine Interface:

      Subroutine lbc_v_to_p ( v,                                        &
     &                        v_p,                                      &
     &                        row_length,                               &
     &                        rows,                                     &
     &                        v_rows,                                   &
     &                        offx,                                     &
     &                        offy)

      Implicit NONE
!
! Description:
!   Interpolates the v-component from the v-grid to the p-grid.
!
! Method:
!   Linear Interpolation of v-comp from v-grid to p-grid.
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
      Integer :: rows          ! No of p rows
      Integer :: v_rows        ! No of v rows
      Integer :: offx, offy    ! Halo sizes

!     v   : v-component on v-grid
!     v_p : v-component on p-grid

      Real    :: v  (1-offx:row_length+offx, 1-offy:v_rows+offy)
      Real    :: v_p(1-offx:row_length+offx, 1-offy:  rows+offy)

! Local scalars:

      Integer :: i,j      ! Loop indices

!- End of header

!     -----------------------------------
!     Interpolate v from v-grid to p-grid
!     -----------------------------------

      Do i = 1-offx, row_length+offx

        Do j = 1-offy+1, v_rows+offy
          v_p(i,j) = ( v(i,j-1) +  v(i,j) ) * 0.5
        End Do

        ! Copy for bottom row of v_p
        v_p(i,1-offy)    = v(i,1-offy)

        ! Copy for top row of v_p
        If (v_rows < rows) Then
          v_p(i,rows+offy) = v(i,v_rows+offy)
        End If

      End Do

      Return
      END SUBROUTINE lbc_v_to_p
