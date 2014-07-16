
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Interpolates v-component from the p-grid to the v-grid.
!
! Subroutine Interface:

      Subroutine lbc_p_to_v ( v_p,                                      &
     &                        v,                                        &
     &                        row_length,                               &
     &                        rows,                                     &
     &                        v_rows,                                   &
     &                        offx,                                     &
     &                        offy )

      Implicit NONE
!
! Description:
!   Interpolates v-component from the p-grid to the v-grid.
!
! Method:
!   Linear Interpolation of v-comp from p-grid to v-grid.
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

      Integer :: row_length    !  No of points
      Integer :: rows          !  No of p rows
      Integer :: v_rows        !  No of v rows
      Integer :: offx, offy    !  Halo sizes

!     v_p : v-component on p-grid
!     v   : v-component on v-grid

      Real    :: v_p(1-offx:row_length+offx, 1-offy:  rows+offy)
      Real    :: v  (1-offx:row_length+offx, 1-offy:v_rows+offy)

! Local scalars:

      Integer :: i,j           !  Loop indices
!
!- End of header

!     ----------------------------------------------
!     Interpolate v from p-grid to v-grid ; v_p to v
!     ----------------------------------------------

      If (v_rows < rows ) Then

        Do i = 1-offx, row_length+offx
          Do j = 1-offy, v_rows+offy
            v(i,j) = ( v_p(i,j) +  v_p(i,j+1) ) * 0.5
          End Do
        End Do

      Else  !  v_rows == rows

        Do i = 1-offx, row_length+offx
          Do j = 1-offy, v_rows+offy-1
            v(i,j) = ( v_p(i,j) +  v_p(i,j+1) ) * 0.5
          End Do
          v(i,v_rows+offy) = v_p(i,rows+offy)
        End Do

      End If

      Return
      END SUBROUTINE lbc_p_to_v
