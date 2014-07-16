
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!     subroutine U_TO_p  for calculating variables held at u points
!     at p points, on the new dynamics grid
!
!     This routine does interior points of array not halos,
!     but requires halo information to be set.
!LL   Model Date      Comments/Programmer
!LL   Vn
!LL   5.1   09/02/00  Use DOMTYP parameters                    P.Burton
! 5.3  17/10/01 Changes required for Single Column Model
!                                             Z. Gardner
! 6.2  03/02/06 Moved to c92_2a. P.Selwood

!
      Subroutine U_TO_P(array_on_u_points,row_length,rows,levels,       &
     &                  off_x, off_y, model_domain,                     &
     &                  at_extremity, array_on_p_points)
!
      Integer                                                           &
               ! INTENT (IN)
     & row_length                                                       &
     &,rows                                                             &
     &,levels                                                           &
     &,model_domain                                                     &
     &,off_x                                                            &
     &,off_y

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Parameters
      Integer                                                           &
     &   PNorth,                                                        &
                      ! North processor address in the neighbor array
     &   PEast,                                                         &
                      ! East processor address in the neighbor array
     &   PSouth,                                                        &
                      ! South processor address in the neighbor array
     &   PWest,                                                         &
                      ! West processor address in the neighbor array
     &   NoDomain     ! Value in neighbor array if the domain has
                      !  no neighbor in this direction. Otherwise
                      !  the value will be the tid of the neighbor
      Parameter (                                                       &
     &   PNorth   = 1,                                                  &
     &   PEast    = 2,                                                  &
     &   PSouth   = 3,                                                  &
     &   PWest    = 4,                                                  &
     &   NoDomain = -1)

      Real                                                              &
               ! INTENT (IN)
     & array_on_u_points(1-off_x:row_length+off_x,                      &
     &                   1-off_y:rows+off_y, levels)

      Real                                                              &
               ! INTENT (OUT)
     & array_on_p_points(row_length, rows, levels)

      Integer                                                           &
               ! local variables
     & i,j,k                                                            &
     &, j0, j1

! DOMTYP contains different model domain types
!
! Author : P.Burton
! History:
! Version  Date      Comment.
! 5.0      15/04/99  New comdeck
! 5.2      15/11/00  add bi_cyclic_lam domain   A. Malcolm

      INTEGER,PARAMETER:: mt_global        = 1
      INTEGER,PARAMETER:: mt_lam           = 2
      INTEGER,PARAMETER:: mt_cyclic_lam    = 3
      INTEGER,PARAMETER:: mt_bi_cyclic_lam = 4
      INTEGER,PARAMETER:: mt_single_column = 5
! DOMTYP end
      j0 = 1
      j1 = rows
      If (model_domain  ==  mt_global) Then
! Do not do poles as this will be done by polar vector wind.
        If (at_extremity(PSouth) ) Then
          j0 = 2
        End If
        If (at_extremity(PNorth) ) Then
          j1 = rows - 1
        End If
      End If

      Do k=1,levels
        Do j= j0, j1

          Do i= 1,row_length
            array_on_p_points(i,j,k)= 0.5 *                             &
     &    ( array_on_u_points(i,j,k) + array_on_u_points(i-1,j,k) )

          End Do

        End Do
      End Do

      return
      END SUBROUTINE U_TO_P
