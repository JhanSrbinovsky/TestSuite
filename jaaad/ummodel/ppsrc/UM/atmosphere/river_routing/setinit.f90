
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************



      subroutine setinit(nx, ny, area, rinit, sto, jmax)

      IMPLICIT NONE
!
!     assume water equivalent to rinit [mm] over 1 grid box is stored
!     as in the initial stage
!     [mm] x [m^2] = [10^(-3) m^3] ===> [kg]
!
      integer nx, ny, i, j, jmax
      real area(nx, ny), sto(nx, ny), rinit
!
      do j = 1, jmax
        do i = 1, nx
          sto(i,j) = rinit * area(i,j)
        end do
      end do
      END SUBROUTINE setinit







!
!  returns latituds at (iy) in (nla)
!
!     NOTICE: this is for GPCP/Xie and Arkin data set
!     which data is located at the center of 2.5 grid box
!
!     iy from SOUTH to NORTH.
!
!     from 23.Feb.1996, by Taikan OKI
!

!
!  returns longitude at (ix) in (nlo)
!
!     NOTICE: this is for GPCP/Xie and Arkin data set
!     which data is located at the center of 2.5 grid box
!
!     ix from west to east.
!
!     from 25.Feb.1996, by Taikan OKI
!

!


!


!


!


!

!
!  returns latituds at (iy) in (nla)
!
!     NOTICE: this is for GPCP/Xie and Arkin data set
!     which data is located at the center of 2.5 grid box
!
!     iy from SOUTH to NORTH.
!
!     from 23.Feb.1996, by Taikan OKI
!

!

!

! ******************************COPYRIGHT******************************
