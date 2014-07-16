
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************


      subroutine setcoef(nx, ny, rlen, rvel, ratmed, rc, jmax)

      IMPLICIT NONE
!
!     calculate coefficient 'c' in the routing model
!     basically by c = u/d, but actual d is assumed to d*ratmed
!     considering the meandering of the river.
!
!     c has dimension of [1/s]
!     sea grids: c = 0.0
!
      integer nx, ny, i, j, jmax
      real rlen(nx, ny), rvel(nx, ny), ratmed, rc(nx, ny)
!
      do i= 1, nx
        do j = 1, ny
          if (rlen(i,j) >  0.0) then
            rc(i, j) = rvel(i,j) / (rlen(i, j) * ratmed)
          else
            rc(i, j) = 0.0
          end if
        end do
      end do
      END SUBROUTINE setcoef








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
