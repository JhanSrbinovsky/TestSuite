
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************










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

      subroutine mmd2kgs(rin, igrcn, din, area, nx, ny, rmiss, jmax)

      IMPLICIT NONE
!
!     convert rin [mm/day] fo din [kg] using area [m^2]
!     [mm/day] x [m^2] = [10^(-3) m^3/day]
!                       ===> [kg/day] / (3600*24) --> [kg/s]
!
      integer nx, ny, i, j, igrcn(nx, ny), jmax
      real rin(nx, ny)
      real din(nx, ny), area(nx, ny), rc, rmiss
      data rc/86400.0/

!      write (6, *) 'mmd2kgs : rc=', rc
      do j = 1, jmax
        do i = 1, nx
!          if ((igrcn(i,j) >= 1).and.(igrcn(i,j) <= 9).and.
!     $         (rin(i,j) /= rmiss)) then
           if ((rin(i,j) /= rmiss)) then
            din(i,j) = rin(i,j) * area(i,j) / rc
!            write (6, *) 'mmd2kgs: ', i, j
!     $           , din(i,j), rin(i,j), area(i,j)
          else
            din(i,j) = 0.0
          end if
        end do
      end do
      END SUBROUTINE mmd2kgs
!

!

! ******************************COPYRIGHT******************************
