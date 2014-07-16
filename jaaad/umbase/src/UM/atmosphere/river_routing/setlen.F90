#if defined(A26_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************



      subroutine setlen (nx, ny, igrcn, inextx, inexty, rlen            &
     &, jmax, rmiss)

      IMPLICIT NONE
!
!     length from (i, j) to the destination in [m]
!
!     river mouth : distance to 1 grid north
!     sea         : 0.0
!
      integer nx, ny, jmax
      integer inextx(nx, ny), inexty(nx, ny), igrcn(nx, ny), i, j
      real getlat0, getlon0, rx, ry, rx2, ry2
      real rlen(nx, ny), rmin, rmax, givelen, rmiss
!
      rx2 = 0.0
      ry2 = 0.0
      rmin = 1.0E10
      rmax = 0.0
!
      do j = 1, jmax
! DEPENDS ON: getlat0
        ry = getlat0(ny-j+1, ny)
        do i = 1, nx
! DEPENDS ON: getlon0
          rx = getlon0(i, nx)
!
          if ((igrcn(i,j) >= 1).and.(igrcn(i,j) <= 8)) then
! DEPENDS ON: getlon0
            rx2 = getlon0(inextx(i,j), nx)
! DEPENDS ON: getlat0
            ry2 = getlat0(ny-inexty(i,j)+1, ny)
! DEPENDS ON: givelen
            rlen(i, j) = givelen(rx, ry, rx2, ry2) * 1000.0
          else if (igrcn(i,j) == 9) then
! DEPENDS ON: getlat0
            ry2 = getlat0(ny-(j-1)+1, ny)
! DEPENDS ON: givelen
            rlen(i, j) = givelen(rx, ry, rx, ry2) * 1000.0
          else
            rlen(i,j) = 0.0
          end if
          if (igrcn(i,j) /= rmiss) then
            if (rlen(i,j) <  rmin) rmin = rlen(i,j)
            if (rlen(i,j) >  rmax) rmax = rlen(i,j)
          end if
!
        end do
      end do
!
      write (6, *) 'setlen: rmax, rmin = ', rmax, rmin
!
      END SUBROUTINE setlen







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
#endif
