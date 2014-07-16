#if defined(A26_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************


      subroutine setarea(nx, ny, area, jmax)

      IMPLICIT NONE
!
!     set area [m^2] of each grid box
!
      integer nx, ny, i, j, jmax
      real arealat1
      real area(nx, ny), atmp
!
      do j = 1, jmax
        if (j <= 90) then
! DEPENDS ON: arealat1
          atmp = arealat1(int(abs(91.0-j)))*1.0E06
        else
! DEPENDS ON: arealat1
          atmp = arealat1(int(abs(j-90.0)))*1.0E06
        end if
!
!        write (6, *) 'setarea: ', j, atmp
        do i = 1, nx
          area(i,j) = atmp
        end do
      end do
!
      END SUBROUTINE setarea








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
