#if defined(A26_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************





      subroutine setrval(nx, ny, rdim, rdat, jmax)

      IMPLICIT NONE
!
!     set uniform value of rdat into rdim(i,j)
!
      integer nx, ny, i, j, jmax
      real rdim(nx, ny), rdat
!
      do i = 1, nx
        do j = 1, jmax
          rdim(i,j) = rdat
        end do
      end do
      END SUBROUTINE setrval





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
