
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


      subroutine initial0(rin, nx, ny, jmax, rmiss                      &
     &     , igrcn, iseq, nseqmax, inextx, inexty, rlen, sto, area)

      IMPLICIT NONE
!
!     Setting Initial Conditions
!
      integer nx, ny, jmax, i, j                                        &
     &     , nseqmax                                                    &
     &     , igrcn(nx, ny), iseq(nx, ny), inextx(nx, ny)                &
     &     , inexty(nx, ny)
      real rin(nx, ny), rmiss
      real rlen(nx, ny), sto(nx, ny), area(nx, ny)
!
      nseqmax = 0
      do i = 1, nx
        do j= 1, ny
          if ( iseq(i,j) >  nseqmax) nseqmax = iseq(i,j)
        end do
      end do
      write(6, '(" initial0: maximum down flow= ",i4)') nseqmax
!
!      write(6, *) 'initial0: set down stream.'

! DEPENDS ON: setnext
      call setnext(nx, ny, igrcn, inextx, inexty)
!
!     set the distance between grids [m]
!
!      write(6, *) 'initial0: set length.'
! DEPENDS ON: setlen
      call setlen (nx, ny, igrcn, inextx, inexty, rlen, jmax            &
     &, rmiss)
!
!      write(6, *) 'initial0: set area size'
! DEPENDS ON: setarea
      call setarea(nx, ny, area, jmax)
!
!     Initial Condition
!a
      END SUBROUTINE initial0
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
