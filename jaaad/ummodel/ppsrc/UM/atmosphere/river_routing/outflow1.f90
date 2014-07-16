
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

!

      subroutine outflow1(sto, drunin, rc, dt                           &
     &     , igrcn, iseq, inextx, inexty, nseqmax, nx, ny, jmax         &
     &     , sto2, din, dout                                            &
     &     , drunall, drivall, stoall, sto2all, doutall, dinall)

      IMPLICIT NONE
!
!  Calculate the storage in the next time step based on the storage;
!  of current time step; sto
!
!  Runoff generated at the grid box is constant during the time step;
!  transfer coefficient is also given; rc
!  For 'runoff drunin' in TRIP seapoints just add to TRIP seapoint
!
      integer nx, ny, i, j, inextx(nx, ny), inexty(nx, ny)              &
     &     , igrcn(nx, ny), jmax, iseq(nx, ny), nseq, nseqmax
      real sto(nx, ny), drunin(nx, ny), rc(nx, ny), dt                  &
     &     , sto2(nx, ny), din(nx, ny), dout(nx, ny)                    &
     &     , drunall, drivall, stoall, sto2all, doutall, dinall
!
! DEPENDS ON: setrval
      call setrval(nx, ny, din, 0.0, jmax)
      drivall = 0.0
      drunall = 0.0
      stoall = 0.0
      dinall = 0.0
      sto2all = 0.0
      doutall = 0.0
      sto2 = 0.0
!********************************************************************
!
      do nseq = 1, nseqmax
        do j = 1, jmax
          do i = 1, nx
            if (iseq(i,j) == nseq) then
              din(i,j) = din(i, j) + drunin(i,j)
              drunall = drunall + drunin(i,j)
              if (rc(i, j) >  0.0) then
!
                sto2(i, j) = sto(i,j) * exp(-(rc(i,j)*dt))              &
     &               + (1.0 - exp(-(rc(i,j)*dt))) * din(i,j) / rc(i,j)
                dout(i, j) = (sto(i,j)-sto2(i,j))/dt + din(i,j)
!
                stoall = stoall + sto(i,j)
                sto2all = sto2all + sto2(i,j)
                doutall = doutall + dout(i,j)*dt
                dinall  = dinall  + din(i,j)*dt
!
                if ((igrcn(i,j) >= 1).and.(igrcn(i,j) <= 8)) then
                  din(inextx(i,j),inexty(i,j))                          &
     &                 = din(inextx(i,j),inexty(i,j)) + dout(i,j)
                  drivall = drivall + dout(i,j)
                end if
!
              else
                write (6, *) '!! rc is not positive ', i, j,rc(i,j)
              end if
            end if
          end do
        end do
      end do
!
      END SUBROUTINE outflow1
!

! ******************************COPYRIGHT******************************
