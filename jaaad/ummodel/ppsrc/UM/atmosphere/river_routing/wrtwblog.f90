
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************





      subroutine wrtWBlog(iofile, iy, im, id, ih, ndev                  &
     &     , stoall, sto2all, dinall, doutall, drunall, drivall         &
     &     , dt)

      IMPLICIT NONE
!
!     write water balance monitoring
!
!     stoall    : Total river channel storage
!     sto2all   : Total river channel storage at the next time step
!     dinall    : Total inflow to the next grid for the next time step
!     doutall   : Total outflow to the next grid for the next time step
!   [ drunall*dt  : Total inflow to the grid from LSM]
!   [ drivall*dt  : Total inflow to the grid from surrounding grids]
!   [(stoall - sto2all + drunall*dt) / (135.3E12) : Mean runoff to the
!                                                   sea
      integer iofile, iy, im, id, ih, ndev
      real rorder, stoall, sto2all, dinall, doutall                     &
     &     , drunall, drivall, dt

      data rorder/1.0E12/
!
      write(iofile, *)' River Routing: ',                               &
     &'Error = Boxinflow - (sum of river inflow + box runoff)'
      write (iofile, 671)
 671  format ('          YY/MM/DD/HH ', '  S(t)  S(t+1)'                &
     &     , '     Din    Dout   Run  Rivflow  Error  Outflow(mm/y) ')
      write(iofile, 600) iy, im, id, ih, stoall/rorder, sto2all/rorder  &
     &           , dinall/rorder, doutall/rorder                        &
     &           , (drunall*dt) / rorder                                &
     &             , (drivall*dt) / rorder                              &
     &           , (dinall - (drivall+drunall)*dt)                      &
     &           , (stoall - sto2all + drunall*dt) / (135.3E12)         &
     &           * 365.0 * real(ndev)
 600  format("10^12 kg: ", i2,"/",i2,"/",i2, "/", i2                    &
!     $           , 6(' ',f7.1), f8.2)
     &           , 6(' ',f7.1), ' ', g8.2, f8.2)
      write (iofile, 602)                                               &
     &(stoall - sto2all + drunall*dt)* real(ndev)/86400,                &
     &(stoall - sto2all + drunall*dt) / (135.3E12)                      &
     &            * real(ndev)
 602  format(' Total outflow to sea or inland points is ',g9.2,         &
     &' Kg/s or ',g9.2,' mm/day')
      write(iofile, *)' River routing: Water balance monitoring'
      write(iofile, *)' (Tot water in - tot water out - diff. water',   &
     &' storage) summed for all gridboxes'
      write (iofile, 672)
 672  format ('          YY/MM/DD/HH ', '  S(t)  S(t+1)'                &
     &     , '     Din    Dout   Din-Dout-(S(t+1)- S(t))Kg')
      write(iofile, 601) iy, im, id, ih, stoall/rorder, sto2all/rorder  &
     &           , dinall/rorder, doutall/rorder                        &
     &           , (dinall - doutall -(sto2all - stoall))

 601  format("10^12 kg: ", i2,"/",i2,"/",i2, "/", i2                    &
     &           , 4(' ',f7.1), ' ', g8.2, f8.2)


      END SUBROUTINE wrtWBlog





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
