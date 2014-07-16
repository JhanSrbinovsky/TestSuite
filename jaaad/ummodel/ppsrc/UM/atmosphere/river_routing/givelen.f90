
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
      real*8 function givelen(rx1, ry1, rx2, ry2)
!
!     give the length (km) between (rx1, ry1) to (rx2, ry2)
!     sphere approximation is applied.
!
!     see page 621 of Rika-Nenpyo (1995)
!
!     from 26/August/1996, by Taikan OKI
!
      real*8 givelat, givelon, giverade, rlat, dx, dy, re
      real rx1, rx2, ry1, ry2, dlon

      dlon = abs(rx1 - rx2)
      if (dlon >= 180.0) dlon = abs(360.0 - dlon)

      if (rx1 == rx2) then
        rlat = (ry1+ry2) / 2.0
! DEPENDS ON: givelat
        givelen = givelat(rlat) * abs(ry1-ry2)
      else if (ry1 == ry2) then
        rlat = ry1
! DEPENDS ON: givelon
        givelen = givelon(rlat) * dlon

      else
        rlat = (ry1+ry2) / 2.0
! DEPENDS ON: giverade
        re = giverade(rlat)
! DEPENDS ON: givelon
        dx = givelon(rlat) * dlon / re
! DEPENDS ON: givelat
        dy = givelat(rlat) * abs(ry1 - ry2) / re
!        write (6, '(5f15.7)') dx*re, dy*re, dx, dy, re
        givelen = acos(cos(dx)*cos(dy)) * re
! Remove dbl as not used on T3E
!        givelen = acos(dble(cos(dx)*cos(dy))) * re
!        write (6, '(5f15.7)') givelen, dsqrt(dx*dx+dy*dy)*re
      end if

      END FUNCTION givelen


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
