
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
      real*8 function givelat(rlat)
!
!     give the length (km) of 1 degree latitude (dy)
!     at the latitude of rlat(degree)
!
!     see page 621 of Rika-Nenpyo (1995)
!
!     from 23/April/1996, by Taikan OKI
!
      real*8 rlat, e2, pi, ra
      real sind
      external sind
      data ra, e2, pi/6378.136, 0.006694470,3.14159265358979/

      givelat = pi / 180.0 * ra * ( 1 - e2)                             &
! DEPENDS ON: sind
     &     / sqrt((1 - e2*sind(rlat)*sind(rlat))**3)

      END FUNCTION givelat


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
