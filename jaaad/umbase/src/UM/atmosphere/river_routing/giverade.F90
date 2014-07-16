#if defined(A26_1A)
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
      real*8 function giverade(rlat)
!
!     give the equivalent radius of the earth (Re)
!     at the latitude of rlat(degree)
!
!     see page 621 of Rika-Nenpyo (1995)
!
!     from 26/August/1996, by Taikan OKI
!
      real*8 rlat, e2, ra, rn
      real sind
      external sind
      data ra, e2/6378.136, 0.006694470/


! DEPENDS ON: sind
      rn = ra / sqrt(1 - e2*sind(rlat)*sind(rlat))
! DEPENDS ON: sind
      giverade = rn * sqrt(1.0 -2*e2*sind(rlat)                         &
! DEPENDS ON: sind
     &     + e2*e2*sind(rlat)*sind(rlat))
!      giverade = ra / sqrt(1 - e2*sind(rlat)*sind(rlat))

      END FUNCTION giverade


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
