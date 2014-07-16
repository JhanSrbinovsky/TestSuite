
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************








      real function arealat1(ilat)
!
!     returns an area of 1 degree x 1 degree mesh in km^2
!
!     from 05/Jul/1996 by Taikan  OKI
!
!     use the results from calcarea (function arealat.f)
!
!     arealat1(89) : 88.0--89.0degree
!     arealat1(1)  : 0.0 -- 1.0 degree
!
      integer ilat
      real area(90)

      data area/                                                        &
     &    12308.46,                                                     &
     &    12304.81,                                                     &
     &    12297.51,                                                     &
     &    12286.57,                                                     &
     &    12271.98,                                                     &
     &    12253.75,                                                     &
     &    12231.89,                                                     &
     &    12206.39,                                                     &
     &    12177.27,                                                     &
     &    12144.53,                                                     &
     &    12108.18,                                                     &
     &    12068.23,                                                     &
     &    12024.68,                                                     &
     &    11977.56,                                                     &
     &    11926.85,                                                     &
     &    11872.59,                                                     &
     &    11814.79,                                                     &
     &    11753.44,                                                     &
     &    11688.58,                                                     &
     &    11620.22,                                                     &
     &    11548.37,                                                     &
     &    11473.05,                                                     &
     &    11394.28,                                                     &
     &    11312.08,                                                     &
     &    11226.47,                                                     &
     &    11137.46,                                                     &
     &    11045.09,                                                     &
     &    10949.37,                                                     &
     &    10850.33,                                                     &
     &    10748.00,                                                     &
     &    10642.39,                                                     &
     &    10533.54,                                                     &
     &    10421.47,                                                     &
     &    10306.22,                                                     &
     &    10187.81,                                                     &
     &    10066.27,                                                     &
     &     9941.64,                                                     &
     &     9813.95,                                                     &
     &     9683.23,                                                     &
     &     9549.52,                                                     &
     &     9412.85,                                                     &
     &     9273.26,                                                     &
     &     9130.79,                                                     &
     &     8985.48,                                                     &
     &     8837.37,                                                     &
     &     8686.49,                                                     &
     &     8532.90,                                                     &
     &     8376.63,                                                     &
     &     8217.73,                                                     &
     &     8056.24,                                                     &
     &     7892.22,                                                     &
     &     7725.70,                                                     &
     &     7556.73,                                                     &
     &     7385.37,                                                     &
     &     7211.67,                                                     &
     &     7035.67,                                                     &
     &     6857.43,                                                     &
     &     6677.01,                                                     &
     &     6494.45,                                                     &
     &     6309.80,                                                     &
     &     6123.14,                                                     &
     &     5934.51,                                                     &
     &     5743.97,                                                     &
     &     5551.58,                                                     &
     &     5357.39,                                                     &
     &     5161.48,                                                     &
     &     4963.90,                                                     &
     &     4764.71,                                                     &
     &     4563.97,                                                     &
     &     4361.75,                                                     &
     &     4158.12,                                                     &
     &     3953.13,                                                     &
     &     3746.85,                                                     &
     &     3539.35,                                                     &
     &     3330.69,                                                     &
     &     3120.95,                                                     &
     &     2910.18,                                                     &
     &     2698.46,                                                     &
     &     2485.86,                                                     &
     &     2272.44,                                                     &
     &     2058.27,                                                     &
     &     1843.43,                                                     &
     &     1627.99,                                                     &
     &     1412.01,                                                     &
     &     1195.56,                                                     &
     &      978.73,                                                     &
     &      761.56,                                                     &
     &      544.15,                                                     &
     &      326.56,                                                     &
     &      108.87/

      if ((ilat <= 90).and.(ilat >= 1)) then
        arealat1 = area(ilat)
!      else if ((ilat <= 180).and.(ilat >= 91)) then
!        arealat1 = area(ilat-90)
      else
        arealat1 = 0.0
      end if

      END FUNCTION arealat1


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
