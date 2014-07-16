
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************









      character*32 function c0fmt(indat, ncol)
!
!     convert integer into ncol characters with 0 head.
!
      character cdum*32, cd(0:9), cdig(32)*1
      integer indat, idat, ncol, i, idig(32)
      data cd/'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/
!
      idat = indat
!      write (6, *) 'idat, ncol=', idat, ncol
      cdum = ''
      if (ncol <= 31) then
        idig(ncol+1) = 0
        do i = ncol, 1, -1
          idig(i) = int(aint(real(idat)/(10**(i-1))))
          cdig(i) = cd(int(aint(real(idat)/(10**(i-1)))))
!          write(6, *) cdig(i)
          idat = idat - idig(i)*(10**(i-1))
        end do
      end if
!
      write (cdum,'(32a1)') (cdig(i), i= ncol, 1, -1)
!      write(6, *) "cdum = ", cdum
      c0fmt = cdum
!      write(6, *) "c0fmt = ", c0fmt
      END FUNCTION c0fmt

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
