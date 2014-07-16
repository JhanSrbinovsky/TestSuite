
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************



      subroutine setnext(nx, ny, igrcn, inextx, inexty)

      IMPLICIT NONE
!
!     set the destination grid point
!
!     (i, j) ===>  (inextx(i,j), inexty(i,j))
!     at river mouth : pointing itself
!     at sea         : 0
!
      integer nx, ny
      integer igrcn(nx, ny), inextx(nx, ny), inexty(nx, ny)
      integer i, j, irnxtx, irnxty, inow
!
      INTEGER DX(9),DY(9)
      DATA DX/0, 1, 1, 1, 0, -1, -1, -1, 0/
      DATA DY/-1, -1, 0, 1, 1, 1, 0, -1, 0/
      DO J = 1, NY
        DO I = 1, NX
          IF (IGRCN(I,J)  >= 1 .AND. IGRCN(I,J)  <= 9) THEN
            inextx(I,J) = I + DX(IGRCN(I,J) )
            IF (inextx(I,J)  ==  0) inextx(I,J) = NX
            IF (inextx(I,J)  ==  NX+1) inextx(I,J) = 1

            inexty(I,J) = J + DY(IGRCN(I,J) )
            IF (inexty(I,J)  ==  0) THEN
              WRITE(6,*)' SETNEXT:ERROR in RIVER DIRECTION FILE'
              inexty(I,J) = 1
            ENDIF
            IF (inexty(I,J)  ==  NY+1) THEN
              WRITE(6,*)' SETNEXT:ERROR in RIVER DIRECTION FILE'
              inexty(I,J) = NY
            ENDIF
          ELSE
            inextx(I,J) = 0
            inexty(I,J) = 0
          ENDIF

        ENDDO
      ENDDO
      END SUBROUTINE setnext







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
