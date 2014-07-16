#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      REAL FUNCTION HINTERP(x,d1,d2,i,l)
!-----------------------------------------------------------------------
!-
!-   Purpose and Methods : Do 2-D interpolation of Orography
!-
!-   Returned value  :
!-   Inputs  : X,D1,D2,I,L
!-   Outputs :
!-   Controls:
!-
!-   Created  12-JUL-2001   Colin Johnson from Interp
!
! History:
! Version   Date                    Comment
!
!  6.2    28/03/06  Minor changes for vn6.2  M.G. Sanderson
!
!-
!VVV  V5.2  HINTERP 12/VII/01
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!-----------------------------------------------------------------------
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: x  ! Interp array
      REAL,                           INTENT(IN) :: d1 ! Rel position
      REAL,                           INTENT(IN) :: d2 ! Rel position
      INTEGER,     INTENT(IN)                    :: i  ! Indicies
      INTEGER,     INTENT(IN)                    :: l  ! Indicies

      REAL               :: x1
      REAL               :: x2
      REAL               :: x3
      REAL, DIMENSION(4) :: xr

      xr(1)=x(i,l)
      xr(2)=x(i,l+1)
      xr(3)=x(i+1,l)
      xr(4)=x(i+1,l+1)

! Do linear interpolation in horizontal directions.
      x1=xr(1)+((xr(3)-xr(1))/dlongm)*d1
      x2=xr(2)+((xr(4)-xr(2))/dlongm)*d1
      x3=x1+((x2-x1)/dlatm)*d2
      hinterp=x3

      END FUNCTION HINTERP
#endif
