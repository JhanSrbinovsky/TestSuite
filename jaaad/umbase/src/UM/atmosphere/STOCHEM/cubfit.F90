#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE CUBFIT(x,y,xout,yout)
!-----------------------------------------------------------------------
!
!     PURPOSE:  Polynomial cubic interpolation in vertical direction.
!               Fit cubic through 2 points and two derivatives.
!
!     INPUTS:   X,Y,Xout     X are the positions of the know Y values,
!                          Xout the required position.
!     OUTPUTS:  Yout         Yout is the required value.
!
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  4.4    28/01/97  Created.  W.J. Collins
!  6.1    20/10/04  Minor improvement to code. M.G. Sanderson
!  6.2    28/03/06  Minor changes for vn6.2  M.G. Sanderson
!-
!VVV  V5.0  CUBFIT 20/vi/01
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER, PARAMETER             :: n=4  ! 4 points to fit cubic
      REAL, DIMENSION(n), INTENT(IN) :: x    ! x coords of points
      REAL, DIMENSION(n), INTENT(IN) :: y    ! y values at points
      REAL, INTENT(IN)               :: xout ! x coord or required point
      REAL, INTENT(OUT)              :: yout ! y value at required point

      REAL :: x2    ! x at point 2
      REAL :: x3    ! x difference between points 2 and 3
      REAL :: x3_2   ! x3 squared
      REAL :: x3_3   ! x3 cubed
      REAL :: a     ! x cubed coefficient
      REAL :: b     ! x squared coefficient
      REAL :: c     ! x coefficient
      REAL :: d     ! constant
      REAL :: dydx2 ! dy/dx at point 2
      REAL :: dydx3 ! dy/dx at point 3
      REAL :: dx    ! x difference between xout and point 2

      dydx2=(y(3)-y(1))/(x(3)-x(1))     ! derivative at point 2
      dydx3=(y(4)-y(2))/(x(4)-x(2))     ! derivative at point 3
      ! Define coordinate system such that x=0 coincides with x2.
      x2=x(2)
      x3=x(3)-x2
      ! fit to point and derivative
      ! at point 2 y   = d =y(2)
      ! and        dydx= c =dydx2
      d=y(2)
      c=dydx2
      x3_2=x3*x3
      x3_3=x3_2*x3
      ! fit to point and derivative
      ! at point 3 y    =  ax3_3+ bx3_2+cx3+d=y(3)
      ! and        dy/dx= 3ax3_2+2bx3  +c    =dydx3
      ! solve for 2 unknowns to give:
      b=(3.0*y(3)-x3*dydx3-2.0*c*x3-3.0*d)/x3_2
      a=(y(3)-b*x3_2-c*x3-d)/x3_3
      dx=xout-x2
      yout=d+dx*(c+dx*(b+dx*a))
      END SUBROUTINE CUBFIT
#endif
