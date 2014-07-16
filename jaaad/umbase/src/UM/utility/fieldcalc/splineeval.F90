#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to calculate Max Wind value and height


!=======================================================================


! ====================================================================



! ====================================================================



!=======================================================================


!=======================================================================


!=======================================================================
SUBROUTINE SplineEval (n, u, x, y, b, c, d, seval)

! Description:
!
! Method:
!
! Owner: Dave Robinson
!
! History:
! Version Date     Comment
! ------- ----     -------
! 1.0     02/05/03 Original Code.  Sara James
! 6.0     12/09/03 Code implemented into UM. Dave Robinson
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: n
REAL,    INTENT(IN) :: u
REAL,    INTENT(IN) :: x(n), y(n), b(n), c(n), d(n)
REAL,    INTENT(OUT) :: seval

! Local constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "SplineEval"

! Local variables:
INTEGER :: jl,jm,ju
REAL :: dx

!-----------------------------------------------------------------------
!------- ypint is the derivative of the function at position u
!------- this subroutine evaluates the cubic spline function
!
!  seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
!    where  x(i) < u < x(i+1), using horner's rule
!  if  u < x(1) then  i = 1  is used.
!  if  u >= x(n) then  i = n  is used.
!  input..
!    n = the number of data points
!    u = the abscissa at which the spline is to be evaluated
!    x,y = the arrays of data abscissas and ordinates
!    b,c,d = arrays of spline coefficients computed by spline
!  a binary search is performed to determine the proper interval.
!-----------------------------------------------------------------------

jl = 1
ju = n+1
DO WHILE ((ju-jl) .gt. 1)
  jm = (jl+ju)/2
  IF ( u <  x(jm) ) ju = jm
  IF ( u >= x(jm) ) jl = jm
ENDDO
!------------------------------------------------  evaluate spline
dx = u - x(jl)
seval = y(jl) + dx*(b(jl) + dx*(c(jl) + dx*d(jl)))
!-----------------------------------------------------------------

END SUBROUTINE SplineEval

#endif
