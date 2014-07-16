#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculate weights for a temporal filtering scheme.



!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
!+ Chebyshev polynomial function.

      REAL FUNCTION ChebPoly ( x, n ) ! in, in

      IMPLICIT NONE

! Method:
!
!   Use a recursion relation.
!
!   The Chebyshev polynomial T_n of order n is given by
!
!     T_0(x) = 1
!     T_1(x) = x
!
!     T_n(x) = 2xT_{n-1}(x) - T_{n-2}(x),  (n > 1)

! Function arguments:

      REAL,    INTENT(IN) :: x
      INTEGER, INTENT(IN) :: n  ! Order (non-negative).

! Local variables:

      INTEGER :: i

      REAL :: LastButOne_val,                                           &
     &        Last_val,                                                 &
     &        Latest_val

!- End of header ------------------------------------------------------

      IF (n == 0) ChebPoly = 1.0
      IF (n == 1) ChebPoly = x

      IF (n > 1) THEN

        LastButOne_val = 1.0
        Last_val       = x

        DO i = 2, n
          Latest_val     = 2.0 * x * Last_val - LastButOne_val
          LastButOne_val = Last_val
          Last_val       = Latest_val
        END DO

        ChebPoly = Latest_val

      END IF

      RETURN
      END FUNCTION ChebPoly
#endif
