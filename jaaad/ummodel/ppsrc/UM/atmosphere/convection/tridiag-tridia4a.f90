
! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE TRIDIAG------------------------------------
!LL
!LL  PURPOSE:  SOLVES THE EQUATIONS A.X=Y, WHERE A IS A
!LL            TRIDIAGONAL MATRIX
!LL
!LL  CALLED FROM: DEEP_GRAD_STRESS
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL
!LL  MODEL            MODIFICATION HISTORY:
!LL VERSION  DATE
!LL   5.4  6/8/2002   New deck created for convection scheme 4A
!LL                                     A.L.M. Grant
!     5.5  20/02/03   Replaced #ENDIF with #endif.      P.Dando
!
!LL   6.2  03/02/05  Added section 5A. R A Stratton
!     6.2  05/11/04   Remove hardwiring of gam array dimension to 38
!                     R. A. Stratton
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL
!LL  DOCUMENTATION : NONE
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE TRIDIAG(A,B,C,R,U,N)
!
      IMPLICIT NONE
!
! VARIABLES INPUT
!
      INTEGER N               ! SIZE OF VECTORS
!
      REAL A(N),                                                        &
                              ! COMPONENTS OF TRIDIAGONAL MATRIX
     &     B(N),                                                        &
     &     C(N),                                                        &
     &     R(N)               ! RHS OF LINEAR EQUATIONS
!
! VARIABLES THAT ARE OUTPUT
!
      REAL U(N)               ! SOLUTION VECTOR
!
! LOCAL VARIABLES
!
      INTEGER J
!
      REAL GAM(N),BET      ! DIMENSION OF GAM SHOULD be the same as N

!
      BET=B(1)
      U(1)=R(1)/BET
      DO J=2,N
       GAM(J)=C(J-1)/BET
       BET=B(J)-A(J)*GAM(J)
!       if (bet q 0) stop
       U(J)=(R(J)-A(J)*U(J-1))/BET
      END DO
      DO J=N-1,1,-1
       U(J)=U(J)-GAM(J+1)*U(J+1)
      END DO
      RETURN
      END SUBROUTINE TRIDIAG
