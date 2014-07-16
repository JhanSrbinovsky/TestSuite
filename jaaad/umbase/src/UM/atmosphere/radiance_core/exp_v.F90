#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to mimic a vector exponential function.
!
! Method:
!       The normal exponential function is called on an array
!       of input values. This is provided for systems where
!       there is no intrinsic EXP_V.
!
! Current owner of code: James Manners
!
! Description of code:
!   FORTRAN 77 with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE EXP_V(N, A, B)
!
!
!
!
      IMPLICIT NONE
!
! Include header files
#include "c_kinds.h"
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    N
!           Number of values to exponentiate
      REAL  (Real64), INTENT(IN) ::                                     &
     &    A(N)
!           Arguments to exponentials
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    B(N)
!           Output arguments
!
!
!     Local variables:
      INTEGER                                                           &
     &    I
!           Loop variable
!
!
!
      DO I=1, N
        B(I)=EXP(A(I))
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE EXP_V
#endif
