#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to solve a set of banded matrix equations.
!
! Method:
!       A set of bands matrix equations is solved using the
!       standard method of Gaussian elimination. Diagonals are
!       numbered downward (i.e. upper diagonals first).
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE BAND_SOLVER(N_MATRIX, N_EQUATION                       &
     &  , IU, IL                                                        &
     &  , A, B                                                          &
     &  , X                                                             &
     &  , RHO                                                           &
     &  , ND_MATRIX, ND_DIAGONAL, ND_EQUATION                           &
     &  )
!
!
!
      IMPLICIT NONE
!
! Include header files
#include "c_kinds.h"
!
!     Sizes of arrays
      INTEGER, INTENT(IN) ::                                            &
     &    ND_MATRIX                                                     &
!           Size alloacted for matrices
     &  , ND_DIAGONAL                                                   &
!           Size allocated for diagonals
     &  , ND_EQUATION
!           Size allocated for equations
!
!     Dummy arguments
      INTEGER, INTENT(IN) ::                                            &
     &    N_MATRIX                                                      &
!           Number of matrices
     &  , N_EQUATION                                                    &
!           Number of equations
     &  , IU                                                            &
!           Number of superdiagonals
     &  , IL
!           Number of subdiagonals
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    A(ND_MATRIX, ND_DIAGONAL, ND_EQUATION)                        &
!           Matrices of coefficients
     &  , B(ND_MATRIX, ND_EQUATION)
!           Righthand sides
      REAL  (Real64), INTENT(OUT) ::                                    &
     &     X(ND_MATRIX, ND_EQUATION)
!           Solution vector
      REAL  (Real64) ::                                                 &
                         !, INTENT(WORK)
     &     RHO(ND_MATRIX)
!           Temporary array
!
!     Local variables
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , J                                                             &
!           Loop variable
     &  , K                                                             &
!           Loop variable
     &  , L                                                             &
!           Loop variable
     &  , IU1
!           Local scalar
!
!
      IU1=IU+1
!     Eliminative phase.
      DO I=N_EQUATION, 2, -1
        DO J=1, MIN(IU, I-1)
          DO L=1, N_MATRIX
            RHO(L)=A(L, IU1-J, I-J)/A(L, IU1, I)
            B(L, I-J)=B(L, I-J)-RHO(L)*B(L, I)
          ENDDO
          DO K=1, MIN(IL, I-1)
            DO L=1, N_MATRIX
              A(L, IU1+K-J, I-J)=A(L, IU1+K-J, I-J)                     &
     &          -RHO(L)*A(L, IU1+K, I)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     Solution and back-substitution:
!
      IF ( (IU == 2).AND.(IL == 2) ) THEN
!       A special version is used for the pentadiagonal case to allow
!       us to chain operations together for efficiency on the CRAY
!       vector machines, as this particular case arises quite often.
!
!       First equation:
        DO L=1, N_MATRIX
          X(L, 1)=B(L, 1)/A(L, 3, 1)
        ENDDO
!       Second equation:
        DO L=1, N_MATRIX
          X(L, 2)=(B(L, 2)-A(L, 4, 2)*X(L, 1))/A(L, 3, 2)
        ENDDO
!       Remaining equations:
        DO I=3, N_EQUATION
          DO L=1, N_MATRIX
            X(L, I)=(B(L, I)-A(L, 4, I)*X(L, I-1)                       &
     &        -A(L, 5, I)*X(L, I-2))/A(L, 3, I)
          ENDDO
        ENDDO
      ELSE
!
!       General case:
        DO I=1, N_EQUATION
          DO L=1, N_MATRIX
               X(L, I)=B(L, I)
          ENDDO
          DO K=1, MIN(IL, I-1)
            DO L=1, N_MATRIX
              X(L, I)=X(L, I)-A(L, IU1+K, I)*X(L, I-K)
            ENDDO
          ENDDO
          DO L=1, N_MATRIX
            X(L, I)=X(L, I)/A(L, IU1, I)
          ENDDO
        ENDDO
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE BAND_SOLVER
#endif
