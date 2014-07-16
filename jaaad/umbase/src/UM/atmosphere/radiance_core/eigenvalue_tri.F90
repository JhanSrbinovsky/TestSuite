#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to find the eigenvalues of a symmetric tridiagonal matrix.
!
! Purpose:
!   To caulate the eigenvalues of a symmetric tridiagonal matrix.
!
! Method:
!   The standard QR-algorithm with shifts is used, though this routine
!   is not a fully general implementation. The algorithm is based on the
!   pseudo-code and description given in "Numerical Analysis" by
!   R. L. Burden and D. J. Faires (PWS-Kent 1989).
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE EIGENVALUE_TRI(N_MATRIX, N_IN, D, E                    &
     &   , TOL, N_MAX_ITERATION                                         &
     &   , ND_MATRIX)
!
!
      IMPLICIT NONE
!
!     Sizes of arrays
      INTEGER, INTENT(IN) ::                                            &
     &   ND_MATRIX
!          Size allocated for matrices treated together.
!
!     Include header files.
#include "c_kinds.h"
#include "def_std_io_icf3z.h"
!
!     Dummy arguments
!
      INTEGER, INTENT(IN) ::                                            &
     &   N_MATRIX                                                       &
!          Number of matrices treated together
     & , N_IN                                                           &
!          Order of input matrix
     & , N_MAX_ITERATION
!          Maximum number of iterations
!
      REAL  (Real64), INTENT(IN) ::                                     &
     &   TOL
!          Tolerance for setting the subdiagonal elements to 0.
!
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &   D(ND_MATRIX, N_IN)                                             &
!          Main diagonal of the matrix: this will hold the eigenvalues
!          on output.
     & , E(ND_MATRIX, N_IN)
!          Subdiagonal of the matrix: E(1) is not used. E is reduced
!          to below the tolerance by the routine.
!
!
!     Local Variables:
!
      INTEGER                                                           &
     &   N                                                              &
!          Current working size of the problem
     & , L                                                              &
!          Loop variable
     & , J                                                              &
!          Loop variable
     & , ITERATION
!          Current iteration
      REAL  (Real64) ::                                                 &
     &   SHIFT(N_MATRIX)                                                &
!          Accumulated `shift'
     & , D_SHIFT(N_MATRIX)                                              &
!          Increment in `shift'
     & , B                                                              &
!          Temporary variable used in solving quadratic
     & , C                                                              &
!          Temporary variable used in solving quadratic
     & , DISCR                                                          &
!          Discriminant used in solving quadratic
     & , KAPPA_1                                                        &
!          First root of quadratic
     & , KAPPA_2
!          Second root of quadratic
      REAL  (Real64) ::                                                 &
     &   ABS_E                                                          &
!          Maximum absolute value of diagonal elements of the
!          current rows of the matrices
     & , SINR(N_MATRIX)                                                 &
!          Sine of current rotation
     & , COSR(N_MATRIX)                                                 &
!          Cosine of current rotations
     & , COSR_TEMP                                                      &
!          Temporary cosine
     & , SQ                                                             &
!          Temporary square root
     & , SUP_DIAG(N_MATRIX)                                             &
!          Element of first superdiagonal of matrix on the J-1st row
     & , SUP_DIAG_OLD(N_MATRIX)
!          Previous value of SUP_DIAG
!
!
!
!     The algorithm proceeds iteratively. The matrix supplied, A, is
!     decomposed as A=QR where Q is orthogonal and R is upper
!     triangular. A'=RQ is then formed and the process is repeated with
!     A'. This leads to a sequence of matrices which converge to one
!     with the eigenvalues down the diagonal.
!
!     Initialization:
!     Reduce the working size of the matrix if the off-diagonal
!     elements are small enough.
      N=N_IN
      ABS_E=0.0E+00_Real64
      DO L=1, N_MATRIX
        ABS_E=MAX(ABS_E, ABS(E(L, N)))
      ENDDO
      DO WHILE ( (N >  1).AND.(ABS_E <  TOL) )
        N=N-1
        DO L=1, N_MATRIX
          ABS_E=MAX(ABS_E, ABS(E(L, N)))
        ENDDO
      ENDDO
!
      ITERATION=0
      DO L=1, N_MATRIX
        SHIFT(L)=0.0E+00_Real64
      ENDDO
!
!
      DO WHILE ( (N >  1).AND.(ITERATION <  N_MAX_ITERATION) )
!
!
        ITERATION=ITERATION+1
!
!       Form an estimate of the first eigenvalue to be found by
!       finding the eigenvalues of the 2x2 matrix at the bottom
!       right-hand corner.
        DO L=1, N_MATRIX
          B=D(L, N-1)+D(L, N)
          C=D(L, N-1)*D(L, N)-E(L, N)*E(L, N)
          DISCR=SQRT(B*B-4.0E+00_Real64*C)
!         For reasons of conditioning we calculate the root of largest
!         magnitude and determine the other from the product of the
!         roots.
          KAPPA_1=0.5E+00_Real64*(B+SIGN(DISCR, B))
          KAPPA_2=C/KAPPA_1
!
!         Calculate the `shift' so as to accelerate convergence to the
!         last eigenvalue. A simple two-branch IF-test should be
!         amenable to vectorization if the vector CPU has a vector
!         mask register.
          IF ( ABS(KAPPA_1-D(L, N)) <                                   &
     &         ABS(KAPPA_2-D(L, N)) ) THEN
            D_SHIFT(L)=KAPPA_1
          ELSE
            D_SHIFT(L)=KAPPA_2
          ENDIF
          SHIFT(L)=SHIFT(L)+D_SHIFT(L)
        ENDDO
!
!       Shift the diagonal elements.
        DO J=1, N
          DO L=1, N_MATRIX
            D(L, J)=D(L, J)-D_SHIFT(L)
          ENDDO
        ENDDO
!
!
!       Form the QR-decompostion of the matrix by constructing
!       rotations to eliminate the sub-diagonal elements. COSR(J)
!       and SINR(J) are the cosine and sine of the rotations to
!       eliminate the element (J, J-1) of the input matrix: these
!       values specify the transpose of Q as we really construct
!       R=Qt.A by this procedure. The upper triangular matrix, R,
!       has two superdiagonals, but in practice only the first
!       is required. As the resulting matrix, RQ, will be a
!       symmetric tridaigonal matrix only its diagonal, D, and
!       the sub-diagonal, E, need be formed.
!
!       Inintialize:
        DO L=1, N_MATRIX
          SUP_DIAG(L)=E(L, 2)
          COSR(L)=1.0E+00_Real64
          SINR(L)=0.0E+00_Real64
        ENDDO
!
        DO J=2, N
!
          DO L=1, N_MATRIX
!
!           This block of code is a little opaque as the variables
!           SINR and COSR are re-used to avoid the need to declare
!           them explicitly as vectors. We form the rotation to
!           elminate E(J) and also calculate E(J-1) of the new matrix
!           RQ using SINR(J-1) for the last time. The new cosine of
!           the rotation must be stored because we still need
!           the old one.
            SQ=SQRT(D(L, J-1)*D(L, J-1)+E(L, J)*E(L, J))
            E(L, J-1)=SINR(L)*SQ
            SINR(L)=E(L, J)/SQ
            COSR_TEMP=D(L, J-1)/SQ
!
!           Adjust the superdiagonal of the previous row of the matrix
!           as required by the elimination. The calculation of D(J-1)
!           actually belongs to the formation of RQ, but is done here
!           before we overwrite COSR.
            SUP_DIAG_OLD(L)=SUP_DIAG(L)
            SUP_DIAG(L)=COSR_TEMP*SUP_DIAG(L)+SINR(L)*D(L, J)
            D(L, J-1)=COSR(L)*D(L, J-1)+SINR(L)*SUP_DIAG(L)
            COSR(L)=COSR_TEMP
!
!           Adjustments to the current row:
            D(L, J)=-SINR(L)*SUP_DIAG_OLD(L)+COSR(L)*D(L, J)
            IF (J <  N) SUP_DIAG(L)=COSR(L)*E(L, J+1)
!
          ENDDO
!
        ENDDO
!
        DO L=1, N_MATRIX
          E(L, N)=SINR(L)*D(L, N)
          D(L, N)=COSR(L)*D(L, N)
        ENDDO
!
!
!       Test for convergence and `shift' the converged eigenvalues.
!       back to their true values.
        ABS_E=0.0E+00_Real64
        DO L=1, N_MATRIX
          ABS_E=MAX(ABS_E, ABS(E(L, N)))
        ENDDO
        DO WHILE ( (N >  1).AND.(ABS_E <  TOL) )
          DO L=1, N_MATRIX
            D(L, N)=D(L, N)+SHIFT(L)
          ENDDO
          N=N-1
          DO L=1, N_MATRIX
            ABS_E=MAX(ABS_E, ABS(E(L, N)))
          ENDDO
        ENDDO
!
!
      ENDDO
!
!
!     Check that convergence has occurred.
      IF (N >  1) THEN
        WRITE(IU_ERR, '(/A)')                                           &
     &    '*** Warning: Convergence has not occurred while '            &
     &    //'calculating eigenvalues.'                                  &
     &    , 'The calculation continues.'
      ELSE
!       Shift the first eigenvalue back to its true value.
        DO L=1, N_MATRIX
          D(L, 1)=D(L, 1)+SHIFT(L)
        ENDDO
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE EIGENVALUE_TRI
#endif
