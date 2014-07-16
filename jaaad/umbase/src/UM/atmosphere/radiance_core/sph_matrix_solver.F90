#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to solve a set of stepped block equations.
!
! Method:
!       A set of linear equations of the stepped block form
!       is solved using Gaussian elimination with pivoting by rows
!       over a restricted range. In this application it should
!       not be necessary to consider all potential equations
!       when choosing pivots and this helps to reduce the
!       band-width of the system.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SPH_MATRIX_SOLVER(N_MATRIX, N_STEP, N_BLOCK            &
     &  , A, B                                                          &
     &  , X                                                             &
     &  , ND_MATRIX, ND_EQUATION, ND_DIAGONAL                           &
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
     &  , ND_EQUATION                                                   &
!           Size allocated for equations
     &  , ND_DIAGONAL
!           Size allocated for diagonals
!
!     Dummy arguments
      INTEGER, INTENT(IN) ::                                            &
     &    N_MATRIX                                                      &
!           Number of matrices
     &  , N_STEP                                                        &
!           Number of steps in the matrix
     &  , N_BLOCK
!           Size of each block
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    A(ND_MATRIX, ND_EQUATION, ND_DIAGONAL)                        &
!           Matrices of coefficients
     &  , B(ND_MATRIX, ND_EQUATION)
!           Righthand sides
      REAL  (Real64), INTENT(OUT) ::                                    &
     &     X(ND_MATRIX, ND_EQUATION)
!           Solution vector
!
!
!     Local variables
      INTEGER                                                           &
     &    I_STEP                                                        &
!           Counter for steps in the matrix
     &  , I_PHASE                                                       &
!           Counter for phase of elimination or back-substitution
     &  , IE                                                            &
!           Number of equation
     &  , IC                                                            &
!           Column of the compressed matrix
     &  , RIGHT                                                         &
!           Rightmost column of the compressed matrix
     &  , ROW_FIRST                                                     &
!           First row of the matrix considered in the second phase
!           of elimination
     &  , ROW_LAST
!           Last row of the matrix considered in the second phase
!           of elimination
      INTEGER                                                           &
     &    I_PIVOT(ND_MATRIX)                                            &
!           Index of pivot
     &  , OFFSET_PIVOT(ND_MATRIX)
!           Offset of the pivoting element relative to the current row
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , J                                                             &
!           Loop variable
     &  , K                                                             &
!           Loop variable
     &  , L
!           Loop variable
      REAL  (Real64) ::                                                 &
     &    PIVOT(ND_MATRIX)                                              &
!           Absolute value of the pivoting element
     &  , RHO(ND_MATRIX)                                                &
!           Scaling applied to the pivotal row in elimination
     &  , AABS                                                          &
!           Absolute value of the current potential pivot
     &  , TMP
!           Temporary variable used in swapping rows
!
!
!
!     Eliminative phase:
      I_STEP=1
      DO WHILE (I_STEP <= N_STEP)
!
!
!       The elminination falls naturally into two phases given
!       the structure of the matrix.
        DO I_PHASE=1, 2
!
          DO I=1, N_BLOCK
!
!           Set the number of the equation for elimination
!           and its column in the reduced matrix.
            IE=N_BLOCK*(2*I_STEP+I_PHASE-3)+I
            IC=N_BLOCK*(3-I_PHASE)+I
            IF (I_STEP <  N_STEP) THEN
              RIGHT=N_BLOCK*(8-2*I_PHASE)
            ELSE
              RIGHT=N_BLOCK*(6-2*I_PHASE)
            ENDIF
!
!           Choose the row for pivoting.
            DO L=1, N_MATRIX
              PIVOT(L)=ABS(A(L, IE, IC))
              I_PIVOT(L)=IE
              OFFSET_PIVOT(L)=0
            ENDDO
!           In the next line (and a little later), we limit the
!           range to the size of the matrix because in the last
!           layer there is no further boundary to consider.
            DO J=IE+1, MIN(IE+N_BLOCK*I_PHASE-I, 2*N_STEP*N_BLOCK)
              DO L=1, N_MATRIX
                AABS=ABS(A(L, J, IC))
                IF (AABS >  PIVOT(L)) THEN
                  PIVOT(L)=AABS
                  I_PIVOT(L)=J
                ENDIF
              ENDDO
            ENDDO
!           In the first phase we also need to consider rows
!           which will have a different offset as they involve
!           conditions on the next boundary below.
            IF (I_PHASE == 1) THEN
              DO J=N_BLOCK*(2*I_STEP-1)+1                               &
     &          , N_BLOCK*(MIN(2*I_STEP+1, 2*N_STEP))
                DO L=1, N_MATRIX
                  AABS=ABS(A(L, J, IC-2*N_BLOCK))
                  IF (AABS >  PIVOT(L)) THEN
                    PIVOT(L)=AABS
                    I_PIVOT(L)=J
                    OFFSET_PIVOT(L)=2*N_BLOCK
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
!
!           Swap the rows regardlessly to allow vectorization.
            DO J=IC, RIGHT
              DO L=1, N_MATRIX
                TMP=A(L, IE, J)
                A(L, IE, J)=A(L, I_PIVOT(L), J-OFFSET_PIVOT(L))
                A(L, I_PIVOT(L), J-OFFSET_PIVOT(L))=TMP
              ENDDO
            ENDDO
            DO L=1, N_MATRIX
              TMP=B(L, IE)
              B(L, IE)=B(L, I_PIVOT(L))
              B(L, I_PIVOT(L))=TMP
            ENDDO
!
!           Now eliminate. In both phases we have equations dealing
!           with the same boundary as the pivoting equation, but in
!           the first phase there is also an additional elimination
!           to be performed referring to the next boundary below.
            ROW_FIRST=IE+1
            IF (I_STEP <  N_STEP) THEN
              ROW_LAST=IE+N_BLOCK*I_PHASE-I
            ELSE
              ROW_LAST=IE+N_BLOCK-I
            ENDIF
            DO J=ROW_FIRST, ROW_LAST
              DO L=1, N_MATRIX
                RHO(L)=A(L, J, IC)/A(L, IE, IC)
                B(L, J)=B(L, J)-RHO(L)*B(L, IE)
              ENDDO
              DO K=IC+1, RIGHT
                DO L=1, N_MATRIX
                  A(L, J, K)=A(L, J, K)-RHO(L)*A(L, IE, K)
                ENDDO
              ENDDO
            ENDDO
!           This is the extra elimination required during the first
!           phase.
            IF (I_PHASE == 1) THEN
              ROW_FIRST=N_BLOCK*(2*I_STEP-1)+1
              IF (I_STEP <  N_STEP) THEN
                ROW_LAST=N_BLOCK*(2*I_STEP+1)
              ELSE
                ROW_LAST=N_BLOCK*2*I_STEP
              ENDIF
              DO J=ROW_FIRST, ROW_LAST
                DO L=1, N_MATRIX
                  RHO(L)=A(L, J, IC-2*N_BLOCK)/A(L, IE, IC)
                  B(L, J)=B(L, J)-RHO(L)*B(L, IE)
                ENDDO
                DO K=IC+1, RIGHT
                  DO L=1, N_MATRIX
                    A(L, J, K-2*N_BLOCK)                                &
     &                =A(L, J, K-2*N_BLOCK)-RHO(L)*A(L, IE, K)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
!
          ENDDO
        ENDDO
!
        I_STEP=I_STEP+1
!
      ENDDO
!
!
!
!     Back-subsititution:
      I_STEP=N_STEP
      DO WHILE (I_STEP >= 1)
!
!
        DO I_PHASE=2, 1, -1
!
          IF (I_STEP <  N_STEP) THEN
            RIGHT=N_BLOCK*(8-2*I_PHASE)
          ELSE
            RIGHT=N_BLOCK*(6-2*I_PHASE)
          ENDIF
!
          DO I=N_BLOCK, 1, -1
!
            IE=N_BLOCK*(2*I_STEP+I_PHASE-3)+I
            IC=(3-I_PHASE)*N_BLOCK+I
!
            DO L=1, N_MATRIX
              X(L, IE)=B(L, IE)
            ENDDO
            DO J=1, RIGHT-IC
              DO L=1, N_MATRIX
                X(L, IE)=X(L, IE)-A(L, IE, J+IC)*X(L, IE+J)
              ENDDO
            ENDDO
            DO L=1, N_MATRIX
              X(L, IE)=X(L, IE)/A(L, IE, IC)
            ENDDO
!
          ENDDO
        ENDDO
!
        I_STEP=I_STEP-1
!
!
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE SPH_MATRIX_SOLVER
#endif
