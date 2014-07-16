#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set up and solve the eigensystem.
!
! Purpose:
!   For a given value of the azimuthal quantum number, MS, this
!   routine returns the positive eigenvalues imposed by the trunctaion
!   in one layer and the corresponsing eigenvectors.
!
! Method:
!   The sub-diagonal of the full matrix is calculated and then reduced
!   to the diagonal and subdiagonal of the reduced matrix. The
!   eigenvalues are then found by calling the QR-algorithm and the
!   eigenvectors are obtained from a recurrence relation.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE EIG_SYS(N_PROFILE, LS_TRUNC, MS, N_RED_EIGENSYSTEM     &
     &  , CG_COEFF, SQS                                                 &
     &  , MU, EIG_VEC                                                   &
     &  , ND_PROFILE, ND_RED_EIGENSYSTEM, ND_MAX_ORDER                  &
     &  )
!
!
      IMPLICIT NONE
!
!
!     Sizes of arrays
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for atmospheric profiles
     &  , ND_RED_EIGENSYSTEM                                            &
!           Size allocated for the reduced eigensystem
     &  , ND_MAX_ORDER
!           Size allocated for the order of the calculation
!
!     Include header files
#include "c_kinds.h"
#include "sph_qr_iter_acf3z.h"
!
!     Dummy variables
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of atmospheric profiles
     &  , LS_TRUNC                                                      &
!           Order of L for truncation
     &  , MS                                                            &
!           Azimuthal quantum number
     &  , N_RED_EIGENSYSTEM
!           Size of the reduced eigenproblem
      REAL  (Real64), INTENT(IN) ::                                     &
     &    CG_COEFF(LS_TRUNC+1-MS)                                       &
!           Clebsch-Gordan coefficients
     &  , SQS(ND_PROFILE, 0: ND_MAX_ORDER)
!           Square roots of S-coefficients
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    MU(ND_PROFILE, ND_RED_EIGENSYSTEM)                            &
!           Eigenvalues
     &  , EIG_VEC(ND_PROFILE, 2*ND_RED_EIGENSYSTEM, ND_RED_EIGENSYSTEM)
!           Eigenvectors
!
!
!     Local variables
      INTEGER                                                           &
     &    J                                                             &
!           Loop variable
     &  , K                                                             &
!           Loop variable
     &  , L                                                             &
!           Loop variable
     &  , LS                                                            &
!           Order of spherical harmonic
     &  , N_MAX_QR_ITERATION
!           Maximum number of QR iterations
      REAL  (Real64) ::                                                 &
     &    TOL                                                           &
!           Tolerance for assessing convergence of the QR-algorithm
     &  , EC(ND_PROFILE, 2*N_RED_EIGENSYSTEM)                           &
!           Sub-diagonal of full matrix
     &  , E(ND_PROFILE, N_RED_EIGENSYSTEM)                              &
!           Sub-diagonal of reduced matrix
     &  , NORMALIZATION(ND_PROFILE)                                     &
!           Normalization factor for the eigenvector
     &  , C_RATIO(ND_PROFILE)                                           &
!           Common ratio in the geometric progression used to rescale
!           the recurrence relation to avoid overflows
     &  , RESCALE(ND_PROFILE)
!           Multiplier to convert the terms of the scaled recurrence
!           to the the final quantities
!
!
!
!     Set the tolerance for convergence of the algorithm from the
!     precision of the machine.
      TOL=RP_TOL_FACTOR_SPH_QR*EPSILON(RP_TOL_FACTOR_SPH_QR)
!
!     Calculate the reduced matrix to yield the eigenvalues. EC_...
!     represent elements of the sub-diagonal of the full matrix:
!     D and E are the diagonal and sub-diagonal of the reduced matrix.
!
!
!     Calculate the sub-diagonal of the full matrix.
      DO J=2, LS_TRUNC+1-MS
        LS=MS-1+J
        DO L=1, N_PROFILE
          EC(L, J)=CG_COEFF(J-1)/(SQS(L, LS)*SQS(L, LS-1))
        ENDDO
      ENDDO

!
!     Retain odd rows and columns of the square of the preceeding
!     matrix. The diagonal terms are stored in MU as this will be
!     reduced to the eigenvalues later.
      DO L=1, N_PROFILE
        MU(L, 1)=EC(L, 2)**2
      ENDDO
      DO J=2, N_RED_EIGENSYSTEM
        DO L=1, N_PROFILE
          MU(L, J)=EC(L, 2*J-1)**2+EC(L, 2*J)**2
          E(L, J)=EC(L, 2*J-2)*EC(L, 2*J-1)
        ENDDO
      ENDDO
!
!     Determine the eigenvalues of the reduced matrix, which
!     are the squares of the (positive) eigenvalues of the
!     full matrix. If the eigensystem is of size 1 no calculation
!     is required.
      IF (N_RED_EIGENSYSTEM >  1) THEN
!       The number of iterations required for convergence increases
!       as the order of truncation rises. A small allowance is made
!       for extra iterations.
        N_MAX_QR_ITERATION=LS_TRUNC+25
! DEPENDS ON: eigenvalue_tri
        CALL EIGENVALUE_TRI(N_PROFILE, N_RED_EIGENSYSTEM                &
     &    , MU, E, TOL, N_MAX_QR_ITERATION                              &
     &    , ND_PROFILE                                                  &
     &    )
      ENDIF
      DO K=1, N_RED_EIGENSYSTEM
        DO L=1, N_PROFILE
          MU(L, K)=SQRT(MU(L, K))
          IF (MU(L, K) >  1.0E+00_Real64) THEN
            C_RATIO(L)=5.0E-01_Real64/MU(L, K)
          ELSE
            C_RATIO(L)=1.0E+00_Real64
          ENDIF
        ENDDO
!
!       Use the recurrence relation to find the eigenvectors of the
!       full matrix. For large values of MU there will be an
!       eigenvector like MU^J and one like MU^-J. The latter (minimal)
!       solution is required, but for |MU|>1 the recurrence is
!       unstable, so the growing solution will will swamp the required
!       solution. Conversely, with downward recurrence, the desired
!       solution grows and will dominate in the recurrence. When
!       |MU|<1 the recurrence is stable in either direction so downward
!       recurrence is used consistently. On further complication must
!       be taken into account: if MU is very large (as can happen with
!       almost conservative scattering) the elements of the eigenvector
!       may be of so large a range of magnitudes that the recurrence
!       overflows. A scaling factor, c, is therefore introduced so that
!       the j'th element of the eigenvector, e_j=c^j.e_j'. s may not
!       be less than 1 for small eigenvalues or the same problem will
!       be introduced with them; the vector e' has elements of order 1.
!
        J=2*N_RED_EIGENSYSTEM
        DO L=1, N_PROFILE
          EIG_VEC(L, J, K)=1.0E+00_Real64
        ENDDO
        J=J-1
        DO L=1, N_PROFILE
          EIG_VEC(L, J, K)=C_RATIO(L)*MU(L, K)/EC(L, J+1)
        ENDDO
        DO WHILE(J >  1)
          J=J-1
          DO L=1, N_PROFILE
            EIG_VEC(L, J, K)                                            &
     &        =(MU(L, K)*EIG_VEC(L, J+1, K)                             &
     &        -C_RATIO(L)*EC(L, J+2)*EIG_VEC(L, J+2, K))                &
     &        *C_RATIO(L)/EC(L, J+1)
          ENDDO
        ENDDO
!
!       Remove the scaling factor, renormalize the eigenvector
!       and rescale by the s-coefficients for later efficiency.
        DO L=1, N_PROFILE
          RESCALE(L)=C_RATIO(L)
          EIG_VEC(L, 1, K)=EIG_VEC(L, 1, K)*RESCALE(L)
          NORMALIZATION(L)=EIG_VEC(L, 1, K)*EIG_VEC(L, 1, K)
        ENDDO
        DO J=2, 2*N_RED_EIGENSYSTEM
          DO L=1, N_PROFILE
            RESCALE(L)=RESCALE(L)*C_RATIO(L)
            EIG_VEC(L, J, K)=EIG_VEC(L, J, K)*RESCALE(L)
            NORMALIZATION(L)=NORMALIZATION(L)                           &
     &        +EIG_VEC(L, J, K)*EIG_VEC(L, J, K)
          ENDDO
        ENDDO
        DO L=1, N_PROFILE
          NORMALIZATION(L)=SQRT(1.0E+00_Real64/NORMALIZATION(L))
        ENDDO
        DO J=1, 2*N_RED_EIGENSYSTEM
          DO L=1, N_PROFILE
            EIG_VEC(L, J, K)=EIG_VEC(L, J, K)*NORMALIZATION(L)          &
     &        /SQS(L, J+MS-1)
          ENDDO
        ENDDO
!
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE EIG_SYS
#endif
