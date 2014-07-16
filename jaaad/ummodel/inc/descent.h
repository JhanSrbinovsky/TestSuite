#if defined(A19_2A)
! DESCENT start

! Number of TRIFFID iterations for gradient descent to equilibrium.
      INTEGER,PARAMETER:: ITER_EQ = 10

! Minimum value for the denominator of the update equation. Ensures
! that gradient descent does not lead to an unstable solution.
      REAL,PARAMETER:: DENOM_MIN=1.0E-6

! Inverse timestep for gradient  descent to equilibrium (/360days).
      REAL,PARAMETER:: GAMMA_EQ = 1.0E-1

! DESCENT emd
#endif
