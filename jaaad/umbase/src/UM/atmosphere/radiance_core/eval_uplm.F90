#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate spherical harmonics excluding expoential.
!
! Purpose:
!   Spherical harmonics, Upsilon_lm, are calculated for given directions
!   for all values of l at a fixed value of m.
!
! Method:
!   Y_mm is known so upward recurrence on l is used.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE EVAL_UPLM(MS, N_MAX_ORDER, N_DIRECTION, X              &
     &   , UP_LM                                                        &
     &   , ND_DIRECTION                                                 &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     Sizes of dummy arrays
      INTEGER, INTENT(IN) ::                                            &
     &     ND_DIRECTION
!             Maximum number of directions
!
!     Include header files.
#include "c_kinds.h"
#include "c_pi.h"
!
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    MS                                                            &
!           Azimuthal quantum number of spherical harmonic
     &  , N_MAX_ORDER
!           Maximum order of harmonics to calculate
      INTEGER, INTENT(IN) ::                                            &
     &    N_DIRECTION
!           Number of directions
      REAL  (Real64), INTENT(IN) ::                                     &
     &    X(ND_DIRECTION)
!           Cosines of polar angels of viewing directions
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    UP_LM(ND_DIRECTION, N_MAX_ORDER+1-MS)
!           Non-azimuthal parts of spherical harmonics
!
!
!     Local variables
      INTEGER                                                           &
     &    LS                                                            &
!           Loop variable
     &  , J                                                             &
!           Loop variable
     &  , K
!           Loop variable
      REAL  (Real64) ::                                                 &
     &    L                                                             &
!           Principal quantum number of harmonic
     &  , M                                                             &
!           Azimuthal quantum number of harmonic
     &  , PRODUCT
!           Factorial terms in Y_lm
!
!
!
!     Start the recurrence for Y_mm.
      PRODUCT=1.0E+00_Real64
      M=REAL(MS, Real64)
      IF (MS >  0) THEN
        DO J=1, MS
          PRODUCT=(1.0E+00_Real64-5.0E-01_Real64/REAL(J, Real64))       &
     &                                                     *PRODUCT
        ENDDO
        DO K=1, N_DIRECTION
          UP_LM(K, 1)=(-1.0E+00_Real64)**MS                             &
     &      *SQRT((1.0E+00_Real64-X(K)*X(K))**MS*PRODUCT                &
     &      *(2.0E+00_Real64*M+1.0E+00_Real64)/(4.0E+00_Real64*PI))
        ENDDO
      ELSE
        DO K=1, N_DIRECTION
          UP_LM(K, 1)=1.0E+00_Real64/SQRT(4.0E+00_Real64*PI)
        ENDDO
      ENDIF
!
!
!     Calculate Y_(m+1),m if it is within bounds.
      IF (MS <  N_MAX_ORDER) THEN
        DO K=1, N_DIRECTION
          UP_LM(K, 2)=X(K)*SQRT(2.0E+00_Real64*M+3.0E+00_Real64)        &
     &      *UP_LM(K, 1)
        ENDDO
      ENDIF
!
!
!     Complete the recurrence on l.
      DO LS=MS+2, N_MAX_ORDER
        L=REAL(LS, Real64)
        DO K=1, N_DIRECTION
          UP_LM(K, LS+1-MS)=X(K)                                        &
     &      *SQRT(((2.0E+00_Real64*L-1.0E+00_Real64)                    &
     &      *(2.0E+00_Real64*L+1.0E+00_Real64))                         &
     &      /((L+M)*(L-M)))*UP_LM(K, LS-MS)                             &
     &      -SQRT(((2.0E+00_Real64*L+1.0E+00_Real64)                    &
     &      *(L-1.0E+00_Real64-M)*(L-1.0E+00_Real64+M))                 &
     &      /((2.0E+00_Real64*L-3.0E+00_Real64)*(L+M)*(L-M)))           &
     &      *UP_LM(K, LS-1-MS)
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE EVAL_UPLM
#endif
