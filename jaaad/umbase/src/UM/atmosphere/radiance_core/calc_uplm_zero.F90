#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate Upsilon_l^m(0) or its derivative.
!
! Purpose:
!   This routine is called to determine the value of a spherical
!   harmonic with theta=pi/2 and phi=0 or the derivative for
!   alternate orders. This minimizes storage.
!
! Method:
!   As this routine is called only once speed is not too critical
!   so direct calculation is used.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE CALC_UPLM_ZERO(MS_MIN, MS_MAX, IA_SPH_MM               &
     &  , LS_LOCAL_TRUNC, UPLM_ZERO                                     &
     &  , ND_MAX_ORDER, ND_SPH_COEFF)
!
!
      IMPLICIT NONE
!
!
!     Include header files
#include "c_kinds.h"
#include "c_pi.h"
!
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_SPH_COEFF                                                  &
!           Number of spherical coefficients
     &  , ND_MAX_ORDER
!           Maximum order of calculation
      INTEGER, INTENT(IN) ::                                            &
     &    MS_MIN                                                        &
!           Lowest azimuthal order calculated
     &  , MS_MAX                                                        &
!           Highest azimuthal order calculated
     &  , IA_SPH_MM(0: ND_MAX_ORDER)                                    &
!           Address of spherical coefficient for (m, m) for each m
     &  , LS_LOCAL_TRUNC(0: ND_MAX_ORDER)
!           Local truncation at this order
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    UPLM_ZERO(ND_SPH_COEFF)
!           Array of Upsilon_l^m and derivatives at polar angles of pi/2
!
!
!     Local variables
      INTEGER                                                           &
     &    LS                                                            &
!           Order of harmonic
     &  , MS                                                            &
!           Azimuthal quantum number of harmonic
     &  , J                                                             &
!           Temporary address
     &  , K
!           Loop variable
!
!
!
      DO MS=MS_MIN, MS_MAX
!
!       Calculate Upsilon_m^m(0) to start the recurrence.
        J=IA_SPH_MM(MS)
        UPLM_ZERO(J)=1.0E+00_Real64/(4.0E+00_Real64*PI)
        DO K=3, 2*MS+1, 2
           UPLM_ZERO(J)=UPLM_ZERO(J)*REAL(K, Real64)/REAL(K-1, Real64)
        ENDDO
        UPLM_ZERO(J)=REAL(1-2*MOD(MS, 2), Real64)*SQRT(UPLM_ZERO(J))
!
!
!       Calculate DUpsilon_{m+1}^m(0) to start the recurrence for the
!       derivatives.
        J=J+1
        UPLM_ZERO(J)=3.0E+00_Real64/(4.0E+00_Real64*PI)
        DO K=5, 2*MS+3, 2
           UPLM_ZERO(J)=UPLM_ZERO(J)*REAL(K, Real64)/REAL(K-3, Real64)
        ENDDO
        UPLM_ZERO(J)=REAL(1-2*MOD(MS, 2), Real64)*SQRT(UPLM_ZERO(J))
!
!       Now apply the recurrence formulae:
        DO LS=MS+2, LS_LOCAL_TRUNC(MS)-1, 2
          J=IA_SPH_MM(MS)+LS-MS
!
!         Recurrence for Upsilon_l^m.
          UPLM_ZERO(J)=-UPLM_ZERO(J-2)                                  &
     &      *SQRT(REAL((2*LS+1)*(LS+MS-1)*(LS-MS-1), Real64)            &
     &      /REAL((2*LS-3)*(LS+MS)*(LS-MS), Real64))
!
!         Recurrence for the derivative Upsilon_(l+1)^m.
          UPLM_ZERO(J+1)=-UPLM_ZERO(J-1)                                &
     &      *SQRT(REAL((2*LS+3)*(LS+1+MS)*(LS+1-MS), Real64)            &
     &      /REAL((2*LS-1)*(LS+MS)*(LS-MS), Real64))
!
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE CALC_UPLM_ZERO
#endif
