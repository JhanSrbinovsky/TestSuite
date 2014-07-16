#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate spherical Clebsch-Gordan coefficients.
!
! Purpose:
!   The routine yields the Clebsch-Gordan coefficients between the
!   spherical harmonics Y_l^m and Y_1^0, c_{lm}^+ in the notation of
!   the description of the algorithm, or <l+1,m|1,0,l,m> in standard
!   notation. These are stored in one array with addressing determined
!   by the truncation.
!
! Method:
!   As this routine is called only once speed is not too critical
!   so direct calculation is used. Only values for m>0 are required.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE CALC_CG_COEFF(LS_MAX_ORDER                             &
     &  , IA_SPH_MM, MS_MIN, MS_TRUNC                                   &
     &  , CG_COEFF                                                      &
     &  , ND_MAX_ORDER, ND_SPH_COEFF)
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
     &    ND_MAX_ORDER                                                  &
!           Size allocated for orders of spherical harmonics
     &  , ND_SPH_COEFF
!           Size of array of spherical coefficients
      INTEGER, INTENT(IN) ::                                            &
     &    LS_MAX_ORDER                                                  &
!           Maximum order of harmonics required
     &  , MS_MIN                                                        &
!           Lowest azimuthal order calculated
     &  , MS_TRUNC(0: ND_MAX_ORDER)                                     &
!           Truncation in MS for this order
     &  , IA_SPH_MM(0: ND_MAX_ORDER)
!           Position of Clebsh-Gordan coefficient with m=0 for the
!           given value of l.
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    CG_COEFF(ND_SPH_COEFF)
!
!
!     Local variables
      INTEGER                                                           &
     &    LS                                                            &
!          Order of harmonic
     &  , MS
!           Azimuthal order of harmonic
      REAL  (Real64) ::                                                 &
     &    INV
!           l-dependent denominator
!
!
!
      DO LS=0, LS_MAX_ORDER
        INV=1.0E+00_Real64/REAL((2*LS+1)*(2*LS+3), Real64)
        DO MS=MS_MIN, MS_TRUNC(LS)
          CG_COEFF(IA_SPH_MM(MS)+LS-MS)                                 &
     &      =SQRT(REAL((LS+1-MS)*(LS+1+MS), Real64)*INV)
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE CALC_CG_COEFF
#endif
