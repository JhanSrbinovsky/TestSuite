#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate hemispheric spherical integrals.
!
! Purpose:
!   This routine calculates hemispheric integrals of products
!   of spherical harmonics for a fixed azimuthal order for use
!   in Marshak's boundary conditions.
!
! Method:
!
!   We require the integral of Y_l'^m* Y_l^m over the downward
!   hemisphere for those l' such that l'+m is odd. If l=l' the
!   integral is trivially 1/2, but otherwise it will be zero
!   unless l+l' is odd. To reduce storage we omit the case l=l'
!   here and then map
!       (l', l) --> ( (l'-m+1)/2, (l-m+2)/2)
!   in the stored array.
!
!   The integrals are evaluated from the values of spherical
!   harmonics and their derivatives at a polar angle of pi/2.
!
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE HEMI_SPH_INTEG(LS_TRUNC, MS, UPLM_ZERO                 &
     &  , KAPPA                                                         &
     &  , ND_MAX_ORDER                                                  &
     &  )
!
!
      IMPLICIT NONE
!
!
!     Sizes of arrays
      INTEGER, INTENT(IN) ::                                            &
     &    ND_MAX_ORDER
!           Size allocated for orders of spherical harmonics
!
!     Include header files.
#include "c_kinds.h"
#include "c_pi.h"
!
!     Dummy arguments
      INTEGER, INTENT(IN) ::                                            &
     &    LS_TRUNC                                                      &
!           The truncating order of the system of equations
     &  , MS
!           Azimuthal order
      REAL  (Real64), INTENT(IN) ::                                     &
     &    UPLM_ZERO(LS_TRUNC+1-MS)
!           Spherical harmonics and derivatives at a polar angle of
!           pi/2
!
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    KAPPA(ND_MAX_ORDER/2, ND_MAX_ORDER/2)
!           Integrals of pairs of spherical harmonics over the downward
!           hemisphere
!
!     Local variables:
      INTEGER                                                           &
     &    LSR_P                                                         &
!           Reduced primed polar order
     &  , LSR
!           Reduced polar order
!
!
!
!     The outer loop is over l' where l'+m is odd. Indexing is done
!     using the reduced indices l'+1-m and l+1-m.
!
      DO LSR_P=2, LS_TRUNC+1-MS, 2
        DO LSR=1, LS_TRUNC-MS, 2
          KAPPA(LSR_P/2, (LSR+1)/2)=2.0E+00_Real64*PI                   &
     &      *REAL(1-2*MOD(LSR_P, 2), Real64)                            &
     &      *UPLM_ZERO(LSR)*UPLM_ZERO(LSR_P)                            &
     &      /REAL((LSR-LSR_P)*(LSR+LSR_P-1+2*MS), Real64)
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE HEMI_SPH_INTEG
#endif
