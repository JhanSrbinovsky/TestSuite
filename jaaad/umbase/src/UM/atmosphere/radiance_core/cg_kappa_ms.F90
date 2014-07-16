#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the product of CG and KAPPA.
!
! Purpose:
!   This routine calculates a sum of products of hemispheric
!   integrals and Clebsch-Gordan coefficients used in the
!   specification of BRDFs. These terms might be stored, but
!   this could involve the use of too much memory.
!
! Method:
!
!
!   Indexing is a bit complicated. The BRDF is truncated at an
!   even order and l'+m and l+m must both be even, so the effective
!   truncation when m is odd is one order lower. Nominally,
!   CGK is CGK(l',l) where l' takes alternate values though l takes
!   consecutive values, so we map into the actual array as
!       (l', l) --> ( (l'-m+2)/2, (l-m+1) )
!
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE CG_KAPPA_MS(MS, LS_TRUNC, LS_BRDF_TRUNC                &
     &  , CG_COEFF, KAPPA                                               &
     &  , CGK                                                           &
     &  , ND_MAX_ORDER, ND_BRDF_TRUNC                                   &
     &  )
!
!
!
      IMPLICIT NONE
!
! Include Header Files
#include "c_kinds.h"
!
!     Sizes of arrays
      INTEGER, INTENT(IN) ::                                            &
     &    ND_MAX_ORDER                                                  &
!           Size allocated for orders of spherical harmonics
     &  , ND_BRDF_TRUNC
!           Size allocated for orders of spherical harmonics
!           in BRDFs
!
!     Dummy arguments
      INTEGER, INTENT(IN) ::                                            &
     &    MS                                                            &
!           Azimuthal order
     &  , LS_TRUNC                                                      &
!           The order of truncation applied to spherical harmonics
     &  , LS_BRDF_TRUNC
!           The order of truncation applied to the BRDF
      REAL  (Real64), INTENT(IN) ::                                     &
     &    CG_COEFF(LS_TRUNC-MS+1)                                       &
!           Clebsch-Gordan coefficients
     &  , KAPPA(ND_MAX_ORDER/2, ND_MAX_ORDER/2)
!           Integrals of pairs of spherical harmonics over the downward
!           hemisphere
!
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    CGK(ND_BRDF_TRUNC/2+1, ND_MAX_ORDER)
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
!     Consider first the case where l+m is even. In this case the
!     documented formula is applied directly, with the omission
!     of a term in the first element of the array where l' would
!     be out of range.
      DO LSR=1, LS_TRUNC+1-MS, 2
        CGK(1, LSR)=CG_COEFF(1)*KAPPA(1, (LSR+1)/2)
        DO LSR_P=3, LS_BRDF_TRUNC-MS+1-MOD(MS, 2), 2
          CGK((LSR_P+1)/2, LSR)                                         &
     &      =CG_COEFF(LSR_P)*KAPPA((LSR_P+1)/2, (LSR+1)/2)              &
     &      +CG_COEFF(LSR_P-1)*KAPPA((LSR_P-1)/2, (LSR+1)/2)
        ENDDO
      ENDDO
!     In the case where l+m is odd the array is generally zero, so
!     we initialize all such values and calculate exceptional cases
!     later. Note that KAPPA does not appear in these loops because
!     in the compressed format these trivially evaluated values are
!     not stored.
      DO LSR=2, LS_TRUNC+1-MS, 2
        DO LSR_P=1, LS_BRDF_TRUNC-MS+1-MOD(MS, 2), 2
          CGK((LSR_P+1)/2, LSR)=0.0E+00_Real64
        ENDDO
      ENDDO
!     The case l=l'+1:
      DO LSR=2, LS_BRDF_TRUNC-MS-MOD(MS, 2)+2, 2
        CGK(LSR/2, LSR)=CG_COEFF(LSR-1)*0.5E+00_Real64
      ENDDO
!     The case l=l'-1:
      DO LSR=2, LS_BRDF_TRUNC-MS-MOD(MS, 2), 2
        CGK(LSR/2+1, LSR)=CG_COEFF(LSR)*0.5E+00_Real64
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE CG_KAPPA_MS
#endif
