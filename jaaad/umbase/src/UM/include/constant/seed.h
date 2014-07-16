#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C) \
 || defined(A19_1A) || defined(A19_2A)
! Start seed
! Description:
!   This file sets the values of the variables FRAC_MIN and FRAC_SEED
!
! Current Code Owner:
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   5.3   25/09/01  Portability changes.  Z. Gardner
!   5.5   17/04/03  Remove reference to obsolete section
!                   A03_7A. T.White
!
      ! Minimum areal fraction for PFTs.
      REAL, PARAMETER:: FRAC_MIN  = 1.0E-6

      ! "Seed" fraction for PFTs.
      REAL, PARAMETER:: FRAC_SEED = 0.01
#endif
! End Seed
