#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      FUNCTION GAUSS_RAN(seed,ransize,nfill)
! ----------------------------------------------------------------------
! Purpose:
! Returns normally distributed random numbers with
! zero mean and unit variance. See Numerical Recipes.
!
! Method:  Box-Muller
!
! Original Programmer: Colin Johnson
!
! Current code owner: Colin Johnson
!
!
! Current Owner of Code: C.E. Johnson
!
! History:
! Version   Date                    Comment
!  5.3    26/07/01  Created.  C.E. Johnson
!  5.3    23/11/01  Now returns array of numbers.  C.E. Johnson
!  5.5    23/02/04  Vectorised version. K. Ketelsen
!  6.1    20/10/04  Now uses integer random number seeds.
!
!VVV  V1.1  GAUSS_RAN 23/11/01 - Array output
! ----------------------------------------------------------------------
      IMPLICIT NONE
! ----------------------------------------------------------------------

      INTEGER, INTENT(IN) :: ransize
      INTEGER, INTENT(IN) :: nfill
      INTEGER, DIMENSION(ransize), INTENT(INOUT) :: seed

      REAL, DIMENSION(nfill) :: GAUSS_RAN

      INTEGER :: i, j
      INTEGER :: anz

      REAL, DIMENSION(nfill+1) :: gset
      REAL :: fac
      REAL :: rsq
      REAL :: v1
      REAL :: v2
      REAL, DIMENSION(2)        :: xran
      REAL, DIMENSION(nfill+1)  :: xran_a
      INTEGER, DIMENSION(nfill) :: ind_g

      j = 0
      DO i=1, nfill, 2
        j = j + 1
        ind_g(j) = i
      END DO
      anz = j

      DO
        DO j = 1, 2*anz
          xran_a(j) = FRANV(seed)
        END DO
!CDIR nodep
        DO i=1,anz
          xran(1) = xran_a(i)
          xran(2) = xran_a(i+1)
          v1 = 2.0 * xran(1) - 1.0
          v2 = 2.0 * xran(2) - 1.0
          rsq = v1**2 + v2**2
          IF (rsq < 1.0 .AND. rsq /= 0.0) THEN
            fac = sqrt(-2.0 * LOG(rsq) / rsq)
            gset(ind_g(i)) = v2 * fac
            gset(ind_g(i)+1) = v1 * fac
            ind_g(i) = -1
          END IF
        END DO
        j = 0
!CDIR nodep
        DO i=1,anz
          IF (ind_g(i) > 0) THEN
            j = j + 1
            ind_g(j)= ind_g(i)
          END IF
        END DO
        anz = j
        IF (anz == 0) EXIT
      END DO
      DO i=1, nfill
        GAUSS_RAN(i) = gset(i)
      END DO

      END FUNCTION GAUSS_RAN
#endif
