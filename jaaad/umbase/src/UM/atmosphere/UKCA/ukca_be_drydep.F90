#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  To assign dry deposition rates in s-1 to array odryrt
!
!           Called from UKCA_CHEMISTRY_CTL
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Current Code Owner:       Colin Johnson/Fiona O'Connor
!
!  Code Description:
!   Language:  FORTRAN 90 (formatted)
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
      SUBROUTINE UKCA_BE_DRYDEP(k, n_pnts, nlev_in_bl2, dryrt, odryrt)

      USE ASAD_MOD,          ONLY:  ndepd, nldepd
      IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "parvars.h"
#include "typsize.h"

      INTEGER, INTENT(IN) :: k                  ! Level number
      INTEGER, INTENT(IN) :: n_pnts             ! No of spatial points
      INTEGER, INTENT(IN) :: nlev_in_bl2(theta_field_size) ! No levs in BL

      REAL, INTENT(IN) :: dryrt(theta_field_size,jpdd)   ! Dry dep rates in s-1
      REAL, INTENT(OUT):: odryrt(theta_field_size,jpspec)! Dry dep rates in s-1

!     Local variables

      INTEGER :: i                              ! Loop variable
      INTEGER :: js                             ! Loop variable
      INTEGER :: nspec                          ! Pointer for species

      DO js = 1, jpspec
        DO i = 1, theta_field_size
          odryrt(i,js) = 0.0
        END DO
      END DO

      DO i = 1,n_pnts
        IF (k <= nlev_in_bl2(i)) THEN
          DO js = 1,ndepd
            nspec = nldepd(js)
            odryrt(i,nspec) = dryrt(i,js)
          END DO
        END IF
      END DO

      RETURN
      END SUBROUTINE UKCA_BE_DRYDEP
#endif
