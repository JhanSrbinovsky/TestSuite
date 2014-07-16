#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description: 
!  Subroutine to assign wet deposition rates in s-1 to
!  array owetrt.
!
!  Called from UKCA routine UKCA_CHEMISTRY_CTL.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Current code owner: Colin Johnson/Fiona O'Connor
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
      SUBROUTINE UKCA_BE_WETDEP(n_pnts, wetrt, owetrt)

      USE ASAD_MOD,      ONLY: ndepw, nldepw 
      IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "parvars.h"
#include "typsize.h"

        INTEGER, INTENT(IN) :: n_pnts                      ! No of spatial pts

        REAL, INTENT(IN) :: wetrt(theta_field_size,jpdw)   ! Wet dep rates (s-1)

        REAL, INTENT(OUT):: owetrt(theta_field_size,jpspec)! Wet dep rates (s-1)

!       Local variables

        INTEGER :: i                         ! Loop variable
        INTEGER :: js                        ! Loop variable
        INTEGER :: nspec                     ! Pointer for species

        DO js = 1, jpspec
          DO i = 1, theta_field_size
            owetrt(i,js) = 0.0
          END DO
        END DO

        DO js = 1,ndepw
          nspec = nldepw(js)
          DO i = 1,n_pnts
            owetrt(i,nspec)= wetrt(i,js)
          ENDDO
        ENDDO

        RETURN
        END SUBROUTINE UKCA_BE_WETDEP
#endif
