#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Main routine for TDF scheme.



!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
!+ Remove problems and inconsistencies from cloud variables.

      SUBROUTINE ModifyCloudVars ( qCL,                                 &
                                                           ! inout
     &                             qCF,                                 &
                                                           ! inout
     &                             area_cloud_fraction,                 &
                                                           ! inout
     &                             bulk_cloud_fraction,                 &
                                                           ! inout
     &                             cloud_fraction_liquid,               &
                                                           ! inout
     &                             cloud_fraction_frozen ) ! inout

! Description:
!
!   Remove problems and inconsistencies from cloud variables to
!   prevent crashes on subsequent calls to certain physics routines.
!
!   1. Remove negative qCLs and qCFs.
!   2. Zero cloud fields where qCL and/or qCF are zero.
!
!
! Current Code Owner: Adam Clayton.
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.2   13/10/00   Original code replacing old routine ZeroCloud.
!                    Adam Clayton
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! Declarations:

      IMPLICIT NONE

! Common blocks:

#include "parvars.h"
#include "typsize.h"
#include "ctfilt.h"

! Subroutine arguments:

      REAL, INTENT(INOUT) ::                                            &
      
     &  qCL ( 1 - halo_i : row_length + halo_i,                         &
     &        1 - halo_j : rows       + halo_j,                         &
     &        1          : wet_levels ),                                &
      
     &  qCF ( 1 - halo_i : row_length + halo_i,                         &
     &        1 - halo_j : rows       + halo_j,                         &
     &        1          : wet_levels ),                                &
      
     &  area_cloud_fraction   ( row_length, rows, wet_levels ),         &
     &  bulk_cloud_fraction   ( row_length, rows, wet_levels ),         &
     &  cloud_fraction_liquid ( row_length, rows, wet_levels ),         &
     &  cloud_fraction_frozen ( row_length, rows, wet_levels )

! Local variables:

      INTEGER :: i, j, k

!- End of header ------------------------------------------------------

      DO k = 1, wet_levels
        DO j = 1, rows
          DO i = 1, row_length

            IF ( qCL(i,j,k) < 0.0 ) qCL(i,j,k) = 0.0
            IF ( qCF(i,j,k) < 0.0 ) qCF(i,j,k) = 0.0

            IF ( qCL(i,j,k) <= q_CC_tol )                               &
     &        cloud_fraction_liquid(i,j,k) = 0.0

            IF ( qCF(i,j,k) <= q_CC_tol )                               &
     &        cloud_fraction_frozen(i,j,k) = 0.0

            IF ( qCL(i,j,k) <= q_CC_tol .AND.                           &
     &           qCF(i,j,k) <= q_CC_tol ) THEN
              area_cloud_fraction(i,j,k) = 0.0
              bulk_cloud_fraction(i,j,k) = 0.0
            END IF

          END DO
        END DO
      END DO


      RETURN
      END SUBROUTINE ModifyCloudVars
#endif
