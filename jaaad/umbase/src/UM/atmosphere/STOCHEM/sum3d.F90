#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE SUM3D(x)
!-----------------------------------------------------------------------
!-
!-   Purpose and Methods : Sums 3D data over all processors.
!-
!-   Inputs  : X        - Input array with restricted latitudes.
!-             DBOUNDS  - Lat and Long indicies for each processor.
!-
!-   Outputs : GLOBAL   - 3-D global output array (exchanged with a
!-                        Common block).
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  5.0    14/05/99  Created.  C.E. Johnson
!  5.5    30/03/04  Removed common block SUM_3D. Now does summation in
!                   one step. K. Ketelsen
!  6.1    21/10/04  No change.
!  6.2    21/10/05  Replaced GSYNC with SSYNC. P.Selwood
!-
!VVV  V2.2  SUM3D 20/X/99
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!----------------------------------------------------------------------

      REAL, DIMENSION(nlnpe,nlpe,nlev), INTENT(IN) :: x

      INTEGER :: info
      REAL, DIMENSION(nlong,mnlat,nlev) :: global

      COMMON /SUM3_D/ global

      global = 0.0
      global(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1,:) = x(:,:,:)

      CALL GC_SSYNC(nproc,info)
      CALL GC_RSUM(nlong*mnlat*nlev,nproc,info,global)
      CALL GC_SSYNC(nproc,info)
      END SUBROUTINE SUM3D
#endif
