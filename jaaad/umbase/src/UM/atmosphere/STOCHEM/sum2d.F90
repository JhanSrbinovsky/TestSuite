#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE SUM2D(x)
!-----------------------------------------------------------------------
!-
!-   Purpose and Methods : Sums 2D data over all processors.
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
!  5.4    23/08/02  Created.  W.J. Collins
!  6.1    29/10/04  Removed incorrect comments and common block SUM_3D
!                   as not required. M.G. Sanderson
!  6.2    21/10/05  Replace GSYNC with SSYNC. P.Selwood
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!----------------------------------------------------------------------

      REAL, DIMENSION(nlnpe,nlpe),   INTENT(IN)   :: x

      INTEGER :: k
      INTEGER :: info
      REAL, DIMENSION(nlong,mnlat,nlev) :: global
      REAL, DIMENSION(nlong,mnlat) :: global_level

      COMMON /SUM3_D/global

      global=0.
      global_level=0.0
      global_level(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1)=x
      CALL GC_SSYNC(nproc,info)
      CALL GC_RSUM(nlong*mnlat,nproc,info,global_level)
      CALL GC_SSYNC(nproc,info)
      global(:,:,1)=global_level

      END SUBROUTINE SUM2D
#endif
