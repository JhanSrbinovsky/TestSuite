#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE AINDEX(nfill,pos,ipos)
!-----------------------------------------------------------------------
!-
!-   Purpose and Methods : RETURNS GRID INDICIES FOR A GRID POINT
!-
!-   Inputs  : POS
!-   Outputs : IPOS
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  4.0    27/07/95  Created.  W.J. Collins
!  5.0    29/05/01  New Dynamics version. C.E. Johnson
!  5.5    10/10/02  Now operates on whole arrays. M.G. Sanderson
!  6.1    20/10/04  Vectorised code. M.G. Sanderson
!  6.2    28/03/06  Minor changes for vn6.2  M.G. Sanderson
!
!-
!VVV  V5.0  AINDEX 29/V/01 - New Dynamics version
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER, INTENT(IN)                       :: nfill
      REAL,    DIMENSION(4,nclprc), INTENT(IN)  :: pos
      INTEGER, DIMENSION(3,nclprc), INTENT(OUT) :: ipos

      INTEGER :: j

      DO j = 1, nfill
        ipos(1,j) = 1 + INT(pos(1,j)/dlong)  ! Longitudinal grid index
        ipos(2,j) = 1 + INT(pos(2,j)/dlat)   ! Latitudinal grid index
        ipos(3,j) = 1                        ! Initialise vertical index
      END DO

      DO j = 1, nlev
        WHERE (pos(3,:) > eta_stochem(j)) ipos(3,:) = j + 1
      END DO

      END SUBROUTINE AINDEX
#endif
