#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE STORED(dstore,ddepo,nbl)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Store depositions if no cells are present
!-
!-   Inputs  : DDEPO,NBL
!-   Outputs : DSTORE
!-   Controls:
!-
!
! Current Code Owner: W.J. Collins
!
!
!  Model            Modification history from model version 6.0:
! version  Date
!  4.5   24/09/98   Created. W.J. Collins
!  6.0   12/09/03   Delete repeated declaration of NBL(NLNPE,NLPE).
!                   Introduce standard UM modification history.
!                                                         P.Dando
!  6.1   23/10/04   No change.
!  6.2   28/03/06   Minor changes for vn6.2  M.G. Sanderson
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER, DIMENSION(nlnpe,nlpe), INTENT(IN) :: nbl

      REAL, DIMENSION(nc,nlnpe,nlpe), INTENT(IN) :: ddepo
      REAL, DIMENSION(nc,nlnpe,nlpe), INTENT(OUT) :: dstore

      INTEGER :: k

! FIND grid boxes IN boundary LAYER with no cells

      DO k=1,nc
        WHERE(nbl == 0)
          dstore(k,:,:) = dstore(k,:,:) + ddepo(k,:,:)
        ELSE WHERE
          dstore(k,:,:) = 0.0
        END WHERE
      END DO

      END SUBROUTINE STORED
#endif
