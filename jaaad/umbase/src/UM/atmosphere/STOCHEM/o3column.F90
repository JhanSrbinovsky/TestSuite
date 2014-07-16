#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE O3COLUMN(o3_ls,o3col)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Calculate o3 column 0-100mb from li+shine o3
!-                         and surface pressure
!-
!-   Inputs  : o3_ls
!-   Outputs : o3col
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  4.5    17/09/98  Created. D.S. Stevenson
!  5.3    17/09/01  O3_LS assumed to be on eta_theta levels. C.E.Johnson
!  6.1    21/10/04  Function HEIGHT renamed as ST_HEIGHT. M.G. Sanderson
!
!-
!VVV  v3.0  O3COLUMN 24/8/00 More general for different model tops.
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!----------------------------------------------------------------------
      REAL, DIMENSION(nlnpe,nlpe,nmetlev), INTENT(IN) :: o3_ls
      REAL, DIMENSION(nlnpe,nlpe), INTENT(OUT) :: o3col

      INTEGER :: k
      REAL :: f

! DEPENDS ON: st_height
      k = ST_HEIGHT(eta_stochem(nlev),'Eta_theta') ! ETA level below mod
      f = (eta_stochem(nlev)-eta_theta(k))/(eta_theta(k+1)-eta_theta(k))
      o3col = (o3_ls(:,:,k)+(o3_ls(:,:,k+1)-o3_ls(:,:,k))*f) / dob

      END SUBROUTINE O3COLUMN
#endif
