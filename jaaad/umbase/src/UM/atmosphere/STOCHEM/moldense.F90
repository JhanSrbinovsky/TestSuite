#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE MOLDENSE(m,pos,lnp,tl,nfill)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Return molecular densities in molecules/cm3
!-
!-   Inputs  : POS,LNP,TL,NFILL
!-   Outputs : M
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  5.2    13/12/00  Created. W.J. Collins
!  5.3    31/08/01  New dynamics version. C.E. Johnson
!  5.5    10/03/04  Calls ETA2P_V for improved optimisation. K. Ketelsen
!  6.1    21/10/04  No change
!
!-
!VVV  V5.2  MOLDENSE 31/8/01 - New dynamics version
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!----------------------------------------------------------------------
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev),  INTENT(IN) :: lnp
      REAL, DIMENSION(4,nclprc),                 INTENT(IN) :: pos
      REAL, DIMENSION(nclprc),                   INTENT(IN) :: tl
      REAL, DIMENSION(nclprc),                  INTENT(OUT) :: m
      INTEGER,                                   INTENT(IN) :: nfill

      INTEGER :: j

! DEPENDS ON: eta2p_v
      CALL ETA2P_V(nfill,pos,lnp,m)

! Factor of 1.0e-6 to convert from m-3 to cm-3
      DO j = 1, nfill
        m(j ) = m(j) * 1.0e-6 * na / (rmol * tl(j))
      END DO

      END SUBROUTINE MOLDENSE
#endif
