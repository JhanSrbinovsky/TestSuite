#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE ETA2P_V(nfill,pos,lnp,val)
!
! Current Owner of Code: M.G. Sanderson
!
! History:
! Version   Date                    Comment
!  5.5    22/01/04  Created. Vector version of ETA2P. K. Ketelsen
!  6.1    06/08/04  No change
!
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE LINTERP_MOD
      IMPLICIT NONE
!----------------------------------------------------------------------

      INTEGER,                                  INTENT(IN) :: nfill
      REAL, DIMENSION(4,nfill),                 INTENT(IN) :: pos ! Posi
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN) :: lnp
      REAL, DIMENSION(nfill),                  INTENT(OUT) :: val

      CALL LINTERP_V(nfill,pos,lnp,val,.TRUE.,.FALSE.)

      val = EXP(val)

      RETURN

      END SUBROUTINE ETA2P_V
#endif
