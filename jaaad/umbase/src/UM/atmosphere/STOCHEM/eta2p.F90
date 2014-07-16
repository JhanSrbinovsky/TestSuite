#if defined(A25_1A)
      REAL FUNCTION ETA2P(pos,lnp)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Converts from eta coordinates to pressure
!-
!-   Returned value  : ETA2P
!-   Inputs  : POS1,POS2,ETA,P_TH
!-   Outputs :
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  3.5    08/03/95  Created. W.J. Collins
!  5.3    29/08/01  New Dynamics version. C.E. Johnson
!  5.5    22/01/04  Vectorised. Now uses LINTERP as a module procedure
!                   K. Ketelsen
!  6.1    06/08/04  No change
!
!-
!VVV  V5.2  ETA2P 29/VIII/01 - New dynamics version
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE LINTERP_MOD
      IMPLICIT NONE
!----------------------------------------------------------------------
      REAL, DIMENSION(4),                       INTENT(IN) :: pos ! Posi
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN) :: lnp
      REAL :: x
      CHARACTER(LEN=72) :: cmessage

! No default settings for now, exit if eta is out of range
      if(pos(3)>1.0 .or. POS(3)<0.0) THEN
        cmessage='ETA out of range'
        WRITE(6,*) cmessage,'POS(3)=',pos(3)
! DEPENDS ON: ereport
        CALL EREPORT('ETA2P',1,cmessage)
      ENDIF

      x=LINTERP(pos,lnp,.TRUE.,.FALSE.)
      ETA2P=EXP(x)

      END FUNCTION ETA2P
#endif
