#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE WATER(pos,q,clw,rh,ql,liq,rhl,nfill)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Return water content of cell (mixing ratio)
!-
!-   Inputs  : POS,LONGM,LATM,Q,CLW
!-   Outputs : QL,LIQ
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  4.1    15/07/96  Created. W.J. Collins
!  5.5    22/10/04  Vectorised code. Now calls LINTERP_V. K. Ketelsen
!  6.1    22/10/04  Now only prints low humidity error if has actually
!                   reset some values. M.G. Sanderson
!  6.2    01/03/06  Interpolates relative humidity to parcel positions.
!                                               M.G. Sanderson
!
!-
!VVV  V5.2  WATER 3/IX/01  - New calcn. of indicies
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_INTF
      USE LINTERP_MOD
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: NFILL

      REAL, DIMENSION(4,nclprc),              INTENT(IN) :: pos
      REAL, DIMENSION(nlonpe,nlatpe,nmetlev), INTENT(IN) :: q
      REAL, DIMENSION(nlonpe,nlatpe,nmetlev), INTENT(IN) :: clw
      REAL, DIMENSION(nlonpe,nlatpe,nmetlev), INTENT(IN) :: rh
      REAL, DIMENSION(nclprc),               INTENT(OUT) :: ql
      REAL, DIMENSION(nclprc),               INTENT(OUT) :: liq
      REAL, DIMENSION(nclprc),               INTENT(OUT) :: rhl

      INTEGER :: j
      INTEGER :: n

      REAL    :: q7
      REAL    :: liq7
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev)  :: x

      x(:,:,0) = q(:,:,1)
      x(:,:,1:nmetlev) = q(:,:,1:nmetlev)

! Linear interpolation
      CALL LINTERP_V(nfill,pos,x,ql,.TRUE.,.FALSE.)

! Set to a minimum value
      n = 0
      DO j=1,nfill
        ql(j) = ql(j) * mair / mh2o ! convert g/g to volume mixing ratio
        IF (ql(j) < 1.0E-10) THEN
          ql(j) = 1.0E-10
          n = n + 1
        END IF
      END DO

      IF (n > 0) THEN
        WRITE(6,*) ' *** WATER: low or -ve humidity reset to 1.0E-10'
        WRITE(6,*) n,' times on PE ',mype
      END IF

      x(:,:,0) = clw(:,:,1)
      x(:,:,1:nmetlev) = clw(:,:,1:nmetlev)
      CALL LINTERP_V(nfill,pos,x,liq,.TRUE.,.FALSE.)

      x(:,:,0) = rh(:,:,1)
      x(:,:,1:nmetlev) = rh(:,:,1:nmetlev)
      CALL LINTERP_V(nfill,pos,x,rhl,.TRUE.,.FALSE.)

      END SUBROUTINE WATER
#endif
