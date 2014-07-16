#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      MODULE P2ETA_MOD

      INTERFACE P2ETA_V
        MODULE PROCEDURE P2ETA_V
      END INTERFACE ! P2ETA_V

      CONTAINS

      SUBROUTINE P2ETA_V(nfill,todo,press,pos,p,pos_out)
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.5   13/02/04   Created. K. Ketelsen
!   6.1   24/09/04   Minor changes for F90 compilation. M.G. Sanderson
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTERFACE
! DEPENDS ON: getmetpoint
        REAL FUNCTION GETMETPOINT(POS,FIELD,Lhalf,LU)
          USE IN_STOCHEM_GRD
          REAL, DIMENSION(4),             INTENT(IN) :: pos
          REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: field
          LOGICAL,                        INTENT(IN) :: Lhalf  ! T for 1
          LOGICAL,              INTENT(IN), OPTIONAL :: lu
END FUNCTION GETMETPOINT
      END INTERFACE

      INTEGER,                                  INTENT(IN) :: nfill
      LOGICAL, DIMENSION(nfill),                INTENT(IN) :: todo
      REAL, DIMENSION(nfill),                   INTENT(IN) :: press
      REAL, DIMENSION(4,nfill),                 INTENT(IN) :: pos
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN) :: p
      REAL, DIMENSION(:,:),                    INTENT(OUT) :: pos_out

      INTEGER                   :: l,k
      CHARACTER(LEN=72)         :: cmessage
      REAL                      :: xxx
      REAL, DIMENSION(nfill)    :: eta_low
      REAL, DIMENSION(nfill)    :: eta_high
      REAL, DIMENSION(nfill)    :: plow
      REAL, DIMENSION(nfill)    :: phigh
      REAL, DIMENSION(nfill)    :: pre_lo
      LOGICAL, DIMENSION(nfill) :: flag,error_flag
      INTEGER, DIMENSION(nfill) :: ind
      REAL, DIMENSION(4,nfill)  :: posr
      INTEGER                   :: nv

!kk   In the test case, count(todo) was less than 10%
!kk   therefore a compress loop was performed to reduce the
!kk   loop count of the search loop.

      nv = 0
!CDIR nodep
      DO l=1,nfill
        IF (todo(l)) THEN
          nv = nv + 1
          ind(nv) = l
          posr(1,nv) = pos(1,l)
          posr(2,nv) = pos(2,l)
          posr(3,nv) = pos(3,l)
          posr(4,nv) = pos(4,l)
          pre_lo(nv) = press(l)
        END IF
      END DO

      DO l=1,nv
! DEPENDS ON: getmetpoint
        phigh(l) = GETMETPOINT(posr(:,l),p(:,:,0),.TRUE.)
      END DO

!kk   vectorized search loop

      flag(1:nv) = .true.
      DO k=1,nmetlev
        DO l=1,nv
          IF (flag(l)) THEN
            plow(l) = phigh(l)
! DEPENDS ON: getmetpoint
            phigh(l) = GETMETPOINT(posr(:,l),p(:,:,k),.TRUE.)
            IF (pre_lo(l) > phigh(l)) THEN
              eta_high(l) = eta_theta(k)
              eta_low(l) = eta_theta(k-1)
              flag(l) = .false.
            END IF
          END IF
        END DO
      END DO

!CDIR nodep
      DO l=1,nv
        IF (pre_lo(l) < plow(l)) THEN
! Logarithmic interpolation in pressure
          xxx = eta_low(l) + (eta_high(l)-eta_low(l)) *                 &
     &      LOG(pre_lo(l)/plow(l)) / LOG(phigh(l)/plow(l))
        ELSE
          xxx = eta_low(l) !If PRESS is slightly below P0
        END IF
        pos_out(3,ind(l)) = xxx
        error_flag(l) =                                                 &
     &    (xxx < eta_stochem(0) .OR. xxx > eta_stochem(nlev))
      END DO

!kk   In case of Error, Rerun complete loop and abort

      IF (ANY(error_flag(1:nv))) THEN
        DO l=1,nfill
          IF (todo(l)) THEN
            IF (pos(3,l) < eta_stochem(0) .OR.                          &
     &        pos(3,l) > eta_stochem(nlev)) THEN
              cmessage = 'Position out of bounds'
              WRITE(6,*) cmessage,'pos(3): ',pos(3,l),' L: ',l
! DEPENDS ON: ereport
              CALL EREPORT('CLMIX2',1,cmessage)
            END IF
          END IF
        END DO
      END IF

      RETURN

      END SUBROUTINE P2ETA_V

      END MODULE P2ETA_MOD

#endif
