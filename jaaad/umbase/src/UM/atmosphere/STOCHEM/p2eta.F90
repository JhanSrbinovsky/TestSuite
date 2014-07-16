#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      REAL FUNCTION P2ETA(press,pos,p)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Converts from pressure coordinates to eta
!-
!-   Returned value  : P2ETA
!-   Inputs  : PRESS,P,POS
!-   Outputs :
!-   Controls:
!-
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  4.0  08/03/95   Created. W.J. Collins
!  4.4  17/10/97   Converted to F90 and parallelised. W.J. Collins
!  5.5  30/11/01   Major revisions for new dynamics. W.J. Collins
!  5.5  16/02/04   Redefined as a module procedure and vectorised
!                  form included for operation on SX6. K. Ketelsen
!  6.1  24/09/04   Reformatted for F90 compilation. M.G. Sanderson
!  6.2  31/01/06   Placed names of interface statement end blocks
!                  behind comments for portability. T. Edwards
!
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

      REAL,                                     INTENT(IN) :: press
      REAL, DIMENSION(4),                INTENT(IN) :: pos     !
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN) :: p

      INTEGER           :: k
      REAL              :: eta_low
      REAL              :: eta_high
      REAL              :: plow
      REAL              :: phigh
      CHARACTER(LEN=72) :: cmessage

! Find pressure at standard eta_theta levels

! DEPENDS ON: getmetpoint
      phigh = GETMETPOINT(pos,p(:,:,0),.TRUE.)
      eta_high = eta_theta(0)
!kk      IF(press>Phigh+tolerance) THEN
!kk        cmessage='**** ERROR: P>P0 ****'
!kk        WRITE(6,*) cmessage,'P=',PRESS,'P0=',Phigh
!kk        CALL EREPORT('P2ETA',1,cmessage)
!kk      ENDIF
      DO k=1,nmetlev
        plow = phigh
        eta_low = eta_high
! DEPENDS ON: getmetpoint
        phigh = GETMETPOINT(pos,p(:,:,k),.TRUE.)
        eta_high = eta_theta(k)
        IF (press > phigh) EXIT
      END DO
!kk      IF(press<Phigh) THEN
!kk        WRITE(cmessage,*) '**** ERROR: P<',Phigh,' ****'
!kk        WRITE(6,*) cmessage,'P=',press
!kk        CALL EREPORT('P2ETA',1,cmessage)
!kk      ENDIF
      IF (press < plow) THEN
! Logarithmic interpolation in pressure
        p2eta = eta_low + (eta_high-eta_low) *                          &
     &    LOG(press/plow) / LOG(phigh/plow)
      ELSE
        p2eta = eta_low !If PRESS is slightly below P0, set Eta to Eta_t
      END IF

      END FUNCTION P2ETA


#endif
