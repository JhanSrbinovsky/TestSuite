#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE FOLLOW(pos,xx,cellflux,clist,flist,follist,            &
     &    cellno,nchem,nflux,nfill,time,day,month,year)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : follow cells
!-
!-   Inputs  : POS,XX,CELLFLUX,CLIST,FLIST,CELLNO,NCHEM,NFLUX,
!-             NFILL,TIME,DAY,MONTH
!-   Outputs :
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  5.1    17/10/00  Created.  W.J. Collins
!  5.1    24/11/00  Revised. and now calls FINDNAME. C.E. Johnson
!  5.5    29/01/04  Uses num3dflux_dim to define cellflux array.
!                   M.G. Sanderson
!  6.1    21/10/04  No change.
!  6.2    21/10/05  Replace GSYNC with SSYNC. P.Selwood
!-
!VVV  V3.0  FOLLOW 23/XI/00 - FINDNAME called.
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_OUT
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: month
      INTEGER, INTENT(IN) :: year
      INTEGER, INTENT(IN) :: nchem
      INTEGER, INTENT(IN) :: nflux
      INTEGER, INTENT(IN) :: nfill
      INTEGER, DIMENSION(numchem), INTENT(IN)   :: clist
      INTEGER, DIMENSION(2,numflux), INTENT(IN) :: flist
      INTEGER, DIMENSION(numfol), INTENT(IN)    :: follist
      INTEGER, DIMENSION(nclprc), INTENT(IN)    :: cellno
      REAL, INTENT(IN) :: time
      REAL, INTENT(IN) :: day
      REAL, DIMENSION(4,nclprc), INTENT(IN)     :: pos
      REAL, DIMENSION(nc,nclprc), INTENT(IN)    :: xx
      REAL, DIMENSION(num3dflux_dim,nclprc), INTENT(IN) :: cellflux

      INTEGER :: i
      INTEGER :: j
      INTEGER :: k
      INTEGER :: info
      CHARACTER(LEN=14) :: filename

! DEPENDS ON: findname
      CALL FINDNAME('f','m','c',0,0,filename,month,year)

      DO k=0,nproc-1
        info=GC_SHM_PUT
        CALL GC_SSYNC(nproc,info)
        IF (mype==k) THEN
          OPEN(61,FILE=filename,STATUS='UNKNOWN',                       &
     &      POSITION='APPEND')
          IF (mype==0)                                                  &
     &      WRITE(61,'(//A,1p,E15.8,A,I7,A,0p,F5.2,A,I2//)')            &
     &      'TIME',TIME,'   NSTEP',nstep,'    DAY ',day,'    MONTH',    &
     &      month
          DO j=1,nfill
            DO i=1,numfol
              IF (follist(i)==cellno(j)) THEN
                WRITE(61,'(I6,4F14.5,I8,I4)') i,pos(:,j),follist(i),mype
!                WRITE(61,'(I6)') nchem
!                WRITE(61,'(4E18.10)') xx(clist(1:nchem),j)
!                WRITE(61,'(I6)') MAXVAL(flist(2,1:nflux))
!                WRITE(61,'(4E18.10)')
!     &            cellflux(1:MAXVAL(flist(2,1:nflux)),j)
              END IF
            END DO
          END DO
          CLOSE(61)
        END IF
      END DO

      END SUBROUTINE FOLLOW
#endif
