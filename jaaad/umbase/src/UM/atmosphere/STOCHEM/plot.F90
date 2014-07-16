#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE PLOT(conc,xx,ipos,nnn,d1,d2,nfill)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Transfer lagrangian concs to eulerian grid
!-                         and do mixing
!-
!-   Inputs  : XX,IPOS,NNN,D1,D2,NFILL
!-   Outputs : CONC,XX
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  3.4    09/12/93  Created. W.J. Collins
!  5.5    23/02/04  Commented out check on ozone value as inhibits
!                   vectorisation. M.G. Sanderson
!  6.1    21/10/04  No change
!
!VVV  V2.2  PLOT 20/X/99
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: nfill

      INTEGER, DIMENSION(3,nclprc),        INTENT(IN) :: ipos
      INTEGER, DIMENSION(nlnpe,nlpe,nlev), INTENT(IN) :: nnn

      REAL, INTENT(IN) :: d1
      REAL, INTENT(IN) :: d2
      REAL, DIMENSION(nc,nclprc),        INTENT(INOUT) :: xx
      REAL, DIMENSION(nc,nlnpe,nlpe,nlev), INTENT(OUT) :: conc

      INTEGER :: i
      INTEGER :: j
      INTEGER :: k
      INTEGER :: l
      INTEGER :: n
      INTEGER, DIMENSION(nc) :: indices=(/(l,l=1,nc)/) ! numbers 1 to NC
      REAL :: d

! Plot back to eulerian grid
! Grid concentrations calculated from cells within each grid
! Mixing of lagrangian cells within each grid

      conc = 0.0

      DO n=1,nfill
        i=ipos(1,n)-lndat+1
        j=ipos(2,n)-ltdat+1
        k=ipos(3,n)
        WHERE(indices/=i_ph) conc(:,i,j,k)=conc(:,i,j,k)+               &
                                                          ! Not pH
     &    xx(:,n)/nnn(i,j,k)
! Count number of non-zero pHs in this volume
!        NNNPH=COUNT(IPOS(1,:)==I.AND.IPOS(2,:)==J+ltdat-1.AND.
!     &    IPOS(3,:)==K.AND.XX(I_PH,:)/=0.)
!        IF(NNNPH>0) CONC(I_PH,I,J,K)=CONC(I_PH,I,J,K)+XX(I_PH,N)/NNNPH

!       IF(xx(i_o3,n)>20.0e-06)
!    &    WRITE(6,'(A,I3,A,I6,A,3I3,A,1PD10.4)')
!    &          ' *** PLOT: O3 >20000 ppb mype= ',mype,
!    &          ' CELL= ',n,' IPOS= ',i,j,k,' O3= ',xx(i_o3,n)
      END DO

! Mixing :
! Set lagrangian cells equal to eulerian concentration
! or exchange mass between eulerian and lagrangian

      DO n=1,nfill
        i=ipos(1,n)-lndat+1
        j=ipos(2,n)-ltdat+1
        k=ipos(3,n)
! Test for sensitivity to D in upper model.
! This magic no. should be coded as a parameter somewhere
        IF (Eta_stochem(k)>.23) THEN        ! i.e. < ~300 hPa
          d=d2
        ELSE
          d=d1
        END IF
        xx(:,n)=xx(:,n)+d*(conc(:,i,j,k)-xx(:,n))
      END DO

      END SUBROUTINE PLOT
#endif
