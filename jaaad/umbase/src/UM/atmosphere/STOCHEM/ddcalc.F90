#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE DDCALC(dstore,dd,ddepo,ipos,pos,bl,orog,nbl,nfill)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Calculate dry depostion rates
!-
!-   Inputs  : DSTORE,DDEPO,IPOS,NC
!-   Outputs : DD
!-   Controls:
!-
!
! Current Code Owner: W.J. Collins
!
! History:
! Version   Date                    Comment
!  3.5    06/01/94  Created.  W.J. Collins
!  5.5    04/09/01  New dynamics version. W.J. Collins
!  5.5    29/11/01  More accurate BL test using orography. W.J. Collins
!  6.1    22/10/04  No change.
!
!-   Updated  29-Nov-2001   Bill Collins Added orog as an argument, more

!VVV  V5.2  DDCALC 4/IX/01 - New dynamics version
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: nfill
      INTEGER, DIMENSION(3,nclprc),   INTENT(IN) :: ipos
      INTEGER, DIMENSION(nlnpe,nlpe), INTENT(IN) :: nbl

      REAL, DIMENSION(nc,nlnpe,nlpe), INTENT(IN) :: ddepo
      REAL, DIMENSION(nc,nlnpe,nlpe), INTENT(IN) :: dstore
      REAL, DIMENSION(4,nclprc),      INTENT(IN) :: pos
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: bl
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: orog
      REAL, DIMENSION(nc,nclprc),     INTENT(OUT) :: dd

      INTEGER :: i
      INTEGER :: l
      INTEGER :: j

      REAL :: Zorog
      REAL :: Zbl

      DO j=1,nfill
        i=ipos(1,j)-lndat+1
        l=ipos(2,j)-ltdat+1
! DEPENDS ON: getmetpoint
        Zorog=GETMETPOINT(pos(:,j),orog,.TRUE.)
! DEPENDS ON: getmetpoint
        Zbl=GETMETPOINT(pos(:,j),bl,.TRUE.)

        IF(pos(4,j)-Earth_Radius-Zorog <= Zbl) THEN
          dd(:,j)=(ddepo(:,i,l)+dstore(:,i,l))/nbl(i,l)
        ELSE
          dd(:,j)=0.
        END IF
      END DO

      END SUBROUTINE DDCALC
#endif
