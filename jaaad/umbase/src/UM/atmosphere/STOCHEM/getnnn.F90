#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE GETNNN(ipos,pos,nnn,nbl,bl,orog,nfill)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Fills IPOS,NNN and NBL
!-
!-   Inputs  : POS,BL
!-   Outputs : IPOS,NNN,NBL
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  4.5    18/11/98  Created.  W.J. Collins
!  5.3    14/08/01  New Dynamics version. C.E. Johnson
!  5.5    14/03/04  Changed to use new version of AINDEX. M.G. Sanderson
!  6.1    22/08/04  Minor tidying of code. M.G. Sanderson
!  6.2    28/03/06  Minor changes for vn6.2  M.G. Sanderson
!-
!VVV  V5.2.1 GETNNN 14/VIII/)1 - vn5.2.1
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: nfill

      REAL, DIMENSION(4,nclprc),            INTENT(IN) :: pos
      REAL, DIMENSION(nlonpe,nlatpe),       INTENT(IN) :: bl
      REAL, DIMENSION(nlonpe,nlatpe),       INTENT(IN) :: orog

      INTEGER, DIMENSION(3,nclprc),        INTENT(OUT) :: ipos
      INTEGER, DIMENSION(nlnpe,nlpe,nlev), INTENT(OUT) :: nnn
      INTEGER, DIMENSION(nlnpe,nlpe),      INTENT(OUT) :: nbl

      INTEGER :: i
      INTEGER :: j
      INTEGER :: k
      INTEGER :: l
      INTEGER, DIMENSION(NLEV) :: Total_in_level

      REAL :: Zorog
      REAL :: Zbl

      nnn=0
      nbl=0

! DEPENDS ON: aindex
      CALL AINDEX(nfill,pos,ipos)    ! calc new indices

      DO j=1,nfill
        i=ipos(1,j)-lndat+1
        k=ipos(2,j)-ltdat+1
        l=ipos(3,j)
        nnn(i,k,l)=nnn(i,k,l)+1
! Increment no. cells in b.l.
! DEPENDS ON: getmetpoint
        Zorog=GETMETPOINT(pos(:,j),orog,.TRUE.)
! DEPENDS ON: getmetpoint
        Zbl=GETMETPOINT(pos(:,j),bl,.TRUE.)
        IF(pos(4,j)-Earth_Radius-Zorog <= Zbl) nbl(i,k)=nbl(i,k)+1
      END DO

      DO k=1,nlev
        Total_in_level(k)=SUM(nnn(:,:,k))
      ENDDO

      WRITE(6,'(A18,10I5,/,10I5)') 'Total in Level: ',Total_in_level

      END SUBROUTINE GETNNN
#endif
