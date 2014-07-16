#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE TEMP(pos,t,t0,tl,nfill)
!-----------------------------------------------------------------------
!-
!-   Purpose and Methods : Return temperature of cell
!-
!-   Inputs  : POS,LONGM,LATM,T
!-   Outputs : TL
!-   Controls:
!-
!
! Current Code Owner: W.J. Collins
!
! History:
! Version   Date                    Comment
!  3.5    11/04/94  Created.  W.J. Collins
!  5.5    13/05/04  Vectorised code, using INTERP_V. M.G. Sanderson
!  6.1    22/10/04  No change.
!
!-
!VVV  V5.2  TEMP 10/IX/01  - New interp call
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_INTF
      USE INTERP_MOD
      IMPLICIT NONE
!-----------------------------------------------------------------------
      REAL, DIMENSION(nlonpe,nlatpe,nmetlev), INTENT(IN) :: t
      REAL, DIMENSION(nlonpe,nlatpe),         INTENT(IN) :: t0
      REAL, DIMENSION(4,nclprc),              INTENT(IN) :: pos
      INTEGER,                                INTENT(IN) :: nfill
      REAL, DIMENSION(nclprc),               INTENT(OUT) :: tl

      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev)           :: tdash
      LOGICAL, DIMENSION(nfill)                          :: todo

      tdash(:,:,0) = t0
      tdash(:,:,1:nmetlev) = t
      todo = .TRUE.

! cubic interpolation in vertical
      CALL INTERP_V(nfill,todo,pos(:,1:nfill),tdash,tl(1:nfill),        &
     &  .TRUE.,.FALSE.)

      END SUBROUTINE TEMP
#endif
