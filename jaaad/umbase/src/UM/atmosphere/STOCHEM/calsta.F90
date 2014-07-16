#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE CALSTA(sdconc,mconc,np,nph,clist,nchem)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : CALCULATE MEANS AND STANDARD DEVIATIONS
!-
!-   Inputs  : SDCONC,MCONC,NP,NPH,CLIST,NCHEM
!-   Outputs : SDCONC,MCONC
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins.
!
! History:
! Version   Date                    Comment
!  3.4    11/04/94  Created.  W.J. Collins
!  4.5    11/06/98  Uses 3D NP array instead of NSTAT. C.E. Johnson
!  4.5    11/11/98  Now handles PH. W.J. Collins
!  5.5    13/02/04  Loop order changed to avoid memory bank conflicts.
!                   (j now innermost loop). K. Ketelsen
!  6.1    21/10/04  Minor reformatting. M.G. Sanderson
!
!-
!VVV  V2.2  CALSTA 20/X/99
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_OUT
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: nchem
      INTEGER, DIMENSION(numchem), INTENT(IN) :: clist
      INTEGER, DIMENSION(nlnpe,nlpe,nlev), INTENT(IN) :: np
      INTEGER, DIMENSION(nlnpe,nlpe,nlev), INTENT(IN) :: nph
      REAL, DIMENSION(numchem,nlnpe,nlpe,nlev), INTENT(INOUT) :: mconc
      REAL, DIMENSION(numchem,nlnpe,nlpe,nlev), INTENT(INOUT) :: sdconc

      INTEGER :: i
      INTEGER :: j
      INTEGER :: k
      INTEGER :: l
      INTEGER :: n
      REAL    :: px
#include "c_mdi.h"

! Calculate mean and standard deviation

      DO i=1,nchem
        DO k=1,nlpe
          DO l=1,nlev
            DO j=1,nlnpe                !kk exchange loops (J inner)
              IF(clist(i)==i_ph) THEN ! for pH only count non-zero value
                n=nph(j,k,l)
              ELSE
                n=np(j,k,l)
              END IF
              IF (n > 1) THEN
                px=sdconc(i,j,k,l)-(mconc(i,j,k,l)**2.0)/n
                IF(px < 1.0e-40) THEN
                  sdconc(i,j,k,l)=0.0
                ELSE
                  sdconc(i,j,k,l)=sqrt(px/(n-1))
                END IF
              ELSE
                sdconc(i,j,k,l)=rmdi
              END IF
              IF (n > 0) THEN
                mconc(i,j,k,l)=mconc(i,j,k,l)/n
              ELSE
                mconc(i,j,k,l)=rmdi
              END IF
            END DO
          END DO
        END DO
      END DO

      END SUBROUTINE CALSTA
#endif
