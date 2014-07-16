#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE STATS(Mconc,sdconc,nm,np,nph,nstat,nnn,conc,clist,     &
     &  nchem)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : ACCUMULATE STATISTICS
!-
!-   Inputs  : CONC,NNN
!-   Outputs : MCONC,SDCONC,NM,NSTAT
!-   Controls:
!-
!-   Created   9-DEC-1993   W.J. Collins
!-   Updated   7-JUN-1995   Bill Collins
!-                           NM is accumulated no. of cells per
!-                           grid volume, NSTAT is no. of timesteps
!-                                        of accumulation.
!-   Updated   6-AUG-1996   Bill Collins  Parameters now in INCLUDE
!-   Updated  17-OCT-1997   Bill Collins  Converted to Fortran 90
!-   Updated  10-MAR-1998   Bill Collins
!-                           Reduced size of Eulerian arrays
!-   Updated  11-JUN-1998   Colin Johnson Added 3-D counter
!-   Updated  11-NOV-1998   Bill Collins  Now copes with PH.
!-   Updated  03-Dec-2001   Bill Collins  Converted to New Dynamics
!-
!VVV  V2.2  STATS 20/X/99
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_OUT
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: nchem
      INTEGER, DIMENSION(nlnpe,nlpe,nlev), INTENT(IN) :: nnn
      INTEGER, DIMENSION(numchem), INTENT(IN) :: clist
      INTEGER, DIMENSION(nlnpe,nlpe,nlev), INTENT(INOUT) :: nm
      INTEGER, DIMENSION(nlnpe,nlpe,nlev), INTENT(INOUT) :: np
      INTEGER, DIMENSION(nlnpe,nlpe,nlev), INTENT(INOUT) :: nph
      INTEGER, INTENT(INOUT) :: nstat

      REAL, DIMENSION(nc,nlnpe,nlpe,nlev), INTENT(IN) :: conc
      REAL, DIMENSION(numchem,nlnpe,nlpe,nlev), INTENT(INOUT) :: mconc
      REAL, DIMENSION(numchem,nlnpe,nlpe,nlev), INTENT(INOUT) :: sdconc

      INTEGER :: k

      DO k=1,nchem
        mconc(k,:,:,:)=mconc(k,:,:,:)+conc(clist(k),:,:,:)
        sdconc(k,:,:,:)=sdconc(k,:,:,:)+conc(clist(k),:,:,:)**2
      END DO
!      Count NP only when at least one cell in grid volume
      WHERE(nnn /= 0) np=np+1
      nm=nm+nnn
      nstat=nstat+1
      END SUBROUTINE STATS
#endif
