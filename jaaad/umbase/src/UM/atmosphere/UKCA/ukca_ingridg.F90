#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject 
! to the terms and conditions set out therein.
! [Met Office Ref SC138] 
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Generate sectional particle mass grid in molecules per particle.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
!  Current Code Owner:       Graham Mann
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
! Subroutine Interface:
      SUBROUTINE UKCA_INGRIDG(NBINS,MSMALL,MLARGE,MBMID,MBLO,MBHI)
!-----------------------------------------------------------------------
!
!     Purpose
!     -------
!     Generate sectional particle mass grid (in molecules per particle)
!
!     Inputs
!     ------
!     NBINS   : Number of aerosol ptcl size bins
!     MSMALL  : Number of molecules for ptcl at smallest bin mid-pt
!     MLARGE  : Number of molecules for ptcl at largest  bin mid-pt
!
!     Outputs
!     -------
!     MBLO, MBMID, MBHI : Mass of lower, mid and upper bin edge
!
!     Local variables
!     ---------------
!     LGMRANGE : Diff. in logarithms of masses of largest/smallest bins
!     LGMGRID  : Diff. in logarithms of masses of limits of each bin
!-----------------------------------------------------------------------
      IMPLICIT NONE
!
!     Subroutine interface
      INTEGER, INTENT(IN)  :: NBINS
      REAL, INTENT(IN)     :: MSMALL
      REAL, INTENT(IN)     :: MLARGE
      REAL, INTENT(OUT)    :: MBMID(NBINS)
      REAL, INTENT(OUT)    :: MBLO(NBINS)
      REAL, INTENT(OUT)    :: MBHI(NBINS)
!
!     Local variables
      INTEGER :: JV
      REAL    :: LGMRANGE
      REAL    :: LGMGRID
!
      LGMRANGE=LOG(MLARGE)-LOG(MSMALL)
      DO JV=1,NBINS+1
       LGMGRID=LOG(MSMALL)+LGMRANGE*FLOAT(JV-1)/FLOAT(NBINS)
       IF(JV < (NBINS+1)) MBLO(JV)=EXP(LGMGRID)
       IF(JV > 1) MBHI(JV-1)=EXP(LGMGRID)
      ENDDO
      DO JV=1,NBINS
       MBMID(JV)=EXP(0.5*LOG(MBLO(JV)*MBHI(JV))) ! geometric mean
      ENDDO

      RETURN
      END SUBROUTINE UKCA_INGRIDG
#endif
