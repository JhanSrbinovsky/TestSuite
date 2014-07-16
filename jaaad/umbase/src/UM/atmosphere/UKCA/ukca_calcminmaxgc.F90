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
!    Used to calculate and print out max,min,mean of
!    each condensable gas phase concentration for
!    error checking in box model.
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
      SUBROUTINE UKCA_CALCMINMAXGC(STR_AT,NBOX,GC)
!----------------------------------------------------------------------
!
!     Purpose
!     -------
!     Used to calculate and print out max,min,mean of
!     each condensable gas phase concentration for
!     error checking in box model
!
!     Inputs
!     ------
!     STR_AT     : String indicating which process the code is up to
!     NBOX       : Number of grid boxes
!     GC         : Condensable cpt number density (molecules cm-3)
!
!     Outputs
!     -------
!     None
!
!     Local variables
!     ---------------
!     GCMIN,GCMAX,GCMEAN : min,max,mean of gas phase condensable conc.
!
!--------------------------------------------------------------------
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES
      IMPLICIT NONE

! Arguments
      INTEGER :: NBOX
      REAL    :: GC(NBOX,NCHEMG)
      CHARACTER*30 :: STR_AT

! Local variables
      INTEGER :: JL
      INTEGER :: JV
      REAL    :: GCMIN
      REAL    :: GCMAX
      REAL    :: GCMEAN
!
      DO JV=1,NCHEMG
       IF(CONDENSABLE(JV)) THEN
        WRITE(6,*) '**********************************'
        WRITE(6,*) STR_AT
        GCMIN=1.0e9
        GCMAX=-1.0e9
        GCMEAN=0.0e0
        DO JL=1,NBOX
         GCMIN=MIN(GC(JL,JV),GCMIN)
         GCMAX=MAX(GC(JL,JV),GCMAX)
         GCMEAN=GCMEAN+GC(JL,JV)/FLOAT(NBOX)
        ENDDO
        WRITE(6,*) 'GC:JV,Min,max,mean=',JV,GCMIN,GCMAX,GCMEAN
        WRITE(6,*) '**********************************'
       ENDIF
      ENDDO
!
      RETURN
      END SUBROUTINE UKCA_CALCMINMAXGC
#endif
