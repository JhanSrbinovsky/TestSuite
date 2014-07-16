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
!    number conc. (ND) & total mass/ptcl (MDT) for each aerosol mode
!    Used for error checking in box model.
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
      SUBROUTINE UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT)
!----------------------------------------------------------------------
!
!     Purpose
!     -------
!     Used to calculate and print out max,min,mean of
!     number conc. (ND) & total mass/ptcl (MDT) for each aerosol mode
!     Used for error checking in box model
!
!     Inputs
!     ------
!     NBOX       : Number of grid boxes
!     ND         : Aerosol ptcl number density for mode (cm^-3)
!     MDT        : Avg tot mass of aerosol ptcl in mode (particle^-1)
!
!     Outputs
!     -------
!     None
!
!     Local variables
!     ---------------
!     XMINN,XMAXN,XMEANN : min,max,mean of particle number conc.
!     XMINM,XMAXM,XMEANM : min,max,mean of avg total mass per particle
!
!--------------------------------------------------------------------
      USE UKCA_MODE_SETUP
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NBOX
      REAL, INTENT(IN) :: ND(NBOX,NMODES)
      REAL, INTENT(IN) :: MDT(NBOX,NMODES)

! Local variables
      INTEGER :: IMODE
      INTEGER :: JL
      REAL    :: XMINN
      REAL    :: XMAXN
      REAL    :: XMEANN
      REAL    :: XMINM
      REAL    :: XMAXM
      REAL    :: XMEANM

      DO IMODE=1,NMODES
       IF(MODE(IMODE)) THEN
        XMINN=1.0e16
        XMAXN=-1.0e16
        XMEANN=0.0e0
        XMINM=1.0e16
        XMAXM=-1.0e16
        XMEANM=0.
        DO JL=1,NBOX
         XMINN=MIN(ND(JL,IMODE),XMINN)
         XMAXN=MAX(ND(JL,IMODE),XMAXN)
         XMEANN=XMEANN+ND(JL,IMODE)/FLOAT(NBOX)
         XMINM=MIN(MDT(JL,IMODE),XMINM)
         XMAXM=MAX(MDT(JL,IMODE),XMAXM)
         XMEANM=XMEANM+MDT(JL,IMODE)/FLOAT(NBOX)
        ENDDO
        write(6,'(1i4,6e15.6)') IMODE,XMINN,XMAXN,                      &
               XMEANN,XMINM,XMAXM,XMEANM
       ENDIF
      ENDDO

      RETURN
      END SUBROUTINE UKCA_CALCMINMAXNDMDT
#endif
