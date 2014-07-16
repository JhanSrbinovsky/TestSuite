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
!    Calculates the average (over kth moment) gravitational
!    settling velocity for a log-normally distributed aerosol
!    population with geometric mean diameter DPG and geometric
!    mean standard deviation SSIGMA.
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
      FUNCTION UKCA_VGRAV_AV_K(K,DPG,SSIGMA,DVISC,MFPA,RHOP)
!------------------------------------------------------------------
!
! Purpose
! -------
! Calculates the average (over kth moment) gravitational
! settling velocity for a log-normally distributed aerosol
! population with geometric mean diameter DPG and geometric
! mean standard deviation SSIGMA following method in Regional
! Particulate Model as described by Binkowski & Shankar (1995).
!
! Parameters
! ----------
! None
!
! Inputs
! ------
! K          : Index of moment for calculation
! DPG        : Geometric mean particle diameter for mode (m)
! SSIGMA     : Geometric standard deviation for mode
! DVISC      : Dynamic viscosity of air (kg m-1 s-1)
! MFPA       : Mean free path of air (m)
! GG         : Gravitational acceleration = 9.80665 ms^-2
! RHOP       : Density of aerosol particle (kgm^-3)
!
! Outputs
! -------
! UKCA_VGRAV_AV_K : Avg. grav. settling velocity (m s-1)
!
! Local variables
! ---------------
! KNG        : Knudsen number for geo. mean sized particle
! UKCA_PREF  : Prefactor term to expression
! LNSQSG     : ln(SSIGMA)*ln(SSIGMA)
! TERM1,TERM2: Terms in average diff coeff expression
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! GG        : Gravitational acceleration = 9.80665 ms^-2
!
!--------------------------------------------------------------------
      USE UKCA_CONSTANTS
      IMPLICIT NONE
!
!     Subroutine interface
      INTEGER, INTENT(IN) :: k
      REAL,    INTENT(IN) :: dpg
      REAL,    INTENT(IN) :: ssigma
      REAL,    INTENT(IN) :: dvisc
      REAL,    INTENT(IN) :: mfpa
      REAL,    INTENT(IN) :: rhop
!
      REAL :: ukca_vgrav_av_k
!
!     Local variables
      REAL :: kng
      REAL :: ukca_pref
      REAL :: lnsqsg
      REAL :: term1
      REAL :: term2
!
      kng=2.0*mfpa/dpg
      ukca_pref=rhop*dpg*dpg*gg/(18.0*dvisc)
      lnsqsg=LOG(ssigma)*LOG(ssigma)
      term1=(4.0*FLOAT(k)+4.0)/2.0
      term2=(2.0*FLOAT(k)+1.0)/2.0
      ukca_vgrav_av_k=ukca_pref*(EXP(term1*lnsqsg)+                     &
                       1.246*kng*EXP(term2*lnsqsg))

      END FUNCTION UKCA_VGRAV_AV_K
#endif
