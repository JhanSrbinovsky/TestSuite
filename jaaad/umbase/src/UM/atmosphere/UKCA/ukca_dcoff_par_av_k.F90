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
!    Calculates the average (over kth moment) diffusion coefficient
!    for a log-normally distributed aerosol population with
!    geometric mean diameter DG, geometric mean standard
!    deviation SSIGMA.
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
      FUNCTION UKCA_DCOFF_PAR_AV_K(K,DPG,SSIGMA,T,DVISC,MFPA)
!--------------------------------------------------------------------
!
! Purpose
! -------
! Calculates the average (over kth moment) diffusion coefficient
! for a log-normally distributed aerosol population with
! geometric mean diameter DG, geometric mean standard
! deviation SSIGMA following method in Regional Particulate
! Model as described by Binkowski & Shankar (1995).
! Also follows expression in Seinfeld & Pandis pg 474
! (Stokes-Einstein relation with slip-flow correction)
!
! Parameters
! ----------
! None
!
! Inputs
! ------
! K              : Index of moment for calculation
! DPG            : Geometric mean particle diameter for mode (m)
! SSIGMA         : Geometric standard deviation for mode
! T              : Temperature of air (K)
! DVISC          : Dynamic viscosity of air (kg m-1 s-1)
! MFPA           : Mean free path of air (m)
!
! Outputs
! -------
! UKCA_DCOFF_PAR_AV_K : Avg. ptcl diffusion coefficient (m^2 s-1)
!
! Local variables
! ---------------
! KNG            : Knudsen number for geo. mean sized particle
! LNSQSG         : ln(SSIGMA)*ln(SSIGMA)
! UKCA_PREF      : Prefactor term to expression
! TERM1,TERM2    : Terms in average diff coeff expression
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! ZBOLTZ         : Boltzmann's constant (kg m2 s-2 K-1 molec-1)
! PPI            : 3.1415927.....
!
!--------------------------------------------------------------------
      USE UKCA_CONSTANTS, ONLY: zboltz, ppi
!
      IMPLICIT NONE
!
! Subroutine interface
      INTEGER :: K
      REAL    :: UKCA_DCOFF_PAR_AV_K
      REAL    :: DPG
      REAL    :: SSIGMA
      REAL    :: T
      REAL    :: DVISC
      REAL    :: MFPA
!
! Local variables
      REAL    :: KNG
      REAL    :: LNSQSG
      REAL    :: UKCA_PREF
      REAL    :: TERM1
      REAL    :: TERM2

      TERM1=(-2.0*FLOAT(K)+1.0)/2.0
      TERM2=(-4.0*FLOAT(K)+4.0)/2.0
      KNG=2.0*MFPA/DPG
      LNSQSG=LOG(SSIGMA)*LOG(SSIGMA)
      UKCA_PREF=ZBOLTZ*T/(3.0*PPI*DVISC*DPG)
      UKCA_DCOFF_PAR_AV_K=UKCA_PREF*(EXP(TERM1*LNSQSG)+                 &
                           1.246*KNG*EXP(TERM2*LNSQSG))

      END FUNCTION UKCA_DCOFF_PAR_AV_K
#endif
