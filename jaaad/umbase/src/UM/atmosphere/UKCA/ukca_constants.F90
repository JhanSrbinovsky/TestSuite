#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module to contain constants used in UKCA
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Current Code Owner:       Colin Johnson/Olaf Morgenstern
!                           Fiona O'Connor
!
!  Code Description:
!   Language:  FORTRAN 90 (formatted)
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
      MODULE UKCA_CONSTANTS

#include "cphyscon.h"      
#include "c_sulchm.h"    
#include "c_rmol.h"  
#include "c_v_m.h"  

! cphyscon includes:
! c_a:                 Earth_radius
! c_g:                 g
! c_lheat:         Lc, Lf
! c_r_cp:        R, Cp, Kappa, Pref
! c_vkman:        vkman
! c_epsilon:        epsilon (ratio of Molec. Wts. of water and air)
! c_omega:        angular rotation
! c_pi:                pi, pi_over_180, recip_pi_over_180
! c_kt_ft:        knot and foot conversion factors
! c_mdi:        Missing data indicators, rmdi, imdi etc.
! c_lapse:        lapse, lapse_trop


! UKCA_ MODE constants

! MODE constants already defined in the UM
! 1) Original values
!      REAL, PARAMETER :: PPI=3.14159265358979323846
!      REAL, PARAMETER :: AVC=6.022E23
!      REAL, PARAMETER :: ZBOLTZ=1.3807E-23
!      REAL, PARAMETER :: VKARMN=0.4
!      REAL, PARAMETER :: RA=287.05
!      REAL, PARAMETER :: GG=9.80665
!      REAL, PARAMETER :: RAD_E=6371229.0
!      REAL, PARAMETER :: MM_DA=AVC*ZBOLTZ/RA
!      REAL, PARAMETER :: RHOSUL=1800.0E0      ! UM is slightly different
!      REAL, PARAMETER :: MMSUL=0.098E0
!      REAL, PARAMETER :: RR=8.314
!      PARAMETER(BCONST=3.44E13)        ! Volume_mode
!      PARAMETER(NU_H2SO4=3.0)          ! Volume_mode
!      PARAMETER(CONVERT=1.0e-18)       ! Volume_mode
!      PARAMETER(RHOW=1000.0,MMW=0.018) ! Volume_mode
!      PARAMETER(RHOSUL=1800.0)`        ! Volume_mode

! 2) UM definitions of MODE constants
      REAL, PARAMETER :: PPI=pi
      REAL, PARAMETER :: AVC=Avogadro
      REAL, PARAMETER :: ZBOLTZ=Boltzmann
      REAL, PARAMETER :: VKARMN=vKman
      REAL, PARAMETER :: RA=R
      REAL, PARAMETER :: RR=rmol
      REAL, PARAMETER :: GG=g
      REAL, PARAMETER :: RAD_E=Earth_Radius
      REAL, PARAMETER :: RHOSUL=rho_so4       ! 1769 or 1800 ?

! 3) MODE constants not already in UM
      REAL, PARAMETER :: MMSUL=0.09808E0     ! M. Wt. H2SO4 kg/mol
      REAL, PARAMETER :: MMW=0.018           ! M. Wt. H2O   kg/mol
      REAL, PARAMETER :: MM_DA=AVC*ZBOLTZ/RA !
      REAL, PARAMETER :: NMOL=1.0E2          !
      REAL, PARAMETER :: TDAYS=0.0           !
      REAL, PARAMETER :: EMS_EPS=1.0E-8      !
      REAL, PARAMETER :: CONC_EPS=1.0E-8     !
      REAL, PARAMETER :: DN_EPS=1.0E-8       !
      REAL, PARAMETER :: BCONST=3.44E13      ! Volume_mode
      REAL, PARAMETER :: NU_H2SO4=3.0        ! Volume_mode
      REAL, PARAMETER :: CONVERT=1.0E-18     ! Volume_mode
      REAL, PARAMETER :: RHOW=1.0E3          ! Volume_mode  kg/l

! 4) ALLOCATABLES
      REAL, ALLOCATABLE, SAVE :: d0(:)       ! Diffusion coefficients
      REAL, ALLOCATABLE, SAVE :: rsurf(:,:)  ! Standard surface resistances

!     molecular masses in kg/mol

      REAL, PARAMETER :: m_air        =  0.02897

      END MODULE UKCA_CONSTANTS
#endif
