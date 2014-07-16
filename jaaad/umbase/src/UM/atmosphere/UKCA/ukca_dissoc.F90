#if defined(A34_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Module to contain photolysis rate arrays
!   Contains subroutines: ukca_strat_photol_init
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
!  Current Code Owner:       Olaf Morgenstern
!                            Colin Johnson
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
     MODULE UKCA_DISSOC

     IMPLICIT NONE

     REAL, ALLOCATABLE :: AJHNO3 (:)
     REAL, ALLOCATABLE :: AJPNA  (:)
     REAL, ALLOCATABLE :: AJH2O2 (:)
     REAL, ALLOCATABLE :: AJ2A   (:)
     REAL, ALLOCATABLE :: AJ2B   (:)
     REAL, ALLOCATABLE :: AJ3    (:)
     REAL, ALLOCATABLE :: AJ3A   (:)
     REAL, ALLOCATABLE :: AJCNITA(:)
     REAL, ALLOCATABLE :: AJCNITB(:)
     REAL, ALLOCATABLE :: AJBRNO3(:)
     REAL, ALLOCATABLE :: AJBRCL (:)
     REAL, ALLOCATABLE :: AJOCLO (:)
     REAL, ALLOCATABLE :: AJCL2O2(:)
     REAL, ALLOCATABLE :: AJHOCL (:)
     REAL, ALLOCATABLE :: AJNO   (:)
     REAL, ALLOCATABLE :: AJNO2  (:)
     REAL, ALLOCATABLE :: AJN2O5 (:)
     REAL, ALLOCATABLE :: AJNO31 (:)
     REAL, ALLOCATABLE :: AJNO32 (:)
     REAL, ALLOCATABLE :: AJBRO  (:)
     REAL, ALLOCATABLE :: AJHCL  (:)
     REAL, ALLOCATABLE :: AJN2O  (:)
     REAL, ALLOCATABLE :: AJHOBR (:)
     REAL, ALLOCATABLE :: AJF11  (:)
     REAL, ALLOCATABLE :: AJF12  (:)
     REAL, ALLOCATABLE :: AJH2O  (:)
     REAL, ALLOCATABLE :: AJCCL4 (:)
     REAL, ALLOCATABLE :: AJF113 (:)
     REAL, ALLOCATABLE :: AJF22  (:)
     REAL, ALLOCATABLE :: AJCH3CL(:)
     REAL, ALLOCATABLE :: AJC2OA (:)
     REAL, ALLOCATABLE :: AJC2OB (:)
     REAL, ALLOCATABLE :: AJMHP  (:)
     REAL, ALLOCATABLE :: AJCH3BR(:)
     REAL, ALLOCATABLE :: AJMCFM (:)
     REAL, ALLOCATABLE :: AJCH4  (:)
     REAL, ALLOCATABLE :: AJF12B1(:)
     REAL, ALLOCATABLE :: AJF13B1(:)
     REAL, ALLOCATABLE :: AJCOF2 (:)
     REAL, ALLOCATABLE :: AJCOFCL(:)
     REAL, ALLOCATABLE :: AJCO2  (:)
     REAL, ALLOCATABLE :: AJCOS  (:)

     CONTAINS

!  UKCA stratospheric photolysis module
!  Test version
! ######################################################################
!
! Subroutine Interface:
!
!---------------------------------------------------------------------------
! Subroutine UKCA_STRAT_PHOTOL_INIT
!------------------------------------------------------------------------
!
! This routine computes stratospheric photolysis rates and merges the
! rates, where necessary, with the tropospheric rates. This is done for
! one level at a time. The stratospheric photolysis routines are taken
! from SLIMCAT.
!
! Version history:
! ----------------
!
! v1.1 Original code          Olaf Morgenstern 22/1/2004
! v1.2 Introduce pointers for photolysis reactions instead of searching
!      matching ASAD names at every call of the subroutine
!                             Olaf Morgenstern 7/7/2005
! v1.2 Modified to include CO2 photolysis
!                             Olaf Morgenstern 5/1/2006
! v1.3 Modified to include branching from HO2NO2 and BrONO2 photolysis.
!                             Olaf Morgenstern 12/1/2006
! v1.4 Modified to include branching from O2 photolysis
! Current code owner: Olaf Morgenstern.
!
      SUBROUTINE UKCA_STRAT_PHOTOL_INIT
      IMPLICIT NONE

#include "parvars.h"
#include "typsize.h"

      ALLOCATE(AJHNO3 (theta_field_size))
      ALLOCATE(AJPNA  (theta_field_size))
      ALLOCATE(AJH2O2 (theta_field_size))
      ALLOCATE(AJ2A   (theta_field_size))
      ALLOCATE(AJ2B   (theta_field_size))
      ALLOCATE(AJ3    (theta_field_size))
      ALLOCATE(AJ3A   (theta_field_size))
      ALLOCATE(AJCNITA(theta_field_size))
      ALLOCATE(AJCNITB(theta_field_size))
      ALLOCATE(AJBRNO3(theta_field_size))
      ALLOCATE(AJBRCL (theta_field_size))
      ALLOCATE(AJOCLO (theta_field_size))
      ALLOCATE(AJCL2O2(theta_field_size))
      ALLOCATE(AJHOCL (theta_field_size))
      ALLOCATE(AJNO   (theta_field_size))
      ALLOCATE(AJNO2  (theta_field_size))
      ALLOCATE(AJN2O5 (theta_field_size))
      ALLOCATE(AJNO31 (theta_field_size))
      ALLOCATE(AJNO32 (theta_field_size))
      ALLOCATE(AJBRO  (theta_field_size))
      ALLOCATE(AJHCL  (theta_field_size))
      ALLOCATE(AJN2O  (theta_field_size))
      ALLOCATE(AJHOBR (theta_field_size))
      ALLOCATE(AJF11  (theta_field_size))
      ALLOCATE(AJF12  (theta_field_size))
      ALLOCATE(AJH2O  (theta_field_size))
      ALLOCATE(AJCCL4 (theta_field_size))
      ALLOCATE(AJF113 (theta_field_size))
      ALLOCATE(AJF22  (theta_field_size))
      ALLOCATE(AJCH3CL(theta_field_size))
      ALLOCATE(AJC2OA (theta_field_size))
      ALLOCATE(AJC2OB (theta_field_size))
      ALLOCATE(AJMHP  (theta_field_size))
      ALLOCATE(AJCH3BR(theta_field_size))
      ALLOCATE(AJMCFM (theta_field_size))
      ALLOCATE(AJCH4  (theta_field_size))
      ALLOCATE(AJF12B1(theta_field_size))
      ALLOCATE(AJF13B1(theta_field_size))
      ALLOCATE(AJCOF2 (theta_field_size))
      ALLOCATE(AJCOFCL(theta_field_size))
      ALLOCATE(AJCO2  (theta_field_size))
      ALLOCATE(AJCOS  (theta_field_size))

      RETURN
      END SUBROUTINE UKCA_STRAT_PHOTOL_INIT

      END MODULE UKCA_DISSOC
#endif
