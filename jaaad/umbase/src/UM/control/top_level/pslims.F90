#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Returns the limits of the pseudo levels
! Subroutine Interface:
      SUBROUTINE PSLIMS(IPFIRST,IPLAST,IFIRST,ILAST)
      USE CSENARIO_MOD
      IMPLICIT NONE

! Description:
!
! Method:
!
! Current code owner:  S.J.Swarbrick
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
! Subroutine arguments:
!   Scalar arguments with intent(in)
      INTEGER IPFIRST
      INTEGER IPLAST

!   Scalar arguments with intent(out):
      INTEGER IFIRST
      INTEGER ILAST

! Global variables:
#include "csubmodl.h"
#include "version.h"
#include "parparm.h"
#include "typsize.h"
#include "nstypes.h"
#include "model.h"

!- End of Header --------------------------------------------------

      IF(IPFIRST == 1) THEN
        IFIRST=1
      ELSE IF(IPFIRST == 21) THEN
        IFIRST=AASPF(1)
      ELSE IF(IPFIRST == 22) THEN
        IFIRST=AASPF(2)
      ELSE IF(IPFIRST == 23) THEN
        IFIRST=AASPF(3)
      ELSE IF(IPFIRST == 24) THEN
        IFIRST=AASPF(4)
      ELSE IF(IPFIRST == 25) THEN
        IFIRST=AASPF(5)
      ELSE
        write(6,*) 'S: PSLIMS(PRELIM).INVALID IPFIRST. NO CHECKING'
        IFIRST=1
      END IF

      IF(IPLAST == 1) THEN
        ILAST=SWBND
      ELSE IF(IPLAST == 2) THEN
        ILAST=LWBND
      ELSE IF(IPLAST == 6) THEN
!Sulphate loading patterns
        ILAST=NSULPAT
      ELSE IF(IPLAST == 7) THEN
!Direct vegetation parametrization: all surface types
        ILAST=NTYPE
      ELSE IF(IPLAST == 8) THEN
!Direct vegetation parametrization: only plant functional types
        ILAST=NPFT
      ELSE IF(IPLAST == 9) THEN
!Direct vegetation parametrization: all tiles
        ILAST=NTILES
      ELSE IF(IPLAST == 21) THEN
        ILAST=AASPL(1)
      ELSE IF(IPLAST == 22) THEN
        ILAST=AASPL(2)
      ELSE IF(IPLAST == 23) THEN
        ILAST=AASPL(3)
      ELSE IF(IPLAST == 24) THEN
        ILAST=AASPL(4)
      ELSE IF(IPLAST == 25) THEN
        ILAST=AASPL(5)
      ELSE
        write(6,*) 'S: PSLIMS(PRELIM). INVALID IPLAST. NO CHECKING'
        ILAST=1000
      END IF

      RETURN
      END SUBROUTINE PSLIMS
#endif
