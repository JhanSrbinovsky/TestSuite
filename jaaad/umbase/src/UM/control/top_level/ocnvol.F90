#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+
! Subroutine Interface:

      SUBROUTINE OCNVOL(LEN,IL1,IL2)

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
! Global variables:

#include "csubmodl.h"
#include "version.h"
#include "parvars.h"
#include "typsize.h"
#include "model.h"

! Subroutine arguments:

!   Scalar arguments with intent(in):

      INTEGER IL1
      INTEGER IL2

!   Scalar arguments with intent(out):

      INTEGER LEN

! Local scalars

      INTEGER ILVS

! - End of Header ----------------------------------------------------

      ILVS=IL2-IL1+1
      IF(ILVS /= NLEVSO) THEN
        WRITE(6,*)'OCNVOL: ERROR. OCEAN COMPRESED FIELD NOT ON',        &
     &  ' FULL MODEL LEVELS'
      ELSE

! Halo types are currently universally set to 1 for ocean prognostics
! so we use halo_type_single in the calculation of LEN.
        LEN=lasize(1,fld_type_p,halo_type_single)*                      &
     &                          lasize(2,fld_type_p,halo_type_single)

      END IF
      RETURN
      END SUBROUTINE OCNVOL
#endif
