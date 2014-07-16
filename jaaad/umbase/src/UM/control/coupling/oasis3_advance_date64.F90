#if defined(OASIS3)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************

      SUBROUTINE OASIS3_ADVANCE_DATE64(ierr64)

! Description: This subroutine is a 64 bit wrapper interface
!              around the call to OASIS3_ADVANCE_DATE32.
!
! Author: R. Hill
!
! Current Code Owner : R. Hill
!=============================================================

      USE OASIS3_atmos_init

      IMPLICIT NONE

      INTEGER :: ierr64

      CALL OASIS3_ADVANCE_DATE32(ierr64)

      RETURN
      END SUBROUTINE OASIS3_ADVANCE_DATE64
#endif
