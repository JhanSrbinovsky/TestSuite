#if defined(C80_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE PR_REHDA---------------------------------------
!LL
!LL  Purpose: Prints out real constants record and checks
!LL           validity of information.
!LL
!LL  Written by A. Dickinson 28/12/89
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL
!LL  Programming standard: Unified Model Documentation Paper No 3
!LL                        Version No 1 15/1/90
!LL
!LL  System component: C25
!LL
!LL  System task: F3
!LL
!LL  Documentation: Unified Model Documentation Paper No F3
!LL           Version No 5 9/2/90
!LL
!LL------------------------------------------------------------
!*L Arguments:-------------------------------------------------
      SUBROUTINE PR_REHDA(REALHD,LEN_REALHD)

      IMPLICIT NONE

      INTEGER                                                           &
     & LEN_REALHD !IN Length of real header

      REAL                                                              &
     & REALHD(LEN_REALHD) !IN Real header

! -------------------------------------------------------------
! Workspace usage:---------------------------------------------
! None
! -------------------------------------------------------------
!*L External subroutines called:-------------------------------
! None
!*-------------------------------------------------------------
! Local variables:---------------------------------------------
      INTEGER I
!--------------------------------------------------------------

!L Internal structure: None

      WRITE(6,'('' '')')
      WRITE(6,'('' REAL CONSTANTS'')')
      WRITE(6,'('' --------------'')')
      WRITE(6,'('' E-W grid spacing in degrees -'',e12.4)')             &
     &REALHD(1)
      WRITE(6,'('' N-S grid spacing in degress -'',e12.4)')             &
     &REALHD(2)
      WRITE(6,'('' Latitude of first row in degrees -'',e12.4)')        &
     &REALHD(3)
      WRITE(6,'('' Longitude of first point in a row in degrees -'',    &
     &e12.4)')REALHD(4)
      WRITE(6,'('' Real latitude of pseudo North Pole in degrees - '',  &
     &e12.4)')REALHD(5)
      WRITE(6,'('' Real longitude of pseudo North Pole in degrees - '', &
     &e12.4)')REALHD(6)
      WRITE(6,'('' Grid orientation in degrees - '',                    &
     &e12.4)')REALHD(7)
      WRITE(6,'(8X,''                   Year        Day       Hour      &
     &Minute     Second'')')
      WRITE(6,'('' Atmosphere time = '',                                &
     &5e12.4)')(REALHD(I),I=8,12)

      WRITE(6,'('' Mass, energy, energy drift = '',3e12.4)')            &
     &REALHD(19),REALHD(20),REALHD(21)

      RETURN
      END SUBROUTINE PR_REHDA
#endif
