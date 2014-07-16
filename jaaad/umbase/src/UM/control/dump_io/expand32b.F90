#if defined(C80_1A) || defined(RECON) || defined(VAROPSVER) \
 || defined(UTILIO)
#if !defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE EXPAND32B--------------------------------------
!LL
!LL  Purpose: Expands from 32 to 64 bit for dump reading routines.
!LL
!LL MC          <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.3  08/04/94  Added to avoid problems in readdump. M.Carter
!LL   4.5  28/10/98  Introduce Single Column Model. J-C Thil.
!     5.3  22/11/01  Enable MPP as the only option for
!                    small executables         E.Leung
!LL
!LL  Programming standard: Unified Model Documentation Paper No 3
!LL                        Version No 1 15/1/90
!LL
!LL  Logical component: R30
!LL
!LL  System task: F3
!LL
!LL  Documentation: Unified Model Documentation Paper No F3
!LL                 Version No 5 9/2/90
!LLEND---------------------------------------------------------
!
!*L Arguments:-------------------------------------------------
      SUBROUTINE EXPAND32B(LENGTH, ARRAY, VERSION)

      IMPLICIT NONE

      INTEGER                                                           &
     & LENGTH,                                                          &
                     !IN length of the field to be expanded
     & VERSION       !IN model version

      REAL                                                              &
     & ARRAY(LENGTH)  !IN/OUT array to be expanded in place

! -------------------------------------------------------------
! Local variables: --------------------------------------------
      REAL HOLD(LENGTH)     ! space for expanded array
      INTEGER I             ! Loop index
! -------------------------------------------------------------
!*L External subroutines called:-------------------------------


! DEPENDS ON: expand21
      CALL EXPAND21(LENGTH,ARRAY,HOLD)
      DO I=1,LENGTH
        ARRAY(I)=HOLD(I)
      ENDDO

      RETURN
      END SUBROUTINE EXPAND32B
#endif
#endif
