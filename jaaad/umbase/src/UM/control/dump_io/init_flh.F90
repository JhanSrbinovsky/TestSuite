#if defined(C80_1A) || defined(PPTOANC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE INIT_FLH --------------------------------------
!LL
!LL  Written by D. Robinson
!LL
!LL  Model   Date     Modification history
!LL version
!LL   3.4   08/09/94  New routine.
!LL
!LL  Programming standard: Unified Model Documentation Paper No 3
!LL                        Version No 4  5/2/92
!LL
!LL  System component: R30
!LL
!LL  System task: F3
!LL
!LL  Purpose:
!LL           Initialises the fixed length header to IMDI except
!LL           for all array dimensions which are set to 1.
!LL
!LL  Documentation:
!LL           Unified Model Documentation Paper No F3
!LL           Version No 5 9/2/90
!LL
!LL------------------------------------------------------------
!*L Arguments:-------------------------------------------------
      SUBROUTINE INIT_FLH (FIXHD,LEN_FIXHD)

      IMPLICIT NONE

      INTEGER                                                           &
     & LEN_FIXHD                                                        &
                        ! IN    Length of fixed length header
     &,FIXHD(LEN_FIXHD) ! INOUT Fixed length header

! Local arrays:------------------------------------------------
! None
! -------------------------------------------------------------
#include "c_mdi.h"
! External subroutines called:---------------------------------
! None
! Local variables:---------------------------------------------
      INTEGER J
! -------------------------------------------------------------


! 1.0 Initialise to IMDI
      DO J = 1,LEN_FIXHD
        FIXHD(J) = IMDI
      ENDDO

! 2.0 Set all array dimensions to 1
      FIXHD(101) = 1     !  Integer Constants
      FIXHD(106) = 1     !  Real Constants
      FIXHD(111) = 1     !  1st dim - Level dependent constants
      FIXHD(112) = 1     !  2nd dim - Level dependent constants
      FIXHD(116) = 1     !  1st dim - Row dependent constants
      FIXHD(117) = 1     !  2nd dim - Row dependent constants
      FIXHD(121) = 1     !  1st dim - Column dependent constants
      FIXHD(122) = 1     !  2nd dim - Column dependent constants
      FIXHD(126) = 1     !  1st dim - Field of constants
      FIXHD(127) = 1     !  2nd dim - Field of constants
      FIXHD(131) = 1     !  Extra constants
      FIXHD(136) = 1     !  Temp History file
      FIXHD(141) = 1     !  Compressed field Index 1
      FIXHD(143) = 1     !  Compressed field Index 2
      FIXHD(145) = 1     !  Compressed field Index 3
      FIXHD(151) = 1     !  1st dim - Lookup Table
      FIXHD(152) = 1     !  2nd dim - Lookup Table
      FIXHD(161) = 1     !  Data

      RETURN
      END SUBROUTINE INIT_FLH
#endif
