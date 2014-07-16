#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Set the STASH addresses for D1
! Subroutine Interface:

!- End of subroutine code -------------------------------------------


!+Compute data lengths and addresses for primary fields
! Subroutine Interface:

!+Test whether level type is discrete (model) or continuous (non-model)
! Function Interface:
      LOGICAL FUNCTION DISCT_LEV(LEV_CODE,ErrorStatus,CMESSAGE)
      IMPLICIT NONE

! Description:
!
! Method:
!
! Current code owner:  S.J.Swarbrick
!
! History:
! Version   Date       Comment
! =======   ====       =======
!   4.1     Apr. 96    Original code.  S.J.Swarbrick
!   5.2     24/10/00   LBC fields (level type 0) stored as one array
!                      for all levels, but addressing requires
!                      seperate levels to be taken into account.
!                                                           P.Burton
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
#include "model.h"

! Function arguments:
!   Scalar arguments with intent(in):
      INTEGER LEV_CODE !Level code from STASHmaster

! ErrorStatus
      INTEGER ErrorStatus
      CHARACTER*80 CMESSAGE

!- End of Header ----------------------------------------------

      IF (LEV_CODE == 0 .OR. LEV_CODE == 1 .OR. LEV_CODE == 2 .OR.      &
     &    LEV_CODE == 6 .OR. LEV_CODE == 10) THEN
        DISCT_LEV=.TRUE.
      ELSE IF (LEV_CODE  >=  0 .AND. LEV_CODE  <=  10) THEN
        DISCT_LEV=.FALSE.
      ELSE
        DISCT_LEV=.FALSE.
        ErrorStatus=1
        CMESSAGE='DISCT_LEV : Invalid level type in STASHmaster'
      END IF
      END FUNCTION DISCT_LEV
!- End of Function code --------------------------------------------

!+Decode the STASH pseudo level code
! Subroutine Interface:
#endif
