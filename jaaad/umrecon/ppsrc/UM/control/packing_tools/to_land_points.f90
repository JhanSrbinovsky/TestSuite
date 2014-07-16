

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine TO_LAND_POINTS-----------------------------------------
!LL
!LL  Purpose:  Selects land points values from input horizontal field DA
!LL            and writes then as contiguous values to array
!LL            DATA_LAND_POINTS.
!LL
!LL  Written by A. Dickinson 11/9/90
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   5.4    06/09/02 SX now uses this deck - add def UTILIO.
!LL                                                  E.Leung
!LL   5.5    25/04/03 Add def FLDIO for future mpp use.   E.Leung
!     6.0    11/09/03 Added def PPTOANC.                  P.Dando
!LL
!LL Programming standard :
!LL
!LL Logical components covered : S171
!LL
!LL Project task :
!LL
!LL  Documentation:  None
!LL
!LL  ------------------------------------------------------------------
!
!*L  Arguments:--------------------------------------------------------

      SUBROUTINE TO_LAND_POINTS                                         &
     & (DATA,DATA_LAND_POINTS,LAND_SEA_MASK,POINTS,LAND_POINTS)

      IMPLICIT NONE

      INTEGER                                                           &
     & POINTS                                                           &
              !IN Total no of both land and sea points to be processed
     &,LAND_POINTS !OUT No of land points

      LOGICAL                                                           &
     & LAND_SEA_MASK(POINTS)    !IN Land-sea mask

      REAL                                                              &
     & DATA_LAND_POINTS(POINTS)                                         &
                                !OUT Data on land points only
     &,DATA(POINTS)             !IN Data on land and sea points

! Workspace usage:-----------------------------------------------------
      INTEGER INDEX_LAND_POINTS(POINTS) !Gather index for land points
!----------------------------------------------------------------------
! External subroutines called:-----------------------------------------
!*---------------------------------------------------------------------
! Local varables:------------------------------------------------------
      INTEGER I ! Integer index
!----------------------------------------------------------------------

!L Compute gather index for land points

      LAND_POINTS = 0
      DO I=1,POINTS
        IF(LAND_SEA_MASK(I))THEN
          LAND_POINTS=LAND_POINTS + 1
          INDEX_LAND_POINTS(LAND_POINTS) = I
        END IF
      END DO

!L Gather land points from input array DATA

      DO I=1,LAND_POINTS
        DATA_LAND_POINTS(I)=DATA(INDEX_LAND_POINTS(I))
      ENDDO

      RETURN
      END SUBROUTINE TO_LAND_POINTS
