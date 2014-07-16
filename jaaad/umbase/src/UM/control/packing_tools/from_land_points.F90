#if defined(C94_1A) || defined(RECON) || defined(UTILIO) \
 || defined(FLDIO)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine FROM_LAND_POINTS---------------------------------------
!LL
!LL  Written by A. Dickinson
!LL
!LL  Purpose:  Selects successive values from input array
!LL            DATA_LAND_POINTS and writes them to land points
!LL            in horizontal field array DATA.
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!     5.4    06/09/02 SX now uses this deck - add def UTILIO.
!                                                    E.Leung
!     5.5    25/04/03 Add def FLDIO for future MPP use.   E.Leung
!LL
!LL  Documentation: None
!LL
!LL  ------------------------------------------------------------------
!
!*L  Arguments:--------------------------------------------------------

      SUBROUTINE FROM_LAND_POINTS                                       &
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
                                !IN Data on land points only
     &,DATA(POINTS)             !OUT Data on land and sea points

! Workspace usage:-----------------------------------------------------
      INTEGER INDEX_LAND_POINTS(POINTS) !Scatter index for land points
!----------------------------------------------------------------------
! External subroutines called:-----------------------------------------
!*---------------------------------------------------------------------
! Local varables:------------------------------------------------------
      INTEGER I ! Integer index
!----------------------------------------------------------------------
! Constants from comdecks:---------------------------------------------
#include "c_mdi.h"
#include "parvars.h"
!----------------------------------------------------------------------

!L Initialise sea points to MDI

      DO I=1,POINTS
        DATA(I)=RMDI
      ENDDO

!L Calculate scatter index for land points

      LAND_POINTS = 0
      DO I=1,POINTS
        IF(LAND_SEA_MASK(I))THEN
          LAND_POINTS=LAND_POINTS + 1
          INDEX_LAND_POINTS(LAND_POINTS) = I
        END IF
      END DO

!L Scatter land points to array DATA
!DIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO I=1,LAND_POINTS
        DATA(INDEX_LAND_POINTS(I))=DATA_LAND_POINTS(I)
      ENDDO

      RETURN
      END SUBROUTINE FROM_LAND_POINTS
#endif
