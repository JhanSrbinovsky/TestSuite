#if defined(A19_1A) || defined(A19_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Counts the number of points containing each surface type and creates
! a TILE_INDEX array specifying the location of these points on the land
! grid.
!
! Subroutine Interface:

      SUBROUTINE TILEPTS(LAND_PTS,FRAC,TILE_PTS,TILE_INDEX)

! module for land-surface namelist
      USE LAND_SURF_MOD, ONLY :                                         &
     & ALL_TILES

      IMPLICIT NONE
!
! Description:
!
! Method:
!
! Current Code Owner:  Richard Betts
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   4.4    16/10/97   Original code.  Peter Cox
!   5.2    15/11/00   Re-Written for New Dynamics      M. Best
!
      INTEGER, INTENT(IN) :: LAND_PTS  !num land points to process

#include "nstypes.h"

      REAL, INTENT(IN) :: FRAC(LAND_PTS,NTYPE)  !fractions of surface types
      
      INTEGER, INTENT(OUT) ::  TILE_PTS(NTYPE)                      &
                                                ! Number of land points which
                                                ! include the nth surface type
     &,                        TILE_INDEX(LAND_PTS,NTYPE)
                                                ! Indices of land points which
                                                ! include the nth surface type
      
      INTEGER :: N,L,C   !local counters: type, land pts, count
!-----------------------------------------------------------------------
! Local parameters
!-----------------------------------------------------------------------
#include "seed.h"

!-----------------------------------------------------------------------
! Create the TILE_INDEX array of land points with each surface type
!-----------------------------------------------------------------------
      TILE_PTS (:) = 0
      TILE_INDEX(:,:) = 0
      
      DO N=1,NTYPE
        C=0
!CDIR NODEP        
        DO L=1,LAND_PTS
          IF ( ((ALL_TILES == 0).AND.                                  &
     &          (FRAC(L,N) >  0.0))                                    &
     &        .OR.                                                     &
     &         ((ALL_TILES == 1).AND.                                  &
     &          ( (N .LT. NTYPE .AND. FRAC(L,NTYPE) .LT. 0.5)          &
     &           .OR.                                                  &
     &            (N .EQ. NTYPE .AND. FRAC(L,NTYPE) .GT. 0.5))) )THEN
            C = C + 1
            TILE_INDEX(C,N) = L
          ENDIF
        ENDDO
        TILE_PTS(N) = C
      ENDDO

      RETURN
      END SUBROUTINE TILEPTS
#endif
