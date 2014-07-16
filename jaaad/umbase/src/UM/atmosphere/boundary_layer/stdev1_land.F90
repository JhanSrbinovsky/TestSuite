#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!!!  SUBROUTINES STDEV1_SEA and STDEV1_LAND ----------------------------
!!!
!!!  Purpose: Calculate the standard deviations of layer 1 turbulent
!!!           fluctuations of temperature and humidity using approximate
!!!           formulae from first order closure.
!!!
!!!    Model            Modification history
!!!   version  date
!!!    5.2   15/11/00   New Deck         M. Best
!    5.3  25/04/01  Add coastal tiling. Nic Gedney
!!!
!!!    Programming standard:
!!!
!!!  -------------------------------------------------------------------
!

!!!  SUBROUTINE STDEV1_SEA ---------------------------------------------
!!!  Layer 1 standard deviations for sea and sea-ice
!!!  -------------------------------------------------------------------

!!!  SUBROUTINE STDEV1_LAND --------------------------------------------
!!!  Layer 1 standard deviations for land tiles
!!!  -------------------------------------------------------------------
      SUBROUTINE STDEV1_LAND (                                          &
     & ROW_LENGTH,ROWS,LAND_PTS,TILE_PTS,LAND_INDEX,TILE_INDEX,FLAND,   &
     & BQ_1,BT_1,FQW_1,FTL_1,RHOKM_1,RHOSTAR,VSHR,Z0M,Z1_TQ,TILE_FRAC,  &
     & Q1_SD,T1_SD,LTIMER                                               &
     & )

      IMPLICIT NONE

      INTEGER                                                           &
     & ROW_LENGTH                                                       &
                             ! IN Number of X points?
     &,ROWS                                                             &
                             ! IN Number of Y points?
     &,LAND_PTS                                                         &
                             ! IN Total number of land points.
     &,TILE_PTS                                                         &
                             ! IN Number of tile points.
     &,LAND_INDEX(LAND_PTS)                                             &
                             ! IN Index of land points.
     &,TILE_INDEX(LAND_PTS)  ! IN Index of tile points.

      LOGICAL                                                           &
     & LTIMER                ! IN logical for TIMER

      REAL                                                              &
     & FLAND(LAND_PTS)                                                  &
     &,BQ_1(ROW_LENGTH,ROWS)                                            &
                             ! IN Buoyancy parameter.
     &,BT_1(ROW_LENGTH,ROWS)                                            &
                             ! IN Buoyancy parameter.
     &,FQW_1(LAND_PTS)                                                  &
                             ! IN Surface flux of QW.
     &,FTL_1(LAND_PTS)                                                  &
                             ! IN Surface flux of TL.
     &,RHOKM_1(LAND_PTS)                                                &
                             ! IN Surface momentum exchange coefficient.
     &,RHOSTAR(ROW_LENGTH,ROWS)                                         &
                               ! IN Surface air density.
     &,VSHR(ROW_LENGTH,ROWS)                                            &
                             ! IN Magnitude of surface-to-lowest-level
!                            !    wind shear.
     &,Z0M(LAND_PTS)                                                    &
                             ! IN Roughness length for momentum.
     &,Z1_TQ(ROW_LENGTH,ROWS)                                           &
                             ! IN Height of lowest tq level.
     &,TILE_FRAC(LAND_PTS)   ! IN Tile fraction.

      REAL                                                              &
     & Q1_SD(ROW_LENGTH,ROWS)                                           &
                             ! INOUT Standard deviation of turbulent
!                            !       fluctuations of surface layer
!                            !       specific humidity (kg/kg).
     &,T1_SD(ROW_LENGTH,ROWS)! INOUT Standard deviation of turbulent
!                            !       fluctuations of surface layer
!                            !       temperature (K).


!  External routines called :-
      EXTERNAL TIMER

#include "c_g.h"

!  Workspace --------------------------------------------------------
      INTEGER                                                           &
     & I,J                                                              &
                             ! Horizontal field index.
     &,K                                                                &
                             ! Tile index.
     &,L                     ! Land field inde.
      REAL                                                              &
     & VS                                                               &
                             ! Surface layer friction velocity
     &,VSF1_CUBED                                                       &
                             ! Cube of surface layer free convective
!                            ! scaling velocity
     &,WS1                   ! Turbulent velocity scale for surface
!                            ! layer

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('STDEV1  ',3)
      ENDIF

      DO K=1,TILE_PTS
        L = TILE_INDEX(K)
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH

        VS = SQRT ( RHOKM_1(L)/RHOSTAR(I,J) * VSHR(I,J) )
        VSF1_CUBED = 1.25*G*(Z1_TQ(I,J) + Z0M(L)) *                     &
     &             ( BT_1(I,J)*FTL_1(L) + BQ_1(I,J)*FQW_1(L) ) /        &
     &                 RHOSTAR(I,J)
        IF ( VSF1_CUBED  >   0.0 ) THEN
          WS1 = ( VSF1_CUBED + VS*VS*VS ) ** (1.0/3.0)
          T1_SD(I,J) = T1_SD(I,J) + MAX ( 0.0 ,                         &
     &      FLAND(L)*TILE_FRAC(L)*1.93*FTL_1(L) / (RHOSTAR(I,J)*WS1) )
          Q1_SD(I,J) = Q1_SD(I,J) + MAX ( 0.0 ,                         &
     &      FLAND(L)*TILE_FRAC(L)*1.93*FQW_1(L) / (RHOSTAR(I,J)*WS1) )
        ENDIF

      ENDDO

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('STDEV1  ',4)
      ENDIF

      RETURN
      END SUBROUTINE STDEV1_LAND
#endif
