
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!!!  SUBROUTINES SF_RIB_LAND and SF_RIB_SEA ---------------------------
!!!
!!!  Purpose: Calculate bulk Richardson number for surface layer
!!!
!!!
!!!    Model            Modification history
!!!   version  date
!!!    5.2   15/11/00   New Deck         M. Best
!    5.3  25/04/01  Add coastal tiling. Nic Gedney
!!!
!!!    Programming standard:
!!!
!!!  ------------------------------------------------------------------

!    SUBROUTINE SF_RIB_LAND--------------------------------------------
!
!    Calculate RIB for land tiles
!
!    ------------------------------------------------------------------
      SUBROUTINE SF_RIB_LAND (                                          &
     & ROW_LENGTH,ROWS,LAND_PTS,TILE_PTS,LAND_INDEX,TILE_INDEX,         &
     & BQ_1,BT_1,QSTAR,QW_1,RESFT,TL_1,TSTAR,VSHR,Z0H,Z0M,Z1_TQ,Z1_UV,  &
     & RIB,DB,LTIMER                                                    &
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
     &,TILE_INDEX(LAND_PTS) ! IN Index of tile points.

      LOGICAL                                                           &
     & LTIMER              ! IN logical for TIMER

      REAL                                                              &
     & BQ_1(ROW_LENGTH,ROWS)                                            &
                             ! IN A buoyancy parameter for lowest atm
!                            !    level. ("beta-q twiddle").
     &,BT_1(ROW_LENGTH,ROWS)                                            &
                             ! IN A buoyancy parameter for lowest atm
!                            !    level. ("beta-T twiddle").
     &,QSTAR(LAND_PTS)                                                  &
                             ! IN Surface saturated sp humidity.
     &,QW_1(ROW_LENGTH,ROWS)                                            &
                             ! IN Total water content of lowest
!                            !    atmospheric layer (kg per kg air).
     &,RESFT(LAND_PTS)                                                  &
                             ! IN Total resistance factor.
     &,TL_1(ROW_LENGTH,ROWS)                                            &
                             ! IN Liquid/frozen water temperature for
!                            !    lowest atmospheric layer (K).
     &,TSTAR(LAND_PTS)                                                  &
                             ! IN Surface temperature (K).
     &,VSHR(ROW_LENGTH,ROWS)                                            &
                             ! IN Magnitude of surface-to-lowest-level
!                            !    wind shear.
     &,Z0H(LAND_PTS)                                                    &
                             ! IN Roughness length for heat and
!                            !    moisture m
     &,Z0M(LAND_PTS)                                                    &
                             ! IN Effective roughness length for
!                            !    momentum
     &,Z1_TQ(ROW_LENGTH,ROWS)                                           &
                             ! IN Height of lowest TQ level (m).
     &,Z1_UV(ROW_LENGTH,ROWS)! IN Height of lowest UV level (m).

      REAL                                                              &
     & RIB(LAND_PTS)                                                    &
                           ! OUT Bulk Richardson number for lowest layer
     &,DB(LAND_PTS)        ! OUT Buoyancy difference between surface
!                          !     and lowest atmospheric level.


!  External routines called :-
      EXTERNAL TIMER


!  Symbolic constants -----------------------------------------------

!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------

!  Workspace --------------------------------------------------------
      INTEGER                                                           &
     & I,J                                                              &
                           ! Horizontal field index.
     &,K                                                                &
                           ! Tile field index.
     &,L                   ! Land field index.

      REAL                                                              &
     & DQ(LAND_PTS)                                                     &
                           ! Sp humidity difference between surface
!                          ! and lowest atmospheric level (Q1 - Q*).
     &,DTEMP(LAND_PTS)     ! Modified temperature difference between
!                            surface and lowest atmospheric level.

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SF_RIB  ',3)
      ENDIF

!-----------------------------------------------------------------------
!!  1 Calculate temperature (strictly, liquid/ice static energy) and
!!    humidity jumps across the surface layer.
!-----------------------------------------------------------------------
      DO K=1,TILE_PTS
        L = TILE_INDEX(K)
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
        DTEMP(L) = TL_1(I,J) - TSTAR(L) +                               &
     &                  (G/CP)*(Z1_TQ(I,J)+Z0M(L)-Z0H(L))
!                                                             ! P243.118
        DQ(L) = QW_1(I,J) - QSTAR(L)                          ! P243.119
      ENDDO

!-----------------------------------------------------------------------
!!  2 Calculate bulk Richardson numbers for the surface layer.
!-----------------------------------------------------------------------
      DO K=1,TILE_PTS
        L = TILE_INDEX(K)
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
        DB(L) = G*(BT_1(I,J)*DTEMP(L) + BQ_1(I,J)*RESFT(L)*DQ(L))
        RIB(L) = Z1_UV(I,J)*DB(L) / ( VSHR(I,J)*VSHR(I,J) )
      ENDDO

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SF_RIB  ',4)
      ENDIF

      RETURN
      END SUBROUTINE SF_RIB_LAND

!    SUBROUTINE SF_RIB_SEA---------------------------------------------
!
!    Calculate RIB for sea, sea-ice and sea-ice leads
!
!    ------------------------------------------------------------------
