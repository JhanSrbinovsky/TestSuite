
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!*LL  SUBROUTINE SF_RESIST----------------------------------------------
!LL
!LL  Purpose: Calculate surface moisture flux resistance factors.
!LL
!LL
!LL    Model            Modification history
!LL   version  date
!LL    5.2   15/11/00   New Deck         M. Best
!  6.1  08/12/03  Add !CDIR NODEP to force vectorisation. R Barnes
!LL
!LL    Programming standard:
!LL
!LL
!LLEND-----------------------------------------------------------------
!*
!*L  Arguments --------------------------------------------------------
      SUBROUTINE SF_RESIST (                                            &
     & ROW_LENGTH,ROWS,LAND_PTS,TILE_PTS,LAND_INDEX,TILE_INDEX,         &
     & CANOPY,CATCH,CH,DQ,EPDT,FLAKE,GC,SNOW,VSHR,                      &
     & FRACA,RESFS,RESFT,LTIMER                                         &
     & )

! module for land-surface namelist
      USE LAND_SURF_MOD, ONLY :                                         &
     & FRAC_SNOW_SUBL_MELT                                              &
     &,MASKD

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
     & LTIMER                ! IN Logical switch for TIMER diags

      REAL                                                              &
     & CANOPY(LAND_PTS)                                                 &
                           ! IN Surface water (kg per sq metre).  F642.
     &,CATCH(LAND_PTS)                                                  &
                           ! IN Surface capacity (max. surface water)
!                          !    (kg per sq metre).  F6416.
     &,CH(LAND_PTS)                                                     &
                           ! IN Transport coefficient for heat and
!                          !    moisture transport
     &,DQ(LAND_PTS)                                                     &
                           ! IN Sp humidity difference between surface
!                          !    and lowest atmospheric level (Q1 - Q*).
     &,EPDT(LAND_PTS)                                                   &
                           ! IN "Potential" Evaporation * Timestep.
!                          !    Dummy variable for first call to routine
     &,FLAKE(LAND_PTS)                                                  &
                           ! IN Lake fraction.
     &,GC(LAND_PTS)                                                     &
                           ! IN Interactive canopy conductance
!                          !    to evaporation (m/s)
     &,SNOW(LAND_PTS)                                                   &
                           ! IN Lying snow amount (kg per sq metre).
     &,VSHR(ROW_LENGTH,ROWS)! IN Magnitude of surface-to-lowest-level
!                          !     windshear

      REAL                                                              &
     & FRACA(LAND_PTS)                                                  &
                           ! OUT Fraction of surface moisture flux with
!                          !     only aerodynamic resistance.
     &,RESFS(LAND_PTS)                                                  &
                           ! OUT Combined soil, stomatal and aerodynamic
!                          !     resistance factor for fraction 1-FRACA.
     &,RESFT(LAND_PTS)     ! OUT Total resistance factor
!                          !     FRACA+(1-FRACA)*RESFS.


!  External routines called :-
      EXTERNAL TIMER


! Workspace -----------------------------------------------------------
      INTEGER                                                           &
     & I,J                                                              &
                   ! Horizontal field index.
     &,K                                                                &
                   ! Tile field index.
     &,L           ! Land field index.

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SFRESIST',3)
      ENDIF

!-----------------------------------------------------------------------
!     Evaporation over land surfaces without snow is limited by
!     soil moisture availability and stomatal resistance.
!     Set FRACA (= fA in the documentation) according to P243.68,
!     and RESFS (= fS) according to P243.75 and P243.61.
!-----------------------------------------------------------------------
!CDIR NODEP
      DO K=1,TILE_PTS
        L = TILE_INDEX(K)
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH

!-----------------------------------------------------------------------
! Calculate the fraction of the flux with only aerodynamic resistance
! (canopy evaporation).
! Set to 1 for negative moisture flux or snow-covered land
! (no surface/stomatal resistance to condensation).
!-----------------------------------------------------------------------
        FRACA(L) = 1.0
        IF (DQ(L) <  0. .AND. SNOW(L) <= 0.) FRACA(L) = 0.0
        IF (DQ(L) <  0. .AND. SNOW(L) <= 0. .AND. CATCH(L) >  0.)       &
     &    FRACA(L) = CANOPY(L) / ( EPDT(L) + CATCH(L) )
        IF (SNOW(L) > 0.) THEN
          IF (FRAC_SNOW_SUBL_MELT == 1) THEN
            FRACA(L) = 1.0 - EXP(-MASKD*SNOW(L))
          ENDIF
        ENDIF
        FRACA(L) = MIN(FRACA(L),1.0)

!-----------------------------------------------------------------------
! Calculate resistance factors for transpiration from vegetation tiles
! and bare soil evaporation from soil tiles.
!-----------------------------------------------------------------------
        RESFS(L) = GC(L) / ( GC(L) + CH(L)*VSHR(I,J) )
        RESFT(L) = FLAKE(L) + (1. - FLAKE(L)) *                         &
     &                        ( FRACA(L) + (1. - FRACA(L))*RESFS(L) )

      ENDDO

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SFRESIST',4)
      ENDIF

      RETURN
      END SUBROUTINE SF_RESIST
