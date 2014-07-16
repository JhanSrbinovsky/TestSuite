
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!*LL  SUBROUTINE SF_EVAP------------------------------------------------
!LL
!LL  Purpose: Calculate surface evaporation and sublimation amounts
!LL           (without applying them to the surface stores).
!LL
!LL
!LL  Suitable for single column usage.
!LL
!LL  Model            Modification history:
!LL version  Date
!LL
!LL   5.2   15/11/00   New Deck         M. Best
!    5.3  25/04/01  Add coastal tiling. Nic Gedney
!LL
!LL  Programming standard: Unified Model Documentation Paper No 4,
!LL                        version 2, dated 18/1/90.
!LL
!LL  Logical component covered: P245.
!LL
!LL  System task:
!LL
!LL  Documentation: UMDP 24
!LL
!LL---------------------------------------------------------------------
!*
!*L Arguments :---------------------------------------------------------
      SUBROUTINE SF_EVAP (                                              &
     & ROW_LENGTH,ROWS,LAND_PTS,NTILES,                                 &
     & LAND_INDEX,TILE_INDEX,TILE_PTS,NSHYD,LTIMER,FLAND,               &
     & ASHTF_TILE,CANOPY,DTRDZ_1,FLAKE,FRACA,SNOW_TILE,RESFS,           &
     & RESFT,RHOKH_1,TILE_FRAC,SMC,WT_EXT_TILE,TIMESTEP,                &
     & FQW_1,FQW_TILE,FTL_1,FTL_TILE,TSTAR_TILE,                        &
     & ECAN,ECAN_TILE,ELAKE_TILE,ESOIL,ESOIL_TILE,EI_TILE,EXT           &
     & )

      IMPLICIT NONE

      INTEGER                                                           &
     & ROW_LENGTH                                                       &
                             ! IN Number of X points?
     &,ROWS                                                             &
                             ! IN Number of Y points?
     &,LAND_PTS                                                         &
                             ! IN Number of land points to be processed.
     &,NTILES                                                           &
                             ! IN Number of tiles per land point.
     &,LAND_INDEX(LAND_PTS)                                             &
                             ! IN Index of land points.
     &,TILE_INDEX(LAND_PTS,NTILES)                                      &
!                            ! IN Index of tile points.
     &,TILE_PTS(NTILES)                                                 &
                             ! IN Number of tile points.
     &,NSHYD                 ! IN Number of soil moisture levels.

      LOGICAL                                                           &
     & LTIMER                ! IN Logical for TIMER.

      REAL                                                              &
     & FLAND(LAND_PTS)                                                  &
                             ! IN Fraction of gridbox which is land.
     &,ASHTF_TILE(LAND_PTS,NTILES)                                      &
!                            ! IN Coefficient to calculate surface
!                            !    heat flux into soil.
     &,CANOPY(LAND_PTS,NTILES)                                          &
!                            ! IN Surface/canopy water on land
!                            !    tiles (kg/m2).
     &,DTRDZ_1(ROW_LENGTH,ROWS)                                         &
!                            ! IN -g.dt/dp for surface layer
     &,FLAKE(LAND_PTS,NTILES)                                           &
                             ! IN Lake fraction.
     &,FRACA(LAND_PTS,NTILES)                                           &
                             ! IN Fraction of surface moisture flux
!                            !    with only aerodynamic resistance
!                            !    for land tiles.
     &,SNOW_TILE(LAND_PTS,NTILES)                                       &
!                            ! IN Lying snow amount on tiles (kg/m2).
     &,RESFS(LAND_PTS,NTILES)                                           &
                             ! IN Combined soil, stomatal and aerodynam.
!                            !    resistance factor for fraction 1-FRACA
!                            !    of land tiles.
     &,RESFT(LAND_PTS,NTILES)                                           &
                             ! IN Total resistance factor
!                            !    FRACA+(1-FRACA)*RESFS.
     &,RHOKH_1(LAND_PTS,NTILES)                                         &
!                            ! IN Surface exchange coefficients.
     &,TILE_FRAC(LAND_PTS,NTILES)                                       &
!                            ! IN Tile fractions.
     &,SMC(LAND_PTS)                                                    &
                             ! IN Available soil moisture (kg/m2).
     &,WT_EXT_TILE(LAND_PTS,NSHYD,NTILES)                               &
!                            ! IN Fraction of transpiration
!                            !    extracted from each soil layer
!                            !    by each tile.
     &,TIMESTEP              ! IN Timestep in seconds.

      REAL                                                              &
     & FQW_1(ROW_LENGTH,ROWS)                                           &
                             ! INOUT Surface moisture flux (kg/m2/s).
     &,FQW_TILE(LAND_PTS,NTILES)                                        &
!                            ! INOUT Local FQW_1 for tiles.
     &,FTL_1(ROW_LENGTH,ROWS)                                           &
                             ! INOUT Surface sensible heat flux (W/m2).
     &,FTL_TILE(LAND_PTS,NTILES)                                        &
!                            ! INOUT Local FTL_1 for tiles.
     &,TSTAR_TILE(LAND_PTS,NTILES)
!                            ! INOUT Tile surface temperatures (K).

      REAL                                                              &
     & ECAN(ROW_LENGTH,ROWS)                                            &
                             ! OUT Gridbox mean evaporation from canopy/
!                            !     surface store (kg per sq m per s).
!                            !     Zero over sea and sea-ice.
     &,ECAN_TILE(LAND_PTS,NTILES)                                       &
!                            ! OUT ECAN for land tiles.
     &,ELAKE_TILE(LAND_PTS,NTILES)                                      &
!                            ! OUT Lake evaporation.
     &,ESOIL(ROW_LENGTH,ROWS)                                           &
                             ! OUT Gridbox mean evapotranspiration from
!                            !     soil moisture (kg per sq m per s).
!                            !     Zero over sea and sea-ice.
     &,ESOIL_TILE(LAND_PTS,NTILES)                                      &
!                            ! OUT ESOIL for land tiles.
     &,EI_TILE(LAND_PTS,NTILES)                                         &
!                            ! OUT Sublimation from snow or land-ice
!                            !     (kg per sq m per s).
     &,EXT(LAND_PTS,NSHYD)   ! OUT Extraction of water from each
!                            !     soil layer (kg/m2/s).


!  External routines called :-
      EXTERNAL TIMER


!  Local and other symbolic constants :-
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
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

      REAL                                                              &
     & DFQW(LAND_PTS)                                                   &
                             ! Increment in GBM moisture flux.
     &,DFTL(LAND_PTS)                                                   &
                             ! Increment in GBM sensible heat flux.
     &,E_TILE_OLD(LAND_PTS,NTILES)                                      &
!                            ! Surface moisture flux before adjustment.
     &,LE_TILE_OLD(LAND_PTS,NTILES)
!                            ! Surf latent heat flux before adjustment.

      REAL                                                              &
     & DIFF_LAT_HTF                                                     &
                             ! Increment in local latent heat flux.
     &,DIFF_SENS_HTF                                                    &
                             ! Increment in local sensible heat flux.
     &,DTSTAR                                                           &
                             ! Increment in local surface temperature.
     &,EDT                   ! Moisture flux x timestep

      INTEGER                                                           &
     & I,J                                                              &
                   ! Loop counter (horizontal field index).
     &,K                                                                &
                   ! Loop counter (land, snow or land-ice field index).
     &,M                                                                &
                   ! Loop counter (soil level index).
     &,L                                                                &
                   ! Loop counter (land point field index).
     &,N           ! Loop counter (tile index).

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SFEVAP  ',3)
      ENDIF

      DO N=1,NTILES
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          E_TILE_OLD(L,N) = FQW_TILE(L,N)
          IF (SNOW_TILE(L,N)  >   0.) THEN
            LE_TILE_OLD(L,N) = (LC + LF)*FQW_TILE(L,N)
          ELSE
            LE_TILE_OLD(L,N) = LC*FQW_TILE(L,N)
          ENDIF
        ENDDO
      ENDDO

      DO N=1,NTILES
        DO L=1,LAND_PTS
          ECAN_TILE(L,N) = 0.
          ESOIL_TILE(L,N) = 0.
          ELAKE_TILE(L,N) = 0.
          EI_TILE(L,N) = 0.
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Sublimation from snow-covered land tiles
!-----------------------------------------------------------------------
      DO N=1,NTILES
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          IF (SNOW_TILE(L,N)  >   0.) THEN
            EI_TILE(L,N) =  FQW_TILE(L,N)
            EDT = EI_TILE(L,N)*TIMESTEP
            IF ( EDT  >   SNOW_TILE(L,N) )                              &
     &        EI_TILE(L,N) = SNOW_TILE(L,N) / TIMESTEP
            FQW_TILE(L,N) = FQW_TILE(L,N) -  EI_TILE(L,N)
          ENDIF
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Surface evaporation from and condensation onto snow-free land
!-----------------------------------------------------------------------
      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        ECAN(I,J) = 0.
        ESOIL(I,J) = 0.
       ENDDO
      ENDDO

      DO N=1,NTILES
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          IF ( FQW_TILE(L,N)  >   0.0 ) THEN
            ECAN_TILE(L,N) = (1. - FLAKE(L,N)) *                        &
     &                       FRACA(L,N) * FQW_TILE(L,N) / RESFT(L,N)
            ESOIL_TILE(L,N) = (1. - FLAKE(L,N)) *                       &
     &                        (1. - FRACA(L,N))*RESFS(L,N)*FQW_TILE(L,N)&
     &                                                      / RESFT(L,N)
            ELAKE_TILE(L,N) = FLAKE(L,N)*FQW_TILE(L,N) / RESFT(L,N)
            EDT = ECAN_TILE(L,N)*TIMESTEP
            IF ( EDT  >   CANOPY(L,N) ) THEN
              ESOIL_TILE(L,N) =  (1. - FLAKE(L,N)) *                    &
     &                           (1. - FRACA(L,N)*CANOPY(L,N)/EDT) *    &
     &                               RESFS(L,N)*FQW_TILE(L,N)/RESFT(L,N)
              ECAN_TILE(L,N) = CANOPY(L,N) / TIMESTEP
            ENDIF
          ELSEIF (SNOW_TILE(L,N) <= 0.) THEN
            IF (TSTAR_TILE(L,N) >= TM) THEN
              ECAN_TILE(L,N) = (1. - FLAKE(L,N))*FQW_TILE(L,N)
              ELAKE_TILE(L,N) = FLAKE(L,N)*FQW_TILE(L,N)
            ELSE
              EI_TILE(L,N) =  FQW_TILE(L,N)
            ENDIF
          ENDIF
          ECAN(I,J) = ECAN(I,J) + TILE_FRAC(L,N)*ECAN_TILE(L,N)
          ESOIL(I,J) = ESOIL(I,J) + TILE_FRAC(L,N)*ESOIL_TILE(L,N)
        ENDDO
      ENDDO

!      print *,'evaptest1',FQW_1,FQW_TILE(:,1),FQW_TILE(:,8),FQW_TILE
!      print *,'evaptest2',ECAN_TILE(:,1),ECAN_TILE(:,8)
!      print *,'evaptest3',ESOIL_TILE(:,1),ESOIL_TILE(:,8)
!      print *,'evaptest4',FRACA(1,1),FRACA(1,8)
!      print *,'evaptest5',RESFT(1,1),RESFT(1,8)
!      print *,'evaptest6',RESFS(1,1),RESFs(1,8)
!      print *,'evaptest71',FRACA(1,1)/ RESFT(1,1)
!      print *,'evaptest72',FRACA(1,8)/ RESFT(1,8)
!      print *,'evaptest91',(1.-FRACA(1,1))*RESFS(1,1)                  &
!     &                                               / RESFT(1,1)
!      print *,'evaptest92',(1.-FRACA(1,8))*RESFS(1,8)                  &
!     &                                               / RESFT(1,8)
!-----------------------------------------------------------------------
! Soil evapotranspiration
!-----------------------------------------------------------------------
      DO L=1,LAND_PTS
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
        EDT = ESOIL(I,J)*TIMESTEP
        IF ( EDT  >   SMC(L) ) THEN
          DO N=1,NTILES
            ESOIL_TILE(L,N) = SMC(L)*ESOIL_TILE(L,N) / EDT
          ENDDO
          ESOIL(I,J) = SMC(L) / TIMESTEP
        ENDIF
      ENDDO

      DO M=1,NSHYD
        DO L=1,LAND_PTS
          EXT(L,M) = 0.
        ENDDO
      ENDDO

      DO M=1,NSHYD
        DO N=1,NTILES
          DO K=1,TILE_PTS(N)
            L = TILE_INDEX(K,N)
            EXT(L,M) = EXT(L,M) + TILE_FRAC(L,N)*WT_EXT_TILE(L,M,N)     &
     &                                          *ESOIL_TILE(L,N)
          ENDDO
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Calculate increments to surface heat fluxes, moisture fluxes and
! temperatures
!-----------------------------------------------------------------------
      DO L=1,LAND_PTS
        DFTL(L) = 0.
        DFQW(L) = 0.
      ENDDO

      DO N=1,NTILES
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          DIFF_LAT_HTF = (LC + LF)*EI_TILE(L,N) + LC*ECAN_TILE(L,N)     &
     &                    + LC*ESOIL_TILE(L,N) + LC*ELAKE_TILE(L,N)     &
     &                    - LE_TILE_OLD(L,N)
          DIFF_SENS_HTF = - DIFF_LAT_HTF /                              &
     &                        ( 1. + ASHTF_TILE(L,N)/(CP*RHOKH_1(L,N)) )
          FTL_TILE(L,N) = FTL_TILE(L,N) + DIFF_SENS_HTF
          DTSTAR = - (DIFF_LAT_HTF + DIFF_SENS_HTF) / ASHTF_TILE(L,N)
          TSTAR_TILE(L,N) = TSTAR_TILE(L,N) + DTSTAR
          DFTL(L) = DFTL(L) + TILE_FRAC(L,N)*DIFF_SENS_HTF
          DFQW(L) = DFQW(L) + TILE_FRAC(L,N)*( ECAN_TILE(L,N) +         &
     &                  ESOIL_TILE(L,N) + EI_TILE(L,N) + ELAKE_TILE(L,N)&
     &                  - E_TILE_OLD(L,N) )
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Update level 1 temperature and humidity and GBM heat and moisture
! fluxes due to limited moisture availability
!-----------------------------------------------------------------------
      DO L=1,LAND_PTS
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
        FTL_1(I,J) = FTL_1(I,J) + FLAND(L)*DFTL(L)
        FQW_1(I,J) = FQW_1(I,J) + FLAND(L)*DFQW(L)
      ENDDO

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SFEVAP  ',4)
      ENDIF

      RETURN
      END SUBROUTINE SF_EVAP
