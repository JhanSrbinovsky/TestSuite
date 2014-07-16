
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!LL  SUBROUTINE SURF_HYD-----------------------------------------------
!LL
!LL  PURPOSE : TO CARRY OUT CANOPY AND SURFACE HYDROLOGY CALCULATIONS
!LL
!LL            CANOPY WATER CONTENT IS DEPRECIATED BY EVAPORATION
!LL
!LL            SNOWMELT IS RUNOFF THE SURFACE WITHOUT INTERACTING
!LL            WITH THE CANOPY
!LL
!LL            THE CANOPY INTERCEPTION AND SURFACE RUNOFF OF
!LL            LARGE-SCALE RAIN IS CALCUALTED
!LL
!LL            THE CANOPY INTERCEPTION AND SURFACE RUNOFF OF
!LL            CONVECTIVE RAIN IS CALCUALTED
!LL
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  WRITTEN FOR CRAY-YMP BY S.ALLEN-BETT AND D.GREGORY
!LL
!LL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 5.2:
!LL VERSION  DATE
!LL  5.2   15/11/00   New Deck         M. Best
!  5.5  12/02/03   Code added for large-scale hydrology.
!                                          Nic Gedney.
!  6.1  17/08/04  Add SSFM code                           Ian Pearman
!  6.2  02/02/06  Remove SSFM. P.Selwood
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
!LL  VERSION NO. 1 18/1/90
!LL
!LL  SYSTEM TASK : P252
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER NO 25
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE SURF_HYD (NPNTS,NTILES,TILE_PTS,TILE_INDEX,            &
     &                     CAN_CPY,E_CANOPY,FRAC,INFIL,CON_RAIN,        &
     &                     LS_RAIN,MELT_TILE,SNOW_MELT,TIMESTEP,        &
     &                     CAN_WCNT,                                    &
     &                     CAN_WCNT_GB,DSMC_DT,                         &
     &                     L_TOP,L_PDM,NSHYD,SOIL_PTS,SOIL_INDEX,       &
     &                     SURF_ROFF,TOT_TFALL,DUN_ROFF,FSAT,           &
     &                     V_SAT,STHU,STHF)

      IMPLICIT NONE

      INTEGER                                                           &
     & NPNTS                                                            &
                            ! IN Total number of land points.
     &,NTILES                                                           &
                            ! IN Number of tiles.
     &,NSHYD                                                            &
                           ! IN Number of soil moisture levels.
     &,SOIL_PTS                                                         &
                           ! IN Number of soil points.
     &,TILE_PTS(NTILES)                                                 &
                            ! IN Number of tile points.
     &,TILE_INDEX(NPNTS,NTILES)                                         &
!                           ! IN Index of tile points.
     &,SOIL_INDEX(NPNTS)   ! IN Array of soil points.

      REAL                                                              &
     & CAN_CPY(NPNTS,NTILES)                                            &
                            ! IN Canopy capacity for land tiles (kg/m2).
     &,E_CANOPY(NPNTS,NTILES)                                           &
                             !IN Canopy evaporation (kg/m2/s).
     &,FRAC(NPNTS,NTILES)                                               &
                            ! IN Tile fractions.
     &,INFIL(NPNTS,NTILES)                                              &
                            ! IN Infiltration rate (kg/m2/s).
     &,CON_RAIN(NPNTS)                                                  &
                            ! IN Convective rain (kg/m2/s).
     &,LS_RAIN(NPNTS)                                                   &
                            ! IN Large-scale rain (kg/m2/s).
     &,MELT_TILE(NPNTS,NTILES)                                          &
!                           ! IN Snow melt on tiles (kg/m2/s).
     &,SNOW_MELT(NPNTS)                                                 &
                            ! IN GBM snow melt (kg/m2/s).
     &,FSAT(NPNTS)                                                      &
                           ! IN Surface saturation fraction.
     &,STHF(NPNTS,NSHYD)                                                &
                           ! INOUT Frozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
     &,STHU(NPNTS,NSHYD)                                                &
                           ! INOUT Unfrozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
     &,V_SAT(NPNTS)                                                     &
                           ! IN Volumetric soil moisture
!                          !    concentration at saturation
!                          !    (m3 H2O/m3 soil).

     &,TIMESTEP             ! IN Timestep (s).


      LOGICAL                                                           &
     & L_TOP                                                            &
                            ! IN TOPMODEL-based hydrology logical.
     &,L_PDM                ! IN PDM logical.
      REAL                                                              &
     & CAN_WCNT(NPNTS,NTILES)!INOUT Tile canopy water contents (kg/m2).

      REAL                                                              &
     & CAN_WCNT_GB(NPNTS)                                               &
                            ! OUT Gridbox canopy water content (kg/m2).
     &,DSMC_DT(NPNTS)                                                   &
                            ! OUT Rate of change of soil moisture
                            !     content (kg/m2/s).
     &,SURF_ROFF(NPNTS)                                                 &
                            ! OUT Cumulative surface runoff (kg/m2/s).
     &,TOT_TFALL(NPNTS)                                                 &
                            ! OUT Cumulative canopy throughfall
                            !     (kg/m2/s).
     &,DUN_ROFF(NPNTS)      ! OUT Cumulative Dunne sfc runoff (kg/m2/s).

!---------------------------------------------------------------------
!  External routines called :-

      EXTERNAL FRUNOFF,SIEVE,PDM


!-----------------------------------------------------------------------

!  Workspace -----------------------------------------------------------
      REAL                                                              &
     & CAN_COND(NPNTS)      ! Canopy condensation (kg/m2/s).

      INTEGER                                                           &
     & I                                                                &
                            ! Loop counter (land points).
     &,J                                                                &
                            ! Loop counter (tile points).
     &,N                    ! Loop counter (tiles).


! Zero cumulative stores
      DO I=1,NPNTS
        CAN_WCNT_GB(I) = 0.0
        TOT_TFALL(I) = 0.0
        SURF_ROFF(I) = 0.0
        DSMC_DT(I)   = 0.0
        DUN_ROFF(I) = 0.0
      ENDDO

! Reduce canopy water content by evaporation
      DO N=1,NTILES
        DO J=1,TILE_PTS(N)
        I = TILE_INDEX(J,N)
          IF (E_CANOPY(I,N)  >   0.0)                                   &
     &      CAN_WCNT(I,N) =                                             &
     &                 MAX( CAN_WCNT(I,N) - E_CANOPY(I,N)*TIMESTEP, 0. )
        ENDDO
      ENDDO

      DO N=1,NTILES

! Surface runoff of snowmelt, assumed to cover 100% of tile
! DEPENDS ON: frunoff
        CALL FRUNOFF (NPNTS,TILE_PTS(N),TILE_INDEX(1,N),1.,             &
     &                CAN_CPY(1,N),CAN_CPY(1,N),INFIL(1,N),             &
     &                MELT_TILE(1,N),FRAC(1,N),TIMESTEP,                &
     &                SURF_ROFF)

! Define canopy condensation when evaporation is negative
        DO J=1,TILE_PTS(N)
        I = TILE_INDEX(J,N)
          IF ( E_CANOPY(I,N)  <   0. ) THEN
           CAN_COND(I) = - E_CANOPY(I,N)
          ELSE
           CAN_COND(I) = 0.
          ENDIF
        ENDDO

! Canopy interception, throughfall and surface runoff for condensation,
! assumed to cover 100% of gridbox
! DEPENDS ON: sieve
        CALL SIEVE (NPNTS,TILE_PTS(N),TILE_INDEX(1,N),1.,               &
     &              CAN_CPY(1,N),CAN_COND,FRAC(1,N),TIMESTEP,           &
     &              CAN_WCNT(1,N),TOT_TFALL)
! DEPENDS ON: frunoff
        CALL FRUNOFF (NPNTS,TILE_PTS(N),TILE_INDEX(1,N),1.,             &
     &                CAN_CPY(1,N),CAN_WCNT(1,N),INFIL(1,N),CAN_COND,   &
     &                FRAC(1,N),TIMESTEP,                               &
     &                SURF_ROFF)

! Canopy interception, throughfall and surface runoff for large-scale
! rain, assumed to cover 100% of gridbox
! DEPENDS ON: sieve
        CALL SIEVE (NPNTS,TILE_PTS(N),TILE_INDEX(1,N),1.,               &
     &              CAN_CPY(1,N),LS_RAIN,FRAC(1,N),TIMESTEP,            &
     &              CAN_WCNT(1,N),TOT_TFALL)
! DEPENDS ON: frunoff
        CALL FRUNOFF (NPNTS,TILE_PTS(N),TILE_INDEX(1,N),1.,             &
     &                CAN_CPY(1,N),CAN_WCNT(1,N),INFIL(1,N),LS_RAIN,    &
     &                FRAC(1,N),TIMESTEP,                               &
     &                SURF_ROFF)

! Canopy interception, throughfall and surface runoff for convective
! rain, assumed to cover 30% of gridbox
! DEPENDS ON: sieve
        CALL SIEVE (NPNTS,TILE_PTS(N),TILE_INDEX(1,N),0.3,              &
     &              CAN_CPY(1,N),CON_RAIN,FRAC(1,N),TIMESTEP,           &
     &              CAN_WCNT(1,N),TOT_TFALL)
! DEPENDS ON: frunoff
        CALL FRUNOFF (NPNTS,TILE_PTS(N),TILE_INDEX(1,N),0.3,            &
     &                CAN_CPY(1,N),CAN_WCNT(1,N),INFIL(1,N),CON_RAIN,   &
     &                FRAC(1,N),TIMESTEP,                               &
     &                SURF_ROFF)

        DO J=1,TILE_PTS(N)
          I = TILE_INDEX(J,N)
          CAN_WCNT_GB(I) = CAN_WCNT_GB(I) + FRAC(I,N)*CAN_WCNT(I,N)
        ENDDO

      ENDDO

!-----------------------------------------------------------------------
! Calculate Saturation excess runoff through PDM:
!-----------------------------------------------------------------------

      IF(SOIL_PTS >  0.AND.L_PDM)                                       &
! DEPENDS ON: pdm
     &  CALL PDM(                                                       &
     &           NPNTS,SOIL_PTS,SOIL_INDEX,NSHYD,                       &
     &          TOT_TFALL,SNOW_MELT,SURF_ROFF,TIMESTEP,                 &
     &          V_SAT,DUN_ROFF,STHU,STHF)

      IF(L_TOP)THEN
        DO I=1,NPNTS
          DUN_ROFF(I)=FSAT(I)*(TOT_TFALL(I) + SNOW_MELT(I)              &
     &                  - SURF_ROFF(I))
        ENDDO
      ENDIF
      DO I=1,NPNTS
        IF(L_TOP.OR.L_PDM)SURF_ROFF(I) = SURF_ROFF(I) + DUN_ROFF(I)
        DSMC_DT(I) = TOT_TFALL(I) + SNOW_MELT(I) - SURF_ROFF(I)
      ENDDO

      RETURN
      END SUBROUTINE SURF_HYD
