
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!*LL  SUBROUTINE RES_TILE ----------------------------------------------
!LL
!LL  Purpose: Calculate gridbox-mean resistance factor used by TR_MIX
!LL           for 7A tiled land surface.
!LL
!LL  Model           Modification history:
!LL version  Date
!LL  5.2   15/11/00   New Deck         M. Best
!LL
!*----------------------------------------------------------------------
      SUBROUTINE RES_TILE (ROW_LENGTH,ROWS,LAND_PTS,LAND_INDEX,         &
     &                     NTILES,TILE_INDEX,TILE_PTS,SOLUBLE,          &
     &                     ARESIST,ARESIST_TILE,CANOPY,CATCH,GS_TILE,   &
     &                     RESIST_B_TILE,SNOW_TILE,TILE_FRAC,RESB,RESS, &
     &                     RES_FACTOR)

      IMPLICIT NONE

      INTEGER                                                           &
     & ROW_LENGTH                                                       &
                             ! IN Number of X points?
     &,ROWS                                                             &
                             ! IN Number of Y points?
     &,LAND_PTS                                                         &
                             ! IN Total number of land points.
     &,LAND_INDEX(LAND_PTS)                                             &
                             ! IN Index of land points.
     &,NTILES                                                           &
                             ! IN Number of land tiles.
     &,TILE_INDEX(LAND_PTS,NTILES)                                      &
!                            ! IN Index of tile points.
     &,TILE_PTS(NTILES)      ! IN Number of tile points.

      LOGICAL                                                           &
     & SOLUBLE               ! IN .TRUE. for soluble aerosols.

      REAL                                                              &
     & ARESIST(ROW_LENGTH,ROWS)                                         &
                                ! IN GBM aerodynamic resistance (s/m).
     &,ARESIST_TILE(LAND_PTS,NTILES)                                    &
!                               ! IN 1/(CD_STD*VSHR) on land tiles.
     &,CANOPY(LAND_PTS,NTILES)                                          &
!                               ! IN Surface water on land tiles (kg/m2)
     &,CATCH(LAND_PTS,NTILES)                                           &
!                               ! IN Surface capacity (max. surface
!                               !    water) of land tiles (kg/m2).
     &,GS_TILE(LAND_PTS,NTILES)                                         &
!                               ! IN Surface conductance for land tiles.
     &,RESIST_B_TILE(LAND_PTS,NTILES)                                   &
!                               ! IN (1/CH-1/CD_STD)/VSHR on land tiles.
     &,SNOW_TILE(LAND_PTS,NTILES)                                       &
!                               ! IN Snow mass on land tiles (kg/m2).
     &,TILE_FRAC(LAND_PTS,NTILES)                                       &
                                ! IN Tile fractions.
     &,RESB                                                             &
                                ! IN Rb(aerosol) / Rb(H2O).
     &,RESS                                                             &
                                ! IN Rs(aerosol) / Rs(H2O).
     &,RES_FACTOR(ROW_LENGTH,ROWS)! OUT Ra/(Ra+Rb+Rs) for dry deposition

      REAL                                                              &
     & DAMP_FACTOR(LAND_PTS,NTILES)                                     &
!                               ! Canopy moistening factor
     &,RS_TILE(LAND_PTS,NTILES)                                         &
!                               ! Surface reistance for land tiles.
     &,STR_RESIST_B                                                     &
                                ! Rb for aerosol.
     &,STR_RESIST_S                                                     &
                                ! Rs for aerosol.
     &,ASNOW                                                            &
                                ! Parameter for snow fraction
!                               ! calculation.
     &,COND_LIM                                                         &
                                ! Low limit for canopy conductance.
     &,R_SNOW                                                           &
                                ! Resistance to dry deposition over snow
     &,SNOW_F                   ! Snow cover fraction.
      PARAMETER (ASNOW=0.2, COND_LIM=1.0E-3, R_SNOW=1.0E3)

      INTEGER                                                           &
     & I,J                                                              &
                   ! Loop counter (horizontal field index).
     &,K                                                                &
                   ! Loop counter (tile field index).
     &,L                                                                &
                   ! Loop counter (land point field index).
     &,N           ! Loop counter (tile index).


      DO N=1,NTILES
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          DAMP_FACTOR(L,N) = 1.0
          IF (GS_TILE(L,N)  >   COND_LIM) THEN
            RS_TILE(L,N) = 1. / GS_TILE(L,N)
          ELSE
            RS_TILE(L,N) = 1. / COND_LIM
          ENDIF
        ENDDO
      ENDDO

      IF (SOLUBLE) THEN
        DO N=1,NTILES
          DO K=1,TILE_PTS(N)
            L = TILE_INDEX(K,N)
            IF( (CATCH(L,N)  >   0.01) .AND.                            &
     &          (CANOPY(L,N)  >   0.0) ) THEN
              IF( CANOPY(L,N)  <=  CATCH(L,N) ) THEN
                DAMP_FACTOR(L,N) = 1. - 0.66667*CANOPY(L,N)/CATCH(L,N)
              ELSE
                DAMP_FACTOR(L,N) = 0.33333
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      DO L=1,LAND_PTS
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
        RES_FACTOR(I,J) = 0.
      ENDDO

      DO N=1,NTILES
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          STR_RESIST_B = RESB*RESIST_B_TILE(L,N)
          STR_RESIST_S = RESS*RS_TILE(L,N)*DAMP_FACTOR(L,N)
          IF (SNOW_TILE(L,N) >  0.) THEN
            SNOW_F = 1. - EXP(-ASNOW*SNOW_TILE(L,N))
            STR_RESIST_S = 1. /                                         &
     &      (SNOW_F/R_SNOW + (1.-SNOW_F)/STR_RESIST_S)
          ENDIF
          RES_FACTOR(I,J) = RES_FACTOR(I,J) +                           &
     &                           ARESIST(I,J)*TILE_FRAC(L,N) /          &
     &                   (ARESIST_TILE(L,N)+STR_RESIST_B+STR_RESIST_S)
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE RES_TILE
