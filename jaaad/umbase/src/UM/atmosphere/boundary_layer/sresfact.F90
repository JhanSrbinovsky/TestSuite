#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      SUBROUTINE SRESFACT (ROW_LENGTH, ROWS, LAND_PTS, LAND_INDEX,      &
     &                     NTILES, TILE_INDEX, TILE_PTS, SOLUBLE,       &
     &                     CANOPY, CATCH,                               &
     &                     GS_TILE, TILE_FRAC,                          &
     &                     SNOW_TILE,                                   &
     &                     ARESIST, ARESIST_TILE, RESIST_B_TILE,        &
     &                     RESB, RESS, RES_FACTOR_LAND)
!
!----------------------------------------------------------------------
!
!  Purpose: Calculate gridbox-mean resistance factor used by
!           BL_TRMIX_DD to calculate dry deposition of tracers.
!           For 8A (MOSES II tiled land surface) bdy lyr versions,
!           adapted from RESTILE8A.
!
! Current code owner: M Woodage
!
! History:
! Version   Date    Comment
! -------   ----    -------
!
!  5.3   15/10/01   New Deck             M. Woodage
!
!
! Code description:
!   Language: FORTRAN 77  + common extensions
!   This code is written to UMDP3 v6 programming standards.
!
! System Components covered:
!
! System task:
!
! Documentation:  UMDP 20, HCTN 30
!
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
! Arguments with intent IN:
!
      INTEGER                                                           &
     & ROW_LENGTH                                                       &
     &,ROWS                                                             &
     &,LAND_PTS                                                         &
                                    !Total number of land points
     &,LAND_INDEX(LAND_PTS)                                             &
                                    !Index of land points.
! For MOSES II
     &,NTILES                                                           &
                                    !Number of land tiles.
     &,TILE_INDEX(LAND_PTS,NTILES)                                      &
                                    !Index of tile points.
     &,TILE_PTS(NTILES)             !Number of tile points.
!
      LOGICAL                                                           &
     & SOLUBLE                      !.TRUE. for soluble aerosols
!
      REAL                                                              &
     & ARESIST(ROW_LENGTH,ROWS)                                         &
                                    !Aerodynamic resistance Ra (s/m)
     &,RESB                                                             &
                                    !Rb(aerosol) / Rb(H2O).
     &,RESS                                                             &
                                    !Rs(aerosol) / Rs(H2O).
! For MOSES II
     &,ARESIST_TILE(LAND_PTS,NTILES)                                    &
                                    ! 1/(CD_STD*VSHR) on land tiles.
     &,CANOPY(LAND_PTS,NTILES)                                          &
                                    !Surface water on land tiles (kg/m2)
     &,CATCH(LAND_PTS,NTILES)                                           &
                                    !Surface capacity (max. surface
!                                   !    water) of land tiles (kg/m2)
     &,GS_TILE(LAND_PTS,NTILES)                                         &
                                    !Surface conductance for land tiles
     &,RESIST_B_TILE(LAND_PTS,NTILES)                                   &
!                                   !(1/CH-1/CD_STD)/VSHR on land tiles.
     &,SNOW_TILE(LAND_PTS,NTILES)                                       &
                                    !Snow mass on land tiles (kg/m2).
     &,TILE_FRAC(LAND_PTS,NTILES)   !Tile fractions.
!
! Arguments with intent OUT:
!
      REAL                                                              &
     & RES_FACTOR_LAND(ROW_LENGTH,ROWS) !Ra/(Ra+Rb+Rs) mean over land
!                                         part of grid box
!
! Local variables:
!
      INTEGER                                                           &
     & I, J                                                             &
                   !Loop counters (horizontal field index).
     &,K                                                                &
                   !Loop counter (tile field index).
     &,L                                                                &
                   !Loop counter (land point field index).
     &,N           !Loop counter (tile index).
!
      REAL                                                              &
     & DAMP_FACTOR(LAND_PTS,NTILES)                                     &
                                     !Canopy moistening factor
     &,RS_TILE(LAND_PTS,NTILES)                                         &
                                     !Surface resistance for land tiles
!
     &,STR_RESIST_B                                                     &
                                     !Rb for aerosol.
     &,STR_RESIST_S                                                     &
                                     !Rs for aerosol.
     &,SNOW_F                        !Snow cover fraction.
!
! Parameters:
      REAL                                                              &
     & ASNOW                                                            &
                                     !Parameter for snow fraction calcn
     &,COND_LIM                                                         &
                                     !Low limit for canopy conductance.
     &,R_SNOW                        !Resistance to dry dep. over snow

      PARAMETER (ASNOW=0.2,                                             &
     &           COND_LIM=1.0E-3,                                       &
     &           R_SNOW=1.0E3)
!
!
!   Initialise DAMP_FACTOR to 1.0
!   Note that for RESIST_B_TILE, negative values have already been
!   set to in SFEXCH8A, but we repeat it here for safety.
!   For STR_RESIST_S values depend on surface type (land, sea,
!   snow, ice) as well as tracer identity. First calculate stomatal
!   resistance (=1/conductance, avoiding dividing by zero).
!
      DO N=1,NTILES
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          DAMP_FACTOR(L,N) = 1.0
            If (RESIST_B_TILE(L,N)  <   0.0)  Then
              RESIST_B_TILE(L,N) = 0.0
            End If
            IF (GS_TILE(L,N)  >   COND_LIM) THEN
              RS_TILE(L,N) = 1. / GS_TILE(L,N)
            ELSE
              RS_TILE(L,N) = 1. / COND_LIM
            ENDIF
        ENDDO
      ENDDO
!
!  For SOLUBLE species (SO2 and NH3) reduce the surface resistance by up
!  to two-thirds if the canopy is damp (the value of 2/3 is empirical).
!  Two special cases need to be trapped here. The canopy capacity
!  (CATCH) is zero at land ice points, so exclude these from the
!  calculation. Also, there is a possibility that canopy water may
!  exceed canopy capacity due to leaves having fallen, so take care
!  of this too.
!
      IF (SOLUBLE) THEN
! Loop over all land tiles:
        DO N=1,NTILES
          DO K=1,TILE_PTS(N)
            L = TILE_INDEX(K,N)
!
            IF( (CATCH(L,N)  >   0.01) .AND.                            &
     &          (CANOPY(L,N)  >   0.0) ) THEN
              If( CANOPY(L,N)  <=  CATCH(L,N) ) THEN
                DAMP_FACTOR(L,N) = 1. - 0.66667*CANOPY(L,N)/CATCH(L,N)
              Else
                DAMP_FACTOR(L,N) = 0.33333
              End If
            ENDIF
!
          ENDDO
        ENDDO
!
      ENDIF
!
! Initialise RES_FACTOR_LAND to zero:
!
      Do J=1, rows
        Do I=1, row_length
          RES_FACTOR_LAND(I,J) = 0.
        End Do
      End Do
!
!  Need to set STR_RESIST_S to suitable values over snow and ice.
!  Where there is snow cover, calculate an approximate snow fraction
!  for the tile using the formula 1-exp(-ASNOW*SNODEP)
!  Note that for atmospheric model run there should not be any sea
!  points with SNODEP >  0, and land_ice points (including Antarctica)
!  should all have large values of SNOW_TILE.
!
!  Note that the value of ARESIST used here in the calculation of
!  RES_FACTOR_LAND must be the same as that incorporated in RHO_ARESIST
!  (i.e. a grid-box mean including land-sea averaging) used in TR_MIX.
!
      DO N=1,NTILES
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J = (LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
! Loop over all land tiles:
!  (Note: Routine TILEPTS sets array elements TILE_PTS(N) to the no. of
!   gridboxes including surface type N, and TILE_INDEX(K,N) to the land
!   array index of the k-th gridbox containing surface type N.
!   See HCTN 30 p 25.)
!
          STR_RESIST_B = RESB*RESIST_B_TILE(L,N)
          STR_RESIST_S = RESS*RS_TILE(L,N)*DAMP_FACTOR(L,N)
!
            IF (SNOW_TILE(L,N) >  0.0 .AND.                             &
     &                             STR_RESIST_S >  0.0) THEN
              SNOW_F = 1. - EXP(-ASNOW*SNOW_TILE(L,N))
              STR_RESIST_S = 1./                                        &
     &            (SNOW_F/R_SNOW + (1.-SNOW_F)/STR_RESIST_S)
            ENDIF
!
          RES_FACTOR_LAND(I,J) = RES_FACTOR_LAND(I,J) +                 &
     &                           ARESIST(I,J)*TILE_FRAC(L,N) /          &
     &                   (ARESIST_TILE(L,N)+STR_RESIST_B+STR_RESIST_S)
!
        ENDDO
      ENDDO
!
!
      RETURN
      END SUBROUTINE SRESFACT
#endif
