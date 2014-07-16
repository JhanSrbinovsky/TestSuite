#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! SUBROUTINE SF_MELT----------------------------------------------------
! Purpose : Calculates surface melting (snow and sea-ice) and increments
!           surface fluxes to satisfy energy balance.
!           Sub-surface snowmelt is calculated and snowdepth incremented
!           by melt and sublimation in P251.
!
!
!    Model            Modification history
!   version  date
!    5.2   15/11/00   New Deck         M. Best
!    5.3  25/04/01  Add coastal tiling. Nic Gedney
!    5.4   11/04/02   Canopy snow model for needleleaf tree tile.
!                     R. Essery
!    5.4   28/08/02 Bug fix for sea ice.  M. Best
!    5.5   31/01/03 Add seaice catagory scheme J. Ridley
!!!  6.2  21/03/05  Remove hardwiring of implicit scheme weight.
!!!                                              M. Diamantakis
!
!    Programming standard:
!
!-----------------------------------------------------------------------
      SUBROUTINE SF_MELT (                                              &
     & ROW_LENGTH,ROWS,LAND_PTS,NTILES,NICE,LAND_INDEX                  &
     &,TILE_INDEX,TILE_PTS,LTIMER,SIMLT,SMLT,FLANDG                     &
     &,ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,DTRDZ_1,ICE_FRACT            &
     &,RESFT,RHOKH_1,RHOKH_1_SICE,TILE_FRAC,TIMESTEP,GAMMA              &
     &,EI_TILE,FQW_1,FQW_ICE,FTL_1,FTL_ICE,FTL_TILE                     &
     &,TSTAR_SEA,TSTAR_SIC,TSTAR_TILE,SNOW_TILE                         &
     &,EI_LAND,EI_SICE,SICE_MELT,ICE_FRACT_NCAT                         &
     &,SICE_MLT_HTF,SNOMLT_SURF_HTF,SNOWMELT,MELT_TILE                  &
     &,TSTAR_SSI,SICE_MELT0)

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
                            ! IN Total number of land points..
     &,LAND_INDEX(LAND_PTS)                                             &
                            ! IN Index of land points.
     &,NTILES                                                           &
                            ! IN Number of tiles per land point.
     &,NICE                                                             &
                            ! IN Number of seaice catagories
     &,TILE_INDEX(LAND_PTS,NTILES)                                      &
!                           ! IN Index of tile points.
     &,TILE_PTS(NTILES)     ! IN Number of tile points.

      LOGICAL                                                           &
     & LTIMER                                                           &
                            ! IN Logical for TIMER.
     &,SIMLT                                                            &
                            ! IN STASH flag for sea-ice melting ht flux.
     &,SMLT                 ! IN STASH flag for snow melting ht flux.

       REAL                                                             &
     & FLANDG(ROW_LENGTH,ROWS)                                          &
!                           ! IN Fraction of gridbox which is land.
     &,ALPHA1(LAND_PTS,NTILES)                                          &
!                           ! IN Gradients of saturated specific
!                           !    humidity with respect to temp.
!                           !    between the bottom model layer
!                           !    and land tile surfaces.
     &,ALPHA1_SICE(ROW_LENGTH,ROWS)                                     &
!                           ! IN ALPHA1 for sea-ice.
     &,ASHTF(ROW_LENGTH,ROWS)                                           &
!                           ! IN Coefficient to calculate surface
!                           !     heat flux into sea-ice (W/m2/K).
     &,ASHTF_TILE(LAND_PTS,NTILES)                                      &
!                           ! IN Coefficient to calculate surface
!                           !    heat flux into soil.
     &,DTRDZ_1(ROW_LENGTH,ROWS)                                         &
!                           ! IN -g.dt/dp for surface layer
     &,ICE_FRACT(ROW_LENGTH,ROWS)                                       &
!                           ! IN Fraction of gridbox which is covered
!                           !    by sea-ice.
     &,ICE_FRACT_NCAT(ROW_LENGTH,ROWS,NICE)                             &
!                           ! IN Fraction of ice in each ice catagory
     &,RESFT(LAND_PTS,NTILES)                                           &
                             !IN Resistance factor.
     &,RHOKH_1(LAND_PTS,NTILES)                                         &
!                           ! IN Surface exchange coeffs for land tiles.
     &,RHOKH_1_SICE(ROW_LENGTH,ROWS)                                    &
                                    ! IN Surface exchange coefficient
!                           !    for sea-ice.
     &,TILE_FRAC(LAND_PTS,NTILES)                                       &
!                           ! IN Tile fractions.
     &,TIMESTEP                                                         &
                            ! IN Timestep (sec).
     &,GAMMA                ! IN implicit weight in level 1

      REAL                                                              &
     & EI_TILE(LAND_PTS,NTILES)                                         &
!                           ! INOUT Sublimation for land tiles (kg/m2/s)
     &,FQW_1(ROW_LENGTH,ROWS)                                           &
!                           ! INOUT GBM surface moisture flux (kg/m2/s).
     &,FQW_ICE(ROW_LENGTH,ROWS)                                         &
!                           ! INOUT FQW for sea-ice.
     &,FTL_1(ROW_LENGTH,ROWS)                                           &
!                           ! INOUT GBM surface sens. heat flux (W/m2).
     &,FTL_ICE(ROW_LENGTH,ROWS)                                         &
!                           ! INOUT FTL for sea-ice.
     &,FTL_TILE(LAND_PTS,NTILES)                                        &
!                           ! INOUT FTL for land tiles.
     &,TSTAR_SEA(ROW_LENGTH,ROWS)                                       &
!                           ! IN Open sea surface temperature (K).
     &,TSTAR_SIC(ROW_LENGTH,ROWS,NICE)                                  &
!                           ! INOUT ice catagory surface temperature (K)
     &,TSTAR_SSI(ROW_LENGTH,ROWS)                                       &
     &,TSTAR_TILE(LAND_PTS,NTILES)                                      &
!                           ! INOUT Land tile surface temperatures (K).
     &,SNOW_TILE(LAND_PTS,NTILES)
!                           ! INOUT Lying snow on land tiles (kg/m2).

      REAL                                                              &
     & EI_LAND(ROW_LENGTH,ROWS)                                         &
!                           ! OUT Sublimation from lying snow (kg/m2/s).
     &,EI_SICE(ROW_LENGTH,ROWS)                                         &
!                           ! OUT Sublimation from sea-ice (kg/m2/s).
     &,MELT_TILE(LAND_PTS,NTILES)                                       &
!                           ! OUT Surface snowmelt on tiles (kg/m2/s).
     &,SICE_MELT(ROW_LENGTH,ROWS,NICE)                                  &
!                           ! OUT Surface melt for sea-ice (W/m2).
     &,SICE_MLT_HTF(ROW_LENGTH,ROWS,NICE)                               &
                                         ! OUT Heat flux due to melt of
!                           !     sea-ice (W/m2).
     &,SNOMLT_SURF_HTF(ROW_LENGTH,ROWS)                                 &
!                           ! OUT Heat flux due to surface melting
!                           !     of snow (W/m2).
     &,SNOWMELT(ROW_LENGTH,ROWS)
!                           ! OUT GBM surface snowmelt (kg/m2/s).
      REAL SICE_MELT0(ROW_LENGTH,ROWS)



!  External routines called :-
      EXTERNAL TIMER


#include "c_0_dg_c.h"
#include "c_lheat.h"
#include "c_r_cp.h"

      REAL                                                              &
     & DFQW                                                             &
                            ! Moisture flux increment.
     &,DFTL                                                             &
                            ! Sensible heat flux increment.
     &,DTSTAR                                                           &
                            ! Surface temperature increment.
     &,LCMELT                                                           &
                            ! Temporary in melt calculations.
     &,LSMELT                                                           &
                            ! Temporary in melt calculations.
     &,RHOKH1_PRIME                                                     &
                            ! Modified forward time-weighted
!                           ! transfer coefficient.
     &,SNOW_MAX                                                         &
                            ! Snow available for melting.
     &,TSTARMAX             ! Maximum gridbox mean surface temperature
!                           ! at sea points with ice.
      INTEGER                                                           &
     & I,J                                                              &
                            ! Loop counter - full horizontal field.
     &,K                                                                &
                            ! Loop counter - tile field
     &,L                                                                &
                            ! Loop counter - land field.
     &,N                    ! Loop counter - tile index.
!
      IF (LTIMER) THEN
! DEPENDS ON: timer
      CALL TIMER('SFMELT  ',3)
      ENDIF

      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        SICE_MELT0(I,J) = 0.0
        DO N=1,NICE
          SICE_MELT(I,J,N) = 0.0
          IF (SIMLT) SICE_MLT_HTF(I,J,N) = 0.0
        ENDDO
        IF (SMLT) SNOMLT_SURF_HTF(I,J) = 0.0
        SNOWMELT(I,J) = 0.0
        EI_LAND(I,J) = 0.0
        EI_SICE(I,J) = 0.0
       ENDDO
      ENDDO

      DO N=1,NTILES
        DO L=1,LAND_PTS
          MELT_TILE(L,N) = 0.
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
!  Melt snow on land tiles if TSTAR_TILE is greater than TM.
!-----------------------------------------------------------------------
      DO N=1,NTILES
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          SNOW_MAX = MAX( 0.0, SNOW_TILE(L,N) - EI_TILE(L,N)*TIMESTEP )
          IF ( SNOW_MAX >  0.0 .AND. TSTAR_TILE(L,N) >  TM ) THEN
            LCMELT = (CP + LC*ALPHA1(L,N)*RESFT(L,N))*RHOKH_1(L,N)      &
     &               + ASHTF_TILE(L,N)
            LSMELT = LCMELT + LF*ALPHA1(L,N)*RHOKH_1(L,N)
            IF (FRAC_SNOW_SUBL_MELT == 0) THEN
              DTSTAR = - MIN( TSTAR_TILE(L,N) - TM ,                    &
     &                        LF*SNOW_MAX / (LCMELT*TIMESTEP) )
            ELSE IF (FRAC_SNOW_SUBL_MELT == 1) THEN
              DTSTAR = - MIN( (TSTAR_TILE(L,N) - TM) *                  &
     &                        (1.0 - EXP(-MASKD*SNOW_MAX)),             &
     &                        LF*SNOW_MAX / (LCMELT*TIMESTEP) )
            ENDIF
            MELT_TILE(L,N) = - LSMELT*DTSTAR / LF
            DFTL = CP*RHOKH_1(L,N)*DTSTAR
            DFQW = ALPHA1(L,N)*RESFT(L,N)*RHOKH_1(L,N)*DTSTAR
            FTL_TILE(L,N) = FTL_TILE(L,N) + DFTL
            EI_TILE(L,N) = EI_TILE(L,N) + DFQW
            TSTAR_TILE(L,N) = TSTAR_TILE(L,N) + DTSTAR
!-----------------------------------------------------------------------
!  Update gridbox-mean quantities
!-----------------------------------------------------------------------
            DFTL = TILE_FRAC(L,N)*DFTL
            DFQW = TILE_FRAC(L,N)*DFQW
            FTL_1(I,J) = FTL_1(I,J) + FLANDG(I,J)*DFTL
            FQW_1(I,J) = FQW_1(I,J) + FLANDG(I,J)*DFQW
          ENDIF
          EI_LAND(I,J) = EI_LAND(I,J) + TILE_FRAC(L,N)*EI_TILE(L,N)
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
!  Increment snow by sublimation and melt
!-----------------------------------------------------------------------
      DO N=1,NTILES
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          SNOW_TILE(L,N) = SNOW_TILE(L,N) -                             &
     &                     (EI_TILE(L,N) + MELT_TILE(L,N))*TIMESTEP
          SNOWMELT(I,J) = SNOWMELT(I,J) + TILE_FRAC(L,N)*MELT_TILE(L,N)
        ENDDO
      ENDDO
      IF (SMLT) THEN
        DO J=1,ROWS
         DO I=1,ROW_LENGTH
          SNOMLT_SURF_HTF(I,J) = LF*SNOWMELT(I,J)
         ENDDO
        ENDDO
      ENDIF

      IF(NICE  ==  1)THEN
      DO J=1,ROWS
       DO I=1,ROW_LENGTH
         DTSTAR=0.0
        IF ( FLANDG(I,J) <  1.0 .AND. ICE_FRACT(I,J) >  0.0 ) THEN
!-----------------------------------------------------------------------
!   Melt sea-ice if TSTAR > TSTARMAX
!-----------------------------------------------------------------------
          EI_SICE(I,J) = FQW_ICE(I,J)
          TSTARMAX = ICE_FRACT(I,J)*TM                                  &
     &         + (1.0 - ICE_FRACT(I,J))*TSTAR_SEA(I,J)
          IF ( TSTAR_SSI(I,J)  >   TSTARMAX ) THEN
            RHOKH1_PRIME = 1. / ( 1. / RHOKH_1_SICE(I,J)                &
     &                       + ICE_FRACT(I,J)*GAMMA*DTRDZ_1(I,J) )
            DTSTAR = TSTARMAX - TSTAR_SSI(I,J)
            LSMELT = (CP + (LC + LF)*ALPHA1_SICE(I,J))*RHOKH1_PRIME     &
     &                                                 + ASHTF(I,J)
            DFTL = CP * RHOKH1_PRIME * DTSTAR
            DFQW = ALPHA1_SICE(I,J) * RHOKH1_PRIME * DTSTAR
            TSTAR_SSI(I,J) = TSTARMAX
             SICE_MELT0(I,J) = - LSMELT * DTSTAR
             IF (SIMLT) SICE_MLT_HTF(I,J,1) = SICE_MELT0(I,J)
            FTL_1(I,J) = FTL_1(I,J) + (1.0-FLANDG(I,J))*DFTL
            FQW_1(I,J) = FQW_1(I,J) + (1.0-FLANDG(I,J))*DFQW
            EI_SICE(I,J) = EI_SICE(I,J) + DFQW
            FTL_ICE(I,J) = FTL_ICE(I,J) + DFTL
            FQW_ICE(I,J) = FQW_ICE(I,J) + DFQW

          ENDIF

        ENDIF
       ENDDO
      ENDDO
      ELSE
       DO J=1,ROWS
        DO I=1,ROW_LENGTH
         IF ( FLANDG(I,J) <  1.0 .AND. ICE_FRACT(I,J) >  0.0 ) THEN
!-----------------------------------------------------------------------
!   Melt sea-ice if TSTAR > TSTARMAX
!-----------------------------------------------------------------------
           EI_SICE(I,J) = FQW_ICE(I,J)
           DTSTAR=0.0
           DO N=1,NICE
             IF (ICE_FRACT_NCAT(I,J,N) >  0.0) THEN
               TSTARMAX = ICE_FRACT_NCAT(I,J,N)*TM
               IF (TSTAR_SIC(I,J,N) >  TSTARMAX) THEN
                 RHOKH1_PRIME = 1. / ( 1. / RHOKH_1_SICE(I,J)           &
     &                  + ICE_FRACT_NCAT(I,J,N)*GAMMA*DTRDZ_1(I,J) )
                 DTSTAR = TSTARMAX - TSTAR_SIC(I,J,N)
                 LSMELT = (CP + (LC + LF)*ALPHA1_SICE(I,J))*RHOKH1_PRIME&
     &                                                 + ASHTF(I,J)
                 DFTL = CP * RHOKH1_PRIME * DTSTAR
                 DFQW = ALPHA1_SICE(I,J) * RHOKH1_PRIME * DTSTAR
                 TSTAR_SIC(I,J,N) = TSTARMAX
                 SICE_MELT(I,J,N) = - LSMELT * DTSTAR
                 IF (SIMLT) SICE_MLT_HTF(I,J,N) = SICE_MELT(I,J,N)
                 FTL_1(I,J)   = FTL_1(I,J) + (1.0-FLANDG(I,J))*DFTL
                 FQW_1(I,J)   = FQW_1(I,J) + (1.0-FLANDG(I,J))*DFQW
                 EI_SICE(I,J) = EI_SICE(I,J) + DFQW
                 FTL_ICE(I,J) = FTL_ICE(I,J) + DFTL
                 FQW_ICE(I,J) = FQW_ICE(I,J) + DFQW

               ENDIF

             ENDIF
           ENDDO
         ENDIF
        ENDDO
       ENDDO
      ENDIF


      IF (LTIMER) THEN
! DEPENDS ON: timer
      CALL TIMER('SFMELT  ',4)
      ENDIF

      RETURN
      END SUBROUTINE SF_MELT
#endif
