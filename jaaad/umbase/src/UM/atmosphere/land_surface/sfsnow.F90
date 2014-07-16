#if defined(A08_7A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE SFSNOW ------------------------------------------------
!LL
!LL  Purpose:  adds the large-scale and convective snowfall to the
!LL            snowdepth;
!LL
!LL  Model            Modification history from model version 5.2:
!LL version  Date
!LL  5.2   15/11/00   New Deck         M. Best
!LL  5.4   11/04/02   Canopy snow model for needleleaf tree tile.
!LL                   R. Essery
!LL
!LL  Programming standard: Unified Model Documentation Paper No.4
!LL                        version no. 2, dated 18/1/90.
!LL
!LL  Logical component covered: P251.
!LL
!LL  System task:
!LL
!LL  Documentation: um documentation paper no 25
!LLEND------------------------------------------------------------------
!
!*L  ARGUMENTS:---------------------------------------------------------
      SUBROUTINE SFSNOW(                                                &
     & NPNTS,NTILES,TILE_PTS,TILE_INDEX,CAN_MODEL,                      &
     & CONV_SNOW,LS_SNOW,CATCH_SNOW,DZ,HCAP,MELT_TILE,SMCL,STHF,        &
     & TILE_FRAC,TSTAR_TILE,TIMESTEP,V_SAT,                             &
     & RGRAIN,SNOW_TILE,SNOW_GRND,TSOIL,                                &
     & LYING_SNOW,SNOW_MELT,SNOMLT_SUB_HTF,L_SNOW_ALBEDO,               &
     & LTIMER)

      IMPLICIT NONE

      INTEGER                                                           &
     & NPNTS                                                            &
                            ! IN Number of gridpoints.
     &,NTILES                                                           &
                            ! IN Number of tiles.
     &,TILE_PTS(NTILES)                                                 &
                            ! IN Number of tile points.
     &,TILE_INDEX(NPNTS,NTILES)                                         &
!                           ! IN Index of tile points.
     &,CAN_MODEL            ! IN Switch for thermal vegetation canopy.

      REAL                                                              &
     & CONV_SNOW(NPNTS)                                                 &
                            ! IN Convective snowfall (kg/m2/s).
     &,LS_SNOW(NPNTS)                                                   &
                            ! IN Large-scale snowfall (kg/m2/s).
     &,CATCH_SNOW(NPNTS)                                                &
                            ! IN Coniferous canopy snow capacity (kg/m2)
     &,DZ                                                               &
                            ! IN Thicknesses of surface soil layer (m).
     &,HCAP(NPNTS)                                                      &
                            ! IN Soil heat capacity (J/K/m3).
     &,MELT_TILE(NPNTS,NTILES)                                          &
!                           ! IN Surface snowmelt on tiles (kg/m2/s).
     &,SMCL(NPNTS)                                                      &
                            ! IN Soil moisture content of surface soil
!                           !    layer (kg/m2).
     &,STHF(NPNTS)                                                      &
                            ! IN Frozen SMC of surface layer as a
!                           !    fraction of saturation.
     &,TILE_FRAC(NPNTS,NTILES)                                          &
                            ! IN Tile fractions.
     &,TSTAR_TILE(NPNTS,NTILES)                                         &
!                           ! IN Tile surface temperatures (K).
     &,TIMESTEP                                                         &
                            ! IN Timestep (s).
     &,V_SAT(NPNTS)         ! IN SMC at saturation (m3 H2O/m3 soil).

      REAL                                                              &
     & RGRAIN(NPNTS,NTILES)                                             &
                            ! INOUT Snow grain size (microns).
     &,SNOW_TILE(NPNTS,NTILES)                                          &
!                           ! INOUT Snow on the ground (kg/m2).
     &,SNOW_GRND(NPNTS)                                                 &
                            ! INOUT Snow below canopy (kg/m2).
     &,TSOIL(NPNTS)         ! INOUT Surface layer temperature (K).

      REAL                                                              &
     & LYING_SNOW(NPNTS)                                                &
                            ! OUT Gridbox snowmass (kg/m2).
     &,SNOW_MELT(NPNTS)                                                 &
                            ! OUT Snowmelt (kg/m2/s).
     &,SNOMLT_SUB_HTF(NPNTS)! OUT Sub-surface snowmelt heat flux
!                           !     (W/m2).

      LOGICAL                                                           &
     & L_SNOW_ALBEDO                                                    &
                            ! IN Flag for prognostic snow albedo.
     &,LTIMER               ! IN Logical for TIMER.

!---------------------------------------------------------------------
!  External routines called :-

      EXTERNAL TIMER

!-----------------------------------------------------------------------

#include "c_0_dg_c.h"
#include "c_densty.h"
#include "c_lheat.h"
#include "c_perma.h"

! Local variables
      REAL                                                              &
     & ASOIL(NPNTS)                                                     &
                           ! 1 / (dz*hcap) for surface soil layer
     &,INTERCEPT(NPNTS)                                                 &
                           ! Conifer canopy snow interception (kg/m2).
     &,SMCLF(NPNTS)                                                     &
                           ! Frozen soil moisture content of surface
!                          ! layer (kg/m2).
     &,SNOWFALL(NPNTS)                                                  &
                           ! Snowfall in timestep (kg/m2).
     &,SUBMELT(NPNTS)                                                   &
                           ! Melt of snow beneath canopy (kg/m2/s).
     &,UNLOAD(NPNTS)       ! Snow unloading from canopy (kg/m2).
      REAL                                                              &
     & R0                                                               &
                           ! Grain size for fresh snow (microns).
     &,RMAX                                                             &
                           ! Maximum snow grain size (microns).
     &,RATE                ! Grain area growth rate (microns**2 / s).
      PARAMETER (R0 = 50., RMAX = 2000.)
      INTEGER I,J,N        ! Loop counters.


      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SFSNOW  ',103)
      ENDIF

!-----------------------------------------------------------------------
! Increment snowdepth by snowfall.
! Snowfall on NLT tile (N=2) partitioned into interception, throughfall
! and unloading if canopy snow model selected (CAN_MODEL=4).
!-----------------------------------------------------------------------
      DO I=1,NPNTS
        SNOWFALL(I) = TIMESTEP*(LS_SNOW(I) + CONV_SNOW(I))
      ENDDO

      DO N=1,NTILES
        IF (CAN_MODEL == 4 .AND. N == 2) THEN
!CDIR NODEP
          DO J=1,TILE_PTS(N)
            I = TILE_INDEX(J,N)
            INTERCEPT(I) = 0.7*(CATCH_SNOW(I) - SNOW_TILE(I,N))         &
     &                     *(1. - EXP(-SNOWFALL(I)/CATCH_SNOW(I)))
            UNLOAD(I) = 0.4*MELT_TILE(I,N)*TIMESTEP
            UNLOAD(I) = MIN(UNLOAD(I), SNOW_TILE(I,N))
            SNOW_TILE(I,N) = SNOW_TILE(I,N) + INTERCEPT(I) - UNLOAD(I)
            SNOW_GRND(I) = SNOW_GRND(I) + UNLOAD(I) +                   &
     &                     SNOWFALL(I) - INTERCEPT(I)
          ENDDO
        ELSE
!CDIR NODEP
        DO J=1,TILE_PTS(N)
          I = TILE_INDEX(J,N)
          SNOW_TILE(I,N) = SNOW_TILE(I,N) + SNOWFALL(I)
        ENDDO
        ENDIF
      ENDDO

!-----------------------------------------------------------------------
! Calculate melt of snow on ground beneath canopy
!-----------------------------------------------------------------------
      IF (CAN_MODEL == 4 .AND. NTILES >  1) THEN
        N = 2
!CDIR NODEP
        DO J=1,TILE_PTS(N)
          I = TILE_INDEX(J,N)
          IF (TSOIL(I) >  TM) THEN
            SMCLF(I) = RHO_WATER*DZ*V_SAT(I)*STHF(I)
            ASOIL(I) = 1. / ( DZ*HCAP(I) + HCAPW*(SMCL(I) - SMCLF(I))   &
     &                                   + HCAPI*SMCLF(I) )
            SUBMELT(I) = MIN( SNOW_GRND(I)/TIMESTEP ,                   &
     &                        (TSOIL(I) - TM)/(LF*ASOIL(I)*TIMESTEP) )
            SNOW_GRND(I) = SNOW_GRND(I) - SUBMELT(I)*TIMESTEP
            TSOIL(I) = TSOIL(I) -                                       &
     &                 TILE_FRAC(I,N)*ASOIL(I)*TIMESTEP*LF*SUBMELT(I)
            SNOW_MELT(I) = SNOW_MELT(I) + TILE_FRAC(I,N)*SUBMELT(I)
            SNOMLT_SUB_HTF(I) = TILE_FRAC(I,N)*LF*SUBMELT(I)
          ENDIF
        ENDDO
      ENDIF

!----------------------------------------------------------------------
! Calculate the gridbox-mean snow mass and surface melt rate
!----------------------------------------------------------------------
      DO I=1,NPNTS
        LYING_SNOW(I) = 0.
        SNOW_MELT(I) = 0.
      ENDDO
      DO N=1,NTILES
!CDIR NODEP
        DO J=1,TILE_PTS(N)
          I = TILE_INDEX(J,N)
          LYING_SNOW(I) = LYING_SNOW(I) + TILE_FRAC(I,N)*SNOW_TILE(I,N)
          SNOW_MELT(I) = SNOW_MELT(I) + TILE_FRAC(I,N)*MELT_TILE(I,N)
        ENDDO
      ENDDO
      IF (CAN_MODEL == 4 .AND. NTILES >  1) THEN
        N = 2
!CDIR NODEP
        DO J=1,TILE_PTS(N)
          I = TILE_INDEX(J,N)
          LYING_SNOW(I) = LYING_SNOW(I) + TILE_FRAC(I,N)*SNOW_GRND(I)
        ENDDO
      ENDIF

!-----------------------------------------------------------------------
! Increment snow grain size used in albedo calculations
!-----------------------------------------------------------------------
      IF ( L_SNOW_ALBEDO ) THEN
      DO N=1,NTILES
!CDIR NODEP
        DO J=1,TILE_PTS(N)
          I = TILE_INDEX(J,N)
          IF ( SNOW_TILE(I,N)  >   0.) THEN
            RATE = 0.6
            IF (TSTAR_TILE(I,N)  <   TM) THEN
              IF (RGRAIN(I,N)  <   150.) THEN
                RATE = 0.06
              ELSE
                RATE = 0.23E6*EXP(-3.7E4/(8.13451*TSTAR_TILE(I,N)))
              ENDIF
            ENDIF
            RGRAIN(I,N) = SQRT( RGRAIN(I,N)**2                          &
     &                          + (RATE/3.14159)*TIMESTEP )             &
     &                              - (RGRAIN(I,N) - R0)*SNOWFALL(I)/2.5
            RGRAIN(I,N) = MIN( RMAX, RGRAIN(I,N) )
            RGRAIN(I,N) = MAX( R0, RGRAIN(I,N) )
          ELSE
            RGRAIN(I,N) = R0
          ENDIF
        ENDDO
      ENDDO
      ENDIF

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SFSNOW  ',104)
      ENDIF

      RETURN
      END SUBROUTINE SFSNOW
#endif
