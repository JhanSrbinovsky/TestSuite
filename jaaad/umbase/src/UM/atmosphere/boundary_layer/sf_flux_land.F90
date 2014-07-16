#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!-----------------------------------------------------------------------
!
! Subroutines SF_FLUX_LAND and SF_FLUX_SEA to calculate explicit surface
! fluxes of heat and moisture
!
!
!    Model            Modification history
!   version  date
!    5.2   15/11/00   New Deck         M. Best
!    5.3  25/04/01  Add coastal tiling. Nic Gedney
!  6.1  01/09/04  Calculate potential evaporation related variables.
!                                                          Nic Gedney
!    6.2  07/11/05  Allow for the salinity of sea water in
!                     the evaporative flux.
!                     J. M. Edwards
!
!    Programming standard:
!
!-----------------------------------------------------------------------

!     SUBROUTINE SF_FLUX_LAND-------------------------------------------
!
!     Calculate explicit surface fluxes of heat and moisture over
!     land tiles
!
!     ------------------------------------------------------------------
      SUBROUTINE SF_FLUX_LAND (                                         &
     & ROW_LENGTH,ROWS,LAND_PTS,TILE_PTS,FLAND,LAND_INDEX,TILE_INDEX,   &
     & CANHC,DZSOIL,HCONS,QS1,QSTAR,QW_1,RADNET,RESFT,RHOKH_1,SMVCST,   &
     & SNOW,TILE_FRAC,TIMESTEP,TL_1,TS1,TSTAR,VFRAC,Z0H,Z0M_EFF,Z1_TQ,  &
     & FQW_1_GB,FTL_1_GB,                                               &
     & ALPHA1,ASHTF,FQW_1,EPOT,FTL_1,RHOKPM,RHOKPM_POT,LTIMER,          &
     & ANTHROP_HEAT)

      USE rad_switches_mod, ONLY: LRAD_EMIS_LAND_GEN,RAD_EMIS_LAND_GEN

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
     & LTIMER              ! IN Logical for TIMER

      REAL                                                              &
     & FLAND(LAND_PTS)                                                  &
     &,CANHC(LAND_PTS)                                                  &
                           ! IN Areal heat capacity of canopy (J/K/m2).
     &,DZSOIL                                                           &
                           ! IN Soil or land-ice surface layer
!                          !    thickness (m).
     &,HCONS(LAND_PTS)                                                  &
                           ! IN Soil thermal conductivity (W/m/K).
     &,QS1(ROW_LENGTH,ROWS)                                             &
                           ! IN Sat. specific humidity qsat(TL_1,PSTAR)
     &,QSTAR(LAND_PTS)                                                  &
                           ! IN Surface qsat.
     &,QW_1(ROW_LENGTH,ROWS)                                            &
                            ! IN Total water content of lowest
!                          !    atmospheric layer (kg per kg air).
     &,RADNET(LAND_PTS)                                                 &
                           ! IN Net surface radiation (W/m2) positive
!                          !    downwards
     &,ANTHROP_HEAT                                                     &
!                          ! IN Anthropogenic contribution to surface
!                          !    heat flux (W/m2). Zero except for
!                          !    urban (N=6) and L_ANTHROP_HEAT=.true.
     &,RESFT(LAND_PTS)                                                  &
                           ! IN Total resistance factor.
     &,RHOKH_1(LAND_PTS)                                                &
                           ! IN Surface exchange coefficient.
     &,SMVCST(LAND_PTS)                                                 &
                           ! IN Volumetric saturation point
!                          !    - zero at land-ice points.
     &,SNOW(LAND_PTS)                                                   &
                           ! IN Lying snow amount (kg/m2).
     &,TILE_FRAC(LAND_PTS)                                              &
                           ! IN Tile fraction.
     &,TIMESTEP                                                         &
                           ! IN Timestep (s).
     &,TL_1(ROW_LENGTH,ROWS)                                            &
                            ! IN Liquid/frozen water temperature for
!                          !     lowest atmospheric layer (K).
     &,TS1(LAND_PTS)                                                    &
                           ! IN Temperature of surface layer (K).
     &,TSTAR(LAND_PTS)                                                  &
                           ! IN Surface temperature (K).
     &,VFRAC(LAND_PTS)                                                  &
                           ! IN Fractional canopy coverage.
     &,Z0H(LAND_PTS)                                                    &
                           ! IN Roughness length for heat and moisture
     &,Z0M_EFF(LAND_PTS)                                                &
                           ! IN Effective roughness length for momentum
     &,Z1_TQ(ROW_LENGTH,ROWS)
!                          ! IN Height of lowest atmospheric level (m).

      REAL                                                              &
     & FQW_1_GB(ROW_LENGTH,ROWS)                                        &
                                 ! INOUT GBM surface flux of
!                                !       QW (kg/m2/s).
     &,FTL_1_GB(ROW_LENGTH,ROWS) ! INOUT GBM surface flux of TL.

      REAL                                                              &
     & ASHTF(LAND_PTS)                                                  &
                           ! OUT Coefficient to calculate surface
!                          !     heat flux into soil (W/m2/K).
     &,ALPHA1(LAND_PTS)                                                 &
                           ! OUT Gradient of saturated specific humidity
!                          !     with respect to temperature between the
!                          !     bottom model layer and the surface.
     &,FQW_1(LAND_PTS)                                                  &
                           ! OUT Local surface flux of QW (kg/m2/s).
     &,EPOT(LAND_PTS)                                                   &
                           ! OUT Potential evaporation (kg/m2/s).
     &,FTL_1(LAND_PTS)                                                  &
                           ! OUT Local surface flux of TL.
     &,RHOKPM(LAND_PTS)                                                 &
                           ! OUT Modified surface exchange coefficient.
     &,RHOKPM_POT(LAND_PTS)! OUT Surface exchange coeff. for
!                            potential evaporation.

#include "c_epslon.h"
#include "c_0_dg_c.h"
#include "c_g.h"
#include "c_lheat.h"
#include "c_r_cp.h"
#include "csigma.h"
#include "c_soilh.h"

! Derived local parameters
      REAL GRCP,LS
      PARAMETER (                                                       &
     & GRCP=G/CP                                                        &
     &,LS=LF+LC                                                         &
                           ! Latent heat of sublimation.
     & )

! Scalars
      INTEGER                                                           &
     & I,J                                                              &
                           ! Horizontal field index.
     &,K                                                                &
                           ! Tile field index.
     &,L                   ! Land point field index.

      REAL                                                              &
     & DQ1                                                              &
                           ! (qsat(TL_1,PSTAR)-QW_1) + g/cp*alpha1*Z1
     &,DS_RATIO                                                         &
                           ! 2 * snowdepth / depth of top soil layer.
     &,D_T                                                              &
                           ! Temporary in calculation of alpha1.
     &,LH                                                               &
                           ! Latent heat (J/K/kg).
     &,RAD_REDUC           ! Radiation term required for surface flux

      EXTERNAL TIMER

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SF_FLUX ',3)
      ENDIF

!-----------------------------------------------------------------------
!!  1 Calculate gradient of saturated specific humidity for use in
!!    calculation of surface fluxes
!-----------------------------------------------------------------------
      DO K=1,TILE_PTS
        L = TILE_INDEX(K)
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
        D_T = TSTAR(L) - TL_1(I,J)
        IF (D_T  >   0.05 .OR. D_T  <   -0.05) THEN
          ALPHA1(L) = (QSTAR(L) - QS1(I,J)) / D_T
        ELSEIF (TL_1(I,J)  >   TM) THEN
          ALPHA1(L) = EPSILON*LC*QS1(I,J)*                              &
     &                (1. + C_VIRTUAL*QS1(I,J)) /                       &
     &             ( R*TL_1(I,J)*TL_1(I,J))
        ELSE
          ALPHA1(L) = EPSILON*LS*QS1(I,J)*                              &
     &                (1. + C_VIRTUAL*QS1(I,J)) /                       &
     &             ( R*TL_1(I,J)*TL_1(I,J))
        ENDIF
      ENDDO
 
!-----------------------------------------------------------------------
!! The following calculation on land tiles is repeated for 
!! LRAD_EMIS_LAND_GEN equals false for bit reproducibility.
!! It would be tidier to move the if statement within the do loop
!! but this affects vectorisation.
!-----------------------------------------------------------------------

      If (LRAD_EMIS_LAND_GEN) Then

        DO K=1,TILE_PTS
          L = TILE_INDEX(K)
          ASHTF(L) = 2.0 * HCONS(L) / DZSOIL
          IF (SNOW(L) >  0.0 .AND. SMVCST(L) /= 0.) THEN
            DS_RATIO = 2.0 * SNOW(L) / (RHO_SNOW * DZSOIL)
            IF (DS_RATIO <= 1.0) THEN
              ASHTF(L) =  ASHTF(L) /                                    &
     &                         (1. + DS_RATIO*(HCONS(L)/SNOW_HCON - 1.))
            ELSE
              ASHTF(L) =  ASHTF(L)*SNOW_HCON / HCONS(L)
            ENDIF
          ENDIF
          ASHTF(L) = (1. - VFRAC(L))*ASHTF(L) + CANHC(L)/TIMESTEP +     &
     &     4*(1. + RAD_EMIS_LAND_GEN*VFRAC(L))*RAD_EMIS_LAND_GEN*SBCON  &
     &      *TS1(L)**3
        ENDDO
      
      Else
      
        DO K=1,TILE_PTS
          L = TILE_INDEX(K)
          ASHTF(L) = 2.0 * HCONS(L) / DZSOIL
          IF (SNOW(L) >  0.0 .AND. SMVCST(L) /= 0.) THEN
            DS_RATIO = 2.0 * SNOW(L) / (RHO_SNOW * DZSOIL)
            IF (DS_RATIO <= 1.0) THEN
              ASHTF(L) =  ASHTF(L) /                                    &
     &                         (1. + DS_RATIO*(HCONS(L)/SNOW_HCON - 1.))
            ELSE
              ASHTF(L) =  ASHTF(L)*SNOW_HCON / HCONS(L)
            ENDIF
          ENDIF
          ASHTF(L) = (1. - VFRAC(L))*ASHTF(L) + CANHC(L)/TIMESTEP +     &
     &     4*(1. + VFRAC(L))*SBCON*TS1(L)**3
        ENDDO
      Endif
      
      DO K=1,TILE_PTS
        L = TILE_INDEX(K)
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH

        LH = LC
        IF (SNOW(L)  >   0.) LH = LS
        RHOKPM(L) = RHOKH_1(L) / ( ASHTF(L)  +                          &
     &                     RHOKH_1(L)*(LH*ALPHA1(L)*RESFT(L) + CP) )
        RHOKPM_POT(L)=RHOKH_1(L) / ( ASHTF(L)  +                        &
     &                     RHOKH_1(L)*(LH*ALPHA1(L) + CP) )
        RAD_REDUC = RADNET(L) + ANTHROP_HEAT                            &
                         - ASHTF(L) * ( TL_1(I,J) - TS1(L)              &
     &                   + GRCP*(Z1_TQ(I,J) + Z0M_EFF(L) - Z0H(L)) )    &
     &                   + CANHC(L)*(TSTAR(L) - TS1(L)) / TIMESTEP
        DQ1 = QS1(I,J) - QW_1(I,J) +                                    &
     &                GRCP*ALPHA1(L)*(Z1_TQ(I,J) + Z0M_EFF(L) - Z0H(L))
        FQW_1(L) = RESFT(L)*RHOKPM(L)*( ALPHA1(L)*RAD_REDUC             &
     &                               + (CP*RHOKH_1(L) + ASHTF(L))*DQ1 )
        EPOT(L) = RHOKPM_POT(L)*( ALPHA1(L)*RAD_REDUC                   &
     &                               + (CP*RHOKH_1(L) + ASHTF(L))*DQ1 )
        FTL_1(L) = RHOKPM(L)*(RAD_REDUC - LH*RESFT(L)*RHOKH_1(L)*DQ1)

        FTL_1_GB(I,J) = FTL_1_GB(I,J) + FLAND(L)*TILE_FRAC(L)*FTL_1(L)
        FQW_1_GB(I,J) = FQW_1_GB(I,J) + FLAND(L)*TILE_FRAC(L)*FQW_1(L)
!        if(l.eq.100.or.l.eq.1000.or.l.eq.1500) print *,'sffluxl',l,k,i,j
      ENDDO ! DO K=1,TILE_PTS
!      print *,'sfluxland',FTL_1_GB,FQW_1_GB
!      print *,'sfluxland',EPOT
!       print *,RHOKPM(L),FQW_1(L),EPOT(L),FTL_1(L),FLAND(L), &
!     &         TILE_FRAC(L),FTL_1_GB(I,J),FQW_1_GB(I,J)
!
      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SF_FLUX ',4)
      ENDIF

      RETURN
      END SUBROUTINE SF_FLUX_LAND

!     SUBROUTINE SF_FLUX_SEA--------------------------------------------
!
!     Calculate explicit surface fluxes of heat and moisture over sea
!     and sea-ice
!
!     ------------------------------------------------------------------

#endif
