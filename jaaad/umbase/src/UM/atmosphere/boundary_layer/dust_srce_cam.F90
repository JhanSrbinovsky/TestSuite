#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
!
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! subroutine DUST_SRCE_CAM
!
! This subroutine contains an old version of dust_srce that is
! currently used in the NWP dust setup in the Crisis Area Models,
! which is based on an old configuration of the HadCM3 dust code. 
! In the future, we hope to unify the code used in  HadGEM and NWP.
!
! Purpose:
!   To calculate flux of mineral dust entrained into atmosphere
!   (currently only calculates flux from bare soil tiles)
!
! Called by bdylyr
!
! Current owners of code: S.Woodward
!
!----------------------------------------------------------------------
!
       SUBROUTINE DUST_SRCE_CAM(                                        &
     & LAND_PTS,NTILES,SM_LEVELS,TILE_PTS,TILE_INDEX,TILE_FRAC,FLAND,   &
     & PSTAR_LAND,TSTAR_TILE,RHOSTAR_LAND,SOIL_LAYER_MOISTURE,SNOW_TILE,&
     & U_S_STD_TILE,MREL_LAND,CLAY_LAND,                                &
     & HO2R2_OROG,                                                      &
     & DUST_FLUX_TILE,U_S_T_TILE,U_S_T_DRY_TILE,                        &
     & SAND_LAND                                                        &
     & )

       IMPLICIT NONE

#include "c_dust_ndiv.h"
#include "c_sulchm.h"
#include "c_dustgen_cam.h"
#include "nstypes.h"
#include "c_g.h"

       INTEGER                                                          &
               !IN
     & LAND_PTS                                                         &
                 ! IN No.of land points in whole grid
     &,NTILES                                                           &
                ! IN No. of land-surface tile types
     &,SM_LEVELS !IN No. of soil levels


       INTEGER                                                          &
               !IN
     & TILE_PTS(NTYPE)                                                  &
                       ! IN Total number of tiles
     &,TILE_INDEX(LAND_PTS,NTYPE) !IN index of tiles on landpts

       REAL                                                             &
            !IN
     & TILE_FRAC(LAND_PTS,NTILES)                                       &
                                  !IN tile fraction
     &,FLAND(LAND_PTS)                                                  &
                       !IN land fraction on land points
     &,PSTAR_LAND(LAND_PTS)                                             &
                            !IN surf. pressure on land pts
     &,TSTAR_TILE(LAND_PTS,NTILES)                                      &
                                   !IN surf. temp on tiles
     &,RHOSTAR_LAND(LAND_PTS)                                           &
                              !IN surf. air density on land pts
     &,SOIL_LAYER_MOISTURE(LAND_PTS,SM_LEVELS)                          &
                                               !IN Soil moisture(kg m-2)
     &,SNOW_TILE(LAND_PTS,NTILES)                                       &
                                  !IN Lying snow on tiles (kg m-2)
     &,U_S_STD_TILE(LAND_PTS,NTILES)                                    &
                                     !IN friction velocity on tiles
     &,MREL_LAND(LAND_PTS,NDIV)                                         &
                                !IN relative soil mass per size division
     &,CLAY_LAND(LAND_PTS)                                              &
                           !IN soil clay fraction
     & ,HO2R2_OROG(LAND_PTS)                                            &
                             !half peak to trough ht of orography
     & ,SAND_LAND(LAND_PTS) !IN soil sand fraction



      REAL                                                              &
           !OUT
     & DUST_FLUX_TILE(LAND_PTS,NTILES,NDIV)                             &
                                            !OUT dust flux (kg m-2 s-1)
     &,U_S_T_TILE(LAND_PTS,NTILES,NDIV)                                 &
                                        !thresh. friction vel. per tile
     &,U_S_T_DRY_TILE(LAND_PTS,NTILES,NDIV)!thresh. frict. vel. per tile
!                                       !excluding soil moisture effects


      INTEGER, PARAMETER :: NDIVH = 9 ! LOCAL: nr divs for horiz flux

      INTEGER                                                           &
              !LOCAL
     & I                                                                &
         !loop counter
     &,J                                                                &
         !loop counter
     &,L                                                                &
         !index of land pt
     &,M                                                                &
         !loop counter, tile types
     &,N                                                                &
         !loop counter, tile points
     &,IDIV !loop counter, dust divisions

      REAL                                                              &
           !LOCAL
     & RATIO                                                            &
             ! U_T_S/U_S_STD
     & ,MREL7(LAND_PTS)                                                 &
                        !local temp: mrel for 31.6 - 100 um radius
     & ,MREL8(LAND_PTS)                                                 &
                        !local temp: mrel for 100 - 316 um
     & ,MREL9(LAND_PTS)                                                 &
                        !local temp: mrel for 316 - 1000 um
     & ,MREL_LAND_ALL(LAND_PTS,NDIVH) !all the mrels



      REAL, PARAMETER :: U_S_MIN =1.E-5 ! LOCAL: minimum val of U_S_STD,
!                                       !below which no dust is produced
!                                       !(avoids divide by zero problem)
      REAL, PARAMETER :: CLAY_MAX =0.2 ! LOC max clay fraction.
      REAL, PARAMETER :: SNOWMIN =1.E-6 ! LOC min snow depth for Dust
      REAL, PARAMETER :: HORIZ_D = 2. !LOC  factor in Hflux eqn (value
!                                           for use with I&W U*td)
      REAL, PARAMETER :: VERT_C = 0.01 !LOC cgs to SI conversion
      REAL, PARAMETER :: H_OROG_LIMIT = 150. !LOC
!                                     1/2pk-trough height above which
!                                     no dust produced
!    U*td from J Iversen and B White, 1982: "Saltation threshold on
!    earth, Mars and Venus". Sedimentology, 29, 111119.
!    These were standard vn6.3 and before.
      REAL, PARAMETER, DIMENSION(NDIVH) :: USTD_BAS = (/                &
     & 1.61,1.33,1.05,0.77,0.50,0.24,0.21,0.32,0.57 /) 
       REAL, PARAMETER, DIMENSION(NDIVH) :: MTOS = (/                   &
     & 2.939E2,9.294E1,2.939E1,9.294E0,2.939E0,9.294E-1,2.939E-1,       &
     & 9.294E-2,2.939E-2   /) !conversion from Mrel to Srel


      REAL                                                              &
           !LOCAL
     & HORIZ_FLUX_TILE(LAND_PTS,NTILES,NDIV)                            &
                                             !horizontal flux per tile
     &, TOT_HORIZ_FLUX_TILE(LAND_PTS,NTILES)                            &
                                             !total over all divs
     &, HORIZ_FLUX_ALL(LAND_PTS,NTILES,NDIV+3)                          &
                                               !per div
     &, HORIZ_FLUX_789(LAND_PTS,NTILES)                                 &
                                        !total over div 7,8,9
     &, MTOT_LAND(LAND_PTS)                                             &
                            !total of mrel_land over all divs
     &, UST_ALL(LAND_PTS,NTILES,NDIVH)  !U*t over all divs                          
                                       




!
!initialisation
!
      DO IDIV = 1,NDIV
        DO M = 1,NTILES
          DO N = 1,TILE_PTS(M)
            L = TILE_INDEX(N,M)
            U_S_T_TILE(L,M,IDIV) = 0.
            U_S_T_DRY_TILE(L,M,IDIV) = 0.
            HORIZ_FLUX_TILE(L,M,IDIV) = 0.
            DUST_FLUX_TILE(L,M,IDIV) = 0.
          ENDDO !TILE_PTS
        ENDDO !NTILES
      ENDDO !NDIV

       DO M = 1,NTILES
          DO N = 1,TILE_PTS(M)
             L = TILE_INDEX(N,M)
             TOT_HORIZ_FLUX_TILE(L,M) = 0.
              HORIZ_FLUX_789(L,M) = 0.
          ENDDO
        ENDDO

        DO IDIV=1,NDIVH
         DO M = 1,NTILES
          DO N = 1,TILE_PTS(M)
             L = TILE_INDEX(N,M)
             HORIZ_FLUX_ALL(L,M,IDIV) = 0.
             UST_ALL(L,M,IDIV) = 0.
           ENDDO
          ENDDO
        ENDDO

        DO N = 1,TILE_PTS(SOIL)
           L = TILE_INDEX(N,SOIL)
           MTOT_LAND(L) = 0.
           DO IDIV = 1,NDIV
             IF (MREL_LAND(L,IDIV)  >   0.)                             &
     &        MTOT_LAND(L) = MTOT_LAND(L) + MREL_LAND(L,IDIV)
           ENDDO
         ENDDO

!
! calculate realtive mass for all divs and put into single array
!
        DO N = 1,TILE_PTS(SOIL)
          L = TILE_INDEX(N,SOIL)
           MREL7(L)=SAND_LAND(L)*0.312
           MREL8(L)=SAND_LAND(L)*0.312
           MREL9(L)=SAND_LAND(L)*0.312
           MREL_LAND_ALL(L,7)=MREL7(L)
           MREL_LAND_ALL(L,8)=MREL8(L)
           MREL_LAND_ALL(L,9)=MREL9(L)
        ENDDO !TILE_PTS
       DO IDIV=1,NDIV
        DO N = 1,TILE_PTS(SOIL)
          L = TILE_INDEX(N,SOIL)
          MREL_LAND_ALL(L,IDIV)=MREL_LAND(L,IDIV)
         ENDDO
        ENDDO

!
!calculate threshold friction velocity for bare soil tiles
!
      DO IDIV = 1,NDIV
        DO N = 1,TILE_PTS(SOIL)
          L = TILE_INDEX(N,SOIL)
          U_S_T_TILE(L,SOIL,IDIV) = UST_A*SOIL_LAYER_MOISTURE(L,1) +    &
     &     UST_B*LOG10(DREP(IDIV)) + UST_C
          IF (U_S_T_TILE(L,SOIL,IDIV)  <   0.) THEN
            U_S_T_TILE(L,SOIL,IDIV) = 0.
          ENDIF
        ENDDO !TILE_PTS
      ENDDO !NDIV
!
!calculate dust flux, if any

      DO IDIV = 1,NDIVH
        DO N = 1,TILE_PTS(SOIL)
          L = TILE_INDEX(N,SOIL)
!
! horizontal flux
!
          UST_ALL(L,SOIL,IDIV)=                                         &
     &          USTD_BAS(IDIV)+UST_A*SOIL_LAYER_MOISTURE(L,1)+UST_C

          IF ((UST_ALL(L,SOIL,IDIV)  <   U_S_STD_TILE(L,SOIL)) .AND.    &
     &    (SNOW_TILE(L,SOIL)  <   SNOWMIN) .AND.                        &
     &    (MREL_LAND_ALL(L,IDIV) > 0.) .AND.                            &
     &    (HO2R2_OROG(L)  <=  H_OROG_LIMIT) .AND.                       &
     &    (U_S_STD_TILE(L,SOIL)  >   U_S_MIN) ) THEN
             RATIO = UST_ALL(L,SOIL,IDIV) / U_S_STD_TILE(L,SOIL)
             HORIZ_FLUX_ALL(L,SOIL,IDIV) = HORIZ_C * RHOSTAR_LAND(L) *  &
     &      HORIZ_D * U_S_STD_TILE(L,SOIL)**3. * (1.0 + RATIO) *        &
     &      (1.0 - RATIO*RATIO) * MREL_LAND_ALL(L,IDIV) / G
             TOT_HORIZ_FLUX_TILE(L,SOIL)=                               &
     &      TOT_HORIZ_FLUX_TILE(L,SOIL)+HORIZ_FLUX_ALL(L,SOIL,IDIV)
            ENDIF
        ENDDO !TILE_PTS
      ENDDO !NDIV

        DO N = 1,TILE_PTS(SOIL)
          L = TILE_INDEX(N,SOIL)
            HORIZ_FLUX_789(L,SOIL)=HORIZ_FLUX_ALL(L,SOIL,7)+            &
     &       HORIZ_FLUX_ALL(L,SOIL,8)+HORIZ_FLUX_ALL(L,SOIL,9)
         ENDDO
!
! vertical flux
!
      DO IDIV = 1,NDIV
        DO N = 1,TILE_PTS(SOIL)
          L = TILE_INDEX(N,SOIL)
          IF ((MREL_LAND(L,IDIV)  >   0.) .AND.                         &
     &  (SNOW_TILE(L,SOIL)  <   SNOWMIN)) THEN

           DUST_FLUX_TILE(L,SOIL,IDIV) =                                &
     &   VERT_C *                                                       &
     &   TOT_HORIZ_FLUX_TILE(L,SOIL) *                                  &
     &   10**(VERT_A * MIN(CLAY_LAND(L),CLAY_MAX) + VERT_B) *           &
     &    MREL_LAND(L,IDIV)/MTOT_LAND(L)  
         ENDIF
        ENDDO !TILE_PTS
      ENDDO !NDIV

!
! correct fluxes for land fraction - for coastal points
!
      DO IDIV = 1,NDIV
       DO M = 1,NTILES
        DO N = 1,TILE_PTS(M)
         L = TILE_INDEX(N,M)
         DUST_FLUX_TILE(L,M,IDIV) = DUST_FLUX_TILE(L,M,IDIV)*FLAND(L)
        ENDDO !TILE_PTS
       ENDDO !NTILES
      ENDDO !NDIV
!
!
!
      RETURN
      END SUBROUTINE DUST_SRCE_CAM
#endif
