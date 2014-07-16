#if defined(A70_1B) || defined(A70_1C) || defined(A70_1Z)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine to calculate albedos of land-surface tiles and gridbox-mean
! albedo for MOSES II.
!

      SUBROUTINE TILE_ALBEDO (                                          &
     & P_FIELD,LAND_FIELD,LAND_INDEX,NTILES,TILE_PTS,                   &
     & TILE_INDEX,L_SPEC_ALBEDO,ALBSOIL,                                &
     & COSZ,FRAC,LAI_IN,RGRAIN,SNOW_TILE,SOOT,TSTAR_TILE,Z0_TILE,       &
     & ALB_TILE,LAND_ALBEDO,CAN_RAD_MOD                                 &
     & )


! module for land-surface namelist
      USE LAND_SURF_MOD, ONLY :                                         &
     & MASKD                       ! masking depth

      IMPLICIT NONE

#include "nstypes.h"

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER                                                           &
     & P_FIELD                                                          &
                                   ! Total number of grid points.
     &,LAND_FIELD                                                       &
                                   ! No. of land points.
     &,NTILES                      ! Number of surface tiles.

      LOGICAL                                                           &
     & L_SPEC_ALBEDO               ! .TRUE. for spectral albedos
!                                  ! and prognostic snow albedo.

!   Array arguments with intent(in):
      INTEGER                                                           &
     & LAND_INDEX(LAND_FIELD)                                           &
                                   ! Index of land points.
     &,TILE_PTS(NTYPE)                                                  &
                                   ! Number of land points which
!                                  ! include the nth surface type.
     &,TILE_INDEX(LAND_FIELD,NTYPE)! Indices of land points which
!                                  ! include the nth surface type.

      REAL                                                              &
     & ALBSOIL(LAND_FIELD)                                              &
                                   ! Soil albedo.
     &,COSZ(P_FIELD)                                                    &
                                   ! Cosine of the zenith angle.
     &,FRAC(LAND_FIELD,NTYPE)                                           &
                                   ! Fractional cover of each
!                                  ! surface type.
     &,LAI_IN(LAND_FIELD,NPFT)                                          &
                                   ! Leaf area index.
     &,RGRAIN(LAND_FIELD,NTILES)                                        &
                                   ! Snow grain size on tiles
!                                  ! (microns).
     &,SNOW_TILE(LAND_FIELD,NTILES)                                     &
                                   ! Lying snow on tiles (kg/m2).
     &,SOOT(P_FIELD)                                                    &
                                   ! Snow soot content (kg/kg).
     &,TSTAR_TILE(LAND_FIELD,NTILES)                                    &
                                    !Tile surface temperatures (K).
     &,Z0_TILE(LAND_FIELD,NTILES)  ! Surface roughness on tiles (m).

!   Array arguments with intent(out):
      REAL                                                              &
     & ALB_TILE(LAND_FIELD,NTILES,4)                                    &
                                    !Albedos for surface tiles.
!                                  !   (*,*,1) - Direct beam visible
!                                  !   (*,*,2) - Diffuse visible
!                                  !   (*,*,3) - Direct beam near-IR
!                                  !   (*,*,4) - Diffuse near-IR
     &,LAND_ALBEDO(P_FIELD,4)      ! GBM albedos.

! Local arrays:
      REAL                                                              &
     & ALBSNC(LAND_FIELD,NTYPE)                                         &
                                   ! Snow-covered albedo of surf types.
     &,ALBSNF(LAND_FIELD,NTYPE)                                         &
                                   ! Snow-free albedo of surf types.
     &,ALB_TYPE(LAND_FIELD,NTYPE,4)                                     &
                                   ! Albedos of surface types.
     &,ALB_SNOW(LAND_FIELD,NTYPE,4)                                     &
                                   ! Snow albedos.
     &,LAI(LAND_FIELD,NPFT)                                             &
                                   ! Adjusted leaf area index.
     &,SNOW(LAND_FIELD)                                                 &
                                   ! Copy of SNOW_TILE.
     &,TSTAR(LAND_FIELD)                                                &
                                   ! Copy of TSTAR_TILE.
     &,Z0(LAND_FIELD)              ! Copy of Z0_TILE.

      INTEGER, PARAMETER ::       ILAYERS_DUMMY=1                       

      INTEGER                                                           &
     & CAN_RAD_MOD                       ! Which canopy radiation model
                                         ! we're using



      REAL                                                              &
     & FAPAR_DIR_DUMMY(LAND_FIELD,NPFT,ILAYERS_DUMMY)                   &
!                                 ! Profile of absorbed PAR -
!                                 ! Direct beam - DUMMY
     &,FAPAR_DIF_DUMMY(LAND_FIELD,NPFT,ILAYERS_DUMMY)
!                                 ! Profile of absorbed PAR -
!                                 ! Diffuse beam -DUMMY

! Local scalars:
      REAL                                                              &
     & DSA                                                              &
                                   ! Deep-snow albedo.
     &,FLIT                                                             &
                                   ! Weighting factor for albedo.
     &,FSNOW                       ! Weighting factor for albedo.

      INTEGER                                                           &
     & BAND,I,J,L,N                ! Loop counters

#include "c_soilh.h"
#include "c_0_dg_c.h"
      REAL                                                              &
     & DTLAND,KLAND,TCLAND
      PARAMETER( DTLAND = 2., KLAND = 0.3/DTLAND, TCLAND = TM-DTLAND)

#include "pftparm.h"
#include "nvegparm.h"

      DO N=1,NTILES
        DO BAND=1,4
          DO L=1,LAND_FIELD
            ALB_TILE(L,N,BAND) = 0.
          ENDDO
        ENDDO
      ENDDO
      DO N=1,NTYPE
        DO BAND=1,4
          DO L=1,LAND_FIELD
            ALB_TYPE(L,N,BAND) = 0.
            ALB_SNOW(L,N,BAND) = 0.
          ENDDO
        ENDDO
      ENDDO

! Impose minimum LAI for bare vegetation
      DO N=1,NPFT
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          LAI(L,N) = MAX( LAI_IN(L,N), 0.5 )
        ENDDO
      ENDDO

      IF (L_SPEC_ALBEDO) THEN
!----------------------------------------------------------------------
! Spectral albedo scheme with prognostic snow albedo
!----------------------------------------------------------------------

! Set albedos of vegetated surface types
! DEPENDS ON: albpft
        CALL ALBPFT       (P_FIELD,LAND_FIELD,                          &
     &                     LAND_INDEX,TILE_INDEX,TILE_PTS,              &
     &                     ILAYERS_DUMMY,                               &
     &                     ALBSOIL,COSZ,LAI,ALB_TYPE,                   &
     &                     FAPAR_DIR_DUMMY,FAPAR_DIF_DUMMY,             &
     &                     CAN_RAD_MOD)



! Set albedos of non-vegetated surface types
        DO BAND=1,4
          DO N=NPFT+1,NTYPE
            DO J=1,TILE_PTS(N)
              L = TILE_INDEX(J,N)
              ALB_TYPE(L,N,BAND) = ALBSNF_NVG(N-NPFT)
              IF ( ALBSNF_NVG(N-NPFT) <  0. )                           &
                                                     ! Soil tile
     &          ALB_TYPE(L,N,BAND) = ALBSOIL(L)
            ENDDO
          ENDDO
        ENDDO

! Calculate snow albedos
! DEPENDS ON: albsnow
        CALL ALBSNOW(P_FIELD,LAND_FIELD,LAND_INDEX,                     &
     &               NTILES,TILE_INDEX,TILE_PTS,                        &
     &               COSZ,RGRAIN,SNOW_TILE,SOOT,ALB_SNOW)

! Adjust surface type albedos for snow cover
        DO L=1,LAND_FIELD
          SNOW(L) = SNOW_TILE(L,1)
          Z0(L) = Z0_TILE(L,1)
        ENDDO
        DO N=1,NTYPE
          IF (NTILES /= 1) THEN
            DO J=1,TILE_PTS(N)
              L = TILE_INDEX(J,N)
              SNOW(L) = SNOW_TILE(L,N)
              Z0(L) = Z0_TILE(L,N)
            ENDDO
          ENDIF
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            IF ( SNOW(L)  >   0.) THEN
              FSNOW = SNOW(L) / (SNOW(L) + 10.*RHO_SNOW*Z0(L))
              DO BAND=1,4
                ALB_TYPE(L,N,BAND) = FSNOW*ALB_SNOW(L,N,BAND) +         &
     &                               (1. - FSNOW)*ALB_TYPE(L,N,BAND)
              ENDDO
            ENDIF
          ENDDO
        ENDDO

      ELSE
!----------------------------------------------------------------------
! Non-spectral albedo scheme with diagnosed snow albedo
!----------------------------------------------------------------------

! Set albedos of vegetated surface types
        DO N=1,NPFT
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            FLIT = 1.0 - EXP(-KEXT(N)*LAI(L,N))
            ALBSNC(L,N) = ALBSNC_MIN(N)*(1 - FLIT) + ALBSNC_MAX(N)*FLIT
            ALBSNF(L,N) = ALBSOIL(L)*(1 - FLIT) + ALBSNF_MAX(N)*FLIT
          ENDDO
        ENDDO

! Set albedos of non-vegetated surface types
        DO N=NPFT+1,NTYPE
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            ALBSNC(L,N) = ALBSNC_NVG(N-NPFT)
            ALBSNF(L,N) = ALBSNF_NVG(N-NPFT)
            IF ( ALBSNF_NVG(N-NPFT) <  0. ) ALBSNF(L,N) = ALBSOIL(L)
          ENDDO
        ENDDO

! Adjust surface type albedos for snow cover
        DO L=1,LAND_FIELD
          TSTAR(L) = TSTAR_TILE(L,1)
          SNOW(L) = SNOW_TILE(L,1)
        ENDDO
        DO N=1,NTYPE
          IF (NTILES /= 1) THEN
            DO J=1,TILE_PTS(N)
              L = TILE_INDEX(J,N)
              TSTAR(L) = TSTAR_TILE(L,N)
              SNOW(L) = SNOW_TILE(L,N)
            ENDDO
          ENDIF
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            IF ( TSTAR(L)  <   TCLAND ) THEN
              DSA = ALBSNC(L,N)
            ELSEIF ( TSTAR(L)  <   TM ) THEN
              DSA = ALBSNC(L,N) + KLAND*(ALBSNF(L,N) - ALBSNC(L,N))     &
     &                                 *(TSTAR(L) - TCLAND)
            ELSE
              DSA = ALBSNC(L,N) + KLAND*(ALBSNF(L,N) - ALBSNC(L,N))     &
     &                                 *(TM - TCLAND)
            ENDIF
            ALB_TYPE(L,N,1) = ALBSNF(L,N) + (DSA - ALBSNF(L,N)) *       &
     &                                    ( 1. - EXP(-MASKD*SNOW(L)) )
          ENDDO
        ENDDO

! Copy albedo to all bands
        DO BAND=2,4
          DO N=1,NTYPE
            DO J=1,TILE_PTS(N)
              L = TILE_INDEX(J,N)
              ALB_TYPE(L,N,BAND) = ALB_TYPE(L,N,1)
            ENDDO
          ENDDO
        ENDDO

      ENDIF       ! Spectral or non-spectral albedo schemes

!----------------------------------------------------------------------
! Calculate GBM surface albedo
!----------------------------------------------------------------------

      DO BAND=1,4
        DO I=1,P_FIELD
          LAND_ALBEDO(I,BAND) = 0.
        ENDDO
        DO N=1,NTYPE
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            I = LAND_INDEX(L)
            LAND_ALBEDO(I,BAND) = LAND_ALBEDO(I,BAND) +                 &
     &                            FRAC(L,N)*ALB_TYPE(L,N,BAND)
          ENDDO
        ENDDO
      ENDDO

!----------------------------------------------------------------------
! Copy albedos as required for aggregate or distinct tiles
!----------------------------------------------------------------------

      IF (NTILES == 1) THEN
        DO BAND=1,4
          DO L=1,LAND_FIELD
            I = LAND_INDEX(L)
            ALB_TILE(L,1,BAND) = LAND_ALBEDO(I,BAND)
          ENDDO
        ENDDO
      ELSE
        DO BAND=1,4
          DO N=1,NTYPE
            DO J=1,TILE_PTS(N)
              L = TILE_INDEX(J,N)
              ALB_TILE(L,N,BAND) = ALB_TYPE(L,N,BAND)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE TILE_ALBEDO
#endif
#endif
