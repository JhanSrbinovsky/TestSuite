
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!**********************************************************************
! Routine to calculate the gridbox mean land surface parameters from
! the areal fractions of the surface types and the structural
! properties of the plant functional types.
!
! Written by Peter Cox (June 1997)
!!!  5.2    15/11/00   Re-Written for New Dynamics      M. Best
!!!  5.4    11/04/02   Canopy snow model for needleleaf tree tile.
!!!                    R. Essery
!**********************************************************************
      SUBROUTINE SPARM (LAND_PTS,NTILES,CAN_MODEL                       &
     &,                 TILE_PTS,TILE_INDEX,FRAC,HT,LAI,SATCON          &
     &,                 CATCH_SNOW,CATCH_TILE,INFIL_TILE,Z0_TILE)

      IMPLICIT NONE

!------------------------ nstypes.h ----------------------------------
!jhan:further renovation of ths file may be necessary params are dependent on dataset
!jhan: ALSO nstypes_cable.h should be unecessary nsoil/soil is only difference
      !--- Number of non-vegetation surface types
      Integer, Parameter :: NNVG  = 4

      !--- Number of plant functional types.
      Integer, Parameter :: NPFT  = 13
      
      !--- Number of surface types.
      Integer, Parameter :: NTYPE =17 
      
      !--- Index of the surface type 'Soil'
      !Integer, Parameter :: SOIL  = 16 
      !dhb599, 20110615: change made as per Peter Vohralik, item 1:
      Integer, Parameter :: SOIL  = 14

!--- Land surface types :
!--- original veg. tiles 
!     1 - Broadleaf Tree
!     2 - Needleleaf Tree
!     3 - C3 Grass
!     4 - C4 Grass
!     5 - Shrub
!--- for testing these tiles are set = 1:5 
!     6 - Broadleaf Tree
!     7 - Needleleaf Tree
!     8 - C3 Grass
!     9 - C4 Grass
!    10 - Shrub
!--- for testing these tiles are set = 0
!    11 - 0 
!    11 - 0
!    11 - 0
!--- original non-veg tiles moved to these indices
!     14 - Urban
!     15 - Water
!     16 - Soil
!     17 - Ice







      INTEGER                                                           &
     & LAND_PTS                                                         &
                             ! IN Number of land points to be processed.
     &,NTILES                                                           &
                             ! IN Number of surface tiles.
     &,CAN_MODEL                                                        &
                                    ! IN Swith for thermal vegetation
     &,TILE_PTS(NTYPE)                                                  &
                                    ! IN Number of land points which
!                                   !    include the nth surface type.
     &,TILE_INDEX(LAND_PTS,NTYPE)   ! IN Indices of land points which
!                                   !    include the nth surface type.

      REAL                                                              &
     & FRAC(LAND_PTS,NTYPE)                                             &
                                  ! IN Fractional cover of each
!                                 !    surface type.
     &,HT(LAND_PTS,NPFT)                                                &
                                  ! IN Vegetation height (m).
     &,LAI(LAND_PTS,NPFT)                                               &
                                  ! IN Leaf area index.
     &,SATCON(LAND_PTS)           ! IN Saturated hydraulic
!                                 !    conductivity (kg/m2/s).

      REAL                                                              &
     & CATCH_TILE(LAND_PTS,NTILES)                                      &
                                    ! OUT Canopy capacity for each tile
!                                   !     (kg/m2).
     &,INFIL_TILE(LAND_PTS,NTILES)                                      &
                                    ! OUT Maximum surface infiltration
!                                   !     rate for each tile (kg/m2/s).
     &,Z0_TILE(LAND_PTS,NTILES)                                         &
                                    ! OUT Roughness length for each
!                                   !     tile (m).
     &,CATCH_SNOW(LAND_PTS)         ! OUT Snow capacity for NLT tile
!                                   !     (kg/m2).
      REAL                                                              &
     & CATCH(LAND_PTS)                                                  &
                                  ! WORK GBM canopy capacity (kg/m2).
     &,CATCH_T(LAND_PTS,NTYPE)                                          &
                                  ! WORK Capacities for types.
     &,FZ0(LAND_PTS)                                                    &
                                  ! WORK Aggregation function of Z0.
     &,INFIL(LAND_PTS)                                                  &
                                  ! WORK GBM infiltration rate(kg/m2/s).
     &,INFIL_T(LAND_PTS,NTYPE)                                          &
                                  ! WORK Infiltration rates for types.
     &,Z0(LAND_PTS)                                                     &
                                  ! WORK GBM roughness length (m).
     &,Z0_T(LAND_PTS,NTYPE)       ! WORK Roughness lengths for types.

      INTEGER                                                           &
     & J,L,N                      ! WORK Loop counters

! PFTPARM defines Surface parameters for each Plant Functional Type

      REAL:: ALBSNC_MAX(NPFT)  ! Snow-covered albedo for large LAI.
      REAL:: ALBSNC_MIN(NPFT)  ! Snow-covered albedo for zero LAI.
      REAL:: ALBSNF_MAX(NPFT)  ! Snow-free albedo for large LAI.
      REAL:: DZ0V_DH(NPFT)     ! Rate of change of vegetation
                               ! roughness length with height.
      REAL:: CATCH0(NPFT)      ! Minimum canopy capacity (kg/m2).
      REAL:: DCATCH_DLAI(NPFT) ! Rate of change of canopy capacity
                               ! with LAI.
      REAL:: INFIL_F(NPFT)     ! Infiltration enhancement factor.
      REAL:: KEXT(NPFT)        ! Light extinction coefficient.
      REAL:: ROOTD_FT(NPFT)    ! Rootdepth (m).
      !----------------------------------------------------------------
      !                     BT    NT   C3G   C4G    S
      !----------------------------------------------------------------
      COMMON  /RUN_PFT/ALBSNC_MAX,ALBSNC_MIN,ALBSNF_MAX,DZ0V_DH,        &
     &  CATCH0,DCATCH_DLAI,INFIL_F,KEXT,ROOTD_FT

! PFTPARM end
! NVEGPARM start

! Surface and vegetation parameters
      REAL :: ALBSNC_NVG(NNVG)  ! Snow-covered albedo.
      REAL :: ALBSNF_NVG(NNVG)  ! Snow-free albedo.
      REAL :: CATCH_NVG(NNVG)   ! Canopy capacity (kg/m2).
      REAL :: GS_NVG(NNVG)          ! Surface conductance (m/s).
      REAL :: INFIL_NVG(NNVG)       ! Infiltration enhancement factor.
      REAL :: ROOTD_NVG(NNVG)   ! Rootdepth (m).
      REAL :: Z0_NVG(NNVG)          ! Roughness length (m).
      REAL :: CH_NVG(NNVG)         ! "Canopy" heat capacity (J/K/m2)
      REAL :: VF_NVG(NNVG)         ! Fractional "canopy" coverage

      COMMON  /RUN_BLVEG/ALBSNC_NVG,ALBSNF_NVG,CATCH_NVG,GS_NVG,        &
     &  INFIL_NVG,ROOTD_NVG,Z0_NVG,CH_NVG,VF_NVG

! NVEGPARM end
! Description:
!   This deck sets up the parameter LB
!
! Current Code Owner: Z Gardner
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 5.3      16/08/01 Added in header and changed definitions. Z Gardner
! 5.5      17/04/03 Remove reference to obsolete section
!                   A03_7A. T.White
!
! Declarations:
! Start blend_h
! Description:
!   This file sets the value of the variable LB
!
! Current Code Owner:
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   5.3   25/09/01  Portability changes.  Z. Gardner

      REAL,PARAMETER:: LB = 20.0 ! Blending height (m).
! End Blend_h

!----------------------------------------------------------------------
! Set parameters for vegetated surface types
!----------------------------------------------------------------------
      DO N=1,NPFT
! DEPENDS ON: pft_sparm
        CALL PFT_SPARM (LAND_PTS,N,TILE_INDEX(1,N),TILE_PTS(N)          &
     &,                 HT(1,N),LAI(1,N),SATCON                         &
     &,                 CATCH_T(1,N),INFIL_T(1,N),Z0_T(1,N))
      ENDDO

      IF (CAN_MODEL  ==  4) THEN
        N=2
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          CATCH_SNOW(L) = 4.4*LAI(L,N)
        ENDDO
      ENDIF

!----------------------------------------------------------------------
! Set parameters for non-vegetated surface types
!----------------------------------------------------------------------
      DO N=NPFT+1,NTYPE
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          CATCH_T(L,N) = CATCH_NVG(N-NPFT)
          INFIL_T(L,N) = INFIL_NVG(N-NPFT)*SATCON(L)
          Z0_T(L,N) = Z0_NVG(N-NPFT)
        ENDDO
      ENDDO

      IF (NTILES  ==  1) THEN
!----------------------------------------------------------------------
! Form means and copy to tile arrays if required for aggregate tiles
!----------------------------------------------------------------------
        DO L=1,LAND_PTS
        CATCH(L) = 0.0
        FZ0(L) = 0.0
          INFIL(L) = 0.0
          Z0(L) = 0.0
      ENDDO

      DO N=1,NTYPE
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          FZ0(L) = FZ0(L) + FRAC(L,N) / (LOG(LB / Z0_T(L,N)))**2
        ENDDO
      ENDDO
        DO L=1,LAND_PTS
          Z0(L) = LB * EXP(-SQRT(1. / FZ0(L)))
        ENDDO

        DO N=1,NTYPE
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          CATCH(L) = CATCH(L) + FRAC(L,N) * CATCH_T(L,N)
            INFIL(L) = INFIL(L) + FRAC(L,N) * INFIL_T(L,N)
        ENDDO
      ENDDO

        DO L=1,LAND_PTS
!         Canopy capacity is average over non-lake surface types
          CATCH_TILE(L,1) = 0.
          IF (FRAC(L,7) <  1.)                                          &
     &      CATCH_TILE(L,1) = CATCH(L) / (1. - FRAC(L,7))
          INFIL_TILE(L,1) = INFIL(L)
          Z0_TILE(L,1) = Z0(L)
        ENDDO

      ELSE
!----------------------------------------------------------------------
! Copy surface-type arrays to tiles if separate tiles used
!----------------------------------------------------------------------
        DO N=1,NTYPE
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            CATCH_TILE(L,N) = CATCH_T(L,N)
            INFIL_TILE(L,N) = INFIL_T(L,N)
            Z0_TILE(L,N) = Z0_T(L,N)
      ENDDO
        ENDDO

      ENDIF

      RETURN
      END SUBROUTINE SPARM
