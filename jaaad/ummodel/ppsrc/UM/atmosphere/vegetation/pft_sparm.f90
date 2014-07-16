
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Routine to calculate the land surface parameters of a given PFT from
! its areal fraction and structural properties.
!
! Written by Peter Cox (June 1997)
!!!  5.2    15/11/00   Re-Written for New Dynamics      M. Best
!**********************************************************************
      SUBROUTINE PFT_SPARM  (LAND_PTS,N,TILE_INDEX,TILE_PTS             &
     &,                      HT,LAI,SATCON,CATCH_T,INFIL_T,Z0_T)


      IMPLICIT NONE

      INTEGER                                                           &
     & LAND_PTS                                                         &
                                  ! IN Number of land points.
     &,N                                                                &
                                  ! IN Plant functional type.
     &,TILE_PTS                                                         &
                                  ! IN Number of land points which
!                                 !    include the surface type.
     &,TILE_INDEX(LAND_PTS)                                             &
                                  ! IN Indices of land points which
!                                 !    include the surface type.
     &,J,L                        ! WORK Loop counters.

      REAL                                                              &
     & HT(LAND_PTS)                                                     &
                                  ! IN Vegetation height (m).
     &,LAI(LAND_PTS)                                                    &
                                  ! IN Leaf area index.
     &,SATCON(LAND_PTS)                                                 &
                                  ! IN Saturated hydraulic conductivity
!                                 !    (kg/m2/s).
     &,CATCH_T(LAND_PTS)                                                &
                                  ! OUT Canopy capacity (kg/m2).
     &,INFIL_T(LAND_PTS)                                                &
                                  ! OUT Maximum surface infiltration
!                                 !     rate (kg/m2/s).
     &,Z0_T(LAND_PTS)             ! OUT Roughness length (m).

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

      DO J=1,TILE_PTS
        L = TILE_INDEX(J)
        Z0_T(L) = DZ0V_DH(N) * HT(L)
        CATCH_T(L) = CATCH0(N) + DCATCH_DLAI(N) * LAI(L)
        INFIL_T(L) = INFIL_F(N) * SATCON(L)
      ENDDO


      RETURN
      END SUBROUTINE PFT_SPARM
