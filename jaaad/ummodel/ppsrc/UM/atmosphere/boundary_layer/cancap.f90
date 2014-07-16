
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!! Subroutine CANCAP ------------------------------------------------
!!!
!!! Purpose : Calculate the heat capacity of a given PFT from its LAI
!!!
!!! version for CRAY YMP
!!!
!!!  Model            Modification history:
!!! version  Date
!!!  5.2  15/11/00 New deck.  M. Best
!!!  5.4  11/04/02   Add can_model=4 option for canopy snow model with
!!!                  needleleaf tree tile.  R. Essery
!!!
!!!END -----------------------------------------------------------------
!**********************************************************************
      SUBROUTINE CANCAP (LAND_PTS,VEG_PTS,VEG_INDEX,CAN_MODEL,FT        &
     &,                  HT,LAI,CANHC,VFRAC)

      IMPLICIT NONE

       INTEGER                                                          &
     & LAND_PTS                                                         &
                                  ! IN Total number of land points.
     &,VEG_PTS                                                          &
                                  ! IN Number of vegetated points.
     &,VEG_INDEX(LAND_PTS)                                              &
                                  ! IN Index of vegetated points.
     &,CAN_MODEL                  ! IN Swith for thermal vegetation
!                                 !    canopy

      INTEGER                                                           &
     & FT                         ! IN Plant functional type.

      REAL                                                              &
     & HT(LAND_PTS)                                                     &
                                  ! IN Vegetation height (m).
     &,LAI(LAND_PTS)                                                    &
                                  ! IN Leaf area index.
     &,CANHC(LAND_PTS)                                                  &
                                  ! OUT Areal heat capacity of
!                                 !     vegetation canopy (J/K/m2).
     &,VFRAC(LAND_PTS)            ! OUT Fractional canopy coverage.

      REAL                                                              &
     & LAI_BAL(LAND_PTS)                                                &
                                  ! WORK Leaf area index in balanced
!                                 !      growth state.
     &,LEAF(LAND_PTS)                                                   &
                                  ! WORK Leaf biomass (kg C/m2).
     &,WOOD(LAND_PTS)             ! WORK Woody biomass (kg C/m2).

      INTEGER                                                           &
     & J,L                        ! WORK Loop counters.

!-----------------------------------------------------------------------
! Parameters
!-----------------------------------------------------------------------
      REAL                                                              &
     & HLEAF                                                            &
                                  ! Specific heat capacity of leaves
!                                 ! (J / K / kg Carbon).
     &,HWOOD                      ! Specific heat capacity of wood
!                                 ! (J / K / kg Carbon).
      PARAMETER ( HLEAF=5.7E4, HWOOD=1.1E4 )

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
!-----------------------------------------------------------------------
!Functional Type dependent parameters
!-----------------------------------------------------------------------
      INTEGER                                                           &
     & C3(NPFT)                                                         &
                                  ! 1 for C3 Plants, 0 for C4 Plants.
     &,CROP(NPFT)                                                       &
                                  ! 1 for crop type, 0 for non-crop.
     &,ORIENT(NPFT)               ! 1 for horizontal, 0 for spherical.

      REAL                                                              &
     & ALPHA(NPFT)                                                      &
                                  ! Quantum efficiency
!                                ! (mol CO2/mol PAR photons).
     &,ALNIR(NPFT)                                                      &
                                  ! Leaf reflection coefficient for
!                                ! near infra-red.
     &,ALPAR(NPFT)                                                      &
                                  ! Leaf reflection coefficient for
!                                ! PAR.
     &,A_WL(NPFT)                                                       &
                                  ! Allometric coefficient relating
!                                ! the target woody biomass to
!                                ! the leaf area index (kg C/m2).
     &,A_WS(NPFT)                                                       &
                                  ! Woody biomass as a multiple of
!                                ! live stem biomass.
     &,B_WL(NPFT)                                                       &
                                  ! Allometric exponent relating
!                                ! the target woody biomass to
!                                ! the leaf area index.
     &,DGL_DM(NPFT)                                                     &
                                  ! Rate of change of leaf turnover
!                                ! rate with moisture availability.
     &,DGL_DT(NPFT)                                                     &
                                  ! Rate of change of leaf turnover
!                                ! rate with temperature (/K)
     &,DQCRIT(NPFT)                                                     &
                                  ! Critical humidity deficit
!                                ! (kg H2O/kg air).
     &,ETA_SL(NPFT)                                                     &
                                  ! Live stemwood coefficient
!                                ! (kg C/m/LAI).
     &,FSMC_OF(NPFT)                                                    &
                                  ! Moisture availability below
!                                ! which leaves are dropped.
     &,F0(NPFT)                                                         &
                                  ! CI/CA for DQ = 0.
     &,GLMIN(NPFT)                                                      &
                                  ! Minimum leaf conductance for H2O
     &,G_AREA(NPFT)                                                     &
                                  ! Disturbance rate (/360days).
     &,G_GROW(NPFT)                                                     &
                                  ! Rate of leaf growth (/360days).
     &,G_LEAF_0(NPFT)                                                   &
                                  ! Minimum turnover rate for leaves
!                                 ! (/360days).
     &,G_ROOT(NPFT)                                                     &
                                  ! Turnover rate for root biomass
!                                 ! (/360days).
     &,G_WOOD(NPFT)                                                     &
                                  ! Turnover rate for woody biomass
!                                 ! (/360days).
     &,KPAR(NPFT)                                                       &
                                  ! PAR Extinction coefficient
!                                ! (m2 leaf/m2 ground).
     &,LAI_MAX(NPFT)                                                    &
                                  ! Maximum projected LAI.
     &,LAI_MIN(NPFT)                                                    &
                                  ! Minimum projected LAI.
     &,NL0(NPFT)                                                        &
                                  ! Top leaf nitrogen concentration
!                                ! (kg N/kg C).
     &,NR_NL(NPFT)                                                      &
                                  ! Ratio of root nitrogen
!                                ! concentration to leaf
!                                ! nitrogen concentration.
     &,NS_NL(NPFT)                                                      &
                                  ! Ratio of stem nitrogen
!                                ! concentration to leaf
!                                ! nitrogen concentration.
     &,OMEGA(NPFT)                                                      &
                                  ! Leaf scattering coefficient
!                                ! for PAR.
     &,OMNIR(NPFT)                                                      &
                                  ! Leaf scattering coefficient for
!                                ! near infra-red.
     &,R_GROW(NPFT)                                                     &
                                  ! Growth respiration fraction.
     &,SIGL(NPFT)                                                       &
                                  ! Specific density of leaf carbon
!                                ! (kg C/m2 leaf).
     &,TLEAF_OF(NPFT)                                                   &
                                  ! Temperature below which leaves are
!                                ! dropped.
     &,TLOW(NPFT)                                                       &
                                  ! Lower temperature for
!                                ! photosynthesis (deg C)
     &,TUPP(NPFT)                 ! Upper temperature for
!                                ! photosynthesis (deg C)

!----------------------------------------------------------------------
!                       BT     NT    C3G    C4G     S
!----------------------------------------------------------------------
      DATA C3      /      1,     1,     1,     0,     1 /
      DATA CROP    /      0,     0,     1,     1,     0 /
      DATA ORIENT  /      0,     0,     0,     0,     0 /
      DATA ALPHA   /   0.08,  0.08,  0.08, 0.040,  0.08 /
      DATA ALNIR   /   0.45,  0.35,  0.58,  0.58,  0.58 /
      DATA ALPAR   /   0.10,  0.07,  0.10,  0.10,  0.10 /
      DATA A_WL    /   0.65,  0.65, 0.005, 0.005,  0.10 /
      DATA A_WS    /  10.00, 10.00,  1.00,  1.00, 10.00 /
      DATA B_WL    /  1.667, 1.667, 1.667, 1.667, 1.667 /
      DATA DGL_DM  /    0.0,   0.0,   0.0,   0.0,   0.0 /
      DATA DGL_DT  /    9.0,   9.0,   0.0,   0.0,   9.0 /
      DATA DQCRIT  /  0.090, 0.060, 0.100, 0.075, 0.100 /
      DATA ETA_SL  /   0.01,  0.01,  0.01,  0.01,  0.01 /
      DATA F0      /  0.875, 0.875, 0.900, 0.800, 0.900 /
      DATA FSMC_OF /   0.00,  0.00,  0.00,  0.00,  0.00 /
      DATA GLMIN   / 1.0E-6,1.0E-6,1.0E-6,1.0E-6,1.0E-6 /
      DATA G_AREA  /  0.005, 0.004,  0.25,  0.25,  0.05 /
      DATA G_GROW  /  20.00, 20.00, 20.00, 20.00, 20.00 /
      DATA G_LEAF_0/   0.25,  0.25,  0.25,  0.25,  0.25 /
      DATA G_ROOT  /   0.25,  0.25,  0.25,  0.25,  0.25 /
      DATA G_WOOD  /   0.01,  0.01,  0.20,  0.20,  0.05 /
      DATA KPAR    /   0.50,  0.50,  0.50,  0.50,  0.50 /
      DATA LAI_MAX /   9.00,  9.00,  4.00,  4.00,  4.00 /
      DATA LAI_MIN /   3.00,  3.00,  1.00,  1.00,  1.00 /
      DATA NL0     /  0.040, 0.030, 0.060, 0.030, 0.030 /
      DATA NR_NL   /   1.00,  1.00,  1.00,  1.00,  1.00 /
      DATA NS_NL   /   0.10,  0.10,  1.00,  1.00,  0.10 /
      DATA OMEGA   /   0.15,  0.15,  0.15,  0.17,  0.15 /
      DATA OMNIR   /   0.70,  0.45,  0.83,  0.83,  0.83 /
      DATA R_GROW  /   0.25,  0.25,  0.25,  0.25,  0.25 /
      DATA SIGL    / 0.0375,0.1000,0.0250,0.0500,0.0500 /
      DATA TLEAF_OF/ 273.15,243.15,258.15,258.15,243.15 /
      DATA TLOW    /    0.0,  -5.0,   0.0,  13.0,   0.0 /
      DATA TUPP    /   36.0,  31.0,  36.0,  45.0,  36.0 /

      DO J=1,VEG_PTS
        L = VEG_INDEX(J)
        CANHC(L) = 0.
        VFRAC(L) = 0.
      ENDDO

      IF (CAN_MODEL  ==  2) THEN
!     Radiative canopy without heat capacity
      DO J=1,VEG_PTS
        L = VEG_INDEX(J)
        CANHC(L) = 0.
        VFRAC(L) = 1. - EXP(-KEXT(FT)*LAI(L))
      ENDDO

      ELSEIF (CAN_MODEL == 3 .OR. CAN_MODEL == 4) THEN
!     Radiative canopy with heat capacity
        DO J=1,VEG_PTS
          L = VEG_INDEX(J)
          LAI_BAL(L) = ( A_WS(FT)*ETA_SL(FT)*HT(L) /                    &
     &                   A_WL(FT) )**(1.0/(B_WL(FT)-1))
          LEAF(L) = SIGL(FT)*LAI_BAL(L)
          WOOD(L) = A_WL(FT)*(LAI_BAL(L)**B_WL(FT))
          CANHC(L) = HLEAF*LEAF(L) + HWOOD*WOOD(L)
          VFRAC(L) = 1. - EXP(-KEXT(FT)*LAI(L))
        ENDDO

      ENDIF

      RETURN
      END SUBROUTINE CANCAP
