
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!**********************************************************************
! Routine to calculate the bulk stomatal resistance and the canopy
! CO2 fluxes
!
!
!    Model            Modification history
!   version  date
!    5.2   15/11/00   New Deck         M. Best
!
!    Programming standard:
!
!**********************************************************************

!***********************************************************************
! Calculates the canopy resistance, net photosynthesis and transpiration
! by scaling-up the leaf level response using the "Big-Leaf" approach
! of Sellers et al. (1994)
!
! Written by Peter Cox (May 1995)
!***********************************************************************
      SUBROUTINE CANOPY (ROW_LENGTH,ROWS,LAND_PTS,LAND_INDEX            &
     &,                  VEG_PTS,VEG_INDEX                              &
     &,                  FT,DQC,IPAR,TSTAR,CO2C,O2,PSTAR                &
     &,                  FPAR,FSMC,LAI                                  &
     &,                  GC,ANETC,CI,RDC)

      IMPLICIT NONE

      INTEGER                                                           &
     & ROW_LENGTH                                                       &
                                  ! IN Number of points on a row
     &,ROWS                                                             &
                                  ! IN Number of rows in a theta field
     &,LAND_PTS                                                         &
                                  ! IN Total number of land points.
     &,LAND_INDEX(LAND_PTS)                                             &
                                  ! IN Index of land points on the
!                                 !    P-grid.
     &,VEG_PTS                                                          &
                                  ! IN Number of vegetated points.
     &,VEG_INDEX(LAND_PTS)        ! IN Index of vegetated points
!                                 !    on the land grid.

      INTEGER                                                           &
     & FT                         ! IN Plant functional type.

      REAL                                                              &
     & CO2C(LAND_PTS)                                                   &
                                  ! IN Canopy level CO2 concentration
!                                 !    (kg CO2/kg air).
     &,DQC(LAND_PTS)                                                    &
                                  ! IN Canopy level specific humidity
!                                 !    deficit (kg H2O/kg air).
     &,O2                                                               &
                                  ! IN Atmospheric O2 concentration
!                                 !    (kg O2/kg air).
     &,PSTAR(LAND_PTS)                                                  &
                                  ! IN Surface pressure (Pa).
     &,IPAR(ROW_LENGTH,ROWS)                                            &
                                  ! IN Incident PAR (W/m2).
     &,TSTAR(LAND_PTS)                                                  &
                                  ! IN Surface temperature (K).
     &,FPAR(LAND_PTS)                                                   &
                                  ! IN PAR absorption factor.
     &,FSMC(LAND_PTS)                                                   &
                                  ! IN Soil water factor.
     &,LAI(LAND_PTS)              ! IN Leaf area index
!                                 !    (m2 leaf/m2 ground).


      REAL                                                              &
     & ANETC(LAND_PTS)                                                  &
                                  ! OUT Net canopy photosynthesis
!                                 !     (mol CO2/m2/s).
     &,CI(LAND_PTS)                                                     &
                                  ! OUT Internal CO2 concentration
!                                 !     (mol CO2/m3).
     &,GC(LAND_PTS)                                                     &
                                  ! OUT Canopy conductance for H2O
!                                 !     (m/s).
     &,RDC(LAND_PTS)                                                    &
                                  ! OUT Canopy dark respiration
!                                 !     (mol CO2/m2/s).
     &,ANETL(LAND_PTS)                                                  &
                                  ! WORK Net leaf photosynthesis
!                                 !      (mol CO2/m2/s/LAI).
     &,APAR(LAND_PTS)                                                   &
                                  ! WORK PAR absorbed by the top leaf
!                                 !      (W/m2).
     &,CA(LAND_PTS)                                                     &
                                  ! WORK Canopy level CO2 pressure
!                                 !      (Pa).
     &,DQM(LAND_PTS)                                                    &
                                  ! WORK Canopy level humidity
!                                 !      deficit (mol H2O/m3).
     &,GL(LAND_PTS)                                                     &
                                  ! WORK Leaf conductance for H2O
!                                 !      (m/s).
     &,OA(LAND_PTS)                                                     &
                                  ! WORK Atmospheric O2 pressure
!                                 !      (Pa).
     &,RD(LAND_PTS)               ! WORK Dark respiration of top leaf
!                                 !      (mol CO2/m2/s).

      INTEGER                                                           &
     & I,J,K,L                      ! WORK Loop counters.

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

!-----------------------------------------------------------------------
! Parameters
!-----------------------------------------------------------------------
      REAL                                                              &
     & R                          ! Gas constant (J/K/mol)
      PARAMETER (R = 8.3144)
!*L------------------COMDECK CCARBON------------------------------------
! Purpose: declares variables and parameters for the carbon cycle
! History:
! version  date         change
! 5.5      26/02/03     add M_CARBON. C Jones.
!----------------------------------------------------------------------
!carbon cycle and vegetation parameters
      REAL                                                              &
     & M_CO2                                                            &
                                  ! molecular weight of CO2
     &,M_AIR                                                            &
                                  ! molecular weight of dry air
     &,M_CARBON                                                         &
                                  ! molecular weight of carbon
     &,EPSILON                                                          &
                                  ! Ratio of molecular weights of water
!                                 !  and dry air.
     &,EPCO2                                                            &
                                  ! Ratio of molecular weights of CO2
!                                 !  and dry air.
     &,EPO2                                                             &
                                  ! Ratio of molecular weights of O2
!                                 !  and dry air.
     &,CO2CONV_A2O                                                      &
                                  ! conversion factor for atmos to
!                                 !  ocean passing of CO2 (mmr to ppmv)
     &,CO2CONV_O2A                ! conversion factor for ocean to
!                                 !  atmos passing of CO2 flux
!                                 !  (mol C/m2/yr to Kg CO2/m2/s)

      PARAMETER (M_AIR=28.966, EPCO2=1.5194, M_CO2=M_AIR*EPCO2,         &
     &           M_CARBON = 12.0, EPSILON = 0.62198, EPO2 = 1.106)

      PARAMETER (CO2CONV_A2O = M_AIR * 1E6 / M_CO2,                     &
     &           CO2CONV_O2A = M_CO2 * 1e-3 / (360.0 * 24.0 * 3600.0))
!*----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Calculate the atmospheric pressures of CO2 and O2
!-----------------------------------------------------------------------
      DO K=1,VEG_PTS
        L = VEG_INDEX(K)
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH

        CA(L) = CO2C(L) / EPCO2 * PSTAR(L)
        OA(L) = O2 / EPO2 * PSTAR(L)
        DQM(L) = DQC(L) / EPSILON * PSTAR(L) / (R * TSTAR(L))

!-----------------------------------------------------------------------
! Calculate the PAR absorbed by the top leaf
!-----------------------------------------------------------------------
        APAR(L) = (1 - OMEGA(FT)) * IPAR(I,J)

      ENDDO

!-----------------------------------------------------------------------
! Call the leaf level model for the top leaf of the C3 and C4 plants
!-----------------------------------------------------------------------

      IF ( C3(FT)  ==  1 ) THEN

! DEPENDS ON: leaf_c3
        CALL LEAF_C3 (LAND_PTS,VEG_PTS,VEG_INDEX,FT                     &
     &,               DQC,APAR,TSTAR,CA,OA,PSTAR,FSMC                   &
     &,               GL,ANETL,CI,RD)

      ELSE

! DEPENDS ON: leaf_c4
        CALL LEAF_C4 (LAND_PTS,VEG_PTS,VEG_INDEX,FT                     &
     &,               DQC,APAR,TSTAR,CA,OA,PSTAR,FSMC                   &
     &,               GL,ANETL,CI,RD)

      ENDIF

!-----------------------------------------------------------------------
! Scale-up to the canopy level
!-----------------------------------------------------------------------
      DO K=1,VEG_PTS
        L = VEG_INDEX(K)

        ANETC(L) = ANETL(L) * FPAR(L)
        GC(L) = FPAR(L) * GL(L)
        RDC(L) = RD(L) * FPAR(L)

      ENDDO

      RETURN

      END SUBROUTINE CANOPY
