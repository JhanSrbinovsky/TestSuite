
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!**********************************************************************
! Subroutine to calculate gridbox mean values of surface conductance
! and carbon fluxes. Also returns net primary productivity, leaf
! turnover and wood respiration of each plant functional type for
! driving TRIFFID.
!
!
! Calculates the leaf resistance and net photosynthesis using:
!  (i) Collatz et al. (1992) C3 photosynthesis model, and the
!      Collatz et al. (1991) C4 photosynthesis model.
! (ii) Jacobs (1994) CI/CA closure.
! Written by Peter Cox (February 1996)
! Adapted for MOSES II tile model by Richard Essery (July 1997)
! Coded into UMVN 6.2  by Pete Falloon (July 2006)
!**********************************************************************
      SUBROUTINE LEAF (LAND_FIELD,VEG_PTS,VEG_INDEX,FT                  &
     &,                CLOS_PTS,OPEN_PTS,CLOS_INDEX,OPEN_INDEX          &
     &,                FSMC,TL,CA,CI,RD,WCARB,WEXPT,WLITE               &
     &,                GL,AL)
     
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

      INTEGER                                                           &
     & LAND_FIELD                                                       &
                                  ! IN Total number of land points.
     &,VEG_PTS                                                          &
                                  ! IN Number of vegetated points.
     &,VEG_INDEX(LAND_FIELD)                                            &
                                  ! IN Index of vegetated points
                                  !    on the land grid.
     &,FT                                                               &
                                  ! IN Plant functional type.
     &,CLOS_INDEX(LAND_FIELD)                                           &
                                  ! IN Index of land points
                                  !    with closed stomata.
     &,CLOS_PTS                                                         &
                                  ! IN Number of land points
                                  !    with closed stomata.
     &,OPEN_INDEX(LAND_FIELD)                                           &
                                  ! IN Index of land points
                                  !    with open stomata.
     &,OPEN_PTS                   ! IN Number of land points
                                  !    with open stomata.

      REAL                                                              &
     & FSMC(LAND_FIELD)                                                 &
                                  ! IN Soil water factor.
     &,TL(LAND_FIELD)                                                   &
                                  ! IN Leaf temperature (K).
     &,RD(LAND_FIELD)                                                   &
                                  ! IN Dark respiration (mol CO2/m2/s).
     &,CA(LAND_FIELD)                                                   &
                                  ! IN Canopy CO2 pressure (Pa).
     &,CI(LAND_FIELD)                                                   &
                                  ! IN Internal CO2 pressure (Pa).
     &,WCARB(LAND_FIELD)                                                &      
     &,WLITE(LAND_FIELD)                                                &
     &,WEXPT(LAND_FIELD)                                                &
                                    ! IN Carboxylation,
                                    !    Light, and
                                    !    export limited gross
                                    !    photosynthetic rates
                                    !    (mol CO2/m2/s).

     &,GL(LAND_FIELD)                                                   &
                                  ! OUT Leaf conductance for H2O (m/s).
     &,AL(LAND_FIELD)                                                   &
                                  ! OUT Net Leaf photosynthesis
!                                 !     (mol CO2/m2/s).
     &,B1(LAND_FIELD)                                                   &
     &,B2(LAND_FIELD)                                                   &
     &,B3(LAND_FIELD)                                                   &
                                  ! WORK Coefficients of the quadratic.
     &,CONV(LAND_FIELD)                                                 &
                                  ! WORK Factor for converting mol/m3
!                                 !      into Pa (J/mol).
     &,GLCO2(LAND_FIELD)                                                &
                                  ! WORK Leaf conductance for CO2 (m/s).
     &,WL(LAND_FIELD)                                                   &
                                  ! WORK Gross leaf phtosynthesis
!                                 !      (mol CO2/m2/s).
     &,WP(LAND_FIELD)             ! WORK Smoothed minimum of
!                                 !      Carboxylation and Light
!                                 !      limited gross photosynthesis
!                                 !      (mol CO2/m2/s).


      INTEGER                                                           &
     & J,L                        ! WORK Loop counters.

!-------------------------------------------------------------------
! Parameters
!-------------------------------------------------------------------
      REAL                                                              &
     & BETA1, BETA2                                                     &
                      ! Coupling coefficients for co-limitation.
     &,R                                                                &
                      ! Gas constant (J/K/mol).
     &,RATIO          ! Ratio of leaf resistance for CO2 to leaf
                     ! resistance for H2O.
      PARAMETER (BETA1 = 0.83, BETA2 = 0.93                             &
     &,          R = 8.3144, RATIO = 1.6)

!----------------------------------------------------------------------
! Calculate the co-limited rate of gross photosynthesis
!----------------------------------------------------------------------
!DIR$ IVDEP


      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        B1(L) = BETA1
        B2(L) = - (WCARB(L) + WLITE(L))
        B3(L) = WCARB(L) * WLITE(L)

        WP(L) = -B2(L)/(2*B1(L))                                        &
     &         - SQRT(B2(L)*B2(L)/(4*B1(L)*B1(L)) - B3(L)/B1(L))

      ENDDO

!DIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        B1(L) = BETA2
        B2(L) = - (WP(L) + WEXPT(L))
        B3(L) = WP(L) * WEXPT(L)

        WL(L) = -B2(L)/(2*B1(L))                                        &
     &         - SQRT(B2(L)*B2(L)/(4*B1(L)*B1(L)) - B3(L)/B1(L))

      ENDDO

!----------------------------------------------------------------------
! Carry out calculations for points with open stomata
!----------------------------------------------------------------------
!DIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

!----------------------------------------------------------------------
! Calculate the net rate of photosynthesis
!----------------------------------------------------------------------
        AL(L) = (WL(L) - RD(L)) * FSMC(L)


!----------------------------------------------------------------------
! Calculate the factor for converting mol/m3 into Pa (J/m3).
!----------------------------------------------------------------------
        CONV(L) = R * TL(L)

!----------------------------------------------------------------------
! Diagnose the leaf conductance
!----------------------------------------------------------------------
        GLCO2(L) = (AL(L) * CONV(L)) / (CA(L) - CI(L))
        GL(L) = RATIO * GLCO2(L)

      ENDDO

!----------------------------------------------------------------------
! Close stomata at points with negative or zero net photosynthesis
! or where the leaf resistance exceeds its maximum value.
!----------------------------------------------------------------------
!DIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        IF (GL(L).LE.GLMIN(FT) .OR. AL(L).LE.0.0) THEN
          GL(L) = GLMIN(FT)
          GLCO2(L) = GL(L) / RATIO
          AL(L) = -RD(L)*FSMC(L)
        ENDIF

      ENDDO

!----------------------------------------------------------------------
! Define fluxes and conductances for points with closed stomata
!----------------------------------------------------------------------
!DIR$ IVDEP
      DO J=1,CLOS_PTS
        L = VEG_INDEX(CLOS_INDEX(J))

        GL(L) = GLMIN(FT)
        GLCO2(L) = GL(L) / RATIO
        AL(L) = -RD(L) *FSMC(L)
      ENDDO

      RETURN
      END SUBROUTINE LEAF

