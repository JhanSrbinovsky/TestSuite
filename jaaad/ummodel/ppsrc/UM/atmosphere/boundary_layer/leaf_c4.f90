
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!**********************************************************************
! Calculates the leaf resistance and net photosynthesis using:
!  (i) Collatz et al. (1992) C3 photosynthesis model
! (ii) Jacobs (1994) CI/CA closure.
!
!    Model            Modification history
!   version  date
!    5.2   15/11/00   New Deck         M. Best
!
!    Programming standard:
!
!**********************************************************************

!**********************************************************************
! Calculates the leaf resistance and net photosynthesis using:
!  (i) Collatz et al. (1991) C4 photosynthesis model
! (ii) Jacobs (1994) CI/CA closure.
!
! Written by Peter Cox (February 1996)
!**********************************************************************
      SUBROUTINE LEAF_C4 (LAND_PTS,VEG_PTS,VEG_INDEX,FT                 &
     &,                   DQ,APAR,TL,CA,OA,PSTAR,FSMC                   &
     &,                   GL,AL,CI,RD)

      IMPLICIT NONE

      INTEGER                                                           &
     & LAND_PTS                                                         &
                                  ! IN Total number of land points.
     &,VEG_PTS                                                          &
                                  ! IN Number of vegetated points.
     &,VEG_INDEX(LAND_PTS)                                              &
                                  ! IN Index of vegetated points
!                                 !    on the land grid.
     &,FT                         ! IN Plant functional type.

      REAL                                                              &
     & DQ(LAND_PTS)                                                     &
                                  ! IN Canopy level specific humidity
!                                 !    deficit (kg H2O/kg air).
     &,APAR(LAND_PTS)                                                   &
                                  ! IN Absorbed PAR (W/m2)
     &,TL(LAND_PTS)                                                     &
                                  ! IN Leaf temperature (K).
     &,CA(LAND_PTS)                                                     &
                                  ! IN Canopy CO2 pressure (Pa).
     &,OA(LAND_PTS)                                                     &
                                  ! IN Atmospheric O2 pressure (Pa).
     &,PSTAR(LAND_PTS)                                                  &
                                  ! IN Atmospheric pressure (Pa).
     &,FSMC(LAND_PTS)                                                   &
                                  ! IN Soil water factor.
     &,GL(LAND_PTS)                                                     &
                                  ! OUT Leaf conductance for H2O (m/s).
     &,AL(LAND_PTS)                                                     &
                                  ! OUT Net Leaf photosynthesis
!                                 !     (mol CO2/m2/s).
     &,RD(LAND_PTS)                                                     &
                                  ! OUT Dark respiration (mol CO2/m2/s).
     &,CI(LAND_PTS)                                                     &
                                  ! OUT Internal CO2 pressure (Pa).
     &,ACR(LAND_PTS)                                                    &
                                  ! WORK Absorbed PAR
!                                 !      (mol photons/m2/s).
     &,B1(LAND_PTS)                                                     &
                                  !
     &,B2(LAND_PTS)                                                     &
                                  !
     &,B3(LAND_PTS)                                                     &
                                  ! WORK Coefficients of the quadratic.
     &,CCP(LAND_PTS)                                                    &
                                  ! WORK Photorespiratory compensatory
!                                 !      point (mol/m3).
     &,CONV(LAND_PTS)                                                   &
                                  ! WORK Factor for converting mol/m3
!                                 !      into Pa (J/mol).
     &,DENOM(LAND_PTS)                                                  &
                                  ! WORK Denominator in equation for VCM
     &,GLCO2(LAND_PTS)                                                  &
                                  ! WORK Leaf conductance for CO2 (m/s).
     &,QTENF(LAND_PTS)                                                  &
                                  ! WORK Q10 function.
     &,TDEGC(LAND_PTS)                                                  &
                                  ! WORK Leaf temperature (deg C).
     &,VCM(LAND_PTS)                                                    &
                                  ! WORK Maximum rate of carboxylation
!                                 !      of Rubisco (mol CO2/m2/s).
     &,VCMAX(LAND_PTS)                                                  &
                                  ! WORK Maximum rate of carboxylation
!                                 !      of Rubisco - without the
!                                 !      temperature factor
!                                 !      (mol CO2/m2/s).
     &,WL(LAND_PTS)                                                     &
                                  ! WORK Gross leaf phtosynthesis
!                                 !      (mol CO2/m2/s).
     &,WCARB(LAND_PTS)                                                  &
                                  ! WORK Carboxylation,
     &,WLITE(LAND_PTS)                                                  &
                                  !      Light, and
     &,WEXPT(LAND_PTS)                                                  &
                                  !      export limited gross
!                                 !      photosynthetic rates
!                                 !      (mol CO2/m2/s).
     &,WP(LAND_PTS)               ! WORK Smoothed minimum of
!                                 !      Carboxylation and Light
!                                 !      limited gross photosynthesis
!                                 !      (mol CO2/m2/s).

      INTEGER                                                           &
     & CLOS_INDEX(LAND_PTS)                                             &
                                  ! WORK Index of land points
!                                 !      with closed stomata.
     &,CLOS_PTS                                                         &
                                  ! WORK Number of land points
!                                 !      with closed stomata.
     &,OPEN_INDEX(LAND_PTS)                                             &
                                  ! WORK Index of land points
!                                 !      with open stomata.
     &,OPEN_PTS                   ! WORK Number of land points
!                                 !      with open stomata.
      INTEGER                                                           &
     & J,L                        ! WORK Loop counters.

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

!-------------------------------------------------------------------
! Parameters
!-------------------------------------------------------------------
      REAL                                                              &
     & BETA1, BETA2                                                     &
                      ! Coupling coefficients for co-limitation.
     &,FDC3, FDC4                                                       &
                      ! Dark respiration coefficients for C3, C4
     &,NEFFC3, NEFFC4                                                   &
                      ! Constant relating VCMAX and leaf N
!                     ! from Schulze et al. 1994 (AMAX = 0.4E-3 * NL
!                     ! - assuming dry matter is 40% carbon by mass)
!                     ! and Jacobs 1994:
!                     ! C3 : VCMAX = 2 * AMAX ; C4 : VCMAX = AMAX
!                     ! (mol/m2/s)
     &,R                                                                &
                      ! Gas constant (J/K/mol).
     &,RATIO                                                            &
                      ! Ratio of leaf resistance for CO2 to leaf
!                     ! resistance for H2O.
     &,ZERODEGC       ! Zero Celsius (K).

      PARAMETER (BETA1 = 0.83, BETA2 = 0.93                             &
     &,          FDC3 = 0.015,  FDC4 = 0.025                            &
     &,          NEFFC3 = 0.8E-3, NEFFC4 = 0.4E-3                       &
     &,          R = 8.3144  , RATIO = 1.6                              &
     &,          ZERODEGC = 273.15)

!----------------------------------------------------------------------
! Initialise counters
!----------------------------------------------------------------------
      CLOS_PTS = 0
      OPEN_PTS = 0

      DO J=1,VEG_PTS
        L = VEG_INDEX(J)

!----------------------------------------------------------------------
! Calculate the points with closed stomata
!----------------------------------------------------------------------
        IF (FSMC(L) == 0.0 .OR. DQ(L) >= DQCRIT(FT)                     &
     &                     .OR. APAR(L) == 0.0) THEN
          CLOS_PTS = CLOS_PTS + 1
          CLOS_INDEX(CLOS_PTS) = J
        ELSE
          OPEN_PTS = OPEN_PTS + 1
          OPEN_INDEX(OPEN_PTS) = J
        ENDIF

!----------------------------------------------------------------------
! Calculate the factor for converting mol/m3 into Pa (J/m3).
!----------------------------------------------------------------------
        CONV(L) = R * TL(L)

      ENDDO

!----------------------------------------------------------------------
! Calculate the photosynthetic parameters
!----------------------------------------------------------------------
!DIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        VCMAX(L) = NEFFC4 * NL0(FT)
        TDEGC(L) = TL(L) - ZERODEGC

        CCP(L) = 0.0

        QTENF(L) = VCMAX(L) * (2.0 ** (0.1 * (TDEGC(L) - 25.0)))
        DENOM(L) = (1 + EXP (0.3 * (TDEGC(L) - TUPP(FT))))              &
     &           * (1 + EXP (0.3 * (TLOW(FT) - TDEGC(L))))
        VCM(L) = QTENF(L) / DENOM(L)

        RD(L) = FDC4 * VCM(L)

      ENDDO

!DIR$ IVDEP
      DO J=1,CLOS_PTS
        L = VEG_INDEX(CLOS_INDEX(J))

        VCMAX(L) = NEFFC4 * NL0(FT)
        TDEGC(L) = TL(L) - ZERODEGC

        QTENF(L) = VCMAX(L) * (2.0 ** (0.1 * (TDEGC(L) - 25.0)))
        DENOM(L) = (1 + EXP (0.3 * (TDEGC(L) - TUPP(FT))))              &
     &           * (1 + EXP (0.3 * (TLOW(FT) - TDEGC(L))))
        VCM(L) = QTENF(L) / DENOM(L)

        RD(L) = FDC4 * VCM(L)

      ENDDO

!----------------------------------------------------------------------
! Carry out calculations for points with open stomata
!----------------------------------------------------------------------
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

!----------------------------------------------------------------------
! Calculate the internal CO2 pressure (Jacobs, 1994).
!----------------------------------------------------------------------
        CI(L) = (CA(L) - CCP(L)) * F0(FT)                               &
     &        * (1 - DQ(L) / DQCRIT(FT)) + CCP(L)

!----------------------------------------------------------------------
! Convert absorbed PAR into mol PAR photons/m2/s
!----------------------------------------------------------------------
        ACR(L) = APAR(L) / 2.19E5

      ENDDO

!----------------------------------------------------------------------
! Calculate the gross photosynthesis for RuBP-Carboxylase, Light and
! Export limited photosynthesis (Collatz et al., 1992).
!----------------------------------------------------------------------
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        WCARB(L) = VCM(L)

        WLITE(L) = ALPHA(FT) * ACR(L)

        WEXPT(L) = 20000.0 * VCM(L) * CI(L) / PSTAR(L)

      ENDDO

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
! Diagnose the leaf conductance
!----------------------------------------------------------------------
        GLCO2(L) = (AL(L) * CONV(L)) / (CA(L) - CI(L))
        GL(L) = GLCO2(L) * RATIO

      ENDDO

!----------------------------------------------------------------------
! Close stomata at points with negative or zero net photosynthesis
! or where the leaf resistance exceeds its maximum value.
!----------------------------------------------------------------------
!DIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        IF (GL(L) <= GLMIN(FT) .OR. AL(L) <= 0.0) THEN
          GL(L) = GLMIN(FT)
          GLCO2(L) = GL(L) / RATIO
          AL(L) = -RD(L) * FSMC(L)
          CI(L) = CA(L) - AL(L) * CONV(L) / GLCO2(L)
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
        AL(L) = -RD(L) * FSMC(L)
        CI(L) = CA(L) - AL(L) * CONV(L) / GLCO2(L)

      ENDDO

      RETURN
      END SUBROUTINE LEAF_C4
