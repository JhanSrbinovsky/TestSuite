
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!**********************************************************************
! Calculates the leaf resistance and net photosynthesis using:
!  (i) Collatz et al. (1992) C3 photosynthesis model
! (ii) Jacobs (1994) CI/CA closure.
!**********************************************************************
      SUBROUTINE LEAF_LIMITS(LAND_FIELD,VEG_PTS,VEG_INDEX,FT            &
     &,                      NL,DQ,APAR,TL,CA,OA,PSTAR,FSMC             &
     &,                 CLOS_PTS,OPEN_PTS,CLOS_INDEX,OPEN_INDEX         &
     &,                 CI,RD,WCARB,WEXPT,WLITE ,FAPARV_layer,VCMAX)


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
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------


!  Factors in expressions for limitation of photosynthesis by transport
!  of products, for C3 and C4 respectively. Clim.Dyn. 15, 183-203 (1999)
!  Appendix A.
      real, PARAMETER      ::  fwe_c3 = 0.5
      real, PARAMETER      ::  fwe_c4 = 20000.0
      
!  q10 factor for plant respiration.
      real, PARAMETER      ::  q10_leaf = 2.0
!  (mol/sec) / (watts) conversion for PAR:
      real, PARAMETER      ::  CONPAR = 2.19e5
      
      REAL                                                              &
     & FD(NPFT)                                                         &
                              ! Dark respiration coefficient.
     & ,NEFF(NPFT)            ! Constant relating VCMAX and leaf N

      DATA FD      /   0.015,  0.015,   0.015,     0.025, 0.015 /
      DATA NEFF    /  0.8E-3 , 0.8E-3 , 0.8E-3 ,  0.4E-3, 0.8E-3  /


      INTEGER                                                           &
     & LAND_FIELD                                                       &
                                  ! IN Total number of land points.
     &,VEG_PTS                                                          &
                                  ! IN Number of vegetated points.
     &,VEG_INDEX(LAND_FIELD)                                            &
                                  ! IN Index of vegetated points
                                  !    on the land grid.
     &,FT                         ! IN Plant functional type.

      INTEGER                                                           &
     & CLOS_INDEX(LAND_FIELD)                                           &
                                  ! OUT Index of land points
                                  !     with closed stomata.
     &,CLOS_PTS                                                         &
                                  ! OUT Number of land points
                                  !     with closed stomata.
     &,OPEN_INDEX(LAND_FIELD)                                           &
                                  ! OUT Index of land points
                                  !     with open stomata.
     &,OPEN_PTS                   ! OUT Number of land points
!                                 !     with open stomata.

      REAL                                                              &
     & NL(LAND_FIELD)                                                   &
                                  ! IN Leaf nitrogen
!                                 !    concentration (kg N/kg C).
     &,DQ(LAND_FIELD)                                                   &
                                  ! IN Canopy level specific humidity
!                                 !    deficit (kg H2O/kg air).
     &,APAR(LAND_FIELD)                                                 &
                                  ! IN Absorbed PAR (W/m2)
     &,TL(LAND_FIELD)                                                   &
                                  ! IN Leaf temperature (K).
     &,CA(LAND_FIELD)                                                   &
                                  ! IN Canopy CO2 pressure (Pa).
     &,OA(LAND_FIELD)                                                   &
                                  ! IN Atmospheric O2 pressure (Pa).
     &,PSTAR(LAND_FIELD)                                                &
                                  ! IN Atmospheric pressure (Pa).
     &,FSMC(LAND_FIELD)                                                 &
                                  ! IN Soil water factor.
     &,GL(LAND_FIELD)                                                   &
                                  ! OUT Leaf conductnace for H2O (m/s).
     &,AL(LAND_FIELD)                                                   &
                                  ! OUT Net Leaf photosynthesis
!                                 !     (mol CO2/m2/s).
     &,CI(LAND_FIELD)                                                   &
                                  ! OUT Internal CO2 pressure (Pa).
     &,RD(LAND_FIELD)                                                   &
                                  ! OUT Dark respiration (mol CO2/m2/s).
     &,WCARB(LAND_FIELD)                                                &
                                  ! OUT Carboxylation, ...
     &,WLITE(LAND_FIELD)                                                &
                                  !     ... Light, and ...
     &,WEXPT(LAND_FIELD)                                                &
                                  !     ... export limited gross ...
!                                 !     ... photosynthetic rates ...
!                                 !     ... (mol CO2/m2/s).
     &,ACR(LAND_FIELD)                                                  &
                                  ! WORK Absorbed PAR
!                                 !      (mol photons/m2/s).
     &,CCP(LAND_FIELD)                                                  &
                                  ! WORK Photorespiratory compensatory
!                                 !      point (mol/m3).
     &,CONV(LAND_FIELD)                                                 &
                                  ! WORK Factor for converting mol/m3
!                                 !      into Pa (J/mol).
     &,DENOM(LAND_FIELD)                                                &
                                  ! WORK Denominator in equation for VCM
     &,GLCO2(LAND_FIELD)                                                &
                                  ! WORK Leaf conductnace for CO2 (m/s).
     &,KC(LAND_FIELD)                                                   &
                                  ! WORK Michaelis constant for CO2 (Pa)
     &,KO(LAND_FIELD)                                                   &
                                  ! WORK Michaelis constant for O2 (Pa).
     &,QTENF(LAND_FIELD)                                                &
                                  ! WORK Q10 function.
     &,TAU(LAND_FIELD)                                                  &
                                  ! WORK CO2/O2 specificity ratio.
     &,TDEGC(LAND_FIELD)                                                &
                                  ! WORK Leaf temperature (deg C).
     &,VCM(LAND_FIELD)                                                  &
                                  ! WORK Maximum rate of carboxylation
!                                 !      of Rubisco (mol CO2/m2/s).
     &,VCMAX(LAND_FIELD)                                                &
                                  ! WORK Maximum rate of carboxylation
!                                 !      of Rubisco - without the
!                                 !      temperature factor
!                                 !      (mol CO2/m2/s).

!     LINA
     &, Vm(LAND_FIELD)                                                  &
     &, FAPARV_layer(LAND_FIELD)                                        &
                                 ! IN Profile of absorbed PAR.
     &, RES(LAND_FIELD)                                                 &
                                 ! OUT light inhibited leaf resp.(mol CO2/m2/s)
     &, Leaf_RD(LAND_FIELD)

      INTEGER                                                           &
     & J,L                        ! WORK Loop counters.


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

        IF (FSMC(L)==0.0 .OR. DQ(L) >= DQCRIT(FT)                       &
     &                     .OR. APAR(L)==0.0) THEN
          CLOS_PTS = CLOS_PTS + 1
          CLOS_INDEX(CLOS_PTS) = J
        ELSE
          OPEN_PTS = OPEN_PTS + 1
          OPEN_INDEX(OPEN_PTS) = J
        ENDIF


      ENDDO

!----------------------------------------------------------------------
! Calculate the photosynthetic parameters
!----------------------------------------------------------------------
!DIR$ IVDEP
      DO J=1,OPEN_PTS
        L = VEG_INDEX(OPEN_INDEX(J))

        VCMAX(L) = NEFF(FT) * NL(L)
         TDEGC(L) = TL(L) - ZeroDegC

! TAU is the Rubisco specificity for CO2 relative to O2. The numbers
! in this equation are from Cox, HCTN 24, "Description ... Vegetation
! Model", equation 53.
        TAU(L) = 2600.0 * (0.57 ** (0.1 * (TDEGC(L) - 25.0)))
        CCP(L) = 0.5 * OA(L) / TAU(L) * C3(FT)


        QTENF(L) = VCMAX(L) * (q10_leaf ** (0.1 * (TDEGC(L) - 25.0)))
        DENOM(L) = (1 + EXP (0.3 * (TDEGC(L) - TUPP(FT))))              &
     &           * (1 + EXP (0.3 * (TLOW(FT) - TDEGC(L))))
        VCM(L) = QTENF(L) / DENOM(L)  ! Cox, HCTN 24, equation 49.
        !RD(L) = FD(FT) * VCM(L)
        LEAf_RD(L) = FD(FT) * VCM(L)
        !lina

      ENDDO

!DIR$ IVDEP
      DO J=1,CLOS_PTS
        L = VEG_INDEX(CLOS_INDEX(J))

        VCMAX(L) = NEFF(FT) * NL(L)
         TDEGC(L) = TL(L) - ZeroDegC

! TAU is the Rubisco specificity for CO2 relative to O2. The numbers
! in this equation are from Cox, HCTN 24, "Description ... Vegetation
! Model", equation 53.
         TAU(L) = 2600.0 * (0.57 ** (0.1 * (TDEGC(L) - 25.0)))
         CCP(L) = 0.5 * OA(L) / TAU(L) * REAL( C3(FT) )

         QTENF(L) = VCMAX(L) * (q10_leaf ** (0.1 * (TDEGC(L) - 25.0)))
         DENOM(L) = (1 + EXP (0.3 * (TDEGC(L) - TUPP(FT))))             &
     &           * (1 + EXP (0.3 * (TLOW(FT) - TDEGC(L))))
        VCM(L) = QTENF(L) / DENOM(L)  ! Cox, HCTN 24, equation 49.

        RD(L) = FD(FT) * VCM(L)

!----------------------------------------------------------------------
! Calculate the internal CO2 pressure (Jacobs, 1994).
!----------------------------------------------------------------------
        CI(L) = (CA(L) - CCP(L)) * F0(FT)                               &
     &        * (1 - DQ(L) / DQCRIT(FT)) + CCP(L)

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
        ACR(L) = APAR(L) / CONPAR

! Mercado et al, Tellus 59B, 553-565 (2007), page 560:
        IF (ACR(L)*1.E6 *FAPARV_layer(L) >  10.) then
          RD(L)=( 0.5-0.05 *log(ACR(L)*FAPARV_layer(L)*1.e6))*Leaf_RD(l)
        ELSE
          RD(L)=Leaf_RD(L)
        ENDIF

      ENDDO

!----------------------------------------------------------------------
! Calculate the gross photosynthesis for RuBP-Carboxylase, Light and
! Export limited photosynthesis (Collatz et al., 1992).
!----------------------------------------------------------------------
!DIR$ IVDEP
      IF (C3(FT)==1) THEN

        DO J=1,OPEN_PTS
          L = VEG_INDEX(OPEN_INDEX(J))
! The numbers
! in these 2 equations are from Cox, HCTN 24, "Description ... Vegetation
! Model", equations 54 and 55.
          KC(L) = 30.0 * (2.1 ** (0.1 * (TDEGC(L) - 25.0)))
          KO(L) = 30000.0 * (1.2 ** (0.1 * (TDEGC(L) - 25.0)))

          WCARB(L) = VCM(L) * (CI(L) - CCP(L))                          &
     &             / (CI(L) + KC(L) * (1. + OA(L) / KO(L)))

          WLITE(L) = ALPHA(FT) * ACR(L) * (CI(L) - CCP(L))              &
     &             / (CI(L) + 2 * CCP(L))

          WEXPT(L) = fwe_c3 * VCM(L)
        ENDDO

      ELSE

        DO J=1,OPEN_PTS

          L = VEG_INDEX(OPEN_INDEX(J))

          WCARB(L) = VCM(L)

          WLITE(L) = ALPHA(FT) * ACR(L)

          WEXPT(L) = fwe_c4 * VCM(L) * CI(L) / PSTAR(L)

        ENDDO

      ENDIF

      RETURN
      END SUBROUTINE LEAF_LIMITS
