
!
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! subroutine DUST_SRCE
!
!
! Purpose:
!   To calculate flux of mineral dust entrained into atmosphere
!   (currently only calculates flux from bare soil tiles)
!
!   Called by bdylyr
!
! Current owners of code: S.Woodward
!
! History:
! Version     Date     Comment
! -------     ----     -------
!
!   5.5      12/02/03  Original code.  S Wooward
!
!   6.2      15/02/06  Fix to production from coastal points;
!                           add to argument list;
!                           set max clay fraction in dust_flux calc.;
!                           switch off production on steep slopes;
!                           change u*t, horiz_flux and dust_flux
!                           calculations for improved physical realism
!                           and tuning to HadGEM1a-like climate
!                                        Stephanie Woodward
!   6.4      03/01/07  Fix SREL initialisation.
!                      Alternative U*t (Bagnold U*td +A SM + C)
!                      and vertical flux ( A.Hi.(1+H789/H1to6) )
!                                                  Stephanie Woodward

!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!
! Documentation: "Modelling the atmospheric lifecycle..."
!                 Woodward, JGR106, D16, pp18155-18166
!
!----------------------------------------------------------------------
!
       SUBROUTINE DUST_SRCE(                                            &
     & LAND_PTS,NTILES,SM_LEVELS,TILE_PTS,TILE_INDEX,TILE_FRAC,FLAND,   &
     & PSTAR_LAND,TSTAR_TILE,RHOSTAR_LAND,SOIL_LAYER_MOISTURE,SNOW_TILE,&
     & U_S_STD_TILE,MREL_LAND,CLAY_LAND,                                &
     & HO2R2_OROG,                                                      &
     & DUST_FLUX_TILE,U_S_T_TILE,U_S_T_DRY_TILE                         &
     & ,SAND_LAND                                                       &
     & )

       IMPLICIT NONE

!C_DUST_NDIV.............................................................
! Description: Contains parameters for mineral dust code
! Current Code Owner: Stephanie Woodward
!
! History:
! Version  Date     Comment
! -------  ----     -------
!  5.5      12/02/03  Original Code.   Stephanie Woodward
!
! Declarations:
!
      INTEGER NDIV        ! number of particle size divisions
      PARAMETER (NDIV = 6)
!.....................................................................
!-------------------COMDECK C_SULCHM--------------------------------
! Parameters for Sulphur Cycle Chemistry
      REAL                                                              &
     &     EVAPTAU,                                                     &
                          ! timescale for dissolved SO4 to evaporate
     &     NUCTAU,                                                      &
                          ! timescale for accumulation mode particles
!                           to nucleate once they enter a cloud.
     &     DIFFUSE_AIT,                                                 &
                          ! diffusion coefficient of Aitken particles
     &     K_SO2OH_HI,                                                  &
                                  ! high pressure reaction rate limit
     &     K_DMS_OH,                                                    &
                                  ! reaction rate for DMS+OH  cc/mcl/s
     &      K4_CH3SO2_O3,                                               &
                             ! Rate coeff for CH3SO2+O3 -> CH3SO3+O2
     &      K5_CH3SO3_HO2,                                              &
                             ! Rate coeff for CH3SO3+HO2 -> MSA+O2
     &      RMM_O3,                                                     &
                             ! relative molecular mass O3
     &     BRAT_SO2,                                                    &
                                  ! branching ratio for SO2 in DMS oxidn
     &     BRAT_MSA,                                                    &
                                  ! branching ratio for MSA in DMS oxidn
     &     AVOGADRO,                                                    &
                                 ! no. of molecules in 1 mole
     &     RMM_H2O2,                                                    &
                                 ! relative molecular mass H2O2 kg/mole
     &     RMM_AIR,                                                     &
                                 ! relative molecular mass dry air
     &     RMM_W,                                                       &
                                 ! relative molecular mass water
     &     RELM_S_H2O2,                                                 &
                                 ! rel atomic mass sulphur/RMM_H2O2
     &     RELM_S_2N,                                                   &
                              ! rel atomic mass Sulphur/2*Nitrogen
     &     PARH,                                                        &
                                ! power of temp dependence of K_SO2OH_LO
     &     K1,                                                          &
                                ! parameters for calcn of K_SO2OH_LO
     &     T1,                                                          &
                                !
     &     FC,                                                          &
                                ! parameters for interpolation between
     &     FAC1,                                                        &
                                !   LO and HI reaction rate limits
     &     K2,K3,K4,                                                    &
                                ! parameters for calcn of K_HO2_HO2
     &     T2,T3,T4,                                                    &
                                !
     &     CLOUDTAU,                                                    &
                                  ! air parcel lifetime in cloud
     &     CHEMTAU,                                                     &
                                  ! chem lifetime in cloud before oxidn
     &     O3_MIN,                                                      &
                              ! min mmr of O3 required for oxidn
     &     THOLD                  ! threshold for cloud liquid water
!
!
      PARAMETER (                                                       &
     &           EVAPTAU = 300.0,                                       &
                                              ! secs  (=5 mins)
     &             NUCTAU = 30.0,                                       &
                                          ! secs
     &       DIFFUSE_AIT = 1.7134E-9,                                   &
                                             ! sq m/s
     &        K_SO2OH_HI = 2.0E-12,                                     &
                                       ! cc/mcl/s from STOCHEM model
     &           K_DMS_OH = 9.1E-12,                                    &
                                          ! cc/mcl/s
     &       K4_CH3SO2_O3 = 1.0E-14,                                    &
                                        ! cc/mcl/s
     &      K5_CH3SO3_HO2 = 4.0E-11,                                    &
     &             RMM_O3 = 4.8E-2,                                     &
                                        ! kg/mole
     &          BRAT_SO2 = 0.9,                                         &
     &           BRAT_MSA = 1.0-BRAT_SO2,                               &
     &           AVOGADRO = 6.022E23,                                   &
                                          ! per mole
     &           RMM_H2O2 = 3.40E-2,                                    &
                                          ! kg/mole
     &            RMM_AIR = 2.896E-2,                                   &
                                          ! kg/mole
     &              RMM_W = 1.8E-2,                                     &
                                          ! kg/mole
     &        RELM_S_H2O2 = 3.206/3.40,                                 &
     &           RELM_S_2N = 3.206/2.80,                                &
     &               PARH = 3.3,                                        &
     &                K1 = 4.0E-31,                                     &
                                       ! (cc/mcl)2/s from STOCHEM
     &                 T1 = 300.0,                                      &
                                          ! K
     &                FC = 0.45,                                        &
                                        ! from STOCHEM model
     &              FAC1 = 1.1904,                                      &
                                    ! 0.75-1.27*LOG10(FC) from STOCHEM
     &                 K2 = 2.2E-13,                                    &
                                          ! cc/mcl/s
     &                 K3 = 1.9E-33,                                    &
                                          ! (cc/mcl)2/s
     &                 K4 = 1.4E-21,                                    &
                                          ! cc/mcl
     &                 T2 = 600.0,                                      &
                                          ! K
     &                 T3 = 890.0,                                      &
                                          ! K
     &                 T4 = 2200.0,                                     &
                                          ! K
     &           CLOUDTAU = 1.08E4,                                     &
                                          ! secs (=3 hours)
     &            CHEMTAU = 9.0E2,                                      &
                                          ! secs (=15 mins)
     &              O3_MIN = 1.6E-8,                                    &
                                        !(kg/kg, equiv. 10ppbv)
     &              THOLD = 1.0E-8                                      &
                                          ! kg/kg
     &          )
!
      REAL RAD_AIT,                                                     &
                            ! median radius of Aitken mode particles
     &     DIAM_AIT,                                                    &
                            !   "    diameter    "
     &     RAD_ACC,                                                     &
                            ! median radius of acccumulation mode
     &     DIAM_ACC,                                                    &
                            !   "    diameter    "
     &     CHI,                                                         &
                            ! mole fraction of S in particle
     &     RHO_SO4,                                                     &
                            ! density of  SO4 particle
     &     SIGMA,                                                       &
                            ! standard devn of particle size distn
!                                 for accumulation mode
     &     E_PARM,                                                      &
                            ! param relating size distns of Ait & Acc
     &     NUM_STAR         ! threshold concn of accu mode particles
                            !  below which PSI=1
!
      PARAMETER (                                                       &
     &           RAD_AIT = 6.5E-9,                                      &
                                             ! m
     &          DIAM_AIT = 2.0*RAD_AIT,                                 &
     &           RAD_ACC = 95.0E-9,                                     &
                                             ! m
     &          DIAM_ACC = 2.0*RAD_ACC,                                 &
     &               CHI = 32.0/132.0,                                  &
     &           RHO_SO4 = 1769.0,                                      &
                                              ! kg/m3
     &             SIGMA = 1.4,                                         &
     &            E_PARM = 0.9398,                                      &
     &          NUM_STAR = 1.0E6                                        &
                                             ! m-3
     &          )
!
      REAL BOLTZMANN       !Boltzmanns constant.
      REAL MFP_REF         !Reference value of mean free path.
      REAL TREF_MFP        !Reference temperature for mean free path.
      REAL PREF_MFP        !Reference pressure for mean free path.
      REAL SIGMA_AIT       !Geometric standard deviation of the Aitken
!                             mode distribution.
!
      PARAMETER (BOLTZMANN = 1.3804E-23)  ! J K-1
      PARAMETER (MFP_REF = 6.6E-8                                       &
                                          ! m
     &        ,  TREF_MFP = 293.15                                      &
                                          ! K
     &        ,  PREF_MFP = 1.01325E5)    ! Pa
      PARAMETER (SIGMA_AIT = 1.30)
!
!*---------------------------------------------------------------------
!C_DUSTGEN..............................................................
! Description: Contains parameters for mineral dust generation
! Current Code Owner: Stephanie Woodward
!
! History:
! Version  Date     Comment
! -------  ----     -------
!  5.5      12/02/03  Original Code.   Stephanie Woodward
!  6.2      14/02/06  Alternative emissions terms for HadGEM1A.
!                                                     Stephanie Woodward
!  6.4      03/01/07  Alternative U*t terms.   Stephanie Woodward
! Parameters for mineral dust generation
!
      REAL, PARAMETER :: HORIZ_C = 2.61 ! C in horizontal flux calc.     
      REAL, PARAMETER :: VERT_A = 13.4 ! A in vertical flux calc
      REAL, PARAMETER :: VERT_B = -6. ! B in vertical flux calc

      REAL,PARAMETER :: RHOP = 2.65E+3  ! density of a dust particle
      REAL,PARAMETER :: Z0B = 0.0003    !roughness length for bare soil
      REAL Z0S(NDIV)  ! smooth roughness len (calc.d from part size)
      REAL DREP(NDIV) ! representative particle diameter
      REAL DMAX(NDIV) ! max diameter of particles in each div.
      REAL DMIN(NDIV) ! min diameter of particles in each div.
                         ! note that by using two arrays here we can set
                         ! up overlapping divisions, however this means
                         ! we have to be careful to make them consistent
!
      DATA Z0S/ .374894E-08, .118552E-07, .374894E-07, .118552E-06,     &
     &           .374894E-06, .118552E-05/
      DATA DREP/ .112468E-06, .355656E-06, .112468E-05, .355656E-05,    &
     &           .112468E-04, .355656E-04/
      DATA DMAX/2.0E-7,6.32456E-7,2.0E-6,                               &
     &          6.32456E-6,2.0E-5,6.32456E-5/
      DATA DMIN/6.32456E-8,2.0E-7,6.32456E-7,                           &
     &          2.0E-6,6.32456E-6,2.0E-5/
!.......................................................................
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


!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------

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


!     ----------------------------------------------------------------------
!     Original hg3c (access 1.3) code (as at 23apr12a, pfv)
!     ----------------------------------------------------------------------
!     REAL, PARAMETER :: U_S_MIN =1.E-5 ! LOCAL: minimum val of U_S_STD,
!                                       !below which no dust is produced
!                                       !(avoids divide by zero problem)
!     REAL, PARAMETER :: CLAY_MAX =0.2 ! LOC max clay fraction.
!     REAL, PARAMETER :: SNOWMIN =1.E-6 ! LOC min snow depth for Dust
!     REAL, PARAMETER :: HORIZ_D = 1. !LOC  factor in Hflux eqn
!     REAL, PARAMETER :: VERT_C = 0.01 !LOC cgs to SI conversion
!     REAL, PARAMETER :: US_AA = 0.0 !LOC Ustar correction (additional)
!     REAL, PARAMETER :: US_AM = 1.5 !LOC Ustar correction (multiplic.)
!     REAL, PARAMETER :: UST_AA = 0.0 !LOC Ustar_t correction (add)
!     REAL, PARAMETER :: UST_AM = 1.0 !LOC Ustar_t correction (multi.)
!     REAL, PARAMETER :: SM_CORR = 1.7 !soil moist. correction factor 
!     REAL, PARAMETER :: H_OROG_LIMIT = 150. !LOC
!                                     1/2pk-trough height above which
!                                     no dust produced
!     REAL, PARAMETER, DIMENSION(NDIVH) :: USTD_BAS = (/               &
!    & 0.85,0.72,0.59,0.46,0.33,0.16,0.14,0.18,0.28 /) 
!                              !Impact U*t derived from Bagnold data    
!     ----------------------------------------------------------------------

!     PFV: following used for xakae:
!     ----------------------------------------------------------------------
!     PFV: HadGEM3 UM 7.6 settings of Stephanie Woodward (23apr12a, pfv)               ! CMIP5
!     ----------------------------------------------------------------------
!     REAL, PARAMETER :: U_S_MIN      = 1.0E-5     ! LOCAL: minimum val of U_S_STD,    ! 1.0E-5
!                                                  ! below which no dust is produced
!                                                  ! (avoids divide by zero problem)
!     REAL, PARAMETER :: CLAY_MAX     = 0.1        ! LOC max clay fraction.            ! 0.2
!     REAL, PARAMETER :: SNOWMIN      = 1.0E-6     ! LOC min snow depth for Dust       ! 1.0E-6
!     REAL, PARAMETER :: HORIZ_D      = 2.5        ! LOC  factor in Hflux eqn          ! 1.0
!     REAL, PARAMETER :: VERT_C       = 0.01       ! LOC cgs to SI conversion          ! 0.01
!     REAL, PARAMETER :: US_AA        = 0.0        ! LOC Ustar correction (additional) ! 0.0
!     REAL, PARAMETER :: US_AM        = 1.6        ! LOC Ustar correction (multiplic.) ! 1.5
!     REAL, PARAMETER :: UST_AA       = 0.0        ! LOC Ustar_t correction (add)      ! 0.0
!     REAL, PARAMETER :: UST_AM       = 1.0        ! LOC Ustar_t correction (multi.)   ! 1.0
!     REAL, PARAMETER :: SM_CORR      = 0.5        ! soil moist. correction factor     ! 1.7
!     REAL, PARAMETER :: H_OROG_LIMIT = 150.       ! LOC max 1/2pk-trough ht for dust  ! 150.0
!     REAL, PARAMETER, DIMENSION(NDIVH) :: &
!    & USTD_BAS = (/0.85,0.72,0.59,0.46,0.33,0.16,0.14,0.18,0.28/) ! U*t from Bagnold data    
!     ----------------------------------------------------------------------

!     PFV: following used for xakaf: (xakae + 8x HORIZ_D, 2.5 to 20.0)
!     ----------------------------------------------------------------------
!     PFV: Modified HadGEM3 UM 7.6 settings of Stephanie Woodward (2may12a)            ! CMIP5
!     ----------------------------------------------------------------------
!     REAL, PARAMETER :: U_S_MIN      = 1.0E-5     ! LOCAL: minimum val of U_S_STD,    ! 1.0E-5
!                                                  ! below which no dust is produced
!                                                  ! (avoids divide by zero problem)
!     REAL, PARAMETER :: CLAY_MAX     = 0.1        ! LOC max clay fraction.            ! 0.2
!     REAL, PARAMETER :: SNOWMIN      = 1.0E-6     ! LOC min snow depth for Dust       ! 1.0E-6
!     REAL, PARAMETER :: HORIZ_D      = 20.0     !*! LOC  factor in Hflux eqn          ! 1.0
!     REAL, PARAMETER :: VERT_C       = 0.01       ! LOC cgs to SI conversion          ! 0.01
!     REAL, PARAMETER :: US_AA        = 0.0        ! LOC Ustar correction (additional) ! 0.0
!     REAL, PARAMETER :: US_AM        = 1.6        ! LOC Ustar correction (multiplic.) ! 1.5
!     REAL, PARAMETER :: UST_AA       = 0.0        ! LOC Ustar_t correction (add)      ! 0.0
!     REAL, PARAMETER :: UST_AM       = 1.0        ! LOC Ustar_t correction (multi.)   ! 1.0
!     REAL, PARAMETER :: SM_CORR      = 0.5        ! soil moist. correction factor     ! 1.7
!     REAL, PARAMETER :: H_OROG_LIMIT = 150.       ! LOC max 1/2pk-trough ht for dust  ! 150.0
!     REAL, PARAMETER, DIMENSION(NDIVH) :: &
!    & USTD_BAS = (/0.85,0.72,0.59,0.46,0.33,0.16,0.14,0.18,0.28/) ! U*t from Bagnold data    
!     ----------------------------------------------------------------------

!     PFV: following used for xakag: (xakae + 5x HORIZ_D, 2.5 to 12.5)
!     ----------------------------------------------------------------------
!     PFV: Modified HadGEM3 UM 7.6 settings of Stephanie Woodward (2may12a)            ! CMIP5
!     ----------------------------------------------------------------------
!     REAL, PARAMETER :: U_S_MIN      = 1.0E-5     ! LOCAL: minimum val of U_S_STD,    ! 1.0E-5
!                                                  ! below which no dust is produced
!                                                  ! (avoids divide by zero problem)
!     REAL, PARAMETER :: CLAY_MAX     = 0.1        ! LOC max clay fraction.            ! 0.2
!     REAL, PARAMETER :: SNOWMIN      = 1.0E-6     ! LOC min snow depth for Dust       ! 1.0E-6
!     REAL, PARAMETER :: HORIZ_D      = 12.5     !*! LOC  factor in Hflux eqn          ! 1.0
!     REAL, PARAMETER :: VERT_C       = 0.01       ! LOC cgs to SI conversion          ! 0.01
!     REAL, PARAMETER :: US_AA        = 0.0        ! LOC Ustar correction (additional) ! 0.0
!     REAL, PARAMETER :: US_AM        = 1.6        ! LOC Ustar correction (multiplic.) ! 1.5
!     REAL, PARAMETER :: UST_AA       = 0.0        ! LOC Ustar_t correction (add)      ! 0.0
!     REAL, PARAMETER :: UST_AM       = 1.0        ! LOC Ustar_t correction (multi.)   ! 1.0
!     REAL, PARAMETER :: SM_CORR      = 0.5        ! soil moist. correction factor     ! 1.7
!     REAL, PARAMETER :: H_OROG_LIMIT = 150.       ! LOC max 1/2pk-trough ht for dust  ! 150.0
!     REAL, PARAMETER, DIMENSION(NDIVH) :: &
!    & USTD_BAS = (/0.85,0.72,0.59,0.46,0.33,0.16,0.14,0.18,0.28/) ! U*t from Bagnold data    
!     ----------------------------------------------------------------------

!     PFV: for xaldc, xaldj: (SW-UM-7.6 settings with 0.2xSW HORIZ_D, from 2.5 to 0.5)
!     ----------------------------------------------------------------------
!     PFV: Modified HadGEM3 UM 7.6 settings of Stephanie Woodward                      ! CMIP5
!     ----------------------------------------------------------------------
!     REAL, PARAMETER :: U_S_MIN      = 1.0E-5     ! LOCAL: minimum val of U_S_STD,    ! 1.0E-5
!                                                  ! below which no dust is produced
!                                                  ! (avoids divide by zero problem)
!     REAL, PARAMETER :: CLAY_MAX     = 0.1        ! LOC max clay fraction.            ! 0.2
!     REAL, PARAMETER :: SNOWMIN      = 1.0E-6     ! LOC min snow depth for Dust       ! 1.0E-6
!     REAL, PARAMETER :: HORIZ_D      = 0.5      !*! LOC  factor in Hflux eqn          ! 1.0
!     REAL, PARAMETER :: VERT_C       = 0.01       ! LOC cgs to SI conversion          ! 0.01
!     REAL, PARAMETER :: US_AA        = 0.0        ! LOC Ustar correction (additional) ! 0.0
!     REAL, PARAMETER :: US_AM        = 1.6        ! LOC Ustar correction (multiplic.) ! 1.5
!     REAL, PARAMETER :: UST_AA       = 0.0        ! LOC Ustar_t correction (add)      ! 0.0
!     REAL, PARAMETER :: UST_AM       = 1.0        ! LOC Ustar_t correction (multi.)   ! 1.0
!     REAL, PARAMETER :: SM_CORR      = 0.5        ! soil moist. correction factor     ! 1.7
!     REAL, PARAMETER :: H_OROG_LIMIT = 150.       ! LOC max 1/2pk-trough ht for dust  ! 150.0
!     REAL, PARAMETER, DIMENSION(NDIVH) :: &
!    & USTD_BAS = (/0.85,0.72,0.59,0.46,0.33,0.16,0.14,0.18,0.28/) ! U*t from Bagnold data    
!     ----------------------------------------------------------------------

!     PFV: for xaldl: (SW-UM-7.6 settings with 1.0xSW HORIZ_D, 24oct13a)
!     ----------------------------------------------------------------------           ! ac13
!     PFV: Modified HadGEM3 UM 7.6 settings of Stephanie Woodward                      ! CMIP5   SW-hg3-7.6
!     ----------------------------------------------------------------------
      REAL, PARAMETER :: U_S_MIN      = 1.0E-5     ! LOCAL: minimum val of U_S_STD,    ! 1.0E-5  1.0E-5
!                                                  ! below which no dust is produced
!                                                  ! (avoids divide by zero problem)
      REAL, PARAMETER :: CLAY_MAX     = 0.1        ! LOC max clay fraction.            ! 0.2     0.1
      REAL, PARAMETER :: SNOWMIN      = 1.0E-6     ! LOC min snow depth for Dust       ! 1.0E-6  1.0E-6
      REAL, PARAMETER :: HORIZ_D      = 2.5        ! LOC  factor in Hflux eqn          ! 1.0     2.5
      REAL, PARAMETER :: VERT_C       = 0.01       ! LOC cgs to SI conversion          ! 0.01    0.01
      REAL, PARAMETER :: US_AA        = 0.0        ! LOC Ustar correction (additional) ! 0.0     0.0
      REAL, PARAMETER :: US_AM        = 1.6        ! LOC Ustar correction (multiplic.) ! 1.5     1.6
      REAL, PARAMETER :: UST_AA       = 0.0        ! LOC Ustar_t correction (add)      ! 0.0     0.0
      REAL, PARAMETER :: UST_AM       = 1.0        ! LOC Ustar_t correction (multi.)   ! 1.0     1.0
      REAL, PARAMETER :: SM_CORR      = 0.5        ! soil moist. correction factor     ! 1.7     0.5
      REAL, PARAMETER :: H_OROG_LIMIT = 150.       ! LOC max 1/2pk-trough ht for dust  ! 150.0   150.0
      REAL, PARAMETER, DIMENSION(NDIVH) :: &
     & USTD_BAS = (/0.85,0.72,0.59,0.46,0.33,0.16,0.14,0.18,0.28/) ! U*t from Bagnold data    
!     ----------------------------------------------------------------------

      REAL, PARAMETER :: FLAND_LIM = 0.99999 ! LOCAL: used to omit coastal points (PFV 25oct13, via SW)

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
     &, HORIZ_FLUX_1TO6(LAND_PTS,NTILES)                                &
                                         !total over divs 1 to 6
     &, MTOT_LAND(LAND_PTS)                                             &
                            !total of mrel_land over all divs
     &, UST_ALL(LAND_PTS,NTILES,NDIVH)                                  &
                                       !U*t over all divs             
     &, SREL(LAND_PTS,NDIVH)                                            &
                                      !Relative surface area by div             
     &, US(LAND_PTS)                                                    &
                                     ! U* on soil tiles 
     &, SMT(LAND_PTS)                                                   &
                               ! Soil moisture term, used in U*t calc
     &, WORK




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
              HORIZ_FLUX_1TO6(L,M) = 0.
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

         DO IDIV=1,NDIVH
          DO N = 1,TILE_PTS(SOIL)
             L = TILE_INDEX(N,SOIL)
          ENDDO
        ENDDO


        DO N = 1,TILE_PTS(SOIL)
           L = TILE_INDEX(N,SOIL)
           US(L)=0.
           MTOT_LAND(L) = 0.
           DO IDIV = 1,NDIV
             IF (MREL_LAND(L,IDIV)  >   0.)                             &
     &        MTOT_LAND(L) = MTOT_LAND(L) + MREL_LAND(L,IDIV)

           ENDDO
         ENDDO

!
! calculate relative mass for all divs and put into single array
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
!
! pre-calculate soil moisture term for U*t
! work = w-w' (Fecan 99),
! w is % volmetric SM = S_L_M * 100(for %age) * (100cm/10cm) / rho=998
! flim= lim of sm. above which no dust prod, ffac is correction factor
!
        DO N = 1,TILE_PTS(SOIL)
          L = TILE_INDEX(N,SOIL)
!          WORK = SOIL_LAYER_MOISTURE(L,1)  * SM_CORR -                  &
!     &            14.*CLAY_LAND(L)*CLAY_LAND(L) - 17.*CLAY_LAND(L)
!dhb599, 20110615, change made for CABLE case as per Peter Vohralik's suggestion, item 2.
          WORK = SOIL_LAYER_MOISTURE(L,1)  * SM_CORR * (10.0/2.2) -  & !CABLE, dz1=2.2cm
     &            14.*CLAY_LAND(L)*CLAY_LAND(L) - 17.*CLAY_LAND(L)     !CABLE
          IF (WORK .GT. 0.) THEN
            SMT(L) =  ((1.0 + 1.21* WORK **.68 )**.5)
          ELSE
            SMT(L) = 1.0
           ENDIF
         ENDDO

!
!lopp over pts in each div calculating dust flux, if any
!
      DO IDIV = 1,NDIVH
        DO N = 1,TILE_PTS(SOIL)
          L = TILE_INDEX(N,SOIL)
!
! horizontal flux
!
          UST_ALL(L,SOIL,IDIV)=                                         &
     &          USTD_BAS(IDIV)*SMT(L)*UST_AM                            
          US(L)=U_S_STD_TILE(L,SOIL)*US_AM
          IF ( (UST_ALL(L,SOIL,IDIV) < US(L)) .AND.                     &
     &     (FLAND(L) >  FLAND_LIM) .AND.                                & ! PFV, 25oct13, omit coastal pts (via SW)
     &     (TSTAR_TILE(L,SOIL) > ZERODEGC) .AND.                        &
     &     (SNOW_TILE(L,SOIL)  <   SNOWMIN) .AND.                       &
     &     (MREL_LAND_ALL(L,IDIV)  >   0.) .AND.                        &
     &     (HO2R2_OROG(L)  <=  H_OROG_LIMIT) .AND.                      &
!    &     ((SOIL_LAYER_MOISTURE(L,1)*SM_CORR <                         &
!    &     (CLAY_LAND(L)+.12)/.03)) .AND.                               &
!dhb599, 20110615, change made for CABLE case as per Peter Vohralik's suggestion, item 2.
     &     ((SOIL_LAYER_MOISTURE(L,1)*SM_CORR*(10.0/2.2) <  &  !CABLE, dz1=2.2cm
     &     (CLAY_LAND(L)+.12)/.03)) .AND.                   &  !CABLE
     &     (US(L) > U_S_MIN) ) THEN                                      
             RATIO = UST_ALL(L,SOIL,IDIV) / US(L)                       
             HORIZ_FLUX_ALL(L,SOIL,IDIV) = HORIZ_C * RHOSTAR_LAND(L) *  &
     &       HORIZ_D * US(L)**3. * (1.0 + RATIO) *                      &
     &       (1.0 - RATIO*RATIO) * MREL_LAND_ALL(L,IDIV) / G
             TOT_HORIZ_FLUX_TILE(L,SOIL)=                               &
     &       TOT_HORIZ_FLUX_TILE(L,SOIL)+HORIZ_FLUX_ALL(L,SOIL,IDIV)
           ENDIF
        ENDDO !TILE_PTS
      ENDDO !NDIV

        DO N = 1,TILE_PTS(SOIL)
          L = TILE_INDEX(N,SOIL)
            HORIZ_FLUX_789(L,SOIL)=HORIZ_FLUX_ALL(L,SOIL,7)+            &
     &       HORIZ_FLUX_ALL(L,SOIL,8)+HORIZ_FLUX_ALL(L,SOIL,9)
            HORIZ_FLUX_1TO6(L,SOIL)=HORIZ_FLUX_ALL(L,SOIL,1)+           &
     &        HORIZ_FLUX_ALL(L,SOIL,2)+HORIZ_FLUX_ALL(L,SOIL,3)+        &
     &        HORIZ_FLUX_ALL(L,SOIL,4)+HORIZ_FLUX_ALL(L,SOIL,5)+        &
     &        HORIZ_FLUX_ALL(L,SOIL,6)
         ENDDO
!
! vertical flux
!
      DO IDIV = 1,NDIV
        DO N = 1,TILE_PTS(SOIL)
          L = TILE_INDEX(N,SOIL)
          IF ((MREL_LAND(L,IDIV)  >   0.) .AND.                         &
     &   (HORIZ_FLUX_1TO6(L,SOIL )   >   0.) .AND.                      &
     &  (SNOW_TILE(L,SOIL)  <   SNOWMIN)) THEN
           DUST_FLUX_TILE(L,SOIL,IDIV) =                                &
     &   VERT_C *                                                       &
     &   ( 1.0 + HORIZ_FLUX_789(L,SOIL)/HORIZ_FLUX_1TO6(L,SOIL ) ) *    &
     &    HORIZ_FLUX_ALL(L,SOIL,IDIV) *                                 &
     &   10**(VERT_A * MIN(CLAY_LAND(L),CLAY_MAX) + VERT_B)
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
!         HORIZ_FLUX_TILE(L,M,IDIV)=HORIZ_FLUX_TILE(L,M,IDIV)*FLAND(L)
         DUST_FLUX_TILE(L,M,IDIV) = DUST_FLUX_TILE(L,M,IDIV)*FLAND(L)
        ENDDO !TILE_PTS
       ENDDO !NTILES
      ENDDO !NDIV
!
!   diagnostics (note only have U_S_T for first 6 divs here)
!
      DO IDIV = 1,NDIV
        DO N = 1,TILE_PTS(SOIL)
         L = TILE_INDEX(N,SOIL)
          U_S_T_TILE(L,SOIL,IDIV)=UST_ALL(L,SOIL,IDIV)
        ENDDO !TILE_PTS
      ENDDO !NDIV
!
!
!
      RETURN
      END SUBROUTINE DUST_SRCE
