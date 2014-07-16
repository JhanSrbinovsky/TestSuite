
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
!C_DUSTGEN_CAM.........................................................
! Description: Contains parameters for mineral dust generation
!
! This version includes values of the constants that are currently used 
! in the NWP dust setup in the Crisis Area Models, which is based on
! an old configuration of the HadCM3 dust code. In the future, we hope
! to unify the code used in  HadGEM and NWP.
!
! Current Code Owner: Stephanie Woodward
! 
!----------------------------------------------------------------------
!
!
!     ! A,B and C in calculation of threshold friction velocity
!     ! Use value of A suitable for Iverson and White T*t(dry) 
!     ! and value of C tuned for use in the Southern Asia CAM 
      REAL, PARAMETER :: UST_A = 0.05
      REAL, PARAMETER :: UST_B = 0
      REAL, PARAMETER :: UST_C = -0.10

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
