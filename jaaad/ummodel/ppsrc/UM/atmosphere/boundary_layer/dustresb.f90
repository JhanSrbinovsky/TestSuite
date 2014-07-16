
!
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! subroutine DUSTRESB
!
!
! Purpose:
!   To calculate the surface layer resistance for mineral dust
!
! Called by sfexch
!
! Current owners of code: S.Woodward
!
! History:
! Version     Date     Comment
! -------     ----     -------
!
!   5.5      12/02/03  Original code   S Woodward
!   6.2      11/05/06  Fix for bit-non-reproducibility   A. Malcolm
!   6.4      15/01/07  Malcolm McVean's portability fix   S.Woodward
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!
! Documentation: "Modelling the atmospheric lifecycle..."
!                 Woodward, JGR106, D16, pp18155-18166
!---------------------------------------------------------------------
!
       SUBROUTINE DUSTRESB(                                             &
     &  ROW_LENGTH,ROWS,                                                &
     &  PSTAR,TSTAR,RHOSTAR,ARESIST,VSHR,CD_STD_DUST,                   &
     &  R_B_DUST                                                        &
     &  )
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

      INTEGER                                                           &
              !IN
     & ROW_LENGTH                                                       &
                  !IN
     &,ROWS       !IN

      REAL                                                              &
           !IN
     & PSTAR(ROW_LENGTH,ROWS)                                           &
                                    !IN surface pressure
     &,TSTAR(ROW_LENGTH,ROWS)                                           &
                                    !IN surface temperature
     &,RHOSTAR(ROW_LENGTH,ROWS)                                         &
                                    !IN surface air density
     &,ARESIST(ROW_LENGTH,ROWS)                                         &
                                    !IN aerodynamic resistance
     &,VSHR(ROW_LENGTH,ROWS)                                            &
                                    !IN surface to lowest lev windspeed
!                                   !   difference
     &,CD_STD_DUST(ROW_LENGTH,ROWS) !IN surface transfer coeffient for
!                                   !   momentum, excluding orographic
!                                   !   form drag

      REAL                                                              &
           !OUT
     & R_B_DUST(ROW_LENGTH,ROWS,NDIV) !OUT surface layer resistance for
!                                     !    mineral dust

!     local variables

      INTEGER                                                           &
     & IDIV                                                             &
            !loop counter, dust divisions
     &,I                                                                &
            !loop counter
     &,J                                                                &
            !loop counter
     &,LEV1 !number of levels for vstokes calculation

      REAL                                                              &
     & NU(ROW_LENGTH,ROWS)                                              &
                                 !kinematic viscosity
     &,ETAA(ROW_LENGTH,ROWS)                                            &
                                 !dynamic viscosity of air
     &,LAMDAA(ROW_LENGTH,ROWS)                                          &
                                 !mean free path of air molecules
     &,VSTOKES1(ROW_LENGTH,ROWS)                                        &
                                 !gravitational settling velocity, lev1
     &,NSTOKES(ROW_LENGTH,ROWS)                                         &
                                 !stokes number = VstokesVshrVshr/nu g
     &,NSCHMIDT(ROW_LENGTH,ROWS)                                        &
                                 !schmidt number = nu/diffusivit
     &,TC(ROW_LENGTH,ROWS)                                              &
                                 !temperature in deg C
     &,ALPHACCF(ROW_LENGTH,ROWS)                                        &
                                 !alpha in cunningham correction factor
     &,CCF(ROW_LENGTH,ROWS)                                             &
                                 !Cunningham correction factor
     &,FVSQ(ROW_LENGTH,ROWS)                                            &
                                 !friction velocity squared
     &,WORK(ROW_LENGTH,ROWS)                                            &
                                 !workspace
     &,STOKES_EXP                                                       &
                                 !stokes term in R_B_DUST equation
     &,SMALLP                    !small +ve number, negligible compared to 1

!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
!*L------------------COMDECK C_PI---------------------------------------
!LL
!LL 4.0 19/09/95  New value for PI. Old value incorrect
!LL               from 12th decimal place. D. Robinson
!LL 5.1 7/03/00   Fixed/Free format P.Selwood
!LL

      ! Pi
      Real, Parameter :: Pi                 = 3.14159265358979323846

      ! Conversion factor degrees to radians
      Real, Parameter :: Pi_Over_180        = Pi/180.0

      ! Conversion factor radians to degrees
      Real, Parameter :: Recip_Pi_Over_180  = 180.0/Pi

!*----------------------------------------------------------------------
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------
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
!C_DUSTGRAV.............................................................
! Description:
! Contains parameters for mineral dust gravitational settling
! Current Code Owner: Stephanie Woodward
!
! History:
! Version  Date     Comment
! -------  ----     -------
!  1      12/02/03  Original Code.   Stephanie Woodward
!
      REAL,PARAMETER :: ACCF=1.257 ! Cunningham correction factor term A
      REAL,PARAMETER :: BCCF=0.4 ! Cunningham correction factor term B
      REAL,PARAMETER :: CCCF=-1.1 ! Cunningham correction factor term C
!.......................................................................
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
!
       EXTERNAL VGRAV
!
!... epsilon() is defined as almost negligible, so eps/100 is negligible
!
      SMALLP = epsilon(1.0) / 100.0
!
!...calc stokes number, schmidt number and finally resistance
!
      LEV1=1

      DO IDIV=1,NDIV
! DEPENDS ON: vgrav
        CALL VGRAV(                                                     &
     &  ROW_LENGTH,ROWS,LEV1,DREP(IDIV),RHOP,PSTAR,TSTAR,               &
     &  VSTOKES1,CCF,ETAA                                               &
     &  )

!CDIR NOVECTOR
        DO J = 1,ROWS
          DO I= 1,ROW_LENGTH
            NSCHMIDT(I,J)=3.*PI*ETAA(I,J)*ETAA(I,J)*DREP(IDIV)/         &
     &       (RHOSTAR(I,J)*BOLTZMANN*TSTAR(I,J)*CCF(I,J))
            NSTOKES(I,J)=VSTOKES1(I,J)*CD_STD_DUST(I,J)*RHOSTAR(I,J)*   &
     &       VSHR(I,J)*VSHR(I,J)/(ETAA(I,J)*G)
            ! Avoid underflow in Stokes term by setting to zero if 
            ! negligible compared to Schmidt term, i.e., if NSTOKES
            ! is too small.
            IF ( 3.0 / NSTOKES(I,J) <                                   &
                 - LOG10( SMALLP *NSCHMIDT(I,J)**(-2./3.) ) ) THEN
               STOKES_EXP = 10.**(-3./NSTOKES(I,J))
            ELSE
               STOKES_EXP = 0.0
            ENDIF
            R_B_DUST(I,J,IDIV)=1./( SQRT(CD_STD_DUST(I,J)) *            &
     &       (NSCHMIDT(I,J)**(-2./3.)+STOKES_EXP) )
          ENDDO !ROW_LENGTH
        ENDDO !ROWS
      ENDDO !NDIV
      
      RETURN
      END SUBROUTINE DUSTRESB
