
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    Subroutine VGRAV -------------------------------------------------
!
! Purpose: To calculate the gravitational sedimentation velocity of
!          tracer particles according to Stoke's law, including the
!          Cunningham correction factor.
!
! Current owners of code:                 S Woodward, M Woodage
!
! History:
! Version    Date     Comment
! -------    ----     -------
!   4.4    03/10/97   Original code        S Woodward, M Woodage
!   5.5    12/02/03   Updated for vn 5.5   S Woodward
!
!
! Code description:
!  Language: FORTRAN77 + extensions
!  Programming standard: UMDP 3 Vn 6
!
! System components covered:
!
! System task:
!
!Documentation: Ref. Pruppacher & Klett
!                    Microphysics of clouds & ppn    1978,1980 edns.
!
!-----------------------------------------------------------------------
!
      SUBROUTINE VGRAV(                                                 &
     & ROW_LENGTH,ROWS,NLEVS,DIAM,RHOP,P,T,                             &
     & VSTOKES,CCF,ETAA)




!
      implicit none
!
!
      INTEGER ROW_LENGTH         !IN row length
      INTEGER ROWS               !IN number of rows
      INTEGER NLEVS              !IN number of levels
!
      REAL DIAM                  !IN particle diameter
      REAL RHOP                  !IN particles density
      REAL P(ROW_LENGTH,ROWS,NLEVS)!IN pressure
      REAL T(ROW_LENGTH,ROWS,NLEVS)!IN temperature
!
      REAL VSTOKES(ROW_LENGTH,ROWS,NLEVS) !OUT sedimentation velocity
      REAL ETAA(ROW_LENGTH,ROWS,NLEVS)!OUT viscosity of air
      REAL CCF(ROW_LENGTH,ROWS,NLEVS) !OUT cunningham correction factor
!
!
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
!
! local variables
!
      INTEGER ILEV               !LOC loop counter for levels
      INTEGER I                  !LOC loop counter
      INTEGER J                  !LOC loop counter
      INTEGER K                  !LOC loop counter
!
      REAL TC(ROW_LENGTH,ROWS)   !LOC temperature in deg C
      REAL LAMDAA(ROW_LENGTH,ROWS)!LOC mean free path of particle
      REAL ALPHACCF(ROW_LENGTH,ROWS)!LOC
!
! Calculate viscosity of air (Pruppacher & Klett p.323)
      DO ILEV=1,NLEVS
        DO J=1,ROWS
          DO I = 1,ROW_LENGTH
           TC(I,J)=T(I,J,ILEV)-ZERODEGC
           IF (TC(I,J)  >=  0.) THEN
            ETAA(I,J,ILEV)=(1.718+0.0049*TC(I,J))*1.E-5
           ELSE
            ETAA(I,J,ILEV)=                                             &
     &            (1.718+0.0049*TC(I,J)-1.2E-5*TC(I,J)*TC(I,J))*1.E-5
           ENDIF
!
          ENDDO !ROW_LENGTH
        ENDDO !ROWS
      ENDDO !NLEVS
!
      DO ILEV=1,NLEVS
        DO J=1,ROWS
          DO I=1,ROW_LENGTH
!
! Calculate mean free path of particle (Pruppacher & Klett p.323)
           LAMDAA(I,J)=MFP_REF*PREF_MFP*T(I,J,ILEV)/                    &
     &      (P(I,J,ILEV)*TREF_MFP)
! Calculate Cunningham correction factor(Pruppacher & Klett p.361)
           ALPHACCF(I,J)=ACCF+BCCF*EXP(CCCF*DIAM*.5/LAMDAA(I,J))
           CCF(I,J,ILEV)=(1.+ALPHACCF(I,J)*LAMDAA(I,J)/(.5*DIAM))
! Calculate sedimentation velocity (Pruppacher & Klett p.362)
           VSTOKES(I,J,ILEV)=CCF(I,J,ILEV)*(DIAM*DIAM*G*RHOP)/          &
     &             (18.*ETAA(I,J,ILEV))
!
          ENDDO !ROW_LENGTH
        ENDDO !ROWS
      ENDDO !NLEV
!
      RETURN
      END SUBROUTINE VGRAV
